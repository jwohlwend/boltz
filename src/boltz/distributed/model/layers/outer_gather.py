# SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
# SPDX-License-Identifier: MIT
#
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.


import math
from typing import Dict, List, Tuple

import torch
import torch.distributed as dist
from torch.distributed._tensor import DTensor, Shard
from torch.distributed.tensor import Partial, Replicate

from boltz.distributed.utils import update_exhaustive_strides


def outer_gather_backward(grad_output, z_shape, idx_q, idx_k, axis, idx_q_mask=None, idx_k_mask=None):
    """
    Backward pass for outer_gather.

    Args:
        grad_output: (..., K, W, H, ...)
        z_shape: Shape of original input z
        idx_q: (..., K, W)
        idx_k: (..., K, H)
        axis: Axis where N starts in z
        idx_q_mask: Optional mask for idx_q of shape (..., K, W). True = valid, False = invalid.
        idx_k_mask: Optional mask for idx_k of shape (..., K, H). True = valid, False = invalid.

    Returns:
        Gradient w.r.t. z
    """
    # Early exit for empty z
    if math.prod(z_shape) == 0:
        # z was empty, gradient is zeros of z_shape
        return torch.zeros(z_shape, dtype=grad_output.dtype, device=grad_output.device)

    batch_shape_z = z_shape[:axis]
    feature_shape_z = z_shape[axis + 2 :]
    if grad_output.shape[:axis] != batch_shape_z:
        # in the forward pass, z's 'axis' is replaced in the output by K, W, H so grad_output.shape[axis:axis+3] is (K, W, H)
        raise ValueError(
            f"grad_output.shape[:axis] must match z_shape[:axis] but got {grad_output.shape[:axis]} vs {z_shape[:axis]}"
        )

    if grad_output.shape[axis + 3 :] != feature_shape_z:
        raise ValueError(
            f"grad_output.shape[axis + 3 :] must match z_shape[axis + 2 :] but got {grad_output.shape[axis + 3 :]} vs {feature_shape_z}"
        )

    if grad_output.shape[: axis + 2] != idx_q.shape:
        # (..., K, W) must match
        raise ValueError(
            f"grad_output.shape[:axis + 2] must match idx_q.shape but got {grad_output.shape[: axis + 2]} vs {idx_q.shape}"
        )

    if grad_output.shape[: axis + 1] + (grad_output.shape[axis + 2],) != idx_k.shape:
        raise ValueError(
            f"grad_output.shape[:axis + 1] + (grad_output.shape[axis + 2],) must match idx_k.shape but got "
            f"{grad_output.shape[:axis] + (grad_output.shape[axis + 2],)} vs {idx_k.shape}"
        )

    # Validate mask shapes if provided
    if idx_q_mask is not None and idx_q_mask.shape != idx_q.shape:
        raise ValueError(f"idx_q_mask shape {idx_q_mask.shape} must match idx_q shape {idx_q.shape}")
    if idx_k_mask is not None and idx_k_mask.shape != idx_k.shape:
        raise ValueError(f"idx_k_mask shape {idx_k_mask.shape} must match idx_k shape {idx_k.shape}")

    ndim_z = len(z_shape)
    flatten_leading_dims = len(batch_shape_z) >= 2  # w at least two leading axes
    flatten_trailing_dims = len(feature_shape_z) >= 2  # w at least two trailing axes
    has_leading_dims = len(batch_shape_z) > 0
    has_trailing_dims = len(feature_shape_z) > 0

    # 1. Normalize Axis
    ndim_z = len(z_shape)
    if axis < 0:
        axis += ndim_z

    # Re-construct broadcasted indices (avoid flattening if possible)
    # idx_q: (B, K, W)
    # idx_k: (B, K, H)

    # Broadcast to grid (B, K, W, H)
    K = idx_q.shape[-2]
    W = idx_q.shape[-1]
    H = idx_k.shape[-1]

    # Reshape indices to (B, K, W) and (B, K, H)
    idx_q_flat = idx_q
    idx_k_flat = idx_k
    idx_q_mask_flat = idx_q_mask
    idx_k_mask_flat = idx_k_mask
    if flatten_leading_dims:
        idx_q_flat = idx_q_flat.flatten(0, -3)  # (B, K, W)
        idx_k_flat = idx_k_flat.flatten(0, -3)  # (B, K, H)
        if idx_q_mask_flat is not None:
            idx_q_mask_flat = idx_q_mask_flat.flatten(0, -3)
        if idx_k_mask_flat is not None:
            idx_k_mask_flat = idx_k_mask_flat.flatten(0, -3)

    if not has_leading_dims:
        idx_q_flat = idx_q_flat.unsqueeze(0)
        idx_k_flat = idx_k_flat.unsqueeze(0)
        if idx_q_mask_flat is not None:
            idx_q_mask_flat = idx_q_mask_flat.unsqueeze(0)
        if idx_k_mask_flat is not None:
            idx_k_mask_flat = idx_k_mask_flat.unsqueeze(0)

    B = idx_q_flat.shape[0]

    # For masked indices, clamp to valid range (0) to avoid index errors
    if idx_q_mask_flat is not None:
        idx_q_flat = torch.where(idx_q_mask_flat, idx_q_flat, torch.zeros_like(idx_q_flat))
    if idx_k_mask_flat is not None:
        idx_k_flat = torch.where(idx_k_mask_flat, idx_k_flat, torch.zeros_like(idx_k_flat))

    # Broadcast to (B, K, W, H)
    q_broad = idx_q_flat.unsqueeze(-1).expand(-1, -1, -1, H)
    k_broad = idx_k_flat.unsqueeze(-2).expand(-1, -1, W, -1)

    # Batch indices (B, K, W, H)
    batch_idx = torch.arange(B, device=grad_output.device).reshape(B, 1, 1, 1).expand(-1, K, W, H)

    # Compute linear indices for (B, N, M) -> (B*N*M)
    # We are updating (B, N, M, D)
    N = z_shape[axis]
    M = z_shape[axis + 1]
    linear_idx = batch_idx * (N * M) + q_broad * M + k_broad  # (B, K, W, H)
    linear_idx_flat = linear_idx.reshape(-1)  # (B*K*W*H)

    # (..., K, W, H, ...) -> (..., K, W, H, D) -> (B* K * W * H, D)
    grad_source_flat = grad_output
    if flatten_trailing_dims:
        # must flatten trailing dims before leading dims to avoid axis offset changes
        grad_source_flat = grad_source_flat.flatten(axis + 3, -1)
    if not has_trailing_dims:
        grad_source_flat = grad_source_flat.unsqueeze(-1)
    grad_source_flat = grad_source_flat.flatten(0, -2)

    # Compute combined mask and zero out gradients at invalid positions
    # combined_mask: (B, K, W, H) - position is valid only if both q and k indices are valid
    if idx_q_mask_flat is not None or idx_k_mask_flat is not None:
        if idx_q_mask_flat is not None and idx_k_mask_flat is not None:
            # Broadcast masks: q_mask (B,K,W) -> (B,K,W,H), k_mask (B,K,H) -> (B,K,W,H)
            combined_mask = idx_q_mask_flat.unsqueeze(-1) & idx_k_mask_flat.unsqueeze(-2)
        elif idx_q_mask_flat is not None:
            combined_mask = idx_q_mask_flat.unsqueeze(-1).expand(-1, -1, -1, H)
        else:
            combined_mask = idx_k_mask_flat.unsqueeze(-2).expand(-1, -1, W, -1)
        combined_mask_flat = combined_mask.reshape(-1, 1).to(grad_source_flat.dtype)
        grad_source_flat = grad_source_flat * combined_mask_flat

    D = grad_source_flat.shape[-1]

    # Initialize flat grad_z
    grad_z_flat = torch.zeros((B * N * M, D), device=grad_output.device, dtype=grad_output.dtype)

    # Scatter add
    grad_z_flat.index_add_(0, linear_idx_flat, grad_source_flat)

    # Unflatten grad_z
    # (B, N, M, D) -> (..., N, M, ...)
    grad_z = grad_z_flat.reshape(z_shape)

    return grad_z


class OuterGather(torch.autograd.Function):
    @staticmethod
    def forward(ctx, z, idx_q, idx_k, axis=1, idx_q_mask=None, idx_k_mask=None):
        """
        Perform outer gather operation: z[b, q, k] for all q in idx_q, k in idx_k.

        Args:
            z: (..., N, M, ...)
            idx_q: (..., K, W)
            idx_k: (..., K, H)
            axis: The dimension index of the first N in z.
            idx_q_mask: Optional mask for idx_q of shape (..., K, W). True = valid, False = invalid.
            idx_k_mask: Optional mask for idx_k of shape (..., K, H). True = valid, False = invalid.

        Returns:
            Tensor of shape (..., K, W, H, ...)
        """
        # 1. Normalize Axis
        ndim_z = z.ndim

        if ndim_z < 2:
            raise ValueError(f"z must have at least 2 dimensions but got {ndim_z}")

        if idx_q.ndim < 2:
            raise ValueError(f"idx_q must have at least 2 dimensions but got {idx_q.ndim}")

        if idx_k.ndim < 2:
            raise ValueError(f"idx_k must have at least 2 dimensions but got {idx_k.ndim}")

        if axis < 0:
            axis += ndim_z

        if not (0 <= axis < z.ndim):
            raise ValueError(f"Axis must be in range [0, {z.ndim - 1}] but got {axis}")

        # 2. Check shapes (Strict Equality)
        # idx_q: (..., K, W)
        # idx_k: (..., K, H)
        # z: (..., N, M, ...)

        batch_shape_idx = idx_q.shape[:-2]
        batch_shape_z = z.shape[:axis]
        feature_shape_z = z.shape[axis + 2 :]

        if batch_shape_z != batch_shape_idx:
            raise ValueError(
                f"Leading dimensions must match exactly but got: z {batch_shape_z} vs idx_q {batch_shape_idx}"
            )
        if idx_k.shape[:-1] != idx_q.shape[:-1]:
            raise ValueError(
                f"All dimensions but the last must match exactly but got: idx_k {idx_k.shape[:-2]} vs idx_q {batch_shape_idx}"
            )

        # Validate masks if provided
        if idx_q_mask is not None and idx_q_mask.shape != idx_q.shape:
            raise ValueError(f"idx_q_mask shape {idx_q_mask.shape} must match idx_q shape {idx_q.shape}")
        if idx_k_mask is not None and idx_k_mask.shape != idx_k.shape:
            raise ValueError(f"idx_k_mask shape {idx_k_mask.shape} must match idx_k shape {idx_k.shape}")
        if idx_q_mask is not None and idx_q_mask.dtype != torch.bool:
            raise TypeError(
                f"idx_q_mask must have dtype torch.bool, got {idx_q_mask.dtype}. Use mask.bool() to convert."
            )
        if idx_k_mask is not None and idx_k_mask.dtype != torch.bool:
            raise TypeError(
                f"idx_k_mask must have dtype torch.bool, got {idx_k_mask.dtype}. Use mask.bool() to convert."
            )

        K = idx_q.shape[-2]
        W = idx_q.shape[-1]
        H = idx_k.shape[-1]

        has_leading_dims = len(batch_shape_z) > 0

        # Reshape z to (B, N, M, D)
        flatten_leading_dims = len(batch_shape_z) >= 2  # w at least two leading axes
        flatten_trailing_dims = len(feature_shape_z) >= 2  # w at least two trailing axes
        z_flat = z
        if flatten_trailing_dims:
            # must flatten trailing dims before leading dims to avoid axis offset changes
            z_flat = z_flat.flatten(axis + 2, -1)
        if flatten_leading_dims:
            z_flat = z_flat.flatten(0, axis - 1)

        # Reshape indices to (B, K, W) and (B, K, H)
        idx_q_flat = idx_q
        idx_k_flat = idx_k
        idx_q_mask_flat = idx_q_mask
        idx_k_mask_flat = idx_k_mask
        if flatten_leading_dims:
            idx_q_flat = idx_q_flat.flatten(0, -3)
            idx_k_flat = idx_k_flat.flatten(0, -3)
            if idx_q_mask_flat is not None:
                idx_q_mask_flat = idx_q_mask_flat.flatten(0, -3)
            if idx_k_mask_flat is not None:
                idx_k_mask_flat = idx_k_mask_flat.flatten(0, -3)

        if not has_leading_dims:
            z_flat = z_flat.unsqueeze(0)
            idx_q_flat = idx_q_flat.unsqueeze(0)
            idx_k_flat = idx_k_flat.unsqueeze(0)
            if idx_q_mask_flat is not None:
                idx_q_mask_flat = idx_q_mask_flat.unsqueeze(0)
            if idx_k_mask_flat is not None:
                idx_k_mask_flat = idx_k_mask_flat.unsqueeze(0)
        B = z_flat.shape[0]

        # Early exit for empty z tensor
        if z_flat.numel() == 0:
            # Validate: the combined mask (outer-AND of q and k masks) must be all-False when z is empty.
            # A position (w, h) is valid only if BOTH idx_q_mask[w] AND idx_k_mask[h] are True.
            # We must check the combined mask, not each mask separately, because:
            # - If idx_q_mask is all-False, no output positions are valid regardless of idx_k_mask
            # - If idx_k_mask is all-False, no output positions are valid regardless of idx_q_mask
            has_valid_output_positions = False
            if idx_q_mask_flat is not None or idx_k_mask_flat is not None:
                if idx_q_mask_flat is not None and idx_k_mask_flat is not None:
                    # Joint validation: outer-AND of the two masks
                    combined_mask = idx_q_mask_flat.unsqueeze(-1) & idx_k_mask_flat.unsqueeze(-2)
                    has_valid_output_positions = combined_mask.any().item()
                elif idx_q_mask_flat is not None:
                    has_valid_output_positions = idx_q_mask_flat.any().item()
                else:
                    has_valid_output_positions = idx_k_mask_flat.any().item()
            else:
                # No masks provided means all positions are implicitly valid
                has_valid_output_positions = idx_q_flat.numel() > 0 and idx_k_flat.numel() > 0

            if has_valid_output_positions:
                raise ValueError(
                    "z is empty but combined mask (idx_q_mask & idx_k_mask) contains valid entries. "
                    "This is a logical error - cannot gather from empty tensor with valid indices."
                )
            # z is empty - return zeros of correct shape
            out_shape = batch_shape_z + (K, W, H) + feature_shape_z
            out = torch.zeros(out_shape, dtype=z.dtype, device=z.device)
            # Save context for backward (same as normal path)
            tensors_to_save = [idx_q, idx_k]
            if idx_q_mask is not None:
                tensors_to_save.append(idx_q_mask)
            if idx_k_mask is not None:
                tensors_to_save.append(idx_k_mask)
            ctx.save_for_backward(*tensors_to_save)
            ctx.has_q_mask = idx_q_mask is not None
            ctx.has_k_mask = idx_k_mask is not None
            ctx.z_shape = z.shape
            ctx.axis = axis
            return out

        # For masked indices, clamp to valid range (0) to avoid index errors
        if idx_q_mask_flat is not None:
            idx_q_flat = torch.where(idx_q_mask_flat, idx_q_flat, torch.zeros_like(idx_q_flat))
        if idx_k_mask_flat is not None:
            idx_k_flat = torch.where(idx_k_mask_flat, idx_k_flat, torch.zeros_like(idx_k_flat))

        # Create broadcasted indices for gather
        # q_broad: (B, K, W, H)
        q_broad = idx_q_flat.unsqueeze(-1).expand(B, K, W, H)
        # k_broad: (B, K, W, H)
        k_broad = idx_k_flat.unsqueeze(-2).expand(B, K, W, H)
        # batch_idx: (B, K, W, H)
        batch_idx = torch.arange(B, device=z.device).reshape(B, 1, 1, 1).expand(B, K, W, H)

        # Gather: z_flat[b, q, k]
        # Result: (B, K, W, H, D)
        out_flat = z_flat[batch_idx, q_broad, k_broad]

        # Zero out invalid positions
        # combined_mask: (B, K, W, H) - position is valid only if both q and k indices are valid
        if idx_q_mask_flat is not None or idx_k_mask_flat is not None:
            if idx_q_mask_flat is not None and idx_k_mask_flat is not None:
                # Broadcast masks: q_mask (B,K,W) -> (B,K,W,H), k_mask (B,K,H) -> (B,K,W,H)
                combined_mask = idx_q_mask_flat.unsqueeze(-1) & idx_k_mask_flat.unsqueeze(-2)
            elif idx_q_mask_flat is not None:
                combined_mask = idx_q_mask_flat.unsqueeze(-1).expand(-1, -1, -1, H)
            else:
                combined_mask = idx_k_mask_flat.unsqueeze(-2).expand(-1, -1, W, -1)
            # Expand to match out_flat: (B, K, W, H) or (B, K, W, H, D) depending on trailing dims
            has_trailing_dims = len(feature_shape_z) > 0
            if has_trailing_dims:
                combined_mask = combined_mask.unsqueeze(-1)
            out_flat = out_flat * combined_mask.to(out_flat.dtype)

        # Reshape to final output
        out = out_flat
        if flatten_trailing_dims:
            out = out.unflatten(-1, z.shape[axis + 2 :])
        if flatten_leading_dims:
            out = out.unflatten(0, batch_shape_z)
        if not has_leading_dims:
            out = out.squeeze(0)

        # Save for backward
        tensors_to_save = [idx_q, idx_k]
        if idx_q_mask is not None:
            tensors_to_save.append(idx_q_mask)
        if idx_k_mask is not None:
            tensors_to_save.append(idx_k_mask)
        ctx.save_for_backward(*tensors_to_save)
        ctx.has_q_mask = idx_q_mask is not None
        ctx.has_k_mask = idx_k_mask is not None
        ctx.z_shape = z.shape
        ctx.axis = axis

        return out

    @staticmethod
    def backward(ctx, grad_output):
        saved = ctx.saved_tensors
        idx_q = saved[0]
        idx_k = saved[1]
        idx = 2
        idx_q_mask = saved[idx] if ctx.has_q_mask else None
        if ctx.has_q_mask:
            idx += 1
        idx_k_mask = saved[idx] if ctx.has_k_mask else None

        grad_z = outer_gather_backward(grad_output, ctx.z_shape, idx_q, idx_k, ctx.axis, idx_q_mask, idx_k_mask)
        return grad_z, None, None, None, None, None


def outer_gather(z, one_hot_q, one_hot_k, axis=1, one_hot_q_mask=None, one_hot_k_mask=None):
    """
    Efficient gather-based equivalent to einsum for window batching index selection.

    Args:
        z: Input tensor (..., N, M, ...)
        one_hot_q: One-hot query indices (..., K, W, N)
        one_hot_k: One-hot key indices (..., K, H, M)
        axis: The dimension index of the first N in z.
        one_hot_q_mask: Optional mask for one_hot_q of shape (..., K, W). True = valid, False = invalid.
        one_hot_k_mask: Optional mask for one_hot_k of shape (..., K, H). True = valid, False = invalid.

    Returns:
        Tensor of shape (..., K, W, H, ...)
    """
    if z.shape[axis] != one_hot_q.shape[-1] or z.shape[axis + 1] != one_hot_k.shape[-1]:
        raise ValueError(
            f"z.shape[axis] must match one_hot_q.shape[-1] and z.shape[axis + 1] must match one_hot_k.shape[-1] but got "
            f"{z.shape[axis : axis + 2]} vs {one_hot_q.shape[-1:]} and {one_hot_k.shape[-1:]}"
        )
    # Condense one-hot to indices
    idx_q = one_hot_q.argmax(dim=-1)
    idx_k = one_hot_k.argmax(dim=-1)
    return OuterGather.apply(z, idx_q, idx_k, axis, one_hot_q_mask, one_hot_k_mask)


def compute_interval_overlap(intervals_a: torch.Tensor, intervals_b: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
    """
    Computes the intersection (overlap) of two sets of n-dimensional intervals where
    the last axis represents [start, end) coordinates (exclusive end).

    Args:
        intervals_a (torch.Tensor): A tensor of shape (..., n_dim, 2).
            The last dimension contains [start, end).
        intervals_b (torch.Tensor): A tensor of shape (..., n_dim, 2).
            Must be broadcastable to intervals_a.shape.

    Returns:
        tuple[torch.Tensor, torch.Tensor]:
            - intervals_overlap (torch.Tensor): The computed intersection intervals with shape
              matching the broadcasted shape of inputs (..., n_dim, 2).
              Contains [overlap_start, overlap_end).
            - mask (torch.Tensor): A boolean mask of shape (...,) corresponding to the leading
              dimensions. It is True if the n-dimensional overlap is valid and non-empty.
              A valid overlap is defined by: all(overlap_end > overlap_start) across the n_dim axis.
    """
    start_a = intervals_a[..., 0]
    end_a = intervals_a[..., 1]
    start_b = intervals_b[..., 0]
    end_b = intervals_b[..., 1]

    overlap_start = torch.maximum(start_a, start_b)
    overlap_end = torch.minimum(end_a, end_b)

    intervals_overlap = torch.stack([overlap_start, overlap_end], dim=-1)

    valid_dims = overlap_end > overlap_start
    mask = torch.all(valid_dims, dim=-1)

    return intervals_overlap, mask


def get_overlap_from_peers(
    rank_peers: torch.Tensor, intervals_a: torch.Tensor, intervals_b: torch.Tensor
) -> List[Dict[str, torch.Tensor | int]]:
    """
    Given a rank table (rank_peers) and two interval tensors, compute the overlapping
    intervals and associated peer ranks.

    Args:
        rank_peers: Tensor of shape matching intervals_a leading dims (...). Contains peer/global ranks
            (e.g., typically from device_mesh.mesh), but no device_mesh semantics are assumed here.
        intervals_a: Tensor of shape (..., n_dim, 2) where last dim is [start, end).
        intervals_b: Tensor broadcastable to intervals_a.shape (..., n_dim, 2).

    Returns:
        List[dict]: Each dict contains:
            - "peer": int global rank
            - "interval": Tensor of shape (n_dim, 2) with [start, end) along each axis.

    Notes:
        - Leading dims of intervals_a/intervals_b must match rank_peers.shape so that boolean
          masking (rank_peers[valid_mask]) is valid.
        - Both interval tensors must use trailing (n_dim, 2) layout with exclusive ends.
    """
    if intervals_a.shape[-1] != 2:
        raise ValueError(f"intervals_a last dim must be 2 (start,end), got {intervals_a.shape}")
    n_dim = intervals_a.shape[-2]
    if intervals_b.shape[-1] != 2 or intervals_b.shape[-2] != n_dim:
        raise ValueError(
            f"intervals_b trailing shape must be ({n_dim}, 2) to match intervals_a; got {intervals_b.shape[-2:]}"
        )

    # Broadcast intervals_b to intervals_a shape
    try:
        intervals_b_broadcast = torch.broadcast_to(intervals_b, intervals_a.shape)
    except RuntimeError as e:
        raise ValueError(f"intervals_b is not broadcastable to intervals_a shape {intervals_a.shape}") from e

    # Validate submesh shape matches leading dims
    leading_shape = intervals_a.shape[:-2]
    if rank_peers.shape != leading_shape:
        raise ValueError(f"rank_peers.shape {rank_peers.shape} must match intervals leading shape {leading_shape}")

    intervals_a_view = intervals_a
    intervals_b_view = intervals_b_broadcast

    overlap_intervals, valid_mask = compute_interval_overlap(intervals_a_view, intervals_b_view)
    if not valid_mask.any().item():
        return []

    if valid_mask.dim() != rank_peers.dim():
        raise ValueError(f"valid_mask dims {valid_mask.dim()} must match rank_peers.dim {rank_peers.dim()}")

    peer_ranks = rank_peers[valid_mask]
    interval_selected = overlap_intervals[valid_mask]

    needed = []
    for idx in range(len(peer_ranks)):
        interval = interval_selected[idx]
        needed.append({"peer": peer_ranks[idx].item(), "interval": interval})
    return needed


class DistributedOuterGather(torch.autograd.Function):
    @staticmethod
    def forward(
        ctx,
        z_dtensor: DTensor,
        idx_n_dtensor: DTensor,
        idx_m_dtensor: DTensor,
        axis: int,
        are_ids_contiguous: bool,
        idx_n_mask: DTensor | None = None,
        idx_m_mask: DTensor | None = None,
    ) -> DTensor:
        """Forward pass for DistributedOuterGather.

        "Outer" refers to gathering a Cartesian product block: for each pair of
        index sets (n, m) you gather the rectangular sub-block z[..., n, m, ...],
        producing an output of shape (..., K, H, *features) where H is the width
        of the m index set and K is the shared leading index count.

        Args:
            z_dtensor (DTensor): Tensor of shape ``(*batch, N, M, *features)``. The
                ``axis`` argument identifies the start of the ``(N, M)`` block in
                ``z_dtensor`` (i.e., ``z_dtensor.shape[axis] == N`` and
                ``z_dtensor.shape[axis + 1] == M``).
            idx_n_dtensor (DTensor): Tensor of shape ``(*batch, K, W)`` containing
                gather indices into the ``N`` dimension. Must be sharded on its
                ``-2`` dimension across one of the two mesh dims that shard
                ``z_dtensor``'s ``axis``/``axis+1``.
            idx_m_dtensor (DTensor): Tensor of shape ``(*batch, K, H)`` containing
                gather indices into the ``M`` dimension. Must share device mesh and
                placements with ``idx_n_dtensor``.
            axis (int): Index in ``z_dtensor`` where the ``(N, M)`` block begins.
            are_ids_contiguous (bool, optional): This is a heuristic for selecting the underlying
                send/recv strategy for performance purpose. Currently only True is supported,
                which means that the idx_n and idx_m tensors map to a contiguous block of z
                along (axis, axis + 1) dimensions for all the shards and for all the leading
                (batch) dimensions. When True, the underlying strategy will use the min/max
                of idx_n and idx_m to compute the needed interval, assuming the resulting
                buffer to be communicated across the ranks is fully (or approximately so)
                utilized. In the case the data inside idx_n and idx_m doesn't mapped to
                contiguous blocks, the result will still be correct by setting
                are_ids_contiguous=True but the buffer of z chunks communicated will contain
                a lot of unused elements, making the sed/recv inefficient.
            idx_n_mask (DTensor, optional): Mask for idx_n_dtensor of shape ``(*batch, K, W)``.
                Same device_mesh and placements as idx_n_dtensor. True = valid, False = invalid.
            idx_m_mask (DTensor, optional): Mask for idx_m_dtensor of shape ``(*batch, K, H)``.
                Same device_mesh and placements as idx_m_dtensor. True = valid, False = invalid.


        Device mesh / placements:
            - ``z_dtensor`` must be sharded along both ``axis`` and ``axis+1``.
            - No ``Partial`` placements are allowed.
            - ``idx_n_dtensor``/``idx_m_dtensor`` must share device mesh and
              placements. Their mesh can be either the same as ``z_dtensor`` or the
              flatten of ``z_dtensor``'s mesh over the two sharding axes
              ``(mesh_dim_axis, mesh_dim_axis_plus_1)``.
            - All co-sharding/co-replication checks enforced in the function apply.

        Returns:
            DTensor: Output of shape ``(*batch, K, H, *features)`` with the same
            device mesh and placements as ``idx_n_dtensor``/``idx_m_dtensor``.
        """

        if not are_ids_contiguous:
            raise NotImplementedError(
                "DistributedOuterGather currently only supports are_ids_contiguous=True. "
                "The current implementation is not efficient when the ids are not contiguous "
                "and especially so if they mapped to distal parts of z, "
                "even though the result will still be correct."
            )

        # 1. Validations
        if (
            not isinstance(z_dtensor, DTensor)
            or not isinstance(idx_n_dtensor, DTensor)
            or not isinstance(idx_m_dtensor, DTensor)
        ):
            raise TypeError("All inputs must be DTensors")

        batch_dims_z = z_dtensor.shape[:axis]

        if batch_dims_z != idx_n_dtensor.shape[:-2]:
            raise ValueError(
                f"Batch dimensions of z ({batch_dims_z}) must match batch dimensions of idx_n ({idx_n_dtensor.shape[:-2]})"
            )

        if batch_dims_z != idx_m_dtensor.shape[:-2]:
            raise ValueError(
                f"Batch dimensions of z ({batch_dims_z}) must match batch dimensions of idx_m ({idx_m_dtensor.shape[:-2]})"
            )

        if idx_n_dtensor.shape[-2] != idx_m_dtensor.shape[-2]:
            raise ValueError(
                f"Number of n indices ({idx_n_dtensor.shape[-2]}) must match number of m indices ({idx_m_dtensor.shape[-2]})"
            )

        mesh = z_dtensor.device_mesh
        idx_mesh = idx_n_dtensor.device_mesh

        # idx_n and idx_m must share device_mesh/placements with each other
        if idx_m_dtensor.device_mesh != idx_mesh:
            raise ValueError("idx_n and idx_m must be on the same DeviceMesh")

        if mesh.ndim < idx_mesh.ndim:
            raise ValueError(
                f"z_dtensor.device_mesh.ndim must be no less than"
                f"idx_n/m_dtensor.device_mesh.ndim but got: z_dtensor: {mesh.ndim}, idx_n: {idx_mesh.ndim}, idx_m: {idx_mesh.ndim}"
            )

        if idx_n_dtensor.placements != idx_m_dtensor.placements:
            raise ValueError("idx_n and idx_m must have the same placements")

        ndim_z = z_dtensor.ndim
        if axis < 0:
            axis += ndim_z
        if axis < 0 or axis + 1 >= ndim_z:
            raise ValueError(f"axis {axis} must satisfy 0 <= axis and axis+1 < z.ndim (got z.ndim={ndim_z})")

        z_placements = z_dtensor.placements
        idx_n_placements = idx_n_dtensor.placements

        # Identify sharding dims
        mesh_dim_axis = None
        mesh_dim_axis_plus_1 = None
        mesh_dim_shard_k = None
        for i_mesh_dim, p_z in enumerate(z_placements):
            check_idx_placements = i_mesh_dim < len(idx_n_placements)
            if isinstance(p_z, Partial) or (check_idx_placements and isinstance(idx_n_placements[i_mesh_dim], Partial)):
                raise ValueError("Partial placements are not supported")
            elif isinstance(p_z, Shard):
                if p_z.dim == axis:
                    mesh_dim_axis = i_mesh_dim
                elif p_z.dim == axis + 1:
                    mesh_dim_axis_plus_1 = i_mesh_dim
                else:
                    # "co-sharding" requirement: we require the device_mesh axis
                    # that shards any z_dtensor axis other
                    # than (axis, axis + 1) to also shard the same tensor axis
                    # in idx_n_dtensor and idx_m_dtensor because the following p2p
                    # communication ops are restricted within the sub-device_mesh
                    # of (mesh_dim_axis, mesh_dim_axis_plus_1)
                    if check_idx_placements and idx_n_placements[i_mesh_dim].dim != p_z.dim:
                        raise ValueError(
                            f"z_dtensor's sharded axis {p_z.dim} is outside of {(axis, axis + 1)} but "
                            f"the same tensor axis in idx_n_dtensor and idx_m_dtensor is not sharded "
                            f"by the same device mesh axis {i_mesh_dim}"
                        )
                # any Shard placement must evenly shard any axis
                if z_dtensor.shape[p_z.dim] % mesh.size(i_mesh_dim) != 0:
                    raise ValueError(
                        f"z_dtensor's sharded axis {p_z.dim} of size {z_dtensor.shape[p_z.dim]} is not "
                        f"evenly divisible by mesh dim {i_mesh_dim} (size {mesh.size(i_mesh_dim)})"
                    )
            elif isinstance(p_z, Replicate):
                if check_idx_placements and not isinstance(idx_n_placements[i_mesh_dim], Replicate):
                    # "co-replicating" requirement: for the same reason as co-sharding,
                    # we require all orthogonal Replicate placements to be consistent across
                    # the 3 DTensors to limit the number of P2P communication ops in the
                    # (mesh_dim_axis, mesh_dim_axis_plus_1) sub-device_mesh
                    raise ValueError(
                        f"z_dtensor's replicate placement at mesh axis {i_mesh_dim} is expected to "
                        f"be Replicate placement in idx_n_dtensor but the latter's placement is {idx_n_placements[i_mesh_dim]}"
                    )
            if check_idx_placements and isinstance(idx_n_placements[i_mesh_dim], Shard):
                i_axis_idx_n = idx_n_placements[i_mesh_dim].dim
                i_axis_idx_n = i_axis_idx_n if i_axis_idx_n >= 0 else i_axis_idx_n + idx_n_dtensor.ndim
                if i_axis_idx_n == idx_n_dtensor.ndim - 2:
                    mesh_dim_shard_k = i_mesh_dim
                else:
                    # Due to the same "co-sharding" requirement as above, any other sharding mesh axis must
                    # shard the same tensor axis between z_dtensor and idx_n_dtensor and idx_m_dtensor
                    if not (
                        isinstance(z_placements[i_mesh_dim], Shard) and z_placements[i_mesh_dim].dim == i_axis_idx_n
                    ):
                        raise ValueError(
                            f"idx_n_dtensor and idx_m_dtensor's sharded axis {i_axis_idx_n} by mesh axis {i_mesh_dim} "
                            f"is not a co-sharded axis: got z_dtensor's placement at mesh axis {i_mesh_dim} as {z_placements[i_mesh_dim]}"
                        )
                # any Shard placement must evenly shard any axis
                if idx_n_dtensor.shape[i_axis_idx_n] % mesh.size(i_mesh_dim) != 0:
                    raise ValueError(
                        f"idx_n_dtensor's sharded axis {i_axis_idx_n} of size {idx_n_dtensor.shape[i_axis_idx_n]} is not "
                        f"evenly divisible by mesh dim {i_mesh_dim} (size {mesh.size(i_mesh_dim)})"
                    )

        if mesh_dim_axis is None or mesh_dim_axis_plus_1 is None:
            raise ValueError(f"z must be sharded along axis {axis} and {axis + 1}")

        if mesh_dim_shard_k is None:
            raise ValueError("idx_n_dtensor must be sharded along axis -2")

        # Allow idx_n/idx_m to reside on a flattened mesh over (mesh_dim_axis, mesh_dim_axis_plus_1)
        mesh_tensor = mesh.mesh
        idx_mesh_tensor = idx_mesh.mesh
        mesh_flat_tensor = torch.flatten(mesh_tensor, start_dim=mesh_dim_axis, end_dim=mesh_dim_axis_plus_1)
        mesh_compatible = torch.equal(idx_mesh_tensor, mesh_tensor) or torch.equal(idx_mesh_tensor, mesh_flat_tensor)
        if not mesh_compatible:
            raise ValueError(
                "idx_n/idx_m device_mesh must match z device_mesh or its flatten over (mesh_dim_axis, mesh_dim_axis_plus_1)"
            )

        if mesh_dim_shard_k not in (mesh_dim_axis, mesh_dim_axis_plus_1):
            # one of (mesh_dim_axis, mesh_dim_axis_plus_1) must shard idx_n_dtensor and idx_m_dtensor
            # along their axis -2
            raise ValueError(
                f"mesh_dim_shard_k {mesh_dim_shard_k} must be one of (mesh_dim_axis, mesh_dim_axis_plus_1) "
                f"but got {mesh_dim_shard_k}"
            )

        # Validate idx_n_mask if provided
        if idx_n_mask is not None:
            if not isinstance(idx_n_mask, DTensor):
                raise TypeError("idx_n_mask must be a DTensor")
            if idx_n_mask.shape != idx_n_dtensor.shape:
                raise ValueError(
                    f"idx_n_mask shape {idx_n_mask.shape} must match idx_n_dtensor shape {idx_n_dtensor.shape}"
                )
            if idx_n_mask.device_mesh != idx_n_dtensor.device_mesh:
                raise ValueError("idx_n_mask must have the same device_mesh as idx_n_dtensor")
            if idx_n_mask.placements != idx_n_dtensor.placements:
                raise ValueError("idx_n_mask must have the same placements as idx_n_dtensor")
            if idx_n_mask.dtype != torch.bool:
                raise TypeError(
                    f"idx_n_mask must have dtype torch.bool, got {idx_n_mask.dtype}. Use mask.bool() to convert."
                )

        # Validate idx_m_mask if provided
        if idx_m_mask is not None:
            if not isinstance(idx_m_mask, DTensor):
                raise TypeError("idx_m_mask must be a DTensor")
            if idx_m_mask.shape != idx_m_dtensor.shape:
                raise ValueError(
                    f"idx_m_mask shape {idx_m_mask.shape} must match idx_m_dtensor shape {idx_m_dtensor.shape}"
                )
            if idx_m_mask.device_mesh != idx_m_dtensor.device_mesh:
                raise ValueError("idx_m_mask must have the same device_mesh as idx_m_dtensor")
            if idx_m_mask.placements != idx_m_dtensor.placements:
                raise ValueError("idx_m_mask must have the same placements as idx_m_dtensor")
            if idx_m_mask.dtype != torch.bool:
                raise TypeError(
                    f"idx_m_mask must have dtype torch.bool, got {idx_m_mask.dtype}. Use mask.bool() to convert."
                )

        if z_dtensor.device_mesh.ndim == idx_n_dtensor.device_mesh.ndim:
            # When idx_n_dtensor.device_mesh.ndim == z_dtensor.device_mesh.ndim,
            # the one of idx_n_dtensor.placements along (mesh_dim_axis, mesh_dim_axis_plus_1)
            # must be Replicate() and the other must be Shard(-2) because of the "co-sharding"
            # and "co-replicating" requirements and the fact that we exclude Partial placements.
            if mesh_dim_shard_k == mesh_dim_axis:
                mesh_dim_replicate_idx = mesh_dim_axis_plus_1
            else:
                mesh_dim_replicate_idx = mesh_dim_axis
            if not isinstance(idx_n_placements[mesh_dim_replicate_idx], Replicate):
                raise ValueError(
                    f"idx_n_dtensor's Replicate placement at mesh axis {mesh_dim_replicate_idx} is expected to "
                    f"be Replicate but the latter's placement is {idx_n_placements[mesh_dim_replicate_idx]}"
                )
        else:
            # NOTE: there is however an exception where idx_n_dtensor.device_mesh is a flattened
            # mesh of z_dtensor along (mesh_dim_axis, mesh_dim_axis_plus_1) then the would-be Replicate()
            # devices participate in sharding along idx_n_dtensor's axis -2 but then the two meshes
            # ndim are different.
            mesh_dim_replicate_idx = None

        # Assert even sharding as per assumption
        if mesh.size(mesh_dim_axis) != mesh.size(mesh_dim_axis_plus_1):
            # TODO: the p2p comm doesn't actually require this
            raise ValueError("Mesh dimensions for sharding z must have equal size")

        # 2. Local Setup
        z_local = z_dtensor.to_local()
        idx_n_local = idx_n_dtensor.to_local()
        idx_m_local = idx_m_dtensor.to_local()
        idx_n_mask_local = idx_n_mask.to_local() if idx_n_mask is not None else None
        idx_m_mask_local = idx_m_mask.to_local() if idx_m_mask is not None else None

        device = z_local.device
        cpu_device = torch.device("cpu")

        # 3. Compute Local Bounding Box (Needed Ranges)
        # TODO: min/max should not take across batch dims
        # The current approach assumes the max z buffer size
        # to be send/recv across the leading batch dimensions,
        # e.g., if:
        #     idx_n_dtensor[0, :, :] and idx_m_dtensor[0, :, :]
        # define/need a z buffer bounded by:
        #     z[0, min_0:max_0, min_1:max_1, ...],
        # while:
        #     idx_n_dtensor[1, :, :] and idx_m_dtensor[1, :, :]
        # need a z buffer bounded by:
        #     z[1, (min_0-delta0):(max_0+delta0), (min_1-delta1):(max_1+delta1), ...],
        # where delta0 > 0 and delta1 > 0 are True, then the actual z buffer
        # will be of shape:
        #     z[:, (min_0-delta0):(max_0+delta0), (min_1-delta1):(max_1+delta1), ...],
        # i.e., the maximal z_buffer is assumed across all leading batch dimensions.
        # This still makes valid computation because the local OuterGather will still
        # correctly figure out the indices from idx_n_local and idx_m_local, but
        # we actually send/recv more data than necessary, making communication
        # suboptimal in general.
        # NOTE: however, in our Boltz application, idx_n_dtensor and idx_m_dtensor are
        # exclusively the atom to token indices, which are homogeneous within a sample
        # across the multiplicity dimension inside AtomAttnEncoder/Decoder so the above "delta"
        # is always 0 as soon as we don't actually have sample batching in the leading dimensions
        # NOTE: supporting more efficient communication with minimal z buffer for send/recv
        # requires a 3D version of the "get_flattened_range_indices()" from utils.py
        # where we need to deal with (B, N, M) indices instead of (D, 2) indices, which
        # actually incur higher memory overhead because then we requires 3 * n_elements_in_z_buf
        # instead of 2 * n_elements_in_z_buf (because of the extra dimension)
        # Build needed interval tensor shape (2, 2): dim0 -> n axis, dim1 -> m axis, last dim is [start, end)
        if idx_n_local.numel() > 0:
            if idx_n_mask_local is not None:
                # Only consider valid indices for interval computation
                if idx_n_mask_local.any():
                    valid_idx_n = idx_n_local[idx_n_mask_local]
                    need_interval_n = torch.stack(valid_idx_n.aminmax())
                else:
                    # All indices are masked out
                    need_interval_n = torch.tensor([0, -1], device=device, dtype=idx_n_local.dtype)
            else:
                # inclusive interval of shape (2,)
                need_interval_n = torch.stack(idx_n_local.aminmax())
        else:
            need_interval_n = torch.tensor([0, -1], device=device, dtype=idx_n_local.dtype)

        if idx_m_local.numel() > 0:
            if idx_m_mask_local is not None:
                # Only consider valid indices for interval computation
                if idx_m_mask_local.any():
                    valid_idx_m = idx_m_local[idx_m_mask_local]
                    need_interval_m = torch.stack(valid_idx_m.aminmax())
                else:
                    # All indices are masked out
                    need_interval_m = torch.tensor([0, -1], device=device, dtype=idx_m_local.dtype)
            else:
                # inclusive interval of shape (2,)
                need_interval_m = torch.stack(idx_m_local.aminmax())
        else:
            need_interval_m = torch.tensor([0, -1], device=device, dtype=idx_m_local.dtype)

        # inclusive interval of both (axis, axis + 1) of shape (2, 2)
        need_interval = torch.stack([need_interval_n, need_interval_m]).to(dtype=torch.long)
        # left-inclusive and right-exclusive interval of shape (2, 2)
        need_interval[:, -1] += 1
        need_start_vec = need_interval[:, 0]
        need_end_vec = need_interval[:, 1]
        # need_interval must be on same device as the device_mesh for later all_gather but
        # save one cpu copy for local indexing computation
        need_start_vec_cpu = need_start_vec.to(cpu_device)

        # 4. Exchange Metadata

        # Compute my owned z range as interval tensor shape (2, 2)
        my_coord_vec = torch.tensor(
            [mesh.get_local_rank(mesh_dim_axis), mesh.get_local_rank(mesh_dim_axis_plus_1)], device=cpu_device
        )
        chunk_sizes = torch.tensor(
            [
                z_dtensor.shape[axis] // mesh.size(mesh_dim_axis),
                z_dtensor.shape[axis + 1] // mesh.size(mesh_dim_axis_plus_1),
            ],
            device=cpu_device,
        )

        own_start = my_coord_vec * chunk_sizes
        own_end = own_start + chunk_sizes
        own_interval = torch.stack([own_start, own_end], dim=-1)  # (2, 2)

        # Gather along axis+1 (M) first, then along axis (N) to keep N leading
        # TODO: due to the limitation of torch DeviceMesh not being able to retrieve
        # group of its submesh, we need to manually all_gather along each mesh axis.
        # An alternative is to create the submesh an flatten it then retrieve the group
        # but the overhead of creating submesh and new groups will recur upon all invocations
        # of this function and the resulting groups and submesh will not be managed, leading to
        # waste of distributed resources.
        group_m = mesh.get_group(mesh_dim_axis_plus_1)
        need_range_m = [torch.zeros_like(need_interval) for _ in range(mesh.size(mesh_dim_axis_plus_1))]
        dist.all_gather(need_range_m, need_interval, group=group_m)
        need_range_m = torch.stack(need_range_m)  # (Grid_M, 2, 2)

        group_n = mesh.get_group(mesh_dim_axis)
        need_range_nm = [torch.zeros_like(need_range_m) for _ in range(mesh.size(mesh_dim_axis))]
        dist.all_gather(need_range_nm, need_range_m, group=group_n)
        need_range_nm = torch.stack(need_range_nm)  # (Grid_N, Grid_M, 2, 2) on device
        need_range_nm_cpu = need_range_nm.cpu()

        ranks_global_on_mesh = mesh.mesh

        # Get my coords in the mesh
        my_coords = mesh.get_coordinate()

        # Slice out the 2D submesh over (mesh_dim_axis, mesh_dim_axis_plus_1) anchored at my_coords on other dims.
        index_list_submesh = []
        for dim in range(ranks_global_on_mesh.ndim):
            if dim == mesh_dim_axis or dim == mesh_dim_axis_plus_1:
                index_list_submesh.append(slice(None))
            else:
                index_list_submesh.append(torch.tensor(my_coords[dim], device=cpu_device))
        ranks_global_on_submesh = ranks_global_on_mesh[tuple(index_list_submesh)]  # shape (size_group_n, size_group_m)

        # 5. P2P Logic
        ops = []
        recv_bufs = {}  # key: peer_rank, value: buffer
        recv_metadata_for_bwd = []

        # Mesh dimensions
        size_group_n = mesh.size(mesh_dim_axis)
        size_group_m = mesh.size(mesh_dim_axis_plus_1)

        # --- RECEIVE PLAN ---
        if need_start_vec[0] >= need_end_vec[0] or need_start_vec[1] >= need_end_vec[1]:
            needed_chunks = []
        else:
            # We need chunks that overlap with need_interval
            # Construct bounds for all chunks in the mesh
            # chunk_i goes 0..size_group_n-1, chunk_j goes 0..size_group_m-1
            # We use meshgrid to generate indices for all chunks
            # Note: This is similar to the Send Plan logic where we flatten the mesh view

            i_ranks_n = torch.arange(size_group_n, device=cpu_device)
            i_ranks_m = torch.arange(size_group_m, device=cpu_device)
            grid_n, grid_m = torch.meshgrid(i_ranks_n, i_ranks_m, indexing="ij")
            grid_coords = torch.stack([grid_n, grid_m], dim=-1)  # (size_group_n, size_group_m, 2)

            # Calculate bounds for each chunk
            starts = grid_coords * chunk_sizes  # (size_group_n, size_group_m, 2)
            ends = starts + chunk_sizes  # (size_group_n, size_group_m, 2)
            peers_own_intervals = torch.stack([starts, ends], dim=-1)  # (size_group_n, size_group_m, 2, 2)

            # Needed bounds (broadcasted)
            # need_start_vec / need_end_vec are scalars (exclusive end) per axis

            # Compute overlaps between peers' owned intervals and the current rank's needed intervals
            need_interval_cpu = need_interval.to(cpu_device)
            # needed_chunks: list of dicts, each dict contains:
            # - "peer": int global rank
            # - "interval": Tensor of shape (2, 2) with [i_start, i_end) in the last axis indicating the
            #   interval of z along its (axis, axis + 1) dimensions.
            needed_chunks = get_overlap_from_peers(ranks_global_on_submesh, peers_own_intervals, need_interval_cpu)

        for item in needed_chunks:
            peer = item["peer"]
            interval = item["interval"]
            starts_global = interval[:, 0]
            lens = interval[:, 1] - interval[:, 0]

            shape = list(z_local.shape)
            shape[axis : axis + 2] = lens.tolist()

            buf = torch.empty(shape, dtype=z_local.dtype, device=device)

            if peer == dist.get_rank():
                # Self-copy
                starts_local = starts_global - own_start
                chunk = z_local.narrow(axis, starts_local[0].item(), lens[0].item()).narrow(
                    axis + 1, starts_local[1].item(), lens[1].item()
                )
                buf.copy_(chunk)
                recv_bufs[peer] = buf
                recv_metadata_for_bwd.append((peer, interval, shape))
            else:
                ops.append(dist.P2POp(dist.irecv, buf, peer))
                recv_bufs[peer] = buf
                recv_metadata_for_bwd.append((peer, interval, shape))

        # --- SEND PLAN ---
        send_metadata_for_bwd = []

        # need_range_nm: (Grid_N, Grid_M, 4) -> (size_group_n, size_group_m, 4)
        # Assuming need_range_nm[i, j] corresponds to peer at coords (i, j) relative to (mesh_dim_axis, mesh_dim_axis_plus_1)
        # Flatten need_range_nm to process all at once

        send_chunks = get_overlap_from_peers(
            ranks_global_on_submesh,
            need_range_nm_cpu.view(size_group_n, size_group_m, 2, 2),
            own_interval.view(1, 1, 2, 2),
        )

        my_rank = dist.get_rank()

        for item in send_chunks:
            peer_rank = item["peer"]
            if peer_rank == my_rank:
                continue

            interval = item["interval"]
            starts_global = interval[:, 0]
            lens = interval[:, 1] - interval[:, 0]

            starts_local = starts_global - own_start

            chunk = z_local.narrow(axis, starts_local[0].item(), lens[0].item()).narrow(
                axis + 1, starts_local[1].item(), lens[1].item()
            )
            chunk = chunk.contiguous()
            ops.append(dist.P2POp(dist.isend, chunk, peer_rank))

            send_metadata_for_bwd.append(
                (
                    peer_rank,
                    starts_local.to(device=cpu_device, dtype=torch.long),
                    lens.to(device=cpu_device, dtype=torch.long),
                )
            )

        # Execute P2P
        if ops:
            reqs = dist.batch_isend_irecv(ops)
            for req in reqs:
                req.wait()

        # 6. Local Computation
        if need_start_vec[0] >= need_end_vec[0] or need_start_vec[1] >= need_end_vec[1]:
            buffer_shape = list(z_local.shape)
            buffer_shape[axis] = 0
            buffer_shape[axis + 1] = 0
            z_buffer = torch.empty(buffer_shape, dtype=z_local.dtype, device=device)
        else:
            buffer_shape = list(z_local.shape)
            buffer_shape[axis] = (need_end_vec[0] - need_start_vec[0]).item()
            buffer_shape[axis + 1] = (need_end_vec[1] - need_start_vec[1]).item()
            z_buffer = torch.zeros(buffer_shape, dtype=z_local.dtype, device=device)

            for item in needed_chunks:
                peer = item["peer"]
                interval = item["interval"]
                if peer not in recv_bufs:
                    raise RuntimeError(
                        f"Peer {peer} not in recv_bufs, which should not happen because "
                        f"the recv planning should guarantee all peer ranks' contribution "
                        f"to the recv buffer arrive by now"
                    )
                buf = recv_bufs[peer]
                starts = interval[:, 0]
                lens = interval[:, 1] - interval[:, 0]

                starts_local_z_buffer = starts - need_start_vec_cpu

                target = z_buffer.narrow(axis, starts_local_z_buffer[0].item(), lens[0].item()).narrow(
                    axis + 1, starts_local_z_buffer[1].item(), lens[1].item()
                )
                target.copy_(buf)

        local_idx_n = idx_n_local - need_start_vec[0].item()
        local_idx_m = idx_m_local - need_start_vec[1].item()

        # Use OuterGather instead of RectangularOuterGather
        out_local = OuterGather.apply(z_buffer, local_idx_n, local_idx_m, axis, idx_n_mask_local, idx_m_mask_local)

        # 7. Output DTensor
        out_global_shape = list(idx_n_dtensor.shape)
        H = idx_m_dtensor.shape[-1]
        feature_shape = z_dtensor.shape[axis + 2 :]

        final_global_shape = list(out_global_shape) + [H] + list(feature_shape)
        strides_out = update_exhaustive_strides(out_local.shape, out_local.stride(), tuple(final_global_shape))

        out_dtensor = DTensor.from_local(
            out_local, idx_n_dtensor.device_mesh, idx_n_placements, shape=tuple(final_global_shape), stride=strides_out
        )

        tensors_to_save = [idx_n_local, idx_m_local]
        if idx_n_mask_local is not None:
            tensors_to_save.append(idx_n_mask_local)
        if idx_m_mask_local is not None:
            tensors_to_save.append(idx_m_mask_local)
        ctx.save_for_backward(*tensors_to_save)
        ctx.has_n_mask = idx_n_mask_local is not None
        ctx.has_m_mask = idx_m_mask_local is not None
        ctx.comm_meta = {
            "recv_metadata_for_bwd": recv_metadata_for_bwd,
            "send_metadata_for_bwd": send_metadata_for_bwd,
            "z_local_shape": z_local.shape,
            "z_buffer_shape": z_buffer.shape,
            "need_interval": need_interval,
            "axis": axis,
            "device_mesh_z": z_dtensor.device_mesh,
            "z_placements": z_placements,
            "z_global_shape": z_dtensor.shape,
            "output_placements": out_dtensor.placements,
            "own_interval": own_interval,
            "device_mesh_output": out_dtensor.device_mesh,
            "mesh_dim_replicate_idx": mesh_dim_replicate_idx,
        }

        return out_dtensor

    @staticmethod
    def backward(ctx, grad_output: DTensor) -> Tuple[DTensor, None, None, None, None, None, None]:
        """Backward pass for DistributedOuterGather.

        Uses forward-phase send/recv metadata (saved in ctx) to drive the P2P
        communication plan that redistributes gradient contributions back to the
        owning shards of ``z_dtensor``.

        Args:
            grad_output (DTensor): Gradient of the output with shape
                ``(*batch, K, H, *features)`` and the same device mesh/placements
                as the forward output.

        Returns:
            Tuple[DTensor, None, None, None, None, None, None]: Gradient for ``z_dtensor`` and
            ``None`` for non-differentiable inputs.
        """
        saved = ctx.saved_tensors
        idx_n_local = saved[0]
        idx_m_local = saved[1]
        idx = 2
        idx_n_mask_local = saved[idx] if ctx.has_n_mask else None
        if ctx.has_n_mask:
            idx += 1
        idx_m_mask_local = saved[idx] if ctx.has_m_mask else None
        meta = ctx.comm_meta
        recv_meta = meta["recv_metadata_for_bwd"]
        send_meta = meta["send_metadata_for_bwd"]
        z_local_shape = meta["z_local_shape"]
        z_buffer_shape = meta["z_buffer_shape"]
        need_interval = meta["need_interval"]
        need_start_vec = need_interval[:, 0]
        axis = meta["axis"]
        device_mesh_z = meta["device_mesh_z"]
        z_placements = meta["z_placements"]
        output_placements = meta["output_placements"]
        z_global_shape = meta["z_global_shape"]
        own_interval = meta["own_interval"]
        device_mesh_output = meta["device_mesh_output"]
        mesh_dim_replicate_idx = meta["mesh_dim_replicate_idx"]

        if device_mesh_output != grad_output.device_mesh:
            raise ValueError(
                f"grad_output device_mesh mismatch: expected same device_mesh from fwd's output: {device_mesh_output}, "
                f"but got {grad_output.device_mesh}"
            )

        # TODO: we can actually support grad_output.placements[mesh_dim_replicate_idx] == Partial("sum")
        # where then the latter division of grad_z_local by device_mesh.size(mesh_dim_replicate_idx) is
        # not more necessary because the grad_output virtually would perform the all-reduce along that mesh axis
        if output_placements != grad_output.placements:
            raise ValueError(
                f"grad_output placements mismatch: expected same placements from fwd's output: {output_placements}, "
                f"but got {grad_output.placements}"
            )

        grad_local = grad_output.to_local()
        # Ensure contiguous for backward ops
        grad_local = grad_local.contiguous()

        local_idx_n = idx_n_local - need_start_vec[0].item()
        local_idx_m = idx_m_local - need_start_vec[1].item()

        # Use OuterGather.backward logic (actually uses outer_gather_backward helper)

        grad_z_buffer = outer_gather_backward(
            grad_local, z_buffer_shape, local_idx_n, local_idx_m, axis, idx_n_mask_local, idx_m_mask_local
        )

        ops = []
        grad_z_local = torch.zeros(z_local_shape, dtype=grad_local.dtype, device=grad_local.device)

        # 1. Backward Send (Reverse of Fwd Recv)
        for peer, interval, shape in recv_meta:
            # grad_z_buffer corresponds to the z_buffer in the forward pass
            # so we need to convert the interval to the local indices of grad_z_buffer
            # by subtracting the need_start_vec, which is the global indices of the z_buffer
            need_start_vec_cpu = need_start_vec.to(interval.device)
            starts_local_z_buffer = interval[:, 0] - need_start_vec_cpu
            lens = interval[:, 1] - interval[:, 0]

            grad_chunk = grad_z_buffer.narrow(axis, starts_local_z_buffer[0].item(), lens[0].item()).narrow(
                axis + 1, starts_local_z_buffer[1].item(), lens[1].item()
            )

            if peer == dist.get_rank():
                # Self-accumulate into grad_z_local
                # Here we need to the offset the z_local as in the forward pass input
                # to copy upstream adjoints into z_local's gradients
                starts_local = interval[:, 0] - own_interval[:, 0]

                target = grad_z_local.narrow(axis, starts_local[0].item(), lens[0].item()).narrow(
                    axis + 1, starts_local[1].item(), lens[1].item()
                )
                target.add_(grad_chunk)
            else:
                ops.append(dist.P2POp(dist.isend, grad_chunk.contiguous(), peer))

        # 2. Backward Recv (Reverse of Fwd Send)
        bwd_recv_bufs = []

        for peer, starts_local, lens in send_meta:
            # send_meta only contains peers != self
            shape = list(z_local_shape)
            shape[axis : axis + 2] = lens.tolist()

            buf = torch.empty(shape, dtype=grad_local.dtype, device=grad_local.device)
            ops.append(dist.P2POp(dist.irecv, buf, peer))
            bwd_recv_bufs.append((buf, starts_local, lens))

        if ops:
            reqs = dist.batch_isend_irecv(ops)
            for req in reqs:
                req.wait()

        for buf, starts_local, lens in bwd_recv_bufs:
            target = grad_z_local.narrow(axis, starts_local[0].item(), lens[0].item()).narrow(
                axis + 1, starts_local[1].item(), lens[1].item()
            )
            target.add_(buf)

        if mesh_dim_replicate_idx is not None:
            # When grad_output.placements[mesh_dim_replicate_idx] == Replicate(),
            # the grad_output contribution among that mesh dimensions are duplicated
            # so we need to divide by the number of devices in that mesh dimension to get the
            # correct gradient.
            grad_z_local = grad_z_local / device_mesh_output.size(mesh_dim_replicate_idx)

        grad_z_dtensor = DTensor.from_local(
            grad_z_local,
            device_mesh_z,
            z_placements,
            shape=z_global_shape,
            stride=update_exhaustive_strides(grad_z_local.shape, grad_z_local.stride(), z_global_shape),
        )

        return grad_z_dtensor, None, None, None, None, None, None


def distributed_outer_gather(
    z_dtensor: DTensor,
    idx_n_dtensor: DTensor,
    idx_m_dtensor: DTensor,
    axis: int = 1,
    are_ids_contiguous: bool = False,
    idx_n_mask: DTensor | None = None,
    idx_m_mask: DTensor | None = None,
) -> DTensor:
    """Distributed outer gather convenience wrapper.

    "Outer" means taking the Cartesian product of two index sets: for every pair
    of indices drawn from ``idx_n_dtensor`` (along the ``N`` axis) and
    ``idx_m_dtensor`` (along the ``M`` axis), the corresponding rectangular
    block ``z_dtensor[..., n, m, ...]`` is gathered, producing an output of
    shape ``(*batch, K, H, *features)`` where ``K`` is the shared leading index
    count and ``H`` is the width of the ``m`` index set.

    Args:
        z_dtensor (DTensor): Shape ``(*batch, N, M, *features)``. Must be sharded
            on the ``(N, M)`` block along two mesh dims corresponding to
            ``axis``/``axis+1``; no ``Partial`` placements allowed.
        idx_n_dtensor (DTensor): Shape ``(*batch, K, W)``. Shares device mesh and
            placements with ``idx_m_dtensor``. Must be sharded on ``-2`` over one
            of the two mesh dims that shard ``z_dtensor``'s ``(N, M)`` block.
            Device mesh can be the same as ``z_dtensor`` or the flatten over the
            two sharding mesh dims.
        idx_m_dtensor (DTensor): Shape ``(*batch, K, H)``, same mesh/placements as
            ``idx_n_dtensor``.
        axis (int, optional): Start axis of the ``(N, M)`` block in ``z_dtensor``.
            Defaults to 1.
        are_ids_contiguous (bool, optional): This is a heuristic for selecting the underlying
            send/recv strategy for performance purpose. Currently only True is supported,
            which means that the idx_n and idx_m tensors map to a contiguous block of z
            along (axis, axis + 1) dimensions for all the shards and for all the leading
            (batch) dimensions. When True, the underlying strategy will use the min/max
            of idx_n and idx_m to compute the needed interval, assuming the resulting
            buffer to be communicated across the ranks is fully (or approximately so)
            utilized. In the case the data inside idx_n and idx_m doesn't mapped to
            contiguous blocks, the result will still be correct by setting
            are_ids_contiguous=True but the buffer of z chunks communicated will contain
            a lot of unused elements, making the sed/recv inefficient.
            Defaults to False as a reminder to the user to understand the performance implications.
        idx_n_mask (DTensor, optional): Mask for idx_n_dtensor of shape ``(*batch, K, W)``.
            Same device_mesh and placements as idx_n_dtensor. True = valid, False = invalid.
        idx_m_mask (DTensor, optional): Mask for idx_m_dtensor of shape ``(*batch, K, H)``.
            Same device_mesh and placements as idx_m_dtensor. True = valid, False = invalid.

    Returns:
        DTensor: Gathered output with shape ``(*batch, K, H, *features)``, using
        the mesh/placements of the index tensors.
    """
    return DistributedOuterGather.apply(
        z_dtensor, idx_n_dtensor, idx_m_dtensor, axis, are_ids_contiguous, idx_n_mask, idx_m_mask
    )
