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
from typing import Optional

import torch
import torch.distributed as dist
from torch.distributed.device_mesh import DeviceMesh
from torch.distributed.tensor import DTensor, Partial, Placement, Shard

from boltz.distributed.model.layers.flatten_and_unflatten import (
    shardwise_flatten,
    shardwise_flatten_sharded,
    shardwise_unflatten_sharded,
)
from boltz.distributed.model.modules.utils import validate_window_batching_parameters
from boltz.distributed.utils import update_exhaustive_strides


def get_query_window_key_range(W: int, H: int, K: int, ids_query_window: torch.Tensor) -> torch.Tensor:
    """
    Get the range of half-window indices (j) that query windows attend to.

    Vectorized version that computes ranges for multiple query windows simultaneously.

    Parameters
    ----------
    W : int
        Atoms per query window (must be even)
    H : int
        Keys per query window (must be divisible by W//2)
    K : int
        Total number of query windows (must be > 0)
    ids_query_window : torch.Tensor
        Query window indices of any shape, with values in range [0, K-1]

    Returns
    -------
    torch.Tensor
        Shape [2, *ids_query_window.shape] where:
        - result[0] contains j_min values
        - result[1] contains j_max values

    Raises
    ------
    AssertionError
        If input parameters don't satisfy constraints

    Examples
    --------
    >>> ids = torch.tensor([0, 2, 9])
    >>> ranges = get_query_window_key_range(32, 128, 10, ids)
    >>> ranges.shape
    torch.Size([2, 3])
    >>> ranges[0]  # j_min values
    tensor([0, 1, 15])
    >>> ranges[1]  # j_max values
    tensor([4, 8, 19])

    >>> ids = torch.tensor([[0, 1], [2, 3]])
    >>> ranges = get_query_window_key_range(32, 128, 10, ids)
    >>> ranges.shape
    torch.Size([2, 2, 2])
    """
    # Validate inputs
    assert W > 0 and W % 2 == 0, "W must be positive and even"
    assert H > 0 and H % (W // 2) == 0, "H must be divisible by W//2"
    assert K > 0, "K must be positive"
    assert torch.all(ids_query_window >= 0) and torch.all(
        ids_query_window < K
    ), f"Query window indices must be in range [0, {K - 1}]"

    # Calculate h (half-windows per query window)
    h = H // (W // 2)

    # Calculate range using proven formula (vectorized)
    # Note: j_max >= j_min is guaranteed by the formula (proven mathematically)

    # j_min: Left-clamped to 0 for early windows when (2*i + 1 - h//2) < 0
    #        This occurs when i < (h//2 - 1)/2, typically first few query windows
    j_min = torch.maximum(torch.zeros_like(ids_query_window), 2 * ids_query_window + 1 - h // 2)

    # j_max: Right-clamped to 2*K-1 for late windows when (2*i + h//2) > 2*K-1
    #        This occurs when i > K - 1 - h//4, typically last few query windows
    j_max = torch.minimum(torch.full_like(ids_query_window, 2 * K - 1), 2 * ids_query_window + h // 2)

    # Co-occurrence: Both clamps can occur simultaneously for very small K (K < h//2 + 1)
    #   In such cases, the query window sees all available half-windows [0, 2K-1]
    #   For typical h=8, this happens when K ≤ 4 (rare in practice)

    # Stack j_min and j_max along new dimension at the front
    return torch.stack([j_min, j_max], dim=0)


def gather_sliding_windows_backward(
    grad_output: torch.Tensor,
    window_start_offsets: torch.Tensor,
    window_size: int,
    axis: int,
    input_shape: tuple,
) -> torch.Tensor:
    """
    Backward pass for sliding window gathering operation.

    Computes gradient w.r.t. input given gradient w.r.t. output.

    Parameters
    ----------
    grad_output : torch.Tensor
        Gradient w.r.t. output, shape (..., n_windows, window_size, ...)
    window_start_offsets : torch.Tensor
        Window starting positions used in forward pass, shape (n_windows,)
    window_size : int
        Size of each window (h)
    axis : int
        Axis along which windowing was applied
    input_shape : tuple
        Shape of the original input tensor

    Returns
    -------
    torch.Tensor
        Gradient w.r.t. input, shape matching input_shape

    Notes
    -----
    Uses index_add_ to accumulate overlapping gradients from multiple windows
    that read the same input positions.
    """
    # Validate input types
    if not isinstance(grad_output, torch.Tensor):
        raise TypeError(f"grad_output must be a torch.Tensor, got {type(grad_output)}")
    if not isinstance(window_start_offsets, torch.Tensor):
        raise TypeError(f"window_start_offsets must be a torch.Tensor, got {type(window_start_offsets)}")

    # Validate window_start_offsets shape
    if window_start_offsets.ndim != 1:
        raise ValueError(
            f"window_start_offsets must be 1D, got {window_start_offsets.ndim}D with shape {window_start_offsets.shape}"
        )

    n_windows = window_start_offsets.shape[0]

    # Normalize and validate axis
    ndim_input = len(input_shape)
    if axis < 0:
        axis_normalized = ndim_input + axis
    else:
        axis_normalized = axis

    if not (0 <= axis_normalized < ndim_input):
        raise ValueError(
            f"axis {axis} out of range for input with {ndim_input} dims "
            f"(normalized to {axis_normalized}, valid range [0, {ndim_input}))"
        )

    # Build expected grad_output shape from input_shape
    # Mapping: (..., in_len, ...) -> (..., n_windows, window_size, ...)
    expected_shape = list(input_shape)
    expected_shape[axis_normalized] = n_windows  # Replace in_len with n_windows
    expected_shape.insert(axis_normalized + 1, window_size)  # Insert window_size after n_windows

    # Validate complete grad_output shape
    if grad_output.shape != tuple(expected_shape):
        raise ValueError(
            f"grad_output shape mismatch:\n"
            f"  Expected: {tuple(expected_shape)}\n"
            f"  Got:      {grad_output.shape}\n"
            f"  (Derived from input_shape={input_shape}, n_windows={n_windows}, "
            f"window_size={window_size}, axis={axis_normalized})"
        )

    device = grad_output.device

    # Use normalized axis for subsequent operations
    axis = axis_normalized
    in_len = input_shape[axis]

    # ============================================================
    # Step A: Normalize grad_output shape for Scatter-Add
    # ============================================================
    # Explanation of the Backward Logic:
    #
    # 1. Preparation (Flattening): Because the input can have arbitrary dimensions
    #    (e.g., Batch, Time, Features or Channels, Time, Height, Width), performing
    #    operations on specific dimensions is tricky. In Step A, we permute and
    #    flatten the tensor into a 2D matrix: [n_windows * window_size, Flat_Features].
    #    This standardizes the problem regardless of the input shape.

    # Forward output was: (..., n_windows, window_size, ...)
    # We want to isolate (n_windows, window_size) and flatten the rest.

    # 1. Move 'window_size' (currently at axis+1) back to end
    #    Current: (..., n_windows, window_size, ...)
    #    Result:  (n_windows, window_size, Flattend_Features)
    g = grad_output.moveaxis([axis, axis + 1], [0, 1]).flatten(2, -1)

    # 4. Final Flatten for Index Add:
    #    Current: (n_windows, window_size, Flattend_Features)
    #    Result: (Total_Elements, Flattened_Features)
    #    Total_Elements = n_windows * window_size
    grad_source_flat = g.flatten(0, 1)

    # ============================================================
    # Step B: Generate Target Indices (Where to accumulate gradients)
    # ============================================================
    # 2. Mapping (target_indices): This is the inverse of the forward index_select.
    #
    #     Forward: "For window w, read index i."
    #
    #     Backward: "For window w, gradient i belongs to index i+window_start_offsets[w]."
    #               We generate a full grid of these destination indices.

    # We need to map every element in (n_windows, window_size) back to
    # the padded input vector index.
    # Formula: index[w, i] = window_start_offsets[w] + pad_top + i
    # where the rhs is exactly the indices of the padded input vector
    # padded_vector in the forward pass, i.e., index[w, i] is the fwd
    # output[..., w, i, ...]'s index in the original source padded_vector.

    # Determine padding (same logic as forward)
    min_k = window_start_offsets.min().item()
    pad_top = max(0, -min_k)

    # 1. Create Base Offsets for each window
    #    Shape: (n_windows, 1)
    base_indices = (window_start_offsets + pad_top).unsqueeze(1)

    # 2. Create Window Steps (0, 1, ..., window_size-1)
    #    Shape: (1, window_size)
    window_steps = torch.arange(window_size, device=device).unsqueeze(0)

    # 3. Broadcast to get full index map
    #    Shape: (n_windows, window_size)
    target_indices = base_indices + window_steps

    # 4. Flatten indices to match the flattened gradients
    #    Shape: (Total_Elements,)
    target_indices_flat = target_indices.view(-1)

    # ============================================================
    # Step C: Accumulate Gradients (Scatter Add)
    # ============================================================
    # 3. Accumulation (index_add_): This is the crucial step. Since multiple windows
    #    might overlap (read from the same input index), their gradients must be summed.
    #    index_add_ handles this atomically.
    #
    #     Note on Padding: We allocate a buffer that includes the padding size.
    #     Gradients computed for "padded zeros" are accumulated into the padding
    #     regions of this buffer.

    # Calculate padded length
    max_k = window_start_offsets.max().item()
    needed_length = max_k + window_size
    pad_bottom = max(0, needed_length - in_len)

    # 1. Create a zero-filled buffer for the PADDED input gradient
    #    Length = pad_top + in_len + pad_bottom
    total_padded_len = pad_top + in_len + pad_bottom
    num_features = grad_source_flat.shape[1]

    grad_padded_buffer = torch.zeros((total_padded_len, num_features), device=device, dtype=grad_output.dtype)

    # 2. Perform the Accumulation
    #    This is the "Overlap-Add" magic. Gradients from overlapping windows
    #    are summed up automatically.
    #    Here grad_source_flat[i] is the upstream adjoint of the fwd output
    #    while target_indices_flat[i] is the index of the fwd padded_vector
    #    and the gradient thereof
    grad_padded_buffer.index_add_(0, target_indices_flat, grad_source_flat)

    # ============================================================
    # Step D: Handle Padding & Reshape
    # ============================================================
    # 4. Slicing: In the final step, we simply slice out the valid middle region
    #    (pad_top to pad_top + in_len), effectively discarding the gradients that
    #    accumulated in the "virtual" padded zones.

    # 1. Slice off the padding (Discard gradients that fell into the pad zones)
    #    We only keep indices [pad_top : pad_top + in_len]
    #    NOTE that in the actual usage case corresponding to the get_indexing_matrix
    #    and single_to_keys, there should be an input mask that would go thru
    #    the same single_to_keys operation as do the sequence data so the mask
    #    would also result in zeros corresponding to the pad zones in the padded_vector,
    #    which implies that it's safe to discard the gradients that fell into the pad zones here
    grad_input_flat = grad_padded_buffer[pad_top : pad_top + in_len]

    # 2. Reshape back to the original input geometry
    #    We flattened (..., In_Len, ...) into (In_Len, Features).
    #    We need to reverse this.

    #    A. Calculate dimensions before and after 'axis'
    #       input_shape = (D1, D2, In_Len, D3, D4)
    #       We need to reshuffle grad_input_flat (In_Len, D1*D2*D3*D4)
    #       back to that shape.

    #    It is cleaner to use the original shape directly but permuted.
    #    Target Layout for reshape: (In_Len, ...)

    #    Construct the permuted shape where 'axis' is at dim 0
    permuted_shape = list(input_shape)
    permuted_shape.pop(axis)
    permuted_shape.insert(0, in_len)

    #    Reshape
    grad_input = grad_input_flat.reshape(permuted_shape)

    #    Inverse Permutation: Move dim 0 back to 'axis'
    grad_input = grad_input.moveaxis(0, axis)

    return grad_input


class GatherSlidingWindows(torch.autograd.Function):
    @staticmethod
    def forward(ctx, input, window_start_offsets, window_size, axis):
        """
        Gather overlapping sliding windows from input at specified starting positions.

        This operation implements efficient windowed attention by extracting windows
        from the input sequence. The underlying mathematical structure is a block
        Toeplitz matrix (see Theorems 1-6 in documentation).

        Example & Intuition:
        --------------------
        Consider gathering 3 windows of size 8 from a sequence of length 10.
        The operation can be viewed as a sparse matrix 'M' (n_windows=3, window_size=8, seq_len=10)
        where each window gathers contiguous elements from the input.

        Semantic mapping to get_indexing_matrix:
        - axis 0: query window index
        - axis 1: slot within window, i.e., index in [0, h-1] where h == H // (W // 2)
        - axis 2: input sequence position, index to the "2 * K" half-windows dimension

        1. Visualizing as Sparse Matrix 'M':
           NOTE: This corresponds to the transposed onehot tensor from get_indexing_matrix:
           M == onehot.transpose(1, 0).transpose(-2, -1)

           The '1's represent which input index is gathered to the output.
           Notice the diagonal pattern shifts by 2 between windows.

           Window 0 (starting at position -3): "Pad Left"
           [0 0 0 0 0 0 0 0 0 0] <- Slot 0 (padded zero)
           [0 0 0 0 0 0 0 0 0 0] <- Slot 1 (padded zero)
           [0 0 0 0 0 0 0 0 0 0] <- Slot 2 (padded zero)
           [1 0 0 0 0 0 0 0 0 0] <- Slot 3 (gathers Input[0])
           [0 1 0 0 0 0 0 0 0 0] <- Slot 4 (gathers Input[1])
           [0 0 1 0 0 0 0 0 0 0] ...
           [0 0 0 1 0 0 0 0 0 0]
           [0 0 0 0 1 0 0 0 0 0]

           Window 1 (starting at position -1): "Shifted +2 from Window 0"
           [0 0 0 0 0 0 0 0 0 0] <- Slot 0 (padded zero)
           [1 0 0 0 0 0 0 0 0 0] <- Slot 1 (gathers Input[0])
           [0 1 0 0 0 0 0 0 0 0] <- Slot 2 (gathers Input[1])
           [0 0 1 0 0 0 0 0 0 0] ...
           [0 0 0 1 0 0 0 0 0 0]
           [0 0 0 0 1 0 0 0 0 0]
           [0 0 0 0 0 1 0 0 0 0]
           [0 0 0 0 0 0 1 0 0 0]

           Window 2 (starting at position +1): "Shifted +2 from Window 1"
           [0 1 0 0 0 0 0 0 0 0] <- Slot 0 (gathers Input[1])
           [0 0 1 0 0 0 0 0 0 0] <- Slot 1 (gathers Input[2])
           [0 0 0 1 0 0 0 0 0 0] ...
           [0 0 0 0 1 0 0 0 0 0]
           [0 0 0 0 0 1 0 0 0 0]
           [0 0 0 0 0 0 1 0 0 0]
           [0 0 0 0 0 0 0 1 0 0]
           [0 0 0 0 0 0 0 0 1 0]

        2. The "Unfold" Implementation:
           Instead of explicit matrix multiplication, we use torch.unfold to create
           sliding windows efficiently. If input is 2D (seq_len, features), gathering
           windows via this operation is equivalent to applying the sparse matrix M.

           For computational efficiency, we slide a window over the (padded) input:
           - Window 0 (offset=-3): reads positions [-3] to [4] (with padding)
           - Window 2 (offset=+1): reads positions [1] to [8]

        Args:
            input: Tensor of arbitrary shape (..., seq_len, ...)
            window_start_offsets: (n_windows,) tensor of integer starting positions
            window_size: int, the size of each output window (h)
            axis: int, the dimension corresponding to sequence length

        Returns:
            Tensor of shape (..., n_windows, window_size, ...)
            1. The 'n_windows' dimension replaces the original sequence dimension at 'axis'
            2. The 'window_size' dimension is placed immediately after 'n_windows' (at axis + 1)
        """

        # 1. Normalize and Save Axis/Shape info
        ndim = input.ndim
        if axis < 0:
            axis += ndim

        in_len = input.shape[axis]

        # 2. Analyze Padding
        min_k = window_start_offsets.min().item()
        max_k = window_start_offsets.max().item()

        pad_top = max(0, -min_k)
        needed_length = max_k + window_size
        pad_bottom = max(0, needed_length - in_len)

        ctx.mark_non_differentiable(window_start_offsets)

        # Save context for backward
        # We save shapes and integers, and the window_start_offsets tensor.
        # We DO NOT save the input (saves memory).
        ctx.save_for_backward(window_start_offsets)
        ctx.params = {
            "input_shape": input.shape,
            "pad_top": pad_top,
            "pad_bottom": pad_bottom,
            "output_len": window_size,
            "axis": axis,
            "in_len": in_len,
        }

        # 3. Forward Logic: Create sliding windows via unfold
        if pad_top == 0 and pad_bottom == 0:
            padded_vector = input
        else:
            pad_arg = [0] * (2 * ndim)
            pad_idx_left = (ndim - 1 - axis) * 2
            pad_idx_right = pad_idx_left + 1
            pad_arg[pad_idx_left] = pad_top
            pad_arg[pad_idx_right] = pad_bottom
            padded_vector = torch.nn.functional.pad(input, pad_arg)

        # Shape: (..., num_windows, ..., window_size)
        # - Position 'axis' now contains num_windows (replaces padded_len)
        # - New dimension window_size added at position -1 (end)
        windows = padded_vector.unfold(axis, window_size, 1)

        # Shape: (n_windows,) - translate window_start_offsets to padded coordinate system
        slice_indices = window_start_offsets + pad_top

        # Shape: (..., n_windows, ..., window_size) - select specific windows along axis
        selected_windows = windows.index_select(axis, slice_indices)

        # Permute: (..., n_windows, ..., window_size) -> (..., n_windows, window_size, ...)
        result = selected_windows.moveaxis(-1, axis + 1)

        return result

    @staticmethod
    def backward(ctx, grad_output):
        """
        Explanation of the Backward Logic

        1. Preparation (Flattening): Because the input can have arbitrary dimensions
        (e.g., Batch, Time, Features or Channels, Time, Height, Width), performing
        operations on specific dimensions is tricky. In Step A, we permute and
        flatten the tensor into a 2D matrix: [N_Windows * Output_Len, Flat_Features].
        This standardizes the problem regardless of the input shape.

        2. Mapping (target_indices): This is the inverse of the forward index_select.

            Forward: "For window w, read index i."

            Backward: "For window w, gradient i belongs to index i+offset[w]."
                      We generate a full grid of these destination indices.

        3. Accumulation (index_add_): This is the crucial step. Since multiple windows
           might overlap (read from the same input index), their gradients must be summed.
           index_add_ handles this atomically.

            Note on Padding: We allocate a buffer that includes the padding size.
            Gradients computed for "padded zeros" are accumulated into the padding
            regions of this buffer.

        4. Slicing: In the final step, we simply slice out the valid middle region
           (pad_top to pad_top + in_len), effectively discarding the gradients that
           accumulated in the "virtual" padded zones.
        """
        # Retrieve context
        (window_start_offsets,) = ctx.saved_tensors
        params = ctx.params
        input_shape = params["input_shape"]
        window_size = params["output_len"]  # Stored as output_len in context for backward compat
        axis = params["axis"]

        # Call standalone backward function (includes validation)
        grad_input = gather_sliding_windows_backward(grad_output, window_start_offsets, window_size, axis, input_shape)

        return grad_input, None, None, None


def compute_query_window_ownership(W: int, H: int, K: int, qw_start: int, qw_end: int) -> dict:
    """
    Compute halo requirements for a given query window ownership.

    Parameters
    ----------
    W : int
        Atoms per query window (must be even)
    H : int
        Keys per query window (must be divisible by W//2)
    K : int
        Total number of query windows (global)
    qw_start : int
        First owned query window (inclusive)
    qw_end : int
        Last owned query window + 1 (exclusive)

    Returns
    -------
    dict
        {
            'hw_owned': (int, int),              # Owned half-windows [2*qw_start, 2*qw_end)
            'hw_needed': (int, int),             # All half-windows needed [start, end)
            'left_halo_size': int,               # Half-windows needed from left neighbor
            'right_halo_size': int,              # Half-windows needed from right neighbor
        }

    Examples
    --------
    >>> # Rank owns QW[4,8) for K=12
    >>> ownership = compute_query_window_ownership(32, 128, 12, 4, 8)
    >>> ownership['hw_owned']
    (8, 16)  # Owns HW[8-15] (inferred from QW range)
    >>> ownership['hw_needed']
    (5, 19)  # Needs HW[5-18]
    >>> ownership['left_halo_size']
    3  # Needs HW[5,6,7] from left neighbor
    """
    assert W > 0 and W % 2 == 0
    assert H > 0 and H % (W // 2) == 0
    assert K > 0
    assert 0 <= qw_start <= qw_end <= K

    # Infer half-window ownership from query window ownership
    # Query window i owns half-windows [2i, 2i+1]
    hw_start = 2 * qw_start
    hw_end = 2 * qw_end

    # Determine which half-windows are needed for owned query windows
    if qw_start < qw_end:
        owned_qw_ids = torch.arange(qw_start, qw_end)
        ranges = get_query_window_key_range(W, H, K, owned_qw_ids)

        hw_need_start = ranges[0].min().item()
        hw_need_end = ranges[1].max().item() + 1  # Exclusive end

        # Compute halo sizes
        left_halo_size = max(0, hw_start - hw_need_start)
        right_halo_size = max(0, hw_need_end - hw_end)
    else:
        # No owned query windows
        hw_need_start = hw_start
        hw_need_end = hw_start
        left_halo_size = 0
        right_halo_size = 0

    return {
        "hw_owned": (hw_start, hw_end),
        "hw_needed": (hw_need_start, hw_need_end),
        "left_halo_size": left_halo_size,
        "right_halo_size": right_halo_size,
    }


def get_halo_from_neighbors(
    rank: int,
    size_group: int,
    n_half_windows_local: int,
    W: int,
    H: int,
    K: int,
) -> tuple[list, list]:
    """
    Compute send/recv metadata for halo exchange (supports multi-hop).

    Returns:
        tuple(recv_meta, send_meta)
        - recv_meta: list of (peer_rank, halo_type, offset_in_halo, length)
          where halo_type is 'left' or 'right'.
        - send_meta: list of (peer_rank, offset_in_local, length)
    """
    # 1. Validate inputs
    assert W > 0 and W % 2 == 0, "W must be positive and even"
    assert H > 0 and H % (W // 2) == 0, "H must be divisible by W//2"
    assert H // (W // 2) % 2 == 0, "H // (W // 2) must be even"
    assert K > 0, "K must be positive"
    if K % size_group != 0:
        raise ValueError(f"K {K} must be an integer multiple of the number of ranks {size_group}.")

    # 2. Vectorized ownership computation
    # Rank ownership: [hw_start, hw_end)
    rank_ids = torch.arange(size_group)
    hw_owned_starts = rank_ids * n_half_windows_local
    hw_owned_ends = (rank_ids + 1) * n_half_windows_local

    # Query window ownership: [qw_start, qw_end)
    qw_starts = hw_owned_starts // 2
    qw_ends = hw_owned_ends // 2

    # Needed half-windows: [hw_need_start, hw_need_end)
    # Note: K is large, so we don't want to compute range for all query windows individually.
    # But compute_query_window_ownership works per rank range.
    # We can vectorize get_query_window_key_range over the start/end QWs of all ranks.

    # For a range of QWs [qs, qe), the needed range is union of needs of all q in [qs, qe).
    # Since QWs are monotonic, needed range is [min(need(qs)), max(need(qe-1))].
    # Exception: if qs >= qe (empty rank), need range is empty/irrelevant.

    # Compute needs for first QW of each rank
    ranges_start = get_query_window_key_range(W, H, K, qw_starts)
    hw_need_starts = ranges_start[0]  # j_min of first QW

    # Compute needs for last QW of each rank
    # Use (qw_ends - 1) but clamp to >= 0 to avoid index -1 for empty ranks
    last_qw_ids = (qw_ends - 1).clamp(min=0)
    ranges_end = get_query_window_key_range(W, H, K, last_qw_ids)
    hw_need_ends = ranges_end[1] + 1  # j_max + 1 of last QW (exclusive end)

    # Handle empty ranks: if qw_start >= qw_end, they own nothing and need nothing (locally)
    # We set need = owned to result in 0 halo size
    is_empty = qw_starts >= qw_ends
    hw_need_starts = torch.where(is_empty, hw_owned_starts, hw_need_starts)
    hw_need_ends = torch.where(is_empty, hw_owned_starts, hw_need_ends)

    # Get current rank's values
    my_hw_start = hw_owned_starts[rank].item()
    my_hw_end = hw_owned_ends[rank].item()
    my_need_start = hw_need_starts[rank].item()
    my_need_end = hw_need_ends[rank].item()

    recv_meta = []  # (peer, halo_type, offset_in_halo, length)
    send_meta = []  # (peer, offset_in_local, length)

    # 3. Identify Neighbors via SearchSorted

    # --- RECV Left ---
    # Need neighbors covering [my_need_start, my_hw_start)
    if my_need_start < my_hw_start:
        # Find first rank whose owned range ends after my_need_start
        # i.e. hw_owned_ends[p] > my_need_start
        p_start_idx = torch.searchsorted(hw_owned_ends, my_need_start, side="right")
        # Find last rank whose owned range starts before my_hw_start
        # i.e. hw_owned_starts[p] < my_hw_start
        p_end_idx = torch.searchsorted(hw_owned_starts, my_hw_start, side="left")

        for peer in range(p_start_idx, p_end_idx):
            if peer == rank:
                continue
            # Overlap: [max(my_need, peer_start), min(my_start, peer_end))
            l_start = max(my_need_start, hw_owned_starts[peer].item())
            l_end = min(my_hw_start, hw_owned_ends[peer].item())
            if l_start < l_end:
                recv_meta.append((peer, "left", l_start - my_need_start, l_end - l_start))

    # --- RECV Right ---
    # Need neighbors covering [my_hw_end, my_need_end)
    if my_hw_end < my_need_end:
        # First rank whose owned range ends after my_hw_end
        p_start_idx = torch.searchsorted(hw_owned_ends, my_hw_end, side="right")
        # Last rank whose owned range starts before my_need_end
        p_end_idx = torch.searchsorted(hw_owned_starts, my_need_end, side="left")

        for peer in range(p_start_idx, p_end_idx):
            if peer == rank:
                continue
            # Overlap: [max(my_end, peer_start), min(my_need_end, peer_end))
            r_start = max(my_hw_end, hw_owned_starts[peer].item())
            r_end = min(my_need_end, hw_owned_ends[peer].item())
            if r_start < r_end:
                recv_meta.append((peer, "right", r_start - my_hw_end, r_end - r_start))

    # --- SEND Left (Peers needing me for their left halo) ---
    # Condition: peer > rank AND peer_need_start < my_hw_end
    # Range of peers: (rank + 1, ...) such that peer_need_start < my_hw_end
    # Since hw_need_starts is sorted (monotonic with rank), we can search
    if rank + 1 < size_group:
        # Find last peer where peer_need_start < my_hw_end
        # searchsorted on hw_need_starts to find insertion point of my_hw_end
        limit_idx = torch.searchsorted(hw_need_starts, my_hw_end, side="left")
        # Valid peers are in range [rank + 1, limit_idx)
        # Note: limit_idx is where value >= my_hw_end starts, so up to limit_idx-1 are < my_hw_end

        for peer in range(rank + 1, limit_idx):
            # Overlap: Peer needs [p_need_start, p_start), I own [my_start, my_end)
            # Intersection: [max(p_need_start, my_start), min(p_start, my_end))
            p_need_start = hw_need_starts[peer].item()
            p_hw_start = hw_owned_starts[peer].item()

            l_start = max(p_need_start, my_hw_start)
            l_end = min(p_hw_start, my_hw_end)
            if l_start < l_end:
                send_meta.append((peer, l_start - my_hw_start, l_end - l_start))

    # --- SEND Right (Peers needing me for their right halo) ---
    # Condition: peer < rank AND peer_need_end > my_hw_start
    # Range of peers: (..., rank - 1) such that peer_need_end > my_hw_start
    # Since hw_need_ends is sorted, we can search
    if rank > 0:
        # Find first peer where peer_need_end > my_hw_start
        # searchsorted on hw_need_ends with my_hw_start
        start_idx = torch.searchsorted(hw_need_ends, my_hw_start, side="right")
        # Valid peers are in range [start_idx, rank)

        for peer in range(start_idx, rank):
            # Overlap: Peer needs [p_end, p_need_end), I own [my_start, my_end)
            # Intersection: [max(p_end, my_start), min(p_need_end, my_end))
            p_hw_end = hw_owned_ends[peer].item()
            p_need_end = hw_need_ends[peer].item()

            r_start = max(p_hw_end, my_hw_start)
            r_end = min(p_need_end, my_hw_end)
            if r_start < r_end:
                send_meta.append((peer, r_start - my_hw_start, r_end - r_start))

    return recv_meta, send_meta


def pack_and_pad(
    input: torch.Tensor, mask: torch.Tensor, axis: int, W: int, keep_input_padding: bool = False
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """
    Left-pack valid elements from input and pad to the next multiple of W.

    This utility prepares variable-length inputs (with padding masks) for downstream
    operations like gather_sliding_windows by:
    1. Left-packing valid elements (moving all True-masked values to the front)
    2. Padding the sequence length to a multiple of W

    Note: This function does NOT reshape the output. Any reshaping (e.g., to
    (2*K, W//2) for half-window layout) should be done by the caller after
    this function returns.

    Process:
    1. Validate inputs (shapes, broadcastability, W is positive and even)
    2. Determine target_len based on keep_input_padding flag
    3. Pad input/mask to target_len if needed
    4. Sort mask descending (stable) to left-pack valid elements
    5. Gather input elements using the sort indices
    6. Zero out invalid positions in the packed output

    Parameters
    ----------
    input : torch.Tensor
        Input tensor of arbitrary shape (..., seq_len, ...)
    mask : torch.Tensor
        Boolean mask with mask.shape[axis] == input.shape[axis].
        Other dimensions must be broadcastable with input.
        True=valid element, False=padding to ignore.
    axis : int
        Dimension containing the sequence
    W : int
        Padding factor (must be positive and even). The output length along
        axis will be padded to the next multiple of W.
    keep_input_padding : bool, optional
        If False (default), the output length is based on the maximum number
        of valid elements across all slices orthogonal to axis, padded to
        the next multiple of W. This removes as much padding as possible.
        If True, the output length is based on input.shape[axis], padded
        to the next multiple of W. This preserves the original sequence length.

    Returns
    -------
    tuple[torch.Tensor, torch.Tensor, torch.Tensor]
        (packed_output, gather_indices, packed_mask)

        packed_output : torch.Tensor
            Left-packed and padded tensor of shape (..., target_len, ...).
            Valid elements are at the front, followed by zeros.
            target_len is a multiple of W.

        gather_indices : torch.Tensor
            Indices used for gathering, shape (..., target_len, ...).
            Required for pack_and_pad_backward to scatter gradients back.

        packed_mask : torch.Tensor
            Boolean mask for packed_output, shape (..., target_len, ...).
            True for valid elements, False for padding.
    """
    # 1. Validate inputs
    if not isinstance(mask, torch.Tensor):
        raise TypeError(f"mask must be a torch.Tensor, got {type(mask)}")
    if mask.dtype != torch.bool:
        raise TypeError(f"mask must have dtype torch.bool, got {mask.dtype}")
    if W <= 0 or W % 2 != 0:
        raise ValueError(f"W must be positive and even, got {W}")
    if not isinstance(keep_input_padding, bool):
        raise TypeError(f"Expected bool for keep_input_padding, got {type(keep_input_padding)}")

    # Normalize axis
    ndim = input.ndim
    if axis < 0:
        axis += ndim
    if not (0 <= axis < ndim):
        raise ValueError(f"axis {axis} out of range for {ndim}D input")

    # Validate mask shape compatibility along axis
    if input.ndim != mask.ndim:
        raise ValueError(f"mask ndim {mask.ndim} must match input ndim {input.ndim}")

    if input.shape[axis] != mask.shape[axis]:
        raise ValueError(
            f"mask and input must match along axis={axis}: "
            f"input.shape[{axis}]={input.shape[axis]}, mask.shape[{axis}]={mask.shape[axis]}"
        )

    # Check broadcastability (other dims)
    try:
        broadcasted_shape = torch.broadcast_shapes(input.shape, mask.shape)
    except RuntimeError as e:
        raise ValueError(f"mask shape {mask.shape} not broadcastable to input shape {input.shape}: {e}")

    if broadcasted_shape != input.shape:
        raise ValueError(
            f"mask shape {mask.shape} broadcasts to {broadcasted_shape}, which mismatches input {input.shape}"
        )

    # 2. Determine target length based on keep_input_padding flag
    if keep_input_padding:
        # Use original input length as basis for padding calculation
        len_basis = input.shape[axis]
    else:
        # Use max valid count to minimize output length (removes excess padding)
        valid_counts = mask.sum(dim=axis, dtype=torch.long)
        max_valid = valid_counts.max().item()
        len_basis = max_valid

    # Round up to next multiple of W
    target_len = ((len_basis + W - 1) // W) * W

    # 3. Pad input and mask to target_len if needed
    current_len = input.shape[axis]
    if target_len > current_len:
        pad_len = target_len - current_len

        # Build padding arguments for torch.nn.functional.pad
        # Format: (left_dim_N, right_dim_N, ..., left_dim_0, right_dim_0)
        pad_arg = [0] * (2 * input.ndim)
        pad_idx = (input.ndim - 1 - axis) * 2 + 1  # Index for right-padding along axis
        pad_arg[pad_idx] = pad_len

        input_padded = torch.nn.functional.pad(input, pad_arg)
        mask_padded = torch.nn.functional.pad(mask, pad_arg)
    else:
        input_padded = input
        mask_padded = mask

    # 4. Left-pack valid elements using stable descending sort on mask
    # Sorting True (1) before False (0) in descending order moves valid elements to the front
    mask_padded_sorted, argsort_mask_padded = torch.sort(mask_padded, dim=axis, descending=True, stable=True)

    # Slice to target_len (handles case where input was longer than target_len)
    slices = [slice(None)] * mask_padded.ndim
    slices[axis] = slice(0, target_len)
    argsort_mask_padded = argsort_mask_padded[tuple(slices)]
    mask_padded_sorted = mask_padded_sorted[tuple(slices)]

    # 5. Expand indices to match input dimensions for gathering
    # torch.gather requires index tensor to have same ndim as input.
    # Expanding broadcasts the sort indices across non-axis dimensions (e.g., features).
    target_gather_shape = list(input_padded.shape)
    target_gather_shape[axis] = target_len
    argsort_mask_padded_expanded = argsort_mask_padded.expand(target_gather_shape)
    mask_padded_sorted_expanded = mask_padded_sorted.expand(target_gather_shape)

    # Gather input elements according to the left-packing order
    input_packed_padded = torch.gather(input_padded, axis, argsort_mask_padded_expanded)

    # 6. Zero out invalid positions
    # The gather operation may have pulled "garbage" values from original invalid positions.
    # Multiplying by the sorted mask ensures only valid elements are non-zero.
    input_packed_padded = input_packed_padded * mask_padded_sorted_expanded.to(input.dtype)

    return input_packed_padded, argsort_mask_padded_expanded, mask_padded_sorted_expanded


def pack_and_pad_backward(
    grad_output: torch.Tensor,
    mask_output: torch.Tensor,
    indices: torch.Tensor,
    input_shape: tuple,
    axis: int,
) -> torch.Tensor:
    """
    Backward pass for pack_and_pad.

    Scatters gradients from the packed/padded output back to the original input shape.
    This reverses the gather operation by using scatter_add with the same indices.

    Parameters
    ----------
    grad_output : torch.Tensor
        Gradient w.r.t. packed output, shape (..., target_len, ...).
        This is the gradient flowing back from operations applied to pack_and_pad's output.
    mask_output : torch.Tensor
        Mask from pack_and_pad's forward pass, shape (..., target_len, ...).
        Used to zero out gradients for invalid (padding) positions.
    indices : torch.Tensor
        Gather indices from pack_and_pad's forward pass (argsort_mask_padded_expanded),
        shape (..., target_len, ...). Used to scatter gradients back to original positions.
    input_shape : tuple
        Shape of the original input tensor to pack_and_pad.
    axis : int
        Dimension containing the sequence (same as in forward pass).

    Returns
    -------
    torch.Tensor
        Gradient w.r.t. original input, shape matching input_shape.
    """
    # Mask out gradients for invalid positions (padding)
    grad_masked = grad_output * mask_output

    # Determine buffer size - may need extra space if target_len > original input length
    target_len = indices.shape[axis]
    current_len = input_shape[axis]
    max_len = max(target_len, current_len)

    # Create gradient buffer for scatter operation
    if max_len > current_len:
        buffer_shape = list(input_shape)
        buffer_shape[axis] = max_len
        grad_input_padded = torch.zeros(buffer_shape, dtype=grad_masked.dtype, device=grad_masked.device)
    else:
        grad_input_padded = torch.zeros(input_shape, dtype=grad_masked.dtype, device=grad_masked.device)

    # Scatter gradients back to their original positions
    # This reverses the gather operation from the forward pass
    grad_input_padded.scatter_add_(axis, indices, grad_masked)

    # Slice back to original input shape if we used an extended buffer
    if max_len > current_len:
        slices = [slice(None)] * len(input_shape)
        slices[axis] = slice(0, current_len)
        grad_input = grad_input_padded[tuple(slices)]
    else:
        grad_input = grad_input_padded

    return grad_input


def gather_sliding_windows(input, window_start_offsets, window_size, axis):
    """
    Gather sliding windows from input using specified offsets.

    This operation implements windowed attention by extracting overlapping windows
    from the input sequence. The underlying mathematical structure is a block Toeplitz
    matrix (see Theorems 1-6 in documentation).

    Args:
        input: Input tensor of shape (..., sequence_len, ...)
        window_start_offsets: Starting positions for each window, shape (n_windows,)
        window_size: Size of each window (h)
        axis: Dimension along which to gather windows

    Returns:
        Tensor of shape (..., n_windows, window_size, ...)
    """
    return GatherSlidingWindows.apply(input, window_start_offsets, window_size, axis)


def get_flattened_range_indices(start_ends: torch.Tensor, device: torch.device) -> tuple[torch.Tensor, torch.Tensor]:
    """
    Generate flattened row and column indices for a set of ranges per row.

    Args:
        start_ends: (D, 2) tensor where col 0 is start, col 1 is end (exclusive).
        device: torch device.

    Returns:
        (row_indices, col_indices) tuple of 1D tensors.
    """
    starts = start_ends[:, 0].to(torch.long)
    ends = start_ends[:, 1].to(torch.long)
    lengths = (ends - starts).clamp(min=0)

    if lengths.sum() == 0:
        return torch.empty(0, dtype=torch.long, device=device), torch.empty(0, dtype=torch.long, device=device)

    # 1. Generate Row Indices
    row_indices = torch.repeat_interleave(torch.arange(start_ends.shape[0], device=device), lengths)

    # 2. Generate Column Indices
    # Global flat counter
    total_length = lengths.sum()
    flat_range = torch.arange(total_length, device=device)

    # Offsets for each row in the flat sequence
    cum_lengths = torch.cumsum(lengths, dim=0)
    shifts = torch.zeros_like(cum_lengths)
    shifts[1:] = cum_lengths[:-1]

    # Expand shifts to match flat range
    shifts_expanded = torch.repeat_interleave(shifts, lengths)
    starts_expanded = torch.repeat_interleave(starts, lengths)

    # col_idx = start + relative_idx
    col_indices = starts_expanded + (flat_range - shifts_expanded)

    return row_indices, col_indices


def _distributed_pack_and_pad(
    input: Optional[DTensor], mask: DTensor, axis: int, W: int, keep_input_padding: bool = False
) -> tuple[Optional[DTensor], DTensor, torch.Tensor, torch.Tensor]:
    """Distributed left-pack and pad operation for DTensor inputs.

    Left-packs valid elements (as indicated by mask) to the front along the specified axis,
    then pads to the right (end) to the next multiple of (W * size_group) where size_group
    is the number of ranks sharding the tensor along axis. Communication is performed to redistribute
    elements across ranks so that the packed result is evenly sharded.

    Note: This function does NOT reshape the output. Any reshaping (e.g., to
    (2*K, W//2) for half-window layout) should be done by the caller.

    Args:
        input (Optional[DTensor]): Input DTensor of shape (..., seq_len, ...) sharded along axis,
            or None. If None, input is set to mask (useful for computing metadata only).
        mask (DTensor): Boolean mask DTensor with shape broadcastable to input. True indicates
            valid elements, False indicates padding to ignore.
        axis (int): Dimension containing the sequence to pack and pad.
        W (int): Padding factor. The output length along axis will be padded to the next
            multiple of (W * size_group).
        keep_input_padding (bool, optional): If False (default), output length is based on
            the maximum number of valid elements across all slices. If True, output length
            is based on input.shape[axis]. Defaults to False.

    Returns:
        tuple[Optional[DTensor], DTensor, torch.Tensor, torch.Tensor]:
            output: Packed and padded DTensor, or None if input was None. Shape is
                (..., target_len, ...) where target_len is a multiple of (W * size_group).
            mask_output: Boolean mask DTensor for output, same shape as output.
            argsort_mask_flat_local: Local argsort indices (2D tensor of shape
                (shape_leading_flat_local, shape_axis_local)) used for backward pass.
            valid_counts_all_ranks: Valid counts per rank, shape (size_group, shape_leading_flat_local).
    """
    # 0. sanity checks
    if input is not None and not isinstance(input, DTensor):
        raise TypeError(f"Expected DTensor, got {type(input)}")
    if not isinstance(mask, DTensor):
        raise TypeError(f"Expected DTensor, got {type(mask)}")
    if not isinstance(axis, int):
        raise TypeError(f"Expected int for axis, got {type(axis)}")
    if not isinstance(W, int):
        raise TypeError(f"Expected int for W, got {type(W)}")
    if not isinstance(keep_input_padding, bool):
        raise TypeError(f"Expected bool for keep_input_padding, got {type(keep_input_padding)}")

    # Mask must be boolean to ensure correct valid count computation
    if mask.dtype != torch.bool:
        raise TypeError(
            f"mask must have dtype torch.bool to avoid precision issues in valid count computation, "
            f"got {mask.dtype}. Use mask.bool() to convert."
        )

    if input is None:
        input = mask
        has_input = False
    else:
        has_input = True
        if input.device_mesh != mask.device_mesh:
            raise ValueError(
                f"input and mask must have the same device mesh but got {input.device_mesh} and {mask.device_mesh}"
            )
        if input.placements != mask.placements:
            raise ValueError(
                f"input and mask must have the same placements but got {input.placements} and {mask.placements}"
            )

    placements = input.placements
    device_mesh = input.device_mesh

    i_dim_device_mesh_shard_axis = None
    for i_dim_device_mesh, placement in enumerate(placements):
        if isinstance(placement, Partial):
            raise ValueError("Partial placements are not supported")
        elif isinstance(placement, Shard):
            if input.shape[placement.dim] % device_mesh.shape[i_dim_device_mesh] != 0:
                raise ValueError(
                    f"Uneven sharding tensor dimension {placement.dim} of size {input.shape[placement.dim]} "
                    f"along device mesh dimension {i_dim_device_mesh} of size {device_mesh.shape[i_dim_device_mesh]} is not supported"
                )
            if placement.dim == axis:
                i_dim_device_mesh_shard_axis = i_dim_device_mesh

    if i_dim_device_mesh_shard_axis is None:
        raise ValueError(f"input is not sharded along axis {axis}")

    try:
        shape_broadcast = torch.broadcast_shapes(input.shape, mask.shape)
    except RuntimeError as e:
        raise ValueError("input and mask shapes cannot be broadcasted") from e

    if shape_broadcast != input.shape:
        raise ValueError(f"broadcasted shape {shape_broadcast} is not equal to input.shape {input.shape}")

    # 1. Get rank info from passed mesh
    rank = device_mesh.get_local_rank(i_dim_device_mesh_shard_axis)
    size_group = device_mesh.size(i_dim_device_mesh_shard_axis)
    group = device_mesh.get_group(i_dim_device_mesh_shard_axis)
    global_rank_peers = torch.distributed.get_process_group_ranks(group)

    local_input = input.to_local()
    local_mask = mask.to_local()
    device = local_input.device

    # Explicitly broadcast mask to match input shape locally
    if local_mask.shape != local_input.shape:
        local_mask = local_mask.expand(local_input.shape)

    # --- Phase 1: Forward Logic ---

    # 1.1 Reshape to 2D (shape_leading_flat_local, N)
    # Move axis to last dim
    local_input_moved = local_input.movedim(axis, -1)
    local_mask_moved = local_mask.movedim(axis, -1)

    shape_leading_flat_local = local_input_moved.numel() // local_input_moved.shape[-1]
    shape_axis_local = local_input_moved.shape[-1]

    local_input_2d = local_input_moved.reshape(shape_leading_flat_local, shape_axis_local)
    local_mask_2d = local_mask_moved.reshape(shape_leading_flat_local, shape_axis_local)

    # 1.2 Local Sort/Pack
    local_valid_counts = local_mask_2d.sum(dim=1, dtype=torch.long)  # (shape_leading_flat_local,)

    # Sort mask to get valid elements to the left
    # local_mask_sorted.shape == argsort_mask_flat_local.shape == local_mask_2d.shape == (shape_leading_flat_local, shape_axis_local)
    local_mask_sorted, argsort_mask_flat_local = torch.sort(local_mask_2d, dim=1, descending=True, stable=True)

    # Gather valid data
    local_valid_data = torch.gather(local_input_2d, 1, argsort_mask_flat_local)

    # Zero out invalid elements
    local_valid_data = local_valid_data * local_mask_sorted.to(local_input_2d.dtype)

    # 1.3 Global Planning (Per Row)
    valid_counts_all_ranks = torch.zeros(size_group, shape_leading_flat_local, device=device, dtype=torch.long)
    dist.all_gather_into_tensor(
        valid_counts_all_ranks, local_valid_counts.unsqueeze(0).to(dtype=torch.long), group=group
    )

    global_ends = valid_counts_all_ranks.cumsum(dim=0)
    global_starts = global_ends - valid_counts_all_ranks
    total_valid = global_ends[-1]  # (shape_leading_flat_local,)

    my_global_start = global_starts[rank]  # (shape_leading_flat_local,)
    my_global_end = my_global_start + local_valid_counts  # (shape_leading_flat_local,)

    # 1.4 Target Partitioning (Global Uniform)

    if keep_input_padding:
        # NOTE: we use input.shape[axis] instead of max_total_valid =mask.sum(dim=axis).max()
        # to stay consistent with the Boltz implementation.
        # This is slightly inefficient because the trailing invalid elements can be reduced
        # by directly pad towards (W * size_group) based on max_total_valid
        # instead of the original input length
        len_basis = input.shape[axis]
    else:
        # Use single global target length derived from MAX row length to ensure continuity.
        # This prevents gaps in valid sequence data across rank boundaries for short rows.
        max_total_valid = total_valid.max().item()
        # while the other comm ops are within the "group", the (shape_leading_flat_local,) axis is virtually
        # sharded if other axes than "axis" in the input are also sharded so we need to
        # do a global all_reduce to get the max_total_valid across the sharded (shape_leading_flat_local,) axis
        tensor_max_total_valid = torch.tensor(max_total_valid, device=device, dtype=torch.long)
        for i_subgroup, subgroup in enumerate(device_mesh.get_all_groups()):
            # we can't use the default world group because the input device_mesh can be
            # a subgroup of the world group, e.g., device_mesh is a submesh. Also, DeviceMesh
            # has no API to return the union of its subgroups so we need to iterate over all subgroups
            if i_subgroup == i_dim_device_mesh_shard_axis:
                # already reduce via all_gather and local max()
                continue
            torch.distributed.all_reduce(tensor_max_total_valid, op=torch.distributed.ReduceOp.MAX, group=subgroup)
        len_basis = tensor_max_total_valid.item()

    target_len_global = math.ceil(len_basis / (W * size_group)) * (W * size_group)
    target_len_local = target_len_global // size_group

    # Vectorized target ranges (Broadcast scalar to all rows)
    all_target_starts = torch.arange(size_group, device=device) * target_len_local  # (WS,)
    all_target_ends = all_target_starts + target_len_local  # (WS,)

    # 1.5 Vectorized Communication

    if has_input:
        output_local = torch.zeros(shape_leading_flat_local, target_len_local, dtype=local_input.dtype, device=device)
        ops = []
        recv_bufs = {}

        # We will save these indices for backward
        # send_indices_map[peer] = (rows, cols) into local_valid_data
        # recv_indices_map[peer] = (rows, cols) into output_local
        send_indices_map = {}
        recv_indices_map = {}

        for peer in range(size_group):
            # --- SEND Logic ---
            # Intersection: [my_start, my_end) AND [target_start, target_end)
            p_target_start = all_target_starts[peer]  # Scalar
            p_target_end = all_target_ends[peer]  # Scalar

            # Broadcast scalar target to (shape_leading_flat_local,)
            overlap_start = torch.maximum(my_global_start, p_target_start)
            overlap_end = torch.minimum(my_global_end, p_target_end)

            # Convert to local indices relative to local_valid_data (starts at 0)
            # local_s = overlap_start - my_global_start
            local_start = (overlap_start - my_global_start).clamp(min=0)
            local_end = (overlap_end - my_global_start).clamp(min=0)

            # (shape_leading_flat_local, 2)
            send_ranges = torch.stack([local_start, local_end], dim=1)
            send_rows, send_cols = get_flattened_range_indices(send_ranges, device)

            send_indices_map[peer] = (send_rows, send_cols)
            send_buf = local_valid_data[send_rows, send_cols]

            if peer != rank:
                if send_buf.numel() > 0:
                    ops.append(dist.P2POp(dist.isend, send_buf, global_rank_peers[peer], group=group))

            # --- RECV Logic ---
            # Intersection: [my_target_start, my_target_end) AND [src_start, src_end)
            # my_target_start is rank * target_len_local (scalar)
            my_t_start = rank * target_len_local
            my_t_end = (rank + 1) * target_len_local

            p_src_start = global_starts[peer]  # (shape_leading_flat_local,)
            p_src_end = global_ends[peer]  # (shape_leading_flat_local,)

            recv_overlap_start = torch.maximum(torch.tensor(my_t_start, device=device), p_src_start)
            recv_overlap_end = torch.minimum(torch.tensor(my_t_end, device=device), p_src_end)

            # Convert to local output indices (relative to my_t_start)
            out_start = (recv_overlap_start - my_t_start).clamp(min=0)
            out_end = (recv_overlap_end - my_t_start).clamp(min=0)

            # (shape_leading_flat_local, 2)
            recv_ranges = torch.stack([out_start, out_end], dim=1)
            recv_rows, recv_cols = get_flattened_range_indices(recv_ranges, device)

            recv_indices_map[peer] = (recv_rows, recv_cols)
            recv_len = recv_rows.numel()

            if peer != rank:
                if recv_len > 0:
                    recv_buf = torch.empty(recv_len, dtype=local_input.dtype, device=device)
                    ops.append(dist.P2POp(dist.irecv, recv_buf, global_rank_peers[peer], group=group))
                    recv_bufs[peer] = recv_buf
            else:
                # Self-copy
                output_local[recv_rows, recv_cols] = send_buf

        if ops:
            reqs = dist.batch_isend_irecv(ops)
            for req in reqs:
                req.wait()

        for peer, buf in recv_bufs.items():
            r_rows, r_cols = recv_indices_map[peer]
            output_local[r_rows, r_cols] = buf

        # 1.6 Reshape to Output Format
        if local_input_moved.ndim > 1:
            # output_local is (shape_leading_flat_local, target_len_local)
            # Reconstruct non-axis dimensions
            output_reshaped = output_local.unflatten(0, local_input_moved.shape[:-1])
            output_final_flat = output_reshaped.movedim(-1, axis)
        else:
            # output_local.shape is (shape_leading_flat_local=1, target_len_local)
            # since (shape_leading_flat_local,) is an temporary axis we added in this function, we need to squeeze it out
            output_final_flat = output_local.squeeze(0)
    else:
        output_final_flat = None

    # 1.7 Generate output masks
    # Compute valid range per row for this rank
    # total_valid is (shape_leading_flat_local,) - global valid count per row
    # rank is local rank in group
    # target_len_local is scalar - max length per rank
    # Formula: i_end_valid_local = min(total_valid, (rank + 1) * target_len_local) - rank * target_len_local
    # This tells us how many valid elements this rank "owns" in the global sequence, starting from 0.
    i_end_valid_local = (
        torch.minimum(total_valid, torch.tensor((rank + 1) * target_len_local, device=device)) - rank * target_len_local
    )
    i_end_valid_local = i_end_valid_local.clamp(min=0)  # (shape_leading_flat_local,)

    # Create mask (shape_leading_flat_local, target_len_local)
    idx_cols = torch.arange(target_len_local, device=device)  # (target_len_local,)
    mask_local_2d = idx_cols.unsqueeze(0) < i_end_valid_local.unsqueeze(
        1
    )  # (shape_leading_flat_local, target_len_local)

    # Reshape mask to match output_final structure
    if local_input_moved.ndim > 1:
        mask_reshaped = mask_local_2d.unflatten(0, local_input_moved.shape[:-1])
        mask_final_flat = mask_reshaped.movedim(-1, axis)
    else:
        mask_final_flat = mask_local_2d.squeeze(0)

    if output_final_flat is None:
        # has_input==False means mask input and mask output
        output_final_flat = mask_final_flat

    # compute global output shape
    # This function doesn't modify any other input axes except for "axis" and for inserting
    # an extra axis at position "axis + 1". The modification for "axis" is guaranteed to be
    # evenly sharded by the i_dim_device_mesh_shard_axis.
    shape_output = list(input.shape)
    shape_output[axis] = output_final_flat.shape[axis] * size_group
    shape_output = tuple(shape_output)

    strides_output = update_exhaustive_strides(output_final_flat.shape, output_final_flat.stride(), shape_output)

    output = DTensor.from_local(output_final_flat, device_mesh, placements, shape=shape_output, stride=strides_output)

    strides_mask_output = update_exhaustive_strides(mask_final_flat.shape, mask_final_flat.stride(), shape_output)
    mask_output = DTensor.from_local(
        mask_final_flat, device_mesh, placements, shape=shape_output, stride=strides_mask_output
    )

    return (output if has_input else None, mask_output, argsort_mask_flat_local, valid_counts_all_ranks)


def _distributed_unpad_and_unpack(
    input: DTensor,
    axis: int,
    argsort_mask_flat_local: torch.Tensor,
    valid_counts_all_ranks_unpacked: torch.Tensor,
    shape_input_expected: torch.Size | None = None,
    device_mesh_expected: DeviceMesh | None = None,
    placements_expected: tuple[Placement, ...] | None = None,
) -> DTensor:
    """Distributed unpad and unpack operation for DTensor inputs.

    Inverse of _distributed_pack_and_pad. Scatters elements from the packed/padded layout
    back to their original positions using the argsort indices from the forward pass.
    Communication is performed to redistribute elements across ranks.

    Args:
        input (DTensor): Packed DTensor of shape (..., target_len, ...) sharded along axis.
        axis (int): Dimension containing the packed sequence.
        argsort_mask_flat_local (torch.Tensor): Local argsort indices from _distributed_pack_and_pad,
            2D tensor of shape (shape_leading_flat_local, shape_axis_local).
        valid_counts_all_ranks_unpacked (torch.Tensor): Valid counts per rank from
            _distributed_pack_and_pad, shape (size_group, shape_leading_flat_local).
        shape_input_expected (torch.Size | None, optional): Expected input shape for validation.
        device_mesh_expected (DeviceMesh | None, optional): Expected device mesh for validation.
        placements_expected (tuple[Placement, ...] | None, optional): Expected placements for validation.

    Returns:
        DTensor: Unpacked DTensor with elements scattered back to their original positions.
            Shape is (..., original_len, ...) where original_len = argsort_mask_flat_local.shape[1] * size_group.
    """
    device = input.device

    # sanity checks
    if not isinstance(input, DTensor):
        raise TypeError(f"Expected DTensor, got {type(input)}")
    if device_mesh_expected is not None and input.device_mesh != device_mesh_expected:
        raise ValueError(f"input device_mesh mismatch: expected {device_mesh_expected}, got {input.device_mesh}")
    if placements_expected is not None and input.placements != placements_expected:
        raise ValueError(f"input placements mismatch: expected {placements_expected}, got {input.placements}")
    if shape_input_expected is not None and input.shape != shape_input_expected:
        raise ValueError(f"input shape mismatch: expected {shape_input_expected}, got {input.shape}")

    if not isinstance(argsort_mask_flat_local, torch.Tensor):
        raise TypeError(f"Expected torch.Tensor, got {type(argsort_mask_flat_local)}")

    if argsort_mask_flat_local.ndim != 2:
        raise ValueError(f"argsort_mask_flat_local must be a 2D tensor, got {argsort_mask_flat_local.ndim}D tensor")

    device_mesh = input.device_mesh
    placements = input.placements

    i_dim_device_mesh_shard_axis = None
    for i_dim_device_mesh, placement in enumerate(placements):
        if isinstance(placement, Partial):
            raise ValueError("Partial placements are not supported")
        elif isinstance(placement, Shard):
            if input.shape[placement.dim] % device_mesh.shape[i_dim_device_mesh] != 0:
                raise ValueError(
                    f"Uneven sharding tensor dimension {placement.dim} of size {input.shape[placement.dim]} "
                    f"along device mesh dimension {i_dim_device_mesh} of size {device_mesh.shape[i_dim_device_mesh]} is not supported"
                )
            if placement.dim == axis:
                i_dim_device_mesh_shard_axis = i_dim_device_mesh

    if i_dim_device_mesh_shard_axis is None:
        raise ValueError(f"input is not sharded along axis {axis}")

    size_group = device_mesh.size(i_dim_device_mesh_shard_axis)
    if valid_counts_all_ranks_unpacked.shape != (size_group, argsort_mask_flat_local.shape[0]):
        raise ValueError(
            f"valid_counts_all_ranks_unpacked.shape {valid_counts_all_ranks_unpacked.shape} != "
            f"(size_group, argsort_mask_flat_local.shape[0]) {size_group, argsort_mask_flat_local.shape[0]}"
        )

    rank = device_mesh.get_local_rank(i_dim_device_mesh_shard_axis)
    group = device_mesh.get_group(i_dim_device_mesh_shard_axis)
    global_rank_peers = torch.distributed.get_process_group_ranks(group)

    # 2. Reshape input to 2D (D, max_target_len)
    input_local = input.to_local()
    input_moved = input_local.movedim(axis, -1)
    input_2d = input_moved.reshape(-1, input_moved.shape[-1])

    # except for the potential difference in padding along the 'axis',
    # argsort_mask_flat_local.shape[0] == input_2d.shape[0]
    if argsort_mask_flat_local.shape[0] != input_2d.shape[0]:
        raise ValueError(
            f"argsort_mask_flat_local.shape[0] {argsort_mask_flat_local.shape[0]} != input_2d.shape[0] {input_2d.shape[0]}"
        )

    # 3. Reconstruct Masks (same as _distributed_pack_and_pad logic)
    # NOTE: indexing the input with these global_{start,end} indices
    # will automatically exclude the invalid elements' contributions to the output
    global_ends = valid_counts_all_ranks_unpacked.cumsum(dim=0)
    global_starts = global_ends - valid_counts_all_ranks_unpacked
    my_global_start = global_starts[rank]

    # The 'target' length from the pad_and_pack logic correspond to the
    # length of the input along 'axis' because this is the inverse of pad_and_pack
    target_len_local = input_local.shape[axis]
    all_target_starts = torch.arange(size_group, device=device) * target_len_local
    all_target_ends = all_target_starts + target_len_local

    # 4. Reverse Communication
    output_valid_2d = torch.zeros(argsort_mask_flat_local.shape, dtype=input_2d.dtype, device=device)
    local_valid_counts = valid_counts_all_ranks_unpacked[rank]

    ops = []
    recv_bufs = {}

    for peer in range(size_group):
        # --- REVERSE RECV (Corresponds to Forward SEND) ---
        # Reconstruct send_mask intervals from forward
        p_target_start = all_target_starts[peer]
        p_target_end = all_target_ends[peer]

        overlap_start = torch.maximum(my_global_start, p_target_start)
        # Reconstruct my_global_end
        my_global_end = my_global_start + local_valid_counts
        overlap_end = torch.minimum(my_global_end, p_target_end)

        local_start = (overlap_start - my_global_start).clamp(min=0)
        local_end = (overlap_end - my_global_start).clamp(min=0)

        send_ranges = torch.stack([local_start, local_end], dim=1)
        send_rows, send_cols = get_flattened_range_indices(send_ranges, device)

        expected_output_len = send_rows.numel()

        if peer != rank:
            if expected_output_len > 0:
                output_recv_buf = torch.empty(expected_output_len, dtype=input_2d.dtype, device=device)
                ops.append(dist.P2POp(dist.irecv, output_recv_buf, global_rank_peers[peer], group=group))
                # Store indices to scatter later
                recv_bufs[peer] = (send_rows, send_cols, output_recv_buf)

        # --- REVERSE SEND (Corresponds to Forward RECV) ---
        # Reconstruct recv_mask intervals from forward
        my_t_start = rank * target_len_local
        my_t_end = (rank + 1) * target_len_local

        p_src_start = global_starts[peer]
        p_src_end = global_ends[peer]

        recv_overlap_start = torch.maximum(torch.tensor(my_t_start, device=device), p_src_start)
        recv_overlap_end = torch.minimum(torch.tensor(my_t_end, device=device), p_src_end)

        out_start = (recv_overlap_start - my_t_start).clamp(min=0)
        out_end = (recv_overlap_end - my_t_start).clamp(min=0)

        recv_ranges = torch.stack([out_start, out_end], dim=1)
        recv_rows, recv_cols = get_flattened_range_indices(recv_ranges, device)

        input_to_send = input_2d[recv_rows, recv_cols]

        if peer != rank:
            if input_to_send.numel() > 0:
                ops.append(dist.P2POp(dist.isend, input_to_send, global_rank_peers[peer], group=group))

        if peer == rank:
            # Self-copy: scatter directly
            # output_valid_2d[send_rows, send_cols] = input_to_send
            output_valid_2d[send_rows, send_cols] = input_to_send

    if ops:
        reqs = dist.batch_isend_irecv(ops)
        for req in reqs:
            req.wait()

    for peer, (r, c, buf) in recv_bufs.items():
        output_valid_2d[r, c] = buf

    # 5. Scatter to Original Shape
    # output_valid_2d is sorted and of same shape as argsort_mask_flat_local
    # Scatter back to unsorted
    output_2d = torch.zeros_like(output_valid_2d)
    output_2d.scatter_(1, argsort_mask_flat_local, output_valid_2d)

    # Reshape D back to non-axis dims
    if input.ndim > 1:
        shape_no_axis = list(input_local.shape)
        shape_no_axis.pop(axis)
        output_reshaped = output_2d.unflatten(0, tuple(shape_no_axis))
        output_final = output_reshaped.movedim(-1, axis)
    else:
        # output_2d.shape is (D=1, argsort_mask_flat_local.shape[1])
        # since (D,) is an temporary axis we added in this function, we need to squeeze it out
        output_final = output_2d.squeeze(0)

    # 'axis' is guaranteed sharded along device mesh dimension i_dim_device_mesh_shard_axis
    shape_axis_output = output_final.shape[axis] * size_group
    shape_output = input.shape[:axis] + (shape_axis_output,) + input.shape[axis + 1 :]

    strides_output = update_exhaustive_strides(output_final.shape, output_final.stride(), shape_output)

    output = DTensor.from_local(output_final, device_mesh, placements, shape=shape_output, stride=strides_output)

    return output


class DistributedPackAndPad(torch.autograd.Function):
    """Autograd function for distributed pack and pad operation.

    Forward: Left-packs valid elements and pads to the right (end) to a multiple of (W * size_group).
    Backward: Unpacks and scatters gradients back to original positions.

    See _distributed_pack_and_pad for detailed documentation.
    """

    @staticmethod
    def forward(
        ctx, input: DTensor, mask: DTensor, axis: int, W: int, keep_input_padding: bool = False
    ) -> tuple[DTensor, DTensor]:
        output, mask_output, argsort_mask_flat_local, valid_counts_all_ranks = _distributed_pack_and_pad(
            input, mask, axis, W, keep_input_padding
        )

        # Context Saving
        ctx.save_for_backward(argsort_mask_flat_local, valid_counts_all_ranks)
        ctx.axis = axis
        ctx.device_mesh = input.device_mesh
        ctx.placements = input.placements
        ctx.shape_output_fwd = output.shape
        ctx.mark_non_differentiable(mask_output, argsort_mask_flat_local, valid_counts_all_ranks)

        return output, mask_output

    @staticmethod
    def backward(ctx, grad_output: DTensor, grad_mask: DTensor) -> tuple[DTensor, None, None, None, None]:
        # 1. Unpack
        (argsort_mask_flat_local, valid_counts_all_ranks) = ctx.saved_tensors

        grad_input = _distributed_unpad_and_unpack(
            grad_output,
            ctx.axis,
            argsort_mask_flat_local,
            valid_counts_all_ranks,
            shape_input_expected=ctx.shape_output_fwd,
            device_mesh_expected=ctx.device_mesh,
            placements_expected=ctx.placements,
        )

        return grad_input, None, None, None, None


class DistributedUnpadAndUnpack(torch.autograd.Function):
    """Autograd function for distributed unpad and unpack operation.

    Forward: Unpacks elements from packed/padded layout back to original positions.
    Backward: Packs gradients using the forward pass of _distributed_pack_and_pad.

    This is the inverse operation of DistributedPackAndPad.
    """

    @staticmethod
    def forward(
        ctx,
        input: DTensor,
        mask: DTensor,
        mask_original: DTensor,
        axis: int,
        keep_input_padding: bool,
    ) -> DTensor:
        # masks must have same ndim as input
        # Non-axis trailing dimensions can be 1 (for broadcasting) or match input
        if mask.ndim != input.ndim:
            raise RuntimeError(
                f"mask ndim {mask.ndim} must equal input ndim {input.ndim}. "
                f"For 3D inputs, use 3D masks with shape (B, N, 1) for broadcasting."
            )

        if mask_original.ndim != input.ndim:
            raise RuntimeError(
                f"mask_original ndim {mask_original.ndim} must equal input ndim {input.ndim}. "
                f"For 3D inputs, use 3D masks with shape (B, N, 1) for broadcasting."
            )

        # Masks must be boolean to ensure correct valid count computation
        if mask.dtype != torch.bool:
            raise TypeError(
                f"mask must have dtype torch.bool to avoid precision issues in valid count computation, "
                f"got {mask.dtype}. Use mask.bool() to convert."
            )
        if mask_original.dtype != torch.bool:
            raise TypeError(
                f"mask_original must have dtype torch.bool to avoid precision issues in valid count computation, "
                f"got {mask_original.dtype}. Use mask_original.bool() to convert."
            )

        # Check shapes are broadcast-compatible
        try:
            shape_input = torch.broadcast_shapes(input.shape, mask.shape)
        except RuntimeError as e:
            raise RuntimeError(f"Shapes of input {input.shape} and mask {mask.shape} are not broadcastable.") from e

        if shape_input != input.shape:
            raise RuntimeError(f"Broadcasted shape {shape_input} is not equal to input shape {input.shape}")

        if input.device_mesh != mask.device_mesh:
            raise RuntimeError(
                f"Input and mask must have the same device mesh but got {input.device_mesh} and {mask.device_mesh}"
            )

        if input.placements != mask.placements:
            raise RuntimeError(
                f"Input and mask must have the same placements but got {input.placements} and {mask.placements}"
            )

        if input.device_mesh != mask_original.device_mesh:
            raise RuntimeError(
                f"Input and mask must have the same device mesh but got {input.device_mesh} and {mask_original.device_mesh}"
            )

        if input.placements != mask_original.placements:
            raise RuntimeError(
                f"Input and mask must have the same placements but got {input.placements} and {mask_original.placements}"
            )

        placements = input.placements
        device_mesh = input.device_mesh

        i_dim_device_mesh_shard_axis = None
        for i_dim_device_mesh, placement in enumerate(placements):
            if isinstance(placement, Partial):
                raise ValueError("Partial placements are not supported")
            elif isinstance(placement, Shard):
                if input.shape[placement.dim] % device_mesh.shape[i_dim_device_mesh] != 0:
                    raise ValueError(
                        f"Uneven sharding tensor dimension {placement.dim} of size {input.shape[placement.dim]} "
                        f"along device mesh dimension {i_dim_device_mesh} of size {device_mesh.shape[i_dim_device_mesh]} is not supported"
                    )
                if placement.dim == axis:
                    i_dim_device_mesh_shard_axis = i_dim_device_mesh

        if i_dim_device_mesh_shard_axis is None:
            raise ValueError(f"input is not sharded along axis {axis}")

        # 0. Normalize axis
        ndim = input.ndim
        if axis < 0:
            axis += ndim

        if not (0 <= axis < ndim):
            raise ValueError(f"axis {axis} out of range for {ndim}D input")

        # 1. While this function does not concern about the 'window' size in application,
        # the usage requirement of _distributed_pack_and_pad requires input padding factor
        # Here we just take out the group_size from the input shape along 'axis' to get the padding factor
        W = input.shape[axis] // device_mesh.shape[i_dim_device_mesh_shard_axis]

        # 2. Expand mask_original in case the mask has been broadcasted from mask_original
        shape_mask_original_expanded = input.shape[:axis] + (mask_original.shape[axis],) + input.shape[axis + 1 :]

        try:
            shape_mask_original_expanded_broadcast = torch.broadcast_shapes(
                shape_mask_original_expanded, mask_original.shape
            )
        except RuntimeError as e:
            raise RuntimeError(
                f"Shapes of input {input.shape} and mask {mask_original.shape} "
                f"are not broadcastable excluding the sequence axis {axis}."
            ) from e

        if shape_mask_original_expanded_broadcast != shape_mask_original_expanded:
            raise ValueError(
                f"mask_original shape {mask_original.shape} is not broadcastable to input shape {input.shape} excluding the sequence axis {axis}"
            )

        if mask_original.shape != shape_mask_original_expanded:
            mask_original_local = mask_original.to_local()
            input_local = input.to_local()

            target_local_shape = (
                input_local.shape[:axis] + (mask_original_local.shape[axis],) + input_local.shape[axis + 1 :]
            )

            mask_original_local_expanded = mask_original_local.expand(target_local_shape)

            strides_mask_original_expanded = tuple(
                0 if mask_original_local_expanded.stride()[i] == 0 else mask_original.stride()[i]
                for i in range(mask_original.ndim)
            )

            mask_original = DTensor.from_local(
                mask_original_local_expanded,
                mask_original.device_mesh,
                mask_original.placements,
                shape=shape_mask_original_expanded,
                stride=strides_mask_original_expanded,
            )

        # 3. Re-compute metadata from mask_original
        # We run forward pass using mask_original as both input and mask.
        # This is cheap(er) and gives us the correct metadata for the backward pass.
        # We need mask_original to be a DTensor.

        _, mask_pack_and_pad, argsort_mask_flat_local, valid_counts_all_ranks = _distributed_pack_and_pad(
            None, mask_original, axis, W, keep_input_padding
        )

        ctx.mark_non_differentiable(mask_original, mask_pack_and_pad, argsort_mask_flat_local, valid_counts_all_ranks)

        # 4. Check mask consistency (if mask_qw_dtensor is provided)
        # We assume mask passed in corresponds to original packed and padded mask.
        if not torch.equal(mask.to_local(), mask_pack_and_pad.to_local()):
            raise ValueError("mask_original does not correspond to mask_qw_dtensor when sorted/packed.")

        # 5. Invert the data
        output = _distributed_unpad_and_unpack(input, axis, argsort_mask_flat_local, valid_counts_all_ranks)

        ctx.save_for_backward(mask_original)
        ctx.params = {"axis": axis, "W": W, "keep_input_padding": keep_input_padding}

        return output

    @staticmethod
    def backward(ctx, grad_output: DTensor) -> tuple[DTensor, None, None, None, None]:
        # grad_output is Square layout (gradient of our output).
        # We want gradient of our input (Window layout).
        # This corresponds to Forward pass of DistributedUnmaskReshape.

        (mask_original,) = ctx.saved_tensors
        params = ctx.params
        axis = params["axis"]
        W = params["W"]
        keep_input_padding = params["keep_input_padding"]

        grad_input, _, _, _ = _distributed_pack_and_pad(grad_output, mask_original, axis, W, keep_input_padding)

        assert grad_input is not None

        return grad_input, None, None, None, None


def distributed_pack_and_pad(
    input_dtensor: DTensor, mask_dtensor: DTensor, W: int, axis: int, keep_input_padding: bool = False
) -> tuple[DTensor, DTensor]:
    """Distributed left-pack and pad operation for DTensor inputs.

    Left-packs valid elements (as indicated by mask) to the front along the specified axis,
    then pads to the right (end) to the next multiple of (W * size_group) where size_group
    is the number of ranks sharding the tensor along axis.

    Args:
        input_dtensor (DTensor): Input DTensor of shape (..., seq_len, ...) sharded along axis.
        mask_dtensor (DTensor): Boolean mask DTensor with shape broadcastable to input.
            True indicates valid elements, False indicates padding to ignore.
        W (int): Padding factor. The output length along axis will be padded to the next
            multiple of (W * size_group).
        axis (int): Dimension containing the sequence to pack and pad.
        keep_input_padding (bool, optional): If False (default), output length is based on
            the maximum number of valid elements across all slices. If True, output length
            is based on input_dtensor.shape[axis]. Defaults to False.

    Returns:
        tuple[DTensor, DTensor]:
            output: Packed and padded DTensor. Shape is (..., target_len, ...) where
                target_len is a multiple of (W * size_group).
            mask_output: Boolean mask DTensor for output, same shape as output.
    """
    return DistributedPackAndPad.apply(input_dtensor, mask_dtensor, axis, W, keep_input_padding)


def distributed_unpad_and_unpack(
    input: DTensor,
    mask: DTensor,
    mask_original: DTensor,
    axis: int,
    keep_input_padding: bool,
) -> DTensor:
    """Distributed unpad and unpack operation for DTensor inputs.

    Inverse of distributed_pack_and_pad. Scatters elements from the packed/padded layout
    back to their original positions.

    Args:
        input (DTensor): Packed DTensor of shape (..., target_len, ...) sharded along axis.
        mask (DTensor): Boolean mask DTensor for input, indicating valid (True) vs
            invalid (False) elements in the packed layout.
        mask_original (DTensor): Boolean mask DTensor indicating valid elements in the
            original unpacked layout. Used to reconstruct the argsort indices.
        axis (int): Dimension containing the packed sequence.
        keep_input_padding (bool): Whether input padding was kept in the forward pass
            of distributed_pack_and_pad.

    Returns:
        DTensor: Unpacked DTensor with elements scattered back to their original positions.
            Shape matches mask_original's shape along axis.
    """
    return DistributedUnpadAndUnpack.apply(input, mask, mask_original, axis, keep_input_padding)


class DistributedGatherSlidingWindows(torch.autograd.Function):
    @staticmethod
    def forward(ctx, dense_dtensor: DTensor, window_size: int, axis: int) -> DTensor:
        """
        Distributed Forward Pass using ownership-based halo exchange.

        Args:
            dense_dtensor: Input DTensor sharded along axis
            window_size: h = H // (W//2) in the original window batching parameters from Boltz
            axis: dense_dtensor's axis to apply windowing
        """
        # 0. sanity checks
        if not isinstance(dense_dtensor, DTensor):
            raise TypeError(f"Expected DTensor, got {type(dense_dtensor)}")

        if not isinstance(axis, int):
            raise TypeError(f"Expected int for axis, got {type(axis)}")
        # Normalize axis
        ndim = dense_dtensor.ndim

        if ndim < 2:
            raise ValueError(f"dense_dtensor must have at least 2 dimensions, got {ndim}D")

        if axis < 0:
            axis += ndim
        if not (0 <= axis < ndim):
            raise ValueError(f"axis {axis} out of range for {ndim}D input")
        if not (isinstance(window_size, int) and window_size > 0 and window_size % 2 == 0):
            # h := window_size must be an even integer per the original get_indexing_matrix function from Boltz,
            # i.e., h / 2 must be an integer
            raise TypeError(f"Expected positive even integer for window_size, got {type(window_size)}")

        # the halo size computation from compute_query_window_ownership assumes the
        # window batching parameters W, H and K and their math relationship from the
        # original get_indexing_matrix function from Boltz, so the input dense_dtensor's
        # shape must satisfy the requirements:
        # 1. dense_dtensor.shape[axis] == 2 * K
        # 2. dense_dtensor.shape[axis + 1] == W // 2
        # In addition, the window_size must be even integer and that
        # window_size * dense_dtensor.shape[axis + 1] gives the resulting H
        if dense_dtensor.shape[axis] % 2 != 0:
            raise ValueError(f"dense_dtensor.shape[{axis}] must be even, got {dense_dtensor.shape[axis]}")
        K = dense_dtensor.shape[axis] // 2
        W = dense_dtensor.shape[axis + 1] * 2
        H = window_size * dense_dtensor.shape[axis + 1]

        device_mesh = dense_dtensor.device_mesh
        placements = dense_dtensor.placements

        i_dim_device_mesh_shard_axis = None
        i_axes_sharded_by_mesh_dim = [None] * device_mesh.ndim
        for i_dim_device_mesh, placement in enumerate(placements):
            if isinstance(placement, Partial):
                raise ValueError("Partial placements are not supported")
            elif isinstance(placement, Shard):
                if dense_dtensor.shape[placement.dim] % device_mesh.shape[i_dim_device_mesh] != 0:
                    raise ValueError(
                        f"Uneven sharding tensor dimension {placement.dim} of size {dense_dtensor.shape[placement.dim]} "
                        f"along device mesh dimension {i_dim_device_mesh} of size {device_mesh.shape[i_dim_device_mesh]} is not supported"
                    )
                if placement.dim == axis:
                    i_dim_device_mesh_shard_axis = i_dim_device_mesh
                if placement.dim == axis + 1:
                    # NOTE: technically, axis + 1 of shape (W // 2) is already supported because
                    # there is no special handling of axis + 1 inside this function but it's only
                    # treated as another axes orthogonal to axis. But the upstream DistributedUnmaskReshape
                    # can not produce sharded placement of axis + 1 so we exclude it from execution
                    # to fence off upstream bugs
                    raise NotImplementedError(f"Sharding along axis {axis + 1} is not supported")
                i_axes_sharded_by_mesh_dim[i_dim_device_mesh] = placement.dim

        if i_dim_device_mesh_shard_axis is None:
            raise ValueError(f"input dense_dtensor is not sharded along axis {axis}")

        # 1. Unpack DTensor and get rank info
        local_tensor = dense_dtensor.to_local()
        rank_in_group = device_mesh.get_local_rank(i_dim_device_mesh_shard_axis)
        size_group = device_mesh.size(i_dim_device_mesh_shard_axis)
        group = device_mesh.get_group(i_dim_device_mesh_shard_axis)
        ranks_global_group = torch.distributed.get_process_group_ranks(group)
        rank_global = ranks_global_group[rank_in_group]

        # 2. Determine ownership from DTensor sharding
        # The DTensor is sharded along axis (sequence of half-windows)
        # Global shape[axis] = 2*K, local_tensor.shape[axis] = (2 * K) / size_group
        # The following code, esp. the usage of compute_query_window_ownership,
        # assumes the underlying sequence length owned by each rank, i.e.,
        # local_tensor.shape[axis] * local_tensor.shape[axis + 1], is a multiple of "W".
        # This is equivalent to: (2 * K // size_group * (W // 2)) % W == 0, i.e.,
        # (K // size_group * W) % W == 0, i.e., K // size_group is an integer
        if K % size_group != 0:
            raise ValueError(
                f"K {K} is not a integer multiple of the number of ranks {size_group} sharding the dense_dtensor along axis {axis}"
            )
        local_hw_len = local_tensor.shape[axis]
        hw_start = rank_in_group * local_hw_len
        hw_end = (rank_in_group + 1) * local_hw_len

        # Query windows: each QW i owns half-windows [2i, 2i+1]
        # So if we own HW [hw_start, hw_end), we own QW [hw_start//2, hw_end//2)
        qw_start = hw_start // 2
        qw_end = hw_end // 2

        # Validate: each rank must own at least one query window
        assert qw_start < qw_end, (
            f"Rank {rank_global} has no query windows to process: QW[{qw_start},{qw_end}). "
            f"This typically means size_group > K. Either reduce size_group or increase K."
        )

        ownership = compute_query_window_ownership(W, H, K, qw_start, qw_end)

        hw_need_start, hw_need_end = ownership["hw_needed"]
        left_halo_size = ownership["left_halo_size"]
        right_halo_size = ownership["right_halo_size"]

        # 3. Halo Exchange (multi-hop)
        device = local_tensor.device

        # Create halo buffers
        left_halo = None
        right_halo = None

        if left_halo_size > 0:
            # Get shape with left_halo_size at axis
            left_shape = list(local_tensor.shape)
            left_shape[axis] = left_halo_size
            left_halo = torch.zeros(left_shape, dtype=local_tensor.dtype, device=device)

        if right_halo_size > 0:
            right_shape = list(local_tensor.shape)
            right_shape[axis] = right_halo_size
            right_halo = torch.zeros(right_shape, dtype=local_tensor.dtype, device=device)

        recv_meta, send_meta = get_halo_from_neighbors(rank_in_group, size_group, local_hw_len, W, H, K)

        ops = []

        # Execute Recvs
        recv_temps = []
        for peer, htype, offset, length in recv_meta:
            target_buffer = left_halo if htype == "left" else right_halo
            # Recv needs contiguous buffer. target_buffer slice might not be.
            # Create temp buffer
            shape = list(target_buffer.shape)
            shape[axis] = length
            recv_buf = torch.empty(shape, dtype=local_tensor.dtype, device=device)
            ops.append(dist.P2POp(dist.irecv, recv_buf, ranks_global_group[peer], group=group))
            recv_temps.append((recv_buf, target_buffer, offset, length))

        # Execute Sends
        for peer, offset, length in send_meta:
            send_buf = local_tensor.narrow(axis, offset, length).contiguous()
            ops.append(dist.P2POp(dist.isend, send_buf, ranks_global_group[peer], group=group))

        if ops:
            reqs = dist.batch_isend_irecv(ops)
            for req in reqs:
                req.wait()

        # Copy temp buffers to actual halo tensors
        for buf, target, off, ln in recv_temps:
            target.narrow(axis, off, ln).copy_(buf)

        # 4. Construct extended local view: [left_halo, local_data, right_halo]
        parts = []
        if left_halo is not None:
            parts.append(left_halo)
        parts.append(local_tensor)
        if right_halo is not None:
            parts.append(right_halo)

        extended_local = torch.cat(parts, dim=axis) if len(parts) > 1 else local_tensor

        # 5. Apply local unfold on extended data
        # Compute local offsets relative to extended_local
        # In extended_local, owned half-windows start at index left_halo_size
        # Needed half-windows span [hw_need_start, hw_need_end)
        # Offset in extended_local = (hw_index - hw_need_start)
        # This leverages the translational equivalence property of the gather_sliding_windows function:
        # T(x[j_s:j_e], offsets[k_s:k_e] - j_s) = T(x, offsets)[k_s:k_e]
        # where "T" is equivalent to the gather_sliding_windows function
        # (or equivalently the underlying Toeplitz matrix operation).

        # Generate offsets for owned query windows
        h = window_size
        offset_start = 1 - h // 2
        owned_qw_offsets_global = torch.arange(offset_start + 2 * qw_start, offset_start + 2 * qw_end, 2, device=device)

        # Convert to local indices in extended_local
        # Global offset points to a half-window index
        # In extended_local: HW[hw_need_start] is at index 0
        local_offsets = owned_qw_offsets_global - hw_need_start

        # Apply efficient unfold on extended data
        # local_result has shape (..., K_local, window_size, W/2, ...) where local_result.shape[axis] == K_local
        local_result = gather_sliding_windows(extended_local, local_offsets, window_size, axis)
        # This function requires the input dense_dtensor.shape[axis] is evenly sharded by the device_mesh
        # and guarantees local_result is evenly sharded along the same axis
        shape_output = list(local_result.shape)
        for i_dim_device_mesh, i_axis in enumerate(i_axes_sharded_by_mesh_dim):
            if i_axis is None:
                continue
            shape_output[i_axis] = shape_output[i_axis] * device_mesh.size(i_dim_device_mesh)
        shape_output = tuple(shape_output)
        strides_output = update_exhaustive_strides(local_result.shape, local_result.stride(), shape_output)

        # 6. Save context for backward
        ctx.save_for_backward(local_offsets)
        ctx.params = {
            "fwd_local_input_shape": local_tensor.shape,
            "fwd_input_shape": dense_dtensor.shape,
            "fwd_output_shape": shape_output,
            "ownership": ownership,
            "axis": axis,
            "mesh": device_mesh,
            "placements": placements,
            "i_dim_device_mesh_shard_axis": i_dim_device_mesh_shard_axis,
            "window_size": window_size,
            "recv_meta": recv_meta,
            "send_meta": send_meta,
        }

        # 7. Wrap in DTensor (sharded on query window dimension)
        return DTensor.from_local(
            local_result, device_mesh, placements, shape=torch.Size(shape_output), stride=tuple(strides_output)
        )

    @staticmethod
    def backward(ctx, grad_output: DTensor) -> tuple[DTensor, None, None]:
        """
        Distributed Backward Pass: Compute gradients with neighbor exchange.

        Steps:
        1. Compute gradient w.r.t. extended_local (includes halos)
        2. Split into [left_halo_grad, local_grad, right_halo_grad]
        3. Send halo grads back to neighbors who own that data
        4. Receive grads from neighbors who used our data as halos
        5. Accumulate and return
        """
        # 1. Unpack context
        (local_offsets,) = ctx.saved_tensors
        params = ctx.params
        fwd_local_input_shape = params["fwd_local_input_shape"]
        fwd_input_shape = params["fwd_input_shape"]
        fwd_output_shape = params["fwd_output_shape"]
        ownership = params["ownership"]
        axis = params["axis"]
        mesh = params["mesh"]
        placements = params["placements"]
        i_dim_device_mesh_shard_axis = params["i_dim_device_mesh_shard_axis"]
        window_size = params["window_size"]
        recv_meta = params["recv_meta"]
        send_meta = params["send_meta"]

        hw_start, hw_end = ownership["hw_owned"]
        hw_need_start, hw_need_end = ownership["hw_needed"]
        left_halo_size = ownership["left_halo_size"]
        right_halo_size = ownership["right_halo_size"]

        # sanity checks
        if not isinstance(grad_output, DTensor):
            raise TypeError(f"Expected DTensor, got {type(grad_output)}")
        if grad_output.device_mesh != mesh:
            raise ValueError(f"grad_output device_mesh mismatch: expected {mesh}, got {grad_output.device_mesh}")
        if grad_output.placements != placements:
            raise ValueError(f"grad_output placements mismatch: expected {placements}, got {grad_output.placements}")
        if grad_output.shape != fwd_output_shape:
            raise ValueError(f"grad_output shape mismatch: expected {fwd_output_shape}, got {grad_output.shape}")

        local_grad_output = grad_output.to_local()

        # 2. Compute gradient w.r.t. extended_local using standalone backward function
        # Build extended_local's shape (same shape as local_input but with extended_len at axis)
        extended_len = hw_need_end - hw_need_start
        extended_shape = list(fwd_local_input_shape)
        extended_shape[axis] = extended_len

        # Call standalone backward to get gradient w.r.t. extended_local
        grad_extended = gather_sliding_windows_backward(
            local_grad_output, local_offsets, window_size, axis, tuple(extended_shape)
        )

        # 3. Split extended gradient
        offset = 0
        grad_left_halo = grad_extended.narrow(axis, offset, left_halo_size) if left_halo_size > 0 else None
        offset += left_halo_size

        local_owned_len = hw_end - hw_start
        grad_local = grad_extended.narrow(axis, offset, local_owned_len)
        offset += local_owned_len

        grad_right_halo = grad_extended.narrow(axis, offset, right_halo_size) if right_halo_size > 0 else None

        # 4. Exchange halo gradients
        group = mesh.get_group(i_dim_device_mesh_shard_axis)
        ranks_global_group = torch.distributed.get_process_group_ranks(group)

        ops = []

        # 1. Send gradients for data I received in fwd (halo grads)
        # fwd recv_meta: (peer, type, offset, length)
        # I received 'length' from 'peer' into my 'type' halo at 'offset'.
        # Now I send that slice of grad back to 'peer'.
        for peer, h_type, offset, length in recv_meta:
            if h_type == "left":
                assert grad_left_halo is not None
                grad_chunk = grad_left_halo.narrow(axis, offset, length)
            else:
                assert grad_right_halo is not None
                grad_chunk = grad_right_halo.narrow(axis, offset, length)

            ops.append(dist.P2POp(dist.isend, grad_chunk.contiguous(), ranks_global_group[peer], group=group))

        # 2. Recv gradients for data I sent in fwd (accumulate to local)
        # fwd send_meta: (peer, offset, length)
        # I sent 'length' from my local at 'offset' to 'peer'.
        # Now I recv that grad from 'peer' and add to my local grad.
        recv_grads = []  # (tensor, offset, length)
        for peer, offset, length in send_meta:
            grad_buf = torch.empty(
                grad_local.shape[:axis] + (length,) + grad_local.shape[axis + 1 :],
                dtype=grad_local.dtype,
                device=grad_local.device,
            )
            ops.append(dist.P2POp(dist.irecv, grad_buf, ranks_global_group[peer], group=group))
            recv_grads.append((grad_buf, offset, length))

        if ops:
            reqs = dist.batch_isend_irecv(ops)
            for req in reqs:
                req.wait()

        # 5. Accumulate received gradients
        for grad_buf, offset, length in recv_grads:
            target = grad_local.narrow(axis, offset, length)
            target.add_(grad_buf)

        # 6. Return DTensor gradient
        strides_grad_input = update_exhaustive_strides(grad_local.shape, grad_local.stride(), fwd_input_shape)
        grad_input = DTensor.from_local(
            grad_local, mesh, placements, shape=torch.Size(fwd_input_shape), stride=tuple(strides_grad_input)
        )

        return grad_input, None, None


def distributed_gather_sliding_windows(dense_dtensor: DTensor, window_size: int, axis: int) -> DTensor:
    """
    Distributed version of gather_sliding_windows for DTensor inputs.

    Args:
        dense_dtensor: Input DTensor sharded along axis dimension
        window_size: h = H // (W//2)
        axis: Dimension to apply windowing

    Returns:
        DTensor sharded on query window dimension
    """
    return DistributedGatherSlidingWindows.apply(dense_dtensor, window_size, axis)


def convert_single_repr_to_window_batched_key(x: DTensor, W: int, H: int) -> DTensor:
    """Converts a single representation tensor to a window-batched key tensor.

    Reshapes and processes the input tensor to create overlapping windows suitable for
    attention keys in windowed attention mechanisms. The input is unflattened into
    half-windows and then gathered using sliding windows.

    Args:
        x: Input tensor of shape (B, N, ...), where B is batch size and N is sequence length.
        W: Query window size.
        H: Key window size.

    Returns:
        A DTensor of shape (B, K, H, ...), where K = N // W is the number of windows.

    Raises:
        TypeError: If ``x`` is not a DTensor.
        ValueError: If ``x`` has fewer than 2 dimensions.
        ValueError: If ``x.shape[1]`` is not divisible by ``W``.
    """
    # input is assumed to be in shape (B, N, ...)
    if not isinstance(x, DTensor):
        raise TypeError(f"x must be a DTensor, but got {type(x)}")
    if x.ndim < 2:
        raise ValueError(f"x must have at least 2 dimensions, but got x.ndim={x.ndim}")

    validate_window_batching_parameters(W, H, True)

    if x.shape[1] % W != 0:
        raise ValueError(f"x.shape[1] must be divisible by W, but got x.shape[1]={x.shape[1]} and W={W}")

    K = x.shape[1] // W
    h = H // (W // 2)

    # (B, K*W, D) -> (B, 2*K, W//2, D)
    x_unflat_hw = shardwise_unflatten_sharded(x, axis=1, sizes=(2 * K, W // 2))
    # (B, 2*K, W//2, D) -> (B, K, h, W // 2, D)
    x_unflat_key = distributed_gather_sliding_windows(x_unflat_hw, window_size=h, axis=1)
    # (B, K, h, W // 2, D) -> (B, K, H, D)
    x_key = shardwise_flatten(x_unflat_key, start_dim=2, end_dim=3)
    return x_key


def convert_single_repr_window_batched_query_to_key(x: DTensor, W: int, H: int) -> DTensor:
    """Converts a window-batched query tensor to a window-batched key tensor.

    First flattens the query-batched input and then converts it to a key tensor using
    ``convert_single_repr_to_window_batched_key``.

    Args:
        x: Input tensor of shape (B, K, W, ...), where B is batch size, K is number of windows,
            and W is query window size.
        W: Query window size.
        H: Key window size.

    Returns:
        A DTensor of shape (B, K, H, ...).

    Raises:
        TypeError: If ``x`` is not a DTensor.
        ValueError: If ``x`` has fewer than 3 dimensions.
        ValueError: If ``x.shape[2]`` is not equal to ``W``.
    """
    # input is assumed to be in shape (B, K, W, ...)
    if not isinstance(x, DTensor):
        raise TypeError(f"x must be a DTensor, but got {type(x)}")
    if x.ndim < 3:
        raise ValueError(f"x must have at least 3 dimensions, but got x.ndim={x.ndim}")

    if x.shape[2] != W:
        raise ValueError(f"x.shape[2] must be equal to W, but got x.shape[2]={x.shape[2]} and W={W}")

    validate_window_batching_parameters(W, H, True)

    # (B, K, W, ...) -> (B, K*W, ...)
    x_flat = shardwise_flatten_sharded(x, start_dim=1, end_dim=2)
    x_key = convert_single_repr_to_window_batched_key(x_flat, W, H)
    return x_key
