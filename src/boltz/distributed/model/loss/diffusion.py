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


import warnings

import torch
from einops import einsum
from torch.distributed.tensor import DTensor, Replicate, Shard

try:
    from boltz.distributed.model.loss.triton.smooth_lddt_loss import (
        grid_launch_config,
        smooth_lddt_loss_bwd_kernel,
        smooth_lddt_loss_fwd_kernel,
    )

    has_smooth_lddt_loss_triton_kernels = True
except ImportError:
    has_smooth_lddt_loss_triton_kernels = False


from boltz.distributed.comm import TransposeComm
from boltz.distributed.model.layers.clip import clip
from boltz.distributed.model.layers.elementwise_op import (
    ElementwiseOp,
    elementwise_op,
    scalar_tensor_op,
    single_tensor_op,
)
from boltz.distributed.model.layers.outer_op import OuterOp, replicate_to_shard_outer_op
from boltz.distributed.model.layers.redistribute_transpose import redistribute_transpose
from boltz.distributed.model.layers.repeat_interleave import shardwise_repeat_interleave
from boltz.distributed.model.layers.replicate_op import ReplicateOp, replicate_op
from boltz.distributed.model.layers.sharded_op import sharded_sum
from boltz.distributed.model.layers.where import where
from boltz.distributed.utils import LayoutRightMap


def weighted_rigid_align(
    true_coords: DTensor,
    pred_coords: DTensor,
    weights: DTensor,
    mask: DTensor,
) -> DTensor:
    """Compute weighted alignment and return the aligned true coordinates using DTensor.

    Implements the same algorithm as boltz.model.loss.diffusionv2.weighted_rigid_align
    (Algorithm 28). Supports only 3D inputs (B, N, 3) for coords and 2D (B, N) for
    weights/mask, with placements (Shard(0), Shard(1), Replicate()).
    SVD is computed in float64 for numerical stability.
    This function is NOT differentiable (uses torch.no_grad() internally).

    The computation is structured to ensure binary-identical results across
    cp_axis_1 column groups by computing on column-0 then broadcasting.
    SVD is computed on rank (:, 0, 0) then broadcast to avoid numerical
    divergence from parallel SVD.

    Parameters
    ----------
    true_coords : DTensor
        Ground truth atom coordinates, shape (B, N, 3).
        Placements: (Shard(0), Shard(1), Replicate()).
    pred_coords : DTensor
        Predicted atom coordinates, shape (B, N, 3).
        Placements: (Shard(0), Shard(1), Replicate()).
    weights : DTensor
        Alignment weights, shape (B, N).
        Placements: (Shard(0), Shard(1), Replicate()).
    mask : DTensor
        Atom mask, shape (B, N).
        Placements: (Shard(0), Shard(1), Replicate()).

    Returns
    -------
    DTensor
        Aligned true coordinates with same placements as input true_coords.

    """
    # Ndim checks (3D coords, 2D weights/mask)
    if true_coords.ndim != 3:
        raise ValueError(f"true_coords must be 3D (B, N, 3), got ndim={true_coords.ndim}")
    if pred_coords.ndim != 3:
        raise ValueError(f"pred_coords must be 3D (B, N, 3), got ndim={pred_coords.ndim}")
    if weights.ndim != 2:
        raise ValueError(f"weights must be 2D (B, N), got ndim={weights.ndim}")
    if mask.ndim != 2:
        raise ValueError(f"mask must be 2D (B, N), got ndim={mask.ndim}")

    # Shape checks
    if true_coords.shape != pred_coords.shape:
        raise ValueError(f"true_coords shape {true_coords.shape} != pred_coords shape {pred_coords.shape}")
    if weights.shape != mask.shape:
        raise ValueError(f"weights shape {weights.shape} != mask shape {mask.shape}")
    if weights.shape != true_coords.shape[:2]:
        raise ValueError(f"weights shape {weights.shape} != expected {true_coords.shape[:2]}")

    # Device mesh checks
    if true_coords.device_mesh != pred_coords.device_mesh:
        raise ValueError("true_coords and pred_coords must be on the same device_mesh")
    if true_coords.device_mesh != weights.device_mesh:
        raise ValueError("true_coords and weights must be on the same device_mesh")
    if true_coords.device_mesh != mask.device_mesh:
        raise ValueError("true_coords and mask must be on the same device_mesh")

    # Placement checks
    placements = (Shard(0), Shard(1), Replicate())
    if true_coords.placements != placements:
        raise ValueError(f"true_coords placements {true_coords.placements} != expected {placements}")
    if pred_coords.placements != placements:
        raise ValueError(f"pred_coords placements {pred_coords.placements} != expected {placements}")
    if weights.placements != placements:
        raise ValueError(f"weights placements {weights.placements} != expected {placements}")
    if mask.placements != placements:
        raise ValueError(f"mask placements {mask.placements} != expected {placements}")

    with torch.no_grad():
        # Convert to local tensors
        true_coords_local = true_coords.to_local()
        pred_coords_local = pred_coords.to_local()
        weights_local = weights.to_local()
        mask_local = mask.to_local()

        device_mesh = true_coords.device_mesh
        # mesh axis 1: reduce along atom dimension (cp_axis_0)
        # mesh axis 2: broadcast along coordinate dimension (cp_axis_1, where coords are replicated)
        group_reduce_atoms = device_mesh.get_group(1)
        group_broadcast = device_mesh.get_group(2)

        rank_coord = device_mesh.get_coordinate()
        assert rank_coord is not None

        batch_size, num_points, dim = true_coords_local.shape
        weights_expanded = (mask_local * weights_local).unsqueeze(-1)

        # Scalar degenerate check (all ranks; compatible with per-batch check below)
        total_num_points = num_points * device_mesh.shape[1]
        if total_num_points < (dim + 1):
            warnings.warn(
                "The size of one of the point clouds is <= dim+1. "
                "`WeightedRigidAlign` cannot return a unique rotation.",
                UserWarning,
                stacklevel=1,
            )

        is_first_column_group = rank_coord[-1] == 0
        rank_broadcast = torch.distributed.get_global_rank(group_broadcast, 0)

        if is_first_column_group:
            # Per-batch degenerate check first (same as diffusionv2: mask.sum(dim=-1) < (dim + 1))
            mask_count_global = mask_local.sum(dim=1).clone()
            torch.distributed.all_reduce(
                mask_count_global,
                op=torch.distributed.ReduceOp.SUM,
                group=group_reduce_atoms,
            )
            degenerate_batch_indices = torch.where(mask_count_global < (dim + 1))[0]
            if degenerate_batch_indices.numel() > 0:
                warnings.warn(
                    f"[rank_coord:{rank_coord}] "
                    "The size of one of the point clouds is <= dim+1. "
                    "`WeightedRigidAlign` cannot return a unique rotation. "
                    f"Batch indices (subset): {degenerate_batch_indices.tolist()}",
                    UserWarning,
                    stacklevel=1,
                )

            # Compute on column-0 then broadcast for binary-identical results
            # Overlapped async reductions for centroids (dim=1 = points, equiv. dim=-2 in serial v2)
            weights_sum_local = weights_expanded.sum(dim=1, keepdim=True)
            req_reduce_weights = torch.distributed.all_reduce(
                weights_sum_local, op=torch.distributed.ReduceOp.SUM, group=group_reduce_atoms, async_op=True
            )

            true_coords_weighted_sum_local = (true_coords_local * weights_expanded).sum(
                dim=1, keepdim=True
            )  # points dim
            req_reduce_true_coords = torch.distributed.all_reduce(
                true_coords_weighted_sum_local,
                op=torch.distributed.ReduceOp.SUM,
                group=group_reduce_atoms,
                async_op=True,
            )

            pred_coords_weighted_sum_local = (pred_coords_local * weights_expanded).sum(
                dim=1, keepdim=True
            )  # points dim
            req_reduce_pred_coords = torch.distributed.all_reduce(
                pred_coords_weighted_sum_local,
                op=torch.distributed.ReduceOp.SUM,
                group=group_reduce_atoms,
                async_op=True,
            )

            req_reduce_weights.wait()
            req_reduce_true_coords.wait()
            true_centroid = true_coords_weighted_sum_local / weights_sum_local

            req_reduce_pred_coords.wait()
            pred_centroid = pred_coords_weighted_sum_local / weights_sum_local
            torch.distributed.broadcast(true_centroid, rank_broadcast, group=group_broadcast)
            torch.distributed.broadcast(pred_centroid, rank_broadcast, group=group_broadcast)
        else:
            true_centroid = torch.empty_like(true_coords_local[:, 0:1, :])
            pred_centroid = torch.empty_like(pred_coords_local[:, 0:1, :])
            torch.distributed.broadcast(true_centroid, rank_broadcast, group=group_broadcast)
            torch.distributed.broadcast(pred_centroid, rank_broadcast, group=group_broadcast)

        # Center the coordinates
        true_coords_centered = true_coords_local - true_centroid
        pred_coords_centered = pred_coords_local - pred_centroid

        # Compute the weighted covariance matrix
        cov_matrix_local = einsum(
            weights_expanded * pred_coords_centered, true_coords_centered, "b n i, b n j -> b i j"
        )
        original_dtype = cov_matrix_local.dtype

        if is_first_column_group:
            # Reduce covariance matrix, compute SVD on rank (:, 0, 0), broadcast
            rank_reduce_cov_matrix = torch.distributed.get_global_rank(group_reduce_atoms, 0)
            torch.distributed.reduce(
                cov_matrix_local,
                op=torch.distributed.ReduceOp.SUM,
                dst=rank_reduce_cov_matrix,
                group=group_reduce_atoms,
            )
            if rank_coord[1] == 0:
                # SVD in float64 for numerical stability
                cov_matrix_64 = cov_matrix_local.to(dtype=torch.float64)
                U, S, V = torch.linalg.svd(cov_matrix_64, driver="gesvd" if cov_matrix_64.is_cuda else None)
                V = V.mH

                # Same logic as diffusionv2: scalar num_points check (v2 uses num_points, not per-batch mask)
                if (S.abs() <= 1e-15).any() and not (total_num_points < (dim + 1)):
                    warnings.warn(
                        f"[rank_coord:{rank_coord}] "
                        "Excessively low rank of "
                        "cross-correlation between aligned point clouds. "
                        "`WeightedRigidAlign` cannot return a unique rotation.",
                        UserWarning,
                        stacklevel=1,
                    )

                # Rotation matrix with proper determinant
                rot_matrix = torch.einsum("b i j, b k j -> b i k", U, V)
                F = torch.eye(dim, dtype=cov_matrix_64.dtype, device=cov_matrix_64.device)[None].repeat(
                    batch_size, 1, 1
                )
                F[:, -1, -1] = torch.det(rot_matrix)
                rot_matrix = einsum(U, F, V, "b i j, b j k, b l k -> b i l")
                rot_matrix = rot_matrix.to(dtype=original_dtype).contiguous()
                torch.distributed.broadcast(rot_matrix, rank_reduce_cov_matrix, group=group_reduce_atoms)
            else:
                rot_matrix = torch.empty(
                    (batch_size, dim, dim), dtype=original_dtype, device=true_coords_local.device
                ).contiguous()
                torch.distributed.broadcast(rot_matrix, rank_reduce_cov_matrix, group=group_reduce_atoms)
            # Broadcast within each row
            torch.distributed.broadcast(rot_matrix, rank_broadcast, group=group_broadcast)
        else:
            rot_matrix = torch.empty(
                (batch_size, dim, dim), dtype=original_dtype, device=true_coords_local.device
            ).contiguous()
            torch.distributed.broadcast(rot_matrix, rank_broadcast, group=group_broadcast)

        # Apply rotation and translation
        aligned_coords_local = einsum(true_coords_centered, rot_matrix, "b n i, b j i -> b n j") + pred_centroid

        # Convert back to DTensor
        aligned_coords = DTensor.from_local(
            aligned_coords_local,
            device_mesh=device_mesh,
            placements=true_coords.placements,
            shape=true_coords.shape,
            stride=true_coords.stride(),
        )

        return aligned_coords


def smooth_lddt_loss(
    pred_coords: DTensor,
    true_coords: DTensor,
    is_nucleotide: DTensor,
    coords_mask: DTensor,
    comm: TransposeComm,
    nucleic_acid_cutoff: float = 30.0,
    other_cutoff: float = 15.0,
    multiplicity: int = 1,
    v2: bool = True,
) -> DTensor:
    """Compute the smooth LDDT loss using DTensor.

    NOTE There is potential memory optimization in diffusionv2.py in Boltz2.

    Parameters
    ----------
    pred_coords: DTensor
        The predicted atom coordinates with placements (Shard(0), Shard(1), Replicate())
    true_coords: DTensor
        The ground truth atom coordinates with placements (Shard(0), Shard(1), Replicate())
    is_nucleotide: DTensor
        The weights for alignment with placements (Shard(0), Shard(1), Replicate())
    coords_mask: DTensor
        The atoms mask with placements (Shard(0), Shard(1), Replicate())
    comm: TransposeComm
        The communication object
    nucleic_acid_cutoff: float
        The cutoff for nucleic acid
    other_cutoff: float
        The cutoff for other atoms
    multiplicity: int
        The multiplicity of the atoms
    v2: bool
        Whether to use the v2 version of the smooth LDDT loss, where the denominator is added with 1e-5.

    Returns
    -------
    DTensor
        The smooth LDDT loss with placement (Replicate(), Replicate(), Replicate())

    """
    is_nucleotide = is_nucleotide.to(torch.bool)

    coords_mask_pairwise_section = redistribute_transpose(
        coords_mask,
        comm,
        output_placements=(Shard(0), Replicate(), Shard(1)),
        dim0=None,
        dim1=None,
    )
    true_dists = replicate_to_shard_outer_op(true_coords, OuterOp.CDIST, 1, comm)
    dtype = true_dists.dtype

    is_nucleotide = shardwise_repeat_interleave(
        is_nucleotide, multiplicity, 0
    )  # (batch_size * multiplicity, num_atoms)
    coords_mask = shardwise_repeat_interleave(coords_mask, multiplicity, 0)
    coords_mask_pairwise_section = shardwise_repeat_interleave(coords_mask_pairwise_section, multiplicity, 0)

    # broadcast is_nucleotide over the second cp axis
    # serial code:is_nucleotide.unsqueeze(-1).expand(-1, -1, is_nucleotide.shape[-1])
    if is_nucleotide.placements != (Shard(0), Shard(1), Replicate()):
        raise ValueError(
            f"is_nucleotide placements {is_nucleotide.placements} != expected (Shard(0), Shard(1), Replicate())"
        )

    # [B, N]
    is_nucleotide_local = is_nucleotide.to_local()
    # [B, N, N]
    is_nucleotide_pair_local = is_nucleotide_local.unsqueeze(-1).expand(-1, -1, is_nucleotide_local.shape[-1])
    shape_is_nucleotide_pair = (is_nucleotide.shape[0], is_nucleotide.shape[1], is_nucleotide.shape[1])
    # torch.Tensor.expand sets the expanded axes' stride to 0. See official doc:
    # https://docs.pytorch.org/docs/stable/generated/torch.Tensor.expand.html#torch-tensor-expand
    stride_is_nucleotide_pair = is_nucleotide.stride() + (0,)
    is_nucleotide_pair = DTensor.from_local(
        is_nucleotide_pair_local,
        device_mesh=is_nucleotide.device_mesh,
        placements=(Shard(0), Shard(1), Shard(2)),
        shape=shape_is_nucleotide_pair,
        stride=stride_is_nucleotide_pair,
    )

    mask = where(
        is_nucleotide_pair,
        scalar_tensor_op(
            nucleic_acid_cutoff,
            true_dists,
            ElementwiseOp.GT,
        ),
        scalar_tensor_op(
            other_cutoff,
            true_dists,
            ElementwiseOp.GT,
        ),
    )
    mask = mask.to(dtype=dtype)

    # Zero out the diagonal. If in CP mode, this means only diagonal ranks participate.
    local_num_samples, local_num_atoms = pred_coords.to_local().shape[:2]
    if comm.is_self_comm:
        diag_mask_local = 1 - torch.eye(local_num_atoms, device=pred_coords.device)
    else:
        diag_mask_local = torch.ones(local_num_atoms, local_num_atoms, device=pred_coords.device)
    diag_mask_local = diag_mask_local.unsqueeze(0).expand(local_num_samples, -1, -1)
    shape_diag_mask = (pred_coords.shape[0], pred_coords.shape[1], pred_coords.shape[1])
    # diag_mask is created from scratch in LayoutRight. The expanded leading axis
    # has stride 0 -- see official doc:
    # https://docs.pytorch.org/docs/stable/generated/torch.Tensor.expand.html#torch-tensor-expand
    stride_diag_mask = (0,) + LayoutRightMap(shape=shape_diag_mask[1:]).strides
    diag_mask = DTensor.from_local(
        diag_mask_local,
        device_mesh=mask.device_mesh,
        placements=mask.placements,
        shape=shape_diag_mask,
        stride=stride_diag_mask,
    )
    mask = elementwise_op(mask, diag_mask, ElementwiseOp.PROD)

    # Apply coordinate mask
    mask = replicate_op(mask, coords_mask, dim_to_unsqueeze_rhs=2, op=ReplicateOp.PROD)
    mask = replicate_op(mask, coords_mask_pairwise_section, dim_to_unsqueeze_rhs=1, op=ReplicateOp.PROD)

    # Compute distances between all pairs of atoms
    pred_dists = replicate_to_shard_outer_op(pred_coords, OuterOp.CDIST, 1, comm)

    dist_diff = single_tensor_op(
        elementwise_op(true_dists, pred_dists, ElementwiseOp.SUB),
        ElementwiseOp.ABS,
    )
    # Compute epsilon values
    eps = single_tensor_op(scalar_tensor_op(0.5, dist_diff, ElementwiseOp.SUB), ElementwiseOp.SIGMOID)
    for cutoff in (1.0, 2.0, 4.0):
        eps = elementwise_op(
            eps,
            single_tensor_op(
                scalar_tensor_op(cutoff, dist_diff, ElementwiseOp.SUB),
                ElementwiseOp.SIGMOID,
            ),
            ElementwiseOp.SUM,
        )
    eps = scalar_tensor_op(0.25, eps, ElementwiseOp.PROD)

    assert mask.requires_grad is False
    num = sharded_sum(elementwise_op(eps, mask, ElementwiseOp.PROD), dim=(-1, -2))
    den = sharded_sum(mask, dim=(-1, -2))  # mask have no gradient. thus no need to use torch.no_grad()
    if v2:
        den = scalar_tensor_op(1e-5, den, ElementwiseOp.SUM)
    else:
        den = clip(den, min_val=1.0, max_val=None)

    lddt = elementwise_op(
        num,
        den,
        ElementwiseOp.DIV,
    )
    lddt = scalar_tensor_op(1.0, lddt, ElementwiseOp.SUB)

    lddt = scalar_tensor_op(
        1 / lddt.shape[0],
        sharded_sum(lddt, dim=0),
        ElementwiseOp.PROD,
    )  # mean along DP axis; placements: (Replicate(), Replicate(), Replicate())

    return lddt


def _smooth_lddt_loss_forward_local(
    pred_coords_local: torch.Tensor,
    true_coords_local: torch.Tensor,
    pred_coords_t_local: torch.Tensor,
    true_coords_t_local: torch.Tensor,
    is_nucleotide_local: torch.Tensor,
    coords_mask_local: torch.Tensor,
    coords_mask_t_local: torch.Tensor,
    is_self_comm: bool,
    nucleic_acid_cutoff: float,
    other_cutoff: float,
    multiplicity: int,
) -> tuple[torch.Tensor, torch.Tensor]:
    """Local computation for smooth LDDT loss.

    This function computes the local numerator and denominator for the smooth LDDT loss
    on shardwise tensors.
    """
    #### Compute forward pass locally to get pair masks ####
    dtype = true_coords_local.dtype

    # Expand is_nucleotide and masks locally
    # (B, N) -> (B * multiplicity, N)
    is_nucleotide_local = is_nucleotide_local.repeat_interleave(multiplicity, dim=0)
    coords_mask_local = coords_mask_local.repeat_interleave(multiplicity, dim=0)
    coords_mask_t_local = coords_mask_t_local.repeat_interleave(multiplicity, dim=0)

    # Construct pairwise nucleotide mask
    num_cols = pred_coords_t_local.shape[1]
    is_nucleotide_pair_local = is_nucleotide_local.unsqueeze(-1).expand(-1, -1, num_cols)  # O(N^2) tensor

    # Compute true distances
    true_dists = torch.cdist(true_coords_local, true_coords_t_local)  # O(N^2) tensor

    # Compute mask based on cutoffs
    mask = torch.where(
        is_nucleotide_pair_local.bool(),
        (true_dists < nucleic_acid_cutoff),
        (true_dists < other_cutoff),
    ).to(dtype=dtype)  # O(N^2) tensor

    # Zero out diagonal
    local_num_samples, local_num_atoms = pred_coords_local.shape[:2]
    if is_self_comm:
        diag_mask_local = 1 - torch.eye(local_num_atoms, num_cols, device=pred_coords_local.device)
    else:
        diag_mask_local = torch.ones(local_num_atoms, num_cols, device=pred_coords_local.device)

    diag_mask_local = diag_mask_local.unsqueeze(0).expand(local_num_samples, -1, -1)  # O(N^2) tensor
    mask = mask * diag_mask_local

    # Apply coordinate mask
    mask = mask * coords_mask_local.unsqueeze(-1)
    if coords_mask_t_local.ndim == 3:
        mask = mask * coords_mask_t_local.transpose(1, 2)
    else:
        mask = mask * coords_mask_t_local.unsqueeze(1)

    #### Compute forward pass ####
    # Compute predicted distances
    pred_dists = torch.cdist(pred_coords_local, pred_coords_t_local)  # O(N^2) tensor
    dist_diff = (true_dists - pred_dists).abs()  # O(N^2) tensor

    # Compute epsilon: O(N^2) tensors
    eps = torch.sigmoid(0.5 - dist_diff)
    for cutoff in (1.0, 2.0, 4.0):
        eps = eps + torch.sigmoid(cutoff - dist_diff)
    eps *= 0.25

    # Compute numerators and denominators
    # Sum over local atoms (rows and cols)
    num_local = (eps * mask).sum(dim=(1, 2))
    den_local = mask.sum(dim=(1, 2))

    return num_local, den_local


def _smooth_lddt_loss_backward_local(
    grad_num_reduced,
    grad_den_reduced,
    pred_coords_local,
    true_coords_local,
    pred_coords_t_local,
    true_coords_t_local,
    is_nucleotide_local,
    coords_mask_local,
    coords_mask_t_local,
    is_self_comm,
    nucleic_acid_cutoff,
    other_cutoff,
    multiplicity,
):
    #### Recompute forward pass locally to get pair masks ####
    dtype = true_coords_local.dtype

    # Expand is_nucleotide and masks locally
    # (B, N) -> (B * multiplicity, N)
    is_nucleotide_local = is_nucleotide_local.repeat_interleave(multiplicity, dim=0)
    coords_mask_local = coords_mask_local.repeat_interleave(multiplicity, dim=0)
    coords_mask_t_local = coords_mask_t_local.repeat_interleave(multiplicity, dim=0)

    # Construct pairwise nucleotide mask
    num_cols = pred_coords_t_local.shape[1]
    is_nucleotide_pair_local = is_nucleotide_local.unsqueeze(-1).expand(-1, -1, num_cols)  # O(N^2) tensor

    # Compute true distances
    true_dists = torch.cdist(true_coords_local, true_coords_t_local)  # O(N^2) tensor

    # Compute mask based on cutoffs
    mask = torch.where(
        is_nucleotide_pair_local.bool(),
        (true_dists < nucleic_acid_cutoff),
        (true_dists < other_cutoff),
    ).to(dtype=dtype)  # O(N^2) tensor

    # Zero out diagonal
    local_num_samples, local_num_atoms = pred_coords_local.shape[:2]
    if is_self_comm:
        diag_mask_local = 1 - torch.eye(local_num_atoms, num_cols, device=pred_coords_local.device)
    else:
        diag_mask_local = torch.ones(local_num_atoms, num_cols, device=pred_coords_local.device)

    diag_mask_local = diag_mask_local.unsqueeze(0).expand(local_num_samples, -1, -1)  # O(N^2) tensor
    mask = mask * diag_mask_local

    # Apply coordinate mask
    mask = mask * coords_mask_local.unsqueeze(-1)
    if coords_mask_t_local.ndim == 3:
        mask = mask * coords_mask_t_local.transpose(1, 2)
    else:
        mask = mask * coords_mask_t_local.unsqueeze(1)

    #### Compute backward pass ####
    # Compute pred diffs and dists
    # diff_vec = P_row - P_col
    diff_vec = pred_coords_local.unsqueeze(2) - pred_coords_t_local.unsqueeze(1)  # O(N^2) tensor
    pred_dists = diff_vec.norm(dim=-1)  # O(N^2) tensor
    dist_diff = (true_dists - pred_dists).abs()  # O(N^2) tensor

    # Compute d_eps_d_diff: O(N^2) tensors
    d_eps_d_diff = torch.zeros_like(dist_diff)
    for cutoff in (0.5, 1.0, 2.0, 4.0):
        val = cutoff - dist_diff
        sig = torch.sigmoid(val)
        d_eps_d_diff -= sig * (1 - sig)
    d_eps_d_diff *= 0.25

    # Compute d_L_d_pred_dists
    # d(diff)/d(pred) = sign(pred - true)
    grad_num_broadcast = grad_num_reduced.view_as(pred_dists[:, 0, 0]).view(-1, 1, 1)
    # O(N^2) tensor
    d_L_d_pred_dists = grad_num_broadcast * mask * d_eps_d_diff * torch.sign(pred_dists - true_dists)

    # Compute gradients w.r.t coords
    # safe normalization
    pred_dists_safe = pred_dists.unsqueeze(-1) + 1e-8  # O(N^2) tensor
    diff_dir = diff_vec / pred_dists_safe  # O(N^2) tensor

    d_L_d_diff_vec = d_L_d_pred_dists.unsqueeze(-1) * diff_dir  # O(N^2) tensor

    grad_pred_local = d_L_d_diff_vec.sum(dim=2)
    grad_pred_t_local = -d_L_d_diff_vec.sum(dim=1)

    return grad_pred_local, grad_pred_t_local


def _smooth_lddt_loss_local_triton_forward(
    pred_coords_local: torch.Tensor,
    true_coords_local: torch.Tensor,
    pred_coords_t_local: torch.Tensor,
    true_coords_t_local: torch.Tensor,
    is_nucleotide_local: torch.Tensor,
    coords_mask_local: torch.Tensor,
    coords_mask_t_local: torch.Tensor,
    is_self_comm: bool,
    nucleic_acid_cutoff: float,
    other_cutoff: float,
    multiplicity: int,
) -> tuple[torch.Tensor, torch.Tensor]:
    """Triton-based local computation for smooth LDDT loss forward pass."""
    if not has_smooth_lddt_loss_triton_kernels:
        raise ImportError("Smooth LDDT loss Triton kernels are not available.")

    if pred_coords_local.dtype == torch.bfloat16 or pred_coords_local.dtype == torch.float16:
        raise ValueError(
            f"Triton kernel for smooth LDDT loss does not support {pred_coords_local.dtype} "
            "due to precision issues. Please use float32."
        )

    # Handle coords_mask_t shape
    if coords_mask_t_local.ndim == 3:
        coords_mask_t_local = (
            coords_mask_t_local.squeeze(-1) if coords_mask_t_local.shape[-1] == 1 else coords_mask_t_local.squeeze(1)
        )

    assert pred_coords_local.shape == pred_coords_t_local.shape
    assert pred_coords_local.shape == true_coords_local.shape
    assert pred_coords_t_local.shape == true_coords_t_local.shape

    # different multiplicity input can share the same mask and is_nucleotide
    # input so we only require B * M % B == 0 here, with pred_coords_local.shape[0] = B * M
    # and is_nucleotide_local.shape[0] = B
    assert pred_coords_local.shape[0] % is_nucleotide_local.shape[0] == 0
    assert pred_coords_local.shape[0] % coords_mask_local.shape[0] == 0
    assert pred_coords_local.shape[0] % coords_mask_t_local.shape[0] == 0

    assert is_nucleotide_local.shape[1] == pred_coords_local.shape[1]
    assert coords_mask_local.shape[1] == pred_coords_local.shape[1]
    assert coords_mask_t_local.shape[1] == pred_coords_local.shape[1]

    # TODO Ensure inputs are contiguous where needed or handle strides
    # Note: Triton handles strides, so we don't strictly need contiguous, but it helps performance

    B = pred_coords_local.shape[0]

    num_output = torch.zeros(B, device=pred_coords_local.device, dtype=pred_coords_local.dtype)
    den_output = torch.zeros(B, device=pred_coords_local.device, dtype=pred_coords_local.dtype)

    # Cast inputs
    is_nucleotide_local = is_nucleotide_local.to(dtype=torch.int8)  # Use int8 for bool

    # Define grid lambda that autotune will use with the selected BLOCK size

    smooth_lddt_loss_fwd_kernel[grid_launch_config](
        pred_coords_local,
        true_coords_local,
        pred_coords_t_local,
        true_coords_t_local,
        is_nucleotide_local,
        coords_mask_local,
        coords_mask_t_local,
        num_output,
        den_output,
        pred_coords_local.stride(0),
        pred_coords_local.stride(1),
        pred_coords_local.stride(2),
        true_coords_local.stride(0),
        true_coords_local.stride(1),
        true_coords_local.stride(2),
        pred_coords_t_local.stride(0),
        pred_coords_t_local.stride(1),
        pred_coords_t_local.stride(2),
        true_coords_t_local.stride(0),
        true_coords_t_local.stride(1),
        true_coords_t_local.stride(2),
        is_nucleotide_local.stride(0),
        is_nucleotide_local.stride(1),
        coords_mask_local.stride(0),
        coords_mask_local.stride(1),
        coords_mask_t_local.stride(0),
        coords_mask_t_local.stride(1),
        nucleic_acid_cutoff,
        other_cutoff,
        is_self_comm,
        pred_coords_local.shape[0],
        pred_coords_local.shape[1],
        coords_mask_local.shape[0],
    )

    return num_output, den_output


def _smooth_lddt_loss_local_triton_backward(
    grad_num_reduced: torch.Tensor,
    grad_den_reduced: torch.Tensor,
    pred_coords_local: torch.Tensor,
    true_coords_local: torch.Tensor,
    pred_coords_t_local: torch.Tensor,
    true_coords_t_local: torch.Tensor,
    is_nucleotide_local: torch.Tensor,
    coords_mask_local: torch.Tensor,
    coords_mask_t_local: torch.Tensor,
    is_self_comm: bool,
    nucleic_acid_cutoff: float,
    other_cutoff: float,
    multiplicity: int,
) -> tuple[torch.Tensor, torch.Tensor]:
    """Triton-based local computation for smooth LDDT loss backward pass."""
    if not has_smooth_lddt_loss_triton_kernels:
        raise ImportError("Smooth LDDT loss Triton kernels are not available.")

    if pred_coords_local.dtype == torch.bfloat16 or pred_coords_local.dtype == torch.float16:
        raise ValueError(
            f"Triton kernel for smooth LDDT loss does not support {pred_coords_local.dtype} "
            "due to precision issues. Please use float32."
        )

    if coords_mask_t_local.ndim == 3:
        coords_mask_t_local = (
            coords_mask_t_local.squeeze(-1) if coords_mask_t_local.shape[-1] == 1 else coords_mask_t_local.squeeze(1)
        )

    assert pred_coords_local.shape == pred_coords_t_local.shape
    assert pred_coords_local.shape == true_coords_local.shape
    assert pred_coords_t_local.shape == true_coords_t_local.shape

    assert pred_coords_local.shape[0] == grad_num_reduced.shape[0]
    assert pred_coords_local.shape[0] == grad_den_reduced.shape[0]

    # different multiplicity input can share the same mask and is_nucleotide
    # input so we only require B * M % B == 0 here, with pred_coords_local.shape[0] = B * M
    # and is_nucleotide_local.shape[0] = B
    assert pred_coords_local.shape[0] % is_nucleotide_local.shape[0] == 0
    assert pred_coords_local.shape[0] % coords_mask_local.shape[0] == 0
    assert pred_coords_local.shape[0] % coords_mask_t_local.shape[0] == 0

    assert is_nucleotide_local.shape[1] == pred_coords_local.shape[1]
    assert coords_mask_local.shape[1] == pred_coords_local.shape[1]
    assert coords_mask_t_local.shape[1] == pred_coords_local.shape[1]

    # Output gradients
    grad_pred_local = torch.zeros_like(pred_coords_local)
    grad_pred_t_local = torch.zeros_like(pred_coords_t_local)

    is_nucleotide_local = is_nucleotide_local.to(dtype=torch.int8)

    smooth_lddt_loss_bwd_kernel[grid_launch_config](
        grad_num_reduced,
        grad_den_reduced,
        pred_coords_local,
        true_coords_local,
        pred_coords_t_local,
        true_coords_t_local,
        is_nucleotide_local,
        coords_mask_local,
        coords_mask_t_local,
        grad_pred_local,
        grad_pred_t_local,
        pred_coords_local.stride(0),
        pred_coords_local.stride(1),
        pred_coords_local.stride(2),
        true_coords_local.stride(0),
        true_coords_local.stride(1),
        true_coords_local.stride(2),
        pred_coords_t_local.stride(0),
        pred_coords_t_local.stride(1),
        pred_coords_t_local.stride(2),
        true_coords_t_local.stride(0),
        true_coords_t_local.stride(1),
        true_coords_t_local.stride(2),
        is_nucleotide_local.stride(0),
        is_nucleotide_local.stride(1),
        coords_mask_local.stride(0),
        coords_mask_local.stride(1),
        coords_mask_t_local.stride(0),
        coords_mask_t_local.stride(1),
        grad_pred_local.stride(0),
        grad_pred_local.stride(1),
        grad_pred_local.stride(2),
        grad_pred_t_local.stride(0),
        grad_pred_t_local.stride(1),
        grad_pred_t_local.stride(2),
        nucleic_acid_cutoff,
        other_cutoff,
        is_self_comm,
        pred_coords_local.shape[0],
        pred_coords_local.shape[1],
        coords_mask_local.shape[0],
    )

    return grad_pred_local, grad_pred_t_local


class SmoothLDDTLossTritonFunction(torch.autograd.Function):
    @staticmethod
    def forward(
        ctx,
        pred_coords: DTensor,
        true_coords: DTensor,
        is_nucleotide: DTensor,
        coords_mask: DTensor,
        comm: TransposeComm,
        nucleic_acid_cutoff: float = 30.0,
        other_cutoff: float = 15.0,
        multiplicity: int = 1,
        v2: bool = True,
    ) -> DTensor:
        # 1. Comm Block: Transpose inputs
        # Output placement (Shard(0), Replicate(), Shard(1)) puts the atom dim on Mesh 2 (Cols)
        # This matches the requirement for pairwise operations where one operand is on Mesh 1 (Rows)
        # and the other on Mesh 2 (Cols).
        target_placements = (Shard(0), Replicate(), Shard(1))

        true_coords_t = redistribute_transpose(true_coords, comm, target_placements, dim0=None, dim1=None)
        pred_coords_t = redistribute_transpose(pred_coords, comm, target_placements, dim0=None, dim1=None)
        coords_mask_t = redistribute_transpose(coords_mask, comm, target_placements, dim0=None, dim1=None)

        # 2. Comp Block: Local computation with autograd
        # We use torch.enable_grad() to allow gradients to flow through the local computation
        # which will be used in backward pass.
        pred_coords_local = pred_coords.to_local()
        true_coords_local = true_coords.to_local()
        is_nucleotide_local = is_nucleotide.to_local()
        coords_mask_local = coords_mask.to_local()

        pred_coords_t_local = pred_coords_t.to_local()
        true_coords_t_local = true_coords_t.to_local()
        coords_mask_t_local = coords_mask_t.to_local()

        num_local, den_local = _smooth_lddt_loss_local_triton_forward(
            pred_coords_local,
            true_coords_local,
            pred_coords_t_local,
            true_coords_t_local,
            is_nucleotide_local,
            coords_mask_local,
            coords_mask_t_local,
            bool(comm.is_self_comm),
            nucleic_acid_cutoff,
            other_cutoff,
            multiplicity,
        )

        # 3. Comm Block: Reduction
        # Reduce over atom groups (Mesh 1 and Mesh 2)
        device_mesh = pred_coords.device_mesh
        group_row = device_mesh.get_group(1)
        group_col = device_mesh.get_group(2)

        # All reduce num and den
        # Note: We can pack them into one tensor for fewer comms
        metrics_local = torch.stack([num_local, den_local])

        # Reduce over row group
        torch.distributed.all_reduce(metrics_local, op=torch.distributed.ReduceOp.SUM, group=group_row)
        # Reduce over col group
        torch.distributed.all_reduce(metrics_local, op=torch.distributed.ReduceOp.SUM, group=group_col)

        num_reduced, den_reduced = metrics_local[0], metrics_local[1]
        if v2:
            den_reduced = den_reduced + 1e-5
        else:
            den_reduced = torch.clamp(den_reduced, min=1.0, max=None)

        # Compute LDDT per sample
        lddt_per_sample = num_reduced / den_reduced
        lddt_per_sample = 1.0 - lddt_per_sample

        # Final average over batch (Mesh 0) if needed
        # The original code does: lddt = scalar_tensor_op(1/N, sharded_sum(lddt, dim=0), ...)
        # Here lddt_per_sample is (B*mult,).
        # We need to sum over batch dimension locally and reduce over Mesh 0.

        # Note: We save tensors for backward
        ctx.comm = comm
        ctx.nucleic_acid_cutoff = nucleic_acid_cutoff
        ctx.other_cutoff = other_cutoff
        ctx.multiplicity = multiplicity
        ctx.is_self_comm = bool(comm.is_self_comm)
        ctx.device_mesh = device_mesh
        ctx.global_batch_size = pred_coords.shape[0]
        ctx.pred_coords_shape = pred_coords.shape
        ctx.pred_coords_stride = pred_coords.stride()

        ctx.save_for_backward(
            pred_coords_local,
            true_coords_local,
            pred_coords_t_local,
            true_coords_t_local,
            is_nucleotide_local,
            coords_mask_local,
            coords_mask_t_local,
            num_reduced,
            den_reduced,
        )

        # Final aggregation
        lddt_sum = lddt_per_sample.sum()
        # Reduce over batch group (Mesh 0)
        group_batch = device_mesh.get_group(0)
        torch.distributed.all_reduce(lddt_sum, op=torch.distributed.ReduceOp.SUM, group=group_batch)

        lddt_final = lddt_sum / pred_coords.shape[0]

        # Return as DTensor (Replicate, Replicate, Replicate)
        # It's a scalar wrapped in DTensor
        return DTensor.from_local(
            lddt_final,
            device_mesh,
            (Replicate(), Replicate(), Replicate()),
            shape=torch.Size(()),
            stride=(),
        )

    @staticmethod
    def backward(ctx, grad_output):
        # grad_output is w.r.t lddt_final (scalar)
        (
            pred_coords_local,
            true_coords_local,
            pred_coords_t_local,
            true_coords_t_local,
            is_nucleotide_local,
            coords_mask_local,
            coords_mask_t_local,
            num_reduced,
            den_reduced,
        ) = ctx.saved_tensors

        # Backprop logic
        # lddt = 1 - num / den
        # Loss L = lddt_final.
        # dL/d(num_reduced) = dL/d(lddt_final) * d(lddt_final)/d(lddt_per_sample) * d(lddt_per_sample)/d(num_reduced)
        # d(lddt_final)/d(lddt_per_sample) = 1 / GlobalBatchSize
        # d(lddt_per_sample)/d(num_reduced) = -1 / den_reduced
        # d(lddt_per_sample)/d(den_reduced) = num_reduced / den_reduced^2

        grad_output_local = grad_output.to_local().item()
        scale = grad_output_local / ctx.global_batch_size

        inv_den = (
            1.0 / den_reduced
        )  # torch.clamp(den_reduced, min=1.0, max=None) or (den_reduced + 1e-5) is done in forward pass
        grad_num_reduced = scale * (-inv_den)
        grad_den_reduced = scale * (num_reduced * inv_den**2)

        grad_pred_local, grad_pred_t_local = _smooth_lddt_loss_local_triton_backward(
            grad_num_reduced,
            grad_den_reduced,
            pred_coords_local,
            true_coords_local,
            pred_coords_t_local,
            true_coords_t_local,
            is_nucleotide_local,
            coords_mask_local,
            coords_mask_t_local,
            ctx.is_self_comm,
            ctx.nucleic_acid_cutoff,
            ctx.other_cutoff,
            ctx.multiplicity,
        )

        # Accumulate partial gradients (reduce over missing dimensions)
        # grad_pred_local is partial sum over Cols (Dim 2)
        # grad_pred_t_local is partial sum over Rows (Dim 1)
        group_row = ctx.device_mesh.get_group(1)
        group_col = ctx.device_mesh.get_group(2)
        torch.distributed.all_reduce(grad_pred_local, op=torch.distributed.ReduceOp.SUM, group=group_col)
        torch.distributed.all_reduce(grad_pred_t_local, op=torch.distributed.ReduceOp.SUM, group=group_row)

        # Comm Block: Reverse transpose
        # We need to move grad_pred_t_local back to pred_coords layout
        # pred_coords_t was (S0, R, S1). grad has same layout.
        # We want (S0, S1, R).

        grad_pred_t_local_transposed = ctx.comm.enqueue_to_dispatch(grad_pred_t_local.contiguous())
        ctx.comm.wait_until_finished()

        # Sum gradients
        total_grad_pred_local = grad_pred_local + grad_pred_t_local_transposed

        # Return None for true_coords (no grad), etc.
        return (
            DTensor.from_local(
                total_grad_pred_local,
                ctx.device_mesh,
                (Shard(0), Shard(1), Replicate()),
                shape=ctx.pred_coords_shape,
                stride=ctx.pred_coords_stride,
            ),  # pred_coords
            None,  # true_coords
            None,  # is_nucleotide
            None,  # coords_mask
            None,  # comm
            None,  # nucleic_acid_cutoff
            None,  # other_cutoff
            None,  # multiplicity
            None,  # v2
        )


def smooth_lddt_loss_triton(
    pred_coords: DTensor,
    true_coords: DTensor,
    is_nucleotide: DTensor,
    coords_mask: DTensor,
    comm: TransposeComm,
    nucleic_acid_cutoff: float = 30.0,
    other_cutoff: float = 15.0,
    multiplicity: int = 1,
    v2: bool = True,
) -> DTensor:
    """Compute the smooth LDDT loss using DTensor and Triton kernel."""
    if not has_smooth_lddt_loss_triton_kernels:
        raise ImportError("Smooth LDDT loss Triton kernels are not available.")

    return SmoothLDDTLossTritonFunction.apply(
        pred_coords,
        true_coords,
        is_nucleotide,
        coords_mask,
        comm,
        nucleic_acid_cutoff,
        other_cutoff,
        multiplicity,
        v2,
    )
