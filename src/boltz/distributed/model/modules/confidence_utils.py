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


"""DTensor-compatible confidence utility functions.

This module provides DTensor implementations of confidence metric computations
used in the Boltz confidence module. Some helpers are shardwise, while others
redistribute inputs to replicated placements for global reductions.
"""

import torch
import torch.nn.functional as F
from torch.autograd.function import FunctionCtx
from torch.distributed.tensor import DTensor, Partial, Replicate, Shard

from boltz.data import const
from boltz.distributed.comm import TransposeComm
from boltz.distributed.model.layers.outer_op import OuterOp, distributed_outer_op
from boltz.distributed.model.loss.confidencev2 import compute_frame_pred
from boltz.distributed.utils import LayoutRightMap, update_exhaustive_strides
from boltz.model.modules.confidence_utils import (
    tm_function as serial_tm_function,
)

# Sentinel value for chain_pair_iptm entries where the chain pair does not exist
# on this rank's batch.  Valid iPTM values are in [0, 1], so -1.0 is unambiguous.
CHAIN_IPTM_SENTINEL = -1.0

# Small constant added to denominators to avoid division by zero in TM/iPTM-style metrics.
_EPS = 1e-5


class _ComputeAggregatedMetricImpl(torch.autograd.Function):
    """Autograd function for computing aggregated metric from logits.

    This implements the forward and backward passes for converting binned logits
    to expected metric values via softmax and weighted sum. The computation is
    shardwise (no communication required) since it operates along the replicated
    bins dimension.

    The metric computation follows:
        probs = softmax(logits, dim=-1)
        metric = sum(probs * bounds, dim=-1)

    Where bounds are bin centers computed as:
        bounds[i] = (i + 0.5) * bin_width, for i in [0, num_bins)
        bin_width = end / num_bins

    The backward pass uses PyTorch autograd on the local computation graph,
    avoiding the need for manual gradient derivation.

    See Also
    --------
    compute_aggregated_metric : The public API function that calls this.
    """

    @staticmethod
    @torch.amp.custom_fwd(device_type="cuda")
    def forward(
        ctx: FunctionCtx,
        logits: DTensor,
        end: float,
    ) -> DTensor:
        """Forward pass for computing aggregated metric.

        Parameters
        ----------
        ctx : FunctionCtx
            The autograd context object for saving tensors for backward.
        logits : DTensor
            Input logits tensor with bins as the last dimension.
            Typical shapes:
            - pLDDT: (B*mult, N_token, 50) with placements (Shard(0), Shard(1),Replicate())
            - pDE/pAE: (B*mult, N_token, N_token, 64) with placements (Shard(0), Shard(1), Shard(2))
        end : float
            Maximum value of the metric range. Default 1.0 for pLDDT, 32.0 for pAE.

        Returns
        -------
        DTensor
            Output metric tensor with the last dimension reduced.
            Shape is logits.shape[:-1].

        Raises
        ------
        TypeError
            If logits is not a DTensor.
        ValueError
            If Partial placements are present or the last dimension is sharded.
        """
        # Type checking
        if not isinstance(logits, DTensor):
            raise TypeError(f"Expected DTensor for logits, got {type(logits)}")

        device_mesh = logits.device_mesh
        placements = logits.placements

        # Validate placements
        last_dim = len(logits.shape) - 1
        for i_dim_device_mesh, placement in enumerate(placements):
            if isinstance(placement, Partial):
                raise ValueError("Partial placements are not supported")
            elif isinstance(placement, Shard):
                # Check that the last dimension (bins) is not sharded
                if placement.dim == last_dim:
                    raise ValueError(
                        f"The bins dimension (dim={last_dim}) must not be sharded for compute_aggregated_metric"
                    )
                # Check that sharded dimensions are evenly divided
                if logits.shape[placement.dim] % device_mesh.shape[i_dim_device_mesh] != 0:
                    raise ValueError(
                        f"Uneven sharding of tensor dimension {placement.dim} of size {logits.shape[placement.dim]} "
                        f"along device mesh dimension {i_dim_device_mesh} of size {device_mesh.shape[i_dim_device_mesh]} is not supported"
                    )

        # Detach and set requires_grad to build a local computation graph
        logits_local_orig = logits.to_local().detach().requires_grad_(logits.requires_grad)

        with torch.enable_grad():
            # Promote to at least float32 to match serial compute_aggregated_metric
            # which uses default-dtype (float32) torch.arange for bounds and
            # sum(probs * bounds) — element-wise ops that stay float32 under autocast.
            compute_dtype = torch.promote_types(logits_local_orig.dtype, torch.float32)
            num_bins = logits_local_orig.shape[-1]
            bin_width = end / num_bins
            bounds = (torch.arange(num_bins, device=logits_local_orig.device, dtype=compute_dtype) + 0.5) * bin_width

            probs = F.softmax(logits_local_orig.to(compute_dtype), dim=-1)

            # Use sum(probs * bounds) instead of matmul(probs, bounds) to match
            # serial code. Under autocast, matmul is on the "lower precision"
            # list and would downcast float32 probs to BF16, while element-wise
            # multiply and sum are not affected by autocast.
            metric_local = torch.sum(
                probs * bounds.view(*((1,) * (probs.ndim - 1)), num_bins),
                dim=-1,
            )

        # Compute output shape and stride (remove last dimension)
        output_shape = tuple(logits.shape[:-1])
        output_stride = LayoutRightMap(output_shape).strides

        # Save tensors for backward pass
        ctx.save_for_backward(logits_local_orig, metric_local)
        ctx.device_mesh = device_mesh
        ctx.placements = placements
        ctx.logits_shape = logits.shape
        ctx.logits_stride = logits.stride()
        ctx.output_shape = output_shape
        ctx.output_stride = output_stride

        # Create output DTensor
        result = DTensor.from_local(
            metric_local.detach(),
            device_mesh=device_mesh,
            placements=placements,
            shape=output_shape,
            stride=output_stride,
        )

        return result

    @staticmethod
    @torch.amp.custom_bwd(device_type="cuda")
    def backward(
        ctx: FunctionCtx,
        grad_output: DTensor,
    ) -> tuple[DTensor | None, None]:
        """Backward pass for computing aggregated metric.

        Computes gradients by backpropagating through the local computation graph
        that was built during the forward pass. This leverages PyTorch's autograd
        rather than manual gradient computation.

        Parameters
        ----------
        ctx : FunctionCtx
            The autograd context containing saved tensors from forward:
            - logits_local_orig: Input logits tensor (local)
            - metric_local: Output metric tensor that holds the computation graph
            - device_mesh, placements: DTensor metadata
            - logits_shape, logits_stride: Shape/stride info for DTensor reconstruction
        grad_output : DTensor
            Gradient of loss with respect to output metric.

        Returns
        -------
        tuple[DTensor | None, None]
            Gradients for each forward input in order:
            - d_logits: DTensor or None, gradient for logits
            - None: end parameter (non-differentiable)

        Notes
        -----
        The gradient computation follows the chain rule:

        Forward:
            probs = softmax(logits, dim=-1)
            metric = sum(probs * bounds, dim=-1)

        Backward (using autograd):
            d_logits = d_metric @ d_metric/d_probs @ d_probs/d_logits
                     = d_metric @ bounds @ softmax_jacobian

        Where softmax_jacobian for dim=-1 is block-diagonal per position,
        making this a purely local operation.
        """
        logits_local, metric_local = ctx.saved_tensors

        if not logits_local.requires_grad:
            return None, None

        grad_output_local = grad_output.to_local()

        # Backprop via the local graph
        (d_logits_local,) = torch.autograd.grad(
            outputs=[metric_local],
            inputs=[logits_local],
            grad_outputs=[grad_output_local],
            retain_graph=False,  # Frees the local graph immediately
        )

        # Wrap gradient in DTensor
        d_logits = DTensor.from_local(
            d_logits_local,
            device_mesh=ctx.device_mesh,
            placements=ctx.placements,
            shape=ctx.logits_shape,
            stride=ctx.logits_stride,
        )

        return d_logits, None


def compute_aggregated_metric(logits: DTensor, end: float = 1.0) -> DTensor:
    """Compute the metric from logits via softmax and weighted sum.

    This is the DTensor-compatible version of the serial compute_aggregated_metric
    function. It converts binned logits to expected metric values by computing
    the softmax-weighted average of bin centers.

    All operations are shardwise (no inter-rank communication required) since
    they operate along the replicated bins dimension (last dimension).

    Parameters
    ----------
    logits : DTensor
        The input logits tensor with bins as the last dimension.
        Typical shapes and placements:
        - pLDDT: (B*mult, N_token, 50) with placements (Shard(0), Replicate())
        - pDE/pAE: (B*mult, N_token, N_token, 64) with placements (Shard(0), Shard(1), Replicate())
    end : float, optional
        Maximum value of the metric range, by default 1.0.
        Use 1.0 for pLDDT, 32.0 for pAE.

    Returns
    -------
    DTensor
        The computed metric tensor with shape logits.shape[:-1].
        Placements are preserved from input.

    Examples
    --------
    >>> # pLDDT computation
    >>> plddt_logits = ...  # DTensor with shape (B*mult, N_token, 50)
    >>> plddt = compute_aggregated_metric(plddt_logits, end=1.0)
    >>> # plddt has shape (B*mult, N_token)

    >>> # pAE computation
    >>> pae_logits = ...  # DTensor with shape (B*mult, N_token, N_token, 64)
    >>> pae = compute_aggregated_metric(pae_logits, end=32.0)
    >>> # pae has shape (B*mult, N_token, N_token)

    See Also
    --------
    boltz.model.modules.confidence_utils.compute_aggregated_metric : Serial version.
    """
    return _ComputeAggregatedMetricImpl.apply(logits, end)


class _LocalShardedSum(torch.autograd.Function):
    """Sum over a sharded dimension with all-reduce within the placement group.

    Forward: local sum over reduced_dim, then all_reduce(SUM) over the process
    group that shards that dimension. Backward: gradient is replicated along
    reduced_dim (evenly to all shards).

    Parameters
    ----------
    x_local : torch.Tensor
        Local shard of the DTensor (i.e. ``dtensor.to_local()``). Must have
        the same shape and layout as the local piece implied by the global
        DTensor shape and input_placements on this rank.
    reduced_dim : int
        Dimension to reduce (0-based). Must be one of the dimensions that has
        ``Shard(d)`` in input_placements; the all-reduce runs over the process
        group for that placement.
    input_placements : tuple[object, ...]
        Placements of the original DTensor (e.g. ``(Shard(0), Shard(1), Replicate())``).
        Used to select which mesh dimension to all-reduce and to iterate
        ``device_mesh.get_all_groups()``.
    device_mesh : DeviceMesh
        The DTensor device mesh. Must match the mesh used to shard the tensor
        that produced x_local.

    Returns
    -------
    torch.Tensor
        Local shard of the sum result. The dimension reduced_dim is removed;
        on ranks that shard that dimension, the local chunk is the same (replicated).
    """

    @staticmethod
    def forward(  # type: ignore[override]
        ctx,
        x_local: torch.Tensor,
        reduced_dim: int,
        input_placements: tuple[object, ...],
        device_mesh,
    ) -> torch.Tensor:
        """Sum over reduced_dim and all-reduce across ranks that shard it."""
        output_local = torch.sum(x_local, dim=reduced_dim, keepdim=False)
        for placement, placement_group in zip(input_placements, device_mesh.get_all_groups()):
            if isinstance(placement, Shard) and placement.dim == reduced_dim:
                torch.distributed.all_reduce(
                    output_local,
                    op=torch.distributed.ReduceOp.SUM,
                    group=placement_group,
                )
        ctx.input_local_shape = x_local.shape
        ctx.reduced_dim = reduced_dim
        return output_local

    @staticmethod
    def backward(ctx, grad_output_local: torch.Tensor) -> tuple[torch.Tensor | None, None, None, None]:
        """Replicate grad_output along reduced_dim for upstream."""
        dx_local = grad_output_local.unsqueeze(ctx.reduced_dim)
        dx_local = dx_local.expand(ctx.input_local_shape).clone(memory_format=torch.contiguous_format)
        return dx_local, None, None, None


class _LocalShardedMax(torch.autograd.Function):
    """Max over a sharded dimension with all-reduce within the placement group.

    Forward: local max over reduced_dim, then all_reduce(MAX) over the process
    group that shards that dimension. Backward: gradient flows only to the
    local argmax elements.

    Parameters
    ----------
    x_local : torch.Tensor
        Local shard of the DTensor. Same semantics as _LocalShardedSum: must
        be the local piece implied by the global shape and input_placements.
    reduced_dim : int
        Dimension to reduce. Must correspond to a Shard in input_placements.
    input_placements : tuple[object, ...]
        Placements of the original DTensor; used to find the group for all_reduce(MAX).
    device_mesh : DeviceMesh
        Device mesh for the DTensor.

    Returns
    -------
    torch.Tensor
        Local shard of the max result (reduced_dim removed). After all-reduce,
        all ranks in the group that shard reduced_dim hold the same values.
    """

    @staticmethod
    def forward(  # type: ignore[override]
        ctx,
        x_local: torch.Tensor,
        reduced_dim: int,
        input_placements: tuple[object, ...],
        device_mesh,
    ) -> torch.Tensor:
        """Max over reduced_dim and all-reduce across ranks that shard it."""
        output_local_keepdim = torch.amax(x_local, dim=reduced_dim, keepdim=True)
        for placement, placement_group in zip(input_placements, device_mesh.get_all_groups()):
            if isinstance(placement, Shard) and placement.dim == reduced_dim:
                torch.distributed.all_reduce(
                    output_local_keepdim,
                    op=torch.distributed.ReduceOp.MAX,
                    group=placement_group,
                )
        ctx.reduced_dim = reduced_dim
        ctx.save_for_backward(x_local, output_local_keepdim)
        return output_local_keepdim.squeeze(reduced_dim)

    @staticmethod
    def backward(ctx, grad_output_local: torch.Tensor) -> tuple[torch.Tensor | None, None, None, None]:
        """Route gradient only to elements that matched the max."""
        x_local, output_local_keepdim = ctx.saved_tensors
        grad_output_local = grad_output_local.unsqueeze(ctx.reduced_dim)
        mask = x_local == output_local_keepdim
        dx_local = grad_output_local * mask
        dx_local = dx_local.expand(x_local.shape).clone(memory_format=torch.contiguous_format)
        return dx_local, None, None, None


def _reduced_placements(
    input_placements: tuple[object, ...], input_shape: torch.Size, reduced_dim: int
) -> tuple[object, ...]:
    """Placements for a tensor after reducing (removing) one dimension.

    Used to build the output DTensor layout when a reduction (e.g. sum or max)
    over reduced_dim is performed: the reduced dimension is squeezed out, so
    placements must be updated accordingly.

    Parameters
    ----------
    input_placements : tuple[object, ...]
        Placements of the tensor before reduction (e.g. (Shard(0), Shard(1), Replicate())).
    input_shape : torch.Size
        Shape of the tensor before reduction. Used to get ndim and to map
        placement dim indices after removing reduced_dim.
    reduced_dim : int
        The dimension that was reduced (0-based). This dimension is removed
        in the output.

    Returns
    -------
    tuple[object, ...]
        Placements for the reduced tensor. The placement that was Shard(reduced_dim)
        becomes Replicate(); any Shard(d) with d > reduced_dim becomes
        Shard(d - 1); Replicate() is unchanged.
    """
    ndim = len(input_shape)
    shift = torch.zeros(ndim, dtype=torch.int64)
    shift[reduced_dim] = 1
    map_dims = (torch.arange(ndim, dtype=torch.int64) - shift.cumsum(0)).tolist()
    output_placements = []
    for placement in input_placements:
        if isinstance(placement, Shard) and placement.dim == reduced_dim:
            output_placements.append(Replicate())
        elif isinstance(placement, Shard):
            output_placements.append(Shard(map_dims[placement.dim]))
        elif isinstance(placement, Replicate):
            output_placements.append(placement)
    return tuple(output_placements)


def _reduced_shape_stride(
    input_shape: torch.Size, input_stride: tuple[int, ...], reduced_dim: int
) -> tuple[torch.Size, tuple[int, ...]]:
    """Shape and stride for a tensor after removing (squeezing) one dimension.

    Used together with _reduced_placements to construct the correct shape and
    stride for DTensor.from_local() when building a reduced output (e.g. after
    sum/max over reduced_dim).

    Parameters
    ----------
    input_shape : torch.Size
        Shape before reduction.
    input_stride : tuple[int, ...]
        Stride of the tensor before reduction (must match input_shape layout).
    reduced_dim : int
        Dimension that was reduced and is to be removed (0-based).

    Returns
    -------
    tuple[torch.Size, tuple[int, ...]]
        output_shape: input_shape with reduced_dim removed (length ndim - 1).
        output_stride: strides with the reduced_dim entry removed, consistent
        with the new shape.
    """
    shape_output = list(input_shape)
    shape_output[reduced_dim] = 1
    shape_output = tuple(shape_output)
    strides_output = update_exhaustive_strides(input_shape, input_stride, shape_output)
    shape_output = tuple(dim for i, dim in enumerate(shape_output) if i != reduced_dim)
    strides_output = tuple(stride for i, stride in enumerate(strides_output) if i != reduced_dim)
    return torch.Size(shape_output), strides_output


class _ComputePtmsImpl(torch.autograd.Function):
    """Fused pTM/ipTM computation from PAE logits and token/chain masks.

    Computes pTM, ipTM, ligand ipTM, protein ipTM, and per-chain-pair ipTM
    using distributed outer ops (TransposeComm) and local sharded sum/max
    (_LocalShardedSum, _LocalShardedMax). Aggregation is performed only within
    each CP group so that results match serial semantics when the serial
    reference is run on the same DP chunk.

    See Also
    --------
    compute_ptms : Public API that builds masks and calls this.
    """

    @staticmethod
    @torch.amp.custom_fwd(device_type="cuda")
    def forward(  # type: ignore[override]
        ctx: FunctionCtx,
        mask_collinear_pred: DTensor,
        token_pad_mask: DTensor,
        asym_id_base: DTensor,
        mol_type: DTensor,
        logits: DTensor,
        multiplicity: int,
        transpose_comm: TransposeComm,
    ) -> tuple[DTensor, DTensor, DTensor, DTensor, dict[int, dict[int, DTensor]]]:
        """Compute pTM, ipTM, ligand/protein ipTM, and chain_pair_iptm from local shards.

        Parameters
        ----------
        mask_collinear_pred : DTensor
            Mask of shape (B, mult, N_token) or (B*mult, N_token), True where
            the token is valid for PTM (non-collinear frame). Typically
            (Shard(0), Shard(1), Replicate()). Must be on the same device_mesh
            as other inputs.
        token_pad_mask : DTensor
            Token padding mask, shape (B, N_token). True = valid token.
            Placements (Shard(0), Shard(1), Replicate()). Repeated internally
            by multiplicity for the computation.
        asym_id_base : DTensor
            Chain ID per token, shape (B, N_token). Used to build inter-chain
            pair masks for ipTM and chain_pair_iptm. Same placements as
            token_pad_mask. Unique values are gathered within cp_axis_0 to form
            chain_pair_iptm keys.
        mol_type : DTensor
            Token molecule type (e.g. PROTEIN, NONPOLYMER), shape (B, N_token).
            Used for ligand_iptm and protein_iptm masks. Same mesh and
            placement convention as token_pad_mask.
        logits : DTensor
            PAE logits, shape (B*mult, N_token, N_token, num_bins). Placements
            (Shard(0), Shard(1), Shard(2)). Bins dimension must be replicated.
        multiplicity : int
            Number of diffusion samples per batch element. Batch dimension
            of logits is B * multiplicity.
        transpose_comm : TransposeComm
            Communication helper for the CP subgroup; used by distributed_outer_op
            to build full pair masks from sharded rows/columns.

        Returns
        -------
        tuple[DTensor, DTensor, DTensor, DTensor, dict[int, dict[int, DTensor]]]
            ptm, iptm, ligand_iptm, protein_iptm (each shape (B,) DTensors), and
            chain_pair_iptm mapping (idx1, idx2) -> DTensor of shape (B,).
            Dict keys are the union of chain IDs across all ranks (world-level),
            so all DP ranks have identical key sets. Entries where a chain pair
            does not exist on this DP rank's batch are filled with sentinel
            value -1.0.

        Requirements
        ------------
        - All DTensor inputs must use the same device_mesh.
        - token_pad_mask, asym_id_base, mol_type must have the same placements
          (typically (Shard(0), Shard(1), Replicate())).
        - logits placements must have the last dimension (bins) replicated.

        Raises
        ------
        TypeError
            If any of mask_collinear_pred, token_pad_mask, asym_id_base,
            mol_type, or logits is not a DTensor.
        ValueError
            If any input's device_mesh differs from token_pad_mask.device_mesh,
            or if transpose_comm is invalid for the current mesh.
        """
        if not isinstance(mask_collinear_pred, DTensor):
            raise TypeError(f"Expected DTensor for mask_collinear_pred, got {type(mask_collinear_pred)}")
        if not isinstance(token_pad_mask, DTensor):
            raise TypeError(f"Expected DTensor for token_pad_mask, got {type(token_pad_mask)}")
        if not isinstance(asym_id_base, DTensor):
            raise TypeError(f"Expected DTensor for asym_id, got {type(asym_id_base)}")
        if not isinstance(mol_type, DTensor):
            raise TypeError(f"Expected DTensor for mol_type, got {type(mol_type)}")
        if not isinstance(logits, DTensor):
            raise TypeError(f"Expected DTensor for logits, got {type(logits)}")

        device_mesh = token_pad_mask.device_mesh
        if mask_collinear_pred.device_mesh != device_mesh:
            raise ValueError(
                "mask_collinear_pred must be on the same device mesh as token_pad_mask, "
                f"got {mask_collinear_pred.device_mesh} and {device_mesh}"
            )
        if asym_id_base.device_mesh != device_mesh:
            raise ValueError(
                f"asym_id must be on the same device mesh as token_pad_mask, got {asym_id_base.device_mesh} and {device_mesh}"
            )
        if mol_type.device_mesh != device_mesh:
            raise ValueError(
                f"mol_type must be on the same device mesh as token_pad_mask, got {mol_type.device_mesh} and {device_mesh}"
            )
        if logits.device_mesh != device_mesh:
            raise ValueError(
                f"logits must be on the same device mesh as token_pad_mask, got {logits.device_mesh} and {device_mesh}"
            )

        group_replicate = device_mesh.get_group("cp_axis_1")

        maski_local = mask_collinear_pred.to_local().bool()
        maski_local = maski_local.reshape(-1, maski_local.shape[-1])
        mask_pad_local = token_pad_mask.to_local().bool().repeat_interleave(multiplicity, dim=0)
        asym_id_local = asym_id_base.to_local().repeat_interleave(multiplicity, dim=0)

        pair_mask_row_local = mask_pad_local & maski_local
        pair_mask_ptm_local = distributed_outer_op(
            pair_mask_row_local,
            op=OuterOp.BITAND,
            axis=1,
            input_t=mask_pad_local,
            transpose_comm=transpose_comm,
            group_replicate=group_replicate,
        )
        pair_mask_iptm_equal_local = distributed_outer_op(
            asym_id_local,
            op=OuterOp.EQUAL,
            axis=1,
            transpose_comm=transpose_comm,
            group_replicate=group_replicate,
        )
        pair_mask_iptm_local = pair_mask_ptm_local & (~pair_mask_iptm_equal_local)

        token_type_local = mol_type.to_local().repeat_interleave(multiplicity, dim=0)
        is_ligand_token_local = token_type_local == const.chain_type_ids["NONPOLYMER"]
        is_protein_token_local = token_type_local == const.chain_type_ids["PROTEIN"]
        ligand_iptm_mask_row_local = distributed_outer_op(
            is_ligand_token_local.bool(),
            op=OuterOp.BITAND,
            axis=1,
            transpose_comm=transpose_comm,
            group_replicate=group_replicate,
            input_t=is_protein_token_local,
        )
        ligand_iptm_mask_col_local = distributed_outer_op(
            is_protein_token_local.bool(),
            op=OuterOp.BITAND,
            axis=1,
            transpose_comm=transpose_comm,
            group_replicate=group_replicate,
            input_t=is_ligand_token_local,
        )
        ligand_iptm_mask_local = (ligand_iptm_mask_row_local | ligand_iptm_mask_col_local) & pair_mask_iptm_local

        protein_iptm_mask_local = distributed_outer_op(
            is_protein_token_local.bool(),
            op=OuterOp.BITAND,
            axis=1,
            transpose_comm=transpose_comm,
            group_replicate=group_replicate,
        )
        protein_iptm_mask_local = protein_iptm_mask_local & pair_mask_iptm_local

        reduced_dim = token_pad_mask.ndim - 1
        n_res_local = _LocalShardedSum.apply(
            mask_pad_local,
            reduced_dim,
            token_pad_mask.placements,
            device_mesh,
        ).unsqueeze(reduced_dim)
        logits_local = logits.to_local().detach()
        n_res_local = n_res_local.detach()

        num_bins = logits_local.shape[-1]
        bin_width = 32.0 / num_bins
        # Use at least float32 for bin centers to match serial code's default-dtype torch.arange
        compute_dtype = torch.promote_types(logits_local.dtype, torch.float32)
        pae_value = (torch.arange(num_bins, device=logits_local.device, dtype=compute_dtype) + 0.5) * bin_width
        pae_value = pae_value.unsqueeze(0)
        tm_value = serial_tm_function(pae_value, n_res_local).unsqueeze(1).unsqueeze(2)
        probs = F.softmax(logits_local.to(compute_dtype), dim=-1)
        tm_expected_value_local = torch.sum(probs * tm_value, dim=-1)

        reduced_dim = tm_expected_value_local.ndim - 1
        ptm_shape1, ptm_stride1 = _reduced_shape_stride(
            logits.shape[:-1], LayoutRightMap(tuple(logits.shape[:-1])).strides, reduced_dim
        )
        ptm_placements1 = _reduced_placements(logits.placements, logits.shape[:-1], reduced_dim)
        ptm_shape2, ptm_stride2 = _reduced_shape_stride(ptm_shape1, ptm_stride1, 1)
        output_placements = _reduced_placements(ptm_placements1, ptm_shape1, 1)
        mask_placements = (Shard(0), Shard(1), Shard(2))

        ptm_mask_local = pair_mask_ptm_local.bool()
        ptm_numerator_local = _LocalShardedSum.apply(
            tm_expected_value_local.masked_fill(~ptm_mask_local, 0),
            reduced_dim,
            logits.placements,
            device_mesh,
        )
        ptm_denominator_local = _LocalShardedSum.apply(
            ptm_mask_local,
            reduced_dim,
            mask_placements,
            device_mesh,
        )
        ptm_local = ptm_numerator_local / (ptm_denominator_local.to(tm_expected_value_local.dtype) + _EPS)
        ptm_local = _LocalShardedMax.apply(ptm_local, 1, ptm_placements1, device_mesh)

        iptm_mask_local = pair_mask_iptm_local.bool()
        iptm_numerator_local = _LocalShardedSum.apply(
            tm_expected_value_local.masked_fill(~iptm_mask_local, 0),
            reduced_dim,
            logits.placements,
            device_mesh,
        )
        iptm_denominator_local = _LocalShardedSum.apply(
            iptm_mask_local,
            reduced_dim,
            mask_placements,
            device_mesh,
        )
        iptm_local = iptm_numerator_local / (iptm_denominator_local.to(tm_expected_value_local.dtype) + _EPS)
        iptm_local = _LocalShardedMax.apply(iptm_local, 1, ptm_placements1, device_mesh)

        ligand_mask_local = ligand_iptm_mask_local.bool()
        ligand_num_local = _LocalShardedSum.apply(
            tm_expected_value_local.masked_fill(~ligand_mask_local, 0),
            reduced_dim,
            logits.placements,
            device_mesh,
        )
        ligand_den_local = _LocalShardedSum.apply(
            ligand_mask_local,
            reduced_dim,
            mask_placements,
            device_mesh,
        )
        ligand_local = ligand_num_local / (ligand_den_local.to(tm_expected_value_local.dtype) + _EPS)
        ligand_local = _LocalShardedMax.apply(ligand_local, 1, ptm_placements1, device_mesh)

        protein_mask_local = protein_iptm_mask_local.bool()
        protein_num_local = _LocalShardedSum.apply(
            tm_expected_value_local.masked_fill(~protein_mask_local, 0),
            reduced_dim,
            logits.placements,
            device_mesh,
        )
        protein_den_local = _LocalShardedSum.apply(
            protein_mask_local,
            reduced_dim,
            mask_placements,
            device_mesh,
        )
        protein_local = protein_num_local / (protein_den_local.to(tm_expected_value_local.dtype) + _EPS)
        protein_local = _LocalShardedMax.apply(protein_local, 1, ptm_placements1, device_mesh)

        chain_pair_iptm: dict[int, dict[int, DTensor]] = {}
        local_asym_ids = set(torch.unique(asym_id_local).tolist())
        cp_axis_0_group = device_mesh.get_group("cp_axis_0")
        cp_obj_list = [None] * torch.distributed.get_world_size(group=cp_axis_0_group)
        torch.distributed.all_gather_object(cp_obj_list, local_asym_ids, group=cp_axis_0_group)
        cp_asym_ids = set().union(*cp_obj_list)

        dp_group = device_mesh.get_group("dp")
        dp_obj_list = [None] * torch.distributed.get_world_size(group=dp_group)
        torch.distributed.all_gather_object(dp_obj_list, cp_asym_ids, group=dp_group)
        world_asym_ids_list = sorted(set().union(*dp_obj_list))
        for idx1 in world_asym_ids_list:
            chain_iptm: dict[int, DTensor] = {}
            for idx2 in world_asym_ids_list:
                if idx1 not in cp_asym_ids or idx2 not in cp_asym_ids:
                    iptm_chain_local = torch.full(
                        (mask_pad_local.size(0),),
                        CHAIN_IPTM_SENTINEL,
                        device=mask_pad_local.device,
                        dtype=tm_expected_value_local.dtype,
                    )
                else:
                    mask_pair_chain_row = maski_local & (asym_id_local == idx2) & mask_pad_local
                    mask_pair_chain_col = (asym_id_local == idx1) & mask_pad_local
                    mask_pair_chain_local = distributed_outer_op(
                        mask_pair_chain_row,
                        op=OuterOp.BITAND,
                        axis=1,
                        transpose_comm=transpose_comm,
                        group_replicate=group_replicate,
                        input_t=mask_pair_chain_col,
                    )
                    mask_pair_chain_local = mask_pair_chain_local.bool()
                    numerator_local = _LocalShardedSum.apply(
                        tm_expected_value_local.masked_fill(~mask_pair_chain_local, 0),
                        reduced_dim,
                        logits.placements,
                        device_mesh,
                    )
                    denominator_local = _LocalShardedSum.apply(
                        mask_pair_chain_local,
                        reduced_dim,
                        mask_placements,
                        device_mesh,
                    )
                    iptm_chain_local = numerator_local / (denominator_local.to(tm_expected_value_local.dtype) + _EPS)
                    iptm_chain_local = _LocalShardedMax.apply(iptm_chain_local, 1, ptm_placements1, device_mesh)

                chain_iptm[idx2] = DTensor.from_local(
                    iptm_chain_local,
                    device_mesh=device_mesh,
                    placements=output_placements,
                    shape=ptm_shape2,
                    stride=ptm_stride2,
                )
            chain_pair_iptm[idx1] = chain_iptm

        ptm = DTensor.from_local(
            ptm_local,
            device_mesh=device_mesh,
            placements=output_placements,
            shape=ptm_shape2,
            stride=ptm_stride2,
        )
        iptm = DTensor.from_local(
            iptm_local,
            device_mesh=device_mesh,
            placements=output_placements,
            shape=ptm_shape2,
            stride=ptm_stride2,
        )
        ligand_iptm = DTensor.from_local(
            ligand_local,
            device_mesh=device_mesh,
            placements=output_placements,
            shape=ptm_shape2,
            stride=ptm_stride2,
        )
        protein_iptm = DTensor.from_local(
            protein_local,
            device_mesh=device_mesh,
            placements=output_placements,
            shape=ptm_shape2,
            stride=ptm_stride2,
        )

        return ptm, iptm, ligand_iptm, protein_iptm, chain_pair_iptm

    @staticmethod
    @torch.amp.custom_bwd(device_type="cuda")
    def backward(  # type: ignore[override]
        ctx: FunctionCtx,
        grad_ptm: DTensor,
        grad_iptm: DTensor,
        grad_ligand_iptm: DTensor,
        grad_protein_iptm: DTensor,
        grad_chain_pair: dict[int, dict[int, DTensor]],
    ) -> tuple[None, None, None, None, None, None, None]:
        """Backward is not supported; returns None for all inputs."""
        return None, None, None, None, None, None, None


def compute_ptms(
    logits: DTensor,
    x_preds: DTensor,
    feats: dict[str, DTensor],
    multiplicity: int,
    transpose_comm: TransposeComm,
) -> tuple[DTensor, DTensor, DTensor, DTensor, dict[int, dict[int, DTensor]]]:
    """Compute pTM and ipTM scores for DTensor inputs.

    This redistributes PAE logits and token-level features to replicated placements
    to compute global reductions (max/sum) across tokens and chains.

    Args:
        logits: DTensor of shape (batch * multiplicity, num_tokens, num_tokens, num_bins)
            with placements (Shard(0), Shard(1), Shard(2)) containing PAE prediction logits.
        x_preds: DTensor of shape (batch * multiplicity, num_atoms, 3) with placements
            (Shard(0), Shard(1), Replicate()) containing predicted atom coordinates.
        feats: DTensor feature dict with required keys (frames_idx, asym_id, atom_to_token,
            atom_pad_mask, atom_resolved_mask, mol_type, token_pad_mask).
        multiplicity: Number of copies per sample in the batch dimension.
        transpose_comm: TransposeComm object for distributed outer operations.

    Returns:
        Tuple containing:
        - ptm: DTensor of shape (batch,) with confidence scores for predicted templates
        - iptm: DTensor of shape (batch,) with interface confidence scores
        - ligand_iptm: DTensor of shape (batch,) with ligand-protein interface scores
        - protein_iptm: DTensor of shape (batch,) with protein-protein interface scores
        - chain_pair_iptm: Dict mapping chain pairs to their interface confidence
          DTensors. Keys are the world-level union of chain IDs (homogeneous across
          all DP ranks). Entries where a chain pair does not exist on this DP rank's
          batch are filled with sentinel value CHAIN_IPTM_SENTINEL (-1.0).
    """
    feats_keys = {
        "frames_idx",
        "asym_id",
        "atom_to_token",
        "atom_pad_mask",
        "atom_resolved_mask",
        "mol_type",
        "token_pad_mask",
    }
    if any(k not in feats for k in feats_keys):
        raise ValueError(f"feats must contain the following keys: {feats_keys}, got {feats.keys()}")
    if not isinstance(logits, DTensor):
        raise TypeError(f"Expected DTensor for logits, got {type(logits)}")
    if not isinstance(x_preds, DTensor):
        raise TypeError(f"Expected DTensor for x_preds, got {type(x_preds)}")
    if feats["frames_idx"].ndim == 4:
        raise ValueError(
            f"frames_idx has unsqueezed ensemble dim (ndim=4, shape={feats['frames_idx'].shape}). "
            "Only E=1 is supported; squeeze the ensemble dim before calling compute_ptms."
        )

    device_mesh = logits.device_mesh
    if x_preds.device_mesh != device_mesh:
        raise ValueError(
            f"x_preds must be on the same device mesh as logits, got {x_preds.device_mesh} and {device_mesh}"
        )
    for key in feats_keys:
        if feats[key].device_mesh != device_mesh:
            raise ValueError(
                f"feats[{key}] must be on the same device mesh as logits, got {feats[key].device_mesh} and {device_mesh}"
            )

    _, mask_collinear_pred = compute_frame_pred(
        x_preds,
        feats["frames_idx"],
        feats,
        multiplicity,
        inference=True,
    )

    ptm, iptm, ligand_iptm, protein_iptm, chain_pair_iptm = _ComputePtmsImpl.apply(
        mask_collinear_pred,
        feats["token_pad_mask"],
        feats["asym_id"],
        feats["mol_type"],
        logits,
        multiplicity,
        transpose_comm,
    )

    return ptm, iptm, ligand_iptm, protein_iptm, chain_pair_iptm
