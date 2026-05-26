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

from copy import deepcopy
from typing import Optional

import torch
import torch.distributed as dist
import torch.nn.functional as F
from torch import Tensor
from torch.autograd.function import FunctionCtx
from torch.distributed.tensor import DTensor, Partial, Replicate, Shard, distribute_tensor

from boltz.data import const
from boltz.distributed.comm import One2OneComm, TransposeComm
from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.layers.atom_to_token import single_repr_token_to_atom
from boltz.distributed.model.layers.clip import clip
from boltz.distributed.model.layers.elementwise_op import (
    ElementwiseOp,
    elementwise_op,
    scalar_tensor_op,
)
from boltz.distributed.model.layers.redistribute_transpose import redistribute_transpose
from boltz.distributed.model.layers.repeat_interleave import shardwise_repeat_interleave
from boltz.distributed.model.layers.sharded_op import sharded_sum
from boltz.distributed.model.loss.triton.cdist_lddt import cdist_lddt
from boltz.distributed.model.loss.triton.cdist_pde import cdist_pde
from boltz.distributed.utils import LayoutMap, LayoutRightMap, get_group_rank_from_axial_shift
from boltz.model.layers.confidence_utils import compute_collinear_mask


class _ResolvedNegativeLogLikelihoodImpl(torch.autograd.Function):
    """Shardwise computation of resolved negative log-likelihood.

    This implements the forward and backward passes for computing the binary
    cross-entropy loss for predicting whether atoms are resolved. The computation
    is shardwise (no communication required) due to the block-diagonal structure
    of token_to_rep_atom with intersperse padding.

    The NLL computation follows:
        ref_mask = bmm(token_to_rep_atom, resolved_mask)
        log_probs = log_softmax(pred_resolved, dim=-1)
        errors = -ref_mask * log_probs[:,:,0] - (1 - ref_mask) * log_probs[:,:,1]

    The backward pass uses PyTorch autograd on the local computation graph,
    avoiding the need for manual gradient derivation.

    See Also
    --------
    resolved_negative_log_likelihood : The public API function that calls this.
    """

    @staticmethod
    @torch.amp.custom_fwd(device_type="cuda")
    def forward(
        ctx: FunctionCtx,
        pred_resolved: DTensor,
        token_to_rep_atom: DTensor,
        true_coords_resolved_mask: DTensor,
    ) -> DTensor:
        """Forward pass for computing resolved negative log-likelihood.

        Parameters
        ----------
        ctx : FunctionCtx
            The autograd context object for saving tensors for backward.
        pred_resolved : DTensor
            Predicted resolved logits with shape (B*mult, N_token, 2).
            Placements: (Shard(0), Shard(1), Replicate())
        token_to_rep_atom : DTensor
            One-hot mapping from tokens to representative atoms (non-multiplexed).
            Shape (B, N_token, N_atom) with block-diagonal structure.
            Placements: (Shard(0), Shard(1), Replicate())
        true_coords_resolved_mask : DTensor
            Resolved mask for atoms. Shape (B*mult, N_atom).
            Placements: (Shard(0), Shard(1), Replicate())

        Returns
        -------
        DTensor
            Error tensor with shape (B*mult, N_token).
            Placements: (Shard(0), Shard(1), Replicate())

        Raises
        ------
        TypeError
            If inputs are not DTensors.
        ValueError
            If Partial placements are present or placements don't match expected.
        """
        # Type checking
        if not isinstance(pred_resolved, DTensor):
            raise TypeError(f"Expected DTensor for pred_resolved, got {type(pred_resolved)}")
        if not isinstance(token_to_rep_atom, DTensor):
            raise TypeError(f"Expected DTensor for token_to_rep_atom, got {type(token_to_rep_atom)}")
        if not isinstance(true_coords_resolved_mask, DTensor):
            raise TypeError(f"Expected DTensor for true_coords_resolved_mask, got {type(true_coords_resolved_mask)}")

        device_mesh = pred_resolved.device_mesh
        pred_placements = pred_resolved.placements

        # Validate placements - check no Partial and expected structure
        for i_dim_device_mesh, placement in enumerate(pred_placements):
            if isinstance(placement, Partial):
                raise ValueError("Partial placements are not supported for pred_resolved")
            elif isinstance(placement, Shard):
                # Check that sharded dimensions are evenly divided
                if pred_resolved.shape[placement.dim] % device_mesh.shape[i_dim_device_mesh] != 0:
                    raise ValueError(
                        f"Uneven sharding of tensor dimension {placement.dim} of size "
                        f"{pred_resolved.shape[placement.dim]} along device mesh dimension "
                        f"{i_dim_device_mesh} of size {device_mesh.shape[i_dim_device_mesh]} is not supported"
                    )

        # Validate device_mesh consistency across all inputs
        if token_to_rep_atom.device_mesh != device_mesh:
            raise ValueError(
                f"Device mesh mismatch: pred_resolved has {device_mesh}, "
                f"token_to_rep_atom has {token_to_rep_atom.device_mesh}"
            )
        if true_coords_resolved_mask.device_mesh != device_mesh:
            raise ValueError(
                f"Device mesh mismatch: pred_resolved has {device_mesh}, "
                f"true_coords_resolved_mask has {true_coords_resolved_mask.device_mesh}"
            )

        # Validate placements consistency across all inputs
        # All inputs should have the same placements on the 3D mesh: (Shard(0), Shard(1), Replicate())
        if token_to_rep_atom.placements != pred_placements:
            raise ValueError(
                f"Placements mismatch: pred_resolved has {pred_placements}, "
                f"token_to_rep_atom has {token_to_rep_atom.placements}"
            )
        if true_coords_resolved_mask.placements != pred_placements:
            raise ValueError(
                f"Placements mismatch: pred_resolved has {pred_placements}, "
                f"true_coords_resolved_mask has {true_coords_resolved_mask.placements}"
            )

        # Validate shape consistency across all inputs
        # Extract dimensions from token_to_rep_atom: (B, N_token, N_atom_padded)
        if len(token_to_rep_atom.shape) != 3:
            raise ValueError(f"token_to_rep_atom must be 3D, got shape {token_to_rep_atom.shape}")
        batch_size = token_to_rep_atom.shape[0]
        n_token = token_to_rep_atom.shape[1]

        # Validate true_coords_resolved_mask shape: (B*mult, N_atom_padded)
        if len(true_coords_resolved_mask.shape) != 2 or true_coords_resolved_mask.shape[0] % batch_size != 0:
            raise ValueError(
                f"true_coords_resolved_mask must be 2D with shape[0] divisible by batch_size ({batch_size}), "
                f"got shape {tuple(true_coords_resolved_mask.shape)}"
            )
        multiplicity = true_coords_resolved_mask.shape[0] // batch_size

        # Validate pred_resolved shape: (B*mult, N_token, 2)
        expected_pred_shape = (batch_size * multiplicity, n_token, 2)
        if tuple(pred_resolved.shape) != expected_pred_shape:
            raise ValueError(
                f"Shape mismatch: pred_resolved has shape {tuple(pred_resolved.shape)}, expected {expected_pred_shape}"
            )

        # Detach and set requires_grad to build a local computation graph.
        # Use promote_types to match serial resolved_loss which casts via .float();
        # promote_types promotes to at least float32 while preserving float64.
        compute_dtype = torch.promote_types(pred_resolved.dtype, torch.float32)
        pred_local = (
            pred_resolved.to_local().detach().to(dtype=compute_dtype).requires_grad_(pred_resolved.requires_grad)
        )
        token_to_rep_atom_local = token_to_rep_atom.to_local().to(dtype=compute_dtype)
        resolved_mask_local = true_coords_resolved_mask.to_local().to(dtype=compute_dtype)

        # Validate n_atom consistency on local shards (both padded to max_atoms_per_shard)
        if token_to_rep_atom_local.shape[-1] != resolved_mask_local.shape[-1]:
            raise ValueError(
                f"Local shard atom dimension mismatch: token_to_rep_atom has {token_to_rep_atom_local.shape[-1]}, "
                f"true_coords_resolved_mask has {resolved_mask_local.shape[-1]}"
            )

        # Infer multiplicity from local shapes
        # token_to_rep_atom_local: (B_local, N_token_local, N_atom_padded)
        # resolved_mask_local: (B_local*mult, N_atom_padded)
        b_local = token_to_rep_atom_local.shape[0]
        multiplicity = resolved_mask_local.shape[0] // b_local

        with torch.enable_grad():
            # Build a local computation graph for the shardwise operations
            # Reshape resolved_mask to (B_local, mult, N_atom_padded) for einsum
            resolved_mask_reshaped = resolved_mask_local.view(b_local, multiplicity, -1)

            # Use einsum to compute ref_mask without repeat_interleave on token_to_rep_atom
            # token_to_rep_atom_local: (B_local, N_token_local, N_atom_padded) -> "btj"
            # resolved_mask_reshaped: (B_local, mult, N_atom_padded) -> "bmj"
            # ref_mask: (B_local, mult, N_token_local) -> "bmt"
            ref_mask = torch.einsum("btj,bmj->bmt", token_to_rep_atom_local, resolved_mask_reshaped)
            # Flatten to (B_local*mult, N_token_local)
            ref_mask = ref_mask.flatten(0, 1)

            # Compute log softmax probabilities
            log_softmax_resolved = F.log_softmax(pred_local, dim=-1)

            # Compute binary cross-entropy errors
            # errors = -ref_mask * log_probs[resolved] - (1 - ref_mask) * log_probs[unresolved]
            errors = -ref_mask * log_softmax_resolved[:, :, 0] - (1 - ref_mask) * log_softmax_resolved[:, :, 1]

        # Compute output shape and stride (same batch and token dims, remove bins dim)
        output_shape = tuple(pred_resolved.shape[:-1])
        output_stride = LayoutRightMap(output_shape).strides

        # Output placements: same as input but without the last (bins) dimension
        # Since pred_resolved has (Shard(0), Shard(1), Replicate()), output is (Shard(0), Shard(1), Replicate())
        # but we need to handle the 2D case
        output_placements = pred_placements

        # Save tensors for backward pass
        ctx.save_for_backward(pred_local, errors)
        ctx.device_mesh = device_mesh
        ctx.pred_placements = pred_placements
        ctx.pred_shape = pred_resolved.shape
        ctx.pred_stride = pred_resolved.stride()
        ctx.output_shape = output_shape
        ctx.output_stride = output_stride
        ctx.output_placements = output_placements

        # Create output DTensor
        result = DTensor.from_local(
            errors.detach(),
            device_mesh=device_mesh,
            placements=output_placements,
            shape=output_shape,
            stride=output_stride,
        )

        return result

    @staticmethod
    @torch.amp.custom_bwd(device_type="cuda")
    def backward(
        ctx: FunctionCtx,
        grad_output: DTensor,
    ) -> tuple[DTensor | None, None, None]:
        """Backward pass for computing resolved negative log-likelihood.

        Computes gradients by backpropagating through the local computation graph
        that was built during the forward pass. This leverages PyTorch's autograd
        rather than manual gradient computation.

        Parameters
        ----------
        ctx : FunctionCtx
            The autograd context containing saved tensors from forward.
        grad_output : DTensor
            Gradient of loss with respect to output errors.

        Returns
        -------
        tuple[DTensor | None, None, None]
            Gradients for each forward input in order:
            - d_pred_resolved: DTensor or None, gradient for pred_resolved
            - None: token_to_rep_atom (non-differentiable, one-hot)
            - None: true_coords_resolved_mask (non-differentiable, ground truth)
        """
        pred_local, errors_local = ctx.saved_tensors

        if not pred_local.requires_grad:
            return None, None, None

        grad_output_local = grad_output.to_local()

        # Backprop via the local graph
        (d_pred_local,) = torch.autograd.grad(
            outputs=[errors_local],
            inputs=[pred_local],
            grad_outputs=[grad_output_local],
            retain_graph=False,  # Frees the local graph immediately
        )

        # Wrap gradient in DTensor
        d_pred = DTensor.from_local(
            d_pred_local,
            device_mesh=ctx.device_mesh,
            placements=ctx.pred_placements,
            shape=ctx.pred_shape,
            stride=ctx.pred_stride,
        )

        return d_pred, None, None


def resolved_negative_log_likelihood(
    pred_resolved: DTensor,
    token_to_rep_atom: DTensor,
    true_coords_resolved_mask: DTensor,
) -> DTensor:
    """Compute shardwise negative log-likelihood for resolved prediction.

    This is the DTensor-compatible version of the resolved loss NLL computation.
    All operations are shardwise (no inter-rank communication required) due to
    the block-diagonal structure of token_to_rep_atom with intersperse padding.

    The multiplicity is inferred from the shapes:
    multiplicity = true_coords_resolved_mask.shape[0] // token_to_rep_atom.shape[0]

    Parameters
    ----------
    pred_resolved : DTensor
        Predicted resolved logits with shape (B*mult, N_token, 2).
        Placements: (Shard(0), Shard(1), Replicate())
    token_to_rep_atom : DTensor
        One-hot mapping from tokens to representative atoms (non-multiplexed).
        Shape (B, N_token, N_atom) with block-diagonal structure.
        Placements: (Shard(0), Shard(1), Replicate())
    true_coords_resolved_mask : DTensor
        Resolved mask for atoms. Shape (B*mult, N_atom).
        Placements: (Shard(0), Shard(1), Replicate())

    Returns
    -------
    DTensor
        Error tensor with shape (B*mult, N_token).
        Placements: (Shard(0), Shard(1), Replicate())

    See Also
    --------
    resolved_loss : Full loss function that uses this for NLL computation.
    boltz.model.loss.confidencev2.resolved_loss : Serial version.
    """
    return _ResolvedNegativeLogLikelihoodImpl.apply(pred_resolved, token_to_rep_atom, true_coords_resolved_mask)


def resolved_loss(
    pred_resolved: DTensor,
    feats: dict[str, DTensor],
    true_coords_resolved_mask: DTensor,
    multiplicity: int = 1,
) -> DTensor:
    """Compute resolved loss using DTensor operations.

    This is the DTensor-compatible version of the resolved_loss function.
    It computes binary cross-entropy loss for predicting whether atoms are resolved
    in the structure.

    The computation is split into:
    1. Part a) Shardwise NLL computation via resolved_negative_log_likelihood
    2. Part b) Distributed weighted sum along token axis, then mean over batch

    Parameters
    ----------
    pred_resolved : DTensor
        Predicted resolved logits with shape (B*mult, N_token, 2).
        Placements: (Shard(0), Shard(1), Replicate())
    feats : dict[str, DTensor]
        Feature dictionary containing:
        - token_to_rep_atom: One-hot mapping (B, N_token, N_atom)
        - token_pad_mask: Padding mask (B, N_token)
    true_coords_resolved_mask : DTensor
        Resolved mask for atoms. Shape (B*mult, N_atom).
        Placements: (Shard(0), Shard(1), Replicate())
    multiplicity : int, optional
        Diffusion batch multiplier, by default 1

    Returns
    -------
    DTensor
        Scalar loss with placements (Replicate(), Replicate(), Replicate())

    See Also
    --------
    resolved_negative_log_likelihood : Shardwise NLL computation.
    boltz.model.loss.confidencev2.resolved_loss : Serial version.
    """
    # Part a) - Shardwise NLL computation (token_to_rep_atom is non-multiplexed)
    errors = resolved_negative_log_likelihood(pred_resolved, feats["token_to_rep_atom"], true_coords_resolved_mask)

    # Part b) - Weighted sum along token axis (dim=-1), then mean over batch
    # Expand pad_mask with multiplicity
    pad_mask = shardwise_repeat_interleave(feats["token_pad_mask"], multiplicity, dim=0)
    # Following diffusion.py pattern (lines 413-420)
    # num = sum(errors * pad_mask, dim=-1)
    num = sharded_sum(elementwise_op(errors, pad_mask, ElementwiseOp.PROD), dim=-1)
    # den = sum(pad_mask, dim=-1)
    den = sharded_sum(pad_mask, dim=-1)
    # loss_per_sample = num / max(den, 1e-7)
    loss_per_sample = elementwise_op(num, clip(den, min_val=1e-7, max_val=None), ElementwiseOp.DIV)

    # Mean over batch dimension (following diffusion.py lines 423-427)
    loss = scalar_tensor_op(
        1.0 / loss_per_sample.shape[0],
        sharded_sum(loss_per_sample, dim=0),
        ElementwiseOp.PROD,
    )

    return loss


# PAE loss numerical stability constants
FRAME_NORM_EPS = 1e-5  # prevents division by zero in frame basis normalization
PAE_DIST_EPS = 1e-8  # prevents sqrt(0) in PAE target distance computation
PAE_LOSS_DENOM_EPS = 1e-7  # prevents division by zero when normalizing by mask sum


def _check_pae_input_consistency(
    pred_pae: DTensor,
    pred_atom_coords: DTensor,
    true_atom_coords: DTensor,
    true_coords_resolved_mask: DTensor,
    feats: dict[str, DTensor],
    multiplicity: int,
) -> None:
    """Validate input tensors for PAE loss computation.

    Checks type, device mesh consistency, placement correctness, and shape
    compatibility for all inputs to pae_loss.

    Parameters
    ----------
    pred_pae : DTensor
        Predicted PAE logits. Expected shape: (B*mult, N_token, N_token, bins).
        Expected placements: (Shard(0), Shard(1), Shard(2)).
    pred_atom_coords : DTensor
        Predicted atom coordinates. Expected shape: (B*mult, N_atom, 3).
        Expected placements: (Shard(0), Shard(1), Replicate()).
    true_atom_coords : DTensor
        True atom coordinates. Expected shape: (B*mult, N_atom, 3).
        Expected placements: (Shard(0), Shard(1), Replicate()).
    true_coords_resolved_mask : DTensor
        Resolved mask for atoms. Expected shape: (B*mult, N_atom).
        Expected placements: (Shard(0), Shard(1), Replicate()).
    feats : dict[str, DTensor]
        Feature dictionary containing frames_idx, frame_resolved_mask, token_pad_mask, etc.
    multiplicity : int
        Diffusion batch multiplier.

    Raises
    ------
    TypeError
        If any input is not a DTensor.
    ValueError
        If device meshes are inconsistent, placements are incorrect, shapes don't match,
        or sharding is uneven.
    """
    # --- Type checks ---
    if not isinstance(pred_pae, DTensor):
        raise TypeError(f"Expected DTensor for pred_pae, got {type(pred_pae)}")
    if not isinstance(pred_atom_coords, DTensor):
        raise TypeError(f"Expected DTensor for pred_atom_coords, got {type(pred_atom_coords)}")
    if not isinstance(true_atom_coords, DTensor):
        raise TypeError(f"Expected DTensor for true_atom_coords, got {type(true_atom_coords)}")
    if not isinstance(true_coords_resolved_mask, DTensor):
        raise TypeError(f"Expected DTensor for true_coords_resolved_mask, got {type(true_coords_resolved_mask)}")

    # --- Device mesh consistency ---
    device_mesh = pred_pae.device_mesh
    if pred_atom_coords.device_mesh != device_mesh:
        raise ValueError(
            f"Device mesh mismatch: pred_pae has {device_mesh}, pred_atom_coords has {pred_atom_coords.device_mesh}"
        )
    if true_atom_coords.device_mesh != device_mesh:
        raise ValueError(
            f"Device mesh mismatch: pred_pae has {device_mesh}, true_atom_coords has {true_atom_coords.device_mesh}"
        )
    if true_coords_resolved_mask.device_mesh != device_mesh:
        raise ValueError(
            f"Device mesh mismatch: pred_pae has {device_mesh}, "
            f"true_coords_resolved_mask has {true_coords_resolved_mask.device_mesh}"
        )

    # --- Placement validation for pred_pae: (Shard(0), Shard(1), Shard(2)) ---
    expected_pae_placements = (Shard(0), Shard(1), Shard(2))
    if pred_pae.placements != expected_pae_placements:
        raise ValueError(f"pred_pae must have placements {expected_pae_placements}, got {pred_pae.placements}")

    # Check sharding divisibility for pred_pae
    # Shard(0) -> dim 0 (B*mult) sharded by mesh dim 0 (dp)
    # Shard(1) -> dim 1 (N_token) sharded by mesh dim 1 (cp_axis_0)
    # Shard(2) -> dim 2 (N_token) sharded by mesh dim 2 (cp_axis_1)
    for mesh_dim, placement in enumerate(pred_pae.placements):
        if isinstance(placement, Shard):
            tensor_dim_size = pred_pae.shape[placement.dim]
            mesh_dim_size = device_mesh.shape[mesh_dim]
            if tensor_dim_size % mesh_dim_size != 0:
                raise ValueError(
                    f"pred_pae dimension {placement.dim} (size {tensor_dim_size}) "
                    f"is not evenly divisible by mesh dimension {mesh_dim} (size {mesh_dim_size})"
                )

    # --- Placement validation for coords: (Shard(0), Shard(1), Replicate()) ---
    expected_coords_placements = (Shard(0), Shard(1), Replicate())
    if pred_atom_coords.placements != expected_coords_placements:
        raise ValueError(
            f"pred_atom_coords must have placements {expected_coords_placements}, got {pred_atom_coords.placements}"
        )
    if true_atom_coords.placements != expected_coords_placements:
        raise ValueError(
            f"true_atom_coords must have placements {expected_coords_placements}, got {true_atom_coords.placements}"
        )

    # Check sharding divisibility for coords
    for mesh_dim, placement in enumerate(pred_atom_coords.placements):
        if isinstance(placement, Shard):
            tensor_dim_size = pred_atom_coords.shape[placement.dim]
            mesh_dim_size = device_mesh.shape[mesh_dim]
            if tensor_dim_size % mesh_dim_size != 0:
                raise ValueError(
                    f"pred_atom_coords dimension {placement.dim} (size {tensor_dim_size}) "
                    f"is not evenly divisible by mesh dimension {mesh_dim} (size {mesh_dim_size})"
                )

    # --- Shape validation ---
    # pred_pae: (B*mult, N_token, N_token, bins)
    if pred_pae.ndim != 4:
        raise ValueError(f"pred_pae must be 4D, got {pred_pae.ndim}D with shape {pred_pae.shape}")

    batch_mult_size = pred_pae.shape[0]

    if pred_pae.shape[1] != pred_pae.shape[2]:
        raise ValueError(f"pred_pae must have equal N_token dimensions (dims 1 and 2), got shape {pred_pae.shape}")

    if batch_mult_size % multiplicity != 0:
        raise ValueError(
            f"pred_pae batch dimension (shape[0]={batch_mult_size}) must be divisible by multiplicity ({multiplicity})"
        )

    # pred_atom_coords / true_atom_coords: (B*mult, N_atom, 3)
    if pred_atom_coords.ndim != 3:
        raise ValueError(
            f"pred_atom_coords must be 3D, got {pred_atom_coords.ndim}D with shape {pred_atom_coords.shape}"
        )
    if pred_atom_coords.shape[0] != batch_mult_size:
        raise ValueError(
            f"pred_atom_coords batch dimension (shape[0]={pred_atom_coords.shape[0]}) "
            f"does not match pred_pae batch dimension ({batch_mult_size})"
        )
    if pred_atom_coords.shape[2] != 3:
        raise ValueError(f"pred_atom_coords must have 3 coordinates, got shape {pred_atom_coords.shape}")

    if true_atom_coords.shape != pred_atom_coords.shape:
        raise ValueError(
            f"true_atom_coords shape {true_atom_coords.shape} "
            f"does not match pred_atom_coords shape {pred_atom_coords.shape}"
        )

    # true_coords_resolved_mask: (B*mult, N_atom)
    if true_coords_resolved_mask.ndim != 2:
        raise ValueError(
            f"true_coords_resolved_mask must be 2D, got {true_coords_resolved_mask.ndim}D "
            f"with shape {true_coords_resolved_mask.shape}"
        )
    if true_coords_resolved_mask.shape[0] != batch_mult_size:
        raise ValueError(
            f"true_coords_resolved_mask batch dimension (shape[0]={true_coords_resolved_mask.shape[0]}) "
            f"does not match pred_pae batch dimension ({batch_mult_size})"
        )

    n_atom = pred_atom_coords.shape[1]
    if true_coords_resolved_mask.shape[1] != n_atom:
        raise ValueError(
            f"true_coords_resolved_mask N_atom dimension ({true_coords_resolved_mask.shape[1]}) "
            f"does not match pred_atom_coords ({n_atom})"
        )

    if "frames_idx" in feats and feats["frames_idx"].ndim == 4:
        raise ValueError(
            f"frames_idx has unsqueezed ensemble dim (ndim=4, shape={feats['frames_idx'].shape}). "
            "Only E=1 is supported; squeeze the ensemble dim before calling pae_loss."
        )
    if "frame_resolved_mask" in feats and feats["frame_resolved_mask"].ndim == 3:
        raise ValueError(
            f"frame_resolved_mask has unsqueezed ensemble dim (ndim=3, shape={feats['frame_resolved_mask'].shape}). "
            "Only E=1 is supported; squeeze the ensemble dim before calling pae_loss."
        )


def pae_loss(
    pred_pae: DTensor,
    pred_atom_coords: DTensor,
    true_atom_coords: DTensor,
    true_coords_resolved_mask: DTensor,
    feats: dict[str, DTensor],
    comm: One2OneComm,
    dist_manager: DistributedManager,
    group_layout: LayoutMap,
    multiplicity: int = 1,
    max_dist: float = 32.0,
) -> DTensor:
    """Compute PAE (Predicted Aligned Error) loss using DTensor operations.

    Tensor-compatible version of the pae_loss function.
    It computes cross-entropy loss for predicting the alignment error between
    predicted and true atom coordinates when expressed in local reference frames.

    Sharding Strategy
    -----------------
    The 3D device mesh has shape (dp, cp_axis_0, cp_axis_1). Placements specify
    which tensor dimension to shard across each mesh dimension:

    ::

        pred_pae shape: (B*mult, N_token, N_token, bins)
        Placements:     (Shard(0), Shard(1), Shard(2))
                           │         │         │
                           ▼         ▼         ▼
        Mesh dims:        dp    cp_axis_0  cp_axis_1

        Example with B=2, mult=2, N=32, bins=64 on dp=2, cp=(2,2):

        Global pred_pae: (4, 32, 32, 64)

        Device Mesh (cp_axis_0 × cp_axis_1):
                    cp_axis_1
                  ┌──────┬──────┐
        cp_axis_0 │  R0  │  R1  │  R0: (2, 16, 16, 64) tokens[0:16, 0:16]
                  ├──────┼──────┤  R1: (2, 16, 16, 64) tokens[0:16, 16:32]
                  │  R2  │  R3  │  R2: (2, 16, 16, 64) tokens[16:32, 0:16]
                  └──────┴──────┘  R3: (2, 16, 16, 64) tokens[16:32, 16:32]

    Parameters
    ----------
    pred_pae : DTensor
        Predicted PAE logits with shape (B*mult, N_token, N_token, num_bins).
        Placements: (Shard(0), Shard(1), Shard(2))
    pred_atom_coords : DTensor
        Predicted atom coordinates with shape (B*mult, N_atom, 3).
        Placements: (Shard(0), Shard(1), Replicate())
    true_atom_coords : DTensor
        True atom coordinates with shape (B*mult, N_atom, 3).
        Placements: (Shard(0), Shard(1), Replicate())
    true_coords_resolved_mask : DTensor
        Resolved mask for atoms. Shape (B*mult, N_atom).
        Placements: (Shard(0), Shard(1), Replicate())
    feats : dict[str, DTensor]
        Feature dictionary containing:

        - frames_idx: Frame atom indices (B, N_token, 3)
          Placements: (Shard(0), Shard(1), Replicate())
        - frame_resolved_mask: Frame validity mask (B, N_token)
          Placements: (Shard(0), Shard(1), Replicate())
        - asym_id: Asymmetric unit IDs (B, N_token)
          Placements: (Shard(0), Shard(1), Replicate())
        - atom_to_token: Atom to token mapping (B, N_atom, N_token)
          Placements: (Shard(0), Shard(1), Replicate()) with intersperse padding
        - atom_pad_mask: Atom padding mask (B, N_atom)
          Placements: (Shard(0), Shard(1), Replicate()) with intersperse padding
        - mol_type: Molecule type (B, N_token)
          Placements: (Shard(0), Shard(1), Replicate())
        - token_pad_mask: Token padding mask (B, N_token)
          Placements: (Shard(0), Shard(1), Replicate())
        - atom_resolved_mask: Atom resolved mask (B, N_atom)
          Placements: (Shard(0), Shard(1), Replicate()) with intersperse padding
        - is_nonpolymer_with_frame: Non-polymer frame indicator (B, N_token)
          Placements: (Shard(0), Shard(1), Replicate())

    comm : One2OneComm
        Communication object for coordinate transpose operations in frame
        computation for non-polymers.
    dist_manager : DistributedManager
        Distributed manager for process group information.
    group_layout : LayoutMap
        Layout map for the 2D CP grid.
    multiplicity : int, optional
        Diffusion batch multiplier, by default 1
    max_dist : float, optional
        Maximum distance for PAE binning, by default 32.0

    Returns
    -------
    DTensor
        Scalar loss with placements (Replicate(), Replicate(), Replicate())

    See Also
    --------
    compute_frame_pred : Distributed frame computation for non-polymers.
    boltz.model.loss.confidencev2.pae_loss : Serial version.
    """
    _check_pae_input_consistency(
        pred_pae,
        pred_atom_coords,
        true_atom_coords,
        true_coords_resolved_mask,
        feats,
        multiplicity,
    )

    return _PAELossImpl.apply(
        pred_pae,
        pred_atom_coords,
        true_atom_coords,
        true_coords_resolved_mask,
        feats,
        comm,
        dist_manager,
        group_layout,
        multiplicity,
        max_dist,
    )


def _express_coordinate_in_frame_distributed(
    atom_coords: Tensor,
    frame_atom_a: Tensor,
    frame_atom_b: Tensor,
    frame_atom_c: Tensor,
    dist_manager: DistributedManager,
    group_layout: LayoutMap,
    transpose_comm: TransposeComm,
) -> Tensor:
    """Distributed express_coordinate_in_frame for 2D-sharded output.

    For rank (shard_i, shard_j) in the 2D CP mesh, computes the block
    [i_start:i_end, j_start:j_end] of the full N_token × N_token output.

    Args:
        atom_coords: [B, mult, N_atom_local, 3] local atom coordinates
        frame_atom_a/b/c: [B, mult, N_token_local] local frame indices (already
            converted to local atom indices)
        dist_manager: DistributedManager for communication
        group_layout: LayoutMap for the CP 2D mesh

    Returns:
        x_transformed: [B, mult, N_token_local_i, N_token_local_j, 3] transformed coordinates
        mask_collinear: [B, mult, N_token_local_i] collinear mask for row tokens
    """
    n_atoms_local = atom_coords.shape[2]
    n_tokens_local = frame_atom_a.shape[-1]

    cp_group = dist_manager.group["cp"]
    rank_coords = group_layout.unravel(dist.get_rank(cp_group))
    atom_offset = rank_coords[0] * n_atoms_local
    global_min = atom_offset
    global_max = atom_offset + n_atoms_local

    # Identify if any rank will rely on nonlocal frame coordinates. If so, all ranks will use the global gather path
    # to avoid deadlock in cooperative operations.
    frame_requires_global_coords = torch.any(
        (frame_atom_a < global_min)
        | (frame_atom_a >= global_max)
        | (frame_atom_b < global_min)
        | (frame_atom_b >= global_max)
        | (frame_atom_c < global_min)
        | (frame_atom_c >= global_max)
    ).to(torch.int32)
    dist.all_reduce(frame_requires_global_coords, op=dist.ReduceOp.MAX, group=cp_group)

    if frame_requires_global_coords.item() == 1:
        a = _gather_frame_coords(
            atom_coords,
            frame_atom_a,
            local_only=False,
            atom_offset=atom_offset,
            n_tokens_local=n_tokens_local,
            dist_manager=dist_manager,
            group_layout=group_layout,
        )
        b = _gather_frame_coords(
            atom_coords,
            frame_atom_b,
            local_only=False,
            atom_offset=atom_offset,
            n_tokens_local=n_tokens_local,
            dist_manager=dist_manager,
            group_layout=group_layout,
        )
        c = _gather_frame_coords(
            atom_coords,
            frame_atom_c,
            local_only=False,
            atom_offset=atom_offset,
            n_tokens_local=n_tokens_local,
            dist_manager=dist_manager,
            group_layout=group_layout,
        )
    else:
        a = _gather_frame_coords(
            atom_coords,
            frame_atom_a,
            local_only=True,
            atom_offset=atom_offset,
            n_tokens_local=n_tokens_local,
            dist_manager=dist_manager,
            group_layout=group_layout,
        )
        b = _gather_frame_coords(
            atom_coords,
            frame_atom_b,
            local_only=True,
            atom_offset=atom_offset,
            n_tokens_local=n_tokens_local,
            dist_manager=dist_manager,
            group_layout=group_layout,
        )
        c = _gather_frame_coords(
            atom_coords,
            frame_atom_c,
            local_only=True,
            atom_offset=atom_offset,
            n_tokens_local=n_tokens_local,
            dist_manager=dist_manager,
            group_layout=group_layout,
        )

    # Exchange b coordinates with transpose peer early to overlap with local frame basis work.
    b_j = transpose_comm.enqueue_to_dispatch(b.contiguous())

    # Build orthonormal frame from local a, b, c
    # a, b, c: [B, mult, N_token_local, 3]
    ab = a - b
    cb = c - b
    w1 = ab / (torch.norm(ab, dim=-1, keepdim=True) + FRAME_NORM_EPS)
    w2 = cb / (torch.norm(cb, dim=-1, keepdim=True) + FRAME_NORM_EPS)
    e1 = (w1 + w2) / (torch.norm(w1 + w2, dim=-1, keepdim=True) + FRAME_NORM_EPS)
    e2 = (w2 - w1) / (torch.norm(w2 - w1, dim=-1, keepdim=True) + FRAME_NORM_EPS)
    e3 = torch.linalg.cross(e1, e2)

    # Collinear mask from correctly gathered frame coordinates
    # Flatten (B, mult, N_token_local) to (...) for compute_collinear_mask
    orig_shape = ab.shape[:-1]  # (B, mult, N_token_local)
    mask_collinear = compute_collinear_mask(
        ab.reshape(-1, 3),
        cb.reshape(-1, 3),
    ).reshape(orig_shape)

    # Ensure transpose exchange completed before using b_j.
    transpose_comm.wait_until_finished()

    # Pairwise displacement: d[i,j] = b_j[j] - b[i]
    # b: [B, mult, N_token_local_i, 3] (local row tokens)
    # b_j: [B, mult, N_token_local_j, 3] (gathered column tokens)
    d = b_j[:, :, None, :, :] - b[:, :, :, None, :]  # [B, mult, N_i, N_j, 3]

    # Project onto local frame basis via batched matmul
    basis = torch.stack([e1, e2, e3], dim=-1)  # [B, mult, N_i, 3, 3]
    x_transformed = torch.matmul(d.unsqueeze(-2), basis[:, :, :, None, :, :]).squeeze(-2)

    return x_transformed, mask_collinear


class _PAELossImpl(torch.autograd.Function):
    """Distributed PAE loss computation with autograd support.

    Handles the forward pass with local autograd for gradient computation,
    and the backward pass with proper gradient routing for DTensors.

    The gradient w.r.t. pred_pae flows through log_softmax → gather → errors.
    The gradient w.r.t. pred_atom_coords is zero because bin_index computation
    (torch.floor) is non-differentiable.
    """

    @staticmethod
    @torch.amp.custom_fwd(device_type="cuda")
    def forward(
        ctx,
        pred_pae: DTensor,
        pred_atom_coords: DTensor,
        true_atom_coords: DTensor,
        true_coords_resolved_mask: DTensor,
        feats: dict[str, DTensor],
        comm: One2OneComm,
        dist_manager: DistributedManager,
        group_layout: LayoutMap,
        multiplicity: int,
        max_dist: float,
    ) -> DTensor:
        device_mesh = pred_pae.device_mesh
        num_bins = pred_pae.shape[-1]
        dp_size = device_mesh.shape[0]

        # Extract local tensors from DTensors
        # pred_pae is (B*mult, N_token, N_token, bins); unflatten to (B_batch, mult, N_i, N_j, bins)
        pred_pae_local = pred_pae.to_local().unflatten(0, (-1, multiplicity))
        pred_atom_coords_local = pred_atom_coords.to_local()
        true_atom_coords_local = true_atom_coords.to_local()
        true_coords_resolved_mask_local = true_coords_resolved_mask.to_local()

        frame_resolved_mask_local = feats["frame_resolved_mask"].to_local()
        token_pad_mask_local = feats["token_pad_mask"].to_local()

        # Get rank info for 2D sharding
        B_local, N_atom_local, _ = true_atom_coords_local.shape
        group_axis_0 = dist_manager.subgroups["cp"][0]
        group_rank_0 = dist.get_rank(group_axis_0)
        atom_offset = N_atom_local * group_rank_0
        if B_local % multiplicity != 0:
            raise ValueError(
                f"true_atom_coords local batch dim ({B_local}) must be " f"divisible by multiplicity ({multiplicity})"
            )
        B_batch = B_local // multiplicity
        if pred_atom_coords_local.shape[0] != true_atom_coords_local.shape[0]:
            raise ValueError(
                f"pred_atom_coords batch dim ({pred_atom_coords_local.shape[0]}) "
                f"!= true_atom_coords batch dim ({true_atom_coords_local.shape[0]})"
            )
        if true_coords_resolved_mask_local.shape != true_atom_coords_local.shape[:2]:
            raise ValueError(
                f"true_coords_resolved_mask shape {tuple(true_coords_resolved_mask_local.shape)} "
                f"must match true_atom_coords[:2] {tuple(true_atom_coords_local.shape[:2])}"
            )

        # --- Step 1: Compute target values and masks ---
        with torch.no_grad():
            mask_frame_true = frame_resolved_mask_local

            # Compute frames for true coords (use DTensor wrapper for serial-consistent results)
            frames_idx_true_dt, _ = compute_frame_pred(
                true_atom_coords,
                feats["frames_idx"],
                feats,
                multiplicity,
                resolved_mask=true_coords_resolved_mask,
            )
            frames_idx_true = frames_idx_true_dt.to_local()

            true_atom_coords_reshaped = true_atom_coords_local.reshape(B_batch, multiplicity, -1, 3)
            transpose_comm = TransposeComm(dist_manager.group["cp"], group_layout)
            true_coords_transformed, mask_collinear_true = _express_coordinate_in_frame_distributed(
                true_atom_coords_reshaped,
                frames_idx_true[:, :, :, 0],
                frames_idx_true[:, :, :, 1],
                frames_idx_true[:, :, :, 2],
                dist_manager,
                group_layout,
                transpose_comm=transpose_comm,
            )
            mask_collinear_true = mask_collinear_true * token_pad_mask_local[:, None, :]

            # Compute frames for pred coords (use DTensor wrapper for serial-consistent results)
            frames_idx_pred_dt, _ = compute_frame_pred(
                pred_atom_coords,
                feats["frames_idx"],
                feats,
                multiplicity,
            )
            frames_idx_pred = frames_idx_pred_dt.to_local()

            pred_atom_coords_reshaped = pred_atom_coords_local.reshape(B_batch, multiplicity, -1, 3)
            pred_coords_transformed, mask_collinear_pred = _express_coordinate_in_frame_distributed(
                pred_atom_coords_reshaped,
                frames_idx_pred[:, :, :, 0],
                frames_idx_pred[:, :, :, 1],
                frames_idx_pred[:, :, :, 2],
                dist_manager,
                group_layout,
                transpose_comm=transpose_comm,
            )
            mask_collinear_pred = mask_collinear_pred * token_pad_mask_local[:, None, :]

            # Compute target PAE distances
            target_pae = torch.sqrt(((true_coords_transformed - pred_coords_transformed) ** 2).sum(-1) + PAE_DIST_EPS)

            # Compute bin indices for cross-entropy
            bin_index = torch.clamp(torch.floor(target_pae * num_bins / max_dist).long(), max=(num_bins - 1))

            # Build pair mask: gather resolved status for all 3 frame atoms.
            # Each diffusion sample has its own resolved mask (symmetry_correction
            # can produce different masks per sample), so preserve the per-sample
            # variation rather than collapsing to sample 0.
            resolved_reshaped = true_coords_resolved_mask_local.reshape(B_batch, multiplicity, -1)
            token_pad_mask_bool = token_pad_mask_local[:, None, :].bool()
            N_token_local = frames_idx_true.shape[-2]

            frames_masked = frames_idx_true.masked_fill(~token_pad_mask_bool.unsqueeze(-1), atom_offset)
            requires_global_gather = torch.any(
                (frames_masked < atom_offset) | (frames_masked >= atom_offset + N_atom_local)
            ).to(dtype=torch.int32)
            if N_token_local != N_atom_local:
                requires_global_gather = torch.ones_like(requires_global_gather)
            dist.all_reduce(requires_global_gather, op=dist.ReduceOp.MAX, group=dist_manager.group["cp"])

            if requires_global_gather.item() == 1:
                resolved_flat = resolved_reshaped.reshape(B_batch * multiplicity, N_atom_local, 1)
                index_flat = frames_idx_true.reshape(B_batch * multiplicity, N_token_local, 3)
                gathered = ring_gather_coordinate(resolved_flat, index_flat, dist_manager, group_layout)
                # gathered: (B*mult, N_token_local, 1, 3) → squeeze → (B_batch, mult, N_token_local, 3)
                frame_resolved_abc = gathered.squeeze(-2).reshape(B_batch, multiplicity, N_token_local, 3)
            else:
                frames_local = frames_idx_true - atom_offset
                frame_resolved_abc = torch.stack(
                    [torch.gather(resolved_reshaped, dim=2, index=frames_local[:, :, :, k]) for k in range(3)],
                    dim=-1,
                )

            b_true_resolved_mask_local = frame_resolved_abc[:, :, :, 1]

            # Exchange masks with transpose peer for column (j) dimension
            b_true_resolved_mask_j = transpose_comm.enqueue_to_dispatch(b_true_resolved_mask_local.contiguous())
            transpose_comm.wait_until_finished()
            token_pad_mask_j = transpose_comm.enqueue_to_dispatch(token_pad_mask_local.contiguous())
            transpose_comm.wait_until_finished()

            pair_mask = (
                mask_frame_true[:, None, :, None]
                * mask_collinear_true[:, :, :, None]
                * mask_collinear_pred[:, :, :, None]
                * b_true_resolved_mask_j[:, :, None, :]
                * token_pad_mask_local[:, None, :, None]
                * token_pad_mask_j[:, None, None, :]
            )

            # Compute local mask sum and reduce across CP ranks.
            sum_mask_local = pair_mask.sum(dim=(-2, -1))
            sum_mask_global_cp = sum_mask_local.clone()
            # Reduce across CP ranks first
            dist.all_reduce(sum_mask_global_cp, op=dist.ReduceOp.SUM, group=dist_manager.group["cp"])

        # --- Step 2: Compute loss with gradients using global normalization ---
        with torch.enable_grad():
            pred_pae_local_grad = pred_pae_local.detach().requires_grad_(pred_pae.requires_grad)
            log_softmax_pae = F.log_softmax(pred_pae_local_grad, dim=-1)
            target_log_prob = torch.gather(log_softmax_pae, dim=-1, index=bin_index.unsqueeze(-1)).squeeze(-1)
            errors = -target_log_prob

            # Local sum of masked errors
            masked_errors = errors * pair_mask
            sum_errors_local = masked_errors.sum(dim=(-2, -1))

            # Normalize by CP-global mask sum per sample for gradient computation.
            loss_per_sample = sum_errors_local / (sum_mask_global_cp + PAE_LOSS_DENOM_EPS)
            loss_local = loss_per_sample.mean() / dp_size

        # --- Step 3: Compute  global loss---
        # Strategy: Reduce errors across CP, normalize per sample, then average across DP ranks.
        #
        # This matches the serial definition:
        #   loss = mean_bm(sum_ij(errors) / sum_ij(mask))
        #
        # Reduce errors across CP ranks
        sum_errors_global_cp = sum_errors_local.detach().clone()
        dist.all_reduce(sum_errors_global_cp, op=dist.ReduceOp.SUM, group=dist_manager.group["cp"])

        # Normalize per sample, then average across local samples
        loss_per_sample_final = sum_errors_global_cp / (sum_mask_global_cp + PAE_LOSS_DENOM_EPS)
        loss_scalar = loss_per_sample_final.mean()

        # Average across DP ranks to match global batch mean
        dist.all_reduce(loss_scalar, op=dist.ReduceOp.SUM, group=dist_manager.group["dp"])
        loss_scalar = loss_scalar / dp_size

        # Save for backward
        ctx.save_for_backward(pred_pae_local_grad, loss_local)
        ctx.device_mesh = device_mesh
        ctx.pred_pae_shape = pred_pae.shape
        ctx.pred_pae_stride = pred_pae.stride()
        ctx.pred_pae_placements = pred_pae.placements
        # Create replicated DTensor for loss
        loss_dtensor = DTensor.from_local(
            loss_scalar,
            device_mesh=device_mesh,
            placements=(Replicate(), Replicate(), Replicate()),
            shape=loss_scalar.shape,
            stride=loss_scalar.stride(),
        )

        return loss_dtensor

    @staticmethod
    @torch.amp.custom_bwd(device_type="cuda")
    def backward(ctx, grad_loss: DTensor):
        pred_pae_local_grad, loss_local = ctx.saved_tensors
        if not pred_pae_local_grad.requires_grad:
            return (
                None,  # pred_pae
                None,  # pred_atom_coords (non-differentiable through floor)
                None,  # true_atom_coords
                None,  # true_coords_resolved_mask
                None,  # feats
                None,  # comm
                None,  # dist_manager
                None,  # group_layout
                None,  # multiplicity
                None,  # max_dist
            )

        grad_loss_scalar = grad_loss.to_local()
        (grad_pred_pae_local,) = torch.autograd.grad(
            outputs=[loss_local],
            inputs=[pred_pae_local_grad],
            grad_outputs=[grad_loss_scalar],
            retain_graph=False,
        )

        # Flatten (B_batch, mult, ...) back to (B_batch*mult, ...) to match DTensor layout
        grad_pred_pae_local = grad_pred_pae_local.flatten(0, 1)

        # Create DTensor gradient for pred_pae
        grad_pred_pae = DTensor.from_local(
            grad_pred_pae_local,
            device_mesh=ctx.device_mesh,
            placements=ctx.pred_pae_placements,
            shape=ctx.pred_pae_shape,
            stride=ctx.pred_pae_stride,
        )

        return (
            grad_pred_pae,  # pred_pae
            None,  # pred_atom_coords (non-differentiable through floor)
            None,  # true_atom_coords
            None,  # true_coords_resolved_mask
            None,  # feats
            None,  # comm
            None,  # dist_manager
            None,  # group_layout
            None,  # multiplicity
            None,  # max_dist
        )


def all_reduce_dist_mat_argmin(
    dist_mat: Tensor,
    group_reduce_dist_mat_row: dist.ProcessGroup,
    order: int = 3,
    inplace: bool = False,
) -> Tensor:
    """Perform distributed argmin operation along given dimension up to the k-th smallest element.

    This function expects dist_mat to be sharded across the process group, such that:
    device mesh = [[  0,   1],
                   [  2,   3]]
    dist_mat    = [[ d00, d01],
                   [ d10, d11]]

    The returned global_argmin will use a different sharding strategy, where it is sharded column-wise and replicated row-wise.
    global_argmin = [[ a0,  a0],
                     [ a1,  a1]]

    Args:
        dist_mat (torch.Tensor): Local shard of distance matrix to perform argmin operation on. Shape = (..., n, n).
        group_reduce_dist_mat_row (dist.ProcessGroup): cp axis 0 process group to reduce the argmin operation to.
        order (int): Number of the smallest elements to find iteratively. Default is 3 to find the closest 3 atoms to construct a frame.
        inplace (bool): If True, perform the operation on dist_mat in place. Default is False.
    Returns:
        torch.Tensor: global argmin indices of the first n smallest elements in shape = (..., order)
    """
    if order > dist_mat.shape[-1] * group_reduce_dist_mat_row.size():
        raise ValueError(
            "order must be less than or equal to the number of elements in the group but got order = {}, dist_mat.shape[dim] = {}, group_reduce_dist_mat_row.size() = {}".format(
                order, dist_mat.shape[-1], group_reduce_dist_mat_row.size()
            )
        )
    if dist_mat.shape[-2] != dist_mat.shape[-1]:
        raise ValueError("distance matrix must be square but got shape = {}".format(dist_mat.shape))

    if not inplace:  # save memory through inplace operation
        dist_mat = dist_mat.clone()

    out_tensor = []
    for _ in range(order):
        # find local min and argmin
        local_min, local_argmin = torch.min(dist_mat, dim=-1)  # shape = (B, N)

        # reduce to global min and check if it is the same as the local min
        global_min = local_min
        dist.all_reduce(global_min, op=dist.ReduceOp.MIN, group=group_reduce_dist_mat_row)

        # locate global min locally
        is_global_min = global_min.unsqueeze(-1) == dist_mat  # shape = (B, N, N)
        has_global_min = is_global_min.any(dim=-1)  # shape = (B, N)

        # offset local argmin by group rank to get global argmin
        global_argmin = local_argmin.clone() + dist.get_rank(group_reduce_dist_mat_row) * dist_mat.shape[-1]

        # mask non-global-min argmin from broadcasting
        global_argmin[~has_global_min] = 0

        # broadcast the output tensor to the process group
        dist.all_reduce(global_argmin, op=dist.ReduceOp.SUM, group=group_reduce_dist_mat_row)  # shape = (B, N)

        # aggregate indices
        out_tensor.append(global_argmin)

        # mask the global min with inf for the next order
        dist_mat[is_global_min] = torch.inf

    out_tensor = torch.stack(out_tensor, dim=-1)
    return out_tensor


def ring_dist_mat_argmin(
    dist_mat: Tensor,
    group_cp: dist.ProcessGroup,
    group_reduce_dist_mat_row: dist.ProcessGroup,
    group_layout: LayoutMap,
    order: int = 3,
) -> Tensor:
    """Perform distributed argmin operation up to the k-th smallest element. The distance matrix is assumed to be square and symmetric.

    This function expects dist_mat to be sharded across the process group, such that:
    device mesh   = [[  0,   1],
                     [  2,   3]]
    dist_mat      = [[ d00, d01],
                     [ d10, d11]]

    The returned global_argmin will use a different sharding strategy, where it is sharded column-wise and replicated row-wise.
    global_argmin = [[ a0,  a0],
                     [ a1,  a1]]

    Args:
        dist_mat (torch.Tensor): Local shard of distance matrix to perform argmin operation on. Shape = (..., n, n).
        group_cp (dist.ProcessGroup): Context parallelism process group
        group_reduce_dist_mat_row (dist.ProcessGroup): Row-wise reduction process group
        group_layout (LayoutMap): Layout map of the process group
        order (int): Number of the smallest elements to find iteratively. Default is 3 to find the closest 3 atoms to construct a frame.
    Returns:
        torch.Tensor: global argmin indices of the first k smallest elements in shape = (..., n, k)
    """
    cp_axis_0_group = group_reduce_dist_mat_row
    if order > dist_mat.shape[-1] * cp_axis_0_group.size():
        raise ValueError(
            "order must be less than or equal to the number of elements in the group but got order = {}, dist_mat.shape[dim] = {}, cp_axis_0_group.size() = {}".format(
                order, dist_mat.shape[-1], cp_axis_0_group.size()
            )
        )

    if dist_mat.shape[-2] != dist_mat.shape[-1]:
        raise ValueError("distance matrix must be square but got shape = {}".format(dist_mat.shape))

    # setup for communication
    rank_coords = group_layout.unravel(dist.get_rank(group_cp))
    topk_comm = One2OneComm(
        group_cp,
        rank_send_to=get_group_rank_from_axial_shift(rank_coords, 1, -1, group_layout),
        rank_recv_from=get_group_rank_from_axial_shift(rank_coords, 1, 1, group_layout),
    )
    topk_idx_comm = deepcopy(topk_comm)

    # find local topk and topk_idx
    max_order = min(order, dist_mat.shape[-1])
    local_topk, local_topk_idx = torch.topk(dist_mat, k=max_order, dim=-1, largest=False)  # shape = (..., n, k)

    # offset local topk_idx by group rank to get global topk_idx
    global_topk = local_topk
    global_topk_idx = local_topk_idx + dist.get_rank(group_reduce_dist_mat_row) * dist_mat.shape[-1]

    # send out to the next rank for the first time
    current_topk = topk_comm.enqueue_to_dispatch(global_topk)
    current_topk_idx = topk_idx_comm.enqueue_to_dispatch(global_topk_idx)
    topk_comm.wait_until_finished()
    topk_idx_comm.wait_until_finished()

    for step in range(cp_axis_0_group.size() - 1):
        # overlap communication with computation by sending out to the next rank
        if step != cp_axis_0_group.size() - 2:
            next_topk = topk_comm.enqueue_to_dispatch(current_topk)
            next_topk_idx = topk_idx_comm.enqueue_to_dispatch(current_topk_idx)

        # concatenate the received topk and topk_idx
        global_topk = torch.cat([global_topk, current_topk], dim=-1)  # shape = (..., n, 2k)
        global_topk_idx = torch.cat([global_topk_idx, current_topk_idx], dim=-1)  # shape = (..., n, 2k)

        # find the largest k values in global_topk and the corresponding indices in global_topk_idx
        max_order = min(order, global_topk.shape[-1])
        topk_values, topk_indices = torch.topk(global_topk, k=max_order, dim=-1, largest=False)  # shape = (..., n, k)

        # select the topk_indices by value from global_topk_idx
        global_topk = topk_values
        global_topk_idx = global_topk_idx.gather(dim=-1, index=topk_indices)  # shape = (..., n, k)

        # receive from the previous rank
        if step != cp_axis_0_group.size() - 2:
            topk_comm.wait_until_finished()
            topk_idx_comm.wait_until_finished()
            current_topk = next_topk
            current_topk_idx = next_topk_idx

    return global_topk_idx


def ring_gather_coordinate(
    coordinate: Tensor,
    global_argmin: Tensor,
    dist_manager: DistributedManager,
    group_layout: LayoutMap,
) -> Tensor:
    """Distributed version of torch.gather in ring topology for coordinate gathering.

    example sharding strategy on N_atoms=2, world_size=4
    device mesh = [[  0,   1],
                   [  2,   3]]

    the gather operation is aggregated through a ring topology by rolling up the coordinates while offsetting the global_argmin by the axial shift.

    step = 0
    coordinate = [[ c0,  c0],   global_argmin = [[ g0,  g0],
                  [ c1,  c1]]                    [ g1,  g1]]

    step = 1 (roll up)
    coordinate = [[ c1,  c1],   global_argmin = [[ g0,  g0],
                  [ c0,  c0]]                    [ g1,  g1]]

    outputs follow same sharding:
    gathered_coords = [[ c0, c0 ],
                       [ c1, c1 ]]

    Args:
        coordinate (torch.Tensor): Sharded coordinates of shape = (B, n_atoms_local, D) to be gathered from.
        global_argmin (torch.Tensor): Sharded global argmin atom indices of shape = (B, n_tokens_local, order).
        dist_manager (DistributedManager): Distributed manager
        group_layout (LayoutMap): Layout map of the process group

    Returns:
        torch.Tensor: gathered coordinates; shape = (B, n_tokens_local, D, order)
    """
    ring_size = group_layout.shape[0]
    bs, n_atoms_local, coord_dim = coordinate.shape
    _, n_tokens_local, order = global_argmin.shape

    batch_idx = (
        torch.arange(bs, device=coordinate.device).view(bs, 1, 1).expand(-1, n_tokens_local, order)
    )  # shape = (B, n_tokens_local, order)

    rank_coords = group_layout.unravel(dist.get_rank(dist_manager.group["cp"]))
    comm = One2OneComm(
        dist_manager.group["cp"],
        rank_send_to=get_group_rank_from_axial_shift(rank_coords, 0, -1, group_layout),
        rank_recv_from=get_group_rank_from_axial_shift(rank_coords, 0, 1, group_layout),
    )

    gathered_coords = torch.zeros(
        bs, n_tokens_local, order, coord_dim, device=coordinate.device, dtype=coordinate.dtype
    )
    for step in range(ring_size):
        if step + 1 != ring_size:
            next_coordinate = comm.enqueue_to_dispatch(coordinate.contiguous())

        idx_range = (rank_coords[0] + step) % ring_size
        shard_start = idx_range * n_atoms_local
        is_argmin_local = (shard_start <= global_argmin) & (global_argmin < (idx_range + 1) * n_atoms_local)

        local_argmin = global_argmin - shard_start
        local_argmin = local_argmin.masked_fill(~is_argmin_local, 0)

        gathered = coordinate[batch_idx, local_argmin]
        gathered = gathered.masked_fill(~is_argmin_local[..., None], 0)
        gathered_coords += gathered

        if step + 1 != ring_size:
            comm.wait_until_finished()
            coordinate = next_coordinate

    return gathered_coords.permute(0, 1, 3, 2)  # shape = (B, n_tokens_local, D, order)


def _gather_frame_coords(
    atom_coords: Tensor,
    frame_atoms: Tensor,
    *,
    local_only: bool,
    atom_offset: int,
    n_tokens_local: int,
    dist_manager: DistributedManager,
    group_layout: LayoutMap,
) -> Tensor:
    """Gather frame atom coordinates for local or global frame indices.

    Args:
        atom_coords: Local atom coordinates with shape [B, mult, N_atom_local, 3].
        frame_atoms: Frame atom indices with shape [B, mult, N_token_local, order].
        local_only: If True, gather directly from local shard using atom_offset.
        atom_offset: Global atom index offset for this shard.
        n_tokens_local: Number of local tokens for this shard.
        dist_manager: DistributedManager for collective communications.
        group_layout: LayoutMap for CP mesh ring gather.

    Returns:
        Gathered frame coordinates with shape [B, mult, N_token_local, 3].
    """
    batch, multiplicity = atom_coords.shape[0], atom_coords.shape[1]
    n_atoms_local = atom_coords.shape[2]
    device = atom_coords.device

    if local_only:
        local_idx = frame_atoms - atom_offset
        batch_indices0 = torch.arange(batch, device=device)[:, None, None]
        batch_indices1 = torch.arange(multiplicity, device=device)[None, :, None]
        return atom_coords[batch_indices0, batch_indices1, local_idx]

    # Flatten (batch, multiplicity) into a single batch dimension for ring gather.
    # ring_gather_coordinate expects (B, n_atoms_local, D) and (B, n_tokens_local, order).
    coords_flat = atom_coords.reshape(batch * multiplicity, n_atoms_local, 3)
    idx_flat = frame_atoms.reshape(batch * multiplicity, n_tokens_local, 1)
    gathered = ring_gather_coordinate(coords_flat, idx_flat, dist_manager, group_layout)
    # Unflatten back to (batch, multiplicity) to match the caller's expectations.
    return gathered.squeeze(-1).reshape(batch, multiplicity, n_tokens_local, 3)


def _fully_distributed_compute_frame_pred(
    pred_atom_coords: Tensor,
    frames_idx_true: Tensor,
    feats: dict[str, Tensor],
    multiplicity: int,
    comm: One2OneComm,
    dist_manager: DistributedManager,
    group_layout: LayoutMap,
    resolved_mask: Optional[Tensor] = None,
    inference: bool = False,
    return_frames_expanded: bool = False,
) -> tuple[Tensor, Tensor, Tensor] | tuple[Tensor, Tensor]:
    """Recompute the frames for non-polymer over 3 atoms given the predicted atom coordinates.

    .. deprecated::
        Unused — superseded by ``compute_frame_pred`` (the DTensor wrapper that
        gathers inputs and delegates to ``_compute_frame_pred``).  Retained for
        reference only.  Contains known resolved_mask indexing bugs that are
        NOT fixed here; see _compute_frame_pred for the corrected version.

    example sharding strategy on N_atoms=2, world_size=4
    device mesh = [[  0,   1],
                   [  2,   3]]

    inputs should follow a sharding strategy below.
    pred_coords = [[ c0,  c0],   frames_idx = [[ f0,  f0],   resolved_mask = [[m0, m0],
                   [ c1,  c1]]                 [ f1,  f1]]                    [m1, m1]]

    outputs follow same sharding:
    frames_idx_pred = [[ f0,  f0],   mask_collinear = [[m0, m0],   frames_expanded = [[c0, c0],
                       [ f1,  f1]]                     [m1, m1]]                      [c1, c1]]

    Args:
        pred_atom_coords (torch.Tensor): Predicted atom coordinates of shape = (B, n_atoms_per_shard, 3).
        frames_idx_true (torch.Tensor): True frames indices of shape = (B, n_atoms_per_shard, order).
        feats (dict): Dictionary of feature tensors
        multiplicity (int): Multiplicity of the predicted atom coordinates
        comm (One2OneComm): Communication class for sending and receiving tensors
        dist_manager (DistributedManager): Distributed manager
        group_layout (LayoutMap): Layout map for context parallelism
        resolved_mask (torch.Tensor, optional): Resolved mask; shape = (B, n_atoms_per_shard). Defaults to None to use atom_resolved_mask and atom_pad_mask in feats.
        inference (bool, optional): Whether to use inference mode which skips resolved_mask. Defaults to False.
        return_frames_expanded (bool, optional): Whether to return the expanded frames for unittest purposes. Defaults to False.

    Returns:
        torch.Tensor: Updated frames indices; shape = (B, N, order)
        torch.Tensor: Mask for collinear or overlapping atoms in the frame; shape = (B, N, order)
        optional torch.Tensor: The closest 3 atom coordinates for the frames; shape = (B, N, order, 3). Returned only if return_frames_expanded is True.
    """
    # group settings
    group_axis_1 = dist_manager.subgroups["cp"][0]
    group_rank_1 = dist.get_rank(group_axis_1)

    # extract necessary features
    asym_id_token = feats["asym_id"]
    asym_id_atom = torch.bmm(feats["atom_to_token"].float(), asym_id_token.unsqueeze(-1).float()).squeeze(-1)
    B, N, _ = pred_atom_coords.shape

    pred_atom_coords = pred_atom_coords.reshape(B // multiplicity, multiplicity, -1, 3)
    frames_idx_pred = (
        frames_idx_true.clone().repeat_interleave(multiplicity, 0).reshape(B // multiplicity, multiplicity, -1, 3)
    )

    frames_expanded = []

    # Iterate through the batch and update the frames for non-polymers
    for i, pred_atom_coord in enumerate(pred_atom_coords):
        # Gather reference atom coordinates per token per order, i.e. frame per token
        pred_atom_coords_sample = pred_atom_coords[i, :, :, :]  # atom coordinates; shape = (multiplicity, N_atoms, 3)
        idx = (
            frames_idx_pred[i] - pred_atom_coords.shape[-2] * group_rank_1
        )  # token to **local** atom indices; shape = (multiplicity, N_tokens)
        idx = idx.masked_fill(
            feats["is_nonpolymer_with_frame"][i][:, None], 0
        )  # reset frames from non-polymer chains with fewer than 3 atoms to 0

        N_tokens = idx.shape[-2]
        batch_idx = (
            torch.arange(multiplicity).view(multiplicity, 1, 1).expand(-1, N_tokens, 3)
        )  # shape = (multiplicity, N_tokens, 3)
        frame_expanded_sample = pred_atom_coords_sample[batch_idx, idx].transpose(
            -1, -2
        )  # shape = (multiplicity, N_atoms, 3, orders)

        assert frame_expanded_sample.shape == (multiplicity, N_tokens, 3, 3)
        frame_expanded_sample = frame_expanded_sample.masked_fill(
            feats["is_nonpolymer_with_frame"][i][None, :, None, None], 0
        )  # reset frames from non-polymer chains with fewer than 3 atoms to 0

        # Gather unique asym_ids from all ranks
        asym_ids_unique = [set() for _ in range(group_axis_1.size())]
        dist.all_gather_object(asym_ids_unique, set(asym_id_token[i].tolist()), group=group_axis_1)
        asym_ids_unique = sorted(set.union(*asym_ids_unique))

        for id in asym_ids_unique:
            mask_chain_token = (asym_id_token[i] == id) * feats["token_pad_mask"][i]
            mask_chain_atom = (asym_id_atom[i] == id) * feats["atom_pad_mask"][i]
            mask_chain_token = mask_chain_token.bool()
            mask_chain_atom = mask_chain_atom.bool()

            # Check if the chain satisfies the criteria for frame recomputation
            #  1. is a non-polymer
            #  2. has at least 3 atoms
            # TODO: streamline this with is_nonpolymer_with_frame
            num_tokens = mask_chain_token.sum()
            num_atoms = mask_chain_atom.sum()
            dist.all_reduce(num_tokens, op=dist.ReduceOp.SUM, group=group_axis_1)
            dist.all_reduce(num_atoms, op=dist.ReduceOp.SUM, group=group_axis_1)

            mol_type = feats["mol_type"][i, mask_chain_token].unique()
            assert len(mol_type) <= 1, "all chains in the batch must have the same mol_type"

            is_target_mol_type = (
                (mol_type.item() == const.chain_type_ids["NONPOLYMER"]) and (num_atoms.item() > 3)
                if len(mol_type) > 0
                else False
            )
            is_target_mol_type = torch.tensor(is_target_mol_type, device=mol_type.device)
            dist.all_reduce(is_target_mol_type, op=dist.ReduceOp.SUM, group=group_axis_1)
            if ~is_target_mol_type:
                continue
            assert (
                num_tokens.item() == num_atoms.item()
            ), "num_tokens and num_atoms must be the same for non-polymers, got num_tokens = {}, num_atoms = {}".format(
                num_tokens.item(), num_atoms.item()
            )

            # Compute all-to-all atom distance matrix, including those that are not part of the chain
            pred_atom_coord_i = pred_atom_coord
            pred_atom_coord_j = comm.enqueue_to_dispatch(pred_atom_coord_i)
            comm.wait_until_finished()
            dist_mat = torch.cdist(pred_atom_coord_i, pred_atom_coord_j)  # shape = (multiplicity, N, N)

            # Restrict neighborhood frame atom search to
            #  1. atoms that are not padding/are resolved, and
            #  2. atoms that are part of the chain
            if inference:
                resolved_mask_i = feats["atom_pad_mask"][i]
            elif resolved_mask is None:
                resolved_mask_i = feats["atom_resolved_mask"][i]
                resolved_mask_i = (
                    resolved_mask_i * feats["atom_pad_mask"][i]
                )  # apply atom_pad_mask for padding in context parallelism
            else:
                resolved_mask_i = resolved_mask[i]

            resolved_mask_j = comm.enqueue_to_dispatch(resolved_mask_i)
            comm.wait_until_finished()
            resolved_pair = (1 - resolved_mask_i[:, None] * resolved_mask_j[None, :]).to(
                torch.float32
            )  # shape = (N, N)
            resolved_pair[resolved_pair == 1] = torch.inf

            mask_chain_atom_i = mask_chain_atom
            mask_chain_atom_j = comm.enqueue_to_dispatch(mask_chain_atom_i)
            comm.wait_until_finished()
            mask_chain_atom_pair = 1 - (mask_chain_atom_i[:, None] * mask_chain_atom_j[None, :]).to(torch.float)
            mask_chain_atom_pair[mask_chain_atom_pair == 1] = torch.inf

            # Sort the atoms by distance
            masked_dist_mat = dist_mat + resolved_pair + mask_chain_atom_pair
            global_argmin = ring_dist_mat_argmin(
                masked_dist_mat,
                order=3,
                group_cp=dist_manager.group["cp"],
                group_reduce_dist_mat_row=dist_manager.subgroups["cp"][1],
                group_layout=group_layout,
            )  # shape = (multiplicity, N_token, order)
            atom_to_ref_coords = ring_gather_coordinate(
                pred_atom_coord_i,
                global_argmin,
                dist_manager,
                group_layout=group_layout,
            )  # shape = (multiplicity, N_token, 3, order)

            # Map reference atom repr to token repr
            # non-polymer has one token per atom so this is a mapping instead of aggregation
            atom_to_token = feats["atom_to_token"][i].float()
            token_to_ref_coords = torch.einsum(
                "miab,ij->mjab", atom_to_ref_coords, atom_to_token
            )  # shape = (multiplicity, N_token, 3, order)
            frames = torch.einsum(
                "mib,ij->mjb", global_argmin.float(), atom_to_token
            ).long()  # shape = (multiplicity, N_token, order)

            # pass reference frames for non-polymer chains
            frame_reorder_index = torch.tensor([1, 0, 2], device=token_to_ref_coords.device)
            frames_idx_pred[i, :, :, :] = (
                frames[:, :, frame_reorder_index] * mask_chain_token[None, :, None]
                + frames_idx_pred[i, :, :, :] * ~mask_chain_token[None, :, None]
            )  # shape = (multiplicity, N_token, order)

            # NOTE: discrepancy in padding will result in different indexing compared to
            #  single-device implementation; compare frames_expanded instead of frames_idx_pred

            # pass reference atom coordinates for non-polymer chains
            frame_expanded_sample = (
                token_to_ref_coords[:, :, :, frame_reorder_index] * mask_chain_token[None, :, None, None]
                + frame_expanded_sample * ~mask_chain_token[None, :, None, None]
            )  # shape = (multiplicity, N_token, 3, order)

        # append per sample
        frames_expanded.append(frame_expanded_sample)

    # concatenate per sample
    frames_expanded = torch.cat(frames_expanded, dim=0).reshape(-1, 3, 3)
    frames_expanded = frames_expanded.transpose(1, 2)  # shape = (..., order, 3)

    # Compute masks for collinear or overlapping atoms in the frame
    mask_collinear_pred = compute_collinear_mask(
        frames_expanded[:, 1] - frames_expanded[:, 0],
        frames_expanded[:, 1] - frames_expanded[:, 2],
    ).reshape(B // multiplicity, multiplicity, -1)

    if return_frames_expanded:
        return frames_idx_pred, mask_collinear_pred * feats["token_pad_mask"][:, None, :], frames_expanded
    else:
        return frames_idx_pred, mask_collinear_pred * feats["token_pad_mask"][:, None, :]


def _compute_frame_pred(
    pred_atom_coords: Tensor,
    frames_idx_true: Tensor,
    feats: dict[str, Tensor],
    asym_id_atom: Tensor,
    multiplicity: int,
    resolved_mask: Optional[Tensor] = None,
    inference: bool = False,
) -> tuple[Tensor, Tensor]:
    """Private tensor implementation compatible with sparse/padded atom indexing of the reference compute_frame_pred in src/boltz/model/layers/confidence_utils.py."""
    # Disable autocast to match serial compute_frame_pred which runs inside
    # torch.amp.autocast("cuda", enabled=False).  Without this, bmm and pow
    # are affected by autocast, producing different intermediate precision
    # than serial and potentially different frame assignments for borderline
    # non-polymer ligand chains.
    with torch.amp.autocast("cuda", enabled=False):
        return _compute_frame_pred_impl(
            pred_atom_coords,
            frames_idx_true,
            feats,
            asym_id_atom,
            multiplicity,
            resolved_mask,
            inference,
        )


def _compute_frame_pred_impl(
    pred_atom_coords: Tensor,
    frames_idx_true: Tensor,
    feats: dict[str, Tensor],
    asym_id_atom: Tensor,
    multiplicity: int,
    resolved_mask: Optional[Tensor] = None,
    inference: bool = False,
) -> tuple[Tensor, Tensor]:
    """Implementation of _compute_frame_pred, called with autocast disabled."""
    asym_id_token = feats["asym_id"]
    B, _, _ = pred_atom_coords.shape
    if B % multiplicity != 0:
        raise ValueError(f"pred_atom_coords batch dim ({B}) must be divisible by multiplicity ({multiplicity})")
    if resolved_mask is not None and resolved_mask.shape != pred_atom_coords.shape[:2]:
        raise ValueError(
            f"resolved_mask shape {tuple(resolved_mask.shape)} must match "
            f"pred_atom_coords[:2] {tuple(pred_atom_coords.shape[:2])}"
        )
    pred_atom_coords = pred_atom_coords.reshape(B // multiplicity, multiplicity, -1, 3)
    # resolved_mask arrives as (B*mult, N_atom).  Reshape to (B_batch, mult,
    # N_atom) so that each diffusion sample's per-sample resolved mask is
    # preserved (symmetry_correction can produce different masks per sample).
    if resolved_mask is not None:
        resolved_mask = resolved_mask.reshape(B // multiplicity, multiplicity, -1)
    frames_idx_pred = (
        frames_idx_true.clone().repeat_interleave(multiplicity, 0).reshape(B // multiplicity, multiplicity, -1, 3)
    )

    for i, pred_atom_coord in enumerate(pred_atom_coords):
        for id in torch.unique(asym_id_token[i]):
            mask_chain_token = (asym_id_token[i] == id) * feats["token_pad_mask"][i]
            mask_chain_atom = (asym_id_atom[i] == id) * feats["atom_pad_mask"][i]
            num_tokens = int(mask_chain_token.sum().item())
            num_atoms = int(mask_chain_atom.sum().item())

            mol_types = feats["mol_type"][i, mask_chain_token.bool()]

            # sanity check: all chains in the batch must have the same mol_type
            mol_type_unique = mol_types.unique()
            assert (
                mol_type_unique.numel() <= 1
            ), f"all chains in the batch must have the same mol_type but got {mol_type_unique}"  # sanity check

            # skip frame reassignment if the chain is not a non-polymer or has fewer than 3 atoms
            if mol_type_unique.item() != const.chain_type_ids["NONPOLYMER"] or num_atoms < 3:
                continue

            # sanity check: num_atoms = num_tokens for non-polymers
            assert (
                num_atoms == num_tokens
            ), "num_atoms and num_tokens must be the same for non-polymers, got num_atoms = {}, num_tokens = {}".format(
                num_atoms, num_tokens
            )

            chain_atom_indices = torch.nonzero(mask_chain_atom.bool(), as_tuple=False).squeeze(-1)
            chain_atom_coords = pred_atom_coord[:, chain_atom_indices]
            dist_mat = ((chain_atom_coords[:, None, :, :] - chain_atom_coords[:, :, None, :]) ** 2).sum(-1) ** 0.5

            if inference:
                resolved_pair = 1 - (
                    feats["atom_pad_mask"][i][chain_atom_indices][None, :]
                    * feats["atom_pad_mask"][i][chain_atom_indices][:, None]
                ).to(torch.float32)
                resolved_pair[resolved_pair == 1] = torch.inf
                indices = torch.sort(dist_mat + resolved_pair, axis=2).indices
            else:
                if resolved_mask is None:
                    # atom_resolved_mask is (B_batch, N_atom); expand to
                    # (B_batch, mult, N_atom) so indexing is uniform.
                    resolved_mask = feats["atom_resolved_mask"][:, None, :].expand(-1, multiplicity, -1)
                # resolved_mask[i]: (mult, N_atom)
                rm_chain = resolved_mask[i][:, chain_atom_indices]  # (mult, N_chain)
                resolved_pair = 1 - (rm_chain[:, None, :] * rm_chain[:, :, None]).to(torch.float32)
                resolved_pair[resolved_pair == 1] = torch.inf
                indices = torch.sort(dist_mat + resolved_pair, axis=2).indices

            frames_local = torch.cat(
                [
                    indices[:, :, 1:2],
                    indices[:, :, 0:1],
                    indices[:, :, 2:3],
                ],
                dim=2,
            )
            frames = chain_atom_indices[frames_local]
            frames_idx_pred[i, :, mask_chain_token.bool(), :] = frames

    frames_expanded = pred_atom_coords[
        torch.arange(0, B // multiplicity, 1)[:, None, None, None].to(frames_idx_pred.device),
        torch.arange(0, multiplicity, 1)[None, :, None, None].to(frames_idx_pred.device),
        frames_idx_pred,
    ].reshape(-1, 3, 3)

    mask_collinear_pred = compute_collinear_mask(
        frames_expanded[:, 1] - frames_expanded[:, 0],
        frames_expanded[:, 1] - frames_expanded[:, 2],
    ).reshape(B // multiplicity, multiplicity, -1)
    return frames_idx_pred, mask_collinear_pred * feats["token_pad_mask"][:, None, :]


def compute_frame_pred(
    pred_atom_coords: DTensor,
    frames_idx_true: DTensor,
    feats: dict[str, DTensor],
    multiplicity: int,
    resolved_mask: Optional[DTensor] = None,
    inference: bool = False,
) -> tuple[DTensor, DTensor]:
    """Distributed wrapper around `compute_frame_pred` for DTensors.

    Gathers DTensor inputs to local tensors, runs the serial implementation,
    then reshapes and redistributes outputs to match DTensor placements.

    Args:
        pred_atom_coords: DTensor of shape (batch * multiplicity, num_atoms, 3).
        frames_idx_true: DTensor of shape (batch, num_tokens, 3).
        feats: DTensor feature dict with required keys (`asym_id`, `atom_to_token`,
            `atom_pad_mask`, `atom_resolved_mask`, `mol_type`, `token_pad_mask`).
        multiplicity: Number of copies per sample in the batch dimension.
        resolved_mask: Optional DTensor (batch, num_atoms) to prefer resolved atoms.
        inference: If True, uses pad mask instead of resolved mask.

    Returns:
        DTensor frame indices and DTensor collinearity mask.
    """
    feats_keys = {
        "asym_id",
        "atom_to_token",
        "atom_pad_mask",
        "atom_resolved_mask",
        "mol_type",
        "token_pad_mask",
    }
    if any(k not in feats for k in feats_keys):
        raise ValueError(f"feats must contain the following keys: {feats_keys}, got {feats.keys()}")

    if frames_idx_true.ndim == 4:
        raise ValueError(
            f"frames_idx_true has unsqueezed ensemble dim (ndim=4, shape={frames_idx_true.shape}). "
            "Only E=1 is supported; squeeze the ensemble dim before calling compute_frame_pred."
        )

    # Check device mesh, placements, and shapes
    device_mesh = pred_atom_coords.device_mesh
    single_repr_placements = (Shard(0), Shard(1), Replicate())
    replicate_placements = (Shard(0), Replicate(), Replicate())

    global_batch_size, num_atoms = feats["atom_pad_mask"].shape
    _, num_tokens = feats["asym_id"].shape
    assert (
        pred_atom_coords.shape[0] == global_batch_size * multiplicity
    ), f"pred_atom_coords must have shape {global_batch_size * multiplicity}, got {pred_atom_coords.shape[0]}"

    expected_placements = {
        "pred_atom_coords": single_repr_placements,
        "frames_idx_true": single_repr_placements,
        "asym_id": single_repr_placements,
        "atom_to_token": single_repr_placements,
        "atom_pad_mask": single_repr_placements,
        "atom_resolved_mask": single_repr_placements,
        "mol_type": single_repr_placements,
        "token_pad_mask": single_repr_placements,  # context parallelism specific
    }
    expected_shape = {
        "pred_atom_coords": (global_batch_size * multiplicity, num_atoms, 3),  # 3D coordinates of the atoms
        "frames_idx_true": (global_batch_size, num_tokens, 3),  # 3 atoms to form a frame per token
        "asym_id": (global_batch_size, num_tokens),  # asym_id of the tokens
        "atom_to_token": (
            global_batch_size,
            num_atoms,
            num_tokens // device_mesh.size(1),
        ),  # mapping from atoms to tokens
        "atom_pad_mask": (global_batch_size, num_atoms),  # padding mask of the atoms
        "atom_resolved_mask": (global_batch_size, num_atoms),  # resolved mask of the atoms
        "mol_type": (global_batch_size, num_tokens),  # mol_type of the tokens
        "token_pad_mask": (global_batch_size, num_tokens),  # padding mask of the tokens (context parallelism specific)
    }

    for k in expected_placements:
        match k:
            case "pred_atom_coords":
                if pred_atom_coords.placements != expected_placements[k]:
                    raise ValueError(
                        f"pred_atom_coords must have placements {expected_placements[k]}, got {pred_atom_coords.placements}"
                    )
                if pred_atom_coords.shape != expected_shape[k]:
                    raise ValueError(
                        f"pred_atom_coords must have shape {expected_shape[k]}, got {pred_atom_coords.shape}"
                    )
            case "frames_idx_true":
                if frames_idx_true.device_mesh != device_mesh:
                    raise ValueError(
                        f"frames_idx_true must be on the same device mesh as pred_atom_coords, got {frames_idx_true.device_mesh} and {device_mesh}"
                    )
                if frames_idx_true.placements != expected_placements[k]:
                    raise ValueError(
                        f"frames_idx_true must have placements {expected_placements[k]}, got {frames_idx_true.placements}"
                    )
                if frames_idx_true.shape != expected_shape[k]:
                    raise ValueError(
                        f"frames_idx_true must have shape {expected_shape[k]}, got {frames_idx_true.shape}"
                    )
            case "resolved_mask":
                if resolved_mask is not None and resolved_mask.device_mesh != device_mesh:
                    raise ValueError(
                        f"resolved_mask must be on the same device mesh as pred_atom_coords, got {resolved_mask.device_mesh} and {device_mesh}"
                    )
                if resolved_mask is not None and resolved_mask.placements != expected_placements[k]:
                    raise ValueError(
                        f"resolved_mask must have placements {expected_placements[k]}, got {resolved_mask.placements}"
                    )
                if resolved_mask is not None and resolved_mask.shape != expected_shape[k]:
                    raise ValueError(f"resolved_mask must have shape {expected_shape[k]}, got {resolved_mask.shape}")
            case _:
                if feats[k].device_mesh != device_mesh:
                    raise ValueError(
                        f"feats[{k}] must be on the same device mesh as pred_atom_coords, got {feats[k].device_mesh} and {device_mesh}"
                    )
                if feats[k].placements != expected_placements[k]:
                    raise ValueError(
                        f"feats[{k}] must have placements {expected_placements[k]}, got {feats[k].placements}"
                    )
                if feats[k].shape != expected_shape[k]:
                    raise ValueError(f"feats[{k}] must have shape {expected_shape[k]}, got {feats[k].shape}")

    # All-gather all inputs/features
    pred_atom_coords_gathered = pred_atom_coords.redistribute(device_mesh, placements=replicate_placements).to_local()
    frames_idx_true_gathered = frames_idx_true.redistribute(device_mesh, placements=replicate_placements).to_local()
    asym_id_gathered = feats["asym_id"].redistribute(device_mesh, placements=replicate_placements).to_local()
    asym_id_atom_gathered = (
        single_repr_token_to_atom(feats["asym_id"].float(), feats["atom_to_token"])
        .redistribute(device_mesh, placements=replicate_placements)
        .to_local()
        .to(torch.int64)
    )
    atom_pad_mask_gathered = (
        feats["atom_pad_mask"].redistribute(device_mesh, placements=replicate_placements).to_local()
    )
    atom_resolved_mask_gathered = (
        feats["atom_resolved_mask"].redistribute(device_mesh, placements=replicate_placements).to_local()
    )
    mol_type_gathered = feats["mol_type"].redistribute(device_mesh, placements=replicate_placements).to_local()
    if resolved_mask is not None:
        resolved_mask_gathered = resolved_mask.redistribute(device_mesh, placements=replicate_placements).to_local()
    else:
        resolved_mask_gathered = None
    token_pad_mask_gathered = (
        feats["token_pad_mask"].redistribute(device_mesh, placements=replicate_placements).to_local()
    )

    feats_gathered = {
        "asym_id": asym_id_gathered,
        "atom_pad_mask": atom_pad_mask_gathered,
        "atom_resolved_mask": atom_resolved_mask_gathered,
        "mol_type": mol_type_gathered,
        "token_pad_mask": token_pad_mask_gathered,
    }

    frames_idx_pred_local, mask_collinear_pred_local = _compute_frame_pred(
        pred_atom_coords_gathered,
        frames_idx_true_gathered,
        feats_gathered,
        asym_id_atom_gathered,
        multiplicity,
        resolved_mask=resolved_mask_gathered,
        inference=inference,
    )

    # Redistribute frames and mask
    shape = torch.Size([global_batch_size, multiplicity, num_tokens, 3])
    cp_submesh = device_mesh["cp_axis_0", "cp_axis_1"]
    frames_idx_pred_cp = distribute_tensor(
        frames_idx_pred_local,
        cp_submesh,
        (Shard(2), Replicate()),
        src_data_rank=0,  # group rank not global rank
    )  # broadcast to rest of cp group to consistency amid potential numerical discrepancies
    frames_idx_pred = DTensor.from_local(
        frames_idx_pred_cp.to_local(),
        device_mesh=device_mesh,
        placements=(Shard(0), Shard(2), Replicate()),
        shape=shape,
        stride=LayoutRightMap(shape).strides,
    )

    shape = torch.Size([global_batch_size, multiplicity, num_tokens])
    mask_collinear_pred_cp = distribute_tensor(
        mask_collinear_pred_local,
        device_mesh["cp_axis_0", "cp_axis_1"],
        (Shard(2), Replicate()),
        src_data_rank=0,  # group rank not global rank
    )  # broadcast to rest of cp group to consistency amid potential numerical discrepancies
    mask_collinear_pred = DTensor.from_local(
        mask_collinear_pred_cp.to_local(),
        device_mesh=device_mesh,
        placements=(Shard(0), Shard(2), Replicate()),
        shape=shape,
        stride=LayoutRightMap(shape).strides,
    )

    return frames_idx_pred, mask_collinear_pred


def lddt_resolved_token(
    pred_atom_coords: DTensor,
    true_atom_coords: DTensor,
    true_coords_resolved_mask: DTensor,
    feats: dict[str, DTensor],
    comm: TransposeComm,
    multiplicity: int = 1,
    cutoff: float = 15.0,
    eps: float = 1e-10,
) -> tuple[DTensor, DTensor]:
    """Compute per-token lDDT scores in distributed setup.

    This function computes the lDDT (local Distance Difference Test) scores for
    each token, which measures the local structural accuracy by comparing pairwise
    distances between predicted and true coordinates.

    The computation uses cdist_lddt with factorized masks and redistribute_transpose
    to handle the distributed pairwise distance computation efficiently.

    Co-sharding Requirements
    ------------------------
    This function assumes specific co-sharding of features along the N_token and N_atom axes:

    1. **Token-Atom Co-sharding**: The N_token axis of token features (mol_type, token_to_rep_atom)
       and the N_atom axis of atom features (pred_atom_coords, true_atom_coords, atom_to_token)
       are co-sharded as diagonal blocks. This means:
       - Atoms of shard i belong ONLY to tokens of shard i
       - The atom-to-token mapping is block-diagonal in the global matrix

    2. **R-set-Token Co-sharding**: Each R-set element corresponds to exactly one token via its
       representative atom. R-set elements are placed in the same shard as their corresponding token:
       - R-set elements of shard i are a SUBSET of tokens in shard i
       - This enables local matmul for r_set_to_rep_atom @ atom_coords

    3. **N_atom_max_per_shard Semantics**: The "atom" axis of token_to_rep_atom, r_set_to_rep_atom,
       and atom_to_token represents the local shard's atom indices (0 to N_atom_max_per_shard-1),
       NOT global atom indices. This is because these tensors are diagonal blocks of the global
       reference versions. See src/boltz/distributed/data/feature/featurizer.py for details.

    These co-sharding properties enable all coordinate projections and mask computations to be
    performed locally without communication (until the final all-reduce for lDDT aggregation).

    Parameters
    ----------
    pred_atom_coords : DTensor
        Predicted atom coordinates with shape [B*mult, N_atom, 3].
        Placements: (Shard(batch_dim), Shard(atom_dim), Replicate())
    true_atom_coords : DTensor
        True atom coordinates with shape [B*mult, N_atom, 3].
        Placements: (Shard(batch_dim), Shard(atom_dim), Replicate())
    true_coords_resolved_mask : DTensor
        Resolved mask for atoms with shape [B*mult, N_atom].
        Placements: (Shard(batch_dim), Shard(atom_dim), Replicate())
    feats : dict[str, DTensor]
        Feature dictionary containing:
        - token_to_rep_atom: One-hot mapping [B, N_token, N_atom_max_per_shard]
        - r_set_to_rep_atom: One-hot mapping [B, N_R, N_atom_max_per_shard]
        - atom_to_token: One-hot mapping [B, N_atom_max_per_shard, N_token]
        - mol_type: Token types [B, N_token]
        All with same device_mesh and placements as pred_atom_coords.
    comm : TransposeComm
        Communication object for redistribute_transpose operations.
    multiplicity : int, optional
        Diffusion batch multiplier, by default 1
    cutoff : float, optional
        Base cutoff distance for lDDT computation, by default 15.0.
        For nucleotide tokens, the cutoff is doubled (cutoff + cutoff * is_nucleotide).
    eps : float, optional
        Small epsilon for numerical stability, by default 1e-10

    Returns
    -------
    target_lddt : DTensor
        Per-token lDDT scores [B*mult, N_token] with same placements as input
    combined_mask : DTensor
        Combined mask (token_resolved_mask * mask_no_match) [B*mult, N_token] with same placements.
        - token_resolved_mask: Whether each token has a resolved representative atom
        - mask_no_match: Whether each token has valid pairs for lDDT computation
    """
    # === Extract device_mesh and placements from input DTensor ===
    device_mesh = pred_atom_coords.device_mesh
    input_placements = pred_atom_coords.placements
    # Validate input placements must be exactly (Shard(0), Shard(1), Replicate())
    expected_placements = (Shard(0), Shard(1), Replicate())
    if input_placements != expected_placements:
        raise ValueError(f"pred_atom_coords placements {input_placements} must be {expected_placements}")

    # Extract features
    token_to_rep_atom = feats["token_to_rep_atom"]  # [B, N_token, N_atom_max_per_shard]
    r_set_to_rep_atom = feats["r_set_to_rep_atom"]  # [B, N_R, N_atom_max_per_shard]
    atom_to_token = feats["atom_to_token"]  # [B, N_atom_max_per_shard, N_token]
    mol_type = feats["mol_type"]  # [B, N_token]

    # === Sanity checks for device_mesh and placements consistency ===
    all_dtensors = [
        ("pred_atom_coords", pred_atom_coords),
        ("true_atom_coords", true_atom_coords),
        ("true_coords_resolved_mask", true_coords_resolved_mask),
        ("token_to_rep_atom", token_to_rep_atom),
        ("r_set_to_rep_atom", r_set_to_rep_atom),
        ("atom_to_token", atom_to_token),
        ("mol_type", mol_type),
    ]
    for name, dtensor in all_dtensors:
        if dtensor.device_mesh != device_mesh:
            raise ValueError(f"{name} has different device_mesh than pred_atom_coords")
        # Check placements match exactly
        if dtensor.placements != expected_placements:
            raise ValueError(f"{name} has placements {dtensor.placements}, expected {expected_placements}")

    # === Extract and validate global shape dimensions ===
    # NOTE: "Global" here refers to DTensor's global semantics (i.e., DTensor.full_tensor() shape),
    # which may differ from the serial equivalent dimensions due to intersperse padding applied
    # during CP data processing and dataloader. For example, N_token_global and N_atom_global
    # include padding to ensure even sharding across CP ranks.
    #
    # Note on DTensor global shape semantics (Shard dim -> global size, Replicate dim -> local size):
    # All have placements (Shard(0), Shard(1), Replicate()) for 3D or (Shard(0), Shard(1), Replicate()) for 2D
    # - token_to_rep_atom: [B, N_token, N_atom_max_per_shard] - dim1 Shard->global, dim2 Replicate->local
    # - r_set_to_rep_atom: [B, N_R, N_atom_max_per_shard] - dim1 Shard->global, dim2 Replicate->local
    # - atom_to_token: [B, N_atom, N_token_max_per_shard] - dim1 Shard->global, dim2 Replicate->local
    # - mol_type: [B, N_token] - dim1 Shard->global
    B_mult_global = pred_atom_coords.shape[0]  # B * multiplicity
    N_atom_global = pred_atom_coords.shape[1]
    B_global = token_to_rep_atom.shape[0]  # B (without multiplicity)
    N_token_global = token_to_rep_atom.shape[1]
    N_R_global = r_set_to_rep_atom.shape[1]
    N_atom_max_per_shard = token_to_rep_atom.shape[2]  # Local atom axis (diagonal block semantics)
    N_token_max_per_shard = atom_to_token.shape[2]  # Local token axis (diagonal block semantics)

    # Validate multiplicity consistency
    if B_mult_global != B_global * multiplicity:
        raise ValueError(
            f"pred_atom_coords batch dim ({B_mult_global}) != B ({B_global}) * multiplicity ({multiplicity})"
        )

    # Validate coordinate shapes
    if true_atom_coords.shape != pred_atom_coords.shape:
        raise ValueError(
            f"true_atom_coords shape {true_atom_coords.shape} != pred_atom_coords shape {pred_atom_coords.shape}"
        )
    if true_coords_resolved_mask.shape != (B_mult_global, N_atom_global):
        raise ValueError(
            f"true_coords_resolved_mask shape {true_coords_resolved_mask.shape} != expected ({B_mult_global}, {N_atom_global})"
        )

    # Validate feature shapes
    # Note: For DTensor, .shape returns global shape. For Shard(dim), global != local.
    # - token_to_rep_atom: dim1 is Shard(1) so shape[1]=N_token_global; dim2 is Replicate so shape[2]=N_atom_max_per_shard
    # - atom_to_token: dim1 is Shard(1) so shape[1]=N_atom_global; dim2 is Replicate so shape[2]=N_token_max_per_shard
    if token_to_rep_atom.shape != (B_global, N_token_global, N_atom_max_per_shard):
        raise ValueError(f"token_to_rep_atom shape {token_to_rep_atom.shape} is invalid")
    if r_set_to_rep_atom.shape != (B_global, N_R_global, N_atom_max_per_shard):
        raise ValueError(f"r_set_to_rep_atom shape {r_set_to_rep_atom.shape} is invalid")
    # atom_to_token: dim1 (atom) is Shard(1)->global, dim2 (token) is Replicate->local
    if atom_to_token.shape != (B_global, N_atom_global, N_token_max_per_shard):
        raise ValueError(f"atom_to_token shape {atom_to_token.shape} is invalid")
    if mol_type.shape != (B_global, N_token_global):
        raise ValueError(f"mol_type shape {mol_type.shape} != expected ({B_global}, {N_token_global})")

    # === Forward-only computation (no gradient support) ===
    # The lDDT metric involves step functions (thresholding distance differences at 0.5, 1, 2, 4 Å),
    # which are not differentiable. Neither the original lddt_dist nor cdist_lddt is mathematically
    # defined with gradients. Therefore, this function does not support backward pass by definition.
    with torch.no_grad():
        # === Get local tensors and consolidate dtype casting ===
        # Promote to at least float32 to match serial get_target_lddt which uses
        # .float() inside autocast(enabled=False).  promote_types preserves
        # float64 for test paths while ensuring BF16 inputs compute in float32.
        pred_coords_local = pred_atom_coords.to_local()  # [local_B*mult, local_N_atom, 3]
        true_coords_local = true_atom_coords.to_local()  # [local_B*mult, local_N_atom, 3]
        coord_dtype = torch.promote_types(pred_coords_local.dtype, torch.float32)
        pred_coords_local = pred_coords_local.to(dtype=coord_dtype)
        true_coords_local = true_coords_local.to(dtype=coord_dtype)

        # Cast all feature tensors to coord_dtype
        token_to_rep_local = token_to_rep_atom.to_local().to(
            dtype=coord_dtype
        )  # [local_B, local_N_token, local_N_atom]
        r_set_to_rep_local = r_set_to_rep_atom.to_local().to(dtype=coord_dtype)  # [local_B, local_N_R, local_N_atom]
        atom_to_token_local = atom_to_token.to_local().to(dtype=coord_dtype)  # [local_B, local_N_atom, local_N_token]
        mol_type_local = mol_type.to_local()  # [local_B, local_N_token]
        resolved_mask_local = true_coords_resolved_mask.to_local().to(dtype=coord_dtype)  # [local_B*mult, local_N_atom]

        # Get local batch size for einsum reshaping
        local_B = token_to_rep_local.shape[0]

        # === Project atom coords to token space (row) using einsum with multiplicity broadcasting ===
        # This avoids repeat_interleave memory overhead by using einsum's implicit broadcasting.
        # Co-sharding assumption: token_to_rep_local operates on local atoms only (diagonal block).
        #
        # token_to_rep_local: [local_B, local_N_token, local_N_atom] -> "bta"
        # pred_coords_local reshaped: [local_B, mult, local_N_atom, 3] -> "bmac"
        # Output: [local_B, mult, local_N_token, 3] -> "bmtc", then reshape to [local_B*mult, local_N_token, 3]
        pred_coords_reshaped = pred_coords_local.view(local_B, multiplicity, -1, 3)  # [local_B, mult, local_N_atom, 3]
        true_coords_reshaped = true_coords_local.view(local_B, multiplicity, -1, 3)  # [local_B, mult, local_N_atom, 3]

        pred_token_coords_row_local = torch.einsum("bta,bmac->bmtc", token_to_rep_local, pred_coords_reshaped).reshape(
            -1, token_to_rep_local.shape[1], 3
        )  # [local_B*mult, local_N_token, 3]
        true_token_coords_row_local = torch.einsum("bta,bmac->bmtc", token_to_rep_local, true_coords_reshaped).reshape(
            -1, token_to_rep_local.shape[1], 3
        )  # [local_B*mult, local_N_token, 3]

        # === Project atom coords to R-set space (col) using einsum with multiplicity broadcasting ===
        # Co-sharding assumption: r_set_to_rep_local operates on local atoms only (diagonal block).
        #
        # r_set_to_rep_local: [local_B, local_N_R, local_N_atom] -> "bra"
        # coords reshaped: [local_B, mult, local_N_atom, 3] -> "bmac"
        # Output: [local_B, mult, local_N_R, 3] -> "bmrc"
        pred_R_coords_col_local = torch.einsum("bra,bmac->bmrc", r_set_to_rep_local, pred_coords_reshaped).reshape(
            -1, r_set_to_rep_local.shape[1], 3
        )  # [local_B*mult, local_N_R, 3]
        true_R_coords_col_local = torch.einsum("bra,bmac->bmrc", r_set_to_rep_local, true_coords_reshaped).reshape(
            -1, r_set_to_rep_local.shape[1], 3
        )  # [local_B*mult, local_N_R, 3]

        # === Compute factorized masks (row and col) using einsum with multiplicity broadcasting ===
        # Masks can vary along the multiplicity axis, so we use einsum to broadcast properly.
        # Co-sharding assumption: the mapping tensors operate on local atoms only.
        #
        # resolved_mask_local: [local_B*mult, local_N_atom] -> reshaped to [local_B, mult, local_N_atom] -> "bma"
        # token_to_rep_local: [local_B, local_N_token, local_N_atom] -> "bta"
        # mask_row: [local_B, mult, local_N_token] -> "bmt", then reshape to [local_B*mult, local_N_token]
        resolved_mask_reshaped = resolved_mask_local.view(local_B, multiplicity, -1)  # [local_B, mult, local_N_atom]

        mask_row_local = torch.einsum("bta,bma->bmt", token_to_rep_local, resolved_mask_reshaped).reshape(
            -1, token_to_rep_local.shape[1]
        )  # [local_B*mult, local_N_token]
        mask_col_local = torch.einsum("bra,bma->bmr", r_set_to_rep_local, resolved_mask_reshaped).reshape(
            -1, r_set_to_rep_local.shape[1]
        )  # [local_B*mult, local_N_R]

        # === Compute cutoff_col based on nucleotide type ===
        # is_nucleotide_token: [local_B, local_N_token]
        # Use atom_to_token_local.dtype for consistency with subsequent bmm operations
        is_nucleotide_token_local = (mol_type_local == const.chain_type_ids["DNA"]).to(
            dtype=atom_to_token_local.dtype
        ) + (mol_type_local == const.chain_type_ids["RNA"]).to(
            dtype=atom_to_token_local.dtype
        )  # [local_B, local_N_token]

        # is_nucleotide_R_element = r_set_to_rep_atom @ (atom_to_token @ is_nucleotide_token)
        # Co-sharding assumption: atom_to_token operates on local atoms/tokens only (diagonal block).
        is_nucleotide_atom_local = torch.bmm(atom_to_token_local, is_nucleotide_token_local.unsqueeze(-1)).squeeze(
            -1
        )  # [local_B, local_N_atom]
        is_nucleotide_R_element_local = torch.bmm(r_set_to_rep_local, is_nucleotide_atom_local.unsqueeze(-1)).squeeze(
            -1
        )  # [local_B, local_N_R]

        # cutoff_col = cutoff + cutoff * is_nucleotide_R_element
        cutoff_col_local = cutoff + cutoff * is_nucleotide_R_element_local  # [local_B, local_N_R]

        # === Get rep_atom indices for diagonal masking (local indices within shard) ===
        # These indices use N_atom_max_per_shard semantics (local atom indices 0 to N_atom_max_per_shard-1).
        # Due to co-sharding, each token/R-element's representative atom is guaranteed to be in the same shard,
        # so argmax on the local diagonal block yields valid local indices for diagonal masking.
        rep_atom_token_local = token_to_rep_local.argmax(dim=-1)  # [local_B, local_N_token]
        rep_atom_r_set_local = r_set_to_rep_local.argmax(dim=-1)  # [local_B, local_N_R]

        # === Derive target placements from input (avoid hardcoding) ===
        # input_placements is e.g. (Shard(0), Shard(1), Replicate()) for the 3D device mesh.
        # Note: DTensor placements tuple length matches mesh dimensions, not tensor dimensions.
        # Target placements after redistribute_transpose: swap Shard(1) <-> Replicate()
        # From (Shard(0), Shard(1), Replicate()) to (Shard(0), Replicate(), Shard(1))
        target_placements = (input_placements[0], input_placements[2], input_placements[1])

        # === Create DTensors for column tensors that need transpose ===
        # Compute shapes and strides for contiguous tensors using LayoutRightMap
        coords_3d_shape = (B_mult_global, N_R_global, 3)
        coords_3d_stride = LayoutRightMap(coords_3d_shape).strides

        pred_R_coords_col_dtensor = DTensor.from_local(
            pred_R_coords_col_local,
            device_mesh,
            input_placements,
            shape=torch.Size(coords_3d_shape),
            stride=coords_3d_stride,
        )
        true_R_coords_col_dtensor = DTensor.from_local(
            true_R_coords_col_local,
            device_mesh,
            input_placements,
            shape=torch.Size(coords_3d_shape),
            stride=coords_3d_stride,
        )

        # Create DTensors for mask_col: [B*mult, N_R] - masks have multiplicity
        mask_col_2d_shape = (B_mult_global, N_R_global)
        mask_col_2d_stride = LayoutRightMap(mask_col_2d_shape).strides

        mask_col_dtensor = DTensor.from_local(
            mask_col_local,
            device_mesh,
            input_placements,
            shape=torch.Size(mask_col_2d_shape),
            stride=mask_col_2d_stride,
        )

        # Create DTensors for cutoff_col and rep_atom_r_set: [B, N_R] - no multiplicity
        feat_2d_shape = (B_global, N_R_global)
        feat_2d_stride = LayoutRightMap(feat_2d_shape).strides

        cutoff_col_dtensor = DTensor.from_local(
            cutoff_col_local,
            device_mesh,
            input_placements,
            shape=torch.Size(feat_2d_shape),
            stride=feat_2d_stride,
        )
        rep_atom_r_set_dtensor = DTensor.from_local(
            rep_atom_r_set_local,
            device_mesh,
            input_placements,
            shape=torch.Size(feat_2d_shape),
            stride=feat_2d_stride,
        )

        # === redistribute_transpose for column tensors ===
        # Transform placements from (S(0), S(1), R) to (S(0), R, S(1)) via all-to-all communication.
        # This distributes the N_R axis across the cp_axis_1 dimension of the device mesh.
        pred_R_coords_col_t = redistribute_transpose(
            pred_R_coords_col_dtensor, comm, target_placements, dim0=None, dim1=None
        )
        true_R_coords_col_t = redistribute_transpose(
            true_R_coords_col_dtensor, comm, target_placements, dim0=None, dim1=None
        )
        mask_col_t = redistribute_transpose(mask_col_dtensor, comm, target_placements, dim0=None, dim1=None)
        cutoff_col_t = redistribute_transpose(cutoff_col_dtensor, comm, target_placements, dim0=None, dim1=None)
        rep_atom_r_set_t = redistribute_transpose(rep_atom_r_set_dtensor, comm, target_placements, dim0=None, dim1=None)

        # === Factorized Pair-Mask Algorithm with cdist_lddt ===
        #
        # The lDDT computation requires pairwise distance comparisons between all (token, R-element) pairs.
        # Instead of materializing the full [N_token, N_R] pair_mask, we use factorized masks:
        #   pair_mask[i,j] = mask_row[i] * mask_col[j]
        #
        # This factorization is valid because:
        # - mask_row[i] = 1 iff token i has a resolved representative atom
        # - mask_col[j] = 1 iff R-element j has a resolved representative atom
        # - A pair (i,j) is valid iff BOTH atoms are resolved
        #
        # For diagonal masking (excluding self-pairs where a token's rep_atom equals an R-element's rep_atom):
        # - atom_indices_row = rep_atom_token: local atom index of each token's representative atom
        # - atom_indices_col = rep_atom_r_set: local atom index of each R-element's representative atom
        #
        # Why local indices work for diagonal masking despite N_atom_max_per_shard semantics:
        # Due to co-sharding, diagonal device_mesh ranks (where cp_axis_0 == cp_axis_1) have:
        # - Row tokens from shard i with local atom indices [0, N_atom_max_per_shard)
        # - Column R-elements ALSO from shard i with the SAME local atom index range
        # Thus, when rep_atom_token[t] == rep_atom_r_set[r], it genuinely means the same physical atom,
        # and the diagonal mask correctly excludes self-pairs.
        #
        # Off-diagonal ranks (cp_axis_0 != cp_axis_1) have row/col from different shards, so their
        # local atom indices never match (different index spaces), and we skip diagonal masking entirely.

        # Determine if this rank is on diagonal (for do_mask_diagonal)
        # Convert to native Python bool for Triton kernel compatibility
        is_diagonal_rank = bool(comm.is_self_comm)

        out_num_local, out_denom_local, mask_no_match_local = cdist_lddt(
            pred_coords_row=pred_token_coords_row_local,  # [local_B*mult, local_N_token, 3]
            pred_coords_col=pred_R_coords_col_t.to_local(),  # [local_B*mult, local_N_R_t, 3]
            true_coords_row=true_token_coords_row_local,  # [local_B*mult, local_N_token, 3]
            true_coords_col=true_R_coords_col_t.to_local(),  # [local_B*mult, local_N_R_t, 3]
            mask_row=mask_row_local,  # [local_B*mult, local_N_token] - factorized row mask
            mask_col=mask_col_t.to_local(),  # [local_B*mult, local_N_R_t] - factorized col mask (transposed)
            multiplicity=multiplicity,
            atom_indices_row=rep_atom_token_local if is_diagonal_rank else None,  # Local indices for diagonal masking
            atom_indices_col=rep_atom_r_set_t.to_local() if is_diagonal_rank else None,  # Transposed local indices
            cutoff_col=cutoff_col_t.to_local(),  # [local_B, local_N_R_t] - per-column cutoff (transposed)
            do_mask_diagonal=is_diagonal_rank,  # Only diagonal ranks need self-pair exclusion
            return_unnormalized_score=True,  # Return partial sums for distributed aggregation
            per_atom=True,  # Per-token output for token-level lDDT
        )
        # Output shapes: out_num, out_denom, mask_no_match are [local_B*mult, local_N_token]

        # === All-reduce across N_R axis (cp_axis_1 group) ===
        # Each rank computed partial lDDT contributions from its local N_R shard.
        # Sum across all N_R shards to get the full lDDT numerator and denominator.
        # This transforms partial sums from (S(0), S(1), partial_N_R) to (S(0), S(1), R).
        group_col = device_mesh.get_group(2)  # cp_axis_1
        dist.all_reduce(out_num_local, op=dist.ReduceOp.SUM, group=group_col)
        dist.all_reduce(out_denom_local, op=dist.ReduceOp.SUM, group=group_col)

        # All-reduce mask_no_match with logical OR (any rank having valid pairs means token has matches)
        dist.all_reduce(mask_no_match_local, op=dist.ReduceOp.MAX, group=group_col)

        # === Compute combined_mask = token_resolved_mask * mask_no_match ===
        # mask_row_local is already token_resolved_mask (computed above via einsum)
        # Both masks don't require gradients
        combined_mask_local = mask_row_local * mask_no_match_local

        # === Normalize to get final lDDT scores ===
        # Preserve input coordinate dtype (e.g., float64 for precision)
        norm = 1.0 / (eps + out_denom_local)
        target_lddt_local = norm * (eps + out_num_local)

        # === Wrap outputs as DTensors ===
        output_shape = (B_mult_global, N_token_global)
        output_stride = LayoutRightMap(output_shape).strides

        target_lddt = DTensor.from_local(
            target_lddt_local,
            device_mesh,
            input_placements,
            shape=torch.Size(output_shape),
            stride=output_stride,
        )
        combined_mask = DTensor.from_local(
            combined_mask_local,
            device_mesh,
            input_placements,
            shape=torch.Size(output_shape),
            stride=output_stride,
        )

        return target_lddt, combined_mask


class _PLDDTLossImpl(torch.autograd.Function):
    """Fused pLDDT loss computation with gradient flow to pred_lddt.

    This fuses the entire pLDDT loss computation into a single autograd Function:
    1. Cross-entropy errors: errors = -sum(one_hot * log_softmax(pred_lddt), dim=-1)
    2. Masked errors: errors * combined_mask
    3. Sum over token dim (all_reduce over CP axis)
    4. Normalize: numerator / clamp(denominator, min=eps)
    5. Sum over batch dim (all_reduce over DP axis)
    6. Mean: loss_sum / batch_size

    WHY THIS WORKS:
    ---------------
    dist.all_reduce is an in-place op that is invisible to PyTorch autograd.
    Autograd records the computation graph as if all_reduce never happened,
    but the actual tensor values ARE the all_reduced values.

    This is correct because all_reduce(SUM) has IDENTITY GRADIENT:
        Forward:  y = sum_over_ranks(x_i), all ranks get the same y
        Backward: ∂L/∂x_i = ∂L/∂y * ∂y/∂x_i = ∂L/∂y * 1 = ∂L/∂y

    So autograd "accidentally" computes the correct gradient by ignoring
    all_reduce, since passing the gradient through unchanged is exactly
    what all_reduce(SUM) backward should do.

    IMPORTANT: This ONLY works for ReduceOp.SUM. Other ops (MEAN, MAX, etc.)
    have non-identity gradients and would produce wrong results. However,
    this is guaranteed by the pLDDT loss semantics: we're computing
        loss = mean_b(sum_t(errors * mask) / sum_t(mask))
    which requires SUM reduction to accumulate partial sums across ranks.

    Gradient flows back to pred_lddt through the log_softmax operation.
    target_lddt and combined_mask are non-differentiable.

    See Also
    --------
    plddt_loss : The DTensor wrapper API that calls this function.
    lddt_resolved_token : Computes target_lddt and combined_mask.
    """

    @staticmethod
    @torch.amp.custom_fwd(device_type="cuda")
    def forward(
        ctx: FunctionCtx,
        pred_lddt: DTensor,
        target_lddt: DTensor,
        combined_mask: DTensor,
    ) -> DTensor:
        """Forward pass for fused pLDDT loss computation.

        Parameters
        ----------
        ctx : FunctionCtx
            The autograd context object for saving tensors for backward.
        pred_lddt : DTensor
            Predicted lDDT logits with shape (B * multiplicity, N_token, num_bins).
            Placements: (Shard(0), Shard(1), Replicate()).
            This is the only tensor that requires gradients.
        target_lddt : DTensor
            Target lDDT scores with shape (B * multiplicity, N_token), values in [0, 1].
            Placements: (Shard(0), Shard(1), Replicate()).
            No gradients (computed from lddt_resolved_token with step functions).
        combined_mask : DTensor
            Combined mask (token_resolved_mask * mask_no_match) with shape (B * multiplicity, N_token).
            Placements: (Shard(0), Shard(1), Replicate()).
            No gradients.

        Returns
        -------
        DTensor
            Scalar loss with placements (Replicate(), Replicate(), Replicate()).
        """
        # === Validate input dimensions ===
        if pred_lddt.ndim != 3:
            raise ValueError(f"pred_lddt must be 3D (B*mult, N_token, num_bins), got {pred_lddt.ndim}D")
        if target_lddt.ndim != 2:
            raise ValueError(f"target_lddt must be 2D (B*mult, N_token), got {target_lddt.ndim}D")
        if combined_mask.ndim != 2:
            raise ValueError(f"combined_mask must be 2D (B*mult, N_token), got {combined_mask.ndim}D")

        # === Validate device_mesh consistency ===
        device_mesh = pred_lddt.device_mesh
        for name, dtensor in [
            ("target_lddt", target_lddt),
            ("combined_mask", combined_mask),
        ]:
            if dtensor.device_mesh != device_mesh:
                raise ValueError(f"{name} has different device_mesh than pred_lddt")

        # === Validate placements ===
        # All inputs should have placements (S(0), S(1), R) for 3D or (S(0), S(1), R) conceptually for 2D
        pred_placements = pred_lddt.placements
        expected_pred_placements = (Shard(0), Shard(1), Replicate())
        if pred_placements != expected_pred_placements:
            raise ValueError(f"pred_lddt placements {pred_placements} must be {expected_pred_placements}")

        # For 2D tensors, placements should match first two dimensions
        expected_2d_placements = (Shard(0), Shard(1), Replicate())
        for name, dtensor in [
            ("target_lddt", target_lddt),
            ("combined_mask", combined_mask),
        ]:
            if dtensor.placements != expected_2d_placements:
                raise ValueError(f"{name} placements {dtensor.placements} must be {expected_2d_placements}")

        # === Validate shape consistency ===
        B_mult_global = pred_lddt.shape[0]
        N_token_global = pred_lddt.shape[1]

        if target_lddt.shape != (B_mult_global, N_token_global):
            raise ValueError(f"target_lddt shape {target_lddt.shape} != expected ({B_mult_global}, {N_token_global})")
        if combined_mask.shape != (B_mult_global, N_token_global):
            raise ValueError(
                f"combined_mask shape {combined_mask.shape} != expected ({B_mult_global}, {N_token_global})"
            )

        # === Get process groups for all_reduce ===
        group_cp = device_mesh.get_group(1)  # CP axis (token dimension)
        group_dp = device_mesh.get_group(0)  # DP axis (batch dimension)

        # === Get local tensors ===
        pred_lddt_local = pred_lddt.to_local().detach().requires_grad_(pred_lddt.requires_grad)
        target_lddt_local = target_lddt.to_local().detach()  # No gradient
        combined_mask_local = combined_mask.to_local().detach()  # No gradient

        # Compute bin indices from target_lddt (no gradient flow through this)
        num_bins = pred_lddt_local.shape[-1]
        bin_index = torch.floor(target_lddt_local * num_bins).long()
        bin_index = torch.clamp(bin_index, max=(num_bins - 1))

        # One-hot encode target bins (no gradient)
        lddt_one_hot = F.one_hot(bin_index, num_classes=num_bins).to(pred_lddt_local.dtype)

        # === Fused subgraph: errors + mask + sum + normalize + sum + mean ===
        with torch.enable_grad():
            # Compute cross-entropy errors (gradient flows through log_softmax)
            log_probs = F.log_softmax(pred_lddt_local, dim=-1)
            errors_local = -torch.sum(lddt_one_hot * log_probs, dim=-1)  # [local_B*mult, local_N_token]

            # Apply combined_mask
            masked_errors_local = errors_local * combined_mask_local  # [local_B*mult, local_N_token]

            # Sum over token dimension (local sum first)
            numerator_local = masked_errors_local.sum(dim=-1)  # [local_B*mult]
            denominator_local = combined_mask_local.sum(dim=-1)  # [local_B*mult]

            # All-reduce over CP axis (token dimension was sharded)
            # Clone to protect against potential upstream saved_for_backward
            # (though sum() output is not typically saved by its upstream, clone is safe and cheap)
            numerator = numerator_local.clone()
            denominator = denominator_local.clone()
            with torch.no_grad():
                dist.all_reduce(numerator, op=dist.ReduceOp.SUM, group=group_cp)
                dist.all_reduce(denominator, op=dist.ReduceOp.SUM, group=group_cp)

            # Normalize: numerator / clamp(denominator, min=eps)
            eps = 1e-7
            denominator_safe = torch.clamp(denominator, min=eps)
            per_sample_loss_local = numerator / denominator_safe  # [local_B*mult]

            # Sum over batch dimension (local sum first)
            loss_sum_local = per_sample_loss_local.sum()  # scalar

            # All-reduce over DP axis (batch dimension was sharded)
            # Clone needed: sum() may save input for backward
            loss_sum = loss_sum_local.clone()
            with torch.no_grad():
                dist.all_reduce(loss_sum, op=dist.ReduceOp.SUM, group=group_dp)

            # Mean over global batch size
            loss_local = loss_sum / B_mult_global

        # === Save for backward ===
        ctx.save_for_backward(pred_lddt_local, loss_local)
        ctx.device_mesh = device_mesh
        ctx.pred_placements = pred_placements
        ctx.pred_lddt_shape = pred_lddt.shape
        ctx.pred_lddt_stride = pred_lddt.stride()

        # === Wrap output as DTensor ===
        # Output is a scalar, fully replicated across all mesh dimensions
        output_placements = (Replicate(), Replicate(), Replicate())

        loss = DTensor.from_local(
            loss_local.detach(),
            device_mesh=device_mesh,
            placements=output_placements,
            shape=torch.Size(()),
            stride=(),
        )

        return loss

    @staticmethod
    @torch.amp.custom_bwd(device_type="cuda")
    def backward(
        ctx: FunctionCtx,
        grad_loss: DTensor,
    ) -> tuple[DTensor | None, None, None]:
        """Backward pass using subgraph trick via autograd.grad.

        The forward pass fuses the entire pLDDT loss computation into one subgraph.
        Autograd ignores all_reduce (which has identity gradient for SUM),
        so we can backprop through the entire fused computation in one call.

        Parameters
        ----------
        ctx : FunctionCtx
            Context containing saved tensors and metadata.
        grad_loss : DTensor
            Gradient of upstream loss w.r.t. this loss output (scalar).

        Returns
        -------
        tuple
            (grad_pred_lddt, None, None) - only pred_lddt has gradients.
        """
        pred_lddt_local, loss_local = ctx.saved_tensors

        if not pred_lddt_local.requires_grad:
            return None, None, None

        # === Validate grad_loss ===
        if grad_loss.shape != torch.Size(()):
            raise ValueError(f"grad_loss must be scalar, got shape {grad_loss.shape}")

        expected_grad_placements = (Replicate(), Replicate(), Replicate())
        if grad_loss.placements != expected_grad_placements:
            raise ValueError(f"grad_loss placements {grad_loss.placements} must be {expected_grad_placements}")

        if grad_loss.device_mesh != ctx.device_mesh:
            raise ValueError("grad_loss has different device_mesh than forward inputs")

        # Get local gradient (scalar, replicated)
        grad_loss_local = grad_loss.to_local()

        # Backprop through entire fused subgraph: loss -> pred_lddt
        # Autograd treats all_reduce as invisible, which is correct since
        # all_reduce(SUM) has identity gradient (see forward pass comment).
        (grad_pred_lddt_local,) = torch.autograd.grad(
            outputs=[loss_local],
            inputs=[pred_lddt_local],
            grad_outputs=[grad_loss_local],
            retain_graph=False,
        )

        # Wrap as DTensor
        grad_pred_lddt = DTensor.from_local(
            grad_pred_lddt_local,
            device_mesh=ctx.device_mesh,
            placements=ctx.pred_placements,
            shape=ctx.pred_lddt_shape,
            stride=ctx.pred_lddt_stride,
        )

        return grad_pred_lddt, None, None


def plddt_loss(
    pred_lddt: DTensor,
    pred_atom_coords: DTensor,
    true_atom_coords: DTensor,
    true_coords_resolved_mask: DTensor,
    feats: dict[str, DTensor],
    comm: TransposeComm,
    multiplicity: int = 1,
    cutoff: float = 15.0,
) -> DTensor:
    """Compute the pLDDT loss using DTensor.

    This is the DTensor version of boltz.model.loss.confidencev2.plddt_loss.
    It computes the cross-entropy loss between predicted lDDT bins and target
    lDDT scores computed from coordinates.

    The entire computation is fused into _PLDDTLossImpl using native PyTorch ops:
    1. lddt_resolved_token: Compute per-token target lDDT scores (no gradient)
    2. _PLDDTLossImpl fuses:
       - Cross-entropy errors: -sum(one_hot * log_softmax(pred_lddt), dim=-1)
       - Masked errors: errors * combined_mask
       - Sum over token dim + all_reduce(SUM) across CP
       - Normalize: numerator / clamp(denominator, min=eps)
       - Sum over batch dim + all_reduce(SUM) across DP
       - Mean: loss_sum / batch_size

    All all_reduce ops use SUM which has identity gradient, allowing the entire
    computation to be captured in a single autograd subgraph.

    Parameters
    ----------
    pred_lddt : DTensor
        Predicted lDDT logits with shape (B * multiplicity, N_token, num_bins).
        Placements: (Shard(0), Shard(1), Replicate()).
    pred_atom_coords : DTensor
        Predicted atom coordinates with shape (B * multiplicity, N_atom, 3).
        Placements: (Shard(0), Shard(1), Replicate()).
    true_atom_coords : DTensor
        Ground truth atom coordinates with shape (B * multiplicity, N_atom, 3).
        Placements: (Shard(0), Shard(1), Replicate()).
    true_coords_resolved_mask : DTensor
        Mask for resolved coordinates with shape (B * multiplicity, N_atom).
        Placements: (Shard(0), Shard(1), Replicate()).
    feats : dict[str, DTensor]
        Dictionary containing feature tensors:
        - "token_to_rep_atom": [B, N_token, N_atom_max_per_shard]
        - "r_set_to_rep_atom": [B, N_R, N_atom_max_per_shard]
        - "atom_to_token": [B, N_atom, N_token_max_per_shard]
        - "mol_type": [B, N_token]
    comm : TransposeComm
        Communication object for redistribute_transpose operations.
    multiplicity : int, optional
        Diffusion batch multiplier, by default 1.
    cutoff : float, optional
        Base cutoff distance for lDDT computation, by default 15.0.

    Returns
    -------
    DTensor
        Scalar loss with placements (Replicate(), Replicate(), Replicate()).
    """
    # Compute target lDDT and combined_mask (no gradients, uses step functions)
    # combined_mask = token_resolved_mask * mask_no_match (computed inside lddt_resolved_token)
    target_lddt, combined_mask = lddt_resolved_token(
        pred_atom_coords=pred_atom_coords,
        true_atom_coords=true_atom_coords,
        true_coords_resolved_mask=true_coords_resolved_mask,
        feats=feats,
        comm=comm,
        multiplicity=multiplicity,
        cutoff=cutoff,
    )

    # _PLDDTLossImpl fuses the entire loss computation using native PyTorch ops:
    # errors -> mask -> sum(CP) -> normalize -> sum(DP) -> mean
    loss = _PLDDTLossImpl.apply(pred_lddt, target_lddt, combined_mask)

    return loss


class _PDELossImpl(torch.autograd.Function):
    """Shardwise computation of PDE loss with gradient flow to pred_pde.

    This computes the PDE cross-entropy loss per token row, with all-reduce across
    the column dimension (cp_axis_1). The row-wise aggregation is done by sharded_sum
    in the wrapper API pde_loss, which has proper autograd support.

    The computation uses cdist_pde which fuses:
        - Distance computation: true_d = cdist(true_coords_row, true_coords_col)
        - Distance computation: pred_d = cdist(pred_coords_row, pred_coords_col)
        - Target PDE: target_pde = abs(true_d - pred_d)
        - Binning: bin_index = clamp(floor(target_pde * num_bins / max_dist), max=num_bins-1)
        - Cross-entropy: errors = -sum(one_hot(bin_index) * log_softmax(pred_pde), dim=-1)
        - Masked sum along column: out_loss_num = sum(errors * mask, dim=-1)

    Gradient flows back to pred_pde through the log_softmax operation in cdist_pde.
    Coordinates and masks are non-differentiable (no gradient flow).

    See Also
    --------
    pde_loss : The DTensor wrapper API that calls this function and does sharded_sum.
    cdist_pde : The fused Triton kernel for PDE cross-entropy computation.
    """

    @staticmethod
    @torch.amp.custom_fwd(device_type="cuda")
    def forward(
        ctx: FunctionCtx,
        pred_pde: DTensor,
        pred_atom_coords: DTensor,
        true_atom_coords: DTensor,
        true_coords_resolved_mask: DTensor,
        token_to_rep_atom: DTensor,
        comm: TransposeComm,
        multiplicity: int,
        max_dist: float,
    ) -> tuple[DTensor, DTensor]:
        """Forward pass for PDE loss computation with transpose and all-reduce.

        Parameters
        ----------
        ctx : FunctionCtx
            The autograd context object for saving tensors for backward.
        pred_pde : DTensor
            Predicted PDE logits with shape (B * multiplicity, N_token, N_token, num_bins).
            Placements: (Shard(0), Shard(1), Shard(2)).
            This is the only tensor that requires gradients.
        pred_atom_coords : DTensor
            Predicted atom coordinates with shape (B * multiplicity, N_atom, 3).
            Placements: (Shard(0), Shard(1), Replicate()).
        true_atom_coords : DTensor
            Ground truth atom coordinates with shape (B * multiplicity, N_atom, 3).
            Placements: (Shard(0), Shard(1), Replicate()).
        true_coords_resolved_mask : DTensor
            Resolved mask with shape (B * multiplicity, N_atom).
            Placements: (Shard(0), Shard(1), Replicate()).
        token_to_rep_atom : DTensor
            Token to representative atom mapping with shape (B, N_token, N_atom).
            Placements: (Shard(0), Shard(1), Replicate()).
        comm : TransposeComm
            Communication object for redistribute_transpose operations.
        multiplicity : int
            Diffusion batch multiplier.
        max_dist : float
            Maximum distance for binning.

        Returns
        -------
        out_loss_num : DTensor
            Partial sum of cross-entropy loss per row, shape [B*mult, N_token_row].
            Placements: (Shard(0), Shard(1), Replicate()).
        out_mask_denom : DTensor
            Partial sum of mask per row, shape [B*mult, N_token_row].
            Placements: (Shard(0), Shard(1), Replicate()).
        """
        # === Validate input dimensions ===
        if pred_pde.ndim != 4:
            raise ValueError(f"pred_pde must be 4D (B*mult, N_token, N_token, num_bins), got {pred_pde.ndim}D")
        if pred_atom_coords.ndim != 3:
            raise ValueError(f"pred_atom_coords must be 3D (B*mult, N_atom, 3), got {pred_atom_coords.ndim}D")
        if true_atom_coords.ndim != 3:
            raise ValueError(f"true_atom_coords must be 3D (B*mult, N_atom, 3), got {true_atom_coords.ndim}D")
        if true_coords_resolved_mask.ndim != 2:
            raise ValueError(
                f"true_coords_resolved_mask must be 2D (B*mult, N_atom), got {true_coords_resolved_mask.ndim}D"
            )
        if token_to_rep_atom.ndim != 3:
            raise ValueError(
                f"token_to_rep_atom must be 3D (B, N_token, N_atom_max_per_shard), got {token_to_rep_atom.ndim}D"
            )

        # === Validate device_mesh consistency ===
        device_mesh = pred_pde.device_mesh
        for name, dtensor in [
            ("pred_atom_coords", pred_atom_coords),
            ("true_atom_coords", true_atom_coords),
            ("true_coords_resolved_mask", true_coords_resolved_mask),
            ("token_to_rep_atom", token_to_rep_atom),
        ]:
            if dtensor.device_mesh != device_mesh:
                raise ValueError(f"{name} has different device_mesh than pred_pde")

        # === Validate placements ===
        # pred_pde is pair representation: (S(0), S(1), S(2)) - sharded on batch and both token axes
        pred_pde_placements = pred_pde.placements
        expected_pred_pde_placements = (Shard(0), Shard(1), Shard(2))
        if pred_pde_placements != expected_pred_pde_placements:
            raise ValueError(f"pred_pde placements {pred_pde_placements} must be {expected_pred_pde_placements}")

        # Other inputs are single representation: (S(0), S(1), R) - sharded on batch and one axis
        input_placements = pred_atom_coords.placements
        expected_input_placements = (Shard(0), Shard(1), Replicate())
        for name, dtensor in [
            ("pred_atom_coords", pred_atom_coords),
            ("true_atom_coords", true_atom_coords),
            ("true_coords_resolved_mask", true_coords_resolved_mask),
            ("token_to_rep_atom", token_to_rep_atom),
        ]:
            if dtensor.placements != expected_input_placements:
                raise ValueError(f"{name} placements {dtensor.placements} must be {expected_input_placements}")

        # === Validate shape consistency ===
        B_mult_global = pred_pde.shape[0]
        N_token_global = pred_pde.shape[1]
        N_atom_global = pred_atom_coords.shape[1]
        B_global = token_to_rep_atom.shape[0]

        # pred_pde must be square in token dimensions
        if pred_pde.shape[1] != pred_pde.shape[2]:
            raise ValueError(f"pred_pde token dimensions must be equal, got shape {pred_pde.shape}")

        # Validate multiplicity consistency
        if B_mult_global != B_global * multiplicity:
            raise ValueError(f"pred_pde batch dim ({B_mult_global}) != B ({B_global}) * multiplicity ({multiplicity})")

        # Validate coordinate shapes match
        if true_atom_coords.shape != pred_atom_coords.shape:
            raise ValueError(
                f"true_atom_coords shape {true_atom_coords.shape} != pred_atom_coords shape {pred_atom_coords.shape}"
            )

        # Validate resolved_mask shape
        if true_coords_resolved_mask.shape != (B_mult_global, N_atom_global):
            raise ValueError(
                f"true_coords_resolved_mask shape {true_coords_resolved_mask.shape} != expected ({B_mult_global}, {N_atom_global})"
            )

        # Validate token_to_rep_atom token dimension matches pred_pde
        if token_to_rep_atom.shape[1] != N_token_global:
            raise ValueError(
                f"token_to_rep_atom N_token ({token_to_rep_atom.shape[1]}) != pred_pde N_token ({N_token_global})"
            )

        # Get local tensors
        pred_pde_local = pred_pde.to_local().detach().requires_grad_(pred_pde.requires_grad)
        pred_coords_local = pred_atom_coords.to_local().detach()  # [local_B*mult, local_N_atom, 3]
        true_coords_local = true_atom_coords.to_local().detach()  # [local_B*mult, local_N_atom, 3]
        token_to_rep_local = (
            token_to_rep_atom.to_local().detach()
        )  # [local_B, local_N_token, local_N_atom_max_per_shard]
        resolved_mask_local = true_coords_resolved_mask.to_local().detach()  # [local_B*mult, local_N_atom]

        # Validate local atom dimension consistency (co-sharding requirement)
        # token_to_rep_local.shape[2] is N_atom_max_per_shard (Replicate dim -> local)
        # pred_coords_local.shape[1] is local N_atom (Shard(1) dim -> local)
        # These must match for the einsum to work correctly
        local_N_atom = pred_coords_local.shape[1]
        N_atom_max_per_shard = token_to_rep_local.shape[2]
        if N_atom_max_per_shard != local_N_atom:
            raise ValueError(
                f"Co-sharding violation: token_to_rep_atom local atom dim ({N_atom_max_per_shard}) "
                f"!= pred_atom_coords local atom dim ({local_N_atom}). "
                "These must match due to diagonal block co-sharding semantics."
            )

        # Promote to at least float32 to match serial get_target_pde which uses
        # .float() inside autocast(enabled=False).
        coord_dtype = torch.promote_types(pred_coords_local.dtype, torch.float32)
        pred_coords_local = pred_coords_local.to(dtype=coord_dtype)
        true_coords_local = true_coords_local.to(dtype=coord_dtype)
        token_to_rep_local = token_to_rep_local.to(dtype=coord_dtype)

        # Get local batch size for einsum reshaping
        local_B = token_to_rep_local.shape[0]
        local_N_token = token_to_rep_local.shape[1]

        # === Project atom coords to token space using einsum with multiplicity broadcasting ===
        # Co-sharding assumption: token_to_rep_local operates on local atoms only (diagonal block).
        pred_coords_reshaped = pred_coords_local.view(local_B, multiplicity, -1, 3)
        true_coords_reshaped = true_coords_local.view(local_B, multiplicity, -1, 3)

        # token_to_rep_local: [local_B, local_N_token, local_N_atom] -> "bta"
        # coords reshaped: [local_B, mult, local_N_atom, 3] -> "bmac"
        # Output: [local_B, mult, local_N_token, 3] -> "bmtc"
        pred_token_coords_row_local = torch.einsum("bta,bmac->bmtc", token_to_rep_local, pred_coords_reshaped).reshape(
            -1, local_N_token, 3
        )  # [local_B*mult, local_N_token, 3]
        true_token_coords_row_local = torch.einsum("bta,bmac->bmtc", token_to_rep_local, true_coords_reshaped).reshape(
            -1, local_N_token, 3
        )  # [local_B*mult, local_N_token, 3]

        # === Compute factorized mask using einsum with multiplicity broadcasting ===
        resolved_mask_reshaped = resolved_mask_local.view(local_B, multiplicity, -1).to(
            dtype=coord_dtype
        )  # [local_B, mult, local_N_atom]
        mask_row_local = torch.einsum("bta,bma->bmt", token_to_rep_local, resolved_mask_reshaped).reshape(
            -1, local_N_token
        )  # [local_B*mult, local_N_token]

        # === Get global shapes for DTensor metadata ===
        B_mult_global = pred_pde.shape[0]
        N_token_global = pred_pde.shape[1]  # Row dimension
        num_bins = pred_pde.shape[-1]

        # === Derive target placements for transpose ===
        # From (S(0), S(1), R) to (S(0), R, S(1))
        target_placements = (input_placements[0], input_placements[2], input_placements[1])

        # === Create DTensors for column tensors that need transpose ===
        coords_3d_shape = (B_mult_global, N_token_global, 3)
        coords_3d_stride = LayoutRightMap(coords_3d_shape).strides

        pred_token_coords_col_dtensor = DTensor.from_local(
            pred_token_coords_row_local.clone(),  # Clone since we need separate row/col
            device_mesh,
            input_placements,
            shape=torch.Size(coords_3d_shape),
            stride=coords_3d_stride,
        )
        true_token_coords_col_dtensor = DTensor.from_local(
            true_token_coords_row_local.clone(),
            device_mesh,
            input_placements,
            shape=torch.Size(coords_3d_shape),
            stride=coords_3d_stride,
        )

        # Create DTensor for mask_col: [B*mult, N_token]
        mask_2d_shape = (B_mult_global, N_token_global)
        mask_2d_stride = LayoutRightMap(mask_2d_shape).strides

        mask_col_dtensor = DTensor.from_local(
            mask_row_local.clone(),
            device_mesh,
            input_placements,
            shape=torch.Size(mask_2d_shape),
            stride=mask_2d_stride,
        )

        # === redistribute_transpose for column tensors ===
        pred_token_coords_col_t = redistribute_transpose(
            pred_token_coords_col_dtensor, comm, target_placements, dim0=None, dim1=None
        )
        true_token_coords_col_t = redistribute_transpose(
            true_token_coords_col_dtensor, comm, target_placements, dim0=None, dim1=None
        )
        mask_col_t = redistribute_transpose(mask_col_dtensor, comm, target_placements, dim0=None, dim1=None)

        # === Fused subgraph: cdist_pde + all_reduce (CP) + normalize + all_reduce (DP) + mean ===
        #
        # WHY THIS WORKS:
        # ---------------
        # dist.all_reduce is an in-place op that is invisible to PyTorch autograd.
        # Autograd records the computation graph as if all_reduce never happened,
        # but the actual tensor values ARE the all_reduced values.
        #
        # This is correct because all_reduce(SUM) has IDENTITY GRADIENT:
        #   Forward:  y = sum_over_ranks(x_i), all ranks get the same y
        #   Backward: ∂L/∂x_i = ∂L/∂y * ∂y/∂x_i = ∂L/∂y * 1 = ∂L/∂y
        #
        # So autograd "accidentally" computes the correct gradient by ignoring
        # all_reduce, since passing the gradient through unchanged is exactly
        # what all_reduce(SUM) backward should do.
        #
        # IMPORTANT: This ONLY works for ReduceOp.SUM. Other ops (MEAN, MAX, etc.)
        # have non-identity gradients and would produce wrong results. However,
        # this is guaranteed by the PDE loss semantics: we're computing
        #   loss = mean_b(sum_{i,j}(errors * mask) / sum_{i,j}(mask))
        # which requires SUM reduction to accumulate partial sums across ranks.
        #
        # Get process group for DP all_reduce (batch dimension)
        group_dp = device_mesh.get_group(0)

        with torch.enable_grad():
            # Kernel returns fully summed outputs [B_mul] (sum over both row and col axes)
            out_loss_num_local, out_mask_denom_local = cdist_pde(
                pred_pde=pred_pde_local,
                true_coords_row=true_token_coords_row_local,
                true_coords_col=true_token_coords_col_t.to_local(),
                pred_coords_row=pred_token_coords_row_local,
                pred_coords_col=pred_token_coords_col_t.to_local(),
                mask_row=mask_row_local,
                mask_col=mask_col_t.to_local(),
                multiplicity=multiplicity,
                num_bins=num_bins,
                max_dist=max_dist,
            )

            # All-reduce across full CP group for numerator and denominator.
            # Shape: [B_mul_local] where B_mul_local = (B * multiplicity) / dp_size.
            # Typically small (e.g., 1-4 elements), so clone cost is negligible.
            #
            # Clone prevents in-place all_reduce from modifying tensors that upstream
            # ops (cdist_pde) might have saved for backward. This is not strictly
            # necessary here because we know cdist_pde only saves its inputs (pred_pde,
            # coords, masks), not its outputs (out_loss_num_local, out_mask_denom_local).
            # We keep the clone for safety and clarity since the cost is negligible.
            #
            # torch.no_grad() is required to prevent autograd from recording all_reduce
            # in the computation graph. Without it, autograd.grad() encounters the
            # all_reduce node during backward and emits a warning because all_reduce
            # has no registered autograd kernel. This is safe because all_reduce(SUM)
            # has identity gradient (grad just passes through unchanged).
            numerator = out_loss_num_local.clone()
            denominator = out_mask_denom_local.clone()
            with torch.no_grad():
                dist.all_reduce(numerator, op=dist.ReduceOp.SUM, group=comm.group)
                dist.all_reduce(denominator, op=dist.ReduceOp.SUM, group=comm.group)

            # Elementwise ops: clamp and divide to get per-sample loss
            eps = 1e-7
            denominator_safe = torch.clamp(denominator, min=eps)
            per_sample_loss_local = numerator / denominator_safe  # [B_mul_local]

            # Sum over local batch samples
            loss_sum_local = per_sample_loss_local.sum()  # scalar

            # All-reduce across DP group to get global sum.
            # Shape: scalar (0-dim tensor), so clone cost is negligible.
            #
            # Clone needed: the .sum() op is a standard PyTorch operation that may
            # internally save tensors for backward. Clone ensures in-place all_reduce
            # won't corrupt any saved state.
            #
            # torch.no_grad() prevents autograd from recording all_reduce in the graph.
            # Same reasoning as above: all_reduce(SUM) has identity gradient.
            loss_sum = loss_sum_local.clone()
            with torch.no_grad():
                dist.all_reduce(loss_sum, op=dist.ReduceOp.SUM, group=group_dp)

            # Mean over global batch size
            loss_local = loss_sum / B_mult_global

        # === Save for backward ===
        ctx.save_for_backward(pred_pde_local, loss_local)
        ctx.device_mesh = device_mesh
        ctx.pred_pde_placements = pred_pde_placements
        ctx.pred_pde_shape = pred_pde.shape
        ctx.pred_pde_stride = pred_pde.stride()

        # === Wrap output as DTensor ===
        # Output is a scalar, fully replicated across all mesh dimensions
        output_placements = (Replicate(), Replicate(), Replicate())

        loss = DTensor.from_local(
            loss_local.detach(),
            device_mesh=device_mesh,
            placements=output_placements,
            shape=torch.Size(()),
            stride=(),
        )

        return loss

    @staticmethod
    @torch.amp.custom_bwd(device_type="cuda")
    def backward(
        ctx: FunctionCtx,
        grad_loss: DTensor,
    ) -> tuple[DTensor | None, None, None, None, None, None, None, None]:
        """Backward pass using subgraph trick via autograd.grad.

        The forward pass fuses the entire PDE loss computation into one subgraph:
        cdist_pde -> all_reduce(CP) -> normalize -> sum -> all_reduce(DP) -> mean.
        Autograd ignores all_reduce (which has identity gradient for SUM),
        so we can backprop through the entire fused computation in one call.

        Parameters
        ----------
        ctx : FunctionCtx
            Context containing saved tensors and metadata.
        grad_loss : DTensor
            Gradient of upstream loss w.r.t. this loss output (scalar).

        Returns
        -------
        tuple
            (grad_pred_pde, None, None, None, None, None, None, None)
            Only pred_pde has gradients.
        """
        pred_pde_local, loss_local = ctx.saved_tensors

        if not pred_pde_local.requires_grad:
            return None, None, None, None, None, None, None, None

        # === Validate grad_loss ===
        # grad_loss should be a scalar DTensor with fully replicated placements
        if grad_loss.shape != torch.Size(()):
            raise ValueError(f"grad_loss must be scalar, got shape {grad_loss.shape}")

        expected_grad_placements = (Replicate(), Replicate(), Replicate())
        if grad_loss.placements != expected_grad_placements:
            raise ValueError(f"grad_loss placements {grad_loss.placements} must be {expected_grad_placements}")

        if grad_loss.device_mesh != ctx.device_mesh:
            raise ValueError("grad_loss has different device_mesh than forward inputs")

        # Get local gradient (scalar, replicated)
        grad_loss_local = grad_loss.to_local()

        # Backprop through entire fused subgraph: loss -> pred_pde
        # Autograd treats all_reduce as invisible, which is correct since
        # all_reduce(SUM) has identity gradient (see forward pass comment).
        (grad_pred_pde_local,) = torch.autograd.grad(
            outputs=[loss_local],
            inputs=[pred_pde_local],
            grad_outputs=[grad_loss_local],
            retain_graph=False,
        )

        # Wrap as DTensor
        grad_pred_pde = DTensor.from_local(
            grad_pred_pde_local,
            device_mesh=ctx.device_mesh,
            placements=ctx.pred_pde_placements,
            shape=ctx.pred_pde_shape,
            stride=ctx.pred_pde_stride,
        )

        return grad_pred_pde, None, None, None, None, None, None, None


def pde_loss(
    pred_pde: DTensor,
    pred_atom_coords: DTensor,
    true_atom_coords: DTensor,
    true_coords_resolved_mask: DTensor,
    feats: dict[str, DTensor],
    comm: TransposeComm,
    multiplicity: int = 1,
    max_dist: float = 32.0,
) -> DTensor:
    """Compute the PDE loss using DTensor.

    This is the DTensor version of boltz.model.loss.confidencev2.pde_loss.
    It computes the cross-entropy loss between predicted PDE bins and target
    PDE scores computed from pairwise coordinate distances.

    The entire computation is fused into _PDELossImpl using native PyTorch ops:
    1. Coordinate mapping via einsum (token_to_rep_atom @ coords)
    2. redistribute_transpose for column tensors
    3. cdist_pde Triton kernel (fused distance + binning + cross-entropy + masking + sum)
    4. all_reduce(SUM) across CP group for numerator/denominator
    5. Normalization: numerator / clamp(denominator, min=eps)
    6. Batch sum + all_reduce(SUM) across DP group
    7. Mean scaling by 1/batch_size

    All all_reduce ops use SUM which has identity gradient, allowing the entire
    computation to be captured in a single autograd subgraph.

    Parameters
    ----------
    pred_pde : DTensor
        Predicted PDE logits with shape (B * multiplicity, N_token, N_token, num_bins).
        Placements: (Shard(0), Shard(1), Shard(2)).
    pred_atom_coords : DTensor
        Predicted atom coordinates with shape (B * multiplicity, N_atom, 3).
        Placements: (Shard(0), Shard(1), Replicate()).
    true_atom_coords : DTensor
        Ground truth atom coordinates with shape (B * multiplicity, N_atom, 3).
        Placements: (Shard(0), Shard(1), Replicate()).
    true_coords_resolved_mask : DTensor
        Mask for resolved coordinates with shape (B * multiplicity, N_atom).
        Placements: (Shard(0), Shard(1), Replicate()).
    feats : dict[str, DTensor]
        Dictionary containing feature tensors:
        - "token_to_rep_atom": [B, N_token, N_atom_max_per_shard]
    comm : TransposeComm
        Communication object for redistribute_transpose operations.
    multiplicity : int, optional
        Diffusion batch multiplier, by default 1.
    max_dist : float, optional
        Maximum distance for binning, by default 32.0.

    Returns
    -------
    DTensor
        Scalar loss with placements (Replicate(), Replicate(), Replicate()).
    """
    # Extract token_to_rep_atom from features
    token_to_rep_atom = feats["token_to_rep_atom"]

    # _PDELossImpl fuses the entire loss computation using native PyTorch ops:
    # cdist_pde -> all_reduce(CP) -> normalize -> sum -> all_reduce(DP) -> mean
    loss = _PDELossImpl.apply(
        pred_pde,
        pred_atom_coords,
        true_atom_coords,
        true_coords_resolved_mask,
        token_to_rep_atom,
        comm,
        multiplicity,
        max_dist,
    )

    return loss


def confidence_loss(
    model_out: dict[str, DTensor],
    feats: dict[str, DTensor],
    true_coords: DTensor,
    true_coords_resolved_mask: DTensor,
    comm: TransposeComm,
    token_level_confidence: bool = True,
    multiplicity: int = 1,
    alpha_pae: float = 0.0,
    mask_loss: Optional[DTensor] = None,
    relative_supervision_weight: float = 0.0,
    dist_manager: Optional[DistributedManager] = None,
    group_layout: Optional[LayoutMap] = None,
) -> dict[str, DTensor | dict[str, DTensor]]:
    """Compute confidence loss using DTensor operations.

    This is the DTensor-compatible version of boltz.model.loss.confidencev2.confidence_loss.
    It aggregates plddt, pde, resolved, and (optionally) pae losses.

    The sub-loss implementations (plddt_loss, pde_loss, resolved_loss) operate at
    token level, matching the Boltz-2 ``token_level_confidence=True`` setting.

    Parameters
    ----------
    model_out : dict[str, DTensor]
        Dictionary containing the model output DTensors:
        - "plddt_logits": Shape [B*mult, N_token, num_bins], Placements (Shard(0), Shard(1), Replicate())
        - "pde_logits": Shape [B*mult, N_token, N_token, num_bins], Placements (Shard(0), Shard(1), Shard(2))
        - "resolved_logits": Shape [B*mult, N_token, 2], Placements (Shard(0), Shard(1), Replicate())
        - "sample_atom_coords": Shape [B*mult, N_atom, 3], Placements (Shard(0), Shard(1), Replicate())
        - "pae_logits" (when alpha_pae > 0): Shape [B*mult, N_token, N_token, num_bins],
          Placements (Shard(0), Shard(1), Shard(2))
    feats : dict[str, DTensor]
        Dictionary containing the model input DTensors:
        - "token_to_rep_atom": Shape [B, N_token, N_atom_max_per_shard], Placements (Shard(0), Shard(1), Replicate())
        - "token_pad_mask": Shape [B, N_token], Placements (Shard(0), Shard(1), Replicate())
        - "r_set_to_rep_atom": Shape [B, N_R, N_atom_max_per_shard], Placements (Shard(0), Shard(1), Replicate())
        - "atom_to_token": Shape [B, N_atom, N_token_max_per_shard], Placements (Shard(0), Shard(1), Replicate())
        - "mol_type": Shape [B, N_token], Placements (Shard(0), Shard(1), Replicate())
    true_coords : DTensor
        The atom coordinates after symmetry correction.
        Shape [B*mult, N_atom, 3], Placements (Shard(0), Shard(1), Replicate())
    true_coords_resolved_mask : DTensor
        The resolved mask after symmetry correction.
        Shape [B*mult, N_atom], Placements (Shard(0), Shard(1), Replicate())
    comm : TransposeComm
        Communication object for redistribute_transpose operations.
    token_level_confidence : bool, optional
        Must be True (default). The atom-level path (False) is not implemented.
    multiplicity : int, optional
        The diffusion batch size, by default 1
    alpha_pae : float, optional
        The weight of the pae loss, by default 0.0.
    mask_loss : DTensor, optional
        Per-sample loss mask. Not yet implemented; must be None.
    relative_supervision_weight : float, optional
        Weight for relative confidence supervision. Not yet implemented; must be 0.0.
    dist_manager : DistributedManager, optional
        Required when alpha_pae > 0.0 for pae_loss communication.
    group_layout : LayoutMap, optional
        Required when alpha_pae > 0.0 for pae_loss 2D CP grid layout.

    Returns
    -------
    dict[str, DTensor | dict[str, DTensor]]
        Dictionary containing:
        - "loss": Scalar DTensor with total loss, Placements (Replicate(), Replicate(), Replicate())
        - "loss_breakdown": dict with individual loss DTensors:
            - "plddt_loss": Scalar DTensor
            - "pde_loss": Scalar DTensor
            - "resolved_loss": Scalar DTensor
            - "pae_loss": Scalar DTensor

    See Also
    --------
    boltz.model.loss.confidencev2.confidence_loss : Serial version.
    plddt_loss : DTensor pLDDT loss computation.
    pde_loss : DTensor PDE loss computation.
    resolved_loss : DTensor resolved loss computation.
    pae_loss : DTensor PAE loss computation.
    """
    if not token_level_confidence:
        raise NotImplementedError(
            "confidence_loss only supports token_level_confidence=True. "
            "The atom-level confidence path is not implemented for DTensor."
        )
    if mask_loss is not None:
        raise NotImplementedError("confidence_loss does not yet support mask_loss. " "Pass mask_loss=None (default).")
    if relative_supervision_weight != 0.0:
        raise NotImplementedError(
            "confidence_loss does not yet support relative_supervision_weight != 0.0. "
            f"Got {relative_supervision_weight}."
        )
    if alpha_pae > 0.0:
        pae = pae_loss(
            model_out["pae_logits"],
            model_out["sample_atom_coords"],
            true_coords,
            true_coords_resolved_mask,
            feats,
            comm,
            dist_manager,
            group_layout,
            multiplicity,
        )
    else:
        device_mesh = model_out["plddt_logits"].device_mesh
        pae = DTensor.from_local(
            torch.tensor(0.0, device=model_out["plddt_logits"].device),
            device_mesh=device_mesh,
            placements=(Replicate(), Replicate(), Replicate()),
            shape=torch.Size(()),
            stride=(),
        )

    # Compute plddt loss
    plddt = plddt_loss(
        pred_lddt=model_out["plddt_logits"],
        pred_atom_coords=model_out["sample_atom_coords"],
        true_atom_coords=true_coords,
        true_coords_resolved_mask=true_coords_resolved_mask,
        feats=feats,
        comm=comm,
        multiplicity=multiplicity,
    )

    # Compute pde loss
    pde = pde_loss(
        pred_pde=model_out["pde_logits"],
        pred_atom_coords=model_out["sample_atom_coords"],
        true_atom_coords=true_coords,
        true_coords_resolved_mask=true_coords_resolved_mask,
        feats=feats,
        comm=comm,
        multiplicity=multiplicity,
    )

    # Compute resolved loss
    resolved = resolved_loss(
        pred_resolved=model_out["resolved_logits"],
        feats=feats,
        true_coords_resolved_mask=true_coords_resolved_mask,
        multiplicity=multiplicity,
    )

    # Sum the losses: loss = plddt + pde + resolved + alpha_pae * pae
    loss = elementwise_op(plddt, pde, ElementwiseOp.SUM)
    loss = elementwise_op(loss, resolved, ElementwiseOp.SUM)
    if alpha_pae > 0.0:
        pae_scaled = scalar_tensor_op(alpha_pae, pae, ElementwiseOp.PROD)
        loss = elementwise_op(loss, pae_scaled, ElementwiseOp.SUM)

    # Build output dictionary
    dict_out = {
        "loss": loss,
        "loss_breakdown": {
            "plddt_loss": plddt,
            "pde_loss": pde,
            "resolved_loss": resolved,
            "pae_loss": pae,
        },
    }

    return dict_out
