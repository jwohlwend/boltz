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

"""DTensor-based context-parallel B-factor loss for Boltz-2.

Implements the B-factor loss as a single torch.autograd.Function, extracting
local tensors once and performing all computation locally with explicit
all_reduce calls at the required communication points.

The B-factor loss is per-token (no pairwise interaction), making it simpler
than the distogram loss.  The only differentiable input is ``pred`` (the
predicted B-factor logits from BFactorModule).

Equivalence to serial code (src/boltz/model/loss/bfactor.py):
  1. Map atom-level B-factors to tokens via token_to_rep_atom matrix
  2. Bin the per-token B-factors into a histogram (one-hot target)
  3. Compute cross-entropy loss between predicted logits and target
  4. Mask by valid tokens (bfactor > 1e-5), average over valid tokens

The serial code computes a single global fraction::

    loss = sum_{b,n}(errors * mask) / (sum_{b,n}(mask) + eps)

This is NOT per-batch-then-averaged.  The distributed version must match
this exactly by reducing both numerator and denominator globally across
all CP ranks and DP ranks before dividing.

Because ``pred`` has single-representation placements
(Shard(0), Shard(1), Replicate()), only the dp and cp0 mesh dimensions
carry unique data.  cp1 is Replicate — all cp1 ranks hold identical
local shards.  We reduce over dp and cp0 only (NOT cp1) to avoid
double-counting the replicated data.

Communication budget:
  Forward (2–3 all_reduce calls):
    1. all_reduce(SUM) over dp group for packed [loss_sum, mask_sum] (1 call)
    2. all_reduce(SUM) over cp0 group for the same packed tensor (1 call)
    Together equivalent to a single all_reduce over the combined dp×cp0
    group: sum_{dp×cp0}(x) = sum_{cp0}(sum_{dp}(x)).
    3. (optional) all_reduce(AVG) over cp1 group for the scalar loss,
    enforcing identical values across Replicate ranks (1 call, when cp1_group
    is provided).
  Backward (0 collective calls):
    The backward of all_reduce(SUM) is identity.
"""

import torch
import torch.distributed as dist
from torch.autograd.function import FunctionCtx
from torch.distributed import ProcessGroup
from torch.distributed.tensor import DTensor, Partial, Replicate, Shard
from torch.distributed.tensor.device_mesh import DeviceMesh


class _BFactorLossCP(torch.autograd.Function):
    """Single autograd.Function for the full B-factor loss.

    Forward: to_local() → local math with explicit all_reduces → from_local()
    Backward: local math only (no communication)
    """

    @staticmethod
    @torch.amp.custom_fwd(device_type="cuda")
    def forward(
        ctx: FunctionCtx,
        pred: DTensor,
        token_to_rep_atom: DTensor | torch.Tensor,
        bfactor: DTensor | torch.Tensor,
        device_mesh: DeviceMesh,
        dp_group: ProcessGroup,
        cp0_group: ProcessGroup,
        cp1_group: ProcessGroup | None,
    ) -> DTensor:
        """Forward pass.

        Parameters
        ----------
        pred : DTensor
            Predicted B-factor logits [B, N, bins], placements
            (Shard(0), Shard(1), Replicate()).
        token_to_rep_atom : DTensor | Tensor
            Token-to-representative-atom mapping [B, N_tokens, max_atoms_per_shard].
            Represents the diagonal block of the global one-hot tensor;
            the last axis is the local atom shard and is NOT sharded across CP.
        bfactor : DTensor | Tensor
            Per-atom B-factors [B, A].
        device_mesh : DeviceMesh
            3D device mesh (dp, cp0, cp1).
        dp_group : ProcessGroup
            Process group for the dp mesh dimension (dim 0).
        cp0_group : ProcessGroup
            Process group for the cp_axis_0 mesh dimension (dim 1).
        cp1_group : ProcessGroup | None
            Process group for the cp1 (Replicate) mesh dimension (dim 2).
            When not None, a mean all-reduce is applied to the scalar loss
            to enforce identical values across cp1 ranks.

        Returns
        -------
        global_loss : DTensor
            Scalar loss, placements (Replicate(), Replicate(), Replicate()).
        """
        # --- Validate differentiable input ---
        if not isinstance(pred, DTensor):
            raise TypeError(f"pred must be DTensor, got {type(pred)}")

        expected_single = (Shard(0), Shard(1), Replicate())
        if pred.placements != expected_single:
            raise ValueError(f"pred placements {pred.placements} must be {expected_single}")
        for i_dim, placement in enumerate(pred.placements):
            if isinstance(placement, Partial):
                raise ValueError(f"Partial placement on pred mesh dim {i_dim} is not supported")
            if isinstance(placement, Shard) and pred.shape[placement.dim] % device_mesh.shape[i_dim] != 0:
                raise ValueError(
                    f"Uneven sharding pred tensor dim {placement.dim} of size "
                    f"{pred.shape[placement.dim]} along mesh dim {i_dim} of size "
                    f"{device_mesh.shape[i_dim]}"
                )

        # Non-differentiable feature inputs: accept DTensor or plain tensor.
        # The data pipeline may provide custom per-shard atom slicing that
        # differs from standard distribute_tensor placements.
        for name, tensor in [("token_to_rep_atom", token_to_rep_atom), ("bfactor", bfactor)]:
            if not isinstance(tensor, (DTensor, torch.Tensor)):
                raise TypeError(f"{name} must be DTensor or Tensor, got {type(tensor)}")

        # --- Extract local tensors ---
        compute_dtype = torch.promote_types(pred.dtype, torch.float32)
        pred_local = pred.to_local().to(compute_dtype)  # [B_local, N_local, bins]
        t2ra_local = (token_to_rep_atom.to_local() if isinstance(token_to_rep_atom, DTensor) else token_to_rep_atom).to(
            compute_dtype
        )
        bf_local = (bfactor.to_local() if isinstance(bfactor, DTensor) else bfactor).to(compute_dtype)

        bins = pred_local.shape[2]

        # --- Construct target (non-differentiable) ---
        # Map atom-level bfactors to tokens
        bfactor_token = torch.bmm(t2ra_local, bf_local.unsqueeze(-1))  # [B_local, N_local, 1]

        # Bin into histogram
        boundaries = torch.linspace(0, 100, bins - 1, device=pred_local.device, dtype=compute_dtype)
        bfactor_token_bin = (bfactor_token > boundaries).sum(dim=-1).long()  # [B_local, N_local]
        bfactor_target = torch.nn.functional.one_hot(bfactor_token_bin, num_classes=bins).to(
            compute_dtype
        )  # [B_local, N_local, bins]

        # Token validity mask
        token_mask = (bfactor_token > 1e-5).squeeze(-1).to(compute_dtype)  # [B_local, N_local]

        # --- Cross-entropy loss ---
        log_softmax_local = torch.nn.functional.log_softmax(pred_local, dim=-1)
        softmax_local = log_softmax_local.exp()  # save for backward

        errors = -(bfactor_target * log_softmax_local).sum(dim=-1)  # [B_local, N_local]
        masked_errors = errors * token_mask  # [B_local, N_local]

        # --- Global reduction matching serial semantics ---
        # Serial: loss = sum_{b,n}(errors * mask) / (sum_{b,n}(mask) + eps)
        # We sum over ALL local dims, then all_reduce across dp and cp0
        # in two sequential calls.  cp1 is Replicate for single-representation
        # data, so cp1 ranks hold identical values — reducing over cp1 would
        # double-count.
        #
        # Two sequential all_reduces are equivalent to a single all_reduce
        # over the combined dp×cp0 group:
        #   sum_{dp×cp0}(x) = sum_{cp0}(sum_{dp}(x))
        packed = torch.stack([masked_errors.sum(), token_mask.sum()])
        dist.all_reduce(packed, op=dist.ReduceOp.SUM, group=dp_group)
        dist.all_reduce(packed, op=dist.ReduceOp.SUM, group=cp0_group)

        global_denom = packed[1] + 1e-5
        global_loss_local = packed[0] / global_denom

        # Average across cp1 (Replicate) ranks to enforce identical scalar
        # loss values.  This relies on pred having Replicate() on the cp1 mesh
        # dimension (validated above as expected_single[2]).  If the Replicate
        # axis were on a different mesh dimension, the group here would need
        # to match that dimension instead.
        if cp1_group is not None:
            dist.all_reduce(global_loss_local, op=dist.ReduceOp.AVG, group=cp1_group)

        # --- Save for backward ---
        if pred.requires_grad:
            ctx.save_for_backward(
                softmax_local,
                bfactor_target,
                token_mask,
                global_denom.unsqueeze(0),  # wrap scalar in 1D for save_for_backward
            )
            ctx.device_mesh = device_mesh
            ctx.pred_placements = pred.placements
            ctx.pred_shape = pred.shape
            ctx.pred_stride = pred.stride()

        # --- Wrap result as DTensor ---
        global_loss_placements = (Replicate(), Replicate(), Replicate())
        return DTensor.from_local(
            global_loss_local,
            device_mesh,
            global_loss_placements,
            shape=(),
            stride=(),
        )

    @staticmethod
    @torch.amp.custom_bwd(device_type="cuda")
    def backward(ctx: FunctionCtx, d_global_loss: DTensor) -> tuple[DTensor | None, None, None, None, None, None, None]:
        """Backward pass — entirely local, no collective communication.

        The gradient of all_reduce(SUM) is identity: each rank computes
        its own local gradient contribution.
        """
        if not ctx.needs_input_grad[0]:
            return None, None, None, None, None, None, None

        softmax_local, bfactor_target, token_mask, (global_denom,) = ctx.saved_tensors
        device_mesh = ctx.device_mesh

        d_gl = (d_global_loss.to_local() if isinstance(d_global_loss, DTensor) else d_global_loss).to(
            softmax_local.dtype
        )

        # Chain rule: loss = sum_{b,n}(errors * mask) / global_denom
        # errors = -sum(target * log_softmax(pred), dim=bins)
        # d_pred[b,n,c] = d_loss * mask[b,n] / global_denom * (softmax[b,n,c] - target[b,n,c])
        scale = d_gl / global_denom

        d_pred_local = (softmax_local - bfactor_target) * token_mask.unsqueeze(-1) * scale

        d_pred = DTensor.from_local(
            d_pred_local,
            device_mesh=device_mesh,
            placements=ctx.pred_placements,
            shape=ctx.pred_shape,
            stride=ctx.pred_stride,
        )

        return d_pred, None, None, None, None, None, None


def bfactor_loss(
    output: dict[str, DTensor],
    feats: dict[str, DTensor],
    device_mesh: DeviceMesh,
    dp_group: ProcessGroup,
    cp0_group: ProcessGroup,
    cp1_group: ProcessGroup | None = None,
) -> DTensor:
    """Compute the B-factor loss using a single fused autograd.Function.

    Parameters
    ----------
    output : dict[str, DTensor]
        Model outputs containing:
        - "pbfactor": [B, N, bins] predicted B-factor logits (DTensor).
    feats : dict[str, DTensor]
        Input features containing:
        - "token_to_rep_atom": [B, N_tokens, max_atoms_per_shard]
          token-to-atom mapping (DTensor).
        - "bfactor": [B, A] per-atom B-factors (DTensor).
    device_mesh : DeviceMesh
        3D device mesh (dp, cp0, cp1).
    dp_group : ProcessGroup
        Process group for the dp mesh dimension (dim 0).
    cp0_group : ProcessGroup
        Process group for the cp_axis_0 mesh dimension (dim 1).
    cp1_group : ProcessGroup | None
        Process group for the cp1 (Replicate) mesh dimension.
        When provided, a mean all-reduce is applied to the scalar loss
        to enforce identical values across cp1 ranks.

    Returns
    -------
    DTensor
        The globally averaged B-factor loss (scalar DTensor).
    """
    with torch.autocast("cuda", enabled=False):
        return _BFactorLossCP.apply(
            output["pbfactor"],
            feats["token_to_rep_atom"],
            feats["bfactor"],
            device_mesh,
            dp_group,
            cp0_group,
            cp1_group,
        )
