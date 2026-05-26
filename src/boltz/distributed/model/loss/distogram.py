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

"""DTensor-based context-parallel distogram loss for Boltz-2.

Implements the distogram loss as a single torch.autograd.Function, extracting
local tensors once and performing all computation locally with explicit
all_reduce calls at the required communication points.

Design rationale for a single autograd.Function:
- Single to_local()/from_local() pair instead of one per DTensor op
- One autograd graph node instead of one per op
- All local math is plain PyTorch and could be torch.compiled
- Minimal dependencies (only TransposeComm for the mask outer product)

Communication budget:
  Forward (4 collective calls):
    1. transpose_then_redistribute for mask outer product (1 call)
    2. all_reduce over CP group for denom (async, overlaps compute) (1 call)
    3. all_reduce over CP group for total (1 call)
    4. all_reduce over DP group for batch mean (1 call)
  Backward (0 collective calls):
    The backward of all_reduce(SUM) is identity — each rank computes
    gradients for its own local spatial chunk with no communication.

Equivalence to serial code (src/boltz/model/loss/distogramv2.py):
- Aggregate mode: sum+normalize target over K → K_eff=1, compute loss against
  each of D predictions, take min over D. For D=1 this is identical to the
  serial aggregate path (min/mean over size-1 dims are identity ops).
- Non-aggregate mode: keep K conformers (K_eff=K), compute all K_eff×D
  cross-entropies vectorized, min over D, mean over K.

Named distogram.py (not distogramv2) because v2 only differs from v1 by the
extra num_distogram axis, which falls back to v1 via unsqueeze(3).
"""

import torch
import torch.distributed as dist
from torch.autograd.function import FunctionCtx
from torch.distributed.tensor import DTensor, Partial, Replicate, Shard

from boltz.distributed.comm import TransposeComm
from boltz.distributed.model.layers.redistribute_transpose_without_dtensor import transpose_then_redistribute


def _build_pairwise_mask_local(
    mask_local: torch.Tensor,
    comm: TransposeComm,
) -> torch.Tensor:
    """Build pairwise [B, N_row, N_col] boolean mask from token mask [B, N_local].

    This performs the R2S outer-BITAND: mask_i & mask_j where mask_j comes from
    the transposed rank. Also zeros the diagonal on self-comm ranks.

    No gradients needed — the mask is treated as a constant.
    """
    # mask_local: [B_local, N_local]
    # Expand for outer product: [B_local, N_local, 1]
    mask_expanded = mask_local.unsqueeze(2)
    # Get the transposed chunk: [B_local, 1, N_local_col]
    mask_transposed = transpose_then_redistribute(mask_expanded, 1, 2, comm)
    # Outer BITAND: [B_local, N_local_row, N_local_col]
    mask_2d = mask_expanded & mask_transposed

    # Zero diagonal on self-comm ranks (where row and column ranges overlap)
    if comm.is_self_comm:
        local_n = mask_local.shape[1]
        diag_mask = ~torch.eye(local_n, dtype=torch.bool, device=mask_local.device)
        mask_2d = mask_2d & diag_mask.unsqueeze(0)

    return mask_2d


class _DistogramLossCP(torch.autograd.Function):
    """Single autograd.Function for the full distogram loss.

    Forward: to_local() → local math with explicit all_reduces → from_local()
    Backward: local math only (no communication)
    """

    @staticmethod
    def forward(
        ctx: FunctionCtx,
        pred: DTensor,
        target: DTensor,
        mask_token: DTensor,
        comm: TransposeComm,
        aggregate_distogram: bool,
    ) -> tuple[DTensor, DTensor]:
        """Forward pass.

        Parameters
        ----------
        pred : DTensor
            Prediction logits [B, N, N, D, bins], placements (Shard(0), Shard(1), Shard(2)).
        target : DTensor
            Target distributions [B, N, N, K, bins], placements (Shard(0), Shard(1), Shard(2)).
        mask_token : DTensor
            Token validity mask [B, N], placements (Shard(0), Shard(1), Replicate()).
        comm : TransposeComm
            Communication handle for CP transpose operations.
        aggregate_distogram : bool
            Whether to aggregate target over K conformers.

        Returns
        -------
        global_loss : DTensor
            Scalar loss, placements (Replicate(), Replicate(), Replicate()).
        batch_loss : DTensor
            Per-example loss [B], placements (Shard(0), Replicate(), Replicate()).
        """
        # --- Validate inputs ---
        if not isinstance(pred, DTensor):
            raise TypeError(f"pred must be DTensor, got {type(pred)}")
        if not isinstance(target, DTensor):
            raise TypeError(f"target must be DTensor, got {type(target)}")
        if not isinstance(mask_token, DTensor):
            raise TypeError(f"mask_token must be DTensor, got {type(mask_token)}")

        device_mesh = pred.device_mesh

        for name, dtensor in [("target", target), ("mask_token", mask_token)]:
            if dtensor.device_mesh != device_mesh:
                raise ValueError(f"{name} has different device_mesh than pred")

        expected_pairlike = (Shard(0), Shard(1), Shard(2))
        if pred.placements != expected_pairlike:
            raise ValueError(f"pred placements {pred.placements} must be {expected_pairlike}")
        if target.placements != expected_pairlike:
            raise ValueError(f"target placements {target.placements} must be {expected_pairlike}")
        expected_mask = (Shard(0), Shard(1), Replicate())
        if mask_token.placements != expected_mask:
            raise ValueError(f"mask_token placements {mask_token.placements} must be {expected_mask}")
        for name, dtensor in [("pred", pred), ("target", target), ("mask_token", mask_token)]:
            # Check placements and handle sharded dimensions
            for i_dim, placement in enumerate(dtensor.placements):
                if isinstance(placement, Partial):
                    raise ValueError(f"Partial placement on {name} mesh dim {i_dim} is not supported")
                elif isinstance(placement, Shard):
                    # Check that sharded dimensions are evenly divided
                    if dtensor.shape[placement.dim] % device_mesh.shape[i_dim] != 0:
                        raise ValueError(
                            f"Uneven sharding {name} tensor dimension {placement.dim} of size {dtensor.shape[placement.dim]} "
                            f"along device mesh dimension {i_dim} of size {device_mesh.shape[i_dim]} is not supported"
                        )
        if pred.ndim != 5:  # noqa: PLR2004
            raise ValueError(f"pred must be 5D [B, N, N, D, bins], got {pred.ndim}D")
        if target.ndim != 5:  # noqa: PLR2004
            raise ValueError(f"target must be 5D [B, N, N, K, bins], got {target.ndim}D")
        if mask_token.ndim != 2:  # noqa: PLR2004
            raise ValueError(f"mask_token must be 2D [B, N], got {mask_token.ndim}D")
        if pred.shape[0] != target.shape[0] or pred.shape[1] != target.shape[1] or pred.shape[2] != target.shape[2]:
            raise ValueError(f"pred shape {pred.shape} and target shape {target.shape} must match on dims 0,1,2")
        if pred.shape[4] != target.shape[4]:
            raise ValueError(f"pred bins {pred.shape[4]} != target bins {target.shape[4]}")
        if mask_token.shape[0] != pred.shape[0] or mask_token.shape[1] != pred.shape[1]:
            raise ValueError(f"mask_token shape {mask_token.shape} inconsistent with pred shape {pred.shape}")

        # --- Extract local tensors (single to_local() per input) ---
        # Compute in at least float32 for numerical stability, while respecting
        # float64 if either input is float64.
        compute_dtype = torch.promote_types(torch.promote_types(pred.dtype, target.dtype), torch.float32)
        pred_local = pred.to_local().to(compute_dtype)  # [B_local, N_row, N_col, D, bins]
        target_local = target.to_local().to(compute_dtype)  # [B_local, N_row, N_col, K, bins]
        mask_token_local = mask_token.to_local().to(torch.bool)  # [B_local, N_local]

        D = pred_local.shape[3]  # noqa: N806
        K = target_local.shape[3]  # noqa: N806

        # --- Build pairwise mask (involves R2S transpose, no grad needed) ---
        mask_local = _build_pairwise_mask_local(mask_token_local, comm).to(compute_dtype)
        # mask_local: [B_local, N_row, N_col]

        # --- Denom: launch async all_reduce so latency overlaps with compute ---
        cp_group = comm.group
        denom_local = mask_local.sum(dim=(-1, -2))  # [B_local]
        denom_work = dist.all_reduce(denom_local, op=dist.ReduceOp.SUM, group=cp_group, async_op=True)

        # --- Target preparation ---
        if aggregate_distogram:
            # Sum over K conformers, normalize → single effective conformer
            P_local = target_local.sum(dim=3)  # [B,N_r,N_c,bins]
            P_denom = P_local.sum(dim=-1, keepdim=True).clamp(min=1)  # [B,N_r,N_c,1]
            P_local = P_local / P_denom  # [B,N_r,N_c,bins]
            P_local = P_local.unsqueeze(3)  # [B,N_r,N_c,1,bins]
            K_eff = 1  # noqa: N806
        else:
            P_local = target_local  # [B,N_r,N_c,K,bins]
            K_eff = K  # noqa: N806

        # --- Vectorized cross-entropy for all (k, d) pairs ---
        log_Q_local = torch.nn.functional.log_softmax(pred_local, dim=-1)  # [B,N_r,N_c,D,bins]
        softmax_local = log_Q_local.exp()  # save for backward

        # Expand P: [B,N_r,N_c,K_eff,bins] → [B,N_r,N_c,K_eff,D,bins]
        P_expanded = P_local.unsqueeze(4).expand(-1, -1, -1, -1, D, -1)
        # Expand log_Q: [B,N_r,N_c,D,bins] → [B,N_r,N_c,K_eff,D,bins]
        log_Q_expanded = log_Q_local.unsqueeze(3).expand(-1, -1, -1, K_eff, -1, -1)

        # Cross-entropy: -sum(P * log_Q, dim=bins) → [B,N_r,N_c,K_eff,D]
        errors_local = -(P_expanded * log_Q_expanded).sum(dim=-1)

        # --- Flatten K_eff*D, apply mask, spatial reduction ---
        errors_flat = errors_local.reshape(
            errors_local.shape[0], errors_local.shape[1], errors_local.shape[2], K_eff * D
        )  # [B,N_r,N_c,K_eff*D]
        mask_exp = mask_local.unsqueeze(-1).expand(-1, -1, -1, K_eff * D)  # [B,N_r,N_c,K_eff*D]
        masked = errors_flat * mask_exp  # [B,N_r,N_c,K_eff*D]

        # --- Separate CP reductions for total and denom ---
        total_local = masked.sum(dim=(1, 2))  # [B_local, K_eff*D]
        dist.all_reduce(total_local, op=dist.ReduceOp.SUM, group=cp_group)
        denom_work.wait()
        denom_local = denom_local + 1e-5  # [B_local]

        # --- Divide by denom ---
        total_local = total_local / denom_local.unsqueeze(-1)  # [B_local, K_eff*D]

        # --- Min over D, mean over K_eff ---
        batch_loss_kd = total_local.reshape(total_local.shape[0], K_eff, D)  # [B_local, K_eff, D]
        min_result = torch.min(batch_loss_kd, dim=-1)  # values: [B_local, K_eff], indices saved
        batch_loss_k = min_result.values
        min_indices = min_result.indices  # [B_local, K_eff] — needed for backward

        batch_loss_local = batch_loss_k.sum(dim=-1) / K_eff  # [B_local]

        # --- Global loss: mean over batch (all_reduce over DP dim) ---
        B_global = pred.shape[0]  # noqa: N806
        global_loss_local = batch_loss_local.sum(dim=0, keepdim=False)  # scalar
        dp_group = device_mesh.get_group(0)
        dist.all_reduce(global_loss_local, op=dist.ReduceOp.SUM, group=dp_group)
        global_loss_local = global_loss_local / B_global

        # --- Save for backward ---
        if pred.requires_grad:
            ctx.save_for_backward(
                softmax_local,  # [B,N_r,N_c,D,bins]
                P_local,  # [B,N_r,N_c,K_eff,bins]
                mask_local,  # [B,N_r,N_c]
                denom_local,  # [B_local]
                min_indices,  # [B_local, K_eff]
            )
            ctx.D = D
            ctx.K_eff = K_eff
            ctx.B_global = B_global
            ctx.device_mesh = device_mesh
            ctx.pred_placements = pred.placements
            ctx.pred_shape = pred.shape
            ctx.pred_stride = pred.stride()

        # --- Wrap results as DTensors ---
        # batch_loss: [B] sharded on DP dim
        batch_loss_placements = (Shard(0), Replicate(), Replicate())
        bl_shape = (B_global,)
        bl_stride = (1,)
        batch_loss_dt = DTensor.from_local(
            batch_loss_local, device_mesh, batch_loss_placements, shape=bl_shape, stride=bl_stride
        )

        # global_loss: scalar, replicated
        global_loss_placements = (Replicate(), Replicate(), Replicate())
        gl_shape = ()
        gl_stride = ()
        global_loss_dt = DTensor.from_local(
            global_loss_local, device_mesh, global_loss_placements, shape=gl_shape, stride=gl_stride
        )

        return global_loss_dt, batch_loss_dt

    @staticmethod
    def backward(
        ctx: FunctionCtx, d_global_loss: DTensor, d_batch_loss: DTensor
    ) -> tuple[DTensor | None, None, None, None, None]:
        """Backward pass — entirely local, no collective communication.

        The forward's all_reduces produce replicated results. Their backward is
        broadcast (expand), which is local since each rank only fills its own chunk.
        """
        if not ctx.needs_input_grad[0]:
            return None, None, None, None, None

        softmax_local, P_local, mask_local, denom_local, min_indices = ctx.saved_tensors
        D = ctx.D  # noqa: N806
        K_eff = ctx.K_eff  # noqa: N806
        B_global = ctx.B_global  # noqa: N806
        device_mesh = ctx.device_mesh

        B_local = softmax_local.shape[0]  # noqa: N806
        N_row = softmax_local.shape[1]  # noqa: N806
        N_col = softmax_local.shape[2]  # noqa: N806

        # d_global_loss is Replicate scalar (DTensor), d_batch_loss may be
        # DTensor, plain Tensor, or None depending on what the caller backward'd through.
        compute_dtype = softmax_local.dtype
        d_gl = (d_global_loss.to_local() if isinstance(d_global_loss, DTensor) else d_global_loss).to(compute_dtype)
        if d_batch_loss is None:
            d_bl = torch.zeros(B_local, device=softmax_local.device, dtype=compute_dtype)
        elif isinstance(d_batch_loss, DTensor):
            d_bl = d_batch_loss.to_local().to(compute_dtype)
        else:
            d_bl = d_batch_loss.to(compute_dtype)

        # -----------------------------------------------------------------
        # Chain rule from global_loss = sum(batch_loss) / B_global
        # d_batch_loss_from_global = d_gl / B_global  (broadcast to [B_local])
        # total d_batch_loss_local = d_bl + d_gl / B_global
        # -----------------------------------------------------------------
        d_batch_local = d_bl + d_gl / B_global  # [B_local]

        # -----------------------------------------------------------------
        # Backward through mean over K_eff: batch_loss = sum(batch_loss_k) / K_eff
        # d_batch_loss_k = d_batch_local / K_eff  → [B_local, K_eff]
        # -----------------------------------------------------------------
        d_batch_loss_k = (d_batch_local / K_eff).unsqueeze(-1).expand(-1, K_eff)  # [B_local, K_eff]

        # -----------------------------------------------------------------
        # Backward through min over D: scatter gradient to argmin index
        # d_batch_loss_kd[b, k, d] = d_batch_loss_k[b, k] if d == min_indices[b, k] else 0
        # → [B_local, K_eff, D]
        # -----------------------------------------------------------------
        d_batch_loss_kd = torch.zeros(B_local, K_eff, D, device=d_batch_local.device, dtype=d_batch_local.dtype)
        d_batch_loss_kd.scatter_(-1, min_indices.unsqueeze(-1), d_batch_loss_k.unsqueeze(-1))

        # -----------------------------------------------------------------
        # Backward through reshape [B,K_eff,D] → [B,K_eff*D]
        # and division by denom: total = total_pre_div / denom
        # d_total_pre_div = d_total / denom
        # -----------------------------------------------------------------
        d_total = d_batch_loss_kd.reshape(B_local, K_eff * D)  # [B_local, K_eff*D]
        d_total = d_total / denom_local.unsqueeze(-1)  # [B_local, K_eff*D]

        # -----------------------------------------------------------------
        # Backward through spatial sum + all_reduce:
        # forward: total = all_reduce(sum(masked, dim=(1,2)))
        # backward of sum is broadcast: d_masked = d_total expanded to [B,N_r,N_c,K_eff*D]
        # backward of all_reduce(SUM) is identity (each rank gets the full gradient
        # and needs to compute its local contribution — no extra comm needed)
        # -----------------------------------------------------------------
        d_masked = d_total.unsqueeze(1).unsqueeze(2).expand(-1, N_row, N_col, -1)  # [B,N_r,N_c,K_eff*D]

        # Backward through mask multiply: d_errors_flat = d_masked * mask_exp
        mask_exp = mask_local.unsqueeze(-1).expand(-1, -1, -1, K_eff * D)
        d_errors_flat = d_masked * mask_exp  # [B,N_r,N_c,K_eff*D]

        # Reshape to [B,N_r,N_c,K_eff,D]
        d_errors = d_errors_flat.reshape(B_local, N_row, N_col, K_eff, D)

        # -----------------------------------------------------------------
        # Backward through cross-entropy: errors = -sum(P * log_Q, dim=-1)
        # d_log_Q_expanded = -P_expanded * d_errors.unsqueeze(-1)
        # Sum over K_eff to get d_log_Q
        # -----------------------------------------------------------------
        P_expanded = P_local.unsqueeze(4).expand(-1, -1, -1, -1, D, -1)  # [B,N_r,N_c,K_eff,D,bins]
        d_log_Q_expanded = -P_expanded * d_errors.unsqueeze(-1)  # [B,N_r,N_c,K_eff,D,bins]
        d_log_Q = d_log_Q_expanded.sum(dim=3)  # [B,N_r,N_c,D,bins]

        # -----------------------------------------------------------------
        # Backward through log_softmax:
        # If y = log_softmax(x), then dy/dx = diag(1) - softmax(x)
        # So d_pred = d_log_Q - softmax * sum(d_log_Q, dim=-1, keepdim=True)
        # -----------------------------------------------------------------
        d_pred_local = d_log_Q - softmax_local * d_log_Q.sum(dim=-1, keepdim=True)

        # --- Wrap gradient as DTensor ---
        d_pred = DTensor.from_local(
            d_pred_local,
            device_mesh=device_mesh,
            placements=ctx.pred_placements,
            shape=ctx.pred_shape,
            stride=ctx.pred_stride,
        )

        return d_pred, None, None, None, None


def distogram_loss(
    output: dict[str, DTensor],
    feats: dict[str, DTensor],
    comm: TransposeComm,
    aggregate_distogram: bool = True,
) -> tuple[DTensor, DTensor]:
    """Compute the distogram loss using a single fused autograd.Function.

    Both aggregate and non-aggregate modes share a unified pipeline:
      1. Prepare target P and set K_eff (differs by mode)
      2. Expand P and log_Q to [B, N, N, K_eff, D, bins]
      3. Compute cross-entropy for all (k, d) pairs
      4. Mask + spatial reduce with async denom + total all_reduce → [B, K_eff*D]
      5. Divide by denom → [B, K_eff*D]
      6. Unflatten → [B, K_eff, D], min over D, mean over K_eff → [B]
      7. Global loss (mean over DP-sharded batch dim)

    All local math is done in a single autograd.Function with explicit
    all_reduce calls, avoiding the overhead of ~15 separate DTensor
    autograd.Function round-trips.

    Boltz-1 (v1) compatibility:
      The v1 serial loss (src/boltz/model/loss/distogram.py) uses 4D tensors
      [B, N, N, bins] with no D or K axes. To use this function as a v1 loss,
      unsqueeze dim 3 of both pred and target before passing them in::

          output["pdistogram"] = pred_4d.unsqueeze(3)      # [B,N,N,bins] → [B,N,N,1,bins]
          feats["disto_target"] = target_4d.unsqueeze(3)    # [B,N,N,bins] → [B,N,N,1,bins]
          distogram_loss(output, feats, comm, aggregate_distogram=True)

      With D=1 and K=1, the min-over-D and mean-over-K steps are identity ops,
      producing results identical to v1.

    Parameters
    ----------
    output : dict[str, DTensor]
        Output of the model containing:
        - "pdistogram": [B, N, N, D, bins] prediction logits (DTensor)
    feats : dict[str, DTensor]
        Input features containing:
        - "disto_target": [B, N, N, K, bins] target distributions (DTensor)
        - "token_disto_mask": [B, N] token validity mask (DTensor)
    comm : TransposeComm
        Communication object for CP transpose and reduction.
    aggregate_distogram : bool
        If True, aggregates target over K conformers into a single normalized
        distribution (K_eff=1). Works with any D (generalizes serial D=1 constraint).
        If False, computes per-conformer loss (K_eff=K) and takes min over D.

    Returns
    -------
    DTensor
        The globally averaged loss (scalar DTensor).
    DTensor
        Per-example loss [B] (DTensor).
    """
    with torch.autocast("cuda", enabled=False):
        return _DistogramLossCP.apply(
            output["pdistogram"],
            feats["disto_target"],
            feats["token_disto_mask"],
            comm,
            aggregate_distogram,
        )
