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
from typing import NamedTuple, Union

import torch
import torch.distributed as dist
from torch import Tensor
from torch.autograd.function import FunctionCtx
from torch.distributed.tensor import DTensor, Replicate, Shard
from torch.nn.attention import SDPBackend, sdpa_kernel

from boltz.distributed.comm import AttentionPairBiasComm
from boltz.distributed.model.layers.dtensor_metadata_tools import (
    raise_if_incorrect_dtensor_metadata_args,
)
from boltz.distributed.model.modules.utils import Precision, SDPAWithBiasBackend, setup_tf32_env
from boltz.distributed.utils import tiled_softmax_attention_update, update_exhaustive_strides

try:
    from torch.nn.attention.flex_attention import flex_attention

    flex_attention_compiled = torch.compile(flex_attention)
    HAS_FLEX_ATTN = True
except ImportError:
    flex_attention_compiled = None
    HAS_FLEX_ATTN = False


class _AttentionPairBiasContextVecParams(NamedTuple):
    """NamedTuple for attention pair bias context vector parameters."""

    ring_comm: AttentionPairBiasComm
    multiplicity: int
    num_heads: int
    head_dim: int
    inf: float
    use_window_batching: bool
    sdpa_with_bias_backend: SDPAWithBiasBackend = SDPAWithBiasBackend.TORCH_FLEX_ATTN


class _AttentionPairBiasShardwiseImpl(torch.autograd.Function):
    """Autograd function for shardwise attention with pair bias.

    This implements the forward and backward passes for window-batched attention
    with pair bias in a DTensor-compatible manner. The computation is performed
    locally on each shard without cross-rank communication (except for the implicit
    DTensor distribution).

    The attention computation follows:
        attn = softmax(q @ k.T / sqrt(head_dim) + z + mask_bias, dim=-2)
        o = attn @ v

    Where dimensions are:
        - q: (B * M, K, W, num_heads, head_dim) - queries per window
        - k: (B * M, K, H, num_heads, head_dim) - keys (H = full attention span)
        - v: (B * M, K, H, num_heads, head_dim) - values
        - z: (B, K, W, H, num_heads) - pair bias (broadcasts over M)
        - mask_key: (B, K, H) - key mask (broadcasts over M)

    The backward pass uses PyTorch autograd on the local computation graph,
    avoiding the need for manual gradient derivation.

    FIXME: bf16 is currently broken in _AttentionPairBiasShardwiseImpl when using activation checkpointing with torch flex attention

    See Also
    --------
    AttentionPairBiasShardwise : The nn.Module wrapper that calls this function.
    SDPAWithBiasBackend : Backend options for the attention computation.
    """

    @staticmethod
    def forward(
        ctx: FunctionCtx,
        q: DTensor,
        k: DTensor,
        v: DTensor,
        z: DTensor,
        mask_key: DTensor | None,
        sdpa_with_bias_backend: SDPAWithBiasBackend,
        num_heads: int,
        head_dim: int,
        inf: float,
    ) -> DTensor:
        """Forward pass for shardwise attention with pair bias.

        Computes multi-head attention with pair bias on window-batched inputs.
        Supports multiple backends for the core SDPA computation.

        Parameters
        ----------
        ctx : FunctionCtx
            The autograd context object for saving tensors for backward.
        q : DTensor
            Query tensor with shape (B * M, K, W, D) where:
            - B is batch size
            - M is multiplicity (diffusion samples)
            - K is number of windows
            - W is window size (typically 32)
            - D is hidden dimension (num_heads * head_dim)
        k : DTensor
            Key tensor with shape (B * M, K, H, D) where:
            - H is the full attention key dimension (typically 128)
        v : DTensor
            Value tensor with shape (B * M, K, H, D), same shape as k.
        z : DTensor
            Pair bias tensor with shape (B, K, W, H, num_heads).
            Note: Does not include multiplicity dimension; broadcasts over M.
        mask_key : DTensor or None
            Key mask tensor with shape (B, K, H) indicating valid key positions.
            None if no masking is needed.
        sdpa_with_bias_backend : SDPAWithBiasBackend
            Backend for computing scaled dot-product attention:
            - REFERENCE: Manual einsum implementation (most compatible)
            - TORCH_SDPA_EFFICIENT_ATTENTION: PyTorch's scaled_dot_product_attention kernel with EFFICIENT_ATTENTION backend
            - TORCH_FLEX_ATTN: PyTorch's FlexAttention with compiled score_mod
        num_heads : int
            Number of attention heads.
        head_dim : int
            Dimension per attention head.
        inf : float
            Large value used for masking invalid positions in attention.

        Returns
        -------
        DTensor
            Output tensor with shape (B * M, K, W, D).

        Raises
        ------
        ValueError
            If tensor dimensions don't match expected shapes.
            If q, k, v, z, mask_key have inconsistent placements or device meshes.
            If q.shape[0] is not divisible by z.shape[0] (multiplicity check).

        Notes
        -----
        - The softmax is computed over the key dimension (dim=-2 in attention matrix)
        - All inputs must have the same DTensor placements and device mesh
        - Computation is promoted to at least FP32 for numerical stability
        - The local computation graph is preserved for backward pass via autograd
        """

        if q.ndim != 4:  # (B * M, K, W, D)
            raise ValueError(f"Input q must have 4 dimensions. Got {q.ndim}.")

        if k.ndim != 4:  # (B * M, K, H, D)
            raise ValueError(f"Input k must have 4 dimensions. Got {k.ndim}.")

        if v.ndim != 4:  # (B * M, K, H, D)
            raise ValueError(f"Input v must have 4 dimensions. Got {v.ndim}.")

        if z.ndim != 5:  # (B, K, W, H, num_heads) with potentially no multiplicity
            raise ValueError(f"Input z must have 5 dimensions. Got {z.ndim}.")

        if mask_key is not None:
            if mask_key.ndim != 3:  # (B, K, H) with potentially no multiplicity
                raise ValueError(f"Input mask_key must have 3 dimensions. Got {mask_key.ndim}.")

            if mask_key.shape != z.shape[:2] + (z.shape[3],):
                raise ValueError(
                    f"Input mask_key must have the same shape as z.shape[:3]. Got {mask_key.shape} and {z.shape[:3]}."
                )

        # Shape checks on the input
        if q.shape[:2] != k.shape[:2] or q.shape[:2] != v.shape[:2]:  # B, K
            raise ValueError(
                f"Input q, k, v must have the same leading two B and K dimensions. Got {q.shape[:2]} and {k.shape[:2]} and {v.shape[:2]}"
            )

        if q.shape[-1] != k.shape[-1] or q.shape[-1] != v.shape[-1]:  # D
            raise ValueError(
                f"Input q, k, v must have the same last dimension. Got {q.shape[-1]} and {k.shape[-1]} and {v.shape[-1]}."
            )

        if q.shape[-1] != num_heads * head_dim:
            raise ValueError(
                f"Input q.shape[-1] and num_heads * head_dim must have the same shape. Got {q.shape[-1]} and {num_heads * head_dim}."
            )

        if k.shape != v.shape:
            raise ValueError(f"Input k and v must have the same shape. Got {k.shape} and {v.shape}.")

        if q.shape[0] % z.shape[0] != 0:  # B * M % B == 0
            raise ValueError(
                f"Input q.shape[0] must be a multiple of z.shape[0]. Got q.shape[0]={q.shape[0]} and z.shape[0]={z.shape[0]}."
            )

        if z.shape[1:3] != q.shape[1:3]:  # K, W
            raise ValueError(
                f"Input q.shape[1:3] and z.shape[1:3] must have the same shape. Got {q.shape[1:3]} and {z.shape[1:3]}."
            )

        if z.shape[3] != k.shape[2]:  # H
            raise ValueError(
                f"Input z.shape[3] and k.shape[2] must have the same shape. Got {z.shape[3]} and {k.shape[2]}."
            )

        if z.shape[-1] != num_heads:
            raise ValueError(
                f"Input z.shape[-1] and num_heads must have the same shape. Got {z.shape[-1]} and {num_heads}."
            )

        if sdpa_with_bias_backend == SDPAWithBiasBackend.TORCH_SDPA_EFFICIENT_ATTENTION and (
            q.shape[-1] % 4 != 0 or k.shape[-1] % 4 != 0 or v.shape[-1] % 4 != 0
        ):
            # torch SDPA errors are shown as warnings instead of errors so we raise for it instead
            raise ValueError(
                f"Torch SDPA Efficient Attention kernel requires q, k, v must have a last dimension that is divisible by 4. "
                f"Got {q.shape[-1]} and {k.shape[-1]} and {v.shape[-1]}."
            )

        # placements and device mesh checks
        if (
            q.placements != k.placements
            or q.placements != v.placements
            or q.placements != z.placements
            or q.placements != mask_key.placements
        ):
            raise ValueError(
                f"Input q, k, v, z, and mask must have the same placements. Got {q.placements} and {k.placements} and {v.placements} and {z.placements} and {mask_key.placements}."
            )
        if (
            q.device_mesh != k.device_mesh
            or q.device_mesh != v.device_mesh
            or q.device_mesh != z.device_mesh
            or q.device_mesh != mask_key.device_mesh
        ):
            raise ValueError(
                f"Input q, k, v, z, and mask must be on the same device mesh. Got {q.device_mesh} and {k.device_mesh} and {v.device_mesh} and {z.device_mesh} and {mask_key.device_mesh}."
            )

        multiplicity = q.shape[0] // z.shape[0]

        q_local_orig = q.to_local().detach().requires_grad_(q.requires_grad)
        k_local_orig = k.to_local().detach().requires_grad_(k.requires_grad)
        v_local_orig = v.to_local().detach().requires_grad_(v.requires_grad)
        z_local_orig = z.to_local().detach().requires_grad_(z.requires_grad)  # (B, K, W, H, num_heads)
        if mask_key is not None:
            mask_key_bias_local_orig = mask_key.to_local().detach().requires_grad_(False)  # (B, K, H)
        else:
            mask_key_bias_local_orig = None

        with torch.enable_grad():
            # enable grad to build a local graph for the shardwise operations
            # We detach inputs to create 'leaf' nodes for our local graph.
            # NOTE: mask_key_bias_local is not differentiable but nonetheless we need to detach it
            q_local = q_local_orig.unflatten(-1, (num_heads, head_dim))  # (B * M, K, W, num_heads, head_dim)
            k_local = k_local_orig.unflatten(-1, (num_heads, head_dim))  # (B * M, K, H, num_heads, head_dim)
            v_local = v_local_orig.unflatten(-1, (num_heads, head_dim))  # (B * M, K, H, num_heads, head_dim)
            z_local = z_local_orig

            if mask_key_bias_local_orig is not None:
                mask_key_bias_local = mask_key_bias_local_orig[:, :, None, :, None]  # (B, K, 1, H, 1)
                mask_key_bias_local = (1 - mask_key_bias_local.to(dtype=q_local.dtype)) * -inf
            else:
                mask_key_bias_local = None

            if multiplicity > 1:
                # unflatten the multiplicity axis so that mask and z are broadcasted along it
                # This has to be (B * multiplicity, ...) -> (B, multiplicity, ...)
                # but never (B * multiplicity, ...) -> (multiplicity, B, ...) due to the upstream
                # order of multiplicity application
                q_local = q_local.unflatten(0, (-1, multiplicity))
                k_local = k_local.unflatten(0, (-1, multiplicity))
                v_local = v_local.unflatten(0, (-1, multiplicity))
                # add singleton axis to mask and z for broadcasting
                z_local = z_local.unsqueeze(1)
                if mask_key_bias_local is not None:
                    mask_key_bias_local = mask_key_bias_local.unsqueeze(1)

            # use at least FP32 for AttnPairBias
            dtype_compute = torch.promote_types(q_local.dtype, torch.float32)
            q_local = q_local.to(dtype_compute)
            k_local = k_local.to(dtype_compute)
            v_local = v_local.to(dtype_compute)
            z_local = z_local.to(dtype_compute)
            if mask_key_bias_local is not None:
                mask_key_bias_local = mask_key_bias_local.to(dtype_compute)

            if (
                sdpa_with_bias_backend == SDPAWithBiasBackend.TORCH_SDPA_EFFICIENT_ATTENTION
                or sdpa_with_bias_backend == SDPAWithBiasBackend.TORCH_FLEX_ATTN
            ):
                # save shape for later reshaping the kernel output back
                if multiplicity > 1:
                    B_local, M_local, K_local, W_local = q_local.shape[:4]
                else:
                    B_local, K_local, W_local = q_local.shape[:3]
                # torch sdpa kernel only supports 4-axes input tensors so we need to
                # move up the head axis then flatten
                # (..., W or H, num_heads, head_dim) -> (..., num_heads, W or H, head_dim)
                q_local = q_local.transpose(-3, -2)
                k_local = k_local.transpose(-3, -2)
                v_local = v_local.transpose(-3, -2)
                # sdpa kernel only accepts bias so we need to sum mask and z into bias
                b_local = z_local
                if mask_key_bias_local is not None:
                    b_local = b_local + mask_key_bias_local
                # (..., W, H, num_heads) -> (..., num_heads, W, H)
                b_local = b_local.moveaxis(-1, -3)
                if multiplicity > 1:
                    # (B, M, K, H, ...) -> (M, B, K, H, ...)
                    q_local = q_local.moveaxis(1, 0).flatten(1, 3)
                    k_local = k_local.moveaxis(1, 0).flatten(1, 3)
                    v_local = v_local.moveaxis(1, 0).flatten(1, 3)
                    b_local = b_local.moveaxis(1, 0).flatten(1, 3)
                else:
                    # (B, K, H, ...) -> (B*K, H, ...)
                    q_local = q_local.flatten(0, 1)
                    k_local = k_local.flatten(0, 1)
                    v_local = v_local.flatten(0, 1)
                    b_local = b_local.flatten(0, 1)
                # run the kernel
                # NOTE: except for SDPBackend.MATH, other kernels can't guarantee consistent backward pass
                # results for the invalid atoms' gradients, which doesn't matter for applications but do matter
                # for testing requirements.
                if sdpa_with_bias_backend == SDPAWithBiasBackend.TORCH_SDPA_EFFICIENT_ATTENTION:
                    # NOTE: dtype_compute is at least FP32 so technically CUDNN_ATTENTION and FLASH_ATTENTION
                    # will not work.
                    with sdpa_kernel(backends=[SDPBackend.EFFICIENT_ATTENTION]):
                        # NOTE: the scale factor is already applied by the kernel, which is default to 1/sqrt(head_dim)
                        # Kernel requires all input to have stride-1 along the last axis
                        o_local = torch.nn.functional.scaled_dot_product_attention(
                            q_local, k_local, v_local, attn_mask=b_local.contiguous()
                        )
                elif sdpa_with_bias_backend == SDPAWithBiasBackend.TORCH_FLEX_ATTN:
                    # flex_attention (compiled Triton/Inductor) requires power-of-2 head_dim and head_dim >= 16.
                    if not (
                        HAS_FLEX_ATTN
                        and q_local.is_cuda
                        and q_local.dtype != torch.float64
                        and is_power_of_2(head_dim)
                        and head_dim >= 16
                    ):
                        raise RuntimeError(
                            f"flex_attention requirements not met: "
                            f"HAS_FLEX_ATTN={HAS_FLEX_ATTN}, is_cuda={q_local.is_cuda}, "
                            f"dtype={q_local.dtype}, head_dim={head_dim}, "
                            f"q_seq_len={q_local.size(-2)}, k_seq_len={k_local.size(-2)}"
                        )

                    if multiplicity > 1:
                        # Squeeze out the M=1 dimension so b_local is (B*K*H, Sq, Sk) to avoid
                        # data-dependent indexing of b[0, ...]
                        # (B * K * num_heads, W, H)
                        b_local = b_local.squeeze(0)

                        def add_bias_to_attn_score(
                            score: torch.Tensor,
                            batch: torch.Tensor,
                            head: torch.Tensor,
                            q_idx: torch.Tensor,
                            k_idx: torch.Tensor,
                        ) -> torch.Tensor:
                            return score + b_local[head, q_idx, k_idx]

                    else:
                        # b_local is (B * K, num_heads, W, H)
                        # with same batch size as q/k/v
                        def add_bias_to_attn_score(
                            score: torch.Tensor,
                            batch: torch.Tensor,
                            head: torch.Tensor,
                            q_idx: torch.Tensor,
                            k_idx: torch.Tensor,
                        ) -> torch.Tensor:
                            return score + b_local[batch, head, q_idx, k_idx]

                    with setup_tf32_env(Precision.FP32):
                        o_local = flex_attention_compiled(q_local, k_local, v_local, score_mod=add_bias_to_attn_score)
                # reshape the tensor back:
                if multiplicity > 1:
                    # (M, B * K * num_heads, ...) -> (B, M, K, ...)
                    o_local = o_local.unflatten(1, (B_local, K_local, num_heads)).moveaxis(0, 1)
                else:
                    # (B * K, ...) -> (B, K, ...)
                    o_local = o_local.unflatten(0, (B_local, K_local))
                # (..., num_heads, W, head_dim) -> (..., W, num_heads, head_dim)
                o_local = o_local.transpose(-3, -2)
            elif sdpa_with_bias_backend == SDPAWithBiasBackend.REFERENCE:
                with setup_tf32_env(Precision.FP32), torch.amp.autocast("cuda", enabled=False):
                    attn = torch.einsum("...wnd,...hnd->...whn", q_local, k_local)
                    attn = attn / head_dim**0.5
                    attn = attn + z_local
                    if mask_key_bias_local is not None:
                        attn = attn + mask_key_bias_local
                    attn = attn.softmax(dim=-2)  # axis = -2 is the key dimension
                    o_local = torch.einsum("...whn,...hnd->...wnd", attn, v_local)  # (B, K, W, num_heads, head_dim)
            # (..., W, num_heads, head_dim) -> (..., W, num_heads * head_dim)
            o_local = o_local.flatten(-2, -1).to(q.dtype)  # (B, K, W, c_s)

            if multiplicity > 1:
                # (B, multiplicity, ...) -> (B * multiplicity, ...)
                o_local = o_local.flatten(0, 1)

        # save the detached tensors for backward pass -- they hold the graph structure
        ctx.save_for_backward(q_local_orig, k_local_orig, v_local_orig, z_local_orig, mask_key_bias_local_orig, o_local)
        ctx.device_mesh = q.device_mesh
        ctx.placements = q.placements
        ctx.q_shape = q.shape
        ctx.q_stride = q.stride()
        ctx.k_shape = k.shape
        ctx.k_stride = k.stride()
        ctx.v_shape = v.shape
        ctx.v_stride = v.stride()
        ctx.z_shape = z.shape
        ctx.z_stride = z.stride()

        o_dtensor = DTensor.from_local(
            o_local.detach(),
            device_mesh=q.device_mesh,
            placements=q.placements,
            shape=q.shape,
            stride=q.stride(),
        )

        return o_dtensor

    @staticmethod
    def backward(
        ctx: FunctionCtx,
        grad_output: DTensor,
    ) -> tuple[DTensor | None, DTensor | None, DTensor | None, DTensor | None, None, None, None, None, None]:
        """Backward pass for shardwise attention with pair bias.

        Computes gradients by backpropagating through the local computation graph
        that was built during the forward pass. This leverages PyTorch's autograd
        rather than manual gradient computation.

        Parameters
        ----------
        ctx : FunctionCtx
            The autograd context containing saved tensors from forward:
            - q_local_orig, k_local_orig, v_local_orig, z_local_orig: Input tensors
            - mask_key_bias_local_orig: Mask bias (non-differentiable)
            - o_local: Output tensor that holds the computation graph
            - device_mesh, placements: DTensor metadata
            - q_shape, q_stride, etc.: Shape/stride info for DTensor reconstruction
        grad_output : DTensor
            Gradient of loss with respect to output, shape (B * M, K, W, D).

        Returns
        -------
        tuple[DTensor | None, ...]
            Gradients for each forward input in order:
            - dq: DTensor or None, gradient for q with shape (B * M, K, W, D)
            - dk: DTensor or None, gradient for k with shape (B * M, K, H, D)
            - dv: DTensor or None, gradient for v with shape (B * M, K, H, D)
            - dz: DTensor or None, gradient for z with shape (B, K, W, H, num_heads)
            - None: mask_key (non-differentiable)
            - None: sdpa_with_bias_backend (non-differentiable)
            - None: num_heads (non-differentiable)
            - None: head_dim (non-differentiable)
            - None: inf (non-differentiable)

        Notes
        -----
        The gradient computation follows the chain rule for attention:

        Forward (with einsum notation, ignoring multiplicity for clarity):
            q_local: (B, K, W, num_heads, head_dim) - "bkwid"
            k_local: (B, K, H, num_heads, head_dim) - "bkhid"
            v_local: (B, K, H, num_heads, head_dim) - "bkhid"
            attn = softmax(einsum("bkwid,bkhid->bkwhi", q, k) / sqrt(d) + z + mask_bias, dim=h)
            o = einsum("bkwhi,bkhid->bkwid", attn, v)

        Backward:
            dv = einsum("bkwhi,bkwid->bkhid", attn, grad_output)
            d_attn = einsum("bkwid,bkhid->bkwhi", grad_output, v)
            d_pre_softmax = attn * (d_attn - sum(attn * d_attn, dim=h, keepdim=True))
            dz = d_pre_softmax
            dq = einsum("bkwhi,bkhid->bkwid", d_pre_softmax, k) / sqrt(d)
            dk = einsum("bkwhi,bkwid->bkhid", d_pre_softmax, q) / sqrt(d)

        The actual implementation uses torch.autograd.grad on the saved local
        computation graph for correctness and maintainability.
        """

        # retrieve the leaf nodes for the local graph
        q_local, k_local, v_local, z_local, _, o_local = ctx.saved_tensors
        inputs_needing_grad = [t for t in (q_local, k_local, v_local, z_local) if t.requires_grad]
        if not inputs_needing_grad:
            # Short-circuit if nothing needed gradients (rare but possible)
            return None, None, None, None, None, None, None, None, None
        grad_output_local = grad_output.to_local()
        # backprop via the local graph -- grads_local only contains the grads for those in inputs_needing_grad
        with setup_tf32_env(Precision.FP32), torch.amp.autocast("cuda", enabled=False):
            grads_local = torch.autograd.grad(
                outputs=[o_local],
                inputs=inputs_needing_grad,
                grad_outputs=[grad_output_local],
                retain_graph=False,  # Frees the local graph immediately
            )

        iter_grads_local = iter(grads_local)
        dq_local = next(iter_grads_local) if q_local.requires_grad else None
        dk_local = next(iter_grads_local) if k_local.requires_grad else None
        dv_local = next(iter_grads_local) if v_local.requires_grad else None
        dz_local = next(iter_grads_local) if z_local.requires_grad else None

        if dq_local is not None:
            dq = DTensor.from_local(
                dq_local, device_mesh=ctx.device_mesh, placements=ctx.placements, shape=ctx.q_shape, stride=ctx.q_stride
            )
        else:
            dq = None

        if dk_local is not None:
            dk = DTensor.from_local(
                dk_local, device_mesh=ctx.device_mesh, placements=ctx.placements, shape=ctx.k_shape, stride=ctx.k_stride
            )
        else:
            dk = None

        if dv_local is not None:
            dv = DTensor.from_local(
                dv_local, device_mesh=ctx.device_mesh, placements=ctx.placements, shape=ctx.v_shape, stride=ctx.v_stride
            )
        else:
            dv = None

        if dz_local is not None:
            dz = DTensor.from_local(
                dz_local, device_mesh=ctx.device_mesh, placements=ctx.placements, shape=ctx.z_shape, stride=ctx.z_stride
            )
        else:
            dz = None

        return dq, dk, dv, dz, None, None, None, None, None


class _AttentionPairBiasContexVecImpl(torch.autograd.Function):
    @staticmethod
    def forward(
        ctx: FunctionCtx,
        q: DTensor,
        k: DTensor,
        v: DTensor,
        z: DTensor,
        mask: DTensor,
        pair_mask: Union[DTensor, None],
        apb_context_vec_params: _AttentionPairBiasContextVecParams,
    ) -> DTensor:
        """

        c_s = num_heads * head_dim checked in vanilla AttentionPairBias.__init__

        Below, N is the global number of tokens, H is the number of heads.

        Parameters
        ----------
        ctx: FunctionCtx
            The context object.
        q : DTensor
            query vectors computed by projection, (B, N, c_s)
        k : DTensor
            key vectors computed by projection, (B, N, c_s)
        v : DTensor
            value vectors computed by projection, (B, N, c_s)
        z : DTensor
            z, (B, H, N, N)
        mask : torch.Tensor
            The pairwise mask tensor (B, N)
        multiplicity : int, optional
            The diffusion batch size, by default 1
        pair_mask: Union[DTensor, None]
            The pairwise mask tensor.
        apb_context_vec_params: tuple
            The parameters for the attention pair bias context vector.

        Key features:
            - Distributed computation across device meshes with various sharding strategies
            - Memory-efficient implementation that operates on local tensor chunks
            - Supports gradient computation through custom backward pass
            - Validates tensor compatibility (type, device mesh, placements, shapes)`

        Raises
        ------
        TypeError
            If dtensor_instance is not a DTensor.
        ValueError
            If the DTensor metadata is incorrect.
        """
        # Check input metadata
        _AttentionPairBiasContexVecImpl.check_forward_input_metadata_and_store(
            ctx, q, k, v, z, mask, pair_mask, apb_context_vec_params
        )
        # Check implementation scope
        _AttentionPairBiasContexVecImpl.check_forward_input_for_impl_state(apb_context_vec_params)

        # ----------------------------------------------
        # Setup inputs to RingAttention.forward()
        # -------------------------------------------------------
        (
            ring_comm,
            multiplicity,
            _,
            _,
            inf,
            use_window_batching,
            sdpa_with_bias_backend,
        ) = apb_context_vec_params

        requires_grad = any(p.requires_grad for p in (q, k, v, z))

        ctx.mark_non_differentiable(mask)
        mask_local: Tensor = mask.to_local()

        # overlay mask comm with qkv projection, pair_mask_ij <- pair_mask_ij + mask_j
        mask_recv: Tensor = ring_comm.comm_transpose_mask.enqueue_to_dispatch(mask_local.contiguous())

        dtype_input = q.dtype
        dtype_compute = torch.promote_types(dtype_input, torch.float32)
        ctx.dtype_compute = dtype_compute
        ctx.dtype_input = dtype_input

        q_local: Tensor = q.to_local().to(dtype=dtype_compute)
        k_local: Tensor = k.to_local().to(dtype=dtype_compute)
        v_local: Tensor = v.to_local().to(dtype=dtype_compute)
        z_local: Tensor = z.to_local().to(dtype=dtype_compute)

        ctx.B_each_chunk = q_local.shape[0]
        ctx.N_each_chunk = q_local.shape[1]

        single_rep_view_shape = ctx.B_each_chunk, ctx.N_each_chunk, ctx.H, ctx.head_dim
        q_local = q_local.view(single_rep_view_shape).requires_grad_(q.requires_grad)
        k_local = k_local.view(single_rep_view_shape).requires_grad_(k.requires_grad)
        v_local = v_local.view(single_rep_view_shape).requires_grad_(v.requires_grad)
        z_local = z_local.permute(0, 3, 1, 2).requires_grad_(z.requires_grad)

        if requires_grad:
            ctx.multiplicity = multiplicity

        ring_comm.comm_transpose_mask.wait_until_finished()

        if use_window_batching or pair_mask is None:  # original behavior
            pair_mask_local = mask_recv[:, None, None, :]
        else:  # only atom-level has pair_mask
            ctx.mark_non_differentiable(pair_mask)
            pair_mask_local: Tensor = pair_mask.to_local()  # shape = (B, I, J)
            pair_mask_local = pair_mask_local[:, None, :, :] * mask_recv[:, None, None, :]

        pair_mask_local = pair_mask_local.to(dtype=dtype_compute)

        with torch.autocast("cuda", enabled=False):
            o_local, ring_attention_simple_data_for_bw = ring_attention_simple_forward(
                q_local,
                k_local,
                v_local,
                z_local,
                pair_mask_local,
                ring_comm,
                inf,
                sdpa_with_bias_backend,
            )
            if requires_grad:
                # Unpack tensors for save_for_backward to enable automatic memory management and hook support
                ctx.save_for_backward(
                    ring_attention_simple_data_for_bw.q_store,
                    ring_attention_simple_data_for_bw.k_t_store,
                    ring_attention_simple_data_for_bw.v_t_store,
                    ring_attention_simple_data_for_bw.z_store,
                    ring_attention_simple_data_for_bw.lse_m,
                    ring_attention_simple_data_for_bw.o_store,
                )
                ctx.ring_comm = ring_attention_simple_data_for_bw.ring_comm
                ctx.sdpa_with_bias_backend = ring_attention_simple_data_for_bw.sdpa_with_bias_backend

        # ---------------------------------------------------------
        # end custom communication
        # ---------------------------------------------------------
        o_local = o_local.reshape(ctx.B_each_chunk, ctx.N_each_chunk, ctx.c_s)  # o_local_b

        # Compute output shape and stride
        shape_output = v.shape[:-1] + (o_local.shape[-1],)

        strides_output = update_exhaustive_strides(o_local.shape, o_local.stride(), shape_output)

        o = DTensor.from_local(
            o_local.to(dtype=dtype_input),
            device_mesh=ctx.device_mesh,
            placements=ctx.single_rep_placements,
            shape=shape_output,
            stride=strides_output,
        )
        return o

    @staticmethod
    def backward(
        ctx: FunctionCtx,
        grad_output: DTensor,
    ) -> tuple[DTensor, DTensor, DTensor, DTensor, None, None, None]:
        """Backward pass implementation.

        Parameters
        ----------
        ctx: FunctionCtx
            The context object.
        grad_output: DTensor
            The gradient of the output tensor.

        Raises
        ------
        TypeError
            If dtensor_instance is not a DTensor.
        ValueError
            If the DTensor metadata is incorrect.
        """
        _AttentionPairBiasContexVecImpl.check_backward_input_metadata(ctx, grad_output)

        # Get ref to local tensor, and reshape to (B, I, H, D)
        grad_output_local: Tensor = grad_output.to_local().reshape(
            ctx.B_each_chunk, ctx.N_each_chunk, ctx.H, ctx.head_dim
        )
        grad_output_local = grad_output_local.to(dtype=ctx.dtype_compute)

        # Call backward separately via refactored-out function
        q_store, k_t_store, v_t_store, z_store, lse_m, o_store = ctx.saved_tensors
        data_for_backward = RingAttentionSimpleDataForBackward(
            q_store=q_store,
            k_t_store=k_t_store,
            v_t_store=v_t_store,
            z_store=z_store,
            lse_m=lse_m,
            ring_comm=ctx.ring_comm,
            sdpa_with_bias_backend=ctx.sdpa_with_bias_backend,
            o_store=o_store,
            multiplicity=ctx.multiplicity,
        )
        del ctx.ring_comm
        del ctx.sdpa_with_bias_backend

        grad_q, grad_k, grad_v, grad_z = ring_attention_simple_backward(
            data_for_backward=data_for_backward,
            do=grad_output_local,
        )
        del data_for_backward  # free up q_store, k_t_store, v_t_store, z_store, lse_m immediately
        grad_q: Tensor = grad_q.to(dtype=ctx.dtype_input)  # (B_each_chunk, N_each_chunk, H, D)
        grad_k: Tensor = grad_k.to(dtype=ctx.dtype_input)  # (B_each_chunk, N_each_chunk, H, D)
        grad_v: Tensor = grad_v.to(dtype=ctx.dtype_input)  # (B_each_chunk, N_each_chunk, H, D)
        grad_z: Tensor = grad_z.to(dtype=ctx.dtype_input)  # (B_each_chunk, N_each_chunk, H, D)

        # Reshape, allocate new memory
        single_rep_target_shape = (ctx.B_each_chunk, ctx.N_each_chunk, ctx.c_s)
        grad_q_flat: Tensor = grad_q.reshape(single_rep_target_shape)
        grad_k_flat: Tensor = grad_k.reshape(single_rep_target_shape)
        grad_v_flat: Tensor = grad_v.reshape(single_rep_target_shape)

        grad_z_flat: Tensor = grad_z.permute((0, 2, 3, 1))

        grad_q_dtensor = DTensor.from_local(
            grad_q_flat,
            device_mesh=ctx.device_mesh,
            placements=ctx.single_rep_placements,
            shape=ctx.shape_q,
            stride=ctx.stride_q,
        )
        grad_k_dtensor = DTensor.from_local(
            grad_k_flat,
            device_mesh=ctx.device_mesh,
            placements=ctx.single_rep_placements,
            shape=ctx.shape_k,
            stride=ctx.stride_k,
        )
        grad_v_dtensor = DTensor.from_local(
            grad_v_flat,
            device_mesh=ctx.device_mesh,
            placements=ctx.single_rep_placements,
            shape=ctx.shape_v,
            stride=ctx.stride_v,
        )
        grad_z_dtensor = DTensor.from_local(
            grad_z_flat,
            device_mesh=ctx.device_mesh,
            placements=ctx.pair_rep_placements,
            shape=ctx.shape_z,
            stride=ctx.stride_z,
        )
        _AttentionPairBiasContexVecImpl.check_backward_output_metadata(
            ctx,
            grad_q_dtensor,
            grad_k_dtensor,
            grad_v_dtensor,
        )
        return grad_q_dtensor, grad_k_dtensor, grad_v_dtensor, grad_z_dtensor, None, None, None

    @staticmethod
    def check_forward_input_metadata_and_store(
        ctx: FunctionCtx,
        q: DTensor,
        k: DTensor,
        v: DTensor,
        z: DTensor,
        mask: DTensor,
        pair_mask: DTensor | None,
        apb_context_vec_params: _AttentionPairBiasContextVecParams,
    ) -> None:
        (
            _,
            _,
            num_heads,
            head_dim,
            _,
            _,
            _,
        ) = apb_context_vec_params

        ctx.H = num_heads
        ctx.head_dim = head_dim

        ctx.B = q.shape[0]
        ctx.N = q.shape[1]
        ctx.c_s = q.shape[-1]

        placements_single_expected = (Shard(0), Shard(1), Replicate())
        placements_pair_expected = (Shard(0), Shard(1), Shard(2))
        ctx.single_rep_placements = q.placements
        ctx.pair_rep_placements = z.placements
        ctx.device_mesh = q.device_mesh
        ctx.shape_q = q.shape
        ctx.stride_q = q.stride()
        ctx.shape_k = k.shape
        ctx.stride_k = k.stride()
        ctx.shape_v = v.shape
        ctx.stride_v = v.stride()
        ctx.shape_z = z.shape
        ctx.stride_z = z.stride()

        check_metadata = raise_if_incorrect_dtensor_metadata_args

        check_metadata(
            q,
            "q",
            check_for_partial_placements=True,
            expected_placements=placements_single_expected,
        )
        check_metadata(
            k,
            "k",
            (ctx.B, ctx.N, ctx.c_s),
            expected_device_mesh=ctx.device_mesh,
            expected_placements=ctx.single_rep_placements,
        )
        check_metadata(
            v,
            "v",
            (ctx.B, ctx.N, ctx.c_s),
            expected_device_mesh=ctx.device_mesh,
            expected_placements=ctx.single_rep_placements,
        )
        check_metadata(
            z,
            "z",
            None,  # shape can be different from single representation(s) due to multiplicity
            expected_device_mesh=ctx.device_mesh,
            check_for_partial_placements=True,
            expected_placements=placements_pair_expected,
        )
        check_metadata(
            mask,
            "mask",
            None,
            expected_device_mesh=ctx.device_mesh,
            expected_placements=ctx.single_rep_placements,
        )
        if pair_mask is not None:
            check_metadata(
                pair_mask,
                "pair_mask",
                None,  # shape can be different from single representation(s) due to multiplicity
                expected_device_mesh=ctx.device_mesh,
                expected_placements=ctx.pair_rep_placements,
            )

    @staticmethod
    def check_backward_input_metadata(ctx: FunctionCtx, grad_output: DTensor) -> None:
        raise_if_incorrect_dtensor_metadata_args(
            grad_output,
            "grad_output",
            expected_shape=(ctx.B, ctx.N, ctx.c_s),
            expected_device_mesh=ctx.device_mesh,
            expected_placements=ctx.single_rep_placements,
        )

    @staticmethod
    def check_backward_output_metadata(
        ctx: FunctionCtx,
        grad_q_dtensor: DTensor,
        grad_k_dtensor: DTensor,
        grad_v_dtensor: DTensor,
    ) -> None:
        """DTensor.from_local(..) requires the specification of stride if
        shape is specified, so the specification of stride is side-stepped
        in this usage by checking the shape determined by DTensor library."""
        metadata_tuple = (
            (grad_q_dtensor, "grad_q", (ctx.B, ctx.N, ctx.c_s)),
            (grad_k_dtensor, "grad_k", (ctx.B, ctx.N, ctx.c_s)),
            (grad_v_dtensor, "grad_v", (ctx.B, ctx.N, ctx.c_s)),
        )
        for dtensor_instance, dtensor_name, expected_shape in metadata_tuple:
            if not dtensor_instance.shape == expected_shape:
                raise ValueError(
                    ", ".join(
                        [
                            f"dtensor '{dtensor_name}' should have shape {expected_shape}",
                            f"but instead has shape {dtensor_instance.shape}.",
                        ]
                    )
                )

    @staticmethod
    def check_forward_input_for_impl_state(apb_context_vec_params: _AttentionPairBiasContextVecParams) -> None:
        """
        Check the implementation scope of the forward pass.

        Parameters
        ----------
        apb_context_vec_params: AttentionPairBiasContextVecParams
            The parameters for the attention pair bias context vector.

        Returns
        -------
        None

        Raises
        ------
        NotImplementedError
        """
        (
            _,
            _,
            _,
            _,
            _,
            use_window_batching,
            _,
        ) = apb_context_vec_params

        if use_window_batching:
            raise NotImplementedError(f"use_window_batching={use_window_batching} is not implemented")


def ring_attention(
    q: Tensor,
    k: Tensor,
    v: Tensor,
    z: Tensor,
    pair_mask: Tensor,
    mask: Tensor,
    ring_comm: AttentionPairBiasComm,
    inf: float = 1e6,
) -> Tensor:
    """Functional interface to RingAttention autograd function.

    Based on vanilla torch tensors.
    """
    return RingAttention.apply(q, k, v, z, pair_mask, mask, ring_comm, inf)


class RingAttention(torch.autograd.Function):
    @staticmethod
    def forward(
        ctx,
        q: Tensor,
        k: Tensor,
        v: Tensor,
        z: Tensor,
        pair_mask: Tensor,
        mask: Tensor,
        ring_comm: AttentionPairBiasComm,
        inf: float = 1e6,
    ) -> Tensor:
        pair_mask = pair_mask[:, None, :, :] * mask[:, None, None, :]
        o, data_for_backward = ring_attention_simple_forward(
            q=q,
            k=k,
            v=v,
            z=z,
            pair_mask=pair_mask,
            ring_comm=ring_comm,
            inf=inf,
        )
        # Unpack tensors for save_for_backward to enable automatic memory management and hook support
        ctx.save_for_backward(
            data_for_backward.q_store,
            data_for_backward.k_t_store,
            data_for_backward.v_t_store,
            data_for_backward.z_store,
            data_for_backward.lse_m,
            data_for_backward.o_store,
        )
        ctx.ring_comm = data_for_backward.ring_comm
        ctx.sdpa_with_bias_backend = data_for_backward.sdpa_with_bias_backend
        return o

    @staticmethod
    def backward(ctx, do: Tensor) -> tuple[Tensor, Tensor, Tensor, Tensor, None, None, None, None]:
        q_store, k_t_store, v_t_store, z_store, lse_m, o_store = ctx.saved_tensors
        data_for_backward = RingAttentionSimpleDataForBackward(
            q_store=q_store,
            k_t_store=k_t_store,
            v_t_store=v_t_store,
            z_store=z_store,
            lse_m=lse_m,
            ring_comm=ctx.ring_comm,
            sdpa_with_bias_backend=ctx.sdpa_with_bias_backend,
            o_store=o_store,
        )
        del ctx.ring_comm
        del ctx.sdpa_with_bias_backend
        dq, dk, dv, dz = ring_attention_simple_backward(data_for_backward=data_for_backward, do=do)
        return dq, dk, dv, dz, None, None, None, None


class RingAttentionSimpleDataForBackward(NamedTuple):
    """Data for backward pass of ring attention."""

    q_store: Tensor
    k_t_store: Tensor
    v_t_store: Tensor
    z_store: Tensor
    lse_m: Tensor
    ring_comm: AttentionPairBiasComm
    sdpa_with_bias_backend: SDPAWithBiasBackend = SDPAWithBiasBackend.REFERENCE
    o_store: Tensor | None = None
    multiplicity: int = 1


def is_power_of_2(n: int) -> bool:
    """Check if n is a power of 2."""
    return (n > 0) and (n & (n - 1) == 0)


def ring_attention_simple_forward(
    q: Tensor,
    k: Tensor,
    v: Tensor,
    z: Tensor,
    pair_mask: Tensor,
    ring_comm: AttentionPairBiasComm,
    inf: float = 1e6,
    sdpa_with_bias_backend: SDPAWithBiasBackend = SDPAWithBiasBackend.TORCH_FLEX_ATTN,
) -> tuple[Tensor, RingAttentionSimpleDataForBackward]:
    """Forward pass of ring attention.

    example sharding strategy on N_tokens=2, world_size=4
    device mesh = [[  0,   1],
                   [  2,   3]]
    q/k/v, z, m = [[ q0,  q0],   [[ k0,  k0],   [[z00, z01],   [[m00, m01],
                   [ q1,  q1]]    [ k1,  k1]]     z10, z11]]     m10, m11]]

    step = 0
    k_t, z      = [[ k0,  k1],   [[z00, z01],
                   [ k0,  k1]]     z10, z11]]

    step = 1 (roll to the left)
    k_t, z      = [[ k1, k0],   [[z01, z00],
                   [ k1, k0]]     z11, z10]]

    now we roll left k, v, z as the ring attention outer loop.

    o           = [[ (q0k0+z00+m00)v0 + (q0k1+z01+m01)v1, (q0k1+z01+m01)v1 + (q0k0+z00+m00)v0 ],
                   [ (q1k0+z10+m10)v0 + (q1k1+z11+m11)v1, (q1k1+z11+m11)v1 + (q1k0+z10+m10)v0 ]]

    Parameters
    ----------
    q : Tensor
        Query tensor (B, I, H, D)
    k : Tensor
        Key tensor (B, J, H, D)
    v : Tensor
        Value tensor (B, J, H, D)
    z : Tensor
        Pair bias projection (B, H, I, J)
    pair_mask : Tensor
        Pair bias mask (B, 1, I, J)
    ring_comm : AttentionPairBiasComm
        Ring communication for async operation
    inf : float
        Infinity value for masking, by default 1e6

    Returns
    -------
    Tensor
        Output tensor (B, I, H, D)
    RingAttentionSimpleDataForBackward
        Data for backward pass
    """
    if sdpa_with_bias_backend == SDPAWithBiasBackend.TORCH_SDPA_EFFICIENT_ATTENTION:
        warnings.warn("torch_sdpa backend is not implemented and will fall back to flex_attention backend")
        sdpa_with_bias_backend = SDPAWithBiasBackend.TORCH_FLEX_ATTN

    # flex_attention (compiled Triton/Inductor) requires power-of-2 head_dim and head_dim >= 16.
    use_flex_attn = (
        sdpa_with_bias_backend == SDPAWithBiasBackend.TORCH_FLEX_ATTN
        and HAS_FLEX_ATTN
        and q.is_cuda
        and q.dtype != torch.float64
        and is_power_of_2(q.size(-1))
        and q.size(-1) >= 16  # head_dim >= 16
    )

    B_q, S, H, D = q.shape
    B_z_orig = z.shape[0]
    multiplicity = B_q // B_z_orig

    embed_dim = q.size(-1)

    # save input for backward
    requires_grad = q.requires_grad or k.requires_grad or v.requires_grad or z.requires_grad
    if requires_grad:
        q_store = q.detach()
    else:
        q_store = None

    # Overlap k, v comm with mask addition
    k_recv = ring_comm.comm_transpose_k.enqueue_to_dispatch(k.contiguous())
    v_recv = ring_comm.comm_transpose_v.enqueue_to_dispatch(v.contiguous())
    z = z.contiguous() + (1 - pair_mask) * -inf
    ring_comm.comm_transpose_k.wait_until_finished()
    ring_comm.comm_transpose_v.wait_until_finished()
    k = k_recv
    v = v_recv

    if requires_grad:
        k_t_store = k.detach()
        v_t_store = v.detach()
    else:
        k_t_store, v_t_store = None, None

    # save input for backward
    if requires_grad:
        z_store = z.detach()
    else:
        z_store = None

    # Ring attention
    o: Union[Tensor, None] = None
    lse_m: Union[Tensor, None] = None
    amax: Union[Tensor, None] = None

    size_device_grid0 = ring_comm.group_layout.shape[0]
    for step in range(size_device_grid0):
        # Overlap k, v, z+m roll-left comm
        if step + 1 != size_device_grid0:
            next_k = ring_comm.comm_k.enqueue_to_dispatch(k)
            next_v = ring_comm.comm_v.enqueue_to_dispatch(v)
            next_z = ring_comm.comm_z.enqueue_to_dispatch(z)

        # Attention by chunk
        if use_flex_attn:
            # Flex attention expects (B, H, S, D)
            q_f = q.transpose(1, 2)
            k_f = k.transpose(1, 2)
            v_f = v.transpose(1, 2)

            if multiplicity > 1:
                # flex_attention doesn't support broadcasting batch dim for bias yet
                # so we use a trick: we reshape B_q * H into the head dimension
                # and use (h // H) // multiplicity to index into z
                q_f = q_f.reshape(1, B_q * H, S, D)
                k_f = k_f.reshape(1, B_q * H, S, D)
                v_f = v_f.reshape(1, B_q * H, S, D)

                def score_mod(score, b, h, q_idx, kv_idx):
                    return score + z[(h // H) // multiplicity, h % H, q_idx, kv_idx]

            else:
                # B_q == B_z
                def score_mod(score, b, h, q_idx, kv_idx):
                    return score + z[b, h, q_idx, kv_idx]

            # flex_attention_compiled with float32 requires full precision
            with setup_tf32_env(Precision.FP32), torch.amp.autocast("cuda", enabled=False):
                block_o, aux_data = flex_attention_compiled(
                    q_f,
                    k_f,
                    v_f,
                    score_mod=score_mod,
                    return_lse=True,
                )

            # block_o: (B_q, H, S, D) or (1, B*H, S, D)
            # aux_data: (B_q, H, S) or (1, B*H, S)
            block_o = block_o.reshape(B_q, H, S, D).transpose(1, 2)
            block_lse_m = aux_data.reshape(B_q, H, S).transpose(1, 2).unsqueeze(-1)
            block_amax = None  # flex_attention doesn't return amax, but we can use lse
        else:
            attn = torch.einsum("bihd,bjhd->bhij", q, k)
            attn = attn / (embed_dim**0.5)

            B_z = z.shape[0]
            B_q = q.shape[0]
            if B_q != B_z:
                attn = (attn.view(B_z, -1, *attn.shape[1:]) + z.unsqueeze(1)).view_as(attn)
            else:
                attn = attn + z

            block_o = torch.einsum("bhij,bjhd->bihd", torch.softmax(attn, dim=-1), v)
            block_amax = attn.amax(dim=-1, keepdim=True)
            block_lse_m = torch.logsumexp(attn - block_amax, dim=-1, keepdim=True)

            block_lse_m = block_lse_m.transpose(-2, -3)
            if block_amax is not None:
                block_amax = block_amax.transpose(-2, -3)

        o, lse_m, amax = tiled_softmax_attention_update(block_o, block_lse_m, block_amax, o, lse_m, amax)

        # Get input for next round
        if step + 1 != size_device_grid0:
            ring_comm.comm_k.wait_until_finished()
            ring_comm.comm_v.wait_until_finished()
            ring_comm.comm_z.wait_until_finished()
            k, v, z = next_k, next_v, next_z

    if requires_grad and not use_flex_attn:
        if multiplicity > 1:
            # Reduce amax across multiplicity
            amax_t = amax.transpose(-2, -3)
            amax_t_view = amax_t.view(B_z_orig, multiplicity, *amax_t.shape[1:])
            amax_reduced = amax_t_view.amax(dim=1)
            z_store -= amax_reduced

            # Adjust lse_m to be relative to amax_reduced
            # We use broadcasting via views to avoid memory-heavy repeat_interleave
            amax_reduced_t = amax_reduced.transpose(-2, -3)  # (B, I, H, 1)
            lse_m_view = lse_m.view(B_z_orig, multiplicity, *lse_m.shape[1:])
            amax_view = amax.view(B_z_orig, multiplicity, *amax.shape[1:])

            lse_m_view += amax_view - amax_reduced_t.unsqueeze(1)  # modify lse_m inplace via view
        else:
            z_store -= amax.transpose(-2, -3)

    data_for_backward = RingAttentionSimpleDataForBackward(
        q_store=q_store,
        k_t_store=k_t_store,
        v_t_store=v_t_store,
        z_store=z_store,
        lse_m=lse_m,
        ring_comm=ring_comm,
        sdpa_with_bias_backend=sdpa_with_bias_backend,
        o_store=o.detach() if requires_grad else None,
        multiplicity=multiplicity,
    )

    return o, data_for_backward


def ring_attention_simple_backward(
    data_for_backward: RingAttentionSimpleDataForBackward,
    do: Tensor,
) -> tuple[Tensor, Tensor, Tensor, Tensor]:
    """Backward pass of ring attention.

    example sharding strategy on N_tokens=2, world_size=4
    device mesh = [[  0,   1],
                   [  2,   3]]

    do          = [[ do0, do0],
                   [ do1, do1]]

    ctx should have saved:
    q    , z    = [[ q0,  q0],   [[z00, z01],
                   [ q1,  q1]]    [z10, z11]]
    k_t, v_t    = [[ k0,  k1],   [[v0,  v1],
                   [ k0,  k1]]    [v0,  v1]]
    lse         = [[ lse0, lse0],
                   [ lse1, lse1]]

    where z has already accounted for the pair mask

    The output gradient should be sharded as follows:
    dq          = [[ dq0, dq0],
                   [ dq1, dq1]]
    dk          = [[ dk0, dk0],
                   [ dk1, dk1]]
    dv          = [[ dv0, dv0],
                   [ dv1, dv1]]
    dz          = [[ dz00, dz01],
                   [ dz10, dz11]]

    Parameters
    ----------
    data_for_backward : RingAttentionSimpleDataForBackward
        Data for backward pass
    do : Tensor
        Gradient of output tensor

    Returns
    -------
    tuple[Tensor, Tensor, Tensor, Tensor]
        Gradients for q, k, v, z
    """
    q, k_t, v_t, z, lse_m, ring_comm, sdpa_with_bias_backend, out_global, multiplicity = data_for_backward
    if sdpa_with_bias_backend == SDPAWithBiasBackend.TORCH_SDPA_EFFICIENT_ATTENTION:
        raise NotImplementedError("torch_sdpa backend is not implemented")

    embed_dim = q.size(-1)

    # Pre-calculate out_global * do for softmax gradient
    # do: (B, S, H, D)
    # out_global: (B, S, H, D)
    # do_o: (B, S, H, 1) -> will be reshaped to (B, H, S, 1) for the reference path
    do_o = torch.sum(do * out_global, dim=-1, keepdim=True)

    lse_m_t = lse_m.transpose(1, 2)  # (B, S, H, 1) -> (B, H, S, 1)

    B_q = q.shape[0]

    # flex_attention (compiled Triton/Inductor) requires power-of-2 head_dim and head_dim >= 16.
    use_flex_attn = (
        sdpa_with_bias_backend == SDPAWithBiasBackend.TORCH_FLEX_ATTN
        and HAS_FLEX_ATTN
        and q.is_cuda
        and q.dtype != torch.float64
        and is_power_of_2(q.size(-1))
        and q.size(-1) >= 16  # head_dim >= 16
    )

    if use_flex_attn:
        # Re-run forward to get aux data
        # In a real implementation we would save this
        # but for Ring Attention we might need to re-run or save it per chunk
        B_q, S, H, D = q.shape

        with torch.enable_grad():
            q_l = q.detach().requires_grad_(True)
            k_l = k_t.detach().requires_grad_(True)
            v_l = v_t.detach().requires_grad_(True)
            z_l = z.detach().requires_grad_(True)

            q_f = q_l.transpose(1, 2)
            k_f = k_l.transpose(1, 2)
            v_f = v_l.transpose(1, 2)

            if multiplicity > 1:
                q_f = q_f.reshape(1, B_q * H, S, D)
                k_f = k_f.reshape(1, B_q * H, S, D)
                v_f = v_f.reshape(1, B_q * H, S, D)

                def score_mod(score, b, h, q_idx, kv_idx):
                    return score + z_l[(h // H) // multiplicity, h % H, q_idx, kv_idx]

            else:
                # B_q == B_z
                def score_mod(score, b, h, q_idx, kv_idx):
                    return score + z_l[b, h, q_idx, kv_idx]

            # flex_attention_compiled with float32 requires full precision
            # Request LSE to perform the global scaling trick
            # This ensures the local gradients match the global context
            with setup_tf32_env(Precision.FP32), torch.amp.autocast("cuda", enabled=False):
                out_l, aux_l = flex_attention_compiled(
                    q_f,
                    k_f,
                    v_f,
                    score_mod=score_mod,
                    return_lse=True,
                )

                out_l = out_l.reshape(B_q, H, S, D).transpose(1, 2)
                lse_l = aux_l.reshape(B_q, H, S, 1).transpose(1, 2)

                size_device_grid0 = ring_comm.group_layout.shape[0]
                if size_device_grid0 == 1:
                    grad_inputs = torch.autograd.grad(out_l, (q_l, k_l, v_l, z_l), do, allow_unused=True)
                else:
                    # Scaling trick: Global contribution of this chunk
                    # Contribution_global = (out_local - out_global) * exp(lse_local - lse_global)
                    # d(Contribution_global)/dx gives the correct softmax gradient including denominator
                    out_scaled = (out_l - out_global) * torch.exp(lse_l - lse_m)
                    grad_inputs = torch.autograd.grad(out_scaled, (q_l, k_l, v_l, z_l), do, allow_unused=True)

            dq, dk, dv, dz = grad_inputs
            dq = dq if dq is not None else torch.zeros_like(q_l)
            dk = dk if dk is not None else torch.zeros_like(k_l)
            dv = dv if dv is not None else torch.zeros_like(v_l)
            dz = dz if dz is not None else torch.zeros_like(z_l)

            if multiplicity > 1 and dz.shape[0] > (q.shape[0] // multiplicity):
                B_z_orig = q.shape[0] // multiplicity
                dz = dz.view(B_z_orig, multiplicity, *dz.shape[1:]).sum(dim=1)
    else:
        # Compute S_ij and A_ij
        s = torch.einsum("bihd,bjhd->bhij", q, k_t)
        s /= embed_dim**0.5

        # Memory efficient in-place softmax reconstruction
        s.sub_(lse_m_t)
        if multiplicity > 1:
            B_z_orig = q.shape[0] // multiplicity
            if z.shape[0] == q.shape[0]:
                s.view(B_z_orig, multiplicity, *s.shape[1:]).add_(z.view(B_z_orig, multiplicity, *z.shape[1:]))
            else:
                s.view(B_z_orig, multiplicity, *s.shape[1:]).add_(z.unsqueeze(1))
        else:
            s.add_(z)

        a = s.exp_()

        # Compute gradient of v and c
        # dV_j = \Sum_{i} A^T_{ij} dO_i
        # c_i = \Sum_{k} v^T_k A_{ik}
        dv = torch.einsum("bihd,bhij->bjhd", do, a)

        # Compute gradient of S (dS = dz)
        # dS_{ij} = A_{ij} * (dO_i v^T_j - (dO_i * out_global_i))
        # do_o has shape (B, S, H, 1), we need (B, H, S, 1) to match a
        do_o_step = do_o.transpose(1, 2)

        # In-place compute dS to save memory
        tmp = torch.einsum("bihd,bjhd->bhij", do, v_t)
        tmp.sub_(do_o_step)
        a.mul_(tmp)
        del tmp

        dS = a
        dz = a
        if multiplicity > 1 and dz.shape[0] > (q.shape[0] // multiplicity):
            B_z_orig = q.shape[0] // multiplicity
            dz = dz.view(B_z_orig, multiplicity, *dz.shape[1:]).sum(dim=1)

        # Compute gradient of q, k
        # dq_i = \Sum_{j} dS_{ij} k^T_j
        # dk_j = \Sum_{i} dS^T_{ij} q_i
        dq = torch.einsum("bhij,bjhd->bihd", dS, k_t) / (
            embed_dim**0.5
        )  # _t here refers to transposition across device mesh
        dk = torch.einsum("bhij,bihd->bjhd", dS, q) / (embed_dim**0.5)

    dv_recv = ring_comm.comm_transpose_v.enqueue_to_dispatch(dv.contiguous())

    # collect and complete dv reduction
    ring_comm.comm_transpose_v.wait_until_finished()
    dv_work = dist.all_reduce(dv_recv, op=dist.ReduceOp.SUM, group=ring_comm.cp_axis_1_group, async_op=True)

    dq = dq.contiguous()
    dq_work = dist.all_reduce(dq, op=dist.ReduceOp.SUM, group=ring_comm.cp_axis_1_group, async_op=True)

    dk_recv = ring_comm.comm_transpose_k.enqueue_to_dispatch(dk.contiguous())
    ring_comm.comm_transpose_k.wait_until_finished()
    dk_work = dist.all_reduce(dk_recv, op=dist.ReduceOp.SUM, group=ring_comm.cp_axis_1_group, async_op=True)

    # Collect all async works
    dq_work.wait()
    dk_work.wait()
    dv_work.wait()

    return dq, dk_recv, dv_recv, dz


class RingAttentionSimple:
    """Ring attention pair bias with context parallelism.

    This class serves as a namespace for the ring attention forward and
    backward functions.  The definitions of these functions outside an
    autograd function subclass, are useful to encapsulate the communication
    and math logic that is used in _AttentionPairBiasContexVecImpl, but
    also may be used elsewhere in the codebase via the simple function
    ring_attention above.
    """

    @staticmethod
    def forward(
        q: Tensor,
        k: Tensor,
        v: Tensor,
        z: Tensor,
        pair_mask: Tensor,
        mask: Tensor,
        ring_comm: AttentionPairBiasComm,
        inf: float = 1e6,
    ) -> tuple[Tensor, RingAttentionSimpleDataForBackward]:
        assert isinstance(q, Tensor), f"q must be a Tensor, got {type(q)}"
        assert isinstance(k, Tensor), f"k must be a Tensor, got {type(k)}"
        assert isinstance(v, Tensor), f"v must be a Tensor, got {type(v)}"
        assert isinstance(z, Tensor), f"z must be a Tensor, got {type(z)}"
        assert isinstance(pair_mask, Tensor), f"pair_mask must be a Tensor, got {type(pair_mask)}"
        assert isinstance(mask, Tensor), f"mask must be a Tensor, got {type(mask)}"
        assert isinstance(
            ring_comm, AttentionPairBiasComm
        ), f"ring_comm must be a AttentionPairBiasComm, got {type(ring_comm)}"
        return ring_attention_simple_forward(q, k, v, z, pair_mask, mask, ring_comm, inf)

    @staticmethod
    def backward(
        data_for_backward: RingAttentionSimpleDataForBackward,
        do: Tensor,
    ) -> tuple[Tensor, Tensor, Tensor, Tensor]:
        return ring_attention_simple_backward(data_for_backward, do)
