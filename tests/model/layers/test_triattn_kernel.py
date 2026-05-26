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

import cuequivariance_torch.primitives.triangle as cueq_triangle
import pytest
import torch
from trifast.torch import _triangle_attention as trifast_triangle_attention
from trifast.torch import triangle_attention_bwd as trifast_triangle_attention_bwd

from boltz.distributed.model.modules.utils import PRECISION_TO_DTYPE, Precision, TriAttnBackend, setup_tf32_env
from boltz.model.layers.triangular_attention.primitives import _attention, cueq_is_installed, trifast_is_installed
from boltz.testing.utils import (
    PRECISION_TO_INF,
    assert_no_percentile_upshift,
    init_tensors_normal,
    init_tensors_uniform,
)


def view_tensor_strided(t: torch.Tensor, argsort_strides: torch.Tensor):
    """
    Create a strided view of a tensor with custom stride ordering.

    This function reorders the strides of a tensor based on the provided argsort
    indices, allowing for custom memory layout configurations without copying data.

    Args:
        t: Input tensor to create a strided view of.
        argsort_strides: A 1D tensor containing indices that specify the desired
            stride ordering. Must have the same number of elements as the number
            of dimensions in `t`.

    Returns:
        A strided view of the input tensor with reordered strides according to
        the specified ordering.

    Raises:
        AssertionError: If the number of dimensions in `t` doesn't match the
            number of elements in `argsort_strides`.
    """
    assert t.ndim == argsort_strides.numel(), "shape and argsort_strides must have the same length"
    shape_sorted_by_strides = torch.tensor(t.shape)[argsort_strides[:-1]]
    strides_sorted = torch.tensor([1] + shape_sorted_by_strides.tolist()).cumprod(dim=0)
    # inverse argsort_strides
    argsort_strides_inv = torch.argsort(argsort_strides)
    strides_new = strides_sorted[argsort_strides_inv]
    ans = torch.as_strided(t.flatten(), t.shape, strides_new.tolist())
    return ans


def make_input(
    B,
    H,
    Q_len,  # Query/output sequence length
    KV_len,  # Key/Value sequence length
    C_hidden,
    device,
    dtype,
    seed: int,
    min_max: tuple[float, float] | None = None,
    inf: float = 1e9,
    use_mask: bool = True,
    argsort_strides: torch.Tensor | None = None,
):
    """Create test input tensors for triangular attention.

    Args:
        B: Batch size
        H: Number of heads
        Q_len: Query sequence length (I in the original formulation)
        KV_len: Key/Value sequence length (J in the original formulation)
        C_hidden: Hidden dimension
        device: Device
        dtype: Data type
        seed: Random seed
        min_max: If provided, randomly sample from [min, max] uniformly for all
            the returned tensors; otherwise, sample from standard normal distribution
        inf: Infinity value for mask
        use_mask: Whether to use mask
        argsort_strides: If provided, reorder the strides of the returned tensors
            according to the provided indices. In effect, the returned
            tensors' stride()[argsort_strides] will be a exclusive cumulative product of
            tensors' shape[argsort_strides].
    """
    N = Q_len
    J = KV_len
    Q = Q_len
    K = KV_len
    V = KV_len
    torch.manual_seed(seed)

    # Create mask
    if use_mask:
        mask = torch.randint(0, 2, (B, N, 1, 1, J), device=device, dtype=dtype, requires_grad=False)
        # Set some regions to zero for testing masking
        if mask.shape[1] > 1:
            mask[0, mask.shape[1] // 2 :, :, :, :] = 0
        if mask.shape[-1] > 1:
            mask[0, :, :, :, mask.shape[-1] // 2 :] = 0
    else:
        mask = None

    q = torch.empty(B, N, H, Q, C_hidden, device=device, dtype=dtype, requires_grad=True)
    k = torch.empty(B, N, H, K, C_hidden, device=device, dtype=dtype, requires_grad=True)
    v = torch.empty(B, N, H, V, C_hidden, device=device, dtype=dtype, requires_grad=True)
    triangle_bias = torch.empty(B, 1, H, N, J, device=device, dtype=dtype, requires_grad=True)

    do = torch.empty_like(q)

    if min_max is None:
        init_tensors_normal([q, k, v, triangle_bias, do])
    else:
        init_tensors_uniform([q, k, v, triangle_bias, do], min_max[0], min_max[1])

    # zero-initialize do for the invalid elements
    if mask is not None:
        with torch.no_grad():
            do = do * mask.any(dim=-1, keepdim=True)

    if argsort_strides is not None:
        assert (
            argsort_strides.ndim == 1 and argsort_strides.numel() == 5
        ), "argsort_strides must be a 1D tensor of length 5"
        q = view_tensor_strided(q, argsort_strides)
        k = view_tensor_strided(k, argsort_strides)
        v = view_tensor_strided(v, argsort_strides)
        mask = view_tensor_strided(mask, argsort_strides)
        triangle_bias = view_tensor_strided(triangle_bias, argsort_strides)
        do = view_tensor_strided(do, argsort_strides)

    return q, k, v, mask, triangle_bias, do


def run_triangle_attention(
    q,
    k,
    v,
    triangle_bias,
    do,
    mask: torch.Tensor | None = None,
    backend: TriAttnBackend = TriAttnBackend.REFERENCE,
    precision: Precision = Precision.FP32,
    check_bwd: bool = True,
    scale: float = 1.0,
    dtype_triangle_bias: torch.dtype | None = None,
):
    """Run triangle attention operation with specified backend and precision.

    This function executes the triangle attention mechanism, which applies attention
    with an additional triangle bias term. It supports a PyTorch reference implementation
    and optimized GPU kernels (CUEQ, TRIFAST), and can run in various precisions.

    Args:
        q: Query tensor of shape [batch, num_heads, seq_len, head_dim].
        k: Key tensor of shape [batch, num_heads, seq_len, head_dim].
        v: Value tensor of shape [batch, num_heads, seq_len, head_dim].
        triangle_bias: Triangle bias tensor of shape [batch, num_heads, seq_len, seq_len].
        do: Gradient output tensor of shape [batch, num_heads, seq_len, head_dim],
            used for backward pass.
        mask: Optional boolean mask tensor of shape [batch, 1, seq_len, seq_len] or
            broadcastable shape. If provided, attention is masked where mask is False.
            Defaults to None.
        backend: Backend implementation to use (TriAttnBackend.REFERENCE, TriAttnBackend.CUEQ,
            or TriAttnBackend.TRIFAST). Defaults to TriAttnBackend.REFERENCE.
        precision: Precision mode for computation (Precision.FP16, Precision.BF16,
            Precision.TF32, Precision.FP32, or Precision.FP64). FP64 is only supported
            with the reference backend. Defaults to Precision.FP32.
        check_bwd: Whether to run the backward pass and compute gradients.
            Defaults to True.
        scale: Scaling factor applied to attention scores. Defaults to 1.0.
        dtype_triangle_bias: Optional dtype override for triangle_bias. If None,
            uses the target dtype from precision. Defaults to None.

    Returns:
        tuple: A 4-tuple containing:
            - output: Attention output tensor of shape [batch, num_heads, seq_len, head_dim].
            - lse_m: Log-sum-exp values (shifted by max) of shape
                [batch, num_heads, seq_len, 1].
            - amax: Maximum attention scores of shape [batch, num_heads, seq_len, 1].
            - input_grads: Dictionary with gradients for 'q', 'k', 'v', and 'triangle_bias'.
                If check_bwd is False, all gradients are None.

    Raises:
        ValueError: If precision is FP64 and backend is not REFERENCE.
        ValueError: If backend is unknown.

    Note:
        - Input tensors are cloned and converted to the target precision internally.
        - The function uses appropriate TF32 environment settings based on precision (except TRIFAST).
        - For CUEQ and TRIFAST backends, LSE/amax values are reshaped/adjusted to match reference format.
    """
    if precision == Precision.FP64 and backend != TriAttnBackend.REFERENCE:
        raise ValueError("FP64 is only supported for reference backend")

    if precision == Precision.TF32 and backend == TriAttnBackend.TRIFAST:
        # trifast hardcodes the "input_precision" to "ieee", i.e., FP32,
        # in tl.dot call, which locks down to doing FP32 matmul when input dtype is FP32.
        raise ValueError("TF32 is not supported for TRIFAST backend")

    device = q.device

    target_dtype = PRECISION_TO_DTYPE[precision]

    # Convert inputs to target precision
    q_work = q.detach().clone().to(dtype=target_dtype, device=device).requires_grad_(True)
    k_work = k.detach().clone().to(dtype=target_dtype, device=device).requires_grad_(True)
    v_work = v.detach().clone().to(dtype=target_dtype, device=device).requires_grad_(True)
    triangle_bias_work = (
        triangle_bias.detach()
        .clone()
        .to(
            dtype=target_dtype if dtype_triangle_bias is None else dtype_triangle_bias,
            device=device,
        )
        .requires_grad_(True)
    )

    if mask is None:
        mask_work = None
    else:
        mask_work = mask.to(dtype=bool, device=device)

    do_work = do.detach().clone().to(dtype=target_dtype, device=device)

    # Must not change the input tensors' memory layout here
    assert q_work.stride() == q.stride(), "q_work.stride() must be the same as q.stride()"
    assert k_work.stride() == k.stride(), "k_work.stride() must be the same as k.stride()"
    assert v_work.stride() == v.stride(), "v_work.stride() must be the same as v.stride()"
    assert (
        triangle_bias_work.stride() == triangle_bias.stride()
    ), "triangle_bias_work.stride() must be the same as triangle_bias.stride()"
    if mask_work is not None:
        (
            mask_work.stride() == mask.stride(),
            "mask_work.stride() must be the same as mask.stride()",
        )
    assert do_work.stride() == do.stride(), "do_work.stride() must be the same as do.stride()"

    # Run forward pass based on backend
    if backend == TriAttnBackend.REFERENCE:
        # reference implementation uses mask bias instead of mask
        inf = PRECISION_TO_INF[precision]
        if mask_work is None:
            mask_bias = torch.zeros(
                (q_work.shape[0], q_work.shape[1], 1, 1, k_work.shape[3]), dtype=target_dtype, device=device
            )
        else:
            mask_bias = inf * (mask_work.to(dtype=target_dtype) - 1.0)
        biases = [mask_bias, triangle_bias_work]
        with setup_tf32_env(precision):
            output, lse_m, amax = _attention(q_work * scale, k_work, v_work, biases, return_lse=True)
            # Run backward pass if requested
            if check_bwd:
                output.backward(do_work)
        # Collect gradients
        input_grads = {
            "q": q_work.grad,
            "k": k_work.grad,
            "v": v_work.grad,
            "triangle_bias": triangle_bias_work.grad,
        }
    elif backend == TriAttnBackend.CUEQ_FWD_TRIFAST_BWD:
        # Forward pass uses CUEQ, backward pass uses TRIFAST
        with setup_tf32_env(precision):
            with torch.no_grad():
                # Run CUEQ forward
                output, lse, amax = cueq_triangle.triangle_attention(
                    q_work,
                    k_work,
                    v_work,
                    triangle_bias_work,
                    mask=mask_work,
                    scale=scale,
                    return_aux=True,
                )
                # add back the singleton K axis resulting from the max reduction
                lse_reshaped = lse.unsqueeze(-1)
                # need to return amax unsqueezed for comparison with reference
                amax = amax.unsqueeze(-1)
                lse_m = lse_reshaped - amax

            # Run backward pass if requested
            if check_bwd:
                # Reshape tensors from CUEQ format to TRIFAST format for backward
                # q: [B, I, H, Q, C_hidden] --> [B, H, I, Q, C_hidden]
                q_trifast = q_work.detach().transpose(-3, -4).contiguous().requires_grad_(True)
                # k: [B, I, H, K, C_hidden] --> [B, H, I, K, C_hidden]
                k_trifast = k_work.detach().transpose(-3, -4).contiguous().requires_grad_(True)
                # v: [B, I, H, V, C_hidden] --> [B, H, I, V, C_hidden]
                v_trifast = v_work.detach().transpose(-3, -4).contiguous().requires_grad_(True)
                # triangle_bias: [B, 1, H, I, J] --> [B, H, I, J]
                triangle_bias_trifast = triangle_bias_work.detach().squeeze(-4).contiguous().requires_grad_(True)
                # output: [B, I, H, V, C_hidden] --> [B, H, I, V, C_hidden]
                o_trifast = output.detach().transpose(-3, -4).contiguous()
                # lse: [B, I, H, Q, 1] --> [B, H, I, Q, 1], and use lse instead of lse_m
                lse_trifast = (lse_m + amax).detach().transpose(-3, -4).contiguous()
                # do: [B, I, H, Q, C_hidden] --> [B, H, I, Q, C_hidden]
                do_trifast = do_work.detach().transpose(-3, -4).contiguous()

                # TRIFAST mask convention: True for invalid positions, False for valid positions
                # mask: [B, I, 1, 1, J] --> [B, I, J]
                mask_trifast = ~(mask_work.detach().squeeze((-2, -3)).contiguous())

                # Call TRIFAST backward
                dq_trifast, dk_trifast, dv_trifast, dtriangle_bias_trifast, _ = trifast_triangle_attention_bwd(
                    do_trifast,
                    q_trifast,
                    k_trifast,
                    v_trifast,
                    triangle_bias_trifast,
                    o_trifast,
                    lse_trifast.squeeze(-1).to(dtype=torch.float32),
                    mask_trifast,
                )

                # Reshape gradients back to CUEQ format
                # dq: [B, H, I, Q, C_hidden] --> [B, I, H, Q, C_hidden]
                dq_cueq = dq_trifast.transpose(-3, -4).contiguous()
                # dk: [B, H, I, K, C_hidden] --> [B, I, H, K, C_hidden]
                dk_cueq = dk_trifast.transpose(-3, -4).contiguous()
                # dv: [B, H, I, V, C_hidden] --> [B, I, H, V, C_hidden]
                dv_cueq = dv_trifast.transpose(-3, -4).contiguous()
                # dtriangle_bias: [B, H, I, J] --> [B, 1, H, I, J]
                dtriangle_bias_cueq = dtriangle_bias_trifast.unsqueeze(-4).contiguous()

                # Manually set gradients
                q_work.grad = dq_cueq
                k_work.grad = dk_cueq
                v_work.grad = dv_cueq
                triangle_bias_work.grad = dtriangle_bias_cueq

                input_grads = {
                    "q": q_work.grad,
                    "k": k_work.grad,
                    "v": v_work.grad,
                    "triangle_bias": triangle_bias_work.grad,
                }
            else:
                input_grads = {"q": None, "k": None, "v": None, "triangle_bias": None}
    elif backend == TriAttnBackend.CUEQ:
        with setup_tf32_env(precision):
            output, lse, amax = cueq_triangle.triangle_attention(
                q_work,
                k_work,
                v_work,
                triangle_bias_work,
                mask=mask_work,
                scale=scale,
                return_aux=True,
            )
            # add back the singleton K axis resulting from the max reduction
            lse_reshaped = lse.unsqueeze(-1)
            # need to return amax unsqueezed for comparison with reference
            amax = amax.unsqueeze(-1)
            lse_m = lse_reshaped - amax
            # manually call backward pass to emulate CP usage if requested
            if check_bwd:
                output.backward(do_work)
                # Collect gradients
                input_grads = {
                    "q": q_work.grad,
                    "k": k_work.grad,
                    "v": v_work.grad,
                    "triangle_bias": triangle_bias_work.grad,
                }
            else:
                input_grads = {"q": None, "k": None, "v": None, "triangle_bias": None}
    elif backend == TriAttnBackend.TRIFAST:
        # No need to setup TF32 environment for TRIFAST
        # as it hardcodes the "input_precision" to "ieee" (FP32) in tl.dot calls.
        q_trifast = q_work.transpose(1, 2)
        k_trifast = k_work.transpose(1, 2)
        v_trifast = v_work.transpose(1, 2)
        triangle_bias_trifast = triangle_bias_work.squeeze(1)
        # TRIFAST mask convention: True for invalid positions, False for valid positions
        if mask_work is None:
            mask_trifast = torch.zeros(
                q_trifast.shape[0], q_trifast.shape[2], q_trifast.shape[3], device=device, dtype=torch.bool
            )
        else:
            # mask: [B, I, 1, 1, K] --> [B, I, K]
            mask_trifast = ~(mask_work.squeeze(dim=(2, 3)).to(dtype=torch.bool))
        output_trifast, lse = trifast_triangle_attention(
            q_trifast, k_trifast, v_trifast, triangle_bias_trifast, mask_trifast
        )
        # output: [B, H, I, J, C_hidden] --> [B, I, H, J, C_hidden]
        output = output_trifast.transpose(1, 2)
        if check_bwd:
            output.backward(do_work)
            # Collect gradients
            input_grads = {
                "q": q_work.grad,
                "k": k_work.grad,
                "v": v_work.grad,
                "triangle_bias": triangle_bias_work.grad,
            }
        # TRIFAST's _triangle_attention API returns lse instead of lse - amax
        # (i.e., lse_m). We return lse_m as-is and set amax to None to indicate this difference.
        # lse: [B, H, I, J] --> [B, I, H, J, 1]
        lse_m = lse.transpose(1, 2).unsqueeze(-1)
        amax = None
    else:
        raise ValueError(f"Unknown backend: {backend}")

    return output, lse_m, amax, input_grads


@pytest.mark.skipif(not torch.cuda.is_available(), reason="CUDA is not available")
@pytest.mark.parametrize(
    "backend",
    [TriAttnBackend.CUEQ, TriAttnBackend.TRIFAST, TriAttnBackend.CUEQ_FWD_TRIFAST_BWD],
    ids=lambda x: f"backend:{x.value}",
)
@pytest.mark.parametrize("precision", [Precision.FP16, Precision.BF16, Precision.TF32, Precision.FP32])
@pytest.mark.parametrize(
    "argsort_strides",
    [
        None,  # LayoutRight
        torch.tensor([4, 2, 3, 1, 0]),  # D->H->S->I->B, commonly seen layout in OF/AF impl
        torch.tensor([4, 3, 2, 1, 0]),  # LayoutLeft
    ],
    ids=lambda x: f"argsort_strides:{'-'.join(map(str, x.tolist())) if x is not None else 'None'}",
)
def test_triangle_attention_kernel(backend, precision, argsort_strides):
    """Test triangle attention kernels against reference implementation with error histogram analysis.

    This test validates triangle attention backend implementations (CUEQ, TRIFAST) against
    the reference PyTorch implementation using error histogram analysis. It compares:
    1. FP64 reference (high precision baseline)
    2. Target precision reference (FP16/BF16/TF32/FP32)
    3. Backend kernel result at target precision (CUEQ or TRIFAST)

    The test also verifies:
    - Forward pass output accuracy
    - Backward pass gradient accuracy (when supported)
    - Memory layout independence (contiguous vs strided)
    - Masking correctness (zero gradients for masked regions)
    """

    if backend == TriAttnBackend.CUEQ and not cueq_is_installed:
        pytest.skip("cuequivariance_torch is not installed")

    if backend == TriAttnBackend.TRIFAST and not trifast_is_installed:
        pytest.skip("trifast is not installed")

    if backend == TriAttnBackend.CUEQ_FWD_TRIFAST_BWD:
        if not cueq_is_installed:
            pytest.skip("cuequivariance_torch is not installed")
        if not trifast_is_installed:
            pytest.skip("trifast is not installed")
        if precision != Precision.FP32:
            pytest.skip("CUEQ_FWD_TRIFAST_BWD only supports FP32 precision")

    if precision == Precision.TF32 and backend in (TriAttnBackend.TRIFAST, TriAttnBackend.CUEQ_FWD_TRIFAST_BWD):
        # TRIFAST hardcodes the "input_precision" to "ieee" (FP32) in tl.dot calls,
        # which locks down to FP32 matmul when input dtype is FP32.
        pytest.skip("TF32 is not supported for TRIFAST backend")

    # Test parameters
    H = 4
    N = 64
    C_hidden = 32
    B = 3
    device = "cuda:0"
    scale = 1 / math.sqrt(C_hidden)

    # Skip backward pass for FP32 (without TF32) with CUEQ backend as it doesn't support it
    # TRIFAST and CUEQ_FWD_TRIFAST_BWD support backward pass for all precisions
    check_bwd = backend in (TriAttnBackend.TRIFAST, TriAttnBackend.CUEQ_FWD_TRIFAST_BWD) or (
        backend == TriAttnBackend.CUEQ
        and (precision == Precision.TF32 or precision == Precision.BF16 or precision == Precision.FP16)
    )

    if precision == Precision.TF32:
        if backend == TriAttnBackend.CUEQ:
            # NOTE: The numerical error from CUEQ triangle attention with TF32 is significant.
            # Use smaller input values to keep gradients in a reasonable range (1e-6 to 1e-5)
            # for accurate gradient checking.
            min_val = -0.02
            max_val = 0.02
        elif backend in (TriAttnBackend.TRIFAST, TriAttnBackend.CUEQ_FWD_TRIFAST_BWD):
            min_val = -0.05
            max_val = 0.05
    elif precision == Precision.BF16:
        min_val = -0.5
        max_val = 0.5
    elif precision == Precision.FP16:
        min_val = -0.5
        max_val = 0.5
    elif precision == Precision.FP32:
        min_val = -0.5
        max_val = 0.5
    else:
        raise ValueError(f"Unsupported precision: {precision}")

    # Create test inputs in FP64 for highest precision, then cast as needed in run_triangle_attention
    seed = 42
    q, k, v, mask, triangle_bias, do = make_input(
        B,
        H,
        N,
        N,
        C_hidden,
        device,
        torch.float64,
        seed,
        min_max=(min_val, max_val),
        inf=1e18,
        use_mask=True,
        argsort_strides=argsort_strides,
    )

    # === RUN COMPUTATIONS WITH DIFFERENT BACKENDS AND PRECISIONS ===

    # Run FP64 reference (high precision baseline)
    o_expected_fp64, lse_m_expected_fp64, amax_expected_fp64, grads_fp64 = run_triangle_attention(
        q,
        k,
        v,
        triangle_bias,
        do,
        mask=mask,
        backend=TriAttnBackend.REFERENCE,
        precision=Precision.FP64,
        check_bwd=check_bwd,
        scale=scale,
    )

    # Run FP32/TF32 reference (alternative precision)
    o_expected_alt, lse_m_expected_alt, amax_expected_alt, grads_alt = run_triangle_attention(
        q,
        k,
        v,
        triangle_bias,
        do,
        mask=mask,
        backend=TriAttnBackend.REFERENCE,
        precision=precision,
        check_bwd=check_bwd,
        scale=scale,
    )

    # Run test backend implementation
    o_result, lse_m_result, amax_result, grads_result = run_triangle_attention(
        q,
        k,
        v,
        triangle_bias,
        do,
        mask=mask,
        backend=backend,
        precision=precision,
        check_bwd=check_bwd,
        scale=scale,
    )

    # Compute masks for proper comparison
    if mask is not None:
        # mask: [B, I, 1, 1, K] --> [B, I, 1, K, 1] for dk and dv masking
        mask_kv = mask.to(dtype=mask.dtype)[:, :, :, 0, :, None]
        # mask_i: [B, I, 1, 1, 1] for output, lse, amax masking
        mask_i = mask.any(dim=-1, keepdim=True)
        if grads_result["triangle_bias"] is not None:
            # mask_j: [B, 1, 1, 1, K] for dtriangle_bias masking
            mask_j = mask.any(dim=1, keepdim=True).to(dtype=grads_result["triangle_bias"].dtype)
        else:
            mask_j = torch.ones((B, 1, 1, 1, N), dtype=o_result.dtype, requires_grad=False, device=o_result.device)
    else:
        mask_kv = torch.ones(
            (B, N, 1, N, 1),
            dtype=k.dtype,
            requires_grad=False,
            device=o_result.device,
        )
        mask_i = torch.ones(
            (B, N, 1, 1, 1),
            dtype=q.dtype,
            requires_grad=False,
            device=o_result.device,
        )
        mask_j = torch.ones(
            (B, 1, 1, 1, N),
            dtype=triangle_bias.dtype,
            requires_grad=False,
            device=o_result.device,
        )

    if precision != Precision.FP32:
        # === ERROR HISTOGRAM ANALYSIS ===
        # Compare kernel implementations (CUEQ/TRIFAST) against PyTorch reference
        # Different backends use different algorithmic approaches, so precision characteristics will differ

        # use mask to select valid elements for comparison
        # to avoid numerical noise from invalid elements
        mask_i = mask_i.to(dtype=bool)
        # convert to FP32 for difference calculation
        assert_no_percentile_upshift(
            o_result[mask_i.expand_as(o_result)],
            o_expected_fp64[mask_i.expand_as(o_expected_fp64)].to(dtype=torch.float32),
            o_expected_alt[mask_i.expand_as(o_expected_alt)],
            names_input=(f"o_{backend.value}_{precision}", "o_ref_fp64", f"o_ref_{precision}"),
        )

        # Test for lse_m and amax (backend-specific due to different return conventions)
        if backend == TriAttnBackend.CUEQ:
            assert_no_percentile_upshift(
                lse_m_result[mask_i.expand_as(lse_m_result)],
                lse_m_expected_fp64[mask_i.expand_as(lse_m_expected_fp64)].to(dtype=torch.float32),
                lse_m_expected_alt[mask_i.expand_as(lse_m_expected_alt)],
                names_input=(f"lse_m_{backend.value}_{precision}", "lse_m_ref_fp64", f"lse_m_ref_{precision}"),
            )
            assert_no_percentile_upshift(
                amax_result[mask_i.expand_as(amax_result)],
                amax_expected_fp64[mask_i.expand_as(amax_expected_fp64)].to(dtype=torch.float32),
                amax_expected_alt[mask_i.expand_as(amax_expected_alt)],
                names_input=(f"amax_{backend.value}_{precision}", "amax_ref_fp64", f"amax_ref_{precision}"),
            )
        elif backend == TriAttnBackend.TRIFAST:
            # TRIFAST returns lse directly (not lse - amax), so compare against lse_m + amax
            assert amax_result is None, "amax should be None for TRIFAST backend"
            assert_no_percentile_upshift(
                lse_m_result[mask_i.expand_as(lse_m_result)],
                (lse_m_expected_fp64 + amax_expected_fp64)[mask_i.expand_as(lse_m_expected_fp64)].to(
                    dtype=torch.float32
                ),
                (lse_m_expected_alt + amax_expected_alt)[mask_i.expand_as(lse_m_expected_alt)],
                names_input=(f"lse_{backend.value}_{precision}", "lse_ref_fp64", f"lse_ref_{precision}"),
            )
        else:
            raise ValueError(f"Unknown backend: {backend}")

        if check_bwd:
            mask_kv = mask_kv.to(dtype=bool)
            mask_j = mask_j.to(dtype=bool)
            # Test gradient error histograms
            for grad_name in ["q", "k", "v", "triangle_bias"]:
                grad_result = grads_result[grad_name]
                grad_expected_fp64 = grads_fp64[grad_name]
                grad_expected_alt = grads_alt[grad_name]
                if grad_name == "q":
                    m = mask_i.expand_as(grad_result)
                elif grad_name == "k" or grad_name == "v":
                    m = mask_kv.expand_as(grad_result)
                elif grad_name == "triangle_bias":
                    m = mask_j.expand_as(grad_result)
                else:
                    raise ValueError(f"Unknown gradient name: {grad_name}")
                assert_no_percentile_upshift(
                    grad_result[m],
                    grad_expected_fp64[m].to(dtype=torch.float32),
                    grad_expected_alt[m],
                    names_input=(
                        f"d_{grad_name}_{backend.value}_{precision}",
                        f"d_{grad_name}_ref_fp64",
                        f"d_{grad_name}_ref_{precision}",
                    ),
                )
    else:
        # === SIMPLE TOLERANCE TESTING (FP32 only) ===
        # For FP32 without TF32, use simple assertion with default tolerances
        torch.testing.assert_close(o_result * mask_i, o_expected_fp64.to(dtype=o_result.dtype) * mask_i)

        if backend in (TriAttnBackend.CUEQ, TriAttnBackend.CUEQ_FWD_TRIFAST_BWD):
            torch.testing.assert_close(lse_m_result * mask_i, lse_m_expected_fp64.to(dtype=lse_m_result.dtype) * mask_i)
            torch.testing.assert_close(amax_result * mask_i, amax_expected_fp64.to(dtype=amax_result.dtype) * mask_i)
        elif backend == TriAttnBackend.TRIFAST:
            assert amax_result is None, "amax should be None for TRIFAST backend"
            # TRIFAST returns lse directly (not lse - amax), so compare against lse_m + amax
            torch.testing.assert_close(
                lse_m_result * mask_i, (lse_m_expected_fp64 + amax_expected_fp64).to(dtype=lse_m_result.dtype) * mask_i
            )
        else:
            raise ValueError(f"Unknown backend: {backend}")

        if check_bwd:
            for grad_name in ["q", "k", "v", "triangle_bias"]:
                grad_result = grads_result[grad_name]
                grad_expected_fp64 = grads_fp64[grad_name]
                torch.testing.assert_close(
                    grad_result, grad_expected_fp64.to(dtype=grad_result.dtype), msg=lambda m: f"d_{grad_name}:\n{m}"
                )

    # check gradients are zero in masked regions to prevent backprop invalid elements' gradient upstream
    if check_bwd:
        for grad_name in ["q", "k", "v", "triangle_bias"]:
            grad_result = grads_result[grad_name]
            grad_expected_fp64 = grads_fp64[grad_name]
            grad_expected_alt = grads_alt[grad_name]
            if grad_name == "q":
                m = mask_i
            elif grad_name == "k" or grad_name == "v":
                m = mask_kv
            elif grad_name == "triangle_bias":
                m = mask_j
            else:
                raise ValueError(f"Unknown gradient name: {grad_name}")
            torch.testing.assert_close(grad_result * ~(m.bool()), torch.zeros_like(grad_result), atol=0, rtol=0)

    # === CONTIGUITY CONSISTENCY TEST ===
    # Test self-consistency between contiguous vs non-contiguous layout
    if argsort_strides is not None:
        # Run the same backend with contiguous layout for comparison
        q_contiguous = q.contiguous().detach().clone()
        k_contiguous = k.contiguous().detach().clone()
        v_contiguous = v.contiguous().detach().clone()
        triangle_bias_contiguous = triangle_bias.contiguous().detach().clone()
        do_contiguous = do.contiguous().detach().clone()
        mask_contiguous = mask.contiguous().detach().clone()

        o_result_contiguous, lse_m_result_contiguous, amax_result_contiguous, grads_result_contiguous = (
            run_triangle_attention(
                q_contiguous,
                k_contiguous,
                v_contiguous,
                triangle_bias_contiguous,
                do_contiguous,
                mask=mask_contiguous,
                backend=backend,
                precision=precision,
                check_bwd=check_bwd,
                scale=scale,
            )
        )

        torch.testing.assert_close(o_result_contiguous, o_result, atol=0, rtol=0)
        torch.testing.assert_close(lse_m_result_contiguous, lse_m_result, atol=0, rtol=0)
        torch.testing.assert_close(amax_result_contiguous, amax_result, atol=0, rtol=0)

        if check_bwd:
            for grad_name in ["q", "k", "v", "triangle_bias"]:
                # triangle_bias's gradient is not binary identical across layouts for CUEQ backend
                # due to atomic operations in the kernel implementation
                atol, rtol = (
                    (None, None) if (grad_name == "triangle_bias" and backend == TriAttnBackend.CUEQ) else (0, 0)
                )
                torch.testing.assert_close(
                    grads_result_contiguous[grad_name], grads_result[grad_name], atol=atol, rtol=rtol
                )
