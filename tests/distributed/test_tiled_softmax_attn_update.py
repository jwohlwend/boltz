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

import pytest
import torch

from boltz.distributed.model.modules.utils import DTYPE_TO_PRECISION
from boltz.distributed.utils import tiled_softmax_attention_update
from boltz.testing.utils import PRECISION_TO_INF


def expected_accum(a, v, dim_softmax, return_amax: bool = True):
    """Compute expected accumulation for validation purposes."""
    dtype_input = a.dtype
    a = a.double()
    v = v.double()
    o = torch.einsum("...i, ...i -> ...", torch.softmax(a, dim=dim_softmax), v)
    amax = a.max(dim=dim_softmax, keepdim=True)[0]
    lse_m = torch.logsumexp(a - amax, dim=dim_softmax, keepdim=True)
    if return_amax:
        return lse_m.to(dtype_input), o.reshape_as(lse_m).to(dtype_input), amax.to(dtype_input)
    else:
        # when return_amax is False, we return lse = lse_m + amax
        # but the lse_m and amax were computed as return_amax is True
        lse = lse_m + amax
        return lse.to(dtype_input), o.reshape_as(lse).to(dtype_input), None


@pytest.fixture
def softmax_test_tensors():
    """Fixture providing test tensors for online softmax accumulation tests."""
    torch.manual_seed(42)
    dtype = torch.float32

    size_chunk = 3
    n_chunks = 21

    n_elems = size_chunk * n_chunks
    size_batch = 3

    device = torch.device("cuda:0")
    a = torch.randn((size_batch, n_elems), dtype=dtype, device=device)
    v = torch.randn_like(a, device=device)

    # generate chunks with specific patterns including inf and -inf
    inf = 1e9
    a = a.reshape(size_batch, n_chunks, size_chunk)

    # Set up specific patterns for batch 0:
    # chunk 0: -inf
    # chunk 1: -inf
    # chunk 2: inf
    # chunk 3: inf
    # chunk 4: -inf
    # chunk 5: -inf
    # chunk 6: inf
    # chunk 7: inf first then -inf
    # chunk 8: -inf first then inf
    # chunk 9: inf first then randn
    # chunk 10: randn first then inf
    # chunk 11: -inf first then randn
    # chunk 12: randn first then -inf
    # chunk 13: randn
    # chunk 14: randn
    # chunk 15: randn
    # chunk 16: inf
    # chunk 17: inf
    # chunk 18: -inf
    # chunk 19: -inf
    # chunk 20: -inf
    a[0, 0:2] = -inf
    a[0, 2:4] = inf
    a[0, 4:6] = -inf
    a[0, 6] = inf
    a[0, 7, : size_chunk // 2] = inf
    a[0, 7, size_chunk // 2 :] = -inf
    a[0, 8, : size_chunk // 2] = -inf
    a[0, 8, size_chunk // 2 :] = inf
    a[0, 9, : size_chunk // 2] = inf
    a[0, 10, size_chunk // 2 :] = inf
    a[0, 11, : size_chunk // 2] = -inf
    a[0, 12, size_chunk // 2 :] = -inf
    a[0, 16:18] = inf
    a[0, 18:21] = -inf

    # inverse the pattern of batch 0 for batch 1
    a[1] = -a[0]

    # a typical masked softmax pattern for batch 2
    a[2, 10:] = -inf

    a = a.flatten(start_dim=1)

    return {
        "a": a,
        "v": v,
        "size_chunk": size_chunk,
        "n_chunks": n_chunks,
        "dim_softmax": 1,
        "device": device,
        "dtype": dtype,
    }


@pytest.mark.parametrize("has_amax", [True, False], ids=lambda x: f"has_amax:{x}")
def test_tiled_softmax_attention_update_correctness(softmax_test_tensors, has_amax):
    """Test that tiled softmax attention update produces correct results."""
    a = softmax_test_tensors["a"]
    v = softmax_test_tensors["v"]
    size_chunk = softmax_test_tensors["size_chunk"]
    n_chunks = softmax_test_tensors["n_chunks"]
    dim_softmax = softmax_test_tensors["dim_softmax"]
    device = softmax_test_tensors["device"]

    if not has_amax:
        # Without amax to take away contribution from extreme values,
        # lse_m is actually lse = lse_m + amax, where the sum often results in
        # catastrophic cancellation. In this case, the tiled_softmax_attention_update
        # can only work if the "-inf" padding in the attention score only shows up after
        # those normal values but not preceding them, i.e., we have a lse pattern of:
        # [... <stretch of normal values>, -inf, -inf, ..., -inf]. See tiled_softmax_attention_update's
        # comments for more details.
        torch.manual_seed(42)
        inf = PRECISION_TO_INF[DTYPE_TO_PRECISION[a.dtype]]
        a = torch.randn_like(a, device=device)
        a[:, a.shape[1] // 2 :] = -inf

    lse_m, o = None, None
    amax = None

    for i_chunk in range(n_chunks):
        i_begin = i_chunk * size_chunk
        i_end = (i_chunk + 1) * size_chunk
        ids_chunk = torch.arange(i_begin, i_end, device=device)
        ids_cum_chunk = torch.arange(0, i_end, device=device)

        a_chunk = a.index_select(dim_softmax, ids_chunk)
        v_chunk = v.index_select(dim_softmax, ids_chunk)

        # Compute expected cumulative results
        lse_m_cum_expected, o_cum_expected, amax_cum_expected = expected_accum(
            a.index_select(dim_softmax, ids_cum_chunk),
            v.index_select(dim_softmax, ids_cum_chunk),
            dim_softmax,
            has_amax,
        )

        # Perform online softmax accumulation
        # Step 1: compute per-block amax, lse_m, and o
        amax_chunk = a_chunk.amax(dim=dim_softmax, keepdim=True)
        a_chunk_delta = a_chunk - amax_chunk

        # Subtract out the current chunk's amax and keep it away from the
        # following computation until absolutely impossible to do so any further
        lse_m_chunk = torch.logsumexp(a_chunk_delta, dim=dim_softmax, keepdim=True)
        s_chunk = torch.exp(a_chunk_delta - lse_m_chunk)
        o_chunk = torch.einsum("...i, ...i -> ...", s_chunk, v_chunk).reshape_as(amax_chunk)

        if not has_amax:
            # when has_amax is False, we use lse instead of lse_m for accumulation across chunks
            lse_m_chunk = lse_m_chunk + amax_chunk
            amax_chunk = None

        # Update accumulated values
        o, lse_m, amax = tiled_softmax_attention_update(o_chunk, lse_m_chunk, amax_chunk, o, lse_m, amax)

        # Verify correctness against expected results
        torch.testing.assert_close(lse_m, lse_m_cum_expected)
        torch.testing.assert_close(o, o_cum_expected)
        if has_amax:
            torch.testing.assert_close(amax, amax_cum_expected)
        else:
            assert amax is None
            assert amax_cum_expected is None


def test_tiled_softmax_attention_update_error_cases():
    """Test error handling for invalid inputs."""
    device = torch.device("cuda:0")
    dtype = torch.float32

    # Create some test tensors
    o_chunk = torch.randn(3, 5, device=device, dtype=dtype)
    lse_m_chunk = torch.randn(3, 1, device=device, dtype=dtype)
    amax_chunk = torch.randn(3, 1, device=device, dtype=dtype)

    # Test case 1: Inconsistent None/not-None parameters for o and lse_m
    with pytest.raises(ValueError, match="o and lse_m must both be None or both be not None"):
        tiled_softmax_attention_update(o_chunk, lse_m_chunk, amax_chunk, o_chunk, None, None)

    # Test case 2: Shape mismatch between lse_m_chunk and amax_chunk
    wrong_shape_amax = torch.randn(3, 2, device=device, dtype=dtype)
    with pytest.raises(ValueError, match="lse_m_chunk and amax_chunk must have the same shape"):
        tiled_softmax_attention_update(o_chunk, lse_m_chunk, wrong_shape_amax, None, None, None)

    # Test case 3: lse_m_chunk doesn't have last dimension of size 1
    wrong_lse_m = torch.randn(3, 2, device=device, dtype=dtype)
    wrong_amax = torch.randn(3, 2, device=device, dtype=dtype)
    with pytest.raises(ValueError, match="lse_m_chunk must have shape \\(\\.\\.\\.\\, 1\\)"):
        tiled_softmax_attention_update(o_chunk, wrong_lse_m, wrong_amax, None, None, None)

    # Test case 4: Different number of dimensions between o_chunk and lse_m_chunk
    wrong_dim_lse_m = torch.randn(3, 1, 1, device=device, dtype=dtype)
    wrong_dim_amax = torch.randn(3, 1, 1, device=device, dtype=dtype)
    with pytest.raises(ValueError, match="o_chunk and lse_m_chunk must have the same number of dimensions"):
        tiled_softmax_attention_update(o_chunk, wrong_dim_lse_m, wrong_dim_amax, None, None, None)

    # Test case 5: Shape mismatch between o_chunk and lse_m_chunk (except last dimension)
    wrong_batch_o = torch.randn(4, 5, device=device, dtype=dtype)  # Different batch size
    with pytest.raises(
        ValueError, match="o_chunk and lse_m_chunk must have the same shape except for the last dimension"
    ):
        tiled_softmax_attention_update(wrong_batch_o, lse_m_chunk, amax_chunk, None, None, None)

    # Test case 6: Shape mismatch between o_chunk and o (non-initial chunk)
    o_accum = torch.randn(3, 5, device=device, dtype=dtype)
    lse_m_accum = torch.randn(3, 1, device=device, dtype=dtype)
    amax_accum = torch.randn(3, 1, device=device, dtype=dtype)
    wrong_shape_o_chunk = torch.randn(3, 6, device=device, dtype=dtype)  # Different feature dimension
    wrong_shape_lse_m_chunk = torch.randn(3, 1, device=device, dtype=dtype)
    with pytest.raises(ValueError, match="o_chunk and o must have the same shape"):
        tiled_softmax_attention_update(
            wrong_shape_o_chunk, wrong_shape_lse_m_chunk, amax_chunk, o_accum, lse_m_accum, amax_accum
        )

    # Test case 7: Shape mismatch between lse_m_chunk and lse_m (non-initial chunk)
    wrong_shape_lse_m_chunk = torch.randn(4, 1, device=device, dtype=dtype)  # Different batch size
    with pytest.raises(ValueError, match="lse_m_chunk and amax_chunk must have the same shape"):
        tiled_softmax_attention_update(o_chunk, wrong_shape_lse_m_chunk, amax_chunk, o_accum, lse_m_accum, amax_accum)

    # Test case 8: Inconsistent amax and amax_chunk (non-initial chunk)
    with pytest.raises(
        ValueError, match="amax and amax_chunk must both be None or both be not None for non-initial chunks"
    ):
        tiled_softmax_attention_update(o_chunk, lse_m_chunk, amax_chunk, o_accum, lse_m_accum, None)

    # Test case 9: Shape mismatch between amax_chunk and amax (non-initial chunk)
    wrong_shape_amax_accum = torch.randn(3, 2, device=device, dtype=dtype)
    with pytest.raises(ValueError, match="amax_chunk and amax must have the same shape"):
        tiled_softmax_attention_update(o_chunk, lse_m_chunk, amax_chunk, o_accum, lse_m_accum, wrong_shape_amax_accum)

    # Test case 10: Shape mismatch between lse_m_chunk and lse_m (non-initial chunk)
    wrong_shape_lse_m_accum = torch.randn(3, 2, device=device, dtype=dtype)
    with pytest.raises(ValueError, match="lse_m_chunk and lse_m must have the same shape"):
        tiled_softmax_attention_update(o_chunk, lse_m_chunk, amax_chunk, o_accum, wrong_shape_lse_m_accum, amax_accum)


@pytest.mark.parametrize(
    "dtype",
    [torch.bfloat16, torch.float16],
    ids=lambda d: str(d).split(".")[-1],
)
@pytest.mark.parametrize("has_amax", [True, False], ids=lambda x: f"has_amax:{x}")
def test_tiled_softmax_attention_update_dtype_preservation(dtype, has_amax):
    """Regression: under autocast, torch.logsumexp promotes BF16/FP16 → FP32.

    torch.logsumexp preserves dtype without autocast but promotes to FP32
    when autocast is active (logsumexp is on autocast's FP32-promotion list).
    In production, _RingMultiHeadTriangleAttentionImpl.forward uses
    @custom_fwd without cast_inputs, which preserves the caller's autocast
    context — so the promotion occurs during BF16-mixed training.

    Without the fix, the has_amax=True path promoted lse_m to FP32 on step 1
    (via logsumexp), which then infected delta_lse → sigmoid → o on step 2+.
    The has_amax=False path (using logsigmoid) was not affected.

    This test wraps the calls in autocast and feeds 5 chunks in the input
    dtype, asserting that o, lse_m, and amax preserve that dtype after every
    step — including step 2+ where the cascade used to occur.
    """
    device = torch.device("cuda:0")
    n_chunks = 5
    batch, feat = 4, 8
    torch.manual_seed(0)

    o, lse_m, amax = None, None, None
    with torch.amp.autocast("cuda", dtype=dtype):
        for _ in range(n_chunks):
            o_chunk = torch.randn(batch, feat, device=device, dtype=dtype)
            lse_m_chunk = torch.randn(batch, 1, device=device, dtype=dtype)
            amax_chunk = torch.randn(batch, 1, device=device, dtype=dtype) if has_amax else None

            o, lse_m, amax = tiled_softmax_attention_update(o_chunk, lse_m_chunk, amax_chunk, o, lse_m, amax)

            assert o.dtype == dtype, f"o promoted to {o.dtype}"
            assert lse_m.dtype == dtype, f"lse_m promoted to {lse_m.dtype}"
            if has_amax:
                assert amax.dtype == dtype, f"amax promoted to {amax.dtype}"
            else:
                assert amax is None

    assert o.isfinite().all(), "o contains non-finite values"
    assert lse_m.isfinite().all(), "lse_m contains non-finite values"
