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

"""Unit tests for window batching utility functions."""

import math

import pytest
import torch

from boltz.distributed.model.layers.utils import (
    gather_sliding_windows,
    gather_sliding_windows_backward,
    get_query_window_key_range,
    pack_and_pad,
    pack_and_pad_backward,
)
from boltz.model.modules.encoders import get_indexing_matrix
from boltz.testing.utils import assert_tensors_identical


def set_batch_diagonal(batch_matrix, k_values, fill_value):
    """
    Sets the k-th diagonal for each matrix in a batch.
    Space complexity: O(B * min(H,W))
    """
    B, H, W = batch_matrix.shape
    device = batch_matrix.device

    # 1. The maximum possible length of any diagonal is min(H, W)
    max_diag_len = min(H, W)

    # 2. Create a base sequence [0, 1, 2, ..., max_len-1]
    # Shape: (1, max_diag_len)
    seq = torch.arange(max_diag_len, device=device).unsqueeze(0)

    # 3. Calculate starting coordinates (r, c) for each k
    # If k > 0: start at (0, k)
    # If k < 0: start at (|k|, 0)
    # Shape: (B, 1)
    start_row = (-k_values).clamp(min=0).unsqueeze(1)
    start_col = k_values.clamp(min=0).unsqueeze(1)

    # 4. Generate the full coordinate grids
    # We broaden the starting points by adding the sequence
    # Shape: (B, max_diag_len)
    rows = start_row + seq
    cols = start_col + seq

    # 5. Create a mask for valid coordinates
    # Because diagonals shift, they might hit the boundary before max_diag_len
    valid_mask = (rows < H) & (cols < W)

    # 6. Create Batch indices to match
    # Shape: (B, max_diag_len)
    batch_idx = torch.arange(B, device=device).unsqueeze(1).expand(-1, max_diag_len)

    # 7. Apply Advanced Indexing
    # We select only the valid coordinates using the boolean mask.
    # PyTorch handles the memory layout mapping here internally.
    batch_matrix[batch_idx[valid_mask], rows[valid_mask], cols[valid_mask]] = fill_value

    return batch_matrix


@pytest.fixture(params=["cpu", "cuda"])
def device(request):
    if request.param == "cuda" and not torch.cuda.is_available():
        pytest.skip("CUDA not available")
    return torch.device(request.param)


@pytest.fixture(
    params=[
        (8, 16, 1),
        (8, 16, 2),
        (8, 16, 5),
        (8, 16, 10),
        (8, 32, 2),
        (8, 32, 5),
        (8, 32, 20),
        (16, 64, 3),
        (16, 64, 10),
        (16, 64, 50),
        (32, 64, 5),
        (32, 64, 20),
        (32, 128, 1),
        (32, 128, 2),
        (32, 128, 3),
        (32, 128, 10),
        (32, 128, 50),
        (32, 128, 100),
        (32, 256, 5),
        (32, 256, 20),
        (64, 256, 10),
        (64, 256, 50),
        (128, 512, 20),
        (128, 512, 100),
    ],
    ids=lambda x: f"W:{x[0]}, H:{x[1]}, K:{x[2]}",
)
def get_toeplitz(request, device):
    """
    Fixture that creates onehot tensor and verifies Toeplitz property.

    Returns (W, H, K, h, batched_toeplitz) where batched_toeplitz has shape (K, h, 2*K).

    Verifies one-hot property: each (query_window, slot) has at most one non-zero.
    """
    W, H, K = request.param
    h = H // (W // 2)

    idx_mat = get_indexing_matrix(K, W, H, device)
    batched_toeplitz = idx_mat.unflatten(-1, (K, h)).transpose(0, 1).transpose(-2, -1).to(dtype=torch.int32)
    # Shape: (K, h, 2*K)

    assert (
        batched_toeplitz.shape == (K, h, 2 * K)
    ), f"get_indexing_matrix({K}, {W}, {H}) produces incorrect shape (post transformation): {batched_toeplitz.shape} != (K, h, 2 * K)"

    # Verify the Toeplitz property
    # All Toeplitz inside batched batched_toeplitz should have one non-zero diagonal of ones,
    # with the first Toeplitz matrix has diagonal at offset[0] = 1 - h // 2 and subsequent
    # Toeplitz matrices have diagonal at offset[i] = offset[i-1] + 2.
    batched_toeplitz_expected = set_batch_diagonal(
        torch.zeros((K, h, 2 * K), device=batched_toeplitz.device, dtype=batched_toeplitz.dtype),
        torch.arange(1 - h // 2, 1 - h // 2 + (K - 1) * 2 + 1, 2, device=batched_toeplitz.device),
        1,
    )

    assert torch.all(
        batched_toeplitz == batched_toeplitz_expected
    ), f"get_indexing_matrix({K}, {W}, {H}) does not produce the expected Toeplitz matrix"

    return W, H, K, h, batched_toeplitz


def test_range_formula_all_windows(get_toeplitz):
    """Test that range formula is correct for ALL query windows using batched call."""
    W, H, K, h, batched_toeplitz = get_toeplitz
    device = batched_toeplitz.device

    # Get ranges for ALL query windows in one batched call
    all_ids = torch.arange(K, device=device)
    ranges = get_query_window_key_range(W, H, K, all_ids)

    # Verify shape
    assert ranges.shape == (2, K)

    # Extract ground truth from onehot using batched sparse COO (no loops!)
    onehot_sparse = batched_toeplitz.to_sparse_coo()
    indices = onehot_sparse.indices()  # Shape: (3, num_nonzeros)
    # indices[0] = query window index (i)
    # indices[1] = slot index
    # indices[2] = half-window index (j)

    qw_idx = indices[0]
    j_idx = indices[2]

    # Use scatter_reduce to compute min/max per query window (PyTorch 1.12+)
    expected_j_min = torch.full((K,), 2 * K, dtype=torch.long, device=device)
    expected_j_max = torch.full((K,), -1, dtype=torch.long, device=device)

    expected_j_min.scatter_reduce_(0, qw_idx, j_idx, reduce="amin", include_self=False)
    expected_j_max.scatter_reduce_(0, qw_idx, j_idx, reduce="amax", include_self=False)

    # Filter to valid windows (those with at least one non-zero)
    valid_mask = expected_j_max >= 0

    assert torch.all(ranges[0, valid_mask] == expected_j_min[valid_mask]), f"W={W},H={H},K={K}: j_min mismatch"
    assert torch.all(ranges[1, valid_mask] == expected_j_max[valid_mask]), f"W={W},H={H},K={K}: j_max mismatch"


@pytest.mark.parametrize("ndim,axis", [(2, 0), (3, 1), (3, -2), (4, 2), (5, 3), (5, -2)])
def test_efficient_unfold_equivalence(get_toeplitz, ndim, axis):
    """
    Test gather_sliding_windows matches einsum with Toeplitz matrix.

    Tests with inputs of varying dimensions and axis positions.
    """
    W, H, K, h, batched_toeplitz = get_toeplitz
    device = batched_toeplitz.device

    # Build shape with 2*K at the specified axis
    shape_list = [2, 3, 4, 5, 6][:ndim]

    # Normalize axis
    norm_axis = axis if axis >= 0 else ndim + axis
    shape_list[norm_axis] = 2 * K

    input_shape = tuple(shape_list)
    dense_input = torch.arange(math.prod(input_shape), device=device, dtype=torch.float32).reshape(input_shape)
    dense_input.requires_grad_(True)

    # Method 1: Generic einsum by moving axis to front
    # Move axis dimension to position 0
    input_permuted = dense_input.moveaxis(norm_axis, 0)  # (2*K, ...)

    # Generic einsum: (K, h, 2*K) × (2*K, ...) → (K, h, ...)
    result_einsum = torch.einsum("kij,j...->ki...", batched_toeplitz.float(), input_permuted)

    # Move K and h dimensions to where axis was
    # result_einsum is (K, h, ...) - move to (...[:axis], K, h, ...[axis:])
    result_einsum = result_einsum.moveaxis([0, 1], [norm_axis, norm_axis + 1])

    # Method 2: efficient_toeplitz_matmul_unfold
    offset_start = 1 - h // 2
    offsets = torch.arange(offset_start, offset_start + 2 * (K - 1) + 1, 2, device=device)

    dense_input_clone = dense_input.detach().clone().requires_grad_(True)
    result_unfold = gather_sliding_windows(dense_input_clone, offsets, h, axis)

    # Verify equivalence
    # Forward is exact copy vs einsum (multiply by 1.0)
    # Should be close, maybe not bitwise identical on GPU
    torch.testing.assert_close(result_einsum, result_unfold)

    # Verify backward pass
    grad_output = torch.arange(result_einsum.numel(), device=device, dtype=result_einsum.dtype).reshape_as(
        result_einsum
    )

    result_einsum.backward(grad_output.detach().clone())
    result_unfold.backward(grad_output.detach().clone())

    torch.testing.assert_close(dense_input.grad, dense_input_clone.grad)


@pytest.mark.parametrize(
    "W,H,K,qw_start,qw_end",
    [
        # Test first windows (includes QW0 with negative offset)
        (32, 128, 10, 0, 3),
        (32, 128, 20, 0, 5),
        # Test last windows (includes final QW with boundary)
        (32, 128, 10, 7, 10),
        (32, 128, 20, 15, 20),
        # Test middle windows (interior, no boundaries)
        (32, 128, 20, 5, 15),
        (32, 128, 50, 20, 30),
        # Test single window
        (32, 128, 10, 3, 4),
        # Different h values
        (32, 64, 10, 2, 6),  # h=4
        (32, 256, 10, 2, 6),  # h=16
        # Larger K
        (32, 128, 100, 40, 60),
    ],
)
def test_translational_symmetry(W, H, K, qw_start, qw_end, device):
    """
    Test Theorem 6: Translational symmetry of Toeplitz multiplication.

    Verifies: T(x[δ:δ+n], offsets - δ) == T(x, offsets)[slice]

    For a subset of query windows computed on a translated input slice,
    the result equals slicing the full computation.
    """
    h = H // (W // 2)

    # Full computation
    torch.manual_seed(42)
    input_full = torch.randn(2 * K, 16, device=device, requires_grad=True)

    offset_start = 1 - h // 2
    offsets_full = torch.arange(offset_start, offset_start + 2 * K, 2, device=device)
    result_full = gather_sliding_windows(input_full, offsets_full, h, axis=0)

    # Determine input span needed for subset
    subset_qw_ids = torch.arange(qw_start, qw_end, device=device)
    ranges = get_query_window_key_range(W, H, K, subset_qw_ids)
    hw_need_start = ranges[0].min().item()
    hw_need_end = ranges[1].max().item() + 1

    # Extract input slice and translate offsets
    input_slice = input_full[hw_need_start:hw_need_end].detach().clone().requires_grad_(True)
    offsets_subset = offsets_full[qw_start:qw_end] - hw_need_start

    # Compute on translated slice
    result_subset = gather_sliding_windows(input_slice, offsets_subset, h, axis=0)

    # Verify: T(x[δ:], offsets-δ) == T(x, offsets)[slice]
    expected_subset = result_full[qw_start:qw_end]
    assert torch.all(result_subset == expected_subset), (
        f"W={W},H={H},K={K},QW[{qw_start},{qw_end}): Forward mismatch: \n"
        f"    {result_subset} \n vs. \n"
        f"    {expected_subset}"
    )

    # Verify backward pass
    grad_output = torch.randn_like(result_subset)

    grad_full = torch.zeros_like(result_full)
    grad_full[qw_start:qw_end] = grad_output
    result_full.backward(grad_full)

    result_subset.backward(grad_output.clone())

    expected_grad_slice = input_full.grad[hw_need_start:hw_need_end]
    # Backward involves gradient accumulation which can have small numerical differences on GPU
    torch.testing.assert_close(input_slice.grad, expected_grad_slice)


def test_backward_validation_errors(device):
    """Test that gather_sliding_windows_backward raises appropriate errors for invalid inputs."""

    # Valid baseline
    grad_output = torch.randn(5, 8, 16, device=device)  # (n_windows=5, window_size=8, features=16)
    window_start_offsets = torch.tensor([-3, -1, 1, 3, 5], device=device)
    window_size = 8
    axis = 0
    input_shape = (12, 16)  # (2*K=12, features=16)

    # Test 1: grad_output not a tensor
    with pytest.raises(TypeError, match="grad_output must be a torch.Tensor"):
        gather_sliding_windows_backward([1, 2, 3], window_start_offsets, window_size, axis, input_shape)

    # Test 2: window_start_offsets not a tensor
    with pytest.raises(TypeError, match="window_start_offsets must be a torch.Tensor"):
        gather_sliding_windows_backward(grad_output, [1, 2, 3], window_size, axis, input_shape)

    # Test 3: window_start_offsets not 1D
    with pytest.raises(ValueError, match="window_start_offsets must be 1D"):
        bad_offsets = torch.tensor([[1, 2], [3, 4]], device=device)
        gather_sliding_windows_backward(grad_output, bad_offsets, window_size, axis, input_shape)

    # Test 4: axis out of range
    with pytest.raises(ValueError, match="axis .* out of range"):
        gather_sliding_windows_backward(grad_output, window_start_offsets, window_size, axis=5, input_shape=input_shape)

    # Test 5: grad_output shape mismatch (wrong n_windows dimension)
    with pytest.raises(ValueError, match="grad_output shape mismatch"):
        bad_grad = torch.randn(7, 8, 16, device=device)  # Wrong n_windows (7 instead of 5)
        gather_sliding_windows_backward(bad_grad, window_start_offsets, window_size, axis, input_shape)

    # Test 6: grad_output shape mismatch (wrong window_size dimension)
    with pytest.raises(ValueError, match="grad_output shape mismatch"):
        bad_grad = torch.randn(5, 10, 16, device=device)  # Wrong window_size (10 instead of 8)
        gather_sliding_windows_backward(bad_grad, window_start_offsets, window_size, axis, input_shape)

    # Test 7: grad_output shape mismatch (wrong feature dimension)
    with pytest.raises(ValueError, match="grad_output shape mismatch"):
        bad_grad = torch.randn(5, 8, 32, device=device)  # Wrong features (32 instead of 16)
        gather_sliding_windows_backward(bad_grad, window_start_offsets, window_size, axis, input_shape)

    # Test 8: grad_output wrong ndim
    with pytest.raises(ValueError, match="grad_output shape mismatch"):
        bad_grad = torch.randn(5, 8, device=device)  # Missing feature dimension
        gather_sliding_windows_backward(bad_grad, window_start_offsets, window_size, axis, input_shape)


@pytest.mark.parametrize(
    "input_shape_extra", [(None,), (4, None), (None, 3), (2, None, 3)], ids=lambda x: f"input_shape_extra:{x}"
)
@pytest.mark.parametrize("keep_input_padding", [False, True], ids=lambda x: f"keep_input_padding:{x}")
def test_pack_and_pad_equivalence(get_toeplitz, input_shape_extra, keep_input_padding):
    """
    Test pack_and_pad utility.

    `None` in fixture `input_shape_extra` is eventually replaced with K * W.
    For example, `input_shape_extra (2, None, 3)` will be replaced with `(2, K * W, 3)` in the test.
    """
    W, H, K, h, batched_toeplitz = get_toeplitz
    device = batched_toeplitz.device

    # Setup parameters
    n_axes_none = sum(1 for x in input_shape_extra if x is None)
    if n_axes_none != 1:
        raise ValueError(f"There can be one and only one 'None' element in the input_shape but got {input_shape_extra}")
    axis = input_shape_extra.index(None)

    input_shape = input_shape_extra[:axis] + (K * W,) + input_shape_extra[axis + 1 :]

    # 1. Generate clean input of shape (K*W, features)
    # This represents the "perfectly padded" sequence
    torch.manual_seed(42)
    input = torch.randn(input_shape, device=device, requires_grad=True)
    input_copy = input.detach().clone().requires_grad_(True)

    mask = torch.randint(0, 2, input_shape, dtype=torch.bool, device=device, requires_grad=False)
    mask_copy = mask.detach().clone().requires_grad_(False)

    mask_sorted, argsort_mask = torch.sort(mask, dim=axis, descending=True, stable=True)
    input_sorted = torch.gather(input, axis, argsort_mask)

    # --- Reference Computation (on clean input) ---
    # Reshape to (2*K, W//2, features) for einsum
    # Note: K*W elements -> reshaped to (2*K, W/2)
    # 2*K * W/2 = K*W. Correct.

    # input_clean: (K*W, F) -> (2*K, W//2, F)
    # We need to ensure the layout matches what the reshape produces.
    # The utility reshapes (K*W) -> (2*K, W//2) along axis.
    # This splits the sequence into chunks of size W/2.
    input_sorted_reshaped = input_sorted.moveaxis(axis, -1).unflatten(-1, (2 * K, W // 2))

    # Einsum: (K, h, 2*K) x (..., 2*K, W//2) -> (..., K, h, W//2)
    # Sum over 2*K dimension
    # batched_toeplitz: (K, h, 2*K)
    ref_output_reshaped = torch.einsum("khj,...jb->...khb", batched_toeplitz.float(), input_sorted_reshaped)
    # (..., K, h, W//2) -> (..., K, h, W//2, ...)
    ref_output = ref_output_reshaped.moveaxis([-3, -2, -1], [axis, axis + 1, axis + 2])
    with torch.no_grad():
        mask_prepared_ref = mask_sorted.unflatten(axis, (2 * K, W // 2))
        mask_output_ref = (
            torch.einsum(
                "khj,...jb->...khb",
                batched_toeplitz.to(dtype=torch.float32),
                # (..., 2 * K, W // 2, ...) -> (..., 2*K, W//2)
                mask_prepared_ref.to(dtype=torch.float32).moveaxis([axis, axis + 1], [-2, -1]),
            ).moveaxis([-3, -2, -1], [axis, axis + 1, axis + 2])
        ).to(dtype=mask_prepared_ref.dtype)

    # Reshape reference result to match gather_sliding_windows output structure
    # Reference einsum output: (..., K, h, W//2, ...)
    # gather_sliding_windows output: (..., n_windows, window_size, ...)
    # For axis=0 input (2*K, W//2, F), it returns (K, h, W//2, F)
    # Structure matches exactly.

    # --- Target Computation ---
    # 1. pack the valid elements and pad to the next multiple of W
    input_prepared, _, mask_prepared = pack_and_pad(
        input_copy, mask_copy, axis, W, keep_input_padding=keep_input_padding
    )

    # check the mask
    assert mask_prepared.shape == input_prepared.shape
    n_valid = mask_copy.expand_as(input_copy).sum(dim=axis)
    # there must be leading n_valid elements and trailing zeros along axis
    # NOTE: cumprod zeros out non-leading True elements
    assert_tensors_identical(mask_prepared.cumprod(dim=axis).sum(dim=axis), n_valid)

    if keep_input_padding:
        mask_prepared_padded = mask_prepared
        input_prepared_padded = input_prepared
    else:
        # when keep_input_padding is False, the mask_prepared will be shorter than the reference
        # along 'axis'. We pad them before along 'axis' with zeros towards the reference length
        pad_len = mask_prepared_ref.flatten(axis, axis + 1).shape[axis] - mask_prepared.shape[axis]
        assert pad_len >= 0, "Padding length should be non-negative for the result when keep_input_padding is False"
        pad_arg = [0] * (2 * mask_prepared.ndim)
        pad_idx = (mask_prepared.ndim - 1 - axis) * 2 + 1
        pad_arg[pad_idx] = pad_len
        input_prepared_padded = torch.nn.functional.pad(input_prepared, pad_arg)
        mask_prepared_padded = torch.nn.functional.pad(mask_prepared, pad_arg)

    # reshape to (..., 2*K, W//2, ...)
    input_prepared_padded = input_prepared_padded.unflatten(axis, (2 * K, W // 2))
    mask_prepared_padded = mask_prepared_padded.unflatten(axis, (2 * K, W // 2))

    assert_tensors_identical(mask_prepared_padded, mask_prepared_ref)

    # 2. Gather Sliding Windows
    # We need offsets for K windows.
    # offset_start = 1 - h // 2
    offset_start = 1 - h // 2
    window_start_offsets = torch.arange(offset_start, offset_start + 2 * K, 2, device=device)

    target_output = gather_sliding_windows(input_prepared_padded, window_start_offsets, h, axis)

    # --- Verification ---

    # 1. Verify forward output
    # target_output: (..., K, h, W//2, ...)
    # Should be binary identical (no fp math involved, just gather/move)
    torch.testing.assert_close(target_output, ref_output * mask_output_ref, atol=0, rtol=0)

    # 2. Verify backward (gradients)
    # We'll compute gradients w.r.t. input_clean vs input_masked

    grad_out = torch.randn_like(ref_output)

    # NOTE: the reference computation doesn't use mask at all so we need to mask out the
    # invalid upstream adjoints
    with torch.no_grad():
        grad_out = grad_out * mask_output_ref

    # Reference backward
    # 1. Backprop through Toeplitz path
    torch.autograd.backward([ref_output], [grad_out])

    # Target backward
    torch.autograd.backward([target_output], [grad_out])

    # Backward involves gradient accumulation (summation) which is not bitwise identical
    # between einsum (reference) and index_add_ (target) due to floating point associativity.
    torch.testing.assert_close(input_copy.grad, input.grad)

    # Verify padded/invalid elements have zero gradient
    assert_tensors_identical(torch.zeros_like(input_copy.grad), input_copy.grad * ~mask_copy)


def test_pack_and_pad_backward_manual(get_toeplitz):
    """Test manual backward pass for pack_and_pad."""
    W, H, K, h, batched_toeplitz = get_toeplitz
    device = batched_toeplitz.device
    axis = 0
    features = 16

    # Setup inputs
    torch.manual_seed(42)
    input_clean = torch.randn(K * W, features, device=device, requires_grad=True)

    num_invalid = 5
    total_len = K * W + num_invalid
    mask_indices = torch.randperm(total_len, device=device)
    valid_indices = mask_indices[: K * W]
    valid_indices_sorted, _ = torch.sort(valid_indices)

    input_masked = torch.zeros(total_len, features, device=device)
    input_masked[valid_indices_sorted] = input_clean.detach()
    input_masked.requires_grad_(True)

    mask = torch.zeros(total_len, dtype=torch.bool, device=device)
    mask[valid_indices_sorted] = True
    mask_reshaped = mask.reshape(total_len, 1)

    # Forward
    output, indices, mask_output = pack_and_pad(input_masked, mask_reshaped, axis, W)

    # Compute grad using autograd
    grad_out = torch.randn_like(output)

    torch.autograd.backward([output], [grad_out])

    grad_autograd = input_masked.grad.clone()
    input_masked.grad.zero_()

    # Compute grad using manual backward
    grad_manual = pack_and_pad_backward(grad_out, mask_output, indices, input_masked.shape, axis)

    # Compare
    torch.testing.assert_close(grad_manual, grad_autograd)
