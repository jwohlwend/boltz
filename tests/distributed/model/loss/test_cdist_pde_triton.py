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

import os
import re
import subprocess
from pathlib import Path

import pytest
import torch
import torch.nn.functional as F

from boltz.distributed.model.loss.triton.cdist_pde import cdist_pde
from boltz.testing.utils import init_tensors_uniform


def cdist_pde_reference(
    pred_pde,
    true_coords_row,
    true_coords_col,
    pred_coords_row,
    pred_coords_col,
    mask_row,
    mask_col,
    multiplicity,
    num_bins=64,
    max_dist=32.0,
):
    """
    Reference implementation for PDE cross-entropy loss.

    This implements the equivalent computation as the Triton kernel but using
    standard PyTorch operations, materializing the full distance matrices.

    Returns fully summed outputs [B_mul] (sum over both row and column dimensions).
    """
    B_mul, N_row, N_col, _ = pred_pde.shape
    B = B_mul // multiplicity

    # Broadcast masks to B_mul if needed
    if mask_row.shape[0] == B:
        mask_row_expanded = mask_row.repeat_interleave(multiplicity, dim=0)
    else:
        mask_row_expanded = mask_row

    if mask_col.shape[0] == B:
        mask_col_expanded = mask_col.repeat_interleave(multiplicity, dim=0)
    else:
        mask_col_expanded = mask_col

    # Compute pair mask [B_mul, N_row, N_col]
    mask = mask_row_expanded.unsqueeze(-1) * mask_col_expanded.unsqueeze(-2)

    # Compute distances
    true_d = torch.cdist(true_coords_row, true_coords_col)
    pred_d = torch.cdist(pred_coords_row, pred_coords_col)
    target_pde = torch.abs(true_d - pred_d)

    # Compute bin indices
    bin_index = torch.floor(target_pde * num_bins / max_dist).long()
    bin_index = torch.clamp(bin_index, max=(num_bins - 1))

    # Compute cross-entropy
    pde_one_hot = F.one_hot(bin_index, num_classes=num_bins).float()
    log_probs = F.log_softmax(pred_pde, dim=-1)
    errors = -1 * torch.sum(pde_one_hot * log_probs, dim=-1)

    # Full sum over both row and column dimensions
    out_loss_num = torch.sum(errors * mask, dim=(-2, -1))  # [B_mul]
    out_mask_denom = torch.sum(mask, dim=(-2, -1))  # [B_mul]

    return out_loss_num, out_mask_denom


@pytest.mark.skipif(not torch.cuda.is_available(), reason="CUDA not available")
@pytest.mark.parametrize("B", [1, 2], ids=lambda x: f"B:{x}")
@pytest.mark.parametrize("multiplicity", [1, 4], ids=lambda x: f"M:{x}")
@pytest.mark.parametrize("N_row", [32, 64], ids=lambda x: f"Nr:{x}")
@pytest.mark.parametrize("N_col", [32, 64], ids=lambda x: f"Nc:{x}")
@pytest.mark.parametrize("mask_row_has_mul", [False, True], ids=lambda x: f"MaskRowMul:{x}")
@pytest.mark.parametrize("mask_col_has_mul", [False, True], ids=lambda x: f"MaskColMul:{x}")
@pytest.mark.parametrize("all_zero_mask", [False, True], ids=lambda x: f"AllZeroMask:{x}")
def test_cdist_pde_correctness(B, multiplicity, N_row, N_col, mask_row_has_mul, mask_col_has_mul, all_zero_mask):
    """Test forward and backward pass correctness against reference implementation."""
    device = torch.device("cuda")
    B_mul = B * multiplicity
    num_bins = 64
    max_dist = 32.0
    min_val, max_val = -1.0, 1.0

    # Generate random coordinates using init_tensors_uniform
    true_coords_row = torch.empty(B_mul, N_row, 3, device=device)
    true_coords_col = torch.empty(B_mul, N_col, 3, device=device)
    pred_coords_row = torch.empty(B_mul, N_row, 3, device=device)
    pred_coords_col = torch.empty(B_mul, N_col, 3, device=device)
    init_tensors_uniform(
        [true_coords_row, true_coords_col, pred_coords_row, pred_coords_col],
        low=min_val,
        high=max_val,
    )

    # pred_pde needs gradient for backward test
    pred_pde_ref = torch.empty(B_mul, N_row, N_col, num_bins, device=device)
    init_tensors_uniform([pred_pde_ref], low=min_val, high=max_val)
    pred_pde_ref.requires_grad_(True)
    pred_pde_triton = pred_pde_ref.detach().clone().requires_grad_(True)

    # Create masks with specified shapes
    if all_zero_mask:
        # All-zero masks for edge case testing
        mask_row_shape = (B_mul, N_row) if mask_row_has_mul else (B, N_row)
        mask_col_shape = (B_mul, N_col) if mask_col_has_mul else (B, N_col)
        mask_row = torch.zeros(mask_row_shape, device=device)
        mask_col = torch.zeros(mask_col_shape, device=device)
    else:
        # Random binary masks
        if mask_row_has_mul:
            mask_row = torch.randint(0, 2, (B_mul, N_row), device=device).float()
        else:
            mask_row = torch.randint(0, 2, (B, N_row), device=device).float()

        if mask_col_has_mul:
            mask_col = torch.randint(0, 2, (B_mul, N_col), device=device).float()
        else:
            mask_col = torch.randint(0, 2, (B, N_col), device=device).float()

        # Ensure at least some masks are 1 for gradient flow
        mask_row.view(-1)[0] = 1.0
        mask_col.view(-1)[0] = 1.0

    # ===== Forward pass test =====
    # Reference computation
    ref_loss_num, ref_mask_denom = cdist_pde_reference(
        pred_pde_ref,
        true_coords_row,
        true_coords_col,
        pred_coords_row,
        pred_coords_col,
        mask_row,
        mask_col,
        multiplicity,
        num_bins,
        max_dist,
    )

    # Triton computation
    triton_loss_num, triton_mask_denom = cdist_pde(
        pred_pde_triton,
        true_coords_row,
        true_coords_col,
        pred_coords_row,
        pred_coords_col,
        mask_row,
        mask_col,
        multiplicity,
        num_bins,
        max_dist,
    )

    # Compare forward outputs
    torch.testing.assert_close(triton_loss_num, ref_loss_num)
    torch.testing.assert_close(triton_mask_denom, ref_mask_denom)

    # ===== Backward pass test =====
    if all_zero_mask:
        # With all-zero mask, outputs should be zero, skip backward test
        assert torch.all(triton_loss_num == 0)
        assert torch.all(triton_mask_denom == 0)
        return

    # Create upstream gradient (mock adjoint) using init_tensors_uniform
    grad_out = torch.empty_like(ref_loss_num)
    init_tensors_uniform([grad_out], low=min_val, high=max_val)

    # Backward pass with upstream gradient
    ref_loss_num.backward(grad_out)
    triton_loss_num.backward(grad_out)

    # Compare gradients
    assert pred_pde_ref.grad is not None
    assert pred_pde_triton.grad is not None
    torch.testing.assert_close(pred_pde_triton.grad, pred_pde_ref.grad)


@pytest.fixture(
    params=[
        # (modifications_dict, expected_error_pattern)
        ({}, None),  # valid inputs
        ({"pred_pde": (4, 32, 32, 32)}, "pred_pde num_bins mismatch"),  # wrong num_bins
        ({"coord_dim": 4}, "Coordinate dimension must be 3"),
        ({"true_coords_row": (4, 64, 3)}, "true_coords_row shape mismatch"),
        ({"true_coords_col": (4, 64, 3)}, "true_coords_col shape mismatch"),
        ({"pred_coords_row": (4, 64, 3)}, "pred_coords_row shape mismatch"),
        ({"pred_coords_col": (4, 64, 3)}, "pred_coords_col shape mismatch"),
        ({"mask_row": (3, 32)}, "mask_row batch dimension"),  # Neither B nor B_mul
        ({"mask_col": (3, 32)}, "mask_col batch dimension"),  # Neither B nor B_mul
        ({"mask_row": (2, 64)}, "mask_row N dimension"),  # Wrong N_row
        ({"mask_col": (2, 64)}, "mask_col N dimension"),  # Wrong N_col
        # requires_grad validation (gradient flow is broken by bin_index computation)
        ({"true_coords_row_requires_grad": True}, "true_coords_row should not require gradients"),
        ({"true_coords_col_requires_grad": True}, "true_coords_col should not require gradients"),
        ({"pred_coords_row_requires_grad": True}, "pred_coords_row should not require gradients"),
        ({"pred_coords_col_requires_grad": True}, "pred_coords_col should not require gradients"),
        ({"mask_row_requires_grad": True}, "mask_row should not require gradients"),
        ({"mask_col_requires_grad": True}, "mask_col should not require gradients"),
    ],
    ids=[
        "valid_inputs",
        "pred_pde_wrong_num_bins",
        "coord_dim_not_3",
        "true_coords_row_wrong_shape",
        "true_coords_col_wrong_shape",
        "pred_coords_row_wrong_shape",
        "pred_coords_col_wrong_shape",
        "mask_row_batch_neither_B_nor_B_mul",
        "mask_col_batch_neither_B_nor_B_mul",
        "mask_row_wrong_n_row",
        "mask_col_wrong_n_col",
        "true_coords_row_requires_grad",
        "true_coords_col_requires_grad",
        "pred_coords_row_requires_grad",
        "pred_coords_col_requires_grad",
        "mask_row_requires_grad",
        "mask_col_requires_grad",
    ],
)
def pde_validation_case(request):
    """Fixture providing (modifications, expected_error) for validation tests"""
    return request.param


@pytest.mark.skipif(not torch.cuda.is_available(), reason="CUDA not available")
def test_cdist_pde_input_validation(pde_validation_case):
    """Test input validation for cdist_pde()"""
    modifications, expected_error = pde_validation_case
    device = torch.device("cuda")

    # Create valid inputs inline
    B, multiplicity, N_row, N_col = 2, 2, 32, 32
    B_mul = B * multiplicity
    num_bins = 64
    max_dist = 32.0
    coord_dim = modifications.get("coord_dim", 3)

    inputs = {
        "pred_pde": torch.randn(B_mul, N_row, N_col, num_bins, device=device),
        "true_coords_row": torch.randn(B_mul, N_row, coord_dim, device=device),
        "true_coords_col": torch.randn(B_mul, N_col, coord_dim, device=device),
        "pred_coords_row": torch.randn(B_mul, N_row, coord_dim, device=device),
        "pred_coords_col": torch.randn(B_mul, N_col, coord_dim, device=device),
        "mask_row": torch.ones(B, N_row, device=device),
        "mask_col": torch.ones(B, N_col, device=device),
    }

    # Apply modifications
    for key, val in modifications.items():
        if key == "coord_dim":
            continue  # Already handled above
        elif key in ("mask_row", "mask_col"):
            inputs[key] = torch.ones(*val, device=device)
        elif key == "pred_pde":
            inputs[key] = torch.randn(*val, device=device)
        elif key in ("true_coords_row", "true_coords_col", "pred_coords_row", "pred_coords_col"):
            inputs[key] = torch.randn(*val, device=device)
        elif key.endswith("_requires_grad"):
            # Handle requires_grad modifications
            tensor_key = key.replace("_requires_grad", "")
            inputs[tensor_key] = inputs[tensor_key].requires_grad_(val)

    if expected_error is None:
        loss_num, mask_denom = cdist_pde(
            **inputs,
            multiplicity=multiplicity,
            num_bins=num_bins,
            max_dist=max_dist,
        )
        # Kernel now outputs fully summed [B_mul] instead of [B_mul, N_row]
        assert loss_num.shape == (B_mul,)
        assert mask_denom.shape == (B_mul,)
    else:
        with pytest.raises(ValueError, match=expected_error):
            cdist_pde(
                **inputs,
                multiplicity=multiplicity,
                num_bins=num_bins,
                max_dist=max_dist,
            )


def assert_no_register_spilling(path_to_ptx_file: Path):
    """Check that a PTX file shows no register spilling."""
    ptx_code = path_to_ptx_file.read_text()

    # Get the target architecture
    sm_arch_match = re.search(r"\.target (sm_\w+)", ptx_code)
    if not sm_arch_match:
        raise RuntimeError(f"No .target directive found in {path_to_ptx_file}")
    sm_arch = sm_arch_match.group(1)

    # Run ptxas
    ptxas_path = os.environ.get("TRITON_PTXAS_PATH", "ptxas")
    cmd = [ptxas_path, "-v", f"--gpu-name={sm_arch}", str(path_to_ptx_file)]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        output = result.stderr
        if "0 bytes spill stores, 0 bytes spill loads" not in output:
            raise RuntimeError(f"Register spilling detected in {path_to_ptx_file}:\n{output}")
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"ptxas failed with error:\n{e.stderr}") from e


@pytest.mark.skipif(not torch.cuda.is_available(), reason="CUDA not available")
def test_no_register_spilling(tmp_path, monkeypatch):
    """Test that kernels don't spill registers."""
    # Setup env for dumping PTX
    monkeypatch.setenv("TRITON_KERNEL_DUMP", "1")
    monkeypatch.setenv("TRITON_DUMP_DIR", str(tmp_path))
    monkeypatch.setenv("TRITON_ALWAYS_COMPILE", "1")
    monkeypatch.setenv("TRITON_CACHE_DIR", str(tmp_path / "cache"))
    monkeypatch.setenv("TRITON_PTXAS_PATH", os.environ.get("TRITON_PTXAS_PATH", "ptxas"))

    device = torch.device("cuda")
    # Use small problem size (like cdist_lddt test) to ensure quick compilation
    B_mul, B = 16, 1
    N = 100
    num_bins = 64

    pred_pde = torch.randn(B_mul, N, N, num_bins, device=device, requires_grad=True)
    true_coords = torch.randn(B_mul, N, 3, device=device)
    pred_coords = torch.randn(B_mul, N, 3, device=device)
    mask = torch.ones(B, N, device=device)

    # Run forward kernel
    loss_num, mask_denom = cdist_pde(
        pred_pde,
        true_coords,
        true_coords,
        pred_coords,
        pred_coords,
        mask,
        mask,
        multiplicity=B_mul // B,
        num_bins=num_bins,
    )

    # Run backward kernel
    loss_num.sum().backward()

    # Check PTX files for forward kernel
    fwd_ptx_files = list(tmp_path.glob("**/_cdist_pde_fwd_kernel.ptx"))
    if not fwd_ptx_files:
        raise RuntimeError(f"No forward kernel PTX file found in {tmp_path}")
    assert_no_register_spilling(fwd_ptx_files[0])

    # Check PTX files for backward kernel
    bwd_ptx_files = list(tmp_path.glob("**/_cdist_pde_bwd_kernel.ptx"))
    if not bwd_ptx_files:
        raise RuntimeError(f"No backward kernel PTX file found in {tmp_path}")
    assert_no_register_spilling(bwd_ptx_files[0])
