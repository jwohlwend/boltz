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

from boltz.distributed.model.loss.triton.cdist_lddt import cdist_lddt
from boltz.model.loss.confidencev2 import lddt_dist
from boltz.testing.utils import init_tensors_uniform


def lddt_dist_reference(
    pred_coords_row,
    pred_coords_col,
    true_coords_row,
    true_coords_col,
    resolved_mask_row,
    resolved_mask_col,
    atom_indices_row=None,  # [B, N_row]
    atom_indices_col=None,  # [B, N_col]
    cutoff=15.0,
    cutoff_col=None,  # [B, N_col] (optional, per-column per-batch cutoff values)
    do_mask_diagonal=True,
    per_atom=False,
    return_denom=False,
):
    """Reference implementation using torch.cdist and lddt_dist"""
    B_mul, N_row, _ = pred_coords_row.shape
    _, N_col, _ = pred_coords_col.shape
    B, _ = resolved_mask_row.shape
    multiplicity = B_mul // B
    device = pred_coords_row.device

    # Broadcast masks
    mask_row = resolved_mask_row.repeat_interleave(multiplicity, dim=0)
    mask_col = resolved_mask_col.repeat_interleave(multiplicity, dim=0)

    # Pair mask: (B_mul, N_row, N_col)
    pair_mask = mask_row.unsqueeze(-1) * mask_col.unsqueeze(-2)

    # Diagonal mask (conditional on do_mask_diagonal)
    if do_mask_diagonal:
        # Default to arange if indices not provided
        # atom_indices are [B, N_row] and [B, N_col], broadcast to B_mul
        if atom_indices_row is not None:
            idx_row = atom_indices_row.repeat_interleave(multiplicity, dim=0)  # [B_mul, N_row]
        else:
            idx_row = torch.arange(N_row, device=device).unsqueeze(0).expand(B_mul, -1)  # [B_mul, N_row]
        if atom_indices_col is not None:
            idx_col = atom_indices_col.repeat_interleave(multiplicity, dim=0)  # [B_mul, N_col]
        else:
            idx_col = torch.arange(N_col, device=device).unsqueeze(0).expand(B_mul, -1)  # [B_mul, N_col]
        # is_diagonal: [B_mul, N_row, N_col]
        is_diagonal = idx_row.unsqueeze(-1) == idx_col.unsqueeze(-2)
        pair_mask = pair_mask * (~is_diagonal)

    dmat_pred = torch.cdist(pred_coords_row, pred_coords_col)
    dmat_true = torch.cdist(true_coords_row, true_coords_col)

    # Compute cutoff tensor: if cutoff_col is provided, broadcast to [B_mul, N_row, N_col]
    if cutoff_col is not None:
        # cutoff_col is [B, N_col], broadcast to [B_mul, 1, N_col] then to [B_mul, N_row, N_col]
        cutoff_expanded = cutoff_col.repeat_interleave(multiplicity, dim=0).unsqueeze(1)  # [B_mul, 1, N_col]
        cutoff_tensor = cutoff_expanded.expand(-1, N_row, -1)  # [B_mul, N_row, N_col]
    else:
        cutoff_tensor = cutoff

    dists_to_score = (dmat_true < cutoff_tensor).float() * pair_mask

    # Use existing reference implementation
    # lddt_dist expects [B, N, N] inputs
    result = lddt_dist(dmat_pred, dmat_true, pair_mask, cutoff=cutoff_tensor, per_atom=per_atom)

    if per_atom:
        score, mask_no_match = result
        if return_denom:
            denom = torch.sum(dists_to_score, dim=-1)
            return score, mask_no_match, denom
        return score, mask_no_match

    score, _total = result
    if return_denom:
        denom = torch.sum(dists_to_score, dim=(-2, -1))
        return score, denom
    return score


@pytest.mark.skipif(not torch.cuda.is_available(), reason="CUDA not available")
# Order: last decorator = outermost loop (slowest varying)
# Put constexpr params (per_atom, do_mask_diagonal, return_unnormalized_score) last to minimize Triton recompilations
@pytest.mark.parametrize("B, multiplicity", [(2, 4)], ids=["B:2, M:4"])
@pytest.mark.parametrize("N_row", [32, 100], ids=lambda x: f"Nr:{x}")
@pytest.mark.parametrize("N_col", [32, 100], ids=lambda x: f"Nc:{x}")
@pytest.mark.parametrize("use_indices_row", [True, False], ids=lambda x: f"IdxR:{x}")
@pytest.mark.parametrize("use_indices_col", [True, False], ids=lambda x: f"IdxC:{x}")
@pytest.mark.parametrize("use_cutoff_col", [True, False], ids=lambda x: f"CutCol:{x}")
@pytest.mark.parametrize("do_mask_diagonal", [True, False], ids=lambda x: f"DiagMask:{x}")
@pytest.mark.parametrize("per_atom", [False, True], ids=lambda x: f"PerAtom:{x}")
@pytest.mark.parametrize("return_unnormalized_score", [False, True], ids=lambda x: f"Unnorm:{x}")
@pytest.mark.parametrize("return_denom", [True, False], ids=lambda x: f"Denom:{x}")
@pytest.mark.parametrize(
    "dtype",
    [torch.bfloat16, torch.float32, torch.float64],
    ids=lambda x: f"Dtype:{x}",
)
def test_cdist_lddt_correctness(
    B,
    multiplicity,
    N_row,
    N_col,
    use_indices_row,
    use_indices_col,
    use_cutoff_col,
    do_mask_diagonal,
    per_atom,
    return_unnormalized_score,
    return_denom,
    dtype,
):
    if torch.promote_types(dtype, torch.float32) != dtype:
        pytest.xfail(f"cdist_lddt requires at least float32 precision but got {dtype}")

    # return_unnormalized_score and per_atom are orthogonal options:
    # - per_atom controls output shape: [B_mul, N_row] vs [B_mul]
    # - return_unnormalized_score controls whether to return raw (out_num, out_denom) vs normalized score

    device = torch.device("cuda")
    B_mul = B * multiplicity

    # Value range for coordinate initialization (controls numerical stability)
    min_val_init = -0.5
    max_val_init = 0.5

    # Generate random data with controlled value range
    pred_coords_row = torch.empty(B_mul, N_row, 3, device=device, dtype=dtype)
    pred_coords_col = torch.empty(B_mul, N_col, 3, device=device, dtype=dtype)
    true_coords_row = torch.empty(B_mul, N_row, 3, device=device, dtype=dtype)
    true_coords_col = torch.empty(B_mul, N_col, 3, device=device, dtype=dtype)
    init_tensors_uniform(
        [pred_coords_row, pred_coords_col, true_coords_row, true_coords_col],
        low=min_val_init,
        high=max_val_init,
    )

    resolved_mask_row = torch.randint(0, 2, (B, N_row), device=device, dtype=dtype)
    resolved_mask_col = torch.randint(0, 2, (B, N_col), device=device, dtype=dtype)

    # Each index can independently be explicit or implicit (arange)
    # atom_indices have shape [B, N_row] and [B, N_col] (batch dimension B, not B_mul)
    # When the other dimension is larger, use a random non-duplicated subset from that
    # dimension to cover the rectangular subset case from minimum_lddt_symmetry_coords
    atom_indices_row = None
    atom_indices_col = None

    if use_indices_row:
        if N_col > N_row:
            # Random non-duplicated subset of arange(N_col) with size N_row, per batch
            # Each batch sample gets a different random permutation
            atom_indices_row = torch.stack(
                [torch.randperm(N_col, device=device)[:N_row].sort().values for _ in range(B)]
            )  # [B, N_row]
        else:
            # N_row >= N_col: use simple arange(N_row) for each batch
            atom_indices_row = torch.arange(N_row, device=device).unsqueeze(0).expand(B, -1).contiguous()  # [B, N_row]

    if use_indices_col:
        if N_row > N_col:
            # Random non-duplicated subset of arange(N_row) with size N_col
            # Each batch sample gets a different random permutation
            atom_indices_col = torch.stack(
                [torch.randperm(N_row, device=device)[:N_col].sort().values for _ in range(B)]
            )  # [B, N_col]
        else:
            # N_col >= N_row: use simple arange(N_col) for each batch
            atom_indices_col = torch.arange(N_col, device=device).unsqueeze(0).expand(B, -1).contiguous()  # [B, N_col]

    # Generate cutoff_col with shape [B, N_col] if use_cutoff_col is True
    # Coordinates are in [-0.5, 0.5], so max distance is ~sqrt(3) ≈ 1.73
    # Use cutoff values in [0.3, 1.2] range to meaningfully filter distances
    cutoff_col = None
    if use_cutoff_col:
        # Generate random cutoff values per batch, per column
        cutoff_col = torch.empty(B, N_col, device=device).uniform_(0.3, 1.2)

    if return_unnormalized_score and return_denom:
        pytest.skip("return_denom is invalid when return_unnormalized_score=True")

    # Reference
    ref_result = lddt_dist_reference(
        pred_coords_row,
        pred_coords_col,
        true_coords_row,
        true_coords_col,
        resolved_mask_row,
        resolved_mask_col,
        atom_indices_row,
        atom_indices_col,
        cutoff_col=cutoff_col,
        do_mask_diagonal=do_mask_diagonal,
        per_atom=per_atom,
        return_denom=return_denom,
    )

    # Triton
    triton_result = cdist_lddt(
        pred_coords_row,
        pred_coords_col,
        true_coords_row,
        true_coords_col,
        resolved_mask_row,
        resolved_mask_col,
        multiplicity,
        atom_indices_row=atom_indices_row,
        atom_indices_col=atom_indices_col,
        cutoff_col=cutoff_col,
        do_mask_diagonal=do_mask_diagonal,
        per_atom=per_atom,
        return_unnormalized_score=return_unnormalized_score,
        return_denom=return_denom,
    )

    eps = 1e-10
    if return_unnormalized_score:
        # return_unnormalized_score returns unnormalized scores before normalization
        # We verify by manually computing the normalized result and comparing to reference
        if per_atom:
            out_num, out_denom, mask_no_match_triton = triton_result
            norm = 1.0 / (eps + out_denom)
            computed_score = norm * (eps + out_num)
            score_ref, mask_no_match_ref = ref_result
            torch.testing.assert_close(computed_score, score_ref)
            torch.testing.assert_close(mask_no_match_triton.to(mask_no_match_ref.dtype), mask_no_match_ref)
        else:
            out_num, out_denom = triton_result
            score_ref = ref_result
            computed_score = out_num / (out_denom + eps)
            computed_score = torch.where(out_denom > 0, computed_score, torch.zeros_like(computed_score))
            torch.testing.assert_close(computed_score, score_ref)
    else:
        if per_atom and return_denom:
            score_ref, mask_no_match_ref, denom_ref = ref_result
            score_triton, mask_no_match_triton, denom_triton = triton_result
            torch.testing.assert_close(score_triton.to(score_ref.dtype), score_ref)
            torch.testing.assert_close(mask_no_match_triton.to(mask_no_match_ref.dtype), mask_no_match_ref)
            torch.testing.assert_close(denom_triton.to(denom_ref.dtype), denom_ref)
        elif per_atom and not return_denom:
            score_ref, mask_no_match_ref = ref_result
            score_triton, mask_no_match_triton = triton_result
            torch.testing.assert_close(score_triton.to(score_ref.dtype), score_ref)
            torch.testing.assert_close(mask_no_match_triton.to(mask_no_match_ref.dtype), mask_no_match_ref)
        elif not per_atom and return_denom:
            score_ref, denom_ref = ref_result
            score_triton, denom_triton = triton_result
            torch.testing.assert_close(score_triton.to(score_ref.dtype), score_ref)
            torch.testing.assert_close(denom_triton.to(denom_ref.dtype), denom_ref)
        else:  # not per_atom and not return_denom
            score_ref = ref_result
            score_triton = triton_result
            torch.testing.assert_close(score_triton.to(score_ref.dtype), score_ref)


@pytest.fixture(
    params=[
        # (modifications_dict, expected_error_pattern)
        ({}, None),  # valid inputs
        ({"mask_row": (3, 32), "mask_col": (3, 48)}, "mask_row batch dimension"),  # Neither B nor B_mul
        ({"coord_dim": 4}, "Coordinate dimension must be 3"),
        ({"pred_coords_col": (4, 48, 3)}, "pred_coords_col shape"),
        ({"true_coords_row": (8, 64, 3)}, "true_coords_row shape"),
        ({"true_coords_col": (8, 64, 3)}, "true_coords_col shape"),
        ({"mask_row": (2, 64)}, "mask_row N dimension"),
        ({"mask_col": (2, 64)}, "mask_col N dimension"),
        (
            {"mask_row": (2, 32), "mask_col": (4, 48)},
            "mask_col batch dimension",
        ),  # mask_col batch is 4, neither B=2 nor B_mul=8
        ({"atom_indices_row": (2, 64)}, "atom_indices_row shape"),
        ({"atom_indices_col": (2, 64)}, "atom_indices_col shape"),
        ({"atom_indices_row": None, "atom_indices_col": None}, None),  # None indices valid
        (
            {"return_unnormalized_score": True, "return_denom": True},
            "return_denom is not valid when return_unnormalized_score=True",
        ),
    ],
    ids=[
        "valid_inputs",
        "mask_batch_neither_B_nor_B_mul",
        "coord_dim_not_3",
        "pred_coords_col_wrong_batch",
        "true_coords_row_wrong_n_row",
        "true_coords_col_wrong_n_col",
        "mask_row_wrong_n_row",
        "mask_col_wrong_n_col",
        "mask_col_batch_neither_B_nor_B_mul",
        "atom_indices_row_wrong_shape",
        "atom_indices_col_wrong_shape",
        "none_indices_valid",
        "return_unnormalized_score_and_denom_invalid",
    ],
)
def validation_case(request):
    """Fixture providing (modifications, expected_error) for validation tests"""
    return request.param


@pytest.mark.skipif(not torch.cuda.is_available(), reason="CUDA not available")
def test_cdist_lddt_validation(validation_case):
    """Test input validation for cdist_lddt()"""
    modifications, expected_error = validation_case
    device = torch.device("cuda")

    # Create valid inputs inline
    B, multiplicity, N_row, N_col = 2, 4, 32, 48
    B_mul = B * multiplicity
    coord_dim = modifications.get("coord_dim", 3)

    rng = torch.Generator(device=device)
    rng.manual_seed(0)
    inputs = {
        "pred_coords_row": torch.randn(B_mul, N_row, coord_dim, device=device, generator=rng),
        "pred_coords_col": torch.randn(B_mul, N_col, coord_dim, device=device, generator=rng),
        "true_coords_row": torch.randn(B_mul, N_row, coord_dim, device=device, generator=rng),
        "true_coords_col": torch.randn(B_mul, N_col, coord_dim, device=device, generator=rng),
        "mask_row": torch.ones(B, N_row, device=device),
        "mask_col": torch.ones(B, N_col, device=device),
        "atom_indices_row": torch.arange(N_row, device=device).unsqueeze(0).expand(B, -1).contiguous(),
        "atom_indices_col": torch.arange(N_col, device=device).unsqueeze(0).expand(B, -1).contiguous(),
        "return_unnormalized_score": False,
        "return_denom": False,
    }

    # Apply modifications
    for key, val in modifications.items():
        if key == "coord_dim":
            continue  # Already handled above
        elif key in ("mask_row", "mask_col"):
            inputs[key] = torch.ones(*val, device=device)
        elif key in ("pred_coords_row", "pred_coords_col", "true_coords_row", "true_coords_col"):
            inputs[key] = torch.randn(*val, device=device)
        elif key in ("atom_indices_row", "atom_indices_col"):
            if val is None:
                inputs[key] = None
            elif isinstance(val, tuple):
                # val is (B, N) shape
                inputs[key] = torch.arange(val[1], device=device).unsqueeze(0).expand(val[0], -1).contiguous()
            else:
                # legacy: val is just N
                inputs[key] = torch.arange(val, device=device).unsqueeze(0).expand(B, -1).contiguous()
        elif key in ("return_unnormalized_score", "return_denom"):
            inputs[key] = val

    if expected_error is None:
        result = cdist_lddt(**inputs, multiplicity=multiplicity)
        assert result.shape == (inputs["pred_coords_row"].shape[0],)
    else:
        with pytest.raises(ValueError, match=expected_error):
            cdist_lddt(**inputs, multiplicity=multiplicity)


def assert_no_register_spilling(path_to_ptx_file: Path):
    ptx_code = path_to_ptx_file.read_text()

    # get the ".target sm_{arch}a" directive from the ptx code
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
    # Setup env for dumping PTX
    monkeypatch.setenv("TRITON_KERNEL_DUMP", "1")
    monkeypatch.setenv("TRITON_DUMP_DIR", str(tmp_path))
    monkeypatch.setenv("TRITON_ALWAYS_COMPILE", "1")
    monkeypatch.setenv("TRITON_CACHE_DIR", str(tmp_path / "cache"))
    monkeypatch.setenv("TRITON_PTXAS_PATH", os.environ.get("TRITON_PTXAS_PATH", "ptxas"))

    device = torch.device("cuda")
    B_mul, B, N = 16, 1, 100

    pred_coords = torch.randn(B_mul, N, 3, device=device)
    true_coords = torch.randn(B_mul, N, 3, device=device)
    mask = torch.ones(B, N, device=device)

    # Run kernel (multiplicity = B_mul // B = 16 // 1 = 16)
    cdist_lddt(pred_coords, pred_coords, true_coords, true_coords, mask, mask, multiplicity=B_mul // B)

    # Check PTX
    ptx_files = list(tmp_path.glob("**/_cdist_lddt_kernel.ptx"))
    if not ptx_files:
        raise RuntimeError(f"No PTX file found in {tmp_path}")

    assert_no_register_spilling(ptx_files[0])
