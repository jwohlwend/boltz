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

"""Tests for compute_plddt_mae_triton.

Verifies that compute_plddt_mae_triton produces scalar results matching the
serial compute_plddt_mae from boltz.model.loss.validation on a single GPU.
"""

from __future__ import annotations

import pytest
import torch

from boltz.distributed.model.loss.validation import compute_plddt_mae_triton
from boltz.model.loss.validation import compute_plddt_mae
from boltz.testing.utils import random_features

EXPECTED_PLDDT_KEYS = {"protein", "ligand", "dna", "rna"}


def _generate_plddt_test_data(multiplicity, resolved_mask_mode, rng):
    """Generate input tensors for plddt MAE tests with batch_size=1."""
    n_tokens = 20
    n_atoms = n_tokens * 20
    feats_host = random_features(
        size_batch=1,
        n_tokens=n_tokens,
        n_atoms=n_atoms,
        n_msa=1,
        atom_counts_per_token_range=(1, 20),
        device=torch.device("cpu"),
        float_value_range=(-5.0, 5.0),
        selected_keys=[
            "token_to_rep_atom",
            "r_set_to_rep_atom",
            "atom_to_token",
            "mol_type",
            "atom_counts_per_token",
        ],
        rng=rng,
    )
    N_atom_actual = feats_host["token_to_rep_atom"].shape[2]
    B_mul = multiplicity
    pred_coords = torch.randn(B_mul, N_atom_actual, 3, dtype=torch.float32)
    true_coords = torch.randn(B_mul, N_atom_actual, 3, dtype=torch.float32)
    pred_lddt = torch.rand(B_mul, n_tokens, dtype=torch.float32)
    if resolved_mask_mode == "ones":
        resolved_mask = torch.ones(B_mul, N_atom_actual, dtype=torch.float32)
    elif resolved_mask_mode == "zeros":
        resolved_mask = torch.zeros(B_mul, N_atom_actual, dtype=torch.float32)
    else:
        resolved_mask = torch.randint(0, 2, (B_mul, N_atom_actual)).float()
    return feats_host, pred_coords, true_coords, pred_lddt, resolved_mask


def _run_comparison(multiplicity, resolved_mask_mode, seed):
    """Run triton vs serial comparison and return both results."""
    torch.manual_seed(seed)
    rng = torch.Generator(device="cpu")
    rng.manual_seed(seed)

    feats_host, pred_coords, true_coords, pred_lddt, resolved_mask = _generate_plddt_test_data(
        multiplicity, resolved_mask_mode, rng
    )

    feats = {
        "r_set_to_rep_atom": feats_host["r_set_to_rep_atom"],
        "mol_type": feats_host["mol_type"],
        "atom_to_token": feats_host["atom_to_token"],
        "token_to_rep_atom": feats_host["token_to_rep_atom"],
    }

    ref_mae, ref_total = compute_plddt_mae(
        pred_atom_coords=pred_coords,
        feats=feats,
        true_atom_coords=true_coords,
        pred_lddt=pred_lddt,
        true_coords_resolved_mask=resolved_mask,
        multiplicity=multiplicity,
    )

    device = torch.device("cuda")
    feats_cuda = {k: v.to(device) for k, v in feats.items()}

    triton_mae, triton_total = compute_plddt_mae_triton(
        pred_atom_coords=pred_coords.to(device),
        feats=feats_cuda,
        true_atom_coords=true_coords.to(device),
        pred_lddt=pred_lddt.to(device),
        true_coords_resolved_mask=resolved_mask.to(device),
        multiplicity=multiplicity,
    )

    return ref_mae, ref_total, triton_mae, triton_total


@pytest.mark.skipif(not torch.cuda.is_available(), reason="CUDA required")
@pytest.mark.parametrize("multiplicity", [1, 2])
@pytest.mark.parametrize("resolved_mask_mode", ["ones", "zeros", "partial"])
def test_compute_plddt_mae_triton(multiplicity, resolved_mask_mode):
    """compute_plddt_mae_triton must match serial compute_plddt_mae."""
    ref_mae, ref_total, triton_mae, triton_total = _run_comparison(
        multiplicity=multiplicity,
        resolved_mask_mode=resolved_mask_mode,
        seed=42,
    )

    assert set(triton_mae.keys()) == EXPECTED_PLDDT_KEYS
    assert set(triton_total.keys()) == EXPECTED_PLDDT_KEYS

    for key in EXPECTED_PLDDT_KEYS:
        if resolved_mask_mode == "zeros":
            assert ref_total[key] == 0.0, f"Expected zero total for '{key}'"
        torch.testing.assert_close(triton_mae[key].cpu(), ref_mae[key])
        torch.testing.assert_close(triton_total[key].cpu(), ref_total[key])
