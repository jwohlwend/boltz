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

from boltz.distributed.model.loss.validation import clash_score
from boltz.distributed.model.loss.validation import factored_lddt_loss as triton_factored_lddt_loss
from boltz.model.loss.validation import factored_lddt_loss
from boltz.testing.utils import random_features


@pytest.mark.skipif(not torch.cuda.is_available(), reason="CUDA not available")
@pytest.mark.parametrize("use_cardinality_weighted", [True, False])
@pytest.mark.parametrize("repeat_atom_mask", [True, False])
def test_factored_lddt_loss_cdist_consistency(use_cardinality_weighted, repeat_atom_mask):
    device = torch.device("cuda")
    rng = torch.Generator(device=device)
    rng.manual_seed(0)
    rng_features = torch.Generator(device=device)
    rng_features.manual_seed(0)

    batch = 2
    num_tokens = 256
    num_atoms = num_tokens * 9
    multiplicity = 5

    true_atom_coords = torch.rand(batch, num_atoms, 3, generator=rng, device=device)
    pred_atom_coords = torch.rand(batch, num_atoms, 3, generator=rng, device=device)

    atom_mask_base = torch.randint(0, 2, (batch, num_atoms), generator=rng, device=device, dtype=torch.float32)

    true_atom_coords = true_atom_coords.repeat_interleave(multiplicity, 0)
    pred_atom_coords = pred_atom_coords.repeat_interleave(multiplicity, 0)
    atom_mask = atom_mask_base.repeat_interleave(multiplicity, 0) if repeat_atom_mask else atom_mask_base

    feats = random_features(
        size_batch=batch,
        n_tokens=num_tokens,
        n_atoms=num_atoms,
        n_msa=1,
        atom_counts_per_token_range=(1, 9),
        device=device,
        float_value_range=(0.0, 1.0),
        selected_keys=["atom_to_token", "asym_id", "mol_type"],
        rng=rng_features,
    )

    feats = {
        "atom_to_token": feats["atom_to_token"],
        "mol_type": feats["mol_type"],
        "asym_id": feats["asym_id"],
    }

    atom_mask_ref = atom_mask_base.repeat_interleave(multiplicity, 0)
    ref_lddt, ref_total = factored_lddt_loss(
        true_atom_coords=true_atom_coords,
        pred_atom_coords=pred_atom_coords,
        feats=feats,
        atom_mask=atom_mask_ref,
        multiplicity=multiplicity,
        cardinality_weighted=use_cardinality_weighted,
    )
    triton_lddt, triton_total = triton_factored_lddt_loss(
        true_atom_coords=true_atom_coords,
        pred_atom_coords=pred_atom_coords,
        feats=feats,
        atom_mask=atom_mask,
        multiplicity=multiplicity,
        cardinality_weighted=use_cardinality_weighted,
    )

    saw_zero_total = False
    for key in ref_lddt:
        torch.testing.assert_close(triton_lddt[key], ref_lddt[key])
        torch.testing.assert_close(triton_total[key], ref_total[key])
        zero_total_mask = ref_total[key] == 0
        if torch.any(zero_total_mask):
            saw_zero_total = True
            torch.testing.assert_close(ref_lddt[key][zero_total_mask], torch.ones_like(ref_lddt[key][zero_total_mask]))
            torch.testing.assert_close(
                triton_lddt[key][zero_total_mask], torch.ones_like(triton_lddt[key][zero_total_mask])
            )
    assert saw_zero_total, "Expected at least one modality to have zero total."


@pytest.mark.skipif(not torch.cuda.is_available(), reason="Test runs on triton kernel but CUDA is not available")
def test_clash_score_counts_and_fraction():
    device = torch.device("cuda")
    clash_cutoff = 2.0
    multiplicity = 2

    coords_repr = torch.tensor(
        [
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [5.0, 0.0, 0.0], [9.0, 0.0, 0.0]],
            [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [6.0, 0.0, 0.0], [9.0, 0.0, 0.0]],
        ],
        device=device,
        dtype=torch.float32,
    )
    token_pad_mask = torch.tensor([[True, True, True, False]], device=device)

    clash_atoms_count, clash_atoms_fraction = clash_score(
        coords_repr=coords_repr,
        token_pad_mask=token_pad_mask,
        multiplicity=multiplicity,
        clash_cutoff=clash_cutoff,
    )

    expected_count = torch.tensor([2, 0], device=device, dtype=clash_atoms_count.dtype)
    expected_fraction = torch.tensor([2.0 / 3.0, 0.0], device=device, dtype=clash_atoms_fraction.dtype)

    torch.testing.assert_close(clash_atoms_count, expected_count)
    torch.testing.assert_close(clash_atoms_fraction, expected_fraction)
