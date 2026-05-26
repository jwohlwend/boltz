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

from boltz.data import const
from boltz.distributed.model.loss.validation import factored_token_lddt_dist_loss_triton
from boltz.model.loss.validation import factored_token_lddt_dist_loss
from boltz.testing.utils import random_features


@pytest.mark.skipif(not torch.cuda.is_available(), reason="CUDA not available")
@pytest.mark.parametrize("use_cardinality_weighted", [True, False])
@pytest.mark.parametrize("pass_pred_d", [True, False])
@pytest.mark.parametrize("pass_true_d", [True, False])
def test_factored_token_lddt_dist_loss_triton_parity(use_cardinality_weighted, pass_pred_d, pass_true_d):
    """Verify factored_token_lddt_dist_loss_triton matches the serial factored_token_lddt_dist_loss."""
    device = torch.device("cuda")
    rng = torch.Generator(device=device)
    rng.manual_seed(42)
    rng_features = torch.Generator(device=device)
    rng_features.manual_seed(42)

    batch = 2
    num_tokens = 256
    multiplicity = 1

    pred_token_coords = torch.rand(batch, num_tokens, 3, generator=rng, device=device)
    true_token_coords = torch.rand(batch, num_tokens, 3, generator=rng, device=device)

    feats = random_features(
        size_batch=batch,
        n_tokens=num_tokens,
        n_atoms=num_tokens * 3,
        n_msa=1,
        atom_counts_per_token_range=(1, 3),
        device=device,
        float_value_range=(0.0, 1.0),
        selected_keys=["mol_type", "asym_id", "token_disto_mask"],
        rng=rng_features,
    )

    mol_type = feats["mol_type"].long()
    token_disto_mask = feats["token_disto_mask"].float()
    asym_id = feats["asym_id"].long()

    pred_d = torch.cdist(pred_token_coords, pred_token_coords)
    true_d = torch.cdist(true_token_coords, true_token_coords)

    serial_feats = {
        "mol_type": mol_type,
        "token_disto_mask": token_disto_mask,
        "asym_id": asym_id,
    }
    ref_lddt, ref_total = factored_token_lddt_dist_loss(
        true_d=true_d,
        pred_d=pred_d,
        feats=serial_feats,
        cardinality_weighted=use_cardinality_weighted,
    )

    triton_lddt, triton_total = factored_token_lddt_dist_loss_triton(
        pred_token_coords=pred_token_coords,
        true_token_coords=true_token_coords,
        mol_type=mol_type,
        token_disto_mask=token_disto_mask,
        asym_id=asym_id,
        multiplicity=multiplicity,
        cardinality_weighted=use_cardinality_weighted,
        pred_d=pred_d if pass_pred_d else None,
        true_d=true_d if pass_true_d else None,
    )

    for key in ref_lddt:
        torch.testing.assert_close(triton_lddt[key], ref_lddt[key], atol=1e-5, rtol=1e-5)
        torch.testing.assert_close(triton_total[key], ref_total[key], atol=1e-5, rtol=1e-5)
        zero_total_mask = ref_total[key] == 0
        if torch.any(zero_total_mask):
            torch.testing.assert_close(ref_lddt[key][zero_total_mask], torch.ones_like(ref_lddt[key][zero_total_mask]))
            torch.testing.assert_close(
                triton_lddt[key][zero_total_mask], torch.ones_like(triton_lddt[key][zero_total_mask])
            )


@pytest.mark.skipif(not torch.cuda.is_available(), reason="CUDA not available")
def test_factored_token_lddt_dist_loss_triton_multiplicity():
    """Verify triton token LDDT works with multiplicity > 1."""
    device = torch.device("cuda")
    rng = torch.Generator(device=device)
    rng.manual_seed(99)
    rng_features = torch.Generator(device=device)
    rng_features.manual_seed(99)

    batch = 2
    num_tokens = 128
    multiplicity = 3

    pred_base = torch.rand(batch, num_tokens, 3, generator=rng, device=device)
    true_base = torch.rand(batch, num_tokens, 3, generator=rng, device=device)

    pred_token_coords = pred_base.repeat_interleave(multiplicity, 0)
    true_token_coords = true_base.repeat_interleave(multiplicity, 0)

    feats = random_features(
        size_batch=batch,
        n_tokens=num_tokens,
        n_atoms=num_tokens * 3,
        n_msa=1,
        atom_counts_per_token_range=(1, 3),
        device=device,
        float_value_range=(0.0, 1.0),
        selected_keys=["mol_type", "asym_id", "token_disto_mask"],
        rng=rng_features,
    )

    mol_type = feats["mol_type"].long()
    token_disto_mask = feats["token_disto_mask"].float()
    asym_id = feats["asym_id"].long()

    triton_lddt, triton_total = factored_token_lddt_dist_loss_triton(
        pred_token_coords=pred_token_coords,
        true_token_coords=true_token_coords,
        mol_type=mol_type,
        token_disto_mask=token_disto_mask,
        asym_id=asym_id,
        multiplicity=multiplicity,
    )

    pred_d = torch.cdist(pred_base, pred_base)
    true_d = torch.cdist(true_base, true_base)

    serial_feats = {
        "mol_type": mol_type,
        "token_disto_mask": token_disto_mask,
        "asym_id": asym_id,
    }
    ref_lddt, ref_total = factored_token_lddt_dist_loss(
        true_d=true_d,
        pred_d=pred_d,
        feats=serial_feats,
    )

    for key in ref_lddt:
        triton_per_sample = triton_lddt[key].reshape(batch, multiplicity)
        for m in range(multiplicity):
            torch.testing.assert_close(triton_per_sample[:, m], ref_lddt[key], atol=1e-5, rtol=1e-5)


@pytest.mark.skipif(not torch.cuda.is_available(), reason="CUDA not available")
def test_factored_token_lddt_dist_loss_triton_zero_total_defaults_to_one():
    """When a modality has no valid pairs, score must be 1.0 (not NaN or 0)."""
    device = torch.device("cuda")
    B = 1
    N = 8

    pred = torch.rand(B, N, 3, device=device)
    true = torch.rand(B, N, 3, device=device)

    mol_type = torch.full((B, N), const.chain_type_ids["PROTEIN"], dtype=torch.long, device=device)
    token_disto_mask = torch.ones(B, N, device=device)
    asym_id = torch.zeros(B, N, dtype=torch.long, device=device)

    lddt_dict, total_dict = factored_token_lddt_dist_loss_triton(
        pred_token_coords=pred,
        true_token_coords=true,
        mol_type=mol_type,
        token_disto_mask=token_disto_mask,
        asym_id=asym_id,
    )

    empty_modalities = [
        "dna_protein",
        "rna_protein",
        "dna_ligand",
        "rna_ligand",
        "ligand_protein",
        "intra_ligand",
        "intra_dna",
        "intra_rna",
        "protein_protein",
    ]
    for key in empty_modalities:
        assert total_dict[key].item() == 0, f"{key} should have zero total"
        assert lddt_dict[key].item() == 1.0, f"{key} should default to 1.0 when no pairs exist"
