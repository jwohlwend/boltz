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

"""Tests for featurizerv2.py (inference featurizer) bug fixes.

Covers:
  1. `visited` set construction: taxonomy-assigned sequences must be excluded
     from the unpaired `available` pool.
  2. In-place mutation: construct_paired_msa must not mutate caller's MSA data.
  3. Deletion indexing: the deletion extraction loop must slice from the
     ORIGINAL full array on every iteration, not from the previous slice.
"""

from types import SimpleNamespace

import numpy as np
import pytest

from boltz.data import const
from boltz.data.feature.featurizerv2 import construct_paired_msa
from boltz.data.types import (
    MSA,
    Chain,
    MSADeletion,
    MSAResidue,
    MSASequence,
    Residue,
    StructureV2,
    Token,
)

# ---------------------------------------------------------------------------
# Helpers -- same synthetic builders as the train tests, but using Token dtype
# (inference featurizer operates on Tokenized which uses Token, not TokenV2)
# ---------------------------------------------------------------------------


def _make_chain(chain_id, res_idx, res_num):
    """Create a single-element Chain structured array.

    Parameters
    ----------
    chain_id : int
        Value for the ``asym_id`` field.
    res_idx : int
        Starting index of this chain's residues in the global residue array.
    res_num : int
        Number of residues in the chain.
    """
    arr = np.zeros(1, dtype=Chain)
    arr[0]["asym_id"] = chain_id
    arr[0]["res_idx"] = res_idx
    arr[0]["res_num"] = res_num
    return arr


def _make_msa(residue_types, taxonomies, deletions_per_seq=None):
    """Build an MSA from explicit residue types, taxonomies, and deletions.

    Parameters
    ----------
    residue_types : list[list[int]]
        Outer = sequences, inner = residue token IDs.
    taxonomies : list[int]
        Taxonomy ID per sequence (-1 for none).
    deletions_per_seq : list[list[tuple[int, int]]], optional
        Per-sequence list of (res_idx, count) deletion entries.
    """
    if deletions_per_seq is None:
        deletions_per_seq = [[] for _ in residue_types]

    all_residues, all_deletions, sequences = [], [], []
    for seq_idx, (res_types, taxon, dels) in enumerate(zip(residue_types, taxonomies, deletions_per_seq)):
        res_start = len(all_residues)
        all_residues.extend(res_types)
        res_end = len(all_residues)

        del_start = len(all_deletions)
        all_deletions.extend(dels)
        del_end = len(all_deletions)

        sequences.append((seq_idx, taxon, res_start, res_end, del_start, del_end))

    return MSA(
        residues=np.array(all_residues, dtype=MSAResidue),
        deletions=np.array(all_deletions, dtype=MSADeletion),
        sequences=np.array(sequences, dtype=MSASequence),
    )


def _make_tokens(asym_ids, res_idxs):
    """Create a Token array with ``asym_id``, ``res_idx``, and ``token_idx`` populated.

    Parameters
    ----------
    asym_ids : list[int]
        Per-token chain identifier.
    res_idxs : list[int]
        Per-token residue index (0-based within its chain).
    """
    n = len(asym_ids)
    tokens = np.zeros(n, dtype=Token)
    tokens["asym_id"] = asym_ids
    tokens["res_idx"] = res_idxs
    tokens["token_idx"] = np.arange(n)
    return tokens


def _make_data(chain_specs, msas):
    """Build a minimal data object accepted by construct_paired_msa.

    Parameters
    ----------
    chain_specs : list[list[int]]
        Per-chain structure residue types (must match MSA query, seq 0).
    msas : dict[int, MSA]
    """
    chains_list, residues_list, asym_ids, res_idxs = [], [], [], []
    offset = 0
    for chain_id, res_types in enumerate(chain_specs):
        n_res = len(res_types)
        chains_list.append(_make_chain(chain_id, res_idx=offset, res_num=n_res))
        res = np.zeros(n_res, dtype=Residue)
        res["res_type"] = res_types
        for i in range(n_res):
            res[i]["res_idx"] = i
        residues_list.append(res)
        for r in range(n_res):
            asym_ids.append(chain_id)
            res_idxs.append(r)
        offset += n_res

    chains = np.concatenate(chains_list)
    residues = np.concatenate(residues_list)
    tokens = _make_tokens(asym_ids, res_idxs)
    structure = SimpleNamespace(chains=chains, residues=residues)
    record = SimpleNamespace(id="test_sample")
    return SimpleNamespace(tokens=tokens, structure=structure, msa=msas, record=record)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture()
def two_chain_taxonomy_data():
    """2 chains x 3 residues, each with 3 MSA sequences.

    seq 0 = query, seq 1 = Human (9606), seq 2 = unique per chain.
    """
    chain0_query = [1, 2, 3]
    chain1_query = [6, 7, 8]

    msa_chain0 = _make_msa(
        residue_types=[chain0_query, [1, 2, 4], [5, 2, 3]],
        taxonomies=[-1, 9606, 7227],
    )
    msa_chain1 = _make_msa(
        residue_types=[chain1_query, [6, 7, 9], [10, 7, 8]],
        taxonomies=[-1, 9606, 7955],
    )
    data = _make_data(
        chain_specs=[chain0_query, chain1_query],
        msas={0: msa_chain0, 1: msa_chain1},
    )
    return data


@pytest.fixture()
def met_unk_mismatch_data():
    """Chain 0: structure [MET, ALA, GLY], MSA query [UNK, ALA, GLY]."""
    met_id = const.token_ids["MET"]
    unk_id = const.token_ids["UNK"]
    ala_id = const.token_ids["ALA"]
    gly_id = const.token_ids["GLY"]

    msa_chain0 = _make_msa(
        residue_types=[[unk_id, ala_id, gly_id]],
        taxonomies=[-1],
    )
    data = _make_data(
        chain_specs=[[met_id, ala_id, gly_id]],
        msas={0: msa_chain0},
    )
    return data, msa_chain0


@pytest.fixture()
def data_with_deletions():
    """1 chain, 3 sequences: query (no deletions), seq 1 (1 del), seq 2 (2 dels)."""
    query = [1, 2, 3, 4, 5]
    msa = _make_msa(
        residue_types=[
            query,
            [1, 2, 4, 4, 5],
            [2, 2, 3, 4, 6],
        ],
        taxonomies=[-1, -1, -1],
        deletions_per_seq=[
            [],
            [(2, 3)],
            [(0, 1), (4, 5)],
        ],
    )
    return _make_data(chain_specs=[query], msas={0: msa})


# ---------------------------------------------------------------------------
# Test 1: visited set -- taxonomy-assigned seqs excluded from available pool
#
# Bug: the old comprehension {(c, s) for c, items in taxonomy_map for s in items}
# produced {(taxon, (chain_id, seq_idx)), ...}.  The downstream membership
# check (chain_id, seq_idx) not in visited never matched, so every sequence
# leaked into the unpaired pool.
# ---------------------------------------------------------------------------


def test_visited_unpaired_fill_excludes_taxonomy_assigned(two_chain_taxonomy_data):
    """Sequences assigned to taxonomy groups should NOT appear as
    unpaired fill in other rows."""
    gap = const.token_ids["-"]
    rng = np.random.default_rng(42)
    msa_data, _, _ = construct_paired_msa(two_chain_taxonomy_data, random=rng, max_seqs=100)

    n_chain0_tokens = 3
    n_rows = msa_data.shape[1]
    human_chain0 = (1, 2, 4)
    human_chain1 = (6, 7, 9)

    for row_idx in range(2, n_rows):
        chain0_row = tuple(msa_data[:n_chain0_tokens, row_idx].tolist())
        chain1_row = tuple(msa_data[n_chain0_tokens:, row_idx].tolist())
        is_gap_chain0 = all(v == gap for v in chain0_row)
        is_gap_chain1 = all(v == gap for v in chain1_row)
        if not is_gap_chain0:
            assert chain0_row != human_chain0, f"Row {row_idx}: chain 0 Human seq reused as unpaired fill"
        if not is_gap_chain1:
            assert chain1_row != human_chain1, f"Row {row_idx}: chain 1 Human seq reused as unpaired fill"


# ---------------------------------------------------------------------------
# Test 2: in-place mutation -- caller's MSA data must not be modified
#
# Bug: msa_residues was a direct view into data.msa[chain_id].residues.
# Modifying msa_residues[...]["res_type"] = ... mutated the caller's data.
# ---------------------------------------------------------------------------


def test_copy_on_write_original_msa_residues_unchanged(met_unk_mismatch_data):
    """After construct_paired_msa, the caller's MSA residue array must
    still contain the original UNK token, not the patched MET."""
    data, original_msa = met_unk_mismatch_data
    original_residues_before = original_msa.residues["res_type"].copy()

    rng = np.random.default_rng(42)
    construct_paired_msa(data, random=rng, max_seqs=10)

    np.testing.assert_array_equal(
        original_msa.residues["res_type"],
        original_residues_before,
        err_msg="construct_paired_msa mutated the caller's MSA residues in-place.",
    )


def test_copy_on_write_idempotent_on_double_call(met_unk_mismatch_data):
    """Calling construct_paired_msa twice on the same data must produce
    identical results, proving the first call didn't corrupt the input."""
    data, _ = met_unk_mismatch_data

    rng1 = np.random.default_rng(42)
    msa1, del1, paired1 = construct_paired_msa(data, random=rng1, max_seqs=10)

    rng2 = np.random.default_rng(42)
    msa2, del2, paired2 = construct_paired_msa(data, random=rng2, max_seqs=10)

    np.testing.assert_array_equal(msa1, msa2, err_msg="Double-call MSA mismatch")
    np.testing.assert_array_equal(del1, del2, err_msg="Double-call deletion mismatch")
    np.testing.assert_array_equal(paired1, paired2, err_msg="Double-call paired mismatch")


# ---------------------------------------------------------------------------
# Test 3: deletion indexing -- must slice from original array every iteration
#
# Bug (inference-specific): chain_deletions = chain_deletions[del_start:del_end]
# progressively shrank the array.  After the query (del_start=0, del_end=0),
# chain_deletions became empty; all subsequent deletions were silently lost.
# With the bug, n_nonzero == 0.  With the fix, n_nonzero == 3.
# ---------------------------------------------------------------------------


def test_deletion_all_present_in_output(data_with_deletions):
    """Every deletion entry must appear in the output tensor."""
    rng = np.random.default_rng(42)
    _, del_data, _ = construct_paired_msa(data_with_deletions, random=rng, max_seqs=10)

    assert del_data[2, 1].item() == 3, (
        f"Seq 1 deletion at res_idx=2: expected 3, got {del_data[2, 1].item()}. "
        "This fails if chain_deletions is progressively shrunk."
    )
    assert del_data[0, 2].item() == 1, f"Seq 2 deletion at res_idx=0: expected 1, got {del_data[0, 2].item()}"
    assert del_data[4, 2].item() == 5, f"Seq 2 deletion at res_idx=4: expected 5, got {del_data[4, 2].item()}"


def test_deletion_query_row_has_no_deletions(data_with_deletions):
    """The query row (row 0) should have zero deletions everywhere."""
    rng = np.random.default_rng(42)
    _, del_data, _ = construct_paired_msa(data_with_deletions, random=rng, max_seqs=10)

    query_dels = del_data[:, 0]
    assert (query_dels == 0).all(), f"Query row should have zero deletions, got {query_dels}"


def test_deletion_nonzero_count_matches_expected(data_with_deletions):
    """With the bug, n_nonzero would be 0 (all deletions lost).
    With the fix, n_nonzero should be exactly 3."""
    rng = np.random.default_rng(42)
    _, del_data, _ = construct_paired_msa(data_with_deletions, random=rng, max_seqs=10)

    n_nonzero = (del_data != 0).sum().item()
    assert n_nonzero == 3, (
        f"Expected 3 nonzero deletion entries, got {n_nonzero}. "
        "A value of 0 indicates the progressive-shrink bug is still present."
    )


# ---------------------------------------------------------------------------
# Real-data integration tests
#
# Load the 8ayv homodimer (2 protein chains sharing one MSA with 4611 seqs,
# 2309 with taxonomy, 716 with deletions) and run construct_paired_msa on
# it to verify the fixes work on production-format data.
# ---------------------------------------------------------------------------


def _load_real_sample(sample_dir):
    """Build a data object from a processed sample directory.

    Uses Token dtype (inference featurizer) instead of TokenV2.
    """
    import json

    manifest = json.loads((sample_dir / "manifest.json").read_text())
    rec = manifest["records"][0]

    structure = StructureV2.load(sample_dir / "structures" / f"{rec['id']}.npz")

    msa_dir = sample_dir / "msa"
    msas = {}
    for chain_info in rec["chains"]:
        msa_id = chain_info.get("msa_id", -1)
        chain_id = chain_info["chain_id"]
        if msa_id != -1:
            msa_path = msa_dir / f"{msa_id}.npz"
            if msa_path.exists():
                msas[chain_id] = MSA.load(msa_path)

    asym_ids, res_idxs = [], []
    for chain in structure.chains:
        cid = int(chain["asym_id"])
        if cid not in msas:
            continue
        for r in range(int(chain["res_num"])):
            asym_ids.append(cid)
            res_idxs.append(r)

    tokens = _make_tokens(asym_ids, res_idxs)
    record = SimpleNamespace(id=rec["id"])
    data = SimpleNamespace(tokens=tokens, structure=structure, msa=msas, record=record)
    return data, msas


@pytest.fixture()
def real_8ayv_data(test_cp_training_base_data_dir_boltz2):
    """Load the 8ayv homodimer from the real test data cache."""
    return _load_real_sample(test_cp_training_base_data_dir_boltz2 / "processed_8ayv")


def test_real_data_no_mutation(real_8ayv_data):
    """MSA residue arrays must not be modified in-place by construct_paired_msa."""
    data, original_msas = real_8ayv_data

    snapshots = {cid: msa.residues["res_type"].copy() for cid, msa in original_msas.items()}

    rng = np.random.default_rng(42)
    construct_paired_msa(data, random=rng, max_seqs=512)

    for cid, before in snapshots.items():
        np.testing.assert_array_equal(
            original_msas[cid].residues["res_type"],
            before,
            err_msg=f"Chain {cid}: MSA residues mutated in-place",
        )


def test_real_data_deletions_preserved(real_8ayv_data):
    """The output deletion matrix must contain nonzero entries when the
    input MSA has sequences with deletions (8ayv has 716)."""
    data, _ = real_8ayv_data

    rng = np.random.default_rng(42)
    _, del_data, _ = construct_paired_msa(data, random=rng, max_seqs=512)

    n_nonzero = (del_data != 0).sum().item()
    assert n_nonzero > 0, (
        "del_data is all zeros despite 716 input sequences having deletions. "
        "This indicates the deletion extraction loop is broken."
    )


def test_real_data_idempotent(real_8ayv_data):
    """Two calls with the same seed must produce identical output."""
    data, _ = real_8ayv_data

    rng1 = np.random.default_rng(42)
    msa1, del1, paired1 = construct_paired_msa(data, random=rng1, max_seqs=512)

    rng2 = np.random.default_rng(42)
    msa2, del2, paired2 = construct_paired_msa(data, random=rng2, max_seqs=512)

    np.testing.assert_array_equal(msa1, msa2, err_msg="Idempotency: MSA mismatch")
    np.testing.assert_array_equal(del1, del2, err_msg="Idempotency: deletion mismatch")
    np.testing.assert_array_equal(paired1, paired2, err_msg="Idempotency: paired mismatch")


def test_real_data_output_shape(real_8ayv_data):
    """Output arrays must have consistent shapes: (n_tokens, n_rows)."""
    data, _ = real_8ayv_data
    n_tokens = len(data.tokens)

    rng = np.random.default_rng(42)
    max_seqs = 128
    msa_data, del_data, paired_data = construct_paired_msa(data, random=rng, max_seqs=max_seqs)

    assert msa_data.shape[0] == n_tokens, f"msa_data dim 0 should be n_tokens={n_tokens}, got {msa_data.shape[0]}"
    assert del_data.shape == msa_data.shape, f"del_data shape {del_data.shape} != msa_data shape {msa_data.shape}"
    n_rows = msa_data.shape[1]
    assert n_rows >= 1, "Must have at least the query row"
    assert n_rows <= max_seqs, f"n_rows={n_rows} exceeds max_seqs={max_seqs}"


def test_real_data_query_row_matches_structure(real_8ayv_data):
    """Row 0 (query) should contain residue types matching the structure."""
    data, _ = real_8ayv_data

    rng = np.random.default_rng(42)
    msa_data, _, _ = construct_paired_msa(data, random=rng, max_seqs=128)

    for chain_id in sorted(data.msa.keys()):
        chain = data.structure.chains[chain_id]
        res_start = int(chain["res_idx"])
        res_num = int(chain["res_num"])
        expected = data.structure.residues[res_start : res_start + res_num]["res_type"]

        token_mask = data.tokens["asym_id"] == chain_id
        query_row = msa_data[token_mask, 0]

        np.testing.assert_array_equal(
            query_row,
            expected,
            err_msg=f"Chain {chain_id}: query row doesn't match structure residues",
        )
