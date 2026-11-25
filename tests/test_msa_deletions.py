"""Regression test for MSA deletion features (jwohlwend/boltz#730).

Commit 04d27c7 re-sliced ``chain_deletions`` from its own previous slice inside the
per-sequence loop. The query row has no deletions, so after the first iteration the
array was empty and every deletion feature came out as zero.
"""

import io
from types import ModuleType, SimpleNamespace

import numpy as np
import pytest

from boltz.data import const
from boltz.data.feature import featurizer, featurizerv2
from boltz.data.parse.a3m import _parse_a3m
from boltz.data.types import MSA

# Lowercase letters are insertions: hom1 has 2 before residue index 3,
# hom2 has 1 before residue index 5.
A3M = """>query
ACDEFGHIKL
>hom1
ACDkkEFGHIKL
>hom2
ACDEFmGHIKL
"""
QUERY = "ACDEFGHIKL"


def _tokenized(msa: MSA) -> SimpleNamespace:
    n = len(QUERY)
    residues = np.zeros(n, dtype=[("res_type", "i1")])
    residues["res_type"] = [
        const.token_ids[const.prot_letter_to_token[c]] for c in QUERY
    ]
    chains = np.zeros(1, dtype=[("res_idx", "i4"), ("res_num", "i4")])
    chains["res_num"] = n
    tokens = np.zeros(n, dtype=[("asym_id", "i4"), ("res_idx", "i4")])
    tokens["res_idx"] = np.arange(n)
    return SimpleNamespace(
        tokens=tokens,
        structure=SimpleNamespace(chains=chains, residues=residues),
        msa={0: msa},
        record=SimpleNamespace(id="test"),
    )


@pytest.mark.parametrize("module", [featurizer, featurizerv2])
def test_deletion_features_are_populated(module: ModuleType) -> None:
    """Each homolog's deletion count lands on its residue; nothing else is set."""
    msa = _parse_a3m(io.StringIO(A3M), taxonomy=None)
    if module is featurizerv2:
        args = (_tokenized(msa), np.random.default_rng(0))
    else:  # the Boltz-1 featurizer takes no random generator
        args = (_tokenized(msa),)
    _, deletion, _ = module.construct_paired_msa(*args, max_seqs=16)
    deletion = np.asarray(deletion)  # (tokens, sequences)

    expected = np.zeros((len(QUERY), 3))
    expected[3, 1] = 2  # hom1
    expected[5, 2] = 1  # hom2
    np.testing.assert_array_equal(deletion, expected)
