"""Explicit intra-chain bonds must reach guidance and survive mmCIF output."""

from types import SimpleNamespace

import numpy as np
import pytest
import torch

from boltz.data import const
from boltz.data.feature.featurizerv2 import process_chain_feature_constraints
from boltz.data.types import (
    BondV2,
    Chain,
)
from boltz.model.potentials.potentials import ConnectionsPotential

COVALENT = const.bond_type_ids["COVALENT"]
SINGLE = const.bond_type_ids["SINGLE"]


@pytest.mark.parametrize(
    "rows,expected_atoms,expected_chains",
    [
        ([], [], []),
        ([(0, 0, 0, 0, 0, 1, SINGLE)], [], []),
        ([(0, 0, 0, 1, 2, 8, COVALENT)], [[2, 8]], []),
        ([(0, 0, 0, 0, 0, 2, COVALENT)], [[0, 2]], []),
        ([(0, 1, 0, 1, 2, 8, COVALENT)], [[2, 8]], [[0, 1]]),
        ([(0, 1, 0, 1, 2, 8, SINGLE)], [[2, 8]], [[0, 1]]),
        (
            [
                (0, 0, 0, 0, 0, 1, SINGLE),
                (0, 0, 0, 1, 2, 8, COVALENT),
                (0, 1, 1, 2, 8, 9, COVALENT),
            ],
            [[2, 8], [8, 9]],
            [[0, 1]],
        ),
    ],
)
def test_explicit_atom_and_chain_connectivity_are_separate(
    rows, expected_atoms, expected_chains
):
    chains = np.zeros(2, dtype=Chain)
    chains["entity_id"] = [0, 1]
    data = SimpleNamespace(
        structure=SimpleNamespace(chains=chains, bonds=np.array(rows, dtype=BondV2))
    )
    feats = process_chain_feature_constraints(data)
    assert feats["connected_atom_index"].T.tolist() == expected_atoms
    assert feats["connected_chain_index"].T.tolist() == expected_chains
    assert feats["connected_atom_index"].shape == (2, len(expected_atoms))
    assert feats["connected_chain_index"].shape == (2, len(expected_chains))


@pytest.mark.parametrize("distance", [1.8, 3.0, 8.0])
def test_same_chain_connection_energy_and_gradient(distance):
    chains = np.zeros(1, dtype=Chain)
    data = SimpleNamespace(
        structure=SimpleNamespace(
            chains=chains, bonds=np.array([(0, 0, 0, 1, 0, 1, COVALENT)], dtype=BondV2)
        )
    )
    feats = {k: v[None] for k, v in process_chain_feature_constraints(data).items()}
    xyz = torch.tensor([[[0.0, 0.0, 0.0], [distance, 0.0, 0.0]]], requires_grad=True)
    potential = ConnectionsPotential()
    parameters = {"buffer": 2.0}
    energy = potential.compute(xyz, feats, parameters)
    assert float(energy.item()) == pytest.approx(max(distance - 2.0, 0.0))
    expected = torch.autograd.grad(energy.sum(), xyz)[0]
    actual = potential.compute_gradient(xyz.detach(), feats, parameters)
    torch.testing.assert_close(actual, expected)
    if distance > 2.0:
        after = potential.compute(xyz.detach() - 0.1 * actual, feats, parameters)
        assert (after < energy).all()
