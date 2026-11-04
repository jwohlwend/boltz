"""Explicit intra-chain bonds must reach guidance without duplicated forces."""

from types import SimpleNamespace

import numpy as np
import pytest
from rdkit.Chem import GetPeriodicTable
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


def test_repeated_and_reversed_internal_links_have_one_force():
    chains = np.zeros(1, dtype=Chain)
    rows = [
        (0, 0, 0, 1, 2, 8, COVALENT),
        (0, 0, 1, 0, 8, 2, COVALENT),
        (0, 0, 0, 1, 2, 8, COVALENT),
    ]
    data = SimpleNamespace(
        structure=SimpleNamespace(chains=chains, bonds=np.array(rows, dtype=BondV2))
    )
    before = data.structure.bonds.copy()
    features = process_chain_feature_constraints(data)
    assert features["connected_atom_index"].T.tolist() == [[2, 8]]
    np.testing.assert_array_equal(data.structure.bonds, before)


def test_self_bond_has_clear_error():
    data = SimpleNamespace(
        structure=SimpleNamespace(
            chains=np.zeros(1, dtype=Chain),
            bonds=np.array([(0, 0, 0, 0, 2, 2, COVALENT)], dtype=BondV2),
        )
    )
    with pytest.raises(ValueError, match="two distinct atoms"):
        process_chain_feature_constraints(data)


@pytest.mark.parametrize("elements", [(6, 16), (16, 16), (6, 6), (6, 7)])
def test_element_based_bounds_for_internal_links_only(elements):
    chains = np.zeros(2, dtype=Chain)
    data = SimpleNamespace(
        structure=SimpleNamespace(
            chains=chains,
            bonds=np.array(
                [
                    (0, 0, 0, 1, 0, 1, COVALENT),
                    (0, 1, 0, 2, 0, 2, COVALENT),
                ],
                dtype=BondV2,
            ),
        )
    )
    ref_elements = torch.nn.functional.one_hot(
        torch.tensor([*elements, 6]), const.num_elements
    )
    feats = process_chain_feature_constraints(data, ref_elements)
    radius = sum(GetPeriodicTable().GetRcovalent(e) for e in elements)
    torch.testing.assert_close(
        feats["connected_atom_lower_bounds"], torch.tensor([radius * 0.875, 0.0])
    )
    torch.testing.assert_close(
        feats["connected_atom_upper_bounds"], torch.tensor([radius * 1.125, 0.0])
    )


def test_unknown_element_has_explicit_fallback():
    data = SimpleNamespace(
        structure=SimpleNamespace(
            chains=np.zeros(1, dtype=Chain),
            bonds=np.array([(0, 0, 0, 1, 0, 1, COVALENT)], dtype=BondV2),
        )
    )
    elements = torch.nn.functional.one_hot(torch.tensor([0, 16]), const.num_elements)
    with pytest.warns(UserWarning, match="legacy connection upper bound"):
        feats = process_chain_feature_constraints(data, elements)
    assert feats["connected_atom_lower_bounds"].tolist() == [0.0]
    assert feats["connected_atom_upper_bounds"].tolist() == [0.0]


@pytest.mark.parametrize("distance", [0.0, 1e-10, 1.2, 1.8, 2.5])
def test_bounded_connection_energy_and_gradient(distance):
    features = {
        "connected_atom_index": torch.tensor([[[0], [1]]]),
        "connected_atom_lower_bounds": torch.tensor([[1.58375]]),
        "connected_atom_upper_bounds": torch.tensor([[2.03625]]),
    }
    xyz = torch.tensor([[[0.0, 0.0, 0.0], [distance, 0.0, 0.0]]], requires_grad=True)
    potential = ConnectionsPotential()
    energy = potential.compute(xyz, features, {"buffer": 2.0})
    actual = potential.compute_gradient(xyz.detach(), features, {"buffer": 2.0})
    expected = torch.autograd.grad(energy.sum(), xyz)[0]
    assert torch.isfinite(actual).all()
    torch.testing.assert_close(actual, expected)
    if 0 < distance < 1.58375:
        assert actual[0, 1, 0] < 0  # Gradient descent separates compressed endpoints.
    elif distance > 2.03625:
        assert actual[0, 1, 0] > 0
