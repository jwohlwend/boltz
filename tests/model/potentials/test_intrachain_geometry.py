"""Geometry guidance must balance closure without affecting unrelated chains."""

from types import SimpleNamespace

import numpy as np
import pytest
import torch

from boltz.data import const
from boltz.data.feature.featurizerv2 import process_residue_constraint_features
from boltz.data.types import (
    BondV2,
    Chain,
    ChiralAtomConstraint,
    PlanarBondConstraint,
    PlanarRing5Constraint,
    PlanarRing6Constraint,
    RDKitBoundsConstraint,
    ResidueConstraints,
    StereoBondConstraint,
)
from boltz.model.potentials.potentials import (
    ConnectionsPotential,
    IntraChainGeometryPotential,
    PoseBustersPotential,
    get_potentials,
)

COVALENT = const.bond_type_ids["COVALENT"]
SINGLE = const.bond_type_ids["SINGLE"]
PARAMETERS = {"bond_buffer": 0.125, "angle_buffer": 0.125, "clash_buffer": 0.1}


def constraint_input(bonds):
    chains = np.zeros(2, dtype=Chain)
    # Nonconsecutive IDs ensure selection uses the chain identity, not its row.
    chains["asym_id"] = [2, 5]
    chains["atom_idx"] = [0, 4]
    chains["atom_num"] = [4, 4]
    constraints = ResidueConstraints(
        rdkit_bounds_constraints=np.array(
            [
                ((0, 1), True, False, 1.52, 1.48),
                ((1, 2), True, False, 1.52, 1.48),
                ((4, 5), True, False, 1.52, 1.48),
                ((0, 2), False, True, 2.6, 2.4),
            ],
            dtype=RDKitBoundsConstraint,
        ),
        chiral_atom_constraints=np.empty(0, dtype=ChiralAtomConstraint),
        stereo_bond_constraints=np.empty(0, dtype=StereoBondConstraint),
        planar_bond_constraints=np.empty(0, dtype=PlanarBondConstraint),
        planar_ring_5_constraints=np.empty(0, dtype=PlanarRing5Constraint),
        planar_ring_6_constraints=np.empty(0, dtype=PlanarRing6Constraint),
    )
    return SimpleNamespace(
        residue_constraints=constraints,
        structure=SimpleNamespace(chains=chains, bonds=np.array(bonds, dtype=BondV2)),
    )


@pytest.mark.parametrize(
    "bonds,expected",
    [
        ([], []),
        ([(2, 2, 0, 0, 0, 1, SINGLE)], []),
        ([(2, 5, 0, 2, 0, 4, COVALENT)], []),
        ([(2, 2, 0, 1, 0, 3, COVALENT)], [0, 1, 3]),
        ([(5, 5, 2, 3, 4, 7, COVALENT)], [2]),
        ([(2, 2, 0, 1, 0, 3, COVALENT), (5, 5, 2, 3, 4, 7, COVALENT)], [0, 1, 2, 3]),
    ],
)
def test_select_existing_bounds_only_in_explicitly_linked_chains(bonds, expected):
    data = constraint_input(bonds)
    before = data.residue_constraints.rdkit_bounds_constraints.copy()
    features = process_residue_constraint_features(data)
    assert features["intrachain_covalent_bound_indices"].tolist() == expected
    np.testing.assert_array_equal(
        data.residue_constraints.rdkit_bounds_constraints, before
    )


def test_missing_residue_constraints_are_empty():
    data = constraint_input([(2, 2, 0, 1, 0, 3, COVALENT)])
    data.residue_constraints = None
    features = process_residue_constraint_features(data)
    assert features["intrachain_covalent_bound_indices"].shape == (0,)
    assert features["rdkit_bounds_index"].shape == (2, 0)


def potential_features():
    data = constraint_input([(2, 2, 0, 1, 0, 3, COVALENT)])
    features = {
        k: v[None] for k, v in process_residue_constraint_features(data).items()
    }
    features["ref_element"] = torch.nn.functional.one_hot(
        torch.full((1, 8), 6), const.num_elements
    )
    return features


DEVICES = ["cpu"] + (["cuda"] if torch.cuda.is_available() else [])


@pytest.mark.parametrize("device", DEVICES)
@pytest.mark.parametrize("coincident", [False, True])
def test_selected_energy_and_analytic_gradient_match_autograd(device, coincident):
    features = {k: v.to(device) for k, v in potential_features().items()}
    xyz = torch.tensor(
        [
            [
                [0.0, 0.0, 0.0],
                [2.2, 0.1, 0.0],
                [3.8, 0.7, 0.0],
                [5.0, 1.0, 0.0],
                [0.0, 5.0, 0.0],
                [2.5, 5.0, 0.0],
                [4.0, 5.0, 0.0],
                [6.0, 5.0, 0.0],
            ]
        ],
        device=device,
    )
    if coincident:
        xyz[:, 1] = xyz[:, 0]
    xyz.requires_grad_()
    before = {k: v.clone() for k, v in features.items()}
    potential = IntraChainGeometryPotential()
    energy = potential.compute(xyz, features, PARAMETERS)
    expected = torch.autograd.grad(energy.sum(), xyz)[0]
    actual = potential.compute_gradient(xyz.detach(), features, PARAMETERS)
    assert torch.isfinite(energy).all() and torch.isfinite(actual).all()
    torch.testing.assert_close(actual, expected)
    assert torch.count_nonzero(actual[:, 4:]) == 0
    for key in features:
        torch.testing.assert_close(features[key], before[key], rtol=0, atol=0)
    if not coincident:
        assert (
            potential.compute(xyz.detach() - 0.01 * actual, features, PARAMETERS)
            < energy
        ).all()


@pytest.mark.parametrize("legacy", [False, True])
def test_no_selected_bounds_is_exactly_inactive(legacy):
    features = potential_features()
    if legacy:
        del features["intrachain_covalent_bound_indices"]
    else:
        features["intrachain_covalent_bound_indices"] = torch.empty(
            (1, 0), dtype=torch.long
        )
    xyz = torch.randn(3, 8, 3)
    potential = IntraChainGeometryPotential()
    assert torch.equal(potential.compute(xyz, features, PARAMETERS), torch.zeros(3))
    assert torch.equal(
        potential.compute_gradient(xyz, features, PARAMETERS), torch.zeros_like(xyz)
    )


@pytest.mark.parametrize(
    "boltz2,physical", [(False, True), (True, False), (True, True)]
)
def test_extra_geometry_is_gated_and_balances_connection_force(boltz2, physical):
    flags = {
        "physical_guidance_update": physical,
        "fk_steering": True,
        "contact_guidance_update": True,
    }
    potentials = get_potentials(flags, boltz2=boltz2)
    extras = [p for p in potentials if isinstance(p, IntraChainGeometryPotential)]
    assert len(extras) == int(boltz2 and physical)
    if extras:
        closure = next(p for p in potentials if type(p) is ConnectionsPotential)
        normal = next(p for p in potentials if type(p) is PoseBustersPotential)
        assert extras[0].parameters["resampling_weight"] == 0
        assert normal.parameters["guidance_weight"] + extras[0].parameters[
            "guidance_weight"
        ] == pytest.approx(closure.parameters["guidance_weight"])
