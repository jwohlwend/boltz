"""Explicit intra-chain bonds must reach guidance and survive mmCIF output."""

from dataclasses import replace

import gemmi
import numpy as np
import pytest
import torch

from boltz.data import const
from boltz.data.types import (
    AtomV2,
    BondV2,
    Chain,
    Coords,
    Ensemble,
    Interface,
    Residue,
    StructureV2,
)
from boltz.data.write.mmcif import _add_covalent_connections, to_mmcif

COVALENT = const.bond_type_ids["COVALENT"]
SINGLE = const.bond_type_ids["SINGLE"]


def linked_structure(interchain=False, duplicate=False):
    names = ["C", "O", "CH3", "N", "CA", "C", "O", "CB", "SG", "N1"]
    atoms = np.zeros(len(names), dtype=AtomV2)
    atoms["name"] = names
    atoms["is_present"] = True
    atoms["coords"] = np.arange(len(names) * 3).reshape(-1, 3) * 0.1
    atoms["coords"][2] = [0.0, 0.0, 0.0]
    atoms["coords"][8] = [1.81, 0.0, 0.0]
    residues = np.zeros(2, dtype=Residue)
    residues["name"] = ["ACE", "CY3"]
    residues["atom_idx"] = [0, 3]
    residues["atom_num"] = [3, 7]
    residues["res_idx"] = [0, 0 if interchain else 1]
    residues["is_present"] = True
    chains = np.zeros(2 if interchain else 1, dtype=Chain)
    chains["name"] = ["longA", "longB"] if interchain else ["B"]
    chains["asym_id"] = np.arange(len(chains))
    chains["entity_id"] = np.arange(len(chains))
    chains["atom_idx"] = [0, 3] if interchain else [0]
    chains["atom_num"] = [3, 7] if interchain else [10]
    chains["res_idx"] = [0, 1] if interchain else [0]
    chains["res_num"] = [1, 1] if interchain else [2]
    chains["mol_type"] = (
        const.chain_type_ids["NONPOLYMER"]
        if interchain
        else const.chain_type_ids["PROTEIN"]
    )
    bonds = [(0, 1 if interchain else 0, 0, 1, 2, 8, COVALENT)]
    if duplicate:
        bonds.append((1 if interchain else 0, 0, 1, 0, 8, 2, COVALENT))
    return StructureV2(
        atoms=atoms,
        bonds=np.array(bonds, dtype=BondV2),
        residues=residues,
        chains=chains,
        interfaces=np.zeros(0, dtype=Interface),
        mask=np.ones(len(chains), dtype=bool),
        coords=np.array([(tuple(x),) for x in atoms["coords"]], dtype=Coords),
        ensemble=np.array([(0, len(atoms))], dtype=Ensemble),
    )


@pytest.mark.parametrize("interchain", [False, True])
@pytest.mark.parametrize("duplicate", [False, True])
def test_mmcif_has_explicit_atoms_and_measured_connection(interchain, duplicate):
    structure = linked_structure(interchain, duplicate)
    plddt = None if interchain else torch.tensor([0.8, 0.7])
    text = to_mmcif(structure, boltz2=True, plddts=plddt)
    block = gemmi.cif.read_string(text).sole_block()
    parsed = gemmi.make_structure_from_block(block)
    assert len(parsed.connections) == 1
    connection = parsed.connections[0]
    assert connection.type == gemmi.ConnectionType.Covale
    assert {connection.partner1.atom_name, connection.partner2.atom_name} == {
        "CH3",
        "SG",
    }
    assert connection.reported_distance == pytest.approx(1.81, abs=1e-5)
    if not interchain:
        assert connection.partner1.res_id.seqid.num == 1
        assert connection.partner2.res_id.seqid.num == 2
        assert len(block.find_values("_ma_qa_metric_local.metric_value")) == 2
    xyz = np.array([a.pos.tolist() for c in parsed[0] for r in c for a in r])
    np.testing.assert_allclose(xyz, structure.atoms["coords"], atol=1e-5, rtol=0)


def test_missing_endpoint_does_not_emit_a_dangling_connection():
    structure = linked_structure()
    atoms = structure.atoms.copy()
    atoms["is_present"][8] = False
    text = to_mmcif(replace(structure, atoms=atoms), boltz2=True)
    assert not gemmi.make_structure_from_block(
        gemmi.cif.read_string(text).sole_block()
    ).connections


def test_non_covalent_residue_internal_bonds_do_not_become_struct_conn():
    structure = linked_structure()
    bonds = structure.bonds.copy()
    bonds["type"] = SINGLE
    text = to_mmcif(replace(structure, bonds=bonds), boltz2=True)
    assert not gemmi.make_structure_from_block(
        gemmi.cif.read_string(text).sole_block()
    ).connections


def test_appending_connections_is_idempotent():
    structure = linked_structure()
    text = to_mmcif(structure, boltz2=True)
    assert _add_covalent_connections(text, structure, True) == text


def test_existing_non_covalent_row_and_quality_categories_are_preserved():
    structure = linked_structure()
    document = gemmi.cif.read_string(
        to_mmcif(structure, boltz2=True, plddts=torch.tensor([0.8, 0.7]))
    )
    block = document.sole_block()
    old = block.get_mmcif_category("_struct_conn.", raw=True)
    old["conn_type_id"][0] = "hydrog"
    block.set_mmcif_category("_struct_conn.", old, raw=True)
    original_quality = block.get_mmcif_category("_ma_qa_metric_local.", raw=True)
    result = _add_covalent_connections(document.as_string(), structure, True)
    updated = gemmi.cif.read_string(result).sole_block()
    category = updated.get_mmcif_category("_struct_conn.", raw=True)
    assert len(category["id"]) == 2
    for column, values in old.items():
        assert category[column][0] == values[0]
    assert category["conn_type_id"] == ["hydrog", "covale"]
    assert (
        updated.get_mmcif_category("_ma_qa_metric_local.", raw=True) == original_quality
    )


def test_legacy_boltz_atom_names_and_connections_are_supported():
    from boltz.data.types import Atom, Bond, Connection, Structure

    source = linked_structure()
    atoms = np.zeros(len(source.atoms), dtype=Atom)
    atoms["coords"] = source.atoms["coords"]
    atoms["is_present"] = True
    for i, name in enumerate(source.atoms["name"]):
        atoms["name"][i] = [ord(c) - 32 for c in name] + [0] * (4 - len(name))
        atoms["element"][i] = {"C": 6, "N": 7, "O": 8, "S": 16}[name[0]]
    structure = Structure(
        atoms=atoms,
        bonds=np.zeros(0, dtype=Bond),
        residues=source.residues,
        chains=source.chains,
        interfaces=source.interfaces,
        mask=source.mask,
        connections=np.array([(0, 0, 0, 1, 2, 8)], dtype=Connection),
    )
    text = to_mmcif(structure)
    parsed = gemmi.make_structure_from_block(gemmi.cif.read_string(text).sole_block())
    assert len(parsed.connections) == 1
    assert {
        parsed.connections[0].partner1.atom_name,
        parsed.connections[0].partner2.atom_name,
    } == {"CH3", "SG"}
