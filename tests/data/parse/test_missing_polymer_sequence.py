"""Coordinate templates without SEQRES must retain their observed residues."""

import warnings

import gemmi
import numpy as np
import pytest
from rdkit import Chem

from boltz.data import const
from boltz.data.parse.mmcif import parse_mmcif
from boltz.data.parse.pdb import parse_pdb


@pytest.fixture
def molecules():
    """Minimal named atom resources; no downloaded CCD assets are needed."""
    resources = {}
    for name in ["ALA", "GLY", "SER", "MET"]:
        mol = Chem.RWMol()
        for atom_name in const.ref_atoms[name]:
            atom = Chem.Atom(atom_name[0])
            atom.SetProp("name", atom_name)
            mol.AddAtom(atom)
        resources[name] = mol.GetMol()
    return resources


def structure(chains, full_sequences=None, models=1):
    """Build a named-coordinate fixture, including arbitrary author numbering."""
    result = gemmi.Structure()
    for model_idx in range(models):
        model = gemmi.Model(str(model_idx + 1))
        for chain_name, entries in chains.items():
            chain = gemmi.Chain(chain_name)
            for j, (res_name, seqnum, icode) in enumerate(entries):
                res = gemmi.Residue()
                res.name, res.seqid, res.het_flag = (
                    res_name,
                    gemmi.SeqId(seqnum, icode),
                    "A",
                )
                for i, atom_name in enumerate(const.ref_atoms[res_name]):
                    atom = gemmi.Atom()
                    atom.name = atom_name
                    atom.element = gemmi.Element(atom_name[0])
                    atom.pos = gemmi.Position(3 * j + i * 0.1, model_idx + 0.2, 0.3)
                    res.add_atom(atom)
                chain.add_residue(res)
            model.add_chain(chain)
        result.add_model(model)
    result.setup_entities()
    if full_sequences:
        for chain_name, sequence in full_sequences.items():
            entity = result.get_entity_of(result[0][chain_name].get_polymer())
            entity.full_sequence = sequence
    return result


@pytest.mark.parametrize("file_format", ["pdb", "cif"])
@pytest.mark.parametrize("models", [1, 2])
@pytest.mark.parametrize(
    "numbering", [((1, " "), (2, " ")), ((42, " "), (47, " ")), ((42, "A"), (42, "B"))]
)
def test_missing_sequence_uses_observed_residues(
    tmp_path, molecules, file_format, models, numbering
):
    entries = [(name, n, ins) for name, (n, ins) in zip(["ALA", "GLY"], numbering)]
    source = structure(
        {"A": entries, "B": [("SER", 8, " "), ("ALA", 9, " ")]}, models=models
    )
    path = tmp_path / ("template." + file_format)
    if file_format == "pdb":
        path.write_text(source.make_pdb_string())
        parser = parse_pdb
    else:
        source.make_mmcif_document().write_file(str(path))
        parser = parse_mmcif
    with pytest.warns(UserWarning, match="observed residues"):
        parsed = parser(
            str(path), mols=molecules, use_assembly=False, compute_interfaces=False
        )
    assert set(parsed.sequences.values()) == {"AG", "SA"}
    assert sorted(parsed.data.chains["res_num"].tolist()) == [2, 2]
    assert len(parsed.data.ensemble) == models
    assert parsed.data.atoms["is_present"].all()
    expected = np.array([a.pos.tolist() for c in source[0] for r in c for a in r])
    np.testing.assert_allclose(
        parsed.data.coords["coords"][: len(expected)], expected, atol=1e-6
    )


@pytest.mark.parametrize("file_format", ["pdb", "cif"])
def test_existing_full_sequence_retains_unresolved_residue(
    tmp_path, molecules, file_format
):
    source = structure(
        {"A": [("ALA", 1, " "), ("GLY", 3, " ")]},
        full_sequences={"A": ["ALA", "SER", "GLY"]},
    )
    path = tmp_path / ("template." + file_format)
    if file_format == "pdb":
        path.write_text(source.make_pdb_string())
        parser = parse_pdb
    else:
        source.make_mmcif_document().write_file(str(path))
        parser = parse_mmcif
    with warnings.catch_warnings(record=True) as emitted:
        parsed = parser(
            str(path), mols=molecules, use_assembly=False, compute_interfaces=False
        )
    assert not [w for w in emitted if "observed residues" in str(w.message)]
    assert list(parsed.sequences.values()) == ["ASG"]
    assert parsed.data.residues["name"].tolist() == ["ALA", "SER", "GLY"]
    assert parsed.data.residues["is_present"].tolist() == [True, False, True]


def test_only_missing_chain_sequence_is_inferred(tmp_path, molecules):
    source = structure(
        {
            "A": [("ALA", 1, " "), ("GLY", 3, " ")],
            "B": [("SER", 1, " "), ("ALA", 2, " ")],
        },
        full_sequences={"A": ["ALA", "SER", "GLY"]},
    )
    path = tmp_path / "mixed.cif"
    source.make_mmcif_document().write_file(str(path))
    with pytest.warns(UserWarning, match="observed residues"):
        parsed = parse_mmcif(
            str(path), mols=molecules, use_assembly=False, compute_interfaces=False
        )
    assert set(parsed.sequences.values()) == {"ASG", "SA"}
