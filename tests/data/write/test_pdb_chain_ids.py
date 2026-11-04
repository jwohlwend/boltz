"""PDB chain labels must never shift its fixed-width coordinate columns."""

import json
import string
from dataclasses import replace

import gemmi
import numpy as np
import pytest
import torch

from boltz.data.types import (
    Atom,
    AtomV2,
    Bond,
    BondV2,
    Chain,
    Connection,
    Coords,
    Ensemble,
    Interface,
    Residue,
    Structure,
    StructureV2,
)
from boltz.data.write.pdb import to_pdb
from boltz.data.write.writer import BoltzWriter


def make_structure(names, boltz2):
    atoms = np.zeros(5 * len(names), dtype=AtomV2 if boltz2 else Atom)
    atoms["is_present"] = True
    atoms["coords"] = np.arange(len(atoms) * 3).reshape(-1, 3) * 0.01357 - 12.34567
    atom_names = ["N", "CA", "C", "O", "CB"]
    for i in range(len(atoms)):
        name = atom_names[i % 5]
        if boltz2:
            atoms["name"][i] = name
        else:
            atoms["name"][i] = [ord(c) - 32 for c in name] + [0] * (4 - len(name))
            atoms["element"][i] = {"N": 7, "C": 6, "O": 8}[name[0]]
    chains = np.zeros(len(names), dtype=Chain)
    residues = np.zeros(len(names), dtype=Residue)
    for i, name in enumerate(names):
        for field, value in {
            "name": name,
            "mol_type": 0,
            "entity_id": 0,
            "asym_id": i,
            "sym_id": i,
            "atom_idx": i * 5,
            "atom_num": 5,
            "res_idx": i,
            "res_num": 1,
        }.items():
            chains[field][i] = value
        for field, value in {
            "name": "ALA",
            "res_idx": 0,
            "atom_idx": i * 5,
            "atom_num": 5,
            "atom_center": i * 5 + 1,
            "atom_disto": i * 5 + 4,
            "is_standard": True,
            "is_present": True,
        }.items():
            residues[field][i] = value
    common = dict(
        atoms=atoms,
        residues=residues,
        chains=chains,
        interfaces=np.zeros(0, dtype=Interface),
        mask=np.ones(len(names), dtype=bool),
    )
    if boltz2:
        return StructureV2(
            **common,
            bonds=np.zeros(0, dtype=BondV2),
            coords=np.array([(tuple(x),) for x in atoms["coords"]], dtype=Coords),
            ensemble=np.array([(0, len(atoms))], dtype=Ensemble),
        )
    return Structure(
        **common,
        bonds=np.zeros(0, dtype=Bond),
        connections=np.zeros(0, dtype=Connection),
    )


@pytest.mark.parametrize("boltz2", [False, True])
@pytest.mark.parametrize(
    "names",
    [
        ["A", "B"],
        ["pro1", "lig1"],
        ["A1", "A2"],
        ["AA", "A", "B"],
        ["_", "a", "0"],
        ["", "X"],
    ],
)
def test_chain_roundtrip_and_fixed_columns(names, boltz2):
    structure = make_structure(names, boltz2)
    before = structure.chains.copy()
    text = to_pdb(structure, boltz2=boltz2)
    parsed = gemmi.read_pdb_string(text)
    assert len(parsed[0]) == len(names)
    ids = [c.name for c in parsed[0]]
    assert len(set(ids)) == len(names)
    assert all(len(c) == 1 and c in string.ascii_letters + string.digits for c in ids)
    for original in names:
        if len(original) == 1 and original in string.ascii_letters + string.digits:
            assert original in ids
    for line in text.splitlines():
        if line.startswith(("ATOM  ", "HETATM", "TER   ")):
            assert len(line) == 80
            assert line[21] in ids
    coords = np.array([a.pos.tolist() for c in parsed[0] for r in c for a in r])
    np.testing.assert_allclose(coords, structure.atoms["coords"], atol=0.000501, rtol=0)
    np.testing.assert_array_equal(structure.chains, before)
    maps = [
        line.strip().split(" ", 4)
        for line in text.splitlines()
        if line.startswith("REMARK 950 CHAIN_ID ")
    ]
    changed = {json.loads(parts[4]): parts[3] for parts in maps}
    assert set(changed) == {
        n for n in names if len(n) != 1 or n not in string.ascii_letters + string.digits
    }


@pytest.mark.parametrize("boltz2", [False, True])
def test_full_alphanumeric_capacity(boltz2):
    names = ["p" + str(i) for i in range(62)]
    text = to_pdb(make_structure(names, boltz2), boltz2=boltz2)
    assert len(gemmi.read_pdb_string(text)[0]) == 62


@pytest.mark.parametrize("boltz2", [False, True])
@pytest.mark.parametrize("names", [["same", "same"], ["p" + str(i) for i in range(63)]])
def test_unrepresentable_chain_set_fails_clearly(names, boltz2):
    with pytest.raises(ValueError, match="mmCIF"):
        to_pdb(make_structure(names, boltz2), boltz2=boltz2)


@pytest.mark.parametrize("boltz2", [False, True])
def test_cross_chain_conect_serials_and_ter_use_same_aliases(boltz2):
    structure = make_structure(["A1", "A2"], boltz2)
    if boltz2:
        structure = replace(
            structure, bonds=np.array([(0, 1, 0, 1, 0, 5, 1)], dtype=BondV2)
        )
    else:
        structure = replace(
            structure, connections=np.array([(0, 1, 0, 1, 0, 5)], dtype=Connection)
        )
    text = to_pdb(structure, boltz2=boltz2)
    atoms = [l for l in text.splitlines() if l.startswith("ATOM  ")]
    conect = [l for l in text.splitlines() if l.startswith("CONECT")]
    assert len(conect) == 1
    assert int(conect[0][6:11]) == int(atoms[0][6:11])
    assert int(conect[0][11:16]) == int(atoms[5][6:11])
    ter = next(l for l in text.splitlines() if l.startswith("TER"))
    assert ter[21] == atoms[0][21] and ter[21] != atoms[5][21]


@pytest.mark.parametrize("boltz2", [False, True])
def test_prediction_callback_exports_one_tensor_in_both_formats(tmp_path, boltz2):
    from boltz.data.types import ChainInfo, Record, StructureInfo

    structure = make_structure(["A1", "A2"], boltz2)
    data_dir = tmp_path / "processed"
    data_dir.mkdir()
    structure.dump(data_dir / "sample.npz")
    record = Record(
        id="sample",
        structure=StructureInfo(),
        interfaces=[],
        chains=[
            ChainInfo(
                chain_id=i,
                chain_name=name,
                mol_type=0,
                cluster_id=-1,
                msa_id=-1,
                num_residues=1,
            )
            for i, name in enumerate(["A1", "A2"])
        ],
    )
    prediction = {
        "exception": False,
        "coords": torch.tensor(structure.atoms["coords"].copy())[None],
        "masks": torch.ones(1, len(structure.atoms)),
    }
    exports = []
    for file_format, extension in [("pdb", "pdb"), ("mmcif", "cif")]:
        output = tmp_path / file_format
        writer = BoltzWriter(data_dir, output, output_format=file_format, boltz2=boltz2)
        writer.write_on_batch_end(
            None, None, prediction, [], {"record": [record]}, 0, 0
        )
        exports.append(
            gemmi.read_structure(
                str(output / "sample" / ("sample_model_0." + extension))
            )
        )
    xyz = [
        np.array([a.pos.tolist() for c in s[0] for r in c for a in r]) for s in exports
    ]
    assert len(exports[0][0]) == len(exports[1][0]) == 2
    np.testing.assert_allclose(xyz[0], xyz[1], rtol=0, atol=0.00051)
