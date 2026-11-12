"""Confidence indices remain traceable to original input chain IDs."""

import json
from dataclasses import asdict

import numpy as np
import pytest
import torch

from boltz.data import const
from boltz.data.types import (
    Atom,
    AtomV2,
    Bond,
    BondV2,
    Chain,
    ChainInfo,
    Connection,
    Coords,
    Ensemble,
    Interface,
    Record,
    Residue,
    Structure,
    StructureInfo,
    StructureV2,
)
from boltz.data.write.writer import BoltzWriter

SCALARS = [
    "confidence_score",
    "ptm",
    "iptm",
    "ligand_iptm",
    "protein_iptm",
    "complex_plddt",
    "complex_iplddt",
    "complex_pde",
    "complex_ipde",
]


def write_case(tmp_path, boltz2, output_format, mask, labels=None, confidence=True):
    count = len(mask)
    if labels is None:
        labels = [chr(65 + i) for i in range(count)]
    atoms = np.zeros(count, dtype=AtomV2 if boltz2 else Atom)
    atoms["is_present"] = True
    if boltz2:
        atoms["name"] = "CA"
    else:
        atoms["name"] = (ord("C") - 32, ord("A") - 32, 0, 0)
        atoms["element"] = 6
    residues = np.zeros(count, dtype=Residue)
    chains = np.zeros(count, dtype=Chain)
    for i in range(count):
        residues[i]["name"] = "GLY"
        residues[i]["res_type"] = const.token_ids["GLY"]
        residues[i]["atom_idx"] = residues[i]["atom_center"] = residues[i][
            "atom_disto"
        ] = i
        residues[i]["atom_num"] = 1
        residues[i]["is_standard"] = residues[i]["is_present"] = True
        chains[i]["name"] = chr(65 + i)
        chains[i]["mol_type"] = const.chain_type_ids["PROTEIN"]
        chains[i]["asym_id"] = i
        chains[i]["entity_id"] = 0  # Identical sequences are copies of one entity.
        chains[i]["res_idx"] = chains[i]["atom_idx"] = i
        chains[i]["res_num"] = chains[i]["atom_num"] = 1
    common = dict(
        atoms=atoms,
        residues=residues,
        chains=chains,
        interfaces=np.zeros(0, dtype=Interface),
        mask=np.array(mask),
    )
    if boltz2:
        structure = StructureV2(
            **common,
            bonds=np.zeros(0, dtype=BondV2),
            coords=np.zeros(count, dtype=Coords),
            ensemble=np.array([(0, count)], dtype=Ensemble),
        )
    else:
        structure = Structure(
            **common,
            bonds=np.zeros(0, dtype=Bond),
            connections=np.zeros(0, dtype=Connection),
        )
    source, output = tmp_path / "structures", tmp_path / "predictions"
    source.mkdir()
    np.savez_compressed(source / "toy.npz", **asdict(structure))
    record = Record(
        id="toy",
        structure=StructureInfo(num_chains=count),
        chains=[ChainInfo(i, label, 0, -1, -1, 1) for i, label in enumerate(labels)],
        interfaces=[],
    )
    valid = sum(mask)
    # Two samples with reversed ranks catch accidental label/ranking coupling.
    coords = torch.arange(2 * valid * 3, dtype=torch.float32).reshape(2, valid, 3)
    if not confidence:
        coords = coords[:1]
    prediction = {"exception": False, "coords": coords, "masks": torch.ones(1, valid)}
    if confidence:
        prediction.update({key: torch.tensor([0.25, 0.75]) for key in SCALARS})
        prediction["plddt"] = torch.full((2, valid), 0.9)
        prediction["pair_chains_iptm"] = {
            i: {
                j: torch.tensor([i / 10 + j / 100, i / 10 + j / 100 + 0.5])
                for j in range(valid)
            }
            for i in range(valid)
        }
    writer = BoltzWriter(source, output, output_format=output_format, boltz2=boltz2)
    writer.write_on_batch_end(None, None, prediction, [], {"record": [record]}, 0, 0)
    return (
        output / "toy",
        prediction,
        [label for label, kept in zip(labels, mask) if kept],
    )


@pytest.mark.parametrize("boltz2", [False, True])
@pytest.mark.parametrize("output_format", ["pdb", "mmcif"])
@pytest.mark.parametrize("mask", [[True, True, True], [True, False, True]])
def test_confidence_maps_reindexed_chains_and_preserves_all_scores(
    tmp_path, boltz2, output_format, mask
):
    output, prediction, names = write_case(tmp_path, boltz2, output_format, mask)
    for sample, rank in [(0, 1), (1, 0)]:
        summary = json.loads((output / f"confidence_toy_model_{rank}.json").read_text())
        assert summary.pop("chain_ids") == {
            str(i): name for i, name in enumerate(names)
        }
        pairs = {
            str(i): {str(j): value[sample].item() for j, value in row.items()}
            for i, row in prediction["pair_chains_iptm"].items()
        }
        expected = {key: prediction[key][sample].item() for key in SCALARS}
        expected.update(
            chains_ptm={i: row[i] for i, row in pairs.items()}, pair_chains_iptm=pairs
        )
        assert summary == expected
        suffix = "cif" if output_format == "mmcif" else "pdb"
        assert (output / f"toy_model_{rank}.{suffix}").is_file()
        plddt = np.load(output / f"plddt_toy_model_{rank}.npz")["plddt"]
        assert np.array_equal(plddt, prediction["plddt"][sample].numpy())


@pytest.mark.parametrize(
    "labels", [["input_receptor_long", "input_ligand_long"], ["same", "same"], ["only"]]
)
def test_original_names_are_not_truncated_or_used_as_unique_keys(tmp_path, labels):
    output, _, _ = write_case(
        tmp_path, True, "mmcif", [True] * len(labels), labels=labels
    )
    summary = json.loads((output / "confidence_toy_model_0.json").read_text())
    assert summary["chain_ids"] == {str(i): name for i, name in enumerate(labels)}


def test_prediction_without_confidence_still_writes_structures(tmp_path):
    output, _, _ = write_case(tmp_path, True, "mmcif", [True, True], confidence=False)
    assert (output / "toy_model_0.cif").is_file()
    assert not list(output.glob("confidence_*.json"))
