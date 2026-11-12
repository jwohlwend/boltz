"""Replay a saved Boltz-2 prediction through the original and patched writer.

No model inference or GPU is needed. Supply a completed affinity-enabled run
containing pre_affinity_<record>.npz and its confidence/PLDDT output files.
"""

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import torch

from boltz.data.types import Record, StructureV2
from boltz.data.write.writer import BoltzWriter


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--record-id", required=True)
    parser.add_argument("--original-writer", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    record_id = args.record_id
    source = args.run / "processed/structures"
    record = Record.load(args.run / "processed/records" / f"{record_id}.json")
    previous = args.run / "predictions" / record_id
    confidence_path = previous / f"confidence_{record_id}_model_0.json"
    summary = json.loads(confidence_path.read_text())
    structure = StructureV2.load(previous / f"pre_affinity_{record_id}.npz")
    coords = torch.from_numpy(structure.coords["coords"].copy()).unsqueeze(0)
    prediction = {
        "exception": False,
        "coords": coords,
        "masks": torch.ones(1, coords.shape[1]),
        "plddt": torch.from_numpy(
            np.load(previous / f"plddt_{record_id}_model_0.npz")["plddt"]
        ).unsqueeze(0),
    }
    for key in [
        "confidence_score",
        "ptm",
        "iptm",
        "ligand_iptm",
        "protein_iptm",
        "complex_plddt",
        "complex_iplddt",
        "complex_pde",
        "complex_ipde",
    ]:
        prediction[key] = torch.tensor([summary[key]], dtype=torch.float64)
    prediction["pair_chains_iptm"] = {
        int(i): {int(j): torch.tensor([v], dtype=torch.float64) for j, v in row.items()}
        for i, row in summary["pair_chains_iptm"].items()
    }
    spec = importlib.util.spec_from_file_location(
        "original_boltz_writer", args.original_writer
    )
    original = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(original)
    for name, writer_type in [
        ("original", original.BoltzWriter),
        ("patched", BoltzWriter),
    ]:
        writer = writer_type(
            source, args.output / name, boltz2=True, output_format="mmcif"
        )
        writer.write_on_batch_end(
            None, None, prediction, [], {"record": [record]}, 0, 0
        )
    relative = Path(record_id) / f"confidence_{record_id}_model_0.json"
    before = json.loads((args.output / "original" / relative).read_text())
    after = json.loads((args.output / "patched" / relative).read_text())
    chain_ids = after.pop("chain_ids")
    assert after == before == summary
    source_structure = StructureV2.load(source / f"{record_id}.npz")
    names = [
        chain.chain_name
        for chain, kept in zip(record.chains, source_structure.mask)
        if kept
    ]
    assert chain_ids == {str(i): name for i, name in enumerate(names)}
    checked = []
    for path in (args.output / "original" / record_id).iterdir():
        other = args.output / "patched" / record_id / path.name
        if path.suffix == ".json":
            continue
        if path.suffix == ".npz":
            # Compare .npy payloads without unpickling optional object arrays;
            # ZIP container timestamps need not match between writer calls.
            with ZipFile(path) as left, ZipFile(other) as right:
                assert left.namelist() == right.namelist()
                assert all(left.read(key) == right.read(key) for key in left.namelist())
        else:
            assert path.read_bytes() == other.read_bytes()
        checked.append(path.name)
    report = {
        "record_id": record_id,
        "source_confidence_sha256": hashlib.sha256(
            confidence_path.read_bytes()
        ).hexdigest(),
        "chain_ids": chain_ids,
        "all_legacy_confidence_fields_equal_saved_prediction": True,
        "unchanged_output_files_or_npz_payloads": sorted(checked),
        "execution": "CPU replay of saved prediction; no new model inference or quality claim",
    }
    (args.output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
