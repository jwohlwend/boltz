"""Score a complete, noncovalent protein/CCD-ligand mmCIF without diffusion.

Experimental example, not a replacement for the prediction CLI. Chain arguments
are mmCIF label_asym_id (subchain) identifiers. No assembly expansion is applied.
Prepare missing atoms/residues externally; this example refuses to invent them.
"""

import argparse
import hashlib
import json
import time
from dataclasses import asdict, replace
from pathlib import Path

import numpy as np
import torch
from rdkit import Chem
from rdkit.Chem import Descriptors

from boltz.data import const
from boltz.data.crop.affinity import AffinityCropper
from boltz.data.feature.featurizerv2 import Boltz2Featurizer
from boltz.data.module.inferencev2 import collate
from boltz.data.mol import load_canonicals, load_molecules
from boltz.data.parse.mmcif import parse_mmcif
from boltz.data.tokenize.boltz2 import Boltz2Tokenizer
from boltz.data.types import AffinityInfo, Input, Record
from boltz.main import (
    Boltz2DiffusionParams,
    BoltzSteeringParams,
    MSAModuleArgs,
    PairformerArgsV2,
)
from boltz.model.models.boltz2 import Boltz2


def prepare_features(path, cache, protein_chains, ligand_chain):
    """Use the released parser, tokenizer, crop and chemical feature pipeline."""
    molecules = load_canonicals(cache / "mols")
    parsed = parse_mmcif(
        str(path),
        mols=molecules,
        moldir=str(cache / "mols"),
        use_assembly=False,
        compute_interfaces=False,
    )
    structure = parsed.data
    if len(structure.ensemble) != 1:
        raise ValueError("Supply one experimental conformation/model.")
    chain_map = {str(c["name"]): c for c in structure.chains}
    if not protein_chains or ligand_chain in protein_chains:
        raise ValueError("Select protein chain(s) and one distinct ligand subchain.")
    if set(protein_chains + [ligand_chain]) - chain_map.keys():
        raise ValueError(
            f"Unknown subchain; available label_asym_id values: {list(chain_map)}"
        )
    ligand = chain_map[ligand_chain]
    if (
        int(ligand["mol_type"]) != const.chain_type_ids["NONPOLYMER"]
        or int(ligand["res_num"]) != 1
    ):
        raise ValueError(
            "The affinity ligand must be a single non-polymer CCD residue."
        )
    if any(
        int(chain_map[c]["mol_type"]) != const.chain_type_ids["PROTEIN"]
        for c in protein_chains
    ):
        raise ValueError("All receptor subchains must be proteins.")
    selected = set(protein_chains + [ligand_chain])
    mask = np.array([str(c["name"]) in selected for c in structure.chains])
    structure = replace(structure, mask=mask)
    for chain in structure.chains[mask]:
        atoms = structure.atoms[
            chain["atom_idx"] : chain["atom_idx"] + chain["atom_num"]
        ]
        if not atoms["is_present"].all():
            raise ValueError(
                f"Subchain {chain['name']} has missing model atoms; prepare it externally first."
            )
    ligand_id = int(ligand["asym_id"])
    for bond in structure.bonds:
        if int(bond["chain_1"]) != int(bond["chain_2"]) and ligand_id in (
            int(bond["chain_1"]),
            int(bond["chain_2"]),
        ):
            raise ValueError(
                "This example does not support covalently attached ligands."
            )
    residue = structure.residues[int(ligand["res_idx"])]
    name = str(residue["name"])
    molecule = Chem.RemoveHs(molecules[name])
    if not 1 < molecule.GetNumAtoms() <= 128:
        raise ValueError("The affinity ligand must contain 2–128 model atoms.")
    mw = Descriptors.MolWt(molecule)
    record = Record(
        id=path.stem,
        structure=parsed.info,
        chains=[],
        interfaces=[],
        affinity=AffinityInfo(chain_id=ligand_id, mw=mw),
    )
    inputs = Input(structure=structure, msa={}, record=record)
    tokenized = Boltz2Tokenizer().tokenize(inputs)
    ligand_tokens = set(
        tokenized.tokens[tokenized.tokens["affinity_mask"]]["token_idx"]
    )
    tokenized = AffinityCropper().crop(tokenized, max_tokens=256, max_atoms=2048)
    if (
        set(tokenized.tokens[tokenized.tokens["affinity_mask"]]["token_idx"])
        != ligand_tokens
    ):
        raise ValueError("The crop did not retain the entire ligand.")
    names = set(tokenized.tokens["res_name"].tolist()) - molecules.keys()
    molecules.update(load_molecules(cache / "mols", names))
    feats = Boltz2Featurizer().process(
        tokenized,
        molecules=molecules,
        random=np.random.default_rng(42),
        training=False,
        max_atoms=None,
        max_tokens=None,
        max_seqs=1,
        pad_to_max_seqs=False,
        single_sequence_prop=1.0,
        compute_frames=True,
        compute_constraint_features=True,
        override_method="other",
        compute_affinity=True,
    )
    # Token and atom order come from the same tokenized object, not an index-only
    # assignment from a separately parsed ligand or PDB file.
    original = np.concatenate(
        [
            structure.coords["coords"][t["atom_idx"] : t["atom_idx"] + t["atom_num"]]
            for t in tokenized.tokens
        ]
    )
    observed = feats["coords"][0, : len(original)].numpy()
    translation = observed - original
    np.testing.assert_allclose(
        translation,
        np.broadcast_to(translation[0], translation.shape),
        atol=1e-5,
        rtol=0,
    )
    metadata = {
        "input_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "protein_subchains": protein_chains,
        "ligand_subchain": ligand_chain,
        "ligand_ccd": name,
        "ligand_mw": mw,
        "omitted_subchains": [
            str(c["name"]) for c in parsed.data.chains if str(c["name"]) not in selected
        ],
        "waters_omitted": True,
        "assembly_expanded": False,
        "msa": "single sequence",
        "tokens": len(tokenized.tokens),
        "atoms": len(original),
        "coordinate_change": "translation to crop centroid only (FP32 rounding)",
    }
    return collate([feats]), metadata


def load_model(cache, device, mw_correction=False):
    """Load the released affinity checkpoint with the prediction CLI defaults."""
    return (
        Boltz2.load_from_checkpoint(
            cache / "boltz2_aff.ckpt",
            map_location="cpu",
            strict=True,
            ema=False,
            use_kernels=False,
            affinity_mw_correction=mw_correction,
            diffusion_process_args=asdict(Boltz2DiffusionParams()),
            pairformer_args=asdict(PairformerArgsV2()),
            msa_args=asdict(MSAModuleArgs(use_paired_feature=True)),
            steering_args=asdict(BoltzSteeringParams(contact_guidance_update=False)),
        )
        .eval()
        .to(device)
    )


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("structure", type=Path)
    parser.add_argument("--protein-chain", action="append", required=True)
    parser.add_argument("--ligand-chain", required=True)
    parser.add_argument("--cache", type=Path, default=Path("~/.boltz").expanduser())
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--recycling-steps", type=int, default=5)
    parser.add_argument("--mw-correction", action="store_true")
    parser.add_argument("--prepare-only", action="store_true")
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Use a fresh output path to preserve previous results.")
    torch.manual_seed(0)
    torch.set_float32_matmul_precision("highest")
    start = time.monotonic()
    feats, metadata = prepare_features(
        args.structure, args.cache, args.protein_chain, args.ligand_chain
    )
    if not args.prepare_only:
        device = torch.device(args.device)
        model = load_model(args.cache, device, args.mw_correction)
        feats = {
            k: v.to(device) if isinstance(v, torch.Tensor) else v
            for k, v in feats.items()
        }
        before = feats["coords"].clone()
        with (
            torch.inference_mode(),
            torch.autocast(
                device.type, dtype=torch.bfloat16, enabled=device.type == "cuda"
            ),
        ):
            out = model(
                feats,
                recycling_steps=args.recycling_steps,
                affinity_from_input_structure=True,
            )
        assert torch.equal(feats["coords"], before)
        assert torch.equal(out["sample_atom_coords"], before[:, 0])
        metadata["predictions"] = {
            k: float(v.item()) for k, v in out.items() if k.startswith("affinity_")
        }
        metadata["input_coordinates_unchanged"] = True
        metadata["mw_correction"] = args.mw_correction
        metadata["recycling_steps"] = args.recycling_steps
        metadata["precision"] = (
            "bf16 autocast trunk / fp32 affinity heads"
            if device.type == "cuda"
            else "fp32"
        )
    metadata["elapsed_s"] = time.monotonic() - start
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(metadata, indent=2, allow_nan=False) + "\n")
    print(json.dumps(metadata, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
