"""Compare direct checkpoint inference with explicit predict-CLI defaults.

Run after preprocessing a small input with ``boltz predict``. This reuses only
its processed features, not cached predictions. No MSA server is contacted.
"""

import argparse
import hashlib
import json
from dataclasses import asdict
from pathlib import Path

import torch

from boltz.data.module.inferencev2 import PredictionDataset, collate
from boltz.data.types import Manifest
from boltz.main import Boltz2DiffusionParams
from boltz.model.layers.pairformer import PairformerModule
from boltz.model.models.boltz2 import Boltz2


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--processed", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    torch.set_num_threads(4)
    torch.set_float32_matmul_precision("highest")
    torch.cuda.reset_peak_memory_stats()
    processed = args.processed
    dataset = PredictionDataset(
        manifest=Manifest.load(processed / "manifest.json"),
        target_dir=processed / "structures",
        msa_dir=processed / "msa",
        mol_dir=args.cache / "mols",
        constraints_dir=processed / "constraints",
        template_dir=processed / "templates",
        extra_mols_dir=processed / "mols",
    )
    feats = {
        k: v.cuda() if isinstance(v, torch.Tensor) else v
        for k, v in collate([dataset[0]]).items()
    }
    checkpoint = args.cache / "boltz2_conf.ckpt"
    saved = torch.load(checkpoint, map_location="cpu", mmap=True, weights_only=False)
    hparams = saved["hyper_parameters"]
    pairformer_before = dict(hparams["pairformer_args"])
    assert "v2" not in pairformer_before and "steering_args" not in hparams

    # Reproduce #679 using the released configuration and actual attention layer.
    legacy_args = dict(pairformer_before)
    legacy_args.update(num_blocks=1, activation_checkpointing=False)
    legacy = (
        PairformerModule(hparams["token_s"], hparams["token_z"], **legacy_args)
        .cuda()
        .eval()
    )
    try:
        legacy(
            torch.zeros(1, 2, hparams["token_s"], device="cuda"),
            torch.zeros(1, 2, 2, hparams["token_z"], device="cuda"),
            mask=torch.ones(1, 2, device="cuda"),
            pair_mask=torch.ones(1, 2, 2, device="cuda"),
        )
    except TypeError as error:
        assert "k_in" in str(error)
        original_error = str(error)
    else:
        raise AssertionError("Expected the upstream default attention to fail")
    del legacy

    # Use predict's inference sampler settings. Historical checkpoint training
    # metadata contains sampler options removed from the current API.
    # Neither pairformer_args nor steering_args is overridden here.
    model = Boltz2.load_from_checkpoint(
        checkpoint,
        map_location="cpu",
        strict=True,
        weights_only=False,
        diffusion_process_args=asdict(Boltz2DiffusionParams()),
    ).eval()
    assert dict(hparams["pairformer_args"]) == pairformer_before
    assert model.steering_args is None
    model.cuda()
    options = dict(recycling_steps=1, num_sampling_steps=20, diffusion_samples=1)
    with torch.inference_mode(), torch.autocast("cuda", dtype=torch.bfloat16):
        torch.manual_seed(91)
        direct = model(feats, **options)
        model.steering_args = {
            "fk_steering": False,
            "physical_guidance_update": False,
            "contact_guidance_update": False,
        }
        torch.manual_seed(91)
        explicit = model(feats, **options)
    checked = []
    for key, value in direct.items():
        if isinstance(value, torch.Tensor):
            assert torch.isfinite(value).all(), key
            assert torch.equal(value, explicit[key]), key
            checked.append(key)
    assert "sample_atom_coords" in checked and "iptm" in checked
    report = {
        "checkpoint": checkpoint.name,
        "checkpoint_sha256": hashlib.sha256(checkpoint.read_bytes()).hexdigest(),
        "original_pairformer_error": original_error,
        "strict_checkpoint_load": True,
        "saved_config_unmodified": True,
        "sampler_settings": asdict(Boltz2DiffusionParams()),
        "none_vs_explicit_disabled_guidance_bit_equal": checked,
        "settings": options,
        "input_id": dataset.manifest.records[0].id,
        "tokens": int(feats["token_pad_mask"].sum()),
        "gpu": torch.cuda.get_device_name(),
        "torch": torch.__version__,
        "peak_allocated_mib": torch.cuda.max_memory_allocated() / 2**20,
        "peak_reserved_mib": torch.cuda.max_memory_reserved() / 2**20,
        "limitations": "One confidence-checkpoint input, shortened diffusion schedule; "
        "correctness smoke test, not a quality or speed benchmark.",
    }
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
