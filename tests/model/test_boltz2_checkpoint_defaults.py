"""Checkpoint defaults must also work outside the predict CLI."""

from copy import deepcopy
from types import MappingProxyType

import pytest
import pytorch_lightning as pl
import torch
from omegaconf import OmegaConf
from torch import nn

from boltz.model.layers.attention import AttentionPairBias as AttentionV1
from boltz.model.layers.attentionv2 import AttentionPairBias as AttentionV2
from boltz.model.models.boltz2 import Boltz2
from boltz.model.modules import diffusionv2


def model_args():
    """Use small real modules without network or checkpoint downloads."""
    return dict(
        atom_s=16,
        atom_z=8,
        token_s=32,
        token_z=16,
        num_bins=8,
        training_args=OmegaConf.create(
            {
                "max_lr": 0.001,
                "weight_decay": 0.0,
                "diffusion_multiplicity": 1,
                "recycling_steps": 0,
            }
        ),
        validation_args={},
        embedder_args={"atom_encoder_depth": 1, "atom_encoder_heads": 4},
        msa_args={"msa_s": 16, "msa_blocks": 0, "msa_dropout": 0.0, "z_dropout": 0.0},
        pairformer_args={"num_blocks": 1, "num_heads": 4},
        score_model_args={
            "dim_fourier": 16,
            "atom_encoder_depth": 1,
            "atom_encoder_heads": 4,
            "token_transformer_depth": 1,
            "token_transformer_heads": 4,
            "atom_decoder_depth": 1,
            "atom_decoder_heads": 4,
            "conditioning_transition_layers": 1,
        },
        diffusion_process_args={},
        diffusion_loss_args={},
        validate_structure=False,
        confidence_prediction=False,
    )


@pytest.mark.parametrize("mapping", [dict, MappingProxyType, OmegaConf.create])
def test_missing_pairformer_version_uses_v2_without_mutating_config(mapping):
    args = model_args()
    pairformer_args = mapping(args["pairformer_args"])
    if OmegaConf.is_config(pairformer_args):
        OmegaConf.set_readonly(pairformer_args, True)
    before = dict(pairformer_args)
    args["pairformer_args"] = pairformer_args
    model = Boltz2(**args).eval()
    assert dict(pairformer_args) == before
    assert model.hparams.pairformer_args["v2"] is True
    assert isinstance(model.pairformer_module.layers[0].attention, AttentionV2)
    with torch.no_grad():
        s, z = model.pairformer_module(
            torch.randn(1, 3, 32),
            torch.randn(1, 3, 3, 16),
            mask=torch.ones(1, 3),
            pair_mask=torch.ones(1, 3, 3),
        )
    assert s.shape == (1, 3, 32)
    assert z.shape == (1, 3, 3, 16)
    assert torch.isfinite(s).all() and torch.isfinite(z).all()


@pytest.mark.parametrize("version", [False, True])
def test_explicit_pairformer_version_is_preserved(version):
    args = model_args()
    args["pairformer_args"]["v2"] = version
    before = deepcopy(args)
    model = Boltz2(**args)
    assert args == before
    assert model.hparams.pairformer_args["v2"] is version
    assert isinstance(
        model.pairformer_module.layers[0].attention,
        AttentionV2 if version else AttentionV1,
    )


def test_checkpoint_round_trip_records_normalized_pairformer_defaults(tmp_path):
    model = Boltz2(**model_args())
    checkpoint = tmp_path / "tiny.ckpt"
    torch.save(
        {
            "state_dict": model.state_dict(),
            "hyper_parameters": dict(model.hparams),
            "pytorch-lightning_version": pl.__version__,
        },
        checkpoint,
    )
    restored = Boltz2.load_from_checkpoint(checkpoint, weights_only=False)
    assert restored.hparams.pairformer_args["v2"] is True
    assert model.state_dict().keys() == restored.state_dict().keys()
    for key, value in model.state_dict().items():
        assert torch.equal(value, restored.state_dict()[key]), key


class ZeroScore(nn.Module):
    def __init__(self, **kwargs):
        super().__init__()
        self.anchor = nn.Parameter(torch.zeros(()))

    def forward(self, r_noisy, **kwargs):
        return torch.zeros_like(r_noisy)


@pytest.fixture
def sampler(monkeypatch):
    monkeypatch.setattr(diffusionv2, "DiffusionModule", ZeroScore)
    return diffusionv2.AtomDiffusion(
        score_model_args={"token_s": 8},
        alignment_reverse_diff=False,
        step_scale_random=None,
    ).eval()


@pytest.mark.parametrize(
    "steering",
    [
        None,
        {},
        {"fk_steering": False},
        {"physical_guidance_update": False},
        {"contact_guidance_update": False},
    ],
)
def test_missing_steering_flags_match_explicit_no_guidance(sampler, steering):
    before = deepcopy(steering)
    explicit = dict(
        fk_steering=False, physical_guidance_update=False, contact_guidance_update=False
    )
    torch.manual_seed(17)
    reference = sampler.sample(
        torch.ones(1, 4), num_sampling_steps=3, steering_args=explicit
    )
    torch.manual_seed(17)
    actual = sampler.sample(
        torch.ones(1, 4), num_sampling_steps=3, steering_args=steering
    )
    assert steering == before
    assert torch.equal(reference["sample_atom_coords"], actual["sample_atom_coords"])
    assert torch.isfinite(actual["sample_atom_coords"]).all()


@pytest.mark.parametrize(
    "flag", ["fk_steering", "physical_guidance_update", "contact_guidance_update"]
)
def test_explicit_guidance_still_runs(sampler, monkeypatch, flag):
    calls = []

    class Potential:
        def compute_parameters(self, t):
            return dict(resampling_weight=1, guidance_weight=1, guidance_interval=1)

        def compute(self, coords, feats, params):
            calls.append("energy")
            return coords.square().sum(dim=(-1, -2))

        def compute_gradient(self, coords, feats, params):
            calls.append("gradient")
            return coords * 0.01

    def potentials(args, boltz2):
        assert args[flag] is True and boltz2
        return [Potential()]

    monkeypatch.setattr(diffusionv2, "get_potentials", potentials)
    steering = {
        flag: True,
        "num_particles": 2,
        "fk_lambda": 0.1,
        "fk_resampling_interval": 1,
        "num_gd_steps": 1,
    }
    before = deepcopy(steering)
    output = sampler.sample(
        torch.ones(1, 4), num_sampling_steps=3, steering_args=steering, feats={}
    )
    assert steering == before
    assert ("energy" if flag == "fk_steering" else "gradient") in calls
    assert output["sample_atom_coords"].shape == (1, 4, 3)
    assert torch.isfinite(output["sample_atom_coords"]).all()
