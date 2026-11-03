"""Exercise the real Boltz2 forward routing without loading model weights."""

from types import SimpleNamespace
from unittest.mock import Mock

import pytest
import torch
from pytorch_lightning import LightningModule
from torch import nn

from boltz.model.models.boltz2 import Boltz2


class Embedder(nn.Module):
    def forward(self, feats, affinity=False):
        return feats["embedding"]


class Head(nn.Module):
    def __init__(self, offset):
        super().__init__()
        self.offset = offset
        self.coordinates = None

    def forward(self, *, x_pred, **kwargs):
        self.coordinates = x_pred
        value = x_pred.square().sum().reshape(1) + self.offset
        return {"affinity_pred_value": value, "affinity_logits_binary": value / 100}


def make_model(ensemble=True, correction=False):
    model = Boltz2.__new__(Boltz2)
    LightningModule.__init__(model)
    model.structure_prediction_training = False
    model.affinity_prediction = True
    model.affinity_ensemble = ensemble
    model.affinity_mw_correction = correction
    model.run_trunk_and_structure = True
    model.skip_run_structure = False
    model.confidence_prediction = True
    model.bond_type_feature = False
    model.use_templates = False
    model.use_kernels = False
    model.is_msa_compiled = False
    model.is_pairformer_compiled = False
    model.predict_bfactor = False
    model.checkpoint_diffusion_conditioning = False
    model.steering_args = {}
    model.input_embedder = Embedder()
    model.s_init = nn.Identity()
    model.z_init_1 = nn.Linear(3, 2, bias=False)
    model.z_init_2 = nn.Linear(3, 2, bias=False)
    model.rel_pos = Mock(return_value=torch.zeros(1, 3, 3, 2))
    model.token_bonds = nn.Linear(1, 2, bias=False)
    model.contact_conditioning = Mock(return_value=torch.zeros(1, 3, 3, 2))
    model.s_norm = model.z_norm = nn.Identity()
    model.s_recycle = model.z_recycle = nn.Identity()
    model.msa_module = Mock(return_value=torch.zeros(1, 3, 3, 2))
    model.pairformer_module = Mock(side_effect=lambda s, z, **kwargs: (s, z))
    model.distogram_module = nn.Identity()
    model.diffusion_conditioning = Mock(return_value=(None,) * 6)
    sampled = torch.arange(24, dtype=torch.float32).reshape(2, 4, 3)
    model.structure_module = SimpleNamespace(
        sample=Mock(return_value={"sample_atom_coords": sampled})
    )
    model.confidence_module = Mock(return_value={"iptm": torch.tensor([0.2, 0.8])})
    if ensemble:
        model.affinity_module1 = Head(1)
        model.affinity_module2 = Head(3)
    else:
        model.affinity_module = Head(1)
    return model.eval()


def features():
    return {
        "embedding": torch.ones(1, 3, 3),
        "token_bonds": torch.zeros(1, 3, 3, 1),
        "token_pad_mask": torch.ones(1, 3),
        "mol_type": torch.tensor([[0, 0, 3]]),
        "asym_id": torch.tensor([[0, 0, 1]]),
        "affinity_token_mask": torch.tensor([[0, 0, 1]]),
        "affinity_mw": [100.0],
        "atom_pad_mask": torch.ones(1, 4),
        "atom_resolved_mask": torch.ones(1, 4),
        "coords": torch.arange(12, dtype=torch.float32)
        .reshape(1, 1, 4, 3)
        .requires_grad_(),
    }


@pytest.mark.parametrize("ensemble", [False, True])
@pytest.mark.parametrize("correction", [False, True])
@pytest.mark.parametrize("skip_structure", [False, True])
def test_input_coordinates_reach_heads_without_sampling(
    ensemble, correction, skip_structure
):
    model = make_model(ensemble, correction)
    model.skip_run_structure = skip_structure
    feats = features()
    original = feats["coords"].clone()
    output = model(
        feats,
        recycling_steps=2,
        diffusion_samples=5,
        affinity_from_input_structure=True,
    )
    model.structure_module.sample.assert_not_called()
    model.diffusion_conditioning.assert_not_called()
    model.confidence_module.assert_not_called()
    assert model.pairformer_module.call_count == 3
    assert "iptm" not in output
    assert torch.equal(feats["coords"], original)
    assert torch.equal(output["sample_atom_coords"], original[:, 0])
    heads = (
        [model.affinity_module1, model.affinity_module2]
        if ensemble
        else [model.affinity_module]
    )
    for head in heads:
        assert torch.equal(head.coordinates, original)
        assert not head.coordinates.requires_grad
        assert head.coordinates.data_ptr() == feats["coords"].data_ptr()
    expected = original.square().sum().reshape(1) + (2 if ensemble else 1)
    if ensemble and correction:
        expected = 1.03525938 * expected - 0.59992683 * 100.0**0.3 + 2.83288489
    torch.testing.assert_close(output["affinity_pred_value"], expected)
    assert torch.isfinite(output["affinity_probability_binary"]).all()


@pytest.mark.parametrize("explicit_default", [False, True])
def test_default_route_still_samples_and_selects_highest_iptm(explicit_default):
    model = make_model()
    kwargs = {"affinity_from_input_structure": False} if explicit_default else {}
    output = model(features(), diffusion_samples=2, **kwargs)
    model.structure_module.sample.assert_called_once()
    model.confidence_module.assert_called_once()
    selected = output["sample_atom_coords"][1][None, None]
    assert torch.equal(model.affinity_module1.coordinates, selected)
    assert torch.equal(model.affinity_module2.coordinates, selected)


@pytest.mark.parametrize(
    "problem",
    [
        "training",
        "no_affinity",
        "no_trunk",
        "missing",
        "rank",
        "batch",
        "conformations",
        "xyz",
        "integer",
        "float64",
        "float16",
        "nan",
        "inf",
        "unresolved",
        "atom_count",
        "resolved_shape",
        "empty",
        "no_ligand",
        "no_protein",
        "cofactor",
        "two_ligands",
    ],
)
def test_invalid_input_fails_before_trunk(problem):
    model, feats = make_model(), features()
    if problem == "training":
        model.train()
    elif problem == "no_affinity":
        model.affinity_prediction = False
    elif problem == "no_trunk":
        model.run_trunk_and_structure = False
    elif problem == "missing":
        del feats["coords"]
    elif problem == "rank":
        feats["coords"] = torch.zeros(1, 4, 3)
    elif problem == "batch":
        feats["coords"] = torch.zeros(2, 1, 4, 3)
    elif problem == "conformations":
        feats["coords"] = torch.zeros(1, 2, 4, 3)
    elif problem == "xyz":
        feats["coords"] = torch.zeros(1, 1, 4, 2)
    elif problem == "integer":
        feats["coords"] = torch.zeros(1, 1, 4, 3, dtype=torch.long)
    elif problem == "float64":
        feats["coords"] = torch.zeros(1, 1, 4, 3, dtype=torch.float64)
    elif problem == "float16":
        feats["coords"] = torch.zeros(1, 1, 4, 3, dtype=torch.float16)
    elif problem in {"nan", "inf"}:
        feats["coords"] = torch.full((1, 1, 4, 3), float(problem))
    elif problem == "unresolved":
        feats["atom_resolved_mask"][0, 1] = 0
    elif problem == "atom_count":
        feats["atom_pad_mask"] = torch.ones(1, 3)
    elif problem == "resolved_shape":
        feats["atom_resolved_mask"] = torch.ones(1, 3)
    elif problem == "empty":
        feats["atom_pad_mask"].zero_()
    elif problem == "no_ligand":
        feats["affinity_token_mask"].zero_()
    elif problem == "no_protein":
        feats["mol_type"].fill_(3)
    elif problem == "cofactor":
        feats["mol_type"][0, 0] = 3
    elif problem == "two_ligands":
        feats["affinity_token_mask"][0, 1] = 1
        feats["mol_type"][0, 1] = 3
        feats["asym_id"][0, 1] = 2
    with pytest.raises(ValueError):
        model(feats, affinity_from_input_structure=True)
    model.pairformer_module.assert_not_called()


def test_unresolved_padding_is_allowed():
    model, feats = make_model(), features()
    feats["atom_pad_mask"][0, -1] = 0
    feats["atom_resolved_mask"][0, -1] = 0
    assert torch.isfinite(
        model(feats, affinity_from_input_structure=True)["affinity_pred_value"]
    ).all()
