# SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
# SPDX-License-Identifier: MIT
#
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

"""Tests for the Boltz2 distributed model wrapper.

Verification checks:
    V1: Construction – all serial parameters are present in the distributed wrapper
    V2: Placeholder modules raise NotImplementedError on forward
    V3: Ready submodules produce DTensor outputs with correct shapes
    V4: configure_optimizers returns valid optimizer configuration
    V5: configure_callbacks returns DistributedEMA when EMA is enabled
    V6: on_after_backward redistributes DTensor gradients to Replicate
    V7: Multi-rank construction and parameter identity across ranks
    V8: on_load_checkpoint adjusts checkpoint hyperparameters
    V9: Non-vacuous guards – distributed modules have DTensor params,
        placeholders have plain params
    V10: bf16 mixed precision – ready submodules produce bf16 outputs under
        autocast, gradients are reduced in >=fp32 (dp=1, cp=1x1)
    V11: EmbeddingParamsReplicated wrapping for token_bonds_type
    V12: BFactorModule wrapping when predict_bfactor=True
    V13: Forward/backward parity – distributed Boltz2 forward matches serial
        (restored from dev-v2 for debug comparison with V16)
    V14: predict_step parity – distributed predict_step matches serial
    V14b: predict_step confidence output – 2-GPU smoke test verifying no
        DTensor values leak into the predict output dict (which would crash
        the BoltzWriter callback)
    V16: Serial vs distributed training_step parity – compares loss, gradients,
        post-optimizer parameters, and CSVLogger-captured logged metrics between
        serial Boltz2.training_step and distributed Boltz2Distributed.training_step
        with all 5 randomness sources controlled (recycling, noise, augmentation,
        sampling, dropout).  Parametrized across dp2-cp1x1 (DP only), dp1-cp2x2
        (CP only), and dp2-cp2x2 (DP + CP) to verify correctness under all
        sharding modes.
    V15: setup – validator wiring, predict no-op, datamodule=None
    V17: Regression tests for get_true_coordinates and loss shape consistency –
        verifies non-symmetry path returns DTensors, symmetry_correction raises
        NotImplementedError, and loss zeros have scalar shape ()
    V18: validation_step and on_validation_epoch_end – real forward pass with
         metric accumulation, aggregation, and logging
"""

import math
import random as stdlib_random
import shutil
import tempfile
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import torch
from pytorch_lightning.loggers import CSVLogger
from torch.distributed.tensor import DTensor, Replicate, Shard, distribute_tensor

import boltz.distributed.model.modules.diffusion as distributed_diffusion_module
import boltz.model.modules.diffusionv2 as serial_diffusion_v2_module
from boltz.data import const
from boltz.data.module.trainingv2 import Boltz2TrainingDataModule as Boltz2TrainingDataModuleSerial
from boltz.data.module.trainingv2 import collate
from boltz.distributed.data.module.trainingv2 import Boltz2TrainingDataModule as BoltzTrainingDataModuleDTensor
from boltz.distributed.data.utils import (
    ATOM_FEATURES_V2,
    LIGAND_GEOMETRY_FEATURES,
    distribute_features,
    map_subgroup_mesh_to_cpu,
)
from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.layers.elementwise_op import ElementwiseOp, elementwise_op, scalar_tensor_op
from boltz.distributed.model.models.boltz2 import Boltz2 as Boltz2Distributed
from boltz.distributed.model.models.boltz2 import _PlaceholderModule
from boltz.distributed.model.validation.rcsb import DistributedRCSBValidator
from boltz.distributed.testing.utils import setup_mock_training_datamodule_config
from boltz.model.layers.attention import AttentionPairBias as AttentionPairBiasV1
from boltz.model.layers.pairformer import PairformerLayer as SerialPairformerLayer
from boltz.model.models.boltz2 import Boltz2 as SerialBoltz2
from boltz.model.validation.rcsb import RCSBValidator
from boltz.testing.utils import (
    SetModuleInfValues,
    concat_data,
    create_boltz2_model_init_params,
    distribute_atom_features,
    get_feature_placements,
    init_module_params_glorot,
    init_tensors_uniform,
    pad_to_length,
    random_features,
    seed_by_rank,
    spawn_multiprocessing,
)


class _DictNamespace:
    """A picklable namespace with both attribute access and .get() support."""

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    def get(self, key, default=None):
        return self.__dict__.get(key, default)


class _LogCapture:
    """Captures LightningModule.log() calls into a dict, backed by a CSVLogger.

    Replaces the usual ``lambda *a, **kw: None`` monkeypatch so that the
    training_step logging code path is exercised rather than silenced.
    After the step, :meth:`flush` persists the captured metrics to the CSV file.
    """

    def __init__(self, csv_logger: CSVLogger):
        self.metrics: dict[str, float] = {}
        self._csv_logger = csv_logger

    def __call__(self, name, value, **kwargs):
        if isinstance(value, torch.Tensor):
            value = value.detach().cpu().item()
        self.metrics[name] = value

    def flush(self, step: int = 0) -> None:
        """Write captured metrics to the backing CSVLogger."""
        self._csv_logger.log_metrics(self.metrics, step=step)
        self._csv_logger.save()


def _setup_training_data_7z64_8b2e(out_dir: Path, base_data_dir: Path) -> Path:
    """Merge 7z64 and 8b2e processed data into a single training directory with records.

    Uses the two smallest samples (8b2e=1062 atoms, 7z64=2278 atoms) to reduce
    GPU memory consumption in FP64 parity tests.
    """
    names = ["7z64", "8b2e"]
    source_dirs = [base_data_dir / f"processed_{name}" for name in names]
    merged = concat_data(out_dir, *source_dirs)
    records_dir = merged / "records"
    records_dir.mkdir(parents=True, exist_ok=True)
    copied: set[str] = set()
    for src in source_dirs:
        for rf in (src / "records").glob("*.json"):
            if rf.name in copied:
                raise ValueError(f"Duplicate record file {rf.name}")
            shutil.copy(rf, records_dir / rf.name)
            copied.add(rf.name)
    return merged


def _deterministic_getitem_monkeypatch(monkeypatch, dataset, base_seed=42):
    """Wrap TrainingDataset.__getitem__ to seed all RNGs per idx.

    Ensures idx=0 picks sample 0 and idx=1 picks sample 1 by intercepting
    np.random.choice calls for dataset and sample selection.  Seeds np.random,
    torch, and Python random with base_seed + idx before each call so that all
    downstream RNG usage (cropper, featurizer, center_random_augmentation) is
    deterministic for a given sample index.
    """
    original_getitem = type(dataset).__getitem__

    def _wrapped_getitem(self, idx):
        np.random.seed(base_seed + idx)
        torch.manual_seed(base_seed + idx)
        stdlib_random.seed(base_seed + idx)

        _original_np_choice = np.random.choice
        _call_count = [0]
        _num_samples = len(self.samples[0])

        def _deterministic_choice(a, p=None, **kwargs):
            _call_count[0] += 1
            result = _original_np_choice(a, p=p, **kwargs)
            if _call_count[0] == 1:
                return 0
            elif _call_count[0] == 2:
                return idx % _num_samples
            return result

        np.random.choice = _deterministic_choice
        try:
            return original_getitem(self, idx)
        finally:
            np.random.choice = _original_np_choice

    monkeypatch.setattr(type(dataset), "__getitem__", _wrapped_getitem)


def _make_training_args(**overrides):
    """Create minimal training_args for Boltz2."""
    defaults = {
        "recycling_steps": 1,
        "sampling_steps": 2,
        "diffusion_multiplicity": 1,
        "diffusion_samples": 1,
        "diffusion_loss_weight": 1.0,
        "distogram_loss_weight": 0.3,
        "confidence_loss_weight": 0.0,
        "bfactor_loss_weight": 0.0,
        "symmetry_correction": False,
        "adam_beta_1": 0.9,
        "adam_beta_2": 0.95,
        "adam_eps": 1e-8,
        "base_lr": 1e-3,
        "max_lr": 1e-3,
        "lr_scheduler": "af3",
        "lr_warmup_no_steps": 10,
        "lr_start_decay_after_n_steps": 100,
        "lr_decay_every_n_steps": 50000,
        "lr_decay_factor": 0.95,
        "weight_decay": 0.0,
    }
    defaults.update(overrides)
    return _DictNamespace(**defaults)


def _make_validation_args(**overrides):
    defaults = {
        "recycling_steps": 0,
        "sampling_steps": 2,
        "diffusion_samples": 1,
        "symmetry_correction": False,
        "run_confidence_sequentially": False,
    }
    defaults.update(overrides)
    return _DictNamespace(**defaults)


TOKEN_S = 32
TOKEN_Z = 16

_BOLTZ2_SELECTED_KEYS = [
    "atom_pad_mask",
    "atom_to_token",
    "pair_mask",
    "token_pad_mask",
    "ref_pos",
    "ref_charge",
    "ref_element",
    "ref_atom_name_chars",
    "ref_space_uid",
    "res_type",
    "profile",
    "deletion_mean",
    "pocket_feature",
    "atom_resolved_mask",
    "mol_type",
    "msa",
    "has_deletion",
    "deletion_value",
    "msa_paired",
    "msa_mask",
    "token_bonds",
    "type_bonds",
    "token_pair_pad_mask",
    "asym_id",
    "residue_index",
    "entity_id",
    "token_index",
    "sym_id",
    "cyclic_period",
    "coords",
    "disto_target",
    "token_disto_mask",
    "atom_counts_per_token",
    "token_to_rep_atom",
    "frames_idx",
    "contact_conditioning",
    "contact_threshold",
    "method_feature",
    "modified",
    "bfactor",
    "plddt",
]


def _create_minimal_serial_boltz2(
    confidence_prediction=False,
    affinity_prediction=False,
    ema=False,
    bond_type_feature=False,
    predict_bfactor=False,
    validate_structure=False,
    validators=None,
    num_val_datasets=1,
):
    """Create a minimal serial Boltz2 model for testing."""
    training_args = _make_training_args()
    validation_args = _make_validation_args()

    pairformer_args = {"num_blocks": 1, "num_heads": 2, "dropout": 0.0}

    model = SerialBoltz2(
        atom_s=16,
        atom_z=8,
        token_s=TOKEN_S,
        token_z=TOKEN_Z,
        num_bins=8,
        training_args=training_args,
        validation_args=validation_args,
        embedder_args={
            "atom_encoder_depth": 1,
            "atom_encoder_heads": 2,
            "activation_checkpointing": False,
        },
        msa_args={"msa_s": 16, "msa_blocks": 1, "msa_dropout": 0.0, "z_dropout": 0.0},
        pairformer_args=pairformer_args,
        score_model_args={
            "sigma_data": 16.0,
            "dim_fourier": 32,
            "atom_encoder_depth": 1,
            "atom_encoder_heads": 2,
            "token_transformer_depth": 1,
            "token_transformer_heads": 2,
            "atom_decoder_depth": 1,
            "atom_decoder_heads": 2,
            "activation_checkpointing": False,
            "conditioning_transition_layers": 1,
        },
        diffusion_process_args={"num_sampling_steps": 2},
        diffusion_loss_args={},
        confidence_prediction=confidence_prediction,
        confidence_model_args=(
            {"pairformer_args": pairformer_args, "confidence_args": {}} if confidence_prediction else None
        ),
        affinity_prediction=affinity_prediction,
        predict_args={"recycling_steps": 0, "sampling_steps": 2, "diffusion_samples": 1, "max_parallel_samples": 1},
        validate_structure=validate_structure,
        structure_prediction_training=True,
        ema=ema,
        use_templates=False,
        predict_bfactor=predict_bfactor,
        bond_type_feature=bond_type_feature,
        validators=validators if validate_structure else None,
        num_val_datasets=num_val_datasets if validate_structure else 1,
    )

    return model


def _prepare_serial_model(ema=False, bond_type_feature=False, predict_bfactor=False):
    """Create a serial model and return its state dict and hparams."""
    model = _create_minimal_serial_boltz2(ema=ema, bond_type_feature=bond_type_feature, predict_bfactor=predict_bfactor)
    return model.state_dict(), dict(model.hparams)


def _init_distributed(rank, grid_group_sizes, device_type, backend, env_map):
    """Common boilerplate: set env vars and initialize DistributedManager."""
    monkeypatch = pytest.MonkeyPatch()
    if env_map is not None:
        for var_name, value in env_map.items():
            monkeypatch.setenv(var_name, f"{rank}" if value == "<INPUT_RANK>" else value)
    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    return monkeypatch, DistributedManager()


def _build_dist_model(serial_hparams, serial_state_dict, dist_manager):
    """Create serial model, load state, wrap as distributed, move to device."""
    serial_model = SerialBoltz2(**serial_hparams)
    serial_model.load_state_dict(serial_state_dict, strict=True)
    serial_model = serial_model.to(dist_manager.device)
    dist_model = Boltz2Distributed(serial_model, dist_manager)
    return dist_model.to(dist_manager.device)


# ====================================================================== #
#  Single-process unit tests (no distributed setup needed)               #
# ====================================================================== #


def test_wrapper_preserves_all_serial_parameters():
    """V1: All serial parameters should appear in the distributed wrapper."""
    serial_model = _create_minimal_serial_boltz2()
    serial_param_names = set(serial_model.state_dict().keys())
    assert len(serial_param_names) > 0, "Serial model has no parameters"

    expected_prefixes = {
        "s_init",
        "z_init_1",
        "z_init_2",
        "s_norm",
        "z_norm",
        "s_recycle",
        "z_recycle",
        "token_bonds",
        "msa_module",
        "pairformer_module",
        "distogram_module",
        "input_embedder",
        "rel_pos",
        "contact_conditioning",
        "diffusion_conditioning",
        "structure_module",
    }
    found_prefixes = {n.split(".")[0] for n in serial_param_names}
    missing = expected_prefixes - found_prefixes
    assert not missing, f"Expected parameter prefixes missing from serial model: {missing}"


def test_wrapper_hparams_saved():
    """V1b: Hyper-parameters should be preserved in save_hyperparameters."""
    serial_model = _create_minimal_serial_boltz2()
    assert hasattr(serial_model, "hparams"), "Serial model should have hparams"
    assert "atom_s" in serial_model.hparams
    assert "token_s" in serial_model.hparams


def test_placeholder_raises_on_forward():
    """V2: Placeholder module should raise NotImplementedError."""
    placeholder = _PlaceholderModule(torch.nn.Linear(4, 8), "TestModule")
    with pytest.raises(NotImplementedError, match="TestModule"):
        placeholder(torch.randn(2, 4))


def test_placeholder_preserves_parameters():
    """V2b: Placeholder should expose the serial module's parameters, frozen."""
    placeholder = _PlaceholderModule(torch.nn.Linear(4, 8), "TestModule")
    params = dict(placeholder.named_parameters())
    assert "_serial.weight" in params
    assert "_serial.bias" in params
    assert params["_serial.weight"].shape == (8, 4)
    for name, p in placeholder.named_parameters():
        assert not p.requires_grad, f"Placeholder param '{name}' should be frozen"


def test_serial_model_ema_callbacks():
    """V5 precondition: serial EMA callback configuration."""
    assert _create_minimal_serial_boltz2(ema=True).use_ema is True
    assert len(_create_minimal_serial_boltz2(ema=True).configure_callbacks()) == 1
    assert _create_minimal_serial_boltz2(ema=False).use_ema is False
    assert len(_create_minimal_serial_boltz2(ema=False).configure_callbacks()) == 0


def test_confidence_model_can_be_created():
    """Confidence serial model can be created with proper args."""
    assert _create_minimal_serial_boltz2(confidence_prediction=True).confidence_prediction is True


def test_bfactor_in_atom_features_v2():
    """Data pipeline: bfactor must be in ATOM_FEATURES_V2 for distributed sharding."""
    assert "bfactor" in ATOM_FEATURES_V2


def test_pairformer_v1_v2_detection_preconditions():
    """V1/V2 detection: serial V1 attention has norm_s, V2 does not.

    The distributed PairformerLayer infers v1/v2 from
    ``hasattr(layer.attention, 'norm_s')``.  This test verifies the
    serial-side preconditions that the heuristic depends on.
    """
    token_s, token_z, num_heads = 32, 16, 2

    v1_layer = SerialPairformerLayer(token_s, token_z, num_heads, v2=False)
    v2_layer = SerialPairformerLayer(token_s, token_z, num_heads, v2=True)

    assert isinstance(v1_layer.attention, AttentionPairBiasV1)
    assert hasattr(v1_layer.attention, "norm_s"), "V1 attention must have norm_s"
    assert not hasattr(v2_layer.attention, "norm_s"), "V2 attention must not have norm_s"

    # The heuristic: not hasattr(layer.attention, "norm_s") → v2
    assert (not hasattr(v1_layer.attention, "norm_s")) is False  # → V1
    assert (not hasattr(v2_layer.attention, "norm_s")) is True  # → V2


def test_checkpoint_lr_is_overwritten():
    """V8: on_load_checkpoint should overwrite lr and weight_decay in checkpoint."""
    serial_model = _create_minimal_serial_boltz2()
    checkpoint = {
        "optimizer_states": [{"param_groups": [{"lr": 0.5, "weight_decay": 0.99}, {"lr": 0.5, "weight_decay": 0.99}]}],
        "lr_schedulers": [{"max_lr": 0.5, "base_lrs": [0.5, 0.5], "_last_lr": [0.5, 0.5]}],
        "hyper_parameters": {
            "training_args": {"max_lr": 0.5, "diffusion_multiplicity": 99, "recycling_steps": 99, "weight_decay": 0.99}
        },
    }
    serial_model.on_load_checkpoint(checkpoint)

    for pg in checkpoint["optimizer_states"][0]["param_groups"]:
        assert pg["lr"] == 1e-3 and pg["weight_decay"] == 0.0

    sched = checkpoint["lr_schedulers"][0]
    assert sched["max_lr"] == 1e-3
    assert all(lr == 1e-3 for lr in sched["base_lrs"])

    hp = checkpoint["hyper_parameters"]["training_args"]
    assert hp["max_lr"] == 1e-3 and hp["diffusion_multiplicity"] == 1 and hp["weight_decay"] == 0.0


# ====================================================================== #
#  Comprehensive multi-rank worker (V3-V7, V9, V11)                      #
# ====================================================================== #


# V13 (forward/backward parity) was removed — V16 (training_step parity) now
# includes dp1-cp2x2 and dp2-cp2x2 parametrizations that provide strictly
# stronger guarantees (loss + gradient + post-optimizer parity) at the same
# topologies.  This matches the Boltz-1 test structure where
# test_boltz1_model_parallel_training_step with cp2x2 is the primary parity test.


# ====================================================================== #
#  V10: bf16 mixed precision (CUDA-only)                                 #
# ====================================================================== #


def _worker_bf16_mixed_precision(
    rank: int,
    serial_state_dict: dict,
    serial_hparams: dict,
    grid_group_sizes: dict,
    device_type: str,
    backend: str,
    env_map: dict[str, str] | None = None,
):
    """Worker: exercise ready submodules under torch.autocast(bf16).

    Verifies forward outputs are bf16, gradients are reduced in >=fp32,
    and clear_autocast_cache() does not error.
    """
    monkeypatch, dm = _init_distributed(rank, grid_group_sizes, device_type, backend, env_map)

    dist_model = _build_dist_model(serial_hparams, serial_state_dict, dm)
    dist_model.train()

    dp = grid_group_sizes["dp"]
    B = max(1, dp)
    N_global = 8 * grid_group_sizes["cp"][0]
    single_pl = (Shard(0), Shard(1), Replicate())
    pair_pl = (Shard(0), Shard(1), Shard(2))

    # Forward under autocast
    with torch.autocast("cuda", dtype=torch.bfloat16):
        s_in = distribute_tensor(
            torch.randn(B, N_global, TOKEN_S, device=dm.device), dm.device_mesh_subgroups, single_pl
        )
        s_out = dist_model.s_init(s_in)
        assert s_out.to_local().dtype == torch.bfloat16

        s_norm_out = dist_model.s_norm(s_in)
        assert s_norm_out.to_local().dtype == torch.float32  # autocast LayerNorm policy

        z1 = dist_model.z_init_1(s_in)
        assert z1.to_local().dtype == torch.bfloat16

        z_in = distribute_tensor(
            torch.randn(B, N_global, N_global, TOKEN_Z, device=dm.device), dm.device_mesh_subgroups, pair_pl
        )
        disto = dist_model.distogram_module(z_in)
        assert disto.to_local().dtype == torch.bfloat16

    # Backward: gradient dtype and placement
    s_grad_in = distribute_tensor(
        torch.randn(B, N_global, TOKEN_S, device=dm.device), dm.device_mesh_subgroups, single_pl
    )
    with torch.autocast("cuda", dtype=torch.bfloat16):
        out = dist_model.s_init(s_grad_in)
        out.to_local().sum().backward()

    w = dist_model.s_init.weight
    assert w.grad is not None and isinstance(w.grad, DTensor)
    for p in w.grad.placements:
        assert isinstance(p, Replicate), f"Weight grad should be Replicate, got {p}"
    assert w.grad.to_local().dtype == w.to_local().dtype  # fp32 via promote_types
    assert w.grad.to_local().abs().sum() > 0

    # clear_autocast_cache branch + recycling path
    dist_model.zero_grad()
    with torch.autocast("cuda", dtype=torch.bfloat16):
        torch.clear_autocast_cache()
        s_zeros = distribute_tensor(
            torch.zeros(B, N_global, TOKEN_S, device=dm.device, dtype=torch.bfloat16),
            dm.device_mesh_subgroups,
            single_pl,
        )
        s_recycled = elementwise_op(
            dist_model.s_init(s_grad_in),
            dist_model.s_recycle(dist_model.s_norm(s_zeros)),
            ElementwiseOp.SUM,
        )
        assert s_recycled.to_local().dtype == torch.bfloat16

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env",
    [
        # Minimal dp=1, cp=1x1: focused autocast dtype check.
        # Full BF16 training+DP+CP is covered by test_boltz2_finetune_from_checkpoint
        # in test_dtensor_boltz2_train.py and test_boltz2_run_predict in
        # test_dtensor_predict.py.
        ((1, (1, 1)), True, "cuda", "ENV"),
    ],
    indirect=("setup_env",),
    ids=["cuda-dp1-cp1x1"],
)
def test_boltz2_bf16_mixed_precision(setup_env):
    """V10: Ready submodules produce bf16 outputs under autocast, grads reduced in fp32.

    Focused autocast check on individual submodules (s_init, z_init_1,
    distogram_module, s_norm) with dp=1, cp=1x1.  Full BF16-mixed training
    across DP and CP topologies is exercised by
    ``test_boltz2_finetune_from_checkpoint`` and ``test_boltz2_run_predict``.
    """
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if not torch.cuda.is_available():
        pytest.skip("CUDA not available")
    if torch.cuda.device_count() < world_size:
        pytest.skip(f"Need {world_size} GPUs, have {torch.cuda.device_count()}")

    sd, hp = _prepare_serial_model(ema=False)
    spawn_multiprocessing(
        _worker_bf16_mixed_precision,
        world_size,
        sd,
        hp,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
    )


# ====================================================================== #
#  V14: predict_step parity (inference sampling)                          #
# ====================================================================== #


def parallel_assert_boltz2_model_predict_step(
    rank,
    grid_group_sizes,
    device_type,
    backend,
    dtype,
    boltz2_model_params,
    module_state_dict,
    predict_args,
    diffusion_samples,
    num_sampling_steps,
    input_feats_global_fp64_host,
    init_noise_global_host,
    step_noise_list_global_host,
    serial_coords_host,
    serial_masks_host,
    env_per_rank=None,
):
    """V14 multi-rank worker: verify distributed predict_step matches serial."""
    monkeypatch = pytest.MonkeyPatch()
    if env_per_rank is not None:
        for var_name, value in env_per_rank.items():
            monkeypatch.setenv(var_name, f"{rank}" if value == "<INPUT_RANK>" else value)

    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    manager = DistributedManager()

    reference_module = SerialBoltz2(**boltz2_model_params)
    reference_module = reference_module.to(dtype=dtype)
    reference_module.load_state_dict(module_state_dict)
    reference_module.structure_module.coordinate_augmentation = False
    reference_module.apply(SetModuleInfValues())
    reference_module = reference_module.to(device=manager.device)
    module = Boltz2Distributed(reference_module, manager)
    module.eval()

    host_tensor_keys = {k for k, v in input_feats_global_fp64_host.items() if isinstance(v, torch.Tensor)}
    _placements = get_feature_placements(
        token_keys=host_tensor_keys,
        msa_keys=host_tensor_keys,
        atom_keys={
            "ref_pos",
            "atom_resolved_mask",
            "ref_element",
            "ref_charge",
            "ref_atom_name_chars",
            "ref_space_uid",
            "coords",
            "atom_pad_mask",
            "atom_to_token",
            "pair_mask",
            "atom_counts_per_token",
            "token_to_rep_atom",
            "bfactor",
            "plddt",
        },
        model_io_keys={"noise"},
        model_io_fp32_keys=set(),
    )
    _placements_token_features = _placements["token_features"]
    _placements_msa_features = _placements["msa_features"]
    _placements_cp_atom_features = _placements["cp_atom_features"]
    _placements_atom_features = _placements["atom_features"]
    _placements_cp_model_io = _placements["cp_model_io"]
    _placements_model_io = _placements["model_io"]

    # ------------------------------------------------------------------
    # Distribute token + MSA features (rank 0 broadcasts)
    # ------------------------------------------------------------------
    if manager.group_rank["world"] == 0:
        input_feats_token_msa_global = {
            k: v.to(device=manager.device, dtype=dtype if v.dtype.is_floating_point else v.dtype)
            for k, v in input_feats_global_fp64_host.items()
            if k in _placements_token_features or k in _placements_msa_features
        }
    else:
        input_feats_token_msa_global = None

    feats_token_msa = distribute_features(
        input_feats_token_msa_global,
        _placements_token_features | _placements_msa_features,
        manager.group["world"],
        manager.group_ranks["world"][0],
        manager.device_mesh_subgroups,
    )

    # ------------------------------------------------------------------
    # Distribute atom features + sampling noise via distribute_atom_features.
    # Noise tensors use the _noise_{i_noise}_{i_mul} naming convention
    # from test_atom_diffusion_sample so that intersperse padding naturally
    # places zeros at padding positions.
    # ------------------------------------------------------------------
    size_batch = input_feats_global_fp64_host["atom_pad_mask"].shape[0]
    inputs_atom = {
        k: v.to(dtype=dtype if v.dtype.is_floating_point else v.dtype)
        for k, v in input_feats_global_fp64_host.items()
        if k in _placements_cp_atom_features
    }

    all_noise = [init_noise_global_host] + list(step_noise_list_global_host)
    for i_noise, noise_host in enumerate(all_noise):
        unflat = noise_host.unflatten(0, (size_batch, diffusion_samples))
        for i_mul in range(diffusion_samples):
            inputs_atom[f"_noise_{i_noise}_{i_mul}"] = unflat[:, i_mul].to(dtype=dtype)

    noise_cp_placements = {}
    noise_placements = {}
    for i_noise in range(len(all_noise)):
        for i_mul in range(diffusion_samples):
            key = f"_noise_{i_noise}_{i_mul}"
            noise_cp_placements[key] = _placements_cp_model_io["noise"]
            noise_placements[key] = _placements_model_io["noise"]

    feats_and_noise = distribute_atom_features(
        inputs=inputs_atom,
        placements_cp=_placements_cp_atom_features | noise_cp_placements,
        placements_dp_cp=_placements_atom_features | noise_placements,
        device_mesh=manager.device_mesh_subgroups,
        cp_group=manager.group["cp"],
        multiplicities={f"_noise_{i}": diffusion_samples for i in range(len(all_noise))},
    )

    noise_dts = []
    for i_noise in range(len(all_noise)):
        noise_dts.append(feats_and_noise.pop(f"_noise_{i_noise}"))
    init_noise_dt = noise_dts[0]
    step_noise_dts = noise_dts[1:]

    feats_dt = {**feats_token_msa, **feats_and_noise}

    # ------------------------------------------------------------------
    # Monkeypatch distributed sample() for determinism
    # ------------------------------------------------------------------
    _orig_center_random_augmentation = distributed_diffusion_module.center_random_augmentation

    def _centering_only_augmentation(atom_coords, atom_mask, **kwargs):
        kwargs["augmentation"] = False
        kwargs["centering"] = True
        return _orig_center_random_augmentation(atom_coords, atom_mask, **kwargs)

    _dt_randn_calls = []
    _dt_randn_sequence = [init_noise_dt] + step_noise_dts

    def _fixed_create_distributed_randn(shape, device_mesh, placements, dtype=torch.float32, scale=1.0):
        idx = len(_dt_randn_calls)
        _dt_randn_calls.append(idx)
        noise_dt = _dt_randn_sequence[idx]
        if scale != 1.0:
            noise_dt = scalar_tensor_op(scale, noise_dt, ElementwiseOp.PROD)
        return noise_dt

    monkeypatch.setattr(distributed_diffusion_module, "center_random_augmentation", _centering_only_augmentation)
    monkeypatch.setattr(distributed_diffusion_module, "create_distributed_randn", _fixed_create_distributed_randn)

    # ------------------------------------------------------------------
    # Run distributed predict_step
    # ------------------------------------------------------------------
    module.predict_args = predict_args
    with torch.no_grad():
        pred_dict = module.predict_step(feats_dt, batch_idx=0)

    assert pred_dict["exception"] is False

    # ------------------------------------------------------------------
    # Compare on gather rank 0 of CP axis 0.
    # The distributed predict_step gathers coords/masks via
    # torch.distributed.gather + concat(dim=1). With dp>1 each DP rank
    # gathers only its own batch element(s), so we slice the serial
    # reference by DP rank.
    # ------------------------------------------------------------------
    tag_group_gather = 0
    if manager.subgroups_rank["cp"][tag_group_gather] == 0:
        gathered_mask = pred_dict["masks"]
        gathered_coords = pred_dict["coords"]

        mask_expanded = gathered_mask.repeat_interleave(diffusion_samples, 0).bool()
        dt_real = gathered_coords[mask_expanded]

        n_dp = grid_group_sizes["dp"]
        dp_rank = manager.group_rank["dp"]
        B_local = size_batch // n_dp
        M = diffusion_samples
        serial_coords_slice = serial_coords_host[dp_rank * B_local * M : (dp_rank + 1) * B_local * M]
        serial_mask_slice = serial_masks_host[dp_rank * B_local : (dp_rank + 1) * B_local]

        serial_coords_device = serial_coords_slice.to(device=manager.device, dtype=dtype)
        serial_mask_expanded = serial_mask_slice.to(device=manager.device).repeat_interleave(M, 0).bool()
        serial_real = serial_coords_device[serial_mask_expanded]

        torch.testing.assert_close(dt_real, serial_real)

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env",
    [
        # dp=2, CP=(2,2) on CUDA — exercises DP slicing + CP gather
        ((2, (2, 2)), True, "cuda", "ENV"),
    ],
    indirect=("setup_env",),
    ids=lambda val: f"{val[2]}-dp:{val[0][0]}-cp{'x'.join(map(str, val[0][1]))}",
)
def test_boltz2_predict_step(setup_env):
    """V14: predict_step parity between distributed and serial Boltz2.

    Tests that the distributed Boltz2.predict_step produces numerically
    identical sampled coordinates compared to the serial implementation
    in eval mode with no backward pass.

    Side-by-side comparison findings (serial vs distributed):

    Bug 1 — Random augmentation always applied in serial sample():
      Serial diffusionv2.py:370-376 calls compute_random_augmentation()
      unconditionally.  Distributed diffusion.py:982-994 gates
      center_random_augmentation() on self.coordinate_augmentation.
      Mitigation: monkeypatch serial compute_random_augmentation to return
      identity rotation and zero translation (same as test_atom_diffusion_sample).

    Bug 2 — alignment_reverse_diff FP32 downcast in serial:
      Serial diffusionv2.py:564-573 forces .float() (FP32) for
      weighted_rigid_align.  Distributed diffusion.py:1076-1082 passes
      DTensors as-is (FP64 in test).
      Mitigation: set alignment_reverse_diff=False.

    Bug 3 — Serial sample() augmentation shape bug for B > 1:
      compute_random_augmentation(multiplicity) returns R of shape (M,3,3)
      but atom_coords has shape (B*M,N,3) — einsum crashes when B > 1.
      Mitigation: the identity augmentation mock returns (B*M,3,3), making
      the einsum a no-op for any batch size.
    """
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("CUDA not available")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"Need {world_size} GPUs, have {torch.cuda.device_count()}")

    dtype = torch.float64
    size_batch_per_rank = 1
    B = size_batch_per_rank * grid_group_sizes["dp"]
    size_cp = grid_group_sizes["cp"][0]

    min_val_init = -0.01
    max_val_init = 0.01
    scale_glorot = 0.05

    num_sampling_steps = 2
    diffusion_samples = 2

    seed = 42
    seed_by_rank(0, seed=seed)

    boltz2_model_params = create_boltz2_model_init_params(use_large_model=False)
    boltz2_model_params["diffusion_process_args"]["alignment_reverse_diff"] = False

    predict_args = {
        "recycling_steps": 0,
        "sampling_steps": num_sampling_steps,
        "diffusion_samples": diffusion_samples,
        "max_parallel_samples": diffusion_samples,
        "write_confidence_summary": False,
        "write_full_pae": False,
    }
    boltz2_model_params["predict_args"] = predict_args

    n_atoms_per_token_min = 8
    n_atoms_per_token_max = 20
    n_tokens = 30 * size_cp
    W = boltz2_model_params["atoms_per_window_queries"]
    n_atoms_raw = n_tokens * n_atoms_per_token_max
    n_atoms = ((n_atoms_raw + W - 1) // W) * W
    n_msa = size_cp * 2

    assert n_atoms % size_cp == 0
    assert n_atoms % W == 0

    input_feats_global_fp64 = random_features(
        size_batch=B,
        n_tokens=n_tokens,
        n_atoms=n_atoms,
        n_msa=n_msa,
        atom_counts_per_token_range=(n_atoms_per_token_min, n_atoms_per_token_max),
        device=device_type,
        float_value_range=(min_val_init, max_val_init),
        selected_keys=_BOLTZ2_SELECTED_KEYS,
        num_disto_bins=boltz2_model_params["num_bins"],
    )
    input_feats_global_fp64["msa"] = torch.randint(
        0, const.num_tokens, (B, n_msa, n_tokens), dtype=torch.int64, device=device_type
    )

    # ------------------------------------------------------------------
    # Build serial model in eval mode
    # ------------------------------------------------------------------
    reference_module = SerialBoltz2(**boltz2_model_params)
    init_module_params_glorot(reference_module, gain=scale_glorot)
    reference_module.apply(SetModuleInfValues())
    reference_module.structure_module.coordinate_augmentation = False
    module_state_dict_fp64 = reference_module.state_dict()
    reference_module = reference_module.to(dtype=torch.float64, device=device_type).eval()

    # ------------------------------------------------------------------
    # Pre-generate deterministic noise for sampling.
    # sample_schedule(N) returns N+1 sigmas (F.pad adds trailing 0), so
    # sigmas_and_gammas has N entries → N denoising steps.  Total
    # torch.randn calls = 1 (init) + N (per-step eps) = N+1.
    # ------------------------------------------------------------------
    _B_M = B * diffusion_samples
    init_noise = torch.empty((_B_M, n_atoms, 3), device=device_type, dtype=dtype)
    step_noise_list = [
        torch.empty((_B_M, n_atoms, 3), device=device_type, dtype=dtype) for _ in range(num_sampling_steps)
    ]
    init_tensors_uniform([init_noise, *step_noise_list], low=min_val_init, high=max_val_init)

    # ------------------------------------------------------------------
    # Monkeypatch serial sample() for determinism (Bug 1/3 mitigations)
    # ------------------------------------------------------------------
    def _identity_compute_random_augmentation(multiplicity_arg, device=None, dtype=None):
        R = torch.eye(3, device=device, dtype=dtype).unsqueeze(0).expand(_B_M, -1, -1)
        tr = torch.zeros(_B_M, 1, 3, device=device, dtype=dtype)
        return R, tr

    _serial_randn_calls = []
    _serial_randn_sequence = [init_noise] + step_noise_list

    def _fixed_randn(*args, **kwargs):
        idx = len(_serial_randn_calls)
        _serial_randn_calls.append(idx)
        return _serial_randn_sequence[idx].clone()

    _monkeypatch = pytest.MonkeyPatch()
    _monkeypatch.setattr(
        serial_diffusion_v2_module, "compute_random_augmentation", _identity_compute_random_augmentation
    )
    _monkeypatch.setattr(serial_diffusion_v2_module.torch, "randn", _fixed_randn)

    # ------------------------------------------------------------------
    # Serial predict_step
    # ------------------------------------------------------------------
    with torch.no_grad():
        serial_pred_dict = reference_module.predict_step(input_feats_global_fp64, batch_idx=0)

    assert serial_pred_dict["exception"] is False
    serial_coords = serial_pred_dict["coords"]
    serial_masks = serial_pred_dict["masks"]

    _monkeypatch.undo()

    # ------------------------------------------------------------------
    # Move everything to CPU for spawn_multiprocessing
    # ------------------------------------------------------------------
    input_feats_host = {k: v.detach().to(device="cpu", copy=True) for k, v in input_feats_global_fp64.items()}
    serial_coords_host = serial_coords.detach().to(device="cpu", copy=True)
    serial_masks_host = serial_masks.detach().to(device="cpu", copy=True)

    spawn_multiprocessing(
        parallel_assert_boltz2_model_predict_step,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        dtype,
        boltz2_model_params,
        module_state_dict_fp64,
        predict_args,
        diffusion_samples,
        num_sampling_steps,
        input_feats_host,
        init_noise.cpu(),
        [n.cpu() for n in step_noise_list],
        serial_coords_host,
        serial_masks_host,
        env_per_rank,
    )


# ====================================================================== #
#  V14b: predict_step confidence output – no DTensor leaks                #
# ====================================================================== #


def _worker_predict_step_confidence(
    rank,
    grid_group_sizes,
    device_type,
    backend,
    dtype,
    boltz2_model_params,
    module_state_dict,
    predict_args,
    diffusion_samples,
    num_sampling_steps,
    input_feats_global_fp64_host,
    init_noise_global_host,
    step_noise_list_global_host,
    env_per_rank=None,
):
    """V14b worker: verify predict_step with confidence produces no DTensor leaks."""
    monkeypatch = pytest.MonkeyPatch()
    if env_per_rank is not None:
        for var_name, value in env_per_rank.items():
            monkeypatch.setenv(var_name, f"{rank}" if value == "<INPUT_RANK>" else value)

    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    manager = DistributedManager()

    reference_module = SerialBoltz2(**boltz2_model_params)
    reference_module = reference_module.to(dtype=dtype)
    reference_module.load_state_dict(module_state_dict)
    reference_module.structure_module.coordinate_augmentation = False
    reference_module.apply(SetModuleInfValues())
    reference_module = reference_module.to(device=manager.device)
    module = Boltz2Distributed(reference_module, manager)
    module.eval()

    host_tensor_keys = {k for k, v in input_feats_global_fp64_host.items() if isinstance(v, torch.Tensor)}
    _placements = get_feature_placements(
        token_keys=host_tensor_keys,
        msa_keys=host_tensor_keys,
        atom_keys={
            "ref_pos",
            "atom_resolved_mask",
            "ref_element",
            "ref_charge",
            "ref_atom_name_chars",
            "ref_space_uid",
            "coords",
            "atom_pad_mask",
            "atom_to_token",
            "pair_mask",
            "atom_counts_per_token",
            "token_to_rep_atom",
            "bfactor",
            "plddt",
        },
        model_io_keys={"noise"},
        model_io_fp32_keys=set(),
    )
    _placements_token_features = _placements["token_features"]
    _placements_msa_features = _placements["msa_features"]
    _placements_cp_atom_features = _placements["cp_atom_features"]
    _placements_cp_model_io = _placements["cp_model_io"]
    _placements_model_io = _placements["model_io"]

    if manager.group_rank["world"] == 0:
        input_feats_token_msa_global = {
            k: v.to(device=manager.device, dtype=dtype if v.dtype.is_floating_point else v.dtype)
            for k, v in input_feats_global_fp64_host.items()
            if k in _placements_token_features or k in _placements_msa_features
        }
    else:
        input_feats_token_msa_global = None

    feats_token_msa = distribute_features(
        input_feats_token_msa_global,
        _placements_token_features | _placements_msa_features,
        manager.group["world"],
        manager.group_ranks["world"][0],
        manager.device_mesh_subgroups,
    )

    size_batch = input_feats_global_fp64_host["atom_pad_mask"].shape[0]
    inputs_atom = {
        k: v.to(dtype=dtype if v.dtype.is_floating_point else v.dtype)
        for k, v in input_feats_global_fp64_host.items()
        if k in _placements_cp_atom_features
    }

    all_noise = [init_noise_global_host] + list(step_noise_list_global_host)
    for i_noise, noise_host in enumerate(all_noise):
        unflat = noise_host.unflatten(0, (size_batch, diffusion_samples))
        for i_mul in range(diffusion_samples):
            inputs_atom[f"_noise_{i_noise}_{i_mul}"] = unflat[:, i_mul].to(dtype=dtype)

    noise_cp_placements = {}
    noise_placements = {}
    for i_noise in range(len(all_noise)):
        for i_mul in range(diffusion_samples):
            key = f"_noise_{i_noise}_{i_mul}"
            noise_cp_placements[key] = _placements_cp_model_io["noise"]
            noise_placements[key] = _placements_model_io["noise"]

    feats_and_noise = distribute_atom_features(
        inputs=inputs_atom,
        placements_cp=_placements_cp_atom_features | noise_cp_placements,
        placements_dp_cp=_placements["atom_features"] | noise_placements,
        device_mesh=manager.device_mesh_subgroups,
        cp_group=manager.group["cp"],
        multiplicities={f"_noise_{i}": diffusion_samples for i in range(len(all_noise))},
    )

    noise_dts = []
    for i_noise in range(len(all_noise)):
        noise_dts.append(feats_and_noise.pop(f"_noise_{i_noise}"))
    init_noise_dt = noise_dts[0]
    step_noise_dts = noise_dts[1:]

    feats_dt = {**feats_token_msa, **feats_and_noise}

    _orig_center_random_augmentation = distributed_diffusion_module.center_random_augmentation

    def _centering_only_augmentation(atom_coords, atom_mask, **kwargs):
        kwargs["augmentation"] = False
        kwargs["centering"] = True
        return _orig_center_random_augmentation(atom_coords, atom_mask, **kwargs)

    _dt_randn_calls = []
    _dt_randn_sequence = [init_noise_dt] + step_noise_dts

    def _fixed_create_distributed_randn(shape, device_mesh, placements, dtype=torch.float32, scale=1.0):
        idx = len(_dt_randn_calls)
        _dt_randn_calls.append(idx)
        noise_dt = _dt_randn_sequence[idx]
        if scale != 1.0:
            noise_dt = scalar_tensor_op(scale, noise_dt, ElementwiseOp.PROD)
        return noise_dt

    monkeypatch.setattr(distributed_diffusion_module, "center_random_augmentation", _centering_only_augmentation)
    monkeypatch.setattr(distributed_diffusion_module, "create_distributed_randn", _fixed_create_distributed_randn)

    module.predict_args = predict_args
    with torch.no_grad():
        pred_dict = module.predict_step(feats_dt, batch_idx=0)

    assert pred_dict["exception"] is False, "predict_step raised an exception"

    expected_confidence_keys = {"pde", "plddt", "complex_plddt", "complex_iplddt", "complex_pde", "complex_ipde"}
    missing = expected_confidence_keys - set(pred_dict.keys())
    assert not missing, f"Missing confidence keys in predict_step output: {missing}"

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env",
    [
        ((1, (1, 1)), True, "cuda", "ENV"),
    ],
    indirect=("setup_env",),
    ids=lambda val: f"{val[2]}-dp:{val[0][0]}-cp{'x'.join(map(str, val[0][1]))}",
)
def test_boltz2_predict_step_confidence(setup_env):
    """V14b: predict_step with confidence_prediction=True produces no DTensor leaks.

    Runs on 1 GPU (dp=1, cp=1x1) and verifies that the predict_step output
    dict contains only plain tensors, catching any DTensor-to-writer leaks
    that would crash the BoltzWriter callback. Even with trivial CP sharding,
    the model produces DTensor outputs; the _assert_no_dtensors_in_output guard
    in predict_step catches any unconverted DTensors at runtime.
    """
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("CUDA not available")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"Need {world_size} GPUs, have {torch.cuda.device_count()}")

    dtype = torch.float64
    size_batch = 1
    size_cp = grid_group_sizes["cp"][0]

    num_sampling_steps = 2
    diffusion_samples = 1

    seed = 42
    seed_by_rank(0, seed=seed)

    boltz2_model_params = create_boltz2_model_init_params(use_large_model=False)
    boltz2_model_params["diffusion_process_args"]["alignment_reverse_diff"] = False
    boltz2_model_params["num_bins"] = 64
    boltz2_model_params["confidence_prediction"] = True
    boltz2_model_params["confidence_model_args"] = {
        "pairformer_args": boltz2_model_params["pairformer_args"],
        "confidence_args": {},
    }

    predict_args = {
        "recycling_steps": 0,
        "sampling_steps": num_sampling_steps,
        "diffusion_samples": diffusion_samples,
        "max_parallel_samples": diffusion_samples,
        "write_confidence_summary": True,
        "write_full_pae": True,
    }
    boltz2_model_params["predict_args"] = predict_args

    n_atoms_per_token_min = 8
    n_atoms_per_token_max = 20
    n_tokens = 30 * size_cp
    W = boltz2_model_params["atoms_per_window_queries"]
    n_atoms_raw = n_tokens * n_atoms_per_token_max
    n_atoms = ((n_atoms_raw + W - 1) // W) * W
    n_msa = size_cp * 2

    assert n_atoms % size_cp == 0
    assert n_atoms % W == 0

    input_feats_global_fp64 = random_features(
        size_batch=size_batch,
        n_tokens=n_tokens,
        n_atoms=n_atoms,
        n_msa=n_msa,
        atom_counts_per_token_range=(n_atoms_per_token_min, n_atoms_per_token_max),
        device=device_type,
        float_value_range=(-0.01, 0.01),
        selected_keys=_BOLTZ2_SELECTED_KEYS,
        num_disto_bins=boltz2_model_params["num_bins"],
    )
    input_feats_global_fp64["msa"] = torch.randint(
        0, const.num_tokens, (size_batch, n_msa, n_tokens), dtype=torch.int64, device=device_type
    )

    reference_module = SerialBoltz2(**boltz2_model_params)
    init_module_params_glorot(reference_module, gain=0.05)
    reference_module.apply(SetModuleInfValues())
    reference_module.structure_module.coordinate_augmentation = False
    module_state_dict_fp64 = reference_module.state_dict()

    _B_M = size_batch * diffusion_samples
    init_noise = torch.empty((_B_M, n_atoms, 3), device=device_type, dtype=dtype)
    step_noise_list = [
        torch.empty((_B_M, n_atoms, 3), device=device_type, dtype=dtype) for _ in range(num_sampling_steps)
    ]
    init_tensors_uniform([init_noise, *step_noise_list], low=-0.01, high=0.01)

    input_feats_host = {k: v.detach().to(device="cpu", copy=True) for k, v in input_feats_global_fp64.items()}

    spawn_multiprocessing(
        _worker_predict_step_confidence,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        dtype,
        boltz2_model_params,
        module_state_dict_fp64,
        predict_args,
        diffusion_samples,
        num_sampling_steps,
        input_feats_host,
        init_noise.cpu(),
        [n.cpu() for n in step_noise_list],
        env_per_rank,
    )


# ====================================================================== #
#  V13: Forward/backward parity (serial vs distributed)                   #
# ====================================================================== #


def _worker_forward_backward_parity(
    rank,
    grid_group_sizes,
    device_type,
    backend,
    dtype,
    boltz2_model_params,
    module_state_dict,
    n_recycles,
    multiplicity_diffusion_train,
    input_feats_global_fp64_host,
    sigmas_expected_global_fp64_host,
    noise_expected_global_fp64_host,
    output_pdistogram_expected_global_host,
    output_denoised_atom_coords_expected_global_host,
    output_pbfactor_expected_global_host,
    output_s_expected_global_host,
    output_z_expected_global_host,
    output_aligned_true_coords_expected_global_host,
    d_output_pdistogram_expected_global_host,
    d_output_denoised_atom_coords_expected_global_host,
    d_output_pbfactor_expected_global_host,
    expected_param_grads_global_host_dict,
    env_per_rank=None,
    use_random_features=True,
    training_data_dir=None,
    canonical_mols_dir=None,
    base_seed=42,
):
    """V13 multi-rank worker: verify distributed forward/backward matches serial."""
    monkeypatch = pytest.MonkeyPatch()
    if env_per_rank is not None:
        for var_name, value in env_per_rank.items():
            monkeypatch.setenv(var_name, f"{rank}" if value == "<INPUT_RANK>" else value)

    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    manager = DistributedManager()

    reference_module = SerialBoltz2(**boltz2_model_params)
    reference_module = reference_module.to(dtype=dtype)
    reference_module.load_state_dict(module_state_dict)
    reference_module.structure_module.coordinate_augmentation = False
    reference_module.apply(SetModuleInfValues())
    reference_module = reference_module.to(device=manager.device)
    module = Boltz2Distributed(reference_module, manager)
    module.train()

    _io_model_io_keys = {"noise", "denoised_atom_coords", "d_denoised_atom_coords", "aligned_true_atom_coords"}

    if use_random_features:
        host_tensor_keys = {k for k, v in input_feats_global_fp64_host.items() if isinstance(v, torch.Tensor)}
        _placements = get_feature_placements(
            token_keys=host_tensor_keys,
            msa_keys=host_tensor_keys,
            atom_keys={
                "ref_pos",
                "atom_resolved_mask",
                "ref_element",
                "ref_charge",
                "ref_atom_name_chars",
                "ref_space_uid",
                "coords",
                "atom_pad_mask",
                "atom_to_token",
                "pair_mask",
                "atom_counts_per_token",
                "token_to_rep_atom",
                "bfactor",
                "plddt",
            },
            model_io_keys=_io_model_io_keys,
            model_io_fp32_keys=set(),
        )

        # Token + MSA features: broadcast from rank 0
        if manager.group_rank["world"] == 0:
            input_feats_token_msa_global = {
                k: v.to(device=manager.device, dtype=dtype if v.dtype.is_floating_point else v.dtype)
                for k, v in input_feats_global_fp64_host.items()
                if k in _placements["token_features"] or k in _placements["msa_features"]
            }
        else:
            input_feats_token_msa_global = None

        feats_token_msa = distribute_features(
            input_feats_token_msa_global,
            _placements["token_features"] | _placements["msa_features"],
            manager.group["world"],
            manager.group_ranks["world"][0],
            manager.device_mesh_subgroups,
        )

        # Atom features: scatter across CP mesh
        size_batch = input_feats_global_fp64_host["atom_pad_mask"].shape[0]
        inputs_atom = {
            k: v.to(dtype=dtype if v.dtype.is_floating_point else v.dtype)
            for k, v in input_feats_global_fp64_host.items()
            if k in _placements["cp_atom_features"]
        }

        noise_unflat = noise_expected_global_fp64_host.unflatten(0, (size_batch, multiplicity_diffusion_train))
        denoised_unflat = output_denoised_atom_coords_expected_global_host.unflatten(
            0, (size_batch, multiplicity_diffusion_train)
        )
        d_denoised_unflat = d_output_denoised_atom_coords_expected_global_host.unflatten(
            0, (size_batch, multiplicity_diffusion_train)
        )
        aligned_coords_unflat = output_aligned_true_coords_expected_global_host.unflatten(
            0, (size_batch, multiplicity_diffusion_train)
        )
        for i_mul in range(multiplicity_diffusion_train):
            inputs_atom[f"noise_{i_mul}"] = noise_unflat[:, i_mul].to(dtype=dtype)
            inputs_atom[f"denoised_atom_coords_{i_mul}"] = denoised_unflat[:, i_mul].to(dtype=dtype)
            inputs_atom[f"d_denoised_atom_coords_{i_mul}"] = d_denoised_unflat[:, i_mul].to(dtype=dtype)
            inputs_atom[f"aligned_true_atom_coords_{i_mul}"] = aligned_coords_unflat[:, i_mul].to(dtype=dtype)

        placements_cp_model_io_mul = {
            f"{k}_{i_mul}": v
            for k, v in _placements["cp_model_io"].items()
            for i_mul in range(multiplicity_diffusion_train)
        }
        placements_cp = _placements["cp_atom_features"] | placements_cp_model_io_mul

        placements_model_io_mul = {
            f"{k}_{i_mul}": v
            for k, v in _placements["model_io"].items()
            for i_mul in range(multiplicity_diffusion_train)
        }
        placements_dp_cp = _placements["atom_features"] | placements_model_io_mul

        multiplicities = {
            "noise": multiplicity_diffusion_train,
            "denoised_atom_coords": multiplicity_diffusion_train,
            "d_denoised_atom_coords": multiplicity_diffusion_train,
            "aligned_true_atom_coords": multiplicity_diffusion_train,
        }

        feats_atom = distribute_atom_features(
            inputs=inputs_atom,
            placements_cp=placements_cp,
            placements_dp_cp=placements_dp_cp,
            device_mesh=manager.device_mesh_subgroups,
            cp_group=manager.group["cp"],
            multiplicities=multiplicities,
        )

        feats = {**feats_token_msa, **feats_atom}

        noise_dt = feats_atom.pop("noise")
        expected_denoised_dt = feats_atom.pop("denoised_atom_coords")
        d_denoised_dt_expected = feats_atom.pop("d_denoised_atom_coords")
        expected_aligned_coords_dt = feats_atom.pop("aligned_true_atom_coords")

    else:
        # --- Dataloader path: load features from distributed data module ---
        DistributedManager.create_group(
            "world_cpu", manager.group_ranks["world"], backend="gloo", use_local_synchronization=True
        )
        DistributedManager.create_group(
            "cp_cpu", manager.group_ranks["cp"], backend="gloo", use_local_synchronization=True
        )
        cp_device_mesh = map_subgroup_mesh_to_cpu(manager)

        _mp_data = pytest.MonkeyPatch()
        cfg = setup_mock_training_datamodule_config(training_data_dir)
        cfg.batch_size = 1
        cfg.samples_per_epoch = grid_group_sizes["dp"]
        cfg.moldir = str(canonical_mols_dir)
        cfg.return_train_symmetries = False
        for ds_cfg in cfg.datasets:
            ds_cfg.filters = None
        seed_by_rank(0, seed=base_seed)
        dm = BoltzTrainingDataModuleDTensor(cfg, manager.device_mesh_subgroups, cp_device_mesh)
        _deterministic_getitem_monkeypatch(_mp_data, dm._serial_module._train_set, base_seed=base_seed)
        dl = dm.train_dataloader()
        batch_cpu = next(iter(dl))
        batch_gpu = dm.transfer_batch_to_device(batch_cpu, manager.device, 0)
        _mp_data.undo()

        feats = {}
        for k, v in batch_gpu.items():
            if isinstance(v, (DTensor, torch.Tensor)) and v.dtype.is_floating_point:
                feats[k] = v.to(dtype=dtype)
            elif isinstance(v, list):
                feats[k] = [
                    item.to(dtype=dtype) if isinstance(item, torch.Tensor) and item.dtype.is_floating_point else item
                    for item in v
                ]
            else:
                feats[k] = v

        # Distribute noise + model I/O reference tensors via intersperse padding
        _placements = get_feature_placements(
            atom_keys=set(),
            model_io_keys=_io_model_io_keys,
            model_io_fp32_keys=set(),
        )
        size_batch = input_feats_global_fp64_host["atom_pad_mask"].shape[0]
        inputs_io = {"atom_counts_per_token": input_feats_global_fp64_host["atom_counts_per_token"].clone()}

        noise_unflat = noise_expected_global_fp64_host.unflatten(0, (size_batch, multiplicity_diffusion_train))
        denoised_unflat = output_denoised_atom_coords_expected_global_host.unflatten(
            0, (size_batch, multiplicity_diffusion_train)
        )
        d_denoised_unflat = d_output_denoised_atom_coords_expected_global_host.unflatten(
            0, (size_batch, multiplicity_diffusion_train)
        )
        aligned_coords_unflat = output_aligned_true_coords_expected_global_host.unflatten(
            0, (size_batch, multiplicity_diffusion_train)
        )
        for i_mul in range(multiplicity_diffusion_train):
            inputs_io[f"noise_{i_mul}"] = noise_unflat[:, i_mul].to(dtype=dtype)
            inputs_io[f"denoised_atom_coords_{i_mul}"] = denoised_unflat[:, i_mul].to(dtype=dtype)
            inputs_io[f"d_denoised_atom_coords_{i_mul}"] = d_denoised_unflat[:, i_mul].to(dtype=dtype)
            inputs_io[f"aligned_true_atom_coords_{i_mul}"] = aligned_coords_unflat[:, i_mul].to(dtype=dtype)

        placements_cp_model_io_mul = {
            f"{k}_{i_mul}": v
            for k, v in _placements["cp_model_io"].items()
            for i_mul in range(multiplicity_diffusion_train)
        }
        placements_cp = _placements["cp_atom_features"] | placements_cp_model_io_mul

        placements_model_io_mul = {
            f"{k}_{i_mul}": v
            for k, v in _placements["model_io"].items()
            for i_mul in range(multiplicity_diffusion_train)
        }
        placements_dp_cp = placements_model_io_mul

        io_feats = distribute_atom_features(
            inputs=inputs_io,
            placements_cp=placements_cp,
            placements_dp_cp=placements_dp_cp,
            device_mesh=manager.device_mesh_subgroups,
            cp_group=manager.group["cp"],
            multiplicities={k: multiplicity_diffusion_train for k in _io_model_io_keys},
        )

        # distribute_atom_features applies intersperse padding per DP rank
        # independently, but CollateDTensor additionally homogenizes local
        # shard shapes across DP ranks via all-reduce MAX.  When samples have
        # different atom counts (7ylz vs 8b2e), the dataloader features are
        # homogenized but the model-I/O DTensors from distribute_atom_features
        # are not.  Pad each DTensor's atom dim to match the batch.
        target_atoms_global = feats["atom_pad_mask"].shape[-1]
        for k in list(io_feats.keys()):
            if io_feats[k].shape[1] < target_atoms_global:
                io_feats[k] = pad_to_length(io_feats[k], dim=1, length=target_atoms_global)

        noise_dt = io_feats.pop("noise")
        expected_denoised_dt = io_feats.pop("denoised_atom_coords")
        d_denoised_dt_expected = io_feats.pop("d_denoised_atom_coords")
        expected_aligned_coords_dt = io_feats.pop("aligned_true_atom_coords")

    # Monkeypatch deterministic noise for distributed forward
    sigmas_device = sigmas_expected_global_fp64_host.to(device=manager.device, dtype=dtype)
    sigmas_dt = distribute_tensor(sigmas_device, manager.device_mesh_subgroups, (Shard(0), Replicate(), Replicate()))

    monkeypatch.setattr(module.structure_module, "noise_distribution", lambda bs, dtype=None: sigmas_dt)
    monkeypatch.setattr(distributed_diffusion_module, "create_distributed_randn", lambda *a, **kw: noise_dt)

    # Distributed forward
    output_dict = module(
        feats,
        recycling_steps=n_recycles,
        multiplicity_diffusion_train=multiplicity_diffusion_train,
    )

    assert "pdistogram" in output_dict
    assert "denoised_atom_coords" in output_dict
    assert "pbfactor" in output_dict
    assert "s" in output_dict
    assert "z" in output_dict
    assert "sigmas" in output_dict
    assert "aligned_true_atom_coords" in output_dict

    token_pad_mask_global = feats["token_pad_mask"].full_tensor()
    token_pair_pad_mask_global = feats["token_pair_pad_mask"].full_tensor()
    atom_pad_mask_global = feats["atom_pad_mask"].full_tensor()
    atom_pad_mask_mul_global = atom_pad_mask_global[:, :, None].repeat_interleave(multiplicity_diffusion_train, 0)

    s_full = output_dict["s"].full_tensor() * token_pad_mask_global[:, :, None]
    expected_s = output_s_expected_global_host.to(device=manager.device, dtype=dtype)
    torch.testing.assert_close(s_full, expected_s)

    z_full = output_dict["z"].full_tensor() * token_pair_pad_mask_global[:, :, :, None]
    expected_z = output_z_expected_global_host.to(device=manager.device, dtype=dtype)
    torch.testing.assert_close(z_full, expected_z)

    pdistogram_full = output_dict["pdistogram"].full_tensor() * token_pair_pad_mask_global[:, :, :, None, None]
    expected_pdistogram = output_pdistogram_expected_global_host.to(device=manager.device, dtype=dtype)
    torch.testing.assert_close(pdistogram_full, expected_pdistogram)

    denoised_full = output_dict["denoised_atom_coords"].full_tensor() * atom_pad_mask_mul_global
    expected_denoised_full = expected_denoised_dt.full_tensor() * atom_pad_mask_mul_global
    torch.testing.assert_close(denoised_full, expected_denoised_full)

    pbfactor_full = output_dict["pbfactor"].full_tensor() * token_pad_mask_global[:, :, None]
    expected_pbfactor = output_pbfactor_expected_global_host.to(device=manager.device, dtype=dtype)
    torch.testing.assert_close(pbfactor_full, expected_pbfactor)

    sigmas_full = output_dict["sigmas"].full_tensor()
    expected_sigmas = sigmas_expected_global_fp64_host.to(device=manager.device, dtype=dtype)
    torch.testing.assert_close(sigmas_full, expected_sigmas)

    aligned_coords_full = output_dict["aligned_true_atom_coords"].full_tensor() * atom_pad_mask_mul_global
    expected_aligned_coords_full = expected_aligned_coords_dt.full_tensor() * atom_pad_mask_mul_global
    torch.testing.assert_close(aligned_coords_full, expected_aligned_coords_full)

    # Backward pass
    d_pdistogram = d_output_pdistogram_expected_global_host.to(device=manager.device, dtype=dtype)
    d_pdistogram_dt = distribute_tensor(
        d_pdistogram, manager.device_mesh_subgroups, output_dict["pdistogram"].placements
    )

    d_pbfactor = d_output_pbfactor_expected_global_host.to(device=manager.device, dtype=dtype)
    d_pbfactor_dt = distribute_tensor(d_pbfactor, manager.device_mesh_subgroups, output_dict["pbfactor"].placements)

    torch.autograd.backward(
        [output_dict["pdistogram"], output_dict["denoised_atom_coords"], output_dict["pbfactor"]],
        [d_pdistogram_dt, d_denoised_dt_expected, d_pbfactor_dt],
    )

    num_grads_checked = 0
    num_nonzero_grads = 0
    for name, param in module.named_parameters():
        canonical_name = name.replace("._serial.", ".")
        if canonical_name in expected_param_grads_global_host_dict:
            expected_grad = expected_param_grads_global_host_dict[canonical_name].to(device=manager.device, dtype=dtype)
            if param.grad is None:
                raise AssertionError(f"Missing gradient for {canonical_name}")
            actual_grad = param.grad.full_tensor() if isinstance(param.grad, DTensor) else param.grad
            num_grads_checked += 1
            if expected_grad.abs().max().item() > 0:
                num_nonzero_grads += 1
            torch.testing.assert_close(
                actual_grad,
                expected_grad,
                msg=lambda msg, cn=canonical_name: f"Gradient mismatch for {cn}: {msg}",
            )

    assert num_grads_checked > 0, "No gradients compared — test is vacuous"
    assert num_nonzero_grads > 0, "All compared gradients are zero — test is vacuous"

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.parametrize(
    ("setup_env", "use_random_features"),
    [
        (((2, (2, 2)), True, "cuda", "ENV"), True),
        (((2, (2, 2)), True, "cuda", "ENV"), False),
    ],
    indirect=["setup_env"],
    ids=["cuda-dp2-cp2x2-random", "cuda-dp2-cp2x2-dataloader"],
)
def test_boltz2_forward_backward_parity(
    setup_env,
    use_random_features,
    test_cp_training_base_data_dir_boltz2,
    canonical_mols_dir,
    tmp_path,
):
    """V13: Forward/backward parity between distributed Boltz2 and serial Boltz2.

    Tests that the distributed Boltz2 wrapper produces numerically identical
    forward outputs and backward gradients compared to the serial implementation,
    with training=True and structure_prediction_training=True.

    This test uses custom upstream gradients (not loss-derived) to isolate
    the forward/backward pipeline from the loss computation.  The model
    configuration, initialization, and features are intentionally aligned
    with test_boltz2_training_step_parity so that a pass here proves the
    forward/backward path is correct under the same numerical regime.

    Parametrized by use_random_features:
    - True: synthetic random features (fast, small dimensions)
    - False: real 7ylz + 8b2e training data via TrainingDataModule
    """
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("CUDA not available")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"Need {world_size} GPUs, have {torch.cuda.device_count()}")

    dtype = torch.float64
    B = grid_group_sizes["dp"]
    size_cp = grid_group_sizes["cp"][0]

    min_val_init = -0.1
    max_val_init = 0.1
    scale_glorot = 0.1

    seed = 42
    seed_by_rank(0, seed=seed)

    boltz2_model_params = create_boltz2_model_init_params(use_large_model=False)
    n_recycles = 0
    multiplicity_diffusion_train = 2

    boltz2_model_params["training_args"].recycling_steps = n_recycles
    boltz2_model_params["training_args"].diffusion_multiplicity = multiplicity_diffusion_train
    boltz2_model_params["no_random_recycling_training"] = True

    if use_random_features:
        n_tokens = 30 * size_cp
        W = boltz2_model_params["atoms_per_window_queries"]
        n_atoms_raw = n_tokens * 20
        n_atoms = ((n_atoms_raw + W - 1) // W) * W
        n_msa = max(size_cp * 2, 2)

        input_feats_global_fp64 = random_features(
            size_batch=B,
            n_tokens=n_tokens,
            n_atoms=n_atoms,
            n_msa=n_msa,
            atom_counts_per_token_range=(8, 20),
            device=device_type,
            float_value_range=(min_val_init, max_val_init),
            selected_keys=_BOLTZ2_SELECTED_KEYS,
            num_disto_bins=boltz2_model_params["num_bins"],
        )
        input_feats_global_fp64["msa"] = torch.randint(
            0, const.num_tokens, (B, n_msa, n_tokens), dtype=torch.int64, device=device_type
        )
        training_data_dir = None
    else:
        boltz2_model_params["num_bins"] = 64
        training_data_dir = _setup_training_data_7z64_8b2e(
            tmp_path / "training_data", test_cp_training_base_data_dir_boltz2
        )
        _mp_data = pytest.MonkeyPatch()
        cfg = setup_mock_training_datamodule_config(training_data_dir)
        cfg.batch_size = B
        cfg.samples_per_epoch = B
        cfg.moldir = str(canonical_mols_dir)
        cfg.return_train_symmetries = False
        for ds_cfg in cfg.datasets:
            ds_cfg.filters = None
        seed_by_rank(0, seed=seed)
        dm = Boltz2TrainingDataModuleSerial(cfg=cfg)
        _deterministic_getitem_monkeypatch(_mp_data, dm._train_set, base_seed=seed)
        dl = dm.train_dataloader()
        raw_batch = next(iter(dl))
        input_feats_global_fp64 = {}
        for k, v in raw_batch.items():
            if isinstance(v, torch.Tensor):
                input_feats_global_fp64[k] = v.to(
                    device=device_type, dtype=dtype if v.dtype.is_floating_point else v.dtype
                )
            elif isinstance(v, list):
                input_feats_global_fp64[k] = [
                    item.to(device=device_type, dtype=dtype if item.dtype.is_floating_point else item.dtype)
                    if isinstance(item, torch.Tensor)
                    else item
                    for item in v
                ]
            else:
                input_feats_global_fp64[k] = v
        n_atoms = input_feats_global_fp64["atom_pad_mask"].shape[-1]
        _mp_data.undo()

        # Real data doesn't include token_pair_pad_mask; compute from token_pad_mask
        if "token_pair_pad_mask" not in input_feats_global_fp64:
            tpm = input_feats_global_fp64["token_pad_mask"]
            input_feats_global_fp64["token_pair_pad_mask"] = tpm[:, :, None] * tpm[:, None, :]

    # Create serial reference module
    reference_module = SerialBoltz2(**boltz2_model_params)
    init_module_params_glorot(reference_module, gain=scale_glorot)
    reference_module.apply(SetModuleInfValues())
    reference_module.structure_module.coordinate_augmentation = False
    module_state_dict_fp64 = {k: v.detach().clone().cpu() for k, v in reference_module.state_dict().items()}
    reference_module = reference_module.to(dtype=torch.float64, device=device_type).train()

    # Pre-generate deterministic sigmas and non-zero noise
    sigmas_expected_global_fp64 = reference_module.structure_module.noise_distribution(
        B * multiplicity_diffusion_train
    ).to(device=device_type, dtype=torch.float64)
    noise_expected_global_fp64 = torch.empty(
        B * multiplicity_diffusion_train, n_atoms, 3, device=device_type, dtype=torch.float64
    )
    init_tensors_uniform([noise_expected_global_fp64], low=min_val_init, high=max_val_init)

    # Monkeypatch serial noise for determinism
    _monkeypatch = pytest.MonkeyPatch()
    _monkeypatch.setattr(
        reference_module.structure_module, "noise_distribution", lambda bs, dtype=None: sigmas_expected_global_fp64
    )
    _monkeypatch.setattr(serial_diffusion_v2_module.torch, "randn_like", lambda t: noise_expected_global_fp64.to(t))

    original_feat_keys = set(input_feats_global_fp64.keys())
    coords_backup = input_feats_global_fp64["coords"].detach().clone()

    # Serial forward
    output_dict_serial = reference_module(
        input_feats_global_fp64,
        recycling_steps=n_recycles,
        multiplicity_diffusion_train=multiplicity_diffusion_train,
    )

    output_pdistogram = output_dict_serial["pdistogram"]
    output_denoised = output_dict_serial["denoised_atom_coords"]
    output_pbfactor = output_dict_serial["pbfactor"]
    output_s = output_dict_serial["s"]
    output_z = output_dict_serial["z"]
    output_aligned_true_coords = output_dict_serial["aligned_true_atom_coords"]

    # Create upstream gradients
    d_output_pdistogram = torch.empty_like(output_pdistogram)
    d_output_denoised = torch.empty_like(output_denoised)
    d_output_pbfactor = torch.empty_like(output_pbfactor)
    init_tensors_uniform(
        [d_output_pdistogram, d_output_denoised, d_output_pbfactor], low=min_val_init, high=max_val_init
    )

    # Mask upstream gradients
    atom_pad_mask = input_feats_global_fp64["atom_pad_mask"]
    atom_pad_mask_mul = atom_pad_mask[:, :, None].repeat_interleave(multiplicity_diffusion_train, 0)
    d_output_denoised = d_output_denoised * atom_pad_mask_mul

    token_pair_pad_mask = input_feats_global_fp64["token_pair_pad_mask"]
    d_output_pdistogram = d_output_pdistogram * token_pair_pad_mask[:, :, :, None, None]

    token_pad_mask = input_feats_global_fp64["token_pad_mask"]
    d_output_pbfactor = d_output_pbfactor * token_pad_mask[:, :, None]

    # Serial backward
    torch.autograd.backward(
        [output_pdistogram, output_denoised, output_pbfactor],
        [d_output_pdistogram, d_output_denoised, d_output_pbfactor],
    )

    grad_params_expected = {
        name: param.grad.detach().to(dtype=dtype, device="cpu", copy=True)
        for name, param in reference_module.named_parameters()
        if param.grad is not None
    }

    # Restore features — serial forward mutates coords in-place and may add keys
    input_feats_global_fp64["coords"] = coords_backup
    for key in list(input_feats_global_fp64.keys()):
        if key not in original_feat_keys:
            del input_feats_global_fp64[key]

    output_pdistogram_host = (
        (output_pdistogram * token_pair_pad_mask[:, :, :, None, None]).detach().to(device="cpu", copy=True)
    )
    output_denoised_host = (output_denoised * atom_pad_mask_mul).detach().to(device="cpu", copy=True)
    output_pbfactor_host = (output_pbfactor * token_pad_mask[:, :, None]).detach().to(device="cpu", copy=True)
    output_s_host = (output_s * token_pad_mask[:, :, None]).detach().to(device="cpu", copy=True)
    output_z_host = (output_z * token_pair_pad_mask[:, :, :, None]).detach().to(device="cpu", copy=True)
    output_aligned_true_coords_host = (
        (output_aligned_true_coords * atom_pad_mask_mul).detach().to(device="cpu", copy=True)
    )

    sigmas_host = sigmas_expected_global_fp64.detach().to(device="cpu", copy=True)
    noise_host = (noise_expected_global_fp64 * atom_pad_mask_mul).detach().to(device="cpu", copy=True)
    d_output_pdistogram_host = d_output_pdistogram.detach().to(device="cpu", copy=True)
    d_output_denoised_host = d_output_denoised.detach().to(device="cpu", copy=True)
    d_output_pbfactor_host = d_output_pbfactor.detach().to(device="cpu", copy=True)

    input_feats_host = {}
    for k, v in input_feats_global_fp64.items():
        if isinstance(v, torch.Tensor):
            input_feats_host[k] = v.detach().to(device="cpu", copy=True)
        elif isinstance(v, list):
            input_feats_host[k] = [
                item.detach().to(device="cpu", copy=True) if isinstance(item, torch.Tensor) else item for item in v
            ]
        else:
            input_feats_host[k] = v

    _monkeypatch.undo()

    spawn_multiprocessing(
        _worker_forward_backward_parity,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        dtype,
        boltz2_model_params,
        module_state_dict_fp64,
        n_recycles,
        multiplicity_diffusion_train,
        input_feats_host,
        sigmas_host,
        noise_host,
        output_pdistogram_host,
        output_denoised_host,
        output_pbfactor_host,
        output_s_host,
        output_z_host,
        output_aligned_true_coords_host,
        d_output_pdistogram_host,
        d_output_denoised_host,
        d_output_pbfactor_host,
        grad_params_expected,
        env_per_rank,
        use_random_features,
        training_data_dir,
        canonical_mols_dir,
        seed,
    )


# ====================================================================== #
#  V16: Serial vs distributed training_step numerical parity              #
# ====================================================================== #


# V15 (training_step smoke test) was removed — V16 (training_step parity)
# provides strictly stronger guarantees, and the recycling path is covered
# by test_boltz2_train_entrypoint / test_boltz2_stop_and_go in
# test_dtensor_boltz2_train.py.  This matches the Boltz-1 test structure
# which has only forward parity + training_step parity (no smoke test).


def _worker_training_step_parity(
    rank,
    grid_group_sizes,
    device_type,
    backend,
    env_per_rank,
    module_state_dict,
    boltz2_model_params,
    input_feats_global_fp64_host,
    sigmas_expected_host,
    noise_expected_host,
    serial_loss_host,
    serial_log_metrics_host,
    serial_grad_dict_host,
    serial_post_opt_dict_host,
    optimizer_step,
    use_random_features=True,
    training_data_dir=None,
    canonical_mols_dir=None,
    base_seed=42,
):
    """V16 multi-rank worker: verify distributed training_step matches serial.

    Compares:
    1. Loss value (serial vs distributed)
    1b. Logged metric values via CSVLogger (serial vs distributed),
        including component-wise grad_norms (compared after backward
        + on_after_backward so grad_norm metrics are available)
    2. Per-parameter gradients after backward + on_after_backward
    3. Post-optimizer parameter values (if optimizer_step=True)
    """
    monkeypatch = pytest.MonkeyPatch()
    if env_per_rank is not None:
        for var_name, value in env_per_rank.items():
            monkeypatch.setenv(var_name, f"{rank}" if value == "<INPUT_RANK>" else value)

    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    manager = DistributedManager()
    dtype = torch.float64

    # Build distributed model from serial state dict
    reference_module = SerialBoltz2(**boltz2_model_params)
    reference_module = reference_module.to(dtype=dtype)
    reference_module.load_state_dict(module_state_dict)
    reference_module.structure_module.coordinate_augmentation = False
    reference_module.apply(SetModuleInfValues())
    reference_module = reference_module.to(device=manager.device)
    module = Boltz2Distributed(reference_module, manager)
    module = module.to(device=manager.device)
    module.train()

    # Capture initial parameter values for non-vacuous optimizer check
    initial_params = {
        name: (p.full_tensor().detach().clone() if isinstance(p, DTensor) else p.detach().clone())
        for name, p in module.named_parameters()
        if p.requires_grad
    }

    # Inject CSVLogger to exercise the logging code path (instead of no-op)
    worker_csv_logger = CSVLogger(save_dir=tempfile.mkdtemp(), name=f"distributed_rank{rank}")
    dist_log = _LogCapture(worker_csv_logger)
    monkeypatch.setattr(module, "log", dist_log)
    monkeypatch.setattr(module, "training_log", lambda *a, **kw: None)

    multiplicity = boltz2_model_params["training_args"].diffusion_multiplicity

    if use_random_features:
        # Distribute features (same pattern as V13)
        host_tensor_keys = {k for k, v in input_feats_global_fp64_host.items() if isinstance(v, torch.Tensor)}
        _placements = get_feature_placements(
            token_keys=host_tensor_keys,
            msa_keys=host_tensor_keys,
            atom_keys={
                "ref_pos",
                "atom_resolved_mask",
                "ref_element",
                "ref_charge",
                "ref_atom_name_chars",
                "ref_space_uid",
                "coords",
                "atom_pad_mask",
                "atom_to_token",
                "pair_mask",
                "atom_counts_per_token",
                "token_to_rep_atom",
                "bfactor",
                "plddt",
            },
            model_io_keys={"noise"},
            model_io_fp32_keys=set(),
        )

        # Token + MSA features: broadcast from rank 0
        if manager.group_rank["world"] == 0:
            input_feats_token_msa_global = {
                k: v.to(device=manager.device, dtype=dtype if v.dtype.is_floating_point else v.dtype)
                for k, v in input_feats_global_fp64_host.items()
                if k in _placements["token_features"] or k in _placements["msa_features"]
            }
        else:
            input_feats_token_msa_global = None

        feats_token_msa = distribute_features(
            input_feats_token_msa_global,
            _placements["token_features"] | _placements["msa_features"],
            manager.group["world"],
            manager.group_ranks["world"][0],
            manager.device_mesh_subgroups,
        )

        # Atom features + noise: scatter across CP mesh with multiplicity
        size_batch = input_feats_global_fp64_host["atom_pad_mask"].shape[0]
        inputs_atom = {
            k: v.to(dtype=dtype if v.dtype.is_floating_point else v.dtype)
            for k, v in input_feats_global_fp64_host.items()
            if k in _placements["cp_atom_features"]
        }

        noise_unflat = noise_expected_host.unflatten(0, (size_batch, multiplicity))
        for i_mul in range(multiplicity):
            inputs_atom[f"noise_{i_mul}"] = noise_unflat[:, i_mul].to(dtype=dtype)

        placements_cp_model_io_mul = {
            f"{k}_{i_mul}": v for k, v in _placements["cp_model_io"].items() for i_mul in range(multiplicity)
        }
        placements_cp = _placements["cp_atom_features"] | placements_cp_model_io_mul

        placements_model_io_mul = {
            f"{k}_{i_mul}": v for k, v in _placements["model_io"].items() for i_mul in range(multiplicity)
        }
        placements_dp_cp = _placements["atom_features"] | placements_model_io_mul

        feats_atom = distribute_atom_features(
            inputs=inputs_atom,
            placements_cp=placements_cp,
            placements_dp_cp=placements_dp_cp,
            device_mesh=manager.device_mesh_subgroups,
            cp_group=manager.group["cp"],
            multiplicities={"noise": multiplicity},
        )

        noise_dt = feats_atom.pop("noise")
        batch = {**feats_token_msa, **feats_atom}

    else:
        # --- Dataloader path: load features from distributed data module ---
        DistributedManager.create_group(
            "world_cpu", manager.group_ranks["world"], backend="gloo", use_local_synchronization=True
        )
        DistributedManager.create_group(
            "cp_cpu", manager.group_ranks["cp"], backend="gloo", use_local_synchronization=True
        )
        cp_device_mesh = map_subgroup_mesh_to_cpu(manager)

        _mp_data = pytest.MonkeyPatch()
        cfg = setup_mock_training_datamodule_config(training_data_dir)
        cfg.batch_size = 1
        cfg.samples_per_epoch = grid_group_sizes["dp"]
        cfg.moldir = str(canonical_mols_dir)
        cfg.return_train_symmetries = False
        for ds_cfg in cfg.datasets:
            ds_cfg.filters = None
        seed_by_rank(0, seed=base_seed)
        dm = BoltzTrainingDataModuleDTensor(cfg, manager.device_mesh_subgroups, cp_device_mesh)
        _deterministic_getitem_monkeypatch(_mp_data, dm._serial_module._train_set, base_seed=base_seed)
        dl = dm.train_dataloader()
        batch_cpu = next(iter(dl))
        batch_gpu = dm.transfer_batch_to_device(batch_cpu, manager.device, 0)
        _mp_data.undo()

        feats = {}
        for k, v in batch_gpu.items():
            if isinstance(v, (DTensor, torch.Tensor)) and v.dtype.is_floating_point:
                feats[k] = v.to(dtype=dtype)
            elif isinstance(v, list):
                feats[k] = [
                    item.to(dtype=dtype) if isinstance(item, torch.Tensor) and item.dtype.is_floating_point else item
                    for item in v
                ]
            else:
                feats[k] = v

        # Distribute noise via intersperse padding + homogenization
        _io_keys = {"noise"}
        _placements = get_feature_placements(
            atom_keys=set(),
            model_io_keys=_io_keys,
            model_io_fp32_keys=set(),
        )
        size_batch = input_feats_global_fp64_host["atom_pad_mask"].shape[0]
        inputs_io = {"atom_counts_per_token": input_feats_global_fp64_host["atom_counts_per_token"].clone()}

        noise_unflat = noise_expected_host.unflatten(0, (size_batch, multiplicity))
        for i_mul in range(multiplicity):
            inputs_io[f"noise_{i_mul}"] = noise_unflat[:, i_mul].to(dtype=dtype)

        placements_cp_model_io_mul = {
            f"{k}_{i_mul}": v for k, v in _placements["cp_model_io"].items() for i_mul in range(multiplicity)
        }
        placements_cp = _placements["cp_atom_features"] | placements_cp_model_io_mul

        placements_model_io_mul = {
            f"{k}_{i_mul}": v for k, v in _placements["model_io"].items() for i_mul in range(multiplicity)
        }
        placements_dp_cp = placements_model_io_mul

        io_feats = distribute_atom_features(
            inputs=inputs_io,
            placements_cp=placements_cp,
            placements_dp_cp=placements_dp_cp,
            device_mesh=manager.device_mesh_subgroups,
            cp_group=manager.group["cp"],
            multiplicities={"noise": multiplicity},
        )

        # Homogenize model I/O DTensors to match dataloader batch atom dim
        target_atoms_global = feats["atom_pad_mask"].shape[-1]
        for k in list(io_feats.keys()):
            if io_feats[k].shape[1] < target_atoms_global:
                io_feats[k] = pad_to_length(io_feats[k], dim=1, length=target_atoms_global)

        noise_dt = io_feats.pop("noise")
        batch = feats

    # Monkeypatch deterministic noise for distributed forward
    sigmas_device = sigmas_expected_host.to(device=manager.device, dtype=dtype)
    sigmas_dt = distribute_tensor(sigmas_device, manager.device_mesh_subgroups, (Shard(0), Replicate(), Replicate()))
    monkeypatch.setattr(module.structure_module, "noise_distribution", lambda bs, dtype=None: sigmas_dt)
    monkeypatch.setattr(distributed_diffusion_module, "create_distributed_randn", lambda *a, **kw: noise_dt)

    # Run distributed training_step
    loss = module.training_step(batch, batch_idx=0)

    # Assert 1: loss matches serial
    loss_local = loss.to_local() if isinstance(loss, DTensor) else loss
    serial_loss_device = serial_loss_host.to(device=manager.device, dtype=dtype)
    torch.testing.assert_close(
        loss_local,
        serial_loss_device,
        msg=lambda msg: f"Rank {rank}: Loss mismatch: {msg}",
    )

    # Backward
    loss_local.backward()

    # on_after_backward redistributes gradients to Replicate and logs grad_norm metrics
    module.on_after_backward()

    # Assert 1b: logged metrics parity (CSVLogger output).
    # Compared after backward + on_after_backward so grad_norm metrics are included.
    assert len(dist_log.metrics) > 0, f"Rank {rank}: No metrics logged — test is vacuous"
    assert set(dist_log.metrics.keys()) == set(serial_log_metrics_host.keys()), (
        f"Rank {rank}: Logged metric keys differ. "
        f"Serial: {sorted(serial_log_metrics_host.keys())}, "
        f"Distributed: {sorted(dist_log.metrics.keys())}"
    )
    for key in sorted(serial_log_metrics_host.keys()):
        torch.testing.assert_close(
            torch.tensor(dist_log.metrics[key], dtype=torch.float64),
            torch.tensor(serial_log_metrics_host[key], dtype=torch.float64),
            msg=lambda msg, k=key: f"Rank {rank}: Logged metric mismatch for {k}: {msg}",
        )
    dist_log.flush(step=0)

    # Assert 2: per-parameter gradient parity (default fp64 tolerances).
    num_grads_checked = 0
    num_nonzero_grads = 0
    for name, param in module.named_parameters():
        canonical_name = name.replace("._serial.", ".")
        if canonical_name not in serial_grad_dict_host:
            continue
        expected_grad = serial_grad_dict_host[canonical_name].to(device=manager.device, dtype=dtype)
        assert param.grad is not None, f"Rank {rank}: Missing gradient for {canonical_name}"
        actual_grad = param.grad.full_tensor() if isinstance(param.grad, DTensor) else param.grad
        num_grads_checked += 1
        if expected_grad.abs().max().item() > 0:
            num_nonzero_grads += 1
        torch.testing.assert_close(
            actual_grad,
            expected_grad,
            msg=lambda msg, cn=canonical_name: (
                f"Rank {rank}: Gradient mismatch for {cn}. "
                f"Serial grad norm: {expected_grad.norm().item():.10f}, "
                f"Distributed grad norm: {actual_grad.norm().item():.10f}. {msg}"
            ),
        )

    assert num_grads_checked > 0, f"Rank {rank}: No gradients compared — test is vacuous"
    assert num_nonzero_grads > 0, f"Rank {rank}: All compared gradients are zero — test is vacuous"

    # Assert 3: optimizer step parity (if requested)
    if optimizer_step:
        optimizer = torch.optim.Adam(module.parameters(), lr=1e-3, betas=(0.9, 0.999))
        optimizer.step()

        num_params_checked = 0
        num_params_changed = 0
        for name, param in module.named_parameters():
            canonical_name = name.replace("._serial.", ".")
            if canonical_name not in serial_post_opt_dict_host:
                continue
            expected_val = serial_post_opt_dict_host[canonical_name].to(device=manager.device, dtype=dtype)
            actual_val = param.full_tensor() if isinstance(param, DTensor) else param.data
            num_params_checked += 1
            if name in initial_params:
                initial_val = initial_params[name].to(device=manager.device, dtype=dtype)
                if not torch.equal(actual_val, initial_val):
                    num_params_changed += 1
            torch.testing.assert_close(
                actual_val,
                expected_val,
                msg=lambda msg, cn=canonical_name: f"Rank {rank}: Post-optimizer mismatch for {cn}. {msg}",
            )
        assert num_params_checked > 0, f"Rank {rank}: No post-optimizer params compared"
        assert num_params_changed > 0, f"Rank {rank}: No parameters changed after optimizer step — test is vacuous"

    # Assert 4: cross-rank loss identity
    torch.distributed.barrier()


@pytest.mark.slow
@pytest.mark.parametrize(
    ("setup_env", "optimizer_step", "use_random_features"),
    [
        (((2, (2, 2)), True, "cuda", "ENV"), True, True),
        (((2, (2, 2)), True, "cuda", "ENV"), True, False),
    ],
    indirect=["setup_env"],
    ids=["cuda-dp2-cp2x2-random", "cuda-dp2-cp2x2-dataloader"],
)
def test_boltz2_training_step_parity(
    setup_env,
    optimizer_step,
    use_random_features,
    test_cp_training_base_data_dir_boltz2,
    canonical_mols_dir,
    tmp_path,
):
    """V16: Serial vs distributed training_step numerical parity.

    Verifies that the distributed Boltz2Distributed.training_step produces
    numerically identical loss, gradients, post-optimizer parameters, and
    logged metrics compared to the serial Boltz2.training_step, with all
    randomness sources controlled:

    1. Recycling: no_random_recycling_training=True (fixed recycling_steps)
    2. Diffusion noise: monkeypatched noise_distribution and randn_like / create_distributed_randn
    3. Coordinate augmentation: disabled
    4. Sampling steps: fixed (no sampling_steps_random)
    5. Dropout: 0.0 in all modules

    Both serial and distributed sessions inject a CSVLogger-backed
    ``_LogCapture`` so that ``self.log()`` calls are exercised (not
    silenced) and the logged metric keys/values are compared.

    Extends V13 (forward/backward parity) to the full training_step control
    flow including loss aggregation, recycling-step broadcasting, and
    gradient redistribution.
    """
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("CUDA not available")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"Need {world_size} GPUs, have {torch.cuda.device_count()}")

    dtype = torch.float64
    min_val_init = -0.1
    max_val_init = 0.1
    scale_glorot = 0.08

    seed = 42
    seed_by_rank(0, seed=seed)

    # Small model config — dropout=0 for determinism
    boltz2_model_params = create_boltz2_model_init_params(use_large_model=False)
    recycling_steps = 0
    multiplicity = 2

    boltz2_model_params["training_args"].recycling_steps = recycling_steps
    boltz2_model_params["training_args"].diffusion_multiplicity = multiplicity
    boltz2_model_params["training_args"].sampling_steps = -1  # not used in training path but read by training_step
    boltz2_model_params["no_random_recycling_training"] = True
    boltz2_model_params["predict_bfactor"] = True
    boltz2_model_params["training_args"].bfactor_loss_weight = 1.0
    boltz2_model_params["validate_structure"] = True

    size_cp = grid_group_sizes["cp"][0]
    B = grid_group_sizes["dp"]

    if use_random_features:
        n_tokens = 30 * size_cp
        W = boltz2_model_params["atoms_per_window_queries"]
        n_atoms_raw = n_tokens * 20
        n_atoms = ((n_atoms_raw + W - 1) // W) * W
        n_msa = max(size_cp * 2, 2)

        input_feats_global_fp64 = random_features(
            size_batch=B,
            n_tokens=n_tokens,
            n_atoms=n_atoms,
            n_msa=n_msa,
            atom_counts_per_token_range=(8, 20),
            device=device_type,
            float_value_range=(min_val_init, max_val_init),
            selected_keys=_BOLTZ2_SELECTED_KEYS,
            num_disto_bins=boltz2_model_params["num_bins"],
        )
        input_feats_global_fp64["msa"] = torch.randint(
            0, const.num_tokens, (B, n_msa, n_tokens), dtype=torch.int64, device=device_type
        )
        input_feats_global_fp64["disto_target"] = input_feats_global_fp64["disto_target"].unsqueeze(3)
        training_data_dir = None
    else:
        boltz2_model_params["num_bins"] = 64
        training_data_dir = _setup_training_data_7z64_8b2e(
            tmp_path / "training_data", test_cp_training_base_data_dir_boltz2
        )
        _mp_data = pytest.MonkeyPatch()
        cfg = setup_mock_training_datamodule_config(training_data_dir)
        cfg.batch_size = B
        cfg.samples_per_epoch = B
        cfg.moldir = str(canonical_mols_dir)
        cfg.return_train_symmetries = False
        for ds_cfg in cfg.datasets:
            ds_cfg.filters = None
        seed_by_rank(0, seed=seed)
        dm = Boltz2TrainingDataModuleSerial(cfg=cfg)
        _deterministic_getitem_monkeypatch(_mp_data, dm._train_set, base_seed=seed)
        dl = dm.train_dataloader()
        raw_batch = next(iter(dl))
        input_feats_global_fp64 = {}
        for k, v in raw_batch.items():
            if isinstance(v, torch.Tensor):
                input_feats_global_fp64[k] = v.to(
                    device=device_type, dtype=dtype if v.dtype.is_floating_point else v.dtype
                )
            elif isinstance(v, list):
                input_feats_global_fp64[k] = [
                    item.to(device=device_type, dtype=dtype if item.dtype.is_floating_point else item.dtype)
                    if isinstance(item, torch.Tensor)
                    else item
                    for item in v
                ]
            else:
                input_feats_global_fp64[k] = v
        n_atoms = input_feats_global_fp64["atom_pad_mask"].shape[-1]
        _mp_data.undo()

        if "token_pair_pad_mask" not in input_feats_global_fp64:
            tpm = input_feats_global_fp64["token_pad_mask"]
            input_feats_global_fp64["token_pair_pad_mask"] = tpm[:, :, None] * tpm[:, None, :]

    # Create serial model with deterministic init
    serial_model = SerialBoltz2(**boltz2_model_params)
    init_module_params_glorot(serial_model, gain=scale_glorot)
    serial_model.apply(SetModuleInfValues())
    serial_model.structure_module.coordinate_augmentation = False
    serial_model = serial_model.to(dtype=dtype, device=device_type)
    serial_model.train()

    # Save state dict for distributed model (before serial forward mutates state)
    module_state_dict = {k: v.detach().clone().cpu() for k, v in serial_model.state_dict().items()}

    # Pre-generate deterministic sigmas and non-zero noise for serial model
    sigmas_serial = serial_model.structure_module.noise_distribution(B * multiplicity).to(dtype=dtype)
    noise_serial = torch.empty(B * multiplicity, n_atoms, 3, device=device_type, dtype=dtype)
    init_tensors_uniform([noise_serial], low=min_val_init, high=max_val_init)
    atom_pad_mask_mul = input_feats_global_fp64["atom_pad_mask"][:, :, None].repeat_interleave(multiplicity, 0)
    noise_serial = noise_serial * atom_pad_mask_mul

    _serial_mp = pytest.MonkeyPatch()
    serial_csv_logger = CSVLogger(save_dir=str(tmp_path), name="serial")
    serial_log = _LogCapture(serial_csv_logger)
    _serial_mp.setattr(serial_model, "log", serial_log)
    _serial_mp.setattr(serial_model, "training_log", lambda *a, **kw: None)
    _serial_mp.setattr(serial_model.structure_module, "noise_distribution", lambda bs, dtype=None: sigmas_serial)
    _serial_mp.setattr(serial_diffusion_v2_module.torch, "randn_like", lambda t: noise_serial.to(t))
    # Monkeypatch serial smooth_lddt_loss to use dense pairwise distances.
    # The original serial code uses sparse indexing (nonzero + F.pairwise_distance)
    # which creates a different autograd backward graph than the distributed dense
    # matrix computation (replicate_to_shard_outer_op CDIST). The different
    # accumulation patterns amplify ~1e-12 forward differences into ~4.5e-7
    # gradient errors. Using dense distances here aligns the backward structure.
    import boltz.model.loss.diffusionv2 as _serial_loss_mod

    def _smooth_lddt_loss_dense(
        pred_coords,
        true_coords,
        is_nucleotide,
        coords_mask=None,
        nucleic_acid_cutoff=30.0,
        other_cutoff=15.0,
        multiplicity=1,
    ):
        compute_dtype = torch.promote_types(pred_coords.dtype, torch.float32)
        N = pred_coords.shape[1]
        lddt = []
        for i in range(true_coords.shape[0]):
            true_dists = torch.cdist(true_coords[i], true_coords[i])

            is_nuc_i = is_nucleotide[i // multiplicity]
            mask_i = coords_mask[i // multiplicity]

            is_nuc_pair = is_nuc_i.unsqueeze(-1).expand(-1, is_nuc_i.shape[-1])

            mask = is_nuc_pair * (true_dists < nucleic_acid_cutoff).to(compute_dtype)
            mask += (1 - is_nuc_pair) * (true_dists < other_cutoff).to(compute_dtype)
            mask *= 1 - torch.eye(N, device=pred_coords.device)
            mask *= mask_i.unsqueeze(-1)
            mask *= mask_i.unsqueeze(-2)

            diff = pred_coords[i].unsqueeze(0) - pred_coords[i].unsqueeze(1)
            pred_dists = (diff * diff).sum(-1).add(1e-30).sqrt()

            dist_diff = (true_dists - pred_dists).abs()

            eps = (
                torch.sigmoid(0.5 - dist_diff)
                + torch.sigmoid(1.0 - dist_diff)
                + torch.sigmoid(2.0 - dist_diff)
                + torch.sigmoid(4.0 - dist_diff)
            ) * 0.25

            lddt_i = (eps * mask).sum() / (mask.sum() + 1e-5)
            lddt.append(lddt_i)

        return 1 - sum(lddt) / len(lddt)

    _serial_mp.setattr(_serial_loss_mod, "smooth_lddt_loss", _smooth_lddt_loss_dense)
    _serial_mp.setattr(serial_diffusion_v2_module, "smooth_lddt_loss", _smooth_lddt_loss_dense)

    # Save coords — serial forward mutates feats["coords"] in-place (flattens ensemble dim)
    coords_backup = input_feats_global_fp64["coords"].detach().clone()
    original_feat_keys = set(input_feats_global_fp64.keys())

    # Run serial training_step
    serial_loss = serial_model.training_step(input_feats_global_fp64, batch_idx=0)
    assert serial_loss is not None, "Serial training_step returned None"
    assert serial_loss.isfinite(), f"Serial loss is not finite: {serial_loss.item()}"

    # Backward on serial model
    serial_loss.backward()

    # on_after_backward logs grad_norm metrics (component-wise and global)
    serial_model.on_after_backward()

    # Flush serial CSVLogger and collect logged metrics (after backward so
    # grad_norm metrics are included alongside the training_step metrics)
    serial_log.flush(step=0)
    serial_log_metrics = dict(serial_log.metrics)
    assert len(serial_log_metrics) > 0, "Serial model logged no metrics — validate_structure may be False"

    # Collect serial gradients
    serial_grad_dict = {}
    for name, param in serial_model.named_parameters():
        if param.grad is not None:
            serial_grad_dict[name] = param.grad.detach().clone().cpu()

    # Collect serial post-optimizer parameters (if needed)
    serial_post_opt_dict = {}
    if optimizer_step:
        # Use a simple Adam (no weight decay) for reproducible parity.
        # configure_optimizers uses AdamW with lr schedulers which complicates comparison.
        serial_optimizer = torch.optim.Adam(serial_model.parameters(), lr=1e-3, betas=(0.9, 0.999))
        serial_optimizer.step()
        for name, param in serial_model.named_parameters():
            serial_post_opt_dict[name] = param.data.detach().clone().cpu()

    serial_loss_host = serial_loss.detach().clone().cpu()
    sigmas_host = sigmas_serial.detach().clone().cpu()
    noise_host = noise_serial.detach().clone().cpu()

    # Restore features — serial forward mutates coords in-place and may add keys
    input_feats_global_fp64["coords"] = coords_backup
    for key in list(input_feats_global_fp64.keys()):
        if key not in original_feat_keys:
            del input_feats_global_fp64[key]
    _serial_mp.undo()

    # Move features to CPU for spawned workers
    input_feats_host = {}
    for k, v in input_feats_global_fp64.items():
        if isinstance(v, torch.Tensor):
            input_feats_host[k] = v.detach().to(device="cpu", copy=True)
        elif isinstance(v, list):
            input_feats_host[k] = [
                item.detach().to(device="cpu", copy=True) if isinstance(item, torch.Tensor) else item for item in v
            ]
        else:
            input_feats_host[k] = v

    # Spawn parallel test
    spawn_multiprocessing(
        _worker_training_step_parity,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        module_state_dict,
        boltz2_model_params,
        input_feats_host,
        sigmas_host,
        noise_host,
        serial_loss_host,
        serial_log_metrics,
        serial_grad_dict,
        serial_post_opt_dict,
        optimizer_step,
        use_random_features,
        training_data_dir,
        canonical_mols_dir,
        seed,
    )


# ====================================================================== #
#  V15: setup() – validator wiring, predict no-op, datamodule=None        #
# ====================================================================== #


def _worker_setup(
    rank,
    grid_group_sizes,
    device_type,
    backend,
    env_map=None,
):
    """Worker: verify setup() wires validators and handles predict stage."""
    monkeypatch, dm = _init_distributed(rank, grid_group_sizes, device_type, backend, env_map)

    validators = [DistributedRCSBValidator(val_names=["RCSB"], confidence_prediction=False, physicalism_metrics=True)]
    serial_model = _create_minimal_serial_boltz2(validate_structure=True, validators=validators, num_val_datasets=1).to(
        dm.device
    )
    dist_model = Boltz2Distributed(serial_model, dm)

    # --- setup("predict") does not populate validator_mapper ---
    dist_model._trainer = SimpleNamespace(
        datamodule=SimpleNamespace(
            val_group_mapper={0: {"label": "RCSB", "symmetry_correction": False}},
        ),
    )
    dist_model.setup("predict")
    assert len(dist_model.validator_mapper) == 0, "predict stage should not wire validators"

    # --- setup("fit") populates validator_mapper from trainer.datamodule ---
    dist_model.setup("fit")
    assert len(dist_model.validator_mapper) == 1, "fit stage should populate validator_mapper"
    assert 0 in dist_model.validator_mapper
    assert isinstance(dist_model.validator_mapper[0], DistributedRCSBValidator)

    # --- val_group_mapper updated from datamodule ---
    assert len(dist_model.val_group_mapper) == 1
    assert dist_model.val_group_mapper[0]["label"] == "RCSB"
    assert dist_model.val_group_mapper[0]["symmetry_correction"] is False

    del dist_model, serial_model

    # --- setup without datamodule is a no-op ---
    validators2 = [DistributedRCSBValidator(val_names=["RCSB"], confidence_prediction=False, physicalism_metrics=True)]
    serial2 = _create_minimal_serial_boltz2(validate_structure=True, validators=validators2, num_val_datasets=1).to(
        dm.device
    )
    dist2 = Boltz2Distributed(serial2, dm)
    dist2._trainer = SimpleNamespace(datamodule=None)
    dist2.setup("fit")
    assert len(dist2.validator_mapper) == 0, "No datamodule should leave validator_mapper empty"
    del dist2, serial2

    # --- setup with validate_structure=False is a no-op regardless of stage ---
    serial3 = _create_minimal_serial_boltz2(validate_structure=False).to(dm.device)
    dist3 = Boltz2Distributed(serial3, dm)
    dist3._trainer = SimpleNamespace(
        datamodule=SimpleNamespace(
            val_group_mapper={0: {"label": "RCSB", "symmetry_correction": False}},
        ),
    )
    dist3.setup("fit")
    assert not hasattr(dist3, "validator_mapper"), "validate_structure=False model should not have validator_mapper"
    del dist3, serial3

    # --- setup("fit") with mismatched num_val_datasets raises AssertionError ---
    validators_mismatch = [
        DistributedRCSBValidator(val_names=["RCSB"], confidence_prediction=False, physicalism_metrics=True)
    ]
    serial_mm = _create_minimal_serial_boltz2(
        validate_structure=True, validators=validators_mismatch, num_val_datasets=1
    ).to(dm.device)
    dist_mm = Boltz2Distributed(serial_mm, dm)
    dist_mm._trainer = SimpleNamespace(
        datamodule=SimpleNamespace(
            val_group_mapper={
                0: {"label": "RCSB", "symmetry_correction": False},
                1: {"label": "EXTRA", "symmetry_correction": False},
            },
        ),
    )
    with pytest.raises(AssertionError, match="num_val_datasets"):
        dist_mm.setup("fit")
    del dist_mm, serial_mm

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env",
    [
        ((2, (2, 2)), True, "cuda", "ENV"),
    ],
    indirect=("setup_env",),
    ids=lambda val: f"{val[2]}-dp:{val[0][0]}-cp{'x'.join(map(str, val[0][1]))}",
)
def test_setup(setup_env):
    """V15: setup() wires validators correctly and predict stage is a no-op.

    Verifies:
      - setup("predict") does not populate validator_mapper
      - setup("fit") populates validator_mapper from trainer.datamodule
      - setup without datamodule is a no-op
      - validate_structure=False means setup("fit") doesn't touch validators
      - Mismatched num_val_datasets raises AssertionError
    """
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if not torch.cuda.is_available():
        pytest.skip("CUDA not available")
    if torch.cuda.device_count() < world_size:
        pytest.skip(f"Need {world_size} GPUs, have {torch.cuda.device_count()}")

    spawn_multiprocessing(
        _worker_setup,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
    )


# ====================================================================== #
#  V18: validation_step and on_validation_epoch_end (real forward pass)   #
# ====================================================================== #


def _worker_validation_step_parity(
    rank,
    grid_group_sizes,
    device_type,
    backend,
    boltz2_model_params,
    module_state_dict,
    input_feats_host,
    noise_host_list,
    serial_per_sample,
    serial_epoch_end_metrics,
    env_per_rank=None,
):
    """V18 multi-rank worker: verify distributed validation_step matches serial.

    Phase 1: Compare raw accumulated metrics after validation_step against
    the serial per-sample reference for this DP rank's sample.

    Phase 2: Compare aggregated metrics after on_validation_epoch_end
    against the serial epoch-end reference.
    """
    monkeypatch = pytest.MonkeyPatch()
    if env_per_rank is not None:
        for var_name, value in env_per_rank.items():
            monkeypatch.setenv(var_name, f"{rank}" if value == "<INPUT_RANK>" else value)

    dtype = torch.float64
    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    manager = DistributedManager()

    dp_rank = manager.group_rank["dp"]

    reference_module = SerialBoltz2(**boltz2_model_params)
    reference_module = reference_module.to(dtype=dtype)
    reference_module.load_state_dict(module_state_dict)
    reference_module.structure_module.coordinate_augmentation = False
    reference_module.apply(SetModuleInfValues())
    reference_module = reference_module.to(device=manager.device)
    module = Boltz2Distributed(reference_module, manager)
    module.eval()

    # Wire validators via setup
    symmetry_correction = boltz2_model_params["validation_args"].symmetry_correction
    num_validators = boltz2_model_params["num_val_datasets"]
    val_names = [v.val_names[0] for v in boltz2_model_params["validators"]]
    worker_val_group_mapper = {
        vi: {"label": val_names[vi], "symmetry_correction": symmetry_correction} for vi in range(num_validators)
    }
    module._trainer = SimpleNamespace(
        datamodule=SimpleNamespace(val_group_mapper=worker_val_group_mapper),
        sanity_checking=False,
    )
    module.setup("fit")
    assert len(module.validator_mapper) == num_validators

    # ------------------------------------------------------------------
    # Distribute token + MSA + atom features for this rank's sample
    # ------------------------------------------------------------------
    diffusion_samples = boltz2_model_params["validation_args"].diffusion_samples
    num_sampling_steps = boltz2_model_params["validation_args"].sampling_steps

    host_tensor_keys = {k for k, v in input_feats_host.items() if isinstance(v, torch.Tensor)}
    _placements = get_feature_placements(
        token_keys=host_tensor_keys,
        msa_keys=host_tensor_keys,
        atom_keys=host_tensor_keys
        & {
            "ref_pos",
            "atom_resolved_mask",
            "ref_element",
            "ref_charge",
            "ref_atom_name_chars",
            "ref_space_uid",
            "coords",
            "atom_pad_mask",
            "atom_to_token",
            "pair_mask",
            "atom_counts_per_token",
            "token_to_rep_atom",
            "bfactor",
            "plddt",
        },
        model_io_keys={"noise"},
        model_io_fp32_keys=set(),
    )

    token_msa_placements = _placements["token_features"] | _placements["msa_features"]

    if manager.group_rank["world"] == 0:
        input_feats_token_msa_global = {
            k: v.to(device=manager.device, dtype=dtype if v.dtype.is_floating_point else v.dtype)
            for k, v in input_feats_host.items()
            if k in token_msa_placements
        }
    else:
        input_feats_token_msa_global = None

    feats_token_msa = distribute_features(
        input_feats_token_msa_global,
        token_msa_placements,
        manager.group["world"],
        manager.group_ranks["world"][0],
        manager.device_mesh_subgroups,
    )

    size_batch = input_feats_host["atom_pad_mask"].shape[0]
    inputs_atom = {
        k: v.to(dtype=dtype if v.dtype.is_floating_point else v.dtype)
        for k, v in input_feats_host.items()
        if k in _placements["cp_atom_features"]
    }

    n_noise = 1 + num_sampling_steps
    for i_noise, noise_host in enumerate(noise_host_list):
        unflat = noise_host.unflatten(0, (size_batch, diffusion_samples))
        for i_mul in range(diffusion_samples):
            inputs_atom[f"_noise_{i_noise}_{i_mul}"] = unflat[:, i_mul].to(dtype=dtype)

    noise_cp_placements = {}
    noise_placements = {}
    for i_noise in range(n_noise):
        for i_mul in range(diffusion_samples):
            key = f"_noise_{i_noise}_{i_mul}"
            noise_cp_placements[key] = _placements["cp_model_io"]["noise"]
            noise_placements[key] = _placements["model_io"]["noise"]

    feats_and_noise = distribute_atom_features(
        inputs=inputs_atom,
        placements_cp=_placements["cp_atom_features"] | noise_cp_placements,
        placements_dp_cp=_placements["atom_features"] | noise_placements,
        device_mesh=manager.device_mesh_subgroups,
        cp_group=manager.group["cp"],
        multiplicities={f"_noise_{i}": diffusion_samples for i in range(n_noise)},
    )

    noise_dts = []
    for i_noise in range(n_noise):
        noise_dts.append(feats_and_noise.pop(f"_noise_{i_noise}"))
    init_noise_dt = noise_dts[0]
    step_noise_dts = noise_dts[1:]

    feats_dt = {**feats_token_msa, **feats_and_noise}

    # Add symmetry features as non-sharded plain tensors/lists (DP-local slice)
    if symmetry_correction:
        _SYM_KEYS = {
            "all_coords",
            "all_resolved_mask",
            "crop_to_all_atom_map",
            "chain_swaps",
            "amino_acids_symmetries",
            "ligand_symmetries",
        }
        for sk in _SYM_KEYS:
            if sk not in input_feats_host:
                continue
            val = input_feats_host[sk]
            if isinstance(val, torch.Tensor):
                feats_dt[sk] = val[dp_rank : dp_rank + 1].to(
                    device=manager.device, dtype=dtype if val.dtype.is_floating_point else val.dtype
                )
            elif isinstance(val, list):
                elem = val[dp_rank]
                if isinstance(elem, torch.Tensor):
                    feats_dt[sk] = elem.unsqueeze(0).to(
                        device=manager.device, dtype=dtype if elem.dtype.is_floating_point else elem.dtype
                    )
                else:
                    feats_dt[sk] = [elem]

    # Non-sharded physicalism keys (clash + PB): DP-sliced, not CP-sharded
    _PHYS_KEYS = LIGAND_GEOMETRY_FEATURES | {"chain_symmetries"}
    for pk in _PHYS_KEYS:
        if pk not in input_feats_host:
            continue
        val = input_feats_host[pk]
        if isinstance(val, torch.Tensor):
            feats_dt[pk] = val[dp_rank : dp_rank + 1].to(
                device=manager.device, dtype=dtype if val.dtype.is_floating_point else val.dtype
            )
        elif isinstance(val, list):
            elem = val[dp_rank]
            if isinstance(elem, torch.Tensor):
                feats_dt[pk] = elem.unsqueeze(0).to(
                    device=manager.device, dtype=dtype if elem.dtype.is_floating_point else elem.dtype
                )
            else:
                feats_dt[pk] = [elem]

    # ------------------------------------------------------------------
    # Monkeypatch distributed sample() for determinism
    # ------------------------------------------------------------------
    _orig_center_random_augmentation = distributed_diffusion_module.center_random_augmentation

    def _centering_only_augmentation(atom_coords, atom_mask, **kwargs):
        kwargs["augmentation"] = False
        kwargs["centering"] = True
        return _orig_center_random_augmentation(atom_coords, atom_mask, **kwargs)

    _dt_randn_calls = []
    _dt_randn_sequence = [init_noise_dt] + step_noise_dts

    def _fixed_create_distributed_randn(shape, device_mesh, placements, dtype=torch.float32, scale=1.0):
        idx = len(_dt_randn_calls)
        _dt_randn_calls.append(idx)
        noise_dt = _dt_randn_sequence[idx]
        if scale != 1.0:
            noise_dt = scalar_tensor_op(scale, noise_dt, ElementwiseOp.PROD)
        return noise_dt

    monkeypatch.setattr(distributed_diffusion_module, "center_random_augmentation", _centering_only_augmentation)
    monkeypatch.setattr(distributed_diffusion_module, "create_distributed_randn", _fixed_create_distributed_randn)

    # ------------------------------------------------------------------
    # Phase 1: Run validation_step (validator 0), compare accumulated metrics
    # ------------------------------------------------------------------
    feats_dt["idx_dataset"] = [torch.tensor([0], device=manager.device)]

    with torch.no_grad():
        module.validation_step(feats_dt, batch_idx=0)

    validator = module.validator_mapper[0]
    fm = validator.folding_metrics
    val_idx = 0

    serial_ref = serial_per_sample[dp_rank]
    compared_phase1 = 0

    disto_loss_metric = fm["disto_loss"][val_idx]["disto_loss"]
    if disto_loss_metric.weight > 0:
        dist_disto_loss = disto_loss_metric.compute().item()
        serial_disto_loss_avg = sum(s["disto_loss"] for s in serial_per_sample) / len(serial_per_sample)
        torch.testing.assert_close(
            torch.tensor(dist_disto_loss, dtype=dtype),
            torch.tensor(serial_disto_loss_avg, dtype=dtype),
            msg=lambda msg: f"Rank {rank}: Phase 1 disto_loss mismatch: {msg}",
        )
        compared_phase1 += 1

    for key in [*const.out_types, "pocket_ligand_protein", "contact_protein_protein"]:
        if key in fm["disto_lddt"][val_idx]:
            metric = fm["disto_lddt"][val_idx][key]
            if metric.weight > 0 and key in serial_ref.get("disto_lddt", {}):
                dist_val = metric.compute().item()
                serial_val = serial_ref["disto_lddt"][key]
                torch.testing.assert_close(
                    torch.tensor(dist_val, dtype=dtype),
                    torch.tensor(serial_val, dtype=dtype),
                    msg=lambda msg, k=key: f"Rank {rank}: Phase 1 disto_lddt_{k} mismatch: {msg}",
                )
                compared_phase1 += 1

    for key in [*const.out_types, "pocket_ligand_protein", "contact_protein_protein"]:
        if key in fm["lddt"][val_idx]:
            metric = fm["lddt"][val_idx][key]
            if metric.weight > 0 and key in serial_ref.get("lddt", {}):
                dist_val = metric.compute().item()
                serial_val = serial_ref["lddt"][key]
                torch.testing.assert_close(
                    torch.tensor(dist_val, dtype=dtype),
                    torch.tensor(serial_val, dtype=dtype),
                    msg=lambda msg, k=key: f"Rank {rank}: Phase 1 lddt_{k} mismatch: {msg}",
                )
                compared_phase1 += 1

    for key in [*const.out_types, "pocket_ligand_protein", "contact_protein_protein"]:
        if key in fm["complex_lddt"][val_idx]:
            metric = fm["complex_lddt"][val_idx][key]
            if metric.weight > 0 and key in serial_ref.get("complex_lddt", {}):
                dist_val = metric.compute().item()
                serial_val = serial_ref["complex_lddt"][key]
                torch.testing.assert_close(
                    torch.tensor(dist_val, dtype=dtype),
                    torch.tensor(serial_val, dtype=dtype),
                    msg=lambda msg, k=key: f"Rank {rank}: Phase 1 complex_lddt_{k} mismatch: {msg}",
                )
                compared_phase1 += 1

    assert compared_phase1 >= 3, f"Rank {rank}: Phase 1 compared only {compared_phase1} metrics — test may be vacuous"

    # Run remaining validators for Phase 2 accumulation
    for vi in range(1, num_validators):
        feats_dt["idx_dataset"] = [torch.tensor([vi], device=manager.device)]

        _dt_randn_calls_vi = []
        _dt_randn_sequence_vi = [init_noise_dt] + step_noise_dts

        def _fixed_randn_vi(
            shape,
            device_mesh,
            placements,
            dtype=torch.float32,
            scale=1.0,
            _seq=_dt_randn_sequence_vi,
            _calls=_dt_randn_calls_vi,
        ):
            idx = len(_calls)
            _calls.append(idx)
            noise_dt = _seq[idx]
            if scale != 1.0:
                noise_dt = scalar_tensor_op(scale, noise_dt, ElementwiseOp.PROD)
            return noise_dt

        monkeypatch.setattr(distributed_diffusion_module, "create_distributed_randn", _fixed_randn_vi)

        # In Phase 1, only 1 validator is used. For tests with multiple validators, validation_step is called for other validators.
        with torch.no_grad():
            module.validation_step(feats_dt, batch_idx=vi)

    # ------------------------------------------------------------------
    # Phase 2: on_validation_epoch_end and compare aggregated metrics
    # ------------------------------------------------------------------
    dist_log = _LogCapture(CSVLogger(save_dir=tempfile.mkdtemp(), name=f"dist_val_rank{rank}"))
    monkeypatch.setattr(module, "log", dist_log)

    module.on_validation_epoch_end()

    _forward_dependent_prefixes = ("val/lddt", "val/complex_lddt", "val/clash", "val/pb", "val/rmsd")
    compared_phase2 = 0
    for key in sorted(serial_epoch_end_metrics):
        if key in dist_log.metrics:
            got = torch.tensor(dist_log.metrics[key], dtype=dtype)
            exp = torch.tensor(serial_epoch_end_metrics[key], dtype=dtype)
            if any(key.startswith(p) for p in _forward_dependent_prefixes):
                torch.testing.assert_close(
                    got,
                    exp,
                    msg=lambda msg, k=key: f"Rank {rank}: Phase 2 epoch-end metric '{k}' mismatch: {msg}",
                )
            else:
                torch.testing.assert_close(
                    got,
                    exp,
                    msg=lambda msg, k=key: f"Rank {rank}: Phase 2 epoch-end metric '{k}' mismatch: {msg}",
                )
            compared_phase2 += 1

    assert compared_phase2 >= 3, f"Rank {rank}: Phase 2 compared only {compared_phase2} metrics — test may be vacuous"

    required_base_metrics = ("val/lddt", "val/disto_lddt", "val/complex_lddt")
    for base_metric in required_base_metrics:
        for vn in val_names:
            suffix = "" if vn == "RCSB" else f"__{vn}"
            required_metric = f"{base_metric}{suffix}"
            assert required_metric in serial_epoch_end_metrics, (
                f"Rank {rank}: serial epoch-end metrics missing '{required_metric}' — "
                f"available: {sorted(serial_epoch_end_metrics)}"
            )
            assert required_metric in dist_log.metrics, (
                f"Rank {rank}: distributed epoch-end metrics missing '{required_metric}' — "
                f"available: {sorted(dist_log.metrics)}"
            )

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.parametrize(
    ("setup_env", "use_random_features", "symmetry_correction", "num_validators"),
    [
        (((2, (2, 2)), True, "cuda", "ENV"), True, False, 2),
        (((2, (2, 2)), True, "cuda", "ENV"), False, True, 1),
    ],
    indirect=["setup_env"],
    ids=["cuda-dp2-cp2x2-random-2val", "cuda-dp2-cp2x2-dataloader-sc-1val"],
)
def test_boltz2_validation_step_parity(
    setup_env,
    use_random_features,
    symmetry_correction,
    num_validators,
    test_cp_training_base_data_dir_boltz2,
    canonical_mols_dir,
    tmp_path,
):
    """V18: validation_step parity between distributed and serial Boltz2.

    Two-phase comparison:
      Phase 1: After validation_step, compare validator MeanMetric values
        (disto_loss, disto_lddt_*, lddt_*, complex_lddt_*) between serial
        (per-sample) and distributed (per DP-rank sample).
      Phase 2: After on_validation_epoch_end (which DP-all-reduces
        MeanMetric internals), compare aggregated logged metrics between
        serial and distributed.

    Uses DP=2, CP=(2,2) -> 8 GPUs, FP64 with default tolerances.
    """
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("CUDA not available")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"Need {world_size} GPUs, have {torch.cuda.device_count()}")

    if use_random_features and symmetry_correction:
        pytest.skip("Symmetry correction not supported with random features")

    if num_validators > 1 and not use_random_features:
        pytest.skip("Multiple validators require multiple val datasets; only supported with use_random_features=True")

    dtype = torch.float64
    min_val_init = -0.01
    max_val_init = 0.01
    scale_glorot = 0.05

    num_sampling_steps = 2
    diffusion_samples = 1

    seed = 42
    seed_by_rank(0, seed=seed)

    boltz2_model_params = create_boltz2_model_init_params(use_large_model=False)
    boltz2_model_params["diffusion_process_args"]["alignment_reverse_diff"] = False
    boltz2_model_params["validate_structure"] = True
    val_names = [f"RCSB_{i}" for i in range(num_validators)] if num_validators > 1 else ["RCSB"]
    boltz2_model_params["validators"] = [
        RCSBValidator(val_names=[vn], confidence_prediction=False, physicalism_metrics=False) for vn in val_names
    ]
    boltz2_model_params["num_val_datasets"] = num_validators
    boltz2_model_params["confidence_prediction"] = False
    boltz2_model_params["validation_args"] = _make_validation_args(
        recycling_steps=0,
        sampling_steps=num_sampling_steps,
        diffusion_samples=diffusion_samples,
        symmetry_correction=symmetry_correction,
    )

    size_cp = grid_group_sizes["cp"][0]
    B = grid_group_sizes["dp"]

    if use_random_features:
        n_atoms_per_token_min = 8
        n_atoms_per_token_max = 20
        n_tokens = 30 * size_cp
        W = boltz2_model_params["atoms_per_window_queries"]
        n_atoms_raw = n_tokens * n_atoms_per_token_max
        n_atoms = ((n_atoms_raw + W - 1) // W) * W
        n_msa = size_cp * 2

        assert n_atoms % size_cp == 0
        assert n_atoms % W == 0

        input_feats_global_fp64 = random_features(
            size_batch=B,
            n_tokens=n_tokens,
            n_atoms=n_atoms,
            n_msa=n_msa,
            atom_counts_per_token_range=(n_atoms_per_token_min, n_atoms_per_token_max),
            device=device_type,
            float_value_range=(min_val_init, max_val_init),
            selected_keys=_BOLTZ2_SELECTED_KEYS,
            num_disto_bins=boltz2_model_params["num_bins"],
        )
        input_feats_global_fp64["msa"] = torch.randint(
            0, const.num_tokens, (B, n_msa, n_tokens), dtype=torch.int64, device=device_type
        )

        input_feats_global_fp64["disto_target"] = input_feats_global_fp64["disto_target"].unsqueeze(3)

        token_to_rep_atom = input_feats_global_fp64["token_to_rep_atom"]
        coords = input_feats_global_fp64["coords"]
        disto_coords_ensemble = torch.bmm(token_to_rep_atom.to(dtype=dtype), coords[:, 0])
        input_feats_global_fp64["disto_coords_ensemble"] = disto_coords_ensemble

        input_feats_global_fp64["connections_edge_index"] = [
            torch.empty(2, 0, dtype=torch.long, device=device_type) for _ in range(B)
        ]
        input_feats_global_fp64["chain_symmetries"] = [[] for _ in range(B)]
    else:
        boltz2_model_params["num_bins"] = 64
        training_data_dir = _setup_training_data_7z64_8b2e(
            tmp_path / "training_data", test_cp_training_base_data_dir_boltz2
        )
        cfg = setup_mock_training_datamodule_config(training_data_dir)
        cfg.overfit = B
        cfg.samples_per_epoch = B
        cfg.moldir = str(canonical_mols_dir)
        cfg.return_train_symmetries = symmetry_correction
        cfg.pad_to_max_tokens = True
        cfg.pad_to_max_atoms = True
        cfg.pad_to_max_seqs = True
        W = boltz2_model_params["atoms_per_window_queries"]
        token_align = size_cp
        atom_align = math.lcm(W, size_cp)
        cfg.max_tokens = ((cfg.max_tokens + token_align - 1) // token_align) * token_align
        cfg.max_atoms = ((cfg.max_atoms + atom_align - 1) // atom_align) * atom_align
        cfg.max_seqs = ((cfg.max_seqs + size_cp - 1) // size_cp) * size_cp
        for ds_cfg in cfg.datasets:
            ds_cfg.filters = None
        seed_by_rank(0, seed=seed)
        dm = Boltz2TrainingDataModuleSerial(cfg=cfg)
        dl = dm.val_dataloader()
        dl_iter = iter(dl)
        raw_samples = [next(dl_iter) for _ in range(B)]

        def _unwrap_bs1(v):
            if isinstance(v, torch.Tensor):
                return v.squeeze(0)
            if isinstance(v, list) and len(v) == 1:
                return v[0]
            return v

        raw_batch = collate([{k: _unwrap_bs1(v) for k, v in s.items()} for s in raw_samples])

        input_feats_global_fp64 = {}
        for k, v in raw_batch.items():
            if isinstance(v, torch.Tensor):
                input_feats_global_fp64[k] = v.to(
                    device=device_type, dtype=dtype if v.dtype.is_floating_point else v.dtype
                )
            elif isinstance(v, list):
                input_feats_global_fp64[k] = [
                    item.to(device=device_type, dtype=dtype if item.dtype.is_floating_point else item.dtype)
                    if isinstance(item, torch.Tensor)
                    else item
                    for item in v
                ]
            else:
                input_feats_global_fp64[k] = v
        n_atoms = input_feats_global_fp64["atom_pad_mask"].shape[-1]

        if "token_pair_pad_mask" not in input_feats_global_fp64:
            tpm = input_feats_global_fp64["token_pad_mask"]
            input_feats_global_fp64["token_pair_pad_mask"] = tpm[:, :, None] * tpm[:, None, :]

    # ------------------------------------------------------------------
    # Slice global batch into individual samples for serial validation
    # ------------------------------------------------------------------
    def _slice_batch(feats, idx):
        batch_i = {}
        for k, v in feats.items():
            if isinstance(v, torch.Tensor):
                batch_i[k] = v[idx : idx + 1].clone()
            elif isinstance(v, list):
                elem = v[idx]
                if isinstance(elem, torch.Tensor):
                    batch_i[k] = elem.unsqueeze(0).clone()
                else:
                    batch_i[k] = [elem]
            else:
                batch_i[k] = v
        batch_i["idx_dataset"] = torch.tensor([0], device=device_type)
        return batch_i

    num_val_samples = B

    # ------------------------------------------------------------------
    # Build serial model
    # ------------------------------------------------------------------
    reference_module = SerialBoltz2(**boltz2_model_params)
    init_module_params_glorot(reference_module, gain=scale_glorot)
    reference_module.apply(SetModuleInfValues())
    reference_module.structure_module.coordinate_augmentation = False
    module_state_dict = reference_module.state_dict()
    reference_module = reference_module.to(dtype=dtype, device=device_type).eval()

    serial_validators = []
    reference_module.val_group_mapper = {}
    reference_module.validator_mapper = {}
    for vi in range(num_validators):
        vn = val_names[vi]
        v = RCSBValidator(val_names=[vn], confidence_prediction=False, physicalism_metrics=True)
        v = v.to(device=device_type, dtype=dtype)
        serial_validators.append(v)
        reference_module.val_group_mapper[vi] = {"label": vn, "symmetry_correction": symmetry_correction}
        reference_module.validator_mapper[vi] = v

    # ------------------------------------------------------------------
    # Pre-generate deterministic noise for sampling
    # ------------------------------------------------------------------
    _B_M = B * diffusion_samples
    init_noise = torch.empty((_B_M, n_atoms, 3), device=device_type, dtype=dtype)
    step_noise_list = [
        torch.empty((_B_M, n_atoms, 3), device=device_type, dtype=dtype) for _ in range(num_sampling_steps)
    ]
    init_tensors_uniform([init_noise, *step_noise_list], low=min_val_init, high=max_val_init)
    all_noise = [init_noise] + step_noise_list

    # ------------------------------------------------------------------
    # Phase 1 serial: per-sample metrics
    # ------------------------------------------------------------------
    serial_per_sample = [{} for _ in range(num_val_samples)]

    _original_torch_randn = torch.randn

    def _run_serial_validation_step_batch(batch_i, sample_idx):
        _serial_randn_calls = []
        noise_for_sample = [n[sample_idx : sample_idx + 1].clone() for n in all_noise]
        _serial_randn_sequence = noise_for_sample

        def _fixed_randn(*args, _seq=_serial_randn_sequence, _calls=_serial_randn_calls, **kwargs):
            idx = len(_calls)
            _calls.append(idx)
            if idx < len(_seq):
                return _seq[idx].clone()
            return _original_torch_randn(*args, **kwargs)

        _serial_mp = pytest.MonkeyPatch()
        _serial_mp.setattr(
            serial_diffusion_v2_module,
            "compute_random_augmentation",
            lambda mult, device=None, dtype=None: (
                torch.eye(3, device=device, dtype=dtype).unsqueeze(0).expand(diffusion_samples, -1, -1),
                torch.zeros(diffusion_samples, 1, 3, device=device, dtype=dtype),
            ),
        )
        _serial_mp.setattr(serial_diffusion_v2_module.torch, "randn", _fixed_randn)
        _serial_mp.setattr(reference_module, "log", lambda *a, **kw: None)

        with torch.no_grad():
            reference_module.validation_step(batch_i, batch_idx=sample_idx)

        _serial_mp.undo()

    def _extract_validator_metrics(validator):
        fm = validator.folding_metrics
        val_idx = 0
        sample_metrics = {}
        disto_loss_metric = fm["disto_loss"][val_idx]["disto_loss"]
        sample_metrics["disto_loss"] = disto_loss_metric.compute().item()
        sample_metrics["disto_lddt"] = {}
        sample_metrics["lddt"] = {}
        sample_metrics["complex_lddt"] = {}
        for m_ in [*const.out_types, "pocket_ligand_protein", "contact_protein_protein"]:
            if m_ in fm["disto_lddt"][val_idx]:
                val = fm["disto_lddt"][val_idx][m_].compute()
                if not torch.isnan(val):
                    sample_metrics["disto_lddt"][m_] = val.item()
            if m_ in fm["lddt"][val_idx]:
                val = fm["lddt"][val_idx][m_].compute()
                if not torch.isnan(val):
                    sample_metrics["lddt"][m_] = val.item()
            if m_ in fm["complex_lddt"][val_idx]:
                val = fm["complex_lddt"][val_idx][m_].compute()
                if not torch.isnan(val):
                    sample_metrics["complex_lddt"][m_] = val.item()
        return sample_metrics

    def _reset_validator_metrics(validator):
        fm = validator.folding_metrics
        val_idx = 0
        for metric_group in ["lddt", "disto_lddt", "complex_lddt", "disto_loss"]:
            for k, metric_obj in fm[metric_group][val_idx].items():
                metric_obj.reset()
        # Reset physicalism metrics so the next sample does not accumulate on top
        if getattr(validator, "physicalism_metrics", None) and hasattr(validator.physicalism_metrics, "keys"):
            for group in ["clash", "pb"]:
                if group in validator.physicalism_metrics:
                    for metric_obj in validator.physicalism_metrics[group][val_idx].values():
                        metric_obj.reset()

    for sample_idx in range(B):
        batch_i = _slice_batch(input_feats_global_fp64, sample_idx)
        batch_i["idx_dataset"] = torch.tensor([0], device=device_type)
        _run_serial_validation_step_batch(batch_i, sample_idx)

        serial_per_sample[sample_idx] = _extract_validator_metrics(serial_validators[0])
        _reset_validator_metrics(serial_validators[0])

    # ------------------------------------------------------------------
    # Phase 2 serial: epoch-end metrics (accumulate both samples)
    # ------------------------------------------------------------------
    for vi in range(num_validators):
        for sample_idx in range(B):
            batch_i = _slice_batch(input_feats_global_fp64, sample_idx)
            batch_i["idx_dataset"] = torch.tensor([vi], device=device_type)
            _run_serial_validation_step_batch(batch_i, sample_idx)

    serial_log = _LogCapture(CSVLogger(save_dir=tempfile.mkdtemp(), name="serial_val"))
    _serial_mp2 = pytest.MonkeyPatch()
    _serial_mp2.setattr(reference_module, "log", serial_log)

    reference_module.on_validation_epoch_end()
    _serial_mp2.undo()

    serial_epoch_end_metrics = dict(serial_log.metrics)

    assert len(serial_per_sample[0]) > 0, "Serial phase 1 produced no metrics for sample 0"
    assert len(serial_per_sample[1]) > 0, "Serial phase 1 produced no metrics for sample 1"
    assert len(serial_epoch_end_metrics) > 0, "Serial phase 2 produced no epoch-end metrics"

    # ------------------------------------------------------------------
    # Move to CPU for spawn_multiprocessing
    # ------------------------------------------------------------------
    input_feats_host = {}
    for k, v in input_feats_global_fp64.items():
        if isinstance(v, torch.Tensor):
            input_feats_host[k] = v.detach().to(device="cpu", copy=True)
        elif isinstance(v, list):
            input_feats_host[k] = [
                item.detach().to(device="cpu", copy=True) if isinstance(item, torch.Tensor) else item for item in v
            ]
        else:
            input_feats_host[k] = v
    noise_host_list = [n.detach().cpu() for n in all_noise]
    serial_per_sample_cpu = list(serial_per_sample)

    boltz2_model_params["validators"] = [
        DistributedRCSBValidator(val_names=[vn], confidence_prediction=False, physicalism_metrics=True)
        for vn in val_names
    ]
    spawn_multiprocessing(
        _worker_validation_step_parity,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        boltz2_model_params,
        module_state_dict,
        input_feats_host,
        noise_host_list,
        serial_per_sample_cpu,
        serial_epoch_end_metrics,
        env_per_rank,
    )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
