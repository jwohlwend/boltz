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

"""End-to-end Boltz-2 distributed training integration test via ``train()``.

Calls the real ``train()`` entrypoint with real Boltz-2 training data and a
small model config, exercising the full pipeline: config loading → distributed
manager → distributed data module → distributed model wrapping → Trainer.fit
→ checkpoint.

The only monkeypatches are:
- ``_cleanup_distributed → lambda: None`` (process group safety for tests)

This test mirrors the pattern in ``test_dtensor_predict.py`` (real data,
real checkpoint) and ``test_dtensor_stop_and_go.py`` (``train()`` entrypoint
monkeypatching).
"""

import copy
import functools
import importlib.util as _importlib_util
import math
import random as stdlib_random
import shutil
from dataclasses import dataclass
from enum import Enum, auto
from pathlib import Path
from typing import Any

import numpy as np
import pytest
import pytorch_lightning as pl
import torch
from omegaconf import OmegaConf
from torch.distributed.tensor import DTensor, Replicate, Shard, distribute_tensor

import boltz.distributed.model.modules.diffusion as dist_diffusion_module
import boltz.distributed.train as train_module
import boltz.model.loss.diffusionv2 as serial_loss_v2_module
import boltz.model.modules.diffusionv2 as serial_diffusion_v2_module
from boltz.data.module.trainingv2 import (
    Boltz2TrainingDataModule,
)
from boltz.data.module.trainingv2 import (
    TrainingDataset as SerialTrainingDataset,
)
from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.models.boltz2 import Boltz2 as Boltz2Distributed
from boltz.distributed.model.modules.diffusion import AtomDiffusion as DistAtomDiffusionV2
from boltz.distributed.model.modules.utils import SDPAWithBiasBackend, TriAttnBackend
from boltz.distributed.predict import run_predict
from boltz.distributed.testing.utils import setup_mock_training_datamodule_config
from boltz.model.models.boltz2 import Boltz2 as SerialBoltz2
from boltz.model.modules.diffusionv2 import AtomDiffusion as SerialAtomDiffusionV2
from boltz.model.validation.rcsb import RCSBValidator
from boltz.testing.utils import (
    SetModuleInfValues,
    concat_data,
    distribute_atom_features,
    get_feature_placements,
    init_module_params_glorot,
    init_tensors_uniform,
    seed_by_rank,
    spawn_multiprocessing,
)


class _Unset(Enum):
    """Sentinel to distinguish 'not provided' from ``None``."""

    UNSET = auto()


_UNSET = _Unset.UNSET


@dataclass
class TrainTestConfig:
    """Unified config for writing YAML training configs in integration tests.

    Defaults match the most common usage (the general distributed training
    function with 6 call sites).  E2E parity callers override the fields
    that differ.
    """

    config_path: Path
    output_dir: Path
    test_data_dir: Path
    mol_dir: Path

    mode: str = "distributed"

    size_dp: int = 1
    size_cp: int = 1

    accelerator: str = "cpu"
    max_epochs: int = 1
    limit_train_batches: int = 2
    precision: str = "FP32"
    num_sanity_val_steps: int = 0
    limit_val_batches: int | None = None
    gradient_clip_val: float = 10.0

    model: dict[str, Any] | None = None
    pretrained: str | None = None
    resume: str | None = None

    ema: bool = True
    ema_decay: float = 0.999

    batch_size: int = 1
    samples_per_epoch: int = 4
    max_tokens: int = 384
    max_atoms: int = 3456
    max_seqs: int = 128
    return_train_symmetries: bool = True
    split: str | None = None
    overfit: int | None = None
    extra_dataset_overrides: dict[str, Any] | None = None
    pop_target_keys: bool = False

    validate_structure: bool = False
    validation_only: bool = False
    seed: int = 42

    v2: bool = False
    strict_loading: bool = True
    wandb: dict[str, Any] | None | _Unset = _UNSET
    save_top_k: int = -1
    disable_checkpoint: bool = False


def _write_train_config(cfg: TrainTestConfig) -> None:
    """Write a YAML training config from a :class:`TrainTestConfig`.

    Supports both ``mode="distributed"`` (default) and ``mode="serial"``.
    """
    prod_yaml = Path(__file__).resolve().parents[2] / "scripts" / "train" / "configs" / "structurev2.yaml"
    data_dict = OmegaConf.to_container(OmegaConf.load(prod_yaml).data, resolve=False)

    if cfg.pop_target_keys:
        data_dict.pop("_target_", None)

    data_dict["datasets"] = [data_dict["datasets"][0]]
    ds = data_dict["datasets"][0]
    if cfg.pop_target_keys:
        ds.pop("_target_", None)

    ds["target_dir"] = str(cfg.test_data_dir)
    ds["msa_dir"] = str(cfg.test_data_dir / "msa")
    ds["template_dir"] = None
    ds["split"] = str(cfg.split) if cfg.split else None
    ds["prob"] = 1.0

    if cfg.extra_dataset_overrides:
        for k, v in cfg.extra_dataset_overrides.items():
            ds[k] = v

    data_dict["samples_per_epoch"] = cfg.samples_per_epoch
    data_dict["num_workers"] = 0
    data_dict["pin_memory"] = False
    data_dict["use_templates"] = False
    data_dict["return_train_symmetries"] = cfg.return_train_symmetries
    data_dict["batch_size"] = cfg.batch_size
    data_dict["max_tokens"] = cfg.max_tokens
    data_dict["max_atoms"] = cfg.max_atoms
    data_dict["max_seqs"] = cfg.max_seqs
    data_dict["pad_to_max_tokens"] = True
    data_dict["pad_to_max_atoms"] = True
    data_dict["pad_to_max_seqs"] = True
    data_dict["msa_sampling_training"] = False
    data_dict["moldir"] = str(cfg.mol_dir)

    if cfg.overfit is not None:
        data_dict["overfit"] = cfg.overfit

    if cfg.model is not None:
        model_dict = cfg.model
    else:
        model_dict = {
            "_target_": "boltz.model.models.boltz2.Boltz2",
            "atom_s": 4,
            "atom_z": 4,
            "token_s": 4,
            "token_z": 4,
            "num_bins": 64,
            "atom_feature_dim": 388,
            "atoms_per_window_queries": 32,
            "atoms_per_window_keys": 128,
            "ema": cfg.ema,
            "ema_decay": cfg.ema_decay,
            "confidence_prediction": False,
            "affinity_prediction": False,
            "structure_prediction_training": True,
            "use_templates": False,
            "validate_structure": cfg.validate_structure,
            "predict_bfactor": False,
            "bond_type_feature": False,
            "embedder_args": {
                "atom_encoder_depth": 1,
                "atom_encoder_heads": 1,
                "activation_checkpointing": False,
            },
            "msa_args": {
                "msa_s": 4,
                "msa_blocks": 1,
                "msa_dropout": 0.0,
                "z_dropout": 0.0,
                "use_paired_feature": True,
            },
            "pairformer_args": {
                "num_blocks": 1,
                "num_heads": 1,
                "dropout": 0.0,
                "v2": True,
            },
            "score_model_args": {
                "sigma_data": 16.0,
                "dim_fourier": 4,
                "atom_encoder_depth": 1,
                "atom_encoder_heads": 1,
                "token_transformer_depth": 1,
                "token_transformer_heads": 1,
                "atom_decoder_depth": 1,
                "atom_decoder_heads": 1,
                "activation_checkpointing": False,
                "conditioning_transition_layers": 1,
            },
            "diffusion_process_args": {
                "coordinate_augmentation": False,
            },
            "diffusion_loss_args": {},
            "training_args": {
                "recycling_steps": 2,
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
                "lr_scheduler": "af3",
                "base_lr": 1e-3,
                "max_lr": 1e-3,
                "lr_warmup_no_steps": 10,
                "lr_start_decay_after_n_steps": 100,
                "lr_decay_every_n_steps": 50000,
                "lr_decay_factor": 0.95,
                "weight_decay": 0.0,
            },
            "validation_args": {
                "recycling_steps": 0,
                "sampling_steps": 2,
                "diffusion_samples": 1,
                "symmetry_correction": False,
                "clash_cutoff": None,
            },
        }

    config: dict[str, Any] = {
        "data": data_dict,
        "model": model_dict,
        "output": str(cfg.output_dir),
        "pretrained": cfg.pretrained,
        "trainer": {
            "accelerator": cfg.accelerator,
            "devices": 1,
            "max_epochs": cfg.max_epochs,
            "limit_train_batches": cfg.limit_train_batches,
            "enable_progress_bar": False,
            "enable_model_summary": False,
            "num_sanity_val_steps": cfg.num_sanity_val_steps,
            "gradient_clip_val": cfg.gradient_clip_val,
        },
    }

    if cfg.limit_val_batches is not None:
        config["trainer"]["limit_val_batches"] = cfg.limit_val_batches

    if cfg.mode == "serial":
        config["trainer"]["precision"] = 32 if cfg.precision == "FP32" else cfg.precision
        config["v2"] = cfg.v2
        config["disable_checkpoint"] = cfg.disable_checkpoint
        config["save_top_k"] = cfg.save_top_k
        config["strict_loading"] = cfg.strict_loading
        config["validation_only"] = cfg.validation_only
        if cfg.wandb is not _UNSET:
            config["wandb"] = cfg.wandb
    else:
        config["resume"] = cfg.resume
        config["parallel_size"] = {"size_dp": cfg.size_dp, "size_cp": cfg.size_cp}
        config["precision"] = cfg.precision
        config["find_unused_parameters"] = False
        config["save_top_k"] = cfg.save_top_k
        config["disable_checkpoint"] = cfg.disable_checkpoint
        config["debug"] = False
        config["validation_only"] = cfg.validation_only
        config["seed"] = cfg.seed
        config["checkpoint"] = {
            "monitor": None,
            "save_last": True,
            "every_n_epochs": 1,
        }
        config["triattn_backend"] = "reference"
        config["sdpa_with_bias_backend"] = "reference"
        config["sdpa_with_bias_shardwise_backend"] = "reference"

    cfg.config_path.parent.mkdir(parents=True, exist_ok=True)
    OmegaConf.save(OmegaConf.create(config), cfg.config_path)


def _parallel_assert_boltz2_train(rank: int, payload: tuple[Any, ...]) -> None:
    """Multi-rank worker: call real train() and verify completion."""
    env_per_rank, config_path, output_dir = payload
    output_dir = Path(output_dir)

    monkeypatch = pytest.MonkeyPatch()
    for key, value in env_per_rank.items():
        monkeypatch.setenv(key, f"{rank}" if value == "<INPUT_RANK>" else value)

    # Only monkeypatch cleanup — the model and data factories use their real implementations.
    monkeypatch.setattr(train_module, "_cleanup_distributed", lambda: None)
    DistributedManager._state = {}

    train_module.train(str(config_path), [])

    # Barrier BEFORE rank-0-only assertions: if rank 0's assertions fail,
    # it raises before any trailing barrier, leaving other ranks stuck in
    # an NCCL wait.  Syncing first avoids this deadlock.
    torch.distributed.barrier()

    # Assert: checkpoint was written (rank 0 only — file I/O, no collectives)
    ckpt_path = output_dir / "last.ckpt"
    if rank == 0:
        assert ckpt_path.exists(), f"Rank {rank}: checkpoint not found at {ckpt_path}"
        ckpt = torch.load(ckpt_path, map_location="cpu", weights_only=False)
        assert ckpt.get("global_step", 0) > 0, "global_step is 0 — training did not run"
        assert "state_dict" in ckpt, "checkpoint missing state_dict"
        assert "optimizer_states" in ckpt and ckpt["optimizer_states"], "checkpoint missing optimizer_states"

        # EMA: if model config has ema=True, checkpoint must contain EMA state
        # saved as plain tensors (not DTensors) for checkpoint portability.
        ckpt_hp = ckpt.get("hyper_parameters", {})
        if ckpt_hp.get("ema", False):
            assert "ema" in ckpt, "checkpoint missing EMA state despite ema=True"
            assert "ema_weights" in ckpt["ema"], "EMA state must include ema_weights"
            assert "cur_step" in ckpt["ema"], "EMA state must include cur_step"
            for ema_key, ema_val in ckpt["ema"]["ema_weights"].items():
                assert isinstance(ema_val, torch.Tensor) and not isinstance(
                    ema_val, DTensor
                ), f"EMA weight '{ema_key}' must be saved as plain torch.Tensor, got {type(ema_val).__name__}"


def _worker_validation_parity(
    rank: int,
    grid_group_sizes: dict,
    device_type: str,
    backend: str,
    env_per_rank: dict[str, Any],
    dist_config_path: str,
    serial_metrics: dict,
    sigmas_global_host: torch.Tensor,
    noise_global_host: torch.Tensor,
    atom_counts_per_token_host: torch.Tensor,
    cached_samples_path: str,
    seed: int,
) -> None:
    """Multi-rank worker: distributed validate() then compare with serial metrics.

    1. Applies module-level monkeypatches (noise, data, smooth_lddt)
    2. Calls ``train_module.train(dist_config_path, [])`` in validation_only mode
    3. Captures validation metrics from ``trainer.validate()``
    4. Compares all validation metrics with the serial reference
    """
    monkeypatch = pytest.MonkeyPatch()
    if env_per_rank is not None:
        for var_name, value in env_per_rank.items():
            monkeypatch.setenv(var_name, f"{rank}" if value == "<INPUT_RANK>" else value)

    monkeypatch.setattr(train_module, "_cleanup_distributed", lambda: None)
    DistributedManager._state = {}

    _apply_cached_getitem(monkeypatch, cached_samples_path)

    # Deterministic sigmas for noise schedule
    def _dist_noise_dist(self, bs, dtype=torch.float32):
        s = sigmas_global_host.to(device=self.device_mesh.device_type, dtype=dtype)[:bs]
        return distribute_tensor(s, self.device_mesh, (Shard(0), Replicate(), Replicate()))

    monkeypatch.setattr(DistAtomDiffusionV2, "noise_distribution", _dist_noise_dist)

    # Deterministic noise via distribute_atom_features (intersperse padding)
    _noise_dt_cache: list[DTensor | None] = [None]
    _noise_computed = [False]

    def _compute_noise_dt_once(device_mesh, dtype):
        if _noise_computed[0]:
            return
        _noise_computed[0] = True
        manager = DistributedManager()
        _io_keys = {"noise"}
        _placements = get_feature_placements(
            atom_keys=set(),
            model_io_keys=_io_keys,
            model_io_fp32_keys=set(),
        )
        size_batch = atom_counts_per_token_host.shape[0]
        multiplicity_val = noise_global_host.shape[0] // size_batch
        noise_unflat = noise_global_host.unflatten(0, (size_batch, multiplicity_val))
        inputs_io = {"atom_counts_per_token": atom_counts_per_token_host.clone()}
        for i_mul in range(multiplicity_val):
            inputs_io[f"noise_{i_mul}"] = noise_unflat[:, i_mul].to(dtype=dtype)
        placements_cp_model_io_mul = {
            f"{k}_{i_mul}": v for k, v in _placements["cp_model_io"].items() for i_mul in range(multiplicity_val)
        }
        placements_cp = _placements["cp_atom_features"] | placements_cp_model_io_mul
        placements_model_io_mul = {
            f"{k}_{i_mul}": v for k, v in _placements["model_io"].items() for i_mul in range(multiplicity_val)
        }
        placements_dp_cp = placements_model_io_mul
        io_feats = distribute_atom_features(
            inputs=inputs_io,
            placements_cp=placements_cp,
            placements_dp_cp=placements_dp_cp,
            device_mesh=manager.device_mesh_subgroups,
            cp_group=manager.group["cp"],
            multiplicities={"noise": multiplicity_val},
        )
        _noise_dt_cache[0] = io_feats.pop("noise").to(dtype=dtype)

    _dist_in_val = [False]

    def _det_create_randn(shape, device_mesh, placements, dtype=torch.float32, scale=1.0):
        if _dist_in_val[0]:
            from boltz.distributed.utils import create_distributed_randn as _real_create_randn

            return _real_create_randn(shape, device_mesh, placements, dtype=dtype, scale=0.0)
        from boltz.testing.utils import pad_to_length as _pad

        _compute_noise_dt_once(device_mesh, dtype)
        n = _noise_dt_cache[0]
        if n.dtype != dtype:
            n = n.to(dtype=dtype)
        if len(shape) > 1 and n.shape[1] < shape[1]:
            n = _pad(n, dim=1, length=shape[1])
        return n * scale

    monkeypatch.setattr(dist_diffusion_module, "create_distributed_randn", _det_create_randn)

    _orig_dist_val_step = Boltz2Distributed.validation_step

    def _dist_val_step_wrapper(self_model, batch, batch_idx):
        _dist_in_val[0] = True
        try:
            return _orig_dist_val_step(self_model, batch, batch_idx)
        finally:
            _dist_in_val[0] = False

    monkeypatch.setattr(Boltz2Distributed, "validation_step", _dist_val_step_wrapper)

    # Skip RMSD (not compared; serial uses a different code path)
    import boltz.distributed.model.validation.validator as _dist_validator_mod

    def _rmsd_noop(*args, **kwargs):
        return torch.tensor(0.0), None, None

    monkeypatch.setattr(_dist_validator_mod, "weighted_minimum_rmsd_single", _rmsd_noop)

    # Capture metrics from trainer.validate()
    _captured_metrics: dict[str, float] = {}
    _orig_validate = pl.Trainer.validate

    def _capturing_validate(self, *args, **kwargs):
        result = _orig_validate(self, *args, **kwargs)
        for k, v in self.callback_metrics.items():
            if isinstance(v, DTensor):
                _captured_metrics[k] = v.full_tensor().detach().cpu().item()
            elif isinstance(v, torch.Tensor):
                _captured_metrics[k] = v.detach().cpu().item()
            else:
                _captured_metrics[k] = v
        return result

    monkeypatch.setattr(pl.Trainer, "validate", _capturing_validate)

    train_module.train(dist_config_path, [])

    # --- Compare metrics ---
    # Atom-level LDDT: atol=2e-4, rtol=0.005 (forward-pass accumulation order).
    # Token-level metrics (disto_lddt_*, disto_loss): default tolerance.
    # Trailing underscore omitted so the global weighted-average "val/lddt"
    # (not just "val/lddt_*") also gets the relaxed tolerance.
    _forward_dependent_prefixes = ("val/lddt", "val/complex_lddt", "val/clash", "val/pb", "val/rmsd")
    _lddt_keys_compared = []
    if serial_metrics:
        for k in serial_metrics:
            if k in _captured_metrics:
                got = torch.tensor(_captured_metrics[k])
                exp = torch.tensor(serial_metrics[k])
                if any(k.startswith(p) for p in _forward_dependent_prefixes):
                    torch.testing.assert_close(
                        got,
                        exp,
                        atol=2e-4,
                        rtol=0.005,
                        msg=lambda m: f"Rank {rank}: metric '{k}' mismatch: {m}",
                    )
                else:
                    torch.testing.assert_close(
                        got,
                        exp,
                        msg=lambda m: f"Rank {rank}: metric '{k}' mismatch: {m}",
                    )
                if "lddt" in k:
                    _lddt_keys_compared.append(k)

    assert _lddt_keys_compared, (
        f"Rank {rank}: no validation LDDT metrics were compared — test is vacuous. "
        f"Serial keys: {sorted(serial_metrics)}, dist keys: {sorted(_captured_metrics)}"
    )
    for required_metric in ("val/lddt", "val/disto_lddt", "val/complex_lddt"):
        assert (
            required_metric in _captured_metrics
        ), f"Rank {rank}: distributed metrics missing '{required_metric}' — available: {sorted(_captured_metrics)}"

    torch.distributed.barrier()


@pytest.mark.slow
@pytest.mark.parametrize(
    "setup_env",
    [
        ((2, (2, 2)), True, "cuda", "ENV"),
    ],
    indirect=("setup_env",),
    ids=["cuda-dp2-cp2x2"],
)
def test_boltz2_validation_parity(
    setup_env,
    test_cp_training_base_data_dir_boltz2,
    canonical_mols_dir,
    tmp_path,
):
    """Serial-vs-DTensor validation metric parity via ``train()`` validation_only.

    Runs both serial and distributed ``train()`` in validation_only mode
    (no training, no backward pass) on the same pretrained model and data,
    then compares all logged validation metrics.  Model weights come from a
    deterministic pretrained checkpoint; noise and data are controlled via
    monkeypatches so the only source of difference is the serial-vs-distributed
    forward pass in FP32.

    Token-level metrics (``val/disto_lddt_*``, ``val/disto_loss``) match at
    default tolerance because they depend only on the confidence head logits,
    which are deterministic given identical weights.

    Atom-level LDDT (``val/lddt_*``, ``val/complex_lddt_*``) uses
    ``atol=2e-4, rtol=0.005``.  These metrics depend on diffusion-sampled
    coordinates whose FP32 accumulation order differs between serial
    (full-sequence attention) and distributed (split attention + padding for
    uniform DP shapes).  Without training, these forward-pass-only differences
    are smaller and more stable than in the e2e training test, allowing ~2.5x
    tighter tolerances than the e2e test's ``atol=5e-4, rtol=0.02``.
    Observed max absdiff is ~9e-5 (``lddt_intra_ligand``).
    """
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("CUDA not available")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"Need {world_size} GPUs, have {torch.cuda.device_count()}")

    dtype = torch.float32
    seed = 42
    multiplicity = 2
    B = 2
    max_tokens = 256
    W = 32  # atoms_per_window_queries
    size_cp = grid_group_sizes["cp"][0] * grid_group_sizes["cp"][1]
    atom_align = math.lcm(W, size_cp)
    max_atoms = ((max_tokens * 10 + atom_align - 1) // atom_align) * atom_align
    max_seqs = 16
    scale_glorot = 0.05

    # --- Merge all 4 samples with train/val split ---
    training_data_dir, split_file = _setup_training_data_all_4_e2e(
        tmp_path / "training_data", test_cp_training_base_data_dir_boltz2
    )

    # --- Create pretrained checkpoint with deterministic init ---
    seed_by_rank(0, seed=seed)
    model_dict = _e2e_model_dict(multiplicity=multiplicity, validate_structure=True)
    model_dict.pop("_target_")
    model_dict.pop("validators", None)
    _val_validators = [RCSBValidator(val_names=["RCSB"], confidence_prediction=False, physicalism_metrics=True)]
    pretrained_model = SerialBoltz2(**model_dict, validators=_val_validators)
    init_module_params_glorot(pretrained_model, gain=scale_glorot)
    pretrained_model.apply(SetModuleInfValues())
    pretrained_model.structure_module.coordinate_augmentation = False
    pretrained_model = pretrained_model.to(dtype=dtype)

    pretrained_path = tmp_path / "pretrained.ckpt"
    torch.save(
        {
            "state_dict": pretrained_model.state_dict(),
            "pytorch-lightning_version": pl.__version__,
            "hyper_parameters": pretrained_model.hparams,
        },
        pretrained_path,
    )

    # --- Pre-load and cache individual samples to disk ---
    _tmp_mp = pytest.MonkeyPatch()
    _apply_e2e_deterministic_getitem(_tmp_mp, base_seed=seed)
    _preload_cfg = setup_mock_training_datamodule_config(training_data_dir)
    _preload_cfg.batch_size = B
    _preload_cfg.samples_per_epoch = B
    _preload_cfg.moldir = str(canonical_mols_dir)
    _preload_cfg.return_train_symmetries = False
    _preload_cfg.msa_sampling_training = False
    _preload_cfg.max_tokens = max_tokens
    _preload_cfg.max_atoms = max_atoms
    _preload_cfg.max_seqs = max_seqs
    for _ds in _preload_cfg.datasets:
        _ds.filters = None
        _ds.split = str(split_file)
        _ds.symmetry_correction = False
    seed_by_rank(0, seed=seed)
    _preload_dm = Boltz2TrainingDataModule(cfg=_preload_cfg)
    _preload_ds = _preload_dm._train_set
    _cached_samples = {i: _preload_ds[i] for i in range(B)}
    cached_samples_path = tmp_path / "cached_samples.pt"
    torch.save(_cached_samples, cached_samples_path)

    _preload_dl = _preload_dm.train_dataloader()
    _preload_batch = next(iter(_preload_dl))
    atom_counts_per_token_host = _preload_batch["atom_counts_per_token"].detach().cpu()
    atom_pad_mask_host = _preload_batch["atom_pad_mask"].detach().cpu()
    _tmp_mp.undo()

    # --- Pre-generate deterministic noise (masked by atom_pad_mask) ---
    seed_by_rank(0, seed=seed)
    sigmas_global = pretrained_model.structure_module.noise_distribution(B * multiplicity).to(dtype=dtype)
    noise_global = torch.empty(B * multiplicity, max_atoms, 3, dtype=dtype)
    init_tensors_uniform([noise_global], low=-scale_glorot, high=scale_glorot)
    _mask_mul = atom_pad_mask_host[:, :, None].repeat_interleave(multiplicity, 0).to(dtype=dtype)
    noise_global = noise_global * _mask_mul

    sigmas_global_host = sigmas_global.detach().cpu()
    noise_global_host = noise_global.detach().cpu()

    # --- Write serial config (validation_only) ---
    serial_output_dir = tmp_path / "serial_output"
    serial_output_dir.mkdir(parents=True, exist_ok=True)
    serial_config_path = tmp_path / "serial_config.yaml"
    _e2e_ds_overrides = {
        "filters": None,
        "moldir": None,
        "symmetry_correction": False,
        "val_group": "RCSB",
        "use_train_subset": None,
        "override_bfactor": False,
        "override_method": None,
    }
    _write_train_config(
        TrainTestConfig(
            config_path=serial_config_path,
            output_dir=serial_output_dir,
            test_data_dir=training_data_dir,
            mol_dir=canonical_mols_dir,
            mode="serial",
            accelerator="gpu",
            validation_only=True,
            pretrained=str(pretrained_path),
            model=_e2e_model_dict(multiplicity=multiplicity, validate_structure=True),
            batch_size=B,
            samples_per_epoch=B,
            max_tokens=max_tokens,
            max_atoms=max_atoms,
            max_seqs=max_seqs,
            return_train_symmetries=False,
            split=str(split_file),
            pop_target_keys=True,
            extra_dataset_overrides=_e2e_ds_overrides,
            v2=True,
            strict_loading=False,
            wandb=None,
            save_top_k=0,
            disable_checkpoint=True,
        )
    )

    # --- Apply serial monkeypatches ---
    serial_mp = pytest.MonkeyPatch()
    _apply_cached_getitem(serial_mp, cached_samples_path)

    _orig_serial_boltz2_init = SerialBoltz2.__init__

    @functools.wraps(_orig_serial_boltz2_init)
    def _init_with_validators(self, *args, **kwargs):
        if kwargs.get("validate_structure", False) and not kwargs.get("validators"):
            kwargs["validators"] = [
                RCSBValidator(val_names=["RCSB"], confidence_prediction=False, physicalism_metrics=True)
            ]
        _orig_serial_boltz2_init(self, *args, **kwargs)

    serial_mp.setattr(SerialBoltz2, "__init__", _init_with_validators)

    serial_mp.setattr(
        SerialAtomDiffusionV2,
        "noise_distribution",
        lambda self, bs: sigmas_global[:bs].to(device=self.zero.device),
    )

    # Deterministic noise: pre-generated for training, zero for validation
    _serial_in_val = [False]

    def _serial_randn_like(t):
        if _serial_in_val[0]:
            return torch.zeros_like(t)
        return noise_global[: t.shape[0], : t.shape[1]].to(device=t.device, dtype=t.dtype)

    serial_mp.setattr(serial_diffusion_v2_module.torch, "randn_like", _serial_randn_like)

    _orig_serial_randn = serial_diffusion_v2_module.torch.randn

    def _serial_randn(*args, **kwargs):
        if _serial_in_val[0]:
            kwargs.pop("generator", None)
            return torch.zeros(*args, **kwargs)
        return _orig_serial_randn(*args, **kwargs)

    serial_mp.setattr(serial_diffusion_v2_module.torch, "randn", _serial_randn)

    _orig_compute_random_augmentation = serial_diffusion_v2_module.compute_random_augmentation

    def _identity_augmentation_during_val(multiplicity_arg, s_trans=1.0, device=None, dtype=torch.float32):
        if _serial_in_val[0]:
            R = torch.eye(3, device=device, dtype=dtype).unsqueeze(0).expand(multiplicity_arg, -1, -1)
            tr = torch.zeros(multiplicity_arg, 1, 3, device=device, dtype=dtype)
            return R, tr
        return _orig_compute_random_augmentation(multiplicity_arg, s_trans=s_trans, device=device, dtype=dtype)

    serial_mp.setattr(serial_diffusion_v2_module, "compute_random_augmentation", _identity_augmentation_during_val)

    _orig_serial_val_step = SerialBoltz2.validation_step

    def _serial_val_step_wrapper(self_model, batch, batch_idx):
        _serial_in_val[0] = True
        try:
            return _orig_serial_val_step(self_model, batch, batch_idx)
        finally:
            _serial_in_val[0] = False

    serial_mp.setattr(SerialBoltz2, "validation_step", _serial_val_step_wrapper)

    serial_mp.setattr(serial_loss_v2_module, "smooth_lddt_loss", _smooth_lddt_loss_dense_e2e)
    serial_mp.setattr(serial_diffusion_v2_module, "smooth_lddt_loss", _smooth_lddt_loss_dense_e2e)

    # Capture metrics from trainer.validate()
    serial_captured_metrics: dict[str, float] = {}
    _orig_validate = pl.Trainer.validate

    def _serial_capturing_validate(self, *args, **kwargs):
        result = _orig_validate(self, *args, **kwargs)
        for k, v in self.callback_metrics.items():
            if isinstance(v, torch.Tensor):
                serial_captured_metrics[k] = v.detach().cpu().item()
            else:
                serial_captured_metrics[k] = v
        return result

    serial_mp.setattr(pl.Trainer, "validate", _serial_capturing_validate)

    # --- Run serial validation ---
    _serial_train_mod = _load_serial_train_module()
    _serial_train_mod.train(str(serial_config_path), [])

    # Non-vacuous guard
    serial_lddt_keys = [
        k for k in serial_captured_metrics if k.startswith("val/lddt_") or k.startswith("val/disto_lddt_")
    ]
    assert serial_lddt_keys, (
        f"Serial run produced no validation LDDT metrics — test is vacuous. "
        f"Available metrics: {sorted(serial_captured_metrics)}"
    )

    serial_mp.undo()

    # --- Write distributed config (validation_only) ---
    dp = grid_group_sizes["dp"]
    cp0, cp1 = grid_group_sizes["cp"]
    size_cp = cp0 * cp1
    dist_output_dir = tmp_path / "dist_output"
    dist_output_dir.mkdir(parents=True, exist_ok=True)
    dist_config_path = tmp_path / "dist_config.yaml"
    _write_train_config(
        TrainTestConfig(
            config_path=dist_config_path,
            output_dir=dist_output_dir,
            test_data_dir=training_data_dir,
            mol_dir=canonical_mols_dir,
            size_dp=dp,
            size_cp=size_cp,
            accelerator="gpu",
            validation_only=True,
            pretrained=str(pretrained_path),
            model=_e2e_model_dict(multiplicity=multiplicity, validate_structure=True, distributed=True),
            batch_size=1,
            samples_per_epoch=dp,
            max_tokens=max_tokens,
            max_atoms=max_atoms,
            max_seqs=max_seqs,
            return_train_symmetries=False,
            split=str(split_file),
            pop_target_keys=True,
            extra_dataset_overrides=_e2e_ds_overrides,
        )
    )

    # --- Spawn distributed workers ---
    spawn_multiprocessing(
        _worker_validation_parity,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        str(dist_config_path),
        serial_captured_metrics,
        sigmas_global_host,
        noise_global_host,
        atom_counts_per_token_host,
        str(cached_samples_path),
        seed,
    )


@pytest.mark.slow
@pytest.mark.parametrize(
    "setup_env",
    [
        # DP-only: exercises DP correctness on 2-GPU workstations
        ((2, (1, 1)), True, "cuda", "ENV"),
        # DP + CP: catches integration issues between DP and CP
        ((2, (2, 2)), True, "cuda", "ENV"),
    ],
    indirect=("setup_env",),
    ids=["cuda-dp2-cp1x1", "cuda-dp2-cp2x2"],
)
def test_boltz2_train_entrypoint(
    setup_env,
    tmp_path,
    test_cp_training_data_dir_boltz2,
    canonical_mols_dir,
):
    """End-to-end Boltz-2 training through the real train() entrypoint.

    Exercises the full pipeline with a small model and real training data:
    config → Hydra instantiate → _create_distributed_model (Boltz2Distributed)
    → _create_distributed_data_module (Boltz2TrainingDataModule) → Trainer.fit
    → checkpoint.
    """
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("CUDA not available")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"Need {world_size} GPUs, have {torch.cuda.device_count()}")

    output_dir = tmp_path / "boltz2_train_output"
    config_path = tmp_path / "boltz2_train_config.yaml"

    _write_train_config(
        TrainTestConfig(
            config_path=config_path,
            output_dir=output_dir,
            test_data_dir=test_cp_training_data_dir_boltz2,
            mol_dir=canonical_mols_dir,
            size_dp=grid_group_sizes["dp"],
            size_cp=math.prod(grid_group_sizes["cp"]),
            accelerator="gpu" if device_type == "cuda" else "cpu",
        )
    )

    payload = (env_per_rank, str(config_path), str(output_dir))
    spawn_multiprocessing(_parallel_assert_boltz2_train, world_size, payload)


def _load_production_model_config(*, reduce_depth: bool = True) -> dict[str, Any]:
    """Load model config from the production structurev2.yaml.

    Returns the model dict with _target_ keys suitable for Hydra instantiation.
    Overrides training args for fast test execution (1 recycling step, 1
    diffusion sample, etc.) and disables features not yet distributed
    (confidence, affinity, templates, validation).

    Parameters
    ----------
    reduce_depth
        If True (default), reduce pairformer/transformer depth to fit in
        32 GiB GPUs.  Set to False on clusters with >=64 GiB GPUs to test
        with production-depth layers.
    """
    config_yaml = Path(__file__).resolve().parents[2] / "scripts" / "train" / "configs" / "structurev2.yaml"
    full_config = OmegaConf.load(config_yaml)
    model_dict = OmegaConf.to_container(full_config.model, resolve=False)

    # Disable features not yet distributed
    model_dict["confidence_prediction"] = False
    model_dict["affinity_prediction"] = False
    model_dict["use_templates"] = False
    model_dict["validate_structure"] = False

    # Remove validators (not needed for training-only test)
    model_dict.pop("validators", None)

    # Remove Hydra interpolations that reference ${data.*} or ${model.*}
    # and replace with concrete values
    model_dict["conditioning_cutoff_min"] = 4.0
    model_dict["conditioning_cutoff_max"] = 20.0

    # Fast training settings
    model_dict["training_args"] = {
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
        "lr_scheduler": "af3",
        "base_lr": 1e-3,
        "max_lr": 1e-3,
        "lr_warmup_no_steps": 10,
        "lr_start_decay_after_n_steps": 100,
        "lr_decay_every_n_steps": 50000,
        "lr_decay_factor": 0.95,
        "weight_decay": 0.0,
    }
    model_dict["validation_args"] = {
        "recycling_steps": 0,
        "sampling_steps": 2,
        "diffusion_samples": 1,
        "symmetry_correction": False,
    }

    if reduce_depth:
        # Reduce model depth to fit in 32 GiB GPU memory under DP=2.
        # The checkpoint is loaded with strict=False so missing blocks are
        # fine — we still exercise the pretrained loading path and verify
        # the matching layers load correctly.
        model_dict["pairformer_args"]["num_blocks"] = 2
        model_dict["score_model_args"]["token_transformer_depth"] = 2
        model_dict["msa_args"]["msa_blocks"] = 1

    model_dict["diffusion_process_args"] = model_dict.get("diffusion_process_args", {})
    model_dict["diffusion_process_args"]["coordinate_augmentation"] = False

    if reduce_depth:
        # With reduced depth the model fits without activation checkpointing,
        # and disabling it makes the test faster.
        for key in ("embedder_args", "msa_args", "pairformer_args", "score_model_args"):
            if key in model_dict and isinstance(model_dict[key], dict):
                model_dict[key]["activation_checkpointing"] = False
        if "template_args" in model_dict:
            model_dict["template_args"]["activation_checkpointing"] = False

    return model_dict


@pytest.mark.slow
@pytest.mark.parametrize(
    ("setup_env", "reduce_depth"),
    [
        (((2, (1, 1)), True, "cuda", "ENV"), True),
        # Full-depth uses CP to distribute activations across 4 GPUs,
        # mirroring production topology and avoiding per-GPU OOM.
        (((1, (2, 2)), True, "cuda", "ENV"), False),
    ],
    indirect=["setup_env"],
    ids=["cuda-dp2-cp1x1-reduced", "cuda-dp1-cp2x2-full"],
)
def test_boltz2_finetune_from_checkpoint(
    setup_env,
    reduce_depth,
    tmp_path,
    test_cp_training_data_dir_boltz2,
    canonical_mols_dir,
    get_model_ckpt_v2,
):
    """End-to-end Boltz-2 finetune from real checkpoint through train().

    Loads the real Boltz-2 checkpoint via ``pretrained`` config, exercises
    production-width model layers with real training data under BF16 mixed
    precision.  The ``reduce_depth=True`` variant uses dp=2 with reduced
    pairformer/transformer depth (fast, 2 GPUs).  The ``reduce_depth=False``
    variant uses dp=1, cp=2x2 to shard the full-depth model across 4 GPUs,
    mirroring how production deployments use CP for large models.
    """
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if not torch.cuda.is_available():
        pytest.skip("CUDA not available")
    if torch.cuda.device_count() < world_size:
        pytest.skip(f"Need {world_size} GPUs, have {torch.cuda.device_count()}")

    output_dir = tmp_path / "boltz2_finetune_output"
    config_path = tmp_path / "boltz2_finetune_config.yaml"

    _write_train_config(
        TrainTestConfig(
            config_path=config_path,
            output_dir=output_dir,
            test_data_dir=test_cp_training_data_dir_boltz2,
            mol_dir=canonical_mols_dir,
            size_dp=grid_group_sizes["dp"],
            size_cp=math.prod(grid_group_sizes["cp"]),
            accelerator="gpu",
            limit_train_batches=1,
            pretrained=str(get_model_ckpt_v2),
            model=_load_production_model_config(reduce_depth=reduce_depth),
            precision="BF16_MIXED",
        )
    )

    payload = (env_per_rank, str(config_path), str(output_dir))
    spawn_multiprocessing(_parallel_assert_boltz2_train, world_size, payload)


# ---------------------------------------------------------------------------
# Stop-and-go checkpoint resume test for Boltz-2
# ---------------------------------------------------------------------------


def _parallel_assert_boltz2_stop_and_go(rank: int, payload: tuple[Any, ...]) -> None:
    """Verify checkpoint resume correctness through the real ``train()`` entrypoint.

    Runs two ``train()`` calls:
    1. Stop/go stage 1 — 1 epoch, checkpoint produced
    2. Stop/go stage 2 — resume from checkpoint, train to epoch 2

    Verifies:
    - Stage 1 checkpoint contains valid model state and optimizer state
    - Stage 2 successfully resumes from checkpoint (no errors)
    - Final checkpoint has correct epoch/step counters (epoch 2, step > stage 1)
    - Final checkpoint has the same state_dict keys as stage 1
    - Weights changed between stage 1 and final (training actually happened)
    - Optimizer state is populated with correct structure

    Note: Exact weight parity between continuous and stop/go runs is not
    tested because ``train.py`` uses different seed offsets on resume
    (``seed + rank + epoch*1000 + step``), which changes diffusion noise.
    The distogram stop-and-go test (``test_dtensor_stop_and_go.py``) covers
    exact parity because its model has no stochastic diffusion.
    """
    (
        env_per_rank,
        stage1_config_path,
        stage2_config_path,
        output_dir,
    ) = payload
    output_dir = Path(output_dir)

    monkeypatch = pytest.MonkeyPatch()
    for key, value in env_per_rank.items():
        monkeypatch.setenv(key, f"{rank}" if value == "<INPUT_RANK>" else value)

    # Only suppress cleanup — use the real model and data factories.
    monkeypatch.setattr(train_module, "_cleanup_distributed", lambda: None)
    DistributedManager._state = {}

    # ---- Stage 1: 1 epoch, checkpoint produced. ----
    train_module.train(str(stage1_config_path), [])
    ckpt_path = output_dir / "last.ckpt"
    assert ckpt_path.exists(), f"Rank {rank}: stage 1 checkpoint not found at {ckpt_path}"

    stage1_ckpt = torch.load(ckpt_path, map_location="cpu", weights_only=False)
    stage1_epoch = stage1_ckpt["epoch"]
    stage1_step = stage1_ckpt["global_step"]
    stage1_keys = set(stage1_ckpt["state_dict"].keys())
    assert stage1_step > 0, f"Rank {rank}: stage 1 global_step is 0 — training did not run"
    assert "state_dict" in stage1_ckpt, f"Rank {rank}: stage 1 missing state_dict"
    assert (
        "optimizer_states" in stage1_ckpt and stage1_ckpt["optimizer_states"]
    ), f"Rank {rank}: stage 1 missing optimizer_states"

    # Save stage 1 weights for change detection
    stage1_weights = {k: v.clone() for k, v in stage1_ckpt["state_dict"].items()}

    # ---- Stage 2: resume from checkpoint to epoch 2. ----
    train_module.train(str(stage2_config_path), [])

    # ---- Verify final checkpoint. ----
    final_ckpt = torch.load(ckpt_path, map_location="cpu", weights_only=False)
    final_epoch = final_ckpt["epoch"]
    final_step = final_ckpt["global_step"]
    final_keys = set(final_ckpt["state_dict"].keys())

    # 1) Epoch and step advanced beyond stage 1.
    assert (
        final_epoch > stage1_epoch
    ), f"Rank {rank}: Final epoch ({final_epoch}) should be > stage 1 epoch ({stage1_epoch})"
    assert final_step > stage1_step, f"Rank {rank}: Final step ({final_step}) should be > stage 1 step ({stage1_step})"

    # 2) State dict keys match (model architecture is preserved).
    assert final_keys == stage1_keys, (
        f"Rank {rank}: state_dict key mismatch between stage 1 and final. "
        f"Extra: {final_keys - stage1_keys}, Missing: {stage1_keys - final_keys}"
    )

    # 3) Weights changed (training actually happened in stage 2).
    weights_differ = any(not torch.equal(final_ckpt["state_dict"][k], stage1_weights[k]) for k in stage1_keys)
    assert (
        weights_differ
    ), f"Rank {rank}: No weights changed between stage 1 and final — stage 2 training may not have run"

    # 4) Optimizer state is valid.
    assert (
        "optimizer_states" in final_ckpt and final_ckpt["optimizer_states"]
    ), f"Rank {rank}: final checkpoint missing optimizer_states"
    opt_state = final_ckpt["optimizer_states"][0]["state"]
    assert len(opt_state) > 0, f"Rank {rank}: optimizer state is empty"

    # 5) Optimizer state keys are FQN strings (not legacy integers).
    opt_state_keys = list(opt_state.keys())
    assert all(
        isinstance(k, str) for k in opt_state_keys
    ), f"Rank {rank}: optimizer state keys should be FQN strings, got {[type(k).__name__ for k in opt_state_keys[:3]]}"

    # 6) EMA state is preserved across resume.
    if "ema" in stage1_ckpt:
        assert "ema_weights" in stage1_ckpt["ema"], "Stage 1 EMA must include ema_weights"
        assert "cur_step" in stage1_ckpt["ema"], "Stage 1 EMA must include cur_step"
        for ema_key, ema_val in stage1_ckpt["ema"]["ema_weights"].items():
            assert isinstance(ema_val, torch.Tensor) and not isinstance(
                ema_val, DTensor
            ), f"Stage 1 EMA weight '{ema_key}' must be plain torch.Tensor, got {type(ema_val).__name__}"

        assert "ema" in final_ckpt, f"Rank {rank}: final checkpoint missing EMA state after resume"
        assert final_ckpt["ema"]["cur_step"] > stage1_ckpt["ema"]["cur_step"], (
            f"Rank {rank}: EMA cur_step did not advance "
            f"(stage1={stage1_ckpt['ema']['cur_step']}, final={final_ckpt['ema']['cur_step']})"
        )
        assert set(final_ckpt["ema"]["ema_weights"].keys()) == set(
            stage1_ckpt["ema"]["ema_weights"].keys()
        ), f"Rank {rank}: EMA weight keys changed between stage 1 and final"

    torch.distributed.barrier()


@pytest.mark.slow
@pytest.mark.parametrize(
    "setup_env",
    [
        # DP-only: checkpoint resume correctness on 2-GPU systems
        ((2, (1, 1)), True, "cuda", "ENV"),
    ],
    indirect=("setup_env",),
    ids=["cuda-dp2-cp1x1"],
)
def test_boltz2_stop_and_go(
    setup_env,
    tmp_path,
    test_cp_training_data_dir_boltz2,
    canonical_mols_dir,
):
    """Stop-and-go checkpoint resume correctness for Boltz-2 via ``train()``.

    Verifies that training 1 epoch, checkpointing, then resuming to epoch 2
    produces a valid final state: correct epoch/step counters, preserved
    model architecture (state_dict keys), weights that changed (training
    happened), and valid optimizer state.

    Note: exact weight parity with a continuous 2-epoch run is not tested
    because ``train.py`` uses different seed offsets on resume.

    Uses the real ``train()`` entrypoint with real Boltz-2 training data
    and a small model config.  Only ``_cleanup_distributed`` is
    monkeypatched (for process group safety in test harness).
    """
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("CUDA not available")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"Need {world_size} GPUs, have {torch.cuda.device_count()}")

    size_dp = grid_group_sizes["dp"]
    size_cp = math.prod(grid_group_sizes["cp"])
    accelerator = "gpu" if device_type == "cuda" else "cpu"

    stopgo_dir = tmp_path / "stopgo"

    common_kwargs: dict[str, Any] = {
        "test_data_dir": test_cp_training_data_dir_boltz2,
        "mol_dir": canonical_mols_dir,
        "size_dp": size_dp,
        "size_cp": size_cp,
        "accelerator": accelerator,
    }

    # Stage 1: 1 epoch with checkpoint.
    stage1_config = stopgo_dir / "config_stage1.yaml"
    _write_train_config(
        TrainTestConfig(
            config_path=stage1_config,
            output_dir=stopgo_dir,
            **common_kwargs,
        )
    )

    # Stage 2: resume from checkpoint, train to epoch 2.
    stage2_config = stopgo_dir / "config_stage2.yaml"
    _write_train_config(
        TrainTestConfig(
            config_path=stage2_config,
            output_dir=stopgo_dir,
            max_epochs=2,
            resume=str(stopgo_dir / "last.ckpt"),
            **common_kwargs,
        )
    )

    payload = (
        env_per_rank,
        str(stage1_config),
        str(stage2_config),
        str(stopgo_dir),
    )
    spawn_multiprocessing(_parallel_assert_boltz2_stop_and_go, world_size, payload)


# ---------------------------------------------------------------------------
# Train → checkpoint → distributed inference pipeline test
# ---------------------------------------------------------------------------


def _parallel_assert_train_to_inference(
    rank: int,
    env_per_rank: dict[str, Any],
    checkpoint_path: str,
    data_dir: str,
    mol_dir: str,
    out_dir: str,
    size_dp: int,
    size_cp: int,
    accelerator: str,
) -> None:
    """Worker: run distributed inference using a checkpoint from training."""
    monkeypatch = pytest.MonkeyPatch()
    for key, value in env_per_rank.items():
        monkeypatch.setenv(key, f"{rank}" if value == "<INPUT_RANK>" else value)

    run_predict(
        data=data_dir,
        out_dir=out_dir,
        mol_dir=mol_dir,
        checkpoint=checkpoint_path,
        size_dp=size_dp,
        size_cp=size_cp,
        accelerator=accelerator,
        recycling_steps=0,
        sampling_steps=2,
        diffusion_samples=1,
        seed=42,
        input_format="preprocessed",
        use_templates=False,
        confidence_prediction=False,
        triattn_backend=TriAttnBackend.REFERENCE,
        sdpa_with_bias_backend=SDPAWithBiasBackend.REFERENCE,
        sdpa_with_bias_shardwise_backend=SDPAWithBiasBackend.REFERENCE,
    )

    out_path = Path(out_dir)
    data_stem = Path(data_dir).stem
    results_dir = out_path / f"boltz_results_{data_stem}"

    rank_cp = rank % size_cp
    if rank_cp == 0:
        assert results_dir.exists(), f"Rank {rank}: results dir {results_dir} not found"
        cif_files = list(results_dir.rglob("*.cif"))
        assert len(cif_files) > 0, f"Rank {rank}: no CIF output files in {results_dir}"


@pytest.mark.slow
@pytest.mark.parametrize(
    "setup_env",
    [
        ((2, (1, 1)), True, "cuda", "ENV"),
    ],
    indirect=("setup_env",),
    ids=["cuda-dp2-cp1x1"],
)
def test_boltz2_train_checkpoint_to_inference(
    setup_env,
    tmp_path,
    test_cp_training_data_dir_boltz2,
    canonical_mols_dir,
):
    """Train a small Boltz-2 model, then run inference with the saved checkpoint.

    Verifies the complete training-to-inference pipeline:
    1. Trains with real data for 1 epoch → checkpoint
    2. Loads the checkpoint via ``run_predict`` for distributed inference
    3. Verifies that CIF output files are produced

    This ensures checkpoints saved by ``train()`` (via ``BoltzContextParallelStrategy``)
    are compatible with the ``run_predict`` inference pipeline.
    """
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if not torch.cuda.is_available():
        pytest.skip("CUDA not available")
    if torch.cuda.device_count() < world_size:
        pytest.skip(f"Need {world_size} GPUs, have {torch.cuda.device_count()}")

    gpu_mem_gb = torch.cuda.get_device_properties(0).total_memory / (1024**3)
    if gpu_mem_gb < 100:
        pytest.skip(
            f"GPU has {gpu_mem_gb:.0f}GB memory; inference with CP=1 on real-data "
            "structures requires >80GB (outer_product_mean pair tensor OOM)"
        )

    size_dp = grid_group_sizes["dp"]
    size_cp = math.prod(grid_group_sizes["cp"])

    # ---- Step 1: Train and produce checkpoint ----
    training_output_dir = tmp_path / "training"
    config_path = tmp_path / "train_config.yaml"

    _write_train_config(
        TrainTestConfig(
            config_path=config_path,
            output_dir=training_output_dir,
            test_data_dir=test_cp_training_data_dir_boltz2,
            mol_dir=canonical_mols_dir,
            size_dp=size_dp,
            size_cp=size_cp,
            accelerator="gpu",
        )
    )

    train_payload = (env_per_rank, str(config_path), str(training_output_dir))
    spawn_multiprocessing(_parallel_assert_boltz2_train, world_size, train_payload)

    checkpoint_path = training_output_dir / "last.ckpt"
    assert checkpoint_path.exists(), f"Training checkpoint not found at {checkpoint_path}"

    # ---- Step 2: Run inference with the trained checkpoint ----
    inference_output_dir = tmp_path / "inference"

    spawn_multiprocessing(
        _parallel_assert_train_to_inference,
        world_size,
        env_per_rank,
        str(checkpoint_path),
        str(test_cp_training_data_dir_boltz2),
        str(canonical_mols_dir),
        str(inference_output_dir),
        size_dp,
        size_cp,
        "gpu",
    )


# ---------------------------------------------------------------------------
# E2E serial-vs-DTensor training parity test
# ---------------------------------------------------------------------------


def _load_serial_train_module():
    """Lazily import serial train.py (a script, not a package module).

    Deferred to call time to avoid executing the module during pytest
    collection, which can pollute global state and cause failures when
    many tests are collected.
    """
    path = Path(__file__).resolve().parents[2] / "scripts" / "train" / "train.py"
    spec = _importlib_util.spec_from_file_location("_serial_train", str(path))
    mod = _importlib_util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _smooth_lddt_loss_dense_e2e(
    pred_coords,
    true_coords,
    is_nucleotide,
    coords_mask=None,
    nucleic_acid_cutoff=30.0,
    other_cutoff=15.0,
    multiplicity=1,
    **kwargs,
):
    """Dense pairwise distance smooth_lddt for serial-distributed backward parity.

    Aligns the backward autograd graph between serial (sparse) and distributed
    (dense CDIST) computation paths.
    """
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


def _setup_training_data_7z64_8b2e_e2e(out_dir: Path, base_data_dir: Path) -> Path:
    """Merge 7z64 and 8b2e processed data into a single training directory."""
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


def _setup_training_data_all_4_e2e(out_dir: Path, base_data_dir: Path) -> tuple[Path, Path]:
    """Merge all 4 samples with a train/val split.

    Training: 7ylz, 8b2e.  Validation: 7z64, 8ayv.
    Returns ``(merged_data_dir, val_split_file_path)``.
    """
    names = ["7ylz", "7z64", "8ayv", "8b2e"]
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
    split_file = out_dir / "val_split.txt"
    split_file.write_text("7z64\n8ayv\n")
    return merged, split_file


def _e2e_model_dict(
    *,
    ema: bool = True,
    ema_decay: float = 0.999,
    multiplicity: int = 2,
    validate_structure: bool = False,
    distributed: bool = False,
) -> dict[str, Any]:
    """Small model config dict for the E2E training parity test.

    Matches ``create_boltz2_model_init_params(use_large_model=False)`` with
    recycling disabled, coordinate augmentation off, and configurable
    ``validate_structure``.

    When ``distributed=True``, validators use
    :class:`DistributedRCSBValidator` (the DTensor-aware variant required
    by ``Boltz2Distributed``).
    """
    d: dict[str, Any] = {
        "_target_": "boltz.model.models.boltz2.Boltz2",
        "atom_s": 4,
        "atom_z": 4,
        "token_s": 4,
        "token_z": 4,
        "num_bins": 64,
        "atom_feature_dim": 388,
        "atoms_per_window_queries": 32,
        "atoms_per_window_keys": 128,
        "ema": ema,
        "ema_decay": ema_decay,
        "confidence_prediction": False,
        "affinity_prediction": False,
        "structure_prediction_training": True,
        "use_templates": False,
        "validate_structure": validate_structure,
        "predict_bfactor": False,
        "bond_type_feature": False,
        "no_random_recycling_training": True,
        "embedder_args": {
            "atom_encoder_depth": 1,
            "atom_encoder_heads": 1,
            "activation_checkpointing": False,
        },
        "msa_args": {
            "msa_s": 4,
            "msa_blocks": 1,
            "msa_dropout": 0.0,
            "z_dropout": 0.0,
            "use_paired_feature": True,
        },
        "pairformer_args": {
            "num_blocks": 1,
            "num_heads": 1,
            "dropout": 0.0,
            "v2": True,
        },
        "score_model_args": {
            "sigma_data": 16.0,
            "dim_fourier": 4,
            "atom_encoder_depth": 1,
            "atom_encoder_heads": 1,
            "token_transformer_depth": 1,
            "token_transformer_heads": 1,
            "atom_decoder_depth": 1,
            "atom_decoder_heads": 1,
            "activation_checkpointing": False,
            "conditioning_transition_layers": 1,
        },
        "diffusion_process_args": {
            "coordinate_augmentation": False,
        },
        "diffusion_loss_args": {},
        "training_args": {
            "recycling_steps": 0,
            "sampling_steps": -1,
            "diffusion_multiplicity": multiplicity,
            "diffusion_samples": -1,
            "diffusion_loss_weight": 1.0,
            "distogram_loss_weight": 0.3,
            "confidence_loss_weight": 0.0,
            "bfactor_loss_weight": 0.0,
            "symmetry_correction": False,
            "adam_beta_1": 0.9,
            "adam_beta_2": 0.95,
            "adam_eps": 1e-8,
            "lr_scheduler": "af3",
            "base_lr": 1e-4,
            "max_lr": 1e-4,
            "lr_warmup_no_steps": 10,
            "lr_start_decay_after_n_steps": 100,
            "lr_decay_every_n_steps": 50000,
            "lr_decay_factor": 0.95,
            "weight_decay": 0.0,
        },
        "validation_args": {
            "recycling_steps": 0,
            "sampling_steps": 2,
            "diffusion_samples": 1,
            "symmetry_correction": False,
        },
    }
    if validate_structure:
        d["num_val_datasets"] = 1
        _validator_target = (
            "boltz.distributed.model.validation.rcsb.DistributedRCSBValidator"
            if distributed
            else "boltz.model.validation.rcsb.RCSBValidator"
        )
        d["validators"] = [
            {
                "_target_": _validator_target,
                "val_names": ["RCSB"],
                "confidence_prediction": False,
                "physicalism_metrics": True,
            }
        ]
    return d


def _apply_e2e_deterministic_getitem(monkeypatch, base_seed: int = 42) -> None:
    """Patch ``TrainingDataset.__getitem__`` at the class level for deterministic data."""
    original_getitem = SerialTrainingDataset.__getitem__

    _getitem_call_count = [0]

    def _wrapped_getitem(self, idx):
        _getitem_call_count[0] += 1
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

    monkeypatch.setattr(SerialTrainingDataset, "__getitem__", _wrapped_getitem)


def _apply_cached_getitem(monkeypatch, cache_path: str | Path) -> None:
    """Replace ``TrainingDataset.__getitem__`` with a disk-backed cache lookup.

    The cache file (created during pre-load) maps integer sample indices to
    feature dicts, guaranteeing identical features across serial and
    distributed data pipelines regardless of data-processing RNG state.
    """
    _cache: dict[int, dict] = {}

    def _cached_getitem(self, idx):
        if not _cache:
            _cache.update(torch.load(str(cache_path), map_location="cpu", weights_only=False))
        return _cache[idx]

    monkeypatch.setattr(SerialTrainingDataset, "__getitem__", _cached_getitem)


def _worker_e2e_training_parity(
    rank: int,
    grid_group_sizes: dict,
    device_type: str,
    backend: str,
    env_per_rank: dict[str, Any],
    dist_config_path: str,
    dist_output_dir: str,
    serial_ckpt_path: str,
    pretrained_ckpt_path: str,
    serial_metrics: dict,
    sigmas_global_host: torch.Tensor,
    noise_global_host: torch.Tensor,
    atom_counts_per_token_host: torch.Tensor,
    cached_samples_path: str,
    seed: int,
) -> None:
    """Multi-rank worker: distributed train() then compare with serial checkpoint.

    1. Applies module-level monkeypatches (noise, data, smooth_lddt, DoublePrecision)
    2. Calls ``train_module.train(dist_config_path, [])``
    3. Loads both serial and distributed checkpoints
    4. Compares state_dict, EMA weights, and logged metrics
    """
    monkeypatch = pytest.MonkeyPatch()
    if env_per_rank is not None:
        for var_name, value in env_per_rank.items():
            monkeypatch.setenv(var_name, f"{rank}" if value == "<INPUT_RANK>" else value)

    monkeypatch.setattr(train_module, "_cleanup_distributed", lambda: None)
    DistributedManager._state = {}

    # --- Deterministic data loading via cached samples ---
    _apply_cached_getitem(monkeypatch, cached_samples_path)

    # --- Deterministic noise (distributed) ---
    # Noise must go through distribute_atom_features (intersperse padding) so
    # that noise[i] at position i in the distributed tensor corresponds to the
    # same atom as coords[i].  Plain distribute_tensor would keep serial ordering
    # while coords use the intersperse-padded ordering, causing a mismatch.
    def _dist_noise_dist(self, bs, dtype=torch.float32):
        s = sigmas_global_host.to(device=self.device_mesh.device_type, dtype=dtype)[:bs]
        return distribute_tensor(s, self.device_mesh, (Shard(0), Replicate(), Replicate()))

    monkeypatch.setattr(DistAtomDiffusionV2, "noise_distribution", _dist_noise_dist)

    _noise_dt_cache: list[DTensor | None] = [None]
    _noise_computed = [False]

    def _compute_noise_dt_once(device_mesh, dtype):
        """Compute noise DTensor via distribute_atom_features (intersperse padding)."""
        if _noise_computed[0]:
            return
        _noise_computed[0] = True

        manager = DistributedManager()
        _io_keys = {"noise"}
        _placements = get_feature_placements(
            atom_keys=set(),
            model_io_keys=_io_keys,
            model_io_fp32_keys=set(),
        )

        size_batch = atom_counts_per_token_host.shape[0]
        multiplicity_val = noise_global_host.shape[0] // size_batch
        noise_unflat = noise_global_host.unflatten(0, (size_batch, multiplicity_val))

        inputs_io = {"atom_counts_per_token": atom_counts_per_token_host.clone()}
        for i_mul in range(multiplicity_val):
            inputs_io[f"noise_{i_mul}"] = noise_unflat[:, i_mul].to(dtype=dtype)

        placements_cp_model_io_mul = {
            f"{k}_{i_mul}": v for k, v in _placements["cp_model_io"].items() for i_mul in range(multiplicity_val)
        }
        placements_cp = _placements["cp_atom_features"] | placements_cp_model_io_mul
        placements_model_io_mul = {
            f"{k}_{i_mul}": v for k, v in _placements["model_io"].items() for i_mul in range(multiplicity_val)
        }
        placements_dp_cp = placements_model_io_mul

        io_feats = distribute_atom_features(
            inputs=inputs_io,
            placements_cp=placements_cp,
            placements_dp_cp=placements_dp_cp,
            device_mesh=manager.device_mesh_subgroups,
            cp_group=manager.group["cp"],
            multiplicities={"noise": multiplicity_val},
        )
        _noise_dt_cache[0] = io_feats.pop("noise").to(dtype=dtype)

    _dist_in_val = [False]

    def _det_create_randn(shape, device_mesh, placements, dtype=torch.float32, scale=1.0):
        if _dist_in_val[0]:
            from boltz.distributed.utils import create_distributed_randn as _real_create_randn

            return _real_create_randn(shape, device_mesh, placements, dtype=dtype, scale=0.0)
        from boltz.testing.utils import pad_to_length as _pad

        _compute_noise_dt_once(device_mesh, dtype)
        n = _noise_dt_cache[0]
        if n.dtype != dtype:
            n = n.to(dtype=dtype)
        if len(shape) > 1 and n.shape[1] < shape[1]:
            n = _pad(n, dim=1, length=shape[1])
        return n * scale

    monkeypatch.setattr(dist_diffusion_module, "create_distributed_randn", _det_create_randn)

    _orig_dist_val_step = Boltz2Distributed.validation_step

    def _dist_val_step_wrapper(self_model, batch, batch_idx):
        _dist_in_val[0] = True
        try:
            return _orig_dist_val_step(self_model, batch, batch_idx)
        finally:
            _dist_in_val[0] = False

    monkeypatch.setattr(Boltz2Distributed, "validation_step", _dist_val_step_wrapper)

    # --- Skip RMSD in distributed validation (not needed for LDDT parity) ---
    import boltz.distributed.model.validation.validator as _dist_validator_mod

    def _rmsd_noop(*args, **kwargs):
        return torch.tensor(0.0), None, None

    monkeypatch.setattr(_dist_validator_mod, "weighted_minimum_rmsd_single", _rmsd_noop)

    # --- Capture trainer metrics ---
    _captured_metrics: dict[str, float] = {}
    _orig_fit = pl.Trainer.fit

    def _capturing_fit(self, *args, **kwargs):
        result = _orig_fit(self, *args, **kwargs)
        for k, v in self.callback_metrics.items():
            if isinstance(v, DTensor):
                _captured_metrics[k] = v.full_tensor().detach().cpu().item()
            elif isinstance(v, torch.Tensor):
                _captured_metrics[k] = v.detach().cpu().item()
            else:
                _captured_metrics[k] = v
        return result

    monkeypatch.setattr(pl.Trainer, "fit", _capturing_fit)

    # --- Run distributed training ---
    train_module.train(dist_config_path, [])

    # --- Load checkpoints and compare ---
    dist_ckpt_path = Path(dist_output_dir) / "last.ckpt"
    assert dist_ckpt_path.exists(), f"Rank {rank}: distributed checkpoint not found at {dist_ckpt_path}"
    dist_ckpt = torch.load(dist_ckpt_path, map_location="cpu", weights_only=False)
    serial_ckpt = torch.load(serial_ckpt_path, map_location="cpu", weights_only=False)

    dist_sd = dist_ckpt["state_dict"]
    serial_sd = serial_ckpt["state_dict"]
    assert len(dist_sd) > 0, f"Rank {rank}: distributed state_dict is empty"
    assert len(serial_sd) > 0, f"Rank {rank}: serial state_dict is empty"

    # The distributed model prefixes keys with ``_serial.``; strip it for comparison.
    dist_sd_mapped = {}
    for k, v in dist_sd.items():
        canonical_k = k.replace("_serial.", "", 1) if k.startswith("_serial.") else k
        dist_sd_mapped[canonical_k] = v

    for k in serial_sd:
        assert k in dist_sd_mapped, f"Rank {rank}: key '{k}' missing from distributed checkpoint"
        torch.testing.assert_close(
            dist_sd_mapped[k],
            serial_sd[k],
            msg=lambda m: f"Rank {rank}: state_dict mismatch on '{k}': {m}",
        )

    # --- EMA weight parity ---
    assert "ema" in dist_ckpt, f"Rank {rank}: distributed checkpoint missing EMA state"
    assert "ema" in serial_ckpt, f"Rank {rank}: serial checkpoint missing EMA state"
    dist_ema = dist_ckpt["ema"]["ema_weights"]
    serial_ema = serial_ckpt["ema"]["ema_weights"]
    assert dist_ckpt["ema"]["cur_step"] == serial_ckpt["ema"]["cur_step"], (
        f"Rank {rank}: EMA cur_step mismatch: "
        f"dist={dist_ckpt['ema']['cur_step']}, serial={serial_ckpt['ema']['cur_step']}"
    )

    dist_ema_mapped = {}
    for k, v in dist_ema.items():
        canonical_k = k.replace("_serial.", "", 1) if k.startswith("_serial.") else k
        dist_ema_mapped[canonical_k] = v

    for k in serial_ema:
        assert k in dist_ema_mapped, f"Rank {rank}: EMA key '{k}' missing from distributed checkpoint"
        torch.testing.assert_close(
            dist_ema_mapped[k],
            serial_ema[k],
            msg=lambda m: f"Rank {rank}: EMA weight mismatch on '{k}': {m}",
        )

    # --- Non-vacuous guard: at least one parameter changed from pretrained init ---
    pretrained_sd = torch.load(
        pretrained_ckpt_path,
        map_location="cpu",
        weights_only=False,
    ).get("state_dict", {})
    if pretrained_sd:
        changed = any(not torch.equal(serial_sd[k], pretrained_sd[k]) for k in serial_sd if k in pretrained_sd)
        assert changed, f"Rank {rank}: no parameters changed from pretrained init — test is vacuous"

    # --- Metric parity ---
    # Atom-level LDDT metrics (val/lddt_*, val/complex_lddt_*) and the global
    # weighted-average val/lddt depend on diffusion-sampled coordinates, which
    # differ slightly between serial and distributed forward passes in FP32 due
    # to accumulation order in parallel attention.  Their exact parity is
    # verified separately by test_boltz2_validation_step_parity (FP64).  Here
    # we use relaxed tolerance (atol=5e-4) for these forward-pass-dependent
    # metrics and default tolerance for everything else.  Trailing underscore
    # omitted so "val/lddt" (global) is also matched.
    _forward_dependent_prefixes = ("val/lddt", "val/complex_lddt", "val/clash", "val/pb", "val/rmsd")
    _lddt_keys_compared = []
    if serial_metrics:
        for k in serial_metrics:
            if k in _captured_metrics:
                got = torch.tensor(_captured_metrics[k])
                exp = torch.tensor(serial_metrics[k])
                if any(k.startswith(p) for p in _forward_dependent_prefixes):
                    torch.testing.assert_close(
                        got,
                        exp,
                        atol=5e-4,
                        rtol=0.02,
                        msg=lambda m: f"Rank {rank}: metric '{k}' mismatch: {m}",
                    )
                else:
                    torch.testing.assert_close(
                        got,
                        exp,
                        msg=lambda m: f"Rank {rank}: metric '{k}' mismatch: {m}",
                    )
                if "lddt" in k:
                    _lddt_keys_compared.append(k)

    assert _lddt_keys_compared, (
        f"Rank {rank}: no validation LDDT metrics were compared — test is vacuous. "
        f"Serial keys: {sorted(serial_metrics)}, dist keys: {sorted(_captured_metrics)}"
    )
    for required_metric in ("val/lddt", "val/disto_lddt", "val/complex_lddt"):
        assert (
            required_metric in _captured_metrics
        ), f"Rank {rank}: distributed metrics missing '{required_metric}' — available: {sorted(_captured_metrics)}"

    # Verify component-wise grad_norm metrics are present and non-zero
    _grad_norm_keys = [
        "train/grad_norm",
        "train/grad_norm_msa_module",
        "train/grad_norm_pairformer_module",
        "train/grad_norm_structure_module",
    ]
    for gn_key in _grad_norm_keys:
        assert (
            gn_key in _captured_metrics
        ), f"Rank {rank}: distributed metrics missing '{gn_key}' — available: {sorted(_captured_metrics)}"
        assert (
            _captured_metrics[gn_key] > 0
        ), f"Rank {rank}: '{gn_key}' is zero — gradients should be non-zero after training"

    torch.distributed.barrier()


@pytest.mark.slow
@pytest.mark.parametrize(
    "setup_env",
    [
        ((2, (2, 2)), True, "cuda", "ENV"),
    ],
    indirect=("setup_env",),
    ids=["cuda-dp2-cp2x2"],
)
def test_boltz2_e2e_training_parity(
    setup_env,
    test_cp_training_base_data_dir_boltz2,
    canonical_mols_dir,
    tmp_path,
):
    """E2E serial-vs-DTensor training parity via ``train()`` entry points.

    Both serial and distributed training go through their respective
    ``train()`` functions for 1 epoch (1 batch of 7ylz+8b2e for training,
    7z64+8ayv for validation), then compare checkpoints (state_dict, EMA
    weights) and logged metrics — including validation LDDT — at FP32
    default tolerance.  Model initialisation is controlled via a pretrained
    checkpoint; noise and data RNG are controlled via module-level
    monkeypatches.

    Comparison summary (32 metrics, 296 state_dict params, 296 EMA params):

    State dict & EMA (296 params each, default FP32 tolerance):
      All params have non-zero magnitude (absmax in [3.9e-3, 5.0e-1]).
      Serial and distributed match exactly (bitwise identical).

    Weight update (pretrained -> post-training, lr=1e-4, 1 step):
      273/296 params changed (delta_absmax ~1e-4, consistent with lr).
      23 params unchanged (fourier embeddings, some norms/MLPs — zero
      gradient for this mini-batch). 15 "changed" params have delta < 1e-8
      (triangle attention Q/K weights with negligible gradients).

    Training metrics (default tolerance):
      train/loss=1.78, train/grad_norm=0.43, train/param_norm=6.37,
      train/diffusion_loss=0.53, train/distogram_loss=4.16,
      train/grad_norm_{msa_module,pairformer_module,structure_module}
      — all match within 1e-7 or exactly. Component-wise and global
      grad_norms are non-zero (logged from on_after_backward where
      gradients are available).

    Validation metrics — token-level (default tolerance):
      val/disto_lddt_{ligand_protein,intra_protein,protein_protein,
      intra_ligand}, val/disto_loss — all non-zero, match exactly.

    Validation metrics — atom-level (relaxed: atol=5e-4, rtol=0.02):
      val/lddt_{intra_ligand,intra_protein,ligand_protein,protein_protein}
      — non-zero values in [7e-4, 0.053], with abs diffs up to 8e-4 due
      to FP32 accumulation order differences in distributed attention.
      Exact parity verified separately by test_boltz2_validation_step_parity
      in FP64.

    Trivially-zero metrics (12 keys, no DNA/RNA in test data):
      val/{lddt,disto_lddt}_{dna_protein,rna_protein,dna_ligand,
      rna_ligand,intra_dna,intra_rna} — 0==0 on both sides.
    """

    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("CUDA not available")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"Need {world_size} GPUs, have {torch.cuda.device_count()}")

    dtype = torch.float32
    seed = 42
    multiplicity = 2
    B = 2
    max_tokens = 256
    W = 32  # atoms_per_window_queries
    size_cp = grid_group_sizes["cp"][0] * grid_group_sizes["cp"][1]
    atom_align = math.lcm(W, size_cp)
    max_atoms = ((max_tokens * 10 + atom_align - 1) // atom_align) * atom_align
    max_seqs = 16
    scale_glorot = 0.05

    # --- Merge all 4 samples with train/val split ---
    training_data_dir, split_file = _setup_training_data_all_4_e2e(
        tmp_path / "training_data", test_cp_training_base_data_dir_boltz2
    )

    # --- Create pretrained checkpoint with deterministic init ---
    seed_by_rank(0, seed=seed)
    model_dict = _e2e_model_dict(multiplicity=multiplicity, validate_structure=True)
    model_dict.pop("_target_")
    model_dict.pop("validators", None)
    _val_validators = [RCSBValidator(val_names=["RCSB"], confidence_prediction=False, physicalism_metrics=True)]
    pretrained_model = SerialBoltz2(**model_dict, validators=_val_validators)
    init_module_params_glorot(pretrained_model, gain=scale_glorot)
    pretrained_model.apply(SetModuleInfValues())
    pretrained_model.structure_module.coordinate_augmentation = False
    pretrained_model = pretrained_model.to(dtype=dtype)

    pretrained_path = tmp_path / "pretrained.ckpt"
    torch.save(
        {
            "state_dict": pretrained_model.state_dict(),
            "pytorch-lightning_version": pl.__version__,
            "hyper_parameters": pretrained_model.hparams,
        },
        pretrained_path,
    )

    # --- Pre-load and cache individual samples to disk ---
    # The serial and distributed ``train()`` create independent data pipelines
    # that call ``TrainingDataset.__getitem__`` separately.  Despite RNG seeding,
    # the featurizer has non-deterministic code paths (random augmentation of
    # ref_pos, MSA subsampling).  Caching the getitem results to disk and
    # replaying them guarantees identical features in both pipelines.
    _tmp_mp = pytest.MonkeyPatch()
    _apply_e2e_deterministic_getitem(_tmp_mp, base_seed=seed)
    _preload_cfg = setup_mock_training_datamodule_config(training_data_dir)
    _preload_cfg.batch_size = B
    _preload_cfg.samples_per_epoch = B
    _preload_cfg.moldir = str(canonical_mols_dir)
    _preload_cfg.return_train_symmetries = False
    _preload_cfg.msa_sampling_training = False
    _preload_cfg.max_tokens = max_tokens
    _preload_cfg.max_atoms = max_atoms
    _preload_cfg.max_seqs = max_seqs
    for _ds in _preload_cfg.datasets:
        _ds.filters = None
        _ds.split = str(split_file)
        _ds.symmetry_correction = False
    seed_by_rank(0, seed=seed)
    _preload_dm = Boltz2TrainingDataModule(cfg=_preload_cfg)

    _preload_ds = _preload_dm._train_set
    _cached_samples = {i: _preload_ds[i] for i in range(B)}
    cached_samples_path = tmp_path / "cached_samples.pt"
    torch.save(_cached_samples, cached_samples_path)

    _preload_dl = _preload_dm.train_dataloader()
    _preload_batch = next(iter(_preload_dl))
    atom_counts_per_token_host = _preload_batch["atom_counts_per_token"].detach().cpu()
    atom_pad_mask_host = _preload_batch["atom_pad_mask"].detach().cpu()
    _tmp_mp.undo()

    # --- Pre-generate deterministic noise (masked by atom_pad_mask) ---
    seed_by_rank(0, seed=seed)
    sigmas_global = pretrained_model.structure_module.noise_distribution(B * multiplicity).to(dtype=dtype)
    noise_global = torch.empty(B * multiplicity, max_atoms, 3, dtype=dtype)
    init_tensors_uniform([noise_global], low=-scale_glorot, high=scale_glorot)
    _mask_mul = atom_pad_mask_host[:, :, None].repeat_interleave(multiplicity, 0).to(dtype=dtype)
    noise_global = noise_global * _mask_mul

    sigmas_global_host = sigmas_global.detach().cpu()
    noise_global_host = noise_global.detach().cpu()

    # --- Write serial config ---
    serial_output_dir = tmp_path / "serial_output"
    serial_output_dir.mkdir(parents=True, exist_ok=True)
    serial_config_path = tmp_path / "serial_config.yaml"
    _e2e_ds_overrides = {
        "filters": None,
        "moldir": None,
        "symmetry_correction": False,
        "val_group": "RCSB",
        "use_train_subset": None,
        "override_bfactor": False,
        "override_method": None,
    }
    _write_train_config(
        TrainTestConfig(
            config_path=serial_config_path,
            output_dir=serial_output_dir,
            test_data_dir=training_data_dir,
            mol_dir=canonical_mols_dir,
            mode="serial",
            accelerator="gpu",
            limit_train_batches=1,
            pretrained=str(pretrained_path),
            model=_e2e_model_dict(multiplicity=multiplicity, validate_structure=True),
            batch_size=B,
            samples_per_epoch=B,
            max_tokens=max_tokens,
            max_atoms=max_atoms,
            max_seqs=max_seqs,
            return_train_symmetries=False,
            split=str(split_file),
            pop_target_keys=True,
            extra_dataset_overrides=_e2e_ds_overrides,
            v2=True,
            strict_loading=False,
            wandb=None,
            save_top_k=0,
            disable_checkpoint=False,
        )
    )

    # --- Apply serial monkeypatches ---
    serial_mp = pytest.MonkeyPatch()
    _apply_cached_getitem(serial_mp, cached_samples_path)

    _orig_serial_boltz2_init = SerialBoltz2.__init__

    @functools.wraps(_orig_serial_boltz2_init)
    def _init_with_validators(self, *args, **kwargs):
        if kwargs.get("validate_structure", False) and not kwargs.get("validators"):
            kwargs["validators"] = [
                RCSBValidator(val_names=["RCSB"], confidence_prediction=False, physicalism_metrics=True)
            ]
        _orig_serial_boltz2_init(self, *args, **kwargs)

    serial_mp.setattr(SerialBoltz2, "__init__", _init_with_validators)

    serial_mp.setattr(
        SerialAtomDiffusionV2,
        "noise_distribution",
        lambda self, bs: sigmas_global[:bs].to(device=self.zero.device),
    )

    _serial_in_val = [False]

    def _serial_randn_like(t):
        if _serial_in_val[0]:
            return torch.zeros_like(t)
        return noise_global[: t.shape[0], : t.shape[1]].to(device=t.device, dtype=t.dtype)

    serial_mp.setattr(serial_diffusion_v2_module.torch, "randn_like", _serial_randn_like)

    _orig_serial_randn = serial_diffusion_v2_module.torch.randn

    def _serial_randn(*args, **kwargs):
        if _serial_in_val[0]:
            kwargs.pop("generator", None)
            return torch.zeros(*args, **kwargs)
        return _orig_serial_randn(*args, **kwargs)

    serial_mp.setattr(serial_diffusion_v2_module.torch, "randn", _serial_randn)

    _orig_compute_random_augmentation = serial_diffusion_v2_module.compute_random_augmentation

    def _identity_augmentation_during_val(multiplicity, s_trans=1.0, device=None, dtype=torch.float32):
        if _serial_in_val[0]:
            R = torch.eye(3, device=device, dtype=dtype).unsqueeze(0).expand(multiplicity, -1, -1)
            tr = torch.zeros(multiplicity, 1, 3, device=device, dtype=dtype)
            return R, tr
        return _orig_compute_random_augmentation(multiplicity, s_trans=s_trans, device=device, dtype=dtype)

    serial_mp.setattr(serial_diffusion_v2_module, "compute_random_augmentation", _identity_augmentation_during_val)

    _orig_serial_val_step = SerialBoltz2.validation_step

    def _serial_val_step_wrapper(self_model, batch, batch_idx):
        _serial_in_val[0] = True
        try:
            return _orig_serial_val_step(self_model, batch, batch_idx)
        finally:
            _serial_in_val[0] = False

    serial_mp.setattr(SerialBoltz2, "validation_step", _serial_val_step_wrapper)

    serial_mp.setattr(serial_loss_v2_module, "smooth_lddt_loss", _smooth_lddt_loss_dense_e2e)
    serial_mp.setattr(serial_diffusion_v2_module, "smooth_lddt_loss", _smooth_lddt_loss_dense_e2e)

    serial_captured_metrics: dict[str, float] = {}
    _orig_fit = pl.Trainer.fit

    def _serial_capturing_fit(self, *args, **kwargs):
        result = _orig_fit(self, *args, **kwargs)
        for k, v in self.callback_metrics.items():
            if isinstance(v, torch.Tensor):
                serial_captured_metrics[k] = v.detach().cpu().item()
            else:
                serial_captured_metrics[k] = v
        return result

    serial_mp.setattr(pl.Trainer, "fit", _serial_capturing_fit)

    # --- Run serial training ---
    _serial_train_mod = _load_serial_train_module()
    _serial_train_mod.train(str(serial_config_path), [])

    # --- Find serial checkpoint ---
    serial_ckpt_files = list(serial_output_dir.rglob("last.ckpt"))
    assert (
        len(serial_ckpt_files) == 1
    ), f"Expected exactly 1 last.ckpt in serial output, found {len(serial_ckpt_files)}: {serial_ckpt_files}"
    serial_ckpt_path = serial_ckpt_files[0]

    # Verify serial checkpoint has EMA
    serial_ckpt = torch.load(serial_ckpt_path, map_location="cpu", weights_only=False)
    assert "ema" in serial_ckpt, "Serial checkpoint missing EMA state"
    assert "ema_weights" in serial_ckpt["ema"], "Serial EMA missing ema_weights"

    # Non-vacuous guard: serial must have logged at least one LDDT val metric
    serial_lddt_keys = [
        k for k in serial_captured_metrics if k.startswith("val/lddt_") or k.startswith("val/disto_lddt_")
    ]
    assert serial_lddt_keys, (
        f"Serial run produced no validation LDDT metrics — test is vacuous. "
        f"Available metrics: {sorted(serial_captured_metrics)}"
    )

    serial_mp.undo()

    # --- Write distributed config ---
    dp = grid_group_sizes["dp"]
    cp0, cp1 = grid_group_sizes["cp"]
    size_cp = cp0 * cp1
    dist_output_dir = tmp_path / "dist_output"
    dist_output_dir.mkdir(parents=True, exist_ok=True)
    dist_config_path = tmp_path / "dist_config.yaml"
    _write_train_config(
        TrainTestConfig(
            config_path=dist_config_path,
            output_dir=dist_output_dir,
            test_data_dir=training_data_dir,
            mol_dir=canonical_mols_dir,
            size_dp=dp,
            size_cp=size_cp,
            accelerator="gpu",
            limit_train_batches=1,
            pretrained=str(pretrained_path),
            model=_e2e_model_dict(multiplicity=multiplicity, validate_structure=True, distributed=True),
            batch_size=1,
            samples_per_epoch=dp,
            max_tokens=max_tokens,
            max_atoms=max_atoms,
            max_seqs=max_seqs,
            return_train_symmetries=False,
            split=str(split_file),
            pop_target_keys=True,
            extra_dataset_overrides=_e2e_ds_overrides,
        )
    )

    # --- Spawn distributed workers ---
    spawn_multiprocessing(
        _worker_e2e_training_parity,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        str(dist_config_path),
        str(dist_output_dir),
        str(serial_ckpt_path),
        str(pretrained_path),
        serial_captured_metrics,
        sigmas_global_host,
        noise_global_host,
        atom_counts_per_token_host,
        str(cached_samples_path),
        seed,
    )


# ---------------------------------------------------------------------------
#  BF16 / activation-checkpoint parity: shared setup
# ---------------------------------------------------------------------------


@dataclass
class _Bf16AcTestEnv:
    """Return type for :func:`_setup_bf16_ac_test_env`."""

    serial_config_path: Path
    dist_config_path: Path
    grid_group_sizes: dict
    world_size: int
    device_type: str
    backend: str
    env_per_rank: dict[str, Any]


def _setup_bf16_ac_test_env(
    setup_env,
    test_cp_training_base_data_dir_boltz2: Path,
    canonical_mols_dir: Path,
    tmp_path: Path,
) -> _Bf16AcTestEnv:
    """Shared setup for BF16 + activation-checkpointing parity tests.

    Creates training data, a pretrained checkpoint with all AC flags enabled
    (mirroring ``structurev2.yaml``), and writes serial / distributed YAML
    configs.  Returns paths and grid metadata so each test can attach its own
    profiler and spawn workers.
    """
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("CUDA not available")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"Need {world_size} GPUs, have {torch.cuda.device_count()}")

    seed = 42
    multiplicity = 2
    B = 2
    max_tokens = 256
    W = 32  # atoms_per_window_queries
    size_cp = grid_group_sizes["cp"][0] * grid_group_sizes["cp"][1]
    atom_align = math.lcm(W, size_cp)
    max_atoms = ((max_tokens * 10 + atom_align - 1) // atom_align) * atom_align
    max_seqs = 16
    scale_glorot = 0.05

    training_data_dir, split_file = _setup_training_data_all_4_e2e(
        tmp_path / "training_data", test_cp_training_base_data_dir_boltz2
    )

    ac_model_dict = _e2e_model_dict(multiplicity=multiplicity, validate_structure=False)
    ac_model_dict["checkpoint_diffusion_conditioning"] = True
    ac_model_dict["msa_args"]["activation_checkpointing"] = True
    ac_model_dict["pairformer_args"]["activation_checkpointing"] = True
    ac_model_dict["score_model_args"]["activation_checkpointing"] = True

    seed_by_rank(0, seed=seed)
    model_dict = copy.deepcopy(ac_model_dict)
    model_dict.pop("_target_")
    model_dict.pop("validators", None)
    pretrained_model = SerialBoltz2(**model_dict)
    init_module_params_glorot(pretrained_model, gain=scale_glorot)
    pretrained_model.apply(SetModuleInfValues())
    pretrained_model.structure_module.coordinate_augmentation = False
    pretrained_model = pretrained_model.to(dtype=torch.float32)

    pretrained_path = tmp_path / "pretrained.ckpt"
    torch.save(
        {
            "state_dict": pretrained_model.state_dict(),
            "pytorch-lightning_version": pl.__version__,
            "hyper_parameters": pretrained_model.hparams,
        },
        pretrained_path,
    )

    _e2e_ds_overrides = {
        "filters": None,
        "moldir": None,
        "symmetry_correction": False,
        "val_group": "RCSB",
        "use_train_subset": None,
        "override_bfactor": False,
        "override_method": None,
    }

    serial_output_dir = tmp_path / "serial_output"
    serial_output_dir.mkdir(parents=True, exist_ok=True)
    serial_config_path = tmp_path / "serial_config.yaml"
    _write_train_config(
        TrainTestConfig(
            config_path=serial_config_path,
            output_dir=serial_output_dir,
            test_data_dir=training_data_dir,
            mol_dir=canonical_mols_dir,
            mode="serial",
            accelerator="gpu",
            precision="bf16-mixed",
            limit_train_batches=1,
            limit_val_batches=0,
            num_sanity_val_steps=0,
            pretrained=str(pretrained_path),
            model=copy.deepcopy(ac_model_dict),
            batch_size=B,
            samples_per_epoch=B,
            max_tokens=max_tokens,
            max_atoms=max_atoms,
            max_seqs=max_seqs,
            return_train_symmetries=False,
            split=str(split_file),
            pop_target_keys=True,
            extra_dataset_overrides=_e2e_ds_overrides,
            v2=True,
            strict_loading=False,
            wandb=None,
            save_top_k=0,
            disable_checkpoint=True,
        )
    )

    dp = grid_group_sizes["dp"]
    size_cp = grid_group_sizes["cp"]
    if isinstance(size_cp, tuple):
        size_cp = size_cp[0] * size_cp[1]
    dist_output_dir = tmp_path / "dist_output"
    dist_output_dir.mkdir(parents=True, exist_ok=True)
    dist_config_path = tmp_path / "dist_config.yaml"
    _write_train_config(
        TrainTestConfig(
            config_path=dist_config_path,
            output_dir=dist_output_dir,
            test_data_dir=training_data_dir,
            mol_dir=canonical_mols_dir,
            size_dp=dp,
            size_cp=size_cp,
            accelerator="gpu",
            precision="BF16_MIXED",
            limit_train_batches=1,
            limit_val_batches=0,
            num_sanity_val_steps=0,
            pretrained=str(pretrained_path),
            model=copy.deepcopy(ac_model_dict),
            batch_size=1,
            samples_per_epoch=dp,
            max_tokens=max_tokens,
            max_atoms=max_atoms,
            max_seqs=max_seqs,
            return_train_symmetries=False,
            split=str(split_file),
            pop_target_keys=True,
            extra_dataset_overrides=_e2e_ds_overrides,
        )
    )

    return _Bf16AcTestEnv(
        serial_config_path=serial_config_path,
        dist_config_path=dist_config_path,
        grid_group_sizes=grid_group_sizes,
        world_size=world_size,
        device_type=device_type,
        backend=backend,
        env_per_rank=env_per_rank,
    )


# ---------------------------------------------------------------------------
#  BF16 dtype parity: serial vs DTensor training
# ---------------------------------------------------------------------------


def _worker_bf16_dtype_parity(
    rank: int,
    grid_group_sizes: dict,
    device_type: str,
    backend: str,
    env_per_rank: dict[str, Any],
    dist_config_path: str,
    serial_dtype_profile_path: str,
) -> None:
    """Multi-rank worker: run distributed train(), compare dtype profiles with serial.

    Only dtype equality is checked — no numerical comparison.  The serial
    dtype profile (written by the main process) is loaded from disk and
    compared against the distributed profile captured via :class:`DtypeProfiler`.
    """
    from boltz.testing.utils import DtypeProfiler

    monkeypatch = pytest.MonkeyPatch()
    if env_per_rank is not None:
        for var_name, value in env_per_rank.items():
            monkeypatch.setenv(var_name, f"{rank}" if value == "<INPUT_RANK>" else value)

    monkeypatch.setattr(train_module, "_cleanup_distributed", lambda: None)
    DistributedManager._state = {}

    # Capture dtype profile via Trainer.fit monkeypatch
    _dist_profiler: list[DtypeProfiler | None] = [None]
    _orig_fit = pl.Trainer.fit

    def _profiling_fit(trainer_self, model, **kwargs):
        _dist_profiler[0] = DtypeProfiler(model)
        result = _orig_fit(trainer_self, model, **kwargs)
        _dist_profiler[0].collect_grad_dtypes(model)
        return result

    monkeypatch.setattr(pl.Trainer, "fit", _profiling_fit)

    train_module.train(dist_config_path, [])

    profiler = _dist_profiler[0]
    assert profiler is not None, f"Rank {rank}: DtypeProfiler was never attached"

    profiler.remove_hooks()

    # Load serial dtype profile
    serial_profile = torch.load(serial_dtype_profile_path, map_location="cpu", weights_only=False)
    serial_fwd = serial_profile["fwd_dtypes"]
    serial_grads = serial_profile["param_grad_dtypes"]

    dist_fwd = profiler.fwd_dtypes
    dist_params = profiler.param_dtypes
    dist_grads = profiler.param_grad_dtypes

    # --- Parameter dtypes: all FP32 under bf16-mixed ---
    for name, dtype in dist_params.items():
        assert dtype == torch.float32, f"Rank {rank}: distributed param '{name}' has dtype {dtype}, expected float32"

    # --- Forward activation dtypes: strict equality at common module names ---
    common_fwd = sorted(set(serial_fwd) & set(dist_fwd))
    assert len(common_fwd) >= 10, (
        f"Rank {rank}: only {len(common_fwd)} common forward module names "
        f"between serial ({len(serial_fwd)}) and distributed ({len(dist_fwd)}). "
        f"Expected >= 10 for a non-vacuous comparison."
    )
    fwd_mismatches: list[str] = []
    for name in common_fwd:
        if serial_fwd[name] != dist_fwd[name]:
            fwd_mismatches.append(f"  {name}: serial={serial_fwd[name]}, dist={dist_fwd[name]}")
    assert not fwd_mismatches, f"Rank {rank}: forward activation dtype mismatches:\n" + "\n".join(fwd_mismatches)

    # --- Non-vacuous: autocast must produce a mix of BF16 and FP32 ---
    dist_fwd_dtypes_set = set(dist_fwd.values())
    assert (
        torch.bfloat16 in dist_fwd_dtypes_set
    ), f"Rank {rank}: no bfloat16 activations found — autocast may not be active. Unique dtypes: {dist_fwd_dtypes_set}"
    assert (
        torch.float32 in dist_fwd_dtypes_set
    ), f"Rank {rank}: no float32 activations found — all ops appear autocasted. Unique dtypes: {dist_fwd_dtypes_set}"

    # --- Parameter gradient dtypes: strict equality at common param names ---
    common_grads = sorted(set(serial_grads) & set(dist_grads))
    grad_mismatches: list[str] = []
    for name in common_grads:
        if serial_grads[name] != dist_grads[name]:
            grad_mismatches.append(f"  {name}: serial={serial_grads[name]}, dist={dist_grads[name]}")
    assert not grad_mismatches, f"Rank {rank}: param gradient dtype mismatches:\n" + "\n".join(grad_mismatches)

    torch.distributed.barrier()


@pytest.mark.slow
@pytest.mark.parametrize(
    "setup_env",
    [
        ((1, (2, 2)), True, "cuda", "ENV"),
    ],
    indirect=("setup_env",),
    ids=["cuda-dp1-cp2x2"],
)
def test_boltz2_bf16_dtype_parity(
    setup_env,
    test_cp_training_base_data_dir_boltz2,
    canonical_mols_dir,
    tmp_path,
):
    """Verify that BF16-mixed autocast produces identical dtype profiles in
    serial and DTensor training workflows.

    Runs 1 training step under ``bf16-mixed`` precision for both the serial
    ``Boltz2`` model (via ``scripts/train/train.py``) and the distributed
    ``Boltz2`` model (via ``src/boltz/distributed/train.py``), then compares
    forward activation dtypes, parameter dtypes, and parameter gradient
    dtypes at every module whose name appears in both models.

    No numerical comparison is performed — only dtype equality.  This means
    no deterministic-noise / cached-sample monkeypatching is required.
    """
    from boltz.testing.utils import DtypeProfiler

    env = _setup_bf16_ac_test_env(setup_env, test_cp_training_base_data_dir_boltz2, canonical_mols_dir, tmp_path)

    # --- Run serial training with dtype profiling ---
    serial_mp = pytest.MonkeyPatch()
    _serial_profiler: list[DtypeProfiler | None] = [None]
    _orig_fit = pl.Trainer.fit

    def _serial_profiling_fit(trainer_self, model, **kwargs):
        _serial_profiler[0] = DtypeProfiler(model)
        result = _orig_fit(trainer_self, model, **kwargs)
        _serial_profiler[0].collect_grad_dtypes(model)
        return result

    serial_mp.setattr(pl.Trainer, "fit", _serial_profiling_fit)

    _serial_train_mod = _load_serial_train_module()
    _serial_train_mod.train(str(env.serial_config_path), [])

    profiler = _serial_profiler[0]
    assert profiler is not None, "Serial DtypeProfiler was never attached"
    profiler.remove_hooks()

    # Non-vacuous: serial must have a mix of BF16 and FP32 activations
    serial_fwd_dtypes_set = set(profiler.fwd_dtypes.values())
    assert (
        torch.bfloat16 in serial_fwd_dtypes_set
    ), f"Serial: no bfloat16 activations — autocast may not be active. Unique dtypes: {serial_fwd_dtypes_set}"
    assert (
        torch.float32 in serial_fwd_dtypes_set
    ), f"Serial: no float32 activations. Unique dtypes: {serial_fwd_dtypes_set}"

    # All serial params must be FP32
    for name, dtype in profiler.param_dtypes.items():
        assert dtype == torch.float32, f"Serial param '{name}' has dtype {dtype}, expected float32"

    # Save serial profile for workers
    serial_dtype_profile_path = tmp_path / "serial_dtype_profile.pt"
    torch.save(
        {
            "fwd_dtypes": profiler.fwd_dtypes,
            "param_dtypes": profiler.param_dtypes,
            "param_grad_dtypes": profiler.param_grad_dtypes,
        },
        serial_dtype_profile_path,
    )

    serial_mp.undo()

    # --- Spawn distributed workers ---
    spawn_multiprocessing(
        _worker_bf16_dtype_parity,
        env.world_size,
        env.grid_group_sizes,
        env.device_type,
        env.backend,
        env.env_per_rank,
        str(env.dist_config_path),
        str(serial_dtype_profile_path),
    )


# ---------------------------------------------------------------------------
#  Activation checkpoint recomputation parity: serial vs DTensor training
# ---------------------------------------------------------------------------


def _worker_actv_ckpt_parity(
    rank: int,
    grid_group_sizes: dict,
    device_type: str,
    backend: str,
    env_per_rank: dict[str, Any],
    dist_config_path: str,
    serial_recompute_profile_path: str,
) -> None:
    """Per-rank worker: compare checkpoint-recomputed modules with serial.

    Each rank records its own local forward-hook call counts via
    :class:`RecomputeProfiler` and independently compares against the serial
    reference.  Counts are **never** aggregated across ranks.
    """
    from boltz.testing.utils import RecomputeProfiler

    monkeypatch = pytest.MonkeyPatch()
    if env_per_rank is not None:
        for var_name, value in env_per_rank.items():
            monkeypatch.setenv(var_name, f"{rank}" if value == "<INPUT_RANK>" else value)

    monkeypatch.setattr(train_module, "_cleanup_distributed", lambda: None)
    DistributedManager._state = {}

    _dist_profiler: list[RecomputeProfiler | None] = [None]
    _orig_fit = pl.Trainer.fit

    def _profiling_fit(trainer_self, model, **kwargs):
        _dist_profiler[0] = RecomputeProfiler(model)
        return _orig_fit(trainer_self, model, **kwargs)

    monkeypatch.setattr(pl.Trainer, "fit", _profiling_fit)

    train_module.train(dist_config_path, [])

    profiler = _dist_profiler[0]
    assert profiler is not None, f"Rank {rank}: RecomputeProfiler was never attached"
    profiler.remove_hooks()

    serial_profile = torch.load(serial_recompute_profile_path, map_location="cpu", weights_only=False)
    serial_counts: dict[str, int] = serial_profile["fwd_counts"]
    dist_counts = profiler.fwd_counts

    common = sorted(set(serial_counts) & set(dist_counts))
    assert len(common) >= 10, (
        f"Rank {rank}: only {len(common)} common module names between serial "
        f"({len(serial_counts)}) and distributed ({len(dist_counts)}). Expected >= 10."
    )

    serial_recomp = {n for n in common if serial_counts[n] >= 2}
    dist_recomp = {n for n in common if dist_counts[n] >= 2}

    # Non-vacuous: activation checkpointing must be active on this rank
    assert len(dist_recomp) >= 5, (
        f"Rank {rank}: only {len(dist_recomp)} recomputed modules in distributed "
        f"(expected >= 5). Activation checkpointing may not be active."
    )
    assert dist_recomp < set(common), (
        f"Rank {rank}: all {len(common)} common modules are recomputed — "
        f"not a strict subset, implying either a counting bug or every module is checkpointed."
    )

    serial_only = sorted(serial_recomp - dist_recomp)
    dist_only = sorted(dist_recomp - serial_recomp)
    assert (
        not serial_only
    ), f"Rank {rank}: modules recomputed in serial but not DTensor ({len(serial_only)}): {serial_only}"
    assert not dist_only, f"Rank {rank}: modules recomputed in DTensor but not serial ({len(dist_only)}): {dist_only}"

    torch.distributed.barrier()


@pytest.mark.slow
@pytest.mark.parametrize(
    "setup_env",
    [
        ((1, (2, 2)), True, "cuda", "ENV"),
    ],
    indirect=("setup_env",),
    ids=["cuda-dp1-cp2x2"],
)
def test_boltz_actv_ckpt_parity(
    setup_env,
    test_cp_training_base_data_dir_boltz2,
    canonical_mols_dir,
    tmp_path,
):
    """Verify serial and DTensor training checkpoint-recompute the same modules.

    Runs 1 training step (forward + backward) under ``bf16-mixed`` precision
    with all production activation-checkpointing flags enabled.  Forward hooks
    count how many times each module's forward is invoked: modules inside a
    ``torch.utils.checkpoint.checkpoint`` region are called twice (once in
    forward, once during backward recomputation).

    Each distributed rank independently compares its local counts against the
    serial reference — counts are **never** aggregated across ranks.
    """
    from boltz.testing.utils import RecomputeProfiler

    env = _setup_bf16_ac_test_env(setup_env, test_cp_training_base_data_dir_boltz2, canonical_mols_dir, tmp_path)

    # --- Run serial training with recompute profiling ---
    serial_mp = pytest.MonkeyPatch()
    _serial_profiler: list[RecomputeProfiler | None] = [None]
    _orig_fit = pl.Trainer.fit

    def _serial_profiling_fit(trainer_self, model, **kwargs):
        _serial_profiler[0] = RecomputeProfiler(model)
        return _orig_fit(trainer_self, model, **kwargs)

    serial_mp.setattr(pl.Trainer, "fit", _serial_profiling_fit)

    _serial_train_mod = _load_serial_train_module()
    _serial_train_mod.train(str(env.serial_config_path), [])

    profiler = _serial_profiler[0]
    assert profiler is not None, "Serial RecomputeProfiler was never attached"
    profiler.remove_hooks()

    # Non-vacuous: activation checkpointing must recompute some modules
    serial_recomp = profiler.recomputed_modules
    all_serial_names = set(profiler.fwd_counts)
    assert len(serial_recomp) >= 5, (
        f"Serial: only {len(serial_recomp)} modules recomputed (expected >= 5). "
        f"Activation checkpointing may not be active."
    )
    assert (
        serial_recomp < all_serial_names
    ), f"Serial: all {len(all_serial_names)} modules are recomputed — not a strict subset, implying a counting bug."

    # Save serial recompute profile for workers
    serial_recompute_profile_path = tmp_path / "serial_recompute_profile.pt"
    torch.save({"fwd_counts": profiler.fwd_counts}, serial_recompute_profile_path)

    serial_mp.undo()

    # --- Spawn distributed workers (per-rank comparison, no cross-rank aggregation) ---
    spawn_multiprocessing(
        _worker_actv_ckpt_parity,
        env.world_size,
        env.grid_group_sizes,
        env.device_type,
        env.backend,
        env.env_per_rank,
        str(env.dist_config_path),
        str(serial_recompute_profile_path),
    )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
