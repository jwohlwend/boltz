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

"""Distributed training entrypoint for Boltz-2 with DTensor context parallelism.

This module wires the full distributed training infrastructure:

- Distributed process-group bootstrap via :class:`DistributedManager`
- :class:`Boltz2Distributed` model wrapping with DTensor CP
- :class:`Boltz2TrainingDataModule` with DTensor feature distribution
- :class:`BoltzContextParallelStrategy` for DTensor checkpoint save/load
- Checkpoint callback with Boltz-2 defaults (``monitor="val/lddt"``, etc.)
- Resume / pretrained-loading plumbing
- Precision configuration (``bf16``, ``bf16-mixed``, ``tf32``, etc.)
- WandB logging and config serialization

Factory functions :func:`_create_distributed_model` and
:func:`_create_distributed_data_module` are extracted as module-level
functions so tests can monkeypatch them with lightweight smoke
implementations (see ``tests/distributed/test_dtensor_stop_and_go.py``).
"""

from __future__ import annotations

import atexit
import os
import random
import string
import sys
import warnings
from collections import OrderedDict
from dataclasses import dataclass
from datetime import timedelta
from math import isqrt
from pathlib import Path
from typing import Any, Optional

import hydra
import omegaconf
import pytorch_lightning as pl
import torch
from hydra import compose, initialize_config_dir
from hydra.core.global_hydra import GlobalHydra
from omegaconf import OmegaConf
from pytorch_lightning import LightningModule, seed_everything
from pytorch_lightning.callbacks import ModelCheckpoint
from pytorch_lightning.loggers import WandbLogger
from pytorch_lightning.utilities import rank_zero_only

try:
    from one_logger_utils.pytorch_lightning import hook_trainer_cls  # type: ignore[import-untyped]

    _one_logger_available = True
except ImportError:
    _one_logger_available = False

from boltz.data.module.trainingv2 import DataConfigV2
from boltz.distributed.data.module.trainingv2 import Boltz2TrainingDataModule
from boltz.distributed.data.utils import map_subgroup_mesh_to_cpu
from boltz.distributed.lightning_strategy import BoltzContextParallelStrategy
from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.models.boltz2 import Boltz2 as Boltz2Distributed
from boltz.distributed.model.modules.utils import (
    PRECISION_TO_LIGHTNING,
    OffloadActvCkptToCPU,
    Precision,
    SDPAWithBiasBackend,
    SetAttnPairBiasBackend,
    SetAttnPairBiasShardwiseBackend,
    SetTriAttnBackend,
    TriAttnBackend,
    setup_tf32_env,
)
from boltz.model.layers.attentionv2 import AttentionPairBias as AttentionPairBiasV2
from boltz.workflow.utils import (
    _DATASET_KEYS_TO_OVERRIDE,
    CUDAMemoryProfile,
    convert_datasets_dict_to_list_config,
)


@dataclass
class DistributedTrainConfig:
    """Configuration dataclass for distributed CP training."""

    data: Any  # DataConfigV2 or OmegaConf dict; converted in _create_distributed_data_module
    model: LightningModule
    output: str
    trainer: Optional[dict[str, Any]] = None
    parallel_size: Optional[dict[str, Any]] = None
    precision: Precision = Precision.FP32
    matmul_precision: Optional[str] = None
    find_unused_parameters: Optional[bool] = False  # Retained for boltz-2 config compat; unused by CP strategy
    save_top_k: Optional[int] = 1
    checkpoint: Optional[dict[str, Any]] = None
    resume: Optional[str] = None
    pretrained: Optional[str] = None
    wandb: Optional[dict[str, Any]] = None
    disable_checkpoint: bool = False
    debug: bool = False
    strict_loading: bool = True
    load_confidence_from_trunk: Optional[bool] = False
    seed: Optional[int] = None
    validation_only: bool = False
    v2: bool = True  # Retained for structurev2.yaml compat; always True for Boltz-2
    triattn_backend: TriAttnBackend = TriAttnBackend.CUEQ
    sdpa_with_bias_backend: SDPAWithBiasBackend = SDPAWithBiasBackend.TORCH_FLEX_ATTN
    sdpa_with_bias_shardwise_backend: SDPAWithBiasBackend = SDPAWithBiasBackend.TORCH_FLEX_ATTN


def _load_and_merge_config(raw_config: str, args: list[str]) -> omegaconf.DictConfig:
    raw_config_dict = omegaconf.OmegaConf.load(raw_config)
    if "defaults" in raw_config_dict:
        config_path = Path(raw_config)
        GlobalHydra.instance().clear()
        with initialize_config_dir(config_dir=str(config_path.parent.absolute()), version_base=None):
            config_dict = compose(config_name=config_path.stem)
    else:
        config_dict = raw_config_dict
    omegaconf.OmegaConf.set_struct(config_dict, False)
    args_dict = omegaconf.OmegaConf.from_dotlist(args)
    if (
        "data" in args_dict
        and "datasets" in args_dict.data
        and "data" in config_dict
        and "datasets" in config_dict.data
    ):
        args_dict["data"]["datasets"] = convert_datasets_dict_to_list_config(
            config_dict.data.datasets,
            args_dict.data.datasets,
            keys_to_override=_DATASET_KEYS_TO_OVERRIDE,
            remove_null_datasets=True,
        )
    return omegaconf.OmegaConf.merge(config_dict, args_dict)


def _parse_precision(value: Any) -> Precision:
    if isinstance(value, Precision):
        return value
    if isinstance(value, str):
        if value in Precision.__members__:
            return Precision[value]
        for precision in Precision:
            if precision.value == value:
                return precision
    raise ValueError(f"Unsupported precision value: {value!r}")


def _parse_backend_enum(value: Any, enum_cls: type) -> Any:
    """Parse a string or enum value into the given enum class."""
    if isinstance(value, enum_cls):
        return value
    if isinstance(value, str):
        try:
            return enum_cls(value)
        except ValueError:
            pass
        if value in enum_cls.__members__:
            return enum_cls[value]
    raise ValueError(f"Unsupported {enum_cls.__name__} value: {value!r}. Valid: {[e.value for e in enum_cls]}")


def _apply_matmul_precision(matmul_precision: Optional[str]) -> None:
    """Apply optional matmul precision setting."""
    if matmul_precision is not None:
        torch.set_float32_matmul_precision(matmul_precision)


def _create_dist_manager(cfg: DistributedTrainConfig) -> DistributedManager:
    trainer_cfg = cfg.trainer or {}
    parallel_cfg = cfg.parallel_size or {}

    accelerator = trainer_cfg.get("accelerator", "gpu")
    accelerator_to_device_type = {"cpu": "cpu", "gpu": "cuda"}
    if accelerator not in accelerator_to_device_type:
        raise ValueError(
            f"Accelerator {accelerator} is not supported; expected one of {sorted(accelerator_to_device_type)}"
        )
    device_type = accelerator_to_device_type[accelerator]

    size_dp = int(parallel_cfg.get("size_dp", 1))
    size_cp = int(parallel_cfg.get("size_cp", 1))
    if size_dp <= 0 or size_cp <= 0:
        raise ValueError(f"size_dp and size_cp must be positive; got size_dp={size_dp}, size_cp={size_cp}")

    # If already initialized (e.g. train() called multiple times in the same
    # process during testing), validate that the requested topology matches the
    # existing singleton and return it.
    if DistributedManager.is_initialized():
        existing = DistributedManager()
        existing_device_type = existing.device.type
        if existing_device_type != device_type:
            raise ValueError(
                f"DistributedManager already initialized with device_type={existing_device_type!r}, "
                f"but this call requests device_type={device_type!r}. "
                f"Cannot change device type without cleanup + reinitialization."
            )
        existing_dp_size = len(existing.group_ranks.get("dp", []))
        existing_cp_size = len(existing.group_ranks.get("cp", []))
        if existing_dp_size and existing_dp_size != size_dp:
            raise ValueError(
                f"DistributedManager already initialized with dp group size {existing_dp_size}, "
                f"but this call requests size_dp={size_dp}. "
                f"Cannot change topology without cleanup + reinitialization."
            )
        if existing_cp_size and existing_cp_size != size_cp:
            raise ValueError(
                f"DistributedManager already initialized with cp group size {existing_cp_size}, "
                f"but this call requests size_cp={size_cp}. "
                f"Cannot change topology without cleanup + reinitialization."
            )
        return existing

    timeout_nccl_minutes = parallel_cfg.get("timeout_nccl")
    timeout_gloo_minutes = parallel_cfg.get("timeout_gloo")
    if timeout_nccl_minutes is not None and timeout_nccl_minutes <= 0:
        raise ValueError("timeout_nccl must be positive when provided")
    if timeout_gloo_minutes is not None and timeout_gloo_minutes <= 0:
        raise ValueError("timeout_gloo must be positive when provided")

    timeout_nccl = timedelta(minutes=timeout_nccl_minutes) if timeout_nccl_minutes is not None else None
    timeout_gloo = timedelta(minutes=timeout_gloo_minutes) if timeout_gloo_minutes is not None else None
    timeout_by_device = {"cuda": timeout_nccl, "cpu": timeout_gloo}

    DistributedManager.initialize(device_type=device_type, timeout=timeout_by_device[device_type])
    atexit.register(DistributedManager.cleanup)
    dist_manager = DistributedManager()
    if not dist_manager.has_dist:
        raise RuntimeError(
            "DistributedManager did not initialize torch.distributed. "
            "Launch this entrypoint under torchrun/slurm with RANK/WORLD_SIZE (or SLURM_* env)."
        )

    if size_dp * size_cp != dist_manager.world_size:
        raise ValueError(
            f"world_size mismatch: process world_size={dist_manager.world_size}, "
            f"expected size_dp*size_cp={size_dp * size_cp}"
        )

    size_cp_axis = isqrt(size_cp)
    if size_cp_axis * size_cp_axis != size_cp:
        raise ValueError(f"size_cp must be a square integer for 2D CP mesh, got {size_cp}")

    grid_group_sizes: OrderedDict[str, int | tuple[int, ...]] = OrderedDict(
        [("dp", size_dp), ("cp", (size_cp_axis, size_cp_axis))]
    )
    DistributedManager.create_grid_group(grid_group_sizes)
    return dist_manager


def _load_pretrained_if_requested(
    model_module: LightningModule,
    cfg: DistributedTrainConfig,
) -> LightningModule:
    """Load pretrained weights into ``model_module`` when ``cfg.pretrained`` is set.

    Returns the model unchanged when no pretrained path is configured or when
    resuming from a training checkpoint (``cfg.resume``).

    When ``cfg.load_confidence_from_trunk`` is True, trunk weights (everything
    except ``structure_module`` and ``distogram_module``) are duplicated under
    the ``confidence_module.`` prefix before loading, so the confidence head
    inherits shared encoder parameters.

    Loading uses ``strict=False`` to support reduced-depth fine-tuning (fewer
    pairformer layers).  ``_validate_checkpoint_architecture`` is called
    post-load to guard against silent V1/V2 attention mismatches that
    ``strict=False`` would otherwise ignore.
    """
    if not cfg.pretrained or cfg.resume:
        return model_module

    if cfg.load_confidence_from_trunk:
        checkpoint = torch.load(cfg.pretrained, map_location="cpu", weights_only=False)
        new_state_dict = {}
        for key, value in checkpoint["state_dict"].items():
            if not key.startswith("structure_module") and not key.startswith("distogram_module"):
                new_state_dict[f"confidence_module.{key}"] = value
        new_state_dict.update(checkpoint["state_dict"])
        checkpoint["state_dict"] = new_state_dict
        random_string = "".join(random.choices(string.ascii_lowercase + string.digits, k=10))
        temp_path = os.path.join(cfg.output, f".tmp_{random_string}.ckpt")
        torch.save(checkpoint, temp_path)
        file_path = temp_path
    else:
        file_path = cfg.pretrained

    hparams = dict(model_module.hparams)
    if getattr(model_module, "validate_structure", False) and hasattr(model_module, "validators"):
        hparams["validators"] = model_module.validators

    loaded = type(model_module).load_from_checkpoint(
        file_path,
        map_location="cpu",
        strict=False,
        **hparams,
    )
    if cfg.load_confidence_from_trunk:
        os.remove(file_path)

    _validate_checkpoint_architecture(loaded)
    return loaded


def _validate_checkpoint_architecture(model: LightningModule) -> None:
    """Verify the loaded model uses the expected attention implementation.

    When ``load_from_checkpoint`` uses ``strict=False``, mismatched hparams
    (e.g. missing ``v2=True`` in ``pairformer_args``) can silently create V1
    attention layers whose extra ``norm_s`` weights are randomly initialized
    instead of loaded from the checkpoint.  This validation catches that.
    """
    pairformer = getattr(model, "pairformer_module", None)
    if pairformer is None:
        return
    layers = getattr(pairformer, "layers", [])
    if not layers:
        return
    layer0_attn = getattr(layers[0], "attention", None)
    if layer0_attn is None:
        return
    if not isinstance(layer0_attn, AttentionPairBiasV2):
        raise RuntimeError(
            f"Pairformer layer 0 attention is {type(layer0_attn).__module__}."
            f"{type(layer0_attn).__name__}, expected AttentionPairBiasV2 "
            f"(boltz.model.layers.attentionv2). This usually means "
            f"pairformer_args is missing 'v2: true' — pass "
            f"pairformer_args=asdict(PairformerArgsV2()) to "
            f"load_from_checkpoint."
        )


def _cleanup_distributed() -> None:
    """Clean up distributed process groups after training.

    Extracted as a module-level function so tests that call ``train()``
    multiple times in the same worker process can monkeypatch it to a
    no-op and defer cleanup to the test's own ``finally`` block.
    """
    if DistributedManager.is_initialized():
        DistributedManager.cleanup()


def _create_distributed_data_module(
    data_config: Any,
    dist_manager: DistributedManager,
) -> pl.LightningDataModule:
    """Construct the distributed Boltz-2 training data module.

    Wraps the serial ``Boltz2TrainingDataModule`` with DTensor context-parallel
    distribution.  Tests may monkeypatch this function to supply a lightweight
    DTensor-producing smoke data module.

    Parameters
    ----------
    data_config
        Data configuration — either a :class:`DataConfigV2` instance or an
        OmegaConf/dict that can be unpacked into one.
    dist_manager
        Initialized :class:`DistributedManager` with grid groups.
    """
    cfg = data_config if isinstance(data_config, DataConfigV2) else DataConfigV2(**data_config)
    device_mesh = dist_manager.device_mesh_subgroups
    device_mesh_cpu = map_subgroup_mesh_to_cpu(dist_manager)
    return Boltz2TrainingDataModule(cfg=cfg, device_mesh=device_mesh, device_mesh_cpu=device_mesh_cpu)


def _create_distributed_model(
    cfg: DistributedTrainConfig,
    dist_manager: DistributedManager,
) -> LightningModule:
    """Construct the distributed Boltz-2 model with DTensor CP wrapping.

    Instantiates the serial model from config, loads pretrained weights if
    requested, moves to device, and wraps with :class:`Boltz2Distributed`.
    Tests may monkeypatch this function to supply a lightweight DTensor-aware
    smoke model.

    Parameters
    ----------
    cfg
        Full training configuration (model, pretrained, strict_loading, etc.).
    dist_manager
        Initialized :class:`DistributedManager` with grid groups.
    """
    model_serial = cfg.model
    model_serial = _load_pretrained_if_requested(model_serial, cfg)
    model_serial = model_serial.to(dist_manager.device)
    dist_model = Boltz2Distributed(model_serial, dist_manager)
    if not cfg.strict_loading:
        dist_model.strict_loading = False
    return dist_model


def train(raw_config: str, args: list[str]) -> None:  # noqa: C901, PLR0912
    """Run distributed training scaffold with strategy/checkpoint wiring."""
    config_dict = _load_and_merge_config(raw_config, args)
    if "precision" in config_dict:
        config_dict.precision = _parse_precision(config_dict.precision)
    for backend_key, backend_cls in (
        ("triattn_backend", TriAttnBackend),
        ("sdpa_with_bias_backend", SDPAWithBiasBackend),
        ("sdpa_with_bias_shardwise_backend", SDPAWithBiasBackend),
    ):
        if backend_key in config_dict and isinstance(config_dict[backend_key], str):
            config_dict[backend_key] = _parse_backend_enum(config_dict[backend_key], backend_cls)

    cuda_memory_profile_cfg = config_dict.pop("CUDAMemoryProfile", None)
    offload_actv_ckpt_cfg = config_dict.pop("OffloadActvCkptToCPU", None)

    cfg = hydra.utils.instantiate(config_dict)
    cfg = DistributedTrainConfig(**cfg)
    if not cfg.v2:
        raise NotImplementedError("DTensor distributed training only supports Boltz-2 (v2=true)")
    Path(cfg.output).mkdir(parents=True, exist_ok=True)
    _apply_matmul_precision(cfg.matmul_precision)

    dist_manager = _create_dist_manager(cfg)

    # Offset RNG seed by rank and, on resume, by epoch + global_step to avoid
    # replaying identical data samples.  Boltz's TrainingDataset ignores the
    # sampler index, so without a resume-aware offset the RNG repeats itself.
    seed_offset = 0
    if cfg.resume:
        ckpt_meta = torch.load(cfg.resume, mmap=True, map_location="cpu", weights_only=False)
        seed_offset = int(ckpt_meta.get("epoch", 0)) + int(ckpt_meta.get("global_step", 0))
        del ckpt_meta
    if cfg.seed is not None:
        seed_everything(dist_manager.group_rank["world"] + seed_offset + int(cfg.seed))

    trainer_cfg = dict(cfg.trainer or {})

    _EXPECTED_DTENSOR_TRAINER = {"devices": 1, "num_nodes": 1}
    for key, expected in _EXPECTED_DTENSOR_TRAINER.items():
        val = trainer_cfg.get(key)
        if val is not None and val != expected:
            raise ValueError(
                f"trainer.{key}={val!r} is incompatible with DTensor context-parallel training "
                f"(expected {expected!r}). The distributed topology is managed by "
                f"DistributedManager via parallel_size, not by Lightning. "
                f"Set parallel_size.size_dp and parallel_size.size_cp instead, "
                f"and use trainer.{key}={expected!r} or omit it."
            )
        trainer_cfg[key] = expected

    num_workers = getattr(getattr(cfg, "data", None), "num_workers", 0)
    if num_workers != 0:
        raise ValueError(
            f"data.num_workers={num_workers} is not supported in DTensor context-parallel training. "
            f"Only num_workers=0 is supported because the DTensor data workflow requires "
            f"main-process collation for distributed tensor construction. "
            f"Set data.num_workers=0 in your config."
        )

    wandb_cfg = cfg.wandb
    if cfg.debug:
        wandb_cfg = None

    data_module = _create_distributed_data_module(cfg.data, dist_manager)
    model_module = _create_distributed_model(cfg, dist_manager)

    model_module.apply(SetTriAttnBackend(cfg.triattn_backend))
    model_module.apply(SetAttnPairBiasBackend(cfg.sdpa_with_bias_backend))
    model_module.apply(SetAttnPairBiasShardwiseBackend(cfg.sdpa_with_bias_shardwise_backend))
    if offload_actv_ckpt_cfg is not None:
        model_module.apply(OffloadActvCkptToCPU(set(offload_actv_ckpt_cfg)))

    if getattr(model_module, "confidence_prediction", False):
        model_module.confidence_prediction = False
        warnings.warn("Confidence prediction is not supported in distributed training mode")

    steering_args = getattr(model_module, "steering_args", None)
    if steering_args is not None:
        for attr in ("fk_steering", "guidance_update"):
            if getattr(steering_args, attr, False):
                setattr(steering_args, attr, False)
        warnings.warn("Steering potentials are not supported in distributed training mode")

    callbacks: list[Any] = []
    if not cfg.disable_checkpoint:
        # Boltz-2 checkpoint defaults; overridable via the ``checkpoint`` config key.
        checkpoint_cfg = dict(cfg.checkpoint or {})
        checkpoint_cfg.setdefault("filename", "{epoch:02d}-{step:05d}")
        checkpoint_cfg.setdefault("monitor", "val/lddt")
        checkpoint_cfg.setdefault("save_top_k", cfg.save_top_k)
        checkpoint_cfg.setdefault("save_last", True)
        checkpoint_cfg.setdefault("save_on_train_epoch_end", True)
        checkpoint_cfg.setdefault("mode", "max")
        checkpoint_cfg.setdefault("every_n_epochs", 1)
        callbacks.append(ModelCheckpoint(dirpath=cfg.output, **checkpoint_cfg))

    if cuda_memory_profile_cfg is not None and cuda_memory_profile_cfg.get("output_path_prefix") is not None:
        output_path = cuda_memory_profile_cfg.output_path_prefix + f"_rank{dist_manager.group_rank['world']}.pickle"
        memory_profile_kwargs = {k: v for k, v in cuda_memory_profile_cfg.items() if k != "output_path_prefix"}
        callbacks.append(CUDAMemoryProfile(output_path=output_path, **memory_profile_kwargs))

    loggers: list[Any] = []
    if wandb_cfg:
        wandb_id = wandb_cfg.get("id")
        wandb_resume = "allow" if wandb_id else None
        wdb_logger = WandbLogger(
            name=wandb_cfg["name"],
            group=wandb_cfg["name"],
            save_dir=cfg.output,
            project=wandb_cfg["project"],
            entity=wandb_cfg["entity"],
            id=wandb_id,
            resume=wandb_resume,
            log_model=False,
        )
        loggers.append(wdb_logger)

        @rank_zero_only
        def save_config_to_wandb() -> None:
            config_out = Path(wdb_logger.experiment.dir) / "run.yaml"
            with config_out.open("w") as file_handle:
                OmegaConf.save(config_dict, file_handle)
            wdb_logger.experiment.save(str(config_out))

        save_config_to_wandb()

    strategy = BoltzContextParallelStrategy(dist_manager=dist_manager)

    if cfg.precision not in PRECISION_TO_LIGHTNING:
        raise ValueError(f"Precision {cfg.precision} is not supported")
    if trainer_cfg.get("precision") is not None:
        raise ValueError(
            "Set precision in the top-level config, not inside trainer. "
            "The trainer.precision key is superseded by the top-level precision setting."
        )
    trainer_cfg["precision"] = PRECISION_TO_LIGHTNING[cfg.precision]

    trainer_kwargs = dict(
        default_root_dir=cfg.output,
        strategy=strategy,
        callbacks=callbacks,
        logger=loggers,
        enable_checkpointing=not cfg.disable_checkpoint,
        reload_dataloaders_every_n_epochs=1,
        use_distributed_sampler=False,  # distributed data module handles its own sharding
        **trainer_cfg,
    )

    if _one_logger_available:
        # Compute global batch size for OneLogger compliance.
        batch_size = getattr(getattr(data_module, "cfg", None), "batch_size", 1)
        one_logger_config: dict[str, Any] = {
            "global_batch_size": dist_manager.group["dp"].size() * batch_size,
        }
        if wandb_cfg:
            one_logger_config.update(
                {
                    "name": wandb_cfg["name"],
                    "group": wandb_cfg["name"],
                    "save_dir": cfg.output,
                    "project": wandb_cfg["project"],
                    "entity": wandb_cfg["entity"],
                    "log_model": False,
                }
            )
        HookedTrainer, one_logger_callback = hook_trainer_cls(pl.Trainer, callback_config=one_logger_config)
        callbacks.append(one_logger_callback)
        trainer = HookedTrainer(**trainer_kwargs)
    else:
        trainer = pl.Trainer(**trainer_kwargs)

    # Suppress expected Lightning warnings in CP mode.
    warnings.filterwarnings(
        "ignore",
        message="It is recommended to use .* when logging on epoch level in "
        "distributed setting to accumulate the metric across devices",
    )
    warnings.filterwarnings(
        "ignore",
        message="The .* does not have many workers which may be a bottleneck. "
        "Consider increasing the value of the `num_workers` .* to improve performance.",
    )

    try:
        with setup_tf32_env(cfg.precision):
            if cfg.validation_only:
                trainer.validate(model_module, datamodule=data_module, ckpt_path=cfg.resume)
            else:
                trainer.fit(model_module, datamodule=data_module, ckpt_path=cfg.resume)
    finally:
        _cleanup_distributed()


if __name__ == "__main__":
    train(sys.argv[1], sys.argv[2:])
