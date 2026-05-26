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

"""Distributed DTensor CP stop/go (checkpoint resume) parity tests.

The primary stop/go test routes through the real ``train()`` entrypoint,
exercising the full save/resume cycle including
:class:`BoltzContextParallelStrategy` checkpoint conversion, checkpoint
callback defaults, and the ``cfg.resume`` auto-resume path.

Cross-mode tests (serial <-> distributed) use direct ``Trainer.fit()`` calls
because the serial leg cannot use ``train()`` (which always creates a
``BoltzContextParallelStrategy``).
"""

from pathlib import Path
from typing import Any

import pytest
import pytorch_lightning as pl
import torch
from omegaconf import OmegaConf
from torch.distributed.checkpoint.state_dict import get_optimizer_state_dict
from torch.distributed.tensor import DTensor, distribute_tensor

import boltz.distributed.train as train_module
from boltz.distributed.lightning_strategy import BoltzContextParallelStrategy
from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.models.boltz2 import Boltz2 as Boltz2Distributed
from boltz.distributed.model.models.boltz2 import _PlaceholderModule
from boltz.distributed.model.modules.utils import has_dtensors
from boltz.model.models.boltz2 import Boltz2 as SerialBoltz2
from boltz.testing.utils import spawn_multiprocessing

from .dtensor_train_harness import (
    DistogramTrainDataModule,
    TinyDistogramCPModel,
    TinyDistogramCPModelWithEMA,
    TinyDistogramSerialModel,
    create_initial_serial_state_dict,
    create_train_dataloader,
)
from .model.models.test_dtensor_boltz2 import _prepare_serial_model


def _to_local(t: torch.Tensor) -> torch.Tensor:
    """Unwrap DTensor to local tensor for comparison."""
    return t.to_local() if isinstance(t, DTensor) else t


def _assert_optimizer_states_match(
    opt_a: torch.optim.Optimizer,
    opt_b: torch.optim.Optimizer,
    label: str,
) -> None:
    """Assert that two optimizers have matching state (exp_avg, exp_avg_sq, step).

    This compares the per-parameter state buffers that Adam maintains.
    If these diverge after resume, future gradient updates will be wrong.
    """
    state_a = opt_a.state_dict()["state"]
    state_b = opt_b.state_dict()["state"]
    assert state_a.keys() == state_b.keys(), f"{label}: optimizer state key mismatch"
    for param_idx in state_a:
        for buf_key in state_a[param_idx]:
            val_a = state_a[param_idx][buf_key]
            val_b = state_b[param_idx][buf_key]
            if isinstance(val_a, torch.Tensor):
                torch.testing.assert_close(
                    _to_local(val_a).cpu(),
                    _to_local(val_b).cpu(),
                    msg=lambda msg, k=buf_key, p=param_idx: (
                        f"{label}: optimizer state mismatch for param {p}, key '{k}'\n{msg}"
                    ),
                )
            else:
                assert val_a == val_b, f"{label}: optimizer scalar mismatch param {param_idx}, key '{buf_key}'"


# ---------------------------------------------------------------------------
# Distogram stop/go via train.py entrypoint
# ---------------------------------------------------------------------------


def _write_distogram_config(
    *,
    config_path: Path,
    output_dir: Path,
    size_dp: int,
    size_cp: int,
    accelerator: str = "cpu",
    max_epochs: int = 2,
    limit_train_batches: int = 2,
    resume: str | None = None,
    weights_seed: int = 37,
    data_seed: int = 53,
    token_z: int = 16,
    num_bins: int = 8,
    num_distograms: int = 2,
    num_conformers: int = 2,
    seq_len: int = 12,
    num_samples: int = 2,
    learning_rate: float = 1e-2,
    ema_decay: float = 0.999,
) -> None:
    """Write a YAML config for a distogram stop/go ``train()`` run."""
    config: dict[str, Any] = {
        "data": {
            "seq_len": seq_len,
            "token_z": token_z,
            "num_bins": num_bins,
            "num_conformers": num_conformers,
            "num_samples": num_samples,
            "seed": data_seed,
        },
        "model": {
            "token_z": token_z,
            "num_bins": num_bins,
            "num_distograms": num_distograms,
            "num_conformers": num_conformers,
            "weights_seed": weights_seed,
            "learning_rate": learning_rate,
            "ema_decay": ema_decay,
        },
        "output": str(output_dir),
        "trainer": {
            "accelerator": accelerator,
            "devices": 1,
            "max_epochs": max_epochs,
            "limit_train_batches": limit_train_batches,
            "enable_progress_bar": False,
            "enable_model_summary": False,
            "num_sanity_val_steps": 0,
        },
        "parallel_size": {"size_dp": size_dp, "size_cp": size_cp},
        "precision": "FP32",
        "find_unused_parameters": False,
        "save_top_k": -1,
        "disable_checkpoint": False,
        "debug": False,
        "validation_only": False,
        "seed": 11,
        "checkpoint": {
            "monitor": None,
            "save_last": True,
            "every_n_epochs": 1,
        },
    }
    if resume is not None:
        config["resume"] = resume
    config_path.parent.mkdir(parents=True, exist_ok=True)
    OmegaConf.save(OmegaConf.create(config), config_path)


def _instantiate_distogram_config(config_dict: Any) -> dict[str, Any]:
    """Replace ``hydra.utils.instantiate`` for distogram stop/go tests.

    Returns the config as a plain dict; ``_create_distogram_distributed_model``
    handles model creation from config + dist_manager.
    """
    cfg = OmegaConf.to_container(config_dict, resolve=True)
    assert isinstance(cfg, dict)
    return cfg


def _create_distogram_distributed_model(
    cfg: Any,
    dist_manager: DistributedManager,
) -> TinyDistogramCPModelWithEMA:
    """Monkeypatched ``_create_distributed_model`` for distogram tests."""
    model_cfg = cfg.model
    serial_state = create_initial_serial_state_dict(
        token_z=model_cfg["token_z"],
        num_bins=model_cfg["num_bins"],
        num_distograms=model_cfg["num_distograms"],
        seed=model_cfg["weights_seed"],
    )
    return TinyDistogramCPModelWithEMA(
        dist_manager=dist_manager,
        token_z=model_cfg["token_z"],
        num_bins=model_cfg["num_bins"],
        num_distograms=model_cfg["num_distograms"],
        num_conformers=model_cfg["num_conformers"],
        serial_state_dict=serial_state,
        learning_rate=model_cfg["learning_rate"],
        ema_decay=model_cfg.get("ema_decay", 0.999),
    )


def _create_distogram_distributed_data_module(
    data_config: Any,
    dist_manager: DistributedManager,
) -> DistogramTrainDataModule:
    """Monkeypatched ``_create_distributed_data_module`` for distogram tests."""
    return DistogramTrainDataModule(
        seq_len=data_config["seq_len"],
        token_z=data_config["token_z"],
        num_bins=data_config["num_bins"],
        num_conformers=data_config["num_conformers"],
        num_samples=data_config["num_samples"],
        seed=data_config["seed"],
        dp_size=dist_manager.group["dp"].size(),
    )


def _assert_checkpoint_optimizer_states_match(
    opt_a: dict[str, Any],
    opt_b: dict[str, Any],
    label: str,
) -> None:
    """Assert optimizer state dicts from two checkpoints match.

    Handles both FQN (string) and legacy integer keys transparently.
    """
    state_a = opt_a["state"]
    state_b = opt_b["state"]
    assert state_a.keys() == state_b.keys(), f"{label}: optimizer state key mismatch"
    for param_key in state_a:
        for buf_key in state_a[param_key]:
            val_a = state_a[param_key][buf_key]
            val_b = state_b[param_key][buf_key]
            if isinstance(val_a, torch.Tensor):
                torch.testing.assert_close(
                    val_a,
                    val_b,
                    msg=lambda msg, k=buf_key, p=param_key: (
                        f"{label}: optimizer state mismatch for param {p}, key '{k}'\n{msg}"
                    ),
                )
            else:
                assert val_a == val_b, f"{label}: optimizer scalar mismatch param {param_key}, key '{buf_key}'"


def _parallel_assert_dtensor_stop_and_go_ema(rank: int, payload: tuple[Any, ...]) -> None:
    """Verify stop/go parity through the real ``train()`` entrypoint.

    Runs three ``train()`` calls:
    1. Continuous baseline — 2 epochs, checkpointing enabled
    2. Stop/go stage 1 — 1 epoch, checkpointing enabled
    3. Stop/go stage 2 — resume from checkpoint, complete to epoch 2

    Compares the final ``last.ckpt`` files for exact parity: model weights,
    optimizer state, EMA shadow weights, epoch, and global step.  This
    validates the full save→resume cycle through ``train.py``, including
    strategy checkpoint conversion and ``cfg.resume``.
    """
    (
        env_per_rank,
        continuous_config_path,
        stage1_config_path,
        stage2_config_path,
        continuous_dir,
        stopgo_dir,
    ) = payload
    continuous_dir = Path(continuous_dir)
    stopgo_dir = Path(stopgo_dir)

    monkeypatch = pytest.MonkeyPatch()
    for key, value in env_per_rank.items():
        monkeypatch.setenv(key, f"{rank}" if value == "<INPUT_RANK>" else value)

    monkeypatch.setattr(train_module.hydra.utils, "instantiate", _instantiate_distogram_config)
    monkeypatch.setattr(train_module, "_create_distributed_model", _create_distogram_distributed_model)
    monkeypatch.setattr(train_module, "_create_distributed_data_module", _create_distogram_distributed_data_module)
    # Suppress per-call cleanup so process groups survive across the 3
    # sequential train() calls; the test's own finally block cleans up.
    monkeypatch.setattr(train_module, "_cleanup_distributed", lambda: None)
    DistributedManager._state = {}

    try:
        # ---- Continuous baseline: 2 epochs in one run. ----
        # _create_dist_manager initializes process groups on the first call;
        # subsequent calls reuse the existing DistributedManager singleton.
        train_module.train(str(continuous_config_path), [])

        # ---- Stop/go stage 1: 1 epoch, checkpoint produced. ----
        train_module.train(str(stage1_config_path), [])
        ckpt_path = stopgo_dir / "last.ckpt"
        assert ckpt_path.exists(), f"Rank {rank}: stage 1 checkpoint not found at {ckpt_path}"

        # Sanity: stage-1 (1 epoch) must differ from the continuous run
        # (2 epochs).  Guards against a vacuous test where both checkpoints
        # are identical before resume even happens.
        continuous_ckpt_early = torch.load(continuous_dir / "last.ckpt", map_location="cpu", weights_only=False)
        stage1_ckpt = torch.load(ckpt_path, map_location="cpu", weights_only=False)
        assert (
            continuous_ckpt_early["epoch"] != stage1_ckpt["epoch"]
            or continuous_ckpt_early["global_step"] != stage1_ckpt["global_step"]
        ), (
            "Stage-1 checkpoint should differ from the 2-epoch continuous checkpoint "
            f"(both have epoch={stage1_ckpt['epoch']}, step={stage1_ckpt['global_step']})"
        )
        weights_differ = any(
            not torch.equal(continuous_ckpt_early["state_dict"][k], stage1_ckpt["state_dict"][k])
            for k in stage1_ckpt["state_dict"]
        )
        assert weights_differ, "Stage-1 weights should differ from 2-epoch continuous weights"

        # ---- Stop/go stage 2: resume from checkpoint to epoch 2. ----
        train_module.train(str(stage2_config_path), [])

        # ---- Compare checkpoint files for parity. ----
        continuous_ckpt = torch.load(continuous_dir / "last.ckpt", map_location="cpu", weights_only=False)
        stopgo_ckpt = torch.load(stopgo_dir / "last.ckpt", map_location="cpu", weights_only=False)

        # Checkpoints should contain only plain tensors (strategy strips DTensors).
        assert not has_dtensors(continuous_ckpt["state_dict"]), "Continuous checkpoint has DTensors"
        assert not has_dtensors(stopgo_ckpt["state_dict"]), "Stop/go checkpoint has DTensors"

        # 1) Model weights must match.
        assert continuous_ckpt["state_dict"].keys() == stopgo_ckpt["state_dict"].keys(), "state_dict key mismatch"
        for key in continuous_ckpt["state_dict"]:
            torch.testing.assert_close(
                continuous_ckpt["state_dict"][key],
                stopgo_ckpt["state_dict"][key],
                msg=lambda msg, k=key: f"Stop/go weight mismatch for {k} on rank {rank}\n{msg}",
            )

        # 2) Epoch and global step must match.
        assert (
            continuous_ckpt["epoch"] == stopgo_ckpt["epoch"]
        ), f"Epoch mismatch: continuous={continuous_ckpt['epoch']}, stopgo={stopgo_ckpt['epoch']}"
        assert (
            continuous_ckpt["global_step"] == stopgo_ckpt["global_step"]
        ), f"Step mismatch: continuous={continuous_ckpt['global_step']}, stopgo={stopgo_ckpt['global_step']}"

        # 3) Optimizer state (Adam exp_avg / exp_avg_sq / step) must match.
        _assert_checkpoint_optimizer_states_match(
            continuous_ckpt["optimizer_states"][0],
            stopgo_ckpt["optimizer_states"][0],
            label=f"rank {rank}",
        )

        # 3b) Optimizer state keys must be FQN strings (not legacy integers).
        opt_state_keys = list(continuous_ckpt["optimizer_states"][0]["state"].keys())
        assert opt_state_keys, f"Rank {rank}: optimizer state is empty"
        assert all(isinstance(k, str) for k in opt_state_keys), (
            f"Rank {rank}: optimizer state keys should be FQN strings, "
            f"got {[type(k).__name__ for k in opt_state_keys[:3]]}"
        )

        # 4) EMA shadow weights and step counter must match.
        assert "ema" in continuous_ckpt, "Continuous checkpoint missing EMA state"
        assert "ema" in stopgo_ckpt, "Stop/go checkpoint missing EMA state"
        assert continuous_ckpt["ema"]["cur_step"] == stopgo_ckpt["ema"]["cur_step"], (
            f"EMA step mismatch: continuous={continuous_ckpt['ema']['cur_step']}, "
            f"stopgo={stopgo_ckpt['ema']['cur_step']}"
        )
        for key in continuous_ckpt["ema"]["ema_weights"]:
            torch.testing.assert_close(
                continuous_ckpt["ema"]["ema_weights"][key],
                stopgo_ckpt["ema"]["ema_weights"][key],
                msg=lambda msg, k=key: f"EMA weight mismatch for {k} on rank {rank}\n{msg}",
            )
    finally:
        DistributedManager.cleanup()
        DistributedManager._state = {}
        monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env",
    [
        ((2, (3, 3)), True, "cpu", "ENV"),
        ((2, (1, 1)), True, "cuda", "ENV"),
        ((1, (2, 2)), True, "cuda", "ENV"),
    ],
    indirect=("setup_env",),
    ids=["cpu-dp2-cp3x3", "cuda-dp2-cp1x1", "cuda-dp1-cp2x2"],
)
def test_stop_and_go_via_train_entrypoint(setup_env, tmp_path):
    """Goals: checkpoint resume parity through the real ``train()`` entrypoint.

    - Continuous 2-epoch run matches stop-at-epoch-1 + resume-to-epoch-2
    - Model weights, optimizer state, EMA shadow weights all match exactly
    - Epoch and global_step counters match
    - Validates ``BoltzContextParallelStrategy`` checkpoint conversion roundtrip
    - Validates ``cfg.resume`` auto-resume path in ``train.py``
    """
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env
    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("CUDA not available")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"Need {world_size} GPUs, have {torch.cuda.device_count()}")

    ema_decay = 0.999
    size_dp = int(grid_group_sizes["dp"])
    cp_group = grid_group_sizes["cp"]
    size_cp = int(cp_group[0] * cp_group[1]) if isinstance(cp_group, tuple) else int(cp_group)
    accelerator = "gpu" if device_type == "cuda" else "cpu"

    continuous_dir = tmp_path / "continuous"
    stopgo_dir = tmp_path / "stopgo"

    common_kwargs: dict[str, Any] = {
        "size_dp": size_dp,
        "size_cp": size_cp,
        "accelerator": accelerator,
        "ema_decay": ema_decay,
    }

    # Continuous baseline: 2 epochs.
    continuous_config = continuous_dir / "config.yaml"
    _write_distogram_config(config_path=continuous_config, output_dir=continuous_dir, max_epochs=2, **common_kwargs)

    # Stage 1: 1 epoch with checkpoint.
    stage1_config = stopgo_dir / "config_stage1.yaml"
    _write_distogram_config(config_path=stage1_config, output_dir=stopgo_dir, max_epochs=1, **common_kwargs)

    # Stage 2: resume to epoch 2.
    stage2_config = stopgo_dir / "config_stage2.yaml"
    _write_distogram_config(
        config_path=stage2_config,
        output_dir=stopgo_dir,
        max_epochs=2,
        resume=str(stopgo_dir / "last.ckpt"),
        **common_kwargs,
    )

    payload = (
        env_per_rank,
        str(continuous_config),
        str(stage1_config),
        str(stage2_config),
        str(continuous_dir),
        str(stopgo_dir),
    )
    spawn_multiprocessing(_parallel_assert_dtensor_stop_and_go_ema, world_size, payload)


# ---------------------------------------------------------------------------
# Cross-mode stop/go: serial <-> distributed checkpoint interop (both dirs)
# ---------------------------------------------------------------------------


def _parallel_assert_cross_mode_stop_and_go(rank: int, payload: tuple[Any, ...]) -> None:
    """Both cross-mode directions in one worker.

    Direction 1 (serial → distributed):
        1-epoch serial train → checkpoint → resume as distributed for epoch 2.
        Compared against a 2-epoch continuous distributed baseline.

    Direction 2 (distributed → serial):
        1-epoch distributed train → checkpoint → resume as serial for epoch 2.
        Compared against a 2-epoch continuous serial baseline (rank 0 only).
    """
    grid_group_sizes, device_type, backend, env_per_rank, output_dir = payload
    output_dir = Path(output_dir)

    token_z, num_bins, num_distograms, num_conformers, seq_len = 16, 8, 2, 2, 12

    monkeypatch = pytest.MonkeyPatch()
    for key, value in env_per_rank.items():
        monkeypatch.setenv(key, f"{rank}" if value == "<INPUT_RANK>" else value)

    try:
        DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
        manager = DistributedManager()
        dp_size = manager.group["dp"].size()

        serial_state = create_initial_serial_state_dict(
            token_z=token_z,
            num_bins=num_bins,
            num_distograms=num_distograms,
            seed=41,
        )

        model_kwargs: dict[str, Any] = {
            "token_z": token_z,
            "num_bins": num_bins,
            "num_distograms": num_distograms,
            "num_conformers": num_conformers,
            "serial_state_dict": serial_state,
            "learning_rate": 1e-2,
        }

        def _make_dl():
            return create_train_dataloader(
                seq_len=seq_len,
                token_z=token_z,
                num_bins=num_bins,
                num_conformers=num_conformers,
                num_samples=2,
                seed=53,
                dp_size=dp_size,
            )

        def _trainer(root, *, strategy=None, ckpt_cb=None, epochs=2):
            kw: dict[str, Any] = {
                "default_root_dir": str(root),
                "accelerator": "cpu" if device_type == "cpu" else "gpu",
                "devices": 1,
                "max_epochs": epochs,
                "limit_train_batches": 2,
                "logger": False,
                "enable_progress_bar": False,
                "enable_model_summary": False,
            }
            if strategy is not None:
                kw["strategy"] = strategy
                kw["use_distributed_sampler"] = False
            if ckpt_cb:
                kw["callbacks"] = [ckpt_cb]
                kw["enable_checkpointing"] = True
            else:
                kw["callbacks"] = []
                kw["enable_checkpointing"] = False
            return pl.Trainer(**kw)

        def _ckpt_cb(dirpath):
            return pl.callbacks.ModelCheckpoint(
                dirpath=str(dirpath),
                filename="epoch-{epoch:02d}",
                every_n_epochs=1,
                save_top_k=-1,
                save_last=True,
            )

        def _check_parity(state_cont, state_res, trainer_cont, trainer_res, model_cont, model_res, label):
            for key in state_cont:
                torch.testing.assert_close(
                    state_res[key], state_cont[key], msg=lambda msg, k=key: f"{label} mismatch for {k}\n{msg}"
                )
            assert trainer_cont.current_epoch == trainer_res.current_epoch
            assert trainer_cont.global_step == trainer_res.global_step
            _assert_optimizer_states_match(trainer_cont.optimizers[0], trainer_res.optimizers[0], label=label)
            e2_c = {k: v for k, v in model_cont._loss_log.items() if k[0] == 1}
            e2_r = {k: v for k, v in model_res._loss_log.items() if k[0] == 1}
            assert e2_c.keys() == e2_r.keys()
            for k in sorted(e2_c):
                assert e2_c[k] == pytest.approx(e2_r[k]), f"{label}: loss mismatch at {k}"

        # ================================================================
        # Direction 1: serial → distributed
        # ================================================================
        s2d = output_dir / "s2d"

        # 2-epoch continuous distributed baseline.
        s2d_cont_model = TinyDistogramCPModel(dist_manager=manager, **model_kwargs)
        s2d_cont_strat = BoltzContextParallelStrategy(dist_manager=manager)
        s2d_cont_tr = _trainer(s2d / "cont", strategy=s2d_cont_strat)
        s2d_cont_tr.fit(s2d_cont_model, train_dataloaders=_make_dl())
        state_s2d_cont = s2d_cont_strat.lightning_module_state_dict()

        # Stage 1: 1-epoch serial (rank 0 only).
        s2d_ckpt = s2d / "serial" / "last.ckpt"
        if manager.rank == 0:
            s2d_s1 = TinyDistogramSerialModel(**model_kwargs)
            s2d_s1_tr = _trainer(s2d / "serial", ckpt_cb=_ckpt_cb(s2d / "serial"), epochs=1)
            s2d_s1_tr.fit(s2d_s1, train_dataloaders=_make_dl())
            assert s2d_ckpt.exists()
            # Guard: serial stage-1 must differ from initial (training was not a no-op).
            for key, init_val in serial_state.items():
                prefixed = f"distogram_module.{key}"
                assert not torch.equal(
                    s2d_s1.state_dict()[prefixed].cpu(), init_val
                ), f"Serial stage-1 weight '{key}' unchanged after 1 epoch"
            # Verify serial checkpoint uses integer keys (standard Lightning).
            s2d_s1_ckpt_data = torch.load(s2d_ckpt, map_location="cpu", weights_only=False)
            s2d_s1_opt_keys = list(s2d_s1_ckpt_data["optimizer_states"][0]["state"].keys())
            assert s2d_s1_opt_keys, "Serial checkpoint optimizer state is empty"
            assert all(isinstance(k, int) for k in s2d_s1_opt_keys), (
                f"Serial checkpoint should have integer optimizer keys, "
                f"got types {[type(k).__name__ for k in s2d_s1_opt_keys[:3]]}"
            )
        torch.distributed.barrier()

        # Stage 2: resume as distributed (loads serial int-key checkpoint via
        # the legacy path in load_optimizer_state_dict).
        s2d_s2 = TinyDistogramCPModel(dist_manager=manager, **model_kwargs)
        s2d_s2_strat = BoltzContextParallelStrategy(dist_manager=manager)
        s2d_s2_tr = _trainer(s2d / "resume", strategy=s2d_s2_strat)
        s2d_s2_tr.fit(s2d_s2, train_dataloaders=_make_dl(), ckpt_path=str(s2d_ckpt))
        state_s2d_res = s2d_s2_strat.lightning_module_state_dict()

        _check_parity(
            state_s2d_cont,
            state_s2d_res,
            s2d_cont_tr,
            s2d_s2_tr,
            s2d_cont_model,
            s2d_s2,
            f"rank {rank} serial→distributed",
        )

        # ================================================================
        # Direction 2: distributed → serial
        # ================================================================
        d2s = output_dir / "d2s"

        # 2-epoch continuous serial baseline (rank 0 only).
        if manager.rank == 0:
            d2s_cont_model = TinyDistogramSerialModel(**model_kwargs)
            d2s_cont_tr = _trainer(d2s / "cont", epochs=2)
            d2s_cont_tr.fit(d2s_cont_model, train_dataloaders=_make_dl())
            state_d2s_cont = {k: v.detach().cpu().clone() for k, v in d2s_cont_model.state_dict().items()}
        torch.distributed.barrier()

        # Stage 1: 1-epoch distributed with checkpoint.
        d2s_s1 = TinyDistogramCPModel(dist_manager=manager, **model_kwargs)
        d2s_s1_strat = BoltzContextParallelStrategy(dist_manager=manager)
        d2s_s1_tr = _trainer(d2s / "dist", strategy=d2s_s1_strat, ckpt_cb=_ckpt_cb(d2s / "dist"), epochs=1)
        d2s_s1_tr.fit(d2s_s1, train_dataloaders=_make_dl())
        # Guard: distributed stage-1 must differ from initial (training was not a no-op).
        d2s_s1_sd = d2s_s1_strat.lightning_module_state_dict()
        for key, init_val in serial_state.items():
            prefixed = f"distogram_module.{key}"
            assert not torch.equal(
                d2s_s1_sd[prefixed].cpu(), init_val
            ), f"Distributed stage-1 weight '{key}' unchanged after 1 epoch"
        # Barrier after all-ranks distributed training, before rank-0-only
        # assertions.  Placing it here avoids deadlock: if rank 0's assertions
        # fail below, other ranks have already passed this sync point and
        # exit cleanly instead of hanging in an NCCL wait.
        torch.distributed.barrier()

        d2s_ckpt = d2s / "dist" / "last.ckpt"
        if manager.rank == 0:
            assert d2s_ckpt.exists()
            ckpt = torch.load(d2s_ckpt, map_location="cpu", weights_only=False)
            assert not has_dtensors(ckpt["state_dict"]), "Distributed ckpt should be plain tensors"
            # Verify distributed checkpoint uses FQN string keys.
            d2s_opt_keys = list(ckpt["optimizer_states"][0]["state"].keys())
            assert d2s_opt_keys, "Distributed checkpoint optimizer state is empty"
            assert all(isinstance(k, str) for k in d2s_opt_keys), (
                f"Distributed checkpoint should have FQN string optimizer keys, "
                f"got types {[type(k).__name__ for k in d2s_opt_keys[:3]]}"
            )
            # Verify FQN keys correspond to actual model parameter names.
            expected_param_names = [n for n, _ in d2s_s1.named_parameters()]
            assert sorted(d2s_opt_keys) == sorted(expected_param_names), (
                f"FQN optimizer keys don't match model parameters.\n"
                f"  Optimizer keys: {sorted(d2s_opt_keys)}\n"
                f"  Model params:   {sorted(expected_param_names)}"
            )
            # Verify param_groups also use FQN keys (not integers).
            pg_params = ckpt["optimizer_states"][0]["param_groups"][0]["params"]
            assert all(isinstance(p, str) for p in pg_params), (
                f"Distributed checkpoint param_groups should use FQN strings, "
                f"got types {[type(p).__name__ for p in pg_params[:3]]}"
            )

        # Stage 2: resume as serial (rank 0 only).
        # The serial model uses Lightning's default load path, which calls
        # optimizer.load_state_dict() — this handles FQN keys transparently
        # via positional mapping in param_groups.
        if manager.rank == 0:
            d2s_s2 = TinyDistogramSerialModel(**model_kwargs)
            d2s_s2_tr = _trainer(d2s / "resume", epochs=2)
            d2s_s2_tr.fit(d2s_s2, train_dataloaders=_make_dl(), ckpt_path=str(d2s_ckpt))
            state_d2s_res = {k: v.detach().cpu().clone() for k, v in d2s_s2.state_dict().items()}
            _check_parity(
                state_d2s_cont,
                state_d2s_res,
                d2s_cont_tr,
                d2s_s2_tr,
                d2s_cont_model,
                d2s_s2,
                "rank 0 distributed→serial",
            )
    finally:
        DistributedManager.cleanup()
        DistributedManager._state = {}
        monkeypatch.undo()


# ---------------------------------------------------------------------------
# Optimizer parameter ordering: serial vs distributed
# ---------------------------------------------------------------------------


def _parallel_assert_optimizer_param_ordering(rank: int, payload: tuple[Any, ...]) -> None:
    """Assert that optimizer parameter ordering matches between serial and distributed models.

    PyTorch optimizer state_dict uses integer keys derived from the iteration
    order of ``model.parameters()`` (which mirrors ``model.named_parameters()``).
    If serial and distributed models yield parameters in different orders, cross-
    topology checkpoint resume silently applies optimizer state (exp_avg, etc.)
    to the wrong parameters.  This test catches that.

    Also verifies that ``get_optimizer_state_dict`` produces FQN keys matching
    ``named_parameters()`` — the mechanism used by ``BoltzContextParallelStrategy``
    for portable, name-keyed optimizer checkpoints.
    """
    grid_group_sizes, device_type, backend, env_per_rank = payload

    monkeypatch = pytest.MonkeyPatch()
    for key, value in env_per_rank.items():
        monkeypatch.setenv(key, f"{rank}" if value == "<INPUT_RANK>" else value)

    try:
        DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
        manager = DistributedManager()

        serial_state = create_initial_serial_state_dict(token_z=16, num_bins=8, num_distograms=2, seed=41)
        model_kwargs = {
            "token_z": 16,
            "num_bins": 8,
            "num_distograms": 2,
            "num_conformers": 2,
            "serial_state_dict": serial_state,
            "learning_rate": 1e-2,
        }

        cp_model = TinyDistogramCPModel(dist_manager=manager, **model_kwargs)
        serial_model = TinyDistogramSerialModel(**model_kwargs)

        # The optimizer receives parameters in iteration order of model.parameters().
        # named_parameters() yields (name, param) in the same order, so comparing
        # the name lists verifies that integer optimizer state keys will align.
        cp_names = [name for name, _ in cp_model.named_parameters()]
        serial_names = [name for name, _ in serial_model.named_parameters()]

        assert cp_names == serial_names, (
            f"Parameter ordering mismatch between serial and distributed models.\n"
            f"  Serial:      {serial_names}\n"
            f"  Distributed: {cp_names}\n"
            f"Optimizer state checkpoint resume will silently apply state to wrong parameters."
        )

        # Also verify that the optimizers themselves see the same number of
        # parameters in their param_groups (guards against extra params from
        # distributed wrappers).
        cp_opt = cp_model.configure_optimizers()
        serial_opt = serial_model.configure_optimizers()
        cp_param_count = sum(len(g["params"]) for g in cp_opt.param_groups)
        serial_param_count = sum(len(g["params"]) for g in serial_opt.param_groups)
        assert cp_param_count == serial_param_count, (
            f"Optimizer param count mismatch: serial has {serial_param_count} params, "
            f"distributed has {cp_param_count} params"
        )

        # Verify get_optimizer_state_dict produces FQN keys matching
        # named_parameters().  This is the save-side mechanism used by
        # BoltzContextParallelStrategy.optimizer_state().

        # Run one optimizer step to populate state buffers (exp_avg, etc.).
        # Use a dummy gradient to avoid going through training_step which
        # requires Trainer integration.
        for p in cp_model.parameters():
            p.grad = torch.randn_like(p.to_local() if isinstance(p, DTensor) else p)
            if isinstance(p, DTensor):
                p.grad = distribute_tensor(p.grad, device_mesh=p.device_mesh, placements=p.placements)
        cp_opt.step()

        fqn_sd = get_optimizer_state_dict(cp_model, cp_opt)
        fqn_state_keys = sorted(fqn_sd["state"].keys())
        fqn_pg_params = sorted(fqn_sd["param_groups"][0]["params"])

        assert all(isinstance(k, str) for k in fqn_state_keys), (
            f"get_optimizer_state_dict should return FQN string keys, "
            f"got types {[type(k).__name__ for k in fqn_state_keys[:3]]}"
        )
        assert fqn_state_keys == sorted(cp_names), (
            f"FQN optimizer state keys don't match named_parameters().\n"
            f"  FQN keys:          {fqn_state_keys}\n"
            f"  named_parameters:  {sorted(cp_names)}"
        )
        assert fqn_pg_params == sorted(cp_names), (
            f"FQN param_group params don't match named_parameters().\n"
            f"  param_groups:      {fqn_pg_params}\n"
            f"  named_parameters:  {sorted(cp_names)}"
        )

        # Verify FQN keys match between serial and distributed models.
        for p in serial_model.parameters():
            p.grad = torch.randn_like(p)
        serial_opt.step()

        serial_fqn_sd = get_optimizer_state_dict(serial_model, serial_opt)
        serial_fqn_keys = sorted(serial_fqn_sd["state"].keys())
        assert fqn_state_keys == serial_fqn_keys, (
            f"FQN keys mismatch between serial and distributed models.\n"
            f"  Distributed: {fqn_state_keys}\n"
            f"  Serial:      {serial_fqn_keys}"
        )
    finally:
        DistributedManager.cleanup()
        DistributedManager._state = {}
        monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env",
    [
        ((2, (3, 3)), True, "cpu", "ENV"),
    ],
    indirect=("setup_env",),
    ids=["cpu-dp2-cp3x3"],
)
def test_optimizer_param_ordering_serial_vs_distributed(setup_env):
    """Goals: optimizer parameter ordering and FQN key consistency.

    - named_parameters() yields identical key lists for both model types
    - Optimizers see the same number of parameters
    - get_optimizer_state_dict produces FQN string keys matching named_parameters()
    - FQN keys are identical between serial and distributed models
    - Guards against silent optimizer state misalignment on cross-topology resume

    See also: test_optimizer_param_ordering_boltz2 for full Boltz-2 model
    coverage, which catches subtle registration-order bugs that the tiny
    harness cannot.
    """
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env
    spawn_multiprocessing(
        _parallel_assert_optimizer_param_ordering,
        world_size,
        (grid_group_sizes, device_type, backend, env_per_rank),
    )


def _parallel_assert_boltz2_optimizer_param_ordering(rank: int, payload: tuple[Any, ...]) -> None:
    """Assert optimizer parameter consistency for the full Boltz-2 model wrapper.

    The Boltz2Distributed wrapper renames placeholder module parameters with
    a ``._serial.`` prefix and may register submodules in a different order
    than the serial model (wrapped modules first, then placeholders).

    Since ``BoltzContextParallelStrategy`` uses FQN-keyed (name-based)
    optimizer state dicts, parameter *ordering* does not need to match.
    This test verifies what does matter for correctness:

    1. The canonical parameter name *sets* (stripping ``._serial.``) are
       identical between serial and distributed models — no missing or extra
       parameters.
    2. The optimizer covers exactly the set of trainable parameters.
    3. ``get_optimizer_state_dict`` produces well-formed FQN string keys
       matching ``named_parameters()``.
    4. Canonical FQN keys match between serial and distributed models,
       ensuring cross-topology checkpoint portability.
    """
    grid_group_sizes, device_type, backend, env_per_rank, serial_state_dict, serial_hparams = payload

    monkeypatch = pytest.MonkeyPatch()
    for key, value in env_per_rank.items():
        monkeypatch.setenv(key, f"{rank}" if value == "<INPUT_RANK>" else value)

    try:
        DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
        manager = DistributedManager()

        serial_model = SerialBoltz2(**serial_hparams)
        serial_model.load_state_dict(serial_state_dict, strict=True)
        serial_model = serial_model.to(manager.device)

        cp_model = Boltz2Distributed(serial_model, manager)
        cp_model = cp_model.to(manager.device)

        # Re-create a fresh serial model for comparison (the original was
        # consumed by the distributed wrapper).
        serial_model2 = SerialBoltz2(**serial_hparams)
        serial_model2.load_state_dict(serial_state_dict, strict=True)

        # Canonical name set comparison: strip ._serial. from distributed names
        cp_names = [name for name, _ in cp_model.named_parameters()]
        serial_names = [name for name, _ in serial_model2.named_parameters()]

        cp_names_canon_set = {n.replace("._serial.", ".") for n in cp_names}
        serial_names_set = set(serial_names)

        missing = serial_names_set - cp_names_canon_set
        extra = cp_names_canon_set - serial_names_set
        assert not missing and not extra, (
            f"Parameter name set mismatch between serial and distributed Boltz2.\n"
            f"  Missing from distributed: {sorted(missing)[:5]}\n"
            f"  Extra in distributed:     {sorted(extra)[:5]}"
        )

        # Build per-name param counts from the serial model and verify
        # the per-name sum equals the optimizer's total.
        cp_result = cp_model.configure_optimizers()
        serial_result = serial_model2.configure_optimizers()

        # Both return ([optimizer], [scheduler_dict]) with lr_scheduler="af3"
        cp_opt = cp_result[0][0] if isinstance(cp_result, tuple) else cp_result
        serial_opt = serial_result[0][0] if isinstance(serial_result, tuple) else serial_result

        serial_trainable_names = {n for n, p in serial_model2.named_parameters() if p.requires_grad}
        serial_opt_total = sum(len(g["params"]) for g in serial_opt.param_groups)
        assert len(serial_trainable_names) == serial_opt_total, (
            f"Serial optimizer param count ({serial_opt_total}) doesn't match "
            f"serial trainable named_parameters count ({len(serial_trainable_names)})"
        )

        # Identify which serial param names correspond to placeholder
        # modules in the distributed model.
        placeholder_prefixes = tuple(
            f"{name}." for name, mod in cp_model.named_modules() if isinstance(mod, _PlaceholderModule)
        )
        serial_placeholder_names = {n for n in serial_trainable_names if n.startswith(placeholder_prefixes)}
        serial_non_placeholder_names = serial_trainable_names - serial_placeholder_names

        # The distributed optimizer should cover exactly the non-placeholder
        # (i.e. trainable) params.
        cp_opt_count = sum(len(g["params"]) for g in cp_opt.param_groups)
        cp_trainable_names = {n for n, p in cp_model.named_parameters() if p.requires_grad}
        assert cp_opt_count == len(cp_trainable_names), (
            f"Distributed optimizer param count ({cp_opt_count}) doesn't match "
            f"distributed trainable param count ({len(cp_trainable_names)})"
        )
        assert cp_opt_count > 0, "Distributed optimizer has zero params"

        # Param count on the non-placeholder subset must match between
        # serial and distributed.
        assert cp_opt_count == len(serial_non_placeholder_names), (
            f"Distributed optimizer has {cp_opt_count} params but serial model "
            f"has {len(serial_non_placeholder_names)} non-placeholder trainable params "
            f"(total serial trainable: {len(serial_trainable_names)}, "
            f"placeholder: {len(serial_placeholder_names)})"
        )

        # Canonicalize distributed trainable names and compare against
        # the serial non-placeholder set.
        cp_trainable_canon = {n.replace("._serial.", ".") for n in cp_trainable_names}
        assert cp_trainable_canon == serial_non_placeholder_names, (
            f"Trainable distributed params don't match serial non-placeholder params.\n"
            f"  Only in distributed: {sorted(cp_trainable_canon - serial_non_placeholder_names)[:5]}\n"
            f"  Only in serial:      {sorted(serial_non_placeholder_names - cp_trainable_canon)[:5]}"
        )

        # Parameter shapes must match between serial and distributed.
        # DTensor params report their global shape, which should equal the
        # serial shape.
        serial_shapes = {n: p.shape for n, p in serial_model2.named_parameters() if p.requires_grad}
        shape_mismatches = []
        for n, p in cp_model.named_parameters():
            if not p.requires_grad:
                continue
            canon = n.replace("._serial.", ".")
            serial_shape = serial_shapes.get(canon)
            if serial_shape is None:
                continue
            cp_shape = p.shape
            if cp_shape != serial_shape:
                shape_mismatches.append((canon, cp_shape, serial_shape))
        assert not shape_mismatches, "Parameter shape mismatches (distributed vs serial):\n" + "\n".join(
            f"  {n}: {cs} vs {ss}" for n, cs, ss in shape_mismatches[:10]
        )

        # FQN keys from get_optimizer_state_dict
        for p in cp_model.parameters():
            p.grad = torch.randn_like(p.to_local() if isinstance(p, DTensor) else p)
            if isinstance(p, DTensor):
                p.grad = distribute_tensor(p.grad, device_mesh=p.device_mesh, placements=p.placements)
        cp_opt.step()

        fqn_sd = get_optimizer_state_dict(cp_model, cp_opt)
        fqn_state_keys = sorted(fqn_sd["state"].keys())

        assert all(isinstance(k, str) for k in fqn_state_keys), (
            f"get_optimizer_state_dict should return FQN string keys, "
            f"got types {[type(k).__name__ for k in fqn_state_keys[:3]]}"
        )
        # FQN state keys only cover trainable params (those in the optimizer)
        cp_trainable_names = sorted(n for n, p in cp_model.named_parameters() if p.requires_grad)
        assert fqn_state_keys == cp_trainable_names, (
            f"FQN optimizer state keys don't match trainable named_parameters().\n"
            f"  FQN keys (first 5):          {fqn_state_keys[:5]}\n"
            f"  trainable params (first 5):  {cp_trainable_names[:5]}"
        )

        # Cross-model FQN key comparison: distributed trainable params
        # (canonicalized) must be a subset of serial optimizer params.
        # Placeholder modules have their params frozen, so only the
        # "ready" distributed submodules appear in the distributed optimizer.
        for p in serial_model2.parameters():
            p.grad = torch.randn_like(p)
        serial_opt.step()

        serial_fqn_sd = get_optimizer_state_dict(serial_model2, serial_opt)
        serial_fqn_keys = set(serial_fqn_sd["state"].keys())

        cp_fqn_canon = {k.replace("._serial.", ".") for k in fqn_state_keys}
        not_in_serial = cp_fqn_canon - serial_fqn_keys
        assert not not_in_serial, (
            f"Distributed optimizer FQN keys not found in serial optimizer:\n" f"  {sorted(not_in_serial)[:5]}"
        )
    finally:
        DistributedManager.cleanup()
        DistributedManager._state = {}
        monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env",
    [
        ((2, (3, 3)), True, "cpu", "ENV"),
        ((1, (2, 2)), True, "cuda", "ENV"),
    ],
    indirect=("setup_env",),
    ids=["cpu-dp2-cp3x3", "cuda-dp1-cp2x2"],
)
def test_optimizer_param_ordering_boltz2(setup_env):
    """Boltz-2 model: optimizer parameter ordering and FQN key consistency.

    Extends the tiny-distogram harness test to the full Boltz-2 wrapper,
    catching registration-order bugs across the many submodules (ready +
    placeholder) that the smaller harness cannot exercise.
    """
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    serial_state_dict, serial_hparams = _prepare_serial_model(ema=False)

    spawn_multiprocessing(
        _parallel_assert_boltz2_optimizer_param_ordering,
        world_size,
        (grid_group_sizes, device_type, backend, env_per_rank, serial_state_dict, serial_hparams),
    )


@pytest.mark.parametrize(
    "setup_env",
    [
        ((2, (3, 3)), True, "cpu", "ENV"),
    ],
    indirect=("setup_env",),
    ids=["cpu-dp2-cp3x3"],
)
def test_cross_mode_stop_and_go(setup_env, tmp_path):
    """Goals: cross-mode checkpoint interop in both directions.

    - Serial→distributed: 1-epoch serial + resume as distributed matches
      2-epoch continuous distributed baseline
    - Distributed→serial: 1-epoch distributed + resume as serial matches
      2-epoch continuous serial baseline
    - Validates BoltzContextParallelStrategy strips DTensor metadata for
      portable checkpoints (the train-with-CP/deploy-with-serial workflow)
    - Verifies serial checkpoints use integer optimizer keys
    - Verifies distributed checkpoints use FQN string optimizer keys
    - Verifies cross-format loading: serial int-key ckpt loaded by distributed
      strategy (legacy path), and distributed FQN-key ckpt loaded by serial
      model (PyTorch positional mapping in optimizer.load_state_dict)
    """
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env
    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("CUDA not available")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"Need {world_size} GPUs, have {torch.cuda.device_count()}")
    output_dir = tmp_path / "cross_mode_output"
    payload = (grid_group_sizes, device_type, backend, env_per_rank, str(output_dir))
    spawn_multiprocessing(_parallel_assert_cross_mode_stop_and_go, world_size, payload)
