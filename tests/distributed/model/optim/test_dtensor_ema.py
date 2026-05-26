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

"""Tests for DistributedEMA – the DTensor-aware EMA callback.

All tests exercise DistributedEMA with DTensor parameters to validate
numerical parity with the base EMA, save/load serialisation, and the
eval-swap lifecycle.
"""

from __future__ import annotations

from typing import Any
from unittest.mock import MagicMock

import pytest
import pytorch_lightning as pl
import torch
from torch import nn
from torch.distributed.tensor import DTensor, distribute_module

from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.optim.ema import DistributedEMA
from boltz.model.optim.ema import EMA
from boltz.testing.utils import assert_all_identical, assert_tensors_identical, spawn_multiprocessing


def _distribute_model(model: pl.LightningModule, manager: DistributedManager) -> None:
    """Distribute all sub-modules' parameters as replicated DTensors on the CP mesh."""
    for _name, child in model.named_children():
        distribute_module(child, manager.device_mesh_subgroups)


class _TinyModule(pl.LightningModule):
    """Minimal model for EMA tests."""

    def __init__(self) -> None:
        super().__init__()
        self.layer = nn.Linear(4, 2, bias=True)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.layer(x)

    def training_step(self, batch: torch.Tensor, batch_idx: int) -> torch.Tensor:
        return self(batch).sum()

    def configure_optimizers(self):
        return torch.optim.SGD(self.parameters(), lr=0.01)


# ---------------------------------------------------------------------------
# Distributed tests (require process group)
# ---------------------------------------------------------------------------


def _parallel_assert_ema_comprehensive(rank: int, payload: tuple[Any, ...]) -> None:
    """Comprehensive EMA lifecycle with DTensors.

    Tests in sequence: init → multi-step parity with serial EMA →
    cross-rank identity → save (DTensor→plain) → load (plain→DTensor
    realignment) → eval swap lifecycle (replace/forward/restore).
    """
    grid_group_sizes, device_type, backend, env_per_rank = payload

    monkeypatch = pytest.MonkeyPatch()
    for key, value in env_per_rank.items():
        monkeypatch.setenv(key, f"{rank}" if value == "<INPUT_RANK>" else value)

    try:
        DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
        manager = DistributedManager()

        # Create distributed + serial models with identical init.
        with torch.random.fork_rng(devices=[]):
            torch.manual_seed(42)
            model = _TinyModule()
        model.to(manager.device)
        _distribute_model(model, manager)

        with torch.random.fork_rng(devices=[]):
            torch.manual_seed(42)
            serial_model = _TinyModule()
        serial_model.to(manager.device)

        for param in model.state_dict().values():
            assert isinstance(param, DTensor)

        # Init EMA on both.
        dist_ema = DistributedEMA(decay=0.9, eval_with_ema=True, warm_start=True)
        serial_ema = EMA(decay=0.9, warm_start=True)
        mock_trainer = MagicMock()
        dist_ema.on_train_start(mock_trainer, model)
        serial_ema.on_train_start(mock_trainer, serial_model)

        for w in dist_ema._ema_weights.values():
            assert isinstance(w, DTensor)

        # Run 5 EMA steps with matching perturbations.
        for step in range(5):
            with torch.random.fork_rng(devices=[]):
                torch.manual_seed(100 + step)
                deltas = [torch.randn_like(p) * 0.3 for p in serial_model.parameters()]

            with torch.no_grad():
                for p, d in zip(serial_model.parameters(), deltas):
                    p.add_(d)
                serial_sd = serial_model.state_dict()
                for n, p in model.named_parameters():
                    (p.data.to_local() if isinstance(p.data, DTensor) else p).copy_(serial_sd[n])

            serial_ema._cur_step = dist_ema._cur_step = step
            mock_trainer.global_step = step + 1
            serial_ema.on_train_batch_end(mock_trainer, serial_model, None, None, 0)
            dist_ema.on_train_batch_end(mock_trainer, model, None, None, 0)

        # Bitwise numerical parity with serial.
        # Use full_tensor() for DTensor-vs-serial comparison (convention).
        for key in serial_ema._ema_weights:
            dw = dist_ema._ema_weights[key]
            assert isinstance(dw, DTensor)
            assert_tensors_identical(dw.full_tensor(), serial_ema._ema_weights[key])

        # EMA actually had an effect: EMA weights must differ from the
        # current model weights (EMA lags behind the perturbed model).
        for key in dist_ema._ema_weights:
            model_w = model.state_dict()[key].to_local()
            ema_w = dist_ema._ema_weights[key].to_local()
            assert not torch.equal(model_w, ema_w), (
                f"EMA weight '{key}' is identical to model weight after 5 steps — " "EMA may not be accumulating"
            )

        # Cross-rank identity.
        world_group = torch.distributed.distributed_c10d._get_default_group()
        for ema_w in dist_ema._ema_weights.values():
            assert_all_identical(ema_w.to_local(), world_group)

        # Save: DTensor EMA weights → plain tensors.
        checkpoint: dict[str, Any] = {}
        dist_ema.on_save_checkpoint(mock_trainer, model, checkpoint)
        assert "ema" in checkpoint
        for w in checkpoint["ema"]["ema_weights"].values():
            assert not isinstance(w, DTensor)

        # Load into fresh EMA: plain → DTensor realignment.
        fresh_ema = DistributedEMA(decay=0.9, eval_with_ema=True)
        fresh_ema.on_load_checkpoint(mock_trainer, model, checkpoint)
        fresh_ema.on_train_start(mock_trainer, model)
        assert fresh_ema._cur_step == dist_ema._cur_step
        for key, loaded_w in fresh_ema._ema_weights.items():
            assert isinstance(loaded_w, DTensor)
            assert loaded_w.placements == model.state_dict()[key].placements
            torch.testing.assert_close(loaded_w.to_local(), dist_ema._ema_weights[key].to_local())

        # Eval swap lifecycle.
        training_weights = {k: v.to_local().clone() for k, v in model.state_dict().items()}
        ema_snapshot = {k: v.to_local().clone() for k, v in fresh_ema._ema_weights.items()}

        fresh_ema.on_validation_start(mock_trainer, model)

        # Backup must be on CPU to avoid doubling GPU memory during validation.
        for k, buf in fresh_ema._weights_buffer.items():
            assert buf.device.type == "cpu", (
                f"Weights buffer '{k}' should be on CPU to save GPU memory, " f"but is on {buf.device}"
            )
            assert not isinstance(buf, DTensor), f"Weights buffer '{k}' should be a plain tensor on CPU, not a DTensor"

        for k, v in model.state_dict().items():
            assert isinstance(v, DTensor), f"DTensor semantics lost during replace for '{k}'"
            assert_tensors_identical(v.to_local(), ema_snapshot[k])

        fresh_ema.on_validation_end(mock_trainer, model)
        for k, v in model.state_dict().items():
            assert isinstance(v, DTensor), f"DTensor semantics lost during restore for '{k}'"
            assert_tensors_identical(v.to_local(), training_weights[k])
    finally:
        DistributedManager.cleanup()
        DistributedManager._state = {}
        monkeypatch.undo()


def _parallel_assert_ema_inference_mode(rank: int, payload: tuple[Any, ...]) -> None:
    """EMA replace/restore under torch.inference_mode(True).

    PyTorch >=2.10 disallows version-counter manipulation on inference
    tensors.  Lightning's trainer.predict() wraps the workflow in
    inference_mode(True), so replace_model_weights and
    restore_original_weights must handle this.
    """
    grid_group_sizes, device_type, backend, env_per_rank = payload

    monkeypatch = pytest.MonkeyPatch()
    for key, value in env_per_rank.items():
        monkeypatch.setenv(key, f"{rank}" if value == "<INPUT_RANK>" else value)

    try:
        DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
        manager = DistributedManager()

        with torch.random.fork_rng(devices=[]):
            torch.manual_seed(42)
            model = _TinyModule()
        model.to(manager.device)
        _distribute_model(model, manager)

        ema = DistributedEMA(decay=0.9, eval_with_ema=True)
        mock_trainer = MagicMock()
        ema.on_train_start(mock_trainer, model)

        # Perturb EMA weights so they differ from model weights
        with torch.no_grad():
            for w in ema._ema_weights.values():
                w.add_(0.5)

        training_weights = {k: v.to_local().clone() for k, v in model.state_dict().items()}
        ema_snapshot = {k: v.to_local().clone() for k, v in ema._ema_weights.items()}

        with torch.inference_mode(True):
            ema.replace_model_weights(model)
            for k, v in model.state_dict().items():
                assert isinstance(v, DTensor), f"DTensor lost for '{k}' under inference_mode replace"
                assert_tensors_identical(v.to_local(), ema_snapshot[k])

            ema.restore_original_weights(model)
            for k, v in model.state_dict().items():
                assert isinstance(v, DTensor), f"DTensor lost for '{k}' under inference_mode restore"
                assert_tensors_identical(v.to_local(), training_weights[k])
    finally:
        DistributedManager.cleanup()
        DistributedManager._state = {}
        monkeypatch.undo()


def _parallel_assert_ema_shape_mismatch_raises(rank: int, payload: tuple[Any, ...]) -> None:
    """_realign_weights_to_model raises ValueError on shape mismatch."""
    grid_group_sizes, device_type, backend, env_per_rank = payload

    monkeypatch = pytest.MonkeyPatch()
    for key, value in env_per_rank.items():
        monkeypatch.setenv(key, f"{rank}" if value == "<INPUT_RANK>" else value)

    try:
        DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
        manager = DistributedManager()

        with torch.random.fork_rng(devices=[]):
            torch.manual_seed(42)
            model = _TinyModule()
        model.to(manager.device)
        _distribute_model(model, manager)

        try:
            DistributedEMA._realign_weights_to_model(
                {"layer.weight": torch.randn(99, 99), "layer.bias": torch.randn(2)}, model
            )
            raise AssertionError("Expected ValueError for shape mismatch")  # noqa: TRY301
        except ValueError as e:
            assert "shape mismatch" in str(e).lower(), f"Unexpected error message: {e}"
    finally:
        DistributedManager.cleanup()
        DistributedManager._state = {}
        monkeypatch.undo()


# ---------------------------------------------------------------------------
# Distributed test parametrisation
# ---------------------------------------------------------------------------


_EMA_TOPOLOGIES = [
    ((1, (2, 2)), True, "cpu", "ENV"),
    ((2, (1, 1)), True, "cuda", "ENV"),
]
_EMA_IDS = ["cpu-dp1-cp2x2", "cuda-dp2-cp1x1"]


@pytest.mark.parametrize("setup_env", _EMA_TOPOLOGIES, indirect=("setup_env",), ids=_EMA_IDS)
def test_ema_comprehensive_lifecycle(setup_env):
    """Goals: full EMA lifecycle — parity, save/load roundtrip, eval swap, cross-rank identity."""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env
    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("CUDA not available")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"Need {world_size} GPUs, have {torch.cuda.device_count()}")
    spawn_multiprocessing(
        _parallel_assert_ema_comprehensive, world_size, (grid_group_sizes, device_type, backend, env_per_rank)
    )


@pytest.mark.parametrize("setup_env", _EMA_TOPOLOGIES, indirect=("setup_env",), ids=_EMA_IDS)
def test_ema_inference_mode(setup_env):
    """Goals: replace/restore under inference_mode(True) — no version_counter error."""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env
    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("CUDA not available")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"Need {world_size} GPUs, have {torch.cuda.device_count()}")
    spawn_multiprocessing(
        _parallel_assert_ema_inference_mode, world_size, (grid_group_sizes, device_type, backend, env_per_rank)
    )


@pytest.mark.parametrize("setup_env", _EMA_TOPOLOGIES, indirect=("setup_env",), ids=_EMA_IDS)
def test_ema_shape_mismatch_raises(setup_env):
    """Goals: _realign_weights_to_model raises ValueError when EMA shape != model shape."""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env
    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("CUDA not available")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"Need {world_size} GPUs, have {torch.cuda.device_count()}")
    spawn_multiprocessing(
        _parallel_assert_ema_shape_mismatch_raises, world_size, (grid_group_sizes, device_type, backend, env_per_rank)
    )
