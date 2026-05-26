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

"""Error-path unit tests for BoltzContextParallelStrategy.

These test input-validation guards (missing model, wrong optimizer count,
etc.) that integration tests would not naturally trigger.

The strategy's DTensor-aware happy-path behaviour — checkpoint conversion,
optimizer-state redistribution, save/load roundtrip — is exercised
end-to-end by:

- ``test_dtensor_train.py`` (one-step optimizer parity with real DTensor model)
- ``test_dtensor_stop_and_go.py`` (full save/resume cycle through ``train()``)
"""

from dataclasses import dataclass

import pytest
import pytorch_lightning as pl
import torch

import boltz.distributed.lightning_strategy as strategy_module


@dataclass
class _DummyDistManager:
    device: torch.device = torch.device("cpu")
    rank: int = 0
    local_rank: int = 0
    world_size: int = 1


class _TinyLightningModule(pl.LightningModule):
    def __init__(self) -> None:
        super().__init__()
        self.layer = torch.nn.Linear(3, 2)

    def configure_optimizers(self):
        return torch.optim.SGD(self.parameters(), lr=0.25)


def _new_strategy() -> strategy_module.BoltzContextParallelStrategy:
    return strategy_module.BoltzContextParallelStrategy(dist_manager=_DummyDistManager())


def _attach_lightning_module(
    strategy: strategy_module.BoltzContextParallelStrategy, module: _TinyLightningModule
) -> None:
    strategy.model = module
    strategy._lightning_module = module


# ---------------------------------------------------------------------------
# lightning_module_state_dict
# ---------------------------------------------------------------------------


def test_module_state_dict_raises_without_model():
    """Goals: RuntimeError when model not attached to strategy."""
    strategy = _new_strategy()
    with pytest.raises(RuntimeError, match="model is not set"):
        strategy.lightning_module_state_dict()


# ---------------------------------------------------------------------------
# load_model_state_dict
# ---------------------------------------------------------------------------


def test_load_model_state_dict_raises_without_lightning_module():
    """Goals: RuntimeError when lightning_module not set."""
    strategy = _new_strategy()
    strategy.model = _TinyLightningModule()
    checkpoint = {"state_dict": {"layer.weight": torch.randn(2, 3), "layer.bias": torch.randn(2)}}
    with pytest.raises(RuntimeError, match="lightning_module is not set"):
        strategy.load_model_state_dict(checkpoint, strict=False)


# ---------------------------------------------------------------------------
# load_checkpoint
# ---------------------------------------------------------------------------


def test_load_checkpoint_uses_checkpoint_io_with_cpu_map_location(monkeypatch):
    """Goals: checkpoint load routes through checkpoint_io with CPU remap."""
    strategy = _new_strategy()
    expected = {"state_dict": {"layer.weight": torch.randn(2, 3)}}
    calls: dict[str, object] = {}

    class _FakeCheckpointIO:
        def load_checkpoint(self, checkpoint_path, map_location=None):
            calls["checkpoint_path"] = checkpoint_path
            calls["map_location"] = map_location
            return expected

    monkeypatch.setattr(strategy_module.torch, "load", lambda *_a, **_k: (_ for _ in ()).throw(AssertionError()))
    strategy.checkpoint_io = _FakeCheckpointIO()

    checkpoint_path = "dummy.ckpt"
    result = strategy.load_checkpoint(checkpoint_path)

    assert result is expected
    assert calls["checkpoint_path"] == checkpoint_path
    assert calls["map_location"] == "cpu"


# ---------------------------------------------------------------------------
# load_optimizer_state_dict
# ---------------------------------------------------------------------------


def test_load_optimizer_state_dict_raises_without_optimizer_states():
    """Goals: ValueError when checkpoint has no optimizer_states key."""
    strategy = _new_strategy()
    module = _TinyLightningModule()
    _attach_lightning_module(strategy, module)
    strategy.optimizers = [torch.optim.SGD(module.parameters(), lr=0.25)]
    with pytest.raises(ValueError, match="no optimizer_states found"):
        strategy.load_optimizer_state_dict({})


def test_load_optimizer_state_dict_raises_on_length_mismatch():
    """Goals: ValueError when checkpoint has wrong number of optimizer states."""
    strategy = _new_strategy()
    module = _TinyLightningModule()
    _attach_lightning_module(strategy, module)
    strategy.optimizers = [
        torch.optim.SGD(module.parameters(), lr=0.1),
        torch.optim.Adam(module.parameters(), lr=0.1),
    ]
    with pytest.raises(ValueError, match="length mismatch"):
        strategy.load_optimizer_state_dict({"optimizer_states": [strategy.optimizers[0].state_dict()]})


def test_load_optimizer_state_dict_raises_on_non_sequence():
    """Goals: TypeError when optimizer_states is not a list/tuple."""
    strategy = _new_strategy()
    module = _TinyLightningModule()
    _attach_lightning_module(strategy, module)
    strategy.optimizers = [torch.optim.SGD(module.parameters(), lr=0.1)]
    with pytest.raises(TypeError, match="must be a list/tuple"):
        strategy.load_optimizer_state_dict({"optimizer_states": {"state": {}}})


def test_load_optimizer_state_dict_no_optimizers():
    """Goals: silently return when no optimizers are attached."""
    strategy = _new_strategy()
    strategy.optimizers = []
    # Should not raise.
    strategy.load_optimizer_state_dict({"optimizer_states": [{"state": {}}]})


# ---------------------------------------------------------------------------
# barrier
# ---------------------------------------------------------------------------


def test_barrier_no_distributed(monkeypatch):
    """Goals: barrier() is a no-op when torch.distributed is not initialized."""
    strategy = _new_strategy()
    monkeypatch.setattr(torch.distributed, "is_initialized", lambda: False)
    # Should not raise.
    strategy.barrier()
