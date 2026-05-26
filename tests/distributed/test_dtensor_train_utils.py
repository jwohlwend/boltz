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

"""Unit tests for distributed train.py utility functions.

These tests cover the pure-logic helpers that don't require a running
distributed process group, catching configuration errors before a full
training launch.
"""

import os
from unittest.mock import MagicMock

import pytest
import torch

import boltz.distributed.train as train_module
from boltz.distributed.model.modules.utils import Precision, setup_tf32_env
from boltz.distributed.train import (
    DistributedTrainConfig,
    _apply_matmul_precision,
    _create_dist_manager,
    _parse_precision,
)


def _make_config(**overrides: object) -> DistributedTrainConfig:
    return DistributedTrainConfig(
        data=MagicMock(),
        model=MagicMock(),
        output="/tmp/test",
        **overrides,
    )


# ---------------------------------------------------------------------------
# _parse_precision
# ---------------------------------------------------------------------------
def test_parse_precision_returns_precision_if_already_precision() -> None:
    """Goals: Precision enum passthrough."""
    assert _parse_precision(Precision.BF16) is Precision.BF16


@pytest.mark.parametrize(
    "name,expected",
    [
        ("FP32", Precision.FP32),
        ("BF16", Precision.BF16),
        ("BF16_MIXED", Precision.BF16_MIXED),
        ("FP16", Precision.FP16),
        ("TF32", Precision.TF32),
        ("FP64", Precision.FP64),
    ],
)
def test_parse_precision_from_string(name: str, expected: Precision) -> None:
    """Goals: string name → Precision enum conversion."""
    assert _parse_precision(name) is expected


@pytest.mark.parametrize("bad_value", ["INVALID", 42, None])
def test_parse_precision_raises_for_unsupported_value(bad_value: object) -> None:
    """Goals: ValueError for invalid precision inputs."""
    with pytest.raises(ValueError, match="Unsupported precision value"):
        _parse_precision(bad_value)


# ---------------------------------------------------------------------------
# _apply_matmul_precision
# ---------------------------------------------------------------------------
def test_apply_matmul_precision_noop_when_unset(monkeypatch: pytest.MonkeyPatch) -> None:
    """Goals: None value does not call set_float32_matmul_precision."""
    calls: list[str] = []

    def _record(value: str) -> None:
        calls.append(value)

    monkeypatch.setattr(torch, "set_float32_matmul_precision", _record)
    _apply_matmul_precision(None)
    assert calls == []


def test_apply_matmul_precision_applies_when_set(monkeypatch: pytest.MonkeyPatch) -> None:
    """Goals: string value is forwarded to set_float32_matmul_precision."""
    calls: list[str] = []

    def _record(value: str) -> None:
        calls.append(value)

    monkeypatch.setattr(torch, "set_float32_matmul_precision", _record)
    _apply_matmul_precision("high")
    assert calls == ["high"]


# ---------------------------------------------------------------------------
# DistributedTrainConfig
# ---------------------------------------------------------------------------
def test_distributed_train_config_defaults() -> None:
    """Goals: DistributedTrainConfig defaults match expected Boltz-2 production values."""
    cfg = _make_config()
    assert cfg.precision is Precision.FP32
    assert cfg.seed is None
    assert cfg.matmul_precision is None
    assert cfg.find_unused_parameters is False
    assert cfg.save_top_k == 1
    assert cfg.validation_only is False


# ---------------------------------------------------------------------------
# _create_dist_manager (already initialized mismatch guards)
# ---------------------------------------------------------------------------
def test_create_dist_manager_raises_on_initialized_device_type_mismatch(monkeypatch: pytest.MonkeyPatch) -> None:
    """Goals: explicit error when reusing singleton with different accelerator."""

    class _FakeDistManager:
        @staticmethod
        def is_initialized() -> bool:
            return True

        def __init__(self) -> None:
            self.device = torch.device("cpu")
            self.group_ranks = {"dp": [0], "cp": [0]}

    monkeypatch.setattr(train_module, "DistributedManager", _FakeDistManager)
    cfg = _make_config(trainer={"accelerator": "gpu"}, parallel_size={"size_dp": 1, "size_cp": 1})
    with pytest.raises(ValueError, match="Cannot change device type"):
        _create_dist_manager(cfg)


def test_create_dist_manager_raises_on_initialized_topology_mismatch(monkeypatch: pytest.MonkeyPatch) -> None:
    """Goals: explicit error when reusing singleton with different dp/cp sizes."""

    class _FakeDistManager:
        @staticmethod
        def is_initialized() -> bool:
            return True

        def __init__(self) -> None:
            self.device = torch.device("cpu")
            self.group_ranks = {"dp": [0], "cp": [0]}

    monkeypatch.setattr(train_module, "DistributedManager", _FakeDistManager)
    cfg = _make_config(trainer={"accelerator": "cpu"}, parallel_size={"size_dp": 2, "size_cp": 1})
    with pytest.raises(ValueError, match="Cannot change topology"):
        _create_dist_manager(cfg)


# ---------------------------------------------------------------------------
# setup_tf32_env
# ---------------------------------------------------------------------------


@pytest.fixture(autouse=False)
def _tf32_snapshot():
    """Save and restore TF32 global state around each test that uses it."""
    orig_matmul = torch.backends.cuda.matmul.allow_tf32
    orig_cudnn = torch.backends.cudnn.allow_tf32
    orig_env = os.environ.get("NVIDIA_TF32_OVERRIDE")
    yield orig_matmul, orig_cudnn, orig_env
    torch.backends.cuda.matmul.allow_tf32 = orig_matmul
    torch.backends.cudnn.allow_tf32 = orig_cudnn
    if orig_env is not None:
        os.environ["NVIDIA_TF32_OVERRIDE"] = orig_env
    else:
        os.environ.pop("NVIDIA_TF32_OVERRIDE", None)


def test_setup_tf32_env_tf32_precision_enables_tf32(_tf32_snapshot) -> None:
    """Goals: TF32 precision sets NVIDIA_TF32_OVERRIDE=1 and enables matmul/cudnn TF32."""
    orig_matmul, orig_cudnn, orig_env = _tf32_snapshot
    with setup_tf32_env(Precision.TF32):
        assert os.environ.get("NVIDIA_TF32_OVERRIDE") == "1"
        assert torch.backends.cuda.matmul.allow_tf32 is True
        assert torch.backends.cudnn.allow_tf32 is True
    # Context manager should restore originals
    assert torch.backends.cuda.matmul.allow_tf32 == orig_matmul
    assert torch.backends.cudnn.allow_tf32 == orig_cudnn
    if orig_env is not None:
        assert os.environ.get("NVIDIA_TF32_OVERRIDE") == orig_env
    else:
        assert "NVIDIA_TF32_OVERRIDE" not in os.environ


def test_setup_tf32_env_fp32_precision_disables_tf32(_tf32_snapshot) -> None:
    """Goals: FP32 precision explicitly disables TF32."""
    with setup_tf32_env(Precision.FP32):
        assert os.environ.get("NVIDIA_TF32_OVERRIDE") == "0"
        assert torch.backends.cuda.matmul.allow_tf32 is False
        assert torch.backends.cudnn.allow_tf32 is False


def test_setup_tf32_env_bf16_precision_leaves_tf32_unchanged(_tf32_snapshot) -> None:
    """Goals: BF16 precision does not modify TF32 state."""
    orig_matmul, orig_cudnn, _ = _tf32_snapshot
    with setup_tf32_env(Precision.BF16):
        assert torch.backends.cuda.matmul.allow_tf32 == orig_matmul
        assert torch.backends.cudnn.allow_tf32 == orig_cudnn


def test_setup_tf32_env_restores_on_exception(_tf32_snapshot) -> None:
    """Goals: TF32 state is restored even when the context manager body raises."""
    orig_matmul, orig_cudnn, _ = _tf32_snapshot
    with pytest.raises(RuntimeError, match="boom"):
        with setup_tf32_env(Precision.TF32):
            raise RuntimeError("boom")
    assert torch.backends.cuda.matmul.allow_tf32 == orig_matmul
    assert torch.backends.cudnn.allow_tf32 == orig_cudnn
