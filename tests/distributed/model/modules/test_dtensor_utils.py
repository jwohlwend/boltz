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

"""Unit tests for DTensor checkpoint conversion helpers."""

import socket
from collections import OrderedDict

import pytest
import torch
from torch.distributed.tensor import DTensor, Replicate, Shard, distribute_tensor

from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.modules.utils import (
    OffloadActvCkptToCPU,
    SDPAWithBiasBackend,
    SetAttnPairBiasBackend,
    SetAttnPairBiasShardwiseBackend,
    SetTriAttnBackend,
    TriAttnBackend,
    _convert_serial_value_to_template_layout,
    convert_distributed_checkpoint_to_serial_state_dict,
    convert_dtensors_to_tensors,
    convert_serial_checkpoint_to_distributed_state_dict,
    has_dtensors,
)
from boltz.testing.utils import create_boltz2_model_init_params, spawn_multiprocessing


def _find_free_port() -> int:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
        sock.bind(("127.0.0.1", 0))
        return int(sock.getsockname()[1])


@pytest.fixture
def single_rank_dist_manager(monkeypatch):
    # DistributedManager cleanup currently destroys process groups but does not
    # always clear singleton state; force reset so each test starts clean.
    DistributedManager._state = {}

    port = str(_find_free_port())
    monkeypatch.setenv("MASTER_ADDR", "127.0.0.1")
    monkeypatch.setenv("MASTER_PORT", port)
    monkeypatch.setenv("WORLD_SIZE", "1")
    monkeypatch.setenv("RANK", "0")
    monkeypatch.setenv("LOCAL_RANK", "0")

    grid_group_sizes = OrderedDict(dp=1, cp=(1, 1))
    backend = DistributedManager.backend_for_device()["cpu"]
    assert backend is not None, "Gloo backend must be available for CPU DTensor tests"
    DistributedManager.initialize(grid_group_sizes, device_type="cpu", backend=backend)
    manager = DistributedManager()
    yield manager
    DistributedManager.cleanup()
    DistributedManager._state = {}


def _as_replicated_dtensor(tensor: torch.Tensor, manager: DistributedManager) -> DTensor:
    return distribute_tensor(
        tensor.to(manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=(Replicate(), Replicate(), Replicate()),
    )


def test_has_dtensors_detects_nested_dtensors(single_rank_dist_manager: DistributedManager):
    """Goals: has_dtensors finds DTensors inside nested dicts/lists."""
    dtensor = _as_replicated_dtensor(torch.randn(2, 3), single_rank_dist_manager)
    nested = {"a": [torch.ones(1), {"b": dtensor}], "c": (1, 2, 3)}
    assert has_dtensors(nested)
    assert not has_dtensors({"a": [1, 2], "b": {"c": torch.tensor([1.0])}})


def test_convert_dtensors_to_tensors_recursively(single_rank_dist_manager: DistributedManager):
    """Goals: convert_dtensors_to_tensors strips DTensor metadata through nested containers."""
    dtensor = _as_replicated_dtensor(torch.randn(2, 3), single_rank_dist_manager)
    nested = {"x": dtensor, "y": [dtensor], "z": ({"k": dtensor},)}

    converted = convert_dtensors_to_tensors(nested)

    assert isinstance(converted["x"], torch.Tensor)
    assert not isinstance(converted["x"], DTensor)
    assert not has_dtensors(converted)
    torch.testing.assert_close(converted["x"], dtensor.to_local())
    torch.testing.assert_close(converted["y"][0], dtensor.to_local())
    torch.testing.assert_close(converted["z"][0]["k"], dtensor.to_local())


def test_checkpoint_roundtrip_serial_to_distributed_to_serial(single_rank_dist_manager: DistributedManager):
    """Goals: serial→distributed→serial roundtrip preserves tensor values."""
    serial_weight = torch.randn(4, 4)
    serial_bias = torch.randn(4)

    state_template = {
        "layer.weight": _as_replicated_dtensor(torch.zeros_like(serial_weight), single_rank_dist_manager),
        "layer.bias": torch.zeros_like(serial_bias),
    }
    checkpoint = {
        "state_dict": {
            "layer.weight": serial_weight.clone(),
            "layer.bias": serial_bias.clone(),
        }
    }

    distributed_state = convert_serial_checkpoint_to_distributed_state_dict(
        checkpoint=checkpoint,
        strict=True,
        state_dict_template=state_template,
    )

    assert isinstance(distributed_state["layer.weight"], DTensor)
    assert not isinstance(distributed_state["layer.bias"], DTensor)
    assert has_dtensors(distributed_state)

    roundtrip_state = convert_distributed_checkpoint_to_serial_state_dict({"state_dict": distributed_state})
    assert not has_dtensors(roundtrip_state)
    torch.testing.assert_close(roundtrip_state["layer.weight"], serial_weight)
    torch.testing.assert_close(roundtrip_state["layer.bias"], serial_bias)


def test_checkpoint_conversion_strict_key_mismatch_raises(single_rank_dist_manager: DistributedManager):
    """Goals: strict mode raises KeyError when checkpoint and template keys differ."""
    checkpoint = {"state_dict": {"foo": torch.tensor([1.0])}}
    template = {"bar": torch.tensor([1.0])}

    with pytest.raises(KeyError, match="State-dict keys do not match template keys"):
        convert_serial_checkpoint_to_distributed_state_dict(
            checkpoint=checkpoint,
            strict=True,
            state_dict_template=template,
        )


# ---------------------------------------------------------------------------
# Error-path tests (parametrized)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "checkpoint,match",
    [
        ({"meta": 1}, "state_dict"),
        ({"state_dict": [1, 2, 3]}, "must be a mapping"),
    ],
    ids=["missing-state-dict", "non-mapping-state-dict"],
)
def test_to_serial_error_paths(checkpoint, match):
    """Goals: convert_distributed_checkpoint_to_serial_state_dict rejects invalid input."""
    with pytest.raises((KeyError, TypeError), match=match):
        convert_distributed_checkpoint_to_serial_state_dict(checkpoint)


@pytest.mark.parametrize(
    "checkpoint,template,error_type,match",
    [
        ({"meta": 1}, {"k": torch.tensor([1.0])}, KeyError, "state_dict"),
        ({"state_dict": {"k": torch.tensor([1.0])}}, None, ValueError, "state_dict_template is required"),
        ({"state_dict": "not_a_dict"}, {"k": torch.tensor([1.0])}, TypeError, "must be a mapping"),
    ],
    ids=["missing-state-dict", "no-template", "non-mapping-state-dict"],
)
def test_to_distributed_error_paths(checkpoint, template, error_type, match):
    """Goals: convert_serial_checkpoint_to_distributed_state_dict rejects invalid input."""
    with pytest.raises(error_type, match=match):
        convert_serial_checkpoint_to_distributed_state_dict(
            checkpoint=checkpoint,
            state_dict_template=template,
        )


# ---------------------------------------------------------------------------
# Non-strict mode: extra keys in checkpoint are preserved
# ---------------------------------------------------------------------------


def test_to_distributed_non_strict_preserves_extra_keys():
    """Extra checkpoint keys not in the template are passed through in non-strict mode."""
    extra_tensor = torch.tensor([42.0])
    checkpoint = {
        "state_dict": {
            "in_template": torch.tensor([1.0]),
            "extra_key": extra_tensor.clone(),
        }
    }
    template = {"in_template": torch.tensor([0.0])}

    result = convert_serial_checkpoint_to_distributed_state_dict(
        checkpoint=checkpoint,
        state_dict_template=template,
        strict=False,
    )

    assert "in_template" in result
    assert "extra_key" in result
    torch.testing.assert_close(result["extra_key"], extra_tensor)


def test_to_distributed_non_strict_ignores_missing_template_keys():
    """Template keys absent from checkpoint are silently skipped in non-strict mode."""
    checkpoint = {"state_dict": {"a": torch.tensor([1.0])}}
    template = {"a": torch.tensor([0.0]), "b": torch.tensor([0.0])}

    result = convert_serial_checkpoint_to_distributed_state_dict(
        checkpoint=checkpoint,
        state_dict_template=template,
        strict=False,
    )

    assert "a" in result
    assert "b" not in result


# ---------------------------------------------------------------------------
# _convert_serial_value_to_template_layout: shape mismatch
# ---------------------------------------------------------------------------


def test_shape_mismatch_raises_value_error(single_rank_dist_manager: DistributedManager):
    """ValueError when serial tensor shape does not match DTensor template shape."""
    template_dtensor = _as_replicated_dtensor(torch.zeros(4, 4), single_rank_dist_manager)
    wrong_shape_tensor = torch.randn(3, 5)

    with pytest.raises(ValueError, match="does not match template shape"):
        _convert_serial_value_to_template_layout(wrong_shape_tensor, template_dtensor)


def test_dtensor_to_dtensor_passthrough(single_rank_dist_manager: DistributedManager):
    """Goals: DTensor value with DTensor template is returned as-is (no re-distribution)."""
    dtensor = _as_replicated_dtensor(torch.randn(3, 3), single_rank_dist_manager)
    template_dtensor = _as_replicated_dtensor(torch.zeros(3, 3), single_rank_dist_manager)

    result = _convert_serial_value_to_template_layout(dtensor, template_dtensor)

    assert isinstance(result, DTensor)
    assert result is dtensor  # should be the exact same object


def _parallel_assert_sharded_template_checkpoint_conversion(rank: int, payload):
    grid_group_sizes, device_type, backend, env_per_rank, serial_weight, serial_bias = payload

    monkeypatch = pytest.MonkeyPatch()
    for var_name, value in env_per_rank.items():
        monkeypatch.setenv(var_name, f"{rank}" if value == "<INPUT_RANK>" else value)

    DistributedManager._state = {}
    try:
        DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
        manager = DistributedManager()

        template_weight = distribute_tensor(
            torch.zeros_like(serial_weight, device=manager.device),
            device_mesh=manager.device_mesh_subgroups,
            placements=(Replicate(), Shard(0), Shard(1)),
        )
        template_bias = torch.zeros_like(serial_bias, device=manager.device)
        template = {
            "layer.weight": template_weight,
            "layer.bias": template_bias,
        }
        checkpoint = {
            "state_dict": {
                "layer.weight": serial_weight.clone(),
                "layer.bias": serial_bias.clone(),
            }
        }

        converted = convert_serial_checkpoint_to_distributed_state_dict(
            checkpoint=checkpoint,
            strict=True,
            state_dict_template=template,
        )

        weight_dtensor = converted["layer.weight"]
        assert isinstance(weight_dtensor, DTensor)
        assert weight_dtensor.placements == template_weight.placements
        torch.testing.assert_close(weight_dtensor.full_tensor().cpu(), serial_weight)

        cp_layout = manager.layout_subgroups["cp"]
        cp_rank = manager.group_rank["cp"]
        i, j = cp_layout.unravel(cp_rank)
        expected_local = torch.chunk(serial_weight, cp_layout.shape[0], dim=0)[i]
        expected_local = torch.chunk(expected_local, cp_layout.shape[1], dim=1)[j]
        torch.testing.assert_close(weight_dtensor.to_local().cpu(), expected_local)
        torch.testing.assert_close(converted["layer.bias"].cpu(), serial_bias)

        serialized = convert_distributed_checkpoint_to_serial_state_dict({"state_dict": converted})
        assert not has_dtensors(serialized)
        # Sharded DTensors must serialize as full global tensors for topology portability.
        torch.testing.assert_close(serialized["layer.weight"], serial_weight)
        torch.testing.assert_close(serialized["layer.bias"], serial_bias)

        # Regression guard: distributed -> serial -> distributed roundtrip for sharded layout.
        roundtrip = convert_serial_checkpoint_to_distributed_state_dict(
            checkpoint={"state_dict": serialized},
            strict=True,
            state_dict_template=template,
        )
        roundtrip_weight = roundtrip["layer.weight"]
        assert isinstance(roundtrip_weight, DTensor)
        assert roundtrip_weight.placements == template_weight.placements
        torch.testing.assert_close(roundtrip_weight.full_tensor().cpu(), serial_weight)
        torch.testing.assert_close(roundtrip_weight.to_local().cpu(), expected_local)
    finally:
        DistributedManager.cleanup()
        DistributedManager._state = {}
        monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env",
    [
        ((1, (2, 2)), True, "cpu", "ENV"),
    ],
    indirect=("setup_env",),
    ids=["cpu-dp1-cp2x2"],
)
def test_sharded_template_checkpoint_conversion_multi_rank(setup_env):
    """Goals: serial→distributed conversion with Shard placements distributes data correctly across ranks."""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env
    serial_weight = torch.arange(16, dtype=torch.float32).reshape(4, 4)
    serial_bias = torch.arange(4, dtype=torch.float32)
    payload = (grid_group_sizes, device_type, backend, env_per_rank, serial_weight, serial_bias)
    spawn_multiprocessing(_parallel_assert_sharded_template_checkpoint_conversion, world_size, payload)


# ---------------------------------------------------------------------------
# SetTriAttnBackend
# ---------------------------------------------------------------------------


def _parallel_assert_set_triattn_backend(rank, env_per_rank, triattn_backend, boltz2_params):
    """Worker: verify SetTriAttnBackend targets only PairformerLayer instances."""
    monkeypatch = pytest.MonkeyPatch()
    if env_per_rank is not None:
        for var_name, value in env_per_rank.items():
            monkeypatch.setenv(var_name, f"{rank}" if value == "<INPUT_RANK>" else value)

    from boltz.distributed.model.layers.pairformer import PairformerLayer
    from boltz.distributed.model.models.boltz2 import Boltz2 as Boltz2Distributed
    from boltz.model.models.boltz2 import Boltz2 as SerialBoltz2

    grid_group_sizes = {"dp": 1, "cp": (1, 1)}
    DistributedManager.initialize(grid_group_sizes, device_type="cuda", backend="nccl")
    manager = DistributedManager()

    serial_model = SerialBoltz2(**boltz2_params).to(device=manager.device).eval()
    dist_model = Boltz2Distributed(serial_model, manager).eval()

    pairformer_layers_before = []
    for name, submodule in dist_model.named_modules():
        if isinstance(submodule, PairformerLayer):
            pairformer_layers_before.append(name)
            assert (
                submodule.triattn_backend == TriAttnBackend.REFERENCE
            ), f"{name}: expected REFERENCE before setter, got {submodule.triattn_backend}"

    assert len(pairformer_layers_before) > 0, "Model must contain at least one PairformerLayer"

    dist_model.apply(SetTriAttnBackend(triattn_backend))

    for name, submodule in dist_model.named_modules():
        if isinstance(submodule, PairformerLayer):
            assert (
                submodule.triattn_backend == triattn_backend
            ), f"{name}: expected {triattn_backend} after setter, got {submodule.triattn_backend}"
        else:
            assert not hasattr(submodule, "triattn_backend") or isinstance(submodule, PairformerLayer), (
                f"Non-PairformerLayer module {name} ({type(submodule).__name__}) "
                f"unexpectedly has triattn_backend attribute"
            )

    DistributedManager.cleanup()


@pytest.mark.parametrize(
    "setup_env",
    [((1, (1, 1)), False, "cuda", "ENV")],
    indirect=True,
    ids=["cuda-dp1-cp1x1"],
)
@pytest.mark.parametrize(
    "triattn_backend",
    [TriAttnBackend.CUEQ, TriAttnBackend.TRIFAST],
    ids=lambda b: b.value,
)
def test_set_triattn_backend(setup_env, triattn_backend):
    """SetTriAttnBackend sets triattn_backend only on PairformerLayer instances."""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if not torch.cuda.is_available():
        pytest.skip("CUDA not available")
    if torch.cuda.device_count() < world_size:
        pytest.skip(f"Need {world_size} GPUs, have {torch.cuda.device_count()}")

    spawn_multiprocessing(
        _parallel_assert_set_triattn_backend,
        world_size,
        env_per_rank,
        triattn_backend,
        create_boltz2_model_init_params(use_large_model=False),
    )


# ---------------------------------------------------------------------------
# SetAttnPairBiasBackend
# ---------------------------------------------------------------------------


def _parallel_assert_set_attn_pair_bias_backend(rank, env_per_rank, sdpa_backend, boltz2_params):
    """Worker: verify SetAttnPairBiasBackend targets only AttentionPairBias instances."""
    monkeypatch = pytest.MonkeyPatch()
    if env_per_rank is not None:
        for var_name, value in env_per_rank.items():
            monkeypatch.setenv(var_name, f"{rank}" if value == "<INPUT_RANK>" else value)

    from boltz.distributed.model.layers.attention import AttentionPairBias, AttentionPairBiasShardwise
    from boltz.distributed.model.models.boltz2 import Boltz2 as Boltz2Distributed
    from boltz.model.models.boltz2 import Boltz2 as SerialBoltz2

    grid_group_sizes = {"dp": 1, "cp": (1, 1)}
    DistributedManager.initialize(grid_group_sizes, device_type="cuda", backend="nccl")
    manager = DistributedManager()

    serial_model = SerialBoltz2(**boltz2_params).to(device=manager.device).eval()
    dist_model = Boltz2Distributed(serial_model, manager).eval()

    attn_modules_before = []
    for name, submodule in dist_model.named_modules():
        if isinstance(submodule, AttentionPairBias):
            attn_modules_before.append(name)
            assert submodule.sdpa_with_bias_backend != sdpa_backend, (
                f"{name}: default backend already matches {sdpa_backend}, "
                "test cannot verify the setter changes anything"
            )

    assert len(attn_modules_before) > 0, "Model must contain at least one AttentionPairBias"

    shardwise_backends_before = {}
    for name, submodule in dist_model.named_modules():
        if isinstance(submodule, AttentionPairBiasShardwise):
            shardwise_backends_before[name] = submodule.sdpa_with_bias_backend

    dist_model.apply(SetAttnPairBiasBackend(sdpa_backend))

    for name, submodule in dist_model.named_modules():
        if isinstance(submodule, AttentionPairBias):
            assert submodule.sdpa_with_bias_backend == sdpa_backend, (
                f"{name}: expected {sdpa_backend} after setter, " f"got {submodule.sdpa_with_bias_backend}"
            )
        if isinstance(submodule, AttentionPairBiasShardwise):
            assert submodule.sdpa_with_bias_backend == shardwise_backends_before[name], (
                f"AttentionPairBiasShardwise {name} was unexpectedly changed by " f"SetAttnPairBiasBackend"
            )

    DistributedManager.cleanup()


@pytest.mark.parametrize(
    "setup_env",
    [((1, (1, 1)), False, "cuda", "ENV")],
    indirect=True,
    ids=["cuda-dp1-cp1x1"],
)
@pytest.mark.parametrize(
    "sdpa_backend",
    [SDPAWithBiasBackend.TORCH_FLEX_ATTN],
    ids=lambda b: b.value,
)
def test_set_attn_pair_bias_backend(setup_env, sdpa_backend):
    """SetAttnPairBiasBackend sets sdpa_with_bias_backend only on AttentionPairBias instances."""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if not torch.cuda.is_available():
        pytest.skip("CUDA not available")
    if torch.cuda.device_count() < world_size:
        pytest.skip(f"Need {world_size} GPUs, have {torch.cuda.device_count()}")

    spawn_multiprocessing(
        _parallel_assert_set_attn_pair_bias_backend,
        world_size,
        env_per_rank,
        sdpa_backend,
        create_boltz2_model_init_params(use_large_model=False),
    )


# ---------------------------------------------------------------------------
# SetAttnPairBiasShardwiseBackend
# ---------------------------------------------------------------------------


def _parallel_assert_set_attn_pair_bias_shardwise_backend(rank, env_per_rank, sdpa_backend, boltz2_params):
    """Worker: verify SetAttnPairBiasShardwiseBackend targets only AttentionPairBiasShardwise instances."""
    monkeypatch = pytest.MonkeyPatch()
    if env_per_rank is not None:
        for var_name, value in env_per_rank.items():
            monkeypatch.setenv(var_name, f"{rank}" if value == "<INPUT_RANK>" else value)

    from boltz.distributed.model.layers.attention import AttentionPairBias, AttentionPairBiasShardwise
    from boltz.distributed.model.models.boltz2 import Boltz2 as Boltz2Distributed
    from boltz.model.models.boltz2 import Boltz2 as SerialBoltz2

    grid_group_sizes = {"dp": 1, "cp": (1, 1)}
    DistributedManager.initialize(grid_group_sizes, device_type="cuda", backend="nccl")
    manager = DistributedManager()

    serial_model = SerialBoltz2(**boltz2_params).to(device=manager.device).eval()
    dist_model = Boltz2Distributed(serial_model, manager).eval()

    shardwise_modules_before = []
    for name, submodule in dist_model.named_modules():
        if isinstance(submodule, AttentionPairBiasShardwise):
            shardwise_modules_before.append(name)
            assert submodule.sdpa_with_bias_backend != sdpa_backend, (
                f"{name}: default backend already matches {sdpa_backend}, "
                "test cannot verify the setter changes anything"
            )

    assert len(shardwise_modules_before) > 0, "Model must contain at least one AttentionPairBiasShardwise"

    ring_backends_before = {}
    for name, submodule in dist_model.named_modules():
        if isinstance(submodule, AttentionPairBias):
            ring_backends_before[name] = submodule.sdpa_with_bias_backend

    dist_model.apply(SetAttnPairBiasShardwiseBackend(sdpa_backend))

    for name, submodule in dist_model.named_modules():
        if isinstance(submodule, AttentionPairBiasShardwise):
            assert submodule.sdpa_with_bias_backend == sdpa_backend, (
                f"{name}: expected {sdpa_backend} after setter, " f"got {submodule.sdpa_with_bias_backend}"
            )
        if isinstance(submodule, AttentionPairBias):
            assert submodule.sdpa_with_bias_backend == ring_backends_before[name], (
                f"AttentionPairBias {name} was unexpectedly changed by " f"SetAttnPairBiasShardwiseBackend"
            )

    DistributedManager.cleanup()


@pytest.mark.parametrize(
    "setup_env",
    [((1, (1, 1)), False, "cuda", "ENV")],
    indirect=True,
    ids=["cuda-dp1-cp1x1"],
)
@pytest.mark.parametrize(
    "sdpa_backend",
    [SDPAWithBiasBackend.TORCH_FLEX_ATTN, SDPAWithBiasBackend.TORCH_SDPA_EFFICIENT_ATTENTION],
    ids=lambda b: b.value,
)
def test_set_attn_pair_bias_shardwise_backend(setup_env, sdpa_backend):
    """SetAttnPairBiasShardwiseBackend sets sdpa_with_bias_backend only on AttentionPairBiasShardwise instances."""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if not torch.cuda.is_available():
        pytest.skip("CUDA not available")
    if torch.cuda.device_count() < world_size:
        pytest.skip(f"Need {world_size} GPUs, have {torch.cuda.device_count()}")

    spawn_multiprocessing(
        _parallel_assert_set_attn_pair_bias_shardwise_backend,
        world_size,
        env_per_rank,
        sdpa_backend,
        create_boltz2_model_init_params(use_large_model=False),
    )


# ---------------------------------------------------------------------------
# OffloadActvCkptToCPU
# ---------------------------------------------------------------------------

_ALL_OFFLOAD_TYPES = ("DiffusionTransformer", "MSAModule", "PairformerModule")


def _parallel_assert_offload_actv_ckpt_to_cpu(rank, env_per_rank, target_names, boltz2_params):
    """Worker: verify OffloadActvCkptToCPU targets only the requested module types."""
    monkeypatch = pytest.MonkeyPatch()
    if env_per_rank is not None:
        for var_name, value in env_per_rank.items():
            monkeypatch.setenv(var_name, f"{rank}" if value == "<INPUT_RANK>" else value)

    from boltz.distributed.model.layers.pairformer import PairformerModule
    from boltz.distributed.model.models.boltz2 import Boltz2 as Boltz2Distributed
    from boltz.distributed.model.modules.transformers import DiffusionTransformer
    from boltz.distributed.model.modules.trunkv2 import MSAModule
    from boltz.model.models.boltz2 import Boltz2 as SerialBoltz2

    name_to_cls = {
        "DiffusionTransformer": DiffusionTransformer,
        "MSAModule": MSAModule,
        "PairformerModule": PairformerModule,
    }
    target_classes = tuple(name_to_cls[n] for n in target_names)
    all_classes = tuple(name_to_cls.values())

    grid_group_sizes = {"dp": 1, "cp": (1, 1)}
    DistributedManager.initialize(grid_group_sizes, device_type="cuda", backend="nccl")
    manager = DistributedManager()

    serial_model = SerialBoltz2(**boltz2_params).to(device=manager.device).eval()
    dist_model = Boltz2Distributed(serial_model, manager).eval()

    for name, submodule in dist_model.named_modules():
        if isinstance(submodule, all_classes):
            assert (
                not submodule.cpu_offloading
            ), f"{name} ({type(submodule).__name__}): cpu_offloading should be False before setter"
            submodule.activation_checkpointing = True

    dist_model.apply(OffloadActvCkptToCPU(set(target_names)))

    found_targeted = {cls: 0 for cls in target_classes}
    for name, submodule in dist_model.named_modules():
        if isinstance(submodule, target_classes):
            assert (
                submodule.cpu_offloading
            ), f"{name} ({type(submodule).__name__}): cpu_offloading should be True after setter"
            for cls in target_classes:
                if isinstance(submodule, cls):
                    found_targeted[cls] += 1
        elif isinstance(submodule, all_classes):
            assert not submodule.cpu_offloading, (
                f"{name} ({type(submodule).__name__}): non-targeted module should still have " f"cpu_offloading=False"
            )

    for cls, count in found_targeted.items():
        assert count > 0, f"Model must contain at least one {cls.__name__} but found none"

    DistributedManager.cleanup()


@pytest.mark.parametrize(
    "setup_env",
    [((1, (1, 1)), False, "cuda", "ENV")],
    indirect=True,
    ids=["cuda-dp1-cp1x1"],
)
@pytest.mark.parametrize(
    "target_names",
    [
        ["DiffusionTransformer"],
        ["MSAModule", "PairformerModule"],
        list(_ALL_OFFLOAD_TYPES),
    ],
    ids=["score_model_only", "msa_and_pairformer", "all_three"],
)
def test_offload_actv_ckpt_to_cpu(setup_env, target_names):
    """OffloadActvCkptToCPU sets cpu_offloading=True only on targeted module types."""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if not torch.cuda.is_available():
        pytest.skip("CUDA not available")
    if torch.cuda.device_count() < world_size:
        pytest.skip(f"Need {world_size} GPUs, have {torch.cuda.device_count()}")

    spawn_multiprocessing(
        _parallel_assert_offload_actv_ckpt_to_cpu,
        world_size,
        env_per_rank,
        target_names,
        create_boltz2_model_init_params(use_large_model=False),
    )


def test_offload_actv_ckpt_to_cpu_invalid_type():
    """OffloadActvCkptToCPU rejects unrecognised module type names."""
    with pytest.raises(ValueError, match="Invalid module type"):
        OffloadActvCkptToCPU({"InvalidModule"})


def test_offload_actv_ckpt_to_cpu_empty_list():
    """OffloadActvCkptToCPU rejects an empty module_types list."""
    with pytest.raises(ValueError, match="must be non-empty"):
        OffloadActvCkptToCPU(set())


def _parallel_assert_offload_actv_ckpt_rejects_no_ckpt(rank, env_per_rank, boltz2_params):
    """Worker: verify OffloadActvCkptToCPU raises when activation_checkpointing is off."""
    monkeypatch = pytest.MonkeyPatch()
    if env_per_rank is not None:
        for var_name, value in env_per_rank.items():
            monkeypatch.setenv(var_name, f"{rank}" if value == "<INPUT_RANK>" else value)

    from boltz.distributed.model.layers.pairformer import PairformerModule
    from boltz.distributed.model.models.boltz2 import Boltz2 as Boltz2Distributed
    from boltz.model.models.boltz2 import Boltz2 as SerialBoltz2

    grid_group_sizes = {"dp": 1, "cp": (1, 1)}
    DistributedManager.initialize(grid_group_sizes, device_type="cuda", backend="nccl")
    manager = DistributedManager()

    serial_model = SerialBoltz2(**boltz2_params).to(device=manager.device).eval()
    dist_model = Boltz2Distributed(serial_model, manager).eval()

    for _name, submodule in dist_model.named_modules():
        if isinstance(submodule, PairformerModule):
            submodule.activation_checkpointing = False
            break

    with pytest.raises(ValueError, match="activation_checkpointing is not enabled"):
        dist_model.apply(OffloadActvCkptToCPU({"PairformerModule"}))

    DistributedManager.cleanup()


@pytest.mark.parametrize(
    "setup_env",
    [((1, (1, 1)), False, "cuda", "ENV")],
    indirect=True,
    ids=["cuda-dp1-cp1x1"],
)
def test_offload_actv_ckpt_to_cpu_rejects_no_ckpt(setup_env):
    """OffloadActvCkptToCPU raises when activation_checkpointing is disabled on a target."""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if not torch.cuda.is_available():
        pytest.skip("CUDA not available")
    if torch.cuda.device_count() < world_size:
        pytest.skip(f"Need {world_size} GPUs, have {torch.cuda.device_count()}")

    spawn_multiprocessing(
        _parallel_assert_offload_actv_ckpt_rejects_no_ckpt,
        world_size,
        env_per_rank,
        create_boltz2_model_init_params(use_large_model=False),
    )
