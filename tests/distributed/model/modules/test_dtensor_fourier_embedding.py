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

"""Tests for DTensor FourierEmbedding module.

Tests both Boltz-1x and Boltz-2 serial FourierEmbedding modules against the unified
DTensor FourierEmbedding implementation, verifying forward-only equivalence.
FourierEmbedding has frozen (non-trainable) parameters, so no backward test is needed.

The V1 and V2 serial FourierEmbedding implementations are identical in structure
and math; the test parametrizes over both to verify the DTensor wrapper accepts
either serial class.

Adapted from Boltz-1x CP test (tests_v1/distributed/model/modules/test_dtensor_encoders.py).
"""

import pytest
import torch
from torch.distributed.tensor import Replicate, Shard, distribute_tensor

from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.modules.encoders import FourierEmbedding as DTensorFourierEmbedding
from boltz.model.modules.encoders import FourierEmbedding as FourierEmbeddingBoltz1
from boltz.model.modules.encodersv2 import FourierEmbedding as FourierEmbeddingBoltz2
from boltz.testing.utils import (
    assert_all_identical,
    assert_tensors_identical,
    seed_by_rank,
    spawn_multiprocessing,
)


def parallel_assert_fourier_embedding(
    rank: int,
    grid_group_sizes,
    device_type: str,
    backend: str,
    env_per_rank,
    serial_class_tag: str,
    dim: int,
    layer_state_dict,
    input_global_host: torch.Tensor,
    output_expected_global_host: torch.Tensor,
):
    """Parallel assertion for FourierEmbedding forward pass."""
    monkeypatch = pytest.MonkeyPatch()
    if env_per_rank is not None:
        for var_name, value in env_per_rank.items():
            if value == "<INPUT_RANK>":
                monkeypatch.setenv(var_name, f"{rank}")
                continue
            monkeypatch.setenv(var_name, value)

    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    manager = DistributedManager()

    # Create serial module from state dict
    serial_class = FourierEmbeddingBoltz1 if serial_class_tag == "boltz1" else FourierEmbeddingBoltz2
    module_serial = serial_class(dim)
    module_serial = module_serial.to(device=manager.device)
    module_serial.load_state_dict(layer_state_dict)

    # Create DTensor module from serial
    module_dt = DTensorFourierEmbedding(module_serial, manager.device_mesh_subgroups)
    module_dt.train()

    # Placements: times is (B,) sharded along DP axis
    placements_times = (Shard(0), Replicate(), Replicate())

    # Distribute input
    input_dt = distribute_tensor(
        input_global_host.to(manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements_times,
    ).requires_grad_(False)

    input_dt_copy = input_dt.detach().clone()

    # Forward pass
    output_dt = module_dt(input_dt)

    # Verify input wasn't modified
    assert_tensors_identical(input_dt_copy.to_local(), input_dt.to_local(), check_grad=False, check_grad_fn=False)

    # Forward compare: local shard
    output_expected_dt = distribute_tensor(
        output_expected_global_host.to(manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements_times,
    )
    torch.testing.assert_close(output_dt.to_local(), output_expected_dt.to_local())

    # Verify all CP ranks produce identical output (frozen params, replicated computation)
    assert_all_identical(output_dt.to_local().detach(), manager.group["cp"])

    # Verify full tensor matches serial reference
    torch.testing.assert_close(output_dt.full_tensor().cpu(), output_expected_global_host)

    # Verify all parameters are frozen
    for name, param in module_dt.named_parameters():
        assert not param.requires_grad, f"Parameter {name} should be frozen but requires_grad=True"

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env",
    [
        ((2, (2, 2)), True, "cuda", "ENV"),
    ],
    indirect=("setup_env",),
    ids=lambda x: f"dp:{x[0][0]}, cp:{x[0][1]}, specify_method:{x[1]}, device_type:{x[2]}, method_init:{x[3]}",
)
@pytest.mark.parametrize(
    "serial_class_tag",
    ["boltz1", "boltz2"],
    ids=["v1", "v2"],
)
def test_fourier_embedding(setup_env, serial_class_tag: str):
    """Test FourierEmbedding DTensor vs serial equivalence for both V1 and V2.

    FourierEmbedding has frozen parameters (non-trainable), so this is a
    forward-only test. V1 and V2 serial implementations are identical; both
    are tested to verify the DTensor wrapper's isinstance check accepts either.
    """
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() < {world_size}")

    B = 2 * grid_group_sizes["dp"]
    dim = 256

    seed_by_rank(0, seed=42)

    serial_class = FourierEmbeddingBoltz1 if serial_class_tag == "boltz1" else FourierEmbeddingBoltz2

    # Create serial module — cast to device before forward
    module_serial = serial_class(dim)
    module_serial = module_serial.to(device=device_type)
    module_serial.train()
    layer_state_dict = module_serial.state_dict()

    # Create input
    input_global = torch.rand((B,), device=device_type, requires_grad=False)

    # Serial forward (no backward — frozen params)
    output_expected_global = module_serial(input_global)

    # Move to host for multiprocessing
    input_global_host = input_global.detach().cpu()
    output_expected_global_host = output_expected_global.detach().cpu()

    spawn_multiprocessing(
        parallel_assert_fourier_embedding,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        serial_class_tag,
        dim,
        layer_state_dict,
        input_global_host,
        output_expected_global_host,
    )
