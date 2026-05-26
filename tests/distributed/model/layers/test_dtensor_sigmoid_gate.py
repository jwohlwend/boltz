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


import itertools

import pytest
import torch
from torch.distributed.tensor import Replicate, Shard, distribute_tensor

from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.layers.sigmoid_gate import sigmoid_gate
from boltz.testing.utils import assert_tensors_identical, spawn_multiprocessing


def serial_sigmoid_gate(x: torch.Tensor, g: torch.Tensor) -> torch.Tensor:
    """Serial implementation of sigmoid gate for comparison."""
    return x * g.sigmoid()


def parallel_assert_sigmoid_gate(
    rank,
    grid_group_sizes,
    device_type,
    backend,
    env_per_rank,
    placements,
    input_x_global_host,
    input_g_global_host,
    output_expected_global_host,
    d_output_expected_global_host,
    d_input_x_expected_global_host,
    d_input_g_expected_global_host,
):
    monkeypatch = pytest.MonkeyPatch()
    if env_per_rank is not None:
        for var_name, value in env_per_rank.items():
            if value == "<INPUT_RANK>":
                monkeypatch.setenv(var_name, f"{rank}")
                continue
            monkeypatch.setenv(var_name, value)

    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    manager = DistributedManager()

    # Distribute input tensors
    input_x_dtensor = distribute_tensor(
        input_x_global_host.to(manager.device), device_mesh=manager.device_mesh_subgroups, placements=placements
    ).requires_grad_(True)

    input_g_dtensor = distribute_tensor(
        input_g_global_host.to(manager.device), device_mesh=manager.device_mesh_subgroups, placements=placements
    ).requires_grad_(True)

    # Distribute expected outputs
    d_output_expected_dtensor = distribute_tensor(
        d_output_expected_global_host.to(manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements,
    )
    output_expected_dtensor = distribute_tensor(
        output_expected_global_host.to(manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements,
        src_data_rank=None,
    )
    d_input_x_expected_dtensor = distribute_tensor(
        d_input_x_expected_global_host.to(manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements,
        src_data_rank=None,
    )
    d_input_g_expected_dtensor = distribute_tensor(
        d_input_g_expected_global_host.to(manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements,
        src_data_rank=None,
    )

    # Create copies to verify inputs aren't modified
    input_x_dtensor_copy = input_x_dtensor.detach().clone().requires_grad_(True)
    input_g_dtensor_copy = input_g_dtensor.detach().clone().requires_grad_(True)

    # Forward pass
    output_dtensor_result = sigmoid_gate(input_x_dtensor, input_g_dtensor)

    # Verify inputs weren't modified
    assert_tensors_identical(
        input_x_dtensor_copy.to_local(), input_x_dtensor.to_local(), check_grad=False, check_grad_fn=False
    )
    assert_tensors_identical(
        input_g_dtensor_copy.to_local(), input_g_dtensor.to_local(), check_grad=False, check_grad_fn=False
    )

    # Test forward pass results
    torch.testing.assert_close(output_dtensor_result.to_local(), output_expected_dtensor.to_local())

    # Backward pass
    d_output_expected_dtensor_copy = d_output_expected_dtensor.detach().clone()
    output_dtensor_result.backward(d_output_expected_dtensor)

    # Verify upstream gradient wasn't modified
    assert_tensors_identical(d_output_expected_dtensor_copy.to_local(), d_output_expected_dtensor.to_local())

    # Test input gradients
    torch.testing.assert_close(input_x_dtensor.grad.to_local(), d_input_x_expected_dtensor.to_local())
    torch.testing.assert_close(input_g_dtensor.grad.to_local(), d_input_g_expected_dtensor.to_local())

    # Test full tensor gathering - verify distributed results match serial results
    output_global_result_host = output_dtensor_result.full_tensor().cpu()
    d_input_x_global_result_host = input_x_dtensor.grad.full_tensor().cpu()
    d_input_g_global_result_host = input_g_dtensor.grad.full_tensor().cpu()

    # Verify full tensors match expected results
    torch.testing.assert_close(output_global_result_host, output_expected_global_host)
    torch.testing.assert_close(d_input_x_global_result_host, d_input_x_expected_global_host)
    torch.testing.assert_close(d_input_g_global_result_host, d_input_g_expected_global_host)

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env",
    itertools.product([(1, (2, 2)), (2, (2, 2))], [True], ["cpu", "cuda"], ["ENV"]),
    indirect=("setup_env",),
    ids=lambda x: f"dp={x[0][0]}, cp={x[0][1]}, specify_method={x[1]}, device_type={x[2]}, method_init={x[3]}",
)
@pytest.mark.parametrize(
    "placements", [(Shard(0), Shard(1), Shard(2)), (Shard(0), Shard(1), Replicate())], ids=["shard", "replicate"]
)
def test_sigmoid_gate_parallel(setup_env, placements):
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    size_ring = grid_group_sizes["cp"][0]
    B = 2 * grid_group_sizes["dp"]
    N = size_ring * 2  # Number of tokens
    D = 8  # Hidden dimension

    seed = 42
    rng = torch.Generator(device=device_type)
    rng.manual_seed(seed)

    # Create input tensors with proper shapes
    input_x_global = torch.rand((B, N, N, D), generator=rng, requires_grad=True, device=device_type)
    input_g_global = torch.rand((B, N, N, D), generator=rng, requires_grad=True, device=device_type)

    # Run serial forward pass
    input_x_global_host = input_x_global.detach().clone().cpu()
    input_g_global_host = input_g_global.detach().clone().cpu()
    output_expected_global = serial_sigmoid_gate(input_x_global, input_g_global)
    output_expected_global_host = output_expected_global.detach().clone().cpu()

    # Create upstream gradient and run backward pass
    d_output_expected_global = torch.rand_like(output_expected_global)
    d_output_expected_global_host = d_output_expected_global.detach().clone().cpu()
    output_expected_global.backward(d_output_expected_global)

    # Launch parallel test across all processes
    spawn_multiprocessing(
        parallel_assert_sigmoid_gate,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        placements,
        input_x_global_host,
        input_g_global_host,
        output_expected_global_host,
        d_output_expected_global_host,
        input_x_global.grad.detach().clone().cpu(),
        input_g_global.grad.detach().clone().cpu(),
    )
