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


import unittest

import pytest
import torch
from torch.distributed.tensor import Replicate, Shard, distribute_tensor

from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.layers.where import where
from boltz.testing.utils import assert_tensors_identical, spawn_multiprocessing


def serial_where(condition: torch.Tensor, x: torch.Tensor, y: torch.Tensor) -> torch.Tensor:
    """Serial implementation of where operation for comparison."""
    return torch.where(condition, x, y)


def parallel_assert_where(
    rank,
    grid_group_sizes,
    device_type,
    backend,
    env_per_rank,
    placements,
    condition_global_host,
    x_global_host,
    y_global_host,
    output_expected_global_host,
    d_output_expected_global_host,
    d_x_expected_global_host,
    d_y_expected_global_host,
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
    condition_dtensor = distribute_tensor(
        condition_global_host.to(manager.device), device_mesh=manager.device_mesh_subgroups, placements=placements
    )
    x_dtensor = distribute_tensor(
        x_global_host.to(manager.device), device_mesh=manager.device_mesh_subgroups, placements=placements
    ).requires_grad_(True)
    y_dtensor = distribute_tensor(
        y_global_host.to(manager.device), device_mesh=manager.device_mesh_subgroups, placements=placements
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
    d_x_expected_dtensor = distribute_tensor(
        d_x_expected_global_host.to(manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements,
        src_data_rank=None,
    )
    d_y_expected_dtensor = distribute_tensor(
        d_y_expected_global_host.to(manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements,
        src_data_rank=None,
    )

    # Create copies to verify inputs aren't modified
    condition_dtensor_copy = condition_dtensor.detach().clone()
    x_dtensor_copy = x_dtensor.detach().clone().requires_grad_(True)
    y_dtensor_copy = y_dtensor.detach().clone().requires_grad_(True)

    # Forward pass
    output_dtensor_result = where(condition_dtensor, x_dtensor, y_dtensor)

    # Verify inputs weren't modified
    assert_tensors_identical(
        condition_dtensor_copy.to_local(), condition_dtensor.to_local(), check_grad=False, check_grad_fn=False
    )
    assert_tensors_identical(x_dtensor_copy.to_local(), x_dtensor.to_local(), check_grad=False, check_grad_fn=False)
    assert_tensors_identical(y_dtensor_copy.to_local(), y_dtensor.to_local(), check_grad=False, check_grad_fn=False)

    # Test forward pass results
    assert (
        output_dtensor_result.shape == output_expected_dtensor.shape
    ), f"Output shape mismatch: {output_dtensor_result.shape} != {output_expected_dtensor.shape}"
    assert (
        output_dtensor_result.stride() == output_expected_dtensor.stride()
    ), f"Output stride mismatch: {output_dtensor_result.stride()} != {output_expected_dtensor.stride()}"
    torch.testing.assert_close(output_dtensor_result.to_local(), output_expected_dtensor.to_local())

    # Backward pass
    d_output_expected_dtensor_copy = d_output_expected_dtensor.detach().clone()
    output_dtensor_result.backward(d_output_expected_dtensor)

    # Verify upstream gradient wasn't modified
    assert_tensors_identical(d_output_expected_dtensor_copy.to_local(), d_output_expected_dtensor.to_local())

    # Test input gradients
    assert (
        x_dtensor.grad.shape == d_x_expected_dtensor.shape
    ), f"Input gradient shape mismatch: {x_dtensor.grad.shape} != {d_x_expected_dtensor.shape}"
    assert (
        x_dtensor.grad.stride() == d_x_expected_dtensor.stride()
    ), f"Input gradient stride mismatch: {x_dtensor.grad.stride()} != {d_x_expected_dtensor.stride()}"
    torch.testing.assert_close(x_dtensor.grad.to_local(), d_x_expected_dtensor.to_local())

    assert (
        y_dtensor.grad.shape == d_y_expected_dtensor.shape
    ), f"Input gradient shape mismatch: {y_dtensor.grad.shape} != {d_y_expected_dtensor.shape}"
    assert (
        y_dtensor.grad.stride() == d_y_expected_dtensor.stride()
    ), f"Input gradient stride mismatch: {y_dtensor.grad.stride()} != {d_y_expected_dtensor.stride()}"
    torch.testing.assert_close(y_dtensor.grad.to_local(), d_y_expected_dtensor.to_local())

    # Test full tensor gathering - verify distributed results match serial results
    output_global_result_host = output_dtensor_result.full_tensor().cpu()
    d_x_global_result_host = x_dtensor.grad.full_tensor().cpu()
    d_y_global_result_host = y_dtensor.grad.full_tensor().cpu()

    # Verify full tensors match expected results
    torch.testing.assert_close(output_global_result_host, output_expected_global_host)
    torch.testing.assert_close(d_x_global_result_host, d_x_expected_global_host)
    torch.testing.assert_close(d_y_global_result_host, d_y_expected_global_host)

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env",
    [
        ((1, (2, 2)), True, "cuda", "ENV"),
        ((2, (2, 2)), True, "cuda", "ENV"),
        ((1, (3, 3)), True, "cuda", "ENV"),
        ((1, (3, 3)), True, "cpu", "ENV"),
    ],
    indirect=("setup_env",),
    ids=lambda x: f"dp:{x[0][0]}, cp:{x[0][1]}, specify_method:{x[1]}, device_type:{x[2]}, method_init:{x[3]}",
)
@pytest.mark.parametrize(
    "placements", [(Shard(0), Shard(1), Shard(2)), (Shard(0), Shard(1), Replicate())], ids=["shard", "replicate"]
)
@pytest.mark.parametrize(
    "condition_type",
    ["random", "threshold"],
)
def test_where_parallel(setup_env, placements, condition_type):
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    size_ring = grid_group_sizes["cp"][0]
    B = 2 * grid_group_sizes["dp"]
    N = size_ring * 100  # Number of tokens
    D = 32  # Hidden dimension

    seed = 42
    rng = torch.Generator(device=device_type)
    rng.manual_seed(seed)

    # Create input tensors
    x_global = torch.empty((B, N, N, D), requires_grad=True, device=device_type)
    y_global = torch.empty((B, N, N, D), requires_grad=True, device=device_type)
    with torch.no_grad():
        x_global.uniform_(-5, 5, generator=rng)
        y_global.uniform_(-3, 7, generator=rng)

    # Create different types of conditions for testing
    if condition_type == "random":
        condition_global = torch.randint(0, 2, (B, N, N, D), dtype=torch.bool, device=device_type, generator=rng)
    elif condition_type == "threshold":
        # Create condition based on x values (test gradient flow)
        condition_global = x_global > 0.0
    else:
        raise ValueError(f"Invalid condition type: {condition_type}")

    # Run serial forward pass
    condition_global_host = condition_global.detach().clone().cpu()
    x_global_host = x_global.detach().clone().cpu()
    y_global_host = y_global.detach().clone().cpu()
    output_expected_global = serial_where(condition_global, x_global, y_global)
    output_expected_global_host = output_expected_global.detach().clone().cpu()

    # Create upstream gradient and run backward pass
    d_output_expected_global = torch.rand_like(output_expected_global)
    d_output_expected_global_host = d_output_expected_global.detach().clone().cpu()
    output_expected_global.backward(d_output_expected_global)

    # Launch parallel test across all processes
    spawn_multiprocessing(
        parallel_assert_where,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        placements,
        condition_global_host,
        x_global_host,
        y_global_host,
        output_expected_global_host,
        d_output_expected_global_host,
        x_global.grad.detach().clone().cpu(),
        y_global.grad.detach().clone().cpu(),
    )


if __name__ == "__main__":
    unittest.main()
