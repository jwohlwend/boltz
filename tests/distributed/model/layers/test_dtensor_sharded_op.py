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


import pytest
import torch
from torch.distributed.tensor import Replicate, Shard, distribute_tensor

from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.layers.sharded_op import sharded_sum
from boltz.testing.utils import assert_tensors_identical, spawn_multiprocessing


def serial_sum(x: torch.Tensor, dims: tuple[int, ...] | int, keepdim: bool = False) -> torch.Tensor:
    """Serial implementation of sum operation for comparison."""
    return torch.sum(x, dim=dims, keepdim=keepdim)


def parallel_assert_sharded_sum(
    rank,
    grid_group_sizes,
    device_type,
    backend,
    env_per_rank,
    placements,
    dims,
    keepdim,
    input_tensor_global_host,
    output_expected_global_host,
    d_output_expected_global_host,
    d_input_expected_global_host,
    output_placements,
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

    device_mesh = manager.device_mesh_subgroups

    # Distribute input tensor
    input_tensor_dtensor = distribute_tensor(
        input_tensor_global_host.to(manager.device), device_mesh=device_mesh, placements=placements
    ).requires_grad_(True)

    # Distribute expected outputs
    d_output_expected_dtensor = distribute_tensor(
        d_output_expected_global_host.to(manager.device),
        device_mesh=device_mesh,
        placements=output_placements,
    )

    # Create copy to verify input isn't modified
    input_tensor_dtensor_copy = input_tensor_dtensor.detach().clone().requires_grad_(True)

    # Forward pass
    output_dtensor_result = sharded_sum(input_tensor_dtensor, dim=dims, keepdim=keepdim)

    # Verify input wasn't modified
    assert_tensors_identical(
        input_tensor_dtensor_copy.to_local(), input_tensor_dtensor.to_local(), check_grad=False, check_grad_fn=False
    )

    # Test forward pass results
    assert (
        output_dtensor_result.placements == output_placements
    ), f"output_placements: {output_placements}, output_dtensor_result.placements: {output_dtensor_result.placements}"
    assert (
        output_dtensor_result.shape == output_expected_global_host.shape
    ), f"Output shape mismatch: {output_dtensor_result.shape} != {output_expected_global_host.shape}"
    assert (
        output_dtensor_result.stride() == output_expected_global_host.stride()
    ), f"Output stride mismatch: {output_dtensor_result.stride()} != {output_expected_global_host.stride()}"
    torch.testing.assert_close(output_dtensor_result.full_tensor().cpu(), output_expected_global_host)

    # Backward pass
    d_output_expected_dtensor_copy = d_output_expected_dtensor.detach().clone()
    output_dtensor_result.backward(d_output_expected_dtensor)

    # Verify upstream gradient wasn't modified
    assert_tensors_identical(d_output_expected_dtensor_copy.to_local(), d_output_expected_dtensor.to_local())

    # Test input gradient
    assert (
        input_tensor_dtensor.grad.placements == placements
    ), f"placements: {placements}, input_tensor_dtensor.grad.placements: {input_tensor_dtensor.grad.placements}"
    assert (
        input_tensor_dtensor.grad.shape == input_tensor_global_host.shape
    ), f"Input gradient shape mismatch: {input_tensor_dtensor.grad.shape} != {input_tensor_global_host.shape}"
    torch.testing.assert_close(input_tensor_dtensor.grad.full_tensor().cpu(), d_input_expected_global_host)

    # Verify full tensors match expected results
    torch.testing.assert_close(output_dtensor_result.full_tensor().cpu(), output_expected_global_host)
    torch.testing.assert_close(input_tensor_dtensor.grad.full_tensor().cpu(), d_input_expected_global_host)

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env",
    [
        ((2, (2, 2)), True, "cuda", "ENV"),
    ],
    indirect=("setup_env",),
    ids=lambda x: f"dp={x[0][0]}, cp={x[0][1]}, specify_method={x[1]}, device_type={x[2]}, method_init={x[3]}",
)
@pytest.mark.parametrize(
    "placements,dims,keepdim,output_placements",
    [
        ((Shard(0), Shard(1), Shard(2)), (0, 1, 2), False, (Replicate(), Replicate(), Replicate())),
        ((Shard(0), Shard(1), Shard(2)), (0, 2), False, (Replicate(), Shard(0), Replicate())),
        ((Shard(0), Shard(1), Shard(2)), (1, 2), False, (Shard(0), Replicate(), Replicate())),
        ((Shard(0), Shard(1), Shard(2)), (1, 2), True, (Shard(0), Replicate(), Replicate())),
        ((Shard(0), Shard(1), Shard(2)), (1,), False, (Shard(0), Replicate(), Shard(1))),
        ((Shard(0), Shard(1), Shard(2)), (1,), True, (Shard(0), Replicate(), Shard(2))),
        ((Shard(0), Shard(1), Replicate()), (1,), False, (Shard(0), Replicate(), Replicate())),
        ((Shard(0), Shard(1), Replicate()), (1,), True, (Shard(0), Replicate(), Replicate())),
        ((Shard(2), Shard(0), Shard(1)), (1,), False, (Shard(1), Shard(0), Replicate())),
        ((Shard(2), Shard(0), Shard(1)), (1,), True, (Shard(2), Shard(0), Replicate())),
    ],
    ids=[
        "pair_dim_0_1_2",
        "pair_dim_0_2",
        "pair_dim_1_2_keepdim",
        "pair_dim_1_2",
        "pair_dim_1",
        "pair_dim_1_keepdim",
        "single_dim_1",
        "single_dim_1_keepdim",
        "s2_s0_s1_dim_1",
        "s2_s0_s1_dim_1_keepdim",
    ],
)
def test_sharded_sum_parallel(setup_env, placements, dims, keepdim, output_placements):
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    size_ring = grid_group_sizes["cp"][0]
    B = 2 * grid_group_sizes["dp"]
    N = size_ring * 3  # Number of tokens
    D = 5  # Hidden dimension

    # Create input tensor with proper shape
    input_tensor_global = torch.randn((B, N, N, D), requires_grad=True, device=device_type)

    # Run serial forward pass
    input_tensor_global_host = input_tensor_global.detach().clone().cpu()
    output_expected_global = serial_sum(input_tensor_global, dims=dims, keepdim=keepdim)
    output_expected_global_host = output_expected_global.detach().clone().cpu()

    # Create upstream gradient and run backward pass
    d_output_expected_global = torch.randn_like(output_expected_global)
    d_output_expected_global_host = d_output_expected_global.detach().clone().cpu()
    output_expected_global.backward(d_output_expected_global)

    # Serial sum grad is not set to upstream gradient
    assert not input_tensor_global.grad.is_set_to(d_output_expected_global)

    # Launch parallel test across all processes
    spawn_multiprocessing(
        parallel_assert_sharded_sum,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        placements,
        dims,
        keepdim,
        input_tensor_global_host,
        output_expected_global_host,
        d_output_expected_global_host,
        input_tensor_global.grad.detach().clone().cpu(),
        output_placements,
    )
