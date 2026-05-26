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
from typing import Optional

import pytest
import torch
import torch.nn.functional as F
from torch.distributed.tensor import Replicate, Shard, distribute_tensor

from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.layers.shardwise_op import (
    ShardwiseOuterOp,
    shardwise_argmax,
    shardwise_log_softmax,
    shardwise_offset,
    shardwise_one_hot,
    shardwise_outer_op,
    shardwise_softmax,
    shardwise_sum,
)
from boltz.testing.utils import assert_tensors_identical, init_tensors_uniform, seed_by_rank, spawn_multiprocessing


def serial_shardwise_sum(x: torch.Tensor, dim: int, keepdim: Optional[bool] = None) -> torch.Tensor:
    """Serial implementation of shardwise sum operation for comparison."""
    if keepdim is None:
        return torch.sum(x, dim=dim)
    else:
        return torch.sum(x, dim=dim, keepdim=keepdim)


def serial_shardwise_one_hot(input: torch.Tensor, num_classes: int = -1) -> torch.Tensor:
    """Serial implementation of shardwise one_hot operation for comparison."""
    return F.one_hot(input, num_classes=num_classes)


def serial_shardwise_softmax(x: torch.Tensor, dim: int = -1) -> torch.Tensor:
    """Serial implementation of shardwise softmax operation for comparison."""
    return F.softmax(x, dim=dim)


def serial_shardwise_log_softmax(x: torch.Tensor, dim: int = -1) -> torch.Tensor:
    """Serial implementation of shardwise log_softmax operation for comparison."""
    return F.log_softmax(x, dim=dim)


def serial_shardwise_argmax(x: torch.Tensor, dim: int, keepdim: Optional[bool] = None) -> torch.Tensor:
    """Serial implementation of shardwise argmax operation for comparison."""
    if keepdim is None:
        return torch.argmax(x, dim=dim)
    return torch.argmax(x, dim=dim, keepdim=keepdim)


def parallel_assert_shardwise_sum(
    rank,
    grid_group_sizes,
    device_type,
    backend,
    env_per_rank,
    placements,
    input_x_global_host,
    output_expected_global_host,
    d_output_expected_global_host,
    d_input_x_expected_global_host,
    dim: int,
    keepdim: Optional[bool] = None,
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

    # Distribute input tensor
    input_x_dtensor = distribute_tensor(
        input_x_global_host.to(manager.device), device_mesh=manager.device_mesh_subgroups, placements=placements
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

    # Create copy to verify input isn't modified
    input_x_dtensor_copy = input_x_dtensor.detach().clone().requires_grad_(True)

    # Forward pass
    output_dtensor_result = shardwise_sum(input_x_dtensor, dim, keepdim)

    # Verify input wasn't modified
    assert_tensors_identical(
        input_x_dtensor_copy.to_local(), input_x_dtensor.to_local(), check_grad=False, check_grad_fn=False
    )

    # Test forward pass results
    assert output_dtensor_result.shape == output_expected_dtensor.shape
    assert output_dtensor_result.stride() == output_expected_dtensor.stride()
    torch.testing.assert_close(output_dtensor_result.to_local(), output_expected_dtensor.to_local())

    # Backward pass
    d_output_expected_dtensor_copy = d_output_expected_dtensor.detach().clone()
    output_dtensor_result.backward(d_output_expected_dtensor)

    # Verify upstream gradient wasn't modified
    assert_tensors_identical(d_output_expected_dtensor_copy.to_local(), d_output_expected_dtensor.to_local())

    # Test input gradient
    assert input_x_dtensor.grad.shape == d_input_x_expected_dtensor.shape
    assert input_x_dtensor.grad.stride() == d_input_x_expected_dtensor.stride()
    torch.testing.assert_close(input_x_dtensor.grad.to_local(), d_input_x_expected_dtensor.to_local())

    # Test full tensor gathering - verify distributed results match serial results
    output_global_result_host = output_dtensor_result.full_tensor().cpu()
    d_input_x_global_result_host = input_x_dtensor.grad.full_tensor().cpu()

    # Verify full tensors match expected results
    torch.testing.assert_close(output_global_result_host, output_expected_global_host)
    torch.testing.assert_close(d_input_x_global_result_host, d_input_x_expected_global_host)

    DistributedManager.cleanup()
    monkeypatch.undo()


def parallel_assert_shardwise_argmax(
    rank,
    grid_group_sizes,
    device_type,
    backend,
    env_per_rank,
    placements,
    input_x_global_host,
    output_expected_global_host,
    dim: int,
    keepdim: Optional[bool],
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

    input_x_dtensor = distribute_tensor(
        input_x_global_host.to(manager.device), device_mesh=manager.device_mesh_subgroups, placements=placements
    )

    output_expected_dtensor = distribute_tensor(
        output_expected_global_host.to(manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements,
        src_data_rank=None,
    )

    input_x_dtensor_copy = input_x_dtensor.detach().clone()

    output_dtensor_result = shardwise_argmax(input_x_dtensor, dim, keepdim)

    assert_tensors_identical(
        input_x_dtensor_copy.to_local(), input_x_dtensor.to_local(), check_grad=False, check_grad_fn=False
    )

    assert output_dtensor_result.shape == output_expected_dtensor.shape
    assert output_dtensor_result.stride() == output_expected_dtensor.stride()
    torch.testing.assert_close(output_dtensor_result.to_local(), output_expected_dtensor.to_local())

    output_global_result_host = output_dtensor_result.full_tensor().cpu()
    torch.testing.assert_close(output_global_result_host, output_expected_global_host)

    assert output_dtensor_result.placements == placements
    assert output_dtensor_result.dtype == torch.long

    with pytest.raises(RuntimeError, match="does not require grad"):
        output_dtensor_result.backward(torch.empty_like(output_dtensor_result))

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
    ids=lambda x: f"dp={x[0][0]}, cp={x[0][1]}, specify_method={x[1]}, device_type={x[2]}, method_init={x[3]}",
)
@pytest.mark.parametrize(
    "placements", [(Shard(0), Shard(1), Shard(2)), (Shard(0), Shard(1), Replicate())], ids=["shard", "replicate"]
)
@pytest.mark.parametrize(
    "sum_config",
    [
        (-1, False),
        (-1, True),
        (2, False),
        (2, True),
    ],
    ids=lambda x: f"dim={x[0]}, keepdim={x[1]}",
)
def test_shardwise_sum_parallel(setup_env, placements, sum_config):
    dim, keepdim = sum_config
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

    # Skip test if trying to sum along a sharded dimension
    input_shape = (B, N, N, D)
    actual_dim = dim if dim >= 0 else len(input_shape) + dim
    for placement in placements:
        if isinstance(placement, Shard) and placement.dim == actual_dim:
            pytest.skip(f"Skipping test: sum along sharded dimension {dim} is not supported")

    seed = 42
    rng = torch.Generator(device=device_type)
    rng.manual_seed(seed)

    # Create input tensor with proper shape
    input_x_global = torch.rand((B, N, N, D), requires_grad=True, device=device_type, generator=rng)

    # Run serial forward pass
    input_x_global_host = input_x_global.detach().clone().cpu()
    output_expected_global = serial_shardwise_sum(input_x_global, dim, keepdim)
    output_expected_global_host = output_expected_global.detach().clone().cpu()

    # Create upstream gradient and run backward pass
    d_output_expected_global = torch.rand_like(output_expected_global)
    d_output_expected_global_host = d_output_expected_global.detach().clone().cpu()
    output_expected_global.backward(d_output_expected_global)

    # Launch parallel test across all processes
    spawn_multiprocessing(
        parallel_assert_shardwise_sum,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        placements,
        input_x_global_host,
        output_expected_global_host,
        d_output_expected_global_host,
        input_x_global.grad.detach().clone().cpu(),
        dim,
        keepdim,
    )


def assert_error_case(rank, grid_group_sizes, device_type, backend, env_per_rank):
    monkeypatch = pytest.MonkeyPatch()
    if env_per_rank is not None:
        for var_name, value in env_per_rank.items():
            if value == "<INPUT_RANK>":
                monkeypatch.setenv(var_name, f"{rank}")
                continue
            monkeypatch.setenv(var_name, value)

    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    manager = DistributedManager()

    seed_by_rank(0)

    B = 2 * grid_group_sizes["dp"]
    N = grid_group_sizes["cp"][0] * 100  # Number of tokens
    D = 32  # Hidden dimension

    # Test case 1: Sum along sharded dimension should raise ValueError
    input_tensor = torch.randn((B, N, D), device=manager.device, requires_grad=True)
    sharded_dtensor = distribute_tensor(
        input_tensor, device_mesh=manager.device_mesh_subgroups, placements=(Shard(0), Shard(1), Replicate())
    )

    # This should raise an error because we're trying to sum along dimension 0, which is sharded
    with pytest.raises(ValueError, match="Sum along sharded dimension 0 is not supported"):
        shardwise_sum(sharded_dtensor, dim=0)
        shardwise_sum(sharded_dtensor, dim=1)

    # Test case 2: Invalid input type
    with pytest.raises(TypeError, match="Expected DTensor"):
        shardwise_sum(input_tensor, dim=1)  # Regular tensor instead of DTensor

    # Test case 3: Invalid dim type
    with pytest.raises(TypeError, match="Expected int for dim"):
        shardwise_sum(sharded_dtensor, dim=1.5)  # Float instead of int

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env",
    [
        ((1, (2, 2)), True, "cuda", "ENV"),
        ((1, (3, 3)), True, "cpu", "ENV"),
    ],
    indirect=("setup_env",),
    ids=lambda x: f"dp={x[0][0]}, cp={x[0][1]}, specify_method={x[1]}, device_type={x[2]}, method_init={x[3]}",
)
def test_shardwise_sum_error_cases(setup_env):
    """Test error cases for shardwise_sum function."""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    spawn_multiprocessing(
        assert_error_case,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
    )


def parallel_assert_shardwise_one_hot(
    rank,
    grid_group_sizes,
    device_type,
    backend,
    env_per_rank,
    placements,
    input_indices_global_host,
    output_expected_global_host,
    num_classes: int = -1,
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

    # Distribute input tensor (indices should be long type, no gradients needed)
    input_indices_dtensor = distribute_tensor(
        input_indices_global_host.to(manager.device), device_mesh=manager.device_mesh_subgroups, placements=placements
    )

    # Distribute expected output
    # For one_hot, the output has the same placements as input for the original dimensions,
    # and the new one-hot dimension will be handled by the shardwise_one_hot implementation
    output_expected_dtensor = distribute_tensor(
        output_expected_global_host.to(manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements,
        src_data_rank=None,
    )

    # Create copy to verify input isn't modified
    input_indices_dtensor_copy = input_indices_dtensor.detach().clone()

    # Forward pass
    output_dtensor_result = shardwise_one_hot(input_indices_dtensor, num_classes)

    # Verify input wasn't modified
    assert_tensors_identical(
        input_indices_dtensor_copy.to_local(), input_indices_dtensor.to_local(), check_grad=False, check_grad_fn=False
    )

    # Test forward pass results
    assert (
        output_dtensor_result.shape == output_expected_dtensor.shape
    ), f"Output shape mismatch: expected {output_expected_dtensor.shape}, got {output_dtensor_result.shape}"
    assert (
        output_dtensor_result.stride() == output_expected_dtensor.stride()
    ), f"Output stride mismatch: expected {output_expected_dtensor.stride()}, got {output_dtensor_result.stride()}"
    torch.testing.assert_close(output_dtensor_result.to_local(), output_expected_dtensor.to_local())

    # Test full tensor gathering - verify distributed results match serial results
    output_global_result_host = output_dtensor_result.full_tensor().cpu()

    # Verify full tensors match expected results
    torch.testing.assert_close(output_global_result_host, output_expected_global_host)

    # Verify placements are correct (same as input placements)
    assert output_dtensor_result.placements == placements

    # Ensure no backward possible
    with pytest.raises(RuntimeError, match="tensors does not require grad and does not have a grad_fn"):
        output_dtensor_result.backward(torch.empty_like(output_dtensor_result))

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
    "placements", [(Shard(0), Shard(1), Replicate()), (Shard(0), Shard(1), Shard(2))], ids=["single", "pair"]
)
@pytest.mark.parametrize(
    "num_classes_config",
    [
        3,  # explicit num_classes
        -1,  # inferred num_classes
    ],
    ids=lambda x: f"num_classes={x}",
)
def test_shardwise_one_hot_parallel(setup_env, placements, num_classes_config):
    num_classes = num_classes_config
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    size_ring = grid_group_sizes["cp"][0]
    B = 2 * grid_group_sizes["dp"]
    N = size_ring * 10  # Number of tokens (smaller for indices)
    D = 8  # Small dimension for indices

    seed = 42
    rng = torch.Generator(device=device_type)
    rng.manual_seed(seed)

    # Create input tensor with integer indices for one-hot encoding (3D to match placements)
    # Use modest range to ensure valid indices
    max_index = 2 if num_classes == -1 else num_classes - 1
    input_indices_global = torch.randint(
        0, max_index + 1, (B, N, N, D), device=device_type, generator=rng, dtype=torch.long
    )

    # Run serial forward pass
    input_indices_global_host = input_indices_global.detach().clone().cpu()
    output_expected_global = serial_shardwise_one_hot(input_indices_global, num_classes)
    output_expected_global_host = output_expected_global.detach().clone().cpu()

    # Ensure no backward possible
    with pytest.raises(RuntimeError, match="tensors does not require grad and does not have a grad_fn"):
        output_expected_global.backward(torch.empty_like(output_expected_global))

    # Launch parallel test across all processes
    spawn_multiprocessing(
        parallel_assert_shardwise_one_hot,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        placements,
        input_indices_global_host,
        output_expected_global_host,
        num_classes,
    )


def parallel_assert_shardwise_softmax(
    rank,
    grid_group_sizes,
    device_type,
    backend,
    env_per_rank,
    placements,
    input_x_global_host,
    output_expected_global_host,
    d_output_expected_global_host,
    d_input_x_expected_global_host,
    dim: int = -1,
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

    # Distribute input tensor
    input_x_dtensor = distribute_tensor(
        input_x_global_host.to(manager.device), device_mesh=manager.device_mesh_subgroups, placements=placements
    ).requires_grad_(True)

    # Distribute expected outputs
    d_output_expected_dtensor = distribute_tensor(
        d_output_expected_global_host.to(manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements,
    )

    # Create copy to verify input isn't modified
    input_x_dtensor_copy = input_x_dtensor.detach().clone().requires_grad_(True)

    # Forward pass
    output_dtensor_result = shardwise_softmax(input_x_dtensor, dim)

    # Verify input wasn't modified
    assert_tensors_identical(
        input_x_dtensor_copy.to_local(), input_x_dtensor.to_local(), check_grad=False, check_grad_fn=False
    )

    # Test forward pass results
    assert output_dtensor_result.shape == output_expected_global_host.shape
    assert output_dtensor_result.stride() == output_expected_global_host.stride()
    torch.testing.assert_close(output_dtensor_result.full_tensor().cpu(), output_expected_global_host.cpu())

    # Backward pass
    d_output_expected_dtensor_copy = d_output_expected_dtensor.detach().clone()
    output_dtensor_result.backward(d_output_expected_dtensor)

    # Verify upstream gradient wasn't modified
    assert_tensors_identical(d_output_expected_dtensor_copy.to_local(), d_output_expected_dtensor.to_local())

    # Test input gradient
    assert input_x_dtensor.grad.shape == d_input_x_expected_global_host.shape
    assert input_x_dtensor.grad.stride() == d_input_x_expected_global_host.stride()
    torch.testing.assert_close(input_x_dtensor.grad.full_tensor().cpu(), d_input_x_expected_global_host.cpu())

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env, dtype",
    (
        params_test := [
            (((2, (2, 2)), True, "cuda", "ENV"), torch.float32),
            (((2, (2, 2)), True, "cuda", "ENV"), torch.bfloat16),
        ]
    ),
    indirect=["setup_env"],
    ids=[
        f"dp={x[0][0][0]}, cp={x[0][0][1]}, specify_method={x[0][1]}, device_type={x[0][2]}, method_init={x[0][3]}, dtype={x[1]}"
        for x in params_test
    ],
)
@pytest.mark.parametrize(
    "dim_config",
    [
        -1,  # last dimension
        2,  # third dimension
    ],
    ids=lambda x: f"dim={x}",
)
def test_shardwise_softmax_parallel(setup_env, dtype, dim_config):
    dim = dim_config
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    placements = (Shard(0), Shard(1), Replicate())
    size_ring = grid_group_sizes["cp"][0]
    B = 2 * grid_group_sizes["dp"]
    N = size_ring * 10  # Number of tokens/atoms/etc for single repr input
    num_bins = 50  # Number of bins for softmax used in plddt_logits, see ConfidenceHeads class

    min_init_val = -0.2
    max_init_val = 0.2

    # Skip test if trying to apply softmax along a sharded dimension
    input_shape = (B, N, num_bins)

    seed = 42
    seed_by_rank(0, seed=seed)

    # Create input tensor with proper shape and dtype
    input_x_global = torch.empty(input_shape, requires_grad=True, device=device_type, dtype=dtype)
    init_tensors_uniform([input_x_global], low=min_init_val, high=max_init_val)

    # Run serial forward pass
    input_x_global_host = input_x_global.detach().clone().cpu()
    output_expected_global = serial_shardwise_softmax(input_x_global, dim)
    output_expected_global_host = output_expected_global.detach().clone().cpu()

    # Create upstream gradient and run backward pass
    d_output_expected_global = torch.empty_like(output_expected_global, dtype=dtype)
    init_tensors_uniform([d_output_expected_global], low=min_init_val, high=max_init_val)
    d_output_expected_global_host = d_output_expected_global.detach().clone().cpu()
    output_expected_global.backward(d_output_expected_global)

    # Launch parallel test across all processes
    spawn_multiprocessing(
        parallel_assert_shardwise_softmax,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        placements,
        input_x_global_host,
        output_expected_global_host,
        d_output_expected_global_host,
        input_x_global.grad.detach().clone().cpu(),
        dim,
    )


def parallel_assert_shardwise_softmax_error_cases(
    rank,
    grid_group_sizes,
    device_type,
    backend,
    env_per_rank,
    placements,
    input_x_global_host,
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

    # Distribute input tensor
    input_x_dtensor = distribute_tensor(
        input_x_global_host.to(manager.device), device_mesh=manager.device_mesh_subgroups, placements=placements
    ).requires_grad_(True)

    # Expect ValueError when softmax is along a sharded dimension - pick last one for simplicity
    with pytest.raises(ValueError, match="Softmax along sharded dimension"):
        shardwise_softmax(input_x_dtensor, dim=-1)

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env",
    (
        params_test := [
            ((1, (1, 1)), True, "cuda", "ENV"),
        ]
    ),
    indirect=["setup_env"],
    ids=[
        f"dp={x[0][0]}, cp={x[0][1]}, specify_method={x[1]}, device_type={x[2]}, method_init={x[3]}"
        for x in params_test
    ],
)
def test_shardwise_softmax_parallel_error_cases(setup_env):
    """Test that shardwise_softmax raises ValueError when softmax is along a sharded dimension."""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    placements = (Shard(0), Shard(1), Shard(2))
    size_ring = grid_group_sizes["cp"][0]
    B = 2 * grid_group_sizes["dp"]
    N = size_ring * 10
    num_bins = 50

    min_init_val = -0.2
    max_init_val = 0.2

    input_shape = (B, N, num_bins)

    seed = 42
    seed_by_rank(0, seed=seed)

    # Create input tensor with proper shape and dtype
    input_x_global = torch.empty(input_shape, requires_grad=True, device=device_type)
    init_tensors_uniform([input_x_global], low=min_init_val, high=max_init_val)

    input_x_global_host = input_x_global.detach().clone().cpu()

    # Launch parallel test across all processes
    spawn_multiprocessing(
        parallel_assert_shardwise_softmax_error_cases,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        placements,
        input_x_global_host,
    )


def parallel_assert_shardwise_log_softmax(
    rank,
    grid_group_sizes,
    device_type,
    backend,
    env_per_rank,
    placements,
    input_x_global_host,
    output_expected_global_host,
    d_output_expected_global_host,
    d_input_x_expected_global_host,
    dim: int = -1,
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

    # Distribute input tensor
    input_x_dtensor = distribute_tensor(
        input_x_global_host.to(manager.device), device_mesh=manager.device_mesh_subgroups, placements=placements
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

    # Create copy to verify input isn't modified
    input_x_dtensor_copy = input_x_dtensor.detach().clone().requires_grad_(True)

    # Forward pass
    output_dtensor_result = shardwise_log_softmax(input_x_dtensor, dim)

    # Verify input wasn't modified
    assert_tensors_identical(
        input_x_dtensor_copy.to_local(), input_x_dtensor.to_local(), check_grad=False, check_grad_fn=False
    )

    # Test forward pass results
    assert output_dtensor_result.shape == output_expected_dtensor.shape
    assert output_dtensor_result.stride() == output_expected_dtensor.stride()
    torch.testing.assert_close(output_dtensor_result.to_local(), output_expected_dtensor.to_local())

    # Backward pass
    d_output_expected_dtensor_copy = d_output_expected_dtensor.detach().clone()
    output_dtensor_result.backward(d_output_expected_dtensor)

    # Verify upstream gradient wasn't modified
    assert_tensors_identical(d_output_expected_dtensor_copy.to_local(), d_output_expected_dtensor.to_local())

    # Test input gradient
    assert input_x_dtensor.grad.shape == d_input_x_expected_dtensor.shape
    assert input_x_dtensor.grad.stride() == d_input_x_expected_dtensor.stride()
    torch.testing.assert_close(input_x_dtensor.grad.to_local(), d_input_x_expected_dtensor.to_local())

    # Test full tensor gathering - verify distributed results match serial results
    output_global_result_host = output_dtensor_result.full_tensor().cpu()
    d_input_x_global_result_host = input_x_dtensor.grad.full_tensor().cpu()

    # Verify full tensors match expected results
    torch.testing.assert_close(output_global_result_host, output_expected_global_host)
    torch.testing.assert_close(d_input_x_global_result_host, d_input_x_expected_global_host)

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
    ids=lambda x: f"dp={x[0][0]}, cp={x[0][1]}, specify_method={x[1]}, device_type={x[2]}, method_init={x[3]}",
)
@pytest.mark.parametrize(
    "placements", [(Shard(0), Shard(1), Shard(2)), (Shard(0), Shard(1), Replicate())], ids=["shard", "replicate"]
)
@pytest.mark.parametrize(
    "dim_config",
    [
        -1,  # last dimension
        2,  # third dimension
    ],
    ids=lambda x: f"dim={x}",
)
def test_shardwise_log_softmax_parallel(setup_env, placements, dim_config):
    dim = dim_config
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    size_ring = grid_group_sizes["cp"][0]
    B = 2 * grid_group_sizes["dp"]
    N = size_ring * 10  # Number of tokens
    D = 8  # Hidden dimension
    min_init_val = -0.2
    max_init_val = 0.2

    # Skip test if trying to apply log_softmax along a sharded dimension
    input_shape = (B, N, N, D)
    actual_dim = dim if dim >= 0 else len(input_shape) + dim
    for placement in placements:
        if isinstance(placement, Shard) and placement.dim == actual_dim:
            pytest.skip(f"Skipping test: log_softmax along sharded dimension {dim} is not supported")

    seed = 42
    seed_by_rank(0, seed=seed)

    # Create input tensor with proper shape
    input_x_global = torch.empty((B, N, N, D), requires_grad=True, device=device_type)
    init_tensors_uniform([input_x_global], low=min_init_val, high=max_init_val)

    # Run serial forward pass
    input_x_global_host = input_x_global.detach().clone().cpu()
    output_expected_global = serial_shardwise_log_softmax(input_x_global, dim)
    output_expected_global_host = output_expected_global.detach().clone().cpu()

    # Create upstream gradient and run backward pass
    d_output_expected_global = torch.empty_like(output_expected_global)
    init_tensors_uniform([d_output_expected_global], low=min_init_val, high=max_init_val)
    d_output_expected_global_host = d_output_expected_global.detach().clone().cpu()
    output_expected_global.backward(d_output_expected_global)

    # Launch parallel test across all processes
    spawn_multiprocessing(
        parallel_assert_shardwise_log_softmax,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        placements,
        input_x_global_host,
        output_expected_global_host,
        d_output_expected_global_host,
        input_x_global.grad.detach().clone().cpu(),
        dim,
    )


@pytest.mark.parametrize(
    "setup_env",
    [
        ((1, (2, 2)), True, "cuda", "ENV"),
        ((2, (2, 2)), True, "cuda", "ENV"),
        ((1, (3, 3)), True, "cuda", "ENV"),
        ((1, (3, 3)), True, "cpu", "ENV"),
    ],
    indirect=("setup_env",),
    ids=lambda x: f"dp={x[0][0]}, cp={x[0][1]}, specify_method={x[1]}, device_type={x[2]}, method_init={x[3]}",
)
@pytest.mark.parametrize(
    "placements", [(Shard(0), Shard(1), Shard(2)), (Shard(0), Shard(1), Replicate())], ids=["shard", "replicate"]
)
@pytest.mark.parametrize(
    "argmax_config",
    [
        (-1, False),
        (-1, True),
        (2, False),
        (2, True),
    ],
    ids=lambda x: f"dim={x[0]}, keepdim={x[1]}",
)
def test_shardwise_argmax_parallel(setup_env, placements, argmax_config):
    dim, keepdim = argmax_config
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    size_ring = grid_group_sizes["cp"][0]
    B = 2 * grid_group_sizes["dp"]
    N = size_ring * 100
    D = 16

    input_shape = (B, N, N, D)
    actual_dim = dim if dim >= 0 else len(input_shape) + dim
    for placement in placements:
        if isinstance(placement, Shard) and placement.dim == actual_dim:
            pytest.skip(f"Skipping test: argmax along sharded dimension {dim} is not supported")

    seed = 123
    rng = torch.Generator(device=device_type)
    rng.manual_seed(seed)

    input_x_global = torch.rand(input_shape, device=device_type, generator=rng)

    input_x_global_host = input_x_global.detach().clone().cpu()
    output_expected_global = serial_shardwise_argmax(input_x_global, dim, keepdim)
    output_expected_global_host = output_expected_global.detach().clone().cpu()

    spawn_multiprocessing(
        parallel_assert_shardwise_argmax,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        placements,
        input_x_global_host,
        output_expected_global_host,
        dim,
        keepdim,
    )


def serial_shardwise_offset(x: torch.Tensor, dim: int, offset_per_rank: float, num_shards: int) -> torch.Tensor:
    """Serial implementation of shardwise offset operation for comparison.

    Simulates what shardwise_offset does across distributed ranks by applying
    rank-dependent offsets to each shard of the tensor.

    Parameters
    ----------
    x : torch.Tensor
        Input tensor.
    dim : int
        The dimension that would be sharded.
    offset_per_rank : float
        The offset value per rank.
    num_shards : int
        Number of shards (ranks) along the specified dimension.

    Returns
    -------
    torch.Tensor
        Tensor with rank-dependent offsets applied to each shard.
    """
    actual_dim = dim if dim >= 0 else x.ndim + dim
    dim_size = x.shape[actual_dim]
    shard_size = dim_size // num_shards

    output = x.clone()
    for rank in range(num_shards):
        # Create slice for this rank's shard
        slices = [slice(None)] * x.ndim
        slices[actual_dim] = slice(rank * shard_size, (rank + 1) * shard_size)

        # Add offset for this rank
        output[tuple(slices)] = output[tuple(slices)] + rank * offset_per_rank

    return output


def parallel_assert_shardwise_offset(
    rank,
    grid_group_sizes,
    device_type,
    backend,
    env_per_rank,
    placements,
    input_x_global_host,
    output_expected_global_host,
    d_output_expected_global_host,
    d_input_x_expected_global_host,
    dim: int,
    offset_per_rank: float,
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

    # Distribute input tensor
    input_x_dtensor = distribute_tensor(
        input_x_global_host.to(manager.device), device_mesh=manager.device_mesh_subgroups, placements=placements
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

    # Create copy to verify input isn't modified
    input_x_dtensor_copy = input_x_dtensor.detach().clone().requires_grad_(True)

    # Forward pass
    output_dtensor_result = shardwise_offset(input_x_dtensor, dim, offset_per_rank)

    # Verify input wasn't modified
    assert_tensors_identical(
        input_x_dtensor_copy.to_local(), input_x_dtensor.to_local(), check_grad=False, check_grad_fn=False
    )

    # Test forward pass results
    assert output_dtensor_result.shape == output_expected_dtensor.shape
    assert output_dtensor_result.stride() == output_expected_dtensor.stride()
    torch.testing.assert_close(output_dtensor_result.to_local(), output_expected_dtensor.to_local())

    # Backward pass
    d_output_expected_dtensor_copy = d_output_expected_dtensor.detach().clone()
    output_dtensor_result.backward(d_output_expected_dtensor)

    # Verify upstream gradient wasn't modified
    assert_tensors_identical(d_output_expected_dtensor_copy.to_local(), d_output_expected_dtensor.to_local())

    # Test input gradient
    assert input_x_dtensor.grad.shape == d_input_x_expected_dtensor.shape
    assert input_x_dtensor.grad.stride() == d_input_x_expected_dtensor.stride()
    torch.testing.assert_close(input_x_dtensor.grad.to_local(), d_input_x_expected_dtensor.to_local())

    # Test full tensor gathering - verify distributed results match serial results
    output_global_result_host = output_dtensor_result.full_tensor().cpu()
    d_input_x_global_result_host = input_x_dtensor.grad.full_tensor().cpu()

    # Verify full tensors match expected results
    torch.testing.assert_close(output_global_result_host, output_expected_global_host)
    torch.testing.assert_close(d_input_x_global_result_host, d_input_x_expected_global_host)

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
    ids=lambda x: f"dp={x[0][0]}, cp={x[0][1]}, specify_method={x[1]}, device_type={x[2]}, method_init={x[3]}",
)
@pytest.mark.parametrize(
    "placements", [(Shard(0), Shard(1), Shard(2)), (Shard(0), Shard(1), Replicate())], ids=["shard", "replicate"]
)
@pytest.mark.parametrize(
    "offset_config",
    [
        (1, 100.0),  # dim 1 with offset 100.0
        (2, 50.0),  # dim 2 with offset 50.0
    ],
    ids=lambda x: f"dim={x[0]}, offset={x[1]}",
)
def test_shardwise_offset_parallel(setup_env, placements, offset_config):
    dim, offset_per_rank = offset_config
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

    # Skip test if the specified dimension is NOT sharded (offset requires sharded dim)
    input_shape = (B, N, N, D)
    actual_dim = dim if dim >= 0 else len(input_shape) + dim
    dim_is_sharded = False
    mesh_axis_for_dim = None
    for i, placement in enumerate(placements):
        if isinstance(placement, Shard) and placement.dim == actual_dim:
            dim_is_sharded = True
            mesh_axis_for_dim = i
            break

    if not dim_is_sharded:
        pytest.skip(f"Skipping test: dimension {dim} is not sharded, but shardwise_offset requires it to be sharded")

    # Get the number of shards for this dimension from the mesh
    # mesh shape is (dp, cp[0], cp[1]) which corresponds to (mesh_dim_0, mesh_dim_1, mesh_dim_2)
    mesh_shape = (grid_group_sizes["dp"], grid_group_sizes["cp"][0], grid_group_sizes["cp"][1])
    num_shards = mesh_shape[mesh_axis_for_dim]

    seed = 42
    rng = torch.Generator(device=device_type)
    rng.manual_seed(seed)

    # Create input tensor with proper shape
    input_x_global = torch.rand((B, N, N, D), requires_grad=True, device=device_type, generator=rng)

    # Run serial forward pass
    input_x_global_host = input_x_global.detach().clone().cpu()
    output_expected_global = serial_shardwise_offset(input_x_global, dim, offset_per_rank, num_shards)
    output_expected_global_host = output_expected_global.detach().clone().cpu()

    # Create upstream gradient and run backward pass
    # For offset, backward is identity (gradient passes through)
    d_output_expected_global = torch.rand_like(output_expected_global)
    d_output_expected_global_host = d_output_expected_global.detach().clone().cpu()
    output_expected_global.backward(d_output_expected_global)

    # Since offset is x + constant, grad_input = grad_output
    # The serial implementation clones input, so grad flows through clone to input
    # grad_input should equal grad_output
    d_input_x_expected_global_host = d_output_expected_global_host.clone()

    # Launch parallel test across all processes
    spawn_multiprocessing(
        parallel_assert_shardwise_offset,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        placements,
        input_x_global_host,
        output_expected_global_host,
        d_output_expected_global_host,
        d_input_x_expected_global_host,
        dim,
        offset_per_rank,
    )


def assert_offset_error_case(rank, grid_group_sizes, device_type, backend, env_per_rank):
    monkeypatch = pytest.MonkeyPatch()
    if env_per_rank is not None:
        for var_name, value in env_per_rank.items():
            if value == "<INPUT_RANK>":
                monkeypatch.setenv(var_name, f"{rank}")
                continue
            monkeypatch.setenv(var_name, value)

    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    manager = DistributedManager()

    seed_by_rank(0)

    B = 2 * grid_group_sizes["dp"]
    N = grid_group_sizes["cp"][0] * 100  # Number of tokens
    D = 32  # Hidden dimension

    # Test case 1: Offset along non-sharded dimension should raise ValueError
    input_tensor = torch.randn((B, N, D), device=manager.device, requires_grad=True)
    # Shard dim 0 and dim 1, but dim 2 is Replicate
    sharded_dtensor = distribute_tensor(
        input_tensor, device_mesh=manager.device_mesh_subgroups, placements=(Shard(0), Shard(1), Replicate())
    )

    # This should raise an error because dim 2 is not sharded
    with pytest.raises(ValueError, match="Dimension 2 must be sharded for shardwise_offset"):
        shardwise_offset(sharded_dtensor, dim=2, offset_per_rank=100.0)

    # Test case 2: Invalid input type
    with pytest.raises(TypeError, match="Expected DTensor"):
        shardwise_offset(input_tensor, dim=1, offset_per_rank=100.0)  # Regular tensor instead of DTensor

    # Test case 3: Invalid dim type
    with pytest.raises(TypeError, match="Expected int for dim"):
        shardwise_offset(sharded_dtensor, dim=1.5, offset_per_rank=100.0)  # Float instead of int

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env",
    [
        ((1, (2, 2)), True, "cuda", "ENV"),
        ((1, (3, 3)), True, "cpu", "ENV"),
    ],
    indirect=("setup_env",),
    ids=lambda x: f"dp={x[0][0]}, cp={x[0][1]}, specify_method={x[1]}, device_type={x[2]}, method_init={x[3]}",
)
def test_shardwise_offset_error_cases(setup_env):
    """Test error cases for shardwise_offset function."""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    spawn_multiprocessing(
        assert_offset_error_case,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
    )


def serial_shardwise_outer_op(x: torch.Tensor, y: torch.Tensor, axis: int, op: ShardwiseOuterOp) -> torch.Tensor:
    """Serial implementation of shardwise outer op for comparison.

    Takes tensors without singletons and computes the outer operation at the specified axis.
    x: (..., L, ...) at axis
    y: (..., R, ...) at axis
    output: (..., L, R, ...) with one more dimension
    """
    # Unsqueeze to create broadcast-compatible shapes
    x_expanded = x.unsqueeze(axis + 1)  # (..., L, 1, ...)
    y_expanded = y.unsqueeze(axis)  # (..., 1, R, ...)

    if op == ShardwiseOuterOp.SUBTRACT:
        return x_expanded - y_expanded
    elif op == ShardwiseOuterOp.ADD:
        return x_expanded + y_expanded
    elif op == ShardwiseOuterOp.LOGICAL_AND:
        return x_expanded & y_expanded
    elif op == ShardwiseOuterOp.EQUAL:
        return x_expanded == y_expanded
    else:
        raise ValueError(f"Unsupported operation: {op}")


def parallel_assert_shardwise_outer_op(
    rank,
    grid_group_sizes,
    device_type,
    backend,
    env_per_rank,
    placements,
    axis,
    op,
    is_differentiable,
    input_x_global_host,
    input_y_global_host,
    output_expected_global_host,
    output_placements_expected,
    d_output_expected_global_host=None,
    d_input_x_expected_global_host=None,
    d_input_y_expected_global_host=None,
):
    """Parallel assertion function for shardwise_outer_op with all ShardwiseOuterOp types."""
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
    )
    input_y_dtensor = distribute_tensor(
        input_y_global_host.to(manager.device), device_mesh=manager.device_mesh_subgroups, placements=placements
    )

    if is_differentiable:
        input_x_dtensor = input_x_dtensor.requires_grad_(True)
        input_y_dtensor = input_y_dtensor.requires_grad_(True)

    # Distribute expected output
    output_expected_dtensor = distribute_tensor(
        output_expected_global_host.to(manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=output_placements_expected,
        src_data_rank=None,
    )

    # Create copies to verify inputs aren't modified
    # Note: detach().clone() drops requires_grad, so we restore it for proper comparison
    input_x_dtensor_copy = input_x_dtensor.detach().clone().requires_grad_(is_differentiable)
    input_y_dtensor_copy = input_y_dtensor.detach().clone().requires_grad_(is_differentiable)

    # Forward pass
    output_dtensor_result = shardwise_outer_op(input_x_dtensor, input_y_dtensor, axis, op)

    # Verify inputs weren't modified
    assert_tensors_identical(
        input_x_dtensor_copy.to_local(), input_x_dtensor.to_local(), check_grad=False, check_grad_fn=False
    )
    assert_tensors_identical(
        input_y_dtensor_copy.to_local(), input_y_dtensor.to_local(), check_grad=False, check_grad_fn=False
    )

    # Test forward pass results
    assert (
        output_dtensor_result.shape == output_expected_dtensor.shape
    ), f"Shape mismatch: {output_dtensor_result.shape} != {output_expected_dtensor.shape}"
    assert (
        output_dtensor_result.stride() == output_expected_dtensor.stride()
    ), f"Stride mismatch: {output_dtensor_result.stride()} != {output_expected_dtensor.stride()}"
    torch.testing.assert_close(output_dtensor_result.to_local(), output_expected_dtensor.to_local())

    # Test full tensor gathering - verify distributed results match serial results
    output_global_result_host = output_dtensor_result.full_tensor().cpu()
    torch.testing.assert_close(output_global_result_host, output_expected_global_host)

    # Verify placements are correct
    assert output_dtensor_result.placements == output_placements_expected

    if is_differentiable:
        # Distribute expected gradients
        d_output_expected_dtensor = distribute_tensor(
            d_output_expected_global_host.to(manager.device),
            device_mesh=manager.device_mesh_subgroups,
            placements=output_placements_expected,
        )
        d_input_x_expected_dtensor = distribute_tensor(
            d_input_x_expected_global_host.to(manager.device),
            device_mesh=manager.device_mesh_subgroups,
            placements=placements,
            src_data_rank=None,
        )
        d_input_y_expected_dtensor = distribute_tensor(
            d_input_y_expected_global_host.to(manager.device),
            device_mesh=manager.device_mesh_subgroups,
            placements=placements,
            src_data_rank=None,
        )

        # Backward pass
        d_output_expected_dtensor_copy = d_output_expected_dtensor.detach().clone()
        output_dtensor_result.backward(d_output_expected_dtensor)

        # Verify upstream gradient wasn't modified
        assert_tensors_identical(d_output_expected_dtensor_copy.to_local(), d_output_expected_dtensor.to_local())

        # Test input gradients
        assert input_x_dtensor.grad.shape == d_input_x_expected_dtensor.shape
        assert input_x_dtensor.grad.stride() == d_input_x_expected_dtensor.stride()
        torch.testing.assert_close(input_x_dtensor.grad.to_local(), d_input_x_expected_dtensor.to_local())

        assert input_y_dtensor.grad.shape == d_input_y_expected_dtensor.shape
        assert input_y_dtensor.grad.stride() == d_input_y_expected_dtensor.stride()
        torch.testing.assert_close(input_y_dtensor.grad.to_local(), d_input_y_expected_dtensor.to_local())

        # Verify full gradient tensors match expected results
        d_input_x_global_result_host = input_x_dtensor.grad.full_tensor().cpu()
        d_input_y_global_result_host = input_y_dtensor.grad.full_tensor().cpu()
        torch.testing.assert_close(d_input_x_global_result_host, d_input_x_expected_global_host)
        torch.testing.assert_close(d_input_y_global_result_host, d_input_y_expected_global_host)
    else:
        # Ensure no backward possible (non-differentiable)
        with pytest.raises(RuntimeError, match="does not require grad"):
            output_dtensor_result.backward(torch.empty_like(output_dtensor_result))

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env",
    [
        ((1, (2, 2)), True, "cuda", "ENV"),
        ((2, (2, 2)), True, "cuda", "ENV"),
    ],
    indirect=("setup_env",),
    ids=lambda x: f"dp={x[0][0]}, cp={x[0][1]}, specify_method={x[1]}, device_type={x[2]}, method_init={x[3]}",
)
@pytest.mark.parametrize(
    # Note: Cannot use Shard on the axis dimension (outer op axis must be local)
    # Input shape is (B, K, L, D) or (B, K, R, D) with axis=2
    # So axis dim (2) must use Replicate
    "placements",
    [(Shard(0), Shard(1), Replicate()), (Replicate(), Shard(1), Replicate())],
    ids=["shard_batch_and_K", "shard_K_only"],
)
@pytest.mark.parametrize(
    "op",
    [ShardwiseOuterOp.SUBTRACT, ShardwiseOuterOp.ADD, ShardwiseOuterOp.LOGICAL_AND, ShardwiseOuterOp.EQUAL],
    ids=["SUBTRACT", "ADD", "LOGICAL_AND", "EQUAL"],
)
def test_shardwise_outer_op_parallel(setup_env, placements, op):
    """Test shardwise_outer_op for all ShardwiseOuterOp types."""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    size_ring = grid_group_sizes["cp"][0]
    B = 2 * grid_group_sizes["dp"]
    K = size_ring * 4  # Number of windows
    L = 8  # Size at axis for x (e.g., queries)
    R = 16  # Size at axis for y (e.g., keys)
    D = 3  # Feature dimension
    axis = 2  # The axis at which to perform outer operation

    seed = 42
    rng = torch.Generator(device=device_type)
    rng.manual_seed(seed)

    is_differentiable = op in (ShardwiseOuterOp.SUBTRACT, ShardwiseOuterOp.ADD)

    # Create input tensors WITHOUT singletons
    # x: (B, K, L, D) and y: (B, K, R, D)
    if op == ShardwiseOuterOp.SUBTRACT:
        input_x_global = torch.rand((B, K, L, D), requires_grad=True, device=device_type, generator=rng)
        input_y_global = torch.rand((B, K, R, D), requires_grad=True, device=device_type, generator=rng)
    elif op == ShardwiseOuterOp.ADD:
        input_x_global = torch.rand((B, K, L, D), requires_grad=True, device=device_type, generator=rng)
        input_y_global = torch.rand((B, K, R, D), requires_grad=True, device=device_type, generator=rng)
    elif op == ShardwiseOuterOp.LOGICAL_AND:
        input_x_global = torch.randint(0, 2, (B, K, L, D), device=device_type, generator=rng).bool()
        input_y_global = torch.randint(0, 2, (B, K, R, D), device=device_type, generator=rng).bool()
    elif op == ShardwiseOuterOp.EQUAL:
        num_unique_values = 5
        input_x_global = torch.randint(0, num_unique_values, (B, K, L, D), device=device_type, generator=rng)
        input_y_global = torch.randint(0, num_unique_values, (B, K, R, D), device=device_type, generator=rng)
    else:
        raise ValueError(f"Unknown op: {op}")

    # Run serial forward pass
    input_x_global_host = input_x_global.detach().clone().cpu()
    input_y_global_host = input_y_global.detach().clone().cpu()
    output_expected_global = serial_shardwise_outer_op(input_x_global, input_y_global, axis, op)
    output_expected_global_host = output_expected_global.detach().clone().cpu()

    # Compute expected output placements (axis+1 is inserted, so Shard dims > axis shift)
    output_placements_expected = list(placements)
    for i, p in enumerate(placements):
        if isinstance(p, Shard) and p.dim > axis:
            output_placements_expected[i] = Shard(p.dim + 1)
    output_placements_expected = tuple(output_placements_expected)

    # Prepare gradient data for differentiable ops
    d_output_expected_global_host = None
    d_input_x_expected_global_host = None
    d_input_y_expected_global_host = None

    if is_differentiable:
        d_output_expected_global = torch.rand_like(output_expected_global)
        d_output_expected_global_host = d_output_expected_global.detach().clone().cpu()
        output_expected_global.backward(d_output_expected_global)
        d_input_x_expected_global_host = input_x_global.grad.detach().clone().cpu()
        d_input_y_expected_global_host = input_y_global.grad.detach().clone().cpu()

    # Launch parallel test across all processes
    spawn_multiprocessing(
        parallel_assert_shardwise_outer_op,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        placements,
        axis,
        op,
        is_differentiable,
        input_x_global_host,
        input_y_global_host,
        output_expected_global_host,
        output_placements_expected,
        d_output_expected_global_host,
        d_input_x_expected_global_host,
        d_input_y_expected_global_host,
    )


def assert_outer_op_error_cases(rank, grid_group_sizes, device_type, backend, env_per_rank):
    """Test error cases for shardwise_outer_op."""
    monkeypatch = pytest.MonkeyPatch()
    if env_per_rank is not None:
        for var_name, value in env_per_rank.items():
            if value == "<INPUT_RANK>":
                monkeypatch.setenv(var_name, f"{rank}")
                continue
            monkeypatch.setenv(var_name, value)

    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    manager = DistributedManager()

    seed_by_rank(0)

    B = 2 * grid_group_sizes["dp"]
    K = grid_group_sizes["cp"][0] * 4
    L = 8
    R = 16
    D = 3
    axis = 2

    placements = (Shard(0), Shard(1), Replicate())

    # Test case 1: Invalid input type (regular tensor instead of DTensor)
    regular_tensor = torch.randn((B, K, L, D), device=manager.device)
    dtensor = distribute_tensor(
        torch.randn((B, K, R, D), device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements,
    )

    with pytest.raises(TypeError, match="Expected DTensor for lhs"):
        shardwise_outer_op(regular_tensor, dtensor, axis, ShardwiseOuterOp.SUBTRACT)

    with pytest.raises(TypeError, match="Expected DTensor for rhs"):
        shardwise_outer_op(dtensor, regular_tensor, axis, ShardwiseOuterOp.SUBTRACT)

    # Test case 2: Mismatched placements
    placements2 = (Shard(0), Replicate(), Shard(2))
    dtensor1 = distribute_tensor(
        torch.randn((B, K, L, D), device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements,
    )
    dtensor2 = distribute_tensor(
        torch.randn((B, K, R, D), device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements2,
    )

    with pytest.raises(ValueError, match="must have the same placements"):
        shardwise_outer_op(dtensor1, dtensor2, axis, ShardwiseOuterOp.SUBTRACT)

    # Test case 3: Trying to shard the axis dimension (outer op must be local)
    placements_axis_shard = (Shard(0), Shard(1), Shard(2))  # Try to shard axis dim (2)
    dtensor3 = distribute_tensor(
        torch.randn((B, K, L, D), device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements_axis_shard,
    )
    dtensor4 = distribute_tensor(
        torch.randn((B, K, R, D), device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements_axis_shard,
    )

    with pytest.raises(ValueError, match="Cannot shard dimension.*outer operation axis"):
        shardwise_outer_op(dtensor3, dtensor4, axis, ShardwiseOuterOp.SUBTRACT)

    # Test case 4: Invalid axis type
    dtensor5 = distribute_tensor(
        torch.randn((B, K, L, D), device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements,
    )
    dtensor6 = distribute_tensor(
        torch.randn((B, K, R, D), device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements,
    )

    with pytest.raises(TypeError, match="Expected int for axis"):
        shardwise_outer_op(dtensor5, dtensor6, "invalid", ShardwiseOuterOp.SUBTRACT)

    # Test case 5: axis out of bounds
    with pytest.raises(ValueError, match="axis.*out of bounds"):
        shardwise_outer_op(dtensor5, dtensor6, 10, ShardwiseOuterOp.SUBTRACT)

    # Test case 6: Mismatched number of dimensions
    dtensor7 = distribute_tensor(
        torch.randn((B, K, L), device=manager.device),  # 3D instead of 4D
        device_mesh=manager.device_mesh_subgroups,
        placements=placements,
    )

    with pytest.raises(ValueError, match="must have the same number of dimensions"):
        shardwise_outer_op(dtensor5, dtensor7, axis, ShardwiseOuterOp.SUBTRACT)

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env",
    [
        ((1, (2, 2)), True, "cuda", "ENV"),
    ],
    indirect=("setup_env",),
    ids=lambda x: f"dp={x[0][0]}, cp={x[0][1]}, specify_method={x[1]}, device_type={x[2]}, method_init={x[3]}",
)
def test_shardwise_outer_op_error_cases(setup_env):
    """Test error cases for shardwise_outer_op."""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    spawn_multiprocessing(
        assert_outer_op_error_cases,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
    )


if __name__ == "__main__":
    unittest.main()
