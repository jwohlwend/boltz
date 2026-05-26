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
from math import isqrt
from typing import Dict, Optional

import pytest
import torch
from torch.distributed.tensor import Replicate, Shard, distribute_module, distribute_tensor

from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.layers.layernorm import LayerNormParamsReplicated
from boltz.testing.utils import assert_all_identical, seed_by_rank


def compute_global_expectation(batch_size, seq_len, features: int | list[int], device):
    """Compute expected results using global tensors.

    Args:
        batch_size: Batch size
        seq_len: Sequence length
        features: Feature dimensions
        device: Device to place tensors on

    Returns:
        tuple: (x_global, layernorm_state_dict, y_global_expected,
                x_global_grad, weight_global_grad, bias_global_grad, dy_global)
    """
    # Create global tensors with deterministic values for reproducibility
    if isinstance(features, int):
        features = [features]
    x_global = torch.randn(batch_size, seq_len, *features, device=device, requires_grad=True)
    layernorm = torch.nn.LayerNorm(features, device=device)
    state_dict_layernorm = layernorm.state_dict()

    # Clone inputs for distribution
    x_global_clone = x_global.detach().clone()

    # Compute on global tensors using standard layernorm operation
    y_global_expected = layernorm(x_global)

    # Create gradients for backward pass
    dy_global = torch.rand_like(y_global_expected)

    # Backward pass on global tensors
    y_global_expected.backward(dy_global)

    return (
        x_global_clone,
        state_dict_layernorm,
        y_global_expected.detach().clone(),
        x_global.grad.detach().clone(),
        layernorm.weight.grad.detach().clone(),
        layernorm.bias.grad.detach().clone(),
        dy_global.detach().clone(),
    )


def parallel_assert_dtensor_layernorm(
    rank: int,
    batch_size: int,
    grid_group_sizes: Dict[str, int],
    device_type: str,
    backend: str,
    env_map: Optional[Dict[str, str]] = None,
):
    """Test distributed layernorm operation in a parallel environment.

    This test validates that the LayerNormParamsReplicated produces identical results to
    standard nn.functional.layer_norm operations with global tensors. It verifies:

    1. Forward pass produces the same results as global tensor computation
    2. Backward pass correctly propagates gradients through the distributed operation
    3. Results and gradients match the equivalent global tensor operations

    Args:
        rank: The process rank in the distributed environment
        grid_group_sizes: Dictionary mapping group names to their sizes for distributed setup
        device_type: Device to run the test on ("cpu" or "cuda")
        backend: The distributed backend to use (e.g., "gloo", "nccl")
        env_map: Optional dictionary of environment variables to set before initialization
    """
    monkeypatch = pytest.MonkeyPatch()
    if env_map is not None:
        for var_name, value in env_map.items():
            if value == "<INPUT_RANK>":
                monkeypatch.setenv(var_name, f"{rank}")
                continue
            monkeypatch.setenv(var_name, value)

    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    manager = DistributedManager()

    size_cp = len(manager.group_ranks["cp"])
    size_ring = isqrt(size_cp)
    if size_ring * size_ring != size_cp:
        raise ValueError(f"cp group size {size_cp} is not a square int")

    # Set test parameters
    seq_len_per_rank = 4
    seq_len_global = size_ring * seq_len_per_rank
    features = [3, 2]

    # Set random seed based on rank for reproducibility
    seed_by_rank(0)

    # Compute global expectations
    (
        x_global,
        state_dict_layernorm,
        y_global_expected,
        x_global_grad,
        weight_global_grad,
        bias_global_grad,
        dy_global,
    ) = compute_global_expectation(batch_size, seq_len_global, features, manager.device)

    # Create distributed tensors
    # Shard the sequence dimension (dim=1) for input tensor
    # this emulates the sharded single representation in the Boltz model
    input_placements = [Shard(dim=0), Shard(dim=1), Replicate()]
    x_dtensor = distribute_tensor(x_global, manager.device_mesh_subgroups, input_placements)
    x_dtensor.requires_grad = True

    layernorm_local = torch.nn.LayerNorm(features, device=manager.device)
    layernorm_local.load_state_dict(state_dict_layernorm)
    layer = LayerNormParamsReplicated(layernorm_local, manager.device_mesh_subgroups)

    layernorm_local_copy = torch.nn.LayerNorm(features, device=manager.device)
    layernorm_local_copy.load_state_dict(state_dict_layernorm)
    layer_dtensor_native = distribute_module(layernorm_local_copy, manager.device_mesh_subgroups)
    x_dtensor_native = x_dtensor.detach().clone().requires_grad_(True)

    # Compute on distributed tensors using LayerNormParamsReplicated
    y_dtensor_result = layer(x_dtensor)
    y_dtensor_native = layer_dtensor_native(x_dtensor_native)

    # Distribute the upstream adjoint for backward pass
    dy_dtensor = distribute_tensor(dy_global, manager.device_mesh_subgroups, input_placements)
    dy_dtensor_native = distribute_tensor(dy_global.detach().clone(), manager.device_mesh_subgroups, input_placements)

    # Perform backward pass
    y_dtensor_result.backward(dy_dtensor)
    y_dtensor_native.backward(dy_dtensor_native)

    # Create distributed tensors from global gradients for comparison
    x_grad_dtensor_expected = distribute_tensor(x_global_grad, manager.device_mesh_subgroups, input_placements)
    y_dtensor_expected = distribute_tensor(y_global_expected, manager.device_mesh_subgroups, input_placements)
    weight_grad_dtensor_expected = distribute_tensor(
        weight_global_grad, manager.device_mesh_subgroups, layer.weight.placements
    )
    bias_grad_dtensor_expected = distribute_tensor(
        bias_global_grad, manager.device_mesh_subgroups, layer.bias.placements
    )

    # Compare results with expected local shards
    torch.testing.assert_close(y_dtensor_expected, y_dtensor_result)
    torch.testing.assert_close(x_grad_dtensor_expected, x_dtensor.grad)
    torch.testing.assert_close(weight_grad_dtensor_expected, layer.weight.grad)
    torch.testing.assert_close(bias_grad_dtensor_expected, layer.bias.grad)

    # Compare results with native DTensor implementation
    assert y_dtensor_result.shape == y_dtensor_native.shape
    assert y_dtensor_result.stride() == y_dtensor_native.stride()
    assert x_dtensor.grad.shape == x_grad_dtensor_expected.shape
    assert x_dtensor.grad.stride() == x_grad_dtensor_expected.stride()
    assert layer.weight.grad.shape == weight_grad_dtensor_expected.shape
    assert layer.weight.grad.stride() == weight_grad_dtensor_expected.stride()
    assert layer.bias.grad.shape == bias_grad_dtensor_expected.shape
    assert layer.bias.grad.stride() == bias_grad_dtensor_expected.stride()

    torch.testing.assert_close(y_dtensor_native, y_dtensor_result)
    torch.testing.assert_close(x_dtensor_native.grad, x_dtensor.grad)
    torch.testing.assert_close(layer_dtensor_native.weight.grad, layer.weight.grad)
    torch.testing.assert_close(layer_dtensor_native.bias.grad, layer.bias.grad)

    # Collect results as global tensors and compare with original global tensors
    y_global_result = y_dtensor_result.full_tensor()
    x_grad_global_result = x_dtensor.grad.full_tensor()
    weight_grad_global_result = layer.weight.grad.full_tensor()
    bias_grad_global_result = layer.bias.grad.full_tensor()

    # Assert output and input gradients match the global computation
    torch.testing.assert_close(y_global_result, y_global_expected)
    torch.testing.assert_close(x_grad_global_result, x_global_grad)
    torch.testing.assert_close(weight_grad_global_result, weight_global_grad)
    torch.testing.assert_close(bias_grad_global_result, bias_global_grad)

    # assert the parameter gradients are identical across all ranks
    assert_all_identical(weight_grad_global_result, manager.group["cp"])
    assert_all_identical(bias_grad_global_result, manager.group["cp"])

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
def test_dtensor_layernorm(setup_env: tuple[dict, int, str, str, str, dict[str, str]]):
    """Test distributed layernorm operation across multiple processes.

    This parametrized test launches multiple processes to test the LayerNormParamsReplicated
    with different configurations. It verifies that layernorm operations work correctly
    in a distributed setting across various:

    - Data parallel (dp) and compute parallel (cp) group sizes
    - Device types (CPU/CUDA)
    - Initialization methods

    The test ensures operations on distributed tensors produce results identical
    to equivalent operations on global tensors, validating the correctness of the
    LayerNormParamsReplicated implementation.

    Args:
        setup_env: Fixture providing the distributed environment configuration
    """
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    batch_size = 2 * grid_group_sizes["dp"]

    torch.multiprocessing.set_start_method("spawn", force=True)
    torch.multiprocessing.spawn(
        fn=parallel_assert_dtensor_layernorm,
        args=(
            batch_size,
            grid_group_sizes,
            device_type,
            backend,
            env_per_rank,
        ),
        nprocs=world_size,
        join=True,
    )


def parallel_assert_dtensor_layernorm_raise_uneven_sharding(
    rank: int,
    batch_size: int,
    grid_group_sizes: Dict[str, int],
    device_type: str,
    backend: str,
    env_map: Optional[Dict[str, str]] = None,
):
    monkeypatch = pytest.MonkeyPatch()
    if env_map is not None:
        for var_name, value in env_map.items():
            if value == "<INPUT_RANK>":
                monkeypatch.setenv(var_name, f"{rank}")
                continue
            monkeypatch.setenv(var_name, value)

    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    manager = DistributedManager()

    size_cp = len(manager.group_ranks["cp"])
    size_ring = isqrt(size_cp)
    if size_ring * size_ring != size_cp:
        raise ValueError(f"cp group size {size_cp} is not a square int")

    # Set test parameters
    seq_len_per_rank = 4
    seq_len_global = size_ring * seq_len_per_rank
    features = [3, 2]

    x_global = torch.ones(batch_size, seq_len_global, *features, device=manager.device)

    # Create distributed tensors
    # Shard the sequence dimension (dim=1) for input tensor
    # this emulates the sharded single representation in the Boltz model
    input_placements = [Shard(dim=0), Shard(dim=1), Replicate()]
    x_dtensor = distribute_tensor(x_global, manager.device_mesh_subgroups, input_placements)
    x_dtensor.requires_grad = True

    layernorm_local = torch.nn.LayerNorm(features, device=manager.device)
    layer = LayerNormParamsReplicated(layernorm_local, manager.device_mesh_subgroups)

    # should raise here
    with pytest.raises(
        ValueError,
        match="Uneven sharding tensor dimension 0 of size 3 along device mesh dimension 0 of size 2 is not supported",
    ):
        _ = layer(x_dtensor)

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env",
    itertools.product([(2, (1, 1))], [True], ["cpu"], ["ENV"]),
    indirect=("setup_env",),
    ids=lambda x: f"dp={x[0][0]}, cp={x[0][1]}, specify_method={x[1]}, device_type={x[2]}, method_init={x[3]}",
)
def test_dtensor_layernorm_raise_on_uneven_sharding(setup_env: tuple[dict, int, str, str, str, dict[str, str]]):
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env
    batch_size = 3

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    torch.multiprocessing.set_start_method("spawn", force=True)
    torch.multiprocessing.spawn(
        fn=parallel_assert_dtensor_layernorm_raise_uneven_sharding,
        args=(
            batch_size,
            grid_group_sizes,
            device_type,
            backend,
            env_per_rank,
        ),
        nprocs=world_size,
        join=True,
    )
