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


from typing import Dict, Optional

import pytest
import torch
from torch.distributed.tensor import DeviceMesh, DTensor, Placement, Replicate, Shard, distribute_tensor

from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.layers.repeat_interleave import shardwise_repeat_interleave
from boltz.testing.utils import assert_tensors_identical, seed_by_rank, spawn_multiprocessing


def compute_global_expectation(shape, repeats, dim_to_repeat, device):
    """Compute global expectation using standard PyTorch operations."""
    # Create input tensor
    x = torch.rand(*shape, device=device, requires_grad=True)

    # Compute on global tensor using standard repeat_interleave operation
    y = torch.repeat_interleave(x, repeats=repeats, dim=dim_to_repeat)

    # Create gradients for backward pass
    dy = torch.rand_like(y)

    # Backward pass on global tensor
    y.backward(dy)

    # Collect input gradient
    input_grad = x.grad.detach().clone()

    return x.detach().clone(), y.detach().clone(), input_grad, dy.detach().clone()


def compute_dtensor_native(
    input_global: torch.Tensor,
    dy_global: torch.Tensor,
    device_mesh: DeviceMesh,
    input_placements: tuple[Placement, ...],
    repeats: int,
    dim_to_repeat: int,
) -> tuple[DTensor, DTensor]:
    """Compute DTensor native operations for comparison."""
    # Create DTensor native input
    input_dtensor = distribute_tensor(input_global.detach().clone(), device_mesh, input_placements).requires_grad_(True)

    # Forward pass with native DTensor repeat_interleave operation
    y_dtensor_result = torch.repeat_interleave(input_dtensor, repeats=repeats, dim=dim_to_repeat)

    # Backward pass with native DTensor op
    dy_dtensor = distribute_tensor(dy_global.detach().clone(), device_mesh, y_dtensor_result.placements)
    y_dtensor_result.backward(dy_dtensor)

    input_grad_dtensor = input_dtensor.grad

    return input_grad_dtensor, y_dtensor_result


def compute_shardwise_repeat_interleave_with_validation(
    input_global: torch.Tensor,
    dy_global: torch.Tensor,
    device_mesh: DeviceMesh,
    input_placements: tuple[Placement, ...],
    repeats: int,
    dim_to_repeat: int,
    label_test_case: str,
) -> tuple[DTensor, DTensor, DTensor]:
    """
    Compute shardwise_repeat_interleave forward and backward pass with input validation checks.

    Returns:
        y_dtensor_result: Forward pass result
        input_dtensor: Input tensor with computed gradient
        dy_dtensor: Distributed upstream gradient
    """
    # Create DTensor input
    input_dtensor = distribute_tensor(input_global.detach().clone(), device_mesh, input_placements).requires_grad_(True)
    input_dtensor_copy = input_dtensor.detach().clone().requires_grad_(True)

    # Compute on distributed tensor using shardwise_repeat_interleave
    y_dtensor_result = shardwise_repeat_interleave(input_dtensor, repeats, dim_to_repeat)

    # verify no change to the fwd input
    assert_tensors_identical(
        input_dtensor.to_local(), input_dtensor_copy.to_local(), check_grad=False, check_grad_fn=False
    )

    # Distribute the upstream adjoint for backward pass
    dy_dtensor = distribute_tensor(dy_global.detach().clone(), device_mesh, y_dtensor_result.placements)

    # Perform backward pass
    dy_dtensor_copy = dy_dtensor.detach().clone()
    y_dtensor_result.backward(dy_dtensor)

    # verify no change to the bwd input
    assert_tensors_identical(dy_dtensor.to_local(), dy_dtensor_copy.to_local(), check_grad=False, check_grad_fn=False)

    # verify input gradient placements are consistent with input placements
    assert (
        input_dtensor.grad.placements == input_placements
    ), f"{label_test_case} inconsistent input gradient placements with input placements"

    return y_dtensor_result, input_dtensor, dy_dtensor


def parallel_assert_dtensor_repeat_interleave(
    rank: int,
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

    # each rank uses the same seed to generate the same input tensors
    seed_by_rank(0, seed=42)

    size_ring = len(manager.subgroups_ranks["cp"][0])

    # Set test parameters
    shape = (3, grid_group_sizes["dp"] * 2, 5, size_ring * 4, 5)
    repeats = 3  # Number of times to repeat each element
    # Shard dimensions for input tensor
    # this emulates the sharded single representation in the Boltz model
    input_placements = (Shard(dim=1), Shard(dim=3), Replicate())

    # Test all dimensions (all are valid since sharded dimensions are evenly divided)
    valid_dims_to_repeat = [0, 1, 2, 3, 4, -1, -2, -3, -4, -5]

    # Test valid repeat_interleave dimensions
    for dim_to_repeat in valid_dims_to_repeat:
        label_test_case = f"for dim {dim_to_repeat}\n"

        # Compute global expectations
        input_global, y_expected_global, input_grad_expected_global, dy_global = compute_global_expectation(
            shape, repeats, dim_to_repeat, manager.device
        )

        # use DTensor native op as an alternative reference
        input_grad_dtensor_native, y_dtensor_result_native = compute_dtensor_native(
            input_global, dy_global, manager.device_mesh_subgroups, input_placements, repeats, dim_to_repeat
        )

        # Compute shardwise_repeat_interleave forward and backward with validation
        y_dtensor_result, input_dtensor, dy_dtensor = compute_shardwise_repeat_interleave_with_validation(
            input_global,
            dy_global,
            manager.device_mesh_subgroups,
            input_placements,
            repeats,
            dim_to_repeat,
            label_test_case,
        )

        # ===================================================================
        # BLOCK 1: Check against DTensor native reference
        # ===================================================================

        # check metadata against DTensor native
        assert (
            y_dtensor_result.placements == y_dtensor_result_native.placements
        ), f"{label_test_case} placements mismatch"
        assert y_dtensor_result.shape == y_dtensor_result_native.shape, f"{label_test_case} shape mismatch"
        assert y_dtensor_result.stride() == y_dtensor_result_native.stride(), f"{label_test_case} stride mismatch"

        # compare forward result with native DTensor op
        torch.testing.assert_close(
            y_dtensor_result.to_local(),
            y_dtensor_result_native.to_local(),
            atol=0,
            rtol=0,
            msg=lambda m: f"{label_test_case} {m}",
        )

        # compare global tensors between shardwise_repeat_interleave and native DTensor results
        y_result_global = y_dtensor_result.full_tensor()
        y_result_global_native = y_dtensor_result_native.full_tensor()

        torch.testing.assert_close(
            y_result_global,
            y_result_global_native,
            atol=0,
            rtol=0,
            msg=lambda m: f"{label_test_case} output vs native: {m}",
        )

        # assert input gradient metadata and values against DTensor native
        assert (
            input_dtensor.grad.placements == input_grad_dtensor_native.placements
        ), f"{label_test_case} input gradient placements mismatch"
        assert (
            input_dtensor.grad.shape == input_grad_dtensor_native.shape
        ), f"{label_test_case} input gradient shape mismatch"
        assert (
            input_dtensor.grad.stride() == input_grad_dtensor_native.stride()
        ), f"{label_test_case} input gradient stride mismatch"

        torch.testing.assert_close(
            input_dtensor.grad.to_local(),
            input_grad_dtensor_native.to_local(),
            atol=0,
            rtol=0,
            msg=lambda m: f"{label_test_case} input gradient mismatch: {m}",
        )

        torch.testing.assert_close(
            input_dtensor.grad.full_tensor(),
            input_grad_dtensor_native.full_tensor(),
            atol=0,
            rtol=0,
            msg=lambda m: f"{label_test_case} input gradient mismatch: {m}",
        )

        # ===================================================================
        # BLOCK 2: Check against global serial expectation
        # ===================================================================
        y_dtensor_expected = distribute_tensor(
            y_expected_global, manager.device_mesh_subgroups, y_dtensor_result.placements
        )

        # Compare results with expected local shards
        torch.testing.assert_close(
            y_dtensor_result.to_local(),
            y_dtensor_expected.to_local(),
            atol=0,
            rtol=0,
            msg=lambda m: f"{label_test_case} {m}",
        )

        # compare forward result with global expectation
        torch.testing.assert_close(
            y_result_global,
            y_expected_global,
            atol=0,
            rtol=0,
            msg=lambda m: f"{label_test_case} output vs global expectation: {m}",
        )

        # create distributed tensor from global result for local shard comparison
        input_grad_expected_dtensor = distribute_tensor(
            input_grad_expected_global, manager.device_mesh_subgroups, input_placements
        )

        # compare local shards with expected
        torch.testing.assert_close(
            input_dtensor.grad.to_local(),
            input_grad_expected_dtensor.to_local(),
            atol=0,
            rtol=0,
            msg=lambda m: f"{label_test_case} input gradient vs global expectation: {m}",
        )

        torch.testing.assert_close(
            input_dtensor.grad.full_tensor(),
            input_grad_expected_global,
            atol=0,
            rtol=0,
            msg=lambda m: f"{label_test_case} input gradient vs global expectation: {m}",
        )

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
def test_dtensor_repeat_interleave(setup_env):
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    spawn_multiprocessing(
        parallel_assert_dtensor_repeat_interleave,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
    )
