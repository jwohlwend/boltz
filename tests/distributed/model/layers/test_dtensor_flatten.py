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


from math import isqrt
from typing import Dict, Optional

import pytest
import torch
from torch.distributed.tensor import DeviceMesh, DTensor, Placement, Replicate, Shard, distribute_tensor

from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.layers.flatten_and_unflatten import shardwise_flatten, shardwise_flatten_sharded
from boltz.testing.utils import assert_tensors_identical, seed_by_rank, spawn_multiprocessing


def compute_global_expectation(shape, start_dim, end_dim, device):
    """Compute global expectation using standard PyTorch operations."""
    # Create tensor for flattening
    x = torch.rand(*shape, device=device, requires_grad=True)

    # Compute on global tensor using standard flatten operation
    y = torch.flatten(x, start_dim=start_dim, end_dim=end_dim)

    # Create gradients for backward pass
    dy = torch.rand_like(y)

    # Backward pass on global tensor
    y.backward(dy)

    # Collect input gradient
    input_grad = x.grad.detach().clone()

    return x.detach().clone(), y.detach().clone(), input_grad, dy.detach().clone()


def compute_dtensor_native(
    x_global: torch.Tensor,
    dy_global: torch.Tensor,
    device_mesh: DeviceMesh,
    input_placements: tuple[Placement, ...],
    start_dim: int,
    end_dim: int,
) -> tuple[DTensor, DTensor]:
    """Compute DTensor native operations for comparison."""
    # Create DTensor native input
    x_dtensor = distribute_tensor(x_global.detach().clone(), device_mesh, input_placements).requires_grad_(True)

    # Forward pass with native DTensor flatten operation
    y_dtensor_result = torch.flatten(x_dtensor, start_dim=start_dim, end_dim=end_dim)

    # Backward pass with native DTensor op
    dy_dtensor = distribute_tensor(dy_global.detach().clone(), device_mesh, y_dtensor_result.placements)
    y_dtensor_result.backward(dy_dtensor)

    x_grad_dtensor = x_dtensor.grad

    return x_grad_dtensor, y_dtensor_result


def compute_shardwise_flatten_with_validation(
    x_global: torch.Tensor,
    dy_global: torch.Tensor,
    device_mesh: DeviceMesh,
    input_placements: tuple[Placement, ...],
    start_dim: int,
    end_dim: int,
    label_test_case: str,
) -> tuple[DTensor, DTensor, DTensor]:
    """
    Compute shardwise_flatten forward and backward pass with input validation checks.

    Returns:
        y_dtensor_result: Forward pass result
        x_dtensor: Input tensor with computed gradient
        dy_dtensor: Distributed upstream gradient
    """
    # Create DTensor input
    x_dtensor = distribute_tensor(x_global.detach().clone(), device_mesh, input_placements).requires_grad_(True)
    x_dtensor_copy = x_dtensor.detach().clone().requires_grad_(True)

    # Compute on distributed tensor using shardwise_flatten
    y_dtensor_result = shardwise_flatten(x_dtensor, start_dim=start_dim, end_dim=end_dim)

    # verify no change to the fwd input
    assert_tensors_identical(x_dtensor.to_local(), x_dtensor_copy.to_local(), check_grad=False, check_grad_fn=False)

    # Distribute the upstream adjoint for backward pass
    dy_dtensor = distribute_tensor(dy_global.detach().clone(), device_mesh, y_dtensor_result.placements)

    # Perform backward pass
    dy_dtensor_copy = dy_dtensor.detach().clone()
    y_dtensor_result.backward(dy_dtensor)

    # verify no change to the bwd input
    assert_tensors_identical(dy_dtensor.to_local(), dy_dtensor_copy.to_local(), check_grad=False, check_grad_fn=False)

    # verify input gradient placements are consistent with input placements
    assert (
        x_dtensor.grad.placements == input_placements
    ), f"{label_test_case} inconsistent input gradient placements with input placements"

    return y_dtensor_result, x_dtensor, dy_dtensor


def parallel_assert_dtensor_flatten(
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

    size_cp = len(manager.group_ranks["cp"])
    size_ring = isqrt(size_cp)
    if size_ring * size_ring != size_cp:
        raise ValueError(f"cp group size {size_cp} is not a square int")

    # Set test parameters - 8D tensor for comprehensive testing
    shape = (3, 5, grid_group_sizes["dp"] * 2, 5, size_ring * 4, 5, 3, 2)
    # Shard the sequence dimension (dim=2) and another dimension (dim=4) for input tensor
    # this emulates the sharded single representation in the Boltz model
    input_placements = (Shard(dim=2), Shard(dim=4), Replicate())

    # Test valid flattening dimensions (not sharded)
    # Sharded dims are 2 and 4, so valid ranges must not include these dimensions
    valid_flatten_params = [
        (0, 1),  # flatten dims 0,1 (no sharded dims)
        (1, 1),  # flatten just dim 1 (no sharded dims)
        (3, 3),  # flatten just dim 3 (no sharded dims)
        (5, 7),  # flatten dims 5,6,7 (no sharded dims)
        (-2, -1),  # flatten dims 6,7 (negative indexing)
        (-3, -1),  # flatten dims 5,6,7 (negative indexing)
        (-8, -7),  # flatten dims 0,1
    ]

    # Test invalid flattening dimensions (include sharded dims 2 and/or 4)
    invalid_flatten_params = [
        (0, 2),  # flatten dims 0,1,2 (includes sharded dim=2)
        (0, 3),  # flatten dims 0,1,2,3 (includes sharded dim=2)
        (2, 3),  # flatten dims 2,3 (includes sharded dim=2)
        (1, 4),  # flatten dims 1,2,3,4 (includes both sharded dims 2,4)
        (3, 5),  # flatten dims 3,4,5 (includes sharded dim=4)
        (4, 5),  # flatten dims 4,5 (includes sharded dim=4)
        (0, 4),  # flatten dims 0,1,2,3,4 (includes both sharded dims)
        (-6, -4),  # flatten dims 2,3,4 (includes both sharded dims 2,4)
        (-4, -1),  # flatten dims 4,5,6,7 (includes sharded dim=4)
    ]

    # Test valid flattening dimensions
    for start_dim, end_dim in valid_flatten_params:
        label_test_case = f"for start_dim={start_dim}, end_dim={end_dim}\n"

        # Compute global expectations
        x_global, y_expected_global, x_grad_expected_global, dy_global = compute_global_expectation(
            shape, start_dim, end_dim, manager.device
        )

        # use DTensor native op as an alternative reference
        x_grad_dtensor_native, y_dtensor_result_native = compute_dtensor_native(
            x_global, dy_global, manager.device_mesh_subgroups, input_placements, start_dim, end_dim
        )

        # Compute shardwise_flatten forward and backward with validation
        y_dtensor_result, x_dtensor, dy_dtensor = compute_shardwise_flatten_with_validation(
            x_global, dy_global, manager.device_mesh_subgroups, input_placements, start_dim, end_dim, label_test_case
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

        # compare global tensors between shardwise_flatten and native DTensor results
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
            x_dtensor.grad.placements == x_grad_dtensor_native.placements
        ), f"{label_test_case} input gradient placements mismatch"
        assert x_dtensor.grad.shape == x_grad_dtensor_native.shape, f"{label_test_case} input gradient shape mismatch"
        assert (
            x_dtensor.grad.stride() == x_grad_dtensor_native.stride()
        ), f"{label_test_case} input gradient stride mismatch"

        torch.testing.assert_close(
            x_dtensor.grad.to_local(),
            x_grad_dtensor_native.to_local(),
            atol=0,
            rtol=0,
            msg=lambda m: f"{label_test_case} input gradient mismatch: {m}",
        )

        torch.testing.assert_close(
            x_dtensor.grad.full_tensor(),
            x_grad_dtensor_native.full_tensor(),
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
        x_grad_expected_dtensor = distribute_tensor(
            x_grad_expected_global, manager.device_mesh_subgroups, input_placements
        )

        # compare local shard with expected
        torch.testing.assert_close(
            x_dtensor.grad.to_local(),
            x_grad_expected_dtensor.to_local(),
            atol=0,
            rtol=0,
            msg=lambda m: f"{label_test_case} input gradient vs global expectation: {m}",
        )

        torch.testing.assert_close(
            x_dtensor.grad.full_tensor(),
            x_grad_expected_global,
            atol=0,
            rtol=0,
            msg=lambda m: f"{label_test_case} input gradient vs global expectation: {m}",
        )

    # Test invalid flattening dimensions (should raise NotImplementedError)
    for start_dim, end_dim in invalid_flatten_params:
        label_test_case = f"for invalid start_dim={start_dim}, end_dim={end_dim}\n"

        # Compute global expectations (this should work fine)
        x_global, _, _, _ = compute_global_expectation(shape, start_dim, end_dim, manager.device)

        # Create DTensor input
        x_dtensor = distribute_tensor(x_global, manager.device_mesh_subgroups, input_placements)
        x_dtensor.requires_grad = True

        # This should raise due to sharded dimension in flatten range
        with pytest.raises(NotImplementedError, match="Flattening dimension .* sharded by device_mesh axis"):
            shardwise_flatten(x_dtensor, start_dim=start_dim, end_dim=end_dim)

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
def test_dtensor_flatten(setup_env):
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    spawn_multiprocessing(
        parallel_assert_dtensor_flatten,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
    )


def compute_shardwise_flatten_sharded_with_validation(
    x_global: torch.Tensor,
    dy_global: torch.Tensor,
    device_mesh: DeviceMesh,
    input_placements: tuple[Placement, ...],
    start_dim: int,
    end_dim: int,
    label_test_case: str,
) -> tuple[DTensor, DTensor, DTensor]:
    """
    Compute shardwise_flatten_sharded forward and backward pass with input validation checks.

    This function is for testing flatten operations that involve the sharded dimension.

    Returns:
        y_dtensor_result: Forward pass result
        x_dtensor: Input tensor with computed gradient
        dy_dtensor: Distributed upstream gradient
    """
    # Create DTensor input
    x_dtensor = distribute_tensor(x_global.detach().clone(), device_mesh, input_placements).requires_grad_(True)
    x_dtensor_copy = x_dtensor.detach().clone().requires_grad_(True)

    # Compute on distributed tensor using shardwise_flatten_sharded
    y_dtensor_result = shardwise_flatten_sharded(x_dtensor, start_dim=start_dim, end_dim=end_dim)

    # verify no change to the fwd input
    assert_tensors_identical(x_dtensor.to_local(), x_dtensor_copy.to_local(), check_grad=False, check_grad_fn=False)

    # Distribute the upstream adjoint for backward pass
    dy_dtensor = distribute_tensor(dy_global.detach().clone(), device_mesh, y_dtensor_result.placements)

    # Perform backward pass
    dy_dtensor_copy = dy_dtensor.detach().clone()
    y_dtensor_result.backward(dy_dtensor)

    # verify no change to the bwd input
    assert_tensors_identical(dy_dtensor.to_local(), dy_dtensor_copy.to_local(), check_grad=False, check_grad_fn=False)

    # verify input gradient placements are consistent with input placements
    assert (
        x_dtensor.grad.placements == input_placements
    ), f"{label_test_case} inconsistent input gradient placements with input placements"

    return y_dtensor_result, x_dtensor, dy_dtensor


def parallel_assert_dtensor_flatten_sharded(
    rank: int,
    grid_group_sizes: Dict[str, int],
    device_type: str,
    backend: str,
    env_map: Optional[Dict[str, str]] = None,
):
    """
    Test shardwise_flatten_sharded which flattens dimensions starting from a sharded axis.

    Unlike shardwise_flatten, this function is designed to flatten dimensions that include
    the sharded dimension. DTensor native op doesn't support this, so we only compare
    against the global serial version as reference.
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

    # each rank uses the same seed to generate the same input tensors
    seed_by_rank(0, seed=42)

    size_dp = grid_group_sizes["dp"]
    size_cp = len(manager.group_ranks["cp"])
    size_ring = isqrt(size_cp)
    if size_ring * size_ring != size_cp:
        raise ValueError(f"cp group size {size_cp} is not a square int")

    # Set test parameters - 6D tensor
    # Shape designed so that sharded dims can be evenly divided
    # dim=1 is sharded by dp, dim=3 is sharded by ring (first dim of cp 2D mesh)
    shape = (2, size_dp * 4, 3, size_ring * 6, 5, 4)

    # Test cases: flatten starting from a sharded dimension
    # Input is sharded on dim=1 (by dp) and dim=3 (by ring)
    # shardwise_flatten_sharded requires start_dim to be the sharded dimension

    # Test Case 1: Shard on dim=1, flatten dims 1,2
    # After flatten: shape becomes (2, size_dp*4*3, size_ring*6, 5, 4)
    # The flattened dim (1) remains sharded by dp
    test_cases_dim1 = [
        # (input_placements, start_dim, end_dim, description)
        ((Shard(dim=1), Replicate(), Replicate()), 1, 2, "flatten dims 1,2 sharded on dim=1"),
        ((Shard(dim=1), Replicate(), Replicate()), 1, 3, "flatten dims 1,2,3 sharded on dim=1"),
        ((Shard(dim=1), Replicate(), Replicate()), 1, -1, "flatten dims 1 to end sharded on dim=1"),
    ]

    # Test Case 2: Shard on dim=3, flatten dims 3,4
    test_cases_dim3 = [
        ((Replicate(), Shard(dim=3), Replicate()), 3, 4, "flatten dims 3,4 sharded on dim=3"),
        ((Replicate(), Shard(dim=3), Replicate()), 3, 5, "flatten dims 3,4,5 sharded on dim=3"),
        ((Replicate(), Shard(dim=3), Replicate()), 3, -1, "flatten dims 3 to end sharded on dim=3"),
    ]

    # Test Case 3: Both dim=1 and dim=3 are sharded
    # dim=1 sharded by dp (mesh dim 0), dim=3 sharded by ring (mesh dim 1)
    # Only test flattening from dim=3 onwards, so dim=1's shard placement is unaffected.
    test_cases_both_sharded = [
        # Flatten from dim=3, dim=1 stays sharded at the same position
        ((Shard(dim=1), Shard(dim=3), Replicate()), 3, 4, "flatten dims 3,4 with both dim=1,3 sharded"),
        ((Shard(dim=1), Shard(dim=3), Replicate()), 3, 5, "flatten dims 3,4,5 with both dim=1,3 sharded"),
        ((Shard(dim=1), Shard(dim=3), Replicate()), 3, -1, "flatten dims 3 to end with both dim=1,3 sharded"),
    ]

    # Test Case 4: Both dim=1 and dim=3 are sharded, flatten from dim=1.
    # Flattening at a lower dim shifts the higher shard's placement index.
    # E.g., flatten dims 1,2 removes 1 dim → Shard(dim=3) must become Shard(dim=2).
    # Format: (input_placements, start_dim, end_dim, expected_output_placements, description)
    test_cases_placement_shift = [
        # flatten dims 1,2: removes 1 dim → Shard(3) shifts to Shard(2)
        (
            (Shard(dim=1), Shard(dim=3), Replicate()),
            1,
            2,
            (Shard(dim=1), Shard(dim=2), Replicate()),
            "flatten dims 1,2 shifting Shard(3)->Shard(2)",
        ),
    ]

    all_test_cases = [
        (pl, sd, ed, pl, desc) for pl, sd, ed, desc in test_cases_dim1 + test_cases_dim3 + test_cases_both_sharded
    ] + test_cases_placement_shift

    for input_placements, start_dim, end_dim, expected_output_placements, description in all_test_cases:
        label_test_case = f"{description} (start_dim={start_dim}, end_dim={end_dim})\n"

        # Compute global expectations using standard PyTorch operations
        x_global, y_expected_global, x_grad_expected_global, dy_global = compute_global_expectation(
            shape, start_dim, end_dim, manager.device
        )

        # NOTE: DTensor native op doesn't support flattening involving sharded dimensions,
        # so we skip DTensor native comparison and only use global serial version as reference.

        # Compute shardwise_flatten_sharded forward and backward with validation
        y_dtensor_result, x_dtensor, dy_dtensor = compute_shardwise_flatten_sharded_with_validation(
            x_global, dy_global, manager.device_mesh_subgroups, input_placements, start_dim, end_dim, label_test_case
        )

        # ===================================================================
        # Check output shape and placements
        # ===================================================================
        # Verify output shape matches expected global shape
        assert (
            y_dtensor_result.shape == y_expected_global.shape
        ), f"{label_test_case} output shape mismatch: got {y_dtensor_result.shape}, expected {y_expected_global.shape}"

        # Verify output placements: Shard dims beyond end_dim shift down by
        # (end_dim - start_dim) because those intermediate dims are merged.
        assert y_dtensor_result.placements == expected_output_placements, (
            f"{label_test_case} output placements mismatch: got {y_dtensor_result.placements}, "
            f"expected {expected_output_placements}"
        )

        # ===================================================================
        # Check against global serial expectation
        # ===================================================================
        # Distribute expected output to compare local shards
        y_dtensor_expected = distribute_tensor(
            y_expected_global, manager.device_mesh_subgroups, y_dtensor_result.placements
        )

        # Compare forward result local shards
        assert_tensors_identical(
            y_dtensor_result.to_local().detach(),
            y_dtensor_expected.to_local().detach(),
        )

        # Compare forward result global tensor
        y_result_global = y_dtensor_result.full_tensor()
        assert_tensors_identical(
            y_result_global.detach(),
            y_expected_global.detach(),
        )

        # ===================================================================
        # Check backward pass against global serial expectation
        # ===================================================================
        # Verify input gradient shape
        assert x_dtensor.grad.shape == x_grad_expected_global.shape, (
            f"{label_test_case} input gradient shape mismatch: got {x_dtensor.grad.shape}, "
            f"expected {x_grad_expected_global.shape}"
        )

        # Distribute expected input gradient for local shard comparison
        x_grad_expected_dtensor = distribute_tensor(
            x_grad_expected_global, manager.device_mesh_subgroups, input_placements
        )

        # Compare input gradient local shards
        assert_tensors_identical(
            x_dtensor.grad.to_local().detach(),
            x_grad_expected_dtensor.to_local().detach(),
        )

        # Compare input gradient global tensor
        assert_tensors_identical(
            x_dtensor.grad.full_tensor().detach(),
            x_grad_expected_global.detach(),
        )

    # ===================================================================
    # Test invalid cases that should raise ValueError
    # ===================================================================

    # Create a test tensor for invalid cases
    x_global_invalid = torch.rand(*shape, device=manager.device)

    # Invalid Case 1: start_dim is NOT sharded
    # Input is sharded on dim=1, but we try to flatten starting from dim=0 (not sharded)
    invalid_not_sharded_cases = [
        # (input_placements, start_dim, end_dim, expected_error_pattern)
        ((Shard(dim=1), Replicate(), Replicate()), 0, 1, "input is not sharded along start_dim"),
        ((Shard(dim=1), Replicate(), Replicate()), 2, 3, "input is not sharded along start_dim"),
        ((Replicate(), Shard(dim=3), Replicate()), 0, 2, "input is not sharded along start_dim"),
        ((Replicate(), Shard(dim=3), Replicate()), 4, 5, "input is not sharded along start_dim"),
    ]

    for input_placements, start_dim, end_dim, error_pattern in invalid_not_sharded_cases:
        x_dtensor = distribute_tensor(x_global_invalid.clone(), manager.device_mesh_subgroups, input_placements)
        with pytest.raises(ValueError, match=error_pattern):
            shardwise_flatten_sharded(x_dtensor, start_dim=start_dim, end_dim=end_dim)

    # Invalid Case 2: start_dim > end_dim
    invalid_dim_order_cases = [
        ((Shard(dim=1), Replicate(), Replicate()), 3, 1, "must be <="),
        ((Replicate(), Shard(dim=3), Replicate()), 5, 3, "must be <="),
    ]

    for input_placements, start_dim, end_dim, error_pattern in invalid_dim_order_cases:
        x_dtensor = distribute_tensor(x_global_invalid.clone(), manager.device_mesh_subgroups, input_placements)
        with pytest.raises(ValueError, match=error_pattern):
            shardwise_flatten_sharded(x_dtensor, start_dim=start_dim, end_dim=end_dim)

    # Invalid Case 3: Dimension out of range
    invalid_out_of_range_cases = [
        ((Shard(dim=1), Replicate(), Replicate()), 10, 11, "out of range"),
        ((Shard(dim=1), Replicate(), Replicate()), 1, 10, "out of range"),
    ]

    for input_placements, start_dim, end_dim, error_pattern in invalid_out_of_range_cases:
        x_dtensor = distribute_tensor(x_global_invalid.clone(), manager.device_mesh_subgroups, input_placements)
        with pytest.raises(ValueError, match=error_pattern):
            shardwise_flatten_sharded(x_dtensor, start_dim=start_dim, end_dim=end_dim)

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
def test_dtensor_flatten_sharded(setup_env):
    """Test shardwise_flatten_sharded for flattening dimensions starting from a sharded axis."""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    spawn_multiprocessing(
        parallel_assert_dtensor_flatten_sharded,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
    )
