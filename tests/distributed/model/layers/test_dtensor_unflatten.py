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
from boltz.distributed.model.layers.flatten_and_unflatten import shardwise_unflatten, shardwise_unflatten_sharded
from boltz.testing.utils import assert_tensors_identical, seed_by_rank, spawn_multiprocessing


def compute_global_expectation(shape, dim, sizes, device):
    """Compute global expectation using standard PyTorch operations."""
    # Create tensor for unflattening
    x = torch.rand(*shape, device=device, requires_grad=True)

    # Compute on global tensor using standard unflatten operation
    y = torch.unflatten(x, dim=dim, sizes=sizes)

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
    dim: int,
    sizes: tuple[int, ...],
) -> tuple[DTensor, DTensor]:
    """Compute DTensor native operations for comparison."""
    # Create DTensor native input
    x_dtensor = distribute_tensor(x_global.detach().clone(), device_mesh, input_placements).requires_grad_(True)

    # Forward pass with native DTensor unflatten operation
    y_dtensor_result = torch.unflatten(x_dtensor, dim=dim, sizes=sizes)

    # Backward pass with native DTensor op
    dy_dtensor = distribute_tensor(dy_global.detach().clone(), device_mesh, y_dtensor_result.placements)
    y_dtensor_result.backward(dy_dtensor)

    x_grad_dtensor = x_dtensor.grad

    return x_grad_dtensor, y_dtensor_result


def compute_shardwise_unflatten_with_validation(
    x_global: torch.Tensor,
    dy_global: torch.Tensor,
    device_mesh: DeviceMesh,
    input_placements: tuple[Placement, ...],
    dim: int,
    sizes: tuple[int, ...],
    label_test_case: str,
) -> tuple[DTensor, DTensor, DTensor]:
    """
    Compute shardwise_unflatten forward and backward pass with input validation checks.

    Returns:
        y_dtensor_result: Forward pass result
        x_dtensor: Input tensor with computed gradient
        dy_dtensor: Distributed upstream gradient
    """
    # Create DTensor input
    x_dtensor = distribute_tensor(x_global.detach().clone(), device_mesh, input_placements).requires_grad_(True)
    x_dtensor_copy = x_dtensor.detach().clone().requires_grad_(True)

    # Compute on distributed tensor using shardwise_unflatten
    y_dtensor_result = shardwise_unflatten(x_dtensor, dim=dim, sizes=sizes)

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


def parallel_assert_dtensor_unflatten(
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

    # Test valid unflattening dimensions (not sharded)
    # Sharded dims are 2 and 4, so valid unflatten dims must not be these dimensions
    # For each dim, we need sizes that multiply to the original dimension size
    valid_unflatten_params = [
        (0, (1, 3)),  # unflatten dim 0 (size 3) into (1, 3)
        (0, (3, 1)),  # unflatten dim 0 (size 3) into (3, 1)
        (1, (1, 5)),  # unflatten dim 1 (size 5) into (1, 5)
        (1, (5, 1)),  # unflatten dim 1 (size 5) into (5, 1)
        (3, (1, 5)),  # unflatten dim 3 (size 5) into (1, 5)
        (5, (1, 5)),  # unflatten dim 5 (size 5) into (1, 5)
        (6, (1, 3)),  # unflatten dim 6 (size 3) into (1, 3)
        (6, (3, 1)),  # unflatten dim 6 (size 3) into (3, 1)
        (7, (1, 2)),  # unflatten dim 7 (size 2) into (1, 2)
        (7, (2, 1)),  # unflatten dim 7 (size 2) into (2, 1)
        (-1, (1, 2)),  # unflatten dim -1 (last dim, size 2) into (1, 2)
        (-2, (1, 3)),  # unflatten dim -2 (dim 6, size 3) into (1, 3)
    ]

    # Test invalid unflattening dimensions (sharded dims 2 and 4)
    invalid_unflatten_params = [
        (2, (1, grid_group_sizes["dp"] * 2)),  # unflatten sharded dim 2
        (2, (grid_group_sizes["dp"], 2)),  # unflatten sharded dim 2 (different split)
        (4, (1, size_ring * 4)),  # unflatten sharded dim 4
        (4, (size_ring, 4)),  # unflatten sharded dim 4 (different split)
        (-6, (1, grid_group_sizes["dp"] * 2)),  # unflatten dim -6 (equivalent to dim 2)
        (-4, (1, size_ring * 4)),  # unflatten dim -4 (equivalent to dim 4)
    ]

    # Test valid unflattening dimensions
    for dim, sizes in valid_unflatten_params:
        label_test_case = f"for dim={dim}, sizes={sizes}\n"

        # Compute global expectations
        x_global, y_expected_global, x_grad_expected_global, dy_global = compute_global_expectation(
            shape, dim, sizes, manager.device
        )

        # use DTensor native op as an alternative reference
        x_grad_dtensor_native, y_dtensor_result_native = compute_dtensor_native(
            x_global, dy_global, manager.device_mesh_subgroups, input_placements, dim, sizes
        )

        # Compute shardwise_unflatten forward and backward with validation
        y_dtensor_result, x_dtensor, dy_dtensor = compute_shardwise_unflatten_with_validation(
            x_global, dy_global, manager.device_mesh_subgroups, input_placements, dim, sizes, label_test_case
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

        # compare global tensors between shardwise_unflatten and native DTensor results
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

    # Test invalid unflattening dimensions (should raise NotImplementedError)
    for dim, sizes in invalid_unflatten_params:
        label_test_case = f"for invalid dim={dim}, sizes={sizes}\n"

        # Compute global expectations (this should work fine)
        x_global, _, _, _ = compute_global_expectation(shape, dim, sizes, manager.device)

        # Create DTensor input
        x_dtensor = distribute_tensor(x_global, manager.device_mesh_subgroups, input_placements)
        x_dtensor.requires_grad = True

        # This should raise due to sharded dimension being unflattened
        with pytest.raises(NotImplementedError, match="Unflattening dimension .* shared by device_mesh axis"):
            shardwise_unflatten(x_dtensor, dim=dim, sizes=sizes)

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
def test_dtensor_unflatten(setup_env):
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    spawn_multiprocessing(
        parallel_assert_dtensor_unflatten,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
    )


# ==============================================================================
# Tests for shardwise_unflatten_sharded
# ==============================================================================


def compute_shardwise_unflatten_sharded_with_validation(
    x_global: torch.Tensor,
    dy_global: torch.Tensor,
    device_mesh: DeviceMesh,
    input_placements: tuple[Placement, ...],
    dim: int,
    sizes: tuple[int, ...],
    label_test_case: str,
) -> tuple[DTensor, DTensor, DTensor]:
    """
    Compute shardwise_unflatten_sharded forward and backward pass with input validation checks.

    Returns:
        y_dtensor_result: Forward pass result
        x_dtensor: Input tensor with computed gradient
        dy_dtensor: Distributed upstream gradient
    """
    # Create DTensor input
    x_dtensor = distribute_tensor(x_global.detach().clone(), device_mesh, input_placements).requires_grad_(True)
    x_dtensor_copy = x_dtensor.detach().clone().requires_grad_(True)

    # Compute on distributed tensor using shardwise_unflatten_sharded
    y_dtensor_result = shardwise_unflatten_sharded(x_dtensor, axis=dim, sizes=sizes)

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


def parallel_assert_dtensor_unflatten_sharded(
    rank: int,
    grid_group_sizes: Dict[str, int],
    device_type: str,
    backend: str,
    env_map: Optional[Dict[str, str]] = None,
):
    """
    Test shardwise_unflatten_sharded which unflattens a sharded dimension.

    Unlike shardwise_unflatten, this function is designed to unflatten a dimension that is
    itself sharded. DTensor native op doesn't support this, so we only compare
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

    # Set test parameters - 5D tensor
    # Shape designed so that sharded dims can be evenly divided
    # dim=1 is sharded by dp with size (size_dp * 4 * 3)
    # dim=2 is sharded by ring with size (size_ring * 6 * 2)
    shape = (2, size_dp * 4 * 3, size_ring * 6 * 2, 5, 4)

    # Test cases: unflatten a sharded dimension
    # shardwise_unflatten_sharded requires dim to be the sharded dimension
    # and sizes[0] must be evenly shardable by the device mesh size

    # Test Case 1: Shard on dim=1, unflatten into 2 dims
    # shape[1] = size_dp * 4 * 3 = size_dp * 12
    # unflatten into (size_dp * 4, 3) -> sizes[0] = size_dp * 4 is evenly shardable by size_dp
    test_cases_dim1 = [
        # (input_placements, dim, sizes, description)
        ((Shard(dim=1), Replicate(), Replicate()), 1, (size_dp * 4, 3), "unflatten dim=1 into (size_dp*4, 3)"),
        ((Shard(dim=1), Replicate(), Replicate()), 1, (size_dp * 2, 6), "unflatten dim=1 into (size_dp*2, 6)"),
        ((Shard(dim=1), Replicate(), Replicate()), 1, (size_dp * 12, 1), "unflatten dim=1 into (size_dp*12, 1)"),
        ((Shard(dim=1), Replicate(), Replicate()), 1, (size_dp * 4, 3, 1), "unflatten dim=1 into (size_dp*4, 3, 1)"),
        ((Shard(dim=1), Replicate(), Replicate()), 1, (size_dp * 2, 2, 3), "unflatten dim=1 into (size_dp*2, 2, 3)"),
    ]

    # Test Case 2: Shard on dim=2, unflatten into 2 dims
    # shape[2] = size_ring * 6 * 2 = size_ring * 12
    # unflatten into (size_ring * 6, 2) -> sizes[0] = size_ring * 6 is evenly shardable by size_ring
    test_cases_dim2 = [
        ((Replicate(), Shard(dim=2), Replicate()), 2, (size_ring * 6, 2), "unflatten dim=2 into (size_ring*6, 2)"),
        ((Replicate(), Shard(dim=2), Replicate()), 2, (size_ring * 4, 3), "unflatten dim=2 into (size_ring*4, 3)"),
        ((Replicate(), Shard(dim=2), Replicate()), 2, (size_ring * 12, 1), "unflatten dim=2 into (size_ring*12, 1)"),
        (
            (Replicate(), Shard(dim=2), Replicate()),
            2,
            (size_ring * 6, 2, 1),
            "unflatten dim=2 into (size_ring*6, 2, 1)",
        ),
        (
            (Replicate(), Shard(dim=2), Replicate()),
            2,
            (size_ring * 2, 3, 2),
            "unflatten dim=2 into (size_ring*2, 3, 2)",
        ),
    ]

    # Test Case 3: Both dim=1 and dim=2 are sharded
    # dim=1 sharded by dp (mesh dim 0), dim=2 sharded by ring (mesh dim 1)
    # unflatten dim=2 so dim=1's shard placement is unaffected
    test_cases_both_sharded = [
        (
            (Shard(dim=1), Shard(dim=2), Replicate()),
            2,
            (size_ring * 6, 2),
            "unflatten dim=2 with both dim=1,2 sharded",
        ),
        (
            (Shard(dim=1), Shard(dim=2), Replicate()),
            2,
            (size_ring * 4, 3),
            "unflatten dim=2 into (size_ring*4, 3) with both sharded",
        ),
        (
            (Shard(dim=1), Shard(dim=2), Replicate()),
            2,
            (size_ring * 2, 3, 2),
            "unflatten dim=2 into (size_ring*2, 3, 2) with both sharded",
        ),
    ]

    # Test Case 4: Both dim=1 and dim=2 are sharded, unflatten dim=1.
    # Unflattening at a lower dim shifts the higher shard's placement index.
    # E.g., unflatten dim=1 into (a, b) adds 1 dim → Shard(dim=2) must become Shard(dim=3).
    # Format: (input_placements, dim, sizes, expected_output_placements, description)
    test_cases_placement_shift = [
        # unflatten dim=1 into 2 parts: adds 1 dim → Shard(2) shifts to Shard(3)
        (
            (Shard(dim=1), Shard(dim=2), Replicate()),
            1,
            (size_dp * 4, 3),
            (Shard(dim=1), Shard(dim=3), Replicate()),
            "unflatten dim=1 into (size_dp*4, 3) shifting Shard(2)->Shard(3)",
        ),
        # unflatten dim=1 into 3 parts: adds 2 dims → Shard(2) shifts to Shard(4)
        (
            (Shard(dim=1), Shard(dim=2), Replicate()),
            1,
            (size_dp * 2, 2, 3),
            (Shard(dim=1), Shard(dim=4), Replicate()),
            "unflatten dim=1 into (size_dp*2, 2, 3) shifting Shard(2)->Shard(4)",
        ),
    ]

    all_test_cases = [
        (pl, d, s, pl, desc) for pl, d, s, desc in test_cases_dim1 + test_cases_dim2 + test_cases_both_sharded
    ] + test_cases_placement_shift

    for input_placements, dim, sizes, expected_output_placements, description in all_test_cases:
        label_test_case = f"{description} (dim={dim}, sizes={sizes})\n"

        # Compute global expectations using standard PyTorch operations
        x_global, y_expected_global, x_grad_expected_global, dy_global = compute_global_expectation(
            shape, dim, sizes, manager.device
        )

        # NOTE: DTensor native op doesn't support unflattening a sharded dimension,
        # so we skip DTensor native comparison and only use global serial version as reference.

        # Compute shardwise_unflatten_sharded forward and backward with validation
        y_dtensor_result, x_dtensor, dy_dtensor = compute_shardwise_unflatten_sharded_with_validation(
            x_global, dy_global, manager.device_mesh_subgroups, input_placements, dim, sizes, label_test_case
        )

        # ===================================================================
        # Check output shape and placements
        # ===================================================================
        # Verify output shape matches expected global shape
        assert (
            y_dtensor_result.shape == y_expected_global.shape
        ), f"{label_test_case} output shape mismatch: got {y_dtensor_result.shape}, expected {y_expected_global.shape}"

        # Verify output placements: Shard dims beyond the unflatten point shift up
        # by (len(sizes) - 1) because that many new dims are introduced.
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

    # Invalid Case 1: dim is NOT sharded
    invalid_not_sharded_cases = [
        # (input_placements, dim, sizes, expected_error_pattern)
        ((Shard(dim=1), Replicate(), Replicate()), 0, (1, 2), "input is not sharded along dim"),
        ((Shard(dim=1), Replicate(), Replicate()), 2, (size_ring * 6, 2), "input is not sharded along dim"),
        ((Replicate(), Shard(dim=2), Replicate()), 0, (1, 2), "input is not sharded along dim"),
        ((Replicate(), Shard(dim=2), Replicate()), 1, (size_dp * 4, 3), "input is not sharded along dim"),
    ]

    for input_placements, dim, sizes, error_pattern in invalid_not_sharded_cases:
        x_dtensor = distribute_tensor(x_global_invalid.clone(), manager.device_mesh_subgroups, input_placements)
        with pytest.raises(ValueError, match=error_pattern):
            shardwise_unflatten_sharded(x_dtensor, axis=dim, sizes=sizes)

    # Invalid Case 2: sizes[0] not evenly shardable
    # Need tensors with shapes that allow sizes product to match but sizes[0] not evenly shardable
    # We need to pick sizes[0] that is NOT divisible by the mesh size
    # For size_dp=2: use 3 (3 % 2 != 0), shape[1] = 2*12 = 24 = 3*8
    # For size_ring=2: use 3 (3 % 2 != 0), shape[2] = 2*12 = 24 = 3*8
    # For size_ring=3: use 4 (4 % 3 != 0), shape[2] = 3*12 = 36 = 4*9

    # Test 1: dim=1 sharded by size_dp, sizes[0] not evenly shardable (when size_dp>1 and 3 % size_dp != 0)
    if size_dp > 1 and 3 % size_dp != 0:
        shape_for_uneven1 = (2, size_dp * 12, size_ring * 12, 5, 4)
        x_for_uneven1 = torch.rand(*shape_for_uneven1, device=manager.device)
        x_dtensor_uneven1 = distribute_tensor(
            x_for_uneven1.clone(), manager.device_mesh_subgroups, (Shard(dim=1), Replicate(), Replicate())
        )
        # sizes = (3, size_dp*4) so 3 * size_dp*4 = size_dp*12 = shape[1], but 3 % size_dp != 0
        with pytest.raises(ValueError, match="must be evenly sharded"):
            shardwise_unflatten_sharded(x_dtensor_uneven1, axis=1, sizes=(3, size_dp * 4))

    # Test 2: dim=2 sharded by size_ring, sizes[0] not evenly shardable
    if size_ring > 1:
        # Choose uneven_val such that uneven_val % size_ring != 0
        # For size_ring=2: use 3, for size_ring=3: use 4
        if size_ring == 2:
            uneven_val = 3
            other_factor = 8  # 3 * 8 = 24 = size_ring * 12
        elif size_ring == 3:
            uneven_val = 4
            other_factor = 9  # 4 * 9 = 36 = size_ring * 12
        else:
            # General case: use (size_ring + 1) if it divides size_ring * 12
            # This test may be skipped for some unusual mesh sizes
            uneven_val = size_ring + 1
            product = size_ring * 12
            if product % uneven_val == 0:
                other_factor = product // uneven_val
            else:
                uneven_val = None  # Skip this test

        if uneven_val is not None and uneven_val % size_ring != 0:
            shape_for_uneven2 = (2, size_dp * 12, size_ring * 12, 5, 4)
            x_for_uneven2 = torch.rand(*shape_for_uneven2, device=manager.device)
            x_dtensor_uneven2 = distribute_tensor(
                x_for_uneven2.clone(), manager.device_mesh_subgroups, (Replicate(), Shard(dim=2), Replicate())
            )
            with pytest.raises(ValueError, match="must be evenly sharded"):
                shardwise_unflatten_sharded(x_dtensor_uneven2, axis=2, sizes=(uneven_val, other_factor))

    # Invalid Case 3: Dimension out of range
    invalid_out_of_range_cases = [
        ((Shard(dim=1), Replicate(), Replicate()), 10, (1, 2), "out of range"),
    ]

    for input_placements, dim, sizes, error_pattern in invalid_out_of_range_cases:
        x_dtensor = distribute_tensor(x_global_invalid.clone(), manager.device_mesh_subgroups, input_placements)
        with pytest.raises(ValueError, match=error_pattern):
            shardwise_unflatten_sharded(x_dtensor, axis=dim, sizes=sizes)

    # Invalid Case 4: sizes has less than 2 elements
    invalid_sizes_cases = [
        ((Shard(dim=1), Replicate(), Replicate()), 1, (size_dp * 12,), "at least two dimensions"),
    ]

    for input_placements, dim, sizes, error_pattern in invalid_sizes_cases:
        x_dtensor = distribute_tensor(x_global_invalid.clone(), manager.device_mesh_subgroups, input_placements)
        with pytest.raises(ValueError, match=error_pattern):
            shardwise_unflatten_sharded(x_dtensor, axis=dim, sizes=sizes)

    # Invalid Case 5: Product of sizes doesn't match dim size
    invalid_product_cases = [
        ((Shard(dim=1), Replicate(), Replicate()), 1, (size_dp * 4, 5), "Expected size"),
    ]

    for input_placements, dim, sizes, error_pattern in invalid_product_cases:
        x_dtensor = distribute_tensor(x_global_invalid.clone(), manager.device_mesh_subgroups, input_placements)
        with pytest.raises(ValueError, match=error_pattern):
            shardwise_unflatten_sharded(x_dtensor, axis=dim, sizes=sizes)

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
def test_dtensor_unflatten_sharded(setup_env):
    """Test shardwise_unflatten_sharded for unflattening a sharded dimension."""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    spawn_multiprocessing(
        parallel_assert_dtensor_unflatten_sharded,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
    )
