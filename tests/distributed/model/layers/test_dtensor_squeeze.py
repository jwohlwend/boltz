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
from math import isqrt
from typing import Dict, Optional

import pytest
import torch
from torch.distributed.device_mesh import DeviceMesh
from torch.distributed.tensor import DTensor, Placement, Replicate, Shard, distribute_tensor

from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.layers.squeeze import shardwise_squeeze, shardwise_unsqueeze
from boltz.testing.utils import assert_tensors_identical, seed_by_rank, spawn_multiprocessing


def compute_global_expectation_unsqueeze(batch_size, seq_len, feature_dim, dim_to_unsqueeze, device):
    x = torch.rand(batch_size, seq_len, feature_dim, device=device, requires_grad=True)

    # Compute on global tensor using standard unsqueeze operation
    y = x.unsqueeze(dim_to_unsqueeze)

    # Create gradients for backward pass
    dy = torch.rand_like(y)

    # Backward pass on global tensor
    y.backward(dy)

    return x.detach().clone(), y.detach().clone(), x.grad.detach().clone(), dy.detach().clone()


def compute_dtensor_native_unsqueeze(
    x_global: torch.Tensor,
    dy_global: torch.Tensor,
    device_mesh: DeviceMesh,
    input_placements: tuple[Placement, ...],
    dim_to_unsqueeze: int,
) -> tuple[DTensor, DTensor]:
    """Compute DTensor native operations for comparison."""
    # Create DTensor native input
    x_dtensor = distribute_tensor(x_global.detach().clone(), device_mesh, input_placements).requires_grad_(True)

    # Forward pass with native DTensor unsqueeze operation
    y_dtensor_result = x_dtensor.unsqueeze(dim_to_unsqueeze)

    # Backward pass with native DTensor op
    dy_dtensor = distribute_tensor(dy_global.detach().clone(), device_mesh, y_dtensor_result.placements)
    y_dtensor_result.backward(dy_dtensor)

    x_grad_dtensor = x_dtensor.grad

    # do the view check on the DTensor native op
    assert (
        y_dtensor_result.to_local().squeeze(dim_to_unsqueeze).is_set_to(x_dtensor.to_local())
    ), f"for dim {dim_to_unsqueeze} output local shard is not a view of the input shard for native DTensor op"
    # do the view check on the DTensor native op for backward pass

    assert x_grad_dtensor.to_local().is_set_to(
        dy_dtensor.to_local().squeeze(dim_to_unsqueeze)
    ), f"for dim {dim_to_unsqueeze} input grad is not a view of the upstream adjoint for native DTensor op"

    return x_grad_dtensor, y_dtensor_result


def parallel_assert_dtensor_unsqueeze(
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

    # Set test parameters
    batch_size = 2
    seq_len_per_rank = 4
    seq_len_global = size_ring * seq_len_per_rank
    feature_dim = 3
    dims_to_unsqueeze = [0, 1, 2, 3, -1, -2, -3, -4]  # Test various dimensions including negative indexing

    for dim_to_unsqueeze in dims_to_unsqueeze:
        label_test_case = f"for dim {dim_to_unsqueeze}\n"
        # Compute global expectations
        x_global, y_expected_global, x_grad_expected_global, dy_global = compute_global_expectation_unsqueeze(
            batch_size, seq_len_global, feature_dim, dim_to_unsqueeze, manager.device
        )

        # Create distributed tensors
        # Shard the sequence dimension (dim=1) for input tensor
        # this emulates the sharded single representation in the Boltz model
        input_placements = (Shard(dim=0), Shard(dim=1), Replicate())

        # use DTensor native op as an alternative reference
        x_grad_dtensor_native, y_dtensor_result_native = compute_dtensor_native_unsqueeze(
            x_global, dy_global, manager.device_mesh_subgroups, input_placements, dim_to_unsqueeze
        )

        # Create DTensor input
        x_dtensor = distribute_tensor(x_global, manager.device_mesh_subgroups, input_placements).requires_grad_(True)
        x_dtensor_copy = x_dtensor.detach().clone().requires_grad_(True)

        # Compute on distributed tensor using shardwise_unsqueeze
        y_dtensor_result = shardwise_unsqueeze(x_dtensor, dim_to_unsqueeze)

        # check if the output local shard is a view of the input. We know
        # that squeeze() and unsqueeze() guarantees a view of the input
        # so we can use them here to do is_set_to() check, which otherwise
        # wouldn't work because is_set_to() also checks the strides, which
        # are different between the pre-squeeze/unsqueeze and post-squeeze/unsqueeze
        assert (
            y_dtensor_result.to_local().squeeze(dim_to_unsqueeze).is_set_to(x_dtensor.to_local())
        ), f"{label_test_case} output local shard is not a view of the input shard"

        # verify no change to the fwd input
        assert_tensors_identical(x_dtensor.to_local(), x_dtensor_copy.to_local(), check_grad=False, check_grad_fn=False)

        # Distribute the upstream adjoint for backward pass
        dy_dtensor = distribute_tensor(dy_global, manager.device_mesh_subgroups, y_dtensor_result.placements)

        # Perform backward pass
        dy_dtensor_copy = dy_dtensor.detach().clone()
        y_dtensor_result.backward(dy_dtensor)

        # check if the input grad is a view of the upstream adjoint
        assert x_dtensor.grad.to_local().is_set_to(
            dy_dtensor.to_local().squeeze(dim_to_unsqueeze)
        ), f"{label_test_case} input grad is not a view of the upstream adjoint"

        # verify no change to the bwd input
        assert_tensors_identical(
            dy_dtensor.to_local(), dy_dtensor_copy.to_local(), check_grad=False, check_grad_fn=False
        )

        # verify input gradient placements are consistent with input placements
        assert (
            x_dtensor.grad.placements == input_placements
        ), f"{label_test_case} inconsistent input gradient placements with input placements"

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

        # assert input gradients' metadata and values against DTensor native
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

        # compare global tensors between shardwise_unsqueeze and native DTensor results
        y_result_global = y_dtensor_result.full_tensor()
        y_result_global_native = y_dtensor_result_native.full_tensor()

        torch.testing.assert_close(
            y_result_global,
            y_result_global_native,
            atol=0,
            rtol=0,
            msg=lambda m: f"{label_test_case} output vs native: {m}",
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

        # create distributed tensor from global results for local shard comparison
        x_grad_dtensor_expected = distribute_tensor(
            x_grad_expected_global, manager.device_mesh_subgroups, input_placements
        )

        # compare local shards with expected
        torch.testing.assert_close(
            x_dtensor.grad.to_local(),
            x_grad_dtensor_expected.to_local(),
            atol=0,
            rtol=0,
            msg=lambda m: f"{label_test_case} input gradient vs global expectation: {m}",
        )

        # compare global gradients with serial expectation
        torch.testing.assert_close(
            x_dtensor.grad.full_tensor(),
            x_grad_expected_global,
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
    ids=lambda x: f"dp={x[0][0]}, cp={x[0][1]}, specify_method={x[1]}, device_type={x[2]}, method_init={x[3]}",
)
def test_dtensor_unsqueeze(setup_env):
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    spawn_multiprocessing(
        parallel_assert_dtensor_unsqueeze,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
    )


def compute_global_expectation_squeeze(batch_size, seq_len, dim_to_squeeze, device):
    # Create tensor with singleton dimensions that can be squeezed
    # Shape will be (batch_size, seq_len) for 2D tensor matching device mesh
    x = torch.rand(batch_size, seq_len, device=device)
    x = x.unsqueeze(dim_to_squeeze)
    x.requires_grad_(True)

    # Compute on global tensor using standard squeeze operation
    y = x.squeeze(dim_to_squeeze)

    # Create gradients for backward pass
    dy = torch.rand_like(y)

    # Backward pass on global tensor
    y.backward(dy)

    return x.detach().clone(), y.detach().clone(), x.grad.detach().clone(), dy.detach().clone()


def compute_dtensor_native_squeeze(
    x_global: torch.Tensor,
    dy_global: torch.Tensor,
    device_mesh: DeviceMesh,
    input_placements: tuple[Placement, ...],
    dim_to_squeeze: int,
) -> tuple[DTensor, DTensor]:
    """Compute DTensor native operations for comparison."""
    # Create DTensor native input
    x_dtensor = distribute_tensor(x_global.detach().clone(), device_mesh, input_placements).requires_grad_(True)

    # Forward pass with native DTensor squeeze operation
    y_dtensor_result = x_dtensor.squeeze(dim_to_squeeze)

    # Backward pass with native DTensor op
    dy_dtensor = distribute_tensor(dy_global.detach().clone(), device_mesh, y_dtensor_result.placements)
    y_dtensor_result.backward(dy_dtensor)

    x_grad_dtensor = x_dtensor.grad

    # do the view check on the DTensor native op
    assert (
        y_dtensor_result.to_local().unsqueeze(dim_to_squeeze).is_set_to(x_dtensor.to_local())
    ), f"for dim {dim_to_squeeze} output local shard is not a view of the input shard for native DTensor op"
    # do the view check on the DTensor native op for backward pass

    assert x_grad_dtensor.to_local().is_set_to(
        dy_dtensor.to_local().unsqueeze(dim_to_squeeze)
    ), f"for dim {dim_to_squeeze} input grad is not a view of the upstream adjoint for native DTensor op"

    return x_grad_dtensor, y_dtensor_result


def parallel_assert_dtensor_squeeze(
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

    # Set test parameters
    batch_size = 2
    seq_len_per_rank = 4
    seq_len_global = size_ring * seq_len_per_rank
    dims_to_squeeze = [0, 1, 2, -1, -2, -3]

    for dim_to_squeeze in dims_to_squeeze:
        label_test_case = f"for dim {dim_to_squeeze}\n"
        # Compute global expectations
        x_global, y_expected_global, x_grad_expected_global, dy_global = compute_global_expectation_squeeze(
            batch_size, seq_len_global, dim_to_squeeze, manager.device
        )

        # Create distributed tensors
        # Shard the batch and sequence dimensions for input tensor
        # Input shape is (batch_size, seq_len, 1)
        input_placements = (Shard(dim=0), Shard(dim=1), Replicate())

        # use DTensor native op as an alternative reference
        wrapped_dim_to_squeeze = dim_to_squeeze if dim_to_squeeze >= 0 else dim_to_squeeze + x_global.ndim
        dim_is_sharded = any(
            placement
            for placement in input_placements
            if isinstance(placement, Shard) and placement.dim == wrapped_dim_to_squeeze
        )

        if not dim_is_sharded:
            x_grad_dtensor_native, y_dtensor_result_native = compute_dtensor_native_squeeze(
                x_global, dy_global, manager.device_mesh_subgroups, input_placements, dim_to_squeeze
            )

        # Create DTensor input
        x_dtensor = distribute_tensor(x_global, manager.device_mesh_subgroups, input_placements).requires_grad_(True)
        x_dtensor_copy = x_dtensor.detach().clone().requires_grad_(True)

        # Compute on distributed tensor using shardwise_squeeze
        if dim_is_sharded:  # short circuit if squeeze on sharded dimensions
            with pytest.raises(ValueError, match=r"Cannot squeeze dimension .* as it is sharded"):
                y_dtensor_result = shardwise_squeeze(x_dtensor, dim_to_squeeze)
            continue

        y_dtensor_result = shardwise_squeeze(x_dtensor, dim_to_squeeze)

        # check if the output local shard is a view of the input. We know
        # that squeeze() and unsqueeze() guarantees a view of the input
        # so we can use them here to do is_set_to() check, which otherwise
        # wouldn't work because is_set_to() also checks the strides, which
        # are different between the pre-squeeze/unsqueeze and post-squeeze/unsqueeze
        assert (
            y_dtensor_result.to_local().unsqueeze(dim_to_squeeze).is_set_to(x_dtensor.to_local())
        ), f"{label_test_case} output local shard is not a view of the input shard"

        # verify no change to the fwd input
        assert_tensors_identical(x_dtensor.to_local(), x_dtensor_copy.to_local(), check_grad=False, check_grad_fn=False)

        # Distribute the upstream adjoint for backward pass
        dy_dtensor = distribute_tensor(dy_global, manager.device_mesh_subgroups, y_dtensor_result.placements)

        # Perform backward pass
        dy_dtensor_copy = dy_dtensor.detach().clone()
        y_dtensor_result.backward(dy_dtensor)

        # check if the input grad is a view of the upstream adjoint
        assert x_dtensor.grad.to_local().is_set_to(
            dy_dtensor.to_local().unsqueeze(dim_to_squeeze)
        ), f"{label_test_case} input grad is not a view of the upstream adjoint"

        # verify no change to the bwd input
        assert_tensors_identical(
            dy_dtensor.to_local(), dy_dtensor_copy.to_local(), check_grad=False, check_grad_fn=False
        )

        # verify input gradient placements are consistent with input placements
        assert (
            x_dtensor.grad.placements == input_placements
        ), f"{label_test_case} inconsistent input gradient placements with input placements"

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

        # assert input gradients' metadata and values against DTensor native
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

        # compare global tensors between shardwise_squeeze and native DTensor results
        y_result_global = y_dtensor_result.full_tensor()
        y_result_global_native = y_dtensor_result_native.full_tensor()

        torch.testing.assert_close(
            y_result_global,
            y_result_global_native,
            atol=0,
            rtol=0,
            msg=lambda m: f"{label_test_case} output vs native: {m}",
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

        # create distributed tensor from global results for local shard comparison
        x_grad_dtensor_expected = distribute_tensor(
            x_grad_expected_global, manager.device_mesh_subgroups, input_placements
        )

        # compare local shards with expected
        torch.testing.assert_close(
            x_dtensor.grad.to_local(),
            x_grad_dtensor_expected.to_local(),
            atol=0,
            rtol=0,
            msg=lambda m: f"{label_test_case} input gradient vs global expectation: {m}",
        )

        # compare global gradients with serial expectation
        torch.testing.assert_close(
            x_dtensor.grad.full_tensor(),
            x_grad_expected_global,
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
    ids=lambda x: f"dp={x[0][0]}, cp={x[0][1]}, specify_method={x[1]}, device_type={x[2]}, method_init={x[3]}",
)
def test_dtensor_squeeze(setup_env):
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    spawn_multiprocessing(
        parallel_assert_dtensor_squeeze,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
    )


if __name__ == "__main__":
    unittest.main()
