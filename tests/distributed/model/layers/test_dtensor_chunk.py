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
from boltz.distributed.model.layers.cat_and_chunk import shardwise_chunk
from boltz.testing.utils import assert_tensors_identical, seed_by_rank, spawn_multiprocessing


def is_slice_of(x: torch.Tensor, chunks: list[torch.Tensor], dim: int) -> list[bool]:
    """Check if each chunk is a view/slice of the original tensor x along the specified dimension."""
    x_chunks = x.split([c.shape[dim] for c in chunks], dim=dim)
    return [chunk.is_set_to(x_chunk) for chunk, x_chunk in zip(chunks, x_chunks)]


def compute_global_expectation(shape, num_chunks, dim_to_chunk, device):
    """Compute global expectation using standard PyTorch operations."""
    # Create input tensor
    x = torch.rand(*shape, device=device, requires_grad=True)

    # Compute on global tensor using standard chunk operation
    y_chunks = x.chunk(num_chunks, dim=dim_to_chunk)

    # Check for forward pass view semantics (are chunks views of the input?)
    is_chunk_view_x = is_slice_of(x, list(y_chunks), dim_to_chunk)

    # Create gradient for backward pass - each chunk gets different gradient
    dy_chunks = [torch.rand_like(chunk) for chunk in y_chunks]

    # Backward pass on global tensors
    # Need to use torch.autograd.backward with multiple tensors
    torch.autograd.backward(y_chunks, dy_chunks)

    # Collect input gradient
    input_grad = x.grad.detach().clone()

    return (
        x.detach().clone(),
        [chunk.detach().clone() for chunk in y_chunks],
        input_grad,
        dy_chunks,
        is_chunk_view_x,
    )


def compute_dtensor_native(
    input_global: torch.Tensor,
    dy_chunks_global: list[torch.Tensor],
    device_mesh: DeviceMesh,
    input_placements: tuple[Placement, ...],
    dim_to_chunk: int,
    num_chunks: int,
) -> tuple[list[DTensor], torch.Tensor, bool]:
    """Compute DTensor native operations for comparison."""
    # Create DTensor native input
    input_dtensor = distribute_tensor(input_global.detach().clone(), device_mesh, input_placements).requires_grad_(True)

    # Forward pass with native DTensor chunk operation
    y_chunks_dtensor_result = torch.chunk(input_dtensor, num_chunks, dim=dim_to_chunk)

    # Convert tuple to list for easier handling
    y_chunks_dtensor_result = list(y_chunks_dtensor_result)

    # Check for forward pass view semantics (are chunks views of the input?)
    is_chunk_view_x = is_slice_of(
        input_dtensor.to_local(), [chunk.to_local() for chunk in y_chunks_dtensor_result], dim_to_chunk
    )

    # Backward pass with native DTensor op
    dy_chunks_dtensor = [
        distribute_tensor(dy_chunk_global.detach().clone(), device_mesh, chunk.placements)
        for dy_chunk_global, chunk in zip(dy_chunks_global, y_chunks_dtensor_result)
    ]

    torch.autograd.backward(y_chunks_dtensor_result, dy_chunks_dtensor)

    input_grad_dtensor = input_dtensor.grad

    return y_chunks_dtensor_result, input_grad_dtensor, is_chunk_view_x


def compute_shardwise_chunk_with_validation(
    input_global: torch.Tensor,
    dy_chunks_global: list[torch.Tensor],
    device_mesh: DeviceMesh,
    input_placements: tuple[Placement, ...],
    dim_to_chunk: int,
    num_chunks: int,
    label_test_case: str,
) -> tuple[list[DTensor], DTensor, list[DTensor], list[bool]]:
    """
    Compute shardwise_chunk forward and backward pass with input validation checks.

    Returns:
        y_chunks_result: Forward pass result (list of chunks)
        input_dtensor: Input tensor with computed gradient
        dy_chunks_dtensor: Distributed upstream gradients
        is_grad_view_dy_result: View semantics check results
    """
    # Create DTensor input
    input_dtensor = distribute_tensor(input_global.detach().clone(), device_mesh, input_placements).requires_grad_(True)
    input_dtensor_copy = input_dtensor.detach().clone().requires_grad_(True)

    # Compute on distributed tensor using shardwise_chunk
    y_chunks_result = shardwise_chunk(input_dtensor, num_chunks, dim_to_chunk)

    # Convert tuple to list for easier handling
    y_chunks_result = list(y_chunks_result)

    # Check for forward pass view semantics (are chunks views of the input?)
    is_chunk_view_x_result = is_slice_of(
        input_dtensor.to_local(), [chunk.to_local() for chunk in y_chunks_result], dim_to_chunk
    )

    # Verify no change to the fwd input
    assert_tensors_identical(
        input_dtensor.to_local(), input_dtensor_copy.to_local(), check_grad=False, check_grad_fn=False
    )

    # Distribute the upstream adjoints for backward pass
    dy_chunks_dtensor = [
        distribute_tensor(dy_chunk_global.detach().clone(), device_mesh, chunk.placements)
        for dy_chunk_global, chunk in zip(dy_chunks_global, y_chunks_result)
    ]

    # Perform backward pass
    dy_chunks_dtensor_copy = [dy_chunk.detach().clone() for dy_chunk in dy_chunks_dtensor]
    torch.autograd.backward(y_chunks_result, dy_chunks_dtensor)

    # Verify no change to the bwd inputs
    for dy_chunk, dy_chunk_copy in zip(dy_chunks_dtensor, dy_chunks_dtensor_copy):
        assert_tensors_identical(dy_chunk.to_local(), dy_chunk_copy.to_local(), check_grad=False, check_grad_fn=False)

    # Verify input gradient placements are consistent with input placements
    assert (
        input_dtensor.grad.placements == input_placements
    ), f"{label_test_case} inconsistent input gradient placements with input placements"

    return y_chunks_result, input_dtensor, dy_chunks_dtensor, is_chunk_view_x_result


def parallel_assert_dtensor_chunk(
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
    shape = (3, 5, grid_group_sizes["dp"] * 2, 5, size_ring * 4, 5, 3, 12)  # Last dim divisible by 2,3,4
    num_chunks = 3  # Number of chunks to split into
    # Shard the sequence dimension (dim=2) and token dimension (dim=4) for input tensor
    # this emulates the sharded representation in the Boltz model
    input_placements = (Shard(2), Shard(4), Replicate())

    # Test valid dimensions (not sharded)
    valid_dims_to_chunk = [0, 1, 3, 5, 6, 7, -1, -2, -3, -5, -7, -8]
    # Test invalid dimensions (sharded): dim 2, 4, -4 (equiv to dim 4), -6 (equiv to dim 2)
    invalid_dims_to_chunk = [2, 4, -4, -6]

    # Test valid chunking dimensions
    for dim_to_chunk in valid_dims_to_chunk:
        label_test_case = f"for dim {dim_to_chunk}\n"

        # Compute global expectations
        input_global, y_chunks_expected_global, input_grad_expected_global, dy_chunks_global, is_chunk_view_x_global = (
            compute_global_expectation(shape, num_chunks, dim_to_chunk, manager.device)
        )

        # Use DTensor native op as an alternative reference
        y_chunks_dtensor_native, input_grad_dtensor_native, is_chunk_view_x_native = compute_dtensor_native(
            input_global, dy_chunks_global, manager.device_mesh_subgroups, input_placements, dim_to_chunk, num_chunks
        )

        # Compute shardwise_chunk forward and backward with validation
        y_chunks_result, input_dtensor, dy_chunks_dtensor, is_chunk_view_x_result = (
            compute_shardwise_chunk_with_validation(
                input_global,
                dy_chunks_global,
                manager.device_mesh_subgroups,
                input_placements,
                dim_to_chunk,
                num_chunks,
                label_test_case,
            )
        )

        # ===================================================================
        # BLOCK 1: Check against DTensor native reference
        # ===================================================================

        # Check metadata against DTensor native - number of chunks
        assert len(y_chunks_result) == len(y_chunks_dtensor_native), f"{label_test_case} number of chunks mismatch"

        # Check each chunk against DTensor native
        for i, (chunk_result, chunk_native) in enumerate(zip(y_chunks_result, y_chunks_dtensor_native)):
            assert (
                chunk_result.placements == chunk_native.placements
            ), f"{label_test_case} chunk {i} placements mismatch"
            assert chunk_result.shape == chunk_native.shape, f"{label_test_case} chunk {i} shape mismatch"
            # In some of the test cases, the DTensor native result will retain the same stride as the input global,
            # which I believe is actually wrong because upon the result.full_tensor(), the would-be padding won't be
            # materialized and it shouldn't be materialized
            # assert chunk_result.stride() == chunk_native.stride(), f"{label_test_case} chunk {i} stride mismatch"

            # Compare forward result with native DTensor op
            torch.testing.assert_close(
                chunk_result.to_local(),
                chunk_native.to_local(),
                atol=0,
                rtol=0,
                msg=lambda m: f"{label_test_case} chunk {i}: {m}",
            )

            # Compare global tensors between shardwise_chunk and native DTensor results
            chunk_result_global = chunk_result.full_tensor()
            chunk_result_global_native = chunk_native.full_tensor()

            torch.testing.assert_close(
                chunk_result_global,
                chunk_result_global_native,
                atol=0,
                rtol=0,
                msg=lambda m: f"{label_test_case} chunk {i} vs native: {m}",
            )

        assert is_chunk_view_x_result == is_chunk_view_x_native, (
            f"{label_test_case} forward pass view semantics mismatch with DTensor native: "
            f"Expected: {is_chunk_view_x_native}, "
            f"Actual: {is_chunk_view_x_result}"
        )

        # Assert input gradient metadata and values against DTensor native
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

        # Compare results with expected local shards
        for i, (chunk_result, chunk_expected) in enumerate(zip(y_chunks_result, y_chunks_expected_global)):
            chunk_dtensor_expected = distribute_tensor(
                chunk_expected, manager.device_mesh_subgroups, chunk_result.placements
            )

            torch.testing.assert_close(
                chunk_result.to_local(),
                chunk_dtensor_expected.to_local(),
                atol=0,
                rtol=0,
                msg=lambda m: f"{label_test_case} chunk {i}: {m}",
            )

            # Compare forward result with global expectation
            torch.testing.assert_close(
                chunk_result.full_tensor(),
                chunk_expected,
                atol=0,
                rtol=0,
                msg=lambda m: f"{label_test_case} chunk {i} vs global expectation: {m}",
            )

        # Create distributed tensor from global result for local shard comparison
        input_grad_expected_dtensor = distribute_tensor(
            input_grad_expected_global, manager.device_mesh_subgroups, input_placements
        )

        # Compare local shards with expected
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

        assert (
            is_chunk_view_x_result == is_chunk_view_x_global
        ), f"{label_test_case} forward pass view semantics mismatch"

    # Test invalid chunking dimensions (should raise NotImplementedError)
    for dim_to_chunk in invalid_dims_to_chunk:
        label_test_case = f"for invalid dim {dim_to_chunk}\n"

        # Compute global expectations (this should work fine)
        input_global, _, _, _, _ = compute_global_expectation(shape, num_chunks, dim_to_chunk, manager.device)

        # Create DTensor input
        input_dtensor = distribute_tensor(input_global, manager.device_mesh_subgroups, input_placements)
        input_dtensor.requires_grad = True

        # This should raise due to sharded dimension
        with pytest.raises(
            NotImplementedError, match=f"Chunking along dimension {dim_to_chunk} shared by device_mesh axis"
        ):
            shardwise_chunk(input_dtensor, num_chunks, dim_to_chunk)

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
def test_dtensor_chunk(setup_env):
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    spawn_multiprocessing(
        parallel_assert_dtensor_chunk,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
    )
