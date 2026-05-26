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
from boltz.distributed.model.layers.cat_and_chunk import shardwise_cat
from boltz.testing.utils import assert_tensors_identical, seed_by_rank, spawn_multiprocessing


def is_slice_of(x: torch.Tensor, chunks: list[torch.Tensor], dim: int) -> list[bool]:
    x_chunks = x.split([c.shape[dim] for c in chunks], dim=dim)
    return [chunk.is_set_to(x_chunk) for chunk, x_chunk in zip(chunks, x_chunks)]


def compute_global_expectation(shape, num_inputs, dim_to_cat, device):
    """Compute global expectation using standard PyTorch operations."""
    inputs = []
    for i in range(num_inputs):
        # Create slightly different tensors for each input to make the test more robust
        x = torch.rand(*shape, device=device, requires_grad=True)
        inputs.append(x)

    # Compute on global tensors using standard cat operation
    y = torch.cat(inputs, dim=dim_to_cat)

    # Create gradients for backward pass
    dy = torch.rand_like(y)

    # Backward pass on global tensors
    y.backward(dy)

    # Collect input gradients
    input_grads = [x.grad.detach().clone() for x in inputs]

    # check for backward pass view semantics
    is_grad_view_dy = is_slice_of(dy, [x.grad for x in inputs], dim_to_cat)

    return [x.detach().clone() for x in inputs], y.detach().clone(), input_grads, dy.detach().clone(), is_grad_view_dy


def compute_dtensor_native(
    inputs_global: list[torch.Tensor],
    dy_global: torch.Tensor,
    device_mesh: DeviceMesh,
    input_placements: tuple[Placement, ...],
    dim_to_cat: int,
) -> tuple[list[DTensor], DTensor]:
    """Compute DTensor native operations for comparison."""
    # Create DTensor native inputs
    inputs_dtensor = [
        distribute_tensor(x_global.detach().clone(), device_mesh, input_placements).requires_grad_(True)
        for x_global in inputs_global
    ]

    # Forward pass with native DTensor cat operation
    y_dtensor_result = torch.cat(inputs_dtensor, dim=dim_to_cat)

    # Backward pass with native DTensor op
    dy_dtensor = distribute_tensor(dy_global.detach().clone(), device_mesh, y_dtensor_result.placements)
    y_dtensor_result.backward(dy_dtensor)

    inputs_grad_dtensor = [x.grad for x in inputs_dtensor]

    # check for backward pass view semantics
    is_grad_view_dy = is_slice_of(dy_dtensor.to_local(), [x.to_local() for x in inputs_grad_dtensor], dim_to_cat)

    return inputs_grad_dtensor, y_dtensor_result, is_grad_view_dy


def compute_shardwise_cat_with_validation(
    inputs_global: list[torch.Tensor],
    dy_global: torch.Tensor,
    device_mesh: DeviceMesh,
    input_placements: tuple[Placement, ...],
    dim_to_cat: int,
    label_test_case: str,
) -> tuple[DTensor, list[DTensor], DTensor, list[bool]]:
    """
    Compute shardwise_cat forward and backward pass with input validation checks.

    Returns:
        y_dtensor_result: Forward pass result
        inputs_dtensor: Input tensors with computed gradients
        dy_dtensor: Distributed upstream gradient
        is_grad_view_dy_result: View semantics check results
    """
    # Create DTensor inputs
    inputs_dtensor = [
        distribute_tensor(x_global.detach().clone(), device_mesh, input_placements).requires_grad_(True)
        for x_global in inputs_global
    ]
    inputs_dtensor_copy = [x_dtensor.detach().clone().requires_grad_(True) for x_dtensor in inputs_dtensor]

    # Compute on distributed tensors using shardwise_cat
    y_dtensor_result = shardwise_cat(inputs_dtensor, dim_to_cat)

    # verify no change to the fwd inputs
    for x_dtensor, x_dtensor_copy in zip(inputs_dtensor, inputs_dtensor_copy):
        assert_tensors_identical(x_dtensor.to_local(), x_dtensor_copy.to_local(), check_grad=False, check_grad_fn=False)

    # Distribute the upstream adjoint for backward pass
    dy_dtensor = distribute_tensor(dy_global.detach().clone(), device_mesh, y_dtensor_result.placements)

    # Perform backward pass
    dy_dtensor_copy = dy_dtensor.detach().clone()
    y_dtensor_result.backward(dy_dtensor)

    # verify no change to the bwd input
    assert_tensors_identical(dy_dtensor.to_local(), dy_dtensor_copy.to_local(), check_grad=False, check_grad_fn=False)

    # verify input gradient placements are consistent with input placements
    for i, inp in enumerate(inputs_dtensor):
        assert (
            inp.grad.placements == input_placements
        ), f"{label_test_case} inconsistent input {i} gradient placements with input placements"

    # check for backward pass view semantics
    is_grad_view_dy_result = is_slice_of(dy_dtensor.to_local(), [x.grad.to_local() for x in inputs_dtensor], dim_to_cat)

    return y_dtensor_result, inputs_dtensor, dy_dtensor, is_grad_view_dy_result


def parallel_assert_dtensor_cat(
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
    shape = (3, 5, grid_group_sizes["dp"] * 2, 5, size_ring * 4, 5, 3, 2)
    num_inputs = 3  # Number of tensors to concatenate
    # Shard the sequence dimension (dim=1) for input tensors
    # this emulates the sharded single representation in the Boltz model
    input_placements = (Shard(dim=2), Shard(dim=4), Replicate())

    # Test valid dimensions (not sharded)
    valid_dims_to_cat = [0, 1, 3, 5, 6, 7, -1, -2, -3, -5, -7, -8]
    # Test invalid dimensions (sharded): dim 2, 4, -4 (equiv to dim 4), -6 (equiv to dim 2)
    invalid_dims_to_cat = [2, 4, -4, -6]

    # Test valid concatenation dimensions
    for dim_to_cat in valid_dims_to_cat:
        label_test_case = f"for dim {dim_to_cat}\n"

        # Compute global expectations
        inputs_global, y_expected_global, inputs_grad_expected_global, dy_global, is_grad_view_dy_global = (
            compute_global_expectation(shape, num_inputs, dim_to_cat, manager.device)
        )

        # use DTensor native op as an alternative reference
        # NOTE: DTensor native cat's backward pass doesn't guarantee view semantics
        # as dim_to_cat == 7 gives a different view semantic result than dim_to_cat == -1,
        # the latter should be the same as the former because ndim = 8
        inputs_grad_dtensor_native, y_dtensor_result_native, _ = compute_dtensor_native(
            inputs_global, dy_global, manager.device_mesh_subgroups, input_placements, dim_to_cat
        )

        # Compute shardwise_cat forward and backward with validation
        y_dtensor_result, inputs_dtensor, dy_dtensor, is_grad_view_dy_result = compute_shardwise_cat_with_validation(
            inputs_global, dy_global, manager.device_mesh_subgroups, input_placements, dim_to_cat, label_test_case
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

        # compare global tensors between shardwise_cat and native DTensor results
        y_result_global = y_dtensor_result.full_tensor()
        y_result_global_native = y_dtensor_result_native.full_tensor()

        torch.testing.assert_close(
            y_result_global,
            y_result_global_native,
            atol=0,
            rtol=0,
            msg=lambda m: f"{label_test_case} output vs native: {m}",
        )

        # assert input gradients' metadata and values against DTensor native
        for i, (inp, inp_grad_native) in enumerate(zip(inputs_dtensor, inputs_grad_dtensor_native)):
            assert (
                inp.grad.placements == inp_grad_native.placements
            ), f"{label_test_case} input {i} gradient placements mismatch"
            assert inp.grad.shape == inp_grad_native.shape, f"{label_test_case} input {i} gradient shape mismatch"
            assert (
                inp.grad.stride() == inp_grad_native.stride()
            ), f"{label_test_case} input {i} gradient stride mismatch"

            torch.testing.assert_close(
                inp.grad.to_local(),
                inp_grad_native.to_local(),
                atol=0,
                rtol=0,
                msg=lambda m: f"{label_test_case} input {i} gradient mismatch: {m}",
            )

            torch.testing.assert_close(
                inp.grad.full_tensor(),
                inp_grad_native.full_tensor(),
                atol=0,
                rtol=0,
                msg=lambda m: f"{label_test_case} input {i} gradient mismatch: {m}",
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

        # create distributed tensors from global results for local shard comparison
        for i, input_grad_expected in enumerate(inputs_grad_expected_global):
            input_grad_expected_dtensor = distribute_tensor(
                input_grad_expected, manager.device_mesh_subgroups, input_placements
            )

            # compare local shards with expected
            torch.testing.assert_close(
                inputs_dtensor[i].grad.to_local(),
                input_grad_expected_dtensor.to_local(),
                atol=0,
                rtol=0,
                msg=lambda m: f"{label_test_case} input {i} gradient vs global expectation: {m}",
            )

            torch.testing.assert_close(
                inputs_dtensor[i].grad.full_tensor(),
                input_grad_expected,
                atol=0,
                rtol=0,
                msg=lambda m: f"{label_test_case} input {i} gradient vs global expectation: {m}",
            )

            # With explicit shape and stride, DTensor.from_local can't guarantee view semantics
            # assert (
            #     is_grad_view_dy_result[i] == is_grad_view_dy_global[i]
            # ), f"{label_test_case} input {i} backward pass view semantics mismatch"

    # Test invalid concatenation dimensions (should raise ValueError)
    for dim_to_cat in invalid_dims_to_cat:
        label_test_case = f"for invalid dim {dim_to_cat}\n"

        # Compute global expectations (this should work fine)
        inputs_global, _, _, _, _ = compute_global_expectation(shape, num_inputs, dim_to_cat, manager.device)

        # Create DTensor inputs
        inputs_dtensor = []
        for x_global in inputs_global:
            x_dtensor = distribute_tensor(x_global, manager.device_mesh_subgroups, input_placements)
            x_dtensor.requires_grad = True
            inputs_dtensor.append(x_dtensor)

        # This should raise due to sharded dimension
        with pytest.raises(
            NotImplementedError, match=f"Concatenation along dimension {dim_to_cat} shared by device_mesh axis"
        ):
            shardwise_cat(inputs_dtensor, dim_to_cat)

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
def test_dtensor_cat(setup_env):
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    spawn_multiprocessing(
        parallel_assert_dtensor_cat,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
    )
