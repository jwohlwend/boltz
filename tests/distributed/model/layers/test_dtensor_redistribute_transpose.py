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

from boltz.distributed.comm import TransposeComm
from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.layers.redistribute_transpose import redistribute_transpose
from boltz.testing.utils import assert_tensors_identical, seed_by_rank, spawn_multiprocessing


def compute_global_expectation(shape, input_placements, output_placements, dim0, dim1, device, device_mesh_shape):
    """Compute global expectation using standard PyTorch operations."""
    # Create tensor for operations
    x = torch.rand(*shape, device=device, requires_grad=True)

    # Determine the type of operation based on placements and transpose dimensions
    has_redistribute = output_placements is not None
    has_local_transpose = dim0 is not None and dim1 is not None
    is_all_replicate = all(isinstance(p, type(Replicate())) for p in input_placements)

    # Compute on global tensor based on operation semantics:
    # 1. redistribute only: no change to global tensor
    # 2. local transpose only
    # 3. redistribute + local transpose: global tensor transpose
    # 4. all-replicate + local transpose: global tensor transpose
    if has_redistribute and has_local_transpose:
        # Case 3: redistribute + local transpose = global transpose when dim{0, 1} are sharded
        assert Shard(dim0) in input_placements and Shard(dim1) in input_placements
        y = torch.transpose(x, dim0=dim0, dim1=dim1)
    elif not has_redistribute and has_local_transpose and is_all_replicate:
        # Case 4: all-replicate + local transpose = global transpose
        y = torch.transpose(x, dim0=dim0, dim1=dim1)
    elif not has_redistribute and has_local_transpose:
        # Case 2: local transpose only -- dim{0, 1} can't be sharded
        assert Shard(dim0) not in input_placements and Shard(dim1) not in input_placements
        y = torch.transpose(x, dim0=dim0, dim1=dim1)
    else:
        # Case 1: redistribute only or no-op = identity
        y = x

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
    output_placements: Optional[tuple[Placement, ...]],
    dim0: Optional[int],
    dim1: Optional[int],
) -> tuple[DTensor, DTensor]:
    """Compute DTensor native operations for comparison."""
    # Create DTensor native input
    x_dtensor = distribute_tensor(x_global.detach().clone(), device_mesh, input_placements).requires_grad_(True)

    # Apply redistribute if output_placements is specified
    if output_placements is not None:
        if output_placements != x_dtensor.placements:
            y_dtensor_result = x_dtensor.redistribute(device_mesh, output_placements)
        else:
            # this can only work if dim{0, 1} are sharded
            assert Shard(dim0) in x_dtensor.placements and Shard(dim1) in x_dtensor.placements
            # swap the shard placements of dim0 and dim1
            output_placements_ = list(output_placements)
            i_axis_mesh_shard_dim0 = output_placements_.index(Shard(dim0))
            i_axis_mesh_shard_dim1 = output_placements_.index(Shard(dim1))
            output_placements_[i_axis_mesh_shard_dim0] = Shard(dim1)
            output_placements_[i_axis_mesh_shard_dim1] = Shard(dim0)
            y_dtensor_result = x_dtensor.redistribute(device_mesh, tuple(output_placements_))
    else:
        y_dtensor_result = x_dtensor

    # Apply local transpose if dim0 and dim1 are specified
    if dim0 is not None and dim1 is not None:
        y_dtensor_result = torch.transpose(y_dtensor_result, dim0=dim0, dim1=dim1)
        # assert view semantics when no redistribute
        if output_placements is None:
            assert y_dtensor_result.to_local().is_set_to(x_dtensor.to_local().transpose(dim0=dim0, dim1=dim1))

    # Backward pass with native DTensor op
    dy_dtensor = distribute_tensor(dy_global.detach().clone(), device_mesh, y_dtensor_result.placements)
    y_dtensor_result.backward(dy_dtensor)

    # assert view semantics of gradients if no redistribute
    if output_placements is None:
        if dim0 is not None and dim1 is not None:
            assert x_dtensor.grad.to_local().is_set_to(dy_dtensor.to_local().transpose(dim0=dim0, dim1=dim1))
        else:
            # DTensor gradient is not  view of the upstream adjoint for no-op case
            assert not x_dtensor.grad.to_local().is_set_to(dy_dtensor.to_local())

    return x_dtensor.grad, y_dtensor_result


def compute_redistribute_transpose_with_validation(
    x_global: torch.Tensor,
    dy_global: torch.Tensor,
    device_mesh: DeviceMesh,
    input_placements: tuple[Placement, ...],
    output_placements: Optional[tuple[Placement, ...]],
    transpose_comm: Optional[TransposeComm],
    dim0: Optional[int],
    dim1: Optional[int],
    label_test_case: str,
) -> tuple[DTensor, DTensor, DTensor]:
    """
    Compute redistribute_transpose forward and backward pass with input validation checks.

    Returns:
        y_dtensor_result: Forward pass result
        x_dtensor: Input tensor with computed gradient
        dy_dtensor: Distributed upstream gradient
    """
    # Create DTensor input
    x_dtensor = distribute_tensor(x_global.detach().clone(), device_mesh, input_placements).requires_grad_(True)
    x_dtensor_copy = x_dtensor.detach().clone().requires_grad_(True)

    # Compute on distributed tensor using redistribute_transpose
    y_dtensor_result = redistribute_transpose(x_dtensor, transpose_comm, output_placements, dim0, dim1)

    # verify no change to the fwd input
    assert_tensors_identical(x_dtensor.to_local(), x_dtensor_copy.to_local(), check_grad=False, check_grad_fn=False)

    # assert view semantics of fwd result if no redistribute
    if output_placements is None:
        if dim0 is not None and dim1 is not None:
            assert y_dtensor_result.to_local().is_set_to(x_dtensor.to_local().transpose(dim0=dim0, dim1=dim1))
        else:
            assert y_dtensor_result.to_local().is_set_to(x_dtensor.to_local())

    # Distribute the upstream adjoint for backward pass
    dy_dtensor = distribute_tensor(dy_global.detach().clone(), device_mesh, y_dtensor_result.placements)

    # Perform backward pass
    dy_dtensor_copy = dy_dtensor.detach().clone()
    y_dtensor_result.backward(dy_dtensor)

    # verify no change to the bwd input
    assert_tensors_identical(dy_dtensor.to_local(), dy_dtensor_copy.to_local(), check_grad=False, check_grad_fn=False)

    # assert view semantics of bwd result if no redistribute
    if output_placements is None:
        if dim0 is not None and dim1 is not None:
            assert x_dtensor.grad.to_local().is_set_to(dy_dtensor.to_local().transpose(dim0=dim0, dim1=dim1))
        else:
            # For some reason, our implementation also follows the native DTensor's semantics
            # for the no-op case
            assert not x_dtensor.grad.to_local().is_set_to(dy_dtensor.to_local())

    # verify input gradient placements are consistent with input placements
    assert (
        x_dtensor.grad.placements == input_placements
    ), f"{label_test_case} inconsistent input gradient placements with input placements"

    return y_dtensor_result, x_dtensor, dy_dtensor


def parallel_assert_dtensor_redistribute_transpose(
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
    shape = (3, 5, grid_group_sizes["dp"] * 2, 5, size_ring * 4, 5, size_ring * 3, 2)

    # Create TransposeComm for tests that need it
    layout_group_cp = manager.layout_subgroups["cp"]
    transpose_comm = TransposeComm(manager.group["cp"], layout_group_cp)

    invalid_test_cases = [
        (
            (Shard(2), Shard(4), Shard(6)),
            (Shard(2), Shard(4), Shard(6)),
            2,
            6,
            ValueError,
            "Inconsistent device mesh coordinate.*",
            manager.subgroups_rank["cp"][0] != manager.group_rank["dp"],  # otherwise assertion trivially pass
        ),  # cross DP - CP mesh transpose will raise if the CP transpose_comm
        (
            (Shard(2), Shard(4), Shard(6)),
            None,
            4,
            0,
            NotImplementedError,
            "Local transpose on sharded dimensions.*",
            True,
        ),  # local transpose on sharded dimensions will raise
        (
            (Shard(2), Shard(4), Shard(6)),
            (Shard(2), Shard(3), Shard(2)),
            4,
            0,
            ValueError,
            "Simultaneous redistribute and local transpose is only supported.*",
            True,
        ),  # redistribute and local transpose without same output placements will raise
        (
            (Shard(2), Shard(4), Shard(6)),
            (Shard(2), Shard(4), Shard(6)),
            4,
            0,
            ValueError,
            "Both dim0 and dim1 must be sharded.*",
            True,
        ),  # redistribute and local transpose on non-sharded dimensions will raise
        (
            (Shard(2), Shard(4), Shard(6)),
            (Replicate(), Replicate(), Replicate()),
            None,
            None,
            ValueError,
            "Input and output placements are not strictly a permutation of each other.*",
            True,
        ),  # redistribute other than device mesh transpose will raise
    ]

    # Define test cases based on user specification
    # Format: (input_placements, output_placements, dim0, dim1, description)
    valid_test_cases = [
        # Global transpose of sharded dtensor along both sharded axes (implies a redistribute)
        (
            (Shard(2), Shard(4), Shard(6)),
            (Shard(2), Shard(4), Shard(6)),
            4,
            6,
            "global transpose along both sharded axes (4,6)",
        ),
        # Redistribute of sharded single representation but no local transpose
        (
            (Shard(2), Shard(4), Replicate()),
            (Shard(2), Replicate(), Shard(4)),
            None,
            None,
            "redistribute only S(2),S(4),R -> S(2),R,S(4)",
        ),
        (
            (Shard(2), Replicate(), Shard(4)),
            (Shard(2), Shard(4), Replicate()),
            None,
            None,
            "redistribute only S(2),R,S(4) -> S(2),S(4),R",
        ),
        # Local transpose only
        ((Shard(2), Replicate(), Replicate()), None, 0, 3, "local transpose only (0,3)"),
        ((Shard(2), Replicate(), Replicate()), None, 4, 6, "local transpose only (4,6)"),
        # No op
        ((Shard(2), Replicate(), Replicate()), None, None, None, "no op"),
        # Local transpose of all-replicate dtensor implying global transpose
        ((Replicate(), Replicate(), Replicate()), None, 0, 3, "local transpose of all-replicate (0,3)"),
    ]

    for input_placements, output_placements, dim0, dim1, description in valid_test_cases:
        label_test_case = f"for {description}\n"

        # Determine if we need transpose_comm for this test case
        needs_transpose_comm = output_placements is not None
        current_transpose_comm = transpose_comm if needs_transpose_comm else None

        # Compute global expectations
        x_global, y_expected_global, x_grad_expected_global, dy_global = compute_global_expectation(
            shape, input_placements, output_placements, dim0, dim1, manager.device, manager.device_mesh_subgroups.shape
        )

        # Use DTensor native op as an alternative reference
        x_grad_dtensor_native, y_dtensor_result_native = compute_dtensor_native(
            x_global, dy_global, manager.device_mesh_subgroups, input_placements, output_placements, dim0, dim1
        )

        # Compute redistribute_transpose forward and backward with validation
        y_dtensor_result, x_dtensor, dy_dtensor = compute_redistribute_transpose_with_validation(
            x_global,
            dy_global,
            manager.device_mesh_subgroups,
            input_placements,
            output_placements,
            current_transpose_comm,
            dim0,
            dim1,
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

        # compare global tensors between redistribute_transpose and native DTensor results
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

    # Test invalid cases to assert raise
    for input_placements, output_placements, dim0, dim1, error, raise_msg, raise_condition in invalid_test_cases:
        label_test_case = f"for {error.__name__} {raise_msg}\n"

        x_global = torch.rand(*shape, device=manager.device, requires_grad=True)

        # Create DTensor input
        x_dtensor = distribute_tensor(x_global, manager.device_mesh_subgroups, input_placements)
        x_dtensor.requires_grad = True

        # Determine if we need transpose_comm for this test case
        needs_transpose_comm = output_placements is not None
        current_transpose_comm = transpose_comm if needs_transpose_comm else None

        # This should raise due to sharded dimension being unflattened
        if raise_condition:
            with pytest.raises(error, match=raise_msg):
                redistribute_transpose(x_dtensor, current_transpose_comm, output_placements, dim0, dim1)
    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env",
    [
        ((1, (1, 1)), True, "cuda", "ENV"),  # 1 GPU (serial equiv)
        ((2, (1, 1)), True, "cuda", "ENV"),  # 2 GPUs, dp=2, cp=1x1
        ((1, (2, 2)), True, "cuda", "ENV"),
        ((2, (2, 2)), True, "cuda", "ENV"),
        ((1, (3, 3)), True, "cuda", "ENV"),
        ((1, (3, 3)), True, "cpu", "ENV"),
    ],
    indirect=("setup_env",),
    ids=lambda x: f"dp={x[0][0]}, cp={x[0][1]}, specify_method={x[1]}, device_type={x[2]}, method_init={x[3]}",
)
def test_dtensor_redistribute_transpose(setup_env):
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    spawn_multiprocessing(
        parallel_assert_dtensor_redistribute_transpose,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
    )
