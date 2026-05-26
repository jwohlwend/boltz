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
from boltz.distributed.model.layers.outer_op import OuterOp, replicate_to_shard_outer_op
from boltz.testing.utils import assert_tensors_identical, seed_by_rank, spawn_multiprocessing


def compute_global_expectation(shape_input, axis, op: OuterOp, device, asymmetric: bool = False):
    """Compute global expectation using standard PyTorch operations."""
    # Create input tensor
    if op == OuterOp.BITAND:
        input_tensor = torch.randint(0, 2, shape_input, device=device, dtype=torch.bool)
    else:
        input_tensor = torch.rand(*shape_input, device=device, requires_grad=True)

    if asymmetric:
        if op == OuterOp.BITAND:
            input_t = torch.randint(0, 2, shape_input, device=device, dtype=torch.bool)
        else:
            input_t = torch.rand(*shape_input, device=device, requires_grad=True)
    else:
        input_t = input_tensor

        # Compute on global tensors using native PyTorch operations
    # Replicate the logic from distributed_outer_op for non-distributed case
    input_expanded = input_tensor.unsqueeze(axis + 1)
    input_t_expanded = input_t.unsqueeze(axis + 1)
    input_t_transposed = input_t_expanded.transpose(axis, axis + 1)

    if op == OuterOp.SUM:
        y = input_expanded + input_t_transposed
    elif op == OuterOp.SUBTRACT:
        y = input_expanded - input_t_transposed
    elif op == OuterOp.EQUAL:
        y = input_expanded == input_t_transposed
    elif op == OuterOp.BITAND:
        y = input_expanded & input_t_transposed
    elif op == OuterOp.PROD:
        y = input_expanded * input_t_transposed
    elif op == OuterOp.CDIST:
        y = torch.cdist(input_tensor, input_t, p=2)

    if op == OuterOp.EQUAL or op == OuterOp.BITAND:
        # Boolean output can't be backpropagated
        return (
            input_tensor.detach().clone(),
            input_t.detach().clone() if asymmetric else None,
            y.detach().clone(),
            None,  # input_grad
            None,  # input_t_grad
            None,  # dy
        )
    else:
        # Create gradients for backward pass
        dy = torch.rand_like(y)

        # Backward pass on global tensors
        y.backward(dy)

        # Collect input gradients
        input_grad = input_tensor.grad.detach().clone()
        input_t_grad = input_t.grad.detach().clone() if asymmetric else None

        return (
            input_tensor.detach().clone(),
            input_t.detach().clone() if asymmetric else None,
            y.detach().clone(),
            input_grad,
            input_t_grad,
            dy.detach().clone(),
        )


def compute_dtensor_native_outer_op(
    input_global: torch.Tensor,
    input_t_global: torch.Tensor | None,
    dy_global: torch.Tensor | None,
    device_mesh: DeviceMesh,
    input_placements: tuple[Placement, ...],
    output_placements: tuple[Placement, ...],
    axis: int,
    op: OuterOp,
) -> tuple[list[DTensor], DTensor]:
    """Compute DTensor native operations for comparison."""
    # Create DTensor native inputs
    input_dtensor = distribute_tensor(input_global.detach().clone(), device_mesh, input_placements).requires_grad_(
        op != OuterOp.BITAND
    )

    if input_t_global is not None:
        input_t_dtensor = distribute_tensor(
            input_t_global.detach().clone(), device_mesh, input_placements
        ).requires_grad_(op != OuterOp.BITAND)
    else:
        input_t_dtensor = input_dtensor  # Symmetric case

    # Forward pass with native DTensor operations (manual outer operation)
    input_expanded = input_dtensor.unsqueeze(axis + 1)
    input_t_expanded = input_t_dtensor.unsqueeze(axis + 1)

    placements_y = (Shard(0), Shard(1), Shard(2))

    # it's necessary to redistribute the input_t_expanded to the placements_y
    # in order to avoid a runtime error raise from the native DTensor backward:
    # redistribute S(1) -> R in backward is not supported
    input_t_expanded_transposed = input_t_expanded.transpose(axis, axis + 1).redistribute(placements=placements_y)

    if op == OuterOp.SUM:
        y_dtensor_result_native = input_expanded + input_t_expanded_transposed
    elif op == OuterOp.SUBTRACT:
        y_dtensor_result_native = input_expanded - input_t_expanded_transposed
    elif op == OuterOp.EQUAL:
        y_dtensor_result_native = input_expanded == input_t_expanded_transposed
    elif op == OuterOp.BITAND:
        y_dtensor_result_native = input_expanded & input_t_expanded_transposed
    elif op == OuterOp.PROD:
        y_dtensor_result_native = input_expanded * input_t_expanded_transposed
    elif op == OuterOp.CDIST:
        y_dtensor_result_with_zeros = torch.sum((input_expanded - input_t_expanded_transposed) ** 2, dim=-1)
        # to avoid the diagonal zeros causing backward pass nan
        y_dtensor_result_native = (
            y_dtensor_result_with_zeros + torch.finfo(y_dtensor_result_with_zeros.dtype).tiny
        ).sqrt()

    if op == OuterOp.EQUAL or op == OuterOp.BITAND or dy_global is None:
        # No backward pass for EQUAL/BITAND operation
        return [], y_dtensor_result_native

    # Backward pass with native DTensor op
    dy_dtensor = distribute_tensor(dy_global.detach().clone(), device_mesh, output_placements)
    y_dtensor_result_native.backward(dy_dtensor)

    # redistribute here to avoid comparing Partial() placements with Replicate() placements, the
    # latter of which is in the replicate_to_shard_outer_op due to the explicit all_reduce op
    inputs_grad_dtensor = [input_dtensor.grad.redistribute(placements=input_placements)]
    if input_t_global is not None:
        inputs_grad_dtensor.append(input_t_dtensor.grad.redistribute(placements=input_placements))

    return inputs_grad_dtensor, y_dtensor_result_native


def compute_replicate_to_shard_outer_op_with_validation(
    input_global: torch.Tensor,
    input_t_global: torch.Tensor | None,
    dy_global: torch.Tensor | None,
    device_mesh: DeviceMesh,
    input_placements: tuple[Placement, ...],
    output_placements: tuple[Placement, ...],
    transpose_comm: TransposeComm,
    axis: int,
    op: OuterOp,
    label_test_case: str,
) -> tuple[DTensor, DTensor, DTensor | None, DTensor | None]:
    """
    Compute replicate_to_shard_outer_op forward and backward pass with input validation checks.

    Returns:
        y_dtensor_result: Forward pass result
        input_dtensor: Input tensor with computed gradient
        input_t_dtensor: Second input tensor with computed gradient (if asymmetric)
        dy_dtensor: Distributed upstream gradient
    """
    # Create DTensor inputs
    input_dtensor = distribute_tensor(input_global.detach().clone(), device_mesh, input_placements).requires_grad_(
        op != OuterOp.BITAND
    )

    if input_t_global is not None:
        input_t_dtensor = distribute_tensor(
            input_t_global.detach().clone(), device_mesh, input_placements
        ).requires_grad_(op != OuterOp.BITAND)
    else:
        input_t_dtensor = None

    input_dtensor_copy = input_dtensor.detach().clone().requires_grad_(op != OuterOp.BITAND)
    input_t_dtensor_copy = (
        input_t_dtensor.detach().clone().requires_grad_(op != OuterOp.BITAND) if input_t_dtensor is not None else None
    )

    # Compute on distributed tensors using replicate_to_shard_outer_op
    y_dtensor_result = replicate_to_shard_outer_op(input_dtensor, op, axis, transpose_comm, input_t_dtensor)

    # Verify no change to the forward inputs
    assert_tensors_identical(
        input_dtensor.to_local(), input_dtensor_copy.to_local(), check_grad=False, check_grad_fn=False
    )
    if input_t_dtensor is not None:
        assert_tensors_identical(
            input_t_dtensor.to_local(), input_t_dtensor_copy.to_local(), check_grad=False, check_grad_fn=False
        )

    # Verify output placements
    assert y_dtensor_result.placements == output_placements, f"{label_test_case} output placements mismatch"

    if op == OuterOp.EQUAL or op == OuterOp.BITAND or dy_global is None:
        # No backward pass for EQUAL/BITAND operation
        return y_dtensor_result, input_dtensor, input_t_dtensor, None

    # Distribute the upstream adjoint for backward pass
    dy_dtensor = distribute_tensor(dy_global.detach().clone(), device_mesh, output_placements)

    # Perform backward pass
    dy_dtensor_copy = dy_dtensor.detach().clone()
    y_dtensor_result.backward(dy_dtensor)

    # Verify no change to the backward input
    assert_tensors_identical(dy_dtensor.to_local(), dy_dtensor_copy.to_local(), check_grad=False, check_grad_fn=False)

    # Verify input gradient placements are consistent with input placements
    assert (
        input_dtensor.grad.placements == input_placements
    ), f"{label_test_case} inconsistent input gradient placements with input placements"

    if input_t_dtensor is not None:
        assert (
            input_t_dtensor.grad.placements == input_placements
        ), f"{label_test_case} inconsistent input_t gradient placements with input placements"

    return y_dtensor_result, input_dtensor, input_t_dtensor, dy_dtensor


def parallel_assert_replicate_to_shard_outer_op(
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

    # Each rank uses the same seed to generate the same input tensors
    seed_by_rank(0, seed=42)

    size_cp = len(manager.group_ranks["cp"])
    size_ring = isqrt(size_cp)
    if size_ring * size_ring != size_cp:
        raise ValueError(f"cp group size {size_cp} is not a square int")

    # Set test parameters
    # Input placements: (Shard(0), Shard(1), Replicate())
    # Output placements: (Shard(0), Shard(1), Shard(2))
    batch_size = 2 * len(manager.group_ranks["dp"])
    seq_len_dim0 = size_ring * 4  # Sharded along dim 0
    embed_dim = 8

    shape_input = (batch_size, seq_len_dim0, embed_dim)
    input_placements = (Shard(0), Shard(1), Replicate())
    output_placements = (Shard(0), Shard(1), Shard(2))

    # Create transpose communication
    layout_map = manager.layout_subgroups["cp"]
    transpose_comm = TransposeComm(manager.group["cp"], layout_map)

    # Test all OuterOp cases and both symmetric/asymmetric
    for op in [OuterOp.SUM, OuterOp.SUBTRACT, OuterOp.EQUAL, OuterOp.BITAND, OuterOp.PROD]:
        for asymmetric in [False, True]:
            axis = 1
            label_test_case = (
                f"op={op}, asymmetric={asymmetric}, axis={axis}, "
                f"input_placements={input_placements}, output_placements={output_placements}\n"
            )

            # Compute global expectations
            (
                input_global,
                input_t_global,
                y_expected_global,
                input_grad_expected_global,
                input_t_grad_expected_global,
                dy_global,
            ) = compute_global_expectation(shape_input, axis, op, manager.device, asymmetric)

            # Use DTensor native op as an alternative reference (always call, even for EQUAL)
            (
                inputs_grad_dtensor_native,
                y_dtensor_result_native,
            ) = compute_dtensor_native_outer_op(
                input_global,
                input_t_global,
                dy_global,  # This will be None for EQUAL operations
                manager.device_mesh_subgroups,
                input_placements,
                output_placements,
                axis,
                op,
            )

            # Compute replicate_to_shard_outer_op forward and backward with validation
            y_dtensor_result, input_dtensor, input_t_dtensor, dy_dtensor = (
                compute_replicate_to_shard_outer_op_with_validation(
                    input_global,
                    input_t_global,
                    dy_global,
                    manager.device_mesh_subgroups,
                    input_placements,
                    output_placements,
                    transpose_comm,
                    axis,
                    op,
                    label_test_case,
                )
            )

            # ===================================================================
            # BLOCK 1: Check against DTensor native reference
            # ===================================================================

            # Check metadata against DTensor native (for all operations)
            assert (
                y_dtensor_result.placements == y_dtensor_result_native.placements
            ), f"{label_test_case} placements mismatch"
            assert y_dtensor_result.shape == y_dtensor_result_native.shape, f"{label_test_case} shape mismatch"
            assert y_dtensor_result.stride() == y_dtensor_result_native.stride(), f"{label_test_case} stride mismatch"

            # Compare forward result with native DTensor op (for all operations)
            torch.testing.assert_close(
                y_dtensor_result.to_local(),
                y_dtensor_result_native.to_local(),
                msg=lambda m: f"{label_test_case} {m}",
            )

            # Compare global tensors between replicate_to_shard_outer_op and native DTensor results
            y_result_global = y_dtensor_result.full_tensor()
            y_result_global_native = y_dtensor_result_native.full_tensor()

            torch.testing.assert_close(
                y_result_global,
                y_result_global_native,
                msg=lambda m: f"{label_test_case} output vs native: {m}",
            )

            # Only check gradients for non-EQUAL and non-BITAND operations
            if op != OuterOp.EQUAL and op != OuterOp.BITAND:
                # Assert input gradients' metadata and values against DTensor native
                # Input gradient comparison
                assert (
                    input_dtensor.grad.placements == inputs_grad_dtensor_native[0].placements
                ), f"{label_test_case} input gradient placements mismatch"
                assert (
                    input_dtensor.grad.shape == inputs_grad_dtensor_native[0].shape
                ), f"{label_test_case} input gradient shape mismatch"
                assert (
                    input_dtensor.grad.stride() == inputs_grad_dtensor_native[0].stride()
                ), f"{label_test_case} input gradient stride mismatch"

                torch.testing.assert_close(
                    input_dtensor.grad.to_local(),
                    inputs_grad_dtensor_native[0].to_local(),
                    msg=lambda m: f"{label_test_case} input gradient mismatch: {m}",
                )

                torch.testing.assert_close(
                    input_dtensor.grad.full_tensor(),
                    inputs_grad_dtensor_native[0].full_tensor(),
                    msg=lambda m: f"{label_test_case} input gradient mismatch: {m}",
                )

                # Input_t gradient comparison (asymmetric case)
                if asymmetric and len(inputs_grad_dtensor_native) > 1:
                    assert (
                        input_t_dtensor.grad.placements == inputs_grad_dtensor_native[1].placements
                    ), f"{label_test_case} input_t gradient placements mismatch"
                    assert (
                        input_t_dtensor.grad.shape == inputs_grad_dtensor_native[1].shape
                    ), f"{label_test_case} input_t gradient shape mismatch"
                    assert (
                        input_t_dtensor.grad.stride() == inputs_grad_dtensor_native[1].stride()
                    ), f"{label_test_case} input_t gradient stride mismatch"

                    torch.testing.assert_close(
                        input_t_dtensor.grad.to_local(),
                        inputs_grad_dtensor_native[1].to_local(),
                        msg=lambda m: f"{label_test_case} input_t gradient mismatch: {m}",
                    )

                    torch.testing.assert_close(
                        input_t_dtensor.grad.full_tensor(),
                        inputs_grad_dtensor_native[1].full_tensor(),
                        msg=lambda m: f"{label_test_case} input_t gradient mismatch: {m}",
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
                msg=lambda m: f"{label_test_case} forward result: {m}",
            )

            # Compare forward result with global expectation
            y_result_global = y_dtensor_result.full_tensor()
            torch.testing.assert_close(
                y_result_global,
                y_expected_global,
                msg=lambda m: f"{label_test_case} forward result vs global expectation: {m}",
            )

            if op != OuterOp.EQUAL and op != OuterOp.BITAND:
                # Check input gradient
                input_grad_expected_dtensor = distribute_tensor(
                    input_grad_expected_global, manager.device_mesh_subgroups, input_placements
                )

                torch.testing.assert_close(
                    input_dtensor.grad.to_local(),
                    input_grad_expected_dtensor.to_local(),
                    msg=lambda m: f"{label_test_case} input gradient vs global expectation: {m}",
                )

                torch.testing.assert_close(
                    input_dtensor.grad.full_tensor(),
                    input_grad_expected_global,
                    msg=lambda m: f"{label_test_case} input gradient vs global expectation: {m}",
                )

                # Check input_t gradient (asymmetric case)
                if asymmetric and input_t_grad_expected_global is not None:
                    input_t_grad_expected_dtensor = distribute_tensor(
                        input_t_grad_expected_global, manager.device_mesh_subgroups, input_placements
                    )

                    torch.testing.assert_close(
                        input_t_dtensor.grad.to_local(),
                        input_t_grad_expected_dtensor.to_local(),
                        msg=lambda m: f"{label_test_case} input_t gradient vs global expectation: {m}",
                    )

                    torch.testing.assert_close(
                        input_t_dtensor.grad.full_tensor(),
                        input_t_grad_expected_global,
                        msg=lambda m: f"{label_test_case} input_t gradient vs global expectation: {m}",
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
def test_replicate_to_shard_outer_op(setup_env):
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    spawn_multiprocessing(
        parallel_assert_replicate_to_shard_outer_op,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
    )
