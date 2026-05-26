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


import pytest
import torch
from torch.distributed.tensor import Replicate, Shard, distribute_tensor

from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.layers.elementwise_op import (
    ElementwiseOp,
    elementwise_op,
    scalar_tensor_op,
    single_tensor_op,
)
from boltz.testing.utils import assert_tensors_identical, spawn_multiprocessing


def serial_elementwise_op(a: torch.Tensor, b: torch.Tensor, op: ElementwiseOp) -> torch.Tensor:
    """Serial implementation of elementwise operation for comparison."""
    if op == ElementwiseOp.SUM:
        return a + b
    elif op == ElementwiseOp.SUB:
        return a - b
    elif op == ElementwiseOp.PROD:
        return a * b
    elif op == ElementwiseOp.DIV:
        return a / b
    elif op == ElementwiseOp.EQUAL:
        return a & b
    elif op == ElementwiseOp.BITAND:
        return a & b
    else:
        raise ValueError(f"Unsupported operation: {op}")


def parallel_assert_elementwise_op(
    rank,
    grid_group_sizes,
    device_type,
    backend,
    env_per_rank,
    placements,
    op,
    input_a_global_host,
    input_b_global_host,
    output_expected_global_host,
    d_output_expected_global_host,
    d_input_a_expected_global_host,
    d_input_b_expected_global_host,
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

    # Distribute input tensors
    input_a_dtensor = distribute_tensor(
        input_a_global_host.to(manager.device), device_mesh=manager.device_mesh_subgroups, placements=placements
    )
    input_b_dtensor = distribute_tensor(
        input_b_global_host.to(manager.device), device_mesh=manager.device_mesh_subgroups, placements=placements
    )
    if op != ElementwiseOp.EQUAL and op != ElementwiseOp.BITAND:
        input_a_dtensor.requires_grad_(True)
        input_b_dtensor.requires_grad_(True)

    # Distribute expected outputs
    output_expected_dtensor = distribute_tensor(
        output_expected_global_host.to(manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements,
        src_data_rank=None,
    )
    if op == ElementwiseOp.EQUAL or op == ElementwiseOp.BITAND:
        d_output_expected_dtensor = None
        d_input_a_expected_dtensor = None
        d_input_b_expected_dtensor = None
    else:
        d_output_expected_dtensor = distribute_tensor(
            d_output_expected_global_host.to(manager.device),
            device_mesh=manager.device_mesh_subgroups,
            placements=placements,
        )
        d_input_a_expected_dtensor = distribute_tensor(
            d_input_a_expected_global_host.to(manager.device),
            device_mesh=manager.device_mesh_subgroups,
            placements=placements,
            src_data_rank=None,
        )
        d_input_b_expected_dtensor = distribute_tensor(
            d_input_b_expected_global_host.to(manager.device),
            device_mesh=manager.device_mesh_subgroups,
            placements=placements,
            src_data_rank=None,
        )

    # Create copies to verify inputs aren't modified
    input_a_dtensor_copy = input_a_dtensor.detach().clone()
    input_b_dtensor_copy = input_b_dtensor.detach().clone()
    if op != ElementwiseOp.EQUAL and op != ElementwiseOp.BITAND:
        input_a_dtensor_copy.requires_grad_(True)
        input_b_dtensor_copy.requires_grad_(True)

    # Forward pass
    output_dtensor_result = elementwise_op(input_a_dtensor, input_b_dtensor, op)

    # Verify inputs weren't modified
    assert_tensors_identical(
        input_a_dtensor_copy.to_local(), input_a_dtensor.to_local(), check_grad=False, check_grad_fn=False
    )
    assert_tensors_identical(
        input_b_dtensor_copy.to_local(), input_b_dtensor.to_local(), check_grad=False, check_grad_fn=False
    )

    # Test forward pass results
    torch.testing.assert_close(output_dtensor_result.to_local(), output_expected_dtensor.to_local())

    # Test shape and stride consistency
    assert (
        output_dtensor_result.shape == output_expected_dtensor.shape
    ), f"Output shape mismatch: {output_dtensor_result.shape} != {output_expected_dtensor.shape}"
    assert (
        output_dtensor_result.stride() == output_expected_dtensor.stride()
    ), f"Output stride mismatch: {output_dtensor_result.stride()} != {output_expected_dtensor.stride()}"

    # Backward pass
    if op == ElementwiseOp.EQUAL or op == ElementwiseOp.BITAND:
        with pytest.raises(RuntimeError, match="tensors does not require grad and does not have a grad_fn"):
            output_dtensor_result.backward(d_output_expected_dtensor)
    else:
        d_output_expected_dtensor_copy = d_output_expected_dtensor.detach().clone()
        output_dtensor_result.backward(d_output_expected_dtensor)

        # Verify upstream gradient wasn't modified
        assert_tensors_identical(d_output_expected_dtensor_copy.to_local(), d_output_expected_dtensor.to_local())

    # Test input gradients
    if op == ElementwiseOp.EQUAL or op == ElementwiseOp.BITAND:
        assert (
            d_output_expected_dtensor is None
            and d_input_a_expected_dtensor is None
            and d_input_b_expected_dtensor is None
        )
        assert input_a_dtensor.grad is None and input_b_dtensor.grad is None
    else:
        torch.testing.assert_close(input_a_dtensor.grad.to_local(), d_input_a_expected_dtensor.to_local())
        torch.testing.assert_close(input_b_dtensor.grad.to_local(), d_input_b_expected_dtensor.to_local())

        # Test gradient shape and stride consistency
        assert (
            input_a_dtensor.grad.shape == d_input_a_expected_dtensor.shape
        ), f"Input A gradient shape mismatch: {input_a_dtensor.grad.shape} != {d_input_a_expected_dtensor.shape}"
        assert (
            input_a_dtensor.grad.stride() == d_input_a_expected_dtensor.stride()
        ), f"Input A gradient stride mismatch: {input_a_dtensor.grad.stride()} != {d_input_a_expected_dtensor.stride()}"
        assert (
            input_b_dtensor.grad.shape == d_input_b_expected_dtensor.shape
        ), f"Input B gradient shape mismatch: {input_b_dtensor.grad.shape} != {d_input_b_expected_dtensor.shape}"
        assert (
            input_b_dtensor.grad.stride() == d_input_b_expected_dtensor.stride()
        ), f"Input B gradient stride mismatch: {input_b_dtensor.grad.stride()} != {d_input_b_expected_dtensor.stride()}"

    # Test full tensor gathering - verify distributed results match serial results
    output_global_result_host = output_dtensor_result.full_tensor().cpu()
    torch.testing.assert_close(output_global_result_host, output_expected_global_host)
    if op != ElementwiseOp.EQUAL and op != ElementwiseOp.BITAND:
        d_input_a_global_result_host = input_a_dtensor.grad.full_tensor().cpu()
        d_input_b_global_result_host = input_b_dtensor.grad.full_tensor().cpu()
        torch.testing.assert_close(d_input_a_global_result_host, d_input_a_expected_global_host)
        torch.testing.assert_close(d_input_b_global_result_host, d_input_b_expected_global_host)

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
    "placements", [(Shard(0), Shard(1), Shard(2)), (Shard(0), Shard(1), Replicate())], ids=["shard", "replicate"]
)
@pytest.mark.parametrize(
    "op",
    [
        ElementwiseOp.SUM,
        ElementwiseOp.SUB,
        ElementwiseOp.PROD,
        ElementwiseOp.DIV,
        ElementwiseOp.EQUAL,
        ElementwiseOp.BITAND,
    ],
    ids=["sum", "sub", "prod", "div", "equal", "bitand"],
)
def test_elementwise_op_parallel(setup_env, placements, op):
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

    seed = 42
    rng = torch.Generator(device=device_type)
    rng.manual_seed(seed)

    # Create input tensors with proper shapes
    if op == ElementwiseOp.EQUAL or op == ElementwiseOp.BITAND:
        input_a_global = torch.randint(0, 2, (B, N, N, D), requires_grad=False, device=device_type, generator=rng)
        input_b_global = torch.randint(0, 2, (B, N, N, D), requires_grad=False, device=device_type, generator=rng)
    else:
        input_a_global = torch.empty((B, N, N, D), requires_grad=True, device=device_type)
        input_b_global = torch.empty((B, N, N, D), requires_grad=True, device=device_type)
        with torch.no_grad():
            input_a_global.uniform_(-1000, 1000, generator=rng)
            input_b_global.uniform_(-1000, 1000, generator=rng)

    # Run serial forward pass
    input_a_global_host = input_a_global.detach().clone().cpu()
    input_b_global_host = input_b_global.detach().clone().cpu()
    output_expected_global = serial_elementwise_op(input_a_global, input_b_global, op)
    output_expected_global_host = output_expected_global.detach().clone().cpu()

    # Create upstream gradient and run backward pass
    if op == ElementwiseOp.EQUAL or op == ElementwiseOp.BITAND:
        with pytest.raises(RuntimeError, match="tensors does not require grad and does not have a grad_fn"):
            output_expected_global.backward(torch.empty_like(output_expected_global))

        d_output_expected_global = None
        d_output_expected_global_host = None
    else:
        d_output_expected_global = torch.rand_like(output_expected_global)
        d_output_expected_global_host = d_output_expected_global.detach().clone().cpu()
        output_expected_global.backward(d_output_expected_global)

    # Launch parallel test across all processes
    spawn_multiprocessing(
        parallel_assert_elementwise_op,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        placements,
        op,
        input_a_global_host,
        input_b_global_host,
        output_expected_global_host,
        d_output_expected_global_host,
        input_a_global.grad.detach().clone().cpu()
        if op != ElementwiseOp.EQUAL and op != ElementwiseOp.BITAND
        else None,
        input_b_global.grad.detach().clone().cpu()
        if op != ElementwiseOp.EQUAL and op != ElementwiseOp.BITAND
        else None,
    )


def serial_single_tensor_op(x: torch.Tensor, op: ElementwiseOp) -> torch.Tensor:
    """Serial implementation of single tensor operation for comparison."""
    if op == ElementwiseOp.COS:
        return torch.cos(x)
    elif op == ElementwiseOp.RELU:
        return torch.relu(x)
    elif op == ElementwiseOp.ROUND:
        return torch.round(x)
    elif op == ElementwiseOp.LOG:
        return torch.log(x)
    elif op == ElementwiseOp.EXP:
        return torch.exp(x)
    elif op == ElementwiseOp.ABS:
        return torch.abs(x)
    elif op == ElementwiseOp.SIGMOID:
        return torch.sigmoid(x)
    else:
        raise ValueError(f"Unsupported single tensor operation: {op}")


def parallel_assert_single_tensor_op(
    rank,
    grid_group_sizes,
    device_type,
    backend,
    env_per_rank,
    placements,
    op,
    input_x_global_host,
    output_expected_global_host,
    d_output_expected_global_host,
    d_input_x_expected_global_host,
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
    output_dtensor_result = single_tensor_op(input_x_dtensor, op)

    # Verify input wasn't modified
    assert_tensors_identical(
        input_x_dtensor_copy.to_local(), input_x_dtensor.to_local(), check_grad=False, check_grad_fn=False
    )

    # Test forward pass results
    torch.testing.assert_close(output_dtensor_result.to_local(), output_expected_dtensor.to_local())

    # Backward pass
    d_output_expected_dtensor_copy = d_output_expected_dtensor.detach().clone()
    if op == ElementwiseOp.ROUND:
        with pytest.raises(RuntimeError, match="tensors does not require grad and does not have a grad_fn"):
            output_dtensor_result.backward(d_output_expected_dtensor)
    else:
        output_dtensor_result.backward(d_output_expected_dtensor)

    # Verify upstream gradient wasn't modified
    assert_tensors_identical(d_output_expected_dtensor_copy.to_local(), d_output_expected_dtensor.to_local())

    # Test input gradient
    if op == ElementwiseOp.ROUND:
        assert input_x_dtensor.grad is None
    else:
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
    "op",
    [
        ElementwiseOp.COS,
        ElementwiseOp.RELU,
        ElementwiseOp.ROUND,
        ElementwiseOp.LOG,
        ElementwiseOp.EXP,
        ElementwiseOp.ABS,
        ElementwiseOp.SIGMOID,
    ],
    ids=["cos", "relu", "round", "log", "exp", "abs", "sigmoid"],
)
def test_single_tensor_op_parallel(setup_env, placements, op):
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

    seed = 42
    rng = torch.Generator(device=device_type)
    rng.manual_seed(seed)

    # Create input tensor with proper shape
    input_x_global = torch.empty((B, N, N, D), requires_grad=True, device=device_type)
    with torch.no_grad():
        if op == ElementwiseOp.LOG:
            # For LOG operation, use positive values only
            # On certain GPU architectures, we need to limit the range of the input
            # to limit the numerical errors
            input_x_global.uniform_(0.1, 50, generator=rng)
        elif op == ElementwiseOp.EXP:
            # For EXP operation, use a reasonable range to avoid overflow
            input_x_global.uniform_(-10, 10, generator=rng)
        else:
            # For other operations, test a wide range of values including negative values and zeros for ReLU testing
            input_x_global.uniform_(-1000, 1000, generator=rng)
            # Ensure we have some zeros for ReLU boundary testing
            input_x_global.view(-1)[::17] = 0.0  # Set every 17th element to zero

    # Run serial forward pass
    input_x_global_host = input_x_global.detach().clone().cpu()
    output_expected_global = serial_single_tensor_op(input_x_global, op)
    output_expected_global_host = output_expected_global.detach().clone().cpu()

    # Create upstream gradient and run backward pass
    d_output_expected_global = torch.rand_like(output_expected_global)
    d_output_expected_global_host = d_output_expected_global.detach().clone().cpu()
    output_expected_global.backward(d_output_expected_global)

    # Launch parallel test across all processes
    spawn_multiprocessing(
        parallel_assert_single_tensor_op,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        placements,
        op,
        input_x_global_host,
        output_expected_global_host,
        d_output_expected_global_host,
        input_x_global.grad.detach().clone().cpu(),
    )


def serial_scalar_tensor_op(scalar: float, tensor: torch.Tensor, op: ElementwiseOp) -> torch.Tensor:
    """Serial implementation of scalar-tensor operation for comparison."""
    if op == ElementwiseOp.SUM:
        return scalar + tensor
    elif op == ElementwiseOp.SUB:
        return scalar - tensor
    elif op == ElementwiseOp.PROD:
        return scalar * tensor
    elif op == ElementwiseOp.DIV:
        return scalar / tensor
    elif op == ElementwiseOp.GT:
        return scalar > tensor
    elif op == ElementwiseOp.LT:
        return scalar < tensor
    elif op == ElementwiseOp.EQUAL:
        return scalar == tensor
    elif op == ElementwiseOp.POW:
        return torch.pow(tensor, scalar)
    elif op == ElementwiseOp.MAX:
        return torch.clamp(tensor, min=scalar)
    else:
        raise ValueError(f"Unsupported scalar-tensor operation: {op}")


def parallel_assert_scalar_tensor_op(
    rank,
    grid_group_sizes,
    device_type,
    backend,
    env_per_rank,
    placements,
    op,
    scalar,
    input_tensor_global_host,
    output_expected_global_host,
    d_output_expected_global_host,
    d_input_tensor_expected_global_host,
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
    input_tensor_dtensor = distribute_tensor(
        input_tensor_global_host.to(manager.device), device_mesh=manager.device_mesh_subgroups, placements=placements
    )
    input_tensor_dtensor.requires_grad_(True)

    # Distribute expected outputs
    if op == ElementwiseOp.GT or op == ElementwiseOp.LT or op == ElementwiseOp.EQUAL:
        d_output_expected_dtensor = None
        d_input_tensor_expected_dtensor = None
    else:
        d_output_expected_dtensor = distribute_tensor(
            d_output_expected_global_host.to(manager.device),
            device_mesh=manager.device_mesh_subgroups,
            placements=placements,
        )
        d_input_tensor_expected_dtensor = distribute_tensor(
            d_input_tensor_expected_global_host.to(manager.device),
            device_mesh=manager.device_mesh_subgroups,
            placements=placements,
            src_data_rank=None,
        )
    output_expected_dtensor = distribute_tensor(
        output_expected_global_host.to(manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements,
        src_data_rank=None,
    )

    # Create copy to verify input isn't modified
    input_tensor_dtensor_copy = input_tensor_dtensor.detach().clone().requires_grad_(input_tensor_dtensor.requires_grad)

    # Forward pass
    output_dtensor_result = scalar_tensor_op(scalar, input_tensor_dtensor, op)

    # Verify input wasn't modified
    assert_tensors_identical(
        input_tensor_dtensor_copy.to_local(), input_tensor_dtensor.to_local(), check_grad=False, check_grad_fn=False
    )

    # Test forward pass results
    torch.testing.assert_close(output_dtensor_result.to_local(), output_expected_dtensor.to_local())

    # Backward pass
    if op == ElementwiseOp.GT or op == ElementwiseOp.LT or op == ElementwiseOp.EQUAL:
        with pytest.raises(RuntimeError, match="tensors does not require grad and does not have a grad_fn"):
            output_dtensor_result.backward(d_output_expected_dtensor)
    else:
        d_output_expected_dtensor_copy = d_output_expected_dtensor.detach().clone()
        output_dtensor_result.backward(d_output_expected_dtensor)

        # Verify upstream gradient wasn't modified
        assert_tensors_identical(d_output_expected_dtensor_copy.to_local(), d_output_expected_dtensor.to_local())

        # Test input gradient for tensor
        torch.testing.assert_close(input_tensor_dtensor.grad.to_local(), d_input_tensor_expected_dtensor.to_local())

    # Test full tensor gathering - verify distributed results match serial results
    output_global_result_host = output_dtensor_result.full_tensor().cpu()

    # Verify full tensors match expected results
    torch.testing.assert_close(output_global_result_host, output_expected_global_host)

    if op != ElementwiseOp.GT and op != ElementwiseOp.LT and op != ElementwiseOp.EQUAL:
        d_input_tensor_global_result_host = input_tensor_dtensor.grad.full_tensor().cpu()
        torch.testing.assert_close(d_input_tensor_global_result_host, d_input_tensor_expected_global_host)
    else:
        # For GT, LT, and EQUAL operations, verify that no gradients were computed
        assert input_tensor_dtensor.grad is None

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env",
    [
        ((1, (2, 2)), True, "cuda", "ENV"),
        ((2, (2, 2)), True, "cuda", "ENV"),
        ((1, (3, 3)), True, "cpu", "ENV"),
    ],
    indirect=("setup_env",),
    ids=lambda x: f"dp={x[0][0]}, cp={x[0][1]}, specify_method={x[1]}, device_type={x[2]}, method_init={x[3]}",
)
@pytest.mark.parametrize(
    "placements", [(Shard(0), Shard(1), Shard(2)), (Shard(0), Shard(1), Replicate())], ids=["shard", "replicate"]
)
@pytest.mark.parametrize(
    "op",
    [
        ElementwiseOp.SUM,
        ElementwiseOp.SUB,
        ElementwiseOp.PROD,
        ElementwiseOp.DIV,
        ElementwiseOp.GT,
        ElementwiseOp.LT,
        ElementwiseOp.EQUAL,
        ElementwiseOp.POW,
        ElementwiseOp.MAX,
    ],
    ids=["sum", "sub", "prod", "div", "gt", "lt", "equal", "pow", "max"],
)
def test_scalar_tensor_op_parallel(setup_env, placements, op):
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

    seed = 42
    rng = torch.Generator(device=device_type)
    rng.manual_seed(seed)

    # Create input tensor
    scalar_value = torch.rand(1, device=device_type, generator=rng).item() * 20.0 - 10.0  # [-10, 10]
    input_tensor_global = torch.empty(
        (B, N, N, D),
        requires_grad=(op != ElementwiseOp.GT and op != ElementwiseOp.LT and op != ElementwiseOp.EQUAL),
        device=device_type,
    )
    with torch.no_grad():
        input_tensor_global.uniform_(-10, 10, generator=rng)
        if op == ElementwiseOp.POW:
            input_tensor_global.abs_()  # fractional power of negative number is complex

    # Run serial forward pass
    input_tensor_global_host = input_tensor_global.detach().clone().cpu()
    output_expected_global = serial_scalar_tensor_op(scalar_value, input_tensor_global, op)
    output_expected_global_host = output_expected_global.detach().clone().cpu()

    # Create upstream gradient and run backward pass
    if op == ElementwiseOp.GT or op == ElementwiseOp.LT or op == ElementwiseOp.EQUAL:
        d_output_expected_global = None
        d_output_expected_global_host = None
    else:
        d_output_expected_global = torch.rand_like(output_expected_global)
        d_output_expected_global_host = d_output_expected_global.detach().clone().cpu()
        output_expected_global.backward(d_output_expected_global)

    # Launch parallel test across all processes
    spawn_multiprocessing(
        parallel_assert_scalar_tensor_op,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        placements,
        op,
        scalar_value,
        input_tensor_global_host,
        output_expected_global_host,
        d_output_expected_global_host,
        input_tensor_global.grad.detach().clone().cpu()
        if op != ElementwiseOp.GT and op != ElementwiseOp.LT and op != ElementwiseOp.EQUAL
        else None,
    )
