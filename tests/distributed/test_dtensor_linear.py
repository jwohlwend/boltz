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
from torch.distributed.tensor import Replicate, Shard, distribute_module, distribute_tensor

from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.layers.linear import LinearParamsReplicated
from boltz.testing.utils import assert_all_identical, seed_by_rank


def compute_global_expectation(batch_size, seq_len, in_features, out_features, device):
    """Compute expected results using global tensors.

    Args:
        batch_size: Batch size
        seq_len: Sequence length
        in_features: Input feature dimension
        out_features: Output feature dimension
        device: Device to place tensors on

    Returns:
        tuple: (x_global, weight_global, bias_global, y_global_expected,
                x_global_grad, weight_global_grad, bias_global_grad, dy_global)
    """
    # Create global tensors with deterministic values for reproducibility
    x_global = torch.randn(batch_size, seq_len, in_features, device=device, requires_grad=True)
    linear = torch.nn.Linear(in_features, out_features, device=device)
    state_dict_linear = linear.state_dict()

    # Clone inputs for distribution
    x_global_clone = x_global.detach().clone()

    # Compute on global tensors using standard linear operation
    y_global_expected = linear(x_global)

    # Create gradients for backward pass
    dy_global = torch.rand_like(y_global_expected)

    # Backward pass on global tensors
    y_global_expected.backward(dy_global)

    return (
        x_global_clone,
        state_dict_linear,
        y_global_expected.detach().clone(),
        x_global.grad.detach().clone(),
        linear.weight.grad.detach().clone(),
        linear.bias.grad.detach().clone(),
        dy_global.detach().clone(),
    )


def parallel_assert_dtensor_linear(
    rank: int,
    grid_group_sizes: Dict[str, int],
    device_type: str,
    backend: str,
    compare_to_native_dtensor: bool,
    env_map: Optional[Dict[str, str]] = None,
):
    """Test distributed linear operation in a parallel environment.

    This test validates that the LinearParamsReplicatedImpl produces identical results to
    standard nn.functional.linear operations with global tensors. It verifies:

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

    try:
        size_cp = len(manager.group_ranks["cp"])
        size_ring = isqrt(size_cp)
        if size_ring * size_ring != size_cp:
            raise ValueError(f"cp group size {size_cp} is not a square int")

        # Set test parameters
        batch_size = 2
        seq_len_per_rank = 4
        seq_len_global = size_ring * seq_len_per_rank
        in_features = 3
        out_features = 5

        # Set random seed based on rank for reproducibility
        seed_by_rank(0)

        # Compute global expectations
        (
            x_global,
            state_dict_linear,
            y_global_expected,
            x_global_grad,
            weight_global_grad,
            bias_global_grad,
            dy_global,
        ) = compute_global_expectation(batch_size, seq_len_global, in_features, out_features, manager.device)

        # Create distributed tensors
        # Shard the sequence dimension (dim=1) for input tensor
        # this emulates the sharded single representation in the Boltz model
        input_placements = [Shard(dim=0), Shard(dim=1), Replicate()]
        x_dtensor = distribute_tensor(x_global, manager.device_mesh_subgroups, input_placements)
        x_dtensor.requires_grad = True

        linear_local = torch.nn.Linear(in_features, out_features, device=manager.device)
        linear_local.load_state_dict(state_dict_linear)
        layer = LinearParamsReplicated(linear_local, manager.device_mesh_subgroups)

        linear_local_copy = torch.nn.Linear(in_features, out_features, device=manager.device)
        linear_local_copy.load_state_dict(state_dict_linear)
        layer_dtensor_native = distribute_module(
            linear_local_copy,
            manager.device_mesh_subgroups,
            output_fn=lambda module, outputs, device_mesh: outputs.redistribute(device_mesh, input_placements),
        )
        x_dtensor_native = x_dtensor.detach().clone().requires_grad_(True)

        # Compute on distributed tensors using LinearParamsReplicatedImpl
        y_dtensor_result = layer(x_dtensor)

        # Distribute the upstream adjoint for backward pass
        dy_dtensor = distribute_tensor(dy_global, manager.device_mesh_subgroups, input_placements)

        # Perform backward pass
        y_dtensor_result.backward(dy_dtensor)

        y_dtensor_native = None
        if compare_to_native_dtensor:
            y_dtensor_native = layer_dtensor_native(x_dtensor_native)
            dy_dtensor_native = distribute_tensor(
                dy_global.detach().clone(), manager.device_mesh_subgroups, input_placements
            )
            y_dtensor_native.backward(dy_dtensor_native)

        # Create distributed tensors from global gradients for comparison
        x_grad_dtensor_expected = distribute_tensor(x_global_grad, manager.device_mesh_subgroups, input_placements)
        y_dtensor_expected = distribute_tensor(y_global_expected, manager.device_mesh_subgroups, input_placements)
        weight_grad_dtensor_expected = distribute_tensor(
            weight_global_grad, manager.device_mesh_subgroups, layer.weight.data.placements
        )
        bias_grad_dtensor_expected = distribute_tensor(
            bias_global_grad, manager.device_mesh_subgroups, layer.bias.data.placements
        )

        # Compare results with expected local shards
        torch.testing.assert_close(y_dtensor_expected, y_dtensor_result)
        torch.testing.assert_close(x_grad_dtensor_expected, x_dtensor.grad)
        torch.testing.assert_close(weight_grad_dtensor_expected, layer.weight.grad)
        torch.testing.assert_close(bias_grad_dtensor_expected, layer.bias.grad)

        if compare_to_native_dtensor:
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

            # assert the gradients are identical across all ranks
            assert_all_identical(weight_grad_global_result, manager.group["cp"])
            assert_all_identical(bias_grad_global_result, manager.group["cp"])
    finally:
        DistributedManager.cleanup()
        monkeypatch.undo()


def parallel_assert_dtensor_linear_bf16_gradient_promotion(
    rank: int,
    grid_group_sizes,
    device_type: str,
    backend: str,
    upstream_grad_dtype: torch.dtype = torch.float32,
    avg_over_replicate_param_grad: bool = True,
    env_map=None,
):
    """Regression for bf16-mixed backward in ``LinearParamsReplicated``.

    Reproduces the mixed-dtype path used by bf16-mixed training:
    - bf16 activations/input shards
    - fp32 replicated parameters
    - upstream gradient dtype is controlled by ``upstream_grad_dtype``:
      * fp32: realistic when loss is computed in fp32 — autograd casts the
        gradient to bf16 (matching the forward output) before custom backward
      * bf16: the common case where the upstream layer also runs in bf16

    Then verifies gradients match explicit promote-types reference math.
    """
    monkeypatch = pytest.MonkeyPatch()
    if env_map is not None:
        for var_name, value in env_map.items():
            if value == "<INPUT_RANK>":
                monkeypatch.setenv(var_name, f"{rank}")
                continue
            monkeypatch.setenv(var_name, value)

    try:
        DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
        manager = DistributedManager()

        size_cp = len(manager.group_ranks["cp"])
        size_ring = isqrt(size_cp)
        if size_ring * size_ring != size_cp:
            raise ValueError(f"cp group size {size_cp} is not a square int")

        batch_size = 2
        seq_len_per_rank = 4
        seq_len_global = size_ring * seq_len_per_rank
        in_features = 8
        out_features = 6

        seed_by_rank(0)
        input_placements = [Shard(dim=0), Shard(dim=1), Replicate()]

        # Simulate bf16 activations arriving at this layer.
        x_global_fp32 = torch.randn(batch_size, seq_len_global, in_features, device=manager.device, dtype=torch.float32)
        x_global_bf16 = x_global_fp32.to(torch.bfloat16)

        # Parameters stay fp32 in mixed precision training.
        linear_ref = torch.nn.Linear(in_features, out_features, device=manager.device, dtype=torch.float32)
        state_dict_ref = {k: v.detach().clone() for k, v in linear_ref.state_dict().items()}
        weight_ref = state_dict_ref["weight"]
        bias_ref = state_dict_ref["bias"]

        linear_local = torch.nn.Linear(in_features, out_features, device=manager.device, dtype=torch.float32)
        linear_local.load_state_dict(state_dict_ref)
        layer = LinearParamsReplicated(
            linear_local,
            manager.device_mesh_subgroups,
            avg_over_replicate_param_grad=avg_over_replicate_param_grad,
        )

        x_dtensor = distribute_tensor(x_global_bf16, manager.device_mesh_subgroups, input_placements)
        x_dtensor.requires_grad_(True)

        dy_upstream = torch.randn(
            batch_size, seq_len_global, out_features, device=manager.device, dtype=upstream_grad_dtype
        )

        with torch.autocast(device_type=device_type, dtype=torch.bfloat16, enabled=True):
            y_dtensor = layer(x_dtensor)
            y_expected = torch.nn.functional.linear(x_global_bf16, weight_ref, bias_ref)
        assert y_dtensor.dtype == torch.bfloat16, f"Expected bf16 output under autocast, got {y_dtensor.dtype}"
        torch.testing.assert_close(y_dtensor.full_tensor(), y_expected)

        dy_dtensor = distribute_tensor(dy_upstream, manager.device_mesh_subgroups, input_placements)
        y_dtensor.backward(dy_dtensor)

        assert layer.weight.grad is not None, "Weight gradient should not be None"
        assert layer.bias.grad is not None, "Bias gradient should not be None"
        assert x_dtensor.grad is not None, "Input gradient should not be None"

        weight_grad_result = layer.weight.grad.full_tensor()
        bias_grad_result = layer.bias.grad.full_tensor()
        x_grad_result = x_dtensor.grad.full_tensor()

        assert weight_grad_result.dtype == torch.float32, f"Weight grad should be fp32, got {weight_grad_result.dtype}"
        assert bias_grad_result.dtype == torch.float32, f"Bias grad should be fp32, got {bias_grad_result.dtype}"
        assert x_grad_result.dtype == torch.bfloat16, f"Input grad should be bf16, got {x_grad_result.dtype}"

        dy_effective = dy_upstream.to(dtype=y_dtensor.dtype)

        weight_grad_expected = torch.einsum("...i,...o->io", dy_effective.float(), x_global_bf16.float())
        bias_grad_expected = dy_effective.float().sum(dim=(0, 1))
        bf16_atol, bf16_rtol = torch.testing._comparison.default_tolerances(torch.bfloat16)
        torch.testing.assert_close(weight_grad_result, weight_grad_expected, atol=2 * bf16_atol, rtol=2 * bf16_rtol)
        torch.testing.assert_close(bias_grad_result, bias_grad_expected, atol=2 * bf16_atol, rtol=2 * bf16_rtol)

        x_grad_expected = torch.einsum("...i,io->...o", dy_effective, weight_ref.to(dy_effective.dtype))
        torch.testing.assert_close(x_grad_result, x_grad_expected)

        assert_all_identical(weight_grad_result, manager.group["cp"])
        assert_all_identical(bias_grad_result, manager.group["cp"])
    finally:
        DistributedManager.cleanup()
        monkeypatch.undo()


@pytest.mark.parametrize("avg_over_replicate", [True, False], ids=["avg_reduce", "no_reduce"])
@pytest.mark.parametrize("upstream_grad_dtype", [torch.float32, torch.bfloat16], ids=["dy_fp32", "dy_bf16"])
@pytest.mark.parametrize(
    "setup_env",
    [
        ((1, (2, 2)), True, "cpu", "ENV"),
        ((1, (2, 2)), True, "cuda", "ENV"),
    ],
    indirect=("setup_env",),
    ids=["cpu-dp1-cp2x2", "cuda-dp1-cp2x2"],
)
def test_dtensor_linear_bf16_gradient_promotion(setup_env, upstream_grad_dtype, avg_over_replicate):
    """Goals: bf16-mixed backward produces fp32 param grads from bf16 activations.

    - bf16 activations + fp32 parameters → bf16 output under autocast
    - Weight/bias grads are fp32 (promoted), input grad stays bf16
    - Grads match explicit reference math within bf16 tolerance
    - Cross-rank gradient consistency
    - Tested with both fp32 and bf16 upstream adjoints
    - Tested with and without avg-over-replicate gradient synchronisation
    """
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("CUDA not available")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"Need {world_size} GPUs, have {torch.cuda.device_count()}")

    torch.multiprocessing.set_start_method("spawn", force=True)
    torch.multiprocessing.spawn(
        fn=parallel_assert_dtensor_linear_bf16_gradient_promotion,
        args=(
            grid_group_sizes,
            device_type,
            backend,
            upstream_grad_dtype,
            avg_over_replicate,
            env_per_rank,
        ),
        nprocs=world_size,
        join=True,
    )


@pytest.mark.parametrize("compare_to_native_dtensor", [True, False])
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
def test_dtensor_linear(setup_env, compare_to_native_dtensor):
    """Goals: LinearParamsReplicated forward/backward matches global-tensor reference.

    - Forward output matches serial F.linear on global tensors
    - Backward produces matching input, weight, and bias gradients
    - Cross-rank gradient consistency for replicated parameters
    """
    if compare_to_native_dtensor:
        pytest.skip(
            "Native PyTorch DTensor bugs introduced in the NGC 25.10 container upgrade (from 25.02) "
            "cause this test to fail."
        )

    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    torch.multiprocessing.set_start_method("spawn", force=True)
    torch.multiprocessing.spawn(
        fn=parallel_assert_dtensor_linear,
        args=(
            grid_group_sizes,
            device_type,
            backend,
            compare_to_native_dtensor,
            env_per_rank,
        ),
        nprocs=world_size,
        join=True,
    )
