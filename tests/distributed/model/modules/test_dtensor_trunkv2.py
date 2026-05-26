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

"""Tests for DTensor-based CP trunk modules for Boltz-2.

This module tests the DTensor context-parallel trunk modules against the
serial implementations imported from trunkv2:
- DistogramModule
- BFactorModule
- ContactConditioning

Verification checks (per module):
    V1: single-proc FW input tensor values unchanged by FW and BW
    V2: single-proc BW input tensor values unchanged by BW
    V4a: multi-proc FW input tensor values unchanged by FW
    V4b: multi-proc FW input tensor values unchanged after BW
    V5: multi-proc BW input tensor values unchanged by BW
    V8: multi-proc FW output tensor values close-to single-proc
    V9: multi-proc FW input gradient values close-to single-proc
    V10: multi-proc parameter gradient values close-to single-proc
    V10b: replicated parameter gradients identical across all CP ranks
"""

import math
from collections import OrderedDict

import pytest
import torch
from torch import Tensor
from torch.distributed.tensor import DTensor, Replicate, Shard, distribute_tensor
from torch.testing import assert_close

from boltz.data import const
from boltz.distributed.comm import TransposeComm
from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.modules.trunkv2 import BFactorModule as BFactorModuleDTensor
from boltz.distributed.model.modules.trunkv2 import ContactConditioning as ContactConditioningDTensor
from boltz.distributed.model.modules.trunkv2 import DistogramModule as DistogramModuleDTensor
from boltz.model.modules.trunkv2 import BFactorModule as SerialBFactorModule
from boltz.model.modules.trunkv2 import ContactConditioning as SerialContactConditioning
from boltz.model.modules.trunkv2 import DistogramModule as SerialDistogramModule
from boltz.testing.utils import (
    assert_all_identical,
    assert_tensors_identical,
    init_module_params_uniform,
    init_tensors_uniform,
    make_random_contact_conditioning_features,
    skip_if_cuda_not_avail_or_device_count_less_than_word_size,
    spawn_multiprocessing,
)

SEED = 42


def _assert_unchanged(actual, expected, *, serial=False):
    """Shorthand for assert_tensors_identical with standard immutability kwargs.

    serial=True uses check_storage_offset=True (serial-side V1/V2 checks).
    serial=False uses check_storage_offset=False (worker-side V4/V5 checks on DTensor locals).
    """
    assert_tensors_identical(
        actual,
        expected,
        check_stride=True,
        check_grad=False,
        check_grad_fn=False,
        check_storage_pointer=False,
        check_storage_offset=serial,
    )


def assert_dtensor_distogram(
    rank: int,
    input_example_on_host: Tensor,
    output_ref_on_host: Tensor,
    output_grad_example_on_host: Tensor,
    input_grad_ref_on_host: Tensor,
    parameter_grads_ref_on_host: OrderedDict[str, Tensor | None],
    module_state_dict: dict,
    token_z: int,
    num_bins: int,
    num_distograms: int,
    grid_group_sizes: dict,
    device_type: str,
    backend: str,
    env_map: dict[str, str] | None = None,
):
    """Worker function for distributed DTensor DistogramModule testing.

    Follows the Boltz-1x CP pattern: uses full_tensor() for comparisons,
    assert_tensors_identical for binary identity checks, and grad.full_tensor()
    for backward gradient comparison.
    """
    monkeypatch = pytest.MonkeyPatch()
    if env_map is not None:
        for var_name, value in env_map.items():
            if value == "<INPUT_RANK>":
                monkeypatch.setenv(var_name, f"{rank}")
                continue
            monkeypatch.setenv(var_name, value)

    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    dist_manager = DistributedManager()

    _assert_dtensor_distogram_impl(
        dist_manager,
        input_example_on_host,
        output_ref_on_host,
        output_grad_example_on_host,
        input_grad_ref_on_host,
        parameter_grads_ref_on_host,
        module_state_dict,
        token_z,
        num_bins,
        num_distograms,
        rank,
    )
    DistributedManager.cleanup()
    monkeypatch.undo()


def _assert_dtensor_distogram_impl(
    dist_manager: DistributedManager,
    input_example_on_host: Tensor,
    output_ref_on_host: Tensor,
    output_grad_example_on_host: Tensor,
    input_grad_ref_on_host: Tensor,
    parameter_grads_ref_on_host: OrderedDict[str, Tensor | None],
    module_state_dict: dict,
    token_z: int,
    num_bins: int,
    num_distograms: int,
    rank: int,
) -> None:
    """Inner implementation of distogram DTensor test (extracted for try/finally)."""
    # Move inputs and reference outputs to device
    input_example_on_device = input_example_on_host.detach().to(dist_manager.device).requires_grad_(True)
    output_ref_on_device = output_ref_on_host.detach().to(dist_manager.device)
    output_grad_example_on_device = output_grad_example_on_host.detach().to(dist_manager.device)
    input_grad_ref_on_device = input_grad_ref_on_host.detach().to(dist_manager.device)

    # Create serial module (needed to wrap in DTensor module)
    serial_module = SerialDistogramModule(
        token_z=token_z,
        num_bins=num_bins,
        num_distograms=num_distograms,
    )
    serial_module.load_state_dict(state_dict=module_state_dict)
    serial_module = serial_module.to(dist_manager.device)

    # Create DTensor module
    dtensor_module = DistogramModuleDTensor(
        module=serial_module,
        dist_manager=dist_manager,
        distogram_comm=TransposeComm(
            dist_manager.group["cp"],
            dist_manager.layout_subgroups["cp"],
        ),
    )
    dtensor_module = dtensor_module.to(dist_manager.device)
    dtensor_module = dtensor_module.train()

    # Create input DTensors with pairlike placements
    pairlike_placements = (Shard(0), Shard(1), Shard(2))

    input_as_dtensor = distribute_tensor(
        input_example_on_device,
        device_mesh=dist_manager.device_mesh_subgroups,
        placements=pairlike_placements,
    ).requires_grad_(input_example_on_device.requires_grad)

    output_grad_as_dtensor = distribute_tensor(
        output_grad_example_on_device,
        device_mesh=dist_manager.device_mesh_subgroups,
        placements=pairlike_placements,
    ).requires_grad_(output_grad_example_on_device.requires_grad)

    # V4a setup: clone input for immutability check
    input_example_clone_as_dtensor = input_as_dtensor.detach().clone().requires_grad_(input_as_dtensor.requires_grad)

    # Forward pass
    output_actual_as_dtensor = dtensor_module(input_as_dtensor)

    # V4a: FW input tensor values unchanged by FW (binary identity - atol=0, rtol=0)
    _assert_unchanged(input_as_dtensor.full_tensor(), input_example_clone_as_dtensor.full_tensor(), serial=True)

    # V8: multi-proc FW output tensor values close-to single-proc
    # Tolerance: fp32 default (atol=1e-5, rtol=1.3e-6). The distogram is
    # always tested in fp32. The only numerical difference is summation order
    # in symmetrize (z + z.T) which introduces ~N*eps ≈ 16*1.2e-7 ≈ 2e-6 error.
    assert_close(output_actual_as_dtensor.full_tensor(), output_ref_on_device, atol=1e-5, rtol=1.3e-6)

    # Verify 5D output shape
    assert (
        output_actual_as_dtensor.shape[3] == num_distograms
    ), f"Expected dim 3 = {num_distograms}, got {output_actual_as_dtensor.shape[3]}"
    assert (
        output_actual_as_dtensor.shape[4] == num_bins
    ), f"Expected dim 4 = {num_bins}, got {output_actual_as_dtensor.shape[4]}"

    # V5 setup: save output for comparison
    output_grad_example_clone_as_dtensor = (
        output_grad_as_dtensor.detach().clone().requires_grad_(output_grad_as_dtensor.requires_grad)
    )
    output_actual_clone_as_dtensor = (
        output_actual_as_dtensor.detach().clone().requires_grad_(output_actual_as_dtensor.requires_grad)
    )

    # Backward pass
    output_actual_as_dtensor.backward(output_grad_as_dtensor)

    # V4b: FW input tensor values unchanged after backward (binary identity - atol=0, rtol=0)
    _assert_unchanged(input_as_dtensor.full_tensor(), input_example_clone_as_dtensor.full_tensor(), serial=True)

    # V5: BW input tensor values unchanged by BW (binary identity - atol=0, rtol=0)
    _assert_unchanged(output_actual_as_dtensor.full_tensor(), output_actual_clone_as_dtensor.full_tensor(), serial=True)
    _assert_unchanged(
        output_grad_as_dtensor.full_tensor(), output_grad_example_clone_as_dtensor.full_tensor(), serial=True
    )

    # V9: multi-proc FW input gradient values close-to single-proc
    # Use grad.full_tensor() because DTensor grad can be in Partial(Sum) placement
    assert input_as_dtensor.grad is not None, "Input DTensor gradient is None - trivial equality guard failed"
    assert input_grad_ref_on_device is not None, "Reference input gradient is None - test setup error"
    input_grad_actual_full: Tensor = input_as_dtensor.grad.full_tensor()
    assert_close(input_grad_actual_full, input_grad_ref_on_device, atol=1e-5, rtol=1.3e-6)

    # Non-vacuous: input gradient must be non-zero
    assert input_grad_actual_full.abs().sum() > 0, "Input gradient is all-zero — backward did not propagate"

    # V10: multi-proc parameter gradient values close-to single-proc
    param_names_checked = []
    for name, param in dtensor_module.named_parameters():
        if name not in parameter_grads_ref_on_host:
            msg = f"Module parameter {name} not in parameter_grads_ref_on_host"
            raise ValueError(msg)

        grad_ref = parameter_grads_ref_on_host[name]

        if param.grad is None and grad_ref is None:
            msg = (
                f"Both actual and reference gradients are None for {name} - "
                "trivial equality, test cannot verify correctness"
            )
            raise ValueError(msg)

        if (param.grad is None) != (grad_ref is None):
            msg = f"Inconsistent grad state for {name} on rank {rank}: result={param.grad}, ref={grad_ref}"
            raise ValueError(msg)

        param_names_checked.append(name)
        if grad_ref is not None:
            # Use full_tensor() to reduce Partial(Sum) gradients
            grad_actual = param.grad.full_tensor() if isinstance(param.grad, DTensor) else param.grad
            param_name = name
            assert_close(
                grad_actual,
                grad_ref.to(dist_manager.device),
                atol=1e-5,
                rtol=1.3e-6,
                msg=lambda m, n=param_name: f"Rank {rank} {n} grad mismatch\n{m}",
            )

            # V10b: replicated parameter gradients identical across all CP ranks
            assert_all_identical(grad_actual.detach(), dist_manager.group["cp"])

            # Non-vacuous: parameter gradient must be non-zero
            assert grad_actual.abs().sum() > 0, f"Rank {rank} {param_name} gradient is all-zero"

    assert (
        len(param_names_checked) >= 2
    ), f"Expected at least 2 parameters (weight, bias), but only checked: {param_names_checked}"


def get_example_input_and_reference_output(
    grid_group_sizes: dict,
    B: int,
    num_tokens_per_device_grid_unit: int,
    token_z: int,
    num_bins: int,
    num_distograms: int,
    dtype_for_test: torch.dtype = torch.float32,
    device_for_test: str = "cpu",
    seed: int = SEED,
):
    """Generate example input and reference output for testing."""
    with torch.random.fork_rng(devices=[], enabled=True):
        torch.manual_seed(seed)

        num_tokens = num_tokens_per_device_grid_unit * grid_group_sizes["cp"][0]

        min_init_val = -0.5
        max_init_val = 0.5

        input_example = torch.empty(
            (B, num_tokens, num_tokens, token_z),
            device=device_for_test,
            dtype=dtype_for_test,
            requires_grad=True,
        )
        init_tensors_uniform([input_example], low=min_init_val, high=max_init_val)
        input_example_copy = input_example.detach().clone().requires_grad_(input_example.requires_grad)

        module = SerialDistogramModule(token_z, num_bins, num_distograms=num_distograms)
        init_module_params_uniform(module, low=min_init_val, high=max_init_val)
        module = module.to(device_for_test)
        module_state_dict = module.state_dict()

        # Run serial forward
        output_ref = module(input_example)

        # V1a: single-proc FW input tensor values unchanged
        _assert_unchanged(input_example, input_example_copy, serial=True)

        output_grad_example = torch.empty_like(output_ref)
        init_tensors_uniform([output_grad_example], low=min_init_val, high=max_init_val)
        output_grad_example_copy = output_grad_example.detach().clone()

        # Serial backward
        torch.autograd.backward([output_ref], [output_grad_example])

        # V1b: single-proc FW input tensor values unchanged after backward
        _assert_unchanged(input_example, input_example_copy, serial=True)

        # V2: single-proc BW input tensor values unchanged
        _assert_unchanged(output_grad_example, output_grad_example_copy, serial=True)

        # Get parameter gradients
        parameter_grads_ref = OrderedDict()
        for name, param in module.named_parameters():
            if param.grad is not None:
                parameter_grads_ref[name] = param.grad.detach().cpu().clone()
            else:
                parameter_grads_ref[name] = None

        # To host for output
        input_example_on_host = input_example.detach().cpu().clone()
        output_ref_on_host = output_ref.detach().cpu().clone()
        output_grad_example_on_host = output_grad_example.detach().cpu().clone()
        input_example_grad_on_host = input_example.grad.detach().cpu().clone()

    return (
        input_example_on_host,
        output_ref_on_host,
        output_grad_example_on_host,
        input_example_grad_on_host,
        parameter_grads_ref,
        module_state_dict,
    )


@pytest.mark.parametrize(
    "setup_env, num_distograms",
    [
        # CUDA dp=1 cp=(1,1): serial-equivalent sanity check with D=1 (v1 fallback)
        (((1, (1, 1)), True, "cuda", "ENV"), 1),
        # CUDA dp=2 cp=(1,1): DP-only path (2 GPUs)
        (((2, (1, 1)), True, "cuda", "ENV"), 1),
        # CUDA dp=2 cp=(2,2): full DP+CP with multi-distogram (D=3)
        (((2, (2, 2)), True, "cuda", "ENV"), 3),
        # CPU dp=2 cp=(3,3): non-power-of-two CP baseline without GPUs
        (((2, (3, 3)), True, "cpu", "ENV"), 1),
    ],
    indirect=("setup_env",),
    ids=[
        "cuda-dp1-cp1x1-D1",
        "cuda-dp2-cp1x1-D1",
        "cuda-dp2-cp2x2-D3",
        "cpu-dp2-cp3x3-D1",
    ],
)
def test_dtensor_distogram_forward_backward(
    num_distograms: int,
    setup_env: dict,
    seed: int = SEED,
):
    """Test that DTensor DistogramModule matches serial Boltz-2 implementation."""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    skip_if_cuda_not_avail_or_device_count_less_than_word_size(
        device_type=device_type,
        world_size=world_size,
    )

    B = 2
    num_tokens_per_device_grid_unit = 8
    token_z = 32
    num_bins = 16

    (
        input_example_on_host,
        output_ref_on_host,
        output_grad_example_on_host,
        input_grad_ref_on_host,
        parameter_grads_ref,
        module_state_dict,
    ) = get_example_input_and_reference_output(
        grid_group_sizes,
        B=B,
        num_tokens_per_device_grid_unit=num_tokens_per_device_grid_unit,
        token_z=token_z,
        num_bins=num_bins,
        num_distograms=num_distograms,
        dtype_for_test=torch.float32,
        device_for_test=device_type,
        seed=seed,
    )

    spawn_multiprocessing(
        assert_dtensor_distogram,
        world_size,
        input_example_on_host,
        output_ref_on_host,
        output_grad_example_on_host,
        input_grad_ref_on_host,
        parameter_grads_ref,
        module_state_dict,
        token_z,
        num_bins,
        num_distograms,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
    )


# ====================================================================== #
#  BFactorModule parity tests                                            #
# ====================================================================== #


def _worker_bfactor_parity(
    rank: int,
    input_on_host: Tensor,
    output_ref_on_host: Tensor,
    grad_output_on_host: Tensor,
    input_grad_ref_on_host: Tensor,
    param_grads_ref: OrderedDict[str, Tensor | None],
    state_dict: dict,
    token_s: int,
    num_bins: int,
    dtype: torch.dtype,
    grid_group_sizes: dict,
    device_type: str,
    backend: str,
    env_map: dict[str, str] | None = None,
):
    """Worker: compare distributed BFactorModule against serial reference.

    Performs V4a, V4b, V5, V8, V9, V10, V10b checks.
    """
    monkeypatch = pytest.MonkeyPatch()
    if env_map is not None:
        for var_name, value in env_map.items():
            if value == "<INPUT_RANK>":
                monkeypatch.setenv(var_name, f"{rank}")
                continue
            monkeypatch.setenv(var_name, value)

    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    dm = DistributedManager()

    serial = SerialBFactorModule(token_s=token_s, num_bins=num_bins)
    serial = serial.to(device=dm.device, dtype=dtype)
    serial.load_state_dict(state_dict)

    dist_mod = BFactorModuleDTensor(serial, device_mesh=dm.device_mesh_subgroups)
    dist_mod.train()

    # Single representation placements: shard B on dp, N on cp0, replicate on cp1
    single_placements = (Shard(0), Shard(1)) + (Replicate(),) * (dm.device_mesh_subgroups.ndim - 2)
    x = distribute_tensor(
        input_on_host.to(device=dm.device, dtype=dtype), dm.device_mesh_subgroups, single_placements
    ).requires_grad_(True)

    # V4a setup: clone input
    x_clone = x.detach().clone().requires_grad_(x.requires_grad)

    out = dist_mod(x)

    # V4a: FW input unchanged
    _assert_unchanged(x.to_local(), x_clone.to_local())

    # V8: forward parity
    fw_atol = 1e-2 if dtype == torch.bfloat16 else 1e-5
    fw_rtol = 1e-2 if dtype == torch.bfloat16 else 1.3e-6
    assert_close(
        out.full_tensor(),
        output_ref_on_host.to(device=dm.device, dtype=dtype),
        atol=fw_atol,
        rtol=fw_rtol,
    )

    # V5 setup: clone output and grad
    grad_out = distribute_tensor(
        grad_output_on_host.to(device=dm.device, dtype=dtype), dm.device_mesh_subgroups, single_placements
    )
    out_clone = out.detach().clone().requires_grad_(out.requires_grad)
    grad_out_clone = grad_out.detach().clone().requires_grad_(grad_out.requires_grad)

    out.backward(grad_out)

    # V4b: FW input unchanged after backward
    _assert_unchanged(x.to_local(), x_clone.to_local())

    # V5: BW inputs unchanged
    _assert_unchanged(out.to_local(), out_clone.to_local())
    _assert_unchanged(grad_out.to_local(), grad_out_clone.to_local())

    # V9: input gradient parity
    grad_atol = 5e-2 if dtype == torch.bfloat16 else 1e-5
    grad_rtol = 5e-2 if dtype == torch.bfloat16 else 1.3e-6
    assert x.grad is not None, "Input gradient is None"
    assert_close(
        x.grad.full_tensor(),
        input_grad_ref_on_host.to(device=dm.device, dtype=dtype),
        atol=grad_atol,
        rtol=grad_rtol,
    )

    # V10: parameter gradient parity
    for name, param in dist_mod.named_parameters():
        ref = param_grads_ref.get(name)
        assert param.grad is not None, f"Param {name} grad is None"
        assert ref is not None, f"Reference grad for {name} is None"
        actual = param.grad.full_tensor() if isinstance(param.grad, DTensor) else param.grad
        assert_close(
            actual,
            ref.to(device=dm.device, dtype=dtype),
            atol=grad_atol,
            rtol=grad_rtol,
            msg=lambda m, n=name: f"{n} grad mismatch\n{m}",
        )

        # V10b: replicated parameter gradients identical across all CP ranks
        grad_full = actual.detach()
        assert_all_identical(grad_full, dm.group["cp"])

        # Non-vacuous: parameter gradients must be non-zero
        assert actual.abs().sum() > 0, f"Param {name} gradient is all-zero"

    # Non-vacuous: input gradient must be non-zero
    assert x.grad.full_tensor().abs().sum() > 0, "Input gradient is all-zero"

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env, dtype",
    [
        # CUDA dp=1 cp=(1,1): serial-equivalent sanity check (1 GPU)
        (((1, (1, 1)), True, "cuda", "ENV"), torch.float32),
        # CUDA dp=1 cp=(2,2): fp32 baseline
        (((1, (2, 2)), True, "cuda", "ENV"), torch.float32),
        # CUDA dp=2 cp=(2,2): DP + CP
        (((2, (2, 2)), True, "cuda", "ENV"), torch.float32),
        # CUDA dp=2 cp=(1,1): bf16 mixed precision, 2-GPU
        (((2, (1, 1)), True, "cuda", "ENV"), torch.bfloat16),
        # CUDA dp=1 cp=(2,2): bf16 with actual CP, 4-GPU
        (((1, (2, 2)), True, "cuda", "ENV"), torch.bfloat16),
        # CPU dp=2 cp=(3,3): DP + non-power-of-two CP for CPU-only CI
        (((2, (3, 3)), True, "cpu", "ENV"), torch.float32),
    ],
    indirect=("setup_env",),
    ids=[
        "cuda-dp1-cp1x1-fp32",
        "cuda-dp1-cp2x2-fp32",
        "cuda-dp2-cp2x2-fp32",
        "cuda-dp2-cp1x1-bf16",
        "cuda-dp1-cp2x2-bf16",
        "cpu-dp2-cp3x3-fp32",
    ],
)
def test_dtensor_bfactor_forward_backward(setup_env, dtype: torch.dtype):
    """BFactorModule: distributed output and gradients match serial reference."""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    skip_if_cuda_not_avail_or_device_count_less_than_word_size(device_type=device_type, world_size=world_size)

    B, token_s, num_bins = 2, 32, 8
    N = 8 * grid_group_sizes["cp"][0]

    with torch.random.fork_rng(devices=[], enabled=True):
        torch.manual_seed(SEED)
        x = torch.randn(B, N, token_s, requires_grad=True)
        serial = SerialBFactorModule(token_s=token_s, num_bins=num_bins)
        init_module_params_uniform(serial, low=-0.5, high=0.5)
        serial = serial.to(dtype=dtype)
        state_dict = {k: v.cpu().clone() for k, v in serial.state_dict().items()}

        x_typed = x.to(dtype=dtype).detach().requires_grad_(True)
        x_copy = x_typed.detach().clone().requires_grad_(True)

        out_ref = serial(x_typed)

        # V1a: serial FW input unchanged
        _assert_unchanged(x_typed, x_copy, serial=True)

        grad_out = torch.randn_like(out_ref)
        grad_out_copy = grad_out.detach().clone()
        out_ref.backward(grad_out)

        # V1b: serial FW input unchanged after backward
        _assert_unchanged(x_typed, x_copy, serial=True)
        # V2: serial BW input (grad_out) unchanged
        _assert_unchanged(grad_out, grad_out_copy, serial=True)

        param_grads = OrderedDict(
            (n, p.grad.detach().cpu().clone()) for n, p in serial.named_parameters() if p.grad is not None
        )

    spawn_multiprocessing(
        _worker_bfactor_parity,
        world_size,
        x_typed.detach().cpu(),
        out_ref.detach().cpu(),
        grad_out.detach().cpu(),
        x_typed.grad.detach().cpu(),
        param_grads,
        state_dict,
        token_s,
        num_bins,
        dtype,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
    )


# ====================================================================== #
#  ContactConditioning parity tests                                      #
# ====================================================================== #


def _worker_contact_conditioning_parity(
    rank: int,
    cc_input_on_host: Tensor,
    ct_input_on_host: Tensor,
    output_ref_on_host: Tensor,
    grad_output_on_host: Tensor,
    cc_grad_ref_on_host: Tensor,
    ct_grad_ref_on_host: Tensor,
    param_grads_ref: OrderedDict[str, Tensor | None],
    state_dict: dict,
    token_z: int,
    dtype: torch.dtype,
    grid_group_sizes: dict,
    device_type: str,
    backend: str,
    env_map: dict[str, str] | None = None,
):
    """Worker: compare distributed ContactConditioning against serial reference.

    Performs V4a, V4b, V5, V8, V9, V10, V10b checks.
    """
    monkeypatch = pytest.MonkeyPatch()
    if env_map is not None:
        for var_name, value in env_map.items():
            if value == "<INPUT_RANK>":
                monkeypatch.setenv(var_name, f"{rank}")
                continue
            monkeypatch.setenv(var_name, value)

    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    dm = DistributedManager()

    serial = SerialContactConditioning(token_z=token_z, cutoff_min=0.0, cutoff_max=22.0)
    serial = serial.to(device=dm.device, dtype=dtype)
    serial.load_state_dict(state_dict)

    dist_mod = ContactConditioningDTensor(serial, device_mesh=dm.device_mesh_subgroups)
    dist_mod.train()

    pair_placements = (Shard(0), Shard(1), Shard(2))

    cc_dt = distribute_tensor(
        cc_input_on_host.to(device=dm.device, dtype=dtype), dm.device_mesh_subgroups, pair_placements
    ).requires_grad_(True)
    ct_dt = distribute_tensor(
        ct_input_on_host.to(device=dm.device, dtype=dtype), dm.device_mesh_subgroups, pair_placements
    ).requires_grad_(True)

    # V4a setup: clone inputs
    cc_dt_clone = cc_dt.detach().clone().requires_grad_(cc_dt.requires_grad)
    ct_dt_clone = ct_dt.detach().clone().requires_grad_(ct_dt.requires_grad)

    out = dist_mod({"contact_conditioning": cc_dt, "contact_threshold": ct_dt})

    # V4a: FW inputs unchanged
    _assert_unchanged(cc_dt.to_local(), cc_dt_clone.to_local())
    _assert_unchanged(ct_dt.to_local(), ct_dt_clone.to_local())

    # V8: forward parity
    fw_atol = 1e-2 if dtype == torch.bfloat16 else 1e-5
    fw_rtol = 1e-2 if dtype == torch.bfloat16 else 1.3e-6
    assert_close(
        out.full_tensor(),
        output_ref_on_host.to(device=dm.device, dtype=dtype),
        atol=fw_atol,
        rtol=fw_rtol,
    )

    # V5 setup: clone output and grad
    grad_out = distribute_tensor(
        grad_output_on_host.to(device=dm.device, dtype=dtype), dm.device_mesh_subgroups, pair_placements
    )
    out_clone = out.detach().clone().requires_grad_(out.requires_grad)
    grad_out_clone = grad_out.detach().clone().requires_grad_(grad_out.requires_grad)

    out.backward(grad_out)

    # V4b: FW inputs unchanged after backward
    _assert_unchanged(cc_dt.to_local(), cc_dt_clone.to_local())

    # V5: BW inputs unchanged
    _assert_unchanged(out.to_local(), out_clone.to_local())
    _assert_unchanged(grad_out.to_local(), grad_out_clone.to_local())

    # V9: input gradient parity
    #
    # Tolerance rationale (fp32): the encoder.weight gradient is computed as
    # dW = x.T @ grad_out, summing over K = B*N*N elements.  Serial sums all K
    # terms in one matmul; distributed partitions into K/W terms per rank then
    # all-reduce-sums across W ranks.  The different accumulation orders produce
    # a relative error of O(sqrt(K) * eps).  For cpu-dp2-cp3x3, K = 2*24*24 =
    # 1152, giving sqrt(1152) * eps ≈ 4.0e-6.  A 2x safety factor covers
    # intermediate ops (Fourier embedding, masking arithmetic).
    if dtype == torch.bfloat16:
        grad_atol, grad_rtol = 5e-2, 5e-2
    else:
        K = cc_input_on_host.shape[0] * cc_input_on_host.shape[1] * cc_input_on_host.shape[2]
        fp32_eps = torch.finfo(torch.float32).eps
        accum_rtol = 2.0 * math.sqrt(K) * fp32_eps
        grad_rtol = max(1.3e-6, accum_rtol)
        grad_atol = max(1e-5, 10.0 * accum_rtol)
    assert cc_dt.grad is not None, "contact_conditioning input gradient is None"
    assert_close(
        cc_dt.grad.full_tensor(),
        cc_grad_ref_on_host.to(device=dm.device, dtype=dtype),
        atol=grad_atol,
        rtol=grad_rtol,
    )
    if ct_grad_ref_on_host is not None:
        assert ct_dt.grad is not None, "contact_threshold input gradient is None"
        assert_close(
            ct_dt.grad.full_tensor(),
            ct_grad_ref_on_host.to(device=dm.device, dtype=dtype),
            atol=grad_atol,
            rtol=grad_rtol,
        )

    # V10: parameter gradient parity (trainable params only)
    for name, param in dist_mod.named_parameters():
        if not param.requires_grad:
            continue
        ref = param_grads_ref.get(name)
        if ref is None:
            continue
        assert param.grad is not None, f"Param {name} grad is None"
        actual = param.grad.full_tensor() if isinstance(param.grad, DTensor) else param.grad
        assert_close(
            actual,
            ref.to(device=dm.device, dtype=dtype),
            atol=grad_atol,
            rtol=grad_rtol,
            msg=lambda m, n=name: f"{n} grad mismatch\n{m}",
        )

        # V10b: replicated parameter gradients identical across all CP ranks
        grad_full = actual.detach()
        assert_all_identical(grad_full, dm.group["cp"])

    # Non-vacuous: encoding_unspecified and encoding_unselected must have
    # non-zero gradients to prove the UNSPECIFIED and UNSELECTED masking
    # branches were exercised.  A zero grad would mean the test data
    # didn't include that contact type, making V10 trivially pass.
    for enc_name in ("encoding_unspecified", "encoding_unselected"):
        enc_param = dict(dist_mod.named_parameters())[enc_name]
        assert enc_param.grad is not None, f"{enc_name} grad is None"
        enc_grad = enc_param.grad.full_tensor() if isinstance(enc_param.grad, DTensor) else enc_param.grad
        assert enc_grad.abs().sum() > 0, (
            f"{enc_name} gradient is all-zero — UNSPECIFIED/UNSELECTED masking "
            f"branch was not exercised. Test data must include this contact type."
        )

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env, dtype",
    [
        # CUDA dp=1 cp=(1,1): serial-equivalent sanity check (1 GPU)
        (((1, (1, 1)), True, "cuda", "ENV"), torch.float32),
        # CUDA dp=1 cp=(2,2): fp32 baseline
        (((1, (2, 2)), True, "cuda", "ENV"), torch.float32),
        # CUDA dp=2 cp=(2,2): DP + CP
        (((2, (2, 2)), True, "cuda", "ENV"), torch.float32),
        # CUDA dp=2 cp=(1,1): bf16 mixed precision, 2-GPU
        (((2, (1, 1)), True, "cuda", "ENV"), torch.bfloat16),
        # CUDA dp=1 cp=(2,2): bf16 with actual CP, 4-GPU
        (((1, (2, 2)), True, "cuda", "ENV"), torch.bfloat16),
        # CPU dp=2 cp=(3,3): DP + non-power-of-two CP for CPU-only CI
        (((2, (3, 3)), True, "cpu", "ENV"), torch.float32),
    ],
    indirect=("setup_env",),
    ids=[
        "cuda-dp1-cp1x1-fp32",
        "cuda-dp1-cp2x2-fp32",
        "cuda-dp2-cp2x2-fp32",
        "cuda-dp2-cp1x1-bf16",
        "cuda-dp1-cp2x2-bf16",
        "cpu-dp2-cp3x3-fp32",
    ],
)
def test_dtensor_contact_conditioning_forward_backward(setup_env, dtype: torch.dtype):
    """ContactConditioning: distributed output and gradients match serial reference."""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    skip_if_cuda_not_avail_or_device_count_less_than_word_size(device_type=device_type, world_size=world_size)

    B, token_z = 2, 16
    N = 8 * grid_group_sizes["cp"][0]
    num_cc_types = len(const.contact_conditioning_info)

    with torch.random.fork_rng(devices=[], enabled=True):
        torch.manual_seed(SEED)

        serial = SerialContactConditioning(token_z=token_z, cutoff_min=0.0, cutoff_max=22.0)
        init_module_params_uniform(serial, low=-0.5, high=0.5)
        serial = serial.to(dtype=dtype)
        state_dict = {k: v.cpu().clone() for k, v in serial.state_dict().items()}

        cc_data, ct_data = make_random_contact_conditioning_features(
            B,
            N,
            num_cc_types,
            dtype=dtype,
            seed=SEED,
        )
        cc_input = cc_data.requires_grad_(True)
        ct_input = ct_data.requires_grad_(True)

        cc_copy = cc_input.detach().clone().requires_grad_(True)
        ct_copy = ct_input.detach().clone().requires_grad_(True)

        feats = {"contact_conditioning": cc_input, "contact_threshold": ct_input}
        out_ref = serial(feats)

        # V1a: serial FW inputs unchanged
        _assert_unchanged(cc_input, cc_copy, serial=True)
        _assert_unchanged(ct_input, ct_copy, serial=True)

        grad_out = torch.randn_like(out_ref)
        grad_out_copy = grad_out.detach().clone()
        out_ref.backward(grad_out)

        # V1b: serial FW inputs unchanged after backward
        _assert_unchanged(cc_input, cc_copy, serial=True)
        _assert_unchanged(ct_input, ct_copy, serial=True)
        # V2: serial BW input (grad_out) unchanged
        _assert_unchanged(grad_out, grad_out_copy, serial=True)

        param_grads = OrderedDict()
        for n, p in serial.named_parameters():
            if p.requires_grad and p.grad is not None:
                param_grads[n] = p.grad.detach().cpu().clone()

        ct_grad = ct_input.grad.detach().cpu().clone() if ct_input.grad is not None else None

    spawn_multiprocessing(
        _worker_contact_conditioning_parity,
        world_size,
        cc_input.detach().cpu(),
        ct_input.detach().cpu(),
        out_ref.detach().cpu(),
        grad_out.detach().cpu(),
        cc_input.grad.detach().cpu(),
        ct_grad,
        param_grads,
        state_dict,
        token_z,
        dtype,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
    )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
