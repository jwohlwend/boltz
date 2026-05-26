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

"""Tests for DTensor ConditionedTransitionBlock module.

Tests both Boltz-1x and Boltz-2 serial ConditionedTransitionBlock modules against
the unified DTensor implementation, verifying forward and backward equivalence.

"""

import pytest
import torch
from torch.distributed.tensor import DTensor, Replicate, Shard, distribute_tensor

from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.modules.transformers import ConditionedTransitionBlock as DTensorCTB
from boltz.model.modules.transformers import ConditionedTransitionBlock as CTBSerialBoltz1
from boltz.model.modules.transformersv2 import ConditionedTransitionBlock as CTBSerialBoltz2
from boltz.testing.utils import (
    assert_tensors_identical,
    seed_by_rank,
    spawn_multiprocessing,
)


def parallel_assert_conditioned_transition_block(
    rank: int,
    grid_group_sizes,
    device_type: str,
    backend: str,
    env_per_rank,
    serial_module_version: str,
    dim_single: int,
    dim_single_cond: int,
    B: int,
    N: int,
    layer_state_dict,
    a_global_host: torch.Tensor,
    s_global_host: torch.Tensor,
    d_out_global_host: torch.Tensor,
    out_expected_global_host: torch.Tensor,
    d_a_expected_global_host: torch.Tensor,
    d_s_expected_global_host: torch.Tensor,
    expected_param_grads_global_host_dict: dict[str, torch.Tensor],
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

    # Create serial module from state dict
    CTBSerial = CTBSerialBoltz1 if serial_module_version == "boltz1" else CTBSerialBoltz2
    module_serial = CTBSerial(dim_single=dim_single, dim_single_cond=dim_single_cond)
    module_serial.load_state_dict(layer_state_dict)
    module_serial = module_serial.to(device=manager.device).train()

    # Create DTensor module from serial
    module_dt = DTensorCTB(
        conditioned_trans_block=module_serial,
        device_mesh=manager.device_mesh_subgroups,
    ).train()

    # Placements
    placements = (Shard(0), Shard(1), Replicate()) if manager.device_mesh_subgroups.ndim == 3 else (Shard(0), Shard(1))

    a_dt = distribute_tensor(
        a_global_host.to(device=manager.device), manager.device_mesh_subgroups, placements
    ).requires_grad_(True)
    s_dt = distribute_tensor(
        s_global_host.to(device=manager.device), manager.device_mesh_subgroups, placements
    ).requires_grad_(True)

    # Copies to verify inputs aren't modified
    a_dt_copy = a_dt.detach().clone().requires_grad_(True)
    s_dt_copy = s_dt.detach().clone().requires_grad_(True)

    # Forward pass
    out_dt: DTensor = module_dt(a_dt, s_dt)

    # Ensure no input mutation
    assert_tensors_identical(a_dt_copy.to_local(), a_dt.to_local(), check_grad=False, check_grad_fn=False)
    assert_tensors_identical(s_dt_copy.to_local(), s_dt.to_local(), check_grad=False, check_grad_fn=False)

    # Forward compare (full gather)
    torch.testing.assert_close(out_dt.full_tensor().cpu(), out_expected_global_host)

    # Backward pass
    d_out_dt = distribute_tensor(d_out_global_host.to(device=manager.device), manager.device_mesh_subgroups, placements)
    out_dt.backward(d_out_dt)

    # Compare input gradients
    torch.testing.assert_close(a_dt.grad.full_tensor().cpu(), d_a_expected_global_host)
    torch.testing.assert_close(s_dt.grad.full_tensor().cpu(), d_s_expected_global_host)

    # Compare parameter gradients
    for name, param in module_dt.named_parameters():
        assert param.grad is not None, f"Parameter {name} has no gradient"
        expected_grad = expected_param_grads_global_host_dict[name]
        torch.testing.assert_close(
            param.grad.full_tensor().cpu(),
            expected_grad,
            msg=lambda m: f"Parameter gradient mismatch for {name}: {m}",
        )

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env",
    [
        ((1, (2, 2)), True, "cuda", "ENV"),
        ((2, (2, 2)), True, "cuda", "ENV"),
    ],
    indirect=("setup_env",),
    ids=lambda x: f"dp:{x[0][0]}, cp:{x[0][1]}, specify_method:{x[1]}, device_type:{x[2]}, method_init:{x[3]}",
)
@pytest.mark.parametrize("serial_module_version", ["boltz1", "boltz2"])
def test_conditioned_transition_block(setup_env, serial_module_version: str):
    """Test ConditionedTransitionBlock DTensor vs serial equivalence."""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() < {world_size}")

    # Module dimensions
    dim_single = 64
    dim_single_cond = 32

    # Data dimensions
    size_ring = grid_group_sizes["cp"][0]
    B = 2 * grid_group_sizes["dp"]
    N = size_ring * 8

    seed_by_rank(0, seed=42)

    # Create serial module
    CTBSerial = CTBSerialBoltz1 if serial_module_version == "boltz1" else CTBSerialBoltz2
    module_serial = CTBSerial(dim_single=dim_single, dim_single_cond=dim_single_cond)
    module_serial = module_serial.train()
    layer_state_dict = module_serial.state_dict()

    # Create input tensors
    a_global = torch.randn(B, N, dim_single, requires_grad=True)
    s_global = torch.randn(B, N, dim_single_cond, requires_grad=True)

    # Serial forward pass
    out_serial = module_serial(a_global, s_global)

    # Create upstream gradient
    d_out = torch.randn_like(out_serial)

    # Serial backward pass
    out_serial.backward(d_out)

    # Collect expected results
    out_expected = out_serial.detach().clone().cpu()
    d_a_expected = a_global.grad.detach().clone().cpu()
    d_s_expected = s_global.grad.detach().clone().cpu()

    expected_param_grads = {}
    for name, param in module_serial.named_parameters():
        assert param.grad is not None, f"Serial parameter {name} has no gradient"
        expected_param_grads[name] = param.grad.detach().clone().cpu()

    # Launch parallel test
    spawn_multiprocessing(
        parallel_assert_conditioned_transition_block,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        serial_module_version,
        dim_single,
        dim_single_cond,
        B,
        N,
        layer_state_dict,
        a_global.detach().clone().cpu(),
        s_global.detach().clone().cpu(),
        d_out.detach().clone().cpu(),
        out_expected,
        d_a_expected,
        d_s_expected,
        expected_param_grads,
    )
