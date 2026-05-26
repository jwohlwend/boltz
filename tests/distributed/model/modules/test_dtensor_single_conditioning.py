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

"""Tests for DTensor SingleConditioning module.

Tests both Boltz-1x and Boltz-2 serial SingleConditioning modules against the unified
DTensor SingleConditioning implementation, verifying forward and backward equivalence.

The functional differences between V1 and V2 SingleConditioning are:
  - ``disable_times``: V2-only flag. When True, fourier time embedding is absent.
  - ``v1_input_layout``: V1 uses a wider input_dim for norm_single/single_embed
    (``2*token_s + 2*num_tokens + 1 + len(pocket_contact_info)``), while V2 uses
    ``2*token_s``.  When True, the V1 serial class is used.

The serial module version is inferred from these flags:
  - ``v1_input_layout=True``  → ``SingleConditioningBoltz1`` (V1 does not support disable_times)
  - ``v1_input_layout=False`` → ``SingleConditioningBoltz2``

Uses float64 with default tolerance for exact comparison.
"""

import pytest
import torch
from torch.distributed.tensor import Replicate, Shard, distribute_tensor

from boltz.data import const
from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.modules.encoders import SingleConditioning as DTensorSingleConditioning
from boltz.model.modules.encoders import SingleConditioning as SingleConditioningBoltz1
from boltz.model.modules.encodersv2 import SingleConditioning as SingleConditioningBoltz2
from boltz.testing.utils import (
    assert_tensors_identical,
    init_module_params_uniform,
    init_tensors_uniform,
    seed_by_rank,
    spawn_multiprocessing,
)


def parallel_assert_single_conditioning(
    rank: int,
    grid_group_sizes,
    device_type: str,
    backend: str,
    env_per_rank,
    serial_class_tag: str,  # "boltz1" or "boltz2" — inferred from flags
    layer_state_dict,
    dtype: torch.dtype,
    # Input tensors (global, on host)
    times_global_host: torch.Tensor,
    s_trunk_global_host: torch.Tensor,
    s_inputs_global_host: torch.Tensor,
    # Expected outputs
    s_expected_global_host: torch.Tensor,
    normed_fourier_expected_global_host: torch.Tensor | None,
    # Upstream gradients
    d_s_global_host: torch.Tensor,
    d_normed_fourier_global_host: torch.Tensor | None,
    # Expected input grads
    d_s_trunk_expected_global_host: torch.Tensor,
    d_s_inputs_expected_global_host: torch.Tensor,
    # Expected parameter grads
    expected_param_grads_global_host_dict: dict[str, torch.Tensor],
    # Module constructor kwargs
    module_kwargs: dict,
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
    if serial_class_tag == "boltz1":
        module_serial = SingleConditioningBoltz1(**module_kwargs)
    else:
        module_serial = SingleConditioningBoltz2(**module_kwargs)
    module_serial = module_serial.to(device=manager.device, dtype=dtype)
    module_serial.load_state_dict(layer_state_dict)
    module_serial = module_serial.train()

    # Create DTensor module from serial
    module_dt = DTensorSingleConditioning(
        layer=module_serial,
        device_mesh=manager.device_mesh_subgroups,
    ).train()

    # Placements
    placements_times = (Shard(0), Replicate(), Replicate())
    placements_s = (Shard(0), Shard(1), Replicate())

    # Distribute inputs
    times_dt = distribute_tensor(
        times_global_host.to(device=manager.device, dtype=dtype),
        manager.device_mesh_subgroups,
        placements_times,
    ).requires_grad_(False)  # times has no gradient in serial code
    s_trunk_dt = distribute_tensor(
        s_trunk_global_host.to(device=manager.device, dtype=dtype),
        manager.device_mesh_subgroups,
        placements_s,
    ).requires_grad_(True)
    s_inputs_dt = distribute_tensor(
        s_inputs_global_host.to(device=manager.device, dtype=dtype),
        manager.device_mesh_subgroups,
        placements_s,
    ).requires_grad_(True)

    # Copies to verify inputs aren't modified
    times_dt_copy = times_dt.detach().clone()
    s_trunk_dt_copy = s_trunk_dt.detach().clone().requires_grad_(True)
    s_inputs_dt_copy = s_inputs_dt.detach().clone().requires_grad_(True)

    # Forward pass
    s_dt, normed_fourier_dt = module_dt(times_dt, s_trunk_dt, s_inputs_dt)

    # Ensure no input mutation
    assert_tensors_identical(times_dt_copy.to_local(), times_dt.to_local(), check_grad=False, check_grad_fn=False)
    assert_tensors_identical(s_trunk_dt_copy.to_local(), s_trunk_dt.to_local(), check_grad=False, check_grad_fn=False)
    assert_tensors_identical(s_inputs_dt_copy.to_local(), s_inputs_dt.to_local(), check_grad=False, check_grad_fn=False)

    # Forward compare
    torch.testing.assert_close(s_dt.full_tensor().cpu(), s_expected_global_host)
    if normed_fourier_expected_global_host is not None:
        assert normed_fourier_dt is not None
        torch.testing.assert_close(normed_fourier_dt.full_tensor().cpu(), normed_fourier_expected_global_host)
    else:
        assert normed_fourier_dt is None

    # Backward pass
    outputs = [s_dt]
    grad_outputs = [
        distribute_tensor(
            d_s_global_host.to(device=manager.device, dtype=dtype),
            manager.device_mesh_subgroups,
            placements_s,
        )
    ]
    if normed_fourier_dt is not None and d_normed_fourier_global_host is not None:
        outputs.append(normed_fourier_dt)
        grad_outputs.append(
            distribute_tensor(
                d_normed_fourier_global_host.to(device=manager.device, dtype=dtype),
                manager.device_mesh_subgroups,
                placements_times,
            )
        )
    torch.autograd.backward(outputs, grad_outputs)

    # Compare input gradients
    torch.testing.assert_close(s_trunk_dt.grad.full_tensor().cpu(), d_s_trunk_expected_global_host)
    torch.testing.assert_close(s_inputs_dt.grad.full_tensor().cpu(), d_s_inputs_expected_global_host)

    # Compare parameter gradients (skip frozen params like FourierEmbedding.proj)
    for name, param in module_dt.named_parameters():
        if not param.requires_grad:
            assert param.grad is None, f"Frozen parameter {name} should have no gradient"
            continue
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
        ((2, (2, 2)), True, "cuda", "ENV"),
    ],
    indirect=("setup_env",),
    ids=lambda x: f"dp:{x[0][0]}, cp:{x[0][1]}, specify_method:{x[1]}, device_type:{x[2]}, method_init:{x[3]}",
)
@pytest.mark.parametrize(
    "disable_times, v1_input_layout",
    [
        (False, False),  # V2 serial, times enabled, standard 2*token_s input_dim
        (False, True),  # V1 serial, times enabled, wider input_dim
        (True, False),  # V2 serial, times disabled (V2-only feature)
        # (True, True) is invalid — V1 does not support disable_times
    ],
    ids=["times:on-input:v2", "times:on-input:v1", "times:off"],
)
def test_single_conditioning(setup_env, disable_times: bool, v1_input_layout: bool):
    """Test SingleConditioning DTensor vs serial equivalence.

    Parametrized on the functional flags that differ between V1 and V2:
      - ``disable_times``:  whether fourier time embedding is skipped (V2-only)
      - ``v1_input_layout``: whether to use V1's wider input_dim (implies V1 serial class)

    The serial module class is inferred from these flags.
    """
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() < {world_size}")

    # Infer serial class from functional flags
    if v1_input_layout:
        assert not disable_times, "V1 does not support disable_times"
        serial_class, serial_class_tag = SingleConditioningBoltz1, "boltz1"
    else:
        serial_class, serial_class_tag = SingleConditioningBoltz2, "boltz2"

    dtype = torch.float64

    # Module dimensions
    token_s = 64
    dim_fourier = 32
    num_transitions = 2
    sigma_data = 1.0

    # Data dimensions
    size_ring = grid_group_sizes["cp"][0]
    B = 2 * grid_group_sizes["dp"]
    N = size_ring * 8

    val_init_min, val_init_max = -0.08, 0.08

    seed_by_rank(0, seed=42)

    # Build serial module kwargs — the input_dim and available kwargs differ by version
    if v1_input_layout:
        # V1: input_dim = 2 * token_s + 2 * const.num_tokens + 1 + len(const.pocket_contact_info)
        input_dim = 2 * token_s + 2 * const.num_tokens + 1 + len(const.pocket_contact_info)
        s_inputs_dim = input_dim - token_s  # s_trunk is token_s, s_inputs is the rest
        module_kwargs = {
            "sigma_data": sigma_data,
            "token_s": token_s,
            "dim_fourier": dim_fourier,
            "num_transitions": num_transitions,
        }
    else:
        # V2: input_dim = 2 * token_s
        s_inputs_dim = token_s  # s_trunk and s_inputs both have token_s dim
        module_kwargs = {
            "sigma_data": sigma_data,
            "token_s": token_s,
            "dim_fourier": dim_fourier,
            "num_transitions": num_transitions,
            "disable_times": disable_times,
        }

    # Create serial module — cast to target dtype BEFORE param init to prevent precision loss
    module_serial = serial_class(**module_kwargs)
    module_serial = module_serial.to(dtype=dtype).train()
    init_module_params_uniform(module_serial, low=val_init_min, high=val_init_max)
    layer_state_dict = module_serial.state_dict()

    # Create input tensors
    times_global = torch.empty(B, dtype=dtype)
    s_trunk_global = torch.empty(B, N, token_s, dtype=dtype, requires_grad=True)
    s_inputs_global = torch.empty(B, N, s_inputs_dim, dtype=dtype, requires_grad=True)
    init_tensors_uniform([times_global, s_trunk_global, s_inputs_global], low=val_init_min, high=val_init_max)

    # Serial forward pass — V1 uses keyword-only args, V2 uses positional
    if v1_input_layout:
        s_serial, normed_fourier_serial = module_serial(
            times=times_global, s_trunk=s_trunk_global, s_inputs=s_inputs_global
        )
    else:
        s_serial, normed_fourier_serial = module_serial(times_global, s_trunk_global, s_inputs_global)

    # Create upstream gradients
    d_s = torch.empty_like(s_serial)
    init_tensors_uniform([d_s], low=val_init_min, high=val_init_max)

    d_normed_fourier = None
    if normed_fourier_serial is not None:
        d_normed_fourier = torch.empty_like(normed_fourier_serial)
        init_tensors_uniform([d_normed_fourier], low=val_init_min, high=val_init_max)

    # Serial backward pass
    outputs = [s_serial]
    grad_outputs = [d_s]
    if normed_fourier_serial is not None and d_normed_fourier is not None:
        outputs.append(normed_fourier_serial)
        grad_outputs.append(d_normed_fourier)
    torch.autograd.backward(outputs, grad_outputs)

    # Collect expected results
    s_expected = s_serial.detach().clone().cpu()
    normed_fourier_expected = (
        normed_fourier_serial.detach().clone().cpu() if normed_fourier_serial is not None else None
    )
    d_s_trunk_expected = s_trunk_global.grad.detach().clone().cpu()
    d_s_inputs_expected = s_inputs_global.grad.detach().clone().cpu()

    expected_param_grads = {}
    for name, param in module_serial.named_parameters():
        if param.requires_grad and param.grad is not None:
            expected_param_grads[name] = param.grad.detach().clone().cpu()

    # Launch parallel test
    spawn_multiprocessing(
        parallel_assert_single_conditioning,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        serial_class_tag,
        layer_state_dict,
        dtype,
        times_global.detach().clone().cpu(),
        s_trunk_global.detach().clone().cpu(),
        s_inputs_global.detach().clone().cpu(),
        s_expected,
        normed_fourier_expected,
        d_s.detach().clone().cpu(),
        d_normed_fourier.detach().clone().cpu() if d_normed_fourier is not None else None,
        d_s_trunk_expected,
        d_s_inputs_expected,
        expected_param_grads,
        module_kwargs,
    )
