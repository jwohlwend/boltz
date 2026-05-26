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

"""Tests for DTensor PairwiseConditioning module.

Tests both Boltz-1x and Boltz-2 serial PairwiseConditioning modules against the unified
DTensor PairwiseConditioning implementation, verifying forward and backward equivalence.

Uses float64 with default tolerance for exact comparison.
"""

import pytest
import torch
from torch.distributed.tensor import Shard, distribute_tensor

from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.modules.encoders import PairwiseConditioning as DTensorPairwiseConditioning
from boltz.model.modules.encoders import PairwiseConditioning as PairwiseConditioningBoltz1
from boltz.model.modules.encodersv2 import PairwiseConditioning as PairwiseConditioningBoltz2
from boltz.testing.utils import (
    assert_tensors_identical,
    init_module_params_uniform,
    init_tensors_uniform,
    seed_by_rank,
    spawn_multiprocessing,
)


def parallel_assert_pairwise_conditioning(
    rank: int,
    grid_group_sizes,
    device_type: str,
    backend: str,
    env_per_rank,
    serial_module_version: str,
    layer_state_dict,
    dtype: torch.dtype,
    # Input tensors (global, on host)
    z_trunk_global_host: torch.Tensor,
    token_rel_pos_feats_global_host: torch.Tensor,
    # Expected outputs
    z_expected_global_host: torch.Tensor,
    # Upstream gradients
    d_z_global_host: torch.Tensor,
    # Expected input grads
    d_z_trunk_expected_global_host: torch.Tensor,
    d_token_rel_pos_feats_expected_global_host: torch.Tensor,
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
    if serial_module_version == "boltz1":
        module_serial = PairwiseConditioningBoltz1(**module_kwargs)
    else:
        module_serial = PairwiseConditioningBoltz2(**module_kwargs)
    module_serial = module_serial.to(device=manager.device, dtype=dtype)
    module_serial.load_state_dict(layer_state_dict)
    module_serial = module_serial.train()

    # Create DTensor module from serial
    module_dt = DTensorPairwiseConditioning(
        layer=module_serial,
        device_mesh=manager.device_mesh_subgroups,
    ).train()

    # Pair placements: shard along both token dimensions
    placements_pair = (Shard(0), Shard(1), Shard(2))

    # Distribute inputs
    z_trunk_dt = distribute_tensor(
        z_trunk_global_host.to(device=manager.device, dtype=dtype),
        manager.device_mesh_subgroups,
        placements_pair,
    ).requires_grad_(True)
    token_rel_pos_feats_dt = distribute_tensor(
        token_rel_pos_feats_global_host.to(device=manager.device, dtype=dtype),
        manager.device_mesh_subgroups,
        placements_pair,
    ).requires_grad_(True)

    # Copies to verify inputs aren't modified
    z_trunk_dt_copy = z_trunk_dt.detach().clone().requires_grad_(True)
    token_rel_pos_feats_dt_copy = token_rel_pos_feats_dt.detach().clone().requires_grad_(True)

    # Forward pass
    z_dt = module_dt(z_trunk_dt, token_rel_pos_feats_dt)

    # Ensure no input mutation
    assert_tensors_identical(z_trunk_dt_copy.to_local(), z_trunk_dt.to_local(), check_grad=False, check_grad_fn=False)
    assert_tensors_identical(
        token_rel_pos_feats_dt_copy.to_local(),
        token_rel_pos_feats_dt.to_local(),
        check_grad=False,
        check_grad_fn=False,
    )

    # Forward compare
    torch.testing.assert_close(z_dt.full_tensor().cpu(), z_expected_global_host)

    # Backward pass
    d_z_dt = distribute_tensor(
        d_z_global_host.to(device=manager.device, dtype=dtype),
        manager.device_mesh_subgroups,
        placements_pair,
    )
    z_dt.backward(d_z_dt)

    # Compare input gradients
    torch.testing.assert_close(z_trunk_dt.grad.full_tensor().cpu(), d_z_trunk_expected_global_host)
    torch.testing.assert_close(
        token_rel_pos_feats_dt.grad.full_tensor().cpu(), d_token_rel_pos_feats_expected_global_host
    )

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
        ((2, (2, 2)), True, "cuda", "ENV"),
    ],
    indirect=("setup_env",),
    ids=lambda x: f"dp:{x[0][0]}, cp:{x[0][1]}, specify_method:{x[1]}, device_type:{x[2]}, method_init:{x[3]}",
)
@pytest.mark.parametrize("serial_module_version", ["boltz1", "boltz2"])
def test_pairwise_conditioning(setup_env, serial_module_version: str):
    """Test PairwiseConditioning DTensor vs serial equivalence for both Boltz-1x and Boltz-2."""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() < {world_size}")

    dtype = torch.float64

    # Module dimensions
    token_z = 32
    dim_token_rel_pos_feats = 8
    num_transitions = 2

    # Data dimensions — N must be divisible by dp * cp
    size_ring = grid_group_sizes["cp"][0]
    B = 2 * grid_group_sizes["dp"]
    N = size_ring * 4  # tokens (keep small for pair O(N^2))

    val_init_min, val_init_max = -0.08, 0.08

    seed_by_rank(0, seed=42)

    module_kwargs = {
        "token_z": token_z,
        "dim_token_rel_pos_feats": dim_token_rel_pos_feats,
        "num_transitions": num_transitions,
    }

    if serial_module_version == "boltz1":
        module_serial = PairwiseConditioningBoltz1(**module_kwargs)
    else:
        module_serial = PairwiseConditioningBoltz2(**module_kwargs)

    # Cast to target dtype BEFORE param init to prevent precision loss
    module_serial = module_serial.to(dtype=dtype).train()
    init_module_params_uniform(module_serial, low=val_init_min, high=val_init_max)
    layer_state_dict = module_serial.state_dict()

    # Create input tensors
    z_trunk_global = torch.empty(B, N, N, token_z, dtype=dtype, requires_grad=True)
    token_rel_pos_feats_global = torch.empty(B, N, N, dim_token_rel_pos_feats, dtype=dtype, requires_grad=True)
    init_tensors_uniform([z_trunk_global, token_rel_pos_feats_global], low=val_init_min, high=val_init_max)

    # Serial forward pass
    z_serial = module_serial(z_trunk_global, token_rel_pos_feats_global)

    # Create upstream gradient
    d_z = torch.empty_like(z_serial)
    init_tensors_uniform([d_z], low=val_init_min, high=val_init_max)

    # Serial backward pass
    z_serial.backward(d_z)

    # Collect expected results
    z_expected = z_serial.detach().clone().cpu()
    d_z_trunk_expected = z_trunk_global.grad.detach().clone().cpu()
    d_token_rel_pos_feats_expected = token_rel_pos_feats_global.grad.detach().clone().cpu()

    expected_param_grads = {}
    for name, param in module_serial.named_parameters():
        if param.grad is not None:
            expected_param_grads[name] = param.grad.detach().clone().cpu()

    # Launch parallel test
    spawn_multiprocessing(
        parallel_assert_pairwise_conditioning,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        serial_module_version,
        layer_state_dict,
        dtype,
        z_trunk_global.detach().clone().cpu(),
        token_rel_pos_feats_global.detach().clone().cpu(),
        z_expected,
        d_z.detach().clone().cpu(),
        d_z_trunk_expected,
        d_token_rel_pos_feats_expected,
        expected_param_grads,
        module_kwargs,
    )
