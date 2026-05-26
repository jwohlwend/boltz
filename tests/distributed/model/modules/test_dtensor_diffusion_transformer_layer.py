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

"""Tests for DTensor DiffusionTransformerLayer module.

Tests both Boltz-1x and Boltz-2 serial DiffusionTransformerLayer modules against
the unified DTensor implementation, verifying forward and backward equivalence.

Supports two attention modes:
- **Window-batched** (``use_ring_comm=False``): Uses ``AttentionPairBiasShardwise``.
  Inputs are 4D/5D window-batched tensors. Tests both V1 and V2 serial modules.
- **Ring attention** (``use_ring_comm=True``): Uses ``AttentionPairBias``.
  Inputs are 3D/4D token-level tensors. V2 only (token-level transformer use case).
"""

from functools import partial

import pytest
import torch
from torch.distributed.tensor import (
    DTensor,
    Replicate,
    Shard,
    distribute_tensor,
)

from boltz.distributed.comm import AttentionPairBiasComm
from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.layers.utils import convert_single_repr_window_batched_query_to_key
from boltz.distributed.model.modules.transformers import (
    DiffusionTransformerLayer as DistributedDiffusionTransformerLayer,
)
from boltz.model.modules.encoders import get_indexing_matrix as get_indexing_matrix_v1
from boltz.model.modules.encoders import single_to_keys as single_to_keys_v1
from boltz.model.modules.encodersv2 import get_indexing_matrix as get_indexing_matrix_v2
from boltz.model.modules.encodersv2 import single_to_keys as single_to_keys_v2
from boltz.model.modules.transformers import DiffusionTransformerLayer as SerialDTLBoltz1
from boltz.model.modules.transformersv2 import DiffusionTransformerLayer as SerialDTLBoltz2
from boltz.testing.utils import (
    SetModuleInfValues,
    assert_all_identical,
    assert_tensors_identical,
    get_param_by_key,
    init_module_params_uniform,
    init_tensors_uniform,
    seed_by_rank,
    spawn_multiprocessing,
)


def parallel_assert_diffusion_transformer_layer(
    rank: int,
    grid_group_sizes,
    device_type: str,
    backend: str,
    env_per_rank,
    serial_module_version: str,
    use_ring_comm: bool,
    dtype: torch.dtype,
    multiplicity: int,
    heads: int,
    dim: int,
    dim_single_cond: int,
    dim_pairwise: int,
    W: int,
    H: int,
    post_layer_norm: bool,
    layer_state_dict,
    a_global_host: torch.Tensor,
    s_global_host: torch.Tensor,
    z_global_host: torch.Tensor,
    mask_global_host: torch.Tensor,
    d_out_global_host: torch.Tensor,
    out_expected_global_host: torch.Tensor,
    d_a_expected_global_host: torch.Tensor,
    d_s_expected_global_host: torch.Tensor,
    d_z_expected_global_host: torch.Tensor,
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

    # Recreate serial module from state dict
    if serial_module_version == "boltz1":
        module_serial = SerialDTLBoltz1(
            heads=heads,
            dim=dim,
            dim_single_cond=dim_single_cond,
            dim_pairwise=dim_pairwise,
        )
    else:
        module_serial = SerialDTLBoltz2(
            heads=heads,
            dim=dim,
            dim_single_cond=dim_single_cond,
            post_layer_norm=post_layer_norm,
        )

    module_serial.load_state_dict(layer_state_dict)
    module_serial = module_serial.to(device=manager.device, dtype=dtype).train()
    module_serial.apply(SetModuleInfValues())

    # Create ring_comm if testing ring attention path
    ring_comm = None
    if use_ring_comm:
        ring_comm = AttentionPairBiasComm(
            manager.group["cp"],
            manager.layout_subgroups["cp"],
            manager.subgroups["cp"][0],
            manager.subgroups["cp"][1],
        )

    module = DistributedDiffusionTransformerLayer(
        diff_transformer_layer=module_serial,
        device_mesh=manager.device_mesh_subgroups,
        ring_comm=ring_comm,
    ).train()

    # Placements depend on attention mode
    placements_single = (Shard(0), Shard(1), Replicate())
    placements_pair = (Shard(0), Shard(1), Shard(2))

    if use_ring_comm:
        # Ring attention: a, s are 3D (B*M, N, D); z is 4D (B, N, N, heads); mask is 2D (B, N)
        a_dt = distribute_tensor(
            a_global_host.to(device=manager.device, dtype=dtype), manager.device_mesh_subgroups, placements_single
        ).requires_grad_(True)
        s_dt = distribute_tensor(
            s_global_host.to(device=manager.device, dtype=dtype), manager.device_mesh_subgroups, placements_single
        ).requires_grad_(True)
        z_dt = distribute_tensor(
            z_global_host.to(device=manager.device, dtype=dtype), manager.device_mesh_subgroups, placements_pair
        ).requires_grad_(True)
        mask_dt = distribute_tensor(
            mask_global_host.to(device=manager.device, dtype=dtype), manager.device_mesh_subgroups, placements_single
        ).requires_grad_(False)
    else:
        # Window-batched: a, s are 4D (B*M, K, W, D); z is 5D (B, K, W, H, D); mask is 3D (B, K, W)
        placements = placements_single
        a_dt = distribute_tensor(
            a_global_host.to(device=manager.device, dtype=dtype), manager.device_mesh_subgroups, placements
        ).requires_grad_(True)
        s_dt = distribute_tensor(
            s_global_host.to(device=manager.device, dtype=dtype), manager.device_mesh_subgroups, placements
        ).requires_grad_(True)
        z_dt = distribute_tensor(
            z_global_host.to(device=manager.device, dtype=dtype), manager.device_mesh_subgroups, placements
        ).requires_grad_(True)
        mask_dt = distribute_tensor(
            mask_global_host.to(device=manager.device, dtype=dtype), manager.device_mesh_subgroups, placements
        ).requires_grad_(False)

    # Copies to ensure inputs aren't modified in-place
    a_dt_copy = a_dt.detach().clone().requires_grad_(True)
    s_dt_copy = s_dt.detach().clone().requires_grad_(True)
    z_dt_copy = z_dt.detach().clone().requires_grad_(True)
    mask_dt_copy = mask_dt.detach().clone()

    # Forward pass
    if use_ring_comm:
        # Ring attention: no to_keys, use multiplicity
        out_dt: DTensor = module(
            a_dt,
            s_dt,
            z_dt,
            mask=mask_dt,
            multiplicity=multiplicity,
            layer_cache=None,
            pair_mask=None,
        )
    else:
        # Window-batched: to_keys converts query→key space
        to_keys = partial(convert_single_repr_window_batched_query_to_key, W=W, H=H)
        out_dt: DTensor = module(
            a_dt,
            s_dt,
            z_dt,
            mask=mask_dt,
            to_keys=to_keys,
            layer_cache=None,
            pair_mask=None,
        )

    # Ensure no input mutation
    assert_tensors_identical(a_dt_copy.to_local(), a_dt.to_local(), check_grad=False, check_grad_fn=False)
    assert_tensors_identical(s_dt_copy.to_local(), s_dt.to_local(), check_grad=False, check_grad_fn=False)
    assert_tensors_identical(z_dt_copy.to_local(), z_dt.to_local(), check_grad=False, check_grad_fn=False)
    assert_tensors_identical(mask_dt_copy.to_local(), mask_dt.to_local(), check_grad=False, check_grad_fn=False)

    # Forward compare
    out_expected_dt = distribute_tensor(
        out_expected_global_host.to(device=manager.device, dtype=dtype),
        manager.device_mesh_subgroups,
        out_dt.placements,
    )
    # Mask the forward comparison: fully-masked attention windows produce numerically
    # undefined softmax output (softmax of all -inf), so we only compare unmasked positions.
    _mask_local = mask_dt.to_local().unsqueeze(-1)
    if multiplicity > 1:
        _mask_local = _mask_local.repeat_interleave(multiplicity, 0)
    torch.testing.assert_close(out_dt.to_local() * _mask_local, out_expected_dt.to_local() * _mask_local)
    _mask_full = mask_dt.full_tensor().unsqueeze(-1)
    if multiplicity > 1:
        _mask_full = _mask_full.repeat_interleave(multiplicity, 0)
    torch.testing.assert_close(
        out_dt.full_tensor() * _mask_full,
        out_expected_global_host.to(device=manager.device, dtype=dtype) * _mask_full,
    )

    # Backward pass
    d_out_dt = distribute_tensor(
        d_out_global_host.to(device=manager.device, dtype=dtype),
        manager.device_mesh_subgroups,
        out_dt.placements,
    )
    d_out_dt_copy = d_out_dt.detach().clone()
    out_dt.backward(d_out_dt)
    assert_tensors_identical(d_out_dt_copy.to_local(), d_out_dt.to_local(), check_grad=False, check_grad_fn=False)

    # Input grad checks
    torch.testing.assert_close(a_dt.grad.full_tensor().cpu(), d_a_expected_global_host.to(dtype=dtype))
    torch.testing.assert_close(s_dt.grad.full_tensor().cpu(), d_s_expected_global_host.to(dtype=dtype))
    torch.testing.assert_close(z_dt.grad.full_tensor().cpu(), d_z_expected_global_host.to(dtype=dtype))

    # Parameter grads (gather full tensors)
    for name, grad_expected_global in expected_param_grads_global_host_dict.items():
        grad_param = get_param_by_key(module, name).grad
        assert grad_param is not None, f"Missing grad for param {name}"

        if isinstance(grad_param, DTensor):
            grad_global_host = grad_param.full_tensor().cpu()
            grad_to_check = grad_param.full_tensor()
        else:
            grad_global_host = grad_param.detach().cpu()
            grad_to_check = grad_param

        torch.testing.assert_close(grad_global_host, grad_expected_global.to(dtype=dtype))
        assert_all_identical(grad_to_check, manager.group["cp"])

    DistributedManager.cleanup()
    monkeypatch.undo()


def _create_serial_reference(
    serial_module_version: str,
    use_ring_comm: bool,
    heads: int,
    dim: int,
    dim_single_cond: int,
    dim_pairwise: int,
    post_layer_norm: bool,
    device_type: str,
    dtype: torch.dtype,
    val_init_min_max: tuple[float, float],
    B: int,
    N: int,
    K: int,
    W: int,
    H: int,
    multiplicity: int = 1,
):
    """Create serial module and compute reference forward/backward outputs.

    Supports two modes:
    - Window-batched (use_ring_comm=False): 4D/5D inputs, to_keys-based attention
    - Ring attention (use_ring_comm=True): 3D/4D token-level inputs, V2 only

    Returns tuple of (layer_state_dict, a, s, z, mask, d_out_global, out_expected_global,
                       d_a_expected_global, d_s_expected_global, d_z_expected_global,
                       expected_param_grads_global)
    """
    if serial_module_version == "boltz1":
        get_indexing_matrix = get_indexing_matrix_v1
        single_to_keys = single_to_keys_v1
        reference_module = SerialDTLBoltz1(
            heads=heads,
            dim=dim,
            dim_single_cond=dim_single_cond,
            dim_pairwise=dim_pairwise,
        ).to(device=device_type, dtype=dtype)
    else:
        get_indexing_matrix = get_indexing_matrix_v2
        single_to_keys = single_to_keys_v2
        reference_module = SerialDTLBoltz2(
            heads=heads,
            dim=dim,
            dim_single_cond=dim_single_cond,
            post_layer_norm=post_layer_norm,
        ).to(device=device_type, dtype=dtype)

    reference_module.train()

    # Keep values small to reduce numerical noise
    init_module_params_uniform(reference_module, low=val_init_min_max[0], high=val_init_min_max[1])
    reference_module.apply(SetModuleInfValues())
    layer_state_dict = reference_module.state_dict()

    if use_ring_comm:
        # ------------------------------------------------------------------
        # Ring attention path: 3D/4D token-level inputs
        # a, s: (B*M, N, D); z: (B, N, N, z_dim); mask: (B, N)
        # V1: z_dim = dim_pairwise (raw pair repr, proj_z computes bias)
        # V2: z_dim = heads (pre-computed bias)
        # ------------------------------------------------------------------
        z_last_dim = dim_pairwise if serial_module_version == "boltz1" else heads
        a = torch.empty((B * multiplicity, N, dim), device=device_type, dtype=dtype, requires_grad=True)
        s = torch.empty((B * multiplicity, N, dim_single_cond), device=device_type, dtype=dtype, requires_grad=True)
        z = torch.empty((B, N, N, z_last_dim), device=device_type, dtype=dtype, requires_grad=True)
        mask = torch.ones((B, N), device=device_type, dtype=dtype)
        mask[0, N // 2 :] = 0
        init_tensors_uniform([a, s, z], low=val_init_min_max[0], high=val_init_min_max[1])

        # Serial forward: no to_keys.
        # Pre-expand z and mask by multiplicity, then pass multiplicity=1 to avoid
        # double-expansion (serial AttentionPairBias internally repeats z by multiplicity).
        a_serial = a.detach().clone().requires_grad_(True)
        s_serial = s.detach().clone().requires_grad_(True)
        z_serial = z.detach().clone().requires_grad_(True)
        z_serial_mul = z_serial.repeat_interleave(multiplicity, 0)
        mask_mul = mask.repeat_interleave(multiplicity, 0)

        if serial_module_version == "boltz1":
            out_expected = reference_module(
                a_serial,
                s_serial,
                z_serial_mul,
                mask=mask_mul.detach(),
                to_keys=None,
                multiplicity=1,  # pre-applied above
                layer_cache=None,
            )
        else:
            out_expected = reference_module(
                a_serial,
                s_serial,
                bias=z_serial_mul,
                mask=mask_mul.detach(),
                to_keys=None,
                multiplicity=1,  # pre-applied above
            )

        d_out = torch.empty_like(out_expected)
        init_tensors_uniform([d_out], low=val_init_min_max[0], high=val_init_min_max[1])
        d_out = d_out * mask_mul.unsqueeze(-1)

        out_expected.backward(d_out)

        out_expected_global_host = out_expected.detach().cpu()
        d_a_expected_global_host = a_serial.grad.detach().cpu()
        d_s_expected_global_host = s_serial.grad.detach().cpu()
        d_z_expected_global_host = z_serial.grad.detach().cpu()
        d_out_global_host = d_out.detach().cpu()
    else:
        # ------------------------------------------------------------------
        # Window-batched path: 4D/5D inputs (V1 and V2)
        # a, s: (B*M, K, W, D); z: (B, K, W, H, z_dim); mask: (B, K, W)
        # ------------------------------------------------------------------
        z_last_dim = dim_pairwise if serial_module_version == "boltz1" else heads
        a = torch.empty((B * multiplicity, K, W, dim), device=device_type, dtype=dtype, requires_grad=True)
        s = torch.empty((B * multiplicity, K, W, dim_single_cond), device=device_type, dtype=dtype, requires_grad=True)
        z = torch.empty((B, K, W, H, z_last_dim), device=device_type, dtype=dtype, requires_grad=True)
        mask = torch.ones((B, K, W), device=device_type, dtype=dtype)
        mask[0, K // 2 :] = 0
        init_tensors_uniform([a, s, z], low=val_init_min_max[0], high=val_init_min_max[1])

        # Serial forward: flatten (B*M, K, ...) -> (B*M*K, ...)
        a_serial = a.detach().clone().requires_grad_(True)
        a_serial_flattened = a_serial.flatten(0, 1)  # (B*M*K, W, dim)
        s_serial = s.detach().clone().requires_grad_(True)
        s_serial_flattened = s_serial.flatten(0, 1)  # (B*M*K, W, dim_single_cond)
        z_serial = z.detach().clone().requires_grad_(True)

        BM = B * multiplicity
        indexing_matrix = get_indexing_matrix(K=K, W=W, H=H, device=device_type).to(dtype=dtype)

        def to_keys_serial(x: torch.Tensor) -> torch.Tensor:
            BMK, W_in, Dflat = x.shape
            assert BMK == BM * K
            assert W_in == W
            x_full = x.view(BM, K * W, Dflat)
            x_key = single_to_keys(x_full, indexing_matrix, W=W, H=H)  # (B*M, K, H, Dflat)
            return x_key.flatten(0, 1)  # (B*M*K, H, Dflat)

        # Pre-apply multiplicity to z and mask, then pass multiplicity=1
        z_serial_multiplex = z_serial.repeat_interleave(multiplicity, 0)
        z_serial_flattened = z_serial_multiplex.flatten(0, 1)
        mask_serial_multiplex = mask.detach().repeat_interleave(multiplicity, 0)
        mask_serial_flattened = mask_serial_multiplex.flatten(0, 1)

        if serial_module_version == "boltz1":
            mask_query = mask_serial_multiplex.view(BM * K, W)
            out_expected = reference_module(
                a_serial_flattened,
                s_serial_flattened,
                z_serial_flattened,
                mask=mask_query,
                to_keys=to_keys_serial,
                multiplicity=1,
                layer_cache=None,
            )
        else:
            out_expected = reference_module(
                a_serial_flattened,
                s_serial_flattened,
                bias=z_serial_flattened,
                mask=mask_serial_flattened,
                to_keys=to_keys_serial,
                multiplicity=1,
            )

        d_out = torch.empty_like(out_expected)
        init_tensors_uniform([d_out], low=val_init_min_max[0], high=val_init_min_max[1])
        d_out = d_out * mask_serial_flattened.unsqueeze(-1)

        out_expected.backward(d_out)

        # Unflatten (B*M*K, W, ...) -> (B*M, K, W, ...)
        out_expected_global_host = out_expected.detach().unflatten(0, a.shape[:2]).cpu()
        d_a_expected_global_host = a_serial.grad.detach().cpu()
        d_s_expected_global_host = s_serial.grad.detach().cpu()
        d_z_expected_global_host = z_serial.grad.detach().cpu()
        d_out_global_host = d_out.detach().unflatten(0, a.shape[:2]).cpu()

    expected_param_grads_global_host_dict = {
        name: param.grad.detach().cpu() for name, param in reference_module.named_parameters() if param.grad is not None
    }

    return (
        layer_state_dict,
        a.detach().cpu(),
        s.detach().cpu(),
        z.detach().cpu(),
        mask.detach().cpu(),
        d_out_global_host,
        out_expected_global_host,
        d_a_expected_global_host,
        d_s_expected_global_host,
        d_z_expected_global_host,
        expected_param_grads_global_host_dict,
    )


@pytest.mark.slow
@pytest.mark.parametrize("multiplicity", [1, 4], ids=lambda m: f"multiplicity:{m}")
@pytest.mark.parametrize(
    "setup_env",
    [
        ((1, (2, 2)), True, "cuda", "ENV"),
        ((2, (2, 2)), True, "cuda", "ENV"),
    ],
    indirect=("setup_env",),
    ids=lambda x: f"dp:{x[0][0]}, cp:{x[0][1]}, specify_method:{x[1]}, device_type:{x[2]}, method_init:{x[3]}",
)
@pytest.mark.parametrize("use_ring_comm", [False, True], ids=["wb", "ring"])
@pytest.mark.parametrize("serial_module_version", ["boltz1", "boltz2"])
def test_diffusion_transformer_layer(setup_env, multiplicity: int, use_ring_comm: bool, serial_module_version: str):
    """Test DiffusionTransformerLayer DTensor vs serial equivalence.

    Parametrized on:
    - ``use_ring_comm``: False for window-batched (AttentionPairBiasShardwise),
      True for ring attention (AttentionPairBias).
    - ``serial_module_version``: "boltz1" or "boltz2".
    """
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() < {world_size}")

    dtype = torch.float32
    seed = 42
    seed_by_rank(0, seed=seed)

    size_cp = grid_group_sizes["cp"][0]
    B = 2 * grid_group_sizes["dp"]
    N = 10 * size_cp * 32  # token count for ring attention (= K * W)
    K = 10 * size_cp  # windows, divisible by cp size
    W = 32  # queries per window (must be even)
    H = 128  # keys per window
    val_init_min_max = (-0.2, 0.2)

    dim = 32
    dim_single_cond = dim
    dim_pairwise = 32
    heads = 2
    post_layer_norm = False

    (
        layer_state_dict,
        a_host,
        s_host,
        z_host,
        mask_host,
        d_out_global_host,
        out_expected_global_host,
        d_a_expected_global_host,
        d_s_expected_global_host,
        d_z_expected_global_host,
        expected_param_grads_global_host_dict,
    ) = _create_serial_reference(
        serial_module_version=serial_module_version,
        use_ring_comm=use_ring_comm,
        heads=heads,
        dim=dim,
        dim_single_cond=dim_single_cond,
        dim_pairwise=dim_pairwise,
        post_layer_norm=post_layer_norm,
        device_type=device_type,
        dtype=dtype,
        val_init_min_max=val_init_min_max,
        B=B,
        N=N,
        K=K,
        W=W,
        H=H,
        multiplicity=multiplicity,
    )

    spawn_multiprocessing(
        parallel_assert_diffusion_transformer_layer,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        serial_module_version,
        use_ring_comm,
        dtype,
        multiplicity,
        heads,
        dim,
        dim_single_cond,
        dim_pairwise,
        W,
        H,
        post_layer_norm,
        layer_state_dict,
        a_host,
        s_host,
        z_host,
        mask_host,
        d_out_global_host,
        out_expected_global_host,
        d_a_expected_global_host,
        d_s_expected_global_host,
        d_z_expected_global_host,
        expected_param_grads_global_host_dict,
    )
