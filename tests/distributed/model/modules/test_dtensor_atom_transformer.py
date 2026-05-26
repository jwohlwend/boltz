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

"""Tests for DTensor AtomTransformer module (window batching).

Tests both Boltz-1x and Boltz-2 serial AtomTransformer modules against the
unified DTensor AtomTransformer implementation, verifying forward and backward
equivalence.

"""

import pytest
import torch
from torch.distributed.tensor import (
    DTensor,
    Replicate,
    Shard,
    distribute_tensor,
)

from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.layers.utils import convert_single_repr_to_window_batched_key
from boltz.distributed.model.modules.transformers import (
    AtomTransformer as DistributedAtomTransformer,
)
from boltz.model.modules.encoders import get_indexing_matrix as get_indexing_matrix_v1
from boltz.model.modules.encoders import single_to_keys as single_to_keys_v1
from boltz.model.modules.encodersv2 import get_indexing_matrix as get_indexing_matrix_v2
from boltz.model.modules.encodersv2 import single_to_keys as single_to_keys_v2
from boltz.model.modules.transformers import AtomTransformer as SerialAtomTransformerBoltz1
from boltz.model.modules.transformersv2 import AtomTransformer as SerialAtomTransformerBoltz2
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


def parallel_assert_atom_transformer(
    rank: int,
    grid_group_sizes,
    device_type: str,
    backend: str,
    env_per_rank,
    serial_module_version: str,
    dtype: torch.dtype,
    multiplicity: int,
    depth: int,
    heads: int,
    dim: int,
    dim_single_cond: int,
    dim_pairwise: int,
    W: int,
    H: int,
    layer_state_dict,
    q_global_host: torch.Tensor,
    c_global_host: torch.Tensor,
    p_global_host: torch.Tensor,
    mask_global_host: torch.Tensor,
    d_out_global_host: torch.Tensor,
    out_expected_global_host: torch.Tensor,
    d_q_expected_global_host: torch.Tensor,
    d_c_expected_global_host: torch.Tensor,
    d_p_expected_global_host: torch.Tensor,
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
        module_serial = SerialAtomTransformerBoltz1(
            depth=depth,
            heads=heads,
            dim=dim,
            dim_single_cond=dim_single_cond,
            dim_pairwise=dim_pairwise,
            attn_window_queries=W,
            attn_window_keys=H,
        )
    else:
        module_serial = SerialAtomTransformerBoltz2(
            attn_window_queries=W,
            attn_window_keys=H,
            depth=depth,
            heads=heads,
            dim=dim,
            dim_single_cond=dim_single_cond,
        )

    module_serial.load_state_dict(layer_state_dict)
    module_serial = module_serial.to(device=manager.device, dtype=dtype).train()
    module_serial.apply(SetModuleInfValues())

    module = DistributedAtomTransformer(
        layer=module_serial,
        device_mesh=manager.device_mesh_subgroups,
    ).train()

    # AtomTransformer inputs are in single repr view: (B, N, D) where N = K * W
    placements_single = (
        (Shard(0), Shard(1), Replicate()) if manager.device_mesh_subgroups.ndim == 3 else (Shard(0), Shard(1))
    )
    # For pair representation p: (B, K, W, H, D)
    placements_pair = (
        (Shard(0), Shard(1), Replicate()) if manager.device_mesh_subgroups.ndim == 3 else (Shard(0), Shard(1))
    )

    q_dt = distribute_tensor(
        q_global_host.to(device=manager.device, dtype=dtype), manager.device_mesh_subgroups, placements_single
    ).requires_grad_(True)
    c_dt = distribute_tensor(
        c_global_host.to(device=manager.device, dtype=dtype), manager.device_mesh_subgroups, placements_single
    ).requires_grad_(True)
    p_dt = distribute_tensor(
        p_global_host.to(device=manager.device, dtype=dtype), manager.device_mesh_subgroups, placements_pair
    ).requires_grad_(True)
    mask_dt = distribute_tensor(
        mask_global_host.to(device=manager.device, dtype=dtype), manager.device_mesh_subgroups, placements_single
    ).requires_grad_(False)

    # Copies to ensure inputs aren't modified in-place
    q_dt_copy = q_dt.detach().clone().requires_grad_(True)
    c_dt_copy = c_dt.detach().clone().requires_grad_(True)
    p_dt_copy = p_dt.detach().clone().requires_grad_(True)
    mask_dt_copy = mask_dt.detach().clone()

    # multiplicity must be 1 for window batching
    out_dt: DTensor = module(
        q=q_dt,
        c=c_dt,
        p=p_dt,
        mask=mask_dt,
        multiplicity=1,
        model_cache=None,
        pair_mask=None,
    )

    # Ensure no input mutation
    assert_tensors_identical(q_dt_copy.to_local(), q_dt.to_local(), check_grad=False, check_grad_fn=False)
    assert_tensors_identical(c_dt_copy.to_local(), c_dt.to_local(), check_grad=False, check_grad_fn=False)
    assert_tensors_identical(p_dt_copy.to_local(), p_dt.to_local(), check_grad=False, check_grad_fn=False)
    assert_tensors_identical(mask_dt_copy.to_local(), mask_dt.to_local(), check_grad=False, check_grad_fn=False)

    # Forward compare (local shards + full gather)
    out_expected_dt = distribute_tensor(
        out_expected_global_host.to(device=manager.device, dtype=dtype),
        manager.device_mesh_subgroups,
        placements_single,
    )
    torch.testing.assert_close(out_dt.to_local(), out_expected_dt.to_local())
    torch.testing.assert_close(out_dt.full_tensor().cpu(), out_expected_global_host.to(dtype=dtype))

    # Backward compare
    d_out_dt = distribute_tensor(
        d_out_global_host.to(device=manager.device, dtype=dtype), manager.device_mesh_subgroups, placements_single
    )
    d_out_dt_copy = d_out_dt.detach().clone()
    out_dt.backward(d_out_dt)
    assert_tensors_identical(d_out_dt_copy.to_local(), d_out_dt.to_local(), check_grad=False, check_grad_fn=False)

    # Input grad checks
    d_q_expected_dt = distribute_tensor(
        d_q_expected_global_host.to(device=manager.device, dtype=dtype),
        manager.device_mesh_subgroups,
        placements_single,
    )
    d_c_expected_dt = distribute_tensor(
        d_c_expected_global_host.to(device=manager.device, dtype=dtype),
        manager.device_mesh_subgroups,
        placements_single,
    )
    torch.testing.assert_close(q_dt.grad.to_local(), d_q_expected_dt.to_local())
    torch.testing.assert_close(c_dt.grad.to_local(), d_c_expected_dt.to_local())

    p_grad_full = p_dt.grad.full_tensor()
    torch.testing.assert_close(p_grad_full.cpu(), d_p_expected_global_host.to(dtype=dtype))

    # Mask-based gradient checks: gradients at invalid (masked) positions should be zero
    mask_dt_local = (
        mask_dt.to_local().repeat_interleave(multiplicity, 0).unsqueeze(-1).bool()
    )  # (B*mult_local, N_local, 1)
    torch.testing.assert_close(q_dt.grad.to_local() * ~(mask_dt_local), torch.zeros_like(q_dt.grad.to_local()))
    torch.testing.assert_close(c_dt.grad.to_local() * ~(mask_dt_local), torch.zeros_like(c_dt.grad.to_local()))

    # For pair repr grad masking, compute key mask from single repr mask (B, N) -> (B, K, H)
    mask_key_dt = convert_single_repr_to_window_batched_key(mask_dt, W=W, H=H)  # DTensor (B, K, H)
    # (B, K_local, H) -> (B, K_local, 1, H, 1) for masking pair view (B, K, W, H, D)
    mask_key_local = mask_key_dt.to_local().unsqueeze(2).unsqueeze(-1).bool()
    torch.testing.assert_close(p_dt.grad.to_local() * ~(mask_key_local), torch.zeros_like(p_dt.grad.to_local()))

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
    depth: int,
    heads: int,
    dim: int,
    dim_single_cond: int,
    dim_pairwise: int,
    W: int,
    H: int,
    B: int,
    K: int,
    multiplicity: int,
    device_type: str,
    dtype: torch.dtype,
    val_init_min_max: tuple[float, float],
):
    """Create serial module and compute reference forward/backward outputs."""
    N = K * W

    if serial_module_version == "boltz1":
        get_indexing_matrix = get_indexing_matrix_v1
        single_to_keys = single_to_keys_v1
        reference_module = SerialAtomTransformerBoltz1(
            depth=depth,
            heads=heads,
            dim=dim,
            dim_single_cond=dim_single_cond,
            dim_pairwise=dim_pairwise,
            attn_window_queries=W,
            attn_window_keys=H,
        ).to(device=device_type, dtype=dtype)
    else:
        get_indexing_matrix = get_indexing_matrix_v2
        single_to_keys = single_to_keys_v2
        reference_module = SerialAtomTransformerBoltz2(
            attn_window_queries=W,
            attn_window_keys=H,
            depth=depth,
            heads=heads,
            dim=dim,
            dim_single_cond=dim_single_cond,
        ).to(device=device_type, dtype=dtype)

    reference_module.train()

    # Inputs in single repr view: (B*M, N, D) for q, c
    q = torch.empty((B * multiplicity, N, dim), device=device_type, dtype=dtype, requires_grad=True)
    c = torch.empty((B * multiplicity, N, dim_single_cond), device=device_type, dtype=dtype, requires_grad=True)
    mask = torch.ones((B, N), device=device_type, dtype=dtype)
    mask[0, N // 2 :] = 0  # mask out second half for first sample

    # Pair repr in window-batched view: (B, K, W, H, D_z)
    # Boltz-1: D_z = dim_pairwise (projected per layer)
    # Boltz-2: D_z = num_heads * depth (pre-computed bias, split across layers)
    z_last_dim = dim_pairwise if serial_module_version == "boltz1" else heads * depth
    p = torch.empty((B, K, W, H, z_last_dim), device=device_type, dtype=dtype, requires_grad=True)

    init_tensors_uniform([q, c, p], low=val_init_min_max[0], high=val_init_min_max[1])
    init_module_params_uniform(reference_module, low=val_init_min_max[0], high=val_init_min_max[1])
    reference_module.apply(SetModuleInfValues())

    layer_state_dict = reference_module.state_dict()

    # Serial forward
    q_serial = q.detach().clone().requires_grad_(True)
    c_serial = c.detach().clone().requires_grad_(True)
    p_serial = p.detach().clone().requires_grad_(True)
    mask_multiplexed = mask.repeat_interleave(multiplicity, 0)

    # to_keys for serial AtomTransformer
    indexing_matrix = get_indexing_matrix(K=K, W=W, H=H, device=device_type).to(dtype=dtype)

    def to_keys_serial(x: torch.Tensor) -> torch.Tensor:
        return single_to_keys(x, indexing_matrix, W=W, H=H)

    if serial_module_version == "boltz1":
        out_expected = reference_module(
            q=q_serial,
            c=c_serial,
            p=p_serial,
            mask=mask_multiplexed,
            multiplicity=multiplicity,
            to_keys=to_keys_serial,
            model_cache=None,
        )
    else:
        out_expected = reference_module(
            q=q_serial,
            c=c_serial,
            bias=p_serial,
            to_keys=to_keys_serial,
            mask=mask_multiplexed,
            multiplicity=multiplicity,
        )

    d_out = torch.empty_like(out_expected)
    init_tensors_uniform([d_out], low=val_init_min_max[0], high=val_init_min_max[1])
    d_out = d_out * mask_multiplexed.unsqueeze(-1)

    out_expected.backward(d_out)

    return (
        layer_state_dict,
        q.detach().cpu(),
        c.detach().cpu(),
        p.detach().cpu(),
        mask.detach().cpu(),
        d_out.detach().cpu(),
        out_expected.detach().cpu(),
        q_serial.grad.detach().cpu(),
        c_serial.grad.detach().cpu(),
        p_serial.grad.detach().cpu(),
        {
            name: param.grad.detach().cpu()
            for name, param in reference_module.named_parameters()
            if param.grad is not None
        },
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
@pytest.mark.parametrize("serial_module_version", ["boltz1", "boltz2"])
def test_atom_transformer(setup_env, multiplicity: int, serial_module_version: str):
    """Test AtomTransformer DTensor vs serial equivalence (window batching)."""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() < {world_size}")

    dtype = torch.float32
    seed = 42
    seed_by_rank(0, seed=seed)

    B = 2 * grid_group_sizes["dp"]
    K = 10 * grid_group_sizes["cp"][0]  # windows, divisible by cp size
    W = 32
    H = 128
    val_init_min_max = (-0.2, 0.2)

    dim = 32
    dim_single_cond = dim
    dim_pairwise = 32
    heads = 2
    depth = 2

    (
        layer_state_dict,
        q_host,
        c_host,
        p_host,
        mask_host,
        d_out_host,
        out_expected_host,
        d_q_expected_host,
        d_c_expected_host,
        d_p_expected_host,
        expected_param_grads_host,
    ) = _create_serial_reference(
        serial_module_version=serial_module_version,
        depth=depth,
        heads=heads,
        dim=dim,
        dim_single_cond=dim_single_cond,
        dim_pairwise=dim_pairwise,
        W=W,
        H=H,
        B=B,
        K=K,
        multiplicity=multiplicity,
        device_type=device_type,
        dtype=dtype,
        val_init_min_max=val_init_min_max,
    )

    spawn_multiprocessing(
        parallel_assert_atom_transformer,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        serial_module_version,
        dtype,
        multiplicity,
        depth,
        heads,
        dim,
        dim_single_cond,
        dim_pairwise,
        W,
        H,
        layer_state_dict,
        q_host,
        c_host,
        p_host,
        mask_host,
        d_out_host,
        out_expected_host,
        d_q_expected_host,
        d_c_expected_host,
        d_p_expected_host,
        expected_param_grads_host,
    )


# ======================================================================
# Test 2: AtomTransformer under autocast bf16 (dtype-only comparison)
# ======================================================================


def parallel_assert_atom_transformer_autocast_bf16(
    rank: int,
    grid_group_sizes,
    device_type: str,
    backend: str,
    env_per_rank,
    depth: int,
    heads: int,
    dim: int,
    dim_single_cond: int,
    W: int,
    H: int,
    layer_state_dict,
    q_global_host: torch.Tensor,
    c_global_host: torch.Tensor,
    p_global_host: torch.Tensor,
    mask_global_host: torch.Tensor,
    serial_output_dtype: torch.dtype,
    serial_grad_dtypes: dict[str, torch.dtype],
    serial_param_grad_dtypes: dict[str, torch.dtype],
):
    """Parallel worker for bf16 autocast dtype test on AtomTransformer."""
    monkeypatch = pytest.MonkeyPatch()
    if env_per_rank is not None:
        for var_name, value in env_per_rank.items():
            if value == "<INPUT_RANK>":
                monkeypatch.setenv(var_name, f"{rank}")
                continue
            monkeypatch.setenv(var_name, value)

    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    manager = DistributedManager()

    dtype = torch.float32

    module_serial = SerialAtomTransformerBoltz2(
        attn_window_queries=W,
        attn_window_keys=H,
        depth=depth,
        heads=heads,
        dim=dim,
        dim_single_cond=dim_single_cond,
    )
    module_serial.load_state_dict(layer_state_dict)
    module_serial = module_serial.to(device=manager.device, dtype=dtype).train()
    module_serial.apply(SetModuleInfValues())

    module = DistributedAtomTransformer(
        layer=module_serial,
        device_mesh=manager.device_mesh_subgroups,
    ).train()

    placements_single = (
        (Shard(0), Shard(1), Replicate()) if manager.device_mesh_subgroups.ndim == 3 else (Shard(0), Shard(1))
    )
    placements_pair = (
        (Shard(0), Shard(1), Replicate()) if manager.device_mesh_subgroups.ndim == 3 else (Shard(0), Shard(1))
    )

    q_dt = distribute_tensor(
        q_global_host.to(device=manager.device, dtype=dtype),
        manager.device_mesh_subgroups,
        placements_single,
    ).requires_grad_(True)
    c_dt = distribute_tensor(
        c_global_host.to(device=manager.device, dtype=dtype),
        manager.device_mesh_subgroups,
        placements_single,
    ).requires_grad_(True)
    p_dt = distribute_tensor(
        p_global_host.to(device=manager.device, dtype=dtype),
        manager.device_mesh_subgroups,
        placements_pair,
    ).requires_grad_(True)
    mask_dt = distribute_tensor(
        mask_global_host.to(device=manager.device, dtype=dtype),
        manager.device_mesh_subgroups,
        placements_single,
    ).requires_grad_(False)

    with torch.autocast("cuda", dtype=torch.bfloat16):
        out_dt = module(q=q_dt, c=c_dt, p=p_dt, mask=mask_dt, multiplicity=1, model_cache=None, pair_mask=None)

    torch.autograd.backward([out_dt], [torch.ones_like(out_dt)])

    assert (
        out_dt.dtype == serial_output_dtype
    ), f"out dtype mismatch: DTensor {out_dt.dtype} vs serial {serial_output_dtype}"

    for name, dt_tensor in [("q", q_dt), ("c", c_dt)]:
        assert dt_tensor.grad is not None, f"{name} grad is None"
        assert (
            dt_tensor.grad.dtype == serial_grad_dtypes[name]
        ), f"{name} grad dtype mismatch: DTensor {dt_tensor.grad.dtype} vs serial {serial_grad_dtypes[name]}"

    for name, param in module.named_parameters():
        if name in serial_param_grad_dtypes and param.grad is not None:
            grad_dtype = param.grad.full_tensor().dtype if hasattr(param.grad, "full_tensor") else param.grad.dtype
            assert (
                grad_dtype == serial_param_grad_dtypes[name]
            ), f"param '{name}' grad dtype mismatch: DTensor {grad_dtype} vs serial {serial_param_grad_dtypes[name]}"

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.slow
@pytest.mark.parametrize(
    "setup_env",
    [((1, (1, 1)), True, "cuda", "ENV")],
    indirect=("setup_env",),
    ids=lambda x: f"dp:{x[0][0]}, cp:{x[0][1]}, device_type:{x[2]}",
)
def test_atom_transformer_autocast_bf16(setup_env):
    """Test DTensor AtomTransformer output dtypes under autocast bf16."""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() < {world_size}")

    seed = 42
    seed_by_rank(0, seed=seed)
    dtype = torch.float32

    B = 1
    K = 10
    W = 32
    H = 128
    val_init_min_max = (-0.2, 0.2)
    N = K * W
    multiplicity = 1

    dim = 32
    dim_single_cond = dim
    heads = 2
    depth = 2
    z_last_dim = heads * depth

    reference_module = SerialAtomTransformerBoltz2(
        attn_window_queries=W,
        attn_window_keys=H,
        depth=depth,
        heads=heads,
        dim=dim,
        dim_single_cond=dim_single_cond,
    ).to(device=device_type, dtype=dtype)
    reference_module.train()

    q = torch.empty((B * multiplicity, N, dim), device=device_type, dtype=dtype, requires_grad=True)
    c = torch.empty((B * multiplicity, N, dim_single_cond), device=device_type, dtype=dtype, requires_grad=True)
    mask = torch.ones((B, N), device=device_type, dtype=dtype)
    p = torch.empty((B, K, W, H, z_last_dim), device=device_type, dtype=dtype, requires_grad=True)
    init_tensors_uniform([q, c, p], low=val_init_min_max[0], high=val_init_min_max[1])
    init_module_params_uniform(reference_module, low=val_init_min_max[0], high=val_init_min_max[1])
    reference_module.apply(SetModuleInfValues())
    layer_state_dict = reference_module.state_dict()

    q_serial = q.detach().clone().requires_grad_(True)
    c_serial = c.detach().clone().requires_grad_(True)
    p_serial = p.detach().clone().requires_grad_(True)

    indexing_matrix = get_indexing_matrix_v2(K=K, W=W, H=H, device=device_type).to(dtype=dtype)

    with torch.autocast("cuda", dtype=torch.bfloat16):
        out_serial = reference_module(
            q=q_serial,
            c=c_serial,
            bias=p_serial,
            to_keys=lambda x: single_to_keys_v2(x, indexing_matrix, W=W, H=H),
            mask=mask.repeat_interleave(multiplicity, 0),
            multiplicity=multiplicity,
        )

    torch.autograd.backward([out_serial], [torch.ones_like(out_serial)])

    serial_output_dtype = out_serial.dtype
    serial_grad_dtypes = {"q": q_serial.grad.dtype, "c": c_serial.grad.dtype}
    serial_param_grad_dtypes = {
        name: param.grad.dtype for name, param in reference_module.named_parameters() if param.grad is not None
    }

    spawn_multiprocessing(
        parallel_assert_atom_transformer_autocast_bf16,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        depth,
        heads,
        dim,
        dim_single_cond,
        W,
        H,
        {k: v.detach().cpu() for k, v in layer_state_dict.items()},
        q.detach().cpu(),
        c.detach().cpu(),
        p.detach().cpu(),
        mask.detach().cpu(),
        serial_output_dtype,
        serial_grad_dtypes,
        serial_param_grad_dtypes,
    )
