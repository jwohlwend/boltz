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

"""Tests for DTensor AtomAttentionDecoder with window batching."""

from functools import partial

import pytest
import torch
from torch.distributed.tensor import distribute_tensor

from boltz.distributed.data.feature.featurizer import pack_atom_features
from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.layers.flatten_and_unflatten import shardwise_unflatten_sharded
from boltz.distributed.model.layers.utils import convert_single_repr_to_window_batched_key
from boltz.distributed.model.modules.encoders import (
    AtomAttentionDecoder as DistributedAtomAttentionDecoder,
)
from boltz.model.modules.encoders import AtomAttentionDecoder as SerialAtomAttentionDecoderBoltz1
from boltz.model.modules.encoders import get_indexing_matrix as get_indexing_matrix_v1
from boltz.model.modules.encoders import single_to_keys as single_to_keys_v1
from boltz.model.modules.encodersv2 import AtomAttentionDecoder as SerialAtomAttentionDecoderBoltz2
from boltz.model.modules.encodersv2 import get_indexing_matrix as get_indexing_matrix_v2
from boltz.model.modules.encodersv2 import single_to_keys as single_to_keys_v2
from boltz.testing.utils import (
    SetModuleInfValues,
    assert_all_identical,
    assert_tensors_close_with_pad,
    distribute_atom_features,
    get_feature_placements,
    get_param_by_key,
    init_module_params_uniform,
    init_tensors_uniform,
    pad_or_shrink_to_length,
    random_features,
    seed_by_rank,
    spawn_multiprocessing,
)

# Subset of keys needed for AtomAttentionDecoder window batching test
_selected_atom_keys = {
    "atom_pad_mask",
    "atom_to_token",
    "atom_counts_per_token",  # Required by pad_and_scatter_atom_features_dtensor
}

_placements = get_feature_placements(
    token_keys=set(),
    msa_keys=set(),
    atom_keys=_selected_atom_keys,
    model_io_keys=set(),
    model_io_fp32_keys=set(),
)
_placements_single = _placements["single"]
_placements_cp_atom_features = _placements["cp_atom_features"]
_placements_atom_features = _placements["atom_features"]


def parallel_assert_atom_attention_decoder_wb(
    rank: int,
    grid_group_sizes,
    device_type: str,
    backend: str,
    env_per_rank,
    serial_module_version: str,
    dtype: torch.dtype,
    multiplicity: int,
    atom_s: int,
    atom_z: int,
    token_s: int,
    W: int,
    H: int,
    atom_decoder_depth: int,
    atom_decoder_heads: int,
    layer_state_dict,
    feats_global_host: dict[str, torch.Tensor],
    a_global_host: torch.Tensor,
    q_global_host: torch.Tensor,
    c_global_host: torch.Tensor,
    p_global_host: torch.Tensor,
    d_r_update_global_host: torch.Tensor,
    r_update_expected_global_host: torch.Tensor,
    d_a_expected_global_host: torch.Tensor,
    d_q_expected_global_host: torch.Tensor,
    d_c_expected_global_host: torch.Tensor,
    d_p_expected_global_host: torch.Tensor,
    expected_param_grads_global_host_dict: dict[str, torch.Tensor],
):
    """Parallel worker function for testing DTensor AtomAttentionDecoder with window batching."""
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
        module_serial = SerialAtomAttentionDecoderBoltz1(
            atom_s=atom_s,
            atom_z=atom_z,
            token_s=token_s,
            attn_window_queries=W,
            attn_window_keys=H,
            atom_decoder_depth=atom_decoder_depth,
            atom_decoder_heads=atom_decoder_heads,
        )
    else:
        module_serial = SerialAtomAttentionDecoderBoltz2(
            atom_s=atom_s,
            token_s=token_s,
            attn_window_queries=W,
            attn_window_keys=H,
            atom_decoder_depth=atom_decoder_depth,
            atom_decoder_heads=atom_decoder_heads,
        )
    module_serial.load_state_dict(layer_state_dict)
    module_serial = module_serial.to(device=manager.device, dtype=dtype).train()
    module_serial.apply(SetModuleInfValues())

    # Create distributed module
    module = DistributedAtomAttentionDecoder(
        layer=module_serial,
        device_mesh=manager.device_mesh_subgroups,
    ).train()

    # Get global masks
    atom_pad_mask_global = feats_global_host["atom_pad_mask"].to(device=manager.device, dtype=torch.bool)
    atom_pad_mask_expanded_global = atom_pad_mask_global.unsqueeze(-1)
    atom_pad_mask_expanded_global_mul = atom_pad_mask_expanded_global.repeat_interleave(multiplicity, dim=0)

    # Distribute atom features
    inputs_atom = {
        k: v.to(dtype=dtype if v.dtype.is_floating_point else v.dtype)
        for k, v in feats_global_host.items()
        if k in _placements_cp_atom_features
    }
    feats_dt = distribute_atom_features(
        inputs_atom,
        _placements_cp_atom_features,
        _placements_atom_features,
        manager.device_mesh_subgroups,
        manager.group["cp"],
    )

    # Pack atom features
    feats_dt_packed = pack_atom_features(feats_dt, set(feats_dt.keys()), W)
    N_atoms_packed = feats_dt_packed["atom_pad_mask"].shape[1]
    K_packed = N_atoms_packed // W

    # Distribute input tensors
    a_dt = distribute_tensor(
        a_global_host.to(device=manager.device, dtype=dtype),
        manager.device_mesh_subgroups,
        _placements_single,
    ).requires_grad_(True)

    q_adjusted = pad_or_shrink_to_length(
        q_global_host.to(device=manager.device, dtype=dtype), axis=1, target_length=N_atoms_packed
    )
    q_dt = distribute_tensor(q_adjusted, manager.device_mesh_subgroups, _placements_single).requires_grad_(True)

    c_adjusted = pad_or_shrink_to_length(
        c_global_host.to(device=manager.device, dtype=dtype), axis=1, target_length=N_atoms_packed
    )
    c_dt = distribute_tensor(c_adjusted, manager.device_mesh_subgroups, _placements_single).requires_grad_(True)

    p_adjusted = pad_or_shrink_to_length(
        p_global_host.to(device=manager.device, dtype=dtype), axis=1, target_length=K_packed
    )
    p_dt = distribute_tensor(p_adjusted, manager.device_mesh_subgroups, _placements_single).requires_grad_(True)

    # Forward pass
    r_update_dt = module(a=a_dt, q=q_dt, c=c_dt, p=p_dt, feats=feats_dt_packed, multiplicity=multiplicity)

    # Forward comparison
    r_update_expected_device = r_update_expected_global_host.to(device=manager.device, dtype=dtype)
    r_update_dt_full = r_update_dt.full_tensor()

    mask_dt_full = feats_dt_packed["atom_pad_mask"].full_tensor()
    mask_dt_full_mul = mask_dt_full.repeat_interleave(multiplicity, dim=0)
    mask_dt_full_mul_expanded = mask_dt_full_mul.unsqueeze(-1)

    assert_tensors_close_with_pad(
        r_update_dt_full * mask_dt_full_mul_expanded,
        r_update_expected_device * atom_pad_mask_expanded_global_mul,
        axis=1,
        pad_val=0,
    )

    # Backward pass
    d_r_update_padded = pad_or_shrink_to_length(
        d_r_update_global_host.to(device=manager.device, dtype=dtype), axis=1, target_length=r_update_dt.shape[1]
    )
    d_r_update_dtensor = distribute_tensor(d_r_update_padded, manager.device_mesh_subgroups, r_update_dt.placements)
    torch.autograd.backward([r_update_dt], [d_r_update_dtensor])

    # Check a gradient
    d_a_expected_device = d_a_expected_global_host.to(device=manager.device, dtype=dtype)
    d_a_dt_full = a_dt.grad.full_tensor()
    torch.testing.assert_close(d_a_dt_full, d_a_expected_device)

    # Check q gradient
    d_q_expected_device = d_q_expected_global_host.to(device=manager.device, dtype=dtype)
    q_grad_full = q_dt.grad.full_tensor()
    assert_tensors_close_with_pad(
        q_grad_full * mask_dt_full_mul_expanded,
        d_q_expected_device * atom_pad_mask_expanded_global_mul,
        axis=1,
        pad_val=0,
    )

    # Check c gradient
    d_c_expected_device = d_c_expected_global_host.to(device=manager.device, dtype=dtype)
    c_grad_full = c_dt.grad.full_tensor()
    assert_tensors_close_with_pad(
        c_grad_full * mask_dt_full_mul_expanded,
        d_c_expected_device * atom_pad_mask_expanded_global_mul,
        axis=1,
        pad_val=0,
    )

    # Check p gradient
    d_p_expected_device = d_p_expected_global_host.to(device=manager.device, dtype=dtype)
    p_grad_full = p_dt.grad.full_tensor()

    mask_dt_query = shardwise_unflatten_sharded(
        feats_dt_packed["atom_pad_mask"], axis=1, sizes=(feats_dt_packed["atom_pad_mask"].shape[1] // W, W)
    )
    mask_dt_query_full = mask_dt_query.full_tensor()
    mask_dt_query_full_expanded = mask_dt_query_full[:, :, :, None, None]
    mask_dt_key = convert_single_repr_to_window_batched_key(feats_dt_packed["atom_pad_mask"], W, H)
    mask_dt_key_full = mask_dt_key.full_tensor()
    mask_dt_key_full_expanded = mask_dt_key_full[:, :, None, :, None]
    mask_dt_pair_full_expanded = mask_dt_query_full_expanded * mask_dt_key_full_expanded

    N_atoms_serial = feats_global_host["atom_pad_mask"].shape[1]
    K_serial = N_atoms_serial // W
    if serial_module_version == "boltz1":
        index_matrix = get_indexing_matrix_v1(K_serial, W, H, manager.device).to(dtype=dtype)
        to_keys_fn = partial(single_to_keys_v1, indexing_matrix=index_matrix, W=W, H=H)
    else:
        index_matrix = get_indexing_matrix_v2(K_serial, W, H, manager.device).to(dtype=dtype)
        to_keys_fn = partial(single_to_keys_v2, indexing_matrix=index_matrix, W=W, H=H)

    mask_key_expected = to_keys_fn(
        feats_global_host["atom_pad_mask"].to(device=manager.device, dtype=dtype).unsqueeze(-1)
    )
    mask_key_expected_expanded = mask_key_expected[:, :, None, :, :]
    mask_query_expected_expanded = atom_pad_mask_expanded_global.unflatten(
        1, (atom_pad_mask_expanded_global.shape[1] // W, W)
    )[:, :, :, None, :]
    mask_pair_expected_expanded = mask_query_expected_expanded * mask_key_expected_expanded

    assert_tensors_close_with_pad(
        p_grad_full * mask_dt_pair_full_expanded,
        d_p_expected_device * mask_pair_expected_expanded,
        axis=1,
        pad_val=0,
    )

    # Parameter grads comparison
    for name, grad_expected_global in expected_param_grads_global_host_dict.items():
        grad_param = get_param_by_key(module, name).grad
        assert grad_param is not None, f"Missing grad for param {name}"

        if hasattr(grad_param, "full_tensor"):
            grad_global_host = grad_param.full_tensor().cpu()
            grad_to_check = grad_param.full_tensor()
        else:
            grad_global_host = grad_param.detach().cpu()
            grad_to_check = grad_param

        torch.testing.assert_close(grad_global_host, grad_expected_global.to(dtype=dtype))
        assert_all_identical(grad_to_check, manager.group["cp"])

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.slow
@pytest.mark.parametrize(
    "setup_env, dtype, multiplicity",
    (
        params_test := [
            (((1, (2, 2)), True, "cuda", "ENV"), torch.float32, 1),
            (((2, (2, 2)), True, "cuda", "ENV"), torch.float32, 4),
        ]
    ),
    indirect=["setup_env"],
    ids=[
        f"dp:{x[0][0][0]}, cp:{x[0][0][1]}, device_type:{x[0][2]}, dtype:{x[1]}, multiplicity:{x[2]}"
        for x in params_test
    ],
)
@pytest.mark.parametrize("serial_module_version", ["boltz1", "boltz2"])
def test_atom_attention_decoder_window_batching(setup_env, dtype, multiplicity, serial_module_version):
    """Test DTensor AtomAttentionDecoder with window batching for V1 and V2."""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() < {world_size}")

    seed = 42
    seed_by_rank(0, seed=seed)

    size_cp = grid_group_sizes["cp"][0]
    B = 1 * grid_group_sizes["dp"]

    W = 32
    H = 128
    val_init_min_max = (-0.1, 0.1)

    n_atoms_per_token_min = 8
    n_atoms_per_token_max = 20
    N_tokens = 100 * size_cp
    N_atoms_raw = N_tokens * n_atoms_per_token_max
    N_atoms = ((N_atoms_raw + W - 1) // W) * W
    N_msa = 1

    atom_s = 8
    atom_z = 8
    token_s = 2
    atom_decoder_depth = 2
    atom_decoder_heads = 2

    # For Boltz-2, p last dim = num_heads * depth (pre-computed bias)
    # For Boltz-1, p last dim = atom_z (pair representation)
    p_last_dim = atom_z if serial_module_version == "boltz1" else atom_decoder_heads * atom_decoder_depth

    selected_keys = list(_selected_atom_keys)

    assert N_tokens % size_cp == 0

    feats = random_features(
        size_batch=B,
        n_tokens=N_tokens,
        n_atoms=N_atoms,
        n_msa=N_msa,
        atom_counts_per_token_range=(n_atoms_per_token_min, n_atoms_per_token_max),
        device=torch.device(device_type),
        float_value_range=val_init_min_max,
        selected_keys=selected_keys,
    )
    feats = {k: v.to(dtype=dtype) if v.dtype == torch.float64 else v for k, v in feats.items()}

    N_atoms_actual = feats["atom_pad_mask"].shape[1]
    assert N_atoms_actual % W == 0
    K = N_atoms_actual // W

    # Generate input tensors
    a = torch.empty((B * multiplicity, N_tokens, token_s * 2), device=device_type, dtype=dtype, requires_grad=True)
    q = torch.empty((B * multiplicity, N_atoms_actual, atom_s), device=device_type, dtype=dtype, requires_grad=True)
    c = torch.empty((B * multiplicity, N_atoms_actual, atom_s), device=device_type, dtype=dtype, requires_grad=True)
    p = torch.empty((B, K, W, H, p_last_dim), device=device_type, dtype=dtype, requires_grad=True)
    init_tensors_uniform([a, q, c, p], low=val_init_min_max[0], high=val_init_min_max[1])

    # Build serial module
    if serial_module_version == "boltz1":
        get_indexing_matrix = get_indexing_matrix_v1
        single_to_keys = single_to_keys_v1
        reference_module = SerialAtomAttentionDecoderBoltz1(
            atom_s=atom_s,
            atom_z=atom_z,
            token_s=token_s,
            attn_window_queries=W,
            attn_window_keys=H,
            atom_decoder_depth=atom_decoder_depth,
            atom_decoder_heads=atom_decoder_heads,
        ).to(device=device_type, dtype=dtype)
    else:
        get_indexing_matrix = get_indexing_matrix_v2
        single_to_keys = single_to_keys_v2
        reference_module = SerialAtomAttentionDecoderBoltz2(
            atom_s=atom_s,
            token_s=token_s,
            attn_window_queries=W,
            attn_window_keys=H,
            atom_decoder_depth=atom_decoder_depth,
            atom_decoder_heads=atom_decoder_heads,
        ).to(device=device_type, dtype=dtype)

    reference_module.train()
    init_module_params_uniform(reference_module, low=val_init_min_max[0], high=val_init_min_max[1])
    reference_module.apply(SetModuleInfValues())
    layer_state_dict = reference_module.state_dict()

    # Serial forward
    feats_serial = {k: v.detach().clone() for k, v in feats.items()}
    a_serial = a.detach().clone().requires_grad_(True)
    q_serial = q.detach().clone().requires_grad_(True)
    c_serial = c.detach().clone().requires_grad_(True)
    p_serial = p.detach().clone().requires_grad_(True)

    index_matrix = get_indexing_matrix(K, W, H, device_type).to(dtype=dtype)
    to_keys = partial(single_to_keys, indexing_matrix=index_matrix, W=W, H=H)

    if serial_module_version == "boltz1":
        r_update_expected = reference_module(
            a=a_serial,
            q=q_serial,
            c=c_serial,
            p=p_serial,
            feats=feats_serial,
            to_keys=to_keys,
            multiplicity=multiplicity,
            model_cache=None,
        )
    else:
        r_update_expected = reference_module(
            a=a_serial,
            q=q_serial,
            c=c_serial,
            atom_dec_bias=p_serial,
            feats=feats_serial,
            to_keys=to_keys,
            multiplicity=multiplicity,
        )

    # Upstream gradient
    d_r_update = torch.empty_like(r_update_expected)
    init_tensors_uniform([d_r_update], low=val_init_min_max[0], high=val_init_min_max[1])
    d_r_update = d_r_update * feats_serial["atom_pad_mask"].unsqueeze(-1).repeat_interleave(multiplicity, dim=0)

    torch.autograd.backward([r_update_expected], [d_r_update])

    # Collect expected outputs
    spawn_multiprocessing(
        parallel_assert_atom_attention_decoder_wb,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        serial_module_version,
        dtype,
        multiplicity,
        atom_s,
        atom_z,
        token_s,
        W,
        H,
        atom_decoder_depth,
        atom_decoder_heads,
        layer_state_dict,
        {k: v.detach().cpu() for k, v in feats.items()},
        a.detach().cpu(),
        q.detach().cpu(),
        c.detach().cpu(),
        p.detach().cpu(),
        d_r_update.detach().cpu(),
        r_update_expected.detach().cpu(),
        a_serial.grad.detach().cpu(),
        q_serial.grad.detach().cpu(),
        c_serial.grad.detach().cpu(),
        p_serial.grad.detach().cpu(),
        {
            name: param.grad.detach().cpu()
            for name, param in reference_module.named_parameters()
            if param.grad is not None
        },
    )


# ======================================================================
# Test 2: AtomAttentionDecoder under autocast bf16 (dtype-only comparison)
# ======================================================================


def parallel_assert_atom_attention_decoder_wb_autocast_bf16(
    rank: int,
    grid_group_sizes,
    device_type: str,
    backend: str,
    env_per_rank,
    atom_s: int,
    token_s: int,
    W: int,
    H: int,
    atom_decoder_depth: int,
    atom_decoder_heads: int,
    layer_state_dict,
    feats_global_host: dict[str, torch.Tensor],
    a_global_host: torch.Tensor,
    q_global_host: torch.Tensor,
    c_global_host: torch.Tensor,
    p_global_host: torch.Tensor,
    serial_output_dtype: torch.dtype,
    serial_grad_dtypes: dict[str, torch.dtype],
    serial_param_grad_dtypes: dict[str, torch.dtype],
):
    """Parallel worker for bf16 autocast dtype test on AtomAttentionDecoder."""
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
    multiplicity = 1

    module_serial = SerialAtomAttentionDecoderBoltz2(
        atom_s=atom_s,
        token_s=token_s,
        attn_window_queries=W,
        attn_window_keys=H,
        atom_decoder_depth=atom_decoder_depth,
        atom_decoder_heads=atom_decoder_heads,
    ).to(device=manager.device, dtype=dtype)
    module_serial.load_state_dict(layer_state_dict)

    module = DistributedAtomAttentionDecoder(
        layer=module_serial,
        device_mesh=manager.device_mesh_subgroups,
    ).train()

    # Distribute atom features
    inputs_atom = {
        k: v.to(dtype=dtype if v.dtype.is_floating_point else v.dtype)
        for k, v in feats_global_host.items()
        if k in _placements_cp_atom_features
    }
    feats_dt = distribute_atom_features(
        inputs_atom,
        _placements_cp_atom_features,
        _placements_atom_features,
        manager.device_mesh_subgroups,
        manager.group["cp"],
    )
    feats_dt_packed = pack_atom_features(feats_dt, set(feats_dt.keys()), W)
    N_atoms_packed = feats_dt_packed["atom_pad_mask"].shape[1]
    K_packed = N_atoms_packed // W

    # Distribute inputs
    a_dt = distribute_tensor(
        a_global_host.to(device=manager.device, dtype=dtype),
        manager.device_mesh_subgroups,
        _placements_single,
    ).requires_grad_(True)
    q_padded = pad_or_shrink_to_length(
        q_global_host.to(device=manager.device, dtype=dtype), axis=1, target_length=N_atoms_packed
    )
    q_dt = distribute_tensor(q_padded, manager.device_mesh_subgroups, _placements_single).requires_grad_(True)
    c_padded = pad_or_shrink_to_length(
        c_global_host.to(device=manager.device, dtype=dtype), axis=1, target_length=N_atoms_packed
    )
    c_dt = distribute_tensor(c_padded, manager.device_mesh_subgroups, _placements_single).requires_grad_(True)
    p_padded = pad_or_shrink_to_length(
        p_global_host.to(device=manager.device, dtype=dtype), axis=1, target_length=K_packed
    )
    p_dt = distribute_tensor(p_padded, manager.device_mesh_subgroups, _placements_single).requires_grad_(True)

    # Forward under autocast
    with torch.autocast("cuda", dtype=torch.bfloat16):
        r_update_dt = module(a=a_dt, q=q_dt, c=c_dt, p=p_dt, feats=feats_dt_packed, multiplicity=multiplicity)

    torch.autograd.backward([r_update_dt], [torch.ones_like(r_update_dt)])

    # Assert output dtype
    assert (
        r_update_dt.dtype == serial_output_dtype
    ), f"r_update dtype mismatch: DTensor {r_update_dt.dtype} vs serial {serial_output_dtype}"

    # Assert input grad dtypes
    for name, dt_tensor in [("a", a_dt), ("q", q_dt), ("c", c_dt)]:
        assert dt_tensor.grad is not None, f"{name} grad is None"
        assert (
            dt_tensor.grad.dtype == serial_grad_dtypes[name]
        ), f"{name} grad dtype mismatch: DTensor {dt_tensor.grad.dtype} vs serial {serial_grad_dtypes[name]}"

    # Assert param grad dtypes
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
    indirect=["setup_env"],
    ids=lambda x: f"dp:{x[0][0]}, cp:{x[0][1]}, device_type:{x[2]}",
)
def test_atom_attention_decoder_wb_autocast_bf16(setup_env):
    """Test DTensor AtomAttentionDecoder output dtypes under autocast bf16."""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() < {world_size}")

    seed = 42
    seed_by_rank(0, seed=seed)

    B = 1
    W = 32
    H = 128
    val_init_min_max = (-0.1, 0.1)
    dtype = torch.float32
    multiplicity = 1

    n_atoms_per_token_min = 8
    n_atoms_per_token_max = 20
    N_tokens = 30
    N_atoms_raw = N_tokens * n_atoms_per_token_max
    N_atoms = ((N_atoms_raw + W - 1) // W) * W
    N_msa = 1

    atom_s = 8
    token_s = 2
    atom_decoder_depth = 2
    atom_decoder_heads = 2
    p_last_dim = atom_decoder_heads * atom_decoder_depth

    feats = random_features(
        size_batch=B,
        n_tokens=N_tokens,
        n_atoms=N_atoms,
        n_msa=N_msa,
        atom_counts_per_token_range=(n_atoms_per_token_min, n_atoms_per_token_max),
        device=torch.device(device_type),
        float_value_range=val_init_min_max,
        selected_keys=list(_selected_atom_keys),
    )
    feats = {k: v.to(dtype=dtype) if v.dtype.is_floating_point else v for k, v in feats.items()}
    N_atoms_actual = feats["atom_pad_mask"].shape[1]
    K = N_atoms_actual // W

    a = torch.empty((B * multiplicity, N_tokens, token_s * 2), device=device_type, dtype=dtype, requires_grad=True)
    q = torch.empty((B * multiplicity, N_atoms_actual, atom_s), device=device_type, dtype=dtype, requires_grad=True)
    c = torch.empty((B * multiplicity, N_atoms_actual, atom_s), device=device_type, dtype=dtype, requires_grad=True)
    p = torch.empty((B, K, W, H, p_last_dim), device=device_type, dtype=dtype, requires_grad=True)
    init_tensors_uniform([a, q, c, p], low=val_init_min_max[0], high=val_init_min_max[1])

    reference_module = SerialAtomAttentionDecoderBoltz2(
        atom_s=atom_s,
        token_s=token_s,
        attn_window_queries=W,
        attn_window_keys=H,
        atom_decoder_depth=atom_decoder_depth,
        atom_decoder_heads=atom_decoder_heads,
    ).to(device=device_type, dtype=dtype)
    reference_module.train()
    init_module_params_uniform(reference_module, low=val_init_min_max[0], high=val_init_min_max[1])
    reference_module.apply(SetModuleInfValues())
    layer_state_dict = reference_module.state_dict()

    # Serial forward under autocast
    a_serial = a.detach().clone().requires_grad_(True)
    q_serial = q.detach().clone().requires_grad_(True)
    c_serial = c.detach().clone().requires_grad_(True)
    p_serial = p.detach().clone().requires_grad_(True)

    index_matrix = get_indexing_matrix_v2(K, W, H, device_type).to(dtype=dtype)
    to_keys = partial(single_to_keys_v2, indexing_matrix=index_matrix, W=W, H=H)

    with torch.autocast("cuda", dtype=torch.bfloat16):
        r_update_serial = reference_module(
            a=a_serial,
            q=q_serial,
            c=c_serial,
            atom_dec_bias=p_serial,
            feats={k: v.clone() for k, v in feats.items()},
            to_keys=to_keys,
            multiplicity=multiplicity,
        )

    torch.autograd.backward([r_update_serial], [torch.ones_like(r_update_serial)])

    serial_output_dtype = r_update_serial.dtype
    serial_grad_dtypes = {"a": a_serial.grad.dtype, "q": q_serial.grad.dtype, "c": c_serial.grad.dtype}
    serial_param_grad_dtypes = {
        name: param.grad.dtype for name, param in reference_module.named_parameters() if param.grad is not None
    }

    spawn_multiprocessing(
        parallel_assert_atom_attention_decoder_wb_autocast_bf16,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        atom_s,
        token_s,
        W,
        H,
        atom_decoder_depth,
        atom_decoder_heads,
        {k: v.detach().cpu() for k, v in layer_state_dict.items()},
        {k: v.detach().cpu() for k, v in feats.items()},
        a.detach().cpu(),
        q.detach().cpu(),
        c.detach().cpu(),
        p.detach().cpu(),
        serial_output_dtype,
        serial_grad_dtypes,
        serial_param_grad_dtypes,
    )
