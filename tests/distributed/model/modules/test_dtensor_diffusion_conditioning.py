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

"""Tests for DTensor DiffusionConditioning module (V2 only).

Tests the DTensor DiffusionConditioning against the V2 serial reference,
verifying forward and backward numerical equivalence.

Uses float64 with default tolerance for exact comparison.
"""

from functools import partial

import pytest
import torch
from torch.distributed.tensor import Replicate, Shard, distribute_tensor

from boltz.data import const as boltz_const
from boltz.distributed.data.feature.featurizer import pack_atom_features
from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.layers.flatten_and_unflatten import shardwise_unflatten_sharded
from boltz.distributed.model.layers.utils import convert_single_repr_to_window_batched_key
from boltz.distributed.model.modules.diffusion_conditioning import (
    DiffusionConditioning as DistributedDiffusionConditioning,
)
from boltz.model.modules.diffusion_conditioning import DiffusionConditioning as SerialDiffusionConditioning
from boltz.model.modules.encodersv2 import get_indexing_matrix, single_to_keys
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

# Subset of keys needed for DiffusionConditioning
_selected_atom_keys = {
    "atom_pad_mask",
    "ref_pos",
    "ref_space_uid",
    "ref_charge",
    "ref_element",
    "ref_atom_name_chars",
    "atom_to_token",
    "atom_counts_per_token",
}

_placements = get_feature_placements(
    token_keys=set(),
    msa_keys=set(),
    atom_keys=_selected_atom_keys,
    model_io_keys=set(),
    model_io_fp32_keys=set(),
)
_placements_cp_atom_features = _placements["cp_atom_features"]
_placements_atom_features = _placements["atom_features"]


def parallel_assert_diffusion_conditioning(
    rank: int,
    grid_group_sizes,
    device_type: str,
    backend: str,
    env_per_rank,
    dtype: torch.dtype,
    # Module dimensions
    atom_s: int,
    atom_z: int,
    token_s: int,
    token_z: int,
    atom_feature_dim: int,
    W: int,
    H: int,
    atom_encoder_depth: int,
    atom_encoder_heads: int,
    token_transformer_depth: int,
    token_transformer_heads: int,
    atom_decoder_depth: int,
    atom_decoder_heads: int,
    layer_state_dict,
    # Inputs
    feats_global_host: dict[str, torch.Tensor],
    s_trunk_global_host: torch.Tensor,
    z_trunk_global_host: torch.Tensor,
    rel_pos_enc_global_host: torch.Tensor,
    # Expected outputs
    q_expected_global_host: torch.Tensor,
    c_expected_global_host: torch.Tensor,
    atom_enc_bias_expected_global_host: torch.Tensor,
    atom_dec_bias_expected_global_host: torch.Tensor,
    token_trans_bias_expected_global_host: torch.Tensor,
    # Upstream grads
    d_q_global_host: torch.Tensor,
    d_c_global_host: torch.Tensor,
    d_atom_enc_bias_global_host: torch.Tensor,
    d_atom_dec_bias_global_host: torch.Tensor,
    d_token_trans_bias_global_host: torch.Tensor,
    # Expected input grads
    d_s_trunk_expected_global_host: torch.Tensor,
    d_z_trunk_expected_global_host: torch.Tensor,
    d_rel_pos_enc_expected_global_host: torch.Tensor,
    # Expected param grads
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

    # Re-create serial module (on this rank's device)
    module_serial = SerialDiffusionConditioning(
        token_s=token_s,
        token_z=token_z,
        atom_s=atom_s,
        atom_z=atom_z,
        atoms_per_window_queries=W,
        atoms_per_window_keys=H,
        atom_encoder_depth=atom_encoder_depth,
        atom_encoder_heads=atom_encoder_heads,
        token_transformer_depth=token_transformer_depth,
        token_transformer_heads=token_transformer_heads,
        atom_decoder_depth=atom_decoder_depth,
        atom_decoder_heads=atom_decoder_heads,
        atom_feature_dim=atom_feature_dim,
    )
    module_serial = module_serial.to(device=manager.device, dtype=dtype)
    module_serial.load_state_dict(layer_state_dict)
    module_serial = module_serial.train()

    # Create DTensor module from serial
    module = DistributedDiffusionConditioning(
        layer=module_serial,
        device_mesh=manager.device_mesh_subgroups,
    ).train()

    # ------------------------------------------------------------------
    # Distribute atom features (unpacked — module packs internally)
    # ------------------------------------------------------------------
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

    # Compute N_atoms_packed from an explicit pack call (for comparison sizing only)
    feats_dt_packed_for_sizing = pack_atom_features(feats_dt, set(feats_dt.keys()), W)
    N_atoms_packed = feats_dt_packed_for_sizing["atom_pad_mask"].shape[1]

    # Global masks for comparison
    atom_pad_mask_global = feats_global_host["atom_pad_mask"].to(device=manager.device, dtype=dtype)
    atom_pad_mask_expanded_global = atom_pad_mask_global.unsqueeze(-1)

    # ------------------------------------------------------------------
    # Distribute token-level tensors
    # ------------------------------------------------------------------
    placements_single = (Shard(0), Shard(1), Replicate())
    placements_pair = (Shard(0), Shard(1), Shard(2))

    s_trunk_dt = distribute_tensor(
        s_trunk_global_host.to(device=manager.device, dtype=dtype),
        manager.device_mesh_subgroups,
        placements_single,
    ).requires_grad_(True)
    z_trunk_dt = distribute_tensor(
        z_trunk_global_host.to(device=manager.device, dtype=dtype),
        manager.device_mesh_subgroups,
        placements_pair,
    ).requires_grad_(True)
    rel_pos_enc_dt = distribute_tensor(
        rel_pos_enc_global_host.to(device=manager.device, dtype=dtype),
        manager.device_mesh_subgroups,
        placements_pair,
    ).requires_grad_(True)

    # ------------------------------------------------------------------
    # Forward pass
    # ------------------------------------------------------------------
    q_dt, c_dt, atom_enc_bias_dt, atom_dec_bias_dt, token_trans_bias_dt = module(
        s_trunk=s_trunk_dt,
        z_trunk=z_trunk_dt,
        relative_position_encoding=rel_pos_enc_dt,
        feats=feats_dt,
    )

    # ------------------------------------------------------------------
    # Forward comparison: q and c (atom-level, padded)
    # ------------------------------------------------------------------
    mask_dt_full = feats_dt_packed_for_sizing["atom_pad_mask"].full_tensor()
    mask_dt_full_expanded = mask_dt_full.unsqueeze(-1)

    q_expected_device = q_expected_global_host.to(device=manager.device, dtype=dtype)
    c_expected_device = c_expected_global_host.to(device=manager.device, dtype=dtype)

    assert_tensors_close_with_pad(
        q_dt.full_tensor() * mask_dt_full_expanded,
        q_expected_device * atom_pad_mask_expanded_global,
        axis=1,
        pad_val=0,
    )
    assert_tensors_close_with_pad(
        c_dt.full_tensor() * mask_dt_full_expanded,
        c_expected_device * atom_pad_mask_expanded_global,
        axis=1,
        pad_val=0,
    )

    # ------------------------------------------------------------------
    # Forward comparison: atom_enc_bias, atom_dec_bias (window-batched atom pair, S(0) S(1) R)
    # ------------------------------------------------------------------
    K_packed = N_atoms_packed // W
    N_atoms_serial = feats_global_host["atom_pad_mask"].shape[1]
    K_serial = N_atoms_serial // W

    mask_dt_query = shardwise_unflatten_sharded(
        feats_dt_packed_for_sizing["atom_pad_mask"], axis=1, sizes=(K_packed, W)
    )
    mask_dt_query_full = mask_dt_query.full_tensor()
    mask_dt_query_full_expanded = mask_dt_query_full[:, :, :, None, None]
    mask_dt_key = convert_single_repr_to_window_batched_key(feats_dt_packed_for_sizing["atom_pad_mask"], W, H)
    mask_dt_key_full = mask_dt_key.full_tensor()
    mask_dt_key_full_expanded = mask_dt_key_full[:, :, None, :, None]
    mask_dt_pair_full_expanded = mask_dt_query_full_expanded * mask_dt_key_full_expanded

    compute_dtype = torch.promote_types(dtype, torch.float32)
    index_matrix = get_indexing_matrix(K_serial, W, H, manager.device).to(dtype=compute_dtype)
    to_keys_fn = partial(single_to_keys, indexing_matrix=index_matrix, W=W, H=H)

    mask_key_expected = to_keys_fn(
        feats_global_host["atom_pad_mask"].to(device=manager.device, dtype=compute_dtype).unsqueeze(-1)
    )
    mask_key_expected_expanded = mask_key_expected[:, :, None, :, :]
    mask_query_expected_expanded = atom_pad_mask_expanded_global.unflatten(
        1, (atom_pad_mask_expanded_global.shape[1] // W, W)
    )[:, :, :, None, :]
    mask_pair_expected_expanded = mask_query_expected_expanded * mask_key_expected_expanded

    for name, bias_dt, bias_expected_host in [
        ("atom_enc_bias", atom_enc_bias_dt, atom_enc_bias_expected_global_host),
        ("atom_dec_bias", atom_dec_bias_dt, atom_dec_bias_expected_global_host),
    ]:
        bias_expected_device = bias_expected_host.to(device=manager.device, dtype=dtype)
        assert_tensors_close_with_pad(
            bias_dt.full_tensor() * mask_dt_pair_full_expanded,
            bias_expected_device * mask_pair_expected_expanded,
            axis=1,
            pad_val=0,
        )

    # ------------------------------------------------------------------
    # Forward comparison: token_trans_bias (token pair level)
    # ------------------------------------------------------------------
    token_trans_bias_expected_device = token_trans_bias_expected_global_host.to(device=manager.device, dtype=dtype)
    torch.testing.assert_close(token_trans_bias_dt.full_tensor(), token_trans_bias_expected_device)

    # ------------------------------------------------------------------
    # Backward pass
    # ------------------------------------------------------------------
    d_q_padded = pad_or_shrink_to_length(
        d_q_global_host.to(device=manager.device, dtype=dtype), axis=1, target_length=N_atoms_packed
    )
    d_c_padded = pad_or_shrink_to_length(
        d_c_global_host.to(device=manager.device, dtype=dtype), axis=1, target_length=N_atoms_packed
    )
    d_atom_enc_bias_padded = pad_or_shrink_to_length(
        d_atom_enc_bias_global_host.to(device=manager.device, dtype=dtype), axis=1, target_length=K_packed
    )
    d_atom_dec_bias_padded = pad_or_shrink_to_length(
        d_atom_dec_bias_global_host.to(device=manager.device, dtype=dtype), axis=1, target_length=K_packed
    )

    d_q_dt = distribute_tensor(d_q_padded, manager.device_mesh_subgroups, q_dt.placements)
    d_c_dt = distribute_tensor(d_c_padded, manager.device_mesh_subgroups, c_dt.placements)
    d_atom_enc_bias_dt = distribute_tensor(
        d_atom_enc_bias_padded, manager.device_mesh_subgroups, atom_enc_bias_dt.placements
    )
    d_atom_dec_bias_dt = distribute_tensor(
        d_atom_dec_bias_padded, manager.device_mesh_subgroups, atom_dec_bias_dt.placements
    )
    d_token_trans_bias_dt = distribute_tensor(
        d_token_trans_bias_global_host.to(device=manager.device, dtype=dtype),
        manager.device_mesh_subgroups,
        token_trans_bias_dt.placements,
    )

    torch.autograd.backward(
        [q_dt, c_dt, atom_enc_bias_dt, atom_dec_bias_dt, token_trans_bias_dt],
        [d_q_dt, d_c_dt, d_atom_enc_bias_dt, d_atom_dec_bias_dt, d_token_trans_bias_dt],
    )

    # Check token-level input gradients
    torch.testing.assert_close(
        s_trunk_dt.grad.full_tensor(),
        d_s_trunk_expected_global_host.to(device=manager.device, dtype=dtype),
    )
    torch.testing.assert_close(
        z_trunk_dt.grad.full_tensor(),
        d_z_trunk_expected_global_host.to(device=manager.device, dtype=dtype),
    )
    torch.testing.assert_close(
        rel_pos_enc_dt.grad.full_tensor(),
        d_rel_pos_enc_expected_global_host.to(device=manager.device, dtype=dtype),
    )

    # Parameter grads
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
    "setup_env, dtype",
    (
        params_test := [
            (((1, (2, 2)), True, "cuda", "ENV"), torch.float64),
            (((2, (2, 2)), True, "cuda", "ENV"), torch.float64),
        ]
    ),
    indirect=["setup_env"],
    ids=[f"dp:{x[0][0][0]}, cp:{x[0][0][1]}, device_type:{x[0][2]}, dtype:{x[1]}" for x in params_test],
)
def test_diffusion_conditioning(setup_env, dtype):
    """Test DTensor DiffusionConditioning (V2) vs serial equivalence."""
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
    val_init_min_max = (-0.08, 0.08)

    n_atoms_per_token_min = 8
    n_atoms_per_token_max = 20
    N_tokens = 50 * size_cp
    N_atoms_raw = N_tokens * n_atoms_per_token_max
    N_atoms = ((N_atoms_raw + W - 1) // W) * W
    N_msa = 1

    atom_s = 8
    atom_z = 8
    token_s = 4
    token_z = 4

    atom_encoder_depth = 2
    atom_encoder_heads = 2
    token_transformer_depth = 3
    token_transformer_heads = 2
    atom_decoder_depth = 2
    atom_decoder_heads = 2

    atom_feature_dim = 3 + 1 + boltz_const.num_elements + 4 * 64

    selected_keys = list(_selected_atom_keys)

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
    K = N_atoms_actual // W

    # Token-level inputs
    s_trunk = torch.empty((B, N_tokens, token_s), device=device_type, dtype=dtype, requires_grad=True)
    z_trunk = torch.empty((B, N_tokens, N_tokens, token_z), device=device_type, dtype=dtype, requires_grad=True)
    rel_pos_enc = torch.empty((B, N_tokens, N_tokens, token_z), device=device_type, dtype=dtype, requires_grad=True)
    init_tensors_uniform([s_trunk, z_trunk, rel_pos_enc], low=val_init_min_max[0], high=val_init_min_max[1])

    # Build serial reference module
    reference_module = SerialDiffusionConditioning(
        token_s=token_s,
        token_z=token_z,
        atom_s=atom_s,
        atom_z=atom_z,
        atoms_per_window_queries=W,
        atoms_per_window_keys=H,
        atom_encoder_depth=atom_encoder_depth,
        atom_encoder_heads=atom_encoder_heads,
        token_transformer_depth=token_transformer_depth,
        token_transformer_heads=token_transformer_heads,
        atom_decoder_depth=atom_decoder_depth,
        atom_decoder_heads=atom_decoder_heads,
        atom_feature_dim=atom_feature_dim,
    ).to(device=device_type, dtype=dtype)
    reference_module.train()
    init_module_params_uniform(reference_module, low=val_init_min_max[0], high=val_init_min_max[1])
    reference_module.apply(SetModuleInfValues())
    layer_state_dict = reference_module.state_dict()

    # Serial forward pass
    feats_serial = {k: v.detach().clone() for k, v in feats.items()}
    s_trunk_serial = s_trunk.detach().clone().requires_grad_(True)
    z_trunk_serial = z_trunk.detach().clone().requires_grad_(True)
    rel_pos_enc_serial = rel_pos_enc.detach().clone().requires_grad_(True)

    q_expected, c_expected, _to_keys, atom_enc_bias_expected, atom_dec_bias_expected, token_trans_bias_expected = (
        reference_module(
            s_trunk=s_trunk_serial,
            z_trunk=z_trunk_serial,
            relative_position_encoding=rel_pos_enc_serial,
            feats=feats_serial,
        )
    )

    # Upstream gradients
    d_q = torch.empty_like(q_expected)
    d_c = torch.empty_like(c_expected)
    d_atom_enc_bias = torch.empty_like(atom_enc_bias_expected)
    d_atom_dec_bias = torch.empty_like(atom_dec_bias_expected)
    d_token_trans_bias = torch.empty_like(token_trans_bias_expected)
    init_tensors_uniform(
        [d_q, d_c, d_atom_enc_bias, d_atom_dec_bias, d_token_trans_bias],
        low=val_init_min_max[0],
        high=val_init_min_max[1],
    )

    # Apply masks to upstream gradients
    mask_expanded = feats_serial["atom_pad_mask"].unsqueeze(-1)
    d_q = d_q * mask_expanded
    d_c = d_c * mask_expanded

    compute_dtype = torch.promote_types(dtype, torch.float32)
    index_matrix = get_indexing_matrix(K, W, H, device_type).to(dtype=compute_dtype)
    to_keys_fn_serial = partial(single_to_keys, indexing_matrix=index_matrix, W=W, H=H)
    mask_key_serial = to_keys_fn_serial(
        feats_serial["atom_pad_mask"].to(dtype=compute_dtype, device=d_atom_enc_bias.device).unsqueeze(-1)
    )
    # Pair mask: (B, K, W, H, 1)
    pair_mask = mask_key_serial[:, :, None, :, :] * mask_expanded.unflatten(1, (K, W))[:, :, :, None, :]
    d_atom_enc_bias = d_atom_enc_bias * pair_mask
    d_atom_dec_bias = d_atom_dec_bias * pair_mask

    # Serial backward
    torch.autograd.backward(
        [q_expected, c_expected, atom_enc_bias_expected, atom_dec_bias_expected, token_trans_bias_expected],
        [d_q, d_c, d_atom_enc_bias, d_atom_dec_bias, d_token_trans_bias],
    )

    expected_param_grads = {
        name: param.grad.detach().cpu() for name, param in reference_module.named_parameters() if param.grad is not None
    }

    spawn_multiprocessing(
        parallel_assert_diffusion_conditioning,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        dtype,
        atom_s,
        atom_z,
        token_s,
        token_z,
        atom_feature_dim,
        W,
        H,
        atom_encoder_depth,
        atom_encoder_heads,
        token_transformer_depth,
        token_transformer_heads,
        atom_decoder_depth,
        atom_decoder_heads,
        {k: v.detach().cpu() for k, v in layer_state_dict.items()},
        {k: v.detach().cpu() for k, v in feats.items()},
        s_trunk.detach().cpu(),
        z_trunk.detach().cpu(),
        rel_pos_enc.detach().cpu(),
        q_expected.detach().cpu(),
        c_expected.detach().cpu(),
        atom_enc_bias_expected.detach().cpu(),
        atom_dec_bias_expected.detach().cpu(),
        token_trans_bias_expected.detach().cpu(),
        d_q.detach().cpu(),
        d_c.detach().cpu(),
        d_atom_enc_bias.detach().cpu(),
        d_atom_dec_bias.detach().cpu(),
        d_token_trans_bias.detach().cpu(),
        s_trunk_serial.grad.detach().cpu(),
        z_trunk_serial.grad.detach().cpu(),
        rel_pos_enc_serial.grad.detach().cpu(),
        expected_param_grads,
    )
