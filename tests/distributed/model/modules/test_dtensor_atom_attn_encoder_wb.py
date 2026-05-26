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

"""Tests for DTensor AtomAttentionEncoder with window batching.

Tests the DTensor AtomAttentionEncoder against V1 and V2 serial references,
verifying forward and backward numerical equivalence.

Uses float64 to enable exact (default tolerance) comparison between
serial and distributed computation paths.
"""

from functools import partial

import pytest
import torch
from torch.distributed.tensor import distribute_tensor

from boltz.distributed.data.feature.featurizer import pack_atom_features
from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.layers.flatten_and_unflatten import shardwise_unflatten_sharded
from boltz.distributed.model.layers.utils import convert_single_repr_to_window_batched_key
from boltz.distributed.model.modules.encoders import (
    AtomAttentionEncoder as DistributedAtomAttentionEncoder,
)
from boltz.model.modules.encoders import AtomAttentionEncoder as SerialAtomAttentionEncoderBoltz1
from boltz.model.modules.encoders import get_indexing_matrix as get_indexing_matrix_v1
from boltz.model.modules.encoders import single_to_keys as single_to_keys_v1
from boltz.model.modules.encodersv2 import AtomAttentionEncoder as SerialAtomAttentionEncoderBoltz2
from boltz.model.modules.encodersv2 import AtomEncoder as SerialAtomEncoderV2
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

# Subset of keys needed for AtomAttentionEncoder test
_selected_atom_keys = {
    "atom_pad_mask",
    "ref_pos",
    "ref_space_uid",
    "ref_charge",
    "ref_element",
    "ref_atom_name_chars",
    "atom_to_token",
    "atom_counts_per_token",
    "token_pad_mask",
}

_placements = get_feature_placements(
    token_keys=set(),
    msa_keys=set(),
    atom_keys=_selected_atom_keys,
    model_io_keys=set(),
    model_io_fp32_keys=set(),
)
_placements_single = _placements["single"]
_placements_pair = _placements["pair"]
_placements_cp_atom_features = _placements["cp_atom_features"]
_placements_atom_features = _placements["atom_features"]


def parallel_assert_atom_attention_encoder_wb(
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
    token_z: int,
    atom_feature_dim: int,
    W: int,
    H: int,
    atom_encoder_depth: int,
    atom_encoder_heads: int,
    structure_prediction: bool,
    layer_state_dict,
    # V2 only: AtomEncoder state dict for generating q/c/p inside worker
    atom_encoder_state_dict,
    feats_global_host: dict[str, torch.Tensor],
    s_trunk_global_host: torch.Tensor | None,
    z_global_host: torch.Tensor | None,
    r_global_host: torch.Tensor | None,
    # Upstream gradients
    d_a_global_host: torch.Tensor,
    d_q_out_global_host: torch.Tensor,
    d_c_out_global_host: torch.Tensor,
    d_p_out_global_host: torch.Tensor,
    # Expected outputs
    a_expected_global_host: torch.Tensor,
    q_out_expected_global_host: torch.Tensor,
    c_out_expected_global_host: torch.Tensor,
    p_out_expected_global_host: torch.Tensor,
    # Expected input grads
    d_s_trunk_expected_global_host: torch.Tensor | None,
    d_z_expected_global_host: torch.Tensor | None,
    d_r_expected_global_host: torch.Tensor | None,
    expected_param_grads_global_host_dict: dict[str, torch.Tensor],
):
    """Parallel worker for DTensor AtomAttentionEncoder window batching test."""
    monkeypatch = pytest.MonkeyPatch()
    if env_per_rank is not None:
        for var_name, value in env_per_rank.items():
            if value == "<INPUT_RANK>":
                monkeypatch.setenv(var_name, f"{rank}")
                continue
            monkeypatch.setenv(var_name, value)

    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    manager = DistributedManager()

    # Recreate serial module -- move to device/dtype BEFORE load_state_dict
    if serial_module_version == "boltz1":
        module_serial = SerialAtomAttentionEncoderBoltz1(
            atom_s=atom_s,
            atom_z=atom_z,
            token_s=token_s,
            token_z=token_z,
            atoms_per_window_queries=W,
            atoms_per_window_keys=H,
            atom_feature_dim=atom_feature_dim,
            atom_encoder_depth=atom_encoder_depth,
            atom_encoder_heads=atom_encoder_heads,
            structure_prediction=structure_prediction,
        )
    else:
        module_serial = SerialAtomAttentionEncoderBoltz2(
            atom_s=atom_s,
            token_s=token_s,
            atoms_per_window_queries=W,
            atoms_per_window_keys=H,
            atom_encoder_depth=atom_encoder_depth,
            atom_encoder_heads=atom_encoder_heads,
            structure_prediction=structure_prediction,
        )
    module_serial = module_serial.to(device=manager.device, dtype=dtype)
    module_serial.load_state_dict(layer_state_dict)
    module_serial = module_serial.train()
    module_serial.apply(SetModuleInfValues())

    # Create distributed module
    module = DistributedAtomAttentionEncoder(
        layer=module_serial,
        device_mesh=manager.device_mesh_subgroups,
    ).train()

    # Global masks
    token_pad_mask_global = feats_global_host.pop("token_pad_mask").to(device=manager.device, dtype=torch.bool)
    token_pad_mask_expanded_global = token_pad_mask_global.unsqueeze(-1)

    atom_pad_mask_global = feats_global_host["atom_pad_mask"].to(device=manager.device, dtype=torch.bool)
    atom_pad_mask_expanded_global = atom_pad_mask_global.unsqueeze(-1)
    atom_pad_mask_expanded_global_mul = atom_pad_mask_expanded_global.repeat_interleave(multiplicity, dim=0)

    # ========================================================================
    # Distribute atom features
    # ========================================================================
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

    # Distribute token-level tensors
    s_trunk_dt = None
    z_dt = None
    r_dt = None
    if structure_prediction:
        s_trunk_dt = distribute_tensor(
            s_trunk_global_host.to(device=manager.device, dtype=dtype),
            manager.device_mesh_subgroups,
            _placements_single,
        ).requires_grad_(True)

        z_dt = distribute_tensor(
            z_global_host.to(device=manager.device, dtype=dtype),
            manager.device_mesh_subgroups,
            _placements_pair,
        ).requires_grad_(True)

        r_padded = pad_or_shrink_to_length(
            r_global_host.to(device=manager.device, dtype=dtype), axis=1, target_length=N_atoms_packed
        )
        r_dt = distribute_tensor(r_padded, manager.device_mesh_subgroups, _placements_single).requires_grad_(True)

    # For V2: compute q, c from DTensor AtomEncoder; bias is passed in separately
    q_dt = None
    c_dt = None
    atom_enc_bias_dt = None
    if serial_module_version == "boltz2":
        atom_encoder_serial = SerialAtomEncoderV2(
            atom_s=atom_s,
            atom_z=atom_z,
            token_s=token_s,
            token_z=token_z,
            atoms_per_window_queries=W,
            atoms_per_window_keys=H,
            atom_feature_dim=atom_feature_dim,
            structure_prediction=structure_prediction,
        )
        atom_encoder_serial = atom_encoder_serial.to(device=manager.device, dtype=dtype)
        atom_encoder_serial.load_state_dict(atom_encoder_state_dict)
        atom_encoder_serial.train()
        atom_encoder_serial.apply(SetModuleInfValues())

        from boltz.distributed.model.modules.encoders import AtomEncoder as DistributedAtomEncoder

        atom_encoder_dt = DistributedAtomEncoder(
            layer=atom_encoder_serial,
            device_mesh=manager.device_mesh_subgroups,
        ).train()

        q_dt, c_dt, _ = atom_encoder_dt(feats=feats_dt_packed, s_trunk=s_trunk_dt, z=z_dt)

        # Distribute the pre-generated atom_enc_bias (p_out_expected_global_host for V2 is the bias)
        # Pad to packed K windows
        p_padded = pad_or_shrink_to_length(
            p_out_expected_global_host.to(device=manager.device, dtype=dtype), axis=1, target_length=K_packed
        )
        atom_enc_bias_dt = distribute_tensor(
            p_padded, manager.device_mesh_subgroups, _placements_single
        ).requires_grad_(True)

    # ========================================================================
    # Forward pass
    # ========================================================================
    a_dt, q_out_dt, c_out_dt, p_out_dt = module(
        feats=feats_dt_packed,
        q=q_dt,
        c=c_dt,
        atom_enc_bias=atom_enc_bias_dt,
        s_trunk=s_trunk_dt if serial_module_version == "boltz1" else None,
        z=z_dt if serial_module_version == "boltz1" else None,
        r=r_dt if structure_prediction else None,
        multiplicity=multiplicity,
    )

    # ========================================================================
    # Forward comparison
    # ========================================================================
    # a is token feature - compare with token mask
    a_expected_device = a_expected_global_host.to(device=manager.device, dtype=dtype)
    token_pad_mask_expanded_global_mul = token_pad_mask_expanded_global.repeat_interleave(multiplicity, dim=0)
    torch.testing.assert_close(
        a_dt.full_tensor() * token_pad_mask_expanded_global_mul,
        a_expected_device * token_pad_mask_expanded_global_mul,
    )

    # q_out and c_out are atom features with multiplicity
    mask_dt_full = feats_dt_packed["atom_pad_mask"].full_tensor()
    mask_dt_full_mul = mask_dt_full.repeat_interleave(multiplicity, dim=0)
    mask_dt_full_mul_expanded = mask_dt_full_mul.unsqueeze(-1)

    q_out_expected_device = q_out_expected_global_host.to(device=manager.device, dtype=dtype)
    assert_tensors_close_with_pad(
        q_out_dt.full_tensor() * mask_dt_full_mul_expanded,
        q_out_expected_device * atom_pad_mask_expanded_global_mul,
        axis=1,
        pad_val=0,
    )

    c_out_expected_device = c_out_expected_global_host.to(device=manager.device, dtype=dtype)
    assert_tensors_close_with_pad(
        c_out_dt.full_tensor() * mask_dt_full_mul_expanded,
        c_out_expected_device * atom_pad_mask_expanded_global_mul,
        axis=1,
        pad_val=0,
    )

    # Compare only the valid 'key' region of the pair repr.
    # Due to pack_atom_features and the resulting difference in atom length,
    # the two pair repr (DTensor vs serial) can have different number of (W, H) windows
    # and the extra windows in either case should be invalid by definition of pack_atom_features'
    # guaranteeing not removing valid atoms. However, for comparing the two pair repr for numerical
    # consistency, we need to mask both (W, H) axes because otherwise the last window
    # can contain non-zero values for the invalid query atoms, failing assert_tensors_close_with_pad.
    # Example of last two windows' mask (from Boltz-1x CP test with W=32, H=128):
    #   mask_dt_key_full_expanded[0, -2:, 0, :, 0]  -- key mask shows partial validity:
    #     window -2: [1,1,...,1, 0,0,...,0]  (51 valid keys, 77 padding)
    #     window -1: [1,1,...,1, 0,0,...,0]  (19 valid keys, 109 padding)
    #   mask_dt_query_full_expanded[0, -2:, :, 0, 0]  -- query mask shows partial validity:
    #     window -2: [1,1,1, 0,...,0]  (3 valid queries, 29 padding)
    #     window -1: [0,0,...,0]  (all padding -- entirely invalid window)
    # Without masking both axes, the all-padding window -1 would have non-zero pair values
    # from the forward pass (computed on garbage padding data) that don't exist in the serial.
    N_atoms_serial = feats_global_host["atom_pad_mask"].shape[1]
    K_serial = N_atoms_serial // W

    mask_dt_query = shardwise_unflatten_sharded(feats_dt_packed["atom_pad_mask"], axis=1, sizes=(K_packed, W))
    mask_dt_query_full = mask_dt_query.full_tensor()
    mask_dt_query_full_expanded = mask_dt_query_full[:, :, :, None, None]
    mask_dt_key = convert_single_repr_to_window_batched_key(feats_dt_packed["atom_pad_mask"], W, H)
    mask_dt_key_full = mask_dt_key.full_tensor()
    mask_dt_key_full_expanded = mask_dt_key_full[:, :, None, :, None]
    mask_dt_pair_full_expanded = mask_dt_query_full_expanded * mask_dt_key_full_expanded

    if serial_module_version == "boltz1":
        index_matrix = get_indexing_matrix_v1(K_serial, W, H, manager.device).to(dtype=dtype)
        to_keys_fn = partial(single_to_keys_v1, indexing_matrix=index_matrix, W=W, H=H)
    else:
        compute_dtype = torch.promote_types(dtype, torch.float32)
        index_matrix = get_indexing_matrix_v2(K_serial, W, H, manager.device).to(dtype=compute_dtype)
        to_keys_fn = partial(single_to_keys_v2, indexing_matrix=index_matrix, W=W, H=H)

    mask_key_expected = to_keys_fn(
        feats_global_host["atom_pad_mask"].to(device=manager.device, dtype=dtype).unsqueeze(-1)
    )
    mask_key_expected_expanded = mask_key_expected[:, :, None, :, :]
    mask_query_expected_expanded = atom_pad_mask_expanded_global.unflatten(
        1, (atom_pad_mask_expanded_global.shape[1] // W, W)
    )[:, :, :, None, :]
    mask_pair_expected_expanded = mask_query_expected_expanded * mask_key_expected_expanded

    p_out_expected_device = p_out_expected_global_host.to(device=manager.device, dtype=dtype)
    assert_tensors_close_with_pad(
        p_out_dt.full_tensor() * mask_dt_pair_full_expanded,
        p_out_expected_device * mask_pair_expected_expanded,
        axis=1,
        pad_val=0,
    )

    # ========================================================================
    # Backward pass
    # ========================================================================
    d_a_dt = distribute_tensor(
        d_a_global_host.to(device=manager.device, dtype=dtype),
        manager.device_mesh_subgroups,
        a_dt.placements,
    )

    d_q_out_padded = pad_or_shrink_to_length(
        d_q_out_global_host.to(device=manager.device, dtype=dtype), axis=1, target_length=q_out_dt.shape[1]
    )
    d_q_out_dt = distribute_tensor(d_q_out_padded, manager.device_mesh_subgroups, q_out_dt.placements)

    d_c_out_padded = pad_or_shrink_to_length(
        d_c_out_global_host.to(device=manager.device, dtype=dtype), axis=1, target_length=c_out_dt.shape[1]
    )
    d_c_out_dt = distribute_tensor(d_c_out_padded, manager.device_mesh_subgroups, c_out_dt.placements)

    d_p_out_padded = pad_or_shrink_to_length(
        d_p_out_global_host.to(device=manager.device, dtype=dtype), axis=1, target_length=p_out_dt.shape[1]
    )
    d_p_out_dt = distribute_tensor(d_p_out_padded, manager.device_mesh_subgroups, p_out_dt.placements)

    torch.autograd.backward(
        [a_dt, q_out_dt, c_out_dt, p_out_dt],
        [d_a_dt, d_q_out_dt, d_c_out_dt, d_p_out_dt],
    )

    # Check input gradients (only for V1 where s_trunk/z flow through AtomAttentionEncoder)
    if structure_prediction and s_trunk_dt is not None and d_s_trunk_expected_global_host is not None:
        d_s_trunk_expected_device = d_s_trunk_expected_global_host.to(device=manager.device, dtype=dtype)
        torch.testing.assert_close(s_trunk_dt.grad.full_tensor(), d_s_trunk_expected_device)

    if structure_prediction and z_dt is not None and d_z_expected_global_host is not None:
        d_z_expected_device = d_z_expected_global_host.to(device=manager.device, dtype=dtype)
        torch.testing.assert_close(z_dt.grad.full_tensor(), d_z_expected_device)

    if structure_prediction and r_dt is not None:
        d_r_expected_device = d_r_expected_global_host.to(device=manager.device, dtype=dtype)
        r_grad_full = r_dt.grad.full_tensor()
        assert_tensors_close_with_pad(
            r_grad_full * mask_dt_full_mul_expanded[:, :, :3],
            d_r_expected_device * atom_pad_mask_expanded_global_mul[:, :, :3],
            axis=1,
            pad_val=0,
        )

    # Parameter grads
    for name, grad_expected_global in expected_param_grads_global_host_dict.items():
        grad_param = get_param_by_key(module, name).grad
        assert grad_param is not None, f"Missing grad for param {name}"

        if hasattr(grad_param, "full_tensor"):
            grad_to_check = grad_param.full_tensor()
        else:
            grad_to_check = grad_param

        torch.testing.assert_close(grad_to_check.cpu(), grad_expected_global.to(dtype=dtype))
        assert_all_identical(grad_to_check, manager.group["cp"])

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.slow
@pytest.mark.parametrize(
    "setup_env, dtype, multiplicity",
    (
        params_test := [
            (((1, (2, 2)), True, "cuda", "ENV"), torch.float64, 1),
            (((2, (2, 2)), True, "cuda", "ENV"), torch.float64, 4),
        ]
    ),
    indirect=["setup_env"],
    ids=[
        f"dp:{x[0][0][0]}, cp:{x[0][0][1]}, device_type:{x[0][2]}, dtype:{x[1]}, multiplicity:{x[2]}"
        for x in params_test
    ],
)
# TODO: Add "boltz1" to serial_module_version to test the V1 internalized_AtomEncoder path.
# Requirements for boltz1 test:
#   - V1 serial AtomAttentionEncoder (monolithic: embed + pair + r_to_q + transformer + scatter)
#   - V1 feature set includes atom_pad_mask in the atom_feats concat (V2 does not)
#   - V1 r_to_q_trans takes 10-dim input: concat([r, zeros(B*M, N, 7)])
#   - V1 serial forward signature: forward(feats, s_trunk, z, r, multiplicity, model_cache)
#   - V1 returns 5 values: (a, q, c, p, to_keys) vs V2's 4: (a, q, c, to_keys)
#   - The DTensor path exercises _atom_encoder() shared function through the V1 code path
@pytest.mark.parametrize("serial_module_version", ["boltz2"])
def test_atom_attention_encoder_window_batching(setup_env, dtype, multiplicity, serial_module_version):
    """Test DTensor AtomAttentionEncoder with window batching."""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() < {world_size}")

    structure_prediction = multiplicity > 1
    seed = 42
    seed_by_rank(0, seed=seed)

    size_cp = grid_group_sizes["cp"][0]
    B = 1 * grid_group_sizes["dp"]

    W = 32
    H = 128
    # Small init range needed because the serial AtomAttentionEncoder uses 1e-6 epsilon in
    # atom_to_token mean normalization (/ (count + 1e-6)) while DTensor uses exact scatter mean.
    # The ~1e-6 relative error scales with value magnitude; (-0.02, 0.02) keeps it within
    # float64 default tolerance (atol=1e-7) while maintaining non-trivial gradient magnitudes
    # (transition layers ~1e-7, attention layers ~1e-4, scatter ~1e-2).
    val_init_min_max = (-0.03, 0.03)

    n_atoms_per_token_min = 8
    n_atoms_per_token_max = 20
    N_tokens = 50 * size_cp
    N_atoms_raw = N_tokens * n_atoms_per_token_max
    N_atoms = ((N_atoms_raw + W - 1) // W) * W
    N_msa = 1

    atom_s = 8
    atom_z = 8
    token_s = 2
    token_z = 2
    atom_encoder_depth = 2
    atom_encoder_heads = 2

    from boltz.data import const as boltz_const

    atom_feature_dim = 3 + 1 + boltz_const.num_elements + 4 * 64

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
    K = N_atoms_actual // W

    # Build serial modules
    # For Boltz-2, p last dim = num_heads * depth (pre-computed bias split for DiffusionTransformer)
    # For Boltz-1, p last dim = atom_z (pair representation)
    p_last_dim = atom_z if serial_module_version == "boltz1" else atom_encoder_heads * atom_encoder_depth

    atom_encoder_state_dict = None
    if serial_module_version == "boltz2":
        # V2: need AtomEncoder to generate q, c (but NOT p -- bias is separate)
        atom_encoder_module = SerialAtomEncoderV2(
            atom_s=atom_s,
            atom_z=atom_z,
            token_s=token_s,
            token_z=token_z,
            atoms_per_window_queries=W,
            atoms_per_window_keys=H,
            atom_feature_dim=atom_feature_dim,
            structure_prediction=structure_prediction,
        ).to(device=device_type, dtype=dtype)
        atom_encoder_module.train()
        init_module_params_uniform(atom_encoder_module, low=val_init_min_max[0], high=val_init_min_max[1])
        atom_encoder_module.apply(SetModuleInfValues())
        atom_encoder_state_dict = {k: v.detach().cpu() for k, v in atom_encoder_module.state_dict().items()}

        reference_module = SerialAtomAttentionEncoderBoltz2(
            atom_s=atom_s,
            token_s=token_s,
            atoms_per_window_queries=W,
            atoms_per_window_keys=H,
            atom_encoder_depth=atom_encoder_depth,
            atom_encoder_heads=atom_encoder_heads,
            structure_prediction=structure_prediction,
        ).to(device=device_type, dtype=dtype)
    else:
        reference_module = SerialAtomAttentionEncoderBoltz1(
            atom_s=atom_s,
            atom_z=atom_z,
            token_s=token_s,
            token_z=token_z,
            atoms_per_window_queries=W,
            atoms_per_window_keys=H,
            atom_feature_dim=atom_feature_dim,
            atom_encoder_depth=atom_encoder_depth,
            atom_encoder_heads=atom_encoder_heads,
            structure_prediction=structure_prediction,
        ).to(device=device_type, dtype=dtype)

    reference_module.train()
    init_module_params_uniform(reference_module, low=val_init_min_max[0], high=val_init_min_max[1])
    reference_module.apply(SetModuleInfValues())
    layer_state_dict = {k: v.detach().cpu() for k, v in reference_module.state_dict().items()}

    # Prepare inputs
    s_trunk = None
    z = None
    r = None
    if structure_prediction:
        s_trunk = torch.empty((B, N_tokens, token_s), device=device_type, dtype=dtype, requires_grad=True)
        z = torch.empty((B, N_tokens, N_tokens, token_z), device=device_type, dtype=dtype, requires_grad=True)
        r = torch.empty((B * multiplicity, N_atoms_actual, 3), device=device_type, dtype=dtype, requires_grad=True)
        init_tensors_uniform([s_trunk, z, r], low=val_init_min_max[0], high=val_init_min_max[1])

    # Serial forward
    feats_serial = {k: v.detach().clone() for k, v in feats.items()}
    s_trunk_serial = s_trunk.detach().clone().requires_grad_(True) if s_trunk is not None else None
    z_serial = z.detach().clone().requires_grad_(True) if z is not None else None
    r_serial = r.detach().clone().requires_grad_(True) if r is not None else None

    if serial_module_version == "boltz2":
        # Run AtomEncoder first to get q, c (we use a separate random bias)
        q_enc, c_enc, _, to_keys_enc = atom_encoder_module(
            feats=feats_serial,
            s_trunk=s_trunk_serial,
            z=z_serial,
        )
        # Create random atom_enc_bias with correct shape (B, K, W, H, num_heads*depth)
        atom_enc_bias = torch.empty((B, K, W, H, p_last_dim), device=device_type, dtype=dtype, requires_grad=True)
        init_tensors_uniform([atom_enc_bias], low=val_init_min_max[0], high=val_init_min_max[1])
        atom_enc_bias_serial = atom_enc_bias.detach().clone().requires_grad_(True)

        # Run AtomAttentionEncoder
        a_expected, q_out_expected, c_out_expected, _ = reference_module(
            feats=feats_serial,
            q=q_enc,
            c=c_enc,
            atom_enc_bias=atom_enc_bias_serial,
            to_keys=to_keys_enc,
            r=r_serial,
            multiplicity=multiplicity,
        )
        p_out_expected = atom_enc_bias_serial
    else:
        a_expected, q_out_expected, c_out_expected, p_out_expected, _ = reference_module(
            feats=feats_serial,
            s_trunk=s_trunk_serial,
            z=z_serial,
            r=r_serial,
            multiplicity=multiplicity,
            model_cache=None,
        )

    # Generate upstream gradients
    d_a = torch.empty_like(a_expected)
    d_q_out = torch.empty_like(q_out_expected)
    d_c_out = torch.empty_like(c_out_expected)
    d_p_out = torch.empty_like(p_out_expected)
    init_tensors_uniform([d_a, d_q_out, d_c_out, d_p_out], low=val_init_min_max[0], high=val_init_min_max[1])

    # Mask upstream gradients
    compute_dtype = torch.promote_types(dtype, torch.float32)
    if serial_module_version == "boltz1":
        index_matrix = get_indexing_matrix_v1(K, W, H, device_type).to(dtype=dtype)
        to_keys_mask = partial(single_to_keys_v1, indexing_matrix=index_matrix, W=W, H=H)
    else:
        index_matrix = get_indexing_matrix_v2(K, W, H, device_type).to(dtype=compute_dtype)
        to_keys_mask = partial(single_to_keys_v2, indexing_matrix=index_matrix, W=W, H=H)

    mask_key_expected_full = to_keys_mask(
        feats_serial["atom_pad_mask"].to(dtype=compute_dtype, device=d_p_out.device).unsqueeze(-1)
    )
    d_a = d_a * feats_serial["token_pad_mask"].unsqueeze(-1).repeat_interleave(multiplicity, dim=0)
    d_q_out = d_q_out * feats_serial["atom_pad_mask"].unsqueeze(-1).repeat_interleave(multiplicity, dim=0)
    d_c_out = d_c_out * feats_serial["atom_pad_mask"].unsqueeze(-1).repeat_interleave(multiplicity, dim=0)
    d_p_out = d_p_out * mask_key_expected_full[:, :, None, :, :]

    # Serial backward
    torch.autograd.backward(
        [a_expected, q_out_expected, c_out_expected, p_out_expected],
        [d_a, d_q_out, d_c_out, d_p_out],
    )

    expected_param_grads = {
        name: param.grad.detach().cpu() for name, param in reference_module.named_parameters() if param.grad is not None
    }

    spawn_multiprocessing(
        parallel_assert_atom_attention_encoder_wb,
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
        token_z,
        atom_feature_dim,
        W,
        H,
        atom_encoder_depth,
        atom_encoder_heads,
        structure_prediction,
        layer_state_dict,
        atom_encoder_state_dict,
        {k: v.detach().cpu() for k, v in feats.items()},
        s_trunk.detach().cpu() if s_trunk is not None else None,
        z.detach().cpu() if z is not None else None,
        r.detach().cpu() if r is not None else None,
        d_a.detach().cpu(),
        d_q_out.detach().cpu(),
        d_c_out.detach().cpu(),
        d_p_out.detach().cpu(),
        a_expected.detach().cpu(),
        q_out_expected.detach().cpu(),
        c_out_expected.detach().cpu(),
        p_out_expected.detach().cpu(),
        s_trunk_serial.grad.detach().cpu() if s_trunk_serial is not None and s_trunk_serial.grad is not None else None,
        z_serial.grad.detach().cpu() if z_serial is not None and z_serial.grad is not None else None,
        r_serial.grad.detach().cpu() if r_serial is not None and r_serial.grad is not None else None,
        expected_param_grads,
    )


# ======================================================================
# Test 2: AtomAttentionEncoder under autocast bf16 (dtype-only comparison)
# ======================================================================


def parallel_assert_atom_attention_encoder_wb_autocast_bf16(
    rank: int,
    grid_group_sizes,
    device_type: str,
    backend: str,
    env_per_rank,
    atom_s: int,
    atom_z: int,
    token_s: int,
    token_z: int,
    atom_feature_dim: int,
    W: int,
    H: int,
    atom_encoder_depth: int,
    atom_encoder_heads: int,
    layer_state_dict,
    atom_encoder_state_dict,
    feats_global_host: dict[str, torch.Tensor],
    s_trunk_global_host: torch.Tensor,
    z_global_host: torch.Tensor,
    r_global_host: torch.Tensor,
    serial_output_dtypes: dict[str, torch.dtype],
    serial_grad_dtypes: dict[str, torch.dtype],
    serial_param_grad_dtypes: dict[str, torch.dtype],
):
    """Parallel worker for bf16 autocast dtype test on AtomAttentionEncoder."""
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
    structure_prediction = True

    module_serial = SerialAtomAttentionEncoderBoltz2(
        atom_s=atom_s,
        token_s=token_s,
        atoms_per_window_queries=W,
        atoms_per_window_keys=H,
        atom_encoder_depth=atom_encoder_depth,
        atom_encoder_heads=atom_encoder_heads,
        structure_prediction=structure_prediction,
    )
    module_serial = module_serial.to(device=manager.device, dtype=dtype)
    module_serial.load_state_dict(layer_state_dict)

    module = DistributedAtomAttentionEncoder(
        layer=module_serial,
        device_mesh=manager.device_mesh_subgroups,
    ).train()

    # Distribute atom features
    feats_global_host.pop("token_pad_mask", None)
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

    s_trunk_dt = distribute_tensor(
        s_trunk_global_host.to(device=manager.device, dtype=dtype),
        manager.device_mesh_subgroups,
        _placements_single,
    ).requires_grad_(True)
    z_dt = distribute_tensor(
        z_global_host.to(device=manager.device, dtype=dtype),
        manager.device_mesh_subgroups,
        _placements_pair,
    ).requires_grad_(True)
    r_padded = pad_or_shrink_to_length(
        r_global_host.to(device=manager.device, dtype=dtype),
        axis=1,
        target_length=N_atoms_packed,
    )
    r_dt = distribute_tensor(r_padded, manager.device_mesh_subgroups, _placements_single).requires_grad_(True)

    # Create DTensor AtomEncoder for q/c
    from boltz.distributed.model.modules.encoders import AtomEncoder as DistributedAtomEncoder

    atom_encoder_serial = SerialAtomEncoderV2(
        atom_s=atom_s,
        atom_z=atom_z,
        token_s=token_s,
        token_z=token_z,
        atoms_per_window_queries=W,
        atoms_per_window_keys=H,
        atom_feature_dim=atom_feature_dim,
        structure_prediction=structure_prediction,
    ).to(device=manager.device, dtype=dtype)
    atom_encoder_serial.load_state_dict(atom_encoder_state_dict)
    atom_encoder_serial.eval()

    atom_encoder_dt = DistributedAtomEncoder(
        layer=atom_encoder_serial, device_mesh=manager.device_mesh_subgroups
    ).eval()

    with torch.no_grad(), torch.autocast("cuda", dtype=torch.bfloat16):
        q_dt, c_dt, _ = atom_encoder_dt(feats=feats_dt_packed, s_trunk=s_trunk_dt.detach(), z=z_dt.detach())

    p_last_dim = atom_encoder_heads * atom_encoder_depth
    K_packed = N_atoms_packed // W
    atom_enc_bias_dt = distribute_tensor(
        torch.randn(1, K_packed, W, H, p_last_dim, device=manager.device, dtype=dtype).expand(
            s_trunk_dt.shape[0], -1, -1, -1, -1
        ),
        manager.device_mesh_subgroups,
        _placements_single,
    )

    # Forward under autocast
    with torch.autocast("cuda", dtype=torch.bfloat16):
        a_dt, q_out_dt, c_out_dt, _ = module(
            feats=feats_dt_packed,
            q=q_dt,
            c=c_dt,
            atom_enc_bias=atom_enc_bias_dt,
            s_trunk=s_trunk_dt,
            z=z_dt,
            r=r_dt,
            multiplicity=multiplicity,
        )

    outputs_with_grad = [(n, t) for n, t in [("a", a_dt), ("q_out", q_out_dt), ("c_out", c_out_dt)] if t.requires_grad]
    torch.autograd.backward(
        [t for _, t in outputs_with_grad],
        [torch.ones_like(t) for _, t in outputs_with_grad],
    )

    # Assert output dtypes
    for name, dt_tensor in [("a", a_dt), ("q_out", q_out_dt), ("c_out", c_out_dt)]:
        assert (
            dt_tensor.dtype == serial_output_dtypes[name]
        ), f"{name} dtype mismatch: DTensor {dt_tensor.dtype} vs serial {serial_output_dtypes[name]}"

    # Assert input grad dtypes
    for name, dt_tensor in [("s_trunk", s_trunk_dt), ("z", z_dt), ("r", r_dt)]:
        if name not in serial_grad_dtypes:
            continue
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
def test_atom_attention_encoder_wb_autocast_bf16(setup_env):
    """Test DTensor AtomAttentionEncoder output dtypes under autocast bf16."""
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
    structure_prediction = True

    n_atoms_per_token_min = 8
    n_atoms_per_token_max = 20
    N_tokens = 30
    N_atoms_raw = N_tokens * n_atoms_per_token_max
    N_atoms = ((N_atoms_raw + W - 1) // W) * W
    N_msa = 1

    atom_s = 8
    atom_z = 8
    token_s = 2
    token_z = 2
    atom_encoder_depth = 2
    atom_encoder_heads = 2

    from boltz.data import const as boltz_const

    atom_feature_dim = 3 + 1 + boltz_const.num_elements + 4 * 64

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
    p_last_dim = atom_encoder_heads * atom_encoder_depth

    s_trunk = torch.empty((B, N_tokens, token_s), device=device_type, dtype=dtype, requires_grad=True)
    z = torch.empty((B, N_tokens, N_tokens, token_z), device=device_type, dtype=dtype, requires_grad=True)
    r = torch.empty((B * multiplicity, N_atoms_actual, 3), device=device_type, dtype=dtype, requires_grad=True)
    init_tensors_uniform([s_trunk, z, r], low=val_init_min_max[0], high=val_init_min_max[1])

    atom_encoder_module = SerialAtomEncoderV2(
        atom_s=atom_s,
        atom_z=atom_z,
        token_s=token_s,
        token_z=token_z,
        atoms_per_window_queries=W,
        atoms_per_window_keys=H,
        atom_feature_dim=atom_feature_dim,
        structure_prediction=structure_prediction,
    ).to(device=device_type, dtype=dtype)
    atom_encoder_module.eval()
    init_module_params_uniform(atom_encoder_module, low=val_init_min_max[0], high=val_init_min_max[1])
    atom_encoder_module.apply(SetModuleInfValues())
    atom_encoder_state_dict = {k: v.detach().cpu() for k, v in atom_encoder_module.state_dict().items()}

    reference_module = SerialAtomAttentionEncoderBoltz2(
        atom_s=atom_s,
        token_s=token_s,
        atoms_per_window_queries=W,
        atoms_per_window_keys=H,
        atom_encoder_depth=atom_encoder_depth,
        atom_encoder_heads=atom_encoder_heads,
        structure_prediction=structure_prediction,
    ).to(device=device_type, dtype=dtype)
    reference_module.train()
    init_module_params_uniform(reference_module, low=val_init_min_max[0], high=val_init_min_max[1])
    reference_module.apply(SetModuleInfValues())
    layer_state_dict = {k: v.detach().cpu() for k, v in reference_module.state_dict().items()}

    # Serial forward under autocast
    s_trunk_serial = s_trunk.detach().clone().requires_grad_(True)
    z_serial = z.detach().clone().requires_grad_(True)
    r_serial = r.detach().clone().requires_grad_(True)

    with torch.no_grad(), torch.autocast("cuda", dtype=torch.bfloat16):
        q_enc, c_enc, _, _ = atom_encoder_module(
            feats={k: v.clone() for k, v in feats.items()},
            s_trunk=s_trunk.detach(),
            z=z.detach(),
        )

    atom_enc_bias = torch.randn(B, K, W, H, p_last_dim, device=device_type, dtype=dtype)

    with torch.autocast("cuda", dtype=torch.bfloat16):
        a_serial, q_out_serial, c_out_serial, _ = reference_module(
            feats={k: v.clone() for k, v in feats.items()},
            q=q_enc,
            c=c_enc,
            atom_enc_bias=atom_enc_bias,
            to_keys=partial(
                single_to_keys_v2,
                indexing_matrix=get_indexing_matrix_v2(K, W, H, device_type).to(dtype=torch.float32),
                W=W,
                H=H,
            ),
            r=r_serial,
            multiplicity=multiplicity,
        )

    outputs_with_grad = [
        (n, t) for n, t in [("a", a_serial), ("q_out", q_out_serial), ("c_out", c_out_serial)] if t.requires_grad
    ]
    torch.autograd.backward(
        [t for _, t in outputs_with_grad],
        [torch.ones_like(t) for _, t in outputs_with_grad],
    )

    serial_output_dtypes = {"a": a_serial.dtype, "q_out": q_out_serial.dtype, "c_out": c_out_serial.dtype}
    serial_grad_dtypes = {}
    if s_trunk_serial.grad is not None:
        serial_grad_dtypes["s_trunk"] = s_trunk_serial.grad.dtype
    if z_serial.grad is not None:
        serial_grad_dtypes["z"] = z_serial.grad.dtype
    if r_serial.grad is not None:
        serial_grad_dtypes["r"] = r_serial.grad.dtype
    serial_param_grad_dtypes = {
        name: param.grad.dtype for name, param in reference_module.named_parameters() if param.grad is not None
    }

    spawn_multiprocessing(
        parallel_assert_atom_attention_encoder_wb_autocast_bf16,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        atom_s,
        atom_z,
        token_s,
        token_z,
        atom_feature_dim,
        W,
        H,
        atom_encoder_depth,
        atom_encoder_heads,
        layer_state_dict,
        atom_encoder_state_dict,
        {k: v.detach().cpu() for k, v in feats.items()},
        s_trunk.detach().cpu(),
        z.detach().cpu(),
        r.detach().cpu(),
        serial_output_dtypes,
        serial_grad_dtypes,
        serial_param_grad_dtypes,
    )
