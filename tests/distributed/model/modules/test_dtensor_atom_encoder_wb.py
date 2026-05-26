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

"""Tests for DTensor AtomEncoder (V2) with window batching.

Tests the DTensor AtomEncoder against the V2 serial AtomEncoder reference,
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
from boltz.distributed.model.modules.encoders import AtomEncoder as DistributedAtomEncoder
from boltz.model.modules.encodersv2 import AtomEncoder as SerialAtomEncoderV2
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

# Subset of keys needed for AtomEncoder test
_selected_atom_keys = {
    "atom_pad_mask",
    "ref_pos",
    "ref_space_uid",
    "ref_charge",
    "ref_element",
    "ref_atom_name_chars",
    "atom_to_token",  # Needed by serial module and pack_atom_features (creates atom_to_token_ids_global)
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
_placements_pair = _placements["pair"]
_placements_cp_atom_features = _placements["cp_atom_features"]
_placements_atom_features = _placements["atom_features"]


def parallel_assert_atom_encoder_wb(
    rank: int,
    grid_group_sizes,
    device_type: str,
    backend: str,
    env_per_rank,
    dtype: torch.dtype,
    atom_s: int,
    atom_z: int,
    token_s: int,
    token_z: int,
    atom_feature_dim: int,
    W: int,
    H: int,
    structure_prediction: bool,
    layer_state_dict,
    feats_global_host: dict[str, torch.Tensor],
    s_trunk_global_host: torch.Tensor | None,
    z_global_host: torch.Tensor | None,
    # Expected outputs
    q_expected_global_host: torch.Tensor,
    c_expected_global_host: torch.Tensor,
    p_expected_global_host: torch.Tensor,
    # Upstream gradients
    d_q_global_host: torch.Tensor,
    d_c_global_host: torch.Tensor,
    d_p_global_host: torch.Tensor,
    # Expected input grads
    d_s_trunk_expected_global_host: torch.Tensor | None,
    d_z_expected_global_host: torch.Tensor | None,
    expected_param_grads_global_host_dict: dict[str, torch.Tensor],
):
    """Parallel worker function for testing DTensor AtomEncoder."""
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
    module_serial = SerialAtomEncoderV2(
        atom_s=atom_s,
        atom_z=atom_z,
        token_s=token_s,
        token_z=token_z,
        atoms_per_window_queries=W,
        atoms_per_window_keys=H,
        atom_feature_dim=atom_feature_dim,
        structure_prediction=structure_prediction,
    )
    # CRITICAL: Move module to target device/dtype BEFORE loading state dict.
    # Otherwise float64 state dict values get truncated to float32 during copy_()
    # into the default float32 nn.Linear params, then .to(float64) can't recover precision.
    module_serial = module_serial.to(device=manager.device, dtype=dtype)
    module_serial.load_state_dict(layer_state_dict)
    module_serial = module_serial.train()
    module_serial.apply(SetModuleInfValues())

    # Create distributed module
    module = DistributedAtomEncoder(
        layer=module_serial,
        device_mesh=manager.device_mesh_subgroups,
    ).train()

    # Get global masks
    atom_pad_mask_global = feats_global_host["atom_pad_mask"].to(device=manager.device, dtype=torch.bool)
    atom_pad_mask_expanded_global = atom_pad_mask_global.unsqueeze(-1)

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

    # Distribute token-level tensors (s_trunk, z)
    s_trunk_dt = None
    z_dt = None
    if structure_prediction and s_trunk_global_host is not None:
        s_trunk_dt = distribute_tensor(
            s_trunk_global_host.to(device=manager.device, dtype=dtype),
            manager.device_mesh_subgroups,
            _placements_single,
        ).requires_grad_(True)
    if structure_prediction and z_global_host is not None:
        z_dt = distribute_tensor(
            z_global_host.to(device=manager.device, dtype=dtype),
            manager.device_mesh_subgroups,
            _placements_pair,
        ).requires_grad_(True)

    # ========================================================================
    # Forward pass
    # ========================================================================
    q_dt, c_dt, p_dt = module(feats=feats_dt_packed, s_trunk=s_trunk_dt, z=z_dt)

    # ========================================================================
    # Forward comparison
    # ========================================================================
    q_expected_device = q_expected_global_host.to(device=manager.device, dtype=dtype)
    c_expected_device = c_expected_global_host.to(device=manager.device, dtype=dtype)
    p_expected_device = p_expected_global_host.to(device=manager.device, dtype=dtype)

    mask_dt_full = feats_dt_packed["atom_pad_mask"].full_tensor()
    mask_dt_full_expanded = mask_dt_full.unsqueeze(-1)

    # q and c: (B, N_atoms_packed, atom_s)
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
    K_packed = N_atoms_packed // W
    N_atoms_serial = feats_global_host["atom_pad_mask"].shape[1]
    K_serial = N_atoms_serial // W

    mask_dt_query = shardwise_unflatten_sharded(feats_dt_packed["atom_pad_mask"], axis=1, sizes=(K_packed, W))
    mask_dt_query_full = mask_dt_query.full_tensor()
    mask_dt_query_full_expanded = mask_dt_query_full[:, :, :, None, None]
    mask_dt_key = convert_single_repr_to_window_batched_key(feats_dt_packed["atom_pad_mask"], W, H)
    mask_dt_key_full = mask_dt_key.full_tensor()
    mask_dt_key_full_expanded = mask_dt_key_full[:, :, None, :, None]
    mask_dt_pair_full_expanded = mask_dt_query_full_expanded * mask_dt_key_full_expanded

    index_matrix = get_indexing_matrix(K_serial, W, H, manager.device).to(dtype=dtype)
    to_keys_fn = partial(single_to_keys, indexing_matrix=index_matrix, W=W, H=H)

    mask_key_expected = to_keys_fn(
        feats_global_host["atom_pad_mask"].to(device=manager.device, dtype=dtype).unsqueeze(-1)
    )
    mask_key_expected_expanded = mask_key_expected[:, :, None, :, :]
    mask_query_expected_expanded = atom_pad_mask_expanded_global.unflatten(
        1, (atom_pad_mask_expanded_global.shape[1] // W, W)
    )[:, :, :, None, :]
    mask_pair_expected_expanded = mask_query_expected_expanded * mask_key_expected_expanded

    assert_tensors_close_with_pad(
        p_dt.full_tensor() * mask_dt_pair_full_expanded,
        p_expected_device * mask_pair_expected_expanded,
        axis=1,
        pad_val=0,
    )

    # ========================================================================
    # Backward pass
    # ========================================================================
    d_q_padded = pad_or_shrink_to_length(
        d_q_global_host.to(device=manager.device, dtype=dtype), axis=1, target_length=N_atoms_packed
    )
    d_c_padded = pad_or_shrink_to_length(
        d_c_global_host.to(device=manager.device, dtype=dtype), axis=1, target_length=N_atoms_packed
    )
    d_p_padded = pad_or_shrink_to_length(
        d_p_global_host.to(device=manager.device, dtype=dtype), axis=1, target_length=K_packed
    )

    d_q_dt = distribute_tensor(d_q_padded, manager.device_mesh_subgroups, q_dt.placements)
    d_c_dt = distribute_tensor(d_c_padded, manager.device_mesh_subgroups, c_dt.placements)
    d_p_dt = distribute_tensor(d_p_padded, manager.device_mesh_subgroups, p_dt.placements)

    torch.autograd.backward([q_dt, c_dt, p_dt], [d_q_dt, d_c_dt, d_p_dt])

    # Check token-level input gradients
    if structure_prediction and s_trunk_dt is not None:
        d_s_trunk_expected_device = d_s_trunk_expected_global_host.to(device=manager.device, dtype=dtype)
        torch.testing.assert_close(s_trunk_dt.grad.full_tensor(), d_s_trunk_expected_device)

    if structure_prediction and z_dt is not None:
        d_z_expected_device = d_z_expected_global_host.to(device=manager.device, dtype=dtype)
        torch.testing.assert_close(z_dt.grad.full_tensor(), d_z_expected_device)

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
@pytest.mark.parametrize("structure_prediction", [True, False], ids=lambda x: f"sp:{x}")
def test_atom_encoder_wb(setup_env, dtype, structure_prediction):
    """Test DTensor AtomEncoder (V2) with window batching."""
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
    token_s = 4
    token_z = 4

    # Compute atom_feature_dim: ref_pos(3) + ref_charge(1) + ref_element(128) + ref_atom_name_chars(256)
    from boltz.data import const as boltz_const

    atom_feature_dim = 3 + 1 + boltz_const.num_elements + 4 * 64  # 388 with default settings

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

    # Token-level inputs (only for structure_prediction)
    s_trunk = None
    z = None
    if structure_prediction:
        s_trunk = torch.empty((B, N_tokens, token_s), device=device_type, dtype=dtype, requires_grad=True)
        z = torch.empty((B, N_tokens, N_tokens, token_z), device=device_type, dtype=dtype, requires_grad=True)
        init_tensors_uniform([s_trunk, z], low=val_init_min_max[0], high=val_init_min_max[1])

    # Build serial reference module
    reference_module = SerialAtomEncoderV2(
        atom_s=atom_s,
        atom_z=atom_z,
        token_s=token_s,
        token_z=token_z,
        atoms_per_window_queries=W,
        atoms_per_window_keys=H,
        atom_feature_dim=atom_feature_dim,
        structure_prediction=structure_prediction,
    ).to(device=device_type, dtype=dtype)
    reference_module.train()
    init_module_params_uniform(reference_module, low=val_init_min_max[0], high=val_init_min_max[1])
    reference_module.apply(SetModuleInfValues())
    layer_state_dict = reference_module.state_dict()

    # Serial forward
    feats_serial = {k: v.detach().clone() for k, v in feats.items()}
    s_trunk_serial = s_trunk.detach().clone().requires_grad_(True) if s_trunk is not None else None
    z_serial = z.detach().clone().requires_grad_(True) if z is not None else None

    q_expected, c_expected, p_expected, _ = reference_module(
        feats=feats_serial,
        s_trunk=s_trunk_serial,
        z=z_serial,
    )

    # Upstream gradients
    d_q = torch.empty_like(q_expected)
    d_c = torch.empty_like(c_expected)
    d_p = torch.empty_like(p_expected)
    init_tensors_uniform([d_q, d_c, d_p], low=val_init_min_max[0], high=val_init_min_max[1])

    # Apply masks to upstream gradients to zero invalid positions
    mask_expanded = feats_serial["atom_pad_mask"].unsqueeze(-1)
    d_q = d_q * mask_expanded
    d_c = d_c * mask_expanded

    # Mask d_p with pair mask (query AND key masks) -- matches V1x test pattern
    compute_dtype = torch.promote_types(dtype, torch.float32)
    index_matrix = get_indexing_matrix(K, W, H, device_type).to(dtype=compute_dtype)
    to_keys_fn_serial = partial(single_to_keys, indexing_matrix=index_matrix, W=W, H=H)
    mask_key_serial = to_keys_fn_serial(
        feats_serial["atom_pad_mask"].to(dtype=compute_dtype, device=d_p.device).unsqueeze(-1)
    )
    # d_p: (B, K, W, H, atom_z) * mask_key_serial: (B, K, 1, H, 1) → (B, K, W, H, atom_z)
    d_p = d_p * mask_key_serial[:, :, None, :, :]

    torch.autograd.backward([q_expected, c_expected, p_expected], [d_q, d_c, d_p])

    expected_param_grads = {
        name: param.grad.detach().cpu() for name, param in reference_module.named_parameters() if param.grad is not None
    }

    spawn_multiprocessing(
        parallel_assert_atom_encoder_wb,
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
        structure_prediction,
        {
            k: v.detach().cpu() for k, v in layer_state_dict.items()
        },  # CPU state dict avoids cross-process CUDA IPC issues
        {k: v.detach().cpu() for k, v in feats.items()},
        s_trunk.detach().cpu() if s_trunk is not None else None,
        z.detach().cpu() if z is not None else None,
        q_expected.detach().cpu(),
        c_expected.detach().cpu(),
        p_expected.detach().cpu(),
        d_q.detach().cpu(),
        d_c.detach().cpu(),
        d_p.detach().cpu(),
        s_trunk_serial.grad.detach().cpu() if s_trunk_serial is not None and s_trunk_serial.grad is not None else None,
        z_serial.grad.detach().cpu() if z_serial is not None and z_serial.grad is not None else None,
        expected_param_grads,
    )


# ======================================================================
# Test 2: AtomEncoder under autocast bf16 (dtype-only comparison)
# ======================================================================


def parallel_assert_atom_encoder_wb_autocast_bf16(
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
    layer_state_dict,
    feats_global_host: dict[str, torch.Tensor],
    s_trunk_global_host: torch.Tensor,
    z_global_host: torch.Tensor,
    q_serial_dtype: torch.dtype,
    c_serial_dtype: torch.dtype,
    p_serial_dtype: torch.dtype,
    serial_grad_dtypes: dict[str, torch.dtype],
    serial_param_grad_dtypes: dict[str, torch.dtype],
):
    """Parallel worker for bf16 autocast dtype test.

    Runs DTensor AtomEncoder forward + backward under
    torch.autocast("cuda", dtype=torch.bfloat16) and asserts output and
    gradient dtypes match the serial reference (computed in main process).
    """
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

    module_serial = SerialAtomEncoderV2(
        atom_s=atom_s,
        atom_z=atom_z,
        token_s=token_s,
        token_z=token_z,
        atoms_per_window_queries=W,
        atoms_per_window_keys=H,
        atom_feature_dim=atom_feature_dim,
        structure_prediction=True,
    )
    module_serial = module_serial.to(device=manager.device, dtype=dtype)
    module_serial.load_state_dict(layer_state_dict)

    module_dt = DistributedAtomEncoder(
        layer=module_serial,
        device_mesh=manager.device_mesh_subgroups,
    ).train()

    # Distribute atom features (same pattern as existing test)
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

    # DTensor forward + backward under autocast
    with torch.autocast("cuda", dtype=torch.bfloat16):
        q_dt, c_dt, p_dt = module_dt(
            feats=feats_dt_packed,
            s_trunk=s_trunk_dt,
            z=z_dt,
        )

    torch.autograd.backward(
        [q_dt, c_dt, p_dt],
        [torch.ones_like(q_dt), torch.ones_like(c_dt), torch.ones_like(p_dt)],
    )

    # Assert forward output dtypes match serial reference
    assert q_dt.dtype == q_serial_dtype, f"q dtype mismatch: DTensor {q_dt.dtype} vs serial {q_serial_dtype}"
    assert c_dt.dtype == c_serial_dtype, f"c dtype mismatch: DTensor {c_dt.dtype} vs serial {c_serial_dtype}"
    assert p_dt.dtype == p_serial_dtype, f"p dtype mismatch: DTensor {p_dt.dtype} vs serial {p_serial_dtype}"

    # Assert input gradient dtypes match serial reference
    assert s_trunk_dt.grad is not None, "s_trunk_dt.grad is None"
    assert z_dt.grad is not None, "z_dt.grad is None"
    assert (
        s_trunk_dt.grad.dtype == serial_grad_dtypes["s_trunk"]
    ), f"s_trunk grad dtype mismatch: DTensor {s_trunk_dt.grad.dtype} vs serial {serial_grad_dtypes['s_trunk']}"
    assert (
        z_dt.grad.dtype == serial_grad_dtypes["z"]
    ), f"z grad dtype mismatch: DTensor {z_dt.grad.dtype} vs serial {serial_grad_dtypes['z']}"

    # Assert parameter gradient dtypes match serial reference
    for name, param in module_dt.named_parameters():
        if name in serial_param_grad_dtypes:
            grad = param.grad
            if grad is None:
                continue
            if hasattr(grad, "full_tensor"):
                grad_dtype = grad.full_tensor().dtype
            else:
                grad_dtype = grad.dtype
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
def test_atom_encoder_wb_autocast_bf16(setup_env):
    """Test DTensor AtomEncoder output dtypes under autocast bf16.

    Verifies that DTensor AtomEncoder produces the same output dtypes as the
    V2 serial AtomEncoder when both run under torch.autocast("cuda", dtype=torch.bfloat16).
    Uses dp=1, cp=(1,1) (1 GPU) since this is a dtype consistency test, not a CP correctness test.
    """
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() < {world_size}")

    seed = 42
    seed_by_rank(0, seed=seed)

    B = 1 * grid_group_sizes["dp"]
    W = 32
    H = 128
    val_init_min_max = (-0.1, 0.1)
    dtype = torch.float32

    n_atoms_per_token_min = 8
    n_atoms_per_token_max = 20
    size_cp = grid_group_sizes["cp"][0]
    N_tokens = 30 * size_cp
    N_atoms_raw = N_tokens * n_atoms_per_token_max
    N_atoms = ((N_atoms_raw + W - 1) // W) * W
    N_msa = 1

    atom_s = 8
    atom_z = 8
    token_s = 4
    token_z = 4

    from boltz.data import const as boltz_const

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
    feats = {k: v.to(dtype=dtype) if v.dtype.is_floating_point else v for k, v in feats.items()}

    s_trunk = torch.empty((B, N_tokens, token_s), device=device_type, dtype=dtype)
    z = torch.empty((B, N_tokens, N_tokens, token_z), device=device_type, dtype=dtype)
    init_tensors_uniform([s_trunk, z], low=val_init_min_max[0], high=val_init_min_max[1])

    reference_module = SerialAtomEncoderV2(
        atom_s=atom_s,
        atom_z=atom_z,
        token_s=token_s,
        token_z=token_z,
        atoms_per_window_queries=W,
        atoms_per_window_keys=H,
        atom_feature_dim=atom_feature_dim,
        structure_prediction=True,
    ).to(device=device_type, dtype=dtype)
    reference_module.eval()
    init_module_params_uniform(reference_module, low=val_init_min_max[0], high=val_init_min_max[1])
    reference_module.apply(SetModuleInfValues())
    layer_state_dict = reference_module.state_dict()

    # Serial forward + backward under autocast (in main process)
    reference_module.train()
    s_trunk_serial = s_trunk.detach().clone().requires_grad_(True)
    z_serial = z.detach().clone().requires_grad_(True)

    with torch.autocast("cuda", dtype=torch.bfloat16):
        q_serial, c_serial, p_serial, _ = reference_module(
            feats={k: v.clone() for k, v in feats.items()},
            s_trunk=s_trunk_serial,
            z=z_serial,
        )

    torch.autograd.backward(
        [q_serial, c_serial, p_serial],
        [torch.ones_like(q_serial), torch.ones_like(c_serial), torch.ones_like(p_serial)],
    )

    serial_grad_dtypes = {
        "s_trunk": s_trunk_serial.grad.dtype,
        "z": z_serial.grad.dtype,
    }
    serial_param_grad_dtypes = {}
    for name, param in reference_module.named_parameters():
        if param.grad is not None:
            serial_param_grad_dtypes[name] = param.grad.dtype

    spawn_multiprocessing(
        parallel_assert_atom_encoder_wb_autocast_bf16,
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
        {k: v.detach().cpu() for k, v in layer_state_dict.items()},
        {k: v.detach().cpu() for k, v in feats.items()},
        s_trunk.detach().cpu(),
        z.detach().cpu(),
        q_serial.dtype,
        c_serial.dtype,
        p_serial.dtype,
        serial_grad_dtypes,
        serial_param_grad_dtypes,
    )
