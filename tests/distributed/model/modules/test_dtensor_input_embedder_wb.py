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

"""Tests for DTensor InputEmbedder for Boltz-2.

Tests the distributed InputEmbedder forward and backward passes using window-
batching atom attention, with the serial Boltz-2 InputEmbedder as reference.

Verification:
    V8: multi-proc FW output tensor values close-to single-proc
    V10: multi-proc parameter gradient values close-to single-proc
"""

import pytest
import torch

from boltz.data import const
from boltz.distributed.data.utils import distribute_features
from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.modules.trunkv2 import InputEmbedder as DistInputEmbedder
from boltz.model.modules.trunkv2 import InputEmbedder as SerialInputEmbedder
from boltz.testing.utils import (
    SetModuleInfValues,
    assert_all_identical,
    distribute_atom_features,
    get_feature_placements,
    get_param_by_key,
    init_module_params_uniform,
    init_tensors_uniform,
    random_features,
    seed_by_rank,
    spawn_multiprocessing,
)

_selected_token_keys = {
    "token_pad_mask",
    "res_type",
}
_selected_msa_keys = {
    "profile",
    "deletion_mean",
}
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

# Additional token-level keys required when add_extra_feats=True
_extra_token_keys = {
    "method_feature",
    "modified",
    "cyclic_period",
    "mol_type",
}


def _get_placements(add_extra_feats: bool):
    """Build placement dicts, optionally including extra conditioning keys."""
    # mol_type and cyclic_period are known to get_feature_placements;
    # method_feature and modified are not, so we add them manually.
    known_extra = {"mol_type", "cyclic_period"}
    token_keys = _selected_token_keys | (known_extra if add_extra_feats else set())
    placements = get_feature_placements(
        token_keys=token_keys,
        msa_keys=_selected_msa_keys,
        atom_keys=_selected_atom_keys,
        model_io_keys=set(),
        model_io_fp32_keys=set(),
    )
    if add_extra_feats:
        placements_single = placements["single"]
        for key in ("method_feature", "modified"):
            placements["token_features"][key] = placements_single
    return placements


# Boltz-2 atom_feature_dim: 3 (ref_pos) + 1 (ref_charge) + 128 (ref_element) + 256 (ref_atom_name_chars)
_ATOM_FEATURE_DIM = 388


def parallel_assert_input_embedder(
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
    W: int,
    H: int,
    atom_encoder_depth: int,
    atom_encoder_heads: int,
    add_extra_feats: bool,
    layer_state_dict,
    feats_global_host: dict[str, torch.Tensor],
    d_s_global_host: torch.Tensor,
    s_expected_global_host: torch.Tensor,
    expected_param_grads_global_host_dict: dict[str, torch.Tensor],
):
    """Parallel worker for testing DTensor InputEmbedder forward and backward."""
    monkeypatch = pytest.MonkeyPatch()
    if env_per_rank is not None:
        for var_name, value in env_per_rank.items():
            if value == "<INPUT_RANK>":
                monkeypatch.setenv(var_name, f"{rank}")
                continue
            monkeypatch.setenv(var_name, value)

    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    manager = DistributedManager()

    try:
        placements = _get_placements(add_extra_feats)
        placements_token_features = placements["token_features"]
        placements_msa_features = placements["msa_features"]
        placements_cp_atom_features = placements["cp_atom_features"]
        placements_atom_features = placements["atom_features"]
        placements_single = placements["single"]

        selected_token_msa_keys = (
            _selected_token_keys | _selected_msa_keys | (_extra_token_keys if add_extra_feats else set())
        )

        # Recreate serial module from state dict
        module_serial = SerialInputEmbedder(
            atom_s=atom_s,
            atom_z=atom_z,
            token_s=token_s,
            token_z=token_z,
            atoms_per_window_queries=W,
            atoms_per_window_keys=H,
            atom_feature_dim=_ATOM_FEATURE_DIM,
            atom_encoder_depth=atom_encoder_depth,
            atom_encoder_heads=atom_encoder_heads,
            add_method_conditioning=add_extra_feats,
            add_modified_flag=add_extra_feats,
            add_cyclic_flag=add_extra_feats,
            add_mol_type_feat=add_extra_feats,
        )
        module_serial = module_serial.to(device=manager.device, dtype=dtype)
        module_serial.load_state_dict(layer_state_dict)
        module_serial = module_serial.train()
        module_serial.apply(SetModuleInfValues())

        # Create distributed module
        module = DistInputEmbedder(
            module=module_serial,
            device_mesh=manager.device_mesh_subgroups,
        ).train()

        # Get token_pad_mask for valid-region comparison
        token_pad_mask_global = feats_global_host["token_pad_mask"].to(device=manager.device, dtype=torch.bool)
        token_pad_mask_expanded_global = token_pad_mask_global.unsqueeze(-1)

        # ====================================================================
        # Distribute token and MSA features
        # ====================================================================
        if manager.group_rank["world"] == 0:
            input_feats_token_msa_global = {
                k: v.to(device=manager.device, dtype=dtype if v.dtype.is_floating_point else v.dtype)
                for k, v in feats_global_host.items()
                if k in selected_token_msa_keys
            }
        else:
            input_feats_token_msa_global = None

        feats_token_msa = distribute_features(
            input_feats_token_msa_global,
            placements_token_features | placements_msa_features,
            manager.group["world"],
            manager.group_ranks["world"][0],
            manager.device_mesh_subgroups,
        )

        # ====================================================================
        # Distribute atom features via distribute_atom_features
        # ====================================================================
        inputs_atom = {
            k: v.to(dtype=dtype if v.dtype.is_floating_point else v.dtype)
            for k, v in feats_global_host.items()
            if k in placements_cp_atom_features
        }

        feats_atom = distribute_atom_features(
            inputs_atom,
            placements_cp_atom_features,
            placements_atom_features,
            manager.device_mesh_subgroups,
            manager.group["cp"],
        )

        # ====================================================================
        # Merge all features (pack_atom_features is internalized in the module)
        # ====================================================================
        feats_dt = {**feats_token_msa, **feats_atom}

        # The serial AtomAttentionEncoder computes atom_to_token_mean with an
        # epsilon bias (+ 1e-6) while the distributed scatter_reduce("mean")
        # uses exact division.  This causes ~1e-5 absolute error in FP64.
        tol = {}
        if dtype == torch.float64:
            tol = {"atol": 5e-5, "rtol": 1e-4}

        # ====================================================================
        # Forward pass
        # ====================================================================
        s_dt = module(feats_dt)

        s_dt_full = s_dt.full_tensor()
        s_expected_device = s_expected_global_host.to(device=manager.device, dtype=dtype)
        torch.testing.assert_close(
            s_dt_full * token_pad_mask_expanded_global,
            s_expected_device * token_pad_mask_expanded_global,
            **tol,
        )

        # ====================================================================
        # Backward pass
        # ====================================================================
        d_s_expected_dtensor = distribute_features(
            {"d_s": d_s_global_host.to(device=manager.device, dtype=dtype)}
            if manager.group_rank["world"] == 0
            else None,
            {"d_s": placements_single},
            manager.group["world"],
            manager.group_ranks["world"][0],
            manager.device_mesh_subgroups,
        )["d_s"]

        s_dt.backward(d_s_expected_dtensor)

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

            torch.testing.assert_close(grad_global_host, grad_expected_global.to(dtype=dtype), **tol)
            assert_all_identical(grad_to_check, manager.group["cp"])

    finally:
        DistributedManager.cleanup()
        monkeypatch.undo()


@pytest.mark.slow
@pytest.mark.parametrize(
    "setup_env, dtype, add_extra_feats",
    (
        params_test := [
            # 2 GPU Test
            (((2, (1, 1)), True, "cuda", "ENV"), torch.float32, False),
            (((2, (1, 1)), True, "cuda", "ENV"), torch.float32, True),
            # # 8 GPU Test
            (((2, (2, 2)), True, "cuda", "ENV"), torch.float32, False),
            (((2, (2, 2)), True, "cuda", "ENV"), torch.float64, False),
            (((2, (2, 2)), True, "cuda", "ENV"), torch.float64, True),
        ]
    ),
    indirect=["setup_env"],
    ids=[f"dp:{x[0][0][0]}, cp:{x[0][0][1]}, dtype:{x[1]}, extra_feats:{x[2]}" for x in params_test],
)
def test_input_embedder_window_batching(setup_env, dtype, add_extra_feats):
    """Test DTensor InputEmbedder forward and backward for Boltz-2."""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    seed = 42
    seed_by_rank(0, seed=seed)

    size_cp = grid_group_sizes["cp"][0]
    B = 1 * grid_group_sizes["dp"]

    val_init_min_max_dtype = {torch.float64: (-0.08, 0.08), torch.float32: (-0.02, 0.02)}
    val_init_min_max = val_init_min_max_dtype[dtype]

    selected_keys = list(_selected_token_keys | _selected_msa_keys | _selected_atom_keys)
    if add_extra_feats:
        selected_keys.extend(["mol_type", "cyclic_period"])

    W = 32
    H = 128

    # n_atoms_per_token range chosen so that N_atoms > W, ensuring there is
    # interspersed padding going into the parallel data path to actually test
    # the pack_atom_features code path inside the distributed InputEmbedder.
    n_atoms_per_token_min = 8
    n_atoms_per_token_max = 20
    N_tokens = 30 * size_cp
    N_atoms = (N_tokens * n_atoms_per_token_max + W - 1) // W * W
    # MSA features are not used in the InputEmbedder
    # so we set N_msa to 1 (independent of size_cp)
    N_msa = 1

    assert N_tokens % size_cp == 0, f"N_tokens ({N_tokens}) must be divisible by size_cp ({size_cp})"

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
    assert feats["token_pad_mask"].shape == (B, N_tokens)

    if add_extra_feats:
        feats["method_feature"] = torch.randint(
            0, const.num_method_types, (B, N_tokens), device=feats["res_type"].device
        )
        feats["modified"] = torch.randint(0, 2, (B, N_tokens), device=feats["res_type"].device)

    atom_s = 8
    atom_z = 8
    token_s = 2
    token_z = 2
    atom_encoder_depth = 2
    atom_encoder_heads = 2

    feats = {k: v.to(dtype=dtype) if v.dtype.is_floating_point else v for k, v in feats.items()}

    reference_module = SerialInputEmbedder(
        atom_s=atom_s,
        atom_z=atom_z,
        token_s=token_s,
        token_z=token_z,
        atoms_per_window_queries=W,
        atoms_per_window_keys=H,
        atom_feature_dim=_ATOM_FEATURE_DIM,
        atom_encoder_depth=atom_encoder_depth,
        atom_encoder_heads=atom_encoder_heads,
        add_method_conditioning=add_extra_feats,
        add_modified_flag=add_extra_feats,
        add_cyclic_flag=add_extra_feats,
        add_mol_type_feat=add_extra_feats,
    ).to(device=device_type, dtype=dtype)
    reference_module.train()

    init_module_params_uniform(reference_module, low=val_init_min_max[0], high=val_init_min_max[1])
    reference_module.apply(SetModuleInfValues())

    layer_state_dict = reference_module.state_dict()

    # Serial forward pass
    feats_serial = {k: v.detach().clone() for k, v in feats.items()}
    s_expected = reference_module(feats_serial)
    s_expected_global_host = s_expected.detach().cpu()

    # Serial backward pass
    d_s = torch.empty_like(s_expected)
    init_tensors_uniform([d_s], low=val_init_min_max[0], high=val_init_min_max[1])
    d_s = d_s * feats_serial["token_pad_mask"].unsqueeze(-1)
    s_expected.backward(d_s)
    d_s_global_host = d_s.detach().cpu()
    expected_param_grads_global_host_dict = {
        name: param.grad.detach().cpu() for name, param in reference_module.named_parameters() if param.grad is not None
    }

    feats_global_host = {k: v.detach().cpu() for k, v in feats.items()}

    spawn_multiprocessing(
        parallel_assert_input_embedder,
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
        W,
        H,
        atom_encoder_depth,
        atom_encoder_heads,
        add_extra_feats,
        layer_state_dict,
        feats_global_host,
        d_s_global_host,
        s_expected_global_host,
        expected_param_grads_global_host_dict,
    )


# ========================================================================
# Standalone BF16 mixed-precision test (via torch.autocast)
# ========================================================================


def parallel_assert_input_embedder_bf16(
    rank: int,
    grid_group_sizes,
    device_type: str,
    backend: str,
    env_per_rank,
    atom_s: int,
    atom_z: int,
    token_s: int,
    token_z: int,
    W: int,
    H: int,
    atom_encoder_depth: int,
    atom_encoder_heads: int,
    add_extra_feats: bool,
    layer_state_dict,
    feats_global_host: dict[str, torch.Tensor],
    d_s_global_host: torch.Tensor,
    s_expected_global_host: torch.Tensor,
    serial_output_dtype: torch.dtype,
    serial_param_grad_dtypes: dict[str, torch.dtype],
):
    """Parallel worker for BF16 mixed-precision test of DTensor InputEmbedder.

    Simulates production BF16 mixed-precision training by wrapping the
    forward pass in ``torch.autocast("cuda", dtype=torch.bfloat16)``.
    Module weights and input features are FP32; autocast handles precision
    inside eligible operations.  The distributed AtomEncoder internally
    disables autocast (matching serial behavior) for numerical stability.

    Checks:
    - Output dtype matches the serial reference output dtype.
    - Output values are close to serial reference (with mixed-precision tolerance).
    - Parameter gradient dtypes match serial gradient dtypes.
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

    try:
        placements = _get_placements(add_extra_feats)
        placements_token_features = placements["token_features"]
        placements_msa_features = placements["msa_features"]
        placements_cp_atom_features = placements["cp_atom_features"]
        placements_atom_features = placements["atom_features"]
        placements_single = placements["single"]

        selected_token_msa_keys = (
            _selected_token_keys | _selected_msa_keys | (_extra_token_keys if add_extra_feats else set())
        )

        # Recreate serial module from FP32 state dict (production mixed-precision setup)
        module_serial = SerialInputEmbedder(
            atom_s=atom_s,
            atom_z=atom_z,
            token_s=token_s,
            token_z=token_z,
            atoms_per_window_queries=W,
            atoms_per_window_keys=H,
            atom_feature_dim=_ATOM_FEATURE_DIM,
            atom_encoder_depth=atom_encoder_depth,
            atom_encoder_heads=atom_encoder_heads,
            add_method_conditioning=add_extra_feats,
            add_modified_flag=add_extra_feats,
            add_cyclic_flag=add_extra_feats,
            add_mol_type_feat=add_extra_feats,
        )
        module_serial = module_serial.to(device=manager.device, dtype=torch.float32)
        module_serial.load_state_dict(layer_state_dict)
        module_serial = module_serial.train()
        module_serial.apply(SetModuleInfValues())

        # AtomEncoder internally disables autocast, matching serial behavior
        module = DistInputEmbedder(
            module=module_serial,
            device_mesh=manager.device_mesh_subgroups,
        ).train()

        # Get token_pad_mask for valid-region comparison
        token_pad_mask_global = feats_global_host["token_pad_mask"].to(device=manager.device, dtype=torch.bool)
        token_pad_mask_expanded_global = token_pad_mask_global.unsqueeze(-1)

        # ====================================================================
        # Distribute token and MSA features (FP32)
        # ====================================================================
        if manager.group_rank["world"] == 0:
            input_feats_token_msa_global = {
                k: v.to(device=manager.device, dtype=torch.float32 if v.dtype.is_floating_point else v.dtype)
                for k, v in feats_global_host.items()
                if k in selected_token_msa_keys
            }
        else:
            input_feats_token_msa_global = None

        feats_token_msa = distribute_features(
            input_feats_token_msa_global,
            placements_token_features | placements_msa_features,
            manager.group["world"],
            manager.group_ranks["world"][0],
            manager.device_mesh_subgroups,
        )

        # ====================================================================
        # Distribute atom features (FP32)
        # ====================================================================
        inputs_atom = {
            k: v.to(dtype=torch.float32 if v.dtype.is_floating_point else v.dtype)
            for k, v in feats_global_host.items()
            if k in placements_cp_atom_features
        }

        feats_atom = distribute_atom_features(
            inputs_atom,
            placements_cp_atom_features,
            placements_atom_features,
            manager.device_mesh_subgroups,
            manager.group["cp"],
        )

        feats_dt = {**feats_token_msa, **feats_atom}

        # ====================================================================
        # Forward pass under autocast (mirrors production training)
        # ====================================================================
        with torch.autocast("cuda", dtype=torch.bfloat16):
            s_dt = module(feats_dt)

        s_dt_full = s_dt.full_tensor()

        # Dtype check: distributed output dtype must match serial
        assert (
            s_dt_full.dtype == serial_output_dtype
        ), f"Distributed output dtype ({s_dt_full.dtype}) != serial output dtype ({serial_output_dtype})"

        # Value check with mixed-precision tolerance
        s_expected_device = s_expected_global_host.to(device=manager.device)
        compare_dtype = torch.promote_types(s_dt_full.dtype, s_expected_device.dtype)
        torch.testing.assert_close(
            (s_dt_full * token_pad_mask_expanded_global).to(compare_dtype),
            (s_expected_device * token_pad_mask_expanded_global).to(compare_dtype),
            atol=0.05,
            rtol=0.05,
        )

        # ====================================================================
        # Backward pass (outside autocast, matching production training)
        # ====================================================================
        d_s_expected_dtensor = distribute_features(
            {"d_s": d_s_global_host.to(device=manager.device)} if manager.group_rank["world"] == 0 else None,
            {"d_s": placements_single},
            manager.group["world"],
            manager.group_ranks["world"][0],
            manager.device_mesh_subgroups,
        )["d_s"]

        s_dt.backward(d_s_expected_dtensor)

        # Gradient dtype check: each parameter's grad dtype must match serial
        for name, expected_dtype in serial_param_grad_dtypes.items():
            grad_param = get_param_by_key(module, name).grad
            assert grad_param is not None, f"Missing grad for param {name}"
            actual_dtype = grad_param.full_tensor().dtype if hasattr(grad_param, "full_tensor") else grad_param.dtype
            assert (
                actual_dtype == expected_dtype
            ), f"Grad dtype mismatch for {name}: distributed={actual_dtype}, serial={expected_dtype}"

    finally:
        DistributedManager.cleanup()
        monkeypatch.undo()


@pytest.mark.slow
@pytest.mark.parametrize(
    "setup_env, add_extra_feats",
    (
        params_test_bf16 := [
            (((2, (2, 2)), True, "cuda", "ENV"), False),
            (((2, (2, 2)), True, "cuda", "ENV"), True),
        ]
    ),
    indirect=["setup_env"],
    ids=[f"dp:{x[0][0][0]}, cp:{x[0][0][1]}, extra_feats:{x[1]}" for x in params_test_bf16],
)
def test_input_embedder_bf16(setup_env, add_extra_feats):
    """BF16 mixed-precision test for DTensor InputEmbedder.

    Simulates production training with ``torch.autocast("cuda", dtype=torch.bfloat16)``
    wrapping the forward pass.  Module weights and input features are FP32;
    autocast handles precision internally.  Verifies that the distributed and
    serial modules produce matching output dtypes and gradient dtypes.
    """
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("CUDA not available")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"Need {world_size} GPUs, have {torch.cuda.device_count()}")

    seed = 42
    seed_by_rank(0, seed=seed)

    size_cp = grid_group_sizes["cp"][0]
    B = 1 * grid_group_sizes["dp"]
    val_init_min_max = (-0.02, 0.02)

    selected_keys = list(_selected_token_keys | _selected_msa_keys | _selected_atom_keys)
    if add_extra_feats:
        selected_keys.extend(["mol_type", "cyclic_period"])

    W = 32
    H = 128
    n_atoms_per_token_min = 8
    n_atoms_per_token_max = 20
    N_tokens = 30 * size_cp
    N_atoms = (N_tokens * n_atoms_per_token_max + W - 1) // W * W
    N_msa = 1

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
    # random_features returns float64; cast to FP32 for production-like mixed precision
    feats = {k: v.to(torch.float32) if v.is_floating_point() else v for k, v in feats.items()}
    assert feats["token_pad_mask"].shape == (B, N_tokens)

    if add_extra_feats:
        feats["method_feature"] = torch.randint(
            0, const.num_method_types, (B, N_tokens), device=feats["res_type"].device
        )
        feats["modified"] = torch.randint(0, 2, (B, N_tokens), device=feats["res_type"].device)

    atom_s = 8
    atom_z = 8
    token_s = 2
    token_z = 2
    atom_encoder_depth = 2
    atom_encoder_heads = 2

    reference_module = SerialInputEmbedder(
        atom_s=atom_s,
        atom_z=atom_z,
        token_s=token_s,
        token_z=token_z,
        atoms_per_window_queries=W,
        atoms_per_window_keys=H,
        atom_feature_dim=_ATOM_FEATURE_DIM,
        atom_encoder_depth=atom_encoder_depth,
        atom_encoder_heads=atom_encoder_heads,
        add_method_conditioning=add_extra_feats,
        add_modified_flag=add_extra_feats,
        add_cyclic_flag=add_extra_feats,
        add_mol_type_feat=add_extra_feats,
    ).to(device=device_type, dtype=torch.float32)
    reference_module.train()

    init_module_params_uniform(reference_module, low=val_init_min_max[0], high=val_init_min_max[1])
    reference_module.apply(SetModuleInfValues())

    layer_state_dict = reference_module.state_dict()

    # Serial forward with autocast (mirrors production training)
    feats_serial = {k: v.detach().clone() for k, v in feats.items()}
    with torch.autocast("cuda", dtype=torch.bfloat16):
        s_expected = reference_module(feats_serial)
    serial_output_dtype = s_expected.dtype

    # Serial backward (outside autocast)
    d_s = torch.empty_like(s_expected)
    init_tensors_uniform([d_s], low=val_init_min_max[0], high=val_init_min_max[1])
    d_s = d_s * feats_serial["token_pad_mask"].unsqueeze(-1)
    s_expected.backward(d_s)
    d_s_global_host = d_s.detach().cpu()

    serial_param_grad_dtypes = {
        name: param.grad.dtype for name, param in reference_module.named_parameters() if param.grad is not None
    }

    s_expected_global_host = s_expected.detach().cpu()
    feats_global_host = {k: v.detach().cpu() for k, v in feats.items()}

    spawn_multiprocessing(
        parallel_assert_input_embedder_bf16,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        atom_s,
        atom_z,
        token_s,
        token_z,
        W,
        H,
        atom_encoder_depth,
        atom_encoder_heads,
        add_extra_feats,
        layer_state_dict,
        feats_global_host,
        d_s_global_host,
        s_expected_global_host,
        serial_output_dtype,
        serial_param_grad_dtypes,
    )
