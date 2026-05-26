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

"""Tests for DTensor DiffusionModule with window batching.

Tests the DTensor DiffusionModule against both V1 and V2 serial references,
verifying forward and backward numerical equivalence.

Parametrized on ``internalized_conditioning``:
- True (V1): module owns pairwise_conditioner, forward takes z_trunk + relative_position_encoding
- False (V2): forward takes pre-computed diffusion_conditioning dict

Uses float64 with default tolerance for exact comparison.
"""

import pytest
import torch
from torch.distributed.tensor import Replicate, Shard, distribute_tensor

from boltz.data import const as boltz_const
from boltz.distributed.comm import AttentionPairBiasComm
from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.modules.diffusion import DiffusionModule as DistributedDiffusionModule
from boltz.distributed.model.modules.diffusion_conditioning import (
    DiffusionConditioning as DistributedDiffusionConditioning,
)
from boltz.model.modules.diffusion import DiffusionModule as SerialDiffusionModuleV1
from boltz.model.modules.diffusion_conditioning import DiffusionConditioning as SerialDiffusionConditioning
from boltz.model.modules.diffusionv2 import DiffusionModule as SerialDiffusionModuleV2
from boltz.testing.utils import (
    SetModuleInfValues,
    distribute_atom_features,
    get_feature_placements,
    get_param_by_key,
    init_module_params_uniform,
    init_tensors_uniform,
    random_features,
    seed_by_rank,
    spawn_multiprocessing,
)

# Atom features needed
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
# Token features needed
_selected_token_keys = {"token_pad_mask"}

_selected_model_io_keys = {
    "r_noisy_expected",
    "r_update_expected",
    "d_r_update_expected",
    "d_r_noisy_expected",
}

_placements = get_feature_placements(
    token_keys=set(),
    msa_keys=set(),
    atom_keys=_selected_atom_keys,
    model_io_keys=_selected_model_io_keys,
    model_io_fp32_keys=set(),
)
_placements_cp_atom_features = _placements["cp_atom_features"]
_placements_atom_features = _placements["atom_features"]
_placements_model_io = _placements["model_io"]
_placements_cp_model_io = _placements["cp_model_io"]


def parallel_assert_diffusion_module(
    rank: int,
    grid_group_sizes,
    device_type: str,
    backend: str,
    env_per_rank,
    internalized_conditioning: bool,
    dtype: torch.dtype,
    multiplicity: int,
    # Module params
    diffusion_module_state_dict,
    diffusion_conditioning_state_dict,  # None for internalized
    module_kwargs: dict,
    conditioning_kwargs: dict | None,  # None for internalized
    W: int,
    H: int,
    # Inputs
    feats_global_host: dict[str, torch.Tensor],
    s_inputs_global_host: torch.Tensor,
    s_trunk_global_host: torch.Tensor,
    z_trunk_global_host: torch.Tensor,
    rel_pos_enc_global_host: torch.Tensor,
    r_noisy_global_host: torch.Tensor,
    times_global_host: torch.Tensor,
    # Expected outputs
    r_update_expected_global_host: torch.Tensor,
    # Upstream grad
    d_r_update_global_host: torch.Tensor,
    # Expected input grads
    d_s_inputs_expected_global_host: torch.Tensor,
    d_s_trunk_expected_global_host: torch.Tensor,
    d_r_noisy_expected_global_host: torch.Tensor,
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

    # Re-create serial DiffusionModule
    if internalized_conditioning:
        serial_diffusion_module = SerialDiffusionModuleV1(**module_kwargs)
    else:
        serial_diffusion_module = SerialDiffusionModuleV2(**module_kwargs)
    serial_diffusion_module = serial_diffusion_module.to(device=manager.device, dtype=dtype)
    serial_diffusion_module.load_state_dict(diffusion_module_state_dict)
    serial_diffusion_module = serial_diffusion_module.train()

    # Create ring_comm for the token-level transformer
    ring_comm = AttentionPairBiasComm(
        manager.group["cp"],
        manager.layout_subgroups["cp"],
        manager.subgroups["cp"][0],
        manager.subgroups["cp"][1],
    )

    # Create DTensor module
    module = DistributedDiffusionModule(
        layer=serial_diffusion_module,
        device_mesh=manager.device_mesh_subgroups,
        ring_comm=ring_comm,
    ).train()

    # ------------------------------------------------------------------
    # Distribute token-level tensors (common to both paths)
    # ------------------------------------------------------------------
    placements_single = (Shard(0), Shard(1), Replicate())
    placements_pair = (Shard(0), Shard(1), Shard(2))
    placements_times = (Shard(0), Replicate(), Replicate())

    token_pad_mask_dt = distribute_tensor(
        feats_global_host["token_pad_mask"].to(device=manager.device, dtype=dtype),
        manager.device_mesh_subgroups,
        placements_single,
    )

    s_inputs_dt = distribute_tensor(
        s_inputs_global_host.to(device=manager.device, dtype=dtype),
        manager.device_mesh_subgroups,
        placements_single,
    ).requires_grad_(True)
    s_trunk_dt = distribute_tensor(
        s_trunk_global_host.to(device=manager.device, dtype=dtype),
        manager.device_mesh_subgroups,
        placements_single,
    ).requires_grad_(True)

    times_dt = distribute_tensor(
        times_global_host.to(device=manager.device, dtype=dtype),
        manager.device_mesh_subgroups,
        placements_times,
    ).requires_grad_(False)

    z_trunk_device = z_trunk_global_host.to(device=manager.device, dtype=dtype)
    rel_pos_enc_device = rel_pos_enc_global_host.to(device=manager.device, dtype=dtype)

    # ------------------------------------------------------------------
    # Distribute atom features and r_noisy (shared by V1 and V2)
    # ------------------------------------------------------------------
    # Both V1 and V2 pass unpacked feats — DiffusionModule packs internally.
    # All atom-level I/O must share the same intersperse-padded atom ordering.
    inputs_atom = {
        k: v.to(dtype=dtype if v.dtype.is_floating_point else v.dtype)
        for k, v in feats_global_host.items()
        if k in _placements_cp_atom_features
    }
    size_batch = feats_global_host["atom_pad_mask"].shape[0]
    io_tensors = {
        "r_noisy_expected": r_noisy_global_host,
        "r_update_expected": r_update_expected_global_host,
        "d_r_update_expected": d_r_update_global_host,
        "d_r_noisy_expected": d_r_noisy_expected_global_host,
    }
    for base_name, tensor_host in io_tensors.items():
        unflat = tensor_host.unflatten(0, (size_batch, multiplicity))
        for i_mul in range(multiplicity):
            inputs_atom[f"{base_name}_{i_mul}"] = unflat[:, i_mul].to(dtype=dtype)

    placements_cp_io_mul = {
        f"{k}_{i_mul}": v for k, v in _placements_cp_model_io.items() for i_mul in range(multiplicity)
    }
    placements_io_mul = {f"{k}_{i_mul}": v for k, v in _placements_model_io.items() for i_mul in range(multiplicity)}
    multiplicities = dict.fromkeys(io_tensors, multiplicity)

    feats_and_io = distribute_atom_features(
        inputs_atom,
        _placements_cp_atom_features | placements_cp_io_mul,
        _placements_atom_features | placements_io_mul,
        manager.device_mesh_subgroups,
        manager.group["cp"],
        multiplicities=multiplicities,
    )
    r_noisy_dt = feats_and_io.pop("r_noisy_expected").requires_grad_(True)
    r_update_expected_dt = feats_and_io.pop("r_update_expected")
    d_r_update_dt_expected = feats_and_io.pop("d_r_update_expected")
    d_r_noisy_expected_dt = feats_and_io.pop("d_r_noisy_expected")
    feats_dt = feats_and_io
    feats_dt["token_pad_mask"] = token_pad_mask_dt

    z_trunk_dt = distribute_tensor(z_trunk_device, manager.device_mesh_subgroups, placements_pair)
    rel_pos_enc_dt = distribute_tensor(rel_pos_enc_device, manager.device_mesh_subgroups, placements_pair)

    # ------------------------------------------------------------------
    # Forward pass (depends on internalized_conditioning)
    # ------------------------------------------------------------------
    if internalized_conditioning:
        # V1: pass z_trunk and rel_pos_enc directly
        r_update_result = module(
            s_inputs=s_inputs_dt,
            s_trunk=s_trunk_dt,
            r_noisy=r_noisy_dt,
            times=times_dt,
            feats=feats_dt,
            z_trunk=z_trunk_dt,
            relative_position_encoding=rel_pos_enc_dt,
            multiplicity=multiplicity,
        )
        r_update_dt = r_update_result["r_update"]
    else:
        # V2: DTensor DiffusionConditioning takes unpacked feats (packs internally)
        serial_conditioning = SerialDiffusionConditioning(**conditioning_kwargs)
        serial_conditioning = serial_conditioning.to(device=manager.device, dtype=dtype)
        serial_conditioning.load_state_dict(diffusion_conditioning_state_dict)
        serial_conditioning = serial_conditioning.eval()

        dtensor_conditioning = DistributedDiffusionConditioning(
            layer=serial_conditioning,
            device_mesh=manager.device_mesh_subgroups,
        ).eval()

        with torch.no_grad():
            q_cond_dt, c_cond_dt, atom_enc_bias_dt, atom_dec_bias_dt, token_trans_bias_dt = dtensor_conditioning(
                s_trunk=s_trunk_dt.detach(),
                z_trunk=z_trunk_dt,
                relative_position_encoding=rel_pos_enc_dt,
                feats=feats_dt,
            )

        diff_cond_dt = {
            "q": q_cond_dt,
            "c": c_cond_dt,
            "atom_enc_bias": atom_enc_bias_dt,
            "atom_dec_bias": atom_dec_bias_dt,
            "token_trans_bias": token_trans_bias_dt,
        }

        r_update_dt = module(
            s_inputs=s_inputs_dt,
            s_trunk=s_trunk_dt,
            r_noisy=r_noisy_dt,
            times=times_dt,
            feats=feats_dt,
            diffusion_conditioning=diff_cond_dt,
            multiplicity=multiplicity,
        )

    # ------------------------------------------------------------------
    # Forward comparison (both V1 and V2 outputs are in intersperse-padded layout)
    # ------------------------------------------------------------------
    torch.testing.assert_close(r_update_dt.full_tensor(), r_update_expected_dt.full_tensor())

    # ------------------------------------------------------------------
    # Backward pass (upstream grad in intersperse-padded layout via distribute_atom_features)
    # ------------------------------------------------------------------
    r_update_dt.backward(d_r_update_dt_expected)

    # Check input gradients (token-level)
    torch.testing.assert_close(
        s_inputs_dt.grad.full_tensor(),
        d_s_inputs_expected_global_host.to(device=manager.device, dtype=dtype),
    )
    torch.testing.assert_close(
        s_trunk_dt.grad.full_tensor(),
        d_s_trunk_expected_global_host.to(device=manager.device, dtype=dtype),
    )

    # r_noisy grad comparison (atom-level, intersperse-padded layout)
    torch.testing.assert_close(r_noisy_dt.grad.full_tensor(), d_r_noisy_expected_dt.full_tensor())

    # Parameter grads
    for name, grad_expected_global in expected_param_grads_global_host_dict.items():
        grad_param = get_param_by_key(module, name).grad
        if grad_param is None:
            continue
        if hasattr(grad_param, "full_tensor"):
            grad_global_host = grad_param.full_tensor().cpu()
        else:
            grad_global_host = grad_param.detach().cpu()
        torch.testing.assert_close(
            grad_global_host,
            grad_expected_global.to(dtype=dtype),
            msg=lambda m: f"Parameter gradient mismatch for {name}: {m}",
        )

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.slow
@pytest.mark.parametrize(
    "setup_env",
    [
        ((2, (2, 2)), True, "cuda", "ENV"),
    ],
    indirect=["setup_env"],
    ids=lambda x: f"dp:{x[0][0]}, cp:{x[0][1]}, device_type:{x[2]}",
)
@pytest.mark.parametrize("multiplicity", [1, 4], ids=lambda x: f"mul:{x}")
@pytest.mark.parametrize("internalized_conditioning", [False, True], ids=["extern", "intern"])
def test_diffusion_module(setup_env, multiplicity, internalized_conditioning: bool):
    """Test DTensor DiffusionModule with window batching.

    Parametrized on ``internalized_conditioning``:
    - False (V2 / externalized): uses DiffusionConditioning to pre-compute q/c/bias, float64
    - True (V1 / internalized): passes z_trunk + relative_position_encoding directly, float32
      (V1 serial code uses hardcoded .float() casts incompatible with float64)
    """
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env
    # V1 serial uses .float() casts internally → must test with float32
    # V2 serial uses promote_types → can test with float64 for exact comparison
    dtype = torch.float32 if internalized_conditioning else torch.float64

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
    val_init_min_max = (-0.5, 0.5)

    n_atoms_per_token_min = 8
    n_atoms_per_token_max = 20
    N_tokens = 30 * size_cp
    N_atoms_raw = N_tokens * n_atoms_per_token_max
    N_atoms = ((N_atoms_raw + W - 1) // W) * W
    N_msa = 1

    atom_s = 8
    token_s = 4
    token_z = 4
    atom_z = 8

    atom_encoder_depth = 2
    atom_encoder_heads = 2
    token_transformer_depth = 2
    token_transformer_heads = 2
    atom_decoder_depth = 2
    atom_decoder_heads = 2
    conditioning_transition_layers = 1

    # V1: ref_pos(3) + ref_charge(1) + atom_pad_mask(1) + ref_element + ref_atom_name_chars(4*64)
    # V2: ref_pos(3) + ref_charge(1) + ref_element + ref_atom_name_chars(4*64)  (no atom_pad_mask)
    atom_feature_dim = 3 + 1 + (1 if internalized_conditioning else 0) + boltz_const.num_elements + 4 * 64

    selected_keys = list(_selected_atom_keys | _selected_token_keys)

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

    # Token-level inputs
    # V1 s_inputs has wider dim: input_dim - token_s, where
    # input_dim = 2 * token_s + 2 * num_tokens + 1 + len(pocket_contact_info)
    # V2 s_inputs has token_s dim (s_trunk and s_inputs are same width)
    if internalized_conditioning:
        v1_input_dim = 2 * token_s + 2 * boltz_const.num_tokens + 1 + len(boltz_const.pocket_contact_info)
        s_inputs_dim = v1_input_dim - token_s
    else:
        s_inputs_dim = token_s
    s_inputs = torch.empty((B, N_tokens, s_inputs_dim), device=device_type, dtype=dtype, requires_grad=True)
    s_trunk = torch.empty((B, N_tokens, token_s), device=device_type, dtype=dtype, requires_grad=True)
    z_trunk = torch.empty((B, N_tokens, N_tokens, token_z), device=device_type, dtype=dtype)
    rel_pos_enc = torch.empty((B, N_tokens, N_tokens, token_z), device=device_type, dtype=dtype)
    init_tensors_uniform([s_inputs, s_trunk, z_trunk, rel_pos_enc], low=val_init_min_max[0], high=val_init_min_max[1])

    # r_noisy: (B*M, N_atoms, 3)
    N_atoms_actual = feats["atom_pad_mask"].shape[1]
    r_noisy = torch.empty((B * multiplicity, N_atoms_actual, 3), device=device_type, dtype=dtype, requires_grad=True)
    times = torch.empty((B * multiplicity,), device=device_type, dtype=dtype)
    init_tensors_uniform([r_noisy, times], low=val_init_min_max[0], high=val_init_min_max[1])

    # ------------------------------------------------------------------
    # Build serial modules and compute reference (depends on internalized_conditioning)
    # ------------------------------------------------------------------
    if internalized_conditioning:
        # V1: module owns pairwise_conditioner, encoder computes q/c/p internally
        module_kwargs = {
            "token_s": token_s,
            "token_z": token_z,
            "atom_s": atom_s,
            "atom_z": atom_z,
            "atoms_per_window_queries": W,
            "atoms_per_window_keys": H,
            "sigma_data": 16,
            "dim_fourier": 32,
            "atom_encoder_depth": atom_encoder_depth,
            "atom_encoder_heads": atom_encoder_heads,
            "token_transformer_depth": token_transformer_depth,
            "token_transformer_heads": token_transformer_heads,
            "atom_decoder_depth": atom_decoder_depth,
            "atom_decoder_heads": atom_decoder_heads,
            "atom_feature_dim": atom_feature_dim,
            "conditioning_transition_layers": conditioning_transition_layers,
        }
        serial_module = SerialDiffusionModuleV1(**module_kwargs).to(device=device_type, dtype=dtype)
        serial_module.train()
        init_module_params_uniform(serial_module, low=val_init_min_max[0], high=val_init_min_max[1])
        serial_module.apply(SetModuleInfValues())
        module_state_dict = serial_module.state_dict()

        conditioning_kwargs = None
        conditioning_state_dict = None

        # Serial forward: V1 takes z_trunk + relative_position_encoding directly
        feats_serial = {k: v.detach().clone() for k, v in feats.items()}
        s_inputs_serial = s_inputs.detach().clone().requires_grad_(True)
        s_trunk_serial = s_trunk.detach().clone().requires_grad_(True)
        r_noisy_serial = r_noisy.detach().clone().requires_grad_(True)

        result_serial = serial_module(
            s_inputs=s_inputs_serial,
            s_trunk=s_trunk_serial,
            z_trunk=z_trunk.detach(),
            r_noisy=r_noisy_serial,
            times=times.detach(),
            relative_position_encoding=rel_pos_enc.detach(),
            feats=feats_serial,
            multiplicity=multiplicity,
        )
        r_update_serial = result_serial["r_update"]

    else:
        # V2: uses DiffusionConditioning to pre-compute conditioning
        module_kwargs = {
            "token_s": token_s,
            "atom_s": atom_s,
            "atoms_per_window_queries": W,
            "atoms_per_window_keys": H,
            "sigma_data": 16,
            "dim_fourier": 32,
            "atom_encoder_depth": atom_encoder_depth,
            "atom_encoder_heads": atom_encoder_heads,
            "token_transformer_depth": token_transformer_depth,
            "token_transformer_heads": token_transformer_heads,
            "atom_decoder_depth": atom_decoder_depth,
            "atom_decoder_heads": atom_decoder_heads,
            "conditioning_transition_layers": conditioning_transition_layers,
        }
        serial_module = SerialDiffusionModuleV2(**module_kwargs).to(device=device_type, dtype=dtype)
        serial_module.train()
        init_module_params_uniform(serial_module, low=val_init_min_max[0], high=val_init_min_max[1])
        serial_module.apply(SetModuleInfValues())
        module_state_dict = serial_module.state_dict()

        conditioning_kwargs = {
            "token_s": token_s,
            "token_z": token_z,
            "atom_s": atom_s,
            "atom_z": atom_z,
            "atoms_per_window_queries": W,
            "atoms_per_window_keys": H,
            "atom_encoder_depth": atom_encoder_depth,
            "atom_encoder_heads": atom_encoder_heads,
            "token_transformer_depth": token_transformer_depth,
            "token_transformer_heads": token_transformer_heads,
            "atom_decoder_depth": atom_decoder_depth,
            "atom_decoder_heads": atom_decoder_heads,
            "atom_feature_dim": atom_feature_dim,
            "conditioning_transition_layers": conditioning_transition_layers,
        }
        serial_conditioning = SerialDiffusionConditioning(**conditioning_kwargs).to(device=device_type, dtype=dtype)
        serial_conditioning.train()
        init_module_params_uniform(serial_conditioning, low=val_init_min_max[0], high=val_init_min_max[1])
        serial_conditioning.apply(SetModuleInfValues())
        conditioning_state_dict = serial_conditioning.state_dict()

        # Serial forward: first conditioning, then diffusion module
        feats_serial = {k: v.detach().clone() for k, v in feats.items()}
        s_inputs_serial = s_inputs.detach().clone().requires_grad_(True)
        s_trunk_serial = s_trunk.detach().clone().requires_grad_(True)
        r_noisy_serial = r_noisy.detach().clone().requires_grad_(True)

        with torch.no_grad():
            q_cond, c_cond, to_keys, atom_enc_bias_cond, atom_dec_bias_cond, token_trans_bias_cond = (
                serial_conditioning(
                    s_trunk=s_trunk.detach(),
                    z_trunk=z_trunk.detach(),
                    relative_position_encoding=rel_pos_enc.detach(),
                    feats={k: v.detach() for k, v in feats.items()},
                )
            )

        diff_cond_serial = {
            "q": q_cond.detach(),
            "c": c_cond.detach(),
            "to_keys": to_keys,
            "atom_enc_bias": atom_enc_bias_cond.detach(),
            "atom_dec_bias": atom_dec_bias_cond.detach(),
            "token_trans_bias": token_trans_bias_cond.detach(),
        }

        r_update_serial = serial_module(
            s_inputs=s_inputs_serial,
            s_trunk=s_trunk_serial,
            r_noisy=r_noisy_serial,
            times=times.detach(),
            feats=feats_serial,
            diffusion_conditioning=diff_cond_serial,
            multiplicity=multiplicity,
        )

    # Upstream gradient
    d_r_update = torch.empty_like(r_update_serial)
    init_tensors_uniform([d_r_update], low=val_init_min_max[0], high=val_init_min_max[1])
    atom_mask_mul = feats_serial["atom_pad_mask"].repeat_interleave(multiplicity, 0).unsqueeze(-1)
    d_r_update = d_r_update * atom_mask_mul

    r_update_serial.backward(d_r_update)

    expected_param_grads = {
        name: param.grad.detach().cpu()
        for name, param in serial_module.named_parameters()
        if param.requires_grad and param.grad is not None
    }

    spawn_multiprocessing(
        parallel_assert_diffusion_module,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        internalized_conditioning,
        dtype,
        multiplicity,
        {k: v.detach().cpu() for k, v in module_state_dict.items()},
        {k: v.detach().cpu() for k, v in conditioning_state_dict.items()} if conditioning_state_dict else None,
        module_kwargs,
        conditioning_kwargs,
        W,
        H,
        {k: v.detach().cpu() for k, v in feats.items()},
        s_inputs.detach().cpu(),
        s_trunk.detach().cpu(),
        z_trunk.detach().cpu(),
        rel_pos_enc.detach().cpu(),
        r_noisy.detach().cpu(),
        times.detach().cpu(),
        r_update_serial.detach().cpu(),
        d_r_update.detach().cpu(),
        s_inputs_serial.grad.detach().cpu(),
        s_trunk_serial.grad.detach().cpu(),
        r_noisy_serial.grad.detach().cpu(),
        expected_param_grads,
    )
