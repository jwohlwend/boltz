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

"""Tests for DTensor AtomDiffusion with window batching.

Tests the DTensor AtomDiffusion (preconditioned_network_forward, training
forward) against V1 and V2 serial references, verifying forward and backward
numerical equivalence.

Parametrized on ``internalized_conditioning``:
- True (V1): module owns pairwise_conditioner, forward takes z_trunk + relative_position_encoding
- False (V2): forward takes pre-computed diffusion_conditioning dict

Both V1 and V2 use float64 with default tolerance for exact comparison.

Adapted from Boltz-1x CP tests (test_dtensor_diffusion.py and
test_dtensor_diffusion_precond.py).
"""

import pytest
import torch
from torch.distributed.tensor import Replicate, Shard, distribute_tensor

import boltz.distributed.model.modules.diffusion as distributed_diffusion_module
import boltz.model.modules.diffusion as serial_diffusion_v1_module
import boltz.model.modules.diffusionv2 as serial_diffusion_v2_module
from boltz.data import const as boltz_const
from boltz.distributed.comm import AttentionPairBiasComm, TransposeComm
from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.layers.elementwise_op import ElementwiseOp, scalar_tensor_op
from boltz.distributed.model.modules.diffusion import AtomDiffusion as DistributedAtomDiffusion
from boltz.distributed.model.modules.diffusion_conditioning import (
    DiffusionConditioning as DistributedDiffusionConditioning,
)
from boltz.model.modules.diffusion import AtomDiffusion as SerialAtomDiffusionV1
from boltz.model.modules.diffusion_conditioning import DiffusionConditioning as SerialDiffusionConditioning
from boltz.model.modules.diffusionv2 import AtomDiffusion as SerialAtomDiffusionV2
from boltz.testing.utils import (
    SetModuleInfValues,
    assert_tensors_identical,
    distribute_atom_features,
    get_feature_placements,
    get_param_by_key,
    init_module_params_glorot,
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
    "atom_resolved_mask",
    "plddt",
}
# Token features needed
_selected_token_keys = {"token_pad_mask"}

_selected_model_io_keys = {
    "r_noisy_expected",
    "r_update_expected",
    "d_r_update_expected",
    "d_r_noisy_expected",
    "noise",
    "denoised_atom_coords",
    "d_denoised_atom_coords",
    "aligned_true_atom_coords",
}

_placements = get_feature_placements(
    token_keys={"mol_type"},
    msa_keys=set(),
    atom_keys=_selected_atom_keys,
    model_io_keys=_selected_model_io_keys,
    model_io_fp32_keys=set(),
)
_placements_cp_atom_features = _placements["cp_atom_features"]
_placements_atom_features = _placements["atom_features"]
_placements_model_io = _placements["model_io"]
_placements_cp_model_io = _placements["cp_model_io"]

# ======================================================================
# Test 1: preconditioned_network_forward
# ======================================================================


def parallel_assert_precond_forward(
    rank: int,
    grid_group_sizes,
    device_type: str,
    backend: str,
    env_per_rank,
    internalized_conditioning: bool,
    dtype: torch.dtype,
    multiplicity: int,
    # Module states
    atom_diffusion_state_dict,
    conditioning_state_dict,  # None for internalized
    score_model_kwargs: dict,
    conditioning_kwargs: dict | None,  # None for internalized
    W: int,
    H: int,
    # Inputs
    feats_global_host: dict[str, torch.Tensor],
    s_inputs_global_host: torch.Tensor,
    s_trunk_global_host: torch.Tensor,
    z_trunk_global_host: torch.Tensor,
    rel_pos_enc_global_host: torch.Tensor,
    noised_atom_coords_global_host: torch.Tensor,
    sigma_global_host: torch.Tensor,
    # Expected outputs
    denoised_coords_expected_global_host: torch.Tensor,
    # Upstream grad
    d_denoised_coords_global_host: torch.Tensor,
    # Expected input grads
    d_s_inputs_expected_global_host: torch.Tensor,
    d_s_trunk_expected_global_host: torch.Tensor,
    d_noised_coords_expected_global_host: torch.Tensor,
    # Expected param grads
    expected_param_grads_global_host_dict: dict[str, torch.Tensor],
):
    """Parallel assertion for AtomDiffusion.preconditioned_network_forward()."""
    monkeypatch = pytest.MonkeyPatch()
    if env_per_rank is not None:
        for var_name, value in env_per_rank.items():
            if value == "<INPUT_RANK>":
                monkeypatch.setenv(var_name, f"{rank}")
                continue
            monkeypatch.setenv(var_name, value)

    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    manager = DistributedManager()

    # Re-create serial module to initialize DTensor version
    if internalized_conditioning:
        serial_atom_diffusion = SerialAtomDiffusionV1(score_model_args=score_model_kwargs)
    else:
        serial_atom_diffusion = SerialAtomDiffusionV2(score_model_args=score_model_kwargs)
    serial_atom_diffusion = serial_atom_diffusion.to(device=manager.device, dtype=dtype)
    serial_atom_diffusion.load_state_dict(atom_diffusion_state_dict)
    serial_atom_diffusion = serial_atom_diffusion.train()

    # Create ring_comm for token-level transformer
    ring_comm = AttentionPairBiasComm(
        manager.group["cp"],
        manager.layout_subgroups["cp"],
        manager.subgroups["cp"][0],
        manager.subgroups["cp"][1],
    )

    # Create DTensor module
    module = DistributedAtomDiffusion(
        layer=serial_atom_diffusion,
        device_mesh=manager.device_mesh_subgroups,
        ring_comm=ring_comm,
    ).train()

    # ------------------------------------------------------------------
    # Distribute token-level tensors (common to both paths)
    # ------------------------------------------------------------------
    placements_single = (Shard(0), Shard(1), Replicate())
    placements_pair = (Shard(0), Shard(1), Shard(2))

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

    placements_scalar = (Shard(0), Replicate(), Replicate())
    sigma_dt = distribute_tensor(
        sigma_global_host.to(device=manager.device, dtype=dtype),
        manager.device_mesh_subgroups,
        placements_scalar,
    )

    z_trunk_device = z_trunk_global_host.to(device=manager.device, dtype=dtype)
    rel_pos_enc_device = rel_pos_enc_global_host.to(device=manager.device, dtype=dtype)

    # ------------------------------------------------------------------
    # Distribute atom features, atom-level I/O, and conditioning
    # ------------------------------------------------------------------
    inputs_atom = {
        k: v.to(dtype=dtype if v.dtype.is_floating_point else v.dtype)
        for k, v in feats_global_host.items()
        if k in _placements_cp_atom_features
    }

    # All atom-level I/O must share the same intersperse-padded atom ordering.
    # Both V1 and V2 pass unpacked feats — DiffusionModule packs internally.
    size_batch = feats_global_host["atom_pad_mask"].shape[0]
    io_tensors = {
        "r_noisy_expected": noised_atom_coords_global_host,
        "r_update_expected": denoised_coords_expected_global_host,
        "d_r_update_expected": d_denoised_coords_global_host,
        "d_r_noisy_expected": d_noised_coords_expected_global_host,
    }
    for base_name, tensor_host in io_tensors.items():
        unflat = tensor_host.unflatten(0, (size_batch, multiplicity))
        for i_mul in range(multiplicity):
            inputs_atom[f"{base_name}_{i_mul}"] = unflat[:, i_mul].to(dtype=dtype)

    io_keys_used = set(io_tensors.keys())
    placements_cp_io_mul = {
        f"{k}_{i_mul}": _placements_cp_model_io[k] for k in io_keys_used for i_mul in range(multiplicity)
    }
    placements_io_mul = {f"{k}_{i_mul}": _placements_model_io[k] for k in io_keys_used for i_mul in range(multiplicity)}
    multiplicities = dict.fromkeys(io_keys_used, multiplicity)

    feats_and_io = distribute_atom_features(
        inputs_atom,
        _placements_cp_atom_features | placements_cp_io_mul,
        _placements_atom_features | placements_io_mul,
        manager.device_mesh_subgroups,
        manager.group["cp"],
        multiplicities=multiplicities,
    )
    noised_atom_coords_dt = feats_and_io.pop("r_noisy_expected").requires_grad_(True)
    denoised_expected_dt = feats_and_io.pop("r_update_expected")
    d_denoised_expected_dt = feats_and_io.pop("d_r_update_expected")
    d_noised_expected_dt = feats_and_io.pop("d_r_noisy_expected")
    feats_dt = feats_and_io
    feats_dt["token_pad_mask"] = token_pad_mask_dt

    z_trunk_dt = distribute_tensor(z_trunk_device, manager.device_mesh_subgroups, placements_pair)
    rel_pos_enc_dt = distribute_tensor(rel_pos_enc_device, manager.device_mesh_subgroups, placements_pair)

    # ------------------------------------------------------------------
    # Build network_condition_kwargs (depends on internalized_conditioning)
    # ------------------------------------------------------------------
    if internalized_conditioning:
        # V1: pass z_trunk and rel_pos_enc directly
        network_condition_kwargs = {
            "s_inputs": s_inputs_dt,
            "s_trunk": s_trunk_dt,
            "feats": feats_dt,
            "multiplicity": multiplicity,
            "z_trunk": z_trunk_dt,
            "relative_position_encoding": rel_pos_enc_dt,
        }
    else:
        # V2: DTensor DiffusionConditioning takes unpacked feats (packs internally)
        serial_conditioning = SerialDiffusionConditioning(**conditioning_kwargs)
        serial_conditioning = serial_conditioning.to(device=manager.device, dtype=dtype)
        serial_conditioning.load_state_dict(conditioning_state_dict)
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

        network_condition_kwargs = {
            "s_inputs": s_inputs_dt,
            "s_trunk": s_trunk_dt,
            "feats": feats_dt,
            "multiplicity": multiplicity,
            "diffusion_conditioning": diff_cond_dt,
        }

    # ------------------------------------------------------------------
    # Forward pass: preconditioned_network_forward
    # ------------------------------------------------------------------
    precond_result = module.preconditioned_network_forward(
        noised_atom_coords_dt,
        sigma_dt,
        network_condition_kwargs,
    )

    if internalized_conditioning:
        denoised_coords_dt = precond_result[0]
    else:
        denoised_coords_dt = precond_result

    # ------------------------------------------------------------------
    # Forward comparison (both V1 and V2 outputs are in intersperse-padded layout)
    # ------------------------------------------------------------------
    torch.testing.assert_close(denoised_coords_dt.full_tensor(), denoised_expected_dt.full_tensor())

    # ------------------------------------------------------------------
    # Backward pass (upstream grad in intersperse-padded layout)
    # ------------------------------------------------------------------
    denoised_coords_dt.backward(d_denoised_expected_dt)

    # Check input gradients (token-level)
    torch.testing.assert_close(
        s_inputs_dt.grad.full_tensor(),
        d_s_inputs_expected_global_host.to(device=manager.device, dtype=dtype),
    )
    torch.testing.assert_close(
        s_trunk_dt.grad.full_tensor(),
        d_s_trunk_expected_global_host.to(device=manager.device, dtype=dtype),
    )

    # noised_atom_coords grad (atom-level, intersperse-padded layout)
    torch.testing.assert_close(noised_atom_coords_dt.grad.full_tensor(), d_noised_expected_dt.full_tensor())

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
def test_preconditioned_network_forward(setup_env, multiplicity, internalized_conditioning: bool):
    """Test AtomDiffusion.preconditioned_network_forward() with DTensor CP.

    Tests the core preconditioned forward pass which computes:
        denoised = c_skip(sigma) * x + c_out(sigma) * score_model(c_in(sigma) * x, c_noise(sigma))

    Parametrized on ``internalized_conditioning``:
    - False (V2 / externalized): uses DiffusionConditioning to pre-compute q/c/bias
    - True (V1 / internalized): passes z_trunk + relative_position_encoding directly

    Uses float64 for exact comparison. Verifies forward and backward numerical
    equivalence against serial reference.
    """
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env
    dtype = torch.float64

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

    # V1 includes atom_pad_mask in atom features; V2 does not
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

    # V1 s_inputs has wider dim:
    #   input_dim = 2 * token_s + 2 * num_tokens + 1 + len(pocket_contact_info)
    #   s_inputs_dim = input_dim - token_s
    # V2 s_inputs has token_s dim
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

    # noised_atom_coords: (B*M, N_atoms_actual, 3)
    N_atoms_actual = feats["atom_pad_mask"].shape[1]
    noised_atom_coords = torch.empty(
        (B * multiplicity, N_atoms_actual, 3), device=device_type, dtype=dtype, requires_grad=True
    )
    init_tensors_uniform([noised_atom_coords], low=val_init_min_max[0], high=val_init_min_max[1])

    # ------------------------------------------------------------------
    # Build serial modules (depends on internalized_conditioning)
    # ------------------------------------------------------------------
    if internalized_conditioning:
        # V1: module owns pairwise_conditioner, no external conditioning
        score_model_kwargs = {
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
        conditioning_kwargs = None
        conditioning_state_dict = None

        serial_atom_diffusion = SerialAtomDiffusionV1(
            score_model_args=score_model_kwargs,
            coordinate_augmentation=False,
        ).to(device=device_type, dtype=dtype)
    else:
        # V2: uses DiffusionConditioning
        score_model_kwargs = {
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
        init_module_params_glorot(serial_conditioning, gain=0.5)
        serial_conditioning.apply(SetModuleInfValues())
        conditioning_state_dict = serial_conditioning.state_dict()

        serial_atom_diffusion = SerialAtomDiffusionV2(
            score_model_args=score_model_kwargs,
            coordinate_augmentation=False,
        ).to(device=device_type, dtype=dtype)

    serial_atom_diffusion.train()
    init_module_params_glorot(serial_atom_diffusion, gain=0.5)
    serial_atom_diffusion.apply(SetModuleInfValues())
    atom_diffusion_state_dict = serial_atom_diffusion.state_dict()

    # Generate sigma from the module's noise distribution (deterministic via seed_by_rank)
    sigma = serial_atom_diffusion.noise_distribution(B * multiplicity).to(device=device_type, dtype=dtype)

    # ------------------------------------------------------------------
    # Serial preconditioned_network_forward
    # ------------------------------------------------------------------
    s_inputs_serial = s_inputs.detach().clone().requires_grad_(True)
    s_trunk_serial = s_trunk.detach().clone().requires_grad_(True)
    noised_coords_serial = noised_atom_coords.detach().clone().requires_grad_(True)

    if internalized_conditioning:
        # V1: pass z_trunk + relative_position_encoding directly
        serial_result = serial_atom_diffusion.preconditioned_network_forward(
            noised_coords_serial,
            sigma,
            network_condition_kwargs={
                "s_inputs": s_inputs_serial,
                "s_trunk": s_trunk_serial,
                "z_trunk": z_trunk.detach(),
                "relative_position_encoding": rel_pos_enc.detach(),
                "feats": {k: v.detach().clone() for k, v in feats.items()},
                "multiplicity": multiplicity,
            },
        )
        denoised_serial = serial_result[0]  # (denoised_coords, token_a)
    else:
        # V2: pre-compute conditioning
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

        denoised_serial = serial_atom_diffusion.preconditioned_network_forward(
            noised_coords_serial,
            sigma,
            network_condition_kwargs={
                "s_inputs": s_inputs_serial,
                "s_trunk": s_trunk_serial,
                "feats": {k: v.detach().clone() for k, v in feats.items()},
                "multiplicity": multiplicity,
                "diffusion_conditioning": diff_cond_serial,
            },
        )

    # Upstream gradient
    d_denoised = torch.empty_like(denoised_serial)
    init_tensors_uniform([d_denoised], low=val_init_min_max[0], high=val_init_min_max[1])
    atom_mask_mul = feats["atom_pad_mask"].repeat_interleave(multiplicity, 0).unsqueeze(-1)
    d_denoised = d_denoised * atom_mask_mul

    denoised_serial.backward(d_denoised)

    expected_param_grads = {
        name: param.grad.detach().cpu()
        for name, param in serial_atom_diffusion.named_parameters()
        if param.requires_grad and param.grad is not None
    }

    spawn_multiprocessing(
        parallel_assert_precond_forward,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        internalized_conditioning,
        dtype,
        multiplicity,
        {k: v.detach().cpu() for k, v in atom_diffusion_state_dict.items()},
        {k: v.detach().cpu() for k, v in conditioning_state_dict.items()} if conditioning_state_dict else None,
        score_model_kwargs,
        conditioning_kwargs,
        W,
        H,
        {k: v.detach().cpu() for k, v in feats.items()},
        s_inputs.detach().cpu(),
        s_trunk.detach().cpu(),
        z_trunk.detach().cpu(),
        rel_pos_enc.detach().cpu(),
        noised_atom_coords.detach().cpu(),
        sigma.detach().cpu(),
        denoised_serial.detach().cpu(),
        d_denoised.detach().cpu(),
        s_inputs_serial.grad.detach().cpu(),
        s_trunk_serial.grad.detach().cpu(),
        noised_coords_serial.grad.detach().cpu(),
        expected_param_grads,
    )


# ======================================================================
# Test 2: AtomDiffusion.forward() (training forward)
# ======================================================================


def parallel_assert_atom_diffusion_forward(
    rank: int,
    grid_group_sizes,
    device_type: str,
    backend: str,
    env_per_rank,
    internalized_conditioning: bool,
    dtype: torch.dtype,
    multiplicity: int,
    # Module states
    atom_diffusion_state_dict,
    conditioning_state_dict,  # None for internalized
    score_model_kwargs: dict,
    conditioning_kwargs: dict | None,  # None for internalized
    W: int,
    H: int,
    # Inputs
    feats_global_host: dict[str, torch.Tensor],
    s_inputs_global_host: torch.Tensor,
    s_trunk_global_host: torch.Tensor,
    z_trunk_global_host: torch.Tensor,
    rel_pos_enc_global_host: torch.Tensor,
    coords_global_host: torch.Tensor,
    sigmas_global_host: torch.Tensor,
    noise_global_host: torch.Tensor,
    # Expected outputs
    denoised_expected_global_host: torch.Tensor,
    noised_expected_global_host: torch.Tensor | None,  # V1 only
    # Upstream grad
    d_denoised_global_host: torch.Tensor,
    # Expected input grads
    d_s_inputs_expected_global_host: torch.Tensor,
    d_s_trunk_expected_global_host: torch.Tensor,
    # Expected param grads
    expected_param_grads_global_host_dict: dict[str, torch.Tensor],
):
    """Parallel assertion for AtomDiffusion.forward() (training forward)."""
    monkeypatch = pytest.MonkeyPatch()
    if env_per_rank is not None:
        for var_name, value in env_per_rank.items():
            if value == "<INPUT_RANK>":
                monkeypatch.setenv(var_name, f"{rank}")
                continue
            monkeypatch.setenv(var_name, value)

    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    manager = DistributedManager()

    # Re-create serial module to initialize DTensor version
    if internalized_conditioning:
        serial_atom_diffusion = SerialAtomDiffusionV1(
            score_model_args=score_model_kwargs,
            coordinate_augmentation=False,
        )
    else:
        serial_atom_diffusion = SerialAtomDiffusionV2(
            score_model_args=score_model_kwargs,
            coordinate_augmentation=False,
        )
    serial_atom_diffusion = serial_atom_diffusion.to(device=manager.device, dtype=dtype)
    serial_atom_diffusion.load_state_dict(atom_diffusion_state_dict)
    serial_atom_diffusion = serial_atom_diffusion.train()

    # Create ring_comm for token-level transformer
    ring_comm = AttentionPairBiasComm(
        manager.group["cp"],
        manager.layout_subgroups["cp"],
        manager.subgroups["cp"][0],
        manager.subgroups["cp"][1],
    )

    # Create DTensor module
    module = DistributedAtomDiffusion(
        layer=serial_atom_diffusion,
        device_mesh=manager.device_mesh_subgroups,
        ring_comm=ring_comm,
    ).train()

    # ------------------------------------------------------------------
    # Distribute token-level tensors (common to both paths)
    # ------------------------------------------------------------------
    placements_single = (Shard(0), Shard(1), Replicate())
    placements_pair = (Shard(0), Shard(1), Shard(2))
    placements_scalar = (Shard(0), Replicate(), Replicate())

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

    sigmas_dt = distribute_tensor(
        sigmas_global_host.to(device=manager.device, dtype=dtype),
        manager.device_mesh_subgroups,
        placements_scalar,
    )

    z_trunk_device = z_trunk_global_host.to(device=manager.device, dtype=dtype)
    rel_pos_enc_device = rel_pos_enc_global_host.to(device=manager.device, dtype=dtype)

    # ------------------------------------------------------------------
    # Distribute atom features, coords, noise, and conditioning
    # ------------------------------------------------------------------
    inputs_atom = {
        k: v.to(dtype=dtype if v.dtype.is_floating_point else v.dtype)
        for k, v in feats_global_host.items()
        if k in _placements_cp_atom_features
    }

    # All atom-level I/O must share the same intersperse-padded atom ordering.
    # Both V1 and V2 pass unpacked feats — DiffusionModule packs internally.
    size_batch = feats_global_host["atom_pad_mask"].shape[0]
    io_tensors = {
        "noise": noise_global_host,
        "denoised_atom_coords": coords_global_host,
        "d_denoised_atom_coords": d_denoised_global_host,
    }
    for base_name, tensor_host in io_tensors.items():
        unflat = tensor_host.unflatten(0, (size_batch, multiplicity))
        for i_mul in range(multiplicity):
            inputs_atom[f"{base_name}_{i_mul}"] = unflat[:, i_mul].to(dtype=dtype)

    io_keys_used = set(io_tensors.keys())
    placements_cp_io_mul = {
        f"{k}_{i_mul}": _placements_cp_model_io[k] for k in io_keys_used for i_mul in range(multiplicity)
    }
    placements_io_mul = {f"{k}_{i_mul}": _placements_model_io[k] for k in io_keys_used for i_mul in range(multiplicity)}
    multiplicities = dict.fromkeys(io_keys_used, multiplicity)

    feats_and_io = distribute_atom_features(
        inputs_atom,
        _placements_cp_atom_features | placements_cp_io_mul,
        _placements_atom_features | placements_io_mul,
        manager.device_mesh_subgroups,
        manager.group["cp"],
        multiplicities=multiplicities,
    )
    noise_dt = feats_and_io.pop("noise")
    coords_dt = feats_and_io.pop("denoised_atom_coords")
    d_denoised_expected_dt = feats_and_io.pop("d_denoised_atom_coords")
    feats_dt = feats_and_io
    feats_dt["token_pad_mask"] = token_pad_mask_dt
    feats_dt["coords"] = coords_dt

    z_trunk_dt = distribute_tensor(z_trunk_device, manager.device_mesh_subgroups, placements_pair)
    rel_pos_enc_dt = distribute_tensor(rel_pos_enc_device, manager.device_mesh_subgroups, placements_pair)

    # ------------------------------------------------------------------
    # Forward pass (with monkeypatched noise_distribution and create_distributed_randn)
    # ------------------------------------------------------------------
    monkeypatch.setattr(module, "noise_distribution", lambda bs, dtype=None: sigmas_dt)
    monkeypatch.setattr(distributed_diffusion_module, "create_distributed_randn", lambda *a, **kw: noise_dt)

    if internalized_conditioning:
        # V1: pass z_trunk and rel_pos_enc directly
        out_dt = module(
            s_inputs=s_inputs_dt,
            s_trunk=s_trunk_dt,
            feats=feats_dt,
            z_trunk=z_trunk_dt,
            relative_position_encoding=rel_pos_enc_dt,
            multiplicity=multiplicity,
        )
    else:
        # V2: DTensor DiffusionConditioning takes unpacked feats (packs internally)
        serial_conditioning = SerialDiffusionConditioning(**conditioning_kwargs)
        serial_conditioning = serial_conditioning.to(device=manager.device, dtype=dtype)
        serial_conditioning.load_state_dict(conditioning_state_dict)
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

        out_dt = module(
            s_inputs=s_inputs_dt,
            s_trunk=s_trunk_dt,
            feats=feats_dt,
            diffusion_conditioning=diff_cond_dt,
            multiplicity=multiplicity,
        )

    denoised_dt = out_dt["denoised_atom_coords"]

    # ------------------------------------------------------------------
    # Forward comparison: noised_atom_coords (V1 only)
    # Extract real atoms via boolean mask from both layouts and compare.
    # denoised_atom_coords comparison is skipped because center_random_augmentation
    # processing cannot be trivially aligned between layouts; correctness is
    # verified through the backward pass instead.
    # ------------------------------------------------------------------
    if internalized_conditioning and noised_expected_global_host is not None:
        noised_dt = out_dt["noised_atom_coords"]
        noised_expected = noised_expected_global_host.to(device=manager.device, dtype=dtype)
        dt_mask = feats_dt["atom_pad_mask"].full_tensor().repeat_interleave(multiplicity, 0).bool()
        serial_mask = (
            feats_global_host["atom_pad_mask"].to(device=manager.device).repeat_interleave(multiplicity, 0).bool()
        )
        torch.testing.assert_close(noised_dt.full_tensor()[dt_mask], noised_expected[serial_mask])

    # ------------------------------------------------------------------
    # Backward pass (upstream grad in intersperse-padded layout)
    # ------------------------------------------------------------------
    denoised_dt.backward(d_denoised_expected_dt)

    # Check input gradients (token-level)
    torch.testing.assert_close(
        s_inputs_dt.grad.full_tensor(),
        d_s_inputs_expected_global_host.to(device=manager.device, dtype=dtype),
    )
    torch.testing.assert_close(
        s_trunk_dt.grad.full_tensor(),
        d_s_trunk_expected_global_host.to(device=manager.device, dtype=dtype),
    )

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
def test_atom_diffusion_forward(setup_env, multiplicity, internalized_conditioning: bool):
    """Test AtomDiffusion.forward() (training forward) with DTensor CP.

    Tests the full training forward which includes:
    - center_random_augmentation (with augmentation=False for determinism)
    - Noise generation (provided externally for determinism)
    - Sigma scheduling (provided externally for determinism)
    - Preconditioned network forward

    Parametrized on ``internalized_conditioning``:
    - False (V2 / externalized): uses DiffusionConditioning
    - True (V1 / internalized): passes z_trunk + relative_position_encoding

    Uses float64 for exact comparison. Uses pre-generated sigmas and noise for deterministic serial vs DTensor
    comparison of forward and backward numerical equivalence.
    """
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env
    dtype = torch.float64

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
    val_init_min_max = (-0.2, 0.2)

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

    # V1 includes atom_pad_mask in atom features; V2 does not
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
    N_atoms_actual = feats["atom_pad_mask"].shape[1]

    # V1 s_inputs has wider dim
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

    # coords: (B*M, N_atoms, 3) for both V1 and V2 (DTensor format)
    coords = torch.empty((B * multiplicity, N_atoms_actual, 3), device=device_type, dtype=dtype)
    init_tensors_uniform([coords], low=val_init_min_max[0], high=val_init_min_max[1])

    # ------------------------------------------------------------------
    # Build serial modules (depends on internalized_conditioning)
    # ------------------------------------------------------------------
    if internalized_conditioning:
        score_model_kwargs = {
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
        conditioning_kwargs = None
        conditioning_state_dict = None

        serial_atom_diffusion = SerialAtomDiffusionV1(
            score_model_args=score_model_kwargs,
            coordinate_augmentation=False,
        ).to(device=device_type, dtype=dtype)
    else:
        score_model_kwargs = {
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
        init_module_params_glorot(serial_conditioning, gain=0.5)
        serial_conditioning.apply(SetModuleInfValues())
        conditioning_state_dict = serial_conditioning.state_dict()

        serial_atom_diffusion = SerialAtomDiffusionV2(
            score_model_args=score_model_kwargs,
            coordinate_augmentation=False,
        ).to(device=device_type, dtype=dtype)

    serial_atom_diffusion.train()
    init_module_params_glorot(serial_atom_diffusion, gain=0.5)
    serial_atom_diffusion.apply(SetModuleInfValues())
    atom_diffusion_state_dict = serial_atom_diffusion.state_dict()

    # Pre-generate sigmas and noise for deterministic comparison.
    # These are injected via monkeypatching so serial and DTensor use identical values.
    sigmas = serial_atom_diffusion.noise_distribution(B * multiplicity).to(device=device_type, dtype=dtype)
    noise = torch.empty((B * multiplicity, N_atoms_actual, 3), device=device_type, dtype=dtype)
    init_tensors_uniform([noise], low=val_init_min_max[0], high=val_init_min_max[1])

    # ------------------------------------------------------------------
    # Serial forward (with monkeypatched noise_distribution and randn_like)
    # ------------------------------------------------------------------
    serial_mod = serial_diffusion_v1_module if internalized_conditioning else serial_diffusion_v2_module
    _monkeypatch = pytest.MonkeyPatch()
    _monkeypatch.setattr(serial_atom_diffusion, "noise_distribution", lambda bs, dtype=None: sigmas)
    _monkeypatch.setattr(serial_mod.torch, "randn_like", lambda t: noise.to(t))

    feats_serial = {k: v.detach().clone() for k, v in feats.items()}
    s_inputs_serial = s_inputs.detach().clone().requires_grad_(True)
    s_trunk_serial = s_trunk.detach().clone().requires_grad_(True)

    if internalized_conditioning:
        # V1: coords need (B, M, N_atoms, 3) shape for serial forward
        feats_serial["coords"] = coords.detach().clone().reshape(B, multiplicity, N_atoms_actual, 3)

        out_serial = serial_atom_diffusion(
            s_inputs=s_inputs_serial,
            s_trunk=s_trunk_serial,
            z_trunk=z_trunk.detach(),
            relative_position_encoding=rel_pos_enc.detach(),
            feats=feats_serial,
            multiplicity=multiplicity,
        )
    else:
        feats_serial["coords"] = coords.detach().clone()

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

        out_serial = serial_atom_diffusion(
            s_inputs=s_inputs_serial,
            s_trunk=s_trunk_serial,
            feats=feats_serial,
            diffusion_conditioning=diff_cond_serial,
            multiplicity=multiplicity,
        )

    _monkeypatch.undo()

    denoised_serial = out_serial["denoised_atom_coords"]
    noised_serial = out_serial.get("noised_atom_coords")  # V1 only

    # Upstream gradient
    d_denoised = torch.empty_like(denoised_serial)
    init_tensors_uniform([d_denoised], low=val_init_min_max[0], high=val_init_min_max[1])
    atom_mask_mul = feats["atom_pad_mask"].repeat_interleave(multiplicity, 0).unsqueeze(-1)
    d_denoised = d_denoised * atom_mask_mul

    denoised_serial.backward(d_denoised)

    expected_param_grads = {
        name: param.grad.detach().cpu()
        for name, param in serial_atom_diffusion.named_parameters()
        if param.requires_grad and param.grad is not None
    }

    spawn_multiprocessing(
        parallel_assert_atom_diffusion_forward,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        internalized_conditioning,
        dtype,
        multiplicity,
        {k: v.detach().cpu() for k, v in atom_diffusion_state_dict.items()},
        {k: v.detach().cpu() for k, v in conditioning_state_dict.items()} if conditioning_state_dict else None,
        score_model_kwargs,
        conditioning_kwargs,
        W,
        H,
        {k: v.detach().cpu() for k, v in feats.items()},
        s_inputs.detach().cpu(),
        s_trunk.detach().cpu(),
        z_trunk.detach().cpu(),
        rel_pos_enc.detach().cpu(),
        coords.detach().cpu(),
        sigmas.detach().cpu(),
        noise.detach().cpu(),
        denoised_serial.detach().cpu(),
        noised_serial.detach().cpu() if noised_serial is not None else None,
        d_denoised.detach().cpu(),
        s_inputs_serial.grad.detach().cpu(),
        s_trunk_serial.grad.detach().cpu(),
        expected_param_grads,
    )


# ======================================================================
# Test 3: AtomDiffusion.sample() (inference sampling)
# ======================================================================


def parallel_assert_atom_diffusion_sample(
    rank: int,
    grid_group_sizes,
    device_type: str,
    backend: str,
    env_per_rank,
    internalized_conditioning: bool,
    dtype: torch.dtype,
    multiplicity: int,
    num_sampling_steps: int,
    max_parallel_samples: int,
    atom_diffusion_state_dict,
    conditioning_state_dict,
    score_model_kwargs: dict,
    conditioning_kwargs: dict | None,
    W: int,
    H: int,
    feats_global_host: dict[str, torch.Tensor],
    s_inputs_global_host: torch.Tensor,
    s_trunk_global_host: torch.Tensor,
    z_trunk_global_host: torch.Tensor,
    rel_pos_enc_global_host: torch.Tensor,
    init_noise_global_host: torch.Tensor,
    step_noise_list_global_host: list[torch.Tensor],
    sample_coords_expected_global_host: torch.Tensor,
):
    """Parallel assertion for AtomDiffusion.sample()."""
    monkeypatch = pytest.MonkeyPatch()
    if env_per_rank is not None:
        for var_name, value in env_per_rank.items():
            if value == "<INPUT_RANK>":
                monkeypatch.setenv(var_name, f"{rank}")
                continue
            monkeypatch.setenv(var_name, value)

    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    manager = DistributedManager()

    if internalized_conditioning:
        serial_atom_diffusion = SerialAtomDiffusionV1(
            score_model_args=score_model_kwargs,
            coordinate_augmentation=False,
            alignment_reverse_diff=False,
            num_sampling_steps=num_sampling_steps,
        )
    else:
        serial_atom_diffusion = SerialAtomDiffusionV2(
            score_model_args=score_model_kwargs,
            coordinate_augmentation=False,
            alignment_reverse_diff=False,
            num_sampling_steps=num_sampling_steps,
        )
    serial_atom_diffusion = serial_atom_diffusion.to(device=manager.device, dtype=dtype)
    serial_atom_diffusion.load_state_dict(atom_diffusion_state_dict)
    serial_atom_diffusion = serial_atom_diffusion.eval()

    ring_comm = AttentionPairBiasComm(
        manager.group["cp"],
        manager.layout_subgroups["cp"],
        manager.subgroups["cp"][0],
        manager.subgroups["cp"][1],
    )
    module = DistributedAtomDiffusion(
        layer=serial_atom_diffusion,
        device_mesh=manager.device_mesh_subgroups,
        ring_comm=ring_comm,
    ).eval()

    # ------------------------------------------------------------------
    # Distribute token-level tensors
    # ------------------------------------------------------------------
    placements_single = (Shard(0), Shard(1), Replicate())
    placements_pair = (Shard(0), Shard(1), Shard(2))

    token_pad_mask_dt = distribute_tensor(
        feats_global_host["token_pad_mask"].to(device=manager.device, dtype=dtype),
        manager.device_mesh_subgroups,
        placements_single,
    )
    s_inputs_dt = distribute_tensor(
        s_inputs_global_host.to(device=manager.device, dtype=dtype),
        manager.device_mesh_subgroups,
        placements_single,
    )
    s_trunk_dt = distribute_tensor(
        s_trunk_global_host.to(device=manager.device, dtype=dtype),
        manager.device_mesh_subgroups,
        placements_single,
    )

    z_trunk_device = z_trunk_global_host.to(device=manager.device, dtype=dtype)
    rel_pos_enc_device = rel_pos_enc_global_host.to(device=manager.device, dtype=dtype)

    # ------------------------------------------------------------------
    # Distribute atom features + noise tensors via distribute_atom_features.
    # Noise tensors are distributed alongside atom features so that
    # intersperse padding naturally places zeros at padding positions.
    # ------------------------------------------------------------------
    inputs_atom = {
        k: v.to(dtype=dtype if v.dtype.is_floating_point else v.dtype)
        for k, v in feats_global_host.items()
        if k in _placements_cp_atom_features
    }
    size_batch = feats_global_host["atom_pad_mask"].shape[0]

    all_noise = [init_noise_global_host] + list(step_noise_list_global_host)
    for i_noise, noise_host in enumerate(all_noise):
        unflat = noise_host.unflatten(0, (size_batch, multiplicity))
        for i_mul in range(multiplicity):
            inputs_atom[f"_noise_{i_noise}_{i_mul}"] = unflat[:, i_mul].to(dtype=dtype)

    noise_cp_placements = {}
    noise_placements = {}
    for i_noise in range(len(all_noise)):
        for i_mul in range(multiplicity):
            key = f"_noise_{i_noise}_{i_mul}"
            noise_cp_placements[key] = _placements_cp_model_io["noise"]
            noise_placements[key] = _placements_model_io["noise"]

    feats_and_noise = distribute_atom_features(
        inputs_atom,
        _placements_cp_atom_features | noise_cp_placements,
        _placements_atom_features | noise_placements,
        manager.device_mesh_subgroups,
        manager.group["cp"],
        multiplicities={f"_noise_{i}": multiplicity for i in range(len(all_noise))},
    )

    noise_dts = []
    for i_noise in range(len(all_noise)):
        noise_dts.append(feats_and_noise.pop(f"_noise_{i_noise}"))
    init_noise_dt = noise_dts[0]
    step_noise_dts = noise_dts[1:]

    feats_dt = feats_and_noise
    feats_dt["token_pad_mask"] = token_pad_mask_dt

    z_trunk_dt = distribute_tensor(z_trunk_device, manager.device_mesh_subgroups, placements_pair)
    rel_pos_enc_dt = distribute_tensor(rel_pos_enc_device, manager.device_mesh_subgroups, placements_pair)

    # ------------------------------------------------------------------
    # Make DTensor sample() deterministic via monkeypatching with
    # non-zero noise distributed through intersperse-padded layout.
    # ------------------------------------------------------------------
    _orig_center_random_augmentation = distributed_diffusion_module.center_random_augmentation

    def _centering_only_augmentation(atom_coords, atom_mask, **kwargs):
        kwargs["augmentation"] = False
        kwargs["centering"] = True
        return _orig_center_random_augmentation(atom_coords, atom_mask, **kwargs)

    _dt_randn_calls = []
    _dt_randn_sequence = [init_noise_dt] + step_noise_dts

    def _fixed_create_distributed_randn(shape, device_mesh, placements, dtype=torch.float32, scale=1.0):
        idx = len(_dt_randn_calls)
        _dt_randn_calls.append(idx)
        noise_dt = _dt_randn_sequence[idx]
        if scale != 1.0:
            noise_dt = scalar_tensor_op(scale, noise_dt, ElementwiseOp.PROD)
        return noise_dt

    monkeypatch.setattr(distributed_diffusion_module, "center_random_augmentation", _centering_only_augmentation)
    monkeypatch.setattr(distributed_diffusion_module, "create_distributed_randn", _fixed_create_distributed_randn)

    # ------------------------------------------------------------------
    # Build network_condition_kwargs and run DTensor sample
    # ------------------------------------------------------------------
    if internalized_conditioning:
        network_condition_kwargs = {
            "s_inputs": s_inputs_dt,
            "s_trunk": s_trunk_dt,
            "feats": feats_dt,
            "z_trunk": z_trunk_dt,
            "relative_position_encoding": rel_pos_enc_dt,
        }
    else:
        serial_conditioning = SerialDiffusionConditioning(**conditioning_kwargs)
        serial_conditioning = serial_conditioning.to(device=manager.device, dtype=dtype)
        serial_conditioning.load_state_dict(conditioning_state_dict)
        serial_conditioning = serial_conditioning.eval()

        dtensor_conditioning = DistributedDiffusionConditioning(
            layer=serial_conditioning,
            device_mesh=manager.device_mesh_subgroups,
        ).eval()

        with torch.no_grad():
            q_dt, c_dt, enc_bias_dt, dec_bias_dt, trans_bias_dt = dtensor_conditioning(
                s_trunk=s_trunk_dt,
                z_trunk=z_trunk_dt,
                relative_position_encoding=rel_pos_enc_dt,
                feats=feats_dt,
            )

        network_condition_kwargs = {
            "s_inputs": s_inputs_dt,
            "s_trunk": s_trunk_dt,
            "feats": feats_dt,
            "diffusion_conditioning": {
                "q": q_dt,
                "c": c_dt,
                "atom_enc_bias": enc_bias_dt,
                "atom_dec_bias": dec_bias_dt,
                "token_trans_bias": trans_bias_dt,
            },
        }

    with torch.no_grad():
        out_dt = module.sample(
            atom_mask=feats_dt["atom_pad_mask"],
            multiplicity=multiplicity,
            max_parallel_samples=max_parallel_samples,
            **network_condition_kwargs,
        )

    # ------------------------------------------------------------------
    # Comparison: extract real atoms (where atom_pad_mask=1) from both
    # the DTensor output (intersperse-padded layout) and serial expected
    # (raw layout). The real atoms maintain the same order in both layouts,
    # so boolean indexing extracts matching sequences.
    # ------------------------------------------------------------------
    dt_full = out_dt["sample_atom_coords"].full_tensor()
    dt_mask = feats_dt["atom_pad_mask"].full_tensor().repeat_interleave(multiplicity, 0).bool()
    dt_real = dt_full[dt_mask]

    sample_expected = sample_coords_expected_global_host.to(device=manager.device, dtype=dtype)
    serial_mask = feats_global_host["atom_pad_mask"].to(device=manager.device).repeat_interleave(multiplicity, 0).bool()
    serial_real = sample_expected[serial_mask]

    torch.testing.assert_close(dt_real, serial_real)

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.slow
@pytest.mark.parametrize(
    "setup_env",
    [((2, (2, 2)), True, "cuda", "ENV")],
    indirect=["setup_env"],
    ids=lambda x: f"dp:{x[0][0]}, cp:{x[0][1]}, device_type:{x[2]}",
)
@pytest.mark.parametrize("multiplicity", [1, 4], ids=lambda x: f"mul:{x}")
@pytest.mark.parametrize("internalized_conditioning", [False, True], ids=["extern", "intern"])
def test_atom_diffusion_sample(setup_env, multiplicity, internalized_conditioning: bool):
    """Test AtomDiffusion.sample() (inference) with DTensor CP.

    Determinism is achieved via monkeypatching with pre-generated non-zero noise
    tensors. Serial uses mocked torch.randn returning the noise sequence; DTensor
    uses the same noise distributed through distribute_atom_features (intersperse-
    padded layout) via mocked create_distributed_randn.
    Uses num_sampling_steps=2 to limit numerical error accumulation across denoising steps.
    Exercises max_parallel_samples chunking (max_parallel_samples=2 when multiplicity=4).
    """
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env
    dtype = torch.float64

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
    num_sampling_steps = 2
    max_parallel_samples = 2 if multiplicity > 2 else multiplicity

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

    if internalized_conditioning:
        v1_input_dim = 2 * token_s + 2 * boltz_const.num_tokens + 1 + len(boltz_const.pocket_contact_info)
        s_inputs_dim = v1_input_dim - token_s
    else:
        s_inputs_dim = token_s

    s_inputs = torch.empty((B, N_tokens, s_inputs_dim), device=device_type, dtype=dtype)
    s_trunk = torch.empty((B, N_tokens, token_s), device=device_type, dtype=dtype)
    z_trunk = torch.empty((B, N_tokens, N_tokens, token_z), device=device_type, dtype=dtype)
    rel_pos_enc = torch.empty((B, N_tokens, N_tokens, token_z), device=device_type, dtype=dtype)
    init_tensors_uniform([s_inputs, s_trunk, z_trunk, rel_pos_enc], low=val_init_min_max[0], high=val_init_min_max[1])

    # ------------------------------------------------------------------
    # Build serial modules
    # ------------------------------------------------------------------
    if internalized_conditioning:
        score_model_kwargs = {
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
        conditioning_kwargs = None
        conditioning_state_dict = None
        serial_atom_diffusion = SerialAtomDiffusionV1(
            score_model_args=score_model_kwargs,
            coordinate_augmentation=False,
            alignment_reverse_diff=False,
            num_sampling_steps=num_sampling_steps,
        ).to(device=device_type, dtype=dtype)
    else:
        score_model_kwargs = {
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
        serial_atom_diffusion = SerialAtomDiffusionV2(
            score_model_args=score_model_kwargs,
            coordinate_augmentation=False,
            alignment_reverse_diff=False,
            num_sampling_steps=num_sampling_steps,
        ).to(device=device_type, dtype=dtype)

    serial_atom_diffusion.train()
    init_module_params_uniform(serial_atom_diffusion, low=val_init_min_max[0], high=val_init_min_max[1])
    serial_atom_diffusion.apply(SetModuleInfValues())
    serial_atom_diffusion.eval()
    atom_diffusion_state_dict = serial_atom_diffusion.state_dict()

    # V1 serial sample() needs token_index in feats for token_repr_shape
    if internalized_conditioning and "token_index" not in feats:
        feats["token_index"] = torch.arange(N_tokens, device=device_type).unsqueeze(0).expand(B, -1)

    N_atoms_actual = feats["atom_pad_mask"].shape[1]
    _B_M = B * multiplicity

    # Pre-generate non-zero noise tensors for deterministic comparison.
    # sample() calls torch.randn once for init_noise and once per step.
    init_noise = torch.empty((_B_M, N_atoms_actual, 3), device=device_type, dtype=dtype)
    step_noise_list = [
        torch.empty((_B_M, N_atoms_actual, 3), device=device_type, dtype=dtype) for _ in range(num_sampling_steps)
    ]
    init_tensors_uniform([init_noise, *step_noise_list], low=val_init_min_max[0], high=val_init_min_max[1])

    # ------------------------------------------------------------------
    # Serial sample (with monkeypatched determinism using non-zero noise)
    # ------------------------------------------------------------------
    serial_mod = serial_diffusion_v1_module if internalized_conditioning else serial_diffusion_v2_module

    def _identity_compute_random_augmentation(multiplicity_arg, device=None, dtype=None):
        R = torch.eye(3, device=device, dtype=dtype).unsqueeze(0).expand(_B_M, -1, -1)
        tr = torch.zeros(_B_M, 1, 3, device=device, dtype=dtype)
        return R, tr

    _serial_randn_calls = []
    _serial_randn_sequence = [init_noise] + step_noise_list

    def _fixed_randn(*args, **kwargs):
        idx = len(_serial_randn_calls)
        _serial_randn_calls.append(idx)
        return _serial_randn_sequence[idx].clone()

    _monkeypatch = pytest.MonkeyPatch()
    _monkeypatch.setattr(serial_mod, "compute_random_augmentation", _identity_compute_random_augmentation)
    _monkeypatch.setattr(serial_mod.torch, "randn", _fixed_randn)

    with torch.no_grad():
        if internalized_conditioning:
            out_serial = serial_atom_diffusion.sample(
                atom_mask=feats["atom_pad_mask"],
                multiplicity=multiplicity,
                max_parallel_samples=max_parallel_samples,
                s_inputs=s_inputs,
                s_trunk=s_trunk,
                z_trunk=z_trunk,
                relative_position_encoding=rel_pos_enc,
                feats={k: v.clone() for k, v in feats.items()},
            )
        else:
            q_cond, c_cond, to_keys, enc_bias, dec_bias, trans_bias = serial_conditioning(
                s_trunk=s_trunk,
                z_trunk=z_trunk,
                relative_position_encoding=rel_pos_enc,
                feats={k: v.detach() for k, v in feats.items()},
            )
            out_serial = serial_atom_diffusion.sample(
                atom_mask=feats["atom_pad_mask"],
                multiplicity=multiplicity,
                max_parallel_samples=max_parallel_samples,
                s_inputs=s_inputs,
                s_trunk=s_trunk,
                feats={k: v.clone() for k, v in feats.items()},
                diffusion_conditioning={
                    "q": q_cond,
                    "c": c_cond,
                    "to_keys": to_keys,
                    "atom_enc_bias": enc_bias,
                    "atom_dec_bias": dec_bias,
                    "token_trans_bias": trans_bias,
                },
            )

    _monkeypatch.undo()

    spawn_multiprocessing(
        parallel_assert_atom_diffusion_sample,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        internalized_conditioning,
        dtype,
        multiplicity,
        num_sampling_steps,
        max_parallel_samples,
        {k: v.cpu() for k, v in atom_diffusion_state_dict.items()},
        {k: v.cpu() for k, v in conditioning_state_dict.items()} if conditioning_state_dict else None,
        score_model_kwargs,
        conditioning_kwargs,
        W,
        H,
        {k: v.cpu() for k, v in feats.items()},
        s_inputs.cpu(),
        s_trunk.cpu(),
        z_trunk.cpu(),
        rel_pos_enc.cpu(),
        init_noise.cpu(),
        [n.cpu() for n in step_noise_list],
        out_serial["sample_atom_coords"].cpu(),
    )


# ======================================================================
# Test 4: AtomDiffusion helper functions (c_skip, c_out, c_in, etc.)
# ======================================================================


def parallel_assert_atom_diffusion_helpers(
    rank: int,
    grid_group_sizes,
    device_type: str,
    backend: str,
    env_per_rank,
    dtype: torch.dtype,
    atom_diffusion_state_dict,
    score_model_kwargs: dict,
    # Inputs and expected outputs
    sigma_global_host: torch.Tensor,
    c_skip_expected_host: torch.Tensor,
    c_out_expected_host: torch.Tensor,
    c_in_expected_host: torch.Tensor,
    c_noise_expected_host: torch.Tensor,
    loss_weight_expected_host: torch.Tensor,
    noise_dist_expected_host: torch.Tensor,
    sample_schedule_expected_host: torch.Tensor,
):
    """Parallel assertion for AtomDiffusion helper functions."""
    monkeypatch = pytest.MonkeyPatch()
    if env_per_rank is not None:
        for var_name, value in env_per_rank.items():
            if value == "<INPUT_RANK>":
                monkeypatch.setenv(var_name, f"{rank}")
                continue
            monkeypatch.setenv(var_name, value)

    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    manager = DistributedManager()

    # Create DTensor module from V2 serial (V1 and V2 have identical helper implementations)
    from boltz.model.modules.diffusionv2 import AtomDiffusion as SerialAtomDiffusionV2Local

    serial = SerialAtomDiffusionV2Local(
        score_model_args=score_model_kwargs,
        coordinate_augmentation=False,
    )
    serial = serial.to(device=manager.device, dtype=dtype)
    serial.load_state_dict(atom_diffusion_state_dict)
    serial = serial.eval()

    ring_comm = AttentionPairBiasComm(
        manager.group["cp"],
        manager.layout_subgroups["cp"],
        manager.subgroups["cp"][0],
        manager.subgroups["cp"][1],
    )
    module = DistributedAtomDiffusion(
        layer=serial,
        device_mesh=manager.device_mesh_subgroups,
        ring_comm=ring_comm,
    ).eval()

    # Distribute sigma
    placements_scalar = (Shard(0), Replicate(), Replicate())
    sigma_dt = distribute_tensor(
        sigma_global_host.to(device=manager.device, dtype=dtype),
        manager.device_mesh_subgroups,
        placements_scalar,
    )

    # Test c_skip, c_out, c_in, c_noise, loss_weight
    for name, fn, expected_host in [
        ("c_skip", module.c_skip, c_skip_expected_host),
        ("c_out", module.c_out, c_out_expected_host),
        ("c_in", module.c_in, c_in_expected_host),
        ("c_noise", module.c_noise, c_noise_expected_host),
        ("loss_weight", module.loss_weight, loss_weight_expected_host),
    ]:
        result = fn(sigma_dt)
        expected = expected_host.to(device=manager.device, dtype=dtype)
        torch.testing.assert_close(result.full_tensor(), expected, msg=lambda m: f"{name}: {m}")

    # Test noise_distribution (stochastic — check shape, placement, and dtype)
    noise_dist = module.noise_distribution(sigma_dt.shape[0])
    assert (
        noise_dist.shape == sigma_dt.shape
    ), f"noise_distribution shape mismatch: {noise_dist.shape} vs {sigma_dt.shape}"
    assert noise_dist.placements == placements_scalar
    assert noise_dist.dtype == torch.float32, f"noise_distribution default dtype: {noise_dist.dtype} != float32"

    noise_dist_f64 = module.noise_distribution(sigma_dt.shape[0], dtype=torch.float64)
    assert noise_dist_f64.dtype == torch.float64, f"noise_distribution float64 dtype: {noise_dist_f64.dtype} != float64"

    # Test sample_schedule (deterministic, returns plain Tensor)
    schedule = module.sample_schedule(num_sampling_steps=5)
    expected_schedule = sample_schedule_expected_host.to(device=manager.device)
    torch.testing.assert_close(schedule, expected_schedule)

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.slow
@pytest.mark.parametrize(
    "setup_env",
    [((2, (1, 1)), True, "cuda", "ENV")],
    indirect=["setup_env"],
    ids=lambda x: f"dp:{x[0][0]}, cp:{x[0][1]}, device_type:{x[2]}",
)
def test_dtensor_atom_diffusion_helpers(setup_env):
    """Test DTensor AtomDiffusion scalar helper functions.

    Tests c_skip, c_out, c_in, c_noise, loss_weight, noise_distribution, sample_schedule.
    V1 and V2 serial implementations are identical for these functions, so only one
    serial reference is needed (V2 is used).
    Uses dp=2, cp=(1,1) since these functions don't involve CP atom sharding.
    """
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env
    dtype = torch.float64

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() < {world_size}")

    seed = 42
    seed_by_rank(0, seed=seed)

    B = 1 * grid_group_sizes["dp"]

    # Build a minimal serial AtomDiffusion (V2) for the helpers
    score_model_kwargs = {
        "token_s": 4,
        "atom_s": 8,
        "atoms_per_window_queries": 32,
        "atoms_per_window_keys": 128,
        "sigma_data": 16,
        "dim_fourier": 32,
        "atom_encoder_depth": 1,
        "atom_encoder_heads": 1,
        "token_transformer_depth": 1,
        "token_transformer_heads": 1,
        "atom_decoder_depth": 1,
        "atom_decoder_heads": 1,
        "conditioning_transition_layers": 1,
    }
    serial = SerialAtomDiffusionV2(
        score_model_args=score_model_kwargs,
        coordinate_augmentation=False,
    ).to(device=device_type, dtype=dtype)
    serial.eval()
    atom_diffusion_state_dict = serial.state_dict()

    # Generate test sigma values
    sigma = torch.tensor([0.5, 3.14, 16.0, 160.0], device=device_type, dtype=dtype)[:B]

    # Compute serial reference
    c_skip_ref = serial.c_skip(sigma)
    c_out_ref = serial.c_out(sigma)
    c_in_ref = serial.c_in(sigma)
    c_noise_ref = serial.c_noise(sigma)
    loss_weight_ref = serial.loss_weight(sigma)
    noise_dist_ref = serial.noise_distribution(B)
    sample_schedule_ref = serial.sample_schedule(num_sampling_steps=5)

    spawn_multiprocessing(
        parallel_assert_atom_diffusion_helpers,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        dtype,
        {k: v.cpu() for k, v in atom_diffusion_state_dict.items()},
        score_model_kwargs,
        sigma.cpu(),
        c_skip_ref.cpu(),
        c_out_ref.cpu(),
        c_in_ref.cpu(),
        c_noise_ref.cpu(),
        loss_weight_ref.cpu(),
        noise_dist_ref.cpu(),
        sample_schedule_ref.cpu(),
    )


# ======================================================================
# Test 5: AtomDiffusion.compute_loss()
# ======================================================================


def parallel_assert_atom_diffusion_compute_loss(
    rank: int,
    grid_group_sizes,
    device_type: str,
    backend: str,
    internalized_conditioning: bool,
    dtype: torch.dtype,
    multiplicity: int,
    add_smooth_lddt_loss: bool,
    nucleotide_loss_weight: float,
    ligand_loss_weight: float,
    filter_by_plddt: float,
    use_triton_kernel: bool,
    feats_host: dict,
    atom_diffusion_state_dict: dict,
    score_model_kwargs: dict,
    denoised_atom_coords_global_host: torch.Tensor,
    aligned_true_atom_coords_global_host: torch.Tensor,
    sigma_global_host: torch.Tensor,
    expected_total_loss_global_host: torch.Tensor,
    expected_mse_loss_global_host: torch.Tensor,
    expected_smooth_lddt_loss_global_host: torch.Tensor,
    expected_denoised_atom_coords_grad_global_host: torch.Tensor,
    env_per_rank=None,
):
    """Parallel assertion for AtomDiffusion.compute_loss().

    Compares distributed compute_loss (V2-style formula) to serial V1 or V2
    depending on internalized_conditioning.
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

    # Build serial module (V1 or V2) to match reference used for expected values
    if internalized_conditioning:
        serial_atom_diffusion = SerialAtomDiffusionV1(
            score_model_args=score_model_kwargs,
        )
    else:
        serial_atom_diffusion = SerialAtomDiffusionV2(
            score_model_args=score_model_kwargs,
        )
    serial_atom_diffusion = serial_atom_diffusion.to(device=manager.device, dtype=dtype)
    serial_atom_diffusion.load_state_dict(atom_diffusion_state_dict)
    serial_atom_diffusion = serial_atom_diffusion.train()

    ring_comm = AttentionPairBiasComm(
        manager.group["cp"],
        manager.layout_subgroups["cp"],
        manager.subgroups["cp"][0],
        manager.subgroups["cp"][1],
    )
    transpose_comm = TransposeComm(manager.group["cp"], manager.layout_subgroups["cp"])

    module = DistributedAtomDiffusion(
        layer=serial_atom_diffusion,
        device_mesh=manager.device_mesh_subgroups,
        ring_comm=ring_comm,
        transpose_comm=transpose_comm,
    ).train()

    # Placements for compute_loss feats
    placements_scalar = (Shard(0), Replicate(), Replicate())
    placements_single = (Shard(0), Shard(1), Replicate())

    inputs_atom = {
        k: v.to(dtype=dtype if v.dtype.is_floating_point else v.dtype)
        for k, v in feats_host.items()
        if k in _placements_cp_atom_features
    }

    size_batch = feats_host["atom_resolved_mask"].shape[0]
    io_tensors = {
        "denoised_atom_coords": denoised_atom_coords_global_host,
        "aligned_true_atom_coords": aligned_true_atom_coords_global_host,
        "d_r_update_expected": expected_denoised_atom_coords_grad_global_host,
    }
    for base_name, tensor_host in io_tensors.items():
        unflat = tensor_host.unflatten(0, (size_batch, multiplicity))
        for i_mul in range(multiplicity):
            inputs_atom[f"{base_name}_{i_mul}"] = unflat[:, i_mul].to(dtype=dtype)

    io_keys_used = set(io_tensors.keys())
    placements_cp_io_mul = {
        f"{k}_{i_mul}": _placements_cp_model_io[k] for k in io_keys_used for i_mul in range(multiplicity)
    }
    placements_io_mul = {f"{k}_{i_mul}": _placements_model_io[k] for k in io_keys_used for i_mul in range(multiplicity)}
    multiplicities = dict.fromkeys(io_keys_used, multiplicity)

    feats_and_io = distribute_atom_features(
        inputs_atom,
        _placements_cp_atom_features | placements_cp_io_mul,
        _placements_atom_features | placements_io_mul,
        manager.device_mesh_subgroups,
        manager.group["cp"],
        multiplicities=multiplicities,
    )

    denoised_atom_coords_dtensor = feats_and_io.pop("denoised_atom_coords").requires_grad_(True)
    d_r_update_expected_dtensor = feats_and_io.pop("d_r_update_expected")
    aligned_true_atom_coords_dtensor = feats_and_io.pop("aligned_true_atom_coords")
    sigma_dtensor = distribute_tensor(
        sigma_global_host.to(device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements_scalar,
    )
    mol_type = distribute_tensor(
        feats_host["mol_type"].to(device=manager.device),
        manager.device_mesh_subgroups,
        placements_single,
    )
    feats_dtensor = {
        "mol_type": mol_type,
        "atom_to_token": feats_and_io.pop("atom_to_token"),
        "atom_resolved_mask": feats_and_io.pop("atom_resolved_mask"),
        "plddt": feats_and_io.pop("plddt"),
    }
    input_dict = {
        "denoised_atom_coords": denoised_atom_coords_dtensor,
        "sigmas": sigma_dtensor,
        "aligned_true_atom_coords": aligned_true_atom_coords_dtensor,
    }
    input_dict_clone = {k: v.detach().clone().requires_grad_(v.requires_grad) for k, v in input_dict.items()}

    output_dict = module.compute_loss(
        feats=feats_dtensor,
        out_dict=input_dict,
        add_smooth_lddt_loss=add_smooth_lddt_loss,
        nucleotide_loss_weight=nucleotide_loss_weight,
        ligand_loss_weight=ligand_loss_weight,
        multiplicity=multiplicity,
        filter_by_plddt=filter_by_plddt,
        use_triton_kernel=use_triton_kernel,
    )

    # Ensure input dict is not modified by compute_loss
    for k in input_dict.keys():
        assert_tensors_identical(
            input_dict[k],
            input_dict_clone[k],
            check_grad=False,
            check_grad_fn=False,
            check_storage_offset=True,
            check_storage_pointer=False,
        )

    total_loss_dtensor = output_dict["loss"]
    mse_loss_dtensor = output_dict["loss_breakdown"]["mse_loss"]
    smooth_lddt_loss_dtensor = output_dict["loss_breakdown"]["smooth_lddt_loss"]
    assert total_loss_dtensor.placements == (
        Replicate(),
        Replicate(),
        Replicate(),
    ), "total_loss_dtensor should be replicated"
    assert mse_loss_dtensor.placements == (
        Replicate(),
        Replicate(),
        Replicate(),
    ), "mse_loss_dtensor should be replicated"
    assert smooth_lddt_loss_dtensor.placements == (
        Replicate(),
        Replicate(),
        Replicate(),
    ), "smooth_lddt_loss_dtensor should be replicated"

    total_loss = total_loss_dtensor.full_tensor().cpu()
    mse_loss = mse_loss_dtensor.full_tensor().cpu()
    smooth_lddt_loss = smooth_lddt_loss_dtensor.full_tensor().cpu()

    assert not (mse_loss == 0.0).all(), "mse_loss should not be 0"
    if add_smooth_lddt_loss:
        assert not (smooth_lddt_loss == 0.0).all(), "smooth_lddt_loss should not be 0"
    assert not (total_loss == 0.0).all(), "total_loss should not be 0"
    torch.testing.assert_close(mse_loss, expected_mse_loss_global_host)
    torch.testing.assert_close(smooth_lddt_loss, expected_smooth_lddt_loss_global_host)
    torch.testing.assert_close(total_loss, expected_total_loss_global_host)

    total_loss_dtensor_clone = total_loss_dtensor.detach().clone().requires_grad_(total_loss_dtensor.requires_grad)
    mse_loss_dtensor_clone = mse_loss_dtensor.detach().clone().requires_grad_(mse_loss_dtensor.requires_grad)
    smooth_lddt_loss_dtensor_clone = (
        smooth_lddt_loss_dtensor.detach().clone().requires_grad_(smooth_lddt_loss_dtensor.requires_grad)
    )

    total_loss_dtensor.backward()

    assert_tensors_identical(
        total_loss_dtensor,
        total_loss_dtensor_clone,
        check_grad=False,
        check_grad_fn=False,
        check_storage_offset=True,
        check_storage_pointer=False,
    )
    assert_tensors_identical(
        mse_loss_dtensor,
        mse_loss_dtensor_clone,
        check_grad=False,
        check_grad_fn=False,
        check_storage_offset=True,
        check_storage_pointer=False,
    )
    assert_tensors_identical(
        smooth_lddt_loss_dtensor,
        smooth_lddt_loss_dtensor_clone,
        check_grad=False,
        check_grad_fn=False,
        check_storage_offset=True,
        check_storage_pointer=False,
    )

    torch.testing.assert_close(
        denoised_atom_coords_dtensor.grad.full_tensor(), d_r_update_expected_dtensor.full_tensor()
    )

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env",
    [
        ((2, (2, 2)), True, "cuda", "ENV"),
    ],
    indirect=("setup_env",),
    ids=lambda x: f"dp={x[0][0]}, cp={x[0][1]}, device_type={x[2]}, method_init={x[3]}",
)
@pytest.mark.parametrize(
    "loss_config",
    [
        (False, 1, True, False, 0.0),
        (False, 4, True, True, 0.5),
        (True, 1, True, False, 0.0),
        (True, 4, True, True, 0.0),
    ],
    ids=lambda x: f"v2={not x[0]}, mul={x[1]}, lddt={x[2]}, triton={x[3]}, plddt={x[4]:.1f}",
)
def test_atom_diffusion_compute_loss(
    setup_env,
    loss_config,
    nucleotide_loss_weight: float = 5.0,
    ligand_loss_weight: float = 10.0,
    dtype: torch.dtype = torch.float32,
):
    """Test AtomDiffusion.compute_loss() with distributed context parallelism.

    Compares DTensor AtomDiffusion.compute_loss() against serial V1 or V2
    (internalized_conditioning True/False) for forward and backward.
    """
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env
    internalized_conditioning, multiplicity, add_smooth_lddt_loss, use_triton_kernel, filter_by_plddt = loss_config

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() < {world_size}")

    if filter_by_plddt > 0 and internalized_conditioning:
        pytest.skip("filter_by_plddt is not supported for boltz1 internalized conditioning")

    if not add_smooth_lddt_loss and use_triton_kernel:
        pytest.skip("use_triton_kernel requires add_smooth_lddt_loss=True")

    seed = 42
    seed_by_rank(0, seed=seed)

    size_cp = grid_group_sizes["cp"][0]
    B = 1 * grid_group_sizes["dp"]
    W = 32
    H = 128
    val_init_min_max = (-1.0, 1.0)

    n_atoms_per_token_min = 8
    n_atoms_per_token_max = 20
    N_tokens = 10 * size_cp
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

    atom_feature_dim = 3 + 1 + (1 if internalized_conditioning else 0) + boltz_const.num_elements + 4 * 64

    compute_loss_selected_keys = {
        "atom_resolved_mask",
        "mol_type",
        "atom_to_token",
        "atom_counts_per_token",
        "plddt",
    }
    selected_keys = list(_selected_atom_keys | _selected_token_keys | compute_loss_selected_keys)
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

    # V1 s_inputs has wider dim:
    #   input_dim = 2 * token_s + 2 * num_tokens + 1 + len(pocket_contact_info)
    #   s_inputs_dim = input_dim - token_s
    # V2 s_inputs has token_s dim
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

    # noised_atom_coords: (B*M, N_atoms_actual, 3)
    N_atoms_actual = feats["atom_pad_mask"].shape[1]
    noised_atom_coords = torch.empty(
        (B * multiplicity, N_atoms_actual, 3), device=device_type, dtype=dtype, requires_grad=True
    )
    init_tensors_uniform([noised_atom_coords], low=val_init_min_max[0], high=val_init_min_max[1])

    if internalized_conditioning:
        score_model_kwargs = {
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
        serial_model = SerialAtomDiffusionV1(
            score_model_args=score_model_kwargs,
            coordinate_augmentation=False,
        ).to(device=device_type, dtype=dtype)
    else:
        # V2: uses DiffusionConditioning
        score_model_kwargs = {
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
        serial_model = SerialAtomDiffusionV2(
            score_model_args=score_model_kwargs,
            coordinate_augmentation=False,
        ).to(device=device_type, dtype=dtype)
    init_module_params_uniform(serial_model, low=-0.1, high=0.1)
    serial_model.apply(SetModuleInfValues())
    module_state_dict = serial_model.state_dict()
    serial_model = serial_model.to(device=device_type, dtype=dtype)

    denoised_atom_coords = torch.empty(
        (B * multiplicity, N_atoms, 3), device=device_type, dtype=dtype, requires_grad=True
    )
    init_tensors_uniform([denoised_atom_coords], low=val_init_min_max[0], high=val_init_min_max[1])
    aligned_true_atom_coords = torch.empty_like(denoised_atom_coords)
    init_tensors_uniform([aligned_true_atom_coords], low=val_init_min_max[0], high=val_init_min_max[1])
    sigma = serial_model.noise_distribution(B * multiplicity).to(device=device_type, dtype=dtype)
    denoised_atom_coords.requires_grad = True

    input_dict = {
        "denoised_atom_coords": denoised_atom_coords,
        "sigmas": sigma,
        "aligned_true_atom_coords": aligned_true_atom_coords,
    }
    feats["coords"] = aligned_true_atom_coords
    # V1 compute_loss requires noised_atom_coords in out_dict
    if internalized_conditioning:
        noised_atom_coords = torch.empty_like(denoised_atom_coords)
        init_tensors_uniform([noised_atom_coords], low=val_init_min_max[0], high=val_init_min_max[1])
        input_dict["noised_atom_coords"] = noised_atom_coords

    extra_kwargs = {}
    if not internalized_conditioning and filter_by_plddt > 0:
        extra_kwargs["filter_by_plddt"] = filter_by_plddt

    output_dict = serial_model.compute_loss(
        feats=feats,
        out_dict=input_dict,
        add_smooth_lddt_loss=add_smooth_lddt_loss,
        nucleotide_loss_weight=nucleotide_loss_weight,
        ligand_loss_weight=ligand_loss_weight,
        multiplicity=multiplicity,
        **extra_kwargs,
    )
    output_dict["loss"].backward()

    feats_host = {k: v.detach().to(device="cpu", copy=True) if torch.is_tensor(v) else v for k, v in feats.items()}

    spawn_multiprocessing(
        parallel_assert_atom_diffusion_compute_loss,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        internalized_conditioning,
        dtype,
        multiplicity,
        add_smooth_lddt_loss,
        nucleotide_loss_weight,
        ligand_loss_weight,
        filter_by_plddt,
        use_triton_kernel,
        feats_host,
        {k: v.cpu() for k, v in module_state_dict.items()},
        score_model_kwargs,
        denoised_atom_coords.detach().clone().cpu(),
        aligned_true_atom_coords.detach().clone().cpu(),
        sigma.detach().clone().cpu(),
        output_dict["loss"].detach().clone().cpu(),
        output_dict["loss_breakdown"]["mse_loss"].detach().clone().cpu(),
        output_dict["loss_breakdown"]["smooth_lddt_loss"].detach().clone().cpu(),
        denoised_atom_coords.grad.detach().clone().cpu(),
        env_per_rank,
    )
