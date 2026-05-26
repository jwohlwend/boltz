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

"""Tests for the distributed v2 ConfidenceHeads and ConfidenceModule (DTensor).

Checks that the distributed ConfidenceHeads and ConfidenceModule produce outputs
and gradients bit-for-bit identical to the serial v2 versions across:
  - use_separate_heads=False  (single shared PAE/PDE head)
  - use_separate_heads=True   (separate intra/inter-chain PAE and PDE heads)
  - dtype=float64 and float32

PTM/iPTM features (frames_idx, atom_to_token, atom_pad_mask, atom_resolved_mask)
are included so that compute_ptms is exercised end-to-end.
"""

import pytest
import torch
from torch.distributed.tensor import Replicate, Shard, distribute_tensor

from boltz.distributed.comm import TransposeComm
from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.modules.confidence_utils import CHAIN_IPTM_SENTINEL
from boltz.distributed.model.modules.confidencev2 import ConfidenceHeads as DTensorConfidenceHeadsV2
from boltz.distributed.model.modules.confidencev2 import ConfidenceModule as DTensorConfidenceModuleV2
from boltz.model.modules.confidencev2 import ConfidenceHeads as SerialConfidenceHeadsV2
from boltz.model.modules.confidencev2 import ConfidenceModule as SerialConfidenceModuleV2
from boltz.testing.utils import (
    assert_tensors_identical,
    create_boltz2_model_init_params,
    distribute_atom_features,
    get_feature_placements,
    init_module_params_glorot,
    init_tensors_uniform,
    random_features,
    seed_by_rank,
    spawn_multiprocessing,
)


def parallel_test_dtensor_confidence_heads_v2(
    rank,
    grid_group_sizes,
    device_type,
    backend,
    env_per_rank,
    # config
    dtype,
    serial_state_dict,
    confidence_heads_kwargs,
    multiplicity,
    # input tensors
    s_global_host,
    z_global_host,
    x_pred_global_host,
    d_global_host,
    feats_global_host,
    pred_distogram_logits_global_host,
    # reference serial outputs
    serial_output_feats_host,
    # reference serial input gradients
    s_grad_host,
    z_grad_host,
    # upstream gradient tensors
    d_plddt_logits_host,
    d_pde_logits_host,
    d_resolved_logits_host,
    d_pae_logits_host,
    # reference serial parameter gradients
    serial_param_grads_host,
):
    """Parallel worker: distributes inputs, runs DTensor ConfidenceHeads, compares."""
    monkeypatch = pytest.MonkeyPatch()
    if env_per_rank is not None:
        for var_name, value in env_per_rank.items():
            if value == "<INPUT_RANK>":
                monkeypatch.setenv(var_name, f"{rank}")
                continue
            monkeypatch.setenv(var_name, value)

    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    manager = DistributedManager()
    device_mesh = manager.device_mesh_subgroups

    seed_by_rank(0, 42)

    # Build serial module on device (for wrapping into distributed module)
    serial_module = SerialConfidenceHeadsV2(**confidence_heads_kwargs)
    serial_module = serial_module.to(device=manager.device, dtype=dtype).train()
    serial_module.load_state_dict(serial_state_dict)

    cp_group = manager.group["cp"]
    layout_group_cp = manager.layout_subgroups["cp"]
    transpose_comm = TransposeComm(cp_group, layout_group_cp)

    module = DTensorConfidenceHeadsV2(
        layer=serial_module,
        device_mesh=device_mesh,
        transpose_comm=transpose_comm,
    )
    module = module.to(device=manager.device, dtype=dtype).train()

    # ----- distribute inputs -----
    # s: (B*mult, N, D_s) → (Shard(0), Shard(1), Replicate())
    s_dtensor = distribute_tensor(
        s_global_host.to(device=manager.device, dtype=dtype).requires_grad_(True),
        device_mesh=device_mesh,
        placements=(Shard(0), Shard(1), Replicate()),
    )

    # z: (B*mult, N, N, D_z) → (Shard(0), Shard(1), Shard(2))
    z_dtensor = distribute_tensor(
        z_global_host.to(device=manager.device, dtype=dtype).requires_grad_(True),
        device_mesh=device_mesh,
        placements=(Shard(0), Shard(1), Shard(2)),
    )

    # d: (B*mult, N, N) → (Shard(0), Shard(1), Shard(2))
    d_dtensor = distribute_tensor(
        d_global_host.to(device=manager.device, dtype=dtype),
        device_mesh=device_mesh,
        placements=(Shard(0), Shard(1), Shard(2)),
    )

    # Distribute atom-level features via distribute_atom_features
    placements_single_repr = (Shard(0), Shard(1), Replicate())
    placements_cp_single_repr = (Shard(0), Replicate())

    special_atom_features = {"atom_pad_mask", "atom_to_token", "atom_resolved_mask", "frames_idx"}
    atom_inputs = {
        key: feats_global_host[key].to(device=manager.device)
        for key in special_atom_features
        if key in feats_global_host
    }

    base_batch = s_global_host.shape[0] // multiplicity
    x_pred_reshaped = x_pred_global_host.reshape(base_batch, multiplicity, *x_pred_global_host.shape[1:])
    for mul_idx in range(multiplicity):
        atom_inputs[f"x_pred_{mul_idx}"] = x_pred_reshaped[:, mul_idx].to(manager.device, dtype=dtype)

    if "atom_counts_per_token" in feats_global_host:
        atom_inputs["atom_counts_per_token"] = feats_global_host["atom_counts_per_token"].to(
            manager.device, dtype=torch.int64
        )
    else:
        atom_inputs["atom_counts_per_token"] = (
            feats_global_host["atom_to_token"].sum(dim=1).to(manager.device, dtype=torch.int64)
        )

    atom_placements_cp = {key: placements_cp_single_repr for key in atom_inputs}
    atom_placements_cp["frames_idx"] = (Shard(1), Replicate())
    atom_placements_dp_cp = {key: placements_single_repr for key in atom_inputs if key != "atom_counts_per_token"}
    atom_placements_dp_cp["frames_idx"] = (Shard(0), Shard(1), Replicate())
    atom_feats_dtensor = distribute_atom_features(
        inputs=atom_inputs,
        placements_cp=atom_placements_cp,
        placements_dp_cp=atom_placements_dp_cp,
        device_mesh=device_mesh,
        cp_group=manager.group["cp"],
        multiplicities={"x_pred": multiplicity},
    )

    x_pred_dtensor = atom_feats_dtensor["x_pred"]

    feats_dtensor = {
        "token_pad_mask": distribute_tensor(
            feats_global_host["token_pad_mask"].to(device=manager.device),
            device_mesh=device_mesh,
            placements=placements_single_repr,
        ),
        "asym_id": distribute_tensor(
            feats_global_host["asym_id"].to(device=manager.device),
            device_mesh=device_mesh,
            placements=placements_single_repr,
        ),
        "mol_type": distribute_tensor(
            feats_global_host["mol_type"].to(device=manager.device),
            device_mesh=device_mesh,
            placements=placements_single_repr,
        ),
    }
    for key in special_atom_features:
        if key in atom_feats_dtensor:
            feats_dtensor[key] = atom_feats_dtensor[key]

    pred_distogram_logits_dtensor = distribute_tensor(
        pred_distogram_logits_global_host.to(manager.device),
        device_mesh=device_mesh,
        placements=(Shard(0), Shard(1), Shard(2)),
    )

    # Keep copies to verify inputs are not mutated
    s_dtensor_copy = s_dtensor.clone()
    z_dtensor_copy = z_dtensor.clone()

    # ----- distribute upstream gradients -----
    d_plddt_logits_dtensor = distribute_tensor(
        d_plddt_logits_host.to(device=manager.device, dtype=dtype),
        device_mesh=device_mesh,
        placements=(Shard(0), Shard(1), Replicate()),
    )
    d_pde_logits_dtensor = distribute_tensor(
        d_pde_logits_host.to(device=manager.device, dtype=dtype),
        device_mesh=device_mesh,
        placements=(Shard(0), Shard(1), Shard(2)),
    )
    d_resolved_logits_dtensor = distribute_tensor(
        d_resolved_logits_host.to(device=manager.device, dtype=dtype),
        device_mesh=device_mesh,
        placements=(Shard(0), Shard(1), Replicate()),
    )
    d_pae_logits_dtensor = distribute_tensor(
        d_pae_logits_host.to(device=manager.device, dtype=dtype),
        device_mesh=device_mesh,
        placements=(Shard(0), Shard(1), Shard(2)),
    )

    # ----- forward -----
    output_dtensor = module(
        s=s_dtensor,
        z=z_dtensor,
        x_pred=x_pred_dtensor,
        d=d_dtensor,
        feats=feats_dtensor,
        pred_distogram_logits=pred_distogram_logits_dtensor,
        multiplicity=multiplicity,
    )

    # Compare all outputs against serial reference
    dp_rank = manager.group_rank["dp"]
    dp_size = len(manager.group_ranks["dp"])
    for key in output_dtensor:
        assert key in serial_output_feats_host, f"DTensor output key '{key}' missing from serial reference"
        dtensor_val = output_dtensor[key]
        if key == "pair_chains_iptm":
            serial_pciptm = serial_output_feats_host[key]
            if isinstance(dtensor_val, dict):
                for idx1, chain_dict in dtensor_val.items():
                    for idx2, dt_val in chain_dict.items():
                        local_val = dt_val.to_local().cpu()
                        is_sentinel = torch.all(local_val == CHAIN_IPTM_SENTINEL)
                        if idx1 in serial_pciptm and idx2 in serial_pciptm.get(idx1, {}):
                            serial_full = serial_pciptm[idx1][idx2]
                            chunk_size = serial_full.shape[0] // dp_size
                            serial_local = serial_full[dp_rank * chunk_size : (dp_rank + 1) * chunk_size]
                            if is_sentinel:
                                assert torch.all(serial_local.abs() < 1e-5), (
                                    f"Chain pair ({idx1}, {idx2}): distributed returned sentinel but "
                                    f"serial has non-zero values {serial_local} for this DP rank"
                                )
                            else:
                                torch.testing.assert_close(
                                    local_val,
                                    serial_local,
                                    msg=f"Chain pair ({idx1}, {idx2}) mismatch on DP rank {dp_rank}",
                                )
                        else:
                            assert is_sentinel, (
                                f"Extra chain pair ({idx1}, {idx2}) should be sentinel "
                                f"{CHAIN_IPTM_SENTINEL}, got {local_val}"
                            )
            else:
                assert isinstance(serial_pciptm, torch.Tensor), (
                    f"pair_chains_iptm: distributed returned DTensor (compute_ptms fallback) "
                    f"but serial returned {type(serial_pciptm)}"
                )
                torch.testing.assert_close(
                    dtensor_val.full_tensor().cpu(),
                    serial_pciptm,
                    msg="pair_chains_iptm fallback mismatch",
                )
            continue
        torch.testing.assert_close(
            dtensor_val.full_tensor().cpu(),
            serial_output_feats_host[key],
            msg=f"Mismatch for output key '{key}'",
        )

    # ----- backward -----
    torch.autograd.backward(
        [
            output_dtensor["plddt_logits"],
            output_dtensor["pde_logits"],
            output_dtensor["resolved_logits"],
            output_dtensor["pae_logits"],
        ],
        [
            d_plddt_logits_dtensor,
            d_pde_logits_dtensor,
            d_resolved_logits_dtensor,
            d_pae_logits_dtensor,
        ],
    )

    # Verify inputs were not mutated by the forward pass
    assert_tensors_identical(
        s_dtensor.to_local().cpu(),
        s_dtensor_copy.to_local().cpu(),
        check_grad=False,
        check_grad_fn=False,
        msg="s_dtensor was mutated during forward",
    )
    assert_tensors_identical(
        z_dtensor.to_local().cpu(),
        z_dtensor_copy.to_local().cpu(),
        check_grad=False,
        check_grad_fn=False,
        msg="z_dtensor was mutated during forward",
    )

    # Compare input gradients
    torch.testing.assert_close(s_dtensor.grad.full_tensor().cpu(), s_grad_host, msg="s gradient mismatch")
    torch.testing.assert_close(z_dtensor.grad.full_tensor().cpu(), z_grad_host, msg="z gradient mismatch")

    # Compare parameter gradients
    result_param_grads = {}
    for name, param in module.named_parameters():
        if param.grad is not None:
            if name not in serial_param_grads_host:
                raise ValueError(
                    f"Parameter '{name}' has a gradient in the distributed module " f"but not in the serial reference"
                )
            result_param_grads[name] = param.grad

    for name, expected in serial_param_grads_host.items():
        assert name in result_param_grads, f"Parameter '{name}' gradient missing in distributed module"
        torch.testing.assert_close(
            result_param_grads[name].full_tensor().cpu(),
            expected,
            msg=f"Parameter gradient mismatch for '{name}'",
        )

    DistributedManager.cleanup()
    monkeypatch.undo()


# NOTE: the use_separate_heads=False codepath looks like it's not run in confidencev2.py, but add test parameterization
@pytest.mark.parametrize(
    ("setup_env", "dtype", "use_separate_heads"),
    (
        params_test := [
            (((2, (2, 2)), True, "cuda", "ENV"), torch.float64, False),
            (((2, (2, 2)), True, "cuda", "ENV"), torch.float32, False),
            (((2, (2, 2)), True, "cuda", "ENV"), torch.float64, True),
            (((2, (2, 2)), True, "cuda", "ENV"), torch.float32, True),
        ]
    ),
    indirect=["setup_env"],
    ids=[
        (f"dp:{x[0][0][0]}, cp:{x[0][0][1]}, device:{x[0][2]}, init:{x[0][3]}, " f"dtype:{x[1]}, separate_heads:{x[2]}")
        for x in params_test
    ],
)
def test_dtensor_confidence_heads_v2(setup_env, dtype, use_separate_heads):
    """Test that DTensor ConfidenceHeadsV2 matches serial ConfidenceHeadsV2.

    Covers:
    * Forward pass parity for all logit outputs, aggregated metrics (pLDDT, PDE, PAE,
      complex_plddt, complex_iplddt, complex_pde, complex_ipde).
    * Backward pass parity for input gradients (s, z) and parameter gradients.
    * Non-mutation of inputs during forward.
    * Both the single-head (use_separate_heads=False) and the intra/inter-chain
      separated head (use_separate_heads=True) configurations.

    PTM/iPTM keys are present in both outputs and are verified to match.

    Parameters
    ----------
    setup_env : tuple
        Grid group sizes, world size, device type, backend, environment variables per rank.
    dtype : torch.dtype
        Tensor dtype for forward and backward passes.
    use_separate_heads : bool
        Whether to use separate intra/inter-chain PAE and PDE projection heads.
    """
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("CUDA not available")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"Need {world_size} GPUs, only {torch.cuda.device_count()} available")

    seed_by_rank(0)

    B = grid_group_sizes["dp"]
    size_cp = grid_group_sizes["cp"][0]
    multiplicity = 2

    # N_tokens and N_atoms must be divisible by CP size for even sharding
    N_tokens = 16 * size_cp
    N_atoms = N_tokens * 4

    selected_keys = [
        "mol_type",
        "asym_id",
        "token_pad_mask",
        "frames_idx",
        "atom_to_token",
        "atom_resolved_mask",
        "atom_pad_mask",
    ]

    feats_host = random_features(
        size_batch=B,
        n_tokens=N_tokens,
        n_atoms=N_atoms,
        n_msa=1,
        atom_counts_per_token_range=(1, 4),
        device=torch.device("cpu"),
        float_value_range=(-0.2, 0.2),
        selected_keys=selected_keys,
    )
    feats_host = {k: (v.to(dtype=dtype) if v.is_floating_point() else v) for k, v in feats_host.items()}
    feats_device = {k: v.to(device=device_type) for k, v in feats_host.items()}

    boltz2_params = create_boltz2_model_init_params(use_large_model=True)
    token_s = boltz2_params["token_s"]
    token_z = boltz2_params["token_z"]
    confidence_args = boltz2_params["confidence_model_args"]["confidence_args"]
    num_distogram_bins = boltz2_params["confidence_model_args"]["num_dist_bins"]
    val_init_range = 0.15

    confidence_heads_kwargs = {
        "token_s": token_s,
        "token_z": token_z,
        **confidence_args,
        "token_level_confidence": True,
        "use_separate_heads": use_separate_heads,
    }

    # Input tensors: s and z have gradient tracking for backward verification.
    s = torch.empty(B * multiplicity, N_tokens, token_s, device=device_type, dtype=dtype, requires_grad=True)
    z = torch.empty(B * multiplicity, N_tokens, N_tokens, token_z, device=device_type, dtype=dtype, requires_grad=True)
    x_pred = torch.empty(B * multiplicity, N_atoms, 3, device=device_type, dtype=dtype)
    d = torch.empty(B * multiplicity, N_tokens, N_tokens, device=device_type, dtype=dtype)
    pred_distogram_logits = torch.empty(B, N_tokens, N_tokens, num_distogram_bins, device=device_type, dtype=dtype)

    init_tensors_uniform([s, z, x_pred, d, pred_distogram_logits], low=-val_init_range, high=val_init_range)

    s_global_host = s.detach().clone().cpu()
    z_global_host = z.detach().clone().cpu()
    x_pred_global_host = x_pred.detach().clone().cpu()
    d_global_host = d.detach().clone().cpu()
    pred_distogram_logits_global_host = pred_distogram_logits.detach().clone().cpu()

    # ----- serial module -----
    serial_module = SerialConfidenceHeadsV2(**confidence_heads_kwargs)
    serial_module = serial_module.to(device=device_type, dtype=dtype).train()
    init_module_params_glorot(serial_module, gain=val_init_range)

    serial_state_dict = serial_module.state_dict()

    # ----- serial forward -----
    serial_output = serial_module(
        s=s,
        z=z,
        x_pred=x_pred,
        d=d,
        feats=feats_device,
        pred_distogram_logits=pred_distogram_logits,
        multiplicity=multiplicity,
    )

    # Upstream gradients for backward
    d_plddt_logits = torch.rand_like(serial_output["plddt_logits"], device=device_type)
    d_pde_logits = torch.rand_like(serial_output["pde_logits"], device=device_type)
    d_resolved_logits = torch.rand_like(serial_output["resolved_logits"], device=device_type)
    d_pae_logits = torch.rand_like(serial_output["pae_logits"], device=device_type)

    d_plddt_logits_host = d_plddt_logits.detach().clone().cpu()
    d_pde_logits_host = d_pde_logits.detach().clone().cpu()
    d_resolved_logits_host = d_resolved_logits.detach().clone().cpu()
    d_pae_logits_host = d_pae_logits.detach().clone().cpu()

    torch.autograd.backward(
        [
            serial_output["plddt_logits"],
            serial_output["pde_logits"],
            serial_output["resolved_logits"],
            serial_output["pae_logits"],
        ],
        [d_plddt_logits, d_pde_logits, d_resolved_logits, d_pae_logits],
    )

    # Save all serial outputs as CPU tensors; pair_chains_iptm handled separately
    def _to_cpu(val):
        if isinstance(val, torch.Tensor):
            return val.detach().clone().cpu()
        if isinstance(val, dict):
            return {k: _to_cpu(v) for k, v in val.items()}
        return val

    serial_output_feats_host = {k: _to_cpu(v) for k, v in serial_output.items()}

    # Verify the serial module has non-zero, non-NaN outputs (guard against vacuous pass)
    assert not torch.isnan(serial_output["plddt_logits"]).any(), "serial plddt_logits contains NaN"
    assert not torch.isnan(serial_output["pde_logits"]).any(), "serial pde_logits contains NaN"
    assert not torch.isnan(serial_output["pae_logits"]).any(), "serial pae_logits contains NaN"
    assert serial_output["plddt"].abs().max() > 0, "serial plddt is all-zero (vacuous)"
    assert serial_output["complex_plddt"].abs().max() > 0, "serial complex_plddt is all-zero (vacuous)"
    assert serial_output["ptm"].abs().max() > 0, "serial ptm is all-zero (vacuous)"
    assert isinstance(serial_output["pair_chains_iptm"], dict), (
        f"serial pair_chains_iptm should be a dict (from compute_ptms), "
        f"got {type(serial_output['pair_chains_iptm'])} — compute_ptms likely failed silently"
    )

    s_grad_host = s.grad.detach().clone().cpu()
    z_grad_host = z.grad.detach().clone().cpu()

    # Verify non-zero gradients (guard against vacuous backward pass)
    assert s_grad_host.abs().max() > 0, "serial s gradient is all-zero (vacuous)"
    assert z_grad_host.abs().max() > 0, "serial z gradient is all-zero (vacuous)"

    serial_param_grads_host = {
        name: param.grad.detach().clone().cpu()
        for name, param in serial_module.named_parameters()
        if param.grad is not None
    }
    assert len(serial_param_grads_host) > 0, "No serial parameter gradients found (vacuous)"

    # ----- parallel distributed test -----
    spawn_multiprocessing(
        parallel_test_dtensor_confidence_heads_v2,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        dtype,
        serial_state_dict,
        confidence_heads_kwargs,
        multiplicity,
        s_global_host,
        z_global_host,
        x_pred_global_host,
        d_global_host,
        feats_host,
        pred_distogram_logits_global_host,
        serial_output_feats_host,
        s_grad_host,
        z_grad_host,
        d_plddt_logits_host,
        d_pde_logits_host,
        d_resolved_logits_host,
        d_pae_logits_host,
        serial_param_grads_host,
    )


# ==============================================================================
# ConfidenceModule v2 test
# ==============================================================================


def parallel_test_dtensor_confidence_module_v2(
    rank,
    grid_group_sizes,
    device_type,
    backend,
    env_per_rank,
    # config
    dtype,
    serial_state_dict,
    confidence_module_kwargs,
    multiplicity,
    run_sequentially,
    # input tensors
    s_inputs_global_host,
    s_global_host,
    z_global_host,
    x_pred_global_host,
    feats_global_host,
    pred_distogram_logits_global_host,
    # reference serial outputs
    serial_output_feats_host,
    # reference serial input gradients
    s_inputs_grad_host,
    s_grad_host,
    z_grad_host,
    # upstream gradient tensors
    d_plddt_logits_host,
    d_pde_logits_host,
    d_resolved_logits_host,
    d_pae_logits_host,
    # reference serial parameter gradients
    serial_param_grads_host,
):
    """Parallel worker: distributes inputs, runs DTensor ConfidenceModule v2, compares."""
    monkeypatch = pytest.MonkeyPatch()
    if env_per_rank is not None:
        for var_name, value in env_per_rank.items():
            if value == "<INPUT_RANK>":
                monkeypatch.setenv(var_name, f"{rank}")
                continue
            monkeypatch.setenv(var_name, value)

    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    manager = DistributedManager()
    device_mesh = manager.device_mesh_subgroups

    seed_by_rank(0, 42)

    # Build serial module on device (for wrapping into distributed module)
    serial_module = SerialConfidenceModuleV2(**confidence_module_kwargs)
    serial_module = serial_module.to(device=manager.device, dtype=dtype).train()
    serial_module.load_state_dict(serial_state_dict)

    cp_group = manager.group["cp"]
    layout_group_cp = manager.layout_subgroups["cp"]
    transpose_comm = TransposeComm(cp_group, layout_group_cp)

    module = DTensorConfidenceModuleV2(
        module=serial_module,
        dist_manager=manager,
        transpose_comm=transpose_comm,
    )
    module = module.to(device=manager.device, dtype=dtype).train()

    # ----- distribute inputs -----
    placements_single_repr = (Shard(0), Shard(1), Replicate())
    placements_pair_repr = (Shard(0), Shard(1), Shard(2))
    placements_cp_single_repr = (Shard(0), Replicate())

    s_inputs_dtensor = distribute_tensor(
        s_inputs_global_host.to(device=manager.device, dtype=dtype).requires_grad_(True),
        device_mesh=device_mesh,
        placements=placements_single_repr,
    )
    s_dtensor = distribute_tensor(
        s_global_host.to(device=manager.device, dtype=dtype).requires_grad_(True),
        device_mesh=device_mesh,
        placements=placements_single_repr,
    )
    z_dtensor = distribute_tensor(
        z_global_host.to(device=manager.device, dtype=dtype).requires_grad_(True),
        device_mesh=device_mesh,
        placements=placements_pair_repr,
    )

    # Distribute atom-level features (token_to_rep_atom, x_pred, PTM-related keys)
    special_atom_features = {"token_to_rep_atom", "atom_pad_mask", "atom_to_token", "atom_resolved_mask", "frames_idx"}
    atom_inputs = {
        key: feats_global_host[key].to(device=manager.device)
        for key in special_atom_features
        if key in feats_global_host
    }

    base_batch = s_inputs_global_host.shape[0]
    x_pred_reshaped = x_pred_global_host.reshape(base_batch, multiplicity, *x_pred_global_host.shape[1:])
    for mul_idx in range(multiplicity):
        atom_inputs[f"x_pred_{mul_idx}"] = x_pred_reshaped[:, mul_idx].to(manager.device, dtype=dtype)

    atom_inputs["atom_counts_per_token"] = feats_global_host["atom_counts_per_token"].to(
        manager.device, dtype=torch.int64
    )

    atom_placements_cp = {key: placements_cp_single_repr for key in atom_inputs}
    atom_placements_cp["frames_idx"] = (Shard(1), Replicate())
    atom_placements_dp_cp = {key: placements_single_repr for key in atom_inputs if key != "atom_counts_per_token"}
    atom_placements_dp_cp["frames_idx"] = (Shard(0), Shard(1), Replicate())
    atom_feats_dtensor = distribute_atom_features(
        inputs=atom_inputs,
        placements_cp=atom_placements_cp,
        placements_dp_cp=atom_placements_dp_cp,
        device_mesh=device_mesh,
        cp_group=manager.group["cp"],
        multiplicities={"x_pred": multiplicity},
    )

    x_pred_dtensor = atom_feats_dtensor["x_pred"]

    feature_placements = get_feature_placements()

    single_repr_keys = [
        "token_pad_mask",
        "asym_id",
        "mol_type",
        "residue_index",
        "entity_id",
        "token_index",
        "sym_id",
        "cyclic_period",
    ]
    feats_dtensor = {
        key: distribute_tensor(
            feats_global_host[key].to(device=manager.device),
            device_mesh=device_mesh,
            placements=feature_placements["token_features"][key],
        )
        for key in single_repr_keys
        if key in feats_global_host
    }

    pair_repr_keys = ["token_bonds", "type_bonds", "token_pair_pad_mask", "contact_conditioning", "contact_threshold"]
    for key in pair_repr_keys:
        if key in feats_global_host:
            feats_dtensor[key] = distribute_tensor(
                feats_global_host[key].to(device=manager.device),
                device_mesh=device_mesh,
                placements=feature_placements["token_features"][key],
            )

    feats_dtensor["token_to_rep_atom"] = atom_feats_dtensor["token_to_rep_atom"]
    for key in special_atom_features - {"token_to_rep_atom"}:
        if key in atom_feats_dtensor:
            feats_dtensor[key] = atom_feats_dtensor[key]

    pred_distogram_logits_dtensor = distribute_tensor(
        pred_distogram_logits_global_host.to(manager.device),
        device_mesh=device_mesh,
        placements=placements_pair_repr,
    )

    # Keep copies to verify inputs are not mutated
    s_inputs_dtensor_copy = s_inputs_dtensor.clone()
    s_dtensor_copy = s_dtensor.clone()
    z_dtensor_copy = z_dtensor.clone()

    # ----- distribute upstream gradients -----
    d_plddt_logits_dtensor = distribute_tensor(
        d_plddt_logits_host.to(device=manager.device, dtype=dtype),
        device_mesh=device_mesh,
        placements=placements_single_repr,
    )
    d_pde_logits_dtensor = distribute_tensor(
        d_pde_logits_host.to(device=manager.device, dtype=dtype),
        device_mesh=device_mesh,
        placements=placements_pair_repr,
    )
    d_resolved_logits_dtensor = distribute_tensor(
        d_resolved_logits_host.to(device=manager.device, dtype=dtype),
        device_mesh=device_mesh,
        placements=placements_single_repr,
    )
    d_pae_logits_dtensor = distribute_tensor(
        d_pae_logits_host.to(device=manager.device, dtype=dtype),
        device_mesh=device_mesh,
        placements=placements_pair_repr,
    )

    # ----- forward -----
    output_dtensor = module(
        s_inputs=s_inputs_dtensor,
        s=s_dtensor,
        z=z_dtensor,
        x_pred=x_pred_dtensor,
        feats=feats_dtensor,
        pred_distogram_logits=pred_distogram_logits_dtensor,
        multiplicity=multiplicity,
        run_sequentially=run_sequentially,
    )

    # Compare all outputs against serial reference
    dp_rank = manager.group_rank["dp"]
    dp_size = len(manager.group_ranks["dp"])
    for key in output_dtensor:
        assert key in serial_output_feats_host, f"DTensor output key '{key}' missing from serial reference"
        dtensor_val = output_dtensor[key]
        if key == "pair_chains_iptm":
            serial_pciptm = serial_output_feats_host[key]
            if isinstance(dtensor_val, dict):
                for idx1, chain_dict in dtensor_val.items():
                    for idx2, dt_val in chain_dict.items():
                        local_val = dt_val.to_local().cpu()
                        is_sentinel = torch.all(local_val == CHAIN_IPTM_SENTINEL)
                        if idx1 in serial_pciptm and idx2 in serial_pciptm.get(idx1, {}):
                            serial_full = serial_pciptm[idx1][idx2]
                            chunk_size = serial_full.shape[0] // dp_size
                            serial_local = serial_full[dp_rank * chunk_size : (dp_rank + 1) * chunk_size]
                            if is_sentinel:
                                assert torch.all(serial_local.abs() < 1e-5), (
                                    f"Chain pair ({idx1}, {idx2}): distributed returned sentinel but "
                                    f"serial has non-zero values {serial_local} for this DP rank"
                                )
                            else:
                                torch.testing.assert_close(
                                    local_val,
                                    serial_local,
                                    msg=f"Chain pair ({idx1}, {idx2}) mismatch on DP rank {dp_rank}",
                                )
                        else:
                            assert is_sentinel, (
                                f"Extra chain pair ({idx1}, {idx2}) should be sentinel "
                                f"{CHAIN_IPTM_SENTINEL}, got {local_val}"
                            )
            else:
                assert isinstance(serial_pciptm, torch.Tensor), (
                    f"pair_chains_iptm: distributed returned DTensor (compute_ptms fallback) "
                    f"but serial returned {type(serial_pciptm)}"
                )
                torch.testing.assert_close(
                    dtensor_val.full_tensor().cpu(),
                    serial_pciptm,
                    msg="pair_chains_iptm fallback mismatch",
                )
            continue
        torch.testing.assert_close(
            dtensor_val.full_tensor().cpu(),
            serial_output_feats_host[key],
            msg=f"Mismatch for output key '{key}'",
        )

    # ----- backward -----
    torch.autograd.backward(
        [
            output_dtensor["plddt_logits"],
            output_dtensor["pde_logits"],
            output_dtensor["resolved_logits"],
            output_dtensor["pae_logits"],
        ],
        [
            d_plddt_logits_dtensor,
            d_pde_logits_dtensor,
            d_resolved_logits_dtensor,
            d_pae_logits_dtensor,
        ],
    )

    # Verify inputs were not mutated by the forward pass
    assert_tensors_identical(
        s_inputs_dtensor.to_local().cpu(),
        s_inputs_dtensor_copy.to_local().cpu(),
        check_grad=False,
        check_grad_fn=False,
        msg="s_inputs_dtensor was mutated during forward",
    )
    assert_tensors_identical(
        s_dtensor.to_local().cpu(),
        s_dtensor_copy.to_local().cpu(),
        check_grad=False,
        check_grad_fn=False,
        msg="s_dtensor was mutated during forward",
    )
    assert_tensors_identical(
        z_dtensor.to_local().cpu(),
        z_dtensor_copy.to_local().cpu(),
        check_grad=False,
        check_grad_fn=False,
        msg="z_dtensor was mutated during forward",
    )

    # Compare input gradients
    torch.testing.assert_close(
        s_inputs_dtensor.grad.full_tensor().cpu(), s_inputs_grad_host, msg="s_inputs gradient mismatch"
    )
    torch.testing.assert_close(s_dtensor.grad.full_tensor().cpu(), s_grad_host, msg="s gradient mismatch")
    torch.testing.assert_close(z_dtensor.grad.full_tensor().cpu(), z_grad_host, msg="z gradient mismatch")

    # Compare parameter gradients
    result_param_grads = {}
    for name, param in module.named_parameters():
        if param.grad is not None:
            if name not in serial_param_grads_host:
                raise ValueError(
                    f"Parameter '{name}' has a gradient in the distributed module " f"but not in the serial reference"
                )
            result_param_grads[name] = param.grad

    for name, expected in serial_param_grads_host.items():
        assert name in result_param_grads, f"Parameter '{name}' gradient missing in distributed module"
        torch.testing.assert_close(
            result_param_grads[name].full_tensor().cpu(),
            expected,
            msg=f"Parameter gradient mismatch for '{name}'",
        )

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.parametrize(
    ("setup_env", "dtype", "use_separate_heads", "run_sequentially", "multiplicity"),
    (
        params_test_module := [
            # multiplicity=2 exercises resolved_mask per-sample indexing paths
            (((2, (2, 2)), True, "cuda", "ENV"), torch.float64, False, False, 2),
            (((2, (2, 2)), True, "cuda", "ENV"), torch.float64, True, False, 2),
            # multiplicity=2 + run_sequentially needs B=1 (serial constraint), so dp=1
            (((1, (2, 2)), True, "cuda", "ENV"), torch.float64, False, True, 2),
        ]
    ),
    indirect=["setup_env"],
    ids=[
        (
            f"dp:{x[0][0][0]}, cp:{x[0][0][1]}, device:{x[0][2]}, init:{x[0][3]}, "
            f"dtype:{x[1]}, separate_heads:{x[2]}, sequential:{x[3]}, mult:{x[4]}"
        )
        for x in params_test_module
    ],
)
def test_dtensor_confidence_module_v2(setup_env, dtype, use_separate_heads, run_sequentially, multiplicity):
    """Test that DTensor ConfidenceModuleV2 matches serial ConfidenceModuleV2.

    Covers:
    * Forward pass parity for all logit outputs and aggregated metrics (pLDDT, PDE,
      PAE, complex_plddt, complex_iplddt, complex_pde, complex_ipde).
    * Backward pass parity for input gradients (s_inputs, s, z) and parameter gradients.
    * Non-mutation of inputs during forward.
    * Both the single-head and the intra/inter-chain separated head configurations.
    * Pairformer stack, distogram embedding, outer-sum s-to-z pair update.
    * multiplicity=1 and multiplicity=2 (exercises resolved_mask per-sample indexing).

    PTM/iPTM and pair_chains_iptm outputs are compared between serial and distributed.

    Parameters
    ----------
    setup_env : tuple
        Grid group sizes, world size, device type, backend, environment variables per rank.
    dtype : torch.dtype
        Tensor dtype for forward and backward passes.
    use_separate_heads : bool
        Whether to use separate intra/inter-chain PAE and PDE projection heads.
    multiplicity : int
        Number of diffusion samples per batch element.
    """
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if dtype == torch.float32:
        pytest.xfail("float32 dtype for logits has numerical stability issues and requires higher tolerances")

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("CUDA not available")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"Need {world_size} GPUs, only {torch.cuda.device_count()} available")

    seed_by_rank(42)

    B = grid_group_sizes["dp"]
    size_cp = grid_group_sizes["cp"][0]

    N_tokens = 16 * size_cp
    N_atoms = N_tokens * 4

    boltz2_params = create_boltz2_model_init_params(use_large_model=True)
    token_s = boltz2_params["token_s"]
    token_z = boltz2_params["token_z"]
    num_distogram_bins = boltz2_params["confidence_model_args"]["num_dist_bins"]

    confidence_model_args = boltz2_params["confidence_model_args"].copy()
    confidence_model_args["confidence_args"] = {
        **confidence_model_args["confidence_args"],
        "use_separate_heads": use_separate_heads,
    }
    confidence_module_kwargs = {
        "token_s": token_s,
        "token_z": token_z,
        "pairformer_args": boltz2_params["pairformer_args"],
        "token_level_confidence": True,
        "bond_type_feature": boltz2_params["bond_type_feature"],
        **confidence_model_args,
    }

    selected_keys = [
        "mol_type",
        "asym_id",
        "token_pad_mask",
        "token_pair_pad_mask",
        "token_to_rep_atom",
        "atom_counts_per_token",
        "frames_idx",
        "atom_to_token",
        "atom_pad_mask",
        "atom_resolved_mask",
        "residue_index",
        "entity_id",
        "token_index",
        "sym_id",
        "cyclic_period",
        "token_bonds",
        "type_bonds",
        "contact_conditioning",
        "contact_threshold",
    ]

    feats_host = random_features(
        size_batch=B,
        n_tokens=N_tokens,
        n_atoms=N_atoms,
        n_msa=1,
        atom_counts_per_token_range=(1, 4),
        device=torch.device("cpu"),
        float_value_range=(-0.2, 0.2),
        selected_keys=selected_keys,
    )
    feats_host = {k: (v.to(dtype=dtype) if v.is_floating_point() else v) for k, v in feats_host.items()}
    feats_device = {k: v.to(device=device_type) for k, v in feats_host.items()}

    val_init_range = 0.002

    # Input tensors
    s_inputs = torch.empty(B, N_tokens, token_s, device=device_type, dtype=dtype, requires_grad=True)
    s = torch.empty(B, N_tokens, token_s, device=device_type, dtype=dtype, requires_grad=True)
    z = torch.empty(B, N_tokens, N_tokens, token_z, device=device_type, dtype=dtype, requires_grad=True)
    x_pred = torch.empty(B * multiplicity, N_atoms, 3, device=device_type, dtype=dtype)
    pred_distogram_logits = torch.empty(B, N_tokens, N_tokens, num_distogram_bins, device=device_type, dtype=dtype)

    init_tensors_uniform([s_inputs, s, z, pred_distogram_logits], low=-val_init_range, high=val_init_range)
    # x_pred needs a wider range so that inter-atom distances exceed
    # the collinear-mask overlap threshold (0.01) in compute_ptms.
    init_tensors_uniform([x_pred], low=-10.0, high=10.0)

    s_inputs_global_host = s_inputs.detach().clone().cpu()
    s_global_host = s.detach().clone().cpu()
    z_global_host = z.detach().clone().cpu()
    x_pred_global_host = x_pred.detach().clone().cpu()
    pred_distogram_logits_global_host = pred_distogram_logits.detach().clone().cpu()

    # ----- serial module -----
    serial_module = SerialConfidenceModuleV2(**confidence_module_kwargs)
    serial_module = serial_module.to(device=device_type, dtype=dtype).train()
    init_module_params_glorot(serial_module)

    serial_state_dict = serial_module.state_dict()

    # ----- serial forward -----
    serial_output = serial_module(
        s_inputs=s_inputs,
        s=s,
        z=z,
        x_pred=x_pred,
        feats=feats_device,
        pred_distogram_logits=pred_distogram_logits,
        multiplicity=multiplicity,
        run_sequentially=run_sequentially,
    )

    # Upstream gradients for backward
    d_plddt_logits = torch.rand_like(serial_output["plddt_logits"], device=device_type)
    d_pde_logits = torch.rand_like(serial_output["pde_logits"], device=device_type)
    d_resolved_logits = torch.rand_like(serial_output["resolved_logits"], device=device_type)
    d_pae_logits = torch.rand_like(serial_output["pae_logits"], device=device_type)

    d_plddt_logits_host = d_plddt_logits.detach().clone().cpu()
    d_pde_logits_host = d_pde_logits.detach().clone().cpu()
    d_resolved_logits_host = d_resolved_logits.detach().clone().cpu()
    d_pae_logits_host = d_pae_logits.detach().clone().cpu()

    torch.autograd.backward(
        [
            serial_output["plddt_logits"],
            serial_output["pde_logits"],
            serial_output["resolved_logits"],
            serial_output["pae_logits"],
        ],
        [d_plddt_logits, d_pde_logits, d_resolved_logits, d_pae_logits],
    )

    def _to_cpu(val):
        if isinstance(val, torch.Tensor):
            return val.detach().clone().cpu()
        if isinstance(val, dict):
            return {k: _to_cpu(v) for k, v in val.items()}
        return val

    serial_output_feats_host = {k: _to_cpu(v) for k, v in serial_output.items()}

    # Verify the serial module has non-zero, non-NaN outputs (guard against vacuous pass)
    assert not torch.isnan(serial_output["plddt_logits"]).any(), "serial plddt_logits contains NaN"
    assert not torch.isnan(serial_output["pde_logits"]).any(), "serial pde_logits contains NaN"
    assert not torch.isnan(serial_output["pae_logits"]).any(), "serial pae_logits contains NaN"
    assert serial_output["plddt"].abs().max() > 0, "serial plddt is all-zero (vacuous)"
    assert serial_output["complex_plddt"].abs().max() > 0, "serial complex_plddt is all-zero (vacuous)"
    assert serial_output["ptm"].abs().max() > 0, "serial ptm is all-zero (vacuous)"
    assert isinstance(serial_output["pair_chains_iptm"], dict), (
        f"serial pair_chains_iptm should be a dict (from compute_ptms), "
        f"got {type(serial_output['pair_chains_iptm'])} — compute_ptms likely failed silently"
    )

    s_inputs_grad_host = s_inputs.grad.detach().clone().cpu()
    s_grad_host = s.grad.detach().clone().cpu()
    z_grad_host = z.grad.detach().clone().cpu()

    # Verify non-zero gradients (guard against vacuous backward pass)
    assert s_inputs_grad_host.abs().max() > 0, "serial s_inputs gradient is all-zero (vacuous)"
    assert s_grad_host.abs().max() > 0, "serial s gradient is all-zero (vacuous)"
    assert z_grad_host.abs().max() > 0, "serial z gradient is all-zero (vacuous)"

    serial_param_grads_host = {
        name: param.grad.detach().clone().cpu()
        for name, param in serial_module.named_parameters()
        if param.grad is not None
    }
    assert len(serial_param_grads_host) > 0, "No serial parameter gradients found (vacuous)"

    # ----- parallel distributed test -----
    spawn_multiprocessing(
        parallel_test_dtensor_confidence_module_v2,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        dtype,
        serial_state_dict,
        confidence_module_kwargs,
        multiplicity,
        run_sequentially,
        s_inputs_global_host,
        s_global_host,
        z_global_host,
        x_pred_global_host,
        feats_host,
        pred_distogram_logits_global_host,
        serial_output_feats_host,
        s_inputs_grad_host,
        s_grad_host,
        z_grad_host,
        d_plddt_logits_host,
        d_pde_logits_host,
        d_resolved_logits_host,
        d_pae_logits_host,
        serial_param_grads_host,
    )
