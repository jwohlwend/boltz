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

from __future__ import annotations

import pytest
import torch
from torch.distributed.tensor import Replicate, Shard, distribute_tensor

import boltz.distributed.model.loss.validation as _validation_module
from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.loss.validation import get_lddt_metrics
from boltz.model.validation.validator import Validator
from boltz.testing.utils import distribute_atom_features, get_feature_placements, random_features, spawn_multiprocessing

_atom_keys = {"atom_pad_mask", "atom_to_token", "atom_counts_per_token"}
_token_keys = {"mol_type", "asym_id", "token_pad_mask"}
_placements = get_feature_placements(atom_keys=_atom_keys, token_keys=_token_keys)
_placements_atom_features = _placements["atom_features"]
_placements_cp_atom_features = _placements["cp_atom_features"]
_placements_token_features = _placements["token_features"]

_placements_pred_coords = {"pred_coords": (Shard(0), Shard(1), Replicate())}
_placements_cp_pred_coords = {"pred_coords": (Shard(0), Replicate())}


def parallel_assert_get_lddt_metrics(rank, payload):
    (
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        n_samples,
        K,
        expand_to_diffusion_samples,
        feats_global_host,
        pred_coords_global_host,
        true_coords_base_host,
        true_coords_global_host,
        true_coords_resolved_mask_base_host,
        true_coords_resolved_mask_host,
        ref_lddt,
        ref_total,
    ) = payload

    monkeypatch = pytest.MonkeyPatch()
    if env_per_rank is not None:
        for var_name, value in env_per_rank.items():
            if value == "<INPUT_RANK>":
                monkeypatch.setenv(var_name, f"{rank}")
                continue
            monkeypatch.setenv(var_name, value)

    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    manager = DistributedManager()

    device = manager.device
    dtype = torch.float32

    size_batch = feats_global_host["atom_pad_mask"].shape[0]
    rank_dp = manager.group_rank["dp"]
    num_dp_ranks = grid_group_sizes["dp"]
    local_batch_size = size_batch // num_dp_ranks
    local_start = rank_dp * local_batch_size
    local_end = local_start + local_batch_size

    def _all_gather_single_repr(single_dtensor):
        single_dtensor = single_dtensor.redistribute(
            placements=[Shard(0), Replicate(), Replicate()],
        )
        return single_dtensor.to_local()

    inputs_atom = {
        k: v.to(dtype=dtype if v.dtype.is_floating_point else v.dtype)
        for k, v in feats_global_host.items()
        if k in _placements_cp_atom_features
    }
    pred_coords_unflat = pred_coords_global_host.unflatten(0, (size_batch, n_samples))
    for i_mul in range(n_samples):
        inputs_atom[f"pred_coords_{i_mul}"] = pred_coords_unflat[:, i_mul].to(dtype=dtype)

    placements_cp = dict(_placements_cp_atom_features)
    placements_dp_cp = dict(_placements_atom_features)
    for i_mul in range(n_samples):
        placements_cp[f"pred_coords_{i_mul}"] = _placements_cp_pred_coords["pred_coords"]
        placements_dp_cp[f"pred_coords_{i_mul}"] = _placements_pred_coords["pred_coords"]

    feats_atom = distribute_atom_features(
        inputs_atom,
        placements_cp,
        placements_dp_cp,
        manager.device_mesh_subgroups,
        manager.group["cp"],
        multiplicities={"pred_coords": n_samples},
    )

    atom_to_token_dtensor = feats_atom["atom_to_token"]

    mol_type_dtensor = distribute_tensor(
        feats_global_host["mol_type"].to(device=device, dtype=torch.int64),
        device_mesh=manager.device_mesh_subgroups,
        placements=_placements_token_features["mol_type"],
    )
    asym_id_dtensor = distribute_tensor(
        feats_global_host["asym_id"].to(device=device, dtype=torch.int64),
        device_mesh=manager.device_mesh_subgroups,
        placements=_placements_token_features["asym_id"],
    )

    mol_type_local = _all_gather_single_repr(mol_type_dtensor)
    asym_id_local = _all_gather_single_repr(asym_id_dtensor)
    pred_coords_local = _all_gather_single_repr(feats_atom["pred_coords"])
    atom_pad_mask_local = _all_gather_single_repr(feats_atom["atom_pad_mask"]).bool()

    local_mul_start = local_start * n_samples
    local_mul_end = local_end * n_samples
    if expand_to_diffusion_samples:
        true_coords_local_unpadded = true_coords_global_host[local_mul_start:local_mul_end].to(
            device=device, dtype=dtype
        )
        mask_local_unpadded = true_coords_resolved_mask_host[local_mul_start:local_mul_end].to(
            device=device, dtype=dtype
        )
        active_unpadded_mask = mask_local_unpadded[0].bool()
    else:
        true_coords_local_unpadded = (
            true_coords_base_host[local_start:local_end].squeeze(0).to(device=device, dtype=dtype)
        )
        mask_local_unpadded = (
            true_coords_resolved_mask_base_host[local_start:local_end].squeeze(0).to(device=device, dtype=dtype)
        )
        active_unpadded_mask = mask_local_unpadded.bool()

    atom_mask_row = atom_pad_mask_local[0].bool()
    n_atoms_padded = atom_mask_row.shape[0]
    n_atoms_active = int(active_unpadded_mask.sum().item())
    if int(atom_mask_row.sum().item()) != n_atoms_active:
        raise ValueError(
            "atom_pad_mask/padded atom-space mismatch: "
            f"sum(atom_pad_mask)={int(atom_mask_row.sum().item())}, n_atoms_active_unpadded={n_atoms_active}"
        )
    if expand_to_diffusion_samples:
        true_coords_local_active = true_coords_local_unpadded[:, :, active_unpadded_mask, :]
        mask_local_active = mask_local_unpadded[:, active_unpadded_mask]

        true_coords_local = torch.zeros(
            true_coords_local_unpadded.shape[0],
            true_coords_local_unpadded.shape[1],
            n_atoms_padded,
            true_coords_local_unpadded.shape[3],
            device=device,
            dtype=dtype,
        )
        true_coords_local[:, :, atom_mask_row, :] = true_coords_local_active

        mask_local = torch.zeros(
            mask_local_unpadded.shape[0],
            n_atoms_padded,
            device=device,
            dtype=dtype,
        )
        mask_local[:, atom_mask_row] = mask_local_active
    else:
        true_coords_local_active = true_coords_local_unpadded[:, active_unpadded_mask, :]
        mask_local_active = mask_local_unpadded[active_unpadded_mask]

        true_coords_local = torch.zeros(
            true_coords_local_unpadded.shape[0],
            n_atoms_padded,
            true_coords_local_unpadded.shape[2],
            device=device,
            dtype=dtype,
        )
        true_coords_local[:, atom_mask_row, :] = true_coords_local_active

        mask_local = torch.zeros(
            n_atoms_padded,
            device=device,
            dtype=dtype,
        )
        mask_local[atom_mask_row] = mask_local_active

    lddt_dict, total_dict = get_lddt_metrics(
        atom_to_token_dtensor=atom_to_token_dtensor,
        num_conformers=K,
        n_samples=n_samples,
        true_coords=true_coords_local,
        true_coords_resolved_mask=mask_local,
        mol_type=mol_type_local,
        asym_id=asym_id_local,
        sample_atom_coords=pred_coords_local,
        expand_to_diffusion_samples=expand_to_diffusion_samples,
    )

    ref_slice = slice(local_mul_start, local_mul_end)
    for key in ref_lddt:
        torch.testing.assert_close(
            lddt_dict[key].cpu(),
            ref_lddt[key][ref_slice],
        )
        torch.testing.assert_close(
            total_dict[key].cpu(),
            ref_total[key][ref_slice],
        )

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.parametrize("expand_to_diffusion_samples", [True, False])
@pytest.mark.parametrize(
    "setup_env",
    [
        ((1, (3, 3)), True, "cuda", "ENV"),
        ((2, (2, 2)), True, "cuda", "ENV"),
    ],
    indirect=("setup_env",),
)
def test_get_lddt_metrics(setup_env, expand_to_diffusion_samples):
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type != "cuda":
        pytest.skip("cdist_lddt requires CUDA")

    if torch.cuda.device_count() < world_size:
        pytest.skip(f"skip cuda test because torch.cuda.device_count() < {world_size}")

    torch.manual_seed(0)
    rng = torch.Generator(device="cpu")
    rng.manual_seed(0)
    size_batch = grid_group_sizes["dp"]
    size_cp = grid_group_sizes["cp"][0]
    n_tokens = size_cp * 20
    n_atoms = n_tokens * 20
    n_samples = 2
    K = 2

    feats_global_host = random_features(
        size_batch=size_batch,
        n_tokens=n_tokens,
        n_atoms=n_atoms,
        n_msa=1,
        atom_counts_per_token_range=(1, 20),
        device=torch.device("cpu"),
        float_value_range=(0.0, 1.0),
        selected_keys=[
            "atom_to_token",
            "atom_pad_mask",
            "atom_counts_per_token",
            "mol_type",
            "asym_id",
            "token_pad_mask",
        ],
        rng=rng,
    )
    token_pad_mask = torch.ones((size_batch, n_tokens), dtype=torch.bool)
    token_pad_mask[:, ::3] = False
    token_pad_mask[:, :2] = True
    feats_global_host["token_pad_mask"] = token_pad_mask.to(dtype=feats_global_host["token_pad_mask"].dtype)

    atom_to_token = feats_global_host["atom_to_token"]
    atom_pad_mask = torch.zeros((size_batch, n_atoms), dtype=torch.bool)
    for batch_idx in range(size_batch):
        token_mask = token_pad_mask[batch_idx]
        atom_mask = atom_to_token[batch_idx][:, token_mask].any(dim=1)
        atom_pad_mask[batch_idx] = atom_mask
        atom_to_token[batch_idx, ~atom_mask, :] = 0
        atom_to_token[batch_idx, :, ~token_mask] = 0

    feats_global_host["atom_pad_mask"] = atom_pad_mask.to(dtype=feats_global_host["atom_pad_mask"].dtype)
    pred_coords_global_host = torch.randn(size_batch * n_samples, n_atoms, 3, dtype=torch.float32)
    true_coords_base_host = torch.randn(size_batch, K, n_atoms, 3, dtype=torch.float32)
    true_coords_global_host = true_coords_base_host.repeat_interleave(n_samples, dim=0)
    true_coords_resolved_mask_base_host = feats_global_host["atom_pad_mask"].to(torch.float32)
    true_coords_resolved_mask_host = true_coords_resolved_mask_base_host.repeat_interleave(n_samples, dim=0)

    ref_lddt: dict[str, torch.Tensor] = {}
    ref_total: dict[str, torch.Tensor] = {}
    feats_serial = {
        "atom_to_token": feats_global_host["atom_to_token"].to(dtype=torch.float32),
        "mol_type": feats_global_host["mol_type"],
        "asym_id": feats_global_host["asym_id"],
        "coords": true_coords_base_host,
    }

    if expand_to_diffusion_samples:
        ref_lddt, ref_total = Validator.get_lddt_metrics(
            None,
            model=None,
            batch=feats_serial,
            out={"sample_atom_coords": pred_coords_global_host},
            idx_dataset=0,
            n_samples=n_samples,
            true_coords_resolved_mask=true_coords_resolved_mask_host,
            true_coords=true_coords_global_host,
            expand_to_diffusion_samples=True,
        )
    else:
        pred_coords_by_batch = pred_coords_global_host.unflatten(0, (size_batch, n_samples))
        for batch_idx in range(size_batch):
            feats_serial_single = {
                "atom_to_token": feats_serial["atom_to_token"][batch_idx : batch_idx + 1],
                "mol_type": feats_serial["mol_type"][batch_idx : batch_idx + 1],
                "asym_id": feats_serial["asym_id"][batch_idx : batch_idx + 1],
                "coords": feats_serial["coords"][batch_idx : batch_idx + 1],
            }
            lddt_single, total_single = Validator.get_lddt_metrics(
                None,
                model=None,
                batch=feats_serial_single,
                out={"sample_atom_coords": pred_coords_by_batch[batch_idx]},
                idx_dataset=0,
                n_samples=n_samples,
                true_coords_resolved_mask=true_coords_resolved_mask_base_host[batch_idx],
                true_coords=true_coords_base_host[batch_idx],
                expand_to_diffusion_samples=False,
            )
            if not ref_lddt:
                for key in lddt_single:
                    ref_lddt[key] = torch.zeros(size_batch * n_samples, K, dtype=lddt_single[key].dtype)
                    ref_total[key] = torch.zeros(size_batch * n_samples, K, dtype=total_single[key].dtype)
            row_slice = slice(batch_idx * n_samples, (batch_idx + 1) * n_samples)
            for key in lddt_single:
                ref_lddt[key][row_slice] = lddt_single[key]
                ref_total[key][row_slice] = total_single[key]

    spawn_multiprocessing(
        parallel_assert_get_lddt_metrics,
        world_size,
        (
            grid_group_sizes,
            device_type,
            backend,
            env_per_rank,
            n_samples,
            K,
            expand_to_diffusion_samples,
            feats_global_host,
            pred_coords_global_host,
            true_coords_base_host,
            true_coords_global_host,
            true_coords_resolved_mask_base_host,
            true_coords_resolved_mask_host,
            dict(ref_lddt),
            dict(ref_total),
        ),
    )


@pytest.mark.parametrize(
    "mutation, expand, match",
    [
        ("batch_size", True, "local batch size 1"),
        ("sample_batch", True, "sample_atom_coords batch must equal"),
        ("sample_ndim", True, "sample_atom_coords must be rank 3"),
        ("true_coords_ndim_expanded", True, "true_coords must be rank 4"),
        ("true_coords_ndim_not_expanded", False, "true_coords must be rank 3"),
        ("mask_rank_expanded", True, "true_coords_resolved_mask must be rank 2"),
        ("mask_rank_not_expanded", False, "true_coords_resolved_mask must be rank 1"),
        ("true_coords_K_expanded", True, "true_coords conformer count"),
        ("true_coords_K_not_expanded", False, "true_coords conformer count"),
        ("true_coords_batch_expanded", True, "true_coords batch dim"),
        ("mol_type_tokens", True, "mol_type N_tokens"),
        ("asym_id_tokens", True, "asym_id N_tokens"),
        ("sample_atoms", True, "sample_atom_coords N_atoms"),
        ("mask_atoms_expanded", True, "true_coords_resolved_mask N_atoms"),
        ("mask_atoms_not_expanded", False, "true_coords_resolved_mask N_atoms"),
        ("true_coords_atoms_expanded", True, "true_coords N_atoms"),
        ("true_coords_atoms_not_expanded", False, "true_coords N_atoms"),
    ],
)
def test_get_lddt_metrics_shape_errors(monkeypatch, mutation, expand, match):
    B, N_tokens, N_atoms, K, n_samples = 1, 4, 8, 2, 2

    atom_to_token = torch.zeros(B, N_atoms, N_tokens)
    mol_type = torch.zeros(B, N_tokens, dtype=torch.long)
    asym_id = torch.zeros(B, N_tokens, dtype=torch.long)
    sample_atom_coords = torch.zeros(B * n_samples, N_atoms, 3)

    if expand:
        true_coords = torch.zeros(B * n_samples, K, N_atoms, 3)
        mask = torch.zeros(B * n_samples, N_atoms)
    else:
        true_coords = torch.zeros(K, N_atoms, 3)
        mask = torch.zeros(N_atoms)

    if mutation == "batch_size":
        atom_to_token = torch.zeros(2, N_atoms, N_tokens)
    elif mutation == "sample_batch":
        sample_atom_coords = torch.zeros(B * n_samples + 1, N_atoms, 3)
    elif mutation == "sample_ndim":
        sample_atom_coords = torch.zeros(B * n_samples * N_atoms, 3)
    elif mutation == "true_coords_ndim_expanded":
        true_coords = torch.zeros(K, N_atoms, 3)
    elif mutation == "true_coords_ndim_not_expanded":
        true_coords = torch.zeros(B * n_samples, K, N_atoms, 3)
    elif mutation == "mask_rank_expanded":
        mask = torch.zeros(N_atoms)
    elif mutation == "mask_rank_not_expanded":
        mask = torch.zeros(B * n_samples, N_atoms)
    elif mutation == "true_coords_K_expanded":
        true_coords = torch.zeros(B * n_samples, K + 1, N_atoms, 3)
    elif mutation == "true_coords_K_not_expanded":
        true_coords = torch.zeros(K + 1, N_atoms, 3)
    elif mutation == "true_coords_batch_expanded":
        true_coords = torch.zeros(B * n_samples + 1, K, N_atoms, 3)
    elif mutation == "mol_type_tokens":
        mol_type = torch.zeros(B, N_tokens + 1, dtype=torch.long)
    elif mutation == "asym_id_tokens":
        asym_id = torch.zeros(B, N_tokens + 1, dtype=torch.long)
    elif mutation == "sample_atoms":
        sample_atom_coords = torch.zeros(B * n_samples, N_atoms + 1, 3)
    elif mutation == "mask_atoms_expanded":
        mask = torch.zeros(B * n_samples, N_atoms + 1)
    elif mutation == "mask_atoms_not_expanded":
        mask = torch.zeros(N_atoms + 1)
    elif mutation == "true_coords_atoms_expanded":
        true_coords = torch.zeros(B * n_samples, K, N_atoms + 1, 3)
    elif mutation == "true_coords_atoms_not_expanded":
        true_coords = torch.zeros(K, N_atoms + 1, 3)

    monkeypatch.setattr(_validation_module, "reconstruct_atom_to_token_global", lambda _: atom_to_token)

    with pytest.raises(ValueError, match=match):
        get_lddt_metrics(
            atom_to_token_dtensor=None,
            num_conformers=K,
            n_samples=n_samples,
            true_coords=true_coords,
            true_coords_resolved_mask=mask,
            mol_type=mol_type,
            asym_id=asym_id,
            sample_atom_coords=sample_atom_coords,
            expand_to_diffusion_samples=expand,
        )
