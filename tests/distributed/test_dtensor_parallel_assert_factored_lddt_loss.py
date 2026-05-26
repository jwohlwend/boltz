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

from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.layers.atom_to_token import reconstruct_atom_to_token_global
from boltz.distributed.model.loss.validation import factored_lddt_loss as dtensor_factored_lddt_loss
from boltz.model.loss.validation import factored_lddt_loss as serial_factored_lddt_loss
from boltz.testing.utils import distribute_atom_features, get_feature_placements, random_features, spawn_multiprocessing

# Get feature placements for the subset of features needed by factored_lddt_loss
_atom_keys = {"atom_pad_mask", "atom_to_token", "atom_counts_per_token"}
_token_keys = {"mol_type", "asym_id"}
_placements = get_feature_placements(atom_keys=_atom_keys, token_keys=_token_keys)
_placements_atom_features = _placements["atom_features"]
_placements_cp_atom_features = _placements["cp_atom_features"]
_placements_token_features = _placements["token_features"]

# Pred/true coords placements: [B*mult, n_atoms, 3]
_placements_pred_coords = {"pred_coords": (Shard(0), Shard(1), Replicate())}
_placements_true_coords = {"true_coords": (Shard(0), Shard(1), Replicate())}
_placements_cp_pred_coords = {"pred_coords": (Shard(0), Replicate())}
_placements_cp_true_coords = {"true_coords": (Shard(0), Replicate())}


def parallel_assert_factored_lddt_loss(rank, payload):
    (
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        multiplicity,
        feats_global_host,
        pred_coords_global_host,
        true_coords_global_host,
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
    pred_coords_unflat = pred_coords_global_host.unflatten(0, (size_batch, multiplicity))
    true_coords_unflat = true_coords_global_host.unflatten(0, (size_batch, multiplicity))
    for i_mul in range(multiplicity):
        inputs_atom[f"pred_coords_{i_mul}"] = pred_coords_unflat[:, i_mul].to(dtype=dtype)
        inputs_atom[f"true_coords_{i_mul}"] = true_coords_unflat[:, i_mul].to(dtype=dtype)

    placements_cp = dict(_placements_cp_atom_features)
    placements_dp_cp = dict(_placements_atom_features)
    for i_mul in range(multiplicity):
        placements_cp[f"pred_coords_{i_mul}"] = _placements_cp_pred_coords["pred_coords"]
        placements_cp[f"true_coords_{i_mul}"] = _placements_cp_true_coords["true_coords"]
        placements_dp_cp[f"pred_coords_{i_mul}"] = _placements_pred_coords["pred_coords"]
        placements_dp_cp[f"true_coords_{i_mul}"] = _placements_true_coords["true_coords"]

    feats_atom = distribute_atom_features(
        inputs_atom,
        placements_cp,
        placements_dp_cp,
        manager.device_mesh_subgroups,
        manager.group["cp"],
        multiplicities={"pred_coords": multiplicity, "true_coords": multiplicity},
    )

    pred_coords_dtensor = feats_atom["pred_coords"]
    true_coords_dtensor = feats_atom["true_coords"]
    atom_pad_mask_dtensor = feats_atom["atom_pad_mask"]
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

    pred_coords_local = _all_gather_single_repr(pred_coords_dtensor)
    true_coords_local = _all_gather_single_repr(true_coords_dtensor)
    atom_mask_base_local = _all_gather_single_repr(atom_pad_mask_dtensor)
    mol_type_local = _all_gather_single_repr(mol_type_dtensor)
    asym_id_local = _all_gather_single_repr(asym_id_dtensor)
    atom_to_token_local = reconstruct_atom_to_token_global(atom_to_token_dtensor)

    feats_local = {
        "atom_to_token": atom_to_token_local,
        "mol_type": mol_type_local,
        "asym_id": asym_id_local,
    }

    lddt_dict, total_dict = dtensor_factored_lddt_loss(
        true_atom_coords=true_coords_local,
        pred_atom_coords=pred_coords_local,
        feats=feats_local,
        atom_mask=atom_mask_base_local,
        multiplicity=multiplicity,
    )

    ref_slice = slice(local_start * multiplicity, local_end * multiplicity)
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


@pytest.mark.parametrize(
    "setup_env",
    [
        ((1, (1, 1)), True, "cuda", "ENV"),
        ((2, (2, 2)), True, "cuda", "ENV"),
    ],
    indirect=("setup_env",),
)
def test_dtensor_parallel_asserrt_factored_lddt_loss(setup_env):
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type != "cuda":
        pytest.skip("cdist_lddt requires CUDA")

    if torch.cuda.device_count() < world_size:
        pytest.skip(f"skip cuda test because torch.cuda.device_count() < {world_size}")

    torch.manual_seed(0)
    rng = torch.Generator(device="cpu")
    rng.manual_seed(0)
    size_batch = grid_group_sizes["dp"]
    cp_shape = grid_group_sizes["cp"]
    size_cp = cp_shape if isinstance(cp_shape, int) else int(torch.tensor(cp_shape).prod().item())
    n_tokens = size_cp * 20
    n_atoms = n_tokens * 20
    multiplicity = 2

    feats_global_host = random_features(
        size_batch=size_batch,
        n_tokens=n_tokens,
        n_atoms=n_atoms,
        n_msa=1,
        atom_counts_per_token_range=(1, 20),
        device=torch.device("cpu"),
        float_value_range=(0.0, 1.0),
        selected_keys=["atom_to_token", "atom_pad_mask", "atom_counts_per_token", "mol_type", "asym_id"],
        rng=rng,
    )
    pred_coords_global_host = torch.randn(size_batch * multiplicity, n_atoms, 3, dtype=torch.float32)
    true_coords_global_host = torch.randn(size_batch * multiplicity, n_atoms, 3, dtype=torch.float32)

    atom_mask_mult = feats_global_host["atom_pad_mask"].to(torch.float32).repeat_interleave(multiplicity, dim=0)
    feats_serial = {
        "atom_to_token": feats_global_host["atom_to_token"].to(dtype=torch.float32),
        "mol_type": feats_global_host["mol_type"],
        "asym_id": feats_global_host["asym_id"],
    }

    ref_lddt, ref_total = serial_factored_lddt_loss(
        true_atom_coords=true_coords_global_host,
        pred_atom_coords=pred_coords_global_host,
        feats=feats_serial,
        atom_mask=atom_mask_mult,
        multiplicity=multiplicity,
    )

    spawn_multiprocessing(
        parallel_assert_factored_lddt_loss,
        world_size,
        (
            grid_group_sizes,
            device_type,
            backend,
            env_per_rank,
            multiplicity,
            feats_global_host,
            pred_coords_global_host,
            true_coords_global_host,
            ref_lddt,
            ref_total,
        ),
    )
