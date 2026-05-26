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

import types

import pytest
import torch
from torch.distributed.tensor import DTensor, Replicate, Shard

from boltz.data.mol import minimum_lddt_symmetry_coords as serial_minimum_lddt_symmetry_coords
from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.models.boltz2 import Boltz2 as DistributedBoltz2
from boltz.testing.utils import (
    distribute_atom_features,
    get_feature_placements,
    random_features,
    spawn_multiprocessing,
)


def _build_symmetry_features_for_batch(
    all_coords_batch: torch.Tensor,
    all_resolved_mask_batch: torch.Tensor,
    crop_to_all_atom_map_batch: torch.Tensor,
    chain_swaps_batch: list | None = None,
):
    """Build symmetry features for a batch of samples (Boltz-2 semantics)."""
    batch_size = all_coords_batch.shape[0]
    if chain_swaps_batch is None:
        chain_swaps_batch = [[[]] for _ in range(batch_size)]
    amino_acids_symmetries_batch = [[] for _ in range(batch_size)]
    ligand_symmetries_batch = [[] for _ in range(batch_size)]

    return {
        "all_coords": all_coords_batch,
        "all_resolved_mask": all_resolved_mask_batch,
        "crop_to_all_atom_map": crop_to_all_atom_map_batch,
        "chain_swaps": chain_swaps_batch,
        "amino_acids_symmetries": amino_acids_symmetries_batch,
        "ligand_symmetries": ligand_symmetries_batch,
    }


def _make_two_chain_swaps(n_atoms: int) -> list:
    """Build chain_swaps for one sample with two equal-length chains that can be swapped."""
    half = n_atoms // 2
    identity = []
    swap_ab = [
        (0, half, half, 2 * half, 0, 1),
        (half, 2 * half, 0, half, 1, 0),
    ]
    return [identity, swap_ab]


_atom_keys = {"atom_pad_mask", "coords", "atom_resolved_mask"}
_placements = get_feature_placements(atom_keys=_atom_keys, token_keys=set())
_placements_atom_features = _placements["atom_features"]
_placements_cp_atom_features = _placements["cp_atom_features"]

_placements_sample_coords = (Shard(0), Shard(1), Replicate())
_placements_cp_sample_coords = (Shard(0), Replicate())

_placements_token_index = (Shard(0), Replicate(), Replicate())
_placements_cp_token_index = (Shard(0), Replicate())


def parallel_assert_get_true_coordinates(rank, payload):
    """Test get_true_coordinates: symmetry path (parity) + non-symmetry path (types/shapes).

    Symmetry path: DP sharding of symmetry features, compares DTensor output against
    serial minimum_lddt_symmetry_coords reference for numerical parity.

    Non-symmetry path: verifies outputs are DTensors (not plain tensors), shapes are
    correct for expanded and unexpanded modes.
    """
    (
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        multiplicity,
        input_feats_global_host,
        sample_coords_global_host,
        feats_symmetry_global_host,
        expected_true_coords_per_mult_sample,
        expected_true_mask_per_mult_sample,
        use_nontrivial_swaps,
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

    size_batch = input_feats_global_host["atom_pad_mask"].shape[0]
    rank_dp = manager.group_rank["dp"]

    inputs_atom = {
        k: v.to(dtype=dtype if v.dtype.is_floating_point else v.dtype)
        for k, v in input_feats_global_host.items()
        if k in _placements_cp_atom_features
    }
    sample_coords_unflat = sample_coords_global_host.unflatten(0, (size_batch, multiplicity))
    for i_mul in range(multiplicity):
        inputs_atom[f"sample_coords_{i_mul}"] = sample_coords_unflat[:, i_mul].to(dtype=dtype)

    placements_cp = _placements_cp_atom_features | {
        f"sample_coords_{i_mul}": _placements_cp_sample_coords for i_mul in range(multiplicity)
    }
    placements_dp_cp = _placements_atom_features | {
        f"sample_coords_{i_mul}": _placements_sample_coords for i_mul in range(multiplicity)
    }

    feats_atom = distribute_atom_features(
        inputs_atom,
        placements_cp,
        placements_dp_cp,
        manager.device_mesh_subgroups,
        manager.group["cp"],
        multiplicities={"sample_coords": multiplicity},
    )

    coords_dtensor = feats_atom["sample_coords"]
    atom_pad_mask_dtensor = feats_atom["atom_pad_mask"]
    n_atoms_padded = atom_pad_mask_dtensor.shape[1]
    assert coords_dtensor.shape == (
        size_batch * multiplicity,
        n_atoms_padded,
        3,
    ), f"coords must be shape ({size_batch * multiplicity}, {n_atoms_padded}, 3) but got {coords_dtensor.shape}"

    # token_index: needed by get_true_coordinates to determine local_batch_size
    from torch.distributed.tensor import distribute_tensor

    token_index_global = input_feats_global_host["token_index"].to(device)
    token_index_dtensor = distribute_tensor(
        token_index_global,
        device_mesh=manager.device_mesh_subgroups,
        placements=_placements_token_index,
    )

    num_dp_ranks = grid_group_sizes["dp"]
    batch_size_global = size_batch
    coords_global_batch = size_batch * multiplicity
    local_batch_size = batch_size_global // num_dp_ranks
    local_start = rank_dp * local_batch_size
    local_end = local_start + local_batch_size

    batch_local = {
        "token_index": token_index_dtensor,
        "all_coords": feats_symmetry_global_host["all_coords"][local_start:local_end].to(device),
        "all_resolved_mask": feats_symmetry_global_host["all_resolved_mask"][local_start:local_end].to(device),
        "crop_to_all_atom_map": feats_symmetry_global_host["crop_to_all_atom_map"][local_start:local_end].to(device),
        "chain_swaps": feats_symmetry_global_host["chain_swaps"][local_start:local_end],
        "amino_acids_symmetries": feats_symmetry_global_host["amino_acids_symmetries"][local_start:local_end],
        "ligand_symmetries": feats_symmetry_global_host["ligand_symmetries"][local_start:local_end],
        "atom_pad_mask": atom_pad_mask_dtensor,
    }
    out_dtensor = {"sample_atom_coords": coords_dtensor}

    dummy = types.SimpleNamespace()

    result = DistributedBoltz2.get_true_coordinates(
        dummy,
        batch=batch_local,
        out=out_dtensor,
        diffusion_samples=multiplicity,
        symmetry_correction=True,
    )

    true_coords_dtensor = result["true_coords"]
    true_mask_dtensor = result["true_coords_resolved_mask"]

    assert (
        true_coords_dtensor.shape[0] == coords_global_batch
    ), f"true_coords_dtensor.shape[0] should be {coords_global_batch}, got {true_coords_dtensor.shape[0]}"
    assert (
        true_mask_dtensor.shape[0] == coords_global_batch
    ), f"true_mask_dtensor.shape[0] should be {coords_global_batch}, got {true_mask_dtensor.shape[0]}"

    actual_coords_full = true_coords_dtensor.full_tensor().cpu()
    # With symmetric correction, true_coords = shardwise_unsqueeze(true_coords, dim=1) is used to fix the shape mismatch, which gives us a 4D tensor.
    if actual_coords_full.ndim == 4:
        actual_coords_full = actual_coords_full.squeeze(1)
    actual_mask_full = true_mask_dtensor.full_tensor().cpu()
    atom_pad_mask_full = atom_pad_mask_dtensor.full_tensor().cpu()

    for mult_idx in range(coords_global_batch):
        expected_coords = expected_true_coords_per_mult_sample[mult_idx]
        expected_mask = expected_true_mask_per_mult_sample[mult_idx]

        batch_idx = mult_idx // multiplicity
        real_atom_mask = atom_pad_mask_full[batch_idx].bool()

        actual_coords_no_pad = actual_coords_full[mult_idx, real_atom_mask, :]
        actual_mask_no_pad = actual_mask_full[mult_idx, real_atom_mask]

        torch.testing.assert_close(
            actual_coords_no_pad,
            expected_coords,
            msg=f"Sample {mult_idx}: true_coords mismatch",
        )
        torch.testing.assert_close(
            actual_mask_no_pad,
            expected_mask.squeeze(0) if expected_mask.ndim > 1 else expected_mask,
            msg=f"Sample {mult_idx}: true_mask mismatch",
        )

    assert result["rmsds"] == 0, f"rmsds should be 0, got {result['rmsds']}"
    assert result["best_rmsd_recall"] == 0, f"best_rmsd_recall should be 0, got {result['best_rmsd_recall']}"

    # ---- Non-symmetry path: type, shape, and expand_to_diffusion_samples ----
    batch_local["coords"] = feats_atom["coords"]
    batch_local["atom_resolved_mask"] = feats_atom["atom_resolved_mask"]

    result_nosym = DistributedBoltz2.get_true_coordinates(
        dummy,
        batch=batch_local,
        out=out_dtensor,
        diffusion_samples=multiplicity,
        symmetry_correction=False,
        expand_to_diffusion_samples=True,
    )
    tc_nosym = result_nosym["true_coords"]
    tm_nosym = result_nosym["true_coords_resolved_mask"]

    assert isinstance(tc_nosym, DTensor), f"Rank {rank}: non-sym true_coords should be DTensor, got {type(tc_nosym)}"
    assert isinstance(tm_nosym, DTensor), f"Rank {rank}: non-sym true_mask should be DTensor, got {type(tm_nosym)}"
    assert tc_nosym.shape == (
        coords_global_batch,
        n_atoms_padded,
        3,
    ), f"Rank {rank}: expanded true_coords shape {tc_nosym.shape} != ({coords_global_batch}, {n_atoms_padded}, 3)"
    assert tm_nosym.shape == (
        coords_global_batch,
        n_atoms_padded,
    ), f"Rank {rank}: expanded true_mask shape {tm_nosym.shape} != ({coords_global_batch}, {n_atoms_padded})"

    result_unexpanded = DistributedBoltz2.get_true_coordinates(
        dummy,
        batch=batch_local,
        out=out_dtensor,
        diffusion_samples=multiplicity,
        symmetry_correction=False,
        expand_to_diffusion_samples=False,
    )
    tc_unexpanded = result_unexpanded["true_coords"]
    assert isinstance(tc_unexpanded, DTensor), f"Rank {rank}: unexpanded true_coords should be DTensor"
    assert (
        tc_unexpanded.shape[0] == batch_size_global
    ), f"Rank {rank}: unexpanded batch dim {tc_unexpanded.shape[0]} should be {batch_size_global}"
    assert (
        tc_unexpanded.shape[0] < tc_nosym.shape[0]
    ), f"Rank {rank}: unexpanded batch {tc_unexpanded.shape[0]} should be < expanded {tc_nosym.shape[0]}"

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.skipif(not torch.cuda.is_available(), reason="CUDA not available")
@pytest.mark.parametrize("use_nontrivial_swaps", [False, True], ids=["identity_swaps", "nontrivial_swaps"])
@pytest.mark.parametrize(
    "setup_env",
    [((2, (2, 2)), True, "cuda", "ENV")],  # dp=2, cp=(2,2), world_size=8
    indirect=("setup_env",),
)
def test_dtensor_get_true_coordinates(setup_env, use_nontrivial_swaps):
    """Test get_true_coordinates: symmetry correction parity + non-symmetry shape/type checks.

    Symmetry path: dp=2, multiplicity=2, compares against serial
    minimum_lddt_symmetry_coords for numerical parity. Parametrized over
    identity and nontrivial chain swaps.

    Non-symmetry path: verifies outputs are DTensors with correct shapes,
    and that expand_to_diffusion_samples=True/False produces correct batch dims.
    """
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if torch.cuda.device_count() < world_size:
        pytest.skip(f"skip cuda test because torch.cuda.device_count() < {world_size}")

    torch.manual_seed(123)

    size_ring = grid_group_sizes["cp"][0]
    num_dp_ranks = grid_group_sizes["dp"]
    batch_size_per_dp_rank = 1
    batch_size_global = batch_size_per_dp_rank * num_dp_ranks
    multiplicity = 2

    n_atoms_per_token = 3
    n_tokens = size_ring * 4
    n_atoms = n_atoms_per_token * n_tokens

    feats_from_random = random_features(
        size_batch=batch_size_global,
        n_tokens=n_tokens,
        n_atoms=n_atoms,
        n_msa=1,
        atom_counts_per_token_range=(1, n_atoms_per_token),
        device=torch.device("cpu"),
        float_value_range=(-1.0, 1.0),
        selected_keys=["atom_pad_mask", "token_index", "atom_counts_per_token", "coords", "atom_resolved_mask"],
    )

    atom_pad_mask_global = feats_from_random["atom_pad_mask"]
    token_index_global = feats_from_random["token_index"]
    atom_counts_per_token_global = feats_from_random["atom_counts_per_token"]

    sample_coords_global = torch.randn((batch_size_global * multiplicity, n_atoms, 3), dtype=torch.float32)

    input_feats_global = {
        "atom_pad_mask": atom_pad_mask_global,
        "token_index": token_index_global,
        "atom_counts_per_token": atom_counts_per_token_global,
        "coords": feats_from_random["coords"],
        "atom_resolved_mask": feats_from_random["atom_resolved_mask"],
    }

    coords_for_symmetry = sample_coords_global[::multiplicity]
    all_coords_global = coords_for_symmetry.clone()
    all_resolved_mask_global = torch.ones((batch_size_global, n_atoms), dtype=torch.bool)
    crop_to_all_atom_map_global = (
        torch.arange(n_atoms, dtype=torch.long).unsqueeze(0).expand(batch_size_global, -1).contiguous()
    )

    chain_swaps_batch = None
    if use_nontrivial_swaps:
        chain_swaps_batch = [_make_two_chain_swaps(n_atoms) for _ in range(batch_size_global)]

    feats_symmetry_global = _build_symmetry_features_for_batch(
        all_coords_global,
        all_resolved_mask_global,
        crop_to_all_atom_map_global,
        chain_swaps_batch=chain_swaps_batch,
    )

    expected_true_coords_per_mult_sample = []
    expected_true_mask_per_mult_sample = []

    for idx in range(batch_size_global):
        for rep in range(multiplicity):
            i = idx * multiplicity + rep
            expected_coords, expected_mask = serial_minimum_lddt_symmetry_coords(
                coords=sample_coords_global[i : i + 1],
                feats=feats_symmetry_global,
                index_batch=idx,
            )
            expected_true_coords_per_mult_sample.append(expected_coords.squeeze(0).detach().clone().cpu())
            expected_true_mask_per_mult_sample.append(expected_mask.detach().clone().cpu())

    input_feats_global_host = {k: v.detach().clone().cpu() for k, v in input_feats_global.items()}
    sample_coords_global_host = sample_coords_global.detach().clone().cpu()

    payload = (
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        multiplicity,
        input_feats_global_host,
        sample_coords_global_host,
        feats_symmetry_global,
        expected_true_coords_per_mult_sample,
        expected_true_mask_per_mult_sample,
        use_nontrivial_swaps,
    )

    spawn_multiprocessing(parallel_assert_get_true_coordinates, world_size, payload)
