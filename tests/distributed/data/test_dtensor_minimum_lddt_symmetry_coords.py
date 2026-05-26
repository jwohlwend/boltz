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
from torch.distributed.tensor import Replicate, Shard

from boltz.data.mol import minimum_lddt_symmetry_coords as serial_minimum_lddt_symmetry_coords
from boltz.distributed.data.feature.symmetry import (
    minimum_lddt_symmetry_coords as dtensor_minimum_lddt_symmetry_coords,
)
from boltz.distributed.manager import DistributedManager
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
    """Build symmetry features for a batch of samples.

    Parameters
    ----------
    chain_swaps_batch : list or None
        If None, each sample gets a single identity (no-op) chain swap ``[[]]``.
        Otherwise, provide a per-sample list of swap combinations.

    """
    batch_size = all_coords_batch.shape[0]
    if chain_swaps_batch is None:
        chain_swaps_batch = [[[]] for _ in range(batch_size)]
    amino_acids_symmetries_batch = [[] for _ in range(batch_size)]
    ligand_symmetries_batch = [[] for _ in range(batch_size)]

    feats = {
        "all_coords": all_coords_batch,
        "all_resolved_mask": all_resolved_mask_batch,
        "crop_to_all_atom_map": crop_to_all_atom_map_batch,
        "chain_swaps": chain_swaps_batch,
        "amino_acids_symmetries": amino_acids_symmetries_batch,
        "ligand_symmetries": ligand_symmetries_batch,
    }
    return feats


def _make_two_chain_swaps(n_atoms: int) -> list:
    """Build chain_swaps for one sample with two equal-length chains that can be swapped.

    Splits the atom range [0, n_atoms) into two halves (chain A and chain B)
    and returns identity + the A<->B swap. Each swap entry is
    (start1, end1, start2, end2, chainidx1, chainidx2).
    """
    half = n_atoms // 2
    identity = []
    swap_ab = [
        (0, half, half, 2 * half, 0, 1),
        (half, 2 * half, 0, half, 1, 0),
    ]
    return [identity, swap_ab]


_atom_keys = {"atom_pad_mask"}
_placements = get_feature_placements(atom_keys=_atom_keys, token_keys=set())
_placements_atom_features = _placements["atom_features"]
_placements_cp_atom_features = _placements["cp_atom_features"]

_placements_sample_coords = {"sample_coords": (Shard(0), Shard(1), Replicate())}
_placements_cp_sample_coords = {"sample_coords": (Shard(0), Replicate())}


def parallel_assert_minimum_lddt_symmetry_coords(rank, payload):
    """Test distributed minimum_lddt_symmetry_coords against serial Boltz-2 version.

    With DP sharding, each DP rank processes its own local batch of symmetry features.
    The coords DTensor is sharded along (DP, CP_0, CP_1) axes.
    """
    (
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        input_feats_global_host,
        feats_symmetry_global_host,
        expected_true_coords_per_sample,
        expected_true_mask_per_sample,
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
    inputs_atom["sample_coords"] = input_feats_global_host["sample_coords"].to(dtype=dtype)

    placements_cp = _placements_cp_atom_features | _placements_cp_sample_coords
    placements_dp_cp = _placements_atom_features | _placements_sample_coords

    feats_atom = distribute_atom_features(
        inputs_atom,
        placements_cp,
        placements_dp_cp,
        manager.device_mesh_subgroups,
        manager.group["cp"],
    )

    coords_dtensor = feats_atom["sample_coords"]
    atom_pad_mask_dtensor = feats_atom["atom_pad_mask"]

    coords_placements = coords_dtensor.placements

    num_dp_ranks = grid_group_sizes["dp"]
    global_batch_size = size_batch
    local_batch_size = global_batch_size // num_dp_ranks
    local_start = rank_dp * local_batch_size
    local_end = local_start + local_batch_size

    feats_local = {
        "all_coords": feats_symmetry_global_host["all_coords"][local_start:local_end].to(device),
        "all_resolved_mask": feats_symmetry_global_host["all_resolved_mask"][local_start:local_end].to(device),
        "crop_to_all_atom_map": feats_symmetry_global_host["crop_to_all_atom_map"][local_start:local_end].to(device),
        "chain_swaps": feats_symmetry_global_host["chain_swaps"][local_start:local_end],
        "amino_acids_symmetries": feats_symmetry_global_host["amino_acids_symmetries"][local_start:local_end],
        "ligand_symmetries": feats_symmetry_global_host["ligand_symmetries"][local_start:local_end],
        "atom_pad_mask": atom_pad_mask_dtensor,
    }

    for i_batch_local in range(local_batch_size):
        global_batch_idx = local_start + i_batch_local

        true_coords_dtensor, true_mask_dtensor = dtensor_minimum_lddt_symmetry_coords(
            coords=coords_dtensor,
            feats=feats_local,
            index_batch_local=i_batch_local,
            i_batch_multiplicity_local=i_batch_local,
        )

        assert (
            true_coords_dtensor.placements == coords_placements
        ), f"Sample {i_batch_local}: true_coords_dtensor.placements mismatch"
        assert (
            true_mask_dtensor.placements == coords_placements
        ), f"Sample {i_batch_local}: true_mask_dtensor.placements mismatch"

        expected_coords = expected_true_coords_per_sample[global_batch_idx]
        expected_mask = expected_true_mask_per_sample[global_batch_idx]

        coords_cp_gathered = true_coords_dtensor.redistribute(
            true_coords_dtensor.device_mesh,
            (true_coords_dtensor.placements[0], Replicate(), Replicate()),
        ).to_local()
        mask_cp_gathered = true_mask_dtensor.redistribute(
            true_mask_dtensor.device_mesh,
            (true_mask_dtensor.placements[0], Replicate(), Replicate()),
        ).to_local()

        atom_pad_mask_gathered = atom_pad_mask_dtensor.redistribute(
            atom_pad_mask_dtensor.device_mesh,
            (atom_pad_mask_dtensor.placements[0], Replicate(), Replicate()),
        ).to_local()

        real_atom_mask = atom_pad_mask_gathered[i_batch_local].bool()
        coords_no_pad = coords_cp_gathered[i_batch_local, real_atom_mask, :]
        mask_no_pad = mask_cp_gathered[i_batch_local, real_atom_mask]

        torch.testing.assert_close(
            coords_no_pad.cpu(),
            expected_coords,
            msg=f"Sample {i_batch_local} (global {global_batch_idx}): true_coords mismatch",
        )
        torch.testing.assert_close(
            mask_no_pad.cpu(),
            expected_mask.squeeze(0) if expected_mask.ndim > 1 else expected_mask,
            msg=f"Sample {i_batch_local} (global {global_batch_idx}): true_mask mismatch",
        )

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.skipif(not torch.cuda.is_available(), reason="CUDA not available")
@pytest.mark.parametrize("use_nontrivial_swaps", [False, True], ids=["identity_swaps", "nontrivial_swaps"])
def test_serial_minimum_lddt_symmetry_coords(use_nontrivial_swaps):
    """Verify serial minimum_lddt_symmetry_coords runs correctly on a single GPU.

    Tests both identity (no-op) and non-trivial (two-chain A<->B) swap cases.
    """
    batch_size = 2
    n_atoms = 24

    torch.manual_seed(123)
    sample_coords = torch.randn((batch_size, n_atoms, 3), dtype=torch.float32)
    all_coords = sample_coords.clone()
    all_resolved_mask = torch.ones((batch_size, n_atoms), dtype=torch.bool)
    crop_to_all_atom_map = torch.arange(n_atoms, dtype=torch.long).unsqueeze(0).expand(batch_size, -1).contiguous()

    chain_swaps_batch = None
    if use_nontrivial_swaps:
        chain_swaps_batch = [_make_two_chain_swaps(n_atoms) for _ in range(batch_size)]

    feats = _build_symmetry_features_for_batch(
        all_coords, all_resolved_mask, crop_to_all_atom_map, chain_swaps_batch=chain_swaps_batch
    )

    for i in range(batch_size):
        true_coords, true_mask = serial_minimum_lddt_symmetry_coords(
            coords=sample_coords[i : i + 1],
            feats=feats,
            index_batch=i,
        )
        assert true_coords.shape[-1] == 3, f"Sample {i}: unexpected coords shape {true_coords.shape}"
        assert true_mask.dtype == torch.bool, f"Sample {i}: unexpected mask dtype {true_mask.dtype}"
        assert true_mask.any(), f"Sample {i}: all-zero mask unexpected with all-resolved input"


@pytest.mark.skipif(not torch.cuda.is_available(), reason="CUDA not available")
@pytest.mark.parametrize("use_nontrivial_swaps", [False, True], ids=["identity_swaps", "nontrivial_swaps"])
@pytest.mark.parametrize(
    "setup_env",
    [
        ((1, (2, 2)), True, "cuda", "ENV"),  # dp=1, cp=(2,2), world_size=4
        ((2, (2, 2)), True, "cuda", "ENV"),  # dp=2, cp=(2,2), world_size=8
    ],
    indirect=("setup_env",),
)
def test_dtensor_minimum_lddt_symmetry_coords(setup_env, use_nontrivial_swaps):
    """Test distributed symmetry correction against Boltz-2 serial implementation.

    Parametrized for 4-GPU and 8-GPU configs, with both identity (no-op) and
    non-trivial (two-chain A<->B) chain swaps.
    """
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if torch.cuda.device_count() < world_size:
        pytest.skip(f"skip cuda test because torch.cuda.device_count() < {world_size}")

    size_ring = grid_group_sizes["cp"][0]
    num_dp_ranks = grid_group_sizes["dp"]
    batch_size_per_dp_rank = 1
    batch_size_global = batch_size_per_dp_rank * num_dp_ranks

    n_atoms_per_token = 3
    n_tokens = size_ring * 4
    n_atoms = n_atoms_per_token * n_tokens

    torch.manual_seed(42)
    feats_from_random = random_features(
        size_batch=batch_size_global,
        n_tokens=n_tokens,
        n_atoms=n_atoms,
        n_msa=1,
        atom_counts_per_token_range=(1, n_atoms_per_token),
        device=torch.device("cpu"),
        float_value_range=(-1.0, 1.0),
        selected_keys=["atom_pad_mask", "coords", "atom_counts_per_token"],
    )

    sample_coords_global = torch.randn((batch_size_global, n_atoms, 3), dtype=torch.float32)

    atom_pad_mask_global = feats_from_random["atom_pad_mask"]
    atom_counts_per_token_global = feats_from_random["atom_counts_per_token"]

    input_feats_global = {
        "atom_pad_mask": atom_pad_mask_global,
        "atom_counts_per_token": atom_counts_per_token_global,
        "sample_coords": sample_coords_global,
    }

    coords_for_symmetry = sample_coords_global
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

    expected_true_coords_per_sample = []
    expected_true_mask_per_sample = []

    for i in range(batch_size_global):
        expected_coords, expected_mask = serial_minimum_lddt_symmetry_coords(
            coords=coords_for_symmetry[i : i + 1],
            feats=feats_symmetry_global,
            index_batch=i,
        )
        expected_true_coords_per_sample.append(expected_coords.squeeze(0).detach().clone().cpu())
        expected_true_mask_per_sample.append(expected_mask.detach().clone().cpu())

    input_feats_global_host = {k: v.detach().clone().cpu() for k, v in input_feats_global.items()}

    payload = (
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        input_feats_global_host,
        feats_symmetry_global,
        expected_true_coords_per_sample,
        expected_true_mask_per_sample,
    )

    spawn_multiprocessing(parallel_assert_minimum_lddt_symmetry_coords, world_size, payload)
