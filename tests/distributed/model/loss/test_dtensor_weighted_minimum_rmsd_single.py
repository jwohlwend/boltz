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

"""Tests for DTensor weighted_minimum_rmsd_single.

Verifies that the distributed (DTensor) version of weighted_minimum_rmsd_single
produces results matching the serial version from boltz.model.loss.validation.
"""

from __future__ import annotations

import pytest
import torch
from torch.distributed.tensor import Replicate, Shard, distribute_tensor

from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.loss.validation import (
    weighted_minimum_rmsd_single as dtensor_weighted_minimum_rmsd_single,
)
from boltz.model.loss.validation import weighted_minimum_rmsd_single as serial_weighted_minimum_rmsd_single
from boltz.testing.utils import distribute_atom_features, get_feature_placements, random_features, spawn_multiprocessing

_atom_keys = {"atom_pad_mask", "atom_to_token"}
_token_keys = {"mol_type"}
_placements = get_feature_placements(atom_keys=_atom_keys, token_keys=_token_keys)
_placements_atom_features = _placements["atom_features"]
_placements_cp_atom_features = _placements["cp_atom_features"]
_placements_token_features = _placements["token_features"]

_placements_coords = (Shard(0), Shard(1), Replicate())
_placements_cp_coords = (Shard(0), Replicate())


def parallel_assert_weighted_minimum_rmsd_single(rank, payload):
    """Test distributed weighted_minimum_rmsd_single against serial version."""
    (
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        input_feats_global_host,
        pred_coords_global_host,
        true_coords_global_host,
        expected_rmsd,
        expected_aligned,
        expected_weights,
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
    device_mesh = manager.device_mesh_subgroups

    inputs_atom = {
        k: v.to(dtype=dtype if v.dtype.is_floating_point else v.dtype)
        for k, v in input_feats_global_host.items()
        if k in _placements_cp_atom_features
    }
    inputs_atom["pred_coords"] = pred_coords_global_host.to(dtype=dtype)
    inputs_atom["true_coords"] = true_coords_global_host.to(dtype=dtype)

    placements_cp = _placements_cp_atom_features | {
        "pred_coords": _placements_cp_coords,
        "true_coords": _placements_cp_coords,
    }
    placements_dp_cp = _placements_atom_features | {
        "pred_coords": _placements_coords,
        "true_coords": _placements_coords,
    }

    feats_atom = distribute_atom_features(
        inputs_atom,
        placements_cp,
        placements_dp_cp,
        device_mesh,
        manager.group["cp"],
    )

    pred_coords_dtensor = feats_atom["pred_coords"]
    true_coords_dtensor = feats_atom["true_coords"]
    atom_pad_mask_dtensor = feats_atom["atom_pad_mask"]
    atom_to_token_dtensor = feats_atom["atom_to_token"]

    mol_type_dtensor = distribute_tensor(
        input_feats_global_host["mol_type"].to(device=device, dtype=dtype),
        device_mesh=device_mesh,
        placements=_placements_token_features["mol_type"],
    )

    atom_mask_float = atom_pad_mask_dtensor

    rmsd, aligned, weights = dtensor_weighted_minimum_rmsd_single(
        pred_atom_coords=pred_coords_dtensor,
        atom_coords=true_coords_dtensor,
        atom_mask=atom_mask_float,
        atom_to_token=atom_to_token_dtensor,
        mol_type=mol_type_dtensor,
    )

    rmsd_full = rmsd.full_tensor().cpu()
    aligned_full = aligned.full_tensor().cpu()

    atom_pad_mask_full = atom_pad_mask_dtensor.full_tensor().cpu()

    torch.testing.assert_close(
        rmsd_full,
        expected_rmsd,
        msg="RMSD mismatch between distributed and serial",
    )

    batch_size = pred_coords_global_host.shape[0]
    for b in range(batch_size):
        real_atom_mask = atom_pad_mask_full[b].bool()
        aligned_no_pad = aligned_full[b, real_atom_mask, :]
        expected_aligned_b = expected_aligned[b]

        torch.testing.assert_close(
            aligned_no_pad,
            expected_aligned_b,
            msg=f"Batch {b}: aligned coords mismatch",
        )

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.skipif(not torch.cuda.is_available(), reason="CUDA not available")
@pytest.mark.parametrize(
    "setup_env",
    [
        ((1, (2, 2)), True, "cuda", "ENV"),  # dp=1, cp=(2,2), world_size=4
        ((2, (2, 2)), True, "cuda", "ENV"),  # dp=2, cp=(2,2), world_size=8
    ],
    indirect=("setup_env",),
)
def test_dtensor_weighted_minimum_rmsd_single(setup_env):
    """Test distributed weighted_minimum_rmsd_single against serial version.

    Uses random features with proper atom-to-token mapping. Compares rmsd,
    aligned coordinates, and weights with the serial implementation.
    """
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if torch.cuda.device_count() < world_size:
        pytest.skip(f"skip cuda test because torch.cuda.device_count() < {world_size}")

    torch.manual_seed(42)

    size_ring = grid_group_sizes["cp"][0]
    num_dp_ranks = grid_group_sizes["dp"]
    batch_size = num_dp_ranks
    n_atoms_per_token = 3
    n_tokens = size_ring * 4
    n_atoms = n_atoms_per_token * n_tokens

    feats_from_random = random_features(
        size_batch=batch_size,
        n_tokens=n_tokens,
        n_atoms=n_atoms,
        n_msa=1,
        atom_counts_per_token_range=(1, n_atoms_per_token),
        device=torch.device("cpu"),
        float_value_range=(-1.0, 1.0),
        selected_keys=["atom_to_token", "mol_type", "atom_pad_mask", "atom_counts_per_token"],
    )

    pred_coords = torch.randn((batch_size, n_atoms, 3), dtype=torch.float32)
    true_coords = torch.randn((batch_size, n_atoms, 3), dtype=torch.float32)

    atom_to_token = feats_from_random["atom_to_token"].float()
    mol_type = feats_from_random["mol_type"]
    atom_pad_mask = feats_from_random["atom_pad_mask"].float()

    serial_device = torch.device("cuda:0")
    expected_rmsd, expected_aligned, expected_weights = serial_weighted_minimum_rmsd_single(
        pred_atom_coords=pred_coords.to(serial_device),
        atom_coords=true_coords.to(serial_device),
        atom_mask=atom_pad_mask.to(serial_device),
        atom_to_token=atom_to_token.to(serial_device),
        mol_type=mol_type.to(serial_device).float(),
    )
    expected_rmsd = expected_rmsd.detach().cpu()
    expected_aligned = expected_aligned.detach().cpu()
    expected_weights = expected_weights.detach().cpu()

    input_feats_global = {
        "atom_to_token": atom_to_token,
        "mol_type": mol_type,
        "atom_pad_mask": atom_pad_mask,
        "atom_counts_per_token": feats_from_random["atom_counts_per_token"],
    }
    input_feats_global_host = {k: v.detach().clone().cpu() for k, v in input_feats_global.items()}

    payload = (
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        input_feats_global_host,
        pred_coords.detach().clone().cpu(),
        true_coords.detach().clone().cpu(),
        expected_rmsd,
        expected_aligned,
        expected_weights,
    )

    spawn_multiprocessing(parallel_assert_weighted_minimum_rmsd_single, world_size, payload)
