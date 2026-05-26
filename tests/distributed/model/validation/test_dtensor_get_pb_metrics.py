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

"""Standalone test for DistributedValidator.get_pb_metrics at CP=(2,2).

Verifies that the distributed override produces identical results to the
serial compute_pb_geometry_metrics, compute_stereo_metrics, and
compute_pb_flatness_metrics by gathering DTensor features and comparing.
Tests both dp=1 (4 GPUs) and dp=2 (8 GPUs).

Test data includes a real PHE (phenylalanine) ligand with CCD-derived
geometry features (66 atom-pair edges, 1 chiral center, 1 aromatic 6-ring)
so that all three PB metric types produce non-trivial results.
"""

from __future__ import annotations

import pytest
import torch
from torch.distributed.tensor import distribute_tensor

from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.validation.rcsb import DistributedRCSBValidator
from boltz.distributed.model.validation.utils import gather_along_cp
from boltz.model.loss.inference import (
    compute_pb_flatness_metrics,
    compute_pb_geometry_metrics,
    compute_stereo_metrics,
)
from boltz.testing.utils import (
    LIGAND_KEYS,
    distribute_atom_features,
    get_feature_placements,
    make_pb_test_data,
    spawn_multiprocessing,
)

_ATOM_KEYS = {"atom_pad_mask", "atom_to_token"}
_TOKEN_KEYS = {"asym_id", "mol_type"}
_placements = get_feature_placements(atom_keys=_ATOM_KEYS, token_keys=_TOKEN_KEYS)
SINGLE_REPR = _placements["single"]

N_SAMPLES = 2


def _distribute_pb_data(batch_serial, out_serial, manager):
    """Distribute batch/out as DTensors for PB metrics.

    Atom features (atom_to_token, atom_pad_mask) and sample_atom_coords
    are distributed via ``distribute_atom_features`` with intersperse
    padding. Token features use ``distribute_tensor``.
    Non-sharded list features are sliced to the local DP rank's batch element.
    """
    mesh = manager.device_mesh_subgroups
    device = manager.device
    dp_rank = manager.group_rank["dp"]
    placements_token = _placements["token_features"]

    B = batch_serial["atom_pad_mask"].shape[0]
    coords = out_serial["sample_atom_coords"]
    coords_unflat = coords.unflatten(0, (B, N_SAMPLES))

    inputs_atom = {
        "atom_counts_per_token": batch_serial["atom_counts_per_token"],
        "atom_to_token": batch_serial["atom_to_token"],
        "atom_pad_mask": batch_serial["atom_pad_mask"],
    }
    for i in range(N_SAMPLES):
        inputs_atom[f"sample_atom_coords_{i}"] = coords_unflat[:, i]

    sample_cp = {f"sample_atom_coords_{i}": _placements["cp_single"] for i in range(N_SAMPLES)}
    sample_dp_cp = {f"sample_atom_coords_{i}": _placements["single"] for i in range(N_SAMPLES)}

    feats_atom = distribute_atom_features(
        inputs=inputs_atom,
        placements_cp=_placements["cp_atom_features"] | sample_cp,
        placements_dp_cp=_placements["atom_features"] | sample_dp_cp,
        device_mesh=mesh,
        cp_group=manager.group["cp"],
        multiplicities={"sample_atom_coords": N_SAMPLES},
    )

    batch_dt: dict = {
        "atom_to_token": feats_atom.pop("atom_to_token"),
        "atom_pad_mask": feats_atom.pop("atom_pad_mask"),
    }

    batch_dt["asym_id"] = distribute_tensor(
        batch_serial["asym_id"].to(device), device_mesh=mesh, placements=placements_token["asym_id"]
    )
    batch_dt["mol_type"] = distribute_tensor(
        batch_serial["mol_type"].to(device), device_mesh=mesh, placements=placements_token["mol_type"]
    )

    for k in LIGAND_KEYS:
        if k in batch_serial:
            batch_dt[k] = [batch_serial[k][dp_rank].to(device)]

    out_dt = {"sample_atom_coords": feats_atom.pop("sample_atom_coords")}
    return batch_dt, out_dt


def _parallel_pb_worker(rank, payload):
    (
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        serial_results_per_batch,
        batch_host,
        out_host,
    ) = payload

    mp = pytest.MonkeyPatch()
    if env_per_rank is not None:
        for k, v in env_per_rank.items():
            mp.setenv(k, f"{rank}" if v == "<INPUT_RANK>" else v)

    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    manager = DistributedManager()
    dp_rank = manager.group_rank["dp"]

    dist_validator = DistributedRCSBValidator(
        val_names=["RCSB"],
        confidence_prediction=False,
        physicalism_metrics=True,
    )
    dist_validator.to(manager.device)

    batch_dt, out_dt = _distribute_pb_data(batch_host, out_host, manager)

    batch_gathered = {
        "asym_id": gather_along_cp(batch_dt["asym_id"]),
        "mol_type": gather_along_cp(batch_dt["mol_type"]),
    }
    out_gathered = {
        "sample_atom_coords": gather_along_cp(out_dt["sample_atom_coords"]),
    }

    pb_failure_dict, pb_total_dict = dist_validator.get_pb_metrics(
        batch_dt,
        out_dt,
        batch_gathered,
        out_gathered,
    )

    expected_failure, expected_total = serial_results_per_batch[dp_rank]
    for key in expected_failure:
        torch.testing.assert_close(
            pb_failure_dict[key].cpu(),
            expected_failure[key],
            msg=f"PB failure mismatch on rank {rank} (dp={dp_rank}), key={key}",
        )
        torch.testing.assert_close(
            pb_total_dict[key].cpu(),
            expected_total[key],
            msg=f"PB total mismatch on rank {rank} (dp={dp_rank}), key={key}",
        )

    DistributedManager.cleanup()
    mp.undo()


@pytest.mark.skipif(not torch.cuda.is_available(), reason="CUDA not available")
@pytest.mark.parametrize(
    "setup_env",
    [
        ((1, (2, 2)), True, "cuda", "ENV"),
        ((2, (2, 2)), True, "cuda", "ENV"),
    ],
    indirect=("setup_env",),
    ids=["dp1-cp2x2", "dp2-cp2x2"],
)
def test_dtensor_get_pb_metrics(setup_env, canonical_mols_dir):
    """Distributed get_pb_metrics matches serial PB computation at CP=(2,2).

    Uses real PHE CCD geometry (66 edges, 1 chiral centre, 1 aromatic 6-ring)
    so that all three PB metric families produce non-trivial totals.
    """
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if torch.cuda.device_count() < world_size:
        pytest.skip(f"Need {world_size} GPUs, have {torch.cuda.device_count()}")

    cp0 = grid_group_sizes["cp"][0]
    num_dp = grid_group_sizes["dp"]
    n_tok = 8 * cp0
    n_atom = 4 * n_tok

    batch, out = make_pb_test_data(
        n_tok,
        n_atom,
        mols_dir=str(canonical_mols_dir),
        batch_size=num_dp,
        n_samples=N_SAMPLES,
        seed=42,
    )

    serial_results_per_batch = []
    for b in range(num_dp):
        b_feats: dict = {
            "atom_to_token": batch["atom_to_token"][b : b + 1],
            "atom_pad_mask": batch["atom_pad_mask"][b : b + 1],
            "asym_id": batch["asym_id"][b : b + 1],
            "mol_type": batch["mol_type"][b : b + 1],
        }
        b_feats.update({k: [batch[k][b]] for k in LIGAND_KEYS})

        b_coords = out["sample_atom_coords"][b * N_SAMPLES : (b + 1) * N_SAMPLES]
        (bl, ba, ic, nl) = compute_pb_geometry_metrics(pred_atom_coords=b_coords, feats=b_feats)
        (cav, ca, sbv, sb) = compute_stereo_metrics(pred_atom_coords=b_coords, feats=b_feats)
        (a5v, a5r, a6v, a6r, dbv, db) = compute_pb_flatness_metrics(pred_atom_coords=b_coords, feats=b_feats)

        assert nl.sum() > 0, "num_ligands should be > 0 (PHE present)"
        assert ca.sum() > 0, "num_chiral_atoms should be > 0 (PHE has 1 chiral center)"
        assert a6r.sum() > 0, "num_aromatic_6_rings should be > 0 (PHE has 1 aromatic ring)"

        failure = {
            "bond_length": bl.cpu(),
            "bond_angle": ba.cpu(),
            "internal_clash": ic.cpu(),
            "atom_chirality": cav.cpu(),
            "bond_stereochemistry": sbv.cpu(),
            "ring_5_flatness": a5v.cpu(),
            "ring_6_flatness": a6v.cpu(),
            "double_bond_flatness": dbv.cpu(),
        }
        total = {
            "bond_length": nl.cpu(),
            "bond_angle": nl.cpu(),
            "internal_clash": nl.cpu(),
            "atom_chirality": ca.cpu(),
            "bond_stereochemistry": sb.cpu(),
            "ring_5_flatness": a5r.cpu(),
            "ring_6_flatness": a6r.cpu(),
            "double_bond_flatness": db.cpu(),
        }
        serial_results_per_batch.append((failure, total))

    payload = (
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        serial_results_per_batch,
        batch,
        out,
    )
    spawn_multiprocessing(_parallel_pb_worker, world_size, payload)
