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

"""Standalone tests for DistributedValidator.get_clash_metrics at CP=(2,2).

Tests:
  - test_clash_score_counts_and_fraction: Triton clash_score kernel sanity check.
  - test_dtensor_get_clash_metrics: Distributed get_clash_metrics vs serial
    compute_chain_clashes, verifying the full DTensor gather path at CP=(2,2)
    with both dp=1 (4 GPUs) and dp=2 (8 GPUs).
"""

from __future__ import annotations

import pytest
import torch
from torch.distributed.tensor import distribute_tensor

from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.loss.validation import clash_score
from boltz.distributed.model.validation.rcsb import DistributedRCSBValidator
from boltz.distributed.model.validation.utils import gather_along_cp
from boltz.model.loss.inference import compute_chain_clashes
from boltz.testing.utils import distribute_atom_features, get_feature_placements, random_features, spawn_multiprocessing

_ATOM_KEYS = {"atom_pad_mask", "atom_to_token", "ref_element"}
_TOKEN_KEYS = {"asym_id"}
_placements = get_feature_placements(atom_keys=_ATOM_KEYS, token_keys=_TOKEN_KEYS)
SINGLE_REPR = _placements["single"]

N_SAMPLES = 2


# ---------------------------------------------------------------------------
# clash_score Triton kernel test
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not torch.cuda.is_available(), reason="CUDA not available")
def test_clash_score_counts_and_fraction():
    """Verify clash_score Triton kernel counts and fraction computation."""
    device = torch.device("cuda")
    clash_cutoff = 2.0
    multiplicity = 2

    coords_repr = torch.tensor(
        [
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [5.0, 0.0, 0.0], [9.0, 0.0, 0.0]],
            [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [6.0, 0.0, 0.0], [9.0, 0.0, 0.0]],
        ],
        device=device,
        dtype=torch.float32,
    )
    token_pad_mask = torch.tensor([[True, True, True, False]], device=device)

    clash_atoms_count, clash_atoms_fraction = clash_score(
        coords_repr=coords_repr,
        token_pad_mask=token_pad_mask,
        multiplicity=multiplicity,
        clash_cutoff=clash_cutoff,
    )

    expected_count = torch.tensor([2, 0], device=device, dtype=clash_atoms_count.dtype)
    expected_fraction = torch.tensor([2.0 / 3.0, 0.0], device=device, dtype=clash_atoms_fraction.dtype)

    torch.testing.assert_close(clash_atoms_count, expected_count)
    torch.testing.assert_close(clash_atoms_fraction, expected_fraction)


# ---------------------------------------------------------------------------
# Distributed get_clash_metrics test
# ---------------------------------------------------------------------------


def _make_clash_test_data(n_tok, n_atom, batch_size=1, seed=42):
    """Create test data for clash metrics."""
    atoms_per_tok = n_atom // n_tok
    rng = torch.Generator().manual_seed(seed)

    feats = random_features(
        size_batch=batch_size,
        n_tokens=n_tok,
        n_atoms=n_atom,
        n_msa=1,
        atom_counts_per_token_range=(atoms_per_tok, atoms_per_tok),
        device=torch.device("cpu"),
        float_value_range=(-1.0, 1.0),
        selected_keys=["atom_to_token", "atom_pad_mask", "atom_counts_per_token", "asym_id", "ref_element"],
        rng=rng,
    )

    batch: dict = {}
    for k, v in feats.items():
        batch[k] = v.to(torch.float32) if v.is_floating_point() else v
    batch["atom_to_token"] = batch["atom_to_token"].to(torch.float32)

    sample_coords = torch.randn(batch_size * N_SAMPLES, n_atom, 3, generator=rng)

    edges = []
    for t in range(n_tok):
        for a in range(atoms_per_tok - 1):
            edges.append([t * atoms_per_tok + a, t * atoms_per_tok + a + 1])
        if t < n_tok - 1:
            edges.append([t * atoms_per_tok + atoms_per_tok - 1, (t + 1) * atoms_per_tok])
    edge_tensor = torch.tensor(edges, dtype=torch.long).T if edges else torch.empty(2, 0, dtype=torch.long)
    batch["connections_edge_index"] = [edge_tensor.clone() for _ in range(batch_size)]

    batch["chain_symmetries"] = []
    for b in range(batch_size):
        unique_asym = batch["asym_id"][b].unique().tolist()
        batch["chain_symmetries"].append([[(int(cid), 0, "A", "PROTEIN", 0)] for cid in unique_asym])

    out = {"sample_atom_coords": sample_coords}
    return batch, out


def _distribute_clash_data(batch_serial, out_serial, manager):
    """Distribute batch/out as DTensors for clash metrics.

    Atom features (atom_to_token, atom_pad_mask, ref_element) and
    sample_atom_coords are distributed via ``distribute_atom_features``
    with intersperse padding. Token features use ``distribute_tensor``.
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
        "ref_element": batch_serial["ref_element"],
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
        "ref_element": feats_atom.pop("ref_element"),
    }

    batch_dt["asym_id"] = distribute_tensor(
        batch_serial["asym_id"].to(device), device_mesh=mesh, placements=placements_token["asym_id"]
    )
    batch_dt["chain_symmetries"] = [batch_serial["chain_symmetries"][dp_rank]]
    batch_dt["connections_edge_index"] = [batch_serial["connections_edge_index"][dp_rank].to(device)]

    out_dt = {"sample_atom_coords": feats_atom.pop("sample_atom_coords")}
    return batch_dt, out_dt


def _parallel_clash_worker(rank, payload):
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

    batch_dt, out_dt = _distribute_clash_data(batch_host, out_host, manager)

    batch_gathered = {
        "asym_id": gather_along_cp(batch_dt["asym_id"]),
        "atom_pad_mask": gather_along_cp(batch_dt["atom_pad_mask"]),
    }
    out_gathered = {
        "sample_atom_coords": gather_along_cp(out_dt["sample_atom_coords"]),
    }

    pair_clash_dict, pair_total_dict = dist_validator.get_clash_metrics(
        batch_dt,
        out_dt,
        batch_gathered,
        out_gathered,
    )

    expected_clash, expected_total = serial_results_per_batch[dp_rank]
    for key in expected_clash:
        torch.testing.assert_close(
            pair_clash_dict[key].cpu(),
            expected_clash[key],
            msg=f"Clash mismatch on rank {rank} (dp={dp_rank}), key={key}",
        )
        torch.testing.assert_close(
            pair_total_dict[key].cpu(),
            expected_total[key],
            msg=f"Clash total mismatch on rank {rank} (dp={dp_rank}), key={key}",
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
def test_dtensor_get_clash_metrics(setup_env):
    """Distributed get_clash_metrics matches serial compute_chain_clashes at CP=(2,2)."""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if torch.cuda.device_count() < world_size:
        pytest.skip(f"Need {world_size} GPUs, have {torch.cuda.device_count()}")

    cp0 = grid_group_sizes["cp"][0]
    num_dp = grid_group_sizes["dp"]
    n_tok = 8 * cp0
    n_atom = 4 * n_tok

    batch, out = _make_clash_test_data(n_tok, n_atom, batch_size=num_dp, seed=42)

    serial_results_per_batch = []
    for b in range(num_dp):
        b_feats = {
            "atom_to_token": batch["atom_to_token"][b : b + 1],
            "atom_pad_mask": batch["atom_pad_mask"][b : b + 1],
            "ref_element": batch["ref_element"][b : b + 1],
            "asym_id": batch["asym_id"][b : b + 1],
            "chain_symmetries": [batch["chain_symmetries"][b]],
            "connections_edge_index": [batch["connections_edge_index"][b]],
        }
        b_coords = out["sample_atom_coords"][b * N_SAMPLES : (b + 1) * N_SAMPLES]
        cd, td = compute_chain_clashes(pred_atom_coords=b_coords, feats=b_feats)
        serial_results_per_batch.append(({k: v.cpu() for k, v in cd.items()}, {k: v.cpu() for k, v in td.items()}))

    payload = (
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        serial_results_per_batch,
        batch,
        out,
    )
    spawn_multiprocessing(_parallel_clash_worker, world_size, payload)
