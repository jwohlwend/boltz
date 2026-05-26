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

"""End-to-end integration tests for DistributedRCSBValidator.

Tests that the distributed validation pipeline produces identical metric
values to the serial RCSBValidator at CP=2 and CP=4.
"""

from __future__ import annotations

from types import SimpleNamespace

import pytest
import torch
from torch import Tensor
from torch.distributed.tensor import DTensor, distribute_tensor

from boltz.data import const
from boltz.distributed.comm import TransposeComm
from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.validation.rcsb import DistributedRCSBValidator
from boltz.model.validation.rcsb import RCSBValidator
from boltz.testing.utils import (
    get_feature_placements,
    random_features,
    skip_if_cuda_not_avail_or_device_count_less_than_word_size,
    spawn_multiprocessing,
)


class _ConfigNamespace(SimpleNamespace):
    """SimpleNamespace with dict-like .get() for OmegaConf compat."""

    def get(self, key, default=None):
        return getattr(self, key, default)


_ATOM_KEYS = {"atom_pad_mask", "atom_to_token", "atom_resolved_mask", "coords", "ref_element", "token_to_rep_atom"}
_TOKEN_KEYS = {"mol_type", "asym_id", "token_pad_mask", "token_index", "token_disto_mask"}
_placements = get_feature_placements(atom_keys=_ATOM_KEYS, token_keys=_TOKEN_KEYS)

SINGLE_REPR = _placements["single"]
PAIR_REPR = _placements["pair"]
ENSEMBLE_REPR = _placements["atom_features"]["coords"]


def _make_mock_model(device="cpu", confidence_prediction=False, diffusion_samples=1, symmetry_correction=True):
    """Build a lightweight mock LightningModule for validator tests."""
    model = SimpleNamespace(
        val_group_mapper={0: {"label": "RCSB", "symmetry_correction": symmetry_correction}},
        validation_args=_ConfigNamespace(
            recycling_steps=0,
            sampling_steps=1,
            diffusion_samples=diffusion_samples,
        ),
        confidence_prediction=confidence_prediction,
        aggregate_distogram=True,
        min_dist=2.0,
        max_dist=22.0,
        num_bins=8,
        num_distograms=1,
        device=device,
        dp_group=None,
        token_level_confidence=True,
        _logged={},
    )

    def _log(name, value, **kwargs):
        model._logged[name] = value

    model.log = _log

    def _get_true_coordinates(
        batch,
        out,
        diffusion_samples,
        symmetry_correction,
        expand_to_diffusion_samples=True,
    ):
        K, L = batch["coords"].shape[1:3]
        mask = batch["atom_resolved_mask"]
        tc = batch["coords"].squeeze(0)
        if expand_to_diffusion_samples:
            tc = tc.repeat((diffusion_samples, 1, 1)).reshape(
                diffusion_samples,
                K,
                L,
                3,
            )
            mask = mask.repeat_interleave(diffusion_samples, dim=0)
        else:
            mask = mask.squeeze(0)
        return {
            "true_coords": tc,
            "true_coords_resolved_mask": mask,
            "rmsds": 0,
            "best_rmsd_recall": 0,
            "best_rmsd_precision": 0,
        }

    model.get_true_coordinates = _get_true_coordinates
    return model


def _make_test_data(
    n_tok,
    n_atom,
    num_bins,
    seed=42,
    n_samples=1,
    confidence=False,
    physicalism=False,
):
    """Generate deterministic batch and output dicts for validator tests."""
    B, K, D = 1, 1, 1
    atoms_per_tok = n_atom // n_tok
    assert n_atom == n_tok * atoms_per_tok

    rng = torch.Generator().manual_seed(seed)

    selected_keys = [
        "atom_pad_mask",
        "atom_to_token",
        "atom_resolved_mask",
        "coords",
        "mol_type",
        "asym_id",
        "token_pad_mask",
        "token_index",
        "token_disto_mask",
        "disto_target",
        "contact_conditioning",
    ]
    if confidence:
        selected_keys += ["token_to_rep_atom", "r_set_to_rep_atom", "frames_idx"]
    if physicalism:
        selected_keys += ["ref_element"]

    feats = random_features(
        size_batch=B,
        n_tokens=n_tok,
        n_atoms=n_atom,
        n_msa=1,
        atom_counts_per_token_range=(atoms_per_tok, atoms_per_tok),
        device=torch.device("cpu"),
        float_value_range=(-1.0, 1.0),
        selected_keys=selected_keys,
        num_disto_bins=num_bins,
        rng=rng,
    )

    batch = {}
    for k, v in feats.items():
        if v.is_floating_point():
            batch[k] = v.to(torch.float32)
        else:
            batch[k] = v

    # atom_to_token must be float for einsum in factored_lddt_loss
    batch["atom_to_token"] = batch["atom_to_token"].to(torch.float32)

    # disto_target: (B, N, N, bins) -> (B, N, N, K, bins) for v2 distogram loss
    batch["disto_target"] = batch["disto_target"].unsqueeze(3)

    # disto_coords_ensemble: needed by serial and distributed compute_disto_lddt
    batch["disto_coords_ensemble"] = torch.randn(B, K * n_tok, 3, generator=rng)

    batch["idx_dataset"] = torch.tensor([0])

    if n_samples == 1:
        sample_coords = batch["coords"][:, 0, :, :].clone()
    else:
        sample_coords = torch.randn(n_samples, n_atom, 3, generator=rng)

    pdisto = torch.randn(B, n_tok, n_tok, D, num_bins, generator=rng)
    out = {"sample_atom_coords": sample_coords, "pdistogram": pdisto}

    if confidence:
        # Pad r_set_to_rep_atom to (B, n_tok, n_atom) for diagonal sharding
        r_set = batch["r_set_to_rep_atom"]
        n_r = r_set.shape[1]
        if n_r < n_tok:
            pad = torch.zeros(B, n_tok - n_r, n_atom, dtype=r_set.dtype, device=r_set.device)
            batch["r_set_to_rep_atom"] = torch.cat([r_set, pad], dim=1)

        batch["frame_resolved_mask"] = torch.ones(B, n_tok, dtype=torch.bool)

        out["plddt"] = torch.randn(n_samples, n_tok, generator=rng).sigmoid()
        out["pde"] = torch.randn(n_samples, n_tok, n_tok, generator=rng)
        out["pae"] = torch.randn(n_samples, n_tok, n_tok, generator=rng)
        for key in (
            "complex_plddt",
            "complex_iplddt",
            "complex_pde",
            "complex_ipde",
            "ptm",
            "iptm",
            "ligand_iptm",
            "protein_iptm",
        ):
            out[key] = torch.randn(n_samples, generator=rng)

    if physicalism:
        edges = []
        for t in range(n_tok):
            for a in range(atoms_per_tok - 1):
                edges.append([t * atoms_per_tok + a, t * atoms_per_tok + a + 1])
            if t < n_tok - 1:
                edges.append([t * atoms_per_tok + atoms_per_tok - 1, (t + 1) * atoms_per_tok])
        if edges:
            batch["connections_edge_index"] = [torch.tensor(edges, dtype=torch.long).T]
        else:
            batch["connections_edge_index"] = [torch.empty(2, 0, dtype=torch.long)]

        batch["chain_symmetries"] = [[[(0, 0, "A", "PROTEIN", 0)]]]

        for key in (
            "ligand_edge_index",
            "ligand_stereo_bond_index",
            "ligand_planar_double_bond_index",
        ):
            batch[key] = [torch.empty(2, 0, dtype=torch.long)]
        for key in (
            "ligand_edge_lower_bounds",
            "ligand_edge_upper_bounds",
        ):
            batch[key] = [torch.empty(0)]
        for key in (
            "ligand_edge_bond_mask",
            "ligand_edge_angle_mask",
            "ligand_chiral_check_mask",
            "ligand_chiral_atom_orientations",
            "ligand_stereo_check_mask",
            "ligand_stereo_bond_orientations",
        ):
            batch[key] = [torch.empty(0, dtype=torch.bool)]
        batch["ligand_chiral_atom_index"] = [torch.empty(4, 0, dtype=torch.long)]
        batch["ligand_aromatic_5_ring_index"] = [torch.empty(5, 0, dtype=torch.long)]
        batch["ligand_aromatic_6_ring_index"] = [torch.empty(6, 0, dtype=torch.long)]

    return batch, out


def _extract_all_metrics(validator):
    """Collect all non-NaN metric values from a validator into a flat dict."""
    results = {}
    for gname, mlist in validator.folding_metrics.items():
        for ds_idx in range(len(mlist)):
            for mname, metric in mlist[ds_idx].items():
                val = metric.compute()
                if not torch.isnan(val):
                    results[f"{gname}/{ds_idx}/{mname}"] = val.item()
    if validator.physicalism_metrics:
        for gname, mlist in validator.physicalism_metrics.items():
            for ds_idx in range(len(mlist)):
                for mname, metric in mlist[ds_idx].items():
                    val = metric.compute()
                    if not torch.isnan(val):
                        results[f"phys/{gname}/{ds_idx}/{mname}"] = val.item()
    if hasattr(validator, "confidence_metrics"):
        for gname, mlist in validator.confidence_metrics.items():
            for ds_idx in range(len(mlist)):
                for mname, metric in mlist[ds_idx].items():
                    val = metric.compute()
                    if not torch.isnan(val):
                        results[f"conf/{gname}/{ds_idx}/{mname}"] = val.item()
    return results


def _parallel_distributed_test(rank, payload):
    (
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        serial_metrics,
        batch_host,
        out_host,
        confidence_prediction,
        physicalism_metrics,
        expand_to_diffusion_samples,
    ) = payload

    n_samples = 2 if confidence_prediction else 1

    mp = pytest.MonkeyPatch()
    if env_per_rank is not None:
        for k, v in env_per_rank.items():
            mp.setenv(k, f"{rank}" if v == "<INPUT_RANK>" else v)

    DistributedManager.initialize(
        grid_group_sizes,
        device_type=device_type,
        backend=backend,
    )
    manager = DistributedManager()
    comm = TransposeComm(manager.group["cp"], manager.layout_subgroups["cp"])

    # --- distribute batch/out as DTensors (inlined) ---
    mesh = manager.device_mesh_subgroups
    device = manager.device
    cp0_size = mesh.get_group("cp_axis_0").size()
    cp0_rank = mesh.get_coordinate()[1]

    placements_token = _placements["token_features"]
    placements_atom = _placements["atom_features"]

    def _dt(tensor, placements):
        return distribute_tensor(
            tensor.to(device),
            device_mesh=mesh,
            placements=placements,
        )

    DIAG_SHARDED_KEYS = {"atom_to_token", "token_to_rep_atom", "r_set_to_rep_atom"}

    batch_dtensor = {}
    for key, val in batch_host.items():
        if key == "idx_dataset":
            batch_dtensor[key] = val.to(device)
            continue
        if key in DIAG_SHARDED_KEYS:
            if key == "r_set_to_rep_atom" and cp0_size > 1:
                # Reorder r_set rows so each shard's row range only contains
                # elements whose representative atom falls within that shard's
                # column range.
                B, nr, nc = val.shape
                nr_l = nr // cp0_size
                nc_l = nc // cp0_size
                valid_mask = val[0].any(dim=-1)
                atom_indices = val[0].argmax(dim=-1)
                shard_of_row = atom_indices // nc_l
                reordered = torch.zeros_like(val)
                for s in range(cp0_size):
                    shard_rows = val[:, (shard_of_row == s) & valid_mask]
                    n = min(shard_rows.shape[1], nr_l)
                    reordered[:, s * nr_l : s * nr_l + n] = shard_rows[:, :n]
                val = reordered

            B, nr, nc = val.shape
            nr_l = nr // cp0_size
            nc_l = nc // cp0_size
            local_block = (
                val[
                    :,
                    cp0_rank * nr_l : (cp0_rank + 1) * nr_l,
                    cp0_rank * nc_l : (cp0_rank + 1) * nc_l,
                ]
                .contiguous()
                .to(device)
            )
            batch_dtensor[key] = DTensor.from_local(
                local_block,
                device_mesh=mesh,
                placements=SINGLE_REPR,
            )
            continue
        if key in placements_atom:
            batch_dtensor[key] = _dt(val, placements_atom[key])
            continue
        if key in placements_token:
            batch_dtensor[key] = _dt(val, placements_token[key])
            continue
        if key in ("disto_target", "contact_conditioning"):
            batch_dtensor[key] = _dt(val, PAIR_REPR)
            continue
        if isinstance(val, Tensor) and val.ndim >= 2:
            batch_dtensor[key] = _dt(val, SINGLE_REPR)
        elif isinstance(val, Tensor):
            batch_dtensor[key] = val.to(device)
        else:
            batch_dtensor[key] = val

    out_dtensor = {
        "sample_atom_coords": _dt(out_host["sample_atom_coords"], SINGLE_REPR),
        "pdistogram": _dt(out_host["pdistogram"], PAIR_REPR),
    }
    for key in ("plddt",):
        if key in out_host:
            out_dtensor[key] = _dt(out_host[key], SINGLE_REPR)
    for key in ("pde", "pae"):
        if key in out_host:
            out_dtensor[key] = _dt(out_host[key], PAIR_REPR)
    for key in (
        "complex_plddt",
        "complex_iplddt",
        "complex_pde",
        "complex_ipde",
        "ptm",
        "iptm",
        "ligand_iptm",
        "protein_iptm",
    ):
        if key in out_host:
            out_dtensor[key] = out_host[key].to(device)
    # --- end distribute ---

    dist_validator = DistributedRCSBValidator(
        val_names=["RCSB"],
        confidence_prediction=confidence_prediction,
        physicalism_metrics=physicalism_metrics,
        rmsd_metrics=True,
        clash_score_metrics=True,
    )
    dist_validator.to(manager.device)

    dist_model = _make_mock_model(
        device=str(manager.device),
        confidence_prediction=confidence_prediction,
        diffusion_samples=n_samples,
    )
    dist_model.dp_group = manager.device_mesh_subgroups.get_group(0)

    dist_validator.common_val_step(
        model=dist_model,
        batch=batch_dtensor,
        out=out_dtensor,
        idx_dataset=0,
        expand_to_diffusion_samples=expand_to_diffusion_samples,
        transpose_comm=comm,
    )

    dist_metrics = _extract_all_metrics(dist_validator)

    # Metrics that only exist in the distributed validator (not in serial).
    _DISTRIBUTED_ONLY_PREFIXES = ("rmsd/", "clash_score/")

    for key in serial_metrics:
        assert key in dist_metrics, f"Metric {key!r} missing in distributed result (rank {rank})"
        torch.testing.assert_close(
            torch.tensor(dist_metrics[key]),
            torch.tensor(serial_metrics[key]),
            msg=lambda m: f"Metric {key!r} mismatch on rank {rank}: {m}",
        )

    for key in dist_metrics:
        if any(key.startswith(p) for p in _DISTRIBUTED_ONLY_PREFIXES):
            continue
        assert key in serial_metrics, f"Unexpected metric {key!r} in distributed (rank {rank})"

    DistributedManager.cleanup()
    mp.undo()


def _parallel_epoch_end_test(rank, payload):
    (
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        confidence_prediction,
        physicalism_metrics,
    ) = payload

    mp = pytest.MonkeyPatch()
    if env_per_rank is not None:
        for k, v in env_per_rank.items():
            mp.setenv(k, f"{rank}" if v == "<INPUT_RANK>" else v)

    DistributedManager.initialize(
        grid_group_sizes,
        device_type=device_type,
        backend=backend,
    )
    manager = DistributedManager()

    dv = DistributedRCSBValidator(
        val_names=["RCSB"],
        confidence_prediction=confidence_prediction,
        physicalism_metrics=physicalism_metrics,
        rmsd_metrics=True,
        clash_score_metrics=True,
    )
    dv.to(manager.device)

    model = _make_mock_model(
        device=str(manager.device),
        confidence_prediction=confidence_prediction,
    )
    model.dp_group = manager.device_mesh_subgroups.get_group(0)

    logged_kwargs = []
    orig_log = model.log

    def _capture(*args, **kwargs):
        logged_kwargs.append(dict(kwargs))
        return orig_log(*args, **kwargs)

    model.log = _capture
    dv.on_epoch_end(model=model)

    DistributedManager.cleanup()
    mp.undo()


@pytest.mark.parametrize(
    "setup_env, confidence_prediction, physicalism_metrics",
    [
        (((1, (2, 2)), True, "cuda", "ENV"), True, True),
    ],
    indirect=("setup_env",),
    ids=[
        "dp1-cp2x2-conf-phys",
    ],
)
@pytest.mark.parametrize("expand_to_diffusion_samples", [True, False])
def test_distributed_rcsb_validator_matches_serial(
    setup_env, confidence_prediction, physicalism_metrics, expand_to_diffusion_samples
):
    """Distributed validator metrics match serial."""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    skip_if_cuda_not_avail_or_device_count_less_than_word_size(
        device_type=device_type,
        world_size=world_size,
    )

    cp_axis_0 = grid_group_sizes["cp"][0]
    n_tok = 8 * cp_axis_0
    n_atom = 32 * cp_axis_0
    num_bins = 8

    # --- run serial reference (inlined) ---
    n_samples = 2 if confidence_prediction else 1
    batch, out = _make_test_data(
        n_tok,
        n_atom,
        num_bins,
        seed=42,
        n_samples=n_samples,
        confidence=confidence_prediction,
        physicalism=physicalism_metrics,
    )
    model = _make_mock_model(
        device="cpu",
        confidence_prediction=confidence_prediction,
        diffusion_samples=n_samples,
    )
    validator = RCSBValidator(
        val_names=["RCSB"],
        confidence_prediction=confidence_prediction,
        physicalism_metrics=physicalism_metrics,
    )
    validator.common_val_step(
        model=model,
        batch=batch,
        out=out,
        idx_dataset=0,
        expand_to_diffusion_samples=expand_to_diffusion_samples,
    )
    serial_metrics = _extract_all_metrics(validator)
    # --- end serial reference ---

    batch_host, out_host = _make_test_data(
        n_tok,
        n_atom,
        num_bins,
        seed=42,
        n_samples=n_samples,
        confidence=confidence_prediction,
        physicalism=physicalism_metrics,
    )

    payload = (
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        serial_metrics,
        batch_host,
        out_host,
        confidence_prediction,
        physicalism_metrics,
        expand_to_diffusion_samples,
    )
    spawn_multiprocessing(_parallel_distributed_test, world_size, payload)


@pytest.mark.parametrize(
    "setup_env, confidence_prediction, physicalism_metrics",
    [
        (((1, (2, 2)), True, "cpu", "ENV"), True, True),
        (((2, (1, 1)), True, "cuda", "ENV"), True, True),
    ],
    indirect=("setup_env",),
    ids=[
        "cpu-cp2x2-conf-phys",
        "cuda-dp2-cp1x1-conf-phys",
    ],
)
def test_distributed_rcsb_validator_epoch_end_sync(setup_env, confidence_prediction, physicalism_metrics):
    """common_on_epoch_end wraps model.log with sync_dist_group."""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    skip_if_cuda_not_avail_or_device_count_less_than_word_size(
        device_type=device_type,
        world_size=world_size,
    )

    payload = (
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        confidence_prediction,
        physicalism_metrics,
    )
    spawn_multiprocessing(_parallel_epoch_end_test, world_size, payload)


def _parallel_dp_all_reduce_test(rank, payload):
    """Worker for test_dp_all_reduce_epoch_end.

    Each DP rank updates the disto_loss MeanMetric with a rank-specific
    value.  After common_on_epoch_end (which calls _dp_all_reduce_metrics),
    the logged disto_loss should be the global weighted mean across DP ranks.
    """
    grid_group_sizes, device_type, backend, env_per_rank = payload

    mp = pytest.MonkeyPatch()
    if env_per_rank is not None:
        for k, v in env_per_rank.items():
            mp.setenv(k, f"{rank}" if v == "<INPUT_RANK>" else v)

    DistributedManager.initialize(
        grid_group_sizes,
        device_type=device_type,
        backend=backend,
    )
    manager = DistributedManager()
    dp_rank = manager.device_mesh_subgroups.get_coordinate()[0]

    dv = DistributedRCSBValidator(
        val_names=["RCSB"],
        confidence_prediction=True,
        physicalism_metrics=True,
        rmsd_metrics=True,
        clash_score_metrics=True,
    )
    dv.to(manager.device)

    rank_values = [10.0, 20.0]
    rank_weights = [1.0, 3.0]
    local_val = torch.tensor(rank_values[dp_rank], device=manager.device)
    local_weight = torch.tensor(rank_weights[dp_rank], device=manager.device)
    dv.folding_metrics["disto_loss"][0]["disto_loss"].update(local_val, local_weight)

    model = _make_mock_model(device=str(manager.device))
    model.dp_group = manager.device_mesh_subgroups.get_group(0)

    dv.on_epoch_end(model=model)

    expected_global_mean = sum(v * w for v, w in zip(rank_values, rank_weights)) / sum(rank_weights)

    logged_disto = model._logged.get("val/disto_loss")
    assert logged_disto is not None, f"Rank {rank}: val/disto_loss not logged"

    logged_val = logged_disto.item() if isinstance(logged_disto, Tensor) else float(logged_disto)
    torch.testing.assert_close(
        torch.tensor(logged_val),
        torch.tensor(expected_global_mean),
        msg=lambda m: f"Rank {rank}: DP all-reduce disto_loss mismatch: {m}",
    )

    DistributedManager.cleanup()
    mp.undo()


@pytest.mark.parametrize(
    "setup_env",
    [
        ((2, (1, 1)), True, "cuda", "ENV"),
    ],
    indirect=("setup_env",),
    ids=["cuda-dp2-cp1x1"],
)
def test_dp_all_reduce_epoch_end(setup_env):
    """DP all-reduce in common_on_epoch_end produces correct global means."""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    skip_if_cuda_not_avail_or_device_count_less_than_word_size(
        device_type=device_type,
        world_size=world_size,
    )

    payload = (grid_group_sizes, device_type, backend, env_per_rank)
    spawn_multiprocessing(_parallel_dp_all_reduce_test, world_size, payload)


def test_distributed_rcsb_validator_mro():
    """MRO gives DistributedValidator precedence over RCSBValidator."""
    from boltz.distributed.model.validation.validator import DistributedValidator
    from boltz.model.validation.validator import Validator

    mro = DistributedRCSBValidator.__mro__
    dv_idx = mro.index(DistributedValidator)
    rcsb_idx = mro.index(RCSBValidator)
    v_idx = mro.index(Validator)

    assert dv_idx < rcsb_idx < v_idx, f"MRO wrong: DV@{dv_idx}, RCSB@{rcsb_idx}, V@{v_idx}"

    assert DistributedRCSBValidator.common_val_step is DistributedValidator.common_val_step
    assert DistributedRCSBValidator.common_on_epoch_end is DistributedValidator.common_on_epoch_end
    assert DistributedRCSBValidator._dp_all_reduce_metrics is DistributedValidator._dp_all_reduce_metrics
    assert DistributedRCSBValidator.on_epoch_end is RCSBValidator.on_epoch_end


def _parallel_epoch_end_logged_metric_names_test(rank, payload):
    (
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        confidence_prediction,
        physicalism_metrics,
        rmsd_metrics,
        clash_score_metrics,
    ) = payload

    mp = pytest.MonkeyPatch()
    if env_per_rank is not None:
        for k, v in env_per_rank.items():
            mp.setenv(k, f"{rank}" if v == "<INPUT_RANK>" else v)

    DistributedManager.initialize(
        grid_group_sizes,
        device_type=device_type,
        backend=backend,
    )
    manager = DistributedManager()

    dv = DistributedRCSBValidator(
        val_names=["RCSB"],
        confidence_prediction=confidence_prediction,
        physicalism_metrics=physicalism_metrics,
        rmsd_metrics=rmsd_metrics,
        clash_score_metrics=clash_score_metrics,
    )
    dv.to(manager.device)

    model = _make_mock_model(
        device=str(manager.device),
        confidence_prediction=confidence_prediction,
    )
    model.dp_group = manager.device_mesh_subgroups.get_group(0)

    logged_names = []

    def _capture_log(name, value, **kwargs):
        logged_names.append(name)

    model.log = _capture_log
    dv.on_epoch_end(model=model)

    # --- build expected logged names (inlined) ---
    expected = set()

    for m_ in [*const.out_types, "pocket_ligand_protein", "contact_protein_protein"]:
        expected.add(f"val/lddt_{m_}")
        expected.add(f"val/disto_lddt_{m_}")
        expected.add(f"val/complex_lddt_{m_}")
    expected.add("val/lddt")
    expected.add("val/disto_lddt")
    expected.add("val/complex_lddt")
    expected.add("val/disto_loss")
    if rmsd_metrics:
        expected.add("val/rmsd")
    if clash_score_metrics:
        expected.add("val/clash_atoms_count")
        expected.add("val/clash_atoms_fraction")

    if confidence_prediction:
        for m in const.out_single_types:
            expected.add(f"val/MAE_plddt_{m}")

        out_types_no_modified = [m for m in const.out_types if m != "modified"]
        for m_ in out_types_no_modified:
            expected.add(f"val/MAE_pde_{m_}")
            expected.add(f"val/MAE_pae_{m_}")

        conf_lddt_keys = [
            "top1_lddt",
            "iplddt_top1_lddt",
            "ipde_top1_lddt",
            "pde_top1_lddt",
            "ptm_top1_lddt",
            "iptm_top1_lddt",
            "ligand_iptm_top1_lddt",
            "protein_iptm_top1_lddt",
            "avg_lddt",
        ]
        for conf_key in conf_lddt_keys:
            for m_ in out_types_no_modified:
                expected.add(f"val/{conf_key}_{m_}")

    if physicalism_metrics:
        clash_keys = [f"asym_{m_}" for m_ in const.clash_types] + [f"sym_{m_}" for m_ in const.out_single_types]
        pb_keys = [
            "bond_length",
            "bond_angle",
            "internal_clash",
            "atom_chirality",
            "bond_stereochemistry",
            "ring_5_flatness",
            "ring_6_flatness",
            "double_bond_flatness",
        ]
        for m in clash_keys:
            expected.add(f"val/clash_{m}")
        for m in pb_keys:
            expected.add(f"val/pb_{m}")

        if confidence_prediction:
            prefixes = [
                "top1",
                "iplddt_top1",
                "pde_top1",
                "ipde_top1",
                "ptm_top1",
                "iptm_top1",
                "ligand_iptm_top1",
                "protein_iptm_top1",
                "avg",
            ]
            for prefix in prefixes:
                for m in clash_keys:
                    expected.add(f"val/{prefix}_clash_{m}")
                for m in pb_keys:
                    expected.add(f"val/{prefix}_pb_{m}")
    # --- end build expected ---

    logged_set = set(logged_names)

    missing = expected - logged_set
    extra = logged_set - expected
    assert not missing, f"Missing metrics: {sorted(missing)}"
    assert not extra, f"Unexpected metrics: {sorted(extra)}"

    DistributedManager.cleanup()
    mp.undo()


@pytest.mark.parametrize(
    "setup_env, confidence_prediction, physicalism_metrics, rmsd_metrics, clash_score_metrics",
    [
        (((1, (1, 1)), True, "cuda", "ENV"), False, False, False, False),
        (((1, (1, 1)), True, "cuda", "ENV"), True, False, False, False),
        (((1, (1, 1)), True, "cuda", "ENV"), False, True, False, False),
        (((1, (1, 1)), True, "cuda", "ENV"), False, False, True, False),
        (((1, (1, 1)), True, "cuda", "ENV"), False, False, False, True),
        (((1, (1, 1)), True, "cuda", "ENV"), True, True, True, True),
    ],
    indirect=("setup_env",),
    ids=["fold-only", "confidence", "physicalism", "rmsd", "clash-score", "all"],
)
def test_epoch_end_logged_metric_names(
    setup_env, confidence_prediction, physicalism_metrics, rmsd_metrics, clash_score_metrics
):
    """All expected metric names appear in model.log during on_epoch_end."""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    skip_if_cuda_not_avail_or_device_count_less_than_word_size(
        device_type=device_type,
        world_size=world_size,
    )

    payload = (
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        confidence_prediction,
        physicalism_metrics,
        rmsd_metrics,
        clash_score_metrics,
    )
    spawn_multiprocessing(_parallel_epoch_end_logged_metric_names_test, world_size, payload)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
