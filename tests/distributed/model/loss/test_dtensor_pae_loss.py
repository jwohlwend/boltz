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


"""Tests for DTensor pae_loss implementation."""

import unittest

import pytest
import torch
from torch.distributed.tensor import Replicate, Shard, distribute_tensor

from boltz.distributed.comm import One2OneComm
from boltz.distributed.data.utils import distribute_features
from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.loss.confidencev2 import pae_loss
from boltz.distributed.utils import get_group_rank_from_axial_shift
from boltz.model.loss.confidencev2 import pae_loss as serial_pae_loss
from boltz.testing.utils import (
    distribute_atom_features,
    init_tensors_uniform,
    random_features,
    seed_by_rank,
    spawn_multiprocessing,
)


def create_heterogeneous_pae_features(
    B: int,
    N_token: int,
    N_atom: int,
    multiplicity: int,
    device: torch.device,
    dtype: torch.dtype,
    base_feats: dict,
) -> tuple[dict, torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
    """Create heterogeneous features for PAE loss testing.

    Generates raw unpadded data with natural heterogeneity from base_feats.
    The heterogeneity comes from:
    - Different atom_counts_per_token distributions (from random_features)
    - Mixed polymer/nonpolymer mol_type (from random_features)
    - Different resolution masks per batch

    Returns:
        feats_global: Feature dictionary (raw, unpadded)
        pred_pae: Predicted PAE logits
        pred_atom_coords: Predicted coordinates (raw, unpadded)
        true_atom_coords: True coordinates (raw, unpadded)
        true_coords_resolved_mask: Resolution mask (raw, unpadded)
    """
    num_bins = 64

    feats_global = {
        "frames_idx": base_feats["frames_idx"].to(device=device),
        "frame_resolved_mask": base_feats["frame_resolved_mask"].to(device=device, dtype=dtype),
        "asym_id": base_feats["asym_id"].to(device=device),
        "atom_to_token": base_feats["atom_to_token"].to(device=device, dtype=dtype),
        "atom_pad_mask": base_feats["atom_pad_mask"].to(device=device, dtype=dtype),
        "mol_type": base_feats["mol_type"].to(device=device),
        "token_pad_mask": base_feats["token_pad_mask"].to(device=device, dtype=dtype),
        "atom_resolved_mask": base_feats["atom_resolved_mask"].to(device=device, dtype=dtype),
        "is_nonpolymer_with_frame": base_feats["is_nonpolymer_with_frame"].to(device=device),
        "atom_counts_per_token": base_feats["atom_counts_per_token"].to(device=device),
    }

    # Main tensors - raw unpadded
    pred_pae = torch.empty(B, multiplicity, N_token, N_token, num_bins, device=device, dtype=dtype)
    init_tensors_uniform([pred_pae], low=-0.5, high=0.5)

    pred_atom_coords = torch.empty(B * multiplicity, N_atom, 3, device=device, dtype=dtype)
    true_atom_coords = torch.empty(B * multiplicity, N_atom, 3, device=device, dtype=dtype)
    init_tensors_uniform([pred_atom_coords, true_atom_coords], low=-10.0, high=10.0)

    # Heterogeneous resolution masks per batch, repeated for each multiplicity copy.
    # pae_loss indexes with arange(0, B*mult, mult) to pick one mask per sample.
    true_coords_resolved_mask = torch.ones(B * multiplicity, N_atom, device=device, dtype=dtype)
    for b in range(B):
        unresolved_fraction = 0.2 + 0.1 * (b % 3)
        n_unresolved = int(N_atom * unresolved_fraction)
        unresolved_indices = torch.randperm(N_atom, device=device)[:n_unresolved]
        for m in range(multiplicity):
            true_coords_resolved_mask[b * multiplicity + m, unresolved_indices] = 0

    return feats_global, pred_pae, pred_atom_coords, true_atom_coords, true_coords_resolved_mask


def parallel_assert_pae_loss(
    rank: int,
    payload: tuple,
):
    """Parallel test function for pae_loss.

    This function runs on each rank in the distributed setup and verifies that
    the DTensor implementation matches the serial reference.
    """
    test_config, inputs_global, feats_global, expected = payload

    # Unpack test configuration
    grid_group_sizes = test_config["grid_group_sizes"]
    device_type = test_config["device_type"]
    backend = test_config["backend"]
    env_per_rank = test_config["env_per_rank"]
    multiplicity = test_config["multiplicity"]
    max_dist = test_config["max_dist"]

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
    group_layout = manager.layout_subgroups["cp"]
    dtype = inputs_global["pred_pae"].dtype

    # Setup communication for coordinate transpose
    rank_coords = group_layout.unravel(manager.group_rank["cp"])
    comm = One2OneComm(
        manager.group["cp"],
        rank_send_to=get_group_rank_from_axial_shift(rank_coords, 0, -1, group_layout),
        rank_recv_from=get_group_rank_from_axial_shift(rank_coords, 0, 1, group_layout),
    )

    # --- Distribute pred_pae (B*mult, N_token, N_token, bins) ---
    pred_pae_dtensor = distribute_tensor(
        inputs_global["pred_pae"].to(manager.device),
        device_mesh=device_mesh,
        placements=(Shard(0), Shard(1), Shard(2)),
    ).requires_grad_(True)

    # --- Distribute atom features with intersperse padding ---
    # Use distribute_atom_features for all atom-indexed tensors including coords
    size_batch = feats_global["atom_counts_per_token"].shape[0]
    pred_coords_unflat = inputs_global["pred_atom_coords"].unflatten(0, (size_batch, multiplicity))
    true_coords_unflat = inputs_global["true_atom_coords"].unflatten(0, (size_batch, multiplicity))
    true_mask_unflat = inputs_global["true_coords_resolved_mask"].unflatten(0, (size_batch, multiplicity))

    inputs_atom = {
        "atom_counts_per_token": feats_global["atom_counts_per_token"].to(dtype=torch.int64),
        "atom_to_token": feats_global["atom_to_token"].to(dtype=dtype),
        "atom_pad_mask": feats_global["atom_pad_mask"].to(dtype=dtype),
        "atom_resolved_mask": feats_global["atom_resolved_mask"].to(dtype=dtype),
        "frames_idx": feats_global["frames_idx"],
    }
    for i_mul in range(multiplicity):
        inputs_atom[f"pred_atom_coords_{i_mul}"] = pred_coords_unflat[:, i_mul].to(dtype=dtype)
        inputs_atom[f"true_atom_coords_{i_mul}"] = true_coords_unflat[:, i_mul].to(dtype=dtype)
        inputs_atom[f"true_coords_resolved_mask_{i_mul}"] = true_mask_unflat[:, i_mul].to(dtype=dtype)

    placements_cp = {
        "atom_counts_per_token": (Shard(0), Replicate()),
        "atom_to_token": (Shard(0), Replicate()),
        "atom_pad_mask": (Shard(0), Replicate()),
        "atom_resolved_mask": (Shard(0), Replicate()),
        "frames_idx": (Shard(1), Replicate()),
    }
    placements_dp_cp = {
        "atom_to_token": (Shard(0), Shard(1), Replicate()),
        "atom_pad_mask": (Shard(0), Shard(1), Replicate()),
        "atom_resolved_mask": (Shard(0), Shard(1), Replicate()),
        "frames_idx": (Shard(0), Shard(1), Replicate()),
    }
    for i_mul in range(multiplicity):
        placements_cp[f"pred_atom_coords_{i_mul}"] = (Shard(0), Replicate())
        placements_cp[f"true_atom_coords_{i_mul}"] = (Shard(0), Replicate())
        placements_cp[f"true_coords_resolved_mask_{i_mul}"] = (Shard(0), Replicate())
        placements_dp_cp[f"pred_atom_coords_{i_mul}"] = (Shard(0), Shard(1), Replicate())
        placements_dp_cp[f"true_atom_coords_{i_mul}"] = (Shard(0), Shard(1), Replicate())
        placements_dp_cp[f"true_coords_resolved_mask_{i_mul}"] = (Shard(0), Shard(1), Replicate())

    feats_atom = distribute_atom_features(
        inputs_atom,
        placements_cp,
        placements_dp_cp,
        device_mesh,
        manager.group["cp"],
        multiplicities={
            "pred_atom_coords": multiplicity,
            "true_atom_coords": multiplicity,
            "true_coords_resolved_mask": multiplicity,
        },
    )

    pred_atom_coords_dtensor = feats_atom["pred_atom_coords"].requires_grad_(True)
    true_atom_coords_dtensor = feats_atom["true_atom_coords"]
    true_coords_resolved_mask_dtensor = feats_atom["true_coords_resolved_mask"]

    # --- Distribute token features (not atom-indexed) ---
    token_placements = (Shard(0), Shard(1), Replicate())
    if manager.group_rank["world"] == 0:
        token_feats_global = {
            "frame_resolved_mask": feats_global["frame_resolved_mask"].to(manager.device),
            "asym_id": feats_global["asym_id"].to(manager.device),
            "mol_type": feats_global["mol_type"].to(manager.device),
            "token_pad_mask": feats_global["token_pad_mask"].to(manager.device),
            "is_nonpolymer_with_frame": feats_global["is_nonpolymer_with_frame"].to(manager.device),
        }
    else:
        token_feats_global = None

    token_placements_map = {
        "frame_resolved_mask": token_placements,
        "asym_id": token_placements,
        "mol_type": token_placements,
        "token_pad_mask": token_placements,
        "is_nonpolymer_with_frame": token_placements,
    }
    feats_token_dtensor = distribute_features(
        token_feats_global,
        token_placements_map,
        manager.group["world"],
        manager.group_ranks["world"][0],
        device_mesh,
    )

    feats_dtensor = {
        "frames_idx": feats_atom["frames_idx"],
        **feats_token_dtensor,
        # Atom features from distribute_atom_features
        "atom_to_token": feats_atom["atom_to_token"],
        "atom_pad_mask": feats_atom["atom_pad_mask"],
        "atom_resolved_mask": feats_atom["atom_resolved_mask"],
    }

    # --- Forward pass ---
    loss_dtensor = pae_loss(
        pred_pae_dtensor,
        pred_atom_coords_dtensor,
        true_atom_coords_dtensor,
        true_coords_resolved_mask_dtensor,
        feats_dtensor,
        comm=comm,
        dist_manager=manager,
        group_layout=group_layout,
        multiplicity=multiplicity,
        max_dist=max_dist,
    )

    # --- Verify forward pass results ---
    loss_actual = loss_dtensor.full_tensor().cpu()
    loss_expected = expected["loss"]
    torch.testing.assert_close(loss_actual, loss_expected)

    # --- Backward pass ---
    loss_dtensor.backward()

    # Verify pred_pae gradient
    assert pred_pae_dtensor.grad is not None, "Gradient not computed for pred_pae"
    torch.testing.assert_close(
        pred_pae_dtensor.grad.full_tensor().cpu(),
        expected["pred_pae_grad"],
    )

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env",
    [
        ((1, (3, 3)), True, "cpu", "ENV"),
        ((2, (2, 2)), True, "cpu", "ENV"),
    ],
    indirect=("setup_env",),
    ids=lambda x: f"dp={x[0][0]}, cp={x[0][1]}, device_type={x[2]}",
)
@pytest.mark.parametrize("multiplicity", [1, 2])
def test_pae_loss(setup_env, multiplicity):
    """Test pae_loss DTensor implementation against serial reference."""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() < {world_size}")

    seed_by_rank(0)

    # Test dimensions
    dp_size = grid_group_sizes["dp"]
    cp_size = grid_group_sizes["cp"][0]

    B = dp_size  # Batch size must equal DP size for per-sample processing
    n_tokens_per_shard = 16
    N_token = n_tokens_per_shard * cp_size
    n_atoms_per_token_min, n_atoms_per_token_max = 1, 18
    avg_atoms_per_token = (n_atoms_per_token_min + n_atoms_per_token_max) / 2
    N_atom = int(avg_atoms_per_token * N_token * 1.25)  # 25% buffer
    N_atom = ((N_atom + cp_size - 1) // cp_size) * cp_size
    max_dist = 32.0
    dtype = torch.float32
    device = torch.device(device_type)

    # Generate base features with random atom-token structure
    base_feats = random_features(
        size_batch=B,
        n_tokens=N_token,
        n_atoms=N_atom,
        n_msa=1,
        atom_counts_per_token_range=(n_atoms_per_token_min, n_atoms_per_token_max),
        device=device,
        float_value_range=(-0.1, 0.1),
        selected_keys=[
            "token_pad_mask",
            "atom_pad_mask",
            "atom_counts_per_token",
            "atom_to_token",
            "asym_id",
            "mol_type",
            "atom_resolved_mask",
            "frames_idx",
            "frame_resolved_mask",
            "is_nonpolymer_with_frame",
        ],
    )

    # Ensure frame_resolved_mask is not killed by random atom_resolved_mask:
    # set atom_resolved_mask to all-ones so frame_resolved_mask reflects only
    # frame geometry (collinear mask), not random atom resolution.
    # true_coords_resolved_mask (created below) still exercises the resolved path.
    base_feats["atom_resolved_mask"] = torch.ones_like(base_feats["atom_resolved_mask"])
    base_feats["frame_resolved_mask"] = torch.ones_like(base_feats["frame_resolved_mask"])

    # Create heterogeneous features - raw unpadded data
    feats_global, pred_pae, pred_atom_coords, true_atom_coords, true_coords_resolved_mask = (
        create_heterogeneous_pae_features(
            B=B,
            N_token=N_token,
            N_atom=N_atom,
            multiplicity=multiplicity,
            device=device,
            dtype=dtype,
            base_feats=base_feats,
        )
    )
    pred_pae.requires_grad_(True)
    pred_atom_coords.requires_grad_(True)

    # Compute serial reference on raw unpadded data
    # Serial pae_loss expects true_coords_resolved_mask with shape (B, N_atom)
    expected_loss, _rel_loss = serial_pae_loss(
        pred_pae=pred_pae,
        pred_atom_coords=pred_atom_coords,
        feats=feats_global,
        true_atom_coords=true_atom_coords,
        true_coords_resolved_mask=true_coords_resolved_mask,
        multiplicity=multiplicity,
        max_dist=max_dist,
    )
    expected_loss.backward()
    assert pred_pae.grad is not None and pred_pae.grad.abs().sum() > 0, "Serial grad is zero"
    # Pack payload as dicts
    test_config = {
        "grid_group_sizes": grid_group_sizes,
        "device_type": device_type,
        "backend": backend,
        "env_per_rank": env_per_rank,
        "multiplicity": multiplicity,
        "max_dist": max_dist,
    }

    # Reshape pred_pae from serial (B, mult, N, N, bins) to distributed (B*mult, N, N, bins)
    pred_pae_flat = pred_pae.detach().clone().flatten(0, 1).cpu()
    pred_pae_grad_flat = pred_pae.grad.detach().clone().flatten(0, 1).cpu()

    inputs_global = {
        "pred_pae": pred_pae_flat,
        "pred_atom_coords": pred_atom_coords.detach().clone().cpu(),
        "true_atom_coords": true_atom_coords.clone().cpu(),
        "true_coords_resolved_mask": true_coords_resolved_mask.clone().cpu(),
    }

    feats_global_cpu = {k: v.clone().cpu() for k, v in feats_global.items()}

    expected = {
        "loss": expected_loss.detach().clone().cpu(),
        "pred_pae_grad": pred_pae_grad_flat,
    }

    payload = (test_config, inputs_global, feats_global_cpu, expected)

    # Launch parallel test
    spawn_multiprocessing(parallel_assert_pae_loss, world_size, payload)


if __name__ == "__main__":
    unittest.main()
