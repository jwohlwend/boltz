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


"""Tests for DTensor pde_loss computation.

This module tests:
1. Whether DTensor pde_loss matches serial pde_loss
2. Whether gradients flow correctly through pred_pde

The pde_loss computation:
    token_to_rep_atom = feats["token_to_rep_atom"]
    token_mask = bmm(token_to_rep_atom, resolved_mask)
    mask = token_mask[:,:,None] * token_mask[:,None,:]

    true_token_coords = bmm(token_to_rep_atom, true_atom_coords)
    pred_token_coords = bmm(token_to_rep_atom, pred_atom_coords)

    true_d = cdist(true_token_coords, true_token_coords)
    pred_d = cdist(pred_token_coords, pred_token_coords)
    target_pde = abs(true_d - pred_d)

    bin_index = clamp(floor(target_pde * num_bins / max_dist), max=num_bins-1)
    errors = cross_entropy(pred_pde, bin_index)
    loss = sum(errors * mask) / sum(mask)
"""

import pytest
import torch
from torch.distributed.tensor import Replicate, Shard

from boltz.distributed.comm import TransposeComm
from boltz.distributed.data.utils import distribute_features
from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.loss.confidencev2 import pde_loss
from boltz.model.loss.confidencev2 import pde_loss as serial_pde_loss
from boltz.testing.utils import (
    distribute_atom_features,
    random_features,
    spawn_multiprocessing,
)


def parallel_assert_pde_loss(
    rank: int,
    payload: tuple,
):
    """Worker function that runs on each rank to test pde_loss DTensor implementation.

    Uses the same setup pattern as parallel_assert_plddt_loss.
    """
    (
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        pred_pde_global_host,
        pred_atom_coords_global_host,
        true_atom_coords_global_host,
        true_coords_resolved_mask_global_host,
        token_to_rep_atom_global_host,
        atom_counts_per_token_host,
        expected_loss_host,
        expected_grad_pred_pde_host,
        multiplicity,
    ) = payload

    # Setup environment variables for this rank
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
    dtype = pred_atom_coords_global_host.dtype

    # Create TransposeComm for redistribute_transpose
    comm = TransposeComm(manager.group["cp"], manager.layout_subgroups["cp"])

    # --- Distribute atom features using distribute_atom_features utility ---
    size_batch = token_to_rep_atom_global_host.shape[0]
    inputs_atom = {
        "atom_counts_per_token": atom_counts_per_token_host.to(dtype=torch.int64),
        "token_to_rep_atom": token_to_rep_atom_global_host.to(dtype=dtype),
    }

    # Add per-multiplicity coordinates and masks
    pred_coords_unflat = pred_atom_coords_global_host.unflatten(0, (size_batch, multiplicity))
    true_coords_unflat = true_atom_coords_global_host.unflatten(0, (size_batch, multiplicity))
    resolved_mask_unflat = true_coords_resolved_mask_global_host.unflatten(0, (size_batch, multiplicity))

    for i_mul in range(multiplicity):
        inputs_atom[f"pred_atom_coords_{i_mul}"] = pred_coords_unflat[:, i_mul].to(dtype=dtype)
        inputs_atom[f"true_atom_coords_{i_mul}"] = true_coords_unflat[:, i_mul].to(dtype=dtype)
        inputs_atom[f"true_coords_resolved_mask_{i_mul}"] = resolved_mask_unflat[:, i_mul].to(dtype=dtype)

    # Define placements for CP submesh and full mesh
    placements_cp = {
        "atom_counts_per_token": (Shard(0), Replicate()),
        "token_to_rep_atom": (Shard(0), Replicate()),
    }
    placements_dp_cp = {
        "token_to_rep_atom": (Shard(0), Shard(1), Replicate()),
    }
    for i_mul in range(multiplicity):
        placements_cp[f"pred_atom_coords_{i_mul}"] = (Shard(0), Replicate())
        placements_cp[f"true_atom_coords_{i_mul}"] = (Shard(0), Replicate())
        placements_cp[f"true_coords_resolved_mask_{i_mul}"] = (Shard(0), Replicate())
        placements_dp_cp[f"pred_atom_coords_{i_mul}"] = (Shard(0), Shard(1), Replicate())
        placements_dp_cp[f"true_atom_coords_{i_mul}"] = (Shard(0), Shard(1), Replicate())
        placements_dp_cp[f"true_coords_resolved_mask_{i_mul}"] = (Shard(0), Shard(1), Replicate())

    # Distribute atom features with intersperse padding
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

    # --- Distribute pred_pde (pair representation) ---
    # pred_pde: [B*mult, N_token, N_token, num_bins] with placements (S(0), S(1), S(2))
    if manager.group_rank["world"] == 0:
        pred_pde_features = {
            "pred_pde": pred_pde_global_host.to(device=manager.device, dtype=torch.float32),
        }
    else:
        pred_pde_features = None
    pred_pde_placements = {
        "pred_pde": (Shard(0), Shard(1), Shard(2)),
    }
    pred_pde_dtensor_dict = distribute_features(
        pred_pde_features,
        pred_pde_placements,
        manager.group["world"],
        manager.group_ranks["world"][0],
        device_mesh,
    )

    # Extract distributed tensors
    pred_atom_coords_dtensor = feats_atom["pred_atom_coords"]
    true_atom_coords_dtensor = feats_atom["true_atom_coords"]
    true_coords_resolved_mask_dtensor = feats_atom["true_coords_resolved_mask"]

    # Create feature dictionary
    feats_dtensor = {
        "token_to_rep_atom": feats_atom["token_to_rep_atom"],
    }

    # Get pred_pde DTensor with gradient tracking
    pred_pde_dtensor = pred_pde_dtensor_dict["pred_pde"]
    pred_pde_dtensor_grad = pred_pde_dtensor.detach().requires_grad_(True)

    # Compute pde_loss
    loss = pde_loss(
        pred_pde=pred_pde_dtensor_grad,
        pred_atom_coords=pred_atom_coords_dtensor,
        true_atom_coords=true_atom_coords_dtensor,
        true_coords_resolved_mask=true_coords_resolved_mask_dtensor,
        feats=feats_dtensor,
        comm=comm,
        multiplicity=multiplicity,
    )

    # Verify loss value
    loss_local = loss.to_local()
    # Match dtype of DTensor output
    expected_loss = expected_loss_host.to(device=loss_local.device, dtype=loss_local.dtype)
    torch.testing.assert_close(loss_local, expected_loss)

    # Verify gradients
    loss_local.backward()
    grad_pred_pde = pred_pde_dtensor_grad.grad

    # Full gather the gradient to compare
    grad_pred_pde_full = grad_pred_pde.full_tensor()
    # Match dtype of DTensor gradient output
    expected_grad = expected_grad_pred_pde_host.to(device=grad_pred_pde_full.device, dtype=grad_pred_pde_full.dtype)
    torch.testing.assert_close(grad_pred_pde_full, expected_grad)

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env",
    [
        ((1, (2, 2)), True, "cuda", "ENV"),
        ((2, (2, 2)), True, "cuda", "ENV"),
    ],
    indirect=("setup_env",),
    ids=lambda x: f"dp:{x[0][0]}, cp:{x[0][1]}, device_type:{x[2]}",
)
@pytest.mark.parametrize("multiplicity", [1, 2], ids=lambda x: f"multiplicity:{x}")
def test_pde_loss(setup_env, multiplicity: int):
    """Test that DTensor pde_loss matches serial reference.

    This test verifies:
    1. Forward pass: DTensor pde_loss matches serial pde_loss
    2. Backward pass: Gradients w.r.t. pred_pde match serial gradients

    The test uses realistic feature generation with proper block-diagonal structure.
    """
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() < {world_size}")

    torch.manual_seed(42)

    dp_size = grid_group_sizes["dp"]
    cp_size = grid_group_sizes["cp"][0] * grid_group_sizes["cp"][1]

    # Generate test data
    B = dp_size  # Batch size equals DP size
    N_token = 32
    N_atom = 140  # Large enough to accommodate token atom counts
    n_atoms_per_token_min = 1
    n_atoms_per_token_max = 4
    num_bins = 64  # Number of PDE bins
    max_dist = 32.0  # Max distance for binning
    dtype = torch.float32  # Use FP32 for testing

    # Make N_atom divisible by cp_size for even sharding
    N_atom = ((N_atom + cp_size - 1) // cp_size) * cp_size

    # Use random_features to generate features with proper block-diagonal structure
    feats = random_features(
        size_batch=B,
        n_tokens=N_token,
        n_atoms=N_atom,
        n_msa=1,
        atom_counts_per_token_range=(n_atoms_per_token_min, n_atoms_per_token_max),
        device=torch.device(device_type),
        float_value_range=(-0.1, 0.1),
        selected_keys=[
            "token_to_rep_atom",
            "atom_counts_per_token",
        ],
    )

    token_to_rep_atom_global = feats["token_to_rep_atom"].to(dtype=dtype)  # [B, N_token, N_atom]
    atom_counts_per_token = feats["atom_counts_per_token"]

    N_atom_actual = token_to_rep_atom_global.shape[2]

    # Generate coordinates using uniform distribution
    # Use [-10, 10] range so average pairwise distance is meaningful
    pred_atom_coords_global = torch.empty(B * multiplicity, N_atom_actual, 3, device=device_type, dtype=dtype).uniform_(
        -10.0, 10.0
    )
    true_atom_coords_global = torch.empty(B * multiplicity, N_atom_actual, 3, device=device_type, dtype=dtype).uniform_(
        -10.0, 10.0
    )

    # Create true_coords_resolved_mask: [B*mult, N_atom]
    true_coords_resolved_mask_global = torch.randint(
        0, 2, (B * multiplicity, N_atom_actual), device=device_type, dtype=dtype
    )
    # Ensure at least half atoms are resolved for meaningful test
    true_coords_resolved_mask_global[:, : N_atom_actual // 2] = 1.0

    # Generate pred_pde: (B*mult, N_token, N_token, num_bins)
    pred_pde_global = torch.randn(
        B * multiplicity, N_token, N_token, num_bins, device=device_type, dtype=torch.float32
    ).requires_grad_(True)

    # Compute serial reference (serial uses float32 internally)
    feats_serial = {
        "token_to_rep_atom": token_to_rep_atom_global.clone().float(),
    }

    expected_loss, _rel_loss = serial_pde_loss(
        pred_pde=pred_pde_global,
        pred_atom_coords=pred_atom_coords_global.float(),
        feats=feats_serial,
        true_atom_coords=true_atom_coords_global.float(),
        true_coords_resolved_mask=true_coords_resolved_mask_global.float(),
        multiplicity=multiplicity,
        max_dist=max_dist,
    )

    # Compute gradients for serial reference
    expected_loss.backward()
    expected_grad_pred_pde = pred_pde_global.grad.clone()

    # Verify that serial reference produces gradients (the gradient flows through log_softmax)
    assert expected_grad_pred_pde is not None, "Serial pde_loss should produce gradients"
    assert not torch.allclose(
        expected_grad_pred_pde, torch.zeros_like(expected_grad_pred_pde)
    ), "Gradients should be non-zero"

    # Prepare payload for parallel test
    payload = (
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        pred_pde_global.detach().clone().cpu(),
        pred_atom_coords_global.clone().cpu(),
        true_atom_coords_global.clone().cpu(),
        true_coords_resolved_mask_global.clone().cpu(),
        token_to_rep_atom_global.clone().cpu(),
        atom_counts_per_token.clone().cpu(),
        expected_loss.detach().clone().cpu(),
        expected_grad_pred_pde.clone().cpu(),
        multiplicity,
    )

    # Launch parallel test
    spawn_multiprocessing(parallel_assert_pde_loss, world_size, payload)


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
