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


"""Tests for DTensor resolved_loss implementation.

This module tests the distributed implementation of resolved_loss which computes
binary cross-entropy loss for predicting whether atoms are resolved. The tests
verify numerical correctness against the serial implementation and proper gradient
computation.

The key challenge is handling the block-diagonal structure of token_to_rep_atom
with intersperse padding, where each CP rank only sees its local token-atom
correspondences. This requires using pad_and_scatter_atom_features_dtensor
instead of simple distribute_tensor for atom features.
"""

import unittest

import pytest
import torch
from torch.distributed.tensor import Replicate, Shard, distribute_tensor

from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.loss.confidencev2 import (
    resolved_loss,
    resolved_negative_log_likelihood,
)
from boltz.model.loss.confidencev2 import resolved_loss as serial_resolved_loss
from boltz.testing.utils import (
    distribute_atom_features,
    init_tensors_uniform,
    random_features,
    spawn_multiprocessing,
)


def parallel_assert_resolved_loss(
    rank: int,
    payload: tuple,
):
    """Parallel test function for resolved_loss.

    This function runs on each rank in the distributed setup and verifies that
    the DTensor implementation matches the serial reference.

    Uses pad_and_scatter_atom_features_dtensor for proper atom feature sharding
    with intersperse padding.
    """
    (
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        pred_resolved_global_host,
        token_to_rep_atom_global_host,
        atom_counts_per_token_host,
        token_pad_mask_global_host,
        true_coords_resolved_mask_global_host,
        expected_loss_host,
        expected_pred_grad_host,
        multiplicity,
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

    device_mesh = manager.device_mesh_subgroups
    dtype = pred_resolved_global_host.dtype

    # Distribute pred_resolved: (B*mult, N_token, 2) with (Shard(0), Shard(1), Replicate())
    pred_resolved_dtensor = distribute_tensor(
        pred_resolved_global_host.to(manager.device),
        device_mesh=device_mesh,
        placements=(Shard(0), Shard(1), Replicate()),
    ).requires_grad_(True)

    # Distribute token_pad_mask: (B, N_token) with (Shard(0), Shard(1), Replicate())
    token_pad_mask_dtensor = distribute_tensor(
        token_pad_mask_global_host.to(manager.device),
        device_mesh=device_mesh,
        placements=(Shard(0), Shard(1), Replicate()),
    )

    # --- Distribute atom features using distribute_atom_features utility ---
    # Prepare inputs dict with all atom features (including per-multiplicity keys)
    size_batch = token_to_rep_atom_global_host.shape[0]
    inputs_atom = {
        "atom_counts_per_token": atom_counts_per_token_host.to(dtype=torch.int64),
        "token_to_rep_atom": token_to_rep_atom_global_host.to(dtype=dtype),
    }
    # Add per-multiplicity resolved masks: unflatten [B*mult, N_atom] -> [B, mult, N_atom]
    resolved_mask_unflat = true_coords_resolved_mask_global_host.unflatten(0, (size_batch, multiplicity))
    for i_mul in range(multiplicity):
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
        placements_cp[f"true_coords_resolved_mask_{i_mul}"] = (Shard(0), Replicate())
        placements_dp_cp[f"true_coords_resolved_mask_{i_mul}"] = (Shard(0), Shard(1), Replicate())

    # Distribute atom features with intersperse padding
    feats_atom = distribute_atom_features(
        inputs_atom,
        placements_cp,
        placements_dp_cp,
        device_mesh,
        manager.group["cp"],
        multiplicities={"true_coords_resolved_mask": multiplicity},
    )

    true_coords_resolved_mask_dtensor = feats_atom["true_coords_resolved_mask"]
    token_to_rep_atom_dtensor = feats_atom["token_to_rep_atom"]

    # Create feature dictionary
    feats_dtensor = {
        "token_to_rep_atom": token_to_rep_atom_dtensor,
        "token_pad_mask": token_pad_mask_dtensor,
    }

    # Create copies to verify inputs aren't modified
    pred_resolved_copy = pred_resolved_dtensor.to_local().detach().clone()

    # Forward pass
    loss_dtensor = resolved_loss(
        pred_resolved_dtensor,
        feats_dtensor,
        true_coords_resolved_mask_dtensor,
        multiplicity=multiplicity,
    )

    # Verify input wasn't modified (values only, not requires_grad)
    torch.testing.assert_close(
        pred_resolved_copy,
        pred_resolved_dtensor.to_local().detach(),
    )

    # Verify forward pass results
    expected_loss_dtensor = distribute_tensor(
        expected_loss_host.to(manager.device),
        device_mesh=device_mesh,
        placements=(Replicate(), Replicate(), Replicate()),
        src_data_rank=None,
    )

    assert (
        loss_dtensor.shape == expected_loss_dtensor.shape
    ), f"Loss shape mismatch: expected {expected_loss_dtensor.shape}, got {loss_dtensor.shape}"
    torch.testing.assert_close(
        loss_dtensor.to_local(),
        expected_loss_dtensor.to_local(),
    )

    # Backward pass
    loss_dtensor.backward()

    # Verify gradient
    assert pred_resolved_dtensor.grad is not None, "Gradient not computed for pred_resolved"
    assert (
        pred_resolved_dtensor.grad.shape == pred_resolved_global_host.shape
    ), f"Grad shape mismatch: expected {pred_resolved_global_host.shape}, got {pred_resolved_dtensor.grad.shape}"

    # Gather and compare gradients
    grad_global_result = pred_resolved_dtensor.grad.full_tensor().cpu()
    torch.testing.assert_close(
        grad_global_result,
        expected_pred_grad_host,
    )

    # Verify full tensor matches expected
    loss_global_result = loss_dtensor.full_tensor().cpu()
    torch.testing.assert_close(loss_global_result, expected_loss_host)

    DistributedManager.cleanup()
    monkeypatch.undo()


def parallel_assert_resolved_nll(
    rank: int,
    payload: tuple,
):
    """Parallel test function for resolved_negative_log_likelihood.

    This tests the shardwise NLL computation in isolation with multiplicity support.
    Uses pad_and_scatter_atom_features_dtensor for proper atom feature sharding.
    """
    (
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        pred_resolved_global_host,
        token_to_rep_atom_global_host,
        atom_counts_per_token_host,
        true_coords_resolved_mask_global_host,
        expected_errors_host,
        expected_d_errors_host,
        expected_pred_grad_host,
        multiplicity,
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

    device_mesh = manager.device_mesh_subgroups
    dtype = pred_resolved_global_host.dtype

    # Distribute pred_resolved: (B*mult, N_token, 2) with (Shard(0), Shard(1), Replicate())
    pred_resolved_dtensor = distribute_tensor(
        pred_resolved_global_host.to(manager.device),
        device_mesh=device_mesh,
        placements=(Shard(0), Shard(1), Replicate()),
    ).requires_grad_(True)

    # --- Distribute atom features using distribute_atom_features utility ---
    size_batch = token_to_rep_atom_global_host.shape[0]
    inputs_atom = {
        "atom_counts_per_token": atom_counts_per_token_host.to(dtype=torch.int64),
        "token_to_rep_atom": token_to_rep_atom_global_host.to(dtype=dtype),
    }
    # Add per-multiplicity resolved masks: unflatten [B*mult, N_atom] -> [B, mult, N_atom]
    resolved_mask_unflat = true_coords_resolved_mask_global_host.unflatten(0, (size_batch, multiplicity))
    for i_mul in range(multiplicity):
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
        placements_cp[f"true_coords_resolved_mask_{i_mul}"] = (Shard(0), Replicate())
        placements_dp_cp[f"true_coords_resolved_mask_{i_mul}"] = (Shard(0), Shard(1), Replicate())

    # Distribute atom features with intersperse padding
    feats_atom = distribute_atom_features(
        inputs_atom,
        placements_cp,
        placements_dp_cp,
        device_mesh,
        manager.group["cp"],
        multiplicities={"true_coords_resolved_mask": multiplicity},
    )

    true_coords_resolved_mask_dtensor = feats_atom["true_coords_resolved_mask"]
    token_to_rep_atom_dtensor = feats_atom["token_to_rep_atom"]

    # Forward pass
    errors_dtensor = resolved_negative_log_likelihood(
        pred_resolved_dtensor,
        token_to_rep_atom_dtensor,
        true_coords_resolved_mask_dtensor,
    )

    # Verify forward pass
    expected_errors_dtensor = distribute_tensor(
        expected_errors_host.to(manager.device),
        device_mesh=device_mesh,
        placements=(Shard(0), Shard(1), Replicate()),
    )

    torch.testing.assert_close(
        errors_dtensor.to_local(),
        expected_errors_dtensor.to_local(),
    )

    # Backward pass with custom gradient
    d_errors = distribute_tensor(
        expected_d_errors_host.to(manager.device),
        device_mesh=device_mesh,
        placements=(Shard(0), Shard(1), Replicate()),
    )
    errors_dtensor.backward(d_errors)

    # Verify gradient
    assert pred_resolved_dtensor.grad is not None, "Gradient not computed"
    grad_global_result = pred_resolved_dtensor.grad.full_tensor().cpu()
    torch.testing.assert_close(
        grad_global_result,
        expected_pred_grad_host,
    )

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env",
    [
        ((1, (2, 2)), True, "cuda", "ENV"),
        ((2, (2, 2)), True, "cuda", "ENV"),
    ],
    indirect=("setup_env",),
    ids=lambda x: f"dp={x[0][0]}, cp={x[0][1]}, device_type={x[2]}",
)
@pytest.mark.parametrize("multiplicity", [1, 2])
def test_resolved_loss(setup_env, multiplicity):
    """Test resolved_loss DTensor implementation against serial reference."""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() < {world_size}")

    seed = 42
    torch.manual_seed(seed)

    # Test dimensions
    dp_size = grid_group_sizes["dp"]
    cp_size = grid_group_sizes["cp"][0]

    # Batch size must equal DP size for per-sample processing in pad_and_scatter
    B = dp_size
    n_tokens_per_shard = 16
    N_token = n_tokens_per_shard * cp_size
    # Atoms: use variable atoms per token (1-3 atoms per token)
    n_atoms_per_token_min, n_atoms_per_token_max = 1, 3
    # Estimate total atoms: avg 2 atoms/token * N_token, rounded to be divisible by cp_size
    avg_atoms_per_token = (n_atoms_per_token_min + n_atoms_per_token_max) / 2
    N_atom = int(avg_atoms_per_token * N_token)
    # Make N_atom divisible by cp_size for even sharding
    N_atom = ((N_atom + cp_size - 1) // cp_size) * cp_size
    dtype = torch.float64

    # Use random_features to generate token_to_rep_atom with proper block-diagonal structure
    # where each token randomly picks one of its owned atoms as representative
    feats = random_features(
        size_batch=B,
        n_tokens=N_token,
        n_atoms=N_atom,
        n_msa=1,  # Not used for this test
        atom_counts_per_token_range=(n_atoms_per_token_min, n_atoms_per_token_max),
        device=torch.device(device_type),
        float_value_range=(-0.1, 0.1),
        selected_keys=["token_to_rep_atom", "token_pad_mask", "atom_pad_mask", "atom_counts_per_token"],
    )

    token_to_rep_atom_global = feats["token_to_rep_atom"].to(dtype=dtype)
    token_pad_mask_global = feats["token_pad_mask"].to(dtype=dtype)
    atom_counts_per_token = feats["atom_counts_per_token"]

    # Add some padding at the end of each shard for testing
    for shard in range(cp_size):
        end_idx = (shard + 1) * n_tokens_per_shard
        token_pad_mask_global[:, end_idx - 2 : end_idx] = 0

    # Create pred_resolved: (B*mult, N_token, 2)
    pred_resolved_global = torch.empty(B * multiplicity, N_token, 2, device=device_type, dtype=dtype)
    init_tensors_uniform([pred_resolved_global], low=-0.5, high=0.5)
    pred_resolved_global.requires_grad_(True)

    # Create true_coords_resolved_mask: (B*mult, N_atom)
    true_coords_resolved_mask_global = torch.randint(0, 2, (B * multiplicity, N_atom), device=device_type, dtype=dtype)

    # Create feature dictionary for serial computation
    feats_global = {
        "token_to_rep_atom": token_to_rep_atom_global,
        "token_pad_mask": token_pad_mask_global,
    }

    # Compute serial reference
    expected_loss = serial_resolved_loss(
        pred_resolved_global,
        feats_global,
        true_coords_resolved_mask_global,
        token_level_confidence=True,
        multiplicity=multiplicity,
    )
    expected_loss.backward()
    expected_loss_host = expected_loss.detach().clone().cpu()
    expected_pred_grad_host = pred_resolved_global.grad.detach().clone().cpu()

    # Prepare payload for parallel test
    # Note: atom_counts_per_token is needed for pad_and_scatter_atom_features_dtensor
    payload = (
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        pred_resolved_global.detach().clone().cpu(),
        token_to_rep_atom_global.clone().cpu(),
        atom_counts_per_token.clone().cpu(),
        token_pad_mask_global.clone().cpu(),
        true_coords_resolved_mask_global.clone().cpu(),
        expected_loss_host,
        expected_pred_grad_host,
        multiplicity,
    )

    # Launch parallel test
    spawn_multiprocessing(parallel_assert_resolved_loss, world_size, payload)


@pytest.mark.parametrize(
    "setup_env",
    [
        ((1, (2, 2)), True, "cuda", "ENV"),
        ((2, (2, 2)), True, "cuda", "ENV"),
    ],
    indirect=("setup_env",),
    ids=lambda x: f"dp={x[0][0]}, cp={x[0][1]}, device_type={x[2]}",
)
def test_resolved_negative_log_likelihood(setup_env):
    """Test resolved_negative_log_likelihood in isolation with multiplicity=2."""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() < {world_size}")

    seed = 123
    torch.manual_seed(seed)

    # Test dimensions
    dp_size = grid_group_sizes["dp"]
    cp_size = grid_group_sizes["cp"][0]
    multiplicity = 2  # Hardcoded multiplicity for this test

    # Batch size must equal DP size for per-sample processing
    B = dp_size
    n_tokens_per_shard = 8
    N_token = n_tokens_per_shard * cp_size
    # Atoms: use variable atoms per token (1-3 atoms per token)
    n_atoms_per_token_min, n_atoms_per_token_max = 1, 3
    avg_atoms_per_token = (n_atoms_per_token_min + n_atoms_per_token_max) / 2
    N_atom = int(avg_atoms_per_token * N_token)
    # Make N_atom divisible by cp_size for even sharding
    N_atom = ((N_atom + cp_size - 1) // cp_size) * cp_size
    dtype = torch.float64

    # Use random_features to generate token_to_rep_atom with proper structure
    feats = random_features(
        size_batch=B,
        n_tokens=N_token,
        n_atoms=N_atom,
        n_msa=1,
        atom_counts_per_token_range=(n_atoms_per_token_min, n_atoms_per_token_max),
        device=torch.device(device_type),
        float_value_range=(-0.1, 0.1),
        selected_keys=["token_to_rep_atom", "atom_counts_per_token"],
    )

    # token_to_rep_atom is NOT multiplexed: shape (B, N_token, N_atom)
    token_to_rep_atom_global = feats["token_to_rep_atom"].to(dtype=dtype)
    atom_counts_per_token = feats["atom_counts_per_token"]

    # Create pred_resolved with multiplicity: (B*mult, N_token, 2)
    pred_resolved_global = torch.empty(B * multiplicity, N_token, 2, device=device_type, dtype=dtype)
    init_tensors_uniform([pred_resolved_global], low=-0.5, high=0.5)
    pred_resolved_global.requires_grad_(True)

    # true_coords_resolved_mask with multiplicity: (B*mult, N_atom)
    true_coords_resolved_mask_global = torch.randint(0, 2, (B * multiplicity, N_atom), device=device_type, dtype=dtype)

    # Compute serial reference for NLL using einsum (matches the DTensor implementation)
    # token_to_rep_atom: (B, N_token, N_atom) -> "btj"
    # resolved_mask reshaped: (B, mult, N_atom) -> "bmj"
    # ref_mask: (B, mult, N_token) -> "bmt" then flatten to (B*mult, N_token)
    resolved_mask_reshaped = true_coords_resolved_mask_global.view(B, multiplicity, N_atom)
    ref_mask = torch.einsum("btj,bmj->bmt", token_to_rep_atom_global, resolved_mask_reshaped)
    ref_mask = ref_mask.flatten(0, 1)  # (B*mult, N_token)

    log_softmax_resolved = torch.nn.functional.log_softmax(pred_resolved_global, dim=-1)
    expected_errors = -ref_mask * log_softmax_resolved[:, :, 0] - (1 - ref_mask) * log_softmax_resolved[:, :, 1]

    # Create gradient for backward pass
    expected_d_errors = torch.empty_like(expected_errors)
    init_tensors_uniform([expected_d_errors], low=-0.5, high=0.5)

    # Backward with custom gradient
    expected_errors.backward(expected_d_errors)
    expected_errors_host = expected_errors.detach().clone().cpu()
    expected_d_errors_host = expected_d_errors.detach().clone().cpu()
    expected_pred_grad_host = pred_resolved_global.grad.detach().clone().cpu()

    # Prepare payload with atom_counts_per_token for pad_and_scatter
    payload = (
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        pred_resolved_global.detach().clone().cpu(),
        token_to_rep_atom_global.clone().cpu(),
        atom_counts_per_token.clone().cpu(),
        true_coords_resolved_mask_global.clone().cpu(),
        expected_errors_host,
        expected_d_errors_host,
        expected_pred_grad_host,
        multiplicity,
    )

    # Launch parallel test
    spawn_multiprocessing(parallel_assert_resolved_nll, world_size, payload)


if __name__ == "__main__":
    unittest.main()
