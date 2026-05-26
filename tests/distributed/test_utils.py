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


import unittest
from random import randint
from unittest.mock import patch

import pytest
import torch
from torch.distributed.tensor import DTensor, Replicate, Shard, distribute_tensor

from boltz.distributed.manager import DistributedManager
from boltz.distributed.utils import create_distributed_randn
from boltz.testing.utils import (
    assert_all_identical,
    distribute_atom_features,
    homogenize_shard_shapes,
    seed_by_rank,
    spawn_multiprocessing,
)


def assert_homogenize_shard_shapes_worker(rank: int, grid_group_sizes, world_size, device_type, backend, env_per_rank):
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

    # Define test cases inside the worker function
    test_cases = [
        # (global_shape, placements, value_to_pad)
        (
            (2 * device_mesh.size(0), 5 * device_mesh.size(1), 4 * device_mesh.size(2)),
            (Shard(0), Shard(1), Shard(2)),
            0.0,
        ),
        (
            (3 * device_mesh.size(0), 4 * device_mesh.size(1)),
            (Shard(0), Shard(1), Replicate()),
            None,
        ),
        (
            (3 * device_mesh.size(0), 4 * device_mesh.size(1)),
            (Shard(0), Replicate(), Shard(1)),
            None,
        ),
        (
            (5 * device_mesh.size(0), 7 * device_mesh.size(1)),
            (Replicate(), Shard(0), Shard(1)),
            0.0,
        ),
        (
            (5 * device_mesh.size(0), 7 * device_mesh.size(1)),
            (Replicate(), Replicate(), Replicate()),
            0.0,
        ),
    ]

    seed_by_rank(manager.group_rank["world"], seed=42)

    for global_shape, placements, value_to_pad in test_cases:
        label_test = f"global_shape={global_shape}, placements={placements}, value_to_pad={value_to_pad}"
        # Create a global tensor with known values for verification
        global_tensor = torch.randn(global_shape, dtype=torch.float32, device=manager.device)

        global_dtensor = distribute_tensor(
            global_tensor, device_mesh=device_mesh, placements=placements, src_data_rank=0
        )

        # slice local shard to make an expected tensor
        shard = global_dtensor.to_local()

        slice_local = []
        for i in range(shard.ndim):
            # randomly generate a slice towards the beginning along each tensor axis
            # so that the resulting equivalent padding is towards the end
            slice_local.append(slice(None, max(1, randint(1, shard.shape[i] - 1)), None))
        shard_sliced = shard[slice_local]

        if value_to_pad is None:
            shard_expected = torch.zeros_like(shard)
        else:
            shard_expected = torch.full_like(shard, torch.tensor(value_to_pad, dtype=shard.dtype))
        # consistent with the target function homogenize_shard_shapes's padding pattern:
        # always pad towards the end (or the last element) along each tensor axis
        shard_expected[slice_local] = shard_sliced

        global_heterogeneous_dtensor = DTensor.from_local(
            shard_sliced,
            device_mesh=global_dtensor.device_mesh,
            placements=global_dtensor.placements,
            shape=global_dtensor.shape,
            stride=global_dtensor.stride(),
        )

        expected = DTensor.from_local(
            shard_expected,
            device_mesh=global_dtensor.device_mesh,
            placements=global_dtensor.placements,
            shape=global_dtensor.shape,
            stride=global_dtensor.stride(),
        )

        results = homogenize_shard_shapes(global_heterogeneous_dtensor, value_to_pad=value_to_pad)

        torch.testing.assert_close(
            results.to_local(),
            expected.to_local(),
            atol=0,
            rtol=0,
            msg=lambda msg: f"{label_test}\n{msg}",
        )

        expected_full = expected.full_tensor()
        results_full = results.full_tensor()
        torch.testing.assert_close(
            expected_full,
            results_full,
            atol=0,
            rtol=0,
            msg=lambda msg: f"{label_test}\n{msg}",
        )

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env",
    [
        ((2, (2, 2)), True, "cuda", "ENV"),
        ((3, (3, 3)), True, "cpu", "ENV"),
    ],
    indirect=("setup_env",),
    ids=lambda x: f"dp={x[0][0]}, cp={x[0][1]}, device_type={x[2]}",
)
def test_homogenize_shard_shapes(setup_env):
    """Test homogenize_shard_shapes function with various sharding configurations."""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    spawn_multiprocessing(
        assert_homogenize_shard_shapes_worker,
        world_size,
        grid_group_sizes,
        world_size,
        device_type,
        backend,
        env_per_rank,
    )


def assert_create_distributed_randn(
    global_rank: int,
    shape: tuple[int, ...],
    placements: tuple[Shard | Replicate, ...],
    grid_group_sizes: dict[str, int],
    device_type: str,
    backend: str,
    env_per_rank: dict[str, str],
    dtype: torch.dtype,
    scale: float,
    seed: int,
    expected_source_ranks: set[int],
):
    """Assert correctness of create_distributed_randn on each rank."""
    monkeypatch = pytest.MonkeyPatch()
    for var_name, value in env_per_rank.items():
        if value == "<INPUT_RANK>":
            monkeypatch.setenv(var_name, f"{global_rank}")
            continue
        monkeypatch.setenv(var_name, value)

    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    manager = DistributedManager()

    # Use device_mesh_subgroups which has 3D structure: (dp, cp_axis_0, cp_axis_1)
    device_mesh = manager.device_mesh_subgroups

    # Set seed for reproducibility
    torch.manual_seed(seed + global_rank)

    # Use the expected source rank from parametrization
    is_source_rank = global_rank in expected_source_ranks

    # Track whether torch.randn was called
    randn_was_called = False
    original_randn = torch.randn

    def mock_randn(*args, **kwargs):
        nonlocal randn_was_called  # does not leak to other subprocesses
        randn_was_called = True
        return original_randn(*args, **kwargs)

    # Mock torch.randn and create distributed random tensor
    with patch("torch.randn", side_effect=mock_randn):
        dtensor = create_distributed_randn(
            shape=shape,
            device_mesh=device_mesh,
            placements=placements,
            dtype=dtype,
            scale=scale,
        )

    # Verify that torch.randn was called only on source ranks
    if is_source_rank:
        assert randn_was_called, f"Source rank {global_rank} should have called torch.randn"
    else:
        assert not randn_was_called, f"Non-source rank {global_rank} should not have called torch.randn"

    # Get local tensor
    local_tensor = dtensor.to_local()

    # Verify DTensor properties
    assert dtensor.device_mesh == device_mesh
    assert dtensor.placements == placements
    assert dtensor.shape == shape
    assert dtensor.dtype == dtype

    # Verify local shape is correct based on placements
    expected_local_shape = list(shape)
    for i_dim_mesh, placement in enumerate(placements):
        if placement.is_shard():
            expected_local_shape[placement.dim] = shape[placement.dim] // device_mesh.shape[i_dim_mesh]
    assert local_tensor.shape == tuple(expected_local_shape)

    # Verify replication: for each mesh dimension with Replicate placement,
    # all ranks along that mesh dimension should have identical values
    for mesh_dim, placement in enumerate(placements):
        if not placement.is_shard():
            # This mesh dimension is replicated, verify all ranks have identical local tensors
            assert_all_identical(
                local_tensor,
                device_mesh.get_group(mesh_dim),
                check_stride=False,
                check_grad=False,
                check_grad_fn=False,
                check_storage_offset=False,
                check_storage_pointer=False,
            )

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env,placements,shape,expected_source_ranks",
    [
        # ===== 3D tensors with 3D mesh =====
        # DP=1, CP=(2,2) tests - 3D mesh (1, 2, 2)
        # Mesh coordinates: rank 0:(0,0,0), 1:(0,0,1), 2:(0,1,0), 3:(0,1,1)
        ([(1, (2, 2)), True, "cpu", "ENV"], (Replicate(), Replicate(), Replicate()), (8, 16, 32), {0}),
        ([(1, (2, 2)), True, "cpu", "ENV"], (Replicate(), Shard(1), Replicate()), (8, 16, 32), {0, 2}),
        ([(1, (2, 2)), True, "cpu", "ENV"], (Replicate(), Replicate(), Shard(2)), (8, 16, 32), {0, 1}),
        ([(1, (2, 2)), True, "cpu", "ENV"], (Replicate(), Shard(0), Shard(1)), (8, 16, 32), {0, 1, 2, 3}),
        # DP=2, CP=(2,2) tests - 3D mesh (2, 2, 2)
        # Mesh coordinates: rank 0:(0,0,0), 1:(0,0,1), 2:(0,1,0), 3:(0,1,1), 4:(1,0,0), 5:(1,0,1), 6:(1,1,0), 7:(1,1,1)
        ([(2, (2, 2)), True, "cpu", "ENV"], (Shard(0), Replicate(), Replicate()), (8, 16, 32), {0, 4}),
        ([(2, (2, 2)), True, "cpu", "ENV"], (Shard(0), Shard(1), Shard(2)), (8, 16, 32), {0, 1, 2, 3, 4, 5, 6, 7}),
        ([(2, (2, 2)), True, "cpu", "ENV"], (Shard(0), Replicate(), Replicate()), (2,), {0, 4}),
        ([(2, (2, 2)), True, "cpu", "ENV"], (Shard(0), Replicate(), Replicate()), (16,), {0, 4}),
        # CUDA tests
        ([(1, (2, 2)), True, "cuda", "ENV"], (Replicate(), Shard(0), Shard(1)), (8, 16, 32), {0, 1, 2, 3}),
        ([(2, (2, 2)), True, "cuda", "ENV"], (Shard(0), Shard(1), Shard(2)), (8, 16, 32), {0, 1, 2, 3, 4, 5, 6, 7}),
        ([(2, (2, 2)), True, "cuda", "ENV"], (Shard(0), Replicate(), Replicate()), (2,), {0, 4}),
    ],
    indirect=("setup_env",),
    ids=[
        "dp:1_3d_all_replicated_cpu",
        "dp:1_3d_shard_dim1_cp0_cpu",
        "dp:1_3d_shard_dim2_cp1_cpu",
        "dp:1_3d_shard_dim01_cp01_cpu",
        "dp:2_3d_shard_dp_cpu",
        "dp:2_3d_shard_all_cpu",
        "dp:2_1d_batch2_shard_dp_cpu",
        "dp:2_1d_batch16_shard_dp_cpu",
        "dp:1_3d_shard_dim01_cp01_cuda",
        "dp:2_3d_shard_all_cuda",
        "dp:2_1d_batch2_shard_dp_cuda",
    ],
)
def test_create_distributed_randn(
    setup_env: dict[str, int],
    placements: tuple[Shard | Replicate, ...],
    shape: tuple[int, ...],
    expected_source_ranks: set[int],
    dtype: torch.dtype = torch.float32,
    scale: float = 1.0,
    seed: int = 42,
):
    """Test create_distributed_randn with various configurations.

    This test covers:
    - 3D tensors with various sharding patterns
    - 1D tensors (diffusion noise use case) with batch_size=1 edge case
    - Per-rank seeds to mimic training scenario
    - Verification of replication correctness across mesh dimensions
    - Verification that only source ranks call torch.randn
    """
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda" and torch.cuda.device_count() < world_size:
        pytest.skip("Not enough GPUs available")

    # Run test across all ranks
    spawn_multiprocessing(
        assert_create_distributed_randn,
        world_size,
        shape,
        placements,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        dtype,
        scale,
        seed,
        expected_source_ranks,
    )


def parallel_assert_distribute_atom_features(
    rank: int,
    grid_group_sizes: dict[str, int],
    world_size: int,
    device_type: str,
    backend: str,
    env_per_rank: dict[str, str],
    # Test data passed from main test function
    feats_global_host: dict[str, torch.Tensor],
    placements_cp: dict[str, tuple],
    placements_dp_cp: dict[str, tuple],
    multiplicities: dict[str, int] | None,
):
    """Parallel worker for testing distribute_atom_features.

    Tests that:
    1. DTensors are created correctly
    2. full_tensor() reconstructs data that matches original (accounting for interspersed padding)
    3. Non-padded values are preserved correctly
    """
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

    # Get test parameters
    size_batch = feats_global_host["atom_pad_mask"].shape[0]
    n_tokens = feats_global_host["token_pad_mask"].shape[1]
    atom_counts_per_token = feats_global_host["atom_counts_per_token"]

    # Convert inputs to device
    dtype = torch.float64
    inputs = {
        k: v.to(device=manager.device, dtype=dtype if v.dtype.is_floating_point else v.dtype)
        for k, v in feats_global_host.items()
        if k in placements_cp
    }

    # Call distribute_atom_features
    feats_dtensor = distribute_atom_features(
        inputs=inputs,
        placements_cp=placements_cp,
        placements_dp_cp=placements_dp_cp,
        device_mesh=device_mesh,
        cp_group=manager.group["cp"],
        multiplicities=multiplicities,
    )

    # Compute shard metadata for extracting non-padded values
    n_rows = device_mesh.size(1)  # cp_axis_0
    n_tokens_per_shard = n_tokens // n_rows
    token_atom_count_cumsum = torch.cat([torch.tensor([0]), atom_counts_per_token[0].cumsum(dim=0)])
    shard_atom_counts = atom_counts_per_token[0].unflatten(0, (n_rows, n_tokens_per_shard)).sum(dim=1)
    max_atoms_per_shard = shard_atom_counts.max().item()
    n_atoms_padded = max_atoms_per_shard * n_rows

    # Build mapping from padded indices to original indices for atom dimension
    # Interspersed padding: shard i has atoms at padded indices [i*max_atoms_per_shard : i*max_atoms_per_shard + actual_atoms_in_shard]
    # and these correspond to original indices [token_atom_count_cumsum[i*n_tokens_per_shard] : token_atom_count_cumsum[(i+1)*n_tokens_per_shard]]
    padded_to_original = torch.full((n_atoms_padded,), -1, dtype=torch.long)
    for i_shard in range(n_rows):
        token_start = n_tokens_per_shard * i_shard
        token_end = n_tokens_per_shard * (i_shard + 1)
        atom_start_orig = token_atom_count_cumsum[token_start].item()
        atom_end_orig = token_atom_count_cumsum[token_end].item()
        actual_atoms = atom_end_orig - atom_start_orig

        padded_start = i_shard * max_atoms_per_shard
        for j in range(actual_atoms):
            padded_to_original[padded_start + j] = atom_start_orig + j

    # Extract valid (non-padding) indices
    valid_atom_mask = padded_to_original >= 0
    valid_padded_indices = torch.where(valid_atom_mask)[0]
    valid_original_indices = padded_to_original[valid_atom_mask]

    # Test each output DTensor
    for k, dtensor in feats_dtensor.items():
        if k in {"atom_to_token", "token_to_rep_atom"}:
            # Block-diagonal matrices need special handling
            # These matrices have placement (Shard(0), Shard(1), Replicate()) but their structure
            # means each shard (i, j) contains the diagonal block (i, i).
            # The token dimension in the local shard is local (N_tokens_per_shard), not global (n_tokens).
            # full_tensor() reconstructs [B, n_atoms_padded, N_tokens_per_shard] (NOT n_tokens).
            #
            # For proper comparison, we extract local shards and compare with the corresponding
            # diagonal blocks from the reference.

            ref = feats_global_host[k].cpu()
            local_shard = dtensor.to_local().cpu()

            # Get this rank's position in the mesh
            rank_dp = manager.group_rank["dp"]
            rank_cp_0 = manager.group_rank["cp_axis_0"]

            # For block-diagonal: the token range is determined by cp_axis_0 (row index)
            token_start = n_tokens_per_shard * rank_cp_0
            token_end = n_tokens_per_shard * (rank_cp_0 + 1)
            atom_start_orig = token_atom_count_cumsum[token_start].item()
            atom_end_orig = token_atom_count_cumsum[token_end].item()
            actual_atoms = atom_end_orig - atom_start_orig

            # local_shard shape: [1, max_atoms_per_shard, N_tokens_per_shard] or [1, N_tokens_per_shard, max_atoms_per_shard]
            if k == "atom_to_token":
                # Shape: [1, max_atoms_per_shard, N_tokens_per_shard]
                # Extract non-padded part
                local_valid = local_shard[0, :actual_atoms, :]

                # Reference block
                ref_block = ref[rank_dp, atom_start_orig:atom_end_orig, token_start:token_end]

                torch.testing.assert_close(
                    local_valid,
                    ref_block.to(local_valid.dtype),
                    atol=0,
                    rtol=0,
                    msg=lambda m, k=k: f"feature {k}:\n  {m}",
                )

            elif k == "token_to_rep_atom":
                # Shape: [1, N_tokens_per_shard, max_atoms_per_shard]
                # Extract non-padded part
                local_valid = local_shard[0, :, :actual_atoms]

                # Reference block
                ref_block = ref[rank_dp, token_start:token_end, atom_start_orig:atom_end_orig]

                torch.testing.assert_close(
                    local_valid,
                    ref_block.to(local_valid.dtype),
                    atol=0,
                    rtol=0,
                    msg=lambda m, k=k: f"feature {k}:\n  {m}",
                )

        elif k == "pair_mask":
            full_tensor = dtensor.full_tensor().cpu()
            # pair_mask has shape [B, n_atoms_padded, n_atoms_padded] after full_tensor
            # Need to extract valid [atom_i, atom_j] pairs
            ref = feats_global_host[k].cpu()
            for b in range(size_batch):
                for i_shard in range(n_rows):
                    for j_shard in range(n_rows):
                        # Row shard atoms
                        row_token_start = n_tokens_per_shard * i_shard
                        row_token_end = n_tokens_per_shard * (i_shard + 1)
                        shard_atom_start = token_atom_count_cumsum[row_token_start].item()
                        shard_atom_end = token_atom_count_cumsum[row_token_end].item()
                        shard_actual_atoms = shard_atom_end - shard_atom_start
                        row_padded_start = i_shard * max_atoms_per_shard

                        # Col shard atoms
                        col_token_start = n_tokens_per_shard * j_shard
                        col_token_end = n_tokens_per_shard * (j_shard + 1)
                        col_atom_start = token_atom_count_cumsum[col_token_start].item()
                        col_atom_end = token_atom_count_cumsum[col_token_end].item()
                        col_actual_atoms = col_atom_end - col_atom_start
                        col_padded_start = j_shard * max_atoms_per_shard

                        # Extract block
                        block_full = full_tensor[
                            b,
                            row_padded_start : row_padded_start + shard_actual_atoms,
                            col_padded_start : col_padded_start + col_actual_atoms,
                        ]
                        block_ref = ref[b, shard_atom_start:shard_atom_end, col_atom_start:col_atom_end]

                        torch.testing.assert_close(
                            block_full,
                            block_ref.to(block_full.dtype),
                            atol=0,
                            rtol=0,
                            msg=lambda m, k=k: f"feature {k}:\n  {m}",
                        )

        else:
            # Single representation features (atom_pad_mask, ref_pos, ref_charge, etc.)
            # Shape: [B, n_atoms_padded, ...] - extract valid atoms using valid_padded_indices
            full_tensor = dtensor.full_tensor().cpu()
            ref = feats_global_host[k].cpu()

            # Handle different tensor dimensions
            if full_tensor.ndim == 2:
                # Shape: [B, n_atoms_padded]
                extracted = full_tensor[:, valid_padded_indices]
                # Reorder to original order
                reordered = torch.empty_like(ref, dtype=full_tensor.dtype)
                for idx, orig_idx in enumerate(valid_original_indices):
                    reordered[:, orig_idx] = extracted[:, idx]

                torch.testing.assert_close(
                    reordered,
                    ref.to(reordered.dtype),
                    atol=0,
                    rtol=0,
                    msg=lambda m, k=k: f"feature {k}:\n  {m}",
                )

            elif full_tensor.ndim == 3:
                # Shape: [B, n_atoms_padded, D] (e.g., ref_pos with D=3)
                extracted = full_tensor[:, valid_padded_indices, :]
                reordered = torch.empty_like(ref, dtype=full_tensor.dtype)
                for idx, orig_idx in enumerate(valid_original_indices):
                    reordered[:, orig_idx, :] = extracted[:, idx, :]

                torch.testing.assert_close(
                    reordered,
                    ref.to(reordered.dtype),
                    atol=0,
                    rtol=0,
                    msg=lambda m, k=k: f"feature {k}:\n  {m}",
                )

            elif full_tensor.ndim == 4:
                # Shape: [B, n_atoms_padded, D1, D2] (e.g., ref_atom_name_chars)
                extracted = full_tensor[:, valid_padded_indices, :, :]
                reordered = torch.empty_like(ref, dtype=full_tensor.dtype)
                for idx, orig_idx in enumerate(valid_original_indices):
                    reordered[:, orig_idx, :, :] = extracted[:, idx, :, :]

                torch.testing.assert_close(
                    reordered,
                    ref.to(reordered.dtype),
                    atol=0,
                    rtol=0,
                    msg=lambda m, k=k: f"feature {k}:\n  {m}",
                )

            elif full_tensor.ndim == 5:
                # Shape: [B, n_atoms_padded, D1, D2, D3] (e.g., ref_element one-hot)
                extracted = full_tensor[:, valid_padded_indices, :, :, :]
                reordered = torch.empty_like(ref, dtype=full_tensor.dtype)
                for idx, orig_idx in enumerate(valid_original_indices):
                    reordered[:, orig_idx, :, :, :] = extracted[:, idx, :, :, :]

                torch.testing.assert_close(
                    reordered,
                    ref.to(reordered.dtype),
                    atol=0,
                    rtol=0,
                    msg=lambda m, k=k: f"feature {k}:\n  {m}",
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
def test_distribute_atom_features(setup_env):
    """Test distribute_atom_features with various atom feature types.

    Tests:
    1. Single representation features (atom_pad_mask, ref_pos, ref_charge, etc.)
    2. Block-diagonal features (atom_to_token, token_to_rep_atom)
    3. Pair representation features (pair_mask)

    Verifies that after distribute_atom_features + full_tensor():
    - Non-padded values match the original input
    - Interspersed padding structure is correct
    """
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("CUDA not available")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"Not enough GPUs: need {world_size}, have {torch.cuda.device_count()}")

    # Use dp_size as batch size (required by distribute_atom_features)
    dp_size = grid_group_sizes["dp"]
    cp_tuple = grid_group_sizes["cp"]

    # Test parameters - tokens must be divisible by CP row size
    n_tokens = 8 * cp_tuple[0]  # e.g., 16 for cp=(2,2)
    # n_atoms must be large enough to accommodate the random atom counts per token
    # With (min=1, max=4) atoms per token, we need at least n_tokens * max_atoms + some buffer for last token
    n_msa = 4
    atom_counts_per_token_range = (1, 3)  # 1-3 atoms per token
    max_atoms_per_token = atom_counts_per_token_range[1]
    # Ensure enough atoms: (n_tokens - 1) * max + at least 1 for last token
    n_atoms = n_tokens * max_atoms_per_token + max_atoms_per_token  # extra buffer

    # Generate random features
    from boltz.testing.utils import random_features

    # Seed for reproducibility
    torch.manual_seed(42)

    feats_global_host = random_features(
        size_batch=dp_size,
        n_tokens=n_tokens,
        n_atoms=n_atoms,
        n_msa=n_msa,
        atom_counts_per_token_range=atom_counts_per_token_range,
        device=torch.device("cpu"),
        float_value_range=(-1.0, 1.0),
        selected_keys=[
            "atom_counts_per_token",
            "atom_pad_mask",
            "atom_to_token",
            "token_to_rep_atom",
            "pair_mask",
            "ref_pos",
            "ref_charge",
            "ref_element",
            "ref_atom_name_chars",
            "ref_space_uid",
            "token_pad_mask",  # Needed for reference
        ],
    )

    # Get actual n_atoms after random_features adjustment
    n_atoms = feats_global_host["atom_pad_mask"].shape[1]

    # Define placements for CP submesh (2-tuple)
    placements_cp_single = (Shard(0), Replicate())
    placements_cp_pair = (Shard(0), Shard(1))

    placements_cp = {
        "atom_counts_per_token": placements_cp_single,
        "atom_pad_mask": placements_cp_single,
        "atom_to_token": placements_cp_single,  # Block-diagonal, stored as single
        "token_to_rep_atom": placements_cp_single,  # Block-diagonal, stored as single
        "pair_mask": placements_cp_pair,
        "ref_pos": placements_cp_single,
        "ref_charge": placements_cp_single,
        "ref_element": placements_cp_single,
        "ref_atom_name_chars": placements_cp_single,
        "ref_space_uid": placements_cp_single,
    }

    # Define placements for full mesh (3-tuple: dp, cp_0, cp_1)
    placements_single = (Shard(0), Shard(1), Replicate())
    placements_pair = (Shard(0), Shard(1), Shard(2))

    placements_dp_cp = {
        "atom_pad_mask": placements_single,
        "atom_to_token": placements_single,
        "token_to_rep_atom": placements_single,
        "pair_mask": placements_pair,
        "ref_pos": placements_single,
        "ref_charge": placements_single,
        "ref_element": placements_single,
        "ref_atom_name_chars": placements_single,
        "ref_space_uid": placements_single,
    }

    # No multiplicity for this basic test
    multiplicities = None

    spawn_multiprocessing(
        parallel_assert_distribute_atom_features,
        world_size,
        grid_group_sizes,
        world_size,
        device_type,
        backend,
        env_per_rank,
        feats_global_host,
        placements_cp,
        placements_dp_cp,
        multiplicities,
    )


if __name__ == "__main__":
    unittest.main()
