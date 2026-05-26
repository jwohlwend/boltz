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

"""Tests for pack_atom_features function.

This module tests the pack_atom_features function which removes per-shard trailing
padding from pad_and_scatter_atom_features_dtensor and creates a packed DTensor
with global trailing padding (multiple of W * size_cp).

Tests verify that:
1. Packed features match serial reference in the valid (non-padding) region
2. The packing operation correctly handles variable atoms per token
3. atom_to_token is converted to global indices (atom_to_token_ids_global)
"""

import warnings

import pytest
import torch
from torch.distributed.tensor import DTensor
from torch.distributed.tensor._utils import compute_global_tensor_info

from boltz.distributed.data.feature.featurizer import pack_atom_features, pad_and_scatter_atom_features_dtensor
from boltz.distributed.manager import DistributedManager
from boltz.testing.utils import (
    assert_tensors_close_with_pad,
    assert_tensors_identical,
    get_feature_placements,
    random_features,
    seed_by_rank,
    spawn_multiprocessing,
)

# Subset of keys needed for pack_atom_features test
_selected_atom_keys = {
    "atom_pad_mask",
    "ref_pos",
    "ref_space_uid",
    "ref_charge",
    "ref_element",
    "ref_atom_name_chars",
    "atom_to_token",
    "atom_counts_per_token",  # Required by pad_and_scatter_atom_features_dtensor
}

# Get feature placements from centralized utility function with atom key subset
# Pass empty sets for unused categories to suppress irrelevant placements
_placements = get_feature_placements(
    token_keys=set(),
    msa_keys=set(),
    atom_keys=_selected_atom_keys,
    model_io_keys=set(),
    model_io_fp32_keys=set(),
)
_placements_cp_atom_features = _placements["cp_atom_features"]
_placements_atom_features = _placements["atom_features"]


def parallel_assert_pack_and_pad_atom_features(
    rank: int,
    grid_group_sizes,
    device_type: str,
    backend: str,
    env_per_rank,
    dtype: torch.dtype,
    W: int,
    # Inputs on host (global tensors)
    feats_global_host: dict[str, torch.Tensor],
):
    """Parallel worker function for testing pack_atom_features."""
    monkeypatch = pytest.MonkeyPatch()
    if env_per_rank is not None:
        for var_name, value in env_per_rank.items():
            if value == "<INPUT_RANK>":
                monkeypatch.setenv(var_name, f"{rank}")
                continue
            monkeypatch.setenv(var_name, value)

    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    manager = DistributedManager()

    # ========================================================================
    # Distribute atom features using pad_and_scatter_atom_features_dtensor
    # ========================================================================
    # This follows the pattern from test_dtensor_model.py:
    # 1. Each dp rank gets one sample from the batch
    # 2. pad_and_scatter_atom_features_dtensor distributes within cp group
    # 3. Results are collated to form full DTensors with dp+cp sharding

    size_batch = feats_global_host["atom_pad_mask"].shape[0]
    assert size_batch == len(manager.group_ranks["dp"]), "size_batch must equal number of dp ranks"
    size_batch_per_dp = size_batch // len(manager.group_ranks["dp"])
    rank_dp = manager.group_rank["dp"]
    i_sample_begin = rank_dp * size_batch_per_dp

    # Prepare inputs for pad_and_scatter_atom_features_dtensor
    # Only cp rank 0 provides the input (it gets scattered to all cp ranks)
    if manager.group_rank["cp"] == 0:
        inputs = {
            k: v[i_sample_begin].to(device=manager.device, dtype=dtype if v.dtype.is_floating_point else v.dtype)
            for k, v in feats_global_host.items()
        }
    else:
        inputs = None

    # Distribute atom features using pad_and_scatter_atom_features_dtensor
    feats_atom_dtensor_cp = pad_and_scatter_atom_features_dtensor(
        inputs,
        _placements_cp_atom_features,
        manager.group["cp"],
        manager.group_ranks["cp"][0],
        manager.device_mesh_subgroups["cp_axis_0", "cp_axis_1"],
    )

    # Collate along batch dimension: per-dp rank atom features -> single DTensor with dp+cp sharding
    feats_atom_shape_stride_global = {
        k: tuple(
            map(
                tuple,
                compute_global_tensor_info(
                    v.to_local().unsqueeze(0), manager.device_mesh_subgroups, _placements_atom_features[k]
                ),
            )
        )
        for k, v in feats_atom_dtensor_cp.items()
    }
    feats_dt = {
        k: DTensor.from_local(
            v.to_local().unsqueeze(0),
            manager.device_mesh_subgroups,
            _placements_atom_features[k],
            shape=feats_atom_shape_stride_global[k][0],
            stride=feats_atom_shape_stride_global[k][1],
        )
        for k, v in feats_atom_dtensor_cp.items()
    }

    # ========================================================================
    # Pack and pad atom features using pack_atom_features
    # This removes per-shard trailing padding from pad_and_scatter_atom_features_dtensor
    # and creates a packed DTensor with global trailing padding (multiple of W * size_cp)
    # ========================================================================
    # Suppress the expected warning about atom_to_token not being packed
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="pack_atom_features: 'atom_to_token'")
        feats_dt_packed = pack_atom_features(feats_dt, set(feats_dt.keys()), W)

    # ========================================================================
    # Verify consistency between packed features and serial reference
    # Packed features should match serial except for trailing padding along N_atoms axis
    # ========================================================================
    N_atoms_actual = feats_global_host["atom_pad_mask"].shape[1]

    # Check that packed features match serial in valid region
    # Note: atom_to_token is excluded because pack_atom_features only outputs
    # atom_to_token_ids_global (global indices), not the original one-hot matrix.
    # The global indices are computed via shardwise_argmax + shardwise_offset before packing.
    atom_feature_keys_to_check = [
        "ref_pos",
        "ref_space_uid",
        "ref_charge",
        "ref_element",
        "ref_atom_name_chars",
        "atom_pad_mask",
    ]
    for key in atom_feature_keys_to_check:
        if key not in feats_dt_packed:
            continue
        packed_full = feats_dt_packed[key].full_tensor()
        serial_ref = feats_global_host[key].to(device=manager.device, dtype=packed_full.dtype)
        assert_tensors_close_with_pad(
            packed_full,
            serial_ref,
            axis=1,
            pad_val=0,
            msg=lambda m, k=key: f"Packed feature {k} mismatch: {m}",
        )

    # Verify atom_to_token_ids_global is present and has correct shape
    assert "atom_to_token_ids_global" in feats_dt_packed, "atom_to_token_ids_global should be in packed features"
    atom_to_token_ids_global_packed = feats_dt_packed["atom_to_token_ids_global"]
    # Shape should be (B, N_atoms_packed) where N_atoms_packed >= N_atoms_actual
    assert atom_to_token_ids_global_packed.shape[0] == size_batch
    assert atom_to_token_ids_global_packed.shape[1] >= N_atoms_actual

    # 'atom_to_token' is returned as it is with a different feature name
    assert "atom_to_token_local_onehot" in feats_dt_packed, "atom_to_token_local_onehot should be in packed features"
    atom_to_token_local_onehot_packed = feats_dt_packed["atom_to_token_local_onehot"]
    atom_to_token_local_onehot_packed_full = atom_to_token_local_onehot_packed.full_tensor()
    atom_to_token_expected_full = feats_dt["atom_to_token"].full_tensor()
    assert_tensors_identical(atom_to_token_local_onehot_packed_full, atom_to_token_expected_full)

    # Verify the global indices match the serial reference in valid region
    # Serial atom_to_token is one-hot (B, N_atoms, N_tokens), extract indices via argmax
    atom_to_token_serial = feats_global_host["atom_to_token"].to(device=manager.device)
    atom_to_token_ids_serial = atom_to_token_serial.argmax(dim=-1)  # (B, N_atoms)
    atom_to_token_ids_global_packed_full = atom_to_token_ids_global_packed.full_tensor()
    assert_tensors_close_with_pad(
        atom_to_token_ids_global_packed_full,
        atom_to_token_ids_serial,
        axis=1,
        pad_val=0,
        msg=lambda m: f"atom_to_token_ids_global mismatch: {m}",
    )

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.slow
@pytest.mark.parametrize(
    "setup_env",
    [
        ((2, (2, 2)), True, "cuda", "ENV"),
    ],
    indirect=("setup_env",),
    ids=lambda x: f"dp:{x[0][0]}, cp:{x[0][1]}, specify_method:{x[1]}, device_type:{x[2]}, method_init:{x[3]}",
)
@pytest.mark.parametrize(
    "n_atoms_per_token_range",
    [(8, 20)],  # can use (1, 1) for debugging purpose
    ids=lambda x: f"atoms_per_token:{x[0]}-{x[1]}",
)
def test_pack_and_pad_atom_features(setup_env, n_atoms_per_token_range: tuple[int, int]):
    """Test pack_atom_features function.

    Verifies that:
    1. Packed features match serial reference in valid region
    2. Variable atoms per token are handled correctly
    3. atom_to_token_ids_global contains correct global token indices
    """
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    # Configuration
    dtype = torch.float32
    seed = 42
    seed_by_rank(0, seed=seed)

    size_cp = grid_group_sizes["cp"][0]
    B = 1 * grid_group_sizes["dp"]  # batch size per rank = 1

    # Test parameters
    n_atoms_per_token_min, n_atoms_per_token_max = n_atoms_per_token_range
    W = 32  # atoms per window for queries
    # N_tokens must be divisible by size_cp for even token sharding
    N_tokens = 1000 * size_cp
    # With max atoms per token, N_atoms = N_tokens * n_atoms_per_token_max
    N_atoms = N_tokens * n_atoms_per_token_max
    N_msa = 1  # minimal MSA

    val_init_min_max = (-0.5, 0.5)

    # Verify constraints
    assert N_tokens % size_cp == 0, f"N_tokens ({N_tokens}) must be divisible by size_cp ({size_cp})"

    # ========================================================================
    # Generate features using random_features
    # This subset of features is for AtomAttentionEncoder usage
    # ========================================================================
    selected_keys = [
        "atom_pad_mask",
        "ref_pos",
        "ref_space_uid",
        "ref_charge",
        "ref_element",
        "ref_atom_name_chars",
        "atom_to_token",
        "atom_counts_per_token",  # Required by pad_and_scatter_atom_features_dtensor
    ]

    feats = random_features(
        size_batch=B,
        n_tokens=N_tokens,
        n_atoms=N_atoms,
        n_msa=N_msa,
        atom_counts_per_token_range=(n_atoms_per_token_min, n_atoms_per_token_max),
        device=torch.device(device_type),
        float_value_range=val_init_min_max,
        selected_keys=selected_keys,
    )

    # Convert float64 to float32 for consistency
    feats = {k: v.to(dtype=dtype) if v.dtype == torch.float64 else v for k, v in feats.items()}

    # Prepare inputs for distributed test
    feats_global_host = {k: v.detach().cpu() for k, v in feats.items()}

    # Launch multiprocess test
    spawn_multiprocessing(
        parallel_assert_pack_and_pad_atom_features,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        dtype,
        W,
        feats_global_host,
    )
