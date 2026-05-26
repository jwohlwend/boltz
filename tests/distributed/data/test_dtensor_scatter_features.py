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


from typing import Dict, Optional

import pytest
import torch
from torch.distributed.tensor import Replicate, Shard

from boltz.distributed.data.feature.featurizer import pad_and_scatter_atom_features_dtensor
from boltz.distributed.data.feature.featurizer_utils import remap_atom_indices_repad
from boltz.distributed.manager import DistributedManager
from boltz.testing.utils import seed_by_rank, spawn_multiprocessing


@pytest.mark.parametrize(
    "old_stride,new_stride,n_shards",
    [
        (5, 8, 3),
        (4, 4, 2),
        (6, 10, 4),
    ],
)
def test_remap_atom_indices_repad(old_stride, new_stride, n_shards):
    """Verify remap_atom_indices_repad correctly remaps padded atom indices."""
    # Build indices that cover multiple shards:
    # for each shard s, place a few valid offsets within [0, old_stride).
    indices = []
    expected = []
    for s in range(n_shards):
        for offset in [0, 1, old_stride - 1]:
            indices.append(s * old_stride + offset)
            expected.append(s * new_stride + offset)

    indices_t = torch.tensor(indices, dtype=torch.int64)
    expected_t = torch.tensor(expected, dtype=torch.int64)

    result = remap_atom_indices_repad(indices_t, old_stride, new_stride)
    torch.testing.assert_close(result, expected_t)


def _make_sample_features(n_tokens, atom_counts_per_token, device):
    """Build minimal feature/placement dicts for a single sample."""
    total_atoms = atom_counts_per_token.sum().item()
    features = {
        "atom_counts_per_token": atom_counts_per_token.to(device),
        "atom_pad_mask": torch.ones(total_atoms, device=device),
        "frames_idx": torch.randint(0, total_atoms, (1, n_tokens, 3), device=device, dtype=torch.int64),
    }
    placements = {
        "atom_counts_per_token": (Shard(0), Replicate()),
        "atom_pad_mask": (Shard(0), Replicate()),
        "frames_idx": (Shard(1), Replicate()),
    }
    return features, placements


def parallel_assert_collate_dtensor_atom_index_remap(
    rank: int,
    grid_group_sizes: Dict[str, int],
    device_type: str,
    backend: str,
    env_map: Optional[Dict[str, str]] = None,
):
    """Verify CollateDTensor remaps atom-index features when samples differ in max_atoms_per_shard."""
    from boltz.distributed.data.utils import CollateDTensor, map_subgroup_mesh_to_cpu

    monkeypatch = pytest.MonkeyPatch()
    if env_map is not None:
        for var_name, value in env_map.items():
            if value == "<INPUT_RANK>":
                monkeypatch.setenv(var_name, f"{rank}")
                continue
            monkeypatch.setenv(var_name, value)

    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    manager = DistributedManager()

    seed_by_rank(0, seed=99)

    cp_submesh = manager.device_mesh_subgroups["cp_axis_0", "cp_axis_1"]
    n_cp = cp_submesh.shape[0]
    cp_group = manager.group["cp"]
    src_rank = cp_submesh.mesh.flatten()[0].item()
    is_src = manager.group_rank["world"] == src_rank
    device = manager.device

    # Two samples with different atom counts → different max_atoms_per_shard
    n_tokens = 4 * n_cp
    counts_a = torch.tensor([2] * n_tokens, device=device)
    counts_b = torch.tensor([4] * n_tokens, device=device)

    feats_a, place_a = _make_sample_features(n_tokens, counts_a, device)
    feats_b, place_b = _make_sample_features(n_tokens, counts_b, device)

    dtensors_a = pad_and_scatter_atom_features_dtensor(
        feats_a if is_src else None, place_a, cp_group, src_rank, cp_submesh
    )
    dtensors_b = pad_and_scatter_atom_features_dtensor(
        feats_b if is_src else None, place_b, cp_group, src_rank, cp_submesh
    )

    # Record per-sample atom dim before collation
    atoms_per_shard_a = dtensors_a["atom_pad_mask"].to_local().shape[0]
    atoms_per_shard_b = dtensors_b["atom_pad_mask"].to_local().shape[0]
    assert atoms_per_shard_a != atoms_per_shard_b, "samples must differ in max_atoms_per_shard"
    smaller, larger = sorted([atoms_per_shard_a, atoms_per_shard_b])

    # Save frames_idx local values before collation
    fidx_local_a = dtensors_a["frames_idx"].to_local().clone()
    fidx_local_b = dtensors_b["frames_idx"].to_local().clone()

    # Collate
    dp_cp_mesh = map_subgroup_mesh_to_cpu(manager)
    collator = CollateDTensor(dp_cp_mesh)
    batch = collator([dtensors_a, dtensors_b])

    # After collation, atom dim should be the larger of the two
    final_atoms_per_shard = batch["atom_pad_mask"].to_local().shape[1]
    assert final_atoms_per_shard == larger

    # frames_idx for the smaller sample should have been remapped
    fidx_batch = batch["frames_idx"].to_local()
    fidx_collated_a = fidx_batch[0]
    fidx_collated_b = fidx_batch[1]

    expected_a = remap_atom_indices_repad(fidx_local_a, atoms_per_shard_a, final_atoms_per_shard)
    expected_b = remap_atom_indices_repad(fidx_local_b, atoms_per_shard_b, final_atoms_per_shard)

    torch.testing.assert_close(fidx_collated_a, expected_a)
    torch.testing.assert_close(fidx_collated_b, expected_b)

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env",
    [
        ((1, (2, 2)), True, "cpu", "ENV"),
    ],
    indirect=("setup_env",),
    ids=lambda x: f"dp:{x[0][0]}, cp:{x[0][1]}, device_type:{x[2]}",
)
def test_collate_dtensor_atom_index_remap(setup_env):
    """Verify CollateDTensor remaps atom-index features when batch samples have different atom counts."""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    spawn_multiprocessing(
        parallel_assert_collate_dtensor_atom_index_remap,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
    )
