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

"""Tests for DTensor center_random_augmentation.

Adapted from Boltz-1x CP tests. Verifies centering, random augmentation,
and consistency across CP ranks.
"""

import pytest
import torch
from torch.distributed.device_mesh import DeviceMesh
from torch.distributed.tensor import Replicate, Shard, distribute_tensor

from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.modules.utils import (
    center_random_augmentation,
)
from boltz.model.modules.utils import (
    center_random_augmentation as center_random_augmentation_serial,
)
from boltz.testing.utils import assert_all_identical, assert_tensors_identical, seed_by_rank, spawn_multiprocessing


def assert_center_random_augmentation(rank, payload):
    (
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        centering,
        augmentation,
        s_trans,
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

    device_mesh: DeviceMesh = manager.device_mesh_subgroups
    replicate_group = device_mesh.get_group("cp_axis_1")

    n_samples_per_rank = 2
    batch_size = n_samples_per_rank * device_mesh.get_group("dp").size()
    n_atoms = 5 * device_mesh.get_group("cp_axis_0").size()

    seed_by_rank(0, seed=42)

    atom_coords_gen_global = torch.randn(
        (batch_size, n_atoms, 3),
        dtype=torch.float32,
        device=manager.device,
    )

    atom_coords = distribute_tensor(
        atom_coords_gen_global,
        device_mesh=device_mesh,
        placements=(Shard(0), Shard(1), Replicate()),
    )

    atom_mask_gen_global = torch.randint(
        0,
        2,
        (batch_size, n_atoms),
        dtype=torch.bool,
        device=manager.device,
    )

    atom_mask = distribute_tensor(
        atom_mask_gen_global,
        device_mesh=device_mesh,
        placements=(Shard(0), Shard(1), Replicate()),
    )

    atom_coords_copy = atom_coords.detach().clone()
    atom_mask_copy = atom_mask.detach().clone()

    # all cp ranks should have same seed but
    # different dp ranks should have different seeds
    seed_by_rank(manager.group_rank["dp"])

    results = center_random_augmentation(
        atom_coords,
        atom_mask,
        augmentation=augmentation,
        centering=centering,
        s_trans=s_trans,
        return_roto=augmentation,
    )
    # no modification to the original tensor
    assert_tensors_identical(atom_coords_copy.to_local(), atom_coords.to_local())
    assert_tensors_identical(atom_mask_copy.to_local(), atom_mask.to_local())

    if augmentation:
        atom_coords_augmented, random_R = results
    else:
        atom_coords_augmented = results

    # check if consistent across replicate ranks
    atom_coords_augmented_local = atom_coords_augmented.to_local()
    assert_all_identical(atom_coords_augmented_local, group=replicate_group)

    # check if mean is 0
    if centering and (not augmentation or s_trans == 0.0):
        atom_mask_global = atom_mask.full_tensor().unsqueeze(-1)
        assert_all_identical(atom_mask_global, group=manager.group["world"])

        atom_coords_augmented_full = atom_coords_augmented.full_tensor()
        assert_all_identical(atom_coords_augmented_full, group=manager.group["world"])

        centroids = (atom_coords_augmented_full * atom_mask_global).sum(dim=1) / atom_mask_global.sum(dim=1)

        torch.testing.assert_close(centroids, torch.zeros_like(centroids))

    if augmentation:
        # Verify DTensor augmentation matches serial augmentation on global data.
        # The V2 serial center_random_augmentation does not support return_roto,
        # so we compare coordinates only. The rotation matrix consistency is
        # verified by the assert_all_identical check above.
        seed_by_rank(manager.group_rank["dp"])
        i_sample_begin = manager.group_rank["dp"] * n_samples_per_rank
        i_sample_end = i_sample_begin + n_samples_per_rank
        atom_coords_global = atom_coords.full_tensor()[i_sample_begin:i_sample_end]
        atom_mask_global = atom_mask.full_tensor()[i_sample_begin:i_sample_end]
        atom_coords_augmented_global_expected = center_random_augmentation_serial(
            atom_coords_global,
            atom_mask_global,
            augmentation=augmentation,
            centering=centering,
            s_trans=s_trans,
        )
        assert_all_identical(atom_coords_augmented_global_expected, group=manager.group["cp"])

        atom_coords_augmented_global_result = atom_coords_augmented.full_tensor()[i_sample_begin:i_sample_end]
        torch.testing.assert_close(atom_coords_augmented_global_result, atom_coords_augmented_global_expected)

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env",
    [
        ((2, (2, 2)), True, "cuda", "ENV"),
    ],
    indirect=("setup_env",),
    ids=lambda x: f"dp:{x[0][0]}, cp:{x[0][1]}, specify_method:{x[1]}, device_type:{x[2]}, method_init:{x[3]}",
)
@pytest.mark.parametrize(
    "config",
    [
        (True, False, 1.0),
        (False, False, 1.0),
        (True, True, 2.0),
    ],
    ids=lambda x: f"centering:{x[0]}, augmentation:{x[1]}, s_trans:{x[2]}",
)
def test_center_random_augmentation(
    setup_env,
    config: tuple[bool, bool, float],
):
    """Test DTensor center_random_augmentation vs serial equivalence."""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env
    centering, augmentation, s_trans = config

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() < {world_size}")

    payload = (
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        centering,
        augmentation,
        s_trans,
    )
    spawn_multiprocessing(assert_center_random_augmentation, world_size, payload)
