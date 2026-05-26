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


import itertools
import os
from typing import Dict, Optional

import pytest
import torch

from boltz.distributed.manager import DistributedManager


def test_manager_singleton(monkeypatch):
    # Test distributed manager singleton functions as expected
    monkeypatch.setenv("MASTER_ADDR", "localhost")
    monkeypatch.setenv("MASTER_PORT", "45678")
    monkeypatch.setenv("RANK", "0")
    monkeypatch.setenv("WORLD_SIZE", "1")
    DistributedManager.initialize({"dp": 1, "cp": 1}, "cpu", "gloo")

    manager_1 = DistributedManager()
    manager_1.random_property = "random_string"
    manager_2 = DistributedManager()

    # Compare attributes
    for attr in manager_1.__dict__.keys():
        assert getattr(manager_1, attr) == getattr(manager_2, attr)
    assert manager_1.random_property == manager_2.random_property
    DistributedManager.cleanup()


def create_manager_and_assert(
    rank_expected: int,
    world_size_expected: int,
    grid_group_sizes_expected: Dict[str, int],
    device_type_expected: str,
    backend_expected: str,
    method_init_expected: str,
    env_map: Optional[Dict[str, str]] = None,
):
    monkeypatch = pytest.MonkeyPatch()
    if env_map is not None:
        for var_name, value in env_map.items():
            if value == "<INPUT_RANK>":
                monkeypatch.setenv(var_name, f"{rank_expected}")
                continue
            monkeypatch.setenv(var_name, value)

    DistributedManager.initialize(grid_group_sizes_expected, device_type=device_type_expected, backend=backend_expected)
    manager = DistributedManager()
    assert (
        manager.has_dist == torch.distributed.is_available()
    ), "DistributedManager.has_dist inconsistent with torch.distributed's availability"
    assert manager.rank == rank_expected
    assert manager.world_size == world_size_expected
    assert manager.group["world"] == torch.distributed.group.WORLD
    assert manager.group_rank["world"] == rank_expected
    assert manager.group_ranks["world"] == list(range(world_size_expected))
    # by default, the underlying DeviceMesh should use layout-right for the grid groups
    layoutMap = manager.layout_device_mesh
    grid_coords_expected = layoutMap.unravel(rank_expected)
    has_subgroups = False
    for i_group, (name_group, size_group) in enumerate(grid_group_sizes_expected.items()):
        if isinstance(size_group, tuple) and all(isinstance(size_group_i, int) for size_group_i in size_group):
            has_subgroups = True
            layoutMap_subgroups = manager.layout_device_mesh_subgroups
            grid_coords_expected_subgroups = layoutMap_subgroups.unravel(rank_expected)
            # check if the layout of the subgroups is correct
            slices_subgroup = list(grid_coords_expected_subgroups)
            for i_subgroup, size_subgroup in enumerate(size_group):
                name_subgroup = f"{name_group}_axis_{i_subgroup}"
                # check if the subgroups' ranks are set consistently with the layout
                assert len(manager.group_ranks[name_subgroup]) == size_subgroup
                assert manager.group_rank[name_subgroup] == grid_coords_expected_subgroups[i_group + i_subgroup]
                # check if the parent group is mapped correctly to the subgroups
                assert manager.subgroups[name_group][i_subgroup] is manager.group[name_subgroup]
                assert manager.subgroups_ranks[name_group][i_subgroup] == manager.group_ranks[name_subgroup]
                assert manager.subgroups_rank[name_group][i_subgroup] == manager.group_rank[name_subgroup]
                slices_subgroup[i_group + i_subgroup] = slice(None)
            # check if the layout of the subgroups is correct
            layoutMap_subgroup = layoutMap_subgroups[*slices_subgroup]
            assert manager.layout_subgroups[name_group].shape == layoutMap_subgroup.shape
            assert manager.layout_subgroups[name_group].strides == layoutMap_subgroup.strides
            assert manager.layout_subgroups[name_group].offset == 0
        elif isinstance(size_group, int):
            assert len(manager.group_ranks[name_group]) == size_group
            assert manager.group_rank[name_group] == grid_coords_expected[i_group]
        else:
            raise ValueError(f"Invalid group size type: {type(size_group)}")
    assert manager.has_subgroups == has_subgroups

    assert manager.backend == backend_expected
    assert manager.device.type == device_type_expected
    assert manager.method_init == method_init_expected
    DistributedManager.cleanup()

    monkeypatch.undo()


def create_default_manager_can_raise(
    rank_expected: int,
    grid_group_sizes_expected: Dict[str, int],
    device_type_expected: str,
    backend_expected: str,
    env_map: Optional[Dict[str, str]] = None,
):
    monkeypatch = pytest.MonkeyPatch()
    if env_map is not None:
        for var_name, value in env_map.items():
            if value == "<INPUT_RANK>":
                monkeypatch.setenv(var_name, f"{rank_expected}")
                continue
            monkeypatch.setenv(var_name, value)
    if "BOLTZ_DISTRIBUTED_INIT_METHOD" in os.environ:
        # setting BOLTZ_DISTRIBUTED_INIT_METHOD without the other
        # relevant env vars for world_size and rank etc will trigger
        # a RuntimeError
        with pytest.raises(RuntimeError):
            DistributedManager.initialize(
                grid_group_sizes_expected, device_type=device_type_expected, backend=backend_expected
            )
    else:
        # default initialization should happen
        DistributedManager.initialize(
            grid_group_sizes_expected, device_type=device_type_expected, backend=backend_expected
        )
        manager = DistributedManager()
        assert manager.initialized
        assert not manager.has_dist
        assert DistributedManager().rank == 0
        assert DistributedManager().world_size == 1
        assert DistributedManager().local_rank == 0
        assert DistributedManager().device == torch.device("cpu")
        assert DistributedManager().backend is None
        assert DistributedManager().method_init is None
        assert DistributedManager().group == {}
        assert DistributedManager().group_rank == {}
        assert DistributedManager().group_ranks == {}
        assert DistributedManager().device_mesh is None
        assert DistributedManager().device_mesh_subgroups is None
        assert DistributedManager().layout_device_mesh is None
        assert DistributedManager().layout_device_mesh_subgroups is None
        assert DistributedManager().has_subgroups is False
        assert DistributedManager().subgroups == {}
        assert DistributedManager().subgroups_ranks == {}
        assert DistributedManager().subgroups_rank == {}
        assert DistributedManager().layout_subgroups == {}

    monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env",
    itertools.product([(1, 2)], [False], ["cpu", "cuda"], [None]),
    indirect=("setup_env",),
    ids=lambda x: f"dp={x[0][0]}, cp={x[0][1]}, specify_method={x[1]}, device_type={x[2]} method_init={x[3]}",
)
def test_manager_default(setup_env):
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env
    # device_type and backend don't matter since there is no valid set of distributed environment variables
    torch.multiprocessing.set_start_method("spawn", force=True)
    torch.multiprocessing.spawn(
        fn=create_default_manager_can_raise,
        args=(grid_group_sizes, device_type, backend, env_per_rank),
        nprocs=world_size,
        join=True,
    )


@pytest.mark.parametrize(
    "setup_env",
    itertools.product([(1, 1), (1, 2), (2, (2, 2)), (1, (4, 4))], [True, False], ["cpu", "cuda"], ["ENV", "SLURM"]),
    indirect=("setup_env",),
    ids=lambda x: f"dp={x[0][0]}, cp={x[0][1]}, specify_method={x[1]}, device_type={x[2]} method_init={x[3]}",
)
def test_manager(setup_env):
    grid_group_sizes, world_size, device_type, backend, method_init, env_per_rank = setup_env
    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip("skip cuda test because torch.cuda.device_count() != world_size")

    torch.multiprocessing.set_start_method("spawn", force=True)
    torch.multiprocessing.spawn(
        fn=create_manager_and_assert,
        args=(world_size, grid_group_sizes, device_type, backend, method_init, env_per_rank),
        nprocs=world_size,
        join=True,
    )


def create_manager_and_group(
    rank: int,
    grid_group_sizes: Dict[str, int],
    device_type: str,
    backend: str,
    env_map: Optional[Dict[str, str]] = None,
):
    if env_map is not None:
        for var_name, value in env_map.items():
            if value == "<INPUT_RANK>":
                os.environ[var_name] = f"{rank}"
                continue
            os.environ[var_name] = value
    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    manager = DistributedManager()

    assert "cp" in grid_group_sizes, "Must have group 'cp' in the input grid_group_sizes"
    name_group = "new_group_test"
    DistributedManager.create_group(name_group, manager.group_ranks["cp"], use_local_synchronization=True)
    assert name_group in manager.group
    assert manager.group_ranks[name_group] == manager.group_ranks["cp"]
    assert manager.group_rank[name_group] == manager.group_rank["cp"]
    DistributedManager.cleanup()


@pytest.mark.parametrize(
    "setup_env",
    itertools.product([(2, 4)], [False], ["cpu"], ["ENV"]),
    indirect=("setup_env",),
    ids=lambda x: f"dp={x[0][0]}, cp={x[0][1]}, specify_method={x[1]}, device_type={x[2]} method_init={x[3]}",
)
def test_manager_create_group(setup_env):
    grid_group_sizes, world_size, device_type, backend, method_init, env_per_rank = setup_env
    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip("skip cuda test because torch.cuda.device_count() != world_size")

    torch.multiprocessing.set_start_method("spawn", force=True)
    torch.multiprocessing.spawn(
        fn=create_manager_and_group,
        args=(grid_group_sizes, device_type, backend, env_per_rank),
        nprocs=world_size,
        join=True,
    )
