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

"""Tests the simple function raise_if_incorrect_dtensor_metadata_args

Verification requirements

    V1a: TypeError not raised when object is a DTensor
    V1b: TypeError raised when object is not a DTensor
    V2a: ValueError not raised when shape matches expected_shape
    V2b: ValueError raised when shape doesn't match expected_shape

    V3a: ValueError raised when device_mesh (d) match expected_device_mesh
    V3b: ValueError raised when device_mesh (d) match expected_device_mesh

    V4a: ValueError raised when placements don't match expected_placements
    V4a: ValueError raised when placements don't match expected_placements

    V5a: ValueError not raised when check_for_partial_placements=True and Partial placement absent
    V5b: ValueError raised when check_for_partial_placements=True and Partial placement present

Implementation status
    V1a/b: done
    V2a/b: done
    V3a/b: done
    V4a/b: done
    V5a/b: done
"""

import itertools
from copy import deepcopy
from typing import Any, Dict, OrderedDict

import pytest
import torch
from torch import Tensor
from torch.distributed import DeviceMesh
from torch.distributed.tensor import DTensor, Partial, Replicate, Shard, distribute_tensor
from torch.distributed.tensor._utils import compute_global_tensor_info

from boltz.distributed.manager import DistributedManager, _GridGroupSizesType
from boltz.distributed.model.layers.dtensor_metadata_tools import raise_if_incorrect_dtensor_metadata_args
from boltz.testing.utils import (
    seed_by_rank,
    skip_if_cuda_not_avail_or_device_count_less_than_word_size,
    spawn_multiprocessing,
)

SEED = 42


def parallel_assert_raise_if_incorrect_dtensor_metadata_args(
    rank: int,
    input_example: OrderedDict[str, Tensor],
    fn_kwargs: dict,  # noqa
    grid_group_sizes: _GridGroupSizesType,
    device_type: str,
    backend: str,
    env_per_rank: Dict[str, str],
):
    monkeypatch = pytest.MonkeyPatch()
    if env_per_rank is not None:
        for var_name, value in env_per_rank.items():
            if value == "<INPUT_RANK>":
                monkeypatch.setenv(var_name, f"{rank}")
                continue
            monkeypatch.setenv(var_name, value)
    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    dist_manager = DistributedManager()

    # Non-boillerplate
    placements_for_single_rep_nonparam = (Shard(0), Shard(1), Replicate())
    input_meta: OrderedDict[str, Dict[str, Any]] = OrderedDict(
        [
            ("x", dict(placements=placements_for_single_rep_nonparam, requires_grad=False)),  # noqa C408
            ("y", dict(placements=placements_for_single_rep_nonparam, requires_grad=False)),  # noqa C408
        ]
    )

    # -------------------------------------------------------------
    # Move inputs and ref outputs to device
    # --------------------------------------------------------------
    input_example_device = OrderedDict(
        [
            (k, v.detach().to(dist_manager.device, copy=True) if isinstance(v, Tensor) else deepcopy(v))
            for (k, v) in input_example.items()
        ]
    )
    # -----------------------------------------------------
    # Create input DTensors
    #   - after move to device
    # ----------------------------------------------------
    input_example_as_dtensor = OrderedDict()
    for name, meta in input_meta.items():
        input_example_as_dtensor[name] = distribute_tensor(
            input_example_device[name], dist_manager.device_mesh_subgroups, meta["placements"]
        ).requires_grad_(meta["requires_grad"])

    # --------------------------------------
    # Verifications
    # --------------------------------------

    # V1a: TypeError not raised if inputs is a DTensor
    raise_if_incorrect_dtensor_metadata_args(
        dtensor_instance=input_example_as_dtensor["x"],
        dtensor_name="x",
    )
    # V1b: TypeError raised if inputs is not a DTensor
    with pytest.raises(TypeError, match="should have type"):
        raise_if_incorrect_dtensor_metadata_args(
            dtensor_instance=input_example_device["x"],
            dtensor_name="x",
        )

    # V2a: ValueError not raised when shapes match
    raise_if_incorrect_dtensor_metadata_args(
        dtensor_instance=input_example_as_dtensor["x"],
        dtensor_name="x",
        expected_shape=input_example_device["y"].shape,
    )
    # V2b: ValueError raised when shape doesn't match
    with pytest.raises(ValueError, match="should have shape"):
        raise_if_incorrect_dtensor_metadata_args(
            dtensor_instance=input_example_as_dtensor["x"],
            dtensor_name="x",
            expected_shape=(
                input_example_device["x"].shape[0] + 1,
                input_example_device["x"].shape[1],
                input_example_device["x"].shape[2],
            ),
        )
    # V3a: ValueError not raised when device_mesh does match
    raise_if_incorrect_dtensor_metadata_args(
        dtensor_instance=input_example_as_dtensor["x"],
        dtensor_name="x",
        expected_device_mesh=input_example_as_dtensor["x"].device_mesh,
    )
    # V3b: ValueError raised when device_mesh doesn't match
    wrong_mesh = DeviceMesh(device_type, torch.arange(2))
    with pytest.raises(ValueError, match="should have device mesh"):
        raise_if_incorrect_dtensor_metadata_args(
            dtensor_instance=input_example_as_dtensor["x"],
            dtensor_name="x",
            expected_device_mesh=wrong_mesh,
        )

    # V4a: ValueError not raised when placements do match
    raise_if_incorrect_dtensor_metadata_args(
        dtensor_instance=input_example_as_dtensor["x"],
        dtensor_name="x",
        expected_placements=input_example_as_dtensor["x"].placements,
    )
    # V4b: ValueError raised when placements don't match
    wrong_placements = (Replicate(), Replicate())
    with pytest.raises(ValueError, match="placements"):
        raise_if_incorrect_dtensor_metadata_args(
            dtensor_instance=input_example_as_dtensor["x"],
            dtensor_name="x",
            expected_placements=wrong_placements,
        )

    # V5a: ValueError not raised placements do not contains Partial
    raise_if_incorrect_dtensor_metadata_args(
        dtensor_instance=input_example_as_dtensor["x"],
        dtensor_name="x",
        check_for_partial_placements=True,
    )
    # ---------------------------------------------------------------------
    # V5b: ValueError raised when Partial placement exists and check_for_partial_placements=True
    #   - Use-case:
    #       - Some custom autograd function f
    #         f.apply(a) -> (x,y)
    #     - Some custom autograd function g
    #         g.apply(x, y) -> z
    #     - g.forward(x, y) might check that x and y have the same placements
    #     but x.placement or y.placement might contain Partial
    # ---------------------------------------------------------------------
    # emulate the forward pass of f
    shape_x_global, stride_x_global = map(
        tuple,
        compute_global_tensor_info(
            input_example_device["x"], dist_manager.device_mesh_subgroups, (Shard(0), Shard(1), Partial())
        ),
    )
    x_with_partial = DTensor.from_local(
        input_example_device["x"],
        device_mesh=dist_manager.device_mesh_subgroups,
        placements=(Shard(0), Shard(1), Partial()),
        shape=shape_x_global,
        stride=stride_x_global,
    )

    shape_y_global, stride_y_global = map(
        tuple,
        compute_global_tensor_info(
            input_example_device["y"], dist_manager.device_mesh_subgroups, (Shard(0), Shard(1), Partial())
        ),
    )
    y_with_partial = DTensor.from_local(
        input_example_device["y"],
        device_mesh=dist_manager.device_mesh_subgroups,
        placements=(Shard(0), Shard(1), Partial()),
        shape=shape_y_global,
        stride=stride_y_global,
    )
    # emulate one version of the metadata check for g.forward()
    with pytest.raises(ValueError, match="placement of type"):
        raise_if_incorrect_dtensor_metadata_args(
            dtensor_instance=x_with_partial,
            dtensor_name="x_with_partial",
        )
        raise_if_incorrect_dtensor_metadata_args(
            dtensor_instance=y_with_partial,
            dtensor_name="y_with_partial",
            expected_placements=x_with_partial.placements,
            check_for_partial_placements=True,
        )

    # emulate a second version of the metadata check for g.forward()
    with pytest.raises(ValueError, match="placement of type"):
        raise_if_incorrect_dtensor_metadata_args(
            dtensor_instance=x_with_partial,
            dtensor_name="x_with_partial",
            check_for_partial_placements=True,
        )
        raise_if_incorrect_dtensor_metadata_args(
            dtensor_instance=y_with_partial,
            dtensor_name="y_with_partial",
            expected_placements=x_with_partial.placements,
        )

    # cleanup
    DistributedManager.cleanup()
    monkeypatch.undo()


def get_example_input(
    len_dim_0: int,
    len_dim_1: int,
    len_dim_2: int,
    seed: int = SEED,
):
    """Generate example input"""
    seed_by_rank(seed)  # specified seed for each rank
    x = torch.randn(len_dim_0, len_dim_1, len_dim_2, requires_grad=False)
    y = torch.randn(len_dim_0, len_dim_1, len_dim_2, requires_grad=False)
    input_example = OrderedDict([("x", x), ("y", y)])
    return input_example


@pytest.mark.parametrize(
    "setup_env",
    itertools.product([(1, (2, 2))], [True], ["cpu"], ["ENV"]),
    indirect=("setup_env",),
    ids=lambda x: f"dp={x[0][0]}, cp={x[0][1]}, specify_method={x[1]}, device_type={x[2]}, method_init={x[3]}",
)
def test_raise_if_incorrect_dtensor_metadata_args(
    setup_env: dict[str, int],
    len_dim_0: int = 2,
    len_dim_1: int = 3,
    len_dim_2: int = 5,
    seed: int = SEED,
):
    """Test raise_if_incorrect_dtensor_metadata_args."""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env
    skip_if_cuda_not_avail_or_device_count_less_than_word_size(device_type, world_size)

    # Get example input and reference output
    input_example = get_example_input(
        len_dim_0=len_dim_0,
        len_dim_1=len_dim_1,
        len_dim_2=len_dim_2,
        seed=seed,
    )

    # Launch parallel test across all processes
    spawn_multiprocessing(
        parallel_assert_raise_if_incorrect_dtensor_metadata_args,
        world_size,
        input_example,
        {},
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
    )
