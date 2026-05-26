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
from typing import Optional

import pytest
import torch
from torch.distributed.tensor import Replicate, Shard, distribute_tensor

from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.layers.clip import clip
from boltz.testing.utils import assert_tensors_identical, spawn_multiprocessing


def serial_clip(tensor: torch.Tensor, min_val: Optional[float] = None, max_val: Optional[float] = None) -> torch.Tensor:
    """Serial implementation of clip operation for comparison."""
    return torch.clip(tensor, min=min_val, max=max_val)


def parallel_assert_clip(
    rank,
    grid_group_sizes,
    device_type,
    backend,
    env_per_rank,
    placements,
    min_val,
    max_val,
    input_global_host,
    output_expected_global_host,
    d_output_expected_global_host,
    d_input_expected_global_host,
):
    monkeypatch = pytest.MonkeyPatch()
    if env_per_rank is not None:
        for var_name, value in env_per_rank.items():
            if value == "<INPUT_RANK>":
                monkeypatch.setenv(var_name, f"{rank}")
                continue
            monkeypatch.setenv(var_name, value)

    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    manager = DistributedManager()

    # Distribute input tensor
    input_dtensor = distribute_tensor(
        input_global_host.to(manager.device), device_mesh=manager.device_mesh_subgroups, placements=placements
    ).requires_grad_(True)

    # Distribute expected outputs
    d_output_expected_dtensor = distribute_tensor(
        d_output_expected_global_host.to(manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements,
    )
    output_expected_dtensor = distribute_tensor(
        output_expected_global_host.to(manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements,
        src_data_rank=None,
    )
    d_input_expected_dtensor = distribute_tensor(
        d_input_expected_global_host.to(manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements,
        src_data_rank=None,
    )

    # Create copy to verify input isn't modified
    input_dtensor_copy = input_dtensor.detach().clone().requires_grad_(True)

    # Forward pass
    output_dtensor_result = clip(input_dtensor, min_val=min_val, max_val=max_val)

    # Verify input wasn't modified
    assert_tensors_identical(
        input_dtensor_copy.to_local(), input_dtensor.to_local(), check_grad=False, check_grad_fn=False
    )

    # Test forward pass results
    assert (
        output_dtensor_result.placements == placements
    ), f"placements: {placements}, output_dtensor_result.placements: {output_dtensor_result.placements}"
    assert (
        output_dtensor_result.shape == output_expected_dtensor.shape
    ), f"Output shape mismatch: {output_dtensor_result.shape} != {output_expected_dtensor.shape}"
    assert (
        output_dtensor_result.stride() == output_expected_dtensor.stride()
    ), f"Output stride mismatch: {output_dtensor_result.stride()} != {output_expected_dtensor.stride()}"
    torch.testing.assert_close(output_dtensor_result.to_local(), output_expected_dtensor.to_local())

    # Backward pass
    d_output_expected_dtensor_copy = d_output_expected_dtensor.detach().clone()
    output_dtensor_result.backward(d_output_expected_dtensor)

    # Verify upstream gradient wasn't modified
    assert_tensors_identical(d_output_expected_dtensor_copy.to_local(), d_output_expected_dtensor.to_local())

    # Test input gradient
    assert (
        input_dtensor.grad.placements == placements
    ), f"placements: {placements}, input_dtensor.grad.placements: {input_dtensor.grad.placements}"
    assert (
        input_dtensor.grad.shape == d_input_expected_dtensor.shape
    ), f"Input gradient shape mismatch: {input_dtensor.grad.shape} != {d_input_expected_dtensor.shape}"
    assert (
        input_dtensor.grad.stride() == d_input_expected_dtensor.stride()
    ), f"Input gradient stride mismatch: {input_dtensor.grad.stride()} != {d_input_expected_dtensor.stride()}"
    torch.testing.assert_close(input_dtensor.grad.to_local(), d_input_expected_dtensor.to_local())

    # Test full tensor gathering - verify distributed results match serial results
    output_global_result_host = output_dtensor_result.full_tensor().cpu()
    d_input_global_result_host = input_dtensor.grad.full_tensor().cpu()

    # Verify full tensors match expected results
    torch.testing.assert_close(output_global_result_host, output_expected_global_host)
    torch.testing.assert_close(d_input_global_result_host, d_input_expected_global_host)

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env",
    [
        ((1, (2, 2)), True, "cuda", "ENV"),
        ((2, (2, 2)), True, "cuda", "ENV"),
        ((1, (3, 3)), True, "cuda", "ENV"),
        ((1, (3, 3)), True, "cpu", "ENV"),
    ],
    indirect=("setup_env",),
    ids=lambda x: f"dp:{x[0][0]}, cp:{x[0][1]}, specify_method:{x[1]}, device_type:{x[2]}, method_init:{x[3]}",
)
@pytest.mark.parametrize(
    "placements", [(Shard(0), Shard(1), Shard(2)), (Shard(0), Shard(1), Replicate())], ids=["shard", "replicate"]
)
@pytest.mark.parametrize(
    "clip_params",
    [
        (0.0, None),  # min_val only (ReLU-like)
        (None, 5.0),  # max_val only
        (-2.0, 1.0),  # both min and max
    ],
    ids=["min_only", "max_only", "min_max"],
)
def test_clip_parallel(setup_env, placements, clip_params):
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env
    min_val, max_val = clip_params

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    size_ring = grid_group_sizes["cp"][0]
    B = 2 * grid_group_sizes["dp"]
    N = size_ring * 100  # Number of tokens
    D = 32  # Hidden dimension

    seed = 42
    rng = torch.Generator(device=device_type)
    rng.manual_seed(seed)

    # Create input tensor with values that will test clipping
    input_global = torch.empty((B, N, N, D), requires_grad=True, device=device_type)
    with torch.no_grad():
        input_global.uniform_(-10, 10, generator=rng)  # Wide range to ensure clipping occurs

    # Run serial forward pass
    input_global_host = input_global.detach().clone().cpu()
    output_expected_global = serial_clip(input_global, min_val=min_val, max_val=max_val)
    output_expected_global_host = output_expected_global.detach().clone().cpu()

    # Create upstream gradient and run backward pass
    d_output_expected_global = torch.rand_like(output_expected_global)
    d_output_expected_global_host = d_output_expected_global.detach().clone().cpu()
    output_expected_global.backward(d_output_expected_global)

    # Launch parallel test across all processes
    spawn_multiprocessing(
        parallel_assert_clip,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        placements,
        min_val,
        max_val,
        input_global_host,
        output_expected_global_host,
        d_output_expected_global_host,
        input_global.grad.detach().clone().cpu(),
    )


if __name__ == "__main__":
    unittest.main()
