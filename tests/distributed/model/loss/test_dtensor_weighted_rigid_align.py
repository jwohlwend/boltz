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

"""Tests for DTensor weighted_rigid_align.

Adapted from Boltz-1x CP tests. Verifies that the DTensor weighted_rigid_align
produces identical results to the serial version, and that outputs are
binary-identical across replicate ranks.
"""

from math import isqrt
from typing import Optional

import pytest
import torch
from torch.distributed.tensor import DeviceMesh, DTensor, Replicate, Shard, distribute_tensor

from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.loss.diffusion import weighted_rigid_align as dtensor_weighted_rigid_align
from boltz.model.loss.diffusionv2 import weighted_rigid_align as serial_weighted_rigid_align
from boltz.testing.utils import assert_all_identical, assert_tensors_identical, seed_by_rank, spawn_multiprocessing


def compute_serial_expectation(
    true_coords_global: torch.Tensor,
    pred_coords_global: torch.Tensor,
    weights_global: torch.Tensor,
    mask_global: torch.Tensor,
    device: torch.device,
) -> torch.Tensor:
    """Compute expected result using serial weighted_rigid_align."""
    true_coords_device = true_coords_global.to(device)
    pred_coords_device = pred_coords_global.to(device)
    weights_device = weights_global.to(device)
    mask_device = mask_global.to(device)

    result = serial_weighted_rigid_align(true_coords_device, pred_coords_device, weights_device, mask_device)
    return result.detach().clone()


def compute_dtensor_weighted_rigid_align_with_validation(
    true_coords_global: torch.Tensor,
    pred_coords_global: torch.Tensor,
    weights_global: torch.Tensor,
    mask_global: torch.Tensor,
    device_mesh: DeviceMesh,
    label_test_case: str,
) -> DTensor:
    """Compute DTensor weighted_rigid_align with input validation checks."""
    coords_placements = (Shard(0), Shard(1), Replicate())

    true_coords_dtensor = distribute_tensor(true_coords_global.detach().clone(), device_mesh, coords_placements)
    pred_coords_dtensor = distribute_tensor(pred_coords_global.detach().clone(), device_mesh, coords_placements)
    weights_dtensor = distribute_tensor(weights_global.detach().clone(), device_mesh, coords_placements)
    mask_dtensor = distribute_tensor(mask_global.detach().clone(), device_mesh, coords_placements)

    # Create copies for validation
    true_coords_copy = true_coords_dtensor.detach().clone()
    pred_coords_copy = pred_coords_dtensor.detach().clone()
    weights_copy = weights_dtensor.detach().clone()
    mask_copy = mask_dtensor.detach().clone()

    result_dtensor = dtensor_weighted_rigid_align(
        true_coords_dtensor, pred_coords_dtensor, weights_dtensor, mask_dtensor
    )

    # Verify no change to inputs
    assert_tensors_identical(
        true_coords_dtensor.to_local(), true_coords_copy.to_local(), check_grad=False, check_grad_fn=False
    )
    assert_tensors_identical(
        pred_coords_dtensor.to_local(), pred_coords_copy.to_local(), check_grad=False, check_grad_fn=False
    )
    assert_tensors_identical(weights_dtensor.to_local(), weights_copy.to_local(), check_grad=False, check_grad_fn=False)
    assert_tensors_identical(mask_dtensor.to_local(), mask_copy.to_local(), check_grad=False, check_grad_fn=False)

    # Verify output placements
    assert (
        result_dtensor.placements == coords_placements
    ), f"{label_test_case} output placements mismatch with input placements"

    # Verify binary identical output across the Replicate device_mesh axis
    assert_all_identical(result_dtensor.to_local(), device_mesh.get_group(2))

    return result_dtensor


def parallel_assert_dtensor_weighted_rigid_align(
    rank: int,
    grid_group_sizes: dict[str, int],
    device_type: str,
    backend: str,
    env_map: Optional[dict[str, str]] = None,
):
    monkeypatch = pytest.MonkeyPatch()
    if env_map is not None:
        for var_name, value in env_map.items():
            if value == "<INPUT_RANK>":
                monkeypatch.setenv(var_name, f"{rank}")
                continue
            monkeypatch.setenv(var_name, value)

    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    manager = DistributedManager()

    seed_by_rank(0, seed=42)

    size_cp = len(manager.group_ranks["cp"])
    size_ring = isqrt(size_cp)
    if size_ring * size_ring != size_cp:
        raise ValueError(f"cp group size {size_cp} is not a square int")

    # Random input test case
    size_batch = 2
    n_atoms_per_rank = 6
    n_atoms_padding_per_rank = 2
    n_atoms = size_ring * n_atoms_per_rank

    true_coords_global = torch.randn((size_batch, n_atoms, 3), dtype=torch.float32) * 2.0 + 10.0
    pred_coords_global = torch.randn((size_batch, n_atoms, 3), dtype=torch.float32) * 5.0 + 15.0

    # Mask with padding per rank shard
    mask_global = torch.zeros((size_batch, size_ring, n_atoms_per_rank), dtype=torch.float32)
    mask_global[:, :, :(-n_atoms_padding_per_rank)] = 1.0
    mask_global = mask_global.reshape(size_batch, n_atoms)
    weights_global = mask_global.clone()

    label = "random_input"

    # Serial expectation
    expected_result = compute_serial_expectation(
        true_coords_global, pred_coords_global, weights_global, mask_global, manager.device
    )

    # DTensor result with validation
    result_dtensor = compute_dtensor_weighted_rigid_align_with_validation(
        true_coords_global,
        pred_coords_global,
        weights_global,
        mask_global,
        manager.device_mesh_subgroups,
        label,
    )

    # Compare local shards
    expected_dtensor = distribute_tensor(expected_result, manager.device_mesh_subgroups, result_dtensor.placements)
    torch.testing.assert_close(
        result_dtensor.to_local(),
        expected_dtensor.to_local(),
        msg=lambda m: f"{label} local shard mismatch: {m}",
    )

    # Compare global tensors
    result_global = result_dtensor.full_tensor()
    torch.testing.assert_close(
        result_global,
        expected_result,
        msg=lambda m: f"{label} global result mismatch: {m}",
    )

    DistributedManager.cleanup()
    monkeypatch.undo()


def _build_degenerate_inputs(case: str, device: torch.device):
    """Build inputs that trigger the degenerate-case warnings in DTensor weighted_rigid_align."""
    if case == "scalar_degenerate":
        # total_num_points <= dim (2 <= 3) → scalar warning
        true_coords = torch.randn(1, 2, 3, dtype=torch.float32, device=device) * 2.0 + 10.0
        pred_coords = torch.randn(1, 2, 3, dtype=torch.float32, device=device) * 5.0 + 15.0
        mask = torch.ones(1, 2, dtype=torch.float32, device=device)
        weights = mask.clone()
    elif case == "per_batch_degenerate":
        # batch size 2, 4 points; mask so batch index 1 has only 2 valid (< dim+1=4) → per-batch warning
        true_coords = torch.randn(2, 4, 3, dtype=torch.float32, device=device) * 2.0 + 10.0
        pred_coords = torch.randn(2, 4, 3, dtype=torch.float32, device=device) * 5.0 + 15.0
        mask = torch.ones(2, 4, dtype=torch.float32, device=device)
        mask[1, 2:] = 0.0  # batch 1: only 2 valid points
        weights = mask.clone()
    elif case == "svd_low_rank":
        # 4 points but collinear → covariance rank-deficient → SVD low-rank warning
        true_coords = torch.zeros(1, 4, 3, dtype=torch.float32, device=device)
        true_coords[:, :, 0] = torch.tensor([0.0, 1.0, 2.0, 3.0], device=device)  # on x-axis
        pred_coords = true_coords.clone() + 0.1 * torch.randn(1, 4, 3, device=device)
        mask = torch.ones(1, 4, dtype=torch.float32, device=device)
        weights = mask.clone()
    else:
        raise ValueError(f"Unknown case: {case}")
    return true_coords, pred_coords, weights, mask


def _expected_warning_substrings(case: str):
    """Substrings that must appear in the warning message for the given degenerate case."""
    if case == "scalar_degenerate":
        return ["The size of one of the point clouds is <= dim+1."]
    if case == "per_batch_degenerate":
        return ["[rank_coord:", "Batch indices (subset):"]
    if case == "svd_low_rank":
        return ["[rank_coord:", "Excessively low rank"]
    raise ValueError(f"Unknown case: {case}")


DEGENERATE_CASES = ("scalar_degenerate", "per_batch_degenerate", "svd_low_rank")


def parallel_assert_dtensor_weighted_rigid_align_degenerate(
    rank: int,
    grid_group_sizes: dict,
    device_type: str,
    backend: str,
    env_map: Optional[dict[str, str]],
):
    """Worker: run DTensor weighted_rigid_align on degenerate inputs and assert expected warnings (all cases in one spawn)."""
    monkeypatch = pytest.MonkeyPatch()
    if env_map is not None:
        for var_name, value in env_map.items():
            if value == "<INPUT_RANK>":
                monkeypatch.setenv(var_name, f"{rank}")
                continue
            monkeypatch.setenv(var_name, value)

    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    manager = DistributedManager()

    for case in DEGENERATE_CASES:
        true_coords_global, pred_coords_global, weights_global, mask_global = _build_degenerate_inputs(
            case, manager.device
        )
        coords_placements = (Shard(0), Shard(1), Replicate())
        true_coords_dt = distribute_tensor(true_coords_global, manager.device_mesh_subgroups, coords_placements)
        pred_coords_dt = distribute_tensor(pred_coords_global, manager.device_mesh_subgroups, coords_placements)
        weights_dt = distribute_tensor(weights_global, manager.device_mesh_subgroups, coords_placements)
        mask_dt = distribute_tensor(mask_global, manager.device_mesh_subgroups, coords_placements)

        with pytest.warns(UserWarning) as record:
            dtensor_weighted_rigid_align(true_coords_dt, pred_coords_dt, weights_dt, mask_dt)

        combined_message = " ".join(str(w.message) for w in record)
        for substring in _expected_warning_substrings(case):
            assert (
                substring in combined_message
            ), f"Case {case}: expected substring {substring!r} in warning message(s), got: {combined_message!r}"

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env",
    [
        ((1, (2, 2)), True, "cuda", "ENV"),
        ((2, (2, 2)), True, "cuda", "ENV"),
    ],
    indirect=("setup_env",),
    ids=lambda x: f"dp:{x[0][0]}, cp:{x[0][1]}, specify_method:{x[1]}, device_type:{x[2]}, method_init:{x[3]}",
)
def test_dtensor_weighted_rigid_align(setup_env):
    """Test DTensor weighted_rigid_align vs serial equivalence."""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() < {world_size}")

    spawn_multiprocessing(
        parallel_assert_dtensor_weighted_rigid_align,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
    )


@pytest.mark.parametrize(
    "setup_env",
    [((1, (1, 1)), True, "cuda", "ENV")],
    indirect=("setup_env",),
    ids=lambda x: f"dp:{x[0][0]}, cp:{x[0][1]}",
)
def test_dtensor_weighted_rigid_align_degenerate_warnings(setup_env):
    """Test that DTensor weighted_rigid_align raises expected warnings for degenerate inputs (dp=1, cp=(1,1))."""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() < {world_size}")

    spawn_multiprocessing(
        parallel_assert_dtensor_weighted_rigid_align_degenerate,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
    )
