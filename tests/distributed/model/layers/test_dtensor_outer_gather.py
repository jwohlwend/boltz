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


import math

import pytest
import torch
import torch.nn.functional as F
from torch.distributed.tensor import Replicate, Shard, distribute_tensor

from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.layers.outer_gather import (
    OuterGather,
    compute_interval_overlap,
    distributed_outer_gather,
    get_overlap_from_peers,
    outer_gather,
)
from boltz.distributed.model.layers.shardwise_op import shardwise_argmax
from boltz.testing.utils import assert_tensors_identical, init_tensors_uniform, seed_by_rank, spawn_multiprocessing


@pytest.mark.parametrize("W", [3, 4, 7, 32], ids=lambda x: f"W:{x}")
@pytest.mark.parametrize("H", [5, 8, 16, 128], ids=lambda x: f"H:{x}")
@pytest.mark.parametrize("K", [1, 2, 10, 50], ids=lambda x: f"K:{x}")
@pytest.mark.parametrize("n", [10, 20], ids=lambda x: f"n:{x}")
@pytest.mark.parametrize("m", [7, 11], ids=lambda x: f"m:{x}")
@pytest.mark.parametrize(
    "shape_extra",
    [
        (None,),
        (4, None),
        (None, 3),
        (2, None, 3),
        (2, 3, None),
        (None, 2, 3),
        (2, 3, None, 4),
        (2, None, 3, 4),
        (2, 3, None, 4, 5),
    ],
    ids=lambda x: f"shape_extra:{'x'.join(map(str, x))}",
)
@pytest.mark.parametrize("device", ["cuda"])
def test_outer_gather(W, H, K, n, m, shape_extra, device):
    """
    Test outer_gather against reference einsum.
    """
    assert (
        sum(1 for x in shape_extra if x is None) == 1
    ), "There can be one and only one 'None' element in the shape_extra"
    axis = shape_extra.index(None)
    shape_z = shape_extra[:axis] + (n, m) + shape_extra[axis + 1 :]
    shape_idx_q = shape_extra[:axis] + (K, W)
    shape_idx_k = shape_extra[:axis] + (K, H)

    torch.manual_seed(42)

    # Test both without mask (None) and with random mask
    for use_mask in [False, True]:
        # Generate random masks
        if use_mask:
            idx_q_mask = torch.rand(shape_idx_q, device=device) > 0.5
            idx_k_mask = torch.rand(shape_idx_k, device=device) > 0.5
        else:
            idx_q_mask = None
            idx_k_mask = None

        # Inputs
        z = torch.empty(shape_z, device=device, requires_grad=True)
        init_tensors_uniform([z], low=-0.2, high=0.2)

        # Random one-hots
        idx_q = torch.randint(0, n, shape_idx_q, device=device).requires_grad_(False)
        idx_k = torch.randint(0, m, shape_idx_k, device=device).requires_grad_(False)

        one_hot_q = torch.nn.functional.one_hot(idx_q, num_classes=n).float().requires_grad_(False)
        one_hot_k = torch.nn.functional.one_hot(idx_k, num_classes=m).float().requires_grad_(False)

        # Zero out one-hot at invalid positions when mask is provided
        if idx_q_mask is not None:
            one_hot_q = one_hot_q * idx_q_mask.unsqueeze(-1).float()
        if idx_k_mask is not None:
            one_hot_k = one_hot_k * idx_k_mask.unsqueeze(-1).float()

        n_axes_trailing_z = len(shape_extra[axis + 1 :])
        assert n_axes_trailing_z <= 26, "There can be at most 26 trailing axes for z"

        symbol_trailing_z = "".join(chr(ord("A") + i) for i in range(n_axes_trailing_z))
        out_expected = torch.einsum(
            f"...ij{symbol_trailing_z},...kwi,...khj->...kwh{symbol_trailing_z}", z, one_hot_q, one_hot_k
        )

        # Test Forward
        z_result = z.detach().clone().requires_grad_(True)
        out_result = outer_gather(
            z_result, one_hot_q, one_hot_k, axis=axis, one_hot_q_mask=idx_q_mask, one_hot_k_mask=idx_k_mask
        )

        assert_tensors_identical(out_result, out_expected, check_grad_fn=False, check_stride=False)

        # Test Backward
        grad_out = torch.empty_like(out_expected)
        init_tensors_uniform([grad_out], low=-0.2, high=0.2)

        out_expected.backward(grad_out)

        out_result.backward(grad_out.detach().clone())

        # NOTE: scatter_add involves atomics in the CUDA backend, which leads to
        # abs. error that scales with number of elements so setting tolerance is needed
        torch.testing.assert_close(z.grad, z_result.grad, atol=5e-5, rtol=5e-5)


@pytest.mark.parametrize("device", ["cuda"])
def test_outer_gather_empty_z_all_masked(device):
    """Test outer_gather with empty z and all-False masks."""
    # Test case: z has zero elements along gather dimensions
    shape_z = (2, 0, 0, 3)  # Empty along axis=1 and axis+1=2
    shape_idx_q = (2, 4, 5)  # K=4, W=5
    shape_idx_k = (2, 4, 8)  # K=4, H=8
    axis = 1

    z = torch.zeros(shape_z, device=device, requires_grad=True)
    idx_q = torch.zeros(shape_idx_q, dtype=torch.long, device=device)
    idx_k = torch.zeros(shape_idx_k, dtype=torch.long, device=device)
    idx_q_mask = torch.zeros(shape_idx_q, dtype=torch.bool, device=device)  # All False
    idx_k_mask = torch.zeros(shape_idx_k, dtype=torch.bool, device=device)  # All False

    out = OuterGather.apply(z, idx_q, idx_k, axis, idx_q_mask, idx_k_mask)

    # Output should be zeros with shape (2, 4, 5, 8, 3)
    expected_shape = (2, 4, 5, 8, 3)
    assert out.shape == expected_shape
    assert (out == 0).all()

    # Test backward
    grad_out = torch.randn_like(out)
    out.backward(grad_out)
    assert z.grad.shape == shape_z


@pytest.mark.parametrize("device", ["cuda"])
@pytest.mark.parametrize(
    "q_mask_mode,k_mask_mode,should_raise",
    [
        # Both masks all-False: combined mask is all-False, no error
        ("all_false", "all_false", False),
        # q_mask all-False, k_mask has valid: combined mask is all-False, no error
        ("all_false", "has_valid", False),
        # q_mask has valid, k_mask all-False: combined mask is all-False, no error
        ("has_valid", "all_false", False),
        # q_mask all-False, k_mask is None: combined mask is all-False, no error
        ("all_false", "none", False),
        # q_mask is None, k_mask all-False: combined mask is all-False, no error
        ("none", "all_false", False),
        # Both masks have valid entries: combined mask has valid entries, should raise
        ("has_valid", "has_valid", True),
        # q_mask has valid, k_mask is None (implicitly all-True): should raise
        ("has_valid", "none", True),
        # q_mask is None (implicitly all-True), k_mask has valid: should raise
        ("none", "has_valid", True),
        # Both masks are None (implicitly all-True): should raise
        ("none", "none", True),
    ],
    ids=lambda x: str(x),
)
def test_outer_gather_empty_z_joint_mask_validation(device, q_mask_mode, k_mask_mode, should_raise):
    """Test that OuterGather validates masks jointly (outer-AND) when z is empty.

    The sanity check should only raise an error when the combined mask
    (idx_q_mask outer-AND idx_k_mask) has valid entries. Separate validation
    would incorrectly raise errors in asymmetric cases where one mask is all-False
    but the other has valid entries.
    """
    shape_z = (2, 0, 0, 3)  # Empty z along axis=1 and axis+1=2
    shape_idx_q = (2, 4, 5)  # K=4, W=5
    shape_idx_k = (2, 4, 8)  # K=4, H=8
    axis = 1

    z = torch.zeros(shape_z, device=device, requires_grad=True)
    idx_q = torch.zeros(shape_idx_q, dtype=torch.long, device=device)
    idx_k = torch.zeros(shape_idx_k, dtype=torch.long, device=device)

    # Configure q_mask based on mode
    if q_mask_mode == "none":
        idx_q_mask = None
    elif q_mask_mode == "all_false":
        idx_q_mask = torch.zeros(shape_idx_q, dtype=torch.bool, device=device)
    else:  # "has_valid"
        idx_q_mask = torch.zeros(shape_idx_q, dtype=torch.bool, device=device)
        idx_q_mask[0, 0, 0] = True  # At least one valid entry

    # Configure k_mask based on mode
    if k_mask_mode == "none":
        idx_k_mask = None
    elif k_mask_mode == "all_false":
        idx_k_mask = torch.zeros(shape_idx_k, dtype=torch.bool, device=device)
    else:  # "has_valid"
        idx_k_mask = torch.zeros(shape_idx_k, dtype=torch.bool, device=device)
        idx_k_mask[0, 0, 0] = True  # At least one valid entry

    if should_raise:
        with pytest.raises(ValueError, match="combined mask.*contains valid entries"):
            OuterGather.apply(z, idx_q, idx_k, axis, idx_q_mask, idx_k_mask)
    else:
        # Should not raise - the combined mask is all-False
        out = OuterGather.apply(z, idx_q, idx_k, axis, idx_q_mask, idx_k_mask)
        expected_shape = (2, 4, 5, 8, 3)
        assert out.shape == expected_shape
        assert (out == 0).all()


@pytest.mark.parametrize("leading_shape", [(), (2,), (2, 3)])
@pytest.mark.parametrize("n_dim", [1, 2, 3])
@pytest.mark.parametrize("make_empty", [False, True])
def test_compute_interval_overlap(leading_shape, n_dim, make_empty):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    # intervals_a: base starts at 0, ends at 5
    start_a = torch.zeros(leading_shape + (n_dim,), device=device, dtype=torch.long)
    end_a = start_a + 5
    intervals_a = torch.stack([start_a, end_a], dim=-1)  # (..., n_dim, 2)

    # intervals_b: optionally shift start beyond end_a in one dim to make empty
    start_b = torch.zeros_like(start_a)
    if make_empty:
        start_b = start_b + 6  # greater than end_a -> empty overlap
    end_b = start_b + 2
    intervals_b = torch.stack([start_b, end_b], dim=-1)

    overlap, mask = compute_interval_overlap(intervals_a, intervals_b)

    # Expected by direct max/min
    expected_start = torch.maximum(start_a, start_b)
    expected_end = torch.minimum(end_a, end_b)
    expected_overlap = torch.stack([expected_start, expected_end], dim=-1)
    expected_mask = torch.all(expected_end > expected_start, dim=-1)

    assert_tensors_identical(overlap, expected_overlap, check_grad=False, check_grad_fn=False)
    assert_tensors_identical(mask, expected_mask, check_grad=False, check_grad_fn=False)


@pytest.mark.parametrize("leading_shape", [(), (2,), (2, 3)])
@pytest.mark.parametrize("n_dim", [1, 2])
@pytest.mark.parametrize("overlap", [True, False])
def test_get_overlap_from_peers(leading_shape, n_dim, overlap):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    # rank_peers carries unique ids per position
    total = math.prod(leading_shape) if leading_shape else 1
    rank_peers = torch.arange(total, device=device).reshape(leading_shape if leading_shape else ())

    # intervals_a: base start = position value, length = 4 in all dims
    base = torch.arange(total, device=device).reshape(leading_shape if leading_shape else ())
    starts = torch.stack([base for _ in range(n_dim)], dim=-1)  # (..., n_dim)
    ends = starts + 4
    intervals_a = torch.stack([starts, ends], dim=-1)  # (..., n_dim, 2)

    # intervals_b: either overlaps or not
    if overlap:
        # Overlap window length 2 starting at base+1: overlap is [base+1, base+3)
        starts_b = torch.stack([base + 1 for _ in range(n_dim)], dim=-1)
        ends_b = starts_b + 2
        expected_mask = torch.ones_like(base, dtype=torch.bool)
        expected_intervals = torch.stack([starts_b, ends_b], dim=-1)
    else:
        # Disjoint: start beyond end_a
        starts_b = ends + 1
        ends_b = starts_b + 2
        expected_mask = torch.zeros_like(base, dtype=torch.bool)
        expected_intervals = torch.stack([starts_b, ends_b], dim=-1)
    intervals_b = torch.stack([starts_b, ends_b], dim=-1)  # (..., n_dim, 2)

    result = get_overlap_from_peers(rank_peers, intervals_a, intervals_b)

    expected_peers = rank_peers[expected_mask]
    expected_intervals = expected_intervals[expected_mask]

    assert len(result) == expected_peers.numel()
    for idx, item in enumerate(result):
        assert item["peer"] == expected_peers.view(-1)[idx].item()
        assert_tensors_identical(item["interval"], expected_intervals.view(-1, n_dim, 2)[idx], check_grad=False)


def parallel_assert_outer_gather(rank, grid_group_sizes, device_type, backend, env_map, dtype):
    """Run distributed outer gather on each rank."""

    monkeypatch = pytest.MonkeyPatch()
    if env_map is not None:
        for var_name, value in env_map.items():
            if value == "<INPUT_RANK>":
                monkeypatch.setenv(var_name, f"{rank}")
                continue
            monkeypatch.setenv(var_name, value)

    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    manager = DistributedManager()
    device = manager.device
    device_mesh = manager.device_mesh_subgroups
    device_mesh_flat = manager.device_mesh

    seed_by_rank(0, 42)

    # We test with different configurations of W, H, K, N
    # Following test_distributed_window_batch_attention style

    for input_shape_extra in [
        (2 * manager.group["dp"].size(), None, 3),
        (4 * manager.group["dp"].size(), None, 128),
    ]:  # 'None' will become the axis to be sharded in z
        if sum(1 for x in input_shape_extra if x is None) != 1:
            raise ValueError(
                f"There can be one and only one 'None' element in the input_shape but got {input_shape_extra}"
            )
        axis = input_shape_extra.index(None)

        for W, H, K_per_rank, N_per_rank in [(4, 8, 1, 8), (8, 16, 2, 16), (32, 128, 3, 100)]:
            N = N_per_rank * manager.group["cp_axis_0"].size()
            M = N_per_rank * manager.group["cp_axis_1"].size()

            # Construct shapes
            z_shape = input_shape_extra[:axis] + (N, M) + input_shape_extra[axis + 1 :]

            # Construct placements
            z_placements = [Shard(0), Shard(axis), Shard(axis + 1)]

            for idx_use_flat_device_mesh in [True, False]:
                if idx_use_flat_device_mesh:
                    K = K_per_rank * manager.group["cp"].size()
                    device_mesh_idx = device_mesh_flat
                    # equivalent to len(idx_n_shape) - 2
                    placements_idx = [Shard(0), Shard(len(input_shape_extra[:axis]))]
                else:
                    K = K_per_rank * manager.group["cp_axis_0"].size()
                    device_mesh_idx = device_mesh
                    # equivalent to len(idx_n_shape) - 2
                    placements_idx = [Shard(0), Shard(len(input_shape_extra[:axis])), Replicate()]

                idx_n_shape = input_shape_extra[:axis] + (K, W)
                idx_m_shape = input_shape_extra[:axis] + (K, H)

                # Test both without mask (None) and with random mask
                for use_mask in [False, True]:
                    label = (
                        f"z_shape:{z_shape}, idx_n_shape:{idx_n_shape}, idx_m_shape:{idx_m_shape}, "
                        f"z_placements:{z_placements}, placements_idx:{placements_idx}, axis:{axis}, "
                        f"idx_use_flat_device_mesh:{idx_use_flat_device_mesh}"
                    )
                    if use_mask:
                        label += " (masked)"
                        idx_n_mask_global = torch.rand(idx_n_shape, device=device) > 0.5
                        idx_m_mask_global = torch.rand(idx_m_shape, device=device) > 0.5

                        # Set second half of K dimension to all-False to simulate
                        # ranks with all-empty data when K is sharded across mesh
                        idx_n_mask_global[..., idx_n_shape[-2] // 2 :, :] = False
                        idx_m_mask_global[..., idx_m_shape[-2] // 2 :, :] = False
                    else:
                        idx_n_mask_global = None
                        idx_m_mask_global = None

                    # Create global data
                    z_global = torch.randn(z_shape, dtype=dtype, device=device, requires_grad=True)

                    idx_n_global = torch.randint(0, N, idx_n_shape, device=device, requires_grad=False)

                    idx_m_global = torch.randint(0, M, idx_m_shape, device=device, requires_grad=False)

                    # Reference
                    out_ref = OuterGather.apply(
                        z_global, idx_n_global, idx_m_global, axis, idx_n_mask_global, idx_m_mask_global
                    )
                    grad_out = torch.randn_like(out_ref)

                    # Use autograd to compute ref grad
                    # We clone z_global to avoid in-place modification issues if any (though here we just read)
                    out_ref.backward(grad_out)
                    grad_z_ref = z_global.grad

                    z_dtensor = distribute_tensor(z_global.detach().clone(), device_mesh, z_placements).requires_grad_(
                        True
                    )
                    idx_n_dtensor = distribute_tensor(idx_n_global, device_mesh_idx, placements_idx)
                    idx_m_dtensor = distribute_tensor(idx_m_global, device_mesh_idx, placements_idx)
                    idx_n_mask_dtensor = (
                        distribute_tensor(idx_n_mask_global, device_mesh_idx, placements_idx)
                        if idx_n_mask_global is not None
                        else None
                    )
                    idx_m_mask_dtensor = (
                        distribute_tensor(idx_m_mask_global, device_mesh_idx, placements_idx)
                        if idx_m_mask_global is not None
                        else None
                    )

                    # Forward
                    # The data in idx_n_global and idx_m_global are not actually contiguous due to randint
                    # but the code still works despite being inefficient for this test setting
                    out_dtensor = distributed_outer_gather(
                        z_dtensor,
                        idx_n_dtensor,
                        idx_m_dtensor,
                        axis=axis,
                        are_ids_contiguous=True,
                        idx_n_mask=idx_n_mask_dtensor,
                        idx_m_mask=idx_m_mask_dtensor,
                    )

                    # Verify Forward
                    out_local = out_dtensor.full_tensor().requires_grad_(True)
                    assert_tensors_identical(
                        out_local,
                        out_ref,
                        check_stride=False,
                        check_grad=False,
                        check_grad_fn=False,
                        msg=lambda m: f"{label} fwd output mismatch:\n {m}",
                    )

                    # Backward
                    grad_out_dtensor = distribute_tensor(
                        grad_out.detach().clone(), out_dtensor.device_mesh, out_dtensor.placements
                    )

                    out_dtensor.backward(grad_out_dtensor)

                    # Verify Backward
                    grad_z_local = z_dtensor.grad.full_tensor()
                    assert_tensors_identical(
                        grad_z_local,
                        grad_z_ref,
                        check_grad=False,
                        check_grad_fn=False,
                        rtol=1e-10,  # FP64 should be very precise but small accumulation errors can happen
                        atol=1e-10,
                        msg=lambda m: f"{label} bwd input gradient mismatch:\n {m}",
                    )

    DistributedManager.cleanup()


@pytest.mark.parametrize(
    "setup_env",
    [
        ((1, (2, 2)), True, "cuda", "ENV"),
        ((2, (2, 2)), True, "cuda", "ENV"),
    ],
    indirect=("setup_env",),
    ids=lambda x: f"dp:{x[0][0]}, cp:{x[0][1]}, specify_method:{x[1]}, device_type:{x[2]}, method_init:{x[3]}",
)
def test_distributed_outer_gather(setup_env):
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    spawn_multiprocessing(
        parallel_assert_outer_gather, world_size, grid_group_sizes, device_type, backend, env_per_rank, torch.float64
    )


def parallel_assert_distributed_outer_gather_w_one_hot(rank, grid_group_sizes, device_type, backend, env_map, dtype):
    """Run distributed outer gather with one-hot indices on each rank."""

    monkeypatch = pytest.MonkeyPatch()
    if env_map is not None:
        for var_name, value in env_map.items():
            if value == "<INPUT_RANK>":
                monkeypatch.setenv(var_name, f"{rank}")
                continue
            monkeypatch.setenv(var_name, value)

    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    manager = DistributedManager()
    device = manager.device
    device_mesh = manager.device_mesh_subgroups
    device_mesh_flat = manager.device_mesh

    seed_by_rank(0, 42)

    for input_shape_extra in [
        (2 * manager.group["dp"].size(), None, 3),
        (4 * manager.group["dp"].size(), None, 128),
    ]:  # 'None' will become the axis to be sharded in z
        if sum(1 for x in input_shape_extra if x is None) != 1:
            raise ValueError(
                f"There can be one and only one 'None' element in the input_shape but got {input_shape_extra}"
            )
        axis = input_shape_extra.index(None)

        for W, H, K_per_rank, N_per_rank in [(4, 8, 1, 8), (8, 16, 2, 16), (32, 128, 3, 100)]:
            N = N_per_rank * manager.group["cp_axis_0"].size()
            M = N_per_rank * manager.group["cp_axis_1"].size()

            # Construct shapes
            z_shape = input_shape_extra[:axis] + (N, M) + input_shape_extra[axis + 1 :]

            # Construct placements
            z_placements = [Shard(0), Shard(axis), Shard(axis + 1)]

            for idx_use_flat_device_mesh in [True, False]:
                if idx_use_flat_device_mesh:
                    K = K_per_rank * manager.group["cp"].size()
                    device_mesh_idx = device_mesh_flat
                    placements_idx = [Shard(0), Shard(len(input_shape_extra[:axis]))]
                else:
                    K = K_per_rank * manager.group["cp_axis_0"].size()
                    device_mesh_idx = device_mesh
                    placements_idx = [Shard(0), Shard(len(input_shape_extra[:axis])), Replicate()]

                idx_n_shape = input_shape_extra[:axis] + (K, W)
                idx_m_shape = input_shape_extra[:axis] + (K, H)

                # Test both without mask (None) and with random mask
                for use_mask in [False, True]:
                    label = (
                        f"z_shape:{z_shape}, idx_n_shape:{idx_n_shape}, idx_m_shape:{idx_m_shape}, "
                        f"z_placements:{z_placements}, placements_idx:{placements_idx}, axis:{axis}, "
                        f"idx_use_flat_device_mesh:{idx_use_flat_device_mesh}"
                    )
                    if use_mask:
                        label += " (masked)"
                        idx_n_mask_global = torch.rand(idx_n_shape, device=device) > 0.5
                        idx_m_mask_global = torch.rand(idx_m_shape, device=device) > 0.5

                        # Set second half of K dimension to all-False to simulate
                        # ranks with all-empty data when K is sharded across mesh
                        idx_n_mask_global[..., idx_n_shape[-2] // 2 :, :] = False
                        idx_m_mask_global[..., idx_m_shape[-2] // 2 :, :] = False
                    else:
                        idx_n_mask_global = None
                        idx_m_mask_global = None

                    # Create global data
                    z_global = torch.randn(z_shape, dtype=dtype, device=device, requires_grad=True)

                    # Create sorted index tensors then convert to one-hot
                    idx_n_indices = torch.randint(0, N, idx_n_shape, device=device)
                    idx_m_indices = torch.randint(0, M, idx_m_shape, device=device)

                    one_hot_n = F.one_hot(idx_n_indices, num_classes=N).to(dtype=z_global.dtype)
                    one_hot_m = F.one_hot(idx_m_indices, num_classes=M).to(dtype=z_global.dtype)

                    # Zero out one-hot at invalid positions when mask is provided
                    if idx_n_mask_global is not None:
                        one_hot_n = one_hot_n * idx_n_mask_global.unsqueeze(-1).to(dtype=z_global.dtype)
                    if idx_m_mask_global is not None:
                        one_hot_m = one_hot_m * idx_m_mask_global.unsqueeze(-1).to(dtype=z_global.dtype)

                    # Reference using einsum with one-hot
                    feature_shape = z_shape[axis + 2 :]
                    feature_flat = int(math.prod(feature_shape)) if feature_shape else 1
                    z_flat = z_global.view(*z_shape[: axis + 2], feature_flat)
                    ref_flat = torch.einsum("...nmf,...kwn,...khm->...kwhf", z_flat, one_hot_n, one_hot_m)
                    out_ref = ref_flat.view(*z_shape[:axis], K, W, H, *feature_shape)

                    grad_out = torch.randn_like(out_ref)

                    # Use autograd to compute ref grad
                    out_ref.backward(grad_out)
                    grad_z_ref = z_global.grad

                    z_dtensor = distribute_tensor(z_global.detach().clone(), device_mesh, z_placements).requires_grad_(
                        True
                    )
                    one_hot_n_dtensor = distribute_tensor(one_hot_n, device_mesh_idx, placements_idx)
                    one_hot_m_dtensor = distribute_tensor(one_hot_m, device_mesh_idx, placements_idx)

                    # Convert one-hot DTensors to index DTensors via shardwise_argmax
                    idx_n_dtensor = shardwise_argmax(one_hot_n_dtensor, dim=-1, keepdim=False)
                    idx_m_dtensor = shardwise_argmax(one_hot_m_dtensor, dim=-1, keepdim=False)

                    idx_n_mask_dtensor = (
                        distribute_tensor(idx_n_mask_global, device_mesh_idx, placements_idx)
                        if idx_n_mask_global is not None
                        else None
                    )
                    idx_m_mask_dtensor = (
                        distribute_tensor(idx_m_mask_global, device_mesh_idx, placements_idx)
                        if idx_m_mask_global is not None
                        else None
                    )

                    # Forward
                    # The data in idx_n_global and idx_m_global are not actually contiguous due to randint
                    # but the code still works despite being inefficient for this test setting
                    out_dtensor = distributed_outer_gather(
                        z_dtensor,
                        idx_n_dtensor,
                        idx_m_dtensor,
                        axis=axis,
                        are_ids_contiguous=True,
                        idx_n_mask=idx_n_mask_dtensor,
                        idx_m_mask=idx_m_mask_dtensor,
                    )

                    # Verify Forward
                    out_local = out_dtensor.full_tensor().requires_grad_(True)
                    assert_tensors_identical(
                        out_local,
                        out_ref,
                        check_stride=False,
                        check_grad=False,
                        check_grad_fn=False,
                        msg=lambda m: f"{label} fwd output mismatch:\n {m}",
                    )

                    # Backward
                    grad_out_dtensor = distribute_tensor(
                        grad_out.detach().clone(), out_dtensor.device_mesh, out_dtensor.placements
                    )

                    out_dtensor.backward(grad_out_dtensor)

                    # Verify Backward
                    grad_z_local = z_dtensor.grad.full_tensor()
                    assert_tensors_identical(
                        grad_z_local,
                        grad_z_ref,
                        check_grad=False,
                        check_grad_fn=False,
                        rtol=1e-10,
                        atol=1e-10,
                        msg=lambda m: f"{label} bwd input gradient mismatch:\n {m}",
                    )

    DistributedManager.cleanup()


@pytest.mark.parametrize(
    "setup_env",
    [
        ((1, (2, 2)), True, "cuda", "ENV"),
        ((2, (2, 2)), True, "cuda", "ENV"),
    ],
    indirect=("setup_env",),
    ids=lambda x: f"dp:{x[0][0]}, cp:{x[0][1]}, specify_method:{x[1]}, device_type:{x[2]}, method_init:{x[3]}",
)
def test_outer_distributed_gather_w_one_hot(setup_env):
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    spawn_multiprocessing(
        parallel_assert_distributed_outer_gather_w_one_hot,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        torch.float64,
    )
