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

import pytest
import torch
from torch.distributed.tensor import Replicate, Shard, distribute_tensor

from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.layers.scatter import distributed_scatter_reduce
from boltz.testing.utils import assert_tensors_identical, seed_by_rank, spawn_multiprocessing


def einsum_scatter_reduce(idx_global, src_global, N_output, axis, reduce, idx_mask_global=None):
    """Reference scatter-reduce using one-hot encoding and einsum.

    This computes the same result as:
        output = zeros(output_shape)
        output.scatter_reduce_(axis, idx_expanded, src, reduce=reduce)

    But uses one-hot encoding and einsum which is more explicit and easier to verify.

    Args:
        idx_global: Index tensor (*batch, N_src) with values in [0, N_output)
        src_global: Source tensor (*batch, N_src, *features)
        N_output: Size of the output's scatter axis
        axis: The scatter axis position
        reduce: "sum" or "mean"
        idx_mask_global: Optional mask tensor (*batch, N_src), True=valid, False=invalid
    """
    # Get shapes
    batch_shape = idx_global.shape[:axis]
    N_src = idx_global.shape[axis]
    feature_shape = src_global.shape[axis + 1 :]
    dtype = src_global.dtype

    # Flatten batch dimensions for easier processing
    B = torch.Size(batch_shape).numel() if batch_shape else 1
    F = torch.Size(feature_shape).numel() if feature_shape else 1

    # Reshape to (B, N_src) and (B, N_src, F)
    idx_flat = idx_global.reshape(B, N_src)
    src_flat = src_global.reshape(B, N_src, F)

    # Create one-hot tensor: (B, N_src, N_output)
    # one_hot[b, i, j] = 1 if idx_flat[b, i] == j
    one_hot = torch.nn.functional.one_hot(idx_flat, num_classes=N_output).to(dtype)

    # Apply mask to one-hot (zeroing out masked entries means they don't contribute to output)
    if idx_mask_global is not None:
        mask_flat = idx_mask_global.reshape(B, N_src, 1).to(dtype)
        one_hot = one_hot * mask_flat

    # Compute output using einsum
    # out[b, j, f] = sum_i(one_hot[b, i, j] * src[b, i, f])
    # einsum: "bin,bif->bnf" where i=N_src, n=N_output, f=F
    out_flat = torch.einsum("bin,bif->bnf", one_hot, src_flat)

    if reduce == "mean":
        # count[b, j] = sum_i(one_hot[b, i, j]) = number of elements scattered to position j
        count = one_hot.sum(dim=1, keepdim=True)  # (B, 1, N_output)
        count = count.permute(0, 2, 1)  # (B, N_output, 1)
        # Clamp to avoid division by zero; positions with count=0 already have out=0 from einsum
        out_flat = out_flat / count.clamp(min=1)

    # Reshape back to original structure
    output_shape = batch_shape + (N_output,) + feature_shape
    out = out_flat.reshape(output_shape)

    return out


def parallel_assert_scatter_reduce(rank, grid_group_sizes, device_type, backend, env_map, dtype, reduce_op):
    """Test distributed scatter_reduce against reference implementation."""
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

    seed_by_rank(0, 42)

    # Test configurations with different shape patterns
    # Format: (batch_shape_before_axis, N_per_rank, N_src_per_rank, feature_shape_after_N)
    # output: (*batch, N, *features)
    # idx: (*batch, N_src)
    # src: (*batch, N_src, *features)
    shape_configs = [
        # Simple cases: output (N, F), idx (N_src,), src (N_src, F)
        ((), 8, 12, (3,)),  # output: (N,3), idx: (N_src,), src: (N_src,3)
        ((), 16, 8, ()),  # output: (N,), idx: (N_src,), src: (N_src,)
        # With batch: output (B, N, F), idx (B, N_src), src (B, N_src, F)
        ((4,), 8, 12, ()),  # output: (4,N), idx: (4,N_src), src: (4,N_src)
        ((4,), 8, 12, (3,)),  # output: (4,N,3), idx: (4,N_src), src: (4,N_src,3)
        ((2,), 16, 8, (5,)),  # output: (2,N,5), idx: (2,N_src), src: (2,N_src,5)
        # Multi-dim features
        ((2,), 8, 10, (3, 4)),  # output: (2,N,3,4), idx: (2,N_src), src: (2,N_src,3,4)
        # Larger test case
        ((2,), 100, 1000, (3,)),  # output: (2,N,3), idx: (2,N_src), src: (2,N_src,3)
    ]

    for batch_shape, N_per_rank, N_src_per_rank, feature_shape in shape_configs:
        for device_mesh in [manager.device_mesh_subgroups, manager.device_mesh]:
            mesh_ndim = device_mesh.ndim
            axis = len(batch_shape)  # axis is right after batch dims

            # Determine placement strategy
            size_group_shard_axis = None
            if axis >= 1 and mesh_ndim >= 2:
                # shard leading tensor axes as well as 'axis'
                placements = (Shard(0), Shard(axis)) + (Replicate(),) * (mesh_ndim - 2)
                size_group_shard_axis = device_mesh.size(1)
            elif mesh_ndim >= 2:
                # axis == 0: only shard 'axis'
                placements = (Replicate(),) * (mesh_ndim - 1) + (Shard(axis),)
                size_group_shard_axis = device_mesh.size(-1)
            else:
                # 1D mesh: shard on axis
                placements = (Shard(axis),) + (Replicate(),) * (mesh_ndim - 1)
                size_group_shard_axis = device_mesh.size(0)

            if size_group_shard_axis is None:
                raise ValueError(f"size_group_shard_axis is None for axis {axis} and device_mesh {device_mesh}")

            N = N_per_rank * size_group_shard_axis
            N_src = N_src_per_rank * size_group_shard_axis

            # Build shapes
            output_shape = batch_shape + (N,) + feature_shape
            idx_shape = batch_shape + (N_src,)
            src_shape = batch_shape + (N_src,) + feature_shape

            # Test both without mask (None) and with random mask
            for use_mask in [False, True]:
                label = f"output:{output_shape}, idx:{idx_shape}, src:{src_shape}, axis:{axis}, reduce:{reduce_op}"
                if use_mask:
                    label += " (masked)"
                    idx_mask_global = torch.rand(idx_shape, device=device) > 0.5
                else:
                    idx_mask_global = None

                # Create test data
                idx_global = torch.randint(0, N, idx_shape, device=device)
                src_global = torch.randn(src_shape, dtype=dtype, device=device, requires_grad=True)

                # Reference computation using einsum with one-hot encoding
                out_ref = einsum_scatter_reduce(idx_global, src_global.detach(), N, axis, reduce_op, idx_mask_global)

                # Create grad_out for backward test
                grad_out = torch.randn_like(out_ref)

                # Reference backward using autograd
                src_global_for_bwd = src_global.detach().clone().requires_grad_(True)
                out_ref_bwd = einsum_scatter_reduce(idx_global, src_global_for_bwd, N, axis, reduce_op, idx_mask_global)
                out_ref_bwd.backward(grad_out)
                grad_src_ref = src_global_for_bwd.grad

                # Distributed computation
                idx_dtensor = distribute_tensor(idx_global, device_mesh, placements)
                src_dtensor = distribute_tensor(src_global.detach().clone(), device_mesh, placements).requires_grad_(
                    True
                )
                idx_mask_dtensor = (
                    distribute_tensor(idx_mask_global, device_mesh, placements) if idx_mask_global is not None else None
                )

                out_dtensor = distributed_scatter_reduce(
                    output_size_per_rank=N_per_rank,
                    axis=axis,
                    idx=idx_dtensor,
                    src=src_dtensor,
                    reduce=reduce_op,
                    idx_mask=idx_mask_dtensor,
                    are_ids_contiguous=True,
                )

                # Compare forward output
                out_local = out_dtensor.full_tensor()
                assert_tensors_identical(
                    out_local.detach(),
                    out_ref,
                    check_stride=False,
                    check_grad=False,
                    check_grad_fn=False,
                    rtol=1e-10,
                    atol=1e-10,
                    msg=lambda m: f"{label} fwd output mismatch:\n {m}",
                )

                # Backward pass
                grad_out_dtensor = distribute_tensor(
                    grad_out.detach().clone(), out_dtensor.device_mesh, out_dtensor.placements
                )

                out_dtensor.backward(grad_out_dtensor)

                # Compare gradients
                grad_src_local = src_dtensor.grad.full_tensor()

                assert_tensors_identical(
                    grad_src_local,
                    grad_src_ref,
                    check_grad=False,
                    check_grad_fn=False,
                    rtol=1e-10,
                    atol=1e-10,
                    msg=lambda m: f"{label} bwd src gradient mismatch:\n {m}",
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
@pytest.mark.parametrize("reduce_op", ["sum", "mean"])
def test_distributed_scatter_reduce(setup_env, reduce_op):
    """Test distributed scatter_reduce operation."""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() < {world_size}")

    spawn_multiprocessing(
        parallel_assert_scatter_reduce,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        torch.float64,
        reduce_op,
    )
