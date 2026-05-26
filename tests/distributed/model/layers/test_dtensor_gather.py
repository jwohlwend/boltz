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
import torch.nn.functional as F
from torch.distributed.tensor import Replicate, Shard, distribute_tensor

from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.layers.gather import distributed_gather
from boltz.testing.utils import assert_tensors_identical, seed_by_rank, spawn_multiprocessing


def parallel_assert_gather(rank, grid_group_sizes, device_type, backend, env_map, dtype):
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

    shape_extras = [
        (None,),
        (4, None),
        (None, 3),
        (2, None, 3),
    ]

    for shape_extra in shape_extras:
        if sum(1 for x in shape_extra if x is None) != 1:
            raise ValueError(f"There can be one and only one 'None' element in shape_extra but got {shape_extra}")
        axis = shape_extra.index(None)

        for N_per_rank, K_per_rank, W in [(8, 6, 2), (16, 4, 3)]:
            for device_mesh in [manager.device_mesh_subgroups, manager.device_mesh]:
                # placements: shard along axis (on mesh dim 0)
                mesh_ndim = device_mesh.ndim
                size_group_shard_axis = None
                if axis >= 1 and mesh_ndim >= 2:
                    # shard leading tensor axes as well as 'axis
                    placements = (Shard(0), Shard(axis)) + (Replicate(),) * (mesh_ndim - 2)
                    size_group_shard_axis = device_mesh.size(1)
                elif mesh_ndim >= 2:
                    # axis == 0 only shard 'axis'
                    placements = (Replicate(),) * (mesh_ndim - 1) + (Shard(axis),)
                    size_group_shard_axis = device_mesh.size(-1)
                else:
                    # not enough mesh dim to shard other than 'axis'
                    placements = (Shard(axis),) + (Replicate(),) * (mesh_ndim - 1)
                    size_group_shard_axis = device_mesh.size(0)

                if size_group_shard_axis is None:
                    raise ValueError(f"size_group_shard_axis is None for axis {axis} and device_mesh {device_mesh}")

                N = N_per_rank * size_group_shard_axis
                K = K_per_rank * size_group_shard_axis

                x_shape = shape_extra[:axis] + (N,) + shape_extra[axis + 1 :]
                idx_shape = shape_extra[:axis] + (K, W)

                # Test both without mask (None) and with random mask
                for use_mask in [False, True]:
                    label = f"x_shape:{x_shape}, idx_shape:{idx_shape}, axis:{axis}"
                    if use_mask:
                        label += " (masked)"
                        idx_mask_global = torch.rand(idx_shape, device=device) > 0.5
                    else:
                        idx_mask_global = None

                    x_global = torch.randn(x_shape, dtype=dtype, device=device, requires_grad=True)
                    idx_global = torch.randint(0, N, idx_shape, device=device)

                    # Reference using one-hot + einsum
                    idx_onehot = F.one_hot(idx_global, num_classes=N).to(dtype=dtype)
                    # Zero out one-hot at invalid positions when mask is provided
                    if idx_mask_global is not None:
                        idx_onehot = idx_onehot * idx_mask_global.unsqueeze(-1).to(dtype=dtype)

                    x_flat = x_global.reshape(
                        *x_global.shape[:axis], x_global.shape[axis], x_global.shape[axis + 1 :].numel()
                    )
                    out_ref_flat = torch.einsum("...nd,...kwn->...kwd", x_flat, idx_onehot)
                    out_ref = out_ref_flat.reshape(
                        *x_global.shape[:axis], idx_shape[-2], idx_shape[-1], *x_global.shape[axis + 1 :]
                    )
                    grad_out = torch.randn_like(out_ref)

                    out_ref.backward(grad_out)
                    grad_x_ref = x_global.grad

                    x_dtensor = distribute_tensor(x_global.detach().clone(), device_mesh, placements).requires_grad_(
                        True
                    )
                    idx_dtensor = distribute_tensor(idx_global, device_mesh, placements)
                    idx_mask_dtensor = (
                        distribute_tensor(idx_mask_global, device_mesh, placements)
                        if idx_mask_global is not None
                        else None
                    )

                    out_dtensor = distributed_gather(
                        x_dtensor, idx_dtensor, axis=axis, are_ids_contiguous=True, idx_mask=idx_mask_dtensor
                    )

                    out_local = out_dtensor.full_tensor().requires_grad_(True)
                    assert_tensors_identical(
                        out_local,
                        out_ref,
                        check_stride=False,
                        check_grad=False,
                        check_grad_fn=False,
                        msg=lambda m: f"{label} fwd output mismatch:\n {m}",
                    )

                    grad_out_dtensor = distribute_tensor(
                        grad_out.detach().clone(), out_dtensor.device_mesh, out_dtensor.placements
                    )

                    out_dtensor.backward(grad_out_dtensor)

                    grad_x_local = x_dtensor.grad.full_tensor()
                    assert_tensors_identical(
                        grad_x_local,
                        grad_x_ref,
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
def test_distributed_gather(setup_env):
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    spawn_multiprocessing(
        parallel_assert_gather,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        torch.float64,
    )
