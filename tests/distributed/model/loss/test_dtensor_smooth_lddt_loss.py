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

from boltz.distributed.comm import TransposeComm
from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.loss.diffusion import smooth_lddt_loss
from boltz.model.loss.diffusion import smooth_lddt_loss as serial_smooth_lddt_loss_v1
from boltz.model.loss.diffusionv2 import smooth_lddt_loss as serial_smooth_lddt_loss_v2
from boltz.testing.utils import spawn_multiprocessing


def parallel_assert_smooth_lddt_loss(
    rank: int,
    payload: tuple,
):
    (
        multiplicity,
        nucleic_acid_cutoff,
        other_cutoff,
        pred_coords,
        true_coords,
        is_nucleotide,
        coords_mask,
        expected_loss_host,
        expected_pred_coords_grad_host,
        grid_group_sizes,
        device_type,
        backend,
        env_map,
        v2,
    ) = payload

    monkeypatch = pytest.MonkeyPatch()
    if env_map is not None:
        for var_name, value in env_map.items():
            if value == "<INPUT_RANK>":
                monkeypatch.setenv(var_name, f"{rank}")
                continue
            monkeypatch.setenv(var_name, value)

    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    manager = DistributedManager()

    layout_map = manager.layout_subgroups["cp"]
    device_mesh = manager.device_mesh_subgroups
    comm = TransposeComm(manager.group["cp"], layout_map)

    pred_coords = distribute_tensor(
        pred_coords.to(manager.device),
        device_mesh=device_mesh,
        placements=(Shard(0), Shard(1), Replicate()),
    )
    true_coords = distribute_tensor(
        true_coords.to(manager.device),
        device_mesh=device_mesh,
        placements=(Shard(0), Shard(1), Replicate()),
    )
    is_nucleotide = distribute_tensor(
        is_nucleotide.to(manager.device),
        device_mesh=device_mesh,
        placements=(Shard(0), Shard(1), Replicate()),
    )
    coords_mask = distribute_tensor(
        coords_mask.to(manager.device),
        device_mesh=device_mesh,
        placements=(Shard(0), Shard(1), Replicate()),
    )

    pred_coords.requires_grad_(True)

    loss = smooth_lddt_loss(
        pred_coords,
        true_coords,
        is_nucleotide,
        coords_mask,
        multiplicity=multiplicity,
        comm=comm,
        nucleic_acid_cutoff=nucleic_acid_cutoff,
        other_cutoff=other_cutoff,
        v2=v2,
    )

    assert (
        loss.shape == expected_loss_host.shape
    ), f"Loss shape mismatch: expected {expected_loss_host.shape}, got {loss.shape}"
    assert (
        loss.stride() == expected_loss_host.stride()
    ), f"Loss stride mismatch: expected {expected_loss_host.stride()}, got {loss.stride()}"

    loss.backward()

    torch.testing.assert_close(loss.full_tensor().cpu(), expected_loss_host)

    assert (
        pred_coords.grad.shape == expected_pred_coords_grad_host.shape
    ), f"Pred coords grad shape mismatch: expected {expected_pred_coords_grad_host.shape}, got {pred_coords.grad.shape}"
    assert (
        pred_coords.grad.stride() == expected_pred_coords_grad_host.stride()
    ), f"Pred coords grad stride mismatch: expected {expected_pred_coords_grad_host.stride()}, got {pred_coords.grad.stride()}"
    torch.testing.assert_close(
        pred_coords.grad.full_tensor().cpu(),
        expected_pred_coords_grad_host,
    )

    # Clean up
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
    ids=lambda x: f"dp={x[0][0]}, cp={x[0][1]}, specify_method={x[1]}, device_type={x[2]}, method_init={x[3]}",
)
@pytest.mark.parametrize("multiplicity", [1, 2])
@pytest.mark.parametrize("v2", [True, False], ids=["v2", "v1"])
def test_smooth_lddt_loss(
    setup_env: tuple,
    multiplicity: int,
    v2: bool,
    nucleic_acid_cutoff: float = 30.0,
    other_cutoff: float = 15.0,
    dtype: torch.dtype = torch.float64,
):
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    seed = 42
    rng = torch.Generator(device_type)
    rng.manual_seed(seed)

    B = 2 * grid_group_sizes["dp"]
    N_per_shard = 32
    N = N_per_shard * grid_group_sizes["cp"][0]

    # Make coordinates large enough that some distances exceed cutoff (default of 30)
    pred_coords = (
        torch.randn(B * multiplicity, N, 3, generator=rng, device=device_type, dtype=dtype) * nucleic_acid_cutoff
    )
    true_coords = (
        torch.randn(B * multiplicity, N, 3, generator=rng, device=device_type, dtype=dtype) * nucleic_acid_cutoff
    )
    pred_coords.requires_grad_(True)

    # Multiplicity is called within smooth_lddt_loss for features
    is_nucleotide = torch.randint(0, 2, (B, N), generator=rng, device=device_type, dtype=dtype)
    coords_mask = torch.randint(0, 2, (B, N), generator=rng, device=device_type, dtype=dtype)

    # mask the last 2 atoms in each shard to emulate inserted virtual atoms in the middle of the atom sequence
    for cp_rank in range(grid_group_sizes["cp"][0]):
        end_idx = (cp_rank + 1) * N_per_shard
        start_idx = end_idx - 2
        coords_mask[:, start_idx:end_idx] = 0

    serial_smooth_lddt_loss = serial_smooth_lddt_loss_v2 if v2 else serial_smooth_lddt_loss_v1
    reference_loss = serial_smooth_lddt_loss(
        pred_coords,
        true_coords,
        is_nucleotide,
        coords_mask,
        multiplicity=multiplicity,
        nucleic_acid_cutoff=nucleic_acid_cutoff,
        other_cutoff=other_cutoff,
    )
    reference_loss.backward()

    # Call subprocesses for parallel testing
    payload = (
        multiplicity,
        nucleic_acid_cutoff,
        other_cutoff,
        pred_coords.detach().clone().cpu(),
        true_coords.clone().cpu(),
        is_nucleotide.clone().cpu(),
        coords_mask.clone().cpu(),
        reference_loss.detach().clone().cpu(),
        pred_coords.grad.detach().clone().cpu(),
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        v2,
    )
    spawn_multiprocessing(
        parallel_assert_smooth_lddt_loss,
        world_size,
        payload,
    )
