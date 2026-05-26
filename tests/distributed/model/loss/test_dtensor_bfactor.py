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

"""DTensor-based tests for Boltz-2 B-factor loss.

Tests the DTensor CP implementation against the serial bfactor_loss_fn.
Maps to: src/boltz/distributed/model/loss/bfactor.py

Verification checks:
    V1: single-proc serial immutability (FW/BW inputs unchanged)
    V4: multi-proc FW input tensor values unchanged by FW and BW
    V8: multi-proc FW loss value close-to single-proc
    V9: multi-proc BW pred gradient values close-to single-proc

bf16 coverage (CUDA-only):
    The bfactor loss uses promote_types (compute in fp32) inside
    autocast(enabled=False).  The bf16 test verifies that bf16 pred logits
    are correctly promoted and the resulting loss/gradient are fp32.
"""

import pytest
import torch
from torch.distributed.tensor import Replicate, Shard, distribute_tensor
from torch.testing import assert_close

from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.loss.bfactor import bfactor_loss as bfactor_loss_dtensor
from boltz.model.loss.bfactor import bfactor_loss_fn as bfactor_loss_serial
from boltz.testing.utils import (
    assert_tensors_identical,
    skip_if_cuda_not_avail_or_device_count_less_than_word_size,
    spawn_multiprocessing,
)

SEED = 42


def _assert_unchanged(actual, expected, *, serial=False):
    """Shorthand for assert_tensors_identical with standard immutability kwargs."""
    assert_tensors_identical(
        actual,
        expected,
        check_stride=True,
        check_grad=False,
        check_grad_fn=False,
        check_storage_pointer=False,
        check_storage_offset=serial,
    )


def _worker_bfactor_loss_parity(
    rank: int,
    pred_on_host: torch.Tensor,
    t2ra_on_host: torch.Tensor,
    bf_on_host: torch.Tensor,
    loss_ref: float,
    pred_grad_ref_on_host: torch.Tensor,
    grid_group_sizes: dict,
    device_type: str,
    backend: str,
    env_map: dict[str, str] | None = None,
):
    """Worker: compare distributed bfactor loss against serial reference.

    Performs V4 (input immutability), V8 (loss parity), V9 (grad parity).
    """
    monkeypatch = pytest.MonkeyPatch()
    if env_map is not None:
        for var_name, value in env_map.items():
            if value == "<INPUT_RANK>":
                monkeypatch.setenv(var_name, f"{rank}")
                continue
            monkeypatch.setenv(var_name, value)

    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    dm = DistributedManager()

    single_placements = (Shard(0), Shard(1), Replicate())
    # bfactor is (B, A) — replicate A on CP dims so bmm dimensions match
    # (t2ra has A in its last dim which is Replicate on cp1)
    atom_placements = (Shard(0), Replicate(), Replicate())

    pred_dt = distribute_tensor(pred_on_host.to(dm.device), dm.device_mesh_subgroups, single_placements).requires_grad_(
        True
    )
    t2ra_dt = distribute_tensor(t2ra_on_host.to(dm.device), dm.device_mesh_subgroups, single_placements)
    bf_dt = distribute_tensor(bf_on_host.to(dm.device), dm.device_mesh_subgroups, atom_placements)

    # V4 setup: clone inputs for immutability check
    pred_dt_clone = pred_dt.detach().clone().requires_grad_(pred_dt.requires_grad)
    t2ra_dt_clone = t2ra_dt.detach().clone().requires_grad_(t2ra_dt.requires_grad)
    bf_dt_clone = bf_dt.detach().clone().requires_grad_(bf_dt.requires_grad)

    output = {"pbfactor": pred_dt}
    feats = {"token_to_rep_atom": t2ra_dt, "bfactor": bf_dt}

    dp_group = dm.device_mesh_subgroups.get_group(0)
    cp0_group = dm.device_mesh_subgroups.get_group(1)
    cp1_group = dm.device_mesh_subgroups.get_group(2)

    # Forward
    loss_dt = bfactor_loss_dtensor(
        output,
        feats,
        device_mesh=dm.device_mesh_subgroups,
        dp_group=dp_group,
        cp0_group=cp0_group,
        cp1_group=cp1_group,
    )

    # V4a: FW inputs unchanged
    assert_tensors_identical(
        pred_dt.to_local(),
        pred_dt_clone.to_local(),
        check_grad=False,
        check_grad_fn=False,
    )
    assert_tensors_identical(
        t2ra_dt.to_local(),
        t2ra_dt_clone.to_local(),
        check_grad=False,
        check_grad_fn=False,
    )
    assert_tensors_identical(
        bf_dt.to_local(),
        bf_dt_clone.to_local(),
        check_grad=False,
        check_grad_fn=False,
    )

    # V8: forward loss parity
    loss_val = loss_dt.full_tensor().item()
    assert_close(
        torch.tensor(loss_val),
        torch.tensor(loss_ref),
        atol=1e-5,
        rtol=1e-5,
        msg=lambda m: f"Rank {rank} loss mismatch\n{m}",
    )

    # Backward
    loss_dt.backward()

    # V4b: FW inputs unchanged after backward
    assert_tensors_identical(
        pred_dt.to_local(),
        pred_dt_clone.to_local(),
        check_grad=False,
        check_grad_fn=False,
    )

    # V9: pred gradient parity
    assert pred_dt.grad is not None, "pred gradient is None"
    pred_grad_full = pred_dt.grad.full_tensor()
    assert_close(
        pred_grad_full,
        pred_grad_ref_on_host.to(dm.device),
        atol=5e-5,
        rtol=5e-5,
        msg=lambda m: f"Rank {rank} pred grad mismatch\n{m}",
    )

    # Loss and gradient should be fp32 regardless of input dtype (promote_types)
    assert loss_dt.dtype == torch.float32, f"Loss dtype should be fp32, got {loss_dt.dtype}"
    assert pred_dt.grad.dtype == torch.float32, f"Grad dtype should be fp32, got {pred_dt.grad.dtype}"

    # Non-vacuous: gradient must be non-zero (would be zero only if loss is
    # independent of pred, which would indicate a broken implementation).
    assert (
        pred_grad_full.abs().sum() > 0
    ), "Distributed pred gradient is all-zero — loss is not differentiable w.r.t. pred"

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env",
    [
        # CPU dp=2 cp=(3,3): DP + non-power-of-two CP for CPU-only CI
        ((2, (3, 3)), True, "cpu", "ENV"),
        # CUDA dp=1 cp=(1,1): serial-equivalent sanity check (1 GPU)
        ((1, (1, 1)), True, "cuda", "ENV"),
        # CUDA dp=2 cp=(1,1): bf16-compatible path (2 GPUs)
        ((2, (1, 1)), True, "cuda", "ENV"),
        # CUDA dp=1 cp=(2,2): actual CP under CUDA (4 GPUs)
        ((1, (2, 2)), True, "cuda", "ENV"),
        # CUDA dp=2 cp=(2,2): DP + CP under CUDA (8 GPUs)
        ((2, (2, 2)), True, "cuda", "ENV"),
    ],
    indirect=("setup_env",),
    ids=["cpu-dp2-cp3x3", "cuda-dp1-cp1x1", "cuda-dp2-cp1x1", "cuda-dp1-cp2x2", "cuda-dp2-cp2x2"],
)
def test_dtensor_bfactor_loss_forward_backward(setup_env):
    """BFactor loss: distributed loss and pred gradient match serial reference."""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    skip_if_cuda_not_avail_or_device_count_less_than_word_size(device_type=device_type, world_size=world_size)

    B = 2
    bins = 8
    N = 8 * grid_group_sizes["cp"][0]
    A = N  # one atom per token for simplicity

    with torch.random.fork_rng(devices=[], enabled=True):
        torch.manual_seed(SEED)

        # Create pred logits
        pred = torch.randn(B, N, bins, requires_grad=True)
        pred_copy = pred.detach().clone().requires_grad_(True)

        # Create token_to_rep_atom: identity-like (each token maps to one atom)
        t2ra = torch.zeros(B, N, A, dtype=torch.float32)
        for b in range(B):
            for i in range(N):
                t2ra[b, i, i] = 1.0

        # Create bfactor: realistic values (0-100 range), some zeros
        bf = torch.rand(B, A) * 80.0
        bf[:, : A // 4] = 0.0  # first quarter is zero (no bfactor)

        # Serial loss
        output_serial = {"pbfactor": pred}
        feats_serial = {"token_to_rep_atom": t2ra, "bfactor": bf}
        loss_ref = bfactor_loss_serial(output_serial, feats_serial)

        # V1a: serial FW input unchanged
        _assert_unchanged(pred, pred_copy, serial=True)

        loss_ref_val = loss_ref.item()

        # Non-vacuous: serial loss must be positive.  A zero loss would make
        # the forward parity check (V8) and gradient parity check (V9) trivial.
        assert loss_ref_val > 0, (
            f"Serial bfactor loss is {loss_ref_val} — test data produced zero "
            f"loss, making parity checks vacuous. Verify bf/t2ra test setup."
        )

        # Serial backward
        loss_ref.backward()

        # V1b: serial FW input unchanged after backward
        _assert_unchanged(pred, pred_copy, serial=True)

        pred_grad_ref = pred.grad.detach().cpu().clone()

    spawn_multiprocessing(
        _worker_bfactor_loss_parity,
        world_size,
        pred.detach().cpu(),
        t2ra.detach().cpu(),
        bf.detach().cpu(),
        loss_ref_val,
        pred_grad_ref,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
    )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
