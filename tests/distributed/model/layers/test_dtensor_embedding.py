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

"""DTensor parity tests for EmbeddingParamsReplicated.

Tests the distributed embedding wrapper against serial nn.Embedding.

Verification checks:
    V4a: multi-proc FW input tensor values unchanged by FW
    V4b: multi-proc FW input tensor values unchanged after BW
    V5:  multi-proc BW input tensor values unchanged by BW
    V8:  multi-proc FW output tensor values close-to single-proc
    V10: multi-proc weight gradient values close-to single-proc
    V10b: replicated weight gradients identical across all CP ranks
"""

import pytest
import torch
import torch.nn as nn
from torch.distributed.tensor import DTensor, Shard, distribute_tensor
from torch.testing import assert_close

from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.layers.embedding import EmbeddingParamsReplicated
from boltz.testing.utils import (
    assert_all_identical,
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


def _worker_embedding_parity(
    rank: int,
    input_on_host: torch.Tensor,
    output_ref_on_host: torch.Tensor,
    grad_output_on_host: torch.Tensor,
    weight_grad_ref_on_host: torch.Tensor,
    state_dict: dict,
    num_embeddings: int,
    embedding_dim: int,
    padding_idx: int | None,
    grid_group_sizes: dict,
    device_type: str,
    backend: str,
    env_map: dict[str, str] | None = None,
):
    """Worker: compare distributed EmbeddingParamsReplicated against serial nn.Embedding."""
    monkeypatch = pytest.MonkeyPatch()
    if env_map is not None:
        for var_name, value in env_map.items():
            monkeypatch.setenv(var_name, f"{rank}" if value == "<INPUT_RANK>" else value)

    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    dm = DistributedManager()

    try:
        serial_emb = nn.Embedding(num_embeddings, embedding_dim, padding_idx=padding_idx)
        serial_emb.load_state_dict(state_dict)
        serial_emb = serial_emb.to(dm.device)

        dist_emb = EmbeddingParamsReplicated(serial_emb, dm.device_mesh_subgroups)

        # Pair-like placements for 3D input [B, N, N]
        pair_placements = (Shard(0), Shard(1), Shard(2))
        x_dt = distribute_tensor(input_on_host.to(dm.device), dm.device_mesh_subgroups, pair_placements)

        # V4a setup
        x_dt_clone = x_dt.detach().clone()

        out = dist_emb(x_dt)

        # V4a: FW input unchanged
        _assert_unchanged(x_dt.to_local(), x_dt_clone.to_local())

        # V8: forward parity
        assert_close(
            out.full_tensor(),
            output_ref_on_host.to(dm.device),
            atol=0,
            rtol=0,
            msg=lambda m: f"Rank {rank} forward output mismatch\n{m}",
        )

        # Backward
        grad_out_dt = distribute_tensor(grad_output_on_host.to(dm.device), dm.device_mesh_subgroups, pair_placements)
        out_clone = out.detach().clone().requires_grad_(out.requires_grad)
        grad_out_dt_clone = grad_out_dt.detach().clone().requires_grad_(grad_out_dt.requires_grad)

        out.backward(grad_out_dt)

        # V4b: FW input unchanged after backward
        _assert_unchanged(x_dt.to_local(), x_dt_clone.to_local())

        # V5: BW inputs (values only) unchanged
        assert_close(out.to_local(), out_clone.to_local(), atol=0, rtol=0)
        assert_close(grad_out_dt.to_local(), grad_out_dt_clone.to_local(), atol=0, rtol=0)

        # V10: weight gradient parity
        assert dist_emb.weight.grad is not None, "Weight gradient is None"
        weight_grad = dist_emb.weight.grad
        assert isinstance(weight_grad, DTensor), f"Weight grad should be DTensor, got {type(weight_grad)}"
        weight_grad_full = weight_grad.full_tensor()
        assert_close(
            weight_grad_full,
            weight_grad_ref_on_host.to(dm.device),
            atol=1e-5,
            rtol=1e-5,
            msg=lambda m: f"Rank {rank} weight grad mismatch\n{m}",
        )

        # V10b: replicated weight gradients identical across all CP ranks
        assert_all_identical(weight_grad_full.detach(), dm.group["cp"])

        # Non-vacuous: weight gradient must be non-zero
        assert weight_grad_full.abs().sum() > 0, "Weight gradient is all-zero"

        # Non-vacuous: verify sharding is active (local < global on at least one dim)
        assert any(
            out.to_local().shape[d] < out.shape[d] for d in range(out.ndim)
        ), f"Sharding not active: local shape {out.to_local().shape} == global shape {out.shape}"

    finally:
        DistributedManager.cleanup()
        monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env, padding_idx",
    [
        # CPU dp=1 cp=(2,2): basic parity
        (((1, (2, 2)), True, "cpu", "ENV"), None),
        # CPU dp=2 cp=(2,2): DP + CP with padding_idx
        (((2, (2, 2)), True, "cpu", "ENV"), 0),
        # CUDA dp=2 cp=(1,1): DP-only, 2-GPU
        (((2, (1, 1)), True, "cuda", "ENV"), None),
        # CUDA dp=1 cp=(2,2): actual CP, 4-GPU
        (((1, (2, 2)), True, "cuda", "ENV"), 0),
    ],
    indirect=("setup_env",),
    ids=["cpu-dp1-cp2x2", "cpu-dp2-cp2x2-pad0", "cuda-dp2-cp1x1", "cuda-dp1-cp2x2-pad0"],
)
def test_dtensor_embedding_forward_backward(setup_env, padding_idx: int | None):
    """EmbeddingParamsReplicated: distributed output and weight gradient match serial nn.Embedding."""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    skip_if_cuda_not_avail_or_device_count_less_than_word_size(device_type=device_type, world_size=world_size)

    num_embeddings = 32
    embedding_dim = 16
    B = 2 * grid_group_sizes["dp"]
    N = 8 * grid_group_sizes["cp"][0]

    with torch.random.fork_rng(devices=[], enabled=True):
        torch.manual_seed(SEED)

        serial_emb = nn.Embedding(num_embeddings, embedding_dim, padding_idx=padding_idx)
        nn.init.uniform_(serial_emb.weight, -0.5, 0.5)
        state_dict = {k: v.cpu().clone() for k, v in serial_emb.state_dict().items()}

        # Pair-like integer input [B, N, N] with diverse indices
        x = torch.randint(0, num_embeddings, (B, N, N))
        if padding_idx is not None:
            x[:, :2, :2] = padding_idx

        out_ref = serial_emb(x)
        grad_out = torch.randn_like(out_ref)
        out_ref.backward(grad_out)

        weight_grad_ref = serial_emb.weight.grad.detach().cpu().clone()

    spawn_multiprocessing(
        _worker_embedding_parity,
        world_size,
        x.cpu(),
        out_ref.detach().cpu(),
        grad_out.detach().cpu(),
        weight_grad_ref,
        state_dict,
        num_embeddings,
        embedding_dim,
        padding_idx,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
    )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
