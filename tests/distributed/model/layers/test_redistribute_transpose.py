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
import random
from typing import Dict, Optional

import pytest
import torch

from boltz.distributed.comm import TransposeComm
from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.layers.redistribute_transpose_without_dtensor import redistribute_transpose
from boltz.testing.utils import assert_tensors_identical


def seed_by_rank(rank: int, seed: int = 42) -> None:
    seed = rank + seed
    random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)


def parallel_assert_redistribute_transpose(
    rank: int,
    grid_group_sizes: Dict[str, int],
    device_type: str,
    backend: str,
    n_tokens_per_rank: int,
    input_host: torch.Tensor,
    output_expected_host: torch.Tensor,
    d_output_expected_host: torch.Tensor,
    d_input_host: torch.Tensor,
    env_map: Optional[Dict[str, str]] = None,
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

    assert input_host.device.type == "cpu"

    layout_map = manager.layout_subgroups["cp"]
    transpose_comm = TransposeComm(manager.group["cp"], layout_map)

    rank_coords = layout_map.unravel(manager.group_rank["cp"])
    i_chunk_begins = [rank_coords[i] * n_tokens_per_rank for i in range(len(rank_coords))]
    i_chunk_ends = [(rank_coords[i] + 1) * n_tokens_per_rank for i in range(len(rank_coords))]

    # Extract the chunk of the input tensor for this rank
    # The .contiguous() is necessary because the equivalent op is single-device scatter
    # of the input, which requires contiguity
    input_chunk = input_host[
        :, i_chunk_begins[0] : i_chunk_ends[0], i_chunk_begins[1] : i_chunk_ends[1], :
    ].contiguous()
    input_chunk = input_chunk.to(device=manager.device)
    input_chunk.requires_grad = True

    d_output_chunk = d_output_expected_host[
        :, i_chunk_begins[0] : i_chunk_ends[0], i_chunk_begins[1] : i_chunk_ends[1], :
    ].contiguous()
    d_output_chunk = d_output_chunk.to(device=manager.device)

    input_chunk_clone = input_chunk.clone()
    # Perform distributed transpose operation
    result = redistribute_transpose(input_chunk, 1, 2, transpose_comm)

    # must not modify the input tensor
    assert_tensors_identical(input_chunk_clone, input_chunk, check_grad_fn=False)

    # Perform backward pass
    d_output_chunk_clone = d_output_chunk.clone()
    torch.autograd.backward([result], [d_output_chunk])
    # no modification to any input
    assert_tensors_identical(d_output_chunk, d_output_chunk_clone)
    # no modification to any input, except that input_chunk now have gradients
    assert_tensors_identical(input_chunk, input_chunk_clone, check_grad_fn=False, check_grad=False)

    # Check forward pass output
    # Extract the expected chunk for this rank
    # The .contiguous() is necessary because the equivalent op is single-device transpose
    # then scatter, where the scatter op requires contiguity
    output_expected_chunk = output_expected_host[
        :, i_chunk_begins[0] : i_chunk_ends[0], i_chunk_begins[1] : i_chunk_ends[1], :
    ].contiguous()

    # Move result to host for comparison
    result_host = result.detach().to(device=output_expected_chunk.device)
    assert_tensors_identical(result_host, output_expected_chunk, check_stride=False)

    # check backward pass
    # The .contiguous() is necessary because the equivalent op is single-device transpose
    # backward then scatter, where the scatter op requires contiguity
    d_input_expected_chunk = d_input_host[
        :, i_chunk_begins[0] : i_chunk_ends[0], i_chunk_begins[1] : i_chunk_ends[1], :
    ].contiguous()
    assert_tensors_identical(input_chunk.grad.detach().cpu(), d_input_expected_chunk)

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.skip_jet_nightly
@pytest.mark.parametrize(
    "setup_env",
    itertools.product(
        [(1, (1, 1)), (2, (1, 1)), (1, (2, 2)), (2, (2, 2)), (1, (3, 3))], [True], ["cpu", "cuda"], ["ENV"]
    ),
    indirect=("setup_env",),
    ids=lambda x: f"dp={x[0][0]}, cp={x[0][1]}, specify_method={x[1]}, device_type={x[2]}, method_init={x[3]}",
)
def test_redistribute_transpose_without_dtensor(setup_env):
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    size_ring = grid_group_sizes["cp"][0]

    n_tokens_global = size_ring * size_ring
    embed_dim = 8
    batch_size = 4
    dtype = torch.float32
    device_host = torch.device("cpu")

    assert (
        n_tokens_global % size_ring == 0
    ), f"n_tokens_global {n_tokens_global} is not a multiple of size_ring {size_ring}"

    seed_by_rank(0)  # same seed for all ranks

    # Create input tensor of shape [B, N, M, D]
    input = torch.randn(
        (batch_size, n_tokens_global, n_tokens_global, embed_dim),
        dtype=dtype,
        device=device_type,
    )
    input.requires_grad = True

    # For comparison, compute the single-device transpose of dimensions 1 and 2
    output_expected = redistribute_transpose(input, 1, 2)

    d_output_expected = torch.rand_like(output_expected)
    torch.autograd.backward([output_expected], [d_output_expected])

    input_host = input.detach().clone().to(device=device_host)
    output_expected_host = output_expected.detach().clone().to(device=device_host)

    d_output_expected_host = d_output_expected.detach().clone().to(device=device_host)
    d_input_host = input.grad.detach().clone().to(device=device_host)

    n_tokens_per_rank = n_tokens_global // size_ring

    torch.multiprocessing.set_start_method("spawn", force=True)
    torch.multiprocessing.spawn(
        fn=parallel_assert_redistribute_transpose,
        args=(
            grid_group_sizes,
            device_type,
            backend,
            n_tokens_per_rank,
            input_host,
            output_expected_host,
            d_output_expected_host,
            d_input_host,
            env_per_rank,
        ),
        nprocs=world_size,
        join=True,
    )
