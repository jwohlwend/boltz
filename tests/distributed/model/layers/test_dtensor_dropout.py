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
from boltz.distributed.model.layers.dropout import apply_dropout_mask_msa_or_pair
from boltz.model.layers.dropout import get_dropout_mask
from boltz.testing.utils import (
    assert_tensors_identical,
    seed_by_rank,
    spawn_multiprocessing,
)


def parallel_assert_dropout(
    rank,
    grid_group_sizes,
    device_type,
    backend,
    env_per_rank,
    dropout,
    training,
    columnwise,
    samples_dropout_global_host,
    seed,
    src_global_host,
    output_expected_global_host,
    d_output_expected_global_host,
    d_src_expected_global_host,
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

    # Input tensors have shape (B, S, N, D) - sharded on dims 0, 1, and 2 (B, S, N)
    placements_input = (Shard(0), Shard(1), Shard(2))

    # Distribute input tensors
    src_dtensor = distribute_tensor(
        src_global_host.to(manager.device), device_mesh=manager.device_mesh_subgroups, placements=placements_input
    ).requires_grad_(True)

    if columnwise:
        placements_samples_dropout = (Shard(0), Replicate(), Shard(2))
    else:
        placements_samples_dropout = (Shard(0), Shard(1), Replicate())

    samples_dropout_dtensor = distribute_tensor(
        samples_dropout_global_host.to(manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements_samples_dropout,
    )

    # Distribute expected outputs
    output_expected_dtensor = distribute_tensor(
        output_expected_global_host.to(manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements_input,
        src_data_rank=None,
    )
    d_output_expected_dtensor = distribute_tensor(
        d_output_expected_global_host.to(manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements_input,
    )
    d_src_expected_dtensor = distribute_tensor(
        d_src_expected_global_host.to(manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements_input,
        src_data_rank=None,
    )

    # Create copies to verify inputs aren't modified
    src_dtensor_copy = src_dtensor.detach().clone().requires_grad_(True)
    samples_dropout_dtensor_copy = samples_dropout_dtensor.detach().clone()

    torch.cuda.manual_seed_all(seed)

    # Forward pass
    output_dtensor_result = apply_dropout_mask_msa_or_pair(
        src_dtensor, dropout, training, columnwise, samples_dropout_dtensor
    )

    # just so this runs but unfortunately we have no way to verify the results due to the lack of way to generate
    # consistent RNG sequences between the serial and distributed versions
    src_dtensor_no_samples_dropout = src_dtensor.detach().clone().requires_grad_(True)
    output_dtensor_result_no_samples_dropout = apply_dropout_mask_msa_or_pair(
        src_dtensor_no_samples_dropout, dropout, training, columnwise
    )
    assert output_dtensor_result_no_samples_dropout.shape == output_dtensor_result.shape

    # Verify inputs weren't modified
    assert_tensors_identical(src_dtensor_copy.to_local(), src_dtensor.to_local(), check_grad=False, check_grad_fn=False)
    assert_tensors_identical(
        samples_dropout_dtensor_copy.to_local(),
        samples_dropout_dtensor.to_local(),
        check_grad=False,
        check_grad_fn=False,
    )

    # Test forward pass results
    assert (
        output_dtensor_result.placements == placements_input
    ), f"placements_input: {placements_input}, output_dtensor_result.placements: {output_dtensor_result.placements}"
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

    # again, no way to verify the results due to the lack of way to generate
    # consistent RNG sequences between the serial and distributed versions
    output_dtensor_result_no_samples_dropout.backward(d_output_expected_dtensor)
    assert src_dtensor_no_samples_dropout.grad.shape == src_dtensor.grad.shape

    # Verify upstream gradient wasn't modified
    assert_tensors_identical(d_output_expected_dtensor_copy.to_local(), d_output_expected_dtensor.to_local())

    # Test input gradients
    assert (
        src_dtensor.grad.placements == placements_input
    ), f"placements_input: {placements_input}, src_dtensor.grad.placements: {src_dtensor.grad.placements}"
    assert (
        src_dtensor.grad.shape == d_src_expected_dtensor.shape
    ), f"Input gradient shape mismatch: {src_dtensor.grad.shape} != {d_src_expected_dtensor.shape}"
    assert (
        src_dtensor.grad.stride() == d_src_expected_dtensor.stride()
    ), f"Input gradient stride mismatch: {src_dtensor.grad.stride()} != {d_src_expected_dtensor.stride()}"
    torch.testing.assert_close(src_dtensor.grad.to_local(), d_src_expected_dtensor.to_local())

    # Verify that samples_dropout_dtensor has no gradients (should be None)
    assert samples_dropout_dtensor.grad is None, "Reference dropout samples_dropout should not have gradients"

    # Test full tensor gathering - verify distributed results match serial results
    src_global_result_host = src_dtensor.full_tensor().cpu()
    output_global_result_host = output_dtensor_result.full_tensor().cpu()
    d_src_global_result_host = src_dtensor.grad.full_tensor().cpu()

    # Verify full tensors match expected results
    torch.testing.assert_close(src_global_result_host, src_global_host)
    torch.testing.assert_close(output_global_result_host, output_expected_global_host)
    torch.testing.assert_close(d_src_global_result_host, d_src_expected_global_host)

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
@pytest.mark.parametrize("dropout", [0.0, 0.5], ids=lambda x: f"dropout={x}")
@pytest.mark.parametrize("training", [True, False], ids=lambda x: f"training={x}")
@pytest.mark.parametrize("columnwise", [True, False], ids=lambda x: f"columnwise={x}")
def test_apply_dropout_mask_msa_or_pair_parallel(setup_env, dropout, training, columnwise):
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    size_ring = grid_group_sizes["cp"][0]
    B = 2 * grid_group_sizes["dp"]
    S = size_ring * 4  # Sequence length
    N = size_ring * 4  # Number of tokens/positions
    D = 64  # Feature dimension

    seed = 42
    seed_by_rank(0, seed)

    # Create input tensors with proper 4D shape (B, S, N, D)
    src_global = torch.rand((B, S, N, D), requires_grad=True, device=device_type)
    # requires_grad but it won't be set
    z_global = torch.rand((B, S, N, D), requires_grad=True, device=device_type)

    # Run serial reference computation
    seed_by_rank(0, seed)  # Reset seed for consistent dropout mask
    src_global_copy = src_global.detach().clone().requires_grad_(True)
    z_global_copy = z_global.detach().clone().requires_grad_(True)

    mask = get_dropout_mask(dropout, z_global_copy, training, columnwise)
    output_expected_global = src_global_copy * mask.to(src_global_copy.dtype)
    output_expected_global_host = output_expected_global.detach().clone().cpu()

    # Create upstream gradient and run backward pass
    d_output_expected_global = torch.rand_like(output_expected_global)
    d_output_expected_global_host = d_output_expected_global.detach().clone().cpu()
    output_expected_global.backward(d_output_expected_global)

    # Verify that z_global_copy has no gradients in the serial version
    assert z_global_copy.grad is None, "Reference tensor z should not have gradients in serial version"

    # Get expected gradients
    d_src_expected_global_host = src_global_copy.grad.detach().clone().cpu()

    # Prepare input data for parallel test
    src_global_host = src_global.detach().clone().cpu()

    # emulate the dropout mask creation to testing the distributed version
    seed_by_rank(0, seed)
    if columnwise:
        samples_dropout_global_host = torch.rand((B, 1, N, 1), device=device_type, dtype=src_global.dtype).cpu()
    else:
        samples_dropout_global_host = torch.rand((B, S, 1, 1), device=device_type, dtype=src_global.dtype).cpu()

    # Launch parallel test across all processes
    spawn_multiprocessing(
        parallel_assert_dropout,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        dropout,
        training,
        columnwise,
        samples_dropout_global_host,
        seed,
        src_global_host,
        output_expected_global_host,
        d_output_expected_global_host,
        d_src_expected_global_host,
    )
