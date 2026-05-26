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

"""Tests for DTensor compute_aggregated_metric and compute_ptms functions.

Ported from boltz-1x-cp with import path updates for the Boltz-2 branch.
Tests verify both forward and backward passes across different device mesh
configurations and placements.
"""

from math import gcd

import pytest
import torch
from torch.distributed.tensor import Replicate, Shard, distribute_tensor

from boltz.data.feature.featurizer import BoltzFeaturizer
from boltz.data.module.inference import load_input
from boltz.data.tokenize.boltz import BoltzTokenizer
from boltz.distributed.comm import TransposeComm
from boltz.distributed.data.feature.featurizer_utils import get_num_atoms_tokens
from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.modules.confidence_utils import (
    CHAIN_IPTM_SENTINEL,
    compute_aggregated_metric,
    compute_ptms,
)
from boltz.model.layers.confidence_utils import (
    compute_aggregated_metric as serial_compute_aggregated_metric,
)
from boltz.model.layers.confidence_utils import (
    compute_ptms as serial_compute_ptms,
)
from boltz.testing.utils import (
    assert_tensors_identical,
    distribute_atom_features,
    init_tensors_uniform,
    random_features,
    spawn_multiprocessing,
)


def _assert_nontrivial_metric(metric: torch.Tensor, metric_name: str) -> None:
    """Guard against degenerate metrics that can make parity checks trivial."""
    metric_cpu = metric.detach().cpu().to(torch.float32)
    if metric_cpu.numel() == 0:
        raise AssertionError(f"{metric_name} is empty")
    if not torch.isfinite(metric_cpu).all():
        raise AssertionError(f"{metric_name} contains non-finite values")
    if torch.all(metric_cpu == 0):
        raise AssertionError(
            f"{metric_name} is trivially all zeros (likely empty effective mask support in compute_ptms)"
        )
    if torch.all(metric_cpu == 1):
        raise AssertionError(f"{metric_name} is trivially all ones")
    if metric_cpu.numel() > 1 and torch.allclose(metric_cpu, metric_cpu.reshape(-1)[0]):
        raise AssertionError(f"{metric_name} is a trivial constant value: {metric_cpu.reshape(-1)[0].item():.6f}")


# ---------------------------------------------------------------------------
# compute_aggregated_metric tests
# ---------------------------------------------------------------------------


def parallel_assert_compute_aggregated_metric(
    rank,
    grid_group_sizes,
    device_type,
    backend,
    env_per_rank,
    placements,
    input_logits_global_host,
    output_expected_global_host,
    d_output_expected_global_host,
    d_input_logits_expected_global_host,
    end: float,
):
    """Compare DTensor compute_aggregated_metric with serial reference."""
    monkeypatch = pytest.MonkeyPatch()
    if env_per_rank is not None:
        for var_name, value in env_per_rank.items():
            monkeypatch.setenv(var_name, f"{rank}" if value == "<INPUT_RANK>" else value)

    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    manager = DistributedManager()

    input_logits_dtensor = distribute_tensor(
        input_logits_global_host.to(manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements,
    ).requires_grad_(True)

    output_placements = placements
    d_output_expected_dtensor = distribute_tensor(
        d_output_expected_global_host.to(manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=output_placements,
    )
    output_expected_dtensor = distribute_tensor(
        output_expected_global_host.to(manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=output_placements,
        src_data_rank=None,
    )
    d_input_logits_expected_dtensor = distribute_tensor(
        d_input_logits_expected_global_host.to(manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements,
        src_data_rank=None,
    )

    input_logits_dtensor_copy = input_logits_dtensor.detach().clone().requires_grad_(True)

    output_dtensor_result = compute_aggregated_metric(input_logits_dtensor, end=end)

    assert_tensors_identical(
        input_logits_dtensor_copy.to_local(),
        input_logits_dtensor.to_local(),
        check_grad=False,
        check_grad_fn=False,
    )

    assert (
        output_dtensor_result.shape == output_expected_dtensor.shape
    ), f"Output shape mismatch: {output_dtensor_result.shape} vs {output_expected_dtensor.shape}"
    assert (
        output_dtensor_result.stride() == output_expected_dtensor.stride()
    ), f"Output stride mismatch: {output_dtensor_result.stride()} vs {output_expected_dtensor.stride()}"
    torch.testing.assert_close(
        output_dtensor_result.to_local(),
        output_expected_dtensor.to_local(),
    )

    d_output_expected_dtensor_copy = d_output_expected_dtensor.detach().clone()
    output_dtensor_result.backward(d_output_expected_dtensor)

    assert_tensors_identical(
        d_output_expected_dtensor_copy.to_local(),
        d_output_expected_dtensor.to_local(),
    )

    assert input_logits_dtensor.grad is not None, "Input gradient should not be None"
    assert (
        input_logits_dtensor.grad.shape == d_input_logits_expected_dtensor.shape
    ), f"Gradient shape mismatch: {input_logits_dtensor.grad.shape} vs {d_input_logits_expected_dtensor.shape}"
    assert (
        input_logits_dtensor.grad.stride() == d_input_logits_expected_dtensor.stride()
    ), f"Gradient stride mismatch: {input_logits_dtensor.grad.stride()} vs {d_input_logits_expected_dtensor.stride()}"
    torch.testing.assert_close(
        input_logits_dtensor.grad.to_local(),
        d_input_logits_expected_dtensor.to_local(),
    )

    output_global_result_host = output_dtensor_result.full_tensor().cpu()
    d_input_logits_global_result_host = input_logits_dtensor.grad.full_tensor().cpu()

    torch.testing.assert_close(
        output_global_result_host,
        output_expected_global_host,
    )
    torch.testing.assert_close(
        d_input_logits_global_result_host,
        d_input_logits_expected_global_host,
    )

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env",
    [
        ((2, (2, 2)), True, "cuda", "ENV"),
    ],
    indirect=("setup_env",),
    ids=["cuda-dp2-cp2x2"],
)
@pytest.mark.parametrize(
    "placements,input_shape_type",
    [
        ((Shard(0), Shard(1), Shard(2)), "pair"),
        ((Shard(0), Shard(1), Replicate()), "single"),
    ],
    ids=["pair-shard-all", "single-shard-B-N"],
)
@pytest.mark.parametrize("end", [1.0, 32.0], ids=["end:1.0", "end:32.0"])
def test_compute_aggregated_metric_parallel(setup_env, placements, input_shape_type, end):
    """Test compute_aggregated_metric with DTensor across multiple configurations."""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("CUDA not available")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"Need {world_size} GPUs, have {torch.cuda.device_count()}")

    size_ring = grid_group_sizes["cp"][0]
    B = 2 * grid_group_sizes["dp"]
    N = size_ring * 10

    num_bins = 50 if end == 1.0 else 64

    if input_shape_type == "pair":
        input_shape = (B, N, N, num_bins)
    else:
        input_shape = (B, N, num_bins)

    init_min, init_max = -1.0, 1.0
    output_shape = input_shape[:-1]
    input_logits_global = torch.empty(input_shape, requires_grad=True, device=device_type)
    d_output_expected_global = torch.empty(output_shape, device=device_type)
    torch.manual_seed(42)
    init_tensors_uniform([input_logits_global, d_output_expected_global], low=init_min, high=init_max)

    input_logits_global_host = input_logits_global.detach().clone().cpu()
    output_expected_global = serial_compute_aggregated_metric(input_logits_global, end=end)
    output_expected_global_host = output_expected_global.detach().clone().cpu()
    d_output_expected_global_host = d_output_expected_global.detach().clone().cpu()
    output_expected_global.backward(d_output_expected_global)

    spawn_multiprocessing(
        parallel_assert_compute_aggregated_metric,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        placements,
        input_logits_global_host,
        output_expected_global_host,
        d_output_expected_global_host,
        input_logits_global.grad.detach().clone().cpu(),
        end,
    )


def parallel_assert_compute_aggregated_metric_error_cases(
    rank,
    grid_group_sizes,
    device_type,
    backend,
    env_per_rank,
):
    """Test error cases for compute_aggregated_metric."""
    monkeypatch = pytest.MonkeyPatch()
    if env_per_rank is not None:
        for var_name, value in env_per_rank.items():
            monkeypatch.setenv(var_name, f"{rank}" if value == "<INPUT_RANK>" else value)

    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    manager = DistributedManager()

    B = 2 * grid_group_sizes["dp"]
    N = grid_group_sizes["cp"][0] * 10
    num_bins = 50

    regular_tensor = torch.empty((B, N, num_bins), device=manager.device, requires_grad=True)
    input_tensor = torch.empty((B, N, num_bins), device=manager.device, requires_grad=True)
    torch.manual_seed(42)
    init_tensors_uniform([regular_tensor, input_tensor], low=-1.0, high=1.0)

    with pytest.raises(TypeError, match="Expected DTensor"):
        compute_aggregated_metric(regular_tensor, end=1.0)

    sharded_bins_dtensor = distribute_tensor(
        input_tensor,
        device_mesh=manager.device_mesh_subgroups,
        placements=(Shard(0), Shard(1), Shard(2)),
    )
    with pytest.raises(ValueError, match="bins dimension.*must not be sharded"):
        compute_aggregated_metric(sharded_bins_dtensor, end=1.0)

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env",
    [
        ((1, (2, 2)), True, "cuda", "ENV"),
    ],
    indirect=("setup_env",),
    ids=["cuda-dp1-cp2x2"],
)
def test_compute_aggregated_metric_error_cases(setup_env):
    """Test error cases for compute_aggregated_metric function."""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("CUDA not available")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"Need {world_size} GPUs, have {torch.cuda.device_count()}")

    spawn_multiprocessing(
        parallel_assert_compute_aggregated_metric_error_cases,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
    )


def parallel_assert_compute_aggregated_metric_no_grad(
    rank,
    grid_group_sizes,
    device_type,
    backend,
    env_per_rank,
):
    """Test compute_aggregated_metric with requires_grad=False."""
    monkeypatch = pytest.MonkeyPatch()
    if env_per_rank is not None:
        for var_name, value in env_per_rank.items():
            monkeypatch.setenv(var_name, f"{rank}" if value == "<INPUT_RANK>" else value)

    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    manager = DistributedManager()

    B = 2 * grid_group_sizes["dp"]
    N = grid_group_sizes["cp"][0] * 10
    num_bins = 50

    torch.manual_seed(42)
    input_tensor = torch.empty((B, N, num_bins), device=manager.device, requires_grad=False)
    init_tensors_uniform([input_tensor], low=-1.0, high=1.0)
    input_dtensor = distribute_tensor(
        input_tensor,
        device_mesh=manager.device_mesh_subgroups,
        placements=(Shard(0), Shard(1), Replicate()),
    )

    output = compute_aggregated_metric(input_dtensor, end=1.0)

    assert output.shape == (B, N), f"Expected shape {(B, N)}, got {output.shape}"
    assert not output.requires_grad, "Output should not require grad when input doesn't"

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env",
    [
        ((2, (2, 2)), True, "cuda", "ENV"),
    ],
    indirect=("setup_env",),
    ids=["cuda-dp2-cp2x2"],
)
def test_compute_aggregated_metric_no_grad(setup_env):
    """Test compute_aggregated_metric with requires_grad=False."""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("CUDA not available")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"Need {world_size} GPUs, have {torch.cuda.device_count()}")

    spawn_multiprocessing(
        parallel_assert_compute_aggregated_metric_no_grad,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
    )


# ---------------------------------------------------------------------------
# compute_ptms tests
# ---------------------------------------------------------------------------


def parallel_assert_compute_ptms(
    rank,
    grid_group_sizes,
    device_type,
    backend,
    env_per_rank,
    input_logits_global_host,
    x_preds_global_host,
    feats_global_host,
    multiplicity: int,
    assert_nontrivial: bool,
):
    """Compare DTensor compute_ptms with serial reference on this rank's DP chunk only."""
    monkeypatch = pytest.MonkeyPatch()
    if env_per_rank is not None:
        for var_name, value in env_per_rank.items():
            monkeypatch.setenv(var_name, f"{rank}" if value == "<INPUT_RANK>" else value)

    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    manager = DistributedManager()
    device_mesh = manager.device_mesh_subgroups
    transpose_comm = TransposeComm(manager.group["cp"], manager.layout_subgroups["cp"])

    placements_pair = (Shard(0), Shard(1), Shard(2))
    placements_single = (Shard(0), Shard(1), Replicate())

    size_batch = feats_global_host["atom_pad_mask"].shape[0]
    dp_size = manager.group["dp"].size()
    local_batch_size = size_batch // dp_size
    dp_rank = manager.group_rank["dp"]
    dp_idx_str = dp_rank * local_batch_size
    dp_idx_end = dp_idx_str + local_batch_size

    chunk_logits = input_logits_global_host.to(manager.device)[dp_idx_str * multiplicity : dp_idx_end * multiplicity]
    chunk_x_preds = x_preds_global_host.to(manager.device)[dp_idx_str * multiplicity : dp_idx_end * multiplicity]
    chunk_feats = {k: v.to(manager.device)[dp_idx_str:dp_idx_end].clone() for k, v in feats_global_host.items()}
    (
        expected_ptm_chunk,
        expected_iptm_chunk,
        expected_ligand_iptm_chunk,
        expected_protein_iptm_chunk,
        expected_chain_pair_chunk,
    ) = serial_compute_ptms(chunk_logits, chunk_x_preds, chunk_feats, multiplicity)
    assert not expected_ptm_chunk.requires_grad
    assert not expected_iptm_chunk.requires_grad
    assert not expected_ligand_iptm_chunk.requires_grad
    assert not expected_protein_iptm_chunk.requires_grad
    for _idx1, chain_dict in expected_chain_pair_chunk.items():
        for _idx2, t in chain_dict.items():
            assert not t.requires_grad
    if assert_nontrivial:
        _assert_nontrivial_metric(expected_iptm_chunk.cpu(), "expected_iptm_chunk")
        _assert_nontrivial_metric(expected_ligand_iptm_chunk.cpu(), "expected_ligand_iptm_chunk")
        _assert_nontrivial_metric(expected_protein_iptm_chunk.cpu(), "expected_protein_iptm_chunk")
    else:
        if not torch.all(expected_iptm_chunk == 0):
            _assert_nontrivial_metric(expected_iptm_chunk.cpu(), "expected_iptm_chunk")
        if not torch.all(expected_ligand_iptm_chunk == 0):
            _assert_nontrivial_metric(expected_ligand_iptm_chunk.cpu(), "expected_ligand_iptm_chunk")
        if not torch.all(expected_protein_iptm_chunk == 0):
            _assert_nontrivial_metric(expected_protein_iptm_chunk.cpu(), "expected_protein_iptm_chunk")

    logits_dtensor = distribute_tensor(
        input_logits_global_host.to(manager.device),
        device_mesh=device_mesh,
        placements=placements_pair,
    )
    x_preds_unflat = x_preds_global_host.unflatten(0, (size_batch, multiplicity))
    inputs_atom = {
        "atom_counts_per_token": feats_global_host["atom_counts_per_token"].to(dtype=torch.int64),
        "atom_to_token": feats_global_host["atom_to_token"].to(dtype=x_preds_global_host.dtype),
        "atom_pad_mask": feats_global_host["atom_pad_mask"].to(dtype=x_preds_global_host.dtype),
        "atom_resolved_mask": feats_global_host["atom_resolved_mask"].to(dtype=x_preds_global_host.dtype),
        "frames_idx": feats_global_host["frames_idx"].to(dtype=torch.int64),
        "x_preds_0": x_preds_unflat[:, 0].to(dtype=x_preds_global_host.dtype),
    }
    for i_mul in range(1, multiplicity):
        inputs_atom[f"x_preds_{i_mul}"] = x_preds_unflat[:, i_mul].to(dtype=x_preds_global_host.dtype)

    placements_cp = {
        "atom_counts_per_token": (Shard(0), Replicate()),
        "atom_to_token": (Shard(0), Replicate()),
        "atom_pad_mask": (Shard(0), Replicate()),
        "atom_resolved_mask": (Shard(0), Replicate()),
        "frames_idx": (Shard(1), Replicate()),
        "x_preds_0": (Shard(0), Replicate()),
    }
    placements_dp_cp = {
        "atom_to_token": (Shard(0), Shard(1), Replicate()),
        "atom_pad_mask": (Shard(0), Shard(1), Replicate()),
        "atom_resolved_mask": (Shard(0), Shard(1), Replicate()),
        "frames_idx": (Shard(0), Shard(1), Replicate()),
        "x_preds_0": (Shard(0), Shard(1), Replicate()),
    }
    for i_mul in range(1, multiplicity):
        placements_cp[f"x_preds_{i_mul}"] = (Shard(0), Replicate())
        placements_dp_cp[f"x_preds_{i_mul}"] = (Shard(0), Shard(1), Replicate())

    feats_atom = distribute_atom_features(
        inputs_atom,
        placements_cp,
        placements_dp_cp,
        device_mesh,
        manager.group["cp"],
        multiplicities={"x_preds": multiplicity},
    )

    x_preds_dtensor = feats_atom["x_preds"]
    feats_dtensor = {
        "frames_idx": feats_atom["frames_idx"],
        "asym_id": distribute_tensor(
            feats_global_host["asym_id"].to(manager.device),
            device_mesh=device_mesh,
            placements=placements_single,
        ),
        "atom_to_token": feats_atom["atom_to_token"],
        "atom_pad_mask": feats_atom["atom_pad_mask"],
        "atom_resolved_mask": feats_atom["atom_resolved_mask"],
        "mol_type": distribute_tensor(
            feats_global_host["mol_type"].to(manager.device),
            device_mesh=device_mesh,
            placements=placements_single,
        ),
        "token_pad_mask": distribute_tensor(
            feats_global_host["token_pad_mask"].to(manager.device),
            device_mesh=device_mesh,
            placements=placements_single,
        ),
    }

    ptm, iptm, ligand_iptm, protein_iptm, chain_pair_iptm = compute_ptms(
        logits_dtensor,
        x_preds_dtensor,
        feats_dtensor,
        multiplicity,
        transpose_comm,
    )
    assert not ptm.requires_grad
    assert not iptm.requires_grad
    assert not ligand_iptm.requires_grad
    assert not protein_iptm.requires_grad
    for _idx1, chain_dict in chain_pair_iptm.items():
        for _idx2, dt in chain_dict.items():
            assert not dt.requires_grad

    torch.testing.assert_close(
        ptm.to_local().cpu().to(dtype=expected_ptm_chunk.dtype),
        expected_ptm_chunk.cpu(),
    )
    torch.testing.assert_close(
        iptm.to_local().cpu().to(dtype=expected_iptm_chunk.dtype),
        expected_iptm_chunk.cpu(),
    )
    torch.testing.assert_close(
        ligand_iptm.to_local().cpu().to(dtype=expected_ligand_iptm_chunk.dtype),
        expected_ligand_iptm_chunk.cpu(),
    )
    torch.testing.assert_close(
        protein_iptm.to_local().cpu().to(dtype=expected_protein_iptm_chunk.dtype),
        expected_protein_iptm_chunk.cpu(),
    )

    for idx1, chain_dict in expected_chain_pair_chunk.items():
        for idx2, expected_value in chain_dict.items():
            torch.testing.assert_close(
                chain_pair_iptm[idx1][idx2].to_local().cpu().to(dtype=expected_value.dtype),
                expected_value.cpu(),
            )

    for idx1, chain_dict in chain_pair_iptm.items():
        for idx2, dt in chain_dict.items():
            if idx1 not in expected_chain_pair_chunk or idx2 not in expected_chain_pair_chunk.get(idx1, {}):
                assert torch.all(
                    dt.to_local() == CHAIN_IPTM_SENTINEL
                ), f"Extra chain pair ({idx1}, {idx2}) should be sentinel {CHAIN_IPTM_SENTINEL}"

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env",
    [
        ((2, (2, 2)), True, "cuda", "ENV"),
    ],
    indirect=("setup_env",),
    ids=["cuda-dp2-cp2x2"],
)
@pytest.mark.parametrize("multiplicity", [1, 2], ids=["multiplicity=1", "multiplicity=2"])
@pytest.mark.parametrize("seed", [0, 42], ids=["seed=0", "seed=42"])
def test_compute_ptms_parallel(setup_env, multiplicity, seed):
    """Test compute_ptms with DTensor using random features."""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("CUDA not available")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"Need {world_size} GPUs, have {torch.cuda.device_count()}")

    size_ring = grid_group_sizes["cp"][0]
    batch_size = grid_group_sizes["dp"]
    n_tokens_per_rank = 20
    n_tokens = size_ring * n_tokens_per_rank
    max_atoms_per_token = 18
    n_atoms_per_rank = n_tokens_per_rank * max_atoms_per_token
    n_atoms = size_ring * n_atoms_per_rank
    num_bins = 64

    rng = torch.Generator(device=device_type)
    rng.manual_seed(seed)
    logits = torch.randn((batch_size * multiplicity, n_tokens, n_tokens, num_bins), device=device_type, generator=rng)
    x_preds = torch.randn((batch_size * multiplicity, n_atoms, 3), device=device_type, generator=rng)
    rng_features = torch.Generator(device=x_preds.device)
    rng_features.manual_seed(seed)

    feats = random_features(
        size_batch=batch_size,
        n_tokens=n_tokens,
        n_atoms=n_atoms,
        n_msa=1,
        atom_counts_per_token_range=(1, max_atoms_per_token),
        device=x_preds.device,
        float_value_range=(-1.0, 1.0),
        selected_keys=[
            "asym_id",
            "atom_to_token",
            "atom_pad_mask",
            "atom_resolved_mask",
            "atom_counts_per_token",
            "mol_type",
            "token_pad_mask",
            "frames_idx",
        ],
        rng=rng_features,
    )

    spawn_multiprocessing(
        parallel_assert_compute_ptms,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        logits.detach().cpu(),
        x_preds.detach().cpu(),
        {k: v.detach().cpu() for k, v in feats.items()},
        multiplicity,
        True,
    )


@pytest.mark.parametrize(
    "setup_env",
    [
        ((1, (2, 2)), True, "cuda", "ENV"),
    ],
    indirect=("setup_env",),
    ids=["cuda-dp1-cp2x2"],
)
@pytest.mark.parametrize("multiplicity", [1, 2], ids=["multiplicity=1", "multiplicity=2"])
@pytest.mark.parametrize("seed", [0, 42], ids=["seed=0", "seed=42"])
def test_compute_ptms_real_data_parallel(
    setup_env,
    multiplicity,
    seed,
    create_preprocessed_handle_boltz1_v1,
):
    """Test compute_ptms with DTensor using real Boltz-1 data."""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("CUDA not available")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"Need {world_size} GPUs, have {torch.cuda.device_count()}")

    processed = create_preprocessed_handle_boltz1_v1
    record = processed.manifest.records[0]
    input_data = load_input(record, processed.targets_dir, processed.msa_dir)
    tokenized = BoltzTokenizer().tokenize(input_data)
    n_atoms_raw, n_tokens_raw = get_num_atoms_tokens(tokenized)

    ring = grid_group_sizes["cp"][0]
    atoms_per_window = 32
    atom_lcm = ring * atoms_per_window // gcd(ring, atoms_per_window)
    max_atoms = ((n_atoms_raw + atom_lcm - 1) // atom_lcm) * atom_lcm
    max_tokens = ((n_tokens_raw + ring - 1) // ring) * ring
    max_seqs = ring

    featurizer = BoltzFeaturizer()
    feats_single = featurizer.process(
        tokenized,
        training=False,
        max_atoms=max_atoms,
        max_tokens=max_tokens,
        max_seqs=max_seqs,
        pad_to_max_seqs=True,
    )
    if not isinstance(feats_single, dict):
        raise TypeError("Expected non-sharded feature dict from BoltzFeaturizer.process")
    selected_keys = [
        "asym_id",
        "atom_to_token",
        "atom_pad_mask",
        "atom_resolved_mask",
        "mol_type",
        "token_pad_mask",
        "frames_idx",
    ]
    feats_single = {k: feats_single[k] for k in selected_keys}
    # v1 featurizer doesn't emit atom_counts_per_token; derive from one-hot atom_to_token
    feats_single["atom_counts_per_token"] = feats_single["atom_to_token"].sum(dim=0).to(torch.int64)

    batch_size = grid_group_sizes["dp"]
    feats = {k: v.unsqueeze(0).repeat_interleave(batch_size, dim=0) for k, v in feats_single.items()}

    n_tokens = feats["token_pad_mask"].shape[1]
    n_atoms = feats["atom_pad_mask"].shape[1]
    num_bins = 64
    rng = torch.Generator(device=device_type)
    rng.manual_seed(seed)
    logits = torch.randn((batch_size * multiplicity, n_tokens, n_tokens, num_bins), device=device_type, generator=rng)
    x_preds = torch.randn((batch_size * multiplicity, n_atoms, 3), device=device_type, generator=rng)

    spawn_multiprocessing(
        parallel_assert_compute_ptms,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        logits.detach().cpu(),
        x_preds.detach().cpu(),
        {k: v.detach().cpu() for k, v in feats.items()},
        multiplicity,
        False,
    )
