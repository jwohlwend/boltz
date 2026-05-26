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


from math import isqrt
from typing import Dict, Optional

import pytest
import torch
from torch import Tensor
from torch.distributed.device_mesh import DeviceMesh
from torch.distributed.tensor import Replicate, Shard, distribute_tensor

from boltz.distributed.comm import TransposeComm
from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.layers.atom_to_token import (
    pair_repr_token_to_atom,
    reconstruct_atom_to_token_global,
    reconstruct_r_set_to_rep_atom_global,
    reconstruct_token_to_rep_atom_global,
    single_repr_atom_to_token,
    single_repr_token_to_atom,
)
from boltz.distributed.testing.utils import create_atom_to_token_dtensor
from boltz.testing.utils import seed_by_rank, spawn_multiprocessing


def create_mock_atom_to_token_tensor(batch_size: int, n_tokens: int, n_atoms: int, cp_size: int) -> Tensor:
    """Create a mock atom_to_token one-hot tensor with diagonal block structure.

    Each atom maps to exactly one token within the same CP shard, producing a
    (B, N_atoms, N_tokens) one-hot matrix with non-zero entries only on the
    block diagonal.  Multiple atoms may map to the same token.
    """
    assert n_atoms >= n_tokens, "n_atoms must be greater than or equal to n_tokens"
    assert n_tokens % cp_size == 0 and n_atoms % cp_size == 0, "n_tokens and n_atoms must be divisible by cp_size"
    num_tokens_per_shard = n_tokens // cp_size
    num_atoms_per_shard = n_atoms // cp_size

    atom_to_token_global = torch.zeros(batch_size, n_atoms, n_tokens)  # block diagonal one-hot matrix
    for sample_idx in range(batch_size):
        for cp_idx in range(cp_size):
            start_token_idx = cp_idx * num_tokens_per_shard
            num_atoms_per_token = torch.randint(
                1, num_atoms_per_shard // num_tokens_per_shard + 1, (num_tokens_per_shard,)
            )
            atom_indices = torch.cumsum(num_atoms_per_token, dim=0)
            atom_indices = torch.clamp(atom_indices, max=num_atoms_per_shard) + cp_idx * num_atoms_per_shard
            atom_indices = [cp_idx * num_atoms_per_shard] + atom_indices.tolist()
            for token_idx_in_shard, (atom_start_idx, atom_end_idx) in enumerate(
                zip(atom_indices[:-1], atom_indices[1:])
            ):
                atom_to_token_global[sample_idx, atom_start_idx:atom_end_idx, start_token_idx + token_idx_in_shard] = 1

    return atom_to_token_global


def create_mock_token_to_rep_atom_tensor(batch_size: int, n_tokens: int, n_atoms: int, cp_size: int) -> Tensor:
    """Create a mock token_to_rep_atom one-hot tensor with diagonal block structure.

    Each token selects exactly one representative atom from the atoms belonging
    to the same CP shard, producing a (B, N_tokens, N_atoms) one-hot matrix
    with non-zero entries only on the block diagonal.
    """
    assert n_atoms >= n_tokens
    assert n_tokens % cp_size == 0 and n_atoms % cp_size == 0
    num_tokens_per_shard = n_tokens // cp_size
    num_atoms_per_shard = n_atoms // cp_size

    token_to_rep_atom_global = torch.zeros(batch_size, n_tokens, n_atoms)
    for sample_idx in range(batch_size):
        for cp_idx in range(cp_size):
            token_start = cp_idx * num_tokens_per_shard
            atom_start = cp_idx * num_atoms_per_shard
            for t in range(num_tokens_per_shard):
                rep_atom = atom_start + torch.randint(0, num_atoms_per_shard, (1,)).item()
                token_to_rep_atom_global[sample_idx, token_start + t, rep_atom] = 1

    return token_to_rep_atom_global


def create_mock_r_set_to_rep_atom_tensor(batch_size: int, n_r_set: int, n_atoms: int, cp_size: int) -> Tensor:
    """Create a mock r_set_to_rep_atom one-hot tensor with diagonal block structure.

    Each R-set element selects exactly one representative atom from the atoms
    belonging to the same CP shard, producing a (B, N_R, N_atoms) one-hot matrix
    with non-zero entries only on the block diagonal.  n_r_set must be divisible
    by cp_size.
    """
    assert n_atoms >= n_r_set
    assert n_r_set % cp_size == 0 and n_atoms % cp_size == 0
    num_r_per_shard = n_r_set // cp_size
    num_atoms_per_shard = n_atoms // cp_size

    r_set_global = torch.zeros(batch_size, n_r_set, n_atoms)
    for sample_idx in range(batch_size):
        for cp_idx in range(cp_size):
            r_start = cp_idx * num_r_per_shard
            atom_start = cp_idx * num_atoms_per_shard
            for r in range(num_r_per_shard):
                rep_atom = atom_start + torch.randint(0, num_atoms_per_shard, (1,)).item()
                r_set_global[sample_idx, r_start + r, rep_atom] = 1

    return r_set_global


def create_single_repr_token_to_atom_global_expectation(
    batch_size: int, n_tokens: int, n_atoms: int, dim: int, cp_size: int, device: torch.device
) -> tuple[Tensor, Tensor, Tensor, Tensor, Tensor]:
    """Create global tensors for single_repr_token_to_atom operation.

    Args:
        batch_size: Batch size
        n_tokens: Number of tokens
        n_atoms: Number of atoms
        dim: Feature dimension
        device: Device to place tensors on

    Returns:
        tuple: (token_single_repr_global, atom_to_token_global)
    """
    # Create global tensors
    token_repr_global = torch.randn(batch_size, n_tokens, dim, device=device, requires_grad=True)
    atom_to_token_global = create_mock_atom_to_token_tensor(batch_size, n_tokens, n_atoms, cp_size)
    atom_to_token_global = atom_to_token_global.to(device)

    # Clone inputs for distribution
    token_repr_global_clone = token_repr_global.detach().clone()
    atom_to_token_global_clone = atom_to_token_global.detach().clone()

    # Compute expected result
    result_global_expected = torch.bmm(atom_to_token_global, token_repr_global)

    # Create gradients for backward pass
    dy_global = torch.rand_like(result_global_expected)

    # Backward pass on global tensors
    result_global_expected.backward(dy_global)

    return (
        token_repr_global_clone,
        atom_to_token_global_clone,
        result_global_expected.detach().clone(),
        token_repr_global.grad.detach().clone(),
        dy_global.detach().clone(),
    )


def create_single_repr_atom_to_token_global_expectation(
    batch_size: int, n_tokens: int, n_atoms: int, dim: int, cp_size: int, device: torch.device
) -> tuple[Tensor, Tensor, Tensor, Tensor, Tensor]:
    """Create global tensors for single_repr_atom_to_token operation.

    Args:
        batch_size: Batch size
        n_tokens: Number of tokens
        n_atoms: Number of atoms
        dim: Feature dimension
        device: Device to place tensors on

    Returns:
        tuple: (atom_single_repr_global, atom_to_token_global, result_global_expected,
                atom_repr_global_grad, dy_global)
    """
    # Create global tensors
    atom_repr_global = torch.randn(batch_size, n_atoms, dim, device=device, requires_grad=True)
    atom_to_token_global = create_mock_atom_to_token_tensor(batch_size, n_tokens, n_atoms, cp_size)
    atom_to_token_global = atom_to_token_global.to(device)

    # Clone inputs for distribution
    atom_repr_global_clone = atom_repr_global.detach().clone()
    atom_to_token_global_clone = atom_to_token_global.detach().clone()

    # Compute expected result
    atom_to_token_sum = atom_to_token_global.sum(dim=1, keepdim=True) + 1e-6
    atom_to_token_mean = atom_to_token_global / atom_to_token_sum
    result_global_expected = torch.bmm(atom_to_token_mean.transpose(1, 2), atom_repr_global)

    # Create gradients for backward pass
    dy_global = torch.rand_like(result_global_expected)

    # Backward pass on global tensors
    result_global_expected.backward(dy_global)

    return (
        atom_repr_global_clone,
        atom_to_token_global_clone,
        result_global_expected.detach().clone(),
        atom_repr_global.grad.detach().clone(),
        dy_global.detach().clone(),
    )


def compute_pair_repr_token_to_atom_global_expectation(batch_size, n_tokens, n_atoms, dim, cp_size, device):
    """Compute expected results using global tensors.

    Args:
        batch_size: Batch size
        n_tokens: Number of tokens
        n_atoms: Number of atoms
        dim: Feature dimension
        device: Device to place tensors on
    Returns:
        tuple: (token_repr_global, atom_to_token_global, result_global_expected,
                token_repr_global_grad, dy_global)
    """
    # Create global tensors
    token_repr_global = torch.randn(batch_size, n_tokens, n_tokens, dim, device=device, requires_grad=True)
    atom_to_token_global = create_mock_atom_to_token_tensor(batch_size, n_tokens, n_atoms, cp_size)
    atom_to_token_global = atom_to_token_global.to(device)

    # Clone inputs for distribution
    token_repr_global_clone = token_repr_global.detach().clone()
    atom_to_token_global_clone = atom_to_token_global.detach().clone()

    # Compute expected result
    result_global_expected = torch.einsum(
        "bijd,bmi,bnj->bmnd", token_repr_global, atom_to_token_global, atom_to_token_global
    )

    # Create gradients for backward pass
    dy_global = torch.rand_like(result_global_expected)

    # Backward pass on global tensors
    result_global_expected.backward(dy_global)

    return (
        token_repr_global_clone,
        atom_to_token_global_clone,
        result_global_expected.detach().clone(),
        token_repr_global.grad.detach().clone(),
        dy_global.detach().clone(),
    )


def assert_single_repr_token_to_atom(
    rank: int,
    grid_group_sizes: Dict[str, int],
    device_type: str,
    backend: str,
    env_map: Optional[Dict[str, str]] = None,
):
    """Test distributed single_repr_token_to_atom operation in a parallel environment.

    This test validates that the single_repr_token_to_atom function produces identical
    results to the equivalent global tensor computation. It verifies:

    1. Forward pass produces the same results as global tensor computation
    2. Backward pass correctly propagates gradients through the distributed operation
    3. Results and gradients match the equivalent global tensor operations

    Args:
        rank: The process rank in the distributed environment
        grid_group_sizes: Dictionary mapping group names to their sizes for distributed setup
        device_type: Device to run the test on ("cpu" or "cuda")
        backend: The distributed backend to use (e.g., "gloo", "nccl")
        env_map: Optional dictionary of environment variables to set before initialization
    """
    monkeypatch = pytest.MonkeyPatch()
    if env_map is not None:
        for var_name, value in env_map.items():
            if value == "<INPUT_RANK>":
                monkeypatch.setenv(var_name, f"{rank}")

    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    manager = DistributedManager()

    size_cp = len(manager.group_ranks["cp"])
    size_ring = isqrt(size_cp)
    if size_ring * size_ring != size_cp:
        raise ValueError(f"cp group size {size_cp} is not a square int")

    # Set test parameters
    batch_size = 2
    n_tokens_per_rank = 4
    n_tokens_global = size_ring * n_tokens_per_rank
    n_atoms_per_rank = n_tokens_per_rank * 3
    n_atoms_global = size_ring * n_atoms_per_rank
    dim = 5

    # Set random seed based on rank for reproducibility
    seed_by_rank(0)

    # Compute global expectations
    (
        token_repr_global,
        atom_to_token_global,
        result_global_expected,
        token_repr_global_grad,
        dy_global,
    ) = create_single_repr_token_to_atom_global_expectation(
        batch_size, n_tokens_global, n_atoms_global, dim, size_ring, manager.device
    )

    # Create distributed tensors
    # token_repr: Shape (B, n_tokens, D) with placement (Shard(0), Shard(1), Replicate())
    single_repr_placements = [Shard(dim=0), Shard(dim=1), Replicate()]
    device_mesh: DeviceMesh = manager.device_mesh_subgroups
    token_repr_dtensor = distribute_tensor(token_repr_global, device_mesh, single_repr_placements)
    token_repr_dtensor.requires_grad = True

    # atom_to_token: Shape (B, n_tokens, n_atoms) with placement (Shard(0), Shard(1), Replicate())
    atom_to_token_dtensor = create_atom_to_token_dtensor(atom_to_token_global, manager.device_mesh_subgroups)

    # Compute on distributed tensors using single_repr_token_to_atom
    result_dtensor = single_repr_token_to_atom(token_repr_dtensor, atom_to_token_dtensor)

    # Distribute the upstream adjoint for backward pass
    # Expected result shape: (B, n_atoms, D)
    dy_dtensor = distribute_tensor(dy_global, device_mesh, single_repr_placements)

    # Perform backward pass
    result_dtensor.backward(dy_dtensor)

    # Create distributed tensors from global expectations for comparison
    token_repr_grad_dtensor_expected = distribute_tensor(token_repr_global_grad, device_mesh, single_repr_placements)
    result_dtensor_expected = distribute_tensor(result_global_expected, device_mesh, single_repr_placements)

    # Compare results with expected local shards
    torch.testing.assert_close(result_dtensor_expected, result_dtensor)
    torch.testing.assert_close(token_repr_grad_dtensor_expected, token_repr_dtensor.grad)

    # Test shape and stride consistency
    assert (
        result_dtensor.shape == result_dtensor_expected.shape
    ), f"Output shape mismatch: {result_dtensor.shape} != {result_dtensor_expected.shape}"
    assert (
        result_dtensor.stride() == result_dtensor_expected.stride()
    ), f"Output stride mismatch: {result_dtensor.stride()} != {result_dtensor_expected.stride()}"
    assert (
        token_repr_dtensor.grad.shape == token_repr_grad_dtensor_expected.shape
    ), f"Gradient shape mismatch: {token_repr_dtensor.grad.shape} != {token_repr_grad_dtensor_expected.shape}"
    assert (
        token_repr_dtensor.grad.stride() == token_repr_grad_dtensor_expected.stride()
    ), f"Gradient stride mismatch: {token_repr_dtensor.grad.stride()} != {token_repr_grad_dtensor_expected.stride()}"

    # Collect results as global tensors and compare with original global tensors
    result_global_result = result_dtensor.full_tensor()
    token_repr_grad_global_result = token_repr_dtensor.grad.full_tensor()

    # Assert output and input gradients match the global computation
    torch.testing.assert_close(result_global_result, result_global_expected)
    torch.testing.assert_close(token_repr_grad_global_result, token_repr_global_grad)

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
def test_single_repr_token_to_atom(setup_env):
    """Test distributed single_repr_token_to_atom operation across multiple processes."""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    spawn_multiprocessing(
        assert_single_repr_token_to_atom,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
    )


def assert_single_repr_atom_to_token(
    rank: int,
    grid_group_sizes: Dict[str, int],
    device_type: str,
    backend: str,
    env_map: Optional[Dict[str, str]] = None,
):
    """Test distributed single_repr_atom_to_token operation in a parallel environment.

    This test validates that the single_repr_atom_to_token function produces identical
    results to the equivalent global tensor computation. It verifies:

    1. Forward pass produces the same results as global tensor computation
    2. Backward pass correctly propagates gradients through the distributed operation
    3. Results and gradients match the equivalent global tensor operations

    Args:
        rank: The process rank in the distributed environment
        grid_group_sizes: Dictionary mapping group names to their sizes for distributed setup
        device_type: Device to run the test on ("cpu" or "cuda")
        backend: The distributed backend to use (e.g., "gloo", "nccl")
        env_map: Optional dictionary of environment variables to set before initialization
    """
    monkeypatch = pytest.MonkeyPatch()
    if env_map is not None:
        for var_name, value in env_map.items():
            if value == "<INPUT_RANK>":
                monkeypatch.setenv(var_name, f"{rank}")

    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    manager = DistributedManager()

    size_cp = len(manager.group_ranks["cp"])
    size_ring = isqrt(size_cp)
    if size_ring * size_ring != size_cp:
        raise ValueError(f"cp group size {size_cp} is not a square int")

    # Set test parameters
    batch_size = 2
    n_tokens_per_rank = 4
    n_tokens_global = size_ring * n_tokens_per_rank
    n_atoms_per_rank = n_tokens_per_rank * 3
    n_atoms_global = size_ring * n_atoms_per_rank
    dim = 5

    # Set random seed based on rank for reproducibility
    seed_by_rank(0)

    # Compute global expectations
    (
        atom_repr_global,
        atom_to_token_global,
        result_global_expected,
        atom_repr_global_grad,
        dy_global,
    ) = create_single_repr_atom_to_token_global_expectation(
        batch_size, n_tokens_global, n_atoms_global, dim, size_ring, manager.device
    )

    # Create distributed tensors
    # atom_repr: Shape (B, n_atoms, D) with placement (Shard(0), Shard(1), Replicate())
    single_repr_placements = [Shard(dim=0), Shard(dim=1), Replicate()]
    device_mesh: DeviceMesh = manager.device_mesh_subgroups
    atom_repr_dtensor = distribute_tensor(atom_repr_global, device_mesh, single_repr_placements)
    atom_repr_dtensor.requires_grad = True

    # atom_to_token: Shape (B, n_atoms, n_tokens) with placement (Shard(0), Shard(1), Replicate())
    atom_to_token_dtensor = create_atom_to_token_dtensor(atom_to_token_global, manager.device_mesh_subgroups)

    # Compute on distributed tensors using single_repr_atom_to_token
    result_dtensor = single_repr_atom_to_token(atom_repr_dtensor, atom_to_token_dtensor)

    # Distribute the upstream adjoint for backward pass
    # Expected result shape: (B, n_tokens, D)
    dy_dtensor = distribute_tensor(dy_global, device_mesh, single_repr_placements)

    # Perform backward pass
    result_dtensor.backward(dy_dtensor)

    # Create distributed tensors from global expectations for comparison
    atom_repr_grad_dtensor_expected = distribute_tensor(atom_repr_global_grad, device_mesh, single_repr_placements)
    result_dtensor_expected = distribute_tensor(result_global_expected, device_mesh, single_repr_placements)

    # Compare results with expected local shards
    assert (
        result_dtensor.shape == result_dtensor_expected.shape
    ), f"Output shape mismatch: {result_dtensor.shape} != {result_dtensor_expected.shape}"
    assert (
        result_dtensor.stride() == result_dtensor_expected.stride()
    ), f"Output stride mismatch: {result_dtensor.stride()} != {result_dtensor_expected.stride()}"
    torch.testing.assert_close(result_dtensor_expected, result_dtensor)

    assert (
        atom_repr_dtensor.grad.shape == atom_repr_grad_dtensor_expected.shape
    ), f"Gradient shape mismatch: {atom_repr_dtensor.grad.shape} != {atom_repr_grad_dtensor_expected.shape}"
    assert (
        atom_repr_dtensor.grad.stride() == atom_repr_grad_dtensor_expected.stride()
    ), f"Gradient stride mismatch: {atom_repr_dtensor.grad.stride()} != {atom_repr_grad_dtensor_expected.stride()}"
    torch.testing.assert_close(atom_repr_grad_dtensor_expected, atom_repr_dtensor.grad)

    # Collect results as global tensors and compare with original global tensors
    result_global_result = result_dtensor.full_tensor()
    atom_repr_grad_global_result = atom_repr_dtensor.grad.full_tensor()

    # Assert output and input gradients match the global computation
    torch.testing.assert_close(result_global_result, result_global_expected)
    torch.testing.assert_close(atom_repr_grad_global_result, atom_repr_global_grad)

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
def test_single_repr_atom_to_token(setup_env):
    """Test distributed single_repr_atom_to_token operation across multiple processes."""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    spawn_multiprocessing(
        assert_single_repr_atom_to_token,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
    )


def assert_pair_repr_token_to_atom(
    rank: int,
    grid_group_sizes: Dict[str, int],
    device_type: str,
    backend: str,
    env_map: Optional[Dict[str, str]] = None,
):
    """Test distributed atom_to_token operation in a parallel environment.

    This test validates that the pair_repr_token_to_atom function produces identical
    results to the equivalent global tensor computation. It verifies:

    1. Forward pass produces the same results as global tensor computation
    2. Backward pass correctly propagates gradients through the distributed operation
    3. Results and gradients match the equivalent global tensor operations

    Args:
        rank: The process rank in the distributed environment
        grid_group_sizes: Dictionary mapping group names to their sizes for distributed setup
        device_type: Device to run the test on ("cpu" or "cuda")
        backend: The distributed backend to use (e.g., "gloo", "nccl")
        env_map: Optional dictionary of environment variables to set before initialization
    """
    monkeypatch = pytest.MonkeyPatch()
    if env_map is not None:
        for var_name, value in env_map.items():
            if value == "<INPUT_RANK>":
                monkeypatch.setenv(var_name, f"{rank}")
                continue
            monkeypatch.setenv(var_name, value)

    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    manager = DistributedManager()

    size_cp = len(manager.group_ranks["cp"])
    size_ring = isqrt(size_cp)
    if size_ring * size_ring != size_cp:
        raise ValueError(f"cp group size {size_cp} is not a square int")

    # Set test parameters
    batch_size = 2
    n_tokens_per_rank = 4
    n_tokens_global = size_ring * n_tokens_per_rank
    n_atoms_per_rank = n_tokens_per_rank * 3
    n_atoms_global = size_ring * n_atoms_per_rank
    dim = 5

    # Set random seed based on rank for reproducibility
    seed_by_rank(0)

    # Compute global expectations
    (
        token_repr_global,
        atom_to_token_global,
        result_global_expected,
        token_repr_global_grad,
        dy_global,
    ) = compute_pair_repr_token_to_atom_global_expectation(
        batch_size, n_tokens_global, n_atoms_global, dim, size_ring, manager.device
    )

    # Create distributed tensors
    # token_repr: Shape (B, n_tokens, n_tokens, D) with placement (Shard(0), Shard(1), Shard(2))
    pair_repr_placements = [Shard(dim=0), Shard(dim=1), Shard(dim=2)]
    device_mesh: DeviceMesh = manager.device_mesh_subgroups
    token_repr_dtensor = distribute_tensor(token_repr_global, device_mesh, pair_repr_placements)
    token_repr_dtensor.requires_grad = True

    # atom_to_token: Shape (B, n_tokens, n_atoms) with placement (Shard(0), Shard(1), Replicate())
    atom_to_token_dtensor = create_atom_to_token_dtensor(atom_to_token_global, manager.device_mesh_subgroups)

    # Create TransposeComm for communication
    # The function requires a transpose communication object
    cp_group = manager.group["cp"]
    layout_group_cp = manager.layout_subgroups["cp"]
    transpose_comm = TransposeComm(cp_group, layout_group_cp)

    # Compute on distributed tensors using pair_repr_token_to_atom
    result_dtensor = pair_repr_token_to_atom(token_repr_dtensor, atom_to_token_dtensor, transpose_comm)

    # Distribute the upstream adjoint for backward pass
    # Expected result shape: (B, n_atoms, n_atoms, D)
    dy_dtensor = distribute_tensor(dy_global, device_mesh, pair_repr_placements)

    # Perform backward pass
    result_dtensor.backward(dy_dtensor)

    # Create distributed tensors from global expectations for comparison
    token_repr_grad_dtensor_expected = distribute_tensor(token_repr_global_grad, device_mesh, pair_repr_placements)
    result_dtensor_expected = distribute_tensor(result_global_expected, device_mesh, pair_repr_placements)

    # Compare results with expected local shards
    assert (
        result_dtensor.shape == result_dtensor_expected.shape
    ), f"Output shape mismatch: {result_dtensor.shape} != {result_dtensor_expected.shape}"
    # We can't guarantee same layout because two different einsum operations are used
    # in the DTensor version and the serial version
    # assert result_dtensor.stride() == result_dtensor_expected.stride(), (
    #     f"Output stride mismatch: {result_dtensor.stride()} != {result_dtensor_expected.stride()}"
    # )
    torch.testing.assert_close(result_dtensor_expected, result_dtensor)

    assert (
        token_repr_dtensor.grad.shape == token_repr_grad_dtensor_expected.shape
    ), f"Gradient shape mismatch: {token_repr_dtensor.grad.shape} != {token_repr_grad_dtensor_expected.shape}"
    assert (
        token_repr_dtensor.grad.stride() == token_repr_grad_dtensor_expected.stride()
    ), f"Gradient stride mismatch: {token_repr_dtensor.grad.stride()} != {token_repr_grad_dtensor_expected.stride()}"
    torch.testing.assert_close(token_repr_grad_dtensor_expected, token_repr_dtensor.grad)

    # Collect results as global tensors and compare with original global tensors
    result_global_result = result_dtensor.full_tensor()
    token_repr_grad_global_result = token_repr_dtensor.grad.full_tensor()

    # Assert output and input gradients match the global computation
    torch.testing.assert_close(result_global_result, result_global_expected)
    torch.testing.assert_close(token_repr_grad_global_result, token_repr_global_grad)

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
def test_pair_repr_token_to_atom(setup_env):
    """Test distributed pair_repr_token_to_atom operation across multiple processes."""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    spawn_multiprocessing(
        assert_pair_repr_token_to_atom,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
    )


def assert_reconstruct_onehot_diag_block_global(
    rank: int,
    grid_group_sizes: Dict[str, int],
    device_type: str,
    backend: str,
    env_map: Optional[Dict[str, str]] = None,
):
    """Validate all three reconstruct functions for diagonally-sharded one-hot DTensors.

    Tests:
    1. reconstruct_atom_to_token_global — round-trip matches original global tensor and
       produces correct results when used with single_repr_token_to_atom.
    2. reconstruct_token_to_rep_atom_global — round-trip matches original global tensor.
    3. reconstruct_r_set_to_rep_atom_global — round-trip matches original global tensor.
    """
    monkeypatch = pytest.MonkeyPatch()
    if env_map is not None:
        for var_name, value in env_map.items():
            if value == "<INPUT_RANK>":
                monkeypatch.setenv(var_name, f"{rank}")

    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    manager = DistributedManager()

    size_cp = len(manager.group_ranks["cp"])
    size_ring = isqrt(size_cp)
    if size_ring * size_ring != size_cp:
        raise ValueError(f"cp group size {size_cp} is not a square int")

    batch_size = grid_group_sizes["dp"]
    n_tokens_per_rank = 4
    n_tokens_global = size_ring * n_tokens_per_rank
    n_atoms_per_rank = n_tokens_per_rank * 3
    n_atoms_global = size_ring * n_atoms_per_rank
    n_r_set_per_rank = 2
    n_r_set_global = size_ring * n_r_set_per_rank
    dim = 5

    seed_by_rank(0)
    device_mesh: DeviceMesh = manager.device_mesh_subgroups

    dp_rank = device_mesh.get_coordinate()[0]
    dp_size = device_mesh.shape[0]
    local_batch = batch_size // dp_size

    # --- reconstruct_atom_to_token_global ---
    token_repr_global = torch.randn(batch_size, n_tokens_global, dim, device=manager.device)
    atom_to_token_global = create_mock_atom_to_token_tensor(batch_size, n_tokens_global, n_atoms_global, size_ring).to(
        manager.device
    )

    single_repr_placements = [Shard(dim=0), Shard(dim=1), Replicate()]
    token_repr_dtensor = distribute_tensor(token_repr_global, device_mesh, single_repr_placements)
    atom_to_token_dtensor = create_atom_to_token_dtensor(atom_to_token_global, device_mesh)

    atom_to_token_reconstructed = reconstruct_atom_to_token_global(atom_to_token_dtensor)
    atom_to_token_dp_local = atom_to_token_global[dp_rank * local_batch : (dp_rank + 1) * local_batch]
    torch.testing.assert_close(atom_to_token_reconstructed, atom_to_token_dp_local)

    result_dtensor = single_repr_token_to_atom(token_repr_dtensor, atom_to_token_dtensor)
    result_full_local = result_dtensor.redistribute(
        placements=[Shard(dim=0), Replicate(), Replicate()],
    ).to_local()
    token_repr_full_local = token_repr_dtensor.redistribute(
        placements=[Shard(dim=0), Replicate(), Replicate()],
    ).to_local()
    expected_full_local = torch.bmm(
        atom_to_token_reconstructed.to(dtype=token_repr_full_local.dtype),
        token_repr_full_local,
    )
    torch.testing.assert_close(result_full_local, expected_full_local)

    # --- reconstruct_token_to_rep_atom_global ---
    token_to_rep_atom_global = create_mock_token_to_rep_atom_tensor(
        batch_size, n_tokens_global, n_atoms_global, size_ring
    ).to(manager.device)
    token_to_rep_atom_dtensor = create_atom_to_token_dtensor(token_to_rep_atom_global, device_mesh)

    token_to_rep_atom_reconstructed = reconstruct_token_to_rep_atom_global(token_to_rep_atom_dtensor)
    token_to_rep_atom_dp_local = token_to_rep_atom_global[dp_rank * local_batch : (dp_rank + 1) * local_batch]
    torch.testing.assert_close(token_to_rep_atom_reconstructed, token_to_rep_atom_dp_local)

    # --- reconstruct_r_set_to_rep_atom_global ---
    r_set_to_rep_atom_global = create_mock_r_set_to_rep_atom_tensor(
        batch_size, n_r_set_global, n_atoms_global, size_ring
    ).to(manager.device)
    r_set_to_rep_atom_dtensor = create_atom_to_token_dtensor(r_set_to_rep_atom_global, device_mesh)

    r_set_to_rep_atom_reconstructed = reconstruct_r_set_to_rep_atom_global(r_set_to_rep_atom_dtensor)
    r_set_to_rep_atom_dp_local = r_set_to_rep_atom_global[dp_rank * local_batch : (dp_rank + 1) * local_batch]
    torch.testing.assert_close(r_set_to_rep_atom_reconstructed, r_set_to_rep_atom_dp_local)

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
def test_reconstruct_onehot_diag_block_global(setup_env):
    """Test all diagonal-block reconstruction functions: atom_to_token, token_to_rep_atom, r_set_to_rep_atom."""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    spawn_multiprocessing(
        assert_reconstruct_onehot_diag_block_global,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
    )
