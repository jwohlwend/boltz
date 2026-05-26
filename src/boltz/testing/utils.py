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

# fmt: off

import json
import math
import multiprocessing
import os
import random
import shutil
import time
import warnings
from collections import OrderedDict
from contextlib import contextmanager
from dataclasses import asdict, dataclass
from functools import partial, reduce
from pathlib import Path
from typing import Any, Callable, Optional

import numpy as np
import pytest
import requests
import torch
import torch.distributed as dist
from torch import Tensor
from torch.distributed.tensor import DeviceMesh, DTensor, Partial, Replicate, Shard
from torch.distributed.tensor._utils import compute_global_tensor_info

from boltz.data import const
from boltz.data.feature.featurizer import (
    BoltzFeaturizer,
    Tokenized,
)
from boltz.data.load import CACHE_DIR
from boltz.distributed.data.feature.featurizer_utils import (
    get_num_atoms_tokens,
    get_pair_mask,
)
from boltz.distributed.data.types import PairMaskMode

# Try to import from main, fall back to defaults if not available
try:
    from boltz.main import MODEL_URL, BoltzDiffusionParams
except ImportError:
    MODEL_URL = ""
    BoltzDiffusionParams = None

# Import from v2 modules for Boltz-2
from boltz.distributed.model.modules.utils import Precision
from boltz.distributed.utils import LayoutRightMap
from boltz.model.layers.attentionv2 import AttentionPairBias as SerialAttentionPairBias
from boltz.model.layers.pair_averaging import PairWeightedAveraging as SerialPairWeightedAveraging
from boltz.model.layers.triangular_attention.attention import TriangleAttention as SerialTriangleAttention
from boltz.model.layers.triangular_attention.attention import (
    TriangleAttentionEndingNode as SerialTriangleAttentionEndingNode,
)

# Import v2 encoder functions for window batching
from boltz.model.modules.encodersv2 import get_indexing_matrix, single_to_keys
from boltz.model.modules.transformersv2 import AtomTransformer as SerialAtomTransformer

PRECISION_TO_INF = {
    Precision.FP16: 6e4,
    Precision.BF16: 1e9,
    Precision.TF32: 1e9,
    Precision.FP32: 1e9,
    Precision.FP64: 1e18,
}


def is_a6000_gpu() -> bool:
    # Check if any of the visible GPUs is an A6000
    for i in range(torch.cuda.device_count()):
        device_name = torch.cuda.get_device_name(i)
        if "A6000" in device_name:
            return True
    return False


def download_model_ckpt() -> Path:
    """Download the model checkpoint for regression and e2e tests."""
    cache = CACHE_DIR / "regression"
    if not cache.exists():
        cache.mkdir(parents=True, exist_ok=True)
    checkpoint_url = MODEL_URL
    model_name = checkpoint_url.split("/")[-1]
    checkpoint = cache / model_name
    if not checkpoint.exists():
        download_file(checkpoint_url, checkpoint)
    return checkpoint


def map_to_device(*args, device: torch.device | str = "cpu") -> tuple[Tensor, ...]:
    return tuple(arg.to(device) if torch.is_tensor(arg) else arg for arg in args)


def get_chunk_size(N_tokens: int, ring_size: int) -> int:
    chunk_size = N_tokens / ring_size
    assert chunk_size.is_integer(), "number of tokens must be divisible by square root of context parallel size"
    return int(chunk_size)


def chunk_along_dim(*args, dim: int, chunks: int, chunk_i: int) -> Tensor | tuple[Tensor, ...]:
    if len(args) == 1:
        return args[0].chunk(chunks, dim=dim)[chunk_i]
    return tuple(t.chunk(chunks, dim=dim)[chunk_i] for t in args)


def permute_final_dims(tensor: torch.Tensor, inds: list[int]):
    zero_index = -1 * len(inds)
    first_inds = list(range(len(tensor.shape[:zero_index])))
    return tensor.permute(first_inds + [zero_index + i for i in inds])


def get_weighted_lddt(
    all_atom_pred_pos: torch.Tensor,
    all_atom_positions: torch.Tensor,
    all_atom_mask: torch.Tensor,
    cutoff: float = 15.0,
    eps: float = 1e-10,
    per_residue: bool = True,
) -> torch.Tensor:
    all_atom_mask = all_atom_mask.unsqueeze(-1)
    n = all_atom_mask.shape[-2]
    dmat_true = torch.sqrt(
        eps
        + torch.sum(
            (all_atom_positions[..., None, :] - all_atom_positions[..., None, :, :]) ** 2,
            dim=-1,
        )
    )

    dmat_pred = torch.sqrt(
        eps
        + torch.sum(
            (all_atom_pred_pos[..., None, :] - all_atom_pred_pos[..., None, :, :]) ** 2,
            dim=-1,
        )
    )
    dists_to_score = (
        (dmat_true < cutoff)
        * all_atom_mask
        * permute_final_dims(all_atom_mask, (1, 0))
        * (1.0 - torch.eye(n, device=all_atom_mask.device))
    )

    dist_l1 = torch.abs(dmat_true - dmat_pred)

    score = (
        (dist_l1 < 0.5).type(dist_l1.dtype)
        + (dist_l1 < 1.0).type(dist_l1.dtype)
        + (dist_l1 < 2.0).type(dist_l1.dtype)
        + (dist_l1 < 4.0).type(dist_l1.dtype)
    )
    score = score * 0.25

    dims = (-1,) if per_residue else (-2, -1)
    norm = 1.0 / (eps + torch.sum(dists_to_score, dim=dims))
    score = norm * (eps + torch.sum(dists_to_score * score, dim=dims))

    return score


def compute_pairwise_lddt_rmsd_matrices(
    coords_a: list[torch.Tensor],
    coords_b: list[torch.Tensor],
) -> tuple[np.ndarray, np.ndarray]:
    """Compute NxM pairwise lDDT and RMSD matrices between two sets of structures.

    Each element is a coordinate tensor of shape ``(n_atoms, 3)``.

    Returns
    -------
    lddt_matrix : np.ndarray, shape (N, M)
    rmsd_matrix : np.ndarray, shape (N, M)
    """
    from boltz.model.loss.diffusion import weighted_rigid_align

    n, m = len(coords_a), len(coords_b)
    lddt_mat = np.zeros((n, m))
    rmsd_mat = np.zeros((n, m))
    for i in range(n):
        for j in range(m):
            a = coords_a[i].unsqueeze(0)
            b = coords_b[j].unsqueeze(0)
            w = torch.ones_like(a[..., 0])
            b_aligned = weighted_rigid_align(true_coords=b, pred_coords=a, weights=w, mask=w)
            lddt_mat[i, j] = get_weighted_lddt(a, b_aligned, w).mean().item()
            rmsd_mat[i, j] = torch.sum((a - b_aligned) ** 2, dim=-1).sqrt().mean().item()
    return lddt_mat, rmsd_mat


def energy_distance_from_matrices(
    cross: np.ndarray,
    intra_a: np.ndarray,
    intra_b: np.ndarray,
    maximize: bool = False,
) -> float:
    """Compute energy distance from pre-computed pairwise metric matrices.

    For a metric where higher is better (``maximize=True``, e.g. lDDT), the
    pairwise distance is ``1 - metric``.  For lower-is-better (``maximize=False``,
    e.g. RMSD), the metric value is used directly as distance.

    Parameters
    ----------
    cross : (N, M) array -- all pairs between the two distributions.
    intra_a : (N, N) array -- all pairs within distribution A.
    intra_b : (M, M) array -- all pairs within distribution B.

    Uses the upper triangle (excluding diagonal) for intra-distribution means.
    """
    mean_cross = cross.mean()
    triu_a = intra_a[np.triu_indices(intra_a.shape[0], k=1)]
    triu_b = intra_b[np.triu_indices(intra_b.shape[0], k=1)]
    mean_intra_a = triu_a.mean()
    mean_intra_b = triu_b.mean()
    if maximize:
        return 2 * (1 - mean_cross) - (1 - mean_intra_a) - (1 - mean_intra_b)
    return 2 * mean_cross - mean_intra_a - mean_intra_b


def matched_mean_metric(matrix: np.ndarray, maximize: bool = False) -> float:
    """Optimal 1-to-1 matching (Hungarian) and return mean of matched values.

    Parameters
    ----------
    matrix : (N, M) array of metric values.
    maximize : if True, maximise the sum (e.g. lDDT); else minimise (e.g. RMSD).
    """
    from scipy.optimize import linear_sum_assignment

    cost = (1.0 - matrix) if maximize else matrix
    row_ind, col_ind = linear_sum_assignment(cost)
    return float(matrix[row_ind, col_ind].mean())


def intra_rowwise_best(matrix: np.ndarray, maximize: bool = False) -> float:
    """Mean of row-wise best value from a square matrix, excluding the diagonal.

    For lDDT (``maximize=True``) returns mean of row-wise max.
    For RMSD (``maximize=False``) returns mean of row-wise min.
    """
    mat = matrix.copy()
    if maximize:
        np.fill_diagonal(mat, -np.inf)
        return float(mat.max(axis=1).mean())
    np.fill_diagonal(mat, np.inf)
    return float(mat.min(axis=1).mean())


def download_file(url, filepath, verbose=True):
    if verbose:
        print(f"Downloading {url} to {filepath}")
    response = requests.get(url)

    target_dir = os.path.dirname(filepath)
    if target_dir and not os.path.exists(target_dir):
        os.makedirs(target_dir)

    # Check if the request was successful
    if response.status_code == 200:
        with open(filepath, "wb") as file:
            file.write(response.content)
    else:
        print(f"Failed to download file. Status code: {response.status_code}")

    return filepath


def detach_and_clone_tensors(tensors: list[Tensor | None], requires_grad: bool = False) -> list[Tensor | None]:
    return [t.detach().clone().requires_grad_(requires_grad) if t is not None else None for t in tensors]


def assert_tensors_identical(
    tensor1: torch.Tensor,
    tensor2: torch.Tensor,
    check_stride: bool = True,
    check_grad: bool = True,
    check_grad_fn: bool = True,
    check_storage_offset: bool = True,
    check_storage_pointer: bool = False,
    rtol: float = 0.0,
    atol: float = 0.0,
    **kwargs_torch_testing_assert_close: dict[str, Any],
) -> None:
    """Verify that two PyTorch tensors are identical with configurable strictness.

    Performs a multi-phase validation to ensure tensors match across different aspects:
    - Phase 1: Core tensor properties (values, device, dtype, layout, stride)
    - Phase 2: Gradient requirements
    - Phase 3: Gradient content comparison (if check_grad=True)
    - Phase 4: Autograd computation graph (if check_grad_fn=True)
    - Phase 5: Memory layout validation on storage offset (if check_storage_offset=True)
    - Phase 6: Storage pointers (if check_storage_pointer=True)

    Args:
        tensor1: First PyTorch tensor to compare
        tensor2: Second PyTorch tensor to compare
        check_stride: Whether to check that strides match
        check_grad: Whether to check that gradients match (if present)
        check_grad_fn: Whether to check that gradient functions match
        check_storage_offset: Whether to check that storage offset matches
        check_storage_pointer: Whether to check that storage pointers match
        rtol: Relative tolerance for torch.testing.assert_close
        atol: Absolute tolerance for torch.testing.assert_close
        **kwargs_torch_testing_assert_close: Additional keyword arguments for torch.testing.assert_close

    Raises:
        AssertionError: If any validation phase fails
    """
    if tensor1 is tensor2:
        return  # Short-circuit for identical objects

    # Phase 1: Core tensor properties and values
    torch.testing.assert_close(
        tensor1,
        tensor2,
        rtol=rtol,
        atol=atol,
        check_device=True,
        check_dtype=True,
        check_layout=True,
        check_stride=check_stride,
        equal_nan=True,
        **kwargs_torch_testing_assert_close,
    )

    # Phase 2: Gradient requirements
    assert tensor1.requires_grad == tensor2.requires_grad, "Input tensors' requires_grad mismatch"

    if check_grad:
        # Phase 3: Gradient content comparison
        grad1 = tensor1.grad
        grad2 = tensor2.grad

        assert (grad1 is None) == (grad2 is None), "Input tensors' gradient existence mismatch"

        if grad1 is not None and grad2 is not None:
            torch.testing.assert_close(
                grad1,
                grad2,
                rtol=rtol,
                atol=atol,
                check_device=True,
                check_dtype=True,
                check_layout=True,
                check_stride=True,
                equal_nan=True,
                **kwargs_torch_testing_assert_close,
            )

    if check_grad_fn:
        # Verify autograd graph compatibility
        assert (
            tensor1.grad_fn == tensor2.grad_fn
        ), "Autograd computation graph mismatch - Input tensors created through different operations"

    # Phase 4: Memory layout validation
    if check_storage_offset:
        assert tensor1.storage_offset() == tensor2.storage_offset(), "Input tensors' Storage offset mismatch"

    # Phase 5: Optional storage validation
    if check_storage_pointer:
        ptr1 = tensor1.storage().data_ptr()
        ptr2 = tensor2.storage().data_ptr()
        assert ptr1 == ptr2, "Input tensors' Storage pointers mismatch"


def assert_tensors_close_with_pad(
    a: torch.Tensor,
    b: torch.Tensor,
    axis: int,
    pad_val: Any = 0,
    **kwargs,
) -> None:
    """Assert that two tensors are close, handling padding along a specified axis.

    Compares the overlapping region of two tensors along the specified axis,
    and verifies that the longer tensor's trailing elements are all equal to `pad_val`.

    This is useful for comparing tensors where one has been padded to a larger size,
    e.g., after distributed_pack_and_pad operations.

    Args:
        a: First tensor to compare.
        b: Second tensor to compare.
        axis: The axis along which to compare and check padding.
        pad_val: The expected value for padding elements (default: 0).
        **kwargs: Additional keyword arguments forwarded to torch.testing.assert_close
            (e.g., atol, rtol, msg).

    Raises:
        AssertionError: If the overlapping regions don't match or trailing elements
            are not equal to pad_val.

    Example:
        >>> a = torch.tensor([1, 2, 3, 0, 0])  # padded tensor
        >>> b = torch.tensor([1, 2, 3])  # original tensor
        >>> assert_tensors_close_with_pad(a, b, axis=0, pad_val=0)  # passes
    """
    len_a = a.shape[axis]
    len_b = b.shape[axis]
    min_len = min(len_a, len_b)

    # Create slices for the overlapping region
    slices_a = [slice(None)] * a.ndim
    slices_b = [slice(None)] * b.ndim
    slices_a[axis] = slice(0, min_len)
    slices_b[axis] = slice(0, min_len)

    # Compare overlapping region
    torch.testing.assert_close(a[tuple(slices_a)], b[tuple(slices_b)], **kwargs)

    # Verify trailing padding in the longer tensor
    if len_a > len_b:
        slices_trailing = [slice(None)] * a.ndim
        slices_trailing[axis] = slice(min_len, len_a)
        trailing = a[tuple(slices_trailing)]
        expected_pad = torch.full_like(trailing, pad_val)
        torch.testing.assert_close(
            trailing,
            expected_pad,
            msg=lambda m: f"Trailing elements in tensor 'a' are not equal to pad_val={pad_val}: {m}",
        )
    elif len_b > len_a:
        slices_trailing = [slice(None)] * b.ndim
        slices_trailing[axis] = slice(min_len, len_b)
        trailing = b[tuple(slices_trailing)]
        expected_pad = torch.full_like(trailing, pad_val)
        torch.testing.assert_close(
            trailing,
            expected_pad,
            msg=lambda m: f"Trailing elements in tensor 'b' are not equal to pad_val={pad_val}: {m}",
        )


def assert_all_identical(tensor: torch.Tensor, group: torch.distributed.ProcessGroup, *args, **kwargs) -> None:
    """Verify that a tensor is identical across all processes in a distributed setup.

    Gathers the tensor from all processes in the specified process group and verifies
    that they are all identical to the input tensor using assert_tensors_identical.

    Args:
        tensor: The PyTorch tensor to verify across processes
        group: The process group to gather tensors from
        *args: Additional positional arguments passed to assert_tensors_identical
        **kwargs: Additional keyword arguments passed to assert_tensors_identical
            (e.g. check_grad, check_grad_fn, check_storage)

    Raises:
        AssertionError: If any tensor from any process differs from the input tensor
    """
    world_size = torch.distributed.get_world_size(group)
    tensor_list = [torch.empty_like(tensor) for _ in range(world_size)]
    torch.distributed.all_gather(tensor_list, tensor, group=group)
    for i in range(world_size):
        assert_tensors_identical(tensor_list[i], tensor, *args, **kwargs)


def save_gradients(model: torch.nn.Module, detach_host: Optional[bool] = False) -> dict[str, torch.Tensor | None]:
    """Save gradients of a model's parameters into a dictionary: parameter_name -> gradient."""
    grad_dict = {name: param.grad if param.grad is not None else None for name, param in model.named_parameters()}
    if detach_host:
        for name, grad in grad_dict.items():
            if grad is not None:
                grad_dict[name] = grad.detach().cpu()
    return grad_dict


def try_assert_and_collect(assertion_func: Callable, *args, error_name: str, errors_list: list[str], **kwargs) -> None:
    """
    Try to run an assertion function and collect any errors.

    Args:
        assertion_func: The assertion function to call (e.g., torch.testing.assert_close)
        *args: Positional arguments to pass to the assertion function
        error_name: Name to use in the error message
        errors_list: List to collect error messages
        **kwargs: Keyword arguments to pass to the assertion function
    """
    try:
        assertion_func(*args, **kwargs)
    except AssertionError as e:
        error_lines = str(e).strip().split("\n")
        last_three_lines = "\n".join(error_lines[-3:]) if len(error_lines) >= 3 else str(e)
        errors_list.append(f"{error_name} assertion failed: {last_three_lines}")


def repad_tensor(tensor: Tensor, pad_mask: Tensor, dim: int = 1) -> Tensor:
    """Reinsert padding into a tensor based on a padding mask.

    Args:
        tensor: Tensor without padding, shape [..., N_items_non_padded, ...]
        pad_mask: Boolean mask indicating valid (non-padding) positions, shape [N_items_total]
        dim: Dimension along which to reinsert padding

    Returns:
        Padded tensor with shape [..., N_items_total, ...]
    """
    assert pad_mask.ndim == 1

    # Get total length including padding positions
    N_total = len(pad_mask)

    # Get shape of padded tensor
    padded_shape = list(tensor.shape)
    padded_shape[dim] = N_total

    # Create padded tensor filled with zeros
    padded_tensor = torch.zeros(padded_shape, device=tensor.device, dtype=tensor.dtype)

    # Get indices of valid (non-padding) positions
    valid_indices = torch.nonzero(pad_mask).squeeze()

    # Create indexing tuples for dynamic slicing
    idx_specs = [slice(None)] * len(padded_shape)
    idx_specs[dim] = valid_indices

    # Assign the non-padded values to their correct positions in the padded tensor
    padded_tensor[tuple(idx_specs)] = tensor
    return padded_tensor


def all_gather_tensors_along_dim(tensor: Tensor, group: dist.ProcessGroup, dim: int = -1) -> Tensor:
    """All gather a tensor by concatenating along a specified dimension.

    Args:
        tensor (torch.Tensor): Tensor to all gather with a shape of (..., size, ...).
        group (dist.ProcessGroup): Process group to all gather on.
        dim (int): Dimension to concatenate along.

    Returns:
        torch.Tensor: All gathered tensor; shape = (..., size * world_size, ...)

    """
    if tensor is None:
        return None

    if tensor.requires_grad:
        raise ValueError("all_gather_tensors_along_dim breaks gradient tracking and sent tensor requires grad")

    size = dist.get_world_size(group)
    recv = [torch.empty_like(tensor) for _ in range(size)]
    dist.all_gather(recv, tensor, group=group)
    output = torch.cat(recv, dim=dim)
    return output


def all_gather_pair_repr_along_dims(
    tensor: Tensor,
    cp_group: dist.ProcessGroup,
    axis_0_size: int,
    dim0: int,
    dim1: int,
) -> Tensor:
    """All gather pair representation along two dimensions.

    Args:
        tensor: Local pair representation tensor
        cp_group: Process group for communication
        axis_0_size: Size of the first dimension of the pair representation
        dim0: First dimension of the pair representation
        dim1: Second dimension of the pair representation

    Returns:
        torch.Tensor: All gathered tensor; shape = (..., axis_0_size * world_size, axis_0_size * world_size, ...)
    """
    if tensor is None:
        return None

    if tensor.requires_grad:
        raise ValueError("all_gather_pair_repr_along_dims breaks gradient tracking and sent tensor requires grad")

    tensor_list = [torch.empty_like(tensor) for _ in range(cp_group.size())]
    dist.all_gather(tensor_list, tensor, group=cp_group)

    tensor_list = [
        torch.cat(tensor_list[i * axis_0_size : (i + 1) * axis_0_size], dim=dim1) for i in range(axis_0_size)
    ]
    return torch.cat(tensor_list, dim=dim0)


def seed_by_rank(rank: int, seed: int = 42) -> None:
    """Set random seeds based on process rank to ensure reproducible but different randomness per rank.

    This function ensures that each process in a distributed setting has its own
    deterministic random state, which is important for reproducible tests while
    maintaining appropriate randomness across different ranks.

    Args:
        rank: The process rank to use for seeding
        seed: Base seed value to which the rank is added (default: 42)
    """
    seed = rank + seed
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)


def skip_if_cuda_not_avail_or_device_count_less_than_word_size(device_type: str, world_size: int) -> None:
    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() < {world_size}")


def spawn_multiprocessing(fn: Callable[[int, ...], None], world_size: int, *args) -> None:
    """Spawn multiple processes using torch.multiprocessing for distributed testing.

    This function provides a convenient wrapper around torch.multiprocessing.spawn()
    with the spawn start method, which is commonly used for distributed PyTorch testing
    to ensure clean process isolation.

    Args:
        fn: The function to execute in each spawned process. The function must accept
            the process rank as its first argument, followed by any additional arguments
            passed via *args. Signature should be: fn(rank: int, *args) -> None
        world_size: Number of processes to spawn (typically equal to the number of GPUs
            or the desired degree of parallelism)
        *args: Additional positional arguments to pass to the spawned function

    Example:
        ```python
        def test_distributed_function(rank: int, tensor_size: int, device_prefix: str):
            device = f"{device_prefix}:{rank}"
            tensor = torch.randn(tensor_size, device=device)
            # ... distributed testing logic ...

            # Spawn 4 processes for testing
            spawn_multiprocessing(test_distributed_function, 4, 1024, "cuda")
        ```

    Note:
        - The spawn start method creates completely isolated processes, which is safer
          for testing but has more overhead than fork on Unix systems
        - Each spawned process will receive a unique rank (0 to world_size-1) as the
          first argument to the provided function
        - This function blocks until all spawned processes complete
    """
    torch.multiprocessing.set_start_method("spawn", force=True)
    torch.multiprocessing.spawn(
        fn=fn,
        args=args,
        nprocs=world_size,
        join=True,
    )


def assert_close_statistics(
    x: Tensor,
    x_ref: Tensor,
    *args,
    mean_threshold: Optional[float] = None,
    median_threshold: Optional[float] = None,
    **kwargs,
) -> None:
    # shortcircuit for empty tensors
    if x.numel() == 0:
        return

    diff = (x - x_ref).abs().cpu()
    mean = diff.mean().item()
    median = diff.median().item()
    maximum = diff.max().item()

    std_diff = diff.std().item()
    min_diff = diff.min().item()

    msg = f"\nMean: {mean:.2e}, Median: {median:.2e}, Max: {maximum:.2e}, Std: {std_diff:.2e}, Min: {min_diff:.2e}"

    if mean_threshold is not None:
        if mean > mean_threshold:
            raise AssertionError(f"Mean diff: {mean:.2e} is greater than {mean_threshold:.2e}; {msg}")
        return

    if median_threshold is not None:
        if median > median_threshold:
            raise AssertionError(f"Median diff: {median:.2e} is greater than {median_threshold:.2e}; {msg}")
        return

    # fall back to torch.testing.assert_close
    try:
        torch.testing.assert_close(x, x_ref, *args, **kwargs)
    except AssertionError as e:
        raise AssertionError(e.args[0] + msg)


def assert_absolute_or_relative_close(
    x: Tensor,
    x_ref: Tensor,
    atol: float,
    rtol: float,
) -> None:
    if x.numel() == 0:
        assert x_ref.numel() == 0, "x_ref is not empty but x is empty"
        return

    abs_diff = (x - x_ref).abs()
    rel_diff = abs_diff / torch.max(x_ref.abs(), x.abs())

    abs_pass = abs_diff <= atol
    rel_pass = rel_diff <= rtol

    either_pass = abs_pass | rel_pass
    if not either_pass.all():
        failed_abs_diffs = abs_diff[~either_pass]
        failed_rel_diffs = rel_diff[~either_pass]
        raise AssertionError(
            f"Absolute diff: {failed_abs_diffs.max():.2e} is greater than {atol:.2e} or relative diff: {failed_rel_diffs.max():.2e} is greater than {rtol:.2e}\n"
            f"Mean abs diff: {failed_abs_diffs.mean():.2e}, Median abs diff: {failed_abs_diffs.median():.2e}, Max abs diff: {failed_abs_diffs.max():.2e}\n"
            f"Mean rel diff: {failed_rel_diffs.mean():.2e}, Median rel diff: {failed_rel_diffs.median():.2e}, Max rel diff: {failed_rel_diffs.max():.2e}"
        )


def hist_diff_log10_bins(
    actual: torch.Tensor,
    expected: torch.Tensor,
    max_diff: float = 100.0,
    bin_edges: torch.Tensor | None = None,
    **kwargs_histogram,
):
    """
    Compute histogram of absolute differences between two tensors using logarithmic bins.

    This function is useful for analyzing numerical precision differences between different
    implementations (e.g., distributed vs single-device) by creating histograms of error
    magnitudes on a logarithmic scale.

    Args:
        actual (torch.Tensor): The actual tensor values to compare.
        expected (torch.Tensor): The expected tensor values to compare against.
        max_diff (float): The upper bound of the histogram bins
        bin_edges (torch.Tensor | None): The edges of the histogram bins.
            If None, the bins are created automatically based on max_diff and the dtype resolution.
        **kwargs_histogram: Additional keyword arguments passed to torch.histogram.

    Returns:
        tuple[torch.Tensor, torch.Tensor]: A tuple containing:
            - hist: Histogram counts for each bin
            - bin_edges: The edges of the logarithmic bins used for the histogram

    Raises:
        ValueError: If actual or expected are not tensors, or if they are not floating point tensors.

    Note:
        - The histogram bins are created on a logarithmic scale from the dtype resolution
          divided by 10 up to max_diff
        - The function automatically handles CUDA tensors by moving data to CPU for
          histogram computation and then back to the original device
        - The dtype used for computation is determined by the less precise of the two input tensors
    """
    if not isinstance(actual, torch.Tensor):
        raise ValueError("actual must be a tensor")

    if not isinstance(expected, torch.Tensor):
        raise ValueError("expected must be a tensor")

    if not actual.is_floating_point():
        raise ValueError("actual must be a floating point tensor")

    if not expected.is_floating_point():
        raise ValueError("expected must be a floating point tensor")

    if torch.finfo(actual.dtype).resolution > torch.finfo(expected.dtype).resolution:
        # actual is less precise than expected
        dtype = actual.dtype
    else:
        # actual is as precise as expected
        dtype = expected.dtype

    if bin_edges is None:
        # max_diff is ignored if bin_edges is provided
        base = 10
        # the lower bound of the histogram to be 1 decimal place lower than the resolution
        min_diff = torch.finfo(dtype).resolution / base
        min_bin = round(math.log10(min_diff))
        max_bin = max(round(math.log10(max_diff)), min_bin)
        n_bins = max_bin - min_bin + 1
        bin_edges = torch.logspace(start=min_bin, end=max_bin, steps=n_bins, base=base, dtype=dtype)
        # add 0 to the left-most bin
        bin_edges = torch.cat([torch.tensor([0.0], dtype=dtype), bin_edges])

    diff = (actual - expected).to(dtype).abs()
    # torch.histogram doesn't work with CUDA tensor of bins so we need to histogram on CPU
    # and move the result back to the original device. See
    # https://github.com/pytorch/pytorch/issues/69519 for details
    hist, _ = torch.histogram(
        diff.cpu(), bins=bin_edges.to(dtype=dtype, device=torch.device("cpu")), **kwargs_histogram
    )
    hist = hist.to(diff.device)

    return hist, bin_edges


def pretty_prints_hist(
    hist: torch.Tensor, bin_edges: torch.Tensor, do_cumsum: bool = False, convert_to_percentage: bool = True
) -> str:
    """
    Pretty print a histogram in a tabular format. Return the string as if printed to stdout.

    This function prints a formatted table showing histogram data with bin ranges
    and corresponding counts/percentages. It's useful for visualizing error distributions
    and numerical precision analysis.

    Args:
        hist (torch.Tensor): Histogram counts for each bin.
        bin_edges (torch.Tensor): The edges of the histogram bins. Should have one more
            element than hist.
        do_cumsum (bool, optional): If True, display cumulative sum of histogram values.
            Defaults to True.
        convert_to_percentage (bool, optional): If True, convert histogram counts to
            percentages of total. Defaults to True.

    Returns:
        str: The formatted string as if printed to stdout.

    Raises:
        ValueError: If hist or bin_edges are not tensors.

    Example:
        The output format shows bin ranges in the first two rows and histogram
        values in the last row:

        |  1.0e-07 |  1.0e-06 |  1.0e-05 | ...
        | - 1.0e-06| - 1.0e-05| - 1.0e-04| ...
        +----------+----------+----------+...
        | 2.5e+01% | 4.5e+01% | 7.8e+01% | ...

    Note:
        - The function prints directly to stdout using print statements
        - Scientific notation is used for both bin edges and histogram values
        - Bin ranges are shown as "lower_bound - upper_bound" format
        - Values are displayed as percentages when convert_to_percentage=True
    """
    if not isinstance(hist, torch.Tensor):
        raise ValueError("hist must be a tensor")

    if not isinstance(bin_edges, torch.Tensor):
        raise ValueError("bin_edges must be a tensor")

    ans = "\t| " + " | ".join([f"{b.item():8.0e}" for b in bin_edges[:-1]]) + " |"
    ans += "\n"
    ans += "\t| " + " | ".join([f"- {b.item():6.0e}" for b in bin_edges[1:]]) + " |"
    ans += "\n"
    ans += "\t" + "+----------" * (bin_edges.numel() - 1) + "+"
    ans += "\n"
    if convert_to_percentage:
        total = hist.sum()
        hist = hist / total
    if do_cumsum:
        hist = hist.cumsum(dim=0)

    if convert_to_percentage:
        ans += "\t| " + " | ".join([f"{(p * 100.0):7.1e}%" for p in hist]) + " |"
    else:
        ans += "\t| " + " | ".join([f"{p:8.1e}" for p in hist]) + " |"

    return ans


def get_param_by_key(module: torch.nn.Module, key_state_dict: str) -> Any:
    """
    Retrieve a parameter or attribute from a PyTorch module using dot notation.

    This function traverses a module's nested structure using a string key with dot
    notation to access deeply nested parameters or attributes. It's useful for
    programmatically accessing specific layers, weights, or other attributes in
    complex neural network architectures.

    Args:
        module (torch.nn.Module): The PyTorch module to search within.
        key_state_dict (str): The dot-separated path to the desired parameter or
            attribute (e.g., "encoder.layer.0.attention.query.weight").

    Returns:
        Any: The parameter, tensor, or attribute found at the specified path.

    Raises:
        ValueError: If module is not a torch.nn.Module, if key_state_dict is not
            a string, or if key_state_dict is empty.
        AttributeError: If the specified path does not exist in the module.

    Example:
        >>> import torch.nn as nn
        >>> model = nn.Sequential(nn.Linear(10, 5), nn.ReLU(), nn.Linear(5, 1))
        >>> weight = get_param_by_key(model, "0.weight")
        >>> bias = get_param_by_key(model, "2.bias")

    Note:
        - The function uses Python's reduce() with getattr() to traverse the path
        - Each part of the dot-separated key must be a valid attribute name
        - The function can access any attribute, not just parameters (weights/biases)
    """
    if not isinstance(module, torch.nn.Module):
        raise ValueError(f"module is not a torch.nn.Module: {type(module)}")

    if not isinstance(key_state_dict, str):
        raise ValueError(f"key_state_dict is not a string: {type(key_state_dict)}")

    if key_state_dict == "":
        raise ValueError("key_state_dict is empty")

    names = key_state_dict.split(sep=".")
    return reduce(getattr, names, module)


def assert_no_percentile_upshift(
    result: torch.Tensor,
    expected: torch.Tensor,
    alternative: torch.Tensor,
    perc: OrderedDict[float, tuple[float, float]] | None = None,
    names_input: tuple[str, str, str] | None = None,
):
    """
    Assert that result tensor accuracy doesn't significantly degrade compared to alternative.

    This function validates that a computation result (e.g., distributed computation)
    maintains comparable numerical accuracy to an alternative implementation (e.g.,
    single-precision reference) when both are compared against a high-precision expected
    value. It uses percentile-based analysis to detect accuracy degradation.

    Args:
        result (torch.Tensor): The tensor result to validate (e.g., from distributed computation).
        expected (torch.Tensor): The high-precision reference tensor (ground truth).
        alternative (torch.Tensor): The alternative implementation result to compare against
            (e.g., single-precision baseline).
        perc (OrderedDict[float, tuple[float, float]] | None, optional): Dictionary mapping percentile
            values to (atol, rtol) tolerance tuples to be used by torch.testing.assert_close() to
            check the consistency of the input percentiles (keys of the dict) of errors between the
            result and alternative. If None, defaults to tolerances for 50th, 75th,
            and 95th percentiles with assert_close()'s default (atol, rtol)
        names_input (tuple[str, str, str] | None, optional): Names of the input tensors to be used in the error message.
            If None, no names will be included in the error message.

    Raises:
        AssertionError: If the result shows significant upward shift in error percentiles
            compared to the alternative, indicating accuracy degradation.

    Note:
        - The function allows downward shifts (better accuracy)
        - Error histograms are generated using logarithmic bins for detailed analysis
        - Assertion errors include histogram visualizations for debugging
        - This is particularly useful for validating distributed computing implementations
          maintain numerical stability compared to single-device baselines

    Example:
        >>> # Validate distributed computation accuracy
        >>> assert_no_percentile_upshift(
        ...     result=distributed_output,
        ...     expected=fp64_reference,
        ...     alternative=fp32_reference
        ... )
    """
    if perc is None:
        perc = OrderedDict({0.25: (None, None), 0.5: (None, None), 0.75: (None, None), 0.95: (None, None)})

    diff_abs_result = (result - expected).abs()
    diff_abs_alternative = (alternative - expected).abs()

    # torch.quantile requires input dtype to be at least fp32
    diff_abs_result = diff_abs_result.to(dtype=torch.promote_types(torch.float32, diff_abs_result.dtype))
    diff_abs_alternative = diff_abs_alternative.to(dtype=torch.promote_types(torch.float32, diff_abs_alternative.dtype))

    percentages = list(perc.keys())
    quantiles_diff_abs_result = torch.quantile(
        diff_abs_result, torch.tensor(percentages, device=result.device, dtype=diff_abs_result.dtype)
    )
    quantiles_diff_abs_alternative = torch.quantile(
        diff_abs_alternative, torch.tensor(percentages, device=alternative.device, dtype=diff_abs_alternative.dtype)
    )

    max_diff = max(diff_abs_result.max().item(), diff_abs_alternative.max().item())
    if max_diff == 0:  # shortcircuit for zero-diff tensors to avoid math.log domain error
        return

    hist_result, bin_edges = hist_diff_log10_bins(result, expected, max_diff=max_diff)
    hist_alternative, _ = hist_diff_log10_bins(alternative, expected, bin_edges=bin_edges)

    str_hist_result = pretty_prints_hist(hist_result, bin_edges)
    str_hist_alternative = pretty_prints_hist(hist_alternative, bin_edges)

    error_msg = "\n"
    if names_input is not None:
        error_msg += f"input: {names_input[0]}\n"
        error_msg += f"expected: {names_input[1]}\n"
        error_msg += f"alternative: {names_input[2]}\n"
        error_msg += "\n"
        error_msg += f"hist(result, ref):\n\n{str_hist_result}\n"
        error_msg += "\n"
        error_msg += f"hist(alternative, ref):\n\n{str_hist_alternative}\n"

    for i, (p, (atol, rtol)) in enumerate(perc.items()):
        if quantiles_diff_abs_result[i] < quantiles_diff_abs_alternative[i]:
            # Assumption: we don't care down-shifting of because they only implies
            # overall higher consistency of the result than the alternative
            continue
        torch.testing.assert_close(
            quantiles_diff_abs_result[i],
            quantiles_diff_abs_alternative[i],
            atol=atol,
            rtol=rtol,
            msg=lambda m: (
                f"""
                shift at {p * 100} percentile:\n{m}\n
                {error_msg}
                """
            ),
        )


def init_tensors_uniform(tensors: list[torch.Tensor], low: float = 0.0, high: float = 1.0) -> None:
    """Initialize a list of tensors with uniform distribution in-place.

    Args:
        tensors: List of tensors to initialize
        low: Lower bound for uniform distribution (default: 0.0)
        high: Upper bound for uniform distribution (default: 1.0)
    """
    with torch.no_grad():
        for tensor in tensors:
            tensor.uniform_(low, high)


def init_tensors_normal(tensors: list[torch.Tensor], *args, **kwargs) -> None:
    """Initialize a list of tensors with values drawn from a normal distribution.

    Args:
        tensors: List of tensors to initialize in-place.
        *args: Positional arguments forwarded to tensor.normal_().
        **kwargs: Keyword arguments forwarded to tensor.normal_().
                  Common kwargs: mean (default 0), std (default 1).
    """
    with torch.no_grad():
        for tensor in tensors:
            tensor.normal_(*args, **kwargs)


def init_module_params_uniform(module: torch.nn.Module, low: float = 0.0, high: float = 1.0) -> None:
    """Initialize all named parameters of a module with uniform distribution in-place.

    Args:
        module: PyTorch module whose parameters to initialize
        low: Lower bound for uniform distribution (default: 0.0)
        high: Upper bound for uniform distribution (default: 1.0)
    """
    with torch.no_grad():
        for name, param in module.named_parameters():
            param.uniform_(low, high)


def init_module_params_glorot(module: torch.nn.Module, gain: float = 1.0) -> None:
    """Initialize parameters with Xavier/Glorot uniform distribution in-place.

    Weight tensors (dim >= 2) use ``xavier_uniform_`` scaled by *gain*.
    Bias / 1-D tensors use ``uniform_(-b, b)`` with ``b = 1/sqrt(fan_out)``.

    Args:
        module: PyTorch module whose parameters to initialize.
        gain: Multiplicative scaling factor for ``xavier_uniform_``
            (default: 1.0).
    """
    with torch.no_grad():
        for _name, param in module.named_parameters():
            if param.dim() >= 2:
                torch.nn.init.xavier_uniform_(param, gain=gain)
            else:
                bound = 1.0 / math.sqrt(param.shape[0]) if param.shape[0] > 0 else 0.01
                param.uniform_(-bound, bound)


def set_dtype_specific_inf_values(module, dtype: torch.dtype) -> None:
    """Set dtype-specific inf values for attention modules named
    "tri_att_start", "tri_att_end", and "pair_weighted_averaging".

    Parameters
    ----------
    module : torch.nn.Module
        The module containing attention modules (e.g., MSAModule, MSALayer,
        PairformerLayer, PairformerModule, tri_att_start, tri_att_end,
        pair_weighted_averaging, etc.)
    dtype : torch.dtype
        The data type to determine the appropriate inf value
    """
    dtype_to_inf = {torch.float32: 1e9, torch.float64: 1e18}
    inf_value = dtype_to_inf.get(dtype, 1e9)

    # Handle MSAModule (contains multiple MSALayers)
    if hasattr(module, "layers"):
        for layer in module.layers:
            # Handle both checkpoint-wrapped and non-wrapped layers
            l = layer._checkpoint_wrapped_module if hasattr(layer, "_checkpoint_wrapped_module") else layer
            _set_layer_inf_values(l, inf_value)
    # Handle single MSALayer
    else:
        _set_layer_inf_values(module, inf_value)


def _set_layer_inf_values(layer, inf_value: float) -> None:
    """Helper function to set inf values on a single layer.

    Parameters
    ----------
    layer : torch.nn.Module
        The layer containing attention modules (tri_att_start, tri_att_end,
        pair_weighted_averaging, etc.)
    inf_value : float
        The inf value to set
    """
    if hasattr(layer, "tri_att_start"):
        layer.tri_att_start.inf = inf_value
    if hasattr(layer, "tri_att_end"):
        layer.tri_att_end.inf = inf_value
    if hasattr(layer, "pair_weighted_averaging"):
        layer.pair_weighted_averaging.inf = inf_value
    if hasattr(layer, "attention"):
        layer.attention.inf = inf_value


class SetModuleInfValues:
    """A callable class that automatically sets dtype-specific infinity values for attention modules.

    This class is designed to be used as a function that can be applied to PyTorch modules
    to automatically configure their infinity values based on their underlying data types.
    It's particularly useful for attention mechanisms where infinity values are used for
    masking operations and need to be adjusted based on the precision of the computations.

    The class supports several attention module types and automatically detects their
    data types by examining specific linear layers within each module. It then sets
    appropriate infinity values that are compatible with the detected dtype to avoid
    numerical overflow or underflow issues.

    Attributes:
        dtype_to_inf (dict): Mapping from PyTorch data types to appropriate infinity values.
            - torch.float32: 1e9 (to avoid overflow in single precision)
            - torch.float64: 1e18 (higher precision allows larger values)
        module_types_to_layer_for_dtype (dict): Mapping from module types to the name
            of the layer used to determine the module's data type. This is used to
            automatically detect the dtype by examining the weight tensor of a specific
            linear layer within each module type.

    Supported Module Types:
        - SerialAttentionPairBias: Uses 'proj_q' layer for dtype detection
        - SerialPairWeightedAveraging: Uses 'proj_m' layer for dtype detection
        - SerialTriangleAttention: Uses 'linear' layer for dtype detection
        - SerialTriangleAttentionEndingNode: Uses 'linear' layer for dtype detection

    Usage:
        >>> inf_setter = SetModuleInfValues()
        >>> model.apply(inf_setter)  # Apply to all modules in a model
        >>> # Or apply to a specific module
        >>> inf_setter(attention_module)

    Raises:
        AttributeError: If a supported module type doesn't have the expected 'inf' attribute.
        TypeError: If a module contains an unsupported data type.

    Note:
        This class is typically used in testing scenarios where consistent infinity
        values across different data types are required for numerical stability and
        reproducible results.
    """

    def __init__(self):
        self.dtype_to_inf = {torch.float32: 1e9, torch.float64: 1e18}
        # find each module type's layer to determine its inherent dtype
        self.module_types_to_layer_for_dtype = {
            SerialAttentionPairBias: "proj_q",
            SerialPairWeightedAveraging: "proj_m",
            SerialTriangleAttention: "linear",
            SerialTriangleAttentionEndingNode: "linear",
        }

    def __call__(self, module: torch.nn.Module) -> None:
        """Set appropriate infinity value for the given module based on its data type.

        This method examines the module type and, if it's a supported attention module,
        determines its data type by inspecting a specific linear layer. It then sets
        the module's 'inf' attribute to an appropriate value based on the detected dtype.

        Args:
            module (torch.nn.Module): The PyTorch module to process. If the module
                type is not in the supported list, the method returns early without
                making any changes.

        Returns:
            None: This method modifies the module in-place by setting its 'inf' attribute.

        Raises:
            AttributeError: If a supported module type doesn't have the expected 'inf'
                attribute that should be set.
            TypeError: If the module's detected data type is not supported (not in
                dtype_to_inf mapping).

        Note:
            This method is designed to be used with PyTorch's Module.apply() method
            for automatic application across an entire model hierarchy.
        """
        type_module = type(module)
        if type_module not in self.module_types_to_layer_for_dtype:
            return

        if not hasattr(module, "inf"):
            raise AttributeError(f"Module {type_module.__name__} should but does not have an 'inf' attribute")

        dtype_module = getattr(module, self.module_types_to_layer_for_dtype[type_module]).weight.dtype
        if dtype_module not in self.dtype_to_inf:
            raise TypeError(f"Unsupported dtype {dtype_module} found in module {type_module.__name__}")
        inf_value = self.dtype_to_inf[dtype_module]
        setattr(module, "inf", inf_value)


class FixBoltzMultiplicityBug:
    """A callable class that fixes the Boltz multiplicity bug by setting reorder_pair_repr_multiplex=True.

    This class is designed to be used with PyTorch's Module.apply() method to recursively
    traverse a model and set the `reorder_pair_repr_multiplex` attribute to True on
    `AtomTransformer` and `AttentionPairBias` modules.

    This fixes the bug described in Boltz github commit 4fa0d0a0c3090ca09e71073fdd58e4108c517382,
    where the pair representation multiplicity was applied in the wrong order, causing
    incorrect behavior when multiplicity > 1.

    Supported Module Types:
        - AtomTransformer: Sets reorder_pair_repr_multiplex = True
        - AttentionPairBias: Sets reorder_pair_repr_multiplex = True

    Usage:
        >>> bug_fixer = FixBoltzMultiplicityBug()
        >>> model.apply(bug_fixer)  # Apply to all modules in a model

    Note:
        This class is typically used in testing scenarios to ensure the multiplicity
        bug is fixed for numerical correctness validation.
    """

    def __init__(self):
        # Module types that need the fix
        self.module_types_to_fix = (SerialAtomTransformer, SerialAttentionPairBias)

    def __call__(self, module: torch.nn.Module) -> None:
        """Set reorder_pair_repr_multiplex=True for supported module types.

        Args:
            module (torch.nn.Module): The PyTorch module to process. If the module
                type is not in the supported list, the method returns early without
                making any changes.

        Returns:
            None: This method modifies the module in-place.
        """
        if isinstance(module, self.module_types_to_fix):
            if hasattr(module, "reorder_pair_repr_multiplex"):
                if isinstance(module, SerialAttentionPairBias):
                    # The token variant of AttentionPairBias can't apply the fix
                    # because otherwise the pairbias z will have mismatching shape
                    # since DiffusionTransformer (token level) never apply the multiplicity
                    # to "z" unlike AtomTransformer (atom level) does.
                    module.reorder_pair_repr_multiplex = module.use_window_batching
                else:
                    module.reorder_pair_repr_multiplex = True


def get_window_batch_key_indices(n_atoms_no_pad: int, W: int, H: int) -> torch.Tensor:
    """Generate key indices for window-based attention batching.

    Creates an indexing matrix that maps global atom positions to key positions
    within each window for window-based attention. The function pads the sequence
    to the next multiple of W and generates indices for H key positions per window.

    Args:
        n_atoms_no_pad (int): Number of atoms without padding.
        W (int): Window size (number of queries per window).
        H (int): Number of key positions per window.

    Returns:
        torch.Tensor: Index tensor of shape (K, H) where K = max_atoms // W.
            Each row contains the key indices for one window. Padded positions
            are represented by 0.
    """
    if n_atoms_no_pad % W == 0:
        max_atoms = n_atoms_no_pad
    else:
        # pad to the next multiple of W
        max_atoms = ((n_atoms_no_pad // W) + 1) * W
    # construct pair mask through indexing matrices
    # TODO construct pair mask directly from AF3 appendix
    index = torch.arange(1, max_atoms + 1)
    index[n_atoms_no_pad:] = 0
    index = index.unsqueeze(0)

    K = max_atoms // W
    keys_indexing_matrix = get_indexing_matrix(K, W, H, index.device)
    to_keys = partial(single_to_keys, indexing_matrix=keys_indexing_matrix, W=W, H=H)
    index_keys = to_keys(index.unsqueeze(-1).float()).view(K, H).long()
    return index_keys


def _pair_masked_global_to_window_batch(
    pair_masked_global: torch.Tensor, n_atoms_no_pad: int, W: int = 32, H: int = 128
) -> torch.Tensor:
    """Convert a global pair representation to window-batched format.

    Transforms a global pairwise interaction matrix into a window-batched format
    suitable for efficient window-based attention computation. The function uses
    sparse matrix operations to efficiently handle the transformation while
    accounting for padding and window alignment.

    The transformation involves:
    1. Computing key indices for each window
    2. Converting to sparse CSR format for efficient manipulation
    3. Adjusting column indices to account for window-specific padding
    4. Reconstructing as a dense tensor in window-batched format

    Args:
        pair_masked_global (torch.Tensor): Input global pair representation.
            Shape can be (N, N) for 2D or (N, N, D) for 3D with embedding dimension.
            Must be square in the first two dimensions. Already masked with zeros
            for invalid padding positions.
        n_atoms_no_pad (int): Number of atoms without padding.
        W (int, optional): Window size. Defaults to 32.
        H (int, optional): Number of keys per window. Defaults to 128.

    Returns:
        torch.Tensor: Window-batched representation of shape (K, W, H) for 2D input
            or (K, W, H, D) for 3D input, where:
            - K = number of windows
            - W = window size
            - H = number of keys per window
            - D = embedding dimension (if present)

    Raises:
        ValueError: If pair_masked_global is not square in the first two dimensions.
        AssertionError: If pair_masked_global is not 2D or 3D.

    Note:
        For batch processing, call this function separately for each batch element.
    """
    # Assumption 0: input pair_masked_global is already masked with invalid padding represented by zeros
    # valid values of zeros are safe to use because this function doesn't rely on the input element values.
    ids_keys_per_window = get_window_batch_key_indices(n_atoms_no_pad, W, H)

    # The limitation of using CSR matrix is that each element along the batch
    # dimension must have the same nnz. This won't be useful for converting mask
    # in a batch but we must loop over the batch dimension and call each entry
    if pair_masked_global.shape[0] != pair_masked_global.shape[1]:
        raise ValueError("pair_masked_global must be square")
    assert pair_masked_global.ndim == 2 or pair_masked_global.ndim == 3, (
        "pair_masked_global must be 2D with potential one trailing embedding dimension. "
        "For batch dimension, loop over the batch dimension and call this function for each entry"
    )
    # count number of leading zeros per window, which can only be >= 0.
    # Left padding > 0 means that we need to increase the column index
    # while == 0 means we need to reset the column offset so that
    # the corresponding rows are left-aligned in the resulting global matrix.
    # This is to match the window-batching behavior of the attention.
    n_left_padding_per_window = (ids_keys_per_window.cumsum(dim=1) == 0).sum(dim=1)
    n_windows_no_pad = ids_keys_per_window.shape[0]
    # sparse map generation
    masked_csr = pair_masked_global.to_sparse_csr(dense_dim=1 if pair_masked_global.ndim == 3 else None)
    crow_ids = masked_csr.crow_indices()
    col_ids_new = masked_csr.col_indices().detach().clone()
    for i_window in range(n_windows_no_pad):
        i_rows_begin = i_window * W
        i_rows_end = min(i_rows_begin + W, crow_ids.shape[0] - 1)
        inz_begin = crow_ids[i_rows_begin]
        inz_end = crow_ids[i_rows_end]
        n_left_padding_this_window = n_left_padding_per_window[i_window]
        if n_left_padding_this_window > 0:
            col_ids_new[inz_begin:inz_end] += n_left_padding_this_window
        else:
            # equivalent to:
            # col_id_min = ids_keys_per_window[i_window].min().item()
            col_id_min = col_ids_new[inz_begin:inz_end].min().item()
            col_ids_new[inz_begin:inz_end] -= col_id_min
    # dim 1 is always the number of atoms regardless of ndim == 2 or 3
    n_atoms_padded = pair_masked_global.shape[1]
    # pad to the next multiple of W towards the end of both dimensions
    n_windows_padded = (n_atoms_padded + W - 1) // W
    target_length = n_windows_padded * W + 1
    current_length = crow_ids.shape[0]
    if target_length > current_length:
        padding_length = target_length - current_length
        last_value = crow_ids[-1]
        padding = last_value.repeat(padding_length)
        crow_ids_new = torch.cat([crow_ids, padding])
    else:
        crow_ids_new = crow_ids
    masked_csr_new = torch.sparse_csr_tensor(
        crow_ids_new,
        col_ids_new,
        masked_csr.values(),
        size=(n_windows_padded * W, H, pair_masked_global.shape[-1])
        if pair_masked_global.ndim == 3
        else (n_windows_padded * W, H),
        dtype=masked_csr.dtype,
        device=masked_csr.device,
    )
    masked_new = masked_csr_new.to_dense().unflatten(0, (n_windows_padded, W))
    return masked_new


@torch.no_grad()
def pair_global_to_window_batch(
    pair_repr_global: torch.Tensor,
    n_atoms_no_pads: torch.Tensor,
    pair_mask_global: torch.Tensor | None = None,
    W: int = 32,
    H: int = 128,
) -> torch.Tensor:
    """Convert batched global pair representations to window-batched format.

    Applies masking to global pair representations and converts each batch element
    to window-batched format. This is a convenience wrapper around
    pair_masked_global_to_window_batch that handles batching and masking.

    Args:
        pair_repr_global (torch.Tensor): Global pair representations with shape
            (B, N, N, D) where B is batch size, N is sequence length, and D is
            embedding dimension.
        n_atoms_no_pads (torch.Tensor): Number of atoms without padding for each batch element.
        pair_mask_global (torch.Tensor | None, optional): Global pair mask with shape (B, N, N) or
            (B, N, N, 1). Used to mask invalid positions. If None, no masking is applied.
            Defaults to None.
        W (int, optional): Window size. Defaults to 32.
        H (int, optional): Number of keys per window. Defaults to 128.

    Returns:
        torch.Tensor: Window-batched pair representations of shape (B, K, W, H, D)
            where K = number of windows, W = window size, H = keys per window,
            and D = embedding dimension.

    Note:
        When pair_mask_global is provided, pair_repr_global and pair_mask_global
        must be broadcastable for element-wise multiplication.
    """
    size_batch = pair_repr_global.shape[0]
    if n_atoms_no_pads.shape != (size_batch,):
        raise ValueError(f"n_atoms_no_pads must be a 1D tensor of size {size_batch}")
    if not (n_atoms_no_pads.dtype == torch.int32 or n_atoms_no_pads.dtype == torch.int64):
        raise ValueError("n_atoms_no_pads must be an int32 or int64 tensor")
    # the product must be broadcastable
    if pair_mask_global is None:
        ans_per_batch = [
            _pair_masked_global_to_window_batch(pair_repr_global[i], n_atoms_no_pads[i], W, H)
            for i in range(size_batch)
        ]
    else:
        pair_repr_global_masked = pair_repr_global * pair_mask_global
        ans_per_batch = [
            _pair_masked_global_to_window_batch(pair_repr_global_masked[i], n_atoms_no_pads[i], W, H)
            for i in range(size_batch)
        ]
    # the assumption here is that _pair_masked_global_to_window_batch will
    # pad to the window-batched result to have the same number of windows for each batch element
    # so we can always stack the results along the batch dimension
    ans = torch.stack(ans_per_batch, dim=0)
    return ans


def get_features(
    tokenized: Tokenized,
    window_batching: bool,
    shard_dims: Optional[tuple[int, int]] = None,
    selected_keys: Optional[list[str]] = None,
) -> dict[str, Tensor] | list[dict[str, Tensor]]:
    """Get features from a tokenized object with key filtering.

    Args:
        tokenized: Tokenized object
        window_batching: Whether to use window batching
        shard_dims: Shard dimensions to enable window batching.
        selected_keys: Selected keys to return. If None, the following keys are returned:
            ["atom_pad_mask", "atom_to_token", "pair_mask", "token_pad_mask"]

    Returns:
        dict[str, Tensor] | list[dict[str, Tensor]]: Features of list of features if shard_dims is not None.
    """
    max_atoms, max_tokens = get_num_atoms_tokens(tokenized)
    if shard_dims is None:
        max_atoms = None
        max_tokens = None
        max_seqs = 1
    else:
        ring_size = shard_dims[0]
        if window_batching:
            pad_to_multiple_atoms = math.lcm(ring_size, 32)
        else:
            pad_to_multiple_atoms = ring_size
        max_atoms = max_atoms + pad_to_multiple_atoms - max_atoms % pad_to_multiple_atoms
        max_tokens = max_tokens + ring_size - max_tokens % ring_size
        max_seqs = ring_size

    featurizer = BoltzFeaturizer()
    feats = featurizer.process(
        tokenized,
        training=False,
        augmentation=False,
        pair_mask_mode=PairMaskMode.NONE if window_batching else PairMaskMode.SEQUENCE_LOCAL_ATTENTION,
        max_atoms=max_atoms,
        max_tokens=max_tokens,
        max_seqs=max_seqs,
        pad_to_max_seqs=True,
        shard_dims=shard_dims,
    )

    if selected_keys is None:
        selected_keys = ["atom_pad_mask", "atom_to_token", "pair_mask", "token_pad_mask"]

    if shard_dims is None:
        return {k: v for k, v in feats.items() if k in selected_keys}
    else:
        return [{k: v for k, v in feats_shard.items() if k in selected_keys} for feats_shard in feats]


def get_to_keys(s_global: torch.Tensor, W: int = 32, H: int = 128) -> Callable[[torch.Tensor], torch.Tensor]:
    """Get to keys function for window-based attention."""
    B, N, D = s_global.shape
    # assume s_global.shape[1] has been padded to be a multiple of W
    assert N % W == 0, "s_global.shape[1] must be a multiple of W"
    K = N // W
    indexing_matrix = get_indexing_matrix(K, W, H, s_global.device).to(s_global.dtype)
    to_keys = partial(single_to_keys, indexing_matrix=indexing_matrix, W=W, H=H)
    return to_keys


def create_msa_module_init_params(use_large_model: bool = False) -> dict[str, Any]:
    """Create initialization parameters for MSAModule.

    Parameters
    ----------
    use_large_model : bool
        Whether to use large model parameters

    Returns
    -------
    dict[str, Any]
        MSAModule initialization parameters
    """

    # Get parameters from the whole set
    boltz1_params = create_boltz1_model_init_params(use_large_model=use_large_model, use_window_batching=True)

    # Extract MSA-specific parameters
    msa_args = boltz1_params["msa_args"].copy()

    # Calculate s_input_dim based on token_s
    token_s = boltz1_params["token_s"]
    msa_args["s_input_dim"] = (
        token_s + 2 * const.num_tokens + 1 + len(const.pocket_contact_info)
    )  # Input sequence dimension

    # Add token_z for compatibility
    msa_args["token_z"] = boltz1_params["token_z"]

    return msa_args


def create_msa_module_init_params_v2(use_large_model: bool = False) -> dict[str, Any]:
    """Create initialization parameters for Boltz-2 MSAModule (model.modules.trunkv2.MSAModule).

    Parameters
    ----------
    use_large_model : bool
        Whether to use large model parameters

    Returns
    -------
    dict[str, Any]
        MSAModule initialization parameters for Boltz-2 (token_s, msa_s, token_z, etc.)
    """
    boltz1_params = create_boltz1_model_init_params(use_large_model=use_large_model, use_window_batching=True)
    msa_args = boltz1_params["msa_args"].copy()
    # Boltz-2 MSAModule uses token_s (single representation dim), not s_input_dim
    msa_args["token_s"] = boltz1_params["token_s"]
    msa_args["token_z"] = boltz1_params["token_z"]
    msa_args["subsample_msa"] = False  # CP does not support MSA subsampling
    msa_args["num_subsampled_msa"] = 1024
    return msa_args


def create_pairformer_module_init_params(use_large_model: bool = False) -> dict[str, Any]:
    """Create initialization parameters for PairformerModule.

    Parameters
    ----------
    use_large_model : bool
        Whether to use large model parameters

    Returns
    -------
    dict[str, Any]
        PairformerModule initialization parameters
    """
    # Get parameters from the whole set
    boltz1_params = create_boltz1_model_init_params(use_large_model=use_large_model, use_window_batching=True)

    # Extract Pairformer-specific parameters
    pairformer_args = boltz1_params["pairformer_args"].copy()

    # Add required shared parameters
    pairformer_args["token_s"] = boltz1_params["token_s"]
    pairformer_args["token_z"] = boltz1_params["token_z"]

    return pairformer_args


def create_diffusion_module_init_params(
    use_large_model: bool = False, use_window_batching: bool = True
) -> dict[str, Any]:
    """Create initialization parameters for DiffusionModule.

    Parameters
    ----------
    use_large_model : bool
        Whether to use large model parameters
    use_window_batching : bool
        Whether to enable window batching

    Returns
    -------
    dict[str, Any]
        DiffusionModule initialization parameters
    """
    # Get parameters from the whole set
    boltz1_params = create_boltz1_model_init_params(
        use_large_model=use_large_model, use_window_batching=use_window_batching
    )

    # Extract Diffusion-specific parameters
    score_model_args = boltz1_params["score_model_args"].copy()

    # Add required shared parameters
    params = {
        "token_s": boltz1_params["token_s"],
        "token_z": boltz1_params["token_z"],
        "atom_s": boltz1_params["atom_s"],
        "atom_z": boltz1_params["atom_z"],
        "atoms_per_window_queries": boltz1_params["atoms_per_window_queries"],
        "atoms_per_window_keys": boltz1_params["atoms_per_window_keys"],
        "atom_feature_dim": boltz1_params["atom_feature_dim"],
        **score_model_args,
    }

    return params


def create_atom_diffusion_init_params(
    use_large_model: bool = False, use_window_batching: bool = True
) -> dict[str, Any]:
    """Create initialization parameters for AtomDiffusion.

    Parameters
    ----------
    use_large_model : bool
        Whether to use large model parameters
    use_window_batching : bool
        Whether to enable window batching

    Returns
    -------
    dict[str, Any]
        AtomDiffusion initialization parameters
    """
    # Get parameters from the whole set
    boltz1_params = create_boltz1_model_init_params(
        use_large_model=use_large_model, use_window_batching=use_window_batching
    )

    # Extract AtomDiffusion-specific parameters
    params = {
        "score_model_args": boltz1_params["score_model_args"],
        **boltz1_params["diffusion_process_args"],
        "compile_score": False,  # couldn't be set in the whole-model parameter but enforced here
        "accumulate_token_repr": False,  # couldn't be set in the whole-model parameter but enforced here
    }

    return params


@dataclass(frozen=False)
class TrainingArgs:
    recycling_steps: int
    sampling_steps: int | None
    diffusion_multiplicity: int
    diffusion_samples: int
    confidence_loss_weight: float
    diffusion_loss_weight: float
    distogram_loss_weight: float


def create_boltz1_model_init_params(
    use_large_model: bool = False, use_window_batching: bool = True, activation_checkpointing: bool = True
) -> dict[str, Any]:
    """Create initialization parameters for Boltz1 model.

    Parameters
    ----------
    use_large_model : bool
        Whether to use large model parameters
    use_window_batching : bool
        Whether to enable window batching
    activation_checkpointing : bool
        Whether to use activation checkpointing

    Returns
    -------
    dict[str, Any]
        Boltz1 model initialization parameters
    """
    if BoltzDiffusionParams is not None:
        diffusion_params = asdict(BoltzDiffusionParams())
    else:
        diffusion_params = {
            "gamma_0": 0.605,
            "gamma_min": 1.107,
            "noise_scale": 0.901,
            "rho": 8,
            "step_scale": 1.638,
            "sigma_min": 0.0004,
            "sigma_max": 160.0,
            "sigma_data": 16.0,
            "P_mean": -1.2,
            "P_std": 1.5,
            "coordinate_augmentation": True,
            "alignment_reverse_diff": True,
            "synchronize_sigmas": True,
            "use_inference_model_cache": True,
        }

    # can't set the following:
    # "compile_score": False,
    # "accumulate_token_repr": False,
    # because model.py will set them explicitly
    # then we would have duplicated keys in the input kwargs
    # but the default settings of these two should work for testing
    diffusion_params.update(
        {
            "coordinate_augmentation": False,  # Turn off for deterministic testing
            "alignment_reverse_diff": True,
            "synchronize_sigmas": False,
            "use_inference_model_cache": False,
            "num_sampling_steps": None,  # only relevant to sample() test but disabled otherwise to prevent accidentally enabling irrelevant code
        }
    )

    if use_large_model:
        atom_s = 128
        atom_z = 64
        token_s = 384
        token_z = 128
        # Distogram module parameters
        num_bins = 64
        # InputEmbedder parameters
        atom_encoder_depth = 2
        atom_encoder_heads = 4
        # MSAModule parameters
        msa_s = 64
        msa_blocks = 4
        pairwise_head_width = 32
        pairwise_num_heads = 4
        # PairformerModule parameters
        num_pairformer_blocks = 4
        num_pairformer_heads = 16
        # AtomDiffusion parameters
        sigma_data = 16
        dim_fourier = 256
        atom_encoder_depth_diffusion = 3
        atom_encoder_heads_diffusion = 4
        token_transformer_depth_diffusion = 24
        token_transformer_heads_diffusion = 16
        atom_decoder_depth_diffusion = 3
        atom_decoder_heads_diffusion = 4
        conditioning_transition_layers_diffusion = 2
        activation_checkpointing_diffusion = activation_checkpointing
        offload_to_cpu_diffusion = False
        # Training parameters
        recycling_steps = 3
        diffusion_multiplicity = 16
        diffusion_samples = 1  # not used in training step
        confidence_loss_weight = 3e-3
        diffusion_loss_weight = 4.0
        distogram_loss_weight = 3e-2
    else:
        atom_s = 4
        atom_z = 2
        token_s = 4
        token_z = 12
        # Distogram module parameters
        num_bins = 4
        # InputEmbedder parameters
        atom_encoder_depth = 1
        atom_encoder_heads = 2
        # MSAModule parameters
        msa_s = 4
        msa_blocks = 1
        pairwise_head_width = 4
        pairwise_num_heads = 2
        # PairformerModule parameters
        num_pairformer_blocks = 1
        num_pairformer_heads = 4
        # AtomDiffusion parameters
        sigma_data = 16
        dim_fourier = 16
        atom_encoder_depth_diffusion = 1
        atom_encoder_heads_diffusion = 2
        token_transformer_depth_diffusion = 2
        token_transformer_heads_diffusion = 2
        atom_decoder_depth_diffusion = 1
        atom_decoder_heads_diffusion = 2
        conditioning_transition_layers_diffusion = 1
        activation_checkpointing_diffusion = activation_checkpointing
        offload_to_cpu_diffusion = False
        # Training parameters
        recycling_steps = 0
        diffusion_multiplicity = 1
        diffusion_samples = 1  # not used in training step
        confidence_loss_weight = 3e-3
        diffusion_loss_weight = 4.0
        distogram_loss_weight = 3e-2

    # Diffusion loss parameters
    diffusion_loss_args = {
        "add_smooth_lddt_loss": True,
        "nucleotide_loss_weight": 5.0,
        "ligand_loss_weight": 10.0,
    }

    params = {
        "atom_s": atom_s,
        "atom_z": atom_z,
        "token_s": token_s,
        "token_z": token_z,
        "num_bins": num_bins,
        "use_window_batching": use_window_batching,
        "embedder_args": {
            "atom_encoder_depth": atom_encoder_depth,
            "atom_encoder_heads": atom_encoder_heads,
            "activation_checkpointing": activation_checkpointing,
            "activation_checkpointing_pair_repr": activation_checkpointing,
        },
        "msa_args": {
            "msa_s": msa_s,
            "msa_blocks": msa_blocks,
            "msa_dropout": 0.0,
            "z_dropout": 0.0,
            "pairwise_head_width": pairwise_head_width,
            "pairwise_num_heads": pairwise_num_heads,
            "use_paired_feature": True,
            "activation_checkpointing": activation_checkpointing,
            "activation_checkpointing_pair_repr": activation_checkpointing,
            "offload_to_cpu": False,
        },
        "pairformer_args": {
            "num_blocks": num_pairformer_blocks,
            "num_heads": num_pairformer_heads,
            "pairwise_head_width": pairwise_head_width,
            "pairwise_num_heads": pairwise_num_heads,
            "dropout": 0.0,
            "no_update_s": False,
            "no_update_z": False,
            "activation_checkpointing": activation_checkpointing,
            "activation_checkpointing_pair_repr": activation_checkpointing,
            "offload_to_cpu": False,
        },
        "atom_feature_dim": 389,  # hardcoded hidden dim from featurizer stacking multiple features
        "atoms_per_window_queries": 32 if use_window_batching else None,
        "atoms_per_window_keys": 128 if use_window_batching else None,
        "no_msa": False,
        "no_atom_encoder": False,
        "min_dist": 2.0,
        "max_dist": 22.0,
        "do_activation_chunking": False,
        "score_model_args": {
            "sigma_data": sigma_data,
            "dim_fourier": dim_fourier,
            "atom_encoder_depth": atom_encoder_depth_diffusion,
            "atom_encoder_heads": atom_encoder_heads_diffusion,
            "token_transformer_depth": token_transformer_depth_diffusion,
            "token_transformer_heads": token_transformer_heads_diffusion,
            "atom_decoder_depth": atom_decoder_depth_diffusion,
            "atom_decoder_heads": atom_decoder_heads_diffusion,
            "conditioning_transition_layers": conditioning_transition_layers_diffusion,
            "activation_checkpointing": activation_checkpointing_diffusion,
            "activation_checkpointing_pair_repr": activation_checkpointing_diffusion,
            "offload_to_cpu": offload_to_cpu_diffusion,
        },
        "diffusion_process_args": diffusion_params,
        # TODO: support validation, confidence and steering args
        "training_args": TrainingArgs(
            recycling_steps=recycling_steps,
            sampling_steps=None,  # not used in training step
            diffusion_multiplicity=diffusion_multiplicity,
            diffusion_samples=diffusion_samples,  # not used in training step
            confidence_loss_weight=confidence_loss_weight,
            diffusion_loss_weight=diffusion_loss_weight,
            distogram_loss_weight=distogram_loss_weight,
        ),
        "validation_args": {},
        "diffusion_loss_args": diffusion_loss_args,
        "confidence_model_args": {},
        "steering_args": {},
    }
    params["score_model_args"].update(
        {
            "atom_s": atom_s,
            "atom_z": atom_z,
            "token_s": token_s,
            "token_z": token_z,
            "atom_feature_dim": 389,  # hardcoded hidden dim from featurizer stacking multiple features
            "atoms_per_window_queries": 32 if use_window_batching else None,
            "atoms_per_window_keys": 128 if use_window_batching else None,
        }
    )

    return params


class DictNamespace:
    """Picklable namespace with both attribute access and dict-style ``.get()``."""

    def __init__(self, **kwargs: Any) -> None:
        self.__dict__.update(kwargs)

    def get(self, key: str, default: Any = None) -> Any:
        return self.__dict__.get(key, default)


def create_boltz2_model_init_params(
    use_large_model: bool = False,
    activation_checkpointing: bool = False,
) -> dict[str, Any]:
    """Create initialization parameters for Boltz2 model.

    Parameters
    ----------
    use_large_model : bool
        Whether to use large model parameters (closer to production config).
    activation_checkpointing : bool
        Whether to use activation checkpointing.

    Returns
    -------
    dict[str, Any]
        Boltz2 model initialization parameters.
    """
    diffusion_process_args = {
        "sigma_min": 0.0004,
        "sigma_max": 160.0,
        "sigma_data": 16.0,
        "rho": 7,
        "P_mean": -1.2,
        "P_std": 1.5,
        "gamma_0": 0.8,
        "gamma_min": 1.0,
        "noise_scale": 1.0,
        "step_scale": 1.0,
        "coordinate_augmentation": False,
        "alignment_reverse_diff": True,
        "synchronize_sigmas": False,
    }

    if use_large_model:
        atom_s = 128
        atom_z = 16
        token_s = 384
        token_z = 128
        num_bins = 64
        atom_encoder_depth = 3
        atom_encoder_heads = 4
        msa_s = 64
        msa_blocks = 4
        pairwise_head_width = 32
        pairwise_num_heads = 4
        num_pairformer_blocks = 4
        num_pairformer_heads = 16
        sigma_data = 16
        dim_fourier = 256
        atom_encoder_depth_diffusion = 3
        atom_encoder_heads_diffusion = 4
        token_transformer_depth_diffusion = 24
        token_transformer_heads_diffusion = 16
        atom_decoder_depth_diffusion = 3
        atom_decoder_heads_diffusion = 4
        conditioning_transition_layers_diffusion = 2
    else:
        atom_s = 4
        atom_z = 2
        token_s = 4
        token_z = 12
        num_bins = 4
        atom_encoder_depth = 1
        atom_encoder_heads = 2
        msa_s = 4
        msa_blocks = 1
        pairwise_head_width = 4
        pairwise_num_heads = 2
        num_pairformer_blocks = 1
        num_pairformer_heads = 4
        sigma_data = 16
        dim_fourier = 16
        atom_encoder_depth_diffusion = 1
        atom_encoder_heads_diffusion = 2
        token_transformer_depth_diffusion = 2
        token_transformer_heads_diffusion = 2
        atom_decoder_depth_diffusion = 1
        atom_decoder_heads_diffusion = 2
        conditioning_transition_layers_diffusion = 1

    training_args = DictNamespace(
        recycling_steps=3 if use_large_model else 0,
        sampling_steps=None,
        sampling_steps_random=None,
        diffusion_multiplicity=2,
        diffusion_samples=1,
        confidence_loss_weight=0.0,
        diffusion_loss_weight=4.0,
        distogram_loss_weight=3e-2,
        bfactor_loss_weight=0.0,
        symmetry_correction=False,
        adam_beta_1=0.9,
        adam_beta_2=0.95,
        adam_eps=1e-8,
        base_lr=1e-3,
        max_lr=1e-3,
        lr_scheduler="af3",
        lr_warmup_no_steps=10,
        lr_start_decay_after_n_steps=100,
        lr_decay_every_n_steps=50000,
        lr_decay_factor=0.95,
        weight_decay=0.0,
    )

    validation_args = DictNamespace(
        recycling_steps=0,
        sampling_steps=2,
        diffusion_samples=1,
        symmetry_correction=False,
        run_confidence_sequentially=False,
    )

    params: dict[str, Any] = {
        "atom_s": atom_s,
        "atom_z": atom_z,
        "token_s": token_s,
        "token_z": token_z,
        "num_bins": num_bins,
        "atom_feature_dim": 388,
        "atoms_per_window_queries": 32,
        "atoms_per_window_keys": 128,
        "embedder_args": {
            "atom_encoder_depth": atom_encoder_depth,
            "atom_encoder_heads": atom_encoder_heads,
            "activation_checkpointing": activation_checkpointing,
            "add_mol_type_feat": True,
            "add_method_conditioning": True,
            "add_modified_flag": True,
            "add_cyclic_flag": True,
        },
        "msa_args": {
            "msa_s": msa_s,
            "msa_blocks": msa_blocks,
            "msa_dropout": 0.0,
            "z_dropout": 0.0,
            "pairwise_head_width": pairwise_head_width,
            "pairwise_num_heads": pairwise_num_heads,
            "use_paired_feature": True,
            "activation_checkpointing": activation_checkpointing,
        },
        "pairformer_args": {
            "num_blocks": num_pairformer_blocks,
            "num_heads": num_pairformer_heads,
            "dropout": 0.0,
            "v2": True,
            "post_layer_norm": False,
            "activation_checkpointing": activation_checkpointing,
        },
        "score_model_args": {
            "sigma_data": sigma_data,
            "dim_fourier": dim_fourier,
            "atom_encoder_depth": atom_encoder_depth_diffusion,
            "atom_encoder_heads": atom_encoder_heads_diffusion,
            "token_transformer_depth": token_transformer_depth_diffusion,
            "token_transformer_heads": token_transformer_heads_diffusion,
            "atom_decoder_depth": atom_decoder_depth_diffusion,
            "atom_decoder_heads": atom_decoder_heads_diffusion,
            "conditioning_transition_layers": conditioning_transition_layers_diffusion,
            "activation_checkpointing": activation_checkpointing,
            "transformer_post_ln": False,
        },
        "diffusion_process_args": diffusion_process_args,
        "diffusion_loss_args": {
            "add_smooth_lddt_loss": True,
            "nucleotide_loss_weight": 5.0,
            "ligand_loss_weight": 10.0,
        },
        "training_args": training_args,
        "validation_args": validation_args,
        "confidence_prediction": False,
        "affinity_prediction": False,
        "structure_prediction_training": True,
        "validate_structure": False,
        "use_templates": False,
        "predict_bfactor": True,
        "bond_type_feature": True,
        "steering_args": None,
        "confidence_model_args": {
            "num_dist_bins": 64,
            "max_dist": 22,
            "add_s_to_z_prod": True,
            "add_s_input_to_s": True,
            "add_z_input_to_z": True,
            "conditioning_cutoff_min": 4.0,
            "conditioning_cutoff_max": 20.0,
            "confidence_args": {
                "num_plddt_bins": 50,
                "num_pde_bins": 64,
                "num_pae_bins": 64,
            },
        },
        "ema": False,
    }
    return params


def random_features(
    size_batch: int,
    n_tokens: int,
    n_atoms: int,
    n_msa: int,
    atom_counts_per_token_range: tuple[int, int],
    device: torch.device,
    float_value_range: tuple[float, float],
    selected_keys: Optional[list[str]] = None,
    num_disto_bins: int = 64,
    rng: Optional[torch.Generator] = None,
) -> dict[str, torch.Tensor]:
    """Generate random feature tensors matching the shapes and dtypes of selected_keys features.

    NOTE: This function uses all-valid masks (mask with all True) for token, atoms and MSA
    NOTE: the returned tensors dtype for floating point features is torch.float64
    NOTE: constraints and symmetry features are not supported yet

    Parameters
    ----------
    size_batch : int
        Batch size
    n_tokens : int
        Number of tokens
    n_atoms : int
        Number of atoms
    n_msa : int
        Number of MSA sequences
    atom_counts_per_token_range : tuple[int, int]
        Range for number of atoms per token (min, max)
    device : torch.device
        Device to create tensors on
    float_value_range : tuple[float, float]
        Range for float values (min, max)
    selected_keys : Optional[list[str]]
        If provided, only return features for these keys. If None, return all features.
    num_disto_bins : int
        Number of bins for distogram target. Default is 64.
    rng : Optional[torch.Generator]
        Optional random generator for deterministic sampling without modifying global RNG state.

    Returns
    -------
    dict[str, torch.Tensor]
        Dictionary of randomly generated feature tensors
    """
    features = {}

    # Generate atom_counts_per_token first to ensure proper atom_to_token mapping
    min_atoms, max_atoms = atom_counts_per_token_range
    if min_atoms < 1:
        raise ValueError(f"min_atoms must be >= 1, got {min_atoms}")
    if max_atoms < min_atoms:
        raise ValueError(f"max_atoms ({max_atoms}) must be >= min_atoms ({min_atoms})")

    # For now, to avoid collating different samples of different atom counts per token, we
    # generate the same atom counts per token for all samples.
    atom_counts_per_token = (
        torch.randint(min_atoms, max_atoms + 1, (n_tokens,), dtype=torch.int64, device=device, generator=rng)
        .unsqueeze(0)
        .repeat_interleave(size_batch, 0)
    )
    # Ensure rewrite-path coverage for confidence.compute_frame_pred:
    # when possible, force the last 3 tokens to be a contiguous NONPOLYMER segment
    # (each has one atom, so the chain has >=3 atoms in total).
    force_nonpolymer_tail = min_atoms <= 1 and n_tokens >= 4 and n_atoms >= n_tokens
    if force_nonpolymer_tail:
        atom_counts_per_token[:, -3:] = 1

    features["atom_counts_per_token"] = atom_counts_per_token

    # Ensure total atoms match n_atoms by adjusting one anchor token.
    # If we force a nonpolymer tail, keep the last 3 tokens fixed at 1 atom each.
    anchor_idx = n_tokens - 4 if force_nonpolymer_tail else n_tokens - 1
    current_total_except_anchor = atom_counts_per_token.sum(dim=1) - atom_counts_per_token[:, anchor_idx]
    if (current_total_except_anchor >= n_atoms).any():
        raise ValueError(
            f"Total atoms {current_total_except_anchor} excluding anchor token {anchor_idx} exceeds n_atoms {n_atoms}"
        )
    atom_counts_per_token[:, anchor_idx] = n_atoms - current_total_except_anchor

    # Create atom_to_token one-hot mapping based on atom_counts_per_token
    atom_to_token_ccol_ids = torch.zeros((size_batch, n_tokens + 1), dtype=torch.int64, device=device)
    atom_to_token_ccol_ids[:, 1:] = atom_counts_per_token.cumsum(dim=1)
    atom_to_token_row_ids = (
        torch.arange(n_atoms, dtype=torch.int64, device=device).unsqueeze(0).repeat_interleave(size_batch, 0)
    )
    atom_to_token_values = torch.ones_like(atom_to_token_row_ids, dtype=torch.int64, device=device)
    atom_to_token_csc = torch.sparse_csc_tensor(
        atom_to_token_ccol_ids,
        atom_to_token_row_ids,
        atom_to_token_values,
        size=(size_batch, n_atoms, n_tokens),
        dtype=torch.int64,
        device=device,
    )
    atom_to_token = atom_to_token_csc.to_dense()

    # TODO: support heterogeneous atom counts per token across samples in a batch
    assert (atom_to_token[0] == atom_to_token).all(), "atom_to_token is not identical across samples in a batch"

    # Create token_to_rep_atom one-hot mapping (each token picks one representative atom randomly)
    # token_to_rep_atom has shape (size_batch, n_tokens, n_atoms)
    # For token i, we pick a random atom from the atoms it owns: [cumsum[i], cumsum[i+1])
    token_atom_start = atom_to_token_ccol_ids[:, :-1]  # (size_batch, n_tokens) - start index of atoms for each token

    # For each token, pick a random offset within its atom range [0, atom_count)
    # Since atom_counts_per_token is identical across batch, we can vectorize
    # Note: min_atoms >= 1 is enforced above, so max_count is always > 0
    max_count = atom_counts_per_token.max().item()
    random_offsets = torch.randint(0, max_count, (size_batch, n_tokens), device=device, generator=rng)
    # Clamp offsets to be within [0, count-1] for each token
    random_offsets = random_offsets % atom_counts_per_token

    # Representative atom index = start + offset
    rep_atom_indices = token_atom_start + random_offsets  # (size_batch, n_tokens)

    # Create one-hot token_to_rep_atom using one_hot
    token_to_rep_atom = torch.nn.functional.one_hot(rep_atom_indices, num_classes=n_atoms)

    features["token_to_rep_atom"] = token_to_rep_atom

    # Create r_set_to_rep_atom: randomly select subset of tokens as R-set elements
    # For each R-set element, assign a representative atom from within that token's atoms
    #
    # Technical notes:
    # - N_R is typically N_resolved_polymer_tokens; we simulate with a random subset of tokens
    # - r_set_to_rep_atom is stored as one-hot tensor [B, N_R, N_atoms]
    # - The featurizer preserves this format and shards it as diagonal blocks aligned with
    #   token sharding. This enables local einsum for atom-to-R-set coordinate mapping
    #   without cross-shard communication in the distributed plddt_loss.
    # - R-set token indices are identical across batch (intentional for consistent sharding)
    #
    n_r = max(
        1,
        n_tokens - torch.randint(0, max(1, n_tokens // 4), (1,), device=device, generator=rng).item(),
    )  # Random N_R <= n_tokens

    # Randomly select which tokens are in the R-set (sorted for consistency across batch)
    r_set_token_indices = torch.randperm(n_tokens, device=device, generator=rng)[:n_r].sort().values

    # Vectorized: get atom start indices and counts for R-set tokens
    # atom_to_token_ccol_ids shape: [B, n_tokens + 1] (cumulative column indices)
    r_set_atom_start = atom_to_token_ccol_ids[:, r_set_token_indices]  # [B, N_R]
    r_set_atom_end = atom_to_token_ccol_ids[:, r_set_token_indices + 1]  # [B, N_R]
    r_set_atom_counts = r_set_atom_end - r_set_atom_start  # [B, N_R]

    # Generate random offsets within each token's atom range
    max_atoms_in_token = r_set_atom_counts.max().item()
    if max_atoms_in_token > 0:
        r_offsets = torch.randint(
            0, max(1, max_atoms_in_token), (size_batch, n_r), device=device, dtype=torch.int64, generator=rng
        )
        # Clamp to valid range using modulo (handles varying atom counts per token)
        r_offsets = r_offsets % torch.clamp(r_set_atom_counts, min=1)
    else:
        r_offsets = torch.zeros((size_batch, n_r), device=device, dtype=torch.int64)

    # Compute representative atom indices
    r_set_rep_atom_indices = r_set_atom_start + r_offsets  # [B, N_R]

    # Create one-hot tensor using F.one_hot
    r_set_to_rep_atom = torch.nn.functional.one_hot(r_set_rep_atom_indices, num_classes=n_atoms).to(
        dtype=torch.float64
    )  # [B, N_R, N_atoms]

    features["r_set_to_rep_atom"] = r_set_to_rep_atom

    # Extract float range values
    min_val, max_val = float_value_range

    # Core features for InputEmbedder
    features["atom_pad_mask"] = torch.ones((size_batch, n_atoms), dtype=torch.float64, device=device)
    features["atom_to_token"] = atom_to_token
    features["pair_mask"] = (
        get_pair_mask(n_atoms).unsqueeze(0).repeat(size_batch, 1, 1).to(dtype=torch.float64, device=device)
    )
    features["token_pad_mask"] = torch.ones((size_batch, n_tokens), dtype=torch.float64, device=device)
    features["ref_pos"] = torch.empty(size_batch, n_atoms, 3, dtype=torch.float64, device=device).uniform_(
        min_val, max_val, generator=rng
    )
    features["ref_charge"] = torch.randint(-3, 4, (size_batch, n_atoms), dtype=torch.int8, device=device, generator=rng)
    features["ref_element"] = torch.randint(
        0,
        const.num_elements,
        (size_batch, n_atoms, const.num_elements),
        dtype=torch.int64,
        device=device,
        generator=rng,
    )
    features["ref_atom_name_chars"] = torch.randint(
        0, 64, (size_batch, n_atoms, 4, 64), dtype=torch.int64, device=device, generator=rng
    )
    features["ref_space_uid"] = torch.randint(
        0, 100, (size_batch, n_atoms), dtype=torch.int64, device=device, generator=rng
    )
    features["res_type"] = torch.randint(
        0, const.num_tokens, (size_batch, n_tokens, const.num_tokens), dtype=torch.int64, device=device, generator=rng
    )
    features["profile"] = torch.empty(
        size_batch, n_tokens, const.num_tokens, dtype=torch.float64, device=device
    ).uniform_(min_val, max_val, generator=rng)
    features["deletion_mean"] = torch.empty(size_batch, n_tokens, dtype=torch.float64, device=device).uniform_(
        min_val, max_val, generator=rng
    )
    features["pocket_feature"] = torch.randint(
        0, 2, (size_batch, n_tokens, 4), dtype=torch.int64, device=device, generator=rng
    )

    # Additional features for AtomDiffusion
    features["atom_resolved_mask"] = torch.ones((size_batch, n_atoms), dtype=torch.float64, device=device)

    # Additional features for MSA module
    features["msa"] = torch.randint(
        0,
        const.num_tokens,
        (size_batch, n_msa, n_tokens, const.num_tokens),
        dtype=torch.int64,
        device=device,
        generator=rng,
    )
    features["has_deletion"] = torch.randint(
        0, 2, (size_batch, n_msa, n_tokens), dtype=torch.bool, device=device, generator=rng
    )
    features["deletion_value"] = torch.empty(size_batch, n_msa, n_tokens, dtype=torch.float64, device=device).uniform_(
        min_val, max_val, generator=rng
    )
    features["msa_paired"] = torch.empty(size_batch, n_msa, n_tokens, dtype=torch.float64, device=device).uniform_(
        min_val, max_val, generator=rng
    )
    features["msa_mask"] = torch.ones((size_batch, n_msa, n_tokens), dtype=torch.int64, device=device)

    # Additional features for Boltz1
    # NOTE: token_bonds typically is binary but in the model workflow will go thru linear projection so it's float
    # anyway
    features["token_bonds"] = torch.empty(
        size_batch, n_tokens, n_tokens, 1, dtype=torch.float64, device=device
    ).uniform_(min_val, max_val, generator=rng)
    features["type_bonds"] = torch.randint(
        0, len(const.bond_types), (size_batch, n_tokens, n_tokens), dtype=torch.long, device=device, generator=rng
    )
    features["token_pair_pad_mask"] = features["token_pad_mask"][:, :, None] * features["token_pad_mask"][:, None, :]

    # Additional features for RelativePositionEncoder
    features["residue_index"] = torch.randint(
        0, 1000, (size_batch, n_tokens), dtype=torch.int64, device=device, generator=rng
    )
    features["entity_id"] = torch.randint(
        0, 10, (size_batch, n_tokens), dtype=torch.int64, device=device, generator=rng
    )
    features["token_index"] = (
        torch.arange(n_tokens, dtype=torch.int64, device=device).unsqueeze(0).expand(size_batch, -1).contiguous()
    )
    features["cyclic_period"] = torch.randint(
        0, 100, (size_batch, n_tokens), dtype=torch.int32, device=device, generator=rng
    )

    # Additional features for AtomDiffusion
    features["coords"] = torch.randn(size_batch, 1, n_atoms, 3, dtype=torch.float64, device=device, generator=rng)
    disto_target = torch.randint(0, num_disto_bins, (size_batch, n_tokens, n_tokens), device=device, generator=rng)
    features["disto_target"] = torch.nn.functional.one_hot(disto_target, num_classes=num_disto_bins).to(
        dtype=torch.float64
    )
    features["token_disto_mask"] = torch.randint(
        0, 2, (size_batch, n_tokens), dtype=torch.float64, device=device, generator=rng
    )

    # Additional features for confidence module
    features["atom_resolved_mask"] = torch.randint(
        0, 2, (size_batch, n_atoms), dtype=torch.float64, device=device, generator=rng
    )

    # Ensure chain-type consistency:
    # - non-polymer tokens have exactly 1 atom
    # - polymer tokens have > 3 atoms and are PROTEIN/DNA/RNA
    atom_counts_per_token = features["atom_to_token"].sum(dim=1).to(torch.int64)
    nonpolymer_flags = atom_counts_per_token == 1
    polymer_type_ids = torch.tensor(
        [
            const.chain_type_ids["PROTEIN"],
            const.chain_type_ids["DNA"],
            const.chain_type_ids["RNA"],
        ],
        dtype=torch.int64,
        device=device,
    )
    polymer_type_idx = torch.randint(
        0, polymer_type_ids.numel(), (size_batch, n_tokens), dtype=torch.int64, device=device, generator=rng
    )
    features["mol_type"] = polymer_type_ids[polymer_type_idx]
    features["mol_type"][nonpolymer_flags] = const.chain_type_ids["NONPOLYMER"]
    if force_nonpolymer_tail:
        features["mol_type"][:, -3:] = const.chain_type_ids["NONPOLYMER"]
        if n_tokens > 3:
            features["mol_type"][:, -4] = const.chain_type_ids["PROTEIN"]

    # Assign asym_id by contiguous mol_type segments per batch.
    # E.g. [0, 0, 0, 1, 1, 0, 2, 2] -> [0, 0, 0, 1, 1, 2, 3, 3]
    asym_id = torch.zeros((size_batch, n_tokens), device=device, dtype=torch.int64)
    if n_tokens > 1:
        mol_changes = features["mol_type"][:, 1:] != features["mol_type"][:, :-1]
        asym_id[:, 1:] = mol_changes.to(torch.int64).cumsum(dim=1)
    features["asym_id"] = asym_id

    # Fuse consecutive asym_id values into shared sym_id buckets per batch.
    # E.g. asym_id unique [0,1,2,3,4] -> sym_id [0,0,1,1,2]
    sym_id = torch.empty_like(asym_id)
    for batch_idx in range(size_batch):
        unique_asym = torch.unique(asym_id[batch_idx])
        sym_map = torch.arange(unique_asym.numel(), device=device) // 2
        for asym_value, sym_value in zip(unique_asym.tolist(), sym_map.tolist()):
            sym_id[batch_idx][asym_id[batch_idx] == asym_value] = sym_value
    features["sym_id"] = sym_id

    if not (features["mol_type"] == const.chain_type_ids["NONPOLYMER"]).any():
        warnings.warn(
            "No non-polymer token is created in random features generation.",
            stacklevel=2,
        )

    # Build frames_idx from atom_to_token using token start offsets.
    # For each token, sample 3 distinct local offsets within its atom range,
    # then map to global indices.
    max_atom_counts_per_token = int(atom_counts_per_token.max().item())
    offset_ids = torch.arange(max_atom_counts_per_token, device=device).view(1, 1, -1)
    rand_scores = torch.rand(size_batch, n_tokens, max_atom_counts_per_token, device=device, generator=rng)
    invalid_offsets = offset_ids >= atom_counts_per_token.unsqueeze(-1)
    rand_scores = rand_scores.masked_fill(invalid_offsets, float("inf"))
    offsets = torch.argsort(rand_scores, dim=-1)[..., :3]
    frames_idx = token_atom_start.unsqueeze(-1) + offsets

    # Match featurizer behavior for tokens with fewer than 3 atoms:
    # use the first atom of the token for all three frame slots.
    small_token_mask = atom_counts_per_token.unsqueeze(-1) < 3
    if small_token_mask.any():
        token_start_triplet = token_atom_start.unsqueeze(-1).expand(-1, -1, 3)
        frames_idx = torch.where(small_token_mask, token_start_triplet, frames_idx)
    features["frames_idx"] = frames_idx

    # Derive frame_resolved_mask: True when all 3 frame atoms are resolved.
    # Matches the featurizer logic in compute_frames_nonpolymer which sets
    # resolved_frame_data[t] = resolved_mask[frames[t]].all().
    batch_expand = torch.arange(size_batch, device=device).view(-1, 1, 1).expand_as(frames_idx)
    frame_atoms_resolved = features["atom_resolved_mask"][batch_expand, frames_idx]  # (B, T, 3)
    features["frame_resolved_mask"] = frame_atoms_resolved.prod(dim=-1)  # (B, T)

    # is_nonpolymer_with_frame: True for non-polymer tokens in chains with >= 3 atoms
    atoms_per_chain = torch.zeros_like(asym_id, dtype=torch.int64)
    atoms_per_chain.scatter_add_(1, asym_id, atom_counts_per_token)
    chain_total_atoms = atoms_per_chain.gather(1, asym_id)
    features["is_nonpolymer_with_frame"] = nonpolymer_flags & (chain_total_atoms >= 3)

    # Boltz-2 specific features (unconditionally generated, opt-in via selected_keys)
    num_cc_types = len(const.contact_conditioning_info)
    rand_type = torch.rand(size_batch, n_tokens, n_tokens, device=device, generator=rng)
    cc = torch.zeros(size_batch, n_tokens, n_tokens, num_cc_types, dtype=torch.float64, device=device)
    cc[:, :, :, 0] = (rand_type < 0.3).to(torch.float64)
    cc[:, :, :, 1] = ((rand_type >= 0.3) & (rand_type < 0.5)).to(torch.float64)
    if num_cc_types > 4:  # noqa: PLR2004
        cc[:, :, :, 4] = (rand_type >= 0.5).to(torch.float64)
    features["contact_conditioning"] = cc
    features["contact_threshold"] = (
        torch.rand(size_batch, n_tokens, n_tokens, device=device, generator=rng) * 22.0
    ).to(torch.float64)
    features["method_feature"] = torch.randint(
        0, const.num_method_types, (size_batch, n_tokens), dtype=torch.int64, device=device, generator=rng
    )
    features["modified"] = torch.randint(0, 2, (size_batch, n_tokens), dtype=torch.int64, device=device, generator=rng)
    features["bfactor"] = torch.empty(size_batch, n_atoms, dtype=torch.float64, device=device).uniform_(
        min_val, max_val, generator=rng
    )
    features["plddt"] = torch.empty(size_batch, n_atoms, dtype=torch.float64, device=device).uniform_(
        0.0, 1.0, generator=rng
    )

    # Return only selected features if specified
    if selected_keys is not None:
        return {k: v for k, v in features.items() if k in selected_keys}
    else:
        return features


def get_features_shardable(
    tokenized: Tokenized,
    pair_mask_mode: PairMaskMode,
    return_shards: bool,
    shard_dims: tuple[int, int],
    selected_keys: Optional[list[str]] = None,
    **kwargs_feat_process: dict[str, Any],
) -> dict[str, torch.Tensor] | list[dict[str, torch.Tensor]]:
    """Get features from a tokenized object with sharding support.

    Args:
        tokenized: Tokenized object
        pair_mask_mode: Pair mask mode to use
        return_shards: Whether to return shards
        shard_dims: Shard dimensions to enable sharding
        selected_keys: Selected keys to return. If None, the following keys are returned:
            ["atom_pad_mask", "atom_to_token", "pair_mask", "token_pad_mask"]

    Returns:
        dict[str, torch.Tensor] | list[dict[str, torch.Tensor]]: Features or list of features if return_shards is True.
    """
    if return_shards and shard_dims is None:
        raise ValueError("shard_dims must be provided if return_shards is True")

    num_atoms, num_tokens = get_num_atoms_tokens(tokenized)
    # always pad to the next multiple of the ring size
    ring_size = shard_dims[0]
    # max_atoms would need to be padded to the least common multiple of ring_size and atoms_per_window_queries=32
    atom_counts_lcm = math.lcm(ring_size, 32)
    max_atoms = ((num_atoms + atom_counts_lcm - 1) // atom_counts_lcm) * atom_counts_lcm
    max_tokens = ((num_tokens + ring_size - 1) // ring_size) * ring_size
    max_seqs = ring_size

    featurizer = BoltzFeaturizer()
    feats = featurizer.process(
        tokenized,
        training=False,
        augmentation=False,
        pair_mask_mode=pair_mask_mode,
        max_atoms=max_atoms,
        max_tokens=max_tokens,
        max_seqs=max_seqs,
        pad_to_max_seqs=True,
        shard_dims=shard_dims if return_shards else None,
        **kwargs_feat_process,
    )

    if selected_keys is None:
        selected_keys = ["atom_pad_mask", "atom_to_token", "pair_mask", "token_pad_mask"]

    if return_shards:
        return [{k: v for k, v in feats_shard.items() if k in selected_keys} for feats_shard in feats]
    else:
        return {k: v for k, v in feats.items() if k in selected_keys}


def get_feature_placements(
    token_keys: Optional[set[str]] = None,
    msa_keys: Optional[set[str]] = None,
    atom_keys: Optional[set[str]] = None,
    model_io_keys: Optional[set[str]] = None,
    model_io_fp32_keys: Optional[set[str]] = None,
):
    """Get comprehensive feature placement definitions for distributed testing.

    Args:
        token_keys: Subset of token feature keys to include. If None, include all token features.
        msa_keys: Subset of MSA feature keys to include. If None, include all MSA features.
        atom_keys: Subset of atom feature keys to include. If None, include all atom features.
        model_io_keys: Subset of model I/O keys to include. If None, include all model I/O features.
        model_io_fp32_keys: Subset of FP32 model I/O keys to include. If None, include all FP32 model I/O features.

    Returns:
        dict: Dictionary containing all placement definitions organized by category.
            Contains both 3-tuple placements (for main device mesh) and 2-tuple cp placements
            (for cp submesh with preexisting batch dimension).
    """

    # Base placement patterns
    placements_single = (Shard(0), Shard(1), Replicate())
    placements_pair = (Shard(0), Shard(1), Shard(2))
    placements_scalar = (Shard(0), Replicate(), Replicate())

    # Base placement patterns for cp submesh with preexisting batch dimension
    placements_cp_single = (Shard(0), Replicate())
    placements_cp_pair = (Shard(0), Shard(1))

    # Helper function to convert 3-tuple placements to 2-tuple cp placements
    def convert_to_cp_placement(placement):
        if placement == placements_single:
            return placements_cp_single
        elif placement == placements_pair:
            return placements_cp_pair
        elif placement == (Shard(0), Shard(2), Replicate()):
            # Special case for coords in atom_features
            return (Shard(1), Replicate())
        else:
            raise ValueError(f"Unsupported placement pattern: {placement}")

    # 3-tuple placements for main device mesh
    placements_token_features_full = OrderedDict(
        {
            # Core features for InputEmbedder
            "token_pad_mask": placements_single,
            "res_type": placements_single,
            "mol_type": placements_single,
            "pocket_feature": placements_single,
            # Additional features for Boltz1
            "token_bonds": placements_pair,
            "type_bonds": placements_pair,
            "token_pair_pad_mask": placements_pair,  # this isn't returned with "window_batching = True"
            # Additional features for RelativePositionEncoder
            "asym_id": placements_single,
            "residue_index": placements_single,
            "entity_id": placements_single,
            "token_index": placements_single,
            "sym_id": placements_single,
            "cyclic_period": placements_single,
            # Additional features for distogram loss
            "disto_target": placements_pair,
            "token_disto_mask": placements_single,
            "disto_coords_ensemble": placements_single,
            # Boltz-2 token features
            "method_feature": placements_single,
            "modified": placements_single,
            "contact_conditioning": placements_pair,
            "contact_threshold": placements_pair,
            "frames_idx": placements_single,
        }
    )

    placements_msa_features_full = OrderedDict(
        {
            # Additional features for MSA module
            "msa": placements_pair,
            "has_deletion": placements_pair,
            "deletion_value": placements_pair,
            "msa_paired": placements_pair,
            "msa_mask": placements_pair,
            # Core features for InputEmbedder (MSA-derived)
            "profile": placements_single,
            "deletion_mean": placements_single,
        }
    )

    placements_atom_features_full = OrderedDict(
        {
            # Core features for InputEmbedder
            "atom_pad_mask": placements_single,
            "atom_to_token": placements_single,
            "pair_mask": placements_pair,  # this isn't returned with "window_batching = True"
            "ref_pos": placements_single,
            "ref_charge": placements_single,
            "ref_element": placements_single,
            "ref_atom_name_chars": placements_single,
            "ref_space_uid": placements_single,
            "atom_resolved_mask": placements_single,
            # Additional features for AtomDiffusion
            # Original Boltz-1x code processes "coords" to shape (B, 1, n_atoms, 3)
            "coords": (Shard(0), Shard(2), Replicate()),
            # Boltz-2 atom features
            "token_to_rep_atom": placements_single,
            "bfactor": placements_single,
            "plddt": placements_single,
        }
    )

    placements_model_io_full = OrderedDict(
        {
            "noise": placements_single,
            "denoised_atom_coords": placements_single,
            "d_denoised_atom_coords": placements_single,
            "aligned_true_atom_coords": placements_single,
            # Additional model I/O for DiffusionModule
            "r_noisy_expected": placements_single,
            "d_r_noisy_expected": placements_single,
            "r_update_expected": placements_single,
            "d_r_update_expected": placements_single,
            # Additional model I/O for AtomDiffusion preconditioned network
            "noised_atom_coords": placements_single,
            "denoised_atom_coords_expected": placements_single,
            "d_denoised_atom_coords_expected": placements_single,
            # Additional model I/O for AtomDiffusion sample
            "sample_atom_coords_expected": placements_single,
        }
    )

    placements_model_io_fp32_full = OrderedDict(
        {
            "denoised_atom_coords_fp32": placements_single,
            # Additional FP32 model I/O for DiffusionModule
            "r_update_fp32": placements_single,
            "d_r_noisy_fp32": placements_single,
            # Additional FP32 model I/O for AtomDiffusion sample
            "sample_atom_coords_fp32": placements_single,
        }
    )

    # Apply subsetting if specified
    placements_token_features = (
        OrderedDict({k: v for k, v in placements_token_features_full.items() if k in token_keys})
        if token_keys is not None
        else placements_token_features_full
    )

    placements_msa_features = (
        OrderedDict({k: v for k, v in placements_msa_features_full.items() if k in msa_keys})
        if msa_keys is not None
        else placements_msa_features_full
    )

    placements_atom_features = (
        OrderedDict({k: v for k, v in placements_atom_features_full.items() if k in atom_keys})
        if atom_keys is not None
        else placements_atom_features_full
    )

    placements_model_io = (
        OrderedDict({k: v for k, v in placements_model_io_full.items() if k in model_io_keys})
        if model_io_keys is not None
        else placements_model_io_full
    )

    placements_model_io_fp32 = (
        OrderedDict({k: v for k, v in placements_model_io_fp32_full.items() if k in model_io_fp32_keys})
        if model_io_fp32_keys is not None
        else placements_model_io_fp32_full
    )

    # 2-tuple placements for cp submesh with preexisting batch dimension
    # Generated using dictionary comprehension from their non-cp counterparts
    placements_cp_atom_features = OrderedDict(
        [
            # Add the additional key specific to cp variant first
            ("atom_counts_per_token", placements_cp_single),
            # Convert existing atom features maintaining original order
            *[(k, convert_to_cp_placement(v)) for k, v in placements_atom_features.items()],
        ]
    )

    placements_cp_model_io = OrderedDict({k: convert_to_cp_placement(v) for k, v in placements_model_io.items()})

    placements_cp_model_io_fp32 = OrderedDict(
        {k: convert_to_cp_placement(v) for k, v in placements_model_io_fp32.items()}
    )

    return {
        # Base patterns
        "single": placements_single,
        "pair": placements_pair,
        "scalar": placements_scalar,
        "cp_single": placements_cp_single,
        "cp_pair": placements_cp_pair,
        # Feature placements (3-tuple)
        "token_features": placements_token_features,
        "msa_features": placements_msa_features,
        "atom_features": placements_atom_features,
        "model_io": placements_model_io,
        "model_io_fp32": placements_model_io_fp32,
        # Feature placements (2-tuple, cp submesh)
        "cp_atom_features": placements_cp_atom_features,
        "cp_model_io": placements_cp_model_io,
        "cp_model_io_fp32": placements_cp_model_io_fp32,
    }


def pad_to_length(t: DTensor, dim: int, length: int) -> DTensor:
    """Pad a DTensor's local shards along *dim* so the global size equals *length*.

    ``distribute_atom_features`` applies intersperse padding per DP rank
    independently, but ``CollateDTensor`` additionally homogenizes local shard
    shapes across DP ranks via an all-reduce MAX.  When samples have different
    atom counts the dataloader features are homogenized but DTensors produced
    by ``distribute_atom_features`` are not.  This helper pads a DTensor's
    local shard along *dim* to match an externally-known correct global size
    (typically obtained from a homogenized batch feature such as
    ``feats["atom_pad_mask"].shape[-1]``).

    Unlike ``homogenize_shard_shapes`` — which derives the target from the
    DTensor's own (potentially inconsistent) global shape — this function
    accepts an explicit *length* that is authoritative across all ranks.

    Parameters
    ----------
    t : DTensor
        Input distributed tensor whose global ``t.shape[dim]`` may be smaller
        than *length*.
    dim : int
        Tensor dimension to pad (e.g. 1 for the atom dimension in
        ``[B*mult, n_atoms, 3]``).
    length : int
        Desired global size along *dim*.  Must be >= ``t.shape[dim]``.

    Returns
    -------
    DTensor
        A new DTensor with ``shape[dim] == length``.

    Raises
    ------
    ValueError
        If ``t.shape[dim] >= length`` (no padding needed — caller bug).
        If ``length`` or ``t.shape[dim]`` is not divisible by the mesh size.
        If multiple mesh axes shard the same tensor dimension.
    """
    if t.shape[dim] > length:
        raise ValueError(f"t.shape[{dim}]={t.shape[dim]} already exceeds target length={length}")
    if t.shape[dim] == length:
        raise ValueError(
            f"t.shape[{dim}]={t.shape[dim]} already equals target length={length}; "
            f"pad_to_length should not be called when no padding is needed"
        )

    # Find the single mesh axis that shards tensor dimension *dim*.
    mesh_size = 1
    for mesh_dim_idx, p in enumerate(t.placements):
        if isinstance(p, Shard) and p.dim == dim:
            if mesh_size != 1:
                raise ValueError(
                    f"pad_to_length does not support multiple mesh axes sharding "
                    f"the same tensor dimension {dim}. "
                    f"Placements: {t.placements}"
                )
            mesh_size = t.device_mesh.size(mesh_dim_idx)

    if t.shape[dim] % mesh_size != 0:
        raise ValueError(
            f"t.shape[{dim}]={t.shape[dim]} is not divisible by mesh size {mesh_size} along the axis sharding dim {dim}"
        )
    if length % mesh_size != 0:
        raise ValueError(f"length={length} is not divisible by mesh size {mesh_size} along the axis sharding dim {dim}")

    local = t.to_local()
    local_target = length // mesh_size
    pad_amount = local_target - local.shape[dim]
    if pad_amount <= 0:
        raise ValueError(
            f"Local shard shape[{dim}]={local.shape[dim]} already >= "
            f"local target {local_target} (global length={length}, "
            f"mesh_size={mesh_size}); no padding possible"
        )

    pad_spec = [0] * (2 * local.ndim)
    pad_spec[2 * (local.ndim - 1 - dim) + 1] = pad_amount
    local = torch.nn.functional.pad(local, pad_spec)

    gshape = list(t.shape)
    gshape[dim] = length
    gshape = torch.Size(gshape)

    return DTensor.from_local(
        local,
        t.device_mesh,
        t.placements,
        shape=gshape,
        stride=LayoutRightMap(tuple(gshape)).strides,
    )


def homogenize_shard_shapes(input: DTensor, value_to_pad: Any | None = None) -> DTensor:
    """Homogenize shard shapes across all ranks by padding local shards to a consistent size.

    NOTE: the involved padding is always towards the end (or the last element) along each
    tensor axis

    In distributed tensor operations, different ranks may have slightly different local shard
    sizes due to uneven data distribution. This function ensures all ranks have consistent
    local shard shapes by padding smaller shards to match the target size.

    The target shape for each tensor dimension is determined by:
    - For non-sharded dimensions: Same as the global tensor shape
    - For sharded dimensions: global_size // mesh_size_along_sharding_mesh_dimension

    All ranks participate in this operation, even if some don't require padding, to ensure
    proper synchronization for subsequent collective operations like DTensor.from_local().

    Args:
        input (DTensor): The input distributed tensor to homogenize
        value_to_pad (Any | None, optional): Value to use for padding. Defaults to 0.

    Returns:
        DTensor: A new DTensor with homogenized local shard shapes across all ranks

    Raises:
        ValueError: If any local shard dimension is larger than the computed target size

    Example:
        >>> # Ranks have different local shard sizes: [10, 8] and [10, 9]
        >>> # After homogenization: both ranks have [10, 9]
        >>> homogenized_dtensor = homogenize_shard_shapes(input_dtensor)

    Note:
        This function is particularly useful in testing scenarios where you need
        consistent shard shapes across ranks for reliable distributed computations.
    """
    # Get the local shard and its shape
    local_shard = input.to_local()
    local_shape = torch.tensor(local_shard.shape, dtype=torch.int64)

    # Get DTensor properties
    global_shape = list(input.shape)
    placements = input.placements
    device_mesh = input.device_mesh

    # Calculate target shape for the local shard
    # Start with global shape and modify only sharded dimensions
    target_shape = torch.tensor(global_shape, dtype=torch.int64)

    for mesh_dim_idx, placement in enumerate(placements):
        if isinstance(placement, Partial):
            raise ValueError(f"Partial placements are not supported: {placement}")
        if isinstance(placement, Shard):
            # For sharded dimensions, target size is global_size / mesh_size_for_this_mesh_dim
            tensor_dim = placement.dim  # This is the tensor dimension being sharded
            mesh_size = device_mesh.size(mesh_dim_idx)
            target_shape[tensor_dim] = global_shape[tensor_dim] // mesh_size

    # Calculate padding amounts using tensor operations (all ranks must compute this even if they don't need padding)
    padding_amounts = target_shape - local_shape

    # Check for invalid cases where local shard is larger than target
    if (padding_amounts < 0).any():
        raise ValueError(f"Local shard shape {local_shape} has axes larger than target shape {target_shape}")

    # Check if this rank needs padding
    needs_padding = (padding_amounts > 0).any()

    if needs_padding:
        # Apply padding to local shard
        # torch.nn.functional.pad expects padding in reverse order (last dim first)
        # Create pad_values as (ndim, 2) tensor: [pad_left, pad_right] for each dimension
        pad_values = torch.zeros((local_shard.ndim, 2), dtype=torch.int64)
        # always pad towards the end (or the last element) along each tensor axis
        pad_values[:, 1] = padding_amounts  # Set right padding amounts

        # Flatten and reverse for torch.nn.functional.pad format
        pad_values = pad_values.flip(0).flatten().tolist()

        if value_to_pad is None:
            value_to_pad = torch.tensor(0, dtype=local_shard.dtype)

        padded_local_shard = torch.nn.functional.pad(local_shard, pad_values, value=value_to_pad)
    else:
        # No padding needed for this rank, but still participate in collective operation
        padded_local_shard = local_shard

    # By definition, this function ensures homogeneous shape across ranks
    shape_output_global, stride_output_global = map(
        tuple,
        compute_global_tensor_info(padded_local_shard, device_mesh, placements),
    )

    # Create new DTensor from padded local shard
    return DTensor.from_local(
        padded_local_shard,
        device_mesh=device_mesh,
        placements=placements,
        shape=shape_output_global,
        stride=stride_output_global,
        run_check=False,  # Skip validation for performance
    )


def concat_data(out_dir: Path, *datas: Path) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    # manifest.json  msa  structures
    # 1. copy the msa contents into msa, raising an error if there are duplicate filenames
    # 2. copy the structures contents into structures, raising an error if there are duplicate filenames
    # 3. merge the manifest.json files, raising an error if there are duplicate filenames
    # 4. write the merged manifest.json
    msa_dir = out_dir / "msa"
    msa_dir.mkdir(parents=True, exist_ok=True)
    copied = set()
    if isinstance(datas, (Path, str)):
        data_lst: list[Path] = [Path(datas)]
    else:
        data_lst = [Path(data) for data in datas]
    for data in data_lst:
        for file in (data / "msa").glob("*"):
            if file.name in copied:
                raise ValueError(f"Duplicate MSA file {file.name}")
            shutil.copy(file, msa_dir / file.name)
            copied.add(file.name)
    structures_dir = out_dir / "structures"
    structures_dir.mkdir(parents=True, exist_ok=True)
    copied_structures = set()
    manifests = []
    for data in data_lst:
        for file in (data / "structures").glob("*"):
            if file.name in copied_structures:
                raise ValueError(f"Duplicate structure file {file.name}")
            shutil.copy(file, structures_dir / file.name)
            copied_structures.add(file.name)
        with open(data / "manifest.json", "r") as f:
            manifest = json.load(f)
            manifests.append(manifest)
    assert all(manifest.keys() == manifests[0].keys() for manifest in manifests), "Manifest keys do not match"
    assert all(set(manifest.keys()) == {"records"} for manifest in manifests), "Manifest keys do not match"
    records = []
    for manifest in manifests:
        records.extend(manifest["records"])
    manifest = {"records": records}
    manifest_file = out_dir / "manifest.json"
    manifest_file.write_text(json.dumps(manifest))
    return out_dir


@contextmanager
def pytorch_use_deterministic_ops():
    """Context manager to enable PyTorch deterministic algorithms in spawn processes.

    This context manager enables deterministic behavior in PyTorch operations
    for reproducibility in testing or debugging. It sets the CUBLAS_WORKSPACE_CONFIG
    environment variable and enables PyTorch's deterministic algorithms mode.

    Important: This context manager must only be used in spawn parallel processes
    (not the main process) because the CUBLAS_WORKSPACE_CONFIG environment variable
    affects the underlying CUDA context for the entire process lifetime. Restricting
    usage to spawned processes prevents side effects on other tests.

    The context manager automatically restores the previous deterministic setting
    and CUBLAS_WORKSPACE_CONFIG value upon exit.

    Yields:
        None

    Raises:
        RuntimeError: If called in the main process (not a spawn parallel process).

    Example:
        >>> # In a spawned test process
        >>> with pytorch_use_deterministic_ops():
        ...     # All PyTorch operations use deterministic algorithms
        ...     result = model(input_tensor)
    """
    # technically, the CUBLAS_WORKSPACE_CONFIG must be set before the
    # "import torch" statement and its effect on the underlying CUDA
    # context will last until the process ends. For that reason, to
    # exclude this env variable's side effects on other tests, we
    # need to restrict this context manager to a spawn parallel
    # processes
    is_spawn_process = multiprocessing.parent_process() is not None
    if is_spawn_process:
        raise RuntimeError("pytorch_use_deterministic_ops() can only be used in spawn parallel processes")
    deterministic_restore = torch.are_deterministic_algorithms_enabled()
    original_env = os.environ.get("CUBLAS_WORKSPACE_CONFIG", None)
    try:
        os.environ["CUBLAS_WORKSPACE_CONFIG"] = ":4096:8"
        torch.use_deterministic_algorithms(True)
        yield
    finally:
        torch.use_deterministic_algorithms(deterministic_restore)
        if original_env is not None:
            os.environ["CUBLAS_WORKSPACE_CONFIG"] = original_env
        else:
            os.environ.pop("CUBLAS_WORKSPACE_CONFIG", None)


@contextmanager
def benchmark_peak_memory_and_runtime():
    """Context manager to benchmark peak memory usage and runtime of a code block.

    This context manager tracks the peak CUDA memory allocated during the execution
    of the wrapped code block and measures the wall-clock time execution duration.
    It yields a dictionary that will be populated with 'peak_mem' (in MB) and 'time'
    (in ms) keys after the block execution completes.

    The memory measurement tracks the *peak allocated memory* relative to the memory
    allocated at the start of the context, attempting to isolate the memory usage
    of the specific operations within the block.

    Yields:
        dict: A dictionary that will contain results after context exit:
            - "peak_mem" (float): Peak memory usage in MB.
            - "time" (float): Execution time in milliseconds.

    Note:
        - Requires CUDA to be available.
        - Performs `torch.cuda.synchronize()` before stopping the timer to ensure
          accurate GPU timing.
        - Clears cache and resets peak memory stats at the beginning.

    Example:
        >>> with benchmark_peak_memory_and_runtime() as stats:
        ...     model(input)
        >>> print(f"Memory: {stats['peak_mem']} MB, Time: {stats['time']} ms")
    """
    torch.cuda.empty_cache()
    torch.cuda.reset_peak_memory_stats()
    start_mem = torch.cuda.memory_allocated()
    start_time = time.time()

    stats = {}
    yield stats

    torch.cuda.synchronize()
    end_time = time.time()
    peak_mem = torch.cuda.max_memory_allocated()

    # We want the peak memory *induced* by the function, relative to start.
    peak_usage_mb = (peak_mem - start_mem) / 1024 / 1024
    duration_ms = (end_time - start_time) * 1000

    stats["peak_mem"] = peak_usage_mb
    stats["time"] = duration_ms


def pad_or_shrink_to_length(
    tensor: torch.Tensor, axis: int, target_length: int, pad_value: float = 0.0
) -> torch.Tensor:
    """Pad or shrink tensor along the specified axis to target_length.

    When shrinking, the function slices from the beginning of the axis.
    When padding, zeros (or pad_value) are appended at the end.

    Args:
        tensor: Input tensor to resize.
        axis: The dimension along which to pad or shrink.
        target_length: The desired length along the axis.
        pad_value: Value to use for padding (default: 0.0).

    Returns:
        Tensor with shape[axis] == target_length.
    """
    current_length = tensor.shape[axis]
    if current_length == target_length:
        return tensor
    elif current_length > target_length:
        # Shrink: slice to target_length
        slices = [slice(None)] * tensor.ndim
        slices[axis] = slice(0, target_length)
        return tensor[tuple(slices)]
    else:
        # Pad: append zeros
        pad_shape = list(tensor.shape)
        pad_shape[axis] = target_length - current_length
        padding = torch.full(pad_shape, pad_value, dtype=tensor.dtype, device=tensor.device)
        return torch.cat([tensor, padding], dim=axis)


def distribute_atom_features(
    inputs: dict[str, Tensor],
    placements_cp: dict[str, tuple],
    placements_dp_cp: dict[str, tuple],
    device_mesh: DeviceMesh,
    cp_group: dist.ProcessGroup,
    cp_mesh_dim_names: tuple[str, ...] = ("cp_axis_0", "cp_axis_1"),
    dp_mesh_dim_name: str = "dp",
    multiplicities: dict[str, int] | None = None,
) -> dict[str, DTensor]:
    """Distribute atom features across a device mesh with intersperse padding.

    This utility abstracts the common workflow of:
    1. Calling pad_and_scatter_atom_features_dtensor() per DP rank to shard across CP
    2. Collating the batch along DP ranks
    3. Combining per-multiplicity DTensors into single DTensor with flattened batch*mult

    Parameters
    ----------
    inputs : dict[str, Tensor]
        Input tensors on host (CPU). Shape [batch, ...] for each feature.
        For multiplicity features, use separate keys like "feat_0", "feat_1".
        Batch size must equal the DP world size.
        May include auxiliary features (e.g., "atom_counts_per_token") needed
        by pad_and_scatter_atom_features_dtensor but not returned in output.
    placements_cp : dict[str, tuple]
        Placements for CP submesh (e.g., (Shard(0), Replicate())).
        Must have the same keys as inputs.
    placements_dp_cp : dict[str, tuple]
        Placements for full mesh (e.g., (Shard(0), Shard(1), Replicate())).
        Keys must be a subset of inputs.keys(). Only features with keys in
        placements_dp_cp will be returned in the output.
    device_mesh : DeviceMesh
        Full device mesh (e.g., 3D: dp, cp_0, cp_1).
    cp_group : dist.ProcessGroup
        Flattened CP process group for this DP slice.
    cp_mesh_dim_names : tuple[str, ...], optional
        Names of CP dimensions in the mesh. Default ("cp_axis_0", "cp_axis_1").
    dp_mesh_dim_name : str, optional
        Name of DP dimension in the mesh. Default "dp".
    multiplicities : dict[str, int] | None, optional
        Features with multiplicity. Key is base name (without _0, _1 suffix),
        value is the multiplicity count. These will be combined in output.

    Returns
    -------
    dict[str, DTensor]
        DTensors distributed across the full mesh. Only features with keys in
        placements_dp_cp are returned. Multiplicity features are combined with
        shape [batch*mult, ...].

    Raises
    ------
    ValueError
        If keys don't match requirements, batch size != DP world size,
        or local batch size != 1.

    Examples
    --------
    >>> # For features with multiplicity=2:
    >>> inputs = {
    ...     "atom_counts_per_token": counts,  # [B, n_tokens] - auxiliary, not in output
    ...     "token_to_rep_atom": tensor_a,    # [B, n_tokens, n_atoms]
    ...     "resolved_mask_0": tensor_b0,     # [B, n_atoms]
    ...     "resolved_mask_1": tensor_b1,     # [B, n_atoms]
    ... }
    >>> placements_cp = {
    ...     "atom_counts_per_token": (Shard(0), Replicate()),
    ...     "token_to_rep_atom": (Shard(0), Replicate()),
    ...     "resolved_mask_0": (Shard(0), Replicate()),
    ...     "resolved_mask_1": (Shard(0), Replicate()),
    ... }
    >>> placements_dp_cp = {
    ...     "token_to_rep_atom": (Shard(0), Shard(1), Replicate()),
    ...     "resolved_mask_0": (Shard(0), Shard(1), Replicate()),
    ...     "resolved_mask_1": (Shard(0), Shard(1), Replicate()),
    ... }
    >>> multiplicities = {"resolved_mask": 2}
    >>> result = distribute_atom_features(
    ...     inputs, placements_cp, placements_dp_cp, device_mesh, cp_group,
    ...     multiplicities=multiplicities
    ... )
    >>> # result["token_to_rep_atom"] has shape [B, n_tokens_padded, n_atoms_padded]
    >>> # result["resolved_mask"] has shape [B*2, n_atoms_padded]
    >>> # "atom_counts_per_token" is NOT in result
    """
    from boltz.distributed.data.feature.featurizer import pad_and_scatter_atom_features_dtensor

    multiplicities = multiplicities or {}

    # Validate key consistency
    # inputs and placements_cp must have the same keys
    if inputs.keys() != placements_cp.keys():
        raise ValueError(
            f"inputs and placements_cp must have the same keys. "
            f"inputs: {set(inputs.keys())}, placements_cp: {set(placements_cp.keys())}"
        )
    # placements_dp_cp keys must be a subset of inputs keys
    if not placements_dp_cp.keys() <= inputs.keys():
        raise ValueError(
            f"placements_dp_cp keys must be a subset of inputs keys. "
            f"placements_dp_cp: {set(placements_dp_cp.keys())}, inputs: {set(inputs.keys())}"
        )

    # Validate multiplicities keys exist in placements_dp_cp (since they'll be in output)
    for base_name, mult in multiplicities.items():
        for i in range(mult):
            key = f"{base_name}_{i}"
            if key not in placements_dp_cp:
                raise ValueError(f"Multiplicity key '{key}' not found in placements_dp_cp")

    # Get mesh info
    cp_mesh = device_mesh[cp_mesh_dim_names]
    dp_dim_idx = device_mesh.mesh_dim_names.index(dp_mesh_dim_name)
    dp_rank = device_mesh.get_coordinate()[dp_dim_idx]
    dp_world_size = device_mesh.size(dp_dim_idx)

    # Get batch size and validate
    sample_key = next(iter(inputs.keys()))
    batch_size = inputs[sample_key].shape[0]
    if batch_size != dp_world_size:
        raise ValueError(f"Batch size ({batch_size}) must equal DP world size ({dp_world_size})")
    batch_size_per_dp = batch_size // dp_world_size
    if batch_size_per_dp != 1:
        raise ValueError(
            f"Local batch size must be 1, got {batch_size_per_dp}. "
            f"pad_and_scatter_atom_features_dtensor can only process one sample at a time."
        )

    # Get CP group info
    cp_rank = dist.get_rank(cp_group)
    cp_src_rank = dist.get_process_group_ranks(cp_group)[0]

    # Determine device
    if device_mesh.device_type == "cuda":
        local_device_idx = torch.cuda.current_device()
        device = torch.device("cuda", local_device_idx)
    else:
        device = torch.device(device_mesh.device_type)

    # Prepare inputs for scatter: select sample for this DP rank
    # Only CP rank 0 provides inputs; others pass None
    _ENSEMBLE_EXPECTED_NDIM = {"frames_idx": 3, "frame_resolved_mask": 2}
    if cp_rank == 0:
        inputs_for_scatter = {}
        for k, v in inputs.items():
            val = v[dp_rank].to(device=device)
            expected_ndim = _ENSEMBLE_EXPECTED_NDIM.get(k)
            if expected_ndim is not None and val.ndim < expected_ndim:
                val = val.unsqueeze(0)
            inputs_for_scatter[k] = val
    else:
        inputs_for_scatter = None

    # Scatter across CP ranks with intersperse padding
    feats_cp = pad_and_scatter_atom_features_dtensor(
        inputs_for_scatter,
        placements_cp,
        cp_group,
        cp_src_rank,
        cp_mesh,
    )

    # Compute global shape/stride and create DTensors on full mesh
    # Only process keys that are in placements_dp_cp (output features)
    def _local_for_dp_wrap(k, v):
        """Get local tensor ready for DP wrapping (unsqueeze batch dim).

        For ensemble-aware features (frames_idx, frame_resolved_mask), the
        featurizer requires an E=1 ensemble dim on input but downstream
        consumers expect (T, ...) per sample.  Squeeze the ensemble dim
        before adding the DP batch dim.
        """
        local = v.to_local()
        if k in _ENSEMBLE_EXPECTED_NDIM:
            local = local.squeeze(0)  # (1, T_shard, ...) -> (T_shard, ...)
        return local.unsqueeze(0)  # Add batch dim

    feats_shape_stride = {
        k: tuple(
            map(
                tuple,
                compute_global_tensor_info(_local_for_dp_wrap(k, v), device_mesh, placements_dp_cp[k]),
            )
        )
        for k, v in feats_cp.items()
        if k in placements_dp_cp
    }

    feats: dict[str, DTensor] = {
        k: DTensor.from_local(
            _local_for_dp_wrap(k, v),  # Add batch dim (squeeze ensemble dim first if needed)
            device_mesh,
            placements_dp_cp[k],
            shape=feats_shape_stride[k][0],
            stride=feats_shape_stride[k][1],
        )
        for k, v in feats_cp.items()
        if k in placements_dp_cp
    }

    # Combine per-multiplicity DTensors into single DTensor with flattened batch*mult
    for base_name, mult in multiplicities.items():
        # Pop the per-multiplicity DTensors
        mult_tensors = [feats.pop(f"{base_name}_{i}") for i in range(mult)]

        # Validate all tensors have same shape, placements, device_mesh
        if not all(t.shape == mult_tensors[0].shape for t in mult_tensors):
            raise ValueError(f"All multiplicity tensors for '{base_name}' must have the same shape")
        if not all(t.placements == mult_tensors[0].placements for t in mult_tensors):
            raise ValueError(f"All multiplicity tensors for '{base_name}' must have the same placements")
        if not all(t.device_mesh == mult_tensors[0].device_mesh for t in mult_tensors):
            raise ValueError(f"All multiplicity tensors for '{base_name}' must have the same device_mesh")

        # Combine: stack along dim 1, then flatten dims 0 and 1
        # [B, ...] -> stack -> [B, mult, ...] -> [B*mult, ...]
        local_cat = torch.cat([t.to_local().unsqueeze(1) for t in mult_tensors], dim=1)
        local_flat = local_cat.flatten(0, 1)

        feats[base_name] = DTensor.from_local(
            local_flat,
            mult_tensors[0].device_mesh,
            mult_tensors[0].placements,
        )

    return feats


def make_random_contact_conditioning_features(
    B: int,
    N: int,
    num_cc_types: int,
    dtype: torch.dtype = torch.float32,
    device: str = "cpu",
    seed: int = 42,
) -> tuple[Tensor, Tensor]:
    """Create random contact conditioning features for testing.

    Exercises all three masking branches (UNSPECIFIED, UNSELECTED, active)
    so that both ``encoding_unspecified`` and ``encoding_unselected``
    parameters receive non-zero gradients.

    Parameters
    ----------
    B : int
        Batch size.
    N : int
        Number of tokens.
    num_cc_types : int
        Number of contact conditioning types (``len(const.contact_conditioning_info)``).
    dtype : torch.dtype
        Data type for the tensors.
    device : str
        Device for the tensors.
    seed : int
        Random seed for reproducibility.

    Returns
    -------
    cc : Tensor
        Contact conditioning features, shape ``(B, N, N, num_cc_types)``.
    ct : Tensor
        Contact threshold values, shape ``(B, N, N)``.
    """
    with torch.random.fork_rng():
        torch.manual_seed(seed)
        rand_type = torch.rand(B, N, N, device=device)
        cc = torch.zeros(B, N, N, num_cc_types, dtype=dtype, device=device)
        cc[:, :, :, 0] = (rand_type < 0.3).to(dtype)  # UNSPECIFIED (~30%)
        cc[:, :, :, 1] = ((rand_type >= 0.3) & (rand_type < 0.5)).to(dtype)  # UNSELECTED (~20%)
        if num_cc_types > 4:  # noqa: PLR2004
            cc[:, :, :, 4] = (rand_type >= 0.5).to(dtype)  # CONTACT (~50%)
        ct = (torch.rand(B, N, N, device=device) * 22.0).to(dtype)
    return cc, ct


def _extract_output_dtypes(
    output: Any,
    prefix: str,
) -> dict[str, torch.dtype]:
    """Extract dtypes from a module's forward output.

    Handles plain tensors, DTensors, tuples/lists of tensors, and dicts
    mapping to tensors.  Returns ``{qualified_name: dtype}`` entries.
    """
    if isinstance(output, Tensor):
        return {prefix: output.dtype}
    if isinstance(output, dict):
        result: dict[str, torch.dtype] = {}
        for k, v in output.items():
            if isinstance(v, Tensor):
                result[f"{prefix}/{k}"] = v.dtype
        return result
    if isinstance(output, (tuple, list)):
        result = {}
        for i, v in enumerate(output):
            if isinstance(v, Tensor):
                result[f"{prefix}/{i}"] = v.dtype
        return result
    return {}


class DtypeProfiler:
    """Capture dtypes of module outputs, parameters, and parameter gradients.

    Attach to a model before a forward + backward pass.  After the pass,
    the three dictionaries ``fwd_dtypes``, ``param_dtypes``, and
    ``param_grad_dtypes`` contain ``{qualified_name: torch.dtype}`` entries
    for every module output, parameter, and parameter gradient respectively.

    Usage::

        profiler = DtypeProfiler(model)
        loss = model(batch).sum()
        loss.backward()
        profiler.collect_grad_dtypes(model)
        profiler.remove_hooks()

    Works transparently with both plain ``torch.Tensor`` and ``DTensor``
    outputs (``DTensor.dtype`` returns the element dtype).
    """

    def __init__(self, model: torch.nn.Module) -> None:
        self.fwd_dtypes: dict[str, torch.dtype] = {}
        self.param_dtypes: dict[str, torch.dtype] = {}
        self.param_grad_dtypes: dict[str, torch.dtype] = {}
        self._handles: list[torch.utils.hooks.RemovableHook] = []
        self._register(model)

    def _make_fwd_hook(self, name: str):
        def hook(_module: torch.nn.Module, _input: Any, output: Any) -> None:
            self.fwd_dtypes.update(_extract_output_dtypes(output, name))

        return hook

    def _register(self, model: torch.nn.Module) -> None:
        for name, param in model.named_parameters():
            self.param_dtypes[name] = param.dtype
        for name, module in model.named_modules():
            self._handles.append(module.register_forward_hook(self._make_fwd_hook(name)))

    def collect_grad_dtypes(self, model: torch.nn.Module) -> None:
        """Snapshot ``param.grad.dtype`` for every parameter with a gradient."""
        for name, param in model.named_parameters():
            if param.grad is not None:
                self.param_grad_dtypes[name] = param.grad.dtype

    def remove_hooks(self) -> None:
        """Remove all registered forward hooks."""
        for h in self._handles:
            h.remove()
        self._handles.clear()

    def __enter__(self) -> "DtypeProfiler":
        return self

    def __exit__(self, *args: Any) -> None:
        self.remove_hooks()


class RecomputeProfiler:
    """Count how many times each module's forward is invoked.

    When activation checkpointing is active, ``torch.utils.checkpoint.checkpoint``
    re-runs the checkpointed function during backward.  Because
    ``nn.Module.__call__`` fires registered forward **pre**-hooks on every
    invocation, modules inside a checkpointed region will have
    ``fwd_counts[name] >= 2`` (once in the original forward, once during
    recomputation) while modules outside checkpointed regions will have
    ``fwd_counts[name] == 1``.

    Pre-hooks (``register_forward_pre_hook``) are used instead of post-forward
    hooks because ``use_reentrant=False`` checkpointing raises an internal
    ``_StopRecomputationError`` via ``early_stop`` to halt recomputation as
    soon as all needed tensors are regenerated.  This exception interrupts
    ``Module.forward()``, so post-forward hooks never fire for modules at or
    after the stop point.  Pre-hooks fire at the *start* of
    ``Module.__call__`` — before ``forward()`` — so they are unaffected.

    Usage::

        profiler = RecomputeProfiler(model)
        loss = model(batch).sum()
        loss.backward()
        profiler.remove_hooks()
        print(profiler.recomputed_modules)

    Works transparently with both plain ``torch.nn.Module`` and DTensor-wrapped
    models (hooks are registered on the local module hierarchy).
    """

    def __init__(self, model: torch.nn.Module) -> None:
        self.fwd_counts: dict[str, int] = {}
        self._handles: list[torch.utils.hooks.RemovableHook] = []
        for name, module in model.named_modules():
            self._handles.append(module.register_forward_pre_hook(self._make_hook(name)))

    def _make_hook(self, name: str):
        def hook(_module: torch.nn.Module, _input: Any) -> None:
            self.fwd_counts[name] = self.fwd_counts.get(name, 0) + 1

        return hook

    @property
    def recomputed_modules(self) -> frozenset[str]:
        """Module names whose forward was called >= 2 times (recomputed by checkpoint)."""
        return frozenset(n for n, c in self.fwd_counts.items() if c >= 2)

    def remove_hooks(self) -> None:
        """Remove all registered forward hooks."""
        for h in self._handles:
            h.remove()
        self._handles.clear()


# ---------------------------------------------------------------------------
# CCD ligand feature loading for PB metric tests
# ---------------------------------------------------------------------------

LIGAND_KEYS = (
    "ligand_edge_index",
    "ligand_edge_lower_bounds",
    "ligand_edge_upper_bounds",
    "ligand_edge_bond_mask",
    "ligand_edge_angle_mask",
    "ligand_chiral_atom_index",
    "ligand_chiral_check_mask",
    "ligand_chiral_atom_orientations",
    "ligand_stereo_bond_index",
    "ligand_stereo_check_mask",
    "ligand_stereo_bond_orientations",
    "ligand_aromatic_5_ring_index",
    "ligand_aromatic_6_ring_index",
    "ligand_planar_double_bond_index",
)


def load_ligand_features_from_ccd(
    mols_dir: str,
    mol_name: str = "PHE",
    atom_offset: int = 0,
) -> tuple[torch.Tensor, dict]:
    """Load ligand geometry features from a CCD pickle file.

    Uses ``load_molecules`` + ``get_symmetries`` from ``boltz.data.mol``
    to read the pre-computed PB feature arrays (edge index, distance
    bounds, chirality, stereo, aromatics) stored inside the RDKit Mol.

    Parameters
    ----------
    mols_dir : str
        Path to the directory containing per-residue ``.pkl`` files
        (e.g. ``tests/test_data/data/mols``).
    mol_name : str
        CCD component name (default ``"PHE"``).
    atom_offset : int
        Offset added to all atom-index features so they refer to the
        correct position within the full-structure atom array.

    Returns
    -------
    coords : torch.Tensor
        Ideal 3-D coordinates, shape ``(n_atoms, 3)``, ``float32``.
    features : dict[str, torch.Tensor]
        Dictionary keyed by the strings in ``LIGAND_KEYS``.
    """
    from boltz.data.mol import get_symmetries, load_molecules

    mols = load_molecules(mols_dir, [mol_name])
    mol = mols[mol_name]
    syms = get_symmetries(mols)
    (
        _syms_ccd, _names_ccd,
        edge_index, lower_bounds, upper_bounds, bond_mask, angle_mask,
        chiral_atom_index, chiral_check_mask, chiral_atom_orientations,
        stereo_bond_index, stereo_check_mask, stereo_bond_orientations,
        aromatic_5_ring_index, aromatic_6_ring_index, planar_double_bond_index,
    ) = syms[mol_name]

    conf = mol.GetConformer(0)
    coords = torch.tensor(
        [[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y, conf.GetAtomPosition(i).z]
         for i in range(mol.GetNumAtoms())],
        dtype=torch.float32,
    )

    features = {
        "ligand_edge_index": torch.tensor(edge_index, dtype=torch.long) + atom_offset,
        "ligand_edge_lower_bounds": torch.tensor(lower_bounds, dtype=torch.float32),
        "ligand_edge_upper_bounds": torch.tensor(upper_bounds, dtype=torch.float32),
        "ligand_edge_bond_mask": torch.tensor(bond_mask),
        "ligand_edge_angle_mask": torch.tensor(angle_mask),
        "ligand_chiral_atom_index": torch.tensor(chiral_atom_index, dtype=torch.long) + atom_offset,
        "ligand_chiral_check_mask": torch.tensor(chiral_check_mask),
        "ligand_chiral_atom_orientations": torch.tensor(chiral_atom_orientations),
        "ligand_stereo_bond_index": torch.tensor(stereo_bond_index, dtype=torch.long) + atom_offset,
        "ligand_stereo_check_mask": torch.tensor(stereo_check_mask),
        "ligand_stereo_bond_orientations": torch.tensor(stereo_bond_orientations),
        "ligand_aromatic_5_ring_index": torch.tensor(aromatic_5_ring_index, dtype=torch.long) + atom_offset,
        "ligand_aromatic_6_ring_index": torch.tensor(aromatic_6_ring_index, dtype=torch.long) + atom_offset,
        "ligand_planar_double_bond_index": torch.tensor(planar_double_bond_index, dtype=torch.long) + atom_offset,
    }
    return coords, features


def make_pb_test_data(n_tok, n_atom, mols_dir, mol_name="PHE", batch_size=1, n_samples=2, seed=42):
    """Create test data with a real CCD ligand for PB metric tests.

    Loads the ligand from the CCD pickle file at ``mols_dir/mol_name.pkl``
    via :func:`load_ligand_features_from_ccd`.

    Layout per batch element (``atoms_per_tok = n_atom // n_tok``):
        Tokens  0 .. n_tok-4 : PROTEIN (asym_id=0)
        Tokens  n_tok-3 .. n_tok-1 : NONPOLYMER / ligand (asym_id=1)

    The 3 ligand tokens x ``atoms_per_tok`` atoms must be >= the number
    of heavy atoms in the CCD component.  CCD ideal coordinates are
    placed at the ligand atom positions in ``sample_atom_coords``;
    remaining protein positions are filled with random noise.
    """
    atoms_per_tok = n_atom // n_tok
    rng = torch.Generator().manual_seed(seed)

    feats = random_features(
        size_batch=batch_size,
        n_tokens=n_tok,
        n_atoms=n_atom,
        n_msa=1,
        atom_counts_per_token_range=(atoms_per_tok, atoms_per_tok),
        device=torch.device("cpu"),
        float_value_range=(-1.0, 1.0),
        selected_keys=["atom_to_token", "mol_type", "atom_pad_mask", "atom_counts_per_token", "asym_id"],
        rng=rng,
    )

    batch: dict = {}
    for k, v in feats.items():
        batch[k] = v.to(torch.float32) if v.is_floating_point() else v
    batch["atom_to_token"] = batch["atom_to_token"].to(torch.float32)

    lig_start_tok = n_tok - 3
    batch["mol_type"][:, :lig_start_tok] = const.chain_type_ids["PROTEIN"]
    batch["mol_type"][:, lig_start_tok:] = const.chain_type_ids["NONPOLYMER"]

    batch["asym_id"][:, :lig_start_tok] = 0
    batch["asym_id"][:, lig_start_tok:] = 1

    lig_atom_start = lig_start_tok * atoms_per_tok

    lig_coords, lig_feats = load_ligand_features_from_ccd(mols_dir, mol_name, atom_offset=lig_atom_start)
    n_lig_atoms = lig_coords.shape[0]

    sample_coords = torch.randn(batch_size * n_samples, n_atom, 3, generator=rng)
    for s in range(batch_size * n_samples):
        sample_coords[s, lig_atom_start : lig_atom_start + n_lig_atoms] = lig_coords

    for k in LIGAND_KEYS:
        batch[k] = [lig_feats[k].clone() for _ in range(batch_size)]

    out = {"sample_atom_coords": sample_coords}
    return batch, out
