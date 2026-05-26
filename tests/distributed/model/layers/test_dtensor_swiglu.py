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


"""Tests a single instance of SwiGLU with DTensor

Verification requirements

    V1: single-proc FW input tensor values unchanged by FW and BW
    V2: single-proc BW input tensor values unchanged by BW
    V3: single-proc FW input tensor grads are zero at padded locations (virtual atoms)
        - for input tensors that require grads

    V4: multi-proc version of V1
    V5: multi-proc version of V2
    V6: multi-proc version of V3: implied by V3 and V9

    V7: multi-proc FW input tensor values and meta match single-proc inputs
    V8: multi-proc FW output tensor values close-to single-proc
    V9: multi-proc FW input gradient values close-to single-proc
    V10: multi-proc parameter gradient values close-to single-proc
    V11: multi-proc parameter gradient values identical across proc's

Implementation status
    V1: implemented
    V2: implemented
    V3: NA (no padding in SwiGLU)
    V4: implemented
    V5: implemented
    V6: NA (no padding in SwiGLU)
    V7: same data
    V8: implemented
    V9: implemented
    V10: NA (no parameters in SwiGLU)
    V11: NA (no parameters in SwiGLU)

Assertion threshold defaults for pytorch
"""

import itertools

import pytest
import torch
from torch import Tensor
from torch.distributed.tensor import DTensor, Replicate, Shard, distribute_tensor
from torch.testing import assert_close

from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.layers.swiglu import SwiGLU as SwiGLUWithDTensor
from boltz.model.modules.utils import SwiGLU as SwiGluBoltz
from boltz.testing.utils import (
    assert_tensors_identical,
    seed_by_rank,
    skip_if_cuda_not_avail_or_device_count_less_than_word_size,
    spawn_multiprocessing,
)

SEED = 42


def assert_swiglu_with_dtensor_fw_bw(
    rank: int,  # noqa
    input_example: Tensor,
    output_ref: Tensor,
    output_grad_example: Tensor,
    input_grads_ref: Tensor,
    grid_group_sizes: tuple[int, ...],
    device_type: str,
    backend: str,
    env_per_rank: dict[str, str],  # noqa
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

    # -------------------------------------------------------------
    # Move inputs and ref outputs to device
    # --------------------------------------------------------------
    input_example_device = input_example.detach().to(manager.device, copy=True)
    output_ref_device = output_ref.detach().to(manager.device, copy=True)
    input_grads_ref_device = input_grads_ref.detach().to(manager.device, copy=True)
    output_grad_example = output_grad_example.detach().to(manager.device, copy=True)

    # -------------------------------------------------------------
    # Create module to test
    #   - do not need to load state_dict
    # --------------------------------------------------------------
    multi_proc_module = SwiGLUWithDTensor()
    multi_proc_module = multi_proc_module.train()
    multi_proc_module = multi_proc_module.to(manager.device)

    # -----------------------------------------------------
    # Create input DTensors
    # ----------------------------------------------------
    placements_for_single_rep_nonparam = (Shard(0), Shard(1), Replicate())

    input_example_as_dtensor = distribute_tensor(
        input_example_device, manager.device_mesh_subgroups, placements_for_single_rep_nonparam
    ).requires_grad_(True)

    output_grad_example_as_dtensor = distribute_tensor(
        output_grad_example,
        manager.device_mesh_subgroups,
        placements_for_single_rep_nonparam,
    ).requires_grad_(False)

    # -------------------------------------------------
    # Run FW
    # -------------------------------------------------
    input_example_clone_as_dtensor = input_example_as_dtensor.detach().clone().requires_grad_(True)
    output_actual_as_dtensor: DTensor = multi_proc_module(input_example_as_dtensor)

    # -------------------------------------------------------
    # V4a: multi-proc FW input tensor values unchanged by FW and BW
    # ------------------------------------------------------
    assert_tensors_identical(
        input_example_clone_as_dtensor.full_tensor(),
        input_example_as_dtensor.full_tensor(),
        check_grad=False,
        check_grad_fn=False,
    )

    # -------------------------------------------------
    # Run BW
    # -------------------------------------------------
    output_grad_example_clone_as_dtensor = (
        output_grad_example_as_dtensor.detach().clone().requires_grad_(output_grad_example_as_dtensor.requires_grad)
    )
    output_actual_as_dtensor.backward(output_grad_example_clone_as_dtensor)

    # -------------------------------------------------------
    # V8: multi-proc FW output tensor values close-to single-proc
    # ------------------------------------------------------
    assert_close(output_actual_as_dtensor.full_tensor(), output_ref_device)

    # -------------------------------------------------------
    # V4b: multi-proc FW input tensor values unchanged by FW and BW
    #   - check again that input is unchanged by BW
    # ------------------------------------------------------
    assert_tensors_identical(
        input_example_clone_as_dtensor.full_tensor(),
        input_example_as_dtensor.full_tensor(),
        check_grad=False,
        check_grad_fn=False,
    )
    # -------------------------------------------------------
    # V5: multi-proc BW input tensor values unchanged by BW
    # ------------------------------------------------------
    assert_tensors_identical(
        output_grad_example_clone_as_dtensor.full_tensor(),
        output_grad_example_as_dtensor.full_tensor(),
        check_grad=False,
        check_grad_fn=False,
    )
    # -------------------------------------------------------
    # V9: multi-proc FW input gradient values close-to single-proc
    # ------------------------------------------------------
    assert_close(
        input_example_as_dtensor.grad.full_tensor(),
        input_grads_ref_device,
    )

    # cleanup
    DistributedManager.cleanup()
    monkeypatch.undo()


def get_example_input_and_reference_output(
    bs: int,
    N_tokens: int,
    dim: int,
    seed: int = SEED,
) -> tuple[Tensor, Tensor, Tensor, Tensor]:
    """Generate example input and reference output for SwiGLU testing.

    Parameters
    ----------
    bs : int
        Batch size.
    N_tokens : int
        Number of tokens.
    dim : int
        Input dimension (must be even for SwiGLU).
    seed : int, optional
        Random seed, by default SEED.

    Returns
    -------
    tuple[Tensor, Tensor, Tensor, Tensor]
        (input_example, output_ref, output_grad_example, input_grads_ref)
    """
    # ----------------------------------------
    # Set random seed
    # ----------------------------------------
    seed_by_rank(seed)

    # ----------------------------------------
    # Create input tensor
    # ----------------------------------------
    input_example = torch.randn(bs, N_tokens, dim, requires_grad=True)
    input_example_copy = input_example.detach().clone().requires_grad_(input_example.requires_grad)

    # -------------------------------------------
    # Create single-proc module and run serial FW
    # -------------------------------------------
    single_proc_module = SwiGluBoltz()
    single_proc_module = single_proc_module.train()
    output_ref = single_proc_module(input_example)

    # ----------------------------------------------------------------
    # V1a: single-proc FW input tensor values unchanged by FW and BW
    # -----------------------------------------------------------------
    assert_tensors_identical(
        input_example,
        input_example_copy,
        check_grad=False,
        check_grad_fn=False,
        check_storage_pointer=False,
    )

    # ----------------------------------------
    # Create output gradient
    # ----------------------------------------
    output_grad_example = torch.randn_like(output_ref, requires_grad=False)
    output_grad_example_copy = output_grad_example.detach().clone().requires_grad_(output_grad_example.requires_grad)

    # ----------------------------------------
    # Serial BW Compute reference gradients
    # ----------------------------------------
    torch.autograd.backward(output_ref, output_grad_example)

    # ----------------------------------------------------------------
    # V1b: single-proc FW input tensor values unchanged by FW and BW
    # -----------------------------------------------------------------
    assert_tensors_identical(
        input_example,
        input_example_copy,
        check_grad=False,
        check_grad_fn=False,
        check_storage_pointer=False,
    )
    # -----------------------------------------------------------------------
    # V2: single-proc BW input tensor values unchanged by BW
    # -----------------------------------------------------------------------
    assert_tensors_identical(
        output_grad_example,
        output_grad_example_copy,
        check_grad=False,
        check_grad_fn=False,
        check_storage_pointer=False,
    )

    # Get input gradients
    input_grads_ref = input_example.grad.clone()

    return input_example, output_ref.detach(), output_grad_example, input_grads_ref


@pytest.mark.parametrize(
    "setup_env",
    itertools.product([(1, (2, 2)), (2, (2, 2))], [True], ["cpu", "cuda"], ["ENV"]),
    indirect=("setup_env",),
    ids=lambda x: f"dp={x[0][0]}, cp={x[0][1]}, specify_method={x[1]}, device_type={x[2]}, method_init={x[3]}",
)
def test_swiglu_with_dtensor(
    setup_env: dict[str, int],
    bs: int = 2,
    N_tokens: int = 4**2,
    dim: int = 2,
    seed: int = SEED,
):
    """Test SwiGLU with DTensor for various configurations.

    Parameters
    ----------
    setup_env : dict[str, int]
        Environment setup for distributed testing.
    bs : int, optional
        Batch size, by default 2.
    N_tokens : int, optional
        Number of tokens, by default 8**2.
    dim : int, optional
        Input dimension, by default 384.
    seed : int, optional
        Random seed, by default SEED.
    """
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env
    skip_if_cuda_not_avail_or_device_count_less_than_word_size(
        device_type=device_type,
        world_size=world_size,
    )

    # (0) Check use-case requirements
    if dim % 2 != 0:
        raise ValueError(f"Dimension must be even for SwiGLU. Got dim={dim}")

    # (1) Get example input and reference output
    input_example, output_ref, output_grad_example, input_grads_ref = get_example_input_and_reference_output(
        bs=bs,
        N_tokens=N_tokens,
        dim=dim,
        seed=seed,
    )

    spawn_multiprocessing(
        assert_swiglu_with_dtensor_fw_bw,
        world_size,
        input_example,
        output_ref,
        output_grad_example,
        input_grads_ref,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
    )


@pytest.mark.parametrize(
    "setup_env",
    itertools.product([(1, (1, 1))], [True], ["cpu", "cuda"], ["ENV"]),
    indirect=("setup_env",),
    ids=lambda x: f"dp={x[0][0]}, cp={x[0][1]}, specify_method={x[1]}, device_type={x[2]}, method_init={x[3]}",
)
def test_swiglu_with_dtensor_for_metadata_checks(
    setup_env: dict[str, int],
    bs: int = 2,
    N_tokens: int = 4**2,
    dim: int = 2,
    seed: int = SEED,
):
    """Test SwiGLU with DTensor for various configurations.

    Parameters
    ----------
    setup_env : dict[str, int]
        Environment setup for distributed testing.
    bs : int, optional
        Batch size, by default 2.
    N_tokens : int, optional
        Number of tokens, by default 8**2.
    dim : int, optional
        Input dimension, by default 384.
    seed : int, optional
        Random seed, by default SEED.
    """
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env
    skip_if_cuda_not_avail_or_device_count_less_than_word_size(
        device_type=device_type,
        world_size=world_size,
    )

    # (0) Check use-case requirements
    if dim % 2 != 0:
        raise ValueError(f"Dimension must be even for SwiGLU. Got dim={dim}")

    # (1) Get example input and reference output
    input_example, output_ref, output_grad_example, input_grads_ref = get_example_input_and_reference_output(
        bs=bs,
        N_tokens=N_tokens,
        dim=dim,
        seed=seed,
    )

    spawn_multiprocessing(
        assert_swiglu_with_dtensor_fw_bw,
        world_size,
        input_example,
        output_ref,
        output_grad_example,
        input_grads_ref,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
    )


if __name__ == "__main__":
    pytest.main([__file__])
