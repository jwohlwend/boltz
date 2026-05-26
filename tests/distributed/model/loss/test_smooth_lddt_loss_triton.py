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


import os
import re
import subprocess
from pathlib import Path

import pytest
import torch
from torch.distributed.tensor import Replicate, Shard, distribute_tensor

from boltz.distributed.comm import TransposeComm
from boltz.distributed.manager import DistributedManager

try:
    from boltz.distributed.model.loss.triton.smooth_lddt_loss import (
        grid_launch_config,
        smooth_lddt_loss_bwd_kernel,
        smooth_lddt_loss_fwd_kernel,
    )

    has_smooth_lddt_loss_triton_kernels = True
except ImportError:
    has_smooth_lddt_loss_triton_kernels = False

from boltz.distributed.model.loss.diffusion import (
    _smooth_lddt_loss_backward_local,
    _smooth_lddt_loss_forward_local,
    _smooth_lddt_loss_local_triton_backward,
    _smooth_lddt_loss_local_triton_forward,
    smooth_lddt_loss_triton,
)
from boltz.distributed.model.modules.utils import PRECISION_TO_DTYPE, Precision, setup_tf32_env
from boltz.model.loss.diffusion import smooth_lddt_loss as smooth_lddt_loss_ref_impl_v1
from boltz.model.loss.diffusionv2 import smooth_lddt_loss as smooth_lddt_loss_ref_impl_v2
from boltz.testing.utils import spawn_multiprocessing


def assert_smooth_lddt_loss_equivalence(rank, payload):
    (
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        dtype,
        multiplicity,
        pred_coords_global,
        true_coords_global,
        is_nucleotide_global,
        coords_mask_global,
        v2,
    ) = payload

    monkeypatch = pytest.MonkeyPatch()
    if env_per_rank is not None:
        for var_name, value in env_per_rank.items():
            if value == "<INPUT_RANK>":
                monkeypatch.setenv(var_name, f"{rank}")
                continue
            monkeypatch.setenv(var_name, value)

    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    manager = DistributedManager()

    # Setup comm
    transpose_comm = TransposeComm(manager.group["cp"], manager.layout_subgroups["cp"])

    # Prepare inputs
    # pred_coords: (Shard(0), Shard(1), Replicate())
    # true_coords: (Shard(0), Shard(1), Replicate())
    # is_nucleotide: (Shard(0), Shard(1), Replicate())
    # coords_mask: (Shard(0), Shard(1), Replicate())

    placements = (Shard(0), Shard(1), Replicate())

    pred_coords_dtensor = distribute_tensor(
        pred_coords_global.to(device=manager.device, dtype=dtype),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements,
    ).requires_grad_(True)

    true_coords_dtensor = distribute_tensor(
        true_coords_global.to(device=manager.device, dtype=dtype),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements,
    ).requires_grad_(False)  # True coords usually don't need grad

    is_nucleotide_dtensor = distribute_tensor(
        is_nucleotide_global.to(device=manager.device, dtype=dtype),  # cast to float/int
        device_mesh=manager.device_mesh_subgroups,
        placements=placements,
    )

    coords_mask_dtensor = distribute_tensor(
        coords_mask_global.to(device=manager.device, dtype=dtype),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements,
    )

    # Extract global tensors for reference implementation (running serially)
    # This avoids potential bugs in distributed reference implementation and compares against
    # the mathematical ground truth (global execution).
    pred_global = pred_coords_dtensor.full_tensor().detach().clone().requires_grad_(True)
    true_global = true_coords_dtensor.full_tensor().detach()
    is_nuc_global = is_nucleotide_dtensor.full_tensor().detach()
    mask_global = coords_mask_dtensor.full_tensor().detach()

    # Run Reference Function (Global/Serial)
    smooth_lddt_loss_ref_impl = smooth_lddt_loss_ref_impl_v2 if v2 else smooth_lddt_loss_ref_impl_v1
    loss_ref = smooth_lddt_loss_ref_impl(
        pred_coords=pred_global,
        true_coords=true_global,
        is_nucleotide=is_nuc_global,
        coords_mask=mask_global,
        multiplicity=multiplicity,
    )

    # Run Custom Function (DTensor version) - Triton
    # We use the original dtensor (v2 copy not needed if we use fresh global for ref)
    # But pred_coords_dtensor tracks grad.
    if dtype == torch.bfloat16:
        # Expect error for bf16
        with pytest.raises(ValueError, match=f"Triton kernel for smooth LDDT loss does not support {dtype}"):
            smooth_lddt_loss_triton(
                pred_coords=pred_coords_dtensor,
                true_coords=true_coords_dtensor,
                is_nucleotide=is_nucleotide_dtensor,
                coords_mask=coords_mask_dtensor,
                comm=transpose_comm,
                multiplicity=multiplicity,
                v2=v2,
            )
        return

    loss_custom = smooth_lddt_loss_triton(
        pred_coords=pred_coords_dtensor,
        true_coords=true_coords_dtensor,
        is_nucleotide=is_nucleotide_dtensor,
        coords_mask=coords_mask_dtensor,
        comm=transpose_comm,
        multiplicity=multiplicity,
        v2=v2,
    )

    # Compare Forward
    # loss_ref is global scalar.
    # loss_custom is DTensor (scalar).
    torch.testing.assert_close(loss_ref, loss_custom.full_tensor())

    # Backward Reference
    loss_ref.backward()
    grad_ref = pred_global.grad

    # Backward Custom
    loss_custom.backward()
    grad_custom = pred_coords_dtensor.grad

    # Compare Backward
    # Gather custom gradients to global to compare with reference
    torch.testing.assert_close(grad_ref, grad_custom.full_tensor())

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.skipif(not torch.cuda.is_available(), reason="CUDA not available")
@pytest.mark.skipif(not has_smooth_lddt_loss_triton_kernels, reason="Triton kernels not available")
@pytest.mark.parametrize("is_self_comm", [True, False])
@pytest.mark.parametrize("dtype_str", ["fp32", "bf16", "tf32", "fp64"])
def test_smooth_lddt_loss_local_triton_equivalence(is_self_comm, dtype_str):
    match dtype_str:
        case "fp32":
            dtype = torch.float32
            precision = Precision.FP32
        case "bf16":
            dtype = torch.bfloat16
            precision = Precision.BF16
        case "tf32":
            dtype = torch.float32
            precision = Precision.TF32
        case "fp64":
            dtype = torch.float64
            precision = Precision.FP64
        case _:
            raise ValueError(f"Unsupported dtype: {dtype_str}")

    with setup_tf32_env(precision):
        # Setup inputs
        B_local = 1
        multiplicity = 16
        N_atom_local = 100
        device = torch.device("cuda")

        init_val_range = 30.0
        pred_coords_local = (
            torch.randn(B_local * multiplicity, N_atom_local, 3, device=device, dtype=dtype) * init_val_range
        )
        true_coords_local = (
            torch.randn(B_local * multiplicity, N_atom_local, 3, device=device, dtype=dtype) * init_val_range
        )
        pred_coords_t_local = (
            torch.randn(B_local * multiplicity, N_atom_local, 3, device=device, dtype=dtype) * init_val_range
        )
        true_coords_t_local = (
            torch.randn(B_local * multiplicity, N_atom_local, 3, device=device, dtype=dtype) * init_val_range
        )

        is_nucleotide_local = torch.randint(0, 2, (B_local, N_atom_local), device=device).bool()
        coords_mask_local = torch.randint(0, 2, (B_local, N_atom_local), device=device).to(dtype=dtype)
        coords_mask_t_local = torch.randint(0, 2, (B_local, N_atom_local), device=device).to(dtype=dtype)

        nucleic_acid_cutoff = 5.0
        other_cutoff = 3.0

        # Clone inputs for Triton forward pass to check for in-place modifications
        pred_coords_local_triton = pred_coords_local.clone()
        true_coords_local_triton = true_coords_local.clone()
        pred_coords_t_local_triton = pred_coords_t_local.clone()
        true_coords_t_local_triton = true_coords_t_local.clone()
        is_nucleotide_local_triton = is_nucleotide_local.clone()
        coords_mask_local_triton = coords_mask_local.clone()
        coords_mask_t_local_triton = coords_mask_t_local.clone()

        # --- Forward Pass ---
        # Run PyTorch version
        num_ref, den_ref = _smooth_lddt_loss_forward_local(
            pred_coords_local,
            true_coords_local,
            pred_coords_t_local,
            true_coords_t_local,
            is_nucleotide_local,
            coords_mask_local,
            coords_mask_t_local,
            is_self_comm,
            nucleic_acid_cutoff,
            other_cutoff,
            multiplicity,
        )

        # Run Triton version
        if dtype == torch.bfloat16:
            with pytest.raises(ValueError, match=f"Triton kernel for smooth LDDT loss does not support {dtype}"):
                _smooth_lddt_loss_local_triton_forward(
                    pred_coords_local_triton,
                    true_coords_local_triton,
                    pred_coords_t_local_triton,
                    true_coords_t_local_triton,
                    is_nucleotide_local_triton,
                    coords_mask_local_triton,
                    coords_mask_t_local_triton,
                    is_self_comm,
                    nucleic_acid_cutoff,
                    other_cutoff,
                    multiplicity,
                )
            return

        num_triton, den_triton = _smooth_lddt_loss_local_triton_forward(
            pred_coords_local_triton,
            true_coords_local_triton,
            pred_coords_t_local_triton,
            true_coords_t_local_triton,
            is_nucleotide_local_triton,
            coords_mask_local_triton,
            coords_mask_t_local_triton,
            is_self_comm,
            nucleic_acid_cutoff,
            other_cutoff,
            multiplicity,
        )

        # Check equivalence
        # PyTorch sum() over bf16 promotes to fp32, while Triton kernel keeps bf16.
        # We cast ref to match triton output for comparison.
        num_ref = num_ref.to(dtype=num_triton.dtype)
        den_ref = den_ref.to(dtype=den_triton.dtype)

        torch.testing.assert_close(num_triton, num_ref)
        torch.testing.assert_close(den_triton, den_ref)

        # Check that inputs were not modified in-place
        torch.testing.assert_close(pred_coords_local_triton, pred_coords_local)
        torch.testing.assert_close(true_coords_local_triton, true_coords_local)
        torch.testing.assert_close(pred_coords_t_local_triton, pred_coords_t_local)
        torch.testing.assert_close(true_coords_t_local_triton, true_coords_t_local)
        torch.testing.assert_close(is_nucleotide_local_triton, is_nucleotide_local)
        torch.testing.assert_close(coords_mask_local_triton, coords_mask_local)
        torch.testing.assert_close(coords_mask_t_local_triton, coords_mask_t_local)

        # --- Backward Pass ---

        # Dummy gradients for num and den (scalars per batch element)
        grad_num_reduced = torch.randn(B_local * multiplicity, device=device, dtype=dtype)
        grad_den_reduced = torch.randn(B_local * multiplicity, device=device, dtype=dtype)

        # Clone gradients for Triton backward pass
        grad_num_reduced_triton = grad_num_reduced.clone()
        grad_den_reduced_triton = grad_den_reduced.clone()

        # Run PyTorch backward
        grad_pred_local_ref, grad_pred_t_local_ref = _smooth_lddt_loss_backward_local(
            grad_num_reduced,
            grad_den_reduced,
            pred_coords_local,
            true_coords_local,
            pred_coords_t_local,
            true_coords_t_local,
            is_nucleotide_local,
            coords_mask_local,
            coords_mask_t_local,
            is_self_comm,
            nucleic_acid_cutoff,
            other_cutoff,
            multiplicity,
        )

        # Run Triton backward
        grad_pred_local_triton, grad_pred_t_local_triton = _smooth_lddt_loss_local_triton_backward(
            grad_num_reduced_triton,
            grad_den_reduced_triton,
            pred_coords_local_triton,
            true_coords_local_triton,
            pred_coords_t_local_triton,
            true_coords_t_local_triton,
            is_nucleotide_local_triton,
            coords_mask_local_triton,
            coords_mask_t_local_triton,
            is_self_comm,
            nucleic_acid_cutoff,
            other_cutoff,
            multiplicity,
        )

        # Check equivalence
        grad_pred_local_ref = grad_pred_local_ref.to(dtype=grad_pred_local_triton.dtype)
        grad_pred_t_local_ref = grad_pred_t_local_ref.to(dtype=grad_pred_t_local_triton.dtype)

        torch.testing.assert_close(grad_pred_local_triton, grad_pred_local_ref, atol=1e-6, rtol=1e-4)
        torch.testing.assert_close(grad_pred_t_local_triton, grad_pred_t_local_ref, atol=1e-6, rtol=1e-4)

        # Check that gradients were not modified in-place
        torch.testing.assert_close(grad_num_reduced_triton, grad_num_reduced)
        torch.testing.assert_close(grad_den_reduced_triton, grad_den_reduced)


@pytest.mark.skipif(not torch.cuda.is_available(), reason="CUDA not available")
@pytest.mark.skipif(not has_smooth_lddt_loss_triton_kernels, reason="Triton kernels not available")
@pytest.mark.parametrize(
    "setup_env",
    (
        params_test := [
            ((2, (2, 2)), True, "cuda", "ENV"),
        ]
    ),
    indirect=("setup_env",),
    ids=[
        f"dp:{x[0][0]}, cp:{x[0][1]}, specify_method:{x[1]}, device_type:{x[2]}, method_init:{x[3]}"
        for x in params_test
    ],
)
@pytest.mark.parametrize("v2", [True, False], ids=["v2", "v1"])
def test_smooth_lddt_loss_equivalence(
    setup_env,
    v2: bool,
    multiplicity: int = 16,
    dtype: torch.dtype = torch.float64,
):
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda" and not torch.cuda.is_available():
        pytest.skip("CUDA not available")

    if device_type == "cuda" and torch.cuda.device_count() < world_size:
        pytest.skip(f"Not enough GPUs. Required: {world_size}, Available: {torch.cuda.device_count()}")

    # Setup dummy data
    B = 1 * grid_group_sizes["dp"]
    N = 1000 * grid_group_sizes["cp"][0]

    B_expanded = B * multiplicity

    init_val_range = 30.0

    pred_coords_global = torch.randn(B_expanded, N, 3) * init_val_range
    true_coords_global = torch.randn(B_expanded, N, 3) * init_val_range
    is_nucleotide_global = torch.randint(0, 2, (B, N)).float()
    coords_mask_global = torch.randint(0, 2, (B, N)).float()

    payload = (
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        dtype,
        multiplicity,
        pred_coords_global,
        true_coords_global,
        is_nucleotide_global,
        coords_mask_global,
        v2,
    )

    spawn_multiprocessing(assert_smooth_lddt_loss_equivalence, world_size, payload)


def assert_no_register_spilling(path_to_ptx_file: Path):
    ptx_code = path_to_ptx_file.read_text()

    # get the ".target sm_{arch}a" directive from the ptx code
    sm_arch_match = re.search(r"\.target (sm_\w+)", ptx_code)
    if not sm_arch_match:
        raise RuntimeError(f"No .target directive found in {path_to_ptx_file}")
    sm_arch = sm_arch_match.group(1)

    # Run ptxas
    # -v: Verbose (prints register/spill stats)
    # --gpu-name=sm_{arch}: Matches target hardware
    ptxas_path = os.environ["TRITON_PTXAS_PATH"]

    cmd = [ptxas_path, "-v", f"--gpu-name={sm_arch}", str(path_to_ptx_file)]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        output = result.stderr
        # Check for spill stores and loads
        # Expected: "0 bytes spill stores, 0 bytes spill loads"
        if "0 bytes spill stores, 0 bytes spill loads" not in output:
            raise RuntimeError(f"Register spilling detected in {path_to_ptx_file}:\n{output}")
        # otherwise, this test will fail with something like this:
        # RuntimeError: Register spilling detected in /tmp/pytest-of-*/**/test_no_register_spilling_fwd0/**/<name_of_kernel>.ptx:
        # ptxas info    : 28 bytes gmem
        # ptxas info    : Compiling entry function 'smooth_lddt_loss_fwd_kernel' for '{sm_arch}'
        # ptxas info    : Function properties for smooth_lddt_loss_fwd_kernel
        #     14976 bytes stack frame, 64736 bytes spill stores, 81760 bytes spill loads
        # ptxas info    : Used 32 registers, used 1 barriers, 14976 bytes cumulative stack size
        # ptxas info    : Compile time = 34792.172 ms
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"ptxas failed with error:\n{e.stderr}") from e


# The spilling test is expensive to run and triton doesn't always follow the recompilation rules
# so we only run the test for a 1 case fwd and bwd
@pytest.mark.skipif(not has_smooth_lddt_loss_triton_kernels, reason="Triton kernels not available")
@pytest.mark.parametrize("precision", [Precision.FP32], ids=lambda x: f"{x}")
@pytest.mark.parametrize("B", [1], ids=lambda x: f"B:{x}")
@pytest.mark.parametrize("M", [16], ids=lambda x: f"M:{x}")
@pytest.mark.parametrize("N", [4608], ids=lambda x: f"N:{x}")
@pytest.mark.parametrize("fwd_or_bwd", ["fwd", "bwd"])
def test_no_register_spilling(tmp_path, monkeypatch, precision, fwd_or_bwd, B, M, N):
    if not torch.cuda.is_available():
        pytest.skip("CUDA not available")

    # setup the env to dump the ptx code
    # NOTE: the cuobjdump version must be recent enough to support the running GPU architecture
    # either wise the cuobjdump call invoked by triton kernel dumping will fail. The user
    # can set, e.g., TRITON_CUOBJDUMP_PATH=$CONDA_PREFIX/bin/cuobjdump or other cuobjdump available
    # that supports the running GPU architecture
    monkeypatch.setenv("TRITON_KERNEL_DUMP", "1")
    monkeypatch.setenv("TRITON_DUMP_DIR", str(tmp_path))
    # Ensure cache dir is unique to avoid hitting cached kernels without dump
    monkeypatch.setenv("TRITON_ALWAYS_COMPILE", "1")
    monkeypatch.setenv("TRITON_CACHE_DIR", str(tmp_path / "cache"))
    # NOTE: the ptxas version must be recent enough to support the running GPU architecture
    # either wise the ptxas call later will fail
    monkeypatch.setenv("TRITON_PTXAS_PATH", os.environ.get("TRITON_PTXAS_PATH", "ptxas"))

    # invoke the kernel to get the ptx code
    D = 3

    device = torch.device("cuda")
    dtype = PRECISION_TO_DTYPE[precision]

    pred_coords_local = torch.randn(B * M, N, D, device=device, dtype=dtype)
    true_coords_local = torch.randn(B * M, N, D, device=device, dtype=dtype)
    pred_coords_t_local = torch.randn(B * M, N, D, device=device, dtype=dtype)
    true_coords_t_local = torch.randn(B * M, N, D, device=device, dtype=dtype)

    is_nucleotide_local = torch.randint(0, 2, (B * M, N), device=device, dtype=torch.bool)
    coords_mask_local = torch.randint(0, 2, (B * M, N), device=device, dtype=dtype)
    coords_mask_t_local = torch.randint(0, 2, (B * M, N), device=device, dtype=dtype)

    num_result = torch.zeros(B * M, device=device, dtype=dtype)
    den_result = torch.zeros(B * M, device=device, dtype=dtype)

    nucleic_acid_cutoff = 5.0
    other_cutoff = 3.0
    is_self_comm = False

    if fwd_or_bwd == "bwd":
        grad_num = torch.randn_like(num_result)
        grad_den = torch.randn_like(den_result)

        grad_pred_coords_local_result = torch.zeros_like(pred_coords_local)
        grad_pred_coords_t_local_result = torch.zeros_like(pred_coords_t_local)
        smooth_lddt_loss_bwd_kernel[grid_launch_config](
            grad_num,
            grad_den,
            pred_coords_local,
            true_coords_local,
            pred_coords_t_local,
            true_coords_t_local,
            is_nucleotide_local,
            coords_mask_local,
            coords_mask_t_local,
            grad_pred_coords_local_result,
            grad_pred_coords_t_local_result,
            pred_coords_local.stride(0),
            pred_coords_local.stride(1),
            pred_coords_local.stride(2),
            true_coords_local.stride(0),
            true_coords_local.stride(1),
            true_coords_local.stride(2),
            pred_coords_t_local.stride(0),
            pred_coords_t_local.stride(1),
            pred_coords_t_local.stride(2),
            true_coords_t_local.stride(0),
            true_coords_t_local.stride(1),
            true_coords_t_local.stride(2),
            is_nucleotide_local.stride(0),
            is_nucleotide_local.stride(1),
            coords_mask_local.stride(0),
            coords_mask_local.stride(1),
            coords_mask_t_local.stride(0),
            coords_mask_t_local.stride(1),
            grad_pred_coords_local_result.stride(0),
            grad_pred_coords_local_result.stride(1),
            grad_pred_coords_local_result.stride(2),
            grad_pred_coords_t_local_result.stride(0),
            grad_pred_coords_t_local_result.stride(1),
            grad_pred_coords_t_local_result.stride(2),
            nucleic_acid_cutoff,
            other_cutoff,
            is_self_comm,
            pred_coords_local.shape[0],
            pred_coords_local.shape[1],
            coords_mask_local.shape[0],
        )
    else:
        smooth_lddt_loss_fwd_kernel[grid_launch_config](
            pred_coords_local,
            true_coords_local,
            pred_coords_t_local,
            true_coords_t_local,
            is_nucleotide_local,
            coords_mask_local,
            coords_mask_t_local,
            num_result,
            den_result,
            pred_coords_local.stride(0),
            pred_coords_local.stride(1),
            pred_coords_local.stride(2),
            true_coords_local.stride(0),
            true_coords_local.stride(1),
            true_coords_local.stride(2),
            pred_coords_t_local.stride(0),
            pred_coords_t_local.stride(1),
            pred_coords_t_local.stride(2),
            true_coords_t_local.stride(0),
            true_coords_t_local.stride(1),
            true_coords_t_local.stride(2),
            is_nucleotide_local.stride(0),
            is_nucleotide_local.stride(1),
            coords_mask_local.stride(0),
            coords_mask_local.stride(1),
            coords_mask_t_local.stride(0),
            coords_mask_t_local.stride(1),
            nucleic_acid_cutoff,
            other_cutoff,
            is_self_comm,
            pred_coords_local.shape[0],
            pred_coords_local.shape[1],
            coords_mask_local.shape[0],
        )

    # parse the ptx code to check for register spilling
    ptx_files = list(tmp_path.glob(f"**/smooth_lddt_loss_{fwd_or_bwd}_kernel.ptx"))

    if not ptx_files:
        raise RuntimeError(f"No PTX file found in {tmp_path}/**/smooth_lddt_loss_{fwd_or_bwd}_kernel.ptx")

    path_to_ptx_file = ptx_files[0]

    assert_no_register_spilling(path_to_ptx_file)
