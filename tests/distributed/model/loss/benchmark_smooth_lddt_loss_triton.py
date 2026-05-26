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


import argparse
import time
from contextlib import contextmanager

import torch

from boltz.distributed.model.loss.diffusion import (
    _smooth_lddt_loss_backward_local,
    _smooth_lddt_loss_forward_local,
    _smooth_lddt_loss_local_triton_backward,
    _smooth_lddt_loss_local_triton_forward,
)
from boltz.distributed.model.modules.utils import Precision, setup_tf32_env

try:
    import triton  # noqa: F401

    HAS_TRITON = True
except ImportError:
    HAS_TRITON = False


@contextmanager
def benchmark_peak_memory_and_runtime(num_warmup=5, num_iter=10):
    """Benchmark with proper warmup and averaging.

    Args:
        num_warmup: Number of warmup iterations (for kernel compilation/autotuning)
        num_iter: Number of timed iterations for averaging
    """
    torch.cuda.empty_cache()
    torch.cuda.reset_peak_memory_stats()
    start_mem = torch.cuda.memory_allocated()

    stats = {"warmup_iters": num_warmup, "timed_iters": num_iter}
    yield stats

    torch.cuda.synchronize()
    peak_mem = torch.cuda.max_memory_allocated()

    # We want the peak memory *induced* by the function, relative to start.
    peak_usage_mb = (peak_mem - start_mem) / 1024 / 1024
    stats["peak_mem"] = peak_usage_mb


def run_benchmark(B_local, N_atoms, multiplicity=1, dtype_str="fp32"):
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
        case _:
            raise ValueError(f"Unsupported dtype: {dtype_str}")

    print(f"\nBenchmark Config: B_local={B_local}, N_atoms={N_atoms}, Multiplicity={multiplicity}, Dtype={dtype_str}")
    device = torch.device("cuda")

    # Calculate grid sizes for diagnostic
    effective_batch = B_local * multiplicity
    fwd_block = 64  # Default from kernel
    bwd_block = 32  # Default from kernel
    fwd_grid_size = (effective_batch, (N_atoms + fwd_block - 1) // fwd_block, (N_atoms + fwd_block - 1) // fwd_block)
    bwd_grid_size = (effective_batch, (N_atoms + bwd_block - 1) // bwd_block, (N_atoms + bwd_block - 1) // bwd_block)
    fwd_total_programs = fwd_grid_size[0] * fwd_grid_size[1] * fwd_grid_size[2]
    bwd_total_programs = bwd_grid_size[0] * bwd_grid_size[1] * bwd_grid_size[2]

    print(f"Effective batch size: {effective_batch}")
    print(f"Forward grid (BLOCK={fwd_block}): {fwd_grid_size} = {fwd_total_programs:,} programs")
    print(f"Backward grid (BLOCK={bwd_block}): {bwd_grid_size} = {bwd_total_programs:,} programs")

    if fwd_total_programs > 100000:
        print("WARNING: Very large grid size may cause atomic contention!")

    print("NOTE: First Triton run includes kernel compilation + autotuning overhead")

    # Setup inputs
    pred_coords_local = torch.randn(B_local * multiplicity, N_atoms, 3, device=device, dtype=dtype)
    true_coords_local = torch.randn(B_local * multiplicity, N_atoms, 3, device=device, dtype=dtype)
    pred_coords_t_local = torch.randn(B_local * multiplicity, N_atoms, 3, device=device, dtype=dtype)
    true_coords_t_local = torch.randn(B_local * multiplicity, N_atoms, 3, device=device, dtype=dtype)

    # Expand inside function, but inputs here match function signature expected input
    # Function expects (B, M) for these, and expands them inside.
    is_nucleotide_local = torch.randint(0, 2, (B_local, N_atoms), device=device, dtype=torch.bool)
    coords_mask_local = torch.randint(0, 2, (B_local, N_atoms), device=device, dtype=dtype)
    coords_mask_t_local = torch.randint(0, 2, (B_local, N_atoms), device=device, dtype=dtype)

    nucleic_acid_cutoff = 5.0
    other_cutoff = 3.0
    is_self_comm = False  # Simplifies mask logic

    print("-" * 60)
    print("FORWARD PASS")

    with setup_tf32_env(precision):
        # PyTorch Forward
        def run_pytorch_fwd():
            return _smooth_lddt_loss_forward_local(
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

        with benchmark_peak_memory_and_runtime(num_warmup=3, num_iter=10) as stats:
            # Warmup
            for _ in range(stats["warmup_iters"]):
                run_pytorch_fwd()
            torch.cuda.synchronize()

            # Timed runs
            start_time = time.time()
            for _ in range(stats["timed_iters"]):
                run_pytorch_fwd()
            torch.cuda.synchronize()
            end_time = time.time()

            stats["time"] = ((end_time - start_time) / stats["timed_iters"]) * 1000

        print(
            f"PyTorch Forward: Peak Memory = {stats['peak_mem']:.2f} MB, Time = {stats['time']:.2f} ms "
            f"(avg of {stats['timed_iters']} runs)"
        )

        if HAS_TRITON:
            # Triton Forward
            def run_triton_fwd():
                return _smooth_lddt_loss_local_triton_forward(
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

            print("  (First run may be slow due to kernel compilation and autotuning...)")
            with benchmark_peak_memory_and_runtime(num_warmup=10, num_iter=10) as stats:
                # Warmup (includes autotuning on first run)
                for i in range(stats["warmup_iters"]):
                    if i == 0:
                        print("  Warming up (compiling/autotuning)...", end="", flush=True)
                    run_triton_fwd()
                    if i == 0:
                        print(" done")
                torch.cuda.synchronize()

                # Timed runs
                start_time = time.time()
                for _ in range(stats["timed_iters"]):
                    run_triton_fwd()
                torch.cuda.synchronize()
                end_time = time.time()

                stats["time"] = ((end_time - start_time) / stats["timed_iters"]) * 1000

            print(
                f"Triton Forward : Peak Memory = {stats['peak_mem']:.2f} MB, Time = {stats['time']:.2f} ms "
                f"(avg of {stats['timed_iters']} runs)"
            )
        else:
            print("Triton not available, skipping Triton Forward")

        print("-" * 60)
        print("BACKWARD PASS")

        # Inputs for backward
        grad_num_reduced = torch.randn(B_local * multiplicity, device=device, dtype=dtype)
        grad_den_reduced = torch.randn(B_local * multiplicity, device=device, dtype=dtype)

        # PyTorch Backward
        def run_pytorch_bwd():
            return _smooth_lddt_loss_backward_local(
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

        with benchmark_peak_memory_and_runtime(num_warmup=3, num_iter=10) as stats:
            # Warmup
            for _ in range(stats["warmup_iters"]):
                run_pytorch_bwd()
            torch.cuda.synchronize()

            # Timed runs
            start_time = time.time()
            for _ in range(stats["timed_iters"]):
                run_pytorch_bwd()
            torch.cuda.synchronize()
            end_time = time.time()

            stats["time"] = ((end_time - start_time) / stats["timed_iters"]) * 1000

        print(
            f"PyTorch Backward: Peak Memory = {stats['peak_mem']:.2f} MB, Time = {stats['time']:.2f} ms "
            f"(avg of {stats['timed_iters']} runs)"
        )

        if HAS_TRITON:
            # Triton Backward
            def run_triton_bwd():
                return _smooth_lddt_loss_local_triton_backward(
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

            print("  (First run may be slow due to kernel compilation and autotuning...)")
            with benchmark_peak_memory_and_runtime(num_warmup=10, num_iter=10) as stats:
                # Warmup (includes autotuning on first run)
                for i in range(stats["warmup_iters"]):
                    if i == 0:
                        print("  Warming up (compiling/autotuning)...", end="", flush=True)
                    run_triton_bwd()
                    if i == 0:
                        print(" done")
                torch.cuda.synchronize()

                # Timed runs
                start_time = time.time()
                for _ in range(stats["timed_iters"]):
                    run_triton_bwd()
                torch.cuda.synchronize()
                end_time = time.time()

                stats["time"] = ((end_time - start_time) / stats["timed_iters"]) * 1000

            print(
                f"Triton Backward : Peak Memory = {stats['peak_mem']:.2f} MB, Time = {stats['time']:.2f} ms "
                f"(avg of {stats['timed_iters']} runs)"
            )
        else:
            print("Triton not available, skipping Triton Backward")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--B_local", type=int, default=1)
    parser.add_argument("--N_atoms", type=int, default=9 * 512)
    parser.add_argument("--multiplicity", type=int, default=16)
    parser.add_argument("--dtype", type=str, default="fp32", choices=["fp32", "bf16", "tf32"])

    args = parser.parse_args()

    if not torch.cuda.is_available():
        print("CUDA not available, exiting")
        exit(0)

    run_benchmark(args.B_local, args.N_atoms, args.multiplicity, args.dtype)
