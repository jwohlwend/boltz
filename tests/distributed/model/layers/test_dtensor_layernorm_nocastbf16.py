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

"""Tests for LayerNormParamsReplicatedNoAutoCastBF16 distributed layer.

This module tests the distributed LayerNormParamsReplicatedNoAutoCastBF16
(from boltz.distributed.model.layers.triangular_attention) against the serial
LayerNorm (from boltz.model.layers.triangular_attention.primitives).
"""

from collections import OrderedDict
from typing import Dict, Optional

import pytest
import torch
from torch.distributed.tensor import Replicate, Shard, distribute_tensor

from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.layers.triangular_attention import LayerNormParamsReplicatedNoAutoCastBF16
from boltz.model.layers.triangular_attention.primitives import LayerNorm as SerialLayerNormNoAutoCastBF16
from boltz.testing.utils import (
    assert_no_percentile_upshift,
    assert_tensors_identical,
    init_module_params_uniform,
    init_tensors_uniform,
    spawn_multiprocessing,
)


def _compute_references(
    B: int, seq_len: int, c_in: int, min_val_init: float, max_val_init: float, device: str = "cpu"
) -> dict:
    """Compute FP64 and FP32 serial references for layernorm.

    Computes on ``device`` for numerical consistency with the DTensor test
    (CUDA tests should pass ``device="cuda"``).  Results are moved to CPU
    for safe transfer across the mp.spawn boundary.
    """
    input_x_fp64 = torch.empty((B, seq_len, c_in), dtype=torch.float64, device=device, requires_grad=True)
    init_tensors_uniform([input_x_fp64], low=min_val_init, high=max_val_init)

    ref_module = SerialLayerNormNoAutoCastBF16(c_in)
    ref_module = ref_module.to(dtype=torch.float64, device=device).train()
    init_module_params_uniform(ref_module, low=min_val_init, high=max_val_init)
    state_dict_fp64 = {k: v.detach().clone().cpu() for k, v in ref_module.state_dict().items()}

    output_fp64 = ref_module(input_x_fp64)
    d_output_fp64 = torch.rand_like(output_fp64)
    output_fp64.backward(d_output_fp64)

    refs = {
        "input_x": input_x_fp64.detach().clone().cpu(),
        "output": output_fp64.detach().clone().cpu(),
        "d_output": d_output_fp64.detach().clone().cpu(),
        "d_input_x": input_x_fp64.grad.detach().clone().cpu(),
        "grad_params": {name: p.grad.detach().clone().cpu() for name, p in ref_module.named_parameters()},
        "state_dict": state_dict_fp64,
    }

    # FP32 serial reference for three-way error histogram comparison
    input_x_fp32 = refs["input_x"].to(dtype=torch.float32, device=device).requires_grad_(True)
    ref_module_fp32 = SerialLayerNormNoAutoCastBF16(c_in)
    ref_module_fp32.load_state_dict(state_dict_fp64)
    ref_module_fp32 = ref_module_fp32.to(dtype=torch.float32, device=device).train()

    output_fp32 = ref_module_fp32(input_x_fp32)
    output_fp32.backward(refs["d_output"].to(dtype=torch.float32, device=device))

    refs["output_fp32"] = output_fp32.detach().clone().cpu()
    refs["d_input_x_fp32"] = input_x_fp32.grad.detach().clone().cpu()
    refs["grad_params_fp32"] = {name: p.grad.detach().clone().cpu() for name, p in ref_module_fp32.named_parameters()}

    return refs


def parallel_assert_dtensor_layernorm_nocastbf16(
    rank: int,
    grid_group_sizes: Dict[str, int],
    device_type: str,
    backend: str,
    env_per_rank: Optional[Dict[str, str]],
    dtype: torch.dtype,
    c_in: int,
    refs_cpu: dict,
    check_error_hist: bool,
):
    """Test distributed LayerNormParamsReplicatedNoAutoCastBF16 in a parallel environment.

    Reference data is computed on CPU in the main process and passed to workers
    via mp.spawn. Workers move tensors to their device before use.
    """
    monkeypatch = pytest.MonkeyPatch()
    if env_per_rank is not None:
        for var_name, value in env_per_rank.items():
            if value == "<INPUT_RANK>":
                monkeypatch.setenv(var_name, f"{rank}")
                continue
            monkeypatch.setenv(var_name, value)

    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    manager = DistributedManager()

    # Move CPU references to worker's device
    refs = {
        "input_x": refs_cpu["input_x"].to(device=manager.device),
        "output": refs_cpu["output"].to(device=manager.device),
        "d_output": refs_cpu["d_output"].to(device=manager.device),
        "d_input_x": refs_cpu["d_input_x"].to(device=manager.device),
        "grad_params": {k: v.to(device=manager.device) for k, v in refs_cpu["grad_params"].items()},
        "state_dict": {k: v.to(device=manager.device) for k, v in refs_cpu["state_dict"].items()},
        "output_fp32": refs_cpu["output_fp32"].to(device=manager.device),
        "d_input_x_fp32": refs_cpu["d_input_x_fp32"].to(device=manager.device),
        "grad_params_fp32": {k: v.to(device=manager.device) for k, v in refs_cpu["grad_params_fp32"].items()},
    }

    if torch.finfo(dtype).resolution < torch.finfo(refs["output"].dtype).resolution:
        raise ValueError(
            f"Target dtype {dtype} has higher precision than reference output's dtype {refs['output'].dtype}"
        )

    # --- build distributed module ---
    module_serial = SerialLayerNormNoAutoCastBF16(c_in)
    module_serial = module_serial.to(dtype=dtype, device=manager.device)
    module_serial.load_state_dict(refs["state_dict"])

    module = LayerNormParamsReplicatedNoAutoCastBF16(module_serial, manager.device_mesh_subgroups)
    module = module.to(device=manager.device).train()

    # Input: (B, seq_len, c_in) — shard on batch (dim 0) and seq (dim 1), replicate feature dim
    placements_input = (Shard(0), Shard(1), Replicate())

    # distribute_tensor with default src_data_rank=0 keeps NCCL streams
    # synchronised with the broadcasts that happened inside the module ctor.
    input_x_dtensor = distribute_tensor(
        refs["input_x"].to(dtype=dtype),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements_input,
    ).requires_grad_(True)

    d_output_dtensor = distribute_tensor(
        refs["d_output"].to(dtype=dtype),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements_input,
    )

    # Forward and backward pass
    output_dtensor_result = module(input_x_dtensor)
    output_dtensor_result.backward(d_output_dtensor)

    # Use full_tensor() to get fully reduced parameter gradients. This works
    # correctly whether the backward returns Replicate placements (eager
    # reduction, no-op) or Partial placements (triggers all-reduce).
    grad_params_global: dict[str, torch.Tensor] = {}
    for name, param in module.named_parameters():
        if param.grad is not None:
            grad_params_global[name] = param.grad.full_tensor()

    if check_error_hist:
        output_expected_dtensor = distribute_tensor(
            refs["output"].to(dtype=dtype),
            device_mesh=manager.device_mesh_subgroups,
            placements=placements_input,
        )
        d_input_x_expected_dtensor = distribute_tensor(
            refs["d_input_x"].to(dtype=dtype),
            device_mesh=manager.device_mesh_subgroups,
            placements=placements_input,
        )
        output_fp32_dtensor = distribute_tensor(
            refs["output_fp32"],
            device_mesh=manager.device_mesh_subgroups,
            placements=placements_input,
        )
        d_input_x_fp32_dtensor = distribute_tensor(
            refs["d_input_x_fp32"],
            device_mesh=manager.device_mesh_subgroups,
            placements=placements_input,
        )

        # --- All collectives done. Only local assertions below. ---

        assert_no_percentile_upshift(
            output_dtensor_result.to_local(),
            output_expected_dtensor.to_local(),
            output_fp32_dtensor.to_local(),
            names_input=("output_cp_fp32", "output_serial_fp64", "output_serial_fp32"),
        )

        assert_no_percentile_upshift(
            input_x_dtensor.grad.to_local(),
            d_input_x_expected_dtensor.to_local(),
            d_input_x_fp32_dtensor.to_local(),
            names_input=("d_input_x_cp_fp32", "d_input_x_serial_fp64", "d_input_x_serial_fp32"),
        )

        # Parameter gradient tolerance: the distributed sum splits the
        # reduction across DP ranks (different FP32 accumulation order).
        # For c_in=8 the bias gradient has only 8 elements; each is a sum
        # of B*seq_len/dp values (~256).  The FP32 ULP at magnitude 256
        # is 256·2^{-23} ≈ 3e-5, so a 1-ULP accumulation-order shift
        # produces ~1.5e-5 absolute error — above the default atol=1e-5.
        # Use 5e-5 (≈2 ULP at magnitude 256) for comfortable first-
        # principles margin.
        perc_param_grad = OrderedDict({0.25: (5e-5, 1e-4), 0.5: (5e-5, 1e-4), 0.75: (5e-5, 1e-4), 0.95: (5e-5, 1e-4)})

        for name, grad_expected in refs["grad_params"].items():
            if name not in grad_params_global:
                raise ValueError(f"Parameter {name}'s gradient is not found in the distributed module")

            assert_no_percentile_upshift(
                grad_params_global[name],
                grad_expected.to(dtype=grad_params_global[name].dtype),
                refs["grad_params_fp32"][name],
                perc=perc_param_grad,
                names_input=(f"d_{name}_cp_fp32", f"d_{name}_serial_fp64", f"d_{name}_serial_fp32"),
            )
    else:
        output_expected_dtensor = distribute_tensor(
            refs["output"].to(dtype=dtype),
            device_mesh=manager.device_mesh_subgroups,
            placements=placements_input,
        )
        d_input_x_expected_dtensor = distribute_tensor(
            refs["d_input_x"].to(dtype=dtype),
            device_mesh=manager.device_mesh_subgroups,
            placements=placements_input,
        )

        # Gather full tensors (collectives) before any assertions.
        output_global_result = output_dtensor_result.full_tensor().cpu()
        d_input_x_global_result = input_x_dtensor.grad.full_tensor().cpu()

        # all_gather for assert_all_identical — do the collective now,
        # assert on the gathered data later.
        gathered_param_grads: dict[str, list[torch.Tensor]] = {}
        cp_group = manager.group["cp"]
        cp_world_size = torch.distributed.get_world_size(cp_group)
        for name in grad_params_global:
            grad_on_device = grad_params_global[name].to(device=manager.device)
            tensor_list = [torch.empty_like(grad_on_device) for _ in range(cp_world_size)]
            torch.distributed.all_gather(tensor_list, grad_on_device, group=cp_group)
            gathered_param_grads[name] = tensor_list

        # --- All collectives done. Only local assertions below. ---

        assert output_dtensor_result.shape == output_expected_dtensor.shape
        assert output_dtensor_result.stride() == output_expected_dtensor.stride()
        torch.testing.assert_close(output_dtensor_result.to_local(), output_expected_dtensor.to_local())

        assert input_x_dtensor.grad.shape == d_input_x_expected_dtensor.shape
        assert input_x_dtensor.grad.stride() == d_input_x_expected_dtensor.stride()
        torch.testing.assert_close(input_x_dtensor.grad.to_local(), d_input_x_expected_dtensor.to_local())

        torch.testing.assert_close(output_global_result, refs["output"].to(dtype=dtype).cpu())
        torch.testing.assert_close(d_input_x_global_result, refs["d_input_x"].to(dtype=dtype).cpu())

        for name, param in module.named_parameters():
            if param.grad is not None:
                if name not in refs["grad_params"]:
                    raise ValueError(f"Parameter {name} has a gradient but is not in the reference")
                torch.testing.assert_close(grad_params_global[name], refs["grad_params"][name].to(dtype=dtype))
                grad_on_device = grad_params_global[name].to(device=manager.device)
                for gathered in gathered_param_grads[name]:
                    assert_tensors_identical(gathered, grad_on_device)

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env, dtype, check_error_hist",
    (
        params_test := [
            ## CUDA tests (2 GPUs)
            (((2, (1, 1)), True, "cuda", "ENV"), torch.float32, True),
            (((2, (1, 1)), True, "cuda", "ENV"), torch.float64, True),
            ## CUDA tests (4 GPUs)
            (((1, (2, 2)), True, "cuda", "ENV"), torch.float32, True),
            (((1, (2, 2)), True, "cuda", "ENV"), torch.float32, True),
            ## CUDA tests (8 GPUs)
            (((2, (2, 2)), True, "cuda", "ENV"), torch.float64, True),
            (((2, (2, 2)), True, "cuda", "ENV"), torch.float32, True),
            ## CPU tests
            (((2, (3, 3)), True, "cpu", "ENV"), torch.float32, True),
        ]
    ),
    indirect=["setup_env"],
    ids=[
        f"dp:{x[0][0][0]}, cp:{x[0][0][1]}, specify_method:{x[0][1]}, device_type:{x[0][2]}, method_init:{x[0][3]}, "
        f"dtype:{x[1]}, check_error_hist:{x[2]}"
        for x in params_test
    ],
)
@pytest.mark.parametrize("c_in", [8, 128])
def test_dtensor_layernorm_nocastbf16(
    setup_env: tuple[dict, int, str, str, str, dict[str, str]],
    dtype: torch.dtype,
    check_error_hist: bool,
    c_in: int,
):
    """Test distributed LayerNormParamsReplicatedNoAutoCastBF16 across multiple processes.

    When check_error_hist=True, uses the three-way error histogram comparison:
    CP FP32 vs serial FP64 ref, compared against serial FP32 vs serial FP64 ref.
    Uses default tolerances from assert_no_percentile_upshift (same as triangle attention test).

    When check_error_hist=False, uses exact match against the FP64 reference.
    """
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    if check_error_hist:
        if grid_group_sizes["dp"] > 2:
            pytest.skip("skip error histogram check for dp > 1 to save test time")

    # Use larger dimensions for error histogram check to emulate realistic workloads
    test_large_model = check_error_hist or dtype == torch.float64

    size_ring = grid_group_sizes["cp"][0]
    B = 2 * grid_group_sizes["dp"]
    if test_large_model:
        seq_len = size_ring * 128
        min_val_init = -5e-2
        max_val_init = 5e-2
    else:
        seq_len = size_ring * 4
        min_val_init = -0.5
        max_val_init = 0.5

    # Compute serial references on the same backend as the DTensor test for
    # numerical consistency.  Results are moved to CPU for mp.spawn transfer.
    torch.manual_seed(42)
    refs_cpu = _compute_references(B, seq_len, c_in, min_val_init, max_val_init, device=device_type)

    spawn_multiprocessing(
        parallel_assert_dtensor_layernorm_nocastbf16,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        dtype,
        c_in,
        refs_cpu,
        check_error_hist,
    )
