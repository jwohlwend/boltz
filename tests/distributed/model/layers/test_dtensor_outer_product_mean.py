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
from torch.distributed.tensor import Shard, distribute_tensor

from boltz.distributed.comm import Ring2DComm
from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.layers.outer_product_mean import OuterProductMean as DistributedOuterProductMean
from boltz.model.layers.outer_product_mean import OuterProductMean as SerialOuterProductMean
from boltz.testing.utils import (
    assert_all_identical,
    assert_no_percentile_upshift,
    assert_tensors_identical,
    get_param_by_key,
    init_module_params_uniform,
    init_tensors_uniform,
    seed_by_rank,
    spawn_multiprocessing,
)


def parallel_assert_outer_prod_mean(
    rank,
    grid_group_sizes,
    device_type,
    backend,
    env_per_rank,
    dtype,
    C_in,
    C_hidden,
    C_out,
    layer_state_dict,
    input_global_host,
    mask_global_host,
    output_expected_global_host,
    d_output_expected_global_host,
    d_input_expected_global_host,
    grad_params_expected_global_host,
    output_global_fp32_host: torch.Tensor | None = None,
    d_input_global_fp32_host: torch.Tensor | None = None,
    grad_params_fp32_global_host: dict[str, torch.Tensor] | None = None,
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

    if torch.finfo(dtype).resolution < torch.finfo(output_expected_global_host.dtype).resolution:
        raise ValueError(
            f"Target dtype {dtype} has higher precision than reference output's dtype {output_expected_global_host.dtype}"
        )

    if ((output_global_fp32_host is None) != (d_input_global_fp32_host is None)) or (
        (output_global_fp32_host is not None) != (grad_params_fp32_global_host is not None)
    ):
        raise ValueError(
            "output_global_fp32_host, d_input_global_fp32_host, and grad_params_fp32_global_host must be either all None or all not None"
        )

    check_error_hist = output_global_fp32_host is not None

    layout_map = manager.layout_subgroups["cp"]
    ring_comm = Ring2DComm(manager.group["cp"], manager.subgroups["cp"][0], layout_map)

    module_serial = SerialOuterProductMean(C_in, C_hidden, C_out).to(dtype=dtype)
    module_serial.load_state_dict(layer_state_dict)
    module_serial = module_serial.to(device=manager.device)
    module = DistributedOuterProductMean(module_serial, manager.device_mesh_subgroups, ring_comm)
    module.train()

    placements_input = (Shard(0), Shard(1), Shard(2))
    # Omitting the src_data_rank parameter to distribute_tensor means the data from rank 0 is sharded
    input_dtensor = distribute_tensor(
        input_global_host.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements_input,
    ).requires_grad_(True)
    mask_dtensor = distribute_tensor(
        mask_global_host.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements_input,
    )
    d_output_expected_dtensor = distribute_tensor(
        d_output_expected_global_host.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements_input,
    )
    output_expected_dtensor = distribute_tensor(
        output_expected_global_host.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements_input,
        src_data_rank=None,
    )
    d_input_expected_dtensor = distribute_tensor(
        d_input_expected_global_host.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements_input,
        src_data_rank=None,
    )

    input_dtensor_copy = input_dtensor.detach().clone().requires_grad_(True)
    mask_dtensor_copy = mask_dtensor.detach().clone()

    if check_error_hist:
        output_dtensor_result = module(input_dtensor, mask_dtensor)
        output_dtensor_result.backward(d_output_expected_dtensor)

        output_fp32_dtensor = distribute_tensor(
            output_global_fp32_host.to(manager.device),
            device_mesh=manager.device_mesh_subgroups,
            placements=placements_input,
            src_data_rank=None,
        )

        d_input_fp32_dtensor = distribute_tensor(
            d_input_global_fp32_host.to(manager.device),
            device_mesh=manager.device_mesh_subgroups,
            placements=placements_input,
            src_data_rank=None,
        )

        assert_no_percentile_upshift(
            output_dtensor_result.to_local(),
            output_expected_dtensor.to_local(),
            output_fp32_dtensor.to_local(),
            names_input=("output_cp_fp32", "output_serial_fp64", "output_serial_fp32"),
        )

        assert_no_percentile_upshift(
            input_dtensor.grad.to_local(),
            d_input_expected_dtensor.to_local(),
            d_input_fp32_dtensor.to_local(),
            names_input=("d_input_cp_fp32", "d_input_serial_fp64", "d_input_serial_fp32"),
        )

        for name, grad_param_expected_global in grad_params_expected_global_host.items():
            grad_param_result_global = get_param_by_key(module, name).grad.full_tensor().cpu()
            assert_no_percentile_upshift(
                grad_param_result_global,
                grad_param_expected_global.to(dtype=grad_param_result_global.dtype),
                grad_params_fp32_global_host[name],
                names_input=(f"d_{name}_cp_fp32", f"d_{name}_serial_fp64", f"d_{name}_serial_fp32"),
            )
    else:
        output_dtensor_result = module(input_dtensor, mask_dtensor)

        # no modification on the input
        assert_tensors_identical(
            input_dtensor_copy.to_local(), input_dtensor.to_local(), check_grad=False, check_grad_fn=False
        )
        assert_tensors_identical(mask_dtensor_copy.to_local(), mask_dtensor.to_local())

        # test for consistent forward results with the single-device
        assert (
            output_dtensor_result.shape == output_expected_dtensor.shape
        ), f"Output shape mismatch: {output_dtensor_result.shape} != {output_expected_dtensor.shape}"
        assert (
            output_dtensor_result.stride() == output_expected_dtensor.stride()
        ), f"Output stride mismatch: {output_dtensor_result.stride()} != {output_expected_dtensor.stride()}"
        torch.testing.assert_close(output_dtensor_result.to_local(), output_expected_dtensor.to_local())

        # check backward pass
        d_output_expected_dtensor_copy = d_output_expected_dtensor.detach().clone()
        output_dtensor_result.backward(d_output_expected_dtensor)

        # backward pass should not modify the upstream adjoint
        assert_tensors_identical(d_output_expected_dtensor_copy.to_local(), d_output_expected_dtensor.to_local())

        assert (
            input_dtensor.grad.shape == d_input_expected_dtensor.shape
        ), f"Gradient shape mismatch: {input_dtensor.grad.shape} != {d_input_expected_dtensor.shape}"
        assert (
            input_dtensor.grad.stride() == d_input_expected_dtensor.stride()
        ), f"Gradient stride mismatch: {input_dtensor.grad.stride()} != {d_input_expected_dtensor.stride()}"
        torch.testing.assert_close(input_dtensor.grad.to_local(), d_input_expected_dtensor.to_local())

        # check gradient of the weight
        grad_params_result_dtensors = {}
        for name, param in module.named_parameters():
            if param.grad is not None:
                if name not in grad_params_expected_global_host:
                    # do an extra check here to make sure the parallel computation don't result in extra gradients
                    raise ValueError(f"Parameter {name} has a resulting gradient but it is not in the reference module")
                grad_params_result_dtensors[name] = param.grad

        for name, grad_param_expected_global_host in grad_params_expected_global_host.items():
            assert name in grad_params_result_dtensors, f"Parameter {name}'s gradient is not found in result gradients"
            grad_params_result = grad_params_result_dtensors[name]
            assert (
                grad_params_result.shape == grad_param_expected_global_host.shape
            ), f"Gradient shape mismatch: {grad_params_result.shape} != {grad_param_expected_global_host.shape}"
            assert (
                grad_params_result.stride() == grad_param_expected_global_host.stride()
            ), f"Gradient stride mismatch: {grad_params_result.stride()} != {grad_param_expected_global_host.stride()}"
            grad_params_result_global = grad_params_result.full_tensor()
            torch.testing.assert_close(grad_params_result_global.cpu(), grad_param_expected_global_host.to(dtype=dtype))
            assert_all_identical(grad_params_result_global, manager.group["cp"])

        # check the results with the full tensor to make sure the module's output and
        # and gradients can be gathered into the consistent results with the single-device
        input_global_result = input_dtensor.full_tensor()
        mask_global_result = mask_dtensor.full_tensor()
        output_global_result = output_dtensor_result.full_tensor()
        d_input_global_result = input_dtensor.grad.full_tensor()

        torch.testing.assert_close(input_global_result.cpu(), input_global_host.to(dtype=dtype))
        torch.testing.assert_close(mask_global_result.cpu(), mask_global_host.to(dtype=dtype))
        torch.testing.assert_close(output_global_result.cpu(), output_expected_global_host.to(dtype=dtype))
        torch.testing.assert_close(d_input_global_result.cpu(), d_input_expected_global_host.to(dtype=dtype))

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env, dtype, check_error_hist",
    (
        params_test := [
            ## CUDA tests (2 GPUs)
            # (((2, (1, 1)), True, "cuda", "ENV"), torch.float32, True),
            # (((2, (1, 1)), True, "cuda", "ENV"), torch.float64, True),
            ## CUDA tests (8 GPUs)
            (((2, (2, 2)), True, "cuda", "ENV"), torch.float64, True),
            (((2, (2, 2)), True, "cuda", "ENV"), torch.float32, True),
            (((2, (2, 2)), True, "cuda", "ENV"), torch.float32, False),
            ## CPU tests
            (((1, (3, 3)), True, "cuda", "ENV"), torch.float32, False),
        ]
    ),
    indirect=["setup_env"],
    ids=[
        f"dp:{x[0][0][0]}, cp:{x[0][0][1]}, specify_method:{x[0][1]}, device_type:{x[0][2]}, method_init:{x[0][3]}, "
        f"dtype:{x[1]}, check_error_hist:{x[2]}"
        for x in params_test
    ],
)
def test_outer_product_parallel(setup_env, dtype, check_error_hist):
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    # dtype is the dtype used by the parallel computation
    # check_error_hist determine whether to compare the error histograms between
    # (CP_in_FP32, serial_in_FP64) and (serial_in_FP32, serial_in_FP64)
    # Typically, check_error_hist will use large input dimensions to emulate
    # the real-world use cases. Same with dtype==torch.float64.

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    if check_error_hist:
        if grid_group_sizes["dp"] > 2:
            pytest.skip("skip error histogram check for dp > 1 to save test time")

    # For float64 and error histogram check, we use a realistic model and input size
    # with heavier computation to test the numerical stability. On the other hand,
    # a smaller model and input size incur less numerical error accumulation to allow
    # a larger range of input values to detect logical bugs inexpensively by using
    # smaller dimensions.
    test_large_model = check_error_hist or dtype == torch.float64

    size_ring = grid_group_sizes["cp"][0]

    B = 2 * grid_group_sizes["dp"]
    if test_large_model:
        N = size_ring * 128
        S = size_ring * 128
        C_in = 64
        C_hidden = 32
        C_out = 128
        min_val_init = -5e-2
        max_val_init = 5e-2
    else:
        N = size_ring * 2
        S = size_ring * 3
        C_in = 3
        C_hidden = 5
        C_out = 3
        min_val_init = -0.5
        max_val_init = 0.5

    seed = 42
    seed_by_rank(0, seed=seed)

    # compute reference results with FP64
    input_global_fp64 = torch.empty((B, S, N, C_in), dtype=torch.float64, requires_grad=True, device=device_type)
    mask_global_fp64 = torch.ones((B, S, N), dtype=torch.float64, requires_grad=False, device=device_type)
    mask_global_fp64[0, (S // size_ring) :, :] = 0
    mask_global_fp64[0, :, (N // size_ring) :] = 0
    reference_module = SerialOuterProductMean(C_in, C_hidden, C_out).to(dtype=torch.float64)
    # The output activation and gradient of the layer weights typically increase by 2 to 3 orders of magnitude,
    # where the ULP would be too large and numerical error distribution becomes very wide, i.e., we would have
    # very unpredictable numerical errors. That would make the test results very noisy and not very useful to
    # detect logical bugs in the code. To avoid this, we use a smaller range for the input and layer weights.
    init_tensors_uniform([input_global_fp64], low=min_val_init, high=max_val_init)
    init_module_params_uniform(reference_module, low=min_val_init, high=max_val_init)
    layer_state_dict_fp64 = reference_module.state_dict()
    reference_module = reference_module.to(device=device_type).train()

    output_expected_global_fp64 = reference_module(input_global_fp64, mask_global_fp64)
    d_output_expected_global_fp64 = torch.rand_like(output_expected_global_fp64)
    output_expected_global_fp64.backward(d_output_expected_global_fp64)

    grad_params_fp64_expected_global_host = {
        name: param.grad.detach().clone().cpu() for name, param in reference_module.named_parameters()
    }

    if check_error_hist:
        input_global_fp32 = input_global_fp64.detach().clone().to(dtype=torch.float32).requires_grad_(True)
        mask_global_fp32 = mask_global_fp64.detach().clone().to(dtype=torch.float32).requires_grad_(False)
        reference_module_fp32 = SerialOuterProductMean(C_in, C_hidden, C_out).to(dtype=torch.float32)
        reference_module_fp32.load_state_dict(layer_state_dict_fp64)
        reference_module_fp32 = reference_module_fp32.to(device=device_type).train()
        output_global_fp32 = reference_module_fp32(input_global_fp32, mask_global_fp32)
        d_output_expected_global_fp32 = d_output_expected_global_fp64.to(dtype=torch.float32)
        output_global_fp32.backward(d_output_expected_global_fp32)

        output_global_fp32_host = output_global_fp32.detach().clone().cpu()
        d_input_global_fp32_host = input_global_fp32.grad.detach().clone().cpu()
        grad_params_fp32_global_host = {
            name: param.grad.detach().clone().cpu() for name, param in reference_module_fp32.named_parameters()
        }
    else:
        output_global_fp32_host = None
        d_input_global_fp32_host = None
        grad_params_fp32_global_host = None

    spawn_multiprocessing(
        parallel_assert_outer_prod_mean,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        dtype,
        C_in,
        C_hidden,
        C_out,
        layer_state_dict_fp64,
        input_global_fp64.detach().clone().cpu(),
        mask_global_fp64.detach().clone().cpu(),
        output_expected_global_fp64.detach().clone().cpu(),
        d_output_expected_global_fp64.detach().clone().cpu(),
        input_global_fp64.grad.detach().clone().cpu(),
        grad_params_fp64_expected_global_host,
        output_global_fp32_host,
        d_input_global_fp32_host,
        grad_params_fp32_global_host,
    )
