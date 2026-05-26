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
from boltz.distributed.model.layers.triangular_mult import (
    TriangleMultiplicationIncoming as DistributedTriangleMultiplicationIncoming,
)
from boltz.distributed.model.layers.triangular_mult import (
    TriangleMultiplicationOutgoing as DistributedTriangleMultiplicationOutgoing,
)
from boltz.distributed.model.layers.triangular_mult import _Direction
from boltz.model.layers.triangular_mult import TriangleMultiplicationIncoming, TriangleMultiplicationOutgoing
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


def parallel_assert_triangle_multiplication(
    rank,
    grid_group_sizes,
    device_type,
    backend,
    env_per_rank,
    dtype,
    dim,
    direction,
    layer_state_dict,
    input_x_global_host,
    mask_global_host,
    output_expected_global_host,
    d_output_expected_global_host,
    d_input_x_expected_global_host,
    grad_params_expected_global_host,
    output_global_fp32_host: torch.Tensor | None = None,
    d_input_x_global_fp32_host: torch.Tensor | None = None,
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

    if ((output_global_fp32_host is None) != (d_input_x_global_fp32_host is None)) or (
        (output_global_fp32_host is not None) != (grad_params_fp32_global_host is not None)
    ):
        raise ValueError(
            "output_global_fp32_host, d_input_x_global_fp32_host, and grad_params_fp32_global_host must be either all None or all not None"
        )

    check_error_hist = output_global_fp32_host is not None

    layout_map = manager.layout_subgroups["cp"]
    ring_comm = Ring2DComm(manager.group["cp"], manager.subgroups["cp"][0], layout_map)

    if direction == _Direction.Outgoing:
        module_serial = TriangleMultiplicationOutgoing(dim)
    elif direction == _Direction.Incoming:
        module_serial = TriangleMultiplicationIncoming(dim)
    else:
        raise ValueError(f"Invalid direction {direction}")
    module_serial.load_state_dict(layer_state_dict)
    module_serial = module_serial.to(dtype=dtype, device=manager.device)

    if direction == _Direction.Outgoing:
        module = DistributedTriangleMultiplicationOutgoing(module_serial, manager.device_mesh_subgroups, ring_comm)
    elif direction == _Direction.Incoming:
        module = DistributedTriangleMultiplicationIncoming(module_serial, manager.device_mesh_subgroups, ring_comm)
    else:
        raise ValueError(f"Invalid direction {direction}")
    module = module.train()

    # Input tensors have the same sharding pattern:
    # x: (B, N, N, D) - sharded on dims 1 and 2 (N and N)
    # mask: (B, N, N) - sharded on dims 1 and 2 (N and N)
    placements = (Shard(0), Shard(1), Shard(2))

    # Distribute input tensors
    input_x_dtensor = distribute_tensor(
        input_x_global_host.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements,
    ).requires_grad_(True)

    mask_dtensor = distribute_tensor(
        mask_global_host.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements,
    )

    # Distribute expected outputs
    d_output_expected_dtensor = distribute_tensor(
        d_output_expected_global_host.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements,
    )
    output_expected_dtensor = distribute_tensor(
        output_expected_global_host.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements,
        src_data_rank=None,
    )
    d_input_x_expected_dtensor = distribute_tensor(
        d_input_x_expected_global_host.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements,
        src_data_rank=None,
    )

    # Create copies to verify inputs aren't modified
    input_x_dtensor_copy = input_x_dtensor.detach().clone().requires_grad_(True)
    mask_dtensor_copy = mask_dtensor.detach().clone()

    if check_error_hist:
        # Forward and backward pass for error histogram checking
        output_dtensor_result = module(input_x_dtensor, mask_dtensor)
        output_dtensor_result.backward(d_output_expected_dtensor)

        output_fp32_dtensor = distribute_tensor(
            output_global_fp32_host.to(device=manager.device),
            device_mesh=manager.device_mesh_subgroups,
            placements=placements,
            src_data_rank=None,
        )

        d_input_x_fp32_dtensor = distribute_tensor(
            d_input_x_global_fp32_host.to(device=manager.device),
            device_mesh=manager.device_mesh_subgroups,
            placements=placements,
            src_data_rank=None,
        )

        # Check that the output tensor has the correct shape
        assert (
            output_dtensor_result.shape == output_expected_dtensor.shape
        ), f"Output DTensor has shape {output_dtensor_result.shape} but expected shape {output_expected_dtensor.shape}"

        # Check that the output tensor has the correct shape
        assert (
            output_dtensor_result.stride() == output_expected_dtensor.stride()
        ), f"Output DTensor has stride {output_dtensor_result.stride()} but expected stride {output_expected_dtensor.stride()}"

        assert (
            input_x_dtensor.grad.shape == d_input_x_expected_dtensor.shape
        ), f"Input DTensor grad has shape {input_x_dtensor.grad.shape} but expected shape {d_input_x_expected_dtensor.shape}"

        assert (
            input_x_dtensor.grad.stride() == d_input_x_expected_dtensor.stride()
        ), f"Input DTensor grad has stride {input_x_dtensor.grad.stride()} but expected stride {d_input_x_expected_dtensor.stride()}"

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

        # Check parameter gradients error histograms
        for name, grad_param_expected_global in grad_params_expected_global_host.items():
            grad_param_result_global = get_param_by_key(module, name).grad.full_tensor().cpu()
            assert_no_percentile_upshift(
                grad_param_result_global,
                grad_param_expected_global.to(dtype=grad_param_result_global.dtype),
                grad_params_fp32_global_host[name],
                names_input=(f"d_{name}_cp_fp32", f"d_{name}_serial_fp64", f"d_{name}_serial_fp32"),
            )
    else:
        # Forward pass
        output_dtensor_result = module(input_x_dtensor, mask_dtensor)

        # Check that the output tensor has the correct shape
        assert (
            output_dtensor_result.shape == output_expected_dtensor.shape
        ), f"Output DTensor has shape {output_dtensor_result.shape} but expected shape {output_expected_dtensor.shape}"

        # Check that the output tensor has the correct shape
        assert (
            output_dtensor_result.stride() == output_expected_dtensor.stride()
        ), f"Output DTensor has stride {output_dtensor_result.stride()} but expected stride {output_expected_dtensor.stride()}"

        # Verify inputs weren't modified
        assert_tensors_identical(
            input_x_dtensor_copy.to_local(), input_x_dtensor.to_local(), check_grad=False, check_grad_fn=False
        )
        assert_tensors_identical(mask_dtensor_copy.to_local(), mask_dtensor.to_local())

        # Test forward pass results
        torch.testing.assert_close(output_dtensor_result.to_local(), output_expected_dtensor.to_local())

        # Backward pass
        d_output_expected_dtensor_copy = d_output_expected_dtensor.detach().clone()
        output_dtensor_result.backward(d_output_expected_dtensor)

        assert (
            input_x_dtensor.grad.shape == d_input_x_expected_dtensor.shape
        ), f"Input DTensor grad has shape {input_x_dtensor.grad.shape} but expected shape {d_input_x_expected_dtensor.shape}"

        assert (
            input_x_dtensor.grad.stride() == d_input_x_expected_dtensor.stride()
        ), f"Input DTensor grad has stride {input_x_dtensor.grad.stride()} but expected stride {d_input_x_expected_dtensor.stride()}"

        # Verify upstream gradient wasn't modified
        assert_tensors_identical(d_output_expected_dtensor_copy.to_local(), d_output_expected_dtensor.to_local())

        # Test input gradients
        torch.testing.assert_close(input_x_dtensor.grad.to_local(), d_input_x_expected_dtensor.to_local())

        # Test full tensor gathering - verify distributed results match serial results
        output_global_result_host = output_dtensor_result.full_tensor().cpu()
        d_input_x_global_result_host = input_x_dtensor.grad.full_tensor().cpu()

        # Verify full tensors match expected results
        torch.testing.assert_close(output_global_result_host, output_expected_global_host.to(dtype=dtype))
        torch.testing.assert_close(d_input_x_global_result_host, d_input_x_expected_global_host.to(dtype=dtype))

        # Test parameter gradients
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
            grad_params_result_global = grad_params_result.full_tensor()
            torch.testing.assert_close(grad_params_result_global.cpu(), grad_param_expected_global_host.to(dtype=dtype))
            assert_all_identical(grad_params_result_global, manager.group["cp"])

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env, dtype, check_error_hist",
    (
        params_test := [
            (((1, (2, 2)), True, "cuda", "ENV"), torch.float32, True),
            (((1, (2, 2)), True, "cuda", "ENV"), torch.float64, False),
            (((2, (2, 2)), True, "cuda", "ENV"), torch.float32, True),
            (((1, (3, 3)), True, "cuda", "ENV"), torch.float32, False),
            (((1, (3, 3)), True, "cpu", "ENV"), torch.float32, False),
        ]
    ),
    indirect=["setup_env"],
    ids=[
        f"dp:{x[0][0][0]}, cp:{x[0][0][1]}, specify_method:{x[0][1]}, device_type:{x[0][2]}, method_init:{x[0][3]}, "
        f"dtype:{x[1]}, check_error_hist:{x[2]}"
        for x in params_test
    ],
)
@pytest.mark.parametrize("direction", [_Direction.Outgoing, _Direction.Incoming])
def test_triangle_multiplication_parallel(setup_env, dtype, check_error_hist, direction):
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
        if grid_group_sizes["dp"] > 1:
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
        N = size_ring * 128  # Number of tokens
        dim = 128  # Hidden dimension
        min_val_init = -5e-2 if dtype == torch.float64 else -1e-3
        max_val_init = -min_val_init
    else:
        N = size_ring * 4  # Number of tokens
        dim = 8  # Hidden dimension
        min_val_init = -0.5
        max_val_init = 0.5

    seed = 42
    seed_by_rank(0, seed=seed)

    # compute reference results with FP64
    input_x_global_fp64 = torch.empty((B, N, N, dim), dtype=torch.float64, requires_grad=True, device=device_type)
    mask_global_fp64 = torch.randint(0, 2, (B, N, N), dtype=torch.float64, requires_grad=False, device=device_type)

    # emulate blocks of pure padding
    mask_global_fp64[0, N // size_ring :, :] = 0
    mask_global_fp64[0, :, N // size_ring :] = 0

    # Create reference serial module
    if direction == _Direction.Outgoing:
        reference_module = TriangleMultiplicationOutgoing(dim)
    elif direction == _Direction.Incoming:
        reference_module = TriangleMultiplicationIncoming(dim)
    else:
        raise ValueError(f"Invalid direction {direction}")

    # The output activation and gradient of the layer weights typically increase by 2 to 3 orders of magnitude,
    # where the ULP would be too large and numerical error distribution becomes very wide, i.e., we would have
    # very unpredictable numerical errors. That would make the test results very noisy and not very useful to
    # detect logical bugs in the code. To avoid this, we use a smaller range for the input and layer weights.
    init_tensors_uniform([input_x_global_fp64], low=min_val_init, high=max_val_init)
    init_module_params_uniform(reference_module, low=min_val_init, high=max_val_init)

    layer_state_dict_fp64 = reference_module.state_dict()
    reference_module = reference_module.to(dtype=torch.float64, device=device_type).train()

    # Run forward pass
    output_expected_global_fp64 = reference_module(input_x_global_fp64, mask_global_fp64)
    d_output_expected_global_fp64 = torch.rand_like(output_expected_global_fp64)
    output_expected_global_fp64.backward(d_output_expected_global_fp64)

    grad_params_fp64_expected_global_host = {
        name: param.grad.detach().clone().cpu() for name, param in reference_module.named_parameters()
    }

    if check_error_hist:
        input_x_global_fp32 = input_x_global_fp64.detach().clone().to(dtype=torch.float32).requires_grad_(True)
        mask_global_fp32 = mask_global_fp64.detach().clone().to(dtype=torch.float32).requires_grad_(False)

        if direction == _Direction.Outgoing:
            reference_module_fp32 = TriangleMultiplicationOutgoing(dim)
        elif direction == _Direction.Incoming:
            reference_module_fp32 = TriangleMultiplicationIncoming(dim)
        else:
            raise ValueError(f"Invalid direction {direction}")

        reference_module_fp32.load_state_dict(layer_state_dict_fp64)
        reference_module_fp32 = reference_module_fp32.to(dtype=torch.float32, device=device_type).train()

        output_global_fp32 = reference_module_fp32(input_x_global_fp32, mask_global_fp32)
        d_output_expected_global_fp32 = d_output_expected_global_fp64.to(dtype=torch.float32)
        output_global_fp32.backward(d_output_expected_global_fp32)

        output_global_fp32_host = output_global_fp32.detach().clone().cpu()
        d_input_x_global_fp32_host = input_x_global_fp32.grad.detach().clone().cpu()
        grad_params_fp32_global_host = {
            name: param.grad.detach().clone().cpu() for name, param in reference_module_fp32.named_parameters()
        }
    else:
        output_global_fp32_host = None
        d_input_x_global_fp32_host = None
        grad_params_fp32_global_host = None

    # Launch parallel test across all processes
    spawn_multiprocessing(
        parallel_assert_triangle_multiplication,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        dtype,
        dim,
        direction,
        layer_state_dict_fp64,
        input_x_global_fp64.detach().clone().cpu(),
        mask_global_fp64.detach().clone().cpu(),
        output_expected_global_fp64.detach().clone().cpu(),
        d_output_expected_global_fp64.detach().clone().cpu(),
        input_x_global_fp64.grad.detach().clone().cpu(),
        grad_params_fp64_expected_global_host,
        output_global_fp32_host,
        d_input_x_global_fp32_host,
        grad_params_fp32_global_host,
    )
