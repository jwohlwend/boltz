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
from boltz.distributed.model.layers.transition import Transition as DistributedTransition
from boltz.model.layers.transition import Transition
from boltz.testing.utils import (
    assert_all_identical,
    assert_no_percentile_upshift,
    assert_tensors_identical,
    get_param_by_key,
    seed_by_rank,
    spawn_multiprocessing,
)


def parallel_assert_transition(
    rank,
    grid_group_sizes,
    device_type,
    backend,
    env_per_rank,
    dtype,
    dim,
    hidden,
    out_dim,
    layer_state_dict,
    input_global_host,
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

    # Create serial reference module
    module_serial = Transition(dim, hidden, out_dim)
    module_serial.load_state_dict(layer_state_dict)
    module_serial = module_serial.to(dtype=dtype, device=manager.device)
    module_serial.train()

    # Create distributed module
    module = DistributedTransition(module_serial, manager.device_mesh_subgroups)
    module.train()

    # Input tensor has shape (B, S, D) - sharded on dims 0 and 1 (B and S)
    placements_input = (Shard(0), Shard(1), Replicate())

    # Distribute input tensor
    input_dtensor = distribute_tensor(
        input_global_host.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements_input,
    ).requires_grad_(True)

    # Distribute expected outputs
    output_expected_dtensor = distribute_tensor(
        output_expected_global_host.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements_input,
        src_data_rank=None,
    )
    d_output_expected_dtensor = distribute_tensor(
        d_output_expected_global_host.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements_input,
    )
    d_input_expected_dtensor = distribute_tensor(
        d_input_expected_global_host.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements_input,
        src_data_rank=None,
    )

    # Create copies to verify inputs aren't modified
    input_dtensor_copy = input_dtensor.detach().clone().requires_grad_(True)

    if check_error_hist:
        # Forward and backward pass for error histogram checking
        output_dtensor_result = module(input_dtensor)
        output_dtensor_result.backward(d_output_expected_dtensor)

        output_fp32_dtensor = distribute_tensor(
            output_global_fp32_host.to(device=manager.device),
            device_mesh=manager.device_mesh_subgroups,
            placements=placements_input,
            src_data_rank=None,
        )

        d_input_fp32_dtensor = distribute_tensor(
            d_input_global_fp32_host.to(device=manager.device),
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
        output_dtensor_result = module(input_dtensor)

        # Verify inputs weren't modified
        assert_tensors_identical(
            input_dtensor_copy.to_local(), input_dtensor.to_local(), check_grad=False, check_grad_fn=False
        )

        # Test forward pass results
        torch.testing.assert_close(output_dtensor_result.to_local(), output_expected_dtensor.to_local())

        # Backward pass
        d_output_expected_dtensor_copy = d_output_expected_dtensor.detach().clone()
        output_dtensor_result.backward(d_output_expected_dtensor)

        # Verify upstream gradient wasn't modified
        assert_tensors_identical(d_output_expected_dtensor_copy.to_local(), d_output_expected_dtensor.to_local())

        # Test input gradients
        torch.testing.assert_close(input_dtensor.grad.to_local(), d_input_expected_dtensor.to_local())

        # Test full tensor gathering - verify distributed results match serial results
        input_global_result_host = input_dtensor.full_tensor().cpu()
        output_global_result_host = output_dtensor_result.full_tensor().cpu()
        d_input_global_result_host = input_dtensor.grad.full_tensor().cpu()

        # Verify full tensors match expected results
        torch.testing.assert_close(input_global_result_host, input_global_host.to(dtype=dtype))
        torch.testing.assert_close(output_global_result_host, output_expected_global_host.to(dtype=dtype))
        torch.testing.assert_close(d_input_global_result_host, d_input_expected_global_host.to(dtype=dtype))

        # Test parameter gradients
        grad_params_result_dtensors = {}
        for name, param in module.named_parameters():
            if param.grad is not None:
                if name not in grad_params_expected_global_host:
                    # do an extra check here to make sure the parallel computation don't result in extra gradients
                    raise ValueError(f"Parameter {name} has a resulting gradient but it is not in the reference module")
                grad_params_result_dtensors[name] = param.grad

        for name, grad_param_expected_global in grad_params_expected_global_host.items():
            assert name in grad_params_result_dtensors, f"Parameter {name}'s gradient is not found in result gradients"
            grad_params_result = grad_params_result_dtensors[name]
            # Test parameter gradients with full tensor gathering
            param_grad_result = grad_params_result.full_tensor()
            torch.testing.assert_close(param_grad_result.cpu(), grad_param_expected_global.to(dtype=dtype))
            assert_all_identical(param_grad_result, manager.group["cp"])

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
    ids=lambda x: f"dp={x[0][0]}, cp={x[0][1]}, specify_method={x[1]}, device_type={x[2]}, method_init={x[3]}",
)
@pytest.mark.parametrize(
    "dtype_and_check_error_hist",
    [
        (torch.float32, False),
        (torch.float32, True),
        (torch.float64, False),
    ],
    ids=lambda x: f"dtype={x[0]}, check_error_hist={x[1]}",
)
def test_transition_parallel(setup_env, dtype_and_check_error_hist):
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    # dtype is the dtype used by the parallel computation
    # check_error_hist determine whether to compare the error histograms between
    # (CP_in_FP32, serial_in_FP64) and (serial_in_FP32, serial_in_FP64)
    # Typically, check_error_hist will use large input dimensions to emulate
    # the real-world use cases. Same with dtype==torch.float64.
    dtype, check_error_hist = dtype_and_check_error_hist

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    if check_error_hist:
        if grid_group_sizes["dp"] > 1:
            pytest.skip("skip error histogram check for dp > 1 to save test time")

    size_ring = grid_group_sizes["cp"][0]
    B = 2 * grid_group_sizes["dp"]
    if check_error_hist or dtype == torch.float64:
        S = size_ring * 128  # Sequence length
        dim = 128  # Input dimension
        hidden = 512  # Hidden dimension
        out_dim = 128  # Output dimension (same as input)
    else:
        S = size_ring * 4  # Sequence length
        dim = 16  # Input dimension
        hidden = 64  # Hidden dimension
        out_dim = 16  # Output dimension (same as input)

    seed = 42
    seed_by_rank(0, seed=seed)

    # compute reference results with FP64
    input_global_fp64 = torch.empty((B, S, dim), dtype=torch.float64, requires_grad=True, device=device_type)

    # Create reference serial module
    reference_module = Transition(dim, hidden, out_dim)

    # Initialize parameters to ensure reproducible behavior
    with torch.no_grad():
        input_global_fp64.uniform_(-5e-2, 5e-2)
        for name, param in reference_module.named_parameters():
            param.uniform_(-5e-2, 5e-2)

    layer_state_dict_fp64 = reference_module.state_dict()
    reference_module = reference_module.to(dtype=torch.float64, device=device_type).train()

    # Run forward pass
    output_expected_global_fp64 = reference_module(input_global_fp64)
    d_output_expected_global_fp64 = torch.rand_like(output_expected_global_fp64)
    output_expected_global_fp64.backward(d_output_expected_global_fp64)

    grad_params_fp64_expected_global_host = {
        name: param.grad.detach().clone().cpu() for name, param in reference_module.named_parameters()
    }

    if check_error_hist:
        input_global_fp32 = input_global_fp64.detach().clone().to(dtype=torch.float32).requires_grad_(True)
        reference_module_fp32 = Transition(dim, hidden, out_dim)
        reference_module_fp32.load_state_dict(layer_state_dict_fp64)
        reference_module_fp32 = reference_module_fp32.to(dtype=torch.float32, device=device_type).train()

        output_global_fp32 = reference_module_fp32(input_global_fp32)
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

    # Launch parallel test across all processes
    spawn_multiprocessing(
        parallel_assert_transition,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        dtype,
        dim,
        hidden,
        out_dim,
        layer_state_dict_fp64,
        input_global_fp64.detach().clone().cpu(),
        output_expected_global_fp64.detach().clone().cpu(),
        d_output_expected_global_fp64.detach().clone().cpu(),
        input_global_fp64.grad.detach().clone().cpu(),
        grad_params_fp64_expected_global_host,
        output_global_fp32_host,
        d_input_global_fp32_host,
        grad_params_fp32_global_host,
    )
