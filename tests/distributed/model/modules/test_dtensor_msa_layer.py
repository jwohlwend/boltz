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

"""Tests for Boltz-2 CP MSALayer (distributed.model.modules.trunkv2.MSALayer)."""

import pytest
import torch
from torch.distributed.tensor import Shard, distribute_tensor

from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.modules.trunkv2 import MSALayer as DistributedMSALayer
from boltz.model.modules.trunkv2 import MSALayer as SerialMSALayer
from boltz.testing.utils import (
    assert_all_identical,
    assert_no_percentile_upshift,
    assert_tensors_identical,
    get_param_by_key,
    init_module_params_uniform,
    init_tensors_uniform,
    seed_by_rank,
    set_dtype_specific_inf_values,
    spawn_multiprocessing,
)


def parallel_assert_msa_layer(
    rank,
    grid_group_sizes,
    device_type,
    backend,
    env_per_rank,
    dtype,
    msa_s,
    token_z,
    msa_dropout,
    z_dropout,
    pairwise_head_width,
    pairwise_num_heads,
    layer_state_dict,
    input_z_global_host,
    input_m_global_host,
    token_mask_global_host,
    msa_mask_global_host,
    output_z_expected_global_host,
    output_m_expected_global_host,
    d_output_z_expected_global_host,
    d_output_m_expected_global_host,
    d_input_z_expected_global_host,
    d_input_m_expected_global_host,
    expected_param_grads_global_host_dict,
    output_z_global_fp32_host: torch.Tensor | None = None,
    output_m_global_fp32_host: torch.Tensor | None = None,
    d_input_z_global_fp32_host: torch.Tensor | None = None,
    d_input_m_global_fp32_host: torch.Tensor | None = None,
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

    if torch.finfo(dtype).resolution < torch.finfo(output_z_expected_global_host.dtype).resolution:
        raise ValueError(
            f"Target dtype {dtype} has higher precision than reference output's dtype {output_z_expected_global_host.dtype}"
        )

    if (
        (output_z_global_fp32_host is None) != (output_m_global_fp32_host is None)
        or (output_z_global_fp32_host is None) != (d_input_z_global_fp32_host is None)
        or (output_z_global_fp32_host is None) != (d_input_m_global_fp32_host is None)
        or (output_z_global_fp32_host is None) != (grad_params_fp32_global_host is None)
    ):
        raise ValueError(
            "output_z_global_fp32_host, output_m_global_fp32_host, d_input_z_global_fp32_host, d_input_m_global_fp32_host, and grad_params_fp32_global_host must be either all None or all not None"
        )

    check_error_hist = output_z_global_fp32_host is not None

    # Create serial reference module
    module_serial = SerialMSALayer(
        msa_s=msa_s,
        token_z=token_z,
        msa_dropout=msa_dropout,
        z_dropout=z_dropout,
        pairwise_head_width=pairwise_head_width,
        pairwise_num_heads=pairwise_num_heads,
    )
    module_serial = module_serial.to(dtype=dtype, device=manager.device)
    module_serial.load_state_dict(layer_state_dict)
    set_dtype_specific_inf_values(module_serial, dtype)

    # Create distributed module
    module = DistributedMSALayer(module_serial, manager)
    module.train()

    # Input tensors have sharding patterns:
    # z: (B, N, N, token_z) - sharded on dims 0, 1, 2
    # m: (B, S, N, msa_s) - sharded on dims 0, 1, 2
    # token_mask: (B, N, N) - sharded on dims 0, 1, 2
    # msa_mask: (B, S, N) - sharded on dims 0, 1, 2
    placements_z_token_mask = (Shard(0), Shard(1), Shard(2))  # For z and token_mask tensors
    placements_m_msa_mask = (Shard(0), Shard(1), Shard(2))  # For m and msa_mask tensors

    # Distribute input tensors
    input_z_dtensor = distribute_tensor(
        input_z_global_host.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements_z_token_mask,
    ).requires_grad_(True)

    input_m_dtensor = distribute_tensor(
        input_m_global_host.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements_m_msa_mask,
    ).requires_grad_(True)

    token_mask_dtensor = distribute_tensor(
        token_mask_global_host.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements_z_token_mask,
    )

    msa_mask_dtensor = distribute_tensor(
        msa_mask_global_host.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements_m_msa_mask,
    )

    # Distribute expected outputs
    d_output_z_expected_dtensor = distribute_tensor(
        d_output_z_expected_global_host.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements_z_token_mask,
    )
    d_output_m_expected_dtensor = distribute_tensor(
        d_output_m_expected_global_host.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements_m_msa_mask,
    )
    output_z_expected_dtensor = distribute_tensor(
        output_z_expected_global_host.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements_z_token_mask,
        src_data_rank=None,
    )
    output_m_expected_dtensor = distribute_tensor(
        output_m_expected_global_host.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements_m_msa_mask,
        src_data_rank=None,
    )
    d_input_z_expected_dtensor = distribute_tensor(
        d_input_z_expected_global_host.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements_z_token_mask,
        src_data_rank=None,
    )
    d_input_m_expected_dtensor = distribute_tensor(
        d_input_m_expected_global_host.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements_m_msa_mask,
        src_data_rank=None,
    )

    # Create copies to verify inputs aren't modified
    input_z_dtensor_copy = input_z_dtensor.detach().clone().requires_grad_(True)
    input_m_dtensor_copy = input_m_dtensor.detach().clone().requires_grad_(True)
    token_mask_dtensor_copy = token_mask_dtensor.detach().clone()
    msa_mask_dtensor_copy = msa_mask_dtensor.detach().clone()

    if check_error_hist:
        # Forward and backward pass for error histogram checking
        output_z_dtensor_result, output_m_dtensor_result = module(
            input_z_dtensor, input_m_dtensor, token_mask_dtensor, msa_mask_dtensor
        )
        torch.autograd.backward(
            [output_z_dtensor_result, output_m_dtensor_result],
            [d_output_z_expected_dtensor, d_output_m_expected_dtensor],
        )

        # Distribute FP32 reference results for comparison
        output_z_fp32_dtensor = distribute_tensor(
            output_z_global_fp32_host.to(device=manager.device),
            device_mesh=manager.device_mesh_subgroups,
            placements=placements_z_token_mask,
            src_data_rank=None,
        )
        output_m_fp32_dtensor = distribute_tensor(
            output_m_global_fp32_host.to(device=manager.device),
            device_mesh=manager.device_mesh_subgroups,
            placements=placements_m_msa_mask,
            src_data_rank=None,
        )
        d_input_z_fp32_dtensor = distribute_tensor(
            d_input_z_global_fp32_host.to(device=manager.device),
            device_mesh=manager.device_mesh_subgroups,
            placements=placements_z_token_mask,
            src_data_rank=None,
        )
        d_input_m_fp32_dtensor = distribute_tensor(
            d_input_m_global_fp32_host.to(device=manager.device),
            device_mesh=manager.device_mesh_subgroups,
            placements=placements_m_msa_mask,
            src_data_rank=None,
        )

        # Check error histograms for outputs
        assert_no_percentile_upshift(
            output_z_dtensor_result.to_local(),
            output_z_expected_dtensor.to_local(),
            output_z_fp32_dtensor.to_local(),
            names_input=("output_z_cp_fp32", "output_z_serial_fp64", "output_z_serial_fp32"),
        )

        assert_no_percentile_upshift(
            output_m_dtensor_result.to_local(),
            output_m_expected_dtensor.to_local(),
            output_m_fp32_dtensor.to_local(),
            names_input=("output_m_cp_fp32", "output_m_serial_fp64", "output_m_serial_fp32"),
        )

        # Check error histograms for input gradients
        assert_no_percentile_upshift(
            input_z_dtensor.grad.to_local(),
            d_input_z_expected_dtensor.to_local(),
            d_input_z_fp32_dtensor.to_local(),
            names_input=("d_input_z_cp_fp32", "d_input_z_serial_fp64", "d_input_z_serial_fp32"),
        )

        assert_no_percentile_upshift(
            input_m_dtensor.grad.to_local(),
            d_input_m_expected_dtensor.to_local(),
            d_input_m_fp32_dtensor.to_local(),
            names_input=("d_input_m_cp_fp32", "d_input_m_serial_fp64", "d_input_m_serial_fp32"),
        )

        # Check parameter gradients error histograms
        for name, grad_param_expected_global in expected_param_grads_global_host_dict.items():
            grad_param_result_global = get_param_by_key(module, name).grad.full_tensor().cpu()
            assert_no_percentile_upshift(
                grad_param_result_global,
                grad_param_expected_global.to(dtype=grad_param_result_global.dtype),
                grad_params_fp32_global_host[name],
                names_input=(f"d_{name}_cp_fp32", f"d_{name}_serial_fp64", f"d_{name}_serial_fp32"),
            )
    else:
        # Forward pass
        output_z_dtensor_result, output_m_dtensor_result = module(
            input_z_dtensor, input_m_dtensor, token_mask_dtensor, msa_mask_dtensor
        )

        # Verify inputs weren't modified
        assert_tensors_identical(
            input_z_dtensor_copy.to_local(), input_z_dtensor.to_local(), check_grad=False, check_grad_fn=False
        )
        assert_tensors_identical(
            input_m_dtensor_copy.to_local(), input_m_dtensor.to_local(), check_grad=False, check_grad_fn=False
        )
        assert_tensors_identical(token_mask_dtensor_copy.to_local(), token_mask_dtensor.to_local())
        assert_tensors_identical(msa_mask_dtensor_copy.to_local(), msa_mask_dtensor.to_local())

        # Test forward pass results
        torch.testing.assert_close(output_z_dtensor_result.to_local(), output_z_expected_dtensor.to_local())
        torch.testing.assert_close(output_m_dtensor_result.to_local(), output_m_expected_dtensor.to_local())

        # Backward pass
        d_output_z_expected_dtensor_copy = d_output_z_expected_dtensor.detach().clone()
        d_output_m_expected_dtensor_copy = d_output_m_expected_dtensor.detach().clone()
        torch.autograd.backward(
            [output_z_dtensor_result, output_m_dtensor_result],
            [d_output_z_expected_dtensor, d_output_m_expected_dtensor],
        )

        # Verify upstream gradients weren't modified
        assert_tensors_identical(d_output_z_expected_dtensor_copy.to_local(), d_output_z_expected_dtensor.to_local())
        assert_tensors_identical(d_output_m_expected_dtensor_copy.to_local(), d_output_m_expected_dtensor.to_local())

        # Test input gradients
        torch.testing.assert_close(input_z_dtensor.grad.to_local(), d_input_z_expected_dtensor.to_local())
        torch.testing.assert_close(input_m_dtensor.grad.to_local(), d_input_m_expected_dtensor.to_local())

        # Test full tensor gathering - verify distributed results match serial results
        output_z_global_result_host = output_z_dtensor_result.full_tensor().cpu()
        output_m_global_result_host = output_m_dtensor_result.full_tensor().cpu()
        d_input_z_global_result_host = input_z_dtensor.grad.full_tensor().cpu()
        d_input_m_global_result_host = input_m_dtensor.grad.full_tensor().cpu()

        # Verify full tensors match expected results
        torch.testing.assert_close(output_z_global_result_host, output_z_expected_global_host.to(dtype=dtype))
        torch.testing.assert_close(output_m_global_result_host, output_m_expected_global_host.to(dtype=dtype))
        torch.testing.assert_close(d_input_z_global_result_host, d_input_z_expected_global_host.to(dtype=dtype))
        torch.testing.assert_close(d_input_m_global_result_host, d_input_m_expected_global_host.to(dtype=dtype))

        # Gather weight gradients using named_parameters
        # NOTE: the layer weights are all replicated and their gradients are in Partial(Sum) state
        # of their dtensor form so testing the full_tensor() results is equivalent to testing the
        # DTensor versions
        result_param_grads_dict = {}
        for name, param in module.named_parameters():
            if param.grad is not None:
                if name not in expected_param_grads_global_host_dict:
                    raise ValueError(f"Parameter {name} has a resulting gradient but it is not in the reference module")
                result_param_grads_dict[name] = param.grad

        # Compare parameter gradients
        for name, expected_grad_global_host in expected_param_grads_global_host_dict.items():
            assert name in result_param_grads_dict, f"Parameter {name}'s gradient is not found in result gradients"
            result_grad = result_param_grads_dict[name]
            result_grad_global = result_grad.full_tensor()
            torch.testing.assert_close(result_grad_global.cpu(), expected_grad_global_host.to(dtype=dtype))
            assert_all_identical(result_grad_global, manager.group["cp"])

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.slow
@pytest.mark.parametrize(
    "setup_env, dtype, check_error_hist",
    (
        params_test := [
            ## CUDA tests (2 GPUs)
            (((2, (1, 1)), True, "cuda", "ENV"), torch.float32, True),
            (((2, (1, 1)), True, "cuda", "ENV"), torch.float64, True),
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
def test_msa_layer_parallel(setup_env, dtype, check_error_hist):
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
        N = size_ring * 128  # Number of tokens
        S = size_ring * 128  # Number of sequences
        msa_s = 64  # MSA embedding dimension
        token_z = 128  # Pairwise embedding dimension
        pairwise_head_width = 32
        pairwise_num_heads = 4
        min_val_init = -5e-2 if dtype == torch.float64 else -1e-3
        max_val_init = -min_val_init
    else:
        N = size_ring * 2  # Number of tokens
        S = size_ring * 3  # Number of sequences
        msa_s = 8  # MSA embedding dimension
        token_z = 12  # Pairwise embedding dimension
        pairwise_head_width = 4
        pairwise_num_heads = 2
        min_val_init = -0.5
        max_val_init = 0.5
    msa_dropout = 0.0  # disable dropout as we have not way to match the random sequences between serial and CP
    z_dropout = 0.0

    seed = 42
    seed_by_rank(0, seed=seed)

    # Compute reference results with FP64
    input_z_global_fp64 = torch.empty((B, N, N, token_z), dtype=torch.float64, requires_grad=True, device=device_type)
    input_m_global_fp64 = torch.empty((B, S, N, msa_s), dtype=torch.float64, requires_grad=True, device=device_type)

    token_mask_global_fp64 = torch.randint(
        0, 2, (B, N, N), dtype=torch.float64, requires_grad=False, device=device_type
    )
    token_mask_global_fp64[0, N // size_ring :, :] = 0
    token_mask_global_fp64[0, :, N // size_ring :] = 0

    msa_mask_global_fp64 = torch.ones((B, S, N), dtype=torch.float64, requires_grad=False, device=device_type)
    msa_mask_global_fp64[0, (S // size_ring) :, :] = 0
    msa_mask_global_fp64[0, :, (N // size_ring) :] = 0

    # Create reference serial module
    reference_module = SerialMSALayer(
        msa_s=msa_s,
        token_z=token_z,
        msa_dropout=msa_dropout,
        z_dropout=z_dropout,
        pairwise_head_width=pairwise_head_width,
        pairwise_num_heads=pairwise_num_heads,
    )

    # Initialize parameters to ensure reproducible behavior
    # The output activation and gradient of the layer weights typically increase by 3 to 4 orders of magnitude,
    # where the ULP would be too large and numerical error distribution becomes very wide, i.e., we would have
    # very unpredictable numerical errors. That would make the test results very noisy and not very useful to
    # detect logical bugs in the code. To avoid this, we use a smaller range for the input and layer weights.
    init_tensors_uniform([input_z_global_fp64, input_m_global_fp64], low=min_val_init, high=max_val_init)
    reference_module = reference_module.to(dtype=torch.float64, device=device_type).train()
    init_module_params_uniform(reference_module, low=min_val_init, high=max_val_init)

    set_dtype_specific_inf_values(reference_module, torch.float64)

    layer_state_dict_fp64 = reference_module.state_dict()

    # Run forward pass
    output_z_expected_global_fp64, output_m_expected_global_fp64 = reference_module(
        input_z_global_fp64, input_m_global_fp64, token_mask_global_fp64, msa_mask_global_fp64
    )
    d_output_z_expected_global_fp64 = torch.rand_like(output_z_expected_global_fp64)
    d_output_m_expected_global_fp64 = torch.rand_like(output_m_expected_global_fp64)
    torch.autograd.backward(
        [output_z_expected_global_fp64, output_m_expected_global_fp64],
        [d_output_z_expected_global_fp64, d_output_m_expected_global_fp64],
    )

    grad_params_fp64_expected_global_host = {
        name: param.grad.detach().to(dtype=dtype, device="cpu", copy=True)
        for name, param in reference_module.named_parameters()
    }

    if check_error_hist:
        # Run serial FP32 reference for three-way error histogram comparison
        input_z_global_fp32 = input_z_global_fp64.detach().to(dtype=torch.float32, copy=True).requires_grad_(True)
        input_m_global_fp32 = input_m_global_fp64.detach().to(dtype=torch.float32, copy=True).requires_grad_(True)
        token_mask_global_fp32 = (
            token_mask_global_fp64.detach().to(dtype=torch.float32, copy=True).requires_grad_(False)
        )
        msa_mask_global_fp32 = msa_mask_global_fp64.detach().to(dtype=torch.float32, copy=True).requires_grad_(False)

        reference_module_fp32 = SerialMSALayer(
            msa_s=msa_s,
            token_z=token_z,
            msa_dropout=msa_dropout,
            z_dropout=z_dropout,
            pairwise_head_width=pairwise_head_width,
            pairwise_num_heads=pairwise_num_heads,
        )

        reference_module_fp32.load_state_dict(layer_state_dict_fp64)
        reference_module_fp32 = reference_module_fp32.to(dtype=torch.float32, device=device_type).train()
        set_dtype_specific_inf_values(reference_module_fp32, torch.float32)

        output_z_global_fp32, output_m_global_fp32 = reference_module_fp32(
            input_z_global_fp32, input_m_global_fp32, token_mask_global_fp32, msa_mask_global_fp32
        )
        d_output_z_expected_global_fp32 = d_output_z_expected_global_fp64.to(dtype=torch.float32)
        d_output_m_expected_global_fp32 = d_output_m_expected_global_fp64.to(dtype=torch.float32)
        torch.autograd.backward(
            [output_z_global_fp32, output_m_global_fp32],
            [d_output_z_expected_global_fp32, d_output_m_expected_global_fp32],
        )

        output_z_global_fp32_host = output_z_global_fp32.detach().to(device="cpu", copy=True)
        output_m_global_fp32_host = output_m_global_fp32.detach().to(device="cpu", copy=True)
        d_input_z_global_fp32_host = input_z_global_fp32.grad.detach().to(device="cpu", copy=True)
        d_input_m_global_fp32_host = input_m_global_fp32.grad.detach().to(device="cpu", copy=True)
        grad_params_fp32_global_host = {
            name: param.grad.detach().to(device="cpu", copy=True)
            for name, param in reference_module_fp32.named_parameters()
        }

        output_z_for_worker = output_z_expected_global_fp64.detach().to(dtype=dtype, device="cpu", copy=True)
        output_m_for_worker = output_m_expected_global_fp64.detach().to(dtype=dtype, device="cpu", copy=True)
        d_input_z_for_worker = input_z_global_fp64.grad.detach().to(dtype=dtype, device="cpu", copy=True)
        d_input_m_for_worker = input_m_global_fp64.grad.detach().to(dtype=dtype, device="cpu", copy=True)
        grad_params_for_worker = grad_params_fp64_expected_global_host
    elif dtype == torch.float32:
        # check_error_hist=False with FP32: the spawned worker compares the CP output
        # directly against the expected output via assert_close.  Because the FP64
        # reference uses genuine FP64 parameters while the CP module's parameters are
        # truncated to FP32 on load, numerical discrepancies from both the parameter
        # truncation and the lower-precision arithmetic accumulate through the composed
        # MSA layer and exceed assert_close tolerances.  To avoid this, we run a serial
        # FP32 reference so that both sides start from identical parameters.
        ref_fp32 = SerialMSALayer(
            msa_s=msa_s,
            token_z=token_z,
            msa_dropout=msa_dropout,
            z_dropout=z_dropout,
            pairwise_head_width=pairwise_head_width,
            pairwise_num_heads=pairwise_num_heads,
        )
        ref_fp32.load_state_dict(layer_state_dict_fp64)
        ref_fp32 = ref_fp32.to(dtype=torch.float32, device=device_type).train()
        set_dtype_specific_inf_values(ref_fp32, torch.float32)

        inp_z = input_z_global_fp64.detach().to(dtype=torch.float32, device=device_type).requires_grad_(True)
        inp_m = input_m_global_fp64.detach().to(dtype=torch.float32, device=device_type).requires_grad_(True)
        tok_mask = token_mask_global_fp64.detach().to(dtype=torch.float32, device=device_type)
        msa_msk = msa_mask_global_fp64.detach().to(dtype=torch.float32, device=device_type)

        out_z, out_m = ref_fp32(inp_z, inp_m, tok_mask, msa_msk)
        d_out_z = d_output_z_expected_global_fp64.to(dtype=torch.float32)
        d_out_m = d_output_m_expected_global_fp64.to(dtype=torch.float32)
        torch.autograd.backward([out_z, out_m], [d_out_z, d_out_m])

        output_z_for_worker = out_z.detach().cpu()
        output_m_for_worker = out_m.detach().cpu()
        d_input_z_for_worker = inp_z.grad.detach().cpu()
        d_input_m_for_worker = inp_m.grad.detach().cpu()
        grad_params_for_worker = {name: param.grad.detach().cpu() for name, param in ref_fp32.named_parameters()}

        output_z_global_fp32_host = None
        output_m_global_fp32_host = None
        d_input_z_global_fp32_host = None
        d_input_m_global_fp32_host = None
        grad_params_fp32_global_host = None
    else:
        # check_error_hist=False with FP64: use FP64 reference directly
        output_z_for_worker = output_z_expected_global_fp64.detach().to(dtype=dtype, device="cpu", copy=True)
        output_m_for_worker = output_m_expected_global_fp64.detach().to(dtype=dtype, device="cpu", copy=True)
        d_input_z_for_worker = input_z_global_fp64.grad.detach().to(dtype=dtype, device="cpu", copy=True)
        d_input_m_for_worker = input_m_global_fp64.grad.detach().to(dtype=dtype, device="cpu", copy=True)
        grad_params_for_worker = grad_params_fp64_expected_global_host

        output_z_global_fp32_host = None
        output_m_global_fp32_host = None
        d_input_z_global_fp32_host = None
        d_input_m_global_fp32_host = None
        grad_params_fp32_global_host = None

    # Launch parallel test across all processes
    spawn_multiprocessing(
        parallel_assert_msa_layer,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        dtype,
        msa_s,
        token_z,
        msa_dropout,
        z_dropout,
        pairwise_head_width,
        pairwise_num_heads,
        layer_state_dict_fp64,
        input_z_global_fp64.detach().to(dtype=dtype, device="cpu", copy=True),
        input_m_global_fp64.detach().to(dtype=dtype, device="cpu", copy=True),
        token_mask_global_fp64.detach().to(dtype=dtype, device="cpu", copy=True),
        msa_mask_global_fp64.detach().to(dtype=dtype, device="cpu", copy=True),
        output_z_for_worker,
        output_m_for_worker,
        d_output_z_expected_global_fp64.detach().to(dtype=dtype, device="cpu", copy=True),
        d_output_m_expected_global_fp64.detach().to(dtype=dtype, device="cpu", copy=True),
        d_input_z_for_worker,
        d_input_m_for_worker,
        grad_params_for_worker,
        output_z_global_fp32_host,
        output_m_global_fp32_host,
        d_input_z_global_fp32_host,
        d_input_m_global_fp32_host,
        grad_params_fp32_global_host,
    )
