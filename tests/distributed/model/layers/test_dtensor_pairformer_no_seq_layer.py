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

from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.layers.pairformer import (
    PairformerNoSeqLayer as DistributedPairformerNoSeqLayer,
)
from boltz.model.layers.pairformer import (
    PairformerNoSeqLayer as SerialPairformerNoSeqLayer,
)
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


def parallel_assert_pairformer_noseq_layer(
    rank,
    grid_group_sizes,
    device_type,
    backend,
    env_per_rank,
    dtype,
    token_z,
    dropout,
    pairwise_head_width,
    pairwise_num_heads,
    post_layer_norm,
    layer_state_dict,
    input_z_global_host,
    pair_mask_global_host,
    output_z_expected_global_host,
    d_output_z_expected_global_host,
    d_input_z_expected_global_host,
    expected_param_grads_global_host_dict,
    output_z_global_fp32_host: torch.Tensor | None = None,
    d_input_z_global_fp32_host: torch.Tensor | None = None,
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

    if (output_z_global_fp32_host is None) != (d_input_z_global_fp32_host is None) or (
        output_z_global_fp32_host is None
    ) != (grad_params_fp32_global_host is None):
        raise ValueError(
            "output_z_global_fp32_host, d_input_z_global_fp32_host, and grad_params_fp32_global_host "
            "must be either all None or all not None"
        )

    check_error_hist = output_z_global_fp32_host is not None

    # Create serial reference module
    module_serial = SerialPairformerNoSeqLayer(
        token_z=token_z,
        dropout=dropout,
        pairwise_head_width=pairwise_head_width,
        pairwise_num_heads=pairwise_num_heads,
        post_layer_norm=post_layer_norm,
    )
    module_serial.load_state_dict(layer_state_dict)
    set_dtype_specific_inf_values(module_serial, dtype)
    module_serial = module_serial.to(dtype=dtype, device=manager.device)

    # Create distributed module
    module = DistributedPairformerNoSeqLayer(module_serial, manager)
    module.train()

    placements_z_pair_mask = (Shard(0), Shard(1), Shard(2))

    input_z_dtensor = distribute_tensor(
        input_z_global_host.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements_z_pair_mask,
    ).requires_grad_(True)
    pair_mask_dtensor = distribute_tensor(
        pair_mask_global_host.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements_z_pair_mask,
    )

    d_output_z_expected_dtensor = distribute_tensor(
        d_output_z_expected_global_host.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements_z_pair_mask,
    )
    output_z_expected_dtensor = distribute_tensor(
        output_z_expected_global_host.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements_z_pair_mask,
        src_data_rank=None,
    )
    d_input_z_expected_dtensor = distribute_tensor(
        d_input_z_expected_global_host.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements_z_pair_mask,
        src_data_rank=None,
    )

    input_z_dtensor_copy = input_z_dtensor.detach().clone().requires_grad_(True)
    pair_mask_dtensor_copy = pair_mask_dtensor.detach().clone()

    if check_error_hist:
        output_z_dtensor_result = module(z=input_z_dtensor, pair_mask=pair_mask_dtensor)
        output_z_dtensor_result.backward(d_output_z_expected_dtensor)
        output_z_fp32_dtensor = distribute_tensor(
            output_z_global_fp32_host.to(device=manager.device),
            device_mesh=manager.device_mesh_subgroups,
            placements=placements_z_pair_mask,
            src_data_rank=None,
        )
        d_input_z_fp32_dtensor = distribute_tensor(
            d_input_z_global_fp32_host.to(device=manager.device),
            device_mesh=manager.device_mesh_subgroups,
            placements=placements_z_pair_mask,
            src_data_rank=None,
        )
        assert_no_percentile_upshift(
            output_z_dtensor_result.to_local(),
            output_z_expected_dtensor.to_local(),
            output_z_fp32_dtensor.to_local(),
            names_input=("output_z_cp_fp32", "output_z_serial_fp64", "output_z_serial_fp32"),
        )
        assert_no_percentile_upshift(
            input_z_dtensor.grad.to_local(),
            d_input_z_expected_dtensor.to_local(),
            d_input_z_fp32_dtensor.to_local(),
            names_input=("d_input_z_cp_fp32", "d_input_z_serial_fp64", "d_input_z_serial_fp32"),
        )
        for name, grad_param_expected_global in expected_param_grads_global_host_dict.items():
            grad_param_result_global = get_param_by_key(module, name).grad.full_tensor().cpu()
            assert_no_percentile_upshift(
                grad_param_result_global,
                grad_param_expected_global.to(dtype=grad_param_result_global.dtype),
                grad_params_fp32_global_host[name],
                names_input=(f"d_{name}_cp_fp32", f"d_{name}_serial_fp64", f"d_{name}_serial_fp32"),
            )
    else:
        output_z_dtensor_result = module(z=input_z_dtensor, pair_mask=pair_mask_dtensor)
        assert_tensors_identical(
            input_z_dtensor_copy.to_local(), input_z_dtensor.to_local(), check_grad=False, check_grad_fn=False
        )
        assert_tensors_identical(pair_mask_dtensor_copy.to_local(), pair_mask_dtensor.to_local())
        torch.testing.assert_close(output_z_dtensor_result.to_local(), output_z_expected_dtensor.to_local())

        # Clone upstream gradients so we can verify backward does not modify them (match layer test)
        d_output_z_expected_dtensor_copy = d_output_z_expected_dtensor.detach().clone()
        torch.autograd.backward([output_z_dtensor_result], [d_output_z_expected_dtensor])
        # Verify upstream gradients were not modified
        assert_tensors_identical(d_output_z_expected_dtensor_copy.to_local(), d_output_z_expected_dtensor.to_local())

        torch.testing.assert_close(input_z_dtensor.grad.to_local(), d_input_z_expected_dtensor.to_local())

        output_z_global_result_host = output_z_dtensor_result.full_tensor().cpu()
        d_input_z_global_result_host = input_z_dtensor.grad.full_tensor().cpu()
        torch.testing.assert_close(output_z_global_result_host, output_z_expected_global_host.to(dtype=dtype))
        torch.testing.assert_close(d_input_z_global_result_host, d_input_z_expected_global_host.to(dtype=dtype))

        result_param_grads_dict = {}
        for name, param in module.named_parameters():
            if param.grad is not None:
                if name not in expected_param_grads_global_host_dict:
                    raise ValueError(f"Parameter {name} has a resulting gradient but it is not in the reference")
                result_param_grads_dict[name] = param.grad
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
            (((1, (2, 2)), True, "cuda", "ENV"), torch.float32, True),
            (((1, (2, 2)), True, "cuda", "ENV"), torch.float64, False),
            (((2, (2, 2)), True, "cuda", "ENV"), torch.float32, True),
            (((2, (2, 2)), True, "cuda", "ENV"), torch.float64, False),
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
def test_pairformer_noseq_layer_parallel(setup_env, dtype, check_error_hist):
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    if check_error_hist and grid_group_sizes["dp"] > 1:
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
        N = size_ring * 128  # Number of tokens (no_seq: pair dimension only)
        token_z = 128  # Token pairwise embedding dimension
        pairwise_head_width = 32
        pairwise_num_heads = 4
        min_val_init = -0.08 if dtype == torch.float64 else -5e-4
        max_val_init = -min_val_init
    else:
        N = size_ring * 2
        token_z = 12  # Token pairwise embedding dimension
        pairwise_head_width = 4
        pairwise_num_heads = 2
        min_val_init = -0.5
        max_val_init = 0.5
    dropout = 0.0  # disable dropout as we have no way to match the random sequences between serial and CP
    post_layer_norm = False

    seed = 42
    seed_by_rank(0, seed=seed)

    # Compute reference results with FP64
    input_z_global_fp64 = torch.empty((B, N, N, token_z), dtype=torch.float64, requires_grad=True, device=device_type)
    pair_mask_global_fp64 = torch.randint(0, 2, (B, N, N), dtype=torch.float64, requires_grad=False, device=device_type)
    pair_mask_global_fp64[0, N // size_ring :, :] = 0
    pair_mask_global_fp64[0, :, N // size_ring :] = 0

    # Create reference serial module
    reference_module = SerialPairformerNoSeqLayer(
        token_z=token_z,
        dropout=dropout,
        pairwise_head_width=pairwise_head_width,
        pairwise_num_heads=pairwise_num_heads,
        post_layer_norm=post_layer_norm,
    )
    init_tensors_uniform([input_z_global_fp64], low=min_val_init, high=max_val_init)
    init_module_params_uniform(reference_module, low=min_val_init, high=max_val_init)
    set_dtype_specific_inf_values(reference_module, torch.float64)
    layer_state_dict_fp64 = reference_module.state_dict()
    reference_module = reference_module.to(dtype=torch.float64, device=device_type).train()

    output_z_expected_global_fp64 = reference_module(input_z_global_fp64, pair_mask_global_fp64)
    d_output_z_expected_global_fp64 = torch.rand_like(output_z_expected_global_fp64)
    output_z_expected_global_fp64.backward(d_output_z_expected_global_fp64)

    grad_params_fp64_expected_global_host = {
        name: param.grad.detach().to(dtype=dtype, device="cpu", copy=True)
        for name, param in reference_module.named_parameters()
    }

    if check_error_hist:
        input_z_global_fp32 = input_z_global_fp64.detach().to(dtype=torch.float32, copy=True).requires_grad_(True)
        pair_mask_global_fp32 = pair_mask_global_fp64.detach().to(dtype=torch.float32, copy=True).requires_grad_(False)
        reference_module_fp32 = SerialPairformerNoSeqLayer(
            token_z=token_z,
            dropout=dropout,
            pairwise_head_width=pairwise_head_width,
            pairwise_num_heads=pairwise_num_heads,
            post_layer_norm=post_layer_norm,
        )
        reference_module_fp32.load_state_dict(layer_state_dict_fp64)
        reference_module_fp32 = reference_module_fp32.to(dtype=torch.float32, device=device_type).train()
        set_dtype_specific_inf_values(reference_module_fp32, torch.float32)
        output_z_global_fp32 = reference_module_fp32(input_z_global_fp32, pair_mask_global_fp32)
        d_output_z_fp32 = d_output_z_expected_global_fp64.to(dtype=torch.float32)
        output_z_global_fp32.backward(d_output_z_fp32)
        output_z_global_fp32_host = output_z_global_fp32.detach().to(device="cpu", copy=True)
        d_input_z_global_fp32_host = input_z_global_fp32.grad.detach().to(device="cpu", copy=True)
        grad_params_fp32_global_host = {
            name: param.grad.detach().to(device="cpu", copy=True)
            for name, param in reference_module_fp32.named_parameters()
        }
    else:
        output_z_global_fp32_host = None
        d_input_z_global_fp32_host = None
        grad_params_fp32_global_host = None

    # Launch parallel test across all processes
    spawn_multiprocessing(
        parallel_assert_pairformer_noseq_layer,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        dtype,
        token_z,
        dropout,
        pairwise_head_width,
        pairwise_num_heads,
        post_layer_norm,
        layer_state_dict_fp64,
        input_z_global_fp64.detach().to(dtype=dtype, device="cpu", copy=True),
        pair_mask_global_fp64.detach().to(dtype=dtype, device="cpu", copy=True),
        output_z_expected_global_fp64.detach().to(dtype=dtype, device="cpu", copy=True),
        d_output_z_expected_global_fp64.detach().to(dtype=dtype, device="cpu", copy=True),
        input_z_global_fp64.grad.detach().to(dtype=dtype, device="cpu", copy=True),
        grad_params_fp64_expected_global_host,
        output_z_global_fp32_host,
        d_input_z_global_fp32_host,
        grad_params_fp32_global_host,
    )
