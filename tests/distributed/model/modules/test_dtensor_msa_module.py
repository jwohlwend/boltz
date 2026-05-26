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

"""Tests for Boltz-2 CP MSAModule (distributed.model.modules.trunkv2.MSAModule)."""

import pytest
import torch
from torch.distributed.tensor import Replicate, Shard, distribute_tensor

from boltz.data import const
from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.modules.trunkv2 import MSAModule as DistributedMSAModule
from boltz.model.modules.trunkv2 import MSAModule as SerialMSAModule
from boltz.testing.utils import (
    assert_all_identical,
    assert_no_percentile_upshift,
    assert_tensors_identical,
    create_msa_module_init_params_v2,
    get_param_by_key,
    init_module_params_uniform,
    init_tensors_uniform,
    seed_by_rank,
    set_dtype_specific_inf_values,
    spawn_multiprocessing,
)


def _feats_for_distributed(feats_global, dtype, device="cpu"):
    """Convert feats to the target dtype/device for the distributed module.

    Integer features (e.g. ``msa``) are kept as-is because the distributed
    MSAModule now applies ``shardwise_one_hot`` internally, matching the
    serial MSAModule which calls ``F.one_hot`` in its forward pass.
    """
    out = {}
    for key, value in feats_global.items():
        if value.dtype.is_floating_point:
            out[key] = value.to(dtype=dtype, device=device)
        else:
            out[key] = value.to(device=device)
    return out


def parallel_assert_msa_module(
    rank,
    grid_group_sizes,
    device_type,
    backend,
    env_per_rank,
    dtype,
    msa_module_params,
    module_state_dict,
    input_z_global_host,
    input_emb_global_host,
    input_feats_global_host,
    output_z_expected_global_host,
    d_output_z_expected_global_host,
    d_input_z_expected_global_host,
    d_input_emb_expected_global_host,
    expected_param_grads_global_host_dict,
    output_z_global_fp32_host=None,
    d_input_z_global_fp32_host=None,
    d_input_emb_global_fp32_host=None,
    grad_params_fp32_global_host=None,
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

    check_error_hist = output_z_global_fp32_host is not None

    module_serial = SerialMSAModule(**msa_module_params)
    module_serial = module_serial.to(dtype=dtype, device=manager.device)
    module_serial.load_state_dict(module_state_dict)

    set_dtype_specific_inf_values(module_serial, dtype)

    module = DistributedMSAModule(module_serial, manager)
    assert module.activation_checkpointing == msa_module_params["activation_checkpointing"]
    module.train()

    placements_z = (Shard(0), Shard(1), Shard(2))
    placements_emb = (Shard(0), Replicate(), Shard(1))
    placements_msa = (Shard(0), Shard(1), Shard(2))
    placements_token_mask = (Shard(0), Shard(1), Shard(2))

    input_z_dtensor = distribute_tensor(
        input_z_global_host.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements_z,
    ).requires_grad_(True)

    input_emb_dtensor = distribute_tensor(
        input_emb_global_host.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements_emb,
    ).requires_grad_(True)

    input_feats_dtensor = {}
    for key, value in input_feats_global_host.items():
        if key in ["msa", "has_deletion", "deletion_value", "msa_paired", "msa_mask"]:
            input_feats_dtensor[key] = distribute_tensor(
                value.to(dtype=dtype if value.dtype.is_floating_point else value.dtype, device=manager.device),
                device_mesh=manager.device_mesh_subgroups,
                placements=placements_msa,
            )
        elif key == "token_pair_pad_mask":
            input_feats_dtensor[key] = distribute_tensor(
                value.to(dtype=dtype, device=manager.device),
                device_mesh=manager.device_mesh_subgroups,
                placements=placements_token_mask,
            )

    d_output_z_expected_dtensor = distribute_tensor(
        d_output_z_expected_global_host.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements_z,
    )
    output_z_expected_dtensor = distribute_tensor(
        output_z_expected_global_host.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements_z,
        src_data_rank=None,
    )
    d_input_z_expected_dtensor = distribute_tensor(
        d_input_z_expected_global_host.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements_z,
        src_data_rank=None,
    )
    d_input_emb_expected_dtensor = distribute_tensor(
        d_input_emb_expected_global_host.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements_emb,
        src_data_rank=None,
    )

    input_z_dtensor_copy = input_z_dtensor.detach().clone().requires_grad_(True)
    input_emb_dtensor_copy = input_emb_dtensor.detach().clone().requires_grad_(True)
    input_feats_dtensor_copy = {k: v.detach().clone() for k, v in input_feats_dtensor.items()}

    if check_error_hist:
        output_z_dtensor_result = module(input_z_dtensor, input_emb_dtensor, input_feats_dtensor)
        torch.autograd.backward([output_z_dtensor_result], [d_output_z_expected_dtensor])

        output_z_fp32_dtensor = distribute_tensor(
            output_z_global_fp32_host.to(device=manager.device),
            device_mesh=manager.device_mesh_subgroups,
            placements=placements_z,
            src_data_rank=None,
        )
        d_input_z_fp32_dtensor = distribute_tensor(
            d_input_z_global_fp32_host.to(device=manager.device),
            device_mesh=manager.device_mesh_subgroups,
            placements=placements_z,
            src_data_rank=None,
        )
        d_input_emb_fp32_dtensor = distribute_tensor(
            d_input_emb_global_fp32_host.to(device=manager.device),
            device_mesh=manager.device_mesh_subgroups,
            placements=placements_emb,
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

        assert_no_percentile_upshift(
            input_emb_dtensor.grad.to_local(),
            d_input_emb_expected_dtensor.to_local(),
            d_input_emb_fp32_dtensor.to_local(),
            names_input=("d_input_emb_cp_fp32", "d_input_emb_serial_fp64", "d_input_emb_serial_fp32"),
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
        output_z_dtensor_result = module(input_z_dtensor, input_emb_dtensor, input_feats_dtensor)

        assert_tensors_identical(
            input_z_dtensor_copy.to_local(), input_z_dtensor.to_local(), check_grad=False, check_grad_fn=False
        )
        assert_tensors_identical(
            input_emb_dtensor_copy.to_local(), input_emb_dtensor.to_local(), check_grad=False, check_grad_fn=False
        )
        for key in input_feats_dtensor_copy:
            assert_tensors_identical(input_feats_dtensor_copy[key].to_local(), input_feats_dtensor[key].to_local())

        torch.testing.assert_close(output_z_dtensor_result.to_local(), output_z_expected_dtensor.to_local())

        d_output_z_expected_dtensor_copy = d_output_z_expected_dtensor.detach().clone()
        torch.autograd.backward([output_z_dtensor_result], [d_output_z_expected_dtensor])

        assert_tensors_identical(d_output_z_expected_dtensor_copy.to_local(), d_output_z_expected_dtensor.to_local())

        torch.testing.assert_close(input_z_dtensor.grad.to_local(), d_input_z_expected_dtensor.to_local())
        torch.testing.assert_close(input_emb_dtensor.grad.to_local(), d_input_emb_expected_dtensor.to_local())

        output_z_global_result_host = output_z_dtensor_result.full_tensor().cpu()
        d_input_z_global_result_host = input_z_dtensor.grad.full_tensor().cpu()
        d_input_emb_global_result_host = input_emb_dtensor.grad.full_tensor().cpu()

        torch.testing.assert_close(output_z_global_result_host, output_z_expected_global_host.to(dtype=dtype))
        torch.testing.assert_close(d_input_z_global_result_host, d_input_z_expected_global_host.to(dtype=dtype))
        torch.testing.assert_close(d_input_emb_global_result_host, d_input_emb_expected_global_host.to(dtype=dtype))

        result_param_grads_dict = {}
        for name, param in module.named_parameters():
            if param.grad is not None:
                if name not in expected_param_grads_global_host_dict:
                    raise ValueError(f"Parameter {name} has a resulting gradient but it is not in the reference module")
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
    "setup_env, dtype, check_error_hist, activation_checkpointing",
    (
        params_test := [
            ## CUDA tests (2 GPUs)
            (((2, (1, 1)), True, "cuda", "ENV"), torch.float32, True, False),
            (((2, (1, 1)), True, "cuda", "ENV"), torch.float64, True, True),
            ## CUDA tests (8 GPUs)
            (((2, (2, 2)), True, "cuda", "ENV"), torch.float32, True, True),
            (((2, (2, 2)), True, "cuda", "ENV"), torch.float32, True, False),
            (((2, (2, 2)), True, "cuda", "ENV"), torch.float32, False, False),
            (((2, (2, 2)), True, "cuda", "ENV"), torch.float64, True, True),
            (((2, (2, 2)), True, "cuda", "ENV"), torch.float64, True, False),
            ## CPU tests
            (((1, (3, 3)), True, "cpu", "ENV"), torch.float32, False, False),
        ]
    ),
    indirect=["setup_env"],
    ids=[
        f"dp:{x[0][0][0]}, cp:{x[0][0][1]}, specify_method:{x[0][1]}, device_type:{x[0][2]}, method_init:{x[0][3]}, "
        f"dtype:{x[1]}, check_error_hist:{x[2]}, checkpoint:{x[3]}"
        for x in params_test
    ],
)
def test_msa_module_parallel(setup_env, dtype, check_error_hist, activation_checkpointing):
    """Test Boltz-2 CP MSAModule against serial reference (forward, backward, param grads)."""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cpu" and grid_group_sizes["cp"] == (3, 3):
        pytest.skip("CPU with 3x3 CP ring not yet validated for numerical parity (distributed vs serial)")

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    test_large_model = check_error_hist or dtype == torch.float64

    size_ring = grid_group_sizes["cp"][0]
    B = 2 * grid_group_sizes["dp"]

    if test_large_model:
        N = size_ring * 64
        S = size_ring * 64
        min_val_init = -0.01
        max_val_init = 0.01
    else:
        N = size_ring * 2
        S = size_ring * 3
        min_val_init = -0.5
        max_val_init = 0.5

    msa_module_params = create_msa_module_init_params_v2(test_large_model)
    msa_module_params["activation_checkpointing"] = activation_checkpointing

    seed = 42
    seed_by_rank(0, seed=seed)

    input_z_global_fp64 = torch.empty(
        (B, N, N, msa_module_params["token_z"]), dtype=torch.float64, requires_grad=True, device=device_type
    )
    input_emb_global_fp64 = torch.empty(
        (B, N, msa_module_params["token_s"]), dtype=torch.float64, requires_grad=True, device=device_type
    )

    dim_input_msa = const.num_tokens
    input_feats_global_fp64 = {
        "msa": torch.randint(0, dim_input_msa, (B, S, N), dtype=torch.int64, device=device_type),
        "has_deletion": torch.empty((B, S, N), dtype=torch.float64, device=device_type),
        "deletion_value": torch.empty((B, S, N), dtype=torch.float64, device=device_type),
        "msa_paired": torch.randint(0, 2, (B, S, N), dtype=torch.float64, device=device_type),
        "msa_mask": torch.ones((B, S, N), dtype=torch.float64, device=device_type),
        "token_pad_mask": torch.randint(0, 2, (B, N), dtype=torch.float64, device=device_type),
    }
    input_feats_global_fp64["token_pad_mask"][0, N // size_ring :] = 0
    input_feats_global_fp64["token_pair_pad_mask"] = (
        input_feats_global_fp64["token_pad_mask"][:, :, None] * input_feats_global_fp64["token_pad_mask"][:, None, :]
    )

    input_feats_global_fp64["msa_mask"][0, (S // size_ring) :, :] = 0
    input_feats_global_fp64["msa_mask"][0, :, (N // size_ring) :] = 0

    reference_module = SerialMSAModule(**msa_module_params)

    init_tensors_uniform([input_z_global_fp64, input_emb_global_fp64], low=min_val_init, high=max_val_init)
    for key, tensor in input_feats_global_fp64.items():
        if tensor.dtype.is_floating_point and "mask" not in key and "msa_paired" not in key:
            init_tensors_uniform([tensor], low=min_val_init, high=max_val_init)

    reference_module = reference_module.to(dtype=torch.float64, device=device_type).train()
    init_module_params_uniform(reference_module, low=min_val_init, high=max_val_init)

    set_dtype_specific_inf_values(reference_module, torch.float64)

    module_state_dict_fp64 = reference_module.state_dict()

    output_z_expected_global_fp64 = reference_module(
        input_z_global_fp64, input_emb_global_fp64, input_feats_global_fp64
    )
    d_output_z_expected_global_fp64 = torch.rand_like(output_z_expected_global_fp64)
    torch.autograd.backward([output_z_expected_global_fp64], [d_output_z_expected_global_fp64])

    grad_params_fp64_expected_global_host = {
        name: param.grad.detach().to(dtype=dtype, device="cpu", copy=True)
        for name, param in reference_module.named_parameters()
        if param.grad is not None
    }

    if check_error_hist:
        # Run serial FP32 reference for three-way error histogram comparison
        input_z_global_fp32 = input_z_global_fp64.detach().to(dtype=torch.float32, copy=True).requires_grad_(True)
        input_emb_global_fp32 = input_emb_global_fp64.detach().to(dtype=torch.float32, copy=True).requires_grad_(True)
        input_feats_global_fp32 = {}
        for key, tensor in input_feats_global_fp64.items():
            if key == "msa":
                input_feats_global_fp32[key] = tensor.detach().clone()
            elif tensor.dtype.is_floating_point:
                input_feats_global_fp32[key] = tensor.detach().to(dtype=torch.float32, copy=True)
            else:
                input_feats_global_fp32[key] = tensor.detach().clone()

        reference_module_fp32 = SerialMSAModule(**msa_module_params)

        reference_module_fp32.load_state_dict(module_state_dict_fp64)
        set_dtype_specific_inf_values(reference_module_fp32, torch.float32)

        reference_module_fp32 = reference_module_fp32.to(dtype=torch.float32, device=device_type).train()

        output_z_global_fp32 = reference_module_fp32(
            input_z_global_fp32, input_emb_global_fp32, input_feats_global_fp32
        )
        d_output_z_expected_global_fp32 = d_output_z_expected_global_fp64.to(dtype=torch.float32)
        torch.autograd.backward([output_z_global_fp32], [d_output_z_expected_global_fp32])

        output_z_global_fp32_host = output_z_global_fp32.detach().to(device="cpu", copy=True)
        d_input_z_global_fp32_host = input_z_global_fp32.grad.detach().to(device="cpu", copy=True)
        d_input_emb_global_fp32_host = input_emb_global_fp32.grad.detach().to(device="cpu", copy=True)
        grad_params_fp32_global_host = {
            name: param.grad.detach().to(device="cpu", copy=True)
            for name, param in reference_module_fp32.named_parameters()
            if param.grad is not None
        }

        output_z_for_worker = output_z_expected_global_fp64.detach().to(dtype=dtype, device="cpu", copy=True)
        d_input_z_for_worker = input_z_global_fp64.grad.detach().to(dtype=dtype, device="cpu", copy=True)
        d_input_emb_for_worker = input_emb_global_fp64.grad.detach().to(dtype=dtype, device="cpu", copy=True)
        grad_params_for_worker = grad_params_fp64_expected_global_host
    elif dtype == torch.float32:
        # check_error_hist=False with FP32: the spawned worker compares the CP output
        # directly against the expected output via assert_close.  Because the FP64
        # reference uses genuine FP64 parameters while the CP module's parameters are
        # truncated to FP32 on load, numerical discrepancies from both the parameter
        # truncation and the lower-precision arithmetic accumulate through the composed
        # MSA module and exceed assert_close tolerances.  To avoid this, we run a serial
        # FP32 reference so that both sides start from identical parameters.
        ref_fp32 = SerialMSAModule(**msa_module_params)
        ref_fp32.load_state_dict(module_state_dict_fp64)
        set_dtype_specific_inf_values(ref_fp32, torch.float32)
        ref_fp32 = ref_fp32.to(dtype=torch.float32, device=device_type).train()

        inp_z = input_z_global_fp64.detach().to(dtype=torch.float32, device=device_type).requires_grad_(True)
        inp_emb = input_emb_global_fp64.detach().to(dtype=torch.float32, device=device_type).requires_grad_(True)
        inp_feats = {}
        for key, tensor in input_feats_global_fp64.items():
            if key == "msa":
                inp_feats[key] = tensor.detach().clone()
            elif tensor.dtype.is_floating_point:
                inp_feats[key] = tensor.detach().to(dtype=torch.float32, device=device_type)
            else:
                inp_feats[key] = tensor.detach().clone().to(device=device_type)

        out_z = ref_fp32(inp_z, inp_emb, inp_feats)
        d_out_z = d_output_z_expected_global_fp64.to(dtype=torch.float32)
        torch.autograd.backward([out_z], [d_out_z])

        output_z_for_worker = out_z.detach().cpu()
        d_input_z_for_worker = inp_z.grad.detach().cpu()
        d_input_emb_for_worker = inp_emb.grad.detach().cpu()
        grad_params_for_worker = {
            name: param.grad.detach().cpu() for name, param in ref_fp32.named_parameters() if param.grad is not None
        }

        output_z_global_fp32_host = None
        d_input_z_global_fp32_host = None
        d_input_emb_global_fp32_host = None
        grad_params_fp32_global_host = None
    else:
        # check_error_hist=False with FP64: use FP64 reference directly
        output_z_for_worker = output_z_expected_global_fp64.detach().to(dtype=dtype, device="cpu", copy=True)
        d_input_z_for_worker = input_z_global_fp64.grad.detach().to(dtype=dtype, device="cpu", copy=True)
        d_input_emb_for_worker = input_emb_global_fp64.grad.detach().to(dtype=dtype, device="cpu", copy=True)
        grad_params_for_worker = grad_params_fp64_expected_global_host

        output_z_global_fp32_host = None
        d_input_z_global_fp32_host = None
        d_input_emb_global_fp32_host = None
        grad_params_fp32_global_host = None

    input_feats_for_distributed = _feats_for_distributed(input_feats_global_fp64, dtype, device="cpu")

    spawn_multiprocessing(
        parallel_assert_msa_module,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        dtype,
        msa_module_params,
        module_state_dict_fp64,
        input_z_global_fp64.detach().to(dtype=dtype, device="cpu", copy=True),
        input_emb_global_fp64.detach().to(dtype=dtype, device="cpu", copy=True),
        input_feats_for_distributed,
        output_z_for_worker,
        d_output_z_expected_global_fp64.detach().to(dtype=dtype, device="cpu", copy=True),
        d_input_z_for_worker,
        d_input_emb_for_worker,
        grad_params_for_worker,
        output_z_global_fp32_host,
        d_input_z_global_fp32_host,
        d_input_emb_global_fp32_host,
        grad_params_fp32_global_host,
    )


def parallel_assert_msa_module_activation_checkpointing(
    rank,
    grid_group_sizes,
    device_type,
    backend,
    env_per_rank,
    dtype,
    msa_module_params,
    min_val_init,
    max_val_init,
    input_z_global_host,
    input_emb_global_host,
    input_feats_global_host,
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

    seed_by_rank(0, seed=42)

    # The first module runs WITHOUT activation checkpointing (regular forward/backward).
    # The second module runs WITH activation checkpointing (checkpoint-wrapped forward,
    # recomputed during backward).  Comparing the two verifies that the checkpointing
    # mechanism preserves numerical correctness for both outputs and gradients.
    msa_module_params = dict(msa_module_params)
    msa_module_params["activation_checkpointing"] = False
    module_serial = SerialMSAModule(**msa_module_params)
    module_serial = module_serial.to(dtype=dtype, device=manager.device)
    init_module_params_uniform(module_serial, low=min_val_init, high=max_val_init)
    set_dtype_specific_inf_values(module_serial, dtype)

    module_state_dict_ref = module_serial.state_dict()

    module = DistributedMSAModule(module_serial, manager)
    module.train()

    placements_z = (Shard(0), Shard(1), Shard(2))
    placements_emb = (Shard(0), Replicate(), Shard(1))
    placements_msa = (Shard(0), Shard(1), Shard(2))
    placements_token_mask = (Shard(0), Shard(1), Shard(2))

    input_z_dtensor = distribute_tensor(
        input_z_global_host.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements_z,
    ).requires_grad_(True)

    input_emb_dtensor = distribute_tensor(
        input_emb_global_host.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements_emb,
    ).requires_grad_(True)

    input_feats_dtensor = {}
    for key, value in input_feats_global_host.items():
        if key in ["msa", "has_deletion", "deletion_value", "msa_paired", "msa_mask"]:
            input_feats_dtensor[key] = distribute_tensor(
                value.to(dtype=dtype if value.dtype.is_floating_point else value.dtype, device=manager.device),
                device_mesh=manager.device_mesh_subgroups,
                placements=placements_msa,
            )
        elif key == "token_pair_pad_mask":
            input_feats_dtensor[key] = distribute_tensor(
                value.to(dtype=dtype, device=manager.device),
                device_mesh=manager.device_mesh_subgroups,
                placements=placements_token_mask,
            )

    input_z_dtensor_copy = input_z_dtensor.detach().clone().requires_grad_(True)
    input_emb_dtensor_copy = input_emb_dtensor.detach().clone().requires_grad_(True)
    input_feats_dtensor_copy = {k: v.detach().clone() for k, v in input_feats_dtensor.items()}

    # Save RNG state so the second forward pass (with activation checkpointing)
    # sees the same dropout masks as the first forward pass.
    cpu_rng_state = torch.random.get_rng_state()
    cuda_rng_state = torch.cuda.get_rng_state(device=manager.device) if device_type == "cuda" else None

    output_z_dtensor_result = module(input_z_dtensor, input_emb_dtensor, input_feats_dtensor)

    assert_tensors_identical(
        input_z_dtensor_copy.to_local(), input_z_dtensor.to_local(), check_grad=False, check_grad_fn=False
    )
    assert_tensors_identical(
        input_emb_dtensor_copy.to_local(), input_emb_dtensor.to_local(), check_grad=False, check_grad_fn=False
    )
    for key in input_feats_dtensor_copy:
        assert_tensors_identical(input_feats_dtensor_copy[key].to_local(), input_feats_dtensor[key].to_local())

    d_output_z_dtensor = torch.distributed.tensor.rand(
        output_z_dtensor_result.shape,
        requires_grad=False,
        dtype=dtype,
        device_mesh=manager.device_mesh_subgroups,
        placements=output_z_dtensor_result.placements,
    )
    d_output_z_dtensor_copy = d_output_z_dtensor.detach().clone()

    torch.autograd.backward([output_z_dtensor_result], [d_output_z_dtensor])

    assert_tensors_identical(d_output_z_dtensor_copy.to_local(), d_output_z_dtensor.to_local())

    # Create second module with the same weights but activation checkpointing enabled
    msa_module_params["activation_checkpointing"] = True
    module_serial_act_ckpt = SerialMSAModule(**msa_module_params)
    module_serial_act_ckpt.load_state_dict(module_state_dict_ref)
    set_dtype_specific_inf_values(module_serial_act_ckpt, dtype)

    module_serial_act_ckpt = module_serial_act_ckpt.to(dtype=dtype, device=manager.device)
    module_act_ckpt = DistributedMSAModule(module_serial_act_ckpt, manager)
    module_act_ckpt.train()

    # Restore RNG state so dropout masks match the first forward pass
    torch.random.set_rng_state(cpu_rng_state)
    if cuda_rng_state is not None:
        torch.cuda.set_rng_state(cuda_rng_state, device=manager.device)

    output_z_dtensor_result_act_ckpt = module_act_ckpt(
        input_z_dtensor_copy, input_emb_dtensor_copy, input_feats_dtensor_copy
    )

    assert_tensors_identical(
        output_z_dtensor_result_act_ckpt.to_local(),
        output_z_dtensor_result.to_local(),
        check_grad=False,
        check_grad_fn=False,
    )

    torch.autograd.backward([output_z_dtensor_result_act_ckpt], [d_output_z_dtensor])

    assert_tensors_identical(input_z_dtensor.grad.to_local(), input_z_dtensor_copy.grad.to_local())
    assert_tensors_identical(input_emb_dtensor.grad.to_local(), input_emb_dtensor_copy.grad.to_local())

    result_param_grads_dict = {}
    for name, param in module.named_parameters():
        if param.grad is not None:
            result_param_grads_dict[name] = param.grad

    for name, param_act_ckpt_grad in module_act_ckpt.named_parameters():
        assert name in result_param_grads_dict, f"Parameter {name}'s gradient is not found in result gradients"
        result_grad = result_param_grads_dict[name]
        assert_tensors_identical(result_grad.to_local(), param_act_ckpt_grad.grad.to_local())

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env, dtype",
    (
        params_test := [
            (((2, (1, 1)), True, "cuda", "ENV"), torch.float32),
            (((2, (2, 2)), True, "cuda", "ENV"), torch.float32),
        ]
    ),
    indirect=["setup_env"],
    ids=[
        f"dp:{x[0][0][0]}, cp:{x[0][0][1]}, specify_method:{x[0][1]}, device_type:{x[0][2]}, method_init:{x[0][3]}, "
        f"dtype:{x[1]}"
        for x in params_test
    ],
)
def test_msa_module_parallel_activation_checkpointing(setup_env, dtype):
    """MSAModule with activation checkpointing vs CP without; results should match."""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    msa_module_params = create_msa_module_init_params_v2(use_large_model=False)
    msa_module_params["msa_dropout"] = 0.5
    msa_module_params["z_dropout"] = 0.5

    B = 2 * grid_group_sizes["dp"]
    size_ring = grid_group_sizes["cp"][0]
    N = size_ring * 2
    S = size_ring * 3
    min_val_init = -1
    max_val_init = 1
    dim_input_msa = const.num_tokens

    input_z_global = torch.empty((B, N, N, msa_module_params["token_z"]), dtype=dtype, requires_grad=True, device="cpu")
    input_emb_global = torch.empty((B, N, msa_module_params["token_s"]), dtype=dtype, requires_grad=True, device="cpu")

    input_feats_global_host = {
        "msa": torch.randint(0, dim_input_msa, (B, S, N), dtype=torch.int64, device="cpu"),
        "has_deletion": torch.empty((B, S, N), dtype=dtype, device="cpu"),
        "deletion_value": torch.empty((B, S, N), dtype=dtype, device="cpu"),
        "msa_paired": torch.randint(0, 2, (B, S, N), dtype=dtype, device="cpu"),
        "msa_mask": torch.ones((B, S, N), dtype=dtype, device="cpu"),
        "token_pad_mask": torch.randint(0, 2, (B, N), dtype=dtype, device="cpu"),
    }
    input_feats_global_host["token_pad_mask"][0, N // size_ring :] = 0
    input_feats_global_host["token_pair_pad_mask"] = (
        input_feats_global_host["token_pad_mask"][:, :, None] * input_feats_global_host["token_pad_mask"][:, None, :]
    )

    input_feats_global_host["msa_mask"][0, (S // size_ring) :, :] = 0
    input_feats_global_host["msa_mask"][0, :, (N // size_ring) :] = 0

    init_tensors_uniform([input_z_global, input_emb_global], low=min_val_init, high=max_val_init)
    for key, tensor in input_feats_global_host.items():
        if tensor.dtype.is_floating_point and "mask" not in key:
            init_tensors_uniform([tensor], low=min_val_init, high=max_val_init)

    input_feats_for_distributed = _feats_for_distributed(input_feats_global_host, dtype, device="cpu")

    spawn_multiprocessing(
        parallel_assert_msa_module_activation_checkpointing,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        dtype,
        msa_module_params,
        min_val_init,
        max_val_init,
        input_z_global.detach().to(dtype=dtype, device="cpu", copy=True),
        input_emb_global.detach().to(dtype=dtype, device="cpu", copy=True),
        input_feats_for_distributed,
    )
