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

"""Tests a single instance of AttentionPairBiasWithDTensor

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
    V11: multi-proc parameter gradident values identical across proc's

Implementation status
    V1:
    V2:
    V3: NA
    V4: implemented
    V5: implemented
    V6: implied by V3 and V9
    V7: same data
    V8: implemented
    V9: implemented
    V10: implemented
    V11:

Assertion threshold defaults for pytorch

dtype       rtol        atol
--------    -------     ------
float16     1e-3        1e-5
bfloat16    1.6e-2      1e-5
float32     1.3e-6      1e-5

"""

from collections import OrderedDict
from copy import deepcopy
from typing import Any, Union

import pytest
import torch
from torch import Tensor
from torch.distributed.tensor import DTensor, Replicate, Shard, distribute_tensor
from torch.testing import assert_close

from boltz.distributed.comm import AttentionPairBiasComm
from boltz.distributed.manager import DistributedManager, _GridGroupSizesType
from boltz.distributed.model.layers.attention import AttentionPairBias as AttentionPairBiasWithDTensor
from boltz.model.layers.attention import AttentionPairBias as AttentionPairBiasSerialV1
from boltz.model.layers.attentionv2 import AttentionPairBias as AttentionPairBiasSerialV2
from boltz.testing.utils import (
    assert_all_identical,
    assert_tensors_identical,
    seed_by_rank,
    skip_if_cuda_not_avail_or_device_count_less_than_word_size,
    spawn_multiprocessing,
)

SEED = 42


def assert_attention_pair_bias_with_dtensor_fw_bw(
    rank: int,
    input_example: OrderedDict[str, Tensor],
    output_ref: Tensor,
    output_grad_example: Tensor,
    input_grads_ref: OrderedDict[str, Tensor],
    c_s: int,
    c_z: int,
    parameter_grads_ref_as_tensors: OrderedDict[str, Tensor],
    num_heads: int,
    layer_state_dict: OrderedDict[str, Tensor],
    inf: float,
    grid_group_sizes: _GridGroupSizesType,
    device_type: str,
    backend: str,
    env_per_rank: dict[str, str],
    serial_version: str,
    apply_initial_norm: bool,
    use_model_cache: bool,
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
    # Setup comm objects
    # --------------------------------------------------------------
    ring_comm = AttentionPairBiasComm(
        process_group=manager.group["cp"],
        group_layout=manager.layout_subgroups["cp"],
        cp_axis_0_group=manager.subgroups["cp"][0],
        cp_axis_1_group=manager.subgroups["cp"][1],
    )
    # -------------------------------------------------------------
    # Move inputs and ref outputs to device
    # --------------------------------------------------------------
    input_example_device = OrderedDict(
        [
            (k, v.detach().to(manager.device) if isinstance(v, Tensor) else deepcopy(v))
            for (k, v) in input_example.items()
        ]
    )
    output_ref = output_ref.to(manager.device)

    inputs_grad_ref_device = OrderedDict(
        [
            (k, v.detach().to(manager.device) if isinstance(v, Tensor) else deepcopy(v))
            for (k, v) in input_grads_ref.items()
        ]
    )
    s_grad_ref = inputs_grad_ref_device.get("s", None)
    z_grad_ref = inputs_grad_ref_device.get("z", None)

    output_grad_example: Tensor = output_grad_example.detach().to(manager.device)

    # -------------------------------------------------------------
    # Create module to test using the appropriate serial version
    # --------------------------------------------------------------
    if serial_version == "v1":
        single_proc_module = AttentionPairBiasSerialV1(
            c_s=c_s,
            c_z=c_z,
            num_heads=num_heads,
            inf=inf,
            initial_norm=apply_initial_norm,
        )
    else:
        single_proc_module = AttentionPairBiasSerialV2(
            c_s=c_s,
            c_z=c_z,
            num_heads=num_heads,
            inf=inf,
            compute_pair_bias=True,
        )
    single_proc_module.load_state_dict(state_dict=layer_state_dict)
    single_proc_module = single_proc_module.train()
    single_proc_module = single_proc_module.to(manager.device)

    multi_proc_module = AttentionPairBiasWithDTensor(
        attn_pair_bias=single_proc_module,
        device_mesh=manager.device_mesh_subgroups,
        ring_comm=ring_comm,
        apply_initial_norm=apply_initial_norm,
        compute_pair_bias=True,  # PairFormer always computes pair bias
        use_model_cache=use_model_cache,
    )
    # -----------------------------------------------------
    # Create input DTensors
    #   s is on device, the whole example input
    #   mask is on device, the whole example input
    #   s_grad is on device, the whole reference tensor
    #   z_grad is on device
    # ----------------------------------------------------
    placements_for_single_rep_nonparam = (Shard(0), Shard(1), Replicate())
    placements_for_pair_rep_nonparam = (Shard(0), Shard(1), Shard(2))
    # Note: pair_mask is not used in Boltz-2 serial AttentionPairBias - only 1D mask is used
    input_meta: OrderedDict[str, dict] = OrderedDict(
        [
            ("s", dict(placements=placements_for_single_rep_nonparam, requires_grad=True)),  # noqa: C408
            ("mask", dict(placements=placements_for_single_rep_nonparam, requires_grad=False)),  # noqa: C408
            ("z", dict(placements=placements_for_pair_rep_nonparam, requires_grad=True)),  # noqa: C408
        ]
    )
    input_example_as_dtensors = OrderedDict()
    for name, meta in input_meta.items():
        input_example_as_dtensors[name] = distribute_tensor(
            input_example_device[name], manager.device_mesh_subgroups, meta["placements"]
        ).requires_grad_(meta["requires_grad"])

    output_ref_dt = distribute_tensor(
        output_ref,
        manager.device_mesh_subgroups,
        placements_for_single_rep_nonparam,
    ).requires_grad_(False)

    output_grad_example_dt = distribute_tensor(
        output_grad_example,
        manager.device_mesh_subgroups,
        placements_for_single_rep_nonparam,
    ).requires_grad_(False)

    # -------------------------------------------------
    # Run FW
    # -------------------------------------------------
    input_example_clone_as_dtensors = OrderedDict(
        [(k, v.detach().clone().requires_grad_(v.requires_grad)) for k, v in input_example_as_dtensors.items()]
    )
    output_actual_dt: DTensor = multi_proc_module(**input_example_as_dtensors)

    # -------------------------------------------------------
    # V8: multi-proc FW output tensor values close-to single-proc
    # ------------------------------------------------------
    assert (
        output_actual_dt.shape == output_ref_dt.shape
    ), f"Output shape mismatch: {output_actual_dt.shape} != {output_ref_dt.shape}"
    assert (
        output_actual_dt.stride() == output_ref_dt.stride()
    ), f"Output stride mismatch: {output_actual_dt.stride()} != {output_ref_dt.stride()}"
    assert_close(output_actual_dt.full_tensor(), output_ref_dt.full_tensor())

    # -------------------------------------------------
    # Run BW
    # -------------------------------------------------
    output_grad_example_clone_dt = (
        output_grad_example_dt.detach().clone().requires_grad_(output_grad_example_dt.requires_grad)
    )
    output_actual_clone_dt = output_actual_dt.detach().clone().requires_grad_(output_actual_dt.requires_grad)
    output_actual_dt.backward(output_grad_example_dt)

    # -------------------------------------------------------
    # V4: FW input tensor values unchanged by FW and BW
    # ---------------------------------------------------------
    for k in input_example_as_dtensors.keys():
        assert_tensors_identical(
            input_example_as_dtensors[k].full_tensor(),
            input_example_clone_as_dtensors[k].full_tensor(),
            check_grad=False,
            check_grad_fn=False,
        )
    # -----------------------------------------------------------
    # V5: BW input tensor values unchanged by BW
    # -------------------------------------------------------------
    assert_tensors_identical(
        output_grad_example_dt,
        output_grad_example_clone_dt,
        check_grad=False,
        check_grad_fn=False,
    )
    assert_tensors_identical(
        output_actual_dt,
        output_actual_clone_dt,
        check_grad=False,
        check_grad_fn=False,
    )
    # --------------------------------------------------------------------
    # V9: multi-proc FW input gradient values close-to single-proc
    #   - conduct the check by materializing the full distributed tensor
    #   - s_grad, z_grad on device
    # ---------------------------------------------------------------------
    s_grad_actual_full: Tensor = input_example_as_dtensors["s"].grad.full_tensor()
    assert (
        s_grad_actual_full.shape == s_grad_ref.shape
    ), f"Gradient shape mismatch: {s_grad_actual_full.shape} != {s_grad_ref.shape}"
    assert (
        s_grad_actual_full.stride() == s_grad_ref.stride()
    ), f"Gradient stride mismatch: {s_grad_actual_full.stride()} != {s_grad_ref.stride()}"
    assert_close(s_grad_actual_full, s_grad_ref)

    if z_grad_ref is not None:
        z_grad_actual_full: Tensor = input_example_as_dtensors["z"].grad.full_tensor()
        assert (
            z_grad_actual_full.shape == z_grad_ref.shape
        ), f"Gradient shape mismatch: {z_grad_actual_full.shape} != {z_grad_ref.shape}"
        assert (
            z_grad_actual_full.stride() == z_grad_ref.stride()
        ), f"Gradient stride mismatch: {z_grad_actual_full.stride()} != {z_grad_ref.stride()}"
        assert_close(z_grad_actual_full, z_grad_ref)

    # --------------------------------------------------------------------
    # V10: multi-proc parameter gradient values close-to single-proc
    #
    #  (1) Trigger reductions on DTensor gradients before evaluating assert
    # ---------------------------------------------------------------------
    param_grads_actual_as_tensors = OrderedDict()
    for name, param in multi_proc_module.named_parameters():
        if (param.grad is None) != (parameter_grads_ref_as_tensors[name] is None):
            raise ValueError(
                f"Inconsistent grad state for {name} on rank {rank}: "
                f"result grad is {param.grad is None}, "
                f"reference grad is {parameter_grads_ref_as_tensors[name] is None}"
            )
        param_grads_actual_as_tensors[name] = None if param.grad is None else param.grad.full_tensor()

    for name, grad_ref in parameter_grads_ref_as_tensors.items():
        if (grad_ref is None) != (param_grads_actual_as_tensors[name] is None):
            raise ValueError(
                f"Inconsistent grad state for {name} on rank {rank}: "
                f"result grad is {param_grads_actual_as_tensors[name] is None}, "
                f"reference grad is {grad_ref is None}"
            )
        grad_actual = param_grads_actual_as_tensors[name]
        assert grad_actual.shape == grad_ref.shape, f"Gradient shape mismatch: {grad_actual.shape} != {grad_ref.shape}"
        assert (
            grad_actual.stride() == grad_ref.stride()
        ), f"Gradient stride mismatch: {grad_actual.stride()} != {grad_ref.stride()}"
        assert_close(
            grad_actual.cpu(),
            grad_ref,
            msg=lambda msg: f"Rank {rank} {name} grad mismatch\n{msg}\ngot:{grad_actual}\nwant:{grad_ref}",
        )
        assert_all_identical(grad_actual, manager.group["cp"])

    DistributedManager.cleanup()
    monkeypatch.undo()


def get_example_input_and_reference_output(
    bs: int,
    N_tokens: int,
    c_s: int,
    c_z: int,
    num_heads: int,
    inf: float,
    serial_version: str,
    apply_initial_norm: bool,
    multiplicity: int = 1,
    seed: int = SEED,
) -> tuple[
    OrderedDict[str, Union[Tensor, Any]],
    Tensor,
    Tensor,
    OrderedDict[str, Union[Tensor, Any]],
    OrderedDict[str, Union[Tensor, Any]],
    OrderedDict[str, Union[Tensor, Any]],
]:
    # ----------------------------------------
    # (0) Check use-case requirements
    # ----------------------------------------
    if multiplicity != 1:
        raise ValueError("multiplicity must be 1 for this use-case")

    # ----------------------------------------
    # (1) Initialize RNG
    # ----------------------------------------
    seed_by_rank(0, seed=seed)

    # -------------------------------------
    # (2) Create example inputs on host
    # Both V1 and V2 pairformer use k_in=s (no to_keys transformation)
    # -------------------------------------
    s = torch.rand(size=(bs, N_tokens, c_s), dtype=torch.float, requires_grad=True)
    input_example: OrderedDict[str, Union[Tensor, Any]] = OrderedDict(
        [
            ("s", s),
            ("z", torch.rand(size=(bs, N_tokens, N_tokens, c_z), dtype=torch.float, requires_grad=True)),
            ("mask", torch.randint(0, 2, size=(bs, N_tokens), dtype=torch.float)),
            ("multiplicity", multiplicity),
        ]
    )
    input_example_clone: OrderedDict[str, Union[Tensor, int, None]] = OrderedDict(
        [(k, v.detach().clone() if isinstance(v, Tensor) else deepcopy(v)) for k, v in input_example.items()]
    )
    # ----------------------------------------------------------------------
    # (3) Create single proc module on cpu / host
    # Run FW, BW with reference module on example inputs
    # -----------------------------------------------------------------------
    if serial_version == "v1":
        module_ref = AttentionPairBiasSerialV1(
            c_s=c_s,
            c_z=c_z,
            num_heads=num_heads,
            inf=inf,
            initial_norm=apply_initial_norm,
        )
    else:
        module_ref = AttentionPairBiasSerialV2(
            c_s=c_s,
            c_z=c_z,
            num_heads=num_heads,
            inf=inf,
            compute_pair_bias=True,
        )
    module_ref.proj_o.reset_parameters()  # avoid zero initialization
    module_ref = module_ref.train()
    state_dict_ref = module_ref.state_dict()

    # For pairformer, k_in=s (queries equal keys)
    # V1: forward(s, z, mask, multiplicity) -- k_in defaults to s internally
    # V2: forward(s, z, mask, k_in=s, multiplicity)
    if serial_version == "v1":
        output_ref = module_ref(
            s=input_example["s"],
            z=input_example["z"],
            mask=input_example["mask"],
            multiplicity=input_example["multiplicity"],
        )
    else:
        output_ref = module_ref(
            s=input_example["s"],
            z=input_example["z"],
            mask=input_example["mask"],
            k_in=input_example["s"],  # k_in = s for pairformer
            multiplicity=input_example["multiplicity"],
        )
    output_grad_example = torch.rand_like(output_ref)
    output_grad_example_clone = output_grad_example.detach().clone()
    output_ref.backward(output_grad_example)

    input_grads_ref_as_tensors: OrderedDict[str, Union[Tensor, None]] = OrderedDict(
        [
            ("s", input_example["s"].grad),
            ("z", input_example["z"].grad),
        ]
    )

    # ----------------------------------------------------------------
    # V1: single-proc FW input tensor values unchanged by FW and BW
    # ----------------------------------------------------------------
    for k in input_example:
        assert_close(input_example[k], input_example_clone[k])

    # ---------------------------------------------------------------
    #   V2: single-proc BW input tensor values unchanged by BW
    # ---------------------------------------------------------------
    assert_close(output_grad_example, output_grad_example_clone)

    # --------------------------------------
    # (3) Get parameter gradients, on host
    # --------------------------------------
    parameter_grads_ref_as_tensors = OrderedDict()
    for name, parameter in module_ref.named_parameters():
        if parameter.grad is not None:
            parameter_grads_ref_as_tensors[name] = parameter.grad
        else:
            parameter_grads_ref_as_tensors[name] = None

    return (
        input_example_clone,
        output_ref.detach(),
        output_grad_example,
        input_grads_ref_as_tensors,
        parameter_grads_ref_as_tensors,
        state_dict_ref,
    )


@pytest.mark.parametrize(
    "setup_env",
    [
        ((2, (2, 2)), True, "cuda", "ENV"),
    ],
    indirect=("setup_env",),
    ids=lambda x: f"dp={x[0][0]}, cp={x[0][1]}, specify_method={x[1]}, device_type={x[2]}, method_init={x[3]}",
)
@pytest.mark.parametrize(
    "version_config",
    [
        # (serial_version, apply_initial_norm, use_model_cache)
        # PairFormer always uses compute_pair_bias=True
        ("v1", True, False),  # V1 PairFormer: initial_norm=True, no cache
        ("v2", False, False),  # V2 PairFormer: no initial norm, no cache
    ],
    ids=lambda x: f"serial:{x[0]}, init_norm:{x[1]}, cache:{x[2]}",
)
def test_attention_pair_bias_with_dtensor_for_pairformer_use_case(
    setup_env: dict[str, int],
    version_config: tuple[str, bool, bool],
    bs: int = 2,
    c_s: int = 2 * 5,
    num_heads: int = 5,
    c_z: int = 7,
    multiplicity: int = 1,
    seed: int = SEED,
):
    serial_version, apply_initial_norm, use_model_cache = version_config
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env
    skip_if_cuda_not_avail_or_device_count_less_than_word_size(device_type, world_size)
    N_tokens: int = grid_group_sizes["cp"][0] * 32
    inf = 1e6

    # ------------------------------------------------------------
    # (0) Check use-case requirements / implementation scope
    # ------------------------------------------------------------
    if multiplicity != 1:
        raise ValueError("multiplicity must be 1 for this use-case")

    # ----------------------------------------
    # (1) Get example inputs and reference outputs
    # ----------------------------------------
    (
        input_example,
        output_ref,
        output_grad_example,
        input_grads_ref_as_tensors,
        parameter_grads_ref_as_tensors,
        state_dict_ref,
    ) = get_example_input_and_reference_output(
        bs,
        N_tokens,
        c_s,
        c_z,
        num_heads,
        inf,
        serial_version,
        apply_initial_norm,
        multiplicity,
        seed=seed,
    )
    spawn_multiprocessing(
        assert_attention_pair_bias_with_dtensor_fw_bw,
        world_size,
        input_example,
        output_ref,
        output_grad_example,
        input_grads_ref_as_tensors,
        c_s,
        c_z,
        parameter_grads_ref_as_tensors,
        num_heads,
        state_dict_ref,
        inf,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        serial_version,
        apply_initial_norm,
        use_model_cache,
    )


if __name__ == "__main__":
    pytest.main([__file__])
