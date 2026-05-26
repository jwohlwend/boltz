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

"""Tests for DTensor-based CP RelativePositionEncoder.

This module tests the DTensor context-parallel RelativePositionEncoder against
both serial v1 and v2 implementations.

Verification checks:
    V1: single-proc FW input tensor values unchanged by FW and BW
    V2: single-proc BW input tensor values unchanged by BW
    V4a: multi-proc FW input tensor values unchanged by FW
    V4b: multi-proc FW input tensor values unchanged after BW
    V5: multi-proc BW input tensor (output_grad) values unchanged by BW
    V8: multi-proc FW output tensor values close-to single-proc
    V9: (N/A — integer inputs have no gradient)
    V10: multi-proc parameter gradient values close-to single-proc
    V10b: replicated parameter gradients identical across all CP ranks

bf16 coverage (CUDA-only):
    Both serial and distributed forward passes are wrapped in
    torch.autocast("cuda", dtype=torch.bfloat16).  The module weights
    remain fp32; autocast handles the downcast for matmul ops.  This
    verifies that outputs and parameter gradients match under mixed
    precision.
"""

from collections import OrderedDict
from contextlib import nullcontext

import pytest
import torch
from torch import Tensor
from torch.distributed.tensor import DTensor, Replicate, Shard, distribute_tensor
from torch.testing import assert_close

from boltz.distributed.comm import TransposeComm
from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.modules.encoders import (
    RelativePositionEncoder as DistributedRelativePositionEncoder,
)
from boltz.model.modules.encoders import (
    RelativePositionEncoder as SerialRelativePositionEncoderV1,
)
from boltz.model.modules.encodersv2 import (
    RelativePositionEncoder as SerialRelativePositionEncoderV2,
)
from boltz.testing.utils import (
    assert_all_identical,
    assert_tensors_identical,
    init_module_params_uniform,
    random_features,
    skip_if_cuda_not_avail_or_device_count_less_than_word_size,
    spawn_multiprocessing,
)

SEED = 42


def _assert_unchanged(actual, expected, *, serial=False):
    """Shorthand for assert_tensors_identical with standard immutability kwargs."""
    assert_tensors_identical(
        actual,
        expected,
        check_stride=True,
        check_grad=False,
        check_grad_fn=False,
        check_storage_pointer=False,
        check_storage_offset=serial,
    )


_RELPOS_KEYS = ["asym_id", "entity_id", "residue_index", "token_index", "sym_id", "cyclic_period"]


def _make_serial_module(
    serial_version: str,
    token_z: int,
    r_max: int,
    s_max: int,
    fix_sym_check: bool,
    cyclic_pos_enc: bool,
):
    """Construct the appropriate serial RelativePositionEncoder."""
    if serial_version == "v1":
        return SerialRelativePositionEncoderV1(token_z=token_z, r_max=r_max, s_max=s_max)
    return SerialRelativePositionEncoderV2(
        token_z=token_z,
        r_max=r_max,
        s_max=s_max,
        fix_sym_check=fix_sym_check,
        cyclic_pos_enc=cyclic_pos_enc,
    )


def _get_serial_reference(
    B: int,
    N: int,
    token_z: int,
    r_max: int,
    s_max: int,
    fix_sym_check: bool,
    cyclic_pos_enc: bool,
    serial_version: str,
    device: str = "cpu",
    use_autocast: bool = False,
    seed: int = SEED,
):
    """Run serial RelativePositionEncoder and collect all reference data.

    Returns inputs, outputs, gradients, and state_dict for distributed comparison.
    Also performs serial immutability checks.

    When use_autocast=True, the forward pass is wrapped in
    torch.autocast("cuda", dtype=torch.bfloat16) and the module weights
    stay fp32 (autocast handles downcast).
    """
    with torch.random.fork_rng(devices=[], enabled=True):
        torch.manual_seed(seed)

        serial = _make_serial_module(serial_version, token_z, r_max, s_max, fix_sym_check, cyclic_pos_enc)
        init_module_params_uniform(serial, low=-0.5, high=0.5)
        serial = serial.to(device=device)
        state_dict = {k: v.cpu().clone() for k, v in serial.state_dict().items()}

        feats = random_features(
            size_batch=B,
            n_tokens=N,
            n_atoms=N,
            n_msa=1,
            atom_counts_per_token_range=(1, 1),
            device=torch.device(device),
            float_value_range=(-1.0, 1.0),
            selected_keys=_RELPOS_KEYS,
        )

        # Clone inputs for V1 immutability check
        feats_copy = {k: v.detach().clone() for k, v in feats.items()}

        ac_ctx = torch.autocast("cuda", dtype=torch.bfloat16) if use_autocast else nullcontext()

        with ac_ctx:
            out_ref = serial(feats)

        # V1a: serial FW input tensor values unchanged
        for k in feats:
            _assert_unchanged(feats[k], feats_copy[k], serial=True)

        grad_out = torch.randn_like(out_ref)
        grad_out_copy = grad_out.detach().clone()

        out_ref.backward(grad_out)

        # V1b: serial FW input tensor values unchanged after backward
        for k in feats:
            _assert_unchanged(feats[k], feats_copy[k], serial=True)

        # V2: serial BW input tensor (grad_out) values unchanged
        _assert_unchanged(grad_out, grad_out_copy, serial=True)

        param_grads = OrderedDict()
        for n, p in serial.named_parameters():
            if p.grad is not None:
                param_grads[n] = p.grad.detach().cpu().clone()

    feats_host = {k: v.detach().cpu().clone() for k, v in feats.items()}

    return (
        feats_host,
        out_ref.detach().cpu(),
        grad_out.detach().cpu(),
        param_grads,
        state_dict,
    )


def _worker_relpos_parity(
    rank: int,
    feats_on_host: dict[str, Tensor],
    output_ref_on_host: Tensor,
    grad_output_on_host: Tensor,
    param_grads_ref: OrderedDict[str, Tensor | None],
    state_dict: dict,
    token_z: int,
    r_max: int,
    s_max: int,
    fix_sym_check: bool,
    cyclic_pos_enc: bool,
    serial_version: str,
    use_autocast: bool,
    grid_group_sizes: dict,
    device_type: str,
    backend: str,
    env_map: dict[str, str] | None = None,
):
    """Worker: compare distributed RelativePositionEncoder against serial reference.

    Performs V4a, V4b, V5, V8, V10, V10b checks.
    When use_autocast=True, wraps the distributed forward in
    torch.autocast("cuda", dtype=torch.bfloat16).
    """
    monkeypatch = pytest.MonkeyPatch()
    if env_map is not None:
        for var_name, value in env_map.items():
            if value == "<INPUT_RANK>":
                monkeypatch.setenv(var_name, f"{rank}")
                continue
            monkeypatch.setenv(var_name, value)

    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    dm = DistributedManager()

    serial = _make_serial_module(serial_version, token_z, r_max, s_max, fix_sym_check, cyclic_pos_enc)
    serial.load_state_dict(state_dict)
    serial = serial.to(device=dm.device)

    transpose_comm = TransposeComm(
        process_group=dm.group["cp"],
        group_layout=dm.layout_subgroups["cp"],
    )
    dist_mod = DistributedRelativePositionEncoder(serial, dm.device_mesh_subgroups, transpose_comm)
    dist_mod.train()

    # Distribute features: single representation placements (Shard(0), Shard(1), Replicate())
    single_placements = (Shard(0), Shard(1)) + (Replicate(),) * (dm.device_mesh_subgroups.ndim - 2)
    feats_dt = {}
    for key, val in sorted(feats_on_host.items()):
        feats_dt[key] = distribute_tensor(
            val.to(dm.device),
            dm.device_mesh_subgroups,
            single_placements,
        )

    # V4a setup: clone inputs for immutability check
    feats_dt_copy = {k: v.detach().clone().requires_grad_(v.requires_grad) for k, v in feats_dt.items()}

    ac_ctx = torch.autocast("cuda", dtype=torch.bfloat16) if use_autocast else nullcontext()

    # Forward
    with ac_ctx:
        out = dist_mod(feats_dt)

    # V4a: FW input tensor values unchanged by FW (binary identity)
    for key in feats_dt:
        _assert_unchanged(feats_dt[key].to_local(), feats_dt_copy[key].to_local())

    # V8: forward parity
    out_full = out.full_tensor()
    ref_on_device = output_ref_on_host.to(device=dm.device)
    # bf16 autocast: output dtype is bf16; ref was also produced under autocast
    fw_atol = 1e-2 if use_autocast else 1e-5
    fw_rtol = 1e-2 if use_autocast else 1.3e-6
    assert_close(
        out_full,
        ref_on_device,
        atol=fw_atol,
        rtol=fw_rtol,
        msg=lambda m: f"Rank {rank} forward output mismatch\n{m}",
    )

    # Backward setup
    pair_placements = (Shard(0), Shard(1), Shard(2))
    grad_out_dt = distribute_tensor(
        grad_output_on_host.to(device=dm.device),
        dm.device_mesh_subgroups,
        pair_placements,
    )

    # V5 setup: clone output and grad for immutability check
    out_clone = out.detach().clone().requires_grad_(out.requires_grad)
    grad_out_dt_clone = grad_out_dt.detach().clone().requires_grad_(grad_out_dt.requires_grad)

    out.backward(grad_out_dt)

    # V4b: FW input tensor values unchanged after backward
    for key in feats_dt:
        _assert_unchanged(feats_dt[key].to_local(), feats_dt_copy[key].to_local())

    # V5: BW input tensor values unchanged by BW
    _assert_unchanged(out.to_local(), out_clone.to_local())
    _assert_unchanged(grad_out_dt.to_local(), grad_out_dt_clone.to_local())

    # V10: parameter gradient parity
    # Tolerance rationale: the linear weight gradient is accumulated across
    # world_size ranks via reduce-sum, changing float32 accumulation order.
    # For N_accum~1376 elements and up to 8 ranks, the expected absolute
    # error is O(sqrt(N_accum) * N_ranks * eps_f32) ≈ 3.5e-05.
    # bf16 autocast: forward used bf16 intermediates so small values may
    # be flushed to zero.  Backward runs outside autocast (grads are fp32)
    # but inherits the bf16 rounding from saved tensors.  With N_ranks=4
    # (cp=2×2), partial-sum vs full-sum accumulation order differences on
    # bf16-rounded inputs can reach ~8×eps_bf16 ≈ 0.063 per element.
    grad_atol = 8e-2 if use_autocast else 5e-5
    grad_rtol = 8e-2 if use_autocast else 5e-5
    checked = 0
    for name, param in dist_mod.named_parameters():
        ref = param_grads_ref.get(name)
        if ref is None:
            continue
        assert param.grad is not None, f"Param {name} grad is None but serial had gradient"
        actual = param.grad.full_tensor() if isinstance(param.grad, DTensor) else param.grad
        assert_close(
            actual,
            ref.to(device=dm.device),
            atol=grad_atol,
            rtol=grad_rtol,
            msg=lambda m, n=name: f"Rank {rank} {n} grad mismatch\n{m}",
        )

        # V10b: replicated parameter gradients identical across all CP ranks
        grad_for_ident = actual.detach() if not isinstance(param.grad, DTensor) else param.grad.full_tensor().detach()
        assert_all_identical(grad_for_ident, dm.group["cp"])

        checked += 1

    assert checked >= 1, f"Expected at least 1 parameter gradient check, got {checked}"

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env, fix_sym_check, cyclic_pos_enc, serial_version, use_autocast",
    (
        _params := [
            (((2, (2, 2)), True, "cuda", "ENV"), False, False, "v2", False),
            (((2, (2, 2)), True, "cuda", "ENV"), True, True, "v2", True),
            (((2, (2, 2)), True, "cuda", "ENV"), False, False, "v1", False),
            (((2, (3, 3)), True, "cpu", "ENV"), False, False, "v2", False),
        ]
    ),
    indirect=("setup_env",),
    ids=[
        f"dp:{x[0][0][0]}, cp:{x[0][0][1]}, device:{x[0][2]}, "
        f"fix_sym:{x[1]}, cyclic:{x[2]}, ver:{x[3]}, autocast:{x[4]}"
        for x in _params
    ],
)
def test_dtensor_relpos_forward_backward(
    setup_env,
    fix_sym_check: bool,
    cyclic_pos_enc: bool,
    serial_version: str,
    use_autocast: bool,
):
    """RelativePositionEncoder: distributed output and gradients match serial reference.

    Covers fp32 parity (CPU+CUDA) and bf16 autocast precision (CUDA-only)
    across multiple mesh topologies, both v1 and v2 serial modules.
    """
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    skip_if_cuda_not_avail_or_device_count_less_than_word_size(
        device_type=device_type,
        world_size=world_size,
    )

    B = 2
    token_z = 32
    r_max = 8
    s_max = 2
    N = 8 * grid_group_sizes["cp"][0]

    (
        feats_host,
        out_ref,
        grad_out,
        param_grads,
        state_dict,
    ) = _get_serial_reference(
        B=B,
        N=N,
        token_z=token_z,
        r_max=r_max,
        s_max=s_max,
        fix_sym_check=fix_sym_check,
        cyclic_pos_enc=cyclic_pos_enc,
        serial_version=serial_version,
        device=device_type,
        use_autocast=use_autocast,
        seed=SEED,
    )

    spawn_multiprocessing(
        _worker_relpos_parity,
        world_size,
        feats_host,
        out_ref,
        grad_out,
        param_grads,
        state_dict,
        token_z,
        r_max,
        s_max,
        fix_sym_check,
        cyclic_pos_enc,
        serial_version,
        use_autocast,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
    )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
