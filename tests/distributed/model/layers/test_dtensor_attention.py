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


from functools import partial

import pytest
import torch
from torch.distributed.tensor import Replicate, Shard, distribute_tensor

from boltz.distributed.comm import AttentionPairBiasComm
from boltz.distributed.data.feature.featurizer_utils import get_pair_mask
from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.layers.attention import AttentionPairBias, AttentionPairBiasShardwise
from boltz.distributed.model.layers.utils import convert_single_repr_window_batched_query_to_key
from boltz.distributed.model.modules.utils import SDPAWithBiasBackend
from boltz.model.layers.attention import AttentionPairBias as SerialAttentionPairBiasV1
from boltz.model.layers.attentionv2 import AttentionPairBias as SerialAttentionPairBiasV2
from boltz.model.modules.encodersv2 import get_indexing_matrix, single_to_keys
from boltz.testing.utils import (
    assert_tensors_identical,
    get_to_keys,
    init_module_params_uniform,
    init_tensors_uniform,
    is_a6000_gpu,
    pair_global_to_window_batch,
    seed_by_rank,
    spawn_multiprocessing,
)


def assert_attention_pair_bias_for_atom_diffusion(
    rank: int,
    payload: tuple,
):
    (
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        dtype,
        multiplicity,
        sdpa_with_bias_backend,
        c_s,
        c_z,
        num_heads,
        inf,
        state_dict_reference,
        s_global_host_fp64,
        z_global_host_fp64,
        mask_global_host_fp64,
        pair_mask_global_host_fp64,
        o_global_host_fp64,
        d_o_global_host_fp64,
        d_s_expected_global_host_fp64,
        d_z_expected_global_host_fp64,
        grad_params_fp64_expected_global_host,
        serial_version,
        apply_initial_norm,
        compute_pair_bias,
        use_model_cache,
    ) = payload

    monkeypatch = pytest.MonkeyPatch()
    if env_per_rank is not None:
        for var_name, value in env_per_rank.items():
            if value == "<INPUT_RANK>":
                monkeypatch.setenv(var_name, f"{rank}")
                continue
            monkeypatch.setenv(var_name, value)

    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    manager = DistributedManager()

    if torch.finfo(dtype).resolution < torch.finfo(s_global_host_fp64.dtype).resolution:
        raise ValueError(
            f"Target dtype {dtype} has higher precision than reference output's dtype {s_global_host_fp64.dtype}"
        )

    seed_by_rank(rank)

    # create module and copy state dict using the appropriate serial version

    if serial_version == "v1":
        serial_module = SerialAttentionPairBiasV1(
            c_s=c_s,
            c_z=c_z,
            num_heads=num_heads,
            inf=inf,
            initial_norm=apply_initial_norm,
        )
    else:
        serial_module = SerialAttentionPairBiasV2(
            c_s=c_s,
            c_z=c_z if compute_pair_bias else None,
            num_heads=num_heads,
            inf=inf,
            compute_pair_bias=compute_pair_bias,
        )
    serial_module.load_state_dict(state_dict_reference)
    serial_module = serial_module.to(device=manager.device)

    ring_comm = AttentionPairBiasComm(
        process_group=manager.group["cp"],
        group_layout=manager.layout_subgroups["cp"],
        cp_axis_0_group=manager.subgroups["cp"][0],
        cp_axis_1_group=manager.subgroups["cp"][1],
    )
    module = AttentionPairBias(
        attn_pair_bias=serial_module,
        device_mesh=manager.device_mesh_subgroups,
        ring_comm=ring_comm,
        sdpa_with_bias_backend=sdpa_with_bias_backend,
        apply_initial_norm=apply_initial_norm,
        compute_pair_bias=compute_pair_bias,
        use_model_cache=use_model_cache,
    )
    module = module.to(device=manager.device, dtype=dtype)
    module = module.train()

    # Distribute input tensors
    placements_single = [Shard(0), Shard(1), Replicate()]
    placements_pair = [Shard(0), Shard(1), Shard(2)]

    s_dtensor = distribute_tensor(
        s_global_host_fp64.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements_single,
    ).requires_grad_(True)
    z_dtensor = distribute_tensor(
        z_global_host_fp64.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements_pair,
    ).requires_grad_(True)
    mask_dtensor = distribute_tensor(
        mask_global_host_fp64.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements_single,
    ).requires_grad_(False)
    pair_mask_dtensor = distribute_tensor(
        pair_mask_global_host_fp64.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements_pair,
    ).requires_grad_(False)

    # Distribute output gradient
    d_o_dtensor = distribute_tensor(
        d_o_global_host_fp64.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements_single,
    ).requires_grad_(False)

    # Create copies to verify inputs/upstream adjoint aren't modified
    s_dtensor_copy = s_dtensor.detach().clone().requires_grad_(True)
    z_dtensor_copy = z_dtensor.detach().clone().requires_grad_(True)
    mask_dtensor_copy = mask_dtensor.detach().clone().requires_grad_(False)
    pair_mask_dtensor_copy = pair_mask_dtensor.detach().clone().requires_grad_(False)
    d_o_dtensor_copy = d_o_dtensor.detach().clone().requires_grad_(False)

    # Forward pass
    o_dtensor = module(
        s=s_dtensor,
        z=z_dtensor,
        mask=mask_dtensor,
        pair_mask=pair_mask_dtensor,
        multiplicity=multiplicity,
    )

    # Verify inputs/upstream adjoint weren't modified
    assert_tensors_identical(s_dtensor_copy.to_local(), s_dtensor.to_local(), check_grad=False, check_grad_fn=False)
    assert_tensors_identical(z_dtensor_copy.to_local(), z_dtensor.to_local(), check_grad=False, check_grad_fn=False)
    assert_tensors_identical(
        mask_dtensor_copy.to_local(), mask_dtensor.to_local(), check_grad=False, check_grad_fn=False
    )
    assert_tensors_identical(
        pair_mask_dtensor_copy.to_local(), pair_mask_dtensor.to_local(), check_grad=False, check_grad_fn=False
    )
    assert_tensors_identical(d_o_dtensor_copy.to_local(), d_o_dtensor.to_local(), check_grad=False, check_grad_fn=False)

    # Test forward pass results
    assert (
        o_dtensor.shape == o_global_host_fp64.shape
    ), f"Output shape mismatch: {o_dtensor.shape} != {o_global_host_fp64.shape}"
    assert (
        o_dtensor.stride() == o_global_host_fp64.stride()
    ), f"Output stride mismatch: {o_dtensor.stride()} != {o_global_host_fp64.stride()}"
    atom_pad_mask_bool = mask_dtensor.full_tensor().bool()
    atom_pad_mask_bool_expanded = atom_pad_mask_bool
    if atom_pad_mask_bool.shape[0] != o_dtensor.shape[0]:
        atom_pad_mask_bool_expanded = atom_pad_mask_bool.repeat_interleave(multiplicity, 0)

    o_dtensor_full = o_dtensor.full_tensor()
    torch.testing.assert_close(
        (o_dtensor_full * atom_pad_mask_bool_expanded[:, :, None]).cpu(),
        (o_global_host_fp64 * atom_pad_mask_bool_expanded[:, :, None].cpu()).to(dtype=dtype),
    )

    # Backward pass
    o_dtensor.backward(d_o_dtensor)

    # Verify upstream gradient wasn't modified
    assert_tensors_identical(s_dtensor_copy.to_local(), s_dtensor.to_local(), check_grad=False, check_grad_fn=False)

    # Test input gradients
    s_inputs_dtensor_grad = s_dtensor.grad.full_tensor()
    torch.testing.assert_close(
        s_inputs_dtensor_grad[~atom_pad_mask_bool_expanded],
        torch.zeros_like(s_inputs_dtensor_grad[~atom_pad_mask_bool_expanded]),
    )
    torch.testing.assert_close(
        s_inputs_dtensor_grad[atom_pad_mask_bool_expanded].cpu(),
        d_s_expected_global_host_fp64[atom_pad_mask_bool_expanded.cpu()].to(dtype=dtype),
    )

    z_inputs_dtensor_grad = z_dtensor.grad.full_tensor()
    pair_mask_dtensor_full = pair_mask_dtensor.full_tensor().bool()

    # In broadcasting mode, atom_pad_mask_bool and pair_mask_dtensor_full
    # already have the original batch size (no need to undo repeat_interleave)
    if mask_dtensor.shape[0] == o_dtensor.shape[0]:
        pair_mask_dtensor_full_z = pair_mask_dtensor_full[::multiplicity]
    else:
        pair_mask_dtensor_full_z = pair_mask_dtensor_full

    # Test z gradient (window batching)
    bs, num_atoms = z_dtensor.shape[:2]
    z_inputs_dtensor_grad_reshaped = pair_global_to_window_batch(
        z_inputs_dtensor_grad,
        n_atoms_no_pads=torch.tensor([num_atoms] * bs, device=manager.device),
        pair_mask_global=pair_mask_dtensor_full_z[:, :, :, None],
    )
    torch.testing.assert_close(
        z_inputs_dtensor_grad_reshaped.cpu(),
        d_z_expected_global_host_fp64.to(dtype=dtype),
    )

    # Gather weight gradients using named_parameters
    result_param_grads_dict = {}
    for name, param in module.named_parameters():
        if param.grad is not None:
            if name not in grad_params_fp64_expected_global_host:
                raise ValueError(f"Parameter {name} has a resulting gradient but it is not in the reference module")
            result_param_grads_dict[name] = param.grad

    # Compare parameter gradients
    for name, expected_grad_global_host in grad_params_fp64_expected_global_host.items():
        assert name in result_param_grads_dict, f"Parameter {name}'s gradient is not found in result gradients"
        result_grad = result_param_grads_dict[name]
        torch.testing.assert_close(result_grad.full_tensor().cpu(), expected_grad_global_host.to(dtype=dtype))

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.slow
@pytest.mark.parametrize(
    "setup_env",
    (
        params_test := [
            ((2, (2, 2)), True, "cuda", "ENV"),
        ]
    ),
    indirect=["setup_env"],
    ids=[
        f"dp:{x[0][0]}, cp:{x[0][1]}, specify_method:{x[1]}, device_type:{x[2]}, method_init:{x[3]}"
        for x in params_test
    ],
)
@pytest.mark.parametrize(
    "config",
    [
        (1, False),
        (3, False),
        (3, True),
    ],
    ids=lambda x: f"multiplicity:{x[0]}, fix_window_batching:{x[1]}",
)
@pytest.mark.parametrize(
    "sdpa_with_bias_backend",
    [SDPAWithBiasBackend.REFERENCE, SDPAWithBiasBackend.TORCH_FLEX_ATTN],
    ids=lambda x: x.value,
)
@pytest.mark.parametrize(
    "version_config",
    [
        # (serial_version, apply_initial_norm, compute_pair_bias, use_model_cache)
        ("v1", False, True, True),  # V1 DTL: initial_norm=False, compute bias, cache z
        ("v2", False, False, False),  # V2 DTL: no init norm, pre-computed bias, no cache
    ],
    ids=lambda x: f"serial:{x[0]}, init_norm:{x[1]}, cpb:{x[2]}, cache:{x[3]}",
)
def test_attention_pair_bias_for_atom_diffusion(
    setup_env,
    config: tuple[int, bool],
    sdpa_with_bias_backend: SDPAWithBiasBackend,
    version_config: tuple[str, bool, bool, bool],
    dtype: torch.dtype = torch.float32,
    c_s: int = 16 * 2,
    c_z: int = 7,
    num_heads: int = 2,
    inf: float = 1e6,
):
    serial_version, apply_initial_norm, compute_pair_bias, use_model_cache = version_config
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env
    multiplicity, fix_window_batching = config

    if sdpa_with_bias_backend == SDPAWithBiasBackend.TORCH_FLEX_ATTN and device_type != "cuda":
        pytest.skip("torch_flex_attn requires cuda device")

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")
        if is_a6000_gpu() and world_size > 1:
            pytest.skip("skip cuda test because distribute_tensor leads to deadlock on A6000 GPUs")

    if multiplicity > 1 and not fix_window_batching:
        pytest.xfail(
            "There is a bug in Boltz1's code due to the difference in the order of repeat_interleave and view calls on single representations and pair bias. The context parallel version doesn't suffer from such bug and therefore won't produce consistent results with that of the serial Boltz1"
        )

    seed_by_rank(0)

    # Use the appropriate serial module for the version
    # When compute_pair_bias=False (V2 DTL), z last dim is num_heads (pre-computed bias)
    # When compute_pair_bias=True, z last dim is c_z (projected through LayerNorm+Linear)
    z_last_dim = c_z if compute_pair_bias else num_heads
    if serial_version == "v1":
        reference_module_fp64 = SerialAttentionPairBiasV1(
            c_s=c_s,
            c_z=c_z,
            num_heads=num_heads,
            inf=inf,
            initial_norm=apply_initial_norm,
        )
    else:
        reference_module_fp64 = SerialAttentionPairBiasV2(
            c_s=c_s,
            c_z=c_z if compute_pair_bias else None,
            num_heads=num_heads,
            inf=inf,
            compute_pair_bias=compute_pair_bias,
        )
    reference_module_fp64 = reference_module_fp64.to(device_type, dtype=torch.float64)
    reference_module_fp64 = reference_module_fp64.train()

    val_init_range = 0.25
    init_module_params_uniform(reference_module_fp64, low=-val_init_range, high=val_init_range)
    state_dict_reference = {k: v.detach().clone().cpu() for k, v in reference_module_fp64.state_dict().items()}

    # mock inputs and output
    bs = 2 * grid_group_sizes["dp"]
    num_atoms = 128 * 4  # multiple of 128 for window batching

    s = torch.empty(
        size=(bs * multiplicity, num_atoms, c_s),
        dtype=torch.float64,
        requires_grad=True,
        device=device_type,
    )  # repeat_interleave happens in AtomAttentionEncoder
    z = torch.empty(
        size=(bs, num_atoms, num_atoms, z_last_dim),
        dtype=torch.float64,
        requires_grad=False,  # z gradient not tested for window batching
        device=device_type,
    )  # repeat_interleave happens in AttentionPairBias
    mask = torch.ones(bs, num_atoms, dtype=torch.float64, device=device_type)
    mask[:, -5:] = 0  # insert padding at the end of the sequence
    pair_mask = get_pair_mask(num_atoms).to(dtype=torch.float64, device=device_type)
    pair_mask = pair_mask.unsqueeze(0).repeat(bs, 1, 1)

    init_tensors_uniform([s, z], low=-val_init_range, high=val_init_range)

    s_global_host_fp64 = s.detach().clone().cpu()
    z_global_host_fp64 = z.detach().clone().cpu()
    mask_global_host_fp64 = mask.detach().clone().cpu()
    pair_mask_global_host_fp64 = pair_mask.detach().clone().cpu()

    # Run serial forward pass (window batching)
    # reshape in AtomAttentionEncoder
    to_keys = get_to_keys(s)

    # reshape in AtomTransformer
    W, H = 32, 128
    B, N, D = s.shape
    NW = N // W

    s_reshaped = s.view((B * NW, W, -1))
    to_keys_new = lambda x: to_keys(x.view(B, NW * W, -1)).view(B * NW, H, -1)  # noqa: E731
    # In Boltz-2, mask is also transformed by to_keys to match key dimension (H=128)
    mask_reshaped = to_keys_new(mask.repeat_interleave(multiplicity, 0).unsqueeze(-1)).squeeze(-1)

    # remap pair representation from square to window shape
    z_wb = pair_global_to_window_batch(
        z,
        n_atoms_no_pads=torch.tensor([num_atoms] * bs, device=device_type),
        pair_mask_global=pair_mask[:, :, :, None],
    ).requires_grad_(True)
    z_reshaped = z_wb

    # reshape in AtomTransformer
    if fix_window_batching:
        # repeat_interleave -> view
        z_reshaped = z_reshaped.repeat_interleave(multiplicity, 0)
        z_reshaped = z_reshaped.view((B * NW, W, H, -1))
    else:
        # view -> repeat_interleave
        z_reshaped = z_reshaped.view((B * NW // multiplicity, W, H, -1))

    # AttentionPairBias forward pass
    # V1: uses to_keys internally (mask must be query-aligned, module transforms it)
    # V2: uses pre-computed k_in (mask already key-aligned)
    if serial_version == "v1":
        # V1 expects query-aligned mask (B*NW, W); to_keys is applied inside forward
        mask_query = mask.repeat_interleave(multiplicity, 0).view(B * NW, W)
        o_attn_global_fp64 = reference_module_fp64(
            s=s_reshaped,
            z=z_reshaped,
            mask=mask_query,
            to_keys=to_keys_new,
            multiplicity=1 if fix_window_batching else multiplicity,
        )
    else:
        k_in_reshaped = to_keys_new(s_reshaped)
        o_attn_global_fp64 = reference_module_fp64(
            s=s_reshaped,
            z=z_reshaped,
            mask=mask_reshaped,
            k_in=k_in_reshaped,
            multiplicity=1 if fix_window_batching else multiplicity,
        )

    # reshape in AtomTransformer
    o_global_fp64 = o_attn_global_fp64.view((B, NW * W, D))

    o_global_host_fp64 = o_global_fp64.detach().clone().cpu()

    # Create upstream gradients, apply masks, and run backward pass
    d_o_global_fp64 = torch.empty_like(o_global_fp64)  # (B, N, D)
    init_tensors_uniform([d_o_global_fp64], low=-val_init_range, high=val_init_range)
    d_o_global_fp64 = d_o_global_fp64 * mask[:, :, None].repeat_interleave(multiplicity, 0)
    d_o_global_host_fp64 = d_o_global_fp64.detach().clone().cpu()

    o_global_fp64.backward(d_o_global_fp64)

    grad_params_fp64_expected_global_host = {
        k: v.grad.detach().clone().cpu() for k, v in reference_module_fp64.named_parameters() if v.grad is not None
    }

    # Get reference input gradients
    d_s_expected_global_host_fp64 = s.grad.detach().clone().cpu()
    d_z_expected_global_host_fp64 = z_wb.grad.detach().clone().cpu()

    payload = (
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        dtype,
        multiplicity,
        sdpa_with_bias_backend,
        c_s,
        c_z,
        num_heads,
        inf,
        state_dict_reference,
        s_global_host_fp64,
        z_global_host_fp64,
        mask_global_host_fp64,
        pair_mask_global_host_fp64,
        o_global_host_fp64,
        d_o_global_host_fp64,
        d_s_expected_global_host_fp64,
        d_z_expected_global_host_fp64,
        grad_params_fp64_expected_global_host,
        # Version config (last 4 elements, extracted by parallel function)
        serial_version,
        apply_initial_norm,
        compute_pair_bias,
        use_model_cache,
    )

    spawn_multiprocessing(assert_attention_pair_bias_for_atom_diffusion, world_size, payload)


def parallel_assert_shardwise_attention_pair_bias(
    rank: int,
    grid_group_sizes: dict[str, int],
    device_type: str,
    backend: str,
    env_map: dict[str, str],
    dtype: torch.dtype,
    sdpa_with_bias_backend: SDPAWithBiasBackend,
    reference_state_dict: dict,
    c_s: int,
    c_z: int,
    num_heads: int,
    inf: float,
    s_global_host: torch.Tensor,
    z_global_host: torch.Tensor,
    mask_global_host: torch.Tensor,
    k_in_global_host: torch.Tensor,  # V2 API: pre-computed k_in
    o_global_host: torch.Tensor,
    d_o_global_host: torch.Tensor,
    d_s_expected_global_host: torch.Tensor,
    d_z_expected_global_host: torch.Tensor,
    d_k_in_expected_global_host: torch.Tensor,  # V2 API: k_in gradient
    grad_params_expected_global_host: dict[str, torch.Tensor],
    serial_version: str,
    apply_initial_norm: bool,
    compute_pair_bias: bool,
    use_model_cache: bool,
):
    monkeypatch = pytest.MonkeyPatch()
    if env_map is not None:
        for var_name, value in env_map.items():
            if value == "<INPUT_RANK>":
                monkeypatch.setenv(var_name, f"{rank}")
                continue
            monkeypatch.setenv(var_name, value)

    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    manager = DistributedManager()
    device = manager.device
    device_mesh = manager.device_mesh_subgroups

    seed_by_rank(0, 42)

    B, K, W, D = s_global_host.shape

    # Setup module using the appropriate serial version
    if serial_version == "v1":
        serial_apb_module = SerialAttentionPairBiasV1(
            c_s=c_s,
            c_z=c_z,
            num_heads=num_heads,
            inf=inf,
            initial_norm=apply_initial_norm,
        )
    else:
        serial_apb_module = SerialAttentionPairBiasV2(
            c_s=c_s,
            c_z=c_z if compute_pair_bias else None,
            num_heads=num_heads,
            inf=inf,
            compute_pair_bias=compute_pair_bias,
        )

    serial_apb_module.load_state_dict(reference_state_dict)
    serial_apb_module = serial_apb_module.to(device=manager.device)

    module = AttentionPairBiasShardwise(
        attn_pair_bias=serial_apb_module,
        device_mesh=device_mesh,
        sdpa_with_bias_backend=sdpa_with_bias_backend,
        apply_initial_norm=apply_initial_norm,
        compute_pair_bias=compute_pair_bias,
        use_model_cache=use_model_cache,
    )
    module = module.to(device=device, dtype=dtype)
    module = module.train()

    # NOTE: only need single rep placements because the "pair" is just (K=N//W, W=32, H=128), and K is sharded along CP0
    placements = (Shard(0), Shard(1), Replicate())

    # Shard the inputs
    s_dtensor = distribute_tensor(s_global_host.to(dtype=dtype, device=device), device_mesh, placements).requires_grad_(
        True
    )
    z_dtensor = distribute_tensor(z_global_host.to(dtype=dtype, device=device), device_mesh, placements).requires_grad_(
        True
    )
    mask_dtensor = distribute_tensor(
        mask_global_host.to(dtype=dtype, device=device), device_mesh, placements
    ).requires_grad_(False)
    d_o_dtensor = distribute_tensor(
        d_o_global_host.to(dtype=dtype, device=device), device_mesh, placements
    ).requires_grad_(False)

    # Create copies to verify inputs/upstream adjoint aren't modified
    s_dtensor_copy = s_dtensor.detach().clone().requires_grad_(True)
    z_dtensor_copy = z_dtensor.detach().clone().requires_grad_(True)
    mask_dtensor_copy = mask_dtensor.detach().clone().requires_grad_(False)
    d_o_dtensor_copy = d_o_dtensor.detach().clone().requires_grad_(False)

    if serial_version == "v1":
        # V1: use to_keys, mask is query-aligned (B, K, W)
        to_keys_dt = partial(convert_single_repr_window_batched_query_to_key, W=W, H=z_global_host.shape[3])
        o_dtensor = module(s_dtensor, z_dtensor, mask_dtensor, to_keys=to_keys_dt)
        k_in_dtensor = None
        k_in_dtensor_copy = None
    else:
        # V2: use pre-computed k_in, mask is key-aligned (B, K, H)
        k_in_dtensor = distribute_tensor(
            k_in_global_host.to(dtype=dtype, device=device), device_mesh, placements
        ).requires_grad_(True)
        k_in_dtensor_copy = k_in_dtensor.detach().clone().requires_grad_(True)
        o_dtensor = module(s_dtensor, z_dtensor, mask_dtensor, k_in=k_in_dtensor)

    # Verify inputs/upstream adjoint weren't modified
    assert_tensors_identical(s_dtensor_copy.to_local(), s_dtensor.to_local(), check_grad=False, check_grad_fn=False)
    assert_tensors_identical(z_dtensor_copy.to_local(), z_dtensor.to_local(), check_grad=False, check_grad_fn=False)
    assert_tensors_identical(
        mask_dtensor_copy.to_local(), mask_dtensor.to_local(), check_grad=False, check_grad_fn=False
    )
    if k_in_dtensor is not None:
        assert_tensors_identical(
            k_in_dtensor_copy.to_local(), k_in_dtensor.to_local(), check_grad=False, check_grad_fn=False
        )
    assert_tensors_identical(d_o_dtensor_copy.to_local(), d_o_dtensor.to_local(), check_grad=False, check_grad_fn=False)

    # Verify forward pass results
    assert (
        o_dtensor.stride() == o_global_host.stride()
    ), f"Output stride mismatch: {o_dtensor.stride()} != {o_global_host.stride()}"

    o_dtensor_full = o_dtensor.full_tensor()

    torch.testing.assert_close(o_dtensor_full.cpu(), o_global_host.cpu().to(dtype=dtype))

    # Run backward pass of distributed shardwise module
    o_dtensor.backward(d_o_dtensor)

    # Verify upstream input wasn't modified
    assert_tensors_identical(s_dtensor_copy.to_local(), s_dtensor.to_local(), check_grad=False, check_grad_fn=False)
    assert_tensors_identical(z_dtensor_copy.to_local(), z_dtensor.to_local(), check_grad=False, check_grad_fn=False)
    if k_in_dtensor is not None:
        assert_tensors_identical(
            k_in_dtensor_copy.to_local(), k_in_dtensor.to_local(), check_grad=False, check_grad_fn=False
        )

    # Verify input gradients
    s_inputs_dtensor_grad = s_dtensor.grad.full_tensor()
    torch.testing.assert_close(
        s_inputs_dtensor_grad.cpu(),
        d_s_expected_global_host.to(dtype=dtype),
    )

    z_inputs_dtensor_grad = z_dtensor.grad.full_tensor()
    torch.testing.assert_close(
        z_inputs_dtensor_grad.cpu(),
        d_z_expected_global_host.to(dtype=dtype),
    )

    # V2 API: Verify k_in gradient (only when k_in is used as separate input)
    if k_in_dtensor is not None and d_k_in_expected_global_host is not None:
        k_in_inputs_dtensor_grad = k_in_dtensor.grad.full_tensor()
        torch.testing.assert_close(
            k_in_inputs_dtensor_grad.cpu(),
            d_k_in_expected_global_host.to(dtype=dtype),
        )

    # Verify parameter gradients
    result_param_grads_dict = {}
    for name, param in module.named_parameters():
        if param.grad is not None:
            if name not in grad_params_expected_global_host:
                raise ValueError(f"Parameter {name} has a resulting gradient but it is not in the reference module")
            result_param_grads_dict[name] = param.grad

    # Compare parameter gradients
    for name, expected_grad_global_host in grad_params_expected_global_host.items():
        assert name in result_param_grads_dict, f"Parameter {name}'s gradient is not found in result gradients"
        result_grad = result_param_grads_dict[name]
        torch.testing.assert_close(result_grad.full_tensor().cpu(), expected_grad_global_host.to(dtype=dtype))

    # clean up
    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env",
    [
        ((2, (2, 2)), True, "cuda", "ENV"),
    ],
    indirect=("setup_env",),
    ids=lambda x: f"dp:{x[0][0]}, cp:{x[0][1]}, specify_method:{x[1]}, device_type:{x[2]}, method_init:{x[3]}",
)
@pytest.mark.parametrize(
    "sdpa_with_bias_backend",
    [
        SDPAWithBiasBackend.REFERENCE,
        SDPAWithBiasBackend.TORCH_SDPA_EFFICIENT_ATTENTION,
        SDPAWithBiasBackend.TORCH_FLEX_ATTN,
    ],
    ids=lambda x: x.value,
)
@pytest.mark.parametrize(
    "version_config",
    [
        # (serial_version, apply_initial_norm, compute_pair_bias, use_model_cache)
        ("v1", False, True, False),  # V1 DTL: initial_norm=False, compute bias, no cache
        ("v2", False, False, False),  # V2 DTL: no init norm, pre-computed bias, no cache
    ],
    ids=lambda x: f"serial:{x[0]}, init_norm:{x[1]}, cpb:{x[2]}, cache:{x[3]}",
)
def test_shardwise_attention_pair_bias(
    setup_env,
    sdpa_with_bias_backend: SDPAWithBiasBackend,
    version_config: tuple[str, bool, bool, bool],
):
    """Test shardwise attention with V1 and V2 serial modules (pre-computed k_in)."""
    serial_version, apply_initial_norm, compute_pair_bias, use_model_cache = version_config
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    seed_by_rank(0)

    dtype: torch.dtype = torch.float32
    c_s: int = 32  # c_s // num_heads must be at least 16 and must be multiple of 4 to test the kernels
    c_z: int = 32
    num_heads: int = 2
    inf: float = 1e6
    # When compute_pair_bias=False (V2 DTL), z last dim is num_heads (pre-computed bias)
    # When compute_pair_bias=True, z last dim is c_z (projected through LayerNorm+Linear)
    z_last_dim = c_z if compute_pair_bias else num_heads
    # mock inputs and output
    B = 2 * grid_group_sizes["dp"]
    W = 32
    H = 128
    num_atoms = 128 * 4  # multiple of 128 for window batching
    K = num_atoms // W  # number of windows = 16
    D = c_s

    # S is reshaped inside atom transformer already into (B, K, W, D)
    s = torch.empty(
        size=(B, K, W, D),
        dtype=dtype,
        requires_grad=True,
        device=device_type,
    )
    z = torch.empty(
        size=(B, K, W, H, z_last_dim),
        dtype=dtype,
        requires_grad=True,
        device=device_type,
    )

    # V2 API: mask will be key-aligned (B, K, H) after transformation
    mask_query_aligned = torch.randint(0, 2, s.shape[:-1], dtype=torch.float, device=device_type, requires_grad=False)
    d_o = torch.empty(
        size=(B, K, W, D),
        dtype=dtype,
        requires_grad=False,
        device=device_type,
    )

    val_init_range = 0.2
    init_tensors_uniform([s, z, d_o], low=-val_init_range, high=val_init_range)

    # mask gradients where the inputs are masked out. We have to do this post-sort because distributed version is done upstream
    d_o = d_o * mask_query_aligned.unsqueeze(-1)

    # reference module using the appropriate serial version
    if serial_version == "v1":
        reference_module = SerialAttentionPairBiasV1(
            c_s=c_s,
            c_z=c_z,
            num_heads=num_heads,
            inf=inf,
            initial_norm=apply_initial_norm,
        )
    else:
        reference_module = SerialAttentionPairBiasV2(
            c_s=c_s,
            c_z=c_z if compute_pair_bias else None,
            num_heads=num_heads,
            inf=inf,
            compute_pair_bias=compute_pair_bias,
        )
    reference_module = reference_module.to(device_type, dtype=dtype)
    reference_module = reference_module.train()

    init_module_params_uniform(reference_module, low=-val_init_range, high=val_init_range)

    reference_state_dict = {k: v.detach().clone().cpu() for k, v in reference_module.state_dict().items()}

    s_reshaped = s.view((B * K, W, -1))  # Q needs to be in this shape

    # Define single device to_keys function
    # This to_keys function assumes that input comes in as shape (B * K, W, D), so translate to (B, N, D)
    def _serial_to_keys(s: torch.Tensor, B: int, K: int, W: int, H: int, D: int) -> torch.Tensor:
        s = s.view(B, K * W, -1)
        indexing_matrix = get_indexing_matrix(K, W, H, s.device).to(dtype=s.dtype)
        return single_to_keys(s, indexing_matrix, W, H)

    to_keys_new = partial(_serial_to_keys, B=B, K=K, W=W, H=H, D=D)

    def _to_keys_new_reshape(x: torch.Tensor) -> torch.Tensor:
        return to_keys_new(x).view(B * K, H, -1)

    # V2 API: Pre-compute k_in and key-aligned mask
    # IMPORTANT: In V2 API, k_in is a separate input (not derived from s in the attention module).
    # To match distributed behavior where s and k_in are independent inputs, we:
    # 1. Compute k_in values from s (for numerical correctness)
    # 2. Detach and create a new leaf tensor for k_in (so s.grad only has Q-path gradients)
    k_in_values = _to_keys_new_reshape(s_reshaped.detach())  # Compute values without gradient connection
    k_in_reshaped = k_in_values.clone().requires_grad_(True)  # Create leaf tensor for gradient tracking
    mask_reshaped = _to_keys_new_reshape(mask_query_aligned.unsqueeze(-1)).squeeze(-1)

    z_reshaped = z.view((B * K, W, H, -1))

    # Run serial forward with the appropriate API
    if serial_version == "v1":
        # V1: pass to_keys, mask is query-aligned (B*K, W); module transforms internally
        mask_query = mask_query_aligned.view(B * K, W)
        o_serial = reference_module(
            s=s_reshaped,
            z=z_reshaped,
            mask=mask_query,
            to_keys=_to_keys_new_reshape,
        )
    else:
        # V2: pass pre-computed k_in, mask is key-aligned (B*K, H)
        o_serial = reference_module(
            s=s_reshaped,
            z=z_reshaped,
            mask=mask_reshaped,
            k_in=k_in_reshaped,
        )

    # clone forward pass output and match distributed module shape
    o_global_host = o_serial.detach().clone().cpu().view(B, K, W, D)

    d_o = d_o.view(B * K, W, D)

    o_serial.backward(d_o)

    # parameter gradients. The serial version has S in shape B, N, D
    d_s_expected_global_host = s.grad.detach().clone().cpu().view(B, K, W, D)
    d_z_expected_global_host = z.grad.detach().clone().cpu()

    if serial_version == "v1":
        # V1: k_in is computed internally from to_keys(s), no separate k_in gradient
        d_k_in_expected_global_host = None
        # V1: mask is query-aligned for distributed test (B, K, W)
        mask_global_host = mask_query_aligned.detach().clone().cpu()
        k_in_global_host = None
    else:
        # V2: k_in gradient
        d_k_in_expected_global_host = k_in_reshaped.grad.detach().clone().cpu().view(B, K, H, D)
        # V2: mask is key-aligned for distributed test (B, K, H)
        mask_global_host = mask_reshaped.detach().clone().cpu().view(B, K, H)
        # V2: k_in for distributed test (B, K, H, D)
        k_in_global_host = k_in_reshaped.detach().clone().cpu().view(B, K, H, D)

    s_global_host = s.detach().clone().cpu()
    z_global_host = z.detach().clone().cpu()
    d_o_global_host = d_o.detach().clone().cpu().view(B, K, W, D)

    grad_params_expected_global_host = {
        k: v.grad.detach().clone().cpu() for k, v in reference_module.named_parameters() if v.grad is not None
    }

    spawn_multiprocessing(
        parallel_assert_shardwise_attention_pair_bias,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        dtype,
        sdpa_with_bias_backend,
        reference_state_dict,
        c_s,
        c_z,
        num_heads,
        inf,
        s_global_host,
        z_global_host,
        mask_global_host,
        k_in_global_host,  # V2 API: pre-computed k_in
        o_global_host,
        d_o_global_host,
        d_s_expected_global_host,
        d_z_expected_global_host,
        d_k_in_expected_global_host,  # V2 API: k_in gradient
        grad_params_expected_global_host,
        serial_version,
        apply_initial_norm,
        compute_pair_bias,
        use_model_cache,
    )
