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
from torch.distributed.tensor import DTensor, Replicate, Shard, distribute_tensor

from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.layers.flatten_and_unflatten import shardwise_unflatten_sharded
from boltz.distributed.model.layers.utils import (
    convert_single_repr_to_window_batched_key,
    distributed_gather_sliding_windows,
    distributed_pack_and_pad,
    distributed_unpad_and_unpack,
    gather_sliding_windows,
    pack_and_pad,
)
from boltz.model.modules.encoders import get_indexing_matrix, single_to_keys
from boltz.testing.utils import assert_tensors_identical, seed_by_rank, spawn_multiprocessing


def parallel_assert_gather_sliding_windows(rank, grid_group_sizes, device_type, backend, env_map):
    """Run distributed version on each rank."""

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
    device_mesh = manager.device_mesh

    seed_by_rank(0, 42)

    for W, H, K_per_rank in [(8, 24, 3), (32, 128, 8), (32, 128, 1)]:
        assert H % (W // 2) == 0, f"H must be divisible by W // 2 but got {H} % {W // 2} != 0"
        h = H // (W // 2)
        assert h % 2 == 0, f"h := H // (W // 2) must be divisible by 2 but got {h} % 2 != 0"
        K = K_per_rank * grid_group_sizes["cp"]
        for n_per_rank in [4, 100]:
            # always shard the (2 * K, ) axis by the "cp" submesh
            for input_shape_extra in [(None,), (4, None), (None, 3), (2, None, 3)]:
                # The "None" element indicates the axis to be used as "half-window" axis
                # there can be only one "None" element in the input_shape
                if sum(1 for x in input_shape_extra if x is None) != 1:
                    raise ValueError(
                        f"There can be one and only one 'None' element in the input_shape but got {input_shape_extra}"
                    )
                axis = input_shape_extra.index(None)
                input_shape = input_shape_extra[:axis] + (2 * K, W // 2) + input_shape_extra[axis + 1 :]

                label = f"W:{W}, H:{H}, K:{K}, input_shape:{input_shape}"

                # Shard the leading axis if they exist and evenly divisible by the mesh shape
                # otherwise replicate
                ndim = len(input_shape)
                if ndim > 2 and axis != 0 and device_mesh.ndim >= 2 and input_shape[0] % device_mesh.size(0) == 0:
                    placements = (Shard(0),) + (Replicate(),) * (device_mesh.ndim - 2) + (Shard(axis),)
                else:
                    placements = (Replicate(),) * (device_mesh.ndim - 1) + (Shard(axis),)

                input_global = torch.randn(input_shape, dtype=torch.float32, device=device, requires_grad=True)

                offset_start = 1 - h // 2
                offsets = torch.arange(offset_start, offset_start + 2 * K, 2, device=device)

                output_ref = gather_sliding_windows(input_global, offsets, h, axis)

                # Backward
                grad_output_ref = torch.randn_like(output_ref)
                output_ref.backward(grad_output_ref)

                # Create sharded DTensor
                input_dtensor = distribute_tensor(
                    input_global.detach().clone(), device_mesh, placements
                ).requires_grad_(True)

                # Distributed forward
                output_dtensor = distributed_gather_sliding_windows(input_dtensor, h, axis)

                # Verify forward
                output_result_global = output_dtensor.full_tensor()
                # due the complicated reshaping involved and the potential concatenation of halo
                # between ranks, it's difficult to guarantee identical strides as in the global
                # case, which is not very useful in practice, so we just check the values
                # Also, compute_global_expectation returns detached tensors which voids grad existence check
                assert_tensors_identical(
                    output_result_global,
                    output_ref,
                    check_stride=False,
                    check_grad=False,
                    check_grad_fn=False,
                    msg=lambda m: f"{label} fwd output mismatch:\n  {m}",
                )

                # Distributed backward
                grad_output_dtensor = distribute_tensor(
                    grad_output_ref.detach().clone(), device_mesh, output_dtensor.placements
                )
                output_dtensor.backward(grad_output_dtensor)

                # Verify backward (with tolerance for numerical precision)
                grad_input_result_global = input_dtensor.grad.full_tensor()

                torch.testing.assert_close(
                    grad_input_result_global,
                    input_global.grad,
                    msg=lambda m: f"{label} input gradient mismatch:\n  {m}",
                )

    DistributedManager.cleanup()


@pytest.mark.parametrize(
    "setup_env",
    [
        ((1, 4), True, "cuda", "ENV"),
        ((1, 7), True, "cuda", "ENV"),
        ((2, 4), True, "cuda", "ENV"),
    ],
    indirect=("setup_env",),
    ids=lambda x: f"dp:{x[0][0]}, cp:{x[0][1]}, specify_method:{x[1]}, device_type:{x[2]}, method_init:{x[3]}",
)
def test_distributed_gather_sliding_windows(setup_env):
    """Test distributed gather sliding windows"""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    spawn_multiprocessing(
        parallel_assert_gather_sliding_windows,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
    )


def parallel_assert_pack_and_pad(rank, grid_group_sizes, device_type, backend, env_map):
    """Run distributed pack and pad test."""

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
    device_mesh = manager.device_mesh

    # 1. Create global data
    seed_by_rank(0, 42)

    for W in [8, 32]:
        for n_per_rank in [4, 100]:
            # For testing purpose, always shard the "axis" along last dimension of the device mesh
            n_global = n_per_rank * device_mesh.size(-1)
            for input_shape_extra in [(None,), (4, None), (None, 3), (2, None, 3)]:
                if sum(1 for x in input_shape_extra if x is None) != 1:
                    raise ValueError(
                        f"There can be one and only one 'None' element in the input_shape but got {input_shape_extra}"
                    )
                axis = input_shape_extra.index(None)
                ndim = len(input_shape_extra)
                input_shape = input_shape_extra[:axis] + (n_global,) + input_shape_extra[axis + 1 :]

                axes_mask_broadcast = torch.arange(ndim).tolist()
                for axis_mask_broadcast in axes_mask_broadcast:
                    if axis_mask_broadcast == axis:
                        # mask has same shape as input
                        shape_mask = input_shape
                    else:
                        # mask has shape 1 along the broadcast axis
                        shape_mask = input_shape[:axis_mask_broadcast] + (1,) + input_shape[axis_mask_broadcast + 1 :]

                    for keep_input_padding in [False, True]:
                        label = f"W:{W}, input_shape:{input_shape}, shape_mask:{shape_mask}, keep_input_padding:{keep_input_padding}"

                        input_global = torch.randn(input_shape, requires_grad=True, device=device)
                        mask_global = torch.randint(0, 2, shape_mask, dtype=torch.bool, device=device)
                        # 1. reference fwd and bwd
                        output_ref, _, mask_ref = pack_and_pad(input_global, mask_global, axis, W, keep_input_padding)
                        grad_output = torch.randn_like(output_ref)
                        # backprop without padding in the reference, which additionally affirms that the extra padding
                        # in DTensor version doesn't contribute to its input gradients regardless of the padding values

                        # We can only backward one of them or both summed.
                        # Let's verify sum
                        torch.autograd.backward([output_ref], [grad_output.detach().clone()])
                        grad_input_ref = input_global.grad.clone()
                        input_global.grad.zero_()

                        # 2. Distribute
                        # Shard the leading axis if they exist and evenly divisible by the mesh shape
                        # otherwise replicate
                        if (
                            ndim > 1
                            and axis != 0
                            and device_mesh.ndim >= 2
                            and input_shape[0] % device_mesh.size(0) == 0
                            and shape_mask[0] % device_mesh.size(0) == 0
                        ):
                            placements = (Shard(0),) + (Replicate(),) * (device_mesh.ndim - 2) + (Shard(axis),)
                        else:
                            placements = (Replicate(),) * (device_mesh.ndim - 1) + (Shard(axis),)
                        i_dim_mesh_shard_axis = placements.index(Shard(axis))
                        world_size_shard_axis = device_mesh.shape[i_dim_mesh_shard_axis]

                        # Use detach().clone() to ensure leaf tensor
                        input_dtensor = distribute_tensor(
                            input_global.detach().clone(), device_mesh, placements
                        ).requires_grad_(True)
                        mask_dtensor = distribute_tensor(mask_global.detach().clone(), device_mesh, placements)

                        # 3. Distributed Unmask
                        input_dtensor_copy_local = input_dtensor.detach().clone().to_local().requires_grad_(True)
                        mask_dtensor_copy_local = mask_dtensor.detach().clone().to_local()
                        output_dtensor, output_mask_dtensor = distributed_pack_and_pad(
                            input_dtensor, mask_dtensor, W, axis, keep_input_padding=keep_input_padding
                        )

                        # assert no change to the fwd inputs
                        assert_tensors_identical(
                            input_dtensor_copy_local, input_dtensor.to_local(), check_grad=False, check_grad_fn=False
                        )
                        assert_tensors_identical(
                            mask_dtensor_copy_local, mask_dtensor.to_local(), check_grad=False, check_grad_fn=False
                        )

                        # 4. verify fwd output
                        output_gathered = output_dtensor.full_tensor()
                        mask_gathered = output_mask_dtensor.full_tensor()

                        # verify the output.shape[axis] is padded to (W * world_size)
                        n_pad_target = W * world_size_shard_axis
                        if keep_input_padding:
                            n_expected = ((input_dtensor.shape[axis] + n_pad_target - 1) // n_pad_target) * n_pad_target
                        else:
                            n_valid_max = mask_global.sum(dim=axis).max().item()
                            n_expected = ((n_valid_max + n_pad_target - 1) // n_pad_target) * n_pad_target

                        assert (
                            output_dtensor.shape[axis] == n_expected
                        ), f"{label} output shape mismatch: {output_dtensor.shape[axis]} != {n_expected}"

                        # The distributed output might be padded more than the reference output
                        # because distributed logic aligns to (W * world_size)
                        pad_len = output_gathered.shape[axis] - output_ref.shape[axis]
                        assert pad_len >= 0, "Distributed output should be at least as large as reference"

                        if pad_len > 0:
                            # Pad output_ref along axis
                            pad_arg = [0] * (2 * output_ref.ndim)
                            pad_idx = (output_ref.ndim - 1 - axis) * 2 + 1
                            pad_arg[pad_idx] = pad_len
                            output_ref_padded = torch.nn.functional.pad(output_ref, pad_arg)
                            mask_ref_padded = torch.nn.functional.pad(mask_ref, pad_arg)
                            grad_output_padded = torch.nn.functional.pad(grad_output, pad_arg)
                        else:
                            output_ref_padded = output_ref
                            mask_ref_padded = mask_ref
                            grad_output_padded = grad_output

                        # target function involves no FLOPS so we can use strict equality
                        assert_tensors_identical(
                            output_gathered,
                            output_ref_padded,
                            check_stride=False,
                            check_grad=False,
                            check_grad_fn=False,
                            msg=lambda m: f"{label} output mismatch:\n  {m}",
                        )

                        assert_tensors_identical(
                            mask_gathered,
                            mask_ref_padded,
                            check_stride=False,
                            check_grad=False,
                            check_grad_fn=False,
                            msg=lambda m: f"{label} mask mismatch:\n  {m}",
                        )

                        # 5. Backward

                        # Distribute grad
                        grad_dtensor = distribute_tensor(grad_output_padded.detach().clone(), device_mesh, placements)
                        # We need dummy grads for the masks

                        # Backward with both
                        # We need to manually call backward because output_dtensor.backward() only accepts one gradient?
                        # No, torch.autograd.backward accepts tensors and grad_tensors
                        torch.autograd.backward([output_dtensor], [grad_dtensor])

                        grad_gathered = input_dtensor.grad.full_tensor()

                        assert_tensors_identical(
                            grad_gathered,
                            grad_input_ref,
                            check_grad=False,
                            check_grad_fn=False,
                            msg=lambda m: f"{label} input gradient mismatch:\n  {m}",
                        )

                        # 6. Test Inverse Forward
                        # We detach output_dtensor_qw to test Inverse backward in isolation
                        output_dtensor_detached = output_dtensor.detach().clone().requires_grad_(True)
                        input_recovered_dtensor = distributed_unpad_and_unpack(
                            output_dtensor_detached,
                            output_mask_dtensor,
                            mask_dtensor,
                            axis,
                            keep_input_padding=keep_input_padding,
                        )

                        # Verify Forward (Recovered Input should match original Input where mask is True)
                        # The inverse operation recovers the valid elements into their original positions.
                        # Invalid positions (where mask is False) are zeroed out by scatter.
                        input_recovered_global = input_recovered_dtensor.full_tensor()

                        # Apply mask to original input for comparison
                        input_global_masked = input_global * mask_global.to(input_global.dtype)
                        # Inverse output also has zeros at invalid positions naturally

                        assert_tensors_identical(
                            input_recovered_global,
                            input_global_masked,
                            check_stride=False,
                            check_grad=False,
                            check_grad_fn=False,
                            msg=lambda m: f"{label} inverse forward mismatch:\n  {m}",
                        )

                        # 7. Test Inverse Backward
                        grad_input_recovered = torch.randn_like(input_global)
                        # Distribute the gradient
                        grad_input_recovered_dtensor = distribute_tensor(
                            grad_input_recovered.detach().clone(), device_mesh, placements
                        )

                        input_recovered_dtensor.backward(grad_input_recovered_dtensor)

                        grad_output_result = output_dtensor_detached.grad.full_tensor()

                        # Reference Backward: The backward of Inverse is Unmask (Gather)
                        # We use the pack and pad reference implementation on the gradient
                        grad_output_ref_qw, _, _ = pack_and_pad(
                            grad_input_recovered, mask_global, axis, W, keep_input_padding
                        )

                        # Pad reference to match distributed output shape
                        pad_len = grad_output_result.shape[axis] - grad_output_ref_qw.shape[axis]
                        if pad_len > 0:
                            pad_arg = [0] * (2 * grad_output_ref_qw.ndim)
                            pad_idx = (grad_output_ref_qw.ndim - 1 - axis) * 2 + 1
                            pad_arg[pad_idx] = pad_len
                            grad_output_ref_padded = torch.nn.functional.pad(grad_output_ref_qw, pad_arg)
                        else:
                            grad_output_ref_padded = grad_output_ref_qw

                        assert_tensors_identical(
                            grad_output_result,
                            grad_output_ref_padded,
                            check_stride=False,
                            check_grad=False,
                            check_grad_fn=False,
                            msg=lambda m: f"{label} inverse backward mismatch:\n  {m}",
                        )

    DistributedManager.cleanup()


@pytest.mark.parametrize(
    "setup_env",
    [
        ((1, 4), True, "cuda", "ENV"),
        ((1, 7), True, "cuda", "ENV"),
        ((2, 4), True, "cuda", "ENV"),
    ],
    indirect=("setup_env",),
    ids=lambda x: f"dp:{x[0][0]}, cp:{x[0][1]}, specify_method:{x[1]}, device_type:{x[2]}, method_init:{x[3]}",
)
def test_distributed_pack_and_pad(setup_env):
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    spawn_multiprocessing(
        parallel_assert_pack_and_pad, world_size, grid_group_sizes, device_type, backend, env_per_rank
    )


def parallel_assert_pack_and_pad_and_gather_sliding_windows(rank, grid_group_sizes, device_type, backend, env_map):
    """Test integration of distributed pack and pad -> gather sliding windows."""

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
    device_mesh = manager.device_mesh

    # 1. Create global data
    seed_by_rank(0, 42)

    for W, H, K_per_rank in [(8, 24, 3), (32, 128, 8), (32, 128, 1)]:
        assert H % (W // 2) == 0, f"H must be divisible by W // 2 but got {H} % {W // 2} != 0"
        h = H // (W // 2)
        assert h % 2 == 0, f"h := H // (W // 2) must be divisible by 2 but got {h} % 2 != 0"
        K = K_per_rank * grid_group_sizes["cp"]
        h = H // (W // 2)
        for n_total in [K * W, K * W + 10 * grid_group_sizes["cp"]]:
            # always shard the (2 * K, ) axis by the "cp" submesh
            for input_shape_extra in [(None,), (4, None), (None, 3), (2, None, 3)]:
                # The "None" element indicates the axis to be used as "half-window" axis
                # there can be only one "None" element in the input_shape
                if sum(1 for x in input_shape_extra if x is None) != 1:
                    raise ValueError(
                        f"There can be one and only one 'None' element in the input_shape but got {input_shape_extra}"
                    )
                axis = input_shape_extra.index(None)
                input_shape = input_shape_extra[:axis] + (n_total,) + input_shape_extra[axis + 1 :]
                shape_mask = input_shape
                ndim = len(input_shape)
                # Shard the leading axis if they exist and evenly divisible by the mesh shape
                # otherwise replicate
                if (
                    ndim > 1
                    and axis != 0
                    and device_mesh.ndim >= 2
                    and input_shape[0] % device_mesh.size(0) == 0
                    and shape_mask[0] % device_mesh.size(0) == 0
                ):
                    # don't shard the dim after "axis" so that these placements can be reused for
                    # sharding the upstream adjoint for testing the backward because the fwd output
                    # inserts new axes after the "axis" dimension
                    placements = (Shard(0),) + (Replicate(),) * (device_mesh.ndim - 2) + (Shard(axis),)
                else:
                    placements = (Replicate(),) * (device_mesh.ndim - 1) + (Shard(axis),)

                label = f"W:{W}, H:{H}, K:{K}, input_shape:{input_shape}"

                # 1. Setup Global Data
                input_global = torch.randn(input_shape, requires_grad=True, device=device)
                mask_global = torch.randint(0, 2, shape_mask, dtype=torch.bool, device=device)

                # 2. Distribute
                input_dtensor = distribute_tensor(
                    input_global.detach().clone().requires_grad_(True), device_mesh, placements
                )
                mask_dtensor = distribute_tensor(mask_global.detach().clone(), device_mesh, placements)

                # 3. Distributed Chain
                # Unmask
                packed_dtensor, packed_mask = distributed_pack_and_pad(input_dtensor, mask_dtensor, W, axis)

                # Due to keep_input_padding=False, input_dtensor.shape[axis] is trimmed to the maximum number of valid elements
                # which is then padded to the multiple of W * size_group_sharding_axis so we need to recompute
                # K in order to reshape the tensors
                K_global = packed_dtensor.shape[axis] // W

                packed_dtensor_hw = shardwise_unflatten_sharded(packed_dtensor, axis, (2 * K_global, W // 2))

                # Gather
                output_dtensor = distributed_gather_sliding_windows(packed_dtensor_hw, h, axis)

                # 4. Single Device Chain
                packed_ref, _, mask_ref = pack_and_pad(input_global, mask_global, axis, W)
                # Apply masks
                packed_ref = packed_ref * mask_ref.to(packed_ref.dtype)

                # packed_ref_hw: (2K, W/2, F)

                # Need to match the global padding logic of distributed_pack_and_pad!
                # Distributed pads total_valid to multiple of W*world_size.
                # Single device pads to multiple of W.
                # We need to pad packed_ref to match packed_dtensor shape for comparison.

                if packed_dtensor.shape[axis] > packed_ref.shape[axis]:
                    # Pad along axis
                    pad_arg = [0] * (2 * packed_ref.ndim)
                    idx = (packed_ref.ndim - 1 - axis) * 2 + 1
                    pad_arg[idx] = packed_dtensor.shape[axis] - packed_ref.shape[axis]
                    packed_ref_padded = torch.nn.functional.pad(packed_ref, pad_arg)

                else:
                    packed_ref_padded = packed_ref

                packed_ref_hw_padded = packed_ref_padded.unflatten(axis, (2 * K_global, W // 2))

                # Gather Ref
                # We need offsets
                offset_start = 1 - h // 2
                offsets_ref = torch.arange(offset_start, offset_start + 2 * K_global, 2, device=device)
                output_ref = gather_sliding_windows(packed_ref_hw_padded, offsets_ref, h, axis)

                # 5. Verify Forward
                output_gathered = output_dtensor.full_tensor()
                torch.testing.assert_close(
                    output_gathered, output_ref, atol=0, rtol=0, msg=lambda m: f"{label} output mismatch:\n  {m}"
                )

                # 6. Verify Backward
                grad_out = torch.randn_like(output_gathered)
                grad_dtensor = distribute_tensor(grad_out.detach().clone(), device_mesh, placements)

                # dummy gradients for masks (required by autograd since we return them)
                # Note: we don't pass them to backward() because masks are non-differentiable
                # and PyTorch will pass None to the backward method for them.

                torch.autograd.backward([output_dtensor], [grad_dtensor])

                # Single backward
                # We use output_ref and packed_ref_qw as roots
                # We directly backprop with the padding, assuring that the extra padding in DTensor case
                # behave similarly
                torch.autograd.backward([output_ref], [grad_out])

                # Verify input gradients
                grad_input_gathered = input_dtensor.grad.full_tensor()
                torch.testing.assert_close(
                    grad_input_gathered,
                    input_global.grad,
                    msg=lambda m: f"{label} input gradient mismatch:\n  {m}",
                )

    DistributedManager.cleanup()


@pytest.mark.parametrize(
    "setup_env",
    [
        ((1, 4), True, "cuda", "ENV"),
        ((1, 7), True, "cuda", "ENV"),
        ((2, 4), True, "cuda", "ENV"),
    ],
    indirect=("setup_env",),
    ids=lambda x: f"dp:{x[0][0]}, cp:{x[0][1]}, specify_method:{x[1]}, device_type:{x[2]}, method_init:{x[3]}",
)
def test_distributed_pack_and_pad_and_gather_sliding_windows(setup_env):
    """
    Test integration of distributed pack and pad -> gather sliding windows.
    """
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    spawn_multiprocessing(
        parallel_assert_pack_and_pad_and_gather_sliding_windows,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
    )


def parallel_assert_window_batch_attention(rank, grid_group_sizes, device_type, backend, env_map, dtype):
    """Test integration of distributed pack and pad -> gather sliding windows -> attention."""

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
    device_mesh = manager.device_mesh

    # 1. Create global data
    seed_by_rank(0, 42)

    for W, H, K_per_rank in [(8, 16, 1), (32, 128, 8), (32, 128, 1)]:
        assert H % (W // 2) == 0, f"H must be divisible by W // 2 but got {H} % {W // 2} != 0"
        h = H // (W // 2)
        assert h % 2 == 0, f"h := H // (W // 2) must be divisible by 2 but got {h} % 2 != 0"
        K = K_per_rank * grid_group_sizes["cp"]
        h = H // (W // 2)
        for n_total in [K * W]:  # Boltz window batching only supports K * W size input
            for input_shape_extra in [
                (2, None, 3),
                (16, None, 128),
            ]:  # single_to_keys requires (B, n, D) shape
                # The "None" element indicates the axis to be used as "half-window" axis
                # there can be only one "None" element in the input_shape
                if sum(1 for x in input_shape_extra if x is None) != 1:
                    raise ValueError(
                        f"There can be one and only one 'None' element in the input_shape but got {input_shape_extra}"
                    )
                axis = input_shape_extra.index(None)
                input_shape = input_shape_extra[:axis] + (n_total,) + input_shape_extra[axis + 1 :]
                shape_mask = input_shape
                ndim = len(input_shape)
                # Shard the leading axis if they exist and evenly divisible by the mesh shape
                # otherwise replicate
                if (
                    ndim > 1
                    and axis != 0
                    and device_mesh.ndim >= 2
                    and input_shape[0] % device_mesh.size(0) == 0
                    and shape_mask[0] % device_mesh.size(0) == 0
                ):
                    # don't shard the dim after "axis" so that these placements can be reused for
                    # sharding the upstream adjoint for testing the backward because the fwd output
                    # inserts new axes after the "axis" dimension
                    placements = (Shard(0),) + (Replicate(),) * (device_mesh.ndim - 2) + (Shard(axis),)
                else:
                    placements = (Replicate(),) * (device_mesh.ndim - 1) + (Shard(axis),)

                label = f"W:{W}, H:{H}, K:{K}, input_shape:{input_shape}"

                # 1. Setup Global Data
                input_global = torch.randn(input_shape, dtype=dtype, requires_grad=True, device=device)
                mask_global = torch.randint(0, 2, shape_mask, dtype=torch.bool, device=device, requires_grad=False)

                # sort the input and mask so that valid elements are leading the axis
                mask_global_sorted, argsort_mask_global = torch.sort(
                    mask_global, dim=axis, descending=True, stable=True
                )
                input_global_sorted = torch.gather(input_global, axis, argsort_mask_global)

                indexing_matrix = get_indexing_matrix(K, W, H, device).to(dtype=input_global.dtype)

                # Get key and query from Boltz window batching (..., K, H, ...)
                key_wb_expected = single_to_keys(input_global_sorted, indexing_matrix, W, H)
                mask_key_wb_expected = single_to_keys(
                    mask_global_sorted.to(dtype=indexing_matrix.dtype), indexing_matrix, W, H
                ).to(dtype=torch.bool)

                query_wb_expected = input_global_sorted.unflatten(axis, (K, W))
                mask_query_wb_expected = mask_global_sorted.unflatten(axis, (K, W))

                # perform linear attention with single iteration (to avoid numerical complications from softmax etc but enough to
                # prove the DTensor logics are correct)
                a_wb_expected = torch.einsum(
                    "bkid,bkjd->bkij",
                    query_wb_expected * mask_query_wb_expected,
                    key_wb_expected * mask_key_wb_expected,
                )
                o_wb_expected = torch.einsum("bkml,bkld->bkmd", a_wb_expected, key_wb_expected * mask_key_wb_expected)

                # 2. Distribute
                input_dtensor = distribute_tensor(
                    input_global.detach().clone().requires_grad_(True), device_mesh, placements
                )
                mask_dtensor = distribute_tensor(mask_global.detach().clone(), device_mesh, placements)

                # 3. Distributed Chain
                # Unmask
                packed_dtensor, packed_mask = distributed_pack_and_pad(input_dtensor, mask_dtensor, W, axis)

                # Due to keep_input_padding=False, input_dtensor.shape[axis] is trimmed to the maximum number of valid elements
                # which is then padded to the multiple of W * size_group_sharding_axis so we need to recompute
                # K in order to reshape the tensors
                K_global = packed_dtensor.shape[axis] // W

                packed_dtensor_hw = shardwise_unflatten_sharded(packed_dtensor, axis, (2 * K_global, W // 2))
                packed_mask_hw = shardwise_unflatten_sharded(packed_mask, axis, (2 * K_global, W // 2))

                # Gather
                # (..., K, h, W//2) -> (..., K, H, ...)
                key_wb_result_reshaped = distributed_gather_sliding_windows(packed_dtensor_hw, h, axis)
                key_wb_result = key_wb_result_reshaped.flatten(axis + 1, axis + 2)
                mask_key_wb_result_reshaped = distributed_gather_sliding_windows(packed_mask_hw, h, axis)
                mask_key_wb_result = mask_key_wb_result_reshaped.flatten(axis + 1, axis + 2)

                query_wb_result = shardwise_unflatten_sharded(packed_dtensor, axis, (K_global, W))
                mask_query_wb_result = shardwise_unflatten_sharded(packed_mask, axis, (K_global, W))

                # DTensor linear attention
                # The DTensor native einsum doesn't support the involved implicit unflattening of two
                # leading batch dimensions so we need to do via shard-wise local attention
                a_wb_result_local = torch.einsum(
                    "bkid,bkjd->bkij",
                    query_wb_result.to_local() * mask_query_wb_result.to_local(),
                    key_wb_result.to_local() * mask_key_wb_result.to_local(),
                )
                o_wb_result_local = torch.einsum(
                    "bkml,bkld->bkmd",
                    a_wb_result_local,
                    key_wb_result.to_local() * mask_key_wb_result.to_local(),
                )
                a_wb_result = DTensor.from_local(a_wb_result_local, device_mesh, placements)
                o_wb_result = DTensor.from_local(o_wb_result_local, device_mesh, placements)

                # verify key and query consistency
                # this requires potentially padding the reference to the next multiple of  W * world_size
                key_wb_result_global = key_wb_result.full_tensor()
                mask_key_wb_result_global = mask_key_wb_result.full_tensor()
                query_wb_result_global = query_wb_result.full_tensor()
                mask_query_wb_result_global = mask_query_wb_result.full_tensor()
                a_wb_result_global = a_wb_result.full_tensor()
                o_wb_result_global = o_wb_result.full_tensor()
                pad_result = False
                pad_expected = False
                if key_wb_result_global.shape[axis] > key_wb_expected.shape[axis]:
                    pad_result = True
                    # pad expected
                    pad_arg = [0] * (2 * key_wb_expected.ndim)
                    idx = (key_wb_expected.ndim - 1 - axis) * 2 + 1
                    pad_arg[idx] = key_wb_result_global.shape[axis] - key_wb_expected.shape[axis]
                    key_wb_expected_padded = torch.nn.functional.pad(key_wb_expected, pad_arg)
                    mask_key_wb_expected_padded = torch.nn.functional.pad(mask_key_wb_expected, pad_arg)
                    query_wb_expected_padded = torch.nn.functional.pad(query_wb_expected, pad_arg)
                    mask_query_wb_expected_padded = torch.nn.functional.pad(mask_query_wb_expected, pad_arg)
                    a_wb_expected_padded = torch.nn.functional.pad(a_wb_expected, pad_arg)
                    o_wb_expected_padded = torch.nn.functional.pad(o_wb_expected, pad_arg)
                    # use result
                    key_wb_result_global_padded = key_wb_result_global
                    mask_key_wb_result_global_padded = mask_key_wb_result_global
                    query_wb_result_global_padded = query_wb_result_global
                    mask_query_wb_result_global_padded = mask_query_wb_result_global
                    a_wb_result_global_padded = a_wb_result_global
                    o_wb_result_global_padded = o_wb_result_global
                elif key_wb_result_global.shape[axis] < key_wb_expected.shape[axis]:
                    pad_result = True
                    # pad result
                    pad_arg = [0] * (2 * key_wb_result_global.ndim)
                    idx = (key_wb_result_global.ndim - 1 - axis) * 2 + 1
                    pad_arg[idx] = key_wb_expected.shape[axis] - key_wb_result_global.shape[axis]
                    key_wb_result_global_padded = torch.nn.functional.pad(key_wb_result_global, pad_arg)
                    mask_key_wb_result_global_padded = torch.nn.functional.pad(mask_key_wb_result_global, pad_arg)
                    query_wb_result_global_padded = torch.nn.functional.pad(query_wb_result_global, pad_arg)
                    mask_query_wb_result_global_padded = torch.nn.functional.pad(mask_query_wb_result_global, pad_arg)
                    a_wb_result_global_padded = torch.nn.functional.pad(a_wb_result_global, pad_arg)
                    o_wb_result_global_padded = torch.nn.functional.pad(o_wb_result_global, pad_arg)
                    # use expected
                    key_wb_expected_padded = key_wb_expected
                    mask_key_wb_expected_padded = mask_key_wb_expected
                    query_wb_expected_padded = query_wb_expected
                    mask_query_wb_expected_padded = mask_query_wb_expected
                    a_wb_expected_padded = a_wb_expected
                    o_wb_expected_padded = o_wb_expected
                else:
                    pad_arg = None  # prevent accidental reuse from previous for loop iteration
                    key_wb_expected_padded = key_wb_expected
                    mask_key_wb_expected_padded = mask_key_wb_expected
                    query_wb_expected_padded = query_wb_expected
                    mask_query_wb_expected_padded = mask_query_wb_expected
                    a_wb_expected_padded = a_wb_expected
                    o_wb_expected_padded = o_wb_expected
                    key_wb_result_global_padded = key_wb_result_global
                    mask_key_wb_result_global_padded = mask_key_wb_result_global
                    query_wb_result_global_padded = query_wb_result_global
                    mask_query_wb_result_global_padded = mask_query_wb_result_global
                    a_wb_result_global_padded = a_wb_result_global
                    o_wb_result_global_padded = o_wb_result_global

                torch.testing.assert_close(
                    mask_key_wb_result_global_padded,
                    mask_key_wb_expected_padded,
                    msg=lambda m: f"{label} mask_key_wb mismatch:\n  {m}",
                )
                torch.testing.assert_close(
                    key_wb_result_global_padded * mask_key_wb_result_global_padded,
                    key_wb_expected_padded * mask_key_wb_expected_padded,
                    msg=lambda m: f"{label} key_wb mismatch:\n  {m}",
                )

                torch.testing.assert_close(
                    mask_query_wb_result_global_padded,
                    mask_query_wb_expected_padded,
                    msg=lambda m: f"{label} mask_query_wb mismatch:\n  {m}",
                )
                torch.testing.assert_close(
                    query_wb_result_global_padded * mask_query_wb_result_global_padded,
                    query_wb_expected_padded * mask_query_wb_expected_padded,
                    msg=lambda m: f"{label} query_wb mismatch:\n  {m}",
                )
                torch.testing.assert_close(
                    a_wb_result_global_padded,
                    a_wb_expected_padded,
                    msg=lambda m: f"{label} a_wb mismatch:\n  {m}",
                )
                torch.testing.assert_close(
                    o_wb_result_global_padded,
                    o_wb_expected_padded,
                    msg=lambda m: f"{label} o_wb mismatch:\n  {m}",
                )

                # check backward pass
                # To make sharding the upstream adjoint easy, we always use the DTensor full_tensor version
                # to generate the upstream adjoints and make necessary padding to the expected version
                # We also zeros out the upstream adjoints to check if the invalid elements in input.grad are zeros
                grad_o_wb_result_global = torch.randn_like(o_wb_result_global) * mask_query_wb_result_global
                if pad_result:
                    # grad_o_wb_result_global is shorter along 'axis' than o_wb_expected_padded
                    grad_o_wb_expected_padded = torch.nn.functional.pad(
                        grad_o_wb_result_global.detach().clone(), pad_arg
                    )
                elif pad_expected:
                    # grad_o_wb_result_global is longer along 'axis'
                    grad_o_wb_expected_padded = (
                        grad_o_wb_result_global.detach().clone().narrow(axis, 0, o_wb_expected_padded.shape[axis])
                    )
                else:
                    grad_o_wb_expected_padded = grad_o_wb_result_global.detach().clone()
                o_wb_expected_padded.backward(grad_o_wb_expected_padded)

                grad_o_wb_dtensor = distribute_tensor(grad_o_wb_result_global.detach().clone(), device_mesh, placements)
                o_wb_result.backward(grad_o_wb_dtensor)

                # verify input gradients
                grad_input_result_global = input_dtensor.grad.full_tensor()
                mask_result_global = mask_dtensor.full_tensor()
                torch.testing.assert_close(
                    grad_input_result_global * mask_result_global,
                    input_global.grad * mask_global,
                    msg=lambda m: f"{label} input gradient mismatch:\n  {m}",
                )

                assert_tensors_identical(
                    grad_input_result_global * (~mask_result_global),
                    torch.zeros_like(grad_input_result_global),
                    check_grad=False,
                    check_grad_fn=False,
                    msg=lambda m: f"{label} input gradient mismatch for invalid elements:\n  {m}",
                )

    DistributedManager.cleanup()


def parallel_assert_single_to_key(rank, grid_group_sizes, device_type, backend, env_map, dtype):
    """Test convert_single_repr_to_window_batched_key against single_to_keys (fwd+bwd)."""

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
    device_mesh = manager.device_mesh

    seed_by_rank(0, 42)

    # Keep the same window-batching tuples used elsewhere: they satisfy Boltz constraints.
    for W, H, K_per_rank in [(8, 16, 1), (32, 128, 8), (32, 128, 1)]:
        assert H % (W // 2) == 0, f"H must be divisible by W // 2 but got {H} % {W // 2} != 0"
        h = H // (W // 2)
        assert h % 2 == 0, f"h := H // (W // 2) must be divisible by 2 but got {h} % 2 != 0"

        K = K_per_rank * grid_group_sizes["cp"]
        N = K * W

        indexing_matrix = get_indexing_matrix(K, W, H, device).to(dtype=dtype)

        # convert_single_repr_to_window_batched_key assumes the sequence axis is dim=1.
        for input_shape_extra in [
            (2, None, 3),
            (16, None, 128),
        ]:
            axis = input_shape_extra.index(None)
            assert axis == 1, "This test assumes the sequence axis is dim=1."
            input_shape = input_shape_extra[:axis] + (N,) + input_shape_extra[axis + 1 :]
            ndim = len(input_shape)
            label = f"W:{W}, H:{H}, K:{K}, input_shape:{input_shape}"

            input_global = torch.randn(input_shape, dtype=dtype, requires_grad=True, device=device)

            # Reference: Boltz window-batched keys.
            key_expected = single_to_keys(input_global, indexing_matrix, W, H)

            # Distribute input on device mesh; shard the sequence axis (dim=1) on the last mesh dim.
            if ndim > 2 and device_mesh.ndim >= 2 and input_shape[0] % device_mesh.size(0) == 0:
                placements = (Shard(0),) + (Replicate(),) * (device_mesh.ndim - 2) + (Shard(1),)
            else:
                placements = (Replicate(),) * (device_mesh.ndim - 1) + (Shard(1),)

            input_dtensor = distribute_tensor(
                input_global.detach().clone().requires_grad_(True),
                device_mesh,
                placements,
            )

            # Distributed: window-batched keys.
            key_dtensor = convert_single_repr_to_window_batched_key(input_dtensor, W, H)
            key_gathered = key_dtensor.full_tensor()

            torch.testing.assert_close(
                key_gathered,
                key_expected,
                msg=lambda m: f"{label} fwd key mismatch:\n  {m}",
            )

            # Backward: compare input gradients.
            grad_key = torch.randn_like(key_expected)
            key_expected.backward(grad_key)

            grad_key_dtensor = distribute_tensor(grad_key.detach().clone(), device_mesh, key_dtensor.placements)
            key_dtensor.backward(grad_key_dtensor)

            grad_input_gathered = input_dtensor.grad.full_tensor()
            torch.testing.assert_close(
                grad_input_gathered,
                input_global.grad,
                msg=lambda m: f"{label} input gradient mismatch:\n  {m}",
            )

    DistributedManager.cleanup()


@pytest.mark.parametrize(
    "setup_env",
    [
        ((1, 4), True, "cuda", "ENV"),
        ((1, 7), True, "cuda", "ENV"),
        ((2, 4), True, "cuda", "ENV"),
    ],
    indirect=("setup_env",),
    ids=lambda x: f"dp:{x[0][0]}, cp:{x[0][1]}, specify_method:{x[1]}, device_type:{x[2]}, method_init:{x[3]}",
)
def test_distributed_window_batch_attention(setup_env):
    """
    Test integration of distributed pack and pad -> gather sliding windows -> attention.
    This test uses Boltz get_indexing_matrix and single_to_keys() to create the window batched query and key
    to perform the reference attention. It then uses the distributed pack and pad -> gather sliding windows
    to perform the distributed attention. Both forward and backward passes are tested.
    """
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    spawn_multiprocessing(
        parallel_assert_window_batch_attention,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        torch.float64,  # use FP64 to avoid dealing with numerical tolerance
    )


@pytest.mark.parametrize(
    "setup_env",
    [
        ((1, 4), True, "cuda", "ENV"),
        ((1, 7), True, "cuda", "ENV"),
        ((2, 4), True, "cuda", "ENV"),
    ],
    indirect=("setup_env",),
    ids=lambda x: f"dp:{x[0][0]}, cp:{x[0][1]}, specify_method:{x[1]}, device_type:{x[2]}, method_init:{x[3]}",
)
def test_distributed_single_to_key(setup_env):
    """Test distributed convert_single_repr_to_window_batched_key vs single_to_keys."""

    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    spawn_multiprocessing(
        parallel_assert_single_to_key,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        torch.float64,  # use FP64 to avoid dealing with numerical tolerance
    )
