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


from typing import Dict, List, Tuple

import torch
import torch.distributed as dist
from torch.distributed._tensor import DTensor, Shard
from torch.distributed.tensor import Partial

from boltz.distributed.model.layers.outer_gather import get_overlap_from_peers
from boltz.distributed.utils import update_exhaustive_strides


class DistributedGather(torch.autograd.Function):
    @staticmethod
    def forward(
        ctx,
        x_dtensor: DTensor,
        idx_dtensor: DTensor,
        axis: int,
        are_ids_contiguous: bool,
        idx_mask: DTensor | None = None,
    ) -> DTensor:
        """Distributed 1D gather.

        Args:
            x_dtensor: DTensor with shape ``(*batch, N, *features)``.
            idx_dtensor: DTensor with shape ``(*batch, K, W)`` that provides gather
                indices into the ``N`` dimension of ``x_dtensor``.
            axis: Axis in ``x_dtensor`` corresponding to ``N`` (and to ``K`` in
                ``idx_dtensor``).
            are_ids_contiguous: This is a heuristic for selecting the underlying
                send/recv strategy for performance purpose. Currently only True is supported,
                which means that the idx_dtensor maps to a contiguous block of x along (axis)
                dimensions for all the shards and for all the leading (batch) dimensions.
                When True, the underlying strategy will use the min/max of idx_dtensor to
                compute the needed interval, assuming the resulting buffer to be communicated
                across the ranks is fully (or approximately so) utilized. In the case the
                data inside idx_dtensor doesn't mapped to contiguous blocks, the result will
                still be correct by setting are_ids_contiguous=True but the buffer of x chunks
                communicated will contain a lot of unused elements, making the sed/recv inefficient.
            idx_mask: Optional DTensor with shape ``(*batch, K, W)`` and same device_mesh
                and placements as ``idx_dtensor``. Elements with True indicate valid indices,
                elements with False indicate invalid indices that should be ignored.
        """
        if not are_ids_contiguous:
            raise NotImplementedError("DistributedGather currently only supports are_ids_contiguous=True")

        if not isinstance(x_dtensor, DTensor) or not isinstance(idx_dtensor, DTensor):
            raise TypeError("x_dtensor and idx_dtensor must be DTensors")

        batch_dims_x = x_dtensor.shape[:axis]
        if batch_dims_x != idx_dtensor.shape[:axis]:
            raise ValueError(f"Batch dimensions must match: x {batch_dims_x} vs idx {idx_dtensor.shape[:axis]}")

        mesh = x_dtensor.device_mesh
        if idx_dtensor.device_mesh != mesh:
            raise ValueError("x and idx must be on the same DeviceMesh")
        if idx_dtensor.placements != x_dtensor.placements:
            raise ValueError("x and idx must have identical placements")

        # Validate idx_mask if provided
        if idx_mask is not None:
            if not isinstance(idx_mask, DTensor):
                raise TypeError("idx_mask must be a DTensor")
            if idx_mask.shape != idx_dtensor.shape:
                raise ValueError(f"idx_mask shape {idx_mask.shape} must match idx_dtensor shape {idx_dtensor.shape}")
            if idx_mask.device_mesh != idx_dtensor.device_mesh:
                raise ValueError("idx_mask must have the same device_mesh as idx_dtensor")
            if idx_mask.placements != idx_dtensor.placements:
                raise ValueError("idx_mask must have the same placements as idx_dtensor")
            if idx_mask.dtype != torch.bool:
                raise TypeError(
                    f"idx_mask must have dtype torch.bool, got {idx_mask.dtype}. Use mask.bool() to convert."
                )

        x_placements = x_dtensor.placements
        idx_placements = idx_dtensor.placements

        ndim_x = x_dtensor.ndim
        if axis < 0:
            axis += ndim_x
        if axis < 0 or axis >= ndim_x:
            raise ValueError(f"axis {axis} out of range for x.ndim={ndim_x}")

        # Identify shard axis on mesh for the gather dimension
        mesh_dim_axis = None
        for i_mesh_dim, p in enumerate(x_placements):
            if isinstance(p, Partial) or isinstance(idx_placements[i_mesh_dim], Partial):
                raise ValueError("Partial placements are not supported")
            if isinstance(p, Shard):
                if p.dim == axis:
                    mesh_dim_axis = i_mesh_dim
                # Enforce even sharding on axis for both x and idx since we require identical device_mesh
                # and placements between the two
                if x_dtensor.shape[p.dim] % mesh.size(i_mesh_dim) != 0:
                    raise ValueError(
                        f"x_dtensor axis {p.dim} size {x_dtensor.shape[p.dim]} not evenly divisible by mesh dim {i_mesh_dim}"
                    )
                if idx_dtensor.shape[p.dim] % mesh.size(i_mesh_dim) != 0:
                    raise ValueError(
                        f"idx_dtensor axis {p.dim} size {idx_dtensor.shape[p.dim]} not evenly divisible by mesh dim {i_mesh_dim}"
                    )

        if mesh_dim_axis is None:
            raise ValueError(f"x must be sharded along axis {axis}")

        x_local = x_dtensor.to_local()
        idx_local = idx_dtensor.to_local()
        idx_mask_local = idx_mask.to_local() if idx_mask is not None else None
        device = x_local.device
        cpu_device = torch.device("cpu")

        # Determine needed interval along axis from local idx
        # need interval is required to have a singleton axis of 1 (representing 1-d gathering)
        # to be used in get_overlap_from_peers
        if idx_local.numel() > 0:
            if idx_mask_local is not None:
                # Only consider valid indices for interval computation
                if idx_mask_local.any():
                    valid_idx = idx_local[idx_mask_local]
                    need_interval = torch.stack(valid_idx.aminmax()).to(dtype=torch.long).unsqueeze(0)  # (1,2)
                    need_interval[:, -1] += 1
                else:
                    # All indices are masked out
                    need_interval = torch.tensor([[0, 0]], device=device, dtype=torch.long)  # (1,2)
            else:
                # aminmax return end-inclusive interval of shape (2,)
                need_interval = torch.stack(idx_local.aminmax()).to(dtype=torch.long).unsqueeze(0)  # (1,2)
                need_interval[:, -1] += 1
        else:
            need_interval = torch.tensor([[0, 0]], device=device, dtype=torch.long)  # (1,2)
        need_start = need_interval[0, 0]
        need_end = need_interval[0, 1]
        need_start_cpu = need_start.to(cpu_device)

        # Owned chunk interval
        coord_axis = mesh.get_local_rank(mesh_dim_axis)
        chunk_size = x_dtensor.shape[axis] // mesh.size(mesh_dim_axis)
        own_start = torch.tensor(coord_axis * chunk_size, device=cpu_device, dtype=torch.long)
        own_end = own_start + chunk_size
        own_interval = torch.stack([own_start, own_end]).unsqueeze(0)  # (1,2)

        # All-gather need intervals along sharded mesh dim (metadata only)
        group_axis = mesh.get_group(mesh_dim_axis)
        need_range = [torch.zeros_like(need_interval) for _ in range(mesh.size(mesh_dim_axis))]
        dist.all_gather(need_range, need_interval, group=group_axis)
        need_range = torch.stack(need_range)  # (size_group,1,2)
        need_range_cpu = need_range.cpu()

        ranks_global_on_mesh = mesh.mesh
        my_coords = mesh.get_coordinate()
        index_list_submesh = []
        for dim in range(ranks_global_on_mesh.ndim):
            if dim == mesh_dim_axis:
                index_list_submesh.append(slice(None))
            else:
                index_list_submesh.append(torch.tensor(my_coords[dim], device=cpu_device))
        ranks_global_on_submesh = ranks_global_on_mesh[tuple(index_list_submesh)]  # (size_group,)

        size_group = mesh.size(mesh_dim_axis)

        # RECEIVE PLAN
        if need_start >= need_end:
            needed_chunks: List[Dict[str, torch.Tensor | int]] = []
        else:
            # peers own intervals
            start_peers_own = torch.arange(size_group, device=cpu_device, dtype=torch.long) * chunk_size
            end_peers_own = start_peers_own + chunk_size
            interval_peers_own = torch.stack([start_peers_own, end_peers_own], dim=-1).unsqueeze(1)  # (size_group,1,2)

            needed_chunks = get_overlap_from_peers(
                ranks_global_on_submesh, interval_peers_own, need_interval.to(cpu_device)
            )

        ops = []
        recv_bufs = {}
        recv_metadata_for_bwd = []

        for item in needed_chunks:
            peer = item["peer"]
            interval = item["interval"]  # (1,2)
            start_global = interval[0, 0]
            length = interval[0, 1] - interval[0, 0]

            shape = list(x_local.shape)
            shape[axis] = length.item()
            buf = torch.empty(shape, dtype=x_local.dtype, device=device)

            if peer == dist.get_rank():
                start_local = start_global - own_start
                buf.copy_(x_local.narrow(axis, start_local.item(), length.item()))
                recv_bufs[peer] = buf
                recv_metadata_for_bwd.append((peer, interval, shape))
            else:
                ops.append(dist.P2POp(dist.irecv, buf, peer))
                recv_bufs[peer] = buf
                recv_metadata_for_bwd.append((peer, interval, shape))

        # SEND PLAN
        send_metadata_for_bwd = []
        send_chunks = get_overlap_from_peers(
            ranks_global_on_submesh,
            need_range_cpu.view(size_group, 1, 2),
            own_interval.view(1, 1, 2),
        )

        my_rank = dist.get_rank()
        for item in send_chunks:
            peer = item["peer"]
            if peer == my_rank:
                continue
            interval = item["interval"]
            start_global = interval[0, 0]
            length = interval[0, 1] - interval[0, 0]
            start_local = start_global - own_start
            chunk = x_local.narrow(axis, start_local.item(), length.item()).contiguous()
            ops.append(dist.P2POp(dist.isend, chunk, peer))
            send_metadata_for_bwd.append(
                (
                    peer,
                    start_local.to(cpu_device, dtype=torch.long),
                    length.to(cpu_device, dtype=torch.long),
                )
            )

        if ops:
            reqs = dist.batch_isend_irecv(ops)
            for req in reqs:
                req.wait()

        # Assemble buffer
        if need_start >= need_end:
            buffer_shape = list(x_local.shape)
            buffer_shape[axis] = 0
            x_buffer = torch.empty(buffer_shape, dtype=x_local.dtype, device=device)
        else:
            buffer_shape = list(x_local.shape)
            buffer_shape[axis] = (need_end - need_start).item()
            x_buffer = torch.zeros(buffer_shape, dtype=x_local.dtype, device=device)

            for item in needed_chunks:
                peer = item["peer"]
                interval = item["interval"]
                buf = recv_bufs.get(peer)
                if buf is None:
                    raise RuntimeError(f"Missing recv buffer for peer {peer}")

                start_global = interval[0, 0]
                length = interval[0, 1] - interval[0, 0]
                start_local_buf = start_global - need_start_cpu
                target = x_buffer.narrow(axis, start_local_buf.item(), length.item())
                target.copy_(buf)

        # Adjust idx into buffer coords
        idx_local_buffer = idx_local - need_start.item()

        # Local computation using linearized gather over (axis, feature) block
        shape_trailing = x_buffer.shape[axis + 1 :]
        shape_trailing_flat = torch.Size(shape_trailing).numel()

        shape_leading = x_buffer.shape[:axis]
        shape_leading_flat = torch.Size(shape_leading).numel()
        L = x_buffer.shape[axis]
        K = idx_local_buffer.shape[-2]
        W = idx_local_buffer.shape[-1]

        # Handle edge case where buffer is empty (all indices masked out)
        if L == 0:
            out_shape = list(shape_leading) + [K, W] + list(shape_trailing)
            out_local = torch.zeros(out_shape, dtype=x_local.dtype, device=device)
        else:
            # For masked indices, clamp to valid buffer range to avoid index errors
            if idx_mask_local is not None:
                idx_local_buffer = torch.where(idx_mask_local, idx_local_buffer, torch.zeros_like(idx_local_buffer))

            x_flat = x_buffer.reshape(shape_leading_flat, L, shape_trailing_flat)
            idx_flat = idx_local_buffer.reshape(shape_leading_flat, K, W)

            # Advanced indexing:
            # dim 0: arange(B) reshaped to (B, 1, 1) to broadcast against (B, K, W)
            batch_idx = torch.arange(shape_leading_flat, device=device).reshape(shape_leading_flat, 1, 1)

            # This performs the gather
            # x_flat[ (B,1,1), (B,K,W), : ] -> (B, K, W, F)
            out_flat = x_flat[batch_idx, idx_flat, :]

            # Zero out invalid positions
            if idx_mask_local is not None:
                mask_flat = idx_mask_local.reshape(shape_leading_flat, K, W, 1).to(out_flat.dtype)
                out_flat = out_flat * mask_flat

            if shape_trailing:
                out_local = out_flat.reshape(*shape_leading, K, W, *shape_trailing)
            else:
                out_local = out_flat.reshape(*shape_leading, K, W)

        out_global_shape = list(idx_dtensor.shape) + list(shape_trailing)
        final_global_shape = tuple(out_global_shape)

        strides_out = update_exhaustive_strides(out_local.shape, out_local.stride(), final_global_shape)

        out_dtensor = DTensor.from_local(
            out_local, idx_dtensor.device_mesh, idx_dtensor.placements, shape=final_global_shape, stride=strides_out
        )

        if idx_mask_local is not None:
            ctx.save_for_backward(idx_local, idx_mask_local)
        else:
            ctx.save_for_backward(idx_local)
        ctx.has_mask = idx_mask_local is not None
        ctx.comm_meta = {
            "recv_metadata_for_bwd": recv_metadata_for_bwd,
            "send_metadata_for_bwd": send_metadata_for_bwd,
            "x_local_shape": x_local.shape,
            "x_buffer_shape": x_buffer.shape,
            "need_interval": need_interval,
            "axis": axis,
            "x_placements": x_placements,
            "x_global_shape": x_dtensor.shape,
            "output_placements": out_dtensor.placements,
            "own_interval": own_interval,
            "device_mesh_output": out_dtensor.device_mesh,
        }

        return out_dtensor

    @staticmethod
    def backward(ctx, grad_output: DTensor) -> Tuple[DTensor, None, None, None, None]:
        if ctx.has_mask:
            idx_local, idx_mask_local = ctx.saved_tensors
        else:
            (idx_local,) = ctx.saved_tensors
            idx_mask_local = None
        meta = ctx.comm_meta
        recv_meta = meta["recv_metadata_for_bwd"]
        send_meta = meta["send_metadata_for_bwd"]
        x_local_shape = meta["x_local_shape"]
        x_buffer_shape = meta["x_buffer_shape"]
        need_interval = meta["need_interval"]
        need_start = need_interval[0, 0]
        axis = meta["axis"]
        x_placements = meta["x_placements"]
        output_placements = meta["output_placements"]
        x_global_shape = meta["x_global_shape"]
        own_interval = meta["own_interval"]
        device_mesh_output = meta["device_mesh_output"]

        if device_mesh_output != grad_output.device_mesh:
            raise ValueError(
                f"grad_output device_mesh mismatch: expected {device_mesh_output}, got {grad_output.device_mesh}"
            )
        if output_placements != grad_output.placements:
            raise ValueError(
                f"grad_output placements mismatch: expected {output_placements}, got {grad_output.placements}"
            )

        grad_local = grad_output.to_local().contiguous()

        local_idx = idx_local - need_start.item()
        shape_trailing = x_buffer_shape[axis + 1 :]
        shape_trailing_flat = torch.Size(shape_trailing).numel()

        shape_leading = x_buffer_shape[:axis]
        shape_leading_flat = torch.Size(shape_leading).numel()
        L = x_buffer_shape[axis]
        K = local_idx.shape[-2]
        W = local_idx.shape[-1]

        grad_flat = grad_local.reshape(shape_leading_flat, K * W, shape_trailing_flat)

        # Handle edge case where buffer is empty (all indices masked out)
        if L == 0:
            grad_x_buffer = torch.zeros(
                *shape_leading, 0, *shape_trailing, dtype=grad_local.dtype, device=grad_local.device
            )
        else:
            # For masked indices, clamp to valid buffer range and zero out gradients
            if idx_mask_local is not None:
                local_idx = torch.where(idx_mask_local, local_idx, torch.zeros_like(local_idx))
                mask_flat = idx_mask_local.reshape(shape_leading_flat, K * W, 1).to(grad_flat.dtype)
                grad_flat = grad_flat * mask_flat

            idx_flat = local_idx.reshape(shape_leading_flat, K * W, 1).expand(
                shape_leading_flat, K * W, shape_trailing_flat
            )

            grad_buf = torch.zeros(
                shape_leading_flat, L, shape_trailing_flat, dtype=grad_local.dtype, device=grad_local.device
            )
            grad_buf.scatter_add_(1, idx_flat, grad_flat)

            grad_x_buffer = grad_buf.reshape(*shape_leading, L, *shape_trailing)

        ops = []
        grad_x_local = torch.zeros(x_local_shape, dtype=grad_local.dtype, device=grad_local.device)

        # Backward send (reverse of recv)
        for peer, interval, _shape in recv_meta:
            start_local_buf = interval[0, 0] - need_start.to(interval.device)
            length = interval[0, 1] - interval[0, 0]
            grad_chunk = grad_x_buffer.narrow(axis, start_local_buf.item(), length.item())

            if peer == dist.get_rank():
                # self accumulate
                start_local = interval[0, 0] - own_interval[0, 0]
                target = grad_x_local.narrow(axis, start_local.item(), length.item())
                target.add_(grad_chunk)
            else:
                ops.append(dist.P2POp(dist.isend, grad_chunk.contiguous(), peer))

        # Backward recv (reverse of send)
        bwd_recv_bufs = []
        for peer, start_local, length in send_meta:
            shape = list(x_local_shape)
            shape[axis] = length.item()
            buf = torch.empty(shape, dtype=grad_local.dtype, device=grad_local.device)
            ops.append(dist.P2POp(dist.irecv, buf, peer))
            bwd_recv_bufs.append((buf, start_local, length))

        if ops:
            reqs = dist.batch_isend_irecv(ops)
            for req in reqs:
                req.wait()

        for buf, start_local, length in bwd_recv_bufs:
            target = grad_x_local.narrow(axis, start_local.item(), length.item())
            target.add_(buf)

        grad_x_dtensor = DTensor.from_local(
            grad_x_local,
            grad_output.device_mesh,
            x_placements,
            shape=x_global_shape,
            stride=update_exhaustive_strides(grad_x_local.shape, grad_x_local.stride(), x_global_shape),
        )

        return grad_x_dtensor, None, None, None, None


def distributed_gather(
    x_dtensor: DTensor,
    idx_dtensor: DTensor,
    axis: int = 1,
    are_ids_contiguous: bool = False,
    idx_mask: DTensor | None = None,
) -> DTensor:
    """Distributed 1D gather.

    Args:
        x_dtensor: DTensor with shape ``(*batch, N, *features)``.
        idx_dtensor: DTensor with shape ``(*batch, K, W)`` that provides gather
            indices into the ``N`` dimension of ``x_dtensor``.
        axis: Axis in ``x_dtensor`` corresponding to ``N`` (and to ``K`` in
            ``idx_dtensor``).
        are_ids_contiguous: This is a heuristic for selecting the underlying
            send/recv strategy for performance purpose. Currently only True is supported,
            which means that the idx_dtensor maps to a contiguous block of x along (axis)
            dimensions for all the shards and for all the leading (batch) dimensions.
            When True, the underlying strategy will use the min/max of idx_dtensor to
            compute the needed interval, assuming the resulting buffer to be communicated
            across the ranks is fully (or approximately so) utilized. In the case the
            data inside idx_dtensor doesn't mapped to contiguous blocks, the result will
            still be correct by setting are_ids_contiguous=True but the buffer of x chunks
            communicated will contain a lot of unused elements, making the sed/recv inefficient.
        idx_mask: Optional DTensor with shape ``(*batch, K, W)`` and same device_mesh
            and placements as ``idx_dtensor``. Elements with True indicate valid indices,
            elements with False indicate invalid indices that should be ignored.
    """
    return DistributedGather.apply(x_dtensor, idx_dtensor, axis, are_ids_contiguous, idx_mask)
