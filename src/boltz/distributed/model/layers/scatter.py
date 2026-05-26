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


class DistributedScatterReduce(torch.autograd.Function):
    @staticmethod
    def forward(
        ctx,
        output_size_per_rank: int,
        axis: int,
        idx_dtensor: DTensor,
        src_dtensor: DTensor,
        reduce: str,
        idx_mask: DTensor | None,
        are_ids_contiguous: bool,
    ) -> DTensor:
        """Distributed scatter reduce.

        Scatters values from src into an output tensor at positions specified by idx,
        applying a reduction operation for duplicate indices. The output is initialized
        to zeros (not from a dst tensor).

        Args:
            output_size_per_rank: Size of the output's scatter axis per rank. The global
                output shape will be ``(*batch, output_size_per_rank * num_shards, *features)``.
            axis: Axis corresponding to the scatter dimension in output and src.
            idx_dtensor: DTensor with shape ``(*batch, N_src)`` that provides scatter
                indices into the output's scatter dimension. Values must be in range
                [0, output_size_per_rank * num_shards).
            src_dtensor: DTensor with shape ``(*batch, N_src, *features)`` - source data
                to scatter into output. Must have same placements as idx_dtensor.
            reduce: Reduction operation - "sum" or "mean".
            idx_mask: Optional DTensor with shape ``(*batch, N_src)`` and same device_mesh
                and placements as ``idx_dtensor``. Elements with True indicate valid indices,
                elements with False indicate invalid indices that should be ignored.
            are_ids_contiguous: This is a heuristic for selecting the underlying
                send/recv strategy for performance purpose. Currently only True is supported,
                which means that the idx_dtensor maps to a contiguous block of dst along (axis)
                dimensions for all the shards and for all the leading (batch) dimensions.
                When True, the underlying strategy will use the min/max of idx_dtensor to
                compute the needed interval, assuming the resulting buffer to be communicated
                across the ranks is fully (or approximately so) utilized.

        Returns:
            DTensor with shape ``(*batch, output_size_per_rank * num_shards, *features)``
            containing the scattered and reduced values.
        """
        if not are_ids_contiguous:
            raise NotImplementedError("DistributedScatterReduce currently only supports are_ids_contiguous=True")

        if reduce not in ("sum", "mean"):
            raise ValueError(
                f"reduce must be 'sum' or 'mean', got '{reduce}'. "
                "Other reductions (amax, amin, prod) are not supported due to complex backward requirements."
            )

        if not isinstance(idx_dtensor, DTensor) or not isinstance(src_dtensor, DTensor):
            raise TypeError("idx_dtensor and src_dtensor must be DTensors")

        # Validate shapes
        # idx: (*batch, N_src) - indices into output dimension
        # src: (*batch, N_src, *features) - values to scatter
        # output: (*batch, output_size_per_rank * num_shards, *features)
        batch_dims = idx_dtensor.shape[:axis]
        feature_shape = src_dtensor.shape[axis + 1 :]

        # idx should be (*batch, N_src)
        if idx_dtensor.ndim != axis + 1:
            raise ValueError(f"idx_dtensor should have ndim={axis + 1} for axis={axis}, got ndim={idx_dtensor.ndim}")

        N_src = idx_dtensor.shape[axis]

        # src should be (*batch, N_src, *features)
        if src_dtensor.shape[:axis] != batch_dims:
            raise ValueError(f"Batch dimensions must match: idx {batch_dims} vs src {src_dtensor.shape[:axis]}")
        if src_dtensor.shape[axis] != N_src:
            raise ValueError(f"src axis {axis} size {src_dtensor.shape[axis]} must match idx axis size {N_src}")

        # Validate device_mesh and placements
        mesh = idx_dtensor.device_mesh
        if src_dtensor.device_mesh != mesh:
            raise ValueError("idx and src must be on the same DeviceMesh")
        if src_dtensor.placements != idx_dtensor.placements:
            raise ValueError("idx and src must have identical placements")

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

        placements = idx_dtensor.placements

        ndim_src = src_dtensor.ndim
        if axis < 0:
            axis += ndim_src
        if axis < 0 or axis >= ndim_src:
            raise ValueError(f"axis {axis} out of range for src.ndim={ndim_src}")

        # Identify shard axis on mesh for the scatter dimension
        mesh_dim_axis = None
        for i_mesh_dim, p in enumerate(placements):
            if isinstance(p, Partial):
                raise ValueError("Partial placements are not supported")
            if isinstance(p, Shard):
                if p.dim == axis:
                    mesh_dim_axis = i_mesh_dim
                # Enforce even sharding for idx and src
                if idx_dtensor.shape[p.dim] % mesh.size(i_mesh_dim) != 0:
                    raise ValueError(
                        f"idx_dtensor axis {p.dim} size {idx_dtensor.shape[p.dim]} "
                        f"not evenly divisible by mesh dim {i_mesh_dim}"
                    )
                if src_dtensor.shape[p.dim] % mesh.size(i_mesh_dim) != 0:
                    raise ValueError(
                        f"src_dtensor axis {p.dim} size {src_dtensor.shape[p.dim]} "
                        f"not evenly divisible by mesh dim {i_mesh_dim}"
                    )
                if idx_mask is not None and idx_mask.shape[p.dim] % mesh.size(i_mesh_dim) != 0:
                    raise ValueError(
                        f"idx_mask axis {p.dim} size {idx_mask.shape[p.dim]} "
                        f"not evenly divisible by mesh dim {i_mesh_dim}"
                    )

        if mesh_dim_axis is None:
            raise ValueError(f"Tensors must be sharded along axis {axis}")

        # Compute output shapes
        size_group = mesh.size(mesh_dim_axis)
        output_global_size = output_size_per_rank * size_group
        output_global_shape = batch_dims + (output_global_size,) + feature_shape

        idx_local = idx_dtensor.to_local()
        src_local = src_dtensor.to_local()
        idx_mask_local = idx_mask.to_local() if idx_mask is not None else None
        device = idx_local.device
        cpu_device = torch.device("cpu")

        # Local output shape: use local batch dims (from src_local) + output_size_per_rank + local features
        batch_dims_local = src_local.shape[:axis]
        feature_shape_local = src_local.shape[axis + 1 :]
        output_local_shape = batch_dims_local + (output_size_per_rank,) + feature_shape_local

        # Compute write interval (where our local idx values want to write to)
        if idx_local.numel() > 0:
            if idx_mask_local is not None:
                if idx_mask_local.any():
                    valid_idx = idx_local[idx_mask_local]
                    write_interval = torch.stack(valid_idx.aminmax()).to(dtype=torch.long).unsqueeze(0)  # (1,2)
                    write_interval[:, -1] += 1
                else:
                    write_interval = torch.tensor([[0, 0]], device=device, dtype=torch.long)
            else:
                write_interval = torch.stack(idx_local.aminmax()).to(dtype=torch.long).unsqueeze(0)
                write_interval[:, -1] += 1
        else:
            write_interval = torch.tensor([[0, 0]], device=device, dtype=torch.long)

        write_start = write_interval[0, 0]
        write_end = write_interval[0, 1]

        # Owned chunk interval (this rank's portion of output)
        coord_axis = mesh.get_local_rank(mesh_dim_axis)
        own_start = torch.tensor(coord_axis * output_size_per_rank, device=cpu_device, dtype=torch.long)
        own_end = own_start + output_size_per_rank
        own_interval = torch.stack([own_start, own_end]).unsqueeze(0)  # (1,2)

        # All-gather write intervals along sharded mesh dim (metadata only)
        group_axis = mesh.get_group(mesh_dim_axis)
        write_range = [torch.zeros_like(write_interval) for _ in range(size_group)]
        dist.all_gather(write_range, write_interval, group=group_axis)
        write_range = torch.stack(write_range)  # (size_group, 1, 2)
        write_range_cpu = write_range.cpu()

        ranks_global_on_mesh = mesh.mesh
        my_coords = mesh.get_coordinate()
        index_list_submesh = []
        for dim in range(ranks_global_on_mesh.ndim):
            if dim == mesh_dim_axis:
                index_list_submesh.append(slice(None))
            else:
                index_list_submesh.append(torch.tensor(my_coords[dim], device=cpu_device))
        ranks_global_on_submesh = ranks_global_on_mesh[tuple(index_list_submesh)]  # (size_group,)

        my_rank = mesh.mesh[tuple(my_coords)].item()

        # Compute owned intervals for all peers
        start_peers_own = torch.arange(size_group, device=cpu_device, dtype=torch.long) * output_size_per_rank
        end_peers_own = start_peers_own + output_size_per_rank
        interval_peers_own = torch.stack([start_peers_own, end_peers_own], dim=-1).unsqueeze(1)  # (size_group, 1, 2)

        # Shape info (batch dims are same for output/src/idx due to same placements)
        N_src_local = idx_local.shape[axis]  # local N_src dimension size
        shape_leading_flat = torch.Size(batch_dims_local).numel() if batch_dims_local else 1
        shape_trailing_flat = torch.Size(feature_shape_local).numel() if feature_shape_local else 1

        # Flatten for easier processing
        # src: (*batch, N_src, *features) -> (B, N_src, F)
        # idx: (*batch, N_src) -> (B, N_src)
        src_flat = src_local.reshape(shape_leading_flat, N_src_local, shape_trailing_flat)
        idx_flat = idx_local.reshape(shape_leading_flat, N_src_local)
        if idx_mask_local is not None:
            mask_flat = idx_mask_local.reshape(shape_leading_flat, N_src_local)
        else:
            mask_flat = None

        # SEND PLAN: determine which peers we need to send data to
        # For each peer whose owned_interval overlaps with our write_interval
        if write_start >= write_end:
            send_plan: List[Dict] = []
        else:
            send_plan = get_overlap_from_peers(
                ranks_global_on_submesh, interval_peers_own, write_interval.to(cpu_device)
            )

        # RECV PLAN: determine which peers will send data to us
        # For each peer whose write_interval overlaps with our owned_interval
        recv_plan = get_overlap_from_peers(
            ranks_global_on_submesh,
            write_range_cpu.view(size_group, 1, 2),
            own_interval.view(1, 1, 2),
        )

        # Batch indices tensor for indexing into flattened (B, N_src) tensors
        # batch_indices[b, n] = b, used to track which batch element each element belongs to
        batch_indices = torch.arange(shape_leading_flat, device=device).unsqueeze(1).expand(-1, N_src_local)

        # Phase 1: Prepare data to send and count elements
        # send_counts[i] = number of (idx, src) pairs this rank will send to peer i
        # Used for all-gather so each rank knows how much data to expect from others
        send_counts = torch.zeros(size_group, device=device, dtype=torch.long)
        # send_data_per_peer[peer_global_rank] = (batch_idx, dst_local_idx, src_values)
        #   - batch_idx: (count,) which batch element each pair belongs to
        #   - dst_local_idx: (count,) destination index relative to peer's owned chunk start
        #   - src_values: (count, F) flattened source values to scatter
        send_data_per_peer: Dict[int, Tuple[torch.Tensor, torch.Tensor, torch.Tensor]] = {}

        for item in send_plan:
            peer = item["peer"]
            interval = item["interval"]
            peer_start = interval[0, 0].item()
            peer_end = interval[0, 1].item()

            # Create mask for idx values that fall in peer's interval
            in_peer_interval = (idx_flat >= peer_start) & (idx_flat < peer_end)
            if mask_flat is not None:
                in_peer_interval = in_peer_interval & mask_flat

            count = in_peer_interval.sum().item()
            # Get peer's rank within the group (peer is global rank, need group-local rank)
            peer_coord = dist.get_group_rank(group_axis, peer)
            send_counts[peer_coord] = count

            if count > 0:
                peer_batch_idx = batch_indices[in_peer_interval]  # (count,)

                # idx values adjusted to be relative to peer's owned chunk start
                # Note: peer_start/peer_end are the overlap interval used for filtering,
                # but we need peer's actual owned start for computing local indices
                peer_owned_start = peer_coord * output_size_per_rank
                peer_idx = idx_flat[in_peer_interval] - peer_owned_start  # (count,)

                # src values
                peer_src = src_flat[in_peer_interval]  # (count, F)

                send_data_per_peer[peer] = (peer_batch_idx, peer_idx, peer_src)

        # All-gather counts from all ranks
        all_counts = [torch.zeros_like(send_counts) for _ in range(size_group)]
        dist.all_gather(all_counts, send_counts, group=group_axis)
        all_counts = torch.stack(all_counts)  # (size_group, size_group)

        # all_counts[i, j] = count of elements rank i sends to rank j
        my_coord_in_group = coord_axis
        recv_counts = all_counts[:, my_coord_in_group]  # counts I receive from each rank

        # Phase 2: Exchange actual data
        ops = []
        recv_bufs: Dict[int, Tuple[torch.Tensor, torch.Tensor, torch.Tensor]] = {}

        # Prepare receives (including self-data from send_data_per_peer)
        for item in recv_plan:
            peer = item["peer"]
            if peer == my_rank:
                # Self-data: directly use local data from send_data_per_peer instead of P2P
                # This unifies self-data and received-data processing in the scatter loop below
                if my_rank in send_data_per_peer:
                    recv_bufs[my_rank] = send_data_per_peer[my_rank]
                continue
            # Get peer's rank within the group (peer is global rank, need group-local rank)
            peer_coord = dist.get_group_rank(group_axis, peer)
            count = recv_counts[peer_coord].item()
            if count > 0:
                recv_batch_buf = torch.empty(count, device=device, dtype=torch.long)
                recv_idx_buf = torch.empty(count, device=device, dtype=idx_local.dtype)
                recv_src_buf = torch.empty(count, shape_trailing_flat, device=device, dtype=src_local.dtype)
                ops.append(dist.P2POp(dist.irecv, recv_batch_buf, peer))
                ops.append(dist.P2POp(dist.irecv, recv_idx_buf, peer))
                ops.append(dist.P2POp(dist.irecv, recv_src_buf, peer))
                recv_bufs[peer] = (recv_batch_buf, recv_idx_buf, recv_src_buf)

        # Prepare sends (skip self since we handled it above)
        for item in send_plan:
            peer = item["peer"]
            if peer == my_rank:
                continue
            if peer in send_data_per_peer:
                peer_batch_idx, peer_idx, peer_src = send_data_per_peer[peer]
                ops.append(dist.P2POp(dist.isend, peer_batch_idx.contiguous(), peer))
                ops.append(dist.P2POp(dist.isend, peer_idx.contiguous(), peer))
                ops.append(dist.P2POp(dist.isend, peer_src.contiguous(), peer))

        # Execute P2P
        if ops:
            reqs = dist.batch_isend_irecv(ops)
            for req in reqs:
                req.wait()

        # Local scatter_reduce operation
        # Initialize output with zeros (no dst tensor to clone from)
        out_local = torch.zeros(output_local_shape, dtype=src_local.dtype, device=device)

        # Linearize output for scatter_reduce: (B, L, F) -> (B*L, F)
        # This allows us to use a single scatter_reduce_ call instead of looping over batches.
        # The key insight: scatter_reduce_ only indexes along one dimension, but we need to
        # index by (batch, position) pairs. By linearizing to (B*L, F), we can compute
        # linear_idx = batch_idx * L + position_idx to get a single index into the B*L dimension.
        out_linear = out_local.reshape(shape_leading_flat * output_size_per_rank, shape_trailing_flat)

        # For "mean" reduction, we need to track counts for each output position
        # count[j] = number of src values scattered to position j
        # This is needed for both the forward (to compute mean) and backward (grad = grad_out / count)
        if reduce == "mean":
            # Initialize counts to 0 (no include_self since there's no dst)
            count_linear = torch.zeros(
                shape_leading_flat * output_size_per_rank, 1, dtype=out_local.dtype, device=device
            )

        # Process all data uniformly (self-data + received data from peers)
        # We always use scatter_add_ since we're accumulating values and dividing for mean at the end
        # Both are stored in recv_bufs as (batch_idx, local_dst_idx, src_values)
        for peer, (recv_batch_idx, recv_idx, recv_src) in recv_bufs.items():
            recv_idx_long = recv_idx.to(torch.long)

            # Compute linearized index: maps (batch, position) -> single index in B*L
            linear_idx = recv_batch_idx * output_size_per_rank + recv_idx_long
            linear_idx_expanded = linear_idx.unsqueeze(-1).expand(-1, shape_trailing_flat)

            out_linear.scatter_add_(0, linear_idx_expanded, recv_src)

            # Accumulate counts for "mean"
            if reduce == "mean":
                ones = torch.ones(linear_idx.shape[0], 1, dtype=out_local.dtype, device=device)
                count_linear.scatter_add_(0, linear_idx.unsqueeze(-1), ones)

        # For "mean", divide the accumulated sum by counts
        # Positions with count=0 are already 0 from initialization, so 0/1=0 naturally
        count_local = None
        if reduce == "mean":
            out_linear = out_linear / count_linear.clamp(min=1).expand(-1, shape_trailing_flat)
            count_local = count_linear.reshape(shape_leading_flat, output_size_per_rank, 1)
            count_local = count_local.expand(-1, -1, shape_trailing_flat).reshape(output_local_shape)

        # Reshape back to original structure
        out_local = out_linear.reshape(output_local_shape)

        strides_out = update_exhaustive_strides(out_local.shape, out_local.stride(), output_global_shape)

        out_dtensor = DTensor.from_local(out_local, mesh, placements, shape=output_global_shape, stride=strides_out)

        # Save for backward
        # For "mean", we also save count_local to correctly scale gradients
        tensors_to_save = [idx_local, src_local]
        if idx_mask_local is not None:
            tensors_to_save.append(idx_mask_local)
        if count_local is not None:
            tensors_to_save.append(count_local)
        ctx.save_for_backward(*tensors_to_save)
        ctx.has_mask = idx_mask_local is not None
        ctx.has_count = count_local is not None
        ctx.reduce = reduce
        ctx.comm_meta = {
            "axis": axis,
            "placements": placements,
            "output_global_shape": output_global_shape,
            "output_size_per_rank": output_size_per_rank,
            "src_global_shape": src_dtensor.shape,
            "output_global_stride": strides_out,
            "src_global_stride": src_dtensor.stride(),
            "own_interval": own_interval,
            "device_mesh": mesh,
            "mesh_dim_axis": mesh_dim_axis,
        }

        return out_dtensor

    @staticmethod
    def backward(ctx, grad_output: DTensor) -> Tuple[None, None, None, DTensor, None, None, None]:
        """Backward pass for DistributedScatterReduce.

        The backward of scatter_reduce is essentially a gather operation:
        - For "sum": grad_src = gather(grad_output, idx)
        - For "mean": grad_src = gather(grad_output / count, idx)

        Note on reusing distributed_gather:
        While the backward is semantically a gather operation, we don't directly call
        distributed_gather here because:
        1. Shape mismatch: distributed_gather expects idx shape (*batch, K, W), but our
           idx has shape (*batch, N_src) - a 1D index per batch element.
        2. For "mean", we need to gather from (grad_output / count), not raw grad_output.
           The count tensor is local to each rank's output shard, so we need to apply the
           division before sending chunks to other ranks during the gather.

        Instead, we implement the gather communication pattern directly, which also allows
        us to handle the count scaling for "mean" reduction efficiently by using the
        count-scaled gradient (grad_gather_source = grad_local / count_local) as the source.

        Returns gradients for: (output_size_per_rank, axis, idx, src, reduce, idx_mask, are_ids_contiguous)
        Only src requires a gradient; others return None.
        """
        # Unpack saved tensors based on what was saved
        saved = ctx.saved_tensors
        idx_local = saved[0]
        src_local = saved[1]
        tensor_idx = 2

        idx_mask_local = None
        if ctx.has_mask:
            idx_mask_local = saved[tensor_idx]
            tensor_idx += 1

        count_local = None
        if ctx.has_count:
            count_local = saved[tensor_idx]

        reduce = ctx.reduce
        meta = ctx.comm_meta
        axis = meta["axis"]
        placements = meta["placements"]
        output_size_per_rank = meta["output_size_per_rank"]
        src_global_shape = meta["src_global_shape"]
        src_global_stride = meta["src_global_stride"]
        own_interval = meta["own_interval"]
        mesh = meta["device_mesh"]
        mesh_dim_axis = meta["mesh_dim_axis"]

        grad_local = grad_output.to_local().contiguous()
        device = grad_local.device
        cpu_device = torch.device("cpu")

        # For "mean", we need to scale by 1/count because:
        #   output[j] = sum of src[i] where idx[i]==j / count[j]
        #   grad_src[i] = grad_output[idx[i]] * d(output[idx[i]])/d(src[i]) = grad_output[idx[i]] / count[idx[i]]
        # We use the count-scaled gradients so that gathered values are automatically divided
        # by count[idx] from the owning rank
        if reduce == "mean" and count_local is not None:
            # Clamp count to avoid division by zero for positions that received no values
            count_clamped = count_local.clamp(min=1)
            grad_gather_source = grad_local / count_clamped
        else:
            grad_gather_source = grad_local

        # grad_src requires gathering from grad_gather_source at idx positions
        # This is essentially the same communication pattern as DistributedGather

        # Compute need interval (what we need from grad_output for our grad_src)
        if idx_local.numel() > 0:
            if idx_mask_local is not None:
                if idx_mask_local.any():
                    valid_idx = idx_local[idx_mask_local]
                    need_interval = torch.stack(valid_idx.aminmax()).to(dtype=torch.long).unsqueeze(0)
                    need_interval[:, -1] += 1
                else:
                    need_interval = torch.tensor([[0, 0]], device=device, dtype=torch.long)
            else:
                need_interval = torch.stack(idx_local.aminmax()).to(dtype=torch.long).unsqueeze(0)
                need_interval[:, -1] += 1
        else:
            need_interval = torch.tensor([[0, 0]], device=device, dtype=torch.long)

        need_start = need_interval[0, 0]
        need_end = need_interval[0, 1]
        need_start_cpu = need_start.to(cpu_device)

        own_start = own_interval[0, 0]

        # All-gather need intervals
        group_axis = mesh.get_group(mesh_dim_axis)
        size_group = mesh.size(mesh_dim_axis)
        need_range = [torch.zeros_like(need_interval) for _ in range(size_group)]
        dist.all_gather(need_range, need_interval, group=group_axis)
        need_range = torch.stack(need_range)
        need_range_cpu = need_range.cpu()

        ranks_global_on_mesh = mesh.mesh
        my_coords = mesh.get_coordinate()
        index_list_submesh = []
        for dim in range(ranks_global_on_mesh.ndim):
            if dim == mesh_dim_axis:
                index_list_submesh.append(slice(None))
            else:
                index_list_submesh.append(torch.tensor(my_coords[dim], device=cpu_device))
        ranks_global_on_submesh = ranks_global_on_mesh[tuple(index_list_submesh)]

        my_rank = mesh.mesh[tuple(my_coords)].item()

        # RECEIVE PLAN for backward (gather pattern)
        if need_start >= need_end:
            needed_chunks: List[Dict] = []
        else:
            start_peers_own = torch.arange(size_group, device=cpu_device, dtype=torch.long) * output_size_per_rank
            end_peers_own = start_peers_own + output_size_per_rank
            interval_peers_own = torch.stack([start_peers_own, end_peers_own], dim=-1).unsqueeze(1)

            needed_chunks = get_overlap_from_peers(
                ranks_global_on_submesh, interval_peers_own, need_interval.to(cpu_device)
            )

        ops = []
        recv_bufs = {}

        for item in needed_chunks:
            peer = item["peer"]
            interval = item["interval"]
            start_global = interval[0, 0]
            length = interval[0, 1] - interval[0, 0]

            shape = list(grad_gather_source.shape)
            shape[axis] = length.item()
            buf = torch.empty(shape, dtype=grad_gather_source.dtype, device=device)

            if peer == my_rank:
                start_local = start_global - own_start
                buf.copy_(grad_gather_source.narrow(axis, start_local.item(), length.item()))
                recv_bufs[peer] = buf
            else:
                ops.append(dist.P2POp(dist.irecv, buf, peer))
                recv_bufs[peer] = buf

        # SEND PLAN for backward
        send_chunks = get_overlap_from_peers(
            ranks_global_on_submesh,
            need_range_cpu.view(size_group, 1, 2),
            own_interval.view(1, 1, 2),
        )

        for item in send_chunks:
            peer = item["peer"]
            if peer == my_rank:
                continue
            interval = item["interval"]
            start_global = interval[0, 0]
            length = interval[0, 1] - interval[0, 0]
            start_local = start_global - own_start
            chunk = grad_gather_source.narrow(axis, start_local.item(), length.item()).contiguous()
            ops.append(dist.P2POp(dist.isend, chunk, peer))

        if ops:
            reqs = dist.batch_isend_irecv(ops)
            for req in reqs:
                req.wait()

        # Assemble gradient buffer
        if need_start >= need_end:
            buffer_shape = list(grad_gather_source.shape)
            buffer_shape[axis] = 0
            grad_buffer = torch.empty(buffer_shape, dtype=grad_gather_source.dtype, device=device)
        else:
            buffer_shape = list(grad_gather_source.shape)
            buffer_shape[axis] = (need_end - need_start).item()
            grad_buffer = torch.zeros(buffer_shape, dtype=grad_gather_source.dtype, device=device)

            for item in needed_chunks:
                peer = item["peer"]
                interval = item["interval"]
                buf = recv_bufs.get(peer)
                if buf is None:
                    raise RuntimeError(f"Missing recv buffer for peer {peer}")

                start_global = interval[0, 0]
                length = interval[0, 1] - interval[0, 0]
                start_local_buf = start_global - need_start_cpu
                target = grad_buffer.narrow(axis, start_local_buf.item(), length.item())
                target.copy_(buf)

        # Gather grad_src from grad_buffer
        idx_local_buffer = idx_local - need_start.item()

        shape_trailing = grad_buffer.shape[axis + 1 :]
        shape_trailing_flat = torch.Size(shape_trailing).numel() if shape_trailing else 1

        shape_leading = grad_buffer.shape[:axis]
        shape_leading_flat = torch.Size(shape_leading).numel()
        L = grad_buffer.shape[axis]
        N_src_local = idx_local.shape[axis]

        if L == 0:
            grad_src_local = torch.zeros_like(src_local)
        else:
            if idx_mask_local is not None:
                idx_local_buffer = torch.where(idx_mask_local, idx_local_buffer, torch.zeros_like(idx_local_buffer))

            grad_flat = grad_buffer.reshape(shape_leading_flat, L, shape_trailing_flat)
            idx_flat = idx_local_buffer.reshape(shape_leading_flat, N_src_local)

            batch_idx = torch.arange(shape_leading_flat, device=device).reshape(shape_leading_flat, 1)
            grad_src_flat = grad_flat[batch_idx, idx_flat, :]  # (B, N_src, F)

            if idx_mask_local is not None:
                mask_flat = idx_mask_local.reshape(shape_leading_flat, N_src_local, 1).to(grad_src_flat.dtype)
                grad_src_flat = grad_src_flat * mask_flat

            if shape_trailing:
                grad_src_local = grad_src_flat.reshape(*shape_leading, N_src_local, *shape_trailing)
            else:
                grad_src_local = grad_src_flat.reshape(*shape_leading, N_src_local)

        # Construct grad_src DTensor using saved strides from forward pass
        grad_src_dtensor = DTensor.from_local(
            grad_src_local.contiguous(),
            grad_output.device_mesh,
            placements,
            shape=src_global_shape,
            stride=src_global_stride,
        )

        # Return: (output_size_per_rank, axis, idx, src, reduce, idx_mask, are_ids_contiguous)
        # Only src requires gradient
        return None, None, None, grad_src_dtensor, None, None, None


def distributed_scatter_reduce(
    output_size_per_rank: int,
    axis: int,
    idx: DTensor,
    src: DTensor,
    reduce: str,
    idx_mask: DTensor | None = None,
    are_ids_contiguous: bool = False,
) -> DTensor:
    """Distributed scatter reduce.

    Scatters values from src into an output tensor at positions specified by idx,
    applying a reduction operation for duplicate indices. The output is initialized
    to zeros.

    This is a P2P-based alternative to:
    ``torch.zeros(output_shape).scatter_reduce_(axis, idx.full_tensor(), src.full_tensor(), reduce)``

    Args:
        output_size_per_rank: Size of the output's scatter axis per rank. The global
            output shape will be ``(*batch, output_size_per_rank * num_shards, *features)``.
        axis: Axis corresponding to the scatter dimension in output and src.
        idx: DTensor with shape ``(*batch, N_src)`` that provides scatter indices
            into the output's scatter dimension. Values must be in range
            [0, output_size_per_rank * num_shards).
        src: DTensor with shape ``(*batch, N_src, *features)`` - source data to scatter.
            Must have same placements as idx.
        reduce: Reduction operation - "sum" or "mean".
        idx_mask: Optional DTensor with shape ``(*batch, N_src)`` and same device_mesh
            and placements as ``idx``. Elements with True indicate valid indices,
            elements with False indicate invalid indices that should be ignored.
        are_ids_contiguous: Heuristic for send/recv strategy. Currently only True is supported.

    Returns:
        DTensor: Result with shape ``(*batch, output_size_per_rank * num_shards, *features)``
        and same placements as idx.
    """
    return DistributedScatterReduce.apply(output_size_per_rank, axis, idx, src, reduce, idx_mask, are_ids_contiguous)
