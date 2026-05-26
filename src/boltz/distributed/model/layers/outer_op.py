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


from enum import Enum, auto
from typing import Optional

import torch
from torch.distributed.tensor import DTensor, Replicate, Shard

from boltz.distributed.comm import TransposeComm
from boltz.distributed.model.layers.redistribute_transpose_without_dtensor import transpose_then_redistribute
from boltz.distributed.utils import LayoutRightMap


class OuterOp(Enum):
    SUM = auto()
    SUBTRACT = auto()
    EQUAL = auto()
    BITAND = auto()
    PROD = auto()
    CDIST = auto()  # Special case for pairwise distance computation


class DistributedOuterOp(torch.autograd.Function):
    """Custom autograd function for distributed outer-[add, subtract, equal, bitand] operations.

    The outer operation assumes the input tensor is sharded along axis 0 of the
    transpose_comm's 2d group grid and also replicated along axis 1 of the same grid:

    input_expanded + input_expanded.transpose(axis, axis + 1)
    return binary_op(input_expanded, input_expanded.transpose(axis, axis + 1))

    Currently, only add, subtract, equal and bitand operations are supported.
    """

    @staticmethod
    def forward(
        ctx,
        input: torch.Tensor,
        input_t: torch.Tensor | None,
        op: OuterOp,
        axis: int,
        transpose_comm: TransposeComm,
        group_replicate: torch.distributed.ProcessGroup,
    ) -> torch.Tensor:
        """Forward pass for DistributedOuterOp.

        Args:
            ctx: Context object to save information for backward pass
            input: Input tensor for outer operation. This is assumed sharded
                   along axis 0 of the transpose_comm's 2d group grid and also
                   replicated along axis 1 of the same grid.
            input_t: Second input tensor for outer operation to be transposed across 2d group grid. Will use `input` if not provided with None.
            op: the binary operation to perform
            axis: Axis along which to perform the outer op
            transpose_comm: Communication object for distributed operations
            group_replicate: Process group for input's replication across ranks

        Returns:
            Tensor with outer op computed

        Raises:
            ValueError: If ranks are inconsistent with the expected process group configuration
        """
        rank_replicate = torch.distributed.get_rank(group_replicate)
        rank_global = transpose_comm.global_rank
        if rank_replicate < 0:
            raise ValueError(
                f"global rank {rank_global} doesn't belong to group_replicate as "
                f"get_rank(group_replicate) returned {rank_replicate}"
            )
        if rank_replicate != transpose_comm.rank_coords[1]:
            raise ValueError(
                f"global rank {rank_global} is not along the input tensor replicating axis, "
                f"which is assumed axis 1 of the transpose_comm's 2d grid, as its rank in the "
                f"grid is {transpose_comm.rank_coords[1]} but the group_replicate's rank is {rank_replicate}"
            )
        ctx.transpose_comm = transpose_comm
        ctx.group_replicate = group_replicate
        ctx.axis = axis
        ctx.op = op
        ctx.is_symmetric = input_t is None

        if input_t is None:
            input_t = input

        input_expanded = input.unsqueeze(axis + 1)
        input_expanded_t = input_t.unsqueeze(axis + 1)
        transposed = transpose_then_redistribute(input_expanded_t, axis, axis + 1, transpose_comm)
        if op == OuterOp.SUM:
            output = input_expanded + transposed
        elif op == OuterOp.SUBTRACT:
            output = input_expanded - transposed
        elif op == OuterOp.EQUAL:
            # boolean output can't be backpropagated but
            # if we were to output a float equivalent, we
            # save the output as a mask to be used in
            # the backward pass
            output = input_expanded == transposed
            ctx.mark_non_differentiable(output)
        elif op == OuterOp.BITAND:
            if input_expanded.dtype.is_floating_point or transposed.dtype.is_floating_point:
                raise ValueError("input_expanded and transposed must be boolean tensors")
            # bitwise AND operation can't be backpropagated
            output = input_expanded & transposed
            ctx.mark_non_differentiable(output)
        else:
            raise ValueError(f"Unsupported operation: {op}")
        return output

    @staticmethod
    def backward(
        ctx, grad_output: torch.Tensor
    ) -> tuple[torch.Tensor | None, torch.Tensor | None, None, None, None, None]:
        """Backward pass for DistributedOuterOp.

        Args:
            ctx: Context object with saved information from forward pass
            grad_output: Gradient tensor from downstream layers

        Returns:
            Tuple containing the gradient for input tensor and None for other parameters
        """
        transpose_comm = ctx.transpose_comm
        group_replicate = ctx.group_replicate
        axis = ctx.axis
        op = ctx.op
        if op == OuterOp.EQUAL or op == OuterOp.BITAND:
            # If EQUAL op and the forward were to output float instead of bool mask
            # then we can use the saved output as mask applied on grad_output here
            # e.g.,
            # if op == OuterOp.EQUAL:
            #     mask = ctx.saved_tensors
            #     grad_output = grad_output * mask
            # BITAND also produces non-differentiable output
            return None, None, None, None, None, None

        # grad on right summand
        grad_transposed = grad_output.sum(dim=axis, keepdim=True).transpose(axis, axis + 1).contiguous()
        if op == OuterOp.SUBTRACT:
            grad_transposed = -grad_transposed
        grad_transposed_recv = transpose_comm.enqueue_to_dispatch(grad_transposed)
        # grad on left summand, which always retain the positive sign
        grad_input_expanded = grad_output.sum(dim=axis + 1, keepdim=True)
        transpose_comm.wait_until_finished()

        # perform allreduce to get the row- and column-wise contributions
        if ctx.is_symmetric:
            grad_input = (grad_input_expanded + grad_transposed_recv).squeeze(dim=axis + 1)
            torch.distributed.all_reduce(grad_input, op=torch.distributed.ReduceOp.SUM, group=group_replicate)
            return grad_input, None, None, None, None, None
        else:
            torch.distributed.all_reduce(grad_input_expanded, op=torch.distributed.ReduceOp.SUM, group=group_replicate)
            torch.distributed.all_reduce(grad_transposed_recv, op=torch.distributed.ReduceOp.SUM, group=group_replicate)
            grad_transposed_recv = grad_transposed_recv.squeeze(dim=axis + 1)
            grad_input_expanded = grad_input_expanded.squeeze(dim=axis + 1)
            return grad_input_expanded, grad_transposed_recv, None, None, None, None


class DistributedCDist(torch.autograd.Function):
    """Custom CP autograd function for torch.cdist.

    Currently supports the default computation, which is the L2 norm.
    Currently supports self-distance.
    """

    @staticmethod
    def forward(
        ctx,
        input_array: torch.Tensor,
        transpose_comm: TransposeComm,
        group_replicate: torch.distributed.ProcessGroup,
    ) -> torch.Tensor:
        """Forward pass for DistributedCDist.

        Args:
            ctx: Context object to save information for backward pass
            input_array: Input tensor for outer operation. This is assumed sharded
                   along axis 0 of the transpose_comm's 2d group grid and also
                   replicated along axis 1 of the same grid.
                   The input tensors sharding dimension is expected to be (-2).
            transpose_comm: Communication object for distributed operations
            group_replicate: Process group for input's replication across ranks

        Returns:
            Tensor with outer op computed

        Raises:
            ValueError: If ranks are inconsistent with the expected process group configuration
        """
        rank_replicate = torch.distributed.get_rank(group_replicate)
        rank_global = transpose_comm.global_rank
        if rank_replicate < 0:
            raise ValueError(
                f"global rank {rank_global} doesn't belong to group_replicate as "
                f"get_rank(group_replicate) returned {rank_replicate}"
            )
        ctx.transpose_comm = transpose_comm
        ctx.group_replicate = group_replicate
        axis = len(input_array.shape) - 2
        ctx.axis = axis
        input_expanded = input_array.unsqueeze(axis + 1)
        transposed = transpose_then_redistribute(input_expanded, axis, axis + 1, transpose_comm)
        # Pairwise distance calculation
        # diff is size (B, N, N, D)
        diff = transposed - input_expanded
        diff_sq = diff * diff
        diff_sum = diff_sq.sum(dim=-1)
        output = diff_sum.sqrt()
        if input_array.requires_grad:
            ctx.save_for_backward(diff)
        return output

    @staticmethod
    def backward(ctx, grad_output: torch.Tensor) -> tuple[torch.Tensor | None, None, None]:
        """Backward pass for DistributedOuterOp.

        Args:
            ctx: Context object with saved information from forward pass
            grad_output: Gradient tensor from downstream layers

        Returns:
            Tuple containing the gradient for input tensor and None for other parameters
        """
        if ctx.needs_input_grad[0] is False:
            return None, None, None
        transpose_comm = ctx.transpose_comm
        group_replicate = ctx.group_replicate
        (diff,) = ctx.saved_tensors
        # Dists is recomputed to save memory.
        diff_sq = diff * diff
        diff_sum = diff_sq.sum(dim=-1)
        dists = diff_sum.sqrt()
        # grad is (B, N, N)
        # grad transposed is (B, N, N)
        grad_transposed = transpose_then_redistribute(grad_output, ctx.axis, ctx.axis + 1, transpose_comm)
        dists = dists + 1e-8
        # (B, N, N, D)
        diff_over_dists = diff / dists.unsqueeze(-1)
        grad_term = diff_over_dists * (grad_output + grad_transposed).unsqueeze(-1)
        # (B, N, N, D) to (B, N, D)
        grad_input = -grad_term.sum(dim=2)
        # Sum reduction across ranks - this is an extension of the above sum.
        torch.distributed.all_reduce(grad_input, op=torch.distributed.ReduceOp.SUM, group=group_replicate)
        return grad_input, None, None


def distributed_cdist(
    input_array: torch.Tensor,
    transpose_comm: TransposeComm | None = None,
    group_replicate: torch.distributed.ProcessGroup | None = None,
):
    """Performs cdist operation, sharded if communication objects are passed, serial if not.

    Args:
        input_array: Input tensor for outer operation. This is assumed sharded
               along axis 0 of the transpose_comm's 2d group grid and also
               replicated along axis 1 of the same grid.
               The input tensors sharding dimension is expected to be (-2).
        transpose_comm: Communication object for distributed operations
        group_replicate: Process group for input's replication across ranks
    """
    if (transpose_comm is None) != (group_replicate is None):
        raise ValueError("transpose_comm and group_replicate must both be provided or both be None")
    if transpose_comm is None:
        return torch.cdist(input_array, input_array)
    return DistributedCDist.apply(input_array, transpose_comm, group_replicate)


class _ReplicateToShardOuterOp(torch.autograd.Function):
    """DTensor version of DistributedOuterOp for outer-[add, subtract, equal, bitand, cdist] operations with R2S transformation.

    The outer operation assumes the input DTensors have placements of (Shard(0), Shard(1), Replicate())
    and produces output DTensors with placements of (Shard(0), Shard(1), Shard(2)).

    The Replicate() placement in the input is used as an extra buffer for the result's Shard(2) placement,
    enabling the R2S (Replicate -> Shard) transformation during the outer operation.

    For most operations:
    input_expanded + input_expanded.transpose(axis, axis + 1)
    return binary_op(input_expanded, input_expanded.transpose(axis, axis + 1))

    For CDIST operation:
    diff = input_expanded.transpose(axis, axis + 1) - input_expanded
    return (diff * diff).sum(dim=-1).sqrt()

    Currently, add, subtract, equal, bitand, and cdist operations are supported.
    """

    @staticmethod
    def forward(
        ctx,
        input: DTensor,
        input_t: DTensor | None,
        op: OuterOp,
        axis: int,
        transpose_comm: TransposeComm,
    ) -> DTensor:
        """Forward pass for _ReplicateToShardOuterOp.

        Args:
            ctx: Context object to save information for backward pass
            input: Input DTensor for outer operation. Must have placements (Shard(0), Shard(1), Replicate())
            input_t: Second input DTensor for outer operation to be transposed across 2d group grid.
                    Will use `input` if not provided with None.
            op: the binary operation to perform
            axis: Axis along which to perform the outer op. Must be one of the Shard placement dimensions (0 or 1).
            transpose_comm: Communication object for distributed operations

        Returns:
            DTensor with outer op computed

        Raises:
            TypeError: If inputs are not DTensors
            ValueError: If device meshes don't match, placements are incorrect, axis is invalid, or uneven sharding detected
        """
        # Type checking
        if not isinstance(input, DTensor):
            raise TypeError(f"input must be DTensor, got {type(input)}")

        if input_t is not None and not isinstance(input_t, DTensor):
            raise TypeError(f"input_t must be DTensor or None, got {type(input_t)}")

        # Check required placements for input: (Shard(0), Shard(1), Replicate())
        input_placements = (Shard(0), Shard(1), Replicate())

        if input.placements != input_placements:
            raise ValueError(f"input must have placements {input_placements}. Got {input.placements}")

        # Get device mesh from input
        device_mesh = input.device_mesh

        # Check that axis is one of the Shard placement dimensions
        shard_dims = [0, 1]  # From Shard(0) and Shard(1)
        if axis not in shard_dims:
            raise ValueError(f"axis must be one of the Shard placement dimensions {shard_dims}. Got {axis}")

        # Check for uneven sharding in input tensor
        if input.shape[0] % device_mesh.shape[0] != 0:
            raise ValueError(
                f"Uneven sharding detected: input tensor dimension 0 of size {input.shape[0]} "
                f"is not evenly divisible by device mesh dimension 0 of size {device_mesh.shape[0]}"
            )

        if input.shape[1] % device_mesh.shape[1] != 0:
            raise ValueError(
                f"Uneven sharding detected: input tensor dimension 1 of size {input.shape[1]} "
                f"is not evenly divisible by device mesh dimension 1 of size {device_mesh.shape[1]}"
            )

        # Set input_t to input if None (symmetric case)
        if input_t is None:
            input_t = input
            is_symmetric = True
        else:
            is_symmetric = False
            # Check device mesh compatibility
            if input_t.device_mesh != input.device_mesh:
                raise ValueError(
                    f"input_t device_mesh mismatch: expected {input.device_mesh}, got {input_t.device_mesh}"
                )
            if input_t.placements != input_placements:
                raise ValueError(f"input_t must have placements {input_placements}. Got {input_t.placements}")
            if input_t.shape != input.shape:
                raise ValueError(f"input and input_t must have the same shape. Got {input.shape} and {input_t.shape}")

        # Infer group_replicate from device mesh (axis 2 corresponds to Replicate())
        group_replicate = device_mesh.get_group(2)

        # Validate ranks (adapted from original DistributedOuterOp)
        rank_replicate = torch.distributed.get_rank(group_replicate)
        rank_global = transpose_comm.global_rank
        if rank_replicate < 0:
            raise ValueError(
                f"global rank {rank_global} doesn't belong to group_replicate as "
                f"get_rank(group_replicate) returned {rank_replicate}"
            )
        if rank_replicate != transpose_comm.rank_coords[1]:
            raise ValueError(
                f"global rank {rank_global} is not along the input tensor replicating axis, "
                f"which is assumed axis 1 of the transpose_comm's 2d grid, as its rank in the "
                f"grid is {transpose_comm.rank_coords[1]} but the group_replicate's rank is {rank_replicate}"
            )

        # Define output placements: (Shard(0), Shard(1), Shard(2)) - R2S transformation
        output_placements = (Shard(0), Shard(1), Shard(2))

        # Save context for backward pass
        ctx.device_mesh = device_mesh
        ctx.input_placements = input_placements
        ctx.output_placements = output_placements
        ctx.transpose_comm = transpose_comm
        ctx.group_replicate = group_replicate
        ctx.axis = axis
        ctx.op = op
        ctx.is_symmetric = is_symmetric
        ctx.input_shape = input.shape
        ctx.input_stride = input.stride()
        ctx.input_t_shape = input_t.shape
        ctx.input_t_stride = input_t.stride()

        # Extract local tensors
        input_local = input.to_local()
        input_t_local = input_t.to_local()

        # Perform the outer operation computation (adapted from original)
        input_expanded = input_local.unsqueeze(axis + 1)
        input_expanded_t = input_t_local.unsqueeze(axis + 1)
        transposed = transpose_then_redistribute(input_expanded_t, axis, axis + 1, transpose_comm)

        if op == OuterOp.SUM:
            output_local = input_expanded + transposed
        elif op == OuterOp.SUBTRACT:
            output_local = input_expanded - transposed
        elif op == OuterOp.EQUAL:
            # boolean output can't be backpropagated but
            # if we were to output a float equivalent, we
            # save the output as a mask to be used in
            # the backward pass
            output_local = input_expanded == transposed
            ctx.mark_non_differentiable(output_local)
        elif op == OuterOp.BITAND:
            # bitwise AND operation can't be backpropagated
            output_local = input_expanded & transposed
            ctx.mark_non_differentiable(output_local)
        elif op == OuterOp.PROD:
            output_local = input_expanded * transposed
            # Save operands for backward: d(a*b)/da = b, d(a*b)/db = a
            if input.requires_grad:
                ctx.save_for_backward(input_expanded, transposed)
        elif op == OuterOp.CDIST:
            # Pairwise distance calculation: L2 norm of difference
            diff_local = input_expanded - transposed
            output_local = (diff_local * diff_local).sum(dim=-1).sqrt()
            # Save diff for backward pass
            if input.requires_grad:
                # this is the gradient of the output with respect to the difference
                # as a prefactor to be multiplied with the downstream gradient
                d_output_d_diff_local = diff_local / (output_local.unsqueeze(-1) + torch.finfo(output_local.dtype).tiny)
                ctx.save_for_backward(d_output_d_diff_local)

        # Compute output shape and stride
        if op == OuterOp.CDIST:
            # CDIST output has reduced the last dimension
            shape_output = input.shape[:2] + (input.shape[1],)
        else:
            # Outer operation adds a dimension: (B, N, D, ...) -> (B, N, N, D, ...) for most ops
            shape_output = input.shape[:2] + (input.shape[1],) + input.shape[2:]

        # Use LayoutRightMap for the output shape
        layout_right = LayoutRightMap(shape_output)
        strides_output = layout_right.strides

        # Convert back to DTensor with output placements (R2S transformation: Replicate -> Shard)
        output = DTensor.from_local(
            output_local, device_mesh, output_placements, shape=shape_output, stride=strides_output
        )

        return output

    @staticmethod
    def backward(ctx, grad_output: DTensor) -> tuple[DTensor | None, DTensor | None, None, None, None]:
        """Backward pass for _ReplicateToShardOuterOp.

        Args:
            ctx: Context object with saved information from forward pass
            grad_output: Gradient DTensor from downstream layers

        Returns:
            Tuple containing the gradient DTensors for input and input_t, and None for other parameters
        """
        # Sanity check grad_output DTensor metadata
        if not isinstance(grad_output, DTensor):
            raise TypeError(f"grad_output must be DTensor, got {type(grad_output)}")

        if grad_output.device_mesh != ctx.device_mesh:
            raise ValueError(
                f"grad_output device_mesh mismatch: expected {ctx.device_mesh}, got {grad_output.device_mesh}"
            )

        if grad_output.placements != ctx.output_placements:
            raise ValueError(
                f"grad_output placements mismatch: expected {ctx.output_placements}, got {grad_output.placements}"
            )

        # Extract local gradient
        grad_output_local = grad_output.to_local()

        # Backward pass logic (adapted from original DistributedOuterOp)
        transpose_comm = ctx.transpose_comm
        group_replicate = ctx.group_replicate
        axis = ctx.axis
        op = ctx.op

        if op == OuterOp.EQUAL or op == OuterOp.BITAND:
            # If EQUAL op and the forward were to output float instead of bool mask
            # then we can use the saved output as mask applied on grad_output here
            # e.g.,
            # if op == OuterOp.EQUAL:
            #     mask = ctx.saved_tensors
            #     grad_output = grad_output * mask
            # BITAND also produces non-differentiable output
            return None, None, None, None, None
        elif op == OuterOp.PROD:
            # d(a*b)/da = b, d(a*b)/db = a — each gradient term uses a different multiplier,
            # so we handle the full reduction here rather than falling through to the common path.
            (input_expanded_local, transposed_local) = ctx.saved_tensors
            # local grad for left operand (input): multiply by transposed, reduce over broadcast dim
            grad_input_expanded = (grad_output_local * transposed_local).sum(dim=axis + 1, keepdim=True)
            # local grad for right operand (input_t): multiply by input_expanded, reduce then transpose back
            grad_transposed = (
                (grad_output_local * input_expanded_local)
                .sum(dim=axis, keepdim=True)
                .transpose(axis, axis + 1)
                .contiguous()
            )
            grad_transposed_recv = transpose_comm.enqueue_to_dispatch(grad_transposed)
            transpose_comm.wait_until_finished()

            if ctx.is_symmetric:
                grad_input_local = (grad_input_expanded + grad_transposed_recv).squeeze(dim=axis + 1)
                torch.distributed.all_reduce(grad_input_local, op=torch.distributed.ReduceOp.SUM, group=group_replicate)
                grad_input = DTensor.from_local(
                    grad_input_local,
                    ctx.device_mesh,
                    ctx.input_placements,
                    shape=ctx.input_shape,
                    stride=ctx.input_stride,
                )
                return grad_input, None, None, None, None
            else:
                torch.distributed.all_reduce(
                    grad_input_expanded, op=torch.distributed.ReduceOp.SUM, group=group_replicate
                )
                torch.distributed.all_reduce(
                    grad_transposed_recv, op=torch.distributed.ReduceOp.SUM, group=group_replicate
                )
                grad_input_expanded = grad_input_expanded.squeeze(dim=axis + 1)
                grad_transposed_recv = grad_transposed_recv.squeeze(dim=axis + 1)
                grad_input = DTensor.from_local(
                    grad_input_expanded,
                    ctx.device_mesh,
                    ctx.input_placements,
                    shape=ctx.input_shape,
                    stride=ctx.input_stride,
                )
                grad_input_t = DTensor.from_local(
                    grad_transposed_recv,
                    ctx.device_mesh,
                    ctx.input_placements,
                    shape=ctx.input_t_shape,
                    stride=ctx.input_t_stride,
                )
                return grad_input, grad_input_t, None, None, None
        elif op == OuterOp.CDIST:
            # we multiply the d_output_d_diff_local by the upstream adjoint so that
            # the rest of the backward pass is the same as the other symmetric ops
            (d_output_d_diff_local,) = ctx.saved_tensors
            # grad_output_local is now d_loss_d_diff as in the upstream adjoint of OuterOp.SUBTRACT
            grad_output_local = grad_output_local.unsqueeze(-1) * d_output_d_diff_local

        # grad on right summand
        grad_transposed = grad_output_local.sum(dim=axis, keepdim=True).transpose(axis, axis + 1).contiguous()
        if op == OuterOp.SUBTRACT or op == OuterOp.CDIST:
            grad_transposed = -grad_transposed
        grad_transposed_recv = transpose_comm.enqueue_to_dispatch(grad_transposed)
        # grad on left summand, which always retain the positive sign
        grad_input_expanded = grad_output_local.sum(dim=axis + 1, keepdim=True)
        transpose_comm.wait_until_finished()

        # perform allreduce to get the row- and column-wise contributions
        if ctx.is_symmetric:
            grad_input_local = (grad_input_expanded + grad_transposed_recv).squeeze(dim=axis + 1)
            torch.distributed.all_reduce(grad_input_local, op=torch.distributed.ReduceOp.SUM, group=group_replicate)

            # Convert gradients back to DTensors
            grad_input = DTensor.from_local(
                grad_input_local, ctx.device_mesh, ctx.input_placements, shape=ctx.input_shape, stride=ctx.input_stride
            )
            return grad_input, None, None, None, None
        else:
            torch.distributed.all_reduce(grad_input_expanded, op=torch.distributed.ReduceOp.SUM, group=group_replicate)
            torch.distributed.all_reduce(grad_transposed_recv, op=torch.distributed.ReduceOp.SUM, group=group_replicate)
            grad_transposed_recv = grad_transposed_recv.squeeze(dim=axis + 1)
            grad_input_expanded = grad_input_expanded.squeeze(dim=axis + 1)

            # Convert gradients back to DTensors
            grad_input = DTensor.from_local(
                grad_input_expanded,
                ctx.device_mesh,
                ctx.input_placements,
                shape=ctx.input_shape,
                stride=ctx.input_stride,
            )
            grad_input_t = DTensor.from_local(
                grad_transposed_recv,
                ctx.device_mesh,
                ctx.input_placements,
                shape=ctx.input_t_shape,
                stride=ctx.input_t_stride,
            )
            return grad_input, grad_input_t, None, None, None


def distributed_outer_op(
    input: torch.Tensor,
    op: OuterOp,
    axis: int,
    input_t: Optional[torch.Tensor] = None,
    transpose_comm: Optional[TransposeComm] = None,
    group_replicate: Optional[torch.distributed.ProcessGroup] = None,
) -> torch.Tensor:
    """Perform an outer op operation with optional distribution across processes.

    This function computes the outer op of a tensor along a specified axis. When
    transpose_comm and group_replicate are provided, the operation is performed in
    a distributed manner across multiple processes.

    Args:
        input: Input tensor for outer op operation. This is assumed sharded
               along axis 0 of the transpose_comm's 2d group grid and also
               replicated along axis 1 of the same grid.
        op: the binary operation to perform
        axis: Axis along which to perform the outer op
        transpose_comm: Optional communication object for distributed operations
        group_replicate: Optional process group for replication across ranks
        input_t: Optional second input tensor for outer op operation. Will use `input` if not provided with None.

    Returns:
        Tensor with outer op computed

    Raises:
        ValueError: If only one of transpose_comm or group_replicate is provided
    """
    if (transpose_comm is None) != (group_replicate is None):
        raise ValueError("transpose_comm and group_replicate must both be provided or both be None")
    if transpose_comm is None and group_replicate is None:
        if input_t is None:
            input_t = input
        input_expanded = input.unsqueeze(axis + 1)
        input_expanded_t = input_t.unsqueeze(axis + 1)
        if op == OuterOp.SUM:
            return input_expanded + input_expanded_t.transpose(axis, axis + 1)
        elif op == OuterOp.SUBTRACT:
            return input_expanded - input_expanded_t.transpose(axis, axis + 1)
        elif op == OuterOp.EQUAL:
            return input_expanded == input_expanded_t.transpose(axis, axis + 1)
        elif op == OuterOp.BITAND:
            return input_expanded & input_expanded_t.transpose(axis, axis + 1)
        else:
            raise ValueError(f"Unsupported operation: {op}")
    else:
        return DistributedOuterOp.apply(input, input_t, op, axis, transpose_comm, group_replicate)


def replicate_to_shard_outer_op(
    input: DTensor,
    op: OuterOp,
    axis: int,
    transpose_comm: TransposeComm,
    input_t: Optional[DTensor] = None,
) -> DTensor:
    """Perform an outer op operation on DTensors with distributed computation and R2S transformation.

    This function computes the outer op of DTensors along a specified axis with
    distributed processing across multiple devices. It performs a Replicate->Shard (R2S)
    transformation where the input's Replicate() placement becomes Shard(2) in the output.

    Args:
        input: Input DTensor for outer op operation. Must have placements (Shard(0), Shard(1), Replicate())
        op: the binary operation to perform (SUM, SUBTRACT, EQUAL, BITAND, or CDIST)
        axis: Axis along which to perform the outer op. Must be 0 or 1 (corresponding to Shard dimensions).
        transpose_comm: communication handle for distributed operations
        input_t: Optional second input DTensor for outer op operation. Will use `input` if not provided.
                Note: For CDIST operation, input_t is ignored and input is used for both operands.

    Returns:
        DTensor with outer op computed and placements (Shard(0), Shard(1), Shard(2))
        For CDIST operation, the output shape is (B, N, N) instead of (B, N, N, D)

    Raises:
        TypeError: If inputs are not DTensors
        ValueError: If placements are incorrect or transpose_comm is None when DTensors are provided
    """
    return _ReplicateToShardOuterOp.apply(input, input_t, op, axis, transpose_comm)
