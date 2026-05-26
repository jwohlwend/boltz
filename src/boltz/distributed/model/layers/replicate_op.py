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

"""Replicate op: lhs op rhs.unsqueeze(dim) with lhs sharded and rhs replicated on that dim."""

from enum import Enum, auto

import torch
from torch.autograd.function import FunctionCtx
from torch.distributed.tensor import DTensor, Partial, Replicate, Shard


class ReplicateOp(Enum):
    """Supported operations for replicate_op."""

    ADD = auto()
    SUB = auto()
    PROD = auto()
    DIV = auto()


class _ReplicateOpImpl(torch.autograd.Function):
    @staticmethod
    def forward(
        ctx: FunctionCtx,
        lhs: DTensor,
        rhs: DTensor,
        dim_to_unsqueeze_rhs: int,
        op: ReplicateOp,
    ) -> DTensor:
        if not isinstance(lhs, DTensor):
            raise TypeError(f"Input 'lhs' must be of type DTensor. Got type {type(lhs)}.")
        if not isinstance(rhs, DTensor):
            raise TypeError(f"Input 'rhs' must be of type DTensor. Got type {type(rhs)}.")
        if not isinstance(dim_to_unsqueeze_rhs, int):
            raise TypeError(f"Input 'dim_to_unsqueeze_rhs' must be of type int. Got type {type(dim_to_unsqueeze_rhs)}.")
        if not isinstance(op, ReplicateOp):
            raise TypeError(f"Input 'op' must be of type ReplicateOp. Got type {type(op)}.")

        if op not in (ReplicateOp.ADD, ReplicateOp.SUB, ReplicateOp.PROD, ReplicateOp.DIV):
            raise ValueError(f"Unsupported operation: {op}. Only ADD, SUB, PROD, and DIV are supported.")

        dim_to_unsqueeze_rhs = (
            dim_to_unsqueeze_rhs if dim_to_unsqueeze_rhs >= 0 else dim_to_unsqueeze_rhs + rhs.ndim + 1
        )
        shape_lhs_expected = (
            rhs.shape[:dim_to_unsqueeze_rhs] + (lhs.shape[dim_to_unsqueeze_rhs],) + rhs.shape[dim_to_unsqueeze_rhs:]
        )

        if lhs.shape != shape_lhs_expected:
            raise ValueError(
                f"Shape mismatch: lhs is expected to have the same shape as rhs except for "
                f"the unsqueezed dimension {shape_lhs_expected} "
                f"but got {lhs.shape}"
            )

        if lhs.device_mesh != rhs.device_mesh:
            raise ValueError(
                f"Device mesh mismatch: lhs.device_mesh={lhs.device_mesh} != rhs.device_mesh={rhs.device_mesh}"
            )

        placements_pair_expected = (Shard(dim_to_unsqueeze_rhs), Replicate())
        dim_device_mesh_reduce = None
        for i_dim_device_mesh, (p_lhs, p_rhs) in enumerate(zip(lhs.placements, rhs.placements)):
            if p_lhs == placements_pair_expected[0] and p_rhs == placements_pair_expected[1]:
                if dim_device_mesh_reduce is not None:
                    raise ValueError(
                        f"Duplicate placements pair {placements_pair_expected} found "
                        f"in lhs.placements {lhs.placements} and rhs.placements {rhs.placements}"
                    )
                dim_device_mesh_reduce = i_dim_device_mesh
                continue
            if isinstance(p_lhs, Partial) or isinstance(p_rhs, Partial):
                raise ValueError("Partial placements are not supported")
            if isinstance(p_lhs, Shard) and isinstance(p_rhs, Shard):
                if p_lhs.dim < dim_to_unsqueeze_rhs:
                    p_rhs_expected = Shard(p_lhs.dim)
                else:
                    p_rhs_expected = Shard(p_lhs.dim - 1)
                if p_rhs != p_rhs_expected:
                    raise ValueError(
                        f"rhs.placements[{i_dim_device_mesh}] is expected to be {p_rhs_expected} but got {p_rhs}"
                    )
            elif p_lhs != p_rhs:
                raise ValueError(
                    f"lhs.placements[{i_dim_device_mesh}] is expected to be the same as rhs.placements[{i_dim_device_mesh}] "
                    f"but got {p_lhs} and {p_rhs}"
                )
        if dim_device_mesh_reduce is None:
            raise ValueError(
                f"lhs.placements is expected to contain Shard({dim_to_unsqueeze_rhs}) and "
                f"rhs.placements is expected to contain Replicate() along the same device mesh axis "
                f"but got {lhs.placements} and {rhs.placements}"
            )

        if lhs.requires_grad or rhs.requires_grad:
            ctx.placements_output = lhs.placements
            ctx.placements_dLHS = lhs.placements
            ctx.placements_dRHS = rhs.placements
            ctx.device_mesh = lhs.device_mesh
            ctx.dim_to_squeeze_dRHS = dim_to_unsqueeze_rhs
            ctx.group_all_reduce_dRHS = ctx.device_mesh.get_group(dim_device_mesh_reduce)
            ctx.op = op
            ctx.shape_dLHS = lhs.shape
            ctx.stride_dLHS = lhs.stride()
            ctx.shape_dRHS = rhs.shape
            ctx.stride_dRHS = rhs.stride()

        lhs_local = lhs.to_local()
        rhs_local = rhs.to_local()
        rhs_unsqueezed_local = rhs_local.unsqueeze(dim_to_unsqueeze_rhs)

        if op == ReplicateOp.ADD:
            output_local = lhs_local + rhs_unsqueezed_local
        elif op == ReplicateOp.SUB:
            output_local = lhs_local - rhs_unsqueezed_local
        elif op == ReplicateOp.PROD:
            output_local = lhs_local * rhs_unsqueezed_local
            ctx.save_for_backward(
                lhs_local.detach().clone() if rhs.requires_grad else None,
                rhs_local.detach().clone() if lhs.requires_grad else None,
            )
        elif op == ReplicateOp.DIV:
            output_local = lhs_local / rhs_unsqueezed_local
            ctx.save_for_backward(
                lhs_local.detach().clone() if rhs.requires_grad else None,
                rhs_local.detach().clone(),
            )
        else:
            raise ValueError(f"Unsupported operation: {op}")

        output = DTensor.from_local(
            output_local, placements=lhs.placements, device_mesh=lhs.device_mesh, shape=lhs.shape, stride=lhs.stride()
        )
        return output

    @staticmethod
    def backward(ctx: FunctionCtx, grad_output: DTensor) -> tuple:
        grad_lhs = None
        grad_rhs = None

        grad_output_local = grad_output.to_local()

        if ctx.op == ReplicateOp.ADD:
            if ctx.needs_input_grad[0]:
                grad_lhs = DTensor.from_local(
                    grad_output_local.clone(),
                    placements=ctx.placements_dLHS,
                    device_mesh=ctx.device_mesh,
                    shape=ctx.shape_dLHS,
                    stride=ctx.stride_dLHS,
                )
            if ctx.needs_input_grad[1]:
                grad_rhs_reduced_local = grad_output_local.sum(dim=ctx.dim_to_squeeze_dRHS)
        elif ctx.op == ReplicateOp.SUB:
            if ctx.needs_input_grad[0]:
                grad_lhs = DTensor.from_local(
                    grad_output_local.clone(),
                    placements=ctx.placements_dLHS,
                    device_mesh=ctx.device_mesh,
                    shape=ctx.shape_dLHS,
                    stride=ctx.stride_dLHS,
                )
            if ctx.needs_input_grad[1]:
                grad_rhs_reduced_local = -grad_output_local.sum(dim=ctx.dim_to_squeeze_dRHS)
        elif ctx.op == ReplicateOp.PROD:
            lhs_local, rhs_local = ctx.saved_tensors
            if ctx.needs_input_grad[0]:
                rhs_unsqueezed_local = rhs_local.unsqueeze(ctx.dim_to_squeeze_dRHS)
                grad_lhs = DTensor.from_local(
                    grad_output_local * rhs_unsqueezed_local,
                    placements=ctx.placements_dLHS,
                    device_mesh=ctx.device_mesh,
                    shape=ctx.shape_dLHS,
                    stride=ctx.stride_dLHS,
                )
            if ctx.needs_input_grad[1]:
                grad_rhs_reduced_local = (grad_output_local * lhs_local).sum(dim=ctx.dim_to_squeeze_dRHS)
        elif ctx.op == ReplicateOp.DIV:
            lhs_local, rhs_local = ctx.saved_tensors
            rhs_unsqueezed_local = rhs_local.unsqueeze(ctx.dim_to_squeeze_dRHS)
            grad_lhs_local = grad_output_local / rhs_unsqueezed_local
            if ctx.needs_input_grad[0]:
                grad_lhs = DTensor.from_local(
                    grad_lhs_local,
                    placements=ctx.placements_dLHS,
                    device_mesh=ctx.device_mesh,
                    shape=ctx.shape_dLHS,
                    stride=ctx.stride_dLHS,
                )
            if ctx.needs_input_grad[1]:
                grad_rhs_reduced_local = -(grad_lhs_local * lhs_local / rhs_unsqueezed_local).sum(
                    dim=ctx.dim_to_squeeze_dRHS
                )

        if ctx.needs_input_grad[1]:
            torch.distributed.all_reduce(
                grad_rhs_reduced_local,
                op=torch.distributed.ReduceOp.SUM,
                group=ctx.group_all_reduce_dRHS,
                async_op=False,
            )
            grad_rhs = DTensor.from_local(
                grad_rhs_reduced_local,
                placements=ctx.placements_dRHS,
                device_mesh=ctx.device_mesh,
                shape=ctx.shape_dRHS,
                stride=ctx.stride_dRHS,
            )

        return grad_lhs, grad_rhs, None, None


def replicate_op(lhs: DTensor, rhs: DTensor, dim_to_unsqueeze_rhs: int, op: ReplicateOp) -> DTensor:
    """lhs op rhs.unsqueeze(dim_to_unsqueeze_rhs); lhs sharded on that dim, rhs replicated on it."""
    return _ReplicateOpImpl.apply(lhs, rhs, dim_to_unsqueeze_rhs, op)
