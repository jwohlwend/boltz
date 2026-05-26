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


from typing import Optional

import torch

from boltz.distributed.comm import TransposeComm


def transpose_then_redistribute(
    input: torch.Tensor, dim0: int, dim1: int, transpose_comm: TransposeComm
) -> torch.Tensor:
    """Transpose a tensor and redistribute it across processes.

    This function first performs a transpose operation on the input tensor
    and then redistributes the result using the provided communication object.

    Parameters
    ----------
    input : torch.Tensor
        Input tensor to transpose
    dim0 : int
        First dimension to transpose
    dim1 : int
        Second dimension to transpose
    transpose_comm: TransposeComm
        Communication object for distributed operations

    Returns
    -------
    torch.Tensor
        Transposed and redistributed tensor

    """
    inputT = input.transpose(dim0, dim1).contiguous()
    inputT_recv = transpose_comm.enqueue_to_dispatch(inputT)
    transpose_comm.wait_until_finished()
    return inputT_recv


class RedistributeTranspose(torch.autograd.Function):
    """Custom autograd function to perform transpose with redistribution

    This operation performs a tensor transpose across a grid of processes
    encapsulated by the input TransposeComm object.
    """

    @staticmethod
    def forward(ctx, input: torch.Tensor, dim0: int, dim1: int, transpose_comm: TransposeComm) -> torch.Tensor:
        """Forward pass for RedistributeTranspose.

        Args:
            ctx: Context object to save information for backward pass
            input: Input tensor to transpose and redistribute
            dim0: First dimension to transpose
            dim1: Second dimension to transpose
            transpose_comm: Communication object for distributed operations

        Returns:
            Transposed and redistributed tensor
        """
        ctx.dim0 = dim0
        ctx.dim1 = dim1
        ctx.transpose_comm = transpose_comm
        return transpose_then_redistribute(input, dim0, dim1, transpose_comm)

    @staticmethod
    def backward(ctx, grad_output: torch.Tensor) -> tuple[torch.Tensor, None, None, None]:
        """Backward pass for RedistributeTranspose.

        Args:
            ctx: Context object with saved information from forward pass
            grad_output: Gradient tensor from downstream layers

        Returns:
            Tuple containing the gradient for input tensor and None for other parameters
        """
        dim0 = ctx.dim0
        dim1 = ctx.dim1
        transpose_comm = ctx.transpose_comm
        grad_input = transpose_then_redistribute(grad_output, dim1, dim0, transpose_comm)
        return grad_input, None, None, None


def redistribute_transpose(
    input: torch.Tensor, dim0: int, dim1: int, transpose_comm: Optional[TransposeComm] = None
) -> torch.Tensor:
    """Transpose a tensor with optional redistribution for distributed training.

    When the input TransposeComm is not None, the input tensor is redistributed
    across the grid of processes encapsulated by the TransposeComm object. This implies
    the return tensor's memory contiguity (layout right be default). By design, this
    intention is to be consistent with the equivalent operation:
    1) inputT_global = input.transpose(dim0, dim1)
    2) scatter(result, [inputT_global_chunk_0.contiguous(), inputT_global_chunk_1.contiguous(), ...],
               src=0) # scatter from root process (0) to all other processes
    where the the scatter op requires contiguous source tensors. Similarly for the the
    backward pass.

    Args:
        input: Input tensor to transpose
        dim0: First dimension to transpose
        dim1: Second dimension to transpose
        transpose_comm: Optional communication object for distributed operations.
                      If None, performs a regular transpose without redistribution.

    Returns:
        Transposed tensor, potentially redistributed across processes
    """
    if transpose_comm is None:
        return input.transpose(dim0, dim1)
    else:
        return RedistributeTranspose.apply(input, dim0, dim1, transpose_comm)
