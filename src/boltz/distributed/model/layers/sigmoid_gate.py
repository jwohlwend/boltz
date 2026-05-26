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


import torch
from torch.distributed.tensor import DTensor, Partial, Shard


class _SigmoidGateImpl(torch.autograd.Function):
    """Distributed implementation of sigmoid gating using DTensors.

    This autograd function implements a distributed sigmoid gating operation that applies
    a sigmoid-activated gate to an input tensor. The operation is performed element-wise
    across distributed tensors while maintaining proper gradient computation.

    The sigmoid gate computes:
        output = x * sigmoid(g)

    Key features:
    - Distributed computation across device meshes with various sharding strategies
    - Memory-efficient implementation that operates on local tensor chunks
    - Supports gradient computation through custom backward pass
    - Validates tensor compatibility (device mesh, placements, shapes)

    Notes
    -----
    Input tensors must be DTensors with:
    - Identical device mesh and placements
    - Compatible shapes (x and g must have the same shape)
    - No Partial placements (not currently supported)
    """

    @staticmethod
    def forward(ctx, x: DTensor, g: DTensor) -> DTensor:
        """Forward pass of distributed sigmoid gating.

        Parameters
        ----------
        ctx : torch.autograd.function.BackwardCFrame
            Context object for saving information needed in backward pass.
        x : DTensor
            Input tensor to be gated. Can have any shape and sharding strategy.
        g : DTensor
            Gate tensor with pre-sigmoid values. Must have identical shape,
            device mesh, and placements as x.

        Returns
        -------
        DTensor
            Output tensor with shape identical to input tensors.
            Contains the result of x * sigmoid(g).

        Raises
        ------
        TypeError
            If inputs are not DTensors.
        ValueError
            If tensors have incompatible device meshes, placements, or if
            Partial placements are used (not supported).
        """
        if not isinstance(x, DTensor):
            raise TypeError(f"Input 'x' must be of type DTensor. Got type {type(x)}.")
        if not isinstance(g, DTensor):
            raise TypeError(f"Input 'g' must be of type DTensor. Got type {type(g)}.")

        device_mesh_input = x.device_mesh
        if g.device_mesh != device_mesh_input:
            raise ValueError(
                f"Input tensors 'x' and 'g' must have identical device mesh. "
                f"Got device meshes {device_mesh_input} and {g.device_mesh}."
            )

        placements_input = x.placements
        for i_dim_device_mesh, placement in enumerate(placements_input):
            if isinstance(placement, Partial):
                raise ValueError("Partial placements are not supported")
            if isinstance(placement, Shard):
                if x.shape[placement.dim] % device_mesh_input.shape[i_dim_device_mesh] != 0:
                    raise ValueError(
                        f"Uneven sharding tensor dimension {placement.dim} of size {x.shape[placement.dim]} "
                        f"along device mesh dimension {i_dim_device_mesh} of size "
                        f"{device_mesh_input.shape[i_dim_device_mesh]} is not supported"
                    )

        if g.placements != placements_input:
            raise ValueError(
                f"Input tensors 'x' and 'g' must have identical placements. "
                f"Got placements {placements_input} and {g.placements}."
            )

        input_shape = x.shape
        if input_shape != g.shape:
            raise ValueError(
                f"Input tensors 'x' and 'g' must have identical shapes. Got shapes {input_shape} and {g.shape}."
            )

        g_local = g.to_local().sigmoid()
        x_gated_local = x.to_local() * g_local

        ctx.save_for_backward(x_gated_local, g_local)
        ctx.device_mesh_input = device_mesh_input
        ctx.placements_input = placements_input
        ctx.input_shape = input_shape

        out = DTensor.from_local(
            x_gated_local,
            device_mesh=device_mesh_input,
            placements=placements_input,
            shape=x.shape,
            stride=x.stride(),
        )
        return out

    @staticmethod
    def backward(ctx, grad_output: DTensor) -> tuple[DTensor, DTensor]:
        """Backward pass of distributed sigmoid gating.

        Computes gradients with respect to both input tensor x and gate tensor g.

        The gradients are:
        - dx = grad_output * sigmoid(g)
        - dg = grad_output * x * sigmoid(g) * (1 - sigmoid(g))

        Parameters
        ----------
        ctx : torch.autograd.function.BackwardCFrame
            Context object containing saved tensors and metadata from forward pass.
        grad_output : DTensor
            Gradient of the loss with respect to the output tensor.
            Must have identical device mesh and placements as the input tensors.

        Returns
        -------
        tuple[DTensor, DTensor]
            Gradients with respect to x and g respectively.
            Both have the same shape and distribution as their corresponding inputs.

        Raises
        ------
        TypeError
            If grad_output is not a DTensor.
        ValueError
            If grad_output has incompatible device mesh or placements compared
            to the input tensors from the forward pass.
        """
        if not isinstance(grad_output, DTensor):
            raise TypeError(f"Input 'grad_output' must be of type DTensor. Got type {type(grad_output)}.")

        if grad_output.device_mesh != ctx.device_mesh_input:
            raise ValueError(
                f"Input 'grad_output' must have the same device mesh as the input tensor. "
                f"Got device meshes {grad_output.device_mesh} and {ctx.device_mesh_input}."
            )

        if grad_output.placements != ctx.placements_input:
            raise ValueError(
                f"Input 'grad_output' must have the same placements as the input tensor. "
                f"Got placements {grad_output.placements} and {ctx.placements_input}."
            )

        if grad_output.shape != ctx.input_shape:
            raise ValueError(
                f"Input 'grad_output' must have the same shape as the input tensor. "
                f"Got shapes {grad_output.shape} and {ctx.input_shape}."
            )

        x_gated_local, g_local = ctx.saved_tensors
        grad_output_local = grad_output.to_local()

        dx_local = grad_output_local * g_local
        dx = DTensor.from_local(
            dx_local,
            device_mesh=ctx.device_mesh_input,
            placements=ctx.placements_input,
            shape=grad_output.shape,
            stride=grad_output.stride(),
        )

        dg_local = grad_output_local * x_gated_local
        dg_local *= 1 - g_local
        dg = DTensor.from_local(
            dg_local,
            device_mesh=ctx.device_mesh_input,
            placements=ctx.placements_input,
            shape=grad_output.shape,
            stride=grad_output.stride(),
        )

        return dx, dg


def sigmoid_gate(x: DTensor, g: DTensor) -> DTensor:
    """Apply sigmoid gating to a distributed tensor.

    This function performs element-wise sigmoid gating: x * sigmoid(g), where both
    input and gate tensors are distributed across multiple devices. The operation
    is performed efficiently using local tensor operations while maintaining
    gradient computation capabilities.

    Parameters
    ----------
    x : DTensor
        Input tensor to be gated. Can have any shape and sharding strategy.
    g : DTensor
        Gate tensor with pre-sigmoid values. Must have identical shape,
        device mesh, and placements as x.

    Returns
    -------
    DTensor
        Gated output tensor with shape identical to input tensors.
        Contains the result of x * sigmoid(g).

    Examples
    --------
    >>> # Assume we have distributed tensors x and g with shape (B, N, D)
    >>> output = sigmoid_gate(x, g)
    >>> # output = x * torch.sigmoid(g), computed in distributed fashion

    Notes
    -----
    - Both input tensors must be DTensors with compatible device meshes and placements
    - Partial placements are not currently supported
    - The function is differentiable and supports gradient computation
    - The operation is performed on local tensor chunks for efficiency
    """
    return _SigmoidGateImpl.apply(x, g)
