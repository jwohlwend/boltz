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

from boltz.distributed.utils import LayoutRightMap


class _WhereImpl(torch.autograd.Function):
    """Distributed implementation of where operation using DTensors.

    This autograd function implements distributed where operations that select
    elements from two tensors based on a condition. The operation is performed
    element-wise across distributed tensors while maintaining proper gradient computation.

    Supported operations:
    - WHERE: output = torch.where(condition, x, y)

    Key features:
    - Distributed computation across device meshes with various sharding strategies
    - Memory-efficient implementation that operates on local tensor chunks
    - Supports gradient computation through custom backward pass
    - Handles broadcasting between condition, x, and y tensors
    """

    @staticmethod
    def forward(ctx, condition: DTensor, x: DTensor, y: DTensor) -> DTensor:
        """Forward pass of distributed where operation.

        Parameters
        ----------
        ctx : torch.autograd.function.BackwardCFrame
            Context object for saving information needed in backward pass.
        condition : DTensor
            Boolean condition tensor. Must be broadcastable with x and y.
        x : DTensor
            Values to select where condition is True.
        y : DTensor
            Values to select where condition is False.

        Returns
        -------
        DTensor
            Output tensor with same shape as the broadcasted shape of inputs.
            Contains x where condition is True, y where condition is False.

        Raises
        ------
        TypeError
            If inputs are not DTensors.
        ValueError
            If Partial placements are used (not supported), or if tensors have
            incompatible device meshes or placements.
        """
        if not isinstance(condition, DTensor):
            raise TypeError(f"Input 'condition' must be of type DTensor. Got type {type(condition)}.")
        if not isinstance(x, DTensor):
            raise TypeError(f"Input 'x' must be of type DTensor. Got type {type(x)}.")
        if not isinstance(y, DTensor):
            raise TypeError(f"Input 'y' must be of type DTensor. Got type {type(y)}.")

        # Validate that all tensors have same device mesh and placements
        if condition.device_mesh != x.device_mesh or x.device_mesh != y.device_mesh:
            raise ValueError(
                f"All input tensors must have identical device mesh. "
                f"Got device meshes {condition.device_mesh}, {x.device_mesh}, {y.device_mesh}."
            )

        if condition.placements != x.placements or x.placements != y.placements:
            raise ValueError(
                f"All input tensors must have identical placements. "
                f"Got placements {condition.placements}, {x.placements}, {y.placements}."
            )

        device_mesh_input = x.device_mesh
        placements_input = x.placements

        # Validate placements
        for i_dim_device_mesh, placement in enumerate(placements_input):
            if isinstance(placement, Partial):
                raise ValueError("Partial placements are not supported")
            elif isinstance(placement, Shard):
                # Check that all tensors can be evenly sharded
                for tensor_name, tensor in [("condition", condition), ("x", x), ("y", y)]:
                    if tensor.shape[placement.dim] % device_mesh_input.shape[i_dim_device_mesh] != 0:
                        raise ValueError(
                            f"Uneven sharding {tensor_name} tensor dimension {placement.dim} of size {tensor.shape[placement.dim]} "
                            f"along device mesh dimension {i_dim_device_mesh} of size {device_mesh_input.shape[i_dim_device_mesh]} is not supported"
                        )

        # Get local tensors
        condition_local = condition.to_local()
        x_local = x.to_local()
        y_local = y.to_local()

        # Perform the where operation
        output_local = torch.where(condition_local, x_local, y_local)

        if x.requires_grad or y.requires_grad:
            # Save condition for backward pass
            condition_local_copy = condition_local.detach().clone()
            ctx.save_for_backward(condition_local_copy)
            ctx.device_mesh_input = device_mesh_input
            ctx.placements_input = placements_input
            ctx.x_requires_grad = x.requires_grad
            ctx.y_requires_grad = y.requires_grad
            ctx.shape_x = x.shape
            ctx.stride_x = x.stride()
            ctx.shape_y = y.shape
            ctx.stride_y = y.stride()

        # x and y shapes are only constrained to be broadcastable
        # without necessarily being the same shape
        shape_output = torch.broadcast_shapes(x.shape, y.shape)
        stride_output = LayoutRightMap(shape_output).strides
        out = DTensor.from_local(
            output_local,
            device_mesh=device_mesh_input,
            placements=placements_input,
            shape=shape_output,
            stride=stride_output,
        )
        return out

    @staticmethod
    def backward(ctx, grad_output: DTensor) -> tuple[None, DTensor | None, DTensor | None]:
        """Backward pass of distributed where operation.

        Computes gradients with respect to x and y inputs.

        The gradients are:
        - For x: grad_output where condition is True, 0 elsewhere
        - For y: grad_output where condition is False, 0 elsewhere
        - For condition: None (condition is not differentiable)

        Parameters
        ----------
        ctx : torch.autograd.function.BackwardCFrame
            Context object containing saved tensors and metadata from forward pass.
        grad_output : DTensor
            Gradients of the loss with respect to the output tensor.

        Returns
        -------
        tuple[None, DTensor | None, DTensor | None]
            Gradients with respect to x and y parameters.
            condition gradient is always None.
        """
        if not (ctx.x_requires_grad or ctx.y_requires_grad):
            return None, None, None

        if not isinstance(grad_output, DTensor):
            raise TypeError(f"Input 'grad_output' must be of type DTensor. Got type {type(grad_output)}.")

        if grad_output.device_mesh != ctx.device_mesh_input:
            raise ValueError(
                f"Input 'grad_output' must have the same device mesh as the input tensors. "
                f"Got device meshes {grad_output.device_mesh} and {ctx.device_mesh_input}."
            )

        if grad_output.placements != ctx.placements_input:
            raise ValueError(
                f"Input 'grad_output' must have the same placements as the input tensors. "
                f"Got placements {grad_output.placements} and {ctx.placements_input}."
            )

        grad_output_local = grad_output.to_local()
        (condition_local,) = ctx.saved_tensors

        # Compute gradients
        grad_x = None
        grad_y = None
        zeros_local = torch.zeros_like(grad_output_local)

        if ctx.x_requires_grad:
            # Gradient flows to x where condition is True
            grad_x_local = torch.where(condition_local, grad_output_local, zeros_local)
            grad_x = DTensor.from_local(
                grad_x_local,
                device_mesh=ctx.device_mesh_input,
                placements=ctx.placements_input,
                shape=ctx.shape_x,
                stride=ctx.stride_x,
            )

        if ctx.y_requires_grad:
            # Gradient flows to y where condition is False
            grad_y_local = torch.where(condition_local, zeros_local, grad_output_local)
            grad_y = DTensor.from_local(
                grad_y_local,
                device_mesh=ctx.device_mesh_input,
                placements=ctx.placements_input,
                shape=ctx.shape_y,
                stride=ctx.stride_y,
            )

        return None, grad_x, grad_y


def where(condition: DTensor, x: DTensor, y: DTensor) -> DTensor:
    """Apply where operation to distributed tensors.

    This function selects elements from x or y based on condition.
    Where condition is True, elements from x are selected; where condition is False,
    elements from y are selected. The operation is performed efficiently using local
    tensor operations while maintaining gradient computation capabilities.

    Parameters
    ----------
    condition : DTensor
        Boolean condition tensor. Must be broadcastable with x and y.
        Should have placements compatible with x and y.
    x : DTensor
        Values to select where condition is True.
        Can have any shape and sharding strategy compatible with condition and y.
    y : DTensor
        Values to select where condition is False.
        Must have same shape, device mesh, and placements as x.

    Returns
    -------
    DTensor
        Output tensor with same shape as the broadcasted shape of inputs.
        Contains x where condition is True, y where condition is False.

    Examples
    --------
    >>> # Assume we have distributed tensors with shape (B, N, D)
    >>> condition = x > 0.0
    >>> result = where(condition, x, y)
    >>> # result = torch.where(condition, x, y), computed in distributed fashion
    >>>
    >>> # Clip using where (equivalent to clip operation)
    >>> clipped = where(x > 5.0, torch.full_like(x, 5.0), x)
    >>> # clipped = torch.where(x > 5.0, 5.0, x), computed in distributed fashion

    Notes
    -----
    - All input tensors must be DTensors with compatible device meshes and placements
    - Partial placements are not currently supported
    - The function is differentiable and supports gradient computation for x and y
    - The condition tensor is not differentiable
    - The operation is performed on local tensor chunks for efficiency
    - Broadcasting is handled by PyTorch's local where operation

    Raises
    ------
    TypeError
        If inputs are not DTensors.
    ValueError
        If Partial placements are used (not supported), or if tensors have
        incompatible device meshes or placements.
    """
    return _WhereImpl.apply(condition, x, y)  # type: ignore
