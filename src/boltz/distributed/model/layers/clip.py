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
from torch.distributed.tensor import DTensor, Partial, Shard


class _ClipImpl(torch.autograd.Function):
    """Distributed implementation of clipping operation using DTensors.

    This autograd function implements distributed clipping operations that constrain
    tensor values to be within specified bounds. The operation is performed element-wise
    across distributed tensors while maintaining proper gradient computation.

    Supported operations:
    - CLIP: output = torch.clip(tensor, min=min_val, max=max_val)

    Key features:
    - Distributed computation across device meshes with various sharding strategies
    - Memory-efficient implementation that operates on local tensor chunks
    - Supports gradient computation through custom backward pass
    - Supports both min and max clipping bounds (either can be None)
    """

    @staticmethod
    def forward(ctx, tensor: DTensor, min_val: Optional[float] = None, max_val: Optional[float] = None) -> DTensor:
        """Forward pass of distributed clipping operation.

        Parameters
        ----------
        ctx : torch.autograd.function.BackwardCFrame
            Context object for saving information needed in backward pass.
        tensor : DTensor
            Input tensor. Can have any shape and sharding strategy.
        min_val : Optional[float], default None
            Minimum value for clipping. If None, no lower bound is applied.
        max_val : Optional[float], default None
            Maximum value for clipping. If None, no upper bound is applied.

        Returns
        -------
        DTensor
            Output tensor with shape identical to input tensor.
            Contains the result of the clipping operation.

        Raises
        ------
        TypeError
            If inputs are not of expected types.
        ValueError
            If Partial placements are used (not supported), or if both min_val and max_val are None.
        """
        if not isinstance(tensor, DTensor):
            raise TypeError(f"Input 'tensor' must be of type DTensor. Got type {type(tensor)}.")

        if min_val is None and max_val is None:
            raise ValueError("At least one of min_val or max_val must be specified for clipping.")

        device_mesh_input = tensor.device_mesh
        placements_input = tensor.placements

        for i_dim_device_mesh, placement in enumerate(placements_input):
            if isinstance(placement, Partial):
                raise ValueError("Partial placements are not supported")
            elif isinstance(placement, Shard):
                if tensor.shape[placement.dim] % device_mesh_input.shape[i_dim_device_mesh] != 0:
                    raise ValueError(
                        f"Uneven sharding tensor dimension {placement.dim} of size {tensor.shape[placement.dim]} "
                        f"along device mesh dimension {i_dim_device_mesh} of size {device_mesh_input.shape[i_dim_device_mesh]} is not supported"
                    )

        tensor_local = tensor.to_local()

        # Perform the clipping operation
        output_local = torch.clip(tensor_local, min=min_val, max=max_val)

        if tensor.requires_grad:
            # Pre-allocate mask in bool to save memory from saving a float tensor_local copy
            mask_local = torch.ones_like(tensor_local, dtype=torch.bool)
            if min_val is not None:
                mask_local = mask_local & (tensor_local >= min_val)  # inclusive in torch.clip
            if max_val is not None:
                mask_local = mask_local & (tensor_local <= max_val)  # inclusive in torch.clip
            ctx.save_for_backward(mask_local)
            ctx.device_mesh_input = device_mesh_input
            ctx.placements_input = placements_input
            ctx.input_shape = tensor.shape
            ctx.min_val = min_val
            ctx.max_val = max_val

        out = DTensor.from_local(
            output_local,
            device_mesh=device_mesh_input,
            placements=placements_input,
            shape=tensor.shape,
            stride=tensor.stride(),
        )
        return out

    @staticmethod
    def backward(ctx, grad_output: DTensor) -> tuple[DTensor | None, None, None]:
        """Backward pass of distributed clipping operation.

        Computes gradients with respect to the input tensor.

        The gradient is:
        - For CLIP: d_tensor = grad_output * mask
          where mask = 1 for elements within bounds, 0 for clipped elements

        Parameters
        ----------
        ctx : torch.autograd.function.BackwardCFrame
            Context object containing saved tensors and metadata from forward pass.
        grad_outputs : tuple
            Gradients of the loss with respect to the output tensors.

        Returns
        -------
        tuple[DTensor | None, None, None]
            Gradients with respect to tensor, min_val, and max_val parameters.
            Only tensor gradient is computed; min_val and max_val gradients are None.
        """
        if not ctx.needs_input_grad[0]:
            return None, None, None

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

        grad_output_local = grad_output.to_local()
        (mask_local,) = ctx.saved_tensors

        # Compute gradient mask: 1 for elements within bounds, 0 for clipped elements
        d_tensor_local = grad_output_local * mask_local
        d_tensor = DTensor.from_local(
            d_tensor_local,
            device_mesh=ctx.device_mesh_input,
            placements=ctx.placements_input,
            shape=grad_output.shape,
            stride=grad_output.stride(),
        )

        return d_tensor, None, None


def clip(tensor: DTensor, min_val: Optional[float] = None, max_val: Optional[float] = None) -> DTensor:
    """Apply clipping operation to a distributed tensor.

    This function constrains tensor values to be within specified bounds.
    Elements below min_val are set to min_val, and elements above max_val are set to max_val.
    The operation is performed efficiently using local tensor operations while maintaining
    gradient computation capabilities.

    Parameters
    ----------
    tensor : DTensor
        Input tensor. Can have any shape and sharding strategy.
    min_val : Optional[float], default None
        Minimum value for clipping. If None, no lower bound is applied.
    max_val : Optional[float], default None
        Maximum value for clipping. If None, no upper bound is applied.

    Returns
    -------
    DTensor
        Output tensor with shape identical to input tensor.
        Contains the result of the clipping operation.

    Examples
    --------
    >>> # Assume we have distributed tensor x with shape (B, N, D)
    >>> clipped_positive = clip(x, min_val=0.0)
    >>> # clipped_positive = torch.clip(x, min=0.0), computed in distributed fashion
    >>>
    >>> clipped_range = clip(x, min_val=-1.0, max_val=1.0)
    >>> # clipped_range = torch.clip(x, min=-1.0, max=1.0), computed in distributed fashion
    >>>
    >>> clipped_max = clip(x, max_val=10.0)
    >>> # clipped_max = torch.clip(x, max=10.0), computed in distributed fashion

    Notes
    -----
    - Input tensor must be a DTensor with any placement strategy
    - Partial placements are not currently supported
    - The function is differentiable and supports gradient computation
    - The operation is performed on local tensor chunks for efficiency
    - At least one of min_val or max_val must be specified

    Raises
    ------
    TypeError
        If input tensor is not a DTensor.
    ValueError
        If Partial placements are used (not supported), or if both min_val and max_val are None.
    """
    return _ClipImpl.apply(tensor, min_val, max_val)  # type: ignore
