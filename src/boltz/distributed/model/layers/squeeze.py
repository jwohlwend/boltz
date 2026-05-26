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


class _ShardwiseUnsqueezeImpl(torch.autograd.Function):
    """Custom autograd function for performing unsqueeze operations on sharded distributed tensors.

    This function implements a differentiable unsqueeze operation that is compatible with
    PyTorch's distributed tensor (DTensor) framework. It handles the complexities of
    updating tensor placements and dimensions when adding a singleton dimension to
    a sharded tensor across multiple devices.

    The implementation ensures that:
    - Shard placements are correctly adjusted when dimensions shift due to unsqueeze
    - Partial placements are not supported (will raise an error)
    - Gradient computation is properly handled in the backward pass
    """

    @staticmethod
    def forward(ctx, x: DTensor, dim: int) -> DTensor:
        """Forward pass: performs unsqueeze operation on a distributed tensor.

        Args:
            ctx: PyTorch autograd context for saving information needed in backward pass
            x (DTensor): Input distributed tensor to unsqueeze
            dim (int): Dimension at which to insert the singleton dimension.
                      Can be negative (counted from the end)

        Returns:
            DTensor: Output tensor with an additional singleton dimension at the specified position

        Raises:
            TypeError: If x is not a DTensor or dim is not an int
            ValueError: If tensor has Partial placements (not supported) or if there's
                       uneven sharding that would prevent proper distribution

        Note:
            The function follows PyTorch's unsqueeze semantics for shape and stride computation.
            For sharded tensors, it updates the shard dimension indices when they are affected
            by the dimension insertion.
        """
        if not isinstance(x, DTensor):
            raise TypeError(f"Input 'x' must be of type DTensor. Got type {type(x)}.")
        if not isinstance(dim, int):
            raise TypeError(f"Input 'dim' must be of type int. Got type {type(dim)}.")

        device_mesh_input = x.device_mesh
        placements_input = x.placements
        shape_input = x.shape
        stride_input = x.stride()

        dim_to_insert = dim if dim >= 0 else x.ndim + 1 + dim
        placements_output = list(placements_input)
        for i_dim_device_mesh, p in enumerate(placements_input):
            if isinstance(p, Partial):
                raise ValueError("Partial placements are not supported")
            if isinstance(p, Shard):
                if x.shape[p.dim] % device_mesh_input.shape[i_dim_device_mesh] != 0:
                    raise ValueError(
                        f"Uneven sharding tensor dimension {p.dim} of size {x.shape[p.dim]} "
                        f"along device mesh dimension {i_dim_device_mesh} of size "
                        f"{device_mesh_input.shape[i_dim_device_mesh]} is not supported"
                    )
                if p.dim >= dim_to_insert:
                    placements_output[i_dim_device_mesh] = Shard(p.dim + 1)

        # update shape and stride according to pytorch unsqueeze() function:
        # https://github.com/pytorch/pytorch/blob/2c16eb9f3db0ba68520e5832d8bb6d3d875bdaeb/aten/src/ATen/native/TensorShape.cpp#L3879-L3890
        shape_output = list(shape_input)
        shape_output.insert(dim_to_insert, 1)
        stride_output = list(stride_input)
        stride_to_insert = 1 if dim_to_insert >= x.ndim else shape_input[dim_to_insert] * stride_input[dim_to_insert]
        stride_output.insert(dim_to_insert, stride_to_insert)

        # Perform unsqueeze on local tensor
        input_local = x.to_local()
        output_local = input_local.unsqueeze(dim_to_insert)

        # Save necessary information for backward pass
        ctx.device_mesh_input = device_mesh_input
        ctx.placements_input = placements_input
        ctx.placements_output = tuple(placements_output)
        ctx.dim_to_squeeze = dim_to_insert
        ctx.shape_input = shape_input
        ctx.stride_input = stride_input

        # Create output DTensor
        out = DTensor.from_local(
            output_local,
            shape=tuple(shape_output),
            stride=tuple(stride_output),
            device_mesh=device_mesh_input,
            placements=ctx.placements_output,
        )
        return out

    @staticmethod
    def backward(ctx, grad_output: DTensor) -> tuple[DTensor, None]:
        """Backward pass: computes gradients by performing squeeze operation.

        This method implements the reverse operation of unsqueeze for gradient computation.
        It takes the gradient with respect to the output and computes the gradient with
        respect to the input by squeezing the dimension that was added in the forward pass.

        Args:
            ctx: PyTorch autograd context containing saved information from forward pass
            grad_output (DTensor): Gradient with respect to the output tensor

        Returns:
            tuple[DTensor, None]: Tuple containing:
                - Gradient with respect to input tensor (DTensor or None if not needed)
                - None for the dim parameter (int parameters don't need gradients)

        Raises:
            TypeError: If grad_output is not a DTensor
            ValueError: If grad_output has incompatible device mesh or placements
                       compared to the original input tensor
        """
        if not isinstance(grad_output, DTensor):
            raise TypeError(f"Input 'grad_output' must be of type DTensor. Got type {type(grad_output)}.")

        if grad_output.device_mesh != ctx.device_mesh_input:
            raise ValueError(
                f"Input 'grad_output' must have the same device mesh as the input tensor. "
                f"Got device meshes {grad_output.device_mesh} and {ctx.device_mesh_input}."
            )

        if grad_output.placements != ctx.placements_output:
            raise ValueError(
                f"Input 'grad_output' must have the same placements as the input tensor. "
                f"Got placements {grad_output.placements} and {ctx.placements_output}."
            )

        if ctx.needs_input_grad[0]:
            # Perform squeeze on gradient
            grad_output_local = grad_output.to_local()
            grad_input_local = grad_output_local.squeeze(ctx.dim_to_squeeze)

            grad_input = DTensor.from_local(
                grad_input_local,
                shape=ctx.shape_input,
                stride=ctx.stride_input,
                device_mesh=ctx.device_mesh_input,
                placements=ctx.placements_input,
            )
        else:
            grad_input = None

        return grad_input, None


class _ShardwiseSqueezeImpl(torch.autograd.Function):
    """Custom autograd function for performing squeeze operations on sharded distributed tensors.

    This function implements a differentiable squeeze operation that is compatible with
    PyTorch's distributed tensor (DTensor) framework. It handles the complexities of
    updating tensor placements and dimensions when removing singleton dimensions from
    a sharded tensor across multiple devices.

    The implementation ensures that:
    - Shard placements are correctly adjusted when dimensions shift due to squeeze
    - Partial placements are not supported (will raise an error)
    - Gradient computation is properly handled in the backward pass
    - Only singleton dimensions (size 1) can be squeezed
    """

    @staticmethod
    def forward(ctx, x: DTensor, dim: int) -> DTensor:
        """Forward pass: performs squeeze operation on a distributed tensor.

        Args:
            ctx: PyTorch autograd context for saving information needed in backward pass
            x (DTensor): Input distributed tensor to squeeze
            dim (int): Dimension to squeeze (must be a singleton dimension).
                      Can be negative (counted from the end)

        Returns:
            DTensor: Output tensor with the singleton dimension removed

        Raises:
            TypeError: If x is not a DTensor or dim is not an int
            ValueError: If tensor has Partial placements (not supported), if there's
                       uneven sharding, or if the dimension to squeeze is not singleton
        """
        if not isinstance(x, DTensor):
            raise TypeError(f"Input 'x' must be of type DTensor. Got type {type(x)}.")
        if not isinstance(dim, int):
            raise TypeError(f"Input 'dim' must be of type int. Got type {type(dim)}.")

        device_mesh_input = x.device_mesh
        placements_input = x.placements
        shape_input = x.shape
        stride_input = x.stride()

        # Convert negative dim to positive
        dim_to_squeeze = dim if dim >= 0 else x.ndim + dim

        # Check if dimension is valid and singleton
        if dim_to_squeeze < 0 or dim_to_squeeze >= x.ndim:
            raise ValueError(f"Dimension {dim} is out of range for tensor with {x.ndim} dimensions")
        if shape_input[dim_to_squeeze] != 1:
            raise ValueError(f"Cannot squeeze dimension {dim_to_squeeze} with size {shape_input[dim_to_squeeze]}")

        placements_output = list(placements_input)
        for i_dim_device_mesh, p in enumerate(placements_input):
            if isinstance(p, Partial):
                raise ValueError("Partial placements are not supported")
            if isinstance(p, Shard):
                # Check if trying to squeeze a sharded dimension
                if p.dim == dim_to_squeeze:
                    raise ValueError(f"Cannot squeeze dimension {dim_to_squeeze} as it is sharded")
                if x.shape[p.dim] % device_mesh_input.shape[i_dim_device_mesh] != 0:
                    raise ValueError(
                        f"Uneven sharding tensor dimension {p.dim} of size {x.shape[p.dim]} "
                        f"along device mesh dimension {i_dim_device_mesh} of size "
                        f"{device_mesh_input.shape[i_dim_device_mesh]} is not supported"
                    )
                if p.dim > dim_to_squeeze:
                    placements_output[i_dim_device_mesh] = Shard(p.dim - 1)

        # Update shape and stride according to pytorch squeeze() function
        shape_output = list(shape_input)
        shape_output.pop(dim_to_squeeze)
        stride_output = list(stride_input)
        stride_output.pop(dim_to_squeeze)

        # Perform squeeze on local tensor
        input_local = x.to_local()
        output_local = input_local.squeeze(dim_to_squeeze)

        # Save necessary information for backward pass
        ctx.device_mesh_input = device_mesh_input
        ctx.placements_input = placements_input
        ctx.placements_output = tuple(placements_output)
        ctx.dim_to_unsqueeze = dim_to_squeeze
        ctx.shape_input = shape_input
        ctx.stride_input = stride_input

        # Create output DTensor
        out = DTensor.from_local(
            output_local,
            shape=tuple(shape_output),
            stride=tuple(stride_output),
            device_mesh=device_mesh_input,
            placements=ctx.placements_output,
        )
        return out

    @staticmethod
    def backward(ctx, grad_output: DTensor) -> tuple[DTensor, None]:
        """Backward pass: computes gradients by performing unsqueeze operation.

        This method implements the reverse operation of squeeze for gradient computation.
        It takes the gradient with respect to the output and computes the gradient with
        respect to the input by unsqueezing the dimension that was removed in the forward pass.

        Args:
            ctx: PyTorch autograd context containing saved information from forward pass
            grad_output (DTensor): Gradient with respect to the output tensor

        Returns:
            tuple[DTensor, None]: Tuple containing:
                - Gradient with respect to input tensor (DTensor or None if not needed)
                - None for the dim parameter (int parameters don't need gradients)

        Raises:
            TypeError: If grad_output is not a DTensor
            ValueError: If grad_output has incompatible device mesh or placements
                       compared to the original input tensor
        """
        if not isinstance(grad_output, DTensor):
            raise TypeError(f"Input 'grad_output' must be of type DTensor. Got type {type(grad_output)}.")

        if grad_output.device_mesh != ctx.device_mesh_input:
            raise ValueError(
                f"Input 'grad_output' must have the same device mesh as the input tensor. "
                f"Got device meshes {grad_output.device_mesh} and {ctx.device_mesh_input}."
            )

        if grad_output.placements != ctx.placements_output:
            raise ValueError(
                f"Input 'grad_output' must have the same placements as the input tensor. "
                f"Got placements {grad_output.placements} and {ctx.placements_output}."
            )

        if ctx.needs_input_grad[0]:
            # Perform unsqueeze on gradient
            grad_output_local = grad_output.to_local()
            grad_input_local = grad_output_local.unsqueeze(ctx.dim_to_unsqueeze)

            grad_input = DTensor.from_local(
                grad_input_local,
                shape=ctx.shape_input,
                stride=ctx.stride_input,
                device_mesh=ctx.device_mesh_input,
                placements=ctx.placements_input,
            )
        else:
            grad_input = None

        return grad_input, None


def shardwise_squeeze(x: DTensor, dim: int) -> DTensor:
    """Performs a squeeze operation on a sharded distributed tensor.

    This function removes a singleton dimension from a distributed tensor at the specified
    position while maintaining proper sharding across multiple devices. It's designed
    to work seamlessly with PyTorch's autograd system for gradient computation.

    Args:
        x (DTensor): Input distributed tensor to squeeze
        dim (int): Dimension to squeeze (must be a singleton dimension with size 1).
                  Can be negative (counted from the end). Valid range is
                  [-x.ndim, x.ndim-1] where negative values are converted
                  to positive using: x.ndim + dim

    Returns:
        DTensor: New distributed tensor with the singleton dimension removed.
                The output tensor will have one fewer dimension than the input,
                with all other dimensions unchanged.

    Raises:
        TypeError: If x is not a DTensor or dim is not an int
        ValueError: If the tensor has unsupported placement types, incompatible
                   sharding configurations, or if the dimension to squeeze is not singleton

    Examples:
        >>> # Assuming we have a 3D distributed tensor of shape (4, 1, 6)
        >>> x = ...  # DTensor with shape (4, 1, 6)
        >>> y = shardwise_squeeze(x, dim=1)  # Remove singleton dimension at position 1
        >>> print(y.shape)  # (4, 6)

        >>> # Using negative indexing
        >>> z = shardwise_squeeze(x, dim=-2)  # Remove singleton dimension at position 1 (from end)
        >>> print(z.shape)  # (4, 6)

    Note:
        This function is a wrapper around the _ShardwiseSqueezeImpl autograd function,
        making it more convenient to use while maintaining full gradient support.

        The function handles the complexity of updating shard placements when dimensions
        are shifted due to the squeeze operation, ensuring correct distributed behavior.

        Unlike PyTorch's squeeze() which can squeeze all singleton dimensions when no dim
        is specified, this function requires an explicit dimension to be specified.
    """
    return _ShardwiseSqueezeImpl.apply(x, dim)


def shardwise_unsqueeze(x: DTensor, dim: int) -> DTensor:
    """Performs an unsqueeze operation on a sharded distributed tensor.

    This function adds a singleton dimension to a distributed tensor at the specified
    position while maintaining proper sharding across multiple devices. It's designed
    to work seamlessly with PyTorch's autograd system for gradient computation.

    Args:
        x (DTensor): Input distributed tensor to unsqueeze
        dim (int): Dimension at which to insert the singleton dimension.
                  Can be negative (counted from the end). Valid range is
                  [-x.ndim-1, x.ndim] where negative values are converted
                  to positive using: x.ndim + 1 + dim

    Returns:
        DTensor: New distributed tensor with an additional singleton dimension.
                The output tensor will have one more dimension than the input,
                with all other dimensions unchanged.

    Raises:
        TypeError: If x is not a DTensor or dim is not an int
        ValueError: If the tensor has unsupported placement types or incompatible
                   sharding configurations

    Examples:
        >>> # Assuming we have a 2D distributed tensor of shape (4, 6)
        >>> x = ...  # DTensor with shape (4, 6)
        >>> y = shardwise_unsqueeze(x, dim=1)  # Insert at dimension 1
        >>> print(y.shape)  # (4, 1, 6)

        >>> # Using negative indexing
        >>> z = shardwise_unsqueeze(x, dim=-1)  # Insert at last position
        >>> print(z.shape)  # (4, 6, 1)

    Note:
        This function is a wrapper around the _ShardwiseUnsqueezeImpl autograd function,
        making it more convenient to use while maintaining full gradient support.

        The function handles the complexity of updating shard placements when dimensions
        are shifted due to the unsqueeze operation, ensuring correct distributed behavior.
    """
    return _ShardwiseUnsqueezeImpl.apply(x, dim)
