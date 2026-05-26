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


import math

import torch
from torch import Tensor
from torch.autograd.function import FunctionCtx
from torch.distributed.tensor import DTensor, Partial, Placement, Shard

from boltz.distributed.utils import LayoutRightMap, update_exhaustive_strides


def _shardwise_flatten(
    x: DTensor,
    start_dim: int = 0,
    end_dim: int = -1,
    output_placements: tuple | None = None,
    input_placements_expected: tuple | None = None,
) -> DTensor:
    """Generalized shardwise flattening operation with validation.

    This function performs input validation and flattening operation.

    Parameters
    ----------
    x : DTensor
        Input DTensor to flatten.
    start_dim : int, optional
        First dimension to flatten. Default is 0.
    end_dim : int, optional
        Last dimension to flatten. Default is -1 (last dimension).
    output_placements : tuple | None, optional
        If provided, use these placements for the output DTensor instead of
        computing them from the input placements. Default is None.
    input_placements_expected : tuple | None, optional
        If provided, skip validation by assuming input has these expected placements.
        Must be present if output_placements is present, and absent if output_placements
        is absent. Default is None.

    Returns
    -------
    DTensor
        Flattened DTensor.

    Raises
    ------
    ValueError
        Checks on the input x and parameters, or if placement argument constraints are violated.
    NotImplementedError
        If any dimension to be flattened is sharded.
    """
    # Validate placement argument constraints
    if (output_placements is None) != (input_placements_expected is None):
        raise ValueError("input_placements_expected must be present if and only if output_placements is present")

    has_input_output_placements = input_placements_expected is not None and output_placements is not None

    # Normalize dimensions
    start_dim_normalized = start_dim if start_dim >= 0 else start_dim + x.ndim
    end_dim_normalized = end_dim if end_dim >= 0 else end_dim + x.ndim

    # Validate dimension ranges
    if start_dim_normalized < 0 or start_dim_normalized >= x.ndim:
        raise ValueError(f"start_dim {start_dim} is out of range for tensor with {x.ndim} dimensions")
    if end_dim_normalized < 0 or end_dim_normalized >= x.ndim:
        raise ValueError(f"end_dim {end_dim} is out of range for tensor with {x.ndim} dimensions")
    if start_dim_normalized > end_dim_normalized:
        raise ValueError(f"start_dim {start_dim} must be <= end_dim {end_dim}")

    # Check that no dimension to be flattened is sharded
    # and there are no unevenly sharded tensor axes
    # Also calculate new placements accounting for dimension changes (if not provided)
    if has_input_output_placements:
        # Check if input placements match expected placements
        if x.placements != input_placements_expected:
            raise ValueError(
                f"Input placements {x.placements} do not match expected placements {input_placements_expected}"
            )
        # Skip validation and use provided placements directly
        new_placements = output_placements
    else:
        dims_removed = end_dim_normalized - start_dim_normalized
        new_placements = []

        for i_dim_device_mesh, placement in enumerate(x.placements):
            if isinstance(placement, Shard):
                i_dim_tensor = placement.dim
                if start_dim_normalized <= i_dim_tensor <= end_dim_normalized:
                    raise NotImplementedError(
                        f"Flattening dimension {i_dim_tensor} sharded by device_mesh axis {i_dim_device_mesh} is not supported"
                    )
                if x.shape[i_dim_tensor] % x.device_mesh.shape[i_dim_device_mesh] != 0:
                    raise ValueError(
                        f"Uneven sharding tensor dimension {i_dim_tensor} of size {x.shape[i_dim_tensor]} "
                        f"along device mesh dimension {i_dim_device_mesh} of size {x.device_mesh.shape[i_dim_device_mesh]} is not supported"
                    )

                # Calculate new placement for this shard
                if i_dim_tensor < start_dim_normalized:
                    # Dimension before flattened region - unchanged
                    new_placements.append(Shard(i_dim_tensor))
                elif i_dim_tensor > end_dim_normalized:
                    # Dimension after flattened region - shift left by dims_removed
                    new_placements.append(Shard(i_dim_tensor - dims_removed))
                # Dimensions within flattened region are handled by the validation above
            elif isinstance(placement, Partial):
                raise ValueError(f"Placements of type {Partial} are not supported")
            else:
                new_placements.append(placement)

    # Perform operation on local tensors
    x_local = x.to_local()
    output_local: Tensor = torch.flatten(x_local, start_dim=start_dim, end_dim=end_dim)

    # Compute output shape and stride
    flattened_size = math.prod(x.shape[start_dim_normalized : (end_dim_normalized + 1)])
    shape_output = x.shape[:start_dim_normalized] + (flattened_size,) + x.shape[end_dim_normalized + 1 :]
    # Use update_exhaustive_strides to compute new strides
    strides_output = update_exhaustive_strides(output_local.shape, output_local.stride(), shape_output)

    # Create output DTensor using input tensor's device mesh and updated placements
    output: DTensor = DTensor.from_local(
        output_local,
        device_mesh=x.device_mesh,
        placements=tuple(new_placements),
        shape=shape_output,
        stride=strides_output,
    )

    return output


def _shardwise_unflatten(
    x: DTensor,
    dim: int,
    sizes: tuple[int, ...],
    output_placements: tuple | None = None,
    input_placements_expected: tuple | None = None,
) -> DTensor:
    """Generalized shardwise unflattening operation with validation.

    This function performs input validation and unflattening operation following
    the torch.unflatten API.

    Parameters
    ----------
    x : DTensor
        Input flattened DTensor to unflatten.
    dim : int
        Dimension to unflatten.
    sizes : tuple[int, ...]
        Sizes for the new dimensions that will replace the specified dimension.
    output_placements : tuple | None, optional
        If provided, use these placements for the output DTensor instead of
        computing them from the input placements. Default is None.
    input_placements_expected : tuple | None, optional
        If provided, skip validation by assuming input has these expected placements.
        Must be present if output_placements is present, and absent if output_placements
        is absent. Default is None.

    Returns
    -------
    DTensor
        Unflattened DTensor with expanded dimensions.

    Raises
    ------
    ValueError
        Checks on the input x and parameters, or if placement argument constraints are violated.
    """
    # Validate placement argument constraints
    if (output_placements is None) != (input_placements_expected is None):
        raise ValueError("input_placements_expected must be present if and only if output_placements is present")

    has_input_output_placements = input_placements_expected is not None and output_placements is not None

    # Normalize dimension
    dim_normalized = dim if dim >= 0 else dim + x.ndim

    # Validate dimension range
    if dim_normalized < 0 or dim_normalized >= x.ndim:
        raise ValueError(f"dim {dim} is out of range for tensor with {x.ndim} dimensions")

    # Also calculate new placements accounting for dimension changes (if not provided)
    if has_input_output_placements:
        # Check if input placements match expected placements
        if x.placements != input_placements_expected:
            raise ValueError(
                f"Input placements {x.placements} do not match expected placements {input_placements_expected}"
            )
        # Skip validation and use provided placements directly
        new_placements = output_placements
    else:
        dims_added = len(sizes) - 1
        new_placements = []

        for i_dim_device_mesh, placement in enumerate(x.placements):
            # Check that the dimension to unflatten is not sharded
            # and there are no unevenly sharded tensor axes
            if isinstance(placement, Shard):
                i_dim_tensor = placement.dim
                if x.shape[i_dim_tensor] % x.device_mesh.shape[i_dim_device_mesh] != 0:
                    raise ValueError(
                        f"Uneven sharding tensor dimension {i_dim_tensor} of size {x.shape[i_dim_tensor]} "
                        f"along device mesh dimension {i_dim_device_mesh} of size {x.device_mesh.shape[i_dim_device_mesh]} is not supported"
                    )
                if i_dim_tensor == dim_normalized:
                    raise NotImplementedError(
                        f"Unflattening dimension {dim} shared by device_mesh axis {i_dim_device_mesh} is not supported"
                    )
                elif i_dim_tensor < dim_normalized:
                    # Dimension before unflattened region - unchanged
                    new_placements.append(Shard(i_dim_tensor))
                else:
                    # Dimension after unflattened region - shift right by dims_added
                    new_placements.append(Shard(i_dim_tensor + dims_added))
            elif isinstance(placement, Partial):
                raise ValueError(f"Placements of type {Partial} are not supported")
            else:
                new_placements.append(placement)

    # Perform operation on local tensors
    x_local = x.to_local()
    output_local: Tensor = torch.unflatten(x_local, dim=dim, sizes=sizes)

    # compute the output DTensor's shape and stride
    shape_output = list(x.shape)
    dim_to_insert = dim_normalized
    shape_output.pop(dim_to_insert)
    shape_output[dim_to_insert:dim_to_insert] = sizes
    # torch unflatten enforce contiguous layout, which is LayoutRight
    layout_right = LayoutRightMap(shape_output)
    strides_output = layout_right.strides

    # Create output DTensor using input tensor's device mesh and updated placements
    output: DTensor = DTensor.from_local(
        output_local,
        device_mesh=x.device_mesh,
        placements=tuple(new_placements),
        shape=tuple(shape_output),
        stride=strides_output,
    )

    return output


class _ShardWiseFlattenImpl(torch.autograd.Function):
    @staticmethod
    def forward(
        ctx: FunctionCtx,
        x: DTensor,
        start_dim: int = 0,
        end_dim: int = -1,
    ) -> DTensor:
        """Forward pass.

        Parameters
        ----------
        ctx : FunctionCtx
            Context object for saving information needed in backward pass.
        x : DTensor
            Input DTensor to flatten.
        start_dim : int, optional
            First dimension to flatten. Default is 0.
        end_dim : int, optional
            Last dimension to flatten. Default is -1 (last dimension).

        Returns
        -------
        DTensor
            Flattened DTensor.
        """
        # Normalize dimensions before flattening
        start_dim_normalized = start_dim if start_dim >= 0 else start_dim + x.ndim
        end_dim_normalized = end_dim if end_dim >= 0 else end_dim + x.ndim

        # Perform flattening with built-in validation
        result = _shardwise_flatten(x, start_dim, end_dim)

        # Save metadata for backward pass
        # For unflattening, we need the dimension in the flattened tensor and the original sizes
        ctx.unflatten_dim = start_dim_normalized  # This is the flattened dimension position
        ctx.unflatten_sizes = x.shape[start_dim_normalized : end_dim_normalized + 1]  # Original sizes
        ctx.device_mesh_input = x.device_mesh
        ctx.placements_input = x.placements
        ctx.placements_output = result.placements

        return result

    @staticmethod
    def backward(
        ctx: FunctionCtx,
        grad_output: DTensor,
    ) -> tuple[DTensor, None, None]:
        """Backward pass.

        Parameters
        ----------
        ctx : FunctionCtx
            Context object containing saved tensors and metadata from forward pass.
        grad_output : DTensor
            Gradient of the loss with respect to the output.

        Returns
        -------
        tuple[DTensor, None, None]
            Gradient with respect to input, None for start_dim and end_dim parameters.
        """
        # Use unflatten operation (inverse of flatten) for backward pass with built-in validation
        grad_x = _shardwise_unflatten(
            grad_output,
            ctx.unflatten_dim,
            ctx.unflatten_sizes,
            output_placements=ctx.placements_input,
            input_placements_expected=ctx.placements_output,
        )

        return grad_x, None, None


class _ShardWiseUnflattenImpl(torch.autograd.Function):
    @staticmethod
    def forward(
        ctx: FunctionCtx,
        x: DTensor,
        dim: int,
        sizes: tuple[int, ...],
    ) -> DTensor:
        """Forward pass.

        Parameters
        ----------
        ctx : FunctionCtx
            Context object for saving information needed in backward pass.
        x : DTensor
            Input DTensor to unflatten.
        dim : int
            Dimension to unflatten.
        sizes : tuple[int, ...]
            Sizes for the new dimensions that will replace the specified dimension.

        Returns
        -------
        DTensor
            Unflattened DTensor.
        """
        # Normalize dimension before unflattening
        dim_normalized = dim if dim >= 0 else dim + x.ndim

        # Perform unflattening with built-in validation
        result = _shardwise_unflatten(x, dim, sizes)

        # Save metadata for backward pass
        # For flattening, we need the start and end dimensions in the unflattened tensor
        ctx.flatten_start_dim = dim_normalized  # This is the start dimension for flattening
        ctx.flatten_end_dim = dim_normalized + len(sizes) - 1  # This is the end dimension for flattening
        ctx.device_mesh_input = x.device_mesh
        ctx.placements_input = x.placements
        ctx.placements_output = result.placements

        return result

    @staticmethod
    def backward(
        ctx: FunctionCtx,
        grad_output: DTensor,
    ) -> tuple[DTensor, None, None]:
        """Backward pass.

        Parameters
        ----------
        ctx : FunctionCtx
            Context object containing saved tensors and metadata from forward pass.
        grad_output : DTensor
            Gradient of the loss with respect to the output.

        Returns
        -------
        tuple[DTensor, None, None]
            Gradient with respect to input, None for dim and sizes parameters.
        """
        # Use flatten operation (inverse of unflatten) for backward pass with built-in validation
        grad_x = _shardwise_flatten(
            grad_output,
            ctx.flatten_start_dim,
            ctx.flatten_end_dim,
            output_placements=ctx.placements_input,
            input_placements_expected=ctx.placements_output,
        )

        return grad_x, None, None


def shardwise_flatten(x: DTensor, start_dim: int = 0, end_dim: int = -1) -> DTensor:
    """Flatten a DTensor along specified dimensions.

    This function flattens a DTensor from start_dim to end_dim (inclusive).
    The dimensions to be flattened must not be sharded on the device mesh.

    Parameters
    ----------
    x : DTensor
        Input DTensor to flatten.
    start_dim : int, optional
        First dimension to flatten. Default is 0.
    end_dim : int, optional
        Last dimension to flatten. Default is -1 (last dimension).

    Returns
    -------
    DTensor
        Flattened DTensor.

    Raises
    ------
    ValueError
        If any specified dimension is sharded or other validation errors.
    NotImplementedError
        If any dimension to be flattened is sharded.
    """
    return _ShardWiseFlattenImpl.apply(x, start_dim, end_dim)


def shardwise_unflatten(x: DTensor, dim: int, sizes: tuple[int, ...]) -> DTensor:
    """Unflatten a DTensor along a specified dimension.

    This function unflattens a DTensor by expanding the specified dimension
    into multiple dimensions with the given sizes. The dimension to be
    unflattened must not be sharded on the device mesh.

    Parameters
    ----------
    x : DTensor
        Input DTensor to unflatten.
    dim : int
        Dimension to unflatten.
    sizes : tuple[int, ...]
        Sizes for the new dimensions that will replace the specified dimension.

    Returns
    -------
    DTensor
        Unflattened DTensor.

    Raises
    ------
    ValueError
        If validation errors occur during unflattening or if placement argument
        constraints are violated.
    NotImplementedError
        If the dimension to be unflattened is sharded.
    """
    return _ShardWiseUnflattenImpl.apply(x, dim, sizes)


def _shardwise_unflatten_sharded_impl(input: DTensor, dim: int, sizes: tuple[int, ...]) -> DTensor:
    """Unflatten a sharded DTensor along a specified dimension.

    This function splits a single dimension into multiple dimensions, similar to
    torch.Tensor.unflatten, but designed for sharded DTensors. The input must be
    sharded along the specified dimension, and sizes[0] must be evenly divisible
    by the device mesh size so that the resulting DTensor is again sharded along
    the same dimension.

    This is the inverse operation of _shardwise_flatten_sharded_impl.

    Args:
        input: Input DTensor to unflatten. Must be sharded along `dim`.
        dim: The dimension to unflatten. Must correspond to a sharded dimension
            in the input's placements. If negative, wraps around.
        sizes: Tuple of integers specifying the shape to unflatten into.
            The product of sizes must equal input.shape[dim].
            sizes[0] must be evenly divisible by the device mesh size.

    Returns:
        DTensor with the specified dimension unflattened into multiple dimensions.
        The output shape is input.shape[:dim] + sizes + input.shape[dim+1:].

    Raises:
        TypeError: If input is not a DTensor, dim is not an int, or sizes
            is not a tuple.
        ValueError: If input has Partial placements, if sizes[0] is not
            evenly shardable, if the product of sizes doesn't match the
            dim size, or if input is not sharded along dim.
    """
    if not isinstance(input, DTensor):
        raise TypeError(f"Expected DTensor, got {type(input)}")
    if not isinstance(dim, int):
        raise TypeError(f"Expected int for dim, got {type(dim)}")
    if not isinstance(sizes, tuple):
        raise TypeError(f"Expected tuple for sizes, got {type(sizes)}")
    if len(sizes) < 2:
        raise ValueError("Must provide at least two dimensions for unflattening")

    ndim = input.ndim

    # Normalize dimension
    if dim < 0:
        dim = ndim + dim

    if not (0 <= dim < ndim):
        raise ValueError(f"dim {dim} out of range for {ndim}D tensor")

    device_mesh = input.device_mesh
    placements = input.placements

    i_mesh_dim_shard_dim = None
    for i_dim_device_mesh, placement in enumerate(placements):
        if isinstance(placement, Partial):
            raise ValueError("Partial placements are not supported")
        elif isinstance(placement, Shard):
            if input.shape[placement.dim] % device_mesh.shape[i_dim_device_mesh] != 0:
                raise ValueError(
                    f"Uneven sharding tensor dimension {placement.dim} of size {input.shape[placement.dim]} "
                    f"along device mesh dimension {i_dim_device_mesh} of size {device_mesh.shape[i_dim_device_mesh]} is not supported"
                )
            if placement.dim == dim:
                i_mesh_dim_shard_dim = i_dim_device_mesh

    if i_mesh_dim_shard_dim is None:
        raise ValueError(f"input is not sharded along dim {dim}")

    size_expected = math.prod(sizes)
    if size_expected != input.shape[dim]:
        raise ValueError(f"Expected size {size_expected} but got {input.shape[dim]}")

    size_group = device_mesh.size(i_mesh_dim_shard_dim)

    # sizes[0] will become the new shape[dim] and should be evenly sharded
    if sizes[0] % size_group != 0:
        raise ValueError(
            f"sizes[0] {sizes[0]} must be evenly sharded along device mesh dimension {i_mesh_dim_shard_dim} of size {size_group}"
        )

    # input.shape[dim] // size_group flattened and sharded into: (sizes[0] // size_group, sizes[1:])
    output_local = input.to_local().unflatten(dim, (sizes[0] // size_group, *sizes[1:]))

    shape_output = input.shape[:dim] + sizes + input.shape[dim + 1 :]
    strides_output = update_exhaustive_strides(output_local.shape, output_local.stride(), shape_output)

    # Adjust Shard dim indices for dimensions shifted by the unflatten.
    # Splitting dim into len(sizes) parts adds (len(sizes) - 1) new dims,
    # so any Shard(d) with d > dim must shift up by that amount.
    n_dims_added = len(sizes) - 1
    output_placements: tuple[Placement, ...] = tuple(
        Shard(p.dim + n_dims_added) if isinstance(p, Shard) and p.dim > dim else p for p in placements
    )

    output: DTensor = DTensor.from_local(
        output_local, device_mesh, output_placements, shape=shape_output, stride=strides_output
    )

    return output


def _shardwise_flatten_sharded_impl(input: DTensor, start_dim: int, end_dim: int) -> DTensor:
    """Flatten consecutive dimensions of a sharded DTensor.

    This function flattens dimensions from start_dim to end_dim (inclusive) into
    a single dimension. The input must be sharded along start_dim, and the
    sharding is preserved on the flattened output dimension.

    This is the inverse operation of _shardwise_unflatten_sharded_impl.

    Args:
        input: Input DTensor to flatten. Must be sharded along `start_dim`.
        start_dim: First dimension to flatten.
        end_dim: Last dimension to flatten (inclusive). If negative, wraps around.

    Returns:
        DTensor with dimensions [start_dim, end_dim] flattened into a single
        dimension at position start_dim.

    Raises:
        TypeError: If input is not a DTensor, or start_dim/end_dim are not int.
        ValueError: If input has Partial placements, if start_dim is not sharded,
            or if dimension indices are invalid.
    """
    if not isinstance(input, DTensor):
        raise TypeError(f"Expected DTensor, got {type(input)}")
    if not isinstance(start_dim, int):
        raise TypeError(f"Expected int for start_dim, got {type(start_dim)}")
    if not isinstance(end_dim, int):
        raise TypeError(f"Expected int for end_dim, got {type(end_dim)}")

    ndim = input.ndim

    # Normalize dimensions
    if start_dim < 0:
        start_dim = ndim + start_dim
    if end_dim < 0:
        end_dim = ndim + end_dim

    if not (0 <= start_dim < ndim):
        raise ValueError(f"start_dim {start_dim} out of range for {ndim}D tensor")
    if not (0 <= end_dim < ndim):
        raise ValueError(f"end_dim {end_dim} out of range for {ndim}D tensor")
    if start_dim > end_dim:
        raise ValueError(f"start_dim {start_dim} must be <= end_dim {end_dim}")

    device_mesh = input.device_mesh
    placements = input.placements

    i_mesh_dim_shard_start = None
    for i_dim_device_mesh, placement in enumerate(placements):
        if isinstance(placement, Partial):
            raise ValueError("Partial placements are not supported")
        elif isinstance(placement, Shard):
            if input.shape[placement.dim] % device_mesh.shape[i_dim_device_mesh] != 0:
                raise ValueError(
                    f"Uneven sharding tensor dimension {placement.dim} of size {input.shape[placement.dim]} "
                    f"along device mesh dimension {i_dim_device_mesh} of size {device_mesh.shape[i_dim_device_mesh]} is not supported"
                )
            if placement.dim == start_dim:
                i_mesh_dim_shard_start = i_dim_device_mesh

    if i_mesh_dim_shard_start is None:
        raise ValueError(f"input is not sharded along start_dim {start_dim}")

    # Flatten locally
    output_local = input.to_local().flatten(start_dim=start_dim, end_dim=end_dim)

    # Compute global output shape
    flattened_size = math.prod(input.shape[start_dim : end_dim + 1])
    shape_output = input.shape[:start_dim] + (flattened_size,) + input.shape[end_dim + 1 :]
    strides_output = update_exhaustive_strides(output_local.shape, output_local.stride(), shape_output)

    # Adjust Shard dim indices for dimensions shifted by the flatten.
    # Merging dims [start_dim, end_dim] removes (end_dim - start_dim) dims,
    # so any Shard(d) with d > end_dim must shift down by that amount.
    n_dims_removed = end_dim - start_dim
    output_placements: tuple[Placement, ...] = tuple(
        Shard(p.dim - n_dims_removed) if isinstance(p, Shard) and p.dim > end_dim else p for p in placements
    )

    output: DTensor = DTensor.from_local(
        output_local, device_mesh, output_placements, shape=shape_output, stride=strides_output
    )

    return output


class ShardwiseUnflattenShardedImpl(torch.autograd.Function):
    """Autograd function to unflatten a sharded DTensor along a specified axis.

    This function performs an unflatten operation on a DTensor while preserving
    the sharding semantics. The input must be sharded along the specified axis,
    and the first element of `sizes` must be evenly divisible by the device mesh
    size along the sharding dimension.

    Example:
        If input has global shape (B, N) sharded along axis=1 across 2 ranks,
        and sizes=(K, W), the output will have global shape (B, K, W) still
        sharded along axis=1. Each rank holds (B, K//2, W) locally.
    """

    @staticmethod
    def forward(ctx: FunctionCtx, input: DTensor, axis: int, sizes: tuple[int, ...]) -> DTensor:
        """Forward pass: unflatten the DTensor along the specified axis.

        Args:
            ctx: Autograd context for saving tensors and metadata for backward.
            input: Input DTensor to unflatten. Must be sharded along `axis`.
            axis: The axis along which to unflatten. Must correspond to a
                sharded dimension in the input's placements.
            sizes: Tuple of integers specifying the new shape for the unflattened
                dimension. The product of sizes must equal input.shape[axis].
                sizes[0] must be evenly divisible by the device mesh size.

        Returns:
            DTensor with the specified axis unflattened into multiple dimensions.
            The output shape is input.shape[:axis] + sizes + input.shape[axis+1:].
            Sharding is preserved along the first unflattened dimension.

        Raises:
            TypeError: If input is not a DTensor, axis is not an int, or sizes
                is not a tuple.
            ValueError: If input has Partial placements, if sizes[0] is not
                evenly shardable, if the product of sizes doesn't match the
                axis dimension, or if input is not sharded along axis.
        """
        output = _shardwise_unflatten_sharded_impl(input, dim=axis, sizes=sizes)
        ctx.axis = axis
        ctx.sizes = sizes
        return output

    @staticmethod
    def backward(ctx: FunctionCtx, grad_output: DTensor) -> tuple[DTensor, None, None]:
        """Backward pass: flatten the gradient back to the original input shape.

        Args:
            ctx: Autograd context containing saved tensors and metadata from forward.
            grad_output: Gradient DTensor with respect to the forward output.
                Must have the same device_mesh, placements, and shape as the
                forward output.

        Returns:
            Tuple of (grad_input, None, None) where grad_input is the gradient
            with respect to the input DTensor, and the None values correspond
            to the non-differentiable axis and sizes arguments.

        Raises:
            TypeError: If grad_output is not a DTensor.
            ValueError: If grad_output has mismatched device_mesh, placements,
                or shape compared to the forward output.
        """
        axis = ctx.axis
        sizes = ctx.sizes
        end_dim = axis + len(sizes) - 1
        grad_input = _shardwise_flatten_sharded_impl(grad_output, start_dim=axis, end_dim=end_dim)
        return grad_input, None, None


def shardwise_unflatten_sharded(input: DTensor, axis: int, sizes: tuple[int, ...]) -> DTensor:
    """Unflatten a sharded DTensor along a specified axis while preserving sharding.

    This function reshapes a DTensor by splitting a single dimension into multiple
    dimensions, similar to torch.Tensor.unflatten, but designed to work correctly
    with distributed tensors that are sharded along the unflattened axis.

    The key constraint is that the first element of `sizes` must be evenly
    divisible by the number of ranks sharding the axis, ensuring the resulting
    tensor can maintain valid sharding semantics.

    Args:
        input: Input DTensor to unflatten. Must be sharded along `axis` with
            even sharding (no remainder when dividing by device mesh size).
        axis: The axis along which to unflatten. Must be a sharded dimension.
        sizes: Tuple of integers specifying the shape to unflatten into.
            Must satisfy: prod(sizes) == input.shape[axis] and
            sizes[0] % device_mesh_size == 0.

    Returns:
        DTensor with shape input.shape[:axis] + sizes + input.shape[axis+1:].
        The tensor remains sharded along the same device mesh dimension,
        with the sharding now applying to the first element of sizes.

    Example:
        >>> # input: DTensor of global shape (4, 128) sharded on axis=1 across 2 ranks
        >>> # Each rank holds (4, 64) locally
        >>> output = shardwise_unflatten_sharded(input, axis=1, sizes=(16, 8))
        >>> # output: DTensor of global shape (4, 16, 8) sharded on axis=1
        >>> # Each rank holds (4, 8, 8) locally
    """
    return ShardwiseUnflattenShardedImpl.apply(input, axis, sizes)


class ShardwiseFlattenShardedImpl(torch.autograd.Function):
    """Autograd function to flatten consecutive dimensions of a sharded DTensor.

    This function performs a flatten operation on a DTensor while preserving
    the sharding semantics. The input must be sharded along start_dim, and the
    sharding is preserved on the flattened output dimension.

    Example:
        If input has global shape (B, K, W) sharded along axis=1 across 2 ranks,
        and we flatten dims 1 and 2, the output will have global shape (B, K*W)
        still sharded along axis=1. Each rank holds (B, K*W//2) locally.
    """

    @staticmethod
    def forward(ctx: FunctionCtx, input: DTensor, start_dim: int, end_dim: int) -> DTensor:
        """Forward pass: flatten the DTensor along the specified dimensions.

        Args:
            ctx: Autograd context for saving tensors and metadata for backward.
            input: Input DTensor to flatten. Must be sharded along `start_dim`.
            start_dim: First dimension to flatten.
            end_dim: Last dimension to flatten (inclusive). If negative, wraps around.

        Returns:
            DTensor with dimensions [start_dim, end_dim] flattened into a single
            dimension at position start_dim.

        Raises:
            TypeError: If input is not a DTensor, or start_dim/end_dim are not int.
            ValueError: If input has Partial placements, if start_dim is not sharded,
                or if dimension indices are invalid.
        """
        # Normalize end_dim for storing in context
        ndim = input.ndim
        if end_dim < 0:
            end_dim = ndim + end_dim

        output = _shardwise_flatten_sharded_impl(input, start_dim=start_dim, end_dim=end_dim)

        # Save the sizes of the flattened dimensions for backward (unflatten)
        ctx.start_dim = start_dim if start_dim >= 0 else ndim + start_dim
        ctx.sizes = tuple(input.shape[ctx.start_dim : end_dim + 1])
        return output

    @staticmethod
    def backward(ctx: FunctionCtx, grad_output: DTensor) -> tuple[DTensor, None, None]:
        """Backward pass: unflatten the gradient back to the original input shape.

        Args:
            ctx: Autograd context containing saved tensors and metadata from forward.
            grad_output: Gradient DTensor with respect to the forward output.
                Must have the same device_mesh, placements, and shape as the
                forward output.

        Returns:
            Tuple of (grad_input, None, None) where grad_input is the gradient
            with respect to the input DTensor, and the None values correspond
            to the non-differentiable start_dim and end_dim arguments.

        Raises:
            TypeError: If grad_output is not a DTensor.
            ValueError: If grad_output has mismatched device_mesh, placements,
                or shape compared to the forward output.
        """
        start_dim = ctx.start_dim
        sizes = ctx.sizes
        grad_input = _shardwise_unflatten_sharded_impl(grad_output, dim=start_dim, sizes=sizes)
        return grad_input, None, None


def shardwise_flatten_sharded(input: DTensor, start_dim: int, end_dim: int) -> DTensor:
    """Flatten consecutive dimensions of a sharded DTensor while preserving sharding.

    This function reshapes a DTensor by merging multiple dimensions into a single
    dimension, similar to torch.Tensor.flatten, but designed to work correctly
    with distributed tensors that are sharded along the start_dim.

    This is the inverse operation of shardwise_unflatten_sharded.

    Args:
        input: Input DTensor to flatten. Must be sharded along `start_dim` with
            even sharding (no remainder when dividing by device mesh size).
        start_dim: First dimension to flatten. Must be a sharded dimension.
        end_dim: Last dimension to flatten (inclusive). If negative, wraps around.

    Returns:
        DTensor with dimensions [start_dim, end_dim] merged into a single dimension.
        The tensor remains sharded along the same device mesh dimension.

    Example:
        >>> # input: DTensor of global shape (4, 16, 8) sharded on axis=1 across 2 ranks
        >>> # Each rank holds (4, 8, 8) locally
        >>> output = shardwise_flatten_sharded(input, start_dim=1, end_dim=2)
        >>> # output: DTensor of global shape (4, 128) sharded on axis=1
        >>> # Each rank holds (4, 64) locally
    """
    return ShardwiseFlattenShardedImpl.apply(input, start_dim, end_dim)
