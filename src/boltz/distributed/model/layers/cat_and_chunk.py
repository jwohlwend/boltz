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
from torch import Tensor
from torch.autograd.function import FunctionCtx
from torch.distributed.tensor import DTensor, Partial, Shard

from boltz.distributed.model.layers.dtensor_metadata_tools import (
    raise_if_incorrect_dtensor_metadata_args,
)
from boltz.distributed.utils import update_exhaustive_strides


def _shardwise_chunk(
    x: DTensor,
    chunks: int,
    dim: int,
) -> tuple[DTensor, ...]:
    """Generalized shardwise chunking operation with validation.

    This function performs input validation and chunking operation.

    Parameters
    ----------
    x : DTensor
        Input DTensor to chunk.
    chunks : int
        Number of chunks to split the dimension into.
    dim : int
        Dimension to chunk along.

    Returns
    -------
    tuple[DTensor, ...]
        Tuple of DTensor instances after chunking.

    Raises
    ------
    TypeError
        See raise_if_incorrect_dtensor_metadata_args for more details.
    ValueError
        Checks on the input x and parameters.
    """
    dim_normalized = dim if dim >= 0 else dim + x.ndim
    # do not allow chunking on a dimension with a sharded placement
    for i_dim_device_mesh, placement in enumerate(x.placements):
        if isinstance(placement, Shard):
            if placement.dim == dim_normalized:
                raise NotImplementedError(
                    f"Chunking along dimension {dim} shared by device_mesh axis {i_dim_device_mesh} is not supported"
                )
            if x.shape[placement.dim] % x.device_mesh.shape[i_dim_device_mesh] != 0:
                raise ValueError(
                    f"Uneven sharding tensor dimension {placement.dim} of size {x.shape[placement.dim]} "
                    f"along device mesh dimension {i_dim_device_mesh} of size {x.device_mesh.shape[i_dim_device_mesh]} is not supported"
                )
        elif isinstance(placement, Partial):
            raise ValueError(f"Placements of type {Partial} are not supported")

    x_local = x.to_local()

    # Perform operation on local tensors
    x_chunks_local: tuple[Tensor, ...] = x_local.chunk(chunks=chunks, dim=dim)

    shapes_output = []
    strides_output = []

    for chunk in x_chunks_local:
        shape_output = list(x.shape)
        shape_output[dim_normalized] = chunk.shape[dim_normalized]
        # we try to be as consistent with the original stride as possible
        # but in principle there is no way to keep the resulting chunked DTensor.full_tensor
        # as views of the original DTensor.full_tensor because upon calling this function,
        # the latter is not materialized in memory yet
        stride_output = update_exhaustive_strides(x.shape, x.stride(), shape_output)
        shapes_output.append(tuple(shape_output))
        strides_output.append(tuple(stride_output))

    # Create output tuple of DTensor using input tensor's device mesh and placements
    # leave the dim check to the torch.Tensor.chunk
    x_in_chunks: tuple[DTensor, ...] = tuple(
        [
            DTensor.from_local(
                chunk,
                device_mesh=x.device_mesh,
                placements=x.placements,
                shape=shapes_output[i],
                stride=strides_output[i],
            )
            for i, chunk in enumerate(x_chunks_local)
        ]
    )
    return x_in_chunks


def _shardwise_cat(
    *inputs: DTensor, dim: int, shape: tuple[int, ...] | None = None, stride: tuple[int, ...] | None = None
) -> DTensor:
    """Generalized shardwise concatenation operation with validation.

    This function performs input validation and concatenation operation.

    Parameters
    ----------
    *inputs : DTensor
        Variable number of input DTensors to concatenate with dim as the last argument
    dim : int
        Dimension to concatenate along
    shape : tuple[int, ...], optional
        Shape of the output DTensor. If not provided, will infer from the inputs.
    stride : tuple[int, ...], optional
        Stride of the output DTensor. If not provided, will infer from the inputs.

    Returns
    -------
    DTensor
        Concatenated DTensor.

    Raises
    ------
    TypeError, ValueError
        Checks on the input tensors and parameters.
    """
    # Validate inputs
    if len(inputs) == 0:
        raise ValueError("Cannot concatenate empty list of tensors.")

    if (shape is None) != (stride is None):
        raise ValueError("Either both shape and stride must be provided or neither")

    first_input = inputs[0]

    placements = first_input.placements
    device_mesh = first_input.device_mesh

    # DTensor Shard(dim) is always normalized to be non-negative
    dim_normalized = dim if dim >= 0 else dim + first_input.ndim
    # Check that concatenation dimension is not sharded
    # and there are no unevenly sharded tensor axes
    for i_dim_device_mesh, placement in enumerate(placements):
        if isinstance(placement, Shard):
            i_dim_tensor = placement.dim
            if i_dim_tensor == dim_normalized:
                # unevenly sharded inputs are not supported
                raise NotImplementedError(
                    f"Concatenation along dimension {dim} shared by device_mesh axis {i_dim_device_mesh} is not supported"
                )
            if first_input.shape[i_dim_tensor] % device_mesh.shape[i_dim_device_mesh] != 0:
                raise ValueError(
                    f"Uneven sharding tensor dimension {i_dim_tensor} of size {first_input.shape[i_dim_tensor]} "
                    f"along device mesh dimension {i_dim_device_mesh} of size {device_mesh.shape[i_dim_device_mesh]} is not supported"
                )
        elif isinstance(placement, Partial):
            raise ValueError(f"Placements of type {Partial} are not supported")

    # Check that all inputs have compatible metadata
    for i, input_tensor in enumerate(inputs[1:], 1):
        raise_if_incorrect_dtensor_metadata_args(
            dtensor_instance=input_tensor,
            dtensor_name=f"_shardwise_cat inputs[{i}]",
            expected_device_mesh=device_mesh,
            expected_placements=placements,
        )

    if shape is None:
        shape = list(first_input.shape)
        shape[dim_normalized] = sum(input_dtensor.shape[dim_normalized] for input_dtensor in inputs)
        shape = tuple(shape)
        stride = update_exhaustive_strides(first_input.shape, first_input.stride(), shape)

    # passing the previous checks implies: all inputs have same device_mesh and same
    # placements. The following torch.cat will check for the local shards in the list
    # for the same shape. If the torch.cat is successful, it means that the inputs
    # will have same shape along all the sharded axes (given that dim can't be sharded)
    # and hence with evenly sharded axes (because first.input is evenly sharded)

    # Perform operation on local tensors
    # leave the dim and shape check to the torch.Tensor.cat
    output_local: Tensor = torch.cat([x.to_local() for x in inputs], dim=dim)

    # Create output DTensor using first input's device mesh and placements
    output: DTensor = DTensor.from_local(
        output_local,
        device_mesh=device_mesh,
        placements=placements,
        shape=shape,
        stride=stride,
    )

    return output


class _ShardwiseChunkImpl(torch.autograd.Function):
    @staticmethod
    def forward(
        ctx: FunctionCtx,
        x: DTensor,
        chunks: int,
        dim: int = -1,
    ) -> tuple[DTensor, ...]:
        """Forward pass.

        Parameters
        ----------
        ctx : FunctionCtx
            Context object for saving information needed in backward pass.
        x : DTensor
            Input DTensor.
        chunks : int
            Number of chunks to split the dimension into.
        dim : int, optional
            Dimension to chunk along. Default is -1 (last dimension).

        Returns
        -------
        tuple[DTensor, ...]
            Tuple of DTensor instances after chunking.
        """
        # Perform chunking with built-in validation
        result = _shardwise_chunk(x, chunks, dim)

        ctx.dim = dim
        ctx.device_mesh_input = x.device_mesh
        ctx.placements_input = x.placements
        ctx.shape_input = x.shape
        ctx.stride_input = x.stride()

        return result

    @staticmethod
    def backward(
        ctx: FunctionCtx,
        *grad_outputs: tuple[DTensor, ...],
    ) -> tuple[DTensor, None, None]:
        """Backward pass.

        Parameters
        ----------
        ctx : FunctionCtx
            Context object containing saved tensors and metadata from forward pass.
        grad_outputs : tuple[DTensor, ...]
            Gradient of the loss with respect to each component of the output.

        Returns
        -------
        tuple[DTensor, None, None]
            Gradient with respect to input, None for chunks and dim parameters.
        """
        # Use cat operation (inverse of chunk) for backward pass with built-in validation
        # no need for shape consistency check here because:
        # 1. torch.autograd.backward will check output's shape against grad_outputs shape
        # 2. the underlying torch.cat will check for shape along other axes than "dim"
        # 3. mismatching shape in the cat grad will be caught upon attaching to the input tensor
        grad_x = _shardwise_cat(*grad_outputs, dim=ctx.dim, shape=ctx.shape_input, stride=ctx.stride_input)

        return grad_x, None, None


class _ShardWiseCatImpl(torch.autograd.Function):
    @staticmethod
    def forward(ctx: FunctionCtx, *inputs) -> DTensor:
        """Forward pass.

        Parameters
        ----------
        ctx : FunctionCtx
            Context object for saving information needed in backward pass.
        *inputs : DTensor
            Variable number of input DTensors to concatenate with dim as the last argument

        Returns
        -------
        DTensor
            Concatenated DTensor.
        """
        # Perform concatenation with built-in validation
        tensors_to_cat = inputs[:-1]
        dim = inputs[-1]
        result = _shardwise_cat(*tensors_to_cat, dim=dim)

        # Save metadata for backward pass
        # shardwise_cat guarantee same device_mesh and placements as in the inputs
        ctx.device_mesh = result.device_mesh
        ctx.placements = result.placements
        ctx.dim = dim
        ctx.n_chunks = len(tensors_to_cat)
        ctx.shapes_and_strides_input = [
            (input_dtensor.shape, input_dtensor.stride()) for input_dtensor in tensors_to_cat
        ]

        # Save the sizes of each input tensor in the concatenation dimension for backward pass
        ctx.split_sizes = [input_tensor.shape[dim] for input_tensor in tensors_to_cat]

        return result

    @staticmethod
    def backward(
        ctx: FunctionCtx,
        grad_output: DTensor,
    ) -> tuple[DTensor, ...]:
        """Backward pass.

        Parameters
        ----------
        ctx : FunctionCtx
            Context object containing saved tensors and metadata from forward pass.
        grad_output : DTensor
            Gradient of the loss with respect to the output.

        Returns
        -------
        tuple[DTensor, ...]
            Gradients with respect to each input tensor, None for dim parameter.
        """
        # Verify grad_output has the same device_mesh and placements as expected
        raise_if_incorrect_dtensor_metadata_args(
            dtensor_instance=grad_output,
            dtensor_name="_ShardWiseCatImpl.backward grad_output",
            expected_device_mesh=ctx.device_mesh,
            expected_placements=ctx.placements,
        )

        # Use local torch.split for backward pass
        grad_output_local = grad_output.to_local()
        grad_inputs_local = torch.split(grad_output_local, ctx.split_sizes, dim=ctx.dim)

        # Wrap each split back into DTensor using the saved metadata
        grad_inputs = tuple(
            DTensor.from_local(
                grad_local,
                device_mesh=ctx.device_mesh,
                placements=ctx.placements,
                shape=ctx.shapes_and_strides_input[i][0],
                stride=ctx.shapes_and_strides_input[i][1],
            )
            for i, grad_local in enumerate(grad_inputs_local)
        )

        return (*grad_inputs, None)


def shardwise_chunk(x: DTensor, chunks: int, dim: int = -1) -> tuple[DTensor, ...]:
    """Chunk a DTensor along a specified dimension.

    This function splits a DTensor into chunks along the specified dimension.
    The dimension must not be sharded on the device mesh.

    Parameters
    ----------
    x : DTensor
        Input DTensor to chunk.
    chunks : int
        Number of chunks to split the dimension into.
    dim : int, optional
        Dimension to chunk along. Default is -1 (last dimension).

    Returns
    -------
    tuple[DTensor, ...]
        Tuple of DTensor instances after chunking.

    Raises
    ------
    ValueError
        If the specified dimension is sharded or other validation errors.
    """
    return _ShardwiseChunkImpl.apply(x, chunks, dim)


def shardwise_cat(inputs: list[DTensor], dim: int = -1) -> DTensor:
    """Concatenate DTensors along a specified dimension.

    This function concatenates DTensors along the specified dimension.
    The dimension must not be sharded on the device mesh, and all tensors
    must have compatible shapes except in the concatenation dimension.

    Parameters
    ----------
    inputs : list[DTensor]
        List of DTensors to concatenate.
    dim : int, optional
        Dimension to concatenate along. Default is -1 (last dimension).

    Returns
    -------
    DTensor
        Concatenated DTensor.

    Raises
    ------
    ValueError
        If the specified dimension is sharded or other validation errors.
    """
    if not inputs:
        raise ValueError("Cannot concatenate empty list of tensors.")

    return _ShardWiseCatImpl.apply(*inputs, dim)
