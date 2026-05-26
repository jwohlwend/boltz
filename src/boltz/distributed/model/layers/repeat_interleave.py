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


class _ShardwiseRepeatInterleaveImpl(torch.autograd.Function):
    @staticmethod
    def forward(
        ctx: FunctionCtx,
        x: DTensor,
        repeats: int,
        dim: int,
    ) -> DTensor:
        """Forward pass.

        Parameters
        ----------
        ctx : FunctionCtx
            Context object for saving information needed in backward pass.
        x : DTensor
            Input DTensor.
        repeats : int
            Number of repetitions for each element.
        dim : int
            Dimension to repeat_interleave along.

        Returns
        -------
        DTensor
            DTensor after repeat_interleave operation.
        """
        # Type checking
        if not isinstance(x, DTensor):
            raise TypeError(f"Expected DTensor, got {type(x)}")
        if not isinstance(repeats, int):
            raise TypeError(f"Expected int for repeats, got {type(repeats)}")
        if not isinstance(dim, int):
            raise TypeError(f"Expected int for dim, got {type(dim)}")

        dim_normalized = dim if dim >= 0 else dim + x.ndim

        # Check placements and handle sharded dimensions
        for i_dim_device_mesh, placement in enumerate(x.placements):
            if isinstance(placement, Shard):
                # Check that sharded dimensions are evenly divided
                if x.shape[placement.dim] % x.device_mesh.shape[i_dim_device_mesh] != 0:
                    raise ValueError(
                        f"Uneven sharding tensor dimension {placement.dim} of size {x.shape[placement.dim]} "
                        f"along device mesh dimension {i_dim_device_mesh} of size {x.device_mesh.shape[i_dim_device_mesh]} is not supported"
                    )
            elif isinstance(placement, Partial):
                raise ValueError(f"Placements of type {Partial} are not supported")

        x_local = x.to_local()

        # Perform operation on local tensors
        output_local: Tensor = torch.repeat_interleave(x_local, repeats=repeats, dim=dim)

        # Compute output shape and stride
        shape_output = list(x.shape)
        shape_output[dim_normalized] = x.shape[dim_normalized] * repeats
        shape_output = tuple(shape_output)

        # Use update_exhaustive_strides to compute new strides
        strides_output = update_exhaustive_strides(output_local.shape, output_local.stride(), shape_output)

        # Create output DTensor using input tensor's device mesh and placements
        result: DTensor = DTensor.from_local(
            output_local,
            device_mesh=x.device_mesh,
            placements=x.placements,
            shape=shape_output,
            stride=strides_output,
        )

        # Save information for backward pass
        ctx.repeats = repeats
        ctx.dim_normalized = dim_normalized
        ctx.input_device_mesh = x.device_mesh
        ctx.input_placements = x.placements
        ctx.input_shape = x.shape
        ctx.input_stride = x.stride()
        ctx.output_shape = result.shape

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
            Gradient with respect to input, None for repeats and dim parameters.
        """
        # Check that grad_output has the expected shape, device_mesh and placements
        raise_if_incorrect_dtensor_metadata_args(
            dtensor_instance=grad_output,
            dtensor_name="_ShardwiseRepeatInterleaveImpl.backward grad_output",
            expected_shape=ctx.output_shape,
            expected_device_mesh=ctx.input_device_mesh,
            expected_placements=ctx.input_placements,
        )

        # Perform backward pass on local tensors
        grad_output_local = grad_output.to_local()

        # Reshape and sum to reverse the repeat_interleave operation
        # Get the original size along the dimension that was repeated
        original_size = grad_output_local.shape[ctx.dim_normalized] // ctx.repeats

        # Unflatten the repeated dimension and sum along the repeats dimension
        grad_unflattened = torch.unflatten(grad_output_local, ctx.dim_normalized, (original_size, ctx.repeats))
        grad_input_local = grad_unflattened.sum(dim=ctx.dim_normalized + 1)

        # Create output DTensor using the saved metadata
        grad_input = DTensor.from_local(
            grad_input_local,
            device_mesh=ctx.input_device_mesh,
            placements=ctx.input_placements,
            shape=ctx.input_shape,
            stride=ctx.input_stride,
        )

        return grad_input, None, None


def shardwise_repeat_interleave(x: DTensor, repeats: int, dim: int) -> DTensor:
    """Repeat elements of a DTensor along a specified dimension.

    This function repeats elements of a DTensor along the specified dimension.
    Each element along the specified dimension is repeated `repeats` times.

    Parameters
    ----------
    x : DTensor
        Input DTensor to repeat_interleave.
    repeats : int
        Number of repetitions for each element.
    dim : int
        Dimension to repeat_interleave along.

    Returns
    -------
    DTensor
        DTensor after repeat_interleave operation.

    Raises
    ------
    TypeError
        If inputs are not of correct type.
    ValueError
        If validation errors occur.
    """
    return _ShardwiseRepeatInterleaveImpl.apply(x, repeats, dim)
