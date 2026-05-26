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
from torch.distributed.tensor import DTensor, Partial, Replicate, Shard

from boltz.distributed.utils import update_exhaustive_strides


class _ShardedSumImpl(torch.autograd.Function):
    """Distributed implementation of sharded summation aggregation on pair input using DTensors."""

    @staticmethod
    def forward(ctx, x: DTensor, dim: tuple[int, ...] | int, keepdim: bool = False) -> DTensor:
        """Forward pass of distributed sharded summation aggregation.

        Parameters
        ----------
        ctx : torch.autograd.function.BackwardCFrame
            Context object for saving information needed in backward pass.
        x : DTensor
            Input tensor. Can have any shape and sharding strategy.
        dim : tuple[int, ...] | int
            Dimensions to reduce over. All of which must be sharded.
        keepdim : bool, default=False
            Whether to keep the reduced dimensions in the output.

        Returns
        -------
        DTensor
            Output tensor with reduced dimensions, maintaining the same device mesh
            and placement strategy as the input tensor.

        Raises
        ------
        TypeError
            If input is not a DTensor.
        ValueError
            If Partial placements are used (not supported), or if dims are invalid.
        """
        if not isinstance(x, DTensor):
            raise TypeError(f"Input 'x' must be of type DTensor. Got type {type(x)}.")

        device_mesh = x.device_mesh
        input_placements = x.placements
        input_shape = x.shape

        # Validate placements and cache sharded dims
        sharded_dims = set()
        for i_dim_device_mesh, placement in enumerate(input_placements):
            if isinstance(placement, Partial):
                raise ValueError("Partial placements are not supported")
            elif isinstance(placement, Shard):
                if input_shape[placement.dim] % device_mesh.shape[i_dim_device_mesh] != 0:
                    raise ValueError(
                        f"Uneven sharding tensor dimension {placement.dim} of size {input_shape[placement.dim]} "
                        f"along device mesh dimension {i_dim_device_mesh} of size {device_mesh.shape[i_dim_device_mesh]} is not supported"
                    )
                sharded_dims.add(placement.dim)
        if isinstance(dim, int):
            dim = (dim,)
        dims = dim

        if not dims:
            raise ValueError("Received empty dims argument")

        reduced_dims = []
        for dim in dims:
            if dim < 0:
                dim = len(input_shape) + dim
            if dim >= len(input_shape):
                raise ValueError(f"Input tensor has {len(input_shape)} dimensions but got dims {dims}")
            if dim not in sharded_dims:
                raise ValueError(f"Expected all dims are sharded but got {dims} for placements: {input_placements}")
            reduced_dims.append(dim)
        reduced_dims = tuple(sorted(reduced_dims))

        # remap sharded dims due to shape change when keepdim=False
        if not keepdim:
            map_dims = {}
            counter = 0
            for dim in range(x.ndim):
                map_dims[dim] = dim - counter
                counter += dim in reduced_dims

        x_local = x.to_local()
        output_local = torch.sum(x_local, dim=reduced_dims, keepdim=keepdim)  # new copy

        output_placements = []
        for placement, placement_group in zip(input_placements, device_mesh.get_all_groups()):
            if isinstance(placement, Shard) and placement.dim in reduced_dims:
                torch.distributed.all_reduce(
                    output_local,
                    op=torch.distributed.ReduceOp.SUM,
                    group=placement_group,
                )
                output_placements.append(Replicate())

            elif keepdim:  # Shortcut for keepdim=True
                output_placements.append(placement)
                continue

            # Shift placement dimensions when keepdim=False
            elif isinstance(placement, Shard):
                output_placements.append(Shard(map_dims[placement.dim]))
            elif isinstance(placement, Replicate):
                output_placements.append(placement)

        if x.requires_grad:
            ctx.device_mesh = device_mesh
            ctx.input_placements = input_placements
            ctx.input_local_shape = x_local.shape
            ctx.input_shape = x.shape
            ctx.input_stride = x.stride()
            ctx.keepdim = keepdim
            ctx.reduced_dims = reduced_dims

        # Compute output shape and stride
        # Shape stays the same as if keepdim=True, just reduced dimensions become size 1
        shape_output = list(x.shape)
        for dim in reduced_dims:
            shape_output[dim] = 1
        shape_output = tuple(shape_output)
        # Use update_exhaustive_strides to compute new strides
        strides_output = update_exhaustive_strides(x.shape, x.stride(), shape_output)
        if not keepdim:
            shape_output_reduced = []
            strides_output_reduced = []
            for i_dim, dim in enumerate(shape_output):
                if i_dim not in reduced_dims:
                    shape_output_reduced.append(dim)
                    strides_output_reduced.append(strides_output[i_dim])
            shape_output = tuple(shape_output_reduced)
            strides_output = tuple(strides_output_reduced)

        out = DTensor.from_local(
            output_local,
            device_mesh=device_mesh,
            placements=output_placements,
            shape=shape_output,
            stride=strides_output,
        )
        return out

    @staticmethod
    def backward(ctx, grad_output: DTensor) -> tuple[DTensor | None, None, None]:
        """Backward pass of distributed sharded summation aggregation.

        The gradient of sum(x, dims) with respect to x is simple broadcasting:
        - grad_output is broadcasted to the original tensor shape
        - Since we only reduce over sharded dimensions, the gradient is just broadcasted back

        Parameters
        ----------
        ctx : torch.autograd.function.BackwardCFrame
            Context object containing saved tensors and metadata from forward pass.
        grad_output : DTensor
            Gradients of the loss with respect to the output tensor.

        Returns
        -------
        tuple[DTensor | None, None, None]
            Gradients with respect to x, dims, and keepdim.
        """
        if not isinstance(grad_output, DTensor):
            raise TypeError(f"Input 'grad_output' must be of type DTensor but got type {type(grad_output)}.")

        if grad_output.device_mesh != ctx.device_mesh:
            raise ValueError(
                f"Input 'grad_output' must have the same device mesh as the input tensor. "
                f"Got device meshes {grad_output.device_mesh} and {ctx.device_mesh}."
            )

        dx_local = grad_output.to_local()
        if not ctx.keepdim:
            for dim in ctx.reduced_dims:
                dx_local = dx_local.unsqueeze(dim)
        dx_local = dx_local.expand(ctx.input_local_shape).clone(memory_format=torch.contiguous_format)

        dx = DTensor.from_local(
            dx_local,
            device_mesh=ctx.device_mesh,
            placements=ctx.input_placements,
            shape=ctx.input_shape,
            stride=ctx.input_stride,
        )
        return dx, None, None


def sharded_sum(x: DTensor, dim: tuple[int, ...] | int, keepdim: bool = False) -> DTensor:
    """Perform sharded summation aggregation.

    Behave similarly to torch.sum but expect all reduced dimensions to be sharded.

    Parameters
    ----------
    x: DTensor
        Input distributed tensor.
    dim: tuple[int, ...] | int
        Dimensions to reduce over.
    keepdim: bool, default=False
        Whether to keep the reduced dimensions in the output.

    Returns
    -------
    DTensor
        Output distributed tensor with dimensions and placements reduced.
    """
    return _ShardedSumImpl.apply(x, dim, keepdim)  # type: ignore
