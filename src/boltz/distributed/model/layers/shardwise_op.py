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


from enum import Enum
from typing import Any, Optional

import torch
import torch.nn.functional as F
from torch.autograd.function import FunctionCtx
from torch.distributed.tensor import DTensor, Partial, Shard

from boltz.distributed.model.layers.dtensor_metadata_tools import (
    raise_if_incorrect_dtensor_metadata_args,
)
from boltz.distributed.utils import LayoutRightMap, update_exhaustive_strides


class _ShardwiseSumImpl(torch.autograd.Function):
    @staticmethod
    def forward(
        ctx: FunctionCtx,
        x: DTensor,
        dim: int,
        keepdim: Optional[bool] = None,
    ) -> DTensor:
        """Forward pass.

        Parameters
        ----------
        ctx : FunctionCtx
            Context object for saving information needed in backward pass.
        x : DTensor
            Input DTensor.
        dim : int
            The dimension to sum over.
        keepdim : Optional[bool]
            Whether to keep the dimension when summing.

        Returns
        -------
        DTensor
            DTensor after sum operation.
        """
        # Type checking
        if not isinstance(x, DTensor):
            raise TypeError(f"Expected DTensor, got {type(x)}")
        if not isinstance(dim, int):
            raise TypeError(f"Expected int for dim, got {type(dim)}")

        device_mesh_input = x.device_mesh
        placements_input = x.placements

        # Check placements and handle sharded dimensions
        actual_dim = dim if dim >= 0 else len(x.shape) + dim
        for i_dim_device_mesh, placement in enumerate(placements_input):
            if isinstance(placement, Partial):
                raise ValueError("Partial placements are not supported")
            elif isinstance(placement, Shard):
                # Check that sharded dimensions are evenly divided
                if x.shape[placement.dim] % device_mesh_input.shape[i_dim_device_mesh] != 0:
                    raise ValueError(
                        f"Uneven sharding tensor dimension {placement.dim} of size {x.shape[placement.dim]} "
                        f"along device mesh dimension {i_dim_device_mesh} of size {device_mesh_input.shape[i_dim_device_mesh]} is not supported"
                    )
                if placement.dim == actual_dim:
                    raise ValueError(f"Sum along sharded dimension {dim} is not supported")

        x_local = x.to_local()

        # Perform operation on local tensors
        output_local = (
            torch.sum(x_local, dim=dim, keepdim=keepdim) if keepdim is not None else torch.sum(x_local, dim=dim)
        )

        shape_output = list(x.shape)
        shape_output[actual_dim] = 1

        # keep the layout mapping but with a new shape
        strides_output = update_exhaustive_strides(x.shape, x.stride(), shape_output)
        if not keepdim:
            # remove the singleton dimension
            shape_output = shape_output[:actual_dim] + shape_output[actual_dim + 1 :]
            strides_output = strides_output[:actual_dim] + strides_output[actual_dim + 1 :]

        # Create output DTensor using input tensor's device mesh and placements
        result: DTensor = DTensor.from_local(
            output_local,
            device_mesh=device_mesh_input,
            placements=placements_input,
            shape=tuple(shape_output),
            stride=strides_output,
        )

        # Save information for backward pass
        if x.requires_grad:
            ctx.dim = dim
            ctx.keepdim = keepdim
            ctx.input_local_shape = x_local.shape
            ctx.device_mesh_input = device_mesh_input
            ctx.placements_input = placements_input
            ctx.output_shape = result.shape
            ctx.shape_input = x.shape
            ctx.stride_input = x.stride()

        return result

    @staticmethod
    def backward(
        ctx: FunctionCtx,
        *grad_outputs,
    ) -> tuple[DTensor, None, None]:
        """Backward pass.

        Parameters
        ----------
        ctx : FunctionCtx
            Context object containing saved tensors and metadata from forward pass.
        grad_outputs : tuple
            Gradients of the loss with respect to the output.

        Returns
        -------
        tuple[DTensor, None, None]
            Gradient with respect to input, None for dim and keepdim parameters.
        """
        grad_output = grad_outputs[0]

        # Check that grad_output has the expected shape, device_mesh and placements
        raise_if_incorrect_dtensor_metadata_args(
            dtensor_instance=grad_output,
            dtensor_name="_ShardwiseSumImpl.backward grad_output",
            expected_shape=ctx.output_shape,
            expected_device_mesh=ctx.device_mesh_input,
            expected_placements=ctx.placements_input,
        )

        # Perform backward pass on local tensors
        grad_output_local = grad_output.to_local()

        # For sum operation, gradient is broadcasted back to original shape
        input_local_shape = list(ctx.input_local_shape)
        if ctx.keepdim:
            dx_local = grad_output_local.expand(input_local_shape)
        else:
            grad_expanded = grad_output_local.unsqueeze(ctx.dim)
            dx_local = grad_expanded.expand(input_local_shape)

        # Create output DTensor using the saved metadata
        grad_input = DTensor.from_local(
            dx_local,
            device_mesh=ctx.device_mesh_input,
            placements=ctx.placements_input,
            shape=ctx.shape_input,
            stride=ctx.stride_input,
        )

        return grad_input, None, None


class _ShardwiseOneHotImpl(torch.autograd.Function):
    @staticmethod
    def forward(
        ctx: FunctionCtx,
        input: DTensor,
        num_classes: int = -1,
    ) -> DTensor:
        """Forward pass.

        Parameters
        ----------
        ctx : FunctionCtx
            Context object for saving information needed in backward pass.
        input : DTensor
            Input DTensor containing class indices.
        num_classes : int
            Number of classes for one-hot encoding. If -1, inferred from input.

        Returns
        -------
        DTensor
            DTensor after one-hot encoding.
        """
        # Type checking
        if not isinstance(input, DTensor):
            raise TypeError(f"Expected DTensor, got {type(input)}")
        if not isinstance(num_classes, int):
            raise TypeError(f"Expected int for num_classes, got {type(num_classes)}")

        device_mesh_input = input.device_mesh
        placements_input = input.placements

        # Check placements and handle sharded dimensions
        for i_dim_device_mesh, placement in enumerate(placements_input):
            if isinstance(placement, Partial):
                raise ValueError("Partial placements are not supported")
            elif isinstance(placement, Shard):
                # Check that sharded dimensions are evenly divided
                if input.shape[placement.dim] % device_mesh_input.shape[i_dim_device_mesh] != 0:
                    raise ValueError(
                        f"Uneven sharding tensor dimension {placement.dim} of size {input.shape[placement.dim]} "
                        f"along device mesh dimension {i_dim_device_mesh} of size {device_mesh_input.shape[i_dim_device_mesh]} is not supported"
                    )

        input_local = input.to_local()

        # Perform one-hot operation on local tensors
        output_local = F.one_hot(input_local, num_classes=num_classes)

        # Compute output shape and stride (one-hot adds a dimension at the end)
        shape_output = input.shape + (output_local.shape[-1],)

        # For one-hot, we append a new dimension, so we can use LayoutRightMap
        layout_right = LayoutRightMap(shape_output)
        strides_output = layout_right.strides

        # Create output DTensor
        result: DTensor = DTensor.from_local(
            output_local,
            device_mesh=device_mesh_input,
            placements=placements_input,
            shape=shape_output,
            stride=strides_output,
        )
        ctx.mark_non_differentiable(result)

        return result

    @staticmethod
    def backward(
        ctx: FunctionCtx,
        *grad_outputs,
    ) -> tuple[None, None]:
        """Backward pass.

        Parameters
        ----------
        ctx : FunctionCtx
            Context object containing saved tensors and metadata from forward pass.
        grad_outputs : tuple
            Gradients of the loss with respect to the output.

        Returns
        -------
        tuple[None, None]
            None gradients for input and num_classes (one_hot is not differentiable w.r.t. indices).
        """
        # one_hot is not differentiable with respect to the input indices
        # Return None for both input and num_classes gradients
        return None, None


class _ShardwiseDistogramImpl(torch.autograd.Function):
    @staticmethod
    def forward(
        ctx: FunctionCtx,
        d: DTensor,
        boundaries: torch.Tensor,
    ) -> DTensor:
        """Forward pass: bin distances into a distogram.

        Computes ``(d.unsqueeze(-1) > boundaries).sum(dim=-1).long()`` element-wise
        on the local shard and wraps the result back into a DTensor with the same
        placements and shape as the input.

        Parameters
        ----------
        ctx : FunctionCtx
            Context object for saving information needed in backward pass.
        d : DTensor
            Input DTensor of pairwise distances.
        boundaries : torch.Tensor
            1-D tensor of bin boundaries (not a DTensor).

        Returns
        -------
        DTensor
            Long DTensor of bin indices, same shape and placements as ``d``.
        """
        device_mesh = d.device_mesh
        placements = d.placements

        d_local = d.to_local()
        output_local = (d_local.unsqueeze(-1) > boundaries).sum(dim=-1).long()

        stride_output = update_exhaustive_strides(output_local.shape, output_local.stride(), d.shape)

        result: DTensor = DTensor.from_local(
            output_local,
            device_mesh=device_mesh,
            placements=placements,
            shape=d.shape,
            stride=stride_output,
        )
        ctx.mark_non_differentiable(result)

        return result

    @staticmethod
    def backward(
        ctx: FunctionCtx,
        *grad_outputs,
    ) -> tuple[None, None]:
        """Backward pass.

        Parameters
        ----------
        ctx : FunctionCtx
            Context object containing saved tensors and metadata from forward pass.
        grad_outputs : tuple
            Gradients of the loss with respect to the output.

        Returns
        -------
        tuple[None, None]
            None gradients for d and boundaries (distogram binning is not differentiable).
        """
        return None, None


class _ShardwiseSoftmaxImpl(torch.autograd.Function):
    @staticmethod
    def forward(
        ctx: FunctionCtx,
        x: DTensor,
        dim: int = -1,
    ) -> DTensor:
        """Forward pass.

        Parameters
        ----------
        ctx : FunctionCtx
            Context object for saving information needed in backward pass.
        x : DTensor
            Input DTensor.
        dim : int
            The dimension to apply softmax over.
        """
        device_mesh_input = x.device_mesh
        placements_input = x.placements

        # Check placements and handle sharded dimensions
        actual_dim = dim if dim >= 0 else len(x.shape) + dim
        for i_dim_device_mesh, placement in enumerate(placements_input):
            if isinstance(placement, Partial):
                raise ValueError("Partial placements are not supported")
            elif isinstance(placement, Shard):
                # Check that the softmax dim is not sharded - must be Replicate()
                if placement.dim == actual_dim:
                    raise ValueError(f"Softmax along sharded dimension {dim} is not supported")

                # Check that sharded dimensions are evenly divided
                if x.shape[placement.dim] % device_mesh_input.shape[i_dim_device_mesh] != 0:
                    raise ValueError(
                        f"Uneven sharding tensor dimension {placement.dim} of size {x.shape[placement.dim]} "
                        f"along device mesh dimension {i_dim_device_mesh} of size {device_mesh_input.shape[i_dim_device_mesh]} is not supported"
                    )

        x_local = x.to_local().detach().requires_grad_(x.requires_grad)

        with torch.enable_grad():
            output_local = F.softmax(x_local, dim=dim)

        shape_output = x.shape
        stride_output = x.stride()

        # Save information for backward pass
        if x.requires_grad:
            ctx.dim = dim
            ctx.device_mesh_input = device_mesh_input
            ctx.placements_input = placements_input
            ctx.output_shape = shape_output
            ctx.save_for_backward(x_local, output_local)  # need x_local here for autograd.grad() in bwd

        output: DTensor = DTensor.from_local(
            output_local.detach(),
            device_mesh=device_mesh_input,
            placements=placements_input,
            shape=shape_output,
            stride=stride_output,
        )

        return output

    @staticmethod
    def backward(
        ctx: FunctionCtx,
        grad_output: DTensor,
    ) -> tuple[DTensor, None]:
        """Backward pass.

        Parameters
        ----------
        ctx : FunctionCtx
            Context object containing saved tensors and metadata from forward pass.
        grad_output : DTensor
            Gradient of the loss with respect to the output.

        Returns
        -------
        tuple[DTensor, None]
            Gradient with respect to input, None for dim parameter.
        """
        # Check that grad_output has the expected shape, device_mesh and placements
        raise_if_incorrect_dtensor_metadata_args(
            dtensor_instance=grad_output,
            dtensor_name="_ShardwiseSoftmaxImpl.backward grad_output",
            expected_shape=ctx.output_shape,
            expected_device_mesh=ctx.device_mesh_input,
            expected_placements=ctx.placements_input,
        )

        # Perform backward pass on local tensors using saved subgraph
        grad_output_local = grad_output.to_local()
        x_local, softmax_output_local = ctx.saved_tensors

        (d_x_local,) = torch.autograd.grad(
            outputs=[softmax_output_local],
            inputs=[x_local],
            grad_outputs=[grad_output_local],
            retain_graph=False,
        )

        # Create output DTensor using the saved metadata
        grad_input = DTensor.from_local(
            d_x_local,
            device_mesh=ctx.device_mesh_input,
            placements=ctx.placements_input,
            shape=grad_output.shape,
            stride=grad_output.stride(),
        )

        return grad_input, None


class _ShardwiseLogSoftmaxImpl(torch.autograd.Function):
    @staticmethod
    def forward(
        ctx: FunctionCtx,
        x: DTensor,
        dim: int = -1,
    ) -> DTensor:
        """Forward pass.

        Parameters
        ----------
        ctx : FunctionCtx
            Context object for saving information needed in backward pass.
        x : DTensor
            Input DTensor.
        dim : int
            The dimension to apply log_softmax over.

        Returns
        -------
        DTensor
            DTensor after log_softmax operation.
        """
        # Type checking
        if not isinstance(x, DTensor):
            raise TypeError(f"Expected DTensor, got {type(x)}")
        if not isinstance(dim, int):
            raise TypeError(f"Expected int for dim, got {type(dim)}")

        device_mesh_input = x.device_mesh
        placements_input = x.placements

        # Check placements and handle sharded dimensions
        actual_dim = dim if dim >= 0 else len(x.shape) + dim
        for i_dim_device_mesh, placement in enumerate(placements_input):
            if isinstance(placement, Partial):
                raise ValueError("Partial placements are not supported")
            elif isinstance(placement, Shard):
                # Check that sharded dimensions are evenly divided
                if x.shape[placement.dim] % device_mesh_input.shape[i_dim_device_mesh] != 0:
                    raise ValueError(
                        f"Uneven sharding tensor dimension {placement.dim} of size {x.shape[placement.dim]} "
                        f"along device mesh dimension {i_dim_device_mesh} of size {device_mesh_input.shape[i_dim_device_mesh]} is not supported"
                    )
                if placement.dim == actual_dim:
                    raise ValueError(f"Log_softmax along sharded dimension {dim} is not supported, must be Replicate")

        x_local = x.to_local()

        # Perform operation on local tensors
        output_local = F.log_softmax(x_local, dim=dim)

        # Create output DTensor using input tensor's device mesh and placements
        shape_output = x.shape
        stride_output = x.stride()
        output: DTensor = DTensor.from_local(
            output_local,
            device_mesh=device_mesh_input,
            placements=placements_input,
            shape=shape_output,
            stride=stride_output,
        )

        # Save information for backward pass
        if x.requires_grad:
            ctx.dim = dim
            ctx.device_mesh_input = device_mesh_input
            ctx.placements_input = placements_input
            ctx.output_shape = output.shape
            ctx.save_for_backward(output_local)

        return output

    @staticmethod
    def backward(
        ctx: FunctionCtx,
        grad_output: DTensor,
    ) -> tuple[DTensor, None]:
        """Backward pass.

        Parameters
        ----------
        ctx : FunctionCtx
            Context object containing saved tensors and metadata from forward pass.
        grad_outputs : tuple
            Gradients of the loss with respect to the output.

        Returns
        -------
        tuple[DTensor, None]
            Gradient with respect to input, None for dim parameter.
        """
        # Check that grad_output has the expected shape, device_mesh and placements
        raise_if_incorrect_dtensor_metadata_args(
            dtensor_instance=grad_output,
            dtensor_name="_ShardwiseLogSoftmaxImpl.backward grad_output",
            expected_shape=ctx.output_shape,
            expected_device_mesh=ctx.device_mesh_input,
            expected_placements=ctx.placements_input,
        )

        # Perform backward pass on local tensors
        grad_output_local = grad_output.to_local()
        (log_softmax_output_local,) = ctx.saved_tensors

        # For log_softmax(x), the gradient is:
        # grad_input = grad_output - exp(log_softmax(x)) * sum(grad_output, dim=dim, keepdim=True)
        # This is equivalent to: grad_input = grad_output - softmax(x) * sum(grad_output, dim=dim, keepdim=True)
        grad_sum = grad_output_local.sum(dim=ctx.dim, keepdim=True)
        dx_local = grad_output_local - torch.exp(log_softmax_output_local) * grad_sum

        # Create output DTensor using the saved metadata
        grad_input = DTensor.from_local(
            dx_local,
            device_mesh=ctx.device_mesh_input,
            placements=ctx.placements_input,
            shape=grad_output.shape,
            stride=grad_output.stride(),
        )

        return grad_input, None


class _ShardwiseArgmaxImpl(torch.autograd.Function):
    @staticmethod
    def forward(
        ctx: FunctionCtx,
        x: DTensor,
        dim: int,
        keepdim: Optional[bool] = None,
    ) -> DTensor:
        """Forward pass for shardwise argmax."""
        if not isinstance(x, DTensor):
            raise TypeError(f"Expected DTensor, got {type(x)}")
        if not isinstance(dim, int):
            raise TypeError(f"Expected int for dim, got {type(dim)}")

        device_mesh_input = x.device_mesh
        placements_input = x.placements

        actual_dim = dim if dim >= 0 else len(x.shape) + dim

        # Validate placements and sharding
        for i_dim_device_mesh, placement in enumerate(placements_input):
            if isinstance(placement, Partial):
                raise ValueError("Partial placements are not supported")
            elif isinstance(placement, Shard):
                if x.shape[placement.dim] % device_mesh_input.shape[i_dim_device_mesh] != 0:
                    raise ValueError(
                        f"Uneven sharding tensor dimension {placement.dim} of size {x.shape[placement.dim]} "
                        f"along device mesh dimension {i_dim_device_mesh} of size {device_mesh_input.shape[i_dim_device_mesh]} is not supported"
                    )
                if placement.dim == actual_dim:
                    raise ValueError(f"Argmax along sharded dimension {dim} is not supported")

        x_local = x.to_local()
        if keepdim is None:
            output_local = torch.argmax(x_local, dim=dim)
        else:
            output_local = torch.argmax(x_local, dim=dim, keepdim=keepdim)

        shape_output = list(x.shape)
        shape_output[actual_dim] = 1
        strides_output = update_exhaustive_strides(x.shape, x.stride(), shape_output)
        if not keepdim:
            shape_output = shape_output[:actual_dim] + shape_output[actual_dim + 1 :]
            strides_output = strides_output[:actual_dim] + strides_output[actual_dim + 1 :]

        result: DTensor = DTensor.from_local(
            output_local,
            device_mesh=device_mesh_input,
            placements=placements_input,
            shape=tuple(shape_output),
            stride=strides_output,
        )
        ctx.mark_non_differentiable(result)
        return result

    @staticmethod
    def backward(ctx: FunctionCtx, *grad_outputs) -> tuple[None, None, None]:
        return None, None, None


class _ShardwiseOffsetImpl(torch.autograd.Function):
    @staticmethod
    def forward(
        ctx: FunctionCtx,
        x: DTensor,
        dim: int,
        offset_per_rank: Any,
    ) -> DTensor:
        """Forward pass for shardwise offset.

        This function adds a rank-dependent offset to each shard of the input tensor.
        The offset for each shard is computed as: rank_on_mesh_axis * offset_per_rank,
        where rank_on_mesh_axis is the rank of the current process along the device mesh
        axis that shards the specified dimension.

        Parameters
        ----------
        ctx : FunctionCtx
            Context object for saving information needed in backward pass.
        x : DTensor
            Input DTensor with dimension `dim` sharded.
        dim : int
            The dimension that must be sharded.
        offset_per_rank : Any
            The offset value per rank. Can be a scalar or tensor that broadcasts
            with x_local.

        Returns
        -------
        DTensor
            DTensor with offset applied: x + rank * offset_per_rank

        Raises
        ------
        TypeError
            If inputs are not of correct type.
        ValueError
            If the specified dimension is not sharded, partial placements exist,
            or uneven sharding is detected.
        """
        # Type checking
        if not isinstance(x, DTensor):
            raise TypeError(f"Expected DTensor, got {type(x)}")
        if not isinstance(dim, int):
            raise TypeError(f"Expected int for dim, got {type(dim)}")

        device_mesh_input = x.device_mesh
        placements_input = x.placements

        # Normalize negative dim
        actual_dim = dim if dim >= 0 else len(x.shape) + dim

        # Find which device_mesh axis shards the specified dim
        mesh_axis_for_dim = None
        for i_dim_device_mesh, placement in enumerate(placements_input):
            if isinstance(placement, Partial):
                raise ValueError("Partial placements are not supported")
            elif isinstance(placement, Shard):
                # Check that sharded dimensions are evenly divided
                if x.shape[placement.dim] % device_mesh_input.shape[i_dim_device_mesh] != 0:
                    raise ValueError(
                        f"Uneven sharding tensor dimension {placement.dim} of size {x.shape[placement.dim]} "
                        f"along device mesh dimension {i_dim_device_mesh} of size {device_mesh_input.shape[i_dim_device_mesh]} is not supported"
                    )
                if placement.dim == actual_dim:
                    mesh_axis_for_dim = i_dim_device_mesh

        # Check that the specified dimension IS sharded
        if mesh_axis_for_dim is None:
            raise ValueError(f"Dimension {dim} must be sharded for shardwise_offset, but it is not")

        x_local = x.to_local()

        # Get the rank along the mesh axis that shards dim
        rank_on_mesh_axis = device_mesh_input.get_local_rank(mesh_axis_for_dim)

        # Compute offset: x_local + rank * offset_per_rank
        output_local = x_local + rank_on_mesh_axis * offset_per_rank

        # Create output DTensor with same shape and stride as input
        result: DTensor = DTensor.from_local(
            output_local,
            device_mesh=device_mesh_input,
            placements=placements_input,
            shape=x.shape,
            stride=x.stride(),
        )

        return result

    @staticmethod
    def backward(
        ctx: FunctionCtx,
        grad_output: DTensor,
    ) -> tuple[DTensor, None, None]:
        """Backward pass for shardwise offset.

        Since the offset is a constant addition (rank * offset_per_rank is constant
        for each shard), the gradient passes through unchanged.

        Parameters
        ----------
        ctx : FunctionCtx
            Context object containing saved tensors and metadata from forward pass.
        grad_output : DTensor
            Gradient of the loss with respect to the output.

        Returns
        -------
        tuple[DTensor, None, None]
            Gradient with respect to input (same as grad_output), None for dim
            and offset_per_rank parameters.
        """
        # Gradient passes through unchanged since offset is constant
        return grad_output, None, None


def shardwise_offset(x: DTensor, dim: int, offset_per_rank: Any) -> DTensor:
    """Add a rank-dependent offset to a DTensor along a sharded dimension.

    This function adds an offset to each shard based on its rank along the
    device mesh axis that shards the specified dimension. The offset for
    each shard is: rank_on_mesh_axis * offset_per_rank.

    Parameters
    ----------
    x : DTensor
        Input DTensor with dimension `dim` sharded.
    dim : int
        The dimension that must be sharded. The function identifies which
        device mesh axis shards this dimension and uses the rank along that
        axis to compute the offset.
    offset_per_rank : Any
        The offset value per rank. Can be a scalar or tensor that broadcasts
        with x_local. The actual offset added is rank * offset_per_rank.

    Returns
    -------
    DTensor
        DTensor with offset applied: x + rank * offset_per_rank

    Raises
    ------
    TypeError
        If inputs are not of correct type.
    ValueError
        If the specified dimension is not sharded, partial placements exist,
        or uneven sharding is detected.

    Examples
    --------
    >>> # Example: Adding position offsets to sharded sequence indices
    >>> # If dim 1 is sharded across 2 ranks with local size 100:
    >>> # - Rank 0: adds 0 * offset_per_rank to its shard
    >>> # - Rank 1: adds 1 * offset_per_rank to its shard
    >>> indices = ...  # DTensor with dim 1 sharded, local shape (B, 100, D)
    >>> offset_per_rank = 100  # Each rank's shard represents 100 elements
    >>> global_indices = shardwise_offset(indices, dim=1, offset_per_rank=offset_per_rank)
    """
    return _ShardwiseOffsetImpl.apply(x, dim, offset_per_rank)


def shardwise_sum(x: DTensor, dim: int, keepdim: Optional[bool] = None) -> DTensor:
    """Sum elements of a DTensor along a specified dimension.

    This function sums elements of a DTensor along the specified dimension.
    The sum operation is performed on local tensor chunks while maintaining
    gradient computation capabilities.

    Parameters
    ----------
    x : DTensor
        Input DTensor to sum.
    dim : int
        Dimension to sum along.
    keepdim : Optional[bool]
        Whether to keep the dimension when summing.

    Returns
    -------
    DTensor
        DTensor after sum operation.

    Raises
    ------
    TypeError
        If inputs are not of correct type.
    ValueError
        If validation errors occur, such as summing along a sharded dimension.
    """
    return _ShardwiseSumImpl.apply(x, dim, keepdim)


def shardwise_one_hot(input: DTensor, num_classes: int = -1) -> DTensor:
    """One-hot encode a DTensor of class indices.

    This function performs one-hot encoding on a DTensor of class indices.
    The operation is performed on local tensor chunks while maintaining
    the distributed tensor structure. The new one-hot dimension is added
    as the last dimension and is replicated across all devices.

    Parameters
    ----------
    input : DTensor
        Input DTensor containing class indices (integer values).
    num_classes : int
        Number of classes for one-hot encoding. If -1, inferred from input.

    Returns
    -------
    DTensor
        DTensor after one-hot encoding with shape input.shape + (num_classes,).

    Raises
    ------
    TypeError
        If inputs are not of correct type.
    ValueError
        If validation errors occur, such as partial placements or uneven sharding.

    Notes
    -----
    The one-hot operation is not differentiable with respect to the input indices,
    so gradients will be None for the input tensor.
    """
    return _ShardwiseOneHotImpl.apply(input, num_classes)


def shardwise_distogram(d: DTensor, boundaries: torch.Tensor) -> DTensor:
    """Bin pairwise distances into a distogram.

    Computes ``(d.unsqueeze(-1) > boundaries).sum(dim=-1).long()`` element-wise
    on local shards while preserving the DTensor placements and shape.

    Parameters
    ----------
    d : DTensor
        Input DTensor of pairwise distances.
    boundaries : torch.Tensor
        1-D tensor of bin boundaries (regular tensor, not a DTensor).

    Returns
    -------
    DTensor
        Long DTensor of bin indices, same shape and placements as ``d``.

    Notes
    -----
    The distogram binning operation is not differentiable, so gradients
    will be None for both inputs.
    """
    return _ShardwiseDistogramImpl.apply(d, boundaries)


def shardwise_softmax(x: DTensor, dim: int = -1) -> DTensor:
    """Apply softmax to a DTensor along a specified dimension.

    Parameters
    ----------
    x : DTensor
        Input DTensor to apply softmax to.
    dim : int
        Dimension to apply softmax over. Default is -1 (last dimension).

    Returns
    -------
    DTensor
        DTensor after softmax operation.
    """
    return _ShardwiseSoftmaxImpl.apply(x, dim)


def shardwise_log_softmax(x: DTensor, dim: int = -1) -> DTensor:
    """Apply log_softmax to a DTensor along a specified dimension.

    This function applies log_softmax to a DTensor along the specified dimension.
    The log_softmax operation is performed on local tensor chunks while maintaining
    gradient computation capabilities.

    Parameters
    ----------
    x : DTensor
        Input DTensor to apply log_softmax to.
    dim : int
        Dimension to apply log_softmax over. Default is -1 (last dimension).

    Returns
    -------
    DTensor
        DTensor after log_softmax operation.

    Raises
    ------
    TypeError
        If inputs are not of correct type.
    ValueError
        If validation errors occur, such as applying log_softmax along a sharded dimension.

    Notes
    -----
    The log_softmax operation is differentiable and supports gradient computation.
    The operation cannot be applied along sharded dimensions as it requires
    communication across devices to compute the softmax normalization.
    """
    return _ShardwiseLogSoftmaxImpl.apply(x, dim)


def shardwise_argmax(x: DTensor, dim: int, keepdim: Optional[bool] = None) -> DTensor:
    """
    Compute argmax of a DTensor along a specified dimension (per shard).

    This is a shard-local argmax: it does not communicate across shards, so the
    resulting indices correspond to the local shard slices. Use only when the
    reduced dimension is not sharded.

    Parameters
    ----------
    x : DTensor
        Input DTensor.
    dim : int
        Dimension to take argmax over.
    keepdim : Optional[bool]
        Whether to retain the reduced dimension with size 1.

    Returns
    -------
    DTensor
        DTensor of dtype long containing indices of the local argmax.

    Raises
    ------
    TypeError
        If inputs are not of correct type.
    ValueError
        If validation errors occur, such as argmax along a sharded dimension
        or partial placements.
    """
    return _ShardwiseArgmaxImpl.apply(x, dim, keepdim)


class ShardwiseOuterOp(Enum):
    """Supported operations for shardwise outer operations.

    These operations support broadcasting between tensors with singleton
    dimensions for computing pairwise operations efficiently.
    """

    SUBTRACT = "subtract"
    """Element-wise subtraction: x - y. Differentiable."""

    ADD = "add"
    """Element-wise addition: x + y. Differentiable."""

    LOGICAL_AND = "logical_and"
    """Element-wise logical AND: x & y. Non-differentiable."""

    EQUAL = "equal"
    """Element-wise equality: x == y. Non-differentiable."""


class _ShardwiseOuterOpImpl(torch.autograd.Function):
    """Unified shardwise outer operation at a specified axis.

    This autograd function handles all outer operations (subtract, addition, logical_and, equal)
    with a shared code path, differing only in the actual math operation performed.
    Supports gradient computation for differentiable operations (SUBTRACT, ADD).

    The operation computes pairwise combinations at the specified axis:
    - x: (..., L, ...) at axis
    - y: (..., R, ...) at axis
    - Result: (..., L, R, ...) with one additional dimension
    """

    @staticmethod
    def forward(
        ctx: FunctionCtx,
        x: DTensor,
        y: DTensor,
        op: ShardwiseOuterOp,
        axis: int,
    ) -> DTensor:
        """Forward pass for outer operations at specified axis.

        Parameters
        ----------
        ctx : FunctionCtx
            Context object for saving information needed in backward pass.
        x : DTensor
            First input tensor with shape (..., L, ...) where L is at position `axis`.
        y : DTensor
            Second input tensor with shape (..., R, ...) where R is at position `axis`.
        op : ShardwiseOuterOp
            The operation to perform (SUBTRACT, ADD, LOGICAL_AND, or EQUAL).
        axis : int
            The axis at which to perform the outer operation.

        Returns
        -------
        DTensor
            Result of the operation with shape (..., L, R, ...).
            The output has one more dimension than the inputs.

        Raises
        ------
        TypeError
            If inputs are not DTensors or axis is not an int.
        ValueError
            If device_mesh, placements don't match, or an unsupported operation is specified.
        """
        # ========== Type checking ==========
        if not isinstance(x, DTensor):
            raise TypeError(f"shardwise_outer_op: Expected DTensor for x, got {type(x)}")
        if not isinstance(y, DTensor):
            raise TypeError(f"shardwise_outer_op: Expected DTensor for y, got {type(y)}")
        if not isinstance(op, ShardwiseOuterOp):
            raise TypeError(f"shardwise_outer_op: Expected ShardwiseOuterOp for op, got {type(op)}")
        if not isinstance(axis, int):
            raise TypeError(f"shardwise_outer_op: Expected int for axis, got {type(axis)}")

        # ========== Validate device_mesh and placements match ==========
        if x.device_mesh != y.device_mesh:
            raise ValueError("shardwise_outer_op: x and y must have the same device_mesh")
        if x.placements != y.placements:
            raise ValueError("shardwise_outer_op: x and y must have the same placements")

        device_mesh = x.device_mesh
        placements_input = x.placements

        # ========== Validate placements ==========
        for i_dim_device_mesh, placement in enumerate(placements_input):
            if isinstance(placement, Partial):
                raise ValueError("shardwise_outer_op: Partial placements are not supported")
            elif isinstance(placement, Shard):
                # The outer operation axis must not be sharded - this is a shardwise op
                # so the outer product must be computed locally without cross-shard communication
                if placement.dim == axis:
                    raise ValueError(
                        f"shardwise_outer_op: Cannot shard dimension {axis} (the outer operation axis) "
                        f"with Shard({placement.dim}). The outer operation must be computed locally "
                        f"on each shard without cross-shard communication."
                    )

                x_dim_size = x.shape[placement.dim]
                y_dim_size = y.shape[placement.dim]
                mesh_dim_size = device_mesh.shape[i_dim_device_mesh]

                if x_dim_size % mesh_dim_size != 0:
                    raise ValueError(
                        f"shardwise_outer_op: Uneven sharding of x tensor dimension {placement.dim} "
                        f"of size {x_dim_size} along device mesh dimension "
                        f"{i_dim_device_mesh} of size {mesh_dim_size} is not supported"
                    )
                if y_dim_size % mesh_dim_size != 0:
                    raise ValueError(
                        f"shardwise_outer_op: Uneven sharding of y tensor dimension {placement.dim} "
                        f"of size {y_dim_size} along device mesh dimension "
                        f"{i_dim_device_mesh} of size {mesh_dim_size} is not supported"
                    )

        # ========== Get local tensors and unsqueeze ==========
        x_local = x.to_local()
        y_local = y.to_local()

        # Unsqueeze local tensors to create broadcast-compatible shapes
        # x: (..., L, ...) → (..., L, 1, ...)  (insert singleton after axis)
        # y: (..., R, ...) → (..., 1, R, ...)  (insert singleton at axis)
        x_local = x_local.unsqueeze(axis + 1)
        y_local = y_local.unsqueeze(axis)

        # Compute output shape: (..., L, R, ...)
        shape_output = list(x.shape)
        shape_output.insert(axis + 1, y.shape[axis])
        shape_output = tuple(shape_output)

        # Adjust placements for the new dimension
        # Shard dimensions > axis need to be incremented
        placements_output = list(placements_input)
        for i_dim_device_mesh, p in enumerate(placements_input):
            if isinstance(p, Shard) and p.dim > axis:
                placements_output[i_dim_device_mesh] = Shard(p.dim + 1)
        placements_output = tuple(placements_output)

        # ========== Compute output stride ==========
        layout_right = LayoutRightMap(shape_output)
        stride_output = layout_right.strides

        # ========== Perform operation based on op type ==========
        if op == ShardwiseOuterOp.SUBTRACT:
            output_local = x_local - y_local
        elif op == ShardwiseOuterOp.ADD:
            output_local = x_local + y_local
        elif op == ShardwiseOuterOp.LOGICAL_AND:
            output_local = x_local & y_local
        elif op == ShardwiseOuterOp.EQUAL:
            output_local = x_local == y_local
        else:
            raise ValueError(f"shardwise_outer_op: Unsupported operation: {op}")

        # ========== Create output DTensor ==========
        result = DTensor.from_local(
            output_local,
            device_mesh=device_mesh,
            placements=placements_output,
            shape=shape_output,
            stride=stride_output,
        )

        # ========== Handle gradient context ==========
        is_differentiable = op in (ShardwiseOuterOp.SUBTRACT, ShardwiseOuterOp.ADD)
        if is_differentiable and (x.requires_grad or y.requires_grad):
            # Differentiable ops (SUBTRACT, ADD) - save context for backward
            ctx.op = op
            ctx.axis = axis
            ctx.device_mesh = device_mesh
            ctx.placements_input = placements_input
            ctx.placements_output = placements_output
            ctx.x_shape = x.shape
            ctx.y_shape = y.shape
            ctx.x_stride = x.stride()
            ctx.y_stride = y.stride()
            ctx.output_shape = shape_output
            ctx.x_requires_grad = x.requires_grad
            ctx.y_requires_grad = y.requires_grad
        else:
            # Non-differentiable operations or no grad required
            ctx.mark_non_differentiable(result)
            ctx.op = op

        return result

    @staticmethod
    def backward(ctx: FunctionCtx, grad_output: DTensor) -> tuple[Optional[DTensor], Optional[DTensor], None, None]:
        """Backward pass for outer operations.

        Differentiable operations:
        - SUBTRACT: For z = x - y:
            grad_x = grad_z (summed over axis+1 where x was broadcast)
            grad_y = -grad_z (summed over axis where y was broadcast)
        - ADD: For z = x + y:
            grad_x = grad_z (summed over axis+1 where x was broadcast)
            grad_y = grad_z (summed over axis where y was broadcast)

        Parameters
        ----------
        ctx : FunctionCtx
            Context object with saved information from forward pass.
        grad_output : DTensor
            Gradient with respect to the output.

        Returns
        -------
        tuple[Optional[DTensor], Optional[DTensor], None, None]
            Gradients for x, y, None for op parameter, and None for axis parameter.
        """
        grad_x = None
        grad_y = None

        # ========== Op-specific gradient computation ==========
        if ctx.op == ShardwiseOuterOp.SUBTRACT:
            # For z = x - y: grad_x = grad_z, grad_y = -grad_z
            # Both need to be summed over their broadcast dimensions
            raise_if_incorrect_dtensor_metadata_args(
                dtensor_instance=grad_output,
                dtensor_name="_ShardwiseOuterOpImpl.backward grad_output",
                expected_shape=ctx.output_shape,
                expected_device_mesh=ctx.device_mesh,
                expected_placements=ctx.placements_output,
            )

            grad_output_local = grad_output.to_local()

            if ctx.x_requires_grad:
                # Sum over axis+1 (where x was broadcast) and squeeze
                grad_x_local = grad_output_local.sum(dim=ctx.axis + 1, keepdim=False)
                grad_x = DTensor.from_local(
                    grad_x_local,
                    device_mesh=ctx.device_mesh,
                    placements=ctx.placements_input,
                    shape=ctx.x_shape,
                    stride=ctx.x_stride,
                )

            if ctx.y_requires_grad:
                # Sum over axis (where y was broadcast), negate, and squeeze
                grad_y_local = -grad_output_local.sum(dim=ctx.axis, keepdim=False)
                grad_y = DTensor.from_local(
                    grad_y_local,
                    device_mesh=ctx.device_mesh,
                    placements=ctx.placements_input,
                    shape=ctx.y_shape,
                    stride=ctx.y_stride,
                )

        elif ctx.op == ShardwiseOuterOp.ADD:
            # For z = x + y: grad_x = grad_z, grad_y = grad_z
            # Both need to be summed over their broadcast dimensions
            raise_if_incorrect_dtensor_metadata_args(
                dtensor_instance=grad_output,
                dtensor_name="_ShardwiseOuterOpImpl.backward grad_output",
                expected_shape=ctx.output_shape,
                expected_device_mesh=ctx.device_mesh,
                expected_placements=ctx.placements_output,
            )

            grad_output_local = grad_output.to_local()

            if ctx.x_requires_grad:
                # Sum over axis+1 (where x was broadcast) and squeeze
                grad_x_local = grad_output_local.sum(dim=ctx.axis + 1, keepdim=False)
                grad_x = DTensor.from_local(
                    grad_x_local,
                    device_mesh=ctx.device_mesh,
                    placements=ctx.placements_input,
                    shape=ctx.x_shape,
                    stride=ctx.x_stride,
                )

            if ctx.y_requires_grad:
                # Sum over axis (where y was broadcast) and squeeze (no negation for addition)
                grad_y_local = grad_output_local.sum(dim=ctx.axis, keepdim=False)
                grad_y = DTensor.from_local(
                    grad_y_local,
                    device_mesh=ctx.device_mesh,
                    placements=ctx.placements_input,
                    shape=ctx.y_shape,
                    stride=ctx.y_stride,
                )

        elif ctx.op == ShardwiseOuterOp.LOGICAL_AND:
            # Non-differentiable
            pass

        elif ctx.op == ShardwiseOuterOp.EQUAL:
            # Non-differentiable
            pass

        else:
            raise ValueError(f"shardwise_outer_op backward: Unsupported operation: {ctx.op}")

        return grad_x, grad_y, None, None


def shardwise_outer_op(lhs: DTensor, rhs: DTensor, axis: int, op: ShardwiseOuterOp) -> DTensor:
    """Compute outer operation at specified axis.

    This function performs an outer operation (subtract, logical_and, equal) between
    two tensors at a specified axis. The operation creates pairwise combinations
    of elements along the specified axis.

    The function internally unsqueezes the inputs to create broadcast-compatible
    shapes and then performs the operation:
    - lhs: (..., L, ...) → unsqueeze to (..., L, 1, ...)
    - rhs: (..., R, ...) → unsqueeze to (..., 1, R, ...)
    - Result: (..., L, R, ...)

    Parameters
    ----------
    lhs : DTensor
        First input tensor with shape (..., L, ...) where L is at position `axis`.
    rhs : DTensor
        Second input tensor with shape (..., R, ...) where R is at position `axis`.
    axis : int
        The axis at which to perform the outer operation. Must be a valid
        dimension index for both tensors.
    op : ShardwiseOuterOp
        The operation to perform:
        - SUBTRACT: lhs - rhs (differentiable)
        - LOGICAL_AND: lhs & rhs (non-differentiable, boolean inputs)
        - EQUAL: lhs == rhs (non-differentiable)

    Returns
    -------
    DTensor
        Result tensor with shape (..., L, R, ...) where:
        - L (from lhs) is at position `axis`
        - R (from rhs) is at position `axis + 1`
        - All other dimensions match the input dimensions.

    Raises
    ------
    TypeError
        If inputs are not DTensors or op is not a ShardwiseOuterOp.
    ValueError
        If device_mesh, placements don't match, axis is invalid, or shapes
        are incompatible for the outer operation.

    Examples
    --------
    >>> # Window batching example
    >>> lhs = ...  # DTensor with shape (B, K, W, D)
    >>> rhs = ...  # DTensor with shape (B, K, H, D)
    >>> result = shardwise_outer_op(lhs, rhs, axis=2, op=ShardwiseOuterOp.SUBTRACT)
    >>> # result has shape (B, K, W, H, D)

    >>> # Computing pairwise differences
    >>> queries = ...  # DTensor with shape (B, N, D)
    >>> keys = ...     # DTensor with shape (B, M, D)
    >>> diffs = shardwise_outer_op(queries, keys, axis=1, op=ShardwiseOuterOp.SUBTRACT)
    >>> # diffs has shape (B, N, M, D)

    Notes
    -----
    This function provides a cleaner API compared to manually unsqueezing
    tensors before calling the operation. It is the preferred way to compute
    outer operations when the input tensors don't already have singleton
    dimensions for broadcasting.
    """
    # Validate axis bounds
    if not isinstance(axis, int):
        raise TypeError(f"shardwise_outer_op: Expected int for axis, got {type(axis)}")

    if not isinstance(lhs, DTensor):
        raise TypeError(f"shardwise_outer_op: Expected DTensor for lhs, got {type(lhs)}")
    if not isinstance(rhs, DTensor):
        raise TypeError(f"shardwise_outer_op: Expected DTensor for rhs, got {type(rhs)}")

    # Normalize negative axis
    ndim = lhs.ndim
    if axis < 0:
        axis = ndim + axis

    if axis < 0 or axis >= ndim:
        raise ValueError(f"shardwise_outer_op: axis {axis} is out of bounds for tensor with {ndim} dimensions")

    if rhs.ndim != ndim:
        raise ValueError(
            f"shardwise_outer_op: lhs and rhs must have the same number of dimensions, "
            f"got lhs.ndim={lhs.ndim} and rhs.ndim={rhs.ndim}"
        )

    # Pass axis to the autograd function - unsqueezing happens on local tensors
    return _ShardwiseOuterOpImpl.apply(lhs, rhs, op, axis)
