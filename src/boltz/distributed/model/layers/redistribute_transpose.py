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
from torch import Tensor
from torch.autograd.function import FunctionCtx
from torch.distributed.tensor import DTensor, Partial, Placement, Shard

from boltz.distributed.comm import TransposeComm
from boltz.distributed.model.layers.dtensor_metadata_tools import (
    raise_if_incorrect_dtensor_metadata_args,
)


def redistribute_transpose(
    input: DTensor,
    transpose_comm: Optional[TransposeComm],
    output_placements: Optional[tuple[Placement, ...]],
    dim0: Optional[int] = None,
    dim1: Optional[int] = None,
) -> DTensor:
    """Transpose a DTensor across device mesh (and locally).

    Use cases in Boltz:
        (1) boltz.model.modules.trunk.py: DistogramModule.forward [impl'd]
            - redistribute_transpose(z, self.distogram_comm, (Shard(0), Shard(1), Shard(2)), 1, 2),

        (2) boltz.model.modules.trunk.py: MSAModule.forward [impl'd]
            - redistribute_transpose(emb, self.comm_transpose, (Shard(0), Shard(1), Shard(2)), 1, 2)

        (3) boltz.model.loss.distogram.py: distogram_loss
            - redistribute_transpose(mask, comm, (Shard(0), Shard(1), Shard(2)), 1, 2)

        (4) boltz.model.modules.encoders.py: AtomAttentionEncoder.forward
            - redistribute_transpose(c, self.transpose_comm_c, (Shard(0), Replicate(), Shard(1)), None, None)

        (5) boltz.model.model.py: Boltz1.forward
            - redistribute_transpose(s_inputs, self.transpose_comm, (Shard(0), Replicate(), Shard(1)), None, None)

    Parameters
    ----------
    input : DTensor
        Input tensor to transpose
    transpose_comm : Optional[TransposeComm]
        Communication object for distributed operations
    output_placements : Optional[tuple[Placement, ...]]
        Output placements for the DTensor.
    dim0 : Optional[int]
        First dimension to transpose locally.
    dim1 : Optional[int]
        Second dimension to transpose locally.

    Returns
    -------
    DTensor
        DTensor transposed across device mesh (and locally).
    """
    return _RedistributeTransposeImpl.apply(input, output_placements, transpose_comm, dim0, dim1)


class _RedistributeTransposeImpl(torch.autograd.Function):
    """Custom autograd function to transpose a DTensor across device mesh (and locally)."""

    @staticmethod
    def forward(
        ctx: FunctionCtx,
        input: DTensor,
        output_placements: Optional[tuple[Placement, ...]],
        transpose_comm: Optional[TransposeComm],
        dim0: Optional[int] = None,
        dim1: Optional[int] = None,
    ) -> DTensor:
        """Forward pass for _RedistributeTransposeImpl custom autograd function.

        Parameters
        ----------
        ctx : FunctionCtx
            Context object to save information for backward pass
        input : DTensor
            Input tensor to transpose and redistribute
        output_placements : tuple[Placement, ...]
            Output placements for the DTensor.
        transpose_comm : Union[TransposeComm, None]
            Communication object for distributed operations
        dim0 : int
            First tensor dimension to transpose locally.
        dim1 : int
            Second tensor dimension to transpose locally.

        Returns
        -------
        DTensor
            DTensor transposed across device mesh (and locally).
        """
        # check input options
        if (dim0 is None) != (dim1 is None):
            raise ValueError(
                " When using redistribute_transpose, either both dim0 and dim1 must be None if no local transposition, or both must be not None for local transposition"
            )

        ctx.is_local_transpose = dim1 is not None
        ctx.is_device_mesh_transpose = transpose_comm is not None

        # Short circuit if no local or device mesh transpose is performed
        if not ctx.is_local_transpose and not ctx.is_device_mesh_transpose:
            return input

        if (transpose_comm is None) != (output_placements is None):
            raise ValueError(
                "transpose_comm and output_placements must be either both None or both not None for device mesh transpose"
            )

        if ctx.is_local_transpose and ctx.is_device_mesh_transpose and (output_placements != input.placements):
            raise ValueError(
                "Simultaneous redistribute and local transpose is only supported when the two transposing axes are the sharding axes involved in the said redistribute. For other usage cases, consider decompose the operation in a redistribute-only followed by a local transpose-only operations by calling redistribute_transpose() twice with different arguments"
            )

        axis_mesh_shard_dim0 = None
        axis_mesh_shard_dim1 = None

        if ctx.is_device_mesh_transpose:
            axes_mesh_transpose = []
        else:
            axes_mesh_transpose = None

        for i_dim_device_mesh, placement in enumerate(input.placements):
            # Check if partial placements
            if isinstance(placement, Partial):
                raise ValueError(
                    f"Partial placements are not supported for redistribute_transpose but {input.placements} is given"
                )

            # Check if sharding is even
            if (
                isinstance(placement, Shard)
                and input.shape[placement.dim] % input.device_mesh.shape[i_dim_device_mesh] != 0
            ):
                raise ValueError(
                    f"Uneven sharding tensor dimension {placement.dim} of size {input.shape[placement.dim]} "
                    f"along device mesh dimension {i_dim_device_mesh} of size {input.device_mesh.shape[i_dim_device_mesh]} is not supported"
                )

            if ctx.is_local_transpose and isinstance(placement, Shard):
                if placement.dim == dim0:
                    axis_mesh_shard_dim0 = i_dim_device_mesh
                if placement.dim == dim1:
                    axis_mesh_shard_dim1 = i_dim_device_mesh

            if axes_mesh_transpose is not None:
                if placement != output_placements[i_dim_device_mesh]:
                    axes_mesh_transpose.append(i_dim_device_mesh)

        # Check if locally transposed dimensions are sharded if both local and device mesh transposes are performed
        if (ctx.is_local_transpose and ctx.is_device_mesh_transpose) and (
            axis_mesh_shard_dim0 is None or axis_mesh_shard_dim1 is None
        ):
            raise ValueError(
                f"Both dim0 and dim1 must be sharded when doing both local and device mesh transposes "
                f"but dim0={dim0} and dim1={dim1} are given with placements={input.placements}"
            )

        # Check if locally transposed dimensions are sharded if only local transpose is performed
        if (ctx.is_local_transpose and not ctx.is_device_mesh_transpose) and (
            axis_mesh_shard_dim0 is not None or axis_mesh_shard_dim1 is not None
        ):
            raise NotImplementedError(
                "Local transpose on sharded dimensions is not supported when only local transpose is performed"
            )

        if ctx.is_device_mesh_transpose:
            device_mesh_coords = input.device_mesh.get_coordinate()
            if len(axes_mesh_transpose) == 0:
                # this implies dim{0, 1} are sharded and output_placements == input_placements
                # but the underlying device mesh transpose will be performed by the transpose_comm
                # along the two Sharding placement axes
                assert axis_mesh_shard_dim0 is not None and axis_mesh_shard_dim1 is not None
                axes_mesh_transpose = [axis_mesh_shard_dim0, axis_mesh_shard_dim1]
            else:
                # assert output placements is strictly a permutation of input placements
                if not (
                    input.placements[axes_mesh_transpose[0]] == output_placements[axes_mesh_transpose[1]]
                    and input.placements[axes_mesh_transpose[1]] == output_placements[axes_mesh_transpose[0]]
                ):
                    raise ValueError(
                        "Input and output placements are not strictly a permutation of each other along mesh transpose axes:"
                        f"input.placements={input.placements} vs. output_placements={output_placements}"
                    )

            # assert the correspondence of transpose_comm's underlying group to the device mesh axes
            if (
                device_mesh_coords[axes_mesh_transpose[0]],
                device_mesh_coords[axes_mesh_transpose[1]],
            ) != transpose_comm.rank_coords:
                raise ValueError(
                    f"Inconsistent device mesh coordinate {device_mesh_coords} along mesh transpose axes {axes_mesh_transpose} "
                    f"compared to transpose_comm rank_coords {transpose_comm.rank_coords}"
                )

        # transpose input tensor
        output_local: Tensor = input.to_local()

        if ctx.is_device_mesh_transpose:
            output_local_ = transpose_comm.enqueue_to_dispatch(output_local.contiguous())
            transpose_comm.wait_until_finished()
            output_local = output_local_
        if ctx.is_local_transpose:
            output_local = output_local.transpose(dim0, dim1)

        output_shape = torch.Size(
            _swap_tuple_elements(input.shape, dim0, dim1) if ctx.is_local_transpose else input.shape
        )
        output_stride = _swap_tuple_elements(input.stride(), dim0, dim1) if ctx.is_local_transpose else input.stride()

        if input.requires_grad:
            ctx.input_shape = input.shape
            ctx.output_shape = output_shape
            ctx.input_stride = input.stride()
            ctx.output_stride = output_stride
            ctx.input_placements = input.placements
            ctx.output_placements = input.placements if output_placements is None else output_placements
            ctx.device_mesh = input.device_mesh
            ctx.dim0 = dim0
            ctx.dim1 = dim1
            ctx.transpose_comm = transpose_comm

        # Create a new DTensor called output
        output: DTensor = DTensor.from_local(
            local_tensor=output_local,
            shape=output_shape,
            stride=output_stride,
            device_mesh=input.device_mesh,
            placements=input.placements if output_placements is None else output_placements,
        )

        return output

    @staticmethod
    def backward(
        ctx: FunctionCtx,
        grad_output: DTensor,
    ) -> tuple[DTensor, None, None, None, None]:
        """Backward pass for _RedistributeTranspose custom autograd.Function

        Parameters
        ----------
        ctx : FunctionCtx
            Context object with saved information from forward pass
        grad_output : DTensor
            Gradient tensor from downstream layers

        Returns
        -------
        tuple[DTensor, None, None, None]
            Tuple containing the gradient for input tensor and None for other parameters
        """
        # Short circuit if no local or device mesh transpose is performed
        if not ctx.is_local_transpose and not ctx.is_device_mesh_transpose:
            return grad_output, None, None, None, None

        # metadata check on grad_output
        raise_if_incorrect_dtensor_metadata_args(
            dtensor_instance=grad_output,
            dtensor_name="grad_output",
            expected_shape=ctx.output_shape,
            expected_device_mesh=ctx.device_mesh,
            expected_placements=ctx.output_placements,
            check_for_partial_placements=False,
        )

        # transpose gradient tensor
        grad_input_local = grad_output.to_local()

        if ctx.is_device_mesh_transpose:
            grad_input_local_ = ctx.transpose_comm.enqueue_to_dispatch(grad_input_local.contiguous())
            ctx.transpose_comm.wait_until_finished()
            grad_input_local = grad_input_local_
        if ctx.is_local_transpose:
            grad_input_local = grad_input_local.transpose(ctx.dim0, ctx.dim1)

        # Create a new DTensor called output
        grad_input: DTensor = DTensor.from_local(
            grad_input_local,
            shape=ctx.input_shape,
            stride=ctx.input_stride,
            device_mesh=ctx.device_mesh,
            placements=ctx.input_placements,
        )
        return grad_input, None, None, None, None


def _swap_tuple_elements(x: tuple[int, ...], i: int, j: int) -> tuple[int, ...]:
    """Swap two elements of a tuple.

    Parameters
    ----------
    x : tuple[int, ...]
        Tuple to swap elements of
    i : int
        Index of first element to swap
    j : int
    """
    y = list(x)
    y[i], y[j] = y[j], y[i]
    return tuple(y)
