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
from torch.distributed.tensor import DTensor, Replicate, Shard


class _ApplyDropoutMaskMsaOrPairImpl(torch.autograd.Function):
    """Distributed implementation of apply_dropout_mask_msa_or_pair using DTensor."""

    @staticmethod
    def forward(
        ctx, src: DTensor, dropout: float, training: bool, columnwise: bool | None, samples_dropout: DTensor | None
    ) -> DTensor:
        """Forward pass of distributed dropout samples_dropout application.

        Parameters
        ----------
        ctx : torch.autograd.function.BackwardCFrame
            Context object for saving information needed in backward pass.
        src : DTensor
            The source tensor to apply dropout to.
        dropout : float
            The dropout rate between 0.0 and 1.0.
        training : bool
            Whether the model is in training mode.
        columnwise : bool, optional
            If True, applies the same samples_dropout to all elements in each column, by default False.
        samples_dropout : DTensor, optional
            These are the uniform random numbers drawn from [0, 1) and samples_dropout >= dropout will be
            used as the dropout mask. If None, a new samples_dropout is created. Note that currently
            there is no effective way to generate consistent random number sequences between the
            serial and distributed versions so this argument can be used to passed in pre-generated
            samples from the serial version sliced according to the input "src" placements as a mock
            to reproduce the random number sequence in the distributed version. We use this method
            in the tests to verify the distributed version is consistent with the serial version.


        Returns
        -------
        DTensor
            The source tensor with dropout applied during training, or unchanged during inference.
        """
        # Check if inputs are of type DTensor
        if not isinstance(src, DTensor):
            raise TypeError(f"Input 'src' must be of type DTensor. Got type {type(src)}.")

        # Verify that src is 4-dimensional as required for indexing patterns
        if src.ndim != 4:
            raise ValueError(f"Input tensor 'src' must be 4-dimensional. Got {src.ndim} dimensions.")

        if samples_dropout is not None:
            if not isinstance(samples_dropout, DTensor):
                raise TypeError(f"Input 'samples_dropout' must be of type DTensor. Got type {type(samples_dropout)}.")

            if samples_dropout.ndim != 4:
                raise ValueError(
                    f"Input tensor 'samples_dropout' must be 4-dimensional. Got {samples_dropout.ndim} dimensions."
                )

            if samples_dropout.device_mesh != src.device_mesh:
                raise ValueError(
                    f"Input tensor 'samples_dropout' must have the same device mesh as the input tensor. "
                    f"Got device meshes {samples_dropout.device_mesh} and {src.device_mesh}."
                )

            if samples_dropout.requires_grad:
                raise ValueError(
                    "Input tensor 'samples_dropout' must not require gradients and its gradient computation is not supported"
                )

            if columnwise:
                if samples_dropout.shape[1] != 1 or samples_dropout.shape[3] != 1:
                    raise ValueError(
                        f"Input tensor 'samples_dropout' must have shape [*, 1, *, 1] for columnwise dropout. Got {samples_dropout.shape}."
                    )
                if samples_dropout.placements != (Shard(0), Replicate(), Shard(2)):
                    raise ValueError(
                        f"Input tensor 'samples_dropout' must have placements (Shard(0), Replicate(), Shard(2)) for columnwise dropout. Got {samples_dropout.placements}."
                    )
            else:
                if samples_dropout.shape[2] != 1 or samples_dropout.shape[3] != 1:
                    raise ValueError(
                        f"Input tensor 'samples_dropout' must have shape [*, *, 1, 1] for rowwise dropout. Got {samples_dropout.shape}."
                    )
                if samples_dropout.placements != (Shard(0), Shard(1), Replicate()):
                    raise ValueError(
                        f"Input tensor 'samples_dropout' must have placements (Shard(0), Shard(1), Replicate()) for rowwise dropout. Got {samples_dropout.placements}."
                    )

            ctx.mark_non_differentiable(samples_dropout)

        # Verify that src.placements is exactly (Shard(0), Shard(1), Shard(2))
        expected_placements = (Shard(0), Shard(1), Shard(2))
        if src.placements != expected_placements:
            raise ValueError(
                f"Input tensor 'src' must have placements {expected_placements}. Got placements {src.placements}."
            )

        # Save context for backward pass
        ctx.device_mesh_input = src.device_mesh
        ctx.placements_input = src.placements
        ctx.training = training
        ctx.columnwise = columnwise
        ctx.dropout = dropout

        # Extract local tensors
        src_local = src.to_local()

        if training:
            if samples_dropout is None:
                # Create dropout samples_dropout using the same logic as the serial version
                shape = list(src_local.shape)
                if columnwise:
                    # equivalent to torch.rand_like(src_local[:, 0, :, 0]
                    shape[1] = 1
                    shape[3] = 1
                else:
                    # equivalent to torch.rand_like(src_local[:, :, 0, 0])
                    shape[2] = 1
                    shape[3] = 1
                if torch.is_autocast_enabled("cuda"):
                    mask_dtype = torch.promote_types(src_local.dtype, torch.float32)
                else:
                    mask_dtype = src_local.dtype
                samples_dropout_local = torch.rand(shape, device=src_local.device, dtype=mask_dtype)
            else:
                samples_dropout_local = samples_dropout.to_local()
            d = samples_dropout_local >= dropout
            if torch.is_autocast_enabled("cuda"):
                scale_dtype = torch.promote_types(src_local.dtype, torch.float32)
            else:
                scale_dtype = src_local.dtype
            d = (d * 1.0 / (1.0 - dropout)).to(dtype=scale_dtype)

            # Apply dropout mask
            result_local = src_local * d

            # Save dropout mask for backward pass
            ctx.save_for_backward(d)
            # Convert result back to DTensor
            result = DTensor.from_local(
                result_local,
                device_mesh=ctx.device_mesh_input,
                placements=ctx.placements_input,
                shape=src.shape,
                stride=src.stride(),
            )

            return result
        else:
            return src

    @staticmethod
    def backward(ctx, grad_output: DTensor) -> tuple[DTensor, None, None, None, None]:
        """Backward pass of distributed dropout mask application.

        Parameters
        ----------
        ctx : torch.autograd.function.BackwardCFrame
            Context object containing saved information from forward pass.
        grad_output : DTensor
            Gradient of the loss with respect to the output.

        Returns
        -------
        tuple[DTensor, None, None, None, None]
            Gradients with respect to inputs. Only src gets a gradient.
        """
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

        # Extract local gradient
        grad_output_local = grad_output.to_local()

        if ctx.training:
            # Extract saved dropout mask
            (dropout_mask,) = ctx.saved_tensors

            # Apply the same dropout mask to the gradient
            grad_src_local = grad_output_local * dropout_mask
        else:
            # During inference, gradient passes through unchanged
            grad_src_local = grad_output_local

        # Convert gradient back to DTensor
        grad_src = DTensor.from_local(
            grad_src_local,
            device_mesh=ctx.device_mesh_input,
            placements=ctx.placements_input,
            shape=grad_output.shape,
            stride=grad_output.stride(),
        )

        # Return gradients: (src, dropout, training, columnwise, mask)
        # Only src needs gradient, others are None
        return grad_src, None, None, None, None


def apply_dropout_mask_msa_or_pair(
    src: DTensor,
    dropout: float,
    training: bool,
    columnwise: bool = False,
    samples_dropout: DTensor | None = None,
) -> DTensor:
    """Apply dropout directly to the source DTensor for MSA or pair representations.

    This function applies dropout to the source DTensor using the shape of z as a reference.
    It behaves like standard dropout during training, and is a no-op during inference.

    When columnwise=True, the same dropout mask is applied to all elements in the same column,
    meaning that entire columns are either kept or dropped together.

    IMPORTANT: This function makes strong assumptions about tensor shape and indexing.
    The reference tensor z must be indexable by [:, 0:1, :, 0:1] (columnwise=True) or
    [:, :, 0:1, 0:1] (columnwise=False). This is specifically designed for MSA and pair
    representation tensors with expected 4D structure.

    Parameters
    ----------
    src : DTensor
        The source DTensor to apply dropout to. Must have placements (Shard(0), Shard(1), Shard(2)).
    dropout : float
        The dropout rate between 0.0 and 1.0
    training : bool
        Whether the model is in training mode
    columnwise : bool, optional
        If True, applies the same mask to all elements in each column, by default False
    samples_dropout : DTensor, optional
        These are the uniform random numbers drawn from [0, 1) and samples_dropout >= dropout will be
        used as the dropout mask. If None, a new samples_dropout is created. Note that currently
        there is no effective way to generate consistent random number sequences between the
        serial and distributed versions so this argument can be used to passed in pre-generated
        samples from the serial version sliced according to the input "src" placements as a mock
        to reproduce the random number sequence in the distributed version. We use this method
        in the tests to verify the distributed version is consistent with the serial version.

    Returns
    -------
    DTensor
        The source DTensor with dropout applied during training, or unchanged during inference

    Notes
    -----
    During training, the values that are kept are scaled by 1/(1-dropout) to maintain
    the expected value of the tensor. During inference (training=False), the input tensor
    is returned unchanged.

    The implementation uses a custom autograd function to handle DTensor operations
    by working with local tensors and properly managing distributed tensor metadata.

    This function enforces specific placement requirements and tensor shape assumptions,
    making it suitable only for MSA and pair representation tensors in the expected format.
    """
    if not training:
        return src
    return _ApplyDropoutMaskMsaOrPairImpl.apply(src, dropout, training, columnwise, samples_dropout)
