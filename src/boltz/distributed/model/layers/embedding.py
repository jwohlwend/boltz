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


from typing import Optional, cast

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.distributed.device_mesh import DeviceMesh
from torch.distributed.tensor import DTensor, Partial, Replicate, Shard, distribute_tensor

from boltz.distributed.utils import update_exhaustive_strides


class _EmbeddingParamsReplicatedImpl(torch.autograd.Function):
    """
    Custom autograd implementation for embedding with replicated parameters.

    The embedding weight is replicated across all device mesh dimensions, while the input
    indices can be sharded or replicated. The output DTensor follows the same placements
    as the input indices.
    """

    @staticmethod
    @torch.amp.custom_fwd(device_type="cuda")
    def forward(
        ctx,
        x: DTensor,
        weight: DTensor,
        padding_idx: Optional[int],
    ) -> DTensor:
        """
        Forward pass for the distributed embedding operation.

        Assumptions and requirements:
        1. Parameters (weight) must be replicated on all device mesh dimensions
        2. Input tensor and parameters must be on the same device mesh
        3. Partial reduction along any input dimension is not supported
        4. Input indices must be integer dtype and must not require gradients

        Args:
            ctx: Context object to store information for backward pass
            x: Input index DTensor with arbitrary placement strategy, except Partial placement
            weight: Weight DTensor with all-replicate placements
            padding_idx: Optional padding index to be passed into F.embedding

        Returns:
            Output DTensor with same placement strategy as input

        Raises:
            ValueError: If any of the placement requirements are violated
        """
        if not isinstance(x, DTensor):
            raise TypeError(f"Expected x to be a DTensor but got {type(x)}.")
        if x.dtype.is_floating_point or x.dtype.is_complex or x.dtype == torch.bool:
            raise ValueError(f"Expected x to be an integer DTensor but got {x.dtype}.")
        if x.requires_grad:
            raise ValueError("x must not require grad in the forward pass")
        if not isinstance(weight, DTensor):
            raise TypeError(f"Expected weight to be a DTensor but got {type(weight)}.")

        device_mesh = x.device_mesh
        if weight.device_mesh != device_mesh:
            raise ValueError("weight and x must be on the same device mesh")

        ndim_device_mesh = device_mesh.ndim
        all_replicate_placements = tuple([Replicate()] * ndim_device_mesh)
        if weight.placements != all_replicate_placements:
            raise ValueError("weight must be replicated on all device mesh dimensions")

        placements_grad_params = list(weight.placements)

        for i_dim_device_mesh, placement in enumerate(x.placements):
            if isinstance(placement, Partial):
                raise ValueError("Partial reduction along any input dimension is not supported")
            if isinstance(placement, Shard):
                if placement.dim >= x.ndim:
                    raise ValueError("Input placement sharding dimension is out of range")
                if x.shape[placement.dim] % device_mesh.shape[i_dim_device_mesh] != 0:
                    raise ValueError(
                        f"Uneven sharding tensor dimension {placement.dim} of size {x.shape[placement.dim]} "
                        f"along device mesh dimension {i_dim_device_mesh} of size {device_mesh.shape[i_dim_device_mesh]} "
                        "is not supported"
                    )
                placements_grad_params[i_dim_device_mesh] = Partial("sum")
            elif not isinstance(placement, Replicate):
                raise ValueError(
                    f"Unsupported x's placements along {i_dim_device_mesh} axis of the device mesh: {placement}"
                )

        x_local = x.to_local()
        weight_local = weight.to_local()

        needs_grad = weight.requires_grad
        if needs_grad:
            with torch.enable_grad():
                weight_local_detached = weight_local.detach().requires_grad_(True)
                output_local = F.embedding(
                    x_local,
                    weight_local_detached,
                    padding_idx=padding_idx,
                )
            ctx.save_for_backward(output_local, weight_local_detached)
        else:
            output_local = F.embedding(
                x_local,
                weight_local,
                padding_idx=padding_idx,
            )

        shape_output = list(output_local.shape)
        for i_dim_mesh, placement in enumerate(x.placements):
            if isinstance(placement, Shard):
                shape_output[placement.dim] *= device_mesh.shape[i_dim_mesh]

        if needs_grad:
            ctx.device_mesh = device_mesh
            ctx.placements_x = x.placements
            ctx.placements_grad_params = placements_grad_params
            ctx.weight_shape = weight.shape
            ctx.weight_stride = weight.stride()
            ctx.output_shape = torch.Size(shape_output)

        stride_output = update_exhaustive_strides(output_local.shape, output_local.stride(), shape_output)
        output = DTensor.from_local(
            output_local,
            device_mesh,
            x.placements,
            shape=torch.Size(shape_output),
            stride=torch.Size(stride_output),
        )
        return output

    @staticmethod
    @torch.amp.custom_bwd(device_type="cuda")
    def backward(  # type: ignore[override]
        ctx, grad_output: DTensor
    ) -> tuple[None, Optional[DTensor], None]:
        """
        Backward pass for the distributed embedding operation.

        Assumptions and requirements:
        1. Gradient of the output is on the same device mesh as the output

        Args:
            ctx: Context object with stored information from forward pass
            grad_output: Gradient of the loss with respect to the output

        Returns:
            Tuple of gradients for input, weight, and padding_idx.
        """
        if grad_output.device_mesh != ctx.device_mesh:
            raise ValueError(
                "_EmbeddingParamsReplicatedImpl: different device mesh between grad_output and the forward input"
            )

        if grad_output.shape != ctx.output_shape:
            raise ValueError(
                "_EmbeddingParamsReplicatedImpl: different shape between grad_output and the forward output"
            )

        grad_weight = None
        if ctx.needs_input_grad[1]:
            grad_output_local = grad_output.to_local()
            output_local, weight_local_detached = ctx.saved_tensors
            (grad_weight_local,) = torch.autograd.grad(
                outputs=[output_local],
                inputs=[weight_local_detached],
                grad_outputs=[grad_output_local],
                retain_graph=False,
            )
            grad_weight = DTensor.from_local(
                grad_weight_local,
                ctx.device_mesh,
                ctx.placements_grad_params,
                shape=ctx.weight_shape,
                stride=ctx.weight_stride,
            )

        return None, grad_weight, None


class EmbeddingParamsReplicated(nn.Module):
    """
    Distributed embedding layer with parameters replicated across all device mesh dimensions.

    This is almost equivalent to
    ```python
    layer = torch.distributed.tensor.distribute_module(layer_local, device_mesh)
    ```
    with the exception that the torch.distributed.tensor.distribute_module version will incur
    significant overhead due to the unnecessary replication of the output tensor along certain
    device mesh dimensions.

    This class avoids such unnecessary overhead by using the custom _EmbeddingParamsReplicatedImpl
    autograd function for forward and backward pass computation instead of relying on the distributed
    module's forward implementation.

    Key requirements:
        1. Parameters (weight) will be replicated on all device mesh dimensions
        2. Input tensor and parameters must be on the same device mesh
        3. Partial reduction along any input dimension is not supported
        4. Input indices must be integer dtype and must not require gradients
        5. Input and outputs must be on the same device mesh with the same placements
        6. Gradients of the weight have Partial("sum") placements along the input's Shard placements'
           dimension so that the all-reduce will be performed along those device-grid dimensions
    """

    def __init__(self, layer_local: nn.Embedding, device_mesh: DeviceMesh):
        if not isinstance(layer_local, nn.Embedding):
            raise TypeError("layer_local is not an instance of nn.Embedding")
        if layer_local.weight.device.type != device_mesh.device_type:
            raise ValueError(
                f"layer_local.weight and device_mesh are not on the same device type: "
                f"{layer_local.weight.device.type} != {device_mesh.device_type}"
            )
        if layer_local.sparse:
            raise ValueError("sparse option is not supported in EmbeddingParamsReplicated")
        if layer_local.scale_grad_by_freq:
            raise ValueError("scale_grad_by_freq option is not supported in EmbeddingParamsReplicated")
        if layer_local.max_norm is not None:
            raise ValueError("max_norm option is not supported in EmbeddingParamsReplicated")

        super().__init__()
        all_replicate_placements = [Replicate()] * device_mesh.ndim
        self.weight = nn.Parameter(
            distribute_tensor(layer_local.weight.data, device_mesh, all_replicate_placements),
            requires_grad=layer_local.weight.requires_grad,
        )
        self.padding_idx = layer_local.padding_idx
        self.num_embeddings = layer_local.num_embeddings
        self.embedding_dim = layer_local.embedding_dim

    def forward(self, input: DTensor) -> DTensor:
        """
        Forward pass for the distributed embedding layer.

        Uses the custom _EmbeddingParamsReplicatedImpl autograd function to perform the computation
        efficiently while preserving correct autograd behavior for distributed tensors.

        Args:
            input: Input index DTensor

        Returns:
            Output DTensor with same placement strategy as input
        """
        return cast(
            DTensor,
            _EmbeddingParamsReplicatedImpl.apply(
                input,
                self.weight,
                self.padding_idx,
            ),
        )
