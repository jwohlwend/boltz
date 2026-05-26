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


from typing import Optional, Union

import torch
import torch.distributed as dist
import torch.nn as nn
import torch.nn.functional as F
from torch import Tensor
from torch.distributed.tensor import DeviceMesh, DTensor, Partial, Replicate, Shard, distribute_tensor

_shape_t = Union[int, list[int], torch.Size]


class _ContextParallelLayerNormImpl(torch.autograd.Function):
    @staticmethod
    def forward(
        ctx,
        x: Tensor,
        normalized_shape: list[int],
        weight: Optional[Tensor],
        bias: Optional[Tensor],
        eps: float,
        reduce_group: dist.ProcessGroup,
    ) -> Tensor:
        """Forward pass of the layer normalization.

        Args:
            ctx: context
            x (Tensor): input tensor
            normalized_shape (list[int]): shape of the input tensor
            weight (Optional[Tensor]): weight tensor
            bias (Optional[Tensor]): bias tensor
            eps (float): a value added to the denominator for numerical stability
            reduce_group (dist.ProcessGroup): process group for all-reduce

        Returns:
            output tensor (Tensor)
        """
        weight_needs_grad = weight is not None and weight.requires_grad
        bias_needs_grad = bias is not None and bias.requires_grad
        # For unknown reasons, using ctx.needs_input_grad in the forward pass can occasionally
        # cause NCCL hanging. ctx.need_input_grad should not be accessed during the forward pass
        # according to this discussion on pytorch forum:
        # https://discuss.pytorch.org/t/is-there-a-diffrence-between-ctx-needs-input-grad-behaviour-vs-input-tensor-requires-grad/195063/2
        if x.requires_grad or weight_needs_grad or bias_needs_grad:
            ctx.reduce_group = reduce_group
            ctx.eps = eps
            ctx.normalized_shape = normalized_shape

            if not x.requires_grad:
                weight = None

            ctx.save_for_backward(x, weight)

        return F.layer_norm(x, normalized_shape, weight, bias, eps)

    @staticmethod
    def backward(
        ctx, grad_output: Tensor
    ) -> tuple[Optional[Tensor], None, Optional[Tensor], Optional[Tensor], None, None]:
        """Backward pass of the layer normalization.

        Although the output between tokens is independent, the backward pass on the weight and bias tensors involves a summation over all tokens, and thus requires all-reduce due to context parallelism.

        Args:
            ctx: context
            grad_output (Tensor): gradient of the output tensor

        Returns:
            gradient for input, weight, bias, eps, and reduce_group (tuple[Optional[Tensor], None, Optional[Tensor], Optional[Tensor], None, None])
        """
        x, weight = ctx.saved_tensors
        eps = ctx.eps
        normalized_shape = ctx.normalized_shape

        if ctx.needs_input_grad[0] or ctx.needs_input_grad[2]:
            dims = tuple(-(i + 1) for i in range(len(normalized_shape)))
            mean = x.mean(dim=dims, keepdim=True)
            var = x.var(dim=dims, unbiased=False, keepdim=True)
            x_norm = (x - mean) / torch.sqrt(var + eps)

        if ctx.needs_input_grad[0]:
            if weight is not None:
                dy = grad_output * weight.view(*([1] * (grad_output.ndim - len(normalized_shape))), *weight.shape)
            else:
                dy = grad_output

            dims = tuple(-(i + 1) for i in range(len(normalized_shape)))
            dy_mean = dy.mean(dim=dims, keepdim=True)
            dy_x_norm_mean = (dy * x_norm).mean(dim=dims, keepdim=True)
            grad_input = (dy - dy_mean - x_norm * dy_x_norm_mean) / torch.sqrt(var + eps)
        else:
            grad_input = None

        if ctx.needs_input_grad[2]:
            reduce_dims = list(range(grad_output.ndim - len(normalized_shape)))
            grad_weight = (grad_output * x_norm).sum(dim=reduce_dims)
            grad_weight = grad_weight.contiguous()
            grad_weight_work = dist.all_reduce(grad_weight, op=dist.ReduceOp.SUM, group=ctx.reduce_group, async_op=True)
        else:
            grad_weight = None

        if ctx.needs_input_grad[3]:
            reduce_dims = list(range(grad_output.ndim - len(normalized_shape)))
            grad_bias = grad_output.sum(dim=reduce_dims)
            grad_bias_work = dist.all_reduce(grad_bias, op=dist.ReduceOp.SUM, group=ctx.reduce_group, async_op=True)
        else:
            grad_bias = None

        # collect all_reduce results at the end
        if ctx.needs_input_grad[2]:
            grad_weight_work.wait()
        if ctx.needs_input_grad[3]:
            grad_bias_work.wait()

        return grad_input, None, grad_weight, grad_bias, None, None


class ContextParallelLayerNorm(nn.LayerNorm):
    def __init__(
        self,
        normalized_shape: _shape_t,
        reduce_group: dist.ProcessGroup,
        eps: float = 1e-5,
        elementwise_affine: bool = True,
        bias: bool = True,
        device=None,
        dtype=None,
    ):
        """Context parallel layer normalization, a wrapper around nn.LayerNorm that supports distributed training.

        Although the output between tokens is independent, the backward pass on the weight and bias tensors involves a summation over all tokens, and thus requires all-reduce due to context parallelism. This means we need a dedicated ContextParallelLayerNorm class for training.

        Args:
            normalized_shape (int or list or torch.Size): input shape from an expected input
                of size
            reduce_group (dist.ProcessGroup): The process group to use for gradient all-reduce.
            eps (float): a value added to the denominator for numerical stability. Default: 1e-5
            elementwise_affine (bool): a boolean value that when set to ``True``, this module
                has learnable per-element affine parameters initialized to ones (for weights)
                and zeros (for biases). Default: ``True``.
            bias (bool): If set to ``False``, the layer will not learn an additive bias (only relevant if
                :attr:`elementwise_affine` is ``True``). Default: ``True``.
            device (torch.device, optional): device on which the module is allocated. Defaults to None.
            dtype (torch.dtype, optional): dtype of the module. Defaults to None.
        """
        super().__init__(normalized_shape, eps, elementwise_affine, bias, device, dtype)
        assert reduce_group is not None, "reduce_group must be provided"
        self.reduce_group = reduce_group

        if isinstance(normalized_shape, int):
            self._normalized_shape_list = [normalized_shape]
        else:
            self._normalized_shape_list = list(normalized_shape)

    def forward(self, input: Tensor) -> Tensor:
        return _ContextParallelLayerNormImpl.apply(
            input, self._normalized_shape_list, self.weight, self.bias, self.eps, self.reduce_group
        )


def get_cp_layernorm(
    normalized_shape: _shape_t,
    reduce_group: dist.ProcessGroup | None = None,
    eps: float = 1e-5,
    elementwise_affine: bool = True,
    bias: bool = True,
    device: torch.device | None = None,
    dtype: torch.dtype | None = None,
) -> nn.LayerNorm | ContextParallelLayerNorm:
    """Get a layer normalization module that is optimized for distributed training.

    Args:
        normalized_shape (int or list or torch.Size): input shape from an expected input
            of size
        reduce_group (dist.ProcessGroup, optional): The process group to use for gradient all-reduce.
            Defaults to None.
        eps (float, optional): a value added to the denominator for numerical stability. Default: 1e-5
        elementwise_affine (bool, optional): a boolean value that when set to ``True``, this module
            has learnable per-element affine parameters initialized to ones (for weights)
            and zeros (for biases). Default: ``True``.
        bias (bool, optional): If set to ``False``, the layer will not learn an additive bias (only relevant if
            :attr:`elementwise_affine` is ``True``). Default: ``True``.
        device (torch.device, optional): device on which the module is allocated. Defaults to None.
        dtype (torch.dtype, optional): dtype of the module. Defaults to None.

    Returns:
        nn.LayerNorm | ContextParallelLayerNorm
    """
    if reduce_group is None:
        return nn.LayerNorm(
            normalized_shape,
            eps=eps,
            elementwise_affine=elementwise_affine,
            bias=bias,
            device=device,
            dtype=dtype,
        )
    else:
        return ContextParallelLayerNorm(
            normalized_shape,
            reduce_group=reduce_group,
            eps=eps,
            elementwise_affine=elementwise_affine,
            bias=bias,
        )


class _LayerNormParamsReplicatedImpl(torch.autograd.Function):
    """
    A custom implementation of LayerNorm with replicated parameters for distributed training.

    This class provides a forward and backward implementation of LayerNorm, ensuring compatibility
    with distributed tensor placements and device meshes. It supports replicated and sharded
    placements for input tensors and replicated placements for weight and bias tensors.

    NOTE: by default, avg reduce over the Replicate placements of the weight and bias gradients
    is performed. This is to ensure identical parameter updates across all ranks and avoid
    gradual divergence during training. This can be disabled by setting
    avg_over_replicate_param_grad to False.

    Methods:
        forward(ctx, x, normalized_shape, weight, bias, eps):
            Computes the forward pass of LayerNorm.

        backward(ctx, grad_output):
            Computes the backward pass of LayerNorm.
    """

    @staticmethod
    @torch.amp.custom_fwd(device_type="cuda")
    def forward(
        ctx,
        x: DTensor,
        normalized_shape: list[int],
        weight: Optional[DTensor],
        bias: Optional[DTensor],
        eps: float,
        cast_params_dtype_to_x: Optional[bool] = False,
        avg_over_replicate_param_grad: bool = True,
    ) -> DTensor:
        """
        Forward pass of LayerNorm with replicated parameters.

        Args:
            ctx: Context for saving tensors for backward computation.
            x (DTensor): Input tensor.
            normalized_shape (list[int]): Shape of the input tensor to normalize.
            weight (Optional[DTensor]): Weight tensor for affine transformation.
            bias (Optional[DTensor]): Bias tensor for affine transformation.
            eps (float): A small value added for numerical stability.
            cast_params_dtype_to_x (Optional[bool]): whether to cast the dtype of
                the weights and bias to the dtype of the input tensor
            avg_over_replicate_param_grad (bool): Whether to perform avg reduce over the
                Replicate placements of the weight and bias gradients. For example,
                if the input DTensor x.placements = (Shard(0), Replicate()), this layer's
                parameters' gradients.placements = (Partial("sum"), Replicate()) if
                self._avg_over_replicate_param_grad is False; otherwise, it will be
                (Partial("sum"), Partial("avg")). The motivation is to ensure identical
                parameter updates across all ranks and avoid gradual divergence during
                training.

        Returns:
            DTensor: The normalized output tensor.
        """
        if not isinstance(x, DTensor):
            dtensor_instance = x
            raise TypeError(
                ", ".join(
                    [
                        f"DTensor instance '{dtensor_instance}' should have type {DTensor}",
                        f"but instead has type {type(dtensor_instance)}.",
                    ]
                )
            )
        device_mesh = x.device_mesh
        ndim_device_mesh = device_mesh.ndim
        all_replicate_placements = tuple([Replicate()] * ndim_device_mesh)
        if weight is not None:
            if not isinstance(weight, DTensor):
                dtensor_instance = weight
                raise TypeError(
                    ", ".join(
                        [
                            f"DTensor instance '{dtensor_instance}' should have type {DTensor}",
                            f"but instead has type {type(dtensor_instance)}.",
                        ]
                    )
                )
            if weight.device_mesh != device_mesh:
                raise ValueError("weight and x must be on the same device mesh")
            if weight.placements != all_replicate_placements:
                raise ValueError("weight must be replicated on all device mesh dimensions")
        if bias is not None:
            if not isinstance(bias, DTensor):
                dtensor_instance = bias
                raise TypeError(
                    ", ".join(
                        [
                            f"DTensor instance '{dtensor_instance}' should have type {DTensor}",
                            f"but instead has type {type(dtensor_instance)}.",
                        ]
                    )
                )
            if bias.device_mesh != device_mesh:
                raise ValueError("bias and x must be on the same device mesh")
            if bias.placements != all_replicate_placements:
                raise ValueError("bias must be replicated on all device mesh dimensions")
        if weight is not None or bias is not None:
            if avg_over_replicate_param_grad:
                placements_grad_params = [Partial("avg")] * ndim_device_mesh
            else:
                placements_grad_params = list(weight.placements) if weight is not None else None
        else:
            placements_grad_params = None
        n_dim_norm = len(normalized_shape)
        for i_dim_device_mesh, p in enumerate(x.placements):
            if isinstance(p, Partial):
                # partial reduction along any input dimension requires complicated backward pass
                raise ValueError("Partial reduction along any input dimension is not supported")
            if isinstance(p, Shard):
                if p.dim >= x.ndim - n_dim_norm:
                    # the normalized dimensions must not be sharded by the device mesh
                    raise ValueError("LayerNorm's normalizing dimensions must not be sharded by the device mesh")
                if x.shape[p.dim] % device_mesh.shape[i_dim_device_mesh] != 0:
                    raise ValueError(
                        f"Uneven sharding tensor dimension {p.dim} of size {x.shape[p.dim]} "
                        f"along device mesh dimension {i_dim_device_mesh} of size "
                        f"{device_mesh.shape[i_dim_device_mesh]} is not supported"
                    )
                # the only supported placement for the input is Shard, which corresponding
                # to the backward's grad partial sum. Otherwise, we can only support Replicate
                # placements for other device mesh dimensions. Also, by using the Partial("sum")
                # placement on the params, the all_reduce is postponed for the params' gradients
                # until needed
                if weight is not None or bias is not None:
                    placements_grad_params[i_dim_device_mesh] = Partial("sum")
            elif not isinstance(p, Replicate):
                raise ValueError(f"Unsupported x's placements along {i_dim_device_mesh} axis of the device mesh: {p}")
        ctx.device_mesh = device_mesh
        # will use x.placements for the x.grad in the backward pass, i.e., this function
        # enforces consistent placements for the input and its gradient
        ctx.placements_x = x.placements
        ctx.placements_grad_params = placements_grad_params

        # Save weight and bias shapes and strides for backward pass
        if weight is not None:
            ctx.weight_shape = weight.shape
            ctx.weight_stride = weight.stride()
        if bias is not None:
            ctx.bias_shape = bias.shape
            ctx.bias_stride = bias.stride()

        weight_needs_grad = weight is not None and weight.requires_grad
        bias_needs_grad = bias is not None and bias.requires_grad
        # IMPORTANT: no modification on *_local for the rest of the code
        x_local = x.to_local()
        if weight is None:
            weight_local = None
        else:
            weight_local = weight.to_local()
            if cast_params_dtype_to_x:
                weight_local = weight_local.to(x.dtype)
        if bias is None:
            bias_local = None
        else:
            bias_local = bias.to_local()
            if cast_params_dtype_to_x:
                bias_local = bias_local.to(x.dtype)
        # For unknown reasons, using ctx.needs_input_grad in the forward pass can occasionally
        # cause NCCL hanging. ctx.need_input_grad should not be accessed during the forward pass
        # according to this discussion on pytorch forum:
        # https://discuss.pytorch.org/t/is-there-a-diffrence-between-ctx-needs-input-grad-behaviour-vs-input-tensor-requires-grad/195063/2
        if x.requires_grad or weight_needs_grad or bias_needs_grad:
            ctx.eps = eps
            ctx.normalized_shape = normalized_shape

            if not x.requires_grad:
                weight = None

            ctx.save_for_backward(x_local, weight_local)

        output_local = F.layer_norm(x_local, normalized_shape, weight_local, bias_local, eps)
        # LayerNorm does not change input's shape
        output = DTensor.from_local(
            output_local,
            device_mesh=device_mesh,
            placements=x.placements,
            shape=x.shape,
            stride=x.stride(),
        )
        return output

    @staticmethod
    @torch.amp.custom_bwd(device_type="cuda")
    def backward(
        ctx, grad_output: DTensor
    ) -> tuple[Optional[DTensor], None, Optional[DTensor], Optional[DTensor], None, None]:
        """
        Backward pass of LayerNorm with replicated parameters.

        Args:
            ctx: Context containing saved tensors and attributes from the forward pass.
            grad_output (DTensor): Gradient of the output tensor.

        Returns:
            tuple: Gradients for input, weight, bias, and other parameters.
        """
        x_local, weight_local = ctx.saved_tensors
        eps = ctx.eps
        normalized_shape = ctx.normalized_shape

        # IMPORTANT: no modification on *_local for the rest of the code
        grad_output_local = grad_output.to_local()

        ids_dim_norm = tuple(-(i + 1) for i in range(len(normalized_shape)))
        if ctx.needs_input_grad[0] or ctx.needs_input_grad[2]:
            mean_local = x_local.mean(dim=ids_dim_norm, keepdim=True)
            var_local = x_local.var(dim=ids_dim_norm, unbiased=False, keepdim=True)
            x_norm_local = (x_local - mean_local) / torch.sqrt(var_local + eps)

        if ctx.needs_input_grad[0]:
            if weight_local is not None:
                dy_local = grad_output_local * weight_local.view(
                    *([1] * (grad_output_local.ndim - len(normalized_shape))), *weight_local.shape
                )
            else:
                dy_local = grad_output_local

            dy_mean_local = dy_local.mean(dim=ids_dim_norm, keepdim=True)
            dy_x_norm_mean_local = (dy_local * x_norm_local).mean(dim=ids_dim_norm, keepdim=True)
            grad_input_local = (dy_local - dy_mean_local - x_norm_local * dy_x_norm_mean_local) / torch.sqrt(
                var_local + eps
            )
            # LayerNorm does not change input's shape in both forward and backward passes
            grad_input = DTensor.from_local(
                grad_input_local,
                device_mesh=ctx.device_mesh,
                placements=ctx.placements_x,
                shape=grad_output.shape,
                stride=grad_output.stride(),
            )
        else:
            grad_input = None

        reduce_dims = list(range(grad_output_local.ndim - len(normalized_shape)))
        if ctx.needs_input_grad[2]:
            grad_weight_local = (grad_output_local * x_norm_local).sum(dim=reduce_dims)
            # all-replicate weight implies identical shape and stride across all ranks
            grad_weight = DTensor.from_local(
                grad_weight_local,
                device_mesh=ctx.device_mesh,
                placements=ctx.placements_grad_params,
                shape=ctx.weight_shape,
                stride=ctx.weight_stride,
            )
        else:
            grad_weight = None

        if ctx.needs_input_grad[3]:
            grad_bias_local = grad_output_local.sum(dim=reduce_dims)
            # all-replicate weight implies identical shape and stride across all ranks
            grad_bias = DTensor.from_local(
                grad_bias_local,
                device_mesh=ctx.device_mesh,
                placements=ctx.placements_grad_params,
                shape=ctx.bias_shape,
                stride=ctx.bias_stride,
            )
        else:
            grad_bias = None

        return grad_input, None, grad_weight, grad_bias, None, None, None


class LayerNormParamsReplicated(nn.Module):
    """
    A LayerNorm module with replicated parameters for distributed training.

    This module wraps around `_LayerNormParamsReplicatedImpl` to provide a user-friendly interface
    for LayerNorm operations using the DTensor API. It supports distributed training with replicated
    and sharded placements for input tensors and replicated placements for weight and bias tensors.

    NOTE: by default, avg reduce over the Replicate placements of the weight and bias gradients
    is performed. This is to ensure identical parameter updates across all ranks and avoid
    gradual divergence during training. This can be disabled by setting
    avg_over_replicate_param_grad to False.

    Args:
        layer_local (nn.LayerNorm): An already-initialized nn.LayerNorm instance.
        device_mesh (DeviceMesh): The device mesh for distributed training.
        avg_over_replicate_param_grad (bool): Whether to perform avg reduce over the
            Replicate placements of the weight and bias gradients. For example,
            if the input DTensor x.placements = (Shard(0), Replicate()), this layer's
            parameters' gradients.placements = (Partial("sum"), Replicate()) if
            self._avg_over_replicate_param_grad is False; otherwise, it will be
            (Partial("sum"), Partial("avg")). The motivation is to ensure identical
            parameter updates across all ranks and avoid gradual divergence during
            training.
    """

    def __init__(
        self, layer_local: nn.LayerNorm, device_mesh: DeviceMesh, avg_over_replicate_param_grad: bool = True
    ) -> None:
        if not isinstance(layer_local, nn.LayerNorm):
            raise TypeError("layer_local is not an instance of nn.LayerNorm")
        if layer_local.weight is not None and layer_local.weight.device.type != device_mesh.device_type:
            raise ValueError(
                f"layer_local.weight and device_mesh are not on the same device type: "
                f"{layer_local.weight.device.type} != {device_mesh.device_type}"
            )
        if layer_local.bias is not None and layer_local.bias.device.type != device_mesh.device_type:
            raise ValueError(
                f"layer_local.bias and device_mesh are not on the same device type: "
                f"{layer_local.bias.device.type} != {device_mesh.device_type}"
            )

        super().__init__()
        self.normalized_shape = layer_local.normalized_shape
        self.eps = layer_local.eps
        self.device_mesh = device_mesh
        self.elementwise_affine = layer_local.elementwise_affine
        self._avg_over_replicate_param_grad = avg_over_replicate_param_grad

        all_replicate_placements = [Replicate()] * device_mesh.ndim

        if layer_local.weight is None:
            self.register_parameter("weight", None)
        else:
            self.weight = nn.Parameter(
                distribute_tensor(layer_local.weight.data, device_mesh, all_replicate_placements)
            )
        if layer_local.bias is None:
            self.register_parameter("bias", None)
        else:
            self.bias = nn.Parameter(distribute_tensor(layer_local.bias.data, device_mesh, all_replicate_placements))

    def forward(self, x: DTensor) -> DTensor:
        """
        Forward pass of LayerNormParamsReplicated.

        Args:
            x (DTensor): Input tensor.

        Returns:
            DTensor: The normalized output tensor.
        """
        return _LayerNormParamsReplicatedImpl.apply(
            x,
            self.normalized_shape,
            self.weight,
            self.bias,
            self.eps,
            True,
            self._avg_over_replicate_param_grad,
        )
