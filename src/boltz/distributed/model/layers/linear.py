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
import torch.distributed as dist
import torch.nn as nn
import torch.nn.functional as F
from torch import Tensor
from torch.distributed.device_mesh import DeviceMesh
from torch.distributed.tensor import DTensor, Partial, Replicate, Shard, distribute_tensor

from boltz.distributed.utils import update_exhaustive_strides


class _ContextParallelLinearImpl(torch.autograd.Function):
    @staticmethod
    def forward(ctx, x: Tensor, weight: Tensor, bias: Optional[Tensor], reduce_group: dist.ProcessGroup) -> Tensor:
        """Forward pass of the linear layer.

        Args:
            ctx: context
            x (Tensor): input tensor
            weight (Tensor): weight tensor
            bias (Optional[Tensor]): bias tensor
            reduce_group (dist.ProcessGroup): process group for all-reduce

        Returns:
            output tensor (Tensor)
        """
        # For unknown reasons, using ctx.needs_input_grad in the forward pass can occasionally
        # cause NCCL hanging. ctx.need_input_grad should not be accessed during the forward pass
        # according to this discussion on pytorch forum:
        # https://discuss.pytorch.org/t/is-there-a-diffrence-between-ctx-needs-input-grad-behaviour-vs-input-tensor-requires-grad/195063/2
        if x.requires_grad or weight.requires_grad or (bias is not None and bias.requires_grad):
            ctx.reduce_group = reduce_group
            ctx.save_for_backward(
                x if weight.requires_grad else None,
                weight if x.requires_grad else None,
            )
        return F.linear(x, weight, bias)

    @staticmethod
    def backward(ctx, grad_output: Tensor) -> tuple[Optional[Tensor], Optional[Tensor], Optional[Tensor], None]:
        """Backward pass of the linear layer.

        Although the output between tokens is independent, the backward pass on the weight and bias tensors involves a summation over all tokens, and thus requires all-reduce due to context parallelism.

        Args:
            ctx: context
            grad_output (Tensor): gradient of the output tensor

        Returns:
            gradient for input, weight, bias, and reduce_group (tuple[Optional[Tensor], Optional[Tensor], Optional[Tensor], None])
        """
        x, weight = ctx.saved_tensors

        if ctx.needs_input_grad[1]:
            dw = torch.einsum("...i,...o->io", grad_output, x)
            dw = dw.contiguous()
            dw_work = dist.all_reduce(dw, op=dist.ReduceOp.SUM, group=ctx.reduce_group, async_op=True)
        else:
            dw = None
            dw_work = None

        if ctx.needs_input_grad[2]:
            dims = list(range(grad_output.ndim - 1))  # aggregate over all but the last dimension
            db = grad_output.sum(dim=dims)
            db = db.contiguous()
            db_work = dist.all_reduce(db, op=dist.ReduceOp.SUM, group=ctx.reduce_group, async_op=True)
        else:
            db = None
            db_work = None

        if ctx.needs_input_grad[0]:
            grad_input = torch.einsum("...i,io->...o", grad_output, weight)
        else:
            grad_input = None

        # collect all work
        if dw_work is not None:
            dw_work.wait()
        if db_work is not None:
            db_work.wait()

        return grad_input, dw, db, None


class _LinearParamsReplicatedImpl(torch.autograd.Function):
    """
    Custom autograd Function implementation for distributed linear operation with replicated parameters.

    The main purpose of this implementation is to avoid the unnecessary overhead seen the the
    equivalent distribute_module-wrapped linear layer, where the output tensors have nonsensical Replicate
    placements along device mesh dimensions that are not intended

    This implementation handles the forward and backward passes for a distributed linear layer where
    parameters (weight and bias) are replicated across the device mesh. The input tensor can have
    various placement strategies.

    NOTE: by default, avg reduce over the Replicate placements of the weight and bias gradients
    is performed. This is to ensure identical parameter updates across all ranks and avoid
    gradual divergence during training. This can be disabled by setting
    avg_over_replicate_param_grad to False.

    Assumptions and requirements:
        (see the respective docstring for forward and backward)
    """

    @staticmethod
    @torch.amp.custom_fwd(device_type="cuda")
    def forward(
        ctx,
        x: DTensor,
        weight: DTensor,
        bias: Optional[DTensor],
        cast_params_dtype_to_x: bool = False,
        avg_over_replicate_param_grad: bool = True,
    ) -> DTensor:
        """
        Forward pass for the distributed linear operation.

        Assumptions and requirements:
        1. Parameters (weight and bias) must be replicated on all device mesh dimensions
        2. Input tensor and parameters must be on the same device mesh
        3. Feature/hidden dimension of the input must not be sharded across the device mesh
        4. Partial reduction along any input dimension is not supported
        5. Input and outputs must be on the same device mesh with the same placements

        Args:
            ctx: Context object to store information for backward pass
            x: Input tensor with arbitrary placement strategy
            weight: Weight tensor (must be replicated across all device mesh dimensions)
            bias: Optional bias tensor (must be replicated if provided)
            cast_params_dtype_to_x: whether to cast the dtype of the weight and bias
                to the dtype of the input tensor
            avg_over_replicate_param_grad: whether to perform avg reduce over the
                Replicate placements of the weight and bias gradients. For example,
                if the input DTensor x.placements = (Shard(0), Replicate()), this layer's
                parameters' gradients.placements = (Partial("sum"), Replicate()) if
                self._avg_over_replicate_param_grad is False; otherwise, it will be
                (Partial("sum"), Partial("avg")). The motivation is to ensure identical
                parameter updates across all ranks and avoid gradual divergence during
                training.

        Returns:
            Output tensor with same placement strategy as input

        Raises:
            ValueError: If any of the placement requirements are violated
        """
        device_mesh = x.device_mesh
        if weight.device_mesh != device_mesh:
            raise ValueError("weight and x must be on the same device mesh")
        if bias is not None and bias.device_mesh != device_mesh:
            raise ValueError("bias and x must be on the same device mesh")
        ndim_device_mesh = device_mesh.ndim
        all_replicate_placements = tuple([Replicate()] * ndim_device_mesh)
        if weight.placements != all_replicate_placements:
            raise ValueError("weight must be replicated on all device mesh dimensions")
        if bias is not None and bias.placements != all_replicate_placements:
            raise ValueError("bias must be replicated on all device mesh dimensions")
        if avg_over_replicate_param_grad:
            placements_grad_params = [Partial("avg")] * ndim_device_mesh
        else:
            # all-replicate placements
            placements_grad_params = list(weight.placements)
        for i_dim_device_mesh, p in enumerate(x.placements):
            if isinstance(p, Partial):
                # partial reduction along any input dimension requires complicated backward pass
                raise ValueError("Partial reduction along any input dimension is not supported")
            if isinstance(p, Shard):
                if p.dim == x.ndim - 1:
                    # the feature or hidden dimension must not be a part of the device mesh
                    raise ValueError("feature or hidden dimension must not be a part of the device mesh")
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
                placements_grad_params[i_dim_device_mesh] = Partial("sum")
            elif not isinstance(p, Replicate):
                raise ValueError(f"Unsupported x's placements along {i_dim_device_mesh} axis of the device mesh: {p}")
        ctx.device_mesh = device_mesh
        # will use x.placements for the x.grad in the backward pass, i.e., this function
        # enforces consistent placements for the input and its gradient
        ctx.placements_x = x.placements
        ctx.placements_grad_params = placements_grad_params
        ctx.shape_input = x.shape
        ctx.stride_input = x.stride()
        ctx.weight_shape = weight.shape
        ctx.weight_stride = weight.stride()
        ctx.dtype_input = x.dtype
        ctx.dtype_weight = weight.dtype
        if bias is not None:
            ctx.bias_shape = bias.shape
            ctx.bias_stride = bias.stride()
            ctx.dtype_bias = bias.dtype
        else:
            ctx.dtype_bias = None
        x_local = x.to_local()
        weight_local = weight.to_local()
        bias_local = None if bias is None else bias.to_local()

        # Save original-precision locals for backward *before* any dtype cast.
        # Native autocast saves fp32 weights and lets the backward autocast
        # context handle further casts.  Saving the bf16-cast version would
        # bake in bf16 rounding on CPU (where custom_bwd does NOT restore
        # autocast), silently lowering gradient precision.
        if x.requires_grad or weight.requires_grad or (bias is not None and bias.requires_grad):
            ctx.save_for_backward(
                x_local.detach().clone() if weight.requires_grad else None,
                weight_local.detach().clone() if x.requires_grad else None,
            )

        if cast_params_dtype_to_x:
            weight_local = weight_local.to(x.dtype)
            if bias_local is not None:
                bias_local = bias_local.to(x.dtype)
        # Extract the local shard to perform the linear operation.
        # This enforces local matrix multiplication without any communication given that:
        # 1. the linear operation is performed locally on each rank along the hidden dimension,
        #    which is agnostic to the device mesh dimensions
        # 2. the weight and bias are replicated on all device mesh dimensions
        # 3. the output has the same placements as the input
        output_local = torch.nn.functional.linear(x_local, weight_local, bias_local)
        # linear only change the last dimension of the input so we need to
        # modify the output shape and strides accordingly
        shape_output = tuple(x.shape[:-1]) + (output_local.shape[-1],)
        strides_output = update_exhaustive_strides(x.shape, x.stride(), shape_output)
        output = DTensor.from_local(output_local, device_mesh, x.placements, shape=shape_output, stride=strides_output)
        return output

    @staticmethod
    def _all_reduce_grad_gteqfp32(
        grad: torch.Tensor,
        device_mesh: DeviceMesh,
        placements: list,
        target_dtype: torch.dtype,
    ) -> torch.Tensor:
        """All-reduce a parameter gradient in at least fp32 across mesh dims.

        For each mesh dimension with a ``Partial`` placement, performs an
        all-reduce in at least float32 to avoid bf16/fp16 accumulation errors.
        Only the parameter-sized gradient is promoted — not the large
        activation tensors.  If the gradient is already >=fp32, it is reduced
        in its native dtype.
        """
        needs_reduce = any(isinstance(p, Partial) and device_mesh.size(dim) > 1 for dim, p in enumerate(placements))
        if not needs_reduce:
            return grad.to(target_dtype)

        reduce_dtype = torch.promote_types(grad.dtype, torch.float32)
        grad = grad.to(reduce_dtype).contiguous()
        for mesh_dim, p in enumerate(placements):
            if not isinstance(p, Partial) or device_mesh.size(mesh_dim) <= 1:
                continue
            group = device_mesh.get_group(mesh_dim)
            if p.reduce_op == "sum":
                dist.all_reduce(grad, op=dist.ReduceOp.SUM, group=group)
            elif p.reduce_op == "avg":
                dist.all_reduce(grad, op=dist.ReduceOp.AVG, group=group)
            else:
                raise ValueError(f"Unsupported reduce_op {p.reduce_op!r} in _all_reduce_grad_gteqfp32")
        return grad.to(target_dtype)

    @staticmethod
    @torch.amp.custom_bwd(device_type="cuda")
    def backward(
        ctx, grad_output: DTensor
    ) -> tuple[Optional[DTensor], Optional[DTensor], Optional[DTensor], None, None]:
        """Backward pass for the distributed linear operation.

        Local einsum stays in the compute dtype (bf16 MMA accumulates in
        fp32 internally on CUDA, so local results are accurate).  The
        precision hazard is cross-rank reduction: implicit ``Partial("SUM")``
        would reduce in bf16, accumulating errors.  We manually all-reduce
        ``dw``/``db`` in fp32 via ``_all_reduce_grad_gteqfp32`` and return
        them with ``Replicate`` placements.

        On CPU (unit tests), ``custom_bwd`` does not restore autocast, so
        we explicitly cast operands to the compute dtype.
        """
        if grad_output.device_mesh != ctx.device_mesh:
            raise ValueError(
                "_LinearParamsReplicatedImpl: different device mesh between grad_output and the forward input"
            )
        x_local, weight_local = ctx.saved_tensors

        if grad_output.placements != ctx.placements_x:
            # DTensor's backward may spuriously all_gather to Replicate();
            # redistribute back to the input's placements.
            grad_output = grad_output.redistribute(ctx.device_mesh, ctx.placements_x)

        grad_output_local = grad_output.to_local()
        all_replicate = tuple([Replicate()] * ctx.device_mesh.ndim)

        # Compute dtype (e.g. bf16 under mixed precision).  We cast operands
        # to this dtype explicitly for CPU compatibility — on CUDA, custom_bwd
        # restores autocast which handles this automatically.
        go_dtype = grad_output_local.dtype

        if ctx.needs_input_grad[1]:
            dw_local = torch.einsum("...i,...o->io", grad_output_local, x_local.to(go_dtype))
            dw_local = _LinearParamsReplicatedImpl._all_reduce_grad_gteqfp32(
                dw_local, ctx.device_mesh, ctx.placements_grad_params, ctx.dtype_weight
            )
            dw = DTensor.from_local(
                dw_local, ctx.device_mesh, all_replicate, shape=ctx.weight_shape, stride=ctx.weight_stride
            )
        else:
            dw = None

        if ctx.needs_input_grad[2]:
            dims = list(range(grad_output_local.ndim - 1))
            if ctx.dtype_bias is None:
                raise RuntimeError("bias gradient requested but bias dtype metadata is missing")
            db_local = grad_output_local.sum(dim=dims)
            db_local = _LinearParamsReplicatedImpl._all_reduce_grad_gteqfp32(
                db_local, ctx.device_mesh, ctx.placements_grad_params, ctx.dtype_bias
            )
            db = DTensor.from_local(
                db_local, ctx.device_mesh, all_replicate, shape=ctx.bias_shape, stride=ctx.bias_stride
            )
        else:
            db = None

        if ctx.needs_input_grad[0]:
            if weight_local is None:
                raise RuntimeError("input gradient requested but saved weight tensor is missing")
            grad_input_local = torch.einsum("...i,io->...o", grad_output_local, weight_local.to(go_dtype))
            grad_input = DTensor.from_local(
                grad_input_local, ctx.device_mesh, ctx.placements_x, shape=ctx.shape_input, stride=ctx.stride_input
            )
        else:
            grad_input = None

        return grad_input, dw, db, None, None


class LinearParamsReplicated(nn.Module):
    """
    Distributed linear layer with parameters replicated across all device mesh dimensions.

    This is almost equivalent to
    ```python
    layer = torch.distributed.tensor.distribute_module(layer_local, device_mesh)
    ```
    with the exception that the torch.distributed.tensor.distribute_module version will incur
    significant overhead due to the unnecessary replication of the output tensor along certain
    device mesh dimensions.

    This class avoids such unnecessary overhead by using the custom _LinearParamsReplicatedImpl
    autograd function for forward and backward pass computation instead of relying on the distributed
    module's forward implementation.

    NOTE: by default, avg reduce over the Replicate placements of the weight and bias gradients
    is performed. This is to ensure identical parameter updates across all ranks and avoid
    gradual divergence during training. This can be disabled by setting
    avg_over_replicate_param_grad to False.

    Key requirements:
        1. Parameters (weight and bias) will replicated on all device mesh dimensions
        2. Input tensor and parameters must be on the same device mesh
        3. Feature/hidden dimension of the input must not be sharded across the device mesh
        4. Partial reduction along any input dimension is not supported
        5. Input and outputs must be on the same device mesh with the same placements
        6. Gradients of the input have the same placements on the same device mesh as the input
        7. Gradients of the weight and bias have Partial("sum") placements along the input's Shard placements'
           dimension so that the all-reduce will be performed along those device-grid dimensions

    """

    def __init__(self, layer_local: nn.Linear, device_mesh: DeviceMesh, avg_over_replicate_param_grad: bool = True):
        """
        Initialize the distributed linear layer.

        Args:
            layer_local: nn.Linear to be distributed
            device_mesh: Device mesh for distributed computation
            avg_over_replicate_param_grad: whether to perform avg reduce over the
                Replicate placements of the weight and bias gradients. For example,
                if the input DTensor x.placements = (Shard(0), Replicate()), this layer's
                parameters' gradients.placements = (Partial("sum"), Replicate()) if
                self._avg_over_replicate_param_grad is False; otherwise, it will be
                (Partial("sum"), Partial("avg")). The motivation is to ensure identical
                parameter updates across all ranks and avoid gradual divergence during
                training.
        """
        if not isinstance(layer_local, nn.Linear):
            raise ValueError("layer_local is not an instance of nn.Linear")
        if layer_local.weight.device.type != device_mesh.device_type:
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
        all_replicate_placements = [Replicate()] * device_mesh.ndim
        self.weight = nn.Parameter(
            distribute_tensor(layer_local.weight.data, device_mesh, all_replicate_placements),
            requires_grad=layer_local.weight.requires_grad,
        )
        if layer_local.bias is None:
            self.register_parameter("bias", None)
        else:
            self.bias = nn.Parameter(
                distribute_tensor(layer_local.bias.data, device_mesh, all_replicate_placements),
                requires_grad=layer_local.bias.requires_grad,
            )
        self._avg_over_replicate_param_grad = avg_over_replicate_param_grad

    def forward(self, input: DTensor) -> DTensor:
        """
        Forward pass for the distributed linear layer.

        Uses the custom _LinearParamsReplicatedImpl autograd function to perform the computation
        efficiently while preserving correct autograd behavior for distributed tensors.

        Args:
            input: Input DTensor with appropriate placement strategy

        Returns:
            Output DTensor with same placement strategy as input
        """
        return _LinearParamsReplicatedImpl.apply(
            input,
            self.weight,
            self.bias,
            True,  # cast_params_dtype_to_x: under bf16-mixed autocast, upstream
            # ops produce bf16 activations while weights stay fp32.  custom_fwd
            # disables autocast inside the function, so F.linear would get
            # mismatched dtypes.  Casting weight to input dtype matches what
            # native autocast does for F.linear.  No-op when dtypes already match.
            self._avg_over_replicate_param_grad,
        )


class ContextParallelLinear(nn.Linear):
    def __init__(self, in_features, out_features, reduce_group, bias=True):
        """Context parallel linear layer, a wrapper around nn.Linear that supports distributed training.

        Although the output between tokens is independent, the backward pass on the weight and bias tensors involves a summation over all tokens, and thus requires all-reduce due to context parallelism. This means we need a dedicated ContextParallelLinear class for training.

        Args:
            in_features: The number of input features.
            out_features: The number of output features.
            reduce_group: The process group to use for the all-reduce in backward pass.
            bias: Whether to use a bias.

        If group is not provided, the layer will behave like a normal nn.Linear.
        """
        super().__init__(in_features, out_features, bias)
        assert reduce_group is not None, "reduce_group must be provided"
        self.reduce_group = reduce_group

    def forward(self, input: Tensor) -> Tensor:
        return _ContextParallelLinearImpl.apply(input, self.weight, self.bias, self.reduce_group)


def get_cp_linear(
    *args, reduce_group: Optional[dist.ProcessGroup] = None, **kwargs
) -> nn.Linear | ContextParallelLinear:
    """Get a context parallel linear layer.

    If group is not provided, the returned layer will fall back to nn.Linear.
    """
    if reduce_group is None:
        return nn.Linear(*args, **kwargs)
    return ContextParallelLinear(*args, reduce_group=reduce_group, **kwargs)
