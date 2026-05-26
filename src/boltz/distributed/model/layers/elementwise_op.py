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


from enum import Enum, auto

import torch
from torch.distributed.tensor import DTensor, Partial, Shard


class ElementwiseOp(Enum):
    """Enumeration of supported elementwise operations."""

    # n-ary ops
    SUM = auto()
    SUB = auto()
    PROD = auto()
    DIV = auto()
    EQUAL = auto()
    BITAND = auto()

    # unary ops
    COS = auto()
    RELU = auto()
    ROUND = auto()
    EXP = auto()
    ABS = auto()
    SIGMOID = auto()

    # comparison ops
    GT = auto()
    LT = auto()
    LOG = auto()

    # scalar-tensor ops
    POW = auto()
    MAX = auto()


class _SingleTensorOpImpl(torch.autograd.Function):
    """Distributed implementation of single-tensor operations using DTensors.

    This autograd function implements distributed single-tensor operations
    like cosine, ReLU, and logarithm. The operations are performed element-wise across distributed
    tensors while maintaining proper gradient computation.

    Supported operations:
    - COS: output = cos(x)
    - RELU: output = max(0, x)
    - ROUND: output = round(x)
    - LOG: output = log(x)
    - EXP: output = exp(x)
    - ABS: output = |x|
    - SIGMOID: output = 1 / (1 + exp(-x))

    Key features:
    - Distributed computation across device meshes with various sharding strategies
    - Memory-efficient implementation that operates on local tensor chunks
    - Supports gradient computation through custom backward pass
    """

    @staticmethod
    def forward(ctx, x: DTensor, op: ElementwiseOp) -> DTensor:
        """Forward pass of distributed single-tensor operation.

        Parameters
        ----------
        ctx : torch.autograd.function.BackwardCFrame
            Context object for saving information needed in backward pass.
        x : DTensor
            Input tensor. Can have any shape and sharding strategy.
        op : ElementwiseOp
            The operation to perform (COS, RELU, ROUND, LOG, EXP, ABS, or SIGMOID).
        Returns
        -------
        DTensor
            Output tensor with shape identical to input tensor.

        Raises
        ------
        TypeError
            If input is not a DTensor.
        ValueError
            If Partial placements are used (not supported), or if op is invalid.
        """
        if not isinstance(x, DTensor):
            raise TypeError(f"Input 'x' must be of type DTensor. Got type {type(x)}.")
        if not isinstance(op, ElementwiseOp):
            raise TypeError(f"Input 'op' must be of type ElementwiseOp. Got type {type(op)}.")

        device_mesh_input = x.device_mesh
        placements_input = x.placements

        for i_dim_device_mesh, placement in enumerate(placements_input):
            if isinstance(placement, Partial):
                raise ValueError("Partial placements are not supported")
            elif isinstance(placement, Shard):
                if x.shape[placement.dim] % device_mesh_input.shape[i_dim_device_mesh] != 0:
                    raise ValueError(
                        f"Uneven sharding tensor dimension {placement.dim} of size {x.shape[placement.dim]} "
                        f"along device mesh dimension {i_dim_device_mesh} of size {device_mesh_input.shape[i_dim_device_mesh]} is not supported"
                    )

        x_local = x.to_local()

        # Perform the operation
        if op == ElementwiseOp.COS:
            output_local = torch.cos(x_local)
            if x.requires_grad:
                x_local_copy = x_local.detach().clone()
                ctx.save_for_backward(x_local_copy)
        elif op == ElementwiseOp.RELU:
            output_local = torch.relu(x_local)
            if x.requires_grad:
                x_local_copy = x_local.detach().clone()
                ctx.save_for_backward(x_local_copy)
        elif op == ElementwiseOp.ROUND:
            output_local = torch.round(x_local)
        elif op == ElementwiseOp.LOG:
            output_local = torch.log(x_local)
            if x.requires_grad:
                x_local_copy = x_local.detach().clone()
                ctx.save_for_backward(x_local_copy)
        elif op == ElementwiseOp.EXP:
            output_local = torch.exp(x_local)
            if x.requires_grad:
                x_local_copy = x_local.detach().clone()
                ctx.save_for_backward(x_local_copy)
        elif op == ElementwiseOp.ABS:
            output_local = torch.abs(x_local)
            if x.requires_grad:
                x_local_copy = x_local.detach().clone()
                ctx.save_for_backward(x_local_copy)
        elif op == ElementwiseOp.SIGMOID:
            output_local = torch.sigmoid(x_local)
            if x.requires_grad:
                ctx.save_for_backward(output_local.clone())
        else:
            raise ValueError(f"Unsupported single-tensor operation: {op}")

        ctx.device_mesh_input = device_mesh_input
        ctx.placements_input = placements_input
        ctx.input_shape = x.shape
        ctx.op = op

        out = DTensor.from_local(
            output_local,
            device_mesh=device_mesh_input,
            placements=placements_input,
            shape=x.shape,
            stride=x.stride(),
        )
        if op == ElementwiseOp.ROUND:
            ctx.mark_non_differentiable(out)
        return out

    @staticmethod
    def backward(ctx, grad_output: DTensor) -> tuple[DTensor | None, None]:
        """Backward pass of distributed single-tensor operation.

        Parameters
        ----------
        ctx : torch.autograd.function.BackwardCFrame
            Context object containing saved tensors and metadata from forward pass.
        grad_output : DTensor
            Gradients of the loss with respect to the output tensor.

        Returns
        -------
        tuple[DTensor | None, None]
            Gradients with respect to x and op.
        """
        if not isinstance(grad_output, DTensor):
            raise TypeError(f"Input 'grad_output' must be of type DTensor. Got type {type(grad_output)}.")

        if grad_output.device_mesh != ctx.device_mesh_input:
            raise ValueError(
                f"Input 'grad_output' must have the same device mesh as the input tensor. "
                f"Got device meshes {grad_output.device_mesh} and {ctx.device_mesh_input}."
            )

        if grad_output.placements != ctx.placements_input:
            raise ValueError(
                f"Input 'grad_output' must have the same placements as the input tensor. "
                f"Got placements {grad_output.placements} and {ctx.placements_input}."
            )

        grad_output_local = grad_output.to_local()
        dx = None
        x_local = ctx.saved_tensors[0]

        # Compute gradients based on operation
        if ctx.op == ElementwiseOp.COS:
            if ctx.needs_input_grad[0]:
                # Derivative of cos(x) is -sin(x)
                dx_local = grad_output_local * (-torch.sin(x_local))
                dx = DTensor.from_local(
                    dx_local,
                    device_mesh=ctx.device_mesh_input,
                    placements=ctx.placements_input,
                    shape=grad_output.shape,
                    stride=grad_output.stride(),
                )
        elif ctx.op == ElementwiseOp.RELU:
            if ctx.needs_input_grad[0]:
                # Derivative of relu(x) is 1 where x > 0, 0 elsewhere
                dx_local = grad_output_local.clone()
                dx_local[x_local <= 0] = 0
                dx = DTensor.from_local(
                    dx_local,
                    device_mesh=ctx.device_mesh_input,
                    placements=ctx.placements_input,
                    shape=grad_output.shape,
                    stride=grad_output.stride(),
                )
        elif ctx.op == ElementwiseOp.ROUND:
            pass  # no gradient through this op
        elif ctx.op == ElementwiseOp.LOG:
            if ctx.needs_input_grad[0]:
                # Derivative of log(x) is 1/x
                dx_local = grad_output_local / x_local
                dx = DTensor.from_local(
                    dx_local,
                    device_mesh=ctx.device_mesh_input,
                    placements=ctx.placements_input,
                    shape=grad_output.shape,
                    stride=grad_output.stride(),
                )
        elif ctx.op == ElementwiseOp.EXP:
            if ctx.needs_input_grad[0]:
                # Derivative of exp(x) is exp(x)
                dx_local = grad_output_local * torch.exp(x_local)
                dx = DTensor.from_local(
                    dx_local,
                    device_mesh=ctx.device_mesh_input,
                    placements=ctx.placements_input,
                    shape=grad_output.shape,
                    stride=grad_output.stride(),
                )
        elif ctx.op == ElementwiseOp.ABS:
            if ctx.needs_input_grad[0]:
                # Derivative of abs(x) is sign(x) = x/abs(x) for x != 0, undefined at x = 0
                # We use torch.sign(x) which handles the case x = 0 by returning 0
                dx_local = grad_output_local * torch.sign(x_local)
                dx = DTensor.from_local(
                    dx_local,
                    device_mesh=ctx.device_mesh_input,
                    placements=ctx.placements_input,
                    shape=grad_output.shape,
                    stride=grad_output.stride(),
                )
        elif ctx.op == ElementwiseOp.SIGMOID:
            if ctx.needs_input_grad[0]:
                # Derivative of sigmoid(x) is sigmoid(x) * (1 - sigmoid(x))
                sigmoid_output = x_local
                dx_local = grad_output_local * sigmoid_output * (1 - sigmoid_output)
                dx = DTensor.from_local(
                    dx_local,
                    device_mesh=ctx.device_mesh_input,
                    placements=ctx.placements_input,
                    shape=grad_output.shape,
                    stride=grad_output.stride(),
                )
        else:
            raise ValueError(f"Unsupported single-tensor operation: {ctx.op}")

        return dx, None


class _ElementwiseOpImpl(torch.autograd.Function):
    """Distributed implementation of elementwise operations using DTensors.

    This autograd function implements distributed elementwise operations that can perform
    summation, multiplication, division, and logical operations between two input tensors. The operations are
    performed element-wise across distributed tensors while maintaining proper gradient
    computation.

    Supported operations:
    - SUM: output = a + b
    - SUB: output = a - b
    - PROD: output = a * b
    - DIV: output = a / b
    - EQUAL: output = a & b
    - BITAND: output = a & b

    Key features:
    - Distributed computation across device meshes with various sharding strategies
    - Memory-efficient implementation that operates on local tensor chunks
    - Supports gradient computation through custom backward pass
    - Validates tensor compatibility (device mesh, placements, shapes)

    Notes
    -----
    Input tensors must be DTensors with:
    - Identical device mesh and placements
    - Compatible shapes (a and b must have the same shape)
    - No Partial placements (not currently supported)
    """

    @staticmethod
    def forward(ctx, a: DTensor, b: DTensor, op: ElementwiseOp) -> DTensor:
        """Forward pass of distributed elementwise operation.

        Parameters
        ----------
        ctx : torch.autograd.function.BackwardCFrame
            Context object for saving information needed in backward pass.
        a : DTensor
            First input tensor. Can have any shape and sharding strategy.
        b : DTensor
            Second input tensor. Must have identical shape, device mesh,
            and placements as a.
        op : ElementwiseOp
            The operation to perform (SUM, PROD, DIV, EQUAL, or BITAND).

        Returns
        -------
        DTensor
            Output tensor with shape identical to input tensors.
            Contains the result of the specified operation.

        Raises
        ------
        TypeError
            If inputs are not DTensors.
        ValueError
            If tensors have incompatible device meshes, placements, or if
            Partial placements are used (not supported), or if op is invalid.
        """
        if not isinstance(a, DTensor):
            raise TypeError(f"Input 'a' must be of type DTensor. Got type {type(a)}.")
        if not isinstance(b, DTensor):
            raise TypeError(f"Input 'b' must be of type DTensor. Got type {type(b)}.")
        if not isinstance(op, ElementwiseOp):
            raise TypeError(f"Input 'op' must be of type ElementwiseOp. Got type {type(op)}.")

        device_mesh_input = a.device_mesh
        if b.device_mesh != device_mesh_input:
            raise ValueError(
                f"Input tensors 'a' and 'b' must have identical device mesh. "
                f"Got device meshes {device_mesh_input} and {b.device_mesh}."
            )

        placements_input = a.placements
        for i_dim_device_mesh, placement in enumerate(placements_input):
            if isinstance(placement, Partial):
                raise ValueError("Partial placements are not supported")
            elif isinstance(placement, Shard):
                if a.shape[placement.dim] % device_mesh_input.shape[i_dim_device_mesh] != 0:
                    raise ValueError(
                        f"Uneven sharding tensor dimension {placement.dim} of size {a.shape[placement.dim]} "
                        f"along device mesh dimension {i_dim_device_mesh} of size {device_mesh_input.shape[i_dim_device_mesh]} is not supported"
                    )

        if b.placements != placements_input:
            raise ValueError(
                f"Input tensors 'a' and 'b' must have identical placements. "
                f"Got placements {placements_input} and {b.placements}."
            )

        input_shape = a.shape
        if input_shape != b.shape:
            raise ValueError(
                f"Input tensors 'a' and 'b' must have identical shapes. Got shapes {input_shape} and {b.shape}."
            )

        a_local = a.to_local()
        b_local = b.to_local()

        # Perform the operation
        if op == ElementwiseOp.SUM:
            output_local = a_local + b_local
        elif op == ElementwiseOp.SUB:
            output_local = a_local - b_local
        elif op == ElementwiseOp.PROD:
            output_local = a_local * b_local
            # TODO: check if we can afford save_for_backward(a_local, b_local) without explicitly copying
            # pytorch's c++ backend has this code here that can determine the necessity of the copy:
            # https://github.com/pytorch/pytorch/blob/7caf6c801ddfaf556a3ca191173b50002c4261f4/torch/csrc/autograd/saved_variable.cpp#L67-L79
            # so we might not need to explicitly copy the tensors here
            a_local_copy = a_local.detach().clone() if b.requires_grad else None
            b_local_copy = b_local.detach().clone() if a.requires_grad else None
            ctx.save_for_backward(a_local_copy, b_local_copy)
        elif op == ElementwiseOp.DIV:
            output_local = a_local / b_local
            # Save tensors for backward pass:
            # - if a.requires_grad, we need b for: da = grad_output / b
            # - if b.requires_grad, we need b and output for: db = -grad_output * output / b
            b_local_copy = b_local.detach().clone() if a.requires_grad else None
            output_over_b_local = output_local / b_local if b.requires_grad else None
            ctx.save_for_backward(b_local_copy, output_over_b_local)
        elif op == ElementwiseOp.EQUAL:
            if a.requires_grad or b.requires_grad:
                raise ValueError("EQUAL operation is not differentiable but requires_grad is True")
            output_local = a_local & b_local
        elif op == ElementwiseOp.BITAND:
            if a.requires_grad or b.requires_grad:
                raise ValueError("BITAND operation is not differentiable but requires_grad is True")
            output_local = a_local & b_local
        else:
            raise ValueError(f"Unsupported operation: {op}")

        ctx.device_mesh_input = device_mesh_input
        ctx.placements_input = placements_input
        ctx.input_shape = input_shape
        ctx.op = op

        out = DTensor.from_local(
            output_local,
            device_mesh=device_mesh_input,
            placements=placements_input,
            shape=a.shape,
            stride=a.stride(),
        )
        return out

    @staticmethod
    def backward(ctx, grad_output) -> tuple[DTensor | None, DTensor | None, None]:
        """Backward pass of distributed elementwise operation.

        Computes gradients with respect to both input tensors a and b.

        The gradients are:
        - For SUM: da = grad_output, db = grad_output
        - For SUB: da = grad_output, db = -grad_output
        - For PROD: da = grad_output * b, db = grad_output * a
        - For DIV: da = grad_output / b, db = -grad_output * output / b
        - For EQUAL: da = 0, db = 0 (not differentiable)

        Parameters
        ----------
        ctx : torch.autograd.function.BackwardCFrame
            Context object containing saved tensors and metadata from forward pass.
        grad_output : DTensor
            Gradient of the loss with respect to the output tensor.
            Must have identical device mesh and placements as the input tensors.

        Returns
        -------
        tuple[DTensor, DTensor, None]
            Gradients with respect to a, b, and None for the op parameter.
            Both gradients have the same shape and distribution as their corresponding inputs.

        Raises
        ------
        TypeError
            If grad_output is not a DTensor.
        ValueError
            If grad_output has incompatible device mesh or placements compared
            to the input tensors from the forward pass.
        """

        if not isinstance(grad_output, DTensor):
            raise TypeError(f"Input 'grad_output' must be of type DTensor. Got type {type(grad_output)}.")

        if grad_output.device_mesh != ctx.device_mesh_input:
            raise ValueError(
                f"Input 'grad_output' must have the same device mesh as the input tensor. "
                f"Got device meshes {grad_output.device_mesh} and {ctx.device_mesh_input}."
            )

        if grad_output.placements != ctx.placements_input:
            raise ValueError(
                f"Input 'grad_output' must have the same placements as the input tensor. "
                f"Got placements {grad_output.placements} and {ctx.placements_input}."
            )

        if grad_output.shape != ctx.input_shape:
            raise ValueError(
                f"Input 'grad_output' must have the same shape as the input tensor. "
                f"Got shapes {grad_output.shape} and {ctx.input_shape}."
            )

        grad_output_local = grad_output.to_local()

        # Compute gradients based on operation
        da = None
        db = None

        if ctx.op == ElementwiseOp.SUM:
            if ctx.needs_input_grad[0]:
                da_local = grad_output_local
                da = DTensor.from_local(
                    da_local,
                    device_mesh=ctx.device_mesh_input,
                    placements=ctx.placements_input,
                    shape=grad_output.shape,
                    stride=grad_output.stride(),
                )
            if ctx.needs_input_grad[1]:
                db_local = grad_output_local
                db = DTensor.from_local(
                    db_local,
                    device_mesh=ctx.device_mesh_input,
                    placements=ctx.placements_input,
                    shape=grad_output.shape,
                    stride=grad_output.stride(),
                )
        elif ctx.op == ElementwiseOp.SUB:
            if ctx.needs_input_grad[0]:
                da_local = grad_output_local
                da = DTensor.from_local(
                    da_local,
                    device_mesh=ctx.device_mesh_input,
                    placements=ctx.placements_input,
                    shape=grad_output.shape,
                    stride=grad_output.stride(),
                )
            if ctx.needs_input_grad[1]:
                db_local = -grad_output_local
                db = DTensor.from_local(
                    db_local,
                    device_mesh=ctx.device_mesh_input,
                    placements=ctx.placements_input,
                    shape=grad_output.shape,
                    stride=grad_output.stride(),
                )
        elif ctx.op == ElementwiseOp.PROD:
            # Unpack saved_tensors once for checkpoint compatibility - must only be done once
            a_local, b_local = ctx.saved_tensors

            if ctx.needs_input_grad[0]:
                da_local = grad_output_local * b_local
                da = DTensor.from_local(
                    da_local,
                    device_mesh=ctx.device_mesh_input,
                    placements=ctx.placements_input,
                    shape=grad_output.shape,
                    stride=grad_output.stride(),
                )
            if ctx.needs_input_grad[1]:
                db_local = grad_output_local * a_local
                db = DTensor.from_local(
                    db_local,
                    device_mesh=ctx.device_mesh_input,
                    placements=ctx.placements_input,
                    shape=grad_output.shape,
                    stride=grad_output.stride(),
                )
        elif ctx.op == ElementwiseOp.DIV:
            # Unpack saved_tensors once for checkpoint compatibility - must only be done once
            b_local, output_over_b_local = ctx.saved_tensors

            if ctx.needs_input_grad[0]:
                # Gradient w.r.t. numerator: da = grad_output / b
                da_local = grad_output_local / b_local
                da = DTensor.from_local(
                    da_local,
                    device_mesh=ctx.device_mesh_input,
                    placements=ctx.placements_input,
                    shape=grad_output.shape,
                    stride=grad_output.stride(),
                )
            if ctx.needs_input_grad[1]:
                # Gradient w.r.t. denominator: db = -grad_output * output / b
                db_local = -grad_output_local * output_over_b_local
                db = DTensor.from_local(
                    db_local,
                    device_mesh=ctx.device_mesh_input,
                    placements=ctx.placements_input,
                    shape=grad_output.shape,
                    stride=grad_output.stride(),
                )
        elif ctx.op == ElementwiseOp.EQUAL:  # dummy op, not differentiable
            pass
        elif ctx.op == ElementwiseOp.BITAND:  # dummy op, not differentiable
            pass
        else:
            raise ValueError(f"Unsupported operation: {ctx.op}")

        return da, db, None


def elementwise_op(a: DTensor, b: DTensor, op: ElementwiseOp) -> DTensor:
    """Apply elementwise operation to two distributed tensors.

    This function performs element-wise operations between two distributed tensors.
    Supported operations are summation (a + b), subtraction (a - b), multiplication (a * b),
    division (a / b), logical AND (a & b), and bitwise AND (a & b). The operation is performed
    efficiently using local tensor operations while maintaining gradient
    computation capabilities.

    Parameters
    ----------
    a : DTensor
        First input tensor. Can have any shape and sharding strategy.
    b : DTensor
        Second input tensor. Must have identical shape, device mesh,
        and placements as a.
    op : ElementwiseOp
        The operation to perform (ElementwiseOp.SUM, ElementwiseOp.SUB, ElementwiseOp.PROD, ElementwiseOp.DIV, ElementwiseOp.EQUAL, or ElementwiseOp.BITAND).

    Returns
    -------
    DTensor
        Output tensor with shape identical to input tensors.
        Contains the result of the specified operation.

    Examples
    --------
    >>> # Assume we have distributed tensors a and b with shape (B, N, D)
    >>> sum_output = elementwise_op(a, b, ElementwiseOp.SUM)
    >>> # sum_output = a + b, computed in distributed fashion
    >>>
    >>> sub_output = elementwise_op(a, b, ElementwiseOp.SUB)
    >>> # sub_output = a - b, computed in distributed fashion
    >>>
    >>> prod_output = elementwise_op(a, b, ElementwiseOp.PROD)
    >>> # prod_output = a * b, computed in distributed fashion
    >>>
    >>> div_output = elementwise_op(a, b, ElementwiseOp.DIV)
    >>> # div_output = a / b, computed in distributed fashion
    >>>
    >>> equal_output = elementwise_op(a, b, ElementwiseOp.EQUAL)
    >>> # equal_output = a == b, computed in distributed fashion
    >>>
    >>> bitand_output = elementwise_op(a, b, ElementwiseOp.BITAND)
    >>> # bitand_output = a & b, computed in distributed fashion

    Notes
    -----
    - Both input tensors must be DTensors with compatible device meshes and placements
    - Partial placements are not currently supported
    - The function is differentiable and supports gradient computation
    - The operation is performed on local tensor chunks for efficiency
    """
    return _ElementwiseOpImpl.apply(a, b, op)  # type: ignore


def single_tensor_op(x: DTensor, op: ElementwiseOp) -> DTensor:
    """Apply single-tensor operation to a distributed tensor.

    This function performs element-wise operations on a single distributed tensor.
    Supports cosine, ReLU, round, logarithm, and exponential operations.

    Parameters
    ----------
    x : DTensor
        Input tensor. Can have any shape and sharding strategy.
    op : ElementwiseOp
        The operation to perform (ElementwiseOp.COS, ElementwiseOp.RELU, ElementwiseOp.ROUND, ElementwiseOp.LOG, ElementwiseOp.EXP, ElementwiseOp.ABS, or ElementwiseOp.SIGMOID).

    Returns
    -------
    DTensor
        Output tensor with shape identical to input tensor.

    Examples
    --------
    >>> # Assume we have distributed tensor x with shape (B, N, D)
    >>> cos_output = single_tensor_op(x, ElementwiseOp.COS)
    >>> # cos_output = cos(x), computed in distributed fashion
    >>>
    >>> relu_output = single_tensor_op(x, ElementwiseOp.RELU)
    >>> # relu_output = max(0, x), computed in distributed fashion
    >>>
    >>> round_output = single_tensor_op(x, ElementwiseOp.ROUND)
    >>> # round_output = round(x), computed in distributed fashion
    >>>
    >>> log_output = single_tensor_op(x, ElementwiseOp.LOG)
    >>> # log_output = log(x), computed in distributed fashion
    >>>
    >>> exp_output = single_tensor_op(x, ElementwiseOp.EXP)
    >>> # exp_output = exp(x), computed in distributed fashion
    >>>
    >>> abs_output = single_tensor_op(x, ElementwiseOp.ABS)
    >>> # abs_output = |x|, computed in distributed fashion
    >>>
    >>> sigmoid_output = single_tensor_op(x, ElementwiseOp.SIGMOID)
    >>> # sigmoid_output = 1 / (1 + exp(-x)), computed in distributed fashion

    Notes
    -----
    - Input tensor must be a DTensor with any placement strategy
    - Partial placements are not currently supported
    - The function is differentiable and supports gradient computation
    - The operation is performed on local tensor chunks for efficiency
    """
    return _SingleTensorOpImpl.apply(x, op)  # type: ignore


class _ScalarTensorOpImpl(torch.autograd.Function):
    """Distributed implementation of scalar-tensor operations using DTensors.

    This autograd function implements distributed operations between a scalar and a DTensor.
    The operations are performed element-wise across distributed tensors while maintaining
    proper gradient computation.

    Supported operations:
    - SUM: output = scalar + tensor
    - SUB: output = scalar - tensor
    - PROD: output = scalar * tensor
    - DIV: output = scalar / tensor
    - GT: output = scalar > tensor
    - LT: output = scalar < tensor
    - EQUAL: output = scalar == tensor
    - POW: output = tensor ** scalar
    - MAX: output = max(scalar, tensor)  (element-wise clamp from below)

    Key features:
    - Distributed computation across device meshes with various sharding strategies
    - Memory-efficient implementation that operates on local tensor chunks
    - Supports gradient computation through custom backward pass
    - Validates tensor compatibility (no Partial placements)
    """

    @staticmethod
    def forward(ctx, scalar: float | int, tensor: DTensor, op: ElementwiseOp) -> DTensor:
        """Forward pass of distributed scalar-tensor operation.

        Parameters
        ----------
        ctx : torch.autograd.function.BackwardCFrame
            Context object for saving information needed in backward pass.
        scalar : float
            Scalar value to operate with.
        tensor : DTensor
            Input tensor. Can have any shape and sharding strategy.
        op : ElementwiseOp
            The operation to perform (SUM, SUB, PROD, DIV, GT, or LT).

        Returns
        -------
        DTensor
            Output tensor with shape identical to input tensor.

        Raises
        ------
        TypeError
            If inputs are not of expected types.
        ValueError
            If Partial placements are used (not supported), or if op is invalid.
        """
        if not isinstance(scalar, (int, float)):
            raise TypeError(f"Input 'scalar' must be of type int or float. Got type {type(scalar)}.")
        if not isinstance(tensor, DTensor):
            raise TypeError(f"Input 'tensor' must be of type DTensor. Got type {type(tensor)}.")
        if not isinstance(op, ElementwiseOp):
            raise TypeError(f"Input 'op' must be of type ElementwiseOp. Got type {type(op)}.")

        device_mesh_input = tensor.device_mesh
        placements_input = tensor.placements

        for i_dim_device_mesh, placement in enumerate(placements_input):
            if isinstance(placement, Partial):
                raise ValueError("Partial placements are not supported")
            elif isinstance(placement, Shard):
                if tensor.shape[placement.dim] % device_mesh_input.shape[i_dim_device_mesh] != 0:
                    raise ValueError(
                        f"Uneven sharding tensor dimension {placement.dim} of size {tensor.shape[placement.dim]} "
                        f"along device mesh dimension {i_dim_device_mesh} of size {device_mesh_input.shape[i_dim_device_mesh]} is not supported"
                    )

        tensor_local = tensor.to_local()

        # Perform the operation
        if op == ElementwiseOp.SUM:
            output_local = scalar + tensor_local
        elif op == ElementwiseOp.SUB:
            output_local = scalar - tensor_local
        elif op == ElementwiseOp.PROD:
            output_local = scalar * tensor_local
        elif op == ElementwiseOp.DIV:
            output_local = scalar / tensor_local
            if tensor.requires_grad:
                # Save tensor for backward pass
                tensor_local_copy = tensor_local.detach().clone()
                ctx.save_for_backward(tensor_local_copy)
        elif op == ElementwiseOp.GT:
            output_local = scalar > tensor_local
        elif op == ElementwiseOp.LT:
            output_local = scalar < tensor_local
        elif op == ElementwiseOp.EQUAL:
            output_local = torch.eq(
                torch.tensor(scalar, device=tensor_local.device, dtype=tensor_local.dtype), tensor_local
            )
        elif op == ElementwiseOp.POW:
            if tensor_local.min() < 0 and not scalar.is_integer():
                raise ValueError(
                    "Negative tensor values are not supported for DTensor POW operation but got scalar: {scalar}"
                )
            output_local = torch.pow(tensor_local, scalar)
            if tensor.requires_grad:
                # Save tensor for backward pass: d_tensor = grad_output * scalar * tensor^(scalar-1)
                if scalar != 0:  # 0 pow has gradient of zero
                    tensor_local_copy = tensor_local.detach().clone()
                    ctx.save_for_backward(tensor_local_copy)
        elif op == ElementwiseOp.MAX:
            # max(scalar, tensor) element-wise, equivalent to tensor.clamp(min=scalar)
            output_local = torch.clamp(tensor_local, min=scalar)
            if tensor.requires_grad:
                # Save mask for backward: gradient passes through where tensor >= scalar
                mask_local = (tensor_local >= scalar).detach()
                ctx.save_for_backward(mask_local)
        else:
            raise ValueError(f"Unsupported scalar-tensor operation: {op}")

        if tensor.requires_grad:
            ctx.device_mesh_input = device_mesh_input
            ctx.placements_input = placements_input
            ctx.input_shape = tensor.shape
            ctx.op = op
            ctx.scalar = scalar

        out = DTensor.from_local(
            output_local,
            device_mesh=device_mesh_input,
            placements=placements_input,
            shape=tensor.shape,
            stride=tensor.stride(),
        )
        if op == ElementwiseOp.GT or op == ElementwiseOp.LT or op == ElementwiseOp.EQUAL:
            ctx.mark_non_differentiable(out)
        return out

    @staticmethod
    def backward(ctx, grad_output) -> tuple[float | None, DTensor | None, None]:
        """Backward pass of distributed scalar-tensor operation.

        Computes gradients with respect to both scalar and tensor inputs.

        The gradients are:
        - For SUM: d_scalar = grad_output.sum(), d_tensor = grad_output
        - For SUB: d_scalar = grad_output.sum(), d_tensor = -grad_output
        - For PROD: d_scalar = (grad_output * tensor).sum(), d_tensor = grad_output * scalar
        - For DIV: d_scalar = -(grad_output * tensor / scalar^2).sum(), d_tensor = -grad_output * scalar / tensor^2
        - For GT: d_scalar = None, d_tensor = None (not differentiable)
        - For LT: d_scalar = None, d_tensor = None (not differentiable)
        - For EQUAL: d_scalar = None, d_tensor = None (not differentiable)

        Parameters
        ----------
        ctx : torch.autograd.function.BackwardCFrame
            Context object containing saved tensors and metadata from forward pass.
        grad_output : DTensor
            Gradients of the loss with respect to the output tensors.

        Returns
        -------
        tuple[float | None, DTensor | None, None]
            Gradients with respect to scalar, tensor, and None for the op parameter.
        """
        if not ctx.needs_input_grad[1]:
            return None, None, None

        if not isinstance(grad_output, DTensor):
            raise TypeError(f"Input 'grad_output' must be of type DTensor. Got type {type(grad_output)}.")

        if grad_output.device_mesh != ctx.device_mesh_input:
            raise ValueError(
                f"Input 'grad_output' must have the same device mesh as the input tensor. "
                f"Got device meshes {grad_output.device_mesh} and {ctx.device_mesh_input}."
            )

        if grad_output.placements != ctx.placements_input:
            raise ValueError(
                f"Input 'grad_output' must have the same placements as the input tensor. "
                f"Got placements {grad_output.placements} and {ctx.placements_input}."
            )

        grad_output_local = grad_output.to_local()
        d_tensor = None

        # Compute gradients based on operation
        if ctx.op == ElementwiseOp.SUM:
            # Gradient w.r.t. tensor: grad_output
            d_tensor = DTensor.from_local(
                grad_output_local,
                device_mesh=ctx.device_mesh_input,
                placements=ctx.placements_input,
                shape=grad_output.shape,
                stride=grad_output.stride(),
            )
        elif ctx.op == ElementwiseOp.SUB:
            # Gradient w.r.t. tensor: -grad_output
            d_tensor_local = -grad_output_local
            d_tensor = DTensor.from_local(
                d_tensor_local,
                device_mesh=ctx.device_mesh_input,
                placements=ctx.placements_input,
                shape=grad_output.shape,
                stride=grad_output.stride(),
            )
        elif ctx.op == ElementwiseOp.PROD:
            # Gradient w.r.t. tensor: grad_output * scalar
            d_tensor_local = grad_output_local * ctx.scalar
            d_tensor = DTensor.from_local(
                d_tensor_local,
                device_mesh=ctx.device_mesh_input,
                placements=ctx.placements_input,
                shape=grad_output.shape,
                stride=grad_output.stride(),
            )
        elif ctx.op == ElementwiseOp.DIV:
            (tensor_local,) = ctx.saved_tensors
            # Gradient w.r.t. tensor: -grad_output * scalar / tensor^2
            d_tensor_local = -grad_output_local * ctx.scalar / (tensor_local**2)
            d_tensor = DTensor.from_local(
                d_tensor_local,
                device_mesh=ctx.device_mesh_input,
                placements=ctx.placements_input,
                shape=grad_output.shape,
                stride=grad_output.stride(),
            )
        elif ctx.op == ElementwiseOp.GT:  # dummy op, not differentiable
            pass  # no gradient through this op
        elif ctx.op == ElementwiseOp.LT:  # dummy op, not differentiable
            pass  # no gradient through this op
        elif ctx.op == ElementwiseOp.EQUAL:  # dummy op, not differentiable
            pass  # no gradient through this op
        elif ctx.op == ElementwiseOp.POW:
            if ctx.scalar == 0:
                d_tensor_local = torch.zeros_like(grad_output_local)
            else:
                (tensor_local,) = ctx.saved_tensors
                d_tensor_local = grad_output_local * ctx.scalar * tensor_local ** (ctx.scalar - 1)
            d_tensor = DTensor.from_local(
                d_tensor_local,
                device_mesh=ctx.device_mesh_input,
                placements=ctx.placements_input,
                shape=grad_output.shape,
                stride=grad_output.stride(),
            )
        elif ctx.op == ElementwiseOp.MAX:
            # Gradient of max(scalar, tensor): passes through where tensor >= scalar, else 0
            (mask_local,) = ctx.saved_tensors
            d_tensor_local = grad_output_local * mask_local
            d_tensor = DTensor.from_local(
                d_tensor_local,
                device_mesh=ctx.device_mesh_input,
                placements=ctx.placements_input,
                shape=grad_output.shape,
                stride=grad_output.stride(),
            )
        else:
            raise ValueError(f"Unsupported scalar-tensor operation: {ctx.op}")

        return None, d_tensor, None


# TODO reduce human error by switching the order (scalar, tensor) to (tensor, scalar)
def scalar_tensor_op(scalar: float | int, tensor: DTensor, op: ElementwiseOp) -> DTensor:
    """Apply scalar-tensor operation to a scalar and distributed tensor.

    This function performs element-wise operations between a scalar and a distributed tensor.
    Supported operations are summation (scalar + tensor), subtraction (scalar - tensor), multiplication (scalar * tensor),
    division (scalar / tensor), greater than comparison (scalar > tensor), less than comparison (scalar < tensor), equality comparison (scalar == tensor), and power (tensor ** scalar). The operation is performed efficiently using local
    tensor operations while maintaining gradient computation capabilities.

    Parameters
    ----------
    scalar : float | int
        Scalar value to operate with.
    tensor : DTensor
        Input tensor. Can have any shape and sharding strategy.
    op : ElementwiseOp
        The operation to perform (ElementwiseOp.SUM, ElementwiseOp.SUB, ElementwiseOp.PROD, ElementwiseOp.DIV, ElementwiseOp.GT, ElementwiseOp.LT, ElementwiseOp.EQUAL, ElementwiseOp.POW, or ElementwiseOp.MAX).

    Returns
    -------
    DTensor
        Output tensor with shape identical to input tensor.
        Contains the result of the specified operation.

    Examples
    --------
    >>> # Assume we have distributed tensor x with shape (B, N, D)
    >>> sum_output = scalar_tensor_op(2.0, x, ElementwiseOp.SUM)
    >>> # sum_output = 2.0 + x, computed in distributed fashion
    >>>
    >>> sub_output = scalar_tensor_op(2.0, x, ElementwiseOp.SUB)
    >>> # sub_output = 2.0 - x, computed in distributed fashion
    >>>
    >>> prod_output = scalar_tensor_op(0.5, x, ElementwiseOp.PROD)
    >>> # prod_output = 0.5 * x, computed in distributed fashion
    >>>
    >>> div_output = scalar_tensor_op(1.0, x, ElementwiseOp.DIV)
    >>> # div_output = 1.0 / x, computed in distributed fashion
    >>>
    >>> gt_output = scalar_tensor_op(0.5, x, ElementwiseOp.GT)
    >>> # gt_output = 0.5 > x, computed in distributed fashion (boolean tensor)
    >>>
    >>> lt_output = scalar_tensor_op(0.5, x, ElementwiseOp.LT)
    >>> # lt_output = 0.5 < x, computed in distributed fashion (boolean tensor)
    >>>
    >>> equal_output = scalar_tensor_op(0.5, x, ElementwiseOp.EQUAL)
    >>> # equal_output = 0.5 == x, computed in distributed fashion (boolean tensor)
    >>>
    >>> pow_output = scalar_tensor_op(2.0, x, ElementwiseOp.POW)
    >>> # pow_output = x ** 2.0, computed in distributed fashion

    Notes
    -----
    - Input tensor must be a DTensor with any placement strategy
    - Partial placements are not currently supported
    - The function is differentiable and supports gradient computation for both scalar and tensor (except GT, LT, and EQUAL which are not differentiable)
    - The operation is performed on local tensor chunks for efficiency

    Raises
    ------
    TypeError
        If inputs are not of expected types.
    ValueError
        If Partial placements are used (not supported), or if op is invalid.
    """
    return _ScalarTensorOpImpl.apply(scalar, tensor, op)  # type: ignore
