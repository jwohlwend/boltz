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

from torch.distributed.tensor import DTensor
from torch.nn import Module

from boltz.distributed.model.layers.cat_and_chunk import shardwise_chunk
from boltz.distributed.model.layers.elementwise_op import ElementwiseOp, elementwise_op
from boltz.distributed.model.layers.sigmoid_gate import sigmoid_gate


class SwiGLU(Module):
    """SwiGLU implemented as a module, for DTensor inputs.

    Gradient computation is implemented in the backward method for the
    implementations of shardwise_chunk, sigmoid_gate, and
    mult_for_same_placement_and_shape.

    See src/boltz/model/modules/utils.py for reference implementation.
    """

    def forward(self, x: DTensor) -> DTensor:
        """Forward pass of SwiGLU.

        DTensor metadata checking is performed in the implementation of
        each of shardwise_chunk, sigmoid_gate, and
        mult_for_same_placement_and_shape.

        Parameters
        ----------
        x : DTensor
            Input tensor.

        Returns
        -------
        DTensor
            Output tensor.

        Raises
        ------
        ValueError
            See _check_forward_input_for_impl_scope for more details.
        """
        self.check_forward_input_for_impl_scope(x)
        y, z = shardwise_chunk(x, chunks=2, dim=-1)
        a: DTensor = sigmoid_gate(x=z, g=z)
        a: DTensor = elementwise_op(a, y, ElementwiseOp.PROD)

        return a

    def check_forward_input_for_impl_scope(self, x: DTensor):
        """Check that the input tensor is compatible with the SwiGLU operation.

        The SwiGLU operations is defined only if the size of the last dimension
        is a multiple of 2.

        Parameters
        ----------
        x : DTensor
            Input tensor.

        Raises
        ------
        ValueError
            If the size of the last dimension is not a multiple of 2.
        """
        if hasattr(x, "shape") and x.shape[-1] % 2 != 0:
            raise ValueError(
                ", ".join(
                    [
                        "SwiGLU operation defined only if the size of the last dimension",
                        f"is a multiple of 2, whereas x.shape={x.shape}",
                    ]
                )
            )
