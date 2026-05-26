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


from torch import nn
from torch.distributed.device_mesh import DeviceMesh
from torch.distributed.tensor import DTensor

from boltz.distributed.model.layers.elementwise_op import ElementwiseOp, elementwise_op
from boltz.distributed.model.layers.layernorm import LayerNormParamsReplicated
from boltz.distributed.model.layers.linear import LinearParamsReplicated
from boltz.distributed.model.layers.sigmoid_gate import sigmoid_gate
from boltz.model.layers.transition import Transition as SerialTransition


class Transition(nn.Module):
    """Distributed two-layer MLP using DTensor."""

    def __init__(
        self,
        layer: SerialTransition,
        device_mesh: DeviceMesh,
    ) -> None:
        """Initialize the distributed Transition module.

        Parameters
        ----------
        layer : SerialTransition
            The serial transition layer containing weights to be distributed.
        device_mesh : DeviceMesh
            Device mesh defining the distributed computation topology.
        """
        super().__init__()
        self.device_mesh = device_mesh
        self.hidden = layer.hidden

        # Map serial layers to distributed versions
        self.norm = LayerNormParamsReplicated(layer.norm, self.device_mesh)
        self.fc1 = LinearParamsReplicated(layer.fc1, self.device_mesh)
        self.fc2 = LinearParamsReplicated(layer.fc2, self.device_mesh)
        self.fc3 = LinearParamsReplicated(layer.fc3, self.device_mesh)

    def forward(self, x: DTensor) -> DTensor:
        """Perform a forward pass.

        Parameters
        ----------
        x : DTensor
            The input data of shape (..., D)

        Returns
        -------
        DTensor
            The output data of shape (..., D)
        """
        x = self.norm(x)

        fc1_out = self.fc1(x)
        fc2_out = self.fc2(x)

        # SwiGLU activation: silu(fc1_out) * fc2_out
        # Since SiLU(x) = x * sigmoid(x), we use sigmoid_gate(fc1_out, fc1_out) * fc2_out
        # NOTE: self.fc1 and self.fc2 have the same dimensionality mapping: dim -> hidden
        # so fc1_out and fc2_out are of the same shape
        x = sigmoid_gate(fc1_out, fc1_out)
        x = elementwise_op(x, fc2_out, ElementwiseOp.PROD)

        x = self.fc3(x)
        return x
