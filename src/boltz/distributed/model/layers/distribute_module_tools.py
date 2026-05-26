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

from typing import Union

from einops.layers.torch import Rearrange
from torch import nn
from torch.distributed.device_mesh import DeviceMesh
from torch.nn import Module

from boltz.distributed.model.layers.layernorm import LayerNormParamsReplicated
from boltz.distributed.model.layers.linear import LinearParamsReplicated


def _convert_each_child_module_to_dtensor_compatible_version(
    module: Module,
    input_module: Module,
    device_mesh: DeviceMesh,
    reduction: str = "contain",
    module_name: Union[str, None] = None,
):
    """
    This function creates an attribute for module that is a DTensor
    compatible version of each of the child modules of input_module.

    For example, in the module.__init__(), this function can be called to
    create attributes of module.  The user may still have some additional
    steps to take in module.forward() to use the DTensor API.

    Parameters
    ----------
    module : Module
        The module that will be modified in place to have a DTensor API.
    input_module : Module
        The non-dtensor module with child modules that do not use DTensor.
    device_mesh : DeviceMesh
        The device mesh.
    reduction : str, optional
        - "contain": An attribute is created for module that is a DTensor
        compatible version of a non-dtensor child of input_module.
        - "sequential": An attribute is created for module that is a
    module_name : str, optional
        The name of the module to be replaced.
        This is only used when reduction is "sequential".
        If not provided, the module name is inferred from the input module.

    Returns
    -------
    None

    Raises
    ------
    NotImplementedError
        If the dtensor version of the child module is not implemented.
    """
    names_to_be_replaced: list[str] = []
    child_replacements: list[Module] = []
    for name, child in input_module.named_children():
        if isinstance(child, nn.Linear):
            names_to_be_replaced.append(name)
            child_replacements.append(LinearParamsReplicated(child, device_mesh=device_mesh))

        elif isinstance(child, nn.LayerNorm):
            names_to_be_replaced.append(name)
            child_replacements.append(LayerNormParamsReplicated(child, device_mesh=device_mesh))

        elif isinstance(child, nn.Sequential):
            _convert_each_child_module_to_dtensor_compatible_version(
                module, child, device_mesh=device_mesh, reduction="sequential", module_name=name
            )

        elif isinstance(child, Rearrange):
            # Not yet implemented
            pass

        else:
            raise NotImplementedError

    if reduction == "contain":
        for name, child_replacement in zip(names_to_be_replaced, child_replacements):
            setattr(module, name, child_replacement)
    elif reduction == "sequential" and module_name is not None:
        setattr(module, module_name, nn.Sequential(*child_replacements))
    else:
        raise NotImplementedError
