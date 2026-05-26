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

from torch import Tensor
from torch.distributed.tensor import DTensor, Shard
from torch.distributed.tensor.placement_types import Replicate


def gather_along_cp(dtensor: DTensor) -> Tensor:
    """Gather a DTensor over CP dimensions, keeping the DP shard, then unwrap to a plain Tensor.

    Redistributes CP mesh dimensions (all except dim 0) to Replicate while
    preserving Shard(0) on the DP dimension. The returned tensor is the
    local DP slice with full spatial extent.

    Parameters
    ----------
    dtensor : DTensor
        Input distributed tensor with arbitrary placements.
        Expects the first mesh dimension to be DP.

    Returns
    -------
    Tensor
        The CP-gathered plain tensor, local to this DP rank.
    """
    expected_mesh_dim_names = ("dp", "cp_axis_0", "cp_axis_1")
    if dtensor.device_mesh.mesh_dim_names != expected_mesh_dim_names:
        raise ValueError(
            "gather_along_cp expects device mesh dim names "
            f"{expected_mesh_dim_names}, got {dtensor.device_mesh.mesh_dim_names}."
        )
    target_placements = [Shard(0)] + [Replicate()] * (dtensor.device_mesh.ndim - 1)
    gathered = dtensor.redistribute(dtensor.device_mesh, target_placements)
    return gathered.to_local()
