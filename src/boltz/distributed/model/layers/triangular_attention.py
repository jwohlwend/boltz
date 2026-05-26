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


import math
from enum import Enum, auto
from typing import Optional

import torch
from torch import nn
from torch.distributed.device_mesh import DeviceMesh
from torch.distributed.tensor import DTensor, Partial, Replicate, Shard, distribute_tensor

from boltz.distributed.comm import Ring2DCommTriAttn
from boltz.distributed.model.layers.layernorm import _LayerNormParamsReplicatedImpl
from boltz.distributed.model.layers.linear import LinearParamsReplicated, _LinearParamsReplicatedImpl
from boltz.distributed.model.layers.sigmoid_gate import sigmoid_gate
from boltz.distributed.model.modules.utils import TriAttnBackend
from boltz.distributed.utils import tiled_softmax_attention_update, update_exhaustive_strides
from boltz.model.layers.triangular_attention.attention import (
    TriangleAttentionEndingNode as SerialTriangleAttentionEndingNode,
)
from boltz.model.layers.triangular_attention.attention import (
    TriangleAttentionStartingNode as SerialTriangleAttentionStartingNode,
)
from boltz.model.layers.triangular_attention.primitives import (
    Attention,
)
from boltz.model.layers.triangular_attention.primitives import LayerNorm as SerialLayerNormNoAutoCastBF16
from boltz.model.layers.triangular_attention.primitives import Linear as SerialLinearNoAutoCastBF16
from boltz.model.layers.triangular_attention.utils import (
    permute_final_dims,
)

try:
    import cuequivariance_torch.primitives.triangle as cueq_triangle

    cueq_is_installed = True
except ImportError:
    cueq_is_installed = False

try:
    from trifast.torch import _triangle_attention as trifast_triangle_attention
    from trifast.torch import triangle_attention_bwd as trifast_triangle_attention_bwd

    trifast_is_installed = True
except ImportError:
    trifast_is_installed = False


def can_run_cueq_triattn_sm100f(
    device: torch.device,
    dtype: torch.dtype,
    dim_token: int,
    dim_hidden: int,
    is_fwd: bool,
) -> bool:
    """Check whether the cuEq SM100f triangle-attention kernel can run.

    Parameters
    ----------
    device : torch.device
        Target device (must be CUDA with SM100 or SM103 compute capability).
    dtype : torch.dtype
        Data type of q/k tensors.
    dim_token : int
        Token (sequence) dimension — ``q.shape[-2]`` or ``kT.shape[3]``.
    dim_hidden : int
        Per-head hidden dimension — ``q.shape[-1]``, i.e. ``c_hidden``.
    is_fwd : bool
        ``True`` for the forward pass, ``False`` for backward.

    Returns
    -------
    bool
        ``True`` when all SM100f constraints are satisfied.
    """
    if device.type != "cuda":
        return False
    if dtype not in (torch.bfloat16, torch.float16):
        return False
    device_cc = torch.cuda.get_device_capability(device)
    if device_cc not in ((10, 0), (10, 3)):
        return False
    if dim_hidden > 128 or dim_hidden % 8 != 0:
        return False
    if not is_fwd or dim_token % 8 == 0:
        return True
    return False


class LayerNormParamsReplicatedNoAutoCastBF16(nn.Module):
    """
    A LayerNorm module with replicated parameters for distributed training and disabled autocast in BF16.

    This module wraps around `_LayerNormParamsReplicatedImpl` to provide a user-friendly interface
    for LayerNorm operations using the DTensor API. It supports distributed training with replicated
    and sharded placements for input tensors and replicated placements for weight and bias tensors.

    Args:
        layer_local (nn.LayerNorm): An already-initialized nn.LayerNorm instance.
        device_mesh (DeviceMesh): The device mesh for distributed training.
    """

    def __init__(self, layer_local: SerialLayerNormNoAutoCastBF16, device_mesh: DeviceMesh) -> None:
        if not isinstance(layer_local, SerialLayerNormNoAutoCastBF16):
            raise ValueError(
                f"layer_local is not an instance of SerialLayerNormNoAutoCastBF16 but got {type(layer_local)}"
            )
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
        self.c_in = layer_local.c_in
        self.normalized_shape = list(self.c_in)
        self.eps = layer_local.eps
        self.device_mesh = device_mesh

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
        d = x.dtype
        if d is torch.bfloat16:
            with torch.autocast("cuda", enabled=False):
                out = _LayerNormParamsReplicatedImpl.apply(
                    x, self.normalized_shape, self.weight, self.bias, self.eps, True
                )
        else:
            out = _LayerNormParamsReplicatedImpl.apply(x, self.normalized_shape, self.weight, self.bias, self.eps)
        return out


class LinearParamsReplicatedNoAutoCastBF16(nn.Module):
    """
    Distributed linear layer with parameters replicated across all device mesh dimensions and disabled autocast in BF16.

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

    Key requirements:
        1. Parameters (weight and bias) will be replicated on all device mesh dimensions
        2. Input tensor and parameters must be on the same device mesh
        3. Feature/hidden dimension of the input must not be sharded across the device mesh
        4. Partial reduction along any input dimension is not supported
        5. Input and outputs must be on the same device mesh with the same placements
        6. Gradients of the input have the same placements on the same device mesh as the input
        7. Gradients of the weight and bias have Partial("sum") placements along the input's Shard placements'
           dimension so that the all-reduce will be performed along those device-grid dimensions

    """

    def __init__(self, layer_local: SerialLinearNoAutoCastBF16, device_mesh: DeviceMesh):
        """
        Initialize the distributed linear layer.

        Args:
            layer_local: nn.Linear to be distributed
            device_mesh: Device mesh for distributed computation
        """
        if not isinstance(layer_local, SerialLinearNoAutoCastBF16):
            raise ValueError(
                f"layer_local is not an instance of SerialLinearNoAutoCastBF16 but got {type(layer_local)}"
            )
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
        self.weight = nn.Parameter(distribute_tensor(layer_local.weight.data, device_mesh, all_replicate_placements))
        if layer_local.bias is None:
            self.register_parameter("bias", None)
        else:
            self.bias = nn.Parameter(distribute_tensor(layer_local.bias.data, device_mesh, all_replicate_placements))

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
        d = input.dtype
        if d is torch.bfloat16:
            with torch.autocast("cuda", enabled=False):
                return _LinearParamsReplicatedImpl.apply(input, self.weight, self.bias, True)
        else:
            return _LinearParamsReplicatedImpl.apply(input, self.weight, self.bias)


class _RingMultiHeadTriangleAttentionImpl(torch.autograd.Function):
    @staticmethod
    @torch.amp.custom_fwd(device_type="cuda")
    def forward(
        ctx,
        q_x: DTensor,
        kv_x: DTensor,
        mask: Optional[DTensor],
        triangle_bias: DTensor,
        weight_q: DTensor,
        weight_k: DTensor,
        weight_v: DTensor,
        no_heads: int,
        c_hidden: int,
        ring_comm: Ring2DCommTriAttn,
        inf: float,
        triattn_backend: TriAttnBackend,
    ) -> DTensor:
        """This function does the linear projection to prepare the q, k and v tensors
        and use them later to compute triangular attention.

        Linear projection and initial data shard redistribution:
           - Triangle bias is reorganized in two stages to avoid cross-rail traffic
           - Key/value pairs are initially shuffled to offset computation along
             attention matrix diagonal

        Stage 1 Triangle Bias Redistribution: Flatten diagonals onto rows/columns
        Original Data Ownership  After Stage 1 (for axis_cp=1)
        ┌───┬───┬───┐           ┌───┬───┬───┐
        │0,0│0,1│0,2│           │0,0│1,1│2,2│ original lower diagonal 0
        ├───┼───┼───┤           ├───┼───┼───┤
        │1,0│1,1│1,2│    →      │1,0│2,1│0,2│ original lower diagonal 1
        ├───┼───┼───┤           ├───┼───┼───┤
        │2,0│2,1│2,2│           │2,0│0,1│1,2│ original lower diagonal 2
        └───┴───┴───┘           └───┴───┴───┘

        Stage 2 Triangle Bias Redistribution: Rotate elements to meet ring attention requirements
        After Stage 2 (for axis_cp=1)
        ┌───┬───┬───┐
        │0,0│1,1│2,2│ original lower diagonal 0
        ├───┼───┼───┤
        │0,2│1,0│2,1│ original lower diagonal 1
        ├───┼───┼───┤
        │0,1│1,2│2,0│ original lower diagonal 2
        └───┴───┴───┘

        Forward pass for Multi-Head Triangle Attention using tri-axial
        virtual all_gather, reduce and all_gather. It uses a ring
        communication pattern across devices:
           - Each device maintains double buffers for k, v, triangle bias and mask
           - Data is shifted along the context parallelism axis in a ring pattern
           - Communication is overlapped with computation for efficiency

        Data sharding Diagram:
        ```
        The algorithm can be summarized as the 2-tuple (i, j) indexing the input data
        ownership of the triangle bias, where subsequent data ownership of other data
        is constrained by matching the corresponding i index of the triangle bias if
        the data contributes to the rows of the attention matrix (Q index) or the corresponding
        j index of the triangle bias if the data contributes to the columns of the
        attention matrix (K index).

        Initial Distribution - For axis_cp=1: (See Ring2DCommTriAttn for more explanation)

        initialized by _PrepQKVCommBiasImpl:
        ┌───┬───┬───┐
        │0,0│1,1│2,2│
        ├───┼───┼───┤
        │0,2│1,0│2,1│
        ├───┼───┼───┤
        │0,1│1,2│2,0│
        └───┴───┴───┘

        At each step, the 2-tuple are updated by upshifted 1 along each column
        e.g., at the end of step 0:
        ┌───┬───┬───┐
        │0,2│1,0│2,1│  # Column 0: 0,0→0,2→0,1 (wrapped)
        ├───┼───┼───┤
        │0,1│1,2│2,0│  # Column 1: 1,1→1,0→1,2
        ├───┼───┼───┤
        │0,0│1,1│2,2│  # Column 2: 2,2→2,1→2,0
        └───┴───┴───┘
        which gives the corresponding Q and K indices for each rank required for step 1

        At the end of step 2, the data is rotated to its "initialized" state,
        where the corresponding tensor along the Q and K axes are saved for
        backward pass

        Args:
            ctx: Context object for autograd
            q_x: Query input tensor
            kv_x: Key-value input tensor
            mask: Optional mask tensor (None creates all-ones mask)
            triangle_bias: Triangle bias tensor
            weight_q: Query projection weights
            weight_k: Key projection weights
            weight_v: Value projection weights
            no_heads: Number of attention heads
            c_hidden: Hidden dimension size
            ring_comm: Ring2DCommTriAttn object for distributed communication
            inf: Infinity value for mask bias computation
            triattn_backend: Triangular attention backend to use
        Returns:
            DTensor: Output tensor of shape [*, H, Q, C_hidden]
        """
        has_mask = mask is not None
        # Check if inputs are of type DTensor
        if not isinstance(q_x, DTensor):
            raise TypeError(f"Input 'q_x' must be of type DTensor. Got type {type(q_x)}.")
        if not isinstance(kv_x, DTensor):
            raise TypeError(f"Input 'kv_x' must be of type DTensor. Got type {type(kv_x)}.")
        if has_mask and not isinstance(mask, DTensor):
            raise TypeError(f"Input 'mask' must be of type DTensor or None. Got type {type(mask)}.")
        if not isinstance(triangle_bias, DTensor):
            raise TypeError(f"Input 'triangle_bias' must be of type DTensor. Got type {type(triangle_bias)}.")
        if not isinstance(weight_q, DTensor):
            raise TypeError(f"Input 'weight_q' must be of type DTensor. Got type {type(weight_q)}.")
        if not isinstance(weight_k, DTensor):
            raise TypeError(f"Input 'weight_k' must be of type DTensor. Got type {type(weight_k)}.")
        if not isinstance(weight_v, DTensor):
            raise TypeError(f"Input 'weight_v' must be of type DTensor. Got type {type(weight_v)}.")
        if not isinstance(triattn_backend, TriAttnBackend):
            raise TypeError(
                f"Input 'triattn_backend' must be of type TriAttnBackend. Got type {type(triattn_backend)}."
            )

        if triattn_backend not in [
            TriAttnBackend.CUEQ,
            TriAttnBackend.TRIFAST,
            TriAttnBackend.REFERENCE,
            TriAttnBackend.CUEQ_FWD_TRIFAST_BWD,
        ]:
            # to prevent accidental usage of unsupported backend
            # so that we don't have to handle the unsupported backend case in the later code
            # every time we have a backend selection logic
            raise NotImplementedError(
                f"Input 'triattn_backend' must be one of {TriAttnBackend.CUEQ, TriAttnBackend.TRIFAST, TriAttnBackend.REFERENCE, TriAttnBackend.CUEQ_FWD_TRIFAST_BWD}. "
                f"Got {triattn_backend}."
            )

        if triattn_backend in (TriAttnBackend.TRIFAST, TriAttnBackend.CUEQ_FWD_TRIFAST_BWD) and not has_mask:
            raise ValueError(
                "trifast or cueq_fwd_trifast_bwd backend requires a mask but mask is None. "
                "Please provide a all-zeros mask to indicate all-valid elements for trifast"
            )

        if inf > torch.finfo(q_x.dtype).max:
            raise ValueError(
                f"Input 'inf'={inf} is larger than max value of dtype {q_x.dtype}: {torch.finfo(q_x.dtype).max}"
            )

        # Check if inputs have identical device mesh
        device_mesh_input = q_x.device_mesh
        if device_mesh_input != kv_x.device_mesh:
            raise ValueError(
                f"Input tensors 'q_x' and 'kv_x' must have identical device mesh. "
                f"Got device meshes {device_mesh_input} and {kv_x.device_mesh}."
            )
        if has_mask and device_mesh_input != mask.device_mesh:
            raise ValueError(
                f"Input tensors 'q_x' and 'mask' must have identical device mesh. "
                f"Got device meshes {device_mesh_input} and {mask.device_mesh}."
            )
        if device_mesh_input != triangle_bias.device_mesh:
            raise ValueError(
                f"Input tensors 'q_x' and 'triangle_bias' must have identical device mesh. "
                f"Got device meshes {device_mesh_input} and {triangle_bias.device_mesh}."
            )
        if device_mesh_input != weight_q.device_mesh:
            raise ValueError(
                f"Input tensors 'q_x' and 'weight_q' must have identical device mesh. "
                f"Got device meshes {device_mesh_input} and {weight_q.device_mesh}."
            )
        if device_mesh_input != weight_k.device_mesh:
            raise ValueError(
                f"Input tensors 'q_x' and 'weight_k' must have identical device mesh. "
                f"Got device meshes {device_mesh_input} and {weight_k.device_mesh}."
            )
        if device_mesh_input != weight_v.device_mesh:
            raise ValueError(
                f"Input tensors 'q_x' and 'weight_v' must have identical device mesh. "
                f"Got device meshes {device_mesh_input} and {weight_v.device_mesh}."
            )

        # Check if q_x, kv_x, and mask_bias have the expected placements (Shard(0), Shard(1), Shard(2))

        expected_input_placements = (Shard(0), Shard(1), Shard(2))

        if q_x.placements != expected_input_placements:
            raise ValueError(
                f"Input tensor 'q_x' must have placements {expected_input_placements}. Got placements {q_x.placements}."
            )
        if kv_x.placements != expected_input_placements:
            raise ValueError(
                f"Input tensor 'kv_x' must have placements {expected_input_placements}. "
                f"Got placements {kv_x.placements}."
            )

        placements_input = q_x.placements

        # Create placement mapping for weight gradients
        # Weight gradients should have Partial("sum") placements corresponding to input's Shard placements
        placements_dweights = [Partial("sum"), Partial("sum"), Partial("sum")]

        # Check mask placements - should match q_x and kv_x
        if has_mask and mask.placements != expected_input_placements:
            raise ValueError(
                f"Input tensor 'mask' must have the same placements as 'q_x' and 'kv_x'. "
                f"Expected placements {expected_input_placements}, got {mask.placements}."
            )

        # Check triangle_bias placements - should be sharded along batch, I, and J dimensions
        # triangle_bias should have shape [B, I, J, H] and be sharded on dimensions 0, 1, and 2
        if triangle_bias.placements != placements_input:
            raise ValueError(
                f"Input tensor 'triangle_bias' must be sharded along batch, I, and J dimensions. "
                f"Expected placements {placements_input}, got {triangle_bias.placements}."
            )

        # Check weight placements - should be all replicated
        all_replicate_placements = [Replicate()] * device_mesh_input.ndim
        if weight_q.placements != tuple(all_replicate_placements):
            raise ValueError(
                f"Weight tensor 'weight_q' must have all replicated placements. "
                f"Expected {tuple(all_replicate_placements)}, got {weight_q.placements}."
            )
        if weight_k.placements != tuple(all_replicate_placements):
            raise ValueError(
                f"Weight tensor 'weight_k' must have all replicated placements. "
                f"Expected {tuple(all_replicate_placements)}, got {weight_k.placements}."
            )
        if weight_v.placements != tuple(all_replicate_placements):
            raise ValueError(
                f"Weight tensor 'weight_v' must have all replicated placements. "
                f"Expected {tuple(all_replicate_placements)}, got {weight_v.placements}."
            )

        # Check ring_comm consistency
        coord_device_mesh_input = device_mesh_input.get_coordinate()
        if coord_device_mesh_input is None:
            raise ValueError(
                f"ring_comm.coord_2d {ring_comm.coord_2d} is not on device_mesh_input {device_mesh_input}."
            )
        if ring_comm.coord_2d != (coord_device_mesh_input[1], coord_device_mesh_input[2]):
            raise ValueError(
                f"Input ring_comm's coord_2d {ring_comm.coord_2d} does not match the "
                f"device mesh's rank coordinates {coord_device_mesh_input} for the sharded dimensions."
            )

        if q_x.shape != kv_x.shape:
            raise ValueError(
                f"Input tensors 'q_x' and 'kv_x' must have the same shape. Got shapes {q_x.shape} and {kv_x.shape}."
            )

        if has_mask and mask.shape != q_x.shape[:-1]:
            raise ValueError(
                f"Input tensor 'mask' must have the same shape as 'q_x' and 'kv_x' except the last dimension. "
                f"Got shapes {mask.shape} and {q_x.shape[:-1]}."
            )

        if triangle_bias.shape != q_x.shape[:-1] + (no_heads,):
            raise ValueError(
                f"Input tensor 'triangle_bias' must have the same shape as 'q_x' and 'kv_x' "
                f"except the last dimension and the last dimension must be equal to no_heads. "
                f"Got shapes {triangle_bias.shape} and {q_x.shape[:-1] + (no_heads,)}."
            )

        # To accommodate the hybrid TriAttnBackend.CUEQ_FWD_TRIFAST_BWD, we dedicated
        # two working flags for the respective fwd and bwd cases in order to reuse the
        # cueq and trifast logics for the hybrid mode. But we need to modify the tensor shape
        # in the hybrid mode due to the different requirements of the two backends.
        # The shapes for the input between cueq and trifast are:
        # q: [B, I, H, Q, C_hidden] vs [B, H, I, Q, C_hidden]
        # kT: [B, I, H, K, C_hidden] vs [B, H, I, K, C_hidden]
        # v: [B, I, H, V, C_hidden] vs [B, H, I, V, C_hidden]
        # triangle_bias: [B, 1, H, I, J] vs [B, H, I, J]
        # mask: [B, I, 1, 1, J] vs [B, I, J]
        triattn_backend_fwd = (
            TriAttnBackend.CUEQ
            if triattn_backend in (TriAttnBackend.CUEQ, TriAttnBackend.CUEQ_FWD_TRIFAST_BWD)
            else triattn_backend
        )
        triattn_backend_bwd = (
            TriAttnBackend.TRIFAST
            if triattn_backend in (TriAttnBackend.TRIFAST, TriAttnBackend.CUEQ_FWD_TRIFAST_BWD)
            else triattn_backend
        )

        if triattn_backend == TriAttnBackend.REFERENCE:
            # we manage the scale in this function scope using torch op in fwd and bwd
            apply_scale = True
        elif triattn_backend in (TriAttnBackend.CUEQ, TriAttnBackend.TRIFAST, TriAttnBackend.CUEQ_FWD_TRIFAST_BWD):
            # we let the TriAttnBackend handle the scale internally in fwd and bwd
            apply_scale = False

        if has_mask:
            ctx.mark_non_differentiable(mask)
        ctx.device_mesh_input = device_mesh_input
        ctx.placements_input = placements_input
        ctx.placements_dweights = placements_dweights
        ctx.input_shape = q_x.shape
        ctx.input_stride = q_x.stride()
        ctx.weight_q_shape = weight_q.shape
        ctx.weight_q_stride = weight_q.stride()
        ctx.weight_k_shape = weight_k.shape
        ctx.weight_k_stride = weight_k.stride()
        ctx.weight_v_shape = weight_v.shape
        ctx.weight_v_stride = weight_v.stride()
        ctx.triangle_bias_shape = triangle_bias.shape
        ctx.triangle_bias_stride = triangle_bias.stride()
        ctx.ring_comm = ring_comm
        q_scale = 1.0 / math.sqrt(c_hidden)  # multiplicative factor for q
        ctx.q_scale = q_scale
        ctx.apply_scale = apply_scale
        ctx.no_heads = no_heads
        ctx.c_hidden = c_hidden
        ctx.has_mask = has_mask
        ctx.triattn_backend_bwd = triattn_backend_bwd

        # Store the mode based on ring_comm.axis_cp
        ctx.mode = _Mode.Ending if ring_comm.axis_cp == 0 else _Mode.Starting

        # Convert DTensors to local tensors for computation
        q_x_local = q_x.to_local()
        kv_x_local = kv_x.to_local()
        mask_local = mask.to_local() if has_mask else None
        triangle_bias_local = triangle_bias.to_local()

        # Handle transpose for ending mode (when ring_comm.axis_cp == 0)
        if ctx.mode == _Mode.Ending:
            # Transpose input tensors from [*, I, J, C] to [*, J, I, C]
            q_x_local = q_x_local.transpose(-2, -3).contiguous()
            kv_x_local = kv_x_local.transpose(-2, -3).contiguous()
            # Transpose triangle_bias from [*, I, J, H] to [*, J, I, H]
            triangle_bias_local = triangle_bias_local.transpose(-2, -3)
            # Transpose mask from [*, I, J] to [*, J, I] if mask exists
            if has_mask:
                mask_local = mask_local.transpose(-1, -2)

        if has_mask:
            if triattn_backend_fwd == TriAttnBackend.CUEQ:
                # if not casting to bool here, cueq will cast it internally anyway
                # so might as well do it here to save some communication bandwidth.
                # Also convert mask to mask_bias: [*, I, J] -> [*, I, 1, 1, J] for cueq
                # TODO: we should cast mask to bool from the dataloader
                mask_bias_local = mask_local[..., :, None, None, :].to(
                    dtype=torch.bool, memory_format=torch.contiguous_format
                )
            elif triattn_backend_fwd == TriAttnBackend.TRIFAST:
                # TRIFAST mask convention: True for invalid positions, False for valid positions
                # TRIFAST mask is of shape [*, I, J]
                mask_bias_local = ~(mask_local.to(dtype=torch.bool, memory_format=torch.contiguous_format))
            elif triattn_backend_fwd == TriAttnBackend.REFERENCE:
                # REFERENCE mask is of shape [*, I, 1, 1, J] and it's an additive mask bias
                mask_bias_local = inf * (mask_local - 1)
                mask_bias_local = mask_bias_local[..., :, None, None, :].contiguous()
        else:
            mask_bias_local = None

        if triattn_backend_fwd == TriAttnBackend.TRIFAST:
            # Convert triangle_bias from [*, I, J, H] to [*, H, I, J] for TRIFAST
            triangle_bias_local = permute_final_dims(triangle_bias_local, (2, 0, 1)).contiguous()
        else:
            # Convert triangle_bias from [*, I, J, H] to [*, 1, H, I, J] for CUEQ and REFERENCE
            triangle_bias_local = permute_final_dims(triangle_bias_local, (2, 0, 1)).unsqueeze(-4).contiguous()
        weight_q_local = weight_q.to_local()
        weight_k_local = weight_k.to_local()
        weight_v_local = weight_v.to_local()

        # send kv_x before computing k and v to avoid sending the latters
        # because in general no_heads > 1 so the linear projection expands
        # the size of kv_x by no_heads times
        kv_x_recv = ring_comm.comm_k_init.enqueue_to_dispatch(kv_x_local)

        # send mask along the axis_cp
        if has_mask:
            mask_bias_recv = ring_comm.comm_mask_init.enqueue_to_dispatch(mask_bias_local)
        else:
            mask_bias_recv = None

        # initialize triangle_bias comm for stage 1
        triangle_bias_recv = ring_comm.comm_bias_init0.enqueue_to_dispatch(triangle_bias_local)

        # [*, Q/K/V, H * C_hidden]
        q = torch.nn.functional.linear(q_x_local, weight_q_local)
        # [*, Q/K, H, C_hidden]
        q = q.view(q.shape[:-1] + (no_heads, -1))

        if triattn_backend_fwd == TriAttnBackend.TRIFAST:
            # [B, I, Q/K, H, C_hidden] --> [B, H, I, Q/K, C_hidden]
            q = permute_final_dims(q, (2, 0, 1, 3))
        else:
            # Both CUEQ and REFERENCE expect q to be of shape [*, H, Q/K, C_hidden]
            q = q.transpose(-2, -3)

        if apply_scale:
            q *= q_scale

        ring_comm.comm_k_init.wait_until_finished()
        # kv_x_recv is ready

        # compute q, k and v
        # kT == k.T is returned
        # [*, Q/K/V, H * C_hidden]
        k = torch.nn.functional.linear(kv_x_recv, weight_k_local)
        # [*, Q/K, H, C_hidden]
        k = k.view(k.shape[:-1] + (no_heads, -1))
        if triattn_backend_fwd == TriAttnBackend.CUEQ:
            # cueq expects k instead of its transpose
            # get kT (virtually k) of shape [*, H, K, C_hidden]
            kT = k.transpose(-2, -3).contiguous()
        elif triattn_backend_fwd == TriAttnBackend.TRIFAST:
            # [B, I, Q/K, H, C_hidden] --> [B, H, I, Q/K, C_hidden]
            kT = permute_final_dims(k, (2, 0, 1, 3)).contiguous()
        elif triattn_backend_fwd == TriAttnBackend.REFERENCE:
            # get k.T of shape [*, H, C_hidden, K]
            kT = permute_final_dims(k, (1, 2, 0))
            # torch.distributed data transfer requires contiguous tensor
            kT = kT.contiguous()

        # wait and initialize triangle_bias comm for stage 2
        ring_comm.comm_bias_init0.wait_until_finished()

        # Due to the two-stage communication, triangle_bias_recv needs an additional
        # buffer. To make the parent module RingMultiHeadTriangleAttention satisfy
        # the requirements of not modifying the input tensor, an additional copy
        # of triangle_bias is required
        triangle_bias_recv1 = ring_comm.comm_bias_init1.enqueue_to_dispatch(triangle_bias_recv)

        # [*, Q/K/V, H * C_hidden]
        v = torch.nn.functional.linear(kv_x_recv, weight_v_local)
        # [*, Q/K, H, C_hidden]
        v = v.view(v.shape[:-1] + (no_heads, -1))
        if triattn_backend_fwd == TriAttnBackend.TRIFAST:
            # [B, I, Q/K, H, C_hidden] --> [B, H, I, Q/K, C_hidden]
            v = permute_final_dims(v, (2, 0, 1, 3))
        else:
            # Both CUEQ and REFERENCE expect v to be of shape [*, H, Q/K, C_hidden]
            v = v.transpose(-2, -3)
        # torch.distributed data transfer requires contiguous tensor
        v = v.contiguous()

        # initial triangle_bias should be ready by now
        # TODO: move this wait inside the loop right before it's needed
        if has_mask:
            ring_comm.comm_mask_init.wait_until_finished()
        ring_comm.comm_bias_init1.wait_until_finished()

        # triangle_bias_recv1 is ready
        # mask_bias_recv is ready

        i_ready = 0
        i_recv = i_ready ^ 1
        kT_buffer = [kT, torch.empty_like(kT)]
        v_buffer = [v, torch.empty_like(v)]
        triangle_bias_buffer = [triangle_bias_recv1, torch.empty_like(triangle_bias_recv)]
        if has_mask:
            mask_bias_buffer = [mask_bias_recv, torch.empty_like(mask_bias_recv)]
        else:
            mask_bias_buffer = None
        o, lse_m, amax = None, None, None
        n_steps = ring_comm.group_layout.shape[ring_comm.axis_cp]
        for step in range(n_steps):
            # launch send/recv for the next round
            # This is done even for the last step to enable saving the tensors for the backward pass
            kT_buffer[i_recv] = ring_comm.comm_k.enqueue_to_dispatch(kT_buffer[i_ready], kT_buffer[i_recv])
            v_buffer[i_recv] = ring_comm.comm_v.enqueue_to_dispatch(v_buffer[i_ready], v_buffer[i_recv])
            if has_mask:
                mask_bias_buffer[i_recv] = ring_comm.comm_mask.enqueue_to_dispatch(
                    mask_bias_buffer[i_ready], mask_bias_buffer[i_recv]
                )
            triangle_bias_buffer[i_recv] = ring_comm.comm_bias.enqueue_to_dispatch(
                triangle_bias_buffer[i_ready], triangle_bias_buffer[i_recv]
            )
            # proceed with current k, v and triangle_bias
            # NOTE: B is batch size; H is head; I and J are pair repr N_token
            # C_hidden is q/k/v embedding dim; Q/K/V are attention dim (N_token)
            # kT.shape == [*, H, C_hidden, K] (default torch variant) or [*, H, K, C_hidden] (cueq variant)
            # q.shape == [*, H, Q, C_hidden]

            if triattn_backend_fwd == TriAttnBackend.CUEQ:
                o_block, lse_block, amax_block = cueq_triangle.triangle_attention(
                    q,
                    kT_buffer[i_ready],
                    v_buffer[i_ready],
                    triangle_bias_buffer[i_ready],
                    mask_bias_buffer[i_ready] if has_mask else None,
                    scale=1.0 if apply_scale else q_scale,
                    return_aux=True,
                )
                amax_block = amax_block.unsqueeze(-1)
                lse_m_block = lse_block.unsqueeze(-1) - amax_block
                # cueq TriAttn returns lse and amax in FP32 so we need to cast them back
                # but the lse and amax returned can contain 1e9 values, which can overflow
                # fp16 so we need to clamp them to the max value of fp16
                if q.dtype == torch.float16:
                    amax_block = amax_block.clamp(
                        min=torch.finfo(torch.float16).min, max=torch.finfo(torch.float16).max
                    )
                    lse_m_block = lse_m_block.clamp(
                        min=torch.finfo(torch.float16).min, max=torch.finfo(torch.float16).max
                    )
                # TODO: verify if there is actually benefit in using FP32 lse_m and amax
                # in tiled_softmax_attention_update
                amax_block = amax_block.to(dtype=q.dtype)
                lse_m_block = lse_m_block.to(dtype=q.dtype)
            elif triattn_backend_fwd == TriAttnBackend.TRIFAST:
                # has_mask == False would have raised before reaching here with TRIFAST backend
                # o_block is of shape [B, H, I, K, C_hidden]
                # lse_block is of shape [B, H, I, Q]
                o_block, lse_block = trifast_triangle_attention(
                    q, kT_buffer[i_ready], v_buffer[i_ready], triangle_bias_buffer[i_ready], mask_bias_buffer[i_ready]
                )
                amax_block = None
                # TRIFAST returns lse directly (not lse - amax) in FP32. This is known to cause accuracy issues
                # due to lower dynamic range in tiled softmax update.
                # Pad a singleton K axis to lse_block to be used inside tiled_softmax_attention_update
                # Here we don't need to canonicalize the shape of o and lse to the REFERENCE backend's shape
                # because the tiled_softmax_attention_update effectively treats the leading axes as virtual
                # batch axes.
                lse_m_block = lse_block.to(dtype=q.dtype).unsqueeze(-1)
            elif triattn_backend_fwd == TriAttnBackend.REFERENCE:
                # [B, I, H, Q, K]
                a = torch.matmul(q, kT_buffer[i_ready])

                # biases[0].shape is [B, I, 1, 1, J]
                if has_mask:
                    a += mask_bias_buffer[i_ready]

                # triangle_bias.shape is [B, 1, H, I, J]
                a += triangle_bias_buffer[i_ready]

                # The following tries to stabilize pure -1e9 chunk
                # in a, which could happen towards the last few
                # chunks of the padding. This is done by keeping
                # track of the max of a and the lse - max(a) during
                # the accumulation. The tiled_softmax_attention_update will
                # first attempt to compute amax_block - amax, which
                # tends to cancel each other out, before updating
                # the accumulators
                amax_block = a.amax(dim=-1, keepdim=True)
                # [*, H, Q, 1]
                lse_m_block = torch.logsumexp(a - amax_block, dim=-1, keepdim=True).to(dtype=a.dtype)
                # [*, H, Q, K]
                a = torch.softmax(a, dim=-1)

                # [*, H, Q, C_hidden]
                o_block = torch.matmul(a, v_buffer[i_ready])

            o, lse_m, amax = tiled_softmax_attention_update(o_block, lse_m_block, amax_block, o, lse_m, amax)
            # wait until next block is ready
            ring_comm.comm_k.wait_until_finished()
            ring_comm.comm_v.wait_until_finished()
            ring_comm.comm_bias.wait_until_finished()
            if has_mask:
                ring_comm.comm_mask.wait_until_finished()
            i_ready ^= 1
            i_recv ^= 1
        # NOTE: The last step's communication is done to reset the data's ownership to its initial state
        # at the beginning of the forward pass
        # NOTE: Although backward pass doesn't need to do block-wise renormalization with amax but only
        # with lse, we need to subtract the masked attention matrix first by amax then lse_m to avoid
        # numerical instability so we need to save both lse_m and amax terms for backward pass
        if (
            q_x.requires_grad
            or kv_x.requires_grad
            or triangle_bias.requires_grad
            or weight_q.requires_grad
            or weight_k.requires_grad
            or weight_v.requires_grad
        ):
            # This should be enough to avoid saving tensors when no gradient is needed
            # TODO: tailor to the individual gradient in terms which intermediate tensors are saved
            # but the challenge is how to avoid deadlocking in the backward pass in branching
            if triattn_backend_fwd == triattn_backend_bwd:
                # when fwd and bwd are the same backend, the processed tensors should satisfy
                # the same shape requirements in the fwd and bwd pass
                ctx.save_for_backward(
                    q_x_local,
                    kv_x_recv,
                    weight_q_local,
                    weight_k_local,
                    weight_v_local,
                    q,
                    kT_buffer[i_ready],
                    v_buffer[i_ready],
                    triangle_bias_buffer[i_ready],
                    mask_bias_buffer[i_ready] if ctx.has_mask else None,
                    o,
                    amax,
                    lse_m,
                )
            elif triattn_backend_fwd == TriAttnBackend.CUEQ and triattn_backend_bwd == TriAttnBackend.TRIFAST:
                # this implies triattn_backend == TriAttnBackend.CUEQ_FWD_TRIFAST_BWD and has_mask is True
                # need to reshape some tensors as if the fwd pass was done with TRIFAST backend
                ctx.save_for_backward(
                    q_x_local,
                    kv_x_recv,
                    weight_q_local,
                    weight_k_local,
                    weight_v_local,
                    q.transpose(-3, -4).contiguous(),  # [B, I, H, Q, C_hidden] --> [B, H, I, Q, C_hidden]
                    kT_buffer[i_ready]
                    .transpose(-3, -4)
                    .contiguous(),  # [B, I, H, K, C_hidden] --> [B, H, I, K, C_hidden]
                    v_buffer[i_ready]
                    .transpose(-3, -4)
                    .contiguous(),  # [B, I, H, V, C_hidden] --> [B, H, I, V, C_hidden]
                    triangle_bias_buffer[i_ready].squeeze(-4).contiguous(),  # [B, 1, H, I, J] --> [B, H, I, J]
                    ~(  # trifast uses inverse boolean mask convention
                        mask_bias_buffer[i_ready].squeeze((-2, -3)).contiguous()
                    ),  # [B, I, 1, 1, J] --> [B, I, J]
                    o.transpose(-3, -4).contiguous(),  # [B, I, H, V, C_hidden] --> [B, H, I, V, C_hidden]
                    None,  # trifast doesn't need amax
                    (lse_m + amax)
                    # trifast assume lse instead of lse_m
                    .transpose(-3, -4)  # [B, I, H, Q, 1] --> [B, H, I, Q, 1]
                    .contiguous(),
                )
            else:
                raise NotImplementedError(
                    f"Unsupported backend {triattn_backend} with fwd backend {triattn_backend_fwd} and bwd backend {triattn_backend_bwd}."
                )

        if triattn_backend_fwd == TriAttnBackend.TRIFAST:
            # o is of shape [B, H, I, Q, C_hidden] or [B, H, J, Q, C_hidden] (if ending mode)
            # which needs to be transposed into [B, I, J, H * C_hidden] or [B, J, I, H * C_hidden] for
            # consistency with input tensor axis semantics and placements
            # and for downstream linear projection
            o = permute_final_dims(o, (1, 2, 0, 3)).flatten(start_dim=-2)
        else:
            # o is of shape [B, I, H, Q, C_hidden] or [B, J, H, Q, C_hidden] (if ending mode)
            # which needs to be transposed into [B, I, J, H * C_hidden] or [B, J, I, H * C_hidden] for
            # consistency with input tensor axis semantics and placements
            # and for downstream linear projection
            o = o.transpose(-2, -3).flatten(start_dim=-2)

        # Handle transpose back for ending mode
        if ctx.mode == _Mode.Ending:
            # Transpose output from [B, J, I, H * C_hidden] to [B, I, J, H * C_hidden]
            o = o.transpose(-2, -3)

        # Convert result back to DTensor
        shape_output = ctx.input_shape[:-1] + (o.shape[-1],)
        stride_output = update_exhaustive_strides(o.shape, o.stride(), shape_output)
        output = DTensor.from_local(
            o,
            device_mesh=device_mesh_input,
            placements=placements_input,
            shape=shape_output,
            stride=stride_output,
        )
        return output

    @staticmethod
    @torch.amp.custom_bwd(device_type="cuda")
    def backward(
        ctx, do: DTensor
    ) -> tuple[DTensor, DTensor, None, DTensor, DTensor, DTensor, DTensor, None, None, None, None, None, None]:
        """Backward pass for Multi-Head Triangle Attention using tri-axial
        virtual all_gather, reduce and all_reduce

        This implements the backward pass for the ring communication pattern used in the forward pass.
        The data ownership follows the same pattern as the forward pass, with gradients being accumulated
        and communicated in a ring pattern.

        Data Ownership Diagram:
        ```
        The algorithm can be summarized as the 2-tuple (i, j) indexing the input data
        ownership of the triangle bias, where subsequent data ownership of other data
        is constrained by matching the corresponding i index of the triangle bias if
        the data contributes to the rows of the attention matrix (Q index) or the corresponding
        j index of the triangle bias if the data contributes to the columns of the
        attention matrix (K index).

        Initial Distribution (saved from forward) - For axis_cp=1: (See Ring2DCommTriAttn for more explanation)

        initialized:             end of step 0
        ┌───┬───┬───┐            ┌───┬───┬───┐
        │0,0│1,1│2,2│            │0,2│1,0│2,1│
        ├───┼───┼───┤  upshift   ├───┼───┼───┤  upshift
        │0,2│1,0│2,1│  ---->     │0,1│1,2│2,0│  ---->  ...
        ├───┼───┼───┤            ├───┼───┼───┤
        │0,1│1,2│2,0│            │0,0│1,1│2,2│
        └───┴───┴───┘            └───┴───┴───┘

        This indexing scheme is exactly the same as the forward pass to all the associated
        intermediate tensors as well as the dependent gradients, whose Q and/or K indices
        strictly matching the 2-tuple shown above. The only data tensor that is distributed
        as the device grid index is the forward's output's gradient,
        i.e., device (i, j) owns do(i, j)

        ```

        Args:
            ctx: Context object containing saved tensors and ring communication object
            do: Gradient of output tensor of shape [B, I, J, H * C_hidden]

        Returns:
            Tuple of gradients for:
            - dq_x: Gradient of query input
            - dkv_x_recv: Gradient of key-value input
            - None: Placeholder for mask gradient (non-differentiable)
            - dtriangle_bias: Gradient of triangle bias
            - dweight_q: Gradient of query weight
            - dweight_k: Gradient of key weight
            - dweight_v: Gradient of value weight
            - None: Placeholder for unused gradients
        """
        # Check if input is of type DTensor
        if not isinstance(do, DTensor):
            raise TypeError(f"Input 'do' must be of type DTensor. Got type {type(do)}.")

        # Check if input has same device mesh and placements as forward inputs
        if do.device_mesh != ctx.device_mesh_input:
            raise ValueError(
                f"Input 'do' must have the same device mesh as the input tensors. "
                f"Got device meshes {do.device_mesh} and {ctx.device_mesh_input}."
            )
        if do.placements != ctx.placements_input:
            raise ValueError(
                f"Input 'do' must have the same placements as the input tensors. "
                f"Got placements {do.placements} and {ctx.placements_input}."
            )

        # Convert gradient input to local tensor
        do_local = do.to_local()

        # Handle transpose for ending mode gradient input
        if ctx.mode == _Mode.Ending:
            # Transpose gradient from [B, I, J, H * C_hidden] to [B, J, I, H * C_hidden]
            do_local = do_local.transpose(-2, -3)
            # else do_local is of shape [B, I, J, H * C_hidden]

        if ctx.triattn_backend_bwd == TriAttnBackend.TRIFAST:
            # flatten and move the head dimension up to: [B, H, I/J, Q, C_hidden]
            do_local = permute_final_dims(
                do_local.unflatten(-1, (ctx.no_heads, ctx.c_hidden)), (2, 0, 1, 3)
            ).contiguous()
        else:
            # flatten and move the head dimension up to: [B, I/J, H, Q, C_hidden]
            do_local = do_local.unflatten(-1, (ctx.no_heads, ctx.c_hidden)).transpose(-2, -3).contiguous()

        (
            q_x,
            kv_x_recv,
            weight_q,
            weight_k,
            weight_v,
            q,
            kT_ready,
            v_ready,
            triangle_bias_ready,
            mask_bias_ready,
            o,
            amax,
            lse_m,
        ) = ctx.saved_tensors

        if ctx.triattn_backend_bwd == TriAttnBackend.CUEQ:
            can_run_sm100f = can_run_cueq_triattn_sm100f(q.device, q.dtype, kT_ready.shape[3], q.shape[-1], False)
            # cueq uses lse instead of lse_m + amax
            # it also requires the singleton K axis to be removed
            lse = (lse_m + amax).squeeze(-1)
        elif ctx.triattn_backend_bwd == TriAttnBackend.TRIFAST:
            # trifast lse_m is actually lse, which must have shape [B, H, I, Q]
            lse = lse_m.squeeze(-1)

        ring_comm: Ring2DCommTriAttn = ctx.ring_comm
        # o is saved from the forward pass and
        # is of shape [B, I, H, Q, C_hidden]

        if ctx.triattn_backend_bwd == TriAttnBackend.REFERENCE:
            # Only needed by REFERENCE backend: doT is of shape [B, I, H, C_hidden, Q]
            doT = do_local.transpose(-2, -1)
            # Only needed by REFERENCE backend: qT is of shape [B, I, H, Q, C_hidden] for dkT computation
            qT = q.transpose(-2, -1)
        dq = torch.empty_like(q)
        dkT = torch.empty_like(kT_ready)
        if ctx.triattn_backend_bwd in (TriAttnBackend.CUEQ, TriAttnBackend.TRIFAST):
            # this is virtually dv as will be returned from cueq triangle attention
            # For CUEQ: dvT is of shape [B, I, H, K, C_hidden]
            # For TRIFAST: dvT is of shape [B, H, I, K, C_hidden]
            dvT = torch.empty_like(v_ready, memory_format=torch.contiguous_format)
        elif ctx.triattn_backend_bwd == TriAttnBackend.REFERENCE:
            # Instead of transposing the attention matrix, we transpose
            # "do" to compute dvT instead
            # dvT is of shape [B, I, H, C_hidden, K]
            dvT = torch.empty(
                v_ready.shape[:-2] + (v_ready.shape[-1], v_ready.shape[-2]), dtype=v_ready.dtype, device=v_ready.device
            )
        dtriangle_bias = torch.empty_like(triangle_bias_ready)

        if ctx.triattn_backend_bwd in (TriAttnBackend.CUEQ, TriAttnBackend.TRIFAST):
            # d is computed internally in cueq/trifast triangle attention
            d = None
        elif ctx.triattn_backend_bwd == TriAttnBackend.REFERENCE:
            # d.shape is [B, I, H, Q, 1]
            d = torch.linalg.vecdot(do_local, o, dim=-1).unsqueeze(-1)
            # prevent d from promoting to fp32 if do_local is fp32
            d = d.to(dtype=o.dtype)

        i_ready = 0
        i_recv = i_ready ^ 1
        kT_buffer = [kT_ready, torch.empty_like(kT_ready)]
        dkT_buffer = [dkT, torch.empty_like(dkT)]
        v_buffer = [v_ready, torch.empty_like(v_ready)]
        dvT_buffer = [dvT, torch.empty_like(dvT)]
        triangle_bias_buffer = [triangle_bias_ready, torch.empty_like(triangle_bias_ready)]
        dtriangle_bias_buffer = [dtriangle_bias, torch.empty_like(dtriangle_bias)]
        if ctx.has_mask:
            mask_bias_buffer = [mask_bias_ready, torch.empty_like(mask_bias_ready)]
        else:
            mask_bias_buffer = None
        apply_scale = ctx.apply_scale
        q_scale = ctx.q_scale
        n_steps = ring_comm.group_layout.shape[ring_comm.axis_cp]
        for step in range(n_steps):
            is_last_step = step == n_steps - 1
            if not is_last_step:
                # launch send/recv for the next round
                # This is done even for the last step to enable saving the tensors for the backward pass
                kT_buffer[i_recv] = ring_comm.comm_k.enqueue_to_dispatch(kT_buffer[i_ready], kT_buffer[i_recv])
                v_buffer[i_recv] = ring_comm.comm_v.enqueue_to_dispatch(v_buffer[i_ready], v_buffer[i_recv])
                if ctx.has_mask:
                    mask_bias_buffer[i_recv] = ring_comm.comm_mask.enqueue_to_dispatch(
                        mask_bias_buffer[i_ready], mask_bias_buffer[i_recv]
                    )
                triangle_bias_buffer[i_recv] = ring_comm.comm_bias.enqueue_to_dispatch(
                    triangle_bias_buffer[i_ready], triangle_bias_buffer[i_recv]
                )

            # proceed with current k, v and triangle_bias
            # NOTE: B is batch size; H is head; I and J are pair repr N_token
            # C_hidden is q/k/v embedding dim; Q/K/V are attention dim (N_token)
            # kT.shape == [*, H, C_hidden, K] (default torch variant) or [*, H, K, C_hidden] (cueq variant)
            #          or [*, H, I, Q/K, C_hidden] (trifast variant)
            # q.shape == [*, H, Q, C_hidden] (default torch or cueq variant) or [*, H, I, Q/K, C_hidden] (trifast variant)

            if ctx.triattn_backend_bwd == TriAttnBackend.CUEQ:
                # dkT_block.shape is [B, I, H, K, C_hidden]
                # dvT_block.shape is [B, I, H, K, C_hidden]
                if can_run_sm100f:
                    # SM100f kernel accepts bias in the same dtype as q;
                    # keep the buffer's native dtype (no-op cast) and let
                    # cuEq's _convert_bias handle any necessary conversion.
                    bias_dtype = triangle_bias_buffer[i_ready].dtype
                else:
                    # Non-SM100f cuEq backward requires float32 bias and lse
                    bias_dtype = torch.float32
                dq_block, dkT_block, dvT_block, dtriangle_bias_block_fp32 = (
                    torch.ops.cuequivariance.triangle_attention_bwd(
                        do_local,
                        o,
                        q,
                        kT_buffer[i_ready],
                        v_buffer[i_ready],
                        triangle_bias_buffer[i_ready].to(dtype=bias_dtype),
                        mask_bias_buffer[i_ready] if ctx.has_mask else None,
                        lse.to(dtype=torch.float32),
                        1.0 if apply_scale else q_scale,
                    )
                )
                dtriangle_bias_block = dtriangle_bias_block_fp32.to(dtype=triangle_bias_buffer[i_ready].dtype)
            elif ctx.triattn_backend_bwd == TriAttnBackend.TRIFAST:
                # A fake dmask tensor is also returned by trifast_triangle_attention_bwd as the last return value
                dq_block, dkT_block, dvT_block, dtriangle_bias_block, _ = trifast_triangle_attention_bwd(
                    do_local,
                    q,
                    kT_buffer[i_ready],
                    v_buffer[i_ready],
                    triangle_bias_buffer[i_ready],
                    o,
                    lse.to(dtype=torch.float32),
                    mask_bias_buffer[i_ready],
                )
            elif ctx.triattn_backend_bwd == TriAttnBackend.REFERENCE:
                # [B, I, H, Q, K]
                a = torch.matmul(q, kT_buffer[i_ready])

                # biases[0].shape is [B, I, 1, 1, J]
                if ctx.has_mask:
                    a += mask_bias_buffer[i_ready]

                # triangle_bias.shape is [B, 1, H, I, J]
                a += triangle_bias_buffer[i_ready]

                # amax and lse_m shape is [B, I, H, Q, 1]
                a -= amax
                # lse_m is fp32 from logsumexp in fwd autocast, so we need to cast it to match a's dtype
                # to avoid promoting a to fp32
                a -= lse_m.to(dtype=a.dtype)

                a = torch.exp(a)

                # dvT_block.shape is [B, I, H, C_hidden, K]
                dvT_block = torch.matmul(doT, a)

                # da.shape is [B, I, H, Q, K]
                da = torch.matmul(do_local, v_buffer[i_ready].transpose(-1, -2))
                # a is no longer needed so we can repurpose its memory for ds
                # ds.shape is [B, I, H, Q, K]
                ds = a
                ds *= da - d

                # TODO: check if the cublas/cutlass backend is optimal with the
                # non-ideal memory layout of kT_buffer[i_ready].transpose(-2, -1)
                dq_block = torch.matmul(ds, kT_buffer[i_ready].transpose(-2, -1))

                # dkT_block.shape is [B, I, H, C_hidden, K]
                dkT_block = torch.matmul(qT, ds)

                # dtriangle_bias_block.shape is [B, 1, H, Q, K]
                dtriangle_bias_block = ds.sum(dim=-4, keepdim=True, dtype=triangle_bias_buffer[i_ready].dtype)

            if step == 0:
                dvT_buffer[i_ready] = dvT_block
                dq = dq_block
                dkT_buffer[i_ready] = dkT_block
                dtriangle_bias_buffer[i_ready] = dtriangle_bias_block
            else:
                dq += dq_block
                ring_comm.comm_dv.wait_until_finished()
                dvT_buffer[i_ready] += dvT_block
                ring_comm.comm_dk.wait_until_finished()
                dkT_buffer[i_ready] += dkT_block
                ring_comm.comm_dbias.wait_until_finished()
                dtriangle_bias_buffer[i_ready] += dtriangle_bias_block

            dvT_buffer[i_recv] = ring_comm.comm_dv.enqueue_to_dispatch(dvT_buffer[i_ready], dvT_buffer[i_recv])
            dkT_buffer[i_recv] = ring_comm.comm_dk.enqueue_to_dispatch(dkT_buffer[i_ready], dkT_buffer[i_recv])
            dtriangle_bias_buffer[i_recv] = ring_comm.comm_dbias.enqueue_to_dispatch(
                dtriangle_bias_buffer[i_ready], dtriangle_bias_buffer[i_recv]
            )

            if not is_last_step:
                # wait until next block is ready
                ring_comm.comm_k.wait_until_finished()
                ring_comm.comm_v.wait_until_finished()
                ring_comm.comm_bias.wait_until_finished()
                if ctx.has_mask:
                    ring_comm.comm_mask.wait_until_finished()
                i_ready ^= 1
                i_recv ^= 1

        # dv, dkT and dtriangle_bias need the extra round of send/recv so that the
        # data's ownership is transferred to the initial state at the beginning of the forward pass
        i_ready ^= 1
        i_recv ^= 1
        ring_comm.comm_dv.wait_until_finished()
        ring_comm.comm_dk.wait_until_finished()
        ring_comm.comm_dbias.wait_until_finished()

        dkT = dkT_buffer[i_ready]
        if ctx.triattn_backend_bwd == TriAttnBackend.CUEQ:
            # dv is already of shape [B, I, H, K, C_hidden]
            dv = dvT_buffer[i_ready]
        elif ctx.triattn_backend_bwd == TriAttnBackend.TRIFAST:
            # [B, H, I, K, C_hidden] --> [B, I, H, K, C_hidden]
            dv = dvT_buffer[i_ready].transpose(-4, -3)
        elif ctx.triattn_backend_bwd == TriAttnBackend.REFERENCE:
            # dvT_buffer[i_ready] is of shape [B, I, H, C_hidden, K]
            # reshaped to dv of shape [B, I, H, K, C_hidden]
            dv = dvT_buffer[i_ready].transpose(-2, -1).contiguous()
        dtriangle_bias = dtriangle_bias_buffer[i_ready]

        # the input tensors are sharded according to what's returned from _RingMHTAFunctorImpl.backward
        # Here, q_x and dq didn't go through shuffling and its data ownership remain stationary
        # kv_x_recv, dkT and dv are shuffled according to ring_comm's comm_k_init
        # while dtriangle_bias is shuffled according to ring_comm's comm_bias_init0 and comm_bias_init1
        # The strategy here is to perform local computation first to get dkv_x, which is then shuffled,
        # because dkv_x in general is smaller in size compared to dkvT and dvT due to no_heads > 1

        # q_x and kv_x_recv are of shape [B, I, J, C_hidden]

        dtriangle_bias_recv = ring_comm.comm_dbias_final0.enqueue_to_dispatch(dtriangle_bias)
        if ctx.triattn_backend_bwd in (TriAttnBackend.CUEQ, TriAttnBackend.REFERENCE):
            # [B, I, Q, H, C_hidden]
            dq_reshaped = dq.transpose(-2, -3)
        elif ctx.triattn_backend_bwd == TriAttnBackend.TRIFAST:
            # [B, H, I, Q, C_hidden] --> [B, I, Q, H, C_hidden]
            dq_reshaped = permute_final_dims(dq, (1, 2, 0, 3))
        # [B, I, Q, H * C_hidden]
        dq_reshaped = dq_reshaped.flatten(start_dim=-2)
        if apply_scale:
            dq_reshaped = dq_reshaped * q_scale
        dq_x = torch.einsum("...z, zc -> ...c", dq_reshaped, weight_q)
        ring_comm.comm_dbias_final0.wait_until_finished()
        dtriangle_bias = ring_comm.comm_dbias_final1.enqueue_to_dispatch(dtriangle_bias_recv, dtriangle_bias)

        # dweight_q is of shape [*, H * C_hidden, C_hidden]
        dweight_q = torch.einsum("...z, ...c -> zc", dq_reshaped, q_x)

        dweight_q_dtensor = DTensor.from_local(
            dweight_q,
            device_mesh=ctx.device_mesh_input,
            placements=ctx.placements_dweights,
            shape=ctx.weight_q_shape,
            stride=ctx.weight_q_stride,
        )

        # reduce the dimensionality of H * C_hidden to C_hidden
        # and sum up contribution from dv and dk before send/recv dkv_x
        # [B, I, K, H, C_hidden]
        dv_reshaped = dv.transpose(-2, -3)
        # [B, I, K, H * C_hidden]
        dv_reshaped = dv_reshaped.flatten(start_dim=-2)

        if ctx.triattn_backend_bwd == TriAttnBackend.CUEQ:
            # dkT is of shape [B, I, H, K, C_hidden]
            # dk_reshaped is of shape [B, I, K, H, C_hidden]
            dk_reshaped = dkT.transpose(-2, -3)
        elif ctx.triattn_backend_bwd == TriAttnBackend.TRIFAST:
            # [B, H, I, K, C_hidden] --> [B, I, K, H, C_hidden]
            dk_reshaped = permute_final_dims(dkT, (1, 2, 0, 3))
        elif ctx.triattn_backend_bwd == TriAttnBackend.REFERENCE:
            # dkT is of shape [B, I, H, C_hidden, K]
            # dk_reshaped is of shape [B, I, K, H, C_hidden]
            dk_reshaped = permute_final_dims(dkT, (2, 0, 1))

        # [B, I, K, H * C_hidden]
        dk_reshaped = dk_reshaped.flatten(start_dim=-2)

        # kv_x is broadcasted to perform linear layer to get v and k
        # so the gradients of kv_x need to be summed up from both contributions
        dkv_x = torch.einsum("...z, zc -> ...c", dv_reshaped, weight_v) + torch.einsum(
            "...z, zc -> ...c", dk_reshaped, weight_k
        )

        dkv_x_recv = ring_comm.comm_dk_final.enqueue_to_dispatch(dkv_x)

        dweight_v = torch.einsum("...z, ...c -> zc", dv_reshaped, kv_x_recv)
        dweight_v_dtensor = DTensor.from_local(
            dweight_v,
            device_mesh=ctx.device_mesh_input,
            placements=ctx.placements_dweights,
            shape=ctx.weight_v_shape,
            stride=ctx.weight_v_stride,
        )

        dweight_k = torch.einsum("...z, ...c -> zc", dk_reshaped, kv_x_recv)
        dweight_k_dtensor = DTensor.from_local(
            dweight_k,
            device_mesh=ctx.device_mesh_input,
            placements=ctx.placements_dweights,
            shape=ctx.weight_k_shape,
            stride=ctx.weight_k_stride,
        )

        ring_comm.comm_dbias_final1.wait_until_finished()
        ring_comm.comm_dk_final.wait_until_finished()

        # Handle transpose back for ending mode gradients
        if ctx.mode == _Mode.Ending:
            # Transpose gradients from [*, J, I, C] back to [*, I, J, C]
            dq_x = dq_x.transpose(-2, -3).contiguous()
            dkv_x_recv = dkv_x_recv.transpose(-2, -3).contiguous()
            # Transpose dtriangle_bias from [*, 1, H, J, I] to [*, 1, H, I, J]
            dtriangle_bias = dtriangle_bias.transpose(-1, -2)

        # Convert gradients back to DTensors
        dq_x_dtensor = DTensor.from_local(
            dq_x,
            device_mesh=ctx.device_mesh_input,
            placements=ctx.placements_input,
            shape=ctx.input_shape,
            stride=ctx.input_stride,
        )
        dkv_x_recv_dtensor = DTensor.from_local(
            dkv_x_recv,
            device_mesh=ctx.device_mesh_input,
            placements=ctx.placements_input,
            shape=ctx.input_shape,
            stride=ctx.input_stride,
        )
        if ctx.triattn_backend_bwd in (TriAttnBackend.CUEQ, TriAttnBackend.REFERENCE):
            # Convert dtriangle_bias from [*, 1, H, I, J] back to [*, I, J, H] for output
            dtriangle_bias_reshaped = permute_final_dims(dtriangle_bias.squeeze(-4), (1, 2, 0)).contiguous()
        elif ctx.triattn_backend_bwd == TriAttnBackend.TRIFAST:
            # [B, H, I, J] -> [B, I, J, H]
            dtriangle_bias_reshaped = permute_final_dims(dtriangle_bias, (1, 2, 0)).contiguous()
        dtriangle_bias_dtensor = DTensor.from_local(
            dtriangle_bias_reshaped,
            device_mesh=ctx.device_mesh_input,
            placements=ctx.placements_input,
            shape=ctx.triangle_bias_shape,
            stride=ctx.triangle_bias_stride,
        )

        return (
            dq_x_dtensor,
            dkv_x_recv_dtensor,
            None,
            dtriangle_bias_dtensor,
            dweight_q_dtensor,
            dweight_k_dtensor,
            dweight_v_dtensor,
            None,
            None,
            None,
            None,
            None,
            None,
        )


class RingMultiHeadTriangleAttention(nn.Module):
    """
    Multi-head triangle attention using a ring communication pattern
    (see the underlying autograd.Function for detail explanation of the algorithm)
    """

    def __init__(
        self,
        layer: Attention,
        device_mesh: DeviceMesh,
        ring_comm: Ring2DCommTriAttn,
        inf: float = 1e9,
    ):
        """
        Args:
            layer:
                Serial Attention instance to convert to distributed version
            device_mesh:
                Device mesh for distributed computation across multiple GPUs
            ring_comm:
               The communication object to use for distributed computation with ring attention
            inf:
                Infinity value for mask bias computation
        """
        super().__init__()

        self.c_q = layer.c_q
        self.c_k = layer.c_k
        self.c_v = layer.c_v
        self.c_hidden = layer.c_hidden
        self.no_heads = layer.no_heads
        self.gating = layer.gating
        self.ring_comm = ring_comm
        self.device_mesh = device_mesh
        self.inf = inf

        # Convert linear layers to DTensor-based counterparts
        # linear_{q,k,v} mapped to LinearParamsReplicated
        self.linear_q = LinearParamsReplicated(layer.linear_q, device_mesh)
        self.linear_k = LinearParamsReplicated(layer.linear_k, device_mesh)
        self.linear_v = LinearParamsReplicated(layer.linear_v, device_mesh)

        # linear_{o,g} mapped to LinearParamsReplicatedNoAutoCastBF16
        self.linear_o = LinearParamsReplicatedNoAutoCastBF16(layer.linear_o, device_mesh)

        self.linear_g = None
        if self.gating and layer.linear_g is not None:
            self.linear_g = LinearParamsReplicatedNoAutoCastBF16(layer.linear_g, device_mesh)

    def forward(
        self,
        q_x: DTensor,
        kv_x: DTensor,
        biases: list[DTensor],
        triattn_backend: TriAttnBackend = TriAttnBackend.REFERENCE,
    ) -> DTensor:
        """
        Args:
            q_x:
                [*, Q, C_q] query data
            kv_x:
                [*, K, C_k] key data
            biases:
                List containing mask (can be None) and triangle_bias
            triattn_backend:
                Triangular attention backend to use
        Returns
            [*, Q, C_q] attention update
        """
        # Linear layer weights are already DTensors from LinearParamsReplicated
        # compute q, k and v and launch initial shifting
        # of biases
        # kT == k.T is returned
        # Handle optional mask
        mask = biases[0] if biases[0] is not None else None
        triangle_bias = biases[1]

        o = _RingMultiHeadTriangleAttentionImpl.apply(
            q_x,
            kv_x,
            mask,
            triangle_bias,
            self.linear_q.weight,
            self.linear_k.weight,
            self.linear_v.weight,
            self.no_heads,
            self.c_hidden,
            self.ring_comm,
            self.inf,
            triattn_backend,
        )

        if self.linear_g is not None:
            # [B, I, J, H * C_hidden]
            g = self.linear_g(q_x)
            o = sigmoid_gate(o, g)

        # [*, Q, C_q]
        o = self.linear_o(o)

        return o


class _Mode(Enum):
    Starting = auto()
    Ending = auto()


class TriangleAttention(nn.Module):
    """Distributed triangle attention layer.

    This layer implements a distributed version of the triangle attention operation,
    which is used in attention mechanisms for protein structure prediction and other applications
    requiring pairwise feature interactions.

    The layer performs the following operations:
    1. Layer normalization of input pairwise features
    2. Linear projection to create triangle bias
    3. Distributed triangle attention computation using ring communication
    4. Transpose operations for ending node configuration

    Parameters
    ----------
    mode : _Mode
        Whether this is a starting or ending triangle attention node.
    layer : SerialTriangleAttentionStartingNode | SerialTriangleAttentionEndingNode
        The serial triangle attention layer to convert to distributed version.
        Used to initialize weights and normalization parameters.
    device_mesh : DeviceMesh
        The device mesh for distributed computation across multiple GPUs.
    comm : Ring2DCommTriAttn
        Ring communication object for efficient distributed triangle attention computation.
    """

    def __init__(
        self,
        mode: _Mode,
        layer: SerialTriangleAttentionStartingNode | SerialTriangleAttentionEndingNode,
        device_mesh: DeviceMesh,
        comm: Ring2DCommTriAttn,
    ) -> None:
        """Initialize the distributed triangle attention layer."""
        super().__init__()
        self.device_mesh = device_mesh
        self.ring_comm = comm
        self.mode = mode

        # Store layer parameters for distributed computation
        self.c_in = layer.c_in
        self.c_hidden = layer.c_hidden
        self.no_heads = layer.no_heads

        self.inf = layer.inf

        # Replicate parameters across the device mesh
        self.layer_norm = LayerNormParamsReplicatedNoAutoCastBF16(layer.layer_norm, self.device_mesh)
        self.linear = LinearParamsReplicatedNoAutoCastBF16(layer.linear, self.device_mesh)

        # Use the ring-based multi-head attention for distributed computation
        self.mha = RingMultiHeadTriangleAttention(
            layer.mha,
            self.device_mesh,
            self.ring_comm,
            self.inf,
        )

        # Validate mode consistency with ring comm
        if mode == _Mode.Starting:
            if not isinstance(layer, SerialTriangleAttentionStartingNode):
                raise ValueError(f"StartingNode mode is inconsistent with layer type {type(layer)}")
            if self.ring_comm.axis_cp != 1:
                raise ValueError(f"StartingNode mode is inconsistent with ring_comm.axis_cp {self.ring_comm.axis_cp}")
        elif mode == _Mode.Ending:
            if not isinstance(layer, SerialTriangleAttentionEndingNode):
                raise ValueError(f"EndingNode mode is inconsistent with layer type {type(layer)}")
            if self.ring_comm.axis_cp != 0:
                raise ValueError(f"EndingNode mode is inconsistent with ring_comm.axis_cp {self.ring_comm.axis_cp}")
        else:
            raise ValueError(f"Invalid mode {mode}")

    def forward(
        self, x: DTensor, mask: Optional[DTensor] = None, triattn_backend: TriAttnBackend = TriAttnBackend.REFERENCE
    ) -> DTensor:
        """Forward pass of the distributed triangle attention layer.

        Parameters
        ----------
        x : DTensor
            Input pairwise tensor with shape (B, I, J, C_in).
            Must be sharded on dimensions 1 and 2.
        mask : DTensor, optional
            Mask tensor with shape (B, I, J) indicating valid positions.
            Must be sharded on dimensions 1 and 2. If None, creates a mask of all ones.
        triattn_backend : TriAttnBackend
            Triangular attention backend to use
        Returns
        -------
        DTensor
            Output pairwise tensor with shape (B, I, J, C_in).
        """
        # Validate input types
        if not isinstance(x, DTensor):
            raise TypeError(f"Input 'x' must be of type DTensor. Got type {type(x)}.")

        if mask is not None:
            if not isinstance(mask, DTensor):
                raise TypeError(f"Input 'mask' must be of type DTensor or None. Got type {type(mask)}.")
            if mask.shape != x.shape[:-1]:
                raise ValueError(
                    f"Input tensor 'mask' must have the same shape as the first 3 dimensions of 'x'. "
                    f"Got mask shape: {mask.shape} vs x shape[:3]: {x.shape[:3]}"
                )

        if triattn_backend in (TriAttnBackend.CUEQ, TriAttnBackend.CUEQ_FWD_TRIFAST_BWD) and not cueq_is_installed:
            raise ValueError(
                "cuequivariance_torch is not installed. For Triangle Attention support, "
                "install using: pip install cuequivariance_ops_torch_cu13==<version> cuequivariance_torch==<version> "
                "where the 'version' tag can be found in the pyproject.toml file"
            )
        if (
            triattn_backend in (TriAttnBackend.TRIFAST, TriAttnBackend.CUEQ_FWD_TRIFAST_BWD)
            and not trifast_is_installed
        ):
            raise ValueError(
                "trifast is not installed. For Triangle Attention support, install using: pip install trifast"
            )
        if triattn_backend == TriAttnBackend.CUEQ_FWD_TRIFAST_BWD and x.dtype != torch.float32:
            raise ValueError(f"CUEQ_FWD_TRIFAST_BWD is only intended for FP32 usage. Got x.dtype {x.dtype}")

        # Normalize input - mask creation moved to MHA implementation
        x = self.layer_norm(x)

        # Compute triangle bias
        triangle_bias = self.linear(x)

        # Prepare biases for attention computation - mask_bias computation moved to MHA implementation
        # Regardless of triattn_backend, the binary mask is passed to the MHA implementation
        # where the default torch variant will convert it internally to a mask bias while
        # the underlying cueq call will use the binary mask directly
        biases = [mask, triangle_bias]

        # Apply distributed multi-head attention - transpose logic moved inside MHA implementation
        output = self.mha(q_x=x, kv_x=x, biases=biases, triattn_backend=triattn_backend)

        return output


class TriangleAttentionStartingNode(TriangleAttention):
    """Distributed triangle attention starting node layer."""

    def __init__(
        self,
        layer: SerialTriangleAttentionStartingNode,
        device_mesh: DeviceMesh,
        comm: Ring2DCommTriAttn,
    ) -> None:
        """Initialize the distributed triangle attention starting node layer.

        Parameters
        ----------
        layer : SerialTriangleAttentionStartingNode
            The serial triangle attention layer to convert to distributed version.
        device_mesh : DeviceMesh
            The device mesh for distributed computation across multiple GPUs.
        comm : Ring2DCommTriAttn
            Ring communication object for efficient distributed triangle attention computation.
        """
        if not layer.starting:
            raise ValueError("Serial layer must be configured as starting=True for TriangleAttentionStartingNode")
        super().__init__(_Mode.Starting, layer, device_mesh, comm)


class TriangleAttentionEndingNode(TriangleAttention):
    """Distributed triangle attention ending node layer."""

    def __init__(
        self,
        layer: SerialTriangleAttentionEndingNode,
        device_mesh: DeviceMesh,
        comm: Ring2DCommTriAttn,
    ) -> None:
        """Initialize the distributed triangle attention ending node layer.

        Parameters
        ----------
        layer : SerialTriangleAttentionEndingNode
            The serial triangle attention layer to convert to distributed version.
        device_mesh : DeviceMesh
            The device mesh for distributed computation across multiple GPUs.
        comm : Ring2DCommTriAttn
            Ring communication object for efficient distributed triangle attention computation.
        """
        if layer.starting:
            raise ValueError("Serial layer must be configured as starting=False for TriangleAttentionEndingNode")
        super().__init__(_Mode.Ending, layer, device_mesh, comm)
