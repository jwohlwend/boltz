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


from typing import Optional, Tuple, Union

import torch
from torch import nn
from torch.distributed.tensor import DTensor
from torch.utils.checkpoint import checkpoint

from boltz.distributed.comm import AttentionPairBiasComm, Ring2DComm, Ring2DCommTriAttn
from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.layers.attention import AttentionPairBias
from boltz.distributed.model.layers.dropout import apply_dropout_mask_msa_or_pair
from boltz.distributed.model.layers.elementwise_op import ElementwiseOp, elementwise_op
from boltz.distributed.model.layers.layernorm import LayerNormParamsReplicated
from boltz.distributed.model.layers.transition import Transition
from boltz.distributed.model.layers.triangular_attention import (
    TriangleAttentionEndingNode,
    TriangleAttentionStartingNode,
)
from boltz.distributed.model.layers.triangular_mult import (
    TriangleMultiplicationIncoming,
    TriangleMultiplicationOutgoing,
)
from boltz.distributed.model.modules.utils import TriAttnBackend, get_cpu_offload_context
from boltz.model.layers.pairformer import (
    PairformerLayer as SerialPairformerLayer,
)
from boltz.model.layers.pairformer import (
    PairformerModule as SerialPairformerModule,
)
from boltz.model.layers.pairformer import (
    PairformerNoSeqLayer as SerialPairformerNoSeqLayer,
)
from boltz.model.layers.pairformer import (
    PairformerNoSeqModule as SerialPairformerNoSeqModule,
)


class PairformerLayer(nn.Module):
    """Distributed PairformerLayer using DTensor (V2: attention with k_in).
    When self.no_seq=True, only the pairwise stack is applied (no sequence track).
    """

    def __init__(
        self,
        layer: Union[SerialPairformerLayer, SerialPairformerNoSeqLayer],
        dist_manager: DistributedManager,
    ) -> None:
        """Initialize the distributed PairformerLayer module.

        Parameters
        ----------
        layer : SerialPairformerLayer or SerialPairformerNoSeqLayer
            The serial Pairformer layer to be distributed.  Passing a
            ``SerialPairformerNoSeqLayer`` activates pair-only mode automatically.
        dist_manager : DistributedManager
            Distributed manager defining the distributed computation topology and groups.
        """
        if not isinstance(layer, (SerialPairformerLayer, SerialPairformerNoSeqLayer)):
            raise TypeError(
                f"layer must be SerialPairformerLayer or SerialPairformerNoSeqLayer, got {type(layer).__name__}"
            )
        super().__init__()
        self.dist_manager = dist_manager
        self.device_mesh = dist_manager.device_mesh_subgroups
        self.no_seq = isinstance(layer, SerialPairformerNoSeqLayer)

        self.token_z = layer.token_z
        self.dropout = layer.dropout
        self.post_layer_norm = layer.post_layer_norm

        # Mutable backend selection for triangle attention layers.
        # Default is REFERENCE (no fused kernels). To switch backend for the
        # entire model, use ``model.apply(SetTriAttnBackend(backend))``
        # (see boltz.distributed.model.modules.utils.SetTriAttnBackend).
        self.triattn_backend = TriAttnBackend.REFERENCE

        # Ring comms for triangular layers
        ring_comm_2d_trimul_outgoing = Ring2DComm(
            self.dist_manager.group["cp"],
            self.dist_manager.subgroups["cp"][0],
            self.dist_manager.layout_subgroups["cp"],
        )
        ring_comm_2d_trimul_incoming = Ring2DComm(
            self.dist_manager.group["cp"],
            self.dist_manager.subgroups["cp"][0],
            self.dist_manager.layout_subgroups["cp"],
        )
        ring_comm_tri_attn_start = Ring2DCommTriAttn(
            self.dist_manager.group["cp"],
            self.dist_manager.layout_subgroups["cp"],
            1,
        )
        ring_comm_tri_attn_end = Ring2DCommTriAttn(
            self.dist_manager.group["cp"],
            self.dist_manager.layout_subgroups["cp"],
            0,
        )

        self.tri_mul_out = TriangleMultiplicationOutgoing(
            layer.tri_mul_out, self.device_mesh, ring_comm_2d_trimul_outgoing
        )
        self.tri_mul_in = TriangleMultiplicationIncoming(
            layer.tri_mul_in, self.device_mesh, ring_comm_2d_trimul_incoming
        )
        self.tri_att_start = TriangleAttentionStartingNode(
            layer.tri_att_start, self.device_mesh, ring_comm_tri_attn_start
        )
        self.tri_att_end = TriangleAttentionEndingNode(layer.tri_att_end, self.device_mesh, ring_comm_tri_attn_end)
        self.transition_z = Transition(layer.transition_z, self.device_mesh)

        if not self.no_seq:
            self.num_heads = layer.num_heads
            # Detect V1 vs V2 from the serial attention module.
            # V1 AttentionPairBias uses initial_norm=True which creates norm_s
            # (see src/boltz/model/layers/attention.py AttentionPairBias.__init__).
            # V2 does not have norm_s.  Prefer the explicit .v2 attribute when
            # set by the serial PairformerLayer; otherwise infer from norm_s.
            if hasattr(layer, "v2"):
                self.v2 = layer.v2
            else:
                self.v2 = not hasattr(layer.attention, "norm_s")
            self.pre_norm_s = LayerNormParamsReplicated(layer.pre_norm_s, self.device_mesh)
            attention_pair_bias_comm = AttentionPairBiasComm(
                self.dist_manager.group["cp"],
                self.dist_manager.layout_subgroups["cp"],
                self.dist_manager.subgroups["cp"][0],
                self.dist_manager.subgroups["cp"][1],
            )
            if self.v2:
                self.attention = AttentionPairBias(
                    layer.attention,
                    self.device_mesh,
                    attention_pair_bias_comm,
                    apply_initial_norm=False,
                    compute_pair_bias=True,
                    use_model_cache=False,
                )
            else:
                self.attention = AttentionPairBias(
                    layer.attention,
                    self.device_mesh,
                    attention_pair_bias_comm,
                    apply_initial_norm=True,
                    compute_pair_bias=True,
                    use_model_cache=True,
                )
            self.transition_s = Transition(layer.transition_s, self.device_mesh)
            if self.post_layer_norm:
                self.s_post_norm = LayerNormParamsReplicated(layer.s_post_norm, self.device_mesh)
            else:
                self.s_post_norm = None

    def forward(
        self,
        s: Optional[DTensor] = None,
        z: Optional[DTensor] = None,
        mask: Optional[DTensor] = None,
        pair_mask: Optional[DTensor] = None,
    ) -> Union[Tuple[DTensor, DTensor], DTensor]:
        """Forward pass. Pass s= and mask= for the full pairformer; omit them for pair-only."""
        assert z is not None and pair_mask is not None
        z = elementwise_op(
            z,
            apply_dropout_mask_msa_or_pair(
                self.tri_mul_out(z, mask=pair_mask),
                self.dropout,
                self.training,
            ),
            ElementwiseOp.SUM,
        )
        z = elementwise_op(
            z,
            apply_dropout_mask_msa_or_pair(
                self.tri_mul_in(z, mask=pair_mask),
                self.dropout,
                self.training,
            ),
            ElementwiseOp.SUM,
        )
        z = elementwise_op(
            z,
            apply_dropout_mask_msa_or_pair(
                self.tri_att_start(z, mask=pair_mask, triattn_backend=self.triattn_backend),
                self.dropout,
                self.training,
            ),
            ElementwiseOp.SUM,
        )
        z = elementwise_op(
            z,
            apply_dropout_mask_msa_or_pair(
                self.tri_att_end(z, mask=pair_mask, triattn_backend=self.triattn_backend),
                self.dropout,
                self.training,
                columnwise=True,
            ),
            ElementwiseOp.SUM,
        )
        z = elementwise_op(z, self.transition_z(z), ElementwiseOp.SUM)
        if s is None:
            return z
        assert mask is not None
        with torch.autocast("cuda", enabled=False):
            safe_dtype = torch.promote_types(s.dtype, torch.float32)
            s = s.to(dtype=safe_dtype)
            s_normed = self.pre_norm_s(s)
            s = elementwise_op(
                s,
                self.attention(s=s_normed, z=z.to(dtype=safe_dtype), mask=mask.to(dtype=safe_dtype), k_in=s_normed),
                ElementwiseOp.SUM,
            )
            s = elementwise_op(s, self.transition_s(s), ElementwiseOp.SUM)
            s = self.s_post_norm(s) if self.s_post_norm is not None else s
        return s, z


class PairformerModule(nn.Module):
    """Distributed PairformerModule using DTensor (V2).

    Handles both full pairformer (sequence + pairwise stacks) and pair-only
    mode. The mode is inferred from the serial module type: passing a
    ``SerialPairformerNoSeqModule`` activates pair-only mode automatically.
    """

    def __init__(
        self,
        module: Union[SerialPairformerModule, SerialPairformerNoSeqModule],
        dist_manager: DistributedManager,
        cpu_offloading: bool = False,
    ) -> None:
        """Initialize the distributed PairformerModule module.

        Parameters
        ----------
        module : SerialPairformerModule or SerialPairformerNoSeqModule
            The serial Pairformer module containing weights and configuration to be distributed.
        dist_manager : DistributedManager
            Distributed manager defining the distributed computation topology and groups.
        cpu_offloading : bool, optional
            Whether to offload checkpoint-boundary activations to CPU when
            activation checkpointing is enabled.  This is a distributed-only
            option (the serial Boltz-2 MSAModule does not support it).
            Defaults to False.
        """
        super().__init__()
        self.dist_manager = dist_manager
        self.device_mesh = dist_manager.device_mesh_subgroups

        no_seq = isinstance(module, SerialPairformerNoSeqModule)
        self.no_seq = no_seq

        self.token_z = module.token_z
        self.num_blocks = module.num_blocks
        self.dropout = module.dropout
        self.post_layer_norm = module.post_layer_norm
        self.activation_checkpointing = module.activation_checkpointing
        self.cpu_offloading = cpu_offloading

        if not no_seq:
            self.num_heads = module.num_heads

        self.layers = nn.ModuleList()
        for serial_layer in module.layers:
            self.layers.append(PairformerLayer(serial_layer, dist_manager))

    def forward(
        self,
        s: Optional[DTensor] = None,
        z: Optional[DTensor] = None,
        mask: Optional[DTensor] = None,
        pair_mask: Optional[DTensor] = None,
    ) -> Union[Tuple[DTensor, DTensor], DTensor]:
        """Forward pass. Pass s= and mask= for the full pairformer; omit them for pair-only."""
        if self.activation_checkpointing and self.training:
            if self.cpu_offloading:
                with get_cpu_offload_context(optimized=True):
                    for layer in self.layers:
                        result = checkpoint(
                            layer,
                            s,
                            z,
                            mask,
                            pair_mask,
                            use_reentrant=False,
                        )
                        if self.no_seq:
                            z = result
                        else:
                            s, z = result
            else:
                for layer in self.layers:
                    result = checkpoint(
                        layer,
                        s,
                        z,
                        mask,
                        pair_mask,
                        use_reentrant=False,
                    )
                    if self.no_seq:
                        z = result
                    else:
                        s, z = result
        else:
            for layer in self.layers:
                result = layer(
                    s=s,
                    z=z,
                    mask=mask,
                    pair_mask=pair_mask,
                )
                if self.no_seq:
                    z = result
                else:
                    s, z = result
        return z if self.no_seq else (s, z)


class PairformerNoSeqLayer(PairformerLayer):
    """Distributed PairformerNoSeqLayer (pairwise stack only, no sequence track)."""

    pass


class PairformerNoSeqModule(PairformerModule):
    """Distributed PairformerNoSeqModule (pairwise stack only)."""

    pass
