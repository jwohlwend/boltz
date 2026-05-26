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

"""DTensor-compatible DiffusionConditioning module for Context Parallelism.

V2-only module that pre-computes conditioning (pair features, atom encoder,
bias projections) outside the diffusion loop. This is a Boltz-2-specific
refactoring that moves conditioning from inside the diffusion step to a
one-time pre-computation.
"""

from torch import nn
from torch.distributed.device_mesh import DeviceMesh
from torch.distributed.tensor import DTensor
from torch.nn import Module

from boltz.distributed.data.feature.featurizer import pack_atom_features
from boltz.distributed.model.layers.cat_and_chunk import shardwise_cat
from boltz.distributed.model.layers.layernorm import LayerNormParamsReplicated
from boltz.distributed.model.layers.linear import LinearParamsReplicated
from boltz.distributed.model.modules.encoders import AtomEncoder, PairwiseConditioning
from boltz.model.modules.diffusion_conditioning import DiffusionConditioning as SerialDiffusionConditioning


class DiffusionConditioning(Module):
    """DTensor DiffusionConditioning for Context Parallelism (V2 only).

    Pre-computes:
    - Pairwise conditioning (z)
    - Atom encoder (q, c, p)
    - Bias projections for atom encoder, atom decoder, and token transformer layers
    """

    def __init__(
        self,
        layer: SerialDiffusionConditioning,
        device_mesh: DeviceMesh,
    ):
        """Initialize the DTensor DiffusionConditioning.

        Parameters
        ----------
        layer : SerialDiffusionConditioning
            The serial DiffusionConditioning module (V2 only).
        device_mesh : DeviceMesh
            The device mesh for distributed tensor operations.

        """
        super().__init__()
        assert isinstance(
            layer, SerialDiffusionConditioning
        ), f"Expected SerialDiffusionConditioning, got {type(layer)}"
        self.device_mesh = device_mesh
        self.atoms_per_window_queries = layer.atom_encoder.atoms_per_window_queries

        # Wrap child modules
        self.pairwise_conditioner = PairwiseConditioning(
            layer=layer.pairwise_conditioner,
            device_mesh=device_mesh,
        )
        self.atom_encoder = AtomEncoder(
            layer=layer.atom_encoder,
            device_mesh=device_mesh,
        )

        # Bias projection layers: ModuleList of Sequential(LayerNorm, Linear)
        self.atom_enc_proj_z = nn.ModuleList()
        for serial_seq in layer.atom_enc_proj_z:
            self.atom_enc_proj_z.append(
                nn.Sequential(
                    LayerNormParamsReplicated(serial_seq[0], device_mesh),
                    LinearParamsReplicated(serial_seq[1], device_mesh),
                )
            )

        self.atom_dec_proj_z = nn.ModuleList()
        for serial_seq in layer.atom_dec_proj_z:
            self.atom_dec_proj_z.append(
                nn.Sequential(
                    LayerNormParamsReplicated(serial_seq[0], device_mesh),
                    LinearParamsReplicated(serial_seq[1], device_mesh),
                )
            )

        self.token_trans_proj_z = nn.ModuleList()
        for serial_seq in layer.token_trans_proj_z:
            self.token_trans_proj_z.append(
                nn.Sequential(
                    LayerNormParamsReplicated(serial_seq[0], device_mesh),
                    LinearParamsReplicated(serial_seq[1], device_mesh),
                )
            )

    def forward(
        self,
        s_trunk: DTensor,
        z_trunk: DTensor,
        relative_position_encoding: DTensor,
        feats: dict[str, DTensor],
    ) -> tuple[DTensor, DTensor, DTensor, DTensor, DTensor]:
        """Forward pass of the DTensor DiffusionConditioning.

        Parameters
        ----------
        s_trunk : DTensor
            Token single representation with shape (B, N, token_s).
            Placements: (Shard(0), Shard(1), Replicate()).
        z_trunk : DTensor
            Token pair representation with shape (B, N, N, token_z).
            Placements: (Shard(0), Shard(1), Shard(2)).
        relative_position_encoding : DTensor
            Relative position encoding with shape (B, N, N, token_z).
            Placements: (Shard(0), Shard(1), Shard(2)).
        feats : dict[str, DTensor]
            Unpacked atom features (with intersperse padding from the CP DTensor
            data loader). The module calls ``pack_atom_features`` internally.

        Returns
        -------
        tuple[DTensor, DTensor, DTensor, DTensor, DTensor]
            q : DTensor with shape (B, N_atoms_packed, atom_s).
                Placements: (Shard(0), Shard(1), Replicate()).
            c : DTensor with shape (B, N_atoms_packed, atom_s).
                Placements: (Shard(0), Shard(1), Replicate()).
            atom_enc_bias : DTensor with shape (B, K, W, H, total_atom_enc_heads).
                Placements: (Shard(0), Shard(1), Replicate()).
                Window-batched atom pair repr; K is sharded, W/H are local window dims.
            atom_dec_bias : DTensor with shape (B, K, W, H, total_atom_dec_heads).
                Placements: (Shard(0), Shard(1), Replicate()).
                Window-batched atom pair repr; K is sharded, W/H are local window dims.
            token_trans_bias : DTensor with shape (B, N, N, total_token_trans_heads).
                Placements: (Shard(0), Shard(1), Shard(2)).
                Token pair repr; both token dims are sharded.

        """
        # Atom features (feats) are expected in unpacked layout (with intersperse padding
        # from the CP DTensor data loader). Each module calls pack_atom_features internally
        # to form a self-contained pack/unpack closure, so that no external caller needs to
        # pre-pack features. This ensures all modules accept atom features directly as
        # produced by the data loader, and future refactoring of the data loading pipeline
        # will not require changes to these modules.
        _keys_atom_features_packed = {
            "atom_pad_mask",
            "ref_pos",
            "ref_space_uid",
            "ref_charge",
            "ref_element",
            "ref_atom_name_chars",
            "atom_to_token",
        }
        feats_packed = pack_atom_features(feats, _keys_atom_features_packed, self.atoms_per_window_queries)

        # Pairwise conditioning
        z = self.pairwise_conditioner(z_trunk, relative_position_encoding)

        # Atom encoder: q, c, p (all in packed layout)
        q, c, p = self.atom_encoder(feats=feats_packed, s_trunk=s_trunk, z=z)

        # Atom encoder bias projections: project p (window-batched atom pair) through each layer, concatenate
        # p: (B, K, W, H, atom_z) with placements (S(0), S(1), R) — K sharded, W/H local window dims
        atom_enc_bias_list = []
        for proj in self.atom_enc_proj_z:
            atom_enc_bias_list.append(proj(p))
        atom_enc_bias = shardwise_cat(atom_enc_bias_list, dim=-1)

        # Atom decoder bias projections
        atom_dec_bias_list = []
        for proj in self.atom_dec_proj_z:
            atom_dec_bias_list.append(proj(p))
        atom_dec_bias = shardwise_cat(atom_dec_bias_list, dim=-1)

        # Token transformer bias projections: project z (token pair)
        # z: (B, N, N, token_z) with placements (S(0), S(1), S(2))
        token_trans_bias_list = []
        for proj in self.token_trans_proj_z:
            token_trans_bias_list.append(proj(z))
        token_trans_bias = shardwise_cat(token_trans_bias_list, dim=-1)

        return q, c, atom_enc_bias, atom_dec_bias, token_trans_bias
