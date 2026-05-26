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

"""DTensor-based implementation of Boltz-2 trunk modules.

This module provides distributed implementations of the MSAModule, MSALayer,
InputEmbedder, DistogramModule, BFactorModule, and ContactConditioning for Boltz-2.

MSAModule and MSALayer use PairformerNoSeq layers for triangle attention
and multiplication.  DistogramModule handles 5D output with num_distograms.
BFactorModule and ContactConditioning are Boltz-2–specific head modules.
InputEmbedder wraps the serial InputEmbedder for context parallelism using
distributed AtomEncoder and AtomAttentionEncoder.
"""

from math import pi
from typing import Dict, Optional, Tuple

import torch
from torch import nn
from torch.distributed.device_mesh import DeviceMesh
from torch.distributed.tensor import DTensor, Replicate, Shard, distribute_tensor
from torch.utils.checkpoint import checkpoint

from boltz.data import const
from boltz.distributed.comm import Ring2DComm, TransposeComm
from boltz.distributed.data.feature.featurizer import pack_atom_features
from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.layers.cat_and_chunk import shardwise_cat
from boltz.distributed.model.layers.clip import clip
from boltz.distributed.model.layers.dropout import apply_dropout_mask_msa_or_pair
from boltz.distributed.model.layers.elementwise_op import ElementwiseOp, elementwise_op
from boltz.distributed.model.layers.embedding import EmbeddingParamsReplicated
from boltz.distributed.model.layers.flatten_and_unflatten import shardwise_unflatten
from boltz.distributed.model.layers.layernorm import LayerNormParamsReplicated
from boltz.distributed.model.layers.linear import LinearParamsReplicated
from boltz.distributed.model.layers.outer_product_mean import OuterProductMean
from boltz.distributed.model.layers.pair_averaging import PairWeightedAveraging, Ring2DCommPairAveraging
from boltz.distributed.model.layers.pairformer import PairformerNoSeqLayer
from boltz.distributed.model.layers.redistribute_transpose import redistribute_transpose
from boltz.distributed.model.layers.replicate_op import ReplicateOp, replicate_op
from boltz.distributed.model.layers.shardwise_op import shardwise_one_hot
from boltz.distributed.model.layers.squeeze import shardwise_unsqueeze
from boltz.distributed.model.layers.transition import Transition
from boltz.distributed.model.modules.encoders import AtomAttentionEncoder as DistAtomAttentionEncoder
from boltz.distributed.model.modules.encoders import AtomEncoder as DistAtomEncoder
from boltz.distributed.model.modules.utils import get_cpu_offload_context
from boltz.distributed.utils import update_exhaustive_strides
from boltz.model.modules.encodersv2 import FourierEmbedding as SerialFourierEmbedding
from boltz.model.modules.trunkv2 import BFactorModule as SerialBFactorModule
from boltz.model.modules.trunkv2 import ContactConditioning as SerialContactConditioning
from boltz.model.modules.trunkv2 import DistogramModule as SerialDistogramModule
from boltz.model.modules.trunkv2 import InputEmbedder as SerialInputEmbedder
from boltz.model.modules.trunkv2 import MSALayer as SerialMSALayer
from boltz.model.modules.trunkv2 import MSAModule as SerialMSAModule


class InputEmbedder(nn.Module):
    """DTensor InputEmbedder for Boltz-2 using window-batching atom attention.

    Wraps the serial Boltz-2 ``InputEmbedder`` and distributes its sub-modules
    for context parallelism.  Atom-level operations delegate to the existing
    distributed :class:`AtomEncoder` and :class:`AtomAttentionEncoder` (window-
    batching variants) while token-level projections use
    :class:`LinearParamsReplicated` / :class:`EmbeddingParamsReplicated`.

    Atom feature packing (via :func:`pack_atom_features`) is internalized in
    :meth:`forward` so that callers only need to pass the raw distributed atom
    features.  The packed features are discarded after the forward pass.
    """

    _KEYS_ATOM_FEATURES_PACKED = {
        "atom_pad_mask",
        "ref_pos",
        "ref_space_uid",
        "ref_charge",
        "ref_element",
        "ref_atom_name_chars",
        "atom_to_token",
    }

    def __init__(
        self,
        module: SerialInputEmbedder,
        device_mesh: DeviceMesh,
    ) -> None:
        """Initialize the distributed InputEmbedder for Boltz-2.

        Parameters
        ----------
        module : SerialInputEmbedder
            The serial InputEmbedder containing weights and configuration.
        device_mesh : DeviceMesh
            Device mesh defining the distributed computation topology.
        """
        super().__init__()

        if not isinstance(module, SerialInputEmbedder):
            raise TypeError(f"Expected SerialInputEmbedder, got {type(module)}")

        self.add_method_conditioning = module.add_method_conditioning
        self.add_modified_flag = module.add_modified_flag
        self.add_cyclic_flag = module.add_cyclic_flag
        self.add_mol_type_feat = module.add_mol_type_feat

        # Atom-level modules -- delegate to existing distributed implementations
        self.atom_encoder = DistAtomEncoder(layer=module.atom_encoder, device_mesh=device_mesh)
        self.atoms_per_window_queries = self.atom_encoder.atoms_per_window_queries

        # atom_enc_proj_z: Sequential(LayerNorm, Linear) projects pair repr to
        # attention bias.  Operates on the last (replicated) dim of p so
        # LayerNormParamsReplicated + LinearParamsReplicated work directly.
        self.atom_enc_proj_z = nn.Sequential(
            LayerNormParamsReplicated(module.atom_enc_proj_z[0], device_mesh=device_mesh),
            LinearParamsReplicated(layer_local=module.atom_enc_proj_z[1], device_mesh=device_mesh),
        )

        self.atom_attention_encoder = DistAtomAttentionEncoder(
            layer=module.atom_attention_encoder, device_mesh=device_mesh
        )

        # Token-level projections (replicated parameters, local ops)
        self.res_type_encoding = LinearParamsReplicated(layer_local=module.res_type_encoding, device_mesh=device_mesh)
        self.msa_profile_encoding = LinearParamsReplicated(
            layer_local=module.msa_profile_encoding, device_mesh=device_mesh
        )

        # Optional conditioning modules
        if self.add_method_conditioning:
            self.method_conditioning_init = EmbeddingParamsReplicated(
                module.method_conditioning_init, device_mesh=device_mesh
            )
        if self.add_modified_flag:
            self.modified_conditioning_init = EmbeddingParamsReplicated(
                module.modified_conditioning_init, device_mesh=device_mesh
            )
        if self.add_cyclic_flag:
            self.cyclic_conditioning_init = LinearParamsReplicated(
                layer_local=module.cyclic_conditioning_init, device_mesh=device_mesh
            )
        if self.add_mol_type_feat:
            self.mol_type_conditioning_init = EmbeddingParamsReplicated(
                module.mol_type_conditioning_init, device_mesh=device_mesh
            )

    def forward(self, feats: dict[str, DTensor], affinity: bool = False) -> DTensor:
        """Forward pass for the distributed InputEmbedder.

        Atom feature packing/unpacking is internalized: ``pack_atom_features``
        is called here to convert the raw distributed atom features into packed
        format for window-batching.  The packed features are discarded after
        the atom encoder and atom attention encoder calls.

        Parameters
        ----------
        feats : dict[str, DTensor]
            Input features (token-level and atom-level).
        affinity : bool, optional
            When True, use ``profile_affinity`` / ``deletion_mean_affinity``
            instead of ``profile`` / ``deletion_mean``.  Defaults to False.

        Returns
        -------
        DTensor
            Token-level single representation ``s`` with shape
            ``(B, N_tokens, token_s)`` and placement
            ``(Shard(0), Shard(1), Replicate())``.
        """
        # Token-level features
        # res_type is integer one-hot; cast to the layer's weight dtype so it
        # works with any precision (float32, float64, bfloat16, etc.).
        res_type = feats["res_type"].to(self.res_type_encoding.weight.dtype)
        if affinity:
            profile = feats["profile_affinity"]
            deletion_mean = shardwise_unsqueeze(feats["deletion_mean_affinity"], -1)
        else:
            profile = feats["profile"]
            deletion_mean = shardwise_unsqueeze(feats["deletion_mean"], -1)

        # Pack atom features for window batching.
        # pack_atom_features converts pad_and_scatter output (with per-shard
        # trailing padding) into packed format with global atom-to-token indices.
        # The packed features are discarded after the forward pass.
        feats_packed = pack_atom_features(feats, self._KEYS_ATOM_FEATURES_PACKED, self.atoms_per_window_queries)

        # Atom encoding: produces (q, c, p).
        # AtomEncoder.forward() internally wraps its body in
        # torch.autocast("cuda", enabled=False), matching the serial code.
        q, c, p = self.atom_encoder(feats_packed)

        # Project pair representation to attention bias
        atom_enc_bias = self.atom_enc_proj_z(p)

        # Atom attention encoding: produces (a, q, c, p)
        a, _, _, _ = self.atom_attention_encoder(
            feats=feats_packed,
            q=q,
            c=c,
            atom_enc_bias=atom_enc_bias,
        )

        # Token-level embedding: sum of atom attention output + learned projections
        profile_cat = shardwise_cat([profile, deletion_mean], dim=-1)

        s = elementwise_op(a, self.res_type_encoding(res_type), ElementwiseOp.SUM)
        s = elementwise_op(s, self.msa_profile_encoding(profile_cat), ElementwiseOp.SUM)

        # Optional conditioning
        if self.add_method_conditioning:
            s = elementwise_op(s, self.method_conditioning_init(feats["method_feature"]), ElementwiseOp.SUM)
        if self.add_modified_flag:
            s = elementwise_op(s, self.modified_conditioning_init(feats["modified"]), ElementwiseOp.SUM)
        if self.add_cyclic_flag:
            cyclic = feats["cyclic_period"].to(self.cyclic_conditioning_init.weight.dtype)
            cyclic = clip(cyclic, max_val=1.0)
            cyclic = shardwise_unsqueeze(cyclic, -1)
            s = elementwise_op(s, self.cyclic_conditioning_init(cyclic), ElementwiseOp.SUM)
        if self.add_mol_type_feat:
            s = elementwise_op(s, self.mol_type_conditioning_init(feats["mol_type"]), ElementwiseOp.SUM)

        return s


class MSALayer(nn.Module):
    """Distributed MSA layer for Boltz-2 using DTensor.

    This is the Boltz-2 version of MSALayer which uses PairformerNoSeqLayer
    for triangle operations instead of individual triangle multiplication and
    attention layers.

    Input/Output Placements:
    - z: (Shard(0), Shard(1), Shard(2)) - Pair representation
    - m: (Shard(0), Shard(1), Shard(2)) - MSA representation
    - token_mask: (Shard(0), Shard(1), Shard(2)) - Token pair mask
    - msa_mask: (Shard(0), Shard(1), Shard(2)) - MSA mask

    Communication:
    - PairWeightedAveraging: Ring communication for MSA-to-pair
    - OuterProductMean: Ring communication for outer product
    - PairformerNoSeqLayer: Triangle operations with ring communication
    """

    def __init__(
        self,
        layer: SerialMSALayer,
        dist_manager: DistributedManager,
    ) -> None:
        """Initialize the distributed MSALayer for Boltz-2.

        Parameters
        ----------
        layer : SerialMSALayer
            The serial MSA layer containing weights and configuration to be distributed.
        dist_manager : DistributedManager
            Distributed manager defining the distributed computation topology and groups.
        """
        super().__init__()
        self.dist_manager = dist_manager
        self.device_mesh = dist_manager.device_mesh_subgroups

        # Store dropout rates
        self.msa_dropout = layer.msa_dropout

        # Create communication objects for distributed computation
        ring_comm_2d_outer_product = Ring2DComm(
            self.dist_manager.group["cp"],
            self.dist_manager.subgroups["cp"][0],
            self.dist_manager.layout_subgroups["cp"],
        )

        ## PWA Implementation with Ring2DCommPairAveraging
        ring_comm_2d_pair_avg = Ring2DCommPairAveraging(
            self.dist_manager.group["cp"],
            self.dist_manager.subgroups["cp"][0],
            self.dist_manager.layout_subgroups["cp"],
        )
        self.pair_weighted_averaging = PairWeightedAveraging(
            layer.pair_weighted_averaging, self.device_mesh, ring_comm_2d_pair_avg
        )

        # Map serial layers to distributed versions
        self.msa_transition = Transition(layer.msa_transition, self.device_mesh)

        # Map PairformerNoSeqLayer to distributed version
        self.pairformer_layer = PairformerNoSeqLayer(layer.pairformer_layer, dist_manager)
        assert self.pairformer_layer.no_seq, (
            f"Expected no_seq=True for PairformerNoSeqLayer, " f"got no_seq={self.pairformer_layer.no_seq}"
        )

        self.outer_product_mean = OuterProductMean(
            layer.outer_product_mean, self.device_mesh, ring_comm_2d_outer_product
        )

    def forward(
        self,
        z: DTensor,
        m: DTensor,
        token_mask: DTensor,
        msa_mask: DTensor,
    ) -> Tuple[DTensor, DTensor]:
        """Perform the forward pass.

        Parameters
        ----------
        z : DTensor
            The pair representation with placement (Shard(0), Shard(1), Shard(2))
        m : DTensor
            The MSA representation with placement (Shard(0), Shard(1), Shard(2))
        token_mask : DTensor
            The token pair mask with placement (Shard(0), Shard(1), Shard(2))
        msa_mask : DTensor
            The MSA mask with placement (Shard(0), Shard(1), Shard(2))

        Returns
        -------
        Tuple[DTensor, DTensor]
            The updated pair representation and MSA representation.
        """
        # Communication to MSA stack
        m = elementwise_op(
            m,
            apply_dropout_mask_msa_or_pair(
                self.pair_weighted_averaging(m, z, token_mask), self.msa_dropout, self.training
            ),
            ElementwiseOp.SUM,
        )
        m = elementwise_op(m, self.msa_transition(m), ElementwiseOp.SUM)

        # Communication to pairwise stack via outer product
        z = elementwise_op(z, self.outer_product_mean(m, msa_mask), ElementwiseOp.SUM)

        # Compute pairwise stack using PairformerNoSeqLayer
        # Note: PairformerNoSeqLayer returns updated z directly (no residual connection needed)
        z = self.pairformer_layer(z=z, pair_mask=token_mask)

        return z, m


class MSAModule(nn.Module):
    """Distributed MSA module for Boltz-2 using DTensor.

    This is the Boltz-2 version of MSAModule which uses PairformerNoSeqLayer
    for triangle operations within each MSA layer.

    Input/Output Placements:
    - z: (Shard(0), Shard(1), Shard(2)) - Pair representation
    - emb: (Shard(0), Replicate(), Shard(1)) - Single representation
    - msa: (Shard(0), Shard(1), Shard(2)) - MSA sequences
    - msa_mask: (Shard(0), Shard(1), Shard(2)) - MSA mask
    - token_pair_pad_mask: (Shard(0), Shard(1), Shard(2)) - Pair mask

    Output:
    - z: (Shard(0), Shard(1), Shard(2)) - Updated pair representation

    Communication:
    - replicate_op: Broadcast emb to MSA dimension
    - PairWeightedAveraging: Ring communication
    - OuterProductMean: Ring communication
    - PairformerNoSeqLayer: Triangle operations with ring communication
    """

    def __init__(
        self,
        module: SerialMSAModule,
        dist_manager: DistributedManager,
        cpu_offloading: bool = False,
    ) -> None:
        """Initialize the distributed MSAModule for Boltz-2.

        Parameters
        ----------
        module : SerialMSAModule
            The serial MSA module containing weights and configuration to be distributed.
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

        # Store attributes from the serial module
        self.msa_blocks = module.msa_blocks
        self.msa_dropout = module.msa_dropout
        self.z_dropout = module.z_dropout
        self.use_paired_feature = module.use_paired_feature
        self.subsample_msa = module.subsample_msa
        self.num_subsampled_msa = module.num_subsampled_msa

        # CP does not support MSA subsampling at module/layer level; require serial config to have it disabled.
        if self.subsample_msa:
            raise NotImplementedError(
                "Subsampling MSA at module level is not supported with context parallelism. "
                "The serial MSAModule must be built with subsample_msa=False."
            )

        # Activation checkpointing is read from the serial module.
        self.activation_checkpointing = getattr(module, "activation_checkpointing", False)
        # CPU offloading is a distributed-only option (the serial Boltz-2 MSAModule
        # does not have this flag).  When enabled together with activation
        # checkpointing, checkpoint-boundary activations are moved to CPU during
        # the forward pass and restored on the backward pass, trading extra
        # CPU<->GPU transfers for reduced GPU memory.
        self.cpu_offloading = cpu_offloading

        # Map serial projections to distributed versions
        # Note: s_proj and msa_proj are LinearParamsReplicated since they operate on
        # features that will be broadcast/replicated
        self.s_proj = LinearParamsReplicated(module.s_proj, self.device_mesh)
        self.msa_proj = LinearParamsReplicated(module.msa_proj, self.device_mesh)

        # Map MSA layers to distributed versions
        self.layers = nn.ModuleList()
        for serial_layer in module.layers:
            self.layers.append(MSALayer(serial_layer, dist_manager))

    def forward(
        self,
        z: DTensor,
        emb: DTensor,
        feats: Dict[str, DTensor],
    ) -> DTensor:
        """Perform the forward pass.

        Parameters
        ----------
        z : DTensor
            The pairwise embeddings with placement (Shard(0), Shard(1), Shard(2))
        emb : DTensor
            The input embeddings with placement (Shard(0), Replicate(), Shard(1))
        feats : Dict[str, DTensor]
            Input features as DTensors

        Returns
        -------
        DTensor
            The output pairwise embeddings.
        """
        # Expected placements
        expected_msa_placement = (Shard(0), Shard(1), Shard(2))
        expected_emb_placement = (Shard(0), Replicate(), Shard(1))

        # Sanity check for emb placement
        if emb.placements != expected_emb_placement:
            raise ValueError(f"Expected emb placement {expected_emb_placement}, but got {emb.placements}")

        # Sanity check for z placement
        if z.placements != expected_msa_placement:
            raise ValueError(f"Expected z placement {expected_msa_placement}, but got {z.placements}")

        # Load relevant features – apply one-hot encoding to match the serial
        # MSAModule (see src/boltz/model/modules/trunkv2.py), then cast from
        # integer to the working dtype so it can be concatenated with the other
        # floating-point features.
        msa = feats["msa"]
        msa = shardwise_one_hot(msa, num_classes=const.num_tokens).to(dtype=z.dtype)
        has_deletion = shardwise_unsqueeze(feats["has_deletion"], -1)
        deletion_value = shardwise_unsqueeze(feats["deletion_value"], -1)
        msa_mask = feats["msa_mask"]
        token_mask = feats["token_pair_pad_mask"]

        # Compute MSA embeddings
        feats_to_cat = [msa, has_deletion, deletion_value]
        if self.use_paired_feature:
            is_paired = shardwise_unsqueeze(feats["msa_paired"], -1)
            feats_to_cat.append(is_paired)

        # Sanity check for feature DTensor placements
        for feat in feats_to_cat:
            if feat.placements != expected_msa_placement:
                raise ValueError(f"Expected MSA feature placement {expected_msa_placement}, but got {feat.placements}")

        # Concatenate MSA features
        m = shardwise_cat(feats_to_cat, dim=-1)

        # Compute input projections
        m = self.msa_proj(m)
        emb_proj = self.s_proj(emb)

        # Use DTensor replicate_op to add emb to MSA
        # emb_proj has placement (Shard(0), Replicate(), Shard(1))
        # m has placement (Shard(0), Shard(1), Shard(2))
        # We need to broadcast emb_proj along the MSA dimension (dim=1)
        m = replicate_op(m, emb_proj, 1, op=ReplicateOp.ADD)

        # Perform MSA blocks.
        # When activation_checkpointing is enabled, saved activations are recomputed
        # during the backward pass.  When cpu_offloading is additionally enabled,
        # the checkpoint-boundary tensors are moved to CPU (module-level offloading)
        # via get_cpu_offload_context, reducing GPU memory at the cost of extra
        # CPU<->GPU transfers.
        if self.activation_checkpointing and self.training:
            if self.cpu_offloading:
                with get_cpu_offload_context(optimized=True):
                    for i in range(self.msa_blocks):
                        z, m = checkpoint(self.layers[i], z, m, token_mask, msa_mask, use_reentrant=False)
            else:
                for i in range(self.msa_blocks):
                    z, m = checkpoint(self.layers[i], z, m, token_mask, msa_mask, use_reentrant=False)
        else:
            for i in range(self.msa_blocks):
                z, m = self.layers[i](z, m, token_mask, msa_mask)

        return z


class DistogramModule(nn.Module):
    """Distogram Module using DTensor for Boltz-2.

    This module wraps a serial DistogramModule and adds DTensor-based
    context parallelism support with num_distograms for 5D output.

    The 4D->5D reshape uses shardwise_unflatten to avoid exposing DTensor
    native operations to the autograd graph (which would cause problematic
    all-gather operations during backward).
    """

    def __init__(
        self,
        module: SerialDistogramModule,
        dist_manager: DistributedManager,
        distogram_comm: Optional[TransposeComm] = None,
    ) -> None:
        """Initialize the DTensor-based distogram module.

        Parameters
        ----------
        module : SerialDistogramModule
            Serial DistogramModule from trunkv2 to be distributed.
            Must have num_distograms and num_bins attributes.
        dist_manager : DistributedManager
            Distributed manager for device mesh and process groups.
        distogram_comm : TransposeComm, optional
            Communication object for CP transpose operations.
            Default is None for serial mode.
        """
        super().__init__()
        self.dist_manager = dist_manager
        self.device_mesh = dist_manager.device_mesh_subgroups

        self.distogram = LinearParamsReplicated(module.distogram, device_mesh=self.device_mesh)
        self.distogram_comm = distogram_comm

        self.num_distograms = module.num_distograms
        self.num_bins = module.num_bins

    def forward(self, z: DTensor) -> DTensor:
        """Perform the forward pass.

        Parameters
        ----------
        z : DTensor
            The pairwise embeddings as DTensor with shape [B, N, N, token_z].
            Sharded along dimensions 1 and 2 for CP.

        Returns
        -------
        DTensor
            The predicted distogram with shape [B, N, N, num_distograms, num_bins].
            Maintains the same sharding as input along dimensions 1 and 2.
        """
        x: DTensor = redistribute_transpose(
            z,
            transpose_comm=self.distogram_comm,
            output_placements=(Shard(0), Shard(1), Shard(2)),
            dim0=1,
            dim1=2,
        )
        y: DTensor = elementwise_op(z, x, ElementwiseOp.SUM)

        output_4d: DTensor = self.distogram(y)

        output_5d: DTensor = shardwise_unflatten(
            output_4d,
            dim=3,
            sizes=(self.num_distograms, self.num_bins),
        )

        return output_5d


class BFactorModule(nn.Module):
    """DTensor BFactorModule for Boltz-2.

    Wraps a serial BFactorModule's ``nn.Linear`` with ``LinearParamsReplicated``
    so that parameter gradients are correctly all-reduced across the CP mesh.

    The forward pass is elementwise on the single representation ``s``
    (shape ``[B, N, token_s]``, placements ``(Shard(0), Shard(1), Replicate())``).
    No cross-shard communication is required.
    """

    def __init__(self, module: SerialBFactorModule, device_mesh: DeviceMesh) -> None:
        """Initialize the distributed BFactorModule.

        Parameters
        ----------
        module : SerialBFactorModule
            Serial BFactorModule from trunkv2 to be distributed.
        device_mesh : DeviceMesh
            The device mesh for distributed tensor operations.
        """
        if not isinstance(module, SerialBFactorModule):
            raise TypeError(f"Expected SerialBFactorModule, got {type(module)}")
        super().__init__()
        self.bfactor = LinearParamsReplicated(module.bfactor, device_mesh=device_mesh)
        self.num_bins = module.num_bins

    def forward(self, s: DTensor) -> DTensor:
        """Predict per-token B-factor histogram.

        Parameters
        ----------
        s : DTensor
            Single representation, shape ``[B, N, token_s]``,
            placements ``(Shard(0), Shard(1), Replicate())``.

        Returns
        -------
        DTensor
            Predicted B-factor logits, shape ``[B, N, num_bins]``,
            same placements as input.
        """
        return self.bfactor(s)


class ContactConditioning(nn.Module):
    """DTensor ContactConditioning for Boltz-2.

    Wraps the serial ContactConditioning module for context parallelism.
    All operations are elementwise on the last dimension of pair features
    ``(B, N, N, *)``, so no cross-shard communication is needed.

    The serial FourierEmbedding's ``proj`` is frozen (no grad), so we keep
    it as a plain ``nn.Linear`` and operate on local tensor shards directly
    via ``to_local()`` / ``DTensor.from_local()``.  The trainable ``encoder``
    linear is wrapped with ``LinearParamsReplicated`` for correct gradient
    all-reduce.
    """

    def __init__(self, module: SerialContactConditioning, device_mesh: DeviceMesh) -> None:
        """Initialize the distributed ContactConditioning.

        Parameters
        ----------
        module : SerialContactConditioning
            Serial ContactConditioning from trunkv2.
        device_mesh : DeviceMesh
            The device mesh for distributed tensor operations.
        """
        if not isinstance(module, SerialContactConditioning):
            raise TypeError(f"Expected SerialContactConditioning, got {type(module)}")
        if const.contact_conditioning_info["UNSPECIFIED"] != 0:
            raise ValueError(
                f"Expected UNSPECIFIED index 0, got {const.contact_conditioning_info['UNSPECIFIED']}. "
                "ContactConditioning forward slices cc[:,:,:,0:1] for UNSPECIFIED."
            )
        if const.contact_conditioning_info["UNSELECTED"] != 1:
            raise ValueError(
                f"Expected UNSELECTED index 1, got {const.contact_conditioning_info['UNSELECTED']}. "
                "ContactConditioning forward slices cc[:,:,:,1:2] for UNSELECTED."
            )
        super().__init__()
        self.device_mesh = device_mesh

        if not isinstance(module.fourier_embedding, SerialFourierEmbedding):
            raise TypeError(f"Expected SerialFourierEmbedding, got {type(module.fourier_embedding)}")
        self.fourier_embedding = module.fourier_embedding
        if self.fourier_embedding.proj.weight.requires_grad or (
            self.fourier_embedding.proj.bias is not None and self.fourier_embedding.proj.bias.requires_grad
        ):
            raise ValueError("FourierEmbedding proj should not have trainable parameters")

        self.encoder = LinearParamsReplicated(module.encoder, device_mesh=device_mesh)

        all_replicate = [Replicate()] * device_mesh.ndim
        self.encoding_unspecified = nn.Parameter(
            distribute_tensor(module.encoding_unspecified.data, device_mesh, all_replicate),
            requires_grad=module.encoding_unspecified.requires_grad,
        )
        self.encoding_unselected = nn.Parameter(
            distribute_tensor(module.encoding_unselected.data, device_mesh, all_replicate),
            requires_grad=module.encoding_unselected.requires_grad,
        )

        self.cutoff_min = module.cutoff_min
        self.cutoff_max = module.cutoff_max

    def forward(self, feats: dict[str, DTensor]) -> DTensor:
        """Compute contact conditioning pairwise embeddings.

        Parameters
        ----------
        feats : dict[str, DTensor]
            Must contain:
            - ``contact_conditioning``: shape ``(B, N, N, num_contact_types)``,
              placements ``(Shard(0), Shard(1), Shard(2))``
            - ``contact_threshold``: shape ``(B, N, N)``,
              placements ``(Shard(0), Shard(1), Shard(2))``

        Returns
        -------
        DTensor
            Contact conditioning embeddings, shape ``(B, N, N, token_z)``,
            placements ``(Shard(0), Shard(1), Shard(2))``.
        """
        cc_dt: DTensor = feats["contact_conditioning"]
        ct_dt: DTensor = feats["contact_threshold"]

        cc_local = cc_dt.to_local()
        ct_local = ct_dt.to_local()

        ct_norm = (ct_local - self.cutoff_min) / (self.cutoff_max - self.cutoff_min)
        ct_flat = ct_norm.flatten()
        fourier_flat = torch.cos(2 * pi * self.fourier_embedding.proj(ct_flat.unsqueeze(-1)))
        ct_fourier = fourier_flat.reshape(ct_norm.shape + (-1,))

        cc_features = cc_local[:, :, :, 2:]
        combined = torch.cat(
            [cc_features, ct_norm.unsqueeze(-1), ct_fourier],
            dim=-1,
        )

        combined_shape = cc_dt.shape[:-1] + (combined.shape[-1],)
        combined_contiguous = combined.contiguous()
        combined_stride = update_exhaustive_strides(
            combined_contiguous.shape, combined_contiguous.stride(), combined_shape
        )
        combined_dt = DTensor.from_local(
            combined_contiguous,
            self.device_mesh,
            cc_dt.placements,
            shape=combined_shape,
            stride=combined_stride,
        )
        encoded_dt = self.encoder(combined_dt)

        # Native DTensor arithmetic is used here despite the general rule
        # against implicit DTensor ops on differentiable paths (CLAUDE.md).
        # This is safe because all multiplications are Replicate × Shard(0,1,2)
        # = local elementwise — no hidden all-gathers.  The resulting
        # Partial(Sum) gradients for the encoding params are reduced to
        # Replicate by on_after_backward.
        unspec_flag = cc_dt[:, :, :, 0:1]  # (B, N, N, 1)
        unsel_flag = cc_dt[:, :, :, 1:2]  # (B, N, N, 1)
        mask_factor = 1.0 - (unspec_flag + unsel_flag)  # (B, N, N, 1)

        result = (
            encoded_dt * mask_factor + self.encoding_unspecified * unspec_flag + self.encoding_unselected * unsel_flag
        )
        return result
