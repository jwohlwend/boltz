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

"""DTensor-compatible encoder/decoder modules for Context Parallelism.

Compatible with both Boltz-1x and Boltz-2 serial modules. Focuses on the
window-batching variant.

Modules:
- RelativePositionEncoder: pairwise relative position features → linear projection
- AtomAttentionDecoder: token → atom position updates
- AtomAttentionEncoder / _atom_encoder: atom-level encoder with window batching
"""

import copy
from math import pi

import torch
from torch import nn
from torch.distributed.device_mesh import DeviceMesh
from torch.distributed.tensor import DTensor, Replicate, Shard
from torch.distributed.tensor import zeros as dtensor_zeros
from torch.nn import Module
from torch.nn.functional import one_hot

from boltz.distributed.comm import TransposeComm
from boltz.distributed.model.layers.cat_and_chunk import shardwise_cat
from boltz.distributed.model.layers.elementwise_op import (
    ElementwiseOp,
    elementwise_op,
    scalar_tensor_op,
    single_tensor_op,
)
from boltz.distributed.model.layers.flatten_and_unflatten import (
    shardwise_flatten,
    shardwise_flatten_sharded,
    shardwise_unflatten_sharded,
)
from boltz.distributed.model.layers.gather import distributed_gather
from boltz.distributed.model.layers.layernorm import LayerNormParamsReplicated
from boltz.distributed.model.layers.linear import LinearParamsReplicated
from boltz.distributed.model.layers.outer_gather import distributed_outer_gather
from boltz.distributed.model.layers.redistribute_transpose import redistribute_transpose
from boltz.distributed.model.layers.repeat_interleave import shardwise_repeat_interleave
from boltz.distributed.model.layers.replicate_op import ReplicateOp, replicate_op
from boltz.distributed.model.layers.scatter import distributed_scatter_reduce
from boltz.distributed.model.layers.shardwise_op import ShardwiseOuterOp, shardwise_outer_op, shardwise_sum
from boltz.distributed.model.layers.squeeze import shardwise_unsqueeze
from boltz.distributed.model.layers.transition import Transition
from boltz.distributed.model.layers.utils import convert_single_repr_to_window_batched_key
from boltz.distributed.model.modules.transformers import AtomTransformer
from boltz.distributed.model.modules.utils import validate_window_batching_parameters
from boltz.model.modules.encoders import AtomAttentionDecoder as AtomAttentionDecoderBoltz1
from boltz.model.modules.encoders import AtomAttentionEncoder as AtomAttentionEncoderBoltz1
from boltz.model.modules.encoders import FourierEmbedding as FourierEmbeddingBoltz1
from boltz.model.modules.encoders import PairwiseConditioning as PairwiseConditioningBoltz1
from boltz.model.modules.encoders import RelativePositionEncoder as SerialRelativePositionEncoderV1
from boltz.model.modules.encoders import SingleConditioning as SingleConditioningBoltz1
from boltz.model.modules.encodersv2 import AtomAttentionDecoder as AtomAttentionDecoderBoltz2
from boltz.model.modules.encodersv2 import AtomAttentionEncoder as AtomAttentionEncoderBoltz2
from boltz.model.modules.encodersv2 import AtomEncoder as AtomEncoderBoltz2
from boltz.model.modules.encodersv2 import FourierEmbedding as FourierEmbeddingBoltz2
from boltz.model.modules.encodersv2 import PairwiseConditioning as PairwiseConditioningBoltz2
from boltz.model.modules.encodersv2 import RelativePositionEncoder as SerialRelativePositionEncoderV2
from boltz.model.modules.encodersv2 import SingleConditioning as SingleConditioningBoltz2


class RelativePositionEncoder(Module):
    """DTensor RelativePositionEncoder for Boltz-1x and Boltz-2.

    Computes pairwise relative position features from single-representation
    features (``asym_id``, ``residue_index``, ``entity_id``, ``token_index``,
    ``sym_id``, ``cyclic_period``) and projects them through a linear layer.

    Under context parallelism the single features are sharded along the token
    dimension.  Pairwise outer comparisons (``feat_i[:, :, None] op feat_j[:, None, :]``)
    require the "column" shard from the transposed rank, obtained via
    ``redistribute_transpose``.

    All intermediate computation (outer ops, clipping, one-hot encoding) is
    non-differentiable and operates on local tensors.  Only the final linear
    projection (``LinearParamsReplicated``) is differentiable.

    Communication budget (forward):
        - 6 ``redistribute_transpose`` calls (one per feature key) for the
          column shards.  Each is a single P2P send/recv.
        - 0 additional collectives (the ``LinearParamsReplicated`` backward
          handles gradient all-reduce).

    Supports both Boltz-1x and Boltz-2 serial modules:
        - Boltz-1x: ``RelativePositionEncoder(token_z, r_max, s_max)``
        - Boltz-2 adds ``fix_sym_check`` and ``cyclic_pos_enc`` flags.
    """

    # Feature keys that need column-shard transpose
    _KEYS_TO_TRANSPOSE = ("asym_id", "entity_id", "residue_index", "token_index", "sym_id", "cyclic_period")

    def __init__(
        self,
        layer: Module,
        device_mesh: DeviceMesh,
        transpose_comm: TransposeComm,
    ) -> None:
        """Initialize the distributed RelativePositionEncoder.

        Parameters
        ----------
        layer : Module
            Serial RelativePositionEncoder (v1 or v2).
        device_mesh : DeviceMesh
            The device mesh (subgroups mesh with dp, cp_axis_0, cp_axis_1).
        transpose_comm : TransposeComm
            Transpose communication for distributed outer operations.
            Separate deep copies are created for each feature key.
        """
        if not isinstance(layer, (SerialRelativePositionEncoderV1, SerialRelativePositionEncoderV2)):
            raise TypeError(f"Expected SerialRelativePositionEncoderV1 or V2, got {type(layer)}")
        super().__init__()
        self.device_mesh = device_mesh
        self.r_max = layer.r_max
        self.s_max = layer.s_max

        # V2-only flags (default to V1 behaviour).
        # V1 has no cyclic_pos_enc attr but always applies the cyclic correction,
        # so default to True to unify V1 and V2-with-flag-on into the same branch.
        self.fix_sym_check = getattr(layer, "fix_sym_check", False)
        self.cyclic_pos_enc = getattr(layer, "cyclic_pos_enc", True)

        # Wrap the linear layer
        self.linear_layer = LinearParamsReplicated(layer.linear_layer, device_mesh=device_mesh)

        # One TransposeComm per feature key (separate P2P buffers).
        # TransposeComm is not an nn.Module, so store as plain dict + individual attrs.
        self._transpose_comms: dict[str, TransposeComm] = {}
        for i, key in enumerate(self._KEYS_TO_TRANSPOSE):
            tc = transpose_comm if i == 0 else copy.deepcopy(transpose_comm)
            self._transpose_comms[key] = tc
            setattr(self, f"_tc_{key}", tc)  # for pickling visibility

    def forward(self, feats: dict[str, DTensor]) -> DTensor:
        """Compute relative position embeddings.

        Parameters
        ----------
        feats : dict[str, DTensor]
            Must contain keys: ``asym_id``, ``entity_id``, ``residue_index``,
            ``token_index``, ``sym_id``, ``cyclic_period``.
            Each has shape ``(B, N)`` with placements
            ``(Shard(0), Shard(1), Replicate())``.

        Returns
        -------
        DTensor
            Relative position embeddings, shape ``(B, N, N, token_z)``,
            placements ``(Shard(0), Shard(1), Shard(2))``.
        """
        expected_placements = (Shard(0), Shard(1))
        for key in self._KEYS_TO_TRANSPOSE:
            dt = feats[key]
            if not isinstance(dt, DTensor):
                raise TypeError(f"Expected DTensor for '{key}', got {type(dt)}")
            # Check first two placements (third may be Replicate for 3D mesh)
            if dt.placements[:2] != expected_placements:
                raise ValueError(
                    f"Expected '{key}' placements to start with {expected_placements}, got {dt.placements}"
                )

        # Get column shards via redistribute_transpose.
        # Row feats: (B, N) placements (Shard(0), Shard(1), Replicate())
        # Column feats: swap cp_axis_0 ↔ cp_axis_1 → (Shard(0), Replicate(), Shard(1))
        # so each rank gets the transpose peer's token shard for outer ops.
        col_placements = (Shard(0), Replicate(), Shard(1))

        feats_col = {}
        for key in self._KEYS_TO_TRANSPOSE:
            feats_col[key] = redistribute_transpose(
                feats[key],
                transpose_comm=self._transpose_comms[key],
                output_placements=col_placements,
                dim0=None,
                dim1=None,
            )

        # Extract local tensors for non-differentiable feature computation
        row = {k: feats[k].to_local() for k in self._KEYS_TO_TRANSPOSE}
        col = {k: feats_col[k].to_local() for k in self._KEYS_TO_TRANSPOSE}

        # Pairwise comparisons: row[:, :, None] op col[:, None, :]
        b_same_chain = torch.eq(row["asym_id"][:, :, None], col["asym_id"][:, None, :])
        b_same_residue = torch.eq(row["residue_index"][:, :, None], col["residue_index"][:, None, :])
        b_same_entity = torch.eq(row["entity_id"][:, :, None], col["entity_id"][:, None, :])

        d_residue = row["residue_index"][:, :, None] - col["residue_index"][:, None, :]

        # Cyclic period adjustment.
        # The serial code guards with torch.any(feats["cyclic_period"] > 0)
        # over the full tensor, but in CP each rank only sees its local row
        # and column shards.  Rather than broadcasting a flag across ranks,
        # we unconditionally apply the correction — the torch.where with
        # fallback period=10000 makes it a no-op when no token has a
        # positive cyclic period (round(d/10000) == 0 for typical d).
        if self.cyclic_pos_enc:
            period = torch.where(
                col["cyclic_period"] > 0,
                col["cyclic_period"],
                torch.zeros_like(col["cyclic_period"]) + 10000,
            )
            d_residue = (d_residue - period[:, None, :] * torch.round(d_residue / period[:, None, :])).long()
        # cyclic_pos_enc=False (V2 only): skip cyclic correction entirely

        d_residue = torch.clip(d_residue + self.r_max, 0, 2 * self.r_max)
        d_residue = torch.where(b_same_chain, d_residue, torch.zeros_like(d_residue) + 2 * self.r_max + 1)
        a_rel_pos = one_hot(d_residue, 2 * self.r_max + 2)

        d_token = torch.clip(
            row["token_index"][:, :, None] - col["token_index"][:, None, :] + self.r_max,
            0,
            2 * self.r_max,
        )
        d_token = torch.where(
            b_same_chain & b_same_residue,
            d_token,
            torch.zeros_like(d_token) + 2 * self.r_max + 1,
        )
        a_rel_token = one_hot(d_token, 2 * self.r_max + 2)

        d_chain = torch.clip(
            row["sym_id"][:, :, None] - col["sym_id"][:, None, :] + self.s_max,
            0,
            2 * self.s_max,
        )
        if self.fix_sym_check:
            # V2 path: sentinel when NOT same entity
            d_chain = torch.where(
                ~b_same_entity,
                torch.zeros_like(d_chain) + 2 * self.s_max + 1,
                d_chain,
            )
        else:
            # V1 path: sentinel when same chain
            d_chain = torch.where(
                b_same_chain,
                torch.zeros_like(d_chain) + 2 * self.s_max + 1,
                d_chain,
            )
        a_rel_chain = one_hot(d_chain, 2 * self.s_max + 2)

        # Concatenate and cast to linear weight dtype
        dtype = self.linear_layer.weight.to_local().dtype
        features_local = torch.cat(
            [
                a_rel_pos.to(dtype),
                a_rel_token.to(dtype),
                b_same_entity.unsqueeze(-1).to(dtype),
                a_rel_chain.to(dtype),
            ],
            dim=-1,
        )

        # Wrap as DTensor with pair placements for the linear layer
        pair_placements = (Shard(0), Shard(1), Shard(2))
        # Pad placements to match mesh ndim (e.g., 3D mesh: dp, cp0, cp1)
        if self.device_mesh.ndim > 3:
            raise ValueError(f"Expected device mesh ndim <= 3, got {self.device_mesh.ndim}")
        full_pair_placements = pair_placements[: self.device_mesh.ndim]

        # Compute global shape and contiguous strides for DTensor.from_local
        asym_dt = feats["asym_id"]
        B_global = asym_dt.shape[0]
        N_global = asym_dt.shape[1]
        feat_dim = features_local.shape[-1]
        global_shape = (B_global, N_global, N_global, feat_dim)
        global_stride = (
            N_global * N_global * feat_dim,
            N_global * feat_dim,
            feat_dim,
            1,
        )

        features_dt = DTensor.from_local(
            features_local.contiguous(),
            self.device_mesh,
            full_pair_placements,
            shape=global_shape,
            stride=global_stride,
        )

        return self.linear_layer(features_dt)


class AtomAttentionDecoder(Module):
    """DTensor AtomAttentionDecoder for window batching.

    Compatible with both Boltz-1x and Boltz-2 serial AtomAttentionDecoder modules.

    This module converts token representations to atom-level position updates
    using an AtomTransformer. The window batching variant uses distributed_gather
    with global atom-to-token indices instead of torch.bmm with one-hot matrices.

    Key operations:
    1. Transform token repr to atom space via a_to_q_trans (Linear)
    2. Gather from token to atom positions via distributed_gather
    3. Run AtomTransformer (window-batched)
    4. Project atom features to position updates via atom_feat_to_atom_pos_update
    """

    def __init__(
        self,
        layer: nn.Module,
        device_mesh: DeviceMesh,
    ):
        """Initialize the DTensor-distributed atom attention decoder.

        Parameters
        ----------
        layer : nn.Module
            The serial AtomAttentionDecoder module to be distributed.
            Accepts both Boltz-1x and Boltz-2 versions.
        device_mesh : DeviceMesh
            The device mesh for distributed tensor operations.

        Raises
        ------
        TypeError
            If layer is not a recognized type.
        """
        super().__init__()

        if not isinstance(layer, (AtomAttentionDecoderBoltz1, AtomAttentionDecoderBoltz2)):
            raise TypeError(
                ", ".join(
                    [
                        f"Instance {layer} should have type "
                        f"{AtomAttentionDecoderBoltz1} or {AtomAttentionDecoderBoltz2}",
                        f"but instead has type {type(layer)}.",
                    ]
                )
            )

        self.could_use_model_cache = isinstance(layer, AtomAttentionDecoderBoltz1)
        self.attn_window_queries = layer.atom_decoder.attn_window_queries
        self.attn_window_keys = layer.atom_decoder.attn_window_keys
        validate_window_batching_parameters(self.attn_window_queries, self.attn_window_keys, use_window_batching=True)

        # a_to_q_trans: LinearNoBias(2 * token_s, atom_s)
        self.a_to_q_trans = LinearParamsReplicated(layer_local=layer.a_to_q_trans, device_mesh=device_mesh)

        # atom_decoder: DTensor AtomTransformer (window batching)
        self.atom_decoder = AtomTransformer(layer=layer.atom_decoder, device_mesh=device_mesh)

        # atom_feat_to_atom_pos_update:
        # Boltz-1: always Sequential(LayerNorm, LinearNoBias) -- no post_layer_norm support
        # Boltz-2 with transformer_post_layer_norm=False: Sequential(LayerNorm, LinearNoBias)
        # Boltz-2 with transformer_post_layer_norm=True: just LinearNoBias (no LayerNorm)
        #
        # Infer transformer_post_layer_norm from the inner DiffusionTransformerLayer's
        # post_lnorm attribute: nn.LayerNorm means True, nn.Identity means False.
        # Boltz-1 DiffusionTransformerLayer has no post_lnorm attribute (always False).
        first_dtl = layer.atom_decoder.diffusion_transformer.layers[0]
        transformer_post_layer_norm = hasattr(first_dtl, "post_lnorm") and not isinstance(
            first_dtl.post_lnorm, nn.Identity
        )

        if transformer_post_layer_norm:
            # Boltz-2 with transformer_post_layer_norm=True: just LinearNoBias
            self.atom_feat_to_atom_pos_update = LinearParamsReplicated(
                layer_local=layer.atom_feat_to_atom_pos_update, device_mesh=device_mesh
            )
        else:
            # Boltz-1, or Boltz-2 with transformer_post_layer_norm=False:
            # Sequential(LayerNorm, LinearNoBias)
            self.atom_feat_to_atom_pos_update = nn.Sequential(
                LayerNormParamsReplicated(layer.atom_feat_to_atom_pos_update[0], device_mesh=device_mesh),
                LinearParamsReplicated(layer_local=layer.atom_feat_to_atom_pos_update[1], device_mesh=device_mesh),
            )

    def forward(
        self,
        a: DTensor,
        q: DTensor,
        c: DTensor,
        p: DTensor,
        feats: dict[str, DTensor],
        multiplicity: int = 1,
        model_cache: dict[str, dict[str, DTensor]] | None = None,
    ) -> DTensor:
        """Forward pass for the DTensor-distributed atom attention decoder.

        All tensors use device mesh (dp, cp_axis_0, cp_axis_1).
        Placements: Shard(0)=dp batch, Shard(1)=cp atom/token axis, Replicate()=cp_axis_1.

        Parameters
        ----------
        a : DTensor
            Token representation, shape (B * M, N_tokens, 2 * token_s).
            Placements: (Shard(0), Shard(1), Replicate()).
        q : DTensor
            Atom query representation, shape (B * M, N_atoms_packed, atom_s).
            Placements: (Shard(0), Shard(1), Replicate()).
        c : DTensor
            Atom conditioning representation, shape (B * M, N_atoms_packed, atom_s).
            Placements: (Shard(0), Shard(1), Replicate()).
        p : DTensor
            Pair representation / pre-computed bias in window-batched format.
            - Boltz-1: shape (B, K, W, H, c_z)
            - Boltz-2: shape (B, K, W, H, num_heads * depth)
            Placements: (Shard(0), Shard(1), Replicate()).
        feats : dict[str, DTensor]
            Features dict containing:
            - "atom_pad_mask": (B, N_atoms_packed), placements per atom_features config
            - "atom_to_token_ids_global": (B, N_atoms_packed), placements per atom_features config
        multiplicity : int, optional
            Number of diffusion samples, by default 1.
        model_cache : dict or None, optional
            Model cache for inference optimization (V1 internalized path only).

        Returns
        -------
        DTensor
            Position updates, shape (B * M, N_atoms_packed, 3).
            Placements: (Shard(0), Shard(1), Replicate()).
        """
        if model_cache is not None and not self.could_use_model_cache:
            raise ValueError("model_cache is only supported with V1 AtomAttentionDecoder")

        W = self.attn_window_queries
        N = q.shape[1]  # N_atoms_packed

        if N % W != 0:
            raise ValueError(
                f"Packed atom sequence length N={N} must be divisible by window size W={W} "
                f"for window batching, but N % W = {N % W}"
            )

        K = N // W

        # Get atom mask (without multiplicity -- matches p's batch dim)
        atom_mask = feats["atom_pad_mask"].bool()

        # Get global atom-to-token indices
        atom_to_token_ids_global = feats["atom_to_token_ids_global"]

        # Apply multiplicity to indices and mask for gather (must match a's batch dim B*M)
        atom_mask_mul = shardwise_repeat_interleave(atom_mask, multiplicity, 0)
        atom_to_token_ids_global = shardwise_repeat_interleave(atom_to_token_ids_global, multiplicity, 0)

        # Unflatten to window view for gather: (B*M, N) -> (B*M, K, W)
        atom_mask_q = shardwise_unflatten_sharded(atom_mask_mul, axis=1, sizes=(K, W))
        atom_to_token_ids_global_q = shardwise_unflatten_sharded(atom_to_token_ids_global, axis=1, sizes=(K, W))

        # a_to_q transform and gather with autocast disabled (matching Boltz-2 serial behavior
        # at encodersv2.py:544 which uses torch.autocast("cuda", enabled=False) for numerical
        # precision in the atom-to-token gather operation)
        with torch.autocast("cuda", enabled=False):
            # (B*M, N_tokens, 2*token_s) -> (B*M, N_tokens, atom_s)
            a_to_q = self.a_to_q_trans(a)

            # Gather from token to atom: (B*M, N_tokens, atom_s) -> (B*M, K, W, atom_s)
            # Equivalent to torch.bmm(atom_to_token, a_to_q) in serial code
            a_to_q = distributed_gather(
                a_to_q, atom_to_token_ids_global_q, axis=1, are_ids_contiguous=True, idx_mask=atom_mask_q
            )

        # Flatten back to atom sequence: (B*M, K, W, atom_s) -> (B*M, N, atom_s)
        a_to_q = shardwise_flatten_sharded(a_to_q, start_dim=1, end_dim=2)

        # Add to q
        q = elementwise_op(q, a_to_q, ElementwiseOp.SUM)

        # V1 model_cache for decoder transformer
        layer_cache = None
        if model_cache is not None:
            cache_prefix = "atomdecoder"
            if cache_prefix not in model_cache:
                model_cache[cache_prefix] = {}
            layer_cache = model_cache[cache_prefix]

        # Call AtomTransformer with window batching
        # multiplicity=1 because multiplicity is already applied to q, c
        # mask should NOT have multiplicity applied (to match p's batch dim)
        q = self.atom_decoder(
            q=q,
            c=c,
            p=p,
            mask=atom_mask,  # NO multiplicity - matches p's batch dim
            multiplicity=1,  # multiplicity already applied to q, c
            model_cache=layer_cache,
            pair_mask=None,  # window batching doesn't use pair_mask
        )

        r_update = self.atom_feat_to_atom_pos_update(q)
        return r_update


class _MaskPaddingAtoms(torch.autograd.Function):
    """Zero out atom representations at padding positions.

    Intersperse padding creates atoms with zero features; ``embed_atom_features``
    maps them to the layer bias, leaking non-zero values into pair features
    (``c_to_p_trans_q/k``) and downstream attention.  This function multiplies
    ``c`` (or any atom-level DTensor with ``(Shard(0), Shard(1), Replicate())``)
    by the binary ``atom_pad_mask``, broadcasting along the feature dim.

    The mask is non-differentiable — backward simply re-applies the mask to
    the upstream gradient.  No collectives are issued.

    Parameters
    ----------
    c : DTensor
        Atom representations, shape ``(B, N_atoms, D)``.
    atom_pad_mask : DTensor
        Binary mask, shape ``(B, N_atoms)``.  1 = real atom, 0 = padding.
    """

    @staticmethod
    @torch.amp.custom_fwd(device_type="cuda")
    def forward(ctx, c: DTensor, atom_pad_mask: DTensor) -> DTensor:
        if not isinstance(c, DTensor):
            raise TypeError(f"c must be a DTensor, got {type(c)}")
        if not isinstance(atom_pad_mask, DTensor):
            raise TypeError(f"atom_pad_mask must be a DTensor, got {type(atom_pad_mask)}")
        if c.device_mesh != atom_pad_mask.device_mesh:
            raise ValueError("c and atom_pad_mask must share the same device_mesh")

        c_local = c.to_local()
        mask_local = atom_pad_mask.to_local().to(c_local.dtype).unsqueeze(-1)
        result_local = c_local * mask_local

        ctx.save_for_backward(mask_local)
        ctx.placements = list(c.placements)
        ctx.device_mesh = c.device_mesh
        ctx.shape_c = c.shape
        ctx.stride_c = c.stride()

        return DTensor.from_local(
            result_local,
            c.device_mesh,
            list(c.placements),
            shape=c.shape,
            stride=c.stride(),
        )

    @staticmethod
    @torch.amp.custom_bwd(device_type="cuda")
    def backward(ctx, grad_output: DTensor):
        (mask_local,) = ctx.saved_tensors
        grad_local = grad_output.to_local() * mask_local
        grad_c = DTensor.from_local(
            grad_local,
            ctx.device_mesh,
            ctx.placements,
            shape=ctx.shape_c,
            stride=ctx.stride_c,
        )
        return grad_c, None


def _mask_padding_atoms(c: DTensor, atom_pad_mask: DTensor) -> DTensor:
    """Zero out representations at padding atom positions.  See :class:`_MaskPaddingAtoms`."""
    return _MaskPaddingAtoms.apply(c, atom_pad_mask)


def _atom_encoder(
    c: DTensor,
    embed_atompair_ref_pos: nn.Module,
    embed_atompair_ref_dist: nn.Module,
    embed_atompair_mask: nn.Module,
    s_to_c_trans: nn.Module | None,
    z_to_p_trans: nn.Module | None,
    c_to_p_trans_q: nn.Module,
    c_to_p_trans_k: nn.Module,
    p_mlp: nn.Module,
    feats: dict[str, DTensor],
    s_trunk: DTensor | None,
    z: DTensor | None,
    structure_prediction: bool,
    W: int,
    H: int,
) -> tuple[DTensor, DTensor, DTensor]:
    """Shared DTensor pair computation for atom encoding (window batching).

    Migrated from V1x DTensor src_v1/boltz/distributed/model/modules/encoders.py
    lines 776-893. Takes pre-embedded c and performs window-batched pair computation.

    This function does NOT manage autocast context. The caller is responsible for
    wrapping the call in torch.autocast("cuda", enabled=False) when used from the
    V2 AtomEncoder path (matching V2 serial encodersv2.py:297). The V1 path
    (AtomAttentionEncoder) should NOT disable autocast, matching V1 serial behavior.

    Parameters
    ----------
    c : DTensor
        Pre-embedded atom single representation, shape (B, N_atoms, atom_s).
        Placements: (Shard(0), Shard(1), Replicate()).
    embed_atompair_ref_pos : nn.Module
        Linear for pair position embedding.
    embed_atompair_ref_dist : nn.Module
        Linear for pair distance embedding.
    embed_atompair_mask : nn.Module
        Linear for pair mask embedding.
    s_to_c_trans : nn.Module or None
        Sequential(LayerNorm, Linear) for token-to-atom conditioning. None if not structure_prediction.
    z_to_p_trans : nn.Module or None
        Sequential(LayerNorm, Linear) for pair token-to-atom conditioning. None if not structure_prediction.
    c_to_p_trans_q : nn.Module
        Sequential(ReLU, Linear) for query pair contribution.
    c_to_p_trans_k : nn.Module
        Sequential(ReLU, Linear) for key pair contribution.
    p_mlp : nn.Module
        Sequential MLP for pair representation refinement.
    feats : dict[str, DTensor]
        Must contain: atom_pad_mask, ref_pos, ref_space_uid, atom_to_token_ids_global.
    s_trunk : DTensor or None
        Token single representation (B, N_tokens, token_s).
        Placements: (Shard(0), Shard(1), Replicate()). None if not structure_prediction.
    z : DTensor or None
        Token pair representation (B, N_tokens, N_tokens, token_z).
        Placements: (Shard(0), Shard(1), Shard(2)). None if not structure_prediction.
    structure_prediction : bool
        Whether to apply token-to-atom conditioning.
    W : int
        Atoms per window for queries.
    H : int
        Atoms per window for keys.

    Returns
    -------
    tuple[DTensor, DTensor, DTensor]
        (q, c, p):
        - q: initial c before conditioning, shape (B, N_atoms, atom_s).
          Placements: (Shard(0), Shard(1), Replicate()).
        - c: s_to_c-conditioned atom representation, shape (B, N_atoms, atom_s).
          Placements: (Shard(0), Shard(1), Replicate()).
        - p: atom pair representation, shape (B, K, W, H, atom_z) where K = N_atoms // W.
          Placements: (Shard(0), Shard(1), Replicate()).
    """
    # Sanity checks: structure_prediction-dependent arguments must be consistently None or non-None.
    # When structure_prediction=True, s_to_c_trans, z_to_p_trans, s_trunk, z are all required.
    # When structure_prediction=False, they must all be None.
    if structure_prediction:
        if s_to_c_trans is None or z_to_p_trans is None:
            raise ValueError("structure_prediction=True requires s_to_c_trans and z_to_p_trans to be provided.")
        if s_trunk is None or z is None:
            raise ValueError("structure_prediction=True requires s_trunk and z to be provided.")
    else:
        if s_to_c_trans is not None or z_to_p_trans is not None:
            raise ValueError("structure_prediction=False but s_to_c_trans or z_to_p_trans was provided.")
        if s_trunk is not None or z is not None:
            raise ValueError("structure_prediction=False but s_trunk or z was provided.")

    N = c.shape[1]
    if N % W != 0:
        raise ValueError(
            f"Sequence length N={N} must be divisible by window size W={W} for window batching, but N % W = {N % W}"
        )
    K = N // W

    # Mimic serial code's .float() casts: promote to at least float32 for
    # numerical stability, but preserve higher precision (e.g. float64) if available.
    compute_dtype = torch.promote_types(c.dtype, torch.float32)

    atom_ref_pos = feats["ref_pos"]
    atom_mask_bool = feats["atom_pad_mask"].bool()
    atom_uid = feats["ref_space_uid"]

    # Convert atom_ref_pos to query/key views
    atom_ref_pos_q = shardwise_unflatten_sharded(atom_ref_pos, axis=1, sizes=(K, W))
    atom_ref_pos_k = convert_single_repr_to_window_batched_key(atom_ref_pos, W, H)

    # Compute distance: d = keys - queries
    # shardwise_outer_op computes lhs - rhs = queries - keys, so negate
    d = shardwise_outer_op(atom_ref_pos_q, atom_ref_pos_k, axis=2, op=ShardwiseOuterOp.SUBTRACT)
    d = scalar_tensor_op(-1.0, d, ElementwiseOp.PROD)
    d_norm = shardwise_sum(elementwise_op(d, d, ElementwiseOp.PROD), dim=-1, keepdim=True)
    d_norm = scalar_tensor_op(1.0, scalar_tensor_op(1.0, d_norm, ElementwiseOp.SUM), ElementwiseOp.DIV)

    # Compute validity mask
    atom_mask_q = shardwise_unflatten_sharded(atom_mask_bool, axis=1, sizes=(K, W))
    atom_mask_k = convert_single_repr_to_window_batched_key(atom_mask_bool, W, H)
    atom_uid_q = shardwise_unflatten_sharded(atom_uid, axis=1, sizes=(K, W))
    atom_uid_k = convert_single_repr_to_window_batched_key(atom_uid, W, H)

    # v = (atom_mask_q & atom_mask_k & (atom_uid_q == atom_uid_k))
    mask_and = shardwise_outer_op(atom_mask_q, atom_mask_k, axis=2, op=ShardwiseOuterOp.LOGICAL_AND)
    uid_eq = shardwise_outer_op(atom_uid_q, atom_uid_k, axis=2, op=ShardwiseOuterOp.EQUAL)
    v = elementwise_op(mask_and, uid_eq, ElementwiseOp.BITAND)
    # Serial: (...).float().unsqueeze(-1) -- use compute_dtype instead of .float()
    v = shardwise_unsqueeze(v, -1).to(compute_dtype)

    # Compute pair representation p: (B, K, W, H, atom_z)
    # TODO: this DTensor native broadcasting elementwise multiplication should be
    # replaced with custom autograd.Function
    p = embed_atompair_ref_pos(d) * v
    p = elementwise_op(p, embed_atompair_ref_dist(d_norm) * v, ElementwiseOp.SUM)
    p = elementwise_op(p, embed_atompair_mask(v) * v, ElementwiseOp.SUM)

    q = c

    if structure_prediction:
        atom_to_token_ids_global = feats["atom_to_token_ids_global"]
        atom_to_token_ids_global_q = shardwise_unflatten_sharded(atom_to_token_ids_global, axis=1, sizes=(K, W))

        # Token-to-atom gather: s_trunk -> s_to_c -> gather to atom positions
        # Serial: s_to_c_trans(s_trunk.float()) -- use compute_dtype
        s_to_c = s_to_c_trans(s_trunk.to(compute_dtype))
        s_to_c = distributed_gather(
            s_to_c, atom_to_token_ids_global_q, axis=1, are_ids_contiguous=True, idx_mask=atom_mask_q
        )
        s_to_c = shardwise_flatten_sharded(s_to_c, start_dim=1, end_dim=2)
        # Serial: c = c + s_to_c.to(c) -- cast back to c's dtype
        c = elementwise_op(c, s_to_c.to(c.dtype), ElementwiseOp.SUM)

        # Pair token-to-atom gather: z -> z_to_p -> outer gather
        # Serial: z_to_p_trans(z.float()) -- use compute_dtype
        z_to_p = z_to_p_trans(z.to(compute_dtype))
        atom_to_token_ids_global_k = convert_single_repr_to_window_batched_key(atom_to_token_ids_global, W, H)
        z_to_p = distributed_outer_gather(
            z_to_p,
            atom_to_token_ids_global_q,
            atom_to_token_ids_global_k,
            axis=1,
            are_ids_contiguous=True,
            idx_n_mask=atom_mask_q,
            idx_m_mask=atom_mask_k,
        )
        # Serial: p = p + z_to_p.to(p) -- cast back to p's dtype
        p = elementwise_op(p, z_to_p.to(p.dtype), ElementwiseOp.SUM)

    # c_to_p contributions in window-batched form
    c_q = shardwise_unflatten_sharded(c, axis=1, sizes=(K, W))
    c_q = c_to_p_trans_q(c_q)
    c_k = convert_single_repr_to_window_batched_key(c, W, H)
    c_k = c_to_p_trans_k(c_k)
    c_qk = shardwise_outer_op(c_q, c_k, axis=2, op=ShardwiseOuterOp.ADD)

    p = elementwise_op(p, c_qk, ElementwiseOp.SUM)
    p = elementwise_op(p, p_mlp(p), ElementwiseOp.SUM)

    return q, c, p


class AtomEncoder(Module):
    """DTensor AtomEncoder for V2 (window batching).

    Wraps V2 serial AtomEncoder. Flat attributes match V2 serial checkpoint keys.
    Performs V2-specific feature construction then calls _atom_encoder() for
    shared pair computation.
    """

    def __init__(
        self,
        layer: AtomEncoderBoltz2,
        device_mesh: DeviceMesh,
    ):
        """Initialize the DTensor-distributed atom encoder.

        Parameters
        ----------
        layer : AtomEncoderBoltz2
            The serial AtomEncoder module (V2) to be distributed.
        device_mesh : DeviceMesh
            The device mesh for distributed tensor operations.
        """
        super().__init__()

        if not isinstance(layer, AtomEncoderBoltz2):
            raise TypeError(
                ", ".join(
                    [
                        f"Instance {layer} should have type {AtomEncoderBoltz2}",
                        f"but instead has type {type(layer)}.",
                    ]
                )
            )

        # V2-specific feature construction config.
        # use_residue_feats_atoms and use_atom_backbone_feat default to False and are never
        # set to True in the Boltz-2 training/inference workflow (structurev2.yaml, boltz2.py).
        # Fail early at init time rather than at forward time.
        if layer.use_residue_feats_atoms:
            raise NotImplementedError(
                "DTensor AtomEncoder does not yet support use_residue_feats_atoms=True. "
                "This requires DTensor-ified one_hot and token-to-atom gather for residue features. "
                "No Boltz-2 workflow currently enables this option."
            )
        if layer.use_atom_backbone_feat:
            raise NotImplementedError(
                "DTensor AtomEncoder does not yet support use_atom_backbone_feat=True. "
                "No Boltz-2 workflow currently enables this option."
            )
        self.use_no_atom_char = layer.use_no_atom_char
        self.use_atom_backbone_feat = layer.use_atom_backbone_feat
        self.use_residue_feats_atoms = layer.use_residue_feats_atoms
        self.structure_prediction = layer.structure_prediction
        self.atoms_per_window_queries = layer.atoms_per_window_queries
        self.atoms_per_window_keys = layer.atoms_per_window_keys
        validate_window_batching_parameters(
            self.atoms_per_window_queries, self.atoms_per_window_keys, use_window_batching=True
        )

        # Atom feature embedding (V2 uses Linear with bias)
        self.embed_atom_features = LinearParamsReplicated(
            layer_local=layer.embed_atom_features, device_mesh=device_mesh
        )

        # Pair computation modules
        self.embed_atompair_ref_pos = LinearParamsReplicated(
            layer_local=layer.embed_atompair_ref_pos, device_mesh=device_mesh
        )
        self.embed_atompair_ref_dist = LinearParamsReplicated(
            layer_local=layer.embed_atompair_ref_dist, device_mesh=device_mesh
        )
        self.embed_atompair_mask = LinearParamsReplicated(
            layer_local=layer.embed_atompair_mask, device_mesh=device_mesh
        )

        if self.structure_prediction:
            self.s_to_c_trans = nn.Sequential(
                LayerNormParamsReplicated(layer.s_to_c_trans[0], device_mesh=device_mesh),
                LinearParamsReplicated(layer_local=layer.s_to_c_trans[1], device_mesh=device_mesh),
            )
            self.z_to_p_trans = nn.Sequential(
                LayerNormParamsReplicated(layer.z_to_p_trans[0], device_mesh=device_mesh),
                LinearParamsReplicated(layer_local=layer.z_to_p_trans[1], device_mesh=device_mesh),
            )

        self.c_to_p_trans_q = nn.Sequential(
            nn.ReLU(),
            LinearParamsReplicated(layer_local=layer.c_to_p_trans_q[1], device_mesh=device_mesh),
        )
        self.c_to_p_trans_k = nn.Sequential(
            nn.ReLU(),
            LinearParamsReplicated(layer_local=layer.c_to_p_trans_k[1], device_mesh=device_mesh),
        )
        self.p_mlp = nn.Sequential(
            nn.ReLU(),
            LinearParamsReplicated(layer_local=layer.p_mlp[1], device_mesh=device_mesh),
            nn.ReLU(),
            LinearParamsReplicated(layer_local=layer.p_mlp[3], device_mesh=device_mesh),
            nn.ReLU(),
            LinearParamsReplicated(layer_local=layer.p_mlp[5], device_mesh=device_mesh),
        )

    def forward(
        self,
        feats: dict[str, DTensor],
        s_trunk: DTensor | None = None,
        z: DTensor | None = None,
    ) -> tuple[DTensor, DTensor, DTensor]:
        """Forward pass for the DTensor-distributed atom encoder.

        Parameters
        ----------
        feats : dict[str, DTensor]
            Atom features including ref_pos, ref_charge, ref_element, atom_pad_mask,
            ref_space_uid, atom_to_token_ids_global, and optionally ref_atom_name_chars.
            Atom-level tensors have placements (Shard(0), Shard(1), Replicate()).
        s_trunk : DTensor or None
            Token single representation (B, N_tokens, token_s).
            Placements: (Shard(0), Shard(1), Replicate()). None if not structure_prediction.
        z : DTensor or None
            Token pair representation (B, N_tokens, N_tokens, token_z).
            Placements: (Shard(0), Shard(1), Shard(2)). None if not structure_prediction.

        Returns
        -------
        tuple[DTensor, DTensor, DTensor]
            (q, c, p):
            - q: initial atom representation, shape (B, N_atoms, atom_s).
              Placements: (Shard(0), Shard(1), Replicate()).
            - c: conditioned atom representation, shape (B, N_atoms, atom_s).
              Placements: (Shard(0), Shard(1), Replicate()).
            - p: atom pair representation, shape (B, K, W, H, atom_z) where K = N_atoms // W.
              Placements: (Shard(0), Shard(1), Replicate()).
        """
        # V2 feature construction and pair computation run with autocast disabled,
        # matching V2 serial AtomEncoder.forward() (encodersv2.py) which wraps
        # the entire forward in torch.autocast("cuda", enabled=False).
        with torch.autocast("cuda", enabled=False):
            atom_feats_list = [
                feats["ref_pos"],
                shardwise_unsqueeze(feats["ref_charge"], -1),
                feats["ref_element"],
            ]
            if not self.use_no_atom_char:
                # ref_atom_name_chars: (B, N, 4, 64) -> (B, N, 256)
                # Dims 2,3 are replicated (not sharded), use shardwise_flatten for non-sharded dims
                atom_feats_list.append(shardwise_flatten(feats["ref_atom_name_chars"], start_dim=2, end_dim=3))

            atom_feats = shardwise_cat(atom_feats_list, dim=-1)
            c = self.embed_atom_features(atom_feats)
            c = _mask_padding_atoms(c, feats["atom_pad_mask"])

            q, c, p = _atom_encoder(
                c=c,
                embed_atompair_ref_pos=self.embed_atompair_ref_pos,
                embed_atompair_ref_dist=self.embed_atompair_ref_dist,
                embed_atompair_mask=self.embed_atompair_mask,
                s_to_c_trans=self.s_to_c_trans if self.structure_prediction else None,
                z_to_p_trans=self.z_to_p_trans if self.structure_prediction else None,
                c_to_p_trans_q=self.c_to_p_trans_q,
                c_to_p_trans_k=self.c_to_p_trans_k,
                p_mlp=self.p_mlp,
                feats=feats,
                s_trunk=s_trunk,
                z=z,
                structure_prediction=self.structure_prediction,
                W=self.atoms_per_window_queries,
                H=self.atoms_per_window_keys,
            )
        return q, c, p


class AtomAttentionEncoder(Module):
    """DTensor AtomAttentionEncoder for window batching.

    Compatible with both Boltz-1x and Boltz-2 serial AtomAttentionEncoder modules.

    V1 (internalized_AtomEncoder=True):
        The V1 serial AtomAttentionEncoder is monolithic — it contains the AtomEncoder
        logic (feature embedding, pair computation) within itself. This DTensor class
        mirrors that by holding all embed/pair modules and calling _atom_encoder() in
        forward to compute q/c/p from raw features.

    V2 (internalized_AtomEncoder=False):
        The V2 architecture splits atom encoding into a separate AtomEncoder class.
        This DTensor class receives pre-computed q/c/bias from the upstream DTensor
        AtomEncoder and only holds r_to_q, AtomTransformer, and atom_to_token_trans.

    Common operations (both V1 and V2):
    1. Apply multiplicity via shardwise_repeat_interleave
    2. Apply r_to_q conditioning (if structure_prediction)
    3. Run AtomTransformer (window-batched)
    4. Atom-to-token scatter via distributed_scatter_reduce with autocast disabled
    """

    def __init__(
        self,
        layer: AtomAttentionEncoderBoltz1 | AtomAttentionEncoderBoltz2,
        device_mesh: DeviceMesh,
    ):
        """Initialize the DTensor-distributed atom attention encoder.

        Parameters
        ----------
        layer : AtomAttentionEncoderBoltz1 | AtomAttentionEncoderBoltz2
            The serial AtomAttentionEncoder module to be distributed.
            Accepts both Boltz-1x and Boltz-2 versions.
        device_mesh : DeviceMesh
            The device mesh for distributed tensor operations.

        Raises
        ------
        TypeError
            If layer is not a recognized type.
        """
        super().__init__()

        if not isinstance(layer, (AtomAttentionEncoderBoltz1, AtomAttentionEncoderBoltz2)):
            raise TypeError(
                ", ".join(
                    [
                        f"Instance {layer} should have type "
                        f"{AtomAttentionEncoderBoltz1} or {AtomAttentionEncoderBoltz2}",
                        f"but instead has type {type(layer)}.",
                    ]
                )
            )

        # V1's serial AtomAttentionEncoder is monolithic: it internalizes the AtomEncoder
        # logic (feature embedding + pair computation) within the same class. In contrast,
        # V2 splits this into a separate AtomEncoder class whose outputs (q, c, bias) are
        # passed into AtomAttentionEncoder from outside.
        #
        # When internalized_AtomEncoder is True (V1), this DTensor class holds all the
        # embed/pair modules and calls _atom_encoder() in forward to compute q/c/p from
        # raw features. When False (V2), q/c/bias are expected as pre-computed inputs
        # from the upstream DTensor AtomEncoder.
        self.internalized_AtomEncoder = isinstance(layer, AtomAttentionEncoderBoltz1)
        self.structure_prediction = layer.structure_prediction

        # V1 only: holds all pair computation modules for _atom_encoder()
        # (V2 delegates these to the separate DTensor AtomEncoder class)
        if self.internalized_AtomEncoder:
            self.atoms_per_window_queries = layer.atoms_per_window_queries
            self.atoms_per_window_keys = layer.atoms_per_window_keys
            validate_window_batching_parameters(
                self.atoms_per_window_queries, self.atoms_per_window_keys, use_window_batching=True
            )

            # V1 atom feature embedding (LinearNoBias)
            self.embed_atom_features = LinearParamsReplicated(
                layer_local=layer.embed_atom_features, device_mesh=device_mesh
            )

            # Pair computation modules (same structure as DTensor AtomEncoder)
            self.embed_atompair_ref_pos = LinearParamsReplicated(
                layer_local=layer.embed_atompair_ref_pos, device_mesh=device_mesh
            )
            self.embed_atompair_ref_dist = LinearParamsReplicated(
                layer_local=layer.embed_atompair_ref_dist, device_mesh=device_mesh
            )
            self.embed_atompair_mask = LinearParamsReplicated(
                layer_local=layer.embed_atompair_mask, device_mesh=device_mesh
            )

            if self.structure_prediction:
                self.s_to_c_trans = nn.Sequential(
                    LayerNormParamsReplicated(layer.s_to_c_trans[0], device_mesh=device_mesh),
                    LinearParamsReplicated(layer_local=layer.s_to_c_trans[1], device_mesh=device_mesh),
                )
                self.z_to_p_trans = nn.Sequential(
                    LayerNormParamsReplicated(layer.z_to_p_trans[0], device_mesh=device_mesh),
                    LinearParamsReplicated(layer_local=layer.z_to_p_trans[1], device_mesh=device_mesh),
                )

            self.c_to_p_trans_q = nn.Sequential(
                nn.ReLU(),
                LinearParamsReplicated(layer_local=layer.c_to_p_trans_q[1], device_mesh=device_mesh),
            )
            self.c_to_p_trans_k = nn.Sequential(
                nn.ReLU(),
                LinearParamsReplicated(layer_local=layer.c_to_p_trans_k[1], device_mesh=device_mesh),
            )
            self.p_mlp = nn.Sequential(
                nn.ReLU(),
                LinearParamsReplicated(layer_local=layer.p_mlp[1], device_mesh=device_mesh),
                nn.ReLU(),
                LinearParamsReplicated(layer_local=layer.p_mlp[3], device_mesh=device_mesh),
                nn.ReLU(),
                LinearParamsReplicated(layer_local=layer.p_mlp[5], device_mesh=device_mesh),
            )

        # V2 only: r_to_q_trans (V1 also has it but under structure_prediction)
        if self.structure_prediction:
            self.r_to_q_trans = LinearParamsReplicated(layer_local=layer.r_to_q_trans, device_mesh=device_mesh)

        # Common: AtomTransformer (window batching)
        self.atom_encoder = AtomTransformer(layer=layer.atom_encoder, device_mesh=device_mesh)

        # Common: atom_to_token_trans Sequential(LinearNoBias, ReLU)
        # Strip nn.ReLU from serial sequential and wrap the Linear
        self.atom_to_token_trans = nn.Sequential(
            LinearParamsReplicated(layer_local=layer.atom_to_token_trans[0], device_mesh=device_mesh),
            nn.ReLU(),
        )

    def forward(
        self,
        feats: dict[str, DTensor],
        q: DTensor | None = None,
        c: DTensor | None = None,
        atom_enc_bias: DTensor | None = None,
        s_trunk: DTensor | None = None,
        z: DTensor | None = None,
        r: DTensor | None = None,
        multiplicity: int = 1,
        model_cache: dict[str, dict[str, DTensor]] | None = None,
    ) -> tuple[DTensor, DTensor, DTensor, DTensor]:
        """Forward pass for the DTensor-distributed atom attention encoder.

        For V2: q, c, atom_enc_bias must be provided (from upstream DTensor AtomEncoder).
        For V1: s_trunk, z are used to compute q, c, p from raw features via _atom_encoder().

        Parameters
        ----------
        feats : dict[str, DTensor]
            Atom features. Must contain atom_pad_mask, atom_to_token_ids_global,
            atom_to_token_local_onehot. For V1 also: ref_pos, ref_space_uid,
            ref_charge, ref_element, ref_atom_name_chars.
        q : DTensor or None
            Pre-computed atom query representation (V2 only), shape (B, N_atoms, atom_s).
        c : DTensor or None
            Pre-computed atom conditioning representation (V2 only), shape (B, N_atoms, atom_s).
        atom_enc_bias : DTensor or None
            Pre-computed pairwise bias (V2 only), shape (B, K, W, H, num_heads * depth).
        s_trunk : DTensor or None
            Token single representation (V1 only), shape (B, N_tokens, token_s).
        z : DTensor or None
            Token pair representation (V1 only), shape (B, N_tokens, N_tokens, token_z).
        r : DTensor or None
            Atom positions for r_to_q conditioning (both V1 and V2),
            shape (B * multiplicity, N_atoms, 3) for V2 or (B * multiplicity, N_atoms, 3) for V1.
        multiplicity : int, optional
            Number of diffusion samples, by default 1.
        model_cache : dict or None, optional
            Model cache for inference optimization (V1 internalized path only).

        Returns
        -------
        tuple[DTensor, DTensor, DTensor, DTensor]
            (a, q, c, p) where:
            - a: Token representation, shape (B * multiplicity, N_tokens_local, token_s)
            - q: Atom query after transformer, shape (B * multiplicity, N_atoms, atom_s)
            - c: Atom conditioning, shape (B * multiplicity, N_atoms, atom_s)
            - p: Pair representation / bias (B, K, W, H, atom_z or num_heads*depth)
        """
        atom_mask = feats["atom_pad_mask"].bool()

        # model_cache is a V1-only feature (internalized path)
        if model_cache is not None and not self.internalized_AtomEncoder:
            raise ValueError("model_cache is only supported with internalized AtomEncoder (V1 path)")

        # V1 model_cache: cache q/c/p from _atom_encoder on first call
        layer_cache = None
        if model_cache is not None:
            cache_prefix = "atomencoder"
            if cache_prefix not in model_cache:
                model_cache[cache_prefix] = {}
            layer_cache = model_cache[cache_prefix]

        if self.internalized_AtomEncoder:
            if model_cache is None or len(layer_cache) == 0:
                # First call or no cache: compute q, c, p from raw features
                atom_feats = shardwise_cat(
                    [
                        feats["ref_pos"],
                        shardwise_unsqueeze(feats["ref_charge"], -1),
                        shardwise_unsqueeze(atom_mask, -1),
                        feats["ref_element"],
                        shardwise_flatten(feats["ref_atom_name_chars"], start_dim=2, end_dim=3),
                    ],
                    dim=-1,
                )
                c = self.embed_atom_features(atom_feats)
                c = _mask_padding_atoms(c, atom_mask)

                q, c, p = _atom_encoder(
                    c=c,
                    embed_atompair_ref_pos=self.embed_atompair_ref_pos,
                    embed_atompair_ref_dist=self.embed_atompair_ref_dist,
                    embed_atompair_mask=self.embed_atompair_mask,
                    s_to_c_trans=self.s_to_c_trans if self.structure_prediction else None,
                    z_to_p_trans=self.z_to_p_trans if self.structure_prediction else None,
                    c_to_p_trans_q=self.c_to_p_trans_q,
                    c_to_p_trans_k=self.c_to_p_trans_k,
                    p_mlp=self.p_mlp,
                    feats=feats,
                    s_trunk=s_trunk,
                    z=z,
                    structure_prediction=self.structure_prediction,
                    W=self.atoms_per_window_queries,
                    H=self.atoms_per_window_keys,
                )

                if model_cache is not None:
                    layer_cache["q"] = q
                    layer_cache["c"] = c
                    layer_cache["p"] = p
            else:
                # Subsequent calls: use cached q/c/p
                q = layer_cache["q"]
                c = layer_cache["c"]
                p = layer_cache["p"]
        else:
            # V2: q, c, bias provided by upstream DTensor AtomEncoder
            if q is None or c is None:
                raise ValueError("V2 AtomAttentionEncoder requires pre-computed q and c from upstream AtomEncoder.")
            p = atom_enc_bias  # may be None if not structure_prediction

        # ====================================================================
        # Common: Apply multiplicity and r_to_q conditioning
        # ====================================================================
        if self.structure_prediction:
            q = shardwise_repeat_interleave(q, multiplicity, 0)

            if self.internalized_AtomEncoder:
                # V1: r_to_q_trans takes (B*M, N, 10) = concat([r, zeros(B*M, N, 7)])
                r_zeros_shape = list(r.shape)
                r_zeros_shape[-1] = 7  # replace last dim 3 → 7
                r_zeros_dt = dtensor_zeros(
                    r_zeros_shape,
                    device_mesh=r.device_mesh,
                    placements=r.placements,
                    dtype=r.dtype,
                    requires_grad=False,
                )
                r_input = shardwise_cat([r, r_zeros_dt], dim=-1)
                r_to_q = self.r_to_q_trans(r_input)
            else:
                # V2: r_to_q_trans takes (B*M, N, 3) directly
                r_to_q = self.r_to_q_trans(r)

            q = elementwise_op(q, r_to_q, ElementwiseOp.SUM)

        c = shardwise_repeat_interleave(c, multiplicity, 0)

        # ====================================================================
        # Common: Run AtomTransformer (window batching)
        # ====================================================================
        # multiplicity=1 because multiplicity is already applied to q, c
        # mask should NOT have multiplicity applied (to match p's batch dim)
        q = self.atom_encoder(
            q=q,
            c=c,
            p=p,
            mask=atom_mask,
            multiplicity=1,  # Multiplicity must be 1 as otherwise meaningless for DTensor window batching code
            model_cache=None,
            pair_mask=None,
        )

        # ====================================================================
        # Common: Atom-to-token scatter aggregation
        # ====================================================================
        # Equivalent to serial's: atom_to_token_mean.T @ q_to_a (mean aggregation)
        # Uses distributed_scatter_reduce with "mean" reduction instead of bmm.
        with torch.autocast("cuda", enabled=False):
            q_to_a = self.atom_to_token_trans(q)

            atom_to_token_ids_global = feats["atom_to_token_ids_global"]
            atom_mask_bool = atom_mask.bool()

            # Apply multiplicity to scatter indices and mask
            atom_to_token_ids_global_mul = shardwise_repeat_interleave(atom_to_token_ids_global, multiplicity, 0)
            atom_mask_bool_mul = shardwise_repeat_interleave(atom_mask_bool, multiplicity, 0)

            # n_tokens_per_shard from atom_to_token_local_onehot (see plan notes)
            n_tokens_per_shard = feats["atom_to_token_local_onehot"].to_local().shape[2]

            # Both serial and distributed now compute exact mean (sum / count).
            # The serial code previously used biased mean: atom_to_token / (count + 1e-6),
            # which was fixed to use count.clamp(min=1) for parity.
            a = distributed_scatter_reduce(
                n_tokens_per_shard,
                1,
                atom_to_token_ids_global_mul,
                q_to_a,
                "mean",
                idx_mask=atom_mask_bool_mul,
                are_ids_contiguous=True,
            )

        return a, q, c, p


class FourierEmbedding(Module):
    """DTensor FourierEmbedding for Context Parallelism.

    Wraps a frozen nn.Linear projection with LinearParamsReplicated.
    Compatible with both Boltz-1x and Boltz-2 serial FourierEmbedding modules.
    """

    def __init__(
        self,
        layer: FourierEmbeddingBoltz1 | FourierEmbeddingBoltz2,
        device_mesh: DeviceMesh,
    ):
        """Initialize the Fourier Embedding layer.

        Parameters
        ----------
        layer : FourierEmbeddingBoltz1 | FourierEmbeddingBoltz2
            The serial Fourier embedding layer.
        device_mesh : DeviceMesh
            The device mesh.

        """
        super().__init__()
        assert isinstance(
            layer, (FourierEmbeddingBoltz1, FourierEmbeddingBoltz2)
        ), f"Expected FourierEmbeddingBoltz1 or FourierEmbeddingBoltz2, got {type(layer)}"
        self.device_mesh = device_mesh
        self.proj = LinearParamsReplicated(layer.proj, self.device_mesh)
        if self.proj.weight.requires_grad or self.proj.bias.requires_grad:
            raise ValueError("Linear layer in FourierEmbedding should not have trainable parameters")

    def forward(self, times: DTensor) -> DTensor:
        """Forward pass of the Fourier embedding layer.

        Parameters
        ----------
        times : DTensor
            The times tensor with shape (B,).
            Placements: (Shard(0), Replicate(), Replicate()).

        Returns
        -------
        DTensor
            The Fourier embedding tensor with shape (B, D).
            Placements: (Shard(0), Replicate(), Replicate()).

        """
        expected_placements = (Shard(0), Replicate(), Replicate())
        if times.placements != expected_placements:
            raise ValueError(f"Times tensor has incorrect placements: {times.placements} != {expected_placements}")

        if times.ndim != 1:
            raise ValueError(f"Times tensor should have shape (B,) but got {times.shape}")

        times = shardwise_unsqueeze(times, dim=1)
        rand_proj = self.proj(times)
        return single_tensor_op(2 * pi * rand_proj, ElementwiseOp.COS)


class SingleConditioning(Module):
    """DTensor SingleConditioning for Context Parallelism.

    Compatible with both Boltz-1x and Boltz-2 serial SingleConditioning modules.
    Handles V2's ``disable_times`` flag: when the serial module was constructed with
    ``disable_times=True``, the fourier_embed / norm_fourier / fourier_to_single
    child modules are absent, and the time-conditioning branch is skipped.
    """

    def __init__(
        self,
        layer: SingleConditioningBoltz1 | SingleConditioningBoltz2,
        device_mesh: DeviceMesh,
    ):
        """Initialize the single conditioning layer with DTensor API.

        Parameters
        ----------
        layer : SingleConditioningBoltz1 | SingleConditioningBoltz2
            The serial single conditioning layer.
        device_mesh : DeviceMesh
            The device mesh.

        """
        super().__init__()
        assert isinstance(
            layer, (SingleConditioningBoltz1, SingleConditioningBoltz2)
        ), f"Expected SingleConditioningBoltz1 or SingleConditioningBoltz2, got {type(layer)}"
        self.device_mesh = device_mesh

        self.norm_single = LayerNormParamsReplicated(layer.norm_single, self.device_mesh)
        self.single_embed = LinearParamsReplicated(layer.single_embed, self.device_mesh)

        # V1 always has fourier time embedding; V2 conditionally creates it based on disable_times.
        self.disable_times = getattr(layer, "disable_times", False)
        if not self.disable_times:
            self.fourier_embed = FourierEmbedding(layer.fourier_embed, self.device_mesh)
            self.norm_fourier = LayerNormParamsReplicated(layer.norm_fourier, self.device_mesh)
            self.fourier_to_single = LinearParamsReplicated(layer.fourier_to_single, self.device_mesh)

        self.transitions = nn.ModuleList([])
        for serial_transition in layer.transitions:
            transition = Transition(
                layer=serial_transition,
                device_mesh=self.device_mesh,
            )
            self.transitions.append(transition)

    def forward(self, times: DTensor, s_trunk: DTensor, s_inputs: DTensor) -> tuple[DTensor, DTensor | None]:
        """Forward pass of the single conditioning layer.

        Parameters
        ----------
        times : DTensor
            The times tensor with shape (B,).
            Placements: (Shard(0), Replicate(), Replicate()).
        s_trunk : DTensor
            The trunk single representation tensor with shape (B, N, D).
            Placements: (Shard(0), Shard(1), Replicate()).
        s_inputs : DTensor
            The inputs single representation tensor with shape (B, N, D).
            Placements: (Shard(0), Shard(1), Replicate()).

        Returns
        -------
        tuple[DTensor, DTensor | None]
            s : DTensor with shape (B, N, 2*D) and placements (Shard(0), Shard(1), Replicate()).
            normed_fourier : DTensor with shape (B, D_fourier) and placements
                (Shard(0), Replicate(), Replicate()), or None if disable_times.

        """
        expected_placements_times = (Shard(0), Replicate(), Replicate())
        if times.placements != expected_placements_times:
            raise ValueError(
                f"Times tensor has incorrect placements: {times.placements} != {expected_placements_times}"
            )
        expected_placements_s = (Shard(0), Shard(1), Replicate())
        if s_trunk.placements != expected_placements_s:
            raise ValueError(
                f"s_trunk tensor has incorrect placements: {s_trunk.placements} != {expected_placements_s}"
            )
        if s_inputs.placements != expected_placements_s:
            raise ValueError(
                f"s_inputs tensor has incorrect placements: {s_inputs.placements} != {expected_placements_s}"
            )

        s = shardwise_cat([s_trunk, s_inputs], dim=-1)
        s = self.single_embed(self.norm_single(s))

        normed_fourier: DTensor | None = None
        if not self.disable_times:
            fourier_embed = self.fourier_embed(times)
            normed_fourier = self.norm_fourier(fourier_embed)
            fourier_to_single = self.fourier_to_single(normed_fourier)
            # fourier_to_single: (B, D) with (S(0), R, R) — broadcast-add to s: (B, N, 2D) with (S(0), S(1), R)
            s = replicate_op(s, fourier_to_single, dim_to_unsqueeze_rhs=1, op=ReplicateOp.ADD)

        for transition in self.transitions:
            s = elementwise_op(transition(s), s, ElementwiseOp.SUM)

        return s, normed_fourier


class PairwiseConditioning(Module):
    """DTensor PairwiseConditioning for Context Parallelism.

    Compatible with both Boltz-1x and Boltz-2 serial PairwiseConditioning modules.
    """

    def __init__(
        self,
        layer: PairwiseConditioningBoltz1 | PairwiseConditioningBoltz2,
        device_mesh: DeviceMesh,
    ):
        """Initialize the pairwise conditioning layer with DTensor API.

        Parameters
        ----------
        layer : PairwiseConditioningBoltz1 | PairwiseConditioningBoltz2
            The serial pairwise conditioning layer.
        device_mesh : DeviceMesh
            The device mesh.

        """
        super().__init__()
        assert isinstance(
            layer, (PairwiseConditioningBoltz1, PairwiseConditioningBoltz2)
        ), f"Expected PairwiseConditioningBoltz1 or PairwiseConditioningBoltz2, got {type(layer)}"
        self.device_mesh = device_mesh

        self.dim_pairwise_init_proj = nn.Sequential(
            LayerNormParamsReplicated(layer.dim_pairwise_init_proj[0], self.device_mesh),
            LinearParamsReplicated(layer.dim_pairwise_init_proj[1], self.device_mesh),
        )

        self.transitions = nn.ModuleList([])
        for serial_transition in layer.transitions:
            transition = Transition(
                layer=serial_transition,
                device_mesh=self.device_mesh,
            )
            self.transitions.append(transition)

    def forward(self, z_trunk: DTensor, token_rel_pos_feats: DTensor) -> DTensor:
        """Forward pass of the pairwise conditioning layer.

        Parameters
        ----------
        z_trunk : DTensor
            The trunk pair representation tensor with shape (B, N, N, D).
            Placements: (Shard(0), Shard(1), Shard(2)).
        token_rel_pos_feats : DTensor
            The token relative position features tensor with shape (B, N, N, D_rel).
            Placements: (Shard(0), Shard(1), Shard(2)).

        Returns
        -------
        DTensor
            The conditioned pair representation tensor with shape (B, N, N, D).
            Placements: (Shard(0), Shard(1), Shard(2)).

        """
        expected_placements = (Shard(0), Shard(1), Shard(2))
        if z_trunk.placements != expected_placements:
            raise ValueError(f"z_trunk tensor has incorrect placements: {z_trunk.placements} != {expected_placements}")
        if token_rel_pos_feats.placements != expected_placements:
            raise ValueError(
                f"token_rel_pos_feats tensor has incorrect placements:"
                f" {token_rel_pos_feats.placements} != {expected_placements}"
            )

        z = shardwise_cat([z_trunk, token_rel_pos_feats], dim=-1)
        z = self.dim_pairwise_init_proj(z)

        for transition in self.transitions:
            z = elementwise_op(transition(z), z, ElementwiseOp.SUM)

        return z
