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

"""DTensor-compatible DiffusionModule and AtomDiffusion for Context Parallelism.

Supports both Boltz-1 and Boltz-2 serial DiffusionModule, controlled by the
``internalized_conditioning`` flag (auto-detected from the serial layer type):

- **Internalized** (Boltz-1): module owns ``pairwise_conditioner``,
  ``AtomAttentionEncoder`` computes q/c/p internally.  Forward takes raw
  ``z_trunk`` + ``relative_position_encoding``.
- **Externalized** (Boltz-2): conditioning is pre-computed by a separate
  ``DiffusionConditioning`` module.  Forward receives a
  ``diffusion_conditioning`` dict.

The token-level DiffusionTransformer uses ring attention (all-to-all), while the
atom-level attention (inside AtomAttentionEncoder/Decoder) uses window-batched
attention in both modes.
"""

import warnings
from copy import deepcopy
from math import exp, sqrt

import torch
from torch import nn
from torch.distributed.device_mesh import DeviceMesh
from torch.distributed.tensor import DTensor, Replicate, Shard, full, zeros

from boltz.data import const
from boltz.distributed.comm import AttentionPairBiasComm, TransposeComm
from boltz.distributed.data.feature.featurizer import pack_atom_features
from boltz.distributed.model.layers.atom_to_token import single_repr_token_to_atom
from boltz.distributed.model.layers.elementwise_op import (
    ElementwiseOp,
    elementwise_op,
    scalar_tensor_op,
    single_tensor_op,
)
from boltz.distributed.model.layers.layernorm import LayerNormParamsReplicated
from boltz.distributed.model.layers.linear import LinearParamsReplicated
from boltz.distributed.model.layers.repeat_interleave import shardwise_repeat_interleave
from boltz.distributed.model.layers.replicate_op import ReplicateOp, replicate_op
from boltz.distributed.model.layers.sharded_op import sharded_sum
from boltz.distributed.model.layers.shardwise_op import shardwise_sum
from boltz.distributed.model.layers.squeeze import shardwise_unsqueeze
from boltz.distributed.model.layers.utils import distributed_pack_and_pad, distributed_unpad_and_unpack
from boltz.distributed.model.loss.diffusion import (
    smooth_lddt_loss,
    smooth_lddt_loss_triton,
    weighted_rigid_align,
)
from boltz.distributed.model.modules.encoders import (
    AtomAttentionDecoder,
    AtomAttentionEncoder,
    PairwiseConditioning,
    SingleConditioning,
)
from boltz.distributed.model.modules.transformers import DiffusionTransformer
from boltz.distributed.model.modules.utils import center_random_augmentation
from boltz.distributed.utils import LayoutRightMap, create_distributed_randn
from boltz.model.modules.diffusion import AtomDiffusion as SerialAtomDiffusionV1
from boltz.model.modules.diffusion import DiffusionModule as SerialDiffusionModuleV1
from boltz.model.modules.diffusionv2 import AtomDiffusion as SerialAtomDiffusionV2
from boltz.model.modules.diffusionv2 import DiffusionModule as SerialDiffusionModuleV2
from boltz.model.modules.utils import default


class DiffusionModule(nn.Module):
    """DTensor DiffusionModule for Context Parallelism.

    Supports both Boltz-1 and Boltz-2 serial DiffusionModule via the
    ``internalized_conditioning`` flag (auto-detected from the serial layer):

    - **Internalized** (Boltz-1): owns ``pairwise_conditioner``, encoder
      computes q/c/p internally.  Forward takes ``z_trunk`` and
      ``relative_position_encoding``.
    - **Externalized** (Boltz-2): receives pre-computed conditioning from
      ``DiffusionConditioning``.  Forward takes ``diffusion_conditioning`` dict.

    In both modes, atom features (``feats``) and ``r_noisy`` are expected in
    **unpacked** layout (with intersperse padding from the CP DTensor data
    loader). The module calls ``pack_atom_features`` and
    ``distributed_pack_and_pad`` internally to form a self-contained pack/unpack
    closure. The ``diffusion_conditioning`` dict (V2 only) arrives in packed
    layout as an inter-layer output from ``DiffusionConditioning`` — this
    coupling is managed between the two modules and is unaffected by data
    loading pipeline changes.

    The token-level DiffusionTransformer uses ring attention and the atom-level
    attention uses window-batched attention.
    """

    def __init__(
        self,
        layer: SerialDiffusionModuleV1 | SerialDiffusionModuleV2,
        device_mesh: DeviceMesh,
        ring_comm: AttentionPairBiasComm | None = None,
    ):
        """Initialize the DTensor DiffusionModule.

        Parameters
        ----------
        layer : SerialDiffusionModuleV1 | SerialDiffusionModuleV2
            The serial DiffusionModule.
        device_mesh : DeviceMesh
            The device mesh for distributed tensor operations.
        ring_comm : AttentionPairBiasComm or None, optional
            Ring communication object for the token-level DiffusionTransformer.

        """
        super().__init__()
        if not isinstance(layer, (SerialDiffusionModuleV1, SerialDiffusionModuleV2)):
            raise TypeError(f"Expected SerialDiffusionModuleV1 or SerialDiffusionModuleV2, got {type(layer)}")

        # Internalized: module owns pairwise_conditioner, encoder computes q/c/p internally.
        # Externalized: receives pre-computed conditioning from DiffusionConditioning.
        self.internalized_conditioning = isinstance(layer, SerialDiffusionModuleV1)

        if isinstance(layer, SerialDiffusionModuleV2):
            warnings.warn(
                "CPU offloading-based activation checkpointing by default is "
                "not used for Boltz-2 so we do not use it in DTensor DiffusionModule.",
                UserWarning,
                stacklevel=2,
            )
        elif isinstance(layer, SerialDiffusionModuleV1):
            warnings.warn(
                "CPU offloading-based activation checkpointing can't be passed to DTensor DiffusionModule via "
                "the input v1 serial layer. We will implement a custom flag in the future to enable it.",
                UserWarning,
                stacklevel=2,
            )

        # Sanity: serial layer must have pairwise_conditioner iff internalized
        if self.internalized_conditioning and not hasattr(layer, "pairwise_conditioner"):
            raise ValueError("internalized_conditioning=True but serial layer has no pairwise_conditioner")
        if not self.internalized_conditioning and hasattr(layer, "pairwise_conditioner"):
            raise ValueError("internalized_conditioning=False but serial layer has pairwise_conditioner")

        self.device_mesh = device_mesh
        self.sigma_data = layer.sigma_data
        self.atoms_per_window_queries = layer.atoms_per_window_queries
        self.atoms_per_window_keys = layer.atoms_per_window_keys
        self.activation_checkpointing = getattr(layer, "activation_checkpointing", False)

        # Common sub-modules (identical attribute names in V1 and V2 serial)
        self.single_conditioner = SingleConditioning(layer.single_conditioner, device_mesh)
        self.atom_attention_encoder = AtomAttentionEncoder(layer.atom_attention_encoder, device_mesh)
        self.s_to_a_linear = nn.Sequential(
            LayerNormParamsReplicated(layer.s_to_a_linear[0], device_mesh),
            LinearParamsReplicated(layer.s_to_a_linear[1], device_mesh),
        )
        self.token_transformer = DiffusionTransformer(layer.token_transformer, device_mesh, ring_comm=ring_comm)
        self.a_norm = LayerNormParamsReplicated(layer.a_norm, device_mesh)
        self.atom_attention_decoder = AtomAttentionDecoder(layer.atom_attention_decoder, device_mesh)

        # Internalized-only: owns pairwise_conditioner
        if self.internalized_conditioning:
            self.pairwise_conditioner = PairwiseConditioning(layer.pairwise_conditioner, device_mesh)

    def forward(
        self,
        s_inputs: DTensor,
        s_trunk: DTensor,
        r_noisy: DTensor,
        times: DTensor,
        feats: dict[str, DTensor],
        # Externalized conditioning (when internalized_conditioning=False)
        diffusion_conditioning: dict[str, DTensor] | None = None,
        # Internalized conditioning inputs (when internalized_conditioning=True)
        z_trunk: DTensor | None = None,
        relative_position_encoding: DTensor | None = None,
        # Common
        multiplicity: int = 1,
        model_cache: dict[str, dict[str, DTensor]] | None = None,
    ) -> DTensor | dict[str, DTensor]:
        """Forward pass of the DTensor DiffusionModule.

        Parameters
        ----------
        s_inputs : DTensor
            Input single representation, shape (B, N, token_s).
            Placements: (Shard(0), Shard(1), Replicate()).
        s_trunk : DTensor
            Trunk single representation, shape (B, N, token_s).
            Placements: (Shard(0), Shard(1), Replicate()).
        r_noisy : DTensor
            Noisy atom coordinates, shape (B*M, N_atoms, 3).
            Placements: (Shard(0), Shard(1), Replicate()).
        times : DTensor
            Time embeddings, shape (B*M,).
            Placements: (Shard(0), Replicate(), Replicate()).
        feats : dict[str, DTensor]
            Unpacked atom features (with intersperse padding from the CP DTensor
            data loader). Must include ``token_pad_mask`` (token-level).
            ``pack_atom_features`` is called internally for both V1 and V2.
        diffusion_conditioning : dict[str, DTensor] or None
            Externalized conditioning (required when internalized_conditioning=False).
            Produced by ``DiffusionConditioning`` in packed layout (inter-layer
            output, not from data loader):
            - "q": (B, N_atoms_packed, atom_s), placements (S(0), S(1), R)
            - "c": (B, N_atoms_packed, atom_s), placements (S(0), S(1), R)
            - "atom_enc_bias": (B, K, W, H, total_enc_heads), placements (S(0), S(1), R)
            - "atom_dec_bias": (B, K, W, H, total_dec_heads), placements (S(0), S(1), R)
            - "token_trans_bias": (B, N, N, total_trans_heads), placements (S(0), S(1), S(2))
        z_trunk : DTensor or None
            Trunk pair representation (required when internalized_conditioning=True).
            Shape (B, N, N, token_z), placements (S(0), S(1), S(2)).
        relative_position_encoding : DTensor or None
            Relative position encoding (required when internalized_conditioning=True).
            Shape (B, N, N, token_z), placements (S(0), S(1), S(2)).
        multiplicity : int
            Number of diffusion samples per batch element.
        model_cache : dict or None, optional
            Model cache for inference optimization (internalized path only).

        Returns
        -------
        DTensor or dict[str, DTensor]
            Internalized: ``{"r_update": DTensor, "token_a": DTensor}``
            Externalized: ``r_update`` DTensor directly.

        """
        # ------------------------------------------------------------------
        # Input sanity checks
        # ------------------------------------------------------------------
        # Conditioning mode consistency
        if self.internalized_conditioning:
            if diffusion_conditioning is not None:
                raise ValueError(
                    "internalized_conditioning: diffusion_conditioning must be None "
                    "(conditioning is computed internally from z_trunk + relative_position_encoding)"
                )
            if z_trunk is None or relative_position_encoding is None:
                raise ValueError("internalized_conditioning: z_trunk and relative_position_encoding are required")
        else:
            if diffusion_conditioning is None:
                raise ValueError("externalized_conditioning: diffusion_conditioning dict is required")
            if z_trunk is not None or relative_position_encoding is not None:
                raise ValueError(
                    "externalized_conditioning: z_trunk and relative_position_encoding must be None "
                    "(conditioning is pre-computed in diffusion_conditioning dict)"
                )

        # Placement checks
        expected_single = (Shard(0), Shard(1), Replicate())
        expected_pair = (Shard(0), Shard(1), Shard(2))
        expected_times = (Shard(0), Replicate(), Replicate())
        if s_inputs.placements != expected_single:
            raise ValueError(f"s_inputs has incorrect placements: {s_inputs.placements} != {expected_single}")
        if s_trunk.placements != expected_single:
            raise ValueError(f"s_trunk has incorrect placements: {s_trunk.placements} != {expected_single}")
        if r_noisy.placements != expected_single:
            raise ValueError(f"r_noisy has incorrect placements: {r_noisy.placements} != {expected_single}")
        if times.placements != expected_times:
            raise ValueError(f"times has incorrect placements: {times.placements} != {expected_times}")
        if self.internalized_conditioning:
            if z_trunk.placements != expected_pair:
                raise ValueError(f"z_trunk has incorrect placements: {z_trunk.placements} != {expected_pair}")
            if relative_position_encoding.placements != expected_pair:
                raise ValueError(
                    f"relative_position_encoding has incorrect placements:"
                    f" {relative_position_encoding.placements} != {expected_pair}"
                )

        # Shape checks: s_inputs/s_trunk batch should NOT include multiplicity
        if s_inputs.shape[0] != feats["token_pad_mask"].shape[0]:
            raise ValueError(
                f"s_inputs batch {s_inputs.shape[0]} != feats['token_pad_mask'] batch"
                f" {feats['token_pad_mask'].shape[0]} (s_inputs should not include multiplicity)"
            )
        if r_noisy.shape[0] != feats["atom_pad_mask"].shape[0] * multiplicity:
            raise ValueError(
                f"r_noisy batch {r_noisy.shape[0]} != atom_pad_mask batch"
                f" {feats['atom_pad_mask'].shape[0]} * multiplicity {multiplicity}"
            )

        # ------------------------------------------------------------------
        # 1. Single conditioning (identical in both paths)
        # ------------------------------------------------------------------
        s_trunk_mult = shardwise_repeat_interleave(s_trunk, multiplicity, 0)
        s_inputs_mult = shardwise_repeat_interleave(s_inputs, multiplicity, 0)
        if self.activation_checkpointing and not self.internalized_conditioning and self.training:
            s, normed_fourier = torch.utils.checkpoint.checkpoint(
                self.single_conditioner, times, s_trunk_mult, s_inputs_mult, use_reentrant=False
            )
        else:
            s, normed_fourier = self.single_conditioner(times, s_trunk_mult, s_inputs_mult)

        # Promote to at least float32 for numerical stability, preserving higher precision
        compute_dtype = torch.promote_types(r_noisy.dtype, torch.float32)

        # ------------------------------------------------------------------
        # Pack atom features and r_noisy (shared by V1 and V2)
        # ------------------------------------------------------------------
        # Atom features (feats) are expected in unpacked layout (with intersperse padding
        # from the CP DTensor data loader). Each module calls pack_atom_features internally
        # to form a self-contained pack/unpack closure, so that no external caller needs to
        # pre-pack features. This ensures all modules accept atom features directly as
        # produced by the data loader, and future refactoring of the data loading pipeline
        # will not require changes to these modules.
        #
        # The pack/unpack lifecycle for window batching:
        # a) upstream input r_noisy and atom feats have interspersed atom padding
        #    due to CP data sharding requirements
        # b) pack_atom_features and distributed_pack_and_pad convert these inputs
        #    to packed format for the window batching of AtomAttentionEncoder and
        #    AtomAttentionDecoder. The q, c and p returned from AtomAttentionEncoder
        #    are also packed.
        # c) AtomAttentionDecoder takes packed q, c and p and packed atom_pad_mask
        #    and atom_to_token_ids_global and outputs r_update_packed.
        # d) distributed_unpad_and_unpack reverts r_update_packed to r_update
        # e) all packed features are discarded after the forward pass.
        # TODO: pack_atom_features and distributed_pack_and_pad should be moved
        #       to the data featurizing layer to save the extra compute and space
        #       due to interspersed padding.
        W = self.atoms_per_window_queries
        _keys_atom_features_packed = {
            "atom_pad_mask",
            "ref_pos",
            "ref_space_uid",
            "ref_charge",
            "ref_element",
            "ref_atom_name_chars",
            "atom_to_token",
        }
        feats_packed = pack_atom_features(feats, _keys_atom_features_packed, W)

        atom_mask_mul = shardwise_repeat_interleave(feats["atom_pad_mask"].bool(), multiplicity, 0)
        atom_mask_mul_expanded = shardwise_unsqueeze(atom_mask_mul, dim=-1)
        r_noisy_packed, atom_mask_r_noisy_packed = distributed_pack_and_pad(r_noisy, atom_mask_mul_expanded, W, axis=1)

        if self.internalized_conditioning:
            # ---- Internalized path (Boltz-1) ----
            # Compute pairwise conditioning z (skipped if cached)
            if model_cache is None or len(model_cache) == 0:
                z = self.pairwise_conditioner(z_trunk, relative_position_encoding)
            else:
                z = None

            # Atom attention encoder: computes q/c/p internally from s_trunk + z
            a, q_skip, c_skip, p_skip = self.atom_attention_encoder(
                feats=feats_packed,
                s_trunk=s_trunk,
                z=z,
                r=r_noisy_packed,
                multiplicity=multiplicity,
                model_cache=model_cache,
            )

            # Token processing (token_pad_mask is token-level, not atom-level — use raw feats)
            a = elementwise_op(a, self.s_to_a_linear(s), ElementwiseOp.SUM)
            mask = feats["token_pad_mask"]
            a = self.token_transformer(
                a, mask=mask.to(a.dtype), s=s, z=z, multiplicity=multiplicity, model_cache=model_cache
            )
            a = self.a_norm(a)

            # Atom attention decoder with internally-computed p_skip
            r_update = self.atom_attention_decoder(
                a=a,
                q=q_skip,
                c=c_skip,
                p=p_skip,
                feats=feats_packed,
                multiplicity=multiplicity,
                model_cache=model_cache,
            )

            # Unpack r_update
            r_update = distributed_unpad_and_unpack(
                r_update, atom_mask_r_noisy_packed, atom_mask_mul_expanded, axis=1, keep_input_padding=False
            )
            return {"r_update": r_update, "token_a": a.detach()}

        else:
            # ---- Externalized path (Boltz-2) ----
            # The diffusion_conditioning dict (q, c, atom_enc_bias, atom_dec_bias,
            # token_trans_bias) is produced by DiffusionConditioning and arrives in
            # packed layout. Unlike atom features from the data loader, these are
            # inter-layer outputs whose format is managed between the producing and
            # consuming modules. This coupling is intentional: changes to the data
            # loading pipeline do not affect the conditioning interface between
            # DiffusionConditioning and DiffusionModule.
            a, q_skip, c_skip, p_skip = self.atom_attention_encoder(
                feats=feats_packed,
                q=diffusion_conditioning["q"].to(compute_dtype),
                c=diffusion_conditioning["c"].to(compute_dtype),
                atom_enc_bias=diffusion_conditioning["atom_enc_bias"].to(compute_dtype),
                r=r_noisy_packed,
                multiplicity=multiplicity,
            )

            # Token processing with pre-computed token_trans_bias
            a = elementwise_op(a, self.s_to_a_linear(s), ElementwiseOp.SUM)
            mask = feats["token_pad_mask"]
            a = self.token_transformer(
                a,
                mask=mask.to(compute_dtype),
                s=s,
                z=diffusion_conditioning["token_trans_bias"].to(compute_dtype),
                multiplicity=multiplicity,
            )
            a = self.a_norm(a)

            # Atom attention decoder with pre-computed atom_dec_bias
            r_update = self.atom_attention_decoder(
                a=a,
                q=q_skip,
                c=c_skip,
                p=diffusion_conditioning["atom_dec_bias"].to(compute_dtype),
                feats=feats_packed,
                multiplicity=multiplicity,
            )

            # Unpack r_update
            r_update = distributed_unpad_and_unpack(
                r_update, atom_mask_r_noisy_packed, atom_mask_mul_expanded, axis=1, keep_input_padding=False
            )
            return r_update


class AtomDiffusion(nn.Module):
    """DTensor AtomDiffusion for Context Parallelism.

    Wraps DiffusionModule with diffusion scheduling (noise preconditioning,
    training forward, and sampling). Scalar diffusion math (c_skip, c_out, c_in,
    c_noise, loss_weight) is identical between V1 and V2.

    Supports both V1 (internalized conditioning) and V2 (externalized conditioning)
    via the ``internalized_conditioning`` flag inherited from the wrapped DiffusionModule.

    Atom features (``feats``) and atom-level tensors (``r_noisy``, ``coords``,
    ``noise``) are expected in **unpacked** layout (with intersperse padding from
    the CP DTensor data loader). The wrapped ``DiffusionModule`` calls
    ``pack_atom_features`` internally — no external packing is needed.
    The ``diffusion_conditioning`` dict (V2 only) arrives in packed layout as an
    inter-layer output from ``DiffusionConditioning``; this format is managed
    between the producing and consuming modules and is unaffected by data loading
    pipeline changes.
    """

    def __init__(
        self,
        layer: SerialAtomDiffusionV1 | SerialAtomDiffusionV2,
        device_mesh: DeviceMesh,
        ring_comm: AttentionPairBiasComm | None = None,
        transpose_comm: TransposeComm | None = None,
    ):
        """Initialize the DTensor AtomDiffusion.

        Parameters
        ----------
        layer : SerialAtomDiffusionV1 | SerialAtomDiffusionV2
            The serial AtomDiffusion module.
        device_mesh : DeviceMesh
            The device mesh for distributed tensor operations.
        ring_comm : AttentionPairBiasComm or None, optional
            Ring communication for the token-level DiffusionTransformer.
        transpose_comm : TransposeComm or None, optional
            Transpose communication for smooth LDDT loss. Required when
            add_smooth_lddt_loss is True in compute_loss.

        """
        super().__init__()
        if not isinstance(layer, (SerialAtomDiffusionV1, SerialAtomDiffusionV2)):
            raise TypeError(f"Expected SerialAtomDiffusionV1 or SerialAtomDiffusionV2, got {type(layer)}")

        self.device_mesh = device_mesh
        self.transpose_comm = transpose_comm

        # Copy diffusion parameters (identical in V1 and V2)
        self.sigma_min = layer.sigma_min
        self.sigma_max = layer.sigma_max
        self.sigma_data = layer.sigma_data
        self.rho = layer.rho
        self.P_mean = layer.P_mean
        self.P_std = layer.P_std
        self.num_sampling_steps = layer.num_sampling_steps
        self.gamma_0 = layer.gamma_0
        self.gamma_min = layer.gamma_min
        self.noise_scale = layer.noise_scale
        self.step_scale = layer.step_scale
        self.coordinate_augmentation = layer.coordinate_augmentation
        self.alignment_reverse_diff = layer.alignment_reverse_diff
        self.synchronize_sigmas = layer.synchronize_sigmas
        self.token_s = layer.token_s

        # Convert the score model to DTensor version
        self.score_model = DiffusionModule(
            layer.score_model,
            device_mesh,
            ring_comm=ring_comm,
        )

        # Derive conditioning mode from the wrapped DiffusionModule
        self.internalized_conditioning = self.score_model.internalized_conditioning

        # V1-only attributes
        self.use_inference_model_cache = getattr(layer, "use_inference_model_cache", False)
        self.accumulate_token_repr = getattr(layer, "accumulate_token_repr", False)
        if self.accumulate_token_repr:
            # v2 doesn't have accumulate_token_repr
            if isinstance(layer, SerialAtomDiffusionV2):
                raise ValueError("accumulate_token_repr should not exist in AtomDiffusionV2")
            # TODO: wrap out_token_feat_update with DTensor layers if accumulate_token_repr is needed
            warnings.warn("OutTokenFeatUpdate is not implemented in DTensor mode. Skipping.")

        self.register_buffer("zero", torch.tensor(0.0), persistent=False)
        self.transpose_comm = deepcopy(transpose_comm)  # for self.compute_loss
        self.v2 = isinstance(layer, SerialAtomDiffusionV2)

    @property
    def device(self):
        """Get the device type of the model."""
        return self.device_mesh.device_type

    # ------------------------------------------------------------------
    # Diffusion preconditioning (DTensor scalar ops)
    # ------------------------------------------------------------------

    def _check_sigma_placement(self, sigma: DTensor) -> None:
        """Validate that sigma has the expected placements (Shard(0), Replicate(), Replicate()).

        Sigma is a 1-D noise-level tensor sharded across the DP axis and
        replicated across CP axes. All preconditioning helpers (c_skip, c_out,
        c_in, c_noise, loss_weight) call this before operating on sigma.
        """
        expected = (Shard(0), Replicate(), Replicate())
        if sigma.placements != expected:
            raise ValueError(f"Sigma tensor has incorrect placements: {sigma.placements} != {expected}")

    def c_skip(self, sigma: DTensor) -> DTensor:
        """Skip-connection scaling: sigma_data^2 / (sigma^2 + sigma_data^2).

        Weights the direct pass-through of noised coordinates in the
        preconditioning formula: denoised = c_skip * noised + c_out * net_out.

        Parameters
        ----------
        sigma : DTensor
            Noise levels, shape (B*M,). Placements: (Shard(0), Replicate(), Replicate()).

        Returns
        -------
        DTensor
            Same shape and placements as sigma.

        """
        self._check_sigma_placement(sigma)
        sigma_sq = scalar_tensor_op(2, sigma, ElementwiseOp.POW)
        denom = scalar_tensor_op(self.sigma_data**2, sigma_sq, ElementwiseOp.SUM)
        return scalar_tensor_op(self.sigma_data**2, denom, ElementwiseOp.DIV)

    def c_out(self, sigma: DTensor) -> DTensor:
        """Output scaling: sigma * sigma_data / sqrt(sigma^2 + sigma_data^2).

        Weights the network output in the preconditioning formula:
        denoised = c_skip * noised + c_out * net_out.

        Parameters
        ----------
        sigma : DTensor
            Noise levels, shape (B*M,). Placements: (Shard(0), Replicate(), Replicate()).

        Returns
        -------
        DTensor
            Same shape and placements as sigma.

        """
        self._check_sigma_placement(sigma)
        numer = scalar_tensor_op(self.sigma_data, sigma, ElementwiseOp.PROD)
        sigma_sq = scalar_tensor_op(2, sigma, ElementwiseOp.POW)
        denom = scalar_tensor_op(self.sigma_data**2, sigma_sq, ElementwiseOp.SUM)
        denom = scalar_tensor_op(0.5, denom, ElementwiseOp.POW)
        return elementwise_op(numer, denom, ElementwiseOp.DIV)

    def c_in(self, sigma: DTensor) -> DTensor:
        """Input scaling: 1 / sqrt(sigma^2 + sigma_data^2).

        Normalizes noised coordinates before feeding into the score model
        in preconditioned_network_forward.

        Parameters
        ----------
        sigma : DTensor
            Noise levels, shape (B*M,). Placements: (Shard(0), Replicate(), Replicate()).

        Returns
        -------
        DTensor
            Same shape and placements as sigma.

        """
        self._check_sigma_placement(sigma)
        sigma_sq = scalar_tensor_op(2, sigma, ElementwiseOp.POW)
        denom = scalar_tensor_op(self.sigma_data**2, sigma_sq, ElementwiseOp.SUM)
        denom = scalar_tensor_op(0.5, denom, ElementwiseOp.POW)
        return scalar_tensor_op(1, denom, ElementwiseOp.DIV)

    def c_noise(self, sigma: DTensor) -> DTensor:
        """Noise conditioning: log(sigma / sigma_data) * 0.25.

        Produces the time embedding input for the score model's
        SingleConditioning / FourierEmbedding layers.

        Parameters
        ----------
        sigma : DTensor
            Noise levels, shape (B*M,). Placements: (Shard(0), Replicate(), Replicate()).

        Returns
        -------
        DTensor
            Same shape and placements as sigma.

        """
        self._check_sigma_placement(sigma)
        scaled = scalar_tensor_op(1 / self.sigma_data, sigma, ElementwiseOp.PROD)
        scaled_local = scaled.to_local().clamp(min=1e-20)
        scaled = DTensor.from_local(
            scaled_local,
            device_mesh=scaled.device_mesh,
            placements=scaled.placements,
            shape=scaled.shape,
            stride=scaled.stride(),
        )
        log_sigma = single_tensor_op(scaled, ElementwiseOp.LOG)
        return scalar_tensor_op(0.25, log_sigma, ElementwiseOp.PROD)

    def loss_weight(self, sigma: DTensor) -> DTensor:
        """Diffusion loss weighting: (sigma^2 + sigma_data^2) / (sigma * sigma_data)^2.

        Used by compute_loss to weight the MSE loss at each noise level.

        Parameters
        ----------
        sigma : DTensor
            Noise levels, shape (B*M,). Placements: (Shard(0), Replicate(), Replicate()).

        Returns
        -------
        DTensor
            Same shape and placements as sigma.

        """
        self._check_sigma_placement(sigma)
        sigma_sq = scalar_tensor_op(2, sigma, ElementwiseOp.POW)
        numer = scalar_tensor_op(self.sigma_data**2, sigma_sq, ElementwiseOp.SUM)
        denom = scalar_tensor_op(self.sigma_data**2, sigma_sq, ElementwiseOp.PROD)
        return elementwise_op(numer, denom, ElementwiseOp.DIV)

    def noise_distribution(self, batch_size: int, dtype: torch.dtype = torch.float32) -> DTensor:
        """Sample noise levels from the training distribution.

        Generates sigma_data * exp(P_mean + P_std * randn(batch_size)).
        Called by forward() to produce per-sample noise levels for the
        diffusion training step.

        Parameters
        ----------
        batch_size : int
            Number of samples (typically B*M after multiplicity expansion).
        dtype : torch.dtype, optional
            Dtype for the generated noise levels. Should match the model's
            compute dtype (e.g. feats["coords"].dtype). Default torch.float32.

        Returns
        -------
        DTensor
            Noise levels, shape (batch_size,).
            Placements: (Shard(0), Replicate(), Replicate()).

        """
        noise = create_distributed_randn(
            (batch_size,),
            device_mesh=self.device_mesh,
            placements=(Shard(0), Replicate(), Replicate()),
            dtype=dtype,
        )
        noise = scalar_tensor_op(self.P_std, noise, ElementwiseOp.PROD)
        noise = single_tensor_op(noise, ElementwiseOp.EXP)
        noise = scalar_tensor_op(self.sigma_data * exp(self.P_mean), noise, ElementwiseOp.PROD)
        return noise

    # ------------------------------------------------------------------
    # Preconditioned network forward
    # ------------------------------------------------------------------

    def preconditioned_network_forward(
        self,
        noised_atom_coords: DTensor,
        sigma: float | DTensor,
        network_condition_kwargs: dict,
    ) -> DTensor | tuple[DTensor, DTensor]:
        """Preconditioned forward pass: c_skip * x + c_out * score_model(c_in * x, c_noise).

        Parameters
        ----------
        noised_atom_coords : DTensor
            Noisy atom coordinates, shape (B*M, N_atoms, 3).
            Placements: (Shard(0), Shard(1), Replicate()).
        sigma : float or DTensor
            Noise level. If float, broadcast to all batch elements.
            If DTensor, shape (B*M,) with placements (Shard(0), Replicate(), Replicate()).
        network_condition_kwargs : dict
            Conditioning arguments for the score model.

        Returns
        -------
        DTensor or tuple[DTensor, DTensor]
            Internalized: ``(denoised_coords, token_a)``
            Externalized: ``denoised_coords``

        """
        batch_size = noised_atom_coords.shape[0]

        if isinstance(sigma, float):
            sigma = full(
                (batch_size,),
                sigma,
                dtype=noised_atom_coords.dtype,
                device_mesh=self.device_mesh,
                placements=(Shard(0), Replicate(), Replicate()),
            )

        # Expand sigma to (B, 3) for element-wise multiply with (B, N, 3)
        padded_sigma = shardwise_repeat_interleave(shardwise_unsqueeze(sigma, dim=-1), 3, -1)

        # r_noisy = c_in(sigma) * noised_atom_coords
        r_noisy = replicate_op(noised_atom_coords, self.c_in(padded_sigma), 1, ReplicateOp.PROD)
        times = self.c_noise(sigma)

        net_out = self.score_model(
            r_noisy=r_noisy,
            times=times,
            **network_condition_kwargs,
        )

        if self.internalized_conditioning:
            # V1: score_model returns dict {"r_update": ..., "token_a": ...}
            r_update = net_out["r_update"]
            token_a = net_out["token_a"]
        else:
            # V2: score_model returns r_update DTensor directly
            r_update = net_out

        # denoised = c_skip(sigma) * noised_atom_coords + c_out(sigma) * r_update
        skip_term = replicate_op(noised_atom_coords, self.c_skip(padded_sigma), 1, ReplicateOp.PROD)
        out_term = replicate_op(r_update, self.c_out(padded_sigma), 1, ReplicateOp.PROD)
        denoised_coords = elementwise_op(skip_term, out_term, ElementwiseOp.SUM)

        if self.internalized_conditioning:
            return denoised_coords, token_a
        return denoised_coords

    # ------------------------------------------------------------------
    # Sampling schedule
    # ------------------------------------------------------------------

    def sample_schedule(self, num_sampling_steps: int | None = None) -> torch.Tensor:
        """Generate sigma schedule for sampling. Returns plain Tensor (scalar schedule)."""
        num_sampling_steps = default(num_sampling_steps, self.num_sampling_steps)
        if num_sampling_steps < 2:
            raise ValueError(f"Need at least 2 sampling steps, got {num_sampling_steps}")
        inv_rho = 1 / self.rho
        steps = torch.arange(num_sampling_steps, device=self.device, dtype=torch.float32)
        sigmas = (
            self.sigma_max**inv_rho
            + steps / (num_sampling_steps - 1) * (self.sigma_min**inv_rho - self.sigma_max**inv_rho)
        ) ** self.rho
        sigmas = sigmas * self.sigma_data
        sigmas = torch.nn.functional.pad(sigmas, (0, 1), value=0.0)
        return sigmas

    # ------------------------------------------------------------------
    # Training forward
    # ------------------------------------------------------------------

    def forward(
        self,
        s_inputs: DTensor,
        s_trunk: DTensor,
        feats: dict[str, DTensor],
        # Externalized conditioning (when internalized_conditioning=False)
        diffusion_conditioning: dict[str, DTensor] | None = None,
        # Internalized conditioning inputs (when internalized_conditioning=True)
        z_trunk: DTensor | None = None,
        relative_position_encoding: DTensor | None = None,
        # Common
        multiplicity: int = 1,
    ) -> dict[str, DTensor]:
        """Training forward: add noise, run preconditioned network, return denoised coords.

        Parameters
        ----------
        s_inputs : DTensor
            Input single representation, shape (B, N, token_s).
        s_trunk : DTensor
            Trunk single representation, shape (B, N, token_s).
        feats : dict[str, DTensor]
            Pre-packed atom features. Must contain 'coords' and 'atom_pad_mask'.
        diffusion_conditioning : dict[str, DTensor] or None
            Externalized conditioning (required when internalized_conditioning=False).
        z_trunk : DTensor or None
            Trunk pair representation (required when internalized_conditioning=True).
        relative_position_encoding : DTensor or None
            Relative position encoding (required when internalized_conditioning=True).
        multiplicity : int
            Number of diffusion samples per batch element.

        Returns
        -------
        dict[str, DTensor]
            denoised_atom_coords, sigmas, aligned_true_atom_coords.

        """
        # Sanity: ensure forward args and feats shapes match conditioning mode
        coords = feats["coords"]
        atom_pad_mask = feats["atom_pad_mask"]
        B, N_atoms = atom_pad_mask.shape[0], atom_pad_mask.shape[1]
        expected_coords_shape = (B * multiplicity, N_atoms, 3)
        mode = "V1" if self.internalized_conditioning else "V2"
        if self.internalized_conditioning:
            if diffusion_conditioning is not None:
                raise ValueError("internalized_conditioning: diffusion_conditioning must be None")
            if z_trunk is None or relative_position_encoding is None:
                raise ValueError("internalized_conditioning: z_trunk and relative_position_encoding are required")
        else:
            if diffusion_conditioning is None:
                raise ValueError("externalized_conditioning: diffusion_conditioning dict is required")
            if z_trunk is not None or relative_position_encoding is not None:
                raise ValueError("externalized_conditioning: z_trunk and relative_position_encoding must be None")
        if coords.shape != expected_coords_shape:
            raise ValueError(
                f"{mode}: feats['coords'] expected shape {expected_coords_shape} "
                f"(B*M, N_atoms, 3) from atom_pad_mask (B={B}, N={N_atoms}) and multiplicity={multiplicity}, "
                f"got {coords.shape}. Caller must expand coords with multiplicity before calling forward()."
            )
        batch_size = B

        coords_dtype = feats["coords"].dtype
        if self.synchronize_sigmas:
            sigmas = self.noise_distribution(batch_size, dtype=coords_dtype)
            sigmas = shardwise_repeat_interleave(sigmas, multiplicity, 0)
        else:
            sigmas = self.noise_distribution(batch_size * multiplicity, dtype=coords_dtype)

        padded_sigmas = shardwise_repeat_interleave(shardwise_unsqueeze(sigmas, dim=-1), 3, -1)

        # Process atom coordinates
        atom_coords = feats["coords"]
        atom_mask = feats["atom_pad_mask"]
        atom_mask = shardwise_repeat_interleave(atom_mask, multiplicity, 0)

        atom_coords = center_random_augmentation(atom_coords, atom_mask, augmentation=self.coordinate_augmentation)

        noise = create_distributed_randn(
            atom_coords.shape,
            device_mesh=self.device_mesh,
            placements=atom_coords.placements,
            dtype=atom_coords.dtype,
        )

        # Add noise: noised = coords + sigma * noise
        noised_atom_coords = elementwise_op(
            atom_coords,
            replicate_op(noise, padded_sigmas, 1, ReplicateOp.PROD),
            ElementwiseOp.SUM,
        )

        # Build network_condition_kwargs based on conditioning mode
        if self.internalized_conditioning:
            network_condition_kwargs = {
                "s_inputs": s_inputs,
                "s_trunk": s_trunk,
                "feats": feats,
                "multiplicity": multiplicity,
                "z_trunk": z_trunk,
                "relative_position_encoding": relative_position_encoding,
            }
        else:
            network_condition_kwargs = {
                "s_inputs": s_inputs,
                "s_trunk": s_trunk,
                "feats": feats,
                "multiplicity": multiplicity,
                "diffusion_conditioning": diffusion_conditioning,
            }

        # Preconditioned network forward
        precond_result = self.preconditioned_network_forward(
            noised_atom_coords,
            sigmas,
            network_condition_kwargs=network_condition_kwargs,
        )

        # V1 returns (denoised, token_a), V2 returns denoised directly
        if self.internalized_conditioning:
            denoised_atom_coords = precond_result[0]
        else:
            denoised_atom_coords = precond_result

        result = {
            "denoised_atom_coords": denoised_atom_coords,
            "sigmas": sigmas,
            "aligned_true_atom_coords": atom_coords,
        }
        if self.internalized_conditioning:
            result["noised_atom_coords"] = noised_atom_coords
        return result

    # ------------------------------------------------------------------
    # Sampling (inference)
    # ------------------------------------------------------------------

    def sample(
        self,
        atom_mask: DTensor,
        num_sampling_steps: int | None = None,
        multiplicity: int = 1,
        max_parallel_samples: int | None = None,
        train_accumulate_token_repr: bool = False,
        **network_condition_kwargs,
    ) -> dict[str, DTensor | None]:
        """Sample from the diffusion model (inference denoising loop).

        Parameters
        ----------
        atom_mask : DTensor
            Atom mask, shape (B, N_atoms).
            Placements: (Shard(0), Shard(1), Replicate()).
        num_sampling_steps : int or None, optional
            Number of sampling steps. If None, uses default.
        multiplicity : int, optional
            Multiplicity factor, by default 1.
        max_parallel_samples : int or None, optional
            Maximum multiplicity samples processed per chunk. If None,
            all samples are processed in a single call.
        train_accumulate_token_repr : bool, optional
            V1-specific flag for accumulating token representations.
            Not yet implemented in DTensor mode; raises NotImplementedError if True.
        **network_condition_kwargs
            Additional conditioning. For externalized: s_inputs, s_trunk, feats,
            diffusion_conditioning. For internalized: s_inputs, s_trunk, feats,
            z_trunk, relative_position_encoding.

        Returns
        -------
        dict[str, DTensor | None]
            sample_atom_coords: denoised coordinates, diff_token_repr: always None
            (accumulate_token_repr not yet implemented in DTensor mode).

        """
        if train_accumulate_token_repr:
            raise NotImplementedError("train_accumulate_token_repr not implemented in DTensor mode yet")

        if max_parallel_samples is None:
            max_parallel_samples = multiplicity

        num_sampling_steps = default(num_sampling_steps, self.num_sampling_steps)
        atom_mask = shardwise_repeat_interleave(atom_mask, multiplicity, 0)
        shape = (*atom_mask.shape, 3)

        # Sampling schedule (plain Tensor, deterministic)
        sigmas = self.sample_schedule(num_sampling_steps)
        gammas = torch.where(sigmas > self.gamma_min, self.gamma_0, 0.0)
        sigmas_and_gammas = list(zip(sigmas[:-1], sigmas[1:], gammas[1:]))

        # Initial noise
        init_sigma = sigmas[0].item()
        atom_coords = create_distributed_randn(
            shape,
            device_mesh=self.device_mesh,
            placements=atom_mask.placements,
            scale=init_sigma,
        )
        atom_coords_denoised = None

        # V1: model_cache for inference optimization
        model_cache = {} if self.use_inference_model_cache else None

        # V1: token representation tracking for accumulate_token_repr
        token_repr = None
        token_a = None

        # Denoising loop
        for step_idx, (sigma_tm, sigma_t, gamma) in enumerate(sigmas_and_gammas):
            aug = self.coordinate_augmentation
            result = center_random_augmentation(
                atom_coords,
                atom_mask,
                augmentation=aug,
                return_second_coords=True,
                second_coords=atom_coords_denoised,
                return_roto=aug,
            )
            if aug:
                atom_coords, atom_coords_denoised, _ = result
            else:
                atom_coords, atom_coords_denoised = result

            sigma_tm, sigma_t, gamma = sigma_tm.item(), sigma_t.item(), gamma.item()
            t_hat = sigma_tm * (1 + gamma)
            noise_var = self.noise_scale**2 * (t_hat**2 - sigma_tm**2)
            eps = create_distributed_randn(
                shape,
                device_mesh=self.device_mesh,
                placements=atom_mask.placements,
                scale=sqrt(noise_var),
            )
            atom_coords_noisy = elementwise_op(atom_coords, eps, ElementwiseOp.SUM)

            with torch.no_grad():
                placements = atom_coords_noisy.placements
                noisy_local = atom_coords_noisy.to_local()
                denoised_local = torch.zeros_like(noisy_local)
                if noisy_local.shape[0] % multiplicity != 0:
                    # this should only happen if all upstream DTensor modules have removed
                    # the non-even sharding requirements and that we actual have a unevenly
                    # sharded batch dimension
                    raise ValueError(
                        f"noisy_local.shape[0] is not divisible by multiplicity: {noisy_local.shape[0]} % {multiplicity} = {noisy_local.shape[0] % multiplicity}"
                    )
                B_local = noisy_local.shape[0] // multiplicity

                sample_ids = torch.arange(multiplicity, device=self.device)
                n_chunks = (multiplicity + max_parallel_samples - 1) // max_parallel_samples
                sample_ids_chunks = sample_ids.chunk(n_chunks)

                for sample_ids_chunk in sample_ids_chunks:
                    # noisy_local is (B_local*M, N, 3).
                    # Unflatten to (B_local, M, N, 3) and index the M axis so each
                    # chunk selects the correct multiplicity slices, then reflatten
                    # to (B_local*chunk_M, N, 3) and rebuild as DTensor.
                    chunk_M = sample_ids_chunk.numel()
                    noisy_chunk_local = noisy_local.unflatten(0, (B_local, multiplicity))[:, sample_ids_chunk].flatten(
                        0, 1
                    )
                    chunk_global_shape = (
                        atom_coords_noisy.shape[0] * chunk_M // multiplicity,
                        atom_coords_noisy.shape[1],
                        3,
                    )
                    noisy_chunk_dt = DTensor.from_local(
                        noisy_chunk_local,
                        device_mesh=self.device_mesh,
                        placements=placements,
                        shape=chunk_global_shape,
                        stride=LayoutRightMap(chunk_global_shape).strides,
                    )

                    precond_kwargs = dict(multiplicity=chunk_M, **network_condition_kwargs)
                    if model_cache is not None:
                        precond_kwargs["model_cache"] = model_cache

                    precond_result = self.preconditioned_network_forward(
                        noisy_chunk_dt,
                        t_hat,
                        network_condition_kwargs=precond_kwargs,
                    )

                    if self.internalized_conditioning:
                        denoised_chunk_dt, token_a = precond_result
                    else:
                        denoised_chunk_dt = precond_result

                    denoised_local.unflatten(0, (B_local, multiplicity))[:, sample_ids_chunk] = (
                        denoised_chunk_dt.to_local().unflatten(0, (B_local, chunk_M))
                    )

                atom_coords_denoised = DTensor.from_local(
                    denoised_local,
                    device_mesh=self.device_mesh,
                    placements=placements,
                    shape=atom_coords_noisy.shape,
                    stride=atom_coords_noisy.stride(),
                )

            # TODO: accumulate_token_repr support (requires DTensor wrapping of OutTokenFeatUpdate)

            # Alignment reverse diffusion: align noisy coords to denoised coords
            if self.alignment_reverse_diff:
                atom_coords_noisy = weighted_rigid_align(
                    atom_coords_noisy,
                    atom_coords_denoised,
                    atom_mask,
                    atom_mask,
                )

            # Next step: x_{t+1} = x_noisy + step_scale * (sigma_t - t_hat) * (x_noisy - x_denoised) / t_hat
            denoised_over_sigma = scalar_tensor_op(
                1 / t_hat,
                elementwise_op(atom_coords_noisy, atom_coords_denoised, ElementwiseOp.SUB),
                ElementwiseOp.PROD,
            )
            atom_coords = elementwise_op(
                atom_coords_noisy,
                scalar_tensor_op(self.step_scale * (sigma_t - t_hat), denoised_over_sigma, ElementwiseOp.PROD),
                ElementwiseOp.SUM,
            )

        return {"sample_atom_coords": atom_coords, "diff_token_repr": token_repr}

    # ------------------------------------------------------------------
    # Compute loss
    # ------------------------------------------------------------------

    def compute_loss(
        self,
        feats: dict[str, DTensor],
        out_dict: dict[str, DTensor],
        add_smooth_lddt_loss: bool = True,
        nucleotide_loss_weight: float = 5.0,
        ligand_loss_weight: float = 10.0,
        multiplicity: int = 1,
        filter_by_plddt: float = 0.0,
        use_triton_kernel: bool = True,
    ) -> dict[str, DTensor]:
        """Compute loss for the diffusion model.

        Parameters
        ----------
        feats : dict[str, DTensor]
            Features for the diffusion model.
        out_dict : dict[str, DTensor]
            Output dictionary from the diffusion model.
        add_smooth_lddt_loss : bool, optional
            Whether to add smooth LDDT loss.
        nucleotide_loss_weight : float, optional
            Weight for nucleotide loss.
        ligand_loss_weight : float, optional
            Weight for ligand loss.
        multiplicity : int, optional
            Multiplicity factor, by default 1.
        filter_by_plddt : float, optional
            Filter by pLDDT threshold, by default 0.0.
        use_triton_kernel : bool, optional
            Whether to use Triton kernel for smooth LDDT loss.

        Returns
        -------
        dict[str, DTensor]
            Loss dictionary containing "loss" and "loss_breakdown" keys.
            "loss" is the total loss tensor.
            "loss_breakdown" is a dictionary containing "mse_loss" and "smooth_lddt_loss" keys.
            "mse_loss" is the MSE loss tensor.
            "smooth_lddt_loss" is the smooth LDDT loss tensor.
        """
        if not self.v2 and filter_by_plddt != 0.0:
            raise ValueError("filter_by_plddt is only supported for V2")

        with torch.autocast("cuda", enabled=False):
            safe_dtype = torch.promote_types(torch.float32, out_dict["denoised_atom_coords"].dtype)
            denoised_atom_coords = out_dict["denoised_atom_coords"].to(dtype=safe_dtype)
            sigmas = out_dict["sigmas"].to(dtype=safe_dtype)

            resolved_atom_mask_uni = feats["atom_resolved_mask"].to(dtype=safe_dtype)
            resolved_atom_mask = shardwise_repeat_interleave(resolved_atom_mask_uni, multiplicity, 0)

            if self.v2 and filter_by_plddt > 0:
                if "plddt" not in feats:
                    raise RuntimeError("Missing required plddt data in feats for plddt filtering")
                plddt_mask = scalar_tensor_op(filter_by_plddt, feats["plddt"], ElementwiseOp.LT)
                resolved_atom_mask_uni_plddt_masked = elementwise_op(
                    resolved_atom_mask_uni, plddt_mask.to(dtype=safe_dtype), ElementwiseOp.PROD
                )
                resolved_atom_mask_plddt_masked = shardwise_repeat_interleave(
                    resolved_atom_mask_uni_plddt_masked, multiplicity, 0
                )
            else:
                resolved_atom_mask_uni_plddt_masked = resolved_atom_mask_uni
                resolved_atom_mask_plddt_masked = resolved_atom_mask

            atom_type = single_repr_token_to_atom(
                feats["mol_type"].to(dtype=safe_dtype), feats["atom_to_token"].to(dtype=safe_dtype)
            )
            atom_type_mult = shardwise_repeat_interleave(atom_type, multiplicity, 0)

            is_nucleotide_mult = elementwise_op(
                scalar_tensor_op(
                    const.chain_type_ids["DNA"],
                    atom_type_mult,
                    ElementwiseOp.EQUAL,
                ),
                scalar_tensor_op(
                    const.chain_type_ids["RNA"],
                    atom_type_mult,
                    ElementwiseOp.EQUAL,
                ),
                ElementwiseOp.SUM,  # or equivalently OR
            )

            nucleotide_loss_weights = scalar_tensor_op(
                nucleotide_loss_weight,
                is_nucleotide_mult,
                ElementwiseOp.PROD,
            )

            ligand_loss_weights = scalar_tensor_op(
                ligand_loss_weight,
                scalar_tensor_op(
                    const.chain_type_ids["NONPOLYMER"],
                    atom_type_mult,
                    ElementwiseOp.EQUAL,
                ),
                ElementwiseOp.PROD,
            )
            align_weights = scalar_tensor_op(
                1.0,
                elementwise_op(
                    nucleotide_loss_weights,
                    ligand_loss_weights,
                    ElementwiseOp.SUM,
                ),
                ElementwiseOp.SUM,
            )

            with torch.no_grad():
                atom_coords = out_dict["aligned_true_atom_coords"].to(dtype=safe_dtype)
                atom_coords_aligned_ground_truth = weighted_rigid_align(
                    atom_coords,
                    denoised_atom_coords,
                    align_weights.to(dtype=safe_dtype),
                    mask=resolved_atom_mask.to(dtype=safe_dtype),
                )
            # Cast back
            atom_coords_aligned_ground_truth: DTensor = atom_coords_aligned_ground_truth.to(
                dtype=denoised_atom_coords.dtype
            )

            # Weighted MSE loss of denoised atom positions (match serial v2 formula)
            mse_loss = elementwise_op(denoised_atom_coords, atom_coords_aligned_ground_truth, ElementwiseOp.SUB)
            mse_loss = scalar_tensor_op(2.0, mse_loss, ElementwiseOp.POW)
            mse_loss = shardwise_sum(mse_loss, dim=-1)
            mse_loss = elementwise_op(mse_loss, resolved_atom_mask_plddt_masked, ElementwiseOp.PROD)

            resolved_align_weights = elementwise_op(align_weights, resolved_atom_mask_plddt_masked, ElementwiseOp.PROD)
            denom = sharded_sum(
                scalar_tensor_op(3.0, resolved_align_weights, ElementwiseOp.PROD),
                dim=-1,
            )
            if self.v2:
                denom = scalar_tensor_op(1e-5, denom, ElementwiseOp.SUM)

            mse_loss = elementwise_op(mse_loss, resolved_align_weights, ElementwiseOp.PROD)
            mse_loss = sharded_sum(mse_loss, dim=-1)
            mse_loss = elementwise_op(mse_loss, denom, ElementwiseOp.DIV)
            loss_weights = self.loss_weight(sigmas)

            mse_loss = elementwise_op(mse_loss, loss_weights, ElementwiseOp.PROD)
            mse_loss = scalar_tensor_op(
                1.0 / mse_loss.shape[0],
                sharded_sum(mse_loss, dim=0),
                ElementwiseOp.PROD,
            )
            total_loss = mse_loss

            if add_smooth_lddt_loss:
                is_nucleotide = elementwise_op(
                    scalar_tensor_op(const.chain_type_ids["DNA"], atom_type, ElementwiseOp.EQUAL),
                    scalar_tensor_op(const.chain_type_ids["RNA"], atom_type, ElementwiseOp.EQUAL),
                    ElementwiseOp.SUM,
                )
                loss_func = (
                    smooth_lddt_loss_triton
                    if use_triton_kernel and self.device_mesh.device_type == "cuda"
                    else smooth_lddt_loss
                )

                lddt_loss = loss_func(
                    denoised_atom_coords,
                    atom_coords,
                    is_nucleotide=is_nucleotide,
                    coords_mask=resolved_atom_mask_uni_plddt_masked,
                    comm=self.transpose_comm,
                    multiplicity=multiplicity,
                    v2=self.v2,
                )
                total_loss = elementwise_op(total_loss, lddt_loss, ElementwiseOp.SUM)
            else:
                lddt_loss = zeros(
                    total_loss.shape,
                    requires_grad=False,
                    device_mesh=total_loss.device_mesh,
                    placements=total_loss.placements,
                )

            loss_breakdown = {
                "mse_loss": mse_loss,
                "smooth_lddt_loss": lddt_loss,
            }

        return {"loss": total_loss, "loss_breakdown": loss_breakdown}
