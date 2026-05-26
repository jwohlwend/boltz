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

"""Distributed Boltz2 model wrapper for context-parallel training and inference.

This module wraps a serial :class:`~boltz.model.models.boltz2.Boltz2` model
and replaces its submodules with DTensor-aware distributed counterparts,
enabling context parallelism (CP) across sequence/pair dimensions.

Architecture overview
---------------------
The wrapper follows the same pattern as Boltz1Distributed
(``other_versions/boltz-1x-cp/src/boltz/distributed/model/model.py``):

1. Accept a fully-initialised serial ``Boltz2`` instance.
2. Replace each serial submodule with its distributed counterpart
   (e.g. ``LinearParamsReplicated``, distributed ``MSAModule``, etc.).
3. Re-implement ``forward`` to use DTensor operations for outer products,
   recycling, and trunk computation.
4. Provide training/validation/predict steps that handle DTensor
   ↔ plain-tensor conversions for loss computation and logging.

Submodule availability
----------------------
Some serial submodules do **not** yet have distributed implementations.
These are wrapped in ``_PlaceholderModule`` which raises
``NotImplementedError`` when their code path is hit during a forward
pass.  The table below summarises the status:

+----------------------------+----------------------------------------------+
| Serial submodule           | Distributed status                           |
+============================+==============================================+
| s_init, z_init_1, z_init_2 | LinearParamsReplicated (ready)                |
| s_norm, z_norm              | LayerNormParamsReplicated (ready)             |
| s_recycle, z_recycle        | LinearParamsReplicated (ready)                |
| token_bonds                 | LinearParamsReplicated (ready)                |
| msa_module                  | MSAModule (trunkv2 distributed, ready)        |
| pairformer_module           | PairformerModule (ready)                      |
| distogram_module            | DistogramModule (trunkv2.py, ready)            |
| rel_pos                     | RelativePositionEncoder (encoders.py, ready)  |
| contact_conditioning        | ContactConditioning (trunkv2.py, ready)        |
| bfactor_module              | BFactorModule (trunkv2.py, ready)              |
| diffusion_conditioning      | DiffusionConditioning (ready)                 |
| structure_module            | AtomDiffusion (ready)                         |
| input_embedder              | InputEmbedder (trunkv2.py, ready)              |
| confidence_module           | ConfidenceModule (confidencev2.py, ready)      |
| template_module             | TODO: needs distributed TemplateModule        |
| affinity_module(s)          | TODO: needs distributed AffinityModule        |
| token_bonds_type            | EmbeddingParamsReplicated (ready)             |
+----------------------------+----------------------------------------------+

When a distributed counterpart becomes available, update the corresponding
``_wrap_*`` method and remove the placeholder.
"""

import gc
import warnings
from copy import deepcopy
from typing import Any, Optional

import numpy as np
import torch
from pytorch_lightning import LightningModule
from torch.distributed.tensor import DTensor, Replicate, Shard
from torch.distributed.tensor import zeros as dtensor_zeros
from torchmetrics import MeanMetric

from boltz.distributed.comm import AttentionPairBiasComm, TransposeComm
from boltz.distributed.data.feature.symmetry import (
    minimum_lddt_symmetry_coords as minimum_lddt_symmetry_coords_dtensor,
)
from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.layers.elementwise_op import ElementwiseOp, elementwise_op, scalar_tensor_op
from boltz.distributed.model.layers.embedding import EmbeddingParamsReplicated
from boltz.distributed.model.layers.flatten_and_unflatten import shardwise_flatten_sharded
from boltz.distributed.model.layers.layernorm import LayerNormParamsReplicated
from boltz.distributed.model.layers.linear import LinearParamsReplicated
from boltz.distributed.model.layers.outer_op import OuterOp, replicate_to_shard_outer_op
from boltz.distributed.model.layers.pairformer import PairformerModule
from boltz.distributed.model.layers.redistribute_transpose import redistribute_transpose
from boltz.distributed.model.layers.repeat_interleave import shardwise_repeat_interleave
from boltz.distributed.model.layers.squeeze import shardwise_squeeze, shardwise_unsqueeze
from boltz.distributed.model.loss.bfactor import bfactor_loss
from boltz.distributed.model.loss.confidencev2 import confidence_loss
from boltz.distributed.model.loss.distogram import distogram_loss
from boltz.distributed.model.modules.confidencev2 import ConfidenceModule
from boltz.distributed.model.modules.diffusion import AtomDiffusion
from boltz.distributed.model.modules.diffusion_conditioning import DiffusionConditioning
from boltz.distributed.model.modules.encoders import RelativePositionEncoder
from boltz.distributed.model.modules.trunkv2 import (
    BFactorModule,
    ContactConditioning,
    DistogramModule,
    InputEmbedder,
    MSAModule,
)
from boltz.distributed.model.optim.ema import DistributedEMA
from boltz.distributed.utils import update_exhaustive_strides
from boltz.model.models.boltz2 import Boltz2 as SerialBoltz2
from boltz.model.optim.scheduler import AlphaFoldLRScheduler


def _ensure_numpy_compatible_dtype(t: torch.Tensor) -> torch.Tensor:
    """Promote tensor dtype to at least float32 for NumPy compatibility.

    NumPy does not support bfloat16/float16.  ``torch.promote_types`` returns
    the wider of *t.dtype* and float32, so half-precision becomes float32 while
    float32/float64 are preserved unchanged.
    """
    return t.to(dtype=torch.promote_types(t.dtype, torch.float32))


def _assert_no_dtensors_in_output(d: dict[str, Any], prefix: str = "") -> None:
    """Raise TypeError if any DTensor values remain in the predict output dict.

    Walks the dict recursively so nested structures like ``pair_chains_iptm``
    are also checked.  Called at the end of ``predict_step`` to ensure the
    writer callback only receives plain tensors.
    """
    for key, val in d.items():
        path = f"{prefix}.{key}" if prefix else key
        if isinstance(val, DTensor):
            raise TypeError(
                f"predict_step output['{path}'] is a DTensor "
                f"(placements={val.placements}). Convert to a plain tensor "
                f"via full_tensor() or to_local() before returning."
            )
        if isinstance(val, dict):
            _assert_no_dtensors_in_output(val, prefix=path)


class _PlaceholderModule(torch.nn.Module):
    """Placeholder for a serial module that does not yet have a distributed implementation.

    Stores the serial module's parameters (so they appear in
    ``state_dict`` for checkpoint compatibility) but raises
    ``NotImplementedError`` on forward.  All parameters are frozen
    because gradients can never flow through a placeholder.
    """

    def __init__(self, serial_module: torch.nn.Module, name: str) -> None:
        super().__init__()
        self._serial = serial_module
        self._name = name
        for p in self._serial.parameters():
            p.requires_grad_(False)

    def forward(self, *args: Any, **kwargs: Any) -> Any:
        raise NotImplementedError(
            f"Distributed version of '{self._name}' is not yet implemented. "
            f"The serial module is stored as a placeholder for parameter "
            f"accounting and checkpoint compatibility."
        )


class Boltz2(LightningModule):
    """Distributed Boltz2 model with context parallelism via DTensor.

    Wraps a fully-initialised serial :class:`~boltz.model.models.boltz2.Boltz2`
    and replaces submodules with DTensor-aware distributed counterparts where
    available.  Submodules without a distributed implementation are wrapped in
    :class:`_PlaceholderModule` so their parameters are tracked in the state
    dict but ``NotImplementedError`` is raised if the code path is hit.
    """

    def __init__(
        self,
        model_serial: SerialBoltz2,
        dist_manager: DistributedManager,
    ) -> None:
        """Initialise the distributed Boltz2 model from a serial model.

        Parameters
        ----------
        model_serial : SerialBoltz2
            Fully initialised serial Boltz2 instance.
        dist_manager : DistributedManager
            Distributed manager with device meshes and process groups.

        Raises
        ------
        Warning
            When a module required for the requested configuration lacks a
            distributed implementation.
        """
        super().__init__()

        self.dist_manager = dist_manager
        self.device_mesh_subgroups = self.dist_manager.device_mesh_subgroups
        self.rank_coords = self.device_mesh_subgroups.get_coordinate()
        self.is_cp_rank_zero = self.rank_coords[1] == 0 and self.rank_coords[2] == 0

        # Preserve hyper-parameters for checkpoint portability
        self.save_hyperparameters(model_serial.hparams)

        self.has_context_parallelism = True

        # ------------------------------------------------------------------ #
        #  Transfer scalar config from serial model                           #
        # ------------------------------------------------------------------ #
        self.training_args = model_serial.training_args
        self.validation_args = model_serial.validation_args
        self.diffusion_loss_args = model_serial.diffusion_loss_args
        self.predict_args = model_serial.predict_args
        self.steering_args = getattr(model_serial, "steering_args", None)
        self.validate_structure = model_serial.validate_structure

        self.no_random_recycling_training = model_serial.no_random_recycling_training
        self.exclude_ions_from_lddt = model_serial.exclude_ions_from_lddt
        self.num_bins = model_serial.num_bins
        self.min_dist = model_serial.min_dist
        self.max_dist = model_serial.max_dist
        self.num_distograms = model_serial.num_distograms
        self.aggregate_distogram = model_serial.aggregate_distogram
        self.is_pairformer_compiled = False
        self.is_msa_compiled = False
        self.is_template_compiled = False
        self.log_loss_every_steps = model_serial.log_loss_every_steps

        self.bond_type_feature = model_serial.bond_type_feature
        self.run_trunk_and_structure = model_serial.run_trunk_and_structure
        self.skip_run_structure = model_serial.skip_run_structure
        self.predict_bfactor = model_serial.predict_bfactor
        self.checkpoint_diffusion_conditioning = model_serial.checkpoint_diffusion_conditioning

        # Confidence / affinity / structure flags
        self.confidence_prediction = model_serial.confidence_prediction
        self.affinity_prediction = False  # TODO: enable when distributed AffinityModule is ready
        self.affinity_ensemble = getattr(model_serial, "affinity_ensemble", False)
        self.affinity_mw_correction = getattr(model_serial, "affinity_mw_correction", True)
        self.token_level_confidence = getattr(model_serial, "token_level_confidence", True)
        self.alpha_pae = model_serial.alpha_pae
        self.structure_prediction_training = model_serial.structure_prediction_training

        if model_serial.affinity_prediction:
            warnings.warn(
                "Affinity prediction is not yet implemented for context parallelism in Boltz2. "
                "Disabling affinity_prediction in the distributed wrapper."
            )

        # EMA – Boltz2 uses callback-based EMA (boltz.model.optim.ema.EMA),
        # and the distributed counterpart is boltz.distributed.model.optim.ema.DistributedEMA.
        # The wrapper itself does *not* manage EMA – the callback does.
        self.use_ema = model_serial.use_ema
        self.ema_decay = model_serial.ema_decay

        # ------------------------------------------------------------------ #
        #  Communication helpers                                              #
        # ------------------------------------------------------------------ #
        layout_group_cp = self.dist_manager.layout_subgroups["cp"]
        self.transpose_comm = TransposeComm(self.dist_manager.group["cp"], layout_group_cp)

        # Process groups for losses that reduce over unique data dimensions
        # only (skipping cp1 which is Replicate for single-repr).
        #
        # We store the individual dp and cp0 groups rather than a single
        # flattened dp×cp0 group.  Creating a flattened group via
        # ``submesh._flatten().get_group()`` triggers ``new_group()`` calls
        # that are NOT coordinated across all ranks when the submesh differs
        # per cp1 position, causing NCCL deadlocks on 4+ GPU topologies.
        # Two sequential all_reduces (dp then cp0) are mathematically
        # equivalent: sum_{dp×cp0}(x) = sum_{cp0}(sum_{dp}(x)).
        self.dp_group = self.device_mesh_subgroups.get_group(0)
        self.cp0_group = self.device_mesh_subgroups.get_group(1)
        self.cp_group = self.dist_manager.group["cp"]

        self._cp1_group = self.device_mesh_subgroups.get_group(2)

        # ------------------------------------------------------------------ #
        # Validation infrastructure                                          #
        # ------------------------------------------------------------------ #
        self.num_val_datasets = model_serial.num_val_datasets
        if self.validate_structure:
            self.val_group_mapper = {}  # maps a dataset index to a validation group name
            self.validator_mapper = {}  # maps a dataset index to a validator
            self.validators = model_serial.validators

        # ------------------------------------------------------------------ #
        #  Wrap submodules that have distributed implementations              #
        # ------------------------------------------------------------------ #
        self.s_init = LinearParamsReplicated(model_serial.s_init, self.device_mesh_subgroups)
        self.z_init_1 = LinearParamsReplicated(model_serial.z_init_1, self.device_mesh_subgroups)
        self.z_init_2 = LinearParamsReplicated(model_serial.z_init_2, self.device_mesh_subgroups)

        self.s_norm = LayerNormParamsReplicated(model_serial.s_norm, self.device_mesh_subgroups)
        self.z_norm = LayerNormParamsReplicated(model_serial.z_norm, self.device_mesh_subgroups)

        self.s_recycle = LinearParamsReplicated(model_serial.s_recycle, self.device_mesh_subgroups)
        self.z_recycle = LinearParamsReplicated(model_serial.z_recycle, self.device_mesh_subgroups)

        self.token_bonds = LinearParamsReplicated(model_serial.token_bonds, self.device_mesh_subgroups)

        # ------------------------------------------------------------------ #
        #  Trunk submodules with distributed implementations                  #
        # ------------------------------------------------------------------ #
        self.msa_module = MSAModule(model_serial.msa_module, self.dist_manager)

        self.pairformer_module = PairformerModule(model_serial.pairformer_module, self.dist_manager)

        self.distogram_module = DistogramModule(
            model_serial.distogram_module,
            self.dist_manager,
            distogram_comm=deepcopy(self.transpose_comm),
        )

        # ------------------------------------------------------------------ #
        #  Submodules that need distributed implementations (placeholders)    #
        # ------------------------------------------------------------------ #
        # InputEmbedder v2: wraps atom-level encoder, feature embedding, etc.
        self.input_embedder = InputEmbedder(
            model_serial.input_embedder,
            device_mesh=self.device_mesh_subgroups,
        )

        # RelativePositionEncoder v2
        self.rel_pos = RelativePositionEncoder(
            model_serial.rel_pos,
            device_mesh=self.dist_manager.device_mesh_subgroups,
            transpose_comm=deepcopy(self.transpose_comm),
        )

        # ContactConditioning
        self.contact_conditioning = ContactConditioning(
            model_serial.contact_conditioning,
            device_mesh=self.dist_manager.device_mesh_subgroups,
        )

        # bond_type_feature embedding (optional)
        if self.bond_type_feature:
            self.token_bonds_type = EmbeddingParamsReplicated(
                model_serial.token_bonds_type, device_mesh=self.dist_manager.device_mesh_subgroups
            )

        # DiffusionConditioning
        self.diffusion_conditioning = DiffusionConditioning(
            model_serial.diffusion_conditioning,
            device_mesh=self.device_mesh_subgroups,
        )

        # AtomDiffusion v2 (structure module)
        ring_comm_diffusion = AttentionPairBiasComm(
            self.dist_manager.group["cp"],
            self.dist_manager.layout_subgroups["cp"],
            self.dist_manager.subgroups["cp"][0],
            self.dist_manager.subgroups["cp"][1],
        )
        self.structure_module = AtomDiffusion(
            model_serial.structure_module,
            device_mesh=self.device_mesh_subgroups,
            ring_comm=ring_comm_diffusion,
            transpose_comm=deepcopy(self.transpose_comm),
        )

        # TemplateModule — preserve weights but skip in forward (not yet distributed)
        if model_serial.use_templates:
            self.template_module = _PlaceholderModule(model_serial.template_module, "TemplateModule")
        # TODO: set self.use_templates = model_serial.use_templates once a real
        # distributed TemplateModule replaces the _PlaceholderModule above.
        self.use_templates = model_serial.use_templates and not isinstance(
            getattr(self, "template_module", None), _PlaceholderModule
        )

        # BFactorModule (optional)
        if self.predict_bfactor:
            self.bfactor_module = BFactorModule(
                model_serial.bfactor_module,
                device_mesh=self.dist_manager.device_mesh_subgroups,
            )

        # ConfidenceModule v2
        if model_serial.confidence_prediction:
            confidence_transpose_comm = deepcopy(self.transpose_comm)
            self.confidence_module = ConfidenceModule(
                model_serial.confidence_module,
                dist_manager=dist_manager,
                transpose_comm=confidence_transpose_comm,
            )

        # AffinityModule(s) (disabled for now)
        if model_serial.affinity_prediction:
            if model_serial.affinity_ensemble:
                self.affinity_module1 = _PlaceholderModule(model_serial.affinity_module1, "AffinityModule1")
                self.affinity_module2 = _PlaceholderModule(model_serial.affinity_module2, "AffinityModule2")
            else:
                self.affinity_module = _PlaceholderModule(model_serial.affinity_module, "AffinityModule")

        # Freeze parameters not involved in structure prediction training
        if not self.structure_prediction_training:
            for name, param in self.named_parameters():
                if (
                    name.split(".")[0] not in ["confidence_module", "affinity_module"]
                    and "out_token_feat_update" not in name
                ):
                    param.requires_grad = False

        # Validate: every trainable parameter must be a DTensor so that
        # on_after_backward can redistribute its gradient.  Non-DTensor
        # trainable params would accumulate plain-tensor gradients that
        # on_after_backward cannot handle.
        for name, param in self.named_parameters():
            if param.requires_grad and not isinstance(param, DTensor):
                raise ValueError(
                    f"Trainable parameter '{name}' is a plain Tensor, not a DTensor. "
                    f"All trainable parameters must be DTensors so that "
                    f"on_after_backward can redistribute their gradients. "
                    f"Either wrap the owning module with a DTensor-aware wrapper "
                    f"or freeze this parameter (requires_grad=False)."
                )

    # ====================================================================== #
    #  Forward pass                                                          #
    # ====================================================================== #

    def forward(
        self,
        feats: dict[str, DTensor],
        recycling_steps: int = 0,
        num_sampling_steps: Optional[int] = None,
        multiplicity_diffusion_train: int = 1,
        diffusion_samples: int = 1,
        max_parallel_samples: Optional[int] = None,
        run_confidence_sequentially: bool = False,
    ) -> dict[str, DTensor]:
        """Forward pass through the distributed Boltz2 model.

        Performs structure prediction using DTensor context parallelism.
        The trunk (input embedding → recycling → MSA → pairformer) is fully
        distributed.  Diffusion conditioning and structure prediction require
        their distributed counterparts to be implemented.

        Parameters
        ----------
        feats : dict[str, DTensor]
            Input features as DTensors.
        recycling_steps : int
            Number of recycling iterations.
        num_sampling_steps : int, optional
            Number of diffusion sampling steps for inference.
        multiplicity_diffusion_train : int
            Training diffusion multiplicity.
        diffusion_samples : int
            Number of diffusion samples for inference.
        max_parallel_samples : int, optional
            Maximum number of parallel diffusion samples.
        run_confidence_sequentially : bool
            Whether to run the confidence module sequentially.

        Returns
        -------
        dict[str, DTensor]
            Dictionary containing model outputs.
        """
        with torch.set_grad_enabled(self.training and self.structure_prediction_training):
            s_inputs = self.input_embedder(feats)

            # Initialise single and pairwise embeddings
            s_init = self.s_init(s_inputs)

            z2 = self.z_init_2(s_inputs)
            z1 = self.z_init_1(s_inputs)

            # Outer sum: globally equivalent to z1[:, :, None, :] + z2[:, None, :, :]
            # Both z1 and z2 have placements (Shard(0), Shard(1), Replicate()).
            # replicate_to_shard_outer_op handles the transpose_then_redistribute
            # for z2 internally and produces (Shard(0), Shard(1), Shard(2)) output.
            # Its backward correctly all-reduces the column/row gradient sums
            # across the Replicate axis.
            z_init = replicate_to_shard_outer_op(
                z1, OuterOp.SUM, axis=1, transpose_comm=self.transpose_comm, input_t=z2
            )

            relative_position_encoding = self.rel_pos(feats)
            z_init = elementwise_op(z_init, relative_position_encoding, ElementwiseOp.SUM)
            z_init = elementwise_op(
                z_init,
                self.token_bonds(feats["token_bonds"].to(dtype=z_init.dtype)),
                ElementwiseOp.SUM,
            )

            if self.bond_type_feature:
                z_init = elementwise_op(
                    z_init,
                    self.token_bonds_type(feats["type_bonds"].long()),
                    ElementwiseOp.SUM,
                )

            z_init = elementwise_op(z_init, self.contact_conditioning(feats), ElementwiseOp.SUM)

            # Initialise recycling buffers using dtensor_zeros to avoid
            # native DTensor dispatch which may trigger implicit communication.
            s: DTensor = dtensor_zeros(
                s_init.shape,
                dtype=s_init.dtype,
                device_mesh=self.device_mesh_subgroups,
                placements=list(s_init.placements),
            )
            z: DTensor = dtensor_zeros(
                z_init.shape,
                dtype=z_init.dtype,
                device_mesh=self.device_mesh_subgroups,
                placements=list(z_init.placements),
            )

            mask = feats["token_pad_mask"].to(dtype=s.dtype)
            pair_mask = feats["token_pair_pad_mask"].to(dtype=z.dtype)

            # Redistribute s_inputs for MSAModule
            # shape: (B, N, D); placements: (S(0), S(1), R) → (S(0), R, S(1))
            s_inputs_redistributed = redistribute_transpose(
                s_inputs,
                transpose_comm=self.transpose_comm,
                output_placements=(Shard(0), Replicate(), Shard(1)),
                dim0=None,
                dim1=None,
            )

            if self.run_trunk_and_structure:
                for i in range(recycling_steps + 1):
                    with torch.set_grad_enabled(
                        self.training and self.structure_prediction_training and (i == recycling_steps)
                    ):
                        if self.training and (i == recycling_steps) and torch.is_autocast_enabled():
                            torch.clear_autocast_cache()

                        # Apply recycling
                        s = elementwise_op(s_init, self.s_recycle(self.s_norm(s)), ElementwiseOp.SUM)
                        z = elementwise_op(z_init, self.z_recycle(self.z_norm(z)), ElementwiseOp.SUM)

                        # Templates (optional)
                        if self.use_templates:
                            # TODO: use distributed template_module when available
                            z = elementwise_op(
                                z,
                                self.template_module(z, feats, pair_mask),
                                ElementwiseOp.SUM,
                            )

                        # MSA module
                        z = elementwise_op(
                            z,
                            self.msa_module(z, s_inputs_redistributed, feats),
                            ElementwiseOp.SUM,
                        )

                        # Pairformer
                        s, z = self.pairformer_module(s, z, mask=mask, pair_mask=pair_mask)

            pdistogram = self.distogram_module(z)
            dict_out: dict[str, DTensor] = {
                "pdistogram": pdistogram,
                "s": s,
                "z": z,
            }

            if self.run_trunk_and_structure and (not self.skip_run_structure):
                # Diffusion conditioning (distributed returns 5 values; to_keys
                # is handled internally by the distributed AtomAttentionEncoder)
                if self.checkpoint_diffusion_conditioning and self.training:
                    q, c, atom_enc_bias, atom_dec_bias, token_trans_bias = torch.utils.checkpoint.checkpoint(
                        self.diffusion_conditioning,
                        s,
                        z,
                        relative_position_encoding,
                        feats,
                        use_reentrant=False,
                    )
                else:
                    q, c, atom_enc_bias, atom_dec_bias, token_trans_bias = self.diffusion_conditioning(
                        s_trunk=s,
                        z_trunk=z,
                        relative_position_encoding=relative_position_encoding,
                        feats=feats,
                    )
                diffusion_conditioning_dict = {
                    "q": q,
                    "c": c,
                    "atom_enc_bias": atom_enc_bias,
                    "atom_dec_bias": atom_dec_bias,
                    "token_trans_bias": token_trans_bias,
                }

                # Inference: reverse diffusion sampling
                if (not self.training) or self.confidence_prediction:
                    with torch.autocast("cuda", enabled=False):
                        compute_dtype = torch.promote_types(s.dtype, torch.float32)
                        struct_out = self.structure_module.sample(
                            s_trunk=s.to(compute_dtype),
                            s_inputs=s_inputs.to(compute_dtype),
                            feats=feats,
                            num_sampling_steps=num_sampling_steps,
                            atom_mask=feats["atom_pad_mask"].to(compute_dtype),
                            multiplicity=diffusion_samples,
                            max_parallel_samples=max_parallel_samples,
                            diffusion_conditioning=diffusion_conditioning_dict,
                        )
                        dict_out.update(struct_out)

                if self.predict_bfactor:
                    dict_out["pbfactor"] = self.bfactor_module(s)

            # Training: diffusion forward pass
            if self.training and self.structure_prediction_training:
                atom_coords = feats["coords"]
                K = atom_coords.shape[1]
                assert K in (multiplicity_diffusion_train, 1)

                # Expand K → multiplicity if needed, then flatten (B, K, L, 3) → (B*K, L, 3).
                if K < multiplicity_diffusion_train:
                    atom_coords = shardwise_repeat_interleave(atom_coords, multiplicity_diffusion_train // K, dim=1)
                feats["coords"] = shardwise_flatten_sharded(atom_coords, start_dim=0, end_dim=1)

                with torch.autocast("cuda", enabled=False):
                    compute_dtype = torch.promote_types(s.dtype, torch.float32)
                    struct_out = self.structure_module(
                        s_trunk=s.to(compute_dtype),
                        s_inputs=s_inputs.to(compute_dtype),
                        feats=feats,
                        multiplicity=multiplicity_diffusion_train,
                        diffusion_conditioning=diffusion_conditioning_dict,  # noqa: F821
                    )
                    dict_out.update(struct_out)

            elif self.training:
                # squeeze(1) removes the singleton ensemble dim:
                # (B, 1, A, 3) → (B*1, A, 3) = (B, A, 3)
                feats["coords"] = shardwise_flatten_sharded(feats["coords"], start_dim=0, end_dim=1)

        if self.confidence_prediction:
            if "frames_idx" in feats and feats["frames_idx"].ndim == 4:
                feats["frames_idx"] = shardwise_flatten_sharded(feats["frames_idx"], start_dim=0, end_dim=1)
            if "frame_resolved_mask" in feats and feats["frame_resolved_mask"].ndim == 3:
                feats["frame_resolved_mask"] = shardwise_flatten_sharded(
                    feats["frame_resolved_mask"], start_dim=0, end_dim=1
                )
            dict_out.update(
                self.confidence_module(
                    s_inputs=s_inputs.detach(),
                    s=s.detach(),
                    z=z.detach(),
                    x_pred=(
                        dict_out["sample_atom_coords"].detach()
                        if not self.skip_run_structure
                        else shardwise_repeat_interleave(feats["coords"], diffusion_samples, dim=0)
                    ),
                    feats=feats,
                    pred_distogram_logits=dict_out["pdistogram"].detach(),
                    multiplicity=diffusion_samples,
                    run_sequentially=run_confidence_sequentially,
                )
            )

        # Affinity (TODO: enable when distributed AffinityModule is ready)

        return dict_out

    # ====================================================================== #
    #  Training step                                                         #
    # ====================================================================== #

    def training_step(self, batch: dict[str, int | DTensor], batch_idx: int) -> DTensor:
        """Training step with distributed loss computation."""
        # Sample recycling steps
        if self.no_random_recycling_training:
            recycling_steps = self.training_args.recycling_steps
        else:
            rgn = np.random.default_rng(self.global_step)
            recycling_steps = rgn.integers(0, self.training_args.recycling_steps + 1).item()

        # Synchronise recycling steps across CP ranks via the flat CP group
        recycling_steps_tensor = torch.tensor(recycling_steps, device=self.device)
        cp_group_global_rank_zero = torch.distributed.get_global_rank(self.cp_group, 0)
        torch.distributed.broadcast(recycling_steps_tensor, src=cp_group_global_rank_zero, group=self.cp_group)
        recycling_steps = recycling_steps_tensor.item()

        if self.training_args.get("sampling_steps_random", None) is not None:
            rgn_sampling = np.random.default_rng(self.global_step)
            sampling_steps = rgn_sampling.choice(self.training_args.sampling_steps_random)
        else:
            sampling_steps = self.training_args.sampling_steps

        # Broadcast sampling_steps across CP ranks for consistency
        sampling_steps_tensor = torch.tensor(int(sampling_steps), device=self.device)
        torch.distributed.broadcast(sampling_steps_tensor, src=cp_group_global_rank_zero, group=self.cp_group)
        sampling_steps = sampling_steps_tensor.item()

        out = self(
            feats=batch,
            recycling_steps=recycling_steps,
            num_sampling_steps=sampling_steps,
            multiplicity_diffusion_train=self.training_args.diffusion_multiplicity,
            diffusion_samples=self.training_args.diffusion_samples,
        )

        # Compute losses
        if self.structure_prediction_training:
            disto_loss = self._compute_distogram_loss(out, batch)
            diffusion_loss_dict = self.structure_module.compute_loss(
                batch,
                out,
                multiplicity=self.training_args.diffusion_multiplicity,
                **self.diffusion_loss_args,
            )
            bfactor_loss_val = self._compute_bfactor_loss(out, batch)
        else:
            zeros_dtensor = dtensor_zeros(
                (),
                requires_grad=False,
                device_mesh=self.device_mesh_subgroups,
                placements=(Replicate(), Replicate(), Replicate()),
            )
            disto_loss = zeros_dtensor
            bfactor_loss_val = zeros_dtensor
            diffusion_loss_dict = {"loss": zeros_dtensor, "loss_breakdown": {}}

        confidence_loss_dict = self._compute_confidence_loss(out, batch)

        # Aggregate losses
        loss = elementwise_op(
            elementwise_op(
                elementwise_op(
                    scalar_tensor_op(
                        self.training_args.diffusion_loss_weight,
                        diffusion_loss_dict["loss"],
                        ElementwiseOp.PROD,
                    ),
                    scalar_tensor_op(
                        self.training_args.distogram_loss_weight,
                        disto_loss,
                        ElementwiseOp.PROD,
                    ),
                    ElementwiseOp.SUM,
                ),
                scalar_tensor_op(
                    # Default 0.0 matches serial boltz2.py (added after other weights).
                    self.training_args.get("bfactor_loss_weight", 0.0),
                    bfactor_loss_val,
                    ElementwiseOp.PROD,
                ),
                ElementwiseOp.SUM,
            ),
            scalar_tensor_op(
                self.training_args.confidence_loss_weight,
                confidence_loss_dict["loss"],
                ElementwiseOp.PROD,
            ),
            ElementwiseOp.SUM,
        )

        if not (self.global_step % self.log_loss_every_steps):
            self.log("train/distogram_loss", disto_loss.to_local())
            self.log("train/diffusion_loss", diffusion_loss_dict["loss"].to_local())
            for k, v in diffusion_loss_dict["loss_breakdown"].items():
                self.log(f"train/{k}", v.to_local() if isinstance(v, DTensor) else v)
            if self.confidence_prediction:
                self.log("train/confidence_loss", confidence_loss_dict["loss"].to_local())
                for k, v in confidence_loss_dict["loss_breakdown"].items():
                    self.log(f"train/{k}", v.to_local() if isinstance(v, DTensor) else v)
            self.log("train/loss", loss.to_local())
            self.training_log()

        return loss

    def _compute_distogram_loss(self, out: dict, batch: dict) -> DTensor:
        """Compute distogram loss using the distributed implementation."""
        disto_loss, _ = distogram_loss(
            out, batch, comm=self.transpose_comm, aggregate_distogram=self.aggregate_distogram
        )
        return disto_loss

    def _compute_confidence_loss(self, out: dict, batch: dict) -> dict:
        """Compute confidence loss if confidence_prediction is enabled."""
        if not self.confidence_prediction:
            zeros = dtensor_zeros(
                (),
                requires_grad=False,
                device_mesh=self.device_mesh_subgroups,
                placements=(Replicate(), Replicate(), Replicate()),
            )
            return {"loss": zeros, "loss_breakdown": {}}

        return_dict = self.get_true_coordinates(
            batch,
            out,
            diffusion_samples=self.training_args.diffusion_samples,
            symmetry_correction=self.training_args.symmetry_correction,
        )

        true_coords = return_dict["true_coords"]
        true_coords_resolved_mask = return_dict["true_coords_resolved_mask"]

        return confidence_loss(
            out,
            batch,
            true_coords,
            true_coords_resolved_mask,
            comm=self.transpose_comm,
            token_level_confidence=self.token_level_confidence,
            alpha_pae=self.alpha_pae,
            multiplicity=self.training_args.diffusion_samples,
            dist_manager=self.dist_manager,
            group_layout=self.dist_manager.layout_subgroups["cp"],
        )

    def _compute_bfactor_loss(self, out: dict, batch: dict) -> DTensor:
        """Compute bfactor loss if enabled."""
        if self.predict_bfactor:
            return bfactor_loss(
                out,
                batch,
                device_mesh=self.device_mesh_subgroups,
                dp_group=self.dp_group,
                cp0_group=self.cp0_group,
                cp1_group=self._cp1_group,
            )
        return dtensor_zeros(
            (),
            requires_grad=False,
            device_mesh=self.device_mesh_subgroups,
            placements=(Replicate(), Replicate(), Replicate()),
        )

    # ====================================================================== #
    #  Logging helpers                                                       #
    # ====================================================================== #

    def training_log(self) -> None:
        """Log training metrics."""
        self.log("train/param_norm", self.parameter_norm(self), prog_bar=False)

        lr = self.trainer.optimizers[0].param_groups[0]["lr"]
        self.log("lr", lr, prog_bar=False)

        self.log("train/param_norm_msa_module", self.parameter_norm(self.msa_module), prog_bar=False)
        self.log("train/param_norm_pairformer_module", self.parameter_norm(self.pairformer_module), prog_bar=False)
        self.log("train/param_norm_structure_module", self.parameter_norm(self.structure_module), prog_bar=False)

        if self.confidence_prediction:
            self.log(
                "train/param_norm_confidence_module",
                self.parameter_norm(self.confidence_module),
                prog_bar=False,
            )

    def gradient_norm(self, module: torch.nn.Module) -> torch.Tensor:
        """Compute L2 norm of gradients for a distributed module.

        Handles both DTensor and plain Tensor parameters (e.g. from
        placeholder modules whose serial weights have not yet been
        distributed).
        """
        norms_sq: list[torch.Tensor] = []
        for p in module.parameters():
            if p.requires_grad and p.grad is not None:
                grad = p.grad.to_local() if isinstance(p.grad, DTensor) else p.grad
                norms_sq.append(grad.norm(p=2) ** 2)
        if len(norms_sq) == 0:
            return torch.tensor(0.0, device=self.dist_manager.device.type)
        return torch.stack(norms_sq).sum().sqrt().to(device=self.dist_manager.device.type)

    def parameter_norm(self, module: torch.nn.Module) -> torch.Tensor:
        """Compute L2 norm of parameters for a distributed module.

        Handles both DTensor and plain Tensor parameters.
        """
        norms_sq: list[torch.Tensor] = []
        for p in module.parameters():
            if p.requires_grad:
                val = p.to_local() if isinstance(p, DTensor) else p
                norms_sq.append(val.norm(p=2) ** 2)
        if len(norms_sq) == 0:
            return torch.tensor(0.0, device=self.dist_manager.device.type)
        return torch.stack(norms_sq).sum().sqrt().to(device=self.dist_manager.device.type)

    # ====================================================================== #
    #  Gradient redistribution                                               #
    # ====================================================================== #

    def on_after_backward(self) -> None:
        """Redistribute DTensor gradients to Replicate placement after backward.

        Called after ``loss.backward()`` and before ``optimizer.step()``.
        Ensures gradients are properly synchronised across context parallel
        ranks via all-reduce.

        The ``__init__`` validates that all trainable parameters are DTensors,
        so only DTensor or None gradients should appear here.  Plain-tensor
        gradients are skipped with a warning as a defensive measure.

        Note: parameter gradients on the Replicate (cp1) axis are already
        synchronised per-layer via ``Partial("avg")`` in the linear/layernorm
        backward (``avg_over_replicate_param_grad=True``).  The redistribute
        call below handles any remaining ``Partial`` placements from other
        sources.
        """
        for name, p in self.named_parameters():
            if p.grad is None:
                continue
            if isinstance(p.grad, DTensor):
                p.grad = p.grad.redistribute(p.grad.device_mesh, [Replicate()] * p.grad.device_mesh.ndim)
            else:
                warnings.warn(
                    f"Parameter '{name}' has a plain-tensor gradient (type={type(p.grad)}), "
                    f"skipping redistribution. This should not happen — all trainable "
                    f"parameters should be DTensors. Check __init__ validation.",
                    stacklevel=2,
                )

        if not (self.global_step % self.log_loss_every_steps):
            self.log("train/grad_norm", self.gradient_norm(self), prog_bar=False)
            self.log("train/grad_norm_msa_module", self.gradient_norm(self.msa_module), prog_bar=False)
            self.log("train/grad_norm_pairformer_module", self.gradient_norm(self.pairformer_module), prog_bar=False)
            self.log("train/grad_norm_structure_module", self.gradient_norm(self.structure_module), prog_bar=False)
            if self.confidence_prediction:
                self.log(
                    "train/grad_norm_confidence_module",
                    self.gradient_norm(self.confidence_module),
                    prog_bar=False,
                )

    # ====================================================================== #
    #  Epoch-level hooks                                                     #
    # ====================================================================== #

    def on_train_epoch_end(self) -> None:
        """Log epoch-level training metrics."""

    # ====================================================================== #
    #  Validation                                                            #
    # ====================================================================== #

    def validation_step(self, batch: dict[str, Any], batch_idx: int) -> None:
        """Validation step delegating to the distributed validator."""
        if self.validate_structure:
            try:
                # idx_dataset is a non-sharded feature; CollateDTensor collates
                # TRAINING_METADATA_FEATURES as a Python list of tensors.

                msg = "Only batch=1 is supported for validation"
                assert len(batch["idx_dataset"]) == 1, msg
                assert batch["idx_dataset"][0].shape[0] == 1, msg

                idx_dataset = batch["idx_dataset"][0].item()
                validator = self.validator_mapper[idx_dataset]

                out = validator.run_model(model=self, batch=batch, idx_dataset=idx_dataset)
                validator.process(
                    model=self,
                    batch=batch,
                    out=out,
                    idx_dataset=idx_dataset,
                    transpose_comm=self.transpose_comm,
                )
            except RuntimeError as e:
                idx_dataset = batch["idx_dataset"][0].item()
                if "out of memory" in str(e):
                    msg = f"| WARNING: ran out of memory, skipping batch, {idx_dataset}"
                    print(msg)
                    torch.cuda.empty_cache()
                    gc.collect()
                    return
                raise e
        else:
            try:
                out = self(
                    batch,
                    recycling_steps=self.validation_args.recycling_steps,
                    num_sampling_steps=self.validation_args.sampling_steps,
                    diffusion_samples=self.validation_args.diffusion_samples,
                    run_confidence_sequentially=self.validation_args.get("run_confidence_sequentially", False),
                )
            except RuntimeError as e:
                idx_dataset = batch["idx_dataset"][0].item()
                if "out of memory" in str(e):
                    msg = f"| WARNING: ran out of memory, skipping batch, {idx_dataset}"
                    print(msg)
                    torch.cuda.empty_cache()
                    gc.collect()
                    return
                raise e

    def on_validation_epoch_end(self) -> None:
        """Aggregate all metrics for each validator."""
        if not self.validate_structure:
            return

        if self.trainer.sanity_checking:
            for validator in self.validator_mapper.values():
                for m in validator.modules():
                    if isinstance(m, MeanMetric):
                        m.reset()
            return

        for validator in self.validator_mapper.values():
            validator.on_epoch_end(model=self)

    def setup(self, stage: str) -> None:
        """Set the model for training, validation and inference."""

        if (
            stage != "predict"
            and hasattr(self.trainer, "datamodule")
            and self.trainer.datamodule
            and self.validate_structure
        ):
            self.val_group_mapper.update(self.trainer.datamodule.val_group_mapper)

            l1 = len(self.val_group_mapper)
            l2 = self.num_val_datasets
            msg = (
                f"Number of validation datasets num_val_datasets={l2} "
                f"does not match the number of val_group_mapper entries={l1}."
            )
            assert l1 == l2, msg

            all_validator_names = []
            for validator in self.validators:
                for val_name in validator.val_names:
                    msg = f"Validator {val_name} duplicated in validators."
                    assert val_name not in all_validator_names, msg
                    all_validator_names.append(val_name)
                    for val_idx, val_group in self.val_group_mapper.items():
                        if val_name == val_group["label"]:
                            self.validator_mapper[val_idx] = validator

            msg = "Mismatch between validator names and val_group_mapper values."
            assert set(all_validator_names) == {x["label"] for x in self.val_group_mapper.values()}, msg

    def get_true_coordinates(
        self,
        batch: dict[str, Any],
        out: dict[str, Any],
        diffusion_samples: int,
        symmetry_correction: bool,
        expand_to_diffusion_samples: bool = True,
    ) -> dict[str, Any]:
        """Compute true coordinates for validation/confidence loss.

        In the distributed case, coordinates are DTensors sharded across the
        (DP, CP_0, CP_1) mesh.  When ``symmetry_correction`` is True, each
        sample is processed by :func:`minimum_lddt_symmetry_coords_dtensor`
        which internally gathers coordinates along CP axes, runs the symmetry
        search on plain tensors using the ``cdist_lddt`` triton kernel, then
        re-shards the results.

        Parameters
        ----------
        batch : dict[str, Any]
            Input features as DTensors.
        out : dict[str, Any]
            Model outputs including ``sample_atom_coords`` as a DTensor.
        diffusion_samples : int
            Number of diffusion samples per batch element.
        symmetry_correction : bool
            Whether to apply symmetry correction.
        expand_to_diffusion_samples : bool
            Whether to expand coordinates to match diffusion samples.

        Returns
        -------
        dict[str, Any]
            Dictionary with ``true_coords``, ``true_coords_resolved_mask``,
            ``rmsds``, and ``best_rmsd_recall``.

        """
        if symmetry_correction:
            assert expand_to_diffusion_samples, "expand_to_diffusion_samples must be true for symmetry correction."

        return_dict: dict[str, Any] = {}
        sample_atom_coords = out["sample_atom_coords"]

        if symmetry_correction:
            true_coords_list: list[DTensor] = []
            true_mask_list: list[DTensor] = []

            local_batch_size = batch["token_index"].to_local().shape[0]
            for i_batch_local in range(local_batch_size):
                for rep in range(diffusion_samples):
                    i_local = i_batch_local * diffusion_samples + rep
                    best_true_coords, best_true_mask = minimum_lddt_symmetry_coords_dtensor(
                        coords=sample_atom_coords,
                        feats=batch,
                        index_batch_local=i_batch_local,
                        i_batch_multiplicity_local=i_local,
                    )
                    true_coords_list.append(best_true_coords)
                    true_mask_list.append(best_true_mask)

            assert len(true_coords_list) >= 1, "There should be at least 1 true coords processed"
            _coords_local = torch.cat([c.to_local() for c in true_coords_list], dim=0)
            _coords_global_shape = list(_coords_local.shape)
            for _pi, _pl in enumerate(true_coords_list[0].placements):
                if hasattr(_pl, "dim"):
                    _coords_global_shape[_pl.dim] *= true_coords_list[0].device_mesh.size(_pi)
            _coords_global_shape = torch.Size(_coords_global_shape)
            true_coords = DTensor.from_local(
                _coords_local,
                device_mesh=true_coords_list[0].device_mesh,
                placements=true_coords_list[0].placements,
                shape=_coords_global_shape,
                stride=update_exhaustive_strides(_coords_local.shape, _coords_local.stride(), _coords_global_shape),
            )
            _mask_local = torch.cat([m.to_local() for m in true_mask_list], dim=0)
            _mask_global_shape = list(_mask_local.shape)
            for _pi, _pl in enumerate(true_mask_list[0].placements):
                if hasattr(_pl, "dim"):
                    _mask_global_shape[_pl.dim] *= true_mask_list[0].device_mesh.size(_pi)
            _mask_global_shape = torch.Size(_mask_global_shape)
            true_coords_resolved_mask = DTensor.from_local(
                _mask_local,
                device_mesh=true_mask_list[0].device_mesh,
                placements=true_mask_list[0].placements,
                shape=_mask_global_shape,
                stride=update_exhaustive_strides(_mask_local.shape, _mask_local.stride(), _mask_global_shape),
            )

            true_coords = shardwise_unsqueeze(true_coords, dim=1)
            return_dict["true_coords"] = true_coords
            return_dict["true_coords_resolved_mask"] = true_coords_resolved_mask
            return_dict["rmsds"] = 0
            return_dict["best_rmsd_recall"] = 0
        else:
            true_coords_resolved_mask = batch["atom_resolved_mask"]
            true_coords = shardwise_squeeze(batch["coords"], dim=1)
            if expand_to_diffusion_samples:
                true_coords = shardwise_repeat_interleave(true_coords, diffusion_samples, 0)
                true_coords_resolved_mask = shardwise_repeat_interleave(true_coords_resolved_mask, diffusion_samples, 0)

            return_dict["true_coords"] = true_coords
            return_dict["true_coords_resolved_mask"] = true_coords_resolved_mask
            return_dict["rmsds"] = 0
            return_dict["best_rmsd_recall"] = 0
            return_dict["best_rmsd_precision"] = 0

        return return_dict

    # ====================================================================== #
    #  Prediction                                                            #
    # ====================================================================== #

    def predict_step(
        self, batch: dict[str, DTensor], batch_idx: int, dataloader_idx: int = 0
    ) -> dict[str, torch.Tensor]:
        """Prediction step with distributed inference.

        Parameters
        ----------
        batch : dict[str, DTensor]
            Input features as DTensors.
        batch_idx : int
            Index of the current batch.
        dataloader_idx : int
            Index of the current dataloader.

        Returns
        -------
        dict[str, torch.Tensor]
            Prediction results gathered on rank 0 of each CP column.
        """
        try:
            out = self(
                batch,
                recycling_steps=self.predict_args["recycling_steps"],
                num_sampling_steps=self.predict_args["sampling_steps"],
                diffusion_samples=self.predict_args["diffusion_samples"],
                max_parallel_samples=self.predict_args.get("max_parallel_samples", None),
                run_confidence_sequentially=True,
            )
            pred_dict: dict[str, Any] = {"exception": False}

            # Gather coords and masks onto rank 0 of each CP column
            tag_group_gather = 0
            ranks_gather = self.dist_manager.subgroups_ranks["cp"][tag_group_gather]
            group_gather = self.dist_manager.subgroups["cp"][tag_group_gather]
            world_size_gather = len(ranks_gather)

            coords = out["sample_atom_coords"].to_local()
            mask = batch["atom_pad_mask"].to_local()

            if self.dist_manager.subgroups_rank["cp"][tag_group_gather] == 0:
                gather_list_coords = [torch.empty_like(coords) for _ in range(world_size_gather)]
                gather_list_mask = [torch.empty_like(mask) for _ in range(world_size_gather)]
            else:
                gather_list_coords = None
                gather_list_mask = None

            torch.distributed.gather(coords, gather_list_coords, dst=ranks_gather[0], group=group_gather)
            torch.distributed.gather(mask, gather_list_mask, dst=ranks_gather[0], group=group_gather)

            if self.dist_manager.subgroups_rank["cp"][tag_group_gather] == 0:
                pred_dict["masks"] = torch.concat(gather_list_mask, dim=1)
                pred_dict["coords"] = _ensure_numpy_compatible_dtype(torch.concat(gather_list_coords, dim=1))
            else:
                pred_dict["masks"] = mask
                pred_dict["coords"] = _ensure_numpy_compatible_dtype(coords)

            if self.confidence_prediction:
                # Per-token and per-pair confidence outputs are sharded across
                # CP ranks. Use full_tensor() to reassemble the global tensor
                # for the writer (predict path only, not hot).
                for key in ["pde", "plddt"]:
                    val = out[key]
                    val = val.full_tensor() if isinstance(val, DTensor) else val
                    pred_dict[key] = _ensure_numpy_compatible_dtype(val)

                # Scalar metrics — already (Shard(0), Replicate(), Replicate()),
                # so to_local() gives the correct value on every rank.
                for key in [
                    "complex_plddt",
                    "complex_iplddt",
                    "complex_pde",
                    "complex_ipde",
                ]:
                    if key in out:
                        val = out[key]
                        pred_dict[key] = val.to_local() if isinstance(val, DTensor) else val

                # Confidence score (matching serial formula)
                cplddt = pred_dict["complex_plddt"]
                if "iptm" in out:
                    iptm_val = out["iptm"].to_local() if isinstance(out["iptm"], DTensor) else out["iptm"]
                    ptm_val = out["ptm"].to_local() if isinstance(out["ptm"], DTensor) else out["ptm"]
                    use_iptm = not torch.allclose(iptm_val, torch.zeros_like(iptm_val))
                    pred_dict["confidence_score"] = (4 * cplddt + (iptm_val if use_iptm else ptm_val)) / 5
                else:
                    pred_dict["confidence_score"] = cplddt

                if self.alpha_pae > 0 and "pae" in out:
                    # pae is sharded across CP — needs full_tensor() to reassemble.
                    val = out["pae"]
                    val = val.full_tensor() if isinstance(val, DTensor) else val
                    pred_dict["pae"] = _ensure_numpy_compatible_dtype(val)

                    # ptm, iptm, *_iptm are globally reduced scalars with
                    # placements (Shard(0), Replicate(), Replicate()) — already
                    # fully reduced across CP.  to_local() extracts the value
                    # with no communication.
                    for key in ["ptm", "iptm", "ligand_iptm", "protein_iptm"]:
                        if key in out:
                            val = out[key]
                            pred_dict[key] = val.to_local() if isinstance(val, DTensor) else val

                    if "pair_chains_iptm" in out:
                        pci = out["pair_chains_iptm"]
                        if isinstance(pci, dict):
                            # Success path: nested dict of DTensors →
                            # plain tensors.  Values are globally reduced
                            # (Shard(0), Replicate(), Replicate()).
                            pred_dict["pair_chains_iptm"] = {
                                k1: {k2: v.to_local() if isinstance(v, DTensor) else v for k2, v in inner.items()}
                                for k1, inner in pci.items()
                            }
                        else:
                            # Fallback: compute_ptms failed and the
                            # confidence module returned a zero tensor.
                            # The writer iterates pair_chains_iptm as a
                            # nested dict — pass an empty dict so the
                            # comprehension produces nothing.
                            pred_dict["pair_chains_iptm"] = {}

            _assert_no_dtensors_in_output(pred_dict)
            return pred_dict

        except RuntimeError as e:
            if "out of memory" in str(e):
                print("| WARNING: ran out of memory, skipping batch")
                torch.cuda.empty_cache()
                gc.collect()
                return {"exception": True}
            raise

    # ====================================================================== #
    #  Optimiser                                                             #
    # ====================================================================== #

    def configure_optimizers(self) -> torch.optim.Optimizer:
        """Configure the optimizer following the serial Boltz2 pattern."""
        param_dict = dict(self.named_parameters())

        if self.structure_prediction_training:
            all_parameter_names = [pn for pn, p in self.named_parameters() if p.requires_grad]
        else:
            all_parameter_names = [
                pn
                for pn, p in self.named_parameters()
                if p.requires_grad and ("out_token_feat_update" in pn or "confidence_module" in pn)
            ]

        weight_decay = self.training_args.get("weight_decay", 0.0)
        if weight_decay > 0 and self.training_args.get("weight_decay_exclude", False):
            nodecay_params_names = [
                pn
                for pn in all_parameter_names
                if (
                    "norm" in pn
                    or "rel_pos" in pn
                    or ".s_init" in pn
                    or ".z_init_" in pn
                    or "token_bonds" in pn
                    or "embed_atom_features" in pn
                    or "dist_bin_pairwise_embed" in pn
                )
            ]
            nodecay_params = [param_dict[pn] for pn in nodecay_params_names]
            decay_params = [param_dict[pn] for pn in all_parameter_names if pn not in nodecay_params_names]
            optim_groups = [
                {"params": decay_params, "weight_decay": weight_decay},
                {"params": nodecay_params, "weight_decay": 0.0},
            ]
            optimizer = torch.optim.AdamW(
                optim_groups,
                betas=(self.training_args.adam_beta_1, self.training_args.adam_beta_2),
                eps=self.training_args.adam_eps,
                lr=self.training_args.base_lr,
            )
        else:
            optimizer = torch.optim.AdamW(
                [param_dict[pn] for pn in all_parameter_names],
                betas=(self.training_args.adam_beta_1, self.training_args.adam_beta_2),
                eps=self.training_args.adam_eps,
                lr=self.training_args.base_lr,
                weight_decay=weight_decay,
            )

        if self.training_args.lr_scheduler == "af3":
            scheduler = AlphaFoldLRScheduler(
                optimizer,
                base_lr=self.training_args.base_lr,
                max_lr=self.training_args.max_lr,
                warmup_no_steps=self.training_args.lr_warmup_no_steps,
                start_decay_after_n_steps=self.training_args.lr_start_decay_after_n_steps,
                decay_every_n_steps=self.training_args.lr_decay_every_n_steps,
                decay_factor=self.training_args.lr_decay_factor,
            )
            return [optimizer], [{"scheduler": scheduler, "interval": "step"}]

        return optimizer

    # ====================================================================== #
    #  EMA – Boltz2 uses callback-based EMA                                  #
    # ====================================================================== #

    def configure_callbacks(self) -> list:
        """Configure model callbacks.

        When EMA is enabled, returns a :class:`DistributedEMA` callback which
        handles DTensor ↔ plain-tensor conversions automatically.
        """
        if self.use_ema:
            return [DistributedEMA(self.ema_decay)]
        return []

    # ====================================================================== #
    #  Checkpoint hooks                                                      #
    # ====================================================================== #

    def on_load_checkpoint(self, checkpoint: dict[str, Any]) -> None:
        """Adjust checkpoint hyperparameters on load (matching serial Boltz2)."""
        lr = self.training_args.max_lr
        weight_decay = self.training_args.weight_decay
        if "optimizer_states" in checkpoint:
            for state in checkpoint["optimizer_states"]:
                for group in state["param_groups"]:
                    group["lr"] = lr
                    group["weight_decay"] = weight_decay
        if "lr_schedulers" in checkpoint:
            for scheduler in checkpoint["lr_schedulers"]:
                scheduler["max_lr"] = lr
                scheduler["base_lrs"] = [lr] * len(scheduler["base_lrs"])
                scheduler["_last_lr"] = [lr] * len(scheduler["_last_lr"])
        if "hyper_parameters" in checkpoint:
            checkpoint["hyper_parameters"]["training_args"]["max_lr"] = lr
            checkpoint["hyper_parameters"]["training_args"]["diffusion_multiplicity"] = (
                self.training_args.diffusion_multiplicity
            )
            checkpoint["hyper_parameters"]["training_args"]["recycling_steps"] = self.training_args.recycling_steps
            checkpoint["hyper_parameters"]["training_args"]["weight_decay"] = self.training_args.weight_decay
