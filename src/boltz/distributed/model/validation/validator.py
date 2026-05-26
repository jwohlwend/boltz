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


"""Distributed validator for Boltz-2 with context parallelism.

Extends the serial :class:`~boltz.model.validation.validator.Validator` to
handle DTensor inputs.  Each metric computation method is overridden to
selectively all-gather only the features it needs, then delegate to the
serial/triton implementation on plain tensors.

Architecture:
  DistributedValidator(Validator)
    - Overrides: run_model, compute_disto_loss, compute_disto_lddt,
      get_lddt_metrics, get_clash_metrics, get_pb_metrics,
      get_confidence_metrics, common_val_step, common_on_epoch_end
    - Each metric override: gather DTensor features → call serial function →
      return plain-tensor results
    - Metric storage, update, and epoch-end aggregation remain in the serial
      base class operating on plain tensors.
"""

from collections import defaultdict
from typing import Optional

import torch
from pytorch_lightning import LightningModule
from torch import Tensor, nn
from torch.distributed.tensor import DTensor
from torchmetrics import MeanMetric

from boltz.distributed.comm import TransposeComm
from boltz.distributed.model.layers.atom_to_token import (
    reconstruct_atom_to_token_global,
    reconstruct_r_set_to_rep_atom_global,
    reconstruct_token_to_rep_atom_global,
    single_repr_rep_atom_to_token,
)
from boltz.distributed.model.loss.distogram import distogram_loss
from boltz.distributed.model.loss.validation import (
    clash_score,
    compute_disto_lddt,
    compute_plddt_mae_triton,
    get_lddt_metrics,
)
from boltz.distributed.model.validation.utils import gather_along_cp
from boltz.model.loss.validation import (
    compute_pae_mae,
    compute_pde_mae,
    weighted_minimum_rmsd_single,
)
from boltz.model.validation.validator import Validator


class DistributedValidator(Validator):
    """Distributed validator that handles DTensor inputs for Boltz-2 CP.

    Overrides metric computation methods from :class:`Validator` to
    selectively all-gather DTensor features before calling the serial
    implementations.  Metric storage and epoch-end aggregation remain in
    the serial base class.
    """

    def __init__(
        self,
        val_names: list[str],
        confidence_prediction: bool = False,
        physicalism_metrics: bool = False,
        rmsd_metrics: bool = False,
        clash_score_metrics: bool = False,
        override_val_method: Optional[str] = None,
    ) -> None:
        """Initialize the distributed validator.

        Parameters
        ----------
        val_names : list[str]
            The list of validation names.
        confidence_prediction : bool
            Whether to predict confidence.
        physicalism_metrics : bool
            Whether to compute physicalism metrics.
        rmsd_metrics : bool
            Whether to compute rmsd metrics.
        clash_score_metrics : bool
            Whether to compute clash score metrics.
        override_val_method : Optional[str]
            The override validation method.
        """
        super().__init__(
            val_names=val_names,
            confidence_prediction=confidence_prediction,
            physicalism_metrics=physicalism_metrics,
            override_val_method=override_val_method,
        )
        self.rmsd_metrics = rmsd_metrics
        self.clash_score_metrics = clash_score_metrics

        if rmsd_metrics:
            self.folding_metrics["rmsd"] = nn.ModuleList([nn.ModuleDict() for _ in range(self.num_val_datasets)])
            for val_idx in range(self.num_val_datasets):
                self.folding_metrics["rmsd"][val_idx]["rmsd"] = MeanMetric()

        if clash_score_metrics:
            self.folding_metrics["clash_score"] = nn.ModuleList([nn.ModuleDict() for _ in range(self.num_val_datasets)])
            for val_idx in range(self.num_val_datasets):
                self.folding_metrics["clash_score"][val_idx]["clash_atoms_count"] = MeanMetric()
                self.folding_metrics["clash_score"][val_idx]["clash_atoms_fraction"] = MeanMetric()

        # In our CP code, lightning is blind to our distributed computation context
        # so MeanMetric is strictly single device and we should not rely on
        # MeanMetric.compute() to get inter-DP all_reduce mean.
        for m in self.modules():
            if isinstance(m, MeanMetric):
                m.sync_on_compute = False
                m._to_sync = False

    def run_model(
        self,
        model: LightningModule,
        batch: dict[str, DTensor],
        idx_dataset: int,
    ) -> dict[str, DTensor]:
        """Compute the forward pass using the distributed model.

        Parameters
        ----------
        model : LightningModule
            The distributed Boltz2 LightningModule.
        batch : dict[str, DTensor]
            Batch features as DTensors.
        idx_dataset : int
            Dataset index.

        Returns
        -------
        dict[str, DTensor]
            Model outputs as DTensors.
        """
        if self.override_val_method is not None:
            raise NotImplementedError("Override validation method is not supported for distributed validation")
            # from boltz.distributed.model.layers.elementwise_op import ElementwiseOp, scalar_tensor_op
            # new_feature = scalar_tensor_op(0.0, batch["method_feature"], ElementwiseOp.PROD)
            # new_feature = scalar_tensor_op(self.override_val_method, new_feature, ElementwiseOp.SUM)
            # batch["method_feature"] = new_feature

        out = model(
            batch,
            recycling_steps=model.validation_args.recycling_steps,
            num_sampling_steps=model.validation_args.sampling_steps,
            diffusion_samples=model.validation_args.diffusion_samples,
            run_confidence_sequentially=model.validation_args.get("run_confidence_sequentially", False),
        )
        return out

    def compute_disto_loss(
        self,
        model: LightningModule,
        out: dict[str, DTensor],
        batch: dict[str, DTensor],
        idx_dataset: int,
        transpose_comm: TransposeComm,
    ) -> Tensor:
        """Compute distogram loss using DTensor-native implementation.

        Uses the existing distributed distogram loss which operates directly
        on DTensors without requiring all-gather.
        """
        val_disto_loss, _ = distogram_loss(
            out, batch, comm=transpose_comm, aggregate_distogram=model.aggregate_distogram
        )
        return val_disto_loss.to_local()

    def compute_disto_lddt(
        self,
        model: LightningModule,
        batch: dict[str, Tensor],
        out: dict[str, Tensor],
        idx_dataset: int,
    ) -> tuple[dict, dict]:
        """Compute distogram lddt."""
        disto_lddt_dict, disto_total_dict = compute_disto_lddt(model, batch, out)
        return disto_lddt_dict, disto_total_dict

    def get_lddt_metrics(
        self,
        *args,
        **kwargs,
    ) -> None:
        """Override the serial get_lddt_metrics method."""
        raise NotImplementedError("DistributedValidator does not need to implement get_lddt_metrics")

    def get_clash_metrics(
        self,
        batch: dict[str, DTensor],
        out: dict[str, DTensor],
        batch_gathered: dict[str, Tensor],
        out_gathered: dict[str, Tensor],
    ) -> tuple[dict, dict]:
        """Compute clash metrics by gathering features at global atom dimension.

        Reuses unstripped ``asym_id``, ``atom_pad_mask``, and
        ``sample_atom_coords`` from ``batch_gathered``/``out_gathered``.
        Gathers ``atom_to_token`` and ``ref_element`` on demand.
        Non-sharded features pass through directly from ``batch``.
        """
        clash_feats = {
            "asym_id": batch_gathered["asym_id"],
            "atom_to_token": reconstruct_atom_to_token_global(batch["atom_to_token"]),
            "ref_element": gather_along_cp(batch["ref_element"]),
            "atom_pad_mask": batch_gathered["atom_pad_mask"],
            "connections_edge_index": batch["connections_edge_index"],
            "chain_symmetries": batch["chain_symmetries"],
        }
        clash_out = {"sample_atom_coords": out_gathered["sample_atom_coords"]}

        result = super().get_clash_metrics(clash_feats, clash_out)
        del clash_feats, clash_out
        return result

    def get_pb_metrics(
        self,
        batch: dict[str, DTensor],
        out: dict[str, DTensor],
        batch_gathered: dict[str, Tensor],
        out_gathered: dict[str, Tensor],
    ) -> tuple[dict, dict]:
        """Compute PB metrics by gathering features at global atom dimension.

        Reuses unstripped ``asym_id``, ``mol_type``, and
        ``sample_atom_coords`` from ``batch_gathered``/``out_gathered``.
        Gathers ``atom_to_token`` on demand.  Ligand features are not
        sharded and pass through directly from ``batch``.
        """
        pb_feats = {
            "asym_id": batch_gathered["asym_id"],
            "atom_to_token": reconstruct_atom_to_token_global(batch["atom_to_token"]),
            "mol_type": batch_gathered["mol_type"],
        }
        _LIGAND_KEYS = (
            "ligand_edge_index",
            "ligand_edge_lower_bounds",
            "ligand_edge_upper_bounds",
            "ligand_edge_bond_mask",
            "ligand_edge_angle_mask",
            "ligand_chiral_atom_index",
            "ligand_chiral_check_mask",
            "ligand_chiral_atom_orientations",
            "ligand_stereo_bond_index",
            "ligand_stereo_check_mask",
            "ligand_stereo_bond_orientations",
            "ligand_aromatic_5_ring_index",
            "ligand_aromatic_6_ring_index",
            "ligand_planar_double_bond_index",
        )
        for k in _LIGAND_KEYS:
            if k in batch:
                pb_feats[k] = batch[k]
        pb_out = {"sample_atom_coords": out_gathered["sample_atom_coords"]}

        result = super().get_pb_metrics(pb_feats, pb_out)
        del pb_feats, pb_out
        return result

    def get_confidence_metrics(
        self,
        model: LightningModule,
        batch: dict[str, Tensor],
        out: dict[str, Tensor],
        idx_dataset: int,
        n_samples: int,
        true_coords: Tensor,
        true_coords_resolved_mask: Tensor,
        expand_to_diffusion_samples: bool,
        batch_gathered: dict[str, Tensor],
        out_gathered: dict[str, Tensor],
    ):
        """Compute confidence metrics using triton pLDDT and serial PDE/PAE.

        Uses :func:`compute_plddt_mae_triton` for pLDDT MAE (avoids
        materialising the full N_token x N_R_set distance/mask matrices),
        and delegates to the serial ``compute_pde_mae`` / ``compute_pae_mae``
        for PDE and PAE.  ``token_to_rep_atom`` is deleted after pLDDT
        computation to free memory.
        """

        atom_pad_mask_1d = batch_gathered["atom_pad_mask_1d"]
        token_pad_mask_1d = batch_gathered["token_pad_mask_1d"]

        # Strip shared features from global to valid-only dimensions.
        # true_coords / true_coords_resolved_mask are passed as parameters;
        # rebind to stripped local copies so the caller's references stay global.
        true_coords = true_coords[..., atom_pad_mask_1d, :]
        true_coords_resolved_mask = true_coords_resolved_mask[..., atom_pad_mask_1d]
        batch_gathered["mol_type"] = batch_gathered["mol_type"][:, token_pad_mask_1d]
        batch_gathered["asym_id"] = batch_gathered["asym_id"][:, token_pad_mask_1d]
        batch_gathered["token_pad_mask"] = batch_gathered["token_pad_mask"][:, token_pad_mask_1d]
        batch_gathered["atom_pad_mask"] = batch_gathered["atom_pad_mask"][:, atom_pad_mask_1d]
        out_gathered["sample_atom_coords"] = out_gathered["sample_atom_coords"][:, atom_pad_mask_1d]

        K = batch["coords"].shape[1]  # ensemble dim is not sharded
        msg = "Confidence_prediction is not supported for num_ensembles_val > 1"
        assert K == 1, msg

        # Reconstruct atom_to_token
        batch_gathered["atom_to_token"] = reconstruct_atom_to_token_global(batch["atom_to_token"])[
            :, atom_pad_mask_1d, :
        ][:, :, token_pad_mask_1d]

        # Gather token-level confidence output
        out_gathered["plddt"] = gather_along_cp(out["plddt"])[:, token_pad_mask_1d]

        # Gather pair-level confidence outputs (pde, pae)
        for key in ("pde", "pae"):
            out_gathered[key] = gather_along_cp(out[key])[:, token_pad_mask_1d, :][:, :, token_pad_mask_1d]

        # Reconstruct diagonally-sharded mapping matrices
        batch_gathered["token_to_rep_atom"] = reconstruct_token_to_rep_atom_global(batch["token_to_rep_atom"])[
            :, token_pad_mask_1d, :
        ][:, :, atom_pad_mask_1d]

        r_set_to_rep_atom_gathered = reconstruct_r_set_to_rep_atom_global(batch["r_set_to_rep_atom"])
        # r_set_to_rep_atom_gathered dim (dim 1) may have per-shard padding; strip rows that are all zeros
        r_set_to_rep_atom_gathered_valid = r_set_to_rep_atom_gathered.any(dim=-1)  # [B, N_R_global], here B is 1
        batch_gathered["r_set_to_rep_atom"] = r_set_to_rep_atom_gathered[:, r_set_to_rep_atom_gathered_valid[0], :][
            :, :, atom_pad_mask_1d
        ]

        # Gather ensemble-aware frame features.
        frames_idx_gathered = gather_along_cp(batch["frames_idx"])
        if frames_idx_gathered.ndim == 4:
            # (B, E=1, T_padded, 3) → squeeze ensemble dim → (B, T_padded, 3)
            frames_idx_gathered = frames_idx_gathered.squeeze(1)
        batch_gathered["frames_idx"] = frames_idx_gathered[:, token_pad_mask_1d]

        frame_resolved_mask_gathered = gather_along_cp(batch["frame_resolved_mask"])
        if frame_resolved_mask_gathered.ndim == 3:
            # (B, E=1, T_padded) → squeeze ensemble dim → (B, T_padded)
            frame_resolved_mask_gathered = frame_resolved_mask_gathered.squeeze(1)
        batch_gathered["frame_resolved_mask"] = frame_resolved_mask_gathered[:, token_pad_mask_1d]

        mae_plddt_dicts: dict[str, list] = defaultdict(list)
        total_mae_plddt_dicts: dict[str, list] = defaultdict(list)
        mae_pde_dicts: dict[str, list] = defaultdict(list)
        total_mae_pde_dicts: dict[str, list] = defaultdict(list)
        mae_pae_dicts: dict[str, list] = defaultdict(list)
        total_mae_pae_dicts: dict[str, list] = defaultdict(list)

        if not expand_to_diffusion_samples:
            true_coords_resolved_mask = true_coords_resolved_mask.unsqueeze(0).repeat((n_samples, 1))

        for ensemble_idx in range(K):
            if expand_to_diffusion_samples:
                true_coords_k = true_coords[:, ensemble_idx]
            else:
                true_coords_k = true_coords[ensemble_idx].unsqueeze(0).repeat((n_samples, 1, 1))

            # pLDDT MAE via triton cdist_lddt (rectangular, per_atom)
            mae_plddt_dict, total_mae_plddt_dict = compute_plddt_mae_triton(
                pred_atom_coords=out_gathered["sample_atom_coords"],
                feats=batch_gathered,
                true_atom_coords=true_coords_k,
                pred_lddt=out_gathered["plddt"],
                true_coords_resolved_mask=true_coords_resolved_mask,
                multiplicity=n_samples,
            )
            for key in mae_plddt_dict:
                mae_plddt_dicts[key].append(mae_plddt_dict[key])
                total_mae_plddt_dicts[key].append(total_mae_plddt_dict[key])

            if ensemble_idx == K - 1:
                del batch_gathered["r_set_to_rep_atom"], out_gathered["plddt"]

            # PDE MAE via serial implementation
            mae_pde_dict, total_mae_pde_dict = compute_pde_mae(
                pred_atom_coords=out_gathered["sample_atom_coords"],
                feats=batch_gathered,
                true_atom_coords=true_coords_k,
                pred_pde=out_gathered["pde"],
                true_coords_resolved_mask=true_coords_resolved_mask,
                multiplicity=n_samples,
            )
            for key in mae_pde_dict:
                mae_pde_dicts[key].append(mae_pde_dict[key])
                total_mae_pde_dicts[key].append(total_mae_pde_dict[key])

            if ensemble_idx == K - 1:
                del batch_gathered["token_to_rep_atom"], out_gathered["pde"]

            # PAE MAE via serial implementation
            mae_pae_dict, total_mae_pae_dict = compute_pae_mae(
                pred_atom_coords=out_gathered["sample_atom_coords"],
                feats=batch_gathered,
                true_atom_coords=true_coords_k,
                pred_pae=out_gathered["pae"],
                true_coords_resolved_mask=true_coords_resolved_mask,
                multiplicity=n_samples,
            )
            for key in mae_pae_dict:
                mae_pae_dicts[key].append(mae_pae_dict[key])
                total_mae_pae_dicts[key].append(total_mae_pae_dict[key])

            if ensemble_idx == K - 1:
                del (
                    batch_gathered["atom_to_token"],
                    batch_gathered["mol_type"],
                    batch_gathered["asym_id"],
                    batch_gathered["frames_idx"],
                    batch_gathered["frame_resolved_mask"],
                    batch_gathered["token_pad_mask"],
                    batch_gathered["atom_pad_mask"],
                    batch_gathered["token_pad_mask_1d"],
                    batch_gathered["atom_pad_mask_1d"],
                    out_gathered["sample_atom_coords"],
                    out_gathered["pae"],
                )

        # Mean over ensembles
        for key in mae_plddt_dicts:
            mae_plddt_dicts[key] = torch.stack(mae_plddt_dicts[key], dim=0).mean(dim=0)
            total_mae_plddt_dicts[key] = torch.stack(total_mae_plddt_dicts[key], dim=0).mean(dim=0)

        for key in mae_pde_dicts:
            mae_pde_dicts[key] = torch.stack(mae_pde_dicts[key], dim=0).mean(dim=0)
            total_mae_pde_dicts[key] = torch.stack(total_mae_pde_dicts[key], dim=0).mean(dim=0)

        for key in mae_pae_dicts:
            mae_pae_dicts[key] = torch.stack(mae_pae_dicts[key], dim=0).mean(dim=0)
            total_mae_pae_dicts[key] = torch.stack(total_mae_pae_dicts[key], dim=0).mean(dim=0)

        return (
            mae_plddt_dicts,
            total_mae_plddt_dicts,
            mae_pde_dicts,
            total_mae_pde_dicts,
            mae_pae_dicts,
            total_mae_pae_dicts,
        )

    def _dp_all_reduce_metrics(self, dp_group: torch.distributed.ProcessGroup) -> None:
        """All-reduce MeanMetric internal states across DP ranks.

        In our CP code, lightning is blind to our distributed computation context
        so MeanMetric is strictly single device and we should not rely on
        MeanMetric.compute() to get inter-DP all_reduce mean.

        Each DP rank accumulates metrics from its own subset of validation
        batches, so ``MeanMetric.mean_value`` (weighted sum) and
        ``MeanMetric.weight`` (total weight) are local.  Before calling
        ``.compute()`` we must sum both across DP ranks so the resulting
        ratio is the global weighted mean — matching the dev-v2
        ``_DP_all_reduce_mean`` semantics.

        This also handles the key-synchronisation problem: different DP
        ranks may see different molecular modalities, so some
        ``MeanMetric`` instances may have ``weight == 0`` on a rank that
        never saw that modality.  After all-reduce, a metric with
        ``weight == 0`` globally will produce ``NaN`` from ``.compute()``,
        which the serial ``common_on_epoch_end`` already maps to ``0.0``.

        """
        for m in self.modules():
            if isinstance(m, MeanMetric):
                torch.distributed.all_reduce(m.mean_value, op=torch.distributed.ReduceOp.SUM, group=dp_group)
                torch.distributed.all_reduce(m.weight, op=torch.distributed.ReduceOp.SUM, group=dp_group)
                m._computed = None
            elif not isinstance(m, (nn.ModuleDict, nn.ModuleList, DistributedValidator)):
                raise ValueError(f"Only support MeanMetric, got {type(m)}")

    def common_on_epoch_end(self, model: LightningModule) -> None:
        """Aggregate metrics at epoch end with DP all-reduce.

        1. All-reduces every ``MeanMetric``'s internal ``mean_value`` and
           ``weight`` across DP ranks so that ``.compute()`` returns the
           global weighted mean (equivalent to dev-v2's
           ``_DP_all_reduce_mean``).
        2. Wraps ``model.log`` with ``sync_dist=False`` to avoid NCCL
           collectives inside Lightning's logging (which can deadlock when
           validation batches are unevenly distributed across DP ranks).
        3. Delegates to the serial ``common_on_epoch_end`` for the actual
           compute / log / reset cycle.
        """
        # In our CP code, lightning is blind to our distributed computation context
        # so MeanMetric is strictly single device and we should not rely on
        # MeanMetric.compute() to get inter-DP all_reduce mean.
        self._dp_all_reduce_metrics(model.dp_group)

        original_log = model.log

        def _log_no_sync(*args, **kwargs):
            kwargs["sync_dist"] = False
            return original_log(*args, **kwargs)

        model.log = _log_no_sync  # type: ignore[assignment]
        try:
            for idx_dataset in range(self.num_val_datasets):
                dataset_name_ori = self.val_names[idx_dataset]
                dataset_name = "" if dataset_name_ori == "RCSB" else f"__{dataset_name_ori}"
                if self.clash_score_metrics:
                    for m in ("clash_atoms_count", "clash_atoms_fraction"):
                        val = self.folding_metrics["clash_score"][idx_dataset][m].compute()
                        val = 0.0 if torch.isnan(val) else val.item()
                        self.folding_metrics["clash_score"][idx_dataset][m].reset()
                        model.log(f"val/{m}{dataset_name}", val)
                if self.rmsd_metrics:
                    avg_rmsd = self.folding_metrics["rmsd"][idx_dataset]["rmsd"].compute()
                    avg_rmsd = 0.0 if torch.isnan(avg_rmsd) else avg_rmsd.item()
                    self.folding_metrics["rmsd"][idx_dataset]["rmsd"].reset()
                    model.log(f"val/rmsd{dataset_name}", avg_rmsd)
            super().common_on_epoch_end(model)
        finally:
            model.log = original_log  # type: ignore[assignment]

    def common_val_step(
        self,
        model: LightningModule,
        batch: dict[str, DTensor],
        out: dict[str, DTensor],
        idx_dataset: int,
        expand_to_diffusion_samples: bool,
        transpose_comm: TransposeComm,
    ) -> None:
        """Run a common validation step with DTensor inputs.

        Gathers DTensor features once, then delegates metric computation
        to the serial base class methods operating on plain tensors.

        Parameters
        ----------
        model : LightningModule
            The distributed Boltz2 model (used for accessing serial
            ``get_true_coordinates`` via the model_serial attribute).
        batch : dict[str, DTensor]
            Batch features as DTensors.
        out : dict[str, DTensor]
            Model outputs as DTensors.
        idx_dataset : int
            Global dataset index.
        expand_to_diffusion_samples : bool
            Whether to expand coordinates to diffusion samples.
        """
        symmetry_correction = model.val_group_mapper[idx_dataset]["symmetry_correction"]
        idx_dataset = self.get_local_val_index(model, idx_dataset)
        n_samples = model.validation_args.diffusion_samples

        # Compute distogram loss and update metrics
        val_disto_loss = self.compute_disto_loss(model, out, batch, idx_dataset, transpose_comm)
        self.folding_metrics["disto_loss"][idx_dataset]["disto_loss"].update(val_disto_loss)

        # Compute distogram lddt and update metrics
        disto_lddt_dict, disto_total_dict = self.compute_disto_lddt(model, batch, out, idx_dataset)

        token_pad_mask = gather_along_cp(batch["token_pad_mask"].bool())
        atom_pad_mask = gather_along_cp(batch["atom_pad_mask"].bool())
        token_pad_mask_1d = token_pad_mask[0]
        atom_pad_mask_1d = atom_pad_mask[0]

        # Get true coords (DTensors) and gather to plain tensors.
        # All gathered tensors are kept at global (unstripped) dimensions so
        # that lDDT (atom_to_token is global) and clash/PB (atom indices are
        # global) can consume them directly.  Confidence metrics strip
        # internally via the 1D pad masks.
        return_dict = self.get_true_coords(
            model,
            batch,
            out,
            n_samples,
            symmetry_correction,
            expand_to_diffusion_samples=expand_to_diffusion_samples,
        )
        true_coords = gather_along_cp(return_dict["true_coords"])
        true_coords_resolved_mask = gather_along_cp(return_dict["true_coords_resolved_mask"])

        mol_type = gather_along_cp(batch["mol_type"])
        asym_id = gather_along_cp(batch["asym_id"])
        sample_atom_coords = gather_along_cp(out["sample_atom_coords"])

        batch_gathered = {
            "token_pad_mask": token_pad_mask,
            "atom_pad_mask": atom_pad_mask,
            "token_pad_mask_1d": token_pad_mask_1d,
            "atom_pad_mask_1d": atom_pad_mask_1d,
            "mol_type": mol_type,
            "asym_id": asym_id,
        }
        out_gathered = {
            "sample_atom_coords": sample_atom_coords,
        }

        # Get lddt metrics (all inputs at global dimensions)
        K = batch["coords"].shape[1]
        if expand_to_diffusion_samples:
            resolved_mask_for_lddt = true_coords_resolved_mask
        else:
            resolved_mask_for_lddt = true_coords_resolved_mask.squeeze(0)
        all_lddt_dict, all_total_dict = get_lddt_metrics(
            atom_to_token_dtensor=batch["atom_to_token"],
            num_conformers=K,
            n_samples=n_samples,
            true_coords=true_coords,
            true_coords_resolved_mask=resolved_mask_for_lddt,
            mol_type=mol_type,
            asym_id=asym_id,
            sample_atom_coords=sample_atom_coords,
            expand_to_diffusion_samples=expand_to_diffusion_samples,
        )

        # Get physical realism metrics
        if self.physicalism_metrics:
            pair_clash_dict, pair_total_dict = self.get_clash_metrics(batch, out, batch_gathered, out_gathered)
            pb_failure_dict, pb_total_dict = self.get_pb_metrics(batch, out, batch_gathered, out_gathered)
        else:
            pair_clash_dict, pair_total_dict = None, None
            pb_failure_dict, pb_total_dict = None, None

        # Filtering based on confidence
        if model.confidence_prediction and n_samples > 1:
            (
                mae_plddt_dicts,
                total_mae_plddt_dicts,
                mae_pde_dicts,
                total_mae_pde_dicts,
                mae_pae_dicts,
                total_mae_pae_dicts,
            ) = self.get_confidence_metrics(
                model,
                batch,
                out,
                idx_dataset,
                n_samples,
                true_coords,
                true_coords_resolved_mask,
                expand_to_diffusion_samples,
                batch_gathered,
                out_gathered,
            )

        # Compute RMSD on gathered plain tensors (first conformer only)
        if self.rmsd_metrics:
            atom_to_token_global = reconstruct_atom_to_token_global(batch["atom_to_token"])

            if expand_to_diffusion_samples:
                atom_coords_rmsd = true_coords[:, 0]
                resolved_mask_rmsd = true_coords_resolved_mask
            else:
                atom_coords_rmsd = true_coords[0].unsqueeze(0).repeat((n_samples, 1, 1))
                resolved_mask_rmsd = true_coords_resolved_mask.repeat((n_samples, 1))

            rmsd_val, _, _ = weighted_minimum_rmsd_single(
                pred_atom_coords=sample_atom_coords,
                atom_coords=atom_coords_rmsd,
                atom_mask=resolved_mask_rmsd.float(),
                atom_to_token=atom_to_token_global.expand(n_samples, -1, -1)
                if atom_to_token_global.shape[0] == 1 and n_samples > 1
                else atom_to_token_global,
                mol_type=mol_type.expand(n_samples, -1) if mol_type.shape[0] == 1 and n_samples > 1 else mol_type,
            )
            del atom_to_token_global, atom_coords_rmsd, resolved_mask_rmsd

        if self.clash_score_metrics:
            clash_cutoff = model.validation_args.get("clash_cutoff")
            if clash_cutoff:
                coords_repr = gather_along_cp(
                    single_repr_rep_atom_to_token(
                        out["sample_atom_coords"],
                        batch["token_to_rep_atom"],
                    )
                )
                clash_count, clash_frac = clash_score(
                    coords_repr=coords_repr,
                    token_pad_mask=token_pad_mask,
                    multiplicity=n_samples,
                    clash_cutoff=clash_cutoff,
                )
                self.folding_metrics["clash_score"][idx_dataset]["clash_atoms_count"].update(clash_count.mean())
                self.folding_metrics["clash_score"][idx_dataset]["clash_atoms_fraction"].update(clash_frac.mean())

        # Update folding metrics
        self.update_lddt_rmsd_metrics(
            {"contact_conditioning": gather_along_cp(batch["contact_conditioning"]), "coords": batch["coords"]},
            all_lddt_dict,
            all_total_dict,
            disto_lddt_dict,
            disto_total_dict,
            idx_dataset,
        )
        if self.rmsd_metrics:
            self.folding_metrics["rmsd"][idx_dataset]["rmsd"].update(rmsd_val.mean())

        # Update physical realism metrics
        if self.physicalism_metrics:
            self.update_physcialism_metrics(
                pair_clash_dict,
                pair_total_dict,
                pb_failure_dict,
                pb_total_dict,
                idx_dataset,
            )

        # Update confidence metrics
        if model.confidence_prediction and n_samples > 1:
            # Pass through scalar confidence outputs (already replicated) to plain tensors
            confidence_out = {}
            for key in (
                "complex_plddt",
                "complex_iplddt",
                "complex_pde",
                "complex_ipde",
                "ptm",
                "iptm",
                "ligand_iptm",
                "protein_iptm",
            ):
                if key in out:
                    val = out[key]
                    confidence_out[key] = val.to_local() if isinstance(val, DTensor) else val

            self.update_confidence_metrics(
                batch,
                confidence_out,
                idx_dataset,
                n_samples,
                all_lddt_dict,
                all_total_dict,
                mae_plddt_dicts,
                total_mae_plddt_dicts,
                mae_pde_dicts,
                total_mae_pde_dicts,
                mae_pae_dicts,
                total_mae_pae_dicts,
                pair_clash_dict,
                pair_total_dict,
                pb_failure_dict,
                pb_total_dict,
                physicalism_metrics=self.physicalism_metrics,
            )
