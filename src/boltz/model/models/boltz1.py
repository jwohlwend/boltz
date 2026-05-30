import gc
import math
import random
from pathlib import Path
from typing import Any, Optional
from typing import Iterable

import numpy as np

import torch
import torch._dynamo
from pytorch_lightning import LightningModule
from torch import Tensor, nn
from torchmetrics import MeanMetric

import boltz.model.layers.initialize as init
from boltz.data.types import Connection, Structure
from boltz.data import const
from boltz.data.feature.symmetry import (
    minimum_lddt_symmetry_coords,
    minimum_symmetry_coords,
)
from lightning_fabric.utilities.apply_func import move_data_to_device
from boltz.model.loss.confidence import confidence_loss
from boltz.model.loss.distogram import distogram_loss
from boltz.model.loss.validation import (
    SampleMetrics,
    compute_pae_mae,
    compute_pde_mae,
    compute_plddt_mae,
    factored_lddt_loss,
    factored_token_lddt_dist_loss,
    weighted_minimum_rmsd,
    weighted_minimum_rmsd_single,
)
from boltz.model.modules.confidence import ConfidenceModule
from boltz.data.write.writer import BoltzWriter
from boltz.model.modules.diffusion import AtomDiffusion
from boltz.model.modules.encoders import RelativePositionEncoder
from boltz.model.modules.trunk import (
    DistogramModule,
    InputEmbedder,
    MSAModule,
    PairformerModule,
)
from boltz.model.modules.utils import ExponentialMovingAverage
from boltz.model.optim.scheduler import AlphaFoldLRScheduler


class Boltz1(LightningModule):
    """Boltz1 model."""

    def __init__(  # noqa: PLR0915, C901, PLR0912
        self,
        atom_s: int,
        atom_z: int,
        token_s: int,
        token_z: int,
        num_bins: int,
        training_args: dict[str, Any],
        validation_args: dict[str, Any],
        embedder_args: dict[str, Any],
        msa_args: dict[str, Any],
        pairformer_args: dict[str, Any],
        score_model_args: dict[str, Any],
        diffusion_process_args: dict[str, Any],
        diffusion_loss_args: dict[str, Any],
        confidence_model_args: dict[str, Any],
        atom_feature_dim: int = 128,
        confidence_prediction: bool = False,
        confidence_imitate_trunk: bool = False,
        alpha_pae: float = 0.0,
        structure_prediction_training: bool = True,
        atoms_per_window_queries: int = 32,
        atoms_per_window_keys: int = 128,
        compile_pairformer: bool = False,
        compile_structure: bool = False,
        compile_confidence: bool = False,
        nucleotide_rmsd_weight: float = 5.0,
        ligand_rmsd_weight: float = 10.0,
        no_msa: bool = False,
        no_atom_encoder: bool = False,
        ema: bool = False,
        ema_decay: float = 0.999,
        min_dist: float = 2.0,
        max_dist: float = 22.0,
        predict_args: Optional[dict[str, Any]] = None,
        steering_args: Optional[dict[str, Any]] = None,
        use_kernels: bool = False,
    ) -> None:
        super().__init__()

        self.save_hyperparameters()

        self.lddt = nn.ModuleDict()
        self.disto_lddt = nn.ModuleDict()
        self.complex_lddt = nn.ModuleDict()
        if confidence_prediction:
            self.top1_lddt = nn.ModuleDict()
            self.iplddt_top1_lddt = nn.ModuleDict()
            self.ipde_top1_lddt = nn.ModuleDict()
            self.pde_top1_lddt = nn.ModuleDict()
            self.ptm_top1_lddt = nn.ModuleDict()
            self.iptm_top1_lddt = nn.ModuleDict()
            self.ligand_iptm_top1_lddt = nn.ModuleDict()
            self.protein_iptm_top1_lddt = nn.ModuleDict()
            self.avg_lddt = nn.ModuleDict()
            self.plddt_mae = nn.ModuleDict()
            self.pde_mae = nn.ModuleDict()
            self.pae_mae = nn.ModuleDict()
        for m in const.out_types + ["pocket_ligand_protein"]:
            self.lddt[m] = MeanMetric()
            self.disto_lddt[m] = MeanMetric()
            self.complex_lddt[m] = MeanMetric()
            if confidence_prediction:
                self.top1_lddt[m] = MeanMetric()
                self.iplddt_top1_lddt[m] = MeanMetric()
                self.ipde_top1_lddt[m] = MeanMetric()
                self.pde_top1_lddt[m] = MeanMetric()
                self.ptm_top1_lddt[m] = MeanMetric()
                self.iptm_top1_lddt[m] = MeanMetric()
                self.ligand_iptm_top1_lddt[m] = MeanMetric()
                self.protein_iptm_top1_lddt[m] = MeanMetric()
                self.avg_lddt[m] = MeanMetric()
                self.pde_mae[m] = MeanMetric()
                self.pae_mae[m] = MeanMetric()
        for m in const.out_single_types:
            if confidence_prediction:
                self.plddt_mae[m] = MeanMetric()
        self.rmsd = MeanMetric()
        self.best_rmsd = MeanMetric()

        self.train_confidence_loss_logger = MeanMetric()
        self.train_confidence_loss_dict_logger = nn.ModuleDict()
        for m in [
            "plddt_loss",
            "resolved_loss",
            "pde_loss",
            "pae_loss",
        ]:
            self.train_confidence_loss_dict_logger[m] = MeanMetric()

        self.ema = None
        self.use_ema = ema
        self.ema_decay = ema_decay

        self.training_args = training_args
        self.validation_args = validation_args
        self.diffusion_loss_args = diffusion_loss_args
        self.predict_args = predict_args
        self.steering_args = steering_args

        self.use_kernels = use_kernels
        self.validation_writer: Optional[BoltzWriter] = None

        self.nucleotide_rmsd_weight = nucleotide_rmsd_weight
        self.ligand_rmsd_weight = ligand_rmsd_weight

        self.num_bins = num_bins
        self.min_dist = min_dist
        self.max_dist = max_dist
        self.is_pairformer_compiled = False

        # Input projections
        s_input_dim = (
            token_s + 2 * const.num_tokens + 1 + len(const.pocket_contact_info)
        )
        self.s_init = nn.Linear(s_input_dim, token_s, bias=False)
        self.z_init_1 = nn.Linear(s_input_dim, token_z, bias=False)
        self.z_init_2 = nn.Linear(s_input_dim, token_z, bias=False)

        # Input embeddings
        full_embedder_args = {
            "atom_s": atom_s,
            "atom_z": atom_z,
            "token_s": token_s,
            "token_z": token_z,
            "atoms_per_window_queries": atoms_per_window_queries,
            "atoms_per_window_keys": atoms_per_window_keys,
            "atom_feature_dim": atom_feature_dim,
            "no_atom_encoder": no_atom_encoder,
            **embedder_args,
        }
        self.input_embedder = InputEmbedder(**full_embedder_args)
        self.rel_pos = RelativePositionEncoder(token_z)
        self.token_bonds = nn.Linear(1, token_z, bias=False)

        # Normalization layers
        self.s_norm = nn.LayerNorm(token_s)
        self.z_norm = nn.LayerNorm(token_z)

        # Recycling projections
        self.s_recycle = nn.Linear(token_s, token_s, bias=False)
        self.z_recycle = nn.Linear(token_z, token_z, bias=False)
        init.gating_init_(self.s_recycle.weight)
        init.gating_init_(self.z_recycle.weight)

        # Pairwise stack
        self.no_msa = no_msa
        if not no_msa:
            self.msa_module = MSAModule(
                token_z=token_z,
                s_input_dim=s_input_dim,
                **msa_args,
            )
        self.pairformer_module = PairformerModule(token_s, token_z, **pairformer_args)
        if compile_pairformer:
            # Big models hit the default cache limit (8)
            self.is_pairformer_compiled = True
            torch._dynamo.config.cache_size_limit = 512
            torch._dynamo.config.accumulated_cache_size_limit = 512
            self.pairformer_module = torch.compile(
                self.pairformer_module,
                dynamic=False,
                fullgraph=False,
            )

        # Output modules
        use_accumulate_token_repr = (
            confidence_prediction
            and "use_s_diffusion" in confidence_model_args
            and confidence_model_args["use_s_diffusion"]
        )
        self.structure_module = AtomDiffusion(
            score_model_args={
                "token_z": token_z,
                "token_s": token_s,
                "atom_z": atom_z,
                "atom_s": atom_s,
                "atoms_per_window_queries": atoms_per_window_queries,
                "atoms_per_window_keys": atoms_per_window_keys,
                "atom_feature_dim": atom_feature_dim,
                **score_model_args,
            },
            compile_score=compile_structure,
            accumulate_token_repr=use_accumulate_token_repr,
            **diffusion_process_args,
        )
        self.distogram_module = DistogramModule(token_z, num_bins)
        self.confidence_prediction = confidence_prediction
        self.alpha_pae = alpha_pae

        self.structure_prediction_training = structure_prediction_training
        self.confidence_imitate_trunk = confidence_imitate_trunk
        if self.confidence_prediction:
            if self.confidence_imitate_trunk:
                self.confidence_module = ConfidenceModule(
                    token_s,
                    token_z,
                    compute_pae=alpha_pae > 0,
                    imitate_trunk=True,
                    pairformer_args=pairformer_args,
                    full_embedder_args=full_embedder_args,
                    msa_args=msa_args,
                    **confidence_model_args,
                )
            else:
                self.confidence_module = ConfidenceModule(
                    token_s,
                    token_z,
                    compute_pae=alpha_pae > 0,
                    **confidence_model_args,
                )
            if compile_confidence:
                self.confidence_module = torch.compile(
                    self.confidence_module, dynamic=False, fullgraph=False
                )

        # Remove grad from weights they are not trained for ddp
        if not structure_prediction_training:
            for name, param in self.named_parameters():
                if name.split(".")[0] != "confidence_module":
                    param.requires_grad = False

    def setup(self, stage: str) -> None:
        """Set the model for training, validation and inference."""
        if stage == "predict" and not (
            torch.cuda.is_available()
            and torch.cuda.get_device_properties(torch.device("cuda")).major >= 8.0  # noqa: PLR2004
        ):
            self.use_kernels = False

    def transfer_batch_to_device(
        self, batch: dict[str, Any], device: torch.device, dataloader_idx: int
    ) -> dict[str, Any]:
        """Move tensor entries to device while keeping CPU-only metadata intact."""
        if batch is None:
            return batch

        cpu_only_keys = {"record", "base_structure"}
        cpu_entries = {}
        tensor_entries = {}

        for key, value in batch.items():
            if key in cpu_only_keys:
                cpu_entries[key] = value
            else:
                tensor_entries[key] = value

        moved = move_data_to_device(tensor_entries, device)
        moved.update(cpu_entries)
        return moved

    def forward(
        self,
        feats: dict[str, Tensor],
        recycling_steps: int = 0,
        num_sampling_steps: Optional[int] = None,
        multiplicity_diffusion_train: int = 1,
        diffusion_samples: int = 1,
        max_parallel_samples: Optional[int] = 1,
        run_confidence_sequentially: bool = False,
    ) -> dict[str, Tensor]:
        dict_out = {}

        # Compute input embeddings
        with torch.set_grad_enabled(
            self.training and self.structure_prediction_training
        ):
            s_inputs = self.input_embedder(feats)

            # Initialize the sequence and pairwise embeddings
            s_init = self.s_init(s_inputs)
            z_init = (
                self.z_init_1(s_inputs)[:, :, None]
                + self.z_init_2(s_inputs)[:, None, :]
            )
            relative_position_encoding = self.rel_pos(feats)
            z_init = z_init + relative_position_encoding
            z_init = z_init + self.token_bonds(feats["token_bonds"].float())

            # Perform rounds of the pairwise stack
            s = torch.zeros_like(s_init)
            z = torch.zeros_like(z_init)

            # Compute pairwise mask
            mask = feats["token_pad_mask"].float()
            pair_mask = mask[:, :, None] * mask[:, None, :]

            for i in range(recycling_steps + 1):
                with torch.set_grad_enabled(self.training and (i == recycling_steps)):
                    # Fixes an issue with unused parameters in autocast
                    if (
                        self.training
                        and (i == recycling_steps)
                        and torch.is_autocast_enabled()
                    ):
                        torch.clear_autocast_cache()

                    # Apply recycling
                    s = s_init + self.s_recycle(self.s_norm(s))
                    z = z_init + self.z_recycle(self.z_norm(z))

                    # Compute pairwise stack
                    if not self.no_msa:
                        z = z + self.msa_module(
                            z, s_inputs, feats, use_kernels=self.use_kernels
                        )

                    # Revert to uncompiled version for validation
                    if self.is_pairformer_compiled and not self.training:
                        pairformer_module = self.pairformer_module._orig_mod  # noqa: SLF001
                    else:
                        pairformer_module = self.pairformer_module

                    s, z = pairformer_module(
                        s,
                        z,
                        mask=mask,
                        pair_mask=pair_mask,
                        use_kernels=self.use_kernels,
                    )

            pdistogram = self.distogram_module(z)
            dict_out = {
                "pdistogram": pdistogram,
                "s": s,
                "z": z,
            }

        # Compute structure module
        if self.training and self.structure_prediction_training:
            dict_out.update(
                self.structure_module(
                    s_trunk=s,
                    z_trunk=z,
                    s_inputs=s_inputs,
                    feats=feats,
                    relative_position_encoding=relative_position_encoding,
                    multiplicity=multiplicity_diffusion_train,
                )
            )

        if (not self.training) or self.confidence_prediction:
            dict_out.update(
                self.structure_module.sample(
                    s_trunk=s,
                    z_trunk=z,
                    s_inputs=s_inputs,
                    feats=feats,
                    relative_position_encoding=relative_position_encoding,
                    num_sampling_steps=num_sampling_steps,
                    atom_mask=feats["atom_pad_mask"],
                    multiplicity=diffusion_samples,
                    max_parallel_samples=max_parallel_samples,
                    train_accumulate_token_repr=self.training,
                    steering_args=self.steering_args,
                )
            )

        if self.confidence_prediction:
            dict_out.update(
                self.confidence_module(
                    s_inputs=s_inputs.detach(),
                    s=s.detach(),
                    z=z.detach(),
                    s_diffusion=(
                        dict_out["diff_token_repr"]
                        if self.confidence_module.use_s_diffusion
                        else None
                    ),
                    x_pred=dict_out["sample_atom_coords"].detach(),
                    feats=feats,
                    pred_distogram_logits=dict_out["pdistogram"].detach(),
                    multiplicity=diffusion_samples,
                    run_sequentially=run_confidence_sequentially,
                    use_kernels=self.use_kernels,
                )
            )
        if self.confidence_prediction and self.confidence_module.use_s_diffusion:
            dict_out.pop("diff_token_repr", None)
        return dict_out

    def get_true_coordinates(
        self,
        batch,
        out,
        diffusion_samples,
        symmetry_correction,
        lddt_minimization=True,
        alignment_mask: Optional[Tensor] = None,
        rmsd_mask: Optional[Tensor] = None,
    ):
        masked_rmsds = None
        best_masked_rmsds = None

        if symmetry_correction:
            min_coords_routine = (
                minimum_lddt_symmetry_coords
                if lddt_minimization
                else minimum_symmetry_coords
            )
            true_coords = []
            true_coords_resolved_mask = []
            rmsds, best_rmsds = [], []
            for idx in range(batch["token_index"].shape[0]):
                best_rmsd = float("inf")
                for rep in range(diffusion_samples):
                    i = idx * diffusion_samples + rep
                    best_true_coords, rmsd, best_true_coords_resolved_mask = (
                        min_coords_routine(
                            coords=out["sample_atom_coords"][i : i + 1],
                            feats=batch,
                            index_batch=idx,
                            nucleotide_weight=self.nucleotide_rmsd_weight,
                            ligand_weight=self.ligand_rmsd_weight,
                        )
                    )
                    rmsds.append(rmsd)
                    true_coords.append(best_true_coords)
                    true_coords_resolved_mask.append(best_true_coords_resolved_mask)
                    if rmsd < best_rmsd:
                        best_rmsd = rmsd
                best_rmsds.append(best_rmsd)
            true_coords = torch.cat(true_coords, dim=0)
            true_coords_resolved_mask = torch.cat(true_coords_resolved_mask, dim=0)
            masked_rmsds = None
            best_masked_rmsds = None
        else:
            true_coords = (
                batch["coords"].squeeze(1).repeat_interleave(diffusion_samples, 0)
            )

            true_coords_resolved_mask = batch["atom_resolved_mask"].repeat_interleave(
                diffusion_samples, 0
            )
            rmsds, best_rmsds = weighted_minimum_rmsd(
                out["sample_atom_coords"],
                batch,
                multiplicity=diffusion_samples,
                nucleotide_weight=self.nucleotide_rmsd_weight,
                ligand_weight=self.ligand_rmsd_weight,
            )

            align_mask = alignment_mask
            if align_mask is not None:
                target_len = out["sample_atom_coords"].shape[0]
                if align_mask.shape[0] != target_len:
                    if align_mask.shape[0] * diffusion_samples == target_len:
                        align_mask = align_mask.repeat_interleave(diffusion_samples, 0)
                    else:
                        raise ValueError(
                            "Alignment mask has unexpected first dimension "
                            f"{align_mask.shape[0]} (expected {target_len})."
                        )
            calc_mask = rmsd_mask
            if calc_mask is not None:
                target_len = out["sample_atom_coords"].shape[0]
                if calc_mask.shape[0] != target_len:
                    if calc_mask.shape[0] * diffusion_samples == target_len:
                        calc_mask = calc_mask.repeat_interleave(diffusion_samples, 0)
                    else:
                        raise ValueError(
                            "RMSD mask has unexpected first dimension "
                            f"{calc_mask.shape[0]} (expected {target_len})."
                        )

            if align_mask is not None or calc_mask is not None:
                masked_rmsds, best_masked_rmsds = weighted_minimum_rmsd(
                    out["sample_atom_coords"],
                    batch,
                    multiplicity=diffusion_samples,
                    nucleotide_weight=self.nucleotide_rmsd_weight,
                    ligand_weight=self.ligand_rmsd_weight,
                    alignment_mask=align_mask,
                    rmsd_mask=calc_mask,
                )
            else:
                masked_rmsds, best_masked_rmsds = rmsds, best_rmsds

        return (
            true_coords,
            rmsds,
            best_rmsds,
            true_coords_resolved_mask,
            masked_rmsds,
            best_masked_rmsds,
        )

    def training_step(self, batch: dict[str, Tensor], batch_idx: int) -> Tensor:
        # Sample recycling steps
        recycling_steps = random.randint(0, self.training_args.recycling_steps)

        # Compute the forward pass
        out = self(
            feats=batch,
            recycling_steps=recycling_steps,
            num_sampling_steps=self.training_args.sampling_steps,
            multiplicity_diffusion_train=self.training_args.diffusion_multiplicity,
            diffusion_samples=self.training_args.diffusion_samples,
        )

        # Compute losses
        if self.structure_prediction_training:
            disto_loss, _ = distogram_loss(
                out,
                batch,
            )
            try:
                diffusion_loss_dict = self.structure_module.compute_loss(
                    batch,
                    out,
                    multiplicity=self.training_args.diffusion_multiplicity,
                    **self.diffusion_loss_args,
                )
            except Exception as e:
                print(f"Skipping batch {batch_idx} due to error: {e}")
                return None

        else:
            disto_loss = 0.0
            diffusion_loss_dict = {"loss": 0.0, "loss_breakdown": {}}

        if self.confidence_prediction:
            # confidence model symmetry correction
            true_coords, _, _, true_coords_resolved_mask, _, _ = self.get_true_coordinates(
                batch,
                out,
                diffusion_samples=self.training_args.diffusion_samples,
                symmetry_correction=self.training_args.symmetry_correction,
            )

            confidence_loss_dict = confidence_loss(
                out,
                batch,
                true_coords,
                true_coords_resolved_mask,
                alpha_pae=self.alpha_pae,
                multiplicity=self.training_args.diffusion_samples,
            )
        else:
            confidence_loss_dict = {
                "loss": torch.tensor(0.0, device=batch["token_index"].device),
                "loss_breakdown": {},
            }

        # Aggregate losses
        loss = (
            self.training_args.confidence_loss_weight * confidence_loss_dict["loss"]
            + self.training_args.diffusion_loss_weight * diffusion_loss_dict["loss"]
            + self.training_args.distogram_loss_weight * disto_loss
        )
        # Log losses
        self.log("train/distogram_loss", disto_loss)
        self.log("train/diffusion_loss", diffusion_loss_dict["loss"])
        for k, v in diffusion_loss_dict["loss_breakdown"].items():
            self.log(f"train/{k}", v)

        if self.confidence_prediction:
            self.train_confidence_loss_logger.update(
                confidence_loss_dict["loss"].detach()
            )

            for k in self.train_confidence_loss_dict_logger.keys():
                self.train_confidence_loss_dict_logger[k].update(
                    confidence_loss_dict["loss_breakdown"][k].detach()
                    if torch.is_tensor(confidence_loss_dict["loss_breakdown"][k])
                    else confidence_loss_dict["loss_breakdown"][k]
                )
        self.log("train/loss", loss)
        self.training_log()
        return loss

    def training_log(self):
        self.log("train/grad_norm", self.gradient_norm(self), prog_bar=False)
        self.log("train/param_norm", self.parameter_norm(self), prog_bar=False)

        lr = self.trainer.optimizers[0].param_groups[0]["lr"]
        self.log("lr", lr, prog_bar=False)

        self.log(
            "train/grad_norm_msa_module",
            self.gradient_norm(self.msa_module),
            prog_bar=False,
        )
        self.log(
            "train/param_norm_msa_module",
            self.parameter_norm(self.msa_module),
            prog_bar=False,
        )

        self.log(
            "train/grad_norm_pairformer_module",
            self.gradient_norm(self.pairformer_module),
            prog_bar=False,
        )
        self.log(
            "train/param_norm_pairformer_module",
            self.parameter_norm(self.pairformer_module),
            prog_bar=False,
        )

        self.log(
            "train/grad_norm_structure_module",
            self.gradient_norm(self.structure_module),
            prog_bar=False,
        )
        self.log(
            "train/param_norm_structure_module",
            self.parameter_norm(self.structure_module),
            prog_bar=False,
        )

        if self.confidence_prediction:
            self.log(
                "train/grad_norm_confidence_module",
                self.gradient_norm(self.confidence_module),
                prog_bar=False,
            )
            self.log(
                "train/param_norm_confidence_module",
                self.parameter_norm(self.confidence_module),
                prog_bar=False,
            )

    def on_train_epoch_end(self):
        self.log(
            "train/confidence_loss",
            self.train_confidence_loss_logger,
            prog_bar=False,
            on_step=False,
            on_epoch=True,
        )
        for k, v in self.train_confidence_loss_dict_logger.items():
            self.log(f"train/{k}", v, prog_bar=False, on_step=False, on_epoch=True)

    def gradient_norm(self, module) -> float:
        # Only compute over parameters that are being trained
        parameters = filter(lambda p: p.requires_grad, module.parameters())
        parameters = filter(lambda p: p.grad is not None, parameters)
        norm = torch.tensor([p.grad.norm(p=2) ** 2 for p in parameters]).sum().sqrt()
        return norm

    def parameter_norm(self, module) -> float:
        # Only compute over parameters that are being trained
        parameters = filter(lambda p: p.requires_grad, module.parameters())
        norm = torch.tensor([p.norm(p=2) ** 2 for p in parameters]).sum().sqrt()
        return norm

    def _build_alignment_masks(
        self,
        batch: dict[str, Tensor],
        record_ids: list[str],
        base_structures: list[Optional[Structure]],
        batch_idx: int,
        precomputed_align_masks: Optional[list[Optional[np.ndarray]]] = None,
        precomputed_rmsd_masks: Optional[list[Optional[np.ndarray]]] = None,
    ) -> tuple[Tensor, Tensor]:
        batch_size, num_atoms, _ = batch["atom_to_token"].shape
        device = batch["atom_pad_mask"].device

        heavy_calc_mask = torch.zeros((batch_size, num_atoms), dtype=torch.bool, device=device)
        peptide_calc_mask = torch.zeros_like(heavy_calc_mask)

        if precomputed_align_masks is None or precomputed_rmsd_masks is None:
            raise ValueError("Precomputed masks are required but were not provided.")

        for structure_idx in range(batch_size):
            record_id = (
                record_ids[structure_idx]
                if structure_idx < len(record_ids)
                else f"batch_{batch_idx}_{structure_idx}"
            )
            structure = (
                base_structures[structure_idx]
                if structure_idx < len(base_structures)
                else None
            )
            pre_align = precomputed_align_masks[structure_idx]
            pre_rmsd = precomputed_rmsd_masks[structure_idx]

            atom_to_token = batch["atom_to_token"][structure_idx].bool()
            entity_ids = batch["entity_id"][structure_idx]
            present_atom = batch["atom_pad_mask"][structure_idx].bool()

            unique_entities = torch.unique(entity_ids)
            res_counts = {int(e.item()): int((entity_ids == e).sum().item()) for e in unique_entities}
            if not res_counts:
                print(f"Warning: no entities found for structure {structure_idx}.")
                continue

            peptide_entity = min(res_counts, key=res_counts.get)
            heavy_candidates = {k: v for k, v in res_counts.items() if k != peptide_entity}
            heavy_entity = max(heavy_candidates, key=heavy_candidates.get) if heavy_candidates else peptide_entity

            tok_is_pep = (entity_ids == peptide_entity)
            tok_is_heavy = (entity_ids == heavy_entity)

            atom_pep = atom_to_token[:, tok_is_pep].any(dim=1)
            atom_heavy = atom_to_token[:, tok_is_heavy].any(dim=1)

            atom_heavy = atom_heavy & ~atom_pep
            atom_pep = atom_pep & present_atom
            atom_heavy = atom_heavy & present_atom

            if atom_heavy.sum() < 3 or atom_pep.sum() < 3:
                print(
                    f"Insufficient atoms for alignment in record {record_id}"
                )
                continue

            if pre_align is None or pre_rmsd is None:
                raise ValueError(
                    f"Missing precomputed masks for record '{record_id}'. "
                    "Please regenerate mask files."
                )

            valid_mask = present_atom
            valid_count = int(valid_mask.sum().item())
            if structure is None or pre_align.shape[0] != valid_count or pre_rmsd.shape[0] != valid_count:
                raise ValueError(
                    f"Precomputed mask length mismatch for record '{record_id}'. "
                    "Ensure masks were generated after removing invalid chains."
                )

            heavy_mask = torch.zeros_like(valid_mask, dtype=torch.bool, device=device)
            peptide_mask = torch.zeros_like(valid_mask, dtype=torch.bool, device=device)

            heavy_mask[valid_mask] = torch.from_numpy(pre_align.astype(bool, copy=False)).to(device=device)
            peptide_mask[valid_mask] = torch.from_numpy(pre_rmsd.astype(bool, copy=False)).to(device=device)

            if heavy_mask.sum() < 3 or peptide_mask.sum() < 3:
                raise ValueError(
                    f"Precomputed masks for record '{record_id}' contain fewer than 3 atoms. "
                    "Please regenerate masks with enough atoms for alignment and RMSD."
                )

            heavy_calc_mask[structure_idx] = heavy_mask
            peptide_calc_mask[structure_idx] = peptide_mask
            print(
                f"[debug] structure {structure_idx} entities: {res_counts} | heavy={heavy_entity} peptide={peptide_entity} | "
                "align_mode=precomputed calc_mode=precomputed"
            )
            print(
                f"Align atoms: {int(heavy_mask.sum().item())} | Calc atoms: {int(peptide_mask.sum().item())}"
            )

        return heavy_calc_mask, peptide_calc_mask

    def _write_validation_cifs(
        self,
        batch: dict[str, Tensor],
        out: dict[str, Tensor],
        base_structures: list[Optional[Structure]],
        record_ids: list[str],
        sample_metrics: Iterable[SampleMetrics],
        n_samples: int,
        batch_idx: int,
    ) -> None:
        if not getattr(self.validation_args, "write_cif_for_validation", True):
            return

        output_dir = Path(getattr(self.validation_args, "val_cif_out_dir", "validation_outputs"))
        output_dir.mkdir(parents=True, exist_ok=True)

        if (
            self.validation_writer is None
            or Path(self.validation_writer.output_dir) != output_dir
        ):
            self.validation_writer = BoltzWriter(
                data_dir=".",
                output_dir=str(output_dir),
                output_format="mmcif",
            )

        records = batch.get("record", [])
        structure_paths = batch.get("structure_path")
        atom_pad_mask = batch["atom_pad_mask"].detach().cpu().bool()
        metrics_map = {metric.sample_idx: metric for metric in sample_metrics}
        total_samples = out["sample_atom_coords"].shape[0]
        samples_per_structure = n_samples if n_samples > 0 else total_samples

        for struct_idx, record_id in enumerate(record_ids):
            if not records or struct_idx >= len(records):
                print(f"[warning] Missing record metadata for '{record_id}', skipping CIF write.")
                continue
            record = records[struct_idx]

            sample_block = out["sample_atom_coords"][
                struct_idx * samples_per_structure : (struct_idx + 1) * samples_per_structure
            ].detach().cpu()
            pad_mask = atom_pad_mask[struct_idx : struct_idx + 1]

            filenames = []
            for sample_offset in range(samples_per_structure):
                sample_idx = struct_idx * samples_per_structure + sample_offset
                metrics = metrics_map.get(sample_idx)
                if metrics is not None:
                    name = (
                        f"prediction_{record_id}_sample_{sample_idx}"
                        f"_whole{metrics.rmsd_whole:.2f}_mask{metrics.rmsd_masked:.2f}"
                    )
                else:
                    name = f"prediction_{record_id}_sample_{sample_idx}"
                filenames.append(name)

            structure_path_entry = None
            if structure_paths is not None and struct_idx < len(structure_paths):
                structure_path_entry = [structure_paths[struct_idx]]

            structure_entry = None
            if base_structures and struct_idx < len(base_structures):
                structure_entry = [base_structures[struct_idx]]

            prediction_payload = {
                "exception": False,
                "coords": sample_block,
                "masks": pad_mask,
                "structure_paths": structure_path_entry,
                "structures": structure_entry,
                "filenames": [filenames],
            }

            writer_batch = {"record": [record]}

            self.validation_writer.write_on_batch_end(
                trainer=None,
                pl_module=self,
                prediction=prediction_payload,
                batch_indices=[0],
                batch=writer_batch,
                batch_idx=batch_idx,
                dataloader_idx=0,
            )

    def _compute_sample_metrics(
        self,
        batch: dict[str, Tensor],
        out: dict[str, Tensor],
        true_coords: Tensor,
        heavy_calc_mask: Tensor,
        peptide_calc_mask: Tensor,
        n_samples: int,
    ) -> list[SampleMetrics]:
        sample_coords = out["sample_atom_coords"]
        device = sample_coords.device
        total_samples = sample_coords.shape[0]
        samples_per_structure = max(n_samples, 1)

        atom_mask_full = batch["atom_resolved_mask"].to(device=device).float()
        atom_to_token_full = batch["atom_to_token"].float().to(device=device)
        mol_type_full = batch["mol_type"].to(device=device)
        heavy_mask_full = heavy_calc_mask.to(device=device, dtype=atom_mask_full.dtype)
        peptide_mask_full = peptide_calc_mask.to(device=device, dtype=atom_mask_full.dtype)

        sample_metrics: list[SampleMetrics] = []

        for sample_idx in range(total_samples):
            struct_idx = sample_idx // samples_per_structure
            pred_sample = sample_coords[sample_idx : sample_idx + 1]
            ref_sample = true_coords[sample_idx : sample_idx + 1]

            struct_atom_mask = atom_mask_full[struct_idx : struct_idx + 1]
            struct_atom_to_token = atom_to_token_full[struct_idx : struct_idx + 1]
            struct_mol_type = mol_type_full[struct_idx : struct_idx + 1]

            whole_rmsd = float("nan")
            masked_rmsd = float("nan")

            try:
                whole_rmsd_tensor, _, _ = weighted_minimum_rmsd_single(
                    pred_sample,
                    ref_sample,
                    struct_atom_mask,
                    struct_atom_to_token,
                    struct_mol_type,
                    nucleotide_weight=self.nucleotide_rmsd_weight,
                    ligand_weight=self.ligand_rmsd_weight,
                )
                whole_rmsd = whole_rmsd_tensor.item()
            except Exception as exc:  # noqa: BLE001
                print(f"Weighted RMSD (MHC Chain) failed for sample {sample_idx}: {exc}")

            align_mask_row = heavy_mask_full[struct_idx : struct_idx + 1]
            calc_mask_row = peptide_mask_full[struct_idx : struct_idx + 1]
            try:
                if (
                    calc_mask_row.float().sum() >= 3
                    and align_mask_row.float().sum() >= 3
                ):
                    masked_rmsd_tensor, _, _ = weighted_minimum_rmsd_single(
                        pred_sample,
                        ref_sample,
                        struct_atom_mask,
                        struct_atom_to_token,
                        struct_mol_type,
                        nucleotide_weight=self.nucleotide_rmsd_weight,
                        ligand_weight=self.ligand_rmsd_weight,
                        alignment_mask=align_mask_row,
                        rmsd_mask=calc_mask_row,
                    )
                    masked_rmsd = masked_rmsd_tensor.item()
            except Exception as exc:  # noqa: BLE001
                print(f"Weighted RMSD (masked) failed for sample {sample_idx}: {exc}")

            print(f"Sample {sample_idx} weighted RMSD (MHC Chain): {whole_rmsd:.3f}Å")
            if math.isnan(masked_rmsd):
                print(f"Sample {sample_idx} weighted RMSD (masked region): nan")
            else:
                print(
                    f"Sample {sample_idx} weighted RMSD (masked region): {masked_rmsd:.3f}Å"
                )

            sample_metrics.append(
                SampleMetrics(
                    sample_idx=sample_idx,
                    rmsd_whole=whole_rmsd,
                    rmsd_masked=masked_rmsd,
                )
            )

        return sample_metrics

    def _extract_record_ids(self, batch: dict[str, Tensor], batch_idx: int) -> list[str]:
        record_ids = batch.get("record_id", None)
        if record_ids is None:
            return [f"batch_{batch_idx}_{idx}" for idx in range(batch["atom_to_token"].shape[0])]
        if isinstance(record_ids, torch.Tensor):
            return [str(r.item()) for r in record_ids]
        return [str(r) for r in record_ids]

    def debug_peptide_mask_residues(
        self,
        base_structures: list[Optional[Structure]],
        record_ids: list[str],
        peptide_calc_mask: torch.Tensor,
        precomputed_rmsd_masks: list[Optional[np.ndarray]],
    ):
        def _labels(structure: Structure, mask) -> list[str]:
            if structure is None or not hasattr(structure, "residues"):
                return []
            if isinstance(mask, torch.Tensor):
                mask = mask.detach().cpu().numpy()
            if mask is None or mask.size == 0:
                return []

            res_arr = structure.residues
            names = set(getattr(res_arr, "dtype", ()).names or [])
            if not {"atom_idx", "atom_num", "name"}.issubset(names):
                return []

            resnum_field = "res_num" if "res_num" in names else ("res_idx" if "res_idx" in names else None)
            out = []
            L = mask.shape[0]
            for i, r in enumerate(res_arr):
                s = int(r["atom_idx"]); e = s + int(r["atom_num"])
                if s < L and mask[s:e].any():
                    rname = r["name"].decode("utf-8", "ignore") if isinstance(r["name"], (bytes, bytearray)) else str(r["name"])
                    rnum = int(r[resnum_field]) if resnum_field else (i + 1)
                    out.append(f"{rname}{rnum}")
            return out

        for i, (structure, rec_id) in enumerate(zip(base_structures, record_ids)):
            if structure is None:
                continue
            pred_mask = peptide_calc_mask[i].detach().cpu().numpy()
            ref_mask  = precomputed_rmsd_masks[i] if i < len(precomputed_rmsd_masks) else None

            pred_list = _labels(structure, pred_mask)
            ref_list  = _labels(structure, ref_mask)

            pred_str = ", ".join(pred_list) if pred_list else "None"
            ref_str  = ", ".join(ref_list)  if ref_list  else "None"
            print(f"[debug-peptide-residues] record {rec_id}\n REF: {ref_str}\n PRD: {pred_str}")

    def validation_step(self, batch: dict[str, Tensor], batch_idx: int):
        # Compute the forward pass
        n_samples = self.validation_args.diffusion_samples
        try:
            out = self(
                batch,
                recycling_steps=self.validation_args.recycling_steps,
                num_sampling_steps=self.validation_args.sampling_steps,
                diffusion_samples=n_samples,
                run_confidence_sequentially=self.validation_args.run_confidence_sequentially,
            )
        except RuntimeError as e:
            if "out of memory" in str(e):
                print("| WARNING: ran out of memory, skipping batch")
                torch.cuda.empty_cache()
                gc.collect()
                return
            else:
                raise e

        try:
            record_ids = self._extract_record_ids(batch, batch_idx)

            batch_size = batch["atom_to_token"].shape[0]

            base_structures = batch.get("base_structure", None)
            if base_structures is None:
                raise ValueError("Validation batch missing 'base_structure'. Ensure dataloader supplies preloaded structures.")
            precomputed_align_masks: list[Optional[np.ndarray]]
            precomputed_rmsd_masks: list[Optional[np.ndarray]]

            if base_structures is None:
                raise ValueError("Validation requires 'base_structure' in batch to avoid disk reloads.")
            if not isinstance(base_structures, list):
                base_structures = [base_structures]
            precomputed_align_masks = [None] * len(base_structures)
            precomputed_rmsd_masks = [None] * len(base_structures)

            align_tensor = batch.get("alignment_mask")
            if isinstance(align_tensor, torch.Tensor):
                align_tensor = align_tensor.bool()
                for idx in range(min(align_tensor.shape[0], len(base_structures))):
                    precomputed_align_masks[idx] = align_tensor[idx].detach().cpu().numpy()

            rmsd_tensor = batch.get("rmsd_mask")
            if isinstance(rmsd_tensor, torch.Tensor):
                rmsd_tensor = rmsd_tensor.bool()
                for idx in range(min(rmsd_tensor.shape[0], len(base_structures))):
                    precomputed_rmsd_masks[idx] = rmsd_tensor[idx].detach().cpu().numpy()

            heavy_calc_mask, peptide_calc_mask = self._build_alignment_masks(
                batch=batch,
                record_ids=record_ids,
                base_structures=base_structures,
                batch_idx=batch_idx,
                precomputed_align_masks=precomputed_align_masks,
                precomputed_rmsd_masks=precomputed_rmsd_masks,
            )

            if getattr(self.validation_args, "debug_peptide_mask_info", False):
                self.debug_peptide_mask_residues(base_structures, record_ids, peptide_calc_mask, precomputed_rmsd_masks)
                    
            true_coords, rmsds, best_rmsds, true_coords_resolved_mask, masked_rmsds, best_masked_rmsds = self.get_true_coordinates(
                batch=batch,
                out=out,
                diffusion_samples=n_samples,
                symmetry_correction=self.validation_args.symmetry_correction,
                alignment_mask=heavy_calc_mask,
                rmsd_mask=peptide_calc_mask,
            )


            boundaries = torch.linspace(2, 22.0, 63, device=out["pdistogram"].device)
            lower = torch.tensor([1.0], device=boundaries.device)
            upper = torch.tensor([27.0], device=boundaries.device)
            exp_boundaries = torch.cat((lower, boundaries, upper))
            mid_points = ((exp_boundaries[:-1] + exp_boundaries[1:]) / 2).to(out["pdistogram"])

            preds = out["pdistogram"]
            pred_softmax = torch.softmax(preds, dim=-1)
            pred_max = pred_softmax.argmax(dim=-1)
            pred_one_hot = torch.nn.functional.one_hot(pred_max, num_classes=preds.shape[-1])
            pred_dist = (pred_one_hot * mid_points).sum(dim=-1)
            true_center = batch["disto_center"]
            true_dists = torch.cdist(true_center, true_center)

            disto_lddt_dict, disto_total_dict = factored_token_lddt_dist_loss(
                feats=batch,
                true_d=true_dists,
                pred_d=pred_dist,
            )
            batch["token_disto_mask"] = batch["token_disto_mask"]

            all_lddt_dict, all_total_dict = factored_lddt_loss(
                feats=batch,
                atom_mask=true_coords_resolved_mask,
                true_atom_coords=true_coords,
                pred_atom_coords=out["sample_atom_coords"],
                multiplicity=n_samples,
            )

            best_lddt_dict, best_total_dict = {}, {}
            best_complex_lddt_dict, best_complex_total_dict = {}, {}
            B = true_coords.shape[0] // n_samples
            if n_samples > 1:
                complex_total = 0
                complex_lddt = 0
                for key in all_lddt_dict.keys():
                    complex_lddt += all_lddt_dict[key] * all_total_dict[key]
                    complex_total += all_total_dict[key]
                complex_lddt /= complex_total + 1e-7
                best_complex_idx = complex_lddt.reshape(-1, n_samples).argmax(dim=1)
                for key in all_lddt_dict:
                    reshaped_lddt = all_lddt_dict[key].reshape(-1, n_samples)
                    reshaped_total = all_total_dict[key].reshape(-1, n_samples)
                    best_idx = reshaped_lddt.argmax(dim=1)
                    best_lddt_dict[key] = reshaped_lddt[torch.arange(B), best_idx]
                    best_total_dict[key] = reshaped_total[torch.arange(B), best_idx]
                    best_complex_lddt_dict[key] = reshaped_lddt[
                        torch.arange(B), best_complex_idx
                    ]
                    best_complex_total_dict[key] = reshaped_total[
                        torch.arange(B), best_complex_idx
                    ]
            else:
                best_lddt_dict = all_lddt_dict
                best_total_dict = all_total_dict
                best_complex_lddt_dict = all_lddt_dict
                best_complex_total_dict = all_total_dict

            def _allowed_metric(key: str) -> bool:
                low = key.lower()
                return ("dna" not in low) and ("rna" not in low) and ("ligand" not in low)

            if self.confidence_prediction and n_samples > 1:
                mae_plddt_dict, total_mae_plddt_dict = compute_plddt_mae(
                    pred_atom_coords=out["sample_atom_coords"],
                    feats=batch,
                    true_atom_coords=true_coords,
                    pred_lddt=out["plddt"],
                    true_coords_resolved_mask=true_coords_resolved_mask,
                    multiplicity=n_samples,
                )
                mae_pde_dict, total_mae_pde_dict = compute_pde_mae(
                    pred_atom_coords=out["sample_atom_coords"],
                    feats=batch,
                    true_atom_coords=true_coords,
                    pred_pde=out["pde"],
                    true_coords_resolved_mask=true_coords_resolved_mask,
                    multiplicity=n_samples,
                )
                mae_pae_dict, total_mae_pae_dict = compute_pae_mae(
                    pred_atom_coords=out["sample_atom_coords"],
                    feats=batch,
                    true_atom_coords=true_coords,
                    pred_pae=out["pae"],
                    true_coords_resolved_mask=true_coords_resolved_mask,
                    multiplicity=n_samples,
                )

                plddt = out["complex_plddt"].reshape(-1, n_samples)
                top1_idx = plddt.argmax(dim=1)
                iplddt = out["complex_iplddt"].reshape(-1, n_samples)
                iplddt_top1_idx = iplddt.argmax(dim=1)
                pde = out["complex_pde"].reshape(-1, n_samples)
                pde_top1_idx = pde.argmin(dim=1)
                ipde = out["complex_ipde"].reshape(-1, n_samples)
                ipde_top1_idx = ipde.argmin(dim=1)
                ptm = out["ptm"].reshape(-1, n_samples)
                ptm_top1_idx = ptm.argmax(dim=1)
                iptm = out["iptm"].reshape(-1, n_samples)
                iptm_top1_idx = iptm.argmax(dim=1)
                ligand_iptm = out["ligand_iptm"].reshape(-1, n_samples)
                ligand_iptm_top1_idx = ligand_iptm.argmax(dim=1)
                protein_iptm = out["protein_iptm"].reshape(-1, n_samples)
                protein_iptm_top1_idx = protein_iptm.argmax(dim=1)

                for key in all_lddt_dict:
                    if not _allowed_metric(key):
                        continue
                    reshaped_lddt = all_lddt_dict[key].reshape(-1, n_samples)
                    reshaped_total = all_total_dict[key].reshape(-1, n_samples)

                    self.top1_lddt[key].update(
                        reshaped_lddt[torch.arange(B), top1_idx],
                        reshaped_total[torch.arange(B), top1_idx],
                    )
                    self.iplddt_top1_lddt[key].update(
                        reshaped_lddt[torch.arange(B), iplddt_top1_idx],
                        reshaped_total[torch.arange(B), iplddt_top1_idx],
                    )
                    self.pde_top1_lddt[key].update(
                        reshaped_lddt[torch.arange(B), pde_top1_idx],
                        reshaped_total[torch.arange(B), pde_top1_idx],
                    )
                    self.ipde_top1_lddt[key].update(
                        reshaped_lddt[torch.arange(B), ipde_top1_idx],
                        reshaped_total[torch.arange(B), ipde_top1_idx],
                    )
                    self.ptm_top1_lddt[key].update(
                        reshaped_lddt[torch.arange(B), ptm_top1_idx],
                        reshaped_total[torch.arange(B), ptm_top1_idx],
                    )
                    self.iptm_top1_lddt[key].update(
                        reshaped_lddt[torch.arange(B), iptm_top1_idx],
                        reshaped_total[torch.arange(B), iptm_top1_idx],
                    )
                    self.ligand_iptm_top1_lddt[key].update(
                        reshaped_lddt[torch.arange(B), ligand_iptm_top1_idx],
                        reshaped_total[torch.arange(B), ligand_iptm_top1_idx],
                    )
                    self.protein_iptm_top1_lddt[key].update(
                        reshaped_lddt[torch.arange(B), protein_iptm_top1_idx],
                        reshaped_total[torch.arange(B), protein_iptm_top1_idx],
                    )

                    self.avg_lddt[key].update(all_lddt_dict[key], all_total_dict[key])
                    self.pde_mae[key].update(mae_pde_dict[key], total_mae_pde_dict[key])
                    self.pae_mae[key].update(mae_pae_dict[key], total_mae_pae_dict[key])

                for key in mae_plddt_dict:
                    if not _allowed_metric(key):
                        continue
                    self.plddt_mae[key].update(
                        mae_plddt_dict[key], total_mae_plddt_dict[key]
                    )

            def _update_metric_collection(collection, values_dict, totals_dict):
                for key in values_dict:
                    if not _allowed_metric(key):
                        continue
                    collection[key].update(values_dict[key], totals_dict[key])

            _update_metric_collection(self.lddt, best_lddt_dict, best_total_dict)
            _update_metric_collection(self.disto_lddt, disto_lddt_dict, disto_total_dict)
            _update_metric_collection(
                self.complex_lddt, best_complex_lddt_dict, best_complex_total_dict
            )
            self.rmsd.update(rmsds)
            self.best_rmsd.update(best_rmsds)

            sample_metrics = self._compute_sample_metrics(
                batch=batch,
                out=out,
                true_coords=true_coords,
                heavy_calc_mask=heavy_calc_mask,
                peptide_calc_mask=peptide_calc_mask,
                n_samples=n_samples,
            )

            whole_values = [m.rmsd_whole for m in sample_metrics if not math.isnan(m.rmsd_whole)]
            masked_values = [m.rmsd_masked for m in sample_metrics if not math.isnan(m.rmsd_masked)]

            if whole_values:
                whole_mean = torch.tensor(whole_values, device=out["sample_atom_coords"].device, dtype=torch.float32).mean()
                self.log("val/weighted_rmsd_whole", whole_mean, prog_bar=False, sync_dist=True, batch_size=batch_size)
            if masked_values:
                masked_mean = torch.tensor(masked_values, device=out["sample_atom_coords"].device, dtype=torch.float32).mean()
                self.log("val/weighted_rmsd_masked", masked_mean, prog_bar=False, sync_dist=True, batch_size=batch_size)


            self._write_validation_cifs(
                batch=batch,
                out=out,
                base_structures=base_structures,
                record_ids=record_ids,
                sample_metrics=sample_metrics,
                n_samples=n_samples,
                batch_idx=batch_idx,
            )


        except RuntimeError as e:
            if "out of memory" in str(e):
                print("| WARNING: ran out of memory, skipping batch")
                torch.cuda.empty_cache()
                gc.collect()
                return
            else:
                raise e

    def on_validation_epoch_end(self):
        avg_lddt = {}
        avg_disto_lddt = {}
        avg_complex_lddt = {}
        if self.confidence_prediction:
            avg_top1_lddt = {}
            avg_iplddt_top1_lddt = {}
            avg_pde_top1_lddt = {}
            avg_ipde_top1_lddt = {}
            avg_ptm_top1_lddt = {}
            avg_iptm_top1_lddt = {}
            avg_ligand_iptm_top1_lddt = {}
            avg_protein_iptm_top1_lddt = {}

            avg_avg_lddt = {}
            avg_mae_plddt = {}
            avg_mae_pde = {}
            avg_mae_pae = {}

        for m in const.out_types + ["pocket_ligand_protein"]:
            avg_lddt[m] = self.lddt[m].compute()
            avg_lddt[m] = 0.0 if torch.isnan(avg_lddt[m]) else avg_lddt[m].item()
            self.lddt[m].reset()
            self.log(f"val/lddt_{m}", avg_lddt[m], prog_bar=False, sync_dist=True)

            avg_disto_lddt[m] = self.disto_lddt[m].compute()
            avg_disto_lddt[m] = (
                0.0 if torch.isnan(avg_disto_lddt[m]) else avg_disto_lddt[m].item()
            )
            self.disto_lddt[m].reset()
            self.log(
                f"val/disto_lddt_{m}", avg_disto_lddt[m], prog_bar=False, sync_dist=True
            )
            avg_complex_lddt[m] = self.complex_lddt[m].compute()
            avg_complex_lddt[m] = (
                0.0 if torch.isnan(avg_complex_lddt[m]) else avg_complex_lddt[m].item()
            )
            self.complex_lddt[m].reset()
            self.log(
                f"val/complex_lddt_{m}",
                avg_complex_lddt[m],
                prog_bar=False,
                sync_dist=True,
            )
            if self.confidence_prediction:
                avg_top1_lddt[m] = self.top1_lddt[m].compute()
                avg_top1_lddt[m] = (
                    0.0 if torch.isnan(avg_top1_lddt[m]) else avg_top1_lddt[m].item()
                )
                self.top1_lddt[m].reset()
                self.log(
                    f"val/top1_lddt_{m}",
                    avg_top1_lddt[m],
                    prog_bar=False,
                    sync_dist=True,
                )
                avg_iplddt_top1_lddt[m] = self.iplddt_top1_lddt[m].compute()
                avg_iplddt_top1_lddt[m] = (
                    0.0
                    if torch.isnan(avg_iplddt_top1_lddt[m])
                    else avg_iplddt_top1_lddt[m].item()
                )
                self.iplddt_top1_lddt[m].reset()
                self.log(
                    f"val/iplddt_top1_lddt_{m}",
                    avg_iplddt_top1_lddt[m],
                    prog_bar=False,
                    sync_dist=True,
                )
                avg_pde_top1_lddt[m] = self.pde_top1_lddt[m].compute()
                avg_pde_top1_lddt[m] = (
                    0.0
                    if torch.isnan(avg_pde_top1_lddt[m])
                    else avg_pde_top1_lddt[m].item()
                )
                self.pde_top1_lddt[m].reset()
                self.log(
                    f"val/pde_top1_lddt_{m}",
                    avg_pde_top1_lddt[m],
                    prog_bar=False,
                    sync_dist=True,
                )
                avg_ipde_top1_lddt[m] = self.ipde_top1_lddt[m].compute()
                avg_ipde_top1_lddt[m] = (
                    0.0
                    if torch.isnan(avg_ipde_top1_lddt[m])
                    else avg_ipde_top1_lddt[m].item()
                )
                self.ipde_top1_lddt[m].reset()
                self.log(
                    f"val/ipde_top1_lddt_{m}",
                    avg_ipde_top1_lddt[m],
                    prog_bar=False,
                    sync_dist=True,
                )
                avg_ptm_top1_lddt[m] = self.ptm_top1_lddt[m].compute()
                avg_ptm_top1_lddt[m] = (
                    0.0
                    if torch.isnan(avg_ptm_top1_lddt[m])
                    else avg_ptm_top1_lddt[m].item()
                )
                self.ptm_top1_lddt[m].reset()
                self.log(
                    f"val/ptm_top1_lddt_{m}",
                    avg_ptm_top1_lddt[m],
                    prog_bar=False,
                    sync_dist=True,
                )
                avg_iptm_top1_lddt[m] = self.iptm_top1_lddt[m].compute()
                avg_iptm_top1_lddt[m] = (
                    0.0
                    if torch.isnan(avg_iptm_top1_lddt[m])
                    else avg_iptm_top1_lddt[m].item()
                )
                self.iptm_top1_lddt[m].reset()
                self.log(
                    f"val/iptm_top1_lddt_{m}",
                    avg_iptm_top1_lddt[m],
                    prog_bar=False,
                    sync_dist=True,
                )

                avg_ligand_iptm_top1_lddt[m] = self.ligand_iptm_top1_lddt[m].compute()
                avg_ligand_iptm_top1_lddt[m] = (
                    0.0
                    if torch.isnan(avg_ligand_iptm_top1_lddt[m])
                    else avg_ligand_iptm_top1_lddt[m].item()
                )
                self.ligand_iptm_top1_lddt[m].reset()
                self.log(
                    f"val/ligand_iptm_top1_lddt_{m}",
                    avg_ligand_iptm_top1_lddt[m],
                    prog_bar=False,
                    sync_dist=True,
                )

                avg_protein_iptm_top1_lddt[m] = self.protein_iptm_top1_lddt[m].compute()
                avg_protein_iptm_top1_lddt[m] = (
                    0.0
                    if torch.isnan(avg_protein_iptm_top1_lddt[m])
                    else avg_protein_iptm_top1_lddt[m].item()
                )
                self.protein_iptm_top1_lddt[m].reset()
                self.log(
                    f"val/protein_iptm_top1_lddt_{m}",
                    avg_protein_iptm_top1_lddt[m],
                    prog_bar=False,
                    sync_dist=True,
                )

                avg_avg_lddt[m] = self.avg_lddt[m].compute()
                avg_avg_lddt[m] = (
                    0.0 if torch.isnan(avg_avg_lddt[m]) else avg_avg_lddt[m].item()
                )
                self.avg_lddt[m].reset()
                self.log(
                    f"val/avg_lddt_{m}", avg_avg_lddt[m], prog_bar=False, sync_dist=True
                )
                avg_mae_pde[m] = self.pde_mae[m].compute().item()
                self.pde_mae[m].reset()
                self.log(
                    f"val/MAE_pde_{m}",
                    avg_mae_pde[m],
                    prog_bar=False,
                    sync_dist=True,
                )
                avg_mae_pae[m] = self.pae_mae[m].compute().item()
                self.pae_mae[m].reset()
                self.log(
                    f"val/MAE_pae_{m}",
                    avg_mae_pae[m],
                    prog_bar=False,
                    sync_dist=True,
                )

        for m in const.out_single_types:
            if self.confidence_prediction:
                avg_mae_plddt[m] = self.plddt_mae[m].compute().item()
                self.plddt_mae[m].reset()
                self.log(
                    f"val/MAE_plddt_{m}",
                    avg_mae_plddt[m],
                    prog_bar=False,
                    sync_dist=True,
                )

        overall_disto_lddt = sum(
            avg_disto_lddt[m] * w for (m, w) in const.out_types_weights.items()
        ) / sum(const.out_types_weights.values())
        self.log("val/disto_lddt", overall_disto_lddt, prog_bar=True, sync_dist=True)

        overall_lddt = sum(
            avg_lddt[m] * w for (m, w) in const.out_types_weights.items()
        ) / sum(const.out_types_weights.values())
        self.log("val/lddt", overall_lddt, prog_bar=True, sync_dist=True)

        overall_complex_lddt = sum(
            avg_complex_lddt[m] * w for (m, w) in const.out_types_weights.items()
        ) / sum(const.out_types_weights.values())
        self.log(
            "val/complex_lddt", overall_complex_lddt, prog_bar=True, sync_dist=True
        )

        if self.confidence_prediction:
            overall_top1_lddt = sum(
                avg_top1_lddt[m] * w for (m, w) in const.out_types_weights.items()
            ) / sum(const.out_types_weights.values())
            self.log("val/top1_lddt", overall_top1_lddt, prog_bar=True, sync_dist=True)

            overall_iplddt_top1_lddt = sum(
                avg_iplddt_top1_lddt[m] * w
                for (m, w) in const.out_types_weights.items()
            ) / sum(const.out_types_weights.values())
            self.log(
                "val/iplddt_top1_lddt",
                overall_iplddt_top1_lddt,
                prog_bar=True,
                sync_dist=True,
            )

            overall_pde_top1_lddt = sum(
                avg_pde_top1_lddt[m] * w for (m, w) in const.out_types_weights.items()
            ) / sum(const.out_types_weights.values())
            self.log(
                "val/pde_top1_lddt",
                overall_pde_top1_lddt,
                prog_bar=True,
                sync_dist=True,
            )

            overall_ipde_top1_lddt = sum(
                avg_ipde_top1_lddt[m] * w for (m, w) in const.out_types_weights.items()
            ) / sum(const.out_types_weights.values())
            self.log(
                "val/ipde_top1_lddt",
                overall_ipde_top1_lddt,
                prog_bar=True,
                sync_dist=True,
            )

            overall_ptm_top1_lddt = sum(
                avg_ptm_top1_lddt[m] * w for (m, w) in const.out_types_weights.items()
            ) / sum(const.out_types_weights.values())
            self.log(
                "val/ptm_top1_lddt",
                overall_ptm_top1_lddt,
                prog_bar=True,
                sync_dist=True,
            )

            overall_iptm_top1_lddt = sum(
                avg_iptm_top1_lddt[m] * w for (m, w) in const.out_types_weights.items()
            ) / sum(const.out_types_weights.values())
            self.log(
                "val/iptm_top1_lddt",
                overall_iptm_top1_lddt,
                prog_bar=True,
                sync_dist=True,
            )

            overall_avg_lddt = sum(
                avg_avg_lddt[m] * w for (m, w) in const.out_types_weights.items()
            ) / sum(const.out_types_weights.values())
            self.log("val/avg_lddt", overall_avg_lddt, prog_bar=True, sync_dist=True)

        self.log("val/rmsd", self.rmsd.compute(), prog_bar=True, sync_dist=True)
        self.rmsd.reset()

        self.log(
            "val/best_rmsd", self.best_rmsd.compute(), prog_bar=True, sync_dist=True
        )
        self.best_rmsd.reset()

    def predict_step(self, batch: Any, batch_idx: int, dataloader_idx: int = 0) -> Any:
        try:
            out = self(
                batch,
                recycling_steps=self.predict_args["recycling_steps"],
                num_sampling_steps=self.predict_args["sampling_steps"],
                diffusion_samples=self.predict_args["diffusion_samples"],
                max_parallel_samples=self.predict_args["diffusion_samples"],
                run_confidence_sequentially=True,
            )
            pred_dict = {"exception": False}
            pred_dict["masks"] = batch["atom_pad_mask"]
            pred_dict["coords"] = out["sample_atom_coords"]
            pred_dict["s"] = out["s"]
            pred_dict["z"] = out["z"]
            if self.predict_args.get("write_confidence_summary", True):
                pred_dict["confidence_score"] = (
                    4 * out["complex_plddt"]
                    + (
                        out["iptm"]
                        if not torch.allclose(
                            out["iptm"], torch.zeros_like(out["iptm"])
                        )
                        else out["ptm"]
                    )
                ) / 5
                for key in [
                    "ptm",
                    "iptm",
                    "ligand_iptm",
                    "protein_iptm",
                    "pair_chains_iptm",
                    "complex_plddt",
                    "complex_iplddt",
                    "complex_pde",
                    "complex_ipde",
                    "plddt",
                ]:
                    pred_dict[key] = out[key]
            if self.predict_args.get("write_full_pae", True):
                pred_dict["pae"] = out["pae"]
            if self.predict_args.get("write_full_pde", False):
                pred_dict["pde"] = out["pde"]
            return pred_dict

        except RuntimeError as e:  # catch out of memory exceptions
            if "out of memory" in str(e):
                print("| WARNING: ran out of memory, skipping batch")
                torch.cuda.empty_cache()
                gc.collect()
                return {"exception": True}
            else:
                raise

    def configure_optimizers(self):
        """Configure the optimizer."""

        if self.structure_prediction_training:
            parameters = [p for p in self.parameters() if p.requires_grad]
        else:
            parameters = [
                p for p in self.confidence_module.parameters() if p.requires_grad
            ] + [
                p
                for p in self.structure_module.out_token_feat_update.parameters()
                if p.requires_grad
            ]

        optimizer = torch.optim.Adam(
            parameters,
            betas=(self.training_args.adam_beta_1, self.training_args.adam_beta_2),
            eps=self.training_args.adam_eps,
            lr=self.training_args.base_lr,
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

    def on_save_checkpoint(self, checkpoint: dict[str, Any]) -> None:
        if self.use_ema:
            checkpoint["ema"] = self.ema.state_dict()

    def on_load_checkpoint(self, checkpoint: dict[str, Any]) -> None:
        if self.use_ema and "ema" in checkpoint:
            self.ema = ExponentialMovingAverage(
                parameters=self.parameters(), decay=self.ema_decay
            )
            if self.ema.compatible(checkpoint["ema"]["shadow_params"]):
                self.ema.load_state_dict(checkpoint["ema"], device=torch.device("cpu"))
            else:
                self.ema = None
                print(
                    "Warning: EMA state not loaded due to incompatible model parameters."
                )

    def on_train_start(self):
        if self.use_ema and self.ema is None:
            self.ema = ExponentialMovingAverage(
                parameters=self.parameters(), decay=self.ema_decay
            )
        elif self.use_ema:
            self.ema.to(self.device)

    def on_train_epoch_start(self) -> None:
        if self.use_ema:
            self.ema.restore(self.parameters())

    def on_train_batch_end(self, outputs, batch: Any, batch_idx: int) -> None:
        # Updates EMA parameters after optimizer.step()
        if self.use_ema:
            self.ema.update(self.parameters())

    def prepare_eval(self) -> None:
        if self.use_ema and self.ema is None:
            self.ema = ExponentialMovingAverage(
                parameters=self.parameters(), decay=self.ema_decay
            )

        if self.use_ema:
            self.ema.store(self.parameters())
            self.ema.copy_to(self.parameters())

    def on_validation_start(self):
        self.prepare_eval()

    def on_predict_start(self) -> None:
        self.prepare_eval()

    def on_test_start(self) -> None:
        self.prepare_eval()
