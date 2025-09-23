import gc
import random
from typing import Any, Optional

import torch
import torch._dynamo
from pytorch_lightning import LightningModule
from torch import Tensor, nn
from torchmetrics import MeanMetric

import boltz.model.layers.initialize as init
from boltz.data import const
from boltz.data.feature.symmetry import (
    minimum_lddt_symmetry_coords,
    minimum_symmetry_coords,
)
from boltz.model.loss.confidence import confidence_loss
from boltz.model.loss.distogram import distogram_loss
from boltz.model.loss.validation import (
    compute_pae_mae,
    compute_pde_mae,
    compute_plddt_mae,
    factored_lddt_loss,
    factored_token_lddt_dist_loss,
    weighted_minimum_rmsd,
)
from boltz.model.modules.confidence import ConfidenceModule
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

    def forward(
        self,
        feats: dict[str, Tensor],
        recycling_steps: int = 0,
        num_sampling_steps: Optional[int] = None,
        multiplicity_diffusion_train: int = 1,
        diffusion_samples: int = 1,
        max_parallel_samples: Optional[int] = None,
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
    ):
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

        return true_coords, rmsds, best_rmsds, true_coords_resolved_mask

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
            true_coords, _, _, true_coords_resolved_mask = self.get_true_coordinates(
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
                "loss": torch.tensor(0.0).to(batch["token_index"].device),
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

        except RuntimeError as e:  # catch out of memory exceptions
            if "out of memory" in str(e):
                print("| WARNING: ran out of memory, skipping batch")
                torch.cuda.empty_cache()
                gc.collect()
                return
            else:
                raise e

        try:
            # Compute distogram LDDT
            boundaries = torch.linspace(2, 22.0, 63)
            lower = torch.tensor([1.0])
            upper = torch.tensor([22.0 + 5.0])
            exp_boundaries = torch.cat((lower, boundaries, upper))
            mid_points = ((exp_boundaries[:-1] + exp_boundaries[1:]) / 2).to(
                out["pdistogram"]
            )

            # Compute predicted dists
            preds = out["pdistogram"]
            pred_softmax = torch.softmax(preds, dim=-1)
            pred_softmax = pred_softmax.argmax(dim=-1)
            pred_softmax = torch.nn.functional.one_hot(
                pred_softmax, num_classes=preds.shape[-1]
            )
            pred_dist = (pred_softmax * mid_points).sum(dim=-1)
            true_center = batch["disto_center"]
            true_dists = torch.cdist(true_center, true_center)

            # Compute lddt's
            batch["token_disto_mask"] = batch["token_disto_mask"]
            disto_lddt_dict, disto_total_dict = factored_token_lddt_dist_loss(
                feats=batch,
                true_d=true_dists,
                pred_d=pred_dist,
            )

            true_coords, rmsds, best_rmsds, true_coords_resolved_mask = (
                self.get_true_coordinates(
                    batch=batch,
                    out=out,
                    diffusion_samples=n_samples,
                    symmetry_correction=self.validation_args.symmetry_correction,
                )
            )

            all_lddt_dict, all_total_dict = factored_lddt_loss(
                feats=batch,
                atom_mask=true_coords_resolved_mask,
                true_atom_coords=true_coords,
                pred_atom_coords=out["sample_atom_coords"],
                multiplicity=n_samples,
            )
        except RuntimeError as e:  # catch out of memory exceptions
            if "out of memory" in str(e):
                print("| WARNING: ran out of memory, skipping batch")
                torch.cuda.empty_cache()
                gc.collect()
                return
            else:
                raise e
        # if the multiplicity used is > 1 then we take the best lddt of the different samples
        # AF3 combines this with the confidence based filtering
        best_lddt_dict, best_total_dict = {}, {}
        best_complex_lddt_dict, best_complex_total_dict = {}, {}
        B = true_coords.shape[0] // n_samples
        if n_samples > 1:
            # NOTE: we can change the way we aggregate the lddt
            complex_total = 0
            complex_lddt = 0
            for key in all_lddt_dict.keys():
                complex_lddt += all_lddt_dict[key] * all_total_dict[key]
                complex_total += all_total_dict[key]
            complex_lddt /= complex_total + 1e-7
            best_complex_idx = complex_lddt.reshape(-1, n_samples).argmax(dim=1)
            for key in all_lddt_dict:
                best_idx = all_lddt_dict[key].reshape(-1, n_samples).argmax(dim=1)
                best_lddt_dict[key] = all_lddt_dict[key].reshape(-1, n_samples)[
                    torch.arange(B), best_idx
                ]
                best_total_dict[key] = all_total_dict[key].reshape(-1, n_samples)[
                    torch.arange(B), best_idx
                ]
                best_complex_lddt_dict[key] = all_lddt_dict[key].reshape(-1, n_samples)[
                    torch.arange(B), best_complex_idx
                ]
                best_complex_total_dict[key] = all_total_dict[key].reshape(
                    -1, n_samples
                )[torch.arange(B), best_complex_idx]
        else:
            best_lddt_dict = all_lddt_dict
            best_total_dict = all_total_dict
            best_complex_lddt_dict = all_lddt_dict
            best_complex_total_dict = all_total_dict

        # Filtering based on confidence
        if self.confidence_prediction and n_samples > 1:
            # note: for now we don't have pae predictions so have to use pLDDT instead of pTM
            # also, while AF3 differentiates the best prediction per confidence type we are currently not doing it
            # consider this in the future as well as weighing the different pLLDT types before aggregation
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
                top1_lddt = all_lddt_dict[key].reshape(-1, n_samples)[
                    torch.arange(B), top1_idx
                ]
                top1_total = all_total_dict[key].reshape(-1, n_samples)[
                    torch.arange(B), top1_idx
                ]
                iplddt_top1_lddt = all_lddt_dict[key].reshape(-1, n_samples)[
                    torch.arange(B), iplddt_top1_idx
                ]
                iplddt_top1_total = all_total_dict[key].reshape(-1, n_samples)[
                    torch.arange(B), iplddt_top1_idx
                ]
                pde_top1_lddt = all_lddt_dict[key].reshape(-1, n_samples)[
                    torch.arange(B), pde_top1_idx
                ]
                pde_top1_total = all_total_dict[key].reshape(-1, n_samples)[
                    torch.arange(B), pde_top1_idx
                ]
                ipde_top1_lddt = all_lddt_dict[key].reshape(-1, n_samples)[
                    torch.arange(B), ipde_top1_idx
                ]
                ipde_top1_total = all_total_dict[key].reshape(-1, n_samples)[
                    torch.arange(B), ipde_top1_idx
                ]
                ptm_top1_lddt = all_lddt_dict[key].reshape(-1, n_samples)[
                    torch.arange(B), ptm_top1_idx
                ]
                ptm_top1_total = all_total_dict[key].reshape(-1, n_samples)[
                    torch.arange(B), ptm_top1_idx
                ]
                iptm_top1_lddt = all_lddt_dict[key].reshape(-1, n_samples)[
                    torch.arange(B), iptm_top1_idx
                ]
                iptm_top1_total = all_total_dict[key].reshape(-1, n_samples)[
                    torch.arange(B), iptm_top1_idx
                ]
                ligand_iptm_top1_lddt = all_lddt_dict[key].reshape(-1, n_samples)[
                    torch.arange(B), ligand_iptm_top1_idx
                ]
                ligand_iptm_top1_total = all_total_dict[key].reshape(-1, n_samples)[
                    torch.arange(B), ligand_iptm_top1_idx
                ]
                protein_iptm_top1_lddt = all_lddt_dict[key].reshape(-1, n_samples)[
                    torch.arange(B), protein_iptm_top1_idx
                ]
                protein_iptm_top1_total = all_total_dict[key].reshape(-1, n_samples)[
                    torch.arange(B), protein_iptm_top1_idx
                ]

                self.top1_lddt[key].update(top1_lddt, top1_total)
                self.iplddt_top1_lddt[key].update(iplddt_top1_lddt, iplddt_top1_total)
                self.pde_top1_lddt[key].update(pde_top1_lddt, pde_top1_total)
                self.ipde_top1_lddt[key].update(ipde_top1_lddt, ipde_top1_total)
                self.ptm_top1_lddt[key].update(ptm_top1_lddt, ptm_top1_total)
                self.iptm_top1_lddt[key].update(iptm_top1_lddt, iptm_top1_total)
                self.ligand_iptm_top1_lddt[key].update(
                    ligand_iptm_top1_lddt, ligand_iptm_top1_total
                )
                self.protein_iptm_top1_lddt[key].update(
                    protein_iptm_top1_lddt, protein_iptm_top1_total
                )

                self.avg_lddt[key].update(all_lddt_dict[key], all_total_dict[key])
                self.pde_mae[key].update(mae_pde_dict[key], total_mae_pde_dict[key])
                self.pae_mae[key].update(mae_pae_dict[key], total_mae_pae_dict[key])

            for key in mae_plddt_dict:
                self.plddt_mae[key].update(
                    mae_plddt_dict[key], total_mae_plddt_dict[key]
                )

        for m in const.out_types:
            if m == "ligand_protein":
                if torch.any(
                    batch["pocket_feature"][
                        :, :, const.pocket_contact_info["POCKET"]
                    ].bool()
                ):
                    self.lddt["pocket_ligand_protein"].update(
                        best_lddt_dict[m], best_total_dict[m]
                    )
                    self.disto_lddt["pocket_ligand_protein"].update(
                        disto_lddt_dict[m], disto_total_dict[m]
                    )
                    self.complex_lddt["pocket_ligand_protein"].update(
                        best_complex_lddt_dict[m], best_complex_total_dict[m]
                    )
                else:
                    self.lddt["ligand_protein"].update(
                        best_lddt_dict[m], best_total_dict[m]
                    )
                    self.disto_lddt["ligand_protein"].update(
                        disto_lddt_dict[m], disto_total_dict[m]
                    )
                    self.complex_lddt["ligand_protein"].update(
                        best_complex_lddt_dict[m], best_complex_total_dict[m]
                    )
            else:
                self.lddt[m].update(best_lddt_dict[m], best_total_dict[m])
                self.disto_lddt[m].update(disto_lddt_dict[m], disto_total_dict[m])
                self.complex_lddt[m].update(
                    best_complex_lddt_dict[m], best_complex_total_dict[m]
                )
        self.rmsd.update(rmsds)
        self.best_rmsd.update(best_rmsds)

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
