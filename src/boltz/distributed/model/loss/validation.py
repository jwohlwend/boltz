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


from collections import defaultdict

import torch
from torch import Tensor
from torch.distributed.tensor import DTensor

from boltz.data import const
from boltz.distributed.model.layers.atom_to_token import reconstruct_atom_to_token_global, single_repr_token_to_atom
from boltz.distributed.model.layers.elementwise_op import (
    ElementwiseOp,
    elementwise_op,
    scalar_tensor_op,
)
from boltz.distributed.model.layers.sharded_op import sharded_sum
from boltz.distributed.model.layers.shardwise_op import shardwise_sum
from boltz.distributed.model.layers.squeeze import shardwise_squeeze
from boltz.distributed.model.loss.diffusion import (
    weighted_rigid_align as dtensor_weighted_rigid_align,
)
from boltz.distributed.model.loss.triton.cdist_lddt import cdist_lddt
from boltz.distributed.model.validation.utils import gather_along_cp
from boltz.model.loss.confidence import (
    lddt_dist,
)


def clash_score(
    coords_repr: torch.Tensor,
    token_pad_mask: torch.Tensor,
    multiplicity: int,
    clash_cutoff: float,
) -> tuple[torch.Tensor, torch.Tensor]:
    """Compute per-sample token clash count and fraction from representative atom coordinates.

    Parameters
    ----------
    coords_repr : torch.Tensor
        Representative atom coordinates, shape ``[B*mul, N_tokens, 3]``.
    token_pad_mask : torch.Tensor
        Token padding mask, shape ``[B, N_tokens]``.
    multiplicity : int
        Diffusion multiplicity.
    clash_cutoff : float
        Distance cutoff for defining a clash in Angstrom.

    Returns
    -------
    tuple[torch.Tensor, torch.Tensor]
        - clash_atoms_count: clashing tokens per sample, shape ``[B*mul]``
        - clash_atoms_fraction: fraction of clashing valid tokens, shape ``[B*mul]``
    """
    _, _, clash_denom = cdist_lddt(
        pred_coords_row=coords_repr,
        pred_coords_col=coords_repr,
        true_coords_row=coords_repr,
        true_coords_col=coords_repr,
        mask_row=token_pad_mask.float(),
        mask_col=token_pad_mask.float(),
        multiplicity=multiplicity,
        cutoff=clash_cutoff,
        do_mask_diagonal=True,
        per_atom=True,
        return_denom=True,
    )
    B = token_pad_mask.shape[0]
    clash_denom_reshaped = clash_denom.reshape(B, multiplicity, -1)
    token_mask_bc = token_pad_mask[:, None, :]

    clash_atoms_count_2d = ((clash_denom_reshaped > 0) & token_mask_bc).sum(dim=2)
    clash_atoms_total = token_pad_mask.sum(dim=1).clamp(min=1)[:, None]

    clash_atoms_count = clash_atoms_count_2d.reshape(-1)
    clash_atoms_fraction = (clash_atoms_count_2d / clash_atoms_total).reshape(-1)
    return clash_atoms_count, clash_atoms_fraction


def factored_lddt_loss(
    true_atom_coords,
    pred_atom_coords,
    feats,
    atom_mask,
    multiplicity=1,
    cardinality_weighted=False,
):
    """Compute the lddt factorized into the different modalities.

    Uses triton kernel cdist_lddt to compute the lddt.

    Parameters
    ----------
    true_atom_coords : torch.Tensor
        Ground truth atom coordinates after symmetry correction, shape [B*mul, N_atoms, 3]
    pred_atom_coords : torch.Tensor
        Predicted atom coordinates, shape [B*mul, N_atoms, 3]
    feats : Dict[str, torch.Tensor]
        Input features with token-level tensors at base batch size B
    atom_mask : torch.Tensor
        Atom mask, shape [B, N_atoms] or [B*mul, N_atoms]. If [B*mul, N_atoms],
        the mask is downsampled by taking every multiplicity-th row.
    multiplicity : int
        Diffusion batch size, by default 1
    cardinality_weighted : bool
        If True, use the cardinality weighted loss, defaults to False

    Returns
    -------
    Dict[str, torch.Tensor]
        The lddt for each modality, each tensor shape [B*mul]
    Dict[str, torch.Tensor]
        The total number of pairs for each modality, each tensor shape [B*mul]

    """
    # extract necessary features
    atom_type = torch.bmm(feats["atom_to_token"].float(), feats["mol_type"].unsqueeze(-1).float()).squeeze(-1).long()

    # Use base masks and rely on cdist_lddt broadcasting across multiplicity.
    if atom_mask.shape[0] == atom_type.shape[0]:
        atom_mask_base = atom_mask
    elif atom_mask.shape[0] == atom_type.shape[0] * multiplicity:
        atom_mask_base = atom_mask[::multiplicity]  # to match with atom_type shape
    else:
        raise ValueError(
            "atom_mask batch dimension must be B or B*mul "
            f"(got {atom_mask.shape[0]} vs B={atom_type.shape[0]}, mul={multiplicity})"
        )

    input_dtype = pred_atom_coords.dtype
    compute_dtype = torch.promote_types(input_dtype, torch.float32)

    ligand_mask = (atom_type == const.chain_type_ids["NONPOLYMER"]).to(dtype=compute_dtype)
    dna_mask = (atom_type == const.chain_type_ids["DNA"]).to(dtype=compute_dtype)
    rna_mask = (atom_type == const.chain_type_ids["RNA"]).to(dtype=compute_dtype)
    protein_mask = (atom_type == const.chain_type_ids["PROTEIN"]).to(dtype=compute_dtype)

    atom_mask_base = atom_mask_base.to(dtype=compute_dtype)
    pred_atom_coords = pred_atom_coords.to(dtype=compute_dtype)
    true_atom_coords = true_atom_coords.to(dtype=compute_dtype)

    def score_and_total(mask_row, mask_col, cutoff, symmetrize=False):
        score, total = cdist_lddt(
            pred_coords_row=pred_atom_coords,
            pred_coords_col=pred_atom_coords,
            true_coords_row=true_atom_coords,
            true_coords_col=true_atom_coords,
            mask_row=mask_row,
            mask_col=mask_col,
            multiplicity=multiplicity,
            cutoff=cutoff,
            do_mask_diagonal=True,
            per_atom=False,
            return_denom=True,
        )
        if symmetrize:
            total = total * 2
        score = torch.where(total > 0, score, torch.ones_like(score))
        return score, total

    mask_dna = atom_mask_base * dna_mask
    mask_rna = atom_mask_base * rna_mask
    mask_ligand = atom_mask_base * ligand_mask
    mask_protein = atom_mask_base * protein_mask

    dna_protein_lddt, dna_protein_total = score_and_total(mask_dna, mask_protein, cutoff=30.0, symmetrize=True)
    rna_protein_lddt, rna_protein_total = score_and_total(mask_rna, mask_protein, cutoff=30.0, symmetrize=True)
    ligand_protein_lddt, ligand_protein_total = score_and_total(mask_ligand, mask_protein, cutoff=15.0, symmetrize=True)
    dna_ligand_lddt, dna_ligand_total = score_and_total(mask_dna, mask_ligand, cutoff=30.0, symmetrize=True)
    rna_ligand_lddt, rna_ligand_total = score_and_total(mask_rna, mask_ligand, cutoff=30.0, symmetrize=True)

    intra_dna_lddt, intra_dna_total = score_and_total(mask_dna, mask_dna, cutoff=30.0)
    intra_rna_lddt, intra_rna_total = score_and_total(mask_rna, mask_rna, cutoff=30.0)

    chain_id = feats["asym_id"]
    atom_chain_id = torch.bmm(feats["atom_to_token"].float(), chain_id.unsqueeze(-1).float()).squeeze(-1).long()

    chain_ids = torch.unique(atom_chain_id[atom_mask_base.bool()])

    def accumulate_chain_scores(base_mask, cutoff):
        score_sum = torch.zeros(true_atom_coords.shape[0], device=true_atom_coords.device)
        total_sum = torch.zeros_like(score_sum)
        for chain_value in chain_ids.tolist():
            chain_mask = (atom_chain_id == chain_value).float()
            mask_chain = base_mask * chain_mask
            if not torch.any(mask_chain):
                continue
            score, total = score_and_total(mask_chain, mask_chain, cutoff=cutoff)
            score_sum = score_sum + score * total
            total_sum = total_sum + total
        score = torch.where(
            total_sum > 0,
            score_sum / (total_sum + 1e-10),
            torch.ones_like(score_sum),
        )
        return score, total_sum

    intra_ligand_lddt, intra_ligand_total = accumulate_chain_scores(mask_ligand, cutoff=15.0)
    intra_protein_lddt, intra_protein_total = accumulate_chain_scores(mask_protein, cutoff=15.0)

    protein_score_sum = torch.zeros(true_atom_coords.shape[0], device=true_atom_coords.device)
    protein_total_sum = torch.zeros_like(protein_score_sum)
    chain_values = chain_ids.tolist()
    for i, chain_i in enumerate(chain_values):
        mask_i = mask_protein * (atom_chain_id == chain_i).float()
        if not torch.any(mask_i):
            continue
        for chain_j in chain_values[i + 1 :]:
            mask_j = mask_protein * (atom_chain_id == chain_j).float()
            if not torch.any(mask_j):
                continue
            score, total = score_and_total(mask_i, mask_j, cutoff=15.0, symmetrize=True)
            protein_score_sum = protein_score_sum + score * total
            protein_total_sum = protein_total_sum + total
    protein_protein_lddt = torch.where(
        protein_total_sum > 0,
        protein_score_sum / (protein_total_sum + 1e-10),
        torch.ones_like(protein_score_sum),
    )
    protein_protein_total = protein_total_sum

    lddt_dict = {
        "dna_protein": dna_protein_lddt,
        "rna_protein": rna_protein_lddt,
        "ligand_protein": ligand_protein_lddt,
        "dna_ligand": dna_ligand_lddt,
        "rna_ligand": rna_ligand_lddt,
        "intra_ligand": intra_ligand_lddt,
        "intra_dna": intra_dna_lddt,
        "intra_rna": intra_rna_lddt,
        "intra_protein": intra_protein_lddt,
        "protein_protein": protein_protein_lddt,
    }

    total_dict = {
        "dna_protein": dna_protein_total,
        "rna_protein": rna_protein_total,
        "ligand_protein": ligand_protein_total,
        "dna_ligand": dna_ligand_total,
        "rna_ligand": rna_ligand_total,
        "intra_ligand": intra_ligand_total,
        "intra_dna": intra_dna_total,
        "intra_rna": intra_rna_total,
        "intra_protein": intra_protein_total,
        "protein_protein": protein_protein_total,
    }
    if not cardinality_weighted:
        for key in total_dict:
            total_dict[key] = (total_dict[key] > 0.0).to(dtype=input_dtype)

    lddt_dict = {key: value.to(dtype=input_dtype) for key, value in lddt_dict.items()}
    total_dict = {key: value.to(dtype=input_dtype) for key, value in total_dict.items()}

    return lddt_dict, total_dict


def factored_token_lddt_dist_loss_triton(
    pred_token_coords,
    true_token_coords,
    mol_type,
    token_disto_mask,
    asym_id,
    multiplicity=1,
    cardinality_weighted=False,
    pred_d=None,
    true_d=None,
):
    """Compute the distogram lddt factorized into different modalities using cdist_lddt.

    Token-level analogue of factored_lddt_loss. When coordinates are provided for
    both sides, uses the cdist_lddt triton kernel to compute pairwise distances
    on-the-fly, avoiding O(N^2) materialization. When pre-computed distance
    matrices are provided (e.g. from a distogram prediction), uses those directly.

    Parameters
    ----------
    pred_token_coords : torch.Tensor
        Predicted token representative coordinates, shape [B*mul, N_tokens, 3].
    true_token_coords : torch.Tensor
        Ground truth token representative coordinates, shape [B*mul, N_tokens, 3].
    mol_type : torch.Tensor
        Molecule type per token, shape [B, N_tokens].
    token_disto_mask : torch.Tensor
        Token validity mask for distogram, shape [B, N_tokens].
    asym_id : torch.Tensor
        Chain (asymmetric unit) identifier per token, shape [B, N_tokens].
    multiplicity : int
        Diffusion multiplicity (B_mul = B * multiplicity), by default 1.
    cardinality_weighted : bool
        If True, return raw pair counts; if False, binarize totals. Default False.
    pred_d : torch.Tensor, optional
        Pre-computed predicted distance matrix, shape [B, N_tokens, N_tokens].
        When provided, overrides pred_token_coords for the predicted distances.
    true_d : torch.Tensor, optional
        Pre-computed true distance matrix, shape [B, N_tokens, N_tokens].
        When provided, overrides true_token_coords for the true distances.

    Returns
    -------
    dict[str, torch.Tensor]
        LDDT score per modality, each shape [B*mul].
    dict[str, torch.Tensor]
        Total (pair count or binary indicator) per modality, each shape [B*mul].

    """

    input_dtype = pred_token_coords.dtype
    compute_dtype = torch.promote_types(input_dtype, torch.float32)

    use_dists = pred_d is not None or true_d is not None

    token_mask = token_disto_mask.to(dtype=compute_dtype)
    true_token_coords = true_token_coords.to(dtype=compute_dtype)

    ligand_mask = (mol_type == const.chain_type_ids["NONPOLYMER"]).to(dtype=compute_dtype)
    dna_mask = (mol_type == const.chain_type_ids["DNA"]).to(dtype=compute_dtype)
    rna_mask = (mol_type == const.chain_type_ids["RNA"]).to(dtype=compute_dtype)
    protein_mask = (mol_type == const.chain_type_ids["PROTEIN"]).to(dtype=compute_dtype)

    nucleotide_mask = dna_mask + rna_mask

    mask_dna = token_mask * dna_mask
    mask_rna = token_mask * rna_mask
    mask_ligand = token_mask * ligand_mask
    mask_protein = token_mask * protein_mask

    if use_dists:
        pairwise_mask = token_mask[:, :, None] * token_mask[:, None, :]
        pairwise_mask = pairwise_mask * (1 - torch.eye(token_mask.shape[1], device=token_mask.device)[None]).to(
            pairwise_mask
        )
        cutoff_matrix = 15 + 15 * (1 - (1 - nucleotide_mask[:, :, None]) * (1 - nucleotide_mask[:, None, :]))
        eff_pred_d = pred_d if pred_d is not None else torch.cdist(pred_token_coords, pred_token_coords)
        eff_true_d = true_d if true_d is not None else torch.cdist(true_token_coords, true_token_coords)

        def score_and_total(mask_row, mask_col, cutoff, symmetrize=False):
            mask_2d = pairwise_mask * (mask_row[:, :, None] * mask_col[:, None, :])
            if symmetrize:
                mask_2d = mask_2d + pairwise_mask * (mask_col[:, :, None] * mask_row[:, None, :])
            # Keep the same API as the cdist path; this branch intentionally uses
            # per-pair cutoffs from cutoff_matrix instead of the scalar cutoff.
            score, total = lddt_dist(eff_pred_d, eff_true_d, mask_2d, cutoff_matrix)
            score = torch.where(total > 0, score, torch.ones_like(score))
            return score, total
    else:

        def score_and_total(mask_row, mask_col, cutoff, symmetrize=False):
            score, total = cdist_lddt(
                pred_coords_row=pred_token_coords,
                pred_coords_col=pred_token_coords,
                true_coords_row=true_token_coords,
                true_coords_col=true_token_coords,
                mask_row=mask_row,
                mask_col=mask_col,
                multiplicity=multiplicity,
                cutoff=cutoff,
                do_mask_diagonal=True,
                per_atom=False,
                return_denom=True,
            )
            if symmetrize:
                total = total * 2
            score = torch.where(total > 0, score, torch.ones_like(score))
            return score, total

    dna_protein_lddt, dna_protein_total = score_and_total(mask_dna, mask_protein, cutoff=30.0, symmetrize=True)
    rna_protein_lddt, rna_protein_total = score_and_total(mask_rna, mask_protein, cutoff=30.0, symmetrize=True)
    ligand_protein_lddt, ligand_protein_total = score_and_total(mask_ligand, mask_protein, cutoff=15.0, symmetrize=True)
    dna_ligand_lddt, dna_ligand_total = score_and_total(mask_dna, mask_ligand, cutoff=30.0, symmetrize=True)
    rna_ligand_lddt, rna_ligand_total = score_and_total(mask_rna, mask_ligand, cutoff=30.0, symmetrize=True)

    intra_dna_lddt, intra_dna_total = score_and_total(mask_dna, mask_dna, cutoff=30.0)
    intra_rna_lddt, intra_rna_total = score_and_total(mask_rna, mask_rna, cutoff=30.0)

    chain_ids = torch.unique(asym_id[token_mask.bool()])

    def accumulate_chain_scores(base_mask, cutoff):
        score_sum = torch.zeros(true_token_coords.shape[0], device=true_token_coords.device)
        total_sum = torch.zeros_like(score_sum)
        for chain_value in chain_ids.tolist():
            chain_mask = (asym_id == chain_value).float()
            mask_chain = base_mask * chain_mask
            if not torch.any(mask_chain):
                continue
            score, total = score_and_total(mask_chain, mask_chain, cutoff=cutoff)
            score_sum = score_sum + score * total
            total_sum = total_sum + total
        score = torch.where(
            total_sum > 0,
            score_sum / (total_sum + 1e-10),
            torch.ones_like(score_sum),
        )
        return score, total_sum

    intra_ligand_lddt, intra_ligand_total = accumulate_chain_scores(mask_ligand, cutoff=15.0)
    intra_protein_lddt, intra_protein_total = accumulate_chain_scores(mask_protein, cutoff=15.0)

    protein_score_sum = torch.zeros(true_token_coords.shape[0], device=true_token_coords.device)
    protein_total_sum = torch.zeros_like(protein_score_sum)
    chain_values = chain_ids.tolist()
    for i, chain_i in enumerate(chain_values):
        mask_i = mask_protein * (asym_id == chain_i).float()
        if not torch.any(mask_i):
            continue
        for chain_j in chain_values[i + 1 :]:
            mask_j = mask_protein * (asym_id == chain_j).float()
            if not torch.any(mask_j):
                continue
            score, total = score_and_total(mask_i, mask_j, cutoff=15.0, symmetrize=True)
            protein_score_sum = protein_score_sum + score * total
            protein_total_sum = protein_total_sum + total
    protein_protein_lddt = torch.where(
        protein_total_sum > 0,
        protein_score_sum / (protein_total_sum + 1e-10),
        torch.ones_like(protein_score_sum),
    )
    protein_protein_total = protein_total_sum

    lddt_dict = {
        "dna_protein": dna_protein_lddt,
        "rna_protein": rna_protein_lddt,
        "ligand_protein": ligand_protein_lddt,
        "dna_ligand": dna_ligand_lddt,
        "rna_ligand": rna_ligand_lddt,
        "intra_ligand": intra_ligand_lddt,
        "intra_dna": intra_dna_lddt,
        "intra_rna": intra_rna_lddt,
        "intra_protein": intra_protein_lddt,
        "protein_protein": protein_protein_lddt,
    }

    total_dict = {
        "dna_protein": dna_protein_total,
        "rna_protein": rna_protein_total,
        "ligand_protein": ligand_protein_total,
        "dna_ligand": dna_ligand_total,
        "rna_ligand": rna_ligand_total,
        "intra_ligand": intra_ligand_total,
        "intra_dna": intra_dna_total,
        "intra_rna": intra_rna_total,
        "intra_protein": intra_protein_total,
        "protein_protein": protein_protein_total,
    }
    if not cardinality_weighted:
        for key in total_dict:
            total_dict[key] = (total_dict[key] > 0.0).to(dtype=input_dtype)

    lddt_dict = {key: value.to(dtype=input_dtype) for key, value in lddt_dict.items()}
    total_dict = {key: value.to(dtype=input_dtype) for key, value in total_dict.items()}

    return lddt_dict, total_dict


def compute_disto_lddt(
    model,
    batch: dict[str, DTensor],
    out: dict[str, DTensor],
) -> tuple[dict[str, torch.Tensor], dict[str, torch.Tensor]]:
    """Compute distogram LDDT by gathering DTensors and calling the triton kernel.

    Distributed override of the serial Validator.compute_disto_lddt. Gathers the
    required DTensors to plain tensors, converts predicted distograms to distance
    matrices, then evaluates factored token LDDT for each (distogram, conformer) pair.

    Parameters
    ----------
    model
        The model, providing min_dist, max_dist, num_bins, num_distograms attributes.
    batch : dict[str, DTensor]
        Batch features as DTensors. Must contain:
        - "disto_coords_ensemble": [K, N_tokens, 3] or [B, K*N_tokens, 3]
        - "mol_type": [B, N_tokens]
        - "token_disto_mask": [B, N_tokens]
        - "asym_id": [B, N_tokens]
    out : dict[str, DTensor]
        Model outputs as DTensors. Must contain:
        - "pdistogram": [B, N, N, D, bins]

    Returns
    -------
    dict[str, torch.Tensor]
        LDDT score per modality, each shape [1] (min over D, mean over K).
    dict[str, torch.Tensor]
        Total per modality, each shape [1].

    """

    disto_coords_ensemble = gather_along_cp(batch["disto_coords_ensemble"])
    mol_type = gather_along_cp(batch["mol_type"])
    token_disto_mask = gather_along_cp(batch["token_disto_mask"])
    asym_id = gather_along_cp(batch["asym_id"])
    pdistogram = gather_along_cp(out["pdistogram"])

    boundaries = torch.linspace(
        model.min_dist,
        model.max_dist,
        model.num_bins - 1,
        device=pdistogram.device,
        dtype=pdistogram.dtype,
    )
    lower = torch.tensor([1.0], device=pdistogram.device, dtype=pdistogram.dtype)
    upper = torch.tensor([model.max_dist + 5.0], device=pdistogram.device, dtype=pdistogram.dtype)
    exp_boundaries = torch.cat((lower, boundaries, upper))
    mid_points = (exp_boundaries[:-1] + exp_boundaries[1:]) / 2

    if "coords" in batch:
        K = gather_along_cp(batch["coords"]).shape[1]
    elif hasattr(model, "num_conformers"):
        K = model.num_conformers
    else:
        raise ValueError("Unable to infer conformer count: expected `batch['coords']` or `model.num_conformers`.")
    true_center = disto_coords_ensemble.reshape(K, -1, 3)

    D = model.num_distograms
    device = pdistogram.device

    disto_lddt_dict = defaultdict(lambda: torch.zeros(K, D, device=device))
    disto_total_dict = defaultdict(lambda: torch.zeros(K, D, device=device))

    for i in range(D):
        preds = pdistogram[:, :, :, i]
        pred_dist_i = mid_points[preds.argmax(dim=-1)]

        for k in range(K):
            true_center_k = true_center[k].unsqueeze(0)

            lddt_dict_, total_dict_ = factored_token_lddt_dist_loss_triton(
                pred_token_coords=true_center_k,
                true_token_coords=true_center_k,
                mol_type=mol_type,
                token_disto_mask=token_disto_mask,
                asym_id=asym_id,
                pred_d=pred_dist_i,
            )

            for key in lddt_dict_:
                disto_lddt_dict[key][k, i] = lddt_dict_[key].item()
                disto_total_dict[key][k, i] = total_dict_[key].item()

    for key in disto_lddt_dict:
        disto_lddt_dict[key] = disto_lddt_dict[key].min(dim=1).values.mean(dim=0)[None]
        disto_total_dict[key] = disto_total_dict[key].min(dim=1).values.mean(dim=0)[None]

    return disto_lddt_dict, disto_total_dict


def get_lddt_metrics(
    atom_to_token_dtensor: DTensor,
    num_conformers: int,
    n_samples: int,
    true_coords: torch.Tensor,
    true_coords_resolved_mask: torch.Tensor,
    mol_type: torch.Tensor,
    asym_id: torch.Tensor,
    sample_atom_coords: torch.Tensor,
    expand_to_diffusion_samples: bool,
) -> tuple[dict[str, torch.Tensor], dict[str, torch.Tensor]]:
    """Compute factored LDDT metrics by gathering DTensors and calling the triton kernel.

    Distributed override of the serial Validator.get_lddt_metrics. Gathers
    atom_to_token on-demand (it is large) while reusing pre-gathered mol_type,
    asym_id, and sample_atom_coords.

    Parameters
    ----------
    atom_to_token_dtensor : DTensor
        Sharded atom-to-token mapping DTensor with placements
        (Shard(0), Shard(1), Replicate()).
    num_conformers : int
        Number of conformers (K) in the ensemble.
    n_samples : int
        Number of diffusion samples (multiplicity).
    true_coords : torch.Tensor
        Ground truth atom coordinates
        Shape [B*mul, K, N_atoms, 3] if ``expand_to_diffusion_samples=True``,
        else [K, N_atoms, 3].
    true_coords_resolved_mask : torch.Tensor
        Resolved atom mask
        Shape [B*mul, N_atoms] when ``expand_to_diffusion_samples=True``,
        else [N_atoms].
    mol_type : torch.Tensor
        Pre-gathered molecule type per token, shape [B, N_tokens].
    asym_id : torch.Tensor
        Pre-gathered chain ID per token, shape [B, N_tokens].
    sample_atom_coords : torch.Tensor
        Pre-gathered predicted atom coordinates, shape [B*mul, N_atoms, 3].
    expand_to_diffusion_samples : bool
        If True, true coordinates/masks are already expanded to diffusion
        samples. If False, this function repeats the non-expanded resolved mask
        across ``n_samples`` to match ``sample_atom_coords``.

    Returns
    -------
    dict[str, torch.Tensor]
        LDDT score per modality, each shape [B*mul, K].
    dict[str, torch.Tensor]
        Total per modality, each shape [B*mul, K].

    """

    atom_to_token = reconstruct_atom_to_token_global(atom_to_token_dtensor)

    if sample_atom_coords.ndim != 3:
        raise ValueError(
            f"sample_atom_coords must be rank 3 ([B*mul, N_atoms, 3]) (got ndim={sample_atom_coords.ndim})"
        )

    if expand_to_diffusion_samples:
        if true_coords.ndim != 4:
            raise ValueError(
                "true_coords must be rank 4 ([B*mul, K, N_atoms, 3]) when "
                f"expand_to_diffusion_samples=True (got ndim={true_coords.ndim})"
            )
        if true_coords_resolved_mask.ndim != 2:
            raise ValueError(
                "true_coords_resolved_mask must be rank 2 when expand_to_diffusion_samples=True "
                f"(got ndim={true_coords_resolved_mask.ndim})"
            )
        true_coords_K = true_coords.shape[1]
        true_coords_n_atoms = true_coords.shape[2]
        true_coords_batch = true_coords.shape[0]
    else:
        if true_coords.ndim != 3:
            raise ValueError(
                "true_coords must be rank 3 ([K, N_atoms, 3]) when "
                f"expand_to_diffusion_samples=False (got ndim={true_coords.ndim})"
            )
        if true_coords_resolved_mask.ndim != 1:
            raise ValueError(
                "true_coords_resolved_mask must be rank 1 when expand_to_diffusion_samples=False "
                f"(got ndim={true_coords_resolved_mask.ndim})"
            )
        true_coords_K = true_coords.shape[0]
        true_coords_n_atoms = true_coords.shape[1]
        true_coords_batch = None
        true_coords_resolved_mask = true_coords_resolved_mask.unsqueeze(0).repeat((n_samples, 1))

    if atom_to_token.shape[0] != 1:
        raise ValueError(
            "get_lddt_metrics currently expects local batch size 1 after atom_to_token reconstruction "
            f"(got atom_to_token.shape[0]={atom_to_token.shape[0]})"
        )
    if atom_to_token.shape[0] * n_samples != sample_atom_coords.shape[0]:
        raise ValueError(
            "sample_atom_coords batch must equal atom_to_token batch * n_samples "
            f"(got sample_atom_coords.shape[0]={sample_atom_coords.shape[0]}, "
            f"atom_to_token.shape[0]={atom_to_token.shape[0]}, n_samples={n_samples})"
        )

    K = num_conformers
    if true_coords_K != K:
        raise ValueError(f"true_coords conformer count ({true_coords_K}) != num_conformers ({K})")
    if true_coords_batch is not None and true_coords_batch != sample_atom_coords.shape[0]:
        raise ValueError(
            f"true_coords batch dim ({true_coords_batch}) != "
            f"sample_atom_coords batch dim ({sample_atom_coords.shape[0]})"
        )

    N_atoms = atom_to_token.shape[1]
    N_tokens = atom_to_token.shape[2]
    if mol_type.shape[-1] != N_tokens:
        raise ValueError(f"mol_type N_tokens ({mol_type.shape[-1]}) != atom_to_token N_tokens ({N_tokens})")
    if asym_id.shape[-1] != N_tokens:
        raise ValueError(f"asym_id N_tokens ({asym_id.shape[-1]}) != atom_to_token N_tokens ({N_tokens})")
    if sample_atom_coords.shape[1] != N_atoms:
        raise ValueError(
            f"sample_atom_coords N_atoms ({sample_atom_coords.shape[1]}) != atom_to_token N_atoms ({N_atoms})"
        )
    if true_coords_resolved_mask.shape[1] != N_atoms:
        raise ValueError(
            f"true_coords_resolved_mask N_atoms ({true_coords_resolved_mask.shape[1]}) != "
            f"atom_to_token N_atoms ({N_atoms})"
        )
    if true_coords_n_atoms != N_atoms:
        raise ValueError(f"true_coords N_atoms ({true_coords_n_atoms}) != atom_to_token N_atoms ({N_atoms})")

    feats = {
        "atom_to_token": atom_to_token,
        "mol_type": mol_type,
        "asym_id": asym_id,
    }

    all_lddt_dict = defaultdict(list)
    all_total_dict = defaultdict(list)
    for ensemble_idx in range(K):
        if expand_to_diffusion_samples:
            true_coords_k = true_coords[:, ensemble_idx]
        else:
            true_coords_k = true_coords[ensemble_idx].unsqueeze(0).repeat((n_samples, 1, 1))

        lddt_dict_k, total_dict_k = factored_lddt_loss(
            true_atom_coords=true_coords_k,
            pred_atom_coords=sample_atom_coords,
            feats=feats,
            atom_mask=true_coords_resolved_mask,
            multiplicity=n_samples,
        )
        for key in lddt_dict_k:
            all_lddt_dict[key].append(lddt_dict_k[key])
            all_total_dict[key].append(total_dict_k[key])

    for key in all_lddt_dict:
        all_lddt_dict[key] = torch.stack(all_lddt_dict[key], dim=1)
        all_total_dict[key] = torch.stack(all_total_dict[key], dim=1)

    return dict(all_lddt_dict), dict(all_total_dict)


def weighted_minimum_rmsd_single(
    pred_atom_coords: DTensor,
    atom_coords: DTensor,
    atom_mask: DTensor,
    atom_to_token: DTensor,
    mol_type: DTensor,
    nucleotide_weight: float = 5.0,
    ligand_weight: float = 10.0,
) -> tuple[DTensor, DTensor, DTensor]:
    """Compute rmsd of the aligned atom coordinates using DTensor operations.

    This is the distributed version that operates on DTensors with placements
    (Shard(0), Shard(1), Replicate()) for coords and (Shard(0), Shard(1), Replicate())
    for 2D features.

    Parameters
    ----------
    pred_atom_coords : DTensor
        Predicted atom coordinates with shape (B, N_atoms, 3).
        Placements: (Shard(0), Shard(1), Replicate())
    atom_coords : DTensor
        Ground truth atom coordinates with shape (B, N_atoms, 3).
        Placements: (Shard(0), Shard(1), Replicate())
    atom_mask : DTensor
        Resolved atom mask with shape (B, N_atoms).
        Placements: (Shard(0), Shard(1), Replicate())
    atom_to_token : DTensor
        Atom to token mapping with shape (B, N_tokens, N_atoms).
        Placements: (Shard(0), Shard(1), Replicate())
    mol_type : DTensor
        Molecule type per token with shape (B, N_tokens).
        Placements: (Shard(0), Shard(1), Replicate())
    nucleotide_weight : float
        Weight for nucleotide atoms in RMSD computation.
    ligand_weight : float
        Weight for ligand atoms in RMSD computation.

    Returns
    -------
    tuple[DTensor, DTensor, DTensor]
        - rmsd: The RMSD value with shape (B,). Placements: (Shard(0), Replicate(), Replicate())
        - atom_coords_aligned: The aligned coordinates with shape (B, N_atoms, 3).
          Placements: (Shard(0), Shard(1), Replicate())
        - align_weights: The alignment weights with shape (B, N_atoms).
          Placements: (Shard(0), Shard(1), Replicate())

    """
    # Validate inputs are DTensors
    if not isinstance(pred_atom_coords, DTensor):
        raise TypeError(f"pred_atom_coords must be DTensor, got {type(pred_atom_coords)}")
    if not isinstance(atom_coords, DTensor):
        raise TypeError(f"atom_coords must be DTensor, got {type(atom_coords)}")
    if not isinstance(atom_mask, DTensor):
        raise TypeError(f"atom_mask must be DTensor, got {type(atom_mask)}")
    if not isinstance(atom_to_token, DTensor):
        raise TypeError(f"atom_to_token must be DTensor, got {type(atom_to_token)}")
    if not isinstance(mol_type, DTensor):
        raise TypeError(f"mol_type must be DTensor, got {type(mol_type)}")

    device_mesh = pred_atom_coords.device_mesh

    # Convert dtypes as needed
    dtype = pred_atom_coords.dtype

    # Compute atom_type by mapping mol_type (token-level) to atom-level
    # atom_type has shape (B, N_atoms) - placements: (Shard(0), Shard(1), Replicate())
    atom_type = single_repr_token_to_atom(
        mol_type.to(dtype=dtype).unsqueeze(-1),  # (B, N_tokens, 1)
        atom_to_token.to(dtype=dtype),  # (B, N_tokens, N_atoms)
    )  # (B, N_atoms, 1)
    atom_type = shardwise_squeeze(atom_type, dim=-1)  # (B, N_atoms)

    # Compute nucleotide mask: is_DNA OR is_RNA
    is_dna = scalar_tensor_op(float(const.chain_type_ids["DNA"]), atom_type, ElementwiseOp.EQUAL)
    is_rna = scalar_tensor_op(float(const.chain_type_ids["RNA"]), atom_type, ElementwiseOp.EQUAL)
    is_nucleotide = elementwise_op(is_dna, is_rna, ElementwiseOp.SUM)

    # Compute ligand mask
    is_ligand = scalar_tensor_op(float(const.chain_type_ids["NONPOLYMER"]), atom_type, ElementwiseOp.EQUAL)

    # Compute weighted contributions
    nucleotide_contribution = scalar_tensor_op(nucleotide_weight, is_nucleotide, ElementwiseOp.PROD)
    ligand_contribution = scalar_tensor_op(ligand_weight, is_ligand, ElementwiseOp.PROD)

    # align_weights = 1 + nucleotide_weight * is_nucleotide + ligand_weight * is_ligand
    align_weights = scalar_tensor_op(
        1.0,
        elementwise_op(nucleotide_contribution, ligand_contribution, ElementwiseOp.SUM),
        ElementwiseOp.SUM,
    )

    # Ensure atom_mask has correct placements
    atom_mask_float = DTensor.from_local(
        atom_mask.to_local().to(dtype=dtype),
        device_mesh,
        atom_mask.placements,
        shape=atom_mask.shape,
        stride=atom_mask.stride(),
    )

    # Perform weighted rigid alignment
    with torch.no_grad():
        atom_coords_aligned_ground_truth = dtensor_weighted_rigid_align(
            atom_coords.to(dtype=dtype),
            pred_atom_coords.to(dtype=dtype),
            align_weights.to(dtype=dtype),
            mask=atom_mask_float,
        )

    # Compute MSE loss: ((pred - aligned_true) ** 2).sum(dim=-1)
    diff = elementwise_op(pred_atom_coords, atom_coords_aligned_ground_truth, ElementwiseOp.SUB)
    diff_sq = scalar_tensor_op(2.0, diff, ElementwiseOp.POW)
    mse_loss = shardwise_sum(diff_sq, dim=-1)  # sum over xyz -> (B, N_atoms)

    # Compute weighted MSE: mse_loss * align_weights * atom_mask
    weighted_mse = elementwise_op(mse_loss, align_weights, ElementwiseOp.PROD)
    weighted_mse = elementwise_op(weighted_mse, atom_mask_float, ElementwiseOp.PROD)

    # Compute denominator: align_weights * atom_mask
    denom = elementwise_op(align_weights, atom_mask_float, ElementwiseOp.PROD)

    # Reduce along atom dimension
    weighted_mse_sum = sharded_sum(weighted_mse, dim=-1)  # (B,)
    denom_sum = sharded_sum(denom, dim=-1)  # (B,)

    # rmsd = sqrt(weighted_mse_sum / denom_sum)
    ratio = elementwise_op(weighted_mse_sum, denom_sum, ElementwiseOp.DIV)
    rmsd = scalar_tensor_op(0.5, ratio, ElementwiseOp.POW)  # sqrt via x^0.5

    return rmsd, atom_coords_aligned_ground_truth, align_weights


def compute_plddt_mae_triton(
    pred_atom_coords: Tensor,
    feats: dict[str, Tensor],
    true_atom_coords: Tensor,
    pred_lddt: Tensor,
    true_coords_resolved_mask: Tensor,
    multiplicity: int = 1,
) -> tuple[dict[str, Tensor], dict[str, Tensor]]:
    """Compute pLDDT MAE using triton ``cdist_lddt`` in rectangular mode.

    Uses ``cdist_lddt`` with rectangular inputs (N_token rows x N_R_set cols)
    and factored masks instead of materialising full distance matrices and
    pair masks.

    Parameters
    ----------
    pred_atom_coords : Tensor
        Predicted atom coordinates, shape ``[B*mul, N_atom, 3]``.
    feats : dict[str, Tensor]
        Feature dict with ``token_to_rep_atom``, ``r_set_to_rep_atom``,
        ``atom_to_token``, ``mol_type``.  All at base batch B.
    true_atom_coords : Tensor
        Ground-truth atom coordinates, shape ``[B*mul, N_atom, 3]``.
    pred_lddt : Tensor
        Predicted per-token lDDT, shape ``[B*mul, N_token]``.
    true_coords_resolved_mask : Tensor
        Per-atom resolved mask, shape ``[B*mul, N_atom]``.
    multiplicity : int
        Diffusion sample count (B_mul = B * multiplicity).

    Returns
    -------
    mae_plddt_dict : dict[str, Tensor]
        Per-modality MAE, each a scalar.
    total_dict : dict[str, Tensor]
        Per-modality total weight, each a scalar.
    """
    token_to_rep_atom = feats["token_to_rep_atom"].float()
    r_set_to_rep_atom = feats["r_set_to_rep_atom"].float()
    atom_to_token = feats["atom_to_token"].float()
    mol_type = feats["mol_type"]

    if multiplicity > 1:
        t2r_expanded = token_to_rep_atom.repeat_interleave(multiplicity, 0)
        r2r_expanded = r_set_to_rep_atom.repeat_interleave(multiplicity, 0)
    else:
        t2r_expanded = token_to_rep_atom
        r2r_expanded = r_set_to_rep_atom

    pred_token_coords = torch.bmm(t2r_expanded, pred_atom_coords)
    true_token_coords = torch.bmm(t2r_expanded, true_atom_coords)
    pred_R_coords = torch.bmm(r2r_expanded, pred_atom_coords)
    true_R_coords = torch.bmm(r2r_expanded, true_atom_coords)

    # Masks at B*mul level so per-sample mask variation is preserved.
    resolved = true_coords_resolved_mask.float()
    mask_row = torch.bmm(t2r_expanded, resolved.unsqueeze(-1)).squeeze(-1)
    mask_col = torch.bmm(r2r_expanded, resolved.unsqueeze(-1)).squeeze(-1)

    # Per-column cutoff based on nucleotide type (base batch B)
    is_nucleotide_token = (mol_type == const.chain_type_ids["DNA"]).float() + (
        mol_type == const.chain_type_ids["RNA"]
    ).float()
    is_nucleotide_atom = torch.bmm(atom_to_token, is_nucleotide_token.unsqueeze(-1)).squeeze(-1)
    is_nucleotide_R = torch.bmm(r_set_to_rep_atom, is_nucleotide_atom.unsqueeze(-1)).squeeze(-1)
    cutoff_col = 15.0 + 15.0 * is_nucleotide_R

    atom_indices_row = token_to_rep_atom.argmax(dim=-1)
    atom_indices_col = r_set_to_rep_atom.argmax(dim=-1)

    target_lddt, mask_no_match = cdist_lddt(
        pred_coords_row=pred_token_coords,
        pred_coords_col=pred_R_coords,
        true_coords_row=true_token_coords,
        true_coords_col=true_R_coords,
        mask_row=mask_row,
        mask_col=mask_col,
        multiplicity=multiplicity,
        atom_indices_row=atom_indices_row,
        atom_indices_col=atom_indices_col,
        cutoff_col=cutoff_col,
        do_mask_diagonal=True,
        per_atom=True,
    )

    atom_mask = mask_row
    if multiplicity > 1:
        token_type = mol_type.repeat_interleave(multiplicity, 0)
    else:
        token_type = mol_type

    protein_mask = (token_type == const.chain_type_ids["PROTEIN"]).float() * atom_mask * mask_no_match
    ligand_mask = (token_type == const.chain_type_ids["NONPOLYMER"]).float() * atom_mask * mask_no_match
    dna_mask = (token_type == const.chain_type_ids["DNA"]).float() * atom_mask * mask_no_match
    rna_mask = (token_type == const.chain_type_ids["RNA"]).float() * atom_mask * mask_no_match

    abs_err = torch.abs(target_lddt - pred_lddt)

    def _mae_and_total(mask):
        total = torch.sum(mask)
        mae = torch.sum(abs_err * mask) / (total + 1e-5)
        return mae, total

    protein_mae, protein_total = _mae_and_total(protein_mask)
    ligand_mae, ligand_total = _mae_and_total(ligand_mask)
    dna_mae, dna_total = _mae_and_total(dna_mask)
    rna_mae, rna_total = _mae_and_total(rna_mask)

    mae_plddt_dict = {
        "protein": protein_mae,
        "ligand": ligand_mae,
        "dna": dna_mae,
        "rna": rna_mae,
    }
    total_dict = {
        "protein": protein_total,
        "ligand": ligand_total,
        "dna": dna_total,
        "rna": rna_total,
    }
    return mae_plddt_dict, total_dict
