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

from __future__ import annotations

from itertools import chain
from numbers import Integral

import torch
from torch.distributed.tensor import DTensor, Replicate, Shard, distribute_tensor

from boltz.distributed.model.loss.triton.cdist_lddt import cdist_lddt


def minimum_lddt_symmetry_coords(
    coords: DTensor,
    feats: dict,
    index_batch_local: int,
    i_batch_multiplicity_local: int = 0,
):
    """Find coordinates with best lDDT under symmetry transformations (Boltz-2 semantics).

    This function handles the distributed case where:
    - coords is a DTensor sharded across (DP, CP_axis_0, CP_axis_1) device mesh
    - symmetry features are plain Tensors with only the local batch for this DP rank

    Unlike the Boltz-1 version, this does NOT perform RMSD alignment (matching
    the Boltz-2 serial ``minimum_lddt_symmetry_coords`` from ``boltz.data.mol``).

    Local Batch Indexing Semantics
    ------------------------------
    In distributed training with DP (data parallelism) and CP (context parallelism):

    1. **coords (DTensor)**: Has global shape [B_global * multiplicity, N_atoms_padded, 3]
       - Sharded along batch dim (DP axis) and atom dim (CP axes)
       - Contains interspersed padding from pad_and_scatter_atom_features_dtensor
       - After CP gather: local shape is [B_local * multiplicity, N_atoms_padded, 3]
       - B_local = B_global / num_DP_ranks (each DP rank owns different samples)
       - Use ``i_batch_multiplicity_local`` to index: coords_local[i_batch_multiplicity_local]

    2. **symmetry features (plain Tensors and nested iterables broadcasted among all ranks)**:
       - The tensors have shape [B_local, N_atoms_no_pad, ...]
       - NOT sharded - each DP rank only has its own local batch
       - NO interspersed padding - contiguous atoms without CP shard padding
       - Use ``index_batch_local`` to index: feats["all_coords"][index_batch_local]

    3. **atom_pad_mask (DTensor in feats)**: Has global shape [B_global, N_atoms_padded]
       - Indicates valid atoms (True/1.0) vs padding (False/0.0)
       - Used to remove interspersed padding from coords for symmetry resolution
       - Then re-add padding before returning as DTensor

    Index Relationship (with multiplicity M):
    - index_batch_local = i_batch_multiplicity_local // M
    - i_batch_multiplicity_local = index_batch_local * M + rep (where rep in [0, M))

    Parameters
    ----------
    coords : DTensor
        Predicted coordinates with shape [B_global * multiplicity, N_atoms_padded, 3].
        The batch dimension includes diffusion_samples multiplicity and is sharded along DP axis.
        Contains interspersed padding from distributed featurization.
    feats : dict
        Dictionary containing (symmetry features are per-DP-rank local tensors):
        - all_coords: Tensor [B_local, N_all, 3]
        - all_resolved_mask: Tensor [B_local, N_all]
        - crop_to_all_atom_map: Tensor [B_local, N_crop]
        - chain_swaps: List[B_local] of swap combinations
        - amino_acids_symmetries: List[B_local] of symmetry groups
        - ligand_symmetries: List[B_local] of symmetry groups
        - atom_pad_mask: DTensor [B_global, N_atoms_padded] - atom padding mask (required)
    index_batch_local : int
        Local batch index into symmetry features (range: [0, B_local)).
    i_batch_multiplicity_local : int
        Local index into coords (range: [0, B_local * multiplicity)).

    Returns
    -------
    tuple[DTensor, DTensor]
        (true_coords_dtensor, true_resolved_mask_dtensor) as DTensors
        with the same placements as the input coords.

    """
    if not isinstance(coords, DTensor):
        raise TypeError(f"coords must be a DTensor, got {type(coords).__name__}.")

    if "atom_pad_mask" not in feats:
        raise KeyError("feats must contain 'atom_pad_mask' key for handling interspersed padding.")
    atom_pad_mask = feats["atom_pad_mask"]
    if not isinstance(atom_pad_mask, DTensor):
        raise TypeError(f"feats['atom_pad_mask'] must be a DTensor, got {type(atom_pad_mask).__name__}.")

    coords_mesh = coords.device_mesh
    coords_placements = coords.placements
    if coords_placements != (Shard(0), Shard(1), Replicate()):
        raise ValueError(f"Expected coords placements (Shard(0), Shard(1), Replicate()), got {coords_placements}")
    if atom_pad_mask.device_mesh != coords_mesh:
        raise ValueError("atom_pad_mask.device_mesh does not match coords.device_mesh")
    if atom_pad_mask.placements != coords_placements:
        raise ValueError(f"Expected atom_pad_mask placements {coords_placements}, got {atom_pad_mask.placements}")

    all_coords = feats["all_coords"]
    all_resolved_mask = feats["all_resolved_mask"]
    crop_to_all_atom_map = feats["crop_to_all_atom_map"]

    # CollateDTensor collates NON_SHARDED_FEATURES_V2 as Python lists of
    # tensors (one per sample, no batch dim). Stack to add the batch dim.
    if isinstance(all_coords, list):
        all_coords = torch.stack(all_coords, dim=0)
    if isinstance(all_resolved_mask, list):
        all_resolved_mask = torch.stack(all_resolved_mask, dim=0)
    if isinstance(crop_to_all_atom_map, list):
        crop_to_all_atom_map = torch.stack(crop_to_all_atom_map, dim=0)

    for key, val in [
        ("all_coords", all_coords),
        ("all_resolved_mask", all_resolved_mask),
        ("crop_to_all_atom_map", crop_to_all_atom_map),
    ]:
        if not isinstance(val, torch.Tensor):
            raise TypeError(f"feats['{key}'] must be a plain torch.Tensor, got {type(val).__name__}.")

    if coords.ndim != 3 or coords.shape[2] != 3:
        raise ValueError("coords must have shape [B, N, 3].")
    if all_coords.ndim != 3 or all_coords.shape[2] != 3:
        raise ValueError("feats['all_coords'] must have shape [B, N_all, 3].")

    chain_swaps_all = feats.get("chain_swaps")
    if not isinstance(chain_swaps_all, (list, tuple)):
        raise TypeError("feats['chain_swaps'] must be a list/tuple of swap combinations.")
    chain_swaps = chain_swaps_all[index_batch_local]
    if not isinstance(chain_swaps, (list, tuple)):
        raise TypeError("chain_swaps must be a list/tuple of swap combinations.")
    if not all(isinstance(combo, (list, tuple)) for combo in chain_swaps):
        raise TypeError("chain_swaps entries must be list/tuple of swaps.")
    if not all(
        isinstance(swap, (list, tuple)) and len(swap) == 6 and all(isinstance(v, Integral) for v in swap)
        for swap in chain.from_iterable(chain_swaps)
    ):
        raise ValueError("chain_swaps swaps must be 6-int tuples: (start1, end1, start2, end2, chainidx1, chainidx2).")

    amino_acids_symmetries_all = feats.get("amino_acids_symmetries")
    if not isinstance(amino_acids_symmetries_all, (list, tuple)):
        raise TypeError("feats['amino_acids_symmetries'] must be a list/tuple of symmetry groups.")
    amino_acids_symmetries = amino_acids_symmetries_all[index_batch_local]

    ligand_symmetries_all = feats.get("ligand_symmetries")
    if not isinstance(ligand_symmetries_all, (list, tuple)):
        raise TypeError("feats['ligand_symmetries'] must be a list/tuple of symmetry groups.")
    ligand_symmetries = ligand_symmetries_all[index_batch_local]

    # --- Validate structure: sym_groups[residue][sym_op] = [(i,j), ...] ---
    for sym_name, sym_groups in [
        ("amino_acids_symmetries", amino_acids_symmetries),
        ("ligand_symmetries", ligand_symmetries),
    ]:
        for group in sym_groups:
            for option in group:
                if not isinstance(option, (list, tuple)):
                    raise ValueError(
                        f"{sym_name} symmetry operation must be a list of (i, j) pairs, got {type(option)}"
                    )
                for swap_pair in option:
                    if not isinstance(swap_pair, (list, tuple)) or len(swap_pair) != 2:
                        raise ValueError(f"{sym_name} entries must be 2-tuples (i, j), got {swap_pair}")

    # --- Gather coords along CP axes, keep DP sharding ---
    cp_gathered_placements = tuple(
        coords_placements[0] if i == 0 else Replicate() for i in range(len(coords_placements))
    )
    coords_cp_gathered = coords.redistribute(coords_mesh, cp_gathered_placements)
    coords_local = coords_cp_gathered.to_local()  # [B_local * mul, N_atoms_padded, 3]

    coords_single_padded = coords_local[i_batch_multiplicity_local : i_batch_multiplicity_local + 1]

    atom_pad_mask_cp_gathered = atom_pad_mask.redistribute(coords_mesh, cp_gathered_placements)
    atom_pad_mask_local = atom_pad_mask_cp_gathered.to_local()
    mask_single = atom_pad_mask_local[index_batch_local].bool()  # [N_atoms_padded]

    # Remove interspersed padding: [1, N_atoms_padded, 3] -> [1, N_atoms_no_pad, 3]
    coords_single = coords_single_padded[:, mask_single, :]

    # Index symmetry features (plain tensors, no multiplicity)
    all_coords_indexed = all_coords[index_batch_local].unsqueeze(0).to(coords_single)
    all_resolved_mask_indexed = all_resolved_mask[index_batch_local].to(coords_single).to(torch.bool)
    crop_to_all_atom_map_indexed = crop_to_all_atom_map[index_batch_local].to(coords_single).to(torch.long)

    n_crop = int(crop_to_all_atom_map_indexed.numel())
    pred_coords_crop = coords_single[:, :n_crop]

    # --- Chain swap selection (Boltz-2 semantics) ---
    best_true_coords = all_coords_indexed[:, crop_to_all_atom_map_indexed].clone()
    best_true_resolved_mask = all_resolved_mask_indexed[crop_to_all_atom_map_indexed].clone()
    best_lddt = -1.0

    for c in chain_swaps:
        true_all_coords = all_coords_indexed.clone()
        true_all_resolved_mask = all_resolved_mask_indexed.clone()
        for start1, end1, start2, end2, _chainidx1, _chainidx2 in c:
            true_all_coords[:, start1:end1] = all_coords_indexed[:, start2:end2]
            true_all_resolved_mask[start1:end1] = all_resolved_mask_indexed[start2:end2]

        true_coords = true_all_coords[:, crop_to_all_atom_map_indexed]
        true_resolved_mask = true_all_resolved_mask[crop_to_all_atom_map_indexed]

        mask_row = true_resolved_mask.unsqueeze(0)
        mask_col = true_resolved_mask.unsqueeze(0)
        lddt = cdist_lddt(
            pred_coords_row=pred_coords_crop,
            pred_coords_col=pred_coords_crop,
            true_coords_row=true_coords,
            true_coords_col=true_coords,
            mask_row=mask_row,
            mask_col=mask_col,
            multiplicity=1,
            cutoff=15.0,
            per_atom=False,
        )[0].item()

        if lddt > best_lddt and torch.sum(true_resolved_mask) > 3:
            best_lddt = lddt
            best_true_coords = true_coords
            best_true_resolved_mask = true_resolved_mask

    # --- Atom-level symmetries (Boltz-2 semantics: best improvement) ---
    true_coords = best_true_coords.clone()
    true_resolved_mask = best_true_resolved_mask.clone()
    for symmetric_amino_or_lig in amino_acids_symmetries + ligand_symmetries:
        best_lddt_improvement = 0.0

        # Precompute all unique indices across all options in this group
        indices_set: set[int] = set()
        for c in symmetric_amino_or_lig:
            for i, j in c:
                indices_set.add(i)
        if len(indices_set) == 0:
            continue
        indices = sorted(indices_set)
        indices = torch.as_tensor(indices, device=pred_coords_crop.device, dtype=torch.long)
        pred_coords_subset = pred_coords_crop[:, indices]

        for c in symmetric_amino_or_lig:
            new_true_coords = true_coords.clone()
            new_true_resolved_mask = true_resolved_mask.clone()
            for i, j in c:
                new_true_coords[:, i] = true_coords[:, j]
                new_true_resolved_mask[i] = true_resolved_mask[j]

            true_coords_subset = true_coords[:, indices]
            new_true_coords_subset = new_true_coords[:, indices]

            mask_row = true_resolved_mask.unsqueeze(0)
            mask_col = true_resolved_mask[indices].unsqueeze(0)
            indices_batch = indices.unsqueeze(0)
            lddt = cdist_lddt(
                pred_coords_row=pred_coords_crop,
                pred_coords_col=pred_coords_subset,
                true_coords_row=true_coords,
                true_coords_col=true_coords_subset,
                mask_row=mask_row,
                mask_col=mask_col,
                multiplicity=1,
                atom_indices_col=indices_batch,
                cutoff=15.0,
                per_atom=False,
            )[0].item()

            new_mask_row = new_true_resolved_mask.unsqueeze(0)
            new_mask_col = new_true_resolved_mask[indices].unsqueeze(0)
            new_lddt = cdist_lddt(
                pred_coords_row=pred_coords_crop,
                pred_coords_col=pred_coords_subset,
                true_coords_row=new_true_coords,
                true_coords_col=new_true_coords_subset,
                mask_row=new_mask_row,
                mask_col=new_mask_col,
                multiplicity=1,
                atom_indices_col=indices_batch,
                cutoff=15.0,
                per_atom=False,
            )[0].item()

            lddt_improvement = new_lddt - lddt
            if lddt_improvement > best_lddt_improvement:
                best_true_coords = new_true_coords
                best_true_resolved_mask = new_true_resolved_mask
                best_lddt_improvement = lddt_improvement

        true_coords = best_true_coords.clone()
        true_resolved_mask = best_true_resolved_mask.clone()

    # --- Re-add interspersed padding and wrap as DTensors ---
    # Shape consistency check (same as Boltz-1 CP): boolean masking at line 167 removes
    # ALL padding (both trailing batch padding and interspersed CP padding) from
    # coords_single, so coords_single.shape[1] = sum(atom_pad_mask) = n_crop.
    # true_coords.shape[1] = len(crop_to_all_atom_map) = n_crop. They must be equal.
    n_atoms_padded = coords_single_padded.shape[1]
    n_real = int(mask_single.sum())
    if true_coords.shape[1] != n_real:
        raise ValueError(
            f"Shape mismatch: true_coords.shape[1]={true_coords.shape[1]} != "
            f"sum(atom_pad_mask)={n_real}. Both should equal the number of crop atoms."
        )
    if true_resolved_mask.shape[0] != n_real:
        raise ValueError(
            f"Shape mismatch: true_resolved_mask.shape[0]={true_resolved_mask.shape[0]} != "
            f"sum(atom_pad_mask)={n_real}. Both should equal the number of crop atoms."
        )

    true_coords_padded = torch.zeros((1, n_atoms_padded, 3), dtype=true_coords.dtype, device=true_coords.device)
    true_coords_padded[:, mask_single, :] = true_coords

    true_resolved_mask_padded = torch.zeros(
        (n_atoms_padded,), dtype=true_resolved_mask.dtype, device=true_resolved_mask.device
    )
    true_resolved_mask_padded[mask_single] = true_resolved_mask

    device_mesh = coords_mesh
    placements = coords_placements

    # Broadcast from CP rank 0 to ensure bitwise-identical results across CP ranks
    # Distribute with the target CP placements directly
    # so each rank only receives the shard it needs.
    cp_mesh = device_mesh["cp_axis_0", "cp_axis_1"]
    cp_shard_placements = (placements[1], placements[2])  # placements: (dp, cp_axis_0, cp_axis_1)

    true_coords_cp = distribute_tensor(
        true_coords_padded, device_mesh=cp_mesh, placements=cp_shard_placements, src_data_rank=0
    )
    _coords_global_shape = true_coords_padded.shape
    true_coords_dtensor = DTensor.from_local(
        true_coords_cp.to_local(),
        device_mesh,
        placements=placements,
        shape=_coords_global_shape,
        stride=true_coords_padded.stride(),
    )

    true_mask_unsqueezed = true_resolved_mask_padded.unsqueeze(0)
    true_mask_cp = distribute_tensor(
        true_mask_unsqueezed, device_mesh=cp_mesh, placements=cp_shard_placements, src_data_rank=0
    )
    _mask_global_shape = true_mask_unsqueezed.shape
    true_resolved_mask_dtensor = DTensor.from_local(
        true_mask_cp.to_local(),
        device_mesh,
        placements=placements,
        shape=_mask_global_shape,
        stride=true_mask_unsqueezed.stride(),
    )

    return true_coords_dtensor, true_resolved_mask_dtensor
