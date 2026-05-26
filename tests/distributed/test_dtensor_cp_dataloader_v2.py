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

# System imports
from pathlib import Path
from typing import Dict, Optional

# Third party imports
import pytest
import torch
from torch.distributed.tensor import DTensor
from torch.utils.data import DataLoader, DistributedSampler

from boltz.data.module.inferencev2 import Boltz2InferenceDataModule, PredictionDataset, collate
from boltz.data.module.trainingv2 import (
    Boltz2TrainingDataModule,
)
from boltz.data.module.trainingv2 import (
    collate as collate_training,
)
from boltz.data.types import Manifest
from boltz.distributed.data.module.inferencev2 import (
    Boltz2InferenceDataModuleDTensor as BoltzInferenceDataModuleDTensor,
)
from boltz.distributed.data.module.trainingv2 import (
    Boltz2TrainingDataModule as BoltzTrainingDataModuleDTensor,
)
from boltz.distributed.data.types import PairMaskMode
from boltz.distributed.data.utils import (
    ATOM_FEATURES_V2 as ATOM_FEATURES,
)
from boltz.distributed.data.utils import (
    NON_SHARDED_FEATURES_V2,
    map_subgroup_mesh_to_cpu,
)
from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.layers.atom_to_token import (
    single_repr_atom_to_token,
    single_repr_token_to_atom,
)
from boltz.distributed.testing.utils import setup_mock_training_datamodule_config
from boltz.main import BoltzProcessedInput
from boltz.testing.utils import (
    concat_data,
    seed_by_rank,
    spawn_multiprocessing,
)

# TODO support CP constraint, template, affinity, ensemble features for inference
INFERENCE_FEATURES_DIFFERENCE = {
    "atom_counts_per_token",
    "chiral_atom_index",
    "chiral_atom_orientations",
    "chiral_reference_mask",
    "connected_atom_index",
    "connected_chain_index",
    "contact_negation_mask",
    "contact_pair_index",
    "contact_thresholds",
    "contact_union_index",
    "ensemble_ref_idxs",
    "pair_mask",
    "planar_bond_index",
    "planar_ring_5_index",
    "planar_ring_6_index",
    "query_to_template",
    "rdkit_bounds_angle_mask",
    "rdkit_bounds_bond_mask",
    "rdkit_bounds_index",
    "rdkit_lower_bounds",
    "rdkit_upper_bounds",
    "record",
    "r_set_to_rep_atom",
    "stereo_bond_index",
    "stereo_bond_orientations",
    "stereo_reference_mask",
    "symmetric_chain_index",
    "template_ca",
    "template_cb",
    "template_frame_rot",
    "template_frame_t",
    "template_mask",
    "template_mask_cb",
    "template_mask_frame",
    "template_restype",
    "token_to_center_atom",
    "token_pair_pad_mask",
    "visibility_ids",
}

TRAINING_FEATURES_DIFFERENCE = {
    "atom_counts_per_token",
    "ensemble_ref_idxs",
    "idx_dataset",
    "record",
    "token_pair_pad_mask",
}

TOKEN_PAIR_FEATURES = {
    "contact_conditioning",
    "contact_threshold",
    "disto_target",
    "token_bonds",
    "token_pair_pad_mask",
    "type_bonds",
}

MSA_FEATURES = {
    "deletion_value",
    "has_deletion",
    "msa",
    "msa_paired",
}

ENSEMBLE_ATOM_FEATURES = {
    "coords",
}

ENSEMBLE_TOKEN_FEATURES = {
    "disto_coords_ensemble",
}


def _map_padded_frames_idx_to_unpadded(
    frames_idx_dtensor: DTensor,
    dp_idx_str: int,
    dp_idx_end: int,
    atom_to_token_dtensor: DTensor,
    frame_resolved_mask_dtensor: DTensor,
) -> torch.Tensor:
    """Map frames_idx from padded atom indices to unpadded indices using atom_to_token."""
    frames_idx_padded = frames_idx_dtensor.full_tensor()[dp_idx_str:dp_idx_end]
    if frames_idx_padded.ndim != 4:
        raise ValueError(f"frames_idx_padded must have ndim=4 (B, E, T, 3), got {frames_idx_padded.ndim}")
    frame_resolved_mask_padded = frame_resolved_mask_dtensor.full_tensor()[dp_idx_str:dp_idx_end].bool()
    if frame_resolved_mask_padded.ndim != 3:
        raise ValueError(
            f"frame_resolved_mask_padded must have ndim=3 (B, E, T), got {frame_resolved_mask_padded.ndim}"
        )

    atom_to_token_local = atom_to_token_dtensor.to_local()
    if atom_to_token_local.ndim != 3:
        raise ValueError(
            f"atom_to_token must have ndim=3 (B, N_atoms, N_tokens_per_row), got {atom_to_token_local.ndim}"
        )
    if atom_to_token_local.shape[0] != frames_idx_padded.shape[0]:
        raise ValueError(
            "atom_to_token batch size does not match frames_idx_padded: "
            f"{atom_to_token_local.shape[0]} vs {frames_idx_padded.shape[0]}"
        )
    device_mesh = atom_to_token_dtensor.device_mesh
    cp_axis_0_group = device_mesh.get_group(1)
    # Local atom mask for this CP row
    atom_mask_local = atom_to_token_local.sum(dim=2) > 0  # (B, N_atoms_per_shard)
    local_counts = atom_mask_local.sum(dim=1).to(torch.int64)  # (B,)

    # Gather counts and masks across CP axis 0
    counts_list = [torch.empty_like(local_counts) for _ in range(cp_axis_0_group.size())]
    torch.distributed.all_gather(counts_list, local_counts, group=cp_axis_0_group)
    counts_all = torch.stack(counts_list, dim=0)  # (cp_rows, B)

    mask_list = [torch.empty_like(atom_mask_local) for _ in range(cp_axis_0_group.size())]
    torch.distributed.all_gather(mask_list, atom_mask_local, group=cp_axis_0_group)
    atom_mask_all = torch.stack(mask_list, dim=0)  # (cp_rows, B, N_atoms_per_shard)

    # Row offsets for global unpadded indices
    row_offsets = torch.cumsum(counts_all, dim=0) - counts_all  # (cp_rows, B)

    # CollateDTensor remaps frames_idx to the final padded stride, so the
    # padded row width is always the local atom dim after collation.
    padded_row_width = atom_to_token_local.shape[1]

    batch_size = frames_idx_padded.shape[0]
    mapped = torch.empty_like(frames_idx_padded)
    for batch_idx in range(batch_size):
        frames = frames_idx_padded[batch_idx]
        resolved_frames = frame_resolved_mask_padded[batch_idx]
        resolved_frames = resolved_frames.unsqueeze(-1).expand_as(frames)

        row_idx = torch.div(frames, padded_row_width, rounding_mode="floor")
        local_idx = frames % padded_row_width

        if row_idx.max().item() >= atom_mask_all.shape[0]:
            raise ValueError(
                f"frames_idx_padded has row index {row_idx.max().item()} out of range for cp rows "
                f"{atom_mask_all.shape[0]}"
            )
        # Build prefix sums per row for mapping local -> unpadded
        prefix = torch.cumsum(atom_mask_all[:, batch_idx].to(torch.int64), dim=1) - 1
        local_mask = atom_mask_all[:, batch_idx]
        # Validate no padding indices are referenced for resolved frames.
        valid = local_mask[row_idx, local_idx]
        if not torch.all(valid[resolved_frames]):
            raise ValueError("frames_idx points to padded atom indices in atom_to_token")
        unpadded_local = prefix[row_idx, local_idx]
        mapped_vals = row_offsets[row_idx, batch_idx] + unpadded_local
        mapped[batch_idx] = torch.where(resolved_frames, mapped_vals, torch.zeros_like(mapped_vals))
    return mapped


def _compare_frames_idx_valid_tokens(
    frames_idx_non_dtensor: torch.Tensor,
    frames_idx_dtensor_unpadded: torch.Tensor,
    frame_resolved_mask_non_dtensor: torch.Tensor,
    frame_resolved_mask_dtensor: torch.Tensor,
    non_dtensor_mask: torch.Tensor,
    dtensor_mask: torch.Tensor,
) -> None:
    """Compare frames_idx using token and frame-resolved masks on each side."""
    if non_dtensor_mask.ndim != 2 or dtensor_mask.ndim != 2:
        raise ValueError(f"token_pad_mask must have ndim=2, got {non_dtensor_mask.ndim} and {dtensor_mask.ndim}")
    batch_size = frames_idx_dtensor_unpadded.shape[0]
    for batch_idx in range(batch_size):
        non_mask = non_dtensor_mask[batch_idx].bool()
        dt_mask = dtensor_mask[batch_idx].bool()
        if non_mask.sum().item() != dt_mask.sum().item():
            raise AssertionError(
                f"frames_idx token count mismatch: non-dtensor={non_mask.sum().item()}, "
                f"dtensor={dt_mask.sum().item()}"
            )
        non_item = frames_idx_non_dtensor[batch_idx]
        dt_item = frames_idx_dtensor_unpadded[batch_idx]
        non_frame_mask = frame_resolved_mask_non_dtensor[batch_idx].bool()
        dt_frame_mask = frame_resolved_mask_dtensor[batch_idx].bool()
        if non_item.ndim != 3 or dt_item.ndim != 3:
            raise ValueError(
                "frames_idx must be ensemble-aware with ndim=3 per batch item; got "
                f"{non_item.ndim} and {dt_item.ndim}"
            )
        non_mask_ens = non_mask.unsqueeze(0).expand(non_item.shape[0], -1) & non_frame_mask
        dt_mask_ens = dt_mask.unsqueeze(0).expand(dt_item.shape[0], -1) & dt_frame_mask
        non_tokens = non_item[non_mask_ens]
        dt_tokens = dt_item[dt_mask_ens]
        torch.testing.assert_close(non_tokens, dt_tokens, atol=0, rtol=0)


def compare_r_set_to_rep_atom(
    ref_tensor: torch.Tensor,
    dtensor_full: torch.Tensor,
    atom_pad_mask_dtensor: torch.Tensor,
    n_shards: int,
) -> None:
    """Compare r_set_to_rep_atom: reference (one-hot) vs DTensor (one-hot diagonal blocks).

    Simplified comparison using:
    - .any(dim=-1) to find valid R-set rows
    - .argmax(dim=-1) to get LOCAL atom indices
    - atom_pad_mask reshaped per-shard to compute atom offsets
    - Convert local -> global indices using per-shard offsets

    Parameters
    ----------
    ref_tensor : torch.Tensor
        Single-device reference tensor (one-hot), shape [B, size_r_set_ref, N_atom_global]
    dtensor_full : torch.Tensor
        DTensor full_tensor() (one-hot diagonal blocks), shape [B, size_r_set_padded, max_atoms_per_shard]
    atom_pad_mask_dtensor : torch.Tensor
        DTensor atom_pad_mask.full_tensor(), shape [B, N_atoms_padded]. True/1 = valid atom.
    n_shards : int
        Number of shards along the R-set sharding axis (cp_axis_0 size).
    """
    # Handle both 2D and 3D tensors
    if ref_tensor.dim() == 2:
        ref_tensor = ref_tensor.unsqueeze(0)
    if dtensor_full.dim() == 2:
        dtensor_full = dtensor_full.unsqueeze(0)
    if atom_pad_mask_dtensor.dim() == 1:
        atom_pad_mask_dtensor = atom_pad_mask_dtensor.unsqueeze(0)

    assert ref_tensor.dim() == 3, f"Expected 3D ref tensor [B, size_r_set, N_atom], got {ref_tensor.shape}"
    assert (
        dtensor_full.dim() == 3
    ), f"Expected 3D dtensor [B, size_r_set_padded, max_atoms_per_shard], got {dtensor_full.shape}"
    assert (
        atom_pad_mask_dtensor.dim() == 2
    ), f"Expected 2D atom_pad_mask [B, N_atoms_padded], got {atom_pad_mask_dtensor.shape}"

    batch_size = ref_tensor.shape[0]
    size_r_set_padded = dtensor_full.shape[1]
    max_atoms_per_shard = dtensor_full.shape[2]
    max_size_r_set_per_shard = size_r_set_padded // n_shards

    for b in range(batch_size):
        ref_b = ref_tensor[b]  # [size_r_set_ref, N_atom_global] one-hot
        dtensor_b = dtensor_full[b]  # [size_r_set_padded, max_atoms_per_shard] one-hot
        atom_mask_b = atom_pad_mask_dtensor[b].long()  # [N_atoms_padded]

        # Find valid R-set rows (non-zero rows)
        ref_valid_mask = ref_b.any(dim=-1)  # [size_r_set_ref]
        dtensor_valid_mask = dtensor_b.any(dim=-1)  # [size_r_set_padded]

        size_r_set_valid_ref = ref_valid_mask.sum().item()
        size_r_set_valid_dtensor = dtensor_valid_mask.sum().item()

        if size_r_set_valid_ref == 0:
            # No valid R-set elements - DTensor should have no valid rows either
            assert (
                size_r_set_valid_dtensor == 0
            ), f"DTensor batch {b} has {size_r_set_valid_dtensor} valid rows, expected 0"
            continue

        # Verify counts match
        if size_r_set_valid_ref != size_r_set_valid_dtensor:
            raise AssertionError(
                f"r_set_to_rep_atom batch {b}: count mismatch. "
                f"Ref has {size_r_set_valid_ref} valid rows, DTensor has {size_r_set_valid_dtensor}"
            )

        # Get global atom indices from reference (straightforward argmax)
        ref_global_indices = ref_b[ref_valid_mask].argmax(dim=-1)  # [size_r_set_valid]

        # Get LOCAL atom indices from DTensor (argmax gives local index within shard)
        dtensor_atom_ids_local = dtensor_b.argmax(dim=-1)  # [size_r_set_padded] local indices

        # Compute per-shard atom offset using atom_pad_mask
        # atom_pad_mask_per_shard[s] = mask for shard s's atoms
        atom_pad_mask_per_shard = atom_mask_b.reshape(n_shards, max_atoms_per_shard)  # [n_shards, max_atoms_per_shard]
        atom_counts_per_shard = atom_pad_mask_per_shard.sum(dim=-1)  # [n_shards]
        # atom_offset[s] = cumsum of atoms before shard s
        atom_offset_per_shard = atom_counts_per_shard.cumsum(dim=0) - atom_counts_per_shard  # [n_shards]

        # Reshape local atom ids by shard: [n_shards, max_size_r_set_per_shard]
        dtensor_atom_ids_local_per_shard = dtensor_atom_ids_local.reshape(n_shards, max_size_r_set_per_shard)

        # Add per-shard offset to convert local -> global
        # atom_offset_per_shard: [n_shards] -> [n_shards, 1] for broadcasting
        dtensor_atom_ids_global_per_shard = dtensor_atom_ids_local_per_shard + atom_offset_per_shard.unsqueeze(-1)

        # Flatten and select valid rows
        dtensor_global_indices = dtensor_atom_ids_global_per_shard.flatten()[
            dtensor_valid_mask
        ]  # [size_r_set_valid_dtensor]

        # Sort both and compare (order may differ due to token-aligned sharding)
        ref_sorted = ref_global_indices.sort().values.to(torch.int64)
        dtensor_sorted = dtensor_global_indices.sort().values.to(torch.int64)

        if not torch.equal(ref_sorted, dtensor_sorted):
            # Find first mismatch for debugging
            mismatch_mask = ref_sorted != dtensor_sorted
            first_mismatch = mismatch_mask.nonzero(as_tuple=True)[0][0].item() if mismatch_mask.any() else -1
            raise AssertionError(
                f"r_set_to_rep_atom batch {b}: global atom indices mismatch.\n"
                f"First mismatch at sorted index {first_mismatch}\n"
                f"Ref (sorted)[:10]: {ref_sorted[:10].tolist()}\n"
                f"DTensor (sorted)[:10]: {dtensor_sorted[:10].tolist()}"
            )


def compare_token_to_rep_atom(
    ref_tensor: torch.Tensor,
    dtensor_full: torch.Tensor,
    atom_pad_mask_ref: torch.Tensor,
    token_pad_mask_ref: torch.Tensor,
    token_pad_mask_dtensor: torch.Tensor,
    n_shards: int,
    atom_counts_per_token: torch.Tensor,
) -> None:
    """Compare serial token_to_rep_atom against DTensor block-diagonal version.

    In the serial pipeline, token_to_rep_atom is (N_tokens, N_atoms) where entry [t, a]
    is nonzero iff atom a is the representative atom for token t. In the distributed
    pipeline, it is stored as diagonal blocks: shard i holds the slice for its token and
    atom ranges, padded to uniform per-shard dimensions.

    V2 adds two padding layers that affect this comparison:
    1) featurizer padding (per-sample token count rounded up to CP divisibility)
    2) collation padding (batched per-shard token count rounded up across samples/ranks)

    This helper handles those layers directly, so callers can pass unmodified serial
    tensors. In particular, serial row boundaries are clipped against cumsum bounds and
    atom-empty shards are skipped.

    Parameters
    ----------
    ref_tensor : torch.Tensor
        Serial token_to_rep_atom. Shape: (B, N_tokens_serial, N_atoms_serial).
    dtensor_full : torch.Tensor
        DTensor token_to_rep_atom after full_tensor().
        Shape: (B, n_shards * N_tps_collated, max_atoms_per_shard_collated).
    atom_pad_mask_ref : torch.Tensor
        Serial atom padding mask. Shape: (B, N_atoms_serial).
    token_pad_mask_ref : torch.Tensor
        Serial token padding mask. Shape: (B, N_tokens_serial).
    token_pad_mask_dtensor : torch.Tensor
        DTensor token padding mask (after full_tensor).
        Shape: (B, n_shards * N_tps_collated).
    n_shards : int
        Number of CP shards along mesh dim 0 (the token-sharding axis).
    atom_counts_per_token : torch.Tensor
        Serial atom counts per token. Shape: (B, N_tokens_serial).

    """
    batch_size = ref_tensor.shape[0]
    for b in range(batch_size):
        ref_b = ref_tensor[b]
        dtensor_b = dtensor_full[b]
        token_mask_ref = token_pad_mask_ref[b].bool()
        token_mask_dt = token_pad_mask_dtensor[b].bool()
        atom_mask_ref = atom_pad_mask_ref[b].bool()

        N_tps_collated = dtensor_b.shape[0] // n_shards

        N_tps_featurizer = int(token_mask_dt[0:N_tps_collated].sum().item())

        cumsum = torch.cat(
            [torch.tensor([0], device=atom_counts_per_token.device), atom_counts_per_token[b].cumsum(dim=0)]
        )

        for shard_i in range(n_shards):
            ref_row_start = shard_i * N_tps_featurizer
            ref_row_end = (shard_i + 1) * N_tps_featurizer
            ref_row_start_clipped = min(ref_row_start, cumsum.shape[0] - 1)
            ref_row_end_clipped = min(ref_row_end, cumsum.shape[0] - 1)
            ref_atom_start = int(cumsum[ref_row_start_clipped].item())
            ref_atom_end = int(cumsum[ref_row_end_clipped].item())
            n_atoms_in_shard = ref_atom_end - ref_atom_start

            if n_atoms_in_shard == 0:
                continue

            ref_block = ref_b[ref_row_start:ref_row_end, ref_atom_start:ref_atom_end]

            ref_tok_mask = token_mask_ref[ref_row_start:ref_row_end]
            ref_atom_mask = atom_mask_ref[ref_atom_start:ref_atom_end]
            ref_unpadded = ref_block[ref_tok_mask][:, ref_atom_mask]

            dt_row_start = shard_i * N_tps_collated
            dt_row_end = (shard_i + 1) * N_tps_collated
            dt_block = dtensor_b[dt_row_start:dt_row_end, :n_atoms_in_shard]

            dt_tok_mask = token_mask_dt[dt_row_start:dt_row_end]
            dt_unpadded = dt_block[dt_tok_mask][:, ref_atom_mask]

            torch.testing.assert_close(
                dt_unpadded,
                ref_unpadded,
                atol=0,
                rtol=0,
                msg=f"token_to_rep_atom mismatch at batch {b}, shard {shard_i}",
            )


def _compare_features(
    common_keys,
    skip_keys,
    data_batch_serial,
    data_batch_dtensor,
    dp_idx_str,
    dp_idx_end,
    n_shards,
    rank,
    manager,
):
    """Compare serial and DTensor features across all common keys.

    Extracts masks, applies per-feature-type masking, and compares values.
    Collects all errors and reports them together. Guards against vacuous
    comparisons by requiring at least one feature with non-trivial data.
    """
    atom_pad_mask_dtensor_full = data_batch_dtensor["atom_pad_mask"].full_tensor()[dp_idx_str:dp_idx_end].bool()
    token_pad_mask_dtensor_full = data_batch_dtensor["token_pad_mask"].full_tensor()[dp_idx_str:dp_idx_end].bool()
    msa_mask_dtensor_full = data_batch_dtensor["msa_mask"].full_tensor()[dp_idx_str:dp_idx_end].bool()
    atom_pad_mask_serial = data_batch_serial["atom_pad_mask"].bool()
    token_pad_mask_serial = data_batch_serial["token_pad_mask"].bool()
    msa_mask_serial = data_batch_serial["msa_mask"].bool()

    token_pad_pair_mask_dtensor_full = token_pad_mask_dtensor_full[:, :, None] * token_pad_mask_dtensor_full[:, None, :]
    token_pad_pair_mask_serial = token_pad_mask_serial[:, :, None] * token_pad_mask_serial[:, None, :]

    errors = []
    any_nonempty_data = False

    for key in common_keys:
        if key in skip_keys:
            continue

        feature_serial = data_batch_serial[key]
        feature_dtensor = data_batch_dtensor[key]
        if not isinstance(feature_serial, torch.Tensor) or not isinstance(feature_dtensor, DTensor):
            continue

        if key == "r_set_to_rep_atom":
            compare_r_set_to_rep_atom(
                ref_tensor=feature_serial,
                dtensor_full=feature_dtensor.full_tensor()[dp_idx_str:dp_idx_end],
                atom_pad_mask_dtensor=atom_pad_mask_dtensor_full,
                n_shards=n_shards,
            )
            any_nonempty_data = True
            continue
        elif key == "token_to_rep_atom":
            compare_token_to_rep_atom(
                ref_tensor=feature_serial,
                dtensor_full=feature_dtensor.full_tensor()[dp_idx_str:dp_idx_end],
                atom_pad_mask_ref=atom_pad_mask_serial,
                token_pad_mask_ref=token_pad_mask_serial,
                token_pad_mask_dtensor=token_pad_mask_dtensor_full,
                n_shards=n_shards,
                atom_counts_per_token=data_batch_serial["atom_counts_per_token"],
            )
            any_nonempty_data = True
            continue

        feature_dtensor_full = None
        feature_serial_full = None
        try:
            if key in ENSEMBLE_ATOM_FEATURES:
                # (B, E, A, ...) — expand atom mask with ensemble dim
                dt_full = feature_dtensor.full_tensor()[dp_idx_str:dp_idx_end]
                mask_dt = atom_pad_mask_dtensor_full.unsqueeze(1).expand(-1, dt_full.shape[1], -1)
                mask_serial = atom_pad_mask_serial.unsqueeze(1).expand(-1, feature_serial.shape[1], -1)
                feature_dtensor_full = dt_full[mask_dt]
                feature_serial_full = feature_serial[mask_serial]
            elif key in ENSEMBLE_TOKEN_FEATURES:
                # (B, E, T, ...) — expand token mask with ensemble dim
                dt_full = feature_dtensor.full_tensor()[dp_idx_str:dp_idx_end]
                mask_dt = token_pad_mask_dtensor_full.unsqueeze(1).expand(-1, dt_full.shape[1], -1)
                mask_serial = token_pad_mask_serial.unsqueeze(1).expand(-1, feature_serial.shape[1], -1)
                feature_dtensor_full = dt_full[mask_dt]
                feature_serial_full = feature_serial[mask_serial]
            elif key in TOKEN_PAIR_FEATURES:
                feature_dtensor_full = feature_dtensor.full_tensor()[dp_idx_str:dp_idx_end][
                    token_pad_pair_mask_dtensor_full
                ]
                feature_serial_full = feature_serial[token_pad_pair_mask_serial]
            elif key in MSA_FEATURES:
                feature_dtensor_full = feature_dtensor.full_tensor()[dp_idx_str:dp_idx_end][msa_mask_dtensor_full]
                feature_serial_full = feature_serial[msa_mask_serial]
            elif key == "frame_resolved_mask":
                dt_full = feature_dtensor.full_tensor()[dp_idx_str:dp_idx_end]
                if dt_full.ndim != 3:
                    raise ValueError(
                        f"frame_resolved_mask must be ensemble-aware with ndim=3 (B, E, T), got {dt_full.ndim}"
                    )
                mask_dt = token_pad_mask_dtensor_full.unsqueeze(1).expand(-1, dt_full.shape[1], -1)
                mask_serial = token_pad_mask_serial.unsqueeze(1).expand(-1, feature_serial.shape[1], -1)
                feature_dtensor_full = dt_full[mask_dt]
                feature_serial_full = feature_serial[mask_serial]
            elif key == "frames_idx":
                frame_resolved_mask_dtensor_full = data_batch_dtensor["frame_resolved_mask"].full_tensor()[
                    dp_idx_str:dp_idx_end
                ]
                frames_idx_dtensor_unpadded = _map_padded_frames_idx_to_unpadded(
                    feature_dtensor,
                    dp_idx_str,
                    dp_idx_end,
                    data_batch_dtensor["atom_to_token"],
                    data_batch_dtensor["frame_resolved_mask"],
                )
                _compare_frames_idx_valid_tokens(
                    feature_serial,
                    frames_idx_dtensor_unpadded,
                    data_batch_serial["frame_resolved_mask"],
                    frame_resolved_mask_dtensor_full,
                    token_pad_mask_serial,
                    token_pad_mask_dtensor_full,
                )
                continue
            elif key in ATOM_FEATURES:
                feature_dtensor_full = feature_dtensor.full_tensor()[dp_idx_str:dp_idx_end][atom_pad_mask_dtensor_full]
                feature_serial_full = feature_serial[atom_pad_mask_serial]
            else:
                feature_dtensor_full = feature_dtensor.full_tensor()[dp_idx_str:dp_idx_end][token_pad_mask_dtensor_full]
                feature_serial_full = feature_serial[token_pad_mask_serial]

            if feature_dtensor_full.numel() > 0:
                if feature_dtensor_full.is_floating_point():
                    if feature_dtensor_full.any():
                        any_nonempty_data = True
                else:
                    if feature_dtensor_full.unique().numel() > 1:
                        any_nonempty_data = True

            torch.testing.assert_close(
                feature_dtensor_full,
                feature_serial_full,
                atol=0,
                rtol=0,
            )
        except Exception as e:
            dt_shape = feature_dtensor_full.shape if feature_dtensor_full is not None else "N/A"
            serial_shape = feature_serial_full.shape if feature_serial_full is not None else "N/A"
            errors.append(
                f"[{key}] rank {rank} cp_rank {manager.group_rank['cp']}: "
                f"dtensor shape {dt_shape} vs serial shape {serial_shape} "
                f"| {e}"
            )

    if errors:
        raise AssertionError(f"Shape/value mismatches on rank {rank} ({len(errors)} failures):\n" + "\n".join(errors))

    assert any_nonempty_data, (
        f"All compared features were empty or zero on rank {rank} -- "
        f"likely a padding/masking bug producing vacuous comparisons"
    )


def parallel_assert_cp_inference_dataloader(
    rank: int,
    processed,
    canonical_mols_dir: Path,
    local_batch_size: int,
    device_type,
    backend,
    grid_group_sizes: Dict[str, int],
    env_map: Optional[dict[str, str]] = None,
):
    monkeypatch = pytest.MonkeyPatch()
    if env_map is not None:
        for var_name, value in env_map.items():
            if value == "<INPUT_RANK>":
                monkeypatch.setenv(var_name, f"{rank}")
                continue
            monkeypatch.setenv(var_name, value)

    # TODO VI: hardcode device_type and backend for now. Do we need GPU?
    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    manager = DistributedManager()
    # create a CPU process groups for comm python objects
    DistributedManager.create_group(
        "world_cpu", manager.group_ranks["world"], backend="gloo", use_local_synchronization=True
    )
    DistributedManager.create_group("cp_cpu", manager.group_ranks["cp"], backend="gloo", use_local_synchronization=True)

    device_mesh = manager.device_mesh_subgroups
    cp_device_mesh = map_subgroup_mesh_to_cpu(manager)
    assert device_mesh.shape == cp_device_mesh.shape
    dp_rank = device_mesh.get_local_rank(0)

    # Serial dataloader
    seed_by_rank(0, seed=42)
    data_module_serial = Boltz2InferenceDataModule(
        manifest=processed.manifest,
        target_dir=processed.targets_dir,
        msa_dir=processed.msa_dir,
        mol_dir=canonical_mols_dir,
        num_workers=0,
    )
    dataset = PredictionDataset(
        manifest=processed.manifest,
        target_dir=processed.targets_dir,
        msa_dir=processed.msa_dir,
        mol_dir=canonical_mols_dir,
    )
    sampler = DistributedSampler(
        dataset,
        num_replicas=device_mesh.shape[0],
        rank=device_mesh.get_local_rank(0),
        shuffle=False,
        drop_last=False,
    )
    dataloader_serial = DataLoader(
        dataset,
        sampler=sampler,
        batch_size=local_batch_size,
        num_workers=0,
        shuffle=False,
        collate_fn=collate,
    )
    data_batch_list_serial = list(dataloader_serial)

    # DTensor CP dataloader
    seed_by_rank(0, seed=42)
    data_module_dtensor = BoltzInferenceDataModuleDTensor(
        manifest=processed.manifest,
        target_dir=processed.targets_dir,
        msa_dir=processed.msa_dir,
        mol_dir=canonical_mols_dir,
        num_workers=0,
        device_mesh=device_mesh,
        device_mesh_cpu=cp_device_mesh,
        local_batch_size=local_batch_size,
        pair_mask_mode=PairMaskMode.SEQUENCE_LOCAL_ATTENTION,
    )
    data_module_dtensor.setup("predict")
    dataloader_dtensor = data_module_dtensor.predict_dataloader()
    data_batch_list_dtensor = list(dataloader_dtensor)

    skip_keys = {
        "atom_pad_mask",
        "token_pad_mask",
        "msa_mask",
        "token_pair_mask",
        "atom_to_token",
    } | NON_SHARDED_FEATURES_V2

    n_shards = device_mesh.get_group("cp_axis_0").size()

    for data_batch_serial, data_batch_dtensor in zip(data_batch_list_serial, data_batch_list_dtensor):
        serial_keys = set(data_batch_serial.keys())
        dtensor_keys = set(data_batch_dtensor.keys())
        assert (serial_keys - INFERENCE_FEATURES_DIFFERENCE) == (dtensor_keys - INFERENCE_FEATURES_DIFFERENCE), (
            f"Feature keys are different: {serial_keys - INFERENCE_FEATURES_DIFFERENCE} "
            f"!= {dtensor_keys - INFERENCE_FEATURES_DIFFERENCE}"
        )
        common_keys = sorted(serial_keys & dtensor_keys)

        data_batch_serial = data_module_serial.transfer_batch_to_device(data_batch_serial, manager.device, 0)
        data_batch_dtensor = data_module_dtensor.transfer_batch_to_device(data_batch_dtensor, manager.device, 0)

        dp_idx_str = dp_rank * local_batch_size
        dp_idx_end = dp_idx_str + local_batch_size

        _compare_features(
            common_keys=common_keys,
            skip_keys=skip_keys,
            data_batch_serial=data_batch_serial,
            data_batch_dtensor=data_batch_dtensor,
            dp_idx_str=dp_idx_str,
            dp_idx_end=dp_idx_end,
            n_shards=n_shards,
            rank=rank,
            manager=manager,
        )

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env",
    [
        ((1, (2, 2)), True, "cuda", "ENV"),
        ((2, (2, 2)), True, "cuda", "ENV"),
        ((1, (3, 3)), True, "cpu", "ENV"),
    ],
    indirect=("setup_env",),
    ids=lambda x: f"dp={x[0][0]}, cp={x[0][1]}, device_type={x[2]}",
)
def test_dtensor_cp_inference_dataloader(
    setup_env: tuple[dict, int, str, str, str, dict[str, str]],
    test_cp_training_base_data_dir_boltz2: Path,
    canonical_mols_dir: Path,
    tmp_path: Path,
    local_batch_size: int = 2,
):
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    # Create combined dataset with multiple samples to support local_batch_size > 1
    names = ["7ylz", "7z64", "8ayv", "8b2e"]
    combined_dataset_path = concat_data(
        tmp_path / "processed_combined",
        *[test_cp_training_base_data_dir_boltz2 / f"processed_{name}" for name in names],
    )

    # Check if we have enough samples for the requested batch size
    total_samples = len(names)
    if total_samples < local_batch_size * grid_group_sizes["dp"]:
        pytest.skip(
            f"Not enough samples ({total_samples}) for local_batch_size * dp_size ({local_batch_size * grid_group_sizes['dp']})"
        )

    # Create a new BoltzProcessedInput with the combined dataset
    combined_processed_handle = BoltzProcessedInput(
        manifest=Manifest.load(combined_dataset_path / "manifest.json"),
        targets_dir=combined_dataset_path / "structures",
        msa_dir=combined_dataset_path / "msa",
        template_dir=None,
        extra_mols_dir=None,
    )

    spawn_multiprocessing(
        parallel_assert_cp_inference_dataloader,
        world_size,
        combined_processed_handle,
        canonical_mols_dir,
        local_batch_size,
        device_type,
        backend,
        grid_group_sizes,
        env_per_rank,
    )


def parallel_assert_atom_to_token_sharding(
    rank: int,
    processed,
    canonical_mols_dir: Path,
    local_batch_size: int,
    device_type,
    backend,
    grid_group_sizes: Dict[str, int],
    env_map: Optional[dict[str, str]] = None,
):
    """Test that sharded atom_to_token operations produce the same results as non-sharded after removing padding."""
    monkeypatch = pytest.MonkeyPatch()
    if env_map is not None:
        for var_name, value in env_map.items():
            if value == "<INPUT_RANK>":
                monkeypatch.setenv(var_name, f"{rank}")
                continue
            monkeypatch.setenv(var_name, value)

    # Initialize distributed environment
    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    manager = DistributedManager()
    DistributedManager.create_group(
        "world_cpu", manager.group_ranks["world"], backend="gloo", use_local_synchronization=True
    )
    DistributedManager.create_group("cp_cpu", manager.group_ranks["cp"], backend="gloo", use_local_synchronization=True)

    device_mesh = manager.device_mesh_subgroups
    cp_device_mesh = map_subgroup_mesh_to_cpu(manager)
    dp_rank = device_mesh.get_local_rank(0)
    cp_rank = manager.group_rank["cp"]

    # Set up data modules
    seed_by_rank(0, seed=42)
    data_module_serial = Boltz2InferenceDataModule(
        manifest=processed.manifest,
        target_dir=processed.targets_dir,
        msa_dir=processed.msa_dir,
        mol_dir=canonical_mols_dir,
        num_workers=0,
    )
    dataset = PredictionDataset(
        manifest=processed.manifest,
        target_dir=processed.targets_dir,
        msa_dir=processed.msa_dir,
        mol_dir=canonical_mols_dir,
    )

    data_module_dtensor = BoltzInferenceDataModuleDTensor(
        manifest=processed.manifest,
        target_dir=processed.targets_dir,
        msa_dir=processed.msa_dir,
        mol_dir=canonical_mols_dir,
        num_workers=0,
        device_mesh=device_mesh,
        device_mesh_cpu=cp_device_mesh,
        local_batch_size=local_batch_size,
    )
    data_module_dtensor.setup("predict")

    # Get data loaders
    sampler = DistributedSampler(
        dataset,
        num_replicas=device_mesh.shape[0],
        rank=device_mesh.get_local_rank(0),
        shuffle=False,
        drop_last=False,
    )
    dataloader_serial = DataLoader(
        dataset,
        sampler=sampler,
        batch_size=local_batch_size,
        num_workers=0,
        shuffle=False,
        collate_fn=collate,
    )
    dataloader_dtensor = data_module_dtensor.predict_dataloader()

    # Process batches
    for data_batch_serial, data_batch_dtensor in zip(dataloader_serial, dataloader_dtensor):
        data_batch_serial = data_module_serial.transfer_batch_to_device(data_batch_serial, manager.device, 0)
        data_batch_dtensor = data_module_dtensor.transfer_batch_to_device(data_batch_dtensor, manager.device, 0)

        # Get masks for removing padding
        dp_idx_str = dp_rank * local_batch_size
        dp_idx_end = dp_idx_str + local_batch_size

        atom_pad_mask_dtensor_full = data_batch_dtensor["atom_pad_mask"].full_tensor()[dp_idx_str:dp_idx_end].bool()
        token_pad_mask_dtensor_full = data_batch_dtensor["token_pad_mask"].full_tensor()[dp_idx_str:dp_idx_end].bool()
        atom_pad_mask_serial = data_batch_serial["atom_pad_mask"].bool()
        token_pad_mask_serial = data_batch_serial["token_pad_mask"].bool()

        # Get atom_to_token matrices
        atom_to_token_serial = data_batch_serial["atom_to_token"]
        atom_to_token_dtensor = data_batch_dtensor["atom_to_token"]

        # Test 1: Single representation token-to-atom
        token_single_repr_serial = data_batch_serial["residue_index"].float()
        token_single_repr_dtensor = data_batch_dtensor["residue_index"].float()

        atom_single_repr_serial = torch.einsum(
            "bij,bj...->bi...", atom_to_token_serial.float(), token_single_repr_serial
        )
        atom_single_repr_dtensor = single_repr_token_to_atom(token_single_repr_dtensor, atom_to_token_dtensor)

        # Remove padding and compare
        for local_sample_idx in range(local_batch_size):
            atom_serial_sample = atom_single_repr_serial[local_sample_idx]
            atom_serial_sample = atom_serial_sample[atom_pad_mask_serial[local_sample_idx]]

            atom_dtensor_full = atom_single_repr_dtensor.full_tensor()[dp_idx_str:dp_idx_end]
            atom_dtensor_sample = atom_dtensor_full[local_sample_idx].float()
            atom_dtensor_sample = atom_dtensor_sample[atom_pad_mask_dtensor_full[local_sample_idx]]

            torch.testing.assert_close(
                atom_dtensor_sample,
                atom_serial_sample,
                atol=1e-6,
                rtol=1e-6,
                msg=f"Single representation token-to-atom mismatch on rank {cp_rank}, sample {local_sample_idx}",
            )

        # Test 2: Single representation atom-to-token
        atom_single_repr_serial = data_batch_serial["ref_charge"].float()
        atom_single_repr_dtensor = data_batch_dtensor["ref_charge"]

        # Non-sharded operation (with normalization)
        atom_to_token_normalized = atom_to_token_serial.float() / (atom_to_token_serial.sum(dim=1, keepdim=True) + 1e-6)
        token_single_repr_serial = torch.einsum("bji,bj...->bi...", atom_to_token_normalized, atom_single_repr_serial)
        token_single_repr_dtensor = single_repr_atom_to_token(atom_single_repr_dtensor.float(), atom_to_token_dtensor)

        # Remove padding and compare
        for local_sample_idx in range(local_batch_size):
            token_serial_sample = token_single_repr_serial[local_sample_idx]
            token_serial_sample = token_serial_sample[token_pad_mask_serial[local_sample_idx]]

            token_dtensor_full = token_single_repr_dtensor.full_tensor()[dp_idx_str:dp_idx_end]
            token_dtensor_sample = token_dtensor_full[local_sample_idx].float()
            token_dtensor_sample = token_dtensor_sample[token_pad_mask_dtensor_full[local_sample_idx]]

            torch.testing.assert_close(
                token_dtensor_sample,
                token_serial_sample,
                atol=1e-6,
                rtol=1e-6,
                msg=f"Single representation atom-to-token mismatch on rank {cp_rank}, sample {local_sample_idx}",
            )

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env",
    [((1, (2, 2)), True, "cuda", "ENV")],
    indirect=("setup_env",),
    ids=lambda x: f"dp={x[0][0]}, cp={x[0][1]}, specify_method={x[1]}, device_type={x[2]}, method_init={x[3]}",
)
def test_atom_to_token_sharding_consistency(
    setup_env: tuple[dict, int, str, str, str, dict[str, str]],
    test_cp_training_base_data_dir_boltz2: Path,
    canonical_mols_dir: Path,
    tmp_path: Path,
    local_batch_size: int = 2,
):
    """Test that sharded atom_to_token operations produce the same results as non-sharded after removing padding.

    This test verifies the sharding strategy described in atom_to_token.py where:
    1. Padding atoms/tokens are added to make dimensions divisible by cp_size
    2. The atom_to_token matrix becomes block diagonal with only diagonal blocks non-zero
    3. Row-wise broadcasting enables local computation on all ranks
    4. After removing padding, sharded and non-sharded results should be identical

    The test covers all three atom_to_token operations:
    - Single representation token-to-atom: residue_index -> atom representation (torch.einsum)
    - Single representation atom-to-token: ref_charge -> token representation (torch.einsum with normalization)
    """
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    # Create combined dataset with multiple samples to support local_batch_size > 1
    names = ["7ylz", "7z64", "8ayv", "8b2e"]
    combined_dataset_path = concat_data(
        tmp_path / "processed_combined",
        *[test_cp_training_base_data_dir_boltz2 / f"processed_{name}" for name in names],
    )

    # Check if we have enough samples for the requested batch size
    total_samples = len(names)
    if total_samples < local_batch_size:
        pytest.skip(
            f"Not enough samples ({total_samples}) for local_batch_size * dp_size ({local_batch_size * grid_group_sizes['dp']})"
        )

    # Create a new BoltzProcessedInput with the combined dataset
    combined_processed_handle = BoltzProcessedInput(
        manifest=Manifest.load(combined_dataset_path / "manifest.json"),
        targets_dir=combined_dataset_path / "structures",
        msa_dir=combined_dataset_path / "msa",
        template_dir=None,
        extra_mols_dir=None,
    )

    spawn_multiprocessing(
        parallel_assert_atom_to_token_sharding,
        world_size,
        combined_processed_handle,
        canonical_mols_dir,
        local_batch_size,
        device_type,
        backend,
        grid_group_sizes,
        env_per_rank,
    )


def parallel_assert_cp_training_dataloader(
    rank: int,
    training_data_dir: Path,
    canonical_mols_dir: Path,
    device_type,
    backend,
    grid_group_sizes: Dict[str, int],
    env_map: Optional[dict[str, str]] = None,
    dataloader_kind: str = "train",
):
    local_batch_size = 1
    monkeypatch = pytest.MonkeyPatch()
    if env_map is not None:
        for var_name, value in env_map.items():
            if value == "<INPUT_RANK>":
                monkeypatch.setenv(var_name, f"{rank}")
                continue
            monkeypatch.setenv(var_name, value)

    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    manager = DistributedManager()
    DistributedManager.create_group(
        "world_cpu", manager.group_ranks["world"], backend="gloo", use_local_synchronization=True
    )
    DistributedManager.create_group("cp_cpu", manager.group_ranks["cp"], backend="gloo", use_local_synchronization=True)

    device_mesh = manager.device_mesh_subgroups
    cp_device_mesh = map_subgroup_mesh_to_cpu(manager)
    dp_rank = device_mesh.get_local_rank(0)

    # Serial reference dataloader
    cfg_serial = setup_mock_training_datamodule_config(training_data_dir)
    cfg_serial.batch_size = local_batch_size
    cfg_serial.moldir = str(canonical_mols_dir)
    if dataloader_kind == "val":
        cfg_serial.overfit = 4
    else:
        cfg_serial.samples_per_epoch = 2 * local_batch_size * device_mesh.shape[0]
    seed_by_rank(0, seed=42)
    data_module_serial = Boltz2TrainingDataModule(cfg=cfg_serial)

    if dataloader_kind == "val":
        serial_set = data_module_serial._val_set
        serial_batch_size = cfg_serial.val_batch_size
    else:
        serial_set = data_module_serial._train_set
        serial_batch_size = cfg_serial.batch_size

    sampler_serial = DistributedSampler(
        serial_set,
        num_replicas=device_mesh.shape[0],
        rank=device_mesh.get_local_rank(0),
        shuffle=False,
        drop_last=False,
    )
    dataloader_serial = DataLoader(
        serial_set,
        sampler=sampler_serial,
        batch_size=serial_batch_size,
        num_workers=0,
        pin_memory=False,
        collate_fn=collate_training,
    )
    data_batch_list_serial = list(dataloader_serial)

    # DTensor distributed dataloader
    cfg_dtensor = setup_mock_training_datamodule_config(training_data_dir)
    cfg_dtensor.batch_size = local_batch_size
    cfg_dtensor.moldir = str(canonical_mols_dir)
    if dataloader_kind == "val":
        cfg_dtensor.overfit = 4
    else:
        cfg_dtensor.samples_per_epoch = cfg_serial.samples_per_epoch
    seed_by_rank(0, seed=42)
    data_module_dtensor = BoltzTrainingDataModuleDTensor(
        cfg=cfg_dtensor,
        device_mesh=device_mesh,
        device_mesh_cpu=cp_device_mesh,
    )
    if dataloader_kind == "val":
        dataloader_dtensor = data_module_dtensor.val_dataloader()
    else:
        dataloader_dtensor = data_module_dtensor.train_dataloader()
    data_batch_list_dtensor = list(dataloader_dtensor)

    skip_keys = (
        {
            "atom_pad_mask",
            "token_pad_mask",
            "msa_mask",
            "atom_to_token",
            "pair_mask",
        }
        | TRAINING_FEATURES_DIFFERENCE
        | NON_SHARDED_FEATURES_V2
    )

    n_shards = device_mesh.shape[1]

    for data_batch_serial, data_batch_dtensor in zip(data_batch_list_serial, data_batch_list_dtensor):
        serial_keys = set(data_batch_serial.keys())
        dtensor_keys = set(data_batch_dtensor.keys())
        serial_keys_filtered = serial_keys - TRAINING_FEATURES_DIFFERENCE
        dtensor_keys_filtered = dtensor_keys - TRAINING_FEATURES_DIFFERENCE
        if serial_keys_filtered != dtensor_keys_filtered:
            missing_in_dtensor = sorted(serial_keys_filtered - dtensor_keys_filtered)
            extra_in_dtensor = sorted(dtensor_keys_filtered - serial_keys_filtered)
            raise AssertionError(
                "Feature keys are different between serial and DTensor training pipelines. "
                f"rank={rank}, cp_rank={manager.group_rank['cp']}, "
                f"serial_pdb_id={data_batch_serial.get('pdb_id')}, "
                f"dtensor_pdb_id={data_batch_dtensor.get('pdb_id')}, "
                f"missing_in_dtensor={missing_in_dtensor}, "
                f"extra_in_dtensor={extra_in_dtensor}"
            )
        common_keys = sorted(serial_keys & dtensor_keys)

        data_batch_serial = data_module_serial.transfer_batch_to_device(data_batch_serial, manager.device, 0)
        data_batch_dtensor = data_module_dtensor.transfer_batch_to_device(data_batch_dtensor, manager.device, 0)

        dp_idx_str = dp_rank * local_batch_size
        dp_idx_end = dp_idx_str + local_batch_size

        _compare_features(
            common_keys=common_keys,
            skip_keys=skip_keys,
            data_batch_serial=data_batch_serial,
            data_batch_dtensor=data_batch_dtensor,
            dp_idx_str=dp_idx_str,
            dp_idx_end=dp_idx_end,
            n_shards=n_shards,
            rank=rank,
            manager=manager,
        )

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env",
    [
        ((1, (3, 3)), True, "cpu", "ENV"),
        ((2, (2, 2)), True, "cuda", "ENV"),
    ],
    indirect=("setup_env",),
    ids=lambda x: f"dp={x[0][0]}, cp={x[0][1]}, device_type={x[2]}",
)
@pytest.mark.parametrize("dataloader_kind", ["train", "val"])
def test_dtensor_cp_training_dataloader(
    setup_env: tuple[dict, int, str, str, str, dict[str, str]],
    test_cp_training_data_dir_boltz2: Path,
    canonical_mols_dir: Path,
    dataloader_kind: str,
):
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    spawn_multiprocessing(
        parallel_assert_cp_training_dataloader,
        world_size,
        test_cp_training_data_dir_boltz2,
        canonical_mols_dir,
        device_type,
        backend,
        grid_group_sizes,
        env_per_rank,
        dataloader_kind,
    )
