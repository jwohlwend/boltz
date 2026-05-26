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

from boltz.distributed.data.utils import (
    PLACEMENT_TYPE_SHARD0_REPLICATE,
    PLACEMENT_TYPE_SHARD0_SHARD1,
    PLACEMENT_TYPE_SHARD1_REPLICATE,
)

BASE_FEATURE_PLACEMENTS_V2: dict[str, tuple] = {
    # Atom features
    "ref_pos": PLACEMENT_TYPE_SHARD0_REPLICATE,
    "ref_charge": PLACEMENT_TYPE_SHARD0_REPLICATE,
    "atom_resolved_mask": PLACEMENT_TYPE_SHARD0_REPLICATE,
    "ref_element": PLACEMENT_TYPE_SHARD0_REPLICATE,
    "ref_atom_name_chars": PLACEMENT_TYPE_SHARD0_REPLICATE,
    "ref_space_uid": PLACEMENT_TYPE_SHARD0_REPLICATE,
    "coords": PLACEMENT_TYPE_SHARD1_REPLICATE,
    "atom_counts_per_token": PLACEMENT_TYPE_SHARD0_REPLICATE,
    "frame_resolved_mask": PLACEMENT_TYPE_SHARD1_REPLICATE,
    "frames_idx": PLACEMENT_TYPE_SHARD1_REPLICATE,
    "atom_pad_mask": PLACEMENT_TYPE_SHARD0_REPLICATE,
    "atom_to_token": PLACEMENT_TYPE_SHARD0_REPLICATE,
    "token_to_rep_atom": PLACEMENT_TYPE_SHARD0_REPLICATE,
    "r_set_to_rep_atom": PLACEMENT_TYPE_SHARD0_REPLICATE,
    "ref_chirality": PLACEMENT_TYPE_SHARD0_REPLICATE,
    "atom_backbone_feat": PLACEMENT_TYPE_SHARD0_REPLICATE,
    "bfactor": PLACEMENT_TYPE_SHARD0_REPLICATE,
    "plddt": PLACEMENT_TYPE_SHARD0_REPLICATE,
    # Token features
    "token_index": PLACEMENT_TYPE_SHARD0_REPLICATE,
    "residue_index": PLACEMENT_TYPE_SHARD0_REPLICATE,
    "asym_id": PLACEMENT_TYPE_SHARD0_REPLICATE,
    "entity_id": PLACEMENT_TYPE_SHARD0_REPLICATE,
    "sym_id": PLACEMENT_TYPE_SHARD0_REPLICATE,
    "mol_type": PLACEMENT_TYPE_SHARD0_REPLICATE,
    "res_type": PLACEMENT_TYPE_SHARD0_REPLICATE,
    "disto_center": PLACEMENT_TYPE_SHARD0_REPLICATE,
    "disto_target": PLACEMENT_TYPE_SHARD0_SHARD1,
    "disto_coords_ensemble": PLACEMENT_TYPE_SHARD1_REPLICATE,
    "token_bonds": PLACEMENT_TYPE_SHARD0_SHARD1,
    "type_bonds": PLACEMENT_TYPE_SHARD0_SHARD1,
    "token_pad_mask": PLACEMENT_TYPE_SHARD0_REPLICATE,
    "token_resolved_mask": PLACEMENT_TYPE_SHARD0_REPLICATE,
    "token_disto_mask": PLACEMENT_TYPE_SHARD0_REPLICATE,
    "token_pair_pad_mask": PLACEMENT_TYPE_SHARD0_SHARD1,
    "pair_mask": PLACEMENT_TYPE_SHARD0_SHARD1,
    "contact_conditioning": PLACEMENT_TYPE_SHARD0_SHARD1,
    "contact_threshold": PLACEMENT_TYPE_SHARD0_SHARD1,
    "cyclic_period": PLACEMENT_TYPE_SHARD0_REPLICATE,
    "method_feature": PLACEMENT_TYPE_SHARD0_REPLICATE,
    "modified": PLACEMENT_TYPE_SHARD0_REPLICATE,
    # MSA features
    "msa": PLACEMENT_TYPE_SHARD0_SHARD1,
    "msa_paired": PLACEMENT_TYPE_SHARD0_SHARD1,
    "deletion_value": PLACEMENT_TYPE_SHARD0_SHARD1,
    "has_deletion": PLACEMENT_TYPE_SHARD0_SHARD1,
    "msa_mask": PLACEMENT_TYPE_SHARD0_SHARD1,
    "deletion_mean": PLACEMENT_TYPE_SHARD0_REPLICATE,
    "profile": PLACEMENT_TYPE_SHARD0_REPLICATE,
}

TRAINING_FEATURE_PLACEMENTS_V2: dict[str, tuple] = {
    **BASE_FEATURE_PLACEMENTS_V2,
    "temp_feature": PLACEMENT_TYPE_SHARD0_REPLICATE,
    "ph_feature": PLACEMENT_TYPE_SHARD0_REPLICATE,
}

INFERENCE_FEATURE_PLACEMENTS_V2: dict[str, tuple] = {
    **BASE_FEATURE_PLACEMENTS_V2,
    "affinity_token_mask": PLACEMENT_TYPE_SHARD0_REPLICATE,
}
