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


import warnings
from collections import OrderedDict

import torch
from torch.distributed import ProcessGroup
from torch.distributed.device_mesh import DeviceMesh
from torch.distributed.tensor import DTensor, Placement, Replicate, Shard

from boltz.data.pad import pad_dim
from boltz.distributed.data.feature.featurizer_utils import (
    _ATOM_IGNORE_FIELDS_FOR_CONTEXT_PARALLEL_PREDICT,
    remap_atom_indices_unpadded_to_padded,
)
from boltz.distributed.data.utils import TensorMetadata
from boltz.distributed.model.layers.shardwise_op import shardwise_argmax, shardwise_offset
from boltz.distributed.model.layers.squeeze import shardwise_unsqueeze
from boltz.distributed.model.layers.utils import distributed_pack_and_pad
from boltz.distributed.utils import LayoutRightMap


def pad_and_scatter_atom_features_dtensor(
    features: dict[str, torch.Tensor] | None,
    placements: dict[str, Placement],
    group: ProcessGroup,
    src_rank_global: int,
    device_mesh: DeviceMesh,
) -> dict[str, DTensor]:
    """Pad and distribute atom features as DTensors for the AtomTransformer in distributed setting.

    This function performs sophisticated sharding of molecular features across a 2D device mesh,
    handling both 1D (atom-based) and 2D (pair-based) features with appropriate padding and
    scattering. Tokens are sharded along the first dimension (i) and duplicated across the
    second dimension (j), ensuring balanced distribution of computational work.

    Parameters
    ----------
    features : dict[str, torch.Tensor] | None
        Dictionary mapping feature names to tensors containing molecular data. Must be provided
        on the source rank and None on all other ranks. Features may include random prefixes
        for testing ordering robustness. Key features include:
        - atom_counts_per_token: Number of atoms per token for sharding calculations
        - coords: Atom coordinates with leading ensemble dimension (E, A, 3)
        - pair_mask: Pairwise interaction masks requiring 2D sharding
        - atom_to_token, token_to_rep_atom: Mapping features with diagonal replication
        - frames_idx: Frame indices requiring cumulative padding adjustments
    placements : dict[str, Placement]
        Dictionary mapping feature names to their desired tensor placements. Must have the
        same keys as features on the source rank. Expected placements:
        - Single representation features: (Shard(0), Replicate()) or (Shard(1), Replicate()) for coords
        - 2D features: (Shard(0), Shard(1)) for pair_mask, (Shard(0), Replicate()) for mapping features
    group : ProcessGroup
        The distributed process group containing all ranks in the device mesh. Used for
        broadcasting metadata and scattering tensor shards.
    src_rank_global : int
        The global rank that serves as the source for the original feature tensors.
        This rank must be included in the process group.
    device_mesh : DeviceMesh
        A 2D square device mesh defining the distributed tensor layout. Currently only
        supports square meshes (n_rows == n_cols). The mesh ranks must match the process group ranks.

    Returns
    -------
    dict[str, DTensor]
        Dictionary mapping feature names to distributed tensors (DTensors) with appropriate
        sharding and padding applied. Features are processed in deterministic order to ensure
        consistency across ranks. Some features may be filtered out based on usage context.

    Raises
    ------
    ValueError
        - If placements is not a dict
        - If group is not a ProcessGroup
        - If device_mesh is not a DeviceMesh or not 2D/square
        - If device is not a torch.device or doesn't match device_mesh type
        - If features/placements key mismatch on source rank
        - If features is not None on non-source ranks
        - If process group ranks don't match device mesh ranks
        - If source rank is not in the process group
        - If number of tokens is not divisible by shard dimension
        - If placement dimensions don't match expected patterns for each feature type

    Notes
    -----
    - Only supports 2D square device meshes currently
    - Filters out fields in _ATOM_IGNORE_FIELDS_FOR_CONTEXT_PARALLEL_PREDICT (training/confidence-only)
    - Uses OrderedDict internally to ensure consistent feature iteration order across ranks
    - Handles complex padding logic for both 1D and 2D features:
      * 1D features: Padded to max_atoms_per_shard along atom dimension
      * 2D features: Complex 2D padding with different strategies per feature type
    - Special handling for coordinates (leading ensemble dimension) and frames_idx (cumulative padding)
    - Broadcasts tensor metadata before scattering to ensure consistent shard shapes
    - Each feature requires specific placement patterns based on its semantic meaning

    Examples
    --------
    >>> features = {"coords": torch.randn(1, 100, 3), "atom_counts_per_token": torch.ones(100)}
    >>> placements = {"coords": (Shard(1), Replicate()), "atom_counts_per_token": (Shard(0), Replicate())}
    >>> dtensors = pad_and_scatter_atom_features_dtensor(features, placements, group, 0, mesh)
    """
    rank_global = torch.distributed.get_rank()
    is_src_rank = rank_global == src_rank_global

    if not isinstance(placements, dict):
        raise ValueError(f"Placements must be a dict, got {placements}")

    if not isinstance(group, ProcessGroup):
        raise ValueError(f"Group must be a ProcessGroup, got {group}")

    if not isinstance(device_mesh, DeviceMesh):
        raise ValueError(f"Device mesh must be a DeviceMesh, got {device_mesh}")

    if is_src_rank:
        if not isinstance(features, dict):
            raise ValueError(f"Features must be a dict on source rank {rank_global}, got {features}")
        if features.keys() != placements.keys():
            raise ValueError(
                f"Features and placements must have the same keys, got {features.keys()} and {placements.keys()}"
            )
    else:
        if features is not None:
            raise ValueError(f"Features must be None on non-source rank {rank_global}, got {features}")

    # check consistency of process group and device mesh
    ranks_in_group = torch.distributed.get_process_group_ranks(group)

    ranks_in_mesh = device_mesh.mesh.clone()  # for inplace modification later
    if ranks_in_group != ranks_in_mesh.flatten().tolist():
        raise ValueError(
            f"Ranks in group {ranks_in_group} do not match ranks in mesh {ranks_in_mesh}, got {ranks_in_group} and {ranks_in_mesh}"
        )

    if src_rank_global not in ranks_in_group:
        raise ValueError(f"Source rank {src_rank_global} not in group {ranks_in_group}")

    # we hardcode the sharding of token and atom features along a square sub-mesh of the device mesh
    if device_mesh.ndim != 2:
        raise ValueError(f"Only 2D device meshes currently supported, got {device_mesh.ndim}")
    n_rows, n_cols = device_mesh.shape
    if n_rows != n_cols:
        raise ValueError(f"Only square device grids currently supported, got {n_rows},{n_cols}")

    # note that ranks_in_mesh might not be consecutive so we need a dictionary
    # this is equivalent to
    # rank_global_to_idx_scatter_list = {
    #     r_global: torch.distributed.get_group_rank(r_global, group) for r_global in ranks_in_mesh
    # }
    rank_global_to_idx_scatter_list = {r: i for i, r in enumerate(ranks_in_mesh.flatten().tolist())}
    # mesh_coord_to_idx_scatter_list[i_row, j_col] -> idx in the group
    # Later in the loop, while we iterate the shards in a LayoutRight order,
    # the mesh_coord_to_idx_scatter_list mapping always respect the layout of the device mesh
    mesh_coord_to_idx_scatter_list = ranks_in_mesh.apply_(rank_global_to_idx_scatter_list.get)

    # torch init_device_mesh and DeviceMesh ctor will do something like torch.cuda.set_device(rank)
    # so the associated device is rank specific, which we can rely on to set the device for the resulting
    # tensors
    device = torch.device(device_mesh.device_type)

    if is_src_rank:
        # metadata only relevant to source rank
        token_atom_counts = features["atom_counts_per_token"]

        N_tokens_total = token_atom_counts.shape[0]
        if N_tokens_total % n_rows != 0:
            raise ValueError(f"Number of tokens ({N_tokens_total}) is not divisible by shard dimension ({n_rows})")

        # atom_counts_per_token is padded to CP-divisible length by the
        # caller, but atom_to_token's token dimension (dim 1) may lag behind
        # because its placement Shard(0) only pads the atom dimension.
        if "atom_to_token" in features and features["atom_to_token"].shape[-1] < N_tokens_total:
            features["atom_to_token"] = pad_dim(
                features["atom_to_token"],
                features["atom_to_token"].ndim - 1,
                N_tokens_total - features["atom_to_token"].shape[-1],
            )

        token_atom_count_cumsum = torch.concatenate(
            (torch.tensor([0], device=token_atom_counts.device), torch.cumsum(token_atom_counts, dim=0))
        )
        N_tokens_per_shard = N_tokens_total // n_rows

        shard_atom_counts_token = token_atom_counts.unflatten(dim=0, sizes=(n_rows, N_tokens_per_shard)).sum(dim=1)
        max_atoms_per_shard_token = shard_atom_counts_token.amax().item()

        # max_atoms_per_shard for consistent intersperse padding
        max_atoms_per_shard = max_atoms_per_shard_token

        # Pre-compute r_set_to_rep_atom metadata (loop-independent, computed once)
        # r_set_to_rep_atom: [size_r_set, N_atoms] one-hot
        #
        # CRITICAL: R-set sharding must ALIGN with token sharding!
        # Each R-set element corresponds to a specific token (via its rep atom).
        # We must put each R-set element in the same shard as its token.
        #
        # The co-sharding of N_tokens and N_R is bijective:
        # - Forward: R-set elements of shard i are a SUBSET of tokens of shard i.
        #   No R-set element can correspond to a token in a different shard.
        # - Inverse: Tokens of shard i are a SUPERSET of its R-set elements.
        #   This holds by definition since each R element maps to exactly one token.
        #
        # PURPOSE: This co-sharding enables LOCAL matmul for mapping coordinates
        # from atom space to R-set space (r_coords = r_set_to_rep_atom @ atom_coords).
        # This mirrors the atom_to_token co-sharding strategy where atoms and tokens
        # are co-sharded as diagonal blocks, enabling local atom-to-token mappings.
        #
        # Algorithm:
        # 1. Compute r_set_to_token = r_set_to_rep_atom @ atom_to_token
        # 2. Determine which shard each R-set element belongs to
        # 3. Group R-set elements by their token shard
        # 4. Each shard (i, *) gets only R-set elements for tokens in shard i
        r_set_precomputed = None
        if "r_set_to_rep_atom" in features and "atom_to_token" in features:
            r_set_to_rep_atom_v = features["r_set_to_rep_atom"]
            # Get valid R-set elements (filter out padding rows - all zeros)
            r_set_valid_mask = r_set_to_rep_atom_v.any(dim=-1)  # [size_r_set_total]
            r_set_valid = r_set_to_rep_atom_v[r_set_valid_mask]  # [size_r_set_valid, N_atom_global]
            size_r_set_valid = r_set_valid.shape[0]

            if size_r_set_valid == 0:
                r_set_precomputed = {"size_r_set_valid": 0}
            else:
                # Get atom_to_token for computing r_set -> token mapping.
                # The general padding loop may have padded dim 0 (atoms) for
                # Shard(0), but r_set_valid has the original atom count in
                # dim 1.  Slice to match.
                n_atoms_r_set = r_set_valid.shape[-1]
                atom_to_token = features["atom_to_token"][:n_atoms_r_set]
                # Both r_set_valid and atom_to_token are one-hot, so the matmul
                # r_set_valid @ atom_to_token is equivalent to two index lookups:
                #   1. argmax over r_set_valid rows -> rep atom index per r_set element
                #   2. argmax over atom_to_token rows -> token index per atom
                # This avoids an O(N_r * N_atoms * N_tokens) int64 scalar loop
                # (no BLAS path for int64) and replaces it with two O(N * M) argmax passes.
                r_set_rep_atom_idx = r_set_valid.argmax(dim=-1)  # [size_r_set_valid]
                atom_to_token_idx = atom_to_token.argmax(dim=-1)  # [N_atoms]
                r_set_token_idx = atom_to_token_idx[r_set_rep_atom_idx]  # [size_r_set_valid]

                # Determine which shard each R-set element belongs to
                r_set_shard_idx = r_set_token_idx // N_tokens_per_shard  # [size_r_set_valid]

                # Count R-set elements per shard using scatter_add
                ones = torch.ones(size_r_set_valid, dtype=torch.long, device=device)
                shard_counts = torch.zeros(n_rows, dtype=torch.long, device=device)
                shard_counts.scatter_add_(0, r_set_shard_idx, ones)

                # Max R-set size per shard for padding
                max_size_r_set_per_shard = max(1, shard_counts.max().item())

                r_set_precomputed = {
                    "size_r_set_valid": size_r_set_valid,
                    "r_set_valid": r_set_valid,
                    "r_set_shard_idx": r_set_shard_idx,
                    "max_size_r_set_per_shard": max_size_r_set_per_shard,
                }

    _2D_feats = {
        "pair_mask",
        "atom_to_token",
        "token_to_rep_atom",
        "frames_idx",
        "r_set_to_rep_atom",  # 2D one-hot [N_R, N_atoms]
    }

    _2D_features_placement_as_single = {
        "atom_to_token",
        "token_to_rep_atom",
        "frames_idx",
        "r_set_to_rep_atom",  # Diagonal block sharding like token_to_rep_atom
    }

    # loop over each feature, create a list of tensors, scatter them then call DTensor.from_local
    # To guarantee the order of iterating thru the features in the dictionary, we need to
    # convert placements to a OrderedDict first so that all ranks go thru the keys in the same order.
    placements_ordered = OrderedDict(sorted(placements.items()))

    result: dict[str, DTensor] = {}
    for k, placement in placements_ordered.items():
        if k in _ATOM_IGNORE_FIELDS_FOR_CONTEXT_PARALLEL_PREDICT:
            continue

        if k == "atom_counts_per_token":
            # Only used in sharding preprocessing.
            continue

        if len(placement) != device_mesh.ndim:
            raise ValueError(
                f"Placement for {k} has {len(placement)} dimensions, expected {device_mesh.ndim} from the device_mesh"
            )

        is_single_repr = k not in _2D_feats
        if is_single_repr:
            if k == "coords":
                placement_expected = (Shard(1), Replicate())
            else:
                placement_expected = (Shard(0), Replicate())
        else:
            if k == "frames_idx":
                placement_expected = (Shard(1), Replicate())
            elif k in _2D_features_placement_as_single:
                placement_expected = (Shard(0), Replicate())
            else:
                placement_expected = (Shard(0), Shard(1))

        if is_src_rank and k == "frames_idx":
            v_src = features[k]
            if v_src.ndim != 3:
                raise ValueError(
                    "frames_idx must have v2 ensemble-aware shape (E, T, 3) "
                    f"with ndim=3, got ndim={v_src.ndim} with shape={tuple(v_src.shape)}"
                )

        if is_src_rank and k == "frame_resolved_mask":
            v_src = features[k]
            if v_src.ndim != 2:
                raise ValueError(
                    "frame_resolved_mask must have v2 ensemble-aware shape (E, T) "
                    f"with ndim=2, got ndim={v_src.ndim} with shape={tuple(v_src.shape)}"
                )

        placement_valid = placement == placement_expected
        if not placement_valid:
            raise ValueError(f"Placement for {k} is {placement}, expected {placement_expected} from the device_mesh")

        if is_src_rank:
            v = features[k]
            # create a list of tensors on the src rank
            scatter_list = [None] * ranks_in_mesh.numel()

            for i in range(n_rows):
                # Entries duplicated over j
                if k not in _2D_feats:
                    # single representation
                    token_start = N_tokens_per_shard * i
                    token_end = N_tokens_per_shard * (i + 1)
                    atom_start = token_atom_count_cumsum[token_start]
                    atom_end = token_atom_count_cumsum[token_end]
                    num_atoms_in_range = atom_end - atom_start
                    pad_amount = max_atoms_per_shard - num_atoms_in_range

                    if k == "coords":
                        # Leading dimension is ensemble count (E, A, 3); shard and pad along atom dim (1).
                        j_duplicates_val = pad_dim(v[:, atom_start:atom_end, ...], 1, pad_amount)
                    else:
                        j_duplicates_val = pad_dim(v[atom_start:atom_end], 0, pad_amount)
                    for j in range(n_cols):
                        # TODO: see if we can avoid the clone here
                        if j_duplicates_val.dtype in [torch.int8, torch.bool]:
                            j_duplicates_val = j_duplicates_val.clone()
                        scatter_list[mesh_coord_to_idx_scatter_list[i, j].item()] = j_duplicates_val
                else:
                    # pair representation
                    # find token and atom ranges for 2d padding
                    row_token_start = N_tokens_per_shard * i
                    row_token_end = N_tokens_per_shard * (i + 1)
                    shard_atom_start = token_atom_count_cumsum[row_token_start]
                    shard_atom_end = token_atom_count_cumsum[row_token_end]
                    shard_atoms_in_range = shard_atom_end - shard_atom_start

                    # 2D entries need separate calculation for each i, j
                    for j in range(n_cols):
                        if k == "pair_mask":
                            col_token_start = N_tokens_per_shard * j
                            col_token_end = N_tokens_per_shard * (j + 1)

                            # 2D indexing and padding, (atoms * atoms)
                            res = torch.zeros(
                                size=(max_atoms_per_shard, max_atoms_per_shard), dtype=v.dtype, device=device
                            )
                            col_atom_start = token_atom_count_cumsum[col_token_start]
                            col_atom_end = token_atom_count_cumsum[col_token_end]
                            col_atoms_in_range = col_atom_end - col_atom_start
                            res[:shard_atoms_in_range, :col_atoms_in_range] = v[
                                shard_atom_start:shard_atom_end, col_atom_start:col_atom_end
                            ]
                        elif k == "atom_to_token":
                            # 2D indexing and padding, (atoms * tokens), internal padding only needed in atom dim.
                            # NOTE: Each j column gets the diagonal representation (i,i) - so columns are i-based here, j is ignored
                            # except for computing the output shard.
                            col_token_start = N_tokens_per_shard * i
                            col_token_end = N_tokens_per_shard * (i + 1)
                            col_tokens_in_range = col_token_end - col_token_start
                            res = torch.zeros(
                                size=(max_atoms_per_shard, N_tokens_per_shard), dtype=v.dtype, device=device
                            )
                            res[:shard_atoms_in_range, :col_tokens_in_range] = v[
                                shard_atom_start:shard_atom_end, col_token_start:col_token_end
                            ]
                        elif k == "frames_idx":
                            # frames_idx shape is (E, T, 3): E ensembles, T tokens,
                            # and for each token 3 global atom indices that define
                            # its local coordinate frame.  The atom indices refer to
                            # positions in the unpadded, unsharded atom array.
                            # After sharding, each shard is padded to
                            # max_atoms_per_shard, shifting atom positions.
                            frame_idx = v[:, row_token_start:row_token_end, :]
                            res = remap_atom_indices_unpadded_to_padded(
                                frame_idx, shard_atom_counts_token, max_atoms_per_shard
                            )
                        elif k == "token_to_rep_atom":
                            # 2D indexing and padding, (tokens * atoms), internal padding only needed in atom dim. Similar to atom_to_token.
                            # NOTE: Each j column gets the diagonal representation (i,i) - so columns are i-based here, j is ignored
                            # except for computing the output shard.
                            col_token_start = N_tokens_per_shard * i
                            col_token_end = N_tokens_per_shard * (i + 1)
                            col_tokens_in_range = col_token_end - col_token_start
                            res = torch.zeros(
                                size=(N_tokens_per_shard, max_atoms_per_shard), dtype=v.dtype, device=device
                            )
                            res[:col_tokens_in_range, :shard_atoms_in_range] = v[
                                col_token_start:col_token_end,
                                shard_atom_start:shard_atom_end,
                            ]
                        elif k == "r_set_to_rep_atom":
                            # Use pre-computed metadata (computed once before loops)
                            # See r_set_precomputed initialization for full documentation on:
                            # - Co-sharding of N_tokens and N_R (bijective relationship)
                            # - Purpose: enabling LOCAL matmul for atom->R-set coordinate mapping
                            # - Algorithm for token-aligned R-set sharding

                            if r_set_precomputed["size_r_set_valid"] == 0:
                                # No valid R-set elements - create empty shard
                                res = torch.zeros((1, max_atoms_per_shard), dtype=v.dtype, device=device)
                            else:
                                r_set_valid = r_set_precomputed["r_set_valid"]
                                r_set_shard_idx = r_set_precomputed["r_set_shard_idx"]
                                max_size_r_set_per_shard = r_set_precomputed["max_size_r_set_per_shard"]

                                # Get indices of R-set elements for shard i
                                shard_mask = r_set_shard_idx == i  # [size_r_set_valid]
                                r_set_in_shard = r_set_valid[shard_mask]  # [count_i, N_atom_global]
                                count_i = r_set_in_shard.shape[0]

                                # Create output: [max_size_r_set_per_shard, max_atoms_per_shard]
                                res = torch.zeros(
                                    (max_size_r_set_per_shard, max_atoms_per_shard), dtype=v.dtype, device=device
                                )
                                if count_i > 0 and shard_atoms_in_range > 0:
                                    # Slice atoms for this token shard i (diagonal block)
                                    res[:count_i, :shard_atoms_in_range] = r_set_in_shard[
                                        :, shard_atom_start:shard_atom_end
                                    ]
                        scatter_list[mesh_coord_to_idx_scatter_list[i, j].item()] = res

        else:
            scatter_list = None

        # broadcast the metadata
        # Assumption: all shards in scatter_list have the same shape, which is the assumption made by the
        # code blocks above
        l_metadata = [TensorMetadata(dtype=v.dtype, shape=scatter_list[0].shape)] if is_src_rank else [None]
        torch.distributed.broadcast_object_list(l_metadata, src=src_rank_global, group=group, device=device)

        # scatter the tensor
        local_shard = torch.empty(l_metadata[0].shape, dtype=l_metadata[0].dtype, device=device)
        torch.distributed.scatter(local_shard, scatter_list, src=src_rank_global, group=group)

        # create the DTensor from local shard
        # Due to the padding, we need to recompute the global shape with padding applied
        shape_global = list(l_metadata[0].shape)
        for i_dim_mesh, p in enumerate(placement):
            if isinstance(p, Shard):
                shape_global[p.dim] *= device_mesh.shape[i_dim_mesh]
        shape_global = tuple(shape_global)
        # Due to the local_shard from torch.empty and data shared by torch.distributed.scatter,
        # the stride is guaranteed to be LayoutRight, i.e., torch's 'contiguous' memory layout
        stride_global = LayoutRightMap(shape=shape_global).strides
        dtensor = DTensor.from_local(
            local_shard, device_mesh, placements=placement, shape=shape_global, stride=stride_global
        )
        result[k] = dtensor

    return result


def pack_atom_features(
    feats: dict[str, DTensor],
    keys_subset: set[str],
    W: int,
) -> OrderedDict[str, DTensor]:
    """Pack and pad atom features using distributed_pack_and_pad.

    This removes per-shard trailing padding from pad_and_scatter_atom_features_dtensor
    and creates a packed DTensor with global trailing padding (multiple of W * size_cp).

    The function handles keys in keys_subset as follows:
    - "atom_to_token": Special handling - converts shard-local indices to global before packing,
      and stores the global indices as "atom_to_token_ids_global". NOTE: The original
      'atom_to_token' one-hot matrix is NOT packed. Directly packing 'atom_to_token' would
      give inconsistent sharding scheme between atom and token dimensions, making the packed
      'atom_to_token' not useful in practice. Only the global indices are packed and returned.
    - All other keys: Treated as generic atom features and packed directly with the mask

    Parameters
    ----------
    feats : dict[str, DTensor]
        Dictionary of atom features as DTensors. Must contain "atom_pad_mask" key
        with shape (B, N_atoms) to use as the packing mask.
    keys_subset : set[str]
        Set of keys to pack from feats. Must contain "atom_pad_mask".
        Only keys present in both feats and keys_subset will be processed.
        All keys (except "atom_to_token" and "atom_pad_mask") are treated as
        generic atom features with N_atoms axis at position 1.
    W : int
        Window size for packing (atoms per window for queries).
        The packed output will have length that is a multiple of W * size_cp.

    Returns
    -------
    OrderedDict[str, DTensor]
        OrderedDict of packed atom features with keys in sorted order. Contains:
        - All keys from keys_subset that were in feats, with packed values
        - "atom_to_token_ids_global" if "atom_to_token" was in keys_subset
        Keys are sorted to ensure consistent iteration order across distributed ranks.

    Raises
    ------
    ValueError
        If feats does not contain "atom_pad_mask" key.
        If feats["atom_pad_mask"] does not have ndim=2 (expected shape: B, N_atoms).
        If keys_subset does not contain "atom_pad_mask".
        If keys_subset contains keys not present in feats.
    NotImplementedError
        If keys_subset contains "coords" (packing coords is not supported).
    """
    if "atom_pad_mask" not in feats:
        raise ValueError("feats must contain 'atom_pad_mask' key")
    if feats["atom_pad_mask"].ndim != 2:
        raise ValueError(
            f"feats['atom_pad_mask'] must have ndim=2 (B, N_atoms), got ndim={feats['atom_pad_mask'].ndim}"
        )
    if "atom_pad_mask" not in keys_subset:
        raise ValueError("keys_subset must contain 'atom_pad_mask'")
    if "coords" in keys_subset:
        raise NotImplementedError("packing 'coords' is not supported")

    # Verify all keys in keys_subset are present in feats
    missing_keys = keys_subset - feats.keys()
    if missing_keys:
        raise ValueError(f"keys_subset contains keys not in feats: {missing_keys}")

    # Sort keys to ensure consistent iteration order across all ranks
    # This is critical because distributed collective operations must be called in the same order
    keys_sorted = sorted(keys_subset)

    # Pack and pad each atom feature in keys_subset
    # Use no_grad to prevent backprop - these are input features not supposed to receive gradients
    feats_dt_packed: OrderedDict[str, DTensor] = OrderedDict()
    with torch.no_grad():
        # Get atom_mask for pack_and_pad (shape: B, N_atoms)
        # Must convert to bool() to match pack_and_pad behavior and avoid
        # precision issues when summing bfloat16 mask values
        atom_mask_dt = feats["atom_pad_mask"].bool()

        for key in keys_sorted:
            # Verify feature shape matches mask shape in first two dimensions (B, N_atoms)
            if feats[key].shape[:2] != atom_mask_dt.shape:
                raise ValueError(
                    f"feats['{key}'].shape[:2]={feats[key].shape[:2]} must match atom_mask_dt.shape={atom_mask_dt.shape}"
                )

            if key == "atom_to_token":
                # Special handling: convert shard-local indices to global BEFORE packing.
                # This is necessary because distributed_pack_and_pad may move atoms between shards,
                # which would result in different sharding schemes between token and atom.
                # On the other hand, shardwise_offset below assumes consistent sharding scheme
                # between token and atom, i.e., any rank must own all atoms from all its own
                # tokens, or any token's atom collection is not sharded. Therefore,
                # shardwise_offset must be applied before distributed_pack_and_pad.
                # 1. Get shard-local token indices from block-diagonal atom_to_token
                atom_to_token_dt = feats[key]  # (B, N_atoms, N_tokens_per_shard)
                if atom_to_token_dt.ndim != 3:
                    raise ValueError(
                        f"feats['atom_to_token'] must have ndim=3 (B, N_atoms, N_tokens_per_shard), "
                        f"got ndim={atom_to_token_dt.ndim}"
                    )
                n_tokens_per_shard = atom_to_token_dt.to_local().shape[2]
                atom_to_token_ids_local = shardwise_argmax(atom_to_token_dt, dim=-1, keepdim=False)
                # 2. Convert to global indices
                atom_to_token_ids_global = shardwise_offset(
                    atom_to_token_ids_local, dim=1, offset_per_rank=n_tokens_per_shard
                )
                # 3. Pack the global indices (not the one-hot matrix)
                mask_for_ids = atom_mask_dt
                packed_ids, _ = distributed_pack_and_pad(atom_to_token_ids_global, mask_for_ids, W, axis=1)
                feats_dt_packed["atom_to_token_ids_global"] = packed_ids
                feats_dt_packed["atom_to_token_local_onehot"] = atom_to_token_dt
                if not getattr(pack_atom_features, "_warned_a2t", False):
                    pack_atom_features._warned_a2t = True
                    warnings.warn(
                        "pack_atom_features: 'atom_to_token' one-hot matrix is NOT packed but returned "
                        "as-is together with 'atom_to_token_ids_global' for window batching AtomAttentionEncoder. "
                        "Directly packing 'atom_to_token' would give inconsistent sharding scheme "
                        "between atom and token, making it not useful in practice.",
                        stacklevel=2,
                    )
            else:
                # Generic atom feature: expand mask and pack
                feat = feats[key]
                mask_for_feat = atom_mask_dt
                # Add trailing dimensions to mask for broadcasting with feat
                while mask_for_feat.ndim < feat.ndim:
                    mask_for_feat = shardwise_unsqueeze(mask_for_feat, -1)
                packed_feat, _packed_mask = distributed_pack_and_pad(feat, mask_for_feat, W, axis=1)
                feats_dt_packed[key] = packed_feat

    return feats_dt_packed
