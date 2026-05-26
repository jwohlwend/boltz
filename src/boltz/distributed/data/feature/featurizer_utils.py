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

from functools import partial

import torch
from torch import Tensor

from boltz.data.types import Tokenized
from boltz.model.modules.encodersv2 import get_indexing_matrix, single_to_keys

# Fields to ignore during context parallel prediction
_ATOM_IGNORE_FIELDS_FOR_CONTEXT_PARALLEL_PREDICT: set[str] = set()

# Features whose values are atom indices into a padded global atom array.
# Both pad_and_scatter_atom_features_dtensor and CollateDTensor must remap
# these values whenever atom padding changes.
ATOM_INDEX_FEATURES = {"frames_idx"}


def remap_atom_indices_unpadded_to_padded(
    indices: Tensor,
    shard_atom_counts: Tensor,
    padded_atoms_per_shard: int,
) -> Tensor:
    """Remap atom indices from an unpadded global atom array to a padded layout.

    In the padded layout each CP shard occupies exactly ``padded_atoms_per_shard``
    positions, with trailing zeros filling unused slots.  This shifts atom
    positions relative to the dense (unpadded) layout.

    Parameters
    ----------
    indices : Tensor
        Atom indices into the unpadded global atom array (arbitrary shape).
    shard_atom_counts : Tensor
        Actual (unpadded) atom count per shard, shape ``(n_shards,)``.
    padded_atoms_per_shard : int
        Fixed size of each shard in the padded layout.

    Returns
    -------
    Tensor
        Indices remapped to the padded layout, same shape as *indices*.
    """
    dtype = indices.dtype
    # cumsum promotes int32→int64 to avoid overflow; cast back explicitly.
    shard_starts = torch.cat(
        [
            torch.zeros(1, device=shard_atom_counts.device, dtype=dtype),
            shard_atom_counts.cumsum(dim=0).to(dtype),
        ]
    )
    pad_offsets = (
        torch.arange(len(shard_atom_counts), device=shard_atom_counts.device, dtype=dtype) * padded_atoms_per_shard
        - shard_starts[:-1]
    )
    shard_idx = torch.bucketize(indices, shard_starts[1:], right=True)
    return indices + pad_offsets[shard_idx]


def remap_atom_indices_repad(
    indices: Tensor,
    old_atoms_per_shard: int,
    new_atoms_per_shard: int,
) -> Tensor:
    """Remap atom indices when per-shard padding changes.

    After ``pad_and_scatter_atom_features_dtensor``, atom-index features
    reference a padded layout with stride ``old_atoms_per_shard``.  When
    collation pads the atom dimension further (to align a batch or across DP
    ranks), the stride grows to ``new_atoms_per_shard`` and every stored index
    must be adjusted.

    Parameters
    ----------
    indices : Tensor
        Atom indices in the old padded layout (arbitrary shape).
    old_atoms_per_shard : int
        Per-shard atom count in the current (old) padded layout.
    new_atoms_per_shard : int
        Per-shard atom count in the target (new) padded layout.

    Returns
    -------
    Tensor
        Indices remapped to the new padded layout, same shape as *indices*.
    """
    if old_atoms_per_shard == new_atoms_per_shard:
        return indices
    shard_of_atom = indices // old_atoms_per_shard
    offset_in_shard = indices % old_atoms_per_shard
    return shard_of_atom * new_atoms_per_shard + offset_in_shard


def get_pair_mask(N_atoms: int, W: int = 32, H: int = 128) -> Tensor:
    """Get the pair mask for the atom transformer.

    Parameters
    ----------
    N_atoms : int
        The number of atoms.
    W : int, optional
        The attention window queries, by default 32.
    H : int, optional
        The attention window keys, by default 128.

    Returns
    -------
    Tensor
        The pair mask.

    """
    mask = torch.zeros(N_atoms, N_atoms)

    if N_atoms % W == 0:
        max_atoms = N_atoms
    else:
        # pad to the next multiple of W
        max_atoms = ((N_atoms // W) + 1) * W

    # construct pair mask through indexing matrices
    # TODO construct pair mask directly from AF3 appendix
    index = torch.arange(1, max_atoms + 1)
    index[N_atoms:] = 0
    index = index.unsqueeze(0)

    K = max_atoms // W
    keys_indexing_matrix = get_indexing_matrix(K, W, H, index.device)
    to_keys = partial(single_to_keys, indexing_matrix=keys_indexing_matrix, W=W, H=H)

    index_queries = index.view(K, W)
    index_keys = to_keys(index.unsqueeze(-1).float()).view(K, H).long()

    for index_query, index_key in zip(index_queries, index_keys):
        index_query = index_query[index_query != 0]
        index_key = index_key[index_key != 0]
        mask[index_query.min() - 1 : index_query.max(), index_key.min() - 1 : index_key.max()] = 1

    return mask


def tokenized_stats(
    tokenized: Tokenized,
) -> dict[str, int]:
    """Get statistics about the tokenized data.

    Parameters
    ----------
    tokenized : Tokenized
        The tokenized data.

    Returns
    -------
    dict[str, int]
        Dictionary containing:
        - num_atoms_total: Total number of atoms across all tokens
        - num_tokens: Number of tokens
        - num_atoms_max: Maximum atoms in any single token
        - num_atoms_min: Minimum atoms in any single token
    """
    num_atoms_total = sum([token["atom_num"] for token in tokenized.tokens])
    num_tokens = len(tokenized.tokens)
    num_atoms_max = max([token["atom_num"] for token in tokenized.tokens])
    num_atoms_min = min([token["atom_num"] for token in tokenized.tokens])

    return {
        "num_atoms_total": num_atoms_total,
        "num_tokens": num_tokens,
        "num_atoms_max": num_atoms_max,
        "num_atoms_min": num_atoms_min,
    }


def get_num_atoms_tokens(tokenized: Tokenized) -> tuple[int, int]:
    """Get the number of atoms and tokens from tokenized data.

    Parameters
    ----------
    tokenized : Tokenized
        The tokenized data.

    Returns
    -------
    tuple[int, int]
        Tuple of (num_atoms_total, num_tokens).
    """
    stats = tokenized_stats(tokenized)
    return stats["num_atoms_total"], stats["num_tokens"]
