import math
import random
from typing import Optional

import numba
import numpy as np
import numpy.typing as npt
import torch
from numba import types
from torch import Tensor, from_numpy
from torch.nn.functional import one_hot

from boltz.data import const
from boltz.data.feature.symmetry import (
    get_amino_acids_symmetries,
    get_chain_symmetries,
    get_ligand_symmetries,
)
from boltz.data.pad import pad_dim
from boltz.data.types import (
    MSA,
    MSADeletion,
    MSAResidue,
    MSASequence,
    Tokenized,
)
from boltz.model.modules.utils import center_random_augmentation

####################################################################################################
# HELPERS
####################################################################################################


def compute_frames_nonpolymer(
    data: Tokenized,
    coords,
    resolved_mask,
    atom_to_token,
    frame_data: list,
    resolved_frame_data: list,
) -> tuple[list, list]:
    """Get the frames for non-polymer tokens.

    Parameters
    ----------
    data : Tokenized
        The tokenized data.
    frame_data : list
        The frame data.
    resolved_frame_data : list
        The resolved frame data.

    Returns
    -------
    tuple[list, list]
        The frame data and resolved frame data.

    """
    frame_data = np.array(frame_data)
    resolved_frame_data = np.array(resolved_frame_data)
    asym_id_token = data.tokens["asym_id"]
    asym_id_atom = data.tokens["asym_id"][atom_to_token]
    token_idx = 0
    atom_idx = 0
    for id in np.unique(data.tokens["asym_id"]):
        mask_chain_token = asym_id_token == id
        mask_chain_atom = asym_id_atom == id
        num_tokens = mask_chain_token.sum()
        num_atoms = mask_chain_atom.sum()
        if (
            data.tokens[token_idx]["mol_type"] != const.chain_type_ids["NONPOLYMER"]
            or num_atoms < 3
        ):
            token_idx += num_tokens
            atom_idx += num_atoms
            continue
        dist_mat = (
            (
                coords.reshape(-1, 3)[mask_chain_atom][:, None, :]
                - coords.reshape(-1, 3)[mask_chain_atom][None, :, :]
            )
            ** 2
        ).sum(-1) ** 0.5
        resolved_pair = 1 - (
            resolved_mask[mask_chain_atom][None, :]
            * resolved_mask[mask_chain_atom][:, None]
        ).astype(np.float32)
        resolved_pair[resolved_pair == 1] = math.inf
        indices = np.argsort(dist_mat + resolved_pair, axis=1)
        frames = (
            np.concatenate(
                [
                    indices[:, 1:2],
                    indices[:, 0:1],
                    indices[:, 2:3],
                ],
                axis=1,
            )
            + atom_idx
        )
        frame_data[token_idx : token_idx + num_atoms, :] = frames
        resolved_frame_data[token_idx : token_idx + num_atoms] = resolved_mask[
            frames
        ].all(axis=1)
        token_idx += num_tokens
        atom_idx += num_atoms
    frames_expanded = coords.reshape(-1, 3)[frame_data]

    mask_collinear = compute_collinear_mask(
        frames_expanded[:, 1] - frames_expanded[:, 0],
        frames_expanded[:, 1] - frames_expanded[:, 2],
    )
    return frame_data, resolved_frame_data & mask_collinear


def compute_collinear_mask(v1, v2):
    norm1 = np.linalg.norm(v1, axis=1, keepdims=True)
    norm2 = np.linalg.norm(v2, axis=1, keepdims=True)
    v1 = v1 / (norm1 + 1e-6)
    v2 = v2 / (norm2 + 1e-6)
    mask_angle = np.abs(np.sum(v1 * v2, axis=1)) < 0.9063
    mask_overlap1 = norm1.reshape(-1) > 1e-2
    mask_overlap2 = norm2.reshape(-1) > 1e-2
    return mask_angle & mask_overlap1 & mask_overlap2


def dummy_msa(residues: np.ndarray) -> MSA:
    """Create a dummy MSA for a chain.

    Parameters
    ----------
    residues : np.ndarray
        The residues for the chain.

    Returns
    -------
    MSA
        The dummy MSA.

    """
    residues = [res["res_type"] for res in residues]
    deletions = []
    sequences = [(0, -1, 0, len(residues), 0, 0)]
    return MSA(
        residues=np.array(residues, dtype=MSAResidue),
        deletions=np.array(deletions, dtype=MSADeletion),
        sequences=np.array(sequences, dtype=MSASequence),
    )


def construct_paired_msa(  # noqa: C901, PLR0915, PLR0912
    data: Tokenized,
    max_seqs: int,
    max_pairs: int = 8192,
    max_total: int = 16384,
    random_subset: bool = False,
) -> tuple[Tensor, Tensor, Tensor]:
    """Pair the MSA data.

    Parameters
    ----------
    data : Input
        The input data.

    Returns
    -------
    Tensor
        The MSA data.
    Tensor
        The deletion data.
    Tensor
        Mask indicating paired sequences.

    """
    # Get unique chains (ensuring monotonicity in the order)
    assert np.all(np.diff(data.tokens["asym_id"], n=1) >= 0)
    chain_ids = np.unique(data.tokens["asym_id"])

    # Get relevant MSA, and create a dummy for chains without
    msa = {k: data.msa[k] for k in chain_ids if k in data.msa}
    for chain_id in chain_ids:
        if chain_id not in msa:
            chain = data.structure.chains[chain_id]
            res_start = chain["res_idx"]
            res_end = res_start + chain["res_num"]
            residues = data.structure.residues[res_start:res_end]
            msa[chain_id] = dummy_msa(residues)

    # Map taxonomies to (chain_id, seq_idx)
    taxonomy_map: dict[str, list] = {}
    for chain_id, chain_msa in msa.items():
        sequences = chain_msa.sequences
        sequences = sequences[sequences["taxonomy"] != -1]
        for sequence in sequences:
            seq_idx = sequence["seq_idx"]
            taxon = sequence["taxonomy"]
            taxonomy_map.setdefault(taxon, []).append((chain_id, seq_idx))

    # Remove taxonomies with only one sequence and sort by the
    # number of chain_id present in each of the taxonomies
    taxonomy_map = {k: v for k, v in taxonomy_map.items() if len(v) > 1}
    taxonomy_map = sorted(
        taxonomy_map.items(),
        key=lambda x: len({c for c, _ in x[1]}),
        reverse=True,
    )

    # Keep track of the sequences available per chain, keeping the original
    # order of the sequences in the MSA to favor the best matching sequences
    visited = {(c, s) for c, items in taxonomy_map for s in items}
    available = {}
    for c in chain_ids:
        available[c] = [
            i for i in range(1, len(msa[c].sequences)) if (c, i) not in visited
        ]

    # Create sequence pairs
    is_paired = []
    pairing = []

    # Start with the first sequence for each chain
    is_paired.append({c: 1 for c in chain_ids})
    pairing.append({c: 0 for c in chain_ids})

    # Then add up to 8191 paired rows
    for _, pairs in taxonomy_map:
        # Group occurences by chain_id in case we have multiple
        # sequences from the same chain and same taxonomy
        chain_occurences = {}
        for chain_id, seq_idx in pairs:
            chain_occurences.setdefault(chain_id, []).append(seq_idx)

        # We create as many pairings as the maximum number of occurences
        max_occurences = max(len(v) for v in chain_occurences.values())
        for i in range(max_occurences):
            row_pairing = {}
            row_is_paired = {}

            # Add the chains present in the taxonomy
            for chain_id, seq_idxs in chain_occurences.items():
                # Roll over the sequence index to maximize diversity
                idx = i % len(seq_idxs)
                seq_idx = seq_idxs[idx]

                # Add the sequence to the pairing
                row_pairing[chain_id] = seq_idx
                row_is_paired[chain_id] = 1

            # Add any missing chains
            for chain_id in chain_ids:
                if chain_id not in row_pairing:
                    row_is_paired[chain_id] = 0
                    if available[chain_id]:
                        # Add the next available sequence
                        seq_idx = available[chain_id].pop(0)
                        row_pairing[chain_id] = seq_idx
                    else:
                        # No more sequences available, we place a gap
                        row_pairing[chain_id] = -1

            pairing.append(row_pairing)
            is_paired.append(row_is_paired)

            # Break if we have enough pairs
            if len(pairing) >= max_pairs:
                break

        # Break if we have enough pairs
        if len(pairing) >= max_pairs:
            break

    # Now add up to 16384 unpaired rows total
    max_left = max(len(v) for v in available.values())
    for _ in range(min(max_total - len(pairing), max_left)):
        row_pairing = {}
        row_is_paired = {}
        for chain_id in chain_ids:
            row_is_paired[chain_id] = 0
            if available[chain_id]:
                # Add the next available sequence
                seq_idx = available[chain_id].pop(0)
                row_pairing[chain_id] = seq_idx
            else:
                # No more sequences available, we place a gap
                row_pairing[chain_id] = -1

        pairing.append(row_pairing)
        is_paired.append(row_is_paired)

        # Break if we have enough sequences
        if len(pairing) >= max_total:
            break

    # Randomly sample a subset of the pairs
    # ensuring the first row is always present
    if random_subset:
        num_seqs = len(pairing)
        if num_seqs > max_seqs:
            indices = np.random.choice(
                list(range(1, num_seqs)), size=max_seqs - 1, replace=False
            )  # noqa: NPY002
            pairing = [pairing[0]] + [pairing[i] for i in indices]
            is_paired = [is_paired[0]] + [is_paired[i] for i in indices]
    else:
        # Deterministic downsample to max_seqs
        pairing = pairing[:max_seqs]
        is_paired = is_paired[:max_seqs]

    # Map (chain_id, seq_idx, res_idx) to deletion
    deletions = {}
    for chain_id, chain_msa in msa.items():
        chain_deletions = chain_msa.deletions
        for sequence in chain_msa.sequences:
            del_start = sequence["del_start"]
            del_end = sequence["del_end"]
            chain_deletions = chain_msa.deletions[del_start:del_end]
            for deletion_data in chain_deletions:
                seq_idx = sequence["seq_idx"]
                res_idx = deletion_data["res_idx"]
                deletion = deletion_data["deletion"]
                deletions[(chain_id, seq_idx, res_idx)] = deletion

    # Add all the token MSA data
    msa_data, del_data, paired_data = prepare_msa_arrays(
        data.tokens, pairing, is_paired, deletions, msa
    )

    msa_data = torch.tensor(msa_data, dtype=torch.long)
    del_data = torch.tensor(del_data, dtype=torch.float)
    paired_data = torch.tensor(paired_data, dtype=torch.float)

    return msa_data, del_data, paired_data


def prepare_msa_arrays(
    tokens,
    pairing: list[dict[int, int]],
    is_paired: list[dict[int, int]],
    deletions: dict[tuple[int, int, int], int],
    msa: dict[int, MSA],
) -> tuple[npt.NDArray[np.int64], npt.NDArray[np.int64], npt.NDArray[np.int64]]:
    """Reshape data to play nicely with numba jit."""
    token_asym_ids_arr = np.array([t["asym_id"] for t in tokens], dtype=np.int64)
    token_res_idxs_arr = np.array([t["res_idx"] for t in tokens], dtype=np.int64)

    chain_ids = sorted(msa.keys())

    # chain_ids are not necessarily contiguous (e.g. they might be 0, 24, 25).
    # This allows us to look up a chain_id by it's index in the chain_ids list.
    chain_id_to_idx = {chain_id: i for i, chain_id in enumerate(chain_ids)}
    token_asym_ids_idx_arr = np.array(
        [chain_id_to_idx[asym_id] for asym_id in token_asym_ids_arr], dtype=np.int64
    )

    pairing_arr = np.zeros((len(pairing), len(chain_ids)), dtype=np.int64)
    is_paired_arr = np.zeros((len(is_paired), len(chain_ids)), dtype=np.int64)

    for i, row_pairing in enumerate(pairing):
        for chain_id in chain_ids:
            pairing_arr[i, chain_id_to_idx[chain_id]] = row_pairing[chain_id]

    for i, row_is_paired in enumerate(is_paired):
        for chain_id in chain_ids:
            is_paired_arr[i, chain_id_to_idx[chain_id]] = row_is_paired[chain_id]

    max_seq_len = max(len(msa[chain_id].sequences) for chain_id in chain_ids)

    # we want res_start from sequences
    msa_sequences = np.full((len(chain_ids), max_seq_len), -1, dtype=np.int64)
    for chain_id in chain_ids:
        for i, seq in enumerate(msa[chain_id].sequences):
            msa_sequences[chain_id_to_idx[chain_id], i] = seq["res_start"]

    max_residues_len = max(len(msa[chain_id].residues) for chain_id in chain_ids)
    msa_residues = np.full((len(chain_ids), max_residues_len), -1, dtype=np.int64)
    for chain_id in chain_ids:
        residues = msa[chain_id].residues.astype(np.int64)
        idxs = np.arange(len(residues))
        chain_idx = chain_id_to_idx[chain_id]
        msa_residues[chain_idx, idxs] = residues

    deletions_dict = numba.typed.Dict.empty(
        key_type=numba.types.Tuple(
            [numba.types.int64, numba.types.int64, numba.types.int64]
        ),
        value_type=numba.types.int64,
    )
    deletions_dict.update(deletions)

    return _prepare_msa_arrays_inner(
        token_asym_ids_arr,
        token_res_idxs_arr,
        token_asym_ids_idx_arr,
        pairing_arr,
        is_paired_arr,
        deletions_dict,
        msa_sequences,
        msa_residues,
        const.token_ids["-"],
    )


deletions_dict_type = types.DictType(types.UniTuple(types.int64, 3), types.int64)


@numba.njit(
    [
        types.Tuple(
            (
                types.int64[:, ::1],  # msa_data
                types.int64[:, ::1],  # del_data
                types.int64[:, ::1],  # paired_data
            )
        )(
            types.int64[::1],  # token_asym_ids
            types.int64[::1],  # token_res_idxs
            types.int64[::1],  # token_asym_ids_idx
            types.int64[:, ::1],  # pairing
            types.int64[:, ::1],  # is_paired
            deletions_dict_type,  # deletions
            types.int64[:, ::1],  # msa_sequences
            types.int64[:, ::1],  # msa_residues
            types.int64,  # gap_token
        )
    ],
    cache=True,
)
def _prepare_msa_arrays_inner(
    token_asym_ids: npt.NDArray[np.int64],
    token_res_idxs: npt.NDArray[np.int64],
    token_asym_ids_idx: npt.NDArray[np.int64],
    pairing: npt.NDArray[np.int64],
    is_paired: npt.NDArray[np.int64],
    deletions: dict[tuple[int, int, int], int],
    msa_sequences: npt.NDArray[np.int64],
    msa_residues: npt.NDArray[np.int64],
    gap_token: int,
) -> tuple[npt.NDArray[np.int64], npt.NDArray[np.int64], npt.NDArray[np.int64]]:
    n_tokens = len(token_asym_ids)
    n_pairs = len(pairing)
    msa_data = np.full((n_tokens, n_pairs), gap_token, dtype=np.int64)
    paired_data = np.zeros((n_tokens, n_pairs), dtype=np.int64)
    del_data = np.zeros((n_tokens, n_pairs), dtype=np.int64)

    # Add all the token MSA data
    for token_idx in range(n_tokens):
        chain_id_idx = token_asym_ids_idx[token_idx]
        chain_id = token_asym_ids[token_idx]
        res_idx = token_res_idxs[token_idx]

        for pair_idx in range(n_pairs):
            seq_idx = pairing[pair_idx, chain_id_idx]
            paired_data[token_idx, pair_idx] = is_paired[pair_idx, chain_id_idx]

            # Add residue type
            if seq_idx != -1:
                res_start = msa_sequences[chain_id_idx, seq_idx]
                res_type = msa_residues[chain_id_idx, res_start + res_idx]
                k = (chain_id, seq_idx, res_idx)
                if k in deletions:
                    del_data[token_idx, pair_idx] = deletions[k]
                msa_data[token_idx, pair_idx] = res_type

    return msa_data, del_data, paired_data


####################################################################################################
# FEATURES
####################################################################################################


def select_subset_from_mask(mask, p):
    num_true = np.sum(mask)
    v = np.random.geometric(p) + 1
    k = min(v, num_true)

    true_indices = np.where(mask)[0]

    # Randomly select k indices from the true_indices
    selected_indices = np.random.choice(true_indices, size=k, replace=False)

    new_mask = np.zeros_like(mask)
    new_mask[selected_indices] = 1

    return new_mask


def process_token_features(
    data: Tokenized,
    max_tokens: Optional[int] = None,
    binder_pocket_conditioned_prop: Optional[float] = 0.0,
    binder_pocket_cutoff: Optional[float] = 6.0,
    binder_pocket_sampling_geometric_p: Optional[float] = 0.0,
    only_ligand_binder_pocket: Optional[bool] = False,
    inference_binder: Optional[list[int]] = None,
    inference_pocket: Optional[list[tuple[int, int]]] = None,
) -> dict[str, Tensor]:
    """Get the token features.

    Parameters
    ----------
    data : Tokenized
        The tokenized data.
    max_tokens : int
        The maximum number of tokens.

    Returns
    -------
    dict[str, Tensor]
        The token features.

    """
    # Token data
    token_data = data.tokens
    token_bonds = data.bonds

    # Token core features
    token_index = torch.arange(len(token_data), dtype=torch.long)
    residue_index = from_numpy(token_data["res_idx"].copy()).long()
    asym_id = from_numpy(token_data["asym_id"].copy()).long()
    entity_id = from_numpy(token_data["entity_id"].copy()).long()
    sym_id = from_numpy(token_data["sym_id"].copy()).long()
    mol_type = from_numpy(token_data["mol_type"].copy()).long()
    res_type = from_numpy(token_data["res_type"].copy()).long()
    res_type = one_hot(res_type, num_classes=const.num_tokens)
    disto_center = from_numpy(token_data["disto_coords"].copy())

    # Token mask features
    pad_mask = torch.ones(len(token_data), dtype=torch.float)
    resolved_mask = from_numpy(token_data["resolved_mask"].copy()).float()
    disto_mask = from_numpy(token_data["disto_mask"].copy()).float()
    cyclic_period = from_numpy(token_data["cyclic_period"].copy())

    # Token bond features
    if max_tokens is not None:
        pad_len = max_tokens - len(token_data)
        num_tokens = max_tokens if pad_len > 0 else len(token_data)
    else:
        num_tokens = len(token_data)

    tok_to_idx = {tok["token_idx"]: idx for idx, tok in enumerate(token_data)}
    bonds = torch.zeros(num_tokens, num_tokens, dtype=torch.float)
    for token_bond in token_bonds:
        token_1 = tok_to_idx[token_bond["token_1"]]
        token_2 = tok_to_idx[token_bond["token_2"]]
        bonds[token_1, token_2] = 1
        bonds[token_2, token_1] = 1

    bonds = bonds.unsqueeze(-1)

    # Pocket conditioned feature
    pocket_feature = (
        np.zeros(len(token_data)) + const.pocket_contact_info["UNSPECIFIED"]
    )
    if inference_binder is not None:
        assert inference_pocket is not None
        pocket_residues = set(inference_pocket)
        for idx, token in enumerate(token_data):
            if token["asym_id"] in inference_binder:
                pocket_feature[idx] = const.pocket_contact_info["BINDER"]
            elif (token["asym_id"], token["res_idx"]) in pocket_residues:
                pocket_feature[idx] = const.pocket_contact_info["POCKET"]
            else:
                pocket_feature[idx] = const.pocket_contact_info["UNSELECTED"]
    elif (
        binder_pocket_conditioned_prop > 0.0
        and random.random() < binder_pocket_conditioned_prop
    ):
        # choose as binder a random ligand in the crop, if there are no ligands select a protein chain
        binder_asym_ids = np.unique(
            token_data["asym_id"][
                token_data["mol_type"] == const.chain_type_ids["NONPOLYMER"]
            ]
        )

        if len(binder_asym_ids) == 0:
            if not only_ligand_binder_pocket:
                binder_asym_ids = np.unique(token_data["asym_id"])

        if len(binder_asym_ids) > 0:
            pocket_asym_id = random.choice(binder_asym_ids)
            binder_mask = token_data["asym_id"] == pocket_asym_id

            binder_coords = []
            for token in token_data:
                if token["asym_id"] == pocket_asym_id:
                    binder_coords.append(
                        data.structure.atoms["coords"][
                            token["atom_idx"] : token["atom_idx"] + token["atom_num"]
                        ]
                    )
            binder_coords = np.concatenate(binder_coords, axis=0)

            # find the tokens in the pocket
            token_dist = np.zeros(len(token_data)) + 1000
            for i, token in enumerate(token_data):
                if (
                    token["mol_type"] != const.chain_type_ids["NONPOLYMER"]
                    and token["asym_id"] != pocket_asym_id
                    and token["resolved_mask"] == 1
                ):
                    token_coords = data.structure.atoms["coords"][
                        token["atom_idx"] : token["atom_idx"] + token["atom_num"]
                    ]

                    # find chain and apply chain transformation
                    for chain in data.structure.chains:
                        if chain["asym_id"] == token["asym_id"]:
                            break

                    token_dist[i] = np.min(
                        np.linalg.norm(
                            token_coords[:, None, :] - binder_coords[None, :, :],
                            axis=-1,
                        )
                    )

            pocket_mask = token_dist < binder_pocket_cutoff

            if np.sum(pocket_mask) > 0:
                pocket_feature = (
                    np.zeros(len(token_data)) + const.pocket_contact_info["UNSELECTED"]
                )
                pocket_feature[binder_mask] = const.pocket_contact_info["BINDER"]

                if binder_pocket_sampling_geometric_p > 0.0:
                    # select a subset of the pocket, according
                    # to a geometric distribution with one as minimum
                    pocket_mask = select_subset_from_mask(
                        pocket_mask, binder_pocket_sampling_geometric_p
                    )

                pocket_feature[pocket_mask] = const.pocket_contact_info["POCKET"]
    pocket_feature = from_numpy(pocket_feature).long()
    pocket_feature = one_hot(pocket_feature, num_classes=len(const.pocket_contact_info))

    # Pad to max tokens if given
    if max_tokens is not None:
        pad_len = max_tokens - len(token_data)
        if pad_len > 0:
            token_index = pad_dim(token_index, 0, pad_len)
            residue_index = pad_dim(residue_index, 0, pad_len)
            asym_id = pad_dim(asym_id, 0, pad_len)
            entity_id = pad_dim(entity_id, 0, pad_len)
            sym_id = pad_dim(sym_id, 0, pad_len)
            mol_type = pad_dim(mol_type, 0, pad_len)
            res_type = pad_dim(res_type, 0, pad_len)
            disto_center = pad_dim(disto_center, 0, pad_len)
            pad_mask = pad_dim(pad_mask, 0, pad_len)
            resolved_mask = pad_dim(resolved_mask, 0, pad_len)
            disto_mask = pad_dim(disto_mask, 0, pad_len)
            pocket_feature = pad_dim(pocket_feature, 0, pad_len)

    token_features = {
        "token_index": token_index,
        "residue_index": residue_index,
        "asym_id": asym_id,
        "entity_id": entity_id,
        "sym_id": sym_id,
        "mol_type": mol_type,
        "res_type": res_type,
        "disto_center": disto_center,
        "token_bonds": bonds,
        "token_pad_mask": pad_mask,
        "token_resolved_mask": resolved_mask,
        "token_disto_mask": disto_mask,
        "pocket_feature": pocket_feature,
        "cyclic_period": cyclic_period,
    }
    return token_features


def process_atom_features(
    data: Tokenized,
    atoms_per_window_queries: int = 32,
    min_dist: float = 2.0,
    max_dist: float = 22.0,
    num_bins: int = 64,
    max_atoms: Optional[int] = None,
    max_tokens: Optional[int] = None,
) -> dict[str, Tensor]:
    """Get the atom features.

    Parameters
    ----------
    data : Tokenized
        The tokenized data.
    max_atoms : int, optional
        The maximum number of atoms.

    Returns
    -------
    dict[str, Tensor]
        The atom features.

    """
    # Filter to tokens' atoms
    atom_data = []
    ref_space_uid = []
    coord_data = []
    frame_data = []
    resolved_frame_data = []
    atom_to_token = []
    token_to_rep_atom = []  # index on cropped atom table
    r_set_to_rep_atom = []
    disto_coords = []
    atom_idx = 0

    chain_res_ids = {}
    for token_id, token in enumerate(data.tokens):
        # Get the chain residue ids
        chain_idx, res_id = token["asym_id"], token["res_idx"]
        chain = data.structure.chains[chain_idx]

        if (chain_idx, res_id) not in chain_res_ids:
            new_idx = len(chain_res_ids)
            chain_res_ids[(chain_idx, res_id)] = new_idx
        else:
            new_idx = chain_res_ids[(chain_idx, res_id)]

        # Map atoms to token indices
        ref_space_uid.extend([new_idx] * token["atom_num"])
        atom_to_token.extend([token_id] * token["atom_num"])

        # Add atom data
        start = token["atom_idx"]
        end = token["atom_idx"] + token["atom_num"]
        token_atoms = data.structure.atoms[start:end]

        # Map token to representative atom
        token_to_rep_atom.append(atom_idx + token["disto_idx"] - start)
        if (chain["mol_type"] != const.chain_type_ids["NONPOLYMER"]) and token[
            "resolved_mask"
        ]:
            r_set_to_rep_atom.append(atom_idx + token["center_idx"] - start)

        # Get token coordinates
        token_coords = np.array([token_atoms["coords"]])
        coord_data.append(token_coords)

        # Get frame data
        res_type = const.tokens[token["res_type"]]

        if token["atom_num"] < 3 or res_type in ["PAD", "UNK", "-"]:
            idx_frame_a, idx_frame_b, idx_frame_c = 0, 0, 0
            mask_frame = False
        elif (token["mol_type"] == const.chain_type_ids["PROTEIN"]) and (
            res_type in const.ref_atoms
        ):
            idx_frame_a, idx_frame_b, idx_frame_c = (
                const.ref_atoms[res_type].index("N"),
                const.ref_atoms[res_type].index("CA"),
                const.ref_atoms[res_type].index("C"),
            )
            mask_frame = (
                token_atoms["is_present"][idx_frame_a]
                and token_atoms["is_present"][idx_frame_b]
                and token_atoms["is_present"][idx_frame_c]
            )
        elif (
            token["mol_type"] == const.chain_type_ids["DNA"]
            or token["mol_type"] == const.chain_type_ids["RNA"]
        ) and (res_type in const.ref_atoms):
            idx_frame_a, idx_frame_b, idx_frame_c = (
                const.ref_atoms[res_type].index("C1'"),
                const.ref_atoms[res_type].index("C3'"),
                const.ref_atoms[res_type].index("C4'"),
            )
            mask_frame = (
                token_atoms["is_present"][idx_frame_a]
                and token_atoms["is_present"][idx_frame_b]
                and token_atoms["is_present"][idx_frame_c]
            )
        else:
            idx_frame_a, idx_frame_b, idx_frame_c = 0, 0, 0
            mask_frame = False
        frame_data.append(
            [idx_frame_a + atom_idx, idx_frame_b + atom_idx, idx_frame_c + atom_idx]
        )
        resolved_frame_data.append(mask_frame)

        # Get distogram coordinates
        disto_coords_tok = data.structure.atoms[token["disto_idx"]]["coords"]
        disto_coords.append(disto_coords_tok)

        # Update atom data. This is technically never used again (we rely on coord_data),
        # but we update for consistency and to make sure the Atom object has valid, transformed coordinates.
        token_atoms = token_atoms.copy()
        token_atoms["coords"] = token_coords[0]  # atom has a copy of first coords
        atom_data.append(token_atoms)
        atom_idx += len(token_atoms)

    disto_coords = np.array(disto_coords)

    # Compute distogram
    t_center = torch.Tensor(disto_coords)
    t_dists = torch.cdist(t_center, t_center)
    boundaries = torch.linspace(min_dist, max_dist, num_bins - 1)
    distogram = (t_dists.unsqueeze(-1) > boundaries).sum(dim=-1).long()
    disto_target = one_hot(distogram, num_classes=num_bins)

    atom_data = np.concatenate(atom_data)
    coord_data = np.concatenate(coord_data, axis=1)
    ref_space_uid = np.array(ref_space_uid)

    # Compute features
    ref_atom_name_chars = from_numpy(atom_data["name"]).long()
    ref_element = from_numpy(atom_data["element"]).long()
    ref_charge = from_numpy(atom_data["charge"])
    ref_pos = from_numpy(
        atom_data["conformer"].copy()
    )  # not sure why I need to copy here..
    ref_space_uid = from_numpy(ref_space_uid)
    coords = from_numpy(coord_data.copy())
    resolved_mask = from_numpy(atom_data["is_present"])
    pad_mask = torch.ones(len(atom_data), dtype=torch.float)
    atom_to_token = torch.tensor(atom_to_token, dtype=torch.long)
    token_to_rep_atom = torch.tensor(token_to_rep_atom, dtype=torch.long)
    r_set_to_rep_atom = torch.tensor(r_set_to_rep_atom, dtype=torch.long)
    frame_data, resolved_frame_data = compute_frames_nonpolymer(
        data,
        coord_data,
        atom_data["is_present"],
        atom_to_token,
        frame_data,
        resolved_frame_data,
    )  # Compute frames for NONPOLYMER tokens
    frames = from_numpy(frame_data.copy())
    frame_resolved_mask = from_numpy(resolved_frame_data.copy())
    # Convert to one-hot
    ref_atom_name_chars = one_hot(
        ref_atom_name_chars % num_bins, num_classes=num_bins
    )  # added for lower case letters
    ref_element = one_hot(ref_element, num_classes=const.num_elements)
    atom_to_token = one_hot(atom_to_token, num_classes=token_id + 1)
    token_to_rep_atom = one_hot(token_to_rep_atom, num_classes=len(atom_data))
    r_set_to_rep_atom = one_hot(r_set_to_rep_atom, num_classes=len(atom_data))

    # Center the ground truth coordinates
    center = (coords * resolved_mask[None, :, None]).sum(dim=1)
    center = center / resolved_mask.sum().clamp(min=1)
    coords = coords - center[:, None]

    # Apply random roto-translation to the input atoms
    ref_pos = center_random_augmentation(
        ref_pos[None], resolved_mask[None], centering=False
    )[0]

    # Compute padding and apply
    if max_atoms is not None:
        assert max_atoms % atoms_per_window_queries == 0
        pad_len = max_atoms - len(atom_data)
    else:
        pad_len = (
            (len(atom_data) - 1) // atoms_per_window_queries + 1
        ) * atoms_per_window_queries - len(atom_data)

    if pad_len > 0:
        pad_mask = pad_dim(pad_mask, 0, pad_len)
        ref_pos = pad_dim(ref_pos, 0, pad_len)
        resolved_mask = pad_dim(resolved_mask, 0, pad_len)
        ref_element = pad_dim(ref_element, 0, pad_len)
        ref_charge = pad_dim(ref_charge, 0, pad_len)
        ref_atom_name_chars = pad_dim(ref_atom_name_chars, 0, pad_len)
        ref_space_uid = pad_dim(ref_space_uid, 0, pad_len)
        coords = pad_dim(coords, 1, pad_len)
        atom_to_token = pad_dim(atom_to_token, 0, pad_len)
        token_to_rep_atom = pad_dim(token_to_rep_atom, 1, pad_len)
        r_set_to_rep_atom = pad_dim(r_set_to_rep_atom, 1, pad_len)

    if max_tokens is not None:
        pad_len = max_tokens - token_to_rep_atom.shape[0]
        if pad_len > 0:
            atom_to_token = pad_dim(atom_to_token, 1, pad_len)
            token_to_rep_atom = pad_dim(token_to_rep_atom, 0, pad_len)
            r_set_to_rep_atom = pad_dim(r_set_to_rep_atom, 0, pad_len)
            disto_target = pad_dim(pad_dim(disto_target, 0, pad_len), 1, pad_len)
            frames = pad_dim(frames, 0, pad_len)
            frame_resolved_mask = pad_dim(frame_resolved_mask, 0, pad_len)

    return {
        "ref_pos": ref_pos,
        "atom_resolved_mask": resolved_mask,
        "ref_element": ref_element,
        "ref_charge": ref_charge,
        "ref_atom_name_chars": ref_atom_name_chars,
        "ref_space_uid": ref_space_uid,
        "coords": coords,
        "atom_pad_mask": pad_mask,
        "atom_to_token": atom_to_token,
        "token_to_rep_atom": token_to_rep_atom,
        "r_set_to_rep_atom": r_set_to_rep_atom,
        "disto_target": disto_target,
        "frames_idx": frames,
        "frame_resolved_mask": frame_resolved_mask,
    }


def process_msa_features(
    data: Tokenized,
    max_seqs_batch: int,
    max_seqs: int,
    max_tokens: Optional[int] = None,
    pad_to_max_seqs: bool = False,
) -> dict[str, Tensor]:
    """Get the MSA features.

    Parameters
    ----------
    data : Tokenized
        The tokenized data.
    max_seqs : int
        The maximum number of MSA sequences.
    max_tokens : int
        The maximum number of tokens.
    pad_to_max_seqs : bool
        Whether to pad to the maximum number of sequences.

    Returns
    -------
    dict[str, Tensor]
        The MSA features.

    """
    # Created paired MSA
    msa, deletion, paired = construct_paired_msa(data, max_seqs_batch)
    msa, deletion, paired = (
        msa.transpose(1, 0),
        deletion.transpose(1, 0),
        paired.transpose(1, 0),
    )  # (N_MSA, N_RES, N_AA)

    # Prepare features
    msa = torch.nn.functional.one_hot(msa, num_classes=const.num_tokens)
    msa_mask = torch.ones_like(msa[:, :, 0])
    profile = msa.float().mean(dim=0)
    has_deletion = deletion > 0
    deletion = np.pi / 2 * np.arctan(deletion / 3)
    deletion_mean = deletion.mean(axis=0)

    # Pad in the MSA dimension (dim=0)
    if pad_to_max_seqs:
        pad_len = max_seqs - msa.shape[0]
        if pad_len > 0:
            msa = pad_dim(msa, 0, pad_len, const.token_ids["-"])
            paired = pad_dim(paired, 0, pad_len)
            msa_mask = pad_dim(msa_mask, 0, pad_len)
            has_deletion = pad_dim(has_deletion, 0, pad_len)
            deletion = pad_dim(deletion, 0, pad_len)

    # Pad in the token dimension (dim=1)
    if max_tokens is not None:
        pad_len = max_tokens - msa.shape[1]
        if pad_len > 0:
            msa = pad_dim(msa, 1, pad_len, const.token_ids["-"])
            paired = pad_dim(paired, 1, pad_len)
            msa_mask = pad_dim(msa_mask, 1, pad_len)
            has_deletion = pad_dim(has_deletion, 1, pad_len)
            deletion = pad_dim(deletion, 1, pad_len)
            profile = pad_dim(profile, 0, pad_len)
            deletion_mean = pad_dim(deletion_mean, 0, pad_len)

    return {
        "msa": msa,
        "msa_paired": paired,
        "deletion_value": deletion,
        "has_deletion": has_deletion,
        "deletion_mean": deletion_mean,
        "profile": profile,
        "msa_mask": msa_mask,
    }


def process_symmetry_features(
    cropped: Tokenized, symmetries: dict
) -> dict[str, Tensor]:
    """Get the symmetry features.

    Parameters
    ----------
    data : Tokenized
        The tokenized data.

    Returns
    -------
    dict[str, Tensor]
        The symmetry features.

    """
    features = get_chain_symmetries(cropped)
    features.update(get_amino_acids_symmetries(cropped))
    features.update(get_ligand_symmetries(cropped, symmetries))

    return features


def process_residue_constraint_features(
    data: Tokenized,
) -> dict[str, Tensor]:
    residue_constraints = data.residue_constraints
    if residue_constraints is not None:
        rdkit_bounds_constraints = residue_constraints.rdkit_bounds_constraints
        chiral_atom_constraints = residue_constraints.chiral_atom_constraints
        stereo_bond_constraints = residue_constraints.stereo_bond_constraints
        planar_bond_constraints = residue_constraints.planar_bond_constraints
        planar_ring_5_constraints = residue_constraints.planar_ring_5_constraints
        planar_ring_6_constraints = residue_constraints.planar_ring_6_constraints

        rdkit_bounds_index = torch.tensor(
            rdkit_bounds_constraints["atom_idxs"].copy(), dtype=torch.long
        ).T
        rdkit_bounds_bond_mask = torch.tensor(
            rdkit_bounds_constraints["is_bond"].copy(), dtype=torch.bool
        )
        rdkit_bounds_angle_mask = torch.tensor(
            rdkit_bounds_constraints["is_angle"].copy(), dtype=torch.bool
        )
        rdkit_upper_bounds = torch.tensor(
            rdkit_bounds_constraints["upper_bound"].copy(), dtype=torch.float
        )
        rdkit_lower_bounds = torch.tensor(
            rdkit_bounds_constraints["lower_bound"].copy(), dtype=torch.float
        )

        chiral_atom_index = torch.tensor(
            chiral_atom_constraints["atom_idxs"].copy(), dtype=torch.long
        ).T
        chiral_reference_mask = torch.tensor(
            chiral_atom_constraints["is_reference"].copy(), dtype=torch.bool
        )
        chiral_atom_orientations = torch.tensor(
            chiral_atom_constraints["is_r"].copy(), dtype=torch.bool
        )

        stereo_bond_index = torch.tensor(
            stereo_bond_constraints["atom_idxs"].copy(), dtype=torch.long
        ).T
        stereo_reference_mask = torch.tensor(
            stereo_bond_constraints["is_reference"].copy(), dtype=torch.bool
        )
        stereo_bond_orientations = torch.tensor(
            stereo_bond_constraints["is_e"].copy(), dtype=torch.bool
        )

        planar_bond_index = torch.tensor(
            planar_bond_constraints["atom_idxs"].copy(), dtype=torch.long
        ).T
        planar_ring_5_index = torch.tensor(
            planar_ring_5_constraints["atom_idxs"].copy(), dtype=torch.long
        ).T
        planar_ring_6_index = torch.tensor(
            planar_ring_6_constraints["atom_idxs"].copy(), dtype=torch.long
        ).T
    else:
        rdkit_bounds_index = torch.empty((2, 0), dtype=torch.long)
        rdkit_bounds_bond_mask = torch.empty((0,), dtype=torch.bool)
        rdkit_bounds_angle_mask = torch.empty((0,), dtype=torch.bool)
        rdkit_upper_bounds = torch.empty((0,), dtype=torch.float)
        rdkit_lower_bounds = torch.empty((0,), dtype=torch.float)
        chiral_atom_index = torch.empty(
            (
                4,
                0,
            ),
            dtype=torch.long,
        )
        chiral_reference_mask = torch.empty((0,), dtype=torch.bool)
        chiral_atom_orientations = torch.empty((0,), dtype=torch.bool)
        stereo_bond_index = torch.empty((4, 0), dtype=torch.long)
        stereo_reference_mask = torch.empty((0,), dtype=torch.bool)
        stereo_bond_orientations = torch.empty((0,), dtype=torch.bool)
        planar_bond_index = torch.empty((6, 0), dtype=torch.long)
        planar_ring_5_index = torch.empty((5, 0), dtype=torch.long)
        planar_ring_6_index = torch.empty((6, 0), dtype=torch.long)

    return {
        "rdkit_bounds_index": rdkit_bounds_index,
        "rdkit_bounds_bond_mask": rdkit_bounds_bond_mask,
        "rdkit_bounds_angle_mask": rdkit_bounds_angle_mask,
        "rdkit_upper_bounds": rdkit_upper_bounds,
        "rdkit_lower_bounds": rdkit_lower_bounds,
        "chiral_atom_index": chiral_atom_index,
        "chiral_reference_mask": chiral_reference_mask,
        "chiral_atom_orientations": chiral_atom_orientations,
        "stereo_bond_index": stereo_bond_index,
        "stereo_reference_mask": stereo_reference_mask,
        "stereo_bond_orientations": stereo_bond_orientations,
        "planar_bond_index": planar_bond_index,
        "planar_ring_5_index": planar_ring_5_index,
        "planar_ring_6_index": planar_ring_6_index,
    }


def process_chain_feature_constraints(
    data: Tokenized,
) -> dict[str, Tensor]:
    structure = data.structure
    if structure.connections.shape[0] > 0:
        connected_chain_index, connected_atom_index = [], []
        for connection in structure.connections:
            connected_chain_index.append([connection["chain_1"], connection["chain_2"]])
            connected_atom_index.append([connection["atom_1"], connection["atom_2"]])
        connected_chain_index = torch.tensor(connected_chain_index, dtype=torch.long).T
        connected_atom_index = torch.tensor(connected_atom_index, dtype=torch.long).T
    else:
        connected_chain_index = torch.empty((2, 0), dtype=torch.long)
        connected_atom_index = torch.empty((2, 0), dtype=torch.long)

    symmetric_chain_index = []
    for i, chain_i in enumerate(structure.chains):
        for j, chain_j in enumerate(structure.chains):
            if j <= i:
                continue
            if chain_i["entity_id"] == chain_j["entity_id"]:
                symmetric_chain_index.append([i, j])
    if len(symmetric_chain_index) > 0:
        symmetric_chain_index = torch.tensor(symmetric_chain_index, dtype=torch.long).T
    else:
        symmetric_chain_index = torch.empty((2, 0), dtype=torch.long)
    return {
        "connected_chain_index": connected_chain_index,
        "connected_atom_index": connected_atom_index,
        "symmetric_chain_index": symmetric_chain_index,
    }


class BoltzFeaturizer:
    """Boltz featurizer."""

    def process(
        self,
        data: Tokenized,
        training: bool,
        max_seqs: int = 4096,
        atoms_per_window_queries: int = 32,
        min_dist: float = 2.0,
        max_dist: float = 22.0,
        num_bins: int = 64,
        max_tokens: Optional[int] = None,
        max_atoms: Optional[int] = None,
        pad_to_max_seqs: bool = False,
        compute_symmetries: bool = False,
        symmetries: Optional[dict] = None,
        binder_pocket_conditioned_prop: Optional[float] = 0.0,
        binder_pocket_cutoff: Optional[float] = 6.0,
        binder_pocket_sampling_geometric_p: Optional[float] = 0.0,
        only_ligand_binder_pocket: Optional[bool] = False,
        inference_binder: Optional[int] = None,
        inference_pocket: Optional[list[tuple[int, int]]] = None,
        compute_constraint_features: bool = False,
    ) -> dict[str, Tensor]:
        """Compute features.

        Parameters
        ----------
        data : Tokenized
            The tokenized data.
        training : bool
            Whether the model is in training mode.
        max_tokens : int, optional
            The maximum number of tokens.
        max_atoms : int, optional
            The maximum number of atoms
        max_seqs : int, optional
            The maximum number of sequences.

        Returns
        -------
        dict[str, Tensor]
            The features for model training.

        """
        # Compute random number of sequences
        if training and max_seqs is not None:
            max_seqs_batch = np.random.randint(1, max_seqs + 1)  # noqa: NPY002
        else:
            max_seqs_batch = max_seqs

        # Compute token features
        token_features = process_token_features(
            data,
            max_tokens,
            binder_pocket_conditioned_prop,
            binder_pocket_cutoff,
            binder_pocket_sampling_geometric_p,
            only_ligand_binder_pocket,
            inference_binder=inference_binder,
            inference_pocket=inference_pocket,
        )

        # Compute atom features
        atom_features = process_atom_features(
            data,
            atoms_per_window_queries,
            min_dist,
            max_dist,
            num_bins,
            max_atoms,
            max_tokens,
        )

        # Compute MSA features
        msa_features = process_msa_features(
            data,
            max_seqs_batch,
            max_seqs,
            max_tokens,
            pad_to_max_seqs,
        )

        # Compute symmetry features
        symmetry_features = {}
        if compute_symmetries:
            symmetry_features = process_symmetry_features(data, symmetries)

        # Compute residue constraint features
        residue_constraint_features = {}
        chain_constraint_features = {}
        if compute_constraint_features:
            residue_constraint_features = process_residue_constraint_features(data)
            chain_constraint_features = process_chain_feature_constraints(data)

        return {
            **token_features,
            **atom_features,
            **msa_features,
            **symmetry_features,
            **residue_constraint_features,
            **chain_constraint_features,
        }
