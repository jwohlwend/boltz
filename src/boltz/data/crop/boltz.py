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

# fmt: off

from dataclasses import replace
from typing import Optional

import numpy as np
from scipy.spatial.distance import cdist

from boltz.data import const
from boltz.data.crop.cropper import Cropper
from boltz.data.types import TokenizedTraining


def pick_random_token(
    tokens: np.ndarray,
    random: np.random.Generator,
) -> np.ndarray:
    """Pick a random token from the data.

    Parameters
    ----------
    tokens : np.ndarray
        The token data.
    random : np.random.Generator
        The random state for reproducibility.

    Returns
    -------
    np.ndarray
        The selected token.

    """
    return tokens[random.randint(len(tokens))]


def pick_chain_token(
    tokens: np.ndarray,
    chain_id: int,
    random: np.random.Generator,
) -> np.ndarray:
    """Pick a random token from a chain.

    Parameters
    ----------
    tokens : np.ndarray
        The token data.
    chain_id : int
        The chain ID.
    random : np.ndarray
        The random state for reproducibility.

    Returns
    -------
    np.ndarray
        The selected token.

    """
    # Filter to chain
    chain_tokens = tokens[tokens["asym_id"] == chain_id]

    # Pick from chain, fallback to all tokens
    if chain_tokens.size:
        query = pick_random_token(chain_tokens, random)
    else:
        query = pick_random_token(tokens, random)

    return query


def pick_interface_token(
    tokens: np.ndarray,
    interface: np.ndarray,
    random: np.random.Generator,
) -> np.ndarray:
    """Pick a random token from an interface.

    Parameters
    ----------
    tokens : np.ndarray
        The token data.
    interface : int
        The interface ID.
    random : np.random.Generator
        The random state for reproducibility.

    Returns
    -------
    np.ndarray
        The selected token.

    """
    # Sample random interface
    chain_1 = int(interface["chain_1"])
    chain_2 = int(interface["chain_2"])

    tokens_1 = tokens[tokens["asym_id"] == chain_1]
    tokens_2 = tokens[tokens["asym_id"] == chain_2]

    # If no interface, pick from the chains
    if tokens_1.size and (not tokens_2.size):
        query = pick_random_token(tokens_1, random)
    elif tokens_2.size and (not tokens_1.size):
        query = pick_random_token(tokens_2, random)
    elif (not tokens_1.size) and (not tokens_2.size):
        query = pick_random_token(tokens, random)
    else:
        # If we have tokens, compute distances
        tokens_1_coords = tokens_1["center_coords"]
        tokens_2_coords = tokens_2["center_coords"]

        dists = cdist(tokens_1_coords, tokens_2_coords)
        cuttoff = dists < const.interface_cutoff

        # In rare cases, the interface cuttoff is slightly
        # too small, then we slightly expand it if it happens
        if not np.any(cuttoff):
            cuttoff = dists < (const.interface_cutoff + 5.0)

        tokens_1 = tokens_1[np.any(cuttoff, axis=1)]
        tokens_2 = tokens_2[np.any(cuttoff, axis=0)]

        # Select random token
        candidates = np.concatenate([tokens_1, tokens_2])
        query = pick_random_token(candidates, random)

    return query


def pick_initial_crop_token(
    tokens: np.ndarray,
    initial_crop: list[int],
    random: np.random.Generator,
) -> np.ndarray:
    """Pick a random token from the initial crop.

    Parameters
    ----------
    tokens : np.ndarray
        The token data.
    initial_crop : list[int]
        The initial crop.
    random : np.random.Generator
        The random state for reproducibility.

    Returns
    -------
    np.ndarray

    """
    # Compute crop centroid
    crop_centroid = np.mean(tokens[initial_crop]["center_coords"], axis=0)

    # Compute distances to all tokens
    dists = cdist(tokens["center_coords"], crop_centroid[None])

    # Pick the closest token
    return tokens[np.argmin(dists[:, 0])]


class BoltzCropper(Cropper):
    """Interpolate between contiguous and spatial crops."""

    def __init__(
        self,
        min_neighborhood: int = 0,
        max_neighborhood: int = 40,
        dna_double_helix: bool = False,
    ) -> None:
        """Initialize the cropper.

        Modulates the type of cropping to be performed.
        Smaller neighborhoods result in more spatial
        cropping. Larger neighborhoods result in more
        continuous cropping. A mix can be achieved by
        providing a list of sizes from which to sample.

        Parameters
        ----------
        min_neighborhood : int
            The minimum neighborhood size, by default 0.
        max_neighborhood : int
            The maximum neighborhood size, by default 40.
        dna_double_helix : bool
            Whether to use DNA double helix cropping, by default False.

        """
        self.neighborhood_sizes = list(range(min_neighborhood, max_neighborhood + 1, 2))
        self.dna_double_helix = dna_double_helix

    def crop(  # noqa: PLR0915
        self,
        data: TokenizedTraining,
        max_tokens: int,
        random: np.random.Generator,
        chain_id: Optional[int] = None,
        interface_id: Optional[int] = None,
        max_atoms: Optional[int] = None,
        return_indices: bool = False,
        initial_crop: Optional[list[int]] = None,
    ) -> TokenizedTraining:
        """Crop the data to a maximum number of tokens.

        Parameters
        ----------
        data : Tokenized
            The tokenized data.
        max_tokens : int
            The maximum number of tokens to crop.
        random : np.random.Generator
            The random state for reproducibility.
        max_atoms : Optional[int]
            The maximum number of atoms to consider.

        Returns
        -------
        TokenizedTraining
            The cropped data.

        """
        # Check inputs
        if chain_id is not None and interface_id is not None:
            msg = "Only one of chain_id or interface_id can be provided."
            raise ValueError(msg)

        # Randomly select a neighborhood size
        neighborhood_size = random.choice(self.neighborhood_sizes)

        # Get token data
        token_data = data.tokens
        token_bonds = data.bonds
        mask = data.structure.mask
        chains = data.structure.chains
        interfaces = data.structure.interfaces

        # Filter to valid chains
        valid_chains = chains[mask]

        # Filter to valid interfaces
        valid_interfaces = interfaces
        valid_interfaces = valid_interfaces[mask[valid_interfaces["chain_1"]]]
        valid_interfaces = valid_interfaces[mask[valid_interfaces["chain_2"]]]

        # Filter to resolved tokens
        valid_tokens = token_data[token_data["resolved_mask"]]

        # Check if we have any valid tokens
        if not valid_tokens.size:
            msg = "No valid tokens in structure"
            raise ValueError(msg)

        # Pick a random token, chain, or interface
        if initial_crop is not None:
            query = pick_initial_crop_token(token_data, initial_crop, random)
        elif chain_id is not None:
            query = pick_chain_token(valid_tokens, chain_id, random)
        elif interface_id is not None:
            interface = interfaces[interface_id]
            query = pick_interface_token(valid_tokens, interface, random)
        elif valid_interfaces.size:
            idx = random.randint(len(valid_interfaces))
            interface = valid_interfaces[idx]
            query = pick_interface_token(valid_tokens, interface, random)
        else:
            idx = random.randint(len(valid_chains))
            chain_id = valid_chains[idx]["asym_id"]
            query = pick_chain_token(valid_tokens, chain_id, random)

        # Sort all tokens by distance to query_coords
        dists = valid_tokens["center_coords"] - query["center_coords"]
        indices = np.argsort(np.linalg.norm(dists, axis=1))

        # Select cropped indices
        cropped: set[int] = set()
        total_atoms = 0

        if initial_crop is not None:
            cropped.update(initial_crop)
            total_atoms = sum(token_data[idx]["atom_num"] for idx in initial_crop)

        for idx in indices:
            # Get the token
            token = valid_tokens[idx]

            neighborhood_size_to_use = neighborhood_size
            center_tokens_to_use = [token]
            new_tokens_acc = []

            # If it is a DNA double helix we may change this
            if (
                self.dna_double_helix
                and token["mol_type"] == const.chain_type_ids["DNA"]
            ):
                base_coords = data.structure.atoms["coords"][
                    token["atom_idx"] : token["atom_idx"] + token["atom_num"]
                ]
                base_is_present = data.structure.atoms["is_present"][
                    token["atom_idx"] : token["atom_idx"] + token["atom_num"]
                ]
                base_coords = base_coords[base_is_present]

                best_dist = 1e9
                best_other_token = None

                for other_token in valid_tokens:
                    if (
                        other_token["mol_type"] == const.chain_type_ids["DNA"]
                        and other_token["asym_id"] != token["asym_id"]
                    ):
                        other_base_coords = data.structure.atoms["coords"][
                            other_token["atom_idx"] : other_token["atom_idx"]
                            + other_token["atom_num"]
                        ]
                        other_base_is_present = data.structure.atoms["is_present"][
                            other_token["atom_idx"] : other_token["atom_idx"]
                            + other_token["atom_num"]
                        ]
                        other_base_coords = other_base_coords[other_base_is_present]

                        dist = np.min(cdist(base_coords, other_base_coords))
                        if dist < best_dist:
                            best_dist = dist
                            best_other_token = other_token

                if best_dist < 3.0:
                    center_tokens_to_use.append(best_other_token)
                    neighborhood_size_to_use = neighborhood_size_to_use // 2

            for center_token in center_tokens_to_use:
                # Get all tokens from this chain
                chain_tokens = token_data[
                    token_data["asym_id"] == center_token["asym_id"]
                ]

                # Pick the whole chain if possible, otherwise select
                # a contiguous subset centered at the query token
                if len(chain_tokens) <= neighborhood_size_to_use:
                    new_tokens = chain_tokens
                else:
                    # First limit to the maximum set of tokens, with the
                    # neighboorhood on both sides to handle edges. This
                    # is mostly for efficiency with the while loop below.
                    min_idx = center_token["res_idx"] - neighborhood_size_to_use
                    max_idx = center_token["res_idx"] + neighborhood_size_to_use

                    max_token_set = chain_tokens
                    max_token_set = max_token_set[max_token_set["res_idx"] >= min_idx]
                    max_token_set = max_token_set[max_token_set["res_idx"] <= max_idx]

                    # Start by adding just the query token
                    new_tokens = max_token_set[
                        max_token_set["res_idx"] == center_token["res_idx"]
                    ]

                    # Expand the neighborhood until we have enough tokens, one
                    # by one to handle some edge cases with non-standard chains.
                    # We switch to the res_idx instead of the token_idx to always
                    # include all tokens from modified residues or from ligands.
                    min_idx = max_idx = center_token["res_idx"]
                    target_size = min(neighborhood_size_to_use, max_token_set.size)
                    while new_tokens.size < target_size:
                        min_idx = min_idx - 1
                        max_idx = max_idx + 1
                        new_tokens = max_token_set
                        new_tokens = new_tokens[new_tokens["res_idx"] >= min_idx]
                        new_tokens = new_tokens[new_tokens["res_idx"] <= max_idx]

                new_tokens_acc.append(new_tokens)

            # Compute new tokens and new atoms
            new_tokens = np.concatenate(new_tokens_acc)
            new_indices = set(new_tokens["token_idx"]) - cropped
            new_tokens = token_data[list(new_indices)]
            new_atoms = np.sum(new_tokens["atom_num"])

            # Stop if we exceed the max number of tokens or atoms
            if (len(new_indices) > (max_tokens - len(cropped))) or (
                (max_atoms is not None) and ((total_atoms + new_atoms) > max_atoms)
            ):
                break

            # Add new indices
            cropped.update(new_indices)
            total_atoms += new_atoms

        # Get the cropped tokens sorted by index
        token_data = token_data[sorted(cropped)]

        # Only keep bonds within the cropped tokens
        indices = token_data["token_idx"]
        token_bonds = token_bonds[np.isin(token_bonds["token_1"], indices)]
        token_bonds = token_bonds[np.isin(token_bonds["token_2"], indices)]

        # Return the cropped tokens
        if return_indices:
            token_ids_mol = set(
                token_data[token_data["mol_type"] == 3]["token_idx"].tolist()
            )
            return replace(data, tokens=token_data, bonds=token_bonds), sorted(
                cropped - token_ids_mol
            )
        else:
            return replace(data, tokens=token_data, bonds=token_bonds)

    def crop_indices(  # noqa: PLR0915
        self,
        data: TokenizedTraining,
        cropped_indices: list[int],
    ) -> TokenizedTraining:
        token_data = data.tokens
        token_ids_mol = token_data[token_data["mol_type"] == 3]["token_idx"].tolist()  # noqa: PLR2004
        cropped_indices = sorted({*token_ids_mol, *cropped_indices})
        token_data = token_data[cropped_indices]
        indices = token_data["token_idx"]
        token_bonds = data.bonds
        token_bonds = token_bonds[np.isin(token_bonds["token_1"], indices)]
        token_bonds = token_bonds[np.isin(token_bonds["token_2"], indices)]

        return replace(data, tokens=token_data, bonds=token_bonds)
