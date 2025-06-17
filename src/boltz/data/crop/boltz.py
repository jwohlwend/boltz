from dataclasses import replace
from typing import Optional

import numpy as np
from scipy.spatial.distance import cdist

from boltz.data import const
from boltz.data.crop.cropper import Cropper
from boltz.data.types import Tokenized


def pick_random_token(
    tokens: np.ndarray,
    random: np.random.RandomState,
) -> np.ndarray:
    """Pick a random token from the data.

    Parameters
    ----------
    tokens : np.ndarray
        The token data.
    random : np.ndarray
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
    random: np.random.RandomState,
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
    random: np.random.RandomState,
) -> np.ndarray:
    """Pick a random token from an interface.

    Parameters
    ----------
    tokens : np.ndarray
        The token data.
    interface : int
        The interface ID.
    random : np.ndarray
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


class BoltzCropper(Cropper):
    """Interpolate between contiguous and spatial crops."""

    def __init__(self, min_neighborhood: int = 0, max_neighborhood: int = 40) -> None:
        """Initialize the cropper.

        Modulates the type of cropping to be performed.
        Smaller neighborhoods result in more spatial
        cropping. Larger neighborhoods result in more
        continuous cropping. A mix can be achieved by
        providing a range over which to sample.

        Parameters
        ----------
        min_neighborhood : int
            The minimum neighborhood size, by default 0.
        max_neighborhood : int
            The maximum neighborhood size, by default 40.

        """
        sizes = list(range(min_neighborhood, max_neighborhood + 1, 2))
        self.neighborhood_sizes = sizes

    def crop(  # noqa: PLR0915
        self,
        data: Tokenized,
        max_tokens: int,
        random: np.random.RandomState,
        max_atoms: Optional[int] = None,
        chain_id: Optional[int] = None,
        interface_id: Optional[int] = None,
    ) -> Tokenized:
        """Crop the data to a maximum number of tokens.

        Parameters
        ----------
        data : Tokenized
            The tokenized data.
        max_tokens : int
            The maximum number of tokens to crop.
        random : np.random.RandomState
            The random state for reproducibility.
        max_atoms : int, optional
            The maximum number of atoms to consider.
        chain_id : int, optional
            The chain ID to crop.
        interface_id : int, optional
            The interface ID to crop.

        Returns
        -------
        Tokenized
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
        if chain_id is not None:
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
        for idx in indices:
            # Get the token
            token = valid_tokens[idx]

            # Get all tokens from this chain
            chain_tokens = token_data[token_data["asym_id"] == token["asym_id"]]

            # Pick the whole chain if possible, otherwise select
            # a contiguous subset centered at the query token
            if len(chain_tokens) <= neighborhood_size:
                new_tokens = chain_tokens
            else:
                # First limit to the maximum set of tokens, with the
                # neighborhood on both sides to handle edges. This
                # is mostly for efficiency with the while loop below.
                min_idx = token["res_idx"] - neighborhood_size
                max_idx = token["res_idx"] + neighborhood_size

                max_token_set = chain_tokens
                max_token_set = max_token_set[max_token_set["res_idx"] >= min_idx]
                max_token_set = max_token_set[max_token_set["res_idx"] <= max_idx]

                # Start by adding just the query token
                new_tokens = max_token_set[max_token_set["res_idx"] == token["res_idx"]]

                # Expand the neighborhood until we have enough tokens, one
                # by one to handle some edge cases with non-standard chains.
                # We switch to the res_idx instead of the token_idx to always
                # include all tokens from modified residues or from ligands.
                min_idx = max_idx = token["res_idx"]
                while new_tokens.size < neighborhood_size:
                    min_idx = min_idx - 1
                    max_idx = max_idx + 1
                    new_tokens = max_token_set
                    new_tokens = new_tokens[new_tokens["res_idx"] >= min_idx]
                    new_tokens = new_tokens[new_tokens["res_idx"] <= max_idx]

            # Compute new tokens and new atoms
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
        return replace(data, tokens=token_data, bonds=token_bonds)
