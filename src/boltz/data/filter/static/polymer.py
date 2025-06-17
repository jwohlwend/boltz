import itertools
from dataclasses import dataclass

import numpy as np
from sklearn.neighbors import KDTree

from boltz.data import const
from boltz.data.filter.static.filter import StaticFilter
from boltz.data.types import Structure


class MinimumLengthFilter(StaticFilter):
    """Filter polymers based on their length.

    We use the number of resolved residues when considering
    the minimum, and the sequence length for the maximum.

    """

    def __init__(self, min_len: int = 4, max_len: int = 5000) -> None:
        """Initialize the filter.

        Parameters
        ----------
        min_len : float, optional
            The minimum allowed length.
        max_len : float, optional
            The maximum allowed length.

        """
        self._min = min_len
        self._max = max_len

    def filter(self, structure: Structure) -> np.ndarray:
        """Filter a chains based on their length.

        Parameters
        ----------
        structure : Structure
            The structure to filter chains from.

        Returns
        -------
        np.ndarray
            The chains to keep, as a boolean mask.

        """
        valid = np.ones(len(structure.chains), dtype=bool)

        for i, chain in enumerate(structure.chains):
            if chain["mol_type"] == const.chain_type_ids["NONPOLYMER"]:
                continue

            res_start = chain["res_idx"]
            res_end = res_start + chain["res_num"]
            residues = structure.residues[res_start:res_end]
            resolved = residues["is_present"].sum()

            if (resolved < self._min) or (resolved > self._max):
                valid[i] = 0

        return valid


class UnknownFilter(StaticFilter):
    """Filter proteins with all unknown residues."""

    def filter(self, structure: Structure) -> np.ndarray:
        """Filter proteins with all unknown residues.

        Parameters
        ----------
        structure : Structure
            The structure to filter chains from.

        Returns
        -------
        np.ndarray
            The chains to keep, as a boolean mask.

        """
        valid = np.ones(len(structure.chains), dtype=bool)
        unk_toks = {
            const.chain_type_ids["PROTEIN"]: const.unk_token_ids["PROTEIN"],
            const.chain_type_ids["DNA"]: const.unk_token_ids["DNA"],
            const.chain_type_ids["RNA"]: const.unk_token_ids["RNA"],
        }

        for i, chain in enumerate(structure.chains):
            if chain["mol_type"] == const.chain_type_ids["NONPOLYMER"]:
                continue

            res_start = chain["res_idx"]
            res_end = res_start + chain["res_num"]
            residues = structure.residues[res_start:res_end]

            unk_id = unk_toks[chain["mol_type"]]
            if np.all(residues["res_type"] == unk_id):
                valid[i] = 0

        return valid


class ConsecutiveCA(StaticFilter):
    """Filter proteins with consecutive CA atoms above a threshold."""

    def __init__(self, max_dist: int = 10.0) -> None:
        """Initialize the filter.

        Parameters
        ----------
        max_dist : float, optional
            The maximum allowed distance.

        """
        self._max_dist = max_dist

    def filter(self, structure: Structure) -> np.ndarray:
        """Filter protein if consecutive CA atoms above a threshold.

        Parameters
        ----------
        structure : Structure
            The structure to filter chains from.

        Returns
        -------
        np.ndarray
            The chains to keep, as a boolean mask.

        """
        valid = np.ones(len(structure.chains), dtype=bool)

        # Remove chain if consecutive CA atoms are above threshold
        for i, chain in enumerate(structure.chains):
            # Skip non-protein chains
            if chain["mol_type"] != const.chain_type_ids["PROTEIN"]:
                continue

            # Get residues
            res_start = chain["res_idx"]
            res_end = res_start + chain["res_num"]
            residues = structure.residues[res_start:res_end]

            # Get c-alphas
            ca_ids = residues["atom_center"]
            ca_atoms = structure.atoms[ca_ids]

            res_valid = residues["is_present"]
            ca_valid = ca_atoms["is_present"] & res_valid
            ca_coords = ca_atoms["coords"]

            # Compute distances between consecutive atoms
            dist = np.linalg.norm(ca_coords[1:] - ca_coords[:-1], axis=1)
            dist = dist > self._max_dist
            dist = dist[ca_valid[1:] & ca_valid[:-1]]

            # Remove the chain if any valid pair is above threshold
            if np.any(dist):
                valid[i] = 0

        return valid


@dataclass(frozen=True)
class Clash:
    """A clash between two chains."""

    chain: int
    other: int
    num_atoms: int
    num_clashes: int


class ClashingChainsFilter(StaticFilter):
    """A filter that filters clashing chains.

    Clashing chains are defined as those with >30% of atoms
    within 1.7 Å of an atom in another chain. If two chains
    are clashing with each other, the chain with the greater
    percentage of clashing atoms will be removed. If the same
    fraction of atoms are clashing, the chain with fewer total
    atoms is removed. If the chains have the same number of
    atoms, then the chain with the larger chain id is removed.

    """

    def __init__(self, dist: float = 1.7, freq: float = 0.3) -> None:
        """Initialize the filter.

        Parameters
        ----------
        dist : float, optional
            The maximum distance for a clash.
        freq : float, optional
            The maximum allowed frequency of clashes.

        """
        self._dist = dist
        self._freq = freq

    def filter(self, structure: Structure) -> np.ndarray:  # noqa: PLR0912, C901
        """Filter out clashing chains.

        Parameters
        ----------
        structure : Structure
            The structure to filter chains from.

        Returns
        -------
        np.ndarray
            The chains to keep, as a boolean mask.

        """
        num_chains = len(structure.chains)
        if num_chains < 2:  # noqa: PLR2004
            return np.ones(num_chains, dtype=bool)

        # Get unique chain pairs
        pairs = itertools.combinations(range(num_chains), 2)

        # Compute clashes
        clashes: list[Clash] = []
        for i, j in pairs:
            # Get the chains
            c1 = structure.chains[i]
            c2 = structure.chains[j]

            # Get the atoms from each chain
            c1_start = c1["atom_idx"]
            c2_start = c2["atom_idx"]
            c1_end = c1_start + c1["atom_num"]
            c2_end = c2_start + c2["atom_num"]

            atoms1 = structure.atoms[c1_start:c1_end]
            atoms2 = structure.atoms[c2_start:c2_end]
            atoms1 = atoms1[atoms1["is_present"]]
            atoms2 = atoms2[atoms2["is_present"]]

            # Skip if either chain has no atoms
            if len(atoms1) == 0 or len(atoms2) == 0:
                continue

            # Compute the number of clashes
            # Compute the distance matrix
            tree = KDTree(atoms1["coords"], metric="euclidean")
            query = tree.query_radius(atoms2["coords"], self._dist)

            c2_clashes = sum(len(neighbors) > 0 for neighbors in query)
            c1_clashes = len(set(itertools.chain.from_iterable(query)))

            # Save results
            if (c1_clashes / len(atoms1)) > self._freq:
                clashes.append(Clash(i, j, len(atoms1), c1_clashes))
            if (c2_clashes / len(atoms2)) > self._freq:
                clashes.append(Clash(j, i, len(atoms2), c2_clashes))

        # Compute indices to clash map
        removed = set()
        ids_to_clash = {(c.chain, c.other): c for c in clashes}

        # Filter out chains according to ruleset
        for clash in clashes:
            # If either is already removed, skip
            if clash.chain in removed or clash.other in removed:
                continue

            # Check if the two chains clash with each other
            other_clash = ids_to_clash.get((clash.other, clash.chain))
            if other_clash is not None:
                # Remove the chain with the most clashes
                clash1_freq = clash.num_clashes / clash.num_atoms
                clash2_freq = other_clash.num_clashes / other_clash.num_atoms
                if clash1_freq > clash2_freq:
                    removed.add(clash.chain)
                elif clash1_freq < clash2_freq:
                    removed.add(clash.other)

                # If same, remove the chain with fewer atoms
                elif clash.num_atoms < other_clash.num_atoms:
                    removed.add(clash.chain)
                elif clash.num_atoms > other_clash.num_atoms:
                    removed.add(clash.other)

                # If same, remove the chain with the larger chain id
                else:
                    removed.add(max(clash.chain, clash.other))

            # Otherwise, just remove the chain directly
            else:
                removed.add(clash.chain)

        # Remove the chains
        valid = np.ones(len(structure.chains), dtype=bool)
        for i in removed:
            valid[i] = 0

        return valid
