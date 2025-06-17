import numpy as np

from boltz.data import const
from boltz.data.filter.static.filter import StaticFilter
from boltz.data.types import Structure


class ExcludedLigands(StaticFilter):
    """Filter excluded ligands."""

    def filter(self, structure: Structure) -> np.ndarray:
        """Filter excluded ligands.

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
            if chain["mol_type"] != const.chain_type_ids["NONPOLYMER"]:
                continue

            res_start = chain["res_idx"]
            res_end = res_start + chain["res_num"]
            residues = structure.residues[res_start:res_end]
            if any(res["name"] in const.ligand_exclusion for res in residues):
                valid[i] = 0

        return valid
