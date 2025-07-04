from dataclasses import dataclass

import numpy as np

from boltz.data import const
from boltz.data.tokenize.tokenizer import Tokenizer
from boltz.data.types import Input, Token, TokenBond, Tokenized


@dataclass
class TokenData:
    """TokenData datatype."""

    token_idx: int
    atom_idx: int
    atom_num: int
    res_idx: int
    res_type: int
    sym_id: int
    asym_id: int
    entity_id: int
    mol_type: int
    center_idx: int
    disto_idx: int
    center_coords: np.ndarray
    disto_coords: np.ndarray
    resolved_mask: bool
    disto_mask: bool
    cyclic_period: int


def token_astuple(token: TokenData) -> tuple:
    """Convert a TokenData object to a tuple."""
    return (
        token.token_idx,
        token.atom_idx,
        token.atom_num,
        token.res_idx,
        token.res_type,
        token.sym_id,
        token.asym_id,
        token.entity_id,
        token.mol_type,
        token.center_idx,
        token.disto_idx,
        token.center_coords,
        token.disto_coords,
        token.resolved_mask,
        token.disto_mask,
        token.cyclic_period,
    )


class BoltzTokenizer(Tokenizer):
    """Tokenize an input structure for training."""

    def tokenize(self, data: Input) -> Tokenized:
        """Tokenize the input data.

        Parameters
        ----------
        data : Input
            The input data.

        Returns
        -------
        Tokenized
            The tokenized data.

        """
        # Get structure data
        struct = data.structure

        # Create token data
        token_data = []

        # Keep track of atom_idx to token_idx
        token_idx = 0
        atom_to_token = {}

        # Filter to valid chains only
        chains = struct.chains[struct.mask]

        for chain in chains:
            # Get residue indices
            res_start = chain["res_idx"]
            res_end = chain["res_idx"] + chain["res_num"]

            for res in struct.residues[res_start:res_end]:
                # Get atom indices
                atom_start = res["atom_idx"]
                atom_end = res["atom_idx"] + res["atom_num"]

                # Standard residues are tokens
                if res["is_standard"]:
                    # Get center and disto atoms
                    center = struct.atoms[res["atom_center"]]
                    disto = struct.atoms[res["atom_disto"]]

                    # Token is present if centers are
                    is_present = res["is_present"] & center["is_present"]
                    is_disto_present = res["is_present"] & disto["is_present"]

                    # Apply chain transformation
                    c_coords = center["coords"]
                    d_coords = disto["coords"]

                    # Create token
                    token = TokenData(
                        token_idx=token_idx,
                        atom_idx=res["atom_idx"],
                        atom_num=res["atom_num"],
                        res_idx=res["res_idx"],
                        res_type=res["res_type"],
                        sym_id=chain["sym_id"],
                        asym_id=chain["asym_id"],
                        entity_id=chain["entity_id"],
                        mol_type=chain["mol_type"],
                        center_idx=res["atom_center"],
                        disto_idx=res["atom_disto"],
                        center_coords=c_coords,
                        disto_coords=d_coords,
                        resolved_mask=is_present,
                        disto_mask=is_disto_present,
                        cyclic_period=chain["cyclic_period"],
                    )
                    token_data.append(token_astuple(token))

                    # Update atom_idx to token_idx
                    atom_to_token.update(
                        dict.fromkeys(range(atom_start, atom_end), token_idx))

                    token_idx += 1

                # Non-standard are tokenized per atom
                else:
                    # We use the unk protein token as res_type
                    unk_token = const.unk_token["PROTEIN"]
                    unk_id = const.token_ids[unk_token]

                    # Get atom coordinates
                    atom_data = struct.atoms[atom_start:atom_end]
                    atom_coords = atom_data["coords"]

                    # Tokenize each atom
                    for i, atom in enumerate(atom_data):
                        # Token is present if atom is
                        is_present = res["is_present"] & atom["is_present"]
                        index = atom_start + i

                        # Create token
                        token = TokenData(
                            token_idx=token_idx,
                            atom_idx=index,
                            atom_num=1,
                            res_idx=res["res_idx"],
                            res_type=unk_id,
                            sym_id=chain["sym_id"],
                            asym_id=chain["asym_id"],
                            entity_id=chain["entity_id"],
                            mol_type=chain["mol_type"],
                            center_idx=index,
                            disto_idx=index,
                            center_coords=atom_coords[i],
                            disto_coords=atom_coords[i],
                            resolved_mask=is_present,
                            disto_mask=is_present,
                            cyclic_period=chain[
                                "cyclic_period"
                            ],  # Enforced to be False in chain parser
                        )
                        token_data.append(token_astuple(token))

                        # Update atom_idx to token_idx
                        atom_to_token[index] = token_idx
                        token_idx += 1

        # Create token bonds
        token_bonds = []

        # Add atom-atom bonds from ligands
        for bond in struct.bonds:
            if (
                bond["atom_1"] not in atom_to_token
                or bond["atom_2"] not in atom_to_token
            ):
                continue
            token_bond = (
                atom_to_token[bond["atom_1"]],
                atom_to_token[bond["atom_2"]],
            )
            token_bonds.append(token_bond)

        # Add connection bonds (covalent)
        for conn in struct.connections:
            if (
                conn["atom_1"] not in atom_to_token
                or conn["atom_2"] not in atom_to_token
            ):
                continue
            token_bond = (
                atom_to_token[conn["atom_1"]],
                atom_to_token[conn["atom_2"]],
            )
            token_bonds.append(token_bond)

        token_data = np.array(token_data, dtype=Token)
        token_bonds = np.array(token_bonds, dtype=TokenBond)
        tokenized = Tokenized(
            tokens=token_data,
            bonds=token_bonds,
            structure=data.structure,
            msa=data.msa,
            record=data.record,
            residue_constraints=data.residue_constraints,
        )
        return tokenized
