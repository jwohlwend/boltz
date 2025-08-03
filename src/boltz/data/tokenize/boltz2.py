from dataclasses import dataclass
from typing import Optional

import numpy as np

from boltz.data import const
from boltz.data.tokenize.tokenizer import Tokenizer
from boltz.data.types import (
    AffinityInfo,
    Input,
    StructureV2,
    TokenBondV2,
    Tokenized,
    TokenV2,
)


@dataclass
class TokenData:
    """TokenData datatype."""

    token_idx: int
    atom_idx: int
    atom_num: int
    res_idx: int
    res_type: int
    res_name: str
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
    modified: bool
    frame_rot: np.ndarray
    frame_t: np.ndarray
    frame_mask: bool
    cyclic_period: int
    affinity_mask: bool = False


def token_astuple(token: TokenData) -> tuple:
    """Convert a TokenData object to a tuple."""
    return (
        token.token_idx,
        token.atom_idx,
        token.atom_num,
        token.res_idx,
        token.res_type,
        token.res_name,
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
        token.modified,
        token.frame_rot,
        token.frame_t,
        token.frame_mask,
        token.cyclic_period,
        token.affinity_mask,
    )


def compute_frame(
    n: np.ndarray,
    ca: np.ndarray,
    c: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute the frame for a residue.

    Parameters
    ----------
    n : np.ndarray
        The N atom.
    ca : np.ndarray
        The C atom.
    c : np.ndarray
        The CA atom.

    Returns
    -------
    np.ndarray
        The frame.

    """
    v1 = c - ca
    v2 = n - ca
    e1 = v1 / (np.linalg.norm(v1) + 1e-10)
    u2 = v2 - e1 * np.dot(e1.T, v2)
    e2 = u2 / (np.linalg.norm(u2) + 1e-10)
    e3 = np.cross(e1, e2)
    rot = np.column_stack([e1, e2, e3])
    t = ca
    return rot, t


def get_unk_token(chain: np.ndarray) -> int:
    """Get the unk token for a residue.

    Parameters
    ----------
    chain : np.ndarray
        The chain.

    Returns
    -------
    int
        The unk token.

    """
    if chain["mol_type"] == const.chain_type_ids["DNA"]:
        unk_token = const.unk_token["DNA"]
    elif chain["mol_type"] == const.chain_type_ids["RNA"]:
        unk_token = const.unk_token["RNA"]
    else:
        unk_token = const.unk_token["PROTEIN"]

    res_id = const.token_ids[unk_token]
    return res_id


def tokenize_structure(  # noqa: C901, PLR0915
    struct: StructureV2,
    affinity: Optional[AffinityInfo] = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Tokenize a structure.

    Parameters
    ----------
    struct : StructureV2
        The structure to tokenize.
    affinity : Optional[AffinityInfo]
        The affinity information.

    Returns
    -------
    np.ndarray
        The tokenized data.
    np.ndarray
        The tokenized bonds.

    """
    # Create token data
    token_data = []

    # Keep track of atom_idx to token_idx
    token_idx = 0
    atom_to_token = {}

    # Filter to valid chains only
    chains = struct.chains[struct.mask]

    # Ensemble atom id start in coords table.
    # For cropper and other operations, hardcoded to 0th conformer.
    offset = struct.ensemble[0]["atom_coord_idx"]

    for chain in chains:
        # Get residue indices
        res_start = chain["res_idx"]
        res_end = chain["res_idx"] + chain["res_num"]
        is_protein = chain["mol_type"] == const.chain_type_ids["PROTEIN"]
        affinity_mask = (affinity is not None) and (
            int(chain["asym_id"]) == int(affinity.chain_id)
        )

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
                # Apply chain transformation
                c_coords = struct.coords[offset + res["atom_center"]]["coords"]
                d_coords = struct.coords[offset + res["atom_disto"]]["coords"]

                # If protein, compute frame, only used for templates
                frame_rot = np.eye(3).flatten()
                frame_t = np.zeros(3)
                frame_mask = False

                if is_protein:
                    # Get frame atoms
                    atom_st = res["atom_idx"]
                    atom_en = res["atom_idx"] + res["atom_num"]
                    atoms = struct.atoms[atom_st:atom_en]

                    # Atoms are always in the order N, CA, C
                    atom_n = atoms[0]
                    atom_ca = atoms[1]
                    atom_c = atoms[2]

                    # Compute frame and mask
                    frame_mask = atom_ca["is_present"]
                    frame_mask &= atom_c["is_present"]
                    frame_mask &= atom_n["is_present"]
                    frame_mask = bool(frame_mask)
                    if frame_mask:
                        frame_rot, frame_t = compute_frame(
                            atom_n["coords"],
                            atom_ca["coords"],
                            atom_c["coords"],
                        )
                        frame_rot = frame_rot.flatten()

                # Create token
                token = TokenData(
                    token_idx=token_idx,
                    atom_idx=res["atom_idx"],
                    atom_num=res["atom_num"],
                    res_idx=res["res_idx"],
                    res_type=res["res_type"],
                    res_name=res["name"],
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
                    modified=False,
                    frame_rot=frame_rot,
                    frame_t=frame_t,
                    frame_mask=frame_mask,
                    cyclic_period=chain["cyclic_period"],
                    affinity_mask=affinity_mask,
                )
                token_data.append(token_astuple(token))

                # Update atom_idx to token_idx
                atom_to_token.update(
                    dict.fromkeys(range(atom_start, atom_end), token_idx)
                )

                token_idx += 1

            # Non-standard are tokenized per atom
            elif chain["mol_type"] == const.chain_type_ids["NONPOLYMER"]:
                # We use the unk protein token as res_type
                unk_token = const.unk_token["PROTEIN"]
                unk_id = const.token_ids[unk_token]

                # Get atom coordinates
                atom_data = struct.atoms[atom_start:atom_end]
                atom_coords = struct.coords[offset + atom_start : offset + atom_end][
                    "coords"
                ]

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
                        res_name=res["name"],
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
                        modified=chain["mol_type"]
                        != const.chain_type_ids["NONPOLYMER"],
                        frame_rot=np.eye(3).flatten(),
                        frame_t=np.zeros(3),
                        frame_mask=False,
                        cyclic_period=chain["cyclic_period"],
                        affinity_mask=affinity_mask,
                    )
                    token_data.append(token_astuple(token))

                    # Update atom_idx to token_idx
                    atom_to_token[index] = token_idx
                    token_idx += 1

            # Modified residues in Boltz-2 are tokenized at residue level
            else:
                res_type = get_unk_token(chain)

                # Get center and disto atoms
                center = struct.atoms[res["atom_center"]]
                disto = struct.atoms[res["atom_disto"]]

                # Token is present if centers are
                is_present = res["is_present"] & center["is_present"]
                is_disto_present = res["is_present"] & disto["is_present"]

                # Apply chain transformation
                c_coords = struct.coords[offset + res["atom_center"]]["coords"]
                d_coords = struct.coords[offset + res["atom_disto"]]["coords"]

                # Create token
                token = TokenData(
                    token_idx=token_idx,
                    atom_idx=res["atom_idx"],
                    atom_num=res["atom_num"],
                    res_idx=res["res_idx"],
                    res_type=res_type,
                    res_name=res["name"],
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
                    modified=True,
                    frame_rot=np.eye(3).flatten(),
                    frame_t=np.zeros(3),
                    frame_mask=False,
                    cyclic_period=chain["cyclic_period"],
                    affinity_mask=affinity_mask,
                )
                token_data.append(token_astuple(token))

                # Update atom_idx to token_idx
                atom_to_token.update(
                    dict.fromkeys(range(atom_start, atom_end), token_idx)
                )

                token_idx += 1

    # Create token bonds
    token_bonds = []

    # Add atom-atom bonds from ligands
    for bond in struct.bonds:
        if bond["atom_1"] not in atom_to_token or bond["atom_2"] not in atom_to_token:
            continue
        token_bond = (
            atom_to_token[bond["atom_1"]],
            atom_to_token[bond["atom_2"]],
            bond["type"] + 1,
        )
        token_bonds.append(token_bond)

    token_data = np.array(token_data, dtype=TokenV2)
    token_bonds = np.array(token_bonds, dtype=TokenBondV2)

    return token_data, token_bonds


class Boltz2Tokenizer(Tokenizer):
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
        # Tokenize the structure
        token_data, token_bonds = tokenize_structure(
            data.structure, data.record.affinity
        )

        # Tokenize the templates
        if data.templates is not None:
            template_tokens = {}
            template_bonds = {}
            for template_id, template in data.templates.items():
                tmpl_token_data, tmpl_token_bonds = tokenize_structure(template)
                template_tokens[template_id] = tmpl_token_data
                template_bonds[template_id] = tmpl_token_bonds
        else:
            template_tokens = None
            template_bonds = None

        # Create the tokenized data
        tokenized = Tokenized(
            tokens=token_data,
            bonds=token_bonds,
            structure=data.structure,
            msa=data.msa,
            record=data.record,
            residue_constraints=data.residue_constraints,
            templates=data.templates,
            template_tokens=template_tokens,
            template_bonds=template_bonds,
            extra_mols=data.extra_mols,
        )
        return tokenized
