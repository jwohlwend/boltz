import re
from typing import Optional

from rdkit import Chem
from torch import Tensor

from boltz.data import const
from boltz.data.types import Structure


def to_pdb(
    structure: Structure,
    plddts: Optional[Tensor] = None,
    boltz2: bool = False,
) -> str:  # noqa: PLR0915
    """Write a structure into a PDB file.

    Parameters
    ----------
    structure : Structure
        The input structure

    Returns
    -------
    str
        the output PDB file

    """
    pdb_lines = []

    atom_index = 1
    atom_reindex_ter = []

    # Load periodic table for element mapping
    periodic_table = Chem.GetPeriodicTable()

    # Index into plddt tensor for current residue.
    res_num = 0
    # Tracks non-ligand plddt tensor indices,
    # Initializing to -1 handles case where ligand is resnum 0
    prev_polymer_resnum = -1
    # Tracks ligand indices.
    ligand_index_offset = 0

    # Add all atom sites.
    for chain in structure.chains:
        # We rename the chains in alphabetical order
        chain_idx = chain["asym_id"]
        chain_tag = chain["name"]

        res_start = chain["res_idx"]
        res_end = chain["res_idx"] + chain["res_num"]

        residues = structure.residues[res_start:res_end]
        for residue in residues:
            res_name = str(residue["name"])
            atom_start = residue["atom_idx"]
            atom_end = residue["atom_idx"] + residue["atom_num"]
            atoms = structure.atoms[atom_start:atom_end]
            atom_coords = atoms["coords"]
            for i, atom in enumerate(atoms):
                # This should not happen on predictions, but just in case.
                if not atom["is_present"]:
                    continue

                record_type = (
                    "ATOM"
                    if chain["mol_type"] != const.chain_type_ids["NONPOLYMER"]
                    else "HETATM"
                )
                name = str(atom["name"])
                if boltz2:
                    atom_name = str(atom["name"])
                    atom_key = re.sub(r"\d", "", atom_name)
                    if atom_key in const.ambiguous_atoms:
                        if isinstance(const.ambiguous_atoms[atom_key], str):
                            element = const.ambiguous_atoms[atom_key]
                        elif res_name in const.ambiguous_atoms[atom_key]:
                            element = const.ambiguous_atoms[atom_key][res_name]
                        else:
                            element = const.ambiguous_atoms[atom_key]["*"]
                    else:
                        element = atom_key[0]
                else:
                    atom_name = atom["name"]
                    atom_name = [chr(c + 32) for c in atom_name if c != 0]
                    atom_name = "".join(atom_name)
                    name = atom_name
                    element = periodic_table.GetElementSymbol(atom["element"].item())

                name = name if len(name) == 4 else f" {name}"  # noqa: PLR2004
                alt_loc = ""
                insertion_code = ""
                occupancy = 1.00
                element = element.upper()
                charge = ""
                residue_index = residue["res_idx"] + 1
                pos = atom_coords[i]
                res_name_3 = (
                    "LIG" if record_type == "HETATM" else str(residue["name"][:3])
                )

                if record_type != "HETATM":
                    # The current residue plddt is stored at the res_num index unless a ligand has previouly been added.
                    b_factor = (
                        100.00
                        if plddts is None
                        else round(
                            plddts[res_num + ligand_index_offset].item() * 100, 2
                        )
                    )
                    prev_polymer_resnum = res_num
                else:
                    # If not a polymer resnum, we can get index into plddts by adding offset relative to previous polymer resnum.
                    ligand_index_offset += 1
                    b_factor = (
                        100.00
                        if plddts is None
                        else round(
                            plddts[prev_polymer_resnum + ligand_index_offset].item()
                            * 100,
                            2,
                        )
                    )

                # PDB is a columnar format, every space matters here!
                atom_line = (
                    f"{record_type:<6}{atom_index:>5} {name:<4}{alt_loc:>1}"
                    f"{res_name_3:>3} {chain_tag:>1}"
                    f"{residue_index:>4}{insertion_code:>1}   "
                    f"{pos[0]:>8.3f}{pos[1]:>8.3f}{pos[2]:>8.3f}"
                    f"{occupancy:>6.2f}{b_factor:>6.2f}          "
                    f"{element:>2}{charge:>2}"
                )
                pdb_lines.append(atom_line)
                atom_reindex_ter.append(atom_index)
                atom_index += 1

            if record_type != "HETATM":
                res_num += 1

        should_terminate = chain_idx < (len(structure.chains) - 1)
        if should_terminate:
            # Close the chain.
            chain_end = "TER"
            chain_termination_line = (
                f"{chain_end:<6}{atom_index:>5}      "
                f"{res_name_3:>3} "
                f"{chain_tag:>1}{residue_index:>4}"
            )
            pdb_lines.append(chain_termination_line)
            atom_index += 1

    # Dump CONECT records.
    all_bonds = list(structure.bonds)
    if hasattr(structure, "connections"):
        all_bonds.extend(list(structure.connections))

    for bond in all_bonds:
        atom1 = structure.atoms[bond["atom_1"]]
        atom2 = structure.atoms[bond["atom_2"]]
        if not atom1["is_present"] or not atom2["is_present"]:
            continue
        atom1_idx = atom_reindex_ter[bond["atom_1"]]
        atom2_idx = atom_reindex_ter[bond["atom_2"]]
        conect_line = f"CONECT{atom1_idx:>5}{atom2_idx:>5}"
        pdb_lines.append(conect_line)

    pdb_lines.append("END")
    pdb_lines.append("")
    pdb_lines = [line.ljust(80) for line in pdb_lines]
    return "\n".join(pdb_lines)
