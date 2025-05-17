import numpy as np
from pathlib import Path
from rdkit import Chem
from rdkit.Geometry import Point3D
from boltz.data.types import Structure, Chain
from boltz.data import const

def write_ligand_sdf(chain: Chain, atoms: np.ndarray, bonds: np.ndarray, output_path: Path) -> None:
    """Write a ligand chain into an SDF file.

    Parameters
    ----------
    chain : Chain
        The ligand chain object from the structure.
    atoms : np.ndarray
        The atom slice corresponding to the ligand chain.
    bonds : np.ndarray
        The bonds for the ligand chain.
    output_path : Path
        The path to the output SDF file.
    """
    rwmol = Chem.RWMol()
    idx_map = {}  # map from structure atom index to RDKit atom index

    # Add atoms to the molecule
    for idx, atom in enumerate(atoms):
        elem = int(atom["element"])
        rd_idx = rwmol.AddAtom(Chem.Atom(elem))
        idx_map[idx] = rd_idx

    # Add bonds to the molecule
    for bond in bonds:
        a1 = idx_map[bond["atom_1"] - chain["atom_idx"]]
        a2 = idx_map[bond["atom_2"] - chain["atom_idx"]]
        btype = bond["type"]
        if btype == const.bond_type_ids["SINGLE"]:
            rd_btype = Chem.BondType.SINGLE
        elif btype == const.bond_type_ids["DOUBLE"]:
            rd_btype = Chem.BondType.DOUBLE
        elif btype == const.bond_type_ids["TRIPLE"]:
            rd_btype = Chem.BondType.TRIPLE
        elif btype == const.bond_type_ids["AROMATIC"]:
            rd_btype = Chem.BondType.AROMATIC
        else:
            rd_btype = Chem.BondType.UNSPECIFIED
        rwmol.AddBond(a1, a2, rd_btype)

    mol = rwmol.GetMol()

    # Set the atom coordinates
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i, atom in enumerate(atoms):
        x, y, z = atom["coords"]
        conf.SetAtomPosition(i, Point3D(float(x), float(y), float(z)))
    mol.AddConformer(conf, assignId=True)

    # Write the SDF file
    Chem.MolToMolFile(mol, str(output_path))
