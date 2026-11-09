"""Tests for the PDB writer."""

import numpy as np

from boltz.data.types import Atom, Bond, Chain, Residue, Structure
from boltz.data.write.pdb import to_pdb


def _enc_name(name: str) -> np.ndarray:
    """Encode an atom name as a 4-byte ASCII array."""
    arr = [ord(c) - 32 for c in name.ljust(4)[:4]]
    return np.array(arr, dtype="i1")


def test_to_pdb_missing_atoms_conect() -> None:
    """Test that missing atoms do not offset CONECT indices or raise IndexError."""
    # 4 atoms: N (idx 0), CA (idx 1, missing), C (idx 2), O (idx 3)
    atoms = np.zeros(4, dtype=Atom)
    atoms[0]["is_present"] = True
    atoms[0]["name"] = _enc_name("N")
    atoms[0]["element"] = 7
    atoms[1]["is_present"] = False
    atoms[1]["name"] = _enc_name("CA")
    atoms[1]["element"] = 6
    atoms[2]["is_present"] = True
    atoms[2]["name"] = _enc_name("C")
    atoms[2]["element"] = 6
    atoms[3]["is_present"] = True
    atoms[3]["name"] = _enc_name("O")
    atoms[3]["element"] = 8

    residues = np.zeros(1, dtype=Residue)
    residues[0]["name"] = "ALA"
    residues[0]["atom_idx"] = 0
    residues[0]["atom_num"] = 4
    residues[0]["res_idx"] = 0

    chains = np.zeros(1, dtype=Chain)
    chains[0]["name"] = "A"
    chains[0]["asym_id"] = 0
    chains[0]["res_idx"] = 0
    chains[0]["res_num"] = 1
    chains[0]["atom_idx"] = 0
    chains[0]["atom_num"] = 4

    bonds = np.zeros(2, dtype=Bond)
    # Bond between N (idx 0) and C (idx 2)
    bonds[0]["atom_1"] = 0
    bonds[0]["atom_2"] = 2
    # Bond between C (idx 2) and O (idx 3)
    bonds[1]["atom_1"] = 2
    bonds[1]["atom_2"] = 3

    struct = Structure(
        atoms=atoms,
        residues=residues,
        chains=chains,
        bonds=bonds,
        connections=np.zeros(0, dtype=Bond),
        interfaces=None,
        mask=None,
    )

    pdb_str = to_pdb(struct, boltz2=False)
    lines = pdb_str.strip().split("\n")

    # Present atoms get serials 1, 2, 3
    # N -> serial 1
    # C -> serial 2
    # O -> serial 3
    conect_lines = [line for line in lines if line.startswith("CONECT")]
    expected_conect_count = 2
    assert len(conect_lines) == expected_conect_count
    assert conect_lines[0].startswith("CONECT    1    2")
    assert conect_lines[1].startswith("CONECT    2    3")


def test_to_pdb_ter_non_zero_asym_id() -> None:
    """Test that TER lines are emitted even when chain asym_id does not start at 0."""
    # 2 chains with asym_id 1 and 2 (e.g. after filtering)
    atoms = np.zeros(2, dtype=Atom)
    atoms[0]["is_present"] = True
    atoms[0]["name"] = _enc_name("N")
    atoms[0]["element"] = 7
    atoms[1]["is_present"] = True
    atoms[1]["name"] = _enc_name("N")
    atoms[1]["element"] = 7

    residues = np.zeros(2, dtype=Residue)
    residues[0]["name"] = "ALA"
    residues[0]["atom_idx"] = 0
    residues[0]["atom_num"] = 1
    residues[0]["res_idx"] = 0
    residues[1]["name"] = "GLY"
    residues[1]["atom_idx"] = 1
    residues[1]["atom_num"] = 1
    residues[1]["res_idx"] = 1

    chains = np.zeros(2, dtype=Chain)
    chains[0]["name"] = "A"
    chains[0]["asym_id"] = 1
    chains[0]["res_idx"] = 0
    chains[0]["res_num"] = 1
    chains[0]["atom_idx"] = 0
    chains[0]["atom_num"] = 1

    chains[1]["name"] = "B"
    chains[1]["asym_id"] = 2
    chains[1]["res_idx"] = 1
    chains[1]["res_num"] = 1
    chains[1]["atom_idx"] = 1
    chains[1]["atom_num"] = 1

    struct = Structure(
        atoms=atoms,
        residues=residues,
        chains=chains,
        bonds=np.zeros(0, dtype=Bond),
        connections=np.zeros(0, dtype=Bond),
        interfaces=None,
        mask=None,
    )

    pdb_str = to_pdb(struct, boltz2=False)
    lines = pdb_str.strip().split("\n")
    ter_lines = [line for line in lines if line.startswith("TER")]
    assert len(ter_lines) == 1
    assert "ALA A   1" in ter_lines[0]
