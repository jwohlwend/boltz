from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import click
import numpy as np
from Bio import Align
from chembl_structure_pipeline.exclude_flag import exclude_flag
from chembl_structure_pipeline.standardizer import standardize_mol
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem, HybridizationType
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.rdchem import BondStereo, Conformer, Mol
from rdkit.Chem.rdDistGeom import GetMoleculeBoundsMatrix
from rdkit.Chem.rdMolDescriptors import CalcNumHeavyAtoms
from scipy.optimize import linear_sum_assignment

from boltz.data import const
from boltz.data.mol import load_molecules
from boltz.data.parse.mmcif import parse_mmcif
from boltz.data.parse.pdb import parse_pdb

from boltz.data.types import (
    AffinityInfo,
    Atom,
    AtomV2,
    Bond,
    BondV2,
    Chain,
    ChainInfo,
    ChiralAtomConstraint,
    Connection,
    Coords,
    Ensemble,
    InferenceOptions,
    Interface,
    PlanarBondConstraint,
    PlanarRing5Constraint,
    PlanarRing6Constraint,
    RDKitBoundsConstraint,
    Record,
    Residue,
    ResidueConstraints,
    StereoBondConstraint,
    Structure,
    StructureInfo,
    StructureV2,
    Target,
    TemplateInfo,
)

####################################################################################################
# DATACLASSES
####################################################################################################


@dataclass(frozen=True)
class ParsedAtom:
    """A parsed atom object."""

    name: str
    element: int
    charge: int
    coords: tuple[float, float, float]
    conformer: tuple[float, float, float]
    is_present: bool
    chirality: int


@dataclass(frozen=True)
class ParsedBond:
    """A parsed bond object."""

    atom_1: int
    atom_2: int
    type: int


@dataclass(frozen=True)
class ParsedRDKitBoundsConstraint:
    """A parsed RDKit bounds constraint object."""

    atom_idxs: tuple[int, int]
    is_bond: bool
    is_angle: bool
    upper_bound: float
    lower_bound: float


@dataclass(frozen=True)
class ParsedChiralAtomConstraint:
    """A parsed chiral atom constraint object."""

    atom_idxs: tuple[int, int, int, int]
    is_reference: bool
    is_r: bool


@dataclass(frozen=True)
class ParsedStereoBondConstraint:
    """A parsed stereo bond constraint object."""

    atom_idxs: tuple[int, int, int, int]
    is_check: bool
    is_e: bool


@dataclass(frozen=True)
class ParsedPlanarBondConstraint:
    """A parsed planar bond constraint object."""

    atom_idxs: tuple[int, int, int, int, int, int]


@dataclass(frozen=True)
class ParsedPlanarRing5Constraint:
    """A parsed planar bond constraint object."""

    atom_idxs: tuple[int, int, int, int, int]


@dataclass(frozen=True)
class ParsedPlanarRing6Constraint:
    """A parsed planar bond constraint object."""

    atom_idxs: tuple[int, int, int, int, int, int]


@dataclass(frozen=True)
class ParsedResidue:
    """A parsed residue object."""

    name: str
    type: int
    idx: int
    atoms: list[ParsedAtom]
    bonds: list[ParsedBond]
    orig_idx: Optional[int]
    atom_center: int
    atom_disto: int
    is_standard: bool
    is_present: bool
    rdkit_bounds_constraints: Optional[list[ParsedRDKitBoundsConstraint]] = None
    chiral_atom_constraints: Optional[list[ParsedChiralAtomConstraint]] = None
    stereo_bond_constraints: Optional[list[ParsedStereoBondConstraint]] = None
    planar_bond_constraints: Optional[list[ParsedPlanarBondConstraint]] = None
    planar_ring_5_constraints: Optional[list[ParsedPlanarRing5Constraint]] = None
    planar_ring_6_constraints: Optional[list[ParsedPlanarRing6Constraint]] = None


@dataclass(frozen=True)
class ParsedChain:
    """A parsed chain object."""

    entity: str
    type: int
    residues: list[ParsedResidue]
    cyclic_period: int
    sequence: Optional[str] = None
    affinity: Optional[bool] = False
    affinity_mw: Optional[float] = None


@dataclass(frozen=True)
class Alignment:
    """A parsed alignment object."""

    query_st: int
    query_en: int
    template_st: int
    template_en: int


####################################################################################################
# HELPERS
####################################################################################################


def convert_atom_name(name: str) -> tuple[int, int, int, int]:
    """Convert an atom name to a standard format.

    Parameters
    ----------
    name : str
        The atom name.

    Returns
    -------
    Tuple[int, int, int, int]
        The converted atom name.

    """
    name = name.strip()
    name = [ord(c) - 32 for c in name]
    name = name + [0] * (4 - len(name))
    return tuple(name)


def compute_3d_conformer(mol: Mol, version: str = "v3") -> bool:
    """Generate 3D coordinates using EKTDG method.

    Taken from `pdbeccdutils.core.component.Component`.

    Parameters
    ----------
    mol: Mol
        The RDKit molecule to process
    version: str, optional
        The ETKDG version, defaults ot v3

    Returns
    -------
    bool
        Whether computation was successful.

    """
    if version == "v3":
        options = AllChem.ETKDGv3()
    elif version == "v2":
        options = AllChem.ETKDGv2()
    else:
        options = AllChem.ETKDGv2()

    options.clearConfs = False
    conf_id = -1

    try:
        conf_id = AllChem.EmbedMolecule(mol, options)

        if conf_id == -1:
            print(
                f"WARNING: RDKit ETKDGv3 failed to generate a conformer for molecule "
                f"{Chem.MolToSmiles(AllChem.RemoveHs(mol))}, so the program will start with random coordinates. "
                f"Note that the performance of the model under this behaviour was not tested."
            )
            options.useRandomCoords = True
            conf_id = AllChem.EmbedMolecule(mol, options)

        AllChem.UFFOptimizeMolecule(mol, confId=conf_id, maxIters=1000)

    except RuntimeError:
        pass  # Force field issue here
    except ValueError:
        pass  # sanitization issue here

    if conf_id != -1:
        conformer = mol.GetConformer(conf_id)
        conformer.SetProp("name", "Computed")
        conformer.SetProp("coord_generation", f"ETKDG{version}")

        return True

    return False


def get_conformer(mol: Mol) -> Conformer:
    """Retrieve an rdkit object for a deemed conformer.

    Inspired by `pdbeccdutils.core.component.Component`.

    Parameters
    ----------
    mol: Mol
        The molecule to process.

    Returns
    -------
    Conformer
        The desired conformer, if any.

    Raises
    ------
    ValueError
        If there are no conformers of the given tyoe.

    """
    # Try using the computed conformer
    for c in mol.GetConformers():
        try:
            if c.GetProp("name") == "Computed":
                return c
        except KeyError:  # noqa: PERF203
            pass

    # Fallback to the ideal coordinates
    for c in mol.GetConformers():
        try:
            if c.GetProp("name") == "Ideal":
                return c
        except KeyError:  # noqa: PERF203
            pass

    # Fallback to boltz2 format
    conf_ids = [int(conf.GetId()) for conf in mol.GetConformers()]
    if len(conf_ids) > 0:
        conf_id = conf_ids[0]
        conformer = mol.GetConformer(conf_id)
        return conformer

    msg = "Conformer does not exist."
    raise ValueError(msg)


def compute_geometry_constraints(mol: Mol, idx_map):
    if mol.GetNumAtoms() <= 1:
        return []

    # Ensure RingInfo is initialized
    mol.UpdatePropertyCache(strict=False)
    Chem.GetSymmSSSR(mol)  # Compute ring information

    bounds = GetMoleculeBoundsMatrix(
        mol,
        set15bounds=True,
        scaleVDW=True,
        doTriangleSmoothing=True,
        useMacrocycle14config=False,
    )
    bonds = set(
        tuple(sorted(b)) for b in mol.GetSubstructMatches(Chem.MolFromSmarts("*~*"))
    )
    angles = set(
        tuple(sorted([a[0], a[2]]))
        for a in mol.GetSubstructMatches(Chem.MolFromSmarts("*~*~*"))
    )

    constraints = []
    for i, j in zip(*np.triu_indices(mol.GetNumAtoms(), k=1)):
        if i in idx_map and j in idx_map:
            constraint = ParsedRDKitBoundsConstraint(
                atom_idxs=(idx_map[i], idx_map[j]),
                is_bond=tuple(sorted([i, j])) in bonds,
                is_angle=tuple(sorted([i, j])) in angles,
                upper_bound=bounds[i, j],
                lower_bound=bounds[j, i],
            )
            constraints.append(constraint)
    return constraints


def compute_chiral_atom_constraints(mol, idx_map):
    constraints = []
    if all([atom.HasProp("_CIPRank") for atom in mol.GetAtoms()]):
        for center_idx, orientation in Chem.FindMolChiralCenters(
            mol, includeUnassigned=False
        ):
            center = mol.GetAtomWithIdx(center_idx)
            neighbors = [
                (neighbor.GetIdx(), int(neighbor.GetProp("_CIPRank")))
                for neighbor in center.GetNeighbors()
            ]
            neighbors = sorted(
                neighbors, key=lambda neighbor: neighbor[1], reverse=True
            )
            neighbors = tuple(neighbor[0] for neighbor in neighbors)
            is_r = orientation == "R"

            if len(neighbors) > 4 or center.GetHybridization() != HybridizationType.SP3:
                continue

            atom_idxs = (*neighbors[:3], center_idx)
            if all(i in idx_map for i in atom_idxs):
                constraints.append(
                    ParsedChiralAtomConstraint(
                        atom_idxs=tuple(idx_map[i] for i in atom_idxs),
                        is_reference=True,
                        is_r=is_r,
                    )
                )

            if len(neighbors) == 4:
                for skip_idx in range(3):
                    chiral_set = neighbors[:skip_idx] + neighbors[skip_idx + 1 :]
                    if skip_idx % 2 == 0:
                        atom_idxs = chiral_set[::-1] + (center_idx,)
                    else:
                        atom_idxs = chiral_set + (center_idx,)
                    if all(i in idx_map for i in atom_idxs):
                        constraints.append(
                            ParsedChiralAtomConstraint(
                                atom_idxs=tuple(idx_map[i] for i in atom_idxs),
                                is_reference=False,
                                is_r=is_r,
                            )
                        )
    return constraints


def compute_stereo_bond_constraints(mol, idx_map):
    constraints = []
    if all([atom.HasProp("_CIPRank") for atom in mol.GetAtoms()]):
        for bond in mol.GetBonds():
            stereo = bond.GetStereo()
            if stereo in {BondStereo.STEREOE, BondStereo.STEREOZ}:
                start_atom_idx, end_atom_idx = (
                    bond.GetBeginAtomIdx(),
                    bond.GetEndAtomIdx(),
                )
                start_neighbors = [
                    (neighbor.GetIdx(), int(neighbor.GetProp("_CIPRank")))
                    for neighbor in mol.GetAtomWithIdx(start_atom_idx).GetNeighbors()
                    if neighbor.GetIdx() != end_atom_idx
                ]
                start_neighbors = sorted(
                    start_neighbors, key=lambda neighbor: neighbor[1], reverse=True
                )
                start_neighbors = [neighbor[0] for neighbor in start_neighbors]
                end_neighbors = [
                    (neighbor.GetIdx(), int(neighbor.GetProp("_CIPRank")))
                    for neighbor in mol.GetAtomWithIdx(end_atom_idx).GetNeighbors()
                    if neighbor.GetIdx() != start_atom_idx
                ]
                end_neighbors = sorted(
                    end_neighbors, key=lambda neighbor: neighbor[1], reverse=True
                )
                end_neighbors = [neighbor[0] for neighbor in end_neighbors]
                is_e = stereo == BondStereo.STEREOE

                atom_idxs = (
                    start_neighbors[0],
                    start_atom_idx,
                    end_atom_idx,
                    end_neighbors[0],
                )
                if all(i in idx_map for i in atom_idxs):
                    constraints.append(
                        ParsedStereoBondConstraint(
                            atom_idxs=tuple(idx_map[i] for i in atom_idxs),
                            is_check=True,
                            is_e=is_e,
                        )
                    )

                if len(start_neighbors) == 2 and len(end_neighbors) == 2:
                    atom_idxs = (
                        start_neighbors[1],
                        start_atom_idx,
                        end_atom_idx,
                        end_neighbors[1],
                    )
                    if all(i in idx_map for i in atom_idxs):
                        constraints.append(
                            ParsedStereoBondConstraint(
                                atom_idxs=tuple(idx_map[i] for i in atom_idxs),
                                is_check=False,
                                is_e=is_e,
                            )
                        )
    return constraints


def compute_flatness_constraints(mol, idx_map):
    planar_double_bond_smarts = Chem.MolFromSmarts("[C;X3;^2](*)(*)=[C;X3;^2](*)(*)")
    aromatic_ring_5_smarts = Chem.MolFromSmarts("[ar5^2]1[ar5^2][ar5^2][ar5^2][ar5^2]1")
    aromatic_ring_6_smarts = Chem.MolFromSmarts(
        "[ar6^2]1[ar6^2][ar6^2][ar6^2][ar6^2][ar6^2]1"
    )

    planar_double_bond_constraints = []
    aromatic_ring_5_constraints = []
    aromatic_ring_6_constraints = []
    for match in mol.GetSubstructMatches(planar_double_bond_smarts):
        if all(i in idx_map for i in match):
            planar_double_bond_constraints.append(
                ParsedPlanarBondConstraint(atom_idxs=tuple(idx_map[i] for i in match))
            )
    for match in mol.GetSubstructMatches(aromatic_ring_5_smarts):
        if all(i in idx_map for i in match):
            aromatic_ring_5_constraints.append(
                ParsedPlanarRing5Constraint(atom_idxs=tuple(idx_map[i] for i in match))
            )
    for match in mol.GetSubstructMatches(aromatic_ring_6_smarts):
        if all(i in idx_map for i in match):
            aromatic_ring_6_constraints.append(
                ParsedPlanarRing6Constraint(atom_idxs=tuple(idx_map[i] for i in match))
            )

    return (
        planar_double_bond_constraints,
        aromatic_ring_5_constraints,
        aromatic_ring_6_constraints,
    )


def get_global_alignment_score(query: str, template: str) -> float:
    """Align a sequence to a template.

    Parameters
    ----------
    query : str
        The query sequence.
    template : str
        The template sequence.

    Returns
    -------
    float
        The global alignment score.

    """
    aligner = Align.PairwiseAligner(scoring="blastp")
    aligner.mode = "global"
    score = aligner.align(query, template)[0].score
    return score


def get_local_alignments(query: str, template: str) -> list[Alignment]:
    """Align a sequence to a template.

    Parameters
    ----------
    query : str
        The query sequence.
    template : str
        The template sequence.

    Returns
    -------
    Alignment
        The alignment between the query and template.

    """
    aligner = Align.PairwiseAligner(scoring="blastp")
    aligner.mode = "local"
    aligner.open_gap_score = -1000
    aligner.extend_gap_score = -1000

    alignments = []
    for result in aligner.align(query, template):
        coordinates = result.coordinates
        alignment = Alignment(
            query_st=int(coordinates[0][0]),
            query_en=int(coordinates[0][1]),
            template_st=int(coordinates[1][0]),
            template_en=int(coordinates[1][1]),
        )
        alignments.append(alignment)

    return alignments


def get_template_records_from_search(
    template_id: str,
    chain_ids: list[str],
    sequences: dict[str, str],
    template_chain_ids: list[str],
    template_sequences: dict[str, str],
    force: bool = False,
    threshold: Optional[float] = None,
) -> list[TemplateInfo]:
    """Get template records from an alignment."""
    # Compute pairwise scores
    score_matrix = []
    for chain_id in chain_ids:
        row = []
        for template_chain_id in template_chain_ids:
            chain_seq = sequences[chain_id]
            template_seq = template_sequences[template_chain_id]
            score = get_global_alignment_score(chain_seq, template_seq)
            row.append(score)
        score_matrix.append(row)

    # Find optimal mapping
    row_ind, col_ind = linear_sum_assignment(score_matrix, maximize=True)

    # Get alignment records
    template_records = []

    for row_idx, col_idx in zip(row_ind, col_ind):
        chain_id = chain_ids[row_idx]
        template_chain_id = template_chain_ids[col_idx]
        chain_seq = sequences[chain_id]
        template_seq = template_sequences[template_chain_id]
        alignments = get_local_alignments(chain_seq, template_seq)

        for alignment in alignments:
            template_record = TemplateInfo(
                name=template_id,
                query_chain=chain_id,
                query_st=alignment.query_st,
                query_en=alignment.query_en,
                template_chain=template_chain_id,
                template_st=alignment.template_st,
                template_en=alignment.template_en,
                force=force,
                threshold=threshold,
            )
            template_records.append(template_record)

    return template_records


def get_template_records_from_matching(
    template_id: str,
    chain_ids: list[str],
    sequences: dict[str, str],
    template_chain_ids: list[str],
    template_sequences: dict[str, str],
    force: bool = False,
    threshold: Optional[float] = None,
) -> list[TemplateInfo]:
    """Get template records from a given matching."""
    template_records = []

    for chain_id, template_chain_id in zip(chain_ids, template_chain_ids):
        # Align the sequences
        chain_seq = sequences[chain_id]
        template_seq = template_sequences[template_chain_id]
        alignments = get_local_alignments(chain_seq, template_seq)
        for alignment in alignments:
            template_record = TemplateInfo(
                name=template_id,
                query_chain=chain_id,
                query_st=alignment.query_st,
                query_en=alignment.query_en,
                template_chain=template_chain_id,
                template_st=alignment.template_st,
                template_en=alignment.template_en,
                force=force,
                threshold=threshold,
            )
            template_records.append(template_record)

    return template_records


def get_mol(ccd: str, mols: dict, moldir: str) -> Mol:
    """Get mol from CCD code.

    Return mol with ccd from mols if it is in mols. Otherwise load it from moldir,
    add it to mols, and return the mol.
    """
    mol = mols.get(ccd)
    if mol is None:
        mol = load_molecules(moldir, [ccd])[ccd]
    return mol


####################################################################################################
# PARSING
####################################################################################################


def parse_ccd_residue(
    name: str, ref_mol: Mol, res_idx: int, drop_leaving_atoms: bool = False
) -> Optional[ParsedResidue]:
    """Parse an MMCIF ligand.

    First tries to get the SMILES string from the RCSB.
    Then, tries to infer atom ordering using RDKit.

    Parameters
    ----------
    name: str
        The name of the molecule to parse.
    ref_mol: Mol
        The reference molecule to parse.
    res_idx : int
        The residue index.

    Returns
    -------
    ParsedResidue, optional
       The output ParsedResidue, if successful.

    """
    unk_chirality = const.chirality_type_ids[const.unk_chirality_type]

    # Check if this is a single heavy atom CCD residue
    if CalcNumHeavyAtoms(ref_mol) == 1:
        # Remove hydrogens
        ref_mol = AllChem.RemoveHs(ref_mol, sanitize=False)

        pos = (0, 0, 0)
        ref_atom = ref_mol.GetAtoms()[0]
        chirality_type = const.chirality_type_ids.get(
            str(ref_atom.GetChiralTag()), unk_chirality
        )
        atom = ParsedAtom(
            name=ref_atom.GetProp("name"),
            element=ref_atom.GetAtomicNum(),
            charge=ref_atom.GetFormalCharge(),
            coords=pos,
            conformer=(0, 0, 0),
            is_present=True,
            chirality=chirality_type,
        )
        unk_prot_id = const.unk_token_ids["PROTEIN"]
        residue = ParsedResidue(
            name=name,
            type=unk_prot_id,
            atoms=[atom],
            bonds=[],
            idx=res_idx,
            orig_idx=None,
            atom_center=0,  # Placeholder, no center
            atom_disto=0,  # Placeholder, no center
            is_standard=False,
            is_present=True,
        )
        return residue

    # Get reference conformer coordinates
    conformer = get_conformer(ref_mol)

    # Parse each atom in order of the reference mol
    atoms = []
    atom_idx = 0
    idx_map = {}  # Used for bonds later

    for i, atom in enumerate(ref_mol.GetAtoms()):
        # Ignore Hydrogen atoms
        if atom.GetAtomicNum() == 1:
            continue

        # Get atom name, charge, element and reference coordinates
        atom_name = atom.GetProp("name")

        # Drop leaving atoms for non-canonical amino acids.
        if drop_leaving_atoms and int(atom.GetProp("leaving_atom")):
            continue

        charge = atom.GetFormalCharge()
        element = atom.GetAtomicNum()
        ref_coords = conformer.GetAtomPosition(atom.GetIdx())
        ref_coords = (ref_coords.x, ref_coords.y, ref_coords.z)
        chirality_type = const.chirality_type_ids.get(
            str(atom.GetChiralTag()), unk_chirality
        )

        # Get PDB coordinates, if any
        coords = (0, 0, 0)
        atom_is_present = True

        # Add atom to list
        atoms.append(
            ParsedAtom(
                name=atom_name,
                element=element,
                charge=charge,
                coords=coords,
                conformer=ref_coords,
                is_present=atom_is_present,
                chirality=chirality_type,
            )
        )
        idx_map[i] = atom_idx
        atom_idx += 1

    # Load bonds
    bonds = []
    unk_bond = const.bond_type_ids[const.unk_bond_type]
    for bond in ref_mol.GetBonds():
        idx_1 = bond.GetBeginAtomIdx()
        idx_2 = bond.GetEndAtomIdx()

        # Skip bonds with atoms ignored
        if (idx_1 not in idx_map) or (idx_2 not in idx_map):
            continue

        idx_1 = idx_map[idx_1]
        idx_2 = idx_map[idx_2]
        start = min(idx_1, idx_2)
        end = max(idx_1, idx_2)
        bond_type = bond.GetBondType().name
        bond_type = const.bond_type_ids.get(bond_type, unk_bond)
        bonds.append(ParsedBond(start, end, bond_type))

    rdkit_bounds_constraints = compute_geometry_constraints(ref_mol, idx_map)
    chiral_atom_constraints = compute_chiral_atom_constraints(ref_mol, idx_map)
    stereo_bond_constraints = compute_stereo_bond_constraints(ref_mol, idx_map)
    planar_bond_constraints, planar_ring_5_constraints, planar_ring_6_constraints = (
        compute_flatness_constraints(ref_mol, idx_map)
    )

    unk_prot_id = const.unk_token_ids["PROTEIN"]
    return ParsedResidue(
        name=name,
        type=unk_prot_id,
        atoms=atoms,
        bonds=bonds,
        idx=res_idx,
        atom_center=0,
        atom_disto=0,
        orig_idx=None,
        is_standard=False,
        is_present=True,
        rdkit_bounds_constraints=rdkit_bounds_constraints,
        chiral_atom_constraints=chiral_atom_constraints,
        stereo_bond_constraints=stereo_bond_constraints,
        planar_bond_constraints=planar_bond_constraints,
        planar_ring_5_constraints=planar_ring_5_constraints,
        planar_ring_6_constraints=planar_ring_6_constraints,
    )


def parse_polymer(
    sequence: list[str],
    raw_sequence: str,
    entity: str,
    chain_type: str,
    components: dict[str, Mol],
    cyclic: bool,
    mol_dir: Path,
) -> Optional[ParsedChain]:
    """Process a sequence into a chain object.

    Performs alignment of the full sequence to the polymer
    residues. Loads coordinates and masks for the atoms in
    the polymer, following the ordering in const.atom_order.

    Parameters
    ----------
    sequence : list[str]
        The full sequence of the polymer.
    entity : str
        The entity id.
    entity_type : str
        The entity type.
    components : dict[str, Mol]
        The preprocessed PDB components dictionary.

    Returns
    -------
    ParsedChain, optional
        The output chain, if successful.

    Raises
    ------
    ValueError
        If the alignment fails.

    """
    ref_res = set(const.tokens)
    unk_chirality = const.chirality_type_ids[const.unk_chirality_type]

    # Get coordinates and masks
    parsed = []
    for res_idx, res_name in enumerate(sequence):
        # Check if modified residue
        # Map MSE to MET
        res_corrected = res_name if res_name != "MSE" else "MET"

        # Handle non-standard residues
        if res_corrected not in ref_res:
            ref_mol = get_mol(res_corrected, components, mol_dir)
            residue = parse_ccd_residue(
                name=res_corrected,
                ref_mol=ref_mol,
                res_idx=res_idx,
                drop_leaving_atoms=True,
            )
            parsed.append(residue)
            continue

        # Load ref residue
        ref_mol = get_mol(res_corrected, components, mol_dir)
        ref_mol = AllChem.RemoveHs(ref_mol, sanitize=False)
        ref_conformer = get_conformer(ref_mol)

        # Only use reference atoms set in constants
        ref_name_to_atom = {a.GetProp("name"): a for a in ref_mol.GetAtoms()}
        ref_atoms = [ref_name_to_atom[a] for a in const.ref_atoms[res_corrected]]

        # Iterate, always in the same order
        atoms: list[ParsedAtom] = []

        for ref_atom in ref_atoms:
            # Get atom name
            atom_name = ref_atom.GetProp("name")
            idx = ref_atom.GetIdx()

            # Get conformer coordinates
            ref_coords = ref_conformer.GetAtomPosition(idx)
            ref_coords = (ref_coords.x, ref_coords.y, ref_coords.z)

            # Set 0 coordinate
            atom_is_present = True
            coords = (0, 0, 0)

            # Add atom to list
            atoms.append(
                ParsedAtom(
                    name=atom_name,
                    element=ref_atom.GetAtomicNum(),
                    charge=ref_atom.GetFormalCharge(),
                    coords=coords,
                    conformer=ref_coords,
                    is_present=atom_is_present,
                    chirality=const.chirality_type_ids.get(
                        str(ref_atom.GetChiralTag()), unk_chirality
                    ),
                )
            )

        atom_center = const.res_to_center_atom_id[res_corrected]
        atom_disto = const.res_to_disto_atom_id[res_corrected]
        parsed.append(
            ParsedResidue(
                name=res_corrected,
                type=const.token_ids[res_corrected],
                atoms=atoms,
                bonds=[],
                idx=res_idx,
                atom_center=atom_center,
                atom_disto=atom_disto,
                is_standard=True,
                is_present=True,
                orig_idx=None,
            )
        )

    if cyclic:
        cyclic_period = len(sequence)
    else:
        cyclic_period = 0

    # Return polymer object
    return ParsedChain(
        entity=entity,
        residues=parsed,
        type=chain_type,
        cyclic_period=cyclic_period,
        sequence=raw_sequence,
    )


def token_spec_to_ids(
    chain_name, residue_index_or_atom_name, chain_to_idx, atom_idx_map, chains
):
    if chains[chain_name].type == const.chain_type_ids["NONPOLYMER"]:
        # Non-polymer chains are indexed by atom name
        _, _, atom_idx = atom_idx_map[(chain_name, 0, residue_index_or_atom_name)]
        return (chain_to_idx[chain_name], atom_idx)
    else:
        # Polymer chains are indexed by residue index
        return chain_to_idx[chain_name], residue_index_or_atom_name - 1


def parse_boltz_schema(  # noqa: C901, PLR0915, PLR0912
    name: str,
    schema: dict,
    ccd: Mapping[str, Mol],
    mol_dir: Optional[Path] = None,
    boltz_2: bool = False,
) -> Target:
    """Parse a Boltz input yaml / json.

    The input file should be a dictionary with the following format:

    version: 1
    sequences:
        - protein:
            id: A
            sequence: "MADQLTEEQIAEFKEAFSLF"
            msa: path/to/msa1.a3m
        - protein:
            id: [B, C]
            sequence: "AKLSILPWGHC"
            msa: path/to/msa2.a3m
        - rna:
            id: D
            sequence: "GCAUAGC"
        - ligand:
            id: E
            smiles: "CC1=CC=CC=C1"
    constraints:
        - bond:
            atom1: [A, 1, CA]
            atom2: [A, 2, N]
        - pocket:
            binder: E
            contacts: [[B, 1], [B, 2]]
            max_distance: 6
        - contact:
            token1: [A, 1]
            token2: [B, 1]
            max_distance: 6
    templates:
        - cif: path/to/template.cif
    properties:
        - affinity:
            binder: E

    Parameters
    ----------
    name : str
        A name for the input.
    schema : dict
        The input schema.
    components : dict
        Dictionary of CCD components.
    mol_dir: Path
        Path to the directory containing the molecules.
    boltz2: bool
        Whether to parse the input for Boltz2.

    Returns
    -------
    Target
        The parsed target.

    """
    # Assert version 1
    version = schema.get("version", 1)
    if version != 1:
        msg = f"Invalid version {version} in input!"
        raise ValueError(msg)

    # Disable rdkit warnings
    blocker = rdBase.BlockLogs()  # noqa: F841

    # First group items that have the same type, sequence and modifications
    items_to_group = {}
    chain_name_to_entity_type = {}

    for item in schema["sequences"]:
        # Get entity type
        entity_type = next(iter(item.keys())).lower()
        if entity_type not in {"protein", "dna", "rna", "ligand"}:
            msg = f"Invalid entity type: {entity_type}"
            raise ValueError(msg)

        # Get sequence
        if entity_type in {"protein", "dna", "rna"}:
            seq = str(item[entity_type]["sequence"])
        elif entity_type == "ligand":
            assert "smiles" in item[entity_type] or "ccd" in item[entity_type]
            assert "smiles" not in item[entity_type] or "ccd" not in item[entity_type]
            if "smiles" in item[entity_type]:
                seq = str(item[entity_type]["smiles"])
            else:
                seq = str(item[entity_type]["ccd"])

        # Group items by entity
        items_to_group.setdefault((entity_type, seq), []).append(item)

        # Map chain names to entity types
        chain_names = item[entity_type]["id"]
        chain_names = [chain_names] if isinstance(chain_names, str) else chain_names
        for chain_name in chain_names:
            chain_name_to_entity_type[chain_name] = entity_type

    # Check if any affinity ligand is present
    affinity_ligands = set()
    properties = schema.get("properties", [])
    if properties and not boltz_2:
        msg = "Affinity prediction is only supported for Boltz2!"
        raise ValueError(msg)

    for prop in properties:
        prop_type = next(iter(prop.keys())).lower()
        if prop_type == "affinity":
            binder = prop["affinity"]["binder"]
            if not isinstance(binder, str):
                # TODO: support multi residue ligands and ccd's
                msg = "Binder must be a single chain."
                raise ValueError(msg)

            if binder not in chain_name_to_entity_type:
                msg = f"Could not find binder with name {binder} in the input!"
                raise ValueError(msg)

            if chain_name_to_entity_type[binder] != "ligand":
                msg = (
                    f"Chain {binder} is not a ligand! "
                    "Affinity is currently only supported for ligands."
                )
                raise ValueError(msg)

            affinity_ligands.add(binder)

    # Check only one affinity ligand is present
    if len(affinity_ligands) > 1:
        msg = "Only one affinity ligand is currently supported!"
        raise ValueError(msg)

    # Go through entities and parse them
    extra_mols: dict[str, Mol] = {}
    chains: dict[str, ParsedChain] = {}
    chain_to_msa: dict[str, str] = {}
    entity_to_seq: dict[str, str] = {}
    is_msa_custom = False
    is_msa_auto = False
    ligand_id = 1
    for entity_id, items in enumerate(items_to_group.values()):
        # Get entity type and sequence
        entity_type = next(iter(items[0].keys())).lower()

        # Get ids
        ids = []
        for item in items:
            if isinstance(item[entity_type]["id"], str):
                ids.append(item[entity_type]["id"])
            elif isinstance(item[entity_type]["id"], list):
                ids.extend(item[entity_type]["id"])

        # Check if any affinity ligand is present
        if len(ids) == 1:
            affinity = ids[0] in affinity_ligands
        elif (len(ids) > 1) and any(x in affinity_ligands for x in ids):
            msg = "Cannot compute affinity for a ligand that has multiple copies!"
            raise ValueError(msg)
        else:
            affinity = False

        # Ensure all the items share the same msa
        msa = -1
        if entity_type == "protein":
            # Get the msa, default to 0, meaning auto-generated
            msa = items[0][entity_type].get("msa", 0)
            if (msa is None) or (msa == ""):
                msa = 0

            # Check if all MSAs are the same within the same entity
            for item in items:
                item_msa = item[entity_type].get("msa", 0)
                if (item_msa is None) or (item_msa == ""):
                    item_msa = 0

                if item_msa != msa:
                    msg = "All proteins with the same sequence must share the same MSA!"
                    raise ValueError(msg)

            # Set the MSA, warn if passed in single-sequence mode
            if msa == "empty":
                msa = -1
                msg = (
                    "Found explicit empty MSA for some proteins, will run "
                    "these in single sequence mode. Keep in mind that the "
                    "model predictions will be suboptimal without an MSA."
                )
                click.echo(msg)

            if msa not in (0, -1):
                is_msa_custom = True
            elif msa == 0:
                is_msa_auto = True

        # Parse a polymer
        if entity_type in {"protein", "dna", "rna"}:
            # Get token map
            if entity_type == "rna":
                token_map = const.rna_letter_to_token
            elif entity_type == "dna":
                token_map = const.dna_letter_to_token
            elif entity_type == "protein":
                token_map = const.prot_letter_to_token
            else:
                msg = f"Unknown polymer type: {entity_type}"
                raise ValueError(msg)

            # Get polymer info
            chain_type = const.chain_type_ids[entity_type.upper()]
            unk_token = const.unk_token[entity_type.upper()]

            # Extract sequence
            raw_seq = items[0][entity_type]["sequence"]
            entity_to_seq[entity_id] = raw_seq

            # Convert sequence to tokens
            seq = [token_map.get(c, unk_token) for c in list(raw_seq)]

            # Apply modifications
            for mod in items[0][entity_type].get("modifications", []):
                code = mod["ccd"]
                idx = mod["position"] - 1  # 1-indexed
                seq[idx] = code

            cyclic = items[0][entity_type].get("cyclic", False)

            # Parse a polymer
            parsed_chain = parse_polymer(
                sequence=seq,
                raw_sequence=raw_seq,
                entity=entity_id,
                chain_type=chain_type,
                components=ccd,
                cyclic=cyclic,
                mol_dir=mol_dir,
            )

        # Parse a non-polymer
        elif (entity_type == "ligand") and "ccd" in (items[0][entity_type]):
            seq = items[0][entity_type]["ccd"]

            if isinstance(seq, str):
                seq = [seq]

            if affinity and len(seq) > 1:
                msg = "Cannot compute affinity for multi residue ligands!"
                raise ValueError(msg)

            residues = []
            affinity_mw = None
            for res_idx, code in enumerate(seq):
                # Get mol
                ref_mol = get_mol(code, ccd, mol_dir)

                if affinity:
                    affinity_mw = AllChem.Descriptors.MolWt(ref_mol)

                    # Add error and warning messaging when computing affinity with ligands too large
                    if ref_mol.GetNumAtoms() > 128:
                        msg = f"The ligand for affinity is too large, ligands with more than 128 atoms are not " \
                              f"supported in the affinity prediction module"
                        raise ValueError(msg)
                    elif ref_mol.GetNumAtoms() > 56:
                        print("WARNING: the ligand used for affinity calculation is larger than 56 heavy-atoms, which "
                              "was the maximum during training, therefore the affinity output might be inaccurate.")

                # Parse residue
                residue = parse_ccd_residue(
                    name=code,
                    ref_mol=ref_mol,
                    res_idx=res_idx,
                )
                residues.append(residue)

            # Create multi ligand chain
            parsed_chain = ParsedChain(
                entity=entity_id,
                residues=residues,
                type=const.chain_type_ids["NONPOLYMER"],
                cyclic_period=0,
                sequence=None,
                affinity=affinity,
                affinity_mw=affinity_mw,
            )

            assert not items[0][entity_type].get("cyclic", False), (
                "Cyclic flag is not supported for ligands"
            )

        elif (entity_type == "ligand") and ("smiles" in items[0][entity_type]):
            seq = items[0][entity_type]["smiles"]

            if affinity:
                seq = standardize(seq)

            mol = AllChem.MolFromSmiles(seq)
            mol = AllChem.AddHs(mol)

            # Set atom names
            canonical_order = AllChem.CanonicalRankAtoms(mol)
            Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
            for atom, can_idx in zip(mol.GetAtoms(), canonical_order):
                atom_name = atom.GetSymbol().upper() + str(can_idx + 1)
                if len(atom_name) > 4:
                    msg = (
                        f"{seq} has an atom with a name longer than "
                        f"4 characters: {atom_name}."
                    )
                    raise ValueError(msg)
                atom.SetProp("name", atom_name)

            success = compute_3d_conformer(mol)
            if not success:
                msg = f"Failed to compute 3D conformer for {seq}"
                raise ValueError(msg)

            mol_no_h = AllChem.RemoveHs(mol, sanitize=False)

            if affinity:
                # Add error and warning messaging when computing affinity with ligands too large
                if mol_no_h.GetNumAtoms() > 128:
                    msg = f"The ligand for affinity is too large, ligands with more than 128 atoms are not supported in the affinity prediction module"
                    raise ValueError(msg)
                elif mol_no_h.GetNumAtoms() > 56:
                    print("WARNING: the ligand used for affinity calculation is larger than 56 heavy-atoms, "
                          "which was the maximum during training, therefore the affinity output might be inaccurate.")

            affinity_mw = AllChem.Descriptors.MolWt(mol_no_h) if affinity else None
            extra_mols[f"LIG{ligand_id}"] = mol_no_h
            residue = parse_ccd_residue(
                name=f"LIG{ligand_id}",
                ref_mol=mol,
                res_idx=0,
            )

            ligand_id += 1
            parsed_chain = ParsedChain(
                entity=entity_id,
                residues=[residue],
                type=const.chain_type_ids["NONPOLYMER"],
                cyclic_period=0,
                sequence=None,
                affinity=affinity,
                affinity_mw=affinity_mw,
            )

            assert not items[0][entity_type].get("cyclic", False), (
                "Cyclic flag is not supported for ligands"
            )

        else:
            msg = f"Invalid entity type: {entity_type}"
            raise ValueError(msg)

        # Add as many chains as provided ids
        for item in items:
            ids = item[entity_type]["id"]
            if isinstance(ids, str):
                ids = [ids]
            for chain_name in ids:
                chains[chain_name] = parsed_chain
                chain_to_msa[chain_name] = msa

    # Check if msa is custom or auto
    if is_msa_custom and is_msa_auto:
        msg = "Cannot mix custom and auto-generated MSAs in the same input!"
        raise ValueError(msg)

    # If no chains parsed fail
    if not chains:
        msg = "No chains parsed!"
        raise ValueError(msg)

    # Create tables
    atom_data = []
    bond_data = []
    res_data = []
    chain_data = []
    protein_chains = set()
    affinity_info = None

    rdkit_bounds_constraint_data = []
    chiral_atom_constraint_data = []
    stereo_bond_constraint_data = []
    planar_bond_constraint_data = []
    planar_ring_5_constraint_data = []
    planar_ring_6_constraint_data = []

    # Convert parsed chains to tables
    atom_idx = 0
    res_idx = 0
    asym_id = 0
    sym_count = {}
    chain_to_idx = {}

    # Keep a mapping of (chain_name, residue_idx, atom_name) to atom_idx
    atom_idx_map = {}

    for asym_id, (chain_name, chain) in enumerate(chains.items()):
        # Compute number of atoms and residues
        res_num = len(chain.residues)
        atom_num = sum(len(res.atoms) for res in chain.residues)

        # Save protein chains for later
        if chain.type == const.chain_type_ids["PROTEIN"]:
            protein_chains.add(chain_name)

        # Add affinity info
        if chain.affinity and affinity_info is not None:
            msg = "Cannot compute affinity for multiple ligands!"
            raise ValueError(msg)

        if chain.affinity:
            affinity_info = AffinityInfo(
                chain_id=asym_id,
                mw=chain.affinity_mw,
            )

        # Find all copies of this chain in the assembly
        entity_id = int(chain.entity)
        sym_id = sym_count.get(entity_id, 0)
        chain_data.append(
            (
                chain_name,
                chain.type,
                entity_id,
                sym_id,
                asym_id,
                atom_idx,
                atom_num,
                res_idx,
                res_num,
                chain.cyclic_period,
            )
        )
        chain_to_idx[chain_name] = asym_id
        sym_count[entity_id] = sym_id + 1

        # Add residue, atom, bond, data
        for res in chain.residues:
            atom_center = atom_idx + res.atom_center
            atom_disto = atom_idx + res.atom_disto
            res_data.append(
                (
                    res.name,
                    res.type,
                    res.idx,
                    atom_idx,
                    len(res.atoms),
                    atom_center,
                    atom_disto,
                    res.is_standard,
                    res.is_present,
                )
            )

            if res.rdkit_bounds_constraints is not None:
                for constraint in res.rdkit_bounds_constraints:
                    rdkit_bounds_constraint_data.append(  # noqa: PERF401
                        (
                            tuple(
                                c_atom_idx + atom_idx
                                for c_atom_idx in constraint.atom_idxs
                            ),
                            constraint.is_bond,
                            constraint.is_angle,
                            constraint.upper_bound,
                            constraint.lower_bound,
                        )
                    )
            if res.chiral_atom_constraints is not None:
                for constraint in res.chiral_atom_constraints:
                    chiral_atom_constraint_data.append(  # noqa: PERF401
                        (
                            tuple(
                                c_atom_idx + atom_idx
                                for c_atom_idx in constraint.atom_idxs
                            ),
                            constraint.is_reference,
                            constraint.is_r,
                        )
                    )
            if res.stereo_bond_constraints is not None:
                for constraint in res.stereo_bond_constraints:
                    stereo_bond_constraint_data.append(  # noqa: PERF401
                        (
                            tuple(
                                c_atom_idx + atom_idx
                                for c_atom_idx in constraint.atom_idxs
                            ),
                            constraint.is_check,
                            constraint.is_e,
                        )
                    )
            if res.planar_bond_constraints is not None:
                for constraint in res.planar_bond_constraints:
                    planar_bond_constraint_data.append(  # noqa: PERF401
                        (
                            tuple(
                                c_atom_idx + atom_idx
                                for c_atom_idx in constraint.atom_idxs
                            ),
                        )
                    )
            if res.planar_ring_5_constraints is not None:
                for constraint in res.planar_ring_5_constraints:
                    planar_ring_5_constraint_data.append(  # noqa: PERF401
                        (
                            tuple(
                                c_atom_idx + atom_idx
                                for c_atom_idx in constraint.atom_idxs
                            ),
                        )
                    )
            if res.planar_ring_6_constraints is not None:
                for constraint in res.planar_ring_6_constraints:
                    planar_ring_6_constraint_data.append(  # noqa: PERF401
                        (
                            tuple(
                                c_atom_idx + atom_idx
                                for c_atom_idx in constraint.atom_idxs
                            ),
                        )
                    )

            for bond in res.bonds:
                atom_1 = atom_idx + bond.atom_1
                atom_2 = atom_idx + bond.atom_2
                bond_data.append(
                    (
                        asym_id,
                        asym_id,
                        res_idx,
                        res_idx,
                        atom_1,
                        atom_2,
                        bond.type,
                    )
                )

            for atom in res.atoms:
                # Add atom to map
                atom_idx_map[(chain_name, res.idx, atom.name)] = (
                    asym_id,
                    res_idx,
                    atom_idx,
                )

                # Add atom to data
                atom_data.append(
                    (
                        atom.name,
                        atom.element,
                        atom.charge,
                        atom.coords,
                        atom.conformer,
                        atom.is_present,
                        atom.chirality,
                    )
                )
                atom_idx += 1

            res_idx += 1

    # Parse constraints
    connections = []
    pocket_constraints = []
    contact_constraints = []
    constraints = schema.get("constraints", [])
    for constraint in constraints:
        if "bond" in constraint:
            if "atom1" not in constraint["bond"] or "atom2" not in constraint["bond"]:
                msg = f"Bond constraint was not properly specified"
                raise ValueError(msg)

            c1, r1, a1 = tuple(constraint["bond"]["atom1"])
            c2, r2, a2 = tuple(constraint["bond"]["atom2"])
            c1, r1, a1 = atom_idx_map[(c1, r1 - 1, a1)]  # 1-indexed
            c2, r2, a2 = atom_idx_map[(c2, r2 - 1, a2)]  # 1-indexed
            connections.append((c1, c2, r1, r2, a1, a2))
        elif "pocket" in constraint:
            if (
                "binder" not in constraint["pocket"]
                or "contacts" not in constraint["pocket"]
            ):
                msg = f"Pocket constraint was not properly specified"
                raise ValueError(msg)

            if len(pocket_constraints) > 0 and not boltz_2:
                msg = f"Only one pocket binders is supported in Boltz-1!"
                raise ValueError(msg)

            max_distance = constraint["pocket"].get("max_distance", 6.0)
            if max_distance != 6.0 and not boltz_2:
                msg = f"Max distance != 6.0 is not supported in Boltz-1!"
                raise ValueError(msg)

            binder = constraint["pocket"]["binder"]
            binder = chain_to_idx[binder]

            contacts = []
            for chain_name, residue_index_or_atom_name in constraint["pocket"][
                "contacts"
            ]:
                contact = token_spec_to_ids(
                    chain_name,
                    residue_index_or_atom_name,
                    chain_to_idx,
                    atom_idx_map,
                    chains,
                )
                contacts.append(contact)

            force = constraint["pocket"].get("force", False)
            pocket_constraints.append((binder, contacts, max_distance, force))
        elif "contact" in constraint:
            if (
                "token1" not in constraint["contact"]
                or "token2" not in constraint["contact"]
            ):
                msg = f"Contact constraint was not properly specified"
                raise ValueError(msg)

            if not boltz_2:
                msg = f"Contact constraint is not supported in Boltz-1!"
                raise ValueError(msg)

            max_distance = constraint["contact"].get("max_distance", 6.0)

            chain_name1, residue_index_or_atom_name1 = constraint["contact"]["token1"]
            token1 = token_spec_to_ids(
                chain_name1,
                residue_index_or_atom_name1,
                chain_to_idx,
                atom_idx_map,
                chains,
            )
            chain_name2, residue_index_or_atom_name2 = constraint["contact"]["token2"]
            token2 = token_spec_to_ids(
                chain_name2,
                residue_index_or_atom_name2,
                chain_to_idx,
                atom_idx_map,
                chains,
            )
            force = constraint["contact"].get("force", False)

            contact_constraints.append((token1, token2, max_distance, force))
        else:
            msg = f"Invalid constraint: {constraint}"
            raise ValueError(msg)

    # Get protein sequences in this YAML
    protein_seqs = {name: chains[name].sequence for name in protein_chains}

    # Parse templates
    template_schema = schema.get("templates", [])
    if template_schema and not boltz_2:
        msg = "Templates are not supported in Boltz 1.0!"
        raise ValueError(msg)

    templates = {}
    template_records = []
    for template in template_schema:
        if "cif" in template:
            path = template["cif"]
            pdb = False
        elif "pdb" in template:
            path = template["pdb"]
            pdb = True
        else:
            msg = "Template was not properly specified, missing CIF or PDB path!"
            raise ValueError(msg)

        template_id = Path(path).stem
        chain_ids = template.get("chain_id", None)
        template_chain_ids = template.get("template_id", None)

        # Check validity of input
        matched = False

        if chain_ids is not None and not isinstance(chain_ids, list):
            chain_ids = [chain_ids]
        if template_chain_ids is not None and not isinstance(template_chain_ids, list):
            template_chain_ids = [template_chain_ids]

        if (
            template_chain_ids is not None
            and chain_ids is not None
        ):
           
                if len(template_chain_ids) == len(chain_ids):
                     if len(template_chain_ids) > 0 and len(chain_ids) > 0:
                        matched = True
                else:
                    msg = (
                        "When providing both the chain_id and template_id, the number of"
                        "template_ids provided must match the number of chain_ids!"
                    )
                    raise ValueError(msg)

        # Get relevant chains ids
        if chain_ids is None:
            chain_ids = list(protein_chains)

        for chain_id in chain_ids:
            if chain_id not in protein_chains:
                msg = (
                    f"Chain {chain_id} assigned for template"
                    f"{template_id} is not one of the protein chains!"
                )
                raise ValueError(msg)

        # Get relevant template chain ids
        if pdb:
            parsed_template = parse_pdb(
                path,
                mols=ccd,
                moldir=mol_dir,
                use_assembly=False,
                compute_interfaces=False,
            )
        else:
            parsed_template = parse_mmcif(
                path,
                mols=ccd,
                moldir=mol_dir,
                use_assembly=False,
                compute_interfaces=False,
            )
        template_proteins = {
            str(c["name"])
            for c in parsed_template.data.chains
            if c["mol_type"] == const.chain_type_ids["PROTEIN"]
        }
        if template_chain_ids is None:
            template_chain_ids = list(template_proteins)

        for chain_id in template_chain_ids:
            if chain_id not in template_proteins:
                msg = (
                    f"Template chain {chain_id} assigned for template"
                    f"{template_id} is not one of the protein chains!"
                )
                raise ValueError(msg)

        force = template.get("force", False)
        if force:
            if "threshold" in template:
                threshold = template["threshold"]
            else:
                msg = f"Template {template_id} must have threshold specified if force is set to True"
                raise ValueError(msg)
        else:
            threshold = float("inf")
        # Compute template records

        if matched:
            template_records.extend(
                get_template_records_from_matching(
                    template_id=template_id,
                    chain_ids=chain_ids,
                    sequences=protein_seqs,
                    template_chain_ids=template_chain_ids,
                    template_sequences=parsed_template.sequences,
                    force=force,
                    threshold=threshold,
                )
            )
        else:
            template_records.extend(
                get_template_records_from_search(
                    template_id=template_id,
                    chain_ids=chain_ids,
                    sequences=protein_seqs,
                    template_chain_ids=template_chain_ids,
                    template_sequences=parsed_template.sequences,
                    force=force,
                    threshold=threshold,
                )
            )
        # Save template
        templates[template_id] = parsed_template.data

    # Convert into datatypes
    residues = np.array(res_data, dtype=Residue)
    chains = np.array(chain_data, dtype=Chain)
    interfaces = np.array([], dtype=Interface)
    mask = np.ones(len(chain_data), dtype=bool)
    rdkit_bounds_constraints = np.array(
        rdkit_bounds_constraint_data, dtype=RDKitBoundsConstraint
    )
    chiral_atom_constraints = np.array(
        chiral_atom_constraint_data, dtype=ChiralAtomConstraint
    )
    stereo_bond_constraints = np.array(
        stereo_bond_constraint_data, dtype=StereoBondConstraint
    )
    planar_bond_constraints = np.array(
        planar_bond_constraint_data, dtype=PlanarBondConstraint
    )
    planar_ring_5_constraints = np.array(
        planar_ring_5_constraint_data, dtype=PlanarRing5Constraint
    )
    planar_ring_6_constraints = np.array(
        planar_ring_6_constraint_data, dtype=PlanarRing6Constraint
    )

    if boltz_2:
        atom_data = [(a[0], a[3], a[5], 0.0, 1.0) for a in atom_data]
        connections = [(*c, const.bond_type_ids["COVALENT"]) for c in connections]
        bond_data = bond_data + connections
        atoms = np.array(atom_data, dtype=AtomV2)
        bonds = np.array(bond_data, dtype=BondV2)
        coords = [(x,) for x in atoms["coords"]]
        coords = np.array(coords, Coords)
        ensemble = np.array([(0, len(coords))], dtype=Ensemble)
        data = StructureV2(
            atoms=atoms,
            bonds=bonds,
            residues=residues,
            chains=chains,
            interfaces=interfaces,
            mask=mask,
            coords=coords,
            ensemble=ensemble,
        )
    else:
        bond_data = [(b[4], b[5], b[6]) for b in bond_data]
        atom_data = [(convert_atom_name(a[0]), *a[1:]) for a in atom_data]
        atoms = np.array(atom_data, dtype=Atom)
        bonds = np.array(bond_data, dtype=Bond)
        connections = np.array(connections, dtype=Connection)
        data = Structure(
            atoms=atoms,
            bonds=bonds,
            residues=residues,
            chains=chains,
            connections=connections,
            interfaces=interfaces,
            mask=mask,
        )

    # Create metadata
    struct_info = StructureInfo(num_chains=len(chains))
    chain_infos = []
    for chain in chains:
        chain_info = ChainInfo(
            chain_id=int(chain["asym_id"]),
            chain_name=chain["name"],
            mol_type=int(chain["mol_type"]),
            cluster_id=-1,
            msa_id=chain_to_msa[chain["name"]],
            num_residues=int(chain["res_num"]),
            valid=True,
            entity_id=int(chain["entity_id"]),
        )
        chain_infos.append(chain_info)

    options = InferenceOptions(
        pocket_constraints=pocket_constraints, contact_constraints=contact_constraints
    )
    record = Record(
        id=name,
        structure=struct_info,
        chains=chain_infos,
        interfaces=[],
        inference_options=options,
        templates=template_records,
        affinity=affinity_info,
    )

    residue_constraints = ResidueConstraints(
        rdkit_bounds_constraints=rdkit_bounds_constraints,
        chiral_atom_constraints=chiral_atom_constraints,
        stereo_bond_constraints=stereo_bond_constraints,
        planar_bond_constraints=planar_bond_constraints,
        planar_ring_5_constraints=planar_ring_5_constraints,
        planar_ring_6_constraints=planar_ring_6_constraints,
    )
    return Target(
        record=record,
        structure=data,
        sequences=entity_to_seq,
        residue_constraints=residue_constraints,
        templates=templates,
        extra_mols=extra_mols,
    )


def standardize(smiles: str) -> Optional[str]:
    """Standardize a molecule and return its SMILES and a flag indicating whether the molecule is valid.
    This version has exception handling, which the original in mol-finder/data doesn't have. I didn't change the mol-finder/data
    since there are a lot of other functions that depend on it and I didn't want to break them.
    """
    LARGEST_FRAGMENT_CHOOSER = rdMolStandardize.LargestFragmentChooser()

    mol = Chem.MolFromSmiles(smiles, sanitize=False)

    exclude = exclude_flag(mol, includeRDKitSanitization=False)

    if exclude:
        raise ValueError("Molecule is excluded")

    # Standardize with ChEMBL data curation pipeline. During standardization, the molecule may be broken
    # Choose molecule with largest component
    mol = LARGEST_FRAGMENT_CHOOSER.choose(mol)
    # Standardize with ChEMBL data curation pipeline. During standardization, the molecule may be broken
    mol = standardize_mol(mol)
    smiles = Chem.MolToSmiles(mol)

    # Check if molecule can be parsed by RDKit (in rare cases, the molecule may be broken during standardization)
    if Chem.MolFromSmiles(smiles) is None:
        raise ValueError("Molecule is broken")

    return smiles
