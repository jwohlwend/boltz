import contextlib
from dataclasses import dataclass, replace
from typing import Optional

import gemmi
import numpy as np
from rdkit import rdBase
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import Conformer, Mol
from sklearn.neighbors import KDTree

from boltz.data import const
from boltz.data.types import (
    Atom,
    Bond,
    Chain,
    Connection,
    Interface,
    Residue,
    Structure,
    StructureInfo,
)

####################################################################################################
# DATACLASSES
####################################################################################################


@dataclass(frozen=True, slots=True)
class ParsedAtom:
    """A parsed atom object."""

    name: str
    element: int
    charge: int
    coords: tuple[float, float, float]
    conformer: tuple[float, float, float]
    is_present: bool
    chirality: int


@dataclass(frozen=True, slots=True)
class ParsedBond:
    """A parsed bond object."""

    atom_1: int
    atom_2: int
    type: int


@dataclass(frozen=True, slots=True)
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


@dataclass(frozen=True, slots=True)
class ParsedChain:
    """A parsed chain object."""

    name: str
    entity: str
    type: str
    residues: list[ParsedResidue]
    sequence: list[str]


@dataclass(frozen=True, slots=True)
class ParsedConnection:
    """A parsed connection object."""

    chain_1: str
    chain_2: str
    residue_index_1: int
    residue_index_2: int
    atom_index_1: str
    atom_index_2: str


@dataclass(frozen=True, slots=True)
class ParsedStructure:
    """A parsed structure object."""

    data: Structure
    info: StructureInfo
    covalents: list[int]


####################################################################################################
# HELPERS
####################################################################################################


def get_dates(block: gemmi.cif.Block) -> tuple[str, str, str]:
    """Get the deposited, released, and last revision dates.

    Parameters
    ----------
    block : gemmi.cif.Block
        The block to process.

    Returns
    -------
    str
        The deposited date.
    str
        The released date.
    str
        The last revision date.

    """
    deposited = "_pdbx_database_status.recvd_initial_deposition_date"
    revision = "_pdbx_audit_revision_history.revision_date"
    deposit_date = revision_date = release_date = ""
    with contextlib.suppress(Exception):
        deposit_date = block.find([deposited])[0][0]
        release_date = block.find([revision])[0][0]
        revision_date = block.find([revision])[-1][0]

    return deposit_date, release_date, revision_date


def get_resolution(block: gemmi.cif.Block) -> float:
    """Get the resolution from a gemmi structure.

    Parameters
    ----------
    block : gemmi.cif.Block
        The block to process.

    Returns
    -------
    float
        The resolution.

    """
    resolution = 0.0
    for res_key in (
        "_refine.ls_d_res_high",
        "_em_3d_reconstruction.resolution",
        "_reflns.d_resolution_high",
    ):
        with contextlib.suppress(Exception):
            resolution = float(block.find([res_key])[0].str(0))
            break
    return resolution


def get_method(block: gemmi.cif.Block) -> str:
    """Get the method from a gemmi structure.

    Parameters
    ----------
    block : gemmi.cif.Block
        The block to process.

    Returns
    -------
    str
        The method.

    """
    method = ""
    method_key = "_exptl.method"
    with contextlib.suppress(Exception):
        methods = block.find([method_key])
        method = ",".join([m.str(0).lower() for m in methods])

    return method


def convert_atom_name(name: str) -> tuple[int, int, int, int]:
    """Convert an atom name to a standard format.

    Parameters
    ----------
    name : str
        The atom name.

    Returns
    -------
    tuple[int, int, int, int]
        The converted atom name.

    """
    name = name.strip()
    name = [ord(c) - 32 for c in name]
    name = name + [0] * (4 - len(name))
    return tuple(name)


def get_unk_token(dtype: gemmi.PolymerType) -> str:
    """Get the unknown token for a given entity type.

    Parameters
    ----------
    dtype : gemmi.EntityType
        The entity type.

    Returns
    -------
    str
        The unknown token.

    """
    if dtype == gemmi.PolymerType.PeptideL:
        unk = const.unk_token["PROTEIN"]
    elif dtype == gemmi.PolymerType.Dna:
        unk = const.unk_token["DNA"]
    elif dtype == gemmi.PolymerType.Rna:
        unk = const.unk_token["RNA"]
    else:
        msg = f"Unknown polymer type: {dtype}"
        raise ValueError(msg)

    return unk


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
    for c in mol.GetConformers():
        try:
            if c.GetProp("name") == "Computed":
                return c
        except KeyError:  # noqa: PERF203
            pass

    for c in mol.GetConformers():
        try:
            if c.GetProp("name") == "Ideal":
                return c
        except KeyError:  # noqa: PERF203
            pass

    msg = "Conformer does not exist."
    raise ValueError(msg)


def compute_covalent_ligands(
    connections: list[gemmi.Connection],
    subchain_map: dict[tuple[str, int], str],
    entities: dict[str, gemmi.Entity],
) -> set[str]:
    """Compute the covalent ligands from a list of connections.

    Parameters
    ----------
    connections: List[gemmi.Connection]
        The connections to process.
    subchain_map: dict[tuple[str, int], str]
        The mapping from chain, residue index to subchain name.
    entities: dict[str, gemmi.Entity]
        The entities in the structure.

    Returns
    -------
    set
        The covalent ligand subchains.

    """
    # Get covalent chain ids
    covalent_chain_ids = set()
    for connection in connections:
        if connection.type.name != "Covale":
            continue

        # Map to correct subchain
        chain_1_name = connection.partner1.chain_name
        chain_2_name = connection.partner2.chain_name

        res_1_id = connection.partner1.res_id.seqid
        res_1_id = str(res_1_id.num) + str(res_1_id.icode).strip()

        res_2_id = connection.partner2.res_id.seqid
        res_2_id = str(res_2_id.num) + str(res_2_id.icode).strip()

        subchain_1 = subchain_map[(chain_1_name, res_1_id)]
        subchain_2 = subchain_map[(chain_2_name, res_2_id)]

        # If non-polymer or branched, add to set
        entity_1 = entities[subchain_1].entity_type.name
        entity_2 = entities[subchain_2].entity_type.name

        if entity_1 in {"NonPolymer", "Branched"}:
            covalent_chain_ids.add(subchain_1)
        if entity_2 in {"NonPolymer", "Branched"}:
            covalent_chain_ids.add(subchain_2)

    return covalent_chain_ids


def compute_interfaces(atom_data: np.ndarray, chain_data: np.ndarray) -> np.ndarray:
    """Compute the chain-chain interfaces from a gemmi structure.

    Parameters
    ----------
    atom_data : List[tuple]
        The atom data.
    chain_data : List[tuple]
        The chain data.

    Returns
    -------
    List[tuple[int, int]]
        The interfaces.

    """
    # Compute chain_id per atom
    chain_ids = []
    for idx, chain in enumerate(chain_data):
        chain_ids.extend([idx] * chain["atom_num"])
    chain_ids = np.array(chain_ids)

    # Filte to present atoms
    coords = atom_data["coords"]
    mask = atom_data["is_present"]

    coords = coords[mask]
    chain_ids = chain_ids[mask]

    # Compute the distance matrix
    tree = KDTree(coords, metric="euclidean")
    query = tree.query_radius(coords, const.atom_interface_cutoff)

    # Get unique chain pairs
    interfaces = set()
    for c1, pairs in zip(chain_ids, query):
        chains = np.unique(chain_ids[pairs])
        chains = chains[chains != c1]
        interfaces.update((c1, c2) for c2 in chains)

    # Get unique chain pairs
    interfaces = [(min(i, j), max(i, j)) for i, j in interfaces]
    interfaces = list({(int(i), int(j)) for i, j in interfaces})
    interfaces = np.array(interfaces, dtype=Interface)
    return interfaces


####################################################################################################
# PARSING
####################################################################################################


def parse_ccd_residue(  # noqa: PLR0915, C901
    name: str,
    components: dict[str, Mol],
    res_idx: int,
    gemmi_mol: Optional[gemmi.Residue] = None,
    is_covalent: bool = False,
) -> Optional[ParsedResidue]:
    """Parse an MMCIF ligand.

    First tries to get the SMILES string from the RCSB.
    Then, tries to infer atom ordering using RDKit.

    Parameters
    ----------
    name: str
        The name of the molecule to parse.
    components : dict
        The preprocessed PDB components dictionary.
    res_idx : int
        The residue index.
    gemmi_mol : Optional[gemmi.Residue]
        The PDB molecule, as a gemmi Residue object, if any.

    Returns
    -------
    ParsedResidue, optional
       The output ParsedResidue, if successful.

    """
    unk_chirality = const.chirality_type_ids[const.unk_chirality_type]
    # Check if we have a PDB structure for this residue,
    # it could be a missing residue from the sequence
    is_present = gemmi_mol is not None

    # Save original index (required for parsing connections)
    if is_present:
        orig_idx = gemmi_mol.seqid
        orig_idx = str(orig_idx.num) + str(orig_idx.icode).strip()
    else:
        orig_idx = None

    # Get reference component
    ref_mol = components[name]

    # Remove hydrogens
    ref_mol = AllChem.RemoveHs(ref_mol, sanitize=False)

    # Check if this is a single atom CCD residue
    if ref_mol.GetNumAtoms() == 1:
        pos = (0, 0, 0)
        if is_present:
            pos = (
                gemmi_mol[0].pos.x,
                gemmi_mol[0].pos.y,
                gemmi_mol[0].pos.z,
            )
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
            is_present=is_present,
            chirality=chirality_type,
        )
        unk_prot_id = const.unk_token_ids["PROTEIN"]
        residue = ParsedResidue(
            name=name,
            type=unk_prot_id,
            atoms=[atom],
            bonds=[],
            idx=res_idx,
            orig_idx=orig_idx,
            atom_center=0,  # Placeholder, no center
            atom_disto=0,  # Placeholder, no center
            is_standard=False,
            is_present=is_present,
        )
        return residue

    # If multi-atom, start by getting the PDB coordinates
    pdb_pos = {}
    if is_present:
        # Match atoms based on names
        for atom in gemmi_mol:
            atom: gemmi.Atom
            pos = (atom.pos.x, atom.pos.y, atom.pos.z)
            pdb_pos[atom.name] = pos

    # Get reference conformer coordinates
    conformer = get_conformer(ref_mol)

    # Parse each atom in order of the reference mol
    atoms = []
    atom_idx = 0
    idx_map = {}  # Used for bonds later

    for i, atom in enumerate(ref_mol.GetAtoms()):
        # Get atom name, charge, element and reference coordinates
        atom_name = atom.GetProp("name")
        charge = atom.GetFormalCharge()
        element = atom.GetAtomicNum()
        ref_coords = conformer.GetAtomPosition(atom.GetIdx())
        ref_coords = (ref_coords.x, ref_coords.y, ref_coords.z)
        chirality_type = const.chirality_type_ids.get(
            str(atom.GetChiralTag()), unk_chirality
        )

        # If the atom is a leaving atom, skip if not in the PDB and is_covalent
        if (
            int(atom.GetProp("leaving_atom")) == 1
            and is_covalent
            and (atom_name not in pdb_pos)
        ):
            continue

        # Get PDB coordinates, if any
        coords = pdb_pos.get(atom_name)
        if coords is None:
            atom_is_present = False
            coords = (0, 0, 0)
        else:
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

    unk_prot_id = const.unk_token_ids["PROTEIN"]
    return ParsedResidue(
        name=name,
        type=unk_prot_id,
        atoms=atoms,
        bonds=bonds,
        idx=res_idx,
        atom_center=0,
        atom_disto=0,
        orig_idx=orig_idx,
        is_standard=False,
        is_present=is_present,
    )


def parse_polymer(  # noqa: C901, PLR0915, PLR0912
    polymer: gemmi.ResidueSpan,
    polymer_type: gemmi.PolymerType,
    sequence: list[str],
    chain_id: str,
    entity: str,
    components: dict[str, Mol],
) -> Optional[ParsedChain]:
    """Process a gemmi Polymer into a chain object.

    Performs alignment of the full sequence to the polymer
    residues. Loads coordinates and masks for the atoms in
    the polymer, following the ordering in const.atom_order.

    Parameters
    ----------
    polymer : gemmi.ResidueSpan
        The polymer to process.
    polymer_type : gemmi.PolymerType
        The polymer type.
    sequence : str
        The full sequence of the polymer.
    chain_id : str
        The chain identifier.
    entity : str
        The entity name.
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
    # Get unknown chirality token
    unk_chirality = const.chirality_type_ids[const.unk_chirality_type]

    # Ignore microheterogenities (pick first)
    sequence = [gemmi.Entity.first_mon(item) for item in sequence]

    # Align full sequence to polymer residues
    # This is a simple way to handle all the different numbering schemes
    result = gemmi.align_sequence_to_polymer(
        sequence,
        polymer,
        polymer_type,
        gemmi.AlignmentScoring(),
    )

    # Get coordinates and masks
    i = 0
    ref_res = set(const.tokens)
    parsed = []
    for j, match in enumerate(result.match_string):
        # Get residue name from sequence
        res_name = sequence[j]

        # Check if we have a match in the structure
        res = None
        name_to_atom = {}

        if match == "|":
            # Get pdb residue
            res = polymer[i]
            name_to_atom = {a.name.upper(): a for a in res}

            # Double check the match
            if res.name != res_name:
                msg = "Alignment mismatch!"
                raise ValueError(msg)

            # Increment polymer index
            i += 1

        # Map MSE to MET, put the selenium atom in the sulphur column
        if res_name == "MSE":
            res_name = "MET"
            if "SE" in name_to_atom:
                name_to_atom["SD"] = name_to_atom["SE"]

        # Handle non-standard residues
        elif res_name not in ref_res:
            residue = parse_ccd_residue(
                name=res_name,
                components=components,
                res_idx=j,
                gemmi_mol=res,
                is_covalent=True,
            )
            parsed.append(residue)
            continue

        # Load regular residues
        ref_mol = components[res_name]
        ref_mol = AllChem.RemoveHs(ref_mol, sanitize=False)
        ref_conformer = get_conformer(ref_mol)

        # Only use reference atoms set in constants
        ref_name_to_atom = {a.GetProp("name"): a for a in ref_mol.GetAtoms()}
        ref_atoms = [ref_name_to_atom[a] for a in const.ref_atoms[res_name]]

        # Iterate, always in the same order
        atoms: list[ParsedAtom] = []

        for ref_atom in ref_atoms:
            # Get atom name
            atom_name = ref_atom.GetProp("name")
            idx = ref_atom.GetIdx()

            # Get conformer coordinates
            ref_coords = ref_conformer.GetAtomPosition(idx)
            ref_coords = (ref_coords.x, ref_coords.y, ref_coords.z)

            # Get coordinated from PDB
            if atom_name in name_to_atom:
                atom = name_to_atom[atom_name]
                atom_is_present = True
                coords = (atom.pos.x, atom.pos.y, atom.pos.z)
            else:
                atom_is_present = False
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

        # Fix naming errors in arginine residues where NH2 is
        # incorrectly assigned to be closer to CD than NH1
        if (res is not None) and (res_name == "ARG"):
            ref_atoms: list[str] = const.ref_atoms["ARG"]
            cd = atoms[ref_atoms.index("CD")]
            nh1 = atoms[ref_atoms.index("NH1")]
            nh2 = atoms[ref_atoms.index("NH2")]

            cd_coords = np.array(cd.coords)
            nh1_coords = np.array(nh1.coords)
            nh2_coords = np.array(nh2.coords)

            if all(atom.is_present for atom in (cd, nh1, nh2)) and (
                np.linalg.norm(nh1_coords - cd_coords)
                > np.linalg.norm(nh2_coords - cd_coords)
            ):
                atoms[ref_atoms.index("NH1")] = replace(nh1, coords=nh2.coords)
                atoms[ref_atoms.index("NH2")] = replace(nh2, coords=nh1.coords)

        # Add residue to parsed list
        if res is not None:
            orig_idx = res.seqid
            orig_idx = str(orig_idx.num) + str(orig_idx.icode).strip()
        else:
            orig_idx = None

        atom_center = const.res_to_center_atom_id[res_name]
        atom_disto = const.res_to_disto_atom_id[res_name]
        parsed.append(
            ParsedResidue(
                name=res_name,
                type=const.token_ids[res_name],
                atoms=atoms,
                bonds=[],
                idx=j,
                atom_center=atom_center,
                atom_disto=atom_disto,
                is_standard=True,
                is_present=res is not None,
                orig_idx=orig_idx,
            )
        )

    # Get polymer class
    if polymer_type == gemmi.PolymerType.PeptideL:
        chain_type = const.chain_type_ids["PROTEIN"]
    elif polymer_type == gemmi.PolymerType.Dna:
        chain_type = const.chain_type_ids["DNA"]
    elif polymer_type == gemmi.PolymerType.Rna:
        chain_type = const.chain_type_ids["RNA"]

    # Return polymer object
    return ParsedChain(
        name=chain_id,
        entity=entity,
        residues=parsed,
        type=chain_type,
        sequence=gemmi.one_letter_code(sequence),
    )


def parse_connection(
    connection: gemmi.Connection,
    chains: list[ParsedChain],
    subchain_map: dict[tuple[str, int], str],
) -> ParsedConnection:
    """Parse (covalent) connection from a gemmi Connection.

    Parameters
    ----------
    connections : gemmi.ConnectionList
        The connection list to parse.
    chains : List[Chain]
        The parsed chains.
    subchain_map : dict[tuple[str, int], str]
        The mapping from chain, residue index to subchain name.

    Returns
    -------
    List[Connection]
        The parsed connections.

    """
    # Map to correct subchains
    chain_1_name = connection.partner1.chain_name
    chain_2_name = connection.partner2.chain_name

    res_1_id = connection.partner1.res_id.seqid
    res_1_id = str(res_1_id.num) + str(res_1_id.icode).strip()

    res_2_id = connection.partner2.res_id.seqid
    res_2_id = str(res_2_id.num) + str(res_2_id.icode).strip()

    subchain_1 = subchain_map[(chain_1_name, res_1_id)]
    subchain_2 = subchain_map[(chain_2_name, res_2_id)]

    # Get chain indices
    chain_1 = next(chain for chain in chains if (chain.name == subchain_1))
    chain_2 = next(chain for chain in chains if (chain.name == subchain_2))

    # Get residue indices
    res_1_idx, res_1 = next(
        (idx, res)
        for idx, res in enumerate(chain_1.residues)
        if (res.orig_idx == res_1_id)
    )
    res_2_idx, res_2 = next(
        (idx, res)
        for idx, res in enumerate(chain_2.residues)
        if (res.orig_idx == res_2_id)
    )

    # Get atom indices
    atom_index_1 = next(
        idx
        for idx, atom in enumerate(res_1.atoms)
        if atom.name == connection.partner1.atom_name
    )
    atom_index_2 = next(
        idx
        for idx, atom in enumerate(res_2.atoms)
        if atom.name == connection.partner2.atom_name
    )

    conn = ParsedConnection(
        chain_1=subchain_1,
        chain_2=subchain_2,
        residue_index_1=res_1_idx,
        residue_index_2=res_2_idx,
        atom_index_1=atom_index_1,
        atom_index_2=atom_index_2,
    )

    return conn


def parse_mmcif(  # noqa: C901, PLR0915, PLR0912
    path: str,
    components: dict[str, Mol],
    use_assembly: bool = True,
) -> ParsedStructure:
    """Parse a structure in MMCIF format.

    Parameters
    ----------
    mmcif_file : PathLike
        Path to the MMCIF file.
    components: dict[str, Mol]
        The preprocessed PDB components dictionary.
    use_assembly: bool
        Whether to use the first assembly.

    Returns
    -------
    ParsedStructure
        The parsed structure.

    """
    # Disable rdkit warnings
    blocker = rdBase.BlockLogs()  # noqa: F841

    # Parse MMCIF input file
    block = gemmi.cif.read(str(path))[0]

    # Extract medatadata
    deposit_date, release_date, revision_date = get_dates(block)
    resolution = get_resolution(block)
    method = get_method(block)

    # Load structure object
    structure = gemmi.make_structure_from_block(block)

    # Clean up the structure
    structure.merge_chain_parts()
    structure.remove_waters()
    structure.remove_hydrogens()
    structure.remove_alternative_conformations()
    structure.remove_empty_chains()

    # Expand assembly 1
    if use_assembly and structure.assemblies:
        how = gemmi.HowToNameCopiedChain.AddNumber
        assembly_name = structure.assemblies[0].name
        structure.transform_to_assembly(assembly_name, how=how)

    # Parse entities
    # Create mapping from subchain id to entity
    entities: dict[str, gemmi.Entity] = {}
    entity_ids: dict[str, int] = {}
    for entity_id, entity in enumerate(structure.entities):
        entity: gemmi.Entity
        if entity.entity_type.name == "Water":
            continue
        for subchain_id in entity.subchains:
            entities[subchain_id] = entity
            entity_ids[subchain_id] = entity_id

    # Create mapping from chain, residue to subchains
    # since a Connection uses the chains and not subchins
    subchain_map = {}
    for chain in structure[0]:
        for residue in chain:
            seq_id = residue.seqid
            seq_id = str(seq_id.num) + str(seq_id.icode).strip()
            subchain_map[(chain.name, seq_id)] = residue.subchain

    # Find covalent ligands
    covalent_chain_ids = compute_covalent_ligands(
        connections=structure.connections,
        subchain_map=subchain_map,
        entities=entities,
    )

    # Parse chains
    chains: list[ParsedChain] = []
    chain_seqs = []
    for raw_chain in structure[0].subchains():
        # Check chain type
        subchain_id = raw_chain.subchain_id()
        entity: gemmi.Entity = entities[subchain_id]
        entity_type = entity.entity_type.name

        # Parse a polymer
        if entity_type == "Polymer":
            # Skip PeptideD, DnaRnaHybrid, Pna, Other
            if entity.polymer_type.name not in {
                "PeptideL",
                "Dna",
                "Rna",
            }:
                continue

            # Add polymer if successful
            parsed_polymer = parse_polymer(
                polymer=raw_chain,
                polymer_type=entity.polymer_type,
                sequence=entity.full_sequence,
                chain_id=subchain_id,
                entity=entity.name,
                components=components,
            )
            if parsed_polymer is not None:
                chains.append(parsed_polymer)
                chain_seqs.append(parsed_polymer.sequence)

        # Parse a non-polymer
        elif entity_type in {"NonPolymer", "Branched"}:
            # Skip UNL or other missing ligands
            if any(components.get(lig.name) is None for lig in raw_chain):
                continue

            residues = []
            for lig_idx, ligand in enumerate(raw_chain):
                # Check if ligand is covalent
                if entity_type == "Branched":
                    is_covalent = True
                else:
                    is_covalent = subchain_id in covalent_chain_ids

                ligand: gemmi.Residue
                residue = parse_ccd_residue(
                    name=ligand.name,
                    components=components,
                    res_idx=lig_idx,
                    gemmi_mol=ligand,
                    is_covalent=is_covalent,
                )
                residues.append(residue)

            if residues:
                chains.append(
                    ParsedChain(
                        name=subchain_id,
                        entity=entity.name,
                        residues=residues,
                        type=const.chain_type_ids["NONPOLYMER"],
                        sequence=None,
                    )
                )

    # If no chains parsed fail
    if not chains:
        msg = "No chains parsed!"
        raise ValueError(msg)

    # Parse covalent connections
    connections: list[ParsedConnection] = []
    for connection in structure.connections:
        # Skip non-covalent connections
        connection: gemmi.Connection
        if connection.type.name != "Covale":
            continue

        parsed_connection = parse_connection(
            connection=connection,
            chains=chains,
            subchain_map=subchain_map,
        )
        connections.append(parsed_connection)

    # Create tables
    atom_data = []
    bond_data = []
    res_data = []
    chain_data = []
    connection_data = []

    # Convert parsed chains to tables
    atom_idx = 0
    res_idx = 0
    asym_id = 0
    sym_count = {}
    chain_to_idx = {}
    res_to_idx = {}

    for asym_id, chain in enumerate(chains):
        # Compute number of atoms and residues
        res_num = len(chain.residues)
        atom_num = sum(len(res.atoms) for res in chain.residues)

        # Find all copies of this chain in the assembly
        entity_id = entity_ids[chain.name]
        sym_id = sym_count.get(entity_id, 0)
        chain_data.append(
            (
                chain.name,
                chain.type,
                entity_id,
                sym_id,
                asym_id,
                atom_idx,
                atom_num,
                res_idx,
                res_num,
            )
        )
        chain_to_idx[chain.name] = asym_id
        sym_count[entity_id] = sym_id + 1

        # Add residue, atom, bond, data
        for i, res in enumerate(chain.residues):
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
            res_to_idx[(chain.name, i)] = (res_idx, atom_idx)

            for bond in res.bonds:
                atom_1 = atom_idx + bond.atom_1
                atom_2 = atom_idx + bond.atom_2
                bond_data.append((atom_1, atom_2, bond.type))

            for atom in res.atoms:
                atom_data.append(
                    (
                        convert_atom_name(atom.name),
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

    # Convert connections to tables
    for conn in connections:
        chain_1_idx = chain_to_idx[conn.chain_1]
        chain_2_idx = chain_to_idx[conn.chain_2]
        res_1_idx, atom_1_offset = res_to_idx[(conn.chain_1, conn.residue_index_1)]
        res_2_idx, atom_2_offset = res_to_idx[(conn.chain_2, conn.residue_index_2)]
        atom_1_idx = atom_1_offset + conn.atom_index_1
        atom_2_idx = atom_2_offset + conn.atom_index_2
        connection_data.append(
            (
                chain_1_idx,
                chain_2_idx,
                res_1_idx,
                res_2_idx,
                atom_1_idx,
                atom_2_idx,
            )
        )

    # Convert into datatypes
    atoms = np.array(atom_data, dtype=Atom)
    bonds = np.array(bond_data, dtype=Bond)
    residues = np.array(res_data, dtype=Residue)
    chains = np.array(chain_data, dtype=Chain)
    connections = np.array(connection_data, dtype=Connection)
    mask = np.ones(len(chain_data), dtype=bool)

    # Compute interface chains (find chains with a heavy atom within 5A)
    interfaces = compute_interfaces(atoms, chains)

    # Return parsed structure
    info = StructureInfo(
        deposited=deposit_date,
        revised=revision_date,
        released=release_date,
        resolution=resolution,
        method=method,
        num_chains=len(chains),
        num_interfaces=len(interfaces),
    )

    data = Structure(
        atoms=atoms,
        bonds=bonds,
        residues=residues,
        chains=chains,
        connections=connections,
        interfaces=interfaces,
        mask=mask,
    )

    return ParsedStructure(data=data, info=info, covalents=[])
