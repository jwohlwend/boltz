import io
import re
from collections.abc import Iterator
from typing import Optional

import gemmi
import ihm
import modelcif
from modelcif import Assembly, AsymUnit, Entity, System, dumper
from modelcif.model import AbInitioModel, Atom, ModelGroup
from rdkit import Chem
from torch import Tensor

from boltz.data import const
from boltz.data.types import Structure


def _add_covalent_connections(text: str, structure: Structure, boltz2: bool) -> str:
    """Preserve explicit bond endpoints as mmCIF covalent connections."""
    links = (
        [b for b in structure.bonds if b["type"] == const.bond_type_ids["COVALENT"]]
        if boltz2
        else list(structure.connections)
    )
    links = [
        b
        for b in links
        if structure.atoms[b["atom_1"]]["is_present"]
        and structure.atoms[b["atom_2"]]["is_present"]
    ]
    if not links:
        return text
    document = gemmi.cif.read_string(text)
    block = document.sole_block()
    model = gemmi.make_structure_from_block(block)

    def address_key(address):
        return (
            address.chain_name,
            address.res_id.seqid.num,
            address.res_id.seqid.icode,
            address.atom_name,
        )

    existing_pairs = {
        frozenset((address_key(c.partner1), address_key(c.partner2)))
        for c in model.connections
        if c.type == gemmi.ConnectionType.Covale
    }
    existing_ids = {c.name for c in model.connections}
    model.connections.clear()
    for link in links:
        endpoints = []
        for side in (1, 2):
            chain = structure.chains[link[f"chain_{side}"]]
            residue = structure.residues[link[f"res_{side}"]]
            atom = structure.atoms[link[f"atom_{side}"]]
            name = (
                str(atom["name"])
                if boltz2
                else "".join(chr(int(c) + 32) for c in atom["name"] if c != 0)
            )
            endpoints.append(
                gemmi.AtomAddress(
                    str(chain["name"]),
                    gemmi.SeqId(int(residue["res_idx"]) + 1, " "),
                    str(residue["name"]),
                    name,
                )
            )
        pair = frozenset(address_key(a) for a in endpoints)
        if pair in existing_pairs:
            continue
        existing_pairs.add(pair)
        index = len(existing_ids) + 1
        while f"boltz_bond_{index}" in existing_ids:
            index += 1
        connection = gemmi.Connection()
        connection.name = f"boltz_bond_{index}"
        existing_ids.add(connection.name)
        connection.type = gemmi.ConnectionType.Covale
        connection.asu = gemmi.Asu.Same
        connection.partner1, connection.partner2 = endpoints
        xyz1 = structure.atoms[link["atom_1"]]["coords"]
        xyz2 = structure.atoms[link["atom_2"]]["coords"]
        connection.reported_distance = float(
            sum((float(a) - float(b)) ** 2 for a, b in zip(xyz1, xyz2)) ** 0.5
        )
        model.connections.append(connection)
    if not model.connections:
        return text
    generated = model.make_mmcif_block(gemmi.MmcifOutputGroups(False, conn=True))
    for category in ("_struct_conn.", "_struct_conn_type."):
        old = block.get_mmcif_category(category, raw=True)
        new = generated.get_mmcif_category(category, raw=True)
        if not new:
            continue
        old_count = len(next(iter(old.values()), []))
        ids = set(old.get("id", []))
        keep = [i for i, name in enumerate(new["id"]) if name not in ids]
        new_count = len(new["id"])
        merged = {
            key: old.get(key, ["?"] * old_count)
            + [new.get(key, ["?"] * new_count)[i] for i in keep]
            for key in dict.fromkeys([*old, *new])
        }
        block.set_mmcif_category(category, merged, raw=True)
    return document.as_string()


def to_mmcif(
    structure: Structure,
    plddts: Optional[Tensor] = None,
    boltz2: bool = False,
) -> str:  # noqa: C901, PLR0915, PLR0912
    """Write a structure into an MMCIF file.

    Parameters
    ----------
    structure : Structure
        The input structure

    Returns
    -------
    str
        the output MMCIF file

    """
    system = System()

    # Load periodic table for element mapping
    periodic_table = Chem.GetPeriodicTable()

    # Map entities to chain_ids
    entity_to_chains = {}
    entity_to_moltype = {}

    for chain in structure.chains:
        entity_id = chain["entity_id"]
        mol_type = chain["mol_type"]
        entity_to_chains.setdefault(entity_id, []).append(chain)
        entity_to_moltype[entity_id] = mol_type

    # Map entities to sequences
    sequences = {}
    for entity in entity_to_chains:
        # Get the first chain
        chain = entity_to_chains[entity][0]

        # Get the sequence
        res_start = chain["res_idx"]
        res_end = chain["res_idx"] + chain["res_num"]
        residues = structure.residues[res_start:res_end]
        sequence = [str(res["name"]) for res in residues]
        sequences[entity] = sequence

    # Create entity objects
    lig_entity = None
    entities_map = {}
    for entity, sequence in sequences.items():
        mol_type = entity_to_moltype[entity]

        if mol_type == const.chain_type_ids["PROTEIN"]:
            alphabet = ihm.LPeptideAlphabet()
            chem_comp = lambda x: ihm.LPeptideChemComp(id=x, code=x, code_canonical="X")  # noqa: E731
        elif mol_type == const.chain_type_ids["DNA"]:
            alphabet = ihm.DNAAlphabet()
            chem_comp = lambda x: ihm.DNAChemComp(id=x, code=x, code_canonical="N")  # noqa: E731
        elif mol_type == const.chain_type_ids["RNA"]:
            alphabet = ihm.RNAAlphabet()
            chem_comp = lambda x: ihm.RNAChemComp(id=x, code=x, code_canonical="N")  # noqa: E731
        elif len(sequence) > 1:
            alphabet = {}
            chem_comp = lambda x: ihm.SaccharideChemComp(id=x)  # noqa: E731
        else:
            alphabet = {}
            chem_comp = lambda x: ihm.NonPolymerChemComp(id=x)  # noqa: E731

        # Handle smiles
        if len(sequence) == 1 and (sequence[0] == "LIG"):
            if lig_entity is None:
                seq = [chem_comp(sequence[0])]
                lig_entity = Entity(seq)
            model_e = lig_entity
        else:
            seq = [
                alphabet[item] if item in alphabet else chem_comp(item)
                for item in sequence
            ]
            model_e = Entity(seq)

        for chain in entity_to_chains[entity]:
            chain_idx = chain["asym_id"]
            entities_map[chain_idx] = model_e

    # We don't assume that symmetry is perfect, so we dump everything
    # into the asymmetric unit, and produce just a single assembly
    asym_unit_map = {}
    for chain in structure.chains:
        # Define the model assembly
        chain_idx = chain["asym_id"]
        chain_tag = str(chain["name"])
        entity = entities_map[chain_idx]
        if entity.type == "water":
            asym = ihm.WaterAsymUnit(
                entity,
                1,
                details="Model subunit %s" % chain_tag,
                id=chain_tag,
            )
        else:
            asym = AsymUnit(
                entity,
                details="Model subunit %s" % chain_tag,
                id=chain_tag,
            )
        asym_unit_map[chain_idx] = asym
    modeled_assembly = Assembly(asym_unit_map.values(), name="Modeled assembly")

    class _LocalPLDDT(modelcif.qa_metric.Local, modelcif.qa_metric.PLDDT):
        name = "pLDDT"
        software = None
        description = "Predicted lddt"

    class _MyModel(AbInitioModel):
        def get_atoms(self) -> Iterator[Atom]:
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
                het = chain["mol_type"] == const.chain_type_ids["NONPOLYMER"]
                chain_idx = chain["asym_id"]
                res_start = chain["res_idx"]
                res_end = chain["res_idx"] + chain["res_num"]

                record_type = (
                    "ATOM"
                    if chain["mol_type"] != const.chain_type_ids["NONPOLYMER"]
                    else "HETATM"
                )

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
                            element = periodic_table.GetElementSymbol(
                                atom["element"].item()
                            )
                        element = element.upper()
                        residue_index = residue["res_idx"] + 1
                        pos = atom_coords[i]

                        if record_type != "HETATM":
                            # The current residue plddt is stored at the res_num index unless a ligand has previouly been added.
                            biso = (
                                100.00
                                if plddts is None
                                else round(
                                    plddts[res_num + ligand_index_offset].item() * 100,
                                    3,
                                )
                            )
                            prev_polymer_resnum = res_num
                        else:
                            # If not a polymer resnum, we can get index into plddts by adding offset relative to previous polymer resnum.
                            ligand_index_offset += 1
                            biso = (
                                100.00
                                if plddts is None
                                else round(
                                    plddts[
                                        prev_polymer_resnum + ligand_index_offset
                                    ].item()
                                    * 100,
                                    3,
                                )
                            )

                        yield Atom(
                            asym_unit=asym_unit_map[chain_idx],
                            type_symbol=element,
                            seq_id=residue_index,
                            atom_id=atom_name,
                            x=f"{pos[0]:.5f}",
                            y=f"{pos[1]:.5f}",
                            z=f"{pos[2]:.5f}",
                            het=het,
                            biso=biso,
                            occupancy=1,
                        )

                    if record_type != "HETATM":
                        res_num += 1

        def add_plddt(self, plddts):
            res_num = 0
            prev_polymer_resnum = (
                -1
            )  # -1 handles case where ligand is the first residue
            ligand_index_offset = 0
            for chain in structure.chains:
                chain_idx = chain["asym_id"]
                res_start = chain["res_idx"]
                res_end = chain["res_idx"] + chain["res_num"]
                residues = structure.residues[res_start:res_end]

                record_type = (
                    "ATOM"
                    if chain["mol_type"] != const.chain_type_ids["NONPOLYMER"]
                    else "HETATM"
                )

                # We rename the chains in alphabetical order
                for residue in residues:
                    residue_idx = residue["res_idx"] + 1

                    atom_start = residue["atom_idx"]
                    atom_end = residue["atom_idx"] + residue["atom_num"]

                    if record_type != "HETATM":
                        # The current residue plddt is stored at the res_num index unless a ligand has previouly been added.
                        self.qa_metrics.append(
                            _LocalPLDDT(
                                asym_unit_map[chain_idx].residue(residue_idx),
                                round(
                                    plddts[res_num + ligand_index_offset].item() * 100,
                                    3,
                                ),
                            )
                        )
                        prev_polymer_resnum = res_num
                    else:
                        # If not a polymer resnum, we can get index into plddts by adding offset relative to previous polymer resnum.
                        self.qa_metrics.append(
                            _LocalPLDDT(
                                asym_unit_map[chain_idx].residue(residue_idx),
                                round(
                                    plddts[
                                        prev_polymer_resnum
                                        + ligand_index_offset
                                        + 1 : prev_polymer_resnum
                                        + ligand_index_offset
                                        + residue["atom_num"]
                                        + 1
                                    ]
                                    .mean()
                                    .item()
                                    * 100,
                                    2,
                                ),
                            )
                        )
                        ligand_index_offset += residue["atom_num"]

                    if record_type != "HETATM":
                        res_num += 1

    # Add the model and modeling protocol to the file and write them out:
    model = _MyModel(assembly=modeled_assembly, name="Model")
    if plddts is not None:
        model.add_plddt(plddts)

    model_group = ModelGroup([model], name="All models")
    system.model_groups.append(model_group)
    ihm.dumper.set_line_wrap(False)

    fh = io.StringIO()
    dumper.write(fh, [system])
    return _add_covalent_connections(fh.getvalue(), structure, boltz2)
