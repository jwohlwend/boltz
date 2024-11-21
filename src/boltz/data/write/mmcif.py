import io
from typing import Iterator

import ihm
from modelcif import Assembly, AsymUnit, Entity, System, dumper
from modelcif.model import AbInitioModel, Atom, ModelGroup
from rdkit import Chem

from boltz.data import const
from boltz.data.types import Structure
from boltz.data.write.utils import generate_tags


def to_mmcif(structure: Structure) -> str:  # noqa: C901
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
    for k,(entity, sequence) in enumerate(sequences.items()):
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
            chem_comp = lambda x: ihm.SaccharideChemComp(id=x,code_canonical=f"X{k}")  # noqa: E731
        else:
            alphabet = {}
            chem_comp = lambda x: ihm.NonPolymerChemComp(id=x,code_canonical=f"X{k}")  # noqa: E731
        seq = [
            alphabet[item] if item in alphabet else chem_comp(item) for item in sequence
        ]
        model_e = Entity(seq)
        for chain in entity_to_chains[entity]:
            chain_idx = chain["asym_id"]
            entities_map[chain_idx] = model_e

    # We don't assume that symmetry is perfect, so we dump everything
    # into the asymmetric unit, and produce just a single assembly
    chain_tags = generate_tags()
    asym_unit_map = {}
    for chain in structure.chains:
        # Define the model assembly
        chain_idx = chain["asym_id"]
        chain_tag = next(chain_tags)
        asym = AsymUnit(
            entities_map[chain_idx],
            details="Model subunit %s" % chain_tag,
            id=chain_tag,
        )
        asym_unit_map[chain_idx] = asym
    modeled_assembly = Assembly(asym_unit_map.values(), name="Modeled assembly")

    # class _LocalPLDDT(modelcif.qa_metric.Local, modelcif.qa_metric.PLDDT):
    #     name = "pLDDT"
    #     software = None
    #     description = "Predicted lddt"

    # class _GlobalPLDDT(modelcif.qa_metric.Global, modelcif.qa_metric.PLDDT):
    #     name = "pLDDT"
    #     software = None
    #     description = "Global pLDDT, mean of per-residue pLDDTs"

    class _MyModel(AbInitioModel):
        def get_atoms(self) -> Iterator[Atom]:
            # Add all atom sites.
            for chain in structure.chains:
                # We rename the chains in alphabetical order
                het = chain["mol_type"] == const.chain_type_ids["NONPOLYMER"]
                chain_idx = chain["asym_id"]
                res_start = chain["res_idx"]
                res_end = chain["res_idx"] + chain["res_num"]

                residues = structure.residues[res_start:res_end]
                for residue in residues:
                    atom_start = residue["atom_idx"]
                    atom_end = residue["atom_idx"] + residue["atom_num"]
                    atoms = structure.atoms[atom_start:atom_end]
                    atom_coords = atoms["coords"]
                    for i, atom in enumerate(atoms):
                        # This should not happen on predictions, but just in case.
                        if not atom["is_present"]:
                            continue

                        name = atom["name"]
                        name = [chr(c + 32) for c in name if c != 0]
                        name = "".join(name)
                        element = periodic_table.GetElementSymbol(
                            atom["element"].item()
                        )
                        element = element.upper()
                        residue_index = residue["res_idx"] + 1
                        pos = atom_coords[i]
                        yield Atom(
                            asym_unit=asym_unit_map[chain_idx],
                            type_symbol=element,
                            seq_id=residue_index,
                            atom_id=name,
                            x=pos[0],
                            y=pos[1],
                            z=pos[2],
                            het=het,
                            biso=1.00,
                            occupancy=1.00,
                        )

        def add_scores(self):
            return
            # local scores
            # plddt_per_residue = {}
            # for i in range(n):
            #     for mask, b_factor in zip(atom_mask[i], b_factors[i]):
            #         if mask < 0.5:
            #             continue
            #         # add 1 per residue, not 1 per atom
            #         if chain_index[i] not in plddt_per_residue:
            #             # first time a chain index is seen: add the key and start the residue dict
            #             plddt_per_residue[chain_index[i]] = {residue_index[i]: b_factor}
            #         if residue_index[i] not in plddt_per_residue[chain_index[i]]:
            #             plddt_per_residue[chain_index[i]][residue_index[i]] = b_factor
            # plddts = []
            # for chain_idx in plddt_per_residue:
            #     for residue_idx in plddt_per_residue[chain_idx]:
            #         plddt = plddt_per_residue[chain_idx][residue_idx]
            #         plddts.append(plddt)
            #         self.qa_metrics.append(
            #             _LocalPLDDT(
            #                 asym_unit_map[chain_idx].residue(residue_idx), plddt
            #             )
            #         )
            # # global score
            # self.qa_metrics.append((_GlobalPLDDT(np.mean(plddts))))

    # Add the model and modeling protocol to the file and write them out:
    model = _MyModel(assembly=modeled_assembly, name="Model")
    # model.add_scores()

    model_group = ModelGroup([model], name="All models")
    system.model_groups.append(model_group)

    fh = io.StringIO()
    dumper.write(fh, [system])
    return fh.getvalue()
