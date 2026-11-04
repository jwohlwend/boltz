from tempfile import NamedTemporaryFile
from typing import Optional

import gemmi
from rdkit.Chem.rdchem import Mol

from boltz.data.parse.mmcif import ParsedStructure, parse_mmcif


def parse_pdb(
    path: str,
    mols: Optional[dict[str, Mol]] = None,
    moldir: Optional[str] = None,
    use_assembly: bool = True,
    compute_interfaces: bool = True,
) -> ParsedStructure:
    with NamedTemporaryFile(suffix=".cif") as tmp_cif_file:
        tmp_cif_path = tmp_cif_file.name
        structure = gemmi.read_structure(str(path))
        structure.setup_entities()

        subchain_counts, subchain_renaming = {}, {}
        for model in structure:
            for chain in model:
                subchain_counts[chain.name] = 0
                for res in chain:
                    if res.subchain not in subchain_renaming:
                        subchain_renaming[res.subchain] = chain.name + str(
                            subchain_counts[chain.name] + 1
                        )
                        subchain_counts[chain.name] += 1
                    res.subchain = subchain_renaming[res.subchain]
        for entity in structure.entities:
            entity.subchains = [
                subchain_renaming[subchain] for subchain in entity.subchains
            ]

        doc = structure.make_mmcif_document()
        doc.write_file(tmp_cif_path)

        return parse_mmcif(
            path=tmp_cif_path,
            mols=mols,
            moldir=moldir,
            use_assembly=use_assembly,
            compute_interfaces=compute_interfaces,
        )
