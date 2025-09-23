from collections.abc import Mapping
from pathlib import Path

from Bio import SeqIO
from rdkit.Chem.rdchem import Mol

from boltz.data.parse.yaml import parse_boltz_schema
from boltz.data.types import Target


def parse_fasta(  # noqa: C901, PLR0912
    path: Path,
    ccd: Mapping[str, Mol],
    mol_dir: Path,
    boltz2: bool = False,
) -> Target:
    """Parse a fasta file.

    The name of the fasta file is used as the name of this job.
    We rely on the fasta record id to determine the entity type.

    > CHAIN_ID|ENTITY_TYPE|MSA_ID
    SEQUENCE
    > CHAIN_ID|ENTITY_TYPE|MSA_ID
    ...

    Where ENTITY_TYPE is either protein, rna, dna, ccd or smiles,
    and CHAIN_ID is the chain identifier, which should be unique.
    The MSA_ID is optional and should only be used on proteins.

    Parameters
    ----------
    fasta_file : Path
        Path to the fasta file.
    ccd : Dict
        Dictionary of CCD components.
    mol_dir : Path
        Path to the directory containing the molecules.
    boltz2 : bool
        Whether to parse the input for Boltz2.

    Returns
    -------
    Target
        The parsed target.

    """
    # Read fasta file
    with path.open("r") as f:
        records = list(SeqIO.parse(f, "fasta"))

    # Make sure all records have a chain id and entity
    for seq_record in records:
        if "|" not in seq_record.id:
            msg = f"Invalid record id: {seq_record.id}"
            raise ValueError(msg)

        header = seq_record.id.split("|")
        assert len(header) >= 2, f"Invalid record id: {seq_record.id}"

        chain_id, entity_type = header[:2]
        if entity_type.lower() not in {"protein", "dna", "rna", "ccd", "smiles"}:
            msg = f"Invalid entity type: {entity_type}"
            raise ValueError(msg)
        if chain_id == "":
            msg = "Empty chain id in input fasta!"
            raise ValueError(msg)
        if entity_type == "":
            msg = "Empty entity type in input fasta!"
            raise ValueError(msg)

    # Convert to yaml format
    sequences = []
    for seq_record in records:
        # Get chain id, entity type and sequence
        header = seq_record.id.split("|")
        chain_id, entity_type = header[:2]
        if len(header) == 3 and header[2] != "":
            assert entity_type.lower() == "protein", (
                "MSA_ID is only allowed for proteins"
            )
            msa_id = header[2]
        else:
            msa_id = None

        entity_type = entity_type.upper()
        seq = str(seq_record.seq)

        if entity_type == "PROTEIN":
            molecule = {
                "protein": {
                    "id": chain_id,
                    "sequence": seq,
                    "modifications": [],
                    "msa": msa_id,
                },
            }
        elif entity_type == "RNA":
            molecule = {
                "rna": {
                    "id": chain_id,
                    "sequence": seq,
                    "modifications": [],
                },
            }
        elif entity_type == "DNA":
            molecule = {
                "dna": {
                    "id": chain_id,
                    "sequence": seq,
                    "modifications": [],
                }
            }
        elif entity_type.upper() == "CCD":
            molecule = {
                "ligand": {
                    "id": chain_id,
                    "ccd": seq,
                }
            }
        elif entity_type.upper() == "SMILES":
            molecule = {
                "ligand": {
                    "id": chain_id,
                    "smiles": seq,
                }
            }

        sequences.append(molecule)

    data = {
        "sequences": sequences,
        "bonds": [],
        "version": 1,
    }

    name = path.stem
    return parse_boltz_schema(name, data, ccd, mol_dir, boltz2)
