"""Create a mapping from structure and chain ID to MSA indices."""

import argparse
import hashlib
import json
import pickle
import subprocess
from pathlib import Path

import pandas as pd
from Bio import SeqIO


def hash_sequence(seq: str) -> str:
    """Hash a sequence."""
    return hashlib.sha256(seq.encode()).hexdigest()


def main(args: argparse.Namespace) -> None:
    """Create clustering."""
    # Set output directory
    output = Path(args.output)

    # Split the sequences into proteins and nucleotides
    with Path(args.sequences).open("r") as handle:
        data = list(SeqIO.parse(handle, "fasta"))

    proteins = []
    nucleotides = []

    for seq in data:
        if set(seq.seq).issubset({"A", "C", "G", "T", "U"}):
            nucleotides.append(seq)
        else:
            proteins.append(seq)

    # Run mmseqs on the protein data
    proteins = [f">{hash_sequence(seq)}\n{seq}" for seq in proteins]
    with (output / "proteins.fasta").open("w") as handle:
        SeqIO.write(proteins, handle, "fasta")

    subprocess.run(
        f"{args.mmseqs} easy-cluster {output / "proteins.fasta"} {output / 'clust_prot'} {output / 'tmp'} --min-seq-id 0.4",  # noqa: E501
        shell=True,  # noqa: S602
        check=True,
    )

    # Load protein clusters
    clustering_path = output / "clust_prot_cluster.tsv"
    protein_data = pd.read_csv(clustering_path, sep="\t", header=None)
    clusters = protein_data[0]
    items = protein_data[1]
    clustering = dict(zip(list(items), list(clusters)))

    # Each unique rna sequence is given an id
    visited = {}
    for nucl in nucleotides:
        nucl_id = hash_sequence(nucl)
        if nucl not in visited:
            clustering[nucl_id] = nucl_id
            visited[nucl] = nucl_id
        else:
            clustering[nucl_id] = visited[nucl]

    # Load ligand data
    with Path(args.ccd).open("rb") as handle:
        ligand_data = pickle.load(handle)  # noqa: S301

    # Each unique ligand CCD is given an id
    visited = {}
    for ccd_code in ligand_data:
        clustering[ccd_code] = ccd_code

    # Save clustering
    with (output / "clustering.json").open("w") as handle:
        json.dump(clustering, handle)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sequences", type=str, help="Input to protein fasta.")
    parser.add_argument("--ccd", type=str, help="Input to rna fasta.")
    parser.add_argument("--output", type=str, help="Output clustering.")
    parser.add_argument(
        "--mmseqs",
        type=str,
        help="Path to mmseqs program.",
        default="mmseqs",
    )
    args = parser.parse_args()
    main(args)
