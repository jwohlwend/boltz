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
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Split the sequences into proteins and nucleotides
    with Path(args.sequences).open("r") as f:
        data = list(SeqIO.parse(f, "fasta"))

    proteins = set()
    shorts = set()
    nucleotides = set()

    # Separate the sequences into proteins, nucleotides and short sequences
    # Short sequences cause a bug in the clustering, so they are separated
    for seq in data:
        if set(str(seq.seq)).issubset({"A", "C", "G", "T", "U", "N"}):
            nucleotides.add(str(seq.seq).strip())
        elif len(str(seq.seq).strip()) < 10:  # noqa: PLR2004
            shorts.add(str(seq.seq).strip())
        else:
            proteins.add(str(seq.seq).strip())

    # Run mmseqs on the protein data
    proteins = [f">{hash_sequence(seq)}\n{seq}" for seq in proteins]
    with (outdir / "proteins.fasta").open("w") as f:
        f.write("\n".join(proteins))

    subprocess.run(
        f"{args.mmseqs} easy-cluster {outdir / 'proteins.fasta'} {outdir / 'clust_prot'} {outdir / 'tmp'} --min-seq-id 0.4",  # noqa: E501
        shell=True,  # noqa: S602
        check=True,
    )

    # Load protein clusters
    clustering_path = outdir / "clust_prot_cluster.tsv"
    protein_data = pd.read_csv(clustering_path, sep="\t", header=None)
    clusters = protein_data[0]
    items = protein_data[1]
    clustering = dict(zip(list(items), list(clusters)))

    # Each shqrt sequence is given an id
    for short in shorts:
        short_id = hash_sequence(short)
        clustering[short_id] = short_id

    # Each unique rna sequence is given an id
    for nucl in nucleotides:
        nucl_id = hash_sequence(nucl)
        clustering[nucl_id] = nucl_id

    # Load ligand data
    with Path(args.ccd).open("rb") as handle:
        ligand_data = pickle.load(handle)  # noqa: S301

    # Each unique ligand CCD is given an id
    for ccd_code in ligand_data:
        clustering[ccd_code] = ccd_code

    # Save clustering
    with (outdir / "clustering.json").open("w") as handle:
        json.dump(clustering, handle)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--sequences",
        type=str,
        help="Input to protein fasta.",
        required=True,
    )
    parser.add_argument(
        "--ccd",
        type=str,
        help="Input to rna fasta.",
        required=True,
    )
    parser.add_argument(
        "--outdir",
        type=str,
        help="Output directory.",
        required=True,
    )
    parser.add_argument(
        "--mmseqs",
        type=str,
        help="Path to mmseqs program.",
        default="mmseqs",
    )
    args = parser.parse_args()
    main(args)
