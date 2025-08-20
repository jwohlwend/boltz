#!/usr/bin/env python3
"""
Generate Boltz ligand YAML configs from a CSV file.

For each row in the input CSV, this script writes a YAML file with the form:

version: 1
sequences:
  - ligand:
      id: A
      smiles: "<SMILES_FROM_ROW>"

By default, it searches for a SMILES column by common names (case-insensitive):
  smiles, smile, smiles_string, smile_string

If an identifier column is present (e.g., id, molecule_id, ligand_id, sample_id,
index, name), it is used to name the output file as ligand_<ID>.yaml. Otherwise,
the zero-based row index is used.

Usage:
  python scripts/process/generate_ligand_yamls.py /path/to/train.csv \
    --out /path/to/output/dir \
    --yaml-id A

This script does not require any external dependencies beyond the Python standard library.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path
from typing import Iterable, Optional


DEFAULT_SMILES_CANDIDATES = (
    "smiles",
    "smile",
    "smiles_string",
    "smile_string",
)

DEFAULT_ID_CANDIDATES = (
    "id",
    "molecule_id",
    "ligand_id",
    "polymer_id",
    "sample_id",
    "index",
    "name",
)


def find_header(headers: Iterable[str], candidates: Iterable[str]) -> Optional[str]:
    """Return the first header whose lowercase matches one of the candidates."""
    lower_to_original = {h.lower(): h for h in headers if h}
    for candidate in candidates:
        original = lower_to_original.get(candidate.lower())
        if original:
            return original
    return None


def sanitize_for_filename(value: str) -> str:
    """Sanitize an arbitrary string for safe use in filenames."""
    if value is None:
        return ""
    safe = re.sub(r"[^A-Za-z0-9._-]", "_", str(value))
    return safe or "row"


def write_yaml_file(output_path: Path, yaml_id: str, smiles: str) -> None:
    """Write a minimal ligand YAML with JSON-quoted SMILES for YAML compatibility."""
    # Use json.dumps to ensure proper escaping inside double quotes. YAML 1.2 is a superset of JSON.
    smiles_quoted = json.dumps(smiles)
    content = (
        "version: 1\n"
        "sequences:\n"
        "  - ligand:\n"
        f"      id: {yaml_id}\n"
        f"      smiles: {smiles_quoted}\n"
    )
    output_path.write_text(content, encoding="utf-8")


def generate_ligand_yamls(
    csv_path: Path,
    output_dir: Path,
    smiles_column: Optional[str] = None,
    id_column: Optional[str] = None,
    yaml_id: str = "A",
    filename_from_row_index: bool = False,
    row_start_index: int = 0,
) -> int:
    """Read CSV and write one YAML per row. Returns the number of files written."""
    output_dir.mkdir(parents=True, exist_ok=True)

    with csv_path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        headers = reader.fieldnames or []

        smiles_col = smiles_column or find_header(headers, DEFAULT_SMILES_CANDIDATES)
        if not smiles_col:
            raise SystemExit(
                "Could not find a SMILES column. Available columns: "
                + ", ".join(headers)
            )

        row_id_col = None if filename_from_row_index else (id_column or find_header(headers, DEFAULT_ID_CANDIDATES))

        written = 0
        for row_index, row in enumerate(reader):
            smiles_value = (row.get(smiles_col) or "").strip()
            if not smiles_value:
                continue

            if filename_from_row_index:
                file_index = row_start_index + row_index
                filename_suffix = str(file_index)
            else:
                raw_id_value = row.get(row_id_col) if row_id_col else str(row_index)
                if raw_id_value is None or str(raw_id_value).strip() == "":
                    raw_id_value = str(row_index)
                filename_suffix = sanitize_for_filename(str(raw_id_value))

            out_file = output_dir / f"ligand_{filename_suffix}.yaml"
            write_yaml_file(out_file, yaml_id=yaml_id, smiles=smiles_value)
            written += 1

    return written


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate ligand YAMLs from CSV")
    parser.add_argument(
        "csv",
        type=Path,
        help="Path to the input CSV file",
    )
    parser.add_argument(
        "--out",
        dest="out_dir",
        type=Path,
        default=Path("examples/generated_ligands"),
        help="Output directory for generated YAML files (default: examples/generated_ligands)",
    )
    parser.add_argument(
        "--smiles-col",
        dest="smiles_col",
        type=str,
        default=None,
        help="Explicit SMILES column name (case-sensitive). If omitted, a common name is auto-detected.",
    )
    parser.add_argument(
        "--id-col",
        dest="id_col",
        type=str,
        default=None,
        help="Explicit ID column name used for naming files. If omitted, a common name is auto-detected.",
    )
    parser.add_argument(
        "--yaml-id",
        dest="yaml_id",
        type=str,
        default="A",
        help="Value to put in sequences[0].ligand.id in each YAML (default: A)",
    )
    parser.add_argument(
        "--filename-source",
        choices=("row-index", "column"),
        default="row-index",
        help="Use row index or a CSV column to name files (default: row-index)",
    )
    parser.add_argument(
        "--row-start",
        dest="row_start",
        type=int,
        default=0,
        help="Starting index when using row-index for filenames (default: 0)",
    )
    parser.add_argument(
        "--clean",
        action="store_true",
        help="If set, removes existing .yaml files in the output directory before generation",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    # Optional cleanup
    if args.clean:
        args.out_dir.mkdir(parents=True, exist_ok=True)
        for p in args.out_dir.glob("*.yaml"):
            try:
                p.unlink()
            except Exception:
                pass
    num_written = generate_ligand_yamls(
        csv_path=args.csv,
        output_dir=args.out_dir,
        smiles_column=args.smiles_col,
        id_column=args.id_col,
        yaml_id=args.yaml_id,
        filename_from_row_index=(args.filename_source == "row-index"),
        row_start_index=args.row_start,
    )
    print(f"Wrote {num_written} YAML files to {args.out_dir}")


if __name__ == "__main__":
    main()


