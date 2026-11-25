"""Regression tests for CSV MSA sequence parsing without model inference."""

from pathlib import Path

import numpy as np
import pytest

from boltz.data.parse.a3m import parse_a3m
from boltz.data.parse.csv import parse_csv


@pytest.mark.parametrize("sequence", ["NA", "NAN", "NaN", "nan", "NULL"])
def test_na_like_protein_sequence(tmp_path: Path, sequence: str) -> None:
    """Preserve valid residues that pandas otherwise interprets as missing."""
    csv_path = tmp_path / "msa.csv"
    csv_path.write_text(f"sequence,key\n{sequence},\n")
    a3m_path = tmp_path / "msa.a3m"
    a3m_path.write_text(f">query\n{sequence}\n")

    actual = parse_csv(csv_path)
    expected = parse_a3m(a3m_path, taxonomy=None)
    np.testing.assert_array_equal(actual.residues, expected.residues)
    np.testing.assert_array_equal(actual.deletions, expected.deletions)
    np.testing.assert_array_equal(actual.sequences, expected.sequences)


def test_empty_sequence_cells(tmp_path: Path) -> None:
    """Skip empty cells as well as whitespace-only sequences."""
    path = tmp_path / "msa.csv"
    path.write_text('sequence,key\n,123\n"",456\n   ,789\nACD,\n')
    actual = parse_csv(path)
    assert len(actual.sequences) == 1
    assert len(actual.residues) == 3


def test_numeric_and_missing_keys(tmp_path: Path) -> None:
    """Keep taxonomy-key inference and missing-key defaults unchanged."""
    path = tmp_path / "msa.csv"
    path.write_text("sequence,key\nNA,123\nACD,\nEFG,NA\nHIK,-1\n")
    actual = parse_csv(path)
    assert actual.sequences["taxonomy"].tolist() == [123, -1, -1, -1]


def test_duplicate_sequences_and_limit(tmp_path: Path) -> None:
    """Retain the existing deduplication and sequence-limit behavior."""
    path = tmp_path / "msa.csv"
    path.write_text("sequence,key\nACD,\nACD,\nEFG,\nHIK,\n")
    actual = parse_csv(path, max_seqs=2)
    assert len(actual.sequences) == 2
    assert len(actual.residues) == 6


def test_invalid_columns(tmp_path: Path) -> None:
    """Reject CSVs that do not have the required column names."""
    path = tmp_path / "msa.csv"
    path.write_text("residues,key\nACD,\n")
    with pytest.raises(ValueError, match="Invalid CSV format"):
        parse_csv(path)
