# SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
# SPDX-License-Identifier: MIT
#
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

"""Tests for scripts/process/cluster.py OS_CMD_INJECTION fix.

Verifies that main() uses subprocess argument lists (no shell=True) and that
the full clustering pipeline (FASTA parsing, mmseqs call, short/nucleotide/
ligand handling, JSON output) is preserved.
"""

import argparse
import hashlib
import json
import pickle
import sys
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts" / "process"))
import cluster as cluster_module  # noqa: E402


def _hash(seq: str) -> str:
    return hashlib.sha256(seq.encode()).hexdigest()


def _write_fasta(path: Path, sequences: list[tuple[str, str]]) -> None:
    """Write a minimal FASTA file.  Each item is (header, sequence)."""
    lines = []
    for header, seq in sequences:
        lines.append(f">{header}")
        lines.append(seq)
    path.write_text("\n".join(lines))


def _write_cluster_tsv(path: Path, mapping: list[tuple[str, str]]) -> None:
    """Write a fake mmseqs clust_prot_cluster.tsv (tab-separated, no header)."""
    lines = [f"{rep}\t{member}" for rep, member in mapping]
    path.write_text("\n".join(lines))


def _write_ccd_pickle(path: Path, ligand_codes: list[str]) -> None:
    """Write a minimal CCD pickle (a dict keyed by ligand code)."""
    data = dict.fromkeys(ligand_codes)
    with path.open("wb") as f:
        pickle.dump(data, f)


class TestMainCallsMmseqs:
    """Verify subprocess.run invocation and full pipeline output."""

    def test_calls_mmseqs_with_correct_args(self, tmp_path):
        seq_a = "MKTAYIAKQRQISFVKSHFSRQ"
        seq_b = "MLLSALVLLLSESGLSGAGGL"
        fasta_path = tmp_path / "input.fasta"
        _write_fasta(fasta_path, [("seqA", seq_a), ("seqB", seq_b)])

        ccd_path = tmp_path / "ccd.pkl"
        _write_ccd_pickle(ccd_path, ["ATP", "NAG"])

        outdir = tmp_path / "out"

        hash_a = _hash(seq_a)
        hash_b = _hash(seq_b)

        args = argparse.Namespace(
            sequences=str(fasta_path),
            ccd=str(ccd_path),
            outdir=str(outdir),
            mmseqs="/usr/bin/mmseqs",
        )

        def _fake_run(cmd, **kwargs):
            outdir.mkdir(parents=True, exist_ok=True)
            _write_cluster_tsv(
                outdir / "clust_prot_cluster.tsv",
                [(hash_a, hash_a), (hash_a, hash_b)],
            )

        with patch("cluster.subprocess.run", side_effect=_fake_run) as mock_run:
            cluster_module.main(args)

        mock_run.assert_called_once()
        call_args, call_kwargs = mock_run.call_args

        cmd = call_args[0]
        assert isinstance(cmd, list), f"Expected list, got {type(cmd)}"
        assert cmd[0] == "/usr/bin/mmseqs"
        assert cmd[1] == "easy-cluster"
        assert cmd[-2:] == ["--min-seq-id", "0.4"]
        assert "shell" not in call_kwargs, "shell=True must not be passed"
        assert call_kwargs.get("check") is True

        clustering_file = outdir / "clustering.json"
        assert clustering_file.exists()
        with clustering_file.open() as f:
            clustering = json.load(f)

        assert hash_a in clustering
        assert hash_b in clustering
        assert clustering[hash_a] == hash_a
        assert clustering[hash_b] == hash_a
        assert "ATP" in clustering
        assert "NAG" in clustering
        assert clustering["ATP"] == "ATP"

    def test_handles_short_and_nucleotide_sequences(self, tmp_path):
        protein = "MKTAYIAKQRQISFVKSHFSRQ"
        short = "MKTAY"
        nucleotide = "ACGUACGU"

        fasta_path = tmp_path / "input.fasta"
        _write_fasta(
            fasta_path,
            [
                ("prot1", protein),
                ("short1", short),
                ("nucl1", nucleotide),
            ],
        )

        ccd_path = tmp_path / "ccd.pkl"
        _write_ccd_pickle(ccd_path, [])

        outdir = tmp_path / "out"
        hash_prot = _hash(protein)

        args = argparse.Namespace(
            sequences=str(fasta_path),
            ccd=str(ccd_path),
            outdir=str(outdir),
            mmseqs="mmseqs",
        )

        def _fake_run(cmd, **kwargs):
            outdir.mkdir(parents=True, exist_ok=True)
            _write_cluster_tsv(
                outdir / "clust_prot_cluster.tsv",
                [(hash_prot, hash_prot)],
            )

        with patch("cluster.subprocess.run", side_effect=_fake_run) as mock_run:
            cluster_module.main(args)

        mock_run.assert_called_once()

        with (outdir / "clustering.json").open() as f:
            clustering = json.load(f)

        hash_short = _hash(short)
        hash_nucl = _hash(nucleotide)
        assert clustering[hash_short] == hash_short, "Short sequence should get self-referential ID"
        assert clustering[hash_nucl] == hash_nucl, "Nucleotide sequence should get self-referential ID"
        assert hash_prot in clustering, "Protein should be in clustering"
