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


"""Tests for scripts/eval/run_evals.py OS_CMD_INJECTION fix.

Verifies that evaluate_structure uses subprocess argument lists (no shell=True)
and preserves the original docker command semantics.
"""

import sys
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts" / "eval"))
import run_evals  # noqa: E402


class TestEvaluateStructureCalls:
    """Verify subprocess invocations produced by evaluate_structure."""

    def test_calls_docker_for_structure_and_ligand(self, tmp_path):
        outdir = tmp_path / "results"
        outdir.mkdir()

        with (
            patch("run_evals.subprocess.run") as mock_run,
            patch("run_evals.os.getuid", return_value=1000),
            patch("run_evals.os.getgid", return_value=2000),
        ):
            run_evals.evaluate_structure(
                name="test_sample",
                pred="/data/pred.cif",
                reference="/data/ref.cif",
                outdir=str(outdir),
                mount="/data",
            )

        assert mock_run.call_count == 2, f"Expected 2 subprocess.run calls, got {mock_run.call_count}"

        # --- Structure comparison call ---
        struct_args, struct_kwargs = mock_run.call_args_list[0]
        struct_cmd = struct_args[0]

        assert "shell" not in struct_kwargs, "shell=True must not be passed"
        assert struct_kwargs.get("check") is False
        assert struct_kwargs.get("capture_output") is True

        assert struct_cmd[:3] == ["sudo", "docker", "run"]
        assert "-u" in struct_cmd
        u_idx = struct_cmd.index("-u")
        assert struct_cmd[u_idx + 1] == "1000:2000"
        assert "--rm" in struct_cmd
        assert "--volume" in struct_cmd
        vol_idx = struct_cmd.index("--volume")
        assert struct_cmd[vol_idx + 1] == "/data:/data"
        assert run_evals.IMAGE_NAME in struct_cmd
        assert "compare-structures" in struct_cmd

        for flag in [
            "--lddt",
            "--bb-lddt",
            "--qs-score",
            "--dockq",
            "--ics",
            "--ips",
            "--rigid-scores",
            "--patch-scores",
            "--tm-score",
            "--fault-tolerant",
        ]:
            assert flag in struct_cmd, f"Missing flag {flag} in structure command"
        assert "-m" in struct_cmd
        m_idx = struct_cmd.index("-m")
        assert struct_cmd[m_idx + 1] == "/data/pred.cif"
        assert "-r" in struct_cmd
        r_idx = struct_cmd.index("-r")
        assert struct_cmd[r_idx + 1] == "/data/ref.cif"
        pep_idx = struct_cmd.index("--min-pep-length")
        assert struct_cmd[pep_idx + 1] == "4"
        nuc_idx = struct_cmd.index("--min-nuc-length")
        assert struct_cmd[nuc_idx + 1] == "4"

        # --- Ligand comparison call ---
        ligand_args, ligand_kwargs = mock_run.call_args_list[1]
        ligand_cmd = ligand_args[0]

        assert "shell" not in ligand_kwargs, "shell=True must not be passed"
        assert "compare-ligand-structures" in ligand_cmd

        for flag in ["--lddt-pli", "--rmsd", "--substructure-match", "--fault-tolerant"]:
            assert flag in ligand_cmd, f"Missing flag {flag} in ligand command"

        expected_ligand_out = str(outdir / "test_sample_ligand.json")
        o_idx = ligand_cmd.index("-o")
        assert ligand_cmd[o_idx + 1] == expected_ligand_out

    def test_skips_existing_outputs(self, tmp_path):
        outdir = tmp_path / "results"
        outdir.mkdir()
        (outdir / "test_sample.json").touch()
        (outdir / "test_sample_ligand.json").touch()

        with patch("run_evals.subprocess.run") as mock_run:
            run_evals.evaluate_structure(
                name="test_sample",
                pred="/data/pred.cif",
                reference="/data/ref.cif",
                outdir=str(outdir),
                mount="/data",
            )

        mock_run.assert_not_called()

    def test_skips_only_structure_when_structure_exists(self, tmp_path):
        outdir = tmp_path / "results"
        outdir.mkdir()
        (outdir / "test_sample.json").touch()

        with (
            patch("run_evals.subprocess.run") as mock_run,
            patch("run_evals.os.getuid", return_value=0),
            patch("run_evals.os.getgid", return_value=0),
        ):
            run_evals.evaluate_structure(
                name="test_sample",
                pred="/data/pred.cif",
                reference="/data/ref.cif",
                outdir=str(outdir),
                mount="/data",
            )

        assert mock_run.call_count == 1
        cmd = mock_run.call_args_list[0][0][0]
        assert "compare-ligand-structures" in cmd

    def test_uid_gid_values(self, tmp_path):
        outdir = tmp_path / "results"
        outdir.mkdir()

        with (
            patch("run_evals.subprocess.run") as mock_run,
            patch("run_evals.os.getuid", return_value=12345),
            patch("run_evals.os.getgid", return_value=67890),
        ):
            run_evals.evaluate_structure(
                name="sample",
                pred="/mnt/pred.cif",
                reference="/mnt/ref.cif",
                outdir=str(outdir),
                mount="/mnt",
            )

        for c in mock_run.call_args_list:
            cmd = c[0][0]
            u_idx = cmd.index("-u")
            assert cmd[u_idx + 1] == "12345:67890"
