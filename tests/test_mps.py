"""MPS (Apple Silicon) tests: verify boltz works on Metal Performance Shaders.

These tests require an Apple Silicon Mac with MPS available.
They are excluded from CI and must be run manually:

    pytest tests/test_mps.py -m mps -v
"""

import json
import os
import subprocess
import sys
import tempfile

import pytest

mps_available = False
try:
    import torch

    mps_available = torch.backends.mps.is_available()
except Exception:
    pass

pytestmark = [
    pytest.mark.mps,
    pytest.mark.skipif(not mps_available, reason="MPS not available"),
]


def _run_boltz_predict(input_yaml, input_filename, tmpdir, extra_args=None):
    """Run boltz predict with --accelerator mps and return (result, predictions_dir)."""
    input_path = os.path.join(tmpdir, input_filename)
    with open(input_path, "w") as f:
        f.write(input_yaml)

    output_dir = os.path.join(tmpdir, "output")
    cmd = [
        sys.executable, "-c", "from boltz.main import cli; cli()", "predict", input_path,
        "--out_dir", output_dir,
        "--accelerator", "mps",
        "--recycling_steps", "1",
        "--diffusion_samples", "1",
    ]
    if extra_args:
        cmd.extend(extra_args)

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=1200,
    )

    stem = os.path.splitext(input_filename)[0]
    pred_dir = os.path.join(output_dir, f"boltz_results_{stem}", "predictions")
    return result, pred_dir


def _find_files(pred_dir, extension):
    """Find all files with given extension in prediction subdirectories."""
    found = []
    if not os.path.isdir(pred_dir):
        return found
    for subdir in os.listdir(pred_dir):
        subdir_path = os.path.join(pred_dir, subdir)
        if os.path.isdir(subdir_path):
            for f in os.listdir(subdir_path):
                if f.endswith(extension):
                    found.append(os.path.join(subdir_path, f))
    return found


def test_autocast_device_type_mps():
    """autocast_device_type returns a valid device type string for MPS."""
    from boltz.model.modules.utils import autocast_device_type

    result = autocast_device_type("mps")
    # Must be accepted by torch.autocast without raising
    with torch.autocast(result, enabled=False):
        pass


def test_checkpoint_map_location_falls_back_to_cpu_for_mps():
    """MPS checkpoint loads should avoid the broken torch.mps.current_device path."""
    from boltz.main import _resolve_checkpoint_load_map_location

    result = _resolve_checkpoint_load_map_location("mps")
    assert result == "cpu"


def test_predict_peptide_mps():
    """Run boltz predict on MPS with a small peptide."""
    input_yaml = """\
version: 1
sequences:
  - protein:
      id: A
      sequence: ACDEFGHIKL
      msa: empty
"""
    with tempfile.TemporaryDirectory() as tmpdir:
        result, pred_dir = _run_boltz_predict(
            input_yaml, "test_mps.yaml", tmpdir,
        )

        assert result.returncode == 0, (
            f"boltz predict --accelerator mps failed:\n"
            f"stdout: {result.stdout}\nstderr: {result.stderr}"
        )

        cif_files = _find_files(pred_dir, ".cif")
        assert len(cif_files) > 0, "No .cif output files found"
        assert os.path.getsize(cif_files[0]) > 0, "CIF file is empty"

        json_files = _find_files(pred_dir, ".json")
        assert len(json_files) > 0, "No confidence JSON files found"


def test_predict_affinity_mps():
    """Run boltz predict on MPS with protein+ligand affinity."""
    input_yaml = """\
version: 1
sequences:
  - protein:
      id: A
      sequence: ACDEFGHIKL
      msa: empty
  - ligand:
      id: L1
      smiles: 'N[C@@H](Cc1ccc(O)cc1)C(=O)O'
properties:
  - affinity:
      binder: L1
"""
    with tempfile.TemporaryDirectory() as tmpdir:
        result, pred_dir = _run_boltz_predict(
            input_yaml, "test_mps_affinity.yaml", tmpdir,
        )

        assert result.returncode == 0, (
            f"boltz predict --accelerator mps failed:\n"
            f"stdout: {result.stdout}\nstderr: {result.stderr}"
        )

        affinity_files = _find_files(pred_dir, ".json")
        affinity_file = None
        for f in affinity_files:
            if os.path.basename(f).startswith("affinity_"):
                affinity_file = f
                break

        assert affinity_file is not None, (
            f"No affinity_*.json found. Files: "
            f"{[os.path.basename(f) for f in affinity_files]}"
        )

        with open(affinity_file) as f:
            affinity_data = json.load(f)

        prob = affinity_data.get("affinity_probability_binary")
        assert isinstance(prob, (int, float)) and 0.0 <= prob <= 1.0, (
            f"affinity_probability_binary invalid: {prob!r}"
        )

        pred_val = affinity_data.get("affinity_pred_value")
        assert isinstance(pred_val, (int, float)) and -10.0 <= pred_val <= 10.0, (
            f"affinity_pred_value invalid: {pred_val!r}"
        )
