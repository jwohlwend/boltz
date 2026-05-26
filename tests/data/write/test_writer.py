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

"""Tests for BoltzWriter multi-diffusion-sample correctness.

Verifies that the current BoltzWriter produces identical CIF output to the
dev-v2 reference version.  The dev-v2 writer crashes with ``KeyError`` when
``diffusion_samples > 1`` and no confidence scores (``idx_to_rank`` is sized
by ``len(records)`` instead of the number of diffusion samples).  We work
around this by calling the dev-v2 reference one sample at a time, then verify
the current writer (called once with all samples) produces identical output.
"""

from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest
import torch
from torch import Tensor

from boltz.data.types import (
    AtomV2,
    BondV2,
    Chain,
    ChainInfo,
    Coords,
    Ensemble,
    Interface,
    Record,
    Residue,
    StructureInfo,
    StructureV2,
)
from boltz.data.write.mmcif import to_mmcif
from boltz.data.write.writer import BoltzWriter

# Atom names used to populate multi-atom residues without duplicates.
_ATOM_NAMES = ["N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ", "OG"]


# ---------------------------------------------------------------------------
# Reference: verbatim dev-v2 write_on_batch_end logic
# ---------------------------------------------------------------------------
def _dev_v2_write_on_batch_end(
    data_dir: Path,
    output_dir: Path,
    prediction: dict[str, Tensor],
    batch: dict[str, Tensor],
) -> None:
    """Faithful reproduction of the dev-v2 BoltzWriter.write_on_batch_end.

    Extracted as a plain function so we can call it with single-sample slices
    to avoid the ``idx_to_rank`` KeyError.
    """
    if prediction["exception"]:
        return

    records: list[Record] = batch["record"]

    coords = prediction["coords"]
    coords = coords.unsqueeze(0)

    pad_masks = prediction["masks"]

    if "confidence_score" in prediction:
        argsort = torch.argsort(prediction["confidence_score"], descending=True)
        idx_to_rank = {idx.item(): rank for rank, idx in enumerate(argsort)}
    else:
        idx_to_rank = {i: i for i in range(len(records))}

    for record, coord, pad_mask in zip(records, coords, pad_masks):
        path = data_dir / f"{record.id}.npz"
        structure: StructureV2 = StructureV2.load(path)

        chain_map = {}
        for i, mask in enumerate(structure.mask):
            if mask:
                chain_map[len(chain_map)] = i

        structure = structure.remove_invalid_chains()

        for model_idx in range(coord.shape[0]):
            model_coord = coord[model_idx]
            coord_unpad = model_coord[pad_mask.bool()]
            coord_unpad = coord_unpad.cpu().numpy()

            atoms = structure.atoms
            atoms["coords"] = coord_unpad
            atoms["is_present"] = True
            structure: StructureV2
            coord_unpad = [(x,) for x in coord_unpad]
            coord_unpad = np.array(coord_unpad, dtype=Coords)

            residues = structure.residues
            residues["is_present"] = True

            interfaces = np.array([], dtype=Interface)
            new_structure: StructureV2 = replace(
                structure,
                atoms=atoms,
                residues=residues,
                interfaces=interfaces,
                coords=coord_unpad,
            )

            chain_info = []
            for chain in new_structure.chains:
                old_chain_idx = chain_map[chain["asym_id"]]
                old_chain_info = record.chains[old_chain_idx]
                new_chain_info = replace(
                    old_chain_info,
                    chain_id=int(chain["asym_id"]),
                    valid=True,
                )
                chain_info.append(new_chain_info)

            struct_dir = output_dir / record.id
            struct_dir.mkdir(parents=True, exist_ok=True)

            plddts = None
            if "plddt" in prediction:
                plddts = prediction["plddt"][model_idx]

            outname = f"{record.id}_model_{idx_to_rank[model_idx]}"
            cif_path = struct_dir / f"{outname}.cif"
            with cif_path.open("w") as f:
                f.write(to_mmcif(new_structure, plddts=plddts, boltz2=True))


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_structure_v2(n_atoms: int, n_residues: int, rng: np.random.Generator) -> StructureV2:
    """Create a minimal valid StructureV2 with one chain.

    Each residue gets ``n_atoms // n_residues`` atoms with distinct names
    drawn from ``_ATOM_NAMES`` so that ``to_mmcif`` does not reject them.
    """
    atoms_per_res = n_atoms // n_residues
    assert atoms_per_res * n_residues == n_atoms
    assert atoms_per_res <= len(_ATOM_NAMES)

    atoms = np.zeros(n_atoms, dtype=AtomV2)
    atoms["coords"] = rng.standard_normal((n_atoms, 3)).astype(np.float32)
    atoms["is_present"] = True
    for i in range(n_atoms):
        atoms[i]["name"] = _ATOM_NAMES[i % atoms_per_res]

    residues = np.zeros(n_residues, dtype=Residue)
    for i in range(n_residues):
        residues[i]["name"] = "ALA"
        residues[i]["res_type"] = 0
        residues[i]["res_idx"] = i
        residues[i]["atom_idx"] = i * atoms_per_res
        residues[i]["atom_num"] = atoms_per_res
        residues[i]["atom_center"] = i * atoms_per_res
        residues[i]["atom_disto"] = i * atoms_per_res
        residues[i]["is_standard"] = True
        residues[i]["is_present"] = True

    chains = np.zeros(1, dtype=Chain)
    chains[0]["name"] = "A"
    chains[0]["mol_type"] = 0
    chains[0]["entity_id"] = 0
    chains[0]["sym_id"] = 0
    chains[0]["asym_id"] = 0
    chains[0]["atom_idx"] = 0
    chains[0]["atom_num"] = n_atoms
    chains[0]["res_idx"] = 0
    chains[0]["res_num"] = n_residues
    chains[0]["cyclic_period"] = 0

    bonds = np.array([], dtype=BondV2)
    interfaces = np.array([], dtype=Interface)
    mask = np.array([True])
    coords = np.array([(c,) for c in atoms["coords"]], dtype=Coords)
    ensemble = np.zeros(1, dtype=Ensemble)
    ensemble[0]["atom_coord_idx"] = 0
    ensemble[0]["atom_num"] = n_atoms

    return StructureV2(
        atoms=atoms,
        bonds=bonds,
        residues=residues,
        chains=chains,
        interfaces=interfaces,
        mask=mask,
        coords=coords,
        ensemble=ensemble,
    )


def _make_record(record_id: str, n_residues: int) -> Record:
    chain_info = ChainInfo(
        chain_id=0,
        chain_name="A",
        mol_type=0,
        cluster_id=0,
        msa_id=0,
        num_residues=n_residues,
        valid=True,
    )
    return Record(
        id=record_id,
        structure=StructureInfo(),
        chains=[chain_info],
        interfaces=[],
    )


def _save_structure_npz(structure: StructureV2, path: Path) -> None:
    """Save StructureV2 to npz, omitting None-valued optional fields."""
    save_dict = {k: v for k, v in vars(structure).items() if v is not None}
    np.savez(str(path), **save_dict)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_diffusion_samples", [1, 5])
def test_boltz_writer_multi_sample_parity(tmp_path, n_diffusion_samples):
    """Current writer (all samples at once) == dev-v2 reference (per sample).

    For each diffusion sample we call the dev-v2 reference with a single-sample
    ``coords`` tensor (shape ``(1, N, 3)``) so ``unsqueeze(0)`` yields
    ``(1, 1, N, 3)`` and ``idx_to_rank = {0: 0}`` — no crash.  The current
    writer receives all samples at once and writes ``model_0 … model_{k-1}``.
    """
    rng = np.random.default_rng(42)
    n_atoms = 30
    n_residues = 10
    n_pad = 2
    n_atoms_padded = n_atoms + n_pad
    record_id = "test_struct"

    # Persist structure to disk
    data_dir = tmp_path / "structures"
    data_dir.mkdir()
    structure = _make_structure_v2(n_atoms, n_residues, rng)
    _save_structure_npz(structure, data_dir / record_id)

    record = _make_record(record_id, n_residues)

    # Deterministic random coords: (n_diffusion_samples, n_atoms_padded, 3)
    torch.manual_seed(123)
    all_coords = torch.randn(n_diffusion_samples, n_atoms_padded, 3)
    pad_mask = torch.zeros(n_atoms_padded, dtype=torch.bool)
    pad_mask[:n_atoms] = True

    # ---- Reference: call dev-v2 logic once per diffusion sample ----
    ref_dir = tmp_path / "ref_output"
    ref_cifs: dict[int, str] = {}
    for sample_idx in range(n_diffusion_samples):
        sample_out = ref_dir / f"_sample_{sample_idx}"
        _dev_v2_write_on_batch_end(
            data_dir=data_dir,
            output_dir=sample_out,
            prediction={
                "exception": False,
                "coords": all_coords[sample_idx : sample_idx + 1],
                "masks": pad_mask.unsqueeze(0),
            },
            batch={"record": [record]},
        )
        cif_path = sample_out / record_id / f"{record_id}_model_0.cif"
        assert cif_path.exists(), f"Reference CIF missing: {cif_path}"
        ref_cifs[sample_idx] = cif_path.read_text()

    # ---- Current writer: single call with all samples ----
    cur_dir = tmp_path / "cur_output"
    cur_writer = BoltzWriter(
        data_dir=str(data_dir),
        output_dir=str(cur_dir),
        output_format="mmcif",
        boltz2=True,
    )
    cur_writer.write_on_batch_end(
        trainer=None,
        pl_module=None,
        prediction={
            "exception": False,
            "coords": all_coords,
            "masks": pad_mask.unsqueeze(0),
        },
        batch_indices=None,
        batch={"record": [record]},
        batch_idx=0,
        dataloader_idx=0,
    )

    # ---- Compare ----
    cur_struct_dir = cur_dir / record_id
    assert cur_struct_dir.exists(), f"Current output dir missing: {cur_struct_dir}"
    cur_cif_files = sorted(cur_struct_dir.glob("*.cif"))
    assert (
        len(cur_cif_files) == n_diffusion_samples
    ), f"Expected {n_diffusion_samples} CIF files, found {len(cur_cif_files)}"

    for sample_idx in range(n_diffusion_samples):
        cur_cif_path = cur_struct_dir / f"{record_id}_model_{sample_idx}.cif"
        assert cur_cif_path.exists(), f"Current CIF missing: {cur_cif_path}"
        cur_text = cur_cif_path.read_text()
        assert cur_text == ref_cifs[sample_idx], (
            f"CIF mismatch for diffusion sample {sample_idx}:\n"
            f"  current file : {cur_cif_path}\n"
            f"  reference    : dev-v2 per-sample call"
        )
