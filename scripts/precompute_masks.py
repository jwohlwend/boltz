from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


def compute_masks(data: np.lib.npyio.NpzFile) -> tuple[np.ndarray, np.ndarray]:
    """Return (alignment_mask, rmsd_mask) inferred from raw NPZ contents."""
    chains = data["chains"]
    chain_mask = data["mask"].astype(bool)

    alignment_segments: list[np.ndarray] = []
    rmsd_segments: list[np.ndarray] = []

    entity_counts: dict[int, int] = {}
    for idx, chain in enumerate(chains):
        if not chain_mask[idx]:
            continue
        entity_id = int(chain["entity_id"])
        entity_counts[entity_id] = entity_counts.get(entity_id, 0) + int(chain["res_num"])

    if not entity_counts:
        return np.zeros(0, dtype=bool), np.zeros(0, dtype=bool)

    # For now, assume the peptide entity is the smallest one by residue count, and MHC/heavy is the largest.
    # Will need to revist.
    peptide_entity = min(entity_counts, key=entity_counts.get)
    heavy_candidates = {k: v for k, v in entity_counts.items() if k != peptide_entity}
    heavy_entity = (
        max(heavy_candidates, key=heavy_candidates.get)
        if heavy_candidates
        else peptide_entity
    )

    for idx, chain in enumerate(chains):
        if not chain_mask[idx]:
            continue
        entity_id = int(chain["entity_id"])
        atom_num = int(chain["atom_num"])

        align_segment = np.zeros(atom_num, dtype=bool)
        rmsd_segment = np.zeros(atom_num, dtype=bool)

        if entity_id == heavy_entity:
            align_segment[:] = True
        if entity_id == peptide_entity:
            rmsd_segment[:] = True

        alignment_segments.append(align_segment)
        rmsd_segments.append(rmsd_segment)

    alignment_mask = (
        np.concatenate(alignment_segments) if alignment_segments else np.zeros(0, dtype=bool)
    )
    rmsd_mask = (
        np.concatenate(rmsd_segments) if rmsd_segments else np.zeros(0, dtype=bool)
    )

    return alignment_mask, rmsd_mask


def process_structure(npz_path: Path, dry_run: bool = False) -> None:
    with np.load(npz_path) as data:
        alignment_mask, rmsd_mask = compute_masks(data)

        if dry_run:
            print(
                f"{npz_path.name}: align atoms={alignment_mask.sum()} "
                f"(len={alignment_mask.shape[0]}), "
                f"rmsd atoms={rmsd_mask.sum()} (len={rmsd_mask.shape[0]})"
            )
            return

        updated = {key: data[key] for key in data.files}
        updated["alignment_mask"] = alignment_mask
        updated["rmsd_mask"] = rmsd_mask

    tmp_path = npz_path.with_suffix(".tmp.npz")
    np.savez_compressed(tmp_path, **updated)
    tmp_path.replace(npz_path)


def main() -> None:
    parser = argparse.ArgumentParser(description="Inject alignment/rmsd masks into structure NPZ files.")
    parser.add_argument(
        "structures_dir",
        type=Path,
        help="Directory containing processed structure NPZ files.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Report mask statistics without writing files.",
    )
    args = parser.parse_args()

    npz_paths = sorted(args.structures_dir.glob("*.npz"))
    if not npz_paths:
        raise SystemExit(f"No .npz files found under {args.structures_dir}")

    for npz_path in npz_paths:
        process_structure(npz_path, dry_run=args.dry_run)

    if args.dry_run:
        print("Dry run complete; no files modified.")
    else:
        print(f"Processed {len(npz_paths)} structures.")


if __name__ == "__main__":
    main()
