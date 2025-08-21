from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List, Dict, Any

import csv as pycsv
import numpy as np
import pytorch_lightning as pl
import torch
from torch.utils.data import DataLoader

from boltz.data.pad import pad_to_max
from boltz.data.tokenize.boltz2 import Boltz2Tokenizer
from boltz.data.feature.featurizerv2 import Boltz2Featurizer
from boltz.data.parse.schema import parse_boltz_schema
from boltz.data.types import Input

@dataclass
class DataConfig:
    csv_path: str
    tokenizer: Boltz2Tokenizer
    featurizer: Boltz2Featurizer
    batch_size: int
    num_workers: int
    pin_memory: bool
    random_seed: int
    max_atoms: int = 0
    max_tokens: int = 0
    max_seqs: int = 0
    pad_to_max_tokens: bool = False
    pad_to_max_atoms: bool = False
    pad_to_max_seqs: bool = False
    atoms_per_window_queries: int = 32
    min_dist: float = 2.0
    max_dist: float = 22.0
    num_bins: int = 64
    val_batch_size: int = 1
    polymer_property_keys: Optional[List[str]] = None


def collate(data: List[Dict[str, torch.Tensor]]) -> Dict[str, torch.Tensor]:
    keys = data[0].keys()
    collated: Dict[str, torch.Tensor] = {}
    for key in keys:
        values = [d[key] for d in data]
        if key not in [
            "all_coords",
            "all_resolved_mask",
            "crop_to_all_atom_map",
            "chain_symmetries",
            "amino_acids_symmetries",
            "ligand_symmetries",
        ]:
            shape = values[0].shape
            if not all(v.shape == shape for v in values):
                values, _ = pad_to_max(values, 0)
            else:
                values = torch.stack(values, dim=0)
        collated[key] = values
    return collated


class PolymerCSVDataset(torch.utils.data.Dataset):
    def __init__(
        self,
        rows: List[Dict[str, Any]],
        tokenizer: Boltz2Tokenizer,
        featurizer: Boltz2Featurizer,
        max_atoms: int,
        max_tokens: int,
        max_seqs: int,
        pad_to_max_atoms: bool,
        pad_to_max_tokens: bool,
        pad_to_max_seqs: bool,
        atoms_per_window_queries: int,
        min_dist: float,
        max_dist: float,
        num_bins: int,
        random_seed: int,
        polymer_property_keys: Optional[List[str]] = ["Tg", "FFV", "Tc", "Density", "Rg"],
    ) -> None:
        super().__init__()
        self.rows = rows
        self.tokenizer = tokenizer
        self.featurizer = featurizer
        self.max_atoms = max_atoms
        self.max_tokens = max_tokens
        self.max_seqs = max_seqs
        self.pad_to_max_atoms = pad_to_max_atoms
        self.pad_to_max_tokens = pad_to_max_tokens
        self.pad_to_max_seqs = pad_to_max_seqs
        self.atoms_per_window_queries = atoms_per_window_queries
        self.min_dist = min_dist
        self.max_dist = max_dist
        self.num_bins = num_bins
        self.random = np.random.default_rng(random_seed)
        self.property_keys = polymer_property_keys

    def __len__(self) -> int:
        return len(self.rows)

    def __getitem__(self, idx: int) -> Dict[str, torch.Tensor]:
        row = self.rows[idx]
        smiles = row["smiles"]

        # Build minimal schema and parse to Target (constructs structure from SMILES via RDKit)
        schema = {"version": 1, "sequences": [{"ligand": {"id": "L", "smiles": smiles}}]}
        target = parse_boltz_schema(name=str(idx), schema=schema, ccd={}, mol_dir=None, boltz_2=True)

        # Create Input (no MSA); include a minimal Record so tokenizer sees .record.affinity
        from boltz.data.types import Record, StructureInfo
        record = Record(
            id=str(idx),
            structure=StructureInfo(),
            chains=[],
            interfaces=[],
            affinity=None,
            polymer_properties=None,
        )
        input_data = Input(target.structure, {}, record=record, extra_mols=target.extra_mols)

        # Tokenize
        tokenized = self.tokenizer.tokenize(input_data)

        # Featurize
        # Minimal molecules dict (empty) and RNG for v2 featurizer
        rng = np.random.default_rng(self.random.integers(0, 1 << 31))
        features = self.featurizer.process(
            tokenized,
            random=rng,
            molecules=input_data.extra_mols or {},
            training=True,
            max_seqs=self.max_seqs,
            atoms_per_window_queries=self.atoms_per_window_queries,
            min_dist=self.min_dist,
            max_dist=self.max_dist,
            num_bins=self.num_bins,
            max_tokens=self.max_tokens if self.pad_to_max_tokens else None,
            max_atoms=self.max_atoms if self.pad_to_max_atoms else None,
            pad_to_max_seqs=self.pad_to_max_seqs,
            compute_symmetries=False,
        )

        # Attach properties in configured order; allow case-insensitive CSV columns
        values: List[Optional[float]] = []
        row_lower = {k.lower(): v for k, v in row.items()}
        for key in self.property_keys:
            v = row.get(key, None)
            if v is None:
                v = row_lower.get(key.lower(), None)
            try:
                values.append(float(v))
            except (TypeError, ValueError):
                values.append(float("nan"))
        tensor = torch.tensor(values, dtype=torch.float)
        features["polymer_properties"] = tensor

        return features


class PolymerCSVDataModule(pl.LightningDataModule):
    def __init__(self, cfg: DataConfig) -> None:
        super().__init__()
        self.cfg = cfg
        self.rows: List[Dict[str, Any]] = []

    def setup(self, stage: Optional[str] = None) -> None:
        # Load CSV rows
        csv_path = Path(self.cfg.csv_path)
        with csv_path.open("r") as f:
            reader = pycsv.DictReader(f)
            fieldnames = [fn.strip() for fn in (reader.fieldnames or [])]
            lower_to_actual = {fn.lower(): fn for fn in fieldnames}
            for row in reader:
                norm = {k.strip(): v for k, v in row.items()}
                # Unify case for smiles
                if "smiles" not in norm:
                    for key in list(norm.keys()):
                        if key.lower() == "smiles":
                            norm["smiles"] = norm[key]
                            break
                if "smiles" not in norm or not norm["smiles"]:
                    continue
                # Ensure canonical keys for properties exist (copy if only case-mismatch)
                for k in self.cfg.polymer_property_keys:
                    if k not in norm:
                        src = lower_to_actual.get(k.lower())
                        if src and src in norm:
                            norm[k] = norm[src]
                self.rows.append(norm)

        self.train_set = PolymerCSVDataset(
            rows=self.rows,
            tokenizer=self.cfg.tokenizer,
            featurizer=self.cfg.featurizer,
            max_atoms=self.cfg.max_atoms,
            max_tokens=self.cfg.max_tokens,
            max_seqs=self.cfg.max_seqs,
            pad_to_max_atoms=self.cfg.pad_to_max_atoms,
            pad_to_max_tokens=self.cfg.pad_to_max_tokens,
            pad_to_max_seqs=self.cfg.pad_to_max_seqs,
            atoms_per_window_queries=self.cfg.atoms_per_window_queries,
            min_dist=self.cfg.min_dist,
            max_dist=self.cfg.max_dist,
            num_bins=self.cfg.num_bins,
            random_seed=self.cfg.random_seed,
            polymer_property_keys=self.cfg.polymer_property_keys
        )

    def train_dataloader(self) -> DataLoader:
        return DataLoader(
            self.train_set,
            batch_size=self.cfg.batch_size,
            num_workers=self.cfg.num_workers,
            pin_memory=self.cfg.pin_memory,
            shuffle=False,
            collate_fn=collate,
        )

    def val_dataloader(self) -> DataLoader:
        # No validation metric for CSV ingestion by default
        return DataLoader(
            self.train_set,
            batch_size=self.cfg.val_batch_size,
            num_workers=self.cfg.num_workers,
            pin_memory=self.cfg.pin_memory,
            shuffle=False,
            collate_fn=collate,
        )


