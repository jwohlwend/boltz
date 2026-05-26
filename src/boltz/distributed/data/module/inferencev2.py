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

# TODO: v2 features not yet supported for distributed sharding:
# - Constraint features (requires compute_constraint_features=True)
# - Template features
# - Affinity module features
#
# NOTE: The following features are produced by the featurizer but not consumed by the model.
# They are silently dropped during distribution and this is intentional:
# - token_to_center_atom: produced but unused by model
# - ensemble_ref_idxs: consumed during featurization only, unused by model

import math
import warnings
from pathlib import Path
from typing import Dict, Optional

import numpy as np
import pytorch_lightning as pl
import torch
from torch.distributed.device_mesh import DeviceMesh
from torch.distributed.tensor import DTensor
from torch.utils.data import DataLoader, DistributedSampler

from boltz.data import const
from boltz.data.feature.featurizerv2 import Boltz2Featurizer
from boltz.data.module.inferencev2 import load_input
from boltz.data.mol import load_canonicals, load_molecules
from boltz.data.pad import pad_dim
from boltz.data.tokenize.boltz2 import Boltz2Tokenizer
from boltz.data.types import Manifest
from boltz.distributed.data.feature.featurizer import pad_and_scatter_atom_features_dtensor
from boltz.distributed.data.feature.featurizer_utils import get_pair_mask
from boltz.distributed.data.module.placements import INFERENCE_FEATURE_PLACEMENTS_V2
from boltz.distributed.data.types import PairMaskMode
from boltz.distributed.data.utils import (
    ATOM_FEATURES_V2,
    CollateDTensor,
    distribute_features,
    get_flattened_group,
)


class PredictionDatasetCPWithDTensorV2(torch.utils.data.Dataset):
    """Prediction dataset with DTensor context parallelism for Boltz2."""

    def __init__(
        self,
        manifest: Manifest,
        target_dir: Path,
        msa_dir: Path,
        mol_dir: Path,
        device_mesh: DeviceMesh,
        device_mesh_cpu: DeviceMesh,
        constraints_dir: Optional[Path] = None,
        template_dir: Optional[Path] = None,
        extra_mols_dir: Optional[Path] = None,
        max_msa_seqs: int = const.max_msa_seqs,
        msa_pad_to_max_seqs: bool = False,
        max_data_retries: int = 5,
        pair_mask_mode: PairMaskMode = PairMaskMode.NONE,
        atoms_per_window_queries: Optional[int] = 32,
        atoms_per_window_keys: Optional[int] = 128,
        num_ensembles: int = 1,
        per_shard_token_multiple: int = 1,
    ) -> None:
        super().__init__()

        self.manifest = manifest
        self.target_dir = target_dir
        self.msa_dir = msa_dir
        self.mol_dir = mol_dir

        if constraints_dir is not None:
            raise NotImplementedError("Constraints are not supported for CP")

        self.constraints_dir = constraints_dir
        self.template_dir = template_dir
        self.extra_mols_dir = extra_mols_dir
        self.max_msa_seqs = max_msa_seqs
        self.msa_pad_to_max_seqs = msa_pad_to_max_seqs
        self.device_mesh = device_mesh
        self.device_mesh_cpu = device_mesh_cpu
        self.tokenizer = Boltz2Tokenizer()
        self.featurizer = Boltz2Featurizer()
        self.canonicals = load_canonicals(self.mol_dir)
        self.max_data_retries = max_data_retries
        self.num_ensembles = num_ensembles
        self.per_shard_token_multiple = per_shard_token_multiple

        if (atoms_per_window_queries is None) != (atoms_per_window_keys is None):
            raise ValueError(
                "atoms_per_window_queries and atoms_per_window_keys must be either both None or both not None"
            )
        if pair_mask_mode == PairMaskMode.SEQUENCE_LOCAL_ATTENTION and atoms_per_window_queries is None:
            raise ValueError("atoms_per_window_queries must not be None if pair_mask_mode is SequenceLocalAttention")
        self.pair_mask_mode = pair_mask_mode
        self.atoms_per_window_queries = atoms_per_window_queries
        self.atoms_per_window_keys = atoms_per_window_keys

        self._cp_submesh = device_mesh_cpu[("cp_axis_0_cpu", "cp_axis_1_cpu")]
        self._cp_submesh_group = get_flattened_group(self._cp_submesh, backend="gloo")
        self.is_cp_rank_zero = tuple(device_mesh.get_coordinate()[1:]) == (0, 0)

        n_shards_axis_0 = self.device_mesh.shape[1]
        if self.max_msa_seqs % n_shards_axis_0 != 0:
            if self.msa_pad_to_max_seqs is False:
                warnings.warn(
                    f"Number CP ranks along process group grid axis 0 {n_shards_axis_0} is not "
                    f"a integer divisor of max_msa_seqs {self.max_msa_seqs}. Will modify max_msa_seqs "
                    f"to a multiple of {n_shards_axis_0} and pad the MSA number of sequences to it"
                )
                self.msa_pad_to_max_seqs = True

        self.FEATURE_TO_DTENSOR_PLACEMENT = INFERENCE_FEATURE_PLACEMENTS_V2
        self._fallback_depth = 0

    @property
    def use_window_batching(self) -> bool:
        return self.pair_mask_mode == PairMaskMode.NONE

    def _raise_or_return_item_0(self, e: Exception) -> None:
        if self.max_data_retries <= 0:
            raise e
        if self._fallback_depth >= self.max_data_retries:
            raise RuntimeError(
                f"Data loading failed {self.max_data_retries} consecutive times. " f"Last error: {e}"
            ) from e
        self._fallback_depth += 1
        try:
            fallback_idx = np.random.randint(0, len(self))
            return self.__getitem__(fallback_idx)
        finally:
            self._fallback_depth -= 1

    def __getitem__(self, idx: int) -> Dict[str, DTensor]:
        if self.is_cp_rank_zero:
            record = self.manifest.records[idx]

            try:
                input_data = load_input(
                    record,
                    self.target_dir,
                    self.msa_dir,
                    constraints_dir=self.constraints_dir,
                    template_dir=self.template_dir,
                    extra_mols_dir=self.extra_mols_dir,
                )
            except Exception as e:  # noqa: BLE001
                print(f"Data loading failed on {record.id} with error {e}. Skipping.")  # noqa: T201
                return self._raise_or_return_item_0(e)

            try:
                tokenized = self.tokenizer.tokenize(input_data)
            except Exception as e:  # noqa: BLE001
                print(f"Tokenizer failed on {record.id} with error {e}. Skipping.")  # noqa: T201
                return self._raise_or_return_item_0(e)

            try:
                molecules = {}
                molecules.update(self.canonicals)
                if input_data.extra_mols:
                    molecules.update(input_data.extra_mols)
                mol_names = set(tokenized.tokens["res_name"].tolist())
                mol_names = mol_names - set(molecules.keys())
                molecules.update(load_molecules(self.mol_dir, mol_names))
            except Exception as e:  # noqa: BLE001
                print(f"Molecule loading failed for {record.id} with error {e}. Skipping.")  # noqa: T201
                return self._raise_or_return_item_0(e)

            seed = 42
            random = np.random.default_rng(seed)

            # Pad dimensions to be divisible by the CP shard dimension (and atoms_per_window_queries for atoms)
            n_shards_axis_0 = self.device_mesh.shape[1]
            W = self.atoms_per_window_queries or 32

            max_tokens = tokenized.tokens.shape[0]
            max_seqs = self.max_msa_seqs
            pad_to_max_seqs = self.msa_pad_to_max_seqs

            token_align = n_shards_axis_0 * self.per_shard_token_multiple
            if max_tokens % token_align != 0:
                max_tokens = ((max_tokens + token_align - 1) // token_align) * token_align

            if self.use_window_batching:
                max_atoms = None
            else:
                max_atoms = int(np.sum(tokenized.tokens["atom_num"])) if len(tokenized.tokens) > 0 else 0
                # Must be divisible by both atoms_per_window_queries and n_shards_axis_0
                atom_align = math.lcm(W, n_shards_axis_0)
                if max_atoms % atom_align != 0:
                    max_atoms = ((max_atoms + atom_align - 1) // atom_align) * atom_align

            if max_seqs % n_shards_axis_0 != 0:
                max_seqs = max_seqs + n_shards_axis_0 - max_seqs % n_shards_axis_0

            try:
                features = self.featurizer.process(
                    tokenized,
                    molecules=molecules,
                    random=random,
                    training=False,
                    max_atoms=max_atoms,
                    max_tokens=max_tokens,
                    max_seqs=max_seqs,
                    pad_to_max_seqs=pad_to_max_seqs,
                    atoms_per_window_queries=self.atoms_per_window_queries,
                    num_ensembles=self.num_ensembles,
                    fix_single_ensemble=self.num_ensembles == 1,
                    compute_frames=True,
                    compute_constraint_features=False,
                )
                # Distributed-specific features not produced by the base featurizer
                mask = features["token_pad_mask"]
                features["token_pair_pad_mask"] = mask[:, None] * mask[None, :]
                if self.pair_mask_mode == PairMaskMode.GLOBAL_ATOM_ATTENTION:
                    N_atoms = len(features["ref_pos"])
                    features["pair_mask"] = torch.ones(N_atoms, N_atoms, dtype=torch.float)
                elif self.pair_mask_mode == PairMaskMode.SEQUENCE_LOCAL_ATTENTION:
                    features["pair_mask"] = get_pair_mask(
                        N_atoms=len(features["ref_pos"]), W=self.atoms_per_window_queries
                    )
            except Exception as e:  # noqa: BLE001
                print(f"Featurizer failed on {record.id} with error {e}. Skipping.")  # noqa: T201
                return self._raise_or_return_item_0(e)

            if not pad_to_max_seqs:
                num_seqs_actual = features["msa_mask"].shape[0]
                target_seqs = max(1, num_seqs_actual)
                if target_seqs % n_shards_axis_0 != 0:
                    target_seqs = target_seqs + n_shards_axis_0 - target_seqs % n_shards_axis_0
                if num_seqs_actual < target_seqs:
                    pad_len = target_seqs - num_seqs_actual
                    msa_feature_keys = ("msa", "msa_paired", "deletion_value", "has_deletion", "msa_mask")
                    for key in msa_feature_keys:
                        if key in features:
                            features[key] = (
                                pad_dim(features[key], 0, pad_len, const.token_ids["-"])
                                if key == "msa"
                                else pad_dim(features[key], 0, pad_len)
                            )

            record_list = [record]
            atom_placements = {
                key: self.FEATURE_TO_DTENSOR_PLACEMENT[key]
                for key in ATOM_FEATURES_V2
                if key in self.FEATURE_TO_DTENSOR_PLACEMENT
            }
            atom_features = {key: features[key] for key in features if key in atom_placements}
            token_and_msa_features = {
                key: features[key]
                for key in features
                if key not in ATOM_FEATURES_V2 and key in self.FEATURE_TO_DTENSOR_PLACEMENT
            }
            token_and_msa_placements = {
                key: self.FEATURE_TO_DTENSOR_PLACEMENT[key]
                for key in self.FEATURE_TO_DTENSOR_PLACEMENT
                if key not in ATOM_FEATURES_V2
            }

        else:
            features = None
            record_list = [None]
            atom_features = None
            atom_placements = {
                key: self.FEATURE_TO_DTENSOR_PLACEMENT[key]
                for key in ATOM_FEATURES_V2
                if key in self.FEATURE_TO_DTENSOR_PLACEMENT
            }
            token_and_msa_features = None
            token_and_msa_placements = {
                key: self.FEATURE_TO_DTENSOR_PLACEMENT[key]
                for key in self.FEATURE_TO_DTENSOR_PLACEMENT
                if key not in ATOM_FEATURES_V2
            }

        if self.use_window_batching:
            atom_placements.pop("pair_mask")

        cp_submesh = self._cp_submesh
        cp_submesh_group = self._cp_submesh_group
        cp_group_src_rank_global = min(torch.distributed.get_process_group_ranks(cp_submesh_group))

        atom_features_dtensor = pad_and_scatter_atom_features_dtensor(
            features=atom_features,
            placements=atom_placements,
            group=cp_submesh_group,
            src_rank_global=cp_group_src_rank_global,
            device_mesh=cp_submesh,
        )

        token_and_msa_features_dtensor = distribute_features(
            features=token_and_msa_features,
            placements=token_and_msa_placements,
            group=cp_submesh_group,
            src_rank_global=cp_group_src_rank_global,
            device_mesh=cp_submesh,
        )

        features_dtensor = {**token_and_msa_features_dtensor, **atom_features_dtensor}

        torch.distributed.broadcast_object_list(record_list, src=cp_group_src_rank_global, group=cp_submesh_group)
        features_dtensor["record"] = record_list[0]

        return features_dtensor

    def __len__(self) -> int:
        return len(self.manifest.records)


class Boltz2InferenceDataModuleDTensor(pl.LightningDataModule):
    """DataModule for Boltz2 distributed inference with DTensor CP."""

    def __init__(
        self,
        manifest: Manifest,
        target_dir: Path,
        msa_dir: Path,
        mol_dir: Path,
        num_workers: int,
        device_mesh: DeviceMesh,
        device_mesh_cpu: DeviceMesh,
        constraints_dir: Optional[Path] = None,
        template_dir: Optional[Path] = None,
        extra_mols_dir: Optional[Path] = None,
        max_msa_seqs: int = const.max_msa_seqs,
        msa_pad_to_max_seqs: bool = False,
        max_data_retries: int = 5,
        pair_mask_mode: PairMaskMode = PairMaskMode.NONE,
        atoms_per_window_queries: int = 32,
        atoms_per_window_keys: int = 128,
        local_batch_size: int = 1,
        num_ensembles: int = 1,
        per_shard_token_multiple: int = 1,
    ) -> None:
        super().__init__()
        if num_workers != 0:
            raise NotImplementedError("num_workers != 0 is not supported for CP")
        self.num_workers = num_workers
        self.manifest = manifest
        self.target_dir = target_dir
        self.msa_dir = msa_dir
        self.mol_dir = mol_dir
        if constraints_dir is not None:
            raise NotImplementedError("Constraints are not supported for CP")
        self.constraints_dir = constraints_dir
        self.template_dir = template_dir
        self.extra_mols_dir = extra_mols_dir
        self.max_msa_seqs = max_msa_seqs
        self.msa_pad_to_max_seqs = msa_pad_to_max_seqs
        self.device_mesh = device_mesh
        self.device_mesh_cpu = device_mesh_cpu
        self.max_data_retries = max_data_retries
        self.pair_mask_mode = pair_mask_mode
        self.atoms_per_window_queries = atoms_per_window_queries
        self.atoms_per_window_keys = atoms_per_window_keys
        self.dataset: Optional[PredictionDatasetCPWithDTensorV2] = None
        self.local_batch_size = local_batch_size
        self.num_ensembles = num_ensembles
        self.per_shard_token_multiple = per_shard_token_multiple

    def setup(self, stage: Optional[str] = None) -> None:
        if stage != "predict":
            raise ValueError(f"Only predict stage is supported for inference but got {stage}")

        self.dataset = PredictionDatasetCPWithDTensorV2(
            manifest=self.manifest,
            target_dir=self.target_dir,
            msa_dir=self.msa_dir,
            mol_dir=self.mol_dir,
            device_mesh=self.device_mesh,
            device_mesh_cpu=self.device_mesh_cpu,
            constraints_dir=self.constraints_dir,
            template_dir=self.template_dir,
            extra_mols_dir=self.extra_mols_dir,
            max_msa_seqs=self.max_msa_seqs,
            msa_pad_to_max_seqs=self.msa_pad_to_max_seqs,
            max_data_retries=self.max_data_retries,
            pair_mask_mode=self.pair_mask_mode,
            num_ensembles=self.num_ensembles,
            per_shard_token_multiple=self.per_shard_token_multiple,
        )

    def predict_dataloader(self) -> DataLoader:
        sampler = DistributedSampler(
            self.dataset,
            num_replicas=self.device_mesh_cpu.shape[0],
            rank=self.device_mesh_cpu.get_local_rank(0),
            shuffle=False,
            drop_last=False,
        )
        custom_collate = CollateDTensor(self.device_mesh_cpu)

        return DataLoader(
            self.dataset,
            batch_size=self.local_batch_size,
            num_workers=self.num_workers,
            pin_memory=False,
            shuffle=False,
            collate_fn=custom_collate,
            sampler=sampler,
        )

    def transfer_batch_to_device(
        self,
        batch: dict,
        device: torch.device,
        dataloader_idx: int,  # noqa: ARG002
    ) -> dict:
        for key in batch:
            if key not in {"record"}:
                batch_local = batch[key].to_local().to(device)
                batch[key] = DTensor.from_local(
                    batch_local,
                    device_mesh=self.device_mesh,
                    placements=batch[key].placements,
                    shape=batch[key].shape,
                    stride=batch[key].stride(),
                )

        return batch
