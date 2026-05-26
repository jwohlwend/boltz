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

from typing import Any, Optional

import pytorch_lightning as pl
import torch
from torch.distributed.device_mesh import DeviceMesh
from torch.distributed.tensor import DTensor, Shard
from torch.utils.data import DataLoader, DistributedSampler

from boltz.data import const
from boltz.data.module.trainingv2 import (
    Boltz2TrainingDataModule as Boltz2TrainingDataModuleSerial,
)
from boltz.data.module.trainingv2 import DataConfigV2
from boltz.data.pad import pad_dim
from boltz.distributed.data.feature.featurizer import pad_and_scatter_atom_features_dtensor
from boltz.distributed.data.module.placements import TRAINING_FEATURE_PLACEMENTS_V2
from boltz.distributed.data.utils import (
    ATOM_FEATURES_V2,
    NON_SHARDED_FEATURES_V2,
    CollateDTensor,
    distribute_features,
    get_flattened_group,
)


class _BaseDatasetCPWithDTensorV2(torch.utils.data.Dataset):
    """Wrap a serial Boltz2 dataset and distribute features as DTensors."""

    def __init__(
        self,
        serial_dataset: torch.utils.data.Dataset,
        device_mesh: DeviceMesh,
        device_mesh_cpu: DeviceMesh,
    ) -> None:
        """Initialize the distributed dataset wrapper.

        Wraps a serial Boltz2 dataset and sets up the device mesh and placement
        mapping used to distribute features as DTensors across context-parallel
        ranks.

        Parameters
        ----------
        serial_dataset : torch.utils.data.Dataset
            The serial (single-rank) Boltz2 dataset to wrap.
        device_mesh : DeviceMesh
            Device mesh for distributed tensor operations on GPU.
        device_mesh_cpu : DeviceMesh
            Device mesh for distributed tensor operations on CPU, used for
            data-loading collectives that run before GPU transfer.

        """
        super().__init__()
        self.serial_dataset = serial_dataset
        self.device_mesh = device_mesh
        self.device_mesh_cpu = device_mesh_cpu
        self._cp_submesh = device_mesh_cpu[("cp_axis_0_cpu", "cp_axis_1_cpu")]
        self._cp_submesh_group = get_flattened_group(self._cp_submesh, backend="gloo")
        self.is_cp_rank_zero = tuple(device_mesh.get_coordinate()[1:]) == (0, 0)

        self.feature_to_dtensor_placement = TRAINING_FEATURE_PLACEMENTS_V2

    def __len__(self) -> int:
        """Get the length of the dataset.

        Returns
        -------
        int
            The number of samples in the underlying serial dataset.

        """
        return len(self.serial_dataset)

    def _distribute_features(self, features: Optional[dict[str, Any]]) -> dict[str, Any]:
        """Distribute serial features as DTensors across context-parallel ranks.

        CP rank zero pads tensor features to be evenly shardable, then broadcasts
        tensor keys and non-sharded metadata to all CP ranks. Atom features are
        distributed via ``pad_and_scatter_atom_features_dtensor`` and token/MSA
        features via ``distribute_features``.

        Parameters
        ----------
        features : dict[str, Any] or None
            Feature dictionary produced by the serial dataset on CP rank zero.
            Must be ``None`` on non-zero CP ranks.

        Returns
        -------
        dict[str, Any]
            Feature dictionary where tensor values are DTensors distributed
            according to ``self.feature_to_dtensor_placement``, and non-sharded
            values are broadcast copies.

        """
        cp_group_src_rank_global = min(torch.distributed.get_process_group_ranks(self._cp_submesh_group))

        if self.is_cp_rank_zero:
            # Synthesize token_pair_pad_mask if the serial featurizer did not
            # produce it.  The distributed model forward requires this feature
            # (boltz2.py, trunkv2.py) but the serial v2 training featurizer
            # does not generate it — only the inference featurizer does.
            if "token_pair_pad_mask" not in features and "token_pad_mask" in features:
                mask = features["token_pad_mask"]
                features["token_pair_pad_mask"] = mask[:, None] * mask[None, :]

            unknown_tensor_keys = sorted(
                key
                for key, value in features.items()
                if isinstance(value, torch.Tensor)
                and key not in self.feature_to_dtensor_placement
                and key not in NON_SHARDED_FEATURES_V2
            )
            if unknown_tensor_keys:
                raise KeyError(
                    "Found tensor feature keys without DTensor placement mapping. "
                    f"Please add placements for: {unknown_tensor_keys}"
                )

            tensor_features_all = {
                key: value
                for key, value in features.items()
                if isinstance(value, torch.Tensor)
                and key in self.feature_to_dtensor_placement
                and key not in NON_SHARDED_FEATURES_V2
            }

            n_shards_axis_0 = self.device_mesh.shape[1]
            for key, tensor in tensor_features_all.items():
                placements = self.feature_to_dtensor_placement[key]
                padded = tensor
                for placement in placements:
                    if not isinstance(placement, Shard):
                        continue
                    shard_dim = placement.dim
                    if shard_dim >= padded.ndim:
                        continue
                    remainder = padded.shape[shard_dim] % n_shards_axis_0
                    if remainder == 0:
                        continue
                    pad_len = n_shards_axis_0 - remainder
                    pad_value = const.token_ids["-"] if key == "msa" and shard_dim == 0 else 0
                    padded = pad_dim(padded, shard_dim, pad_len, pad_value)
                tensor_features_all[key] = padded
            tensor_feature_keys = sorted(tensor_features_all.keys())
            keys_payload = [tensor_feature_keys]
            torch.distributed.broadcast_object_list(
                keys_payload,
                src=cp_group_src_rank_global,
                group=self._cp_submesh_group,
            )
            tensor_feature_keys_shared = keys_payload[0]

            atom_placements = {
                key: self.feature_to_dtensor_placement[key]
                for key in ATOM_FEATURES_V2
                if key in tensor_feature_keys_shared
            }
            atom_features = {key: tensor_features_all[key] for key in atom_placements}
            token_and_msa_features = {
                key: value for key, value in tensor_features_all.items() if key not in ATOM_FEATURES_V2
            }
            token_and_msa_placements = {key: self.feature_to_dtensor_placement[key] for key in token_and_msa_features}
            non_sharded_features = {key: value for key, value in features.items() if key in NON_SHARDED_FEATURES_V2}
            object_payload = [non_sharded_features]
        else:
            keys_payload = [None]
            torch.distributed.broadcast_object_list(
                keys_payload,
                src=cp_group_src_rank_global,
                group=self._cp_submesh_group,
            )
            tensor_feature_keys_shared = keys_payload[0]
            atom_features = None
            token_and_msa_features = None
            atom_placements = {
                key: self.feature_to_dtensor_placement[key]
                for key in ATOM_FEATURES_V2
                if key in tensor_feature_keys_shared
            }
            token_and_msa_placements = {
                key: placement
                for key, placement in self.feature_to_dtensor_placement.items()
                if key in tensor_feature_keys_shared and key not in ATOM_FEATURES_V2
            }
            object_payload = [None]

        atom_features_dtensor = pad_and_scatter_atom_features_dtensor(
            features=atom_features,
            placements=atom_placements,
            group=self._cp_submesh_group,
            src_rank_global=cp_group_src_rank_global,
            device_mesh=self._cp_submesh,
        )
        token_and_msa_features_dtensor = distribute_features(
            features=token_and_msa_features,
            placements=token_and_msa_placements,
            group=self._cp_submesh_group,
            src_rank_global=cp_group_src_rank_global,
            device_mesh=self._cp_submesh,
        )

        torch.distributed.broadcast_object_list(
            object_payload,
            src=cp_group_src_rank_global,
            group=self._cp_submesh_group,
        )

        features_dtensor = {**token_and_msa_features_dtensor, **atom_features_dtensor}
        features_dtensor.update(object_payload[0] or {})
        return features_dtensor


class TrainingDatasetCPWithDTensorV2(_BaseDatasetCPWithDTensorV2):
    """Training dataset with DTensor context parallelism for Boltz2."""

    def __getitem__(self, idx: int) -> dict[str, Any]:
        """Fetch and distribute a single training sample.

        CP rank zero retrieves the sample from the serial dataset; all other
        CP ranks receive the distributed DTensor features via collectives.

        Parameters
        ----------
        idx : int
            Sample index in the serial dataset.

        Returns
        -------
        dict[str, Any]
            Distributed feature dictionary with DTensor values.

        """
        features = self.serial_dataset[idx] if self.is_cp_rank_zero else None
        return self._distribute_features(features)


class ValidationDatasetCPWithDTensorV2(_BaseDatasetCPWithDTensorV2):
    """Validation dataset with DTensor context parallelism for Boltz2."""

    def __init__(
        self,
        serial_dataset: torch.utils.data.Dataset,
        device_mesh: DeviceMesh,
        device_mesh_cpu: DeviceMesh,
        val_skip_sample_threshold_tokens: Optional[int] = None,
        val_skip_sample_threshold_atoms: Optional[int] = None,
        val_skip_sample_threshold_seqs: Optional[int] = None,
    ) -> None:
        """Initialize the distributed validation dataset.

        Parameters
        ----------
        serial_dataset : torch.utils.data.Dataset
            The serial (single-rank) Boltz2 validation dataset to wrap.
        device_mesh : DeviceMesh
            Device mesh for distributed tensor operations on GPU.
        device_mesh_cpu : DeviceMesh
            Device mesh for distributed tensor operations on CPU.
        val_skip_sample_threshold_tokens : int, optional
            Skip validation samples with more tokens than this threshold to
            prevent OOM.
        val_skip_sample_threshold_atoms : int, optional
            Skip validation samples with more atoms than this threshold to
            prevent OOM.
        val_skip_sample_threshold_seqs : int, optional
            Skip validation samples with more MSA sequences than this threshold
            to prevent OOM.

        """
        super().__init__(serial_dataset=serial_dataset, device_mesh=device_mesh, device_mesh_cpu=device_mesh_cpu)
        self.val_skip_sample_threshold_tokens = val_skip_sample_threshold_tokens
        self.val_skip_sample_threshold_atoms = val_skip_sample_threshold_atoms
        self.val_skip_sample_threshold_seqs = val_skip_sample_threshold_seqs

    def __getitem__(self, idx: int) -> dict[str, Any]:
        """Fetch and distribute a single validation sample.

        On CP rank zero, iterates from ``idx`` through the dataset looking for
        a sample that satisfies all ``val_skip_sample_threshold_*`` constraints.
        If no sample passes, raises ``RuntimeError``. Other CP ranks receive
        the distributed DTensor features via collectives.

        Parameters
        ----------
        idx : int
            Starting sample index in the serial dataset.

        Returns
        -------
        dict[str, Any]
            Distributed feature dictionary with DTensor values.

        Raises
        ------
        RuntimeError
            If every sample in the dataset is filtered out by the thresholds.

        """
        if self.is_cp_rank_zero:
            num_items = len(self.serial_dataset)
            for shift in range(num_items):
                curr_idx = (idx + shift) % num_items
                features = self.serial_dataset[curr_idx]

                if self.val_skip_sample_threshold_tokens is not None and "token_pad_mask" in features:
                    tokens = int(features["token_pad_mask"].sum().item())
                    if tokens > self.val_skip_sample_threshold_tokens:
                        continue
                if self.val_skip_sample_threshold_atoms is not None and "atom_pad_mask" in features:
                    atoms = int(features["atom_pad_mask"].sum().item())
                    if atoms > self.val_skip_sample_threshold_atoms:
                        continue
                if self.val_skip_sample_threshold_seqs is not None and "msa_mask" in features:
                    msa_mask = features["msa_mask"]
                    seqs = int((msa_mask.sum(dim=1) > 0).sum().item()) if msa_mask.ndim == 2 else int(msa_mask.shape[0])
                    if seqs > self.val_skip_sample_threshold_seqs:
                        continue
                break
            else:
                raise RuntimeError("All validation samples were filtered out by val_skip_sample_threshold_*")
        else:
            features = None

        return self._distribute_features(features)


class Boltz2TrainingDataModule(pl.LightningDataModule):
    """DataModule for Boltz2 distributed training with DTensor CP."""

    def __init__(
        self,
        cfg: DataConfigV2,
        device_mesh: DeviceMesh,
        device_mesh_cpu: DeviceMesh,
    ) -> None:
        """Initialize the distributed training data module.

        Wraps the serial ``Boltz2TrainingDataModule`` and its datasets with
        DTensor context-parallel distribution. Internally creates
        :class:`TrainingDatasetCPWithDTensorV2` and
        :class:`ValidationDatasetCPWithDTensorV2` from the serial module's
        datasets.

        Parameters
        ----------
        cfg : DataConfigV2
            The data configuration.
        device_mesh : DeviceMesh
            Device mesh for distributed tensor operations on GPU.
        device_mesh_cpu : DeviceMesh
            Device mesh for distributed tensor operations on CPU.

        Raises
        ------
        NotImplementedError
            If ``cfg.num_workers != 0``, since multi-worker loading is
            incompatible with DTensor CP collectives in the dataset.

        """
        super().__init__()
        if cfg.num_workers != 0:
            raise NotImplementedError("num_workers != 0 is not supported for CP")

        self.cfg = cfg
        self.device_mesh = device_mesh
        self.device_mesh_cpu = device_mesh_cpu
        self._serial_module = Boltz2TrainingDataModuleSerial(cfg=cfg)
        self.val_group_mapper = self._serial_module.val_group_mapper

        self._train_set = TrainingDatasetCPWithDTensorV2(
            serial_dataset=self._serial_module._train_set,
            device_mesh=self.device_mesh,
            device_mesh_cpu=self.device_mesh_cpu,
        )
        self._val_set = ValidationDatasetCPWithDTensorV2(
            serial_dataset=self._serial_module._val_set,
            device_mesh=self.device_mesh,
            device_mesh_cpu=self.device_mesh_cpu,
            val_skip_sample_threshold_tokens=cfg.val_skip_sample_threshold_tokens,
            val_skip_sample_threshold_atoms=cfg.val_skip_sample_threshold_atoms,
            val_skip_sample_threshold_seqs=cfg.val_skip_sample_threshold_seqs,
        )

    def setup(self, stage: Optional[str] = None) -> None:  # noqa: ARG002
        """Run the setup for the DataModule.

        No-op because the serial module and CP-wrapped datasets are fully
        initialized in ``__init__``.

        Parameters
        ----------
        stage : str, optional
            The stage, one of 'fit', 'validate', 'test'.

        """
        return

    def train_dataloader(self) -> DataLoader:
        """Get the training dataloader.

        Returns
        -------
        DataLoader
            The training dataloader with a ``DistributedSampler`` partitioned
            across data-parallel replicas and a ``CollateDTensor`` collate
            function.

        """
        sampler = DistributedSampler(
            self._train_set,
            num_replicas=self.device_mesh_cpu.shape[0],
            rank=self.device_mesh_cpu.get_local_rank(0),
            shuffle=False,
            drop_last=False,
        )
        custom_collate = CollateDTensor(self.device_mesh_cpu)
        return DataLoader(
            self._train_set,
            sampler=sampler,
            batch_size=self.cfg.batch_size,
            num_workers=self.cfg.num_workers,
            pin_memory=self.cfg.pin_memory,
            shuffle=False,
            collate_fn=custom_collate,
        )

    def val_dataloader(self) -> DataLoader:
        """Get the validation dataloader.

        Returns
        -------
        DataLoader
            The validation dataloader with a ``DistributedSampler`` partitioned
            across data-parallel replicas and a ``CollateDTensor`` collate
            function.

        """
        sampler = DistributedSampler(
            self._val_set,
            num_replicas=self.device_mesh_cpu.shape[0],
            rank=self.device_mesh_cpu.get_local_rank(0),
            shuffle=False,
            drop_last=False,
        )
        custom_collate = CollateDTensor(self.device_mesh_cpu)
        return DataLoader(
            self._val_set,
            sampler=sampler,
            batch_size=self.cfg.val_batch_size,
            num_workers=self.cfg.num_workers,
            pin_memory=self.cfg.pin_memory,
            shuffle=False,
            collate_fn=custom_collate,
        )

    def transfer_batch_to_device(
        self,
        batch: dict,
        device: torch.device,
        dataloader_idx: int,  # noqa: ARG002
    ) -> dict:
        """Transfer a batch from CPU DTensors to the target device.

        DTensor values are moved by extracting the local shard, transferring it
        to ``device``, and re-wrapping with the GPU ``device_mesh``. Plain
        tensors and lists of tensors are transferred directly.

        Parameters
        ----------
        batch : dict
            The batch to transfer.
        device : torch.device
            The target device (typically a CUDA device).
        dataloader_idx : int
            The dataloader index (unused).

        Returns
        -------
        dict
            The batch with all tensor values on ``device``.

        """
        for key, value in batch.items():
            if isinstance(value, DTensor):
                batch_local = value.to_local().to(device)
                batch[key] = DTensor.from_local(
                    batch_local,
                    device_mesh=self.device_mesh,
                    placements=value.placements,
                    shape=value.shape,
                    stride=value.stride(),
                )
            elif isinstance(value, list):
                batch[key] = [item.to(device) if isinstance(item, torch.Tensor) else item for item in value]
            elif isinstance(value, torch.Tensor):
                batch[key] = value.to(device)

        return batch
