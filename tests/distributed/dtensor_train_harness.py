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

"""Shared harness for lightweight DTensor distributed training tests."""

from __future__ import annotations

from pathlib import Path

import pytorch_lightning as pl
import torch
from torch.distributed.tensor import Replicate, Shard, distribute_tensor
from torch.utils.data import DataLoader, Dataset

from boltz.distributed.comm import TransposeComm
from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.loss.distogram import distogram_loss as distogram_loss_dtensor
from boltz.distributed.model.modules.trunkv2 import DistogramModule as DistogramModuleDTensor
from boltz.distributed.model.optim.ema import DistributedEMA
from boltz.model.loss.distogramv2 import distogram_loss as distogram_loss_serial
from boltz.model.modules.trunkv2 import DistogramModule as SerialDistogramModule
from boltz.testing.utils import init_module_params_uniform


class SyntheticDistogramDataset(Dataset):
    """Deterministic synthetic dataset for CP distogram training tests."""

    def __init__(
        self,
        *,
        seq_len: int,
        token_z: int,
        num_bins: int,
        num_conformers: int,
        num_samples: int,
        seed: int,
    ) -> None:
        super().__init__()
        generator = torch.Generator()
        generator.manual_seed(seed)

        samples: list[dict[str, torch.Tensor]] = []
        for _ in range(num_samples):
            z = torch.rand((seq_len, seq_len, token_z), generator=generator)
            target = torch.softmax(
                torch.randn((seq_len, seq_len, num_conformers, num_bins), generator=generator),
                dim=-1,
            )
            mask = torch.randint(0, 2, (seq_len,), generator=generator, dtype=torch.bool)
            samples.append({"z": z, "target": target, "mask": mask})
        self.samples = samples

    def __len__(self) -> int:
        return len(self.samples)

    def __getitem__(self, index: int) -> dict[str, torch.Tensor]:
        return self.samples[index]


class TinyDistogramCPModel(pl.LightningModule):
    """Tiny LightningModule that exercises DTensor distogram forward/loss in CP."""

    def __init__(
        self,
        *,
        dist_manager: DistributedManager,
        token_z: int,
        num_bins: int,
        num_distograms: int,
        num_conformers: int,
        serial_state_dict: dict[str, torch.Tensor],
        learning_rate: float,
        adam_betas: tuple[float, float] = (0.9, 0.999),
    ) -> None:
        super().__init__()
        self.dist_manager = dist_manager
        self.num_conformers = num_conformers
        self.learning_rate = learning_rate
        self.adam_betas = adam_betas
        self.last_pred_local: torch.Tensor | None = None

        serial_module = SerialDistogramModule(
            token_z=token_z,
            num_bins=num_bins,
            num_distograms=num_distograms,
        )
        serial_module.load_state_dict(serial_state_dict)
        serial_module.to(dist_manager.device)  # Must be on mesh device before DTensor wrapping
        self.loss_comm = TransposeComm(dist_manager.group["cp"], dist_manager.layout_subgroups["cp"])
        self.distogram_module = DistogramModuleDTensor(
            module=serial_module,
            dist_manager=dist_manager,
            distogram_comm=self.loss_comm,
        )

    def on_train_start(self) -> None:
        # Loss log keyed by (epoch, batch_idx) for stop/go trajectory comparison.
        self._loss_log: dict[tuple[int, int], float] = {}

    def _assert_sharding_active(
        self,
        global_tensor: torch.Tensor,
        dtensor: torch.Tensor,
        label: str,
    ) -> None:
        """Assert local shape is smaller than global for each Shard dim.

        Catches bugs where CP silently falls back to replication.
        """
        from torch.distributed.tensor import DTensor as _DTensor

        if not isinstance(dtensor, _DTensor):
            raise AssertionError(f"{label}: expected DTensor, got {type(dtensor)}")

        local = dtensor.to_local()
        for mesh_dim, placement in enumerate(dtensor.placements):
            if isinstance(placement, Shard):
                shard_dim = placement.dim
                mesh_size = dtensor.device_mesh.size(mesh_dim)
                if mesh_size > 1:
                    assert local.shape[shard_dim] < global_tensor.shape[shard_dim], (
                        f"{label}: Shard({shard_dim}) on mesh dim {mesh_dim} (size={mesh_size}) "
                        f"did not reduce local shape — global {global_tensor.shape[shard_dim]}, "
                        f"local {local.shape[shard_dim]}.  CP may not be active."
                    )

    def training_step(self, batch: dict[str, torch.Tensor], batch_idx: int) -> torch.Tensor:
        z = batch["z"].to(self.dist_manager.device)
        target = batch["target"].to(self.dist_manager.device)
        mask = batch["mask"].to(self.dist_manager.device)

        # Mesh dims: (dp, cp_row, cp_col).  Shard(0) splits the batch across
        # dp ranks; DataLoader supplies global batch (micro_batch * dp_size).
        z_dtensor = distribute_tensor(
            z,
            device_mesh=self.dist_manager.device_mesh_subgroups,
            placements=(Shard(0), Shard(1), Shard(2)),
        )
        target_dtensor = distribute_tensor(
            target,
            device_mesh=self.dist_manager.device_mesh_subgroups,
            placements=(Shard(0), Shard(1), Shard(2)),
        )
        mask_dtensor = distribute_tensor(
            mask,
            device_mesh=self.dist_manager.device_mesh_subgroups,
            placements=(Shard(0), Shard(1), Replicate()),
        )

        # Verify DTensor sharding is active — proves CP is real, not silently
        # replicated.  Each mesh dim shrinks the corresponding tensor dim.
        self._assert_sharding_active(z, z_dtensor, "z")

        pred_dtensor = self.distogram_module(z_dtensor)
        # Capture local output for tests without altering Trainer internals.
        self.last_pred_local = pred_dtensor.to_local().detach().cpu()
        global_loss_dtensor, _ = distogram_loss_dtensor(
            {"pdistogram": pred_dtensor},
            {"disto_target": target_dtensor, "token_disto_mask": mask_dtensor},
            self.loss_comm,
            aggregate_distogram=False,
        )
        loss = global_loss_dtensor.to_local()
        self.log("train/loss", loss, on_step=True, on_epoch=True, prog_bar=False)
        self._loss_log[(self.current_epoch, batch_idx)] = loss.detach().item()
        return loss

    def configure_optimizers(self):
        return torch.optim.Adam(self.parameters(), lr=self.learning_rate, betas=self.adam_betas)


class TinyDistogramCPModelWithEMA(TinyDistogramCPModel):
    """TinyDistogramCPModel with ``DistributedEMA`` via ``configure_callbacks()``."""

    def __init__(
        self,
        *,
        ema_decay: float = 0.999,
        **kwargs: object,
    ) -> None:
        super().__init__(**kwargs)  # type: ignore[arg-type]
        self.ema_decay = ema_decay

    def configure_callbacks(self) -> list[pl.Callback]:
        return [DistributedEMA(decay=self.ema_decay)]


class TinyDistogramSerialModel(pl.LightningModule):
    """Serial (plain-tensor) counterpart of ``TinyDistogramCPModel``.

    Same loss, hyperparameters, and state-dict keys — used for cross-mode
    checkpoint tests (serial ↔ distributed).
    """

    def __init__(
        self,
        *,
        token_z: int,
        num_bins: int,
        num_distograms: int,
        num_conformers: int,
        serial_state_dict: dict[str, torch.Tensor],
        learning_rate: float,
        adam_betas: tuple[float, float] = (0.9, 0.999),
    ) -> None:
        super().__init__()
        self.num_conformers = num_conformers
        self.learning_rate = learning_rate
        self.adam_betas = adam_betas

        self.distogram_module = SerialDistogramModule(
            token_z=token_z,
            num_bins=num_bins,
            num_distograms=num_distograms,
        )
        self.distogram_module.load_state_dict(serial_state_dict)

    def on_train_start(self) -> None:
        self._loss_log: dict[tuple[int, int], float] = {}

    def training_step(self, batch: dict[str, torch.Tensor], batch_idx: int) -> torch.Tensor:
        z = batch["z"].to(self.device)
        target = batch["target"].to(self.device)
        mask = batch["mask"].to(self.device)

        pred = self.distogram_module(z)
        global_loss, _ = distogram_loss_serial(
            {"pdistogram": pred},
            {"disto_target": target, "token_disto_mask": mask},
            aggregate_distogram=False,
        )
        self.log("train/loss", global_loss, on_step=True, on_epoch=True, prog_bar=False)
        self._loss_log[(self.current_epoch, batch_idx)] = global_loss.detach().item()
        return global_loss

    def configure_optimizers(self):
        return torch.optim.Adam(self.parameters(), lr=self.learning_rate, betas=self.adam_betas)


class DistogramTrainDataModule(pl.LightningDataModule):
    """LightningDataModule wrapping :class:`SyntheticDistogramDataset`.

    Used by stop/go tests that route through ``train.py``, which requires a
    ``LightningDataModule`` (not a bare ``DataLoader``).
    """

    def __init__(
        self,
        *,
        seq_len: int,
        token_z: int,
        num_bins: int,
        num_conformers: int,
        num_samples: int,
        seed: int,
        batch_size: int = 1,
        dp_size: int = 1,
    ) -> None:
        super().__init__()
        self._seq_len = seq_len
        self._token_z = token_z
        self._num_bins = num_bins
        self._num_conformers = num_conformers
        self._num_samples = num_samples
        self._seed = seed
        self._batch_size = batch_size
        self._dp_size = dp_size

    def train_dataloader(self) -> DataLoader:
        return create_train_dataloader(
            seq_len=self._seq_len,
            token_z=self._token_z,
            num_bins=self._num_bins,
            num_conformers=self._num_conformers,
            num_samples=self._num_samples,
            seed=self._seed,
            batch_size=self._batch_size,
            dp_size=self._dp_size,
        )


def create_initial_serial_state_dict(
    *,
    token_z: int,
    num_bins: int,
    num_distograms: int,
    seed: int,
) -> dict[str, torch.Tensor]:
    """Create deterministic initial serial distogram parameters."""
    with torch.random.fork_rng(devices=[]):
        torch.manual_seed(seed)
        serial_module = SerialDistogramModule(
            token_z=token_z,
            num_bins=num_bins,
            num_distograms=num_distograms,
        )
        init_module_params_uniform(serial_module, low=-0.25, high=0.25)
        return {key: value.detach().clone().cpu() for key, value in serial_module.state_dict().items()}


def create_train_dataloader(
    *,
    seq_len: int,
    token_z: int,
    num_bins: int,
    num_conformers: int,
    num_samples: int,
    seed: int,
    batch_size: int = 1,
    dp_size: int = 1,
) -> DataLoader:
    """Create deterministic dataloader for distributed train tests.

    Produces global batches of ``batch_size * dp_size``; the test harness
    splits them across dp ranks via ``distribute_tensor(Shard(0))``.
    """
    global_batch = batch_size * dp_size
    dataset = SyntheticDistogramDataset(
        seq_len=seq_len,
        token_z=token_z,
        num_bins=num_bins,
        num_conformers=num_conformers,
        num_samples=num_samples,
        seed=seed,
    )
    return DataLoader(dataset, batch_size=global_batch, shuffle=False, num_workers=0)


def state_dicts_differ(
    baseline_state: dict[str, torch.Tensor],
    candidate_state: dict[str, torch.Tensor],
) -> bool:
    """Return True if at least one key differs."""
    for key in baseline_state:
        if key not in candidate_state:
            return True
        if not torch.equal(baseline_state[key], candidate_state[key].cpu()):
            return True
    return False


def save_rank_state_dict(output_dir: Path, rank: int, state_dict: dict[str, torch.Tensor]) -> Path:
    """Save rank-local serialized state dict for cross-rank checks."""
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / f"rank_{rank}_state.pt"
    torch.save(state_dict, path)
    return path
