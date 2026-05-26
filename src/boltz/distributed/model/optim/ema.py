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

"""DTensor-aware EMA callback for distributed context-parallel training.

This module provides :class:`DistributedEMA`, a drop-in subclass of the
base :class:`~boltz.model.optim.ema.EMA` callback that handles the
DTensor ↔ plain-tensor conversions required when model parameters live
on a distributed device mesh.

Key differences from the base EMA callback:

* **Save** – EMA shadow weights (which are DTensors during training) are
  converted to plain local tensors before being written to the checkpoint,
  keeping checkpoints portable across topologies.
* **Load** – Plain-tensor EMA weights from a checkpoint are re-distributed
  to match the model's current DTensor placements before the first
  training step.
* **Weight swap** (for validation with EMA) – The model-weight backup
  preserves DTensor metadata instead of moving to CPU (which strips
  DTensor placement information).
"""

from __future__ import annotations

import logging
from typing import Any

import torch
from pytorch_lightning import LightningModule, Trainer
from torch.distributed.tensor import DTensor, Replicate

from boltz.distributed.model.modules.utils import convert_dtensors_to_tensors
from boltz.model.optim.ema import EMA

logger = logging.getLogger(__name__)


def _assert_replicate_only(tensor: DTensor, name: str) -> None:
    """Raise if *tensor* has any non-Replicate placements.

    EMA uses in-place ``.data`` arithmetic which bypasses DTensor dispatch,
    so it is only correct when every placement is ``Replicate``.
    """
    for placement in tensor.placements:
        if not isinstance(placement, Replicate):
            raise ValueError(
                f"DistributedEMA requires all placements to be Replicate for "
                f"correct in-place arithmetic, but parameter '{name}' has "
                f"placement {placement!r}.  Shard/Partial placements would "
                f"silently produce incorrect EMA updates."
            )


class DistributedEMA(EMA):
    """DTensor-aware Exponential Moving Average callback.

    For models whose parameters are plain :class:`torch.Tensor` objects this
    callback is behaviourally identical to :class:`EMA`.  When parameters are
    :class:`~torch.distributed.tensor.DTensor` instances (e.g. under context
    parallelism), the additional conversions described in the module docstring
    are applied transparently.
    """

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        super().__init__(*args, **kwargs)
        self._placements_validated: bool = False

    @staticmethod
    def _realign_weights_to_model(
        weights: dict[str, torch.Tensor],
        pl_module: LightningModule,
    ) -> dict[str, torch.Tensor]:
        """Re-distribute plain-tensor weights to match a model's DTensor layout.

        For every key in *weights* whose corresponding model parameter is a
        :class:`DTensor`, the plain tensor is wrapped with the same device-mesh
        and placements.  Keys that are already DTensors or whose templates are
        plain tensors pass through unchanged.

        Args:
            weights: Dict mapping parameter names to tensors (plain or DTensor).
            pl_module: The LightningModule whose ``state_dict()`` defines the
                target DTensor layout.

        Returns:
            A new dict with values distributed to match the model layout.

        Raises:
            ValueError: If an EMA weight shape does not match the model
                template shape, or if a DTensor template has non-Replicate
                placements.
        """
        model_state = pl_module.state_dict()
        ema_keys = set(weights.keys())
        model_keys = set(model_state.keys())
        missing_in_model = ema_keys - model_keys
        missing_in_ema = model_keys - ema_keys
        if missing_in_model:
            logger.warning(
                "DistributedEMA: %d EMA key(s) not found in model state_dict (stale checkpoint?): %s",
                len(missing_in_model),
                sorted(missing_in_model)[:5],
            )
        if missing_in_ema:
            logger.warning(
                "DistributedEMA: %d model key(s) not found in EMA weights (new layers added?): %s",
                len(missing_in_ema),
                sorted(missing_in_ema)[:5],
            )

        realigned: dict[str, torch.Tensor] = {}
        for key, weight in weights.items():
            template = model_state.get(key)
            if template is not None and isinstance(template, DTensor) and not isinstance(weight, DTensor):
                if tuple(weight.shape) != tuple(template.shape):
                    raise ValueError(
                        f"DistributedEMA: shape mismatch for '{key}': "
                        f"EMA weight has shape {tuple(weight.shape)} but model "
                        f"template has shape {tuple(template.shape)}"
                    )
                _assert_replicate_only(template, key)

                weight = weight.to(device=template.device_mesh.device_type, dtype=template.dtype)
                # Placements are validated as Replicate above, and all
                # ranks hold the same EMA checkpoint — from_local avoids
                # the redundant all-gather that distribute_tensor would do.
                weight = DTensor.from_local(
                    weight,
                    device_mesh=template.device_mesh,
                    placements=template.placements,
                    shape=template.shape,
                    stride=template.stride(),
                )
            realigned[key] = weight

        for key in missing_in_ema:
            realigned[key] = model_state[key].detach().clone()

        return realigned

    def apply_ema(self, pl_module: LightningModule) -> None:
        """Apply EMA update, asserting Replicate placements on DTensors.

        The base EMA class performs arithmetic via ``.data`` which bypasses
        DTensor dispatch.  This is only correct when all placements are
        ``Replicate``.  This override validates that invariant on the first
        call and then delegates to the base implementation.
        """
        # Validate on first call only (placements don't change during training)
        if not self._placements_validated:
            for k, ema_w in self._ema_weights.items():
                if isinstance(ema_w, DTensor):
                    _assert_replicate_only(ema_w, f"ema[{k}]")
            for k, param in pl_module.state_dict().items():
                if isinstance(param, DTensor):
                    _assert_replicate_only(param, f"model[{k}]")
            self._placements_validated = True

        super().apply_ema(pl_module)

    def on_save_checkpoint(
        self,
        trainer: Trainer,  # noqa: ARG002
        pl_module: LightningModule,  # noqa: ARG002
        checkpoint: dict[str, Any],
    ) -> None:
        """Save EMA state, converting any DTensors to plain tensors."""
        if self.ema_initialized:
            checkpoint["ema"] = {
                "cur_step": self._cur_step,
                "ema_weights": convert_dtensors_to_tensors(self._ema_weights),
            }

    def on_train_start(self, trainer: Trainer, pl_module: LightningModule) -> None:  # noqa: ARG002
        """Initialise or re-distribute EMA weights to match the model layout.

        * **First run** – clones the model's ``state_dict()`` (which may
          contain DTensors).
        * **Resumed run** – plain-tensor EMA weights loaded from checkpoint
          are redistributed to match the model's DTensor placements.
        """
        # Re-validate placements on next apply_ema (topology may have changed).
        self._placements_validated = False

        if not self.ema_initialized:
            self._ema_weights = {k: p.detach().clone() for k, p in pl_module.state_dict().items()}
        else:
            # Re-distribute plain-tensor EMA weights to match the model's DTensor layout.
            self._ema_weights = self._realign_weights_to_model(self._ema_weights, pl_module)

        # Move to correct device (preserves DTensor placements).
        self._ema_weights = {k: p.to(pl_module.device) for k, p in self._ema_weights.items()}

    def replace_model_weights(self, pl_module: LightningModule) -> None:
        """Replace model weights with EMA weights, backing up originals to CPU.

        DTensor metadata (device mesh and placements) is stored alongside the
        CPU tensors so :meth:`restore_original_weights` can reconstruct the
        original DTensor layout without an extra ``state_dict()`` call.
        Keeping the backup on CPU instead of on-device frees GPU memory for
        the duration of validation.

        The ``inference_mode(False)`` guard is required because Lightning's
        ``trainer.predict()`` wraps the entire workflow in
        ``torch.inference_mode(True)`` by default (since Lightning ≥2.x).
        Operations like ``.detach()``, ``.clone()``, and ``load_state_dict``
        manipulate version counters, which PyTorch ≥2.10 disallows on
        inference tensors.
        """
        with torch.inference_mode(False):
            self._weights_buffer: dict[str, torch.Tensor] = {}
            self._weights_dtensor_meta: dict[str, tuple] = {}
            for k, p in pl_module.state_dict().items():
                if isinstance(p, DTensor):
                    self._weights_dtensor_meta[k] = (p.device_mesh, p.placements, p.shape, p.stride())
                    self._weights_buffer[k] = p.to_local().detach().clone().cpu()
                else:
                    self._weights_buffer[k] = p.detach().clone().cpu()
            ema_weights_to_load = self._realign_weights_to_model(self._ema_weights, pl_module)
            pl_module.load_state_dict(ema_weights_to_load, strict=False)

    def restore_original_weights(self, pl_module: LightningModule) -> None:
        """Restore original model weights from the CPU backup.

        For keys that were originally DTensors, the plain CPU tensor is
        redistributed back to the stored device mesh and placements.

        See :meth:`replace_model_weights` for why the ``inference_mode(False)``
        guard is necessary.
        """
        with torch.inference_mode(False):
            restored: dict[str, torch.Tensor] = {}
            for k, cpu_tensor in self._weights_buffer.items():
                if k in self._weights_dtensor_meta:
                    mesh, placements, shape, stride = self._weights_dtensor_meta[k]
                    restored[k] = DTensor.from_local(
                        cpu_tensor.to(device=mesh.device_type),
                        device_mesh=mesh,
                        placements=placements,
                        shape=shape,
                        stride=stride,
                    )
                else:
                    restored[k] = cpu_tensor.to(device=pl_module.device)
            pl_module.load_state_dict(restored, strict=False)
            del self._weights_buffer
            del self._weights_dtensor_meta
