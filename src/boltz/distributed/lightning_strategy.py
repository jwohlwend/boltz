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

"""Lightning strategy helpers for DTensor context-parallel training."""

import logging
from pathlib import Path
from typing import Any, Mapping, Optional

import torch
from pytorch_lightning.core.optimizer import LightningOptimizer
from pytorch_lightning.strategies import SingleDeviceStrategy
from torch.distributed.checkpoint.state_dict import (
    StateDictOptions,
    get_optimizer_state_dict,
    set_optimizer_state_dict,
)
from torch.distributed.tensor import DTensor, Replicate
from typing_extensions import override

from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.modules.utils import (
    convert_distributed_checkpoint_to_serial_state_dict,
    convert_dtensors_to_tensors,
    convert_serial_checkpoint_to_distributed_state_dict,
)

logger = logging.getLogger(__name__)


def _redistribute_optimizer_state_to_params(optimizer: torch.optim.Optimizer) -> None:
    """Redistribute plain-tensor optimizer state to match DTensor parameters.

    After loading a checkpoint, optimizer state buffers (``exp_avg``,
    ``exp_avg_sq``, etc.) are plain tensors.  If the corresponding parameter
    is a DTensor, the state must be re-distributed to the same mesh and
    placements to avoid mixed-type errors on ``optimizer.step()``.
    """
    for group in optimizer.param_groups:
        for param in group["params"]:
            if not isinstance(param, DTensor):
                continue
            param_state = optimizer.state.get(param)
            if param_state is None:
                continue
            if not all(isinstance(p, Replicate) for p in param.placements):
                raise ValueError(
                    f"Only Replicate placements are supported for optimizer state redistribution, "
                    f"got {param.placements}"
                )
            for state_key, state_val in param_state.items():
                if isinstance(state_val, torch.Tensor) and not isinstance(state_val, DTensor):
                    state_val = state_val.to(device=param.device_mesh.device_type)
                    # All ranks load the same checkpoint, so state_val is
                    # already identical across ranks.  from_local avoids the
                    # redundant all-gather that distribute_tensor would do.
                    param_state[state_key] = DTensor.from_local(
                        state_val,
                        device_mesh=param.device_mesh,
                        placements=param.placements,
                        shape=state_val.shape,
                        stride=state_val.stride(),
                    )


class BoltzContextParallelStrategy(SingleDeviceStrategy):
    """DTensor-aware strategy for context-parallel checkpoint handling.

    This strategy intentionally stays close to single-device semantics while
    customizing checkpoint save/load behavior:

    - save: convert model DTensors to regular tensors for portable checkpoints
    - load: map serial checkpoint tensors into the currently-materialized model
      state_dict template (which may include DTensor parameters/buffers)
    """

    strategy_name = "boltz_context_parallel"

    def __init__(self, dist_manager: DistributedManager, *args: Any, **kwargs: Any) -> None:
        self.dist_manager = dist_manager
        super().__init__(*args, device=self.dist_manager.device, **kwargs)
        self.global_rank = self.dist_manager.rank
        self.local_rank = self.dist_manager.local_rank
        self.world_size = self.dist_manager.world_size

    @property
    @override
    def is_global_zero(self) -> bool:
        return self.global_rank == 0

    @override
    def barrier(self, *args: Any, **kwargs: Any) -> None:
        if torch.distributed.is_initialized():
            torch.distributed.barrier()

    @override
    def save_checkpoint(
        self, checkpoint: dict[str, Any], filepath: str | Path, storage_options: Optional[Any] = None
    ) -> None:
        """Save checkpoint with all DTensors converted to plain tensors.

        This ensures *every* DTensor in the checkpoint (model state_dict,
        optimizer states, EMA shadow weights, etc.) is stored in a portable
        serial format, not just the ``state_dict`` key.
        """
        checkpoint = convert_dtensors_to_tensors(checkpoint)
        super().save_checkpoint(checkpoint, filepath, storage_options=storage_options)

    @override
    def load_checkpoint(self, checkpoint_path: str | Path) -> dict[str, Any]:
        # Route through Lightning's checkpoint_io so strategy-level checkpoint
        # backends and remapping behavior are respected.
        # Use map_location="cpu" (string, not callable) to avoid loading onto a
        # serialized CUDA device that may not exist on the current node.
        return self.checkpoint_io.load_checkpoint(checkpoint_path, map_location="cpu")

    @override
    def model_to_device(self) -> None:
        # Distributed models are expected to be materialized on the target device
        # before strategy setup; calling .to(...) here can cause context issues.
        return None

    @override
    def lightning_module_state_dict(self) -> dict[str, Any]:
        if self.model is None:
            raise RuntimeError(
                "BoltzContextParallelStrategy.model is not set. "
                "Attach the model before exporting a checkpoint state_dict."
            )
        checkpoint = {"state_dict": self.model.state_dict()}
        return convert_distributed_checkpoint_to_serial_state_dict(checkpoint)

    @override
    def load_model_state_dict(self, checkpoint: Mapping[str, Any], strict: bool = False) -> None:
        if self.lightning_module is None:
            raise RuntimeError(
                "BoltzContextParallelStrategy.lightning_module is not set. "
                "Attach the LightningModule before loading model state."
            )
        state_template = self.lightning_module.state_dict()
        distributed_state_dict = convert_serial_checkpoint_to_distributed_state_dict(
            checkpoint=checkpoint,
            strict=strict,
            state_dict_template=state_template,
        )
        self.lightning_module.load_state_dict(distributed_state_dict, strict=strict)

    @override
    def optimizer_state(self, optimizer: Any) -> dict[str, Any]:
        """Return optimizer state dict with FQN (fully qualified name) keys.

        Uses ``torch.distributed.checkpoint.state_dict.get_optimizer_state_dict``
        to produce parameter-name-keyed state (e.g. ``"distogram_module.weight"``
        instead of ``0``), making checkpoints portable across model topologies
        regardless of parameter registration order.
        """
        if isinstance(optimizer, LightningOptimizer):
            optimizer = optimizer._optimizer
        if self.lightning_module is None:
            raise RuntimeError(
                "BoltzContextParallelStrategy.lightning_module is not set. "
                "Attach the LightningModule before saving optimizer state."
            )
        return get_optimizer_state_dict(self.lightning_module, optimizer)

    @override
    def load_optimizer_state_dict(self, checkpoint: Mapping[str, Any]) -> None:
        """Load optimizer state, handling both FQN-keyed and legacy int-keyed formats.

        Checkpoints saved by this strategy use FQN (fully qualified name) string
        keys produced by ``get_optimizer_state_dict``.  Legacy checkpoints (or
        those saved by the default Lightning strategy) use integer keys.

        This method auto-detects the format by inspecting the first key in
        ``state_dict["state"]``:

        - **FQN keys (str):** loaded via ``set_optimizer_state_dict`` which maps
          parameter names back to the live optimizer's internal indices and
          redistributes plain tensors to match DTensor parameter placements.
        - **Integer keys (int):** loaded via ``optimizer.load_state_dict()``
          followed by manual DTensor redistribution (legacy path).
        """
        if not self.optimizers:
            return

        if "optimizer_states" not in checkpoint:
            raise ValueError("Checkpoint is passed into load_optimizer_state_dict but no optimizer_states found")

        optimizer_states = checkpoint["optimizer_states"]
        if not isinstance(optimizer_states, (list, tuple)):
            raise TypeError(f"Checkpoint field 'optimizer_states' must be a list/tuple, got {type(optimizer_states)}")
        if len(optimizer_states) != len(self.optimizers):
            raise ValueError(
                "Optimizer-state length mismatch: "
                f"checkpoint has {len(optimizer_states)} state entries but strategy has {len(self.optimizers)} optimizers"
            )

        for index, optimizer in enumerate(self.optimizers):
            optimizer_state = optimizer_states[index]
            if isinstance(optimizer, LightningOptimizer):
                optimizer = optimizer._optimizer

            state_keys = list(optimizer_state.get("state", {}).keys())
            uses_fqn_keys = state_keys and isinstance(state_keys[0], str)

            if uses_fqn_keys:
                if self.lightning_module is None:
                    raise RuntimeError(
                        "BoltzContextParallelStrategy.lightning_module is not set. "
                        "Attach the LightningModule before loading optimizer state."
                    )
                set_optimizer_state_dict(
                    self.lightning_module,
                    optimizer,
                    optim_state_dict=optimizer_state,
                    options=StateDictOptions(full_state_dict=True),
                )
            else:
                logger.warning(
                    "Loading optimizer state with legacy integer keys. "
                    "This is expected when resuming from a checkpoint saved by "
                    "the default Lightning strategy or an older version of the "
                    "distributed trainer. Future checkpoints will use FQN keys."
                )
                optimizer.load_state_dict(optimizer_state)

            # Ensure optimizer state buffers match DTensor parameter placements.
            # For the FQN path, set_optimizer_state_dict may already redistribute;
            # this is a safe no-op when state is already distributed correctly.
            _redistribute_optimizer_state_to_params(optimizer)
