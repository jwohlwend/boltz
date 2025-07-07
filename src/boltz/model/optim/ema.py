# --------------------------------------------------------------------------------------
# Modified from Bio-Diffusion (https://github.com/BioinfoMachineLearning/bio-diffusion):
# Modified from : https://github.com/BioinfoMachineLearning/bio-diffusion/blob/main/src/utils/__init__.py
# --------------------------------------------------------------------------------------

from typing import Any, Optional

import torch
from pytorch_lightning import Callback, LightningModule, Trainer
from pytorch_lightning.utilities.exceptions import MisconfigurationException
from pytorch_lightning.utilities.types import STEP_OUTPUT


class EMA(Callback):
    """Implements Exponential Moving Averaging (EMA).

    When training a model, this callback maintains moving averages
    of the trained parameters. When evaluating, we use the moving
    averages copy of the trained parameters. When saving, we save
    an additional set of parameters with the prefix `ema`.

    Adapted from:
    https://github.com/NVIDIA/NeMo/blob/main/nemo/collections/common/callbacks/ema.py
    https://github.com/BioinfoMachineLearning/bio-diffusion/blob/main/src/utils/__init__.py

    """

    def __init__(
        self,
        decay: float = 0.999,
        apply_ema_every_n_steps: int = 1,
        start_step: int = 0,
        eval_with_ema: bool = True,
        warm_start: bool = True,
    ) -> None:
        """Initialize the EMA callback.

        Parameters
        ----------
        decay: float
            The exponential decay, has to be between 0-1.
        apply_ema_every_n_steps: int, optional (default=1)
            Apply EMA every n global steps.
        start_step: int, optional (default=0)
            Start applying EMA from ``start_step`` global step onwards.
        eval_with_ema: bool, optional (default=True)
            Validate the EMA weights instead of the original weights.
            Note this means that when saving the model, the
            validation metrics are calculated with the EMA weights.

        """
        if not (0 <= decay <= 1):
            msg = "EMA decay value must be between 0 and 1"
            raise MisconfigurationException(msg)

        self._ema_weights: Optional[dict[str, torch.Tensor]] = None
        self._cur_step: Optional[int] = None
        self._weights_buffer: Optional[dict[str, torch.Tensor]] = None
        self.apply_ema_every_n_steps = apply_ema_every_n_steps
        self.start_step = start_step
        self.eval_with_ema = eval_with_ema
        self.decay = decay
        self.warm_start = warm_start

    @property
    def ema_initialized(self) -> bool:
        """Check if EMA weights have been initialized.

        Returns
        -------
        bool
            Whether the EMA weights have been initialized.

        """
        return self._ema_weights is not None

    def state_dict(self) -> dict[str, Any]:
        """Return the current state of the callback.

        Returns
        -------
        dict[str, Any]
            The current state of the callback.

        """
        return {
            "cur_step": self._cur_step,
            "ema_weights": self._ema_weights,
        }

    def load_state_dict(self, state_dict: dict[str, Any]) -> None:
        """Load the state of the callback.

        Parameters
        ----------
        state_dict: dict[str, Any]
            The state of the callback to load.

        """
        self._cur_step = state_dict["cur_step"]
        self._ema_weights = state_dict["ema_weights"]

    def should_apply_ema(self, step: int) -> bool:
        """Check if EMA should be applied at the current step.

        Parameters
        ----------
        step: int
            The current global step.

        Returns
        -------
        bool
            True if EMA should be applied, False otherwise.

        """
        return (
            step != self._cur_step
            and step >= self.start_step
            and step % self.apply_ema_every_n_steps == 0
        )

    def apply_ema(self, pl_module: LightningModule) -> None:
        """Apply EMA to the model weights.

        Parameters
        ----------
        pl_module: LightningModule
            The LightningModule instance.

        """
        decay = self.decay
        if self.warm_start:
            decay = min(decay, (1 + self._cur_step) / (10 + self._cur_step))

        for k, orig_weight in pl_module.state_dict().items():
            ema_weight = self._ema_weights[k]
            if (
                ema_weight.data.dtype != torch.long  # noqa: PLR1714
                and orig_weight.data.dtype != torch.long  # skip non-trainable weights
            ):
                diff = ema_weight.data - orig_weight.data
                diff.mul_(1.0 - decay)
                ema_weight.sub_(diff)

    def on_load_checkpoint(
        self,
        trainer: Trainer,
        pl_module: LightningModule,  # noqa: ARG002
        checkpoint: dict[str, Any],
    ) -> None:
        """Load the EMA weights from the checkpoint.

        Parameters
        ----------
        trainer: Trainer
            The Trainer instance.
        pl_module: LightningModule
            The LightningModule instance.
        checkpoint: dict[str, Any]
            The checkpoint to load.

        """
        if "ema" in checkpoint:
            print("LOADING CHECKPOINT RUNNING")
            self.load_state_dict(checkpoint["ema"])

    def on_save_checkpoint(
        self,
        trainer: Trainer,  # noqa: ARG002
        pl_module: LightningModule,  # noqa: ARG002
        checkpoint: dict[str, Any],
    ) -> None:
        """Save the EMA weights to the checkpoint.

        Parameters
        ----------
        trainer: Trainer
            The Trainer instance.
        pl_module: LightningModule
            The LightningModule instance.
        checkpoint: dict[str, Any]
            The checkpoint to save.

        """
        if self.ema_initialized:
            checkpoint["ema"] = self.state_dict()

    def on_train_start(self, trainer: Trainer, pl_module: LightningModule) -> None:  # noqa: ARG002
        """Initialize EMA weights and move to device.

        Parameters
        ----------
        trainer: pl.Trainer
            The Trainer instance.
        pl_module: pl.LightningModule
            The LightningModule instance.

        """
        # Create EMA weights if not already initialized
        if not self.ema_initialized:
            self._ema_weights = {
                k: p.detach().clone() for k, p in pl_module.state_dict().items()
            }

        # Move EMA weights to the correct device
        self._ema_weights = {
            k: p.to(pl_module.device) for k, p in self._ema_weights.items()
        }

    def on_train_batch_end(
        self,
        trainer: Trainer,
        pl_module: LightningModule,
        outputs: STEP_OUTPUT,  # noqa: ARG002
        batch: Any,  # noqa: ARG002
        batch_idx: int,  # noqa: ARG002
    ) -> None:
        """Apply EMA to the model weights at the end of each training batch.

        Parameters
        ----------
        trainer: Trainer
            The Trainer instance.
        pl_module: LightningModule
            The LightningModule instance.
        outputs: STEP_OUTPUT
            The outputs of the model.
        batch: Any
            The current batch.
        batch_idx: int
            The index of the current batch.

        """
        if self.should_apply_ema(trainer.global_step):
            self._cur_step = trainer.global_step
            self.apply_ema(pl_module)

    def replace_model_weights(self, pl_module: LightningModule) -> None:
        """Replace model weights with EMA weights.

        Parameters
        ----------
        pl_module: LightningModule
            The LightningModule instance.

        """
        self._weights_buffer = {
            k: p.detach().clone().to("cpu") for k, p in pl_module.state_dict().items()
        }
        pl_module.load_state_dict(self._ema_weights, strict=False)

    def restore_original_weights(self, pl_module: LightningModule) -> None:
        """Restore model weights to original weights.

        Parameters
        ----------
        pl_module: LightningModule
            The LightningModule instance.

        """
        pl_module.load_state_dict(self._weights_buffer, strict=False)
        del self._weights_buffer

    def _on_eval_start(self, pl_module: LightningModule) -> None:
        """Use EMA weights for evaluation.

        Parameters
        ----------
        pl_module: LightningModule
            The LightningModule instance.

        """
        if self.ema_initialized and self.eval_with_ema:
            self.replace_model_weights(pl_module)

    def _on_eval_end(self, pl_module: LightningModule) -> None:
        """Restore original weights after evaluation.

        Parameters
        ----------
        pl_module: LightningModule
            The LightningModule instance.

        """
        if self.ema_initialized and self.eval_with_ema:
            self.restore_original_weights(pl_module)

    def on_validation_start(
        self,
        trainer: Trainer,  # noqa: ARG002
        pl_module: LightningModule,
    ) -> None:
        """Use EMA weights for validation.

        Parameters
        ----------
        trainer: Trainer
            The Trainer instance.
        pl_module: LightningModule
            The LightningModule instance.

        """
        self._on_eval_start(pl_module)

    def on_validation_end(
        self,
        trainer: Trainer,  # noqa: ARG002
        pl_module: LightningModule,
    ) -> None:
        """Restore original weights after validation.

        Parameters
        ----------
        trainer: Trainer
            The Trainer instance.
        pl_module: LightningModule
            The LightningModule instance.

        """
        self._on_eval_end(pl_module)

    def on_test_start(
        self,
        trainer: Trainer,  # noqa: ARG002
        pl_module: LightningModule,
    ) -> None:
        """Use EMA weights for testing.

        Parameters
        ----------
        trainer: Trainer
            The Trainer instance.
        pl_module: LightningModule
            The LightningModule instance.

        """
        self._on_eval_start(pl_module)

    def on_test_end(
        self,
        trainer: Trainer,  # noqa: ARG002
        pl_module: LightningModule,
    ) -> None:
        """Restore original weights after testing.

        Parameters
        ----------
        trainer: Trainer
            The Trainer instance.
        pl_module: LightningModule
            The LightningModule instance.

        """
        self._on_eval_end(pl_module)

    def on_predict_start(
        self,
        trainer: Trainer,  # noqa: ARG002
        pl_module: LightningModule,
    ) -> None:
        """Use EMA weights for prediction.

        Parameters
        ----------
        trainer: Trainer
            The Trainer instance.
        pl_module: LightningModule
            The LightningModule instance.

        """
        self._on_eval_start(pl_module)

    def on_predict_end(
        self,
        trainer: Trainer,  # noqa: ARG002
        pl_module: LightningModule,
    ) -> None:
        """Restore original weights after prediction.

        Parameters
        ----------
        trainer: Trainer
            The Trainer instance.
        pl_module: LightningModule
            The LightningModule instance.

        """
        self._on_eval_end(pl_module)
