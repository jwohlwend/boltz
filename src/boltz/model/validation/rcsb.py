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

# fmt: off


from typing import Optional

import torch
from pytorch_lightning import LightningModule

from boltz.model.validation.validator import Validator


class RCSBValidator(Validator):
    """Validation step implementation for RCSB."""

    def __init__(
        self,
        val_names: list[str],
        confidence_prediction: bool = False,
        physicalism_metrics: bool = False,
        override_val_method: Optional[str] = None,
    ) -> None:
        super().__init__(
            val_names=val_names,
            confidence_prediction=confidence_prediction,
            physicalism_metrics=physicalism_metrics,
            override_val_method=override_val_method,
        )

    def process(
        self,
        model: LightningModule,
        batch: dict[str, torch.Tensor],
        out: dict[str, torch.Tensor],
        idx_dataset: int,
    ) -> None:
        """Compute features.

        Parameters
        ----------
        model : LightningModule
            The LightningModule model.
        batch : Dict[str, torch.Tensor]
            The batch input.
        out : Dict[str, torch.Tensor]
            The output of the model.

        """
        symmetry_correction = model.val_group_mapper[idx_dataset]["symmetry_correction"]
        expand_to_diffusion_samples = (
            symmetry_correction  # True # TODO Mateo why is this set to sym correction?
        )

        # For now all was dumped into the common operation in the parent Validator class
        self.common_val_step(
            model,
            batch,
            out,
            idx_dataset,
            expand_to_diffusion_samples=expand_to_diffusion_samples,
        )

        # TODO: Implement the RCSB specific validation step

    def on_epoch_end(self, model: LightningModule) -> None:
        # For now all was dumped into the common operation in the parent Validator class
        self.common_on_epoch_end(model)
