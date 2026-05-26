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


"""Distributed RCSB validator for Boltz-2 with context parallelism.

Uses diamond inheritance so that:

    DistributedRCSBValidator
        ├── DistributedValidator   (DTensor gathering + metric overrides)
        └── RCSBValidator          (process / on_epoch_end flow)
                └── Validator      (metric storage, serial compute logic)

MRO: DistributedRCSBValidator → DistributedValidator → RCSBValidator → Validator

This gives ``DistributedValidator`` precedence for every method it
overrides (``run_model``, ``common_val_step``, ``compute_disto_loss``, …)
while ``RCSBValidator.process()`` and ``RCSBValidator.on_epoch_end()``
provide the RCSB-specific entry points that delegate to
``self.common_val_step`` / ``self.common_on_epoch_end`` – both of which
resolve to the distributed versions from ``DistributedValidator``.
"""

from typing import Optional

import torch
from pytorch_lightning import LightningModule

from boltz.distributed.comm import TransposeComm
from boltz.distributed.model.validation.validator import DistributedValidator
from boltz.model.validation.rcsb import RCSBValidator


class DistributedRCSBValidator(DistributedValidator, RCSBValidator):
    """Distributed RCSB validator with DTensor-aware metric computation.

    Inherits:
    * Metric overrides & gathering from :class:`DistributedValidator`.
    * ``process`` / ``on_epoch_end`` flow from :class:`RCSBValidator`.
    """

    def __init__(
        self,
        val_names: list[str],
        confidence_prediction: bool = False,
        physicalism_metrics: bool = False,
        rmsd_metrics: bool = False,
        clash_score_metrics: bool = False,
        override_val_method: Optional[str] = None,
    ) -> None:
        """Initialize the distributed RCSB validator.

        Parameters
        ----------
        val_names : list[str]
            The list of validation names.
        confidence_prediction : bool
            Whether to predict confidence.
        physicalism_metrics : bool
            Whether to compute physicalism metrics.
        rmsd_metrics : bool
            Whether to compute rmsd metrics.
        clash_score_metrics : bool
            Whether to compute clash score metrics.
        override_val_method : Optional[str]
            The override validation method.
        """
        # Bypass cooperative super().__init__ to avoid MRO signature
        # mismatches between DistributedValidator and RCSBValidator.
        # Both ultimately initialise the same Validator base; calling it
        # directly avoids passing kwargs that the intermediate classes
        # don't expect.
        DistributedValidator.__init__(
            self,
            val_names=val_names,
            confidence_prediction=confidence_prediction,
            physicalism_metrics=physicalism_metrics,
            rmsd_metrics=rmsd_metrics,
            clash_score_metrics=clash_score_metrics,
            override_val_method=override_val_method,
        )

    def process(
        self,
        model: LightningModule,
        batch: dict[str, torch.Tensor],
        out: dict[str, torch.Tensor],
        idx_dataset: int,
        transpose_comm: TransposeComm,
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
        idx_dataset : int
            Global dataset index.
        transpose_comm : TransposeComm
            The transpose communication object.
        """
        symmetry_correction = model.val_group_mapper[idx_dataset]["symmetry_correction"]
        expand_to_diffusion_samples = symmetry_correction  # True # TODO Mateo why is this set to sym correction?

        # For now all was dumped into the common operation in the parent Validator class
        self.common_val_step(
            model,
            batch,
            out,
            idx_dataset,
            expand_to_diffusion_samples=expand_to_diffusion_samples,
            transpose_comm=transpose_comm,
        )
