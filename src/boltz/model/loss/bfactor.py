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

import torch
from torch import Tensor


def bfactor_loss_fn(
    output: dict[str, Tensor],
    feats: dict[str, Tensor],
) -> Tensor:
    """Compute the bfactor loss.

    Parameters
    ----------
    output : dict[str, Tensor]
        Output of the model
    feats : dict[str, Tensor]
        Input features

    Returns
    -------
    Tensor
        The globally averaged loss.

    """
    with torch.autocast("cuda", enabled=False):
        compute_dtype = torch.promote_types(output["pbfactor"].dtype, torch.float32)
        pred = output["pbfactor"].to(compute_dtype)  # (B, L, bins)
        bins = pred.shape[2]  # num_bins
        token_to_rep_atom = feats["token_to_rep_atom"]

        # Compute target histogram
        bfactor_atom = feats["bfactor"].unsqueeze(-1)  # (B, L)
        bfactor_token = torch.bmm(token_to_rep_atom.to(compute_dtype), bfactor_atom)

        boundaries = torch.linspace(0, 100, bins - 1, device=bfactor_token.device)
        bfactor_token_bin = (bfactor_token > boundaries).sum(dim=-1).long()
        bfactor_target = torch.nn.functional.one_hot(
            bfactor_token_bin, num_classes=bins
        )

        # Combine target mask and padding mask
        token_mask = (bfactor_token > 1e-5).squeeze(-1).to(compute_dtype)

        # Compute the bfactor loss
        errors = -1 * torch.sum(
            bfactor_target * torch.nn.functional.log_softmax(pred, dim=-1),
            dim=-1,
        )
        loss = torch.sum(errors * token_mask) / (torch.sum(token_mask) + 1e-5)
        return loss
