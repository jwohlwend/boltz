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


def distogram_loss(
    output: dict[str, Tensor],
    feats: dict[str, Tensor],
    aggregate_distogram: bool = True,
) -> tuple[Tensor, Tensor]:
    """Compute the  distogram loss.

    Parameters
    ----------
    output : Dict[str, Tensor]
        Output of the model
    feats : Dict[str, Tensor]
        Input features

    Returns
    -------
    Tensor
        The globally averaged loss.
    Tensor
        Per example loss.

    """
    with torch.autocast("cuda", enabled=False):
        # Get predicted distograms
        compute_dtype = torch.promote_types(output["pdistogram"].dtype, torch.float32)
        pred = output["pdistogram"].to(compute_dtype)  # (B, L, L, num_distograms, disto_bins)
        D = pred.shape[3]  # num_distograms  # noqa: N806
        assert len(pred.shape) == 5  # noqa: PLR2004

        # Compute target distogram
        target = feats["disto_target"]  # (B, L, L, K, disto_bins)
        assert len(target.shape) == 5  # noqa: PLR2004

        if aggregate_distogram:
            msg = "Cannot aggregate GT distogram when num_distograms > 1"
            assert pred.shape[3] == 1, msg

            pred = pred.squeeze(3)  # (B, L, L, disto_bins)

            # Aggregate distogram over K conformers
            target = target.sum(dim=3)  # (B, L, L, disto_bins)

            # Normalize distogram
            P = target / target.sum(-1)[..., None].clamp(min=1)  # noqa: N806

            # Combine target mask and padding mask
            mask = feats["token_disto_mask"]
            mask = mask[:, None, :] * mask[:, :, None]
            mask = mask * (1 - torch.eye(mask.shape[1])[None]).to(pred)

            # Compute the distogram loss
            log_Q = torch.nn.functional.log_softmax(pred, dim=-1)  # noqa: N806
            errors = -1 * torch.sum(
                P * log_Q,
                dim=-1,
            )
            denom = 1e-5 + torch.sum(mask, dim=(-1, -2))
            mean = errors * mask
            mean = torch.sum(mean, dim=-1)
            mean = mean / denom[..., None]
            batch_loss = torch.sum(mean, dim=-1)
            global_loss = torch.mean(batch_loss)
        else:
            # We want to compute the loss for each pair of conformer K and predicted
            # distogram

            # Loop through conformers and compute the loss
            batch_loss = []
            for k in range(target.shape[3]):
                # Get the target distogram for conformer k
                # (B, L, L, K, disto_bins) -> (B, L, L, D, disto_bins)
                P_k = target[:, :, :, k : k + 1, :].repeat_interleave(D, dim=3)  # noqa: N806

                # Compute the distogram loss to all predicted distograms
                log_Q = torch.nn.functional.log_softmax(pred, dim=-1)  # noqa: N806
                errors = -1 * torch.sum(
                    P_k * log_Q,
                    dim=-1,
                )  # (B, L, L, D)

                # Compute mask
                mask = feats["token_disto_mask"]
                mask = mask[:, None, :] * mask[:, :, None]
                mask = mask * (1 - torch.eye(mask.shape[1])[None]).to(pred)
                mask = mask.unsqueeze(-1).repeat_interleave(D, -1)  # (B, L, L, D)

                denom = 1e-5 + torch.sum(mask, dim=(-2, -3))  # (B, D)
                mean = errors * mask
                mean = torch.sum(mean, dim=-2)  # (B, L, D)
                mean = mean / denom[..., None, :]
                b_loss = torch.sum(mean, dim=-2)  # (B, D)

                batch_loss.append(b_loss)

            batch_loss = torch.stack(batch_loss, dim=1)  # (B, K, D)

            # Compute the batch loss by taking the min over the predicted distograms
            # and the average across conformers
            batch_loss = torch.min(batch_loss, dim=-1).values.mean(dim=1)
            global_loss = torch.mean(batch_loss)

        return global_loss, batch_loss
