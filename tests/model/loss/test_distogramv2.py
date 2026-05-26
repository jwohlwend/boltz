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

"""Tests for the Boltz-2 distogramv2 serial loss function.

Tests verify:
1. v2 matches v1 for K=1, D=1 (backward-compatible legacy behavior)
2. Aggregate invariant: K identical conformers == K=1
3. Aggregate error: D>1 raises
4. Non-aggregate semantic: min-over-D selects best prediction
5. Masking: diagonal + padded tokens get zero gradient; fully masked batch → ~0 loss
6. Gradient correctness: finite across all (K,D,agg) configs + numerical gradcheck
7. Autocast: loss stays float32
"""

import pytest
import torch

from boltz.model.loss.distogram import distogram_loss as distogram_loss_v1
from boltz.model.loss.distogramv2 import distogram_loss as distogram_loss_v2


def test_aggregate_k1_matches_v1():
    """v2 aggregate with K=1, D=1 must produce identical results to v1 (forward + backward)."""
    B, N, num_bins = 2, 16, 64
    D = 1

    with torch.random.fork_rng():
        torch.manual_seed(42)

        pred_v2 = torch.randn(B, N, N, D, num_bins, requires_grad=True)
        pred_v1 = pred_v2.squeeze(3).detach().clone().requires_grad_(True)

        target_idx = torch.randint(0, num_bins, (B, N, N))
        target_onehot = torch.nn.functional.one_hot(target_idx, num_classes=num_bins).float()
        target_v2 = target_onehot.unsqueeze(3)
        target_v1 = target_onehot

        mask = torch.ones(B, N)
        mask[0, 12:] = 0
        mask[1, 8:] = 0

        global_loss_v2, batch_loss_v2 = distogram_loss_v2(
            {"pdistogram": pred_v2}, {"disto_target": target_v2, "token_disto_mask": mask}, aggregate_distogram=True
        )
        global_loss_v1, batch_loss_v1 = distogram_loss_v1(
            {"pdistogram": pred_v1}, {"disto_target": target_v1, "token_disto_mask": mask}
        )

        torch.testing.assert_close(global_loss_v2, global_loss_v1)
        torch.testing.assert_close(batch_loss_v2, batch_loss_v1)

        global_loss_v2.backward()
        global_loss_v1.backward()
        torch.testing.assert_close(pred_v2.grad.squeeze(3), pred_v1.grad)


def test_aggregate_identical_conformers():
    """K identical conformers must produce the same loss and gradient as K=1."""
    B, N, num_bins = 2, 8, 16
    K, D = 3, 1

    with torch.random.fork_rng():
        torch.manual_seed(42)

        pred = torch.randn(B, N, N, D, num_bins, requires_grad=True)
        pred_single = pred.detach().clone().requires_grad_(True)

        target_idx = torch.randint(0, num_bins, (B, N, N))
        target_onehot = torch.nn.functional.one_hot(target_idx, num_classes=num_bins).float()
        target_multi = target_onehot.unsqueeze(3).repeat(1, 1, 1, K, 1)
        target_single = target_onehot.unsqueeze(3)

        mask = torch.ones(B, N)

        loss_multi, _ = distogram_loss_v2(
            {"pdistogram": pred}, {"disto_target": target_multi, "token_disto_mask": mask}, aggregate_distogram=True
        )
        loss_single, _ = distogram_loss_v2(
            {"pdistogram": pred_single},
            {"disto_target": target_single, "token_disto_mask": mask},
            aggregate_distogram=True,
        )

        torch.testing.assert_close(loss_multi, loss_single)

        loss_multi.backward()
        loss_single.backward()
        torch.testing.assert_close(pred.grad, pred_single.grad)


def test_aggregate_rejects_multi_distogram():
    """Aggregate mode must reject D>1."""
    pred = torch.randn(2, 8, 8, 2, 16)
    target = torch.rand(2, 8, 8, 1, 16)
    mask = torch.ones(2, 8)

    with pytest.raises(AssertionError, match="Cannot aggregate GT distogram when num_distograms > 1"):
        distogram_loss_v2(
            {"pdistogram": pred}, {"disto_target": target, "token_disto_mask": mask}, aggregate_distogram=True
        )


def test_non_aggregate_min_selects_best_prediction():
    """Non-aggregate min-over-D should yield low loss when one prediction matches the target."""
    B, N, num_bins = 1, 4, 8
    D = 2

    with torch.random.fork_rng():
        torch.manual_seed(42)

        target_idx = torch.randint(0, num_bins, (B, N, N))
        target = torch.nn.functional.one_hot(target_idx, num_classes=num_bins).float().unsqueeze(3)

        pred = torch.randn(B, N, N, D, num_bins)
        # First prediction gets high logit for correct bin; second is random
        pred_good = torch.zeros(B, N, N, num_bins)
        pred_good.scatter_(-1, target_idx.unsqueeze(-1), 10.0)
        with torch.no_grad():
            pred[:, :, :, 0, :] = pred_good
            pred[:, :, :, 1, :] = torch.randn(B, N, N, num_bins)
        pred.requires_grad_(True)

        mask = torch.ones(B, N)
        global_loss, _ = distogram_loss_v2(
            {"pdistogram": pred}, {"disto_target": target, "token_disto_mask": mask}, aggregate_distogram=False
        )

        assert global_loss < 1.0, f"Loss should be low when one prediction matches, got {global_loss}"


def test_masking_zeroes_gradients():
    """Diagonal, padded-token, and fully-masked-batch positions must get zero gradient / ~0 loss."""
    B, N, num_bins = 2, 8, 16
    K, D = 1, 1

    with torch.random.fork_rng():
        torch.manual_seed(42)

        # --- Part 1: diagonal + padded tokens ---
        pred = torch.randn(B, N, N, D, num_bins, requires_grad=True)
        target = torch.rand(B, N, N, K, num_bins)

        mask = torch.ones(B, N)
        mask[:, N // 2 :] = 0  # mask out second half of tokens

        global_loss, _ = distogram_loss_v2(
            {"pdistogram": pred}, {"disto_target": target, "token_disto_mask": mask}, aggregate_distogram=True
        )
        global_loss.backward()

        # Diagonal elements should have zero gradient
        for b in range(B):
            for i in range(N):
                assert torch.allclose(
                    pred.grad[b, i, i, :, :], torch.zeros_like(pred.grad[b, i, i, :, :]), atol=1e-6
                ), f"Diagonal gradient at [{b},{i},{i}] should be zero"

        # Masked rows and columns should have zero gradient
        for b in range(B):
            for i in range(N // 2, N):
                assert torch.allclose(
                    pred.grad[b, i, :, :, :], torch.zeros_like(pred.grad[b, i, :, :, :]), atol=1e-6
                ), f"Gradient for masked row {i} should be zero"
                assert torch.allclose(
                    pred.grad[b, :, i, :, :], torch.zeros_like(pred.grad[b, :, i, :, :]), atol=1e-6
                ), f"Gradient for masked column {i} should be zero"

        # --- Part 2: fully masked batch element ---
        pred2 = torch.randn(B, N, N, D, num_bins, requires_grad=True)
        target2 = torch.rand(B, N, N, K, num_bins)
        mask2 = torch.ones(B, N)
        mask2[1, :] = 0

        _, batch_loss = distogram_loss_v2(
            {"pdistogram": pred2}, {"disto_target": target2, "token_disto_mask": mask2}, aggregate_distogram=True
        )

        assert batch_loss[1] < 1e-4, f"Fully masked batch should have ~0 loss, got {batch_loss[1]}"


@pytest.mark.parametrize(
    "loss_config",
    [
        (1, 1, True),
        (3, 1, True),
        (1, 2, False),
        (3, 2, False),
    ],
    ids=lambda c: f"K={c[0]}|D={c[1]}|agg={c[2]}",
)
def test_gradient_finite(loss_config):
    """Gradients must be finite for all (K, D, aggregate) configurations."""
    K, D, aggregate = loss_config
    B, N, num_bins = 2, 16, 64

    with torch.random.fork_rng():
        torch.manual_seed(42)

        pred = torch.randn(B, N, N, D, num_bins, requires_grad=True)
        target = torch.rand(B, N, N, K, num_bins)
        target = target / target.sum(dim=-1, keepdim=True).clamp(min=1e-8)
        mask = torch.ones(B, N)
        mask[0, N // 2 :] = 0

        global_loss, batch_loss = distogram_loss_v2(
            {"pdistogram": pred}, {"disto_target": target, "token_disto_mask": mask}, aggregate_distogram=aggregate
        )

        # Shape sanity (replaces removed test_output_shapes)
        assert global_loss.shape == (), f"Global loss should be scalar, got {global_loss.shape}"
        assert batch_loss.shape == (B,), f"Batch loss should be [B], got {batch_loss.shape}"

        global_loss.backward()
        assert not torch.isnan(pred.grad).any(), f"NaN in gradients for {loss_config}"
        assert not torch.isinf(pred.grad).any(), f"Inf in gradients for {loss_config}"


def test_gradient_numerical_check():
    """Numerical gradient check (torch.autograd.gradcheck) for aggregate mode."""
    B, N, num_bins = 1, 4, 8

    with torch.random.fork_rng():
        torch.manual_seed(42)

        pred = torch.randn(B, N, N, 1, num_bins, requires_grad=True, dtype=torch.float64)
        target = torch.rand(B, N, N, 1, num_bins, dtype=torch.float64)
        target = target / target.sum(dim=-1, keepdim=True).clamp(min=1e-8)
        mask = torch.ones(B, N, dtype=torch.float64)

        def loss_fn(p):
            return distogram_loss_v2(
                {"pdistogram": p}, {"disto_target": target, "token_disto_mask": mask}, aggregate_distogram=True
            )[0]

        result = torch.autograd.gradcheck(loss_fn, pred, eps=1e-5, atol=1e-3, rtol=1e-2, nondet_tol=1e-3)
        assert result, "Gradient check failed for aggregate mode"


@pytest.mark.skipif(not torch.cuda.is_available(), reason="CUDA required for autocast test")
def test_autocast_disabled_in_loss():
    """Loss computation must stay float32 even under autocast."""
    B, N, num_bins = 2, 8, 16

    with torch.random.fork_rng(devices=["cuda"]):
        torch.manual_seed(42)

        pred = torch.randn(B, N, N, 1, num_bins, device="cuda", requires_grad=True)
        target = torch.rand(B, N, N, 1, num_bins, device="cuda")
        mask = torch.ones(B, N, device="cuda")

        with torch.autocast("cuda", dtype=torch.float16):
            global_loss, batch_loss = distogram_loss_v2(
                {"pdistogram": pred}, {"disto_target": target, "token_disto_mask": mask}, aggregate_distogram=True
            )

        assert global_loss.dtype == torch.float32, f"Global loss should be float32, got {global_loss.dtype}"
        assert batch_loss.dtype == torch.float32, f"Batch loss should be float32, got {batch_loss.dtype}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
