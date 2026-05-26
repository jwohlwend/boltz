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

# ruff: noqa
# fmt: off

import torch
import torch.nn as nn

import unittest

from boltz.model.layers.outer_product_mean import OuterProductMean


class OuterProductMeanTest(unittest.TestCase):
    def setUp(self):
        self.c_in = 32
        self.c_hidden = 16
        self.c_out = 64

        # Use torch.random.fork_rng + torch.no_grad context managers instead of
        # pytorch_lightning.seed_everything + torch.set_grad_enabled(False) to
        # avoid polluting global RNG state and leaking disabled gradients into
        # subsequent tests in the suite (the latter caused "element 0 of tensors
        # does not require grad" failures in test_triattn_kernel and
        # test_distogramv2).
        with torch.random.fork_rng(), torch.no_grad():
            torch.manual_seed(1100)
            self.layer = OuterProductMean(self.c_in, self.c_hidden, self.c_out)

            # Initialize layer
            for name, param in self.layer.named_parameters():
                nn.init.normal_(param, mean=1.0, std=1.0)

            # Set to eval mode
            self.layer.eval()

    def test_chunk(self):
        chunk_sizes = [16, 33, 64, 83, 100]
        B, S, N = 1, 49, 84
        m = torch.randn(size=(B, S, N, self.c_in))
        mask = torch.randint(low=0, high=1, size=(B, S, N))

        with torch.no_grad():
            exp_output = self.layer(m=m, mask=mask)
            for chunk_size in chunk_sizes:
                with self.subTest(chunk_size=chunk_size):
                    act_output = self.layer(m=m, mask=mask, chunk_size=chunk_size)
                    assert torch.allclose(exp_output, act_output, atol=1e-8)
