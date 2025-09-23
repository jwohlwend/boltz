import pytorch_lightning
import torch
import torch.nn as nn

import unittest

from boltz.model.layers.triangular_attention.attention import TriangleAttention


class OuterProductMeanTest(unittest.TestCase):
    def setUp(self):
        self.c_in = 128
        self.c_hidden = 32
        self.no_heads = 1

        torch.set_grad_enabled(False)
        pytorch_lightning.seed_everything(1100)
        self.layer = TriangleAttention(self.c_in, self.c_hidden, self.no_heads)

        # Initialize layer
        for name, param in self.layer.named_parameters():
            nn.init.normal_(param, mean=1.0, std=1.0)

    def test_chunk(self):
        chunk_sizes = [16, 33, 64, 100]
        B, N = 1, 99
        m = torch.randn(size=(B, N, N, self.c_in))
        mask = torch.randint(low=0, high=1, size=(B, N, N))

        with torch.no_grad():
            exp_output = self.layer(x=m, mask=mask)
            for chunk_size in chunk_sizes:
                with self.subTest(chunk_size=chunk_size):
                    act_output = self.layer(x=m, mask=mask, chunk_size=chunk_size)
                    assert torch.allclose(exp_output, act_output, atol=1e-8)
