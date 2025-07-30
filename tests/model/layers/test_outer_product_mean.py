import pytorch_lightning
import torch
import torch.nn as nn

import unittest

from boltz.model.layers.outer_product_mean import OuterProductMean


class OuterProductMeanTest(unittest.TestCase):
    def setUp(self):
        self.c_in = 32
        self.c_hidden = 16
        self.c_out = 64

        torch.set_grad_enabled(False)
        pytorch_lightning.seed_everything(1100)
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
