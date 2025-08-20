import unittest

import torch
import torch.nn as nn
import pytorch_lightning

from boltz.model.modules.affinity import PolymerPropertyHeadsTransformer


class PolymerPropertyHeadsTransformerTest(unittest.TestCase):
    def setUp(self):
        pytorch_lightning.seed_everything(1100)
        torch.set_grad_enabled(False)

        self.token_z = 32
        self.input_token_s = 16

        self.head = PolymerPropertyHeadsTransformer(
            token_z=self.token_z,
            input_token_s=self.input_token_s,
            num_blocks=1,
            num_heads=1,
            activation_checkpointing=False,
            use_cross_transformer=False,
            groups={},
        )

        # Initialize deterministically
        for _, param in self.head.named_parameters():
            if param.requires_grad:
                nn.init.normal_(param, mean=0.0, std=0.02)

        self.head.eval()

    def test_forward_shape(self):
        B, L = 3, 10
        z = torch.randn(B, L, L, self.token_z)
        # Ensure at least one valid token to avoid degenerate all-zero mask
        pad_mask = torch.ones(B, L, dtype=torch.bool)
        feats = {"token_pad_mask": pad_mask}

        with torch.no_grad():
            out = self.head(z=z, feats=feats, multiplicity=1)

        self.assertIn("polymer_properties_pred", out)
        y = out["polymer_properties_pred"]
        self.assertEqual(y.shape, (B, 5))
        self.assertTrue(torch.isfinite(y).all())

    def test_forward_with_multiplicity(self):
        B, L, mult = 2, 7, 3
        z = torch.randn(B * mult, L, L, self.token_z)
        # Random padding mask (some zeros) for base batch dim B
        pad_mask = torch.randint(low=0, high=2, size=(B, L)).to(torch.bool)
        feats = {"token_pad_mask": pad_mask}

        with torch.no_grad():
            out = self.head(z=z, feats=feats, multiplicity=mult)

        y = out["polymer_properties_pred"]
        self.assertEqual(y.shape, (B * mult, 5))
        self.assertTrue(torch.isfinite(y).all())


if __name__ == "__main__":
    unittest.main()


