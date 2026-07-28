import tempfile
import unittest
from pathlib import Path

import numpy as np
import torch

from boltz.data.write.writer import BoltzWriter


class WriteDistogramTest(unittest.TestCase):
    def setUp(self):
        self.num_tokens = 7
        self.num_padded = 10
        self.num_bins = 64

        self.writer = BoltzWriter(
            data_dir=tempfile.mkdtemp(),
            output_dir=tempfile.mkdtemp(),
            write_distogram=True,
        )

        pad = self.num_padded - self.num_tokens
        self.batch = {
            "token_pad_mask": torch.tensor([[1.0] * self.num_tokens + [0.0] * pad]),
            "asym_id": torch.tensor([[0] * 4 + [1] * 3 + [0] * pad]),
            "residue_index": torch.arange(self.num_padded)[None],
            "entity_id": torch.zeros(1, self.num_padded, dtype=torch.long),
            "mol_type": torch.zeros(1, self.num_padded, dtype=torch.long),
        }

    def _write(self, logits):
        with tempfile.TemporaryDirectory() as tmp_dir:
            path = Path(tmp_dir) / "distogram.npz"
            self.writer._write_distogram(  # noqa: SLF001
                path, {"pdistogram": logits}, self.batch, 0
            )
            return dict(np.load(path))

    def test_boltz1_shape(self):
        """A single distogram is written unchanged."""
        logits = torch.randn(1, self.num_padded, self.num_padded, self.num_bins)
        arrays = self._write(logits)

        expected = logits[0, : self.num_tokens, : self.num_tokens].numpy()
        self.assertEqual(
            arrays["distogram_logits"].shape,
            (self.num_tokens, self.num_tokens, self.num_bins),
        )
        np.testing.assert_array_equal(
            arrays["distogram_logits"], expected.astype(np.float16)
        )

    def test_boltz2_shape(self):
        """A stack holding a single distogram loses its stacking axis."""
        logits = torch.randn(1, self.num_padded, self.num_padded, 1, self.num_bins)
        arrays = self._write(logits)

        expected = logits[0, : self.num_tokens, : self.num_tokens, 0].numpy()
        self.assertEqual(
            arrays["distogram_logits"].shape,
            (self.num_tokens, self.num_tokens, self.num_bins),
        )
        np.testing.assert_array_equal(
            arrays["distogram_logits"], expected.astype(np.float16)
        )

    def test_distogram_stack_is_kept(self):
        """A stack holding several distograms is written whole."""
        logits = torch.randn(1, self.num_padded, self.num_padded, 3, self.num_bins)
        arrays = self._write(logits)

        self.assertEqual(
            arrays["distogram_logits"].shape,
            (self.num_tokens, self.num_tokens, 3, self.num_bins),
        )

    def test_token_annotation(self):
        """The token annotation covers the unpadded tokens only."""
        logits = torch.randn(1, self.num_padded, self.num_padded, 1, self.num_bins)
        arrays = self._write(logits)

        for key in ("asym_id", "residue_index", "entity_id", "mol_type"):
            self.assertEqual(arrays[key].shape, (self.num_tokens,))
        np.testing.assert_array_equal(arrays["asym_id"], [0, 0, 0, 0, 1, 1, 1])
        np.testing.assert_array_equal(arrays["residue_index"], np.arange(7))


if __name__ == "__main__":
    unittest.main()
