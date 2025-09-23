import os
import pickle
from dataclasses import asdict
import pprint

import torch
import torch.nn as nn

import pytest
import unittest

from lightning_fabric import seed_everything

from boltz.main import MODEL_URL
from boltz.model.model import Boltz1

import test_utils

tests_dir = os.path.dirname(os.path.abspath(__file__))
test_data_dir = os.path.join(tests_dir, "data")


@pytest.mark.regression
class RegressionTester(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        cache = os.path.expanduser("~/.boltz")
        checkpoint_url = MODEL_URL
        model_name = checkpoint_url.split("/")[-1]
        checkpoint = os.path.join(cache, model_name)
        if not os.path.exists(checkpoint):
            test_utils.download_file(checkpoint_url, checkpoint)

        regression_feats_path = os.path.join(
            test_data_dir, "ligand_regression_feats.pkl"
        )
        if not os.path.exists(regression_feats_path):
            regression_feats_url = "https://www.dropbox.com/scl/fi/1avbcvoor5jcnvpt07tp6/ligand_regression_feats.pkl?rlkey=iwtm9gpxgrbp51jbizq937pqf&st=jnbky253&dl=1"
            test_utils.download_file(regression_feats_url, regression_feats_path)

        regression_feats = torch.load(regression_feats_path, map_location=device)
        model_module: nn.Module = Boltz1.load_from_checkpoint(
            checkpoint, map_location=device
        )
        model_module.to(device)
        model_module.eval()

        coords = regression_feats["feats"]["coords"]
        # Coords should be rank 4
        if len(coords.shape) == 3:
            coords = coords.unsqueeze(0)
        regression_feats["feats"]["coords"] = coords
        for key, val in regression_feats["feats"].items():
            if hasattr(val, "to"):
                regression_feats["feats"][key] = val.to(device)

        cls.model_module = model_module.to(device)
        cls.regression_feats = regression_feats

    def test_input_embedder(self):
        exp_s_inputs = self.regression_feats["s_inputs"]
        act_s_inputs = self.model_module.input_embedder(self.regression_feats["feats"])

        assert torch.allclose(exp_s_inputs, act_s_inputs, atol=1e-5)

    def test_rel_pos(self):
        exp_rel_pos_encoding = self.regression_feats["relative_position_encoding"]
        act_rel_pos_encoding = self.model_module.rel_pos(self.regression_feats["feats"])

        assert torch.allclose(exp_rel_pos_encoding, act_rel_pos_encoding, atol=1e-5)

    @pytest.mark.slow
    def test_structure_output(self):
        exp_structure_output = self.regression_feats["structure_output"]
        s = self.regression_feats["s"]
        z = self.regression_feats["z"]
        s_inputs = self.regression_feats["s_inputs"]
        feats = self.regression_feats["feats"]
        relative_position_encoding = self.regression_feats["relative_position_encoding"]
        multiplicity_diffusion_train = self.regression_feats[
            "multiplicity_diffusion_train"
        ]

        self.model_module.structure_module.coordinate_augmentation = False
        self.model_module.structure_module.sigma_data = 0.0

        seed_everything(self.regression_feats["seed"])
        act_structure_output = self.model_module.structure_module(
            s_trunk=s,
            z_trunk=z,
            s_inputs=s_inputs,
            feats=feats,
            relative_position_encoding=relative_position_encoding,
            multiplicity=multiplicity_diffusion_train,
        )

        act_keys = act_structure_output.keys()
        exp_keys = exp_structure_output.keys()
        assert act_keys == exp_keys

        # Other keys have some randomness, so we will only check the keys that
        # we can make deterministic with sigma_data = 0.0 (above).
        check_keys = ["noised_atom_coords", "aligned_true_atom_coords"]
        for key in check_keys:
            exp_val = exp_structure_output[key]
            act_val = act_structure_output[key]
            assert exp_val.shape == act_val.shape, f"Shape mismatch in {key}"
            assert torch.allclose(exp_val, act_val, atol=1e-4)


if __name__ == "__main__":
    unittest.main()
