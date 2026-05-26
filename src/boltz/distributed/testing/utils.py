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


from pathlib import Path

import hydra
import torch
from omegaconf import OmegaConf
from torch import Tensor
from torch.distributed.device_mesh import DeviceMesh
from torch.distributed.tensor import DTensor, Replicate, Shard, distribute_tensor

from boltz.data.module.trainingv2 import DataConfigV2

# Config file paths
ROOT_DIR = Path(__file__).resolve().parents[4]
CONFIG_FILE_BASE = ROOT_DIR / "scripts" / "train" / "configs" / "structurev2.yaml"


def create_atom_to_token_dtensor(atom_to_token_global: Tensor, device_mesh: DeviceMesh) -> DTensor:
    """Create a distributed tensor for atom_to_token with proper placement.

    Args:
        atom_to_token_global: Global atom_to_token tensor of shape (B, n_atoms, n_tokens)
        device_mesh: DeviceMesh instance

    Returns:
        DTensor: Distributed atom_to_token tensor with placement (Shard(0), Shard(1), Replicate())
    """
    # Get block diagonal chunk of atom_to_token_global
    n_atoms, n_tokens = atom_to_token_global.shape[1:]
    cp_axis_0_size = device_mesh.get_group("cp_axis_0").size()

    atom_to_token_local = []
    for cp_idx in range(cp_axis_0_size):
        start_token_idx = cp_idx * n_tokens // cp_axis_0_size
        end_token_idx = (cp_idx + 1) * n_tokens // cp_axis_0_size
        start_atom_idx = cp_idx * n_atoms // cp_axis_0_size
        end_atom_idx = (cp_idx + 1) * n_atoms // cp_axis_0_size
        atom_to_token_local.append(atom_to_token_global[:, start_atom_idx:end_atom_idx, start_token_idx:end_token_idx])

    atom_to_token_local = torch.cat(atom_to_token_local, dim=1)

    placements = (Shard(dim=0), Shard(dim=1), Replicate())
    atom_to_token_dtensor = distribute_tensor(atom_to_token_local, device_mesh, placements)
    return atom_to_token_dtensor


def setup_mock_training_datamodule_config(test_data_dir: Path) -> DataConfigV2:
    """Setup mock training datamodule configuration by loading and merging config files.

    Args:
        test_data_dir: Base path for test data directory

    Returns:
        Configured DataConfigV2 instance with test data paths
    """
    config_dict = OmegaConf.load(CONFIG_FILE_BASE)

    # Provide temporary valid paths for required string fields before instantiate.
    for dataset_cfg in config_dict.data.datasets:
        dataset_cfg.target_dir = "."
        dataset_cfg.msa_dir = "."

    # Instantiate the configuration
    cfg = hydra.utils.instantiate(config_dict)

    data_config = DataConfigV2(**cfg.data)

    # Test data comes from RCSB only; keep a single dataset entry.
    if len(data_config.datasets) > 1:
        data_config.datasets = [data_config.datasets[0]]

    # Override paths to use the prepared Boltz2 training test dataset layout.
    data_config.datasets[0].target_dir = str(test_data_dir)
    data_config.datasets[0].msa_dir = str(test_data_dir / "msa")
    data_config.datasets[0].split = None
    data_config.datasets[0].template_dir = None
    data_config.datasets[0].prob = 1.0

    # Keep tests small and deterministic.
    data_config.samples_per_epoch = 4
    data_config.num_workers = 0
    data_config.pin_memory = False
    data_config.use_templates = False

    # Enable symmetry features for training to test symmetry feature broadcasting
    data_config.return_train_symmetries = True

    return data_config
