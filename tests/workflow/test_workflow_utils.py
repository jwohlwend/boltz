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

"""Unit tests for workflow utilities."""

import pickle

import pytest
import torch
import torch.nn as nn
from omegaconf import DictConfig, ListConfig, OmegaConf
from pytorch_lightning import LightningModule, Trainer
from torch.utils.data import DataLoader, TensorDataset

from boltz.workflow.utils import CUDAMemoryProfile, convert_datasets_dict_to_list_config

# ========== Fixtures ==========


@pytest.fixture
def simple_base_config():
    """Create a simple base ListConfig with one item."""
    return OmegaConf.create([{"target_dir": "/path1", "prob": 0.5}])


@pytest.fixture
def multi_item_base_config():
    """Create a base ListConfig with multiple items."""
    return OmegaConf.create(
        [
            {"target_dir": "/path1", "prob": 0.5, "msa_dir": "/msa1"},
            {"target_dir": "/path2", "prob": 0.3, "msa_dir": "/msa2"},
            {"target_dir": "/path3", "prob": 0.2, "msa_dir": "/msa3"},
        ]
    )


@pytest.fixture
def keys_set_simple():
    """Simple set of valid keys."""
    return {"target_dir", "prob"}


@pytest.fixture
def keys_set_extended():
    """Extended set of valid keys."""
    return {"target_dir", "prob", "msa_dir"}


# ========== Happy Path Tests ==========


def test_single_item_single_key_override(simple_base_config, keys_set_simple):
    """Test overriding a single key in a single-item ListConfig."""
    override = OmegaConf.create({"0": {"target_dir": "/new/path"}})

    result = convert_datasets_dict_to_list_config(simple_base_config, override, keys_set_simple)

    assert isinstance(result, ListConfig)
    assert result[0].target_dir == override["0"]["target_dir"]
    assert result[0].prob == simple_base_config[0].prob  # Unchanged

    # Verify original is not modified
    assert simple_base_config[0].target_dir == "/path1"


def test_single_item_multiple_keys_override(simple_base_config, keys_set_simple):
    """Test overriding multiple keys in a single-item ListConfig."""
    override = OmegaConf.create({"0": {"target_dir": "/new/path", "prob": 0.9}})

    result = convert_datasets_dict_to_list_config(simple_base_config, override, keys_set_simple)

    assert result[0].target_dir == override["0"]["target_dir"]
    assert result[0].prob == override["0"]["prob"]


def test_multiple_items_single_override(multi_item_base_config, keys_set_extended):
    """Test overriding a single item in a multi-item ListConfig."""
    override = OmegaConf.create({"1": {"target_dir": "/new/path"}})

    result = convert_datasets_dict_to_list_config(multi_item_base_config, override, keys_set_extended)

    # First item unchanged
    assert result[0].target_dir == multi_item_base_config[0].target_dir
    assert result[0].prob == multi_item_base_config[0].prob
    assert result[0].msa_dir == multi_item_base_config[0].msa_dir

    # Second item modified
    assert result[1].target_dir == override["1"]["target_dir"]
    assert result[1].prob == multi_item_base_config[1].prob  # Unchanged
    assert result[1].msa_dir == multi_item_base_config[1].msa_dir  # Unchanged

    # Third item unchanged
    assert result[2].target_dir == multi_item_base_config[2].target_dir


def test_multiple_items_multiple_overrides(multi_item_base_config, keys_set_extended):
    """Test overriding multiple items in a multi-item ListConfig."""
    override = OmegaConf.create(
        {
            "0": {"target_dir": "/new/path0", "prob": 0.8},
            "2": {"msa_dir": "/new/msa2"},
        }
    )

    result = convert_datasets_dict_to_list_config(multi_item_base_config, override, keys_set_extended)

    # First item modified
    assert result[0].target_dir == override["0"]["target_dir"]
    assert result[0].prob == override["0"]["prob"]
    assert result[0].msa_dir == multi_item_base_config[0].msa_dir  # Unchanged

    # Second item unchanged
    assert result[1].target_dir == multi_item_base_config[1].target_dir
    assert result[1].prob == multi_item_base_config[1].prob

    # Third item modified
    assert result[2].target_dir == multi_item_base_config[2].target_dir  # Unchanged
    assert result[2].msa_dir == override["2"]["msa_dir"]


def test_all_items_overridden(multi_item_base_config, keys_set_extended):
    """Test overriding all items in a multi-item ListConfig."""
    override = OmegaConf.create(
        {
            "0": {"target_dir": "/new0"},
            "1": {"target_dir": "/new1"},
            "2": {"target_dir": "/new2"},
        }
    )

    result = convert_datasets_dict_to_list_config(multi_item_base_config, override, keys_set_extended)

    assert result[0].target_dir == override["0"]["target_dir"]
    assert result[1].target_dir == override["1"]["target_dir"]
    assert result[2].target_dir == override["2"]["target_dir"]


def test_deep_copy_behavior(simple_base_config, keys_set_simple):
    """Test that the function returns a deep copy and doesn't mutate the original."""
    original_target_dir = simple_base_config[0].target_dir
    override = OmegaConf.create({"0": {"target_dir": "/new/path"}})

    result = convert_datasets_dict_to_list_config(simple_base_config, override, keys_set_simple)

    # Modify the result
    modified_value = "/another/path"
    result[0].target_dir = modified_value

    # Original should be unchanged
    assert simple_base_config[0].target_dir == original_target_dir
    assert result[0].target_dir == modified_value


@pytest.mark.parametrize(
    "field_name,field_value,expected",
    [
        ("str_field", "new_value", "new_value"),
        ("int_field", 100, 100),
        ("float_field", 2.71, 2.71),
        ("bool_field", False, False),
    ],
)
def test_override_with_different_value_types(field_name, field_value, expected):
    """Test overriding with different value types (string, int, float, bool)."""
    base = OmegaConf.create(
        [
            {
                "str_field": "value",
                "int_field": 42,
                "float_field": 3.14,
                "bool_field": True,
            }
        ]
    )
    keys = {"str_field", "int_field", "float_field", "bool_field"}
    override = OmegaConf.create({"0": {field_name: field_value}})

    result = convert_datasets_dict_to_list_config(base, override, keys)

    assert result[0][field_name] == override["0"][field_name]
    assert result[0][field_name] == expected


def test_partial_key_override(multi_item_base_config, keys_set_extended):
    """Test overriding only some keys, leaving others unchanged."""
    override = OmegaConf.create({"1": {"prob": 0.99}})

    result = convert_datasets_dict_to_list_config(multi_item_base_config, override, keys_set_extended)

    assert result[1].target_dir == multi_item_base_config[1].target_dir  # Unchanged
    assert result[1].prob == override["1"]["prob"]  # Changed
    assert result[1].msa_dir == multi_item_base_config[1].msa_dir  # Unchanged


# ========== Error Cases: Invalid base ==========


@pytest.mark.parametrize(
    "base_input,error_match",
    [
        (OmegaConf.create({"0": {"target_dir": "/path"}}), "base must be a ListConfig"),  # DictConfig
        ([{"target_dir": "/path"}], "base must be a ListConfig"),  # Plain list
    ],
)
def test_base_not_listconfig_raises(base_input, error_match):
    """Test that passing non-ListConfig as base raises ValueError."""
    override = OmegaConf.create({"0": {"target_dir": "/new"}})
    keys = {"target_dir"}

    with pytest.raises(ValueError, match=error_match):
        convert_datasets_dict_to_list_config(base_input, override, keys)


# ========== Error Cases: Invalid override ==========


@pytest.mark.parametrize(
    "override_input,error_match",
    [
        (OmegaConf.create([{"target_dir": "/new"}]), "override must be a DictConfig"),  # ListConfig
        ({"0": {"target_dir": "/new"}}, "override must be a DictConfig"),  # Plain dict
    ],
)
def test_override_not_dictconfig_raises(simple_base_config, keys_set_simple, override_input, error_match):
    """Test that passing non-DictConfig as override raises ValueError."""
    with pytest.raises(ValueError, match=error_match):
        convert_datasets_dict_to_list_config(simple_base_config, override_input, keys_set_simple)


def test_empty_override_raises(simple_base_config, keys_set_simple):
    """Test that empty override DictConfig raises ValueError."""
    override = OmegaConf.create({})

    with pytest.raises(ValueError, match="Input DictConfig override is empty"):
        convert_datasets_dict_to_list_config(simple_base_config, override, keys_set_simple)


@pytest.mark.parametrize(
    "invalid_key,base_len,error_match",
    [
        ("3", 3, "Invalid keys in override"),  # Out of range
        ("-1", 1, "Invalid keys in override"),  # Negative
        ("abc", 1, "Invalid keys in override"),  # Non-integer
    ],
)
def test_override_with_invalid_index_raises(invalid_key, base_len, error_match):
    """Test that override with invalid index raises ValueError."""
    base = OmegaConf.create([{"field": f"value{i}"} for i in range(base_len)])
    override = OmegaConf.create({invalid_key: {"field": "/new"}})
    keys = {"field"}

    with pytest.raises(ValueError, match=error_match):
        convert_datasets_dict_to_list_config(base, override, keys)


def test_override_with_invalid_nested_keys_raises(simple_base_config, keys_set_simple):
    """Test that override with keys not in keys_to_override raises ValueError."""
    override = OmegaConf.create({"0": {"target_dir": "/new", "invalid_key": "value"}})

    with pytest.raises(ValueError, match="Invalid keys in override of item 0"):
        convert_datasets_dict_to_list_config(simple_base_config, override, keys_set_simple)


def test_override_with_only_invalid_keys_raises(simple_base_config, keys_set_simple):
    """Test that override with only invalid keys raises ValueError."""
    override = OmegaConf.create({"0": {"invalid_key": "value"}})

    with pytest.raises(ValueError, match="Invalid keys in override of item 0"):
        convert_datasets_dict_to_list_config(simple_base_config, override, keys_set_simple)


# ========== Null dataset removal ==========


def test_remove_null_dataset_single(multi_item_base_config, keys_set_extended):
    """Base has 2 items; override sets item 1 to null with remove_null_datasets=True; result has 1 item."""
    base = OmegaConf.create(
        [
            {"target_dir": "/path1", "prob": 0.5, "msa_dir": "/msa1"},
            {"target_dir": "/path2", "prob": 0.5, "msa_dir": "/msa2"},
        ]
    )
    override = OmegaConf.create({"1": None})

    result = convert_datasets_dict_to_list_config(base, override, keys_set_extended, remove_null_datasets=True)

    assert len(result) == 1
    assert result[0].target_dir == "/path1"
    assert result[0].prob == 0.5
    assert result[0].msa_dir == "/msa1"


def test_remove_null_dataset_multiple(multi_item_base_config, keys_set_extended):
    """Base has 3 items; override nullifies items 0 and 2; result has 1 item (former index 1)."""
    override = OmegaConf.create({"0": None, "2": None})

    result = convert_datasets_dict_to_list_config(
        multi_item_base_config, override, keys_set_extended, remove_null_datasets=True
    )

    assert len(result) == 1
    assert result[0].target_dir == "/path2"
    assert result[0].prob == 0.3
    assert result[0].msa_dir == "/msa2"


def test_null_override_without_flag_raises(multi_item_base_config, keys_set_extended):
    """Override sets an item to null with remove_null_datasets=False (default); raises ValueError."""
    override = OmegaConf.create({"1": None})

    with pytest.raises(ValueError, match="Override for item 1 is null but remove_null_datasets is False"):
        convert_datasets_dict_to_list_config(multi_item_base_config, override, keys_set_extended)


def test_remove_null_with_partial_override(multi_item_base_config, keys_set_extended):
    """Override nullifies item 1 and partially overrides item 0; result has one item with overrides applied."""
    override = OmegaConf.create({"0": {"target_dir": "/new/path0", "prob": 0.8}, "1": None})

    result = convert_datasets_dict_to_list_config(
        multi_item_base_config, override, keys_set_extended, remove_null_datasets=True
    )

    assert len(result) == 2  # item 0 (merged) and item 2 (unchanged)
    assert result[0].target_dir == "/new/path0"
    assert result[0].prob == 0.8
    assert result[0].msa_dir == "/msa1"
    assert result[1].target_dir == "/path3"
    assert result[1].prob == 0.2
    assert result[1].msa_dir == "/msa3"


# ========== Edge Cases ==========


def test_single_item_listconfig():
    """Test with a single-item ListConfig."""
    base = OmegaConf.create([{"field": "value"}])
    override = OmegaConf.create({"0": {"field": "new_value"}})
    keys = {"field"}

    result = convert_datasets_dict_to_list_config(base, override, keys)

    assert len(result) == len(base)
    assert result[0].field == override["0"]["field"]


def test_large_listconfig():
    """Test with a large ListConfig (10 items)."""
    base_size = 10
    base = OmegaConf.create([{"field": f"value{i}"} for i in range(base_size)])
    override = OmegaConf.create(
        {
            "0": {"field": "new0"},
            "5": {"field": "new5"},
            "9": {"field": "new9"},
        }
    )
    keys = {"field"}

    result = convert_datasets_dict_to_list_config(base, override, keys)

    assert len(result) == base_size
    assert result[0].field == override["0"]["field"]
    assert result[1].field == base[1].field  # Unchanged
    assert result[5].field == override["5"]["field"]
    assert result[9].field == override["9"]["field"]


@pytest.mark.parametrize(
    "override_value,expected",
    [
        (None, None),
        ("", ""),
        (0, 0),
    ],
)
def test_override_with_special_values(override_value, expected):
    """Test overriding with special values (None, empty string, zero)."""
    base = OmegaConf.create([{"field": "value"}])
    override = OmegaConf.create({"0": {"field": override_value}})
    keys = {"field"}

    result = convert_datasets_dict_to_list_config(base, override, keys)

    assert result[0].field == override["0"]["field"]
    assert result[0].field == expected


def test_empty_keys_set():
    """Test with empty keys_to_override set."""
    base = OmegaConf.create([{}])
    override = OmegaConf.create({"0": {}})
    keys = set()

    result = convert_datasets_dict_to_list_config(base, override, keys)

    # Should work since both base and override are empty
    assert len(result) == len(base)


def test_nested_dictconfig_values():
    """Test with nested DictConfig as values."""
    base = OmegaConf.create([{"config": {"nested": "value"}}])
    override = OmegaConf.create({"0": {"config": {"nested": "new_value"}}})
    keys = {"config"}

    result = convert_datasets_dict_to_list_config(base, override, keys)

    assert result[0].config.nested == override["0"]["config"]["nested"]


def test_list_values_in_config():
    """Test with list values in the config."""
    base = OmegaConf.create([{"data_list": [1, 2, 3]}])
    override = OmegaConf.create({"0": {"data_list": [4, 5, 6]}})
    keys = {"data_list"}

    result = convert_datasets_dict_to_list_config(base, override, keys)

    assert result[0].data_list == override["0"]["data_list"]


# ========== Relaxed validation: base may have extra keys ==========


def test_base_with_extra_keys_allowed():
    """Base items may have keys beyond keys_to_override; override works and extra keys preserved."""
    base = OmegaConf.create([{"target_dir": "/path1", "prob": 0.5, "extra_key": "extra_value", "another_extra": 42}])
    override = OmegaConf.create({"0": {"target_dir": "/new/path"}})
    keys = {"target_dir", "prob"}

    result = convert_datasets_dict_to_list_config(base, override, keys)

    assert result[0].target_dir == "/new/path"
    assert result[0].prob == 0.5
    assert result[0].extra_key == "extra_value"
    assert result[0].another_extra == 42


def test_base_items_with_different_keys():
    """Base items may have heterogeneous key sets; override of shared keys works."""
    # Mimics real YAML: dataset 0 has split/val_group, dataset 1 has override_method/override_bfactor
    base = OmegaConf.create(
        [
            {"target_dir": "/path1", "prob": 0.5, "split": "train", "val_group": "RCSB"},
            {"target_dir": "/path2", "prob": 0.5, "override_method": "AFDB", "override_bfactor": True},
        ]
    )
    keys = {"target_dir", "prob", "split", "val_group", "override_method", "override_bfactor"}
    override = OmegaConf.create({"0": {"target_dir": "/new0"}, "1": {"target_dir": "/new1"}})

    result = convert_datasets_dict_to_list_config(base, override, keys)

    assert result[0].target_dir == "/new0"
    assert result[0].split == "train"
    assert result[1].target_dir == "/new1"
    assert result[1].override_method == "AFDB"


def test_override_key_not_in_base_item_raises():
    """Override key in keys_to_override but absent from that base item raises ValueError."""
    base = OmegaConf.create(
        [
            {"target_dir": "/path1", "prob": 0.5, "split": "train"},
            {"target_dir": "/path2", "prob": 0.5},  # no "split" key
        ]
    )
    keys = {"target_dir", "prob", "split"}
    override = OmegaConf.create({"1": {"split": "val"}})  # override split for item 1 which has no split

    with pytest.raises(ValueError, match="contain keys not present in base item"):
        convert_datasets_dict_to_list_config(base, override, keys)


# ========== Documentation Example Test ==========


def test_docstring_example():
    """Test the example from the function's docstring."""
    base = OmegaConf.create([{"target_dir": "/path1", "prob": 0.5}])
    override = OmegaConf.create({"0": {"target_dir": "/new/path"}})
    keys = {"target_dir", "prob"}

    result = convert_datasets_dict_to_list_config(base, override, keys)

    assert result[0].target_dir == override["0"]["target_dir"]
    assert result[0].prob == base[0].prob


# ========== Real-world Scenario Tests ==========


def test_realistic_dataset_config():
    """Test with a realistic dataset configuration scenario."""
    base = OmegaConf.create(
        [
            {
                "_target_": "DatasetA",
                "target_dir": "/data/train",
                "msa_dir": "/msa/train",
                "prob": 0.7,
                "sampler": "uniform",
                "cropper": "random",
                "split": "train",
            },
            {
                "_target_": "DatasetB",
                "target_dir": "/data/val",
                "msa_dir": "/msa/val",
                "prob": 0.3,
                "sampler": "weighted",
                "cropper": "center",
                "split": "val",
            },
        ]
    )

    # Command-line override: Change first dataset's directory and probability
    override = OmegaConf.create(
        {
            "0": {
                "target_dir": "/new/train/path",
                "prob": 0.9,
            }
        }
    )

    keys = {"_target_", "target_dir", "msa_dir", "prob", "sampler", "cropper", "split"}

    result = convert_datasets_dict_to_list_config(base, override, keys)

    # First dataset modified
    assert result[0]._target_ == base[0]._target_  # Unchanged
    assert result[0].target_dir == override["0"]["target_dir"]  # Changed
    assert result[0].prob == override["0"]["prob"]  # Changed
    assert result[0].sampler == base[0].sampler  # Unchanged

    # Second dataset unchanged
    assert result[1].target_dir == base[1].target_dir
    assert result[1].prob == base[1].prob


def test_string_index_consistency(multi_item_base_config, keys_set_extended):
    """Test that string indices work correctly and consistently."""
    override = OmegaConf.create(
        {
            "0": {"target_dir": "/new0"},
            "1": {"target_dir": "/new1"},
        }
    )

    result = convert_datasets_dict_to_list_config(multi_item_base_config, override, keys_set_extended)

    assert result[0].target_dir == override["0"]["target_dir"]
    assert result[1].target_dir == override["1"]["target_dir"]
    assert result[2].target_dir == multi_item_base_config[2].target_dir  # Unchanged


def test_preserves_omegaconf_metadata(simple_base_config, keys_set_simple):
    """Test that OmegaConf metadata is preserved."""
    override = OmegaConf.create({"0": {"target_dir": "/new"}})

    result = convert_datasets_dict_to_list_config(simple_base_config, override, keys_set_simple)

    assert isinstance(result, ListConfig)
    assert isinstance(result[0], DictConfig)
    assert OmegaConf.is_config(result)
    assert OmegaConf.is_list(result)


# ========== CUDAMemoryProfile Tests ==========


class SimpleLightningModule(LightningModule):
    """Simple Lightning module for testing CUDAMemoryProfile."""

    def __init__(self):
        super().__init__()
        self.layer1 = nn.Linear(100, 50)
        self.layer2 = nn.Linear(50, 10)

    def forward(self, x):
        x = torch.relu(self.layer1(x))
        return self.layer2(x)

    def training_step(self, batch, batch_idx):
        x, y = batch
        y_pred = self(x)
        return nn.functional.mse_loss(y_pred, y)

    def configure_optimizers(self):
        return torch.optim.Adam(self.parameters(), lr=0.001)


@pytest.mark.skipif(not torch.cuda.is_available(), reason="CUDA not available")
def test_cuda_memory_profile_callback(tmp_path):
    """Test CUDAMemoryProfile callback creates memory snapshot file."""
    output_file = tmp_path / "memory_snapshot.pickle"
    memory_profiler = CUDAMemoryProfile(output_path=output_file, max_entries=100000)

    model = SimpleLightningModule()
    x_data = torch.randn(32, 100)
    y_data = torch.randn(32, 10)
    dataset = TensorDataset(x_data, y_data)
    dataloader = DataLoader(dataset, batch_size=8)

    trainer = Trainer(
        max_epochs=1,
        accelerator="gpu",
        devices=1,
        callbacks=[memory_profiler],
        enable_progress_bar=False,
        enable_model_summary=False,
        logger=False,
    )
    trainer.fit(model, dataloader)

    assert output_file.exists(), f"Memory snapshot file should be created at {output_file}"

    with open(output_file, "rb") as f:
        snapshot_data = pickle.load(f)

    assert isinstance(snapshot_data, dict), "Snapshot should be a dictionary"
    assert "segments" in snapshot_data or "device_traces" in snapshot_data, "Snapshot should contain memory data"

    if "device_traces" in snapshot_data and len(snapshot_data["device_traces"]) > 0:
        device_0_traces = snapshot_data["device_traces"][0]
        assert isinstance(device_0_traces, list), "Device traces should be a list"


@pytest.mark.skipif(not torch.cuda.is_available(), reason="CUDA not available")
def test_cuda_memory_profile_with_kwargs(tmp_path):
    """Test CUDAMemoryProfile callback forwards kwargs correctly."""
    output_file = tmp_path / "memory_snapshot_kwargs.pickle"
    max_entries = 50000
    memory_profiler = CUDAMemoryProfile(output_path=output_file, max_entries=max_entries)

    assert memory_profiler._kwargs == {"max_entries": max_entries}
    assert memory_profiler._output_path == output_file

    model = SimpleLightningModule()
    x_data = torch.randn(16, 100)
    y_data = torch.randn(16, 10)
    dataset = TensorDataset(x_data, y_data)
    dataloader = DataLoader(dataset, batch_size=8)

    trainer = Trainer(
        max_epochs=1,
        accelerator="gpu",
        devices=1,
        callbacks=[memory_profiler],
        enable_progress_bar=False,
        enable_model_summary=False,
        logger=False,
    )
    trainer.fit(model, dataloader)

    assert output_file.exists(), "Memory snapshot file should be created"


@pytest.mark.skipif(not torch.cuda.is_available(), reason="CUDA not available")
def test_cuda_memory_profile_creates_parent_dirs(tmp_path):
    """Test that CUDAMemoryProfile creates parent directories if they don't exist."""
    output_file = tmp_path / "nested" / "dirs" / "memory_snapshot.pickle"
    assert not output_file.parent.exists()

    memory_profiler = CUDAMemoryProfile(output_path=output_file)

    assert output_file.parent.exists(), "Parent directories should be created"
    assert output_file.parent.is_dir(), "Parent path should be a directory"

    model = SimpleLightningModule()
    x_data = torch.randn(8, 100)
    y_data = torch.randn(8, 10)
    dataset = TensorDataset(x_data, y_data)
    dataloader = DataLoader(dataset, batch_size=4)

    trainer = Trainer(
        max_epochs=1,
        accelerator="gpu",
        devices=1,
        callbacks=[memory_profiler],
        enable_progress_bar=False,
        enable_model_summary=False,
        logger=False,
    )
    trainer.fit(model, dataloader)

    assert output_file.exists(), "Memory snapshot should be created in nested directory"
