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

import copy
import warnings
from pathlib import Path

import omegaconf
import pytorch_lightning as pl
import torch


def convert_datasets_dict_to_list_config(
    base: omegaconf.ListConfig,
    override: omegaconf.DictConfig,
    keys_to_override: set[str],
    *,
    remove_null_datasets: bool = False,
) -> omegaconf.ListConfig:
    """Convert a DictConfig with string integer keys to a ListConfig by merging into base.

    This function provides a workaround for OmegaConf's limitation of partially overriding
    entries nested in a ListConfig by enabling the user to specify the partial overrides
    using a DictConfig with string integer keys so that entry under those string integer keys
    are merged into the base ListConfig[index] structure.

    Base items may have extra keys beyond keys_to_override. Override entries may only
    contain keys that are both in keys_to_override and present in the corresponding
    base item.

    When remove_null_datasets is True, an override value of None (e.g. data.datasets.1=null
    on the CLI) removes that list entry from the result instead of merging.

    Args:
        base: The base ListConfig to merge overrides into. Items may have any keys;
            override keys must exist in the target base item.
        override: A DictConfig with string integer keys (e.g., {"0": {...}, "1": null}).
            Each key is an index in the base ListConfig. Values are either a DictConfig
            of keys to merge, or None to remove that entry (only when remove_null_datasets
            is True).
        keys_to_override: Set of keys that are allowed in override items (whitelist).
            Used for validation; override keys must be a subset of this set and of
            the target base item's keys.
        remove_null_datasets: If True, override entries with value None cause that
            list index to be removed from the result. If False, a None value raises
            ValueError.

    Returns:
        A new ListConfig that is a deep copy of base with the override values merged in at
        the specified indices, and with null-marked entries removed when remove_null_datasets
        is True. The user can use OmegaConf.merge to merge the returned ListConfig
        with the base ListConfig to get the effect of partial overrides as described above.

    Raises:
        ValueError: If base is not a ListConfig, override is not a DictConfig, override is
            empty, override contains invalid index keys, override items have keys not in
            keys_to_override, override items have keys not present in the corresponding
            base item, or an override value is None while remove_null_datasets is False.

    Example:
        >>> base = OmegaConf.create([{"target_dir": "/path1", "prob": 0.5}])
        >>> override = OmegaConf.create({"0": {"target_dir": "/new/path"}})
        >>> result = convert_datasets_dict_to_list_config(base, override, {"target_dir", "prob"})
        >>> result[0].target_dir
        '/new/path'
    """
    if not isinstance(base, omegaconf.ListConfig):
        raise ValueError(f"base must be a ListConfig, got {type(base)}")

    if not isinstance(override, omegaconf.DictConfig):
        raise ValueError(f"override must be a DictConfig, got {type(override)}")

    # expecting command line override data.datasets.<some_integer> where <some_integer> is an integer
    # in the range [0, len(base) - 1]
    keys_dataset_ids_override = set(map(str, range(len(base))))

    if len(override.keys()) == 0:
        raise ValueError(
            "Input DictConfig override is empty. "
            "Please specify at least one item using <some_integer> as its key "
            f"where <some_integer> is in the set of {keys_dataset_ids_override}"
        )

    if not (override.keys() <= keys_dataset_ids_override):
        raise ValueError(f"Invalid keys in override: {override.keys()}. Valid keys are: {keys_dataset_ids_override}")

    ans = copy.deepcopy(base)
    indices_to_remove = set()
    for i in range(len(base)):
        i_str = str(i)
        if i_str in override:
            if override[i_str] is None:
                if remove_null_datasets:
                    indices_to_remove.add(i)
                    continue
                raise ValueError(f"Override for item {i_str} is null but remove_null_datasets is False")
            if not (override[i_str].keys() <= keys_to_override):
                raise ValueError(
                    f"Invalid keys in override of item {i_str}: "
                    f"{override[i_str].keys()}. Valid keys are: {keys_to_override}"
                )
            if not (override[i_str].keys() <= base[i].keys()):
                raise ValueError(
                    f"Override keys {override[i_str].keys()} for item {i_str} "
                    f"contain keys not present in base item: "
                    f"{override[i_str].keys() - base[i].keys()}"
                )
            for k, v in override[i_str].items():
                ans[i][k] = v
    if indices_to_remove:
        ans = omegaconf.OmegaConf.create([ans[i] for i in range(len(ans)) if i not in indices_to_remove])
    return ans


# Default whitelist for CLI overrides of data.datasets[*], aligned with DatasetConfig (trainingv2).
# Used by train entrypoints when calling convert_datasets_dict_to_list_config.
_DATASET_KEYS_TO_OVERRIDE = {
    "_target_",
    "target_dir",
    "msa_dir",
    "prob",
    "sampler",
    "cropper",
    "template_dir",
    "filters",
    "split",
    "symmetry_correction",
    "val_group",
    "use_train_subset",
    "moldir",
    "override_bfactor",
    "override_method",
}


class CUDAMemoryProfile(pl.Callback):
    """PyTorch Lightning callback for profiling CUDA memory usage.

    Captures a detailed history of CUDA memory allocations and deallocations
    throughout training or prediction, then dumps a memory snapshot at the end.
    The snapshot can be analyzed with the PyTorch Memory Visualizer.

    Uses ``torch.cuda.memory._record_memory_history`` /
    ``torch.cuda.memory._dump_snapshot`` under the hood.

    Parameters
    ----------
    output_path : Path or str
        Path where the memory snapshot pickle file will be saved.
        Parent directories are created automatically.
    *args
        Forwarded to ``torch.cuda.memory._record_memory_history()``.
    **kwargs
        Forwarded to ``torch.cuda.memory._record_memory_history()``.
        Common kwargs include ``max_entries`` (default 100 000).

    Examples
    --------
    >>> profiler = CUDAMemoryProfile("profiling/mem_rank0.pickle", max_entries=300000)
    >>> trainer = pl.Trainer(callbacks=[profiler])
    """

    def __init__(self, output_path: Path | str, *args, **kwargs):
        super().__init__()
        self._output_path = Path(output_path) if isinstance(output_path, str) else output_path
        self._args = args
        self._kwargs = kwargs
        self._output_path.parent.mkdir(parents=True, exist_ok=True)

    def _start_recording(self) -> None:
        torch.cuda.memory._record_memory_history(*self._args, **self._kwargs)

    def _stop_and_dump(self) -> None:
        try:
            torch.cuda.memory._dump_snapshot(str(self._output_path))
        except Exception as e:
            warnings.warn(f"CUDAMemoryProfile: Failed to capture memory snapshot: {e}")
        torch.cuda.memory._record_memory_history(enabled=None)

    def on_predict_start(self, trainer: pl.Trainer, pl_module: pl.LightningModule) -> None:
        self._start_recording()

    def on_predict_end(self, trainer: pl.Trainer, pl_module: pl.LightningModule) -> None:
        self._stop_and_dump()

    def on_train_start(self, trainer: pl.Trainer, pl_module: pl.LightningModule) -> None:
        self._start_recording()

    def on_train_end(self, trainer: pl.Trainer, pl_module: pl.LightningModule) -> None:
        self._stop_and_dump()
