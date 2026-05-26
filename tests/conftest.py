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

import os
import pickle
import re
import shutil
import tarfile
import tempfile
import urllib.request
import warnings
from collections import OrderedDict
from copy import deepcopy
from dataclasses import asdict
from math import prod
from pathlib import Path

import gdown
import pytest
import torch
from torch import Tensor

from boltz.data.feature.featurizer import BoltzFeaturizer
from boltz.data.load import CACHE_DIR
from boltz.data.load.load import load
from boltz.data.module.inference import BoltzInferenceDataModule
from boltz.data.tokenize.boltz import BoltzTokenizer
from boltz.data.types import Input, Manifest, Target, Tokenized
from boltz.distributed.data.types import PairMaskMode
from boltz.distributed.manager import DistributedManager
from boltz.main import (
    MOL_URL,
    BoltzDiffusionParams,
    BoltzProcessedInput,
    BoltzSteeringParams,
    parse_yaml,
)
from boltz.model.models.boltz1 import Boltz1
from boltz.testing.utils import concat_data, download_model_ckpt


def pytest_addoption(parser):
    """Add custom CLI options for pytest.

    Adds the --name_regex option to filter tests by matching their node IDs
    against a regular expression pattern.

    Args:
        parser: The pytest argument parser to add options to.
    """
    parser.addoption("--name_regex", action="store", default=None, help="Run tests matching a regex")


def pytest_collection_modifyitems(config, items):
    """Filter collected test items based on the --name_regex option.

    If --name_regex is provided, only tests whose node ID matches the given
    regular expression pattern will be selected for execution. Non-matching
    tests are deselected and reported via the pytest_deselected hook.

    Args:
        config: The pytest config object containing CLI options.
        items: List of collected test items to filter (modified in-place).
    """
    regex = config.getoption("--name_regex")
    if regex:
        r = re.compile(regex)
        selected = []
        deselected = []
        for item in items:
            if r.search(item.nodeid):
                selected.append(item)
            else:
                deselected.append(item)
        items[:] = selected
        config.hook.pytest_deselected(items=deselected)


ROOT_DIR = Path(__file__).resolve().parents[1]
EXAMPLE_DIR = ROOT_DIR / "examples"

EXAMPLE_PROT_YAML = EXAMPLE_DIR / "prot.yaml"
EXAMPLE_PROT_CUSTOM_MSA_YAML = EXAMPLE_DIR / "prot_custom_msa.yaml"
EXAMPLE_MULTIMER_YAML = EXAMPLE_DIR / "multimer.yaml"  # can lead to SIGKILL on cp = (4, 4)
EXAMPLE_CYCLIC_PROT_YAML = EXAMPLE_DIR / "cyclic_prot.yaml"
EXAMPLE_YAMLS = [
    EXAMPLE_PROT_YAML,
    EXAMPLE_MULTIMER_YAML,
]

TEST_INFERENCE_DIR = ROOT_DIR / "tests" / "data" / "inference"
TEST_YAML = TEST_INFERENCE_DIR / "test_input.yaml"
TEST_POCKET_YAMLS = [
    TEST_INFERENCE_DIR / "pocket_small.yaml",
]

SEED = 42
CCD_URL = "https://huggingface.co/boltz-community/boltz-1/resolve/main/ccd.pkl"


def download_ccd() -> Path:
    ccd_path = CACHE_DIR / "ccd.pkl"
    if not ccd_path.exists():
        with tempfile.TemporaryDirectory() as temp_dir:
            zip_path = os.path.join(temp_dir, "ccd.pkl")
            gdown.download(CCD_URL, zip_path)
            shutil.move(zip_path, ccd_path)

    return ccd_path


@pytest.fixture(scope="session")
def ccd():
    ccd_path = download_ccd()
    with ccd_path.open("rb") as file:
        return pickle.load(file)  # noqa: S301


@pytest.fixture(scope="function")
def example_cyclic_prot_input(ccd):
    target: Target = parse_yaml(EXAMPLE_CYCLIC_PROT_YAML, ccd)
    input = Input(target.structure, {}, target.record)
    return input


@pytest.fixture(scope="function")
def example_multimer_input(ccd):
    target: Target = parse_yaml(EXAMPLE_MULTIMER_YAML, ccd)
    input = Input(target.structure, {}, target.record)
    return input


@pytest.fixture(scope="function", params=EXAMPLE_YAMLS)
def example_input(ccd, request):
    yaml = request.param
    target: Target = parse_yaml(yaml, ccd)
    input = Input(target.structure, {}, target.record)
    return input


@pytest.fixture(scope="function")
def example_multimer_tokenized(example_multimer_input: Input) -> Tokenized:
    tokenizer = BoltzTokenizer()
    return tokenizer.tokenize(example_multimer_input)


@pytest.fixture(scope="function")
def example_tokenized(example_input: Input) -> Tokenized:
    tokenizer = BoltzTokenizer()
    return tokenizer.tokenize(example_input)


@pytest.fixture(scope="function")
def example_features(example_tokenized: Tokenized) -> dict[str, Tensor]:
    featurizer = BoltzFeaturizer()
    feats = featurizer.process(
        example_tokenized,
        training=False,
        augmentation=False,
        pair_mask_mode=PairMaskMode.SEQUENCE_LOCAL_ATTENTION,
    )
    return {k: v.unsqueeze(0) for k, v in feats.items()}  # create batch dimension


@pytest.fixture(scope="function")
def example_prot_input(ccd):
    target: Target = parse_yaml(EXAMPLE_PROT_YAML, ccd)
    input = Input(target.structure, {}, target.record)
    return input


@pytest.fixture(scope="function")
def example_prot_tokenized(example_prot_input: Input) -> Tokenized:
    tokenizer = BoltzTokenizer()
    return tokenizer.tokenize(example_prot_input)


@pytest.fixture(scope="function", params=TEST_POCKET_YAMLS)
def example_pocket_input(ccd, request):
    yaml = request.param
    target: Target = parse_yaml(yaml, ccd)
    input = Input(target.structure, {}, target.record)
    return input


@pytest.fixture(scope="function")
def setup_env(request, monkeypatch):
    (n_procs_dp, n_procs_cp), specify_method, device_type, method_init = request.param
    if isinstance(n_procs_cp, tuple) and all(isinstance(n, int) for n in n_procs_cp):
        world_size = n_procs_dp * prod(n_procs_cp)
    elif isinstance(n_procs_cp, int):
        world_size = n_procs_dp * n_procs_cp
    else:
        raise ValueError(f"Invalid type for CP ranks: {type(n_procs_cp)}")
    env_per_rank = None
    if specify_method:
        if method_init is None:
            raise ValueError("method_init must be specified if specify_method is True")
        # this emulates the behavior in the DistributedManager in the
        # case where the user explicitly sets the init method. In this
        # case, we don't need to clean up the env variables not in the
        # the scope of the "method_init" because only the "method_init"
        # in question is used to initialize the distributed manager
        monkeypatch.setenv("BOLTZ_DISTRIBUTED_INIT_METHOD", method_init)
    elif not specify_method:
        # this emulates the behavior in the DistributedManager in the
        # case where the user does not specify the init method but rely
        # on available environment variables to initialize the DistributedManager.
        # We need to clean up the existing environment variables so that
        # the tests truly respect the method_init. There are two sub-cases:
        # 1. the user does not have any environment variables set, in which case
        #    the DistributedManager will default initialize, which falls back
        #    to the single-device case
        # 2. the user has set the environment variables related to the
        #    ENV or SLURM init method, in which case the DistributedManager
        #    will use the environment variables to initialize
        # Sub-case 2 is handled by the code in the subsequent if statements
        monkeypatch.delenv("MASTER_ADDR", raising=False)
        monkeypatch.delenv("MASTER_PORT", raising=False)
        monkeypatch.delenv("WORLD_SIZE", raising=False)
        monkeypatch.delenv("RANK", raising=False)
        monkeypatch.delenv("LOCAL_RANK", raising=False)
        monkeypatch.delenv("SLURM_LAUNCH_NODE_IPADDR", raising=False)
        monkeypatch.delenv("SLURM_NPROCS", raising=False)
        monkeypatch.delenv("SLURM_PROCID", raising=False)
        monkeypatch.delenv("SLURM_LOCALID", raising=False)
    if method_init == "ENV":
        monkeypatch.setenv("MASTER_ADDR", "localhost")
        monkeypatch.setenv("MASTER_PORT", "29500")
        monkeypatch.setenv("WORLD_SIZE", f"{world_size}")
        env_per_rank = {"RANK": "<INPUT_RANK>", "LOCAL_RANK": "<INPUT_RANK>"}
    elif method_init == "SLURM":
        monkeypatch.setenv("SLURM_LAUNCH_NODE_IPADDR", "localhost")
        monkeypatch.setenv("SLURM_NPROCS", f"{world_size}")
        env_per_rank = {"SLURM_PROCID": "<INPUT_RANK>", "SLURM_LOCALID": "<INPUT_RANK>"}
    backend = DistributedManager.backend_for_device()[device_type]
    grid_group_sizes = OrderedDict(dp=n_procs_dp, cp=n_procs_cp)
    yield grid_group_sizes, world_size, device_type, backend, method_init, env_per_rank


@pytest.fixture(scope="session")
def get_inference_golden_value_dir_v1() -> Path:
    return load("unittests/test_inference_pipeline_golden_values")


@pytest.fixture(scope="session")
def test_cp_integration_data_dir_v1() -> Path:
    base_data_dir = load("unittests/test_cp_integration")
    return base_data_dir


@pytest.fixture(scope="session")
def test_cp_integration_data_dir_boltz1_v1() -> Path:
    base_data_dir = load("unittests/test_cp_integration_boltz1")
    return base_data_dir


@pytest.fixture(scope="session", params=["7ylz", "7z64", "8ayv", "8b2e"])
def get_preprocessed_boltz1_v1(test_cp_integration_data_dir_boltz1_v1, request):
    name = request.param
    return test_cp_integration_data_dir_boltz1_v1 / f"processed_{name}"


@pytest.fixture(scope="session", params=["7ylz", "7z64", "8ayv", "8b2e"])
def get_preprocessed_v1(test_cp_integration_data_dir_v1, request):
    name = request.param
    return test_cp_integration_data_dir_v1 / f"processed_{name}"


@pytest.fixture(scope="session")
def create_preprocessed_handle_boltz1_v1(
    get_preprocessed_boltz1_v1: Path,
) -> BoltzProcessedInput:
    f_manifest = get_preprocessed_boltz1_v1 / "manifest.json"
    dir_structure = get_preprocessed_boltz1_v1 / "structures"
    dir_msa = get_preprocessed_boltz1_v1 / "msa"
    assert f_manifest.is_file(), f"Manifest file {f_manifest} does not exist"
    assert dir_structure.is_dir(), f"Structure directory {dir_structure} does not exist"
    assert dir_msa.is_dir(), f"MSA directory {dir_msa} does not exist"
    processed = BoltzProcessedInput(
        manifest=Manifest.load(f_manifest),
        targets_dir=dir_structure,
        msa_dir=dir_msa,
    )
    return processed


@pytest.fixture(scope="session")
def create_preprocessed_handle_v1(get_preprocessed_v1: Path) -> BoltzProcessedInput:
    f_manifest = get_preprocessed_v1 / "manifest.json"
    dir_structure = get_preprocessed_v1 / "structures"
    dir_msa = get_preprocessed_v1 / "msa"
    assert f_manifest.is_file(), f"Manifest file {f_manifest} does not exist"
    assert dir_structure.is_dir(), f"Structure directory {dir_structure} does not exist"
    assert dir_msa.is_dir(), f"MSA directory {dir_msa} does not exist"
    processed = BoltzProcessedInput(
        manifest=Manifest.load(f_manifest),
        targets_dir=dir_structure,
        msa_dir=dir_msa,
    )
    return processed


@pytest.fixture(scope="function")
def create_datamodule_serial_v1(create_preprocessed_handle_v1):
    data_module = BoltzInferenceDataModule(
        manifest=create_preprocessed_handle_v1.manifest,
        target_dir=create_preprocessed_handle_v1.targets_dir,
        msa_dir=create_preprocessed_handle_v1.msa_dir,
        augmentation=False,
        num_workers=0,
    )  # default use_cache=False
    return data_module


@pytest.fixture(scope="session")
def get_model_v1_ckpt():
    f_ckpt = download_model_ckpt()
    assert f_ckpt.is_file(), f"Checkpoint file {f_ckpt} does not exist"
    return f_ckpt


@pytest.fixture(scope="session")
def get_predict_args_v1():
    predict_args = {
        "recycling_steps": 10,  # Boltz uses 10 for evaluation (https://github.com/jwohlwend/boltz/blob/main/docs/evaluation.md#evaluation-setup)
        "sampling_steps": 200,
        "diffusion_samples": 1,
        "write_confidence_summary": False,
        "write_full_pae": False,
        "write_full_pde": False,
    }
    return predict_args


@pytest.fixture(scope="session")
def get_diffusion_params_v1():
    diffusion_params = BoltzDiffusionParams()
    return asdict(diffusion_params)


@pytest.fixture(scope="session")
def get_steering_params_no_potentials_v1():
    steering_params = BoltzSteeringParams()
    steering_params.fk_steering = False
    steering_params.guidance_update = False
    return asdict(steering_params)


@pytest.fixture(scope="session")
def get_score_model_args_v1():
    return {
        "sigma_data": 16,
        "dim_fourier": 256,
        "atom_encoder_depth": 3,
        "atom_encoder_heads": 4,
        "token_transformer_depth": 24,
        "token_transformer_heads": 16,
        "atom_decoder_depth": 3,
        "atom_decoder_heads": 4,
        "activation_checkpointing": False,
    }


@pytest.fixture(scope="session")
def _load_model_v1_with_ckpt_serial(
    get_model_v1_ckpt,
    get_predict_args_v1,
    get_diffusion_params_v1,
    get_steering_params_no_potentials_v1,
    get_score_model_args_v1,
):
    kwargs_model = {
        "confidence_prediction": False,
        "predict_args": get_predict_args_v1,
        "diffusion_process_args": get_diffusion_params_v1,
        "steering_args": get_steering_params_no_potentials_v1,
        "score_model_args": get_score_model_args_v1,
        "ema": False,
    }
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        model = Boltz1.load_from_checkpoint(get_model_v1_ckpt, strict=False, map_location="cpu", **kwargs_model)
    ckpt_dict = torch.load(get_model_v1_ckpt, map_location="cpu", weights_only=False)
    # merge the model kwargs into the hyper_parameters entry
    ckpt_dict["hyper_parameters"].update(kwargs_model)
    # drop deprecated args
    ckpt_dict["hyper_parameters"].pop("chain_sampling_args")
    ckpt_dict["hyper_parameters"].pop("recycling_detach")
    ckpt_dict["hyper_parameters"].pop("run_trunk_and_structure")
    # the CP-related settings must be absent to prevent mix-and-matching
    # data sharding and model distributed computation settings
    assert "dist_manager" not in ckpt_dict["hyper_parameters"]
    assert "layout_map_cp" not in ckpt_dict["hyper_parameters"]
    return model, ckpt_dict


@pytest.fixture(scope="function")
def load_model_v1_with_ckpt_serial(_load_model_v1_with_ckpt_serial):
    _, ckpt_dict = _load_model_v1_with_ckpt_serial
    ckpt_dict = deepcopy(ckpt_dict)
    model = Boltz1(**ckpt_dict["hyper_parameters"])
    model.load_state_dict(ckpt_dict["state_dict"], strict=False)
    return model, ckpt_dict


@pytest.fixture(scope="function", params=[False, True], ids=["no_tf32", "tf32"])
def setup_tf32(request):
    """Configure TF32 settings and reset them after the test."""
    use_tf32 = request.param

    original_env = os.environ.get("NVIDIA_TF32_OVERRIDE", None)
    original_matmul_tf32 = torch.backends.cuda.matmul.allow_tf32
    original_cudnn_tf32 = torch.backends.cudnn.allow_tf32

    tf32_value = "1" if use_tf32 else "0"
    os.environ["NVIDIA_TF32_OVERRIDE"] = tf32_value
    torch.backends.cuda.matmul.allow_tf32 = use_tf32
    torch.backends.cudnn.allow_tf32 = use_tf32

    yield use_tf32

    if original_env is not None:
        os.environ["NVIDIA_TF32_OVERRIDE"] = original_env
    else:
        os.environ.pop("NVIDIA_TF32_OVERRIDE", None)
    torch.backends.cuda.matmul.allow_tf32 = original_matmul_tf32
    torch.backends.cudnn.allow_tf32 = original_cudnn_tf32


@pytest.fixture
def get_training_data_v1():
    """Download and return the path to the training data truncated set.

    Returns the path to boltz_training_truncated_set directory which contains:
    - openfold_processed_targets/
    - openfold_processed_msa/
    - train_ids.txt
    - validation_ids.txt

    The parent directory (training_data) contains symmetry.pkl.
    """
    return load("unittests/training_data_truncated_set") / "training_data" / "boltz_training_truncated_set"


####################################################################################################
# Boltz 2 UTILITIES
####################################################################################################


def _build_processed_input_boltz2(base_dir: Path) -> BoltzProcessedInput:
    """Build a BoltzProcessedInput from a processed data directory."""
    manifest = Manifest.load(base_dir / "manifest.json")
    targets_dir = base_dir / "structures"
    msa_dir = base_dir / "msa"
    template_dir = base_dir / "templates" if (base_dir / "templates").exists() else None
    extra_mols_dir = base_dir / "extra_mols" if (base_dir / "extra_mols").exists() else None
    return BoltzProcessedInput(
        manifest=manifest,
        targets_dir=targets_dir,
        msa_dir=msa_dir,
        constraints_dir=None,
        template_dir=template_dir,
        extra_mols_dir=extra_mols_dir,
    )


def _concat_data_with_records(out_dir: Path, *datas: Path) -> Path:
    """Concatenate processed Boltz2 dirs, including records."""
    merged = concat_data(out_dir, *datas)
    records_dir = merged / "records"
    records_dir.mkdir(parents=True, exist_ok=True)

    copied_records = set()
    for data in datas:
        src_records_dir = Path(data) / "records"
        for record_file in src_records_dir.glob("*.json"):
            if record_file.name in copied_records:
                raise ValueError(f"Duplicate record file {record_file.name}")
            shutil.copy(record_file, records_dir / record_file.name)
            copied_records.add(record_file.name)

    return merged


@pytest.fixture(scope="session")
def canonical_mols_dir() -> Path:
    """Download canonical molecules to cache if needed, return the mols directory.

    The mols directory contains per-residue pickle files of RDKit Mol objects derived
    from the Chemical Component Dictionary (CCD). Each file (e.g. ALA.pkl) provides
    the reference 3D structure, atom names, bonds, and other chemical metadata for a
    single CCD component.

    This is the Boltz2 equivalent of the `ccd.pkl` file used by Boltz1. Where
    Boltz1 deserializes the entire CCD dictionary into memory at once, Boltz2 splits it
    into individual per-component files so that only the needed subset is loaded.
    `load_canonicals` eagerly loads the 20 standard amino acid residues plus UNK,
    while `load_molecules` and `get_mol` load additional non-standard components
    lazily on demand during tokenization, featurization, and structure parsing.

    Returns
    -------
    Path
        The path to the mols directory.
    """
    cache_mols = CACHE_DIR / "mols"
    cache_tar = CACHE_DIR / "mols.tar"
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    if not cache_tar.exists():
        with tempfile.TemporaryDirectory() as temp_dir:
            tmp_tar = os.path.join(temp_dir, "mols.tar")
            gdown.download(MOL_URL, tmp_tar)
            shutil.move(tmp_tar, cache_tar)
    if not cache_mols.exists():
        with tarfile.open(str(cache_tar), "r") as tar:
            tar.extractall(CACHE_DIR)  # noqa: S202
    return cache_mols


@pytest.fixture(scope="session")
def test_cp_training_base_data_dir_boltz2() -> Path:
    """Raw Boltz2 training archive root containing processed_{id} dirs."""
    return load("unittests/test_cp_training_data_boltz2")


@pytest.fixture(scope="session", params=["processed_7ylz", "processed_7z64", "processed_8ayv", "processed_8b2e"])
def get_preprocessed_boltz2(test_cp_training_base_data_dir_boltz2: Path, request: pytest.FixtureRequest) -> Path:
    """Per-sample preprocessed directory for Boltz-2 predict tests."""
    return test_cp_training_base_data_dir_boltz2 / request.param


@pytest.fixture(scope="session")
def test_cp_training_data_dir_boltz2(
    tmp_path_factory: pytest.TempPathFactory,
    test_cp_training_base_data_dir_boltz2: Path,
) -> Path:
    """Merged Boltz2 training directory with records for the training data module."""
    names = ["7ylz", "7z64", "8ayv", "8b2e"]
    source_dirs = [test_cp_training_base_data_dir_boltz2 / f"processed_{name}" for name in names]
    out_dir = tmp_path_factory.mktemp("cp_training_data_boltz2") / "processed_training"
    return _concat_data_with_records(out_dir, *source_dirs)


@pytest.fixture(scope="session")
def create_preprocessed_handle_boltz2(test_cp_training_data_dir_boltz2: Path) -> BoltzProcessedInput:
    """Build a BoltzProcessedInput from the merged Boltz2 training data directory."""
    return _build_processed_input_boltz2(test_cp_training_data_dir_boltz2)


@pytest.fixture(scope="session")
def get_inference_golden_value_dir_v2() -> Path:
    return load("unittests/test_inference_pipeline_golden_values_boltz2")


@pytest.fixture(scope="session")
def get_model_ckpt_v2() -> Path:
    from boltz.main import BOLTZ2_URL_WITH_FALLBACK

    cache = CACHE_DIR
    cache.mkdir(parents=True, exist_ok=True)
    checkpoint = cache / "boltz2_conf.ckpt"
    if not checkpoint.exists():
        for i, url in enumerate(BOLTZ2_URL_WITH_FALLBACK):
            try:
                urllib.request.urlretrieve(url, str(checkpoint))  # noqa: S310
                break
            except Exception as e:  # noqa: BLE001
                if i == len(BOLTZ2_URL_WITH_FALLBACK) - 1:
                    msg = f"Failed to download Boltz-2 checkpoint: {e}"
                    raise RuntimeError(msg) from e
    assert checkpoint.is_file(), f"Checkpoint {checkpoint} does not exist"
    return checkpoint
