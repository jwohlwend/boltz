import os
import random
import string
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Any

import hydra
import omegaconf
import pytorch_lightning as pl
import torch
import torch.multiprocessing
from omegaconf import OmegaConf, listconfig
from pytorch_lightning import LightningModule
from pytorch_lightning.callbacks.model_checkpoint import ModelCheckpoint
from pytorch_lightning.loggers import WandbLogger
from pytorch_lightning.strategies import DDPStrategy
from pytorch_lightning.utilities import rank_zero_only

from boltz.data.module.polymer_trainingv2 import (
    PolymerTrainingDataModule,
    DataConfig as PolymerDataConfig,
)
from boltz.data.module.polymer_csv import PolymerCSVDataModule, DataConfig as PolymerCSVDataConfig
from boltz.data.tokenize.boltz2 import Boltz2Tokenizer
from boltz.data.feature.featurizerv2 import Boltz2Featurizer


@dataclass
class TrainConfig:
    """Train configuration.

    Attributes
    ----------
    data : DataConfig
        The data configuration.
    model : ModelConfig
        The model configuration.
    output : str
        The output directory.
    trainer : Optional[dict]
        The trainer configuration.
    resume : Optional[str]
        The resume checkpoint.
    pretrained : Optional[str]
        The pretrained model.
    wandb : Optional[dict]
        The wandb configuration.
    disable_checkpoint : bool
        Disable checkpoint.
    matmul_precision : Optional[str]
        The matmul precision.
    find_unused_parameters : Optional[bool]
        Find unused parameters.
    save_top_k : Optional[int]
        Save top k checkpoints.
    validation_only : bool
        Run validation only.
    debug : bool
        Debug mode.
    strict_loading : bool
        Fail on mismatched checkpoint weights.
    load_confidence_from_trunk: Optional[bool]
        Load pre-trained confidence weights from trunk.

    """

    data: Any
    model: LightningModule
    output: str
    trainer: Optional[dict] = None
    resume: Optional[str] = None
    pretrained: Optional[str] = None
    wandb: Optional[dict] = None
    disable_checkpoint: bool = False
    matmul_precision: Optional[str] = None
    find_unused_parameters: Optional[bool] = False
    save_top_k: Optional[int] = 1
    validation_only: bool = False
    debug: bool = False
    strict_loading: bool = True
    load_confidence_from_trunk: Optional[bool] = False


def train(raw_config: str, args: list[str]) -> None:  # noqa: C901, PLR0912, PLR0915
    """Run training.

    Parameters
    ----------
    raw_config : str
        The input yaml configuration.
    args : list[str]
        Any command line overrides.

    """
    # Load the configuration
    raw_config = omegaconf.OmegaConf.load(raw_config)

    # Apply input arguments
    args = omegaconf.OmegaConf.from_dotlist(args)
    raw_config = omegaconf.OmegaConf.merge(raw_config, args)

    # Instantiate the task
    cfg = hydra.utils.instantiate(raw_config)
    cfg = TrainConfig(**cfg)

    # Set matmul precision
    if cfg.matmul_precision is not None:
        torch.set_float32_matmul_precision(cfg.matmul_precision)

    # Create trainer dict
    trainer = cfg.trainer
    if trainer is None:
        trainer = {}

    # Flip some arguments in debug mode
    devices = trainer.get("devices", 1)

    wandb = cfg.wandb
    if cfg.debug:
        if isinstance(devices, int):
            devices = 1
        elif isinstance(devices, (list, listconfig.ListConfig)):
            devices = [devices[0]]
        trainer["devices"] = devices
        cfg.data.num_workers = 0
        if wandb:
            wandb = None

    # Create objects
    model_module = cfg.model
    polymer_only = getattr(model_module, "polymer_only_training", False)
    if polymer_only:
        # Allow CSV ingestion if user provides csv_path in data
        if "csv_path" in cfg.data:
            ds = cfg.data
            csv_cfg = {
                "csv_path": ds["csv_path"],
                # For CSV/SMILES + StructureV2, use v2 tokenizer/featurizer regardless of YAML
                "tokenizer": Boltz2Tokenizer(),
                "featurizer": Boltz2Featurizer(),
                "samples_per_epoch": ds["samples_per_epoch"],
                "batch_size": ds["batch_size"],
                "num_workers": ds["num_workers"],
                "pin_memory": ds["pin_memory"],
                "random_seed": ds["random_seed"],
                "symmetries": ds.get("symmetries", None),
                "max_atoms": ds["max_atoms"],
                "max_tokens": ds["max_tokens"],
                "max_seqs": ds["max_seqs"],
                "pad_to_max_tokens": ds.get("pad_to_max_tokens", False),
                "pad_to_max_atoms": ds.get("pad_to_max_atoms", False),
                "pad_to_max_seqs": ds.get("pad_to_max_seqs", False),
                "crop_validation": ds.get("crop_validation", False),
                "return_train_symmetries": ds.get("return_train_symmetries", False),
                "return_val_symmetries": ds.get("return_val_symmetries", True),
                "atoms_per_window_queries": ds.get("atoms_per_window_queries", 32),
                "min_dist": ds.get("min_dist", 2.0),
                "max_dist": ds.get("max_dist", 22.0),
                "num_bins": ds.get("num_bins", 64),
                "val_batch_size": ds.get("val_batch_size", 1),
            }
            data_config = PolymerCSVDataConfig(**csv_cfg)
            data_module = PolymerCSVDataModule(data_config)
        else:
            data_config = PolymerDataConfig(**cfg.data)
            data_module = PolymerTrainingDataModule(data_config)
    else:
        data_config = FullDataConfig(**cfg.data)
        data_module = BoltzTrainingDataModule(data_config)

    if cfg.pretrained and not cfg.resume:
        # Load the pretrained weights into the confidence module
        if cfg.load_confidence_from_trunk:
            checkpoint = torch.load(cfg.pretrained, map_location="cpu")

            # Modify parameter names in the state_dict
            new_state_dict = {}
            for key, value in checkpoint["state_dict"].items():
                if not key.startswith("structure_module") and not key.startswith(
                    "distogram_module"
                ):
                    new_key = "confidence_module." + key
                    new_state_dict[new_key] = value
            new_state_dict.update(checkpoint["state_dict"])

            # Update the checkpoint with the new state_dict
            checkpoint["state_dict"] = new_state_dict

            # Save the modified checkpoint
            random_string = "".join(
                random.choices(string.ascii_lowercase + string.digits, k=10)
            )
            file_path = os.path.dirname(cfg.pretrained) + "/" + random_string + ".ckpt"
            print(
                f"Saving modified checkpoint to {file_path} created by broadcasting trunk of {cfg.pretrained} to confidence module."
            )
            torch.save(checkpoint, file_path)
        else:
            file_path = cfg.pretrained

        print(f"Loading model from {file_path}")
        model_module = type(model_module).load_from_checkpoint(
            file_path, map_location="cpu", strict=False, **(model_module.hparams)
        )

        if cfg.load_confidence_from_trunk:
            os.remove(file_path)

    # Create checkpoint callback
    callbacks = []
    dirpath = cfg.output
    if not cfg.disable_checkpoint:
        # In polymer-only training we don't compute lddt. Monitor polymer properties loss instead.
        is_polymer_only = getattr(model_module, "polymer_only_training", False)
        if is_polymer_only:
            mc = ModelCheckpoint(
                monitor="train/polymer_properties_loss",
                save_top_k=cfg.save_top_k,
                save_last=True,
                mode="min",
                every_n_epochs=1,
            )
        else:
            mc = ModelCheckpoint(
                monitor="val/lddt",
                save_top_k=cfg.save_top_k,
                save_last=True,
                mode="max",
                every_n_epochs=1,
            )
        callbacks = [mc]

    # Create wandb logger
    loggers = []
    if wandb:
        wdb_logger = WandbLogger(
            name=wandb["name"],
            group=wandb["name"],
            save_dir=cfg.output,
            project=wandb["project"],
            entity=wandb["entity"],
            log_model=False,
        )
        loggers.append(wdb_logger)
        # Save the config to wandb

        @rank_zero_only
        def save_config_to_wandb() -> None:
            config_out = Path(wdb_logger.experiment.dir) / "run.yaml"
            with Path.open(config_out, "w") as f:
                OmegaConf.save(raw_config, f)
            wdb_logger.experiment.save(str(config_out))

        save_config_to_wandb()

    # Set up trainer
    strategy = "auto"
    if (isinstance(devices, int) and devices > 1) or (
        isinstance(devices, (list, listconfig.ListConfig)) and len(devices) > 1
    ):
        strategy = DDPStrategy(find_unused_parameters=cfg.find_unused_parameters)

    trainer = pl.Trainer(
        default_root_dir=str(dirpath),
        strategy=strategy,
        callbacks=callbacks,
        logger=loggers,
        enable_checkpointing=not cfg.disable_checkpoint,
        reload_dataloaders_every_n_epochs=1,
        **trainer,
    )

    if not cfg.strict_loading:
        model_module.strict_loading = False

    if cfg.validation_only:
        trainer.validate(
            model_module,
            datamodule=data_module,
            ckpt_path=cfg.resume,
        )
    else:
        trainer.fit(
            model_module,
            datamodule=data_module,
            ckpt_path=cfg.resume,
        )


if __name__ == "__main__":
    arg1 = sys.argv[1]
    arg2 = sys.argv[2:]
    train(arg1, arg2)
