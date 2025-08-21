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

from boltz.data.module.training_polymer import (
    PolymerCSVDataModule,
    DataConfig as PolymerCSVDataConfig
)
from boltz.main import download_boltz2, get_cache_path
"""
Tokenizer and featurizer are provided via YAML (Hydra _target_), no direct imports here.
"""


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
    default_pretrained: Optional[bool] = False


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

    # Mirror user-provided polymer_property_keys into model.affinity_model_args if needed
    try:
        if raw_config.get("model", {}).get("polymer_only_training", False) and (
            "csv_path" in raw_config.get("data", {})
        ):
            if raw_config.get("data", {}).get("polymer_property_keys", None):
                model_args = raw_config.get("model", {})
                affinity_args = model_args.get("affinity_model_args", {}) or {}
                if "polymer_property_keys" not in affinity_args:
                    affinity_args["polymer_property_keys"] = (
                        raw_config.data.polymer_property_keys
                    )
                    raw_config.model.affinity_model_args = affinity_args
    except Exception:
        pass

    # Instantiate the task
    cfg = hydra.utils.instantiate(raw_config)
    cfg = TrainConfig(**cfg)

    # If requested, use the default downloaded Boltz2 checkpoint when no pretrained is provided
    if cfg.default_pretrained and not cfg.pretrained:
        cache_dir = Path(get_cache_path()).expanduser()
        cache_dir.mkdir(parents=True, exist_ok=True)
        # Ensure the default checkpoint is available
        download_boltz2(cache_dir)
        cfg.pretrained = str(cache_dir / "boltz2_conf.ckpt")

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
    # Always use CSV DataConfig/DataModule
    if "csv_path" not in cfg.data or not cfg.data["csv_path"]:
        raise ValueError("data.csv_path must be provided for CSV training")
    data_config = PolymerCSVDataConfig(**cfg.data)
    data_module = PolymerCSVDataModule(data_config)

    if cfg.pretrained and not cfg.resume:
        # Load and optionally modify the checkpoint before loading
        checkpoint = torch.load(cfg.pretrained, map_location="cpu")
        modified = False

        # Optionally broadcast trunk weights into confidence module
        if cfg.load_confidence_from_trunk:
            new_state_dict = {}
            for key, value in checkpoint["state_dict"].items():
                if not key.startswith("structure_module") and not key.startswith(
                    "distogram_module"
                ):
                    new_key = "confidence_module." + key
                    new_state_dict[new_key] = value
            new_state_dict.update(checkpoint["state_dict"])
            checkpoint["state_dict"] = new_state_dict
            modified = True

        # If polymer-only training, remap affinity head weights into polymer head if present
        try:
            is_polymer_only = getattr(model_module, "polymer_only_training", False)
            if is_polymer_only:
                sd = checkpoint["state_dict"]
                has_aff = any(k.startswith("affinity_module.") for k in sd)
                # If the pretrained has an affinity_module, map it to polymer_module
                if has_aff:
                    remapped = {}
                    for k, v in sd.items():
                        if k.startswith("affinity_module."):
                            new_key = k.replace("affinity_module.", "polymer_module.", 1)
                            remapped[new_key] = v
                        else:
                            remapped[k] = v
                    checkpoint["state_dict"] = remapped
                    modified = True
        except Exception:
            pass

        # Build a filtered state_dict that matches current model shapes
        src_sd = checkpoint.get("state_dict", {})
        tgt_sd = model_module.state_dict()
        filtered = {}
        skipped = []
        for k, v in src_sd.items():
            if k in tgt_sd and getattr(v, "shape", None) == getattr(tgt_sd[k], "shape", None):
                filtered[k] = v
            else:
                skipped.append(k)
        print(
            f"Loading pretrained weights: matching={len(filtered)} skipped={len(skipped)}"
        )
        missing, unexpected = model_module.load_state_dict(filtered, strict=False)
        if missing:
            print(f"Missing keys (not loaded): {len(missing)}")
        if unexpected:
            print(f"Unexpected keys (ignored): {len(unexpected)}")

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
