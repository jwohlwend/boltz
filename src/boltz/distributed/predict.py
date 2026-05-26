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

"""Distributed inference entrypoint for Boltz-2 with DTensor context parallelism.

Adapted from the Boltz-1x-CP ``distributed/main.py::run_predict`` function.
Differences from the Boltz-1 version:

- Uses :class:`Boltz2` (serial) + :class:`Boltz2Distributed` wrapper.
- Uses :class:`Boltz2InferenceDataModuleDTensor` (v2 featurizer/tokenizer).
- Requires ``mol_dir`` for per-residue CCD molecule files.
- Precision default is ``bf16-mixed`` (matching serial Boltz-2 inference).
"""

from __future__ import annotations

import atexit
import warnings
from collections import OrderedDict
from dataclasses import asdict
from datetime import timedelta
from math import isqrt
from pathlib import Path
from typing import Any, Literal, Optional

import torch
from lightning_fabric.plugins.precision.utils import _convert_fp_tensor
from lightning_utilities import apply_to_collection
from pytorch_lightning import Trainer, seed_everything
from pytorch_lightning.callbacks import Callback
from pytorch_lightning.plugins.precision import HalfPrecision
from pytorch_lightning.strategies import SingleDeviceStrategy

from boltz.data import const
from boltz.data.types import Manifest
from boltz.data.write.writer import BoltzWriter
from boltz.distributed.data.module.inferencev2 import Boltz2InferenceDataModuleDTensor
from boltz.distributed.data.types import PairMaskMode
from boltz.distributed.data.utils import map_subgroup_mesh_to_cpu
from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.layers.triangular_attention import can_run_cueq_triattn_sm100f
from boltz.distributed.model.models.boltz2 import Boltz2 as Boltz2Distributed
from boltz.distributed.model.modules.utils import (
    PRECISION_TO_LIGHTNING,
    Precision,
    SDPAWithBiasBackend,
    SetAttnPairBiasBackend,
    SetAttnPairBiasShardwiseBackend,
    SetTriAttnBackend,
    TriAttnBackend,
    setup_tf32_env,
)
from boltz.main import (
    Boltz2DiffusionParams,
    BoltzProcessedInput,
    check_inputs,
    filter_inputs_structure,
    process_inputs,
)
from boltz.model.models.boltz2 import Boltz2 as Boltz2Serial


class HalfPrecisionAllowFrozen(HalfPrecision):
    """Half precision plugin that handles Boltz-2's frozen dataclasses.

    Lightning's default ``HalfPrecision`` plugin raises ``ValueError`` when
    ``apply_to_collection`` encounters a frozen dataclass in the batch.  This
    subclass overrides ``convert_input`` to pass ``allow_frozen=True``, which
    lets Lightning recurse into frozen dataclass fields without trying to
    ``setattr`` on the frozen instance.
    """

    def convert_input(self, data: Any) -> Any:
        return apply_to_collection(
            data,
            function=_convert_fp_tensor,
            dtype=torch.Tensor,
            dst_type=self._desired_input_dtype,
            allow_frozen=True,
        )


def run_predict(
    data: str | Path,
    out_dir: str | Path,
    mol_dir: str | Path,
    checkpoint: str | Path,
    size_dp: int = 1,
    size_cp: int = 1,
    accelerator: str = "gpu",
    recycling_steps: int = 3,
    sampling_steps: int = 200,
    diffusion_samples: int = 1,
    max_parallel_samples: Optional[int] = None,
    step_scale: float = 1.5,
    output_format: Literal["pdb", "mmcif"] = "mmcif",
    seed: Optional[int] = None,
    max_msa_seqs: int = const.max_msa_seqs,
    msa_pad_to_max_seqs: bool = False,
    input_format: str = "preprocessed",
    timeout_nccl: Optional[float] = 30,
    timeout_gloo: Optional[float] = 30,
    precision: Precision = Precision.BF16_MIXED,
    atoms_per_window_queries_keys: tuple[Optional[int], Optional[int]] = (32, 128),
    pair_mask_mode: PairMaskMode = PairMaskMode.NONE,
    local_batch_size: int = 1,
    num_ensembles: int = 1,
    extra_callbacks: Optional[list[Callback]] = None,
    confidence_prediction: bool = True,
    write_full_pae: bool = False,
    use_templates: bool = True,
    triattn_backend: TriAttnBackend = TriAttnBackend.CUEQ,
    sdpa_with_bias_backend: SDPAWithBiasBackend = SDPAWithBiasBackend.TORCH_FLEX_ATTN,
    sdpa_with_bias_shardwise_backend: SDPAWithBiasBackend = SDPAWithBiasBackend.TORCH_FLEX_ATTN,
    auto_pad_tokens_for_sm100f: bool = True,
    cuda_memory_profile: bool = False,
    override: bool = False,
    max_data_retries: int = 5,
) -> None:
    """Run distributed Boltz-2 structure prediction with DTensor context parallelism.

    Parameters
    ----------
    data : str or Path
        Path to the input data. For ``input_format="preprocessed"``, a
        directory containing manifest.json, structures/, msa/. For
        ``input_format="config_files"``, a YAML/FASTA file or directory
        of such files.
    out_dir : str or Path
        Output directory for predictions.
    mol_dir : str or Path
        Directory containing per-residue CCD molecule pickle files.
    checkpoint : str or Path
        Path to the Boltz-2 model checkpoint.
    size_dp : int
        Number of data-parallel ranks.
    size_cp : int
        Total number of context-parallel ranks (must be a perfect square).
    accelerator : str
        Device accelerator ("gpu" or "cpu").
    recycling_steps : int
        Number of recycling iterations for the trunk.
    sampling_steps : int
        Number of diffusion denoising steps.
    diffusion_samples : int
        Number of independent diffusion samples per input.
    max_parallel_samples : int or None
        Max diffusion samples to run in parallel (None = all at once).
    step_scale : float
        Step scale for the diffusion schedule.
    output_format : str
        Output structure format ("pdb" or "mmcif").
    seed : int or None
        Random seed for reproducibility.
    max_msa_seqs : int
        Maximum number of MSA sequences.
    msa_pad_to_max_seqs : bool
        Whether to pad MSA to max_msa_seqs.
    input_format : str
        Input data format ("preprocessed" or "config_files").
    timeout_nccl : float or None
        NCCL timeout in minutes.
    timeout_gloo : float or None
        Gloo timeout in minutes.
    precision : Precision
        Model precision mode.
    atoms_per_window_queries_keys : tuple
        (queries, keys) window sizes for atom attention batching.
    pair_mask_mode : PairMaskMode
        Pair mask mode (NONE = window batching).
    local_batch_size : int
        Per-rank batch size.
    num_ensembles : int
        Number of ensemble members for structure prediction.
    extra_callbacks : list[Callback] or None
        Additional Lightning callbacks.
    confidence_prediction : bool
        Whether to run the confidence module (pLDDT, pTM, iPTM, PAE, PDE).
    write_full_pae : bool
        Whether to write full PAE matrices (requires confidence module).
    use_templates : bool
        Reserved for future use. Currently ignored — template weights are
        always loaded from the checkpoint but the distributed TemplateModule
        is not yet implemented, so templates are skipped during forward.
    triattn_backend : TriAttnBackend
        Backend for distributed triangle attention in PairformerLayer.
    sdpa_with_bias_backend : SDPAWithBiasBackend
        Backend for ring-attention AttentionPairBias layers.
    sdpa_with_bias_shardwise_backend : SDPAWithBiasBackend
        Backend for window-batched AttentionPairBiasShardwise layers.
    auto_pad_tokens_for_sm100f : bool
        When True and the cuEq TriAttn backend is selected with BF16 or
        BF16_MIXED precision on an SM100/SM103 GPU, pad token counts so
        that each CP shard has a multiple-of-8 sequence length, enabling
        the SM100f cuEq TriAttn kernel.
    cuda_memory_profile : bool
        When True, records CUDA memory history and dumps a per-rank
        snapshot pickle to ``out_dir`` at the end of prediction.
    override : bool
        When True, rerun predictions even if output already exists.
    max_data_retries : int
        Maximum number of retry attempts when a data sample fails to load.
        Set to 0 to raise immediately on the first error.

    """
    atoms_per_window_queries, atoms_per_window_keys = atoms_per_window_queries_keys

    accelerator_to_device = {"cpu": "cpu", "gpu": "cuda"}
    if accelerator not in accelerator_to_device:
        raise ValueError(f"Accelerator {accelerator!r} not recognised; expected one of {sorted(accelerator_to_device)}")
    device_type = accelerator_to_device[accelerator]

    if device_type == "cuda" and not torch.cuda.is_available():
        raise ValueError("accelerator='gpu' requires CUDA, but torch.cuda.is_available() is False")

    if triattn_backend in (TriAttnBackend.CUEQ, TriAttnBackend.TRIFAST) and (
        device_type != "cuda" or not torch.cuda.is_available()
    ):
        raise ValueError(
            f"triattn_backend={triattn_backend.value!r} requires CUDA, "
            f"but accelerator={accelerator!r} (device_type={device_type!r}), "
            f"torch.cuda.is_available()={torch.cuda.is_available()}"
        )

    if not isinstance(size_cp, int) or size_cp <= 0:
        raise TypeError("size_cp must be a positive integer")
    if timeout_nccl is not None and timeout_nccl <= 0:
        raise TypeError("timeout_nccl must be a positive float")
    if timeout_gloo is not None and timeout_gloo <= 0:
        raise TypeError("timeout_gloo must be a positive float")
    if precision not in PRECISION_TO_LIGHTNING:
        raise ValueError(f"Precision {precision} not supported")

    timeout_nccl_td = timedelta(minutes=timeout_nccl) if timeout_nccl is not None else None
    timeout_gloo_td = timedelta(minutes=timeout_gloo) if timeout_gloo is not None else None
    timeout_by_device = {"cpu": timeout_gloo_td, "cuda": timeout_nccl_td}

    # --- Distributed setup ---
    DistributedManager.initialize(device_type=device_type, timeout=timeout_by_device[device_type])
    atexit.register(DistributedManager.cleanup)
    dist_manager = DistributedManager()

    if size_dp * size_cp != dist_manager.world_size:
        raise ValueError(f"world_size mismatch: {dist_manager.world_size} != size_dp*size_cp ({size_dp}*{size_cp})")

    size_cp_axis = isqrt(size_cp)
    if size_cp_axis * size_cp_axis != size_cp:
        raise ValueError(f"size_cp must be a perfect square, got {size_cp}")

    # Load checkpoint hparams to read existing pairformer_args / msa_args,
    # then ensure the critical V2 flags are set without overriding model
    # dimensions (num_heads, num_blocks, msa_s, etc.) with production defaults.
    #
    # Without the v2 flag, Boltz2Serial creates V1 attention layers (which use
    # a different norm layout — 128 extra `norm_s` LayerNorms per pairformer
    # layer).  Under strict=False this was silent: the V1 norm parameters were
    # randomly initialised and the V2 checkpoint weights were simply dropped,
    # producing garbage predictions.  We now merge the flag from the checkpoint
    # and load with strict=True to catch any mismatch.
    ckpt_raw = torch.load(str(checkpoint), map_location="cpu", weights_only=False, mmap=True)
    ckpt_hp = ckpt_raw.get("hyper_parameters", {})
    del ckpt_raw

    pairformer_args = dict(ckpt_hp.get("pairformer_args", {}))
    pairformer_args.setdefault("v2", True)

    msa_args = dict(ckpt_hp.get("msa_args", {}))
    msa_args.setdefault("use_paired_feature", True)

    dim_hidden = pairformer_args.get("pairwise_head_width", 32)
    sm100f_per_shard_token_multiple = 1
    if auto_pad_tokens_for_sm100f and triattn_backend == TriAttnBackend.CUEQ:
        if precision in (Precision.BF16, Precision.BF16_MIXED):
            if can_run_cueq_triattn_sm100f(dist_manager.device, torch.bfloat16, 8, dim_hidden, True):
                sm100f_per_shard_token_multiple = 8

    grid_group_sizes: OrderedDict[str, int | tuple[int, ...]] = OrderedDict(
        [("dp", size_dp), ("cp", (size_cp_axis, size_cp_axis))]
    )
    DistributedManager.create_grid_group(grid_group_sizes)

    DistributedManager.create_group(
        "world_cpu",
        dist_manager.group_ranks["world"],
        backend="gloo",
        use_local_synchronization=True,
        timeout=timeout_by_device["cpu"],
    )
    DistributedManager.create_group(
        "cp_cpu",
        dist_manager.group_ranks["cp"],
        backend="gloo",
        use_local_synchronization=True,
        timeout=timeout_by_device["cpu"],
    )
    device_mesh_cpu = map_subgroup_mesh_to_cpu(dist_manager)

    # --- Data loading ---
    torch.set_grad_enabled(False)
    if seed is not None:
        seed_everything(seed)

    out_dir = Path(out_dir).expanduser()
    mol_dir = Path(mol_dir).expanduser()
    data = Path(data).expanduser()

    if dist_manager.group_rank["world"] == 0:
        processed: Optional[BoltzProcessedInput] = None
        try:
            out_dir.mkdir(parents=True, exist_ok=True)

            if input_format == "config_files":
                results_dir = out_dir / f"boltz_results_{data.stem}"
                results_dir.mkdir(parents=True, exist_ok=True)

                data_files = check_inputs(data)
                process_inputs(
                    data=data_files,
                    out_dir=results_dir,
                    ccd_path=None,
                    mol_dir=mol_dir,
                    use_msa_server=False,
                    msa_server_url="",
                    msa_pairing_strategy="greedy",
                    boltz2=True,
                    max_msa_seqs=max_msa_seqs,
                )

                manifest = Manifest.load(results_dir / "processed" / "manifest.json")
                filtered_manifest = filter_inputs_structure(
                    manifest=manifest,
                    outdir=results_dir,
                    override=override,
                )
                processed_dir = results_dir / "processed"
                processed = BoltzProcessedInput(
                    manifest=filtered_manifest,
                    targets_dir=processed_dir / "structures",
                    msa_dir=processed_dir / "msa",
                    constraints_dir=(
                        (processed_dir / "constraints") if (processed_dir / "constraints").exists() else None
                    ),
                    template_dir=((processed_dir / "templates") if (processed_dir / "templates").exists() else None),
                    extra_mols_dir=((processed_dir / "mols") if (processed_dir / "mols").exists() else None),
                )
            elif input_format == "preprocessed":
                manifest = Manifest.load(data / "manifest.json")
                filtered_manifest = filter_inputs_structure(
                    manifest=manifest,
                    outdir=out_dir,
                    override=override,
                )
                processed = BoltzProcessedInput(
                    manifest=filtered_manifest,
                    targets_dir=data / "structures",
                    msa_dir=data / "msa",
                    constraints_dir=(data / "constraints") if (data / "constraints").exists() else None,
                    template_dir=(data / "templates") if (data / "templates").exists() else None,
                    extra_mols_dir=(data / "extra_mols") if (data / "extra_mols").exists() else None,
                )
            else:
                raise ValueError(
                    f"Unsupported input_format={input_format!r}; expected 'preprocessed' or 'config_files'"
                )
        finally:
            torch.distributed.broadcast_object_list([processed], src=0, group=dist_manager.group["world_cpu"])
    else:
        processed_recv: list[Optional[BoltzProcessedInput]] = [None]
        torch.distributed.broadcast_object_list(processed_recv, src=0, group=dist_manager.group["world_cpu"])
        processed = processed_recv[0]
        if processed is None:
            raise RuntimeError("Rank 0 failed during input processing; see rank 0 logs for the root cause.")

    out_dir = out_dir / f"boltz_results_{data.stem}"
    if dist_manager.group_rank["world"] == 0:
        out_dir.mkdir(parents=True, exist_ok=True)

    data_module = Boltz2InferenceDataModuleDTensor(
        manifest=processed.manifest,
        target_dir=processed.targets_dir,
        msa_dir=processed.msa_dir,
        mol_dir=mol_dir,
        num_workers=0,
        device_mesh=dist_manager.device_mesh_subgroups,
        device_mesh_cpu=device_mesh_cpu,
        constraints_dir=None,
        template_dir=processed.template_dir,
        extra_mols_dir=processed.extra_mols_dir,
        max_msa_seqs=max_msa_seqs,
        msa_pad_to_max_seqs=msa_pad_to_max_seqs,
        max_data_retries=max_data_retries,
        pair_mask_mode=pair_mask_mode,
        atoms_per_window_queries=atoms_per_window_queries,
        atoms_per_window_keys=atoms_per_window_keys,
        local_batch_size=local_batch_size,
        num_ensembles=num_ensembles,
        per_shard_token_multiple=sm100f_per_shard_token_multiple,
    )

    # --- Model loading ---
    predict_args = {
        "recycling_steps": recycling_steps,
        "sampling_steps": sampling_steps,
        "diffusion_samples": diffusion_samples,
        "max_parallel_samples": max_parallel_samples,
        "write_confidence_summary": confidence_prediction,
        "write_full_pae": write_full_pae,
    }

    diffusion_params = Boltz2DiffusionParams()
    diffusion_params.step_scale = step_scale

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=UserWarning)
        model_serial: Boltz2Serial = Boltz2Serial.load_from_checkpoint(
            str(checkpoint),
            strict=True,
            predict_args=predict_args,
            map_location=dist_manager.device,
            diffusion_process_args=asdict(diffusion_params),
            pairformer_args=pairformer_args,
            msa_args=msa_args,
            confidence_prediction=confidence_prediction,
            ema=False,
        )
    model_serial.eval()
    model_distributed = Boltz2Distributed(model_serial, dist_manager).eval()
    model_distributed.apply(SetTriAttnBackend(triattn_backend))
    model_distributed.apply(SetAttnPairBiasBackend(sdpa_with_bias_backend))
    model_distributed.apply(SetAttnPairBiasShardwiseBackend(sdpa_with_bias_shardwise_backend))

    # --- Callbacks ---
    callbacks_lst: list[Callback] = []
    if dist_manager.group_rank["cp"] == 0:
        name_subdir = f"predictions_dp{dist_manager.group_rank['dp']}_cp{dist_manager.group_rank['cp']}"
        pred_writer = BoltzWriter(
            data_dir=processed.targets_dir,
            output_dir=out_dir / name_subdir,
            output_format=output_format,
            boltz2=True,
        )
        callbacks_lst.append(pred_writer)
    if cuda_memory_profile:
        from boltz.workflow.utils import CUDAMemoryProfile

        mem_snapshot_path = out_dir / f"cuda_memory_profile_rank{dist_manager.group_rank['world']}.pickle"
        callbacks_lst.append(CUDAMemoryProfile(output_path=mem_snapshot_path, max_entries=300000))
    if extra_callbacks is not None:
        callbacks_lst.extend(extra_callbacks)

    callbacks: list[Callback] | None = callbacks_lst if callbacks_lst else None

    # --- Trainer ---
    strategy = SingleDeviceStrategy(device=dist_manager.device)
    devices = size_dp if size_cp == 1 else "auto"
    precision_lightning = PRECISION_TO_LIGHTNING.get(precision, "bf16-mixed")

    plugins = None
    if precision == Precision.BF16:
        plugins = [HalfPrecisionAllowFrozen(precision_lightning)]
        precision_lightning = None

    trainer = Trainer(
        default_root_dir=out_dir,
        strategy=strategy,
        callbacks=callbacks,
        accelerator=accelerator,
        devices=devices,
        precision=precision_lightning,
        plugins=plugins,
    )

    if dist_manager.group_rank["world"] == 0:
        print(f"Boltz-2 distributed inference: precision={precision}, dp={size_dp}, cp={size_cp}")

    with setup_tf32_env(precision):
        trainer.predict(
            model_distributed,
            datamodule=data_module,
            return_predictions=False,
        )
    DistributedManager.cleanup()
