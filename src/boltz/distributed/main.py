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

"""CLI entrypoint for distributed Boltz-2 inference with DTensor context parallelism.

Thin click wrapper that resolves checkpoint / CCD-molecule paths and forwards
every option to :func:`boltz.distributed.predict.run_predict`.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import click

from boltz.data import const
from boltz.distributed.data.types import PairMaskMode
from boltz.distributed.model.modules.utils import Precision, SDPAWithBiasBackend, TriAttnBackend
from boltz.main import download_boltz2, get_cache_path


@click.group()
def cli() -> None:
    """Boltz-2 distributed inference."""
    return


@cli.command()
@click.argument("data", type=click.Path(exists=True))
@click.option(
    "--out_dir",
    type=click.Path(exists=False),
    help="The path where to save the predictions.",
    default="./",
)
@click.option(
    "--cache",
    type=click.Path(exists=False),
    help=("The directory where to download the data and model. " "Default is ~/.boltz, or $BOLTZ_CACHE if set."),
    default=get_cache_path,
)
@click.option(
    "--checkpoint",
    type=click.Path(exists=True),
    help="Path to a Boltz-2 model checkpoint. Downloaded automatically if not provided.",
    default=None,
)
@click.option(
    "--mol_dir",
    type=click.Path(exists=True),
    help="Directory containing per-residue CCD molecule pickle files. " "Resolved from --cache if not provided.",
    default=None,
)
@click.option(
    "--size_dp",
    type=int,
    help="Number of data-parallel ranks. Default is 1.",
    default=1,
)
@click.option(
    "--size_cp",
    type=int,
    help="Total number of context-parallel ranks (must be a perfect square). Default is 1.",
    default=1,
)
@click.option(
    "--accelerator",
    type=click.Choice(["gpu", "cpu"]),
    help="The accelerator to use for prediction. Default is gpu.",
    default="gpu",
)
@click.option(
    "--recycling_steps",
    type=int,
    help="The number of recycling steps. Default is 3.",
    default=3,
)
@click.option(
    "--sampling_steps",
    type=int,
    help="The number of diffusion sampling steps. Default is 200.",
    default=200,
)
@click.option(
    "--diffusion_samples",
    type=int,
    help="The number of independent diffusion samples per input. Default is 1.",
    default=1,
)
@click.option(
    "--max_parallel_samples",
    type=int,
    help="Max diffusion samples to run in parallel (None = all at once). Default is None.",
    default=None,
)
@click.option(
    "--step_scale",
    type=float,
    help=(
        "Step scale for the diffusion schedule. Lower values increase diversity "
        "among samples (recommended between 1 and 2). Default is 1.5."
    ),
    default=1.5,
)
@click.option(
    "--output_format",
    type=click.Choice(["pdb", "mmcif"]),
    help="The output structure format. Default is mmcif.",
    default="mmcif",
)
@click.option(
    "--seed",
    type=int,
    help="Random seed for reproducibility. Default is None.",
    default=None,
)
@click.option(
    "--max_msa_seqs",
    type=int,
    help=f"Maximum number of MSA sequences. Default is {const.max_msa_seqs}.",
    default=const.max_msa_seqs,
)
@click.option(
    "--msa_pad_to_max_seqs",
    is_flag=True,
    help="Whether to pad MSA to max_msa_seqs. Default is False.",
)
@click.option(
    "--input_format",
    type=click.Choice(["preprocessed", "config_files"], case_sensitive=False),
    help="Data format for input. If 'preprocessed', expects a folder with "
    "manifest.json, msa/ folder and structures/ folder with preprocessed data. "
    "If 'config_files', expects yaml/fasta files. Default is preprocessed.",
    default="preprocessed",
)
@click.option(
    "--timeout_nccl",
    type=float,
    help="NCCL timeout in minutes. Default is 30.",
    default=30,
)
@click.option(
    "--timeout_gloo",
    type=float,
    help="Gloo timeout in minutes. Default is 30.",
    default=30,
)
@click.option(
    "--precision",
    type=click.Choice([Precision.BF16.value, Precision.BF16_MIXED.value, Precision.TF32.value, Precision.FP32.value]),
    help="Model precision mode. Default is BF16_MIXED.",
    default=Precision.BF16_MIXED.value,
)
@click.option(
    "--atoms_per_window_queries_keys",
    nargs=2,
    type=int,
    help="(queries, keys) window sizes for atom attention batching. Default is 32 128.",
    default=(32, 128),
)
@click.option(
    "--pair_mask_mode",
    type=click.Choice(
        [PairMaskMode.NONE.value, PairMaskMode.GLOBAL_ATOM_ATTENTION.value, PairMaskMode.SEQUENCE_LOCAL_ATTENTION.value]
    ),
    help="Pair mask mode. Default is None (window batching).",
    default=PairMaskMode.NONE.value,
)
@click.option(
    "--local_batch_size",
    type=int,
    help="Per-rank batch size. Default is 1.",
    default=1,
)
@click.option(
    "--num_ensembles",
    type=int,
    help="Number of ensemble members. Default is 1.",
    default=1,
)
@click.option(
    "--write_full_pae",
    is_flag=True,
    help="Whether to write full PAE matrices. Default is False.",
)
@click.option(
    "--use_templates",
    type=bool,
    help="Whether to use template features. Default is True.",
    default=True,
)
@click.option(
    "--triattn_backend",
    type=click.Choice([TriAttnBackend.CUEQ.value, TriAttnBackend.TRIFAST.value, TriAttnBackend.REFERENCE.value]),
    help="Triangle attention backend to use. Default is cueq.",
    default=TriAttnBackend.CUEQ.value,
)
@click.option(
    "--sdpa_with_bias_backend",
    type=click.Choice([SDPAWithBiasBackend.REFERENCE.value, SDPAWithBiasBackend.TORCH_FLEX_ATTN.value]),
    help="SDPA backend for ring-attention AttentionPairBias layers. Default is torch_flex_attn.",
    default=SDPAWithBiasBackend.TORCH_FLEX_ATTN.value,
)
@click.option(
    "--sdpa_with_bias_shardwise_backend",
    type=click.Choice(
        [
            SDPAWithBiasBackend.REFERENCE.value,
            SDPAWithBiasBackend.TORCH_SDPA_EFFICIENT_ATTENTION.value,
            SDPAWithBiasBackend.TORCH_FLEX_ATTN.value,
        ]
    ),
    help="SDPA backend for window-batched AttentionPairBiasShardwise layers. Default is torch_flex_attn.",
    default=SDPAWithBiasBackend.TORCH_FLEX_ATTN.value,
)
@click.option(
    "--auto_pad_tokens_for_sm100f/--no_auto_pad_tokens_for_sm100f",
    default=True,
    help="Pad token counts so each CP shard is a multiple of 8 for SM100f cuEq TriAttn. Default is True.",
)
@click.option(
    "--cuda_memory_profile",
    is_flag=True,
    default=False,
    help="Profile CUDA memory usage and dump a snapshot pickle per rank.",
)
@click.option(
    "--override",
    is_flag=True,
    default=False,
    help="Override existing predictions even if output already exists. Default is False.",
)
def predict(
    data: str,
    out_dir: str,
    cache: str,
    checkpoint: Optional[str],
    mol_dir: Optional[str],
    size_dp: int,
    size_cp: int,
    accelerator: str,
    recycling_steps: int,
    sampling_steps: int,
    diffusion_samples: int,
    max_parallel_samples: Optional[int],
    step_scale: float,
    output_format: str,
    seed: Optional[int],
    max_msa_seqs: int,
    msa_pad_to_max_seqs: bool,
    input_format: str,
    timeout_nccl: float,
    timeout_gloo: float,
    precision: str,
    atoms_per_window_queries_keys: tuple[int, int],
    pair_mask_mode: str,
    local_batch_size: int,
    num_ensembles: int,
    write_full_pae: bool,
    use_templates: bool,
    triattn_backend: str,
    sdpa_with_bias_backend: str,
    sdpa_with_bias_shardwise_backend: str,
    auto_pad_tokens_for_sm100f: bool,
    cuda_memory_profile: bool,
    override: bool,
) -> None:
    """Run distributed Boltz-2 structure prediction.

    DATA is the path to the input data directory.
    """
    cache_path = Path(cache).expanduser()
    cache_path.mkdir(parents=True, exist_ok=True)

    # Resolve checkpoint: download if not provided
    if checkpoint is None:
        download_boltz2(cache_path)
        checkpoint = str(cache_path / "boltz2_conf.ckpt")

    # Resolve mol_dir: use cache default if not provided
    if mol_dir is None:
        download_boltz2(cache_path)
        mol_dir = str(cache_path / "mols")

    from boltz.distributed.predict import run_predict

    run_predict(
        data=data,
        out_dir=out_dir,
        mol_dir=mol_dir,
        checkpoint=checkpoint,
        size_dp=size_dp,
        size_cp=size_cp,
        accelerator=accelerator,
        recycling_steps=recycling_steps,
        sampling_steps=sampling_steps,
        diffusion_samples=diffusion_samples,
        max_parallel_samples=max_parallel_samples,
        step_scale=step_scale,
        output_format=output_format,
        seed=seed,
        max_msa_seqs=max_msa_seqs,
        msa_pad_to_max_seqs=msa_pad_to_max_seqs,
        input_format=input_format,
        timeout_nccl=timeout_nccl,
        timeout_gloo=timeout_gloo,
        precision=Precision(precision),
        atoms_per_window_queries_keys=atoms_per_window_queries_keys,
        pair_mask_mode=PairMaskMode(pair_mask_mode),
        local_batch_size=local_batch_size,
        num_ensembles=num_ensembles,
        write_full_pae=write_full_pae,
        use_templates=use_templates,
        triattn_backend=TriAttnBackend(triattn_backend),
        sdpa_with_bias_backend=SDPAWithBiasBackend(sdpa_with_bias_backend),
        sdpa_with_bias_shardwise_backend=SDPAWithBiasBackend(sdpa_with_bias_shardwise_backend),
        auto_pad_tokens_for_sm100f=auto_pad_tokens_for_sm100f,
        cuda_memory_profile=cuda_memory_profile,
        override=override,
    )


if __name__ == "__main__":
    cli()
