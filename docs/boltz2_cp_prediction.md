# Boltz CP Distributed Prediction with DTensor Context Parallelism

> **Note:** The current implementation supports Boltz-2 only.

This document describes how to run distributed structure prediction
using `src/boltz/distributed/main.py`, which provides a Click CLI that
delegates to `src/boltz/distributed/predict.py::run_predict` for DTensor-based
context parallelism (CP) combined with data parallelism (DP).

## Entrypoint

The distributed prediction CLI is:

```
src/boltz/distributed/main.py predict <DATA> [options]
```

`DATA` is the path to the input data directory. The CLI resolves the model
checkpoint and CCD molecule directory automatically (downloading them to
`~/.boltz` if not provided), then forwards all options to `run_predict`.

### Launching with `torchrun` or `srun`

The script is designed to be launched with **either** `torchrun` (for
single-node or multi-node runs outside SLURM) **or** `srun` (for SLURM
clusters). The `DistributedManager` (`src/boltz/distributed/manager.py`)
auto-detects the launch method:

- **`torchrun`** — detected via `RANK`, `WORLD_SIZE`, `LOCAL_RANK`,
  `MASTER_ADDR`, `MASTER_PORT` environment variables.
- **`srun` (SLURM)** — detected via `SLURM_PROCID`, `SLURM_NPROCS`,
  `SLURM_LOCALID`, `SLURM_LAUNCH_NODE_IPADDR` environment variables.

The detection order can be forced by setting
`BOLTZ_DISTRIBUTED_INIT_METHOD=ENV` or `BOLTZ_DISTRIBUTED_INIT_METHOD=SLURM`.

Example with `torchrun` (single node, 4 GPUs, dp=1, cp=4):

```bash
torchrun \
  --nnodes 1 \
  --nproc_per_node 4 \
  src/boltz/distributed/main.py predict \
  /path/to/preprocessed_data \
  --out_dir ./predictions \
  --size_dp 1 \
  --size_cp 4 \
  --recycling_steps 3 \
  --sampling_steps 200 \
  --diffusion_samples 5
```

Example with `srun` (multi-node SLURM, 8 GPUs total, dp=2, cp=4):

```bash
srun --ntasks-per-node=4 --nodes=2 \
  python src/boltz/distributed/main.py predict \
  /path/to/preprocessed_data \
  --out_dir ./predictions \
  --size_dp 2 \
  --size_cp 4 \
  --checkpoint /path/to/boltz2_conf.ckpt \
  --mol_dir /path/to/mols
```

The constraint `size_dp * size_cp == world_size` must hold, and `size_cp`
must be a perfect square (the CP mesh is 2D).

---

## Input Data

The distributed inference pipeline currently supports only **preprocessed**
input data (`--input_format preprocessed`). The data directory must contain:

- `manifest.json` — describes the samples to predict.
- `structures/` — preprocessed structure files.
- `msa/` — MSA files for each target.
- `templates/` (optional) — template structure files.
- `extra_mols/` (optional) — additional molecule definitions.

Rank 0 loads the manifest and broadcasts it to all other ranks via
`torch.distributed.broadcast_object_list` over a CPU (Gloo) process group.

---

## CLI Options

All options are provided as Click flags on the `predict` subcommand.

### Required Arguments

| Argument | Description                              |
| -------- | ---------------------------------------- |
| `DATA`   | Path to the preprocessed data directory. |

### Common Options

| Option         | Type | Default       | Description                                                               |
| -------------- | ---- | ------------- | ------------------------------------------------------------------------- |
| `--out_dir`    | path | `./`          | Output directory for predictions.                                         |
| `--cache`      | path | `~/.boltz`    | Download cache for checkpoint and CCD molecules. Respects `$BOLTZ_CACHE`. |
| `--checkpoint` | path | auto-download | Path to a Boltz model checkpoint.                                         |
| `--mol_dir`    | path | auto-download | Directory containing per-residue CCD molecule pickle files.               |

### Parallelism

| Option      | Type | Default | Description                                                                |
| ----------- | ---- | ------- | -------------------------------------------------------------------------- |
| `--size_dp` | int  | `1`     | Number of data-parallel ranks.                                             |
| `--size_cp` | int  | `1`     | Total context-parallel ranks (must be a perfect square: 1, 4, 9, 16, ...). |

The product `size_dp * size_cp` must equal the total world size
(`nproc_per_node * nnodes`).

### Diffusion Sampling

| Option                   | Type  | Default | Description                                                                              |
| ------------------------ | ----- | ------- | ---------------------------------------------------------------------------------------- |
| `--recycling_steps`      | int   | `3`     | Number of trunk recycling iterations.                                                    |
| `--sampling_steps`       | int   | `200`   | Number of diffusion denoising steps.                                                     |
| `--diffusion_samples`    | int   | `1`     | Number of independent diffusion samples per input.                                       |
| `--max_parallel_samples` | int   | `None`  | Max samples to run in parallel (`None` = all at once).                                   |
| `--step_scale`           | float | `1.5`   | Diffusion schedule step scale. Lower values increase sample diversity (recommended 1–2). |

### Model and Precision

| Option          | Type   | Default      | Description                                            |
| --------------- | ------ | ------------ | ------------------------------------------------------ |
| `--precision`   | enum   | `BF16_MIXED` | Model precision: `BF16`, `BF16_MIXED`, `TF32`, `FP32`. |
| `--accelerator` | choice | `gpu`        | Device accelerator: `gpu` or `cpu`.                    |
| `--seed`        | int    | `None`       | Random seed for reproducibility.                       |

### Attention Kernel Backends

| Option                               | Values                                                           | Default           | Description                                                       |
| ------------------------------------ | ---------------------------------------------------------------- | ----------------- | ----------------------------------------------------------------- |
| `--triattn_backend`                  | `cueq`, `trifast`, `reference`                                   | `cueq`            | Triangle attention kernel. `cueq` requires CUDA + cuequivariance. |
| `--sdpa_with_bias_backend`           | `reference`, `torch_flex_attn`                                   | `torch_flex_attn` | SDPA backend for ring-attention layers.                           |
| `--sdpa_with_bias_shardwise_backend` | `reference`, `torch_sdpa_efficient_attention`, `torch_flex_attn` | `torch_flex_attn` | SDPA backend for window-batched attention layers.                 |

### Data Processing

| Option                  | Type   | Default        | Description                                                                                                     |
| ----------------------- | ------ | -------------- | --------------------------------------------------------------------------------------------------------------- |
| `--input_format`        | choice | `preprocessed` | Input data format. Only `preprocessed` is currently supported.                                                  |
| `--max_msa_seqs`        | int    | `4096`         | Maximum number of MSA sequences.                                                                                |
| `--msa_pad_to_max_seqs` | flag   | `False`        | Pad MSA to `max_msa_seqs`.                                                                                      |
| `--use_templates`       | bool   | `True`         | Reserved for future use. Template weights are loaded but the distributed TemplateModule is not yet implemented. |

### Window Batching

| Option                            | Type    | Default  | Description                                                                                   |
| --------------------------------- | ------- | -------- | --------------------------------------------------------------------------------------------- |
| `--atoms_per_window_queries_keys` | int int | `32 128` | (queries, keys) window sizes for atom attention batching.                                     |
| `--pair_mask_mode`                | choice  | `None`   | Pair mask mode: `None` (window batching), `GlobalAtomAttention`, or `SequenceLocalAttention`. |

### Output

| Option               | Type   | Default | Description                                           |
| -------------------- | ------ | ------- | ----------------------------------------------------- |
| `--output_format`    | choice | `mmcif` | Output structure format: `pdb` or `mmcif`.            |
| `--write_full_pae`   | flag   | `False` | Write full PAE matrices (requires confidence module). |
| `--local_batch_size` | int    | `1`     | Per-rank batch size.                                  |
| `--num_ensembles`    | int    | `1`     | Number of ensemble members for structure prediction.  |

### Timeouts and Profiling

| Option                  | Type  | Default | Description                                               |
| ----------------------- | ----- | ------- | --------------------------------------------------------- |
| `--timeout_nccl`        | float | `30`    | NCCL timeout in minutes (for CUDA).                       |
| `--timeout_gloo`        | float | `30`    | Gloo timeout in minutes (for CPU).                        |
| `--cuda_memory_profile` | flag  | `False` | Dump a CUDA memory snapshot pickle per rank to `out_dir`. |

---

## Inference Pipeline Stages

The `run_predict` function in `src/boltz/distributed/predict.py` executes
the following stages:

### 1. Distributed Setup

- Initializes `DistributedManager` with the appropriate device type and timeout.
- Validates `size_dp * size_cp == world_size` and that `size_cp` is a perfect
  square.
- Creates a 2D CP grid: `OrderedDict([("dp", size_dp), ("cp", (sqrt_cp, sqrt_cp))])`.
- Creates CPU-backed Gloo process groups (`world_cpu`, `cp_cpu`) for data
  broadcast and metadata exchange.

### 2. Data Broadcast and Loading

- Rank 0 loads the manifest and constructs `BoltzProcessedInput`.
- The processed input is broadcast to all ranks via
  `broadcast_object_list` over the `world_cpu` group.
- A `Boltz2InferenceDataModuleDTensor` data module is created with the device
  mesh, which uses `PredictionDatasetCPWithDTensorV2` internally. Each sample
  is featurized, tokenized, and distributed as DTensors across CP ranks.

### 3. Model Loading

- The serial Boltz model is loaded from the checkpoint using
  `Boltz2Serial.load_from_checkpoint` with `strict=True`.
- Checkpoint hyperparameters (e.g., `pairformer_args.v2`, `msa_args.use_paired_feature`)
  are read from the checkpoint and merged to ensure the correct model
  architecture is instantiated.
- The serial model is wrapped in the `Boltz2Distributed` wrapper, which
  replaces submodules with DTensor-aware distributed counterparts.
- Attention kernel backends are configured via `model.apply(SetTriAttnBackend(...))`,
  `model.apply(SetAttnPairBiasBackend(...))`, and
  `model.apply(SetAttnPairBiasShardwiseBackend(...))`.

### 4. Prediction

- A Lightning `Trainer` with `SingleDeviceStrategy` runs `trainer.predict()`.
- Only CP rank 0 within each DP group writes output files, via `BoltzWriter`.
- Each DP rank writes to its own subdirectory:
  `<out_dir>/boltz_results_<data_stem>/predictions_dp{dp_rank}_cp0/`.
- Prediction runs inside a `setup_tf32_env` context when TF32 precision
  is selected.

---

## CP-Specific Settings (not in serial `predict`)

The following settings exist in the distributed `predict` CLI but have no
counterpart in the serial `predict` command:

### Parallelism topology

| Key         | Type | Description                                             |
| ----------- | ---- | ------------------------------------------------------- |
| `--size_dp` | int  | Data-parallel group size.                               |
| `--size_cp` | int  | Context-parallel group size (must be a perfect square). |

### Precision

A top-level `--precision` enum (`BF16`, `BF16_MIXED`, `TF32`, `FP32`) that
replaces Lightning's `trainer.precision`. For `BF16` mode, a custom
`HalfPrecisionAllowFrozen` plugin is used to handle Boltz's frozen
dataclasses during input conversion.

### Attention kernel backends

| Option                               | Values                                                           | Description                                                       |
| ------------------------------------ | ---------------------------------------------------------------- | ----------------------------------------------------------------- |
| `--triattn_backend`                  | `reference`, `cueq`, `trifast`                                   | Triangular attention kernel. `cueq` does not support FP32 or CPU. |
| `--sdpa_with_bias_backend`           | `reference`, `torch_flex_attn`                                   | SDPA backend for ring-attention layers.                           |
| `--sdpa_with_bias_shardwise_backend` | `reference`, `torch_sdpa_efficient_attention`, `torch_flex_attn` | SDPA backend for window-batched attention layers.                 |

### CUDA memory profiling

When `--cuda_memory_profile` is set, each rank writes a pickle file to
`<out_dir>/cuda_memory_profile_rank<world_rank>.pickle`.

### Timeouts

| Key              | Type  | Default | Description                             |
| ---------------- | ----- | ------- | --------------------------------------- |
| `--timeout_nccl` | float | `30`    | NCCL timeout in minutes (CUDA).         |
| `--timeout_gloo` | float | `30`    | Gloo timeout in minutes (CPU/metadata). |

---

## Differences from Serial Prediction

| Aspect                | Serial (`src/boltz/main.py predict`)        | CP (`src/boltz/distributed/main.py predict`)                                          |
| --------------------- | ------------------------------------------- | ------------------------------------------------------------------------------------- |
| Multi-GPU strategy    | Lightning DDP (`DDPStrategy`)               | `SingleDeviceStrategy` + DTensor CP mesh                                              |
| Device management     | Lightning (`--devices`, `--num_nodes`)      | `DistributedManager` via `--size_dp`, `--size_cp`                                     |
| Launch method         | `python src/boltz/main.py predict`          | `torchrun` or `srun`                                                                  |
| Input formats         | `config_files` (YAML/FASTA), `preprocessed` | `preprocessed` only                                                                   |
| `num_workers`         | Configurable                                | Fixed at `0` (DTensor CP requires main-process collation)                             |
| Precision             | Lightning `--precision` string              | Top-level `--precision` enum                                                          |
| Attention backends    | Not configurable                            | `--triattn_backend`, `--sdpa_with_bias_backend`, `--sdpa_with_bias_shardwise_backend` |
| CUDA memory profiling | Not available                               | `--cuda_memory_profile` flag                                                          |
| Confidence prediction | Supported                                   | Not yet supported (`write_confidence_summary=False`)                                  |
| Steering potentials   | Supported                                   | Not yet supported                                                                     |
| Affinity prediction   | Supported                                   | Not yet supported                                                                     |
| Template features     | Supported                                   | Weights loaded but distributed TemplateModule not yet implemented                     |
| Constraint features   | Supported                                   | Not yet supported                                                                     |
| Checkpoint loading    | Lightning default                           | Reads checkpoint hparams, merges v2 flags, loads with `strict=True`                   |
| Output writing        | All ranks write                             | Only CP rank 0 per DP group writes output                                             |

---
