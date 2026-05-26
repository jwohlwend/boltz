# Boltz CP Distributed Training with DTensor Context Parallelism

> **Note:** The current implementation supports Boltz-2 only.

This document describes how to run distributed training using
`src/boltz/distributed/train.py`, which implements DTensor-based context
parallelism (CP) combined with data parallelism (DP).

## Entrypoint

The distributed training entrypoint is:

```
src/boltz/distributed/train.py <config.yaml> [override1=value1] [override2=value2] ...
```

It accepts a YAML config file as the first positional argument, followed by
zero or more OmegaConf-style dotlist overrides. The config is loaded via
Hydra (respecting `defaults:` chains), struct mode is disabled to allow
adding new keys, and CLI overrides are merged on top.

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

Example with `torchrun` (single node, 8 GPUs):

```bash
torchrun \
  --nnodes 1 \
  --nproc_per_node 8 \
  src/boltz/distributed/train.py \
  scripts/train/configs/structurev2_small_cp.yaml \
  parallel_size.size_dp=2 \
  parallel_size.size_cp=4 \
  output=<output_dir> \
  ...
```

Example with `srun` (multi-node SLURM):

```bash
srun --ntasks-per-node=8 --nodes=2 \
  python src/boltz/distributed/train.py \
  scripts/train/configs/structurev2_cp.yaml \
  parallel_size.size_dp=4 \
  parallel_size.size_cp=4 \
  output=<output_dir> \
  ...
```

The constraint `size_dp * size_cp == world_size` must hold, and `size_cp`
must be a perfect square (the CP mesh is 2D).

---

## Configuration Hierarchy

### Base config: `structurev2.yaml`

`scripts/train/configs/structurev2.yaml` is the base configuration for
Boltz serial training (used by `scripts/train/train.py`). It defines the
full model architecture, data pipeline, training hyperparameters, and
defaults to single-device training:

```yaml
trainer:
  accelerator: cuda
  devices: 1
  num_nodes: 1
  precision: bf16-mixed
```

### CP overlay: `structurev2_cp.yaml`

`scripts/train/configs/structurev2_cp.yaml` extends `structurev2.yaml` with
settings required for DTensor context-parallel training:

```yaml
defaults:
  - structurev2
  - _self_

trainer:
  accelerator: gpu     # must be gpu (not cuda) — CP code manages devices
  devices: 1           # must be 1 — each Lightning Trainer manages one device
  num_nodes: 1         # must be 1 — multi-node is handled by torchrun/SLURM
  precision: null      # superseded by top-level precision setting

parallel_size:
  size_cp: 1           # context-parallel group size (must be a perfect square)
  size_dp: 1           # data-parallel group size
  timeout_nccl: 30     # NCCL timeout in minutes
  timeout_gloo: 30     # Gloo timeout in minutes

precision: BF16_MIXED  # top-level precision enum (FP32, TF32, BF16, BF16_MIXED, FP16, FP64)

triattn_backend: cueq
sdpa_with_bias_backend: torch_flex_attn
sdpa_with_bias_shardwise_backend: torch_flex_attn

data:
  num_workers: 0       # must be 0 — DTensor CP requires main-process collation
  pin_memory: false
```

Key differences from the serial config:

| Setting | Serial (`structurev2.yaml`) | CP (`structurev2_cp.yaml`) |
|---|---|---|
| `trainer.accelerator` | `cuda` | `gpu` |
| `trainer.devices` | configurable (e.g. 8) | must be `1` |
| `trainer.num_nodes` | configurable (e.g. 4) | must be `1` |
| `trainer.precision` | `bf16-mixed` | `null` (use top-level `precision`) |
| `parallel_size` | not present | `size_dp`, `size_cp`, timeouts |
| `precision` (top-level) | not present | `BF16_MIXED` enum |
| `triattn_backend` | not present | triangular attention kernel backend |
| `sdpa_with_bias_backend` | not present | ring-attention SDPA backend |
| `sdpa_with_bias_shardwise_backend` | not present | window-batched SDPA backend |
| `data.num_workers` | `2` | `0` (required) |
| `data.pin_memory` | `false` | `false` |
| `CUDAMemoryProfile` | not present | optional memory profiling |

### Small-model variants

- **`structurev2_small.yaml`** extends `structurev2.yaml` with reduced model
  depth (12 pairformer blocks, 3 MSA blocks), smaller sequence limits
  (`max_tokens: 256`, `max_atoms: 2048`), and no activation checkpointing.
  It also sets `trainer.devices: 8` and `trainer.num_nodes: 4` for multi-GPU
  serial DDP training.

- **`structurev2_small_cp.yaml`** extends `structurev2_cp.yaml` (not the
  serial small variant) with the same reduced model depth and sequence
  limits. It does **not** set `trainer.devices` or `trainer.num_nodes`
  because the CP config already fixes those to `1`. The parallel topology is
  controlled entirely via `parallel_size.size_dp` and
  `parallel_size.size_cp` CLI overrides.

---

## CP-Specific Settings (not in serial `train.py`)

The following settings exist in `DistributedTrainConfig` but have no
counterpart in the serial `TrainConfig`:

### `parallel_size`

Controls the distributed topology.

| Key | Type | Description |
|---|---|---|
| `size_dp` | int | Data-parallel group size |
| `size_cp` | int | Context-parallel group size (must be a perfect square: 1, 4, 9, 16, ...) |
| `timeout_nccl` | int | NCCL timeout in minutes (for CUDA) |
| `timeout_gloo` | int | Gloo timeout in minutes (for CPU) |

The product `size_dp * size_cp` must equal the total world size
(`nproc_per_node * nnodes`).

### `precision` (top-level)

An enum (`FP32`, `TF32`, `BF16`, `BF16_MIXED`, `FP16`, `FP64`) that
replaces `trainer.precision`. Setting `trainer.precision` directly raises an
error in the CP entrypoint.

### Attention kernel backends

| Key | Values | Description |
|---|---|---|
| `triattn_backend` | `reference`, `cueq`, `trifast`, `cueq_fwd_trifast_bwd` | Triangular attention kernel. `cueq` does not support FP32. |
| `sdpa_with_bias_backend` | `reference`, `torch_sdpa_efficient_attention`, `torch_flex_attn` | SDPA backend for ring-attention layers |
| `sdpa_with_bias_shardwise_backend` | `reference`, `torch_sdpa_efficient_attention`, `torch_flex_attn` | SDPA backend for window-batched (shardwise) attention layers |

### `OffloadActvCkptToCPU`

Enables CPU offloading of activation-checkpoint boundary tensors for
selected distributed module types. When activation checkpointing is active,
intermediate activations saved for backward are moved to CPU during the
forward pass and restored on backward, trading extra CPU-GPU transfers for
reduced GPU memory.

```yaml
# Disabled by default.
OffloadActvCkptToCPU: null

# Enable for all three supported module types:
OffloadActvCkptToCPU:
  - DiffusionTransformer
  - MSAModule
  - PairformerModule
```

CLI override:

```bash
OffloadActvCkptToCPU='[DiffusionTransformer,MSAModule,PairformerModule]'
```

Valid module type names: `DiffusionTransformer`, `MSAModule`,
`PairformerModule`. Any subset may be specified.

**Prerequisite: activation checkpointing must be enabled on every targeted
module.** The setter raises `ValueError` if a targeted module has
`activation_checkpointing=False`. The standard activation-checkpointing
overrides cover four config sections:

```bash
model.msa_args.activation_checkpointing=true
model.template_args.activation_checkpointing=true
model.pairformer_args.activation_checkpointing=true
model.score_model_args.activation_checkpointing=true
```

However, there is a fifth `DiffusionTransformer` instance nested inside the
`InputEmbedder` (via `AtomAttentionEncoder` -> `AtomTransformer` ->
`DiffusionTransformer`). This instance is controlled by `embedder_args`, not
`score_model_args`. If `DiffusionTransformer` is included in
`OffloadActvCkptToCPU`, you must also enable:

```bash
model.embedder_args.activation_checkpointing=true
```

Otherwise the setter will raise because the `InputEmbedder`'s
`DiffusionTransformer` still has `activation_checkpointing=False`.

### `CUDAMemoryProfile`

Optional CUDA memory profiling. Each rank writes a pickle file.

```yaml
CUDAMemoryProfile:
  output_path_prefix: null   # set a path prefix to enable, e.g. "profiling/mem"
  max_entries: 300000
```

### `checkpoint`

Overrides for Lightning's `ModelCheckpoint`. The CP entrypoint applies
Boltz defaults (`monitor="val/lddt"`, `save_last=True`,
`save_on_train_epoch_end=True`, `mode="max"`, `every_n_epochs=1`) but any
key can be overridden via CLI:

```bash
checkpoint.monitor=val/disto_lddt_protein_protein \
checkpoint.save_top_k=3
```

### `seed`

When set, seeds are offset by `rank + epoch + global_step` on resume to
avoid replaying identical data samples across ranks and restarts.

### `validation_only`

When `true`, runs `trainer.validate()` instead of `trainer.fit()`. Useful
for evaluating a checkpoint without training.

---

## CLI Overrides

All config keys can be overridden from the command line using OmegaConf
dotlist syntax. The config's struct mode is disabled, so new keys can also
be introduced.

### Dataset overrides via CLI

A recently added utility (`convert_datasets_dict_to_list_config` from
`src/boltz/workflow/utils.py`) enables **partial overrides of individual
dataset entries** within the `data.datasets` list directly from the CLI.
This was not possible before because OmegaConf does not natively support
partial updates to ListConfig entries.

The utility converts dict-style index keys (e.g. `data.datasets.0.key`)
into proper list merges against the base config. The overridable dataset
keys are:

`_target_`, `target_dir`, `msa_dir`, `prob`, `sampler`, `cropper`,
`template_dir`, `filters`, `split`, `symmetry_correction`, `val_group`,
`use_train_subset`, `moldir`, `override_bfactor`, `override_method`

Examples:

```bash
# Override dataset 0's data directory and disable filters
data.datasets.0.target_dir=/path/to/data \
data.datasets.0.msa_dir=/path/to/msa \
'data.datasets.0.filters=[]' \

# Remove dataset 1 entirely (null removes the entry)
data.datasets.1=null \
```

When `filters` needs to be set to an empty list, it must be quoted to
prevent shell interpretation: `'data.datasets.0.filters=[]'`.

### Common CLI overrides

```bash
# Parallelism topology
parallel_size.size_dp=<int>
parallel_size.size_cp=<int>

# Output and checkpointing
output=<output_dir>
resume=<checkpoint_path>
pretrained=<pretrained_path>
checkpoint.monitor=<metric_name>

# Training schedule
trainer.max_epochs=<int>
trainer.accumulate_grad_batches=<int>

# Logging frequency
trainer.log_every_n_steps=<int>     # Lightning: how often to flush to logger (default: 50)
model.log_loss_every_steps=<int>    # Model: how often to compute and log losses/norms (default: 50)

# Data
data.samples_per_epoch=<int>
data.overfit=<int>                  # Activate overfit mode: use first N samples, validate on train data
data.datasets.0.target_dir=<path>
data.datasets.0.split=<path_or_null>

# Precision
precision=BF16_MIXED

# WandB (all keys required when wandb section is present)
wandb.name=<run_name>
wandb.project=<project>
wandb.id=<run_id>
wandb.entity=<entity_or_null>

# Validation only
validation_only=true
```

### Logging frequency

There are two independent gates that control when training metrics appear
in the logger:

1. **`model.log_loss_every_steps`** (default: 50) — The Boltz model only
   calls `self.log()` (and computes parameter/gradient norms) every N
   global steps. Steps where this condition is not met skip all logging
   computation entirely.

2. **`trainer.log_every_n_steps`** (default: 50) — PyTorch Lightning only
   flushes `self.log()` calls to the configured logger (e.g. WandB) every
   N steps. Calls on other steps are buffered but not written.

Both conditions must be satisfied for metrics to appear. For small
datasets with few steps per epoch, set both to `1`.

### Overfit mode

Setting `data.overfit=N` activates overfit mode:

- Training samples are truncated to the first N per dataset.
- Validation uses the **training** datasets instead of the validation split.
- Sample weights are normalized for uniform sampling.

When using overfit mode, set `data.datasets.0.split=null` to prevent the
train/val split from reducing the training set. Set
`data.samples_per_epoch` to control the epoch length (divided by `size_dp`
for per-rank step count).

---

## Differences from Serial Training (`scripts/train/train.py`)

| Aspect | Serial (`scripts/train/train.py`) | CP (`src/boltz/distributed/train.py`) |
|---|---|---|
| Multi-GPU strategy | Lightning DDP (`DDPStrategy`) | `BoltzContextParallelStrategy` (DTensor) |
| Device management | Lightning (`trainer.devices`, `trainer.num_nodes`) | `DistributedManager` via `parallel_size` |
| Launch method | `python scripts/train/train.py` | `torchrun` or `srun` |
| `num_workers` | Configurable (default: 2) | Must be `0` |
| Precision | `trainer.precision` (Lightning string) | Top-level `precision` enum |
| Attention backends | Not configurable | `triattn_backend`, `sdpa_with_bias_backend`, `sdpa_with_bias_shardwise_backend` |
| CUDA memory profiling | Not available | `CUDAMemoryProfile` section |
| Confidence prediction | Supported | Not yet supported (auto-disabled with warning) |
| Steering potentials | Supported | Not yet supported (auto-disabled with warning) |
| Checkpoint strategy | Lightning default | `BoltzContextParallelStrategy` (DTensor-aware save/load) |
