# Fold-CP: A Context Parallelism Framework for Biomolecular Modeling

Context parallelism (CP) for distributed inference and training for
biomolecular folding models across multiple GPUs using a 2D CP mesh combined
with data parallelism, demonstrated with the Boltz model.

⚠️ **Note**<br>
This repository demonstrates a proof-of-concept implementation of Fold-CP with Boltz-2. <br>
Learn more about Fold-CP here: https://research.nvidia.com/labs/dbr/assets/data/manuscripts/fold_cp.pdf

For an introduction to the Boltz family of biomolecular interaction models,
see the [public Boltz repository](https://github.com/jwohlwend/boltz).

## Copyright and License Compliance

- The context parallel code is licensed under the terms and conditions as written in [the license file](../../../licenses/LICENSE)

- The original Boltz code is licensed under their respective MIT license (See the [third-party-attr.txt](../../../licenses/third-party-attr.txt))

- This project will download and install additional third-party open source software projects. Review the license terms of these open source projects before use

## Key Capabilities

- **Distributed inference** with DTensor context parallelism
- **Distributed training** with DTensor context parallelism
- Combined data parallelism (DP) and context parallelism (CP)
- Multiple attention kernel backends: cuEquivariance, trifast, FlexAttention
- Support for BF16, BF16-mixed, TF32, and FP32 precision modes

## Requirements

- Python 3.10+
- PyTorch 2.9+ with CUDA support
- Multiple NVIDIA GPUs (CP requires at least 4 GPUs; CP size must be a
  perfect square)
- `torchrun` or SLURM `srun` for multi-process launching

## Distributed Inference

Distributed inference uses `src/boltz/distributed/main.py predict` to run
structure prediction with DTensor context parallelism.

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

For full documentation of all options, the inference pipeline stages, and
differences from serial prediction, see the
[Distributed Inference Guide](docs/boltz2_cp_prediction.md).

## Distributed Training

Distributed training uses `src/boltz/distributed/train.py` with a YAML
config file to run training with DTensor context parallelism.

```bash
torchrun \
  --nnodes 1 \
  --nproc_per_node 8 \
  src/boltz/distributed/train.py \
  scripts/train/configs/structurev2_small_cp.yaml \
  parallel_size.size_dp=2 \
  parallel_size.size_cp=4 \
  output=<output_dir>
```

For full documentation of the configuration hierarchy, CP-specific settings,
CLI overrides, and differences from serial training, see the
[Distributed Training Guide](docs/boltz2_cp_training.md).

## Contributing

This project is currently not accepting contributions.

## Security

See [SECURITY.md](SECURITY.md) for vulnerability reporting instructions.

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for
details.
