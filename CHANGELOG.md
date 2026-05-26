# Boltz Context Parallelism Changelog

## 0.1.0 (initial release)

### New Features

- Distributed inference with DTensor context parallelism (`src/boltz/distributed/predict.py`)
- Distributed training with DTensor context parallelism (`src/boltz/distributed/train.py`)
- 2D CP mesh support (Shard x Shard context parallelism)
- Data parallelism combined with context parallelism (DP x CP)
- Multiple attention kernel backends: cuEquivariance, trifast, FlexAttention
- CUDA memory profiling support
- Preprocessed input format for distributed inference
