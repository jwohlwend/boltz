# Boltz-2 checkpoint defaults: bounded validation

Addresses the missing attention-version and optional steering settings described
in issues #679 and #680. Those discussions also mention a community-fork fix;
this patch and its tests were implemented against upstream `b1ebfc46`.

## Automated checks

```bash
CUDA_VISIBLE_DEVICES='' OMP_NUM_THREADS=1 PYTHONPATH=src \
  python -m pytest tests/model/test_boltz2_checkpoint_defaults.py tests/model/layers -q
```

16 passed. The 14 focused cases cover dict, immutable mapping and read-only
OmegaConf inputs; explicit attention versions; a strict Lightning checkpoint
round trip; omitted and partial steering settings; and each enabled guidance
branch with a test potential. The same focused tests against unchanged upstream
gave 12 failures and 2 passing explicit-version controls. Test-sized real model
modules are used; the sampler-only cases replace the score network with a small
deterministic module.

## Released weights on H200

After preprocessing a small input with the normal CLI, run:

```bash
PYTHONPATH=src python examples/checkpoint_defaults_smoke.py \
  --cache /path/to/boltz-cache \
  --processed /path/to/boltz_results_input/processed \
  --output checkpoint-defaults-h200.json
```

The adjacent JSON records the checkpoint hash, settings, tensor comparisons and
process-local PyTorch memory peaks. The actual run used physical GPU 0, NVIDIA
H200, PyTorch 2.11.0+cu128 and Lightning 2.5.0, with 117 protein tokens, one
recycle and 20 sampling steps. This cached input's mismatching MSA was replaced
by the featurizer's standard query-only fallback in both compared calls.

The original attention default reproduces the unexpected `k_in` TypeError.
The patched model loads released weights with `strict=True`, without overriding
`pairformer_args` or `steering_args`. All 20 tensor outputs are finite and
bit-identical between omitted steering settings and explicit disabled guidance
after resetting the random seed. The input configuration remains unchanged.

The loader deliberately supplies `Boltz2DiffusionParams`, as `boltz predict`
does: the released checkpoint's training metadata includes a removed
`mse_rotational_alignment` sampler option. This check does not claim that passing
every historical checkpoint hyperparameter directly into the constructor works.
It covers the confidence checkpoint only, a shortened diffusion schedule, and a
single input. It is a correctness smoke test, not a quality or speed benchmark.

Targeted Ruff E4/E7/E9/F/I and `git diff --check` pass. The repository's entire
legacy inference/kernel test suite was not run.
