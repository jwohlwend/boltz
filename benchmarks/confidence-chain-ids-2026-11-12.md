# Confidence chain IDs: compatibility validation

Issue #173 asks for the input chain IDs alongside per-chain confidence scores.
The patch adds `chain_ids`, an index-to-input-ID mapping, and retains all legacy
score fields and numeric keys. It uses the writer's existing chain remapping
after invalid chains are removed.

## Automated checks

```bash
CUDA_VISIBLE_DEVICES='' OMP_NUM_THREADS=1 PYTHONPATH=src \
  python -m pytest tests/data/test_confidence_chain_ids.py tests/model/layers -q
```

12 focused cases and 2 neighboring layer tests pass. Against unchanged upstream
`b1ebfc46`, 11 focused cases fail because the mapping is absent and the
no-confidence control passes. The tests exercise actual PDB/mmCIF writing in
Boltz-1 and Boltz-2, masked-chain removal, two differently ranked predictions,
single chains, long original IDs, repeated names and disabled confidence.
Every legacy confidence field and PLDDT array is checked independently of the
new mapping. No model or tensor computations are changed.

## Replay of a saved prediction

```bash
CUDA_VISIBLE_DEVICES='' PYTHONPATH=src python examples/confidence_chain_ids_replay.py \
  --run /path/to/boltz_results_input --record-id INPUT_ID \
  --original-writer /path/to/unpatched/src/boltz/data/write/writer.py \
  --output /path/to/replay-output
```

An existing L3001 prediction (two chains, A/B) was replayed on CPU through the
unmodified and patched writer. All legacy confidence fields equal the original
saved summary. The mmCIF bytes, PLDDT payload and saved pre-affinity structure
payload are identical; the only confidence-summary addition is
`chain_ids: {"0": "A", "1": "B"}`. The adjacent JSON records the source summary
hash and the checked files. NPZ member bytes are compared without unpickling;
ZIP container timestamps are excluded from the comparison.

Targeted Ruff E4/E7/E9/F/I and `git diff --check` pass. This is output
compatibility validation using synthetic fixtures and one saved real result,
not new inference, a GPU benchmark or evidence about model quality. The full
legacy inference/kernel test suite was not run.
