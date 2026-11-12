# Training record metadata validation

Issue #686 reports mismatched cluster keys and absent MSA links in the RCSB
training preprocessor. Implemented against upstream `b1ebfc46` after reading
the issue and the producer/consumer code. The discussion also mentions a
community-fork fix; no code from that implementation was used.

## Reproduce

Install Boltz and `scripts/process/requirements.txt`, then run:

```bash
CUDA_VISIBLE_DEVICES='' OMP_NUM_THREADS=1 PYTHONPATH=src \
  python -m pytest tests/data/test_training_record_metadata.py tests/model/layers -q
```

13 focused tests and 2 neighboring layer tests pass. On unchanged upstream the
same focused cases give 12 failures and 1 passing legacy-mapping control.
The tests skip when optional processing dependencies are absent.

Coverage includes:

- Cluster keys produced by `cluster.py` consumed directly by `rcsb.py` for
  proteins, DNA, RNA and a single CCD ligand. MMseqs itself is replaced with a
  deterministic TSV-producing test double; clustering quality is not measured.
- Full-sequence keys preferred over legacy PDB/entity keys, legacy fallback,
  missing clusters, repeated/truncated chain names, and multi-residue ligands.
- Existing protein MSA files linked; missing files, directories and nucleotide
  files not linked. A real parsed A3M is saved and loaded through the training
  input loader, including both copies of the protein chain.
- Saved record metadata, numeric entity IDs, interface indices, and forwarding
  the MSA directory through serial and parallel scheduling branches. The latter
  replaces worker scheduling and Redis, rather than launching external services.
- A small mmCIF generated with Gemmi and parsed by the real preprocessing parser:
  a missing MSE residue retains its canonical M letter in the full sequence;
  the chain table includes the current `cyclic_period` field with default zero.
  This also catches an existing 9-versus-10-field table-construction failure.

The run used Python 3.12, NumPy 1.26.4, Gemmi 0.6.5 and RDKit in the existing
Boltz environment, with CUDA hidden. Targeted Ruff E4/E7/E9/F/I and
`git diff --check` pass. No Redis service, MMseqs executable, dataset download,
model training, GPU run or model-quality claim is required for these metadata
checks. The entire legacy inference/kernel test suite was not run.

Existing preprocessed records are skipped by the existing cache policy. Use a
fresh output directory to apply the metadata fix; existing models or datasets
are not silently rewritten.
