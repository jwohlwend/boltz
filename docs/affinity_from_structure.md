# Experimental affinity scoring of a supplied holo structure

`Boltz2.forward(..., affinity_from_input_structure=True)` evaluates the affinity
heads using the coordinates supplied in `feats['coords']`. It retains the
learned trunk, both affinity heads, their ensemble and the optional existing
molecular-weight correction. It bypasses diffusion conditioning, coordinate
sampling, confidence prediction and best-IPTM pose selection.

This is an opt-in Python API and example for discussion in
[issue #702](https://github.com/jwohlwend/boltz/issues/702), not a new flag of the
standard `boltz predict` CLI. The default forward path and checkpoint/state-dict
format are unchanged. There is no structure relaxation or affinity fine-tuning.

## API contract

Load a Boltz-2 **affinity** checkpoint, call `eval()`, and supply ordinary Boltz-2
features prepared for that same complex and crop. Input coordinates must be
float32, finite, and shaped `[1, 1, atoms, 3]` (one complex, one conformation).
Every real atom must have a resolved coordinate; unresolved padding is allowed.
Atom order must match the chemical, token, mapping and mask features. The
receptor must be protein and exactly one small-molecule ligand must be marked
by `affinity_token_mask`. RNA/DNA, ions, additional ligands and other cofactors
are outside this initial API contract.

```python
model.eval()
with torch.inference_mode():
    result = model(
        feats,
        recycling_steps=5,
        affinity_from_input_structure=True,
    )
```

The usual affinity fields are returned. `sample_atom_coords` contains the
unchanged input conformation with shape `[1, atoms, 3]`; it denotes the supplied
pose in this mode. Confidence/PAE/IPTM outputs are absent. Diffusion sample count,
sampling steps and sequential confidence options do not generate additional
structures in this mode. The input feature dictionary is not mutated.

For CUDA, use the normal prediction precision context (bf16 autocast for the
trunk; the forward already disables autocast around the affinity heads):

```python
with torch.inference_mode(), torch.autocast('cuda', dtype=torch.bfloat16):
    result = model(feats, recycling_steps=5, affinity_from_input_structure=True)
```

## mmCIF example

The example uses the existing Boltz mmCIF parser and CCD molecules, tokenizer,
stock affinity crop (200 protein tokens, 256 total tokens, 2048 atoms), and
featurizer. It uses single-sequence features; no MSA server is called. It requires
the existing `boltz2_aff.ckpt` and extracted `mols/` under the specified cache.

```bash
curl -L https://files.rcsb.org/download/1HSG.cif -o 1HSG.cif
python examples/affinity_from_structure.py 1HSG.cif \
  --protein-chain A --protein-chain B --ligand-chain C \
  --cache ~/.boltz --output known-structure-affinity.json
```

Use mmCIF **label_asym_id/subchain** identifiers, not author chain identifiers:
in 1HSG the MK1 ligand is subchain C even though its author chain is B. Both
protein subchains are selected explicitly. No biological assembly is expanded.
`--prepare-only` checks parsing/chemistry and writes input metadata without
loading model weights. `--mw-correction` enables the existing fixed correction;
the example defaults to no correction, as does the prediction CLI. Use a fresh
output file for each invocation.

The example refuses multi-model structures, missing model atoms, incomplete
ligand crops, and covalently attached ligands. A complete structure prepared by
an external tool may be required. For example, a model CCD can retain a hydrogen
that is absent in an experimental X-ray file; the example rejects that case
rather than replacing its coordinate with zero. It does not reconstruct missing
residues, infer coordinates from SMILES, or match a separately ordered ligand
array by index. The released parser chooses alternate conformations, removes
waters/hydrogens and applies its existing residue conventions (including MSE
and symmetric arginine atom naming). Coordinates of the retained atoms undergo
only the featurizer's common centering translation, within FP32 rounding; the
forward does not modify them further. The output records selected/omitted
subchains, input hash, crop dimensions, molecular weight and precision.

## Interpretation and testing

The trained weights are unchanged, but scoring a fixed experimental geometry
changes the standard inference workflow. This feature does **not** establish
better affinity accuracy, calibration, or equivalence to affinity from a
generated pose. Single-sequence context, missing waters and cofactors, the
chosen biological assembly, and structure preparation can all affect results.
The Boltz-2 paper notes limitations for cofactors and multimeric partners.

Focused CPU tests exercise the real forward with small instrumented modules:
both heads receive the supplied tensor, diffusion/confidence are never called,
the ensemble/MW calculation is retained, invalid inputs fail early, and the
default path still selects the highest-IPTM generated pose. Run:

```bash
pytest -q tests/model/test_input_structure_affinity.py
```

These focused checks are not the entire upstream test suite. A real-checkpoint
GPU smoke test additionally needs to verify routing, coordinate preservation,
and default behavior; a scientific accuracy study requires more complexes and
experimental affinity labels.
