# Neural Style Transfer via Diffusion Latent Steering

This document describes the neural style transfer functionality for Boltz, which allows you to optimize protein structures by manipulating the internal latent representations.

## Overview

Neural style transfer for protein structures works by optimizing the pairwise (z) latent representation to achieve two objectives:

1. **Content Objective**: Minimize coordinate RMSD to a reference structure (using CA atoms)
2. **Style Objective**: Maximize a confidence score (iPTM - interface predicted TM-score)

This is analogous to neural style transfer in computer vision, where:
- **Content** = structural coordinates (what the protein looks like spatially)
- **Style** = confidence/quality metrics (how "good" the structure is according to the model)

## How It Works

### Architecture Modifications

The `Boltz2` model has been extended with two new methods:

1. **`get_trunk_latents()`**: Runs the trunk modules (embedder, MSA, pairformer) and returns:
   - `s`: Single representation [batch, num_tokens, token_s]
   - `z`: Pair representation [batch, num_tokens, num_tokens, token_z]
   - Diffusion conditioning
   - Other intermediate values

2. **`run_from_latents()`**: Takes latent representations and runs:
   - Diffusion module (structure prediction)
   - Confidence module (iPTM, pTM, pLDDT, etc.)

### Optimization Process

```
┌─────────────────┐
│  Input Sequence │
└────────┬────────┘
         │
         v
┌─────────────────┐
│  Trunk Modules  │  ← InputEmbedder, MSA, Pairformer
│   (Fixed/NoGrad)│
└────────┬────────┘
         │
         v
┌─────────────────┐
│  Latents: s, z  │  ← z becomes optimizable
└────────┬────────┘
         │
         │  Gradient Descent Loop
         │  ┌──────────────────────────┐
         └─→│  1. Recompute Conditioning│
            │  2. Run Diffusion Sampling│
            │  3. Run Confidence Heads  │
            │  4. Compute Losses:       │
            │     - Content: RMSD       │
            │     - Style: -iPTM        │
            │  5. Backprop to z         │
            │  6. Update z              │
            └──────────────────────────┘
                      │
                      v
            ┌─────────────────┐
            │ Optimized       │
            │ Structure       │
            └─────────────────┘
```

## Usage

### Option 1: Using Pre-processed Features (Recommended)

First, run a standard Boltz prediction to get processed features:

```bash
# Run Boltz prediction (this generates features)
boltz predict input.yaml --checkpoint boltz2_conf.ckpt

# The features will be in the prediction output directory
# You can extract and save them for reuse
```

Then run neural style transfer:

```bash
python scripts/neural_style_transfer_example.py \
    --features processed_features.pt \
    --reference reference_structure.pdb \
    --checkpoint boltz2_conf.ckpt \
    --output results/ \
    --num-iterations 100 \
    --learning-rate 0.01 \
    --content-weight 1.0 \
    --style-weight 1.0 \
    --num-sampling-steps 20
```

### Option 2: Full Pipeline (Advanced)

```bash
python scripts/neural_style_transfer.py \
    --input input.yaml \
    --reference reference.pdb \
    --checkpoint boltz2_conf.ckpt \
    --output results/ \
    --num-iterations 100 \
    --learning-rate 0.01 \
    --content-weight 1.0 \
    --style-weight 1.0
```

## Parameters

### Optimization Parameters

- `--num-iterations`: Number of gradient descent iterations (default: 100)
- `--learning-rate`: Step size for gradient descent (default: 0.01)
- `--content-weight`: Weight for content loss (RMSD) (default: 1.0)
- `--style-weight`: Weight for style loss (iPTM) (default: 1.0)

### Diffusion Parameters

- `--num-sampling-steps`: Number of diffusion denoising steps (default: 20)
  - More steps = higher quality but slower
  - Typical range: 10-50

### Atom Selection

- `--atom-selection`: Atom type for content loss (default: "CA")
  - Options: "CA" (alpha carbon), "N", "C", "CB", etc.
  - CA is recommended as it captures backbone structure

## Loss Functions

### Content Loss

```python
content_loss = RMSD(pred_CA_coords, reference_CA_coords)
```

- Measures structural similarity to reference
- Uses CA atoms by default
- Structures are centered before RMSD computation
- Lower is better

### Style Loss

```python
style_loss = -iPTM
```

- Maximizes interface predicted TM-score
- iPTM measures model confidence in inter-chain contacts
- Can also use pTM (overall TM-score) or pLDDT
- Higher iPTM is better (so we minimize negative iPTM)

### Total Loss

```python
total_loss = content_weight * content_loss + style_weight * style_loss
```

## Output

The script generates the following outputs in the specified output directory:

1. **`metrics.npz`**: NumPy archive containing optimization metrics
   - `iteration`: Iteration numbers
   - `total_loss`: Total loss at each iteration
   - `content_loss`: RMSD values
   - `style_loss`: iPTM values
   - `rmsd`: RMSD values (same as content_loss)
   - `iptm`: iPTM scores

2. **`optimized_latents.pt`**: PyTorch file with optimized latents
   - `z`: Optimized pairwise representation
   - `s`: Single representation
   - `s_inputs`: Input embeddings

3. **`final_coords.pt`**: Final predicted coordinates

4. **`final_output.pt`**: Complete output dictionary including all confidence scores

## Visualization

You can visualize the optimization progress:

```python
import numpy as np
import matplotlib.pyplot as plt

# Load metrics
data = np.load('results/metrics.npz')

# Plot RMSD over iterations
plt.figure(figsize=(12, 4))

plt.subplot(1, 3, 1)
plt.plot(data['iteration'], data['rmsd'])
plt.xlabel('Iteration')
plt.ylabel('RMSD (Å)')
plt.title('Content Loss (RMSD)')

plt.subplot(1, 3, 2)
plt.plot(data['iteration'], data['iptm'])
plt.xlabel('Iteration')
plt.ylabel('iPTM')
plt.title('Style Score (iPTM)')

plt.subplot(1, 3, 3)
plt.plot(data['iteration'], data['total_loss'])
plt.xlabel('Iteration')
plt.ylabel('Total Loss')
plt.title('Total Loss')

plt.tight_layout()
plt.savefig('optimization_progress.png')
```

## Tips and Best Practices

### Learning Rate

- Start with `lr=0.01` for most cases
- Increase to `0.05-0.1` for faster convergence (may be less stable)
- Decrease to `0.001-0.005` if optimization is unstable

### Loss Weights

- **Equal weights** (`content_weight=1.0`, `style_weight=1.0`): Balanced optimization
- **High content weight** (`content_weight=10.0`, `style_weight=1.0`): Prioritize matching reference structure
- **High style weight** (`content_weight=1.0`, `style_weight=10.0`): Prioritize model confidence

### Number of Iterations

- 50-100 iterations usually sufficient for convergence
- Monitor the loss curves to check for convergence
- If loss plateaus, optimization has converged

### Diffusion Sampling Steps

- More steps = better quality predictions but slower
- 20 steps is a good default
- Can reduce to 10 for faster experimentation
- Increase to 40-50 for final high-quality results

## Advanced Usage

### Custom Loss Functions

You can modify the optimization objective by editing the loss computation in the scripts:

```python
# Example: Add a pLDDT objective
plddt = output["complex_plddt"][0]
plddt_loss = -plddt

total_loss = (
    content_weight * content_loss +
    style_weight * style_loss +
    plddt_weight * plddt_loss
)
```

### Optimizing Different Latents

The current implementation optimizes the `z` (pairwise) latent. You can also optimize:

- `s` (single representation): Contains per-residue features
- Both `s` and `z`: Joint optimization
- Conditioning tensors: Direct manipulation of diffusion inputs

Example:

```python
# Enable gradients on both s and z
s = s.detach().clone()
s.requires_grad = True
z = z.detach().clone()
z.requires_grad = True

optimizer = torch.optim.Adam([s, z], lr=learning_rate)
```

### Different Confidence Scores

You can optimize for different confidence metrics:

- **iPTM**: Interface TM-score (inter-chain contacts)
- **pTM**: Overall TM-score (full structure)
- **pLDDT**: Per-residue confidence
- **iPLDDT**: Interface-weighted pLDDT

Simply replace `iptm` with your desired metric in the loss computation.

## Troubleshooting

### "Out of memory" errors

- Reduce `num_sampling_steps`
- Use CPU instead of GPU (slower but uses system RAM)
- Process one sample at a time

### Optimization not converging

- Decrease learning rate
- Increase number of iterations
- Check that reference structure is compatible with input sequence
- Adjust loss weights

### Poor final structure quality

- Increase `num_sampling_steps` for better diffusion sampling
- Try different weight combinations
- Check that reference structure is reasonable

## Citation

If you use this neural style transfer functionality in your research, please cite:

```bibtex
@article{wohlwend2024boltz1,
  title={Boltz-1: Democratizing biomolecular interaction modeling},
  author={Wohlwend, Jeremy and ...},
  journal={bioRxiv},
  year={2024}
}
```

## License

This extension follows the same license as Boltz (see main LICENSE file).
