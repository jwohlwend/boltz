#!/usr/bin/env python3
"""
Simplified Neural Style Transfer Example for Boltz.

This script demonstrates the neural style transfer concept using pre-processed
Boltz features. It optimizes the pairwise (z) latent to achieve:
1. Content objective: Minimize coordinate RMSD to a reference structure
2. Style objective: Maximize iPTM confidence score

Usage:
    python neural_style_transfer_example.py \\
        --features input_features.pt \\
        --reference reference.pdb \\
        --checkpoint boltz2_conf.ckpt \\
        --output results/

The input features should be a PyTorch .pt file containing a dictionary with
processed Boltz features (obtained from the standard Boltz prediction pipeline).
"""

import argparse
import logging
from pathlib import Path

import torch
import numpy as np

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


def load_reference_ca_coords(pdb_path: Path) -> torch.Tensor:
    """Load CA coordinates from a PDB file."""
    from Bio import PDB

    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("reference", pdb_path)

    coords = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if "CA" in residue:
                    coords.append(residue["CA"].get_coord())

    if not coords:
        raise ValueError(f"No CA atoms found in {pdb_path}")

    return torch.tensor(coords, dtype=torch.float32)


def extract_ca_coords_from_atoms(
    atom_coords: torch.Tensor,
    atom_pad_mask: torch.Tensor,
    num_tokens: int,
) -> torch.Tensor:
    """
    Extract CA coordinates from full atom representation.

    Assumes CA is the second atom (index 1) in each residue/token.
    """
    # Reshape atom coords to [num_tokens, atoms_per_token, 3]
    num_atoms = atom_pad_mask.sum().item()
    atoms_per_token = num_atoms // num_tokens

    # Get CA indices (assuming CA is at position 1 in each residue)
    ca_indices = torch.arange(1, num_atoms, atoms_per_token, device=atom_coords.device)
    ca_coords = atom_coords[ca_indices]

    return ca_coords


def compute_aligned_rmsd(
    coords1: torch.Tensor,
    coords2: torch.Tensor,
) -> torch.Tensor:
    """
    Compute RMSD between two coordinate sets after centering.

    Note: For gradient-based optimization, we use simple centering rather
    than full Kabsch alignment to avoid SVD instabilities.
    """
    # Center both structures
    coords1_centered = coords1 - coords1.mean(dim=0, keepdim=True)
    coords2_centered = coords2 - coords2.mean(dim=0, keepdim=True)

    # Compute RMSD
    diff = coords1_centered - coords2_centered
    msd = (diff ** 2).sum(dim=-1).mean()
    rmsd = torch.sqrt(msd + 1e-8)  # Add epsilon for numerical stability

    return rmsd


def optimize_latents(
    model,
    feats: dict,
    reference_ca_coords: torch.Tensor,
    num_iterations: int = 100,
    learning_rate: float = 0.01,
    content_weight: float = 1.0,
    style_weight: float = 1.0,
    num_sampling_steps: int = 20,
    device: str = "cuda",
):
    """
    Optimize the z (pairwise) latent via gradient descent.

    Objective:
        L = content_weight * RMSD(pred, ref) - style_weight * iPTM
    """
    logger.info("=" * 80)
    logger.info("Starting Neural Style Transfer Optimization")
    logger.info("=" * 80)
    logger.info(f"Iterations: {num_iterations}")
    logger.info(f"Learning rate: {learning_rate}")
    logger.info(f"Content weight: {content_weight}")
    logger.info(f"Style weight: {style_weight}")
    logger.info(f"Sampling steps: {num_sampling_steps}")
    logger.info("=" * 80)

    # Move reference to device
    reference_ca_coords = reference_ca_coords.to(device)

    # Get initial trunk latents
    logger.info("Computing initial trunk latents...")
    with torch.no_grad():
        latents = model.get_trunk_latents(feats, recycling_steps=0)

    s = latents["s"]
    z = latents["z"]
    s_inputs = latents["s_inputs"]
    rel_pos_enc = latents["relative_position_encoding"]

    logger.info(f"Latent shapes: s={s.shape}, z={z.shape}")

    # Enable gradients on z
    z = z.detach().clone()
    z.requires_grad = True

    # Setup optimizer
    optimizer = torch.optim.Adam([z], lr=learning_rate)

    # Track metrics
    metrics = {
        "iteration": [],
        "total_loss": [],
        "content_loss": [],
        "style_loss": [],
        "iptm": [],
        "rmsd": [],
    }

    # Get number of tokens for CA extraction
    num_tokens = feats["token_pad_mask"].sum().item()

    # Optimization loop
    logger.info("\nStarting optimization...")
    logger.info("-" * 80)

    for iteration in range(num_iterations):
        optimizer.zero_grad()

        # Run diffusion and confidence from current z
        output = model.run_from_latents(
            s=s,
            z=z,
            s_inputs=s_inputs,
            relative_position_encoding=rel_pos_enc,
            feats=feats,
            num_sampling_steps=num_sampling_steps,
            diffusion_samples=1,
            recompute_conditioning=True,
        )

        # Extract predicted coordinates
        pred_coords = output["sample_atom_coords"][0]  # [num_atoms, 3]

        # Extract CA coordinates
        pred_ca_coords = extract_ca_coords_from_atoms(
            pred_coords,
            feats["atom_pad_mask"][0],
            num_tokens,
        )

        # Match lengths
        min_len = min(len(pred_ca_coords), len(reference_ca_coords))
        pred_ca_coords = pred_ca_coords[:min_len]
        ref_ca_coords = reference_ca_coords[:min_len]

        # Compute content loss (RMSD)
        content_loss = compute_aligned_rmsd(pred_ca_coords, ref_ca_coords)

        # Compute style loss (negative iPTM to maximize iPTM)
        iptm = output.get("iptm", torch.tensor([0.0], device=device))
        if iptm.dim() == 0:
            iptm = iptm.unsqueeze(0)
        iptm_value = iptm[0]
        style_loss = -iptm_value

        # Total loss
        total_loss = content_weight * content_loss + style_weight * style_loss

        # Backpropagate
        total_loss.backward()

        # Gradient clipping for stability
        torch.nn.utils.clip_grad_norm_([z], max_norm=1.0)

        optimizer.step()

        # Log metrics
        metrics["iteration"].append(iteration)
        metrics["total_loss"].append(total_loss.item())
        metrics["content_loss"].append(content_loss.item())
        metrics["style_loss"].append(-style_loss.item())  # Store as positive iPTM
        metrics["iptm"].append(iptm_value.item())
        metrics["rmsd"].append(content_loss.item())

        if iteration % 10 == 0 or iteration == num_iterations - 1:
            logger.info(
                f"Iter {iteration:3d} | "
                f"Loss: {total_loss.item():7.4f} | "
                f"RMSD: {content_loss.item():6.3f} Å | "
                f"iPTM: {iptm_value.item():5.3f}"
            )

    logger.info("-" * 80)
    logger.info("Optimization complete!")

    # Final evaluation
    logger.info("\nGenerating final structure...")
    with torch.no_grad():
        final_output = model.run_from_latents(
            s=s,
            z=z,
            s_inputs=s_inputs,
            relative_position_encoding=rel_pos_enc,
            feats=feats,
            num_sampling_steps=num_sampling_steps,
            diffusion_samples=1,
            recompute_conditioning=True,
        )

    logger.info("=" * 80)
    logger.info("Final Results:")
    logger.info(f"  iPTM: {final_output.get('iptm', [0])[0]:.4f}")
    logger.info(f"  pTM:  {final_output.get('ptm', [0])[0]:.4f}")
    logger.info(f"  pLDDT: {final_output.get('complex_plddt', [0])[0]:.4f}")
    logger.info("=" * 80)

    return {
        "metrics": metrics,
        "final_coords": final_output["sample_atom_coords"],
        "final_output": final_output,
        "optimized_z": z.detach(),
        "optimized_s": s,
        "s_inputs": s_inputs,
    }


def save_results(results: dict, output_dir: Path, feats: dict):
    """Save optimization results and metrics."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # Save metrics
    metrics_file = output_dir / "metrics.npz"
    np.savez(
        metrics_file,
        **{k: np.array(v) for k, v in results["metrics"].items()}
    )
    logger.info(f"Saved metrics to {metrics_file}")

    # Save optimized latents
    latents_file = output_dir / "optimized_latents.pt"
    torch.save({
        "z": results["optimized_z"],
        "s": results["optimized_s"],
        "s_inputs": results["s_inputs"],
    }, latents_file)
    logger.info(f"Saved optimized latents to {latents_file}")

    # Save final coordinates
    coords_file = output_dir / "final_coords.pt"
    torch.save(results["final_coords"], coords_file)
    logger.info(f"Saved final coordinates to {coords_file}")

    # Save final output dict
    output_file = output_dir / "final_output.pt"
    torch.save(results["final_output"], output_file)
    logger.info(f"Saved final output to {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Neural style transfer via diffusion latent steering"
    )
    parser.add_argument(
        "--features",
        type=Path,
        required=True,
        help="Pre-processed Boltz features (.pt file)",
    )
    parser.add_argument(
        "--reference",
        type=Path,
        required=True,
        help="Reference PDB file for content loss",
    )
    parser.add_argument(
        "--checkpoint",
        type=Path,
        required=True,
        help="Boltz2 model checkpoint",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("style_transfer_results"),
        help="Output directory",
    )
    parser.add_argument(
        "--num-iterations",
        type=int,
        default=100,
        help="Number of optimization iterations",
    )
    parser.add_argument(
        "--learning-rate",
        type=float,
        default=0.01,
        help="Learning rate for gradient descent",
    )
    parser.add_argument(
        "--content-weight",
        type=float,
        default=1.0,
        help="Weight for content loss (RMSD)",
    )
    parser.add_argument(
        "--style-weight",
        type=float,
        default=1.0,
        help="Weight for style loss (iPTM)",
    )
    parser.add_argument(
        "--num-sampling-steps",
        type=int,
        default=20,
        help="Number of diffusion sampling steps",
    )
    parser.add_argument(
        "--device",
        type=str,
        default="cuda" if torch.cuda.is_available() else "cpu",
        help="Device to run on",
    )

    args = parser.parse_args()

    # Load model
    logger.info(f"Loading model from {args.checkpoint}...")
    from boltz.model.models.boltz2 import Boltz2

    model = Boltz2.load_from_checkpoint(str(args.checkpoint))
    model = model.to(args.device)
    model.eval()

    # Load features
    logger.info(f"Loading features from {args.features}...")
    feats = torch.load(args.features)

    # Move features to device
    for key, value in feats.items():
        if isinstance(value, torch.Tensor):
            feats[key] = value.to(args.device)

    # Load reference coordinates
    logger.info(f"Loading reference from {args.reference}...")
    reference_coords = load_reference_ca_coords(args.reference)

    # Run optimization
    results = optimize_latents(
        model=model,
        feats=feats,
        reference_ca_coords=reference_coords,
        num_iterations=args.num_iterations,
        learning_rate=args.learning_rate,
        content_weight=args.content_weight,
        style_weight=args.style_weight,
        num_sampling_steps=args.num_sampling_steps,
        device=args.device,
    )

    # Save results
    save_results(results, args.output, feats)

    logger.info(f"\nResults saved to {args.output}")
    logger.info("Done!")


if __name__ == "__main__":
    main()
