#!/usr/bin/env python3
"""
Neural Style Transfer via Diffusion Latent Steering for Boltz.

This script performs neural style transfer on protein structures by optimizing
the pairwise (z) latent representation to achieve two objectives:
1. Content: Minimize coordinate loss to a reference structure (CA atoms)
2. Style: Maximize confidence score (iPTM)

The approach:
1. Run the trunk modules to get initial latents (s, z)
2. Enable gradient computation on z
3. Iteratively:
   - Run diffusion and confidence from z
   - Compute content loss (RMSD to reference)
   - Compute style loss (negative iPTM)
   - Backpropagate and update z
"""

import argparse
import logging
from pathlib import Path
from typing import Optional

import torch
import yaml
from pytorch_lightning import Trainer

from boltz.data import const
from boltz.data.feature.batch import collate_fn
from boltz.data.feature.featurize import featurize
from boltz.data.parse.a3m import parse_a3m
from boltz.data.parse.fasta import parse_fasta
from boltz.main import setup_model

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


def load_reference_coords(
    pdb_path: Path,
    atom_selection: str = "CA",
) -> torch.Tensor:
    """
    Load reference coordinates from a PDB file.

    Parameters
    ----------
    pdb_path : Path
        Path to reference PDB file
    atom_selection : str
        Atom type to select (e.g., "CA" for alpha carbons)

    Returns
    -------
    torch.Tensor
        Reference coordinates [num_atoms, 3]
    """
    from Bio import PDB

    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("reference", pdb_path)

    coords = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if atom_selection in residue:
                    atom = residue[atom_selection]
                    coords.append(atom.get_coord())

    if not coords:
        msg = f"No {atom_selection} atoms found in {pdb_path}"
        raise ValueError(msg)

    return torch.tensor(coords, dtype=torch.float32)


def extract_ca_coords(
    atom_coords: torch.Tensor,
    feats: dict,
    atom_selection: str = "CA",
) -> torch.Tensor:
    """
    Extract specific atom coordinates from full atom representation.

    Parameters
    ----------
    atom_coords : torch.Tensor
        Full atom coordinates [num_atoms, 3]
    feats : dict
        Feature dictionary containing atom type information
    atom_selection : str
        Atom type to select (e.g., "CA")

    Returns
    -------
    torch.Tensor
        Selected atom coordinates [num_selected, 3]
    """
    # Get atom names from features
    atom_names = feats.get("atom_names", None)
    if atom_names is None:
        # Fallback: assume every Nth atom is the desired type
        # This is a simplification and may need adjustment
        logger.warning(
            "atom_names not found in features, using stride-based selection"
        )
        # For proteins, CA is typically the 2nd atom (index 1) in each residue
        # This is a rough approximation
        num_tokens = feats["token_pad_mask"].sum().item()
        stride = len(atom_coords) // num_tokens
        ca_indices = torch.arange(1, len(atom_coords), stride)
        return atom_coords[ca_indices]

    # Find indices of desired atoms
    ca_indices = []
    for i, name in enumerate(atom_names):
        if name == atom_selection:
            ca_indices.append(i)

    if not ca_indices:
        msg = f"No {atom_selection} atoms found in predicted structure"
        raise ValueError(msg)

    return atom_coords[torch.tensor(ca_indices)]


def compute_rmsd(
    coords1: torch.Tensor,
    coords2: torch.Tensor,
    align: bool = True,
) -> torch.Tensor:
    """
    Compute RMSD between two sets of coordinates.

    Parameters
    ----------
    coords1 : torch.Tensor
        First set of coordinates [N, 3]
    coords2 : torch.Tensor
        Second set of coordinates [N, 3]
    align : bool
        Whether to align structures before computing RMSD

    Returns
    -------
    torch.Tensor
        RMSD value (scalar tensor)
    """
    if coords1.shape != coords2.shape:
        msg = f"Coordinate shapes don't match: {coords1.shape} vs {coords2.shape}"
        raise ValueError(msg)

    if align:
        # Center both structures
        coords1_centered = coords1 - coords1.mean(dim=0, keepdim=True)
        coords2_centered = coords2 - coords2.mean(dim=0, keepdim=True)

        # Compute optimal rotation using Kabsch algorithm
        # H = coords2_centered.T @ coords1_centered
        # U, S, Vt = torch.linalg.svd(H)
        # R = U @ Vt

        # For simplicity, just use centered coordinates
        # (full alignment would require SVD which can be unstable with gradients)
        coords1 = coords1_centered
        coords2 = coords2_centered

    # Compute RMSD
    diff = coords1 - coords2
    msd = (diff ** 2).sum(dim=-1).mean()
    rmsd = torch.sqrt(msd)

    return rmsd


def neural_style_transfer(
    model,
    feats: dict,
    reference_coords: torch.Tensor,
    num_iterations: int = 100,
    learning_rate: float = 0.01,
    content_weight: float = 1.0,
    style_weight: float = 1.0,
    num_sampling_steps: int = 20,
    atom_selection: str = "CA",
    device: str = "cuda",
) -> dict:
    """
    Perform neural style transfer via diffusion latent steering.

    Parameters
    ----------
    model : Boltz2
        Loaded Boltz2 model
    feats : dict
        Input features
    reference_coords : torch.Tensor
        Reference coordinates for content loss [num_atoms, 3]
    num_iterations : int
        Number of optimization iterations
    learning_rate : float
        Learning rate for gradient descent
    content_weight : float
        Weight for content loss
    style_weight : float
        Weight for style loss (iPTM)
    num_sampling_steps : int
        Number of diffusion sampling steps
    atom_selection : str
        Atom type for content loss (e.g., "CA")
    device : str
        Device to run on

    Returns
    -------
    dict
        Results containing optimized coordinates and metrics
    """
    logger.info("Starting neural style transfer...")
    logger.info(f"Iterations: {num_iterations}")
    logger.info(f"Learning rate: {learning_rate}")
    logger.info(f"Content weight: {content_weight}")
    logger.info(f"Style weight: {style_weight}")

    # Move reference coords to device
    reference_coords = reference_coords.to(device)

    # Get initial trunk latents
    logger.info("Computing initial trunk latents...")
    with torch.no_grad():
        latents = model.get_trunk_latents(feats, recycling_steps=0)

    s = latents["s"]
    z = latents["z"]
    s_inputs = latents["s_inputs"]
    rel_pos_enc = latents["relative_position_encoding"]
    diffusion_conditioning = latents["diffusion_conditioning"]

    # Enable gradients on z (pairwise latent)
    z = z.detach().clone()
    z.requires_grad = True

    # Setup optimizer
    optimizer = torch.optim.Adam([z], lr=learning_rate)

    # Optimization loop
    results = {
        "iterations": [],
        "content_losses": [],
        "style_losses": [],
        "total_losses": [],
        "iptm_scores": [],
        "rmsd_values": [],
    }

    for iteration in range(num_iterations):
        optimizer.zero_grad()

        # Run diffusion and confidence from current z
        # Note: we recompute conditioning from the modified z
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

        # Extract CA coordinates (or other atom type)
        # This is a simplified version - may need adjustment based on actual feature structure
        try:
            pred_ca_coords = extract_ca_coords(pred_coords, feats, atom_selection)
        except Exception as e:
            logger.warning(f"Error extracting CA coords: {e}. Using all atoms.")
            pred_ca_coords = pred_coords

        # Ensure shapes match
        min_len = min(len(pred_ca_coords), len(reference_coords))
        pred_ca_coords = pred_ca_coords[:min_len]
        ref_coords = reference_coords[:min_len]

        # Compute content loss (RMSD to reference)
        content_loss = compute_rmsd(pred_ca_coords, ref_coords, align=True)

        # Compute style loss (negative iPTM to maximize iPTM)
        iptm = output.get("iptm", torch.tensor([0.0], device=device))[0]
        style_loss = -iptm

        # Total loss
        total_loss = content_weight * content_loss + style_weight * style_loss

        # Backpropagate
        total_loss.backward()
        optimizer.step()

        # Log metrics
        results["iterations"].append(iteration)
        results["content_losses"].append(content_loss.item())
        results["style_losses"].append(-style_loss.item())  # Store as positive iPTM
        results["total_losses"].append(total_loss.item())
        results["iptm_scores"].append(iptm.item())
        results["rmsd_values"].append(content_loss.item())

        if iteration % 10 == 0:
            logger.info(
                f"Iteration {iteration}: "
                f"Total Loss = {total_loss.item():.4f}, "
                f"RMSD = {content_loss.item():.4f}, "
                f"iPTM = {iptm.item():.4f}"
            )

    # Final evaluation
    logger.info("Running final evaluation...")
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

    results["final_coords"] = final_output["sample_atom_coords"]
    results["final_iptm"] = final_output.get("iptm", torch.tensor([0.0]))[0].item()
    results["optimized_z"] = z.detach()

    logger.info(f"Final iPTM: {results['final_iptm']:.4f}")

    return results


def main():
    """Main entry point for neural style transfer script."""
    parser = argparse.ArgumentParser(
        description="Neural style transfer via diffusion latent steering"
    )
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Input YAML file or FASTA file",
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
        help="Model checkpoint path",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("style_transfer_output"),
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
        "--atom-selection",
        type=str,
        default="CA",
        help="Atom type for content loss (e.g., CA, N, C)",
    )
    parser.add_argument(
        "--device",
        type=str,
        default="cuda" if torch.cuda.is_available() else "cpu",
        help="Device to run on",
    )

    args = parser.parse_args()

    # Create output directory
    args.output.mkdir(parents=True, exist_ok=True)

    # Load model
    logger.info(f"Loading model from {args.checkpoint}...")
    model = setup_model(
        checkpoint=str(args.checkpoint),
        devices=1,
        accelerator="gpu" if args.device == "cuda" else "cpu",
    )
    model = model.to(args.device)
    model.eval()

    # Load reference coordinates
    logger.info(f"Loading reference coordinates from {args.reference}...")
    reference_coords = load_reference_coords(args.reference, args.atom_selection)

    # Process input
    logger.info(f"Processing input from {args.input}...")
    if args.input.suffix == ".yaml":
        with open(args.input) as f:
            data = yaml.safe_load(f)
        # Process YAML input (simplified - may need adjustment)
        # This would need proper implementation based on Boltz's data processing
        raise NotImplementedError("YAML input processing not yet implemented")
    elif args.input.suffix in [".fasta", ".fa"]:
        # Process FASTA
        records = parse_fasta(args.input)
        # This is simplified - proper featurization would be needed
        raise NotImplementedError("FASTA input processing not yet implemented")
    else:
        msg = f"Unsupported input format: {args.input.suffix}"
        raise ValueError(msg)

    # For now, we'll need the user to provide processed features
    # TODO: Implement proper input processing
    logger.error(
        "Input processing not yet fully implemented. "
        "Please provide pre-processed features."
    )

    # Perform neural style transfer
    # This would be called once input processing is implemented:
    # results = neural_style_transfer(
    #     model=model,
    #     feats=feats,
    #     reference_coords=reference_coords,
    #     num_iterations=args.num_iterations,
    #     learning_rate=args.learning_rate,
    #     content_weight=args.content_weight,
    #     style_weight=args.style_weight,
    #     num_sampling_steps=args.num_sampling_steps,
    #     atom_selection=args.atom_selection,
    #     device=args.device,
    # )

    # Save results
    # output_path = args.output / "optimized_structure.pdb"
    # logger.info(f"Saving results to {output_path}...")
    # ... save PDB and metrics ...

    logger.info("Done!")


if __name__ == "__main__":
    main()
