# started from code from https://github.com/lucidrains/alphafold3-pytorch, MIT License, Copyright (c) 2024 Phil Wang

import einx
import torch
import torch.nn.functional as F
from einops import einsum, rearrange


def weighted_rigid_align(
    true_coords,  # Float['b n 3'],       #  true coordinates
    pred_coords,  # Float['b n 3'],       # predicted coordinates
    weights,  # Float['b n'],             # weights for each atom
    mask,  # Bool['b n'] | None = None    # mask for variable lengths
):  # -> Float['b n 3']:
    """Algorithm 28 : note there is a problem with the pseudocode in the paper where predicted and
    GT are swapped in algorithm 28, but correct in equation (2)."""

    out_shape = torch.broadcast_shapes(true_coords.shape, pred_coords.shape)
    *batch_size, num_points, dim = out_shape
    weights = (mask * weights).unsqueeze(-1)

    # Compute weighted centroids
    true_centroid = (true_coords * weights).sum(dim=-2, keepdim=True) / weights.sum(
        dim=-2, keepdim=True
    )
    pred_centroid = (pred_coords * weights).sum(dim=-2, keepdim=True) / weights.sum(
        dim=-2, keepdim=True
    )

    # Center the coordinates
    true_coords_centered = true_coords - true_centroid
    pred_coords_centered = pred_coords - pred_centroid

    if torch.any(mask.sum(dim=-1) < (dim + 1)):
        print(
            "Warning: The size of one of the point clouds is <= dim+1. "
            + "`WeightedRigidAlign` cannot return a unique rotation."
        )

    # Compute the weighted covariance matrix
    cov_matrix = einsum(
        weights * pred_coords_centered,
        true_coords_centered,
        "... n i, ... n j -> ... i j",
    )

    # Compute the SVD of the covariance matrix, required float32 for svd and determinant
    original_dtype = cov_matrix.dtype
    cov_matrix_32 = cov_matrix.to(dtype=torch.float32)

    U, S, V = torch.linalg.svd(
        cov_matrix_32, driver="gesvd" if cov_matrix_32.is_cuda else None
    )
    V = V.mH

    # Catch ambiguous rotation by checking the magnitude of singular values
    if (S.abs() <= 1e-15).any() and not (num_points < (dim + 1)):
        print(
            "Warning: Excessively low rank of "
            + "cross-correlation between aligned point clouds. "
            + "`WeightedRigidAlign` cannot return a unique rotation."
        )

    # Compute the rotation matrix
    rot_matrix = torch.einsum("... i j, ... k j -> ... i k", U, V).to(
        dtype=torch.float32
    )

    # Ensure proper rotation matrix with determinant 1
    F = torch.eye(dim, dtype=cov_matrix_32.dtype, device=cov_matrix.device)[
        None
    ].repeat(*batch_size, 1, 1)
    F[..., -1, -1] = torch.det(rot_matrix)
    rot_matrix = einsum(U, F, V, "... i j, ... j k, ... l k -> ... i l")
    rot_matrix = rot_matrix.to(dtype=original_dtype)

    # Apply the rotation and translation
    aligned_coords = (
        einsum(true_coords_centered, rot_matrix, "... n i, ... j i -> ... n j")
        + pred_centroid
    )
    aligned_coords.detach_()

    return aligned_coords


def smooth_lddt_loss(
    pred_coords,  # Float['b n 3'],
    true_coords,  # Float['b n 3'],
    is_nucleotide,  # Bool['b n'],
    coords_mask,  # Bool['b n'] | None = None,
    nucleic_acid_cutoff: float = 30.0,
    other_cutoff: float = 15.0,
    multiplicity: int = 1,
):  # -> Float['']:
    """Algorithm 27
    pred_coords: predicted coordinates
    true_coords: true coordinates
    Note: for efficiency pred_coords is the only one with the multiplicity expanded
    TODO: add weighing which overweight the smooth lddt contribution close to t=0 (not present in the paper)
    """
    lddt = []
    for i in range(true_coords.shape[0]):
        true_dists = torch.cdist(true_coords[i], true_coords[i])

        is_nucleotide_i = is_nucleotide[i // multiplicity]
        coords_mask_i = coords_mask[i // multiplicity]

        is_nucleotide_pair = is_nucleotide_i.unsqueeze(-1).expand(
            -1, is_nucleotide_i.shape[-1]
        )

        mask = is_nucleotide_pair * (true_dists < nucleic_acid_cutoff).float()
        mask += (1 - is_nucleotide_pair) * (true_dists < other_cutoff).float()
        mask *= 1 - torch.eye(pred_coords.shape[1], device=pred_coords.device)
        mask *= coords_mask_i.unsqueeze(-1)
        mask *= coords_mask_i.unsqueeze(-2)

        valid_pairs = mask.nonzero()
        true_dists_i = true_dists[valid_pairs[:, 0], valid_pairs[:, 1]]

        pred_coords_i1 = pred_coords[i, valid_pairs[:, 0]]
        pred_coords_i2 = pred_coords[i, valid_pairs[:, 1]]
        pred_dists_i = F.pairwise_distance(pred_coords_i1, pred_coords_i2)

        dist_diff_i = torch.abs(true_dists_i - pred_dists_i)

        eps_i = (
            F.sigmoid(0.5 - dist_diff_i)
            + F.sigmoid(1.0 - dist_diff_i)
            + F.sigmoid(2.0 - dist_diff_i)
            + F.sigmoid(4.0 - dist_diff_i)
        ) / 4.0

        lddt_i = eps_i.sum() / (valid_pairs.shape[0] + 1e-5)
        lddt.append(lddt_i)

    # average over batch & multiplicity
    return 1.0 - torch.stack(lddt, dim=0).mean(dim=0)
