# SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
# SPDX-License-Identifier: MIT
#
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.


"""Tests for pair_mask factorization and DTensor lddt_resolved_token computation.

This module tests:
1. Whether the pair_mask construction in serial plddt_loss can be factorized
2. Whether cdist_lddt with factorized masks matches torch.cdist + lddt_dist
3. Whether DTensor lddt_resolved_token() matches serial reference

The pair_mask factorization:
    pair_mask = atom_mask[:,:,None] & atom_mask[:,None,:]  # outer product
    pair_mask = pair_mask & ~eye  # remove diagonal
    pair_mask = einsum("bnm,bkm->bnk", pair_mask, r_set_to_rep_atom)
    pair_mask = bmm(token_to_rep_atom, pair_mask)

Can be factorized into:
    mask_row = bmm(token_to_rep_atom, atom_mask)
    mask_col = bmm(r_set_to_rep_atom, atom_mask)
    factorized_mask = mask_row[:,:,None] & mask_col[:,None,:] & diagonal_mask

Where diagonal_mask handles the atom-level self-pair exclusion.
"""

import pytest
import torch
from torch.distributed.tensor import Replicate, Shard

from boltz.data import const
from boltz.distributed.comm import TransposeComm
from boltz.distributed.data.utils import distribute_features
from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.loss.confidencev2 import lddt_resolved_token, plddt_loss
from boltz.distributed.model.loss.triton.cdist_lddt import cdist_lddt
from boltz.model.loss.confidencev2 import lddt_dist
from boltz.model.loss.confidencev2 import plddt_loss as serial_plddt_loss
from boltz.testing.utils import (
    distribute_atom_features,
    random_features,
    spawn_multiprocessing,
)


def test_pair_mask_factorization():
    """Test pair_mask factorization with realistic features from random_features.

    This test checks if the factorized mask approach can reproduce
    the serial pair_mask construction using proper block-diagonal
    token_to_rep_atom and r_set_to_rep_atom matrices.
    """
    if not torch.cuda.is_available():
        pytest.skip("CUDA not available")

    device = torch.device("cuda")
    torch.manual_seed(42)
    dtype = torch.float32  # Use FP32 for testing

    # Test dimensions
    B = 2
    N_token = 160
    n_atoms_per_token_min, n_atoms_per_token_max = 1, 3
    # N_atom should be >= N_token * max_atoms_per_token to fit all atoms
    N_atom = N_token * (n_atoms_per_token_min + n_atoms_per_token_max) // 2  # 32

    # Generate realistic features using random_features
    feats = random_features(
        size_batch=B,
        n_tokens=N_token,
        n_atoms=N_atom,
        n_msa=1,
        atom_counts_per_token_range=(n_atoms_per_token_min, n_atoms_per_token_max),
        device=device,
        float_value_range=(-0.1, 0.1),
        selected_keys=[
            "token_to_rep_atom",
            "r_set_to_rep_atom",
        ],
    )

    # Extract features and convert to boolean where appropriate
    token_to_rep_atom = feats["token_to_rep_atom"].bool()  # [B, N_token, N_atom]
    r_set_to_rep_atom = feats["r_set_to_rep_atom"].bool()  # [B, N_R, N_atom]

    N_atom_actual = token_to_rep_atom.shape[2]

    # Create random atom_mask [B, N_atom] - boolean
    atom_mask = torch.randint(0, 2, (B, N_atom_actual), device=device, dtype=torch.bool)

    # ========== Serial pair_mask construction (using boolean ops) ==========
    # Step 1: Outer product of atom_mask
    pair_mask = atom_mask.unsqueeze(-1) & atom_mask.unsqueeze(-2)  # [B, N_atom, N_atom]

    # Step 2: Remove diagonal (self-pairs in atom space)
    diag_mask = ~torch.eye(N_atom_actual, device=device, dtype=torch.bool)
    pair_mask = pair_mask & diag_mask[None, :, :]

    # Step 3: Project columns via r_set_to_rep_atom (need float for einsum)
    pair_mask_float = pair_mask.to(dtype=dtype)
    r_set_float = r_set_to_rep_atom.to(dtype=dtype)
    pair_mask_float = torch.einsum("bnm,bkm->bnk", pair_mask_float, r_set_float)  # [B, N_atom, N_R]

    # Step 4: Project rows via token_to_rep_atom
    token_float = token_to_rep_atom.to(dtype=dtype)
    pair_mask_serial = torch.bmm(token_float, pair_mask_float)  # [B, N_token, N_R]

    # Convert back to boolean (values are 0 or 1 due to one-hot matrices)
    pair_mask_serial = pair_mask_serial.bool()

    # ========== Factorized mask construction ==========
    # mask_row[b, t] = True if rep_atom of token t is resolved
    mask_row = torch.bmm(token_float, atom_mask.unsqueeze(-1).to(dtype=dtype)).squeeze(-1).bool()  # [B, N_token]

    # mask_col[b, r] = True if rep_atom of R-set r is resolved
    mask_col = torch.bmm(r_set_float, atom_mask.unsqueeze(-1).to(dtype=dtype)).squeeze(-1).bool()  # [B, N_R]

    # Outer product of factorized masks
    factorized_outer = mask_row.unsqueeze(-1) & mask_col.unsqueeze(-2)  # [B, N_token, N_R]

    # Diagonal mask: exclude pairs where rep_atom_token == rep_atom_r_set
    rep_atom_token = token_to_rep_atom.int().argmax(dim=-1)  # [B, N_token]
    rep_atom_r_set = r_set_to_rep_atom.int().argmax(dim=-1)  # [B, N_R]

    # diagonal_mask[b, t, r] = True if rep_atom_token[b,t] != rep_atom_r_set[b,r]
    diagonal_mask = rep_atom_token.unsqueeze(-1) != rep_atom_r_set.unsqueeze(-2)  # [B, N_token, N_R]

    # Factorized mask with diagonal exclusion
    pair_mask_factorized = factorized_outer & diagonal_mask  # [B, N_token, N_R]

    # ========== Compare ==========
    # The test: are they equal?
    assert torch.equal(pair_mask_factorized, pair_mask_serial), "Factorized mask does NOT equal serial pair_mask"


@pytest.mark.parametrize("multiplicity", [1, 2], ids=lambda x: f"multiplicity:{x}")
def test_pair_mask_factorized_cdist_lddt(multiplicity):
    """Test that cdist_lddt with factorized masks matches cdist + lddt_dist.

    This test verifies that:
    1. cdist_lddt forward output matches torch.cdist + lddt_dist
    2. Multiplicity is handled correctly for coordinates and masks
    """
    if not torch.cuda.is_available():
        pytest.skip("CUDA not available")

    device = torch.device("cuda")
    torch.manual_seed(42)
    dtype = torch.float32  # Use FP32 for testing

    # Test dimensions
    B = 2
    N_token = 32
    n_atoms_per_token_min, n_atoms_per_token_max = 1, 3
    N_atom = N_token * (n_atoms_per_token_min + n_atoms_per_token_max) // 2

    # Generate realistic features using random_features
    feats = random_features(
        size_batch=B,
        n_tokens=N_token,
        n_atoms=N_atom,
        n_msa=1,
        atom_counts_per_token_range=(n_atoms_per_token_min, n_atoms_per_token_max),
        device=device,
        float_value_range=(-0.1, 0.1),
        selected_keys=[
            "token_to_rep_atom",
            "r_set_to_rep_atom",
        ],
    )

    token_to_rep_atom = feats["token_to_rep_atom"].to(dtype=dtype)  # [B, N_token, N_atom]
    r_set_to_rep_atom = feats["r_set_to_rep_atom"].to(dtype=dtype)  # [B, N_R, N_atom]

    N_R = r_set_to_rep_atom.shape[1]
    N_atom_actual = token_to_rep_atom.shape[2]

    # Create random atom_mask [B*mult, N_atom] - boolean (with multiplicity)
    atom_mask = torch.randint(0, 2, (B * multiplicity, N_atom_actual), device=device, dtype=torch.bool)
    # Ensure at least half atoms are resolved
    atom_mask[:, : N_atom_actual // 2] = True

    # Create random coordinates for tokens (rows) and R-set (columns) with multiplicity
    # Use uniform [-10, 10] to ensure meaningful distance distribution around cutoff
    # Shape: [B*mult, N_token/N_R, 3]
    pred_coords_row = torch.empty(B * multiplicity, N_token, 3, device=device, dtype=dtype).uniform_(-10.0, 10.0)
    true_coords_row = torch.empty(B * multiplicity, N_token, 3, device=device, dtype=dtype).uniform_(-10.0, 10.0)
    pred_coords_col = torch.empty(B * multiplicity, N_R, 3, device=device, dtype=dtype).uniform_(-10.0, 10.0)
    true_coords_col = torch.empty(B * multiplicity, N_R, 3, device=device, dtype=dtype).uniform_(-10.0, 10.0)

    # Compute average distance to set cutoff at roughly the median distance
    # This ensures meaningful coverage of pairs both within and outside cutoff
    pred_d_test = torch.cdist(pred_coords_row, pred_coords_col)
    true_d_test = torch.cdist(true_coords_row, true_coords_col)
    avg_dist = (pred_d_test.mean().item() + true_d_test.mean().item()) / 2.0
    cutoff = avg_dist

    # ========== Serial approach: torch.cdist + lddt_dist with full pair_mask ==========
    # Expand token_to_rep_atom and r_set_to_rep_atom for multiplicity
    token_to_rep_atom_mult = token_to_rep_atom.repeat_interleave(multiplicity, dim=0)  # [B*mult, N_token, N_atom]
    r_set_to_rep_atom_mult = r_set_to_rep_atom.repeat_interleave(multiplicity, dim=0)  # [B*mult, N_R, N_atom]

    # Construct pair_mask with multiplicity
    atom_mask_float = atom_mask.to(dtype=dtype)
    pair_mask = atom_mask_float.unsqueeze(-1) * atom_mask_float.unsqueeze(-2)  # [B*mult, N_atom, N_atom]
    diag_mask = 1.0 - torch.eye(N_atom_actual, device=device, dtype=dtype)
    pair_mask = pair_mask * diag_mask[None, :, :]
    pair_mask = torch.einsum("bnm,bkm->bnk", pair_mask, r_set_to_rep_atom_mult)  # [B*mult, N_atom, N_R]
    pair_mask_serial = torch.bmm(token_to_rep_atom_mult, pair_mask)  # [B*mult, N_token, N_R]

    # Compute distances using torch.cdist
    pred_d_serial = torch.cdist(pred_coords_row, pred_coords_col)  # [B*mult, N_token, N_R]
    true_d_serial = torch.cdist(true_coords_row, true_coords_col)  # [B*mult, N_token, N_R]

    # Compute lddt using lddt_dist
    lddt_serial, mask_no_match_serial = lddt_dist(
        pred_d_serial, true_d_serial, pair_mask_serial, cutoff=cutoff, per_atom=True
    )

    # ========== Factorized approach: cdist_lddt ==========
    # Factorized masks with multiplicity [B*mult, N_token] and [B*mult, N_R]
    mask_row = torch.bmm(token_to_rep_atom_mult, atom_mask_float.unsqueeze(-1)).squeeze(-1)  # [B*mult, N_token]
    mask_col = torch.bmm(r_set_to_rep_atom_mult, atom_mask_float.unsqueeze(-1)).squeeze(-1)  # [B*mult, N_R]

    # Representative atom indices [B, N_token] and [B, N_R] (no multiplicity)
    rep_atom_token = token_to_rep_atom.argmax(dim=-1)  # [B, N_token]
    rep_atom_r_set = r_set_to_rep_atom.argmax(dim=-1)  # [B, N_R]

    lddt_cdist, mask_no_match_cdist = cdist_lddt(
        pred_coords_row=pred_coords_row,
        pred_coords_col=pred_coords_col,
        true_coords_row=true_coords_row,
        true_coords_col=true_coords_col,
        mask_row=mask_row,
        mask_col=mask_col,
        multiplicity=multiplicity,
        atom_indices_row=rep_atom_token,
        atom_indices_col=rep_atom_r_set,
        cutoff=cutoff,
        do_mask_diagonal=True,
        per_atom=True,
    )

    # ========== Compare forward pass ==========
    # cdist_lddt Triton kernel uses float32 internally for performance,
    # so convert serial reference to match for comparison
    torch.testing.assert_close(lddt_cdist, lddt_serial.to(lddt_cdist.dtype))
    torch.testing.assert_close(mask_no_match_cdist, mask_no_match_serial.to(mask_no_match_cdist.dtype))


def parallel_assert_lddt_resolved_token(
    rank: int,
    payload: tuple,
):
    """Parallel test function for DTensor lddt_resolved_token().

    This function runs on each rank in the distributed setup and verifies that
    the DTensor lddt_resolved_token implementation matches the serial reference.
    """
    (
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        pred_atom_coords_global_host,
        true_atom_coords_global_host,
        true_coords_resolved_mask_global_host,
        token_to_rep_atom_global_host,
        r_set_to_rep_atom_global_host,
        atom_to_token_global_host,
        mol_type_global_host,
        atom_counts_per_token_host,
        expected_target_lddt_host,
        expected_combined_mask_host,
        multiplicity,
        cutoff_value,
    ) = payload

    monkeypatch = pytest.MonkeyPatch()
    if env_per_rank is not None:
        for var_name, value in env_per_rank.items():
            if value == "<INPUT_RANK>":
                monkeypatch.setenv(var_name, f"{rank}")
                continue
            monkeypatch.setenv(var_name, value)

    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    manager = DistributedManager()

    device_mesh = manager.device_mesh_subgroups
    dtype = pred_atom_coords_global_host.dtype

    # Create TransposeComm for redistribute_transpose
    comm = TransposeComm(manager.group["cp"], manager.layout_subgroups["cp"])

    # --- Distribute atom features using distribute_atom_features utility ---
    size_batch = token_to_rep_atom_global_host.shape[0]
    inputs_atom = {
        "atom_counts_per_token": atom_counts_per_token_host.to(dtype=torch.int64),
        "token_to_rep_atom": token_to_rep_atom_global_host.to(dtype=dtype),
        "r_set_to_rep_atom": r_set_to_rep_atom_global_host.to(dtype=dtype),
        "atom_to_token": atom_to_token_global_host.to(dtype=dtype),
    }

    # Add per-multiplicity coordinates and masks
    # Unflatten [B*mult, N_atom, 3] -> [B, mult, N_atom, 3]
    pred_coords_unflat = pred_atom_coords_global_host.unflatten(0, (size_batch, multiplicity))
    true_coords_unflat = true_atom_coords_global_host.unflatten(0, (size_batch, multiplicity))
    resolved_mask_unflat = true_coords_resolved_mask_global_host.unflatten(0, (size_batch, multiplicity))

    for i_mul in range(multiplicity):
        inputs_atom[f"pred_atom_coords_{i_mul}"] = pred_coords_unflat[:, i_mul].to(dtype=dtype)
        inputs_atom[f"true_atom_coords_{i_mul}"] = true_coords_unflat[:, i_mul].to(dtype=dtype)
        inputs_atom[f"true_coords_resolved_mask_{i_mul}"] = resolved_mask_unflat[:, i_mul].to(dtype=dtype)

    # Define placements for CP submesh and full mesh
    # Note: distribute_atom_features expects Replicate() for second dim of atom features
    placements_cp = {
        "atom_counts_per_token": (Shard(0), Replicate()),
        "token_to_rep_atom": (Shard(0), Replicate()),
        "r_set_to_rep_atom": (Shard(0), Replicate()),
        "atom_to_token": (Shard(0), Replicate()),
    }
    placements_dp_cp = {
        "token_to_rep_atom": (Shard(0), Shard(1), Replicate()),
        "r_set_to_rep_atom": (Shard(0), Shard(1), Replicate()),
        "atom_to_token": (Shard(0), Shard(1), Replicate()),
    }
    for i_mul in range(multiplicity):
        placements_cp[f"pred_atom_coords_{i_mul}"] = (Shard(0), Replicate())
        placements_cp[f"true_atom_coords_{i_mul}"] = (Shard(0), Replicate())
        placements_cp[f"true_coords_resolved_mask_{i_mul}"] = (Shard(0), Replicate())
        placements_dp_cp[f"pred_atom_coords_{i_mul}"] = (Shard(0), Shard(1), Replicate())
        placements_dp_cp[f"true_atom_coords_{i_mul}"] = (Shard(0), Shard(1), Replicate())
        placements_dp_cp[f"true_coords_resolved_mask_{i_mul}"] = (Shard(0), Shard(1), Replicate())

    # Distribute atom features with intersperse padding
    feats_atom = distribute_atom_features(
        inputs_atom,
        placements_cp,
        placements_dp_cp,
        device_mesh,
        manager.group["cp"],
        multiplicities={
            "pred_atom_coords": multiplicity,
            "true_atom_coords": multiplicity,
            "true_coords_resolved_mask": multiplicity,
        },
    )

    # --- Distribute token features (mol_type) using distribute_features ---
    # mol_type is a token feature [B, N_token], not an atom feature
    # Only rank 0 in the world group provides the features, others pass None
    if manager.group_rank["world"] == 0:
        token_features = {
            "mol_type": mol_type_global_host.to(device=manager.device, dtype=torch.int64),
        }
    else:
        token_features = None
    token_placements = {
        "mol_type": (Shard(0), Shard(1), Replicate()),
    }
    token_feats_dtensor = distribute_features(
        token_features,
        token_placements,
        manager.group["world"],
        manager.group_ranks["world"][0],
        device_mesh,
    )

    # Extract distributed tensors
    pred_atom_coords_dtensor = feats_atom["pred_atom_coords"]
    true_atom_coords_dtensor = feats_atom["true_atom_coords"]
    true_coords_resolved_mask_dtensor = feats_atom["true_coords_resolved_mask"]

    # Create feature dictionary
    feats_dtensor = {
        "token_to_rep_atom": feats_atom["token_to_rep_atom"],
        "r_set_to_rep_atom": feats_atom["r_set_to_rep_atom"],
        "atom_to_token": feats_atom["atom_to_token"],
        "mol_type": token_feats_dtensor["mol_type"],
    }

    # Call distributed lddt_resolved_token()
    # Returns (target_lddt, combined_mask) where combined_mask = token_resolved_mask * mask_no_match
    target_lddt_dtensor, combined_mask_dtensor = lddt_resolved_token(
        pred_atom_coords_dtensor,
        true_atom_coords_dtensor,
        true_coords_resolved_mask_dtensor,
        feats_dtensor,
        comm,
        multiplicity=multiplicity,
        cutoff=cutoff_value,
    )

    # Verify against serial reference
    target_lddt_global = target_lddt_dtensor.full_tensor().cpu()
    # Match dtype of DTensor output (preserves input coordinate dtype)
    expected_target_lddt_global = expected_target_lddt_host.to(dtype=target_lddt_global.dtype)

    torch.testing.assert_close(
        target_lddt_global,
        expected_target_lddt_global,
    )

    combined_mask_global = combined_mask_dtensor.full_tensor().cpu()
    # Match dtype of DTensor output (inherits from coordinate dtype)
    expected_combined_mask_global = expected_combined_mask_host.to(dtype=combined_mask_global.dtype)

    torch.testing.assert_close(
        combined_mask_global,
        expected_combined_mask_global,
    )

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env",
    [
        ((1, (2, 2)), True, "cuda", "ENV"),
        ((2, (2, 2)), True, "cuda", "ENV"),
    ],
    indirect=("setup_env",),
    ids=lambda x: f"dp:{x[0][0]}, cp:{x[0][1]}, device_type:{x[2]}",
)
@pytest.mark.parametrize("multiplicity", [1, 2], ids=lambda x: f"multiplicity:{x}")
def test_lddt_resolved_token(setup_env, multiplicity):
    """Test DTensor lddt_resolved_token() implementation against serial reference.

    This test verifies that the distributed lddt_resolved_token computation matches the
    serial implementation by:
    1. Generating realistic features using random_features()
    2. Computing serial reference using exact code from plddt_loss
    3. Sharding inputs using distribute_atom_features()
    4. Calling distributed lddt_resolved_token() function
    5. Verifying output matches serial reference
    """
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() < {world_size}")

    seed = 42
    torch.manual_seed(seed)

    # Test dimensions
    dp_size = grid_group_sizes["dp"]
    cp_size = grid_group_sizes["cp"][0]

    # Batch size must equal DP size for per-sample processing in pad_and_scatter
    B = dp_size
    n_tokens_per_shard = 16
    N_token = n_tokens_per_shard * cp_size
    # Atoms: use variable atoms per token (1-3 atoms per token)
    n_atoms_per_token_min, n_atoms_per_token_max = 1, 3
    # Estimate total atoms: avg 2 atoms/token * N_token, rounded to be divisible by cp_size
    avg_atoms_per_token = (n_atoms_per_token_min + n_atoms_per_token_max) / 2
    N_atom = int(avg_atoms_per_token * N_token)
    # Make N_atom divisible by cp_size for even sharding
    N_atom = ((N_atom + cp_size - 1) // cp_size) * cp_size
    dtype = torch.float32  # Use FP32 for testing

    # Use random_features to generate features with proper block-diagonal structure
    feats = random_features(
        size_batch=B,
        n_tokens=N_token,
        n_atoms=N_atom,
        n_msa=1,
        atom_counts_per_token_range=(n_atoms_per_token_min, n_atoms_per_token_max),
        device=torch.device(device_type),
        float_value_range=(-0.1, 0.1),
        selected_keys=[
            "token_to_rep_atom",
            "r_set_to_rep_atom",
            "atom_to_token",
            "mol_type",
            "atom_counts_per_token",
        ],
    )

    token_to_rep_atom_global = feats["token_to_rep_atom"].to(dtype=dtype)  # [B, N_token, N_atom]
    r_set_to_rep_atom_global = feats["r_set_to_rep_atom"].to(dtype=dtype)  # [B, N_R, N_atom]
    atom_to_token_global = feats["atom_to_token"].to(dtype=dtype)  # [B, N_atom, N_token]
    mol_type_global = feats["mol_type"]  # [B, N_token]
    atom_counts_per_token = feats["atom_counts_per_token"]

    N_atom_actual = token_to_rep_atom_global.shape[2]

    # Generate coordinates using uniform distribution
    # Use [-10, 10] range to ensure meaningful distance distribution around cutoff
    # The cutoff will be computed as average distance to ensure coverage above and below
    pred_atom_coords_global = torch.empty(B * multiplicity, N_atom_actual, 3, device=device_type, dtype=dtype).uniform_(
        -10.0, 10.0
    )
    true_atom_coords_global = torch.empty(B * multiplicity, N_atom_actual, 3, device=device_type, dtype=dtype).uniform_(
        -10.0, 10.0
    )

    # Create true_coords_resolved_mask: [B*mult, N_atom]
    true_coords_resolved_mask_global = torch.randint(
        0, 2, (B * multiplicity, N_atom_actual), device=device_type, dtype=dtype
    )
    # Ensure at least half atoms are resolved for meaningful test
    true_coords_resolved_mask_global[:, : N_atom_actual // 2] = 1.0

    # === SERIAL REFERENCE (exact copy from plddt_loss lines 178-231) ===
    atom_mask = true_coords_resolved_mask_global

    R_set_to_rep_atom = r_set_to_rep_atom_global.repeat_interleave(multiplicity, 0).to(dtype=dtype)

    token_type = mol_type_global.repeat_interleave(multiplicity, 0)
    is_nucleotide_token = (token_type == const.chain_type_ids["DNA"]).to(dtype=dtype) + (
        token_type == const.chain_type_ids["RNA"]
    ).to(dtype=dtype)

    atom_to_token = atom_to_token_global.to(dtype=dtype).repeat_interleave(multiplicity, 0)
    token_to_rep_atom = token_to_rep_atom_global.to(dtype=dtype).repeat_interleave(multiplicity, 0)

    true_token_coords = torch.bmm(token_to_rep_atom, true_atom_coords_global)
    pred_token_coords = torch.bmm(token_to_rep_atom, pred_atom_coords_global)

    true_d = torch.cdist(
        true_token_coords,
        torch.bmm(R_set_to_rep_atom, true_atom_coords_global),
    )
    pred_d = torch.cdist(
        pred_token_coords,
        torch.bmm(R_set_to_rep_atom, pred_atom_coords_global),
    )

    # pair_mask construction
    pair_mask = atom_mask.unsqueeze(-1) * atom_mask.unsqueeze(-2)  # [B, N_atom, N_atom]
    pair_mask = pair_mask * (1 - torch.eye(pair_mask.shape[1], device=pair_mask.device))[None, :, :]
    pair_mask = torch.einsum("bnm,bkm->bnk", pair_mask, R_set_to_rep_atom)  # [B, N_atom, N_R]
    pair_mask = torch.bmm(token_to_rep_atom, pair_mask)  # [B, N_token, N_R]

    is_nucleotide_R_element = torch.bmm(
        R_set_to_rep_atom, torch.bmm(atom_to_token, is_nucleotide_token.unsqueeze(-1))
    ).squeeze(-1)

    # Compute average inter-atom distance for cutoff to ensure coverage below and above cutoff
    # This gives better numerical stability in the lDDT computation
    avg_pred_dist = pred_d.mean().item()
    avg_true_dist = true_d.mean().item()
    cutoff_value = (avg_pred_dist + avg_true_dist) / 2.0

    cutoff = cutoff_value + cutoff_value * is_nucleotide_R_element.reshape(B * multiplicity, 1, -1).repeat(
        1, true_d.shape[1], 1
    )

    # lddt_dist (per_atom=True)
    expected_target_lddt, expected_mask_no_match = lddt_dist(pred_d, true_d, pair_mask, cutoff, per_atom=True)

    # Compute token_resolved_mask (whether each token has a resolved representative atom)
    # This matches the computation in lddt_resolved_token
    token_resolved_mask = torch.bmm(
        token_to_rep_atom, atom_mask.unsqueeze(-1).to(dtype=token_to_rep_atom.dtype)
    ).squeeze(-1)
    expected_combined_mask = token_resolved_mask * expected_mask_no_match

    # Verify that lddt_dist does not support gradient computation.
    # The lDDT metric uses step functions (thresholding at 0.5, 1, 2, 4 Å) which break the
    # autograd graph. Even with requires_grad=True on inputs, backward() raises RuntimeError.
    pred_d_grad = pred_d.detach().clone().requires_grad_(True)
    true_d_grad = true_d.detach().clone().requires_grad_(True)
    lddt_out, _ = lddt_dist(pred_d_grad, true_d_grad, pair_mask, cutoff, per_atom=True)
    with pytest.raises(RuntimeError, match="does not require grad and does not have a grad_fn"):
        lddt_out.sum().backward()

    # Prepare payload for parallel test
    payload = (
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        pred_atom_coords_global.clone().cpu(),
        true_atom_coords_global.clone().cpu(),
        true_coords_resolved_mask_global.clone().cpu(),
        token_to_rep_atom_global.clone().cpu(),
        r_set_to_rep_atom_global.clone().cpu(),
        atom_to_token_global.clone().cpu(),
        mol_type_global.clone().cpu(),
        atom_counts_per_token.clone().cpu(),
        expected_target_lddt.detach().clone().cpu(),
        expected_combined_mask.detach().clone().cpu(),
        multiplicity,
        cutoff_value,
    )

    # Launch parallel test
    spawn_multiprocessing(parallel_assert_lddt_resolved_token, world_size, payload)


def parallel_assert_plddt_loss(
    rank: int,
    payload: tuple,
):
    """Worker function that runs on each rank to test plddt_loss DTensor implementation.

    Uses the same setup pattern as parallel_assert_lddt_resolved_token.
    """
    (
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        pred_lddt_global_host,
        pred_atom_coords_global_host,
        true_atom_coords_global_host,
        true_coords_resolved_mask_global_host,
        token_to_rep_atom_global_host,
        r_set_to_rep_atom_global_host,
        atom_to_token_global_host,
        mol_type_global_host,
        atom_counts_per_token_host,
        expected_loss_host,
        expected_grad_pred_lddt_host,
        multiplicity,
    ) = payload

    # Setup environment variables for this rank (same as test_lddt_resolved_token)
    monkeypatch = pytest.MonkeyPatch()
    if env_per_rank is not None:
        for var_name, value in env_per_rank.items():
            if value == "<INPUT_RANK>":
                monkeypatch.setenv(var_name, f"{rank}")
                continue
            monkeypatch.setenv(var_name, value)

    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    manager = DistributedManager()

    device_mesh = manager.device_mesh_subgroups
    dtype = pred_atom_coords_global_host.dtype

    # Create TransposeComm for redistribute_transpose
    comm = TransposeComm(manager.group["cp"], manager.layout_subgroups["cp"])

    # --- Distribute atom features using distribute_atom_features utility ---
    # (same pattern as parallel_assert_lddt_resolved_token)
    size_batch = token_to_rep_atom_global_host.shape[0]
    inputs_atom = {
        "atom_counts_per_token": atom_counts_per_token_host.to(dtype=torch.int64),
        "token_to_rep_atom": token_to_rep_atom_global_host.to(dtype=dtype),
        "r_set_to_rep_atom": r_set_to_rep_atom_global_host.to(dtype=dtype),
        "atom_to_token": atom_to_token_global_host.to(dtype=dtype),
    }

    # Add per-multiplicity coordinates and masks
    pred_coords_unflat = pred_atom_coords_global_host.unflatten(0, (size_batch, multiplicity))
    true_coords_unflat = true_atom_coords_global_host.unflatten(0, (size_batch, multiplicity))
    resolved_mask_unflat = true_coords_resolved_mask_global_host.unflatten(0, (size_batch, multiplicity))

    for i_mul in range(multiplicity):
        inputs_atom[f"pred_atom_coords_{i_mul}"] = pred_coords_unflat[:, i_mul].to(dtype=dtype)
        inputs_atom[f"true_atom_coords_{i_mul}"] = true_coords_unflat[:, i_mul].to(dtype=dtype)
        inputs_atom[f"true_coords_resolved_mask_{i_mul}"] = resolved_mask_unflat[:, i_mul].to(dtype=dtype)

    # Define placements for CP submesh and full mesh
    placements_cp = {
        "atom_counts_per_token": (Shard(0), Replicate()),
        "token_to_rep_atom": (Shard(0), Replicate()),
        "r_set_to_rep_atom": (Shard(0), Replicate()),
        "atom_to_token": (Shard(0), Replicate()),
    }
    placements_dp_cp = {
        "token_to_rep_atom": (Shard(0), Shard(1), Replicate()),
        "r_set_to_rep_atom": (Shard(0), Shard(1), Replicate()),
        "atom_to_token": (Shard(0), Shard(1), Replicate()),
    }
    for i_mul in range(multiplicity):
        placements_cp[f"pred_atom_coords_{i_mul}"] = (Shard(0), Replicate())
        placements_cp[f"true_atom_coords_{i_mul}"] = (Shard(0), Replicate())
        placements_cp[f"true_coords_resolved_mask_{i_mul}"] = (Shard(0), Replicate())
        placements_dp_cp[f"pred_atom_coords_{i_mul}"] = (Shard(0), Shard(1), Replicate())
        placements_dp_cp[f"true_atom_coords_{i_mul}"] = (Shard(0), Shard(1), Replicate())
        placements_dp_cp[f"true_coords_resolved_mask_{i_mul}"] = (Shard(0), Shard(1), Replicate())

    # Distribute atom features with intersperse padding
    feats_atom = distribute_atom_features(
        inputs_atom,
        placements_cp,
        placements_dp_cp,
        device_mesh,
        manager.group["cp"],
        multiplicities={
            "pred_atom_coords": multiplicity,
            "true_atom_coords": multiplicity,
            "true_coords_resolved_mask": multiplicity,
        },
    )

    # --- Distribute token features (mol_type and pred_lddt) using distribute_features ---
    if manager.group_rank["world"] == 0:
        token_features = {
            "mol_type": mol_type_global_host.to(device=manager.device, dtype=torch.int64),
            "pred_lddt": pred_lddt_global_host.to(device=manager.device, dtype=torch.float32),
        }
    else:
        token_features = None
    token_placements = {
        "mol_type": (Shard(0), Shard(1), Replicate()),
        "pred_lddt": (Shard(0), Shard(1), Replicate()),
    }
    token_feats_dtensor = distribute_features(
        token_features,
        token_placements,
        manager.group["world"],
        manager.group_ranks["world"][0],
        device_mesh,
    )

    # Extract distributed tensors
    pred_atom_coords_dtensor = feats_atom["pred_atom_coords"]
    true_atom_coords_dtensor = feats_atom["true_atom_coords"]
    true_coords_resolved_mask_dtensor = feats_atom["true_coords_resolved_mask"]

    # Create feature dictionary
    feats_dtensor = {
        "token_to_rep_atom": feats_atom["token_to_rep_atom"],
        "r_set_to_rep_atom": feats_atom["r_set_to_rep_atom"],
        "atom_to_token": feats_atom["atom_to_token"],
        "mol_type": token_feats_dtensor["mol_type"],
    }

    # Get pred_lddt DTensor with gradient tracking
    pred_lddt_dtensor = token_feats_dtensor["pred_lddt"]
    pred_lddt_dtensor_grad = pred_lddt_dtensor.detach().requires_grad_(True)

    # Compute plddt_loss
    loss = plddt_loss(
        pred_lddt=pred_lddt_dtensor_grad,
        pred_atom_coords=pred_atom_coords_dtensor,
        true_atom_coords=true_atom_coords_dtensor,
        true_coords_resolved_mask=true_coords_resolved_mask_dtensor,
        feats=feats_dtensor,
        comm=comm,
        multiplicity=multiplicity,
    )

    # Verify loss value
    loss_local = loss.to_local()
    # Match dtype of DTensor output (may inherit from coordinate dtype)
    expected_loss = expected_loss_host.to(device=loss_local.device, dtype=loss_local.dtype)
    torch.testing.assert_close(loss_local, expected_loss)

    # Verify gradients
    loss_local.backward()
    grad_pred_lddt = pred_lddt_dtensor_grad.grad

    # Full gather the gradient to compare
    grad_pred_lddt_full = grad_pred_lddt.full_tensor()
    # Match dtype of DTensor gradient output
    expected_grad = expected_grad_pred_lddt_host.to(device=grad_pred_lddt_full.device, dtype=grad_pred_lddt_full.dtype)
    torch.testing.assert_close(grad_pred_lddt_full, expected_grad)

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env",
    [
        ((1, (2, 2)), True, "cuda", "ENV"),
        ((2, (2, 2)), True, "cuda", "ENV"),
    ],
    indirect=("setup_env",),
    ids=lambda x: f"dp:{x[0][0]}, cp:{x[0][1]}, device_type:{x[2]}",
)
@pytest.mark.parametrize("multiplicity", [1, 2], ids=lambda x: f"multiplicity:{x}")
def test_plddt_loss(setup_env, multiplicity: int):
    """Test that DTensor plddt_loss matches serial reference.

    This test verifies:
    1. Forward pass: DTensor plddt_loss matches serial plddt_loss
    2. Backward pass: Gradients w.r.t. pred_lddt match serial gradients

    The test uses the same feature generation as test_lddt_resolved_token.
    """
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() < {world_size}")

    torch.manual_seed(42)

    dp_size = grid_group_sizes["dp"]
    cp_size = grid_group_sizes["cp"][0] * grid_group_sizes["cp"][1]

    # Generate test data
    B = dp_size  # Batch size equals DP size
    N_token = 32
    N_atom = 140  # Large enough to accommodate token atom counts
    n_atoms_per_token_min = 1
    n_atoms_per_token_max = 4
    num_bins = 50  # Number of pLDDT bins
    dtype = torch.float32  # Use FP32 for testing

    # Make N_atom divisible by cp_size for even sharding
    N_atom = ((N_atom + cp_size - 1) // cp_size) * cp_size

    # Use random_features to generate features with proper block-diagonal structure
    feats = random_features(
        size_batch=B,
        n_tokens=N_token,
        n_atoms=N_atom,
        n_msa=1,
        atom_counts_per_token_range=(n_atoms_per_token_min, n_atoms_per_token_max),
        device=torch.device(device_type),
        float_value_range=(-0.1, 0.1),
        selected_keys=[
            "token_to_rep_atom",
            "r_set_to_rep_atom",
            "atom_to_token",
            "mol_type",
            "atom_counts_per_token",
        ],
    )

    token_to_rep_atom_global = feats["token_to_rep_atom"].to(dtype=dtype)  # [B, N_token, N_atom]
    r_set_to_rep_atom_global = feats["r_set_to_rep_atom"].to(dtype=dtype)  # [B, N_R, N_atom]
    atom_to_token_global = feats["atom_to_token"].to(dtype=dtype)  # [B, N_atom, N_token]
    mol_type_global = feats["mol_type"]  # [B, N_token]
    atom_counts_per_token = feats["atom_counts_per_token"]

    N_atom_actual = token_to_rep_atom_global.shape[2]

    # Generate coordinates using uniform distribution
    # Use [-10, 10] range so average pairwise distance is ~15 (matching cutoff in serial plddt_loss)
    # This ensures meaningful coverage of pairs both within and outside cutoff
    pred_atom_coords_global = torch.empty(B * multiplicity, N_atom_actual, 3, device=device_type, dtype=dtype).uniform_(
        -10.0, 10.0
    )
    true_atom_coords_global = torch.empty(B * multiplicity, N_atom_actual, 3, device=device_type, dtype=dtype).uniform_(
        -10.0, 10.0
    )

    # Create true_coords_resolved_mask: [B*mult, N_atom]
    true_coords_resolved_mask_global = torch.randint(
        0, 2, (B * multiplicity, N_atom_actual), device=device_type, dtype=dtype
    )
    # Ensure at least half atoms are resolved for meaningful test
    true_coords_resolved_mask_global[:, : N_atom_actual // 2] = 1.0

    # Generate pred_lddt: (B*mult, N_token, num_bins)
    pred_lddt_global = torch.randn(
        B * multiplicity, N_token, num_bins, device=device_type, dtype=torch.float32
    ).requires_grad_(True)

    # Compute serial reference (serial uses float32 internally)
    feats_serial = {
        "token_to_rep_atom": token_to_rep_atom_global.clone().float(),
        "r_set_to_rep_atom": r_set_to_rep_atom_global.clone().float(),
        "atom_to_token": atom_to_token_global.clone().float(),
        "mol_type": mol_type_global.clone(),
    }

    expected_loss, _rel_loss = serial_plddt_loss(
        pred_lddt=pred_lddt_global,
        pred_atom_coords=pred_atom_coords_global.float(),
        feats=feats_serial,
        true_atom_coords=true_atom_coords_global.float(),
        true_coords_resolved_mask=true_coords_resolved_mask_global.float(),
        token_level_confidence=True,
        multiplicity=multiplicity,
    )

    # Compute gradients for serial reference
    expected_loss.backward()
    expected_grad_pred_lddt = pred_lddt_global.grad.clone()

    # Verify that serial reference produces gradients (the gradient flows through log_softmax)
    assert expected_grad_pred_lddt is not None, "Serial plddt_loss should produce gradients"
    assert not torch.allclose(
        expected_grad_pred_lddt, torch.zeros_like(expected_grad_pred_lddt)
    ), "Gradients should be non-zero"

    # Prepare payload for parallel test
    payload = (
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        pred_lddt_global.detach().clone().cpu(),
        pred_atom_coords_global.clone().cpu(),
        true_atom_coords_global.clone().cpu(),
        true_coords_resolved_mask_global.clone().cpu(),
        token_to_rep_atom_global.clone().cpu(),
        r_set_to_rep_atom_global.clone().cpu(),
        atom_to_token_global.clone().cpu(),
        mol_type_global.clone().cpu(),
        atom_counts_per_token.clone().cpu(),
        expected_loss.detach().clone().cpu(),
        expected_grad_pred_lddt.clone().cpu(),
        multiplicity,
    )

    # Launch parallel test
    spawn_multiprocessing(parallel_assert_plddt_loss, world_size, payload)


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
