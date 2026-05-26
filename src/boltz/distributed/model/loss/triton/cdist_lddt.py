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

import torch
import triton
import triton.language as tl


@triton.jit
def _cdist_lddt_kernel(
    # Pointers to inputs
    pred_coords_row,
    pred_coords_col,
    true_coords_row,
    true_coords_col,
    mask_row,
    mask_col,
    atom_indices_row,  # Can be None (use implicit arange)
    atom_indices_col,  # Can be None (use implicit arange)
    cutoff_col_ptr,  # Can be None (use scalar cutoff)
    out_num,
    out_denom,
    # Shape and strides for pred_coords_row [B_mul, N_row, 3]
    B_mul_pred_coords_row,
    N_pred_coords_row,
    stride_pred_coords_row_b,
    stride_pred_coords_row_n,
    stride_pred_coords_row_d,
    # Shape and strides for pred_coords_col [B_mul, N_col, 3]
    B_mul_pred_coords_col,
    N_pred_coords_col,
    stride_pred_coords_col_b,
    stride_pred_coords_col_n,
    stride_pred_coords_col_d,
    # Shape and strides for true_coords_row [B_mul, N_row, 3]
    B_mul_true_coords_row,
    N_true_coords_row,
    stride_true_coords_row_b,
    stride_true_coords_row_n,
    stride_true_coords_row_d,
    # Shape and strides for true_coords_col [B_mul, N_col, 3]
    B_mul_true_coords_col,
    N_true_coords_col,
    stride_true_coords_col_b,
    stride_true_coords_col_n,
    stride_true_coords_col_d,
    # Shape and strides for mask_row [B, N_row]
    B_mask_row,
    N_mask_row,
    stride_mask_row_b,
    stride_mask_row_n,
    # Shape and strides for mask_col [B, N_col]
    B_mask_col,
    N_mask_col,
    stride_mask_col_b,
    stride_mask_col_n,
    # Shape and strides for atom_indices_row [B, N_row] (only used if USE_EXPLICIT_INDICES_ROW)
    B_atom_indices_row,
    N_atom_indices_row,
    stride_atom_indices_row_b,
    stride_atom_indices_row_n,
    # Shape and strides for atom_indices_col [B, N_col] (only used if USE_EXPLICIT_INDICES_COL)
    B_atom_indices_col,
    N_atom_indices_col,
    stride_atom_indices_col_b,
    stride_atom_indices_col_n,
    # Shape and strides for cutoff_col [B, N_col] (only used if USE_CUTOFF_COL)
    B_cutoff_col,
    N_cutoff_col,
    stride_cutoff_col_b,
    stride_cutoff_col_n,
    # Shape and strides for out_num [B_mul] or [B_mul, N_row] if PER_ATOM
    B_mul_out_num,
    stride_out_num_b,
    stride_out_num_n,  # Only used if PER_ATOM
    # Shape and strides for out_denom [B_mul] or [B_mul, N_row] if PER_ATOM
    B_mul_out_denom,
    stride_out_denom_b,
    stride_out_denom_n,  # Only used if PER_ATOM
    # Constants
    cutoff,
    # Block sizes
    BLOCK_M: tl.constexpr,
    BLOCK_N: tl.constexpr,
    # Coordinate dimension size (e.g., 3 for 3D)
    SIZE_DIM_D: tl.constexpr,
    # Flags for implicit arange (avoids creating torch.arange tensors)
    USE_EXPLICIT_INDICES_ROW: tl.constexpr,  # If False, use offs_m as indices
    USE_EXPLICIT_INDICES_COL: tl.constexpr,  # If False, use offs_n as indices
    # Flag for per-column cutoff
    USE_CUTOFF_COL: tl.constexpr,  # If True, use per-column cutoff values instead of scalar
    # Flag for diagonal masking
    DO_MASK_DIAGONAL: tl.constexpr,  # If True, exclude self-pairs where indices match
    # Flag for per-atom output
    PER_ATOM: tl.constexpr,  # If True, output per-row scores; otherwise total score
    # Memory layout orders for make_block_ptr (computed from argsort of strides)
    # For 3D coord blocks [1, N, D], order is a tuple of 3 indices
    ORDER_PRED_COORDS_ROW_0: tl.constexpr,
    ORDER_PRED_COORDS_ROW_1: tl.constexpr,
    ORDER_PRED_COORDS_ROW_2: tl.constexpr,
    ORDER_PRED_COORDS_COL_0: tl.constexpr,
    ORDER_PRED_COORDS_COL_1: tl.constexpr,
    ORDER_PRED_COORDS_COL_2: tl.constexpr,
    ORDER_TRUE_COORDS_ROW_0: tl.constexpr,
    ORDER_TRUE_COORDS_ROW_1: tl.constexpr,
    ORDER_TRUE_COORDS_ROW_2: tl.constexpr,
    ORDER_TRUE_COORDS_COL_0: tl.constexpr,
    ORDER_TRUE_COORDS_COL_1: tl.constexpr,
    ORDER_TRUE_COORDS_COL_2: tl.constexpr,
    # For 2D mask blocks [1, N], order is a tuple of 2 indices
    ORDER_MASK_ROW_0: tl.constexpr,
    ORDER_MASK_ROW_1: tl.constexpr,
    ORDER_MASK_COL_0: tl.constexpr,
    ORDER_MASK_COL_1: tl.constexpr,
    # For 2D atom_indices blocks [1, N], order is a tuple of 2 indices (only used if USE_EXPLICIT_INDICES_*)
    ORDER_ATOM_INDICES_ROW_0: tl.constexpr,
    ORDER_ATOM_INDICES_ROW_1: tl.constexpr,
    ORDER_ATOM_INDICES_COL_0: tl.constexpr,
    ORDER_ATOM_INDICES_COL_1: tl.constexpr,
    # For 2D cutoff_col blocks [1, N], order is a tuple of 2 indices (only used if USE_CUTOFF_COL)
    ORDER_CUTOFF_COL_0: tl.constexpr,
    ORDER_CUTOFF_COL_1: tl.constexpr,
    # Flags for mask multiplicity (whether masks have B*mul or B batch dimension)
    MASK_ROW_HAS_MUL: tl.constexpr,  # If True, mask_row has [B*mul, N] shape
    MASK_COL_HAS_MUL: tl.constexpr,  # If True, mask_col has [B*mul, N] shape
    # Explicit multiplicity factor
    MULTIPLICITY: tl.constexpr,
):
    # Validate coordinate dimension at compile time
    tl.static_assert(SIZE_DIM_D == 3, "SIZE_DIM_D must be 3 (3D coordinates)")

    # Next power of 2 for coordinate dimension (required by tl.arange)
    # Since SIZE_DIM_D is always 3, BLOCK_D is 4 (next power of 2)
    BLOCK_D: tl.constexpr = 4

    # 1. Grid Identification
    pid_batch = tl.program_id(0)  # Handles flattened batch * multiplicity
    pid_m = tl.program_id(1)  # Row block index
    pid_n = tl.program_id(2)  # Col block index

    # 2. Multiplicity & Broadcasting
    # Use explicit multiplicity parameter
    # Determine which sample in the original [B] batch this corresponds to
    batch_idx = pid_batch // MULTIPLICITY
    # Batch indices for masks (depends on whether masks have multiplicity)
    batch_idx_mask_row = pid_batch if MASK_ROW_HAS_MUL else batch_idx
    batch_idx_mask_col = pid_batch if MASK_COL_HAS_MUL else batch_idx

    # 3. Create block pointers using make_block_ptr (full-dimensional to capture stride-1 axis)
    # For coordinate tensors: 3D block [1, BLOCK_M, BLOCK_D] from [B_mul, N_row, 3]
    # Including batch dimension ensures order tuple captures all strides including if batch has stride=1
    pred_row_block_ptr = tl.make_block_ptr(
        base=pred_coords_row,
        shape=(B_mul_pred_coords_row, N_pred_coords_row, SIZE_DIM_D),
        strides=(stride_pred_coords_row_b, stride_pred_coords_row_n, stride_pred_coords_row_d),
        offsets=(pid_batch, pid_m * BLOCK_M, 0),
        block_shape=(1, BLOCK_M, BLOCK_D),
        order=(ORDER_PRED_COORDS_ROW_0, ORDER_PRED_COORDS_ROW_1, ORDER_PRED_COORDS_ROW_2),
    )
    pred_col_block_ptr = tl.make_block_ptr(
        base=pred_coords_col,
        shape=(B_mul_pred_coords_col, N_pred_coords_col, SIZE_DIM_D),
        strides=(stride_pred_coords_col_b, stride_pred_coords_col_n, stride_pred_coords_col_d),
        offsets=(pid_batch, pid_n * BLOCK_N, 0),
        block_shape=(1, BLOCK_N, BLOCK_D),
        order=(ORDER_PRED_COORDS_COL_0, ORDER_PRED_COORDS_COL_1, ORDER_PRED_COORDS_COL_2),
    )
    true_row_block_ptr = tl.make_block_ptr(
        base=true_coords_row,
        shape=(B_mul_true_coords_row, N_true_coords_row, SIZE_DIM_D),
        strides=(stride_true_coords_row_b, stride_true_coords_row_n, stride_true_coords_row_d),
        offsets=(pid_batch, pid_m * BLOCK_M, 0),
        block_shape=(1, BLOCK_M, BLOCK_D),
        order=(ORDER_TRUE_COORDS_ROW_0, ORDER_TRUE_COORDS_ROW_1, ORDER_TRUE_COORDS_ROW_2),
    )
    true_col_block_ptr = tl.make_block_ptr(
        base=true_coords_col,
        shape=(B_mul_true_coords_col, N_true_coords_col, SIZE_DIM_D),
        strides=(stride_true_coords_col_b, stride_true_coords_col_n, stride_true_coords_col_d),
        offsets=(pid_batch, pid_n * BLOCK_N, 0),
        block_shape=(1, BLOCK_N, BLOCK_D),
        order=(ORDER_TRUE_COORDS_COL_0, ORDER_TRUE_COORDS_COL_1, ORDER_TRUE_COORDS_COL_2),
    )

    # For mask tensors: 2D block [1, BLOCK_M] from [B or B*mul, N_row]
    # Uses batch_idx_mask_row/col which is pid_batch if mask has multiplicity, else batch_idx
    mask_row_block_ptr = tl.make_block_ptr(
        base=mask_row,
        shape=(B_mask_row, N_mask_row),
        strides=(stride_mask_row_b, stride_mask_row_n),
        offsets=(batch_idx_mask_row, pid_m * BLOCK_M),
        block_shape=(1, BLOCK_M),
        order=(ORDER_MASK_ROW_0, ORDER_MASK_ROW_1),
    )
    mask_col_block_ptr = tl.make_block_ptr(
        base=mask_col,
        shape=(B_mask_col, N_mask_col),
        strides=(stride_mask_col_b, stride_mask_col_n),
        offsets=(batch_idx_mask_col, pid_n * BLOCK_N),
        block_shape=(1, BLOCK_N),
        order=(ORDER_MASK_COL_0, ORDER_MASK_COL_1),
    )

    # For atom_indices tensors: 2D block [1, BLOCK_M/N] from [B, N_row/col]
    # Uses batch_idx for broadcasting (same as masks)
    # Only create block pointers if explicit indices are used
    if USE_EXPLICIT_INDICES_ROW:
        atom_indices_row_block_ptr = tl.make_block_ptr(
            base=atom_indices_row,
            shape=(B_atom_indices_row, N_atom_indices_row),
            strides=(stride_atom_indices_row_b, stride_atom_indices_row_n),
            offsets=(batch_idx, pid_m * BLOCK_M),
            block_shape=(1, BLOCK_M),
            order=(ORDER_ATOM_INDICES_ROW_0, ORDER_ATOM_INDICES_ROW_1),
        )
    if USE_EXPLICIT_INDICES_COL:
        atom_indices_col_block_ptr = tl.make_block_ptr(
            base=atom_indices_col,
            shape=(B_atom_indices_col, N_atom_indices_col),
            strides=(stride_atom_indices_col_b, stride_atom_indices_col_n),
            offsets=(batch_idx, pid_n * BLOCK_N),
            block_shape=(1, BLOCK_N),
            order=(ORDER_ATOM_INDICES_COL_0, ORDER_ATOM_INDICES_COL_1),
        )

    # 5. Load Data using block pointers with boundary_check and padding_option
    # boundary_check specifies which dimensions to check for out-of-bounds
    # padding_option="zero" fills out-of-bounds values with 0
    # Only check dims 1 and 2 (N and D) for coords; dim 0 (batch) is always in bounds
    # Reshape to squeeze the batch dimension of size 1: [1, BLOCK_M, BLOCK_D] -> [BLOCK_M, BLOCK_D]
    pred_row = tl.reshape(
        tl.load(pred_row_block_ptr, boundary_check=(1, 2), padding_option="zero"),
        (BLOCK_M, BLOCK_D),
    )
    pred_col = tl.reshape(
        tl.load(pred_col_block_ptr, boundary_check=(1, 2), padding_option="zero"),
        (BLOCK_N, BLOCK_D),
    )
    true_row = tl.reshape(
        tl.load(true_row_block_ptr, boundary_check=(1, 2), padding_option="zero"),
        (BLOCK_M, BLOCK_D),
    )
    true_col = tl.reshape(
        tl.load(true_col_block_ptr, boundary_check=(1, 2), padding_option="zero"),
        (BLOCK_N, BLOCK_D),
    )

    # Only check dim 1 (N) for masks; dim 0 (batch) is always in bounds
    # Reshape to squeeze: [1, BLOCK_M] -> [BLOCK_M]
    m_row = tl.reshape(
        tl.load(mask_row_block_ptr, boundary_check=(1,), padding_option="zero"),
        (BLOCK_M,),
    ).to(tl.int1)
    m_col = tl.reshape(
        tl.load(mask_col_block_ptr, boundary_check=(1,), padding_option="zero"),
        (BLOCK_N,),
    ).to(tl.int1)

    # 5. Distance Computation (Pairwise)
    # dist = sqrt(sum((x[:, None] - y[None, :])^2))
    # Expand dims for broadcasting: [BLOCK_M, 1, 3] - [1, BLOCK_N, 3]
    delta_pred = pred_row[:, None, :] - pred_col[None, :, :]
    d_pred_sq = tl.sum(delta_pred * delta_pred, axis=2)
    d_pred = tl.sqrt(d_pred_sq)

    delta_true = true_row[:, None, :] - true_col[None, :, :]
    d_true_sq = tl.sum(delta_true * delta_true, axis=2)
    d_true = tl.sqrt(d_true_sq)

    # 6. Validity Masking
    # Base validity: both atoms resolved (out-of-bounds already zeroed by masked load)
    valid = m_row[:, None] & m_col[None, :]

    # Offsets for Row/Col dimensions (needed for diagonal masking and per_atom accumulation)
    offs_m = pid_m * BLOCK_M + tl.arange(0, BLOCK_M)
    offs_n = pid_n * BLOCK_N + tl.arange(0, BLOCK_N)
    boundary_mask_m = offs_m < N_pred_coords_row
    boundary_mask_n = offs_n < N_pred_coords_col

    # 7. Diagonal Masking (conditional on DO_MASK_DIAGONAL)
    if DO_MASK_DIAGONAL:
        # Get row indices: either from explicit tensor or use implicit arange (offs_m)
        # atom_indices use batch_idx (same as masks) for [B, N] broadcasting
        if USE_EXPLICIT_INDICES_ROW:
            # Load using block pointer with boundary_check and padding_option
            # padding_option=-1 to distinguish from valid indices (which are >= 0)
            idx_row = tl.reshape(
                tl.load(atom_indices_row_block_ptr, boundary_check=(1,), padding_option="zero"),
                (BLOCK_M,),
            )
            # Set out-of-bounds indices to -1 (won't match any valid column index)
            idx_row = tl.where(boundary_mask_m, idx_row, -1)
        else:
            # Use offs_m directly as indices (equivalent to torch.arange(N_row))
            idx_row = offs_m

        # Get col indices: either from explicit tensor or use implicit arange (offs_n)
        if USE_EXPLICIT_INDICES_COL:
            # Load using block pointer with boundary_check and padding_option
            # padding_option=-2 to distinguish from valid indices and row padding
            idx_col = tl.reshape(
                tl.load(atom_indices_col_block_ptr, boundary_check=(1,), padding_option="zero"),
                (BLOCK_N,),
            )
            # Set out-of-bounds indices to -2 (won't match any valid row index)
            idx_col = tl.where(boundary_mask_n, idx_col, -2)
        else:
            # Use offs_n directly as indices (equivalent to torch.arange(N_col))
            idx_col = offs_n

        # Apply diagonal mask (exclude self-pairs where indices match)
        is_diagonal = idx_row[:, None] == idx_col[None, :]
        valid = valid & (~is_diagonal)

    # 7. Scoring
    # Only score pairs where true distance < cutoff
    # If USE_CUTOFF_COL, use per-column cutoff values [B, N_col]; otherwise use scalar cutoff
    if USE_CUTOFF_COL:
        # Load per-column cutoff values for this block using block pointer
        # Uses batch_idx for broadcasting (B_mul -> B), same as masks
        cutoff_col_block_ptr = tl.make_block_ptr(
            base=cutoff_col_ptr,
            shape=(B_cutoff_col, N_cutoff_col),
            strides=(stride_cutoff_col_b, stride_cutoff_col_n),
            offsets=(batch_idx, pid_n * BLOCK_N),
            block_shape=(1, BLOCK_N),
            order=(ORDER_CUTOFF_COL_0, ORDER_CUTOFF_COL_1),
        )
        cutoff_vals = tl.reshape(
            tl.load(cutoff_col_block_ptr, boundary_check=(1,), padding_option="zero"),
            (BLOCK_N,),
        )
        # Broadcast [BLOCK_N] over [BLOCK_M, BLOCK_N]: cutoff_vals[None, :] is [1, BLOCK_N]
        within_cutoff = d_true < cutoff_vals[None, :]
    else:
        within_cutoff = d_true < cutoff
    active = valid & within_cutoff

    diff = tl.abs(d_pred - d_true)
    score_05 = (diff < 0.5).to(tl.float32)
    score_10 = (diff < 1.0).to(tl.float32)
    score_20 = (diff < 2.0).to(tl.float32)
    score_40 = (diff < 4.0).to(tl.float32)

    score_tile = 0.25 * (score_05 + score_10 + score_20 + score_40)

    # Zero out invalid entries
    score_tile = tl.where(active, score_tile, 0.0)
    denom_tile = tl.where(active, 1.0, 0.0)

    # 8. Accumulation
    if PER_ATOM:
        # Sum over columns (BLOCK_N dimension) for per-row scores
        num_per_row = tl.sum(score_tile, axis=1)  # [BLOCK_M]
        denom_per_row = tl.sum(denom_tile, axis=1)  # [BLOCK_M]

        # Atomic Add to global buffers [B_mul, N_row]
        out_ptr_num = out_num + pid_batch * stride_out_num_b + offs_m * stride_out_num_n
        out_ptr_denom = out_denom + pid_batch * stride_out_denom_b + offs_m * stride_out_denom_n

        tl.atomic_add(out_ptr_num, num_per_row, mask=boundary_mask_m, sem="relaxed")
        tl.atomic_add(out_ptr_denom, denom_per_row, mask=boundary_mask_m, sem="relaxed")
    else:
        # Sum over the entire tile
        num_sum = tl.sum(score_tile)
        denom_sum = tl.sum(denom_tile)

        # Atomic Add to global buffers [B_mul]
        tl.atomic_add(out_num + pid_batch * stride_out_num_b, num_sum, sem="relaxed")
        tl.atomic_add(out_denom + pid_batch * stride_out_denom_b, denom_sum, sem="relaxed")


def cdist_lddt(
    pred_coords_row,  # [B_mul, N_row, 3]  (B_mul = B * multiplicity)
    pred_coords_col,  # [B_mul, N_col, 3]
    true_coords_row,  # [B_mul, N_row, 3]
    true_coords_col,  # [B_mul, N_col, 3]
    mask_row,  # [B, N_row] or [B_mul, N_row]
    mask_col,  # [B, N_col] or [B_mul, N_col]
    multiplicity,  # Required: explicit multiplicity factor (B_mul = B * multiplicity)
    atom_indices_row=None,  # [B, N_row] (optional, for diagonal masking)
    atom_indices_col=None,  # [B, N_col] (optional, for diagonal masking)
    cutoff=15.0,
    cutoff_col=None,  # [B, N_col] (optional, per-column cutoff values per batch)
    eps=1e-10,
    do_mask_diagonal=True,  # If True, exclude self-pairs where indices match
    return_unnormalized_score=False,  # If True, return (out_num, out_denom) before normalization
    per_atom=False,  # If True, return per-row lDDT scores and mask_no_match
    return_denom=False,  # If True, also return the denominator (pair counts)
):
    """
    Computes lDDT score directly from coordinates without materializing distance matrices.
    Handles rectangular inputs and broadcasting of multiplicity.

    Note: Inputs must share a dtype. Computation runs in at least float32
    (via torch.promote_types with float32), and outputs are cast back to the input dtype.

    Parameters
    ----------
    pred_coords_row : torch.Tensor
        Predicted coordinates for row atoms, shape [B_mul, N_row, 3]
    pred_coords_col : torch.Tensor
        Predicted coordinates for column atoms, shape [B_mul, N_col, 3]
    true_coords_row : torch.Tensor
        True coordinates for row atoms, shape [B_mul, N_row, 3]
    true_coords_col : torch.Tensor
        True coordinates for column atoms, shape [B_mul, N_col, 3]
    mask_row : torch.Tensor
        Resolved mask for row atoms, shape [B, N_row] or [B*mul, N_row].
        If [B, N_row], broadcasts to B_mul. If [B*mul, N_row], used directly.
    mask_col : torch.Tensor
        Resolved mask for column atoms, shape [B, N_col] or [B*mul, N_col].
        If [B, N_col], broadcasts to B_mul. If [B*mul, N_col], used directly.
    multiplicity : int
        Required. Explicit multiplicity factor where B_mul = B * multiplicity.
        This determines the base batch size B for validating tensor shapes.
    atom_indices_row : torch.Tensor, optional
        Explicit indices for row atoms, shape [B, N_row]. If None, uses implicit arange.
        The batch dimension B matches mask_row (not B_mul), enabling per-sample indices.
    atom_indices_col : torch.Tensor, optional
        Explicit indices for column atoms, shape [B, N_col]. If None, uses implicit arange.
        The batch dimension B matches mask_col (not B_mul), enabling per-sample indices.
    cutoff : float
        Distance cutoff for lDDT scoring (default 15.0). Used as fallback if cutoff_col is None.
    cutoff_col : torch.Tensor, optional
        Per-column distance cutoff values, shape [B, N_col]. If provided, overrides scalar cutoff.
        The batch dimension B matches mask_col (not B_mul), enabling per-sample cutoffs.
        This enables nucleotide-dependent cutoffs (e.g., 15.0 for protein, 30.0 for nucleotides).
    eps : float
        Small epsilon for numerical stability (default 1e-10)
    do_mask_diagonal : bool
        If True (default), exclude self-pairs where atom indices match.
    return_unnormalized_score : bool
        If True, return raw (out_num, out_denom) before normalization.
        This is useful for distributed aggregation where partial sums from different
        shards need to be allreduced before computing the final normalized score.
        Can be combined with per_atom to control output shape.
        If False (default), return normalized lDDT score.
    per_atom : bool
        If True, return per-row lDDT scores [B_mul, N_row] and mask_no_match [B_mul, N_row].
        If False (default), return total lDDT score [B_mul].
        Can be combined with return_unnormalized_score.
    return_denom : bool
        If True, also return the denominator (pair counts) used for normalization.
        Not valid when return_unnormalized_score=True.

    Returns
    -------
    If return_unnormalized_score=True and per_atom=True:
        out_num : torch.Tensor
            Per-row unnormalized sum of scores, shape [B_mul, N_row].
        out_denom : torch.Tensor
            Per-row sum of valid pair counts, shape [B_mul, N_row].
        mask_no_match : torch.Tensor
            Mask indicating rows with valid pairs, shape [B_mul, N_row].
        Note: For distributed allreduce on out_num and out_denom, then normalize manually:
            norm = 1.0 / (eps + allreduced_denom)
            score = norm * (eps + allreduced_num)
            # mask_no_match should be allreduced with logical OR across shards
    If return_unnormalized_score=True and per_atom=False:
        out_num : torch.Tensor
            Unnormalized sum of scores, shape [B_mul].
        out_denom : torch.Tensor
            Sum of valid pair counts, shape [B_mul].
    If return_unnormalized_score=False and per_atom=False:
        score : torch.Tensor
            lDDT score per batch element, shape [B_mul]
        denom : torch.Tensor, optional
            Pair counts per batch element, shape [B_mul]
    If return_unnormalized_score=False and per_atom=True:
        score : torch.Tensor
            Per-row lDDT score, shape [B_mul, N_row]
        mask_no_match : torch.Tensor
            Boolean mask indicating rows with valid pairs, shape [B_mul, N_row]
        denom : torch.Tensor, optional
            Pair counts per row, shape [B_mul, N_row]
    """
    # Extract reference shapes from primary tensors
    B_mul, N_row, dim_d = pred_coords_row.shape
    _, N_col, _ = pred_coords_col.shape
    B_or_B_mul_mask_row, _ = mask_row.shape
    B_or_B_mul_mask_col, _ = mask_col.shape

    # Validate coordinate dimension (must be 3D coordinates)
    if dim_d != 3:
        raise ValueError(f"Coordinate dimension must be 3, got {dim_d}")

    if return_unnormalized_score and return_denom:
        raise ValueError("return_denom is not valid when return_unnormalized_score=True")

    # Compute B from explicit multiplicity
    if B_mul % multiplicity != 0:
        raise ValueError(f"Coordinate batch dimension ({B_mul}) must be divisible by multiplicity ({multiplicity})")
    B = B_mul // multiplicity

    # Check if masks have multiplicity
    mask_row_has_mul = B_or_B_mul_mask_row == B_mul
    mask_col_has_mul = B_or_B_mul_mask_col == B_mul

    # Validate mask shapes: must be either (B, N) or (B_mul, N)
    if B_or_B_mul_mask_row != B and B_or_B_mul_mask_row != B_mul:
        raise ValueError(f"mask_row batch dimension ({B_or_B_mul_mask_row}) must be either B ({B}) or B_mul ({B_mul})")
    if B_or_B_mul_mask_col != B and B_or_B_mul_mask_col != B_mul:
        raise ValueError(f"mask_col batch dimension ({B_or_B_mul_mask_col}) must be either B ({B}) or B_mul ({B_mul})")

    # Note: return_unnormalized_score and per_atom are orthogonal options:
    # - per_atom controls output shape: [B_mul, N_row] vs [B_mul]
    # - return_unnormalized_score controls whether to return raw (out_num, out_denom) vs normalized score

    # Expected shapes for all tensors
    coord_row_shape = (B_mul, N_row, dim_d)
    coord_col_shape = (B_mul, N_col, dim_d)
    idx_row_shape = (B, N_row)
    idx_col_shape = (B, N_col)

    # Validate coordinate tensors
    if pred_coords_row.shape != coord_row_shape:
        raise ValueError(f"pred_coords_row shape {tuple(pred_coords_row.shape)} must be {coord_row_shape}")
    if pred_coords_col.shape != coord_col_shape:
        raise ValueError(f"pred_coords_col shape {tuple(pred_coords_col.shape)} must be {coord_col_shape}")
    if true_coords_row.shape != coord_row_shape:
        raise ValueError(f"true_coords_row shape {tuple(true_coords_row.shape)} must be {coord_row_shape}")
    if true_coords_col.shape != coord_col_shape:
        raise ValueError(f"true_coords_col shape {tuple(true_coords_col.shape)} must be {coord_col_shape}")

    # Validate mask N dimensions
    if mask_row.shape[1] != N_row:
        raise ValueError(f"mask_row N dimension ({mask_row.shape[1]}) must be {N_row}")
    if mask_col.shape[1] != N_col:
        raise ValueError(f"mask_col N dimension ({mask_col.shape[1]}) must be {N_col}")

    # Validate optional index tensors
    if atom_indices_row is not None and atom_indices_row.shape != idx_row_shape:
        raise ValueError(f"atom_indices_row shape {tuple(atom_indices_row.shape)} must be {idx_row_shape}")
    if atom_indices_col is not None and atom_indices_col.shape != idx_col_shape:
        raise ValueError(f"atom_indices_col shape {tuple(atom_indices_col.shape)} must be {idx_col_shape}")

    # Validate optional cutoff_col tensor
    cutoff_col_shape = (B, N_col)
    if cutoff_col is not None and cutoff_col.shape != cutoff_col_shape:
        raise ValueError(f"cutoff_col shape {tuple(cutoff_col.shape)} must be {cutoff_col_shape}")

    device = pred_coords_row.device

    if (
        pred_coords_row.dtype != pred_coords_col.dtype
        or pred_coords_row.dtype != true_coords_row.dtype
        or pred_coords_row.dtype != true_coords_col.dtype
    ):
        raise ValueError(
            "pred/true coords dtypes must match: "
            f"row={pred_coords_row.dtype}, col={pred_coords_col.dtype}, "
            f"true_row={true_coords_row.dtype}, true_col={true_coords_col.dtype}"
        )
    input_dtype = pred_coords_row.dtype
    compute_dtype = torch.promote_types(input_dtype, torch.float32)

    if pred_coords_row.dtype != compute_dtype:
        pred_coords_row = pred_coords_row.to(compute_dtype)
    if pred_coords_col.dtype != compute_dtype:
        pred_coords_col = pred_coords_col.to(compute_dtype)
    if true_coords_row.dtype != compute_dtype:
        true_coords_row = true_coords_row.to(compute_dtype)
    if true_coords_col.dtype != compute_dtype:
        true_coords_col = true_coords_col.to(compute_dtype)

    # Output buffers
    if per_atom:
        # Per-row outputs: [B_mul, N_row]
        out_num = torch.zeros(B_mul, N_row, device=device, dtype=compute_dtype)
        out_denom = torch.zeros(B_mul, N_row, device=device, dtype=compute_dtype)
    else:
        # Total outputs: [B_mul]
        out_num = torch.zeros(B_mul, device=device, dtype=compute_dtype)
        out_denom = torch.zeros(B_mul, device=device, dtype=compute_dtype)

    # Block sizes
    BLOCK_M = 32
    BLOCK_N = 32

    grid = (B_mul, triton.cdiv(N_row, BLOCK_M), triton.cdiv(N_col, BLOCK_N))

    # Determine whether to use explicit indices or implicit arange
    use_explicit_indices_row = atom_indices_row is not None
    use_explicit_indices_col = atom_indices_col is not None

    # Determine whether to use per-column cutoff
    use_cutoff_col = cutoff_col is not None

    # Compute memory layout order for make_block_ptr
    # order = argsort of strides (ascending), giving fastest-varying dim first
    order_pred_coords_row = tuple(torch.tensor(pred_coords_row.stride()).argsort().tolist())
    order_pred_coords_col = tuple(torch.tensor(pred_coords_col.stride()).argsort().tolist())
    order_true_coords_row = tuple(torch.tensor(true_coords_row.stride()).argsort().tolist())
    order_true_coords_col = tuple(torch.tensor(true_coords_col.stride()).argsort().tolist())
    order_mask_row = tuple(torch.tensor(mask_row.stride()).argsort().tolist())
    order_mask_col = tuple(torch.tensor(mask_col.stride()).argsort().tolist())
    # Compute order for atom_indices if they are provided
    order_atom_indices_row = (
        tuple(torch.tensor(atom_indices_row.stride()).argsort().tolist())
        if use_explicit_indices_row
        else (0, 1)  # Default order, not used when indices are None
    )
    order_atom_indices_col = (
        tuple(torch.tensor(atom_indices_col.stride()).argsort().tolist())
        if use_explicit_indices_col
        else (0, 1)  # Default order, not used when indices are None
    )
    # Compute order for cutoff_col if provided
    order_cutoff_col = (
        tuple(torch.tensor(cutoff_col.stride()).argsort().tolist())
        if use_cutoff_col
        else (0, 1)  # Default order, not used when cutoff_col is None
    )

    _cdist_lddt_kernel[grid](
        # Pointers
        pred_coords_row,
        pred_coords_col,
        true_coords_row,
        true_coords_col,
        mask_row,
        mask_col,
        atom_indices_row,  # Can be None if USE_EXPLICIT_INDICES_ROW=False
        atom_indices_col,  # Can be None if USE_EXPLICIT_INDICES_COL=False
        cutoff_col,  # Can be None if USE_CUTOFF_COL=False
        out_num,
        out_denom,
        # Shape and strides for pred_coords_row [B_mul, N_row, 3]
        pred_coords_row.shape[0],
        pred_coords_row.shape[1],
        pred_coords_row.stride(0),
        pred_coords_row.stride(1),
        pred_coords_row.stride(2),
        # Shape and strides for pred_coords_col [B_mul, N_col, 3]
        pred_coords_col.shape[0],
        pred_coords_col.shape[1],
        pred_coords_col.stride(0),
        pred_coords_col.stride(1),
        pred_coords_col.stride(2),
        # Shape and strides for true_coords_row [B_mul, N_row, 3]
        true_coords_row.shape[0],
        true_coords_row.shape[1],
        true_coords_row.stride(0),
        true_coords_row.stride(1),
        true_coords_row.stride(2),
        # Shape and strides for true_coords_col [B_mul, N_col, 3]
        true_coords_col.shape[0],
        true_coords_col.shape[1],
        true_coords_col.stride(0),
        true_coords_col.stride(1),
        true_coords_col.stride(2),
        # Shape and strides for mask_row [B, N_row]
        mask_row.shape[0],
        mask_row.shape[1],
        mask_row.stride(0),
        mask_row.stride(1),
        # Shape and strides for mask_col [B, N_col]
        mask_col.shape[0],
        mask_col.shape[1],
        mask_col.stride(0),
        mask_col.stride(1),
        # Shape and strides for atom_indices_row [B, N_row] (only used if USE_EXPLICIT_INDICES_ROW)
        atom_indices_row.shape[0] if use_explicit_indices_row else 0,
        atom_indices_row.shape[1] if use_explicit_indices_row else 0,
        atom_indices_row.stride(0) if use_explicit_indices_row else 0,
        atom_indices_row.stride(1) if use_explicit_indices_row else 0,
        # Shape and strides for atom_indices_col [B, N_col] (only used if USE_EXPLICIT_INDICES_COL)
        atom_indices_col.shape[0] if use_explicit_indices_col else 0,
        atom_indices_col.shape[1] if use_explicit_indices_col else 0,
        atom_indices_col.stride(0) if use_explicit_indices_col else 0,
        atom_indices_col.stride(1) if use_explicit_indices_col else 0,
        # Shape and strides for cutoff_col [B, N_col] (only used if USE_CUTOFF_COL is True)
        cutoff_col.shape[0] if use_cutoff_col else 0,
        cutoff_col.shape[1] if use_cutoff_col else 0,
        cutoff_col.stride(0) if use_cutoff_col else 0,
        cutoff_col.stride(1) if use_cutoff_col else 0,
        # Shape and strides for out_num [B_mul] or [B_mul, N_row] if per_atom
        out_num.shape[0],
        out_num.stride(0),
        out_num.stride(1) if per_atom else 0,
        # Shape and strides for out_denom [B_mul] or [B_mul, N_row] if per_atom
        out_denom.shape[0],
        out_denom.stride(0),
        out_denom.stride(1) if per_atom else 0,
        # Constants
        # When cutoff_col is provided, it overrides the scalar cutoff (pass 0.0 as dummy)
        0.0 if use_cutoff_col else float(cutoff),
        # Block sizes
        BLOCK_M=BLOCK_M,
        BLOCK_N=BLOCK_N,
        # Coordinate dimension size (BLOCK_D computed inside kernel as next power of 2)
        SIZE_DIM_D=dim_d,
        # Flags for implicit arange
        USE_EXPLICIT_INDICES_ROW=use_explicit_indices_row,
        USE_EXPLICIT_INDICES_COL=use_explicit_indices_col,
        # Flag for per-column cutoff
        USE_CUTOFF_COL=use_cutoff_col,
        # Flag for diagonal masking
        DO_MASK_DIAGONAL=do_mask_diagonal,
        # Flag for per-atom output
        PER_ATOM=per_atom,
        # Memory layout orders for make_block_ptr (3-tuple for coords, 2-tuple for masks)
        ORDER_PRED_COORDS_ROW_0=order_pred_coords_row[0],
        ORDER_PRED_COORDS_ROW_1=order_pred_coords_row[1],
        ORDER_PRED_COORDS_ROW_2=order_pred_coords_row[2],
        ORDER_PRED_COORDS_COL_0=order_pred_coords_col[0],
        ORDER_PRED_COORDS_COL_1=order_pred_coords_col[1],
        ORDER_PRED_COORDS_COL_2=order_pred_coords_col[2],
        ORDER_TRUE_COORDS_ROW_0=order_true_coords_row[0],
        ORDER_TRUE_COORDS_ROW_1=order_true_coords_row[1],
        ORDER_TRUE_COORDS_ROW_2=order_true_coords_row[2],
        ORDER_TRUE_COORDS_COL_0=order_true_coords_col[0],
        ORDER_TRUE_COORDS_COL_1=order_true_coords_col[1],
        ORDER_TRUE_COORDS_COL_2=order_true_coords_col[2],
        ORDER_MASK_ROW_0=order_mask_row[0],
        ORDER_MASK_ROW_1=order_mask_row[1],
        ORDER_MASK_COL_0=order_mask_col[0],
        ORDER_MASK_COL_1=order_mask_col[1],
        ORDER_ATOM_INDICES_ROW_0=order_atom_indices_row[0],
        ORDER_ATOM_INDICES_ROW_1=order_atom_indices_row[1],
        ORDER_ATOM_INDICES_COL_0=order_atom_indices_col[0],
        ORDER_ATOM_INDICES_COL_1=order_atom_indices_col[1],
        ORDER_CUTOFF_COL_0=order_cutoff_col[0],
        ORDER_CUTOFF_COL_1=order_cutoff_col[1],
        # Flags for mask multiplicity
        MASK_ROW_HAS_MUL=mask_row_has_mul,
        MASK_COL_HAS_MUL=mask_col_has_mul,
        # Explicit multiplicity factor
        MULTIPLICITY=multiplicity,
        num_warps=4,
        num_stages=3,
    )

    # Final reduction / normalization
    if per_atom:
        # Per-row outputs: out_num, out_denom have shape [B_mul, N_row]
        # mask_no_match: True where there are valid pairs for this row
        mask_no_match = (out_denom > 0).to(input_dtype)

        if return_unnormalized_score:
            # Return raw (out_num, out_denom, mask_no_match) for distributed allreduce
            return out_num.to(input_dtype), out_denom.to(input_dtype), mask_no_match

        # Per-row normalization: score = (eps + sum(dists_to_score * score)) / (eps + sum(dists_to_score))
        # This matches lddt_dist's per_atom=True behavior
        norm = 1.0 / (eps + out_denom)
        score = norm * (eps + out_num)
        if return_denom:
            return score.to(input_dtype), mask_no_match, out_denom.to(input_dtype)
        return score.to(input_dtype), mask_no_match

    # Total outputs: out_num, out_denom have shape [B_mul]
    if return_unnormalized_score:
        # Return raw (out_num, out_denom) for distributed allreduce
        return out_num.to(input_dtype), out_denom.to(input_dtype)

    # Total score normalization
    # Avoid division by zero
    result = out_num / (out_denom + eps)

    # If denominator is 0, result should be 0 (no valid atoms)
    score = torch.where(out_denom > 0, result, torch.zeros_like(result))
    if return_denom:
        return score.to(input_dtype), out_denom.to(input_dtype)
    return score.to(input_dtype)
