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
def _cdist_pde_fwd_kernel(
    # Pointers to inputs (passed as tensors for make_block_ptr)
    pred_pde,  # [B_mul, N_row, N_col, num_bins]
    true_coords_row,  # [B_mul, N_row, 3]
    true_coords_col,  # [B_mul, N_col, 3]
    pred_coords_row,  # [B_mul, N_row, 3]
    pred_coords_col,  # [B_mul, N_col, 3]
    mask_row,  # [B, N_row] or [B_mul, N_row]
    mask_col,  # [B, N_col] or [B_mul, N_col]
    # Output pointers (scalar per batch, use atomic_add)
    out_loss_num_ptr,  # [B_mul]
    out_mask_denom_ptr,  # [B_mul]
    # Shape info for pred_pde [B_mul, N_row, N_col, num_bins]
    B_mul,
    N_row,
    N_col,
    NUM_BINS: tl.constexpr,
    stride_pde_b,
    stride_pde_i,
    stride_pde_j,
    stride_pde_k,
    # Strides for true_coords_row [B_mul, N_row, 3]
    stride_tc_row_b,
    stride_tc_row_n,
    stride_tc_row_d,
    # Strides for true_coords_col [B_mul, N_col, 3]
    stride_tc_col_b,
    stride_tc_col_n,
    stride_tc_col_d,
    # Strides for pred_coords_row [B_mul, N_row, 3]
    stride_pc_row_b,
    stride_pc_row_n,
    stride_pc_row_d,
    # Strides for pred_coords_col [B_mul, N_col, 3]
    stride_pc_col_b,
    stride_pc_col_n,
    stride_pc_col_d,
    # Shape info for mask_row [B_mask_row, N_row]
    B_mask_row,
    stride_mask_row_b,
    stride_mask_row_n,
    # Shape info for mask_col [B_mask_col, N_col]
    B_mask_col,
    stride_mask_col_b,
    stride_mask_col_n,
    # Constants
    max_dist,
    # Block sizes
    BLOCK_M: tl.constexpr,
    BLOCK_N: tl.constexpr,
    # Coordinate dimension
    SIZE_DIM_D: tl.constexpr,
    # Mask multiplicity flags
    MASK_ROW_HAS_MUL: tl.constexpr,
    MASK_COL_HAS_MUL: tl.constexpr,
    MULTIPLICITY: tl.constexpr,
    # Memory layout orders for make_block_ptr (computed from argsort of strides)
    ORDER_TC_ROW_0: tl.constexpr,
    ORDER_TC_ROW_1: tl.constexpr,
    ORDER_TC_ROW_2: tl.constexpr,
    ORDER_TC_COL_0: tl.constexpr,
    ORDER_TC_COL_1: tl.constexpr,
    ORDER_TC_COL_2: tl.constexpr,
    ORDER_PC_ROW_0: tl.constexpr,
    ORDER_PC_ROW_1: tl.constexpr,
    ORDER_PC_ROW_2: tl.constexpr,
    ORDER_PC_COL_0: tl.constexpr,
    ORDER_PC_COL_1: tl.constexpr,
    ORDER_PC_COL_2: tl.constexpr,
    ORDER_MASK_ROW_0: tl.constexpr,
    ORDER_MASK_ROW_1: tl.constexpr,
    ORDER_MASK_COL_0: tl.constexpr,
    ORDER_MASK_COL_1: tl.constexpr,
    ORDER_PDE_0: tl.constexpr,
    ORDER_PDE_1: tl.constexpr,
    ORDER_PDE_2: tl.constexpr,
    ORDER_PDE_3: tl.constexpr,
):
    """
    Forward kernel for computing PDE cross-entropy loss.

    For each (i, j) pair in a tile:
    1. Compute true_d = ||true_coords[i] - true_coords[j]||
    2. Compute pred_d = ||pred_coords[i] - pred_coords[j]||
    3. Compute target_pde = |true_d - pred_d|
    4. Compute bin_index = clamp(floor(target_pde * num_bins / max_dist), max=num_bins-1)
    5. Compute log_softmax of pred_pde logits
    6. Compute cross-entropy: ce_loss = -log_softmax[bin_index]
    7. Accumulate: out_loss_num[b] += sum_{i,j}(ce_loss * mask[i,j])
                   out_mask_denom[b] += sum_{i,j}(mask[i,j])
    """
    tl.static_assert(SIZE_DIM_D == 3, "SIZE_DIM_D must be 3 (3D coordinates)")

    # Block dimension for coordinates (next power of 2 of 3)
    BLOCK_D: tl.constexpr = 4

    # Grid identification
    pid_batch = tl.program_id(0)  # batch * multiplicity
    pid_m = tl.program_id(1)  # row block
    pid_n = tl.program_id(2)  # col block

    # Multiplicity handling for mask broadcasting
    batch_idx = pid_batch // MULTIPLICITY
    batch_idx_mask_row = pid_batch if MASK_ROW_HAS_MUL else batch_idx
    batch_idx_mask_col = pid_batch if MASK_COL_HAS_MUL else batch_idx

    # ============================================
    # 1. Create block pointers and load coordinates
    # ============================================
    # Using make_block_ptr lets the compiler handle index types automatically
    # For 3D coordinate tensors: [B_mul, N, 3], load block [1, BLOCK_M/N, BLOCK_D]
    tc_row_block_ptr = tl.make_block_ptr(
        base=true_coords_row,
        shape=(B_mul, N_row, SIZE_DIM_D),
        strides=(stride_tc_row_b, stride_tc_row_n, stride_tc_row_d),
        offsets=(pid_batch, pid_m * BLOCK_M, 0),
        block_shape=(1, BLOCK_M, BLOCK_D),
        order=(ORDER_TC_ROW_0, ORDER_TC_ROW_1, ORDER_TC_ROW_2),
    )
    tc_col_block_ptr = tl.make_block_ptr(
        base=true_coords_col,
        shape=(B_mul, N_col, SIZE_DIM_D),
        strides=(stride_tc_col_b, stride_tc_col_n, stride_tc_col_d),
        offsets=(pid_batch, pid_n * BLOCK_N, 0),
        block_shape=(1, BLOCK_N, BLOCK_D),
        order=(ORDER_TC_COL_0, ORDER_TC_COL_1, ORDER_TC_COL_2),
    )
    pc_row_block_ptr = tl.make_block_ptr(
        base=pred_coords_row,
        shape=(B_mul, N_row, SIZE_DIM_D),
        strides=(stride_pc_row_b, stride_pc_row_n, stride_pc_row_d),
        offsets=(pid_batch, pid_m * BLOCK_M, 0),
        block_shape=(1, BLOCK_M, BLOCK_D),
        order=(ORDER_PC_ROW_0, ORDER_PC_ROW_1, ORDER_PC_ROW_2),
    )
    pc_col_block_ptr = tl.make_block_ptr(
        base=pred_coords_col,
        shape=(B_mul, N_col, SIZE_DIM_D),
        strides=(stride_pc_col_b, stride_pc_col_n, stride_pc_col_d),
        offsets=(pid_batch, pid_n * BLOCK_N, 0),
        block_shape=(1, BLOCK_N, BLOCK_D),
        order=(ORDER_PC_COL_0, ORDER_PC_COL_1, ORDER_PC_COL_2),
    )

    # Load coordinates with boundary_check and reshape to squeeze batch dim
    # boundary_check=(1, 2) checks N and D dims; batch dim is always in bounds
    true_row = tl.reshape(
        tl.load(tc_row_block_ptr, boundary_check=(1, 2), padding_option="zero"),
        (BLOCK_M, BLOCK_D),
    )
    true_col = tl.reshape(
        tl.load(tc_col_block_ptr, boundary_check=(1, 2), padding_option="zero"),
        (BLOCK_N, BLOCK_D),
    )
    pred_row = tl.reshape(
        tl.load(pc_row_block_ptr, boundary_check=(1, 2), padding_option="zero"),
        (BLOCK_M, BLOCK_D),
    )
    pred_col = tl.reshape(
        tl.load(pc_col_block_ptr, boundary_check=(1, 2), padding_option="zero"),
        (BLOCK_N, BLOCK_D),
    )

    # ============================================
    # 2. Compute pairwise distances [BLOCK_M, BLOCK_N]
    # ============================================
    # true_d[i,j] = ||true_row[i] - true_col[j]||
    delta_true = true_row[:, None, :] - true_col[None, :, :]  # [BLOCK_M, BLOCK_N, BLOCK_D]
    d_true_sq = tl.sum(delta_true * delta_true, axis=2)  # [BLOCK_M, BLOCK_N]
    d_true = tl.sqrt(d_true_sq)

    # pred_d[i,j] = ||pred_row[i] - pred_col[j]||
    delta_pred = pred_row[:, None, :] - pred_col[None, :, :]  # [BLOCK_M, BLOCK_N, BLOCK_D]
    d_pred_sq = tl.sum(delta_pred * delta_pred, axis=2)  # [BLOCK_M, BLOCK_N]
    d_pred = tl.sqrt(d_pred_sq)

    # ============================================
    # 3. Compute target_pde and bin_index [BLOCK_M, BLOCK_N]
    # ============================================
    target_pde = tl.abs(d_true - d_pred)

    # bin_index = clamp(floor(target_pde * num_bins / max_dist), max=num_bins-1)
    bin_index_float = target_pde * NUM_BINS / max_dist
    bin_index = tl.minimum(tl.floor(bin_index_float).to(tl.int32), NUM_BINS - 1)

    # ============================================
    # 4. Load masks using block_ptr
    # ============================================
    # For 2D mask tensors: [B or B_mul, N], load block [1, BLOCK_M/N]
    mask_row_block_ptr = tl.make_block_ptr(
        base=mask_row,
        shape=(B_mask_row, N_row),
        strides=(stride_mask_row_b, stride_mask_row_n),
        offsets=(batch_idx_mask_row, pid_m * BLOCK_M),
        block_shape=(1, BLOCK_M),
        order=(ORDER_MASK_ROW_0, ORDER_MASK_ROW_1),
    )
    mask_col_block_ptr = tl.make_block_ptr(
        base=mask_col,
        shape=(B_mask_col, N_col),
        strides=(stride_mask_col_b, stride_mask_col_n),
        offsets=(batch_idx_mask_col, pid_n * BLOCK_N),
        block_shape=(1, BLOCK_N),
        order=(ORDER_MASK_COL_0, ORDER_MASK_COL_1),
    )

    # Load masks and reshape to squeeze batch dim
    m_row = tl.reshape(
        tl.load(mask_row_block_ptr, boundary_check=(1,), padding_option="zero"),
        (BLOCK_M,),
    )
    m_col = tl.reshape(
        tl.load(mask_col_block_ptr, boundary_check=(1,), padding_option="zero"),
        (BLOCK_N,),
    )

    # Combined pair mask [BLOCK_M, BLOCK_N]
    pair_mask = m_row[:, None] * m_col[None, :]

    # ============================================
    # 5. Load pred_pde logits using block_ptr and compute cross-entropy
    # ============================================
    # For 4D pred_pde: [B_mul, N_row, N_col, num_bins], load block [1, BLOCK_M, BLOCK_N, NUM_BINS]
    pde_block_ptr = tl.make_block_ptr(
        base=pred_pde,
        shape=(B_mul, N_row, N_col, NUM_BINS),
        strides=(stride_pde_b, stride_pde_i, stride_pde_j, stride_pde_k),
        offsets=(pid_batch, pid_m * BLOCK_M, pid_n * BLOCK_N, 0),
        block_shape=(1, BLOCK_M, BLOCK_N, NUM_BINS),
        order=(ORDER_PDE_0, ORDER_PDE_1, ORDER_PDE_2, ORDER_PDE_3),
    )

    # Load and reshape to squeeze batch dim: [1, BLOCK_M, BLOCK_N, NUM_BINS] -> [BLOCK_M, BLOCK_N, NUM_BINS]
    logits = tl.reshape(
        tl.load(pde_block_ptr, boundary_check=(1, 2), padding_option="zero"),
        (BLOCK_M, BLOCK_N, NUM_BINS),
    )

    # Compute log_softmax along the last dimension
    # log_softmax(x) = x - log(sum(exp(x)))
    max_logits = tl.max(logits, axis=2)[:, :, None]  # [BLOCK_M, BLOCK_N, 1]
    logits_shifted = logits - max_logits
    exp_logits = tl.exp(logits_shifted)
    sum_exp = tl.sum(exp_logits, axis=2)[:, :, None]  # [BLOCK_M, BLOCK_N, 1]
    log_sum_exp = tl.log(sum_exp)
    log_probs = logits_shifted - log_sum_exp  # [BLOCK_M, BLOCK_N, NUM_BINS]

    # Gather log_probs at bin_index for each (i, j)
    # We need log_probs[i, j, bin_index[i, j]]
    # Use tl.gather along axis=2 (the bins dimension)
    selected_log_prob = tl.gather(log_probs, bin_index[:, :, None], axis=2)  # [BLOCK_M, BLOCK_N, 1]
    selected_log_prob = tl.reshape(selected_log_prob, (BLOCK_M, BLOCK_N))  # [BLOCK_M, BLOCK_N]

    # Cross-entropy loss: -log_prob
    ce_loss = -selected_log_prob  # [BLOCK_M, BLOCK_N]

    # Apply pair mask
    ce_loss_masked = ce_loss * pair_mask  # [BLOCK_M, BLOCK_N]

    # ============================================
    # 6. Accumulate tile sum and atomic add (scalar per batch)
    # ============================================
    # Sum over both row and column dimensions (full tile contribution)
    loss_tile_sum = tl.sum(ce_loss_masked)  # scalar
    denom_tile_sum = tl.sum(pair_mask)  # scalar

    # Atomic add to per-batch output (scalar)
    out_loss_ptr = out_loss_num_ptr + pid_batch
    out_denom_ptr = out_mask_denom_ptr + pid_batch

    tl.atomic_add(out_loss_ptr, loss_tile_sum, sem="relaxed")
    tl.atomic_add(out_denom_ptr, denom_tile_sum, sem="relaxed")


@triton.jit
def _cdist_pde_bwd_kernel(
    # Pointers to inputs (passed as tensors for make_block_ptr)
    grad_out_loss_num,  # [B_mul] - upstream gradient (scalar per batch)
    pred_pde,  # [B_mul, N_row, N_col, num_bins]
    true_coords_row,  # [B_mul, N_row, 3]
    true_coords_col,  # [B_mul, N_col, 3]
    pred_coords_row,  # [B_mul, N_row, 3]
    pred_coords_col,  # [B_mul, N_col, 3]
    mask_row,  # [B, N_row] or [B_mul, N_row]
    mask_col,  # [B, N_col] or [B_mul, N_col]
    # Output tensor (for make_block_ptr)
    grad_pred_pde,  # [B_mul, N_row, N_col, num_bins]
    # Shape info for pred_pde [B_mul, N_row, N_col, num_bins]
    B_mul,
    N_row,
    N_col,
    NUM_BINS: tl.constexpr,
    stride_pde_b,
    stride_pde_i,
    stride_pde_j,
    stride_pde_k,
    # Strides for true_coords_row [B_mul, N_row, 3]
    stride_tc_row_b,
    stride_tc_row_n,
    stride_tc_row_d,
    # Strides for true_coords_col [B_mul, N_col, 3]
    stride_tc_col_b,
    stride_tc_col_n,
    stride_tc_col_d,
    # Strides for pred_coords_row [B_mul, N_row, 3]
    stride_pc_row_b,
    stride_pc_row_n,
    stride_pc_row_d,
    # Strides for pred_coords_col [B_mul, N_col, 3]
    stride_pc_col_b,
    stride_pc_col_n,
    stride_pc_col_d,
    # Shape info for mask_row [B_mask_row, N_row]
    B_mask_row,
    stride_mask_row_b,
    stride_mask_row_n,
    # Shape info for mask_col [B_mask_col, N_col]
    B_mask_col,
    stride_mask_col_b,
    stride_mask_col_n,
    # Stride for grad_out_loss_num [B_mul]
    stride_grad_out_b,
    # Constants
    max_dist,
    # Block sizes
    BLOCK_M: tl.constexpr,
    BLOCK_N: tl.constexpr,
    # Coordinate dimension
    SIZE_DIM_D: tl.constexpr,
    # Mask multiplicity flags
    MASK_ROW_HAS_MUL: tl.constexpr,
    MASK_COL_HAS_MUL: tl.constexpr,
    MULTIPLICITY: tl.constexpr,
    # Memory layout orders for make_block_ptr
    ORDER_TC_ROW_0: tl.constexpr,
    ORDER_TC_ROW_1: tl.constexpr,
    ORDER_TC_ROW_2: tl.constexpr,
    ORDER_TC_COL_0: tl.constexpr,
    ORDER_TC_COL_1: tl.constexpr,
    ORDER_TC_COL_2: tl.constexpr,
    ORDER_PC_ROW_0: tl.constexpr,
    ORDER_PC_ROW_1: tl.constexpr,
    ORDER_PC_ROW_2: tl.constexpr,
    ORDER_PC_COL_0: tl.constexpr,
    ORDER_PC_COL_1: tl.constexpr,
    ORDER_PC_COL_2: tl.constexpr,
    ORDER_MASK_ROW_0: tl.constexpr,
    ORDER_MASK_ROW_1: tl.constexpr,
    ORDER_MASK_COL_0: tl.constexpr,
    ORDER_MASK_COL_1: tl.constexpr,
    ORDER_PDE_0: tl.constexpr,
    ORDER_PDE_1: tl.constexpr,
    ORDER_PDE_2: tl.constexpr,
    ORDER_PDE_3: tl.constexpr,
):
    """
    Backward kernel for PDE cross-entropy loss.

    Computes gradient w.r.t. pred_pde only.

    For cross-entropy loss with log_softmax:
        loss = -log_softmax(logits)[target]
        d_loss/d_logits = softmax(logits) - one_hot(target)

    With upstream gradient (scalar per batch) and mask:
        grad_pred_pde[i,j,:] = (softmax - one_hot(bin_idx)) * mask[i,j] * grad_out[b]
    """
    tl.static_assert(SIZE_DIM_D == 3, "SIZE_DIM_D must be 3 (3D coordinates)")

    BLOCK_D: tl.constexpr = 4

    # Grid identification
    pid_batch = tl.program_id(0)
    pid_m = tl.program_id(1)
    pid_n = tl.program_id(2)

    # Multiplicity handling
    batch_idx = pid_batch // MULTIPLICITY
    batch_idx_mask_row = pid_batch if MASK_ROW_HAS_MUL else batch_idx
    batch_idx_mask_col = pid_batch if MASK_COL_HAS_MUL else batch_idx

    # ============================================
    # 1. Create block pointers and load coordinates
    # ============================================
    tc_row_block_ptr = tl.make_block_ptr(
        base=true_coords_row,
        shape=(B_mul, N_row, SIZE_DIM_D),
        strides=(stride_tc_row_b, stride_tc_row_n, stride_tc_row_d),
        offsets=(pid_batch, pid_m * BLOCK_M, 0),
        block_shape=(1, BLOCK_M, BLOCK_D),
        order=(ORDER_TC_ROW_0, ORDER_TC_ROW_1, ORDER_TC_ROW_2),
    )
    tc_col_block_ptr = tl.make_block_ptr(
        base=true_coords_col,
        shape=(B_mul, N_col, SIZE_DIM_D),
        strides=(stride_tc_col_b, stride_tc_col_n, stride_tc_col_d),
        offsets=(pid_batch, pid_n * BLOCK_N, 0),
        block_shape=(1, BLOCK_N, BLOCK_D),
        order=(ORDER_TC_COL_0, ORDER_TC_COL_1, ORDER_TC_COL_2),
    )
    pc_row_block_ptr = tl.make_block_ptr(
        base=pred_coords_row,
        shape=(B_mul, N_row, SIZE_DIM_D),
        strides=(stride_pc_row_b, stride_pc_row_n, stride_pc_row_d),
        offsets=(pid_batch, pid_m * BLOCK_M, 0),
        block_shape=(1, BLOCK_M, BLOCK_D),
        order=(ORDER_PC_ROW_0, ORDER_PC_ROW_1, ORDER_PC_ROW_2),
    )
    pc_col_block_ptr = tl.make_block_ptr(
        base=pred_coords_col,
        shape=(B_mul, N_col, SIZE_DIM_D),
        strides=(stride_pc_col_b, stride_pc_col_n, stride_pc_col_d),
        offsets=(pid_batch, pid_n * BLOCK_N, 0),
        block_shape=(1, BLOCK_N, BLOCK_D),
        order=(ORDER_PC_COL_0, ORDER_PC_COL_1, ORDER_PC_COL_2),
    )

    # Load coordinates with boundary_check and reshape
    true_row = tl.reshape(
        tl.load(tc_row_block_ptr, boundary_check=(1, 2), padding_option="zero"),
        (BLOCK_M, BLOCK_D),
    )
    true_col = tl.reshape(
        tl.load(tc_col_block_ptr, boundary_check=(1, 2), padding_option="zero"),
        (BLOCK_N, BLOCK_D),
    )
    pred_row = tl.reshape(
        tl.load(pc_row_block_ptr, boundary_check=(1, 2), padding_option="zero"),
        (BLOCK_M, BLOCK_D),
    )
    pred_col = tl.reshape(
        tl.load(pc_col_block_ptr, boundary_check=(1, 2), padding_option="zero"),
        (BLOCK_N, BLOCK_D),
    )

    # ============================================
    # 2. Recompute distances and bin_index
    # ============================================
    delta_true = true_row[:, None, :] - true_col[None, :, :]
    d_true_sq = tl.sum(delta_true * delta_true, axis=2)
    d_true = tl.sqrt(d_true_sq)

    delta_pred = pred_row[:, None, :] - pred_col[None, :, :]
    d_pred_sq = tl.sum(delta_pred * delta_pred, axis=2)
    d_pred = tl.sqrt(d_pred_sq)

    target_pde = tl.abs(d_true - d_pred)
    bin_index_float = target_pde * NUM_BINS / max_dist
    bin_index = tl.minimum(tl.floor(bin_index_float).to(tl.int32), NUM_BINS - 1)  # [BLOCK_M, BLOCK_N]

    # ============================================
    # 3. Load masks using block_ptr
    # ============================================
    mask_row_block_ptr = tl.make_block_ptr(
        base=mask_row,
        shape=(B_mask_row, N_row),
        strides=(stride_mask_row_b, stride_mask_row_n),
        offsets=(batch_idx_mask_row, pid_m * BLOCK_M),
        block_shape=(1, BLOCK_M),
        order=(ORDER_MASK_ROW_0, ORDER_MASK_ROW_1),
    )
    mask_col_block_ptr = tl.make_block_ptr(
        base=mask_col,
        shape=(B_mask_col, N_col),
        strides=(stride_mask_col_b, stride_mask_col_n),
        offsets=(batch_idx_mask_col, pid_n * BLOCK_N),
        block_shape=(1, BLOCK_N),
        order=(ORDER_MASK_COL_0, ORDER_MASK_COL_1),
    )

    m_row = tl.reshape(
        tl.load(mask_row_block_ptr, boundary_check=(1,), padding_option="zero"),
        (BLOCK_M,),
    )
    m_col = tl.reshape(
        tl.load(mask_col_block_ptr, boundary_check=(1,), padding_option="zero"),
        (BLOCK_N,),
    )

    pair_mask = m_row[:, None] * m_col[None, :]  # [BLOCK_M, BLOCK_N]

    # ============================================
    # 4. Load upstream gradient (scalar per batch) using make_block_ptr
    # ============================================
    # grad_out_loss_num is [B_mul], load single scalar for this batch
    grad_out_block_ptr = tl.make_block_ptr(
        base=grad_out_loss_num,
        shape=(B_mul,),
        strides=(stride_grad_out_b,),
        offsets=(pid_batch,),
        block_shape=(1,),
        order=(0,),
    )
    grad_out = tl.reshape(tl.load(grad_out_block_ptr), ())  # scalar

    # ============================================
    # 5. Load pred_pde logits using block_ptr and compute gradient
    # ============================================
    pde_block_ptr = tl.make_block_ptr(
        base=pred_pde,
        shape=(B_mul, N_row, N_col, NUM_BINS),
        strides=(stride_pde_b, stride_pde_i, stride_pde_j, stride_pde_k),
        offsets=(pid_batch, pid_m * BLOCK_M, pid_n * BLOCK_N, 0),
        block_shape=(1, BLOCK_M, BLOCK_N, NUM_BINS),
        order=(ORDER_PDE_0, ORDER_PDE_1, ORDER_PDE_2, ORDER_PDE_3),
    )
    logits = tl.reshape(
        tl.load(pde_block_ptr, boundary_check=(1, 2), padding_option="zero"),
        (BLOCK_M, BLOCK_N, NUM_BINS),
    )

    # Compute softmax
    max_logits = tl.max(logits, axis=2)[:, :, None]  # [BLOCK_M, BLOCK_N, 1]
    logits_shifted = logits - max_logits
    exp_logits = tl.exp(logits_shifted)
    sum_exp = tl.sum(exp_logits, axis=2)[:, :, None]  # [BLOCK_M, BLOCK_N, 1]
    softmax_probs = exp_logits / sum_exp  # [BLOCK_M, BLOCK_N, NUM_BINS]

    # Create one_hot for targets: [BLOCK_M, BLOCK_N, NUM_BINS]
    offs_k = tl.arange(0, NUM_BINS)
    one_hot = (offs_k[None, None, :] == bin_index[:, :, None]).to(tl.float32)

    # Gradient of cross-entropy w.r.t. logits: softmax - one_hot
    grad_logits = softmax_probs - one_hot  # [BLOCK_M, BLOCK_N, NUM_BINS]

    # Apply mask and upstream gradient (scalar per batch)
    # grad_pred_pde[i,j,:] = grad_logits * mask[i,j] * grad_out[b]
    scale = pair_mask * grad_out  # [BLOCK_M, BLOCK_N] (scalar broadcast)
    grad_logits_scaled = grad_logits * scale[:, :, None]  # [BLOCK_M, BLOCK_N, NUM_BINS]

    # ============================================
    # 6. Store gradient using block_ptr
    # ============================================
    # Create block_ptr for gradient output (same layout as pred_pde)
    grad_pde_block_ptr = tl.make_block_ptr(
        base=grad_pred_pde,
        shape=(B_mul, N_row, N_col, NUM_BINS),
        strides=(stride_pde_b, stride_pde_i, stride_pde_j, stride_pde_k),
        offsets=(pid_batch, pid_m * BLOCK_M, pid_n * BLOCK_N, 0),
        block_shape=(1, BLOCK_M, BLOCK_N, NUM_BINS),
        order=(ORDER_PDE_0, ORDER_PDE_1, ORDER_PDE_2, ORDER_PDE_3),
    )
    # Expand to 4D for store: [BLOCK_M, BLOCK_N, NUM_BINS] -> [1, BLOCK_M, BLOCK_N, NUM_BINS]
    grad_logits_4d = tl.reshape(grad_logits_scaled, (1, BLOCK_M, BLOCK_N, NUM_BINS))
    tl.store(grad_pde_block_ptr, grad_logits_4d, boundary_check=(1, 2))


class _CdistPDEImpl(torch.autograd.Function):
    """Autograd function wrapping forward and backward Triton kernels."""

    @staticmethod
    def forward(
        ctx,
        pred_pde,
        true_coords_row,
        true_coords_col,
        pred_coords_row,
        pred_coords_col,
        mask_row,
        mask_col,
        multiplicity,
        num_bins,
        max_dist,
    ):
        """
        Forward pass: compute PDE cross-entropy fully summed per batch.

        Args:
            pred_pde: [B_mul, N_row, N_col, num_bins] - predicted logits
            true_coords_row: [B_mul, N_row, 3] - true coordinates for rows
            true_coords_col: [B_mul, N_col, 3] - true coordinates for columns
            pred_coords_row: [B_mul, N_row, 3] - predicted coordinates for rows
            pred_coords_col: [B_mul, N_col, 3] - predicted coordinates for columns
            mask_row: [B, N_row] or [B_mul, N_row] - row mask
            mask_col: [B, N_col] or [B_mul, N_col] - column mask
            multiplicity: int - B_mul = B * multiplicity
            num_bins: int - number of bins
            max_dist: float - maximum distance for binning

        Returns:
            out_loss_num: [B_mul] - sum of CE loss over all (i,j) pairs per batch
            out_mask_denom: [B_mul] - sum of mask over all (i,j) pairs per batch
        """
        B_mul, N_row, N_col, num_bins_tensor = pred_pde.shape
        device = pred_pde.device

        # Validate num_bins
        if num_bins_tensor != num_bins:
            raise ValueError(f"pred_pde num_bins mismatch: got {num_bins_tensor}, expected {num_bins}")

        # Compute B from multiplicity
        if B_mul % multiplicity != 0:
            raise ValueError(f"B_mul ({B_mul}) must be divisible by multiplicity ({multiplicity})")
        B = B_mul // multiplicity

        # Validate that coordinates and masks don't require gradients
        # (gradient flow is broken by .long() in bin_index computation)
        if true_coords_row.requires_grad:
            raise ValueError(
                "true_coords_row should not require gradients " "(gradient flow is broken by bin_index computation)"
            )
        if true_coords_col.requires_grad:
            raise ValueError(
                "true_coords_col should not require gradients " "(gradient flow is broken by bin_index computation)"
            )
        if pred_coords_row.requires_grad:
            raise ValueError(
                "pred_coords_row should not require gradients " "(gradient flow is broken by bin_index computation)"
            )
        if pred_coords_col.requires_grad:
            raise ValueError(
                "pred_coords_col should not require gradients " "(gradient flow is broken by bin_index computation)"
            )
        if mask_row.requires_grad:
            raise ValueError("mask_row should not require gradients")
        if mask_col.requires_grad:
            raise ValueError("mask_col should not require gradients")

        # Validate coordinate dimensions
        if true_coords_row.shape[-1] != 3:
            raise ValueError(f"Coordinate dimension must be 3, got true_coords_row shape {true_coords_row.shape}")

        # Validate coordinate shapes
        if true_coords_row.shape != (B_mul, N_row, 3):
            raise ValueError(
                f"true_coords_row shape mismatch: got {true_coords_row.shape}, " f"expected ({B_mul}, {N_row}, 3)"
            )
        if true_coords_col.shape != (B_mul, N_col, 3):
            raise ValueError(
                f"true_coords_col shape mismatch: got {true_coords_col.shape}, " f"expected ({B_mul}, {N_col}, 3)"
            )
        if pred_coords_row.shape != (B_mul, N_row, 3):
            raise ValueError(
                f"pred_coords_row shape mismatch: got {pred_coords_row.shape}, " f"expected ({B_mul}, {N_row}, 3)"
            )
        if pred_coords_col.shape != (B_mul, N_col, 3):
            raise ValueError(
                f"pred_coords_col shape mismatch: got {pred_coords_col.shape}, " f"expected ({B_mul}, {N_col}, 3)"
            )

        # Check mask dimensions
        mask_row_has_mul = mask_row.shape[0] == B_mul
        mask_col_has_mul = mask_col.shape[0] == B_mul

        # Validate mask_row shape
        if mask_row.shape[0] not in (B, B_mul):
            raise ValueError(
                f"mask_row batch dimension must be B ({B}) or B_mul ({B_mul}), " f"got {mask_row.shape[0]}"
            )
        if mask_row.shape[1] != N_row:
            raise ValueError(f"mask_row N dimension mismatch: got {mask_row.shape[1]}, expected {N_row}")

        # Validate mask_col shape
        if mask_col.shape[0] not in (B, B_mul):
            raise ValueError(
                f"mask_col batch dimension must be B ({B}) or B_mul ({B_mul}), " f"got {mask_col.shape[0]}"
            )
        if mask_col.shape[1] != N_col:
            raise ValueError(f"mask_col N dimension mismatch: got {mask_col.shape[1]}, expected {N_col}")

        # Don't materialize zero gradients for non-differentiable outputs (out_mask_denom)
        # This makes grad_out_mask_denom be None in backward() instead of zeros
        ctx.set_materialize_grads(False)

        # Output buffers - scalar per batch [B_mul]
        # Use float64 if input is float64, otherwise float32
        output_dtype = pred_pde.dtype if pred_pde.dtype == torch.float64 else torch.float32
        out_loss_num = torch.zeros(B_mul, device=device, dtype=output_dtype)
        out_mask_denom = torch.zeros(B_mul, device=device, dtype=output_dtype)

        # Block sizes (tuned for N=4096, num_bins=64)
        BLOCK_M = 8
        BLOCK_N = 8

        # Compute memory layout order for make_block_ptr
        # order = argsort of strides (ascending), giving fastest-varying dim first
        order_tc_row = tuple(torch.tensor(true_coords_row.stride()).argsort().tolist())
        order_tc_col = tuple(torch.tensor(true_coords_col.stride()).argsort().tolist())
        order_pc_row = tuple(torch.tensor(pred_coords_row.stride()).argsort().tolist())
        order_pc_col = tuple(torch.tensor(pred_coords_col.stride()).argsort().tolist())
        order_mask_row = tuple(torch.tensor(mask_row.stride()).argsort().tolist())
        order_mask_col = tuple(torch.tensor(mask_col.stride()).argsort().tolist())
        order_pde = tuple(torch.tensor(pred_pde.stride()).argsort().tolist())

        # Grid
        grid = (B_mul, triton.cdiv(N_row, BLOCK_M), triton.cdiv(N_col, BLOCK_N))

        _cdist_pde_fwd_kernel[grid](
            pred_pde,
            true_coords_row,
            true_coords_col,
            pred_coords_row,
            pred_coords_col,
            mask_row,
            mask_col,
            out_loss_num,
            out_mask_denom,
            # Shape info
            B_mul,
            N_row,
            N_col,
            num_bins,  # constexpr
            pred_pde.stride(0),
            pred_pde.stride(1),
            pred_pde.stride(2),
            pred_pde.stride(3),
            # Coord strides
            true_coords_row.stride(0),
            true_coords_row.stride(1),
            true_coords_row.stride(2),
            true_coords_col.stride(0),
            true_coords_col.stride(1),
            true_coords_col.stride(2),
            pred_coords_row.stride(0),
            pred_coords_row.stride(1),
            pred_coords_row.stride(2),
            pred_coords_col.stride(0),
            pred_coords_col.stride(1),
            pred_coords_col.stride(2),
            # Mask info
            mask_row.shape[0],
            mask_row.stride(0),
            mask_row.stride(1),
            mask_col.shape[0],
            mask_col.stride(0),
            mask_col.stride(1),
            # Constants
            max_dist,
            # Block sizes (tuned for N=4096, num_bins=64)
            BLOCK_M=BLOCK_M,
            BLOCK_N=BLOCK_N,
            SIZE_DIM_D=3,
            # Flags
            MASK_ROW_HAS_MUL=mask_row_has_mul,
            MASK_COL_HAS_MUL=mask_col_has_mul,
            MULTIPLICITY=multiplicity,
            # Memory layout orders
            ORDER_TC_ROW_0=order_tc_row[0],
            ORDER_TC_ROW_1=order_tc_row[1],
            ORDER_TC_ROW_2=order_tc_row[2],
            ORDER_TC_COL_0=order_tc_col[0],
            ORDER_TC_COL_1=order_tc_col[1],
            ORDER_TC_COL_2=order_tc_col[2],
            ORDER_PC_ROW_0=order_pc_row[0],
            ORDER_PC_ROW_1=order_pc_row[1],
            ORDER_PC_ROW_2=order_pc_row[2],
            ORDER_PC_COL_0=order_pc_col[0],
            ORDER_PC_COL_1=order_pc_col[1],
            ORDER_PC_COL_2=order_pc_col[2],
            ORDER_MASK_ROW_0=order_mask_row[0],
            ORDER_MASK_ROW_1=order_mask_row[1],
            ORDER_MASK_COL_0=order_mask_col[0],
            ORDER_MASK_COL_1=order_mask_col[1],
            ORDER_PDE_0=order_pde[0],
            ORDER_PDE_1=order_pde[1],
            ORDER_PDE_2=order_pde[2],
            ORDER_PDE_3=order_pde[3],
            num_warps=4,
            num_stages=3,
        )

        # Save for backward
        ctx.save_for_backward(
            pred_pde,
            true_coords_row,
            true_coords_col,
            pred_coords_row,
            pred_coords_col,
            mask_row,
            mask_col,
        )
        ctx.multiplicity = multiplicity
        ctx.num_bins = num_bins
        ctx.max_dist = max_dist
        ctx.mask_row_has_mul = mask_row_has_mul
        ctx.mask_col_has_mul = mask_col_has_mul

        # out_mask_denom only depends on mask (not pred_pde), so no gradient
        ctx.mark_non_differentiable(out_mask_denom)

        return out_loss_num, out_mask_denom

    @staticmethod
    def backward(ctx, grad_out_loss_num, grad_out_mask_denom):
        """
        Backward pass: compute gradient w.r.t. pred_pde.

        grad_out_mask_denom should be None since out_mask_denom is marked
        non-differentiable and ctx.set_materialize_grads(False) is set.
        """
        # Validate grad_out_mask_denom is None (out_mask_denom is non-differentiable)
        if grad_out_mask_denom is not None:
            raise ValueError(
                "grad_out_mask_denom should be None since out_mask_denom is "
                "marked non-differentiable (it only depends on mask, not pred_pde)"
            )

        (
            pred_pde,
            true_coords_row,
            true_coords_col,
            pred_coords_row,
            pred_coords_col,
            mask_row,
            mask_col,
        ) = ctx.saved_tensors

        multiplicity = ctx.multiplicity
        num_bins = ctx.num_bins
        max_dist = ctx.max_dist
        mask_row_has_mul = ctx.mask_row_has_mul
        mask_col_has_mul = ctx.mask_col_has_mul

        B_mul, N_row, N_col, _ = pred_pde.shape

        # Validate grad_out_loss_num shape (scalar per batch)
        if grad_out_loss_num.shape != (B_mul,):
            raise ValueError(
                f"grad_out_loss_num shape mismatch: got {grad_out_loss_num.shape}, " f"expected ({B_mul},)"
            )

        # Output gradient buffer
        grad_pred_pde = torch.zeros_like(pred_pde)

        # Block sizes (tuned for N=4096, num_bins=64)
        BLOCK_M = 4
        BLOCK_N = 4

        # Compute memory layout order for make_block_ptr
        order_tc_row = tuple(torch.tensor(true_coords_row.stride()).argsort().tolist())
        order_tc_col = tuple(torch.tensor(true_coords_col.stride()).argsort().tolist())
        order_pc_row = tuple(torch.tensor(pred_coords_row.stride()).argsort().tolist())
        order_pc_col = tuple(torch.tensor(pred_coords_col.stride()).argsort().tolist())
        order_mask_row = tuple(torch.tensor(mask_row.stride()).argsort().tolist())
        order_mask_col = tuple(torch.tensor(mask_col.stride()).argsort().tolist())
        order_pde = tuple(torch.tensor(pred_pde.stride()).argsort().tolist())
        grad_out_contiguous = grad_out_loss_num.contiguous()

        # Grid
        grid = (B_mul, triton.cdiv(N_row, BLOCK_M), triton.cdiv(N_col, BLOCK_N))

        _cdist_pde_bwd_kernel[grid](
            grad_out_contiguous,
            pred_pde,
            true_coords_row,
            true_coords_col,
            pred_coords_row,
            pred_coords_col,
            mask_row,
            mask_col,
            grad_pred_pde,
            # Shape info
            B_mul,
            N_row,
            N_col,
            num_bins,
            pred_pde.stride(0),
            pred_pde.stride(1),
            pred_pde.stride(2),
            pred_pde.stride(3),
            # Coord strides
            true_coords_row.stride(0),
            true_coords_row.stride(1),
            true_coords_row.stride(2),
            true_coords_col.stride(0),
            true_coords_col.stride(1),
            true_coords_col.stride(2),
            pred_coords_row.stride(0),
            pred_coords_row.stride(1),
            pred_coords_row.stride(2),
            pred_coords_col.stride(0),
            pred_coords_col.stride(1),
            pred_coords_col.stride(2),
            # Mask info
            mask_row.shape[0],
            mask_row.stride(0),
            mask_row.stride(1),
            mask_col.shape[0],
            mask_col.stride(0),
            mask_col.stride(1),
            # Stride for grad_out_loss_num
            grad_out_contiguous.stride(0),
            # Constants
            max_dist,
            # Block sizes (tuned for N=4096, num_bins=64)
            BLOCK_M=BLOCK_M,
            BLOCK_N=BLOCK_N,
            SIZE_DIM_D=3,
            # Flags
            MASK_ROW_HAS_MUL=mask_row_has_mul,
            MASK_COL_HAS_MUL=mask_col_has_mul,
            MULTIPLICITY=multiplicity,
            # Memory layout orders
            ORDER_TC_ROW_0=order_tc_row[0],
            ORDER_TC_ROW_1=order_tc_row[1],
            ORDER_TC_ROW_2=order_tc_row[2],
            ORDER_TC_COL_0=order_tc_col[0],
            ORDER_TC_COL_1=order_tc_col[1],
            ORDER_TC_COL_2=order_tc_col[2],
            ORDER_PC_ROW_0=order_pc_row[0],
            ORDER_PC_ROW_1=order_pc_row[1],
            ORDER_PC_ROW_2=order_pc_row[2],
            ORDER_PC_COL_0=order_pc_col[0],
            ORDER_PC_COL_1=order_pc_col[1],
            ORDER_PC_COL_2=order_pc_col[2],
            ORDER_MASK_ROW_0=order_mask_row[0],
            ORDER_MASK_ROW_1=order_mask_row[1],
            ORDER_MASK_COL_0=order_mask_col[0],
            ORDER_MASK_COL_1=order_mask_col[1],
            ORDER_PDE_0=order_pde[0],
            ORDER_PDE_1=order_pde[1],
            ORDER_PDE_2=order_pde[2],
            ORDER_PDE_3=order_pde[3],
            num_warps=4,
            num_stages=2,
        )

        # Return gradients (None for non-differentiable inputs)
        return grad_pred_pde, None, None, None, None, None, None, None, None, None


def cdist_pde(
    pred_pde,
    true_coords_row,
    true_coords_col,
    pred_coords_row,
    pred_coords_col,
    mask_row,
    mask_col,
    multiplicity,
    num_bins=64,
    max_dist=32.0,
):
    """
    Compute PDE cross-entropy loss without materializing O(N_token^2) distance matrices.

    This function computes the cross-entropy portion of the PDE loss directly from
    coordinates, fusing distance computation, binning, and cross-entropy into a
    single kernel that only uses O(tile_size^2) local memory.

    The computation is equivalent to:
        true_d = torch.cdist(true_coords_row, true_coords_col)
        pred_d = torch.cdist(pred_coords_row, pred_coords_col)
        target_pde = torch.abs(true_d - pred_d)
        bin_index = torch.clamp(torch.floor(target_pde * num_bins / max_dist).long(), max=num_bins-1)
        one_hot = F.one_hot(bin_index, num_classes=num_bins)
        errors = -torch.sum(one_hot * F.log_softmax(pred_pde, dim=-1), dim=-1)
        out_loss_num = torch.sum(errors * mask, dim=(-2, -1))  # sum over both dims
        out_mask_denom = torch.sum(mask, dim=(-2, -1))  # sum over both dims

    Parameters
    ----------
    pred_pde : torch.Tensor
        Predicted PDE logits, shape [B_mul, N_row, N_col, num_bins].
        This is the only input that requires gradients.
    true_coords_row : torch.Tensor
        True coordinates for row tokens, shape [B_mul, N_row, 3].
    true_coords_col : torch.Tensor
        True coordinates for column tokens, shape [B_mul, N_col, 3].
    pred_coords_row : torch.Tensor
        Predicted coordinates for row tokens, shape [B_mul, N_row, 3].
    pred_coords_col : torch.Tensor
        Predicted coordinates for column tokens, shape [B_mul, N_col, 3].
    mask_row : torch.Tensor
        Mask for row tokens, shape [B, N_row] or [B_mul, N_row].
        If [B, N_row], broadcasts to B_mul.
    mask_col : torch.Tensor
        Mask for column tokens, shape [B, N_col] or [B_mul, N_col].
        If [B, N_col], broadcasts to B_mul.
    multiplicity : int
        Required. Explicit multiplicity factor where B_mul = B * multiplicity.
    num_bins : int, optional
        Number of distance bins for PDE. Default: 64.
    max_dist : float, optional
        Maximum distance for binning. Default: 32.0.

    Returns
    -------
    out_loss_num : torch.Tensor
        Sum of cross-entropy loss over all (i,j) pairs per batch, shape [B_mul].
        out_loss_num[b] = sum_{i,j}(CE_loss[i,j] * mask[i,j])
    out_mask_denom : torch.Tensor
        Sum of mask over all (i,j) pairs per batch, shape [B_mul].
        out_mask_denom[b] = sum_{i,j}(mask[i,j])

    Notes
    -----
    To compute the final normalized PDE loss:
        loss = out_loss_num / (eps + out_mask_denom)

    For distributed training, allreduce out_loss_num and out_mask_denom separately
    before computing the normalized loss.
    """
    return _CdistPDEImpl.apply(
        pred_pde,
        true_coords_row,
        true_coords_col,
        pred_coords_row,
        pred_coords_col,
        mask_row,
        mask_col,
        multiplicity,
        num_bins,
        max_dist,
    )
