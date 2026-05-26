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


try:
    import triton
    import triton.language as tl
except ImportError:
    raise ImportError("Triton is not available. Will not import smooth_lddt_loss kernels.")


def grid_launch_config(kwargs):
    return (
        kwargs["shape_coords_axis_0"],
        triton.cdiv(kwargs["shape_coords_axis_1"], kwargs["BLOCK"]),
        triton.cdiv(kwargs["shape_coords_axis_1"], kwargs["BLOCK"]),
    )


@triton.heuristics(
    # register pressure is high with BLOCK >= 64 and will spill at BLOCK = 128
    # due to the pair-wise distance computation
    {"BLOCK": lambda args: 32, "num_warps": lambda args: 4},
)
@triton.jit
def smooth_lddt_loss_fwd_kernel(
    pred_coords_ptr,
    true_coords_ptr,
    pred_coords_t_ptr,
    true_coords_t_ptr,
    is_nucleotide_ptr,
    coords_mask_ptr,
    coords_mask_t_ptr,
    num_output_ptr,
    den_output_ptr,
    stride_pred_b,
    stride_pred_n,
    stride_pred_d,
    stride_true_b,
    stride_true_n,
    stride_true_d,
    stride_pred_t_b,
    stride_pred_t_n,
    stride_pred_t_d,
    stride_true_t_b,
    stride_true_t_n,
    stride_true_t_d,
    stride_nuc_b,
    stride_nuc_n,
    stride_mask_b,
    stride_mask_n,
    stride_mask_t_b,
    stride_mask_t_n,
    nucleic_acid_cutoff,
    other_cutoff,
    is_self_comm: tl.constexpr,
    shape_coords_axis_0,
    shape_coords_axis_1,
    shape_mask_axis_0,
    BLOCK: tl.constexpr,
):
    pid = tl.program_id(0)
    pid_m = tl.program_id(1)
    pid_n = tl.program_id(2)

    # Batch offset
    batch_idx = pid

    multiplicity = shape_coords_axis_0 // shape_mask_axis_0

    # we can reuse the same mask batch without pre-repeat_interleave
    # the mask with multiplicity
    batch_idx_mask = batch_idx // multiplicity

    # Pointers to current batch
    # Compute pointers by adding offsets * strides
    # Note: pointers are 64-bit

    pred_coords_cur = pred_coords_ptr + batch_idx * stride_pred_b
    true_coords_cur = true_coords_ptr + batch_idx * stride_true_b
    pred_coords_t_cur = pred_coords_t_ptr + batch_idx * stride_pred_t_b
    true_coords_t_cur = true_coords_t_ptr + batch_idx * stride_true_t_b

    is_nucleotide_cur = is_nucleotide_ptr + batch_idx_mask * stride_nuc_b
    coords_mask_cur = coords_mask_ptr + batch_idx_mask * stride_mask_b
    coords_mask_t_cur = coords_mask_t_ptr + batch_idx_mask * stride_mask_t_b

    # Offsets for M dimension
    offs_m = pid_m * BLOCK + tl.arange(0, BLOCK)
    mask_m = offs_m < shape_coords_axis_1

    # Offsets for shape_coords_axis_1 dimension
    offs_n = pid_n * BLOCK + tl.arange(0, BLOCK)
    mask_n = offs_n < shape_coords_axis_1

    # Load data
    # is_nucleotide: (BLOCK_M,)
    # Use float mask for simpler math
    is_nuc = tl.load(is_nucleotide_cur + offs_m * stride_nuc_n, mask=mask_m, other=0.0).to(tl.int1)

    # coords_mask: (BLOCK_M,)
    c_mask = tl.load(coords_mask_cur + offs_m * stride_mask_n, mask=mask_m, other=0.0)

    # coords_mask_t: (BLOCK_N,)
    c_mask_t = tl.load(coords_mask_t_cur + offs_n * stride_mask_t_n, mask=mask_n, other=0.0)

    # Combined mask (BLOCK_M, BLOCK_N)
    # mask = c_mask[:, None] * c_mask_t[None, :]
    combined_mask = c_mask[:, None] * c_mask_t[None, :]

    # Handle diagonal
    is_diag = (offs_m[:, None] == offs_n[None, :]) & is_self_comm
    combined_mask = tl.where(is_diag, 0.0, combined_mask)

    # Distances and Differences
    d_idx = tl.arange(0, 4)
    mask_d = d_idx < 3

    # True Coords
    tc_ptr = true_coords_cur + offs_m[:, None] * stride_true_n + d_idx[None, :] * stride_true_d
    tc_m = tl.load(tc_ptr, mask=mask_m[:, None] & mask_d[None, :], other=0.0)

    tc_t_ptr = true_coords_t_cur + offs_n[:, None] * stride_true_t_n + d_idx[None, :] * stride_true_t_d
    tc_t_n = tl.load(tc_t_ptr, mask=mask_n[:, None] & mask_d[None, :], other=0.0)

    diff_true = tc_m[:, None, :] - tc_t_n[None, :, :]
    true_dist_sq = tl.sum(diff_true * diff_true, axis=2)

    # Pred Coords
    pc_ptr = pred_coords_cur + offs_m[:, None] * stride_pred_n + d_idx[None, :] * stride_pred_d
    pc_m = tl.load(pc_ptr, mask=mask_m[:, None] & mask_d[None, :], other=0.0)

    pc_t_ptr = pred_coords_t_cur + offs_n[:, None] * stride_pred_t_n + d_idx[None, :] * stride_pred_t_d
    pc_t_n = tl.load(pc_t_ptr, mask=mask_n[:, None] & mask_d[None, :], other=0.0)

    diff_pred = pc_m[:, None, :] - pc_t_n[None, :, :]
    pred_dist_sq = tl.sum(diff_pred * diff_pred, axis=2)

    true_dist = tl.sqrt(true_dist_sq)
    pred_dist = tl.sqrt(pred_dist_sq)

    # Cutoff mask
    # is_nuc is (BLOCK_M,), broadcast to (BLOCK_M, BLOCK_N)
    cutoff = tl.where(is_nuc[:, None], nucleic_acid_cutoff, other_cutoff)
    dist_mask = true_dist < cutoff

    final_mask = combined_mask * dist_mask

    # Epsilon
    dist_diff = tl.abs(true_dist - pred_dist)

    eps = tl.sigmoid(0.5 - dist_diff)
    eps += tl.sigmoid(1.0 - dist_diff)
    eps += tl.sigmoid(2.0 - dist_diff)
    eps += tl.sigmoid(4.0 - dist_diff)
    eps *= 0.25

    # Accumulate
    num_val = eps * final_mask
    den_val = final_mask

    # Sum within block
    # Accumulate in fp32 for precision
    num_sum = tl.sum(num_val.to(tl.float32))
    den_sum = tl.sum(den_val.to(tl.float32))

    # Atomic Add to global
    tl.atomic_add(num_output_ptr + batch_idx, num_sum.to(pred_coords_ptr.dtype.element_ty))
    tl.atomic_add(den_output_ptr + batch_idx, den_sum.to(pred_coords_ptr.dtype.element_ty))


@triton.heuristics(
    # register pressure is high with BLOCK >= 64 and will spill at BLOCK = 128
    # due to the pair-wise distance computation
    {"BLOCK": lambda args: 16, "num_warps": lambda args: 2},
)
@triton.jit
def smooth_lddt_loss_bwd_kernel(
    grad_num_reduced_ptr,
    grad_den_reduced_ptr,
    pred_coords_ptr,
    true_coords_ptr,
    pred_coords_t_ptr,
    true_coords_t_ptr,
    is_nucleotide_ptr,
    coords_mask_ptr,
    coords_mask_t_ptr,
    grad_pred_local_ptr,
    grad_pred_t_local_ptr,
    stride_pred_b,
    stride_pred_n,
    stride_pred_d,
    stride_true_b,
    stride_true_n,
    stride_true_d,
    stride_pred_t_b,
    stride_pred_t_n,
    stride_pred_t_d,
    stride_true_t_b,
    stride_true_t_n,
    stride_true_t_d,
    stride_nuc_b,
    stride_nuc_n,
    stride_mask_b,
    stride_mask_n,
    stride_mask_t_b,
    stride_mask_t_n,
    stride_grad_pred_b,
    stride_grad_pred_n,
    stride_grad_pred_d,
    stride_grad_pred_t_b,
    stride_grad_pred_t_n,
    stride_grad_pred_t_d,
    nucleic_acid_cutoff,
    other_cutoff,
    is_self_comm: tl.constexpr,
    shape_coords_axis_0,
    shape_coords_axis_1,
    shape_mask_axis_0,
    BLOCK: tl.constexpr,
):
    pid = tl.program_id(0)
    pid_m = tl.program_id(1)
    pid_n = tl.program_id(2)

    batch_idx = pid

    multiplicity = shape_coords_axis_0 // shape_mask_axis_0

    # we can reuse the same mask batch without pre-repeat_interleave
    # the mask with multiplicity
    batch_idx_mask = batch_idx // multiplicity

    # Pointers
    pred_coords_cur = pred_coords_ptr + batch_idx * stride_pred_b
    true_coords_cur = true_coords_ptr + batch_idx * stride_true_b
    pred_coords_t_cur = pred_coords_t_ptr + batch_idx * stride_pred_t_b
    true_coords_t_cur = true_coords_t_ptr + batch_idx * stride_true_t_b

    grad_pred_local_cur = grad_pred_local_ptr + batch_idx * stride_grad_pred_b
    grad_pred_t_local_cur = grad_pred_t_local_ptr + batch_idx * stride_grad_pred_t_b

    is_nucleotide_cur = is_nucleotide_ptr + batch_idx_mask * stride_nuc_b
    coords_mask_cur = coords_mask_ptr + batch_idx_mask * stride_mask_b
    coords_mask_t_cur = coords_mask_t_ptr + batch_idx_mask * stride_mask_t_b

    # Gradients are scalars per batch
    grad_num = tl.load(grad_num_reduced_ptr + batch_idx)

    # Offsets
    offs_m = pid_m * BLOCK + tl.arange(0, BLOCK)
    mask_m = offs_m < shape_coords_axis_1
    offs_n = pid_n * BLOCK + tl.arange(0, BLOCK)
    mask_n = offs_n < shape_coords_axis_1

    # --- Recompute Forward Pass Intermediates ---

    # Masks
    is_nuc = tl.load(is_nucleotide_cur + offs_m * stride_nuc_n, mask=mask_m, other=0.0).to(tl.int1)
    c_mask = tl.load(coords_mask_cur + offs_m * stride_mask_n, mask=mask_m, other=0.0)
    c_mask_t = tl.load(coords_mask_t_cur + offs_n * stride_mask_t_n, mask=mask_n, other=0.0)

    combined_mask = c_mask[:, None] * c_mask_t[None, :]
    is_diag = (offs_m[:, None] == offs_n[None, :]) & is_self_comm
    combined_mask = tl.where(is_diag, 0.0, combined_mask)

    # Distances and Differences
    d_idx = tl.arange(0, 4)
    mask_d = d_idx < 3

    # True Coords
    tc_ptr = true_coords_cur + offs_m[:, None] * stride_true_n + d_idx[None, :] * stride_true_d
    tc_m = tl.load(tc_ptr, mask=mask_m[:, None] & mask_d[None, :], other=0.0)

    tc_t_ptr = true_coords_t_cur + offs_n[:, None] * stride_true_t_n + d_idx[None, :] * stride_true_t_d
    tc_t_n = tl.load(tc_t_ptr, mask=mask_n[:, None] & mask_d[None, :], other=0.0)

    diff_true = tc_m[:, None, :] - tc_t_n[None, :, :]
    true_dist_sq = tl.sum(diff_true * diff_true, axis=2)

    # Pred Coords
    pc_ptr = pred_coords_cur + offs_m[:, None] * stride_pred_n + d_idx[None, :] * stride_pred_d
    pc_m = tl.load(pc_ptr, mask=mask_m[:, None] & mask_d[None, :], other=0.0)

    pc_t_ptr = pred_coords_t_cur + offs_n[:, None] * stride_pred_t_n + d_idx[None, :] * stride_pred_t_d
    pc_t_n = tl.load(pc_t_ptr, mask=mask_n[:, None] & mask_d[None, :], other=0.0)

    diff_pred = pc_m[:, None, :] - pc_t_n[None, :, :]
    pred_dist_sq = tl.sum(diff_pred * diff_pred, axis=2)

    mask_0 = (d_idx == 0)[None, None, :]
    mask_1 = (d_idx == 1)[None, None, :]
    mask_2 = (d_idx == 2)[None, None, :]

    diff_pred_x = tl.sum(diff_pred * mask_0, axis=2)
    diff_pred_y = tl.sum(diff_pred * mask_1, axis=2)
    diff_pred_z = tl.sum(diff_pred * mask_2, axis=2)

    true_dist = tl.sqrt(true_dist_sq)
    pred_dist = tl.sqrt(pred_dist_sq)

    # Cutoff Mask
    cutoff = tl.where(is_nuc[:, None], nucleic_acid_cutoff, other_cutoff)
    dist_mask = true_dist < cutoff
    final_mask = combined_mask * dist_mask

    # --- Backward Computation ---

    # Compute d_eps_d_diff in fp32 for precision
    dist_diff = tl.abs(true_dist - pred_dist).to(tl.float32)
    d_eps_d_diff = tl.zeros([BLOCK, BLOCK], dtype=tl.float32)

    # Loop over cutoffs: 0.5, 1.0, 2.0, 4.0
    # Unrolling for simplicity/speed in Triton
    # 0.5
    val = 0.5 - dist_diff
    sig = tl.sigmoid(val)
    d_eps_d_diff -= sig * (1.0 - sig)
    # 1.0
    val = 1.0 - dist_diff
    sig = tl.sigmoid(val)
    d_eps_d_diff -= sig * (1.0 - sig)
    # 2.0
    val = 2.0 - dist_diff
    sig = tl.sigmoid(val)
    d_eps_d_diff -= sig * (1.0 - sig)
    # 4.0
    val = 4.0 - dist_diff
    sig = tl.sigmoid(val)
    d_eps_d_diff -= sig * (1.0 - sig)

    d_eps_d_diff *= 0.25

    # sign(pred - true)
    # if pred > true: 1, else -1. But be careful about 0.
    # Using tl.where
    # sign_diff = tl.where(pred_dist > true_dist, 1.0, -1.0)
    # Actually torch.sign returns 0 for 0.
    diff_dist = pred_dist - true_dist
    sign_diff = tl.where(diff_dist > 0, 1.0, tl.where(diff_dist < 0, -1.0, 0.0))

    # d_L_d_pred_dists
    # grad_num is scalar broadcasted
    d_L_d_pred_dists = grad_num * final_mask * d_eps_d_diff * sign_diff

    # diff_dir = diff_vec / (pred_dist + 1e-8)
    # pred_dist_safe
    pred_dist_safe = pred_dist + 1e-8
    inv_dist = 1.0 / pred_dist_safe

    # Compute grad_vec (shape_coords_axis_1^2, 3) effectively
    # factor = d_L_d_pred_dists * inv_dist
    factor = d_L_d_pred_dists * inv_dist

    d_L_d_diff_x = factor * diff_pred_x
    d_L_d_diff_y = factor * diff_pred_y
    d_L_d_diff_z = factor * diff_pred_z

    # Accumulate gradients locally
    # grad_pred_local (M, 3) = sum_over_N(d_L_d_diff)
    # grad_pred_t_local (shape_coords_axis_1, 3) = sum_over_M(-d_L_d_diff)

    # Sum over shape_coords_axis_1 (cols) for grad_pred_local
    grad_x_m = tl.sum(d_L_d_diff_x, axis=1)
    grad_y_m = tl.sum(d_L_d_diff_y, axis=1)
    grad_z_m = tl.sum(d_L_d_diff_z, axis=1)

    # Sum over M (rows) for grad_pred_t_local
    # Note: grad_pred_t_local is -sum
    grad_x_n = tl.sum(d_L_d_diff_x, axis=0)
    grad_y_n = tl.sum(d_L_d_diff_y, axis=0)
    grad_z_n = tl.sum(d_L_d_diff_z, axis=0)

    # Atomic Add to output buffers
    # grad_pred_local: (B, M, 3)
    dtype = grad_pred_local_ptr.dtype.element_ty
    tl.atomic_add(grad_pred_local_cur + offs_m * stride_grad_pred_n + 0, grad_x_m.to(dtype), mask=mask_m)
    tl.atomic_add(grad_pred_local_cur + offs_m * stride_grad_pred_n + 1, grad_y_m.to(dtype), mask=mask_m)
    tl.atomic_add(grad_pred_local_cur + offs_m * stride_grad_pred_n + 2, grad_z_m.to(dtype), mask=mask_m)

    # grad_pred_t_local: (B, shape_coords_axis_1, 3)
    # Negate
    tl.atomic_add(grad_pred_t_local_cur + offs_n * stride_grad_pred_t_n + 0, (-grad_x_n).to(dtype), mask=mask_n)
    tl.atomic_add(grad_pred_t_local_cur + offs_n * stride_grad_pred_t_n + 1, (-grad_y_n).to(dtype), mask=mask_n)
    tl.atomic_add(grad_pred_t_local_cur + offs_n * stride_grad_pred_t_n + 2, (-grad_z_n).to(dtype), mask=mask_n)
