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


from enum import Enum, auto
from typing import Tuple

import torch
from torch import nn
from torch.distributed.device_mesh import DeviceMesh
from torch.distributed.tensor import DTensor, Shard

from boltz.distributed.comm import Ring2DComm
from boltz.distributed.model.layers.layernorm import LayerNormParamsReplicated
from boltz.distributed.model.layers.linear import LinearParamsReplicated
from boltz.distributed.model.layers.sigmoid_gate import sigmoid_gate
from boltz.distributed.utils import update_exhaustive_strides
from boltz.model.layers.triangular_mult import (
    TriangleMultiplicationIncoming as SerialTriangleMultiplicationIncoming,
)
from boltz.model.layers.triangular_mult import (
    TriangleMultiplicationOutgoing as SerialTriangleMultiplicationOutgoing,
)


class _XposeArgs(Enum):
    lhs = auto()
    rhs = auto()


def _distributed_bmm(
    lhs: torch.Tensor,
    rhs: torch.Tensor,
    comm: Ring2DComm,
    permute_lhs: tuple[int, ...] | None = None,
    permute_rhs: tuple[int, ...] | None = None,
    permute_out: tuple[int, ...] | None = None,
    xpose_args: _XposeArgs | None = None,
) -> torch.Tensor:
    """Perform distributed batch matrix multiplication using ring communication.

    This function implements a memory-efficient distributed batch matrix
    multiply operation across a 2D process grid using ring communication patterns.
    It computes the matrix multiplication of two tensors while minimizing memory usage
    through double buffering and overlapping computation with communication.

    The algorithm works by:
    1. Optionally permuting input tensors to desired layouts
    2. Setting up communication buffers based on transpose requirements
    3. Using ring communication to rotate tensor chunks across processes
    4. Performing overlapped computation and communication with double buffering
    5. Accumulating partial results to compute the final distributed bmm

    Communication Patterns
    ----------------------
    The function uses Ring2DComm to implement sophisticated communication patterns
    across a 2D process grid. Below are ASCII diagrams illustrating the key phases:

    **Phase 1: Initial 2D Grid Setup**

    For a 3x3 process grid, each process (i,j) initially owns tensor chunks:
    ```
    ┌─────┬─────┬─────┐
    │(0,0)│(0,1)│(0,2)│  ← Row 0
    ├─────┼─────┼─────┤
    │(1,0)│(1,1)│(1,2)│  ← Row 1
    ├─────┼─────┼─────┤
    │(2,0)│(2,1)│(2,2)│  ← Row 2
    └─────┴─────┴─────┘
      ↑     ↑     ↑
     Col 0 Col 1 Col 2
    ```

    **Phase 2: Transpose Communication (if xpose_args specified)**

    e.g., when xpose_args=_XposeArgs.rhs, RHS tensor is transposed across the 2D grid:
    ```
    Original RHS Ownership    After Transpose Communication
    ┌─────┬─────┬─────┐      ┌─────┬─────┬─────┐
    │ R00 │ R01 │ R02 │      │ R00 │ R10 │ R20 │
    ├─────┼─────┼─────┤  →   ├─────┼─────┼─────┤
    │ R10 │ R11 │ R12 │      │ R01 │ R11 │ R21 │
    ├─────┼─────┼─────┤      ├─────┼─────┼─────┤
    │ R20 │ R21 │ R22 │      │ R02 │ R12 │ R22 │
    └─────┴─────┴─────┘      └─────┴─────┴─────┘
    ```
    When xpose_args=_XposeArgs.lhs, LHS tensor is similarly transposed across the 2D grid

    **Phase 3: Initial Ring Setup**

    Row initialization (comm_row_init): Each row i shifts left by i positions
    ```
    Before Row Init              After Row Init
    ┌─────┬─────┬─────┐         ┌─────┬─────┬─────┐
    │ L00 │ L01 │ L02 │ ←shift 0│ L00 │ L01 │ L02 │
    ├─────┼─────┼─────┤         ├─────┼─────┼─────┤
    │ L10 │ L11 │ L12 │ ←shift 1│ L11 │ L12 │ L10 │
    ├─────┼─────┼─────┤         ├─────┼─────┼─────┤
    │ L20 │ L21 │ L22 │ ←shift 2│ L22 │ L20 │ L21 │
    └─────┴─────┴─────┘         └─────┴─────┴─────┘
    ```

    Column initialization (comm_col_init): Each column j shifts up by j positions
    ```
    Before Col Init              After Col Init
    ┌─────┬─────┬─────┐         ┌─────┬─────┬─────┐
    │ R00 │ R01 │ R02 │         │ R00 │ R11 │ R22 │
    ├─────┼─────┼─────┤  shift  ├─────┼─────┼─────┤
    │ R10 │ R11 │ R12 │   ↑     │ R10 │ R21 │ R02 │
    ├─────┼─────┼─────┤   0,1,2 ├─────┼─────┼─────┤
    │ R20 │ R21 │ R22 │         │ R20 │ R01 │ R12 │
    └─────┴─────┴─────┘         └─────┴─────┴─────┘
    ```

    **Phase 4: Ring Computation Loop**

    For each iteration k in range(grid_size):
    1. Compute partial matmul: out += matmul(lhs_chunk, rhs_chunk)
    2. Ring shift both tensors for next iteration

    Ring communication pattern (each step shifts by 1):
    ```
    Step 0 → Step 1 → Step 2 (back to original)

    LHS Row Shifts (left by 1):
    ┌─────┬─────┬─────┐    ┌─────┬─────┬─────┐    ┌─────┬─────┬─────┐
    │ L00 │ L01 │ L02 │ →  │ L01 │ L02 │ L00 │ →  │ L02 │ L00 │ L01 │
    ├─────┼─────┼─────┤    ├─────┼─────┼─────┤    ├─────┼─────┼─────┤
    │ L11 │ L12 │ L10 │ →  │ L12 │ L10 │ L11 │ →  │ L10 │ L11 │ L12 │
    ├─────┼─────┼─────┤    ├─────┼─────┼─────┤    ├─────┼─────┼─────┤
    │ L22 │ L20 │ L21 │ →  │ L20 │ L21 │ L22 │ →  │ L21 │ L22 │ L20 │
    └─────┴─────┴─────┘    └─────┴─────┴─────┘    └─────┴─────┴─────┘

    RHS Column Shifts (up by 1):
    ┌─────┬─────┬─────┐    ┌─────┬─────┬─────┐    ┌─────┬─────┬─────┐
    │ R00 │ R11 │ R22 │    │ R10 │ R21 │ R02 │    │ R20 │ R01 │ R12 │
    ├─────┼─────┼─────┤    ├─────┼─────┼─────┤    ├─────┼─────┼─────┤
    │ R10 │ R21 │ R02 │ →  │ R20 │ R01 │ R12 │ →  │ R00 │ R11 │ R22 │
    ├─────┼─────┼─────┤    ├─────┼─────┼─────┤    ├─────┼─────┼─────┤
    │ R20 │ R01 │ R12 │    │ R00 │ R11 │ R22 │    │ R10 │ R21 │ R02 │
    └─────┴─────┴─────┘    └─────┴─────┴─────┘    └─────┴─────┴─────┘
    ```

    **Double Buffering Strategy**

    The algorithm uses double buffering to overlap communication with computation:
    ```
    Time →  │ Compute │ Compute │ Compute │
            │ Buffer0 │ Buffer1 │ Buffer0 │
            │    ↓    │    ↓    │    ↓    │
    Comm →  │   Send  │   Send  │   Send  │
            │  Buffer1│ Buffer0 │ Buffer1 │
            │   Recv  │   Recv  │   Recv  │
            │ Buffer1 │ Buffer0 │ Buffer1 │
    ```

    This ensures that while one buffer is being used for computation, the other
    buffer is being prepared through communication for the next iteration.


    Parameters
    ----------
    lhs : torch.Tensor
        Left-hand side tensor for matrix multiplication.
        Typically has shape (B, ...) where B is batch dimension.
    rhs : torch.Tensor
        Right-hand side tensor for matrix multiplication.
        Must be compatible with lhs for matrix multiplication after permutations.
    comm : Ring2DComm
        Ring communication object configured for 2D process grid communication.
        Provides row and column communication groups for distributed computation.
    permute_lhs : tuple[int, ...] | None, optional
        Permutation indices to apply to lhs tensor before computation. Typically
        the permutation with group the batch-like axes into leading axes and reshape
        the last two axes into "N" and "K" dimensions (in the NMK notation)
        If None, no permutation is applied. Default is None.
    permute_rhs : tuple[int, ...] | None, optional
        Permutation indices to apply to rhs tensor before computation. Typically
        the permutation with group the batch-like axes into leading axes and reshape
        the last two axes into "K" and "M" dimensions (in the NMK notation)
        If None, no permutation is applied. Default is None.
    permute_out : tuple[int, ...] | None, optional
        Permutation indices to apply to output tensor after computation. Typically
        the permutation reverts the resulting permutation of the output matrix
        due to the permutation of the input lhs' and rhs' axes.
        If None, no permutation is applied. Default is None.
    xpose_args : _XposeArgs | None, optional
        Specifies which tensor requires transpose communication:
        - _XposeArgs.lhs: Transpose communication for left-hand side tensor
        - _XposeArgs.rhs: Transpose communication for right-hand side tensor
        - None: No transpose communication required
        Default is None.

    Returns
    -------
    torch.Tensor
        Result of the distributed batch matrix multiplication.
        Shape depends on input shapes and permutation arguments.

    Examples
    --------
    Typical usage in triangle multiplication:

    >>> # For outgoing triangle multiplication
    >>> result = _distributed_bmm(
    ...     lhs=tensor_a,
    ...     rhs=tensor_b,
    ...     comm=ring_comm,
    ...     permute_lhs=(0, 3, 1, 2),  # (B, n, k, D) -> (B, D, n, k)
    ...     permute_rhs=(0, 3, 2, 1),  # (B, m, k, D) -> (B, D, k, m)
    ...     permute_out=(0, 2, 3, 1),  # (B, D, n, m) -> (B, n, m, D)
    ...     xpose_args=_XposeArgs.rhs
    ... )
    """
    if permute_lhs is not None:
        lhs = lhs.permute(permute_lhs)
    # this enforces lhs and rhs to be a clone so that the in-place modification
    # does not affect the input tensor
    lhs = lhs.clone(memory_format=torch.contiguous_format)
    if permute_rhs is not None:
        rhs = rhs.permute(permute_rhs)
    rhs = rhs.clone(memory_format=torch.contiguous_format)

    if xpose_args == _XposeArgs.lhs:
        lhs_recv = comm.comm_2d_trans.enqueue_to_dispatch(lhs)
        rhs_recv = rhs
        rhs = torch.empty_like(rhs_recv)
    elif xpose_args == _XposeArgs.rhs:
        rhs_recv = comm.comm_2d_trans.enqueue_to_dispatch(rhs)
        lhs_recv = lhs
        lhs = torch.empty_like(lhs_recv)
    elif xpose_args is None:
        lhs_recv = lhs
        lhs = torch.empty_like(lhs_recv)
        rhs_recv = rhs
        rhs = torch.empty_like(rhs_recv)
    else:
        raise ValueError(f"Invalid xpose_args: {xpose_args}")

    # post the comm_2d_trans.wait_until_finished() (or no wait if xpose_args is not None),
    # *_recv are the correct tensors to operate on
    i_ready = 0
    i_recv = i_ready ^ 1
    lhs_buffer = [lhs_recv, lhs]
    rhs_buffer = [rhs_recv, rhs]

    if xpose_args is not None:
        comm.comm_2d_trans.wait_until_finished()

    lhs_buffer[i_recv] = comm.comm_row_init.enqueue_to_dispatch(lhs_buffer[i_ready], lhs_buffer[i_recv])
    rhs_buffer[i_recv] = comm.comm_col_init.enqueue_to_dispatch(rhs_buffer[i_ready], rhs_buffer[i_recv])

    i_ready ^= 1
    i_recv ^= 1

    out = torch.zeros_like(lhs_buffer[i_ready])

    comm.comm_row_init.wait_until_finished()
    comm.comm_col_init.wait_until_finished()

    # Double buffering computation
    for k_step in range(comm.group_layout.shape[1]):
        lhs_ready = lhs_buffer[i_ready]
        rhs_ready = rhs_buffer[i_ready]
        if k_step < comm.group_layout.shape[1] - 1:
            lhs_buffer[i_recv] = comm.comm_row.enqueue_to_dispatch(lhs_ready, lhs_buffer[i_recv])
            rhs_buffer[i_recv] = comm.comm_col.enqueue_to_dispatch(rhs_ready, rhs_buffer[i_recv])
        out = out + torch.matmul(lhs_ready, rhs_ready)
        if k_step < comm.group_layout.shape[1] - 1:
            comm.comm_row.wait_until_finished()
            comm.comm_col.wait_until_finished()
            i_ready = i_ready ^ 1
            i_recv = i_recv ^ 1

    if permute_out is not None:
        out = out.permute(permute_out)
    return out


class _Direction(Enum):
    Outgoing = auto()
    Incoming = auto()


class _TriangleMultiplicationImpl(torch.autograd.Function):
    """Distributed implementation of triangle multiplication using ring communication.

    This autograd function implements a memory-efficient distributed triangle multiplication
    operation across a 2D process grid. The computation is parallelized using ring
    communication patterns to reduce memory usage and communication overhead.

    The triangle multiplication computes:

    for Outgoing:
        o = torch.einsum("bnkd,bmkd->bnmd", a * mask, b * mask)

    for Incoming:
        o = torch.einsum("bknd,bkmd->bnmd", a * mask, b * mask)

    Key features:
    - Distributed across a 2D grid with sharding on token dimensions (dim 1 and 2)
    - Uses ring communication to rotate data chunks during computation
    - Memory-efficient implementation that avoids materializing full tensors
    - Supports gradient computation through custom backward pass

    Notes
    -----
    Input tensors must be DTensors with:
    - Shape: (B, N_token1, N_token2, c_hidden) for tensors a and b
    - Shape: (B, N_token1, N_token2, 1) for mask tensor
    - Sharding on dimensions 1 and 2 (Shard(1) and Shard(2) placements)
    - Identical device mesh and placements across all inputs

    The algorithm uses a ring-based communication pattern where:
    - Tensor b is transposed and rotated by row
    - Tensor a is rotated by column
    - Each process computes partial matrix products and accumulates results
    """

    @staticmethod
    @torch.amp.custom_fwd(device_type="cuda")
    def forward(ctx, x: DTensor, mask: DTensor, g: DTensor, comm: Ring2DComm, direction: _Direction) -> DTensor:
        """Forward pass of distributed triangle multiplication computation.

        Parameters
        ----------
        ctx : torch.autograd.function.BackwardCFrame
            Context object for saving information needed in backward pass.
        x : DTensor
            Input tensor with shape (B, N_token1, N_token2, c_hidden * 2).
            Must be sharded on dimensions 1 and 2.
        mask : DTensor
            Mask tensor with shape (B, N_token1, N_token2) indicating valid positions.
            Must be sharded on dimensions 1 and 2.
        g : DTensor
            pre-sigmoid gate tensor with shape (B, N_token1, N_token2, c_hidden * 2) indicating valid positions.
            Must be sharded on dimensions 1 and 2.
        comm : Ring2DComm
            Ring communication object configured for the distributed computation.
        direction : _Direction
            Direction of the triangle multiplication, Outgoing or Incoming.

        Returns
        -------
        DTensor
            Output tensor with shape (B, N_token1, N_token2, c_hidden).
            Contains the distributed triangle multiplication result.
        """
        # Check if inputs are of type DTensor
        if not isinstance(x, DTensor):
            raise TypeError(f"Input 'x' must be of type DTensor. Got type {type(x)}.")
        if not isinstance(mask, DTensor):
            raise TypeError(f"Input 'mask' must be of type DTensor. Got type {type(mask)}.")
        if not isinstance(g, DTensor):
            raise TypeError(f"Input 'g' must be of type DTensor. Got type {type(g)}.")

        # Check if inputs have identical device mesh
        device_mesh_input = x.device_mesh
        if device_mesh_input != mask.device_mesh:
            raise ValueError(
                f"Input tensors 'x' and 'mask' must have identical device mesh. "
                f"Got device meshes {device_mesh_input} and {mask.device_mesh}."
            )
        if device_mesh_input != g.device_mesh:
            raise ValueError(
                f"Input tensors 'x' and 'g' must have identical device mesh. "
                f"Got device meshes {device_mesh_input} and {g.device_mesh}."
            )

        # Check if inputs have identical placements
        placements_input = x.placements
        if placements_input != mask.placements:
            raise ValueError(
                f"Input tensors 'x' and 'mask' must have identical placements. "
                f"Got placements {placements_input} and {mask.placements}."
            )
        if placements_input != g.placements:
            raise ValueError(
                f"Input tensors 'x' and 'g' must have identical placements. "
                f"Got placements {placements_input} and {g.placements}."
            )
        if placements_input != (Shard(0), Shard(1), Shard(2)):
            # For debugging, we requires the placements to be (Shard(0), Shard(1), Shard(2))
            # TODO: remove this to only use the previous check
            raise ValueError(
                f"Input tensor 'x's placements are not (Shard(0), Shard(1), Shard(2)). "
                f"Got placements {placements_input}."
            )

        # Check input shapes
        if x.shape[-1] % 2 != 0:
            raise ValueError(f"Input tensor 'x' must have an even number of hidden dimension size. Got {x.shape[-1]}")

        if x.ndim != 4:
            raise ValueError(f"Input tensor 'x' must have 4 dimensions. Got {x.ndim} dimensions.")

        if mask.ndim != 3:
            raise ValueError(f"Input tensor 'mask' must have 3 dimensions. Got {mask.ndim} dimensions.")

        if mask.shape != x.shape[:3]:
            raise ValueError(
                f"Input tensor 'mask' must have the same shape as the first 3 dimensions of 'x'. "
                f"Got mask shape: {mask.shape} vs x shape[:3]: {x.shape[:3]}"
            )
        if g.shape != x.shape:
            raise ValueError(
                f"Input tensor 'g' must have the same shape as 'x'. Got g shape: {g.shape} vs x shape: {x.shape}"
            )

        # Perform consistency check between the ring_comm and the device_mesh_input
        i_tensor_dim_to_i_grid_axis = [-1] * x.ndim
        for i_grid_axis, placement in enumerate(placements_input):
            if isinstance(placement, Shard):
                i_tensor_dim_to_i_grid_axis[placement.dim] = i_grid_axis
        if i_tensor_dim_to_i_grid_axis[1] == -1 or i_tensor_dim_to_i_grid_axis[2] == -1:
            raise ValueError(f"Input tensors' dimensions 1 and 2 must be sharded. Got placements {placements_input}.")

        # Check ring_comm consistency
        if comm.group_col != device_mesh_input.get_group(i_tensor_dim_to_i_grid_axis[1]):
            raise ValueError(
                "Input ring_comm's group_col process group is not the same as the group sharding the input tensors' axis 1"
            )

        coord_device_mesh_input = device_mesh_input.get_coordinate()
        if coord_device_mesh_input is None:
            raise ValueError(f"ring_comm.coord_2d {comm.coord_2d} is not on device_mesh_input {device_mesh_input}.")
        if comm.coord_2d != (
            coord_device_mesh_input[i_tensor_dim_to_i_grid_axis[1]],
            coord_device_mesh_input[i_tensor_dim_to_i_grid_axis[2]],
        ):
            raise ValueError(
                f"Input ring_comm's coord_2d {comm.coord_2d} does not match the "
                f"device mesh's rank coordinates {coord_device_mesh_input} for the sharded dimensions."
            )

        ctx.mark_non_differentiable(mask)

        # Apply mask and prepare for computation
        mask_local = mask.to_local().unsqueeze(-1)
        g_local = g.to_local().sigmoid()
        x_local = x.to_local() * mask_local
        x_local *= g_local

        # the _distributed_bmm will permute a_local and b_local and make
        # the resulting tensors contiguous so we don't need to clone them here
        a_local, b_local = torch.chunk(x_local, 2, dim=-1)

        # Store tensors for backward pass
        if x.requires_grad:
            # here x_local is masked and gated
            ctx.save_for_backward(a_local, b_local, mask_local, x_local, g_local)
            ctx.comm = comm
            ctx.shape_x_input = x.shape
            ctx.stride_x_input = x.stride()
            ctx.shape_g_input = g.shape
            ctx.stride_g_input = g.stride()
            ctx.placements_input = placements_input
            ctx.device_mesh_input = device_mesh_input
            ctx.direction = direction

        if direction == _Direction.Outgoing:
            permute_lhs = (0, 3, 1, 2)  # from (B, n, k, D) to (B, D, n, k)
            permute_rhs = (0, 3, 2, 1)  # from (B, m, k, D) to (B, D, k, m)
            permute_out = (0, 2, 3, 1)  # from (B, D, n, m) to (B, n, m, D)
            xpose_args = _XposeArgs.rhs
        elif direction == _Direction.Incoming:
            permute_lhs = (0, 3, 2, 1)  # from (B, k, n, D) to (B, D, n, k)
            permute_rhs = (0, 3, 1, 2)  # from (B, k, m, D) to (B, D, k, m)
            permute_out = (0, 2, 3, 1)  # from (B, D, n, m) to (B, n, m, D)
            xpose_args = _XposeArgs.lhs
        else:
            raise ValueError(f"Invalid direction: {direction}")

        out_local = _distributed_bmm(
            a_local,
            b_local,
            comm,
            permute_lhs=permute_lhs,
            permute_rhs=permute_rhs,
            permute_out=permute_out,
            xpose_args=xpose_args,
        ).contiguous()

        shape_output = x.shape[:-1] + (out_local.shape[-1],)
        stride_output = update_exhaustive_strides(x.shape, x.stride(), shape_output)
        # Convert back to DTensor
        out = DTensor.from_local(
            out_local,
            device_mesh=device_mesh_input,
            placements=placements_input,
            shape=shape_output,
            stride=stride_output,
        )
        return out

    @staticmethod
    @torch.amp.custom_bwd(device_type="cuda")
    def backward(ctx, d_loss_d_out: DTensor) -> Tuple[DTensor, None, DTensor, None, None, None]:
        """Backward pass of distributed triangle multiplication computation."""
        if not isinstance(d_loss_d_out, DTensor):
            raise TypeError(f"Input 'd_loss_d_out' must be of type DTensor. Got type {type(d_loss_d_out)}.")

        if d_loss_d_out.device_mesh != ctx.device_mesh_input:
            raise ValueError(
                f"Input 'd_loss_d_out' must have the same device mesh as the input tensors. "
                f"Got device meshes {d_loss_d_out.device_mesh} and {ctx.device_mesh_input}."
            )

        if d_loss_d_out.placements != ctx.placements_input:
            raise ValueError(
                f"Input 'd_loss_d_out' must have the same placements as the input tensors. "
                f"Got placements {d_loss_d_out.placements} and {ctx.placements_input}."
            )

        a, b, mask_local, x_masked_gated_local, g_local = ctx.saved_tensors
        comm = ctx.comm
        direction = ctx.direction

        # cast d_loss_d_out to the same dtype as a (saved tensor) to avoid type promotion to FP32
        # Note: torch.amp.custom_bwd disables autocast, so operations run in the input dtype.
        # If the upstream adjoint (d_loss_d_out) arrives as FP32 (e.g. from loss scaling or downstream FP32 layers),
        # mixed-precision ops with saved BF16 tensors would promote to FP32, causing potential communication
        # buffer mismatches and NCCL hangs. Explicit casting ensures consistent precision.
        d_loss_d_out_local = d_loss_d_out.to_local().to(dtype=a.dtype)

        if direction == _Direction.Outgoing:
            lhs_da = d_loss_d_out_local
            rhs_da = b
            permute_lhs_da = (0, 3, 1, 2)  # from (B, n, m, D) to (B, D, n, m)
            permute_rhs_da = (0, 3, 1, 2)  # from (B, m, k, D) to (B, D, m, k)
            permute_out_da = (0, 2, 3, 1)  # from (B, D, n, k) to (B, n, k, D)
            xpose_args_da = None

            lhs_db = d_loss_d_out_local
            rhs_db = a
            permute_lhs_db = (0, 3, 2, 1)  # from (B, n, m, D) to (B, D, m, n)
            permute_rhs_db = (0, 3, 1, 2)  # from (B, n, k, D) to (B, D, n, k)
            permute_out_db = (0, 2, 3, 1)  # from (B, D, m, k) to (B, m, k, D)
            xpose_args_db = _XposeArgs.lhs

        elif direction == _Direction.Incoming:
            lhs_da = b
            rhs_da = d_loss_d_out_local
            permute_lhs_da = (0, 3, 1, 2)  # from (B, k, m, D) to (B, D, k, m)
            permute_rhs_da = (0, 3, 2, 1)  # from (B, n, m, D) to (B, D, m, n)
            permute_out_da = (0, 2, 3, 1)  # from (B, D, k, n) to (B, k, n, D)
            xpose_args_da = _XposeArgs.rhs

            lhs_db = a
            rhs_db = d_loss_d_out_local
            permute_lhs_db = (0, 3, 1, 2)  # from (B, k, n, D) to (B, D, k, n)
            permute_rhs_db = (0, 3, 1, 2)  # from (B, n, m, D) to (B, D, n, m)
            permute_out_db = (0, 2, 3, 1)  # from (B, D, k, m) to (B, k, m, D)
            xpose_args_db = None
        else:
            raise ValueError(f"Invalid direction: {direction}")

        d_loss_d_a_local = _distributed_bmm(
            lhs_da,
            rhs_da,
            comm,
            permute_lhs=permute_lhs_da,
            permute_rhs=permute_rhs_da,
            permute_out=permute_out_da,
            xpose_args=xpose_args_da,
        ).contiguous()

        # Phase 2: d_loss_d_b
        d_loss_d_b_local = _distributed_bmm(
            lhs_db,
            rhs_db,
            comm,
            permute_lhs=permute_lhs_db,
            permute_rhs=permute_rhs_db,
            permute_out=permute_out_db,
            xpose_args=xpose_args_db,
        ).contiguous()

        # concatenate and apply mask to gradients
        dab_local = torch.cat([d_loss_d_a_local, d_loss_d_b_local], dim=-1)

        x_masked_gated_local *= 1 - g_local
        dg_local = dab_local * x_masked_gated_local

        dg = DTensor.from_local(
            dg_local,
            device_mesh=ctx.device_mesh_input,
            placements=ctx.placements_input,
            shape=ctx.shape_g_input,
            stride=ctx.stride_g_input,
        )

        dx_local = dab_local
        dx_local *= mask_local
        dx_local *= g_local

        # Convert gradients back to DTensors
        dx = DTensor.from_local(
            dx_local,
            device_mesh=ctx.device_mesh_input,
            placements=ctx.placements_input,
            shape=ctx.shape_x_input,
            stride=ctx.stride_x_input,
        )

        return dx, None, dg, None, None


class TriangleMultiplication(nn.Module):
    """Distributed triangle multiplication layer.

    This layer implements a distributed version of the triangle multiplication operation,
    which is used in attention mechanisms for protein structure prediction and other applications
    requiring pairwise feature interactions.

    The layer performs the following operations:
    1. Layer normalization of input pairwise features
    2. Linear projections to create two representation streams (a and b)
    3. Distributed triangle multiplication computation using ring communication
    4. Output gating and final linear projection

    Parameters
    ----------
    layer : SerialTriangleMultiplicationOutgoing | SerialTriangleMultiplicationIncoming
        The serial triangle multiplication layer to convert to distributed version.
        Used to initialize projection weights and normalization parameters.
    device_mesh : DeviceMesh
        The device mesh for distributed computation across multiple GPUs.
    comm : Ring2DComm
        Ring communication object for efficient distributed triangle multiplication computation.
    """

    def __init__(
        self,
        direction: _Direction,
        layer: SerialTriangleMultiplicationOutgoing | SerialTriangleMultiplicationIncoming,
        device_mesh: DeviceMesh,
        comm: Ring2DComm,
    ) -> None:
        """Initialize the distributed triangle multiplication layer."""
        super().__init__()
        self.device_mesh = device_mesh
        self.ring_comm = comm

        self.norm_in = LayerNormParamsReplicated(layer.norm_in, self.device_mesh)
        self.p_in = LinearParamsReplicated(layer.p_in, self.device_mesh)
        self.g_in = LinearParamsReplicated(layer.g_in, self.device_mesh)

        self.norm_out = LayerNormParamsReplicated(layer.norm_out, self.device_mesh)
        self.p_out = LinearParamsReplicated(layer.p_out, self.device_mesh)
        self.g_out = LinearParamsReplicated(layer.g_out, self.device_mesh)

        if direction == _Direction.Outgoing:
            if not isinstance(layer, SerialTriangleMultiplicationOutgoing):
                raise ValueError(f"Invalid layer type for direction {direction}: {type(layer)}")
        elif direction == _Direction.Incoming:
            if not isinstance(layer, SerialTriangleMultiplicationIncoming):
                raise ValueError(f"Invalid layer type for direction {direction}: {type(layer)}")
        else:
            raise ValueError(f"Invalid direction {direction}")
        self._direction = direction

    def forward(self, x: DTensor, mask: DTensor) -> DTensor:
        """Forward pass of the distributed triangle multiplication layer.

        Parameters
        ----------
        x : DTensor
            Input pairwise tensor with shape (B, N, N, D).
            Must be sharded on dimensions 1 and 2.
        mask : DTensor
            Mask tensor with shape (B, N, N) indicating valid positions.
            Must be sharded on dimensions 1 and 2.

        Returns
        -------
        DTensor
            Output pairwise tensor with shape (B, N, N, D).
        """
        # Stabilize pair embedding tensor with layer norm
        x = self.norm_in(x)
        x_in = x
        g_out = self.g_out(x_in)

        # Decompress: D -> 2D
        g = self.g_in(x)
        x = self.p_in(x)

        # Distributed triangular multiplication (mask is applied inside the implementation)
        x = _TriangleMultiplicationImpl.apply(x, mask, g, self.ring_comm, self._direction)

        # Output gating
        x = self.p_out(self.norm_out(x))
        x = sigmoid_gate(x, g_out)

        return x


class TriangleMultiplicationOutgoing(TriangleMultiplication):
    """Distributed triangle multiplication outgoing layer."""

    def __init__(
        self,
        layer: SerialTriangleMultiplicationOutgoing,
        device_mesh: DeviceMesh,
        comm: Ring2DComm,
    ) -> None:
        """Initialize the distributed triangle multiplication outgoing layer.

        Parameters
        ----------
        layer : SerialTriangleMultiplicationOutgoing
            The serial triangle multiplication outgoing layer to convert to distributed version.
        device_mesh : DeviceMesh
            The device mesh for distributed computation across multiple GPUs.
        comm : Ring2DComm
            Ring communication object for efficient distributed triangle multiplication computation.
        """
        super().__init__(_Direction.Outgoing, layer, device_mesh, comm)


class TriangleMultiplicationIncoming(TriangleMultiplication):
    """Distributed triangle multiplication incoming layer."""

    def __init__(
        self,
        layer: SerialTriangleMultiplicationIncoming,
        device_mesh: DeviceMesh,
        comm: Ring2DComm,
    ) -> None:
        """Initialize the distributed triangle multiplication incoming layer.

        Parameters
        ----------
        layer : SerialTriangleMultiplicationIncoming
            The serial triangle multiplication incoming layer to convert to distributed version.
        device_mesh : DeviceMesh
            The device mesh for distributed computation across multiple GPUs.
        comm : Ring2DComm
            Ring communication object for efficient distributed triangle multiplication computation.
        """
        super().__init__(_Direction.Incoming, layer, device_mesh, comm)
