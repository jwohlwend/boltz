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


from copy import deepcopy

import torch
import torch.distributed as dist
from torch.distributed.tensor import DTensor, Shard
from torch.distributed.tensor.device_mesh import DeviceMesh

from boltz.distributed.comm import One2OneComm, TransposeComm, ternary_parity
from boltz.distributed.model.layers.layernorm import LayerNormParamsReplicated
from boltz.distributed.model.layers.linear import LinearParamsReplicated
from boltz.distributed.utils import LayoutMap, get_group_rank_from_axial_shift, tiled_softmax_attention_update
from boltz.model.layers.pair_averaging import PairWeightedAveraging as SerialPairWeightedAveraging


class Ring2DCommPairAveraging:
    """
    Implements communication primitives for distributed operations on a 2D grid of devices.

    This class provides general-purpose ring communication patterns for operations like
    TriangleMultiplication and OuterProductMean across a 2D grid of devices. Unlike
    Ring2DCommTriAttn which is specialized for triangular attention, this class provides
    more general ring communication patterns.

    The communication patterns implemented include:
    1. Transpose communication for matrix operations
    2. Row-wise ring communication (left shifts)
    3. Column-wise ring communication (up shifts)

    Parameters
    ----------
    group_2d : dist.ProcessGroup
        The process group representing the 2D grid of devices. This should include
        all processes in the distributed computation.
    group_col : dist.ProcessGroup
        A subprocess group that provides communication between ranks in the same column.
    group_layout : LayoutMap
        A mapping from the 2D grid index to the flattened index of the devices on the 2D grid.
        Must represent a square grid (same dimensions in both axes).

    Notes
    -----
    The class implements various communication patterns needed for distributed matrix
    operations, including initial communication (with different shift patterns based on
    coordinates) and subsequent iterations (with fixed shifts).

    Communication ordering is carefully managed to prevent deadlocks by using
    ternary_parity to determine consistent send/receive ordering across different ranks.
    """

    def __init__(
        self,
        group_2d: dist.ProcessGroup,
        group_col: dist.ProcessGroup,
        group_layout: LayoutMap,
    ):
        """
        Ring comm over a 2d grid of devices with comm happening along both axes
        Arguments:
            group_2d: Group torch process group that provides communication
                across the full cross-device
            group_col: Subprocess group that provides communication
                between ranks in the same column
            group_layout: mapping from the 2d grid index to the flatten index
            of the devices on the 2d grid
        """
        # TODO: consolidate the ring 2d comm groups with other modules e,g. triangle attn
        self.group_2d = group_2d
        self.group_col = group_col
        self.group_layout = group_layout
        ranks_group_2d = set(dist.get_process_group_ranks(self.group_2d))
        ranks_group_col = set(dist.get_process_group_ranks(self.group_col))

        if not ranks_group_col.issubset(ranks_group_2d):
            raise ValueError("The col ranks are not a subset of ranks_group_2d")

        self.size_2d = dist.get_world_size(self.group_2d)

        if self.size_2d != self.group_layout.numel:
            raise ValueError(
                f"size of group_2d {self.size_2d} differs from the number of elements in group_layout {self.group_layout.numel}"
            )

        if self.group_layout.shape[0] != self.group_layout.shape[1]:
            raise ValueError(f"group_layout.shape {self.group_layout.shape} is not square")

        self.rank_2d = dist.get_rank(self.group_2d)
        self.coord_2d = self.group_layout.unravel(self.rank_2d)

        # all the send/recv ranks must be global in order to use isend/irecv
        # only need transpose at the beginning of the batched GEMM for b or a
        self.comm_2d_trans = TransposeComm(self.group_2d, self.group_layout)

        # always do left shift per row
        # for iteration 0, i'th row left shift by i column
        self.send_rank_row_init = get_group_rank_from_axial_shift(
            self.coord_2d, 1, -self.coord_2d[0], self.group_layout
        )
        self.recv_rank_row_init = get_group_rank_from_axial_shift(self.coord_2d, 1, self.coord_2d[0], self.group_layout)

        self.comm_row_init = One2OneComm(
            self.group_2d,
            self.send_rank_row_init,
            self.recv_rank_row_init,
            parity=ternary_parity(self.rank_2d, self.send_rank_row_init, self.recv_rank_row_init),
        )
        # for other iterations left shift by 1 col
        self.send_rank_row = get_group_rank_from_axial_shift(self.coord_2d, 1, -1, self.group_layout)
        self.recv_rank_row = get_group_rank_from_axial_shift(self.coord_2d, 1, 1, self.group_layout)

        self.comm_row = One2OneComm(
            self.group_2d,
            self.send_rank_row,
            self.recv_rank_row,
            parity=ternary_parity(self.rank_2d, self.send_rank_row, self.recv_rank_row),
        )

        # always do up shift per col
        # for iteration 0, j'th col up shift by j row
        self.send_rank_col_init = get_group_rank_from_axial_shift(
            self.coord_2d, 0, -self.coord_2d[1], self.group_layout
        )
        self.recv_rank_col_init = get_group_rank_from_axial_shift(self.coord_2d, 0, self.coord_2d[1], self.group_layout)
        self.comm_col_init = One2OneComm(
            self.group_2d,
            self.send_rank_col_init,
            self.recv_rank_col_init,
            parity=ternary_parity(self.rank_2d, self.send_rank_col_init, self.recv_rank_col_init),
        )
        # for other iterations, up shift by 1 row
        self.send_rank_col = get_group_rank_from_axial_shift(self.coord_2d, 0, -1, self.group_layout)
        self.recv_rank_col = get_group_rank_from_axial_shift(self.coord_2d, 0, 1, self.group_layout)
        self.comm_col = One2OneComm(
            self.group_2d,
            self.send_rank_col,
            self.recv_rank_col,
            parity=ternary_parity(self.rank_2d, self.send_rank_col, self.recv_rank_col),
        )

        self.comm_d_init = deepcopy(self.comm_row_init)
        self.comm_d = deepcopy(self.comm_row)

        self.comm_db = deepcopy(self.comm_col)

        # down shift j'th col by j row to reset db's data ownership as the input b
        self.send_rank_db_final = get_group_rank_from_axial_shift(self.coord_2d, 0, self.coord_2d[1], self.group_layout)
        self.recv_rank_db_final = get_group_rank_from_axial_shift(
            self.coord_2d, 0, -self.coord_2d[1], self.group_layout
        )
        self.comm_db_final = One2OneComm(
            self.group_2d,
            self.send_rank_db_final,
            self.recv_rank_db_final,
            parity=ternary_parity(self.rank_2d, self.send_rank_db_final, self.recv_rank_db_final),
        )

        self.comm_2d_trans_lse_m = deepcopy(self.comm_2d_trans)
        self.comm_2d_trans_amax = deepcopy(self.comm_2d_trans)


class _PairWeightedAveragingImpl(torch.autograd.Function):
    """Distributed implementation of pair weighted averaging using ring communication.

    This autograd function implements a memory-efficient distributed pair weighted averaging
    operation across a 2D process grid. The computation is parallelized using ring
    communication patterns to reduce memory usage and communication overhead.

    The pair weighted averaging operation computes:
    o[s,i] = sum_j(softmax(b[i,j]) * v[s,j])

    where the softmax is taken over the last dimension j, and is masked by the input mask.

    Key features:
    - Distributed across a 2D grid with sharding on sequence (dim 1) and token (dim 2) dimensions
    - Uses ring communication to rotate data chunks during computation
    - Memory-efficient implementation that avoids materializing full tensors
    - Supports gradient computation through custom backward pass

    Notes
    -----
    Input tensors must be DTensors with:
    - Shape: (B, H, S, N, c_h) for tensor v (value)
    - Shape: (B, H, N, N) for tensor b (bias/attention weights)
    - Shape: (B, N, N) for mask tensor
    - Sharding on dimensions 2 and 3 (Shard(2) and Shard(3) placements)
    - Identical device mesh and placements across all inputs

    The algorithm uses a ring-based communication pattern where:
    - Tensor v is rotated by row
    - Tensor b is rotated by column
    - Each process computes partial weighted averages and accumulates results
    """

    @staticmethod
    @torch.amp.custom_fwd(device_type="cuda")
    def forward(
        ctx, v: DTensor, b: DTensor, mask: DTensor, g: DTensor, comm: Ring2DCommPairAveraging, n_heads: int, inf: float
    ) -> DTensor:
        """Forward pass of distributed pair weighted averaging computation.

        Computes the pair weighted averaging of input tensors v and b using distributed
        ring communication to minimize memory usage and communication overhead.

        Parameters
        ----------
        ctx : torch.autograd.function.BackwardCFrame
            Context object for saving information needed in backward pass.
        v : DTensor
            Value tensor with shape (B, S, N, h_c_h).
            Must be sharded on dimensions 1 and 2.
        b : DTensor
            Bias/attention weights tensor with shape (B, N, N, n_heads).
            Must have compatible device mesh and placements and sharded on dimensions 1 and 2.
        mask : DTensor
            Mask tensor with shape (B, N, N).
            Must have same device mesh and placements as input tensors and sharded on dimensions 1 and 2.
        g : DTensor
            pre-sigmoid gate tensor with shape (B, S, N, h_c_h).
            Must have same device mesh and placements as input tensors and sharded on dimensions 1 and 2.
        comm : Ring2DCommPairAveraging
            Ring communication object configured for the distributed computation.
        n_heads : int
            Number of heads. The input tensor v's last dimension is h_c_h = n_heads * c_h where
            c_h will be derived by h_c_h // n_heads

        Returns
        -------
        DTensor
            gated output tensor with shape (B, S, N, h_c_h).
            Contains the distributed pair weighted averaging result ready for the output projection.

        Raises
        ------
        TypeError
            If inputs are not DTensor type.
        ValueError
            If tensor shapes, device meshes, or placements are incompatible,
            or if ring communication setup is inconsistent.
        """

        if not isinstance(comm, Ring2DCommPairAveraging):
            raise ValueError(f"Input comm must be of type Ring2DCommPairAveraging. Got type {type(comm)}.")

        # Check if inputs are of type DTensor
        if (
            not isinstance(v, DTensor)
            or not isinstance(b, DTensor)
            or not isinstance(mask, DTensor)
            or not isinstance(g, DTensor)
        ):
            raise TypeError(
                f"Inputs 'v', 'b', 'mask' and 'g' must be of type DTensor. Got types {type(v)}, {type(b)}, {type(mask)}, and {type(g)}."
            )

        # Check if inputs have identical device mesh
        device_mesh_input = v.device_mesh
        if device_mesh_input != b.device_mesh:
            raise ValueError(
                f"Input tensors 'v' and 'b' must have identical device mesh. "
                f"Got device meshes {device_mesh_input} and {b.device_mesh}."
            )
        if device_mesh_input != mask.device_mesh:
            raise ValueError(
                f"Input tensor 'mask' must have the same device mesh as the input tensors 'v' and 'b'. "
                f"Got device meshes {mask.device_mesh} and {device_mesh_input}."
            )
        if device_mesh_input != g.device_mesh:
            raise ValueError(
                f"Input tensor 'g' must have the same device mesh as the input tensors 'v' and 'b'. "
                f"Got device meshes {g.device_mesh} and {device_mesh_input}."
            )

        # Check if inputs have compatible placements
        placements_input = v.placements
        if placements_input != b.placements:
            raise ValueError(
                f"Input tensors 'v' and 'b' must have identical placements. "
                f"Got placements {placements_input} and {b.placements}."
            )
        if placements_input != mask.placements:
            raise ValueError(
                f"Input tensor 'mask' must have the same placements as the input tensors 'v' and 'b'. "
                f"Got placements {mask.placements} and {placements_input}."
            )
        if placements_input != g.placements:
            raise ValueError(
                f"Input tensor 'g' must have the same placements as the input tensors 'v' and 'b'. "
                f"Got placements {g.placements} and {placements_input}."
            )
        if placements_input != (Shard(0), Shard(1), Shard(2)):
            # For debugging, we requires the placements to be (Shard(0), Shard(1), Shard(2))
            # TODO: remove this to only use the previous check
            raise ValueError(
                f"Input tensors 'v', 'b' and 'mask's placements are not (Shard(0), Shard(1), Shard(2)). "
                f"Got placements {placements_input}."
            )

        # Check tensor dimensions
        if v.ndim != 4:
            raise ValueError(f"Input tensor 'v' must have 4 dimensions. Got {v.ndim} dimensions.")
        if b.ndim != 4:
            raise ValueError(f"Input tensor 'b' must have 4 dimensions. Got {b.ndim} dimensions.")
        if mask.ndim != 3:
            raise ValueError(f"Input tensor 'mask' must have 3 dimensions. Got {mask.ndim} dimensions.")
        if g.ndim != 4:
            raise ValueError(f"Input tensor 'g' must have 4 dimensions. Got {g.ndim} dimensions.")

        v_local = v.to_local()
        b_local = b.to_local()
        mask_local = mask.to_local()
        g_local = g.to_local()

        # Check shape compatibility
        B, S, N, h_c_h = v_local.shape
        if h_c_h % n_heads != 0:
            raise ValueError(
                f"Input tensor 'v' must have a last dimension divisible by n_heads {n_heads}. Got {h_c_h}."
            )
        c_h = h_c_h // n_heads
        if b_local.shape != (B, N, N, n_heads):
            raise ValueError(
                f"Input tensor 'b' must have shape ({B}, {N}, {N}, {n_heads}). Got shape {b.to_local().shape}."
            )
        if mask_local.shape != (B, N, N):
            raise ValueError(f"Input tensor 'mask' must have shape ({B}, {N}, {N}). Got shape {mask.to_local().shape}.")
        if g_local.shape != (B, S, N, h_c_h):
            raise ValueError(
                f"Input tensor 'g' must have shape ({B}, {S}, {N}, {h_c_h}). Got shape {g.to_local().shape}."
            )

        # Perform consistency check between the comm and the device_mesh_input
        # Check if Shard(1) and Shard(2) exist in placements_input (for S and N dimensions)
        i_tensor_dim_to_i_grid_axis = [-1] * v.ndim
        for i_grid_axis, placement in enumerate(placements_input):
            if isinstance(placement, Shard):
                i_tensor_dim_to_i_grid_axis[placement.dim] = i_grid_axis
        if i_tensor_dim_to_i_grid_axis[1] == -1 or i_tensor_dim_to_i_grid_axis[2] == -1:
            raise ValueError(
                f"Input tensors 'v', 'b' and 'mask's dimensions 1 and 2 must be sharded."
                f"Got placements {placements_input}."
            )

        # Check if ring_comm.group_col match the device_mesh_input's group
        if comm.group_col != device_mesh_input.get_group(i_tensor_dim_to_i_grid_axis[1]):
            raise ValueError(
                "Input comm's group_col process group is not the same as the group sharding the input tensors' axis 1"
            )

        # Check if the rank coordinates are consistent
        coord_device_mesh_input = device_mesh_input.get_coordinate()
        if coord_device_mesh_input is None:
            raise ValueError(f"comm.coord_2d {comm.coord_2d} is not on device_mesh_input {device_mesh_input}.")
        if comm.coord_2d != (
            coord_device_mesh_input[i_tensor_dim_to_i_grid_axis[1]],
            coord_device_mesh_input[i_tensor_dim_to_i_grid_axis[2]],
        ):
            raise ValueError(
                f"Input comm's coord_2d {comm.coord_2d} does not match the "
                f"device mesh's rank coordinates {coord_device_mesh_input} for the sharded dimensions "
                f"{i_tensor_dim_to_i_grid_axis[1]} and {i_tensor_dim_to_i_grid_axis[2]}."
            )

        requires_grad = v.requires_grad or b.requires_grad
        ctx.mark_non_differentiable(mask)

        # # Save device mesh and placements for backward pass
        ctx.device_mesh_input = device_mesh_input
        ctx.placements_input = placements_input
        ctx.input_shape = v.shape
        ctx.input_stride = v.stride()
        ctx.g_shape = g.shape
        ctx.g_stride = g.stride()
        ctx.b_shape = b.shape
        ctx.b_stride = b.stride()
        ctx.comm = comm
        ctx.n_heads = n_heads
        ctx.c_h = c_h

        # reshape v_local from (B, S, N, n_heads * c_h) to (B, n_heads, S, N, c_h)
        # TODO: reshape v_local to (B, n_heads, c_h, S, N) to be more GEMM friendly
        v_local = (
            v_local.unflatten(-1, (n_heads, c_h)).permute(0, 3, 1, 2, 4).clone(memory_format=torch.contiguous_format)
        )

        if requires_grad:
            v_local_copy = v_local.detach().clone()
        else:
            v_local_copy = None

        v_local_recv = comm.comm_row_init.enqueue_to_dispatch(v_local)

        # convert mask_local to bias and reshape it to (B, N0, N1, 1)
        mask_bias_local = 1 - mask_local
        mask_bias_local *= -inf
        mask_bias_local = mask_bias_local.unsqueeze(-1)

        # apply mask_bias to b and reshape b_local from (B, N0, N1, n_heads) to (B, n_heads, N1, N0)
        b_local = b_local + mask_bias_local
        bT_local = b_local.permute(0, 3, 2, 1).contiguous()

        bT_local_recv = comm.comm_2d_trans.enqueue_to_dispatch(bT_local)

        g_local = g_local.sigmoid()

        comm.comm_2d_trans.wait_until_finished()

        bT_local = comm.comm_col_init.enqueue_to_dispatch(bT_local_recv, bT_local)

        # cumulative amax until the current block
        amax = None
        # cumulative lse_m
        lse_m = None

        # bT_local is ready
        comm.comm_col_init.wait_until_finished()

        # v_local_recv is ready
        comm.comm_row_init.wait_until_finished()

        bT_local_buffer = [bT_local, bT_local_recv]
        v_local_buffer = [v_local_recv, v_local]
        i_ready = 0
        i_recv = i_ready ^ 1

        # receive other tensor blocks
        num_steps = comm.group_layout.shape[0]
        o_local = None
        for step in range(num_steps):
            there_is_another_step = (step + 1) < num_steps
            if there_is_another_step:
                v_local_buffer[i_recv] = comm.comm_row.enqueue_to_dispatch(
                    v_local_buffer[i_ready], v_local_buffer[i_recv]
                )
                bT_local_buffer[i_recv] = comm.comm_col.enqueue_to_dispatch(
                    bT_local_buffer[i_ready], bT_local_buffer[i_recv]
                )

            # (B, n_heads, 1, N)
            amax_block = bT_local_buffer[i_ready].amax(dim=-2, keepdim=True)
            lse_m_block = torch.logsumexp(bT_local_buffer[i_ready] - amax_block, dim=-2, keepdim=True)

            p = bT_local_buffer[i_ready].softmax(dim=-2)

            o_block = torch.einsum("bhsjd,bhji->bhsid", v_local_buffer[i_ready], p)

            # reshape o_block from (B, n_heads, S, N, c_h) to (B, n_heads, N, S * c_h)
            # reshape lse_m_block from (B, n_heads, 1, N) to (B, n_heads, N, 1)
            # reshape amax_block from (B, n_heads, 1, N) to (B, n_heads, N, 1)
            o_local, lse_m, amax = tiled_softmax_attention_update(
                o_block.transpose(-3, -2).flatten(start_dim=-2),
                lse_m_block.transpose(-2, -1),
                amax_block.transpose(-2, -1),
                o_local,
                lse_m,
                amax,
            )

            if there_is_another_step:
                comm.comm_row.wait_until_finished()
                comm.comm_col.wait_until_finished()
                i_ready = i_ready ^ 1
                i_recv = i_recv ^ 1

        # reshape o_local from (B, n_heads, N, S *c_h) to (B, S, N, n_heads * c_h)
        o_local = o_local.unflatten(-1, (S, c_h)).permute(0, 3, 2, 1, 4).flatten(start_dim=-2)

        o_local_copy = o_local.detach().clone(memory_format=torch.contiguous_format)

        if requires_grad:
            # transpose lse_m and amax
            amax_recv = comm.comm_2d_trans_amax.enqueue_to_dispatch(amax)
            # normalize b_local to be post-softmax matrix
            # b_local is of shape: (B, n_heads, N0, N1)
            b_local = b_local.permute(0, 3, 1, 2).contiguous()
            lse_m_recv = comm.comm_2d_trans_lse_m.enqueue_to_dispatch(lse_m)
            # subtract amax from b_local first
            # as amax tends to store extreme values from b_local masked
            # lse_m and amax of shape: (B, n_heads, N0, 1)
            # transpose lse across the grid to match b_local's placements
            # b_local = torch.exp(b_local - amax - lse_m)
            comm.comm_2d_trans_amax.wait_until_finished()
            b_local -= amax_recv
            comm.comm_2d_trans_lse_m.wait_until_finished()
            b_local -= lse_m_recv
            b_local.exp_()
            ctx.save_for_backward(v_local_copy, b_local, g_local, o_local_copy)

        o_local *= g_local

        o = DTensor.from_local(
            o_local,
            device_mesh=device_mesh_input,
            placements=placements_input,
            shape=ctx.input_shape,
            stride=ctx.input_stride,
        )

        return o

    @staticmethod
    @torch.amp.custom_bwd(device_type="cuda")
    def backward(ctx, do: DTensor) -> tuple[DTensor, DTensor, None, None, None, None, None]:
        if not isinstance(do, DTensor):
            raise TypeError(f"Input 'do' must be of type DTensor. Got type {type(do)}.")

        if do.device_mesh != ctx.device_mesh_input:
            raise ValueError(
                f"Input 'do' must have the same device mesh as the input tensors. "
                f"Got device meshes {do.device_mesh} and {ctx.device_mesh_input}."
            )

        if do.placements != ctx.placements_input:
            raise ValueError(
                f"Input 'do' must have the same placements as the input tensors. "
                f"Got placements {do.placements} and {ctx.placements_input}."
            )

        v_local, p_local, g_local, o_local = ctx.saved_tensors
        comm = ctx.comm
        num_steps = comm.group_layout.shape[0]

        S = do.to_local().shape[1]

        p_local_recv = comm.comm_col_init.enqueue_to_dispatch(p_local)

        # do.to_local() is of shape (B, S, N, n_heads * c_h)
        # g_local is of shape (B, S, N, n_heads * c_h)
        # cast do to the same dtype as g_local to avoid type promotion to FP32 if do is FP32
        # which can cause NCCL hang due to size mismatch in P2P communication
        do_local = do.to_local().to(dtype=g_local.dtype) * g_local
        dsigmoid = 1 - g_local

        # input o_local is of shape (B, S, N, n_heads * c_h)
        dsigmoid *= do_local
        dsigmoid *= o_local

        dg = DTensor.from_local(
            dsigmoid,
            device_mesh=ctx.device_mesh_input,
            placements=ctx.placements_input,
            shape=ctx.g_shape,
            stride=ctx.g_stride,
        )

        # reshape do from (B, S, N, n_heads * c_h) to (B, n_heads, S * c_h, N)
        do_local = (
            do_local.unflatten(-1, (ctx.n_heads, ctx.c_h))  # (B, S, N, n_heads, c_h)
            .permute(0, 3, 1, 4, 2)  # (B, n_heads, S, c_h, N)
            .flatten(start_dim=-3, end_dim=-2)  # (B, n_heads, S * c_h, N)
            .contiguous()
        )

        do_local_recv = comm.comm_row_init.enqueue_to_dispatch(do_local)

        # reshape o_local from (B, S, N, n_heads * c_h) to (B, n_heads, S * c_h, N)
        o_local = (
            o_local.unflatten(-1, (ctx.n_heads, ctx.c_h))  # (B, S, N, n_heads, c_h)
            .permute(0, 3, 1, 4, 2)  # (B, n_heads, S, c_h, N)
            .flatten(start_dim=-3, end_dim=-2)  # (B, n_heads, S * c_h, N)
            .contiguous()
        )

        d = torch.einsum("bhti,bhti->bhi", do_local, o_local).contiguous()

        d_recv = comm.comm_d_init.enqueue_to_dispatch(d)

        # reshape v_local from (B, n_heads, S, N, c_h) to (B, n_heads, S * c_h, N)
        v_local = v_local.transpose(-2, -1).flatten(start_dim=-3, end_dim=-2).contiguous()

        comm.comm_col_init.wait_until_finished()  # p_local_recv is ready
        comm.comm_row_init.wait_until_finished()  # do_local_recv is ready
        comm.comm_d_init.wait_until_finished()  # d_recv is ready

        i_ready = 0
        i_recv = i_ready ^ 1

        p_local_buffer = [p_local_recv, p_local]
        do_local_buffer = [do_local_recv, do_local]
        d_buffer = [d_recv, d]
        db_local_buffer = [torch.zeros_like(p_local), torch.zeros_like(p_local)]

        dv_local = torch.zeros_like(o_local)

        for step in range(num_steps):
            there_is_another_step = step < num_steps - 1
            if there_is_another_step:
                p_local_buffer[i_recv] = comm.comm_col.enqueue_to_dispatch(
                    p_local_buffer[i_ready], p_local_buffer[i_recv]
                )
                do_local_buffer[i_recv] = comm.comm_row.enqueue_to_dispatch(
                    do_local_buffer[i_ready], do_local_buffer[i_recv]
                )
                d_buffer[i_recv] = comm.comm_d.enqueue_to_dispatch(d_buffer[i_ready], d_buffer[i_recv])

            dv_block = torch.einsum("bhti,bhij->bhtj", do_local_buffer[i_ready], p_local_buffer[i_ready])
            dv_local += dv_block

            # dp
            db_local_block = torch.einsum("bhti,bhtj->bhij", do_local_buffer[i_ready], v_local).contiguous()
            # dp - d
            db_local_block -= d_buffer[i_ready].unsqueeze(-1)
            # p * (dp - d)
            db_local_block *= p_local_buffer[i_ready]
            # virtual all-reduce db_local
            if step > 0:
                comm.comm_db.wait_until_finished()
            db_local_buffer[i_ready] += db_local_block
            # db send/recv will carry through num_steps
            # explicitly cast to the recv buffer's dtype to prevent NCCL hang due to potential dtype mismatch
            # if db_local_buffer[i_ready] was promoted to FP32 during accumulation
            db_local_buffer[i_recv] = comm.comm_db.enqueue_to_dispatch(
                db_local_buffer[i_ready].to(dtype=db_local_buffer[i_recv].dtype), db_local_buffer[i_recv]
            )

            if there_is_another_step:
                comm.comm_col.wait_until_finished()
                comm.comm_row.wait_until_finished()
                comm.comm_d.wait_until_finished()
                i_ready = i_ready ^ 1
                i_recv = i_recv ^ 1

        # reshape dv_local from (B, n_heads, S * c_h, N) to (B, S, N, n_heads * c_h)
        dv_local = dv_local.unflatten(-2, (S, ctx.c_h)).permute(0, 2, 4, 1, 3).flatten(start_dim=-2)

        dv = DTensor.from_local(
            dv_local,
            device_mesh=ctx.device_mesh_input,
            placements=ctx.placements_input,
            shape=ctx.input_shape,
            stride=ctx.input_stride,
        )

        # reshape db_local from (B, n_heads, N0, N1) to (B, N0, N1, n_heads)
        # the last comm of db is ready
        comm.comm_db.wait_until_finished()
        i_ready = i_ready ^ 1
        i_recv = i_recv ^ 1

        # revert the comm_col_init to get the db data ownership as the input b
        db_local_buffer[i_recv] = comm.comm_db_final.enqueue_to_dispatch(
            db_local_buffer[i_ready], db_local_buffer[i_recv]
        )

        comm.comm_db_final.wait_until_finished()
        i_ready = i_ready ^ 1
        i_recv = i_recv ^ 1

        # reshape db_local from (B, n_heads, N0, N1) to (B, N0, N1, n_heads)
        db_local = db_local_buffer[i_ready].permute(0, 2, 3, 1)

        db = DTensor.from_local(
            db_local,
            device_mesh=ctx.device_mesh_input,
            placements=ctx.placements_input,
            shape=ctx.b_shape,
            stride=ctx.b_stride,
        )

        return dv, db, None, dg, None, None, None


class PairWeightedAveraging(torch.nn.Module):
    """Pair weighted averaging layer."""

    def __init__(
        self,
        layer: SerialPairWeightedAveraging,
        device_mesh: DeviceMesh,
        comm: Ring2DCommPairAveraging,
    ) -> None:
        super().__init__()
        self.comm = comm
        self.c_m = layer.c_m
        self.c_z = layer.c_z
        self.c_h = layer.c_h
        self.num_heads = layer.num_heads
        self.inf = layer.inf

        self.device_mesh = device_mesh
        self.comm = comm

        self.norm_m = LayerNormParamsReplicated(layer.norm_m, self.device_mesh)
        self.norm_z = LayerNormParamsReplicated(layer.norm_z, self.device_mesh)

        self.proj_m = LinearParamsReplicated(layer.proj_m, self.device_mesh)
        self.proj_g = LinearParamsReplicated(layer.proj_g, self.device_mesh)
        self.proj_z = LinearParamsReplicated(layer.proj_z, self.device_mesh)
        self.proj_o = LinearParamsReplicated(layer.proj_o, self.device_mesh)

    def forward(self, m: DTensor, z: DTensor, mask: DTensor) -> DTensor:
        # Compute layer norms
        m = self.norm_m(m)
        z = self.norm_z(z)

        g = self.proj_g(m)

        # TODO: fuse the m -> {v, g} projection in one kernel
        v = self.proj_m(m)

        b = self.proj_z(z)

        o = _PairWeightedAveragingImpl.apply(v, b, mask, g, self.comm, self.num_heads, self.inf)

        o = self.proj_o(o)
        return o
