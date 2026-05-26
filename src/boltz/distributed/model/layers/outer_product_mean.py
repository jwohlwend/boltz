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
import torch.distributed as dist
from torch import nn
from torch.distributed.tensor import DTensor, Shard
from torch.distributed.tensor.device_mesh import DeviceMesh

from boltz.distributed.comm import Ring2DComm
from boltz.distributed.model.layers.layernorm import LayerNormParamsReplicated
from boltz.distributed.model.layers.linear import LinearParamsReplicated
from boltz.distributed.utils import LayoutRightMap
from boltz.model.layers.outer_product_mean import OuterProductMean as SerialOuterProductMean


class _OuterProductMeanImpl(torch.autograd.Function):
    """Distributed implementation of outer product mean using ring communication.

    This autograd function implements a memory-efficient distributed outer product mean
    operation across a 2D process grid. The computation is parallelized using ring
    communication patterns to reduce memory usage and communication overhead.

    The outer product mean computes:
    z[i,j] = mean_s(a[s,i] ⊗ b[s,j])

    where ⊗ denotes outer product, and the mean is taken over the sequence dimension s.

    Key features:
    - Distributed across a 2D grid with sharding on sequence (dim 1) and token (dim 2) dimensions
    - Uses ring communication to rotate data chunks during computation
    - Memory-efficient implementation that avoids materializing full tensors
    - Supports gradient computation through custom backward pass

    Notes
    -----
    Input tensors must be DTensors with:
    - Shape: (B, N_seq, N_token, c_hidden) for tensors a and b
    - Shape: (B, N_seq, N_token, 1) for mask tensor
    - Sharding on dimensions 1 and 2 (Shard(1) and Shard(2) placements)
    - Identical device mesh and placements across all inputs

    The algorithm uses a ring-based communication pattern where:
    - Tensor a is transposed and rotated by row
    - Tensor b is rotated by column
    - Each process computes partial outer products and accumulates results
    """

    @staticmethod
    @torch.amp.custom_fwd(device_type="cuda")
    def forward(ctx, a: DTensor, b: DTensor, mask: DTensor, ring_comm: Ring2DComm):
        """Forward pass of distributed outer product mean computation.

        Computes the outer product mean of input tensors a and b using distributed
        ring communication to minimize memory usage and communication overhead.

         a/b sharding is done along the s and N dimensions (dim 1 and 2).
        For an i, j rank, the output is the sum of outer products over a[..., i] and b[..., j].
        Example initial layout for a 2x2 grid::
         [[a00, a01], and [[b00, b01],
          [a10, a11]] and  [b10, b11]]
        a is transposed and rotated by row, b is NOT transposed, and is rotated by column.
        After transpose, the layout is:
         [[a00, a10],
          [a01, a11]]
        An initial offset is added for each. For a, row i is rotated by i elements, and likewise for the columns of b.
        After offset,
         [[a00, a10], and [[b00, b11],
          [a11, a10]] and  [b10, b10]]
         Note that (for example) grid element (1, 0) has a[11] and b[10], corresponding to a[.., i] and b[..., ,j],
         with matching secondary index of 1.
        After 1 rotation
         [[a10, a00], and [[b10, b00],
          [a01, a11]] and  [b00, b11]],
        the same (1, 0) index has a01 and b00, with the same i,j match.


        Parameters
        ----------
        ctx : torch.autograd.function.BackwardCFrame
            Context object for saving information needed in backward pass.
        a : DTensor
            First input tensor with shape (B, N_seq, N_token, c_hidden).
            Must be sharded on dimensions 1 and 2.
        b : DTensor
            Second input tensor with shape (B, N_seq, N_token, c_hidden).
            Must have identical shape, device mesh, and placements as tensor a.
        mask : DTensor
            Mask tensor with shape (B, N_seq, N_token).
            Must have same device mesh and placements as input tensors.
        ring_comm : Ring2DComm
            Ring communication object configured for the distributed computation.

        Returns
        -------
        DTensor
            Output tensor with shape (B, N_token, N_token, c_hidden*c_hidden).
            Contains the distributed outer product mean result.

        Raises
        ------
        TypeError
            If inputs are not DTensor type.
        ValueError
            If tensor shapes, device meshes, or placements are incompatible,
            or if ring communication setup is inconsistent.
        """
        # Check if inputs a and b are of type DTensor
        if not isinstance(a, DTensor) or not isinstance(b, DTensor) or not isinstance(mask, DTensor):
            raise TypeError(
                f"Inputs 'a', 'b', and 'mask' must be of type DTensor. Got types {type(a)}, {type(b)}, and {type(mask)}."
            )

        # Check if inputs a and b have identical device mesh
        device_mesh_input = a.device_mesh
        if device_mesh_input != b.device_mesh:
            raise ValueError(
                f"Input tensors 'a' and 'b' must have identical device mesh. "
                f"Got device meshes {device_mesh_input} and {b.device_mesh}."
            )
        if device_mesh_input != mask.device_mesh:
            raise ValueError(
                f"Input tensor 'mask' must have the same device mesh as the input tensors 'a' and 'b'. "
                f"Got device meshes {mask.device_mesh} and {device_mesh_input}."
            )

        # Check if inputs a and b have identical placements
        placements_input = a.placements
        if placements_input != b.placements:
            raise ValueError(
                f"Input tensors 'a' and 'b' must have identical placements. "
                f"Got placements {placements_input} and {b.placements}."
            )
        if placements_input != mask.placements:
            # TODO: in the future, if a and b are sharded along the hidden dimension, we need to
            # skip the check on the corresponding grid axes
            raise ValueError(
                f"Input tensor 'mask' must have the same placements as the input tensors 'a' and 'b'. "
                f"Got placements {mask.placements} and {placements_input}."
            )
        if placements_input != (Shard(0), Shard(1), Shard(2)):
            # For debugging, we requires the placements to be (Shard(0), Shard(1), Shard(2))
            # TODO: remove this to only use the previous check
            raise ValueError(
                f"Input tensors 'a' and 'b''s placements are not (Shard(0), Shard(1), Shard(2)). "
                f"Got placements {placements_input}."
            )

        # Check if inputs a and b have the same shape
        if a.shape != b.shape:
            raise ValueError(f"Input tensors 'a' and 'b' must have the same shape. Got shapes {a.shape} and {b.shape}.")

        if mask.shape[:3] != a.shape[:3]:
            raise ValueError(
                f"Input tensor 'mask' doesn't have the same size in the 3 leading dimensions "
                f"as the input tensors 'a' and 'b': Got shape mask's shape: {mask.shape} vs. a.shape: {a.shape}"
            )

        # to stay off potential DTensor issue, let's not use negative dim axis when we can avoid it
        # but this requires we assume the semantics of the axes in order to check the placements:
        # a.shape [B, N_seq, N_token, c_hidden]
        if a.ndim != 4:
            raise ValueError(
                f"Input tensors 'a' and 'b' must have 4 dimensions. "
                f"Got {a.ndim} dimensions for tensor 'a' and {b.ndim} dimensions for tensor 'b'."
            )

        # Perform consistency check between the ring_comm and the device_mesh_input
        # TODO: we could also check the device_mesh_input.mesh tensor, but it would be too expensive as the check
        # will likely go through all elements on the mesh
        # NOTE: leading batch dimensions may or may not be sharded, as this algorithm operates orthogonally to them.
        # 1. Check if Shard(1) and Shard(2) exist in placements_input
        i_tensor_dim_to_i_grid_axis = [-1] * a.ndim
        for i_grid_axis, placement in enumerate(placements_input):
            if isinstance(placement, Shard):
                i_tensor_dim_to_i_grid_axis[placement.dim] = i_grid_axis
        if i_tensor_dim_to_i_grid_axis[1] == -1 or i_tensor_dim_to_i_grid_axis[2] == -1:
            raise ValueError(
                f"Input tensors 'a', 'b' and 'mask's dimensions 1 and 2 must be sharded. Got placements {placements_input}."
            )
        # 2. Check if ring_comm.group_col match the device_mesh_input's group
        # NOTE: ring_comm.group_col is the group sharding the input tensors' axis 1
        if ring_comm.group_col != device_mesh_input.get_group(i_tensor_dim_to_i_grid_axis[1]):
            raise ValueError(
                "Input ring_comm's group_col process group is not the same as the group sharding the input tensors' axis 1"
            )
        # 3. Check if the rank coordinates are consistent
        coord_device_mesh_input = device_mesh_input.get_coordinate()
        if coord_device_mesh_input is None:
            raise ValueError(
                f"ring_comm.coord_2d {ring_comm.coord_2d} is not on device_mesh_input {device_mesh_input}."
            )
        if ring_comm.coord_2d != (
            coord_device_mesh_input[i_tensor_dim_to_i_grid_axis[1]],
            coord_device_mesh_input[i_tensor_dim_to_i_grid_axis[2]],
        ):
            raise ValueError(
                f"Input ring_comm's coord_2d {ring_comm.coord_2d} does not match the "
                f"device mesh's rank coordinates {coord_device_mesh_input} for the sharded dimensions "
                f"{i_tensor_dim_to_i_grid_axis[1]} and {i_tensor_dim_to_i_grid_axis[2]}."
            )

        ctx.mark_non_differentiable(mask)
        mask_local = mask.to_local().unsqueeze(-1)
        # DTensor.to_local() returns a view to the shard so we need to clone it to avoid modifying the original DTensor.
        a_local = a.to_local() * mask_local
        b_local = b.to_local() * mask_local
        a_local_copy = a_local.detach().clone()
        b_local_copy = b_local.detach().clone()

        # Sum mask count along columns to get divisor for mean.

        B, _, N, c_hidden = a_local.shape

        # send off A transpose + row init, b column init
        a_recv = ring_comm.comm_transpose_row_init.enqueue_to_dispatch(a_local)
        b_recv = ring_comm.comm_col_init.enqueue_to_dispatch(b_local)

        z_local = torch.zeros((B, N, N, c_hidden, c_hidden), dtype=a_local.dtype, device=a_local.device)

        ring_comm.comm_col_init.wait_until_finished()
        ring_comm.comm_transpose_row_init.wait_until_finished()

        num_mask_local = mask_local[:, :, 0].sum(1)[:, None, None]
        num_mask_work = dist.all_reduce(num_mask_local, group=ring_comm.group_col, async_op=True)
        a_buffer = [a_recv, a_local]
        b_buffer = [b_recv, b_local]
        i_ready = 0
        i_recv = i_ready ^ 1
        for k_step in range(ring_comm.group_layout.shape[1]):
            a_ready = a_buffer[i_ready]
            b_ready = b_buffer[i_ready]
            if k_step < ring_comm.group_layout.shape[1] - 1:
                a_buffer[i_recv] = ring_comm.comm_row.enqueue_to_dispatch(a_ready, a_buffer[i_recv])
                b_buffer[i_recv] = ring_comm.comm_col.enqueue_to_dispatch(b_ready, b_buffer[i_recv])
            z_local = z_local + torch.einsum("bsic,bsjd->bijcd", a_ready, b_ready)
            if k_step < ring_comm.group_layout.shape[1] - 1:
                ring_comm.comm_row.wait_until_finished()
                ring_comm.comm_col.wait_until_finished()
                i_ready = i_ready ^ 1
                i_recv = i_recv ^ 1

        num_mask_work.wait()
        num_mask_local_clamped = num_mask_local.clamp(min=1)
        z_local = z_local.flatten(start_dim=-2) / num_mask_local_clamped

        # Compute output shape and stride
        shape_output = (a.shape[0], a.shape[2], a.shape[2], z_local.shape[-1])  # (B, N, N, c_hidden * c_hidden)

        # Use LayoutRightMap for the output shape
        layout_right = LayoutRightMap(shape_output)
        strides_output = layout_right.strides

        if a.requires_grad or b.requires_grad:
            ctx.save_for_backward(a_local_copy, b_local_copy, mask_local.detach().clone(), num_mask_local_clamped)
            ctx.ring_comm = ring_comm
            ctx.placements_input = placements_input
            ctx.device_mesh_input = device_mesh_input
            ctx.input_shape_a = a.shape
            ctx.input_stride_a = a.stride()
            ctx.input_shape_b = b.shape
            ctx.input_stride_b = b.stride()

        z = DTensor.from_local(
            z_local,
            device_mesh=device_mesh_input,
            placements=placements_input,
            shape=shape_output,
            stride=strides_output,
        )
        return z

    @staticmethod
    @torch.amp.custom_bwd(device_type="cuda")
    def backward(ctx, grad_z: DTensor) -> tuple[DTensor, DTensor, None, None]:
        """Backward pass of distributed outer product mean computation.

        Computes gradients with respect to input tensors a and b using the same
        ring communication pattern as the forward pass but with transposed operations.

        The gradient computation follows:
        - grad_a[s,i] = sum_j(grad_z[i,j] ⊗ b[s,j])
        - grad_b[s,j] = sum_i(grad_z[i,j] ⊗ a[s,i])

        The sharding of a and b is described above in the forward pass.
        z is similarly sharded, but the two sharding dimensions are both (N, N).
        For an i, j rank:
        Gradient of A is the sum of outer products over grad_z[j, ...] and b[j, ...].
        Gradient of B is the sum of outer products over grad_z[..., j] and a[i, ...].
        For the gradient of a, the approach is identical to that of the forward pass, this time with transposition of
        grad_z, rotation of b by row, and grad_z by column. View the above schematic for details.
        For the gradient of b, the approach is similar, but with no transposition of grad_z.
        Starting with:
         [[g00, g01], and [[a00, a01],
          [g10, g11]] and  [a10, a11]]
        After initial rotation:
         [[g00, g11], and [[a00, a01],
          [g10, g01]] and  [a11, a10]]
        where index i,j =(0, 1) has g11 and a01 corresponding to grad_z[..., j] and a[i, ...], with matching secondary
        index of 1.
        The next rotation of grad_z up and b left yields:
         [[g10, g01], and [[a01, a00],
          [g00, g11]] and  [a11, a10]]
        where index i,j = (0, 1) has g01 and a00, corresponding to grad_z[..., j] and a[i, ...], with matching
        secondary index of 0.

        Parameters
        ----------
        ctx : torch.autograd.function.BackwardCFrame
            Context object containing saved tensors and communication setup.
        grad_z : DTensor
            Gradient tensor from upstream with shape (B, N_token, N_token, c_hidden*c_hidden).
            Must have same device mesh and placements as forward pass inputs.

        Returns
        -------
        tuple[DTensor, DTensor, None, None]
            Gradients with respect to inputs:
            - grad_a: Gradient for tensor a with shape (B, N_seq, N_token, c_hidden)
            - grad_b: Gradient for tensor b with shape (B, N_seq, N_token, c_hidden)
            - None: No gradient for mask (marked non-differentiable)
            - None: No gradient for ring_comm

        Raises
        ------
        TypeError
            If grad_z is not a DTensor.
        ValueError
            If grad_z has incompatible shape, device mesh, or placements.
        """

        if not isinstance(grad_z, DTensor):
            raise TypeError(f"Input 'grad_z' must be of type DTensor. Got type {type(grad_z)}.")

        if grad_z.ndim != 4:
            raise ValueError(f"Input 'grad_z' must have 4 dimensions but got {grad_z.ndim} dimensions")

        if grad_z.device_mesh != ctx.device_mesh_input:
            raise ValueError(
                f"Input 'grad_z' must have the same device mesh as the input tensors 'a' and 'b'. "
                f"Got device meshes {grad_z.device_mesh} and {ctx.device_mesh_input}."
            )

        if grad_z.placements != ctx.placements_input:
            raise ValueError(
                f"Input 'grad_z' must have the same placements as the input tensors 'a' and 'b'. "
                f"Got placements {grad_z.placements} and {ctx.placements_input}."
            )

        a_local, b_local, mask_local, num_mask_local_clamped = ctx.saved_tensors

        # reshape grad_z's last axis to perform GEMV
        c_hidden = a_local.shape[-1]
        # apply the mask and clone
        grad_z_local = grad_z.to_local().unflatten(-1, (c_hidden, c_hidden)) / num_mask_local_clamped.unsqueeze(-1)

        if grad_z_local.shape[:3] != (
            a_local.shape[0],
            a_local.shape[2],
            a_local.shape[2],
        ):
            raise ValueError(
                f"grad_z shard shape {grad_z_local.shape} does not match expected shape "
                f"({a_local.shape[0]}, {a_local.shape[2]}, {a_local.shape[2]}) for outer product mean."
            )

        # Save for reset for B grad compute.
        grad_z_save = grad_z_local.clone()

        ring_comm: Ring2DComm = ctx.ring_comm

        # Initialize gradient buffers
        grad_a_local = torch.zeros_like(a_local)

        # For grad_a computation: Transpose z grad, rotate z by column and b by row.
        grad_z_recv = ring_comm.comm_transpose_col_init.enqueue_to_dispatch(grad_z_local)
        b_recv = ring_comm.comm_row_init.enqueue_to_dispatch(b_local)
        ring_comm.comm_transpose_col_init.wait_until_finished()
        ring_comm.comm_row_init.wait_until_finished()

        b_buffer = [b_recv, b_local]
        grad_z_buffer = [grad_z_recv, grad_z_local]
        i_ready = 0
        i_recv = i_ready ^ 1

        # Compute grad_a
        for k_step in range(ring_comm.group_layout.shape[1]):
            b_ready = b_buffer[i_ready]
            grad_z_ready = grad_z_buffer[i_ready]

            if k_step < ring_comm.group_layout.shape[1] - 1:
                b_buffer[i_recv] = ring_comm.comm_row.enqueue_to_dispatch(b_ready, b_buffer[i_recv])
                grad_z_buffer[i_recv] = ring_comm.comm_col.enqueue_to_dispatch(grad_z_ready, grad_z_buffer[i_recv])

            grad_a_local = grad_a_local + torch.einsum("bijcd,bsjd->bsic", grad_z_ready, b_ready)

            if k_step < ring_comm.group_layout.shape[1] - 1:
                ring_comm.comm_row.wait_until_finished()
                ring_comm.comm_col.wait_until_finished()
                i_ready = i_ready ^ 1
                i_recv = i_recv ^ 1

        grad_a_local = grad_a_local * mask_local
        grad_a = DTensor.from_local(
            grad_a_local,
            device_mesh=ctx.device_mesh_input,
            placements=ctx.placements_input,
            shape=ctx.input_shape_a,
            stride=ctx.input_stride_a,
        )

        # For grad_b computation: no transpose, rotate z by column, a by row.
        # NOTE: we're using the stored original grad_output.
        grad_z_recv = ring_comm.comm_col_init.enqueue_to_dispatch(grad_z_save)
        a_recv = ring_comm.comm_row_init.enqueue_to_dispatch(a_local)
        ring_comm.comm_row_init.wait_until_finished()
        ring_comm.comm_col_init.wait_until_finished()

        a_buffer = [a_recv, a_local]
        grad_z_buffer = [grad_z_recv, grad_z_local]
        i_ready = 0
        i_recv = i_ready ^ 1

        # Reuse b for grad_b, since we're done with it.
        grad_b_local = b_local.view(b_local.shape)
        grad_b_local *= 0.0
        # Compute grad_b
        for k_step in range(ring_comm.group_layout.shape[1]):
            a_ready = a_buffer[i_ready]
            grad_z_ready = grad_z_buffer[i_ready]

            if k_step < ring_comm.group_layout.shape[1] - 1:
                a_buffer[i_recv] = ring_comm.comm_row.enqueue_to_dispatch(a_ready, a_buffer[i_recv])
                grad_z_buffer[i_recv] = ring_comm.comm_col.enqueue_to_dispatch(grad_z_ready, grad_z_buffer[i_recv])

            grad_b_local = grad_b_local + torch.einsum("bijcd,bsic->bsjd", grad_z_ready, a_ready)

            if k_step < ring_comm.group_layout.shape[1] - 1:
                ring_comm.comm_row.wait_until_finished()
                ring_comm.comm_col.wait_until_finished()
                i_ready = i_ready ^ 1
                i_recv = i_recv ^ 1

        grad_b_local = grad_b_local * mask_local
        grad_b = DTensor.from_local(
            grad_b_local,
            device_mesh=ctx.device_mesh_input,
            placements=ctx.placements_input,
            shape=ctx.input_shape_b,
            stride=ctx.input_stride_b,
        )

        return grad_a, grad_b, None, None


class OuterProductMean(nn.Module):
    """Distributed outer product mean layer for sequence-to-pair transformations.

    This layer implements a distributed version of the outer product mean operation,
    which transforms sequence representations into pairwise representations. It's
    commonly used in protein structure prediction and other tasks requiring
    sequence-to-pair information propagation.

    The layer performs the following operations:
    1. Layer normalization of input sequences
    2. Linear projections to create two representation streams (a and b)
    3. Distributed outer product mean computation using ring communication
    4. Final linear projection to output dimension

    The outer product mean operation computes:
    z[i,j] = mean_s(proj_a(norm(m))[s,i] ⊗ proj_b(norm(m))[s,j])

    where the mean is taken over the sequence dimension s, masked by the input mask.

    Parameters
    ----------
    layer : SerialOuterProductMean
        The serial outer product mean layer to convert to distributed version.
        Used to initialize projection weights and normalization parameters.
        The weights and biases of the input layer will be replicated across the device mesh.
    device_mesh : DeviceMesh
        The device mesh for distributed computation across multiple GPUs.
    comm : Ring2DComm
        Ring communication object for efficient distributed outer product computation.

    Attributes
    ----------
    c_hidden : int
        Hidden dimension size from the projection layers.
    c_in : int
        Input dimension size.
    norm : LayerNormParamsReplicated
        Distributed layer normalization.
    proj_a : LinearParamsReplicated
        First projection layer (input -> hidden).
    proj_b : LinearParamsReplicated
        Second projection layer (input -> hidden).
    proj_o : LinearParamsReplicated
        Output projection layer (hidden*hidden -> output).

    Notes
    -----
    This implementation requires input tensors to be DTensors with appropriate
    sharding patterns. The layer is designed for large-scale distributed training
    where memory efficiency is critical.
    """

    def __init__(self, layer: SerialOuterProductMean, device_mesh: DeviceMesh, comm: Ring2DComm) -> None:
        """Initialize the distributed outer product mean layer.

        Parameters
        ----------
        layer : SerialOuterProductMean
            The serial outer product mean layer containing weights to be distributed.
        device_mesh : DeviceMesh
            Device mesh defining the distributed computation topology.
        comm : Ring2DComm
            Ring communication handler for distributed outer product operations.

        """
        super().__init__()
        self.device_mesh = device_mesh
        self.ring_comm = comm
        self.c_hidden = layer.c_hidden
        self.c_in = layer.proj_a.in_features
        self.norm = LayerNormParamsReplicated(layer.norm, self.device_mesh)
        self.proj_a = LinearParamsReplicated(layer.proj_a, self.device_mesh)
        self.proj_b = LinearParamsReplicated(layer.proj_b, self.device_mesh)
        self.proj_o = LinearParamsReplicated(layer.proj_o, self.device_mesh)

    def forward(self, m: DTensor, mask: DTensor) -> DTensor:
        """Forward pass of the distributed outer product mean layer.

        Transforms sequence representations into pairwise representations using
        the distributed outer product mean operation with masking support.

        Parameters
        ----------
        m : DTensor
            Input sequence tensor with shape (B, S, N, c_in).
            - B: batch size
            - S: sequence length
            - N: number of tokens
            - c_in: input feature dimension
            It's expected that the sequence and token dimensions are both sharded
            the device mesh, which is further required to be the same as the one
            the weights of the layer are placed on
        mask : DTensor
            Mask tensor with shape (B, S, N) indicating valid positions.
            Values should be 1.0 for valid positions and 0.0 for masked positions.
            It's expected that the sequence and token dimensions are both sharded
            the device mesh, which is further required to be the same as the one
            the weights of the layer are placed on

        Returns
        -------
        DTensor
            Output pairwise tensor with shape (B, N, N, c_out).
            Contains pairwise representations between all token pairs.

        Notes
        -----
        The computation pipeline is:
        1. Apply layer normalization to input sequences
        2. Project to two hidden representations (a and b) and apply masking
        3. Compute distributed outer product mean: z[i,j] = mean_s(a[s,i] ⊗ b[s,j])
        4. Apply final projection to get output dimension
        """
        # No need to expand mask here because it's handled in _OuterProductMeanImpl

        # Compute projections
        m = self.norm(m)
        a = self.proj_a(m)
        b = self.proj_b(m)

        z = _OuterProductMeanImpl.apply(a, b, mask, self.ring_comm)
        z = self.proj_o(z)
        return z
