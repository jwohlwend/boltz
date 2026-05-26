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


########################################################################################
# On sharding strategy and DTensor placements of atom_to_token
########################################################################################
# For performance optimization, we shard atom_to_token to minimize communication overhead for any operations on atom_to_token. These include:
# 1. translation of token- to atom-level single representation through torch.bmm
# 2. translation of atom- to token-level single representation through torch.bmm
# 3. pair representation token-to-atom through einsum
#
# Instead of a distributed matrix multiplication, we pad atom_to_token with padding atoms and
# tokens for context parallelism, such that only the block diagonal locations are non-zero. For
# example, say we have 3 tokens and 4 atoms on cp = (2, 2), atom_to_token without parallelism can be:
#
# [[1, 0, 0],
#  [0, 1, 0],
#  [0, 1, 0],
#  [0, 0, 1]]
#
# In context parallelism, we first pad with virtual tokens such that n_tokens is divisible by cp_size,
# and then pad with virtual atoms such that n_atoms is divisible by cp_size. Let us visualize it on the device mesh:
#
#                                               1 0 | 0 0
# 1 0 | 0                 1 0 | 0 0             0 1 | 0 0
# 0 1 | 0    pad tokens   0 1 | 0 0  pad atoms  0 1 | 0 0
# ---------  --------->   ---------  -------->  ---------
# 0 1 | 0                 0 1 | 0 0             0 0 | 1 0
# 0 0 | 1                 0 0 | 1 0             0 0 | 0 0
#                                               0 0 | 0 0
#
# Note that all tokens inside of each shard now have every of their own atoms in the same shard. As
# such, any mapping between atom to token or vice versa can be done locally on the diagonal ranks.
#
# To enable this for the off-diagonal ranks, we broadcast the block diagonal matrix row-wise, and
# since all single representations are replicated by row, now all ranks can perform these mapping
# locally without communication.
#
# 1 0 | 0 0                     1 0 | 1 0
# 0 1 | 0 0                     0 1 | 0 1
# 0 1 | 0 0  row-wise broadcast 0 1 | 0 1
# ---------  -----------------> ---------
# 0 0 | 1 0                     1 0 | 1 0
# 0 0 | 0 0                     0 0 | 0 0
# 0 0 | 0 0                     0 0 | 0 0
#
# Mapping on pair representation is the exception where a transposition of atom_to_token on the
# device mesh is needed.
########################################################################################


import torch
from torch import Tensor
from torch.distributed.tensor import DTensor, Replicate, Shard

from boltz.distributed.comm import TransposeComm
from boltz.distributed.utils import update_exhaustive_strides


class SingleReprTokenToAtomFunction(torch.autograd.Function):
    """Autograd function for transforming token-level single representation to atom-level single representation."""

    @staticmethod
    def forward(
        ctx,
        token_single_repr: DTensor,
        atom_to_token: DTensor,
    ) -> DTensor:
        """
        Transform a token-level single representation to an atom-level single representation.

        Args:
            token_single_repr: The token-level single representation. Shape: (B, n_tokens, D) and placement: (Shard(0), Shard(1), Replicate())
            atom_to_token: The atom to token one-hot mapping except for padding atoms/tokens. Shape: (B, n_atoms, n_tokens) and placement: (Shard(0), Shard(1), Replicate())

        Returns:
            The atom-level single representation. Shape: (B, n_atoms, D)
        """
        single_repr_placements = (Shard(dim=0), Shard(dim=1), Replicate())  # same as atom_to_token placements
        if atom_to_token.placements != single_repr_placements:
            raise ValueError(
                f"Expect atom_to_token to have placements {single_repr_placements}, but got {atom_to_token.placements}"
            )
        if token_single_repr.placements != single_repr_placements:
            raise ValueError(
                f"Expect token_single_repr to have placements {single_repr_placements}, but got {token_single_repr.placements}"
            )

        # Perform local bmm and distribute
        atom_to_token_local = atom_to_token.to_local().to(
            dtype=token_single_repr.dtype
        )  # NOTE in case atom_to_token is int
        token_single_repr_local = token_single_repr.to_local()
        o = torch.einsum("bij,bj...->bi...", atom_to_token_local, token_single_repr_local)

        # Compute output shape and stride for the global DTensor
        # For einsum "bij,bj...->bi...", output shape is (B, n_atoms, D)
        # where the last axis could be omitted if the input token_single_repr is 2D
        # tensor
        shape_output = atom_to_token.shape[:2] + o.shape[2:]

        # Use LayoutRightMap for the output shape
        strides_output = update_exhaustive_strides(o.shape, o.stride(), shape_output)

        o = DTensor.from_local(
            o, atom_to_token.device_mesh, single_repr_placements, shape=shape_output, stride=strides_output
        )

        if token_single_repr.requires_grad:
            ctx.device_mesh = atom_to_token.device_mesh
            ctx.single_repr_placements = single_repr_placements
            ctx.token_single_repr_shape = token_single_repr.shape
            ctx.token_single_repr_stride = token_single_repr.stride()
            ctx.save_for_backward(atom_to_token_local)

        return o

    @staticmethod
    def backward(ctx, grad_output: DTensor) -> tuple[DTensor, None, None]:
        """
        Backward pass for single_repr_token_to_atom.
        """
        if grad_output.placements != (Shard(dim=0), Shard(dim=1), Replicate()):
            raise ValueError(
                f"Expect grad_output to have placements {(Shard(dim=0), Shard(dim=1), Replicate())}, but got {grad_output.placements}"
            )
        if grad_output.device_mesh != ctx.device_mesh:
            raise ValueError(
                f"Expect grad_output to have device mesh {ctx.device_mesh}, but got {grad_output.device_mesh}"
            )

        do = grad_output.to_local()
        (atom_to_token_local,) = ctx.saved_tensors

        d_token_single_repr = torch.einsum("bji,bj...->bi...", atom_to_token_local, do)
        d_token_single_repr = DTensor.from_local(
            d_token_single_repr,
            ctx.device_mesh,
            ctx.single_repr_placements,
            shape=ctx.token_single_repr_shape,
            stride=ctx.token_single_repr_stride,
        )

        return d_token_single_repr, None, None


class SingleReprAtomToTokenFunction(torch.autograd.Function):
    """Autograd function for transforming atom-level single representation to token-level single representation."""

    @staticmethod
    @torch.amp.custom_fwd(device_type="cuda")
    def forward(ctx, atom_single_repr: DTensor, atom_to_token: DTensor) -> DTensor:
        """
        Transform an atom-level single representation to a token-level single representation.

        Args:
            atom_single_repr: The atom-level single representation.
                Shape: (B, n_atoms, D) and placement: (Shard(0), Shard(1), Replicate())
            atom_to_token: The atom to token one-hot mapping except for padding atoms/tokens.
                Shape: (B, n_atoms, n_tokens_per_rank) and placement: (Shard(0), Shard(1), Replicate())

        Returns:
            The token-level single representation. Shape: (B, n_tokens, D)
        """
        single_repr_placements = (Shard(dim=0), Shard(dim=1), Replicate())  # same as atom_to_token placements
        if atom_to_token.placements != single_repr_placements:
            raise ValueError(
                f"Expect atom_to_token to have placements {single_repr_placements}, but got {atom_to_token.placements}"
            )
        if atom_single_repr.placements != single_repr_placements:
            raise ValueError(
                f"Expect atom_single_repr to have placements {single_repr_placements}, but got {atom_single_repr.placements}"
            )

        # TODO potential performance optimization by moving division after bmm

        # Normalize atom_to_token
        atom_to_token_local = atom_to_token.to_local().to(
            dtype=atom_single_repr.dtype
        )  # NOTE in case atom_to_token is int
        atom_to_token_sum = atom_to_token_local.sum(dim=1, keepdim=True).clamp(min=1)
        atom_to_token_mean = atom_to_token_local / atom_to_token_sum

        # Perform local bmm and distribute
        atom_single_repr_local = atom_single_repr.to_local()
        o = torch.einsum("bji,bj...->bi...", atom_to_token_mean, atom_single_repr_local)

        # Compute output shape and stride for the global DTensor
        # Output should be (B, n_tokens, D). By definition, atom_to_token.shape[2] == n_tokens_per_rank
        # which by definition guarantee atom -> token mapping is uniform across ranks
        # so n_tokens == n_tokens_per_rank * size_cp
        n_tokens = atom_to_token.shape[2] * atom_to_token.device_mesh.get_group(1).size()
        # where the last axis could be omitted if the input atom_single_repr is 2D
        shape_output = (atom_to_token.shape[0], n_tokens) + o.shape[2:]

        # Use LayoutRightMap for the output shape
        strides_output = update_exhaustive_strides(o.shape, o.stride(), shape_output)

        o = DTensor.from_local(
            o, atom_to_token.device_mesh, single_repr_placements, shape=shape_output, stride=strides_output
        )

        if atom_single_repr.requires_grad:
            ctx.device_mesh = atom_to_token.device_mesh
            ctx.single_repr_placements = single_repr_placements
            ctx.atom_single_repr_shape = atom_single_repr.shape
            ctx.atom_single_repr_stride = atom_single_repr.stride()
            ctx.save_for_backward(atom_to_token_mean)

        return o

    @staticmethod
    @torch.amp.custom_bwd(device_type="cuda")
    def backward(ctx, grad_output: DTensor) -> tuple[DTensor, None, None]:
        """
        Backward pass for single_repr_atom_to_token.
        """
        if grad_output.placements != (Shard(dim=0), Shard(dim=1), Replicate()):
            raise ValueError(
                f"Expect grad_output to have placements {(Shard(dim=0), Shard(dim=1), Replicate())}, but got {grad_output.placements}"
            )
        if grad_output.device_mesh != ctx.device_mesh:
            raise ValueError(
                f"Expect grad_output to have device mesh {ctx.device_mesh}, but got {grad_output.device_mesh}"
            )

        do = grad_output.to_local()
        (atom_to_token_mean,) = ctx.saved_tensors

        # Perform local bmm and distribute
        d_atom_single_repr = torch.einsum("bij,bj...->bi...", atom_to_token_mean, do)
        d_atom_single_repr = DTensor.from_local(
            d_atom_single_repr,
            ctx.device_mesh,
            ctx.single_repr_placements,
            shape=ctx.atom_single_repr_shape,
            stride=ctx.atom_single_repr_stride,
        )

        return d_atom_single_repr, None, None


class PairReprTokenToAtomFunction(torch.autograd.Function):
    """Autograd function for transforming token-level pair representation to atom-level pair representation."""

    @staticmethod
    def forward(
        ctx,
        token_repr: DTensor,
        atom_to_token: DTensor,
        transpose_comm: TransposeComm,
    ) -> DTensor:
        """
        Transform a token-level pair representation to an atom-level pair representation.

        Args:
            token_repr: The token-level pair representation. Shape: (B, n_tokens, n_tokens, D) and placement: (Shard(0), Shard(1), Shard(2))
            atom_to_token: The atom to token one-hot mapping except for padding atoms/tokens. Shape: (B, n_atoms, n_tokens) and placement: (Shard(0), Shard(1), Replicate())
            transpose_comm: The transpose communication object.

        Returns:
            The atom-level pair representation. Shape: (B, n_atoms, n_atoms, D)
        """
        single_repr_placements = (Shard(dim=0), Shard(dim=1), Replicate())  # same as atom_to_token placements
        if atom_to_token.placements != single_repr_placements:
            raise ValueError(
                f"Expect atom_to_token to have placements {single_repr_placements}, but got {atom_to_token.placements}"
            )
        pair_repr_placements = (Shard(dim=0), Shard(dim=1), Shard(dim=2))
        if token_repr.placements != pair_repr_placements:
            raise ValueError(
                f"Expect token_repr to have placements {pair_repr_placements}, but got {token_repr.placements}"
            )

        if atom_to_token.requires_grad:
            raise ValueError("atom_to_token should not require grad")

        # Perform transpose communication to get atom_to_token_local_j
        atom_to_token_local = atom_to_token.to_local().to(dtype=token_repr.dtype)  # NOTE in case atom_to_token is int
        atom_to_token_local = atom_to_token_local.contiguous()  # for both forward and backward
        atom_to_token_local_j = transpose_comm.enqueue_to_dispatch(atom_to_token_local)

        # Perform overlapped einsum operation
        token_repr_local = token_repr.to_local()

        # TODO potential performance optimization by op fusion versus communication overlap
        o = torch.einsum("bijd,bmi->bmjd", token_repr_local, atom_to_token_local)
        transpose_comm.wait_until_finished()
        o = torch.einsum("bmjd,bnj->bmnd", o, atom_to_token_local_j)

        # Compute output shape and stride for the global DTensor
        # Output should be (B, n_atoms, n_atoms, D)
        shape_output = (atom_to_token.shape[0], atom_to_token.shape[1], atom_to_token.shape[1]) + o.shape[3:]

        # Use LayoutRightMap for the output shape
        strides_output = update_exhaustive_strides(o.shape, o.stride(), shape_output)

        o = DTensor.from_local(
            o, token_repr.device_mesh, pair_repr_placements, shape=shape_output, stride=strides_output
        )

        # Save tensors needed for backward pass
        if token_repr.requires_grad:
            ctx.transpose_comm = transpose_comm
            ctx.device_mesh = token_repr.device_mesh
            ctx.pair_repr_placements = pair_repr_placements
            ctx.token_repr_shape = token_repr.shape
            ctx.token_repr_stride = token_repr.stride()
            ctx.save_for_backward(atom_to_token_local)

        return o

    @staticmethod
    def backward(ctx, grad_output: DTensor) -> tuple[DTensor, None, None, None]:
        """
        Backward pass for pair_repr_token_to_atom.

        Args:
            grad_output: Gradient w.r.t. output with shape (B, n_atoms, n_atoms, D)

        Returns:
            Tuple of gradients: (grad_token_repr, None, None)
            - grad_token_repr: Gradient w.r.t. token_repr
            - None: No gradient for atom_to_token (as specified)
            - None: No gradient for transpose_comm (not differentiable)
        """
        (atom_to_token_local,) = ctx.saved_tensors
        transpose_comm = ctx.transpose_comm
        atom_to_token_local_j = ctx.transpose_comm.enqueue_to_dispatch(atom_to_token_local)

        # Perform overlapped einsum operation
        do_local = grad_output.to_local()

        d_token_repr = torch.einsum("bmnd,bmi->bind", do_local, atom_to_token_local)
        transpose_comm.wait_until_finished()
        d_token_repr = torch.einsum("bind,bnj->bijd", d_token_repr, atom_to_token_local_j)
        d_token_repr = DTensor.from_local(
            d_token_repr,
            ctx.device_mesh,
            ctx.pair_repr_placements,
            shape=ctx.token_repr_shape,
            stride=ctx.token_repr_stride,
        )

        return d_token_repr, None, None, None


class SingleReprRepAtomToTokenFunction(torch.autograd.Function):
    """Autograd function for token representative-atom projection."""

    @staticmethod
    def forward(ctx, atom_single_repr: DTensor, token_to_rep_atom: DTensor) -> DTensor:
        """Project atom-level single representation to token-level via representative atoms.

        Supports a multiplicity factor: atom_single_repr may have shape (B*mult, n_atoms, D)
        while token_to_rep_atom has shape (B, n_tokens, n_atoms). The token_to_rep_atom map
        is broadcast across the multiplicity dimension.

        Args:
            atom_single_repr: Atom-level representation. Shape: (B*mult, n_atoms, D),
                placements: (Shard(0), Shard(1), Replicate()).
            token_to_rep_atom: Token->representative-atom one-hot map.
                Shape: (B, n_tokens, n_atoms), placements: (Shard(0), Shard(1), Replicate()).

        Returns:
            Token-level representation. Shape: (B*mult, n_tokens, D),
                placements: (Shard(0), Shard(1), Replicate()).
        """
        single_repr_placements = (Shard(0), Shard(1), Replicate())
        if atom_single_repr.placements != single_repr_placements:
            raise ValueError(
                f"Expect atom_single_repr to have placements {single_repr_placements}, but got {atom_single_repr.placements}"
            )
        if token_to_rep_atom.placements != single_repr_placements:
            raise ValueError(
                f"Expect token_to_rep_atom to have placements {single_repr_placements}, but got {token_to_rep_atom.placements}"
            )

        token_to_rep_local = token_to_rep_atom.to_local().to(dtype=atom_single_repr.dtype)
        atom_single_repr_local = atom_single_repr.to_local()

        # atom_single_repr may carry a multiplicity factor: (B*mult, N_atom, D) vs (B, N_token, N_atom).
        # Reshape to (B, mult, N_atom, D) so token_to_rep_local broadcasts over the mult dimension.
        B_local = token_to_rep_local.shape[0]
        mult = atom_single_repr_local.shape[0] // B_local
        atom_reshaped = atom_single_repr_local.reshape(B_local, mult, *atom_single_repr_local.shape[1:])
        # (B, N_token, N_atom) @ (B, mult, N_atom, D) -> (B, mult, N_token, D)
        o = torch.einsum("btj,bmj...->bmt...", token_to_rep_local, atom_reshaped)
        # Flatten back to (B*mult, N_token, D)
        o = o.reshape(B_local * mult, *o.shape[2:])

        shape_output = (token_to_rep_atom.shape[0] * mult,) + token_to_rep_atom.shape[1:2] + o.shape[2:]
        strides_output = update_exhaustive_strides(o.shape, o.stride(), shape_output)
        o = DTensor.from_local(
            o, token_to_rep_atom.device_mesh, single_repr_placements, shape=shape_output, stride=strides_output
        )

        if atom_single_repr.requires_grad:
            ctx.device_mesh = atom_single_repr.device_mesh
            ctx.single_repr_placements = single_repr_placements
            ctx.atom_single_repr_shape = atom_single_repr.shape
            ctx.atom_single_repr_stride = atom_single_repr.stride()
            ctx.mult = mult
            ctx.save_for_backward(token_to_rep_local)

        return o

    @staticmethod
    def backward(ctx, grad_output: DTensor) -> tuple[DTensor, None, None]:
        """Backward pass for representative atom projection."""
        if grad_output.placements != (Shard(0), Shard(1), Replicate()):
            raise ValueError(
                f"Expect grad_output to have placements {(Shard(0), Shard(1), Replicate())}, but got {grad_output.placements}"
            )
        if grad_output.device_mesh != ctx.device_mesh:
            raise ValueError(
                f"Expect grad_output to have device mesh {ctx.device_mesh}, but got {grad_output.device_mesh}"
            )

        do = grad_output.to_local()
        (token_to_rep_local,) = ctx.saved_tensors
        mult = ctx.mult
        B_local = token_to_rep_local.shape[0]
        # Reshape grad from (B*mult, N_token, D) to (B, mult, N_token, D)
        do_reshaped = do.reshape(B_local, mult, *do.shape[1:])
        # (B, N_token, N_atom)^T @ (B, mult, N_token, D) -> (B, mult, N_atom, D)
        d_atom_reshaped = torch.einsum("btj,bmt...->bmj...", token_to_rep_local, do_reshaped)
        # Flatten back to (B*mult, N_atom, D)
        d_atom_single_repr = d_atom_reshaped.reshape(B_local * mult, *d_atom_reshaped.shape[2:])
        d_atom_single_repr = DTensor.from_local(
            d_atom_single_repr,
            ctx.device_mesh,
            ctx.single_repr_placements,
            shape=ctx.atom_single_repr_shape,
            stride=ctx.atom_single_repr_stride,
        )
        return d_atom_single_repr, None, None


def single_repr_token_to_atom(
    token_single_repr: DTensor,
    atom_to_token: DTensor,
) -> DTensor:
    """
    Transform a token-level single representation to an atom-level single representation.

    Args:
        token_single_repr: The token-level single representation. Shape: (B, n_tokens, D) and placement: (Shard(0), Shard(1), Replicate())
        atom_to_token: The atom to token mapping. Shape: (B, n_tokens, n_atoms) and placement: (Shard(0), Shard(1), Replicate())
        device_mesh: The device mesh.

    Returns:
        The atom-level single representation. Shape: (B, n_atoms, D)
    """
    return SingleReprTokenToAtomFunction.apply(token_single_repr, atom_to_token)


def single_repr_atom_to_token(
    atom_single_repr: DTensor,
    atom_to_token: DTensor,
) -> DTensor:
    """
    Transform an atom-level single representation to a token-level single representation.

    Args:
        atom_single_repr: The atom-level single representation. Shape: (B, n_atoms, D) and placement: (Shard(0), Shard(1), Replicate())
        atom_to_token: The atom to token mapping. Shape: (B, n_atoms, n_tokens) and placement: (Shard(0), Shard(1), Replicate())
        device_mesh: The device mesh.

    Returns:
        The token-level single representation. Shape: (B, n_tokens, D)
    """
    return SingleReprAtomToTokenFunction.apply(atom_single_repr, atom_to_token)


def pair_repr_token_to_atom(
    token_repr: DTensor,
    atom_to_token: DTensor,
    transpose_comm: TransposeComm,
) -> DTensor:
    """
    Transform a token-level pair representation to an atom-level pair representation.

    Args:
        token_repr: The token-level pair representation. Shape: (B, n_tokens, n_tokens, D) and placement: (Shard(0), Shard(1), Shard(2))
        atom_to_token: The atom to token mapping. Shape: (B, n_tokens, n_atoms) and placement: (Shard(0), Shard(1), Replicate())
        transpose_comm: The transpose communication object.

    Returns:
        The atom-level pair representation. Shape: (B, n_atoms, n_atoms, D)
    """
    return PairReprTokenToAtomFunction.apply(token_repr, atom_to_token, transpose_comm)


def single_repr_rep_atom_to_token(
    atom_single_repr: DTensor,
    token_to_rep_atom: DTensor,
) -> DTensor:
    """Project atom-level single representation to token-level using representative atoms."""
    return SingleReprRepAtomToTokenFunction.apply(atom_single_repr, token_to_rep_atom)


def _reconstruct_onehot_diag_block_global(
    dtensor: DTensor,
) -> Tensor:
    """Reconstruct a diagonally-sharded one-hot DTensor into a full plain tensor.

    Both ``atom_to_token`` (atoms×tokens) and ``token_to_rep_atom`` (tokens×atoms)
    use diagonal block sharding: shard *i* only contains the non-zero block
    relating atoms of shard *i* to tokens of shard *i*.  After all-reduce the
    complete global matrix is recovered.

    Parameters
    ----------
    dtensor : DTensor
        Diagonally-sharded DTensor with placements ``(Shard(0), Shard(1), Replicate())``.

    Returns
    -------
    Tensor
        Reconstructed global tensor of shape ``(B, N_atoms_global, N_tokens_global)``
        or ``(B, N_tokens_global, N_atoms_global)``.
    """
    device_mesh = dtensor.device_mesh

    expected_placements = (Shard(dim=0), Shard(dim=1), Replicate())
    if dtensor.placements != expected_placements:
        raise ValueError(f"Expected placements {expected_placements}, got {dtensor.placements}")

    for i_dim_mesh, placement in enumerate(expected_placements):
        if isinstance(placement, Shard) and dtensor.shape[placement.dim] % device_mesh.shape[i_dim_mesh] != 0:
            raise ValueError(
                f"Uneven sharding tensor dimension {placement.dim} of size {dtensor.shape[placement.dim]} "
                f"along device mesh dimension {i_dim_mesh} of size {device_mesh.shape[i_dim_mesh]} is not supported"
            )

    local = dtensor.to_local()
    B = local.shape[0]
    assert B == 1, "Only batch size 1 is supported"

    n_per_shard_dim1 = local.shape[1]
    n_per_shard_dim2 = local.shape[2]

    cp_axis_0_size = device_mesh.get_group("cp_axis_0").size()
    cp_axis_0_rank, cp_axis_1_rank = device_mesh.get_coordinate()[1:]

    n_global_dim1 = n_per_shard_dim1 * cp_axis_0_size
    n_global_dim2 = n_per_shard_dim2 * cp_axis_0_size

    result = torch.zeros(
        B,
        n_global_dim1,
        n_global_dim2,
        dtype=local.dtype,
        device=local.device,
    )

    start_dim1 = cp_axis_0_rank * n_per_shard_dim1
    end_dim1 = start_dim1 + n_per_shard_dim1
    start_dim2 = cp_axis_0_rank * n_per_shard_dim2
    end_dim2 = start_dim2 + n_per_shard_dim2

    if cp_axis_1_rank == 0:
        n_non_zeros = local.sum(dim=2)
        if not ((n_non_zeros == 0) | (n_non_zeros == 1)).all():
            raise ValueError(
                f"Input DTensor shard is not one-hot for CP rank ({cp_axis_0_rank}, {cp_axis_1_rank}): "
                f"found rows with sum not in {{0, 1}}"
            )
        result[:, start_dim1:end_dim1, start_dim2:end_dim2] = local

    torch.distributed.all_reduce(result, op=torch.distributed.ReduceOp.SUM, group=device_mesh.get_group("cp_axis_0"))
    torch.distributed.all_reduce(result, op=torch.distributed.ReduceOp.SUM, group=device_mesh.get_group("cp_axis_1"))

    return result


def reconstruct_token_to_rep_atom_global(token_to_rep_atom_dtensor: DTensor) -> Tensor:
    """Reconstruct the full ``token_to_rep_atom`` matrix from a diagonally-sharded DTensor.

    The reconstruction mirrors :func:`reconstruct_atom_to_token_global` but for
    the transposed mapping ``(B, N_tokens, N_atoms)``.

    Parameters
    ----------
    token_to_rep_atom_dtensor : DTensor
        Diagonally-sharded DTensor with placements ``(Shard(0), Shard(1), Replicate())``.
        Local shape: ``(B, N_tokens_per_shard, max_atoms_per_shard)``.

    Returns
    -------
    Tensor
        ``(B, N_tokens_global, N_atoms_global)``
    """
    return _reconstruct_onehot_diag_block_global(token_to_rep_atom_dtensor)


def reconstruct_r_set_to_rep_atom_global(r_set_dtensor: DTensor) -> Tensor:
    """Reconstruct the full ``r_set_to_rep_atom`` matrix from a diagonally-sharded DTensor.

    Parameters
    ----------
    r_set_dtensor : DTensor
        Diagonally-sharded DTensor with placements ``(Shard(0), Shard(1), Replicate())``.
        Local shape: ``(B, max_r_set_per_shard, max_atoms_per_shard)``.

    Returns
    -------
    Tensor
        ``(B, N_R_global, N_atoms_global)``
    """
    return _reconstruct_onehot_diag_block_global(r_set_dtensor)


def reconstruct_atom_to_token_global(atom_to_token_dtensor: DTensor) -> Tensor:
    """
    Reconstruct the original full atom_to_token tensor from a DTensor with (Shard, Shard, Replicate) placements.

    This function reverses the context parallel sharding strategy by:
    1. Gathering the local tensor from each rank
    2. Reconstructing the block diagonal structure
    3. Removing padding to get the original tensor

    Args:
        atom_to_token_dtensor: DTensor with placements (Shard(0), Shard(1), Replicate())
                               Shape: (global_batch_size, n_atoms_per_rank, n_tokens)

    Returns:
        Tensor: The reconstructed global atom_to_token tensor
                Shape: (local_batch_size, n_atoms_global, n_tokens_global)
    """
    result = _reconstruct_onehot_diag_block_global(atom_to_token_dtensor)
    assert torch.max(result) == 1
    return result
