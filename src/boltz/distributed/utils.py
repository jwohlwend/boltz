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


from math import isqrt
from typing import Callable, Iterable, List, Self, Sequence

import numpy as np
import torch
from torch import Tensor
from torch.distributed.tensor import DTensor, Replicate, Shard
from torch.distributed.tensor.device_mesh import DeviceMesh


class LayoutMap:
    """
    A class representing a mapping between multidimensional indices and flat indices.
    This is basing on the C++ std::layout_stride::mapping but adds on top it the mapping
    from the flat index to the multidimensional indices.

    Parameters
    ----------
    strides : tuple of ints
        The strides of the layout.
    shape : tuple of ints, optional
        The shape of the layout. If not provided, the shape will be inferred from the strides.

    Raises
    ------
    ValueError
        If the input strides or shape is invalid.

    Attributes
    ----------
    numel : int
        The total number of elements in the layout.
    shape : tuple of ints
        The shape of the layout.
    """

    def __init__(
        self,
        strides: tuple[int, ...],
        shape: tuple[int, ...],
        offset: int = 0,
    ):
        """
        Initialize the layout mapping.

        Parameters
        ----------
        strides : tuple of ints
            The strides of the layout.
        shape : tuple of ints
            The shape of the layout.
        offset : int, optional
            The offset of the layout.

        Raises
        ------
        ValueError
            If the input strides or shape is invalid.

        Notes
        -----
        The input strides must be a cumulative product of some permutation of the input shape.
        """
        if not all(isinstance(stride, (int, np.int64)) and stride > 0 for stride in strides):
            raise ValueError(f"Input strides contain non-integer or negative values: {strides}")

        self._has_negative_shape = any(s < 0 for s in shape)

        if self._has_negative_shape:
            raise ValueError(f"Input shape contain negative values: {shape}")

        self._has_zero_shape = any(s == 0 for s in shape)

        if self._has_zero_shape:
            # NOTE: the C++ standard does allow the shape to be zero along some axes, which
            # is not useful for our usage case. We nonetheless can relax the condition to
            # allow zero-sized axes but it requires more testing
            raise ValueError(f"Input shape contain zero values: {shape}")

        self._strides = strides
        self._n_axes = len(strides)

        if len(shape) != self._n_axes:
            raise ValueError(f"Shape {shape} and strides {strides} must have the same length")

        self._shape = shape
        self._numel = np.prod(self._shape)
        self._offset = offset

        # singleton axes can confound the uniqueness and exhaustiveness check, e.g.,
        # for layout right of (3, 1, 5), the strides are (5, 5, 1) but direct argsort
        # on the strides will give the permuted exhaustive stride of (1, 5, 15), which
        # corresponds to the strides of (5, 15, 1) (argsort of (2, 0, 1)), which will
        # fail the uniqueness check. This is purely artifact of the stable sorting where
        # the single axis can potentially be arbitrarily placed before or after the
        # other axes with the same stride. The correct thing to do is to handle the ties
        # involving the singleton axes so that we sort by shape if two stride elements are
        # tied.
        shape_and_strides = np.array(
            list(zip(self._shape, self._strides)),
            dtype=np.dtype([("shape", int), ("strides", int)]),
        )
        argsort_ascend_strides_and_shape = np.argsort(shape_and_strides, order=["strides", "shape"])

        self.is_unique = self._is_unique(argsort_ascend_strides_and_shape)
        self.is_exhaustive = self._is_exhaustive(argsort_ascend_strides_and_shape)

        if not self.is_unique:
            raise ValueError(f"Input strides {strides} and shape {shape} do not give unique layout.")

        self._required_span_size = self._compute_required_span_size()
        self._argsort_descend_strides = argsort_ascend_strides_and_shape[::-1]
        self._argsort_ascend_strides = argsort_ascend_strides_and_shape

    def _compute_required_span_size(self) -> int:
        """
        Calculate the minimal span size required to represent the layout, e.g.,
        by storing the elements in a contiguous piece of memory

        See the C++ std::layout_stride's requirements here:
        https://eel.is/c++draft/views.multidim#mdspan.layout.stride.expo-1

        """
        if self._n_axes == 0:
            return 1
        if self._has_zero_shape:
            return 0
        return 1 + sum((self._shape[i] - 1) * self._strides[i] for i in range(self._n_axes))

    def _strides_exhaustive(self, permutation: np.ndarray) -> np.ndarray:
        """
        Calculate the expected exhaustive strides for a given permutation.

        For a valid exhaustive layout, the strides should follow a pattern where
        strides[p[i]] equals strides[p[i-1]] * shape[p[i-1]] for all i > 0,
        and strides[p[0]] equals 1, where p is the permutation of indices.

        Parameters
        ----------
        permutation : np.ndarray
            Permutation of indices in ascending order of strides.

        Returns
        -------
        tuple[np.ndarray, np.ndarray]
            A tuple containing:
            - The permuted strides
            - The expected exhaustive strides for the given permutation
        """
        strides = np.array(self._strides)
        shape = np.array(self._shape)
        shape_permuted = shape[permutation]
        strides_permuted = strides[permutation]
        shape_shifted = np.concatenate([[1], shape_permuted[:-1]])
        strides_shifted = np.concatenate([[1], strides_permuted[:-1]])
        return strides_permuted, strides_shifted * shape_shifted

    def _is_unique(self, permutation: np.ndarray) -> bool:
        """
        Check if the layout mapping is unique, i.e.,
        if the mapping: n-dimensional index -> flat index is injective.

        This implements the necessary and sufficient condition for a
        LayoutMap, i.e., the equivalent std::layout_stride::mapping, to
        be unique.

        See the C++ std::layout_stride's requirements here for the
        uniqueness condition:
        https://eel.is/c++draft/views.multidim#mdspan.layout.stride.cons

        Parameters
        ----------
        permutation : np.ndarray
            Permutation of indices per the standard's requirement. Note that
            the C++ standard doesn't specify what the permutation should be
            as soon as it exists. In practice, one can construct a proof by
            induction that the permutation can be chosen to be the ascending
            ordering of the strides.

        Returns
        -------
        bool
            True if the layout has unique strides, False otherwise.
        """
        if self._n_axes == 0:
            return True
        strides, strides_exhaustive = self._strides_exhaustive(permutation)
        ans = np.all(strides >= strides_exhaustive)
        return ans

    def _is_exhaustive(self, permutation: np.ndarray) -> bool:
        """
        Check if the layout mapping is exhaustive, i.e.,
        if the mapping: n-dimensional index -> flat index is surjective.

        This implements the necessary and sufficient condition for a
        LayoutMap, i.e., the equivalent std::layout_stride::mapping, to
        be exhaustive.

        See the C++ std::layout_stride's requirements here for the
        exhaustiveness condition:
        https://eel.is/c++draft/views.multidim#mdspan.layout.stride.obs-5.2

        Parameters
        ----------
        permutation : np.ndarray
            Permutation of indices per the standard's requirement. Note that
            the C++ standard doesn't specify what the permutation should be
            as soon as it exists. In practice, one can construct a proof by
            induction that the permutation can be chosen to be the ascending
            ordering of the strides.

        Returns
        -------
        bool
            True if the layout is exhaustive, False otherwise.
        """
        if self._n_axes == 0:
            return True
        strides, strides_exhaustive = self._strides_exhaustive(permutation)
        ans = np.all(strides == strides_exhaustive)
        return ans

    @property
    def offset(self) -> int:
        """
        Get the offset of the layout.
        """
        return self._offset

    @property
    def required_span_size(self) -> int:
        """
        Get the required span size of the layout.
        """
        return self._required_span_size

    @property
    def numel(self) -> int:
        """
        Get the total number of elements in the layout.
        Returns
        -------
        int
            The total number of elements.
        """
        return self._numel

    @property
    def shape(self) -> tuple[int, ...]:
        """
        Get the shape of the layout.

        Returns
        -------
        tuple of ints
            The shape of the layout.
        """
        return self._shape

    @property
    def strides(self) -> tuple[int, ...]:
        """
        Get the strides of the layout.

        Returns
        -------
        tuple of ints
            The strides of the layout.

        Notes
        -----
        The strides are a cumulative product of some permutation of the input shape.
        """
        return self._strides

    def __call__(self, ids: tuple[int, ...]) -> int:
        """
        Get the flat index corresponding to the given multidimensional index.

        Parameters
        ----------
        ids : tuple of ints
            The multidimensional index.

        Returns
        -------
        int
            The flat index.

        Raises
        ------
        ValueError
            If the input index is out of range.
        """
        if len(ids) != self._n_axes:
            raise ValueError(f"Expected {self._n_axes} elements in ids but got only {len(ids)}")

        if len(ids) == 0:
            return self._offset

        if self._shape is not None:
            for axis, idx in enumerate(ids):
                if idx < 0 or idx >= self._shape[axis]:
                    raise ValueError(
                        f"Expected ids to satisfy 0 <= ids[{axis}] <= {self._shape[axis] - 1} "
                        f"but found ids[{axis}] == {idx}"
                    )
        return np.dot(ids, self._strides) + self._offset

    def unravel(self, flat_index: int) -> tuple[int, ...]:
        """
        Convert a flat index to a multidimensional index.

        Parameters
        ----------
        flat_index : int
            The flat index.

        Returns
        -------
        tuple of ints
            The multidimensional index.

        Raises
        ------
        TypeError
            If the input is not an integer.
        ValueError
            If the input index is out of range.
        """
        if not self.is_unique:
            # double check the uniqueness of the layout
            raise ValueError(f"Layout is not unique, cannot unravel {flat_index}")

        if not isinstance(flat_index, (int, np.integer)):
            raise TypeError(f"Expected arg to be an int, but instead got type {type(flat_index)}")

        remaining = flat_index - self._offset

        if remaining < 0 or remaining >= self._required_span_size:
            raise ValueError(
                f"Expected flat_index in range [{self._offset}, {self._offset + self._required_span_size - 1}], "
                f"but instead got {flat_index}"
            )

        indices = [0] * self._n_axes  # Initialize indices

        for i_dim in self._argsort_descend_strides:
            stride = self._strides[i_dim]
            size = self._shape[i_dim]
            indices[i_dim] = (remaining // stride) % size
            remaining -= indices[i_dim] * stride

        if remaining != 0:
            msg = f"Input flat_index {flat_index} is out of the valid range of span."
            if not self.is_exhaustive:
                msg += " Given the layout is not exhaustive, the input flat_index can fall into the unmapped region."
            raise ValueError(msg)

        return tuple(indices)

    def __getitem__(self, slices: tuple[slice | int, ...]) -> Self:
        """
        Create a new LayoutMap by slicing the current layout along specified dimensions.

        This method allows for creating sub-layouts by slicing the original layout,
        similar to numpy array slicing. Dimensions can be reduced by using integer
        indices, or transformed by using slices with custom start, stop, and step values.

        Parameters
        ----------
        slices : tuple[slice | int, ...] or slice or int
            Slices to apply along each dimension. Can be a single slice/int or a tuple
            of slices/integers. Integer indices collapse the corresponding dimension,
            while slices transform the dimension according to start, stop, and step.
            If fewer slices are provided than dimensions, remaining dimensions will
            be sliced with full-range slices (equivalent to ':').

        Returns
        -------
        Self
            A new LayoutMap instance with updated shape, strides and offset reflecting
            the applied slicing.

        Raises
        ------
        ValueError
            If slicing with negative or zero steps, or if start is less than or equal to stop.
        TypeError
            If slice elements are not of type slice or int.

        Examples
        --------
        >>> layout = LayoutMap((12, 4, 1), (2, 3, 4))
        >>> # Slice first dimension from indices 1 to 3 with step 2
        >>> sub_layout = layout[slice(1, 3, 2), :, :]
        >>> # Collapse second dimension by selecting index 1
        >>> sub_layout = layout[:, 1, :]
        >>> # Slice only first dimension, remaining dimensions use full range
        >>> sub_layout = layout[1:]  # Equivalent to layout[1:, :, :]
        """
        if not isinstance(slices, tuple) and (isinstance(slices, slice) or isinstance(slices, int)):
            slices = (slices,)

        # Pad slices with full-range slices if needed
        if len(slices) < self._n_axes:
            full_slice = slice(None)  # This is equivalent to ':'
            slices = slices + (full_slice,) * (self._n_axes - len(slices))

        new_shape = []
        new_strides = []
        new_offset = self.offset

        for axis, s in enumerate(slices):
            if isinstance(s, (int, np.int64)):
                # Collapse dimension and adjust offset
                new_offset += s * self.strides[axis]
            elif isinstance(s, slice):
                start, stop, step = s.indices(self.shape[axis])
                if step <= 0:
                    raise ValueError("Unsupported slicing: Negative or zero steps")
                if start >= stop:
                    # NOTE: the start == stop could be supported because we could have
                    # a layout with shape[i] == 0. But it wouldn't be useful for our usage cases.
                    raise ValueError("Unsupported slicing: start not smaller than stop")

                # Calculate new dimension length
                dim_len = (stop - start + step - 1) // step
                dim_len = max(0, dim_len)

                # Update metadata
                new_shape.append(dim_len)
                new_strides.append(self.strides[axis] * step)
                new_offset += start * self.strides[axis]
            else:
                raise TypeError(f"Unsupported slice type: {type(s)}")

        return LayoutMap(tuple(new_strides), tuple(new_shape), new_offset)


class LayoutRightMap(LayoutMap):
    """
    A class representing a right-aligned layout mapping.

    Parameters
    ----------
    shape : tuple of ints
        The shape of the layout.
    """

    def __init__(self, shape: tuple[int, ...]):
        """
        Initialize the layout mapping.

        Parameters
        ----------
        shape : tuple of ints
            The shape of the layout.
        """
        strides = np.ones_like(shape)
        strides[1:] = shape[:0:-1]
        strides = np.cumprod(strides)[::-1]
        super().__init__(tuple(strides), shape=shape)


class LayoutLeftMap(LayoutMap):
    """
    A class representing a left-aligned layout mapping.

    Parameters
    ----------
    shape : tuple of ints
        The shape of the layout.
    """

    def __init__(self, shape: tuple[int, ...]):
        """
        Initialize the layout mapping.

        Parameters
        ----------
        shape : tuple of ints
            The shape of the layout.
        """
        strides = np.ones_like(shape)
        strides[1:] = shape[:-1]
        strides = np.cumprod(strides)
        super().__init__(tuple(strides), shape=shape)


def slice_repr_mask(
    s: torch.Tensor,
    z: torch.Tensor,
    mask: torch.Tensor,
    pair_mask: torch.Tensor,
    n_ranks: int,
    layout_group: LayoutMap,
) -> tuple[List[torch.Tensor], List[torch.Tensor], List[torch.Tensor], List[torch.Tensor]]:
    """
    Slice the single- and pair-representation tensors and masks along the given coordinates.

    This function slices the input tensors into n_ranks segments along the token dimensions,
    with the segment boundaries determined by the layout mapping of the ranks.

    Parameters
    ----------
    s : torch.Tensor
        Single representation tensor.
    z : torch.Tensor
        Pair representation tensor.
    mask : torch.Tensor
        Single representation mask.
    pair_mask : torch.Tensor
        Pair representation mask.
    n_ranks : int
        The number of ranks.
    layout_group : LayoutMap
        The layout mapping of the ranks.
    Returns
    -------
    tuple[List[torch.Tensor], List[torch.Tensor], List[torch.Tensor], List[torch.Tensor]]
        Sliced single representation tensors, pair representation tensors, single representation masks,
        and pair representation masks.

    Raises
    ------
    ValueError
        If the shape of any input tensor is incompatible with the layout mapping.
    """
    if z.shape[-2] != z.shape[-3]:
        raise ValueError(f"z is not square tensor in the middle two axes but of shape {z.shape}")
    if s.shape[-2] != z.shape[-3]:
        raise ValueError(f"Incompatible s shape {s.shape} and z shape {z.shape}")
    if mask.shape != s.shape[:-1]:
        raise ValueError(f"Incompatible s shape {s.shape} and mask shape {mask.shape}")
    if pair_mask.shape != z.shape[:-1]:
        raise ValueError(f"Incompatible z shape {z.shape} and pair_mask shape {pair_mask.shape}")
    n_tokens = s.shape[-2]
    coords = [layout_group.unravel(rank) for rank in range(n_ranks)]
    n_ranks_axis = isqrt(n_ranks)
    if n_ranks_axis * n_ranks_axis != n_ranks:
        raise ValueError(f"Input n_ranks is not a square int: {n_ranks}")
    if n_tokens % n_ranks_axis:
        raise ValueError(
            f"Input tensors size along the token dimensions {n_tokens} not divisible by square root of n_ranks {n_ranks}"
        )
    stride = n_tokens // n_ranks_axis
    s_slices = []
    z_slices = []
    mask_slices = []
    pair_mask_slices = []
    for i_row, j_col in coords:
        i_row_begin = i_row * stride
        i_row_end = (i_row + 1) * stride
        j_col_begin = j_col * stride
        j_col_end = (j_col + 1) * stride
        s_slices.append(s[..., i_row_begin:i_row_end, :].contiguous())
        mask_slices.append(mask[..., i_row_begin:i_row_end].contiguous())
        z_slices.append(
            z[
                ...,
                i_row_begin:i_row_end,
                j_col_begin:j_col_end,
                :,
            ].contiguous()
        )
        pair_mask_slices.append(
            pair_mask[
                ...,
                i_row_begin:i_row_end,
                j_col_begin:j_col_end,
            ].contiguous()
        )
    return s_slices, z_slices, mask_slices, pair_mask_slices


def gather_repr(
    s_slices: List[torch.Tensor],
    z_slices: List[torch.Tensor],
    rank: int,
    group: torch.distributed.ProcessGroup,
    layout_group: LayoutMap,
    s_global: torch.Tensor,
    z_global: torch.Tensor,
) -> tuple[torch.Tensor, torch.Tensor]:
    """
    Gather representation slices from all ranks in the group and fill the global representations.

    This function performs all-gather operations on the representation slices from each rank,
    and fills the global representations `s_global` and `z_global`.
    The single representation `s_global` is assumed to be sliced along the rows of the grid of process ranks,
    i.e., rank[i, :] owns s_global[..., i, :].
    The pair representation `z_global` is sliced along both of its middle axes.

    Parameters
    ----------
    s_slices : List[torch.Tensor]
        Slices of single representation `s` from each rank of the group.
    z_slices : List[torch.Tensor]
        Slices of pair representation `z` from each rank of the group.
    rank : int
        The rank of the current process.
    group : torch.distributed.ProcessGroup
        The process group for all-gather operation.
    layout_group : LayoutMap
        The layout mapping of the group.
    s_global : torch.Tensor
        The global representation `s` to be filled.
    z_global : torch.Tensor
        The global representation `z` to be filled.

    Returns
    -------
    tuple[torch.Tensor, torch.Tensor]
        The gathered global representations `s_global` and `z_global`.

    Raises
    ------
    ValueError
        If the input tensors have incompatible shapes.
    """
    if z_global.shape[-2] != z_global.shape[-3]:
        raise ValueError(f"z_global is not square tensor in the middle two axes but of shape {z_global.shape}")
    if s_global.shape[-2] != z_global.shape[-3]:
        raise ValueError(f"Incompatible s_global shape {s_global.shape} and z_global shape {z_global.shape}")
    n_tokens = s_global.shape[-2]
    ranks = torch.distributed.get_process_group_ranks(group)
    n_ranks = len(ranks)
    coords = [layout_group.unravel(rank) for rank in range(n_ranks)]
    n_ranks_axis = isqrt(n_ranks)
    if n_ranks_axis * n_ranks_axis != n_ranks:
        raise ValueError(f"Input n_ranks is not a square int: {n_ranks}")
    if n_tokens % n_ranks_axis:
        raise ValueError(
            f"Input tensors size along the token dimensions {n_tokens} not divisible by square root of n_ranks {n_ranks}"
        )
    stride = n_tokens // n_ranks_axis
    req_gather_z = torch.distributed.all_gather(z_slices, z_slices[rank], group=group, async_op=True)
    req_gather_s = torch.distributed.all_gather(s_slices, s_slices[rank], group=group, async_op=True)
    req_gather_s.wait()
    for i_rank in range(len(coords)):
        i_row, j_col = coords[i_rank]
        s_global[..., i_row * stride : (i_row + 1) * stride, :] = s_slices[i_rank]
    req_gather_z.wait()
    for i_rank in range(len(coords)):
        i_row, j_col = coords[i_rank]
        z_global[..., i_row * stride : (i_row + 1) * stride, j_col * stride : (j_col + 1) * stride, :] = z_slices[
            i_rank
        ]
    return s_global, z_global


def get_group_rank_from_axial_shift(coord: tuple[int, ...], axis: int, delta: int, layout_group: LayoutMap) -> int:
    """
    Get the rank of a process after shifting its coordinates along an axis.

    Parameters
    ----------
    coord : tuple of ints
        The current coordinates of the process in the group layout.
    axis : int
        The axis along which to shift the coordinates.
    delta : int
        The amount to shift the coordinates by (can be positive or negative).
    layout_group : LayoutMap
        The layout mapping of the process group.

    Returns
    -------
    int
        The rank of the process after shifting its coordinates.

    Raises
    ------
    ValueError
        If the coordinates are incompatible with the layout shape or if the axis is out of range.
    """
    if len(coord) != len(layout_group.shape):
        raise ValueError(f"Incompatible coord {coord} and layout_group shape {layout_group.shape}")
    if axis >= len(coord):
        raise ValueError(f"Axis {axis} is out of range for coord {coord}")
    coord_shifted = list(coord)
    coord_shifted[axis] = (coord_shifted[axis] + delta) % layout_group.shape[axis]
    return layout_group(coord_shifted)


def all_reduce_weighted_mean(
    weights,
    values,
    group_reduce: torch.distributed.ProcessGroup,
    dim: int | tuple[int, ...] = -1,
    eps: float = 0.0,
) -> Tensor:
    """Perform distributed weighted sum operation.

    Args:
        weights: weights tensor
        values: values tensor
        group_reduce: process group for reduction
        dim: dimension to perform weighted mean operation on; default is -1
        eps: epsilon value to avoid division by zero; default is 0.0

    Returns:
        Tensor: weighted mean of values
    """
    values_local = torch.sum(weights * values, dim=dim)
    values_work = torch.distributed.all_reduce(
        values_local,
        op=torch.distributed.ReduceOp.SUM,
        group=group_reduce,
        async_op=True,
    )
    weights_local = torch.sum(weights, dim=dim) + eps
    torch.distributed.all_reduce(weights_local, op=torch.distributed.ReduceOp.SUM, group=group_reduce)
    values_work.wait()
    return values_local / weights_local


def tiled_softmax_attention_update(
    o_chunk: torch.Tensor,
    lse_m_chunk: torch.Tensor,
    amax_chunk: torch.Tensor | None,
    o: torch.Tensor | None = None,
    lse_m: torch.Tensor | None = None,
    amax: torch.Tensor | None = None,
):
    """
    Update online softmax attention accumulation with a new chunk of data.

    This function implements a numerically stable online softmax computation that processes
    data in chunks. It maintains running statistics (output accumulation, log-sum-exp, and
    maximum values) to compute the final softmax result incrementally without storing all
    intermediate values in memory.

    The algorithm is particularly useful for attention mechanisms where the sequence length
    is too large to fit in memory at once, allowing for tiled/chunked processing while
    maintaining mathematical equivalence to full softmax computation.

    Args:
        o_chunk (torch.Tensor): Output accumulation for the current chunk.
            Shape: (..., D) where D is the feature dimension.
        lse_m_chunk (torch.Tensor): Log-sum-exp minus max for the current chunk.
            Shape: (..., 1) - must have last dimension of size 1.
            Note: When amax_chunk is None, this is actually lse (not lse_m = lse - amax).
        amax_chunk (torch.Tensor | None): Maximum value for the current chunk.
            Shape: (..., 1) - must match lse_m_chunk shape.
            If None, the function operates without amax tracking (has_amax == False mode).
            See Note below for implications.
        o (torch.Tensor | None, optional): Accumulated output from previous chunks.
            If None, this is treated as the first chunk. Must have same shape as o_chunk
            if provided. Defaults to None.
        lse_m (torch.Tensor | None, optional): Accumulated log-sum-exp minus max from
            previous chunks. Must have same shape as lse_m_chunk if provided.
            Note: When amax is None, this is actually lse (not lse_m = lse - amax).
            Defaults to None.
        amax (torch.Tensor | None, optional): Maximum value across all previous chunks.
            Must have same shape as amax_chunk if provided. If amax_chunk is None,
            this must also be None, and the returned amax will also be None. Defaults to None.

    Returns:
        tuple[torch.Tensor, torch.Tensor, torch.Tensor]: A tuple containing:
            - o (torch.Tensor): Updated accumulated output with current chunk incorporated.
            - lse_m (torch.Tensor): Updated log-sum-exp minus max value.
            - amax (torch.Tensor): Updated maximum value across all processed chunks.

    Raises:
        ValueError: If input tensors have incompatible shapes or if the None/not-None
            consistency of o, lse_m, and amax is violated.

        Mathematical Derivation:
        The online softmax update rules are derived as follows:

        Given:
        - amax: cumulative max until current block
        - amax_chunk: max of current block
        - amax_next: resulting amax post-update = max(amax, amax_chunk)
        - lse_chunk = log(sum(exp(a_chunk)))
        - lse_m_chunk = log(sum(exp(a_chunk - amax_chunk))) = lse_chunk - amax_chunk
        - lse = log(sum(exp(a)))
        - lse_m = log(sum(exp(a - amax))) = lse - amax
        - out = unormalized_out / exp(lse) = unormalized_out / exp(lse_m + amax)

        Update rule for lse_next:
        lse_next = lse + torch.log(1 + torch.exp(lse_chunk - lse))
                 = lse - F.logsigmoid(lse - lse_chunk)

        Update rule for lse_m_next:
        lse_m_next = lse_next - amax_next
                   = lse - F.logsigmoid(lse - lse_chunk) - amax_next
                   = lse - (lse - lse_chunk + F.logsigmoid(lse_chunk - lse)) - amax_next
                   = lse_chunk - amax_next - F.logsigmoid(lse_chunk - lse)
                   = lse_m_chunk + amax_chunk - amax_next - F.logsigmoid(lse_m_chunk + amax_chunk - lse_m - amax)
                   = lse_m_chunk + amax_chunk - amax_next - F.logsigmoid(lse_m_chunk - lse_m + amax_chunk - amax)
                   = lse_m_chunk + amax_chunk - amax_next - (amax_chunk - amax - logsumexp([-lse_m_chunk + lse_m, amax_chunk - amax]))
                   = lse_m_chunk + amax - amax_next + logsumexp([-lse_m_chunk + lse_m, amax_chunk - amax])
                   = lse_m_chunk + logsumexp([amax - amax_next - lse_m_chunk + lse_m, amax_chunk - amax_next])

        Update rule for o_next (numerically stable form):
        o_next = torch.exp(lse - lse_next) * o + torch.exp(lse_chunk - lse_next) * o_chunk

        The following computation is more numerically stable:
        o_next = o / (1 + torch.exp(lse_chunk - lse)) + torch.exp(lse_chunk - lse + F.logsigmoid(lse - lse_chunk)) * o_chunk
               = o * F.sigmoid(lse - lse_chunk) + torch.exp(lse_chunk - lse) * F.sigmoid(lse - lse_chunk) * o_chunk
               = o * F.sigmoid(lse - lse_chunk) + F.sigmoid(lse_chunk - lse) * o_chunk
               = o * (1 - F.sigmoid(lse_chunk - lse)) + F.sigmoid(lse_chunk - lse) * o_chunk
               = o - F.sigmoid(lse_chunk - lse) * (o - o_chunk)
               = o - F.sigmoid(lse_m_chunk + amax_chunk - lse_m - amax) * (o - o_chunk)

    Note:
        All three optional parameters (o, lse_m, amax) must be either all None (indicating
        this is the first chunk) or all provided with compatible shapes. The function uses
        numerically stable computations with log-space arithmetic to prevent overflow/underflow
        issues common in softmax computations.

        The mathematical operations maintain the invariant that the final result is equivalent
        to computing softmax over the concatenation of all processed chunks.

        **Behavior when amax_chunk is None (has_amax == False mode):**
        When amax_chunk (and correspondingly amax) is None, the function operates in a mode
        where:

        a) lse_m is actually lse (i.e., we assume lse_m = lse and don't track amax separately)
        b) The returned amax will always be None throughout all chunk updates
        c) The update equations have **limited dynamic range** due to catastrophic cancellation
           in the computation of d_lse_m = lse_m - lse_m_chunk. Without amax to handle extreme
           values, this subtraction can lose precision when values differ significantly.

        Despite the limited dynamic range, this mode still works correctly in common use cases
        where "-inf" padding in attention scores appears only after normal values (not before).
        In such cases, the lse pattern is [...normal values..., -inf, -inf, ..., -inf], and
        the trailing -inf values are effectively discarded because sigmoid(-inf) ≈ 0 and
        logsigmoid(inf) ≈ 0, even though -inf + normal_value = -inf arithmetically.

        For better numerical stability across all scenarios, it is recommended to use
        has_amax == True mode by providing amax_chunk.
    """
    if not ((o is None) == (lse_m is None)):
        raise ValueError("o and lse_m must both be None or both be not None")

    # has_amax == False will ignored amax terms entirely but assumes lse_m = lse_m + amax
    has_amax = amax_chunk is not None

    is_initial_chunk = o is None

    if has_amax and lse_m_chunk.shape != amax_chunk.shape:
        raise ValueError("lse_m_chunk and amax_chunk must have the same shape")

    shape_o = o_chunk.shape

    if lse_m_chunk.shape[-1] != 1:
        raise ValueError("lse_m_chunk must have shape (..., 1)")

    if o_chunk.ndim != lse_m_chunk.ndim:
        raise ValueError("o_chunk and lse_m_chunk must have the same number of dimensions")

    if lse_m_chunk.shape[:-1] != shape_o[:-1]:
        raise ValueError("o_chunk and lse_m_chunk must have the same shape except for the last dimension")

    if not is_initial_chunk:
        if o_chunk.shape != o.shape:
            raise ValueError("o_chunk and o must have the same shape")
        if lse_m_chunk.shape != lse_m.shape:
            raise ValueError("lse_m_chunk and lse_m must have the same shape")
        if (amax is None) != (amax_chunk is None):
            raise ValueError("amax and amax_chunk must both be None or both be not None for non-initial chunks")
        if has_amax and amax_chunk.shape != amax.shape:
            raise ValueError("amax_chunk and amax must have the same shape")

    if is_initial_chunk:
        o = o_chunk
        lse_m = lse_m_chunk
        amax = amax_chunk
    else:
        if has_amax:
            d_lse_m = lse_m - lse_m_chunk
            amax_next = torch.maximum(amax_chunk, amax)
            delta_lse = amax_chunk - amax - d_lse_m
            o = o - torch.sigmoid(delta_lse) * (o - o_chunk)
            # torch.logsumexp unconditionally promotes BF16/FP16 → FP32.
            # Cast back to the accumulator dtype to prevent FP32 leaking into
            # lse_m → delta_lse → sigmoid → o on subsequent ring steps.
            lse_m = lse_m_chunk + torch.logsumexp(
                torch.cat([(amax - amax_next) + d_lse_m, amax_chunk - amax_next], dim=-1),
                dim=-1,
                keepdim=True,
            ).to(dtype=lse_m_chunk.dtype)
            # # TODO: the following double-buffer approach can be used equivalently
            # # but save some memory
            # # First create the double buffer to store:
            # # d_lse_m[..., 0] = amax_chunk - amax + lse_m - lse_m_chunk
            # # d_lse_m[..., 1] = (o - o_chunk)
            # d_lse_m = amax_chunk.repeat_interleave(2, dim=-1)
            # amax = amax.squeeze(-1)
            # lse_m_chunk = lse_m_chunk.squeeze(-1)
            # lse_m = lse_m.squeeze(-1)
            # d_lse_m[..., 0] -= amax
            # d_lse_m[..., 0] += lse_m_chunk
            # d_lse_m[..., 0] -= lse_m
            # d_lse_m[..., 1] = o.squeeze(-1)
            # d_lse_m[..., 1] -= o_chunk.squeeze(-1)
            # # then do the sigmoid and update o
            # d_lse_m[..., 0].sigmoid_()
            # d_lse_m[..., 0] *= d_lse_m[..., 1]
            # o = o - d_lse_m[..., 0].reshape_as(o)
            # # reuse the double buffer to update lse_m
            # # amax_next = torch.maximum(amax_chunk, amax)
            # # d_lse_m[..., 0] = -amax_next
            # # d_lse_m[..., 1] = -amax_next
            # d_lse_m[..., 0] = amax_next.squeeze(-1)
            # d_lse_m[..., 0].neg_()
            # d_lse_m[..., 1] = d_lse_m[..., 0]
            # # d_lse_m[..., 0] = amax - amax_next + lse_m - lse_m_chunk
            # d_lse_m[..., 0] += amax
            # d_lse_m[..., 0] += lse_m
            # d_lse_m[..., 0] -= lse_m_chunk
            # # d_lse_m[..., 1] = amax_chunk - amax_next
            # d_lse_m[..., 1] += amax_chunk.squeeze(-1)
            # lse_m = d_lse_m.logsumexp(dim=-1, keepdim=True)
            # lse_m += lse_m_chunk.unsqueeze(-1)
            amax = amax_next
        else:
            # NOTE: without amax taking away contribution from extreme values,
            # d_lse_m can result in catastrophic cancellation. This whole branch
            # of update therefore has much smaller dynamic range as compared to
            # has_max being True. Nonetheless, in commonly encountered usage cases
            # where the "-inf" padding in the attention score only shows up after
            # those normal values but not preceding them, we would have a lse pattern of:
            # [... <stretch of normal values>, -inf, -inf, ..., -inf], which works
            # still despite -inf + normal_value = -inf because sigmoid(-inf) ~= 0
            # and logsigmoid(inf) ~= 0 so the trailing -inf would be virtually
            # discarded
            d_lse_m = lse_m - lse_m_chunk
            delta_lse = -d_lse_m
            o = o - torch.sigmoid(delta_lse) * (o - o_chunk)
            # when amax is None, lse_m is lse, i.e., we assume lse_m = lse_m + amax
            lse_m = lse_m - torch.nn.functional.logsigmoid(d_lse_m)
            amax = None

    return o, lse_m, amax


def create_and_broadcast_tensor_into_placements(
    shape: tuple[int, ...],
    create_local_fn: Callable[[Iterable[int], torch.dtype, torch.device], torch.Tensor],
    device_mesh: DeviceMesh,
    placements: tuple[Shard | Replicate, ...],
    dtype: torch.dtype = torch.float32,
) -> Tensor:
    """Create a local tensor and broadcast it independently along each replicate axis.

    With multiple replicate axes, we first create a local tensor at the source rank and broadcast
    it to all other ranks. Source rank is identified as the intersection of group rank zeros along
    all replicate axes. As such, create_local_fn is only called on a subset of ranks, i.e. the
    source rank(s).

    Parameters
    ----------
    shape : tuple[int, ...]
        Shape of the random tensor.
    create_local_fn : Callable[[Iterable[int], torch.dtype, torch.device], torch.Tensor]
        Function to create a local tensor from the shape of local tensor and dtype.
    device_mesh : DeviceMesh
        The device mesh for DTensor operations.
    placements : tuple[Shard, ...]
        The placements of the random tensor.
    dtype : torch.dtype, optional
        The dtype of the target tensor. Defaults to torch.float32.

    Returns
    -------
    Tensor
        The local tensor after broadcasting.
    """
    axis_groups = (
        device_mesh.get_all_groups()
    )  # reference: https://github.com/pytorch/pytorch/blob/v2.8.0/torch/distributed/device_mesh.py#L793

    shape_local = list(shape)
    is_mesh_dim_replicated = [placement.is_replicate() for placement in placements]

    # check placements
    for i_dim_mesh, placement in enumerate(placements):
        if placement.is_shard():
            # check if sharding is even on cp axes
            if shape[placement.dim] % device_mesh.shape[i_dim_mesh] != 0:
                raise ValueError(
                    f"Uneven sharding tensor dimension {placement.dim} of size {shape[placement.dim]} "
                    f"along device mesh dimension {i_dim_mesh} of size {device_mesh.shape[i_dim_mesh]} is not supported"
                )
            shape_local[placement.dim] = shape[placement.dim] // device_mesh.shape[i_dim_mesh]
        elif placement.is_partial():
            raise ValueError(f"Partial placements are not supported yet but got {placements}")

    is_source_rank = all(
        not is_replicated or group_rank == 0
        for is_replicated, group_rank in zip(is_mesh_dim_replicated, device_mesh.get_coordinate())
    )

    if is_source_rank:
        tensor_local = create_local_fn(shape_local=shape_local, dtype=dtype, device=device_mesh.device_type)
    else:
        tensor_local = torch.empty(shape_local, device=device_mesh.device_type, dtype=dtype)

    for axis_group, is_replicated in zip(axis_groups, is_mesh_dim_replicated):
        if not is_replicated:
            continue
        torch.distributed.broadcast(
            tensor_local,
            torch.distributed.get_global_rank(axis_group, 0),
            group=axis_group,
        )

    return tensor_local


def create_distributed_randn(
    shape: tuple[int, ...],
    device_mesh: DeviceMesh,
    placements: tuple[Shard | Replicate, ...],
    dtype: torch.dtype = torch.float32,
    scale: float = 1.0,
):
    """Create a distributed random normal distributed tensor.

    Parameters
    ----------
    shape : tuple[int, ...]
        Shape of the random tensor.
    device_mesh : DeviceMesh
        The device mesh for DTensor operations.
    placements : tuple[Shard, ...]
        The placements of the random tensor.
    dtype : torch.dtype, optional
        The dtype of the random tensor. Defaults to torch.float32.
    scale : float, optional
        Scale of the normal distribution. Defaults to 1.0.

    Returns
    -------
    DTensor
        The randn DTensor with the corresponding placements.
    """

    def create_randn_fn(shape_local, dtype, device):
        return torch.randn(shape_local, dtype=dtype, device=device) * scale

    tensor_local = create_and_broadcast_tensor_into_placements(
        shape=shape,
        create_local_fn=create_randn_fn,
        device_mesh=device_mesh,
        placements=placements,
        dtype=dtype,
    )
    # by broadcasting inside the create_and_broadcast_tensor_into_placements
    # tensor_local is guaranteed to be of same shape across ranks
    # FIXME: create_and_broadcast_tensor_into_placements should be responsible
    # for creating the DTensor output
    shape_output = list(tensor_local.shape)
    for i_dim_mesh, p in enumerate(placements):
        if isinstance(p, Shard):
            shape_output[p.dim] *= device_mesh.shape[i_dim_mesh]
    shape_output = tuple(shape_output)
    stride_output = update_exhaustive_strides(tensor_local.shape, tensor_local.stride(), shape_output)

    return DTensor.from_local(tensor_local, device_mesh, placements, shape=shape_output, stride=stride_output)


def update_exhaustive_strides(
    shape_original: Sequence[int], strides_original: Sequence[int], shape_new: Sequence[int]
) -> Sequence[int]:
    """
    Update strides to maintain the same memory layout pattern when shape changes.

    This function computes new strides that preserve the same axis ordering and memory
    layout pattern as the original exhaustive layout, but with a new shape. The resulting
    strides will create an exhaustive layout with the same dimension ordering as the
    original layout.

    An exhaustive layout is one where the mapping from multidimensional indices to flat
    indices is surjective, meaning every valid flat index corresponds to at least one
    multidimensional index. Meanwhile, a non-unique layout is not practically useful
    for our application so we further require the input shape and strides to form an
    unique layout, which implies the output layout is also unique. Overall, both the
    input and output layouts are bijective

    Parameters
    ----------
    shape_original : Sequence[int]
        The original shape of the tensor layout.
    strides_original : Sequence[int]
        The original strides of the tensor layout. Must form an exhaustive layout
        with shape_original.
    shape_new : Sequence[int]
        The new shape for which to compute corresponding strides. Must have the
        same number of dimensions as shape_original.

    Returns
    -------
    Sequence[int]
        New strides that maintain the same memory layout pattern as the original
        but are compatible with the new shape. The resulting strides will form
        an exhaustive layout with shape_new.

    Raises
    ------
    ValueError
        If the original layout (shape_original, strides_original) is not exhaustive.

    Examples
    --------
    >>> # Original layout: right-aligned (row-major) for shape (2, 3, 4)
    >>> shape_orig = (2, 3, 4)
    >>> strides_orig = (12, 4, 1)  # exhaustive right-aligned strides
    >>> shape_new = (3, 5, 2)
    >>> new_strides = update_exhaustive_strides(shape_orig, strides_orig, shape_new)
    >>> # Result: (10, 2, 1) - maintains right-aligned pattern

    Notes
    -----
    The algorithm works by:
    1. Creating a LayoutMap from the original shape and strides
    2. Verifying the original layout is exhaustive
    3. Reordering the new shape according to the original layout's stride ordering
    4. Computing exhaustive strides for the reordered new shape
    5. Reordering the computed strides back to match the original dimension order

    This is useful when reshaping tensors while preserving their memory access patterns,
    particularly in distributed computing scenarios where maintaining consistent
    memory layouts across different tensor shapes is important.
    """
    layout_original = LayoutMap(tuple(strides_original), tuple(shape_original))
    if not layout_original.is_exhaustive:
        raise ValueError(f"Input layout with shape {shape_original} and strides {strides_original} is not exhaustive")
    shape_new_ascending = np.array(shape_new)[layout_original._argsort_ascend_strides]
    argsort_output = np.argsort(layout_original._argsort_ascend_strides)
    strides_new_ascending = np.concatenate(([1], shape_new_ascending[:-1])).cumprod()
    strides_new = strides_new_ascending[argsort_output]
    return tuple(strides_new.tolist())
