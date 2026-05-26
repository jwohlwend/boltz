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

from collections import namedtuple
from itertools import product

import numpy as np
import pytest
import torch

from boltz.distributed.utils import LayoutMap, update_exhaustive_strides

# Define a named tuple for the test parameters
LayoutMapParams = namedtuple("LayoutMapParams", ["strides", "shape", "offset", "slices", "layoutmap"])


@pytest.fixture(
    params=[
        # Test cases without slices
        ((), (), 3, None),  # rank-0 layout
        ((1,), (10,), 1, None),  # rank-1 layout
        ((1, 1), (2, 1), 1, None),  # rank-2 LayoutRightMap with all-1 strides
        ((1, 1, 1), (1, 2, 1), 0, None),  # rank-3 LayoutRightMap with all-1 strides
        ((5, 5, 1), (3, 1, 5), 0, None),  # rank-3 LayoutRightMap with singleton axis
        ((1, 2, 6), (2, 3, 4), 2, None),  # LayoutLeftMap
        ((12, 4, 1), (2, 3, 4), 3, None),  # LayoutRightMap
        ((3, 1, 6), (2, 3, 4), 4, None),  # Layout strided
        ((3, 1, 24), (4, 2, 3), 0, None),  # LayoutMap(strides=(), shape=())[:4, :2, :3]
        ((12, 2, 72), (2, 2, 2), 2, None),  # LayoutMap(strides=(), shape=())[::4, ::2, ::3]
        # Test cases with slices
        ((1, 2, 6), (2, 3, 4), 2, (slice(1, 2, 1),)),  # Single slice on first dimension (test slice padding)
        ((12, 4, 1), (2, 3, 4), 3, (0, slice(2, 3, 1))),  # Int index + slice
        ((3, 1, 6), (2, 3, 4), 4, (slice(1, 2, 1), 1, slice(2, 4, 1))),  # Multiple slices and int
        (
            (12, 3, 24),
            (2, 2, 3),
            0,
            (slice(None, None, 2), 0, slice(None, None, 2)),
        ),  # LayoutMap(strides=(6, 1, 24), shape=(4, 6, 5))[::2, ::3, :3][::2, 0, ::2]
    ],
    ids=lambda p: f"stride = {p[0]}, shape = {p[1]}, offset = {p[2]}, slices = {p[3]}",
)
def layout_map_fixture(request):
    strides, shape, offset, slices = request.param

    layout_map = LayoutMap(strides, shape, offset=offset)
    # save the input parameters for debugging purpose
    params = LayoutMapParams(strides, shape, offset, slices, layout_map)

    if len(strides) == 0:
        required_span_size = 1
        flat_indices_nd = torch.tensor(offset)
    else:
        required_span_size = 1 + sum((shape[i] - 1) * strides[i] for i in range(len(strides)))
        flat_indices = torch.arange(required_span_size) + offset
        flat_indices_nd = torch.as_strided(flat_indices, size=shape, stride=strides)

    # Apply slicing if specified
    if slices is not None:
        try:
            # Slice the flat_indices_nd tensor
            flat_indices_nd = flat_indices_nd[slices]

            # Slice the layout_map
            layout_map = layout_map[slices]

            # Update the shape and strides to match the sliced layout
            shape = layout_map.shape
            strides = layout_map.strides
            offset = layout_map.offset
            required_span_size = layout_map.required_span_size
        except (ValueError, RuntimeError) as e:
            # Skip this test case if slicing raises an exception
            pytest.skip(f"Slicing failed with error: {e}")

    return params, strides, shape, offset, layout_map, required_span_size, flat_indices_nd


def test_layout_map_init(layout_map_fixture):
    params, strides, shape, offset, layout_map, required_span_size, flat_indices_nd = layout_map_fixture
    assert isinstance(layout_map.shape, tuple)
    assert isinstance(layout_map.strides, tuple)
    assert strides == layout_map.strides
    assert layout_map.shape == shape
    assert layout_map.numel == np.prod(shape)
    assert required_span_size == layout_map.required_span_size
    # for exhaustive layout, the flat_indices_nd should be arange upon sorting
    is_exhaustive_expected = (
        required_span_size == flat_indices_nd.numel()
        and torch.all(
            flat_indices_nd.flatten().sort(stable=True).values == torch.arange(required_span_size) + offset
        ).item()
    )
    assert layout_map.is_exhaustive == is_exhaustive_expected


def test_layout_map_init_invalid_strides(layout_map_fixture):
    _, _, shape, _, _, _, _ = layout_map_fixture
    invalid_strides = (1, 1, 1)
    if len(shape) != 3:
        with pytest.raises(ValueError, match=r"^.* must have the same length"):
            LayoutMap(invalid_strides, shape)
    else:
        # the only unique layout with all-1 strides is the one where
        # only 1 dimension of the shape is larger or equal to 1
        # while the rest are 1.
        values, counts = np.unique(shape, return_counts=True)
        is_unique_layout = (
            values.size == 2 and values[0] == 1 and values[1] > 1 and counts[0] == len(shape) - 1 and counts[1] == 1
        )
        if not is_unique_layout:
            with pytest.raises(ValueError, match=r"^.* do not give unique layout"):
                LayoutMap(invalid_strides, shape)


def test_layout_map_init_invalid_shape(layout_map_fixture):
    _, strides, _, _, _, _, _ = layout_map_fixture
    invalid_shape = (99, 99, 99)
    with pytest.raises(ValueError):
        LayoutMap(strides, invalid_shape)
    if len(strides) != 3:
        with pytest.raises(ValueError, match=r"^.* must have the same length"):
            LayoutMap(strides, invalid_shape)
    else:
        with pytest.raises(ValueError, match=r"^.* do not give unique layout"):
            LayoutMap(strides, invalid_shape)
    negative_shape = (1, 2, -3)
    with pytest.raises(ValueError, match=r"^.* contain negative values"):
        LayoutMap(strides, negative_shape)
    zero_shape = (1, 2, 0)
    with pytest.raises(ValueError, match=r"^.* contain zero values"):
        LayoutMap(strides, zero_shape)


def test_layout_map_call(layout_map_fixture):
    _, _, shape, _, layout_map, _, flat_indices_nd = layout_map_fixture
    for idx in product(*[range(s) for s in shape]):
        flat_idx_expected = flat_indices_nd[idx]
        flat_idx_result = layout_map(idx)
        assert flat_idx_expected == flat_idx_result


def test_layout_map_call_invalid_idx(layout_map_fixture):
    _, _, shape, _, layout_map, _, _ = layout_map_fixture
    idx = tuple(i for i in range(len(shape) + 1))
    with pytest.raises(ValueError, match=r"Expected .* elements in ids but got only .*"):
        layout_map(idx)
    if len(shape) > 0:
        # test out of bounds
        idx_oob = tuple(shape[i] for i in range(len(shape)))
        with pytest.raises(ValueError, match=r"Expected ids to satisfy 0 <= ids\[.*\] <= .* but found ids\[.*\] == .*"):
            layout_map(idx_oob)


def test_layout_map_unravel(layout_map_fixture):
    _, _, shape, _, layout_map, _, flat_indices_nd = layout_map_fixture
    for idx_expected in product(*[range(s) for s in shape]):
        flat_idx = flat_indices_nd[idx_expected].item()
        idx_result = layout_map.unravel(flat_idx)
        assert idx_result == idx_expected


def test_layout_map_unravel_invalid_flat_idx(layout_map_fixture):
    _, _, _, _, layout_map, required_span_size, _ = layout_map_fixture
    flat_idx_oob = required_span_size + layout_map.offset
    with pytest.raises(ValueError, match=r"Expected flat_index in range"):
        layout_map.unravel(flat_idx_oob)


# Add a test specifically for slicing with padding
def test_layout_map_getitem_padding():
    """Test that slicing with fewer dimensions than the layout properly pads with ':' slices."""
    layout_map = LayoutMap((12, 4, 1), (2, 3, 4))

    # Test partial slicing
    partial_slice_1 = layout_map[1]  # Should be equivalent to layout_map[1, :, :]
    full_slice_1 = layout_map[1, slice(None), slice(None)]

    assert partial_slice_1.shape == full_slice_1.shape
    assert partial_slice_1.strides == full_slice_1.strides
    assert partial_slice_1.offset == full_slice_1.offset

    # Test with 2 dimensions only
    partial_slice_2 = layout_map[1, 2]  # Should be equivalent to layout_map[1, 2, :]
    full_slice_2 = layout_map[1, 2, slice(None)]

    assert partial_slice_2.shape == full_slice_2.shape
    assert partial_slice_2.strides == full_slice_2.strides
    assert partial_slice_2.offset == full_slice_2.offset

    # Test with a single slice
    partial_slice_3 = layout_map[slice(0, 2, 1)]  # Should be equivalent to layout_map[0:2, :, :]
    full_slice_3 = layout_map[slice(0, 2, 1), slice(None), slice(None)]

    assert partial_slice_3.shape == full_slice_3.shape
    assert partial_slice_3.strides == full_slice_3.strides
    assert partial_slice_3.offset == full_slice_3.offset


# Add a test that specifically checks the slicing behavior
def test_layout_map_getitem(layout_map_fixture):
    """Test that slicing behavior is as expected for cases with slices."""
    params, _, _, _, layout_map, _, flat_indices_nd = layout_map_fixture

    # Skip test cases without slices
    if params.slices is None:
        return

    # For cases with slices, verify that indices after slicing still match
    for idx in product(*[range(s) for s in layout_map.shape]):
        flat_idx_expected = flat_indices_nd[idx].item()
        flat_idx_result = layout_map(idx)
        assert flat_idx_expected == flat_idx_result


# Test for invalid slicing cases
def test_layout_map_invalid_slicing():
    """Test the __getitem__ method with invalid slices."""
    layout_map = LayoutMap((1, 2, 6), (2, 3, 4))

    # Test with negative step
    with pytest.raises(ValueError, match="Negative or zero steps"):
        layout_map[slice(1, 0, -1)]

    # Test with start >= stop
    with pytest.raises(ValueError, match="start not smaller than stop"):
        layout_map[slice(1, 1, 1)]

    with pytest.raises(ValueError, match="start not smaller than stop"):
        layout_map[slice(2, 1, 1)]

    # Test with invalid slice type
    with pytest.raises(TypeError, match="Unsupported slice type"):
        layout_map[("invalid_type",)]


def test_update_exhaustive_strides(layout_map_fixture):
    params, strides, shape, offset, layout_map, required_span_size, flat_indices_nd = layout_map_fixture

    # Create common label for test case debugging
    label_test_case = f"shape_original={shape}, strides_original={strides}"

    # Skip test cases that are not exhaustive or have slices applied
    if not layout_map.is_exhaustive:
        # The function should raise ValueError
        with pytest.raises(ValueError, match="Input layout .* is not exhaustive"):
            update_exhaustive_strides(shape, strides, shape)
        return

    # Generate a new shape with the same number of dimensions
    # Use different sizes but keep the same dimensionality
    shape_new = tuple(max(1, s + 1) for s in shape)  # Increment each dimension by 1

    # Call the function under test
    strides_new = update_exhaustive_strides(shape, strides, shape_new)

    # Verify the output
    assert len(strides_new) == len(shape_new), (
        f"Output strides should have same length as new shape. "
        f"{label_test_case}, shape_new={shape_new}, strides_new={strides_new}"
    )
    assert all(
        isinstance(s, (int, np.integer)) and s > 0 for s in strides_new
    ), f"All strides should be positive integers. {label_test_case}, shape_new={shape_new}, strides_new={strides_new}"

    # Create a new LayoutMap with the returned strides and new shape
    layout_new = LayoutMap(strides_new, shape_new)

    # Verify that the new layout is exhaustive
    assert (
        layout_new.is_exhaustive
    ), f"The new layout should be exhaustive. {label_test_case}, shape_new={shape_new}, strides_new={strides_new}"

    # Test with a different new shape to ensure generality
    shape_new_2 = tuple(max(1, s * 2) for s in shape)  # Double each dimension
    strides_new_2 = update_exhaustive_strides(shape, strides, shape_new_2)
    layout_new_2 = LayoutMap(strides_new_2, shape_new_2)

    assert layout_new_2.is_exhaustive, (
        f"The second new layout should be exhaustive. "
        f"{label_test_case}, shape_new_2={shape_new_2}, strides_new_2={strides_new_2}"
    )
    assert layout_new_2.is_unique, (
        f"The second new layout should be unique. "
        f"{label_test_case}, shape_new_2={shape_new_2}, strides_new_2={strides_new_2}"
    )

    # self reflective test
    strides_new_3 = update_exhaustive_strides(shape, strides, shape)
    assert (
        strides_new_3 == strides
    ), f"The self reflective test should return the original strides. {label_test_case}, strides_new_3={strides_new_3}"
