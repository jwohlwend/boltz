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

"""Tests for query window ownership computation in distributed windowed attention."""

import pytest

from boltz.distributed.model.layers.utils import (
    compute_query_window_ownership,
    get_halo_from_neighbors,
)


def compute_rank_ranges(K: int, n_ranks: int, rank: int) -> tuple[tuple[int, int], tuple[int, int]]:
    """
    Test helper to compute rank ranges for simulating DTensor sharding.

    IMPORTANT: Half-window ownership is derived from query window ownership,
    not computed independently. This ensures hw_start = 2*qw_start always holds.
    """
    # Query windows (contiguous assignment)
    qw_per_rank = K // n_ranks
    qw_remainder = K % n_ranks

    if rank < qw_remainder:
        qw_start = rank * (qw_per_rank + 1)
        qw_end = qw_start + qw_per_rank + 1
    else:
        qw_start = rank * qw_per_rank + qw_remainder
        qw_end = qw_start + qw_per_rank

    # Half-windows: DERIVED from query windows (QW i owns HW [2i, 2i+1])
    hw_start = 2 * qw_start
    hw_end = 2 * qw_end

    return (qw_start, qw_end), (hw_start, hw_end)


@pytest.mark.parametrize("W,H", params := [(32, 128), (16, 64), (32, 256)], ids=[f"W:{x[0]}-H:{x[1]}" for x in params])
@pytest.mark.parametrize("K_per_rank", [1, 2, 10, 12, 32, 100], ids=lambda x: f"K_per_rank:{x}")
@pytest.mark.parametrize("n_ranks", [2, 3, 4, 8], ids=lambda x: f"n_ranks:{x}")
def test_ownership_coverage_and_validity(W, H, K_per_rank, n_ranks):
    """
    Test that ownership mapping covers all query windows exactly once and is valid.

    NOTE: this test doesn't cover the case K % n_ranks != 0, i.e.,
    the number of query windows is not a multiple of the number of ranks, which is
    not supported by the DistributedGatherSlidingWindows currently.

    Verifies:
    1. All query windows [0, K) are owned by exactly one rank
    2. No gaps or overlaps in ownership
    3. All ranges are valid (start <= end)
    4. Halos are non-negative
    """
    K = K_per_rank * n_ranks
    all_owned_qws = []

    # Assuming uniform distribution as per DistributedGatherSlidingWindows requirement K % n_ranks == 0
    # The compute_rank_ranges handles non-uniform but we restrict here to uniform to match utils.py logic
    # Actually utils.py requires K % n_ranks == 0.
    local_hw_len = (2 * K) // n_ranks

    for r in range(n_ranks):
        (qw_start, qw_end), (hw_start, hw_end) = compute_rank_ranges(K, n_ranks, r)
        ownership = compute_query_window_ownership(W, H, K, qw_start, qw_end)

        hw_need_start, hw_need_end = ownership["hw_needed"]
        left_halo = ownership["left_halo_size"]
        right_halo = ownership["right_halo_size"]

        # Valid ranges
        assert qw_start <= qw_end, f"Rank {r}: invalid qw_range"
        assert hw_start <= hw_end, f"Rank {r}: invalid hw_owned"
        assert hw_need_start <= hw_need_end, f"Rank {r}: invalid hw_needed"

        # Non-negative halos
        assert left_halo >= 0, f"Rank {r}: negative left_halo"
        assert right_halo >= 0, f"Rank {r}: negative right_halo"

        # Verify sufficiency of neighbors using get_halo_from_neighbors
        recv_meta, _ = get_halo_from_neighbors(r, n_ranks, local_hw_len, W, H, K)

        received_left = sum(length for _, htype, _, length in recv_meta if htype == "left")
        received_right = sum(length for _, htype, _, length in recv_meta if htype == "right")

        assert (
            received_left == left_halo
        ), f"Rank {r}: insufficient left halo from neighbors. Need {left_halo}, received {received_left}"
        assert (
            received_right == right_halo
        ), f"Rank {r}: insufficient right halo from neighbors. Need {right_halo}, received {received_right}"

        # Collect owned windows
        if qw_start < qw_end:
            all_owned_qws.extend(range(qw_start, qw_end))

    # Verify complete coverage
    assert sorted(all_owned_qws) == list(range(K)), f"W={W},H={H},K={K},n_ranks={n_ranks}: Coverage mismatch"


@pytest.mark.parametrize("W,H", params := [(32, 128), (16, 64), (32, 256)], ids=[f"W:{x[0]}-H:{x[1]}" for x in params])
@pytest.mark.parametrize("K_per_rank", [1, 2, 10, 12, 32, 100], ids=lambda x: f"K_per_rank:{x}")
@pytest.mark.parametrize("n_ranks", [2, 3, 4, 8], ids=lambda x: f"n_ranks:{x}")
def test_get_halo_from_neighbors_symmetry(W, H, K_per_rank, n_ranks):
    """
    Test that get_halo_from_neighbors produces symmetric recv/send patterns.
    If Rank A expects to receive from Rank B, Rank B must expect to send to Rank A.
    """
    K = K_per_rank * n_ranks
    local_hw_len = (2 * K) // n_ranks

    for r in range(n_ranks):
        recv_meta, send_meta = get_halo_from_neighbors(r, n_ranks, local_hw_len, W, H, K)

        # Check each receive item
        for peer, htype, offset_in_halo, length in recv_meta:
            # Check peer's send list
            peer_recv, peer_send = get_halo_from_neighbors(peer, n_ranks, local_hw_len, W, H, K)

            # Look for a matching send: to me (r), same length
            # Note: offset logic is trickier to verify cross-rank without full reconstruction,
            # but length and existence are critical.
            matches = [s for s in peer_send if s[0] == r and s[2] == length]
            assert matches, (
                f"Rank {r} expects to recv {length} items (type {htype}) from {peer}, "
                f"but {peer} has no matching send to {r}. Peer sends: {peer_send}"
            )

        # Check each send item
        for peer, offset_in_local, length in send_meta:
            # Check peer's recv list
            peer_recv, peer_send = get_halo_from_neighbors(peer, n_ranks, local_hw_len, W, H, K)

            matches = [r_item for r_item in peer_recv if r_item[0] == r and r_item[3] == length]
            assert matches, (
                f"Rank {r} expects to send {length} items to {peer}, "
                f"but {peer} has no matching recv from {r}. Peer recvs: {peer_recv}"
            )


@pytest.mark.parametrize("W,H", params := [(32, 128), (16, 64), (32, 256)], ids=[f"W:{x[0]}-H:{x[1]}" for x in params])
@pytest.mark.parametrize("K_per_rank", [1, 2, 10, 12, 32, 100], ids=lambda x: f"K_per_rank:{x}")
@pytest.mark.parametrize("n_ranks", [2, 3, 4, 8], ids=lambda x: f"n_ranks:{x}")
def test_halo_sufficiency(W, H, K_per_rank, n_ranks):
    """Test that halos are sufficient to cover all needed half-windows."""
    K = K_per_rank * n_ranks
    for r in range(n_ranks):
        (qw_start, qw_end), (hw_start, hw_end) = compute_rank_ranges(K, n_ranks, r)
        ownership = compute_query_window_ownership(W, H, K, qw_start, qw_end)

        hw_need_start, hw_need_end = ownership["hw_needed"]
        left_halo = ownership["left_halo_size"]
        right_halo = ownership["right_halo_size"]

        # With halos, we should cover all needed half-windows
        available_start = hw_start - left_halo
        available_end = hw_end + right_halo

        assert available_start <= hw_need_start, f"Rank {r}: left halo insufficient"
        assert available_end >= hw_need_end, f"Rank {r}: right halo insufficient"


@pytest.mark.parametrize(
    "W,H,K,n_ranks",
    [
        (32, 128, 10, 2),
        (32, 128, 12, 3),
        (32, 128, 32, 4),
    ],
)
def test_boundary_ranks_one_sided_halos(W, H, K, n_ranks):
    """Test that first and last ranks have one-sided halos."""
    # First rank
    (qw_start, qw_end), (hw_start, hw_end) = compute_rank_ranges(K, n_ranks, 0)
    ownership_0 = compute_query_window_ownership(W, H, K, qw_start, qw_end)
    assert ownership_0["left_halo_size"] == 0, "First rank should have no left halo"

    # Last rank
    (qw_start, qw_end), (hw_start, hw_end) = compute_rank_ranges(K, n_ranks, n_ranks - 1)
    ownership_last = compute_query_window_ownership(W, H, K, qw_start, qw_end)
    assert ownership_last["right_halo_size"] == 0, "Last rank should have no right halo"


def test_ownership_example_k12_n3():
    """Concrete example test for K=12, n_ranks=3."""
    W, H, K, n_ranks = 32, 128, 12, 3

    # Rank 0
    (qw_start, qw_end), (hw_start, hw_end) = compute_rank_ranges(K, n_ranks, 0)
    assert (qw_start, qw_end) == (0, 4)
    assert (hw_start, hw_end) == (0, 8)
    own0 = compute_query_window_ownership(W, H, K, qw_start, qw_end)
    assert own0["left_halo_size"] == 0  # First rank

    # Rank 1
    (qw_start, qw_end), (hw_start, hw_end) = compute_rank_ranges(K, n_ranks, 1)
    assert (qw_start, qw_end) == (4, 8)
    assert (hw_start, hw_end) == (8, 16)

    # Rank 2
    (qw_start, qw_end), (hw_start, hw_end) = compute_rank_ranges(K, n_ranks, 2)
    assert (qw_start, qw_end) == (8, 12)
    assert (hw_start, hw_end) == (16, 24)
    own2 = compute_query_window_ownership(W, H, K, qw_start, qw_end)
    assert own2["right_halo_size"] == 0  # Last rank


def test_n_ranks_greater_than_K():
    """Test edge case where we have more ranks than query windows."""
    W, H, K, n_ranks = 32, 128, 5, 10

    ranks_with_windows = 0
    ranks_without_windows = 0

    for r in range(n_ranks):
        (qw_start, qw_end), (hw_start, hw_end) = compute_rank_ranges(K, n_ranks, r)
        ownership = compute_query_window_ownership(W, H, K, qw_start, qw_end)

        if qw_start < qw_end:
            ranks_with_windows += 1
        else:
            ranks_without_windows += 1
            # Ranks without windows should have zero halos
            assert ownership["left_halo_size"] == 0
            assert ownership["right_halo_size"] == 0

    assert ranks_with_windows == K, "Should have K ranks with windows"
    assert ranks_without_windows == n_ranks - K


def get_halo_from_neighbors_iterative(
    rank: int,
    size_group: int,
    n_half_windows_local: int,
    W: int,
    H: int,
    K: int,
) -> tuple[list, list]:
    """
    Reference iterative implementation of get_halo_from_neighbors for testing purposes.
    Copied from previous implementation.
    """
    # 2. Compute ownership for ALL ranks
    all_ownerships = []
    for r in range(size_group):
        r_hw_start = r * n_half_windows_local
        r_hw_end = (r + 1) * n_half_windows_local
        r_qw_start = r_hw_start // 2
        r_qw_end = r_hw_end // 2
        if r_qw_start < r_qw_end:
            all_ownerships.append(compute_query_window_ownership(W, H, K, r_qw_start, r_qw_end))
        else:
            all_ownerships.append(None)  # Handle empty ranks if K < size_group

    my_own = all_ownerships[rank]
    # Handle case where rank has no windows
    if my_own is None:
        hw_start, hw_end = rank * n_half_windows_local, (rank + 1) * n_half_windows_local
        hw_need_start, hw_need_end = hw_start, hw_start
    else:
        hw_start, hw_end = my_own["hw_owned"]
        hw_need_start, hw_need_end = my_own["hw_needed"]

    recv_meta = []  # (peer, halo_type, offset_in_halo, length)
    send_meta = []  # (peer, offset_in_local, length)

    for peer in range(size_group):
        if peer == rank:
            continue

        peer_own = all_ownerships[peer]
        peer_hw_start = peer * n_half_windows_local
        peer_hw_end = (peer + 1) * n_half_windows_local

        # --- RECV Logic (What I need from peer) ---
        # overlap between two intervals:
        # I need: [hw_need_start, hw_start) and
        # Peer owns: [peer_hw_start, peer_hw_end)
        l_start = max(hw_need_start, peer_hw_start)
        l_end = min(hw_start, peer_hw_end)
        if l_start < l_end:
            recv_meta.append((peer, "left", l_start - hw_need_start, l_end - l_start))

        # overlap between two intervals:
        # I need: [hw_end, hw_need_end) and
        # Peer owns: [peer_hw_start, peer_hw_end)
        r_start = max(hw_end, peer_hw_start)
        r_end = min(hw_need_end, peer_hw_end)
        if r_start < r_end:
            recv_meta.append((peer, "right", r_start - hw_end, r_end - r_start))

        # --- SEND Logic (What peer needs from me) ---
        if peer_own is None:
            continue
        p_need_start, p_need_end = peer_own["hw_needed"]
        p_hw_start, p_hw_end = peer_own["hw_owned"]

        # overlap between two intervals:
        # Peer needs: [p_need_start, p_hw_start) and
        # I own: [hw_start, hw_end)
        l_start = max(p_need_start, hw_start)
        l_end = min(p_hw_start, hw_end)
        if l_start < l_end:
            send_meta.append((peer, l_start - hw_start, l_end - l_start))

        # overlap between two intervals:
        # Peer needs: [p_hw_end, p_need_end) and
        # I own: [hw_start, hw_end)
        r_start = max(p_hw_end, hw_start)
        r_end = min(p_need_end, hw_end)
        if r_start < r_end:
            send_meta.append((peer, r_start - hw_start, r_end - r_start))

    return recv_meta, send_meta


@pytest.mark.parametrize("W,H", params := [(32, 128), (16, 64), (32, 256)], ids=[f"W:{x[0]}-H:{x[1]}" for x in params])
@pytest.mark.parametrize("K_per_rank", [1, 2, 10, 12, 32, 100], ids=lambda x: f"K_per_rank:{x}")
@pytest.mark.parametrize("n_ranks", [2, 3, 4, 8], ids=lambda x: f"n_ranks:{x}")
def test_get_halo_vectorized_vs_iterative(W, H, K_per_rank, n_ranks):
    """
    Compare new vectorized get_halo_from_neighbors against the old iterative implementation.
    """
    K = K_per_rank * n_ranks

    local_hw_len = (2 * K) // n_ranks

    for r in range(n_ranks):
        recv_vec, send_vec = get_halo_from_neighbors(r, n_ranks, local_hw_len, W, H, K)
        recv_iter, send_iter = get_halo_from_neighbors_iterative(r, n_ranks, local_hw_len, W, H, K)

        # Sort to ensure order consistency for comparison
        # recv tuple: (peer, type, offset, length)
        recv_vec.sort()
        recv_iter.sort()

        # send tuple: (peer, offset, length)
        send_vec.sort()
        send_iter.sort()

        assert recv_vec == recv_iter, f"Rank {r}: Recv mismatch.\nVectorized: {recv_vec}\nIterative: {recv_iter}"
        assert send_vec == send_iter, f"Rank {r}: Send mismatch.\nVectorized: {send_vec}\nIterative: {send_iter}"
