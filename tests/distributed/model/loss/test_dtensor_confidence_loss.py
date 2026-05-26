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


"""Tests for DTensor confidence_loss wrapper implementation.

This module tests the distributed implementation of confidence_loss which aggregates
plddt_loss, pde_loss, and resolved_loss. The tests verify numerical correctness
against the serial implementation and proper gradient computation.

The confidence_loss wrapper coordinates DTensor placements across sub-loss functions
and returns aggregated scalar losses.
"""

from math import gcd

import pytest
import torch
from torch.distributed.tensor import Replicate, Shard, distribute_tensor

from boltz.data.feature.featurizer import BoltzFeaturizer
from boltz.data.module.inference import load_input
from boltz.data.tokenize.boltz import BoltzTokenizer
from boltz.distributed.comm import TransposeComm
from boltz.distributed.data.feature.featurizer_utils import get_num_atoms_tokens
from boltz.distributed.data.utils import distribute_features
from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.loss.confidencev2 import (
    compute_frame_pred as dtensor_compute_frame_pred,
)
from boltz.distributed.model.loss.confidencev2 import confidence_loss
from boltz.model.layers.confidence_utils import compute_frame_pred as serial_compute_frame_pred
from boltz.model.loss.confidencev2 import confidence_loss as serial_confidence_loss
from boltz.testing.utils import (
    distribute_atom_features,
    random_features,
    spawn_multiprocessing,
)


def _assert_nontrivial_expected_mask_collinear(
    expected_mask_collinear_host: torch.Tensor,
    token_pad_mask_host: torch.Tensor,
    test_name: str,
) -> None:
    """Ensure expected mask_collinear has valid support and non-trivial positives."""
    token_pad_mask_bool = token_pad_mask_host.bool()
    if expected_mask_collinear_host.numel() == 0:
        raise AssertionError(f"{test_name}: expected_mask_collinear_host is empty")
    if not token_pad_mask_bool.any():
        raise AssertionError(f"{test_name}: token_pad_mask has no valid tokens")

    for batch_idx in range(expected_mask_collinear_host.shape[0]):
        valid = expected_mask_collinear_host[batch_idx, :, token_pad_mask_bool[batch_idx]]
        if valid.numel() == 0:
            raise AssertionError(f"{test_name}: batch {batch_idx} has no valid tokens for expected_mask_collinear")
        if not valid.any():
            raise AssertionError(
                f"{test_name}: batch {batch_idx} expected_mask_collinear on valid tokens is always False"
            )


def parallel_assert_compute_frame_pred(rank: int, payload: tuple) -> None:
    (
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        multiplicity,
        pred_atom_coords_host,
        feats_host,
        atom_counts_per_token_host,
        expected_frames_idx_host,
        expected_mask_collinear_host,
    ) = payload

    monkeypatch = pytest.MonkeyPatch()
    if env_per_rank is not None:
        for var_name, value in env_per_rank.items():
            if value == "<INPUT_RANK>":
                monkeypatch.setenv(var_name, f"{rank}")
                continue
            monkeypatch.setenv(var_name, value)

    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    manager = DistributedManager()
    device_mesh = manager.device_mesh_subgroups

    single_repr_placements = (Shard(0), Shard(1), Replicate())
    replicate_placements = (Shard(0), Replicate(), Replicate())

    size_batch = feats_host["atom_pad_mask"].shape[0]
    pred_coords_unflat = pred_atom_coords_host.unflatten(0, (size_batch, multiplicity))
    inputs_atom = {
        "atom_counts_per_token": atom_counts_per_token_host.to(dtype=torch.int64),
        "pred_atom_coords_0": pred_coords_unflat[:, 0].to(dtype=pred_atom_coords_host.dtype),
        "atom_to_token": feats_host["atom_to_token"].to(dtype=pred_atom_coords_host.dtype),
        "atom_pad_mask": feats_host["atom_pad_mask"].to(dtype=pred_atom_coords_host.dtype),
        "atom_resolved_mask": feats_host["atom_resolved_mask"].to(dtype=pred_atom_coords_host.dtype),
        "frames_idx": feats_host["frames_idx"].to(dtype=torch.int64),
    }
    for i_mul in range(1, multiplicity):
        inputs_atom[f"pred_atom_coords_{i_mul}"] = pred_coords_unflat[:, i_mul].to(dtype=pred_atom_coords_host.dtype)

    placements_cp = {
        "atom_counts_per_token": (Shard(0), Replicate()),
        "atom_to_token": (Shard(0), Replicate()),
        "atom_pad_mask": (Shard(0), Replicate()),
        "atom_resolved_mask": (Shard(0), Replicate()),
        "frames_idx": (Shard(1), Replicate()),
        "pred_atom_coords_0": (Shard(0), Replicate()),
    }
    placements_dp_cp = {
        "atom_to_token": (Shard(0), Shard(1), Replicate()),
        "atom_pad_mask": (Shard(0), Shard(1), Replicate()),
        "atom_resolved_mask": (Shard(0), Shard(1), Replicate()),
        "frames_idx": (Shard(0), Shard(1), Replicate()),
        "pred_atom_coords_0": (Shard(0), Shard(1), Replicate()),
    }
    for i_mul in range(1, multiplicity):
        placements_cp[f"pred_atom_coords_{i_mul}"] = (Shard(0), Replicate())
        placements_dp_cp[f"pred_atom_coords_{i_mul}"] = (Shard(0), Shard(1), Replicate())

    feats_atom = distribute_atom_features(
        inputs_atom,
        placements_cp,
        placements_dp_cp,
        device_mesh,
        manager.group["cp"],
        multiplicities={"pred_atom_coords": multiplicity},
    )
    pred_atom_coords = feats_atom["pred_atom_coords"]
    frames_idx_true = feats_atom["frames_idx"]

    feats = {
        "asym_id": distribute_tensor(
            feats_host["asym_id"].to(manager.device), device_mesh=device_mesh, placements=single_repr_placements
        ),
        "atom_to_token": feats_atom["atom_to_token"],
        "atom_pad_mask": feats_atom["atom_pad_mask"],
        "atom_resolved_mask": feats_atom["atom_resolved_mask"],
        "mol_type": distribute_tensor(
            feats_host["mol_type"].to(manager.device), device_mesh=device_mesh, placements=single_repr_placements
        ),
        "token_pad_mask": distribute_tensor(
            feats_host["token_pad_mask"].to(manager.device), device_mesh=device_mesh, placements=single_repr_placements
        ),
    }

    frames_idx_pred, mask_collinear_pred = dtensor_compute_frame_pred(
        pred_atom_coords,
        frames_idx_true,
        feats,
        multiplicity=multiplicity,
    )

    note = ""
    try:
        frames_idx_pred_local = frames_idx_pred.redistribute(device_mesh, placements=replicate_placements).to_local()
        dp_rank = manager.group_rank["dp"]
        local_batch_size = size_batch // manager.group["dp"].size()
        dp_idx_str = dp_rank * local_batch_size
        dp_idx_end = dp_idx_str + local_batch_size

        # Proxy comparison for frame indices: compare geometry implied by indices.
        # This avoids coupling test correctness to a specific index-space convention.
        pred_atom_coords_local = (
            pred_atom_coords.redistribute(device_mesh, placements=replicate_placements)
            .to_local()
            .unflatten(0, (local_batch_size, multiplicity))
        )
        expected_atom_coords_local = pred_atom_coords_host.to(manager.device).unflatten(0, (size_batch, multiplicity))[
            dp_idx_str:dp_idx_end
        ]
        expected_frames_idx_local = expected_frames_idx_host.to(manager.device)[dp_idx_str:dp_idx_end]

        batch_idx = torch.arange(local_batch_size, device=manager.device)[:, None, None, None]
        mult_idx = torch.arange(multiplicity, device=manager.device)[None, :, None, None]

        dt_frames = pred_atom_coords_local[batch_idx, mult_idx, frames_idx_pred_local]
        host_frames = expected_atom_coords_local[batch_idx, mult_idx, expected_frames_idx_local]
        # Use mean frame center as an index-invariant proxy.
        dt_frame_centers = dt_frames.mean(dim=-2)
        host_frame_centers = host_frames.mean(dim=-2)

        token_pad_mask_non_dtensor = feats_host["token_pad_mask"][dp_idx_str:dp_idx_end].bool().to(manager.device)
        token_pad_mask_dtensor = (
            feats["token_pad_mask"].redistribute(device_mesh, placements=replicate_placements).to_local().bool()
        )
        for batch_i in range(local_batch_size):
            non_mask = token_pad_mask_non_dtensor[batch_i]
            dt_mask = token_pad_mask_dtensor[batch_i]

            # Ensure test isn't trivially passing
            if non_mask.sum().item() != dt_mask.sum().item():
                raise AssertionError(
                    "frames_idx proxy token count mismatch: "
                    f"non-dtensor={non_mask.sum().item()}, dtensor={dt_mask.sum().item()}"
                )
            if not non_mask.any():
                raise AssertionError(f"batch {batch_i} has no valid tokens")

            # Compare collected frame center coordinates
            torch.testing.assert_close(
                dt_frame_centers[batch_i, :, dt_mask].cpu(),
                host_frame_centers[batch_i, :, non_mask].cpu(),
            )

    except AssertionError as e:
        note += "Test failed when comparing frames_idx_pred: " + str(e) + "\n"

    try:
        mask_collinear_full = mask_collinear_pred.full_tensor().cpu()
        token_pad_mask_non_dtensor = feats_host["token_pad_mask"].bool().cpu()
        token_pad_mask_dtensor = feats["token_pad_mask"].full_tensor().bool().cpu()
        expected_mask_collinear_host_cpu = expected_mask_collinear_host.cpu()

        for batch_idx in range(mask_collinear_full.shape[0]):
            non_mask = token_pad_mask_non_dtensor[batch_idx]
            dt_mask = token_pad_mask_dtensor[batch_idx]
            mask_collinear_valid = mask_collinear_full[batch_idx, :, dt_mask]
            expected_mask_collinear_valid = expected_mask_collinear_host_cpu[batch_idx, :, non_mask]

            # Ensure test isn't trivially passing
            if non_mask.sum().item() != dt_mask.sum().item():
                raise AssertionError(
                    "mask_collinear token count mismatch: "
                    f"non-dtensor={non_mask.sum().item()}, dtensor={dt_mask.sum().item()}"
                )
            if not mask_collinear_valid.any() or not expected_mask_collinear_valid.any():
                raise AssertionError(
                    "test can trivially pass on mask_collinear_pred since mask_collinear_pred on valid tokens is always False"
                )

            # Compare mask_collinear
            torch.testing.assert_close(
                mask_collinear_valid,
                expected_mask_collinear_valid,
            )
    except AssertionError as e:
        note += "Test failed when comparing mask_collinear_pred: " + str(e) + "\n"

    if note:
        raise AssertionError(note)

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env",
    [
        ((2, (2, 2)), True, "cuda", "ENV"),
        ((2, (2, 2)), True, "cpu", "ENV"),
    ],
    indirect=("setup_env",),
    ids=lambda x: f"dp={x[0][0]}, cp={x[0][1]}, specify_method={x[1]}, device_type={x[2]}, method_init={x[3]}",
)
@pytest.mark.parametrize("multiplicity", [1, 2], ids=["multiplicity=1", "multiplicity=2"])
@pytest.mark.parametrize("seed", [0, 42], ids=["seed=0", "seed=42"])
def test_dtensor_compute_frame_pred(setup_env: tuple, multiplicity: int, seed: int) -> None:
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    size_ring = grid_group_sizes["cp"][0]
    batch_size = grid_group_sizes["dp"]
    n_tokens_per_rank = 20
    n_tokens = size_ring * n_tokens_per_rank
    max_atoms_per_token = 18
    n_atoms_per_rank = n_tokens_per_rank * max_atoms_per_token
    n_atoms = size_ring * n_atoms_per_rank

    rng = torch.Generator(device=device_type)
    rng.manual_seed(seed)

    pred_atom_coords = torch.randn(
        (batch_size * multiplicity, n_atoms, 3), device=device_type, generator=rng, requires_grad=True
    )
    # enforce collinearity by setting some atoms to zero
    collinear_mask = torch.randint(
        0, 2, (batch_size * multiplicity, n_atoms), device=device_type, generator=rng, dtype=torch.bool
    )
    pred_atom_coords = torch.where(collinear_mask.unsqueeze(-1), pred_atom_coords, torch.zeros_like(pred_atom_coords))

    rng_features = torch.Generator(device=pred_atom_coords.device)
    rng_features.manual_seed(seed)
    feats = random_features(
        size_batch=batch_size,
        n_tokens=n_tokens,
        n_atoms=n_atoms,
        n_msa=1,
        atom_counts_per_token_range=(1, max_atoms_per_token),
        device=pred_atom_coords.device,
        float_value_range=(-1.0, 1.0),
        selected_keys=[
            "asym_id",
            "atom_to_token",
            "atom_pad_mask",
            "atom_resolved_mask",
            "frames_idx",
            "atom_counts_per_token",
            "mol_type",
            "token_pad_mask",
        ],
        rng=rng_features,
    )
    atom_pad_mask_bool = feats["atom_pad_mask"].bool()
    frames_idx = feats["frames_idx"]
    atom_pad_per_token = atom_pad_mask_bool.unsqueeze(1).expand(-1, n_tokens, -1)
    frame_atom_valid = torch.gather(atom_pad_per_token, dim=2, index=frames_idx)
    assert torch.all(frame_atom_valid), "random_features generated frames_idx pointing to masked atoms"
    feats_serial = {
        "asym_id": feats["asym_id"],
        "atom_to_token": feats["atom_to_token"],
        "atom_pad_mask": feats["atom_pad_mask"],
        "atom_resolved_mask": feats["atom_resolved_mask"],
        "mol_type": feats["mol_type"],
        "token_pad_mask": feats["token_pad_mask"],
    }
    expected_frames_idx_host, expected_mask_collinear_host = serial_compute_frame_pred(
        pred_atom_coords,
        feats["frames_idx"],
        feats_serial,
        multiplicity,
    )
    _assert_nontrivial_expected_mask_collinear(
        expected_mask_collinear_host.detach().cpu(),
        feats["token_pad_mask"].detach().cpu(),
        "test_dtensor_compute_frame_pred",
    )

    payload = (
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        multiplicity,
        pred_atom_coords.detach().clone().cpu(),
        {k: v.detach().clone().cpu() for k, v in feats.items()},
        feats["atom_counts_per_token"].detach().clone().cpu(),
        expected_frames_idx_host.detach().clone().cpu(),
        expected_mask_collinear_host.detach().clone().cpu(),
    )

    spawn_multiprocessing(
        parallel_assert_compute_frame_pred,
        world_size,
        payload,
    )


@pytest.mark.parametrize(
    "setup_env",
    [
        ((1, (2, 2)), True, "cuda", "ENV"),
        ((1, (2, 2)), True, "cpu", "ENV"),
    ],
    indirect=("setup_env",),
    ids=lambda x: f"dp={x[0][0]}, cp={x[0][1]}, specify_method={x[1]}, device_type={x[2]}, method_init={x[3]}",
)
@pytest.mark.parametrize("multiplicity", [1, 2], ids=["multiplicity=1", "multiplicity=2"])
@pytest.mark.parametrize("seed", [0, 42], ids=["seed=0", "seed=42"])
def test_dtensor_compute_frame_pred_real_data_parallel(
    setup_env: tuple,
    multiplicity: int,
    seed: int,
    create_preprocessed_handle_boltz1_v1,
) -> None:
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    processed = create_preprocessed_handle_boltz1_v1
    record = processed.manifest.records[0]
    input_data = load_input(record, processed.targets_dir, processed.msa_dir)
    tokenized = BoltzTokenizer().tokenize(input_data)
    n_atoms_raw, n_tokens_raw = get_num_atoms_tokens(tokenized)

    ring = grid_group_sizes["cp"][0]
    atoms_per_window = 32
    atom_lcm = ring * atoms_per_window // gcd(ring, atoms_per_window)
    max_atoms = ((n_atoms_raw + atom_lcm - 1) // atom_lcm) * atom_lcm
    max_tokens = ((n_tokens_raw + ring - 1) // ring) * ring
    max_seqs = ring

    feats_single = BoltzFeaturizer().process(
        tokenized,
        training=False,
        max_atoms=max_atoms,
        max_tokens=max_tokens,
        max_seqs=max_seqs,
        pad_to_max_seqs=True,
    )
    if not isinstance(feats_single, dict):
        raise TypeError("Expected non-sharded feature dict from BoltzFeaturizer.process")

    selected_keys = [
        "asym_id",
        "atom_to_token",
        "atom_pad_mask",
        "atom_resolved_mask",
        "frames_idx",
        "mol_type",
        "token_pad_mask",
    ]
    feats_single = {k: feats_single[k] for k in selected_keys}
    # v1 featurizer doesn't emit atom_counts_per_token; derive from one-hot atom_to_token
    feats_single["atom_counts_per_token"] = feats_single["atom_to_token"].sum(dim=0).to(torch.int64)

    batch_size = grid_group_sizes["dp"]
    feats = {k: v.unsqueeze(0).repeat_interleave(batch_size, dim=0).to(device_type) for k, v in feats_single.items()}

    n_atoms = feats["atom_pad_mask"].shape[1]
    rng = torch.Generator(device=device_type)
    rng.manual_seed(seed)
    pred_atom_coords = torch.randn((batch_size * multiplicity, n_atoms, 3), device=device_type, generator=rng)
    collinear_mask = torch.randint(
        0,
        2,
        (batch_size * multiplicity, n_atoms),
        device=device_type,
        generator=rng,
        dtype=torch.bool,
    )
    pred_atom_coords = torch.where(collinear_mask.unsqueeze(-1), pred_atom_coords, torch.zeros_like(pred_atom_coords))
    feats_serial = {
        "asym_id": feats["asym_id"],
        "atom_to_token": feats["atom_to_token"],
        "atom_pad_mask": feats["atom_pad_mask"],
        "atom_resolved_mask": feats["atom_resolved_mask"],
        "mol_type": feats["mol_type"],
        "token_pad_mask": feats["token_pad_mask"],
    }
    expected_frames_idx_host, expected_mask_collinear_host = serial_compute_frame_pred(
        pred_atom_coords,
        feats["frames_idx"],
        feats_serial,
        multiplicity,
    )
    _assert_nontrivial_expected_mask_collinear(
        expected_mask_collinear_host.detach().cpu(),
        feats["token_pad_mask"].detach().cpu(),
        "test_dtensor_compute_frame_pred_real_data_parallel",
    )

    payload = (
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        multiplicity,
        pred_atom_coords.detach().clone().cpu(),
        {k: v.detach().clone().cpu() for k, v in feats.items()},
        feats["atom_counts_per_token"].detach().clone().cpu(),
        expected_frames_idx_host.detach().clone().cpu(),
        expected_mask_collinear_host.detach().clone().cpu(),
    )

    spawn_multiprocessing(
        parallel_assert_compute_frame_pred,
        world_size,
        payload,
    )


def parallel_assert_confidence_loss(
    rank: int,
    payload: tuple,
):
    """Worker function that runs on each rank to test confidence_loss DTensor implementation.

    This function:
    1. Initializes the distributed environment
    2. Distributes atom and token features using appropriate utilities
    3. Calls the DTensor confidence_loss function
    4. Verifies forward pass matches serial reference
    5. Verifies gradients flow correctly through all logit tensors
    """
    (
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        multiplicity,
        alpha_pae,
        pred_lddt_global_host,
        pred_pde_global_host,
        pred_pae_global_host,
        pred_resolved_global_host,
        pred_atom_coords_global_host,
        true_atom_coords_global_host,
        true_coords_resolved_mask_global_host,
        token_to_rep_atom_global_host,
        r_set_to_rep_atom_global_host,
        atom_to_token_global_host,
        mol_type_global_host,
        token_pad_mask_global_host,
        atom_counts_per_token_host,
        pae_feats_host,
        expected_loss_host,
        expected_loss_breakdown_host,
        expected_grad_pred_lddt_host,
        expected_grad_pred_pde_host,
        expected_grad_pred_pae_host,
        expected_grad_pred_resolved_host,
    ) = payload

    # Setup environment variables for this rank
    monkeypatch = pytest.MonkeyPatch()
    if env_per_rank is not None:
        for var_name, value in env_per_rank.items():
            if value == "<INPUT_RANK>":
                monkeypatch.setenv(var_name, f"{rank}")
                continue
            monkeypatch.setenv(var_name, value)

    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    manager = DistributedManager()

    device_mesh = manager.device_mesh_subgroups
    dtype = pred_atom_coords_global_host.dtype

    # Create TransposeComm for redistribute_transpose
    comm = TransposeComm(manager.group["cp"], manager.layout_subgroups["cp"])

    # --- Distribute atom features using distribute_atom_features utility ---
    size_batch = token_to_rep_atom_global_host.shape[0]
    inputs_atom = {
        "atom_counts_per_token": atom_counts_per_token_host.to(dtype=torch.int64),
        "token_to_rep_atom": token_to_rep_atom_global_host.to(dtype=dtype),
        "r_set_to_rep_atom": r_set_to_rep_atom_global_host.to(dtype=dtype),
        "atom_to_token": atom_to_token_global_host.to(dtype=dtype),
    }
    if pae_feats_host is not None:
        inputs_atom["frames_idx"] = pae_feats_host["frames_idx"]
        inputs_atom["atom_pad_mask"] = pae_feats_host["atom_pad_mask"].to(dtype=dtype)
        inputs_atom["atom_resolved_mask"] = pae_feats_host["atom_resolved_mask"].to(dtype=dtype)

    # Add per-multiplicity coordinates and masks
    pred_coords_unflat = pred_atom_coords_global_host.unflatten(0, (size_batch, multiplicity))
    true_coords_unflat = true_atom_coords_global_host.unflatten(0, (size_batch, multiplicity))
    resolved_mask_unflat = true_coords_resolved_mask_global_host.unflatten(0, (size_batch, multiplicity))

    for i_mul in range(multiplicity):
        inputs_atom[f"pred_atom_coords_{i_mul}"] = pred_coords_unflat[:, i_mul].to(dtype=dtype)
        inputs_atom[f"true_atom_coords_{i_mul}"] = true_coords_unflat[:, i_mul].to(dtype=dtype)
        inputs_atom[f"true_coords_resolved_mask_{i_mul}"] = resolved_mask_unflat[:, i_mul].to(dtype=dtype)

    # Define placements for CP submesh and full mesh
    placements_cp = {
        "atom_counts_per_token": (Shard(0), Replicate()),
        "token_to_rep_atom": (Shard(0), Replicate()),
        "r_set_to_rep_atom": (Shard(0), Replicate()),
        "atom_to_token": (Shard(0), Replicate()),
    }
    placements_dp_cp = {
        "token_to_rep_atom": (Shard(0), Shard(1), Replicate()),
        "r_set_to_rep_atom": (Shard(0), Shard(1), Replicate()),
        "atom_to_token": (Shard(0), Shard(1), Replicate()),
    }
    if pae_feats_host is not None:
        placements_cp["frames_idx"] = (Shard(1), Replicate())
        placements_cp["atom_pad_mask"] = (Shard(0), Replicate())
        placements_cp["atom_resolved_mask"] = (Shard(0), Replicate())
        placements_dp_cp["frames_idx"] = (Shard(0), Shard(1), Replicate())
        placements_dp_cp["atom_pad_mask"] = (Shard(0), Shard(1), Replicate())
        placements_dp_cp["atom_resolved_mask"] = (Shard(0), Shard(1), Replicate())
    for i_mul in range(multiplicity):
        placements_cp[f"pred_atom_coords_{i_mul}"] = (Shard(0), Replicate())
        placements_cp[f"true_atom_coords_{i_mul}"] = (Shard(0), Replicate())
        placements_cp[f"true_coords_resolved_mask_{i_mul}"] = (Shard(0), Replicate())
        placements_dp_cp[f"pred_atom_coords_{i_mul}"] = (Shard(0), Shard(1), Replicate())
        placements_dp_cp[f"true_atom_coords_{i_mul}"] = (Shard(0), Shard(1), Replicate())
        placements_dp_cp[f"true_coords_resolved_mask_{i_mul}"] = (Shard(0), Shard(1), Replicate())

    # Distribute atom features with intersperse padding
    feats_atom = distribute_atom_features(
        inputs_atom,
        placements_cp,
        placements_dp_cp,
        device_mesh,
        manager.group["cp"],
        multiplicities={
            "pred_atom_coords": multiplicity,
            "true_atom_coords": multiplicity,
            "true_coords_resolved_mask": multiplicity,
        },
    )

    # --- Distribute token features using distribute_features ---
    if manager.group_rank["world"] == 0:
        token_features = {
            "mol_type": mol_type_global_host.to(device=manager.device, dtype=torch.int64),
            "token_pad_mask": token_pad_mask_global_host.to(device=manager.device, dtype=dtype),
            "pred_lddt": pred_lddt_global_host.to(device=manager.device, dtype=torch.float32),
            "pred_resolved": pred_resolved_global_host.to(device=manager.device, dtype=torch.float32),
        }
        if pae_feats_host is not None:
            token_features["frame_resolved_mask"] = pae_feats_host["frame_resolved_mask"].to(
                device=manager.device, dtype=dtype
            )
            token_features["asym_id"] = pae_feats_host["asym_id"].to(device=manager.device)
            token_features["is_nonpolymer_with_frame"] = pae_feats_host["is_nonpolymer_with_frame"].to(
                device=manager.device
            )
    else:
        token_features = None
    token_placements = {
        "mol_type": (Shard(0), Shard(1), Replicate()),
        "token_pad_mask": (Shard(0), Shard(1), Replicate()),
        "pred_lddt": (Shard(0), Shard(1), Replicate()),
        "pred_resolved": (Shard(0), Shard(1), Replicate()),
    }
    if pae_feats_host is not None:
        token_placements["frame_resolved_mask"] = (Shard(0), Shard(1), Replicate())
        token_placements["asym_id"] = (Shard(0), Shard(1), Replicate())
        token_placements["is_nonpolymer_with_frame"] = (Shard(0), Shard(1), Replicate())
    token_feats_dtensor = distribute_features(
        token_features,
        token_placements,
        manager.group["world"],
        manager.group_ranks["world"][0],
        device_mesh,
    )

    # --- Distribute pair representations (pred_pde, and pred_pae when alpha_pae > 0) ---
    if manager.group_rank["world"] == 0:
        pair_features = {
            "pred_pde": pred_pde_global_host.to(device=manager.device, dtype=torch.float32),
        }
        if alpha_pae > 0.0:
            pair_features["pred_pae"] = pred_pae_global_host.to(device=manager.device, dtype=torch.float32)
    else:
        pair_features = None
    pair_placements = {
        "pred_pde": (Shard(0), Shard(1), Shard(2)),
    }
    if alpha_pae > 0.0:
        pair_placements["pred_pae"] = (Shard(0), Shard(1), Shard(2))
    pair_dtensor_dict = distribute_features(
        pair_features,
        pair_placements,
        manager.group["world"],
        manager.group_ranks["world"][0],
        device_mesh,
    )

    # Extract distributed tensors
    pred_atom_coords_dtensor = feats_atom["pred_atom_coords"]
    true_atom_coords_dtensor = feats_atom["true_atom_coords"]
    true_coords_resolved_mask_dtensor = feats_atom["true_coords_resolved_mask"]

    # Create feature dictionary for confidence_loss
    feats_dtensor = {
        "token_to_rep_atom": feats_atom["token_to_rep_atom"],
        "r_set_to_rep_atom": feats_atom["r_set_to_rep_atom"],
        "atom_to_token": feats_atom["atom_to_token"],
        "mol_type": token_feats_dtensor["mol_type"],
        "token_pad_mask": token_feats_dtensor["token_pad_mask"],
    }
    if pae_feats_host is not None:
        feats_dtensor["frames_idx"] = feats_atom["frames_idx"]
        feats_dtensor["atom_pad_mask"] = feats_atom["atom_pad_mask"]
        feats_dtensor["atom_resolved_mask"] = feats_atom["atom_resolved_mask"]
        feats_dtensor["frame_resolved_mask"] = token_feats_dtensor["frame_resolved_mask"]
        feats_dtensor["asym_id"] = token_feats_dtensor["asym_id"]
        feats_dtensor["is_nonpolymer_with_frame"] = token_feats_dtensor["is_nonpolymer_with_frame"]

    # Get model_out tensors with gradient tracking
    pred_lddt_dtensor = token_feats_dtensor["pred_lddt"].detach().requires_grad_(True)
    pred_pde_dtensor = pair_dtensor_dict["pred_pde"].detach().requires_grad_(True)
    pred_resolved_dtensor = token_feats_dtensor["pred_resolved"].detach().requires_grad_(True)

    # Build model_out dictionary
    model_out = {
        "plddt_logits": pred_lddt_dtensor,
        "pde_logits": pred_pde_dtensor,
        "resolved_logits": pred_resolved_dtensor,
        "sample_atom_coords": pred_atom_coords_dtensor,
    }

    pred_pae_dtensor = None
    if alpha_pae > 0.0:
        pred_pae_dtensor = pair_dtensor_dict["pred_pae"].detach().requires_grad_(True)
        model_out["pae_logits"] = pred_pae_dtensor

    # Compute confidence_loss
    confidence_loss_kwargs = {
        "model_out": model_out,
        "feats": feats_dtensor,
        "true_coords": true_atom_coords_dtensor,
        "true_coords_resolved_mask": true_coords_resolved_mask_dtensor,
        "comm": comm,
        "multiplicity": multiplicity,
        "alpha_pae": alpha_pae,
    }
    if alpha_pae > 0.0:
        confidence_loss_kwargs["dist_manager"] = manager
        confidence_loss_kwargs["group_layout"] = manager.layout_subgroups["cp"]

    result = confidence_loss(**confidence_loss_kwargs)

    # Verify output structure
    assert "loss" in result, "Result must contain 'loss' key"
    assert "loss_breakdown" in result, "Result must contain 'loss_breakdown' key"
    assert "plddt_loss" in result["loss_breakdown"], "loss_breakdown must contain 'plddt_loss'"
    assert "pde_loss" in result["loss_breakdown"], "loss_breakdown must contain 'pde_loss'"
    assert "resolved_loss" in result["loss_breakdown"], "loss_breakdown must contain 'resolved_loss'"
    assert "pae_loss" in result["loss_breakdown"], "loss_breakdown must contain 'pae_loss'"

    # Verify individual loss values and placements first for better diagnostics
    expected_placements = (Replicate(), Replicate(), Replicate())
    for loss_name in ["plddt_loss", "pde_loss", "resolved_loss", "pae_loss"]:
        subloss_dtensor = result["loss_breakdown"][loss_name]
        assert (
            subloss_dtensor.placements == expected_placements
        ), f"{loss_name} placements {subloss_dtensor.placements} != expected {expected_placements}"
        loss_value = subloss_dtensor.to_local()
        expected_value = expected_loss_breakdown_host[loss_name].to(device=loss_value.device, dtype=loss_value.dtype)
        torch.testing.assert_close(loss_value, expected_value, msg=f"{loss_name} mismatch")

    # Verify total loss value
    loss_local = result["loss"].to_local()
    expected_loss = expected_loss_host.to(device=loss_local.device, dtype=loss_local.dtype)
    torch.testing.assert_close(loss_local, expected_loss)

    # Verify total loss placements are fully replicated
    loss_dtensor = result["loss"]
    assert (
        loss_dtensor.placements == expected_placements
    ), f"Loss placements {loss_dtensor.placements} != expected {expected_placements}"

    # Backward pass on DTensor directly
    loss_dtensor.backward()

    # Verify gradients for pred_lddt
    grad_pred_lddt = pred_lddt_dtensor.grad
    assert grad_pred_lddt is not None, "Gradient not computed for pred_lddt"
    grad_pred_lddt_full = grad_pred_lddt.full_tensor()
    expected_grad_lddt = expected_grad_pred_lddt_host.to(
        device=grad_pred_lddt_full.device, dtype=grad_pred_lddt_full.dtype
    )
    torch.testing.assert_close(grad_pred_lddt_full, expected_grad_lddt)

    # Verify gradients for pred_pde
    grad_pred_pde = pred_pde_dtensor.grad
    assert grad_pred_pde is not None, "Gradient not computed for pred_pde"
    grad_pred_pde_full = grad_pred_pde.full_tensor()
    expected_grad_pde = expected_grad_pred_pde_host.to(device=grad_pred_pde_full.device, dtype=grad_pred_pde_full.dtype)
    torch.testing.assert_close(grad_pred_pde_full, expected_grad_pde)

    # Verify gradients for pred_resolved
    grad_pred_resolved = pred_resolved_dtensor.grad
    assert grad_pred_resolved is not None, "Gradient not computed for pred_resolved"
    grad_pred_resolved_full = grad_pred_resolved.full_tensor()
    expected_grad_resolved = expected_grad_pred_resolved_host.to(
        device=grad_pred_resolved_full.device, dtype=grad_pred_resolved_full.dtype
    )
    torch.testing.assert_close(grad_pred_resolved_full, expected_grad_resolved)

    # Verify gradients for pred_pae (when alpha_pae > 0)
    if alpha_pae > 0.0:
        grad_pred_pae = pred_pae_dtensor.grad
        assert grad_pred_pae is not None, "Gradient not computed for pred_pae"
        grad_pred_pae_full = grad_pred_pae.full_tensor()
        expected_grad_pae = expected_grad_pred_pae_host.to(
            device=grad_pred_pae_full.device, dtype=grad_pred_pae_full.dtype
        )
        torch.testing.assert_close(grad_pred_pae_full, expected_grad_pae)

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env",
    [
        ((1, (2, 2)), True, "cuda", "ENV"),
        ((2, (2, 2)), True, "cuda", "ENV"),
    ],
    indirect=("setup_env",),
    ids=lambda x: f"dp:{x[0][0]}, cp:{x[0][1]}, device_type:{x[2]}",
)
@pytest.mark.parametrize("alpha_pae", [0.0, 0.5], ids=lambda x: f"alpha_pae:{x}")
@pytest.mark.parametrize("multiplicity", [1, 2], ids=lambda x: f"multiplicity:{x}")
def test_confidence_loss(setup_env, alpha_pae: float, multiplicity: int):
    """Test that DTensor confidence_loss matches serial reference.

    This test verifies:
    1. Forward pass: DTensor confidence_loss matches serial confidence_loss
    2. Backward pass: Gradients w.r.t. pred_lddt, pred_pde, pred_resolved (and pred_pae
       when alpha_pae > 0) match serial gradients
    3. Output structure contains plddt_loss, pde_loss, resolved_loss, pae_loss

    The test uses realistic feature generation with proper block-diagonal structure.
    """
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() < {world_size}")

    torch.manual_seed(42)

    dp_size = grid_group_sizes["dp"]
    cp_size = grid_group_sizes["cp"][0] * grid_group_sizes["cp"][1]

    # Generate test data
    B = dp_size  # Batch size equals DP size
    N_token = 32
    N_atom = 140  # Large enough to accommodate token atom counts
    n_atoms_per_token_min = 1
    n_atoms_per_token_max = 4
    num_bins_lddt = 50  # Number of pLDDT bins
    num_bins_pde = 64  # Number of PDE bins
    dtype = torch.float32  # Use FP32 for testing

    # Make N_atom divisible by cp_size for even sharding
    N_atom = ((N_atom + cp_size - 1) // cp_size) * cp_size

    # Use random_features to generate features with proper block-diagonal structure
    selected_keys = [
        "token_to_rep_atom",
        "r_set_to_rep_atom",
        "atom_to_token",
        "mol_type",
        "token_pad_mask",
        "atom_counts_per_token",
    ]
    if alpha_pae > 0.0:
        selected_keys += [
            "frames_idx",
            "frame_resolved_mask",
            "asym_id",
            "atom_pad_mask",
            "atom_resolved_mask",
            "is_nonpolymer_with_frame",
        ]

    feats = random_features(
        size_batch=B,
        n_tokens=N_token,
        n_atoms=N_atom,
        n_msa=1,
        atom_counts_per_token_range=(n_atoms_per_token_min, n_atoms_per_token_max),
        device=torch.device(device_type),
        float_value_range=(-0.5, 0.5),
        selected_keys=selected_keys,
    )

    if alpha_pae > 0.0:
        feats["atom_resolved_mask"] = torch.ones_like(feats["atom_resolved_mask"])
        feats["frame_resolved_mask"] = torch.ones_like(feats["frame_resolved_mask"])

    token_to_rep_atom_global = feats["token_to_rep_atom"].to(dtype=dtype)  # [B, N_token, N_atom]
    r_set_to_rep_atom_global = feats["r_set_to_rep_atom"].to(dtype=dtype)  # [B, N_R, N_atom]
    atom_to_token_global = feats["atom_to_token"].to(dtype=dtype)  # [B, N_atom, N_token]
    mol_type_global = feats["mol_type"]  # [B, N_token]
    token_pad_mask_global = feats["token_pad_mask"].to(dtype=dtype)  # [B, N_token]
    atom_counts_per_token = feats["atom_counts_per_token"]

    N_atom_actual = token_to_rep_atom_global.shape[2]

    # Generate coordinates using uniform distribution
    # Use [-10, 10] range for meaningful pairwise distances
    pred_atom_coords_global = torch.empty(B * multiplicity, N_atom_actual, 3, device=device_type, dtype=dtype).uniform_(
        -10.0, 10.0
    )
    true_atom_coords_global = torch.empty(B * multiplicity, N_atom_actual, 3, device=device_type, dtype=dtype).uniform_(
        -10.0, 10.0
    )

    # Create true_coords_resolved_mask: [B*mult, N_atom]
    true_coords_resolved_mask_global = torch.randint(
        0, 2, (B * multiplicity, N_atom_actual), device=device_type, dtype=dtype
    )
    # Ensure at least half atoms are resolved for meaningful test
    true_coords_resolved_mask_global[:, : N_atom_actual // 2] = 1.0

    # Generate pred_lddt: (B*mult, N_token, num_bins_lddt)
    pred_lddt_global = torch.randn(
        B * multiplicity, N_token, num_bins_lddt, device=device_type, dtype=torch.float32
    ).requires_grad_(True)

    # Generate pred_pde: (B*mult, N_token, N_token, num_bins_pde)
    pred_pde_global = torch.randn(
        B * multiplicity, N_token, N_token, num_bins_pde, device=device_type, dtype=torch.float32
    ).requires_grad_(True)

    # Generate pred_resolved: (B*mult, N_token, 2)
    pred_resolved_global = torch.randn(
        B * multiplicity, N_token, 2, device=device_type, dtype=torch.float32
    ).requires_grad_(True)

    # Generate pred_pae when needed: (B*mult, N_token, N_token, num_bins_pde)
    num_bins_pae = 64
    pred_pae_global = None
    if alpha_pae > 0.0:
        pred_pae_global = torch.randn(
            B * multiplicity, N_token, N_token, num_bins_pae, device=device_type, dtype=torch.float32
        ).requires_grad_(True)

    # Build model_out for serial reference
    model_out_serial = {
        "plddt_logits": pred_lddt_global,
        "pde_logits": pred_pde_global,
        "resolved_logits": pred_resolved_global,
        "sample_atom_coords": pred_atom_coords_global.float(),
    }
    if alpha_pae > 0.0:
        model_out_serial["pae_logits"] = pred_pae_global

    # Build feats for serial reference
    feats_serial = {
        "token_to_rep_atom": token_to_rep_atom_global.clone().float(),
        "r_set_to_rep_atom": r_set_to_rep_atom_global.clone().float(),
        "atom_to_token": atom_to_token_global.clone().float(),
        "mol_type": mol_type_global.clone(),
        "token_pad_mask": token_pad_mask_global.clone().float(),
    }
    if alpha_pae > 0.0:
        feats_serial["frames_idx"] = feats["frames_idx"].clone()
        feats_serial["frame_resolved_mask"] = feats["frame_resolved_mask"].clone().float()
        feats_serial["asym_id"] = feats["asym_id"].clone()
        feats_serial["atom_pad_mask"] = feats["atom_pad_mask"].clone().float()
        feats_serial["atom_resolved_mask"] = feats["atom_resolved_mask"].clone().float()
        feats_serial["is_nonpolymer_with_frame"] = feats["is_nonpolymer_with_frame"].clone()

    # Compute serial reference
    expected_result = serial_confidence_loss(
        model_out=model_out_serial,
        feats=feats_serial,
        true_coords=true_atom_coords_global.float(),
        true_coords_resolved_mask=true_coords_resolved_mask_global.float(),
        token_level_confidence=True,
        multiplicity=multiplicity,
        alpha_pae=alpha_pae,
    )

    expected_loss = expected_result["loss"]
    expected_loss_breakdown = {
        "plddt_loss": expected_result["loss_breakdown"]["plddt_loss"],
        "pde_loss": expected_result["loss_breakdown"]["pde_loss"],
        "resolved_loss": expected_result["loss_breakdown"]["resolved_loss"],
        "pae_loss": torch.tensor(expected_result["loss_breakdown"]["pae_loss"], dtype=dtype),
    }

    # Compute gradients for serial reference
    expected_loss.backward()
    expected_grad_pred_lddt = pred_lddt_global.grad.clone()
    expected_grad_pred_pde = pred_pde_global.grad.clone()
    expected_grad_pred_resolved = pred_resolved_global.grad.clone()
    expected_grad_pred_pae = None
    if alpha_pae > 0.0:
        expected_grad_pred_pae = pred_pae_global.grad.clone()
        assert expected_grad_pred_pae is not None, "Serial confidence_loss should produce gradients for pred_pae"

    # Verify that serial reference produces gradients
    assert expected_grad_pred_lddt is not None, "Serial confidence_loss should produce gradients for pred_lddt"
    assert expected_grad_pred_pde is not None, "Serial confidence_loss should produce gradients for pred_pde"
    assert expected_grad_pred_resolved is not None, "Serial confidence_loss should produce gradients for pred_resolved"

    # Collect PAE-specific features for the parallel worker
    pae_feats = None
    if alpha_pae > 0.0:
        pae_feats = {
            "frames_idx": feats["frames_idx"].clone().cpu(),
            "frame_resolved_mask": feats["frame_resolved_mask"].clone().cpu(),
            "asym_id": feats["asym_id"].clone().cpu(),
            "atom_pad_mask": feats["atom_pad_mask"].clone().cpu(),
            "atom_resolved_mask": feats["atom_resolved_mask"].clone().cpu(),
            "is_nonpolymer_with_frame": feats["is_nonpolymer_with_frame"].clone().cpu(),
        }

    # Prepare payload for parallel test
    payload = (
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        multiplicity,
        alpha_pae,
        pred_lddt_global.detach().clone().cpu(),
        pred_pde_global.detach().clone().cpu(),
        pred_pae_global.detach().clone().cpu() if pred_pae_global is not None else None,
        pred_resolved_global.detach().clone().cpu(),
        pred_atom_coords_global.clone().cpu(),
        true_atom_coords_global.clone().cpu(),
        true_coords_resolved_mask_global.clone().cpu(),
        token_to_rep_atom_global.clone().cpu(),
        r_set_to_rep_atom_global.clone().cpu(),
        atom_to_token_global.clone().cpu(),
        mol_type_global.clone().cpu(),
        token_pad_mask_global.clone().cpu(),
        atom_counts_per_token.clone().cpu(),
        pae_feats,
        expected_loss.detach().clone().cpu(),
        {k: v.detach().clone().cpu() if isinstance(v, torch.Tensor) else v for k, v in expected_loss_breakdown.items()},
        expected_grad_pred_lddt.clone().cpu(),
        expected_grad_pred_pde.clone().cpu(),
        expected_grad_pred_pae.clone().cpu() if expected_grad_pred_pae is not None else None,
        expected_grad_pred_resolved.clone().cpu(),
    )

    # Launch parallel test
    spawn_multiprocessing(parallel_assert_confidence_loss, world_size, payload)


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
