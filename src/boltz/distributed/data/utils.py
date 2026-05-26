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


from collections import OrderedDict
from typing import NamedTuple, Optional

import torch
from torch.distributed import ProcessGroup
from torch.distributed.tensor import DTensor, Placement, Replicate, Shard, distribute_tensor
from torch.distributed.tensor.device_mesh import DeviceMesh
from torch.types import Tensor

from boltz.data.pad import pad_to_max
from boltz.distributed.data.feature.featurizer_utils import ATOM_INDEX_FEATURES, remap_atom_indices_repad
from boltz.distributed.manager import DistributedManager
from boltz.distributed.utils import LayoutRightMap

# Placement types used inside the DataSet class which only sees the NxN CP group
PLACEMENT_TYPE_SHARD0_REPLICATE = (Shard(0), Replicate())  # For single repr
PLACEMENT_TYPE_SHARD0_SHARD1 = (Shard(0), Shard(1))  # For pair repr
PLACEMENT_TYPE_SHARD1_REPLICATE = (Shard(1), Replicate())  # For ensemble-aware features


ATOM_FEATURES = {
    "ref_pos",
    "atom_resolved_mask",
    "ref_element",
    "ref_charge",
    "ref_atom_name_chars",
    "ref_space_uid",
    "coords",
    "atom_pad_mask",
    "atom_to_token",
    "atom_counts_per_token",
    "pair_mask",
    "token_to_rep_atom",
    "r_set_to_rep_atom",
    "frames_idx",
}

ATOM_FEATURES_V2 = ATOM_FEATURES | {
    "ref_chirality",
    "atom_backbone_feat",
    "bfactor",
    "plddt",
}

SYMMETRY_FEATURES = {
    "all_coords",
    "all_resolved_mask",
    "crop_to_all_atom_map",
    "chain_symmetries",
    "amino_acids_symmetries",
    "ligand_symmetries",
}
SYMMETRY_FEATURES_V2 = SYMMETRY_FEATURES | {
    "chain_swaps",
}

RESIDUE_CONSTRAINT_FEATURES = {
    "rdkit_bounds_index",
    "rdkit_bounds_bond_mask",
    "rdkit_bounds_angle_mask",
    "rdkit_upper_bounds",
    "rdkit_lower_bounds",
    "chiral_atom_index",
    "chiral_reference_mask",
    "chiral_atom_orientations",
    "stereo_bond_index",
    "stereo_reference_mask",
    "stereo_bond_orientations",
    "planar_bond_index",
    "planar_ring_5_index",
    "planar_ring_6_index",
}

CHAIN_CONSTRAINT_FEATURES = {
    "connected_chain_index",
    "connected_atom_index",
    "symmetric_chain_index",
}

CONTACT_CONSTRAINT_FEATURES = {
    "contact_pair_index",
    "contact_union_index",
    "contact_negation_mask",
    "contact_thresholds",
}

NON_SHARDED_FEATURES = {"record"} | SYMMETRY_FEATURES | RESIDUE_CONSTRAINT_FEATURES | CHAIN_CONSTRAINT_FEATURES
TRAINING_METADATA_FEATURES = {
    "chain_swaps",
    "activity_name",
    "activity_qualifier",
    "idx_dataset",
    "sid",
    "cid",
    "normalized_protein_accession",
    "pair_id",
    "pdb_id",
}

LIGAND_GEOMETRY_FEATURES = {
    "ligand_edge_index",
    "ligand_edge_lower_bounds",
    "ligand_edge_upper_bounds",
    "ligand_edge_bond_mask",
    "ligand_edge_angle_mask",
    "connections_edge_index",
    "ligand_chiral_atom_index",
    "ligand_chiral_check_mask",
    "ligand_chiral_atom_orientations",
    "ligand_stereo_bond_index",
    "ligand_stereo_check_mask",
    "ligand_stereo_bond_orientations",
    "ligand_aromatic_5_ring_index",
    "ligand_aromatic_6_ring_index",
    "ligand_planar_double_bond_index",
}

NON_SHARDED_FEATURES_V2 = (
    {"record"}
    | SYMMETRY_FEATURES_V2
    | RESIDUE_CONSTRAINT_FEATURES
    | CHAIN_CONSTRAINT_FEATURES
    | {
        "affinity_mw",
        "ensemble_ref_idxs",
        "template_force",
        "template_force_threshold",
    }
    | CONTACT_CONSTRAINT_FEATURES
    | TRAINING_METADATA_FEATURES
    | LIGAND_GEOMETRY_FEATURES
)


class TensorMetadata(NamedTuple):
    dtype: torch.dtype
    shape: torch.Size


def broadcast_feature_tensors_metadata(
    features: dict[str, Tensor] | None, group: ProcessGroup, src_rank_global: int
) -> OrderedDict[str, tuple[torch.dtype, torch.Size]]:
    """Broadcast tensor metadata from source rank to all ranks in a process group.

    This function extracts metadata (dtype and shape) from feature tensors on the source
    rank, broadcasts this information to all other ranks in the process group, and returns
    it as an OrderedDict to preserve feature ordering consistency across ranks. This allows
    non-source ranks to know the structure of tensors without having the actual tensor data,
    and ensures all ranks iterate over features in the same order.

    Parameters
    ----------
    features : dict[str, Tensor] | None
        Dictionary mapping feature names to tensors. Must be a dict on the source rank
        and None on all other ranks. Feature names can include random prefixes to test
        ordering consistency across distributed ranks.
    group : ProcessGroup
        The distributed process group to broadcast within.
    src_rank_global : int
        The global rank that serves as the source for broadcasting metadata.
        This rank must be included in the process group.

    Returns
    -------
    OrderedDict[str, tuple[torch.dtype, torch.Size]]
        Ordered dictionary mapping feature names to tuples of (dtype, shape) where dtype
        is the PyTorch data type and shape is the tensor's size. The order is preserved
        from the source rank to ensure consistent iteration order across all ranks.

    Raises
    ------
    ValueError
        If source rank doesn't provide a dict input, non-source ranks provide non-None
        input, or feature values are not tensors.

    Notes
    -----
    - Returns an OrderedDict to guarantee consistent feature iteration order across ranks
    - This is critical for distributed processing where all ranks must process features
      in the same sequence to maintain synchronization
    - The metadata includes TensorMetadata named tuples containing dtype and shape info
    - Feature names may include random prefixes for testing ordering robustness
    """
    rank_global = torch.distributed.get_rank()
    is_src_rank = rank_global == src_rank_global

    # Check that src_rank_global is in the process group
    ranks_in_group = torch.distributed.get_process_group_ranks(group)
    if src_rank_global not in ranks_in_group:
        raise ValueError(f"Source rank {src_rank_global} not in group {ranks_in_group}")

    if is_src_rank and not isinstance(features, dict):
        raise ValueError(f"Source rank {src_rank_global} must have a dict input, got {type(features)}")
    elif not is_src_rank and features is not None:
        raise ValueError(f"Non-source rank (not {src_rank_global}) must have None input, this rank is {rank_global}")
    if is_src_rank:
        metadata = OrderedDict()
        for k, v in features.items():
            if not isinstance(v, Tensor):
                raise ValueError(f"Feature {k} is not a tensor, got {type(v)}")
            metadata[k] = TensorMetadata(dtype=v.dtype, shape=v.shape)
    else:
        metadata = None
    l_metadata = [metadata]
    torch.distributed.broadcast_object_list(l_metadata, src=src_rank_global, group=group)
    metadata = l_metadata[0]
    return metadata


def distribute_features(
    features: dict[str, Tensor],
    placements: dict[str, Placement],
    group: ProcessGroup,
    src_rank_global: int,
    device_mesh: DeviceMesh,
) -> dict[str, DTensor]:
    """Distribute feature tensors from source rank to all ranks as DTensors.

    This function takes feature tensors from a source rank and distributes them across
    a device mesh according to specified placements. The metadata (dtype and shape) is
    first broadcast to all ranks using broadcast_feature_tensors_metadata (which preserves
    ordering), then each rank creates or uses existing tensors to form DTensors with the
    specified placements. All ranks process features in the same order for consistency.

    Parameters
    ----------
    features : dict[str, Tensor]
        Dictionary mapping feature names to tensors. Must contain tensors on the source
        rank and can be None on non-source ranks (only validated on source rank).
        Feature names may include random prefixes for testing ordering robustness.
    placements : dict[str, Placement]
        Dictionary mapping feature names to their desired tensor placements (e.g., Shard,
        Replicate). Must have the same keys as features on the source rank.
    group : ProcessGroup
        The distributed process group that contains all ranks in the device mesh.
    src_rank_global : int
        The global rank that serves as the source for the original tensors.
        This rank must be included in the process group.
    device_mesh : DeviceMesh
        The device mesh defining the distributed tensor layout. Its ranks must match
        the ranks in the process group.

    Returns
    -------
    dict[str, DTensor]
        Dictionary mapping feature names to distributed tensors (DTensors) with the
        specified placements across the device mesh. While returned as a regular dict,
        the internal processing ensures consistent ordering across ranks.

    Raises
    ------
    ValueError
        If features and placements don't have the same keys (validated only on source rank),
        if the ranks in the process group don't match the ranks in the device mesh,
        or if the source rank is not in the process group.

    Notes
    -----
    - Uses broadcast_feature_tensors_metadata to ensure consistent feature ordering across ranks
    - Non-source ranks create empty tensors with the correct dtype and shape from metadata
    - Each tensor is distributed using PyTorch's distribute_tensor function with specified placements
    - The resulting DTensors are ready for distributed computation with consistent sharding
    - Order preservation is critical for synchronization during distributed training
    - Returns regular dict for backward compatibility, but internal ordering is guaranteed

    Examples
    --------
    >>> features = {"abc123_feature_0": torch.randn(10, 20)}  # Only on source rank
    >>> placements = {"abc123_feature_0": Shard(0)}
    >>> dtensors = distribute_features(features, placements, group, 0, mesh)
    """
    rank_global = torch.distributed.get_rank()
    is_src_rank = rank_global == src_rank_global

    if is_src_rank and features.keys() != placements.keys():
        raise ValueError(
            f"Features and placements must have the same keys, got {sorted(features.keys())} and {sorted(placements.keys())}"
        )
    ranks_in_group = torch.distributed.get_process_group_ranks(group)
    ranks_in_mesh = device_mesh.mesh.flatten().tolist()
    if ranks_in_group != ranks_in_mesh:
        raise ValueError(
            f"Ranks in group {ranks_in_group} do not match ranks in mesh {ranks_in_mesh}, got {ranks_in_group} and {ranks_in_mesh}"
        )

    if src_rank_global not in ranks_in_group:
        raise ValueError(f"Source rank {src_rank_global} not in group {ranks_in_group}")

    # this guarantees only src_rank_global have the features so the returned metadata is ordered by the iteration order
    # of the features dictionary in the src_rank_global. This is important for the later for loop to iterate over the features
    # in the same order among all ranks
    metadata: OrderedDict[str, TensorMetadata] = broadcast_feature_tensors_metadata(features, group, src_rank_global)

    # on the other hand, to stay backward compatible with other usage in the data processing code, we return
    # regular dict and it's up to the caller to make sure the iteration order is consistent among ranks
    ans = {}
    for name_feature, m in metadata.items():
        dtype, shape = m
        placement = placements[name_feature]
        if is_src_rank:
            t = features[name_feature].to(device=device_mesh.device_type)
        else:
            t = torch.empty(shape, dtype=dtype, device=device_mesh.device_type)
        ans[name_feature] = distribute_tensor(t, device_mesh, placements=placement)
    return ans


def broadcast_tensors(
    features: dict[str, Tensor] | None,
    group: ProcessGroup,
    src_rank_global: int,
    device: str = "cpu",
) -> OrderedDict[str, Tensor]:
    """Broadcast tensors from source rank to all ranks in a process group.

    This function broadcasts tensor data from the source rank to all other ranks
    in the process group. The metadata (dtype and shape) is first broadcast using
    broadcast_feature_tensors_metadata, then each tensor is broadcast individually.
    All ranks process features in the same order (determined by the source rank's
    iteration order) to ensure consistent communication and avoid deadlocks.

    Parameters
    ----------
    features : dict[str, Tensor] | None
        Dictionary mapping feature names to tensors. Must be a dict containing tensors
        on the source rank and None on all other ranks.
    group : ProcessGroup
        The distributed process group to broadcast within.
    src_rank_global : int
        The global rank that serves as the source for broadcasting.
        This rank must be included in the process group.
    device : str
        Device to use for broadcasting. For NCCL backend, this should be "cuda".
        For gloo backend, this can be "cpu". Default is "cpu".

    Returns
    -------
    OrderedDict[str, Tensor]
        Ordered dictionary mapping feature names to broadcast tensors on the specified
        device. The order is preserved from the source rank to ensure consistent
        iteration order across all ranks.

    Raises
    ------
    ValueError
        If source rank doesn't provide a dict input, non-source ranks provide non-None
        input, or feature values are not tensors (raised by broadcast_feature_tensors_metadata).

    Notes
    -----
    - Uses broadcast_feature_tensors_metadata to ensure consistent feature ordering across ranks
    - Non-source ranks create empty tensors with the correct dtype and shape from metadata
    - Each tensor is broadcast using torch.distributed.broadcast
    - Output tensors are on the specified device
    - Order preservation is critical for synchronization during distributed processing

    Examples
    --------
    >>> # On source rank (rank 0)
    >>> features = {"coords": torch.randn(10, 3), "mask": torch.ones(10)}
    >>> result = broadcast_tensors(features, group, src_rank_global=0, device="cuda")
    >>> # On non-source ranks
    >>> result = broadcast_tensors(None, group, src_rank_global=0, device="cuda")
    >>> # All ranks now have the same tensors in result
    """
    rank_global = torch.distributed.get_rank()
    is_src_rank = rank_global == src_rank_global

    # broadcast_feature_tensors_metadata validates inputs and returns an OrderedDict
    # with consistent ordering across all ranks (based on source rank's iteration order)
    metadata: OrderedDict[str, TensorMetadata] = broadcast_feature_tensors_metadata(features, group, src_rank_global)

    # Iterate over metadata in consistent order across all ranks
    ans = OrderedDict()
    for name_feature, m in metadata.items():
        dtype, shape = m
        if is_src_rank:
            tensor = features[name_feature].to(device=device).contiguous()
        else:
            tensor = torch.empty(shape, dtype=dtype, device=device)
        torch.distributed.broadcast(tensor, src=src_rank_global, group=group)
        ans[name_feature] = tensor

    return ans


class CollateDTensor:
    def __init__(self, output_device_mesh: DeviceMesh):
        # Check that shape is like (dp, cp_axis_0, cp_axis_1)
        if output_device_mesh.ndim != 3:
            raise ValueError(f"CollateDTensor expects a DP-CP-CP device mesh but got ndim {output_device_mesh.ndim}")

        self._output_device_mesh = output_device_mesh

    def __call__(self, data: list[dict[str, DTensor]]) -> dict[str, DTensor]:
        """Collate the data.

        Parameters
        ----------
        data : List[Dict[str, DTensor]]
            The data to collate.

        Returns
        -------
        Dict[str, DTensor]
            The collated data.

        """
        # Get the keys
        keys = data[0].keys()

        # Pre-scan: determine final atom dim for atom-index remapping.
        # Atom-index features store indices into a padded global atom array
        # whose stride is max_atoms_per_shard. When collation pads the atom
        # dimension (because samples differ or DP ranks differ), that stride
        # changes and every stored index must be adjusted.
        _ATOM_DIM_REFERENCE_KEY = "atom_pad_mask"
        has_atom_index_features = any(k in ATOM_INDEX_FEATURES for k in keys)
        if has_atom_index_features and _ATOM_DIM_REFERENCE_KEY in keys:
            ref_locals = [d[_ATOM_DIM_REFERENCE_KEY].to_local() for d in data]
            old_atoms_per_sample = [v.shape[0] for v in ref_locals]
            local_max_atoms = max(old_atoms_per_sample)

            global_max_atoms_t = torch.tensor([local_max_atoms], device=self._output_device_mesh.device_type)
            group = self._output_device_mesh.get_group(0)
            torch.distributed.all_reduce(global_max_atoms_t, op=torch.distributed.ReduceOp.MAX, group=group)
            final_atoms_per_shard = int(global_max_atoms_t.item())
        else:
            old_atoms_per_sample = None
            final_atoms_per_shard = None

        # Collate the data
        collated = {}
        for key in keys:
            # special handling for non-DTensor features
            if key in NON_SHARDED_FEATURES_V2:
                collated[key] = [d[key] for d in data]
                continue

            # change batch dim to shard, special handling for coords since it has a leading singleton dim
            values = [d[key] for d in data]
            placements = values[0].placements

            values_local = [value.to_local() for value in values]

            # Remap atom-index features before padding so indices reflect the
            # final (post-collation) padded atom layout.
            if key in ATOM_INDEX_FEATURES and old_atoms_per_sample is not None:
                for i in range(len(values_local)):
                    values_local[i] = remap_atom_indices_repad(
                        values_local[i], old_atoms_per_sample[i], final_atoms_per_shard
                    )

            # local collate
            values_local, _ = pad_to_max(values_local, 0)  # internally will stack if shapes are the same
            values_local = values_local.contiguous()  # contiguous implies layout right

            # global collate
            local_shape_max = torch.tensor(
                values_local.shape,
                device=self._output_device_mesh.device_type,
            )
            group = self._output_device_mesh.get_group(0)
            torch.distributed.all_reduce(local_shape_max, op=torch.distributed.ReduceOp.MAX, group=group)

            # Pad local tensor to match global shape if needed
            if values_local.shape != tuple(local_shape_max.tolist()):
                # Calculate padding needed for each dimension
                current_shape = values_local.shape
                target_shape = tuple(local_shape_max.tolist())
                num_dims = len(current_shape)

                # Create padding tuple (reverse order for F.pad)
                padding = []
                for i in range(num_dims):
                    dim_idx = num_dims - 1 - i  # Reverse order for F.pad
                    pad_needed = target_shape[dim_idx] - current_shape[dim_idx]
                    padding.extend([0, pad_needed])  # [left, right] for each dimension

                # Pad the tensor
                values_local = torch.nn.functional.pad(values_local, tuple(padding), value=0)

            # expand to get global shape and strides
            shape_scaling = torch.ones_like(local_shape_max)
            shape_scaling[0] = self._output_device_mesh.shape[0]  # dp dimension

            new_placements = [Shard(0)]
            for placement in placements:
                if isinstance(placement, Shard):
                    shape_scaling[placement.dim + 1] = self._output_device_mesh.shape[
                        placement.dim + 1
                    ]  # +1 because input placement lacks dp dimension
                    new_placements.append(Shard(placement.dim + 1))
                elif isinstance(placement, Replicate):
                    new_placements.append(placement)
                else:
                    raise ValueError(f"Unsupported placement: {placement}")

            global_shape = local_shape_max * shape_scaling
            strides = LayoutRightMap(tuple(global_shape.tolist())).strides  # coherent with contiguous layout

            collated[key] = DTensor.from_local(
                values_local,
                device_mesh=self._output_device_mesh,
                placements=new_placements,
                shape=torch.Size(global_shape.tolist()),
                stride=strides,
            )

        return collated


def map_subgroup_mesh_to_cpu(dist_manager: "DistributedManager") -> DeviceMesh:
    """Map the subgroup mesh to the CPU device mesh.

    Parameters
    ----------
    manager : DistributedManager
        The distributed manager.
    """
    device_mesh = dist_manager.device_mesh_subgroups
    if device_mesh.device_type == "cpu":
        return DeviceMesh.from_group(
            group=[
                dist_manager.group["dp"],
                dist_manager.group["cp_axis_0"],
                dist_manager.group["cp_axis_1"],
            ],
            device_type="cpu",
            mesh=device_mesh.mesh.clone(),
            mesh_dim_names=(
                "dp_cpu",
                "cp_axis_0_cpu",
                "cp_axis_1_cpu",
            ),
        )
    elif device_mesh.device_type == "cuda":
        if "dp_cpu" not in dist_manager.group_ranks:
            dist_manager.create_group(
                "dp_cpu",
                dist_manager.group_ranks["dp"],
                backend="gloo",
                use_local_synchronization=True,
            )
            dist_manager.create_group(
                "cp_axis_0_cpu",
                dist_manager.group_ranks["cp_axis_0"],
                backend="gloo",
                use_local_synchronization=True,
            )
            dist_manager.create_group(
                "cp_axis_1_cpu",
                dist_manager.group_ranks["cp_axis_1"],
                backend="gloo",
                use_local_synchronization=True,
            )

        device_mesh_cpu = DeviceMesh.from_group(
            group=[
                dist_manager.group["dp_cpu"],
                dist_manager.group["cp_axis_0_cpu"],
                dist_manager.group["cp_axis_1_cpu"],
            ],
            device_type="cpu",
            mesh=device_mesh.mesh.clone(),
            mesh_dim_names=(
                "dp_cpu",
                "cp_axis_0_cpu",
                "cp_axis_1_cpu",
            ),
        )

        # Check if cpu and cuda group ranks match
        cuda_group_ranks = (
            dist_manager.group_ranks["dp"],
            dist_manager.subgroups_ranks["cp"][0],
            dist_manager.subgroups_ranks["cp"][1],
        )
        cpu_group_ranks = tuple(
            torch.distributed.get_process_group_ranks(group) for group in device_mesh_cpu.get_all_groups()
        )
        for cuda_group_ranks_, cpu_group_ranks_ in zip(cuda_group_ranks, cpu_group_ranks):
            if set(cuda_group_ranks_) != set(cpu_group_ranks_):
                raise ValueError(
                    f"New CPU group ranks {cpu_group_ranks_} do not match with existing CUDA group ranks {cuda_group_ranks_}"
                )

        return device_mesh_cpu
    else:
        raise ValueError(f"Unknown device type {device_mesh.device_type}")


def get_flattened_group(device_mesh: DeviceMesh, backend: Optional[str] = None) -> ProcessGroup:
    """Get the flattened process group from a device mesh.

    The original _flatten method creates a new group using default parameters, which can lead to
    inconsistent backend with the original mesh. This function creates an additional group with a
    pre-specified backend.

    Examples:
        mesh = DeviceMesh(device_type="cuda", mesh=[[0, 1], [2, 3]])
        group = get_flattened_group(mesh, backend="gloo")
        print(group)
        # ProcessGroup(type: GLOO, backend: gloo, devices: [0, 1, 2, 3])

    Parameters
    ----------
    device_mesh : DeviceMesh
        The device mesh.

    Returns
    -------
    group : ProcessGroup
        The flattened group.
    """
    new_group = device_mesh._flatten().get_group()
    if backend is None:
        return new_group

    return torch.distributed.new_group(
        torch.distributed.get_process_group_ranks(new_group),
        backend=backend,
        use_local_synchronization=True,
    )
