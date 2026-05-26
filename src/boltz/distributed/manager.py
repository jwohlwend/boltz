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


import os
from copy import deepcopy
from math import prod
from typing import Any, Dict, Optional, OrderedDict, Union
from warnings import warn

import torch

from boltz.distributed.utils import LayoutMap, LayoutRightMap

# grid_group_sizes objects must have
#   (1) .values() attribute, (2) .items() attribute
_GridGroupSizesType = OrderedDict[str, Union[int, tuple[int, ...]]]


class DistributedManager:
    """
    Borg-style singleton class for managing distributed state.

    Parameters
    ----------
    None

    Attributes
    ----------
    rank : int
        Rank of the process.
    world_size : int
        Total number of processes.
    local_rank : int
        Local rank of the process.
    device : torch.device
        Device used by the process.
    backend : str
        Backend used for distributed communication.
    device_mesh : torch.distributed.device_mesh
        Device mesh used for distributed communication.
    group : dict
        Dictionary of groups used for distributed communication.
    group_rank : dict
        Dictionary of group ranks used for distributed communication.
    group_ranks : dict
        Dictionary of group ranks used for distributed communication.
    method_init : str
        Method used to initialize the distributed manager.

    Examples
    --------
    >>> DistributedManager.initialize()
    >>> manager = DistributedManager()
    >>> manager.rank
    0
    >>> manager.world_size
    1
    """

    _state = {}

    def __new__(cls):
        """
        Creates a new instance of the DistributedManager class.

        Parameters
        ----------
        None
        Returns
        -------
        instance : DistributedManager
            New instance of the DistributedManager class.
        """
        instance = super().__new__(cls)
        instance.__dict__ = cls._state
        # default the properties so that the default initialize() could work
        if not hasattr(instance, "_initialized"):
            instance._initialized = False
        if not hasattr(instance, "_has_dist"):
            instance._has_dist = False
        if not hasattr(instance, "_rank"):
            instance._rank = 0
        if not hasattr(instance, "_world_size"):
            instance._world_size = 1
        if not hasattr(instance, "_local_rank"):
            instance._local_rank = 0
        if not hasattr(instance, "_device"):
            instance._device = torch.device("cpu")
        if not hasattr(instance, "_backend"):
            instance._backend = None
        if not hasattr(instance, "_device_mesh"):
            instance._device_mesh = None
        if not hasattr(instance, "_layout_device_mesh"):
            instance._layout_device_mesh = None
        if not hasattr(instance, "_has_subgroups"):
            instance._has_subgroups = False
        if not hasattr(instance, "_device_mesh_subgroups"):
            instance._device_mesh_subgroups = None
        if not hasattr(instance, "_layout_device_mesh_subgroups"):
            instance._layout_device_mesh_subgroups = None
        if not hasattr(instance, "_group"):
            instance._group = {}
        if not hasattr(instance, "_group_rank"):
            instance._group_rank = {}
        if not hasattr(instance, "_group_ranks"):
            instance._group_ranks = {}
        if not hasattr(instance, "_subgroups"):
            instance._subgroups = {}
        if not hasattr(instance, "_subgroups_rank"):
            instance._subgroups_rank = {}
        if not hasattr(instance, "_subgroups_ranks"):
            instance._subgroups_ranks = {}
        if not hasattr(instance, "_layout_subgroups"):
            instance._layout_subgroups = {}
        if not hasattr(instance, "_method_init"):
            instance._method_init = None
        return instance

    @classmethod
    def methods_init_available(cls) -> set[str]:
        """
        Returns the set of available initialization methods for the DistributedManager.

        Returns
        -------
        methods : set[str]
            Set of available initialization methods.
        """
        return {"ENV", "SLURM"}

    @classmethod
    def backend_for_device(cls) -> Dict[str, Optional[str]]:
        """
        Returns the mapping of device types to their default backend.

        Returns
        -------
        backend_for_device : dict[str, str or None]
            Mapping of device types to their default backend.
        """
        backend_for_device = {
            "cuda": "nccl" if torch.distributed.is_nccl_available() else None,
            "cpu": "gloo" if torch.distributed.is_gloo_available() else None,
        }
        return backend_for_device

    @classmethod
    def is_initialized(cls) -> bool:
        """
        Checks if the DistributedManager singleton has been initialized.

        Parameters
        ----------
        None

        Returns
        -------
        initialized : bool
            True if the DistributedManager singleton has been initialized, False otherwise.
        """
        return cls._state.get("_initialized", False)

    def __init__(self):
        """
        Initializes the DistributedManager instance.

        Parameters
        ----------
        None

        Raises
        ------
        RuntimeError
            If the DistributedManager instance is being instantiated before the singleton class is initialized.
        """
        if not self._initialized:
            raise RuntimeError(
                "A DistributedManager instance is being instantiated before "
                "the singleton class is initialized, which can lead to communication "
                "failure among processes. Please call DistributedManager.initialize() "
                "before instantiating any `DistributedManager` instance. "
            )
        super().__init__()

    def __getattr__(self, name: str) -> Any:
        """
        Gets the attribute of the DistributedManager instance.

        Parameters
        ----------
        name : str
            Name of the attribute to get.

        Returns
        -------
        attribute : Any
            Attribute of the DistributedManager instance.

        Raises
        ------
        AttributeError
            If the attribute does not exist.
        """
        # to enable read-only access to the shared _state data
        key_state = f"_{name}"
        has_key_shared_state = key_state in self.__dict__
        has_key = name in self.__dict__
        if has_key_shared_state:
            return self.__dict__[key_state]
        elif has_key:
            return self.__dict__[name]
        else:
            raise AttributeError(f'Attribute "{name}" or "_{name}" not found.')

    def __str__(self):
        """
        Returns a string representation of the DistributedManager instance.

        Parameters
        ----------
        None

        Returns
        -------
        str : str
            String representation of the DistributedManager instance.
        """
        output = (
            f"Initialized process {self.rank} of {self.world_size} using "
            f"method '{self.method_init}'. Device set to {str(self.device)}. Backend is {self.backend}"
        )
        return output

    @staticmethod
    def _setup(
        grid_group_sizes: Optional[_GridGroupSizesType] = None,
        device_type: str = "cuda",
        backend: Optional[str] = None,
        rank: int = -1,
        node_rank: int = -1,
        world_size: int = -1,
        local_rank: Optional[int] = None,
        addr: str = "localhost",
        port: str = "29500",
        method_init: str = "ENV",
        **kwargs_init_pg,
    ):
        """
        Sets up the DistributedManager instance.

        Parameters
        ----------
        grid_group_sizes : OrderedDict, optional
            OrderedDict of group sizes used for distributed communication.
            See create_grid_group() for details and examples.
        device_type : str, optional
            Type of device used for distributed communication.
        backend : str, optional
            Backend used for distributed communication.
        rank : int, optional
            Rank of the process.
        node_rank : int, optional
            Node rank of the process.
        world_size : int, optional
            Total number of processes.
        local_rank : int, optional
            Local rank of the process.
        addr : str, optional
            Address used for distributed communication.
        port : str, optional
            Port used for distributed communication.
        method_init : str, optional
            Method used to initialize the distributed manager.
        kwargs_init_pg:
            kwargs to forward to torch.distributed.init_process_group call

        Returns
        -------
        None
        """
        # TODO: could relax this to allow, e.g., "cuda" for "gloo"
        if device_type == "cuda" and not torch.cuda.is_available():
            raise RuntimeError(f"Input device type {device_type} but torch.cuda is not available")

        if world_size != -1 and grid_group_sizes is not None:
            total_size = 1
            assert hasattr(grid_group_sizes, "values")
            for value in grid_group_sizes.values():
                if isinstance(value, tuple) and all(isinstance(v, int) for v in value):
                    total_size *= prod(value)
                elif isinstance(value, int):
                    total_size *= value
                else:
                    raise RuntimeError(
                        f"Values in grid_group_sizes must be either int or tuple[int, ...], got {type(value)}"
                    )

            if world_size != total_size:
                raise RuntimeError(
                    f"Non-default world_size {world_size} != product of grid_group_sizes values ({total_size})"
                )

        backend_for_device = DistributedManager.backend_for_device()

        if backend_for_device["cpu"] is None and backend_for_device["cuda"] is None:
            raise RuntimeError(f"No backend available for the supported device types: {backend_for_device.keys()}")

        if device_type not in backend_for_device:
            raise RuntimeError(f"Invalid input device type {device_type}: only supports {backend_for_device.keys()}")

        if backend is None:
            backend = backend_for_device[device_type]
        elif backend != backend_for_device[device_type]:
            raise RuntimeError(f"Invalid input backend {backend} for input device type {device_type}")

        # set these in order to call torch.distributed.init_process_group
        os.environ["MASTER_ADDR"] = addr
        os.environ["MASTER_PORT"] = str(port)

        # instantiate the singleton
        DistributedManager._state["_initialized"] = True
        manager = DistributedManager()

        manager._has_dist = torch.distributed.is_available()

        manager._rank = rank
        manager._world_size = world_size
        manager._node_rank = node_rank
        if device_type == "cuda":
            if manager.world_size > torch.cuda.device_count() and manager.world_size % torch.cuda.device_count():
                warn(
                    "world_size is not a multiple of torch.cuda.device_count() so cuda devices could be shared by multiple ranks"
                )
            # will try to guess a local_rank from GPU counts
            if local_rank is None:
                manager._local_rank = manager.rank % torch.cuda.device_count()
            else:
                manager._local_rank = local_rank
            manager._device = torch.device(f"cuda:{manager.local_rank}")
        else:
            if local_rank is not None:
                manager._local_rank = local_rank
            manager._device = torch.device("cpu")

        if not manager.has_dist:
            warn("DistributedManager initialized without torch.distributed package")
            # TODO: triage the importance of having a default device according to the
            # input backend
            return

        if manager.device.type == "cuda":
            # set device before init_process_group to avoid unintended
            # cuda context and to avoid potential NCCL issues
            torch.cuda.set_device(manager.device)
            torch.cuda.device(manager.device)
            torch.cuda.empty_cache()

        manager._backend = backend

        # initialize torch.distributed
        if manager.device.type == "cuda" and backend == "nccl":
            try:
                # to prevent nccl hang and other potential issues:
                # see e.g., https://github.com/pytorch/pytorch/issues/142356
                torch.distributed.init_process_group(
                    manager.backend,
                    rank=manager.rank,
                    world_size=manager.world_size,
                    device_id=manager.device,
                    **kwargs_init_pg,
                )
            except TypeError:
                torch.distributed.init_process_group(
                    manager.backend, rank=manager.rank, world_size=manager.world_size, **kwargs_init_pg
                )
        else:
            torch.distributed.init_process_group(
                manager.backend, rank=manager.rank, world_size=manager.world_size, **kwargs_init_pg
            )

        manager._group["world"] = torch.distributed.group.WORLD
        manager._group_rank["world"] = manager.rank
        manager._group_ranks["world"] = torch.distributed.get_process_group_ranks(manager.group["world"])

        manager._method_init = method_init

        if grid_group_sizes is not None:
            DistributedManager.create_grid_group(grid_group_sizes)

    @staticmethod
    def _create_device_mesh_and_groups(name: list[str], shape: list[int], suffix_mesh: Optional[str] = None) -> None:
        """
        Creates a device mesh and associated process groups for distributed communication.

        Parameters
        ----------
        name : list[str]
            Names of the dimensions in the device mesh.
        shape : list[int]
            Shape of the device mesh, representing the sizes of each dimension.
        suffix_mesh : str, optional
            Suffix to append to the device mesh name for identification, defaults to None.

        Returns
        -------
        None

        Raises
        ------
        RuntimeError
            If DistributedManager is not initialized.
        RuntimeError
            If torch.distributed package is not available.
        RuntimeError
            If method_init is invalid or not available.
        RuntimeError
            If backend is invalid or not available.
        RuntimeError
            If device type is invalid or not available.
        RuntimeError
            If world_size does not match the expected world size computed from shape.
        """
        if not DistributedManager.is_initialized():
            raise RuntimeError("DistributedManager is not initialized upon calling _create_device_mesh_and_groups")
        if not DistributedManager._state["_has_dist"] or not torch.distributed.is_available():
            raise RuntimeError(
                "_create_device_mesh_and_groups requires torch.distributed package, which is not available"
            )
        if (
            DistributedManager._state["_method_init"] is None
            or DistributedManager._state["_method_init"] not in DistributedManager.methods_init_available()
        ):
            raise RuntimeError(
                f"Invalid DistributedManager method_init {DistributedManager._state['_method_init']} "
                "(most likely because it was default initialized)"
            )
        if (
            DistributedManager._state["_backend"] is None
            or DistributedManager._state["_backend"] not in DistributedManager.backend_for_device().values()
        ):
            raise RuntimeError(
                f"Invalid DistributedManager backend {DistributedManager._state['_backend']} "
                "(most likely because it was default initialized)"
            )
        if (
            DistributedManager._state["_device"] is None
            or DistributedManager._state["_device"].type not in DistributedManager.backend_for_device().keys()
        ):
            raise RuntimeError(
                f"Invalid DistributedManager device type {DistributedManager._state['_device'].type} "
                "(most likely because it was default initialized)"
            )

        world_size_expected = prod(shape)

        if world_size_expected != DistributedManager._state["_world_size"]:
            raise RuntimeError(
                f"world_size {DistributedManager._state['_world_size']} does not match the expected world size "
                f"{world_size_expected} computed from the input shape {shape}"
            )

        device_type = DistributedManager._state["_device"].type
        name_mesh = f"_device_mesh_{suffix_mesh}" if suffix_mesh is not None else "_device_mesh"

        # TODO: support arbitrary user-input layout
        layout = LayoutRightMap(tuple(shape))
        DistributedManager._state[f"_layout{name_mesh}"] = layout

        grid2rank = torch.as_strided(torch.arange(world_size_expected), size=layout.shape, stride=layout.strides)
        DistributedManager._state[name_mesh] = torch.distributed.device_mesh.DeviceMesh(
            device_type, grid2rank, mesh_dim_names=tuple(name)
        )

        for i_group in range(len(name)):
            name_group = name[i_group]

            if name_group in DistributedManager._state["_group"]:
                # skip those already created, e.g., from another call of this function
                continue

            DistributedManager._state["_group"][name_group] = DistributedManager._state[name_mesh].get_group(name_group)
            DistributedManager._state["_group_rank"][name_group] = torch.distributed.get_group_rank(
                DistributedManager._state["_group"][name_group], DistributedManager._state["_rank"]
            )
            DistributedManager._state["_group_ranks"][name_group] = torch.distributed.get_process_group_ranks(
                DistributedManager._state["_group"][name_group]
            )

    @staticmethod
    def create_grid_group(grid_group_sizes: _GridGroupSizesType) -> None:
        """
        Creates a grid group for distributed communication.

        Parameters
        ----------
        grid_group_sizes : OrderedDict[str, int | tuple[int, ...]]
            Dictionary of group sizes used for distributed communication. The keys of the OrderedDict
            are the group names and the values are the group sizes. The group sizes can be an integer
            or a tuple of integers. If it is a tuple of integers, they partition the ranks of that group
            into a subgrid of the corresponding shape. The layout of the ranks in the groups and subgroups
            follows the LayoutRightMap convention, where the last group's (or its last subgroup) ranks are
            contiguous global rank on the device grid.

        Returns
        -------
        None

        Notes
        -----
        This method updates the following dictionaries in the DistributedManager._state:
        - _group: Maps group names to torch.distributed.ProcessGroup objects
        - _group_rank: Maps group names to the rank of the current process in that group
        - _group_ranks: Maps group names to lists of all ranks in that group
        - _subgroups: Maps parent group names to lists of their subgroups' ProcessGroup objects
        - _subgroups_rank: Maps parent group names to lists of ranks within each subgroup
        - _subgroups_ranks: Maps parent group names to lists of all ranks for each subgroup
        - _layout_subgroups: Maps parent group names to LayoutMap objects for subgroup coordinate mapping
        - _device_mesh: The main PyTorch DeviceMesh object created from shape_groups
        - _device_mesh_subgroups: The DeviceMesh object for subgroups created from shape_subgroups

        This method creates PyTorch DeviceMesh objects to facilitate distributed tensor operations
        and communication. The DeviceMesh provides a logical view of the physical device layout,
        enabling efficient collective operations across process groups and subgroups.

        Examples
        --------
        >>> # Create a grid with data parallel (dp) dimension of 1 and
        >>> # a 2x2 communication parallel (cp) grid
        >>> from collections import OrderedDict
        >>> grid_group_sizes = OrderedDict([("dp", 1), ("cp", (2, 2))])
        >>> DistributedManager.initialize(grid_group_sizes, device_type="cuda")
        >>> manager = DistributedManager()
        >>> # After initialization, the following dictionaries are populated:
        >>> # manager.group contains ProcessGroup objects for 'world', 'dp', 'cp'
        >>> # and also 'cp_axis_0' and 'cp_axis_1' for the cp subgroups
        >>> # manager.group_rank contains the current process's rank in each group
        >>> # manager.group_ranks contains all ranks for each group
        >>>
        >>> # Device mesh objects are created and accessible
        >>> device_mesh = manager._device_mesh  # DeviceMesh for parent groups ['dp', 'cp']
        >>> device_mesh_subgroups = manager._device_mesh_subgroups  # DeviceMesh for ['dp', 'cp_axis_0', 'cp_axis_1']
        >>> print(f"Parent device mesh shape: {device_mesh.shape()}")  # Expected: (1, 4)
        >>> print(f"Subgroups device mesh shape: {device_mesh_subgroups.shape()}")  # Expected: (1, 2, 2)
        >>>
        >>> # Verify that the subgroups are accessible both directly and via the parent group
        >>> assert "cp_axis_0" in manager.group
        >>> assert "cp_axis_1" in manager.group
        >>> assert manager.subgroups["cp"][0] is manager.group["cp_axis_0"]
        >>> assert manager.subgroups["cp"][1] is manager.group["cp_axis_1"]
        >>>
        >>> # The subgroups represent slices along specific axes of the mesh
        >>> # For a rank with coordinates (i, j) in a 2x2 grid:
        >>> # "cp_axis_0" contains all ranks (0:2, j) - varying along the first axis, fixed j
        >>> # "cp_axis_1" contains all ranks (i, 0:2) - fixed i, varying along the second axis
        >>> # This is consistent with torch.distributed.DeviceMesh's group definition
        >>> # For rank 1 (coords: (0, 1)):
        >>> #   - It belongs to cp_axis_0 group containing ranks [1, 3] (all ranks with j=1)
        >>> #   - It belongs to cp_axis_1 group containing ranks [0, 1] (all ranks with i=0)
        >>>
        >>> # Examine the layout of layout_subgroups['cp']
        >>> layout = manager.layout_subgroups['cp']
        >>> print(f"Shape: {layout.shape}")  # Expected: (2, 2)
        >>> print(f"Strides: {layout.strides}")  # Expected: (2, 1) for LayoutRightMap
        >>>
        >>> # Demonstration of how the layout maps coordinates to ranks
        >>> for i in range(layout.shape[0]):
        ...     for j in range(layout.shape[1]):
        ...         rank = layout.ravel((i, j))
        ...         coords = layout.unravel(rank)
        ...         print(f"Coordinates ({i}, {j}) map to rank {rank}, which maps back to {coords}")
        ...
        Coordinates (0, 0) map to rank 0, which maps back to (0, 0)
        Coordinates (0, 1) map to rank 1, which maps back to (0, 1)
        Coordinates (1, 0) map to rank 2, which maps back to (1, 0)
        Coordinates (1, 1) map to rank 3, which maps back to (1, 1)

        Raises
        ------
        RuntimeError
            If DistributedManager is not initialized.
        RuntimeError
            If torch.distributed package is not available.
        RuntimeError
            If method_init is invalid or not available.
        RuntimeError
            If backend is invalid or not available.
        RuntimeError
            If device type is invalid or not available.
        RuntimeError
            If values in grid_group_sizes are not int or tuple[int, ...].
        """
        shape_groups = []
        name_groups = []
        shape_subgroups = []
        name_subgroups = []
        group2subgroup = {}
        group2subgroup_axes = {}
        assert hasattr(grid_group_sizes, "items")
        for k, v in grid_group_sizes.items():
            if isinstance(v, tuple) and all(isinstance(v_i, int) for v_i in v):
                # Create a new dimension of the DeviceMesh for each group
                shape_groups.append(prod(v))
                name_groups.append(k)
                # Create a new dimension of the DeviceMesh for each subgroup
                # to allow torch DTensor placement on the subgroups' DeviceMesh,
                # where each subgroup axis is treated as a separate dimension in the mesh
                shape_subgroups.extend(v)
                names_this_subgroup = [f"{k}_axis_{i}" for i in range(len(v))]
                name_subgroups.extend(names_this_subgroup)
                # map each group to its subgroups along each axis
                group2subgroup[k] = names_this_subgroup
                group2subgroup_axes[k] = list(range(len(name_subgroups) - len(v), len(name_subgroups)))
            elif isinstance(v, int):
                shape_groups.append(v)
                name_groups.append(k)
                shape_subgroups.append(v)
                name_subgroups.append(k)
            else:
                raise RuntimeError(f"Values in grid_group_sizes must be either int or tuple[int, ...], got {type(v)}")

        # TODO: might not always need the device_mesh for parent groups
        # but one could just create them via create_group()
        DistributedManager._create_device_mesh_and_groups(name_groups, shape_groups)
        if (name_groups == name_subgroups) != (shape_groups == shape_subgroups):
            raise RuntimeError(
                f"Inconsistent group ({name_groups}, {shape_groups}) and "
                f"subgroup ({name_subgroups}, {shape_subgroups}) settings"
            )

        DistributedManager._state["_has_subgroups"] = name_groups != name_subgroups
        if DistributedManager._state["_has_subgroups"]:
            if len(group2subgroup) == 0:
                raise RuntimeError("group2subgroup is empty while _has_subgroups is True")
            DistributedManager._create_device_mesh_and_groups(name_subgroups, shape_subgroups, suffix_mesh="subgroups")
            layout = DistributedManager._state["_layout_device_mesh_subgroups"]
            coords = DistributedManager._state["_device_mesh_subgroups"].get_coordinate()
            for name_group, name_subgroups in group2subgroup.items():
                # map the parent process group name to the subgroups
                DistributedManager._state["_subgroups"][name_group] = [
                    DistributedManager._state["_group"][name_subgroup] for name_subgroup in name_subgroups
                ]
                DistributedManager._state["_subgroups_ranks"][name_group] = [
                    DistributedManager._state["_group_ranks"][name_subgroup] for name_subgroup in name_subgroups
                ]
                DistributedManager._state["_subgroups_rank"][name_group] = [
                    DistributedManager._state["_group_rank"][name_subgroup] for name_subgroup in name_subgroups
                ]
                # create the subgroup layout for each parent group
                # TODO: support LayoutMap.reshape to simplify this
                axes_subgroup = group2subgroup_axes[name_group]
                slices = deepcopy(coords)
                for axis in axes_subgroup:
                    slices[axis] = slice(None)
                layout_subgroup = layout[*slices]
                # create a LayoutMap for the subgroups with offset 0 to be used for a bijective mapping
                # between rank within the subgroups and the subgrid
                DistributedManager._state["_layout_subgroups"][name_group] = LayoutMap(
                    layout_subgroup.strides, layout_subgroup.shape, offset=0
                )

    @staticmethod
    def create_group(name: str, ranks: list[int], **kwargs_dist_ng) -> None:
        """
        Creates a new process group for distributed communication.

        Parameters
        ----------
        name : str
            Name of the group.
        ranks : list[int]
            Ranks of the processes in the group.
        **kwargs_dist_ng
            Keyword arguments to pass to torch.distributed.new_group.

        Returns
        -------
        None

        Notes
        -----
        This method creates a new process group with the given name and ranks,
        and stores the group, ranks, and group rank in the DistributedManager state.
        """
        DistributedManager._state["_group"][name] = torch.distributed.new_group(ranks=ranks, **kwargs_dist_ng)
        DistributedManager._state["_group_ranks"][name] = ranks
        DistributedManager._state["_group_rank"][name] = torch.distributed.get_group_rank(
            DistributedManager._state["_group"][name], DistributedManager._state["_rank"]
        )

    @staticmethod
    def _initialize_env(*args, **kwargs):
        """
        Initializes the DistributedManager instance using environment variables.

        Parameters
        ----------
        *args : list
            Variable length argument list.
        **kwargs : dict
            Arbitrary keyword arguments.

        Returns
        -------
        None
        """
        if not ("RANK" in os.environ and "WORLD_SIZE" in os.environ):
            raise RuntimeError(
                "environment variable RANK and WORLD_SIZE must be set to initialize "
                "torch.distributed using the env:// method"
            )
        rank = os.environ.get("RANK")
        world_size = os.environ.get("WORLD_SIZE")
        local_rank = os.environ.get("LOCAL_RANK")
        # From LightningEnvironment.node_rank()
        group_rank = os.environ.get("GROUP_RANK", 0)
        node_rank = int(os.environ.get("NODE_RANK", group_rank))
        try:
            rank = int(rank)
            world_size = int(world_size)
            if local_rank is not None:
                local_rank = int(local_rank)
        except TypeError:
            raise RuntimeError(
                "environment variables RANK, LOCAL_RANK and WORLD_SIZE must be specified as integer "
                f"but got rank={rank}, local_rank={local_rank}, world_size={world_size}"
            )

        DistributedManager._setup(
            *args,
            rank=rank,
            node_rank=node_rank,
            world_size=world_size,
            local_rank=local_rank,
            addr=os.environ.get("MASTER_ADDR"),
            port=os.environ.get("MASTER_PORT"),
            method_init="ENV",
            **kwargs,
        )

    @staticmethod
    def _initialize_slurm(*args, **kwargs):
        """
        Initializes the DistributedManager instance using SLURM environment variables.

        Parameters
        ----------
        *args : list
            Variable length argument list.
        **kwargs : dict
            Arbitrary keyword arguments.

        Returns
        -------
        None
        """
        keys = ("SLURM_PROCID", "SLURM_NPROCS", "SLURM_LOCALID", "SLURM_LAUNCH_NODE_IPADDR")
        if not all(k in os.environ for k in keys):
            raise RuntimeError(
                f"environment variables {keys} must be set to initialize torch.distributed using the slurm"
            )
        rank = os.environ.get("SLURM_PROCID")
        node_rank = int(os.environ.get("SLURM_NODEID", 0))
        world_size = os.environ.get("SLURM_NPROCS")
        local_rank = os.environ.get("SLURM_LOCALID")
        addr = os.environ.get("SLURM_LAUNCH_NODE_IPADDR")
        try:
            rank = int(rank)
            world_size = int(world_size)
            if local_rank is not None:
                local_rank = int(local_rank)
        except TypeError:
            raise RuntimeError(
                "environment variables SLURM_{PROCID,NPROCS,LOCALID} must be specified as integer "
                f"but got PROCID={rank}, LOCALID={local_rank}, NPROCS={world_size}"
            )

        DistributedManager._setup(
            *args,
            rank=rank,
            node_rank=node_rank,
            world_size=world_size,
            local_rank=local_rank,
            addr=addr,
            method_init="SLURM",
            **kwargs,
        )

    @staticmethod
    def initialize(
        grid_group_sizes: Optional[OrderedDict[str, int | tuple[int, ...]]] = None,
        device_type: str = "cuda",
        backend: Optional[str] = None,
        **kwargs_init_pg,
    ):
        """
        Initializes the DistributedManager instance.

        Parameters
        ----------
        grid_group_sizes : OrderedDict[str, int | tuple[int, ...]]
            Dictionary of group sizes used for distributed communication. The keys of the OrderedDict
            are the group names and the values are the group sizes. The group sizes can be an integer
            or a tuple of integers. If it is a tuple of integers, they partition the ranks of that group
            into a subgrid of the corresponding shape. The layout of the ranks in the groups and subgroups
            follows the LayoutRightMap convention, where the last group's (or its last subgroup) ranks are
            contiguous global rank on the device grid. See create_grid_group() for details and examples
        device_type : str, optional
            Type of device used for distributed communication.
        backend : str, optional
            Backend used for distributed communication.
        kwargs_init_pg:
            kwargs to forward to torch.distributed.init_process_group call

        Returns
        -------
        None
        """
        if DistributedManager.is_initialized():
            warn("DistributedManager is already initialized. Skip initialize()")
            return
        if backend == "nccl":
            # https://pytorch.org/docs/master/notes/cuda.html#id5
            os.environ["TORCH_NCCL_ASYNC_ERROR_HANDLING"] = "0"
        method_init = os.getenv("BOLTZ_DISTRIBUTED_INIT_METHOD")
        if method_init is not None and method_init not in DistributedManager.methods_init_available():
            raise ValueError(
                f"Unknown value set for BOLTZ_DISTRIBUTED_INIT_METHOD={method_init}. "
                f"Allowed options are one of {DistributedManager.methods_init_available()}"
            )
        if method_init is None:
            try:
                DistributedManager._initialize_env(
                    grid_group_sizes, device_type=device_type, backend=backend, **kwargs_init_pg
                )
            except RuntimeError as except_env:
                try:
                    DistributedManager._initialize_slurm(
                        grid_group_sizes, device_type=device_type, backend=backend, **kwargs_init_pg
                    )
                except RuntimeError as except_slurm:
                    warn(
                        "Could not initialize DistributedManager with either the env:// method nor the slurm method.\n"
                        f"Error from env:// method: {except_env} \n"
                        f"Error from the slurm method: {except_slurm} \n"
                        "Will default initialize DistributedManager"
                    )
                    DistributedManager._state["_initialized"] = True
        elif method_init == "ENV":
            DistributedManager._initialize_env(
                grid_group_sizes, device_type=device_type, backend=backend, **kwargs_init_pg
            )
        elif method_init == "SLURM":
            DistributedManager._initialize_slurm(
                grid_group_sizes, device_type=device_type, backend=backend, **kwargs_init_pg
            )

    @staticmethod
    def cleanup():
        """
        Cleans up the DistributedManager instance.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        if DistributedManager._state.get("_group", {}) != {}:
            if torch.distributed.is_initialized():
                # somewhere else has already called torch.distributed.destroy_process_group()
                # need to skip to avoid double destruction
                if DistributedManager._state["_device"].type == "cuda" and torch.cuda.is_available():
                    torch.distributed.barrier(device_ids=[DistributedManager._state["_local_rank"]])
                else:
                    torch.distributed.barrier()
                torch.distributed.destroy_process_group()
            else:
                # otherwise, just clean up the state
                DistributedManager._state = {}
