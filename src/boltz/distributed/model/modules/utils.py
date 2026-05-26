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

"""Distributed model module utilities.

This module provides:

- precision-related helpers for DTensor-based distributed training/inference
- DTensor checkpoint conversion helpers used by context-parallel strategy code
"""

import os
from contextlib import contextmanager
from enum import Enum
from typing import Any, Mapping, Optional

import torch
import torch.distributed as dist
from torch import Tensor
from torch.distributed.device_mesh import DeviceMesh
from torch.distributed.tensor import DTensor, Replicate, Shard, distribute_tensor

from boltz.distributed.model.layers.replicate_op import ReplicateOp, replicate_op
from boltz.distributed.utils import (
    LayoutRightMap,
    all_reduce_weighted_mean,
    create_and_broadcast_tensor_into_placements,
    create_distributed_randn,
)
from boltz.model.modules.utils import random_rotations


def get_cpu_offload_hooks(optimized: bool = True):
    """Return hooks for moving tensors to CPU asynchronously during activation checkpointing.

    Handles both regular ``torch.Tensor`` and ``DTensor`` by offloading the
    underlying local shards.  When *optimized* is True, a dedicated CUDA stream
    and pinned memory are used for true asynchronous transfers.

    Parameters
    ----------
    optimized : bool, optional
        Use an optimised async implementation with a dedicated CUDA stream and
        pinned memory.  Defaults to True.

    Returns
    -------
    tuple[Callable, Callable]
        A ``(pack_hook, unpack_hook)`` pair for use with
        ``torch.autograd.graph.saved_tensors_hooks``.
    """
    offload_stream = torch.cuda.Stream() if optimized else None

    def pack_hook(tensor: Tensor):
        orig_cls = tensor.__class__
        is_dtensor = isinstance(tensor, DTensor)

        if is_dtensor:
            local_tensor = tensor.to_local()
            metadata = (tensor.device_mesh, tensor.placements, tensor.shape, tensor.stride())
        else:
            local_tensor = tensor
            metadata = None

        if local_tensor.is_cuda:
            if optimized:
                with torch.cuda.stream(offload_stream):
                    cpu_tensor = torch.empty(
                        local_tensor.shape, dtype=local_tensor.dtype, device="cpu", pin_memory=True
                    )
                    cpu_tensor.copy_(local_tensor, non_blocking=True)
                return (cpu_tensor, orig_cls, metadata)
            else:
                return (local_tensor.to("cpu", non_blocking=True), orig_cls, metadata)

        return tensor

    def unpack_hook(pack_data):
        if not isinstance(pack_data, tuple):
            return pack_data

        cpu_tensor, cls, metadata = pack_data

        if optimized:
            with torch.cuda.stream(offload_stream):
                gpu_tensor = cpu_tensor.to("cuda", non_blocking=True)
            torch.cuda.current_stream().wait_stream(offload_stream)
        else:
            gpu_tensor = cpu_tensor.to("cuda", non_blocking=True)

        if cls is DTensor:
            device_mesh, placements, shape, stride = metadata
            return DTensor.from_local(gpu_tensor, device_mesh, placements, shape=shape, stride=stride)

        return gpu_tensor

    return pack_hook, unpack_hook


def get_cpu_offload_context(optimized: bool = True):
    """Return a context manager that offloads checkpoint-boundary tensors to CPU.

    When used together with ``torch.utils.checkpoint.checkpoint``, saved
    activations (the *plateau*) are moved to CPU inside the context and
    transparently restored to GPU on the backward pass.

    Parameters
    ----------
    optimized : bool, optional
        Use the optimised async offloading path.  Defaults to True.

    Returns
    -------
    torch.autograd.graph.saved_tensors_hooks
        A context manager wrapping the pack/unpack hooks.
    """
    pack, unpack = get_cpu_offload_hooks(optimized=optimized)
    return torch.autograd.graph.saved_tensors_hooks(pack, unpack)


def extract_checkpointing_config(layer: torch.nn.Module) -> tuple[bool, bool]:
    """Extract activation checkpointing configuration from a single layer.

    Detects if the layer has been wrapped with fairscale's checkpoint_wrapper,
    which replaces the forward method with a functools.partial object.

    Parameters
    ----------
    layer : nn.Module
        A single layer module that may have been wrapped with checkpoint_wrapper.

    Returns
    -------
    tuple[bool, bool]
        (activation_checkpointing, cpu_offloading):
        - activation_checkpointing: True if the layer has checkpointing enabled
        - cpu_offloading: True if checkpointing is configured to offload to CPU

    """
    import functools

    forward_func = getattr(layer.forward, "func", None)
    if (
        isinstance(layer.forward, functools.partial)
        and forward_func is not None
        and forward_func.__name__ == "_checkpointed_forward"
    ):
        cpu_offloading = layer.forward.args[-1]
        return True, cpu_offloading

    return False, False


def has_dtensors(obj: Any) -> bool:
    """Recursively check whether an object contains any DTensors.

    Args:
        obj: Value to inspect. Supported container recursion includes dict/list/tuple.

    Returns:
        ``True`` when at least one DTensor is present, otherwise ``False``.
    """
    if isinstance(obj, DTensor):
        return True
    if isinstance(obj, dict):
        return any(has_dtensors(value) for value in obj.values())
    if isinstance(obj, (list, tuple)):
        return any(has_dtensors(value) for value in obj)
    return False


def convert_dtensors_to_tensors(obj: Any) -> Any:
    """Recursively convert DTensors to plain tensors.

    For ``Replicate``-only placements, ``to_local()`` returns the full
    global tensor with no communication. For any ``Shard``/``Partial``
    placement, this function uses ``full_tensor()`` so checkpoints keep
    global tensor semantics and remain topology-portable.

    Args:
        obj: Value potentially containing DTensors.

    Returns:
        Input structure with all DTensors replaced by plain tensors.
    """
    if isinstance(obj, DTensor):
        if all(isinstance(placement, Replicate) for placement in obj.placements):
            return obj.to_local()
        return obj.full_tensor()
    if isinstance(obj, dict):
        # Keep collective ordering deterministic across ranks when nested
        # sharded DTensors are serialized.
        keys_sorted = sorted(obj.keys(), key=repr)
        return {key: convert_dtensors_to_tensors(obj[key]) for key in keys_sorted}
    if isinstance(obj, list):
        return [convert_dtensors_to_tensors(value) for value in obj]
    if isinstance(obj, tuple):
        return tuple(convert_dtensors_to_tensors(value) for value in obj)
    return obj


def convert_distributed_checkpoint_to_serial_state_dict(checkpoint: Mapping[str, Any]) -> dict[str, Any]:
    """Extract and convert a distributed checkpoint state dict to serial tensors.

    Args:
        checkpoint: Mapping containing at least a ``"state_dict"`` entry.

    Returns:
        A plain ``dict`` where any DTensor entries are converted to plain tensors.

    Raises:
        KeyError: If ``"state_dict"`` is missing from ``checkpoint``.
        TypeError: If ``checkpoint["state_dict"]`` is not mapping-like.
    """
    if "state_dict" not in checkpoint:
        raise KeyError("Checkpoint does not contain 'state_dict'")

    state_dict = checkpoint["state_dict"]
    if not isinstance(state_dict, Mapping):
        raise TypeError("'state_dict' must be a mapping")

    converted = convert_dtensors_to_tensors(state_dict)
    if not isinstance(converted, dict):
        return dict(converted)
    return converted


def _convert_serial_value_to_template_layout(value: Any, template_value: Any) -> Any:
    """Convert one checkpoint value to match a template value layout/device/dtype.

    Handles four tensor-to-tensor cases plus a non-tensor passthrough:

    * **DTensor → DTensor**: validate shape and stride, return as-is.
    * **Tensor → DTensor**: validate shape, distribute to template's mesh/placements.
    * **DTensor → Tensor**: unwrap via ``to_local()``, cast to template device/dtype.
    * **Tensor → Tensor**: cast to template device/dtype.
    * **Non-tensor**: return unchanged.
    """
    # --- Common validation for any tensor-to-tensor conversion ---------------
    both_tensors = isinstance(value, torch.Tensor) and isinstance(template_value, torch.Tensor)
    if both_tensors:
        if tuple(value.shape) != tuple(template_value.shape):
            raise ValueError(
                f"Value shape {tuple(value.shape)} does not match template shape {tuple(template_value.shape)}"
            )
        if tuple(value.stride()) != tuple(template_value.stride()):
            raise ValueError(
                f"Value stride {tuple(value.stride())} does not match template stride {tuple(template_value.stride())}"
            )

    # --- DTensor template ----------------------------------------------------
    if isinstance(template_value, DTensor):
        if isinstance(value, DTensor):
            return value
        if not isinstance(value, torch.Tensor):
            raise TypeError(f"Expected tensor value for DTensor template, got {type(value)}")

        value = value.to(device=template_value.device_mesh.device_type, dtype=template_value.dtype)
        if all(isinstance(p, Replicate) for p in template_value.placements):
            # All ranks load the same checkpoint, so the value is already
            # identical across ranks.  from_local avoids the redundant
            # all-gather that distribute_tensor would trigger.
            return DTensor.from_local(
                value,
                device_mesh=template_value.device_mesh,
                placements=template_value.placements,
                shape=value.shape,
                stride=value.stride(),
            )
        return distribute_tensor(
            value,
            device_mesh=template_value.device_mesh,
            placements=template_value.placements,
        )

    # --- Plain tensor template -----------------------------------------------
    if isinstance(template_value, torch.Tensor):
        if isinstance(value, DTensor):
            value = value.to_local()
        if isinstance(value, torch.Tensor):
            return value.to(device=template_value.device, dtype=template_value.dtype)

    # --- Fallback: unwrap DTensor or pass through ----------------------------
    if isinstance(value, DTensor):
        return value.to_local()
    return value


def convert_serial_checkpoint_to_distributed_state_dict(
    checkpoint: Mapping[str, Any],
    strict: bool = False,
    state_dict_template: Optional[Mapping[str, Any]] = None,
) -> dict[str, Any]:
    """Convert a serial checkpoint state dict to match a distributed state template.

    This helper intentionally works from an explicit ``state_dict_template`` rather
    than constructing a full distributed model, so strategy tests can run without the
    full CP model stack.

    Args:
        checkpoint: Mapping containing a serial ``"state_dict"``.
        strict: Enforce key parity between serial state and template when ``True``.
        state_dict_template: A mapping (typically ``lightning_module.state_dict()``)
            that defines desired output layout/type per key.

    Returns:
        A new state dict aligned to ``state_dict_template``.

    Raises:
        KeyError: If required checkpoint fields are missing, or strict key parity fails.
        TypeError: If ``checkpoint["state_dict"]`` is not mapping-like.
        ValueError: If ``state_dict_template`` is not provided.
    """
    if "state_dict" not in checkpoint:
        raise KeyError("Checkpoint does not contain 'state_dict'")
    if state_dict_template is None:
        raise ValueError("state_dict_template is required to convert serial checkpoint to distributed layout")

    state_dict = checkpoint["state_dict"]
    if not isinstance(state_dict, Mapping):
        raise TypeError("'state_dict' must be a mapping")

    template_keys = set(state_dict_template.keys())
    state_keys = set(state_dict.keys())
    missing_keys = template_keys - state_keys
    extra_keys = state_keys - template_keys
    if strict and (missing_keys or extra_keys):
        msg = "State-dict keys do not match template keys."
        if missing_keys:
            msg += f" Missing keys: {sorted(missing_keys)}."
        if extra_keys:
            msg += f" Extra keys: {sorted(extra_keys)}."
        raise KeyError(msg)

    converted_state: dict[str, Any] = {}
    for key, template_value in state_dict_template.items():
        if key not in state_dict:
            continue
        converted_state[key] = _convert_serial_value_to_template_layout(state_dict[key], template_value)

    if not strict:
        for key in extra_keys:
            converted_state[key] = convert_dtensors_to_tensors(state_dict[key])

    return converted_state


def validate_window_batching_parameters(
    attn_window_queries: Optional[int], attn_window_keys: Optional[int], use_window_batching: bool
) -> None:
    """Validates parameters for window batching in attention mechanisms.

    Args:
        attn_window_queries: Size of the query window. Must be a positive even integer if provided.
        attn_window_keys: Size of the key window. Must be a positive integer if provided.
        use_window_batching: Whether window batching is enabled.

    Raises:
        ValueError: If ``attn_window_queries`` and ``attn_window_keys`` are not both None or both not None.
        ValueError: If ``use_window_batching`` is True but ``attn_window_queries`` is None.
        ValueError: If ``attn_window_queries`` is not a positive even integer.
        ValueError: If ``attn_window_keys`` is not a positive integer.
        ValueError: If ``attn_window_keys`` is not divisible by ``attn_window_queries // 2``.
    """
    if (attn_window_queries is None) != (attn_window_keys is None):
        raise ValueError("attn_window_queries and attn_window_keys must be either both None or both not None")

    if (attn_window_queries is None) == use_window_batching:
        raise ValueError(
            f"attn_window_queries and attn_window_keys must be None if use_window_batching is False, otherwise they must be not None, but got attn_window_queries={attn_window_queries}, attn_window_keys={attn_window_keys} and use_window_batching={use_window_batching}"
        )

    if attn_window_queries is not None:
        if not isinstance(attn_window_queries, int) or attn_window_queries <= 0:
            raise ValueError("attn_window_queries must be a positive integer")

        if attn_window_queries % 2 != 0:
            raise ValueError("attn_window_queries must be even")

    if attn_window_keys is not None:
        if not isinstance(attn_window_keys, int) or attn_window_keys <= 0:
            raise ValueError("attn_window_keys must be a positive integer")

        if attn_window_keys % (attn_window_queries // 2) != 0:
            raise ValueError("attn_window_keys must be divisible by attn_window_queries // 2")


class Precision(Enum):
    """Precision modes for model computation."""

    BF16 = "BF16"
    BF16_MIXED = "BF16_MIXED"
    FP16 = "FP16"
    TF32 = "TF32"
    FP32 = "FP32"
    FP64 = "FP64"


class SDPAWithBiasBackend(Enum):
    """Scaled dot-product attention with bias backend implementations."""

    REFERENCE = "reference"
    TORCH_SDPA_EFFICIENT_ATTENTION = "torch_sdpa_efficient_attention"
    TORCH_FLEX_ATTN = "torch_flex_attn"


class TriAttnBackend(Enum):
    """Triangle attention backend implementations (for distributed triangular attention)."""

    REFERENCE = "reference"
    CUEQ = "cueq"
    TRIFAST = "trifast"
    CUEQ_FWD_TRIFAST_BWD = "cueq_fwd_trifast_bwd"


class SetTriAttnBackend:
    """Callable that sets ``triattn_backend`` on every :class:`PairformerLayer` in a model.

    Designed for use with :meth:`torch.nn.Module.apply`::

        from boltz.distributed.model.modules.utils import SetTriAttnBackend, TriAttnBackend
        model.apply(SetTriAttnBackend(TriAttnBackend.CUEQ))

    ``MSALayer`` is **not** targeted directly because it contains a
    ``PairformerNoSeqLayer`` child which is reached by the recursive
    ``apply`` traversal.
    """

    def __init__(self, triattn_backend: TriAttnBackend) -> None:
        # Lazy import: PairformerLayer imports TriAttnBackend from this module,
        # so a top-level import would create a circular dependency.
        from boltz.distributed.model.layers.pairformer import PairformerLayer

        valid = (
            TriAttnBackend.REFERENCE,
            TriAttnBackend.CUEQ,
            TriAttnBackend.TRIFAST,
            TriAttnBackend.CUEQ_FWD_TRIFAST_BWD,
        )
        if triattn_backend not in valid:
            raise ValueError(f"triattn_backend must be one of {valid} but got {triattn_backend}")
        self.triattn_backend = triattn_backend
        self.supported_module_types = (PairformerLayer,)

    def __call__(self, module: torch.nn.Module) -> None:
        if not isinstance(module, self.supported_module_types):
            return
        if not hasattr(module, "triattn_backend"):
            raise AttributeError(
                f"Module {type(module).__name__} should but does not have a 'triattn_backend' attribute"
            )
        module.triattn_backend = self.triattn_backend


class SetAttnPairBiasBackend:
    """Callable that sets ``sdpa_with_bias_backend`` on every :class:`AttentionPairBias` in a model.

    Designed for use with :meth:`torch.nn.Module.apply`::

        from boltz.distributed.model.modules.utils import SDPAWithBiasBackend, SetAttnPairBiasBackend
        model.apply(SetAttnPairBiasBackend(SDPAWithBiasBackend.TORCH_FLEX_ATTN))

    Only ``REFERENCE`` and ``TORCH_FLEX_ATTN`` are valid for ring-attention
    :class:`AttentionPairBias`; see the validation in its ``__init__``.
    """

    def __init__(self, sdpa_with_bias_backend: SDPAWithBiasBackend) -> None:
        # Lazy import: attention.py imports from this module, so a top-level
        # import would create a circular dependency.
        from boltz.distributed.model.layers.attention import AttentionPairBias

        valid = (SDPAWithBiasBackend.REFERENCE, SDPAWithBiasBackend.TORCH_FLEX_ATTN)
        if sdpa_with_bias_backend not in valid:
            raise ValueError(f"sdpa_with_bias_backend must be one of {valid} but got {sdpa_with_bias_backend}")
        self.sdpa_with_bias_backend = sdpa_with_bias_backend
        self._target_type = AttentionPairBias

    def __call__(self, module: torch.nn.Module) -> None:
        if not isinstance(module, self._target_type):
            return
        if not hasattr(module, "sdpa_with_bias_backend"):
            raise AttributeError(
                f"Module {type(module).__name__} should but does not have a " f"'sdpa_with_bias_backend' attribute"
            )
        module.sdpa_with_bias_backend = self.sdpa_with_bias_backend


class SetAttnPairBiasShardwiseBackend:
    """Callable that sets ``sdpa_with_bias_backend`` on every :class:`AttentionPairBiasShardwise` in a model.

    Designed for use with :meth:`torch.nn.Module.apply`::

        from boltz.distributed.model.modules.utils import SDPAWithBiasBackend, SetAttnPairBiasShardwiseBackend
        model.apply(SetAttnPairBiasShardwiseBackend(SDPAWithBiasBackend.TORCH_SDPA_EFFICIENT_ATTENTION))

    All three ``SDPAWithBiasBackend`` members are valid for window-batched
    :class:`AttentionPairBiasShardwise`.
    """

    def __init__(self, sdpa_with_bias_backend: SDPAWithBiasBackend) -> None:
        # Lazy import: attention.py imports from this module, so a top-level
        # import would create a circular dependency.
        from boltz.distributed.model.layers.attention import AttentionPairBiasShardwise

        valid = (
            SDPAWithBiasBackend.REFERENCE,
            SDPAWithBiasBackend.TORCH_SDPA_EFFICIENT_ATTENTION,
            SDPAWithBiasBackend.TORCH_FLEX_ATTN,
        )
        if sdpa_with_bias_backend not in valid:
            raise ValueError(f"sdpa_with_bias_backend must be one of {valid} but got {sdpa_with_bias_backend}")
        self.sdpa_with_bias_backend = sdpa_with_bias_backend
        self._target_type = AttentionPairBiasShardwise

    def __call__(self, module: torch.nn.Module) -> None:
        if not isinstance(module, self._target_type):
            return
        if not hasattr(module, "sdpa_with_bias_backend"):
            raise AttributeError(
                f"Module {type(module).__name__} should but does not have a " f"'sdpa_with_bias_backend' attribute"
            )
        module.sdpa_with_bias_backend = self.sdpa_with_bias_backend


class OffloadActvCkptToCPU:
    """Callable that enables ``cpu_offloading`` on selected distributed module types.

    Designed for use with :meth:`torch.nn.Module.apply`::

        from boltz.distributed.model.modules.utils import OffloadActvCkptToCPU
        model.apply(OffloadActvCkptToCPU(["DiffusionTransformer", "PairformerModule"]))

    Each targeted module must already have ``activation_checkpointing = True``;
    a :class:`ValueError` is raised otherwise.

    Parameters
    ----------
    module_types : set[str]
        Subset of ``{"DiffusionTransformer", "MSAModule", "PairformerModule"}``.
    """

    def __init__(self, module_types: set[str]) -> None:
        from boltz.distributed.model.layers.pairformer import PairformerModule
        from boltz.distributed.model.modules.transformers import DiffusionTransformer
        from boltz.distributed.model.modules.trunkv2 import MSAModule

        valid_map: dict[str, type] = {
            "DiffusionTransformer": DiffusionTransformer,
            "MSAModule": MSAModule,
            "PairformerModule": PairformerModule,
        }
        module_types = set(module_types)
        invalid = module_types - valid_map.keys()
        if invalid:
            raise ValueError(
                f"Invalid module type(s) {sorted(invalid)} for OffloadActvCkptToCPU. "
                f"Valid types: {sorted(valid_map)}"
            )
        if not module_types:
            raise ValueError("module_types must be non-empty")
        self._target_types = tuple(valid_map[n] for n in sorted(module_types))

    def __call__(self, module: torch.nn.Module) -> None:
        if not isinstance(module, self._target_types):
            return
        for attr in ("activation_checkpointing", "cpu_offloading"):
            if not hasattr(module, attr):
                raise AttributeError(f"Module {type(module).__name__} should but does not have a '{attr}' attribute")
        if not module.activation_checkpointing:
            raise ValueError(
                f"Cannot enable cpu_offloading on {type(module).__name__} because "
                f"activation_checkpointing is not enabled. Enable it first "
                f"(e.g. model.msa_args/pairformer_args/score_model_args"
                f".activation_checkpointing=true)."
            )
        module.cpu_offloading = True


PRECISION_TO_DTYPE = {
    Precision.BF16: torch.bfloat16,
    Precision.BF16_MIXED: torch.bfloat16,
    Precision.FP16: torch.float16,
    Precision.TF32: torch.float32,
    Precision.FP32: torch.float32,
    Precision.FP64: torch.float64,
}


DTYPE_TO_PRECISION = {
    # no BF16-MIXED mapping as it's only relevant for training (mostly)
    # also, this util dict is only used to look up dtype-specific attention INF values
    # in the tests
    torch.bfloat16: Precision.BF16,
    torch.float16: Precision.FP16,
    torch.float32: Precision.FP32,
    torch.float64: Precision.FP64,
}


PRECISION_TO_LIGHTNING = {
    Precision.BF16: "bf16-true",
    Precision.BF16_MIXED: "bf16-mixed",
    Precision.FP16: "fp16-true",
    Precision.TF32: "32",
    Precision.FP32: "32",
    Precision.FP64: "64",
}


@contextmanager
def setup_tf32_env(precision: Precision):
    """Context manager to setup TF32 environment based on precision setting.

    This context manager temporarily modifies TF32 settings for CUDA operations
    and automatically restores the original settings when exiting the context.

    Args:
        precision (Precision): Target precision mode

    Example:
        >>> with setup_tf32_env(Precision.TF32):
        ...     # TF32 is enabled for this block
        ...     result = model(input_tensor)
        >>> # Original TF32 settings are restored

        >>> with setup_tf32_env(Precision.FP32):
        ...     # TF32 is explicitly disabled for pure FP32
        ...     result = model(input_tensor)

    Note:
        This affects both CUDA matrix operations and cuDNN operations.
        The original environment is always restored, even if an exception occurs.
    """
    # Store original TF32 settings
    original_env = os.environ.get("NVIDIA_TF32_OVERRIDE", None)
    original_matmul_tf32 = torch.backends.cuda.matmul.allow_tf32
    original_cudnn_tf32 = torch.backends.cudnn.allow_tf32

    # Setup TF32 environment based on precision
    use_tf32 = precision == Precision.TF32

    if use_tf32:
        os.environ["NVIDIA_TF32_OVERRIDE"] = "1"
        torch.backends.cuda.matmul.allow_tf32 = True
        torch.backends.cudnn.allow_tf32 = True
    elif precision == Precision.FP32:
        # Explicitly disable TF32 for pure FP32
        os.environ["NVIDIA_TF32_OVERRIDE"] = "0"
        torch.backends.cuda.matmul.allow_tf32 = False
        torch.backends.cudnn.allow_tf32 = False
    # For fp64 or other precisions, leave TF32 settings unchanged

    try:
        yield
    finally:
        # Restore original TF32 settings
        if original_env is not None:
            os.environ["NVIDIA_TF32_OVERRIDE"] = original_env
        else:
            os.environ.pop("NVIDIA_TF32_OVERRIDE", None)
        torch.backends.cuda.matmul.allow_tf32 = original_matmul_tf32
        torch.backends.cudnn.allow_tf32 = original_cudnn_tf32


def create_and_broadcast_random_rotation(
    shape: tuple[int, ...],
    device_mesh: DeviceMesh,
    dtype: torch.dtype = torch.float32,
) -> Tensor:
    """Create a distributed random rotation matrix.

    Parameters
    ----------
    shape : tuple[int, ...]
        Shape of the random rotation matrix, e.g. (B, 3, 3).
    device_mesh : DeviceMesh
        The device mesh for DTensor operations.
    dtype : torch.dtype, optional
        The dtype of the random rotation matrix. Defaults to torch.float32.

    Returns
    -------
    Tensor
        The local random rotation matrix tensor (shard of the batch dim).

    """

    def create_rand_rot_fn(shape_local, dtype, device):
        return random_rotations(shape_local[0], dtype=dtype, device=device)

    tensor_local = create_and_broadcast_tensor_into_placements(
        shape=shape,
        create_local_fn=create_rand_rot_fn,
        device_mesh=device_mesh,
        placements=(Shard(0), Replicate(), Replicate()),
        dtype=dtype,
    )
    return tensor_local


@torch.no_grad()
def randomly_rotate(
    coords: DTensor,
    second_coords: Optional[DTensor] = None,
    return_roto: bool = False,
) -> tuple[DTensor, Optional[DTensor], Optional[DTensor]]:
    """Randomly rotate coordinates using DTensor operations. Does not support backward.

    Parameters
    ----------
    coords : DTensor
        The coordinates to rotate, shape (B, N, 3).
        Placements: (Shard(0), Shard(1), Replicate()).
    second_coords : Optional[DTensor], optional
        Optional second coordinates to rotate with the same rotation matrix.
    return_roto : bool, optional
        Whether to return the rotation matrix.

    Returns
    -------
    tuple[DTensor, Optional[DTensor], Optional[DTensor]]
        Rotated coords, rotated second_coords (if provided), rotation matrix (if requested).

    """
    if coords.requires_grad:
        raise ValueError("randomly_rotate does not support backward pass but got coords.requires_grad is True")

    device_mesh = coords.device_mesh
    placements = coords.placements

    if placements != (Shard(0), Shard(1), Replicate()):
        raise ValueError(f"Expected placements (Shard(0), Shard(1), Replicate()), got {placements}")

    if second_coords is not None and second_coords.placements != (Shard(0), Shard(1), Replicate()):
        raise ValueError(
            f"Expected second_coords placements (Shard(0), Shard(1), Replicate()), got {second_coords.placements}"
        )

    size_batch = coords.shape[0]

    # Create random rotation matrix (local shard)
    R_local = create_and_broadcast_random_rotation(
        shape=(size_batch, 3, 3),
        device_mesh=device_mesh,
        dtype=coords.to_local().dtype,
    )

    # Apply rotation using einsum on local tensors
    coords_local = coords.to_local()
    coords_rotated_local = torch.einsum("bmd,bds->bms", coords_local, R_local)

    coords_rotated = DTensor.from_local(
        coords_rotated_local,
        device_mesh=device_mesh,
        placements=placements,
        shape=coords.shape,
        stride=coords.stride(),
    )

    # Handle second_coords
    second_coords_rotated = None
    if second_coords is not None:
        second_coords_local = second_coords.to_local()
        second_coords_rotated_local = torch.einsum("bmd,bds->bms", second_coords_local, R_local)
        second_coords_rotated = DTensor.from_local(
            second_coords_rotated_local,
            device_mesh=device_mesh,
            placements=placements,
            shape=second_coords.shape,
            stride=second_coords.stride(),
        )

    # Return rotation matrix if requested
    roto = None
    if return_roto:
        shape_roto = (size_batch, 3, 3)
        stride_roto = LayoutRightMap(shape_roto).strides
        roto = DTensor.from_local(
            R_local,
            device_mesh,
            (Shard(0), Replicate(), Replicate()),
            shape=shape_roto,
            stride=stride_roto,
        )

    return coords_rotated, second_coords_rotated, roto


def center_random_augmentation(
    atom_coords: DTensor,
    atom_mask: DTensor,
    s_trans: float = 1.0,
    augmentation: bool = True,
    centering: bool = True,
    return_second_coords: bool = False,
    second_coords: Optional[DTensor] = None,
    return_roto: bool = False,
) -> tuple[DTensor, DTensor, DTensor] | tuple[DTensor, DTensor] | DTensor:
    """Center and randomly augment coordinates using DTensor operations. Does not support backward.

    Parameters
    ----------
    atom_coords : DTensor
        Atom coordinates, shape (B, N, 3).
        Placements: (Shard(0), Shard(1), Replicate()).
    atom_mask : DTensor
        Atom mask, shape (B, N).
        Placements: (Shard(0), Shard(1), Replicate()).
    s_trans : float, optional
        Translation scale factor, by default 1.0.
    augmentation : bool, optional
        Whether to add random rotation + translation, by default True.
    centering : bool, optional
        Whether to center coordinates to zero mean, by default True.
    return_second_coords : bool, optional
        Whether to return transformed second coordinates, by default False.
    second_coords : Optional[DTensor], optional
        Second coordinates to apply the same transformation.
    return_roto : bool, optional
        Whether to return the rotation matrix, by default False.

    Returns
    -------
    DTensor | tuple[DTensor, ...]
        Augmented coordinates, and optionally second coords and rotation matrix.

    """
    if atom_coords.requires_grad:
        raise ValueError("center_random_augmentation does not support backward pass")

    if second_coords is not None and second_coords.requires_grad:
        raise ValueError("center_random_augmentation does not support backward pass for second_coords")

    if return_roto and not augmentation:
        raise ValueError("cannot return rotation matrix when augmentation is False")

    device_mesh = atom_coords.device_mesh
    input_placements = atom_coords.placements
    if input_placements != (Shard(0), Shard(1), Replicate()):
        raise ValueError(f"Expected placements (Shard(0), Shard(1), Replicate()), got {input_placements}")

    cp_axis_0_group = device_mesh.get_group("cp_axis_0")
    cp_axis_1_group = device_mesh.get_group("cp_axis_1")
    cp_axis_1_rank = device_mesh.get_local_rank("cp_axis_1")

    if centering:
        # Compute mean on cp_axis_1 rank 0, then broadcast to all column ranks
        if cp_axis_1_rank == 0:
            atom_coords_local = atom_coords.to_local().requires_grad_(False)
            atom_mask_local = atom_mask.to_local().requires_grad_(False)

            atom_mean_local = all_reduce_weighted_mean(
                atom_mask_local.unsqueeze(-1),
                atom_coords_local,
                group_reduce=cp_axis_0_group,
                dim=1,
            )
        else:
            atom_coords_local = atom_coords.to_local()
            atom_mean_local = torch.empty_like(atom_coords_local[:, 0, :])

        # Broadcast mean across cp_axis_1
        dist.broadcast(atom_mean_local, dist.get_global_rank(cp_axis_1_group, 0), cp_axis_1_group)

        shape_atom_mean_global = (atom_coords.shape[0], atom_coords.shape[-1])
        stride_atom_mean_global = LayoutRightMap(shape_atom_mean_global).strides
        atom_mean = DTensor.from_local(
            atom_mean_local,
            device_mesh=device_mesh,
            placements=(Shard(0), Replicate(), Replicate()),
            shape=shape_atom_mean_global,
            stride=stride_atom_mean_global,
        )

        atom_coords = replicate_op(atom_coords, atom_mean, 1, ReplicateOp.SUB)
        if second_coords is not None:
            second_coords = replicate_op(second_coords, atom_mean, 1, ReplicateOp.SUB)

    if augmentation:
        atom_coords, second_coords, roto = randomly_rotate(
            atom_coords,
            second_coords=second_coords,
            return_roto=return_roto,
        )

        # Generate and apply random translation
        batch_size = atom_coords.shape[0]
        random_trans = create_distributed_randn(
            (batch_size, 1, 3),
            device_mesh=device_mesh,
            placements=(Shard(0), Replicate(), Replicate()),
            dtype=atom_coords.dtype,
            scale=s_trans,
        )

        with torch.no_grad():
            atom_coords_local = atom_coords.to_local()
            random_trans_local = random_trans.to_local()
            atom_coords = DTensor.from_local(
                atom_coords_local + random_trans_local,
                device_mesh=device_mesh,
                placements=input_placements,
                shape=atom_coords.shape,
                stride=atom_coords.stride(),
            )

        if second_coords is not None:
            with torch.no_grad():
                second_coords_local = second_coords.to_local()
                second_coords = DTensor.from_local(
                    second_coords_local + random_trans_local,
                    device_mesh=device_mesh,
                    placements=input_placements,
                    shape=second_coords.shape,
                    stride=second_coords.stride(),
                )

    if return_second_coords and return_roto:
        return atom_coords, second_coords, roto
    elif return_second_coords:
        return atom_coords, second_coords
    elif return_roto:
        return atom_coords, roto
    else:
        return atom_coords
