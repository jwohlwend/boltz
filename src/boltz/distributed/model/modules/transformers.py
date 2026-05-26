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

"""DTensor-compatible transformer modules for Context Parallelism.

Compatible with both Boltz-1x and Boltz-2 serial modules. Supports both
window-batching (AttentionPairBiasShardwise) and ring attention (AttentionPairBias),
dispatched via the ``ring_comm`` parameter at construction time.
"""

from functools import partial
from typing import Callable, Union

from torch import nn
from torch.distributed.device_mesh import DeviceMesh
from torch.distributed.tensor import DTensor
from torch.nn import Module, ModuleList
from torch.utils.checkpoint import checkpoint

from boltz.distributed.comm import AttentionPairBiasComm
from boltz.distributed.model.layers.attention import AttentionPairBias as AttentionPairBiasRing
from boltz.distributed.model.layers.attention import AttentionPairBiasShardwise
from boltz.distributed.model.layers.cat_and_chunk import shardwise_chunk
from boltz.distributed.model.layers.elementwise_op import ElementwiseOp, elementwise_op
from boltz.distributed.model.layers.flatten_and_unflatten import (
    shardwise_flatten_sharded,
    shardwise_unflatten_sharded,
)
from boltz.distributed.model.layers.layernorm import LayerNormParamsReplicated
from boltz.distributed.model.layers.linear import LinearParamsReplicated
from boltz.distributed.model.layers.sigmoid_gate import sigmoid_gate
from boltz.distributed.model.layers.swiglu import SwiGLU as SwiGLUWithDTensor
from boltz.distributed.model.layers.utils import convert_single_repr_window_batched_query_to_key
from boltz.distributed.model.modules.utils import (
    extract_checkpointing_config,
    get_cpu_offload_context,
    validate_window_batching_parameters,
)
from boltz.model.layers.attention import AttentionPairBias as AttentionPairBiasBoltz1
from boltz.model.modules.transformers import AdaLN as AdaLNBoltz1
from boltz.model.modules.transformers import AtomTransformer as AtomTransformerBoltz1
from boltz.model.modules.transformers import ConditionedTransitionBlock as ConditionedTransitionBlockBoltz1
from boltz.model.modules.transformers import DiffusionTransformer as DiffusionTransformerBoltz1
from boltz.model.modules.transformers import DiffusionTransformerLayer as DiffusionTransformerLayerBoltz1
from boltz.model.modules.transformersv2 import AdaLN as AdaLNBoltz2
from boltz.model.modules.transformersv2 import AtomTransformer as AtomTransformerBoltz2
from boltz.model.modules.transformersv2 import ConditionedTransitionBlock as ConditionedTransitionBlockBoltz2
from boltz.model.modules.transformersv2 import DiffusionTransformer as DiffusionTransformerBoltz2
from boltz.model.modules.transformersv2 import DiffusionTransformerLayer as DiffusionTransformerLayerBoltz2


class AdaLN(Module):
    """Adaptive Layer Normalization for DTensor.

    Compatible with both Boltz-1x and Boltz-2 serial AdaLN modules.

    Both versions have identical child modules: a_norm (LayerNorm),
    s_norm (LayerNorm), s_scale (Linear), s_bias (Linear, no bias).
    """

    def __init__(
        self,
        ada_layer_norm: nn.Module,
        device_mesh: DeviceMesh,
    ):
        """Initialize the DTensor-distributed adaptive layer normalization.

        Parameters
        ----------
        ada_layer_norm : nn.Module
            The serial AdaLN module to be distributed. Accepts both
            boltz.model.modules.transformers.AdaLN (Boltz-1x) and
            boltz.model.modules.transformersv2.AdaLN (Boltz-2).
        device_mesh : DeviceMesh
            The device mesh for distributed tensor operations.

        Raises
        ------
        TypeError
            If ada_layer_norm is not an instance of AdaLNBoltz1 or AdaLNBoltz2.
        """
        super().__init__()

        # (1) Set non-module, non-parameter attributes
        self.device_mesh: DeviceMesh = device_mesh

        # (2) Sanity checks
        if not isinstance(ada_layer_norm, (AdaLNBoltz1, AdaLNBoltz2)):
            raise TypeError(
                ", ".join(
                    [
                        f"Instance {ada_layer_norm} should have type {AdaLNBoltz1} or {AdaLNBoltz2}",
                        f"but instead has type {type(ada_layer_norm)}.",
                    ]
                )
            )
        if not isinstance(self.device_mesh, DeviceMesh):
            raise TypeError(
                ", ".join(
                    [
                        f"Instance device_mesh should have type {DeviceMesh}",
                        f"but instead has type {type(self.device_mesh)}.",
                    ]
                )
            )

        # (3) Initialize child modules explicitly
        self.a_norm = LayerNormParamsReplicated(ada_layer_norm.a_norm, device_mesh=device_mesh)
        self.s_norm = LayerNormParamsReplicated(ada_layer_norm.s_norm, device_mesh=device_mesh)
        self.s_scale = LinearParamsReplicated(layer_local=ada_layer_norm.s_scale, device_mesh=device_mesh)
        self.s_bias = LinearParamsReplicated(layer_local=ada_layer_norm.s_bias, device_mesh=device_mesh)

    def forward(self, a: DTensor, s: DTensor) -> DTensor:
        """Forward pass for the DTensor-distributed adaptive layer normalization.

        All tensors use device mesh (dp, cp_axis_0, cp_axis_1).
        Placements: (Shard(0), Shard(1), Replicate()) — batch over dp, sequence/window
        index over cp_axis_0, features replicated over cp_axis_1.

        Parameters
        ----------
        a : DTensor
            The input tensor, shape (B, N, dim) or (B*M, K, W, dim) for window batching.
            Placements: (Shard(0), Shard(1), Replicate()).
        s : DTensor
            The conditioning tensor, shape (B, N, dim_single_cond) or
            (B*M, K, W, dim_single_cond) for window batching.
            Placements: (Shard(0), Shard(1), Replicate()).

        Returns
        -------
        DTensor
            The output tensor, same shape and placements as a.
        """
        a: DTensor = self.a_norm(a)
        s: DTensor = self.s_norm(s)

        gate_input: DTensor = self.s_scale(s)
        a: DTensor = sigmoid_gate(x=a, g=gate_input)
        b: DTensor = self.s_bias(s)
        c: DTensor = elementwise_op(a, b, op=ElementwiseOp.SUM)

        return c


class ConditionedTransitionBlock(Module):
    """Conditioned Transition Block for DTensor.

    Compatible with both Boltz-1x and Boltz-2 serial ConditionedTransitionBlock modules.

    Both versions have identical child modules: adaln, swish_gate (Sequential of
    LinearNoBias + SwiGLU), a_to_b, b_to_a, output_projection (Sequential of
    Linear + Sigmoid). The Sigmoid is stripped and replaced by sigmoid_gate
    in the forward pass for DTensor compatibility.
    """

    def __init__(
        self,
        conditioned_trans_block: nn.Module,
        device_mesh: DeviceMesh,
    ):
        """Initialize the DTensor-distributed conditioned transition block.

        Parameters
        ----------
        conditioned_trans_block : nn.Module
            The serial ConditionedTransitionBlock module to be distributed.
            Accepts both Boltz-1x and Boltz-2 versions.
        device_mesh : DeviceMesh
            The device mesh for distributed tensor operations.

        Raises
        ------
        TypeError
            If conditioned_trans_block is not a recognized type.
        """
        super().__init__()

        if not isinstance(
            conditioned_trans_block, (ConditionedTransitionBlockBoltz1, ConditionedTransitionBlockBoltz2)
        ):
            raise TypeError(
                ", ".join(
                    [
                        f"Instance {conditioned_trans_block} should have type "
                        f"{ConditionedTransitionBlockBoltz1} or {ConditionedTransitionBlockBoltz2}",
                        f"but instead has type {type(conditioned_trans_block)}.",
                    ]
                )
            )

        self.adaln = AdaLN(
            ada_layer_norm=conditioned_trans_block.adaln,
            device_mesh=device_mesh,
        )
        self.swish_gate = nn.Sequential(
            LinearParamsReplicated(
                layer_local=conditioned_trans_block.swish_gate[0],
                device_mesh=device_mesh,
            ),
            SwiGLUWithDTensor(),
        )

        self.a_to_b = LinearParamsReplicated(
            layer_local=conditioned_trans_block.a_to_b,
            device_mesh=device_mesh,
        )
        self.b_to_a = LinearParamsReplicated(
            layer_local=conditioned_trans_block.b_to_a,
            device_mesh=device_mesh,
        )

        # Strip the sigmoid from output_projection - sigmoid operation is handled
        # via sigmoid_gate in the forward pass for DTensor compatibility.
        # Preserves the parameter initialization from the serial module.
        self.output_projection = nn.Sequential(
            LinearParamsReplicated(
                layer_local=conditioned_trans_block.output_projection[0],
                device_mesh=device_mesh,
            ),
        )

    def forward(
        self,
        a: DTensor,
        s: DTensor,
    ) -> DTensor:
        """Forward pass for the DTensor-distributed conditioned transition block.

        All tensors use placements (Shard(0), Shard(1), Replicate()) on mesh (dp, cp_axis_0, cp_axis_1).

        Parameters
        ----------
        a : DTensor
            The input tensor, shape (B, N, dim) or (B*M, K, W, dim).
            Placements: (Shard(0), Shard(1), Replicate()).
        s : DTensor
            The conditioning tensor, shape (B, N, dim_single_cond) or
            (B*M, K, W, dim_single_cond).
            Placements: (Shard(0), Shard(1), Replicate()).

        Returns
        -------
        DTensor
            The output tensor, same shape and placements as a.
        """
        a: DTensor = self.adaln(a, s)
        c: DTensor = self.swish_gate(a)
        b: DTensor = self.a_to_b(a)
        b: DTensor = elementwise_op(c, b, op=ElementwiseOp.PROD)
        a: DTensor = sigmoid_gate(x=self.b_to_a(b), g=self.output_projection[0](s))
        return a


class DiffusionTransformerLayer(Module):
    """Diffusion Transformer Layer for DTensor.

    Compatible with both Boltz-1x and Boltz-2 serial DiffusionTransformerLayer modules.

    Supports two attention modes, dispatched by ``ring_comm``:
    - **Window-batched** (``ring_comm=None``): Uses ``AttentionPairBiasShardwise``.
      Input z is 5D ``(B, K, W, H, D)`` with ``(S(0), S(1), R)`` placements.
    - **Ring attention** (``ring_comm`` provided): Uses ``AttentionPairBias``.
      Input z is 4D ``(B, N, N, D)`` with ``(S(0), S(1), S(2))`` placements.

    Config flags (apply_initial_norm, compute_pair_bias, use_model_cache) are
    auto-detected from the serial module's AttentionPairBias attributes for both
    attention types.
    """

    def __init__(
        self,
        diff_transformer_layer: nn.Module,
        device_mesh: DeviceMesh,
        ring_comm: AttentionPairBiasComm | None = None,
    ):
        """Initialize the DTensor-distributed diffusion transformer layer.

        Parameters
        ----------
        diff_transformer_layer : nn.Module
            The serial DiffusionTransformerLayer module to be distributed.
            Accepts both Boltz-1x and Boltz-2 versions.
        device_mesh : DeviceMesh
            The device mesh for distributed tensor operations.
        ring_comm : AttentionPairBiasComm or None, optional
            Ring communication object. When provided, uses ring attention
            (AttentionPairBias); when None, uses window-batched attention
            (AttentionPairBiasShardwise). Default None.

        Raises
        ------
        TypeError
            If diff_transformer_layer is not a recognized type.
        """
        super().__init__()

        if not isinstance(diff_transformer_layer, (DiffusionTransformerLayerBoltz1, DiffusionTransformerLayerBoltz2)):
            raise TypeError(
                ", ".join(
                    [
                        f"Instance {diff_transformer_layer} should have type "
                        f"{DiffusionTransformerLayerBoltz1} or {DiffusionTransformerLayerBoltz2}",
                        f"but instead has type {type(diff_transformer_layer)}.",
                    ]
                )
            )

        self.adaln = AdaLN(
            ada_layer_norm=diff_transformer_layer.adaln,
            device_mesh=device_mesh,
        )

        # Auto-detect V1/V2 config flags for the attention module
        serial_attn = diff_transformer_layer.pair_bias_attn
        is_boltz1 = isinstance(serial_attn, AttentionPairBiasBoltz1)
        apply_initial_norm = getattr(serial_attn, "initial_norm", False)
        compute_pair_bias = True if is_boltz1 else getattr(serial_attn, "compute_pair_bias", True)
        use_model_cache = is_boltz1

        if ring_comm is not None:
            # Ring attention (all-to-all) — used for token-level transformer
            self.pair_bias_attn = AttentionPairBiasRing(
                attn_pair_bias=serial_attn,
                device_mesh=device_mesh,
                ring_comm=ring_comm,
                apply_initial_norm=apply_initial_norm,
                compute_pair_bias=compute_pair_bias,
                use_model_cache=use_model_cache,
            )
        else:
            # Window-batched attention — used for atom-level transformer
            self.pair_bias_attn = AttentionPairBiasShardwise(
                attn_pair_bias=serial_attn,
                device_mesh=device_mesh,
                apply_initial_norm=apply_initial_norm,
                compute_pair_bias=compute_pair_bias,
                use_model_cache=use_model_cache,
            )

        # Track attention mode for forward dispatch
        self.use_window_batching = isinstance(self.pair_bias_attn, AttentionPairBiasShardwise)

        # In DiffusionTransformerLayer, output_projection_linear is a class attribute.
        # output_projection wraps it with Sigmoid, which is replaced by sigmoid_gate
        # in the forward pass for DTensor compatibility.
        self.output_projection_linear = LinearParamsReplicated(
            layer_local=diff_transformer_layer.output_projection_linear,
            device_mesh=device_mesh,
        )
        self.output_projection = nn.Sequential(self.output_projection_linear)

        self.transition = ConditionedTransitionBlock(
            conditioned_trans_block=diff_transformer_layer.transition,
            device_mesh=device_mesh,
        )

        # Handle post_layer_norm (Boltz-2 only)
        self.post_lnorm = None
        if hasattr(diff_transformer_layer, "post_lnorm") and not isinstance(
            diff_transformer_layer.post_lnorm, nn.Identity
        ):
            self.post_lnorm = LayerNormParamsReplicated(
                diff_transformer_layer.post_lnorm,
                device_mesh=device_mesh,
            )

    def forward(
        self,
        a: DTensor,
        s: DTensor,
        z: DTensor,
        mask: Union[DTensor, None] = None,
        to_keys: Union[Callable[[DTensor], DTensor], None] = None,
        multiplicity: int = 1,
        layer_cache: Union[dict[str, dict[str, DTensor]], None] = None,
        pair_mask: Union[DTensor, None] = None,
    ) -> DTensor:
        """Forward pass for the DTensor-distributed diffusion transformer layer.

        Supports two modes:
        - Window-batched: a/s are 4D (B*M, K, W, D), z is 5D (B, K, W, H, D),
          mask is 3D (B, K, W). Uses to_keys for query→key conversion.
        - Ring attention: a/s are 3D (B*M, N, D), z is 4D (B, N, N, D),
          mask is 2D (B, N). Uses multiplicity for batch expansion.

        Parameters
        ----------
        a : DTensor
            The input tensor.
        s : DTensor
            The conditioning tensor.
        z : DTensor
            The pair representation / pre-computed bias tensor.
        mask : DTensor or None, optional
            The mask tensor.
        to_keys : Callable or None, optional
            Function to transform tensors from query space to key space.
            Used by AttentionPairBiasShardwise for window batching.
        multiplicity : int, optional
            The multiplicity (number of diffusion samples), by default 1.
        layer_cache : dict or None, optional
            Cache for storing projected z during diffusion rollout.
        pair_mask : DTensor or None, optional
            The pair mask tensor.

        Returns
        -------
        DTensor
            The output tensor, same shape and placements as a.
        """
        b: DTensor = self.adaln(a, s)

        if self.use_window_batching:
            if multiplicity != 1:
                raise NotImplementedError(
                    "DiffusionTransformerLayer: window batching mode does not need multiplicity "
                    "but use memory-efficient algorithm to avoid having to explicitly apply multiplicity. "
                    "Multiplicity must be 1 in this mode."
                )
            # Window-batched attention: to_keys converts query→key space
            b: DTensor = self.pair_bias_attn(
                s=b,
                z=z,
                mask=mask,
                to_keys=to_keys,
                model_cache=layer_cache,
            )
        else:
            # Ring attention: uses multiplicity and pair_mask
            b: DTensor = self.pair_bias_attn(
                s=b,
                z=z,
                mask=mask,
                multiplicity=multiplicity,
                model_cache=layer_cache,
                pair_mask=pair_mask,
            )

        b: DTensor = sigmoid_gate(g=self.output_projection[0](s), x=b)

        # Residual connections
        a: DTensor = elementwise_op(a, b, op=ElementwiseOp.SUM)
        c: DTensor = self.transition(a, s)
        a: DTensor = elementwise_op(a, c, op=ElementwiseOp.SUM)

        # Optional post layer norm (Boltz-2 only)
        if self.post_lnorm is not None:
            a = self.post_lnorm(a)

        return a


class DiffusionTransformer(Module):
    """Multi-layer DiffusionTransformer for DTensor.

    Compatible with both Boltz-1x and Boltz-2 serial DiffusionTransformer modules.

    Key difference: Boltz-2 splits the bias across layers (last dim = num_heads * L),
    while Boltz-1 passes the same z to all layers (each layer projects independently).

    Boltz-2's pair_bias_attn=False is not supported (dead code in serial).
    """

    def __init__(
        self,
        diff_transformer: nn.Module,
        device_mesh: DeviceMesh,
        ring_comm: AttentionPairBiasComm | None = None,
    ):
        """Initialize the DTensor-distributed multi-layer diffusion transformer.

        Parameters
        ----------
        diff_transformer : nn.Module
            The serial DiffusionTransformer module to be distributed.
            Accepts both Boltz-1x and Boltz-2 versions.
        device_mesh : DeviceMesh
            The device mesh for distributed tensor operations.
        ring_comm : AttentionPairBiasComm or None, optional
            Ring communication object. When provided, uses ring attention;
            when None, uses window-batched attention. Default None.

        Raises
        ------
        TypeError
            If diff_transformer is not a recognized type.
        NotImplementedError
            If Boltz-2 serial module has pair_bias_attn=False.
        """
        super().__init__()

        if not isinstance(diff_transformer, (DiffusionTransformerBoltz1, DiffusionTransformerBoltz2)):
            raise TypeError(
                ", ".join(
                    [
                        f"Instance {diff_transformer} should have type "
                        f"{DiffusionTransformerBoltz1} or {DiffusionTransformerBoltz2}",
                        f"but instead has type {type(diff_transformer)}.",
                    ]
                )
            )

        # Boltz-2: raise if pair_bias_attn=False (dead code in serial, not supported)
        if isinstance(diff_transformer, DiffusionTransformerBoltz2):
            if not getattr(diff_transformer, "pair_bias_attn", True):
                raise NotImplementedError(
                    "DTensor DiffusionTransformer does not support pair_bias_attn=False. "
                    "This is dead code in the serial Boltz-2 implementation."
                )

        # Detect Boltz-2 bias-splitting mode:
        # Boltz-2 DiffusionTransformer receives bias with last dim = num_heads * L
        # and splits it across layers. Boltz-1 passes the same z to all layers.
        self.split_bias_across_layers = isinstance(diff_transformer, DiffusionTransformerBoltz2)

        # Track attention mode for forward dispatch
        self.use_window_batching = ring_comm is None

        # Detect activation checkpointing and CPU offloading.
        #
        # Boltz-1x: fairscale checkpoint_wrapper replaces each layer's forward
        # method, so we inspect per-layer.
        #
        # Boltz-2: the parent DiffusionTransformer stores
        # ``activation_checkpointing`` as a module-level attribute and handles
        # it in its own ``forward()``.  The individual layers are plain modules.
        # We check the parent attribute as a fallback when per-layer detection
        # yields False.
        activation_checkpointing = set()
        cpu_offloading = set()
        for serial_layer in diff_transformer.layers:
            has_ckpt, has_offload = extract_checkpointing_config(serial_layer)
            activation_checkpointing.add(has_ckpt)
            cpu_offloading.add(has_offload)

        if len(activation_checkpointing) > 1:
            raise ValueError(
                "All layers must have the same activation checkpointing configuration but got different values: ",
                activation_checkpointing,
            )
        if len(cpu_offloading) > 1:
            raise ValueError(
                "All layers must have the same CPU offloading configuration but got different values: ",
                cpu_offloading,
            )

        layer_level_ckpt = activation_checkpointing.pop() if activation_checkpointing else False
        parent_level_ckpt = getattr(diff_transformer, "activation_checkpointing", False)
        self.activation_checkpointing = layer_level_ckpt or parent_level_ckpt
        self.cpu_offloading = cpu_offloading.pop() if cpu_offloading else False

        self.layers = ModuleList(
            [
                DiffusionTransformerLayer(
                    diff_transformer_layer=layer,
                    device_mesh=device_mesh,
                    ring_comm=ring_comm,
                )
                for layer in diff_transformer.layers
            ]
        )

    def forward(
        self,
        a: DTensor,
        s: DTensor,
        z: DTensor,
        mask: Union[DTensor, None] = None,
        to_keys: Union[Callable[[DTensor], DTensor], None] = None,
        multiplicity: int = 1,
        model_cache: Union[dict[str, dict[str, DTensor]], None] = None,
        pair_mask: Union[DTensor, None] = None,
    ) -> DTensor:
        """Forward pass for the DTensor-distributed multi-layer diffusion transformer.

        Supports two modes:
        - Window-batched: a/s 4D, z 5D, mask 3D, uses to_keys.
        - Ring attention: a/s 3D, z 4D, mask 2D, uses multiplicity.

        Parameters
        ----------
        a : DTensor
            The input tensor.
        s : DTensor
            The conditioning tensor.
        z : DTensor
            The pair representation / pre-computed bias tensor.
        mask : DTensor or None, optional
            The mask tensor.
        to_keys : Callable or None, optional
            Function to transform tensors from query space to key space (window-batched).
        multiplicity : int, optional
            The multiplicity (number of diffusion samples), by default 1.
        model_cache : dict or None, optional
            Cache for storing projected z during diffusion rollout.
        pair_mask : DTensor or None, optional
            The pair mask tensor (ring attention only).

        Returns
        -------
        DTensor
            The output tensor, same shape and placements as a.
        """
        if self.split_bias_across_layers and len(self.layers) > 1:
            # Boltz-2: split z last dim across layers
            L = len(self.layers)
            if z.shape[-1] % L != 0:
                raise ValueError(
                    f"Boltz-2 bias last dimension ({z.shape[-1]}) must be evenly divisible by "
                    f"the number of layers ({L}). The Boltz-2 architecture guarantees this because "
                    f"the upstream bias construction (DiffusionConditioning and InputEmbedder) "
                    f"produces z.shape[-1] = num_heads * depth by design."
                )
            # Window-batched: z is 5D (B, K, W, H, heads*L) → L chunks of (B, K, W, H, heads)
            # Ring attention: z is 4D (B, N, N, heads*L) → L chunks of (B, N, N, heads)
            # Both modes use the same shardwise_chunk operation.
            z_chunks = shardwise_chunk(z, chunks=L, dim=-1)
        else:
            z_chunks = None  # Boltz-1: same z for all layers, or single layer

        for i, layer in enumerate(self.layers):
            layer_cache = None
            if model_cache is not None:
                prefix_cache = "layer_" + str(i)
                if prefix_cache not in model_cache:
                    model_cache[prefix_cache] = {}
                layer_cache = model_cache[prefix_cache]

            z_i = z_chunks[i] if z_chunks is not None else z

            if self.activation_checkpointing and self.training:
                if self.cpu_offloading:
                    with get_cpu_offload_context(optimized=True):
                        a = checkpoint(
                            layer,
                            a,
                            s,
                            z_i,
                            mask,
                            to_keys,
                            multiplicity,
                            layer_cache,
                            pair_mask,
                            use_reentrant=False,
                        )
                else:
                    a = checkpoint(
                        layer,
                        a,
                        s,
                        z_i,
                        mask,
                        to_keys,
                        multiplicity,
                        layer_cache,
                        pair_mask,
                        use_reentrant=False,
                    )
            else:
                a = layer(
                    a,
                    s,
                    z_i,
                    mask=mask,
                    to_keys=to_keys,
                    multiplicity=multiplicity,
                    layer_cache=layer_cache,
                    pair_mask=pair_mask,
                )
        return a


class AtomTransformer(Module):
    """AtomTransformer for DTensor (window batching).

    Compatible with both Boltz-1x and Boltz-2 serial AtomTransformer modules.

    Reshapes single repr (B, N, D) -> window-batched (B, K, W, D) using
    shardwise_unflatten_sharded, delegates to DiffusionTransformer, then
    flattens back using shardwise_flatten_sharded.

    Unlike the serial version which flattens (B, K) -> (B*K), the DTensor version
    keeps B and K as separate axes since both are sharded on the device mesh.
    """

    def __init__(
        self,
        layer: nn.Module,
        device_mesh: DeviceMesh,
    ):
        """Initialize the DTensor-distributed atom transformer.

        Parameters
        ----------
        layer : nn.Module
            The serial AtomTransformer module to be distributed.
            Accepts both Boltz-1x and Boltz-2 versions.
        device_mesh : DeviceMesh
            The device mesh for distributed tensor operations.

        Raises
        ------
        TypeError
            If layer is not a recognized type.
        """
        super().__init__()

        if not isinstance(layer, (AtomTransformerBoltz1, AtomTransformerBoltz2)):
            raise TypeError(
                ", ".join(
                    [
                        f"Instance {layer} should have type {AtomTransformerBoltz1} or {AtomTransformerBoltz2}",
                        f"but instead has type {type(layer)}.",
                    ]
                )
            )

        validate_window_batching_parameters(layer.attn_window_queries, layer.attn_window_keys, use_window_batching=True)

        self.attn_window_queries = layer.attn_window_queries
        self.attn_window_keys = layer.attn_window_keys
        self.diffusion_transformer = DiffusionTransformer(
            diff_transformer=layer.diffusion_transformer,
            device_mesh=device_mesh,
        )

    def forward(
        self,
        q: DTensor,
        c: DTensor,
        p: DTensor,
        mask: Union[DTensor, None] = None,
        multiplicity: int = 1,
        model_cache: Union[dict[str, dict[str, DTensor]], None] = None,
        pair_mask: Union[DTensor, None] = None,
        to_keys: None = None,
    ) -> DTensor:
        """Forward pass for the DTensor-distributed atom transformer (window batching).

        All tensors use device mesh (dp, cp_axis_0, cp_axis_1).
        Placements: Shard(0)=dp batch, Shard(1)=cp atom/window axis, Replicate()=cp_axis_1.
        Internally reshapes q/c from (B*M, N, D) to (B*M, K, W, D) via shardwise_unflatten_sharded,
        delegates to DiffusionTransformer, then flattens back.

        Parameters
        ----------
        q : DTensor
            Query single representation, shape (B*M, N, dim) where N = K * W.
            Placements: (Shard(0), Shard(1), Replicate()).
        c : DTensor
            Conditioning single representation, shape (B*M, N, dim_single_cond).
            Placements: (Shard(0), Shard(1), Replicate()).
        p : DTensor
            Pair representation in window-batched format.
            - Boltz-1: shape (B, K, W, H, c_z)
            - Boltz-2: shape (B, K, W, H, num_heads * depth)
            Placements: (Shard(0), Shard(1), Replicate()).
        mask : DTensor or None, optional
            The mask tensor, shape (B, N) or (B*M, N).
            Placements: (Shard(0), Shard(1), Replicate()).
        multiplicity : int, optional
            The multiplicity (number of diffusion samples), by default 1.
            Must be 1 for window batching mode.
        model_cache : dict or None, optional
            Cache for storing projected z during diffusion rollout.
        pair_mask : DTensor or None, optional
            The pair mask tensor. Not supported in window batching mode.
        to_keys : None, optional
            Not used -- to_keys is constructed internally for window batching.

        Returns
        -------
        DTensor
            The output tensor, shape (B*M, N, dim).
            Placements: (Shard(0), Shard(1), Replicate()).
        """
        W = self.attn_window_queries
        H = self.attn_window_keys

        if pair_mask is not None:
            raise NotImplementedError("pair_mask is not supported in AtomTransformer window batching mode")

        if multiplicity != 1:
            raise NotImplementedError(
                "AtomTransformer window batching mode uses memory-efficient algorithm "
                "to avoid having to explicitly apply multiplicity. Multiplicity must be 1 in this mode."
            )

        if q.shape[1] % W != 0:
            raise ValueError(f"q.shape[1] must be divisible by W, but got q.shape[1]={q.shape[1]} and W={W}")

        if c.shape[1] != q.shape[1]:
            raise ValueError(
                f"c.shape[1] must be equal to q.shape[1], but got c.shape[1]={c.shape[1]} and q.shape[1]={q.shape[1]}"
            )

        if mask is not None and mask.shape[1] != q.shape[1]:
            raise ValueError(
                f"mask.shape[1] must be equal to q.shape[1], "
                f"but got mask.shape[1]={mask.shape[1]} and q.shape[1]={q.shape[1]}"
            )

        B, N, D = q.shape
        K = N // W

        # NOTE: p is already in shape (B, K, W, H, D_z)
        if p.ndim != 5:
            raise ValueError(f"p must have 5 dimensions, but got p.ndim={p.ndim}")

        if p.shape[1:-1] != (K, W, H):
            raise ValueError(f"p.shape[1:-1] must be (K, W, H) = {(K, W, H)}, but got p.shape[1:-1]={p.shape[1:-1]}")

        if B % p.shape[0] != 0:
            raise ValueError(f"B must be divisible by p.shape[0], but got B={B} and p.shape[0]={p.shape[0]}")

        # Reshape the single repr into window-batched query view:
        # (B, N, D) -> (B, K, W, D)
        # Unlike the serial version, we don't flatten the resulting (B, K) axes
        # since both of them are sharded on the device mesh.
        q = shardwise_unflatten_sharded(q, axis=1, sizes=(K, W))
        c = shardwise_unflatten_sharded(c, axis=1, sizes=(K, W))
        if mask is not None:
            mask = shardwise_unflatten_sharded(mask, axis=1, sizes=(K, W))

        to_keys_new = partial(
            convert_single_repr_window_batched_query_to_key, W=self.attn_window_queries, H=self.attn_window_keys
        )

        # Main transformer
        q = self.diffusion_transformer(
            q,
            c,
            p,
            to_keys=to_keys_new,
            mask=mask,
            multiplicity=multiplicity,
            model_cache=model_cache,
            pair_mask=pair_mask,
        )

        # Flatten the window-batched query view back to the original single repr view:
        # (B, K, W, D) -> (B, N, D)
        q = shardwise_flatten_sharded(q, start_dim=1, end_dim=2)

        return q
