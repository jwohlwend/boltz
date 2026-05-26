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
from typing import Callable, Union

from torch import nn
from torch.distributed.device_mesh import DeviceMesh
from torch.distributed.tensor import DTensor

from boltz.distributed.comm import AttentionPairBiasComm
from boltz.distributed.model.layers.attention_impl import (
    _AttentionPairBiasContextVecParams,
    _AttentionPairBiasContexVecImpl,
    _AttentionPairBiasShardwiseImpl,
)
from boltz.distributed.model.layers.layernorm import LayerNormParamsReplicated
from boltz.distributed.model.layers.linear import LinearParamsReplicated
from boltz.distributed.model.layers.sigmoid_gate import sigmoid_gate
from boltz.distributed.model.modules.utils import SDPAWithBiasBackend
from boltz.model.layers.attention import AttentionPairBias as AttentionPairBiasSerialV1
from boltz.model.layers.attentionv2 import AttentionPairBias as AttentionPairBiasSerialV2


class AttentionPairBias(nn.Module):
    """Attention pair bias module based on DTensor with ring attention.

    This module implements global (non-window-batched) attention with pair bias
    using ring communication patterns for context parallelism.

    The __init__() method follows the pattern of distribute_module(), and
    so takes a device mesh as an argument.  See the following link for details:
    https://docs.pytorch.org/docs/stable/distributed.tensor.html#torch.distributed.tensor.distribute_module

    Configuration Flags
    -------------------
    The module supports both V1 (Boltz-1x) and V2 (Boltz-2) API styles through
    configuration flags:

    - apply_initial_norm: V1=True (has norm_s LayerNorm), V2=False (no initial norm)
    - compute_pair_bias: V1=True (project z via LayerNorm+Linear),
      V2=configurable (False for DiffusionTransformerLayer where z is pre-computed bias).
      Mutually exclusive with use_model_cache when False.
    - use_model_cache: V1=True (cache z projection), V2=False (no caching).
      Only valid when compute_pair_bias=True.

    Use Cases
    ---------
    U1: PairFormerModule (global attention, no window batching)
        - multiplicity=1, compute_pair_bias=True
        - k_in=s (queries equal keys)

    U2: AtomDiffusion with multiplicity (non-window-batched, all-to-all)
        - multiplicity >= 1
        - k_in=s or pre-computed

    Note: Window batching use cases should use AttentionPairBiasShardwise instead.
    """

    def __init__(
        self,
        attn_pair_bias: nn.Module,
        device_mesh: DeviceMesh,
        ring_comm: AttentionPairBiasComm,
        sdpa_with_bias_backend: SDPAWithBiasBackend = SDPAWithBiasBackend.REFERENCE,
        # Configuration flags for V1/V2 API compatibility
        apply_initial_norm: bool = False,  # V1=True, V2=False (default to V2)
        compute_pair_bias: bool = True,  # V1=True, V2=configurable (False for DiffusionTransformerLayer)
        use_model_cache: bool = False,  # V1=True, V2=False (default to V2)
    ) -> None:
        """Initialize the attention pair bias layer.

        Parameters
        ----------
        attn_pair_bias : nn.Module
            The serial attention pair bias layer to convert to DTensor.
        device_mesh : DeviceMesh
            The device mesh for distributed tensor operations.
        ring_comm : AttentionPairBiasComm
            The ring communication object for context parallelism.
        sdpa_with_bias_backend : SDPAWithBiasBackend, optional
            The attention backend to use. Default is REFERENCE.
        apply_initial_norm : bool, optional
            Whether to apply LayerNorm to input s. V1=True, V2=False. Default False.
        compute_pair_bias : bool, optional
            Whether to compute pair bias (LayerNorm + Linear on z). V1=True,
            V2=configurable (False for DiffusionTransformerLayer where z is
            pre-computed bias). Mutually exclusive with use_model_cache=True
            when False. Default True.
        use_model_cache : bool, optional
            Whether to cache z projection for diffusion rollout. V1=True, V2=False.
            Only valid when compute_pair_bias=True. Default False.

        Raises
        ------
        TypeError
            If device_mesh is not a DeviceMesh or ring_comm is not an AttentionPairBiasComm.
        ValueError
            If sdpa_with_bias_backend is not supported, or if use_model_cache=True
            with compute_pair_bias=False.
        """
        super().__init__()

        # (0) Type check on serial module
        if not isinstance(attn_pair_bias, (AttentionPairBiasSerialV1, AttentionPairBiasSerialV2)):
            raise TypeError(
                ", ".join(
                    [
                        f"Instance {attn_pair_bias} should have type "
                        f"{AttentionPairBiasSerialV1} or {AttentionPairBiasSerialV2}",
                        f"but instead has type {type(attn_pair_bias)}.",
                    ]
                )
            )

        # (1) Set non-module, non-parameter attributes from serial module
        self.c_s = attn_pair_bias.c_s
        self.num_heads = attn_pair_bias.num_heads
        self.head_dim = attn_pair_bias.head_dim
        self.inf = attn_pair_bias.inf
        # Mutable backend selection for scaled dot-product attention.
        # Default is REFERENCE. To switch backend for the entire model, use
        # ``model.apply(SetAttnPairBiasBackend(backend))``
        # (see boltz.distributed.model.modules.utils.SetAttnPairBiasBackend).
        self.sdpa_with_bias_backend = (
            sdpa_with_bias_backend
            if isinstance(sdpa_with_bias_backend, SDPAWithBiasBackend)
            else SDPAWithBiasBackend(sdpa_with_bias_backend)
        )
        if self.sdpa_with_bias_backend not in [
            SDPAWithBiasBackend.TORCH_FLEX_ATTN,
            SDPAWithBiasBackend.REFERENCE,
        ]:
            raise ValueError(
                f"Unsupported sdpa_with_bias_backend: {self.sdpa_with_bias_backend}. "
                f"Only TORCH_FLEX_ATTN and REFERENCE are supported."
            )

        # Configuration flags — use_model_cache caches the z projection output,
        # which requires compute_pair_bias=True.
        if use_model_cache and not compute_pair_bias:
            raise ValueError(
                "use_model_cache=True requires compute_pair_bias=True because the cache "
                "stores the z projection output. Got compute_pair_bias=False."
            )
        self.apply_initial_norm = apply_initial_norm
        self.compute_pair_bias = compute_pair_bias
        self.use_model_cache = use_model_cache

        # Ring attention does not support window batching
        self.use_window_batching = False

        self.device_mesh = device_mesh
        self.ring_comm = ring_comm

        # (2) Sanity checks on non-module, non-parameter attributes
        if not isinstance(self.device_mesh, DeviceMesh):
            raise TypeError(f"Input '{device_mesh}' must be of type {DeviceMesh}. Got type {type(self.device_mesh)}.")
        if not isinstance(self.ring_comm, AttentionPairBiasComm):
            raise TypeError(
                f"Input '{ring_comm}' must be of type {AttentionPairBiasComm}. Got type {type(self.ring_comm)}."
            )

        # (3) Initialize child modules explicitly from serial module
        if self.apply_initial_norm:
            self.norm_s = LayerNormParamsReplicated(attn_pair_bias.norm_s, device_mesh=device_mesh)

        self.proj_q = LinearParamsReplicated(layer_local=attn_pair_bias.proj_q, device_mesh=device_mesh)
        self.proj_k = LinearParamsReplicated(layer_local=attn_pair_bias.proj_k, device_mesh=device_mesh)
        self.proj_v = LinearParamsReplicated(layer_local=attn_pair_bias.proj_v, device_mesh=device_mesh)
        self.proj_g = LinearParamsReplicated(layer_local=attn_pair_bias.proj_g, device_mesh=device_mesh)
        self.proj_o = LinearParamsReplicated(layer_local=attn_pair_bias.proj_o, device_mesh=device_mesh)

        # (4) proj_z: Strip the Rearrange to avoid changing placements of z.
        #     When compute_pair_bias=True, serial proj_z is Sequential(LayerNorm, Linear, Rearrange)
        #     -> keep only LayerNorm and Linear. The permute is done manually in the forward pass.
        #     When compute_pair_bias=False, serial proj_z is just a Rearrange (no projection
        #     needed; z is already the pre-computed bias).
        if self.compute_pair_bias:
            self.proj_z = nn.Sequential(
                LayerNormParamsReplicated(attn_pair_bias.proj_z[0], device_mesh=device_mesh),
                LinearParamsReplicated(layer_local=attn_pair_bias.proj_z[1], device_mesh=device_mesh),
            )

    def forward(
        self,
        s: DTensor,
        z: DTensor,
        mask: DTensor,
        pair_mask: Union[DTensor, None] = None,
        multiplicity: int = 1,
        k_in: Union[DTensor, None] = None,
        model_cache: Union[OrderedDict, None] = None,
    ) -> DTensor:
        """Forward pass for ring attention with pair bias.

        Parameters
        ----------
        s : DTensor
            The input sequence tensor (queries), with shape (B, N, c_s) or (B*M, N, c_s)
            where M is multiplicity.
        z : DTensor
            The input pairwise tensor, with shape (B, N, N, c_z).
        mask : DTensor
            The token mask tensor with shape (B, N) or (B*M, N).
        pair_mask : DTensor or None, optional
            The pairwise mask tensor with shape (B, N, N). If None, only uses 1D mask.
        multiplicity : int, optional
            The diffusion batch size, by default 1.
        k_in : DTensor or None, optional
            Pre-computed key input tensor. If None, uses s as key input (k_in=s).
            For V2 API, caller should pass k_in explicitly.
        model_cache : OrderedDict or None, optional
            Cache for storing projected z during diffusion rollout. Only used if
            use_model_cache=True was set at init.

        Returns
        -------
        DTensor
            The output tensor, with shape (B*M, N, c_s).

        Raises
        ------
        ValueError
            If mask shape is incompatible with k_in.
        """
        # -------------------------------------------------
        # Begin DTensor ops
        #   DTensor metadata checks done for each operation
        # -------------------------------------------------
        if self.apply_initial_norm:
            s: DTensor = self.norm_s(s)  # Layer norm

        # V2 API: k_in is passed explicitly; V1 API: k_in defaults to s
        if k_in is None:
            k_in = s

        # Sanity check: mask should have same sequence length as k_in
        if mask.shape[-1] != k_in.shape[-2]:
            raise ValueError(
                f"mask sequence length ({mask.shape[-1]}) must match k_in sequence length ({k_in.shape[-2]}). "
                f"For V2 API with to_keys transformation, transform mask before passing."
            )

        # Compute projections
        q_proj_out: DTensor = self.proj_q(s)  # (B, N, c_s)
        k_proj_out: DTensor = self.proj_k(k_in)  # (B, N, c_s) or (B, H, c_s) if transformed
        v_proj_out: DTensor = self.proj_v(k_in)  # (B, N, c_s) or (B, H, c_s) if transformed
        g_proj_out: DTensor = self.proj_g(s)  # (B, N, c_s)

        # ------------------------------------------------------------
        # Project z to num_heads dimensions (V1: compute_pair_bias=True)
        # or use z as-is (V2: compute_pair_bias=False, z is pre-computed bias)
        # ------------------------------------------------------------
        if self.compute_pair_bias:
            #  input z: (B, N, N, c_z)
            #  output z: (B, N, N, num_heads) after proj_z (without Rearrange)
            if self.use_model_cache and model_cache is not None:
                if "z" not in model_cache:
                    z: DTensor = self.proj_z(z)  # (B, N, N, num_heads)
                    model_cache["z"] = z
                else:
                    z = model_cache["z"]
            else:
                z: DTensor = self.proj_z(z)  # (B, N, N, num_heads)
        # else: z is already the pre-computed bias with shape (B, N, N, num_heads)

        # ------------------------------------------------------------
        # Compute context vectors
        # ------------------------------------------------------------
        apb_context_vec_params = _AttentionPairBiasContextVecParams(
            ring_comm=self.ring_comm,
            multiplicity=multiplicity,
            num_heads=self.num_heads,
            head_dim=self.head_dim,
            inf=self.inf,
            use_window_batching=self.use_window_batching,
            sdpa_with_bias_backend=self.sdpa_with_bias_backend,
        )
        o_contex_vec = _AttentionPairBiasContexVecImpl.apply(
            q_proj_out,
            k_proj_out,
            v_proj_out,
            z,  # (B, N, N, H)
            mask,
            pair_mask,
            apb_context_vec_params,
        )
        # ------------------------------------------------------------
        # Gate and project context vectors
        # ------------------------------------------------------------
        gated_context_vec: DTensor = sigmoid_gate(x=o_contex_vec, g=g_proj_out)
        o: DTensor = self.proj_o(gated_context_vec)

        return o


class AttentionPairBiasShardwise(nn.Module):
    """Shardwise attention with pair bias for window-batched context parallelism.

    This module implements multi-head attention with pair bias specifically designed
    for window batching scenarios in context parallelism (CP). Unlike the standard
    `AttentionPairBias` which uses ring communication patterns, this implementation
    operates on sharded windows where each shard can be processed independently.

    The key difference from `AttentionPairBias`:
    - Used for window batching scenarios
    - Uses `to_keys` function OR pre-computed `k_in` to transform queries to key space
    - Operates on 4D single representations (B, K, W, D) and 5D pair representations
      (B, K, W, H, num_heads)
    - Does not apply multiplicity to z/mask, instead broadcasts them

    Configuration Flags
    -------------------
    The module supports both V1 (Boltz-1x) and V2 (Boltz-2) API styles:

    - apply_initial_norm: V1=True (has norm_s LayerNorm), V2=False (no initial norm)
    - compute_pair_bias: V1=True (always compute via LayerNorm+Linear),
      V2=False (z is pre-computed bias, no projection needed).
      Mutually exclusive with use_model_cache when False.
    - use_model_cache: V1=True (cache z projection), V2=False (no caching).
      Only valid when compute_pair_bias=True.

    Use Cases
    ---------
    V1 API (to_keys inside forward):
        AtomTransformer.forward(to_keys=to_keys)
            AttentionPairBiasShardwise.forward(to_keys=to_keys)

    V2 API (k_in pre-computed):
        DiffusionTransformerLayer.forward(to_keys=to_keys)
            k_in = to_keys(s)
            mask = to_keys(mask)
            AttentionPairBiasShardwise.forward(k_in=k_in, mask=mask)

    Attributes
    ----------
    c_s : int
        Hidden dimension of single representation (num_heads * head_dim).
    num_heads : int
        Number of attention heads.
    head_dim : int
        Dimension per attention head.
    inf : float
        Large value used for masking invalid positions.
    apply_initial_norm : bool
        Whether to apply layer normalization to input (V1=True, V2=False).
    compute_pair_bias : bool
        Whether to compute pair bias via LayerNorm+Linear (V1=True, V2=False).
    use_model_cache : bool
        Whether to cache z projection (V1=True, V2=False).
    device_mesh : DeviceMesh
        The device mesh for distributed computation.
    sdpa_with_bias_backend : SDPAWithBiasBackend
        Backend for scaled dot-product attention computation.
    """

    def __init__(
        self,
        attn_pair_bias: nn.Module,
        device_mesh: DeviceMesh,
        sdpa_with_bias_backend: SDPAWithBiasBackend = SDPAWithBiasBackend.REFERENCE,
        # Configuration flags for V1/V2 API compatibility
        apply_initial_norm: bool = False,  # V1=True, V2=False (default to V2)
        compute_pair_bias: bool = True,  # V1=True, V2=False
        use_model_cache: bool = False,  # V1=True, V2=False (default to V2)
    ) -> None:
        """Initialize the shardwise attention pair bias layer.

        Parameters
        ----------
        attn_pair_bias : nn.Module
            The serial attention pair bias layer to convert to DTensor.
        device_mesh : DeviceMesh
            The device mesh for distributed tensor operations.
        sdpa_with_bias_backend : SDPAWithBiasBackend, optional
            Backend for computing scaled dot-product attention with bias.
            Default is REFERENCE.
        apply_initial_norm : bool, optional
            Whether to apply LayerNorm to input s. V1=True, V2=False. Default False.
        compute_pair_bias : bool, optional
            Whether to compute pair bias (LayerNorm + Linear on z). V1=True, V2=False.
            Mutually exclusive with use_model_cache=True when compute_pair_bias=False.
            Default True.
        use_model_cache : bool, optional
            Whether to cache z projection for diffusion rollout. V1=True, V2=False.
            Only valid when compute_pair_bias=True. Default False.

        Raises
        ------
        TypeError
            If device_mesh is not a DeviceMesh instance, or if attn_pair_bias is not
            a recognized serial AttentionPairBias type.
        ValueError
            If use_model_cache=True with compute_pair_bias=False.
        """
        super().__init__()

        # (0) Type check on serial module
        if not isinstance(attn_pair_bias, (AttentionPairBiasSerialV1, AttentionPairBiasSerialV2)):
            raise TypeError(
                ", ".join(
                    [
                        f"Instance {attn_pair_bias} should have type "
                        f"{AttentionPairBiasSerialV1} or {AttentionPairBiasSerialV2}",
                        f"but instead has type {type(attn_pair_bias)}.",
                    ]
                )
            )

        # (1) Set non-module, non-parameter attributes
        self.c_s = attn_pair_bias.c_s
        self.num_heads = attn_pair_bias.num_heads
        self.head_dim = attn_pair_bias.head_dim
        self.inf = attn_pair_bias.inf

        # Configuration flags — compute_pair_bias and use_model_cache are mutually exclusive:
        # use_model_cache caches the z projection, which requires compute_pair_bias=True.
        if use_model_cache and not compute_pair_bias:
            raise ValueError(
                "use_model_cache=True requires compute_pair_bias=True because the cache "
                "stores the z projection output. Got compute_pair_bias=False."
            )
        self.apply_initial_norm = apply_initial_norm
        self.compute_pair_bias = compute_pair_bias
        self.use_model_cache = use_model_cache

        self.device_mesh = device_mesh
        # Mutable backend selection for scaled dot-product attention.
        # Default is REFERENCE. To switch backend for the entire model, use
        # ``model.apply(SetAttnPairBiasShardwiseBackend(backend))``
        # (see boltz.distributed.model.modules.utils.SetAttnPairBiasShardwiseBackend).
        self.sdpa_with_bias_backend = sdpa_with_bias_backend

        # (2) Sanity checks on non-module, non-parameter attributes
        if not isinstance(self.device_mesh, DeviceMesh):
            raise TypeError(f"Input '{device_mesh}' must be of type {DeviceMesh}. Got type {type(self.device_mesh)}.")

        # (3) Initialize child modules explicitly from serial module
        if self.apply_initial_norm:
            self.norm_s = LayerNormParamsReplicated(attn_pair_bias.norm_s, device_mesh=device_mesh)

        self.proj_q = LinearParamsReplicated(layer_local=attn_pair_bias.proj_q, device_mesh=device_mesh)
        self.proj_k = LinearParamsReplicated(layer_local=attn_pair_bias.proj_k, device_mesh=device_mesh)
        self.proj_v = LinearParamsReplicated(layer_local=attn_pair_bias.proj_v, device_mesh=device_mesh)
        self.proj_g = LinearParamsReplicated(layer_local=attn_pair_bias.proj_g, device_mesh=device_mesh)
        self.proj_o = LinearParamsReplicated(layer_local=attn_pair_bias.proj_o, device_mesh=device_mesh)

        # (4) proj_z: Strip the Rearrange to avoid changing placements of z.
        #     When compute_pair_bias=True, serial proj_z is Sequential(LayerNorm, Linear, Rearrange)
        #     -> keep only LayerNorm and Linear.
        #     When compute_pair_bias=False, serial proj_z is just a Rearrange (no projection
        #     needed; z is already the pre-computed bias).
        if self.compute_pair_bias:
            self.proj_z = nn.Sequential(
                LayerNormParamsReplicated(attn_pair_bias.proj_z[0], device_mesh=device_mesh),
                LinearParamsReplicated(layer_local=attn_pair_bias.proj_z[1], device_mesh=device_mesh),
            )

    def forward(
        self,
        s: DTensor,
        z: DTensor,
        mask: DTensor,
        to_keys: Union[Callable[[DTensor], DTensor], None] = None,
        k_in: Union[DTensor, None] = None,
        model_cache: Union[OrderedDict, None] = None,
    ) -> DTensor:
        """Forward pass for shardwise attention with pair bias.

        Computes multi-head attention with pair bias on window-batched inputs.
        The attention is computed within each window shard independently.

        Two API modes are supported:

        V1 API (to_keys provided):
            - to_keys transforms s to k_in internally
            - mask is transformed by to_keys internally
            - mask shape: (B, K, W) - query-aligned

        V2 API (k_in provided):
            - k_in is pre-computed by caller
            - mask is pre-transformed by caller to key-aligned shape
            - mask shape: (B, K, H) - key-aligned

        Parameters
        ----------
        s : DTensor
            Input single representation tensor with shape (B * M, K, W, c_s) where:
            - B is batch size
            - M is multiplicity (diffusion samples)
            - K is number of windows
            - W is window size (typically 32)
            - c_s is hidden dimension
        z : DTensor
            Input pair representation tensor with shape (B, K, W, H, c_z) where:
            - H is the attention key dimension (typically 128)
            - c_z is pair hidden dimension
            Note: z is NOT multiplied by M; it broadcasts along the multiplicity axis.
        mask : DTensor
            Mask tensor indicating valid positions.
            - V1 API (to_keys): shape (B, K, W) - will be transformed to (B, K, H)
            - V2 API (k_in): shape (B, K, H) - already key-aligned
        to_keys : Callable or None, optional
            Function to transform tensors from query space (B, K, W, ...) to
            key space (B, K, H, ...). Mutually exclusive with k_in.
        k_in : DTensor or None, optional
            Pre-computed key input tensor with shape (B*M, K, H, c_s).
            Mutually exclusive with to_keys.
        model_cache : OrderedDict or None, optional
            Cache for storing projected z during diffusion rollout. Only used if
            use_model_cache=True was set at init.

        Returns
        -------
        DTensor
            Output tensor with shape (B * M, K, W, c_s).

        Raises
        ------
        ValueError
            If both to_keys and k_in are provided (mutually exclusive).
            If neither to_keys nor k_in is provided.
            If s does not have 4 dimensions.
            If z does not have 5 dimensions.
            If mask dimensions don't match expected shapes.
            If s.shape[0] is not divisible by z.shape[0].

        Notes
        -----
        This module avoids the multiplicity memory overhead by broadcasting z and mask
        along the multiplicity dimension rather than replicating them.
        """
        # Check mutual exclusivity of to_keys and k_in
        if to_keys is not None and k_in is not None:
            raise ValueError("to_keys and k_in are mutually exclusive. Provide only one.")
        if to_keys is None and k_in is None:
            raise ValueError("Either to_keys or k_in must be provided.")

        # Shape validations
        if s.ndim != 4:
            raise ValueError(f"s must have 4 dimensions (B*M, K, W, D), but got s.ndim={s.ndim}")
        if z.ndim != 5:
            raise ValueError(f"z must have 5 dimensions (B, K, W, H, c_z), but got z.ndim={z.ndim}")

        if s.shape[1:3] != z.shape[1:3]:
            raise ValueError(
                f"s.shape[1:3] must be equal to z.shape[1:3], but got s.shape[1:3]={s.shape[1:3]} "
                f"and z.shape[1:3]={z.shape[1:3]}"
            )

        if s.shape[0] % z.shape[0] != 0:
            # NOTE: this module doesn't apply multiplicity to z because it broadcasts z (and mask)
            # to the attention score by design. This avoids multiplying the memory storage of
            # the pair representation throughout the entire AtomTransformer and its submodules.
            raise ValueError(
                f"s.shape[0] must be divisible by z.shape[0], but got s.shape[0]={s.shape[0]} "
                f"and z.shape[0]={z.shape[0]}"
            )

        # Validate mask shape based on API mode
        if mask is not None:
            if to_keys is not None:
                # V1 API: mask should be query-aligned (B, K, W)
                if mask.ndim != 3:
                    raise ValueError(f"V1 API: mask must have 3 dimensions (B, K, W), but got mask.ndim={mask.ndim}")
                if mask.shape != z.shape[:3]:
                    raise ValueError(
                        f"V1 API: mask.shape must equal z.shape[:3], but got mask.shape={mask.shape} "
                        f"and z.shape[:3]={z.shape[:3]}"
                    )
            else:
                # V2 API: mask should be key-aligned (B, K, H)
                if mask.ndim != 3:
                    raise ValueError(f"V2 API: mask must have 3 dimensions (B, K, H), but got mask.ndim={mask.ndim}")
                # For V2 API, mask.shape[2] should be H (key dimension), not W (query dimension)
                if mask.shape[:2] != z.shape[:2]:
                    raise ValueError(
                        f"V2 API: mask.shape[:2] must equal z.shape[:2], but got mask.shape[:2]={mask.shape[:2]} "
                        f"and z.shape[:2]={z.shape[:2]}"
                    )

        # -------------------------------------------------
        # Begin DTensor ops
        # -------------------------------------------------
        if self.apply_initial_norm:
            s: DTensor = self.norm_s(s)  # Layer norm

        # Compute k_in and mask_key based on API mode
        if to_keys is not None:
            # V1 API: transform s and mask using to_keys
            k_in_computed = to_keys(s)
            mask_key = to_keys(mask)
        else:
            # V2 API: use provided k_in and mask (already key-aligned)
            k_in_computed = k_in
            mask_key = mask

        # Project z to num_heads dimensions (V1: compute_pair_bias=True)
        # or use z as-is (V2: compute_pair_bias=False, z is pre-computed bias)
        if self.compute_pair_bias:
            # V1: project z through LayerNorm+Linear, optionally cache the result
            if self.use_model_cache and model_cache is not None:
                if "z" not in model_cache:
                    z: DTensor = self.proj_z(z)  # (B, K, W, H, num_heads)
                    model_cache["z"] = z
                else:
                    z = model_cache["z"]
            else:
                z: DTensor = self.proj_z(z)  # (B, K, W, H, num_heads)
        # else: V2 — z is already the pre-computed bias with shape
        # (B, K, W, H, num_heads), no projection needed.

        # Compute projections
        q_proj_out: DTensor = self.proj_q(s)  # (B, K, W, c_s)
        k_proj_out: DTensor = self.proj_k(k_in_computed)  # (B, K, H, c_s)
        v_proj_out: DTensor = self.proj_v(k_in_computed)  # (B, K, H, c_s)
        g_proj_out: DTensor = self.proj_g(s)  # (B, K, W, c_s)

        o = _AttentionPairBiasShardwiseImpl.apply(
            q_proj_out,
            k_proj_out,
            v_proj_out,
            z,
            mask_key,
            self.sdpa_with_bias_backend,
            self.num_heads,
            self.head_dim,
            self.inf,
        )

        # ------------------------------------------------------------
        # Gate and project context vectors
        # ------------------------------------------------------------
        o_gated: DTensor = sigmoid_gate(x=o, g=g_proj_out)
        # (B, K, W, c_s)
        o_final: DTensor = self.proj_o(o_gated)

        return o_final
