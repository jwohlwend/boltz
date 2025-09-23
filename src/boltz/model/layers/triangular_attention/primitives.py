# Copyright 2021 AlQuraishi Laboratory
# Copyright 2021 DeepMind Technologies Limited
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import math
from typing import Callable, List, Optional, Tuple

import torch
from einops import rearrange
from torch import nn

from boltz.model.layers import initialize
from boltz.model.layers.triangular_attention.utils import (
    flatten_final_dims,
    permute_final_dims,
)


class Linear(nn.Linear):
    """
    A Linear layer with built-in nonstandard initializations. Called just
    like torch.nn.Linear.

    Implements the initializers in 1.11.4, plus some additional ones found
    in the code.
    """

    def __init__(
        self,
        in_dim: int,
        out_dim: int,
        bias: bool = True,
        init: str = "default",
        init_fn: Optional[Callable[[torch.Tensor, torch.Tensor], None]] = None,
        precision=None,
    ):
        """Initialize the linear layer.

        Parameters
        ----------
        in_dim : int
            The final dimension of inputs to the layer
        out_dim : int
            The final dimension of layer outputs
        bias : bool, default=True
            Whether to learn an additive bias
        init : str, default='default'
            The initializer to use. Choose from:

            - "default": LeCun fan-in truncated normal initialization
            - "relu": He initialization w/ truncated normal distribution
            - "glorot": Fan-average Glorot uniform initialization
            - "gating": Weights=0, Bias=1
            - "normal": Normal initialization with std=1/sqrt(fan_in)
            - "final": Weights=0, Bias=0

            Overridden by init_fn if the latter is not None.
        init_fn : callable, optional
            A custom initializer taking weight and bias as inputs.
            Overrides init if not None.

        """
        super().__init__(in_dim, out_dim, bias=bias)

        if bias:
            with torch.no_grad():
                self.bias.fill_(0)

        with torch.no_grad():
            if init_fn is not None:
                init_fn(self.weight, self.bias)
            else:
                if init == "default":
                    initialize.lecun_normal_init_(self.weight)
                elif init == "relu":
                    initialize.he_normal_init_(self.weight)
                elif init == "glorot":
                    initialize.glorot_uniform_init_(self.weight)
                elif init == "gating":
                    initialize.gating_init_(self.weight)
                    if bias:
                        self.bias.fill_(1.0)
                elif init == "normal":
                    initialize.normal_init_(self.weight)
                elif init == "final":
                    initialize.final_init_(self.weight)
                else:
                    raise ValueError("Invalid init string.")

        self.precision = precision

    def forward(self, input: torch.Tensor) -> torch.Tensor:
        d = input.dtype
        if self.precision is not None:
            with torch.autocast("cuda", enabled=False):
                bias = (
                    self.bias.to(dtype=self.precision)
                    if self.bias is not None
                    else None
                )
                return nn.functional.linear(
                    input.to(dtype=self.precision),
                    self.weight.to(dtype=self.precision),
                    bias,
                ).to(dtype=d)

        if d is torch.bfloat16:
            with torch.autocast("cuda", enabled=False):
                bias = self.bias.to(dtype=d) if self.bias is not None else None
                return nn.functional.linear(input, self.weight.to(dtype=d), bias)

        return nn.functional.linear(input, self.weight, self.bias)


class LayerNorm(nn.Module):
    def __init__(self, c_in, eps=1e-5):
        super(LayerNorm, self).__init__()

        self.c_in = (c_in,)
        self.eps = eps

        self.weight = nn.Parameter(torch.ones(c_in))
        self.bias = nn.Parameter(torch.zeros(c_in))

    def forward(self, x):
        d = x.dtype
        if d is torch.bfloat16:
            with torch.autocast("cuda", enabled=False):
                out = nn.functional.layer_norm(
                    x,
                    self.c_in,
                    self.weight.to(dtype=d),
                    self.bias.to(dtype=d),
                    self.eps,
                )
        else:
            out = nn.functional.layer_norm(
                x,
                self.c_in,
                self.weight,
                self.bias,
                self.eps,
            )

        return out


@torch.jit.ignore
def softmax_no_cast(t: torch.Tensor, dim: int = -1) -> torch.Tensor:
    """
    Softmax, but without automatic casting to fp32 when the input is of
    type bfloat16
    """
    d = t.dtype
    if d is torch.bfloat16:
        with torch.autocast("cuda", enabled=False):
            s = torch.nn.functional.softmax(t, dim=dim)
    else:
        s = torch.nn.functional.softmax(t, dim=dim)

    return s


# @torch.jit.script
def _attention(
    query: torch.Tensor,
    key: torch.Tensor,
    value: torch.Tensor,
    biases: List[torch.Tensor],
) -> torch.Tensor:
    # [*, H, C_hidden, K]
    key = permute_final_dims(key, (1, 0))

    # [*, H, Q, K]
    a = torch.matmul(query, key)

    for b in biases:
        a += b

    a = softmax_no_cast(a, -1)

    # [*, H, Q, C_hidden]
    a = torch.matmul(a, value)

    return a


@torch.compiler.disable
def kernel_triangular_attn(q, k, v, tri_bias, mask, scale):
    from cuequivariance_torch.primitives.triangle import triangle_attention
    return triangle_attention(q, k, v, tri_bias, mask=mask, scale=scale)


class Attention(nn.Module):
    """
    Standard multi-head attention using AlphaFold's default layer
    initialization. Allows multiple bias vectors.
    """

    def __init__(
        self,
        c_q: int,
        c_k: int,
        c_v: int,
        c_hidden: int,
        no_heads: int,
        gating: bool = True,
    ):
        """Initialize the attention layer.

        Parameters
        ----------
        c_q : int
            Input dimension of query data
        c_k : int
            Input dimension of key data
        c_v : int
            Input dimension of value data
        c_hidden : int
            Per-head hidden dimension
        no_heads : int
            Number of attention heads
        gating : bool, default=True
            Whether the output should be gated using query data

        """
        super().__init__()

        self.c_q = c_q
        self.c_k = c_k
        self.c_v = c_v
        self.c_hidden = c_hidden
        self.no_heads = no_heads
        self.gating = gating

        # DISCREPANCY: c_hidden is not the per-head channel dimension, as
        # stated in the supplement, but the overall channel dimension.

        self.linear_q = Linear(
            self.c_q, self.c_hidden * self.no_heads, bias=False, init="glorot"
        )
        self.linear_k = Linear(
            self.c_k, self.c_hidden * self.no_heads, bias=False, init="glorot"
        )
        self.linear_v = Linear(
            self.c_v, self.c_hidden * self.no_heads, bias=False, init="glorot"
        )
        self.linear_o = Linear(
            self.c_hidden * self.no_heads, self.c_q, bias=False, init="final"
        )

        self.linear_g = None
        if self.gating:
            self.linear_g = Linear(
                self.c_q, self.c_hidden * self.no_heads, bias=False, init="gating"
            )

        self.sigmoid = nn.Sigmoid()

    def _prep_qkv(
        self, q_x: torch.Tensor, kv_x: torch.Tensor, apply_scale: bool = True
    ) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        # [*, Q/K/V, H * C_hidden]
        q = self.linear_q(q_x)
        k = self.linear_k(kv_x)
        v = self.linear_v(kv_x)

        # [*, Q/K, H, C_hidden]
        q = q.view(q.shape[:-1] + (self.no_heads, -1))
        k = k.view(k.shape[:-1] + (self.no_heads, -1))
        v = v.view(v.shape[:-1] + (self.no_heads, -1))

        # [*, H, Q/K, C_hidden]
        q = q.transpose(-2, -3)
        k = k.transpose(-2, -3)
        v = v.transpose(-2, -3)

        if apply_scale:
            q /= math.sqrt(self.c_hidden)

        return q, k, v

    def _wrap_up(self, o: torch.Tensor, q_x: torch.Tensor) -> torch.Tensor:
        if self.linear_g is not None:
            g = self.sigmoid(self.linear_g(q_x))

            # [*, Q, H, C_hidden]
            g = g.view(g.shape[:-1] + (self.no_heads, -1))
            o = o * g

        # [*, Q, H * C_hidden]
        o = flatten_final_dims(o, 2)

        # [*, Q, C_q]
        o = self.linear_o(o)

        return o

    def forward(
        self,
        q_x: torch.Tensor,
        kv_x: torch.Tensor,
        tri_bias: torch.Tensor,
        mask_bias: torch.Tensor,
        mask: torch.Tensor,
        use_kernels: bool = False,
    ) -> torch.Tensor:
        """Compute attention.

        Parameters
        ----------
        q_x : torch.Tensor
            [*, Q, C_q] query data
        kv_x : torch.Tensor
            [*, K, C_k] key data
        tri_bias : torch.Tensor
            [*, H, Q, K] triangular bias
        mask_bias : torch.Tensor
            [*, H, Q, K] mask bias
        mask : torch.Tensor
            [*, Q, K] mask
        use_kernels : bool, default=False
            Whether to use optimized CUDA kernels

        Returns
        -------
            [*, Q, C_q] attention update

        """
        # Attention kernel applies scaling internally
        q, k, v = self._prep_qkv(
            q_x,
            kv_x,
            apply_scale=not use_kernels,
        )

        if use_kernels:
            scale = 1.0 / math.sqrt(self.c_hidden)
            o = kernel_triangular_attn(
                q,
                k,
                v,
                tri_bias=tri_bias,
                mask=mask.bool(),
                scale=scale,
            )
            o = o.transpose(-2, -3)
        else:
            biases = [mask_bias, tri_bias]
            o = _attention(q, k, v, biases)
            o = o.transpose(-2, -3)

        o = self._wrap_up(o, q_x)

        return o


def _trifast_attn(q, k, v, biases):
    orig_n_dims = len(q.shape)

    if len(biases) != 2:
        raise ValueError(f"Trifast expects two bias terms, found {len(biases)}")

    mask, b = biases

    if len(b.shape) == 5:
        # Sometimes there is an extra batch dim -- why?
        b = b.squeeze(1)

    if orig_n_dims == 4:
        # add fake batch dim
        q = q.unsqueeze(0)
        k = k.unsqueeze(0)
        v = v.unsqueeze(0)
        # b = b.unsqueeze(0) not sure why this and only this has a batch dim?
        mask = mask.unsqueeze(0)

    if len(q.shape) != 5:
        raise ValueError(f"Trifast expects q/k/v to be 5D, found {len(q.shape)}")

    # Reorder q/k/v
    q = rearrange(q, "b i h j d -> b h i j d")
    k = rearrange(k, "b i h j d -> b h i j d")
    v = rearrange(v, "b i h j d -> b h i j d")

    # Make mask the right shape.
    mask = rearrange(mask, "b i () () j -> b i j").bool()

    # Delay import to here to avoid initializing cuda too early
    from trifast import triangle_attention

    o = triangle_attention(q, k, v, b, mask)
    o = rearrange(o, "b h i j d -> b i j h d")

    # Remove the batch dim if we added it.
    if orig_n_dims == 4:
        o = o.squeeze(0)
    return o
