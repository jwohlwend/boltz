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

from functools import partial, partialmethod
from typing import Optional

import torch
import torch.nn as nn

from boltz.model.layers.triangular_attention.primitives import (
    Attention,
    LayerNorm,
    Linear,
)
from boltz.model.layers.triangular_attention.utils import (
    chunk_layer,
    permute_final_dims,
)


class TriangleAttention(nn.Module):
    """Implement Algorithm 12."""

    def __init__(
        self,
        c_in: int,
        c_hidden: int,
        no_heads: int,
        starting: bool = True,
        inf: float = 1e9,
    ) -> None:
        super().__init__()

        self.c_in = c_in
        self.c_hidden = c_hidden
        self.no_heads = no_heads
        self.starting = starting
        self.inf = inf

        self.layer_norm = LayerNorm(self.c_in)

        self.linear = Linear(c_in, self.no_heads, bias=False, init="normal")

        self.mha = Attention(
            self.c_in, self.c_in, self.c_in, self.c_hidden, self.no_heads
        )

    @torch.jit.ignore
    def _chunk(
        self,
        x: torch.Tensor,
        tri_bias: torch.Tensor,
        mask_bias: torch.Tensor,
        mask: torch.Tensor,
        chunk_size: int,
        use_kernels: bool = False,
    ) -> torch.Tensor:
        """Compute triangle attention.

        Parameters
        ----------
        x : torch.Tensor
            Input tensor of shape [*, I, J, C_in]
        biases : list[torch.Tensor]
            List of bias tensors of shape [*, H, I, J]
        chunk_size : int
            Size of chunks for memory efficient computation
        use_kernels : bool, default=False
            Whether to use optimized CUDA kernels

        Returns
        -------
        torch.Tensor
            Output tensor of shape [*, I, J, C_in]

        """
        mha_inputs = {
            "q_x": x,
            "kv_x": x,
            "tri_bias": tri_bias,
            "mask_bias": mask_bias,
            "mask": mask,
        }

        return chunk_layer(
            partial(
                self.mha,
                use_kernels=use_kernels,
            ),
            mha_inputs,
            chunk_size=chunk_size,
            no_batch_dims=len(x.shape[:-2]),
            _out=None,
        )

    def forward(
        self,
        x: torch.Tensor,
        mask: Optional[torch.Tensor] = None,
        chunk_size: Optional[int] = None,
        use_kernels: bool = False,
    ) -> torch.Tensor:
        """Compute triangle attention.

        Parameters
        ----------
        x : torch.Tensor
            Input tensor of shape [*, I, J, C_in]
        mask : torch.Tensor, optional
            Attention mask of shape [*, I, J]
        chunk_size : int, optional
            Size of chunks for memory efficient computation
        use_kernels : bool, default=False
            Whether to use optimized CUDA kernels

        Returns
        -------
        torch.Tensor
            Output tensor of shape [*, I, J, C_in]

        """
        if mask is None:
            # [*, I, J]
            mask = x.new_ones(
                x.shape[:-1],
            )

        if not self.starting:
            x = x.transpose(-2, -3)
            mask = mask.transpose(-1, -2)

        # [*, I, J, C_in]
        x = self.layer_norm(x)

        # [*, I, 1, 1, J]
        mask = mask[..., :, None, None, :]
        mask_bias = self.inf * (mask - 1)

        # [*, H, I, J]
        triangle_bias = permute_final_dims(self.linear(x), (2, 0, 1))

        # [*, 1, H, I, J]
        triangle_bias = triangle_bias.unsqueeze(-4)

        if chunk_size is not None and not use_kernels:
            x = self._chunk(
                x,
                triangle_bias,
                mask_bias,
                mask,
                chunk_size,
                use_kernels=use_kernels,
            )
        else:
            x = self.mha(
                x,
                x,
                triangle_bias,
                mask_bias,
                mask,
                use_kernels=use_kernels,
            )

        if not self.starting:
            x = x.transpose(-2, -3)

        return x


# Implements Algorithm 13
TriangleAttentionStartingNode = TriangleAttention


class TriangleAttentionEndingNode(TriangleAttention):
    """Implement Algorithm 14."""

    __init__ = partialmethod(TriangleAttention.__init__, starting=False)
