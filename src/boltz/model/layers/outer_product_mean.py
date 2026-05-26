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


import torch
from torch import Tensor, nn

import boltz.model.layers.initialize as init


class OuterProductMean(nn.Module):
    """Outer product mean layer."""

    def __init__(self, c_in: int, c_hidden: int, c_out: int) -> None:
        """Initialize the outer product mean layer.

        Parameters
        ----------
        c_in : int
            The input dimension.
        c_hidden : int
            The hidden dimension.
        c_out : int
            The output dimension.

        """
        super().__init__()
        self.c_hidden = c_hidden
        self.norm = nn.LayerNorm(c_in)
        self.proj_a = nn.Linear(c_in, c_hidden, bias=False)
        self.proj_b = nn.Linear(c_in, c_hidden, bias=False)
        self.proj_o = nn.Linear(c_hidden * c_hidden, c_out)
        init.final_init_(self.proj_o.weight)
        init.final_init_(self.proj_o.bias)

    def forward(self, m: Tensor, mask: Tensor, chunk_size: int = None) -> Tensor:
        """Forward pass.

        Parameters
        ----------
        m : torch.Tensor
            The sequence tensor (B, S, N, c_in).
        mask : torch.Tensor
            The mask tensor (B, S, N).

        Returns
        -------
        torch.Tensor
            The output tensor (B, N, N, c_out).

        """
        # Expand mask
        mask = mask.unsqueeze(-1).to(m)

        # Compute projections
        m = self.norm(m)
        a = self.proj_a(m) * mask
        b = self.proj_b(m) * mask

        # Compute outer product mean
        if chunk_size is not None and not self.training:
            # Compute pairwise mask
            for i in range(0, mask.shape[1], 64):
                if i == 0:
                    num_mask = (mask[:, i : i + 64, None, :] * mask[:, i : i + 64, :, None]).sum(1)
                else:
                    num_mask += (mask[:, i : i + 64, None, :] * mask[:, i : i + 64, :, None]).sum(1)
            num_mask = num_mask.clamp(min=1)

            # Compute squentially in chunks
            for i in range(0, self.c_hidden, chunk_size):
                a_chunk = a[:, :, :, i : i + chunk_size]
                sliced_weight_proj_o = self.proj_o.weight[:, i * self.c_hidden : (i + chunk_size) * self.c_hidden]

                z = torch.einsum("bsic,bsjd->bijcd", a_chunk, b)
                z = z.reshape(*z.shape[:3], -1)
                z = z / num_mask

                # Project to output
                if i == 0:
                    z_out = z.to(m) @ sliced_weight_proj_o.T
                else:
                    z_out = z_out + z.to(m) @ sliced_weight_proj_o.T

            z_out = z_out + self.proj_o.bias  # add bias
            return z_out
        else:
            mask = mask[:, :, None, :] * mask[:, :, :, None]
            num_mask = mask.sum(1).clamp(min=1)
            # Cast to at least float32 for numerical stability, using
            # promote_types to preserve higher-precision dtypes (e.g. float64).
            compute_dtype = torch.promote_types(a.dtype, torch.float32)
            z = torch.einsum("bsic,bsjd->bijcd", a.to(compute_dtype), b.to(compute_dtype))
            z = z.reshape(*z.shape[:3], -1)
            z = z / num_mask

            # Project to output
            z = self.proj_o(z.to(m))
            return z
