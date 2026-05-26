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

from boltz.model.layers import initialize as init


@torch.compiler.disable
def kernel_triangular_mult(
    x,
    direction,
    mask,
    norm_in_weight,
    norm_in_bias,
    p_in_weight,
    g_in_weight,
    norm_out_weight,
    norm_out_bias,
    p_out_weight,
    g_out_weight,
    eps,
):
    from cuequivariance_torch.primitives.triangle import triangle_multiplicative_update

    return triangle_multiplicative_update(
        x,
        direction=direction,
        mask=mask,
        norm_in_weight=norm_in_weight,
        norm_in_bias=norm_in_bias,
        p_in_weight=p_in_weight,
        g_in_weight=g_in_weight,
        norm_out_weight=norm_out_weight,
        norm_out_bias=norm_out_bias,
        p_out_weight=p_out_weight,
        g_out_weight=g_out_weight,
        eps=eps,
    )


class TriangleMultiplicationOutgoing(nn.Module):
    """TriangleMultiplicationOutgoing."""

    def __init__(self, dim: int = 128) -> None:
        """Initialize the TriangularUpdate module.

        Parameters
        ----------
        dim: int
            The dimension of the input, default 128

        """
        super().__init__()

        self.norm_in = nn.LayerNorm(dim, eps=1e-5)
        self.p_in = nn.Linear(dim, 2 * dim, bias=False)
        self.g_in = nn.Linear(dim, 2 * dim, bias=False)

        self.norm_out = nn.LayerNorm(dim)
        self.p_out = nn.Linear(dim, dim, bias=False)
        self.g_out = nn.Linear(dim, dim, bias=False)

        init.bias_init_one_(self.norm_in.weight)
        init.bias_init_zero_(self.norm_in.bias)

        init.lecun_normal_init_(self.p_in.weight)
        init.gating_init_(self.g_in.weight)

        init.bias_init_one_(self.norm_out.weight)
        init.bias_init_zero_(self.norm_out.bias)

        init.final_init_(self.p_out.weight)
        init.gating_init_(self.g_out.weight)

    def forward(self, x: Tensor, mask: Tensor, use_kernels: bool = False) -> Tensor:
        """Perform a forward pass.

        Parameters
        ----------
        x: torch.Tensor
            The input data of shape (B, N, N, D)
        mask: torch.Tensor
            The input mask of shape (B, N, N)
        use_kernels: bool
            Whether to use the kernel

        Returns
        -------
        x: torch.Tensor
            The output data of shape (B, N, N, D)

        """
        if use_kernels:
            return kernel_triangular_mult(
                x,
                direction="outgoing",
                mask=mask,
                norm_in_weight=self.norm_in.weight,
                norm_in_bias=self.norm_in.bias,
                p_in_weight=self.p_in.weight,
                g_in_weight=self.g_in.weight,
                norm_out_weight=self.norm_out.weight,
                norm_out_bias=self.norm_out.bias,
                p_out_weight=self.p_out.weight,
                g_out_weight=self.g_out.weight,
                eps=1e-5,
            )

        # Input gating: D -> D
        x = self.norm_in(x)
        x_in = x
        x = self.p_in(x) * self.g_in(x).sigmoid()

        # Apply mask
        x = x * mask.unsqueeze(-1)

        # Split input and cast to float
        safe_dtype = torch.promote_types(x.dtype, torch.float32)
        a, b = torch.chunk(x.to(dtype=safe_dtype), 2, dim=-1)

        # Triangular projection
        x = torch.einsum("bikd,bjkd->bijd", a, b)

        # Output gating
        x = self.p_out(self.norm_out(x)) * self.g_out(x_in).sigmoid()

        return x


class TriangleMultiplicationIncoming(nn.Module):
    """TriangleMultiplicationIncoming."""

    def __init__(self, dim: int = 128) -> None:
        """Initialize the TriangularUpdate module.

        Parameters
        ----------
        dim: int
            The dimension of the input, default 128

        """
        super().__init__()

        self.norm_in = nn.LayerNorm(dim, eps=1e-5)
        self.p_in = nn.Linear(dim, 2 * dim, bias=False)
        self.g_in = nn.Linear(dim, 2 * dim, bias=False)

        self.norm_out = nn.LayerNorm(dim)
        self.p_out = nn.Linear(dim, dim, bias=False)
        self.g_out = nn.Linear(dim, dim, bias=False)

        init.bias_init_one_(self.norm_in.weight)
        init.bias_init_zero_(self.norm_in.bias)

        init.lecun_normal_init_(self.p_in.weight)
        init.gating_init_(self.g_in.weight)

        init.bias_init_one_(self.norm_out.weight)
        init.bias_init_zero_(self.norm_out.bias)

        init.final_init_(self.p_out.weight)
        init.gating_init_(self.g_out.weight)

    def forward(self, x: Tensor, mask: Tensor, use_kernels: bool = False) -> Tensor:
        """Perform a forward pass.

        Parameters
        ----------
        x: torch.Tensor
            The input data of shape (B, N, N, D)
        mask: torch.Tensor
            The input mask of shape (B, N, N)
        use_kernels: bool
            Whether to use the kernel

        Returns
        -------
        x: torch.Tensor
            The output data of shape (B, N, N, D)

        """
        if use_kernels:
            return kernel_triangular_mult(
                x,
                direction="incoming",
                mask=mask,
                norm_in_weight=self.norm_in.weight,
                norm_in_bias=self.norm_in.bias,
                p_in_weight=self.p_in.weight,
                g_in_weight=self.g_in.weight,
                norm_out_weight=self.norm_out.weight,
                norm_out_bias=self.norm_out.bias,
                p_out_weight=self.p_out.weight,
                g_out_weight=self.g_out.weight,
                eps=1e-5,
            )

        # Input gating: D -> D
        x = self.norm_in(x)
        x_in = x
        x = self.p_in(x) * self.g_in(x).sigmoid()

        # Apply mask
        x = x * mask.unsqueeze(-1)

        # Split input and cast to float
        safe_dtype = torch.promote_types(x.dtype, torch.float32)
        a, b = torch.chunk(x.to(dtype=safe_dtype), 2, dim=-1)

        # Triangular projection
        x = torch.einsum("bkid,bkjd->bijd", a, b)

        # Output gating
        x = self.p_out(self.norm_out(x)) * self.g_out(x_in).sigmoid()

        return x
