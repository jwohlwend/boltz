from typing import Optional, Union

import torch
from torch import Tensor, nn

import boltz.model.layers.initialize as init


class Transition(nn.Module):
    """Perform a two-layer MLP."""

    def __init__(
        self,
        dim: int = 128,
        hidden: int = 512,
        out_dim: Optional[int] = None,
    ) -> None:
        """Initialize the TransitionUpdate module.

        Parameters
        ----------
        dim: int
            The dimension of the input, default 128
        hidden: int
            The dimension of the hidden, default 512
        out_dim: Optional[int]
            The dimension of the output, default None

        """
        super().__init__()
        if out_dim is None:
            out_dim = dim

        self.norm = nn.LayerNorm(dim, eps=1e-5)
        self.fc1 = nn.Linear(dim, hidden, bias=False)
        self.fc2 = nn.Linear(dim, hidden, bias=False)
        self.fc3 = nn.Linear(hidden, out_dim, bias=False)
        self.silu = nn.SiLU()
        self.hidden = hidden

        init.bias_init_one_(self.norm.weight)
        init.bias_init_zero_(self.norm.bias)

        init.lecun_normal_init_(self.fc1.weight)
        init.lecun_normal_init_(self.fc2.weight)
        init.final_init_(self.fc3.weight)

    def forward(
        self,
        x: Tensor,
        chunk_size: Optional[int] = None,
        return_intermediates: bool = False,
    ) -> Union[Tensor, tuple[Tensor, dict[str, Tensor]]]:
        """Perform a forward pass.

        Parameters
        ----------
        x: torch.Tensor
            The input data of shape (..., D)

        Returns
        -------
        x: torch.Tensor
            The output data of shape (..., D)

        """
        x = self.norm(x)
        x_norm = x

        if chunk_size is None or self.training:
            fc1 = self.fc1(x)
            fc2 = self.fc2(x)
            hidden = self.silu(fc1) * fc2
            out = self.fc3(hidden)
            if return_intermediates:
                return out, {
                    "x_norm": x_norm,
                    "fc1": fc1,
                    "fc2": fc2,
                    "hidden": hidden,
                }
            return out

        # Compute in chunks
        fc1_chunks = []
        fc2_chunks = []
        hidden_chunks = []
        for i in range(0, self.hidden, chunk_size):
            fc1_slice = self.fc1.weight[i : i + chunk_size, :]
            fc2_slice = self.fc2.weight[i : i + chunk_size, :]
            fc3_slice = self.fc3.weight[:, i : i + chunk_size]
            fc1_chunk = x @ fc1_slice.T
            fc2_chunk = x @ fc2_slice.T
            x_chunk = self.silu(fc1_chunk) * fc2_chunk
            out = x_chunk @ fc3_slice.T if i == 0 else out + x_chunk @ fc3_slice.T
            if return_intermediates:
                fc1_chunks.append(fc1_chunk)
                fc2_chunks.append(fc2_chunk)
                hidden_chunks.append(x_chunk)
        if return_intermediates:
            fc1 = torch.cat(fc1_chunks, dim=-1)
            fc2 = torch.cat(fc2_chunks, dim=-1)
            hidden = torch.cat(hidden_chunks, dim=-1)
            return out, {
                "x_norm": x_norm,
                "fc1": fc1,
                "fc2": fc2,
                "hidden": hidden,
            }
        return out
