import torch
from torch import Tensor
from torch.nn.functional import pad


def pad_dim(data: Tensor, dim: int, pad_len: float, value: float = 0) -> Tensor:
    """Pad a tensor along a given dimension.

    Parameters
    ----------
    data : Tensor
        The input tensor.
    dim : int
        The dimension to pad.
    pad_len : float
        The padding length.
    value : int, optional
        The value to pad with.

    Returns
    -------
    Tensor
        The padded tensor.

    """
    if pad_len == 0:
        return data

    total_dims = len(data.shape)
    padding = [0] * (2 * (total_dims - dim))
    padding[2 * (total_dims - 1 - dim) + 1] = pad_len
    return pad(data, tuple(padding), value=value)


def pad_to_max(data: list[Tensor], value: float = 0) -> tuple[Tensor, Tensor]:
    """Pad the data in all dimensions to the maximum found.

    Parameters
    ----------
    data : list[Tensor]
        list of tensors to pad.
    value : float
        The value to use for padding.

    Returns
    -------
    Tensor
        The padded tensor.
    Tensor
        The padding mask.

    """
    if isinstance(data[0], str):
        return data, 0

    # Check if all have the same shape
    if all(d.shape == data[0].shape for d in data):
        return torch.stack(data, dim=0), 0

    # Get the maximum in each dimension
    num_dims = len(data[0].shape)
    max_dims = [max(d.shape[i] for d in data) for i in range(num_dims)]

    # Get the padding lengths
    pad_lengths = []
    for d in data:
        dims = []
        for i in range(num_dims):
            dims.append(0)
            dims.append(max_dims[num_dims - i - 1] - d.shape[num_dims - i - 1])
        pad_lengths.append(dims)

    # Pad the data
    padding = [
        pad(torch.ones_like(d), pad_len, value=0)
        for d, pad_len in zip(data, pad_lengths)
    ]
    data = [pad(d, pad_len, value=value) for d, pad_len in zip(data, pad_lengths)]

    # Stack the data
    padding = torch.stack(padding, dim=0)
    data = torch.stack(data, dim=0)

    return data, padding
