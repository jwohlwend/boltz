import torch


def compute_relative_distribution_perfect_correlation(
    binned_distribution_1, binned_distribution_2
):
    """
    Compute the relative distribution between two binned distributions with perfect correlation.

    Parameters
    ----------
    binned_distribution_1 : torch.Tensor
        The first binned distribution, shape (..., K).
    binned_distribution_2 : torch.Tensor
        The second binned distribution, shape (..., K).

    Returns
    -------
    torch.Tensor
        The relative distribution, shape (..., 2K - 1).

    """
    K = binned_distribution_1.shape[-1]
    relative_distribution = torch.zeros(
        binned_distribution_1.shape[:-1] + (2 * K - 1,),
        device=binned_distribution_1.device,
    )
    zero = torch.zeros(
        binned_distribution_1.shape[:-1] + (1,), device=binned_distribution_1.device
    )

    binned_distribution_1 = torch.cat([zero, binned_distribution_1], dim=-1)
    binned_distribution_2 = torch.cat([zero, binned_distribution_2], dim=-1)

    cumulative_1 = torch.cumsum(binned_distribution_1, dim=-1)
    cumulative_2 = torch.cumsum(binned_distribution_2, dim=-1)

    for i in range(K):
        relative_distribution[..., K - 1 + i] = torch.sum(
            torch.relu(
                torch.minimum(
                    cumulative_1[..., 1 + i :], cumulative_2[..., 1 : K + 1 - i]
                )
                - torch.maximum(cumulative_1[..., i:-1], cumulative_2[..., : K - i]),
            )
        )

    for i in range(1, K):
        relative_distribution[..., K - 1 - i] = torch.sum(
            torch.relu(
                torch.minimum(
                    cumulative_2[..., 1 + i :], cumulative_1[..., 1 : K + 1 - i]
                )
                - torch.maximum(cumulative_2[..., i:-1], cumulative_1[..., : K - i]),
            )
        )

    return relative_distribution
