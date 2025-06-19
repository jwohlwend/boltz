import torch
from torch import Tensor


def bfactor_loss_fn(
    output: dict[str, Tensor],
    feats: dict[str, Tensor],
) -> Tensor:
    """Compute the bfactor loss.

    Parameters
    ----------
    output : dict[str, Tensor]
        Output of the model
    feats : dict[str, Tensor]
        Input features

    Returns
    -------
    Tensor
        The globally averaged loss.

    """
    with torch.autocast("cuda", enabled=False):
        # Get predicted distograms
        pred = output["pbfactor"].float()  # (B, L, bins)
        bins = pred.shape[2]  # num_bins
        token_to_rep_atom = feats["token_to_rep_atom"]

        # Compute target histogram
        bfactor_atom = feats["bfactor"].unsqueeze(-1)  # (B, L)
        bfactor_token = torch.bmm(token_to_rep_atom.float(), bfactor_atom)

        boundaries = torch.linspace(0, 100, bins - 1, device=bfactor_token.device)
        bfactor_token_bin = (bfactor_token > boundaries).sum(dim=-1).long()
        bfactor_target = torch.nn.functional.one_hot(
            bfactor_token_bin, num_classes=bins
        )

        # Combine target mask and padding mask
        token_mask = (bfactor_token > 1e-5).squeeze(-1).float()

        # Compute the bfactor loss
        errors = -1 * torch.sum(
            bfactor_target * torch.nn.functional.log_softmax(pred, dim=-1),
            dim=-1,
        )
        loss = torch.sum(errors * token_mask) / (torch.sum(token_mask) + 1e-5)
        return loss
