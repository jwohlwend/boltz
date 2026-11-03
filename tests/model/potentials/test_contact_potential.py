"""Check manual contact guidance against the energy it differentiates."""

import os

import pytest
import torch

from boltz.model.potentials.potentials import ContactPotentital


DEVICE = os.environ.get("BOLTZ_TEST_DEVICE", "cpu")


def contact_features(union_index, negation_mask):
    return {
        "contact_pair_index": torch.tensor(
            [[[0, 0, 2], [1, 2, 3]]], device=DEVICE
        ),
        "contact_union_index": torch.tensor([union_index], device=DEVICE),
        "contact_negation_mask": torch.tensor([negation_mask], device=DEVICE),
        "contact_thresholds": torch.tensor([[0.7, 1.5, 1.0]], device=DEVICE),
    }


@pytest.mark.parametrize("union_lambda", [0.0, 0.5, 1.0, 3.0])
@pytest.mark.parametrize("union_index", [[0, 0, 0], [0, 0, 1], [0, 1, 2]])
@pytest.mark.parametrize("batch_size", [1, 3])
@pytest.mark.parametrize("negation_mask", [[True, True, True], [True, False, True]])
def test_contact_gradient_matches_autograd(
    union_lambda, union_index, batch_size, negation_mask
):
    """Cover unions, singleton controls, negation, and diffusion multiplicity."""
    coords = torch.tensor(
        [[0.2, 0.1, 0.3], [1.4, -0.4, 0.1], [-0.6, 2.2, 0.5], [0.8, 1.6, -1.1]],
        device=DEVICE,
    ).repeat(batch_size, 1, 1)
    # Different samples retain the same atom/constraint indexing.
    coords *= torch.linspace(0.9, 1.1, batch_size, device=DEVICE)[:, None, None]
    feats = contact_features(union_index, negation_mask)
    parameters = {"union_lambda": union_lambda}
    potential = ContactPotentital()
    with torch.enable_grad():
        coords.requires_grad_(True)
        energy = potential.compute(coords, feats, parameters)
        expected = torch.autograd.grad(energy.sum(), coords)[0]
    actual = potential.compute_gradient(coords.detach(), feats, parameters)
    torch.testing.assert_close(actual, expected, atol=2e-6, rtol=2e-5)


def flat_and_violated_contacts():
    coords = torch.tensor(
        [[[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [0.0, 3.0, 0.0], [3.0, 3.0, 0.0]]],
        device=DEVICE,
    )
    feats = {
        "contact_pair_index": torch.tensor([[[0, 2], [1, 3]]], device=DEVICE),
        "contact_union_index": torch.tensor([[0, 0]], device=DEVICE),
        "contact_negation_mask": torch.tensor([[True, True]], device=DEVICE),
        "contact_thresholds": torch.tensor([[1.0, 1.0]], device=DEVICE),
    }
    return coords, feats, {"union_lambda": 1.0}


def test_contact_gradient_matches_finite_difference():
    """An independent numerical derivative, away from flat-bottom corners."""
    coords, feats, parameters = flat_and_violated_contacts()
    potential = ContactPotentital()
    delta = torch.zeros_like(coords)
    delta[0, 3, 0] = 1e-3
    expected = (
        potential.compute(coords + delta, feats, parameters)
        - potential.compute(coords - delta, feats, parameters)
    ) / 2e-3
    actual = potential.compute_gradient(coords, feats, parameters)[0, 3, 0]
    torch.testing.assert_close(actual, expected.squeeze(), atol=2e-5, rtol=2e-3)


def test_contact_gradient_step_lowers_computed_energy():
    """One satisfied alternative must not invert the other energy derivative."""
    coords, feats, parameters = flat_and_violated_contacts()
    potential = ContactPotentital()
    before = potential.compute(coords, feats, parameters)
    gradient = potential.compute_gradient(coords, feats, parameters)
    after = potential.compute(coords - 1e-2 * gradient, feats, parameters)
    assert torch.all(after < before)


def test_empty_contact_gradient_is_zero():
    coords, feats, parameters = flat_and_violated_contacts()
    feats = {key: value[..., :0] for key, value in feats.items()}
    gradient = ContactPotentital().compute_gradient(coords, feats, parameters)
    torch.testing.assert_close(gradient, torch.zeros_like(coords))
