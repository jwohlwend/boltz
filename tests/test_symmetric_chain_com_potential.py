import torch
import pytest

from boltz.model.potentials.potentials import SymmetricChainCOMPotential


@pytest.mark.parametrize("batch_size", [1, 2])
def test_symmetric_chain_com_potential_gradient_matches_autograd(batch_size):
    num_atoms = 9
    # 2 tokens/chains: chain 0 (3 atoms), chain 1 (4 atoms), 2 padding atoms
    atom_to_token = torch.zeros(1, num_atoms, 2)
    atom_to_token[0, 0:3, 0] = 1.0
    atom_to_token[0, 3:7, 1] = 1.0

    atom_pad_mask = torch.tensor(
        [[True, True, True, True, True, True, True, False, False]]
    )

    feats = {
        "atom_to_token": atom_to_token,
        "asym_id": torch.tensor([[0, 1]]),
        "symmetric_chain_index": torch.tensor([[[0], [1]]]),
        "atom_pad_mask": atom_pad_mask,
    }

    torch.manual_seed(42)
    coords = torch.randn(batch_size, num_atoms, 3, requires_grad=True)
    parameters = {"buffer": 10.0}

    pot = SymmetricChainCOMPotential()

    energy = pot.compute(coords, feats, parameters)
    expected_grad = torch.autograd.grad(energy.sum(), coords)[0]

    actual_grad = pot.compute_gradient(coords.detach(), feats, parameters)

    torch.testing.assert_close(actual_grad, expected_grad, atol=1e-6, rtol=1e-5)
    # Ensure padded atoms have zero gradient
    assert torch.all(actual_grad[:, 7:] == 0)
