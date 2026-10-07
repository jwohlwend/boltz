import json
import math

import pytest
import torch

from boltz.data.write.writer import _nan_to_none
from boltz.model.layers.confidence_utils import compute_pae_summaries


def _feats(asym_ids, pad_mask=None):
    asym_id = torch.tensor([asym_ids])
    if pad_mask is None:
        pad_mask = torch.ones_like(asym_id)
    else:
        pad_mask = torch.tensor([pad_mask])
    return {"asym_id": asym_id, "token_pad_mask": pad_mask}


def _block_pae(asym_ids, values):
    """PAE matrix whose value depends only on the (aligned, scored) chain pair."""
    n = len(asym_ids)
    pae = torch.zeros(1, n, n)
    for i in range(n):
        for j in range(n):
            pae[0, i, j] = values[(asym_ids[i], asym_ids[j])]
    return pae


def test_two_chain_block_means():
    asym = [0, 0, 0, 1, 1]
    values = {(0, 0): 1.0, (0, 1): 8.0, (1, 0): 6.0, (1, 1): 2.0}
    pae = _block_pae(asym, values)
    complex_pae, complex_ipae, pair = compute_pae_summaries(pae, _feats(asym), 1)

    # Off-diagonal pairs: 3*2 intra chain 0, 2*1 intra chain 1, 6 + 6 inter.
    expected_complex = (6 * 1.0 + 2 * 2.0 + 6 * 8.0 + 6 * 6.0) / 20
    assert complex_pae.shape == (1,)
    assert torch.isclose(complex_pae[0], torch.tensor(expected_complex))
    assert torch.isclose(complex_ipae[0], torch.tensor((6 * 8.0 + 6 * 6.0) / 12))
    assert torch.isclose(pair[0][0][0], torch.tensor(1.0))
    assert torch.isclose(pair[0][1][0], torch.tensor(8.0))
    assert torch.isclose(pair[1][0][0], torch.tensor(6.0))
    assert torch.isclose(pair[1][1][0], torch.tensor(2.0))


def test_diagonal_excluded():
    asym = [0, 0]
    pae = torch.tensor([[[100.0, 3.0], [5.0, 100.0]]])
    complex_pae, _, pair = compute_pae_summaries(pae, _feats(asym), 1)
    assert torch.isclose(complex_pae[0], torch.tensor(4.0))
    assert torch.isclose(pair[0][0][0], torch.tensor(4.0))


def test_single_chain_interface_is_nan():
    asym = [0, 0, 0]
    pae = torch.full((1, 3, 3), 2.0)
    complex_pae, complex_ipae, pair = compute_pae_summaries(pae, _feats(asym), 1)
    assert torch.isclose(complex_pae[0], torch.tensor(2.0))
    assert torch.isnan(complex_ipae[0])
    assert list(pair) == [0]


def test_padding_and_frame_mask():
    asym = [0, 0, 1, 1]
    pad = [1, 1, 1, 0]  # last token is padding
    values = {(0, 0): 1.0, (0, 1): 4.0, (1, 0): 9.0, (1, 1): 2.0}
    pae = _block_pae(asym, values)
    feats = _feats(asym, pad)

    _, complex_ipae, pair = compute_pae_summaries(pae, feats, 1)
    # Valid inter-chain pairs: (0,2), (1,2) -> 4.0 ; (2,0), (2,1) -> 9.0
    assert torch.isclose(complex_ipae[0], torch.tensor((2 * 4.0 + 2 * 9.0) / 4))
    # Chain 1 has a single valid token -> no intra-chain pair.
    assert torch.isnan(pair[1][1][0])

    # Masking the frame of token 2 removes every row aligned on chain 1.
    frame_mask = torch.tensor([[1, 1, 0, 1]])
    _, complex_ipae, pair = compute_pae_summaries(pae, feats, 1, frame_mask=frame_mask)
    assert torch.isclose(complex_ipae[0], torch.tensor(4.0))
    assert torch.isnan(pair[1][0][0])
    assert torch.isclose(pair[0][1][0], torch.tensor(4.0))


def test_multiplicity_batches_are_independent():
    asym = [0, 1]
    feats = _feats(asym)
    pae = torch.stack(
        [
            torch.tensor([[0.0, 2.0], [4.0, 0.0]]),
            torch.tensor([[0.0, 10.0], [20.0, 0.0]]),
        ]
    )
    complex_pae, complex_ipae, pair = compute_pae_summaries(pae, feats, 2)
    assert complex_pae.shape == (2,)
    assert torch.allclose(complex_pae, torch.tensor([3.0, 15.0]))
    assert torch.allclose(complex_ipae, torch.tensor([3.0, 15.0]))
    assert torch.allclose(pair[0][1], torch.tensor([2.0, 10.0]))
    assert torch.allclose(pair[1][0], torch.tensor([4.0, 20.0]))


def test_nan_serializes_to_null():
    assert _nan_to_none(1.5) == 1.5
    assert _nan_to_none(float("nan")) is None
    payload = json.dumps({"complex_ipae": _nan_to_none(float("nan"))})
    assert json.loads(payload)["complex_ipae"] is None
    with pytest.raises(ValueError):
        json.loads('{"x": NaN}', parse_constant=_reject_constant)


def _reject_constant(name):
    msg = f"invalid JSON constant {name}"
    raise ValueError(msg)


def test_frame_mask_without_padding_matches_plain_mean():
    asym = [0, 0, 1]
    pae = torch.rand(1, 3, 3) * 30
    feats = _feats(asym)
    a, _, _ = compute_pae_summaries(pae, feats, 1)
    b, _, _ = compute_pae_summaries(pae, feats, 1, frame_mask=torch.ones(1, 3))
    off_diag = pae[0][~torch.eye(3, dtype=torch.bool)]
    assert torch.isclose(a[0], off_diag.mean())
    assert torch.isclose(a[0], b[0])
    assert not math.isnan(a[0].item())
