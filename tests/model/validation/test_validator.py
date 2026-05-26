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

"""Tests for Validator metric logging in common_on_epoch_end.

Validates that all metrics computed and accumulated during validation_step
are properly read, reset, and logged when the epoch ends. The key concern
is that confidence-ranked lDDT, PDE MAE, and PAE MAE metrics were being
accumulated but never logged — a silent data loss bug.
"""

from __future__ import annotations

from unittest.mock import MagicMock

import pytest
import torch

from boltz.data import const
from boltz.model.validation.validator import Validator

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Pair-type keys returned by factored_lddt_loss / compute_pde_mae / compute_pae_mae.
# These are const.out_types minus "modified".
PAIR_METRIC_KEYS = [m for m in const.out_types if m != "modified"]

# Folding metric iteration keys (includes pocket/contact variants).
FOLDING_METRIC_KEYS = [*const.out_types, "pocket_ligand_protein", "contact_protein_protein"]

CONFIDENCE_PREFIXES = [
    "top1",
    "iplddt_top1",
    "ipde_top1",
    "pde_top1",
    "ptm_top1",
    "iptm_top1",
    "ligand_iptm_top1",
    "protein_iptm_top1",
    "avg",
]


def _make_mock_model(confidence_prediction: bool = True) -> tuple[MagicMock, dict[str, float]]:
    """Create a mock LightningModule that records log() calls.

    Returns
    -------
    model : MagicMock
        Mock with `.log()` wired to record calls.
    logged : dict
        Mapping from metric name to logged value, populated by `model.log()`.
    """
    model = MagicMock()
    model.confidence_prediction = confidence_prediction
    logged: dict[str, float] = {}

    def _log(name: str, value, **kwargs):
        logged[name] = value

    model.log = _log
    return model, logged


def _make_validator(confidence_prediction: bool = True) -> Validator:
    """Create a Validator with a single 'RCSB' validation dataset."""
    return Validator(
        val_names=["RCSB"],
        confidence_prediction=confidence_prediction,
        physicalism_metrics=False,
    )


def _populate_folding_metrics(validator: Validator) -> None:
    """Feed known values into all folding MeanMetrics so compute() is non-NaN."""
    idx = 0
    for m_ in FOLDING_METRIC_KEYS:
        validator.folding_metrics["lddt"][idx][m_].update(torch.tensor(0.5), torch.tensor(1.0))
        validator.folding_metrics["disto_lddt"][idx][m_].update(torch.tensor(0.4), torch.tensor(1.0))
        validator.folding_metrics["complex_lddt"][idx][m_].update(torch.tensor(0.6), torch.tensor(1.0))
    validator.folding_metrics["disto_loss"][idx]["disto_loss"].update(torch.tensor(0.1))


def _populate_confidence_metrics(validator: Validator) -> None:
    """Feed known values into all confidence MeanMetrics that exist.

    Populates plddt_mae (always initialized) and — if they have been
    initialized — the confidence-ranked lDDT, pde_mae, and pae_mae metrics.
    """
    idx = 0
    # plddt_mae — always initialized
    for m in const.out_single_types:
        validator.confidence_metrics["plddt_mae"][idx][m].update(torch.tensor(0.05), torch.tensor(1.0))

    # Confidence-ranked lDDT metrics
    for prefix in CONFIDENCE_PREFIXES:
        label = f"{prefix}_lddt"
        for key in PAIR_METRIC_KEYS:
            if key in validator.confidence_metrics[label][idx]:
                validator.confidence_metrics[label][idx][key].update(torch.tensor(0.45), torch.tensor(1.0))

    # PDE MAE and PAE MAE
    for mae_label in ["pde_mae", "pae_mae"]:
        for key in PAIR_METRIC_KEYS:
            if key in validator.confidence_metrics[mae_label][idx]:
                validator.confidence_metrics[mae_label][idx][key].update(torch.tensor(0.15), torch.tensor(1.0))


# ---------------------------------------------------------------------------
# Expected metric names
# ---------------------------------------------------------------------------


def _expected_folding_metric_names() -> set[str]:
    """Return the set of metric names that folding metrics should produce."""
    names: set[str] = set()
    for m_ in FOLDING_METRIC_KEYS:
        names.add(f"val/lddt_{m_}")
        names.add(f"val/disto_lddt_{m_}")
        names.add(f"val/complex_lddt_{m_}")
    names.add("val/disto_loss")
    names.add("val/disto_lddt")
    names.add("val/lddt")
    names.add("val/complex_lddt")
    return names


def _expected_plddt_mae_metric_names() -> set[str]:
    """Return the set of metric names for plddt_mae."""
    return {f"val/MAE_plddt_{m}" for m in const.out_single_types}


def _expected_confidence_lddt_metric_names() -> set[str]:
    """Return the set of metric names for confidence-ranked lDDT."""
    names: set[str] = set()
    for prefix in CONFIDENCE_PREFIXES:
        for key in PAIR_METRIC_KEYS:
            names.add(f"val/{prefix}_lddt_{key}")
    return names


def _expected_pde_mae_metric_names() -> set[str]:
    """Return the set of metric names for PDE MAE."""
    return {f"val/MAE_pde_{key}" for key in PAIR_METRIC_KEYS}


def _expected_pae_mae_metric_names() -> set[str]:
    """Return the set of metric names for PAE MAE."""
    return {f"val/MAE_pae_{key}" for key in PAIR_METRIC_KEYS}


def _all_expected_confidence_metric_names() -> set[str]:
    """All confidence metrics that should be logged (union of plddt/pde/pae/lddt)."""
    return (
        _expected_plddt_mae_metric_names()
        | _expected_confidence_lddt_metric_names()
        | _expected_pde_mae_metric_names()
        | _expected_pae_mae_metric_names()
    )


# ---------------------------------------------------------------------------
# Tests — Folding metrics (should pass before fix)
# ---------------------------------------------------------------------------


def test_folding_metrics_logged_without_confidence():
    validator = _make_validator(confidence_prediction=False)
    model, logged = _make_mock_model(confidence_prediction=False)
    _populate_folding_metrics(validator)

    validator.common_on_epoch_end(model)

    expected = _expected_folding_metric_names()
    missing = expected - set(logged.keys())
    assert not missing, f"Missing folding metrics: {sorted(missing)}"


def test_folding_metrics_have_correct_values():
    validator = _make_validator(confidence_prediction=False)
    model, logged = _make_mock_model(confidence_prediction=False)
    _populate_folding_metrics(validator)

    validator.common_on_epoch_end(model)

    assert logged["val/lddt_dna_protein"] == pytest.approx(0.5)
    assert logged["val/disto_lddt_dna_protein"] == pytest.approx(0.4)
    assert logged["val/complex_lddt_dna_protein"] == pytest.approx(0.6)
    assert logged["val/disto_loss"] == pytest.approx(0.1, abs=1e-5)


def test_no_confidence_metrics_when_disabled():
    """When confidence_prediction=False, no confidence metrics should be logged."""
    validator = _make_validator(confidence_prediction=False)
    model, logged = _make_mock_model(confidence_prediction=False)
    _populate_folding_metrics(validator)

    validator.common_on_epoch_end(model)

    confidence_names = {"MAE_plddt", "MAE_pde", "MAE_pae", "top1_lddt", "avg_lddt"}
    for name in logged:
        for cn in confidence_names:
            assert cn not in name, f"Unexpected confidence metric logged: {name}"


# ---------------------------------------------------------------------------
# Tests — Confidence metrics (main TDD tests — should FAIL before fix)
# ---------------------------------------------------------------------------


def test_plddt_mae_logged():
    """plddt_mae should already be logged (pre-existing behavior)."""
    validator = _make_validator(confidence_prediction=True)
    model, logged = _make_mock_model(confidence_prediction=True)
    _populate_folding_metrics(validator)
    _populate_confidence_metrics(validator)

    validator.common_on_epoch_end(model)

    expected = _expected_plddt_mae_metric_names()
    missing = expected - set(logged.keys())
    assert not missing, f"Missing plddt_mae metrics: {sorted(missing)}"


def test_confidence_lddt_metrics_initialized():
    """Confidence-ranked lDDT MeanMetrics must exist in the ModuleDict."""
    validator = _make_validator(confidence_prediction=True)

    for prefix in CONFIDENCE_PREFIXES:
        label = f"{prefix}_lddt"
        for key in PAIR_METRIC_KEYS:
            assert (
                key in validator.confidence_metrics[label][0]
            ), f"MeanMetric not initialized for confidence_metrics['{label}'][0]['{key}']"


def test_pde_mae_metrics_initialized():
    """PDE MAE MeanMetrics must exist in the ModuleDict."""
    validator = _make_validator(confidence_prediction=True)
    for key in PAIR_METRIC_KEYS:
        assert (
            key in validator.confidence_metrics["pde_mae"][0]
        ), f"MeanMetric not initialized for confidence_metrics['pde_mae'][0]['{key}']"


def test_pae_mae_metrics_initialized():
    """PAE MAE MeanMetrics must exist in the ModuleDict."""
    validator = _make_validator(confidence_prediction=True)
    for key in PAIR_METRIC_KEYS:
        assert (
            key in validator.confidence_metrics["pae_mae"][0]
        ), f"MeanMetric not initialized for confidence_metrics['pae_mae'][0]['{key}']"


def test_confidence_lddt_metrics_logged():
    """All confidence-ranked lDDT metrics must appear in logged output."""
    validator = _make_validator(confidence_prediction=True)
    model, logged = _make_mock_model(confidence_prediction=True)
    _populate_folding_metrics(validator)
    _populate_confidence_metrics(validator)

    validator.common_on_epoch_end(model)

    expected = _expected_confidence_lddt_metric_names()
    missing = expected - set(logged.keys())
    assert not missing, f"Missing confidence lDDT metrics: {sorted(missing)}"


def test_pde_mae_logged():
    """PDE MAE metrics must appear in logged output."""
    validator = _make_validator(confidence_prediction=True)
    model, logged = _make_mock_model(confidence_prediction=True)
    _populate_folding_metrics(validator)
    _populate_confidence_metrics(validator)

    validator.common_on_epoch_end(model)

    expected = _expected_pde_mae_metric_names()
    missing = expected - set(logged.keys())
    assert not missing, f"Missing PDE MAE metrics: {sorted(missing)}"


def test_pae_mae_logged():
    """PAE MAE metrics must appear in logged output."""
    validator = _make_validator(confidence_prediction=True)
    model, logged = _make_mock_model(confidence_prediction=True)
    _populate_folding_metrics(validator)
    _populate_confidence_metrics(validator)

    validator.common_on_epoch_end(model)

    expected = _expected_pae_mae_metric_names()
    missing = expected - set(logged.keys())
    assert not missing, f"Missing PAE MAE metrics: {sorted(missing)}"


def test_all_confidence_metrics_logged():
    """Comprehensive check: every confidence metric must be logged."""
    validator = _make_validator(confidence_prediction=True)
    model, logged = _make_mock_model(confidence_prediction=True)
    _populate_folding_metrics(validator)
    _populate_confidence_metrics(validator)

    validator.common_on_epoch_end(model)

    expected = _all_expected_confidence_metric_names()
    missing = expected - set(logged.keys())
    assert not missing, f"Missing confidence metrics ({len(missing)}): {sorted(missing)}"


def test_no_metric_name_collisions():
    """All logged metric names must be unique (no overwrites)."""
    validator = _make_validator(confidence_prediction=True)
    logged_names: list[str] = []

    def _log(name: str, value, **kwargs):
        logged_names.append(name)

    model = MagicMock()
    model.confidence_prediction = True
    model.log = _log

    _populate_folding_metrics(validator)
    _populate_confidence_metrics(validator)
    validator.common_on_epoch_end(model)

    assert len(logged_names) == len(
        set(logged_names)
    ), f"Duplicate metric names: {[n for n in logged_names if logged_names.count(n) > 1]}"


# ---------------------------------------------------------------------------
# Tests — Metric reset lifecycle
# ---------------------------------------------------------------------------


def test_folding_metrics_reset_after_epoch_end():
    """After common_on_epoch_end, folding MeanMetrics should be reset."""
    validator = _make_validator(confidence_prediction=False)
    model, _ = _make_mock_model(confidence_prediction=False)
    _populate_folding_metrics(validator)

    validator.common_on_epoch_end(model)

    for m_ in FOLDING_METRIC_KEYS:
        val = validator.folding_metrics["lddt"][0][m_].compute()
        assert torch.isnan(val), f"lddt[{m_}] not reset: {val}"


def test_confidence_metrics_reset_after_epoch_end():
    """After common_on_epoch_end, confidence MeanMetrics should be reset."""
    validator = _make_validator(confidence_prediction=True)
    model, _ = _make_mock_model(confidence_prediction=True)
    _populate_folding_metrics(validator)
    _populate_confidence_metrics(validator)

    validator.common_on_epoch_end(model)

    for m in const.out_single_types:
        val = validator.confidence_metrics["plddt_mae"][0][m].compute()
        assert torch.isnan(val), f"plddt_mae[{m}] not reset: {val}"

    for mae_label in ["pde_mae", "pae_mae"]:
        for key in PAIR_METRIC_KEYS:
            if key in validator.confidence_metrics[mae_label][0]:
                val = validator.confidence_metrics[mae_label][0][key].compute()
                assert torch.isnan(val), f"{mae_label}[{key}] not reset: {val}"

    for prefix in CONFIDENCE_PREFIXES:
        label = f"{prefix}_lddt"
        for key in PAIR_METRIC_KEYS:
            if key in validator.confidence_metrics[label][0]:
                val = validator.confidence_metrics[label][0][key].compute()
                assert torch.isnan(val), f"{label}[{key}] not reset: {val}"


def test_two_epoch_cycle():
    """Simulate two validation epochs and verify independent accumulation.

    This tests the core lifecycle:
    1. Epoch 1: populate -> epoch_end -> log values -> reset
    2. Epoch 2: populate with different values -> epoch_end -> log new values
    """
    validator = _make_validator(confidence_prediction=False)

    # --- Epoch 1 ---
    model1, logged1 = _make_mock_model(confidence_prediction=False)
    _populate_folding_metrics(validator)
    validator.common_on_epoch_end(model1)
    epoch1_lddt = logged1["val/lddt_dna_protein"]

    # --- Epoch 2: different values ---
    model2, logged2 = _make_mock_model(confidence_prediction=False)
    for m_ in FOLDING_METRIC_KEYS:
        validator.folding_metrics["lddt"][0][m_].update(torch.tensor(0.9), torch.tensor(1.0))
        validator.folding_metrics["disto_lddt"][0][m_].update(torch.tensor(0.8), torch.tensor(1.0))
        validator.folding_metrics["complex_lddt"][0][m_].update(torch.tensor(0.7), torch.tensor(1.0))
    validator.folding_metrics["disto_loss"][0]["disto_loss"].update(torch.tensor(0.01))
    validator.common_on_epoch_end(model2)
    epoch2_lddt = logged2["val/lddt_dna_protein"]

    assert epoch1_lddt == pytest.approx(0.5)
    assert epoch2_lddt == pytest.approx(0.9)
    assert epoch1_lddt != epoch2_lddt, "Epoch 2 should reflect new values, not stale ones"


def test_unpopulated_metrics_log_zero():
    """If no validation batches run (e.g. empty val set), metrics should log 0.0.

    This is relevant for resume scenarios where the validator is freshly
    created but no val batches have run yet.
    """
    validator = _make_validator(confidence_prediction=False)
    model, logged = _make_mock_model(confidence_prediction=False)

    validator.common_on_epoch_end(model)

    for m_ in FOLDING_METRIC_KEYS:
        assert logged[f"val/lddt_{m_}"] == 0.0, f"Unpopulated lddt_{m_} should be 0.0 (NaN->0.0 fallback)"


# ---------------------------------------------------------------------------
# Tests — Dataset name suffix
# ---------------------------------------------------------------------------


def test_rcsb_has_no_suffix():
    validator = _make_validator(confidence_prediction=False)
    model, logged = _make_mock_model(confidence_prediction=False)
    _populate_folding_metrics(validator)
    validator.common_on_epoch_end(model)

    assert "val/lddt_dna_protein" in logged


def test_non_rcsb_has_suffix():
    validator = Validator(
        val_names=["CUSTOM"],
        confidence_prediction=False,
        physicalism_metrics=False,
    )
    model, logged = _make_mock_model(confidence_prediction=False)
    _populate_folding_metrics(validator)
    validator.common_on_epoch_end(model)

    assert "val/lddt_dna_protein__CUSTOM" in logged
    assert "val/lddt_dna_protein" not in logged
