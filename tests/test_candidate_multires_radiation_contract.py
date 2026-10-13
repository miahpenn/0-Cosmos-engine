"""Strict contract tests for the candidate multi-resolution radiation witnesses."""

import pytest

from tools.validate_initial_quadrature_cmc_candidate import (
    REQUIRED_RADIATION_DIAGNOSTICS,
    _required_radiation_summary,
)


def _finite_diagnostics():
    return {
        key: float(index + 1)
        for index, key in enumerate(REQUIRED_RADIATION_DIAGNOSTICS)
    }


def test_radiation_summary_preserves_every_required_witness():
    diagnostics = _finite_diagnostics()

    summary = _required_radiation_summary(diagnostics)

    assert tuple(summary) == REQUIRED_RADIATION_DIAGNOSTICS
    assert summary == diagnostics


def test_radiation_summary_rejects_a_missing_required_witness():
    diagnostics = _finite_diagnostics()
    missing = "radiation_outer_transport_to_source_S"
    del diagnostics[missing]

    with pytest.raises(KeyError, match=missing):
        _required_radiation_summary(diagnostics)


@pytest.mark.parametrize("bad_value", [float("nan"), float("inf"), -float("inf")])
def test_radiation_summary_rejects_nonfinite_witness_values(bad_value):
    diagnostics = _finite_diagnostics()
    diagnostics["radiation_outer_ratio"] = bad_value

    with pytest.raises(FloatingPointError, match="radiation_outer_ratio is non-finite"):
        _required_radiation_summary(diagnostics)
