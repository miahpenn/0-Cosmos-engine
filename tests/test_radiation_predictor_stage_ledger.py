"""Deterministic algebraic checks for predictor-stage attribution."""
import math

import pytest

from engine.radiation_predictor_stage_ledger import reconstruct_metric_cone, validate_ledger_steps


def test_reconstruct_metric_inverse_and_cone_margin_independently():
    values = reconstruct_metric_cone(
        a=0.6089553401215667,
        X=0.9438956524107626,
        U_E=0.0002965253923984548,
        U_r=0.0002433162925154048,
    )
    assert values["gamma_rr_inv_reconstructed"] == pytest.approx(
        1.4630613181946635, abs=2e-15
    )
    assert values["cone_margin_reconstructed"] == pytest.approx(
        2.2171301264504095e-6, abs=2e-15
    )
    assert values["ratio_absS_over_E_reconstructed"] < 1.0


def test_reconstruction_reproduces_captured_predictor_crossing():
    values = reconstruct_metric_cone(
        a=0.6032035002114748,
        X=0.946558319470313,
        U_E=0.0002965253923984548,
        U_r=0.0002433162925154048,
    )
    assert values["gamma_rr_inv_reconstructed"] == pytest.approx(
        1.4853571835116133, abs=2e-15
    )
    assert values["cone_margin_reconstructed"] == pytest.approx(
        -1.6900205171508503e-8, abs=2e-15
    )
    assert values["ratio_absS_over_E_reconstructed"] == pytest.approx(
        1.0000569941246342, abs=2e-12
    )


@pytest.mark.parametrize("a", [0.0, -1.0, math.nan, math.inf])
def test_reconstruction_rejects_invalid_radial_metric(a):
    with pytest.raises(ValueError):
        reconstruct_metric_cone(a, 1.0, 1.0, 0.0)


def test_empty_ledger_fails_closed():
    result = validate_ledger_steps([], range(75, 80))
    assert result["ok"] is False
    assert result["step_count"] == 0
    assert result["errors"]


def test_invalid_inverse_metric_factor_rejected():
    with pytest.raises(ValueError):
        reconstruct_metric_cone(1.0, 0.0, 1.0, 0.1)


def test_nonpositive_radiation_energy_rejected():
    with pytest.raises(ValueError):
        reconstruct_metric_cone(1.0, 1.0, 0.0, 0.1)
