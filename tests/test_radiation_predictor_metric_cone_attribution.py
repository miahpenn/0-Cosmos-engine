"""Captured-state regression for the radiation predictor's metric-induced cone crossing.

Numbers are from the source-pinned radiation trace artifact, workflow run #7.
This is an algebraic diagnostic test, not full trajectory validation.
"""
import math

import pytest

from engine.valencia import spherical_metric_from_bssn


def test_terminal_predictor_cone_crossing_is_metric_inverse_change():
    # Cell i=79, r=79.5: accepted metric then explicit predictor metric.
    accepted = spherical_metric_from_bssn(
        79.5, a=0.6089553401215667, b=1.281466561904975,
        X=0.9438956524107626, alpha=1.0, beta=0.0,
    )
    predictor = spherical_metric_from_bssn(
        79.5, a=0.6032035002114748, b=1.287561770875216,
        X=0.946558319470313, alpha=1.0, beta=0.0,
    )

    assert accepted.gamma_rr_inv == pytest.approx(1.4630613181946635)
    assert predictor.gamma_rr_inv == pytest.approx(1.4853571835116133)
    assert predictor.gamma_rr_inv > accepted.gamma_rr_inv

    # Hold the actual predictor conservative variables fixed. The ordered
    # trace budget attributes the final margin change to the metric update.
    U_E = 0.0002965253923984548
    U_r = 0.0002433162925154048
    C_on_accepted_metric = U_E - math.sqrt(accepted.gamma_rr_inv) * abs(U_r)
    C_on_predictor_metric = U_E - math.sqrt(predictor.gamma_rr_inv) * abs(U_r)

    assert C_on_accepted_metric == pytest.approx(
        2.2171301264504095e-6, abs=2e-15
    )
    assert C_on_predictor_metric == pytest.approx(
        -1.6900205171508503e-8, abs=2e-15
    )
    assert C_on_predictor_metric - C_on_accepted_metric == pytest.approx(
        -2.234030331621918e-6, abs=2e-15
    )

    # The realizability crossing is physical E < |S|, not a sign error in
    # the densitized coordinate momentum itself.
    physical_E = U_E / predictor.sqrt_gamma
    physical_abs_S = (
        math.sqrt(predictor.gamma_rr_inv) * abs(U_r) / predictor.sqrt_gamma
    )
    assert physical_abs_S > physical_E
    assert physical_abs_S / physical_E == pytest.approx(
        1.0000569941246342, abs=2e-12
    )
