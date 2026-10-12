from types import SimpleNamespace

import numpy as np

from tools.audit_radiation_frozen_state_dt_stage_budget import metric_predictor_budget


def test_metric_predictor_budget_closes_and_checks_zero_shift_a_x_identities():
    n = 5
    grid = SimpleNamespace(n=n, centers=np.arange(n, dtype=float) + 0.5)
    a = np.full(n, 1.2)
    b = np.full(n, 0.9)
    X = np.full(n, 0.95)
    alpha = np.ones(n)
    Aa = np.full(n, 0.07)
    K = np.full(n, 0.12)
    accepted = SimpleNamespace(
        a=a.copy(), b=b.copy(), X=X.copy(), alpha=alpha.copy(),
        Aa=Aa.copy(), K=K.copy(),
    )
    rates = {
        "a": -2.0 * alpha * a * Aa,
        "b": np.full(n, 0.02),
        "X": alpha * X * K / 3.0,
    }
    dt = 0.003
    predicted = SimpleNamespace(
        a=a + dt * rates["a"],
        b=b + dt * rates["b"],
        X=X + dt * rates["X"],
    )

    report = metric_predictor_budget(grid, accepted, predicted, rates, dt)

    assert report["max_abs_gamma_identity_closure"] < 1.0e-14
    assert report["max_abs_rhs_a_identity_residual"] == 0.0
    assert report["max_abs_rhs_X_identity_residual"] == 0.0
    assert len(report["outer_cells"]) == 5
