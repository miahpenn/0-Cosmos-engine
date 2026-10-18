"""Regression for the conditional scalar–dark-matter exchange projection.

This is an implementation test of the archived conditional source law, not
an admission of the underlying microphysical coupling premise.
"""
import numpy as np

from engine.matter_system import (
    BSSNMetricSlice,
    MetricDerivativeSet,
    Species,
    _metric_arrays,
    _species_rhs_components,
    initialize_dust,
    primitives,
)


def test_dark_matter_rest_current_source_matches_nonzero_shift_projection():
    # Non-unit lapse and nonzero shift ensure the coordinate/normal derivative
    # conversion cannot accidentally pass as a flat, zero-shift identity.
    metric = BSSNMetricSlice(
        r=np.asarray([0.5, 1.5]),
        a=np.asarray([1.1, 1.2]),
        b=np.asarray([0.95, 1.05]),
        X=np.asarray([0.9, 1.0]),
        alpha=np.asarray([0.7, 0.85]),
        beta=np.asarray([0.12, -0.08]),
        Aa=np.zeros(2),
        K=np.zeros(2),
        Lambda=np.zeros(2),
        B=np.zeros(2),
    )
    zero = np.zeros(2)
    derivs = MetricDerivativeSet(
        time={},
        radial={
            "alpha": zero.copy(),
            "beta": zero.copy(),
            "rr": zero.copy(),
            "thth": zero.copy(),
        },
    )
    rho_initial = np.asarray([1.3, 0.8])
    velocity = np.asarray([0.21, -0.17])
    metrics = _metric_arrays(metric)
    state = initialize_dust(metrics, rho_initial, velocity)
    dphi_t = np.asarray([0.31, -0.22])
    dphi_r = np.asarray([-0.14, 0.27])
    beta_dm = -0.04

    (_flux, _geom_e, _geom_s, q_e, q_s, q_rest) = _species_rhs_components(
        metric, derivs, state, Species.DARK_MATTER,
        dphi_t=dphi_t, dphi_r=dphi_r, beta_dm=beta_dm,
    )
    recovered = primitives(metrics, state, Species.DARK_MATTER)

    expected_rest = np.empty(2)
    expected_e = np.empty(2)
    expected_s = np.empty(2)
    missing_shift = np.empty(2)
    for i, (m, p) in enumerate(zip(metrics, recovered)):
        # Independent contraction u^mu d_mu(phi) from coordinate components:
        # u^t=W/alpha, u^r=W(v^r-beta^r/alpha).
        W = p.lorentz()
        u_grad_phi = W / m.alpha * (
            dphi_t[i] + (m.alpha * p.v_r - m.beta) * dphi_r[i]
        )
        expected_rest[i] = (
            -m.alpha * m.sqrt_gamma * beta_dm * p.rho * u_grad_phi
        )
        q_t = beta_dm * p.rho * dphi_t[i]
        q_r = beta_dm * p.rho * dphi_r[i]
        expected_e[i] = -m.sqrt_gamma * (q_t - m.beta * q_r)
        expected_s[i] = m.alpha * m.sqrt_gamma * q_r
        missing_shift[i] = (
            -m.alpha * m.sqrt_gamma * beta_dm * p.rho
            * W / m.alpha * (dphi_t[i] + m.alpha * p.v_r * dphi_r[i])
        )

    assert np.all(np.abs(q_rest) > 0.0)
    np.testing.assert_allclose(q_rest, expected_rest, rtol=2e-14, atol=0.0)
    np.testing.assert_allclose(q_e, expected_e, rtol=2e-14, atol=0.0)
    np.testing.assert_allclose(q_s, expected_s, rtol=2e-14, atol=0.0)

    # Signed/factor negative controls: either error must be observable for
    # this deliberately non-degenerate fixture.
    assert not np.allclose(q_rest, -expected_rest, rtol=1e-12, atol=0.0)
    assert not np.allclose(q_rest, expected_rest / metric.alpha, rtol=1e-12, atol=0.0)
    assert not np.allclose(q_rest, missing_shift, rtol=1e-12, atol=0.0)
