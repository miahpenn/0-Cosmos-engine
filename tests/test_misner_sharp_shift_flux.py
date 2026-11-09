"""Shift term of the Misner-Sharp ledger (engine.worldtube.misner_sharp_shift_flux).

d_t M|_r = 4 pi R^2 (T^r_t R' - T^r_r Rdot) at fixed coordinate r. The production ledger's
flux_T + work_pR is the beta = 0 part; misner_sharp_shift_flux is the rest.
"""
import math

import numpy as np
import pytest

from engine.worldtube import misner_sharp_shift_flux


def _ledger(R, Rr, Rdot, alpha, beta, gamma_rr, rho, j_r, S_rr):
    pr = S_rr / gamma_rr                      # mixed S^r_r, as production's total.pr
    flux_T = -4.0 * math.pi * R**2 * alpha * (j_r / gamma_rr) * Rr
    work_pR = -4.0 * math.pi * R**2 * pr * Rdot
    shift = misner_sharp_shift_flux(R, Rr, Rdot, alpha, beta, rho, pr, j_r)
    return flux_T + work_pR + shift


def _covariant(R, Rr, Rdot, alpha, beta, gamma_rr, rho, j_r, S_rr):
    """4 pi R^2 (T^r_t R' - T^r_r Rdot) from the full (t, r) metric and T_{mu nu}."""
    g = np.array([[-alpha**2 + gamma_rr * beta**2, gamma_rr * beta],
                  [gamma_rr * beta, gamma_rr]])
    n_low = np.array([-alpha, 0.0])
    j_low = np.array([beta * j_r, j_r])                    # j . n = 0
    S_low = S_rr * np.array([[beta**2, beta], [beta, 1.0]])  # S . n = 0
    T = rho * np.outer(n_low, n_low) + np.outer(n_low, j_low) + np.outer(j_low, n_low) + S_low
    T_up_low = np.linalg.inv(g) @ T                        # T^mu_nu
    return 4.0 * math.pi * R**2 * (T_up_low[1, 0] * Rr - T_up_low[1, 1] * Rdot)


def test_shift_term_completes_the_covariant_mass_balance():
    rng = np.random.default_rng(20261109)
    for _ in range(200):
        args = dict(
            R=rng.uniform(0.5, 20.0), Rr=rng.uniform(0.2, 2.0), Rdot=rng.normal(),
            alpha=rng.uniform(0.2, 2.0), beta=rng.normal(scale=0.7), gamma_rr=rng.uniform(0.3, 3.0),
            rho=rng.uniform(0.0, 2.0), j_r=rng.normal(), S_rr=rng.uniform(-0.5, 2.0),
        )
        assert _ledger(**args) == pytest.approx(_covariant(**args), rel=1e-12, abs=1e-12)


def test_shift_term_vanishes_without_shift():
    assert misner_sharp_shift_flux(3.0, 1.1, 0.4, 0.9, 0.0, 1.0, 0.3, -0.2) == 0.0


def test_flrw_in_a_stretched_radial_coordinate():
    """Flat FLRW with r' = r f(t): same slices (alpha = 1, j = 0), shift beta = -r' fdot/f.

    Exactly, M = (4 pi/3) rho R^3 with R = a r'/f, so d_t M|_{r'} = 4 pi R^3 (-H p - rho fdot/f).
    """
    a, H, rho, p, f, fdot, rp = 2.0, 0.3, 1.7, 0.4, 1.3, 0.25, 4.0
    R = a * rp / f
    Rr = a / f
    Rdot = R * (H - fdot / f)
    beta = -rp * fdot / f
    ledger = _ledger(R, Rr, Rdot, 1.0, beta, (a / f) ** 2, rho, 0.0, p * (a / f) ** 2)
    exact = 4.0 * math.pi * R**3 * (-H * p - rho * fdot / f)
    assert ledger == pytest.approx(exact, rel=1e-13)
    # Without the shift term the balance is wrong by 4 pi R^2 beta (rho + p) R'.
    flux_work = ledger - misner_sharp_shift_flux(R, Rr, Rdot, 1.0, beta, rho, p, 0.0)
    assert abs(flux_work - exact) > 0.1 * abs(exact)


def test_production_diagnostics_report_a_zero_shift_term():
    from engine.production_kernel import V55ProductionKernel

    kernel = V55ProductionKernel()
    state = kernel.initialize(resolution=24, r_max=40.0)
    assert kernel.diagnostics(state)["flux_shift"] == 0.0
