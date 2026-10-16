"""Exact geometric controls for the independent covariant stress divergence."""
from __future__ import annotations

import importlib.util
import os
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("COSMOS_PINNED_SOURCE_ROOT", str(ROOT))
SPEC = importlib.util.spec_from_file_location(
    "_issue4_covariant_total_t_audit",
    ROOT / "tools" / "audit_total_stress_energy_divergence.py",
)
if SPEC is None or SPEC.loader is None:
    raise ImportError("cannot load independent total-T divergence auditor")
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)


def _d1(values: np.ndarray, parity: int) -> np.ndarray:
    """Polynomial-exact five-point derivative for the analytic test profiles."""
    del parity  # The test functions are explicit polynomials on the positive half-line.
    values = np.asarray(values, dtype=np.float64)
    n = values.size
    dr = 0.25
    if values.ndim != 1 or n < 5:
        raise ValueError("test derivative expects a one-dimensional profile with n>=5")
    result = np.empty_like(values)
    target = np.zeros(5, dtype=np.float64)
    target[1] = 1.0
    for i in range(n):
        start = min(max(i - 2, 0), n - 5)
        indices = np.arange(start, start + 5)
        offsets = indices.astype(np.float64) - float(i)
        vandermonde = np.vstack([offsets**power for power in range(5)])
        weights = np.linalg.solve(vandermonde, target)
        result[i] = float(weights @ values[indices]) / dr
    return result


def _zero(n: int) -> np.ndarray:
    return np.zeros(n, dtype=np.float64)


def test_flat_static_isotropic_stress_has_zero_covariant_divergence():
    n = 128
    dr = 0.25
    r = (np.arange(n, dtype=float) + 0.5) * dr
    alpha = np.ones(n)
    gamma_rr = np.ones(n)
    gamma_thth = r**2
    rho = np.full(n, 2.3)
    pressure = np.full(n, 0.7)
    j = _zero(n)
    zeros = _zero(n)
    div_t, div_r = AUDIT.covariant_divergence_from_profiles(
        r, _d1, alpha, gamma_rr, gamma_thth,
        rho, j, pressure, pressure,
        zeros, zeros, zeros, zeros, zeros,
    )
    assert float(np.max(np.abs(div_t))) < 2.0e-12
    assert float(np.max(np.abs(div_r))) < 2.0e-11


def test_flat_static_anisotropic_stress_matches_independent_radial_identity():
    n = 128
    dr = 0.25
    r = (np.arange(n, dtype=float) + 0.5) * dr
    alpha = np.ones(n)
    gamma_rr = np.ones(n)
    gamma_thth = r**2
    rho = np.full(n, 0.2)
    pr = 0.4 + 0.03 * r**2
    pt = 0.4 + 0.02 * r**2
    j = _zero(n)
    zeros = _zero(n)
    div_t, div_r = AUDIT.covariant_divergence_from_profiles(
        r, _d1, alpha, gamma_rr, gamma_thth,
        rho, j, pr, pt,
        zeros, zeros, zeros, zeros, zeros,
    )
    # In flat static spherical coordinates: D_r = p_r' + 2(p_r-p_t)/r.
    expected = 0.06 * r + 2.0 * (pr - pt) / r
    assert np.allclose(div_t, 0.0, rtol=0.0, atol=2.0e-12)
    assert np.allclose(div_r, expected, rtol=0.0, atol=2.0e-11)

def _d1_for_dr(dr: float):
    """Return a local polynomial-exact derivative for a chosen uniform spacing."""
    def derivative(values: np.ndarray, parity: int) -> np.ndarray:
        del parity  # The manufactured profile is evaluated on r>0.
        values = np.asarray(values, dtype=np.float64)
        n = values.size
        if values.ndim != 1 or n < 5:
            raise ValueError("test derivative expects a one-dimensional profile with n>=5")
        out = np.empty_like(values)
        target = np.zeros(5, dtype=np.float64)
        target[1] = 1.0
        for i in range(n):
            start = min(max(i - 2, 0), n - 5)
            indices = np.arange(start, start + 5)
            offsets = indices.astype(np.float64) - float(i)
            vandermonde = np.vstack([offsets**power for power in range(5)])
            weights = np.linalg.solve(vandermonde, target)
            out[i] = float(weights @ values[indices]) / dr
        return out
    return derivative


def _shifted_flat_minkowski_control(n: int) -> dict[str, float]:
    """Flat spacetime in T=a(t)+f(r), R=r coordinates; fluid is physically static.

    This control has nonzero radial shift, non-unit lapse, time-dependent g_tt
    and g_tr, and radial gradients of g_tr and gamma_rr. A constant physical
    perfect-fluid stress has exactly zero covariant divergence.
    """
    length = 20.0
    dr = length / n
    r = (np.arange(n, dtype=np.float64) + 0.5) * dr
    derivative = _d1_for_dr(dr)
    q, scale = 0.35, 3.0
    rho0, pressure0 = 2.3, 0.7
    a_dot, a_ddot = 1.2, 0.08

    # T=a(t)+f(r) applied to Minkowski time T and areal coordinate R=r.
    slope = q * r / np.sqrt(r**2 + scale**2)  # f'(r)
    gamma_rr = 1.0 - slope**2
    alpha = a_dot / np.sqrt(gamma_rr)
    beta_r = -a_dot * slope / gamma_rr
    gamma_thth = r**2

    # Eulerian projections of a constant perfect fluid at rest in (T,R).
    rho = (rho0 + pressure0) / gamma_rr - pressure0
    j = -(rho0 + pressure0) * slope / np.sqrt(gamma_rr)
    pr = pressure0 + (rho0 + pressure0) * slope**2 / gamma_rr
    pt = np.full_like(r, pressure0)
    sqrt_gamma = np.sqrt(gamma_rr) * gamma_thth
    sqrt_minus_g = alpha * sqrt_gamma

    # Exact time derivatives of the densitized mixed components and metric.
    dt_density_t = -rho0 * a_ddot * r**2
    dt_density_r = np.zeros_like(r)
    dt_gtt = np.full_like(r, -2.0 * a_dot * a_ddot)
    dt_gtr = -a_ddot * slope
    dt_gamma_rr = np.zeros_like(r)
    dt_gamma_thth = np.zeros_like(r)

    # Check the mixed-coordinate components against their exact transformed values.
    density_t = sqrt_minus_g * (-rho + beta_r * j / alpha)
    density_r = sqrt_minus_g * (j / alpha)
    flux_t = sqrt_minus_g * (
        beta_r * (rho + pr) - alpha * (j / gamma_rr) - beta_r**2 * j / alpha
    )
    flux_r = sqrt_minus_g * (pr - beta_r * j / alpha)
    assert np.allclose(density_t, -rho0 * a_dot * r**2, rtol=0.0, atol=2e-11)
    assert np.allclose(density_r, -(rho0 + pressure0) * slope * r**2, rtol=0.0, atol=2e-11)
    assert float(np.max(np.abs(flux_t))) < 2e-11
    assert np.allclose(flux_r, pressure0 * a_dot * r**2, rtol=0.0, atol=2e-11)
    assert float(np.max(np.abs(beta_r))) > 0.4
    assert float(np.max(np.abs(alpha - 1.0))) > 0.2

    div_t, div_r = AUDIT.covariant_divergence_from_profiles(
        r, derivative, alpha, gamma_rr, gamma_thth,
        rho, j, pr, pt,
        dt_density_t, dt_density_r, dt_gtt, dt_gamma_rr, dt_gamma_thth,
        beta_r=beta_r, dt_gtr=dt_gtr,
    )
    annulus = (r >= 2.0) & (r <= 10.0)
    return {
        "div_t_rms": float(np.sqrt(np.mean(div_t[annulus] ** 2))),
        "div_t_max": float(np.max(np.abs(div_t[annulus]))),
        "div_r_rms": float(np.sqrt(np.mean(div_r[annulus] ** 2))),
        "div_r_max": float(np.max(np.abs(div_r[annulus]))),
    }


def test_nonzero_shift_flat_minkowski_control_closes_and_converges():
    coarse, medium, fine = (
        _shifted_flat_minkowski_control(n) for n in (80, 160, 320)
    )
    assert coarse["div_t_max"] < 2e-11
    assert medium["div_t_max"] < 2e-11
    assert fine["div_t_max"] < 2e-11

    # D_r contains fourth-order radial metric and flux derivatives. Require the
    # manufactured residual to reduce at better than 3.5th order on refinements.
    order_1 = np.log(coarse["div_r_rms"] / medium["div_r_rms"]) / np.log(2.0)
    order_2 = np.log(medium["div_r_rms"] / fine["div_r_rms"]) / np.log(2.0)
    assert 3.5 < order_1 < 4.5
    assert 3.5 < order_2 < 4.5
    assert fine["div_r_rms"] < 5e-10

