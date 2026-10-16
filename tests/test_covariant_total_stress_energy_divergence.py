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


def test_finite_volume_spherical_face_flux_reconstruction_is_exact_for_quadratic_flux():
    n = 64
    dr = 0.25
    r = (np.arange(n, dtype=np.float64) + 0.5) * dr
    flux = 0.7 * r**2
    derivative = AUDIT.finite_volume_radial_flux_derivative(flux, dr)
    expected = 1.4 * r
    assert np.allclose(derivative, expected, rtol=0.0, atol=2.0e-12)


def test_flat_static_isotropic_stress_has_zero_finite_volume_covariant_divergence():
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
    div_t, div_r = AUDIT.covariant_divergence_fv_from_profiles(
        r, _d1, alpha, gamma_rr, gamma_thth,
        rho, j, pressure, pressure,
        zeros, zeros, zeros, zeros, zeros,
    )
    assert float(np.max(np.abs(div_t))) < 2.0e-12
    assert float(np.max(np.abs(div_r))) < 2.0e-11

def test_finite_volume_face_reconstruction_respects_regular_origin_parity():
    n = 32
    dr = 0.25
    r = (np.arange(n, dtype=np.float64) + 0.5) * dr

    # F_r = r^2 g(r^2), with quadratic g, is reconstructed exactly at r=dr.
    even_flux = 0.7 * r**2 * (1.0 + 0.03 * r**2 + 0.002 * r**4)
    even_faces = AUDIT.spherical_flux_face_values(even_flux, dr, origin_power=2)
    expected_even = 0.7 * dr**2 * (1.0 + 0.03 * dr**2 + 0.002 * dr**4)
    assert np.isclose(even_faces[1], expected_even, rtol=0.0, atol=2.0e-14)

    # F_t = r^3 g(r^2) uses the odd regular-origin power.
    odd_flux = 0.4 * r**3 * (1.0 + 0.02 * r**2 + 0.001 * r**4)
    odd_faces = AUDIT.spherical_flux_face_values(odd_flux, dr, origin_power=3)
    expected_odd = 0.4 * dr**3 * (1.0 + 0.02 * dr**2 + 0.001 * dr**4)
    assert np.isclose(odd_faces[1], expected_odd, rtol=0.0, atol=2.0e-14)

    # The radial metric-source density scales as r at a regular origin.
    source_density = 0.6 * r * (1.0 + 0.015 * r**2 + 0.0005 * r**4)
    source_faces = AUDIT.spherical_flux_face_values(source_density, dr, origin_power=1)
    expected_source = 0.6 * dr * (1.0 + 0.015 * dr**2 + 0.0005 * dr**4)
    assert np.isclose(source_faces[1], expected_source, rtol=0.0, atol=2.0e-14)

