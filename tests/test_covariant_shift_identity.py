"""Manufactured shift controls for the shift-aware covariant divergence."""
from __future__ import annotations

import numpy as np

from tools.covariant_shift_identity import (
    covariant_divergence_from_terms,
    eulerian_to_coordinate_components,
)


def _flat_chart(n: int = 120, beta_amp: float = 0.035, length: float = 4.3):
    """Minkowski chart R=r+t*v(r), T=t, with a smooth odd radial shift."""
    r_max = 12.0
    dr = r_max / n
    r = (np.arange(n, dtype=np.float64) + 0.5) * dr
    beta = beta_amp * r * np.exp(-(r / length) ** 2)
    beta_r = beta_amp * np.exp(-(r / length) ** 2) * (
        1.0 - 2.0 * (r / length) ** 2
    )
    energy = np.full_like(r, 2.3)
    pressure = np.full_like(r, 0.7)
    zeros = np.zeros_like(r)
    components = eulerian_to_coordinate_components(
        energy, zeros, pressure, pressure,
        np.ones_like(r), beta, np.ones_like(r), r**2,
    )

    # Exact chart-derived derivatives. These are analytic controls, not the
    # production stencil. In particular dot(g_tr)=v*v', dot(g_rr)=2*v',
    # dot(g_thth)=2*r*v, and dot(sqrt(-g)T^t_t)=-E*(r^2*v'+2*r*v).
    d_t_density_t = -energy * (r**2 * beta_r + 2.0 * r * beta)
    d_t_density_r = zeros
    d_r_flux_t = (energy + pressure) * (2.0 * r * beta + r**2 * beta_r)
    d_r_flux_r = 2.0 * pressure * r
    d_t_metric = {
        "gtt": zeros,
        "gtr": beta * beta_r,
        "grr": 2.0 * beta_r,
        "gthth": 2.0 * r * beta,
    }
    d_r_metric = {
        "gtt": 2.0 * beta * beta_r,
        "gtr": beta_r,
        "grr": zeros,
        "gthth": 2.0 * r,
    }
    terms = {
        "d_t_density_t": d_t_density_t,
        "d_t_density_r": d_t_density_r,
        "d_r_flux_t": d_r_flux_t,
        "d_r_flux_r": d_r_flux_r,
        "d_t_metric": d_t_metric,
        "d_r_metric": d_r_metric,
    }
    return r, beta, beta_r, energy, pressure, components, terms


def _annulus_rms(values: np.ndarray, r: np.ndarray) -> float:
    mask = (r >= 2.0) & (r <= 10.0)
    return float(np.sqrt(np.mean(np.asarray(values)[mask] ** 2)))


def test_zero_shift_reduces_to_existing_issue4_component_contract():
    n = 16
    x = np.linspace(0.1, 2.0, n)
    alpha = 0.9 + 0.03 * x
    grr = 1.1 + 0.04 * x
    gth = x**2 + 0.5
    energy = 0.2 + 0.01 * x
    j = -0.003 * x
    pr = 0.07 + 0.002 * x
    pt = 0.05 + 0.001 * x
    beta = np.zeros_like(x)
    c = eulerian_to_coordinate_components(
        energy, j, pr, pt, alpha, beta, grr, gth,
    )
    sqrt_gamma = np.sqrt(grr) * gth
    assert np.allclose(c["Ttt"], energy / alpha**2, rtol=0, atol=1e-15)
    assert np.allclose(c["Ttr"], j / (alpha * grr), rtol=0, atol=1e-15)
    assert np.allclose(c["Trr"], pr / grr, rtol=0, atol=1e-15)
    assert np.allclose(c["density_t"], -alpha * sqrt_gamma * energy, rtol=0, atol=1e-15)
    assert np.allclose(c["density_r"], sqrt_gamma * j, rtol=0, atol=1e-15)
    assert np.allclose(c["flux_t"], -alpha**2 * sqrt_gamma * j / grr, rtol=0, atol=1e-15)
    assert np.allclose(c["flux_r"], alpha * sqrt_gamma * pr, rtol=0, atol=1e-15)
    assert np.allclose(c["gtt"], -alpha**2, rtol=0, atol=1e-15)
    assert np.allclose(c["gtr"], 0.0, rtol=0, atol=1e-15)


def test_regular_nonzero_shift_minkowski_chart_closes_to_roundoff():
    r, beta, beta_r, energy, pressure, c, terms = _flat_chart()
    assert float(np.max(np.abs(beta))) > 0.05
    div_t, div_r = covariant_divergence_from_terms(c, **terms)
    assert float(np.max(np.abs(div_t))) < 2.0e-13
    assert float(np.max(np.abs(div_r))) < 2.0e-13
    assert _annulus_rms(div_t, r) < 2.0e-14
    assert _annulus_rms(div_r, r) < 2.0e-14


def test_shift_negative_controls_reject_omitted_cross_terms_and_wrong_flux_sign():
    r, beta, beta_r, energy, pressure, c, terms = _flat_chart()
    correct_t, correct_r = covariant_divergence_from_terms(c, **terms)
    wrong_flux_terms = dict(terms)
    wrong_flux_terms["d_r_flux_t"] = -terms["d_r_flux_t"]
    wrong_t, _ = covariant_divergence_from_terms(c, **wrong_flux_terms)
    assert _annulus_rms(wrong_t, r) > 1.0e-2

    no_radial_cross = dict(terms)
    no_radial_cross["d_r_metric"] = dict(terms["d_r_metric"])
    no_radial_cross["d_r_metric"]["gtr"] = np.zeros_like(r)
    _, wrong_r = covariant_divergence_from_terms(c, **no_radial_cross)
    assert _annulus_rms(wrong_r, r) > 1.0e-4

    no_time_cross = dict(terms)
    no_time_cross["d_t_metric"] = dict(terms["d_t_metric"])
    no_time_cross["d_t_metric"]["gtr"] = np.zeros_like(r)
    wrong_t_cross, _ = covariant_divergence_from_terms(c, **no_time_cross)
    assert _annulus_rms(wrong_t_cross, r) > 1.0e-5
    assert _annulus_rms(correct_t, r) < 2.0e-14
    assert _annulus_rms(correct_r, r) < 2.0e-14
