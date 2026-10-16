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
    del parity  # These analytic controls are smooth even profiles on the sampled positive half-line.
    return np.gradient(np.asarray(values, dtype=float), 0.25, edge_order=2)


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
    pt = 0.4 + 0.06 * r**2
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
