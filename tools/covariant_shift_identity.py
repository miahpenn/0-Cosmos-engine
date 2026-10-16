"""Shift-aware spherical covariant stress-energy divergence identities.

Diagnostic utility only. Inputs use Eulerian energy E, covariant radial current
J_r, and orthonormal radial/tangential pressures. Derivative arrays are
provided by the caller so analytic manufactured controls can use exact
derivatives and production audits can use their declared discrete operators.
"""
from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import numpy as np


_METRIC_KEYS = ("gtt", "gtr", "grr", "gthth")


def _profile(name: str, value: Any, shape: tuple[int, ...]) -> np.ndarray:
    out = np.asarray(value, dtype=np.float64)
    if out.shape != shape:
        raise ValueError(f"{name} has shape {out.shape}, expected {shape}")
    if not np.all(np.isfinite(out)):
        raise FloatingPointError(f"{name} contains non-finite values")
    return out


def eulerian_to_coordinate_components(
    energy: Any,
    j_r: Any,
    p_r: Any,
    p_t: Any,
    alpha: Any,
    beta_r: Any,
    gamma_rr: Any,
    gamma_thth: Any,
) -> dict[str, np.ndarray]:
    """Build contravariant/mixed T components and densitized mixed fluxes.

    The radial current convention matches the zero-shift Issue #4 auditor:
    at beta^r=0, T^t_r=J_r/alpha and T^r_t=-alpha*J^r. The metric is
    ds^2=(-alpha^2+gamma_rr*beta^2)dt^2
          +2 gamma_rr*beta dt dr + gamma_rr dr^2 + gamma_thth dOmega^2.

    Returned Ttt/Ttr/Trr/Tthth are contravariant components. density_t/r
    and flux_t/r are sqrt(-g)*T^t{}_nu and sqrt(-g)*T^r{}_nu respectively.
    """
    E = np.asarray(energy, dtype=np.float64)
    if E.ndim != 1 or E.size < 1:
        raise ValueError("energy must be a nonempty one-dimensional profile")
    shape = E.shape
    E = _profile("energy", E, shape)
    j = _profile("j_r", j_r, shape)
    pr = _profile("p_r", p_r, shape)
    pt = _profile("p_t", p_t, shape)
    a = _profile("lapse", alpha, shape)
    beta = _profile("radial shift", beta_r, shape)
    grr = _profile("gamma_rr", gamma_rr, shape)
    gth = _profile("gamma_thth", gamma_thth, shape)
    if np.any(a <= 0.0) or np.any(grr <= 0.0) or np.any(gth <= 0.0):
        raise ValueError("lapse and spatial metric components must be positive")

    sqrt_gamma = np.sqrt(grr) * gth
    sqrt_minus_g = a * sqrt_gamma

    Ttt = E / a**2
    Ttr = -beta * E / a**2 + j / (a * grr)
    Trr = E * beta**2 / a**2 - 2.0 * beta * j / (a * grr) + pr / grr
    Tthth = pt / gth

    # Mixed components obtained by lowering the second tensor index.
    mixed_t_t = -E + beta * j / a
    mixed_t_r = j / a
    mixed_r_t = beta * (E + pr) - a * j / grr - beta**2 * j / a
    mixed_r_r = pr - beta * j / a

    gtt = -a**2 + grr * beta**2
    gtr = grr * beta

    return {
        "sqrt_gamma": sqrt_gamma,
        "sqrt_minus_g": sqrt_minus_g,
        "Ttt": Ttt,
        "Ttr": Ttr,
        "Trr": Trr,
        "Tthth": Tthth,
        "mixed_t_t": mixed_t_t,
        "mixed_t_r": mixed_t_r,
        "mixed_r_t": mixed_r_t,
        "mixed_r_r": mixed_r_r,
        "density_t": sqrt_minus_g * mixed_t_t,
        "density_r": sqrt_minus_g * mixed_t_r,
        "flux_t": sqrt_minus_g * mixed_r_t,
        "flux_r": sqrt_minus_g * mixed_r_r,
        "gtt": gtt,
        "gtr": gtr,
        "grr": grr,
        "gthth": gth,
    }


def covariant_divergence_from_terms(
    components: Mapping[str, Any],
    *,
    d_t_density_t: Any,
    d_t_density_r: Any,
    d_r_flux_t: Any,
    d_r_flux_r: Any,
    d_t_metric: Mapping[str, Any],
    d_r_metric: Mapping[str, Any],
) -> tuple[np.ndarray, np.ndarray]:
    """Return (nabla_mu T^mu_t, nabla_mu T^mu_r) using supplied derivatives.

    Uses nabla_mu T^mu{}_nu =
      [d_t(sqrt(-g)T^t_nu)+d_r(sqrt(-g)T^r_nu)]/sqrt(-g)
      - 1/2 T^{ab} d_nu g_ab,
    including both g_tr cross terms and the spherical angular contribution.
    Derivatives are explicit inputs by design; this keeps the algebra layer
    separable from any particular finite-difference/parity implementation.
    """
    denom = np.asarray(components["sqrt_minus_g"], dtype=np.float64)
    shape = denom.shape
    if denom.ndim != 1 or np.any(denom <= 0.0) or not np.all(np.isfinite(denom)):
        raise ValueError("sqrt(-g) must be a finite positive one-dimensional profile")

    def arr(name: str, value: Any) -> np.ndarray:
        return _profile(name, value, shape)

    dt = {k: arr(f"d_t_metric[{k}]", d_t_metric[k]) for k in _METRIC_KEYS}
    dr = {k: arr(f"d_r_metric[{k}]", d_r_metric[k]) for k in _METRIC_KEYS}
    dt_den_t = arr("d_t_density_t", d_t_density_t)
    dt_den_r = arr("d_t_density_r", d_t_density_r)
    dr_flux_t = arr("d_r_flux_t", d_r_flux_t)
    dr_flux_r = arr("d_r_flux_r", d_r_flux_r)
    Ttt = arr("Ttt", components["Ttt"])
    Ttr = arr("Ttr", components["Ttr"])
    Trr = arr("Trr", components["Trr"])
    Tthth = arr("Tthth", components["Tthth"])

    source_t = -0.5 * (
        Ttt * dt["gtt"] + 2.0 * Ttr * dt["gtr"]
        + Trr * dt["grr"] + 2.0 * Tthth * dt["gthth"]
    )
    source_r = -0.5 * (
        Ttt * dr["gtt"] + 2.0 * Ttr * dr["gtr"]
        + Trr * dr["grr"] + 2.0 * Tthth * dr["gthth"]
    )
    div_t = (dt_den_t + dr_flux_t) / denom + source_t
    div_r = (dt_den_r + dr_flux_r) / denom + source_r
    if not np.all(np.isfinite(div_t)) or not np.all(np.isfinite(div_r)):
        raise FloatingPointError("covariant divergence contains non-finite values")
    return div_t, div_r
