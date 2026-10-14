"""Derived, research-only observables for saved one-dimensional GEAR slices.

This module post-processes captured states. It does not advance the state,
change the metric, or inject a source into the GEAR/COSMOS machine.

The archived corrected radial source audited on 2026-10-13 has SHA-256
0f2fc685f90b901275b33670c24aee3d2f1dbe8712d6e8de95d426ad2e1f5f4b.
The older proper-space helper has SHA-256
862bc9bb2ed0346450c25f8d75f4b2d435b10a3186bb8755dfc2aa3e5d6ba0c6.
These are provenance anchors, not hashes of this observer module.
"""
from __future__ import annotations

from math import isfinite
from typing import Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

from research.gear_proper_radial_v1 import archive_proper_radius

FloatArray = NDArray[np.float64]
CORRECTED_ARCHIVE_SOURCE_SHA256 = (
    "0f2fc685f90b901275b33670c24aee3d2f1dbe8712d6e8de95d426ad2e1f5f4b"
)
LEGACY_HELPER_SOURCE_SHA256 = (
    "862bc9bb2ed0346450c25f8d75f4b2d435b10a3186bb8755dfc2aa3e5d6ba0c6"
)


def _vec(values: ArrayLike, *, name: str, size: int | None = None) -> FloatArray:
    try:
        arr = np.asarray(values, dtype=np.float64)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric") from exc
    if arr.ndim != 1 or arr.size < 1:
        raise ValueError(f"{name} must be a nonempty one-dimensional array")
    if size is not None and arr.size != size:
        raise ValueError(f"{name} must have length {size}")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must be finite")
    return arr


def signed_mode_flux(
    radius: ArrayLike,
    lapse: ArrayLike,
    radial_scale_factor: ArrayLike,
    phi: ArrayLike,
    pi: ArrayLike,
    *,
    G: float = 1.0,
) -> FloatArray:
    """Evaluate archived signed mode flux on the supplied sample locations.

    F_X = -4*pi*G*r^2*N/a^3*Phi_X*Pi_X. Positive F follows the archive's
    outward-transport sign convention. Callers must state whether these are
    cell-centre diagnostic samples or true face samples; no interpolation is
    silently performed here.
    """
    r = _vec(radius, name="radius")
    n = _vec(lapse, name="lapse", size=r.size)
    a = _vec(radial_scale_factor, name="radial_scale_factor", size=r.size)
    p = _vec(phi, name="phi", size=r.size)
    q = _vec(pi, name="pi", size=r.size)
    if np.any(r <= 0.0) or np.any(a <= 0.0) or np.any(n <= 0.0):
        raise ValueError("radius, lapse, and radial_scale_factor must be positive")
    if not isfinite(G) or G < 0.0:
        raise ValueError("G must be finite and nonnegative")
    out = -4.0 * np.pi * G * r**2 * n / a**3 * p * q
    if not np.all(np.isfinite(out)):
        raise FloatingPointError("flux calculation produced nonfinite values")
    return out


def profile_peak_observables(
    radius: ArrayLike,
    radial_scale_factor: ArrayLike,
    lapse: ArrayLike,
    phi_s: ArrayLike,
    pi_s: ArrayLike,
    phi_d: ArrayLike,
    pi_d: ArrayLike,
    *,
    G: float = 1.0,
) -> dict[str, object]:
    """Return proper-radius and signed S/D flux profiles plus the D peak.

    The peak is the maximum signed D flux, not the maximum absolute flux.
    radius remains the areal-coordinate sample vector; ell is a metric-derived
    diagnostic map using the archived d ell = a dr rule.
    """
    r = _vec(radius, name="radius")
    a = _vec(radial_scale_factor, name="radial_scale_factor", size=r.size)
    n = _vec(lapse, name="lapse", size=r.size)
    ps = _vec(phi_s, name="phi_s", size=r.size)
    pis = _vec(pi_s, name="pi_s", size=r.size)
    pd = _vec(phi_d, name="phi_d", size=r.size)
    pid = _vec(pi_d, name="pi_d", size=r.size)
    ell = archive_proper_radius(r, a)
    fs = signed_mode_flux(r, n, a, ps, pis, G=G)
    fd = signed_mode_flux(r, n, a, pd, pid, G=G)
    idx = int(np.argmax(fd))
    return {
        "areal_radius": r.copy(),
        "proper_radius": ell,
        "flux_s": fs,
        "flux_d": fd,
        "d_peak_index": idx,
        "d_peak_areal_radius": float(r[idx]),
        "d_peak_proper_radius": float(ell[idx]),
        "d_peak_signed_flux": float(fd[idx]),
        "radius_label_note": "areal_radius and proper_radius are distinct; no speed inferred",
    }


def misner_sharp_face_mass(face_radius: float, a_left: float, a_right: float) -> float:
    """Estimate geometric Misner-Sharp mass at an interior face from metric data.

    Uses the arithmetic face value a_f=(a_left+a_right)/2 and
    M_MS=R_f*(1-a_f**(-2))/2. This is intentionally distinct from integrating
    a sampled energy density; callers must record the definition they use.
    """
    vals = (face_radius, a_left, a_right)
    if not all(isfinite(float(v)) for v in vals):
        raise ValueError("face-radius and metric values must be finite")
    if face_radius <= 0.0 or a_left <= 0.0 or a_right <= 0.0:
        raise ValueError("face radius and metric values must be positive")
    a_face = 0.5 * (a_left + a_right)
    return 0.5 * face_radius * (1.0 - 1.0 / (a_face * a_face))


def face_flux_from_centres(
    face_radius: float,
    spacing: float,
    lapse_left: float,
    lapse_right: float,
    a_left: float,
    a_right: float,
    phi_s_left: float,
    phi_s_right: float,
    pi_s_left: float,
    pi_s_right: float,
    phi_d_left: float,
    phi_d_right: float,
    pi_d_left: float,
    pi_d_right: float,
    *,
    G: float = 1.0,
) -> float:
    """Compute face flux from adjacent cell values using a direct face gradient.

    This declares the discrete convention rather than treating it as uniquely
    implied by the continuum formula: Phi_face=(X_right-X_left)/dr; N, a,
    and Pi are arithmetic face averages. It matches the independent audit
    convention for the corrected radial run.
    """
    vals = (
        face_radius, spacing, lapse_left, lapse_right, a_left, a_right,
        phi_s_left, phi_s_right, pi_s_left, pi_s_right,
        phi_d_left, phi_d_right, pi_d_left, pi_d_right, G,
    )
    if not all(isfinite(float(v)) for v in vals):
        raise ValueError("face-flux inputs must be finite")
    if face_radius <= 0.0 or spacing <= 0.0 or min(lapse_left, lapse_right, a_left, a_right) <= 0.0:
        raise ValueError("face radius, spacing, lapse, and metric must be positive")
    if G < 0.0:
        raise ValueError("G must be nonnegative")
    nf = 0.5 * (lapse_left + lapse_right)
    af = 0.5 * (a_left + a_right)
    phi_sf = (phi_s_right - phi_s_left) / spacing
    phi_df = (phi_d_right - phi_d_left) / spacing
    pi_sf = 0.5 * (pi_s_left + pi_s_right)
    pi_df = 0.5 * (pi_d_left + pi_d_right)
    return float(-4.0 * np.pi * G * face_radius**2 * nf / af**3 * (phi_sf*pi_sf + phi_df*pi_df))


def balance_residual(time: ArrayLike, mass_ms: ArrayLike, total_flux: ArrayLike) -> float:
    """Compute C=M_MS(t1)-M_MS(t0)+integral(F_total dt) by trapezoid.

    Use dense accepted-step/stage-start samples for time integration. Applying
    this only to sparse saved profile peaks can give a quadrature-dominated
    answer and is not equivalent to a per-step ledger.
    """
    t = _vec(time, name="time")
    m = _vec(mass_ms, name="mass_ms", size=t.size)
    f = _vec(total_flux, name="total_flux", size=t.size)
    if t.size < 2 or np.any(np.diff(t) <= 0.0):
        raise ValueError("time must contain at least two strictly increasing samples")
    return float(m[-1] - m[0] + np.sum(0.5 * (f[:-1] + f[1:]) * np.diff(t)))


def rk4_proper_time_increment(dt: float, lapse_stages: Sequence[float]) -> float:
    """Integrate d tau=N(0,t)dt consistently with one classical RK4 step."""
    if not isfinite(dt) or dt <= 0.0:
        raise ValueError("dt must be finite and positive")
    n = _vec(lapse_stages, name="lapse_stages", size=4)
    if np.any(n <= 0.0):
        raise ValueError("all lapse stages must be positive")
    return float(dt * (n[0] + 2*n[1] + 2*n[2] + n[3]) / 6.0)


__all__ = [
    "CORRECTED_ARCHIVE_SOURCE_SHA256",
    "LEGACY_HELPER_SOURCE_SHA256",
    "signed_mode_flux",
    "profile_peak_observables",
    "misner_sharp_face_mass",
    "face_flux_from_centres",
    "balance_residual",
    "rk4_proper_time_increment",
]
