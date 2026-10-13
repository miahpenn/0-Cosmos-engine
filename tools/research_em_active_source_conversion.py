"""Research-only Lorentz source conversion for active V5.5 storage conventions.

This is an algebraic adapter contract, not part of the production RHS. It uses
Gaussian units (c=1), signature (-,+,+,+), and the convention
Q_nu = F_{nu lambda} J^lambda for the matter force covector.

Input electric_covector_r is E_r = sqrt(gamma_rr) * Ehat_r, not the
orthonormal field itself. sqrt_gamma excludes the angular sin(theta) and
the 4*pi angular integration.

The active storage contracts differ:
* Eulerian variables (used by radiation): U_E=sqrt(gamma) E,
  U_Sr=sqrt(gamma) S_r. Add the returned Eulerian source pair directly.
* Mixed variables (used by dust/baryons/DM): U_t=sqrt(-g) T^t_t
  = sqrt(gamma) (-alpha E + beta^r S_r), U_r=sqrt(-g) T^t_r
  = sqrt(gamma) S_r. Convert sources with dU_t=-alpha*q_E+beta^r*q_Sr
  and dU_r=q_Sr, equivalently alpha*sqrt(gamma)*Q_t/r.

No charge carrier, current law, field normalization, or stage coupling is
selected here. The caller must supply same-stage current/geometry and apply
equal-and-opposite field/matter exchange at the actual integrator stage.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


def _finite_1d(values, name: str) -> np.ndarray:
    """Return a defensive copy of a finite, real, nonempty one-dimensional array."""
    try:
        raw = np.asarray(values)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite real 1-D array") from exc
    if raw.ndim != 1 or raw.size == 0 or raw.dtype.kind not in "iuf":
        raise ValueError(f"{name} must be a finite real 1-D array")
    result = np.array(raw, dtype=float, copy=True)
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must contain only finite values")
    return result


@dataclass(frozen=True)
class RadialLorentzSourceMap:
    """Conditional force source maps for both active conservative representations."""

    force_covector_t: np.ndarray
    force_covector_r: np.ndarray
    eulerian_energy: np.ndarray
    eulerian_momentum_r: np.ndarray
    mixed_energy_t: np.ndarray
    mixed_momentum_r: np.ndarray


def radial_lorentz_source_map(
    alpha,
    beta_r,
    sqrt_gamma,
    electric_covector_r,
    charge_density,
    spatial_current_r,
) -> RadialLorentzSourceMap:
    """Map a radial electric Lorentz force into active energy/momentum stores.

    Inputs are same-stage cell samples. The current is the Eulerian spatial
    current j_e^r (contravariant radial index); charge_density is the Eulerian
    charge density rho_e. For the declared convention,

        Q_r = rho_e E_r
        Q_t = -alpha E_r j_e^r + beta^r rho_e E_r
        q_E = alpha sqrt(gamma) E_r j_e^r
        q_Sr = alpha sqrt(gamma) rho_e E_r.

    The returned mixed sources equal sqrt(-g) Q_t/r and also the storage
    conversion -alpha*q_E + beta^r*q_Sr, q_Sr. The function does not mutate
    inputs, infer a current law, or inject sources into the production kernel.
    """
    arrays = {
        "alpha": _finite_1d(alpha, "alpha"),
        "beta_r": _finite_1d(beta_r, "beta_r"),
        "sqrt_gamma": _finite_1d(sqrt_gamma, "sqrt_gamma"),
        "electric_covector_r": _finite_1d(
            electric_covector_r, "electric_covector_r"
        ),
        "charge_density": _finite_1d(charge_density, "charge_density"),
        "spatial_current_r": _finite_1d(
            spatial_current_r, "spatial_current_r"
        ),
    }
    shapes = {value.shape for value in arrays.values()}
    if len(shapes) != 1:
        detail = ", ".join(f"{key}={value.shape}" for key, value in arrays.items())
        raise ValueError(f"all source arrays must have the same shape: {detail}")

    a = arrays["alpha"]
    beta = arrays["beta_r"]
    rootg = arrays["sqrt_gamma"]
    e_r = arrays["electric_covector_r"]
    rho_e = arrays["charge_density"]
    j_r = arrays["spatial_current_r"]
    if np.any(a <= 0.0):
        raise ValueError("alpha must be strictly positive")
    if np.any(rootg < 0.0):
        raise ValueError("sqrt_gamma must be nonnegative")

    with np.errstate(over="ignore", invalid="ignore"):
        q_r = rho_e * e_r
        q_t = -a * e_r * j_r + beta * rho_e * e_r
        q_energy = a * rootg * e_r * j_r
        q_momentum = a * rootg * q_r
        mixed_energy = a * rootg * q_t
        mixed_momentum = a * rootg * q_r

    outputs = {
        "force_covector_t": q_t,
        "force_covector_r": q_r,
        "eulerian_energy": q_energy,
        "eulerian_momentum_r": q_momentum,
        "mixed_energy_t": mixed_energy,
        "mixed_momentum_r": mixed_momentum,
    }
    if not all(np.all(np.isfinite(value)) for value in outputs.values()):
        raise ValueError("Lorentz source map produced a non-finite value")
    return RadialLorentzSourceMap(**outputs)
