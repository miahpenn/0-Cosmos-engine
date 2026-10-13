"""Research-only current closure for a separately specified charged-dust sector.

This is a conventional model option, not a charge assignment to existing GEAR
baryons/DM/radiation or real scalars S/D/phi, and not a claim that GEAR derives
electric charge. It is deliberately not imported by the production kernel.

For a pressureless dust species with constant caller-supplied charge-to-rest-mass
ratio eta=q/m, proper rest-mass density rho0 and Eulerian velocity v^r:
    W = (1 - gamma_rr (v^r)^2)^(-1/2)
    rho_e = eta * rho0 * W
    j_e^r = rho_e * v^r

The discrete charge flux is formed from the *same numerical rest-mass face flux*
used by that species:
    I_charge = 4*pi*eta*F_rest
where F_rest is the per-unit-solid-angle Valencia face flux for sqrt(gamma) D.
Using the existing numerical flux, instead of separately reconstructing the
continuum expression at faces, preserves exact discrete charge/rest-mass
continuity for constant eta. This module provides no matter evolution, Maxwell
field evolution, force injection, or production normalization.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

FOUR_PI = 4.0 * math.pi


def _finite_1d(values, name: str, *, minimum_size: int = 1) -> np.ndarray:
    try:
        raw = np.asarray(values)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite real 1-D array") from exc
    if raw.ndim != 1 or raw.size < minimum_size or raw.dtype.kind not in "iuf":
        raise ValueError(f"{name} must be a finite real 1-D array")
    result = np.array(raw, dtype=float, copy=True)
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must contain only finite values")
    return result


def _finite_scalar(value, name: str) -> float:
    if isinstance(value, (bool, np.bool_)) or not isinstance(
        value, (int, float, np.integer, np.floating)
    ):
        raise ValueError(f"{name} must be a finite real scalar")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{name} must be a finite real scalar")
    return result


@dataclass(frozen=True)
class ChargedDustCurrent:
    """Eulerian charge/current density derived from one charged dust primitive."""

    lorentz_factor: np.ndarray
    charge_density: np.ndarray
    spatial_current_r: np.ndarray


def charged_dust_eulerian_current(
    proper_rest_mass_density,
    velocity_r,
    gamma_rr,
    charge_to_rest_mass: float,
) -> ChargedDustCurrent:
    """Derive J's Eulerian charge/current from a separate charged-dust primitive.

    Arrays must be sampled on one common stage/grid. Density is proper rest-mass
    density, velocity_r is contravariant Eulerian v^r, and gamma_rr is the
    physical radial metric component. No speed clipping is applied: |v|>=1
    fails closed. The caller-supplied q/m may be positive, negative, or zero.
    """
    rho0 = _finite_1d(proper_rest_mass_density, "proper_rest_mass_density")
    count = rho0.size
    velocity = _finite_1d(velocity_r, "velocity_r", minimum_size=count)
    gamma = _finite_1d(gamma_rr, "gamma_rr", minimum_size=count)
    if velocity.shape != rho0.shape or gamma.shape != rho0.shape:
        raise ValueError("density, velocity and gamma_rr arrays must have identical shapes")
    eta = _finite_scalar(charge_to_rest_mass, "charge_to_rest_mass")
    if np.any(rho0 < 0.0):
        raise ValueError("proper rest-mass density must be nonnegative")
    if np.any(gamma <= 0.0):
        raise ValueError("gamma_rr must be strictly positive")

    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        v2 = gamma * velocity * velocity
        if np.any(v2 >= 1.0):
            raise ValueError("charged dust must be subluminal in the Eulerian frame")
        lorentz = 1.0 / np.sqrt(1.0 - v2)
        charge = eta * rho0 * lorentz
        current = charge * velocity
    if not all(np.all(np.isfinite(x)) for x in (lorentz, charge, current)):
        raise ValueError("charged-dust current closure produced a non-finite value")
    return ChargedDustCurrent(
        lorentz_factor=lorentz,
        charge_density=charge,
        spatial_current_r=current,
    )


def charged_dust_charge_flux_from_mass_flux(
    rest_mass_flux_faces,
    charge_to_rest_mass: float,
) -> np.ndarray:
    """Convert the actual Valencia rest-mass face flux into integrated charge flux.

    rest_mass_flux_faces must be the output of the same-stage numerical Valencia
    flux for the charged-dust species, per unit solid angle and already including
    the code's radial metric/lapse/shift factors. The returned outward flux is
    integrated over the full sphere. For constant q/m this construction makes
    discrete charge transport exactly proportional to that species' numerical
    rest-mass transport, including at the outer boundary.
    """
    mass_flux = _finite_1d(rest_mass_flux_faces, "rest_mass_flux_faces", minimum_size=2)
    eta = _finite_scalar(charge_to_rest_mass, "charge_to_rest_mass")
    with np.errstate(over="ignore", invalid="ignore"):
        charge_flux = FOUR_PI * eta * mass_flux
    if not np.all(np.isfinite(charge_flux)):
        raise ValueError("charged-dust flux conversion produced a non-finite value")
    return charge_flux
