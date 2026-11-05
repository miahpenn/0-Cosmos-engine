#!/usr/bin/env python3
"""Initializer-free analytic regular-family constructor for Issue #4.

This module constructs arrays only. It does not import the engine, vendor,
production kernel, initializer, RHS, or any evolution step.
"""
from __future__ import annotations

import numpy as np

STATE_ARRAY_NAMES = (
    "grid.centers",
    *(f"geometry.{name}" for name in
      ("a", "b", "X", "Aa", "K", "Lambda", "alpha", "beta", "B")),
    *(f"scalars.{name}" for name in ("S", "PS", "D", "PD", "phi", "Pi")),
    *(f"matter.{species}.{name}" for species in
      ("dark_matter", "baryons", "radiation")
      for name in ("rest", "energy_t", "momentum_r")),
)

PROFILE = {
    "geometry": "exact_flat_areal_slice",
    "r_max": 40.0,
    "scalar_width": 8.0,
    "fluid_width": 10.0,
    "scalar_amplitudes": {"S": 1/10, "D": 1/20, "phi": 1/25},
    "fluid_rest_density_amplitudes": {
        "dark_matter": 1/1000, "baryons": 1/2000, "radiation": 1/3000,
    },
    "family_scalar_momentum_amplitudes": {
        "PS": 1/50, "PD": 1/75, "Pi": 1/100,
    },
    "family_velocity_slopes": {
        "dark_matter": 1/100, "baryons": 1/120, "radiation": 1/140,
    },
    "radiation_equation_of_state": "p=rho/3",
    "dust_equation_of_state": "p=0",
    "family_parameter": "s",
}


def _finite(name: str, values: np.ndarray) -> np.ndarray:
    result = np.asarray(values, dtype=np.float64)
    if result.ndim != 1 or not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must be a finite one-dimensional array")
    return np.ascontiguousarray(result)


def _fluid_conserved(rho: np.ndarray, pressure: np.ndarray,
                     velocity: np.ndarray) -> tuple[np.ndarray, ...]:
    """Flat-slice Valencia map: gamma_rr=alpha=sqrt(gamma)=1, beta=0."""
    v = _finite("fluid velocity", velocity)
    if np.any(np.abs(v) >= 1.0):
        raise ValueError("analytic fluid velocity must remain subluminal")
    W = 1.0 / np.sqrt(1.0 - v*v)
    h = rho + pressure
    energy = h * W*W - pressure
    momentum = h * W*W * v
    # D=rho W; T^t_t=-E; T^t_r=S_r on this exact flat ADM slice.
    return rho*W, -energy, momentum


def construct_regular_family(centers: np.ndarray, s: float = 0.0) -> dict[str, np.ndarray]:
    """Return U(s) as canonical state arrays, with explicit analytic profiles.

    All metric fields are exact Minkowski/flat-slice values. Scalar fields and
    fluid rest densities are even Gaussians. Radial velocities are odd regular
    profiles proportional to r*exp(-r^2/w^2). Scalar conjugate momenta are even.
    The parameter s changes only scalar conjugate momenta and fluid velocities.
    """
    r = _finite("grid centers", centers)
    if r.size < 5 or np.any(r <= 0.0) or np.any(np.diff(r) <= 0.0):
        raise ValueError("centers must be strictly increasing positive cell centers")
    if not np.isfinite(s):
        raise ValueError("family parameter s must be finite")

    scalar_shape = np.exp(-(r / PROFILE["scalar_width"])**2)
    fluid_shape = np.exp(-(r / PROFILE["fluid_width"])**2)
    state: dict[str, np.ndarray] = {"grid.centers": r.copy()}

    # Exact flat geometry; Lambda_metric and evolved Lambda are both zero.
    for name, value in {
        "a": 1.0, "b": 1.0, "X": 1.0, "Aa": 0.0, "K": 0.0,
        "Lambda": 0.0, "alpha": 1.0, "beta": 0.0, "B": 0.0,
    }.items():
        state[f"geometry.{name}"] = np.full_like(r, value)

    for name, amplitude in PROFILE["scalar_amplitudes"].items():
        state[f"scalars.{name}"] = amplitude * scalar_shape
    for name, amplitude in PROFILE["family_scalar_momentum_amplitudes"].items():
        state[f"scalars.{name}"] = s * amplitude * scalar_shape

    for species, amplitude in PROFILE["fluid_rest_density_amplitudes"].items():
        rho = amplitude * fluid_shape
        slope = PROFILE["family_velocity_slopes"][species]
        velocity = s * slope * r * fluid_shape
        pressure = rho / 3.0 if species == "radiation" else np.zeros_like(rho)
        rest, energy_t, momentum_r = _fluid_conserved(rho, pressure, velocity)
        state[f"matter.{species}.rest"] = rest
        state[f"matter.{species}.energy_t"] = energy_t
        state[f"matter.{species}.momentum_r"] = momentum_r

    missing = set(STATE_ARRAY_NAMES).difference(state)
    extra = set(state).difference(STATE_ARRAY_NAMES)
    if missing or extra:
        raise RuntimeError(f"canonical state inventory mismatch: missing={sorted(missing)}, extra={sorted(extra)}")
    for name in STATE_ARRAY_NAMES:
        state[name] = _finite(name, state[name])
        if state[name].shape != r.shape:
            raise ValueError(f"{name} has a shape inconsistent with grid.centers")
    return state


def construct_triplet(centers: np.ndarray, epsilon: float) -> dict[str, dict[str, np.ndarray]]:
    """Independently construct U(0), U(+epsilon), and U(-epsilon)."""
    if not np.isfinite(epsilon) or epsilon <= 0.0:
        raise ValueError("epsilon must be positive and finite")
    return {
        "base": construct_regular_family(centers, 0.0),
        "plus": construct_regular_family(centers, float(epsilon)),
        "minus": construct_regular_family(centers, -float(epsilon)),
    }
