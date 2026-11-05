#!/usr/bin/env python3
"""Equation-consistent matter-direction control on a frozen flat background.

This is an analytic instantaneous matter tangent, not a coupled Einstein-
matter solution. It uses continuum Gaussian derivatives and does not import
engine modules, a production RHS, initializer, kernel, or evolution step.
"""
from __future__ import annotations

import numpy as np

from tools.issue4_manufactured_regular_family import STATE_ARRAY_NAMES, PROFILE

# Locked values copied from the pinned scalar/Valencia source identities.
BETA_DM = -0.04
COSMOS_V0 = 8.242522415500654e-5
COSMOS_LAMBDA = 1.0
COSMOS_D_POT = 0.10


def _finite(name: str, value) -> np.ndarray:
    out = np.asarray(value, dtype=np.float64)
    if out.ndim != 1 or not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must be a finite one-dimensional array")
    return np.ascontiguousarray(out)


def gaussian_laplacian(r: np.ndarray, values: np.ndarray, width: float) -> np.ndarray:
    """Exact continuum 3-D radial Laplacian for A exp(-r^2/width^2)."""
    rr = _finite("radius", r)
    f = _finite("Gaussian profile", values)
    if width <= 0.0 or f.shape != rr.shape:
        raise ValueError("invalid Gaussian Laplacian inputs")
    return f * (4.0 * rr * rr / width**4 - 6.0 / width**2)


def gaussian_radial_derivative(r: np.ndarray, values: np.ndarray, width: float) -> np.ndarray:
    """Exact continuum radial derivative for A exp(-r^2/width^2)."""
    rr = _finite("radius", r)
    f = _finite("Gaussian profile", values)
    if width <= 0.0 or f.shape != rr.shape:
        raise ValueError("invalid Gaussian derivative inputs")
    return -2.0 * rr * f / width**2


def _fluid_conserved(r: np.ndarray, rho: np.ndarray, pressure: np.ndarray,
                     velocity: np.ndarray, *, radiation: bool) -> tuple[np.ndarray, ...]:
    """Canonical flat-areal conservative map: sqrt(gamma)=r^2."""
    v = _finite("fluid velocity", velocity)
    if np.any(np.abs(v) >= 1.0):
        raise ValueError("analytic fluid velocity must remain subluminal")
    W = 1.0 / np.sqrt(1.0 - v*v)
    h = rho + pressure
    energy = h * W*W - pressure
    momentum = h * W*W * v
    rest = np.zeros_like(rho) if radiation else r*r*rho*W
    return rest, r*r*energy, r*r*momentum


def analytic_matter_rates(r: np.ndarray, *, d1=None, d2=None) -> dict[str, np.ndarray]:
    """Return equation rates at the declared zero-velocity base.

    Scalar rates use the flat-slice pinned equations. Optional d1/d2 are the\n    pinned grid operators; absent operators select exact continuum Gaussians.
    Fluid rates satisfy spherical continuity and the Valencia/Euler momentum
    equation at rest, including the archived DM exchange covector.
    """
    rr = _finite("grid centers", r)
    if (d1 is None) != (d2 is None):
        raise ValueError("d1 and d2 must be supplied together")
    def radial_derivative(values, width):
        if d1 is None:
            return gaussian_radial_derivative(rr, values, width)
        return _finite("pinned first derivative", d1(values, parity=1))
    def radial_laplacian(values, width):
        if d1 is None:
            return gaussian_laplacian(rr, values, width)
        first = _finite("pinned first derivative", d1(values, parity=1))
        second = _finite("pinned second derivative", d2(values, parity=1))
        return second + 2.0 * first / rr
    scalar_shape = np.exp(-(rr / PROFILE["scalar_width"])**2)
    fluid_shape = np.exp(-(rr / PROFILE["fluid_width"])**2)

    scalar_fields = {
        name: amplitude * scalar_shape
        for name, amplitude in PROFILE["scalar_amplitudes"].items()
    }
    S, D, phi = scalar_fields["S"], scalar_fields["D"], scalar_fields["phi"]
    rho_dm = PROFILE["fluid_rest_density_amplitudes"]["dark_matter"] * fluid_shape

    S_accel = radial_laplacian(S, PROFILE["scalar_width"]) - S
    D_accel = radial_laplacian(D, PROFILE["scalar_width"]) + D
    phi_lap = radial_laplacian(phi, PROFILE["scalar_width"])
    phi_prime = -COSMOS_V0 * COSMOS_LAMBDA * np.exp(-COSMOS_LAMBDA * phi)
    phi_accel = phi_lap - phi_prime + BETA_DM * rho_dm

    phi_r = radial_derivative(phi, PROFILE["scalar_width"])
    rates = {
        "S_dot": np.zeros_like(rr), "PS_dot": S_accel,
        "D_dot": np.zeros_like(rr), "PD_dot": D_accel,
        "phi_dot": np.zeros_like(rr), "Pi_dot": phi_accel,
        "dark_matter_rho_dot": np.zeros_like(rr),
        "baryons_rho_dot": np.zeros_like(rr),
        "radiation_rho_dot": np.zeros_like(rr),
        "dark_matter_v_dot": BETA_DM * phi_r,
        "baryons_v_dot": np.zeros_like(rr),
    }
    rho_rad = PROFILE["fluid_rest_density_amplitudes"]["radiation"] * fluid_shape
    # At rest, (rho+p) dv/dt = -dp/dr for an uncoupled perfect fluid.
    dp_rad = radial_derivative(rho_rad, PROFILE["fluid_width"]) / 3.0
    rates["radiation_v_dot"] = -dp_rad / (rho_rad + rho_rad / 3.0)
    return rates


def construct_equation_consistent_family(centers: np.ndarray, s: float = 0.0, *, rates=None) -> dict[str, np.ndarray]:
    """Construct U(s)=U(0)+s*Udot_matter to first order at s=0.

    Geometry is deliberately held at exact flat areal values. Scalar fields
    and fluid densities have zero first derivative at the zero-momentum base;
    their conjugate momenta/velocities receive the equation-derived rates.
    """
    r = _finite("grid centers", centers)
    if r.size < 5 or np.any(r <= 0.0) or np.any(np.diff(r) <= 0.0):
        raise ValueError("centers must be strictly increasing positive cell centers")
    if not np.isfinite(s):
        raise ValueError("family parameter must be finite")

    rates = analytic_matter_rates(r)
    scalar_shape = np.exp(-(r / PROFILE["scalar_width"])**2)
    fluid_shape = np.exp(-(r / PROFILE["fluid_width"])**2)
    state: dict[str, np.ndarray] = {"grid.centers": r.copy()}
    for name, value in {
        "a": 1.0, "b": 1.0, "X": 1.0, "Aa": 0.0, "K": 0.0,
        "Lambda": 0.0, "alpha": 1.0, "beta": 0.0, "B": 0.0,
    }.items():
        state[f"geometry.{name}"] = np.full_like(r, value)

    for name, amplitude in PROFILE["scalar_amplitudes"].items():
        state[f"scalars.{name}"] = amplitude * scalar_shape
    for name, rate_name in (("PS", "PS_dot"), ("PD", "PD_dot"), ("Pi", "Pi_dot")):
        state[f"scalars.{name}"] = s * rates[rate_name]

    for species, amplitude in PROFILE["fluid_rest_density_amplitudes"].items():
        rho = amplitude * fluid_shape
        pressure = rho / 3.0 if species == "radiation" else np.zeros_like(rho)
        velocity = s * rates[f"{species}_v_dot"]
        rest, energy_t, momentum_r = _fluid_conserved(
            r, rho, pressure, velocity, radiation=(species == "radiation")
        )
        state[f"matter.{species}.rest"] = rest
        state[f"matter.{species}.energy_t"] = energy_t
        state[f"matter.{species}.momentum_r"] = momentum_r

    if set(state) != set(STATE_ARRAY_NAMES):
        raise RuntimeError("equation-consistent canonical state inventory mismatch")
    for key in STATE_ARRAY_NAMES:
        state[key] = _finite(key, state[key])
        if state[key].shape != r.shape:
            raise ValueError(f"{key} shape mismatch")
    return state


def construct_equation_consistent_triplet(centers: np.ndarray, epsilon: float, *, rates=None) -> dict[str, dict[str, np.ndarray]]:
    if not np.isfinite(epsilon) or epsilon <= 0.0:
        raise ValueError("epsilon must be positive and finite")
    return {
        "base": construct_equation_consistent_family(centers, 0.0, rates=rates),
        "plus": construct_equation_consistent_family(centers, float(epsilon), rates=rates),
        "minus": construct_equation_consistent_family(centers, -float(epsilon), rates=rates),
    }
