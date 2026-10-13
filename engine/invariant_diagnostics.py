"""Invariant and numerical-witness definitions for V5.5.

These routines distinguish geometric invariants from coordinate diagnostics.
No event is converted into a bounce, branch change, or physical stop.
"""
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class InvariantSnapshot:
    hamiltonian_max: float
    momentum_max: float
    connection_max: float
    determinant_min: float
    lapse_min: float
    lapse_max: float
    trapping_min: float
    misner_sharp_mass: float
    misner_sharp_current: float


def finite_snapshot(s: InvariantSnapshot) -> bool:
    return all(
        np.isfinite(v)
        for v in (
            s.hamiltonian_max,
            s.momentum_max,
            s.connection_max,
            s.determinant_min,
            s.lapse_min,
            s.lapse_max,
            s.trapping_min,
            s.misner_sharp_mass,
            s.misner_sharp_current,
        )
    )


def trapping_indicator_from_areal_radius(
    spatial_gradient_R: np.ndarray,
    normal_time_derivative_R: np.ndarray,
    radial_inverse_metric: np.ndarray,
) -> np.ndarray:
    """Return g^(ab) partial_a R partial_b R."""
    return radial_inverse_metric * spatial_gradient_R**2 - normal_time_derivative_R**2


def misner_sharp_mass_from_chi(
    R: np.ndarray, chi: np.ndarray
) -> np.ndarray:
    return 0.5 * R * (1.0 - chi)


def current_balance(dm_dt: np.ndarray, flux: np.ndarray) -> np.ndarray:
    return np.asarray(dm_dt) - np.asarray(flux)


def scaled_max_norm(
    values: np.ndarray,
    scale: np.ndarray | float,
    floor: float = 1e-30,
) -> float:
    """Dimensionless max norm relative to an explicitly supplied scale."""
    v = np.max(np.abs(np.asarray(values)))
    s = np.max(np.abs(np.asarray(scale)))
    return float(v / max(float(s), floor))
