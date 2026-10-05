"""Invariant diagnostics for the unified spherical production engine."""
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
    return all(np.isfinite(v) for v in (
        s.hamiltonian_max, s.momentum_max, s.connection_max,
        s.determinant_min, s.lapse_min, s.lapse_max,
        s.trapping_min, s.misner_sharp_mass, s.misner_sharp_current,
    ))


def trapping_indicator_from_areal_radius(
    spatial_gradient_R: np.ndarray,
    normal_time_derivative_R: np.ndarray,
    radial_inverse_metric: np.ndarray,
) -> np.ndarray:
    """Return g^(ab) d_a R d_b R without dividing by the lapse."""
    return radial_inverse_metric * spatial_gradient_R**2 - normal_time_derivative_R**2


def misner_sharp_mass_from_chi(R: np.ndarray, chi: np.ndarray) -> np.ndarray:
    return 0.5 * R * (1.0 - chi)


def current_balance(dm_dt: np.ndarray, flux: np.ndarray) -> np.ndarray:
    return dm_dt - flux


def relative_norm(x: np.ndarray, floor: float = 1e-30) -> float:
    return float(np.max(np.abs(x)) / max(float(np.max(np.abs(x))), floor))
