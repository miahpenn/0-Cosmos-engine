"""Reusable conservative matter RHS for synchronized time integrators."""
import numpy as np

from .matter_system import (
    BSSNMetricSlice,
    ConservedSpecies,
    MetricDerivativeSet,
    Species,
    _species_rhs_components,
)


def species_rhs(
    metric: BSSNMetricSlice,
    metric_derivatives: MetricDerivativeSet,
    state: ConservedSpecies,
    species: Species,
    dphi_t: np.ndarray | None = None,
    dphi_r: np.ndarray | None = None,
    beta_dm: float = -0.04,
) -> ConservedSpecies:
    """Return the semi-discrete conservative RHS without taking a finite step.

    Incoming radiation primitive recovery remains guarded. Actual finite-time
    and production-stage admissibility checks stay on their existing paths.
    """
    face_flux, source_e, source_s, q_e, q_s = _species_rhs_components(
        metric,
        metric_derivatives,
        state,
        species,
        dphi_t=dphi_t,
        dphi_r=dphi_r,
        beta_dm=beta_dm,
    )
    inv_dr = 1.0 / (metric.r[1] - metric.r[0])
    return ConservedSpecies(
        rest=-inv_dr * (face_flux[1:, 0] - face_flux[:-1, 0]),
        energy_t=(
            -inv_dr * (face_flux[1:, 1] - face_flux[:-1, 1])
            + (source_e + q_e)
        ),
        momentum_r=(
            -inv_dr * (face_flux[1:, 2] - face_flux[:-1, 2])
            + (source_s + q_s)
        ),
    )
