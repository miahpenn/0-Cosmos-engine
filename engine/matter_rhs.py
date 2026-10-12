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
    beta_rad: float = 0.0,
) -> ConservedSpecies:
    """Return the conservative RHS without advancing by an artificial dt.

    Radiation admissibility is checked on the incoming state while primitives
    are recovered, then on actual finite-time updates in evolve_species().
    A differential-rate query must not emulate a unit-time physical step.
    """
    rhs, _, _, _, _ = _species_rhs_components(
        metric,
        metric_derivatives,
        state,
        species,
        dphi_t=dphi_t,
        dphi_r=dphi_r,
        beta_dm=beta_dm,
        beta_rad=beta_rad,
    )
    return rhs
