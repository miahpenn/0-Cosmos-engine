"""Reusable conservative matter RHS for synchronized time integrators."""
import numpy as np

from .matter_system import (
    BSSNMetricSlice,
    ConservedSpecies,
    MetricDerivativeSet,
    Species,
    evolve_species,
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
    """Return d(state)/dt using the conservative transport operator.

    The underlying spatial operator is independent of the caller's dt, so a
    unit-step application is algebraically equivalent to extracting its RHS.
    Keeping this wrapper separate lets the production PIRK stage synchronize
    matter with the geometry without introducing a second physical evolution
    law.
    """
    advanced = evolve_species(
        metric,
        metric_derivatives,
        state,
        species,
        1.0,
        dphi_t=dphi_t,
        dphi_r=dphi_r,
        beta_dm=beta_dm,
    )
    return ConservedSpecies(
        advanced.rest - state.rest,
        advanced.energy_t - state.energy_t,
        advanced.momentum_r - state.momentum_r,
    )
