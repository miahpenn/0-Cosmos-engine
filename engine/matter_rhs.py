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
    recovery_metric: BSSNMetricSlice | None = None,
) -> ConservedSpecies:
    """Return d(state)/dt using the conservative transport operator.

    The underlying spatial operator is independent of the caller's dt, so a
    unit-step application is algebraically equivalent to extracting its RHS.
    Keeping this wrapper separate lets the production PIRK stage synchronize
    matter with the geometry without introducing a second physical evolution
    law.
    """
    # dt=1 is only a bookkeeping device for extracting the RHS. It is
    # not a physical update and must not trigger radiation admissibility
    # rejection. The actual RK stages are checked when they are formed.
    evolve_kwargs = dict(
        dphi_t=dphi_t,
        dphi_r=dphi_r,
        beta_dm=beta_dm,
        validate_physical_state=False,
    )
    # Keep the diagnostic-off invocation identical to the original call.
    if recovery_metric is not None:
        evolve_kwargs["recovery_metric"] = recovery_metric
    advanced = evolve_species(
        metric,
        metric_derivatives,
        state,
        species,
        1.0,
        **evolve_kwargs,
    )
    return ConservedSpecies(
        advanced.rest - state.rest,
        advanced.energy_t - state.energy_t,
        advanced.momentum_r - state.momentum_r,
    )
