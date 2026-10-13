"""Semidiscrete conservative matter RHS, shared with finite-step advancement."""
from __future__ import annotations

import numpy as np

from . import matter_system
from .matter_system import (
    BSSNMetricSlice,
    ConservedSpecies,
    MetricDerivativeSet,
    Species,
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
    """Return the semidiscrete d(state)/dt without a finite-time update.

    Fluxes, hydrodynamic sources, and dark-matter coupling are evaluated by the
    same spatial operator used by evolve_species(). The actual update routine
    applies the flux divergence and sources in the legacy componentwise order;
    this function combines them as rates for synchronized integrators.
    """
    (
        face_flux,
        source_e,
        source_s,
        coupling_e,
        coupling_s,
        inv_dr,
    ) = matter_system._species_operator_components(
        metric,
        metric_derivatives,
        state,
        species,
        dphi_t=dphi_t,
        dphi_r=dphi_r,
        beta_dm=beta_dm,
    )

    rhs_rest = -inv_dr * (face_flux[1:, 0] - face_flux[:-1, 0])
    rhs_energy = -inv_dr * (face_flux[1:, 1] - face_flux[:-1, 1])
    rhs_energy += source_e + coupling_e
    rhs_momentum = -inv_dr * (face_flux[1:, 2] - face_flux[:-1, 2])
    rhs_momentum += source_s + coupling_s
    return ConservedSpecies(rhs_rest, rhs_energy, rhs_momentum)
