"""Derived COSMOS/cycle observables from the unified spacetime."""
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class UnifiedCosmosObservables:
    H_eff: float
    e_folds: float
    phi_outer: float
    Pi_outer: float
    rho_outer: float
    rho_dm_outer: float
    rho_b_outer: float
    rho_r_outer: float


def effective_hubble(geometry, volumes) -> float:
    return -float(
        np.sum(volumes * geometry.K) / np.sum(volumes)
    ) / 3.0


def append_efolds(old_efolds: float, H_eff: float, dt: float) -> float:
    """Integrate d ln a_eff / dt = H_eff without prescribing a scale factor."""
    return old_efolds + H_eff * dt


def outer_observables(grid, geometry, scalars, total_stress,
                      dm_density, baryon_density, radiation_density,
                      e_folds: float) -> UnifiedCosmosObservables:
    outer = grid.centers >= 0.8 * grid.r_max
    return UnifiedCosmosObservables(
        H_eff=effective_hubble(geometry, grid.volumes),
        e_folds=float(e_folds),
        phi_outer=float(np.mean(scalars.phi[outer])),
        Pi_outer=float(np.mean(scalars.Pi[outer])),
        rho_outer=float(np.mean(total_stress.rho[outer])),
        rho_dm_outer=float(np.mean(dm_density[outer])),
        rho_b_outer=float(np.mean(baryon_density[outer])),
        rho_r_outer=float(np.mean(radiation_density[outer])),
    )
