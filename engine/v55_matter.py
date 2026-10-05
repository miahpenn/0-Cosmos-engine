"""V5.5 matter bundle connecting scalar and fluid sectors.

The bundle keeps the physical roles distinct:
- local S and inverted D are archive scalar fields;
- COSMOS phi is the corrected normalized scalar;
- dark matter and baryons are conservative fluids;
- radiation is a p=rho/3 fluid;
- homogeneous Bianchi-I shear is represented only in the homogeneous
  cosmology lane and is explicitly absent from the exact spherical local
  geometry.

No D-to-matter identification or phenomenological interface source exists.
"""
from dataclasses import dataclass
import numpy as np

from . import reference_pirk_unified as q
from .matter_system import (
    BSSNMetricSlice,
    ConservedSpecies,
    Species,
    project_species,
    initialize_dust,
    initialize_radiation,
)


RHO_TOTAL_PRESENT = 9.64116e-5
RHO_DM_PRESENT = 2.5857e-5
RHO_B_PRESENT = 4.0306e-6
OMEGA_R_PRESENT = 9.2e-5
OMEGA_SHEAR_PRESENT = 1.745e-8
RHO_R_PRESENT = OMEGA_R_PRESENT * RHO_TOTAL_PRESENT
RHO_SHEAR_PRESENT = OMEGA_SHEAR_PRESENT * RHO_TOTAL_PRESENT


@dataclass
class V55MatterState:
    dark_matter: ConservedSpecies
    baryons: ConservedSpecies
    radiation: ConservedSpecies

    def copy(self) -> "V55MatterState":
        return V55MatterState(
            self.dark_matter.copy(),
            self.baryons.copy(),
            self.radiation.copy(),
        )


def metric_slice_from_q(g, s) -> BSSNMetricSlice:
    return BSSNMetricSlice(
        r=np.asarray(g.centers),
        a=np.asarray(s.a),
        b=np.asarray(s.b),
        X=np.asarray(s.X),
        alpha=np.asarray(s.alpha),
        beta=np.asarray(s.beta),
        Aa=np.asarray(s.Aa),
        K=np.asarray(s.K),
        Lambda=np.asarray(s.Lambda),
        B=np.asarray(s.B),
    )


def initialize_from_archive(g, s) -> V55MatterState:
    """Initialize the homogeneous COSMOS matter densities in the unified grid."""
    metrics = []
    for r, a, b, X, alpha, beta in zip(
        g.centers, s.a, s.b, s.X, s.alpha, s.beta
    ):
        from .valencia import spherical_metric_from_bssn
        metrics.append(
            spherical_metric_from_bssn(
                float(r), float(a), float(b), float(X),
                float(alpha), float(beta)
            )
        )

    dm = initialize_dust(
        metrics,
        np.full(g.n, RHO_DM_PRESENT, dtype=float),
    )
    baryons = initialize_dust(
        metrics,
        np.full(g.n, RHO_B_PRESENT, dtype=float),
    )
    radiation = initialize_radiation(
        metrics,
        np.full(g.n, RHO_R_PRESENT, dtype=float),
    )
    return V55MatterState(dm, baryons, radiation)


def total_fluid_projection(
    metric: BSSNMetricSlice,
    state: V55MatterState,
) -> dict[str, np.ndarray]:
    """Return rho, radial pressure, tangential pressure, and j for Einstein RHS."""
    pieces = [
        project_species(metric, state.dark_matter, Species.DARK_MATTER),
        project_species(metric, state.baryons, Species.BARYON),
        project_species(metric, state.radiation, Species.RADIATION),
    ]
    return {
        key: sum(piece[key] for piece in pieces)
        for key in ("rho", "pr", "pt", "j")
    }


def total_matter_projection(
    g,
    s,
    scalar_fields,
    matter: V55MatterState,
) -> dict[str, np.ndarray]:
    """Combine the frozen scalar source with the conservative fluid source."""
    metric = metric_slice_from_q(g, s)
    se, spr, spa, sj = q.matter_projection(g, s, scalar_fields)
    fluid = total_fluid_projection(metric, matter)
    return {
        "rho": se + fluid["rho"],
        "pr": spr + fluid["pr"],
        "pt": spa + fluid["pt"],
        "j": sj + fluid["j"],
        "fluid": fluid,
    }


def dm_density(metric: BSSNMetricSlice, matter: V55MatterState) -> np.ndarray:
    return project_species(
        metric, matter.dark_matter, Species.DARK_MATTER
    )["rho"] * (8.0 * np.pi)


def exchange_pair(beta_dm: float, rho_dm: np.ndarray,
                  dphi_t: np.ndarray, dphi_r: np.ndarray) -> dict[str, np.ndarray]:
    """Return the locked DM/scalar covector exchange pair."""
    q_t = beta_dm * rho_dm * dphi_t
    q_r = beta_dm * rho_dm * dphi_r
    return {
        "dm_t": q_t,
        "dm_r": q_r,
        "phi_t": -q_t,
        "phi_r": -q_r,
        "net_t": q_t - q_t,
        "net_r": q_r - q_r,
    }
