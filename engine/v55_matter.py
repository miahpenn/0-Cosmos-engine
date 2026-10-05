"""V5.5 matter bundle connecting scalar and fluid sectors.

The bundle keeps the physical roles distinct:
- local S and inverted D are archive scalar fields;
- COSMOS phi is the corrected normalized scalar;
- dark matter and baryons are conservative fluids;
- radiation is a p=rho/3 fluid;
- homogeneous Bianchi-I shear is represented only in the homogeneous
  cosmology lane and is absent from the exact spherical local geometry.

No D-to-matter identification or phenomenological interface source exists.
"""
from dataclasses import dataclass
import math
import numpy as np

from .matter_system import (
    BSSNMetricSlice,
    ConservedSpecies,
    Species,
    primitives,
    project_species,
    initialize_dust,
    initialize_radiation,
)
from .scalar_system import scalar_projection, ScalarFields
from .valencia import spherical_metric_from_bssn


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


def _metrics(g, s):
    return [
        spherical_metric_from_bssn(
            float(r), float(a), float(b), float(X),
            float(alpha), float(beta)
        )
        for r, a, b, X, alpha, beta in zip(
            g.centers, s.a, s.b, s.X, s.alpha, s.beta
        )
    ]


def initialize_from_archive(g, s, include_radiation=True) -> V55MatterState:
    """Initialize the corrected archive matter operating point."""
    metrics = _metrics(g, s)
    dm = initialize_dust(
        metrics, np.full(g.n, RHO_DM_PRESENT, dtype=float)
    )
    baryons = initialize_dust(
        metrics, np.full(g.n, RHO_B_PRESENT, dtype=float)
    )
    radiation = initialize_radiation(
        metrics,
        np.full(g.n, RHO_R_PRESENT if include_radiation else 0.0, dtype=float),
    )
    return V55MatterState(dm, baryons, radiation)


def normalized_fluid_densities(
    metric: BSSNMetricSlice,
    state: V55MatterState,
) -> tuple[np.ndarray, np.ndarray]:
    metrics = _metrics_from_slice(metric)
    dm = primitives(metrics, state.dark_matter, Species.DARK_MATTER)
    baryons = primitives(metrics, state.baryons, Species.BARYON)
    return (
        np.asarray([q.rho for q in dm]),
        np.asarray([q.rho for q in baryons]),
    )


def _metrics_from_slice(metric: BSSNMetricSlice):
    return [
        spherical_metric_from_bssn(
            float(r), float(a), float(b), float(X),
            float(alpha), float(beta)
        )
        for r, a, b, X, alpha, beta in zip(
            metric.r, metric.a, metric.b, metric.X,
            metric.alpha, metric.beta
        )
    ]


def total_fluid_projection(
    metric: BSSNMetricSlice,
    state: V55MatterState,
) -> dict[str, np.ndarray]:
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
    scalar_fields: ScalarFields,
    matter: V55MatterState,
) -> dict[str, np.ndarray]:
    """Combine scalar and conservative-fluid stress projections once each."""
    metric = metric_slice_from_q(g, s)
    scalar_e, scalar_pr, scalar_pt, scalar_j = scalar_projection(
        g, s, scalar_fields
    )
    fluid = total_fluid_projection(metric, matter)

    return {
        "rho": scalar_e + fluid["rho"],
        "pr": scalar_pr + fluid["pr"],
        "pt": scalar_pt + fluid["pt"],
        "j": scalar_j + fluid["j"],
        "fluid": fluid,
        "scalar": {
            "rho": scalar_e,
            "pr": scalar_pr,
            "pt": scalar_pt,
            "j": scalar_j,
        },
    }


def dm_density(metric: BSSNMetricSlice, matter: V55MatterState) -> np.ndarray:
    """Return archive-normalized DM rest density, not Einstein-normalized rho."""
    return normalized_fluid_densities(
        metric, matter
    )[0]


def exchange_pair(
    beta_dm: float,
    rho_dm: np.ndarray,
    dphi_t: np.ndarray,
    dphi_r: np.ndarray,
) -> dict[str, np.ndarray]:
    """Return the locked DM/scalar covector exchange pair."""
    q_t = beta_dm * rho_dm * dphi_t
    q_r = beta_dm * rho_dm * dphi_r
    return {
        "dm_t": q_t,
        "dm_r": q_r,
        "phi_t": -q_t,
        "phi_r": -q_r,
        "net_t": np.zeros_like(q_t),
        "net_r": np.zeros_like(q_r),
    }


def radiation_density_normalized(
    metric: BSSNMetricSlice, state: V55MatterState
) -> np.ndarray:
    projection = project_species(
        metric, state.radiation, Species.RADIATION
    )
    return projection["rho"] * (8.0 * math.pi)


def species_states(metric: BSSNMetricSlice, state: V55MatterState):
    metrics = _metrics_from_slice(metric)
    return {
        "dark_matter": primitives(
            metrics, state.dark_matter, Species.DARK_MATTER
        ),
        "baryons": primitives(
            metrics, state.baryons, Species.BARYON
        ),
        "radiation": primitives(
            metrics, state.radiation, Species.RADIATION
        ),
    }
