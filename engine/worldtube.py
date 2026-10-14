"""Worldtube observables derived from the solved spherical spacetime.

The worldtube is a diagnostic extraction at a chosen coordinate location; it
is not a shell, boundary condition, or separate dynamical system.
"""
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class WorldtubeObservation:
    coordinate_r: float
    areal_radius: float
    misner_sharp_mass: float
    chi: float
    lapse: float
    proper_clock_rate: float
    H_eff: float
    stress_flux: float
    work_term: float



def proper_radial_coordinate(grid, geometry):
    """Return proper radial distance at cell centres on the solved slice.

    The production BSSN inverse radial metric is gamma^rr = X**2 / a, so
    gamma_rr = a / X**2 and d ell / dr = sqrt(a) / abs(X). On the uniform cell-centred grid, midpoint
    quadrature places the first centre at half a cell and subsequent centres
    at cumulative full-cell increments. This is the main-engine counterpart
    of the archived GEAR proper-radius integral, expressed in the current
    metric variables rather than reusing the archive's polar-areal metric a.
    """
    centers = np.asarray(grid.centers, dtype=float)
    a = np.asarray(geometry.a, dtype=float)
    X = np.asarray(geometry.X, dtype=float)
    if centers.ndim != 1 or centers.size == 0:
        raise ValueError("grid centres must be a non-empty one-dimensional array")
    if a.shape != centers.shape or X.shape != centers.shape:
        raise ValueError("radial metric arrays must match the grid centres")
    dr = float(grid.dr)
    if not np.isfinite(dr) or dr <= 0.0:
        raise ValueError("grid spacing must be finite and positive")
    if (
        not np.all(np.isfinite(a))
        or not np.all(np.isfinite(X))
        or np.any(a <= 0.0)
        or np.any(X == 0.0)
    ):
        raise ValueError("radial metric must be finite and positive definite")
    radial_factor = np.sqrt(a / (X * X))
    if not np.all(np.isfinite(radial_factor)):
        raise ValueError("proper radial metric factor is non-finite")
    return dr * (np.cumsum(radial_factor) - 0.5 * radial_factor)


def areal_radius(grid, geometry):
    r = np.asarray(grid.centers)
    return r * np.sqrt(geometry.b) / geometry.X


def areal_radius_derivative(grid, geometry):
    """Return dR/dr with the regular-centre parity of the areal radius.

    R = r*sqrt(b)/X is odd under radial reflection when b and X are regular
    even fields. The derivative stencil must therefore use odd-field ghost
    parity for R; the resulting derivative is even and regular at the origin.
    """
    R = areal_radius(grid, geometry)
    return grid.cell_derivative_fourth(R, parity=-1)


def invariant_chi(grid, geometry):
    R = areal_radius(grid, geometry)
    Rr = areal_radius_derivative(grid, geometry)
    Ktheta = geometry.K / 3.0 - geometry.Aa / 2.0
    normal_dR = -R * Ktheta
    return (geometry.X**2 / geometry.a) * Rr**2 - normal_dR**2


def misner_sharp(grid, geometry):
    R = areal_radius(grid, geometry)
    chi = invariant_chi(grid, geometry)
    return 0.5 * R * (1.0 - chi)


def marginal_roots(grid, chi):
    r = np.asarray(grid.centers)
    roots = []
    for i in np.where(chi[:-1] * chi[1:] <= 0.0)[0]:
        if chi[i] != chi[i + 1]:
            roots.append(
                float(
                    r[i] - chi[i] * (r[i + 1] - r[i])
                    / (chi[i + 1] - chi[i])
                )
            )
    return roots


def summarize(grid, geometry, *, surface_r=10.0, H_eff=0.0,
              stress_flux=0.0, work_term=0.0):
    R = areal_radius(grid, geometry)
    chi = invariant_chi(grid, geometry)
    mass = 0.5 * R * (1.0 - chi)
    k = int(np.argmin(np.abs(grid.centers - surface_r)))
    return WorldtubeObservation(
        coordinate_r=float(grid.centers[k]),
        areal_radius=float(R[k]),
        misner_sharp_mass=float(mass[k]),
        chi=float(chi[k]),
        lapse=float(geometry.alpha[k]),
        proper_clock_rate=float(geometry.alpha[0]),
        H_eff=float(H_eff),
        stress_flux=float(stress_flux),
        work_term=float(work_term),
    )


def current_residual(times, masses, rhs):
    """Residual of dM/dt = rhs on an already sampled worldtube ledger."""
    t = np.asarray(times, dtype=float)
    m = np.asarray(masses, dtype=float)
    r = np.asarray(rhs, dtype=float)
    if len(t) < 3:
        return np.asarray([])
    return np.gradient(m, t, edge_order=2) - r
