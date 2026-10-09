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


def areal_radius(grid, geometry):
    r = np.asarray(grid.centers)
    return r * np.sqrt(geometry.b) / geometry.X


def invariant_chi(grid, geometry):
    r = np.asarray(grid.centers)
    R = areal_radius(grid, geometry)
    Rr = grid.cell_derivative_fourth(R, parity=1)
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
    """Interior residual of dM/dt = rhs on a sampled worldtube ledger.

    The first and last samples are not scored because their derivatives are
    one-sided. Interior samples adjacent to a sub-nanounit interval are also
    excluded: a terminal remainder that tiny can arise from floating-point
    time accumulation, and a centered derivative across it is ill-conditioned.
    This avoids reporting a numerical differencing artifact as a conservation
    residual. The source history itself is not modified.
    """
    t = np.asarray(times, dtype=float)
    m = np.asarray(masses, dtype=float)
    r = np.asarray(rhs, dtype=float)
    if t.shape != m.shape or t.shape != r.shape:
        raise ValueError("times, masses, and rhs must have identical shapes")
    if len(t) < 3:
        return np.asarray([])
    residual = np.full_like(m, np.nan, dtype=float)
    residual[1:-1] = np.gradient(m, t, edge_order=2)[1:-1] - r[1:-1]

    # The N=320 preregistration explicitly permits a <1e-9 final remainder
    # when the accumulated time is already within tolerance of final_time.
    # Exclude only the derivative samples touching such an interval.
    short_intervals = np.flatnonzero(np.diff(t) < 1.0e-9)
    for left in short_intervals:
        for index in (int(left), int(left) + 1):
            if 0 < index < len(t) - 1:
                residual[index] = np.nan
    return residual
