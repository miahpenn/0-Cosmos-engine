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


def misner_sharp_shift_flux(R, Rr, Rdot, alpha, beta, rho, pr, j):
    """Shift contribution to dM/dt of the Misner-Sharp mass at FIXED COORDINATE r.

    With d_t = alpha n + beta d_r and the 3+1 split of T (rho, covariant j_r, mixed S^r_r = pr),
    d_t M|_r = 4 pi R^2 (T^r_t R' - T^r_r Rdot), where
        T^r_t = beta rho - alpha j^r - beta^2 j_r / alpha + beta S^r_r,
        T^r_r = S^r_r - beta j_r / alpha,
    and Rdot = d_t R at fixed r. The beta = 0 part is the existing ledger,
    flux_T = -4 pi R^2 alpha j^r R' and work_pR = -4 pi R^2 S^r_r Rdot; this returns the remainder,
        4 pi R^2 beta [ (rho + S^r_r - beta j_r / alpha) R' + j_r Rdot / alpha ],
    which is exactly zero for beta = 0 (production CMC). All inputs use the same normalization
    as flux_T and work_pR.
    """
    return 4.0 * np.pi * R**2 * beta * (
        (rho + pr - beta * j / alpha) * Rr + j * Rdot / alpha
    )


def current_residual(times, masses, rhs):
    """Residual of dM/dt = rhs on an already sampled worldtube ledger."""
    t = np.asarray(times, dtype=float)
    m = np.asarray(masses, dtype=float)
    r = np.asarray(rhs, dtype=float)
    if len(t) < 3:
        return np.asarray([])
    return np.gradient(m, t, edge_order=2) - r
