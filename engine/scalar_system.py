"""V5.5 scalar-field evolution and stress projections.

The scalar equations are the frozen spherical BSSN equations from the
archive. The helper accepts any grid exposing fourth-order cell derivatives
with the parity convention used by the pinned reference kernel.
"""
from dataclasses import dataclass
import numpy as np
import math

KAPPA = 8.0 * math.pi
V0 = 8.242522415500654e-5
D_POT = 0.10
LAMBDA = 1.0


@dataclass(frozen=True)
class ScalarFields:
    S: np.ndarray
    PS: np.ndarray
    D: np.ndarray
    PD: np.ndarray
    phi: np.ndarray
    Pi: np.ndarray

    def copy(self) -> "ScalarFields":
        return ScalarFields(
            self.S.copy(), self.PS.copy(), self.D.copy(),
            self.PD.copy(), self.phi.copy(), self.Pi.copy()
        )


def _d1(grid, x, parity):
    if hasattr(grid, "cell_derivative_fourth"):
        return grid.cell_derivative_fourth(x, parity=parity)
    return grid.d1(x, parity)


def _d2(grid, x, parity):
    if hasattr(grid, "cell_second_derivative_fourth"):
        return grid.cell_second_derivative_fourth(x, parity=parity)
    return grid.d2(x, parity)


def cosmos_potential(phi):
    return V0 * (np.exp(-LAMBDA * phi) - D_POT)


def cosmos_potential_prime(phi):
    return -V0 * LAMBDA * np.exp(-LAMBDA * phi)


def scalar_rhs_arrays(grid, geometry, fields, beta_dm=-0.04,
                      rho_dm=None):
    """Return d_t(S,PS,D,PD,phi,Pi) for the frozen V5.5 scalar system."""
    r = np.asarray(grid.centers)
    invr = geometry.X**2 / geometry.a
    ap = _d1(grid, geometry.a, 1)
    bp = _d1(grid, geometry.b, 1)
    Xp = _d1(grid, geometry.X, 1)
    alphap = _d1(grid, geometry.alpha, 1)

    def one(f, p, Vp):
        fp = _d1(grid, f, 1)
        fpp = _d2(grid, f, 1)
        lap = invr * (
            fpp - fp * (
                0.5 * ap / geometry.a
                - bp / geometry.b
                + Xp / geometry.X
                - 2.0 / r
            )
        )
        field_t = geometry.alpha * p + geometry.beta * fp
        p_t = (
            geometry.beta * _d1(grid, p, 1)
            + geometry.alpha * lap
            + geometry.alpha * geometry.K * p
            + geometry.X**-2 * alphap * fp
            - geometry.alpha * Vp
        )
        return field_t, p_t

    St, PSt = one(fields.S, fields.PS, fields.S)
    Dt, Pdt = one(fields.D, fields.PD, -fields.D)
    phit, Pit = one(
        fields.phi, fields.Pi, cosmos_potential_prime(fields.phi)
    )
    if rho_dm is not None:
        # The covariant scalar-DM exchange contributes alpha * beta * rho
        # to the coordinate-time Pi RHS. This is the opposite scalar source
        # to the conservative DM source in the same ADM slice.
        Pit = Pit + geometry.alpha * beta_dm * rho_dm
    return ScalarFields(St, PSt, Dt, Pdt, phit, Pit)


def scalar_projection(grid, geometry, fields):
    """Return the archive-normalized scalar Einstein projections."""
    invr = geometry.X**2 / geometry.a
    e = np.zeros_like(fields.S)
    pr = np.zeros_like(fields.S)
    pt = np.zeros_like(fields.S)
    j = np.zeros_like(fields.S)

    for f, p, V in (
        (fields.S, fields.PS, 0.5 * fields.S**2),
        (fields.D, fields.PD, -0.5 * fields.D**2),
        (fields.phi, fields.Pi, cosmos_potential(fields.phi)),
    ):
        fp = _d1(grid, f, 1)
        gr2 = invr * fp * fp
        scale = KAPPA if f is fields.phi else 1.0
        e += (0.5 * p * p + 0.5 * gr2 + V) / scale
        pr += (0.5 * p * p + 0.5 * gr2 - V) / scale
        pt += (0.5 * p * p - 0.5 * gr2 - V) / scale
        j += (-p * fp) / scale

    return e, pr, pt, j
