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

# GEAR-67 on-site quartic of the S-D pair (locked input, recorded 2026-11-10).
#
# GEAR-67 models the pair as two identical elements with an on-site quartic lambda x^4 / 4. In pair
# coordinates (engine labels: S the surviving oscillator, D the ejected mode) the full potential is
#     V_SD = 1/2 S^2 - 1/2 D^2 + (lambda/8)(S^4 + 6 S^2 D^2 + D^4)      (g = 1, kappa = 0).
# lambda is an unpaid root (GEAR-67 3(e)): it has no derived value. The production machine runs the
# lambda = 0 case, in which S and D interact only through gravity; that was implicit until this
# constant made it explicit. Any lambda > 0 run is research and must say so.
SD_QUARTIC_LAMBDA = 0.0


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


def sd_potential(S, D, lam=None):
    """V_SD(S, D) of GEAR-67 in engine labels (see SD_QUARTIC_LAMBDA)."""
    lam = SD_QUARTIC_LAMBDA if lam is None else lam
    return 0.5 * S**2 - 0.5 * D**2 + (lam / 8.0) * (S**4 + 6.0 * S**2 * D**2 + D**4)


def sd_potential_gradient(S, D, lam=None):
    """(dV_SD/dS, dV_SD/dD). For lambda = 0 these are exactly (S, -D)."""
    lam = SD_QUARTIC_LAMBDA if lam is None else lam
    if lam == 0.0:
        return S, -D
    dS = S + (lam / 2.0) * S**3 + (1.5 * lam) * S * D**2
    dD = -D + (lam / 2.0) * D**3 + (1.5 * lam) * S**2 * D
    return dS, dD


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
            + invr * alphap * fp
            - geometry.alpha * Vp
        )
        return field_t, p_t

    VpS, VpD = sd_potential_gradient(fields.S, fields.D)
    St, PSt = one(fields.S, fields.PS, VpS)
    Dt, Pdt = one(fields.D, fields.PD, VpD)
    phit, Pit = one(
        fields.phi, fields.Pi, cosmos_potential_prime(fields.phi)
    )
    if rho_dm is not None:
        # Pi is a normal derivative, so coordinate-time evolution carries alpha.
        Pit = Pit + geometry.alpha * beta_dm * rho_dm
    return ScalarFields(St, PSt, Dt, Pdt, phit, Pit)


def scalar_projection(grid, geometry, fields):
    """Return the archive-normalized scalar Einstein projections."""
    invr = geometry.X**2 / geometry.a
    e = np.zeros_like(fields.S)
    pr = np.zeros_like(fields.S)
    pt = np.zeros_like(fields.S)
    j = np.zeros_like(fields.S)

    # The S-D potential is joint (GEAR-67): its quadratic parts are attributed to S and D as before
    # and the quartic, including the S-D cross term, is added once with D.
    V_quartic = sd_potential(fields.S, fields.D) - (0.5 * fields.S**2 - 0.5 * fields.D**2)
    if SD_QUARTIC_LAMBDA == 0.0:
        V_quartic = 0.0
    for f, p, V in (
        (fields.S, fields.PS, 0.5 * fields.S**2),
        (fields.D, fields.PD, -0.5 * fields.D**2 + V_quartic),
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
