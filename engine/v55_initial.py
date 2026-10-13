"""Corrected V5.5 production initial data.

The production initializer uses the archive-corrected COSMOS operating point:

    a=33.8983
    phi=0.179055
    phidot=0.00341913
    rho_DM=2.5857e-5
    rho_b =4.0306e-6

Radiation is enabled by default using rho_r = Omega_r*rho_total.
The archived Bianchi-I shear is deliberately excluded from the exact
spherical geometry.

The older reference_pirk_unified.make_initial() is retained unchanged as an
audit/reference artifact; it contains a historical variable-mapping mismatch
and is therefore not used here.
"""
from dataclasses import dataclass
import math
import numpy as np

from .scalar_system import cosmos_potential
from .v55_matter import (
    RHO_TOTAL_PRESENT,
    RHO_DM_PRESENT,
    RHO_B_PRESENT,
    RHO_R_PRESENT,
    V55MatterState,
    initialize_from_archive,
)


@dataclass(frozen=True)
class CorrectedInitialData:
    geometry: object
    scalars: object
    matter: V55MatterState
    H0: float
    include_radiation: bool


def _d1(grid, x, parity=1):
    if hasattr(grid, "cell_derivative_fourth"):
        return grid.cell_derivative_fourth(x, parity=parity)
    return grid.d1(x, parity)


# Origin half-cell weights integrate a smooth even source through r^4 exactly.
_ORIGIN_HALF_CELL_WEIGHTS = (
    401.0 / 720.0,
    -31.0 / 480.0,
    11.0 / 1440.0,
)


def _integrate_initial_origin_half_cell(src, dr):
    """Integrate the origin half-cell using even-series quadrature through r^4."""
    values = np.asarray(src, dtype=float)
    step = float(dr)
    if values.ndim != 1 or values.size < 3:
        raise ValueError("origin half-cell quadrature requires at least three samples")
    if not np.isfinite(step) or step <= 0.0:
        raise ValueError("initial radial quadrature requires finite positive dr")
    if not np.all(np.isfinite(values[:3])):
        raise FloatingPointError("origin half-cell source samples are non-finite")
    result = step * float(np.dot(_ORIGIN_HALF_CELL_WEIGHTS, values[:3]))
    if not np.isfinite(result):
        raise FloatingPointError("origin half-cell quadrature produced non-finite integral")
    return result


def _integrate_initial_radial_source(src, dr):
    """Integrate an even radial source from the origin to every cell centre.

    The origin half-cell uses the analytic even-series quadrature. Remaining
    intervals use cubic Lagrange-integral weights; the first interval uses an
    even ghost and the final interval a one-sided cubic rule. The source law
    is unchanged: only its radial integration is made centre-compatible.
    """
    values = np.asarray(src, dtype=float)
    step = float(dr)
    if values.ndim != 1 or values.size < 4:
        raise ValueError("initial radial quadrature requires at least four cells")
    if not np.isfinite(step) or step <= 0.0:
        raise ValueError("initial radial quadrature requires finite positive dr")
    if not np.all(np.isfinite(values)):
        raise FloatingPointError("initial radial quadrature source is non-finite")

    half_cell = _integrate_initial_origin_half_cell(values, step)
    increments = np.empty(values.size - 1, dtype=float)
    # Interval [r_0, r_1]; the even ghost has src[-1] == src[0].
    increments[0] = step * (
        12.0 * values[0] + 13.0 * values[1] - values[2]
    ) / 24.0
    # Interior intervals [r_(i-1), r_i], integrated cubic interpolation.
    for i in range(2, values.size - 1):
        increments[i - 1] = step * (
            -values[i - 2] + 13.0 * values[i - 1]
            + 13.0 * values[i] - values[i + 1]
        ) / 24.0
    # Final interval: one-sided cubic interpolation at the outer edge.
    increments[-1] = step * (
        values[-4] - 5.0 * values[-3]
        + 19.0 * values[-2] + 9.0 * values[-1]
    ) / 24.0

    integral = np.cumsum(np.concatenate(([half_cell], increments)))
    if not np.all(np.isfinite(integral)):
        raise FloatingPointError("initial radial quadrature produced non-finite integral")
    return integral


def build_initial_data(
    grid,
    vacuum_state_factory,
    amplitude=0.01,
    width=7.0,
    D_amplitude=1.0e-10,
    include_radiation=True,
):
    r = np.asarray(grid.centers)
    n = len(r)
    dr = float(grid.dr)
    rho_r = RHO_R_PRESENT if include_radiation else 0.0

    phi0 = 0.179055
    pi0 = 0.00341913
    dm0 = RHO_DM_PRESENT
    b0 = RHO_B_PRESENT

    # Archive-normalized spherical budget. The exact spherical branch
    # excludes anisotropic shear, so it is not silently folded into rho.
    rho0 = 0.5 * pi0**2 + float(cosmos_potential(np.array([phi0]))[0])
    rho0 += dm0 + b0 + rho_r
    H0 = math.sqrt(max(rho0 / 3.0, 0.0))

    S = amplitude * np.exp(-(r / width) ** 2)
    D = D_amplitude * np.exp(-(r / width) ** 2)
    PS = np.zeros_like(r)
    PD = D.copy()
    phi = np.full(n, phi0)
    Pi = np.full(n, pi0)
    Kt = np.full(n, -H0)

    Sp = _d1(grid, S, 1)
    Dp = _d1(grid, D, 1)
    Ktp = _d1(grid, Kt, 1)
    j = -(PD * Dp)
    B = np.ones_like(r)

    # Same constrained areal initial-data reconstruction used by the archived
    # V5.5 gate, with corrected COSMOS variables and optional radiation.
    for _ in range(6):
        rhoS = 0.5 * (PS**2 + B * Sp**2) + 0.5 * S**2
        rhoD = 0.5 * (PD**2 + B * Dp**2) - 0.5 * D**2
        rhoPhi = (0.5 * Pi**2 + cosmos_potential(phi)) / (8.0 * math.pi)
        rhoMatter = (dm0 + b0 + rho_r) / (8.0 * math.pi)
        rho = rhoS + rhoD + rhoPhi + rhoMatter
        Kr = Kt + r * Ktp + 4.0 * math.pi * r * j
        src = (
            1.0 - 8.0 * math.pi * r * r * rho
            + 2.0 * r * r * Kr * Kt
            + r * r * Kt * Kt
        )
        integ = _integrate_initial_radial_source(src, dr)
        B = integ / r

    if np.any(~np.isfinite(B)) or np.any(B <= 0.0):
        raise FloatingPointError("invalid constrained initial metric")

    A = 1.0 / np.sqrt(B)
    Kr = Kt + r * Ktp + 4.0 * math.pi * r * j
    K = Kr + 2.0 * Kt
    Aa = Kr - K / 3.0
    a = A ** (4.0 / 3.0)
    b = A ** (-2.0 / 3.0)
    X = A ** (-1.0 / 3.0)
    Lambda = (
        _d1(grid, a, 1) / (2.0 * a * a)
        - _d1(grid, b, 1) / (a * b)
        + 2.0 / r * (1.0 / b - 1.0 / a)
    )

    state = vacuum_state_factory(grid)
    state.a[:] = a
    state.b[:] = b
    state.X[:] = X
    state.alpha[:] = 1.0
    state.beta[:] = 0.0
    state.Aa[:] = Aa
    state.K[:] = K
    state.Lambda[:] = Lambda
    state.B[:] = 0.75 * Lambda

    scalars = __import__("engine.scalar_system", fromlist=["ScalarFields"]).ScalarFields(
        S, PS, D, PD, phi, Pi
    )
    matter = initialize_from_archive(
        grid,
        state,
        include_radiation=include_radiation,
    )
    return CorrectedInitialData(
        geometry=state,
        scalars=scalars,
        matter=matter,
        H0=H0,
        include_radiation=include_radiation,
    )
