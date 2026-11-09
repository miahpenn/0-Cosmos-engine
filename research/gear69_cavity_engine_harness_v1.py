"""RESEARCH ONLY. Gear 69/70B nonlinear cavity in the production engine (CMC, beta = 0).

Measures the growth of the uncoupled mode D (V_D = -D^2/2, mu = 1, G = 1 normalization) in a
cavity with a Dirichlet wall at FIXED AREAL RADIUS R = 10, through V55ProductionKernel.step, and
compares with the polar-areal prediction (GEAR-86/87): along the centre worldline,

    Gamma_tau = (d D/d tau_c) / D_c = PD(0) / D(0)      (beta = 0: the centre normal is the worldline)
    (Gamma_tau / gamma_lin - 1) / D_c^2  ->  2 w(0) - alpha2(0) = 4.5826.

What the harness changes, by monkeypatching only (production files untouched):
  - V0 = 0, phi = Pi = 0, S = PS = 0, all matter densities 1e-30 (the fluid solvers need a positive rest current): a pure D cavity. (The COSMOS potential's
    V0 (1 - D_POT) would otherwise add a uniform energy density, GEAR-85.)
  - The wall: the coordinate r_w where the areal radius r sqrt(b)/X equals R is found on every RHS
    call. For D's derivatives, cells beyond r_w hold odd-reflected ghost values about r_w (cubic
    interpolation), so D = 0 at the areal wall to sub-cell accuracy. D's stress-energy and RHS are
    zero beyond r_w, so the exterior is vacuum and the ghosts never source the geometry.
Initial data: the gear's growing mode to O(A^3), D = A u + A^3 w, d_t D = gamma (A u + 3 A^3 w),
PD = d_t D / alpha
(u = j0(k r), k = pi/R), with the areal reconstruction of engine/v55_initial (K_theta = 0,
K_r from the momentum constraint). Not time-symmetric: the observable is a property of the growing-
mode solution along the centre worldline, so it does not need the slices to coincide.
"""
from __future__ import annotations

import dataclasses
import math

import numpy as np
from scipy.integrate import cumulative_simpson, solve_bvp
from scipy.interpolate import CubicSpline
from scipy.special import spherical_jn

import engine.production_kernel as pk
import engine.scalar_system as ss
import engine.v55_matter as vm
from engine import v55_initial
from engine.matter_system import _metric_arrays, initialize_dust, initialize_radiation
from engine.v55_matter import V55MatterState, metric_slice_from_q

R_WALL, MU, G = 10.0, 1.0, 1.0
K_WAVE = math.pi / R_WALL
GAMMA_LIN = math.sqrt(MU * MU - K_WAVE * K_WAVE)
PREDICTED_SLOPE = 4.5826   # 2 w(0) - alpha2(0), GEAR-86 (polar-areal, G = 1, R = 10, mu = 1)
MATTER_FLOOR = 1.0e-30     # dust/radiation density; the fluid solvers need a positive rest current

_ORIG_RHS = pk.scalar_rhs_arrays
_ORIG_PROJ = vm.scalar_projection


def wall_coordinate(grid, geometry) -> float:
    r = np.asarray(grid.centers, float)
    areal = r * np.sqrt(np.asarray(geometry.b)) / np.asarray(geometry.X)
    i = int(np.searchsorted(areal, R_WALL))
    if i == 0 or i >= len(r):
        raise RuntimeError("areal wall outside the grid")
    f = (R_WALL - areal[i - 1]) / (areal[i] - areal[i - 1])
    return float(r[i - 1] + f * (r[i] - r[i - 1]))


def _ghost(r, values, r_w):
    v = np.asarray(values, float).copy()
    inside = r < r_w
    sp = CubicSpline(np.r_[r[inside], r_w], np.r_[v[inside], 0.0])
    out = r > r_w
    mirror = 2.0 * r_w - r[out]
    ghost = np.where(mirror > 0.0, -sp(np.clip(mirror, r[0], r_w)), 0.0)
    ghost[(r[out] - r_w) > 4.0 * (r[1] - r[0])] = 0.0
    v[out] = ghost
    return v


def _wall_fields(grid, geometry, fields):
    r = np.asarray(grid.centers, float)
    r_w = wall_coordinate(grid, geometry)
    return r, r_w, dataclasses.replace(
        fields, D=_ghost(r, fields.D, r_w), PD=_ghost(r, fields.PD, r_w)
    )


def _rhs_with_wall(grid, geometry, fields, *args, **kwargs):
    r, r_w, fg = _wall_fields(grid, geometry, fields)
    out = _ORIG_RHS(grid, geometry, fg, *args, **kwargs)
    beyond = r > r_w
    Dt = np.where(beyond, 0.0, out.D)
    PDt = np.where(beyond, 0.0, out.PD)
    return dataclasses.replace(out, D=Dt, PD=PDt)


def _projection_with_wall(grid, geometry, fields):
    r, r_w, fg = _wall_fields(grid, geometry, fields)
    e, pr, pt, j = _ORIG_PROJ(grid, geometry, fg)
    beyond = r > r_w
    return tuple(np.where(beyond, 0.0, x) for x in (e, pr, pt, j))


def install():
    """Apply the harness patches (V0 = 0, wall in the D RHS and in the stress projection)."""
    ss.V0 = 0.0
    pk.scalar_rhs_arrays = _rhs_with_wall
    vm.scalar_projection = _projection_with_wall


def w_profile(r):
    """GEAR-86/87 O(A^3) mode-shape correction w(r) and alpha2 (polar-areal), from the gear script."""
    g2 = GAMMA_LIN ** 2
    k = K_WAVE
    rr = np.linspace(0, R_WALL, 20001)
    u = spherical_jn(0, k * rr); up = -k * spherical_jn(1, k * rr)
    upp = np.empty_like(rr); upp[1:] = -2 * up[1:] / rr[1:] - k * k * u[1:]; upp[0] = -k * k / 3
    rho2 = 0.5 * (g2 * u * u + up * up - MU * MU * u * u)
    m2 = cumulative_simpson(4 * math.pi * G * rr * rr * rho2, x=rr, initial=0.0)
    a2 = np.zeros_like(rr); a2[1:] = m2[1:] / rr[1:]
    ap = np.zeros_like(rr)
    ap[1:] = a2[1:] / rr[1:] + 2 * math.pi * G * rr[1:] * (g2 * u[1:] ** 2 + up[1:] ** 2 + MU * MU * u[1:] ** 2)
    ia = cumulative_simpson(ap, x=rr, initial=0.0); al2 = -a2[-1] - (ia[-1] - ia)
    q, s = al2 - a2, al2 + a2
    qp = CubicSpline(rr, q).derivative()(rr)
    div = np.empty_like(rr); div[1:] = 2 * q[1:] * up[1:] / rr[1:] + qp[1:] * up[1:] + q[1:] * upp[1:]; div[0] = -k * k * q[0]
    S3 = -div - MU * MU * s * u - 3 * g2 * q * u
    src = CubicSpline(rr, S3); lam = MU * MU - 9 * g2

    def f(x, y):
        c = src(x) - lam * y[0]; xs = np.where(x == 0, 1.0, x)
        return np.vstack((y[1], np.where(x == 0, c / 3, c - 2 * y[1] / xs)))
    mesh = np.linspace(0, R_WALL, 1000)
    sol = solve_bvp(f, lambda ya, yb: np.array([ya[1], yb[0]]), mesh, np.zeros((2, 1000)), tol=1e-8, max_nodes=20000)
    assert sol.success
    rin = np.clip(r, 0, R_WALL)
    return np.where(r < R_WALL, sol.sol(rin)[0], 0.0), float(2 * np.atleast_1d(sol.sol(0.0)[0])[0] - al2[0])


def initial_state(N: int, A: float, r_max: float = 20.0, use_w: bool = True):
    kernel = pk.V55ProductionKernel()
    st = kernel.initialize(resolution=N, r_max=r_max, include_radiation=False)
    grid = st.grid
    r = np.asarray(grid.centers, float); dr = float(grid.dr)
    inside = r < R_WALL
    u = np.where(inside, spherical_jn(0, K_WAVE * r), 0.0)
    w = w_profile(r)[0] if use_w else 0.0 * r
    D = A * u + A ** 3 * w
    PD = GAMMA_LIN * (A * u + 3.0 * A ** 3 * w)
    z = np.zeros_like(r)
    # Areal reconstruction as in engine/v55_initial (K_theta = 0); D sector normalized with G = 1.
    d1 = lambda x: v55_initial._d1(grid, x, 1)  # noqa: E731
    Dp = d1(_ghost(r, D, R_WALL))
    Bm = np.ones_like(r)
    for _ in range(8):
        rho = np.where(inside, 0.5 * (PD ** 2 + Bm * Dp ** 2) - 0.5 * D ** 2, 0.0)
        src = 1.0 - 8.0 * math.pi * r * r * rho
        Bm = v55_initial._cumulative_origin_integral(grid, r, dr, src) / r
    j = np.where(inside, -(PD * Dp), 0.0)
    Kr = 4.0 * math.pi * r * j
    Am = 1.0 / np.sqrt(Bm)
    geo = st.geometry
    geo.a[:] = Am ** (4.0 / 3.0); geo.b[:] = Am ** (-2.0 / 3.0); geo.X[:] = Am ** (-1.0 / 3.0)
    geo.K[:] = Kr; geo.Aa[:] = Kr - Kr / 3.0
    geo.Lambda[:] = (d1(geo.a) / (2.0 * geo.a * geo.a) - d1(geo.b) / (geo.a * geo.b)
                     + 2.0 / r * (1.0 / geo.b - 1.0 / geo.a))
    geo.beta[:] = 0.0; geo.B[:] = 0.0
    st.scalars = dataclasses.replace(st.scalars, S=z.copy(), PS=z.copy(), D=D, PD=PD, phi=z.copy(), Pi=z.copy())
    metrics = _metric_arrays(metric_slice_from_q(grid, geo))
    # The fluid solvers require a positive rest current; 1e-30 is ~25 orders below the D energy.
    tiny = np.full_like(r, MATTER_FLOOR)
    st.matter = V55MatterState(initialize_dust(metrics, tiny.copy()), initialize_dust(metrics, tiny.copy()),
                               initialize_radiation(metrics, tiny.copy()), None)
    geo.alpha[:] = kernel._solve_lapse_for_state(grid, geo, st.scalars, st.matter, None)[0]
    # PD is the NORMAL derivative, PD = d_t D / alpha (beta = 0). The growing-mode rate above is per
    # coordinate time, so divide by the lapse; without this the data carry an O(A^2) error (alpha - 1
    # ~ -149 A^2) that excites the decaying mode. One pass suffices at O(A^2): the lapse's dependence
    # on PD is itself O(A^2) smaller.
    st.scalars = dataclasses.replace(st.scalars, PD=PD / np.asarray(geo.alpha))
    geo.alpha[:] = kernel._solve_lapse_for_state(grid, geo, st.scalars, st.matter, None)[0]
    return kernel, st


def centre_value(grid, values):
    r = np.asarray(grid.centers, float)[:3]
    A = np.vstack([np.ones(3), r ** 2, r ** 4]).T
    return float(np.linalg.solve(A, np.asarray(values, float)[:3])[0])


def run(N: int, A: float, T: float, cfl: float = 0.06, every: float = 0.25, r_max: float = 20.0, use_w: bool = True):
    """Evolve and return rows (t, tau_c, D_c, Gamma_tau, max|M| in cells 0-2)."""
    install()
    kernel, st = initial_state(N, A, r_max=r_max, use_w=use_w)
    dt = cfl * float(st.grid.dr)
    n_every = max(1, int(round(every / dt)))
    rows, tau = [], 0.0
    for n in range(int(round(T / dt)) + 1):
        if n % n_every == 0:
            Dc = centre_value(st.grid, st.scalars.D); PDc = centre_value(st.grid, st.scalars.PD)
            rows.append((float(st.t), tau, Dc, PDc / Dc))
        a0 = centre_value(st.grid, st.geometry.alpha)
        st = kernel.step(st, dt)
        tau += 0.5 * dt * (a0 + centre_value(st.grid, st.geometry.alpha))
    return np.array(rows)


if __name__ == "__main__":
    import sys
    N = int(sys.argv[1]); T = float(sys.argv[2])
    for A in map(float, sys.argv[3:]):
        rows = run(N, A, T)
        for t, tau, Dc, G_tau in rows[:: max(1, len(rows) // 8)]:
            print(f"A={A} N={N} t={t:.2f} tau_c={tau:.3f} D_c={Dc:.4e} Gamma_tau={G_tau:.8f} "
                  f"(G/g-1)/Dc^2={(G_tau / GAMMA_LIN - 1) / Dc ** 2:.3f}", flush=True)
