"""RESEARCH ONLY. Dynamic bag test in the production engine, gravity on (HANDOFF Nov 11, Task 2).

A Thomas-Fermi bag pocket (research/pocket_with_vacuum_tov_v1.py: surface at P_S = DeltaV, x = 3.266) evolved by the
production kernel (CMC, beta = 0), GEAR-67 quartic at lambda = 5000, S at its turning point (PS = 0), D at its ejected
value D* outside the surface and a 1e-10 seed inside. phi, V0 and the fluids are off (V0 = 0, phi = 0, matter 1e-30).

Vacuum zero. With D ejected outside, V_SD = -DeltaV there, and over the domain that negative vacuum would outweigh the
pocket (the collapse seen in earlier runs). mode "bag" therefore adds a uniform +DeltaV to the S/D energy (and -DeltaV to
its pressures): the total vacuum is 0 outside and +DeltaV inside the pocket. This is phi's vacuum frozen at v = 1
("V_total = 0 outside", phi = 2.303 under the identity), stated as a choice; the TOV profile is solved with the same v.
mode "nod" is the control: the same matter, no D outside (seed only), no offset.

Unpaid inputs: lambda = 5000; the frozen vacuum v = 1 stands in for phi (no V0, phi_today or beta_D enters). G_N = 1.
Usage (repo root): python research/sd_bag_engine_v1.py bag|nod x_c N RMAX T out.npz
"""
import dataclasses
import math
import sys

sys.path.insert(0, ".")

import numpy as np  # noqa: E402
from scipy.integrate import solve_ivp  # noqa: E402

import engine.scalar_system as ss  # noqa: E402
import engine.v55_matter as vm  # noqa: E402
from engine import v55_initial  # noqa: E402
from engine.matter_system import _metric_arrays, initialize_dust, initialize_radiation  # noqa: E402
from engine.production_kernel import V55ProductionKernel  # noqa: E402
from engine.v55_matter import V55MatterState, metric_slice_from_q  # noqa: E402

LAM = 5000.0
DV = 1.0 / (2 * LAM)
DSTAR = math.sqrt(2.0 / LAM)
X_SURF = math.sqrt(32.0 / 3.0)


def tov_profile(xc, v):
    """x(r) inside the bag, physical units (r in 1/mu)."""
    def rhs(r, y):
        m, x = y
        P, rho = 3 * x * x / 32, x * (1 + 9 * x / 32)
        dP = -(rho + P) * (m + 4 * math.pi * r**3 * (P - v)) / (r * (r - 2 * m))
        return [4 * math.pi * r * r * (rho + v), dP / (3 * x / 16)]
    surf = lambda r, y: y[1] - X_SURF  # noqa: E731
    surf.terminal = True
    r0 = 1e-5
    s = solve_ivp(rhs, (r0, 5.0), [4 * math.pi / 3 * r0**3 * (xc * (1 + 9 * xc / 32) + v), xc], events=surf, rtol=1e-10, atol=1e-13,
                  max_step=1e-3)
    return s.t / math.sqrt(DV), s.y[1]


def main():
    mode, xc = sys.argv[1], float(sys.argv[2])
    N, RMAX, T, out = int(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), sys.argv[6]
    C = DV if mode == "bag" else 0.0
    ss.V0 = 0.0; ss.SD_QUARTIC_LAMBDA = LAM
    rp, xp = tov_profile(xc, 1.0 if mode == "bag" else 0.0)
    R = rp[-1]

    proj0 = vm.scalar_projection

    def projection(grid, geometry, fields):
        e, pr, pt, j = proj0(grid, geometry, fields)
        return e + C, pr - C, pt - C, j
    vm.scalar_projection = projection

    k = V55ProductionKernel()
    st = k.initialize(resolution=N, r_max=RMAX, include_radiation=False)
    grid = st.grid; r = np.asarray(grid.centers); dr = float(grid.dr); z = np.zeros_like(r)
    A_in = np.sqrt(np.interp(r, rp, xp, left=xp[0], right=xp[-1]) / LAM)
    A_s = math.sqrt(X_SURF / LAM)
    S = np.where(r <= R, A_in, A_s * np.exp(-2.0 * (r - R)))          # outside S has mass 2: skin depth 1/2
    D = 1e-10 * np.exp(-(r / (0.5 * R)) ** 2)
    if mode == "bag":
        D = D + DSTAR * 0.5 * (1 + np.tanh((r - R - 0.5) / 0.7))
    d1 = lambda f: v55_initial._d1(grid, f, 1)  # noqa: E731
    Sp, Dp = d1(S), d1(D); Bm = np.ones_like(r)
    for _ in range(14):
        rho = 0.5 * Bm * (Sp**2 + Dp**2) + ss.sd_potential(S, D) + C
        Bm = v55_initial._cumulative_origin_integral(grid, r, dr, 1 - 8 * math.pi * r * r * rho) / r
    Am = 1 / np.sqrt(Bm); g = st.geometry
    g.a[:] = Am ** (4 / 3); g.b[:] = Am ** (-2 / 3); g.X[:] = Am ** (-1 / 3); g.K[:] = 0; g.Aa[:] = 0; g.beta[:] = 0; g.B[:] = 0
    g.Lambda[:] = d1(g.a) / (2 * g.a**2) - d1(g.b) / (g.a * g.b) + 2 / r * (1 / g.b - 1 / g.a)
    mets = _metric_arrays(metric_slice_from_q(grid, g)); tiny = np.full_like(r, 1e-30)
    st = dataclasses.replace(st, geometry=g,
                             scalars=dataclasses.replace(st.scalars, S=S, PS=z.copy(), D=D, PD=z.copy(), phi=z.copy(), Pi=z.copy()),
                             matter=V55MatterState(initialize_dust(mets, tiny), initialize_dust(mets, tiny),
                                                   initialize_radiation(mets, tiny), None))
    g.alpha[:] = k._solve_lapse_for_state(grid, g, st.scalars, st.matter, None)[0]
    print(f"[{mode}] x_c = {xc}, R = {R:.3f}, M(r_max) = {0.5*r[-1]*(1-Bm[-1]):.4f}, 2m/R at surface = {(1-np.interp(R, r, Bm)):.3f}, "
          f"lambda = {LAM:g}, vacuum offset = {C:.1e}", flush=True)
    dt = 0.06 * dr; fine = max(1, int(round(0.1 / dt)))
    rec = {kk: [] for kk in ("t", "tau", "alpha", "S", "D")}
    nsteps = int(round(T / dt))

    def save():
        np.savez_compressed(out, r=r, R=R, xc=xc, mode=mode, **{kk: np.array(vv) for kk, vv in rec.items()})
    for n in range(nsteps + 1):
        if n % fine == 0:
            for kk, val in (("t", st.t), ("tau", st.tau), ("alpha", st.geometry.alpha), ("S", st.scalars.S), ("D", st.scalars.D)):
                rec[kk].append(np.array(val, dtype=float).copy())
            if n % (fine * 10) == 0:
                Dd = np.abs(st.scalars.D)
                print(f"t={st.t:6.2f} a0={st.geometry.alpha[0]:.4f} amin/max={st.geometry.alpha.min():.3f}/{st.geometry.alpha.max():.3f} "
                      f"x(0)={LAM*st.scalars.S[0]**2:.3f} D(0)/D*={st.scalars.D[0]/DSTAR:+.2e} |D|max/D*={Dd.max()/DSTAR:.3f}", flush=True)
                save()
        if n < nsteps:
            try:
                st = k.step(st, dt)
            except Exception as exc:
                print(f"STOPPED at t={st.t:.3f}: {type(exc).__name__}: {str(exc)[:160]}", flush=True)
                break
    save()


if __name__ == "__main__":
    main()
