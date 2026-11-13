"""RESEARCH ONLY. Dark-matter pocket rerun under the one-budget rule (HANDOFF to projects Nov 11 (2), Task 4;
AMENDMENT-REV20-2 §4, which withdrew the earlier runs of research/dm_pocket_v1.py at production V0).

One-budget rule at lambda = 5000:
  - the hand floor is off (engine.scalar_system.D_POT = 0) and the floor is D's ejected vacuum at the run's own lambda
    (D = D* = sqrt(2/lambda) outside matter, energy -1/(2 lambda) at G_N = 1, i.e. -8 pi/(2 lambda) archive-normalised);
  - V0 and every present-day density (dark matter, baryons, radiation) are multiplied by
        c^2 = 8 pi / (2 lambda * 8.2425e-6) = 304.9,
    so that D's step equals the rescaled hand floor 0.10 c^2 V0 exactly (the identity of PR #469). Pi_today scales with H,
    by c = 17.46; phi_today is unchanged. H_today becomes 0.0990 (Hubble radius 10.1, one e-fold per 10.1 time units);
  - beta_D = 0 by default; -0.04 as a labelled variant (research/gear91_coupling_v1.py, beta_S = 0).
Matter pocket: a ball of S with x(r) = x_c (1 - tanh((r - R)/w))/2, x = lambda A^2, R = 8, w = 1, S = A at its turning
point (PS = 0). D: settled outside, relaxed with friction against the time-averaged S before the run (the cold start of
PR #470, done here in flat space on the production grid, so the energy released where D falls into its well is
removed), plus a 1e-10 seed inside. phi = phi_today, Pi = c Pi_today; dark matter, baryons, radiation uniform at the
rescaled present densities, at rest with the normal (comoving) observers.
Initial slice: production's expanding CMC slice (K = -3 H0 uniform; areal Hamiltonian reconstruction as in
engine/v55_initial.build_initial_data) with these profiles.
Cases:
  full       S pocket, D (cold start), beta_D as given
  noD        D = 0 everywhere; the floor restored as the hand floor at the same value (D_POT = 0.10 at c^2 V0, the
             identity), so the background is identical; S pocket present
  nomatter   no S; D = D* everywhere (one-budget floor); background only
Recorded every 0.5: t, alpha, areal R, S, PS, D, PD, phi, Pi, dark-matter density (archive-normalised), the dark-matter
conserved energy and rest mass per cell (times cell volume), proper volume per cell.

Unpaid inputs: lambda = 5000; V0 (rescaled by c^2, so its ratio to the floor is production's); phi_today = 0.179055
(imported "today", map row 10); beta_D = 0 or -0.04 (beta_S = 0); D's seed 1e-10; x_c; the cold start. G_N = 1;
G_reset not involved.
Usage (repo root): python research/dm_one_budget_v1.py CASE x_c beta_D N RMAX T out.npz
"""
import dataclasses
import math
import sys

sys.path.insert(0, ".")

import numpy as np  # noqa: E402

import engine.scalar_system as ss  # noqa: E402
import engine.v55_matter as vm  # noqa: E402
from engine import v55_initial  # noqa: E402
from engine import worldtube as wt  # noqa: E402
from engine.matter_system import _metric_arrays, initialize_dust, initialize_radiation  # noqa: E402
from engine.production_kernel import V55ProductionKernel  # noqa: E402
from engine.v55_matter import V55MatterState, dm_density, metric_slice_from_q  # noqa: E402
from research import gear91_coupling_v1 as g91  # noqa: E402

LAM = 5000.0
V0_PROD = 8.242522415500654e-5
C2 = 8 * math.pi / (2 * LAM * 8.2425e-6)
PHI_TODAY, PI_TODAY = 0.179055, 0.00341913
R_POCKET, W_EDGE = 8.0, 1.0


def cold_D(r, dr, x, T=200.0, gamma=0.3, seed=1e-10, beta_D=0.0, phi=PHI_TODAY):
    """D relaxed with friction against the time-averaged S (S^2 = A^2/2), flat space, cell-centred grid."""
    eD = math.exp(-2 * beta_D * phi)
    A2 = x / LAM
    D = np.sqrt(np.clip((2 / LAM) * (eD - 0.75 * x), 0, None))
    w = r * D; wp = w.copy(); dt = 0.4 * dr
    for _ in range(int(round(T / dt))):
        lap = np.empty_like(w)
        lap[1:-1] = (w[2:] - 2 * w[1:-1] + w[:-2]) / dr**2
        lap[0] = (w[1] - 3 * w[0]) / dr**2                       # odd ghost: u(-r0) = -u(r0)
        lap[-1] = (w[-2] - 2 * w[-1] + w[-1] * (r[-1] + dr) / r[-1]) / dr**2   # Neumann on D at r_max
        Dv = w / r
        aw = lap + w * (eD - 0.5 * LAM * Dv * Dv - 0.75 * LAM * A2)
        c = 0.5 * gamma * dt
        wn = (2 * w - (1 - c) * wp + dt * dt * aw) / (1 + c)
        wp, w = w, wn
    return w / r + seed * np.exp(-(r / (0.5 * R_POCKET)) ** 2)


def initial(case, x_c, beta_D, N, RMAX):
    ss.SD_QUARTIC_LAMBDA = LAM
    ss.V0 = C2 * V0_PROD
    ss.D_POT = 0.10 if case == "noD" else 0.0
    k = V55ProductionKernel()
    st = k.initialize(resolution=N, r_max=RMAX, include_radiation=False)
    grid = st.grid
    r = np.asarray(grid.centers, float); dr = float(grid.dr); z = np.zeros_like(r)
    x = x_c * 0.5 * (1 - np.tanh((r - R_POCKET) / W_EDGE)) if case != "nomatter" else z.copy()
    S = np.sqrt(x / LAM)
    if case == "noD":
        D = z.copy()
    elif case == "nomatter":
        D = np.full_like(r, math.sqrt(2 / LAM) * math.exp(-beta_D * PHI_TODAY))
    else:
        D = cold_D(r, dr, x, beta_D=beta_D)
    phi = np.full_like(r, PHI_TODAY); Pi = np.full_like(r, math.sqrt(C2) * PI_TODAY)
    dm0, b0, rr0 = C2 * vm.RHO_DM_PRESENT, C2 * vm.RHO_B_PRESENT, C2 * vm.RHO_R_PRESENT
    # background H0 from the far-field budget (archive-normalised; S = 0, D at its far value)
    Dfar = D[-1]
    Vfar = 8 * math.pi * float(g91.potential(0.0, Dfar, PHI_TODAY, 0.0, beta_D, 1.0, LAM))
    rho0 = 0.5 * Pi[0] ** 2 + float(ss.cosmos_potential(np.array([PHI_TODAY]))[0]) + Vfar + dm0 + b0 + rr0
    H0 = math.sqrt(rho0 / 3.0)
    Kt = np.full_like(r, -H0)
    d1 = lambda f: v55_initial._d1(grid, f, 1)  # noqa: E731
    Sp, Dp = d1(S), d1(D)
    B = np.ones_like(r)
    for _ in range(8):
        rho = (0.5 * B * (Sp**2 + Dp**2) + g91.potential(S, D, phi, 0.0, beta_D, 1.0, LAM)
               + (0.5 * Pi**2 + ss.cosmos_potential(phi)) / ss.KAPPA + (dm0 + b0 + rr0) / ss.KAPPA)
        src = 1.0 - 8 * math.pi * r * r * rho + 3.0 * r * r * Kt * Kt
        B = v55_initial._cumulative_origin_integral(grid, r, dr, src) / r
    if np.any(B <= 0) or not np.all(np.isfinite(B)):
        raise FloatingPointError("invalid constrained initial metric")
    Am = 1 / np.sqrt(B); g = st.geometry
    g.a[:] = Am ** (4 / 3); g.b[:] = Am ** (-2 / 3); g.X[:] = Am ** (-1 / 3)
    g.K[:] = 3 * Kt; g.Aa[:] = 0; g.beta[:] = 0; g.B[:] = 0
    g.Lambda[:] = d1(g.a) / (2 * g.a**2) - d1(g.b) / (g.a * g.b) + 2 / r * (1 / g.b - 1 / g.a)
    g.B[:] = 0.75 * g.Lambda
    mets = _metric_arrays(metric_slice_from_q(grid, g))
    st = dataclasses.replace(
        st, geometry=g,
        scalars=dataclasses.replace(st.scalars, S=S, PS=z.copy(), D=D, PD=z.copy(), phi=phi, Pi=Pi),
        matter=V55MatterState(initialize_dust(mets, np.full_like(r, dm0)), initialize_dust(mets, np.full_like(r, b0)),
                              initialize_radiation(mets, np.full_like(r, rr0)), None))
    return k, st, r, H0


def main():
    case, x_c, beta_D = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
    N, RMAX, T, out = int(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6]), sys.argv[7]
    import contextlib
    ctx = g91.coupled(beta_S=0.0, beta_D=beta_D) if beta_D != 0.0 else contextlib.nullcontext()
    with ctx:
        k, st, r, H0 = initial(case, x_c, beta_D, N, RMAX)
        grid = st.grid
        st.geometry.alpha[:] = k._solve_lapse_for_state(grid, st.geometry, st.scalars, st.matter, None)[0]
        print(f"[{case} x_c={x_c} beta_D={beta_D}] lambda={LAM:g} c^2={C2:.2f} V0={ss.V0:.4e} D_POT={ss.D_POT} H0={H0:.4f}",
              flush=True)
        dr = float(grid.dr); dt = 0.06 * dr; every = max(1, int(round(0.5 / dt)))
        keys = ("t", "alpha", "R", "S", "PS", "D", "PD", "phi", "Pi", "K", "rho_dm", "dm_E", "dm_rest", "dV")
        rec = {kk: [] for kk in keys}
        nsteps = int(round(T / dt))
        vols = np.asarray(grid.volumes, float)
        for n in range(nsteps + 1):
            if n % every == 0:
                g = st.geometry
                sqrtg = np.sqrt(g.a) * g.b / g.X**3          # sqrt(gamma)/(r^2 sin) for gamma_rr = a/X^2, gamma_TT = r^2 b/X^2
                rdm = dm_density(metric_slice_from_q(grid, g), st.matter)
                dmst = st.matter.dark_matter
                for kk, v in (("t", st.t), ("alpha", g.alpha), ("R", wt.areal_radius(grid, g)), ("S", st.scalars.S),
                              ("PS", st.scalars.PS), ("D", st.scalars.D), ("PD", st.scalars.PD),
                              ("phi", st.scalars.phi), ("Pi", st.scalars.Pi), ("K", g.K), ("rho_dm", rdm),
                              ("dm_E", np.asarray(dmst.energy_t) * vols), ("dm_rest", np.asarray(dmst.rest) * vols),
                              ("dV", sqrtg * vols)):
                    rec[kk].append(np.array(v, dtype=float).copy())
                if n % (every * 4) == 0:
                    print(f"t={st.t:6.2f} a0={g.alpha[0]:.4f} K0={g.K[0]:+.4f} S(0)={st.scalars.S[0]:+.3e} "
                          f"|D|(0)={abs(st.scalars.D[0]):.2e} D(r=15)={st.scalars.D[np.searchsorted(r,15)]:.4e} "
                          f"rho_dm(0)/far={rdm[0]/rdm[-5]:.4f}", flush=True)
                    np.savez_compressed(out, r=r, H0=H0, **{kk: np.array(vv) for kk, vv in rec.items()})
            if n < nsteps:
                try:
                    st = k.step(st, dt)
                except Exception as exc:
                    print(f"STOPPED at t={st.t:.3f}: {type(exc).__name__}: {str(exc)[:200]}", flush=True)
                    break
        np.savez_compressed(out, r=r, H0=H0, **{kk: np.array(vv) for kk, vv in rec.items()})


if __name__ == "__main__":
    main()
