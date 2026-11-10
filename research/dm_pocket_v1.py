"""RESEARCH ONLY. Dark-matter pocket experiment (HANDOFF-to-projects Nov 10, Task 3).

Uniform dark matter at the present-day density, phi at production values (V0 on: "at production V0"), a self-held
S oscillaton (research/sd_oscillaton_v1.py profile, S at its turning point) with a D seed inside, GEAR-67 quartic,
production kernel (CMC, beta = 0), no wall. Three couplings:
  off        GEAR-91 off (S, D interact with phi and dark matter by gravity only; dark matter still couples to phi, beta)
  on         GEAR-91 on, beta_S = 0, beta_D = beta (research/gear91_coupling_v1.py)
  nosource   GEAR-91 on but without phi's 8 pi dV_SD/dphi source (D's mass shifts with phi; nothing flows to phi):
             separates the D -> phi -> dark matter channel from the change in D's own dynamics.
Initial slice: time-symmetric (K = 0, polar-areal Hamiltonian constraint with all sources), not today's expanding
background; phi uniform at phi_today with Pi_today. Baryons and radiation are negligible (1e-30).

Unpaid inputs: lambda (GEAR-67 3(e)); V0 and phi_today (GEAR-92); beta_D choice. G_N = 1; G_reset not involved.

Usage (repo root): python research/dm_pocket_v1.py OSC.json off|on|nosource D0 wD N RMAX T out.npz
"""
import dataclasses
import json
import math
import sys

sys.path.insert(0, ".")

import numpy as np  # noqa: E402

import engine.scalar_system as ss  # noqa: E402
from engine import v55_initial  # noqa: E402
from engine.matter_system import _metric_arrays, initialize_dust, initialize_radiation  # noqa: E402
from engine.production_kernel import V55ProductionKernel  # noqa: E402
from engine.v55_matter import RHO_DM_PRESENT, V55MatterState, dm_density, metric_slice_from_q  # noqa: E402
from research import gear91_coupling_v1 as g91  # noqa: E402

PHI_TODAY, PI_TODAY = 0.179055, 0.00341913


def main():
    osc = json.load(open(sys.argv[1]))
    mode = sys.argv[2]
    D0, wD = float(sys.argv[3]), float(sys.argv[4])
    N, RMAX, T, out = int(sys.argv[5]), float(sys.argv[6]), float(sys.argv[7]), sys.argv[8]
    lam = osc["lam"]
    ss.SD_QUARTIC_LAMBDA = lam                      # V0 stays at its production value
    cfg = dict(off=None, on=dict(beta_S=0.0, beta_D=g91.BETA),
               nosource=dict(beta_S=0.0, beta_D=g91.BETA, phi_source=False))[mode]

    k = V55ProductionKernel()
    st = k.initialize(resolution=N, r_max=RMAX, include_radiation=False)
    grid = st.grid
    r = np.asarray(grid.centers); dr = float(grid.dr); z = np.zeros_like(r)
    S = np.interp(r, osc["r"], osc["A"], right=0.0)
    D = D0 * np.exp(-(r / wD) ** 2)
    phi = np.full_like(r, PHI_TODAY); Pi = np.full_like(r, PI_TODAY)
    rho_dm = np.full_like(r, RHO_DM_PRESENT)      # archive-normalised (8 pi rho_physical)

    if cfg is None:
        Vsd = ss.sd_potential(S, D)
    else:
        Vsd = g91.potential(S, D, phi, cfg["beta_S"], cfg["beta_D"], 1.0, lam)
    d1 = lambda x: v55_initial._d1(grid, x, 1)  # noqa: E731
    Sp, Dp = d1(S), d1(D)
    rho_phi = (0.5 * Pi**2 + ss.cosmos_potential(phi)) / ss.KAPPA
    Bm = np.ones_like(r)
    for _ in range(12):
        rho = 0.5 * Bm * (Sp**2 + Dp**2) + Vsd + rho_phi + rho_dm / ss.KAPPA
        Bm = v55_initial._cumulative_origin_integral(grid, r, dr, 1 - 8 * math.pi * r * r * rho) / r
    Am = 1 / np.sqrt(Bm); g = st.geometry
    g.a[:] = Am ** (4 / 3); g.b[:] = Am ** (-2 / 3); g.X[:] = Am ** (-1 / 3)
    g.K[:] = 0; g.Aa[:] = 0; g.beta[:] = 0; g.B[:] = 0
    g.Lambda[:] = d1(g.a) / (2 * g.a**2) - d1(g.b) / (g.a * g.b) + 2 / r * (1 / g.b - 1 / g.a)
    mets = _metric_arrays(metric_slice_from_q(grid, g)); tiny = np.full_like(r, 1e-30)
    st = dataclasses.replace(
        st, geometry=g,
        scalars=dataclasses.replace(st.scalars, S=S, PS=z.copy(), D=D, PD=z.copy(), phi=phi, Pi=Pi),
        matter=V55MatterState(initialize_dust(mets, rho_dm), initialize_dust(mets, tiny),
                              initialize_radiation(mets, tiny), None))

    import contextlib
    ctx = g91.coupled(**cfg) if cfg is not None else contextlib.nullcontext()
    with ctx:
        g.alpha[:] = k._solve_lapse_for_state(grid, g, st.scalars, st.matter, None)[0]
        print(f"[{mode}] initial M(r_max) = {0.5 * r[-1] * (1 - Bm[-1]):.4f}; lambda = {lam:g}; at production V0", flush=True)
        dt = 0.06 * dr; fine = max(1, int(round(0.5 / dt)))
        keys = ("t", "tau", "alpha", "S", "D", "phi", "Pi", "rho_dm")
        rec = {kk: [] for kk in keys}
        nsteps = int(round(T / dt))
        for n in range(nsteps + 1):
            if n % fine == 0:
                rdm = dm_density(metric_slice_from_q(grid, st.geometry), st.matter)
                for kk, v in (("t", st.t), ("tau", st.tau), ("alpha", st.geometry.alpha), ("S", st.scalars.S),
                              ("D", st.scalars.D), ("phi", st.scalars.phi), ("Pi", st.scalars.Pi), ("rho_dm", rdm)):
                    rec[kk].append(np.array(v, dtype=float).copy())
                if n % (fine * 4) == 0:
                    Dd = np.abs(st.scalars.D)
                    gg, kap, gk = g91.g_kappa(st.scalars.phi[0], beta_D=g91.BETA, beta_S=0.0)
                    print(f"t={st.t:6.2f} a0={st.geometry.alpha[0]:.4f} |D|max={Dd.max():.2e}@r={r[Dd.argmax()]:.1f} "
                          f"D(0)={st.scalars.D[0]:+.2e} phi(0)={st.scalars.phi[0]:.6f} rho_dm(0)={rdm[0]:.6e} "
                          f"g/kappa(0)={float(gk):.3f}", flush=True)
                    np.savez_compressed(out, r=r, lam=lam, mode=mode, **{kk: np.array(vv) for kk, vv in rec.items()})
            if n < nsteps:
                try:
                    st = k.step(st, dt)
                except Exception as exc:  # keep the history on an engine failure
                    print(f"STOPPED at t={st.t:.3f}: {type(exc).__name__}: {str(exc)[:160]}", flush=True)
                    break
    np.savez_compressed(out, r=r, lam=lam, mode=mode, **{kk: np.array(vv) for kk, vv in rec.items()})


if __name__ == "__main__":
    main()
