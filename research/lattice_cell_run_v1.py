"""RESEARCH ONLY. Driver for runs in the lattice cell (research/lattice_cell_v1.py; HANDOFF Nov 14 2026 (4)).

One-budget rule (lambda = 5000, V0 and present densities x c^2 = 304.92, D's ejected vacuum is the floor), beta_D =
beta_S = 0, seed 1e-10. Records every `every` time units: the fields, the geometry (a, b, X, K, Aa, alpha, areal
radius, gamma_rr), the total stress-energy projections (G = 1 units), the Hamiltonian and momentum constraint residuals
(H = raw - 16 pi rho, M = raw - 8 pi j), dust dark-matter density, the cell's total baryon and dark-matter rest mass (conserved densities, which carry sqrt(gamma), summed x dr; baryons
have no source, so their total tests the wall's zero flux), proper cell volumes (per steradian: r^2 dr weights,
no 4 pi), and Kdot of the lapse closure.
Usage: python research/lattice_cell_run_v1.py CASE x_c N RMAX T out.npz [closure=ii] [start=cold] [cfl=0.2] [every=0.5]
       CASE = full | noD | nomatter
"""
import math
import sys

sys.path.insert(0, ".")

import numpy as np  # noqa: E402

import research.lattice_cell_v1 as lc  # noqa: E402
from engine import v55_pirk_adapter as adapter  # noqa: E402
from engine import worldtube as wt  # noqa: E402
from engine.stress_energy import assemble_total_stress_energy  # noqa: E402
from engine.v55_matter import dm_density, metric_slice_from_q  # noqa: E402

KEYS = ("t", "alpha", "R", "S", "PS", "D", "PD", "phi", "Pi", "K", "Aa", "a", "b", "X", "Lambda", "grr", "dV",
        "rho", "pr", "pt", "j", "Hc", "Mc", "rho_dm", "kdot", "bar_rest", "dm_rest")


def snapshot(k, st):
    grid, g = st.grid, st.geometry
    _, vacuum, _ = adapter.vendor_modules()
    tot = assemble_total_stress_energy(grid, g, st.scalars, st.matter)
    raw = vacuum.constraints(grid, g)
    vols = np.asarray(grid.volumes, float)
    kdot = k._solve_lapse_for_state(grid, g, st.scalars, st.matter, None)[1]
    return dict(t=st.t, alpha=g.alpha, R=wt.areal_radius(grid, g), S=st.scalars.S, PS=st.scalars.PS, D=st.scalars.D,
                PD=st.scalars.PD, phi=st.scalars.phi, Pi=st.scalars.Pi, K=g.K, Aa=g.Aa, a=g.a, b=g.b, X=g.X,
                Lambda=g.Lambda, grr=g.a / g.X**2, dV=np.sqrt(g.a) * g.b / g.X**3 * vols, rho=tot.rho, pr=tot.pr,
                pt=tot.pt, j=tot.j, Hc=raw["hamiltonian"] - 16.0 * math.pi * tot.rho,
                Mc=raw["momentum"] - 8.0 * math.pi * tot.j,
                rho_dm=dm_density(metric_slice_from_q(grid, g), st.matter), kdot=kdot,
                bar_rest=float(np.sum(np.asarray(st.matter.baryons.rest)) * grid.dr),
                dm_rest=float(np.sum(np.asarray(st.matter.dark_matter.rest)) * grid.dr))


def run(case, x_c, N, RMAX, T, out, closure="ii", start="cold", cfl=0.2, every=0.5, verbose=True):
    k, st, r, H0 = lc.cell_initial(case, x_c, 0.0, N, RMAX, closure=closure, start=start)
    dt = cfl * float(st.grid.dr)
    n_every = max(1, int(round(every / dt)))
    nsteps = int(round(T / dt))
    rec = {kk: [] for kk in KEYS}
    status = "ok"
    for n in range(nsteps + 1):
        if n % n_every == 0:
            s = snapshot(k, st)
            for kk in KEYS:
                rec[kk].append(np.array(s[kk], dtype=float).copy())
            if verbose and n % (n_every * 8) == 0:
                g = st.geometry
                print(f"t={st.t:6.2f} alpha0={g.alpha[0]:.4f} alpha_wall={g.alpha[-1]:.6f} K0={g.K[0]:+.5f} S0={st.scalars.S[0]:+.3e}"
                      f" |D|0={abs(st.scalars.D[0]):.2e} a/b-1(wall)={g.a[-1]/g.b[-1]-1:+.1e} maxH={np.abs(s['Hc']).max():.1e}", flush=True)
                np.savez_compressed(out, r=r, H0=H0, status=status, **{kk: np.array(v) for kk, v in rec.items()})
        if n < nsteps:
            try:
                st = k.step(st, dt)
            except Exception as exc:  # noqa: BLE001
                status = f"STOPPED at t={st.t:.3f}: {type(exc).__name__}: {str(exc)[:200]}"
                print(status, flush=True)
                break
    np.savez_compressed(out, r=r, H0=H0, status=status, case=case, x_c=x_c, closure=closure, start=start, cfl=cfl,
                        **{kk: np.array(v) for kk, v in rec.items()})
    return status


if __name__ == "__main__":
    a = sys.argv
    kw = dict(closure=a[7] if len(a) > 7 else "ii", start=a[8] if len(a) > 8 else "cold",
              cfl=float(a[9]) if len(a) > 9 else 0.2, every=float(a[10]) if len(a) > 10 else 0.5)
    print(f"[cell {a[1]} x_c={a[2]} N={a[3]} r_max={a[4]} T={a[5]} {kw}]", flush=True)
    print(run(a[1], float(a[2]), int(a[3]), float(a[4]), float(a[5]), a[6], **kw))
