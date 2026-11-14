"""RESEARCH ONLY. Driver for runs in the lattice cell (research/lattice_cell_v1.py; HANDOFF Nov 14 2026 (4)).

One-budget rule (lambda = 5000, V0 and present densities x c^2 = 304.92, D's ejected vacuum is the floor), beta_D =
beta_S = 0, seed 1e-10. Records every `every` time units: the fields, the geometry (a, b, X, K, Aa, alpha, areal
radius, gamma_rr), the total stress-energy projections (G = 1 units), the Hamiltonian and momentum constraint residuals
(H = raw - 16 pi rho, M = raw - 8 pi j), dust dark-matter density, the cell's total baryon and dark-matter rest mass (conserved densities, which carry sqrt(gamma), summed x dr; baryons
have no source, so their total tests the wall's zero flux), proper cell volumes (per steradian: r^2 dr weights,
no 4 pi), and Kdot of the lapse closure. Every second step: the fine time series `ts` (columns TS_KEYS) for the
settling clock and S's local frequency: S projected on its initial profile, the pumping rates P1, P2, D's kinetic
energy, and the sums of D^2, S^2, PS^2 over the pocket (r < 11) and over the whole cell, with the centre's and the
wall's proper time.
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


TS_KEYS = ("t", "s", "S0", "P1", "P2", "KE_D", "D2", "S2", "PS2", "vol", "P1_cell", "KE_D_cell", "D2_cell", "S2_cell",
           "PS2_cell", "vol_cell", "alpha0", "alpha_wall", "tau0", "tau_wall")


def fine_series(st, prof, pk, vols, tau):
    """Pocket (r < 11) and whole-cell sums for the settling clock and S's local frequency (proper volume, per sr)."""
    g, sc = st.geometry, st.scalars
    dV = np.sqrt(g.a) * g.b / g.X**3 * vols
    S, PS, D, PD = (np.asarray(x) for x in (sc.S, sc.PS, sc.D, sc.PD))
    p1 = -PD * 1.5 * LAM * S * S * D * dV              # rate of energy into D's own sector, per unit proper time
    p2 = 1.5 * LAM * S * PS * D * D * dV               # rate of energy out of S's own sector
    ke = 0.5 * PD * PD * dV
    f = lambda q: float(q[pk].sum())  # noqa: E731
    return (st.t, float(np.sum((S * prof * dV)[pk]) / np.sum((prof * prof * dV)[pk])) if prof is not None else 0.0, S[0],
            f(p1), f(p2), f(ke), f(D * D * dV), f(S * S * dV), f(PS * PS * dV), f(dV),
            float(p1.sum()), float(ke.sum()), float((D * D * dV).sum()), float((S * S * dV).sum()),
            float((PS * PS * dV).sum()), float(dV.sum()), g.alpha[0], g.alpha[-1], tau[0], tau[1])


LAM = 5000.0


def run(case, x_c, N, RMAX, T, out, closure="ii", start="cold", cfl=0.2, every=0.5, verbose=True):
    k, st, r, H0 = lc.cell_initial(case, x_c, 0.0, N, RMAX, closure=closure, start=start)
    dt = cfl * float(st.grid.dr)
    prof = np.array(st.scalars.S, dtype=float).copy() if np.any(np.asarray(st.scalars.S) != 0) else None
    pk = r < 11.0; vols_ = np.asarray(st.grid.volumes, float); ts = []; tau = [0.0, 0.0]
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
                np.savez_compressed(out, r=r, H0=H0, status=status, case=case, x_c=x_c, closure=closure, start=start,
                                    cfl=cfl, ts=np.array(ts), **{kk: np.array(v) for kk, v in rec.items()})
        if n % 2 == 0:
            ts.append(fine_series(st, prof, pk, vols_, tau))
        if n < nsteps:
            try:
                a0, aw = st.geometry.alpha[0], st.geometry.alpha[-1]
                st = k.step(st, dt)
                tau[0] += 0.5 * (a0 + st.geometry.alpha[0]) * dt; tau[1] += 0.5 * (aw + st.geometry.alpha[-1]) * dt
            except Exception as exc:  # noqa: BLE001
                status = f"STOPPED at t={st.t:.3f}: {type(exc).__name__}: {str(exc)[:200]}"
                print(status, flush=True)
                break
    np.savez_compressed(out, r=r, H0=H0, status=status, case=case, x_c=x_c, closure=closure, start=start, cfl=cfl,
                        ts=np.array(ts), **{kk: np.array(v) for kk, v in rec.items()})
    return status


if __name__ == "__main__":
    a = sys.argv
    kw = dict(closure=a[7] if len(a) > 7 else "ii", start=a[8] if len(a) > 8 else "cold",
              cfl=float(a[9]) if len(a) > 9 else 0.2, every=float(a[10]) if len(a) > 10 else 0.5)
    print(f"[cell {a[1]} x_c={a[2]} N={a[3]} r_max={a[4]} T={a[5]} {kw}]", flush=True)
    print(run(a[1], float(a[2]), int(a[3]), float(a[4]), float(a[5]), a[6], **kw))
