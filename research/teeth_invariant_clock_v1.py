"""RESEARCH ONLY. Teeth test on a slicing-independent clock (HANDOFF-to-projects Nov 10, Task 4; backlog item 6).

AMENDMENT-2-GEAR-68 asked whether the CMC lapse is a fixed function of D's local energy (eta^2 ~ 0.85, cross-run
residual ~0.55 in ln alpha). The lapse is a slicing statement. In spherical symmetry the slicing-independent clock rate
is the norm of the Kodama vector, |K|^2 = -(grad R)^2 = -(1 - 2m/R) (engine.worldtube.invariant_chi), with m the
Misner-Sharp mass and R the areal radius: the
redshift factor N_K = sqrt|1 - 2m/R| relative to Kodama time (the Killing time wherever the region is static).
It is built from invariants only (m, R), so it does not depend on the slicing.

Runs (`run`): AMENDMENT-GEAR-68's no-wall set rebuilt on the production engine (CMC, beta = 0, lambda = 0, V0 = 0,
phi = 0, matter at floor): S = A exp(-(r/w)^2), D = eps exp(-(r/w)^2), time-symmetric, N = 160, r_max = 80, T = 45,
w = 7, eps = 1e-10, A = 0.005, 0.01, 0.02. Profiles every 1.0: t, r, S, PS, D, PD, alpha, R, m.
Points: D-dominated (|D| > 1e-3, |D| > 5|S|, r < 60) and untrapped (1 - 2m/R > 0, where the Kodama vector is timelike;
about 4% of D-dominated points are trapped and excluded). Analysis (`teeth`): AMENDMENT-2-GEAR-68's scores (eta^2 on 30 quantile bins, cross-run RMS error) for ln alpha (CMC,
reproduction) and ln N_K (invariant), against local D^2, act_D = PD^2 + D^2, |rho_D|, the enclosed active source,
and the controls r, t; plus m/R itself.

Unpaid inputs: none beyond the production lambda = 0 and V0 = 0 (phi off). G_N = 1; G_reset not involved.
Usage (repo root): python research/teeth_invariant_clock_v1.py run A out.npz  |  teeth ref.npz other.npz ...
"""
import dataclasses
import math
import sys

sys.path.insert(0, ".")

import numpy as np  # noqa: E402


def run(A, out, N=160, r_max=80.0, T=45.0, w=7.0, eps=1e-10):
    import engine.production_kernel as pk
    import engine.scalar_system as ss
    from engine import v55_initial
    from engine import worldtube as wt
    from engine.matter_system import _metric_arrays, initialize_dust, initialize_radiation
    from engine.v55_matter import V55MatterState, metric_slice_from_q

    ss.V0 = 0.0
    kernel = pk.V55ProductionKernel()
    st = kernel.initialize(resolution=N, r_max=r_max, include_radiation=False)
    grid = st.grid; r = np.asarray(grid.centers, float); dr = float(grid.dr)
    gss = np.exp(-(r / w) ** 2); S = A * gss; D = eps * gss; z = np.zeros_like(r)
    d1 = lambda x: v55_initial._d1(grid, x, 1)  # noqa: E731
    Sp, Dp = d1(S), d1(D); Bm = np.ones_like(r)
    for _ in range(10):
        rho = 0.5 * Bm * (Sp**2 + Dp**2) + 0.5 * S**2 - 0.5 * D**2
        Bm = v55_initial._cumulative_origin_integral(grid, r, dr, 1.0 - 8.0 * math.pi * r * r * rho) / r
    Am = 1.0 / np.sqrt(Bm); geo = st.geometry
    geo.a[:] = Am ** (4 / 3); geo.b[:] = Am ** (-2 / 3); geo.X[:] = Am ** (-1 / 3); geo.K[:] = 0; geo.Aa[:] = 0
    geo.Lambda[:] = d1(geo.a) / (2 * geo.a**2) - d1(geo.b) / (geo.a * geo.b) + 2.0 / r * (1 / geo.b - 1 / geo.a)
    geo.beta[:] = 0; geo.B[:] = 0
    st = dataclasses.replace(st, scalars=dataclasses.replace(st.scalars, S=S, PS=z.copy(), D=D, PD=z.copy(),
                                                            phi=z.copy(), Pi=z.copy()))
    met = _metric_arrays(metric_slice_from_q(grid, geo)); tiny = np.full_like(r, 1e-30)
    st = dataclasses.replace(st, matter=V55MatterState(initialize_dust(met, tiny.copy()), initialize_dust(met, tiny.copy()),
                                                       initialize_radiation(met, tiny.copy()), None))
    geo.alpha[:] = kernel._solve_lapse_for_state(grid, geo, st.scalars, st.matter, None)[0]
    dt = 0.06 * dr; nsteps = int(round(T / dt)); every = max(1, int(round(1.0 / dt)))
    keys = ("t", "S", "PS", "D", "PD", "alpha", "R", "m")
    rec = {k: [] for k in keys}
    for n in range(nsteps + 1):
        if n % every == 0:
            R = np.asarray(wt.areal_radius(grid, st.geometry))
            m = np.asarray(wt.misner_sharp(grid, st.geometry))
            for k, v in (("t", st.t), ("S", st.scalars.S), ("PS", st.scalars.PS), ("D", st.scalars.D),
                         ("PD", st.scalars.PD), ("alpha", st.geometry.alpha), ("R", R), ("m", m)):
                rec[k].append(np.array(v, dtype=float).copy())
        if n < nsteps:
            try:
                st = kernel.step(st, dt)
            except Exception as exc:
                print(f"STOPPED at t={st.t:.3f}: {type(exc).__name__}: {str(exc)[:160]}"); break
    np.savez_compressed(out, r=r, A=A, **{k: np.array(v) for k, v in rec.items()})
    print(f"A={A}: saved {len(rec['t'])} profiles to t={rec['t'][-1]:.2f}")


def _points(path):
    d = np.load(path); r = d["r"]; dr = r[1] - r[0]
    rows = []
    for i, t in enumerate(d["t"]):
        S, PS, D, PD, al, R, m = (d[k][i] for k in ("S", "PS", "D", "PD", "alpha", "R", "m"))
        Dp, Sp = np.gradient(D, dr), np.gradient(S, dr)
        rhoD = 0.5 * PD**2 + 0.5 * Dp**2 - 0.5 * D**2
        act = PD**2 + D**2; actS = PS**2 + S**2
        Mact = np.cumsum(4 * np.pi * r**2 * (act + actS)) * dr
        chi = 1.0 - 2.0 * m / R
        NK = np.sqrt(np.abs(chi))
        sel = (r < 60) & (np.abs(D) > 1e-3) & (np.abs(D) > 5 * np.abs(S)) & (chi > 0.0)   # untrapped: Kodama timelike
        for row in zip(np.log(al[sel]), np.log(NK[sel]), np.abs(rhoD[sel]), act[sel], D[sel] ** 2, Mact[sel],
                       r[sel], np.full(sel.sum(), t), chi[sel]):
            rows.append(row)
    return np.array(rows)


def _binned(x, y, nb=30):
    lx = np.log10(x + 1e-9); edges = np.quantile(lx, np.linspace(0, 1, nb + 1)); edges[-1] += 1e-12
    idx = np.clip(np.digitize(lx, edges) - 1, 0, nb - 1)
    mu = np.array([y[idx == k].mean() if np.any(idx == k) else np.nan for k in range(nb)])
    within = sum(((y[idx == k] - mu[k]) ** 2).sum() for k in range(nb) if np.any(idx == k))
    return 1 - within / ((y - y.mean()) ** 2).sum(), edges, mu


def _predict(edges, mu, x):
    lx = np.log10(x + 1e-9); c = 0.5 * (edges[1:] + edges[:-1]); ok = ~np.isnan(mu)
    return np.interp(lx, c[ok], mu[ok])


NAMES = ["|rho_D| (local energy)", "act_D (local active)", "D^2 (local)", "M_act(r) (enclosed, non-local)",
         "r (position; control)", "t (time; control)", "1 - 2m/R (enclosed; N_K's definition)"]
COLS = [2, 3, 4, 5, 6, 7, 8]


def teeth(files):
    data = {f: _points(f) for f in files}
    for f, d in data.items():
        print(f"{f}: {len(d)} points; ln alpha {d[:, 0].min():.2f}..{d[:, 0].max():.2f}; "
              f"ln N_K {d[:, 1].min():.3f}..{d[:, 1].max():.3f}")
    for yc, label in ((0, "ln alpha (CMC lapse; reproduction)"), (1, "ln N_K (Kodama norm; slicing-independent)")):
        print(f"\n== {label}: eta^2 within each run")
        for nm, c in zip(NAMES, COLS):
            print(f"  {nm:34s}" + " ".join(f"{_binned(d[:, c], d[:, yc])[0]:.3f}" for d in data.values()))
        ref = files[0]
        print(f"   cross-run RMS error (curve from {ref}) [spread]:")
        for nm, c in zip(NAMES, COLS):
            _, edges, mu = _binned(data[ref][:, c], data[ref][:, yc])
            row = []
            for f in files[1:]:
                d = data[f]; lx = np.log10(d[:, c] + 1e-9); ok = (lx >= edges[0]) & (lx <= edges[-1])
                err = np.sqrt(np.mean((_predict(edges, mu, d[ok, c]) - d[ok, yc]) ** 2))
                row.append(f"{err:.3f} [{d[:, yc].std():.3f}]")
            print(f"  {nm:34s}" + "  ".join(row))


if __name__ == "__main__":
    if sys.argv[1] == "run":
        run(float(sys.argv[2]), sys.argv[3])
    else:
        teeth(sys.argv[2:])
