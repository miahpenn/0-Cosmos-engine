"""RESEARCH ONLY. GEAR-67 lambda scan in the production engine (CMC, gravity on, no wall).
S Gaussian (amplitude S0, width 5, at rest) and a broad D seed (1e-4, width 15). Records alpha, S, D every 0.05 and
snapshots (with K) every 0.5; on an engine failure the history so far is kept.
Usage (repo root, COSMOS_PINNED_SOURCE_ROOT set): python research/sd_quartic_scan_v1.py lam S0 N T out.json
The saddle's vacuum is 0* (GEAR-91): S/D add nothing there and phi supplies it, at scale V0 (production 8.24e-5).
These runs set V0 = 0 literally (research choice, recorded): phi's supply is off, 0* resolves to plain 0 at the saddle,
and the ejected D region carries -1/(2 lambda)."""
import sys, math, json; sys.path.insert(0, ".")
import numpy as np, dataclasses
import engine.scalar_system as ss
from engine import v55_initial
from engine.matter_system import _metric_arrays, initialize_dust, initialize_radiation
from engine.v55_matter import V55MatterState, metric_slice_from_q
from engine.production_kernel import V55ProductionKernel
lam, S0, wS, D0, wD, N, T, RMAX = (float(sys.argv[1]), float(sys.argv[2]), 5.0, 1e-4, 15.0, int(sys.argv[3]), float(sys.argv[4]), 40.0)
ss.V0 = 0.0; ss.SD_QUARTIC_LAMBDA = lam
k = V55ProductionKernel(); st = k.initialize(resolution=N, r_max=RMAX, include_radiation=False)
grid = st.grid; r = np.asarray(grid.centers); dr = float(grid.dr); z = np.zeros_like(r)
S = S0 * np.exp(-(r / wS) ** 2); D = D0 * np.exp(-(r / wD) ** 2)
d1 = lambda x: v55_initial._d1(grid, x, 1); Sp, Dp = d1(S), d1(D); Bm = np.ones_like(r)
for _ in range(8):
    rho = 0.5 * Bm * (Sp ** 2 + Dp ** 2) + ss.sd_potential(S, D)
    Bm = v55_initial._cumulative_origin_integral(grid, r, dr, 1 - 8 * math.pi * r * r * rho) / r
Am = 1 / np.sqrt(Bm); g = st.geometry
g.a[:] = Am ** (4/3); g.b[:] = Am ** (-2/3); g.X[:] = Am ** (-1/3); g.K[:] = 0; g.Aa[:] = 0; g.beta[:] = 0; g.B[:] = 0
g.Lambda[:] = d1(g.a) / (2 * g.a**2) - d1(g.b) / (g.a * g.b) + 2 / r * (1 / g.b - 1 / g.a)
mets = _metric_arrays(metric_slice_from_q(grid, g)); tiny = np.full_like(r, 1e-30)
st = dataclasses.replace(st, geometry=g, scalars=dataclasses.replace(st.scalars, S=S, PS=z.copy(), D=D, PD=z.copy(), phi=z.copy(), Pi=z.copy()),
    matter=V55MatterState(initialize_dust(mets, tiny), initialize_dust(mets, tiny), initialize_radiation(mets, tiny), None))
g.alpha[:] = k._solve_lapse_for_state(grid, g, st.scalars, st.matter, None)[0]
dt = 0.06 * dr; every = int(round(0.5 / dt)); fine = max(1, int(round(0.05 / dt))); snaps = []; hist = []
for n in range(int(round(T / dt)) + 1):
    if n % fine == 0:
        hist.append({"t": st.t, "tau": st.tau, "alpha": st.geometry.alpha.tolist(), "S": st.scalars.S.tolist(), "D": st.scalars.D.tolist()})
    if n % every == 0:
        snaps.append({"t": st.t, "D": st.scalars.D.tolist(), "S": st.scalars.S.tolist(), "alpha0": float(st.geometry.alpha[0]),
                      "alpha": st.geometry.alpha.tolist(), "K": st.geometry.K.tolist()})
    try:
        st = k.step(st, dt)
    except Exception as exc:
        print(f"STOPPED at t={st.t:.3f}: {type(exc).__name__}: {str(exc)[:160]}", flush=True); break
json.dump({"r": r.tolist(), "lam": lam, "S0": S0, "N": N, "dt": dt, "snaps": snaps, "hist": hist}, open(sys.argv[5], "w"))
for s_ in snaps[:: max(1, len(snaps) // 8)]:
    Dd = np.abs(np.array(s_["D"])); print(f"lam={lam} t={s_['t']:.1f} |D|max {Dd.max():.2e} at r={r[Dd.argmax()]:.1f}  D(0) {s_['D'][0]:+.2e}  D(r=10) {np.interp(10, r, s_['D']):+.2e} D(r=20) {np.interp(20, r, s_['D']):+.2e}  alpha0 {s_['alpha0']:.4f} alpha min/max {min(s_['alpha']):.4f}/{max(s_['alpha']):.4f}", flush=True)
