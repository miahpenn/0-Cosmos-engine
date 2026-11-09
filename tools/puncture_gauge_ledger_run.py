"""Misner-Sharp balance at r = 10 in the puncture-gauge research kernel with a moving S pulse."""
import sys, math, json; sys.path.insert(0, ".")
import numpy as np, dataclasses
from engine import v55_initial
from engine.matter_system import _metric_arrays, initialize_dust, initialize_radiation
from engine.v55_matter import V55MatterState, metric_slice_from_q
import engine.scalar_system as ss
from research.puncture_gauge_kernel_v1 import PunctureGaugeKernel
ss.V0 = 0.0
N, T, A, W, RMAX = int(sys.argv[1]), float(sys.argv[2]), 0.05, 3.0, 60.0
k = PunctureGaugeKernel(eta=2.0, ko_epsilon=float(sys.argv[4]) if len(sys.argv) > 4 else 0.0)
st = k.initialize(resolution=N, r_max=RMAX, include_radiation=False)
grid = st.grid; r = np.asarray(grid.centers); dr = float(grid.dr); z = np.zeros_like(r)
S = A * np.exp(-((r - 0.0) / W) ** 2)
Sp = v55_initial._d1(grid, S, 1); Bm = np.ones_like(r)
for _ in range(8):
    rho = 0.5 * Bm * Sp ** 2 + 0.5 * S ** 2
    Bm = v55_initial._cumulative_origin_integral(grid, r, dr, 1.0 - 8 * math.pi * r * r * rho) / r
Am = 1 / np.sqrt(Bm); g = st.geometry
g.a[:] = Am ** (4 / 3); g.b[:] = Am ** (-2 / 3); g.X[:] = Am ** (-1 / 3); g.K[:] = 0; g.Aa[:] = 0
d1 = lambda x: v55_initial._d1(grid, x, 1)
g.Lambda[:] = d1(g.a) / (2 * g.a * g.a) - d1(g.b) / (g.a * g.b) + 2 / r * (1 / g.b - 1 / g.a)
g.alpha[:] = 1.0; g.beta[:] = 0.0; g.B[:] = 0.0
mets = _metric_arrays(metric_slice_from_q(grid, g)); tiny = np.full_like(r, 1e-30)
st = dataclasses.replace(st, geometry=g, scalars=dataclasses.replace(st.scalars, S=S, PS=z.copy(), D=z.copy(), PD=z.copy(), phi=z.copy(), Pi=z.copy()),
    matter=V55MatterState(initialize_dust(mets, tiny), initialize_dust(mets, tiny), initialize_radiation(mets, tiny), None))
dt = 0.1 * dr; every = max(1, int(round(0.1 / dt))); rows = []
for n in range(int(round(T / dt)) + 1):
    if n % every == 0:
        d = k.diagnostics(st)
        rows.append([st.t, d["M_MS"], d["flux_T"], d["work_pR"], d["flux_shift"], float(st.geometry.beta[int(np.argmin(abs(r - 10)))])])
    st = k.step(st, dt)
a = np.array(rows); t, Mm, fl, wk, sh, bw = a.T
dM = np.gradient(Mm, t, edge_order=2); m = (t >= 1.0) & (t <= T - 1.0)
res_with = np.max(abs(dM - fl - wk - sh)[m]); res_without = np.max(abs(dM - fl - wk)[m])
print(f"N={N}: max|beta(10)| {np.max(abs(bw)):.2e}  max|flux| {np.max(abs(fl[m])):.2e} max|work| {np.max(abs(wk[m])):.2e} max|shift| {np.max(abs(sh[m])):.2e}  "
      f"residual with shift {res_with:.2e}  without {res_without:.2e}", flush=True)
json.dump(rows, open(sys.argv[3], "w"))
