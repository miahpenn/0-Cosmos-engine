"""RESEARCH ONLY. Self-held S pocket (single-harmonic oscillaton profile) in the production engine (CMC, no wall).
Usage: python research/sd_pocket_engine_v1.py osc.json D0 wD N RMAX T out.npz
S(t=0) = A(r) (turning point, PS = 0); D = D0 exp(-(r/wD)^2), PD = 0; metric from the Hamiltonian constraint
(polar-areal, K = 0). Records alpha, S, D, K every 0.1.
V0 = 0 switches off the phi sector's vacuum supply (research choice, recorded): the S/D potential adds nothing at the
saddle, so the vacuum there is whatever phi supplies (GEAR-91, 'V = 0*'); with phi's supply off, the ejected D region
carries -1/(2 lambda)."""
import sys, math, json; sys.path.insert(0, ".")
import numpy as np, dataclasses
import engine.scalar_system as ss
from engine import v55_initial
from engine.matter_system import _metric_arrays, initialize_dust, initialize_radiation
from engine.v55_matter import V55MatterState, metric_slice_from_q
from engine.production_kernel import V55ProductionKernel
osc = json.load(open(sys.argv[1])); D0, wD = float(sys.argv[2]), float(sys.argv[3])
N, RMAX, T, out = int(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6]), sys.argv[7]
lam = osc["lam"]; ss.V0 = 0.0; ss.SD_QUARTIC_LAMBDA = lam
k = V55ProductionKernel(); st = k.initialize(resolution=N, r_max=RMAX, include_radiation=False)
grid = st.grid; r = np.asarray(grid.centers); dr = float(grid.dr); z = np.zeros_like(r)
S = np.interp(r, osc["r"], osc["A"], right=0.0); D = D0 * np.exp(-(r / wD) ** 2)
d1 = lambda x: v55_initial._d1(grid, x, 1); Sp, Dp = d1(S), d1(D); Bm = np.ones_like(r)
for _ in range(12):
    rho = 0.5 * Bm * (Sp ** 2 + Dp ** 2) + ss.sd_potential(S, D)
    Bm = v55_initial._cumulative_origin_integral(grid, r, dr, 1 - 8 * math.pi * r * r * rho) / r
Am = 1 / np.sqrt(Bm); g = st.geometry
g.a[:] = Am ** (4/3); g.b[:] = Am ** (-2/3); g.X[:] = Am ** (-1/3); g.K[:] = 0; g.Aa[:] = 0; g.beta[:] = 0; g.B[:] = 0
g.Lambda[:] = d1(g.a) / (2 * g.a**2) - d1(g.b) / (g.a * g.b) + 2 / r * (1 / g.b - 1 / g.a)
mets = _metric_arrays(metric_slice_from_q(grid, g)); tiny = np.full_like(r, 1e-30)
st = dataclasses.replace(st, geometry=g, scalars=dataclasses.replace(st.scalars, S=S, PS=z.copy(), D=D, PD=z.copy(), phi=z.copy(), Pi=z.copy()),
    matter=V55MatterState(initialize_dust(mets, tiny), initialize_dust(mets, tiny), initialize_radiation(mets, tiny), None))
g.alpha[:] = k._solve_lapse_for_state(grid, g, st.scalars, st.matter, None)[0]
M0 = 0.5 * r[-1] * (1 - Bm[-1]); print(f"initial M(r_max) = {M0:.4f} (oscillaton M = {osc['M']:.4f})", flush=True)
dt = 0.06 * dr; fine = max(1, int(round(0.1 / dt))); rec = {"t": [], "tau": [], "alpha": [], "S": [], "D": [], "K": []}
nsteps = int(round(T / dt))
for n in range(nsteps + 1):
    if n % fine == 0:
        for key, v in (("t", st.t), ("tau", st.tau), ("alpha", st.geometry.alpha), ("S", st.scalars.S), ("D", st.scalars.D), ("K", st.geometry.K)):
            rec[key].append(np.array(v, dtype=float).copy())
        if n % (fine * 5) == 0:
            Dd = np.abs(st.scalars.D); Sd = np.abs(st.scalars.S)
            print(f"t={st.t:6.2f} a0={st.geometry.alpha[0]:.4f} amin/max={st.geometry.alpha.min():.4f}/{st.geometry.alpha.max():.4f} "
                  f"S0={st.scalars.S[0]:+.5f} |S|max={Sd.max():.5f} |D|max={Dd.max():.2e}@r={r[Dd.argmax()]:.1f} D(0)={st.scalars.D[0]:+.2e}", flush=True)
            np.savez_compressed(out, r=r, lam=lam, **{kk: np.array(vv) for kk, vv in rec.items()})
    if n < nsteps:
        st = k.step(st, dt)
np.savez_compressed(out, r=r, lam=lam, **{kk: np.array(vv) for kk, vv in rec.items()})
