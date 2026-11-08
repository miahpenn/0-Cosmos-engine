"""Centre momentum-constraint floor (Issue #4): temporal vs spatial separation.

Run from a checkout of main (CMC centre row #436 and initializer origin quadrature #439 in place).
Default initial data, r_max = 40, evolve to t = 0.60 with real kernel.step calls.
Records total M (M - 8 pi j) at cells 0-5 and the 5-10 annulus, plus H and C, at t = 0 and t = 0.60.
Study T: N fixed, dt = c dr for c in (0.06, 0.03, 0.015).
Study S: dt fixed at 0.005, N in (60, 120, 240, 480).
Usage: python <this file> <study T|S> <N or dt-list> output.json
"""
import json, math, sys, time
import numpy as np
sys.path.insert(0, ".")
from engine.production_kernel import V55ProductionKernel
from engine import v55_pirk_adapter as adapter
_, vac, _ = adapter.vendor_modules()

def measure(k, s):
    c = vac.constraints(s.grid, s.geometry); d = k.diagnostics(s, profiles=True); p = d.get("profiles", d)
    H = np.asarray(c["hamiltonian"]) - 16 * math.pi * np.asarray(p["total_rho"])
    M = np.asarray(c["momentum"]) - 8 * math.pi * np.asarray(p["total_j"])
    C = np.asarray(c["connection"]); r = np.asarray(s.grid.centers); ann = (r >= 5) & (r <= 10)
    return {"M_cells_0_5": [float(x) for x in M[:6]], "r_0_5": [float(x) for x in r[:6]],
            **{f"{q}_{w}": float(np.max(np.abs(v[sl]))) for q, v in (("M", M), ("H", H), ("C", C))
               for w, sl in (("c02", slice(0, 3)), ("ann", ann))}}

def run(n, dt):
    k = V55ProductionKernel(); s = k.initialize(resolution=n, r_max=40.0)
    out = {"N": n, "dt": dt, "t0": measure(k, s)}; t0 = time.time()
    steps = int(round(0.6 / dt))
    for _ in range(steps):
        s = k.step(s, dt)
    out.update({"t_end": float(s.t), "steps": steps, "t0.60": measure(k, s), "wall_s": time.time() - t0})
    m = out["t0.60"]
    print(f"N={n:>3} dt={dt:.5f} steps={steps}: M c0-2 {m['M_c02']:.3e} ann {m['M_ann']:.2e} | H c0-2 {m['H_c02']:.2e} | C c0-2 {m['C_c02']:.2e} | M0..2 {[f'{x:.2e}' for x in m['M_cells_0_5'][:3]]} ({out['wall_s']:.0f}s)", flush=True)
    return out

study = sys.argv[1]; res = []
if study == "T":
    n = int(sys.argv[2]); dr = 40.0 / n
    for c in (0.06, 0.03, 0.015):
        res.append(run(n, c * dr))
else:
    for n in [int(x) for x in sys.argv[2].split(",")]:
        res.append(run(n, 0.005))
json.dump(res, open(sys.argv[3], "w"), indent=1)
