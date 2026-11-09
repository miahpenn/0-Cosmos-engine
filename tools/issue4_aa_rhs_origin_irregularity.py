"""Issue #4: irregular (r^0) part at the origin of each Aa RHS piece, on regular evolved states."""
import sys
import numpy as np
sys.path.insert(0, ".")
from engine.production_kernel import V55ProductionKernel
from engine import v55_pirk_adapter as adapter
_, vac, _ = adapter.vendor_modules()

def c0(f, r, k=4):
    A = np.vstack([np.ones(k), r[:k]**2, r[:k]**4]).T
    return np.linalg.lstsq(A, np.asarray(f)[:k], rcond=None)[0][0]

for N in (120, 240, 480):
    k = V55ProductionKernel(); s = k.initialize(resolution=N, r_max=40.0)
    dt = 0.005
    for _ in range(int(round(0.3 / dt))):
        s = k.step(s, dt)
    g = s.geometry.copy(); g.alpha = k._solve_lapse_for_state(s.grid, g, s.scalars, s.matter, None)[0]
    t = k._geometry_stage_terms(s.grid, g, s.scalars, s.matter, None)
    r = np.asarray(s.grid.centers)
    l2, l3 = t["primary_l2"]["Aa"], t["primary_l3"]["Aa"]
    ex = t["explicit"]
    print(f"N={N}: Aa state c0 {c0(g.Aa, r):+.2e} | RHS c0: L2 {c0(l2, r):+.2e}  L3 {c0(l3, r):+.2e}  sum {c0(l2 + l3, r):+.2e} | "
          f"d_t(a-b) c0 {c0(ex['a'] - ex['b'], r):+.2e} | L2 cells0-3 {' '.join(f'{x:+.1e}' for x in l2[:4])}", flush=True)
