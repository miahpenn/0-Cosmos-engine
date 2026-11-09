"""Issue #4 original symptom: inner momentum residual over t in [0, 12] at N = 40, 60, 80.
Run from a checkout (cwd) of the commit under test. Records max |M_total| on cells 0-2 and on
the 1 <= r <= 5 inner region and the 5-10 annulus every 0.5 time units. dt = 0.06 dr."""
import sys, math, json, time
import numpy as np
sys.path.insert(0, ".")
from engine.production_kernel import V55ProductionKernel
from engine import v55_pirk_adapter as adapter
_, vac, _ = adapter.vendor_modules()
out = {}
for N in [int(x) for x in sys.argv[1].split(",")]:
    k = V55ProductionKernel(); s = k.initialize(resolution=N, r_max=40.0)
    dt = 0.06 * float(s.grid.dr); every = int(round(0.5 / dt)); rows = []; t0 = time.time()
    r = np.asarray(s.grid.centers)
    for n in range(int(round(12.0 / dt)) + 1):
        if n % every == 0:
            c = vac.constraints(s.grid, s.geometry); d_ = k.diagnostics(s, profiles=True); p = d_.get("profiles", d_)
            M = np.abs(np.asarray(c["momentum"]) - 8 * math.pi * np.asarray(p["total_j"]))
            rows.append({"t": float(s.t), "c02": float(M[:3].max()), "inner_1_5": float(M[(r >= 1) & (r <= 5)].max()),
                         "ann_5_10": float(M[(r >= 5) & (r <= 10)].max()), "argmax_r": float(r[int(np.argmax(M[r <= 10]))])})
        s = k.step(s, dt)
    out[N] = rows
    last = rows[-1]
    print(f"N={N}: t={last['t']:.2f} c0-2 {last['c02']:.2e} inner {last['inner_1_5']:.2e} ann {last['ann_5_10']:.2e} argmax r {last['argmax_r']:.2f} ({time.time()-t0:.0f}s)", flush=True)
json.dump(out, open(sys.argv[2], "w"), indent=1)
