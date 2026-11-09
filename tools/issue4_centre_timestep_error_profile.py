"""Issue #4: where does the centre time-step error live? N fixed, dt and dt/2 to t=0.6; profile of
|field(dt) - field(dt/2)| at cells 0..11 and the annulus max, for every evolved field and the lapse."""
import sys, math, json
import numpy as np
sys.path.insert(0, ".")
from engine.production_kernel import V55ProductionKernel
from engine import v55_pirk_adapter as adapter
_, vac, _ = adapter.vendor_modules()
N = int(sys.argv[1]); dts = [float(x) for x in sys.argv[2].split(",")]
runs = []
for dt in dts:
    k = V55ProductionKernel(); s = k.initialize(resolution=N, r_max=40.0)
    for _ in range(int(round(0.6 / dt))):
        s = k.step(s, dt)
    c = vac.constraints(s.grid, s.geometry); d_ = k.diagnostics(s, profiles=True); p = d_.get("profiles", d_)
    M = np.asarray(c["momentum"]) - 8 * math.pi * np.asarray(p["total_j"])
    f = {n: np.asarray(getattr(s.geometry, n), float) for n in ("a", "b", "X", "K", "Aa", "Lambda", "alpha")}
    f.update({n: np.asarray(getattr(s.scalars, n), float) for n in ("S", "PS", "D", "PD", "phi", "Pi")})
    f["M"] = M
    runs.append(f); r = np.asarray(s.grid.centers)
ann = (r >= 5) & (r <= 10)
out = {}
for n in runs[0]:
    d = runs[0][n] - runs[1][n]
    out[n] = {"cells_0_11": [float(x) for x in d[:12]], "ann_max": float(np.max(np.abs(d[ann]))),
              "c02_max": float(np.max(np.abs(d[:3])))}
    print(f"{n:7s} c0-2 {out[n]['c02_max']:.2e} ann {out[n]['ann_max']:.2e} | " + " ".join(f"{x:+.1e}" for x in d[:10]), flush=True)
json.dump(out, open(sys.argv[3], "w"), indent=1)
