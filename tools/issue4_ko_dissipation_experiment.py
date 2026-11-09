"""Issue #4 experiment (research only): does vendor KO dissipation on the geometry remove the
growing centre time-step error? Monkeypatches the kernel's geometry stage terms; production untouched.
Usage: ko_exp.py EPS N dt1,dt2 out.json"""
import sys, math, json
import numpy as np
sys.path.insert(0, ".")
from engine.production_kernel import V55ProductionKernel
from engine import v55_pirk_adapter as adapter
_, vac, _ = adapter.vendor_modules()
ko = vac.ko_dissipation
EPS = float(sys.argv[1]); N = int(sys.argv[2]); dts = [float(x) for x in sys.argv[3].split(",")]
orig = V55ProductionKernel._geometry_stage_terms

def patched(self, grid, g, s, m, em):
    t = orig(self, grid, g, s, m, em)
    if EPS == 0.0:
        return t
    t = dict(t); t["explicit"] = dict(t["explicit"]); t["primary_l3"] = dict(t["primary_l3"])
    for name in ("a", "b", "X"):
        t["explicit"][name] = t["explicit"][name] + ko(grid, getattr(g, name), 1, epsilon=EPS)
    t["primary_l3"]["Aa"] = t["primary_l3"]["Aa"] + ko(grid, g.Aa, 1, epsilon=EPS)
    t["primary_l3"]["K"] = t["primary_l3"]["K"] + ko(grid, g.K, 1, epsilon=EPS)
    t["lambda_l3"] = t["lambda_l3"] + ko(grid, g.Lambda, -1, epsilon=EPS)
    return t
V55ProductionKernel._geometry_stage_terms = patched

res = []
for dt in dts:
    k = V55ProductionKernel(); s = k.initialize(resolution=N, r_max=40.0)
    for _ in range(int(round(0.6 / dt))):
        s = k.step(s, dt)
    c = vac.constraints(s.grid, s.geometry); d_ = k.diagnostics(s, profiles=True); p = d_.get("profiles", d_)
    r = np.asarray(s.grid.centers); ann = (r >= 5) & (r <= 10)
    M = np.asarray(c["momentum"]) - 8 * math.pi * np.asarray(p["total_j"])
    H = np.asarray(c["hamiltonian"]) - 16 * math.pi * np.asarray(p["total_rho"])
    row = {"N": N, "dt": dt, "eps": EPS, "M0": float(M[0]), "M_c02": float(np.max(abs(M[:3]))),
           "M_ann": float(np.max(abs(M[ann]))), "H_c02": float(np.max(abs(H[:3]))), "H_ann": float(np.max(abs(H[ann]))),
           "Aa0": float(s.geometry.Aa[0])}
    res.append(row)
    print(f"eps={EPS} N={N} dt={dt}: M cell0 {row['M0']:+.3e}  M c0-2 {row['M_c02']:.2e}  M ann {row['M_ann']:.2e}  H c0-2 {row['H_c02']:.2e}  H ann {row['H_ann']:.2e}", flush=True)
json.dump(res, open(sys.argv[4], "w"), indent=1)
