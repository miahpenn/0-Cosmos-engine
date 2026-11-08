"""Legacy vs origin-accurate initializer quadrature over real production-kernel steps.

Run from a checkout whose engine/v55_initial.py has _cumulative_origin_integral (the
repaired initializer). The legacy rule (rectangle over [0, r0] plus plain trapezoids)
is swapped in by monkeypatching that helper; nothing else differs. Records the total
Hamiltonian, momentum and connection constraints at cells 0-2 and in 5 <= r <= 10, at
t = 0 and t = 0.60 (dt = 0.06 dr, r_max = 40), for N = 60, 120, 240.
Usage: python <this file> output.json
"""
import json, math, sys
import numpy as np
sys.path.insert(0, ".")
import engine.v55_initial as vi
from engine.production_kernel import V55ProductionKernel
from engine import v55_pirk_adapter as adapter

_, vac, _ = adapter.vendor_modules()
FIXED = vi._cumulative_origin_integral
def legacy(grid, r, dr, f):
    return np.cumsum(np.r_[0.5 * dr * f[0], 0.5 * dr * (f[1:] + f[:-1])])

def measure(k, s):
    c = vac.constraints(s.grid, s.geometry); d = k.diagnostics(s, profiles=True); p = d.get("profiles", d)
    H = np.asarray(c["hamiltonian"]) - 16 * math.pi * np.asarray(p["total_rho"])
    M = np.asarray(c["momentum"]) - 8 * math.pi * np.asarray(p["total_j"])
    C = np.asarray(c["connection"]); r = np.asarray(s.grid.centers); ann = (r >= 5) & (r <= 10)
    return {name: {"cells_0_2": float(np.max(np.abs(v[:3]))), "annulus_5_10": float(np.max(np.abs(v[ann])))}
            for name, v in (("H_total", H), ("M_total", M), ("C", C))}

def run(n, rule):
    vi._cumulative_origin_integral = rule
    try:
        k = V55ProductionKernel(); s = k.initialize(resolution=n, r_max=40.0)
    finally:
        vi._cumulative_origin_integral = FIXED
    out = {"t0": measure(k, s)}
    dt = 0.06 * float(s.grid.dr)
    for _ in range(int(round(0.6 / dt))):
        s = k.step(s, dt)
    out["t0.60"] = measure(k, s); out["t_end"] = float(s.t)
    return out

res = {}
for n in (60, 120, 240):
    res[n] = {"legacy": run(n, legacy), "fixed": run(n, FIXED)}
    for lab in ("legacy", "fixed"):
        x = res[n][lab]
        print(f"N={n:>3} {lab:6s} t=0: H {x['t0']['H_total']['cells_0_2']:.2e}/{x['t0']['H_total']['annulus_5_10']:.1e} | "
              f"t={x['t_end']:.2f}: " + " ".join(f"{q} {x['t0.60'][q]['cells_0_2']:.2e}/{x['t0.60'][q]['annulus_5_10']:.1e}" for q in ("M_total", "H_total", "C")))
json.dump(res, open(sys.argv[1], "w"), indent=1)
