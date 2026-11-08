"""Old vs new CMC centre row over real production-kernel steps (read-only comparison)."""
import importlib.util, json, math, sys, time
import numpy as np
sys.path.insert(0, ".")
import engine.production_kernel as pk
import engine.cmc_gauge as new
from engine import v55_pirk_adapter as adapter

spec = importlib.util.spec_from_file_location("engine.cmc_gauge_old", sys.argv[1])
old = importlib.util.module_from_spec(spec); old.__package__ = "engine"; spec.loader.exec_module(old)
T_END, CFL, R_MAX = 0.60, 0.06, 40.0

def constraints(k, s):
    _, vacuum, _ = adapter.vendor_modules()
    c = vacuum.constraints(s.grid, s.geometry)
    d = k.diagnostics(s, profiles=True)
    prof = d.get("profiles", d)
    rho = np.asarray(prof["total_rho"], float); j = np.asarray(prof["total_j"], float)
    return np.asarray(c["hamiltonian"]) - 16 * math.pi * rho, np.asarray(c["momentum"]) - 8 * math.pi * j, np.asarray(c["connection"])

def run(n, solver):
    pk.solve_cmc_lapse = solver
    k = pk.V55ProductionKernel()
    s = k.initialize(resolution=n, r_max=R_MAX)
    dt = CFL * float(s.grid.dr)
    steps = int(round(T_END / dt))
    hist = []
    for i in range(steps):
        s = k.step(s, dt)
        if (i + 1) % max(1, steps // 5) == 0 or i + 1 == steps:
            H, M, C = constraints(k, s)
            r = np.asarray(s.grid.centers)
            ann = (r >= 5) & (r <= 10)
            hist.append({"t": float(s.t), "M_cells_0_2": float(np.max(np.abs(M[:3]))), "M_annulus": float(np.max(np.abs(M[ann]))),
                         "H_cells_0_2": float(np.max(np.abs(H[:3]))), "H_annulus": float(np.max(np.abs(H[ann]))),
                         "C_cells_0_2": float(np.max(np.abs(C[:3]))), "C_annulus": float(np.max(np.abs(C[ann])))})
    return s, hist

out = {}
for n in (60, 120, 240):
    t0 = time.time()
    s_old, h_old = run(n, old.solve_cmc_lapse)
    s_new, h_new = run(n, new.solve_cmc_lapse)
    r = np.asarray(s_new.grid.centers)
    diff = {}
    for name in ("alpha", "a", "b", "X", "Aa", "K", "Lambda"):
        dd = np.abs(np.asarray(getattr(s_new.geometry, name)) - np.asarray(getattr(s_old.geometry, name)))
        diff[name] = {"cells_0_2": float(dd[:3].max()), "r_ge_5": float(dd[r >= 5].max())}
    al = np.asarray(s_new.geometry.alpha); ao = np.asarray(s_old.geometry.alpha)
    out[n] = {"old": h_old, "new": h_new, "final_state_new_minus_old": diff,
              "alpha1_minus_alpha0": {"old": float(ao[1] - ao[0]), "new": float(al[1] - al[0])}, "seconds": time.time() - t0}
    print(f"N={n} ({out[n]['seconds']:.0f}s) t={h_new[-1]['t']:.3f}")
    for lab, h in (("old", h_old), ("new", h_new)):
        print(f"   {lab}: " + " | ".join(f"t={x['t']:.2f} M0-2={x['M_cells_0_2']:.2e} Mann={x['M_annulus']:.2e} H0-2={x['H_cells_0_2']:.2e} C0-2={x['C_cells_0_2']:.2e} Cann={x['C_annulus']:.1e}" for x in h[-2:]))
    print("   final new-old:", {k: (f"{v['cells_0_2']:.1e}", f"{v['r_ge_5']:.1e}") for k, v in diff.items()}, "a1-a0:", out[n]["alpha1_minus_alpha0"])
json.dump(out, open(sys.argv[2], "w"), indent=1)
