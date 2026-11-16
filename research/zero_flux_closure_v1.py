"""RESEARCH ONLY. Task 2 (HANDOFF Nov 15 (3)): the native outer closure for expanding runs, one-stone scope (confirmed by
Jeremiah: "identical neighbours, in step: one-stone simplification"). Compares production's open edge with the lattice
wall (research/lattice_cell_v1.LatticeCellKernel: W1-W5; K-dot from the wall's compatibility, not from the stored lapse)
on production's own merged machine (#503: coupled floor, lambda matched today, background-only default).
Note: results/zero-flux-closure-2026-11-16/README.md.

Runs (short):
  B  background, production default start: open edge vs wall, both against the production cosmos lane.
  S  the same with a seeded lapse error (alpha = 1 + eps (1 - r^2/R^2) before the first solve): the open-edge mode
     (results/floor-swap-2026-11-15/KERNEL_DRIFT.md) against the wall, which takes no lapse from the past.
Usage: python research/zero_flux_closure_v1.py B|S N RMAX T [eps]
"""
import math
import sys

sys.path.insert(0, ".")
import numpy as np  # noqa: E402
from scipy.integrate import solve_ivp  # noqa: E402

from engine.cosmos import CosmosParams, construct_present_day, cosmos_rhs  # noqa: E402
from engine.production_kernel import V55ProductionKernel  # noqa: E402
from research.lattice_cell_v1 import LatticeCellKernel, as_cell  # noqa: E402


def lane():
    P = CosmosParams(); s0 = construct_present_day(P)
    y0 = [s0.a, s0.H, s0.phi, s0.pi_phi, s0.rho_dm, s0.rho_b, s0.rho_r, s0.rho_shear]
    return solve_ivp(lambda t, y: list(cosmos_rhs(tuple(y), P)), (0, 2000), y0, method="DOP853", rtol=1e-11, atol=1e-18, dense_output=True)


def run(kind, N, R, T, eps=0.0, closure="i", cfl=0.06):
    st = V55ProductionKernel().initialize(resolution=N, r_max=R)
    k = V55ProductionKernel() if kind == "open" else LatticeCellKernel(closure)
    if kind != "open":
        st = as_cell(st)
    r = np.asarray(st.grid.centers)
    st.geometry.alpha[:] = 1.0 + eps * (1 - (r / R) ** 2)
    st.geometry.alpha[:] = k._solve_lapse_for_state(st.grid, st.geometry, st.scalars, st.matter, None)[0]
    dt = cfl * float(st.grid.dr); n = 0; out = []
    every = max(1, int(round(10.0 / dt)))
    while st.t < T:
        if n % every == 0:
            m = r < 0.5 * R
            out.append((st.t, -float(np.mean(st.geometry.K[m])) / 3, float(np.mean(st.scalars.phi[m])), float(st.geometry.alpha[0]),
                        float(np.ptp(st.geometry.alpha))))
        st = k.step(st, dt); n += 1
    return np.array(out)


if __name__ == "__main__" and sys.argv[1] in ("B", "S"):
    mode, N, R, T = sys.argv[1], int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
    eps = float(sys.argv[5]) if len(sys.argv) > 5 else (0.01 if mode == "S" else 0.0)
    L = lane()
    print(f"mode {mode}: N {N}, r_max {R}, to t = {T}, lapse seed eps = {eps}; production defaults (#503); lane = engine.cosmos")
    res = {kind: run(kind, N, R, T, eps) for kind in ("open", "wall")}
    print("   t       open edge: dH/H     dphi        alpha(0)-1   alpha spread |  wall: dH/H        dphi        alpha(0)-1   alpha spread")
    for i in range(len(res["open"])):
        o, w = res["open"][i], res["wall"][i]; y = L.sol(o[0])
        print(f"  {o[0]:6.1f}   {(o[1]-y[1])/y[1]:+.2e}   {o[2]-y[2]:+.2e}   {o[3]-1:+.2e}   {o[4]:.2e}   |  {(w[1]-y[1])/y[1]:+.2e}   {w[2]-y[2]:+.2e}   {w[3]-1:+.2e}   {w[4]:.2e}")


def pocket_compare(wall_path, open_path, r_core=12.0):
    """P  pocket (one-budget lambda = 5000, x_c 2.0, r_max 24): the wall (lattice_cell_run_v1) against the open edge
    (dm_one_budget_v1 main, ALPHA_SEED=one). Same physics; slices differ (conformally flat vs areal), so compare
    physical totals: S-D energy above the vacuum inside r_core and outside, its total, and H in the core."""
    from research.reset_handoff_v1 import sd_energy
    out = {}
    for lab, p in (("wall", wall_path), ("open", open_path)):
        d = np.load(p, allow_pickle=True); t = d["t"]; r = d["r"]; rows = []
        for i in range(len(t)):
            e = sd_energy(d, i); dV = d["dV"][i]; m = r < r_core
            rows.append((t[i], float(np.sum(e * dV * m)), float(np.sum(e * dV * ~m)), -float(np.mean(d["K"][i][m])) / 3))
        out[lab] = np.array(rows)
    return out


if __name__ == "__main__" and sys.argv[1] == "P":
    o = pocket_compare(sys.argv[2], sys.argv[3])
    w, op = o["wall"], o["open"]
    print("pocket, x_c 2.0, r_max 24, lambda 5000 (one budget): S-D energy (scaled, x proper volume per sr) and core H")
    print("    t     wall: E(r<12)   E(r>12)    total      H(core)   |  open: E(r<12)   E(r>12)    total      H(core)   | total open/wall")
    for tq in (0, 2, 4, 6, 8, 10, 12, 15, 20, 25, 30, 35, 40):
        i = int(np.argmin(np.abs(w[:, 0] - tq))); j = int(np.argmin(np.abs(op[:, 0] - tq)))
        if abs(op[j, 0] - tq) > 0.3 or abs(w[i, 0] - tq) > 0.3:
            continue
        print(f"  {tq:4.0f}   {w[i,1]:.4e}   {w[i,2]:.4e}   {w[i,1]+w[i,2]:.4e}   {w[i,3]:.5f}   |  {op[j,1]:.4e}   {op[j,2]:.4e}   {op[j,1]+op[j,2]:.4e}   {op[j,3]:.5f}   | {(op[j,1]+op[j,2])/(w[i,1]+w[i,2]):.4f}")


def same_start(kind, N=120, R=24.0, T=30.0, cfl=0.06, out=None):
    """Q  the same initial state (the cell's conformally flat one-budget pocket, x_c 2.0, lambda 5000) evolved with the
    wall and with production's open edge, identical dt. Records the S-D energy inside 12 and outside, and core H."""
    import dataclasses as _dc
    from engine.v55_pirk_adapter import vendor_modules
    from research.lattice_cell_v1 import cell_initial
    from research.reset_handoff_v1 import sd_energy
    k_cell, st, r, H0 = cell_initial("full", 2.0, 0.0, N, R)
    if kind == "open":
        grid_ops, _, _ = vendor_modules()
        st = _dc.replace(st, grid=grid_ops.SphericalCellGrid(N, R))
        k = V55ProductionKernel()
        st.geometry.alpha[:] = 1.0
        st.geometry.alpha[:] = k._solve_lapse_for_state(st.grid, st.geometry, st.scalars, st.matter, None)[0]
    else:
        k = k_cell
    dt = cfl * float(st.grid.dr); n = 0; rows = []
    every = max(1, int(round(0.5 / dt)))
    while st.t < T + 1e-9:
        if n % every == 0:
            from engine.v55_matter import metric_slice_from_q  # noqa: F401
            g = st.geometry; X = np.asarray(g.X)
            dV = np.sqrt(g.a) * g.b / X**3 * np.asarray(st.grid.volumes, float)
            d = dict(r=r, grr=[g.a / X**2], S=[st.scalars.S], D=[st.scalars.D], PS=[st.scalars.PS], PD=[st.scalars.PD])
            e = sd_energy({kk: (np.asarray(v) if kk == "r" else np.asarray(v)) for kk, v in d.items()}, 0)
            m = r < 12.0
            rows.append((st.t, float(np.sum(e * dV * m)), float(np.sum(e * dV * ~m)), -float(np.mean(g.K[m])) / 3, float(g.alpha[0]), float(g.alpha[-1])))
        try:
            st = k.step(st, dt); n += 1
        except Exception as exc:  # noqa: BLE001
            print(kind, "stopped at", st.t, exc); break
    np.save(out, np.array(rows))


if __name__ == "__main__" and sys.argv[1] == "Q":
    same_start(sys.argv[2], out=sys.argv[3])
