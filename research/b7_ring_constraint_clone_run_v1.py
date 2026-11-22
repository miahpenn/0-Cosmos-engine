"""RESEARCH ONLY. B7 W2 driver: positive and negative controls, then the placeholder-chart ring test (no time evolution).

Usage: python research/b7_ring_constraint_clone_run_v1.py RING.npz OUTDIR     (RING.npz = source C t300 member, hash recorded)
Design and thresholds: Project doc claude/b7-w2-constraint-clone-design.md (fixed before the build). Engine used read-only.
"""
import json
import math
import os
import sys

sys.path.insert(0, ".")
import numpy as np  # noqa: E402

from engine import v55_pirk_adapter as adapter  # noqa: E402
from engine.production_kernel import V55ProductionKernel, build_initial_data  # noqa: E402
from engine.stress_energy import assemble_total_stress_energy  # noqa: E402
from engine import scalar_system as ss  # noqa: E402
from research import b7_ring_constraint_clone_v1 as clone  # noqa: E402
from research import b7_ring_constraint_v1 as con  # noqa: E402
from research import b7_ring_source_v1 as src  # noqa: E402

PC_TOL = 1e-12            # fixed before the build (design Section 3)
G4C_FACTOR = 10.0         # ring residual <= 10 x PC2 residual (design Section 4)
R_MAX = 60.0
GEOM_FIELDS = ("a", "b", "X", "Aa", "K", "Lambda")
LAM = ss.SD_QUARTIC_LAMBDA


def rel(x, y):
    return float(np.max(np.abs(x - y)) / max(float(np.max(np.abs(y))), 1e-300))


def residuals(kernel, grid, init):
    st = clone.research_state(kernel, grid, init)
    d = kernel.diagnostics(st)
    tot = assemble_total_stress_energy(grid, st.geometry, st.scalars, st.matter)
    al = np.asarray(st.geometry.alpha)
    return {"H_full_absmax": d["hamiltonian_full_absmax"], "M_full_absmax": d["momentum_full_absmax"],
            "H_inner_absmax": d["hamiltonian_inner_absmax"], "M_inner_absmax": d["momentum_inner_absmax"],
            "scale_16pi_rho_max": float(16 * math.pi * np.max(np.abs(tot.rho))),
            "scale_8pi_j_max": float(8 * math.pi * np.max(np.abs(tot.j))),
            "alpha_min": float(al.min()), "alpha_max": float(al.max())}


def engine_slice(grid, vac, structure, **kw):
    return build_initial_data(grid, vac.flat_state, initial_structure=structure, **kw)


def pc(grid, vac, structure, sea_floor=True, **kw):
    eng = engine_slice(grid, vac, structure, **kw)
    s = eng.scalars
    cl, _ = clone.reconstruct_slice(grid, vac.flat_state, s.S, s.PS, s.D, s.PD, s.phi, s.Pi, sea_floor=sea_floor)
    diffs = {f: rel(getattr(cl.geometry, f), getattr(eng.geometry, f)) for f in GEOM_FIELDS}
    return eng, cl, {"max_rel_diff": max(diffs.values()), "per_field": diffs, "H0_equal": cl.H0 == eng.H0,
                     "result": "PASS" if max(diffs.values()) <= PC_TOL and cl.H0 == eng.H0 else "FAIL"}


def ring_fields(npz, grid):
    z = np.load(npz)
    r_t = z["r"]
    out = {}
    for k_toy, k in (("t300_s", "S"), ("t300_ps", "PS"), ("t300_d", "D"), ("t300_pd", "PD")):
        out[k] = np.interp(np.asarray(grid.centers), r_t, z[k_toy]) / math.sqrt(LAM)   # PLACEHOLDER chart, u = 1
    n = len(grid.centers)
    out["phi"] = np.full(n, clone.PHI0)
    out["Pi"] = np.full(n, clone.PI0)
    return out


def main(npz, outdir):
    os.makedirs(outdir, exist_ok=True)
    kernel = V55ProductionKernel()
    grid_ops, vac, _ = adapter.vendor_modules()
    rep = {"lambda_ruler": LAM, "ring_member_sha256": src.sha256_file(npz), "runtime": src.environment_report(),
           "r_max": R_MAX, "chart": "PLACEHOLDER S=s/sqrt(lambda), D=d/sqrt(lambda), u=1 (not the engine ring; W1 pending); toy field truncated at r_max",
           "pc": {}, "grids": {}, "nc": {}}
    for dr in (0.2, 0.1, 0.05):
        n = int(round(R_MAX / dr))
        grid = grid_ops.SphericalCellGrid(n, R_MAX)
        g = {"n": n, "dr": dr}
        # positive controls (engine build_initial_data = initializer call, no evolution)
        pcs = {}
        e1, c1, pcs["PC1_none"] = pc(grid, vac, "none")
        e2, c2, pcs["PC2_rescaled_pocket"] = pc(grid, vac, "rescaled_pocket", pocket_x=0.5)
        e3, c3, pcs["PC3_hand_floor_reference"] = pc(grid, vac, "hand_floor_reference", sea_floor=False)
        g["pc"] = pcs
        g["pc_residuals"] = {"PC2_engine": residuals(kernel, grid, e2)}
        s2 = e2.scalars
        pz, _ = clone.reconstruct_slice(grid, vac.flat_state, s2.S, s2.PS, s2.D, s2.PD, s2.phi, s2.Pi, Kt_mode="zero")
        g["pc_residuals"]["PC2_clone_Kt_zero"] = residuals(kernel, grid, pz)
        # ring
        f = ring_fields(npz, grid)
        ring = {}
        for mode in ("hubble", "zero"):
            init, aux = clone.reconstruct_slice(grid, vac.flat_state, f["S"], f["PS"], f["D"], f["PD"], f["phi"], f["Pi"], Kt_mode=mode)
            ring[mode] = residuals(kernel, grid, init)
            ref = g["pc_residuals"]["PC2_engine" if mode == "hubble" else "PC2_clone_Kt_zero"]
            ring[mode]["G4c_vs_PC2"] = {k: ring[mode][k] <= G4C_FACTOR * ref[k] for k in ("H_full_absmax", "M_full_absmax", "H_inner_absmax", "M_inner_absmax")}
            ring[mode]["PC2_ref"] = {k: ref[k] for k in ("H_full_absmax", "M_full_absmax", "H_inner_absmax", "M_inner_absmax")}
        Sp = np.gradient(f["S"], grid.dr)
        Dp = np.gradient(f["D"], grid.dr)
        ring["j_parts_absmax"] = {"PS_Sr": float(np.max(np.abs(f["PS"] * Sp))), "PD_Dr": float(np.max(np.abs(f["PD"] * Dp)))}
        g["ring"] = ring
        # negative controls on the ring (hubble)
        full, _ = clone.reconstruct_slice(grid, vac.flat_state, f["S"], f["PS"], f["D"], f["PD"], f["phi"], f["Pi"])
        nops, _ = clone.reconstruct_slice(grid, vac.flat_state, f["S"], f["PS"], f["D"], f["PD"], f["phi"], f["Pi"], include_ps_momentum=False)
        nokin, _ = clone.reconstruct_slice(grid, vac.flat_state, f["S"], f["PS"], f["D"], f["PD"], f["phi"], f["Pi"], include_ps_kinetic=False)
        nc = {"NC1_K_rel_diff": rel(nops.geometry.K, full.geometry.K), "NC1_residuals_without_PS_Sr": residuals(kernel, grid, nops),
              "NC2_metric_rel_diff": max(rel(getattr(nokin.geometry, q), getattr(full.geometry, q)) for q in ("a", "b", "X")),
              "NC2_residuals_without_PS_kinetic": residuals(kernel, grid, nokin)}
        flat = {"s": f["S"], "d": f["D"], "ps": f["PS"], "pd": f["PD"]}
        drr = float(grid.dr)
        g4a_false = con.g4_gr_initial_data({"s": f["S"], "d": f["D"], "ps": f["PS"], "pd": f["PD"]}, drr,
                                           j_supplied=con.momentum_density(flat, drr, include_ps=False))["checks"][0]["result"]
        g4a_true = con.g4_gr_initial_data(flat, drr)["checks"][0]["result"]
        nc["NC1_G4a_without"], nc["NC1_G4a_with"] = g4a_false, g4a_true
        zS = np.zeros_like(f["S"])
        a_, _ = clone.reconstruct_slice(grid, vac.flat_state, zS, f["PS"], f["D"], f["PD"], f["phi"], f["Pi"])
        b_, _ = clone.reconstruct_slice(grid, vac.flat_state, zS, f["PS"], f["D"], f["PD"], f["phi"], f["Pi"], include_ps_momentum=False)
        nc["NC3_K_rel_diff_S0"] = rel(a_.geometry.K, b_.geometry.K)
        g["nc"] = nc
        rep["grids"]["dr%g" % dr] = g
    # verdicts
    pcs_ok = all(v["result"] == "PASS" for g in rep["grids"].values() for v in g["pc"].values())
    nc_ok = all(g["nc"]["NC1_K_rel_diff"] > PC_TOL and g["nc"]["NC1_G4a_without"] == "FAIL" and g["nc"]["NC1_G4a_with"] == "PASS"
                and g["nc"]["NC2_metric_rel_diff"] > PC_TOL and g["nc"]["NC3_K_rel_diff_S0"] <= PC_TOL for g in rep["grids"].values())
    rep["verdict"] = {"PC1_PC3_all_grids": "PASS" if pcs_ok else "FAIL", "NC1_NC3_as_expected": nc_ok,
                      "G4c_hubble_all_measures_within_10x_PC2": all(all(g["ring"]["hubble"]["G4c_vs_PC2"].values()) for g in rep["grids"].values())}
    json.dump(rep, open(os.path.join(outdir, "report.json"), "w"), indent=1, default=float)
    return rep


if __name__ == "__main__":
    main(*sys.argv[1:3])
