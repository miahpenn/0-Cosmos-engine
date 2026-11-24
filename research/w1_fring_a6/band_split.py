"""Band M split into Aa, K' and j parts on stored data (method fixed in
results/research/w1_fring_a6_band_split/PRERUN.md). No time step.
Usage (repo root): python -m research.w1_fring_a6.band_split OUT.json"""
import json
import sys
from pathlib import Path

import numpy as np

from engine import v55_pirk_adapter as ad
from engine.production_kernel import V55ProductionKernel
from engine.stress_energy import assemble_total_stress_energy
from research.w1_fring_a6.m_profile import load_state

ARMS = Path("results/research/w1_fring_a6/arms")
TRACE = Path("results/research/w1_fring_origin_trace")
KEYS = ("0.015", "0.8235", "1.6455", "3.291")
LO, HI = 0.05, 0.15


def explicit(grid, g, scal, mat, vac):
    d = lambda v, p: grid.cell_derivative_fourth(v, parity=p)
    r = np.asarray(grid.centers)
    parts = {"Aa_prime": d(g.Aa, 1), "three_Aa_over_r": 3*g.Aa/r,
             "Aa_b_prime_over_b": 1.5*g.Aa*d(g.b, 1)/g.b, "six_Aa_chi_prime": 6*g.Aa*(-.5*d(g.X, 1)/g.X),
             "K_part": -(2/3)*d(g.K, 1),
             "j_part": -8*np.pi*assemble_total_stress_energy(grid, g, scal, mat).j}
    parts["Aa_part"] = parts["Aa_prime"] + parts["three_Aa_over_r"] + parts["Aa_b_prime_over_b"] + parts["six_Aa_chi_prime"]
    M = vac.constraints(grid, g)["momentum"] + parts["j_part"]
    return r, M, parts


def stored_split(cells):
    Aa = np.array(cells["M_Aa"]) + np.array(cells["M_other_geometry"])
    K, j, M = np.array(cells["M_K"]), -np.array(cells["M_source"]), np.array(cells["M"])
    gap = float(np.max(np.abs(np.array(cells["M_geometry"]) - np.array(cells["M_source"]) - M)))
    return M, {"Aa_part": Aa, "K_part": K, "j_part": j}, gap


def band(r, dr, M, parts):
    w = (r >= LO) & (r <= HI)
    G = {k: float(np.sum(np.abs(v[w]))*dr) for k, v in parts.items() if k in ("Aa_part", "K_part", "j_part")}
    G["Aa_plus_K"] = float(np.sum(np.abs(parts["Aa_part"][w] + parts["K_part"][w]))*dr)
    G["M"] = float(np.sum(np.abs(M[w]))*dr)
    G["j_share"] = G["j_part"]/(G["Aa_part"] + G["K_part"] + G["j_part"])
    G["cells_r"] = r[w].tolist()
    return G


def main(out):
    res = {"band_r": [LO, HI], "grids": {}, "outer_edge": {}}
    # two grids at t 0.015 from stored sample cells
    srcs = {"dr0.05_dt1.5e-3 (A1 start)": (json.loads((ARMS/"A1.json").read_text())["saves"]["0.015"]["row"], .05),
            "dr0.05_dt7.5e-4 (A2 save)": (json.loads((ARMS/"A2.json").read_text())["saves"]["0.015"]["row"], .05),
            "dr0.025_dt7.5e-4": (json.loads((TRACE/"A-0.025.json").read_text())["rows"][-1]["cells"], .025),
            "dr0.025_dt3.75e-4": (json.loads((TRACE/"A-0.025-half.json").read_text())["rows"][-1]["cells"], .025)}
    for name, (cells, dr) in srcs.items():
        r = (np.array([0, 1, 2, 3, 4, 5]) + .5)*dr
        cells6 = {k: v[:6] for k, v in cells.items()}
        M, parts, gap = stored_split(cells6)
        res["grids"][name] = {"adds_back_gap": gap, **band(r, dr, M, parts)}
    c, f = res["grids"]["dr0.05_dt7.5e-4 (A2 save)"], res["grids"]["dr0.025_dt7.5e-4"]
    grow = {k: f[k]/c[k] for k in ("Aa_part", "K_part", "Aa_plus_K", "M")}
    res["matched_dt_growth"] = grow
    p1 = f["j_share"] < c["j_share"]
    p2 = (4 <= grow["Aa_part"] <= 16 and 4 <= grow["K_part"] <= 16 and grow["M"] > 1
          and (f["Aa_plus_K"] - c["Aa_plus_K"]) >= .5*(f["M"] - c["M"]))
    cl = res["grids"]["dr0.05_dt1.5e-3 (A1 start)"]
    res["matched_courant_growth"] = {k: f[k]/cl[k] for k in ("Aa_part", "K_part", "Aa_plus_K", "M")}
    res["matched_courant_j_share"] = [cl["j_share"], f["j_share"]]
    # outer edge, full grid at the A6 saves
    go, vac, _ = ad.vendor_modules()
    grid = go.SphericalCellGrid(1200, 60.)
    state = V55ProductionKernel().initialize(resolution=1200, r_max=60., outer_edge="open_edge",
                                             eight_node_areal=True, initial_structure="none")
    p3, addback_ok = True, True
    for arm in ("A1", "A2"):
        rec = json.loads((ARMS/f"{arm}.json").read_text())
        for key in KEYS:
            with np.load(ARMS/f"{arm}-t{key}.npz") as z:
                load_state(state, z, vac)
            r, M, parts = explicit(grid, state.geometry, state.scalars, state.matter, vac)
            total = parts["Aa_part"] + parts["K_part"] + parts["j_part"]
            rel = float(np.max(np.abs(total - M))/np.max(np.abs(M)))
            stored_gap = float(np.max(np.abs(M[[0, 1, 2, 3, 4, 5, -3, -2, -1]] - np.array(rec["saves"][key]["row"]["M"]))))
            terms = sum(np.abs(parts[k]) for k in ("Aa_prime", "three_Aa_over_r", "Aa_b_prime_over_b",
                                                   "six_Aa_chi_prime", "K_part", "j_part"))
            in_eps = float(np.max(np.abs(total - M)/(np.finfo(float).eps*terms)))
            addback_ok &= rel <= 1e-12
            w = np.flatnonzero(r >= 30)
            i = int(w[np.argmax(np.abs(M[w]))])
            mags = {k: abs(float(parts[k][i])) for k in ("Aa_part", "K_part", "j_part")}
            share = mags["K_part"]/sum(mags.values())
            p3 &= share >= .5
            res["outer_edge"][f"{arm} t{key}"] = {"cell": i, "r": float(r[i]), "M": float(M[i]),
                "parts": {k: float(parts[k][i]) for k in ("Aa_part", "K_part", "j_part")}, "K_share": share,
                "adds_back_rel": rel, "adds_back_in_eps_x_sum_terms": in_eps, "stored_cells_gap": stored_gap,
                "Aa_pieces": {k: float(parts[k][i]) for k in ("Aa_prime", "three_Aa_over_r", "Aa_b_prime_over_b", "six_Aa_chi_prime")}}
    res["full_grid_adds_back_1e-12_of_max_M"] = bool(addback_ok)
    res["prediction_1_j_share_shrinks"] = bool(p1)
    res["prediction_2_Aa_K_grow_8x_and_carry"] = bool(p2)
    res["prediction_3_K_carries_outer_edge"] = bool(p3)
    Path(out).write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main(sys.argv[1])
