"""Compare the spatial replays with Aster's targets, split band M on two grids, and record what the (B) projection
touches beyond cell 0. Reads replay outputs only (plus one initial-data evaluation for the projection record).
Usage (repo root): python -m research.w1_fring_a6.spatial_compare REPLAY_DIR SLICE_NPZ OUT.json"""
import json
import sys
from pathlib import Path

import numpy as np

from engine import scalar_system as ss, v55_initial as vi, v55_pirk_adapter as ad
from engine.production_kernel import V55ProductionKernel
from research.w1_fring_a6 import spatial_replay as sr

TARGET_A = {(.025, .000375): 3.5574e-12, (.0125, .000375): 2.8602e-11,
            (.025, .0001875): 6.8067e-12, (.0125, .0001875): 8.2393e-11}
TERMS = ("Aa_prime", "three_Aa_over_r", "Aa_b_prime_over_b", "six_Aa_chi_prime", "K_part", "j_part")


def band_G(z, lo=.05, hi=.25):
    r = z["r"]
    dr = float(r[1] - r[0])
    w = (r >= lo) & (r <= hi)
    Aa = sum(z[k] for k in TERMS[:4])
    parts = {"Aa_part": Aa, "K_part": z["K_part"], "j_part": z["j_part"]}
    G = {k: float(np.sum(np.abs(v[w]))*dr) for k, v in parts.items()}
    G["Aa_plus_K"] = float(np.sum(np.abs(Aa[w] + z["K_part"][w]))*dr)
    G["M"] = float(np.sum(np.abs(z["M"][w]))*dr)
    G["j_share"] = G["j_part"]/(G["Aa_part"] + G["K_part"] + G["j_part"])
    o = np.flatnonzero(r >= 30)
    i = int(o[np.argmax(np.abs(z["M"][o]))])
    mags = {k: abs(float(v[i])) for k, v in parts.items()}
    G["outer"] = {"r": float(r[i]), "M": float(z["M"][i]), **{k: float(v[i]) for k, v in parts.items()},
                  "K_share": mags["K_part"]/sum(mags.values())}
    total = sum(z[k] for k in TERMS)
    G["adds_back_in_eps_x_sum_terms"] = float(np.max(np.abs(total - z["M"])/(np.finfo(float).eps*sum(np.abs(z[k]) for k in TERMS))))
    return G


def main(rdir, slice_npz, out):
    rdir = Path(rdir)
    res = {"runs": {}, "split": {}, "predictions": {}}
    rows = {}
    for f in sorted(rdir.glob("[AB]-dr*.json")):
        d = json.loads(f.read_text())
        k = ("B" if d["recon"] else "A", d["dr"], d["dt"])
        rows[k] = d["rows"][-1]
        t = TARGET_A.get(k[1:]) if k[0] == "A" else None
        res["runs"][f.stem] = {"peak_abs_M": rows[k]["max_abs_M"], "peak_cell": rows[k]["max_cell"],
                               "cell0_abs_M": abs(rows[k]["M_cell0"]), "band_max_abs_M": rows[k]["band_max_abs_M"],
                               "target_peak": t, "ratio_to_target": rows[k]["max_abs_M"]/t if t else None,
                               "patched": d["patched"]}
        res["split"][f.stem] = band_G(np.load(rdir/f"{f.stem}-split.npz"))
    for a in "AB":
        for dt in (.000375, .0001875):
            res["runs"][f"band_ratio_{a}_dt{dt:g}"] = rows[(a, .0125, dt)]["band_max_abs_M"]/rows[(a, .025, dt)]["band_max_abs_M"]
    res["runs"]["B_cell0_fall_dr0.0125_dt0.000375"] = abs(rows[("A", .0125, .000375)]["M_cell0"])/abs(rows[("B", .0125, .000375)]["M_cell0"])
    for a in "AB":
        for dt in ("0.000375", "0.0001875"):
            c, f_ = res["split"][f"{a}-dr0.025-dt{dt}"], res["split"][f"{a}-dr0.0125-dt{dt}"]
            g = {k: f_[k]/c[k] for k in ("Aa_part", "K_part", "Aa_plus_K", "M")}
            res["predictions"][f"{a} dt{dt}"] = {
                "growth": g, "j_share": [c["j_share"], f_["j_share"]],
                "1_j_share_shrinks": f_["j_share"] < c["j_share"],
                "2_Aa_K_grow_4_16_and_carry": bool(4 <= g["Aa_part"] <= 16 and 4 <= g["K_part"] <= 16 and g["M"] > 1
                                                   and f_["Aa_plus_K"] - c["Aa_plus_K"] >= .5*(f_["M"] - c["M"])),
                "3_K_carries_outer_edge": c["outer"]["K_share"] >= .5 and f_["outer"]["K_share"] >= .5}
    # what the (B) projection writes, on the dr 0.025 initial data
    go, vac, _ = ad.vendor_modules()
    grid = go.SphericalCellGrid(2400, 60.)
    with np.load(sr.STATIC) as z:
        vals = [np.array(z["dr0.025_"+k], dtype=float) for k in ("S", "PS", "D", "PD")]
    vals.extend((np.full(2400, .179055), np.full(2400, .00341913)))
    st = V55ProductionKernel().initialize(resolution=2400, r_max=60., outer_edge="open_edge", eight_node_areal=True,
        initial_structure="supplied_slice", supplied_slice=vi.SuppliedSlice(ss.ScalarFields(*vals),
        np.asarray(grid.centers).copy(), sr.H0_CHART, sr.EDGE, "projection record, dr 0.025 corrected slice"))
    g0 = st.geometry.copy()
    full = vac.enforce_algebraic_regularity(grid, g0)
    ours = g0.copy()
    sr.project_first_cell(vac, grid, ours)
    rec = {}
    for name, g in (("ours_first_cell_only", ours), ("vendor_function_whole", full)):
        rec[name] = {k: {"cell0": float(getattr(g, k)[0] - getattr(g0, k)[0]),
                         "beyond_cell0_max_abs": float(np.max(np.abs(getattr(g, k)[1:] - getattr(g0, k)[1:])))}
                     for k in ("a", "b", "Aa", "beta", "Lambda", "B")}
        rec[name]["a_b2_minus_1_beyond_cell0_max_abs"] = float(np.max(np.abs(g.a[1:]*g.b[1:]**2 - 1)))
    rec["unmodified_a_b2_minus_1_max_abs"] = float(np.max(np.abs(g0.a*g0.b**2 - 1)))
    rec["unmodified_B_minus_0.75_Lambda_max_abs"] = float(np.max(np.abs(g0.B - .75*g0.Lambda)))
    res["projection_record"] = rec
    Path(out).write_text(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    main(*sys.argv[1:4])
