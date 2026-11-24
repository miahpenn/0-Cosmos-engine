"""Compare the #595 rerun with Aster's material (README-34 5b206037, centre report fed74597, reproduction ZIP
on Drive): fine input, sampled trajectories, band definition, centre-trial readbacks. Stored data only, no step.
Usage (repo root): python research/w1_fring_a6/aster_compare.py ASTER_SNAPSHOT_DIR OUT.json"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline

snap, out = Path(sys.argv[1]), Path(sys.argv[2])
RUNS = Path("results/research/w1_fring_spatial_rerun/runs")
res = {}
a = np.load(snap/"results/research/w1_fring_spatial/static_0p0125.npz")
o = np.load("results/research/w1_fring_spatial_rerun/slice_dr0.0125.npz")
res["input_dr0.0125"] = {k: {"bitwise": a[k].tobytes() == o["dr0.0125_"+k].tobytes(),
                             "max_abs_diff": float(np.abs(a[k] - o["dr0.0125_"+k]).max()),
                             "max_abs_value": float(np.abs(a[k]).max())} for k in ("S", "PS", "D", "PD")}
res["input_dr0.0125"]["hb_x_max_abs_diff"] = float(np.abs(a["x"] - o["dr0.0125_x"]).max())
c = json.loads((snap/"results/research/w1_fring_spatial/comparison.json").read_text())
radii = np.array(c["radii"])
res["trajectories"] = {}
for dr in ("0.025", "0.0125"):
    for dt in ("0.000375", "0.0001875"):
        tag = f"A-dr{float(dr):g}-dt{float(dt):g}"
        ast = c["runs"][dr][dt]["rows"]
        ours = json.loads((RUNS/f"{tag}.json").read_text())["rows"]
        every = round(.00375/float(dt))
        z = np.load(RUNS/f"{tag}-split.npz")
        scale = ast[-1]["first_r"]/z["r"][0]
        R = z["r"]*scale
        m = R < .4
        res["trajectories"][tag] = {
            "samples": [{"t": s["t"], "aster_M_cell0": s["cell_0"]["M"], "ours_M_cell0": ours[i*every]["M_cell0"],
                         "aster_peak_M": s["peak"]["M"], "ours_peak_M": ours[i*every]["max_abs_M"],
                         "aster_peak_H": s["peak"]["H"], "ours_peak_H": ours[i*every]["max_abs_H"]}
                        for i, s in enumerate(ast)],
            "band_aster": ast[-1]["band_max"]["M"], "band_ours_cell_centres": ours[-1]["band_max_abs_M"],
            "band_ours_aster_definition": float(np.abs(CubicSpline(R[m], z["M"][m])(radii)).max()),
            "areal_scale_used": float(scale)}
t = res["trajectories"]
for dt in ("0.000375", "0.0001875"):
    f, g = t[f"A-dr0.0125-dt{float(dt):g}"], t[f"A-dr0.025-dt{float(dt):g}"]
    res[f"band_ratio_dt{dt}"] = {"aster": f["band_aster"]/g["band_aster"],
                                  "ours_cell_centres": f["band_ours_cell_centres"]/g["band_ours_cell_centres"],
                                  "ours_aster_definition": f["band_ours_aster_definition"]/g["band_ours_aster_definition"]}
s = json.loads((snap/"results/research/w1_fring_centre_projection/summary.json").read_text())
res["centre_summary_keys"] = list(s.keys())
out.write_text(json.dumps(res, indent=1))
print(json.dumps({k: v for k, v in res.items() if k.startswith("band_ratio")}, indent=1))
