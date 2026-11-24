"""Diagnostic reanalysis of #584's stored smoke data (no new steps). Usage: reanalyse.py results/research/w1_fring_smoke"""
import json
import sys
from pathlib import Path

import numpy as np

d = Path(sys.argv[1])
A = {a: json.loads((d / f"arm_{a}.json").read_text()) for a in ("A025", "B025", "A05", "B05")}
out = {"note": "Stored per-step data hold only the max-norm H and M and the areal R of the max; per-cell residuals "
               "(cells 1-5) were not stored, so they cannot be read without new steps.", "pairs": {}}
for arm in ("A", "B"):
    f, c = A[arm + "025"]["rows"], A[arm + "05"]["rows"]
    pr = {}
    for q in ("H", "M"):
        # common t: step 2k on dr 0.025 = step k on dr 0.05, k = 0..50
        t = np.array([f[2 * k]["t"] for k in range(51)])
        assert np.allclose(t, [c[k]["t"] for k in range(51)])
        rf = np.array([f[2 * k][q + "_max"] for k in range(51)])
        rc = np.array([c[k][q + "_max"] for k in range(51)])
        Rf = np.array([f[2 * k][q + "_argR"] for k in range(51)])
        Rc = np.array([c[k][q + "_argR"] for k in range(51)])
        orig = (Rf < 0.025) & (Rc < 0.05)
        sel = [1, 2, 5, 10, 25, 50]
        row = {"t": [float(t[k]) for k in sel], "dr025": [float(rf[k]) for k in sel], "dr05": [float(rc[k]) for k in sel],
               "ratio_05_over_025": [float(rc[k] / rf[k]) for k in sel],
               "argR_025": [float(Rf[k]) for k in sel], "argR_05": [float(Rc[k]) for k in sel],
               "both_at_origin_cell": [bool(orig[k]) for k in sel], "t0_values": [float(rf[0]), float(rc[0])]}
        # t-laws over the later half (t >= 0.0375) on each grid's full stored range
        for g, rows in (("dr025", f), ("dr05", c)):
            tt = np.array([x["t"] for x in rows])
            yy = np.array([x[q + "_max"] for x in rows])
            m = tt >= 0.0375
            m &= tt <= 0.075   # common range
            p = np.polyfit(np.log(tt[m]), np.log(yy[m]), 1)
            e = np.polyfit(tt[m], np.log(yy[m]), 1)
            rp = np.std(np.log(yy[m]) - np.polyval(p, np.log(tt[m])))
            re = np.std(np.log(yy[m]) - np.polyval(e, tt[m]))
            row["law_" + g] = {"power_exponent": float(p[0]), "power_resid": float(rp),
                               "exp_rate": float(e[0]), "exp_resid": float(re)}
        pr[q] = row
    out["pairs"][arm] = pr
(d / "reanalysis.json").write_text(json.dumps(out, indent=1))
for arm, pr in out["pairs"].items():
    for q, row in pr.items():
        print(arm, q, "ratio", ["%.2f" % x for x in row["ratio_05_over_025"]], "origin", row["both_at_origin_cell"])
        print("   argR025", row["argR_025"], "argR05", row["argR_05"])
        print("   025", ["%.2e" % x for x in row["dr025"]], "05", ["%.2e" % x for x in row["dr05"]])
        print("   laws", row["law_dr025"], row["law_dr05"])
