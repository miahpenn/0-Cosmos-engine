"""Diagnostic only (stored native.npz, no simulation, no gate change): direction of the W1 input tails and whether
G0.5'a-d (AMENDMENT 2, Drive 1bpHE2NiBAoYmbMiOOakOqYhBrFARJLd1PsaP6baWmb8) catch a reversed (inward) tail.
Signed flux F(r) = -4 pi r^2 f_t f_r with f_t = P - H0 r f_r (chart term removed); positive = outward.
Usage: python3 research/w1_c9_direction/run.py native.npz OUTDIR
"""
import json
import math
import os
import sys
from pathlib import Path

sys.path.insert(0, os.getcwd())
import numpy as np  # noqa: E402

from engine import scalar_system as ss, v55_initial as vi, v55_pirk_adapter as ad  # noqa: E402
from research import b7_w2_origin_check_v1 as oc  # noqa: E402

z = np.load(sys.argv[1])
out = Path(sys.argv[2])
out.mkdir(parents=True, exist_ok=True)
go, vac, _ = ad.vendor_modules()
PHI0 = oc.clone.PHI0
OM2, K2, OM3, K3 = 3.818758, 3.543171, 1.5 * 3.818758, 5.363607
P_D_REF, A_D_PRED = 4.207e-10, 2.2245e-6
ENV = {0.025: (3.708e-8, 9.693e-8), 0.05: (3.709e-8, 9.696e-8)}
Ds = float(ss.sd_rest_value(PHI0))


def fit(r, a2, k):
    M = np.stack([np.ones_like(r), np.cos(2 * k * r), np.sin(2 * k * r)], 1)
    A, c, s = np.linalg.lstsq(M, a2, rcond=None)[0]
    return float(A), float(math.hypot(c, s) / A)


def g05(grid, dr, H0, S, PS, D, PD):  # G0.5'a-d exactly as research/w1_c9_g05p/run.py
    r = np.asarray(grid.centers)
    PDr = PD - H0 * r * grid.cell_derivative_fourth(D, parity=1)
    PSr = PS - H0 * r * grid.cell_derivative_fourth(S, parity=1)
    aD2 = r**2 * ((D - Ds) ** 2 + (PDr / OM2) ** 2)
    aS2 = r**2 * (S**2 + (PSr / OM3) ** 2)
    wD, wS = (r >= 20) & (r <= 60 - 5 * dr), (r >= 35) & (r <= 60 - 5 * dr)
    AD, mD = fit(r[wD], aD2[wD], K2)
    AS, mS = fit(r[wS], aS2[wS], K3)
    P = 2 * math.pi * OM2 * K2 * AD
    cD, cS = float(np.max(np.abs(aD2[-5:] / AD - 1))), float(np.max(np.abs(aS2[-5:] / AS - 1)))
    sc = math.sqrt(AD) / A_D_PRED * 1.03
    eD, eP = float(np.max(np.abs(D[-5:] - Ds))), float(np.max(np.abs(PD[-5:])))
    return {"a": abs(P / P_D_REF - 1) <= 0.03, "b": mD <= 0.01 and mS <= 0.01, "c": cD <= 0.03 + mD and cS <= 0.03 + mS,
            "d": eD <= ENV[dr][0] * sc and eP <= ENV[dr][1] * sc,
            "values": {"P_rel": P / P_D_REF - 1, "mD": mD, "mS": mS, "cD": cD, "cS": cS,
                       "D_last5/bound": eD / (ENV[dr][0] * sc), "PD_last5/bound": eP / (ENV[dr][1] * sc)}}


rep = {}
for dr, n in ((0.025, 2400), (0.05, 1200)):
    key = "dr%g" % dr
    grid = go.SphericalCellGrid(n, 60.0)
    r = np.asarray(grid.centers)
    H0 = float(vi.build_initial_data(grid, vac.flat_state, initial_structure="none", eight_node_areal=True).H0)
    S, PS, D, PD = (np.array(z[key + "_" + k], dtype=np.float64) for k in ("S", "PS", "D", "PD"))
    Dr, Sr = grid.cell_derivative_fourth(D, parity=1), grid.cell_derivative_fourth(S, parity=1)
    chD, chS = H0 * r * Dr, H0 * r * Sr
    wD, wS = (r >= 20) & (r <= 60 - 5 * dr), (r >= 35) & (r <= 60 - 5 * dr)
    FD = -4 * math.pi * r**2 * (PD - chD) * Dr
    FS = -4 * math.pi * r**2 * (PS - chS) * Sr
    flux = {"D": {"window_mean": float(np.mean(FD[wD])), "window_min": float(np.min(FD[wD])),
                  "share_positive": float(np.mean(FD[wD] > 0)), "window": [20.0, 60 - 5 * dr],
                  "P_D_ref": P_D_REF},
            "S": {"window_mean": float(np.mean(FS[wS])), "window_min": float(np.min(FS[wS])),
                  "share_positive": float(np.mean(FS[wS] > 0)), "window": [35.0, 60 - 5 * dr]}}
    # reversed tail: keep the field, reverse its rest momentum (f_t -> -f_t), keep the producer's chart term
    sw = 0.5 * (1 + np.tanh((r - 12.0) / 2))
    cases = {"baseline": (S, PS, D, PD),
             "S tail reversed": (S, PS - 2 * sw * (PS - chS), D, PD),
             "D tail reversed": (S, PS, D, PD - 2 * sw * (PD - chD))}
    rep[key] = {"H0": H0, "signed_flux": flux,
                "outward": {"D": flux["D"]["window_mean"] > 0, "S": flux["S"]["window_mean"] > 0},
                "faults": {nm: g05(grid, dr, H0, *a) for nm, a in cases.items()}}
    for nm in ("S tail reversed", "D tail reversed"):
        f = rep[key]["faults"][nm]
        f["caught_by"] = [g for g in "abcd" if not f[g]]
(out / "report.json").write_text(json.dumps(rep, indent=1, default=float))
for k, v in rep.items():
    print(k, "flux D mean %.4g (+%.3f) S mean %.4g (+%.3f)" % (v["signed_flux"]["D"]["window_mean"], v["signed_flux"]["D"]["share_positive"],
          v["signed_flux"]["S"]["window_mean"], v["signed_flux"]["S"]["share_positive"]), v["outward"])
    for nm, f in v["faults"].items():
        print("  ", nm, f.get("caught_by"), {x: round(y, 4) for x, y in f["values"].items()})
