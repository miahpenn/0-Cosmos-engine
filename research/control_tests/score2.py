"""Stage-2 (and stage-3) Test 2 gates, scored on stored outputs (python3 with mpmath; no engine import, no run).

  python3 -m research.control_tests.score2 STAGE1_DIR STAGE2_DIR [CDOT_FLOOR.json]

STAGE1_DIR: results/research/control_tests_stage1 (arms/ holds the 1a floor arms; readback.json holds the 1a GK.1
floors). STAGE2_DIR: where test2 wrote its arms and mp_readback wrote readback.json. Writes STAGE2_DIR/score2.json.

Gates as predeclared, none changed (design 3.3; A1 sections 2-3; A2 P9-P11 and section 3; A3 sections 2, 3, 6, 7;
A4-A5 GK.1; A5 arm list; A20 and A21 as adopted on Jeremiah's "Yes to both", 08:07:50 UTC):
  G2.2  (H - floor) and (M - floor) ratio T2-3 / T2-6 (amp 1e-2 / 5e-3) in [3.6, 4.4].
  G2.3  M(T2-3) / (M(T2-4) / 9) >= 4 (T2-4 at amp 3e-2, A21). A point-instance is two-sided while its .025 floor is
        <= 1.44e-9, otherwise one-sided (can register FAIL only); PASS needs >= 10 two-sided counted over >= 3 saves.
  G2.4  Q = (4M/r + d_tH) / |d_tH| at the 4 band points, saves t >= 1.2, |Q| <= 0.5 (the #602 rule), T2-3 with
        the T2-6 repeat. d_tH: centred over the steps next to the save (one-sided second order at the last save).
  Edge  T2-3 (r_max 8) against T2-E (12): |Y_8 - Y_12| <= 3 |d alpha_c| |Y_8| + 10 floor, Y = H, M (A3 section 6).
  GK.1  |K residual - prediction| <= 10 x the 1a floor at the same point, save and stage (A4, A5), T2-3 and T2-4.
  GS.2  the semi-discrete M rate (Cdot at the saves): amp^2 ratio T2-3 / T2-6 in [3.6, 4.4]; dr ratio
        T2-3 / (T2-4 / 9) >= 4. Floors from test2 cdot_floor (a reanalysis of the stage-1 saves).
  GM.1  60-digit rho and j from the stored S, P_S: the largest band shift of 8 pi j is <= 1e-6 x the band M rate
        x dt, at every save (T2-3).
  GS.1, GS.3 (stage 3, only once T2-5 exists): Simpson split of A2 P10 on the per-step Cdot series.
  G2.1 (stage 3) is not implemented here: it needs the exact Klein-Gordon reference (design 3.2).
Aggregation (A3 section 3): a point-instance counts only where every compared signal is >= 10 x its floor; PASS
needs >= 10 counted over >= 3 saves, all passing; FAIL needs >= 10 counted with >= 25% failing; else UNDETERMINED.
Saves t 0.3 to 2.1 (t = 0 is not a gate save). Floors: the 1a arm with the same dr, dt and r_max (A5 g), signed.
"""
import json
import math
import sys
from pathlib import Path

import mpmath as mp
import numpy as np

from research.control_tests import mp60 as M
from research.control_tests.mp_readback import fixed_points

R_PTS = [None, 0.075, 0.125, 0.175, 0.225]
BAND = (1, 2, 3, 4)
ALL = (0, 1, 2, 3, 4)
TWO_SIDED_FLOOR = 1.44e-9          # A21
SCALE_T24 = 9.0                    # A21: (3e-2 / 1e-2)^2


def floor_arm(dr, dt, r_max):
    tag = {0.05: "05", 0.025: "025"}[round(dr, 6)]
    dtag = {3.75e-4: "dt0", 1.875e-4: "dt1"}[round(dt, 10)]
    return f"1a_dr{tag}_{dtag}_r{r_max:g}"


def verdict(insts, two_sided_only_for_pass=False):
    c = [x for x in insts if x["counted"]]
    fails = [x for x in c if not x["pass"]]
    base = [x for x in c if x.get("two_sided", True)] if two_sided_only_for_pass else c
    saves = {x["save"] for x in base}
    if len(base) >= 10 and len(saves) >= 3 and not fails:
        v = "PASS"
    elif len(c) >= 10 and len(fails) >= 0.25 * len(c):
        v = "FAIL"
    else:
        v = "UNDETERMINED"
    out = {"verdict": v, "counted": len(c), "failing": len(fails), "saves_counted": len(saves),
           "failing_instances": fails, "instances": insts}
    if two_sided_only_for_pass:
        out["two_sided_counted"] = len(base)
        out["one_sided_counted"] = len(c) - len(base)
    return out


class Run:
    def __init__(self, d, name):
        self.name = name
        self.rec = json.loads((d / f"{name}.json").read_text())
        self.saves = self.rec["saves"]
        self.dt, self.dr = self.rec["dt"], self.rec["dr"]
        self.near = {r["step"]: r for r in self.rec.get("near_save", [])}

    def y(self, key, k, p):
        return float(self.saves[k][key][p])

    def dtH(self, k, p):
        n = self.saves[k]["step"]
        if n + 1 in self.near and n - 1 in self.near:
            return (self.near[n + 1]["H"][p] - self.near[n - 1]["H"][p]) / (2 * self.dt)
        if n - 1 in self.near and n - 2 in self.near:
            return (3 * self.saves[k]["H"][p] - 4 * self.near[n - 1]["H"][p] + self.near[n - 2]["H"][p]) / (2 * self.dt)
        return None


def g22(r3, r6, fl):
    out = {}
    for key in ("H", "M"):
        ins = []
        for k in range(1, len(r3.saves)):
            for p in ALL:
                F = fl.y(key, k, p)
                y3, y6 = r3.y(key, k, p), r6.y(key, k, p)
                ratio = (y3 - F) / (y6 - F) if y6 != F else float("nan")
                ins.append({"save": k, "t": r3.saves[k]["t"], "point": p, "Y_T2-3": y3, "Y_T2-6": y6, "floor": F,
                            "ratio": ratio, "counted": abs(y3) >= 10 * abs(F) and abs(y6) >= 10 * abs(F),
                            "pass": 3.6 <= ratio <= 4.4})
        out[key] = verdict(ins)
    vs = {out[k]["verdict"] for k in ("H", "M")}
    out["verdict"] = "FAIL" if "FAIL" in vs else ("PASS" if vs == {"PASS"} else "UNDETERMINED")
    return out


def g23(r3, r4, f05, f025):
    ins = []
    for k in range(1, len(r3.saves)):
        for p in ALL:
            F5, F25 = f05.y("M", k, p), f025.y("M", k, p)
            y3, y4 = r3.y("M", k, p), r4.y("M", k, p)
            ratio = y3 / (y4 / SCALE_T24) if y4 != 0 else float("nan")
            ins.append({"save": k, "t": r3.saves[k]["t"], "point": p, "M_T2-3": y3, "M_T2-4": y4,
                        "floor_05": F5, "floor_025": F25, "ratio": ratio,
                        "counted": abs(y3) >= 10 * abs(F5) and abs(y4) >= 10 * abs(F25),
                        "two_sided": abs(F25) <= TWO_SIDED_FLOOR, "pass": ratio >= 4})
    return verdict(ins, two_sided_only_for_pass=True)


def g24(run, fl):
    ins = []
    for k in range(len(run.saves)):
        t = run.saves[k]["t"]
        if t < 1.2 - 1e-9:
            continue
        for p in BAND:
            dh = run.dtH(k, p)
            m = run.y("M", k, p)
            F = fl.y("M", k, p)
            if dh is None or dh == 0:
                ins.append({"save": k, "t": t, "point": p, "counted": False, "pass": False, "note": "no d_tH"})
                continue
            Q = (4 * m / R_PTS[p] + dh) / abs(dh)
            ins.append({"save": k, "t": t, "point": p, "M": m, "d_tH": dh, "Q": Q, "floor_M": F,
                        "counted": abs(m) >= 10 * abs(F) and abs(R_PTS[p] * dh / 4) >= 10 * abs(F),
                        "pass": abs(Q) <= 0.5})
    out = verdict(ins)
    out["counting_note"] = ("counted where |M| and |r d_tH / 4| are both >= 10 x the 1a M floor (no 1a d_tH floor "
                            "exists): Wright's reading, for Lumen^2 to confirm")
    return out


def edge(r8, r12, f8, f12):
    out = {}
    for key in ("H", "M"):
        ins = []
        for k in range(1, len(r8.saves)):
            da = r8.saves[k]["alpha_rhs_cell0"] - r12.saves[k]["alpha_rhs_cell0"]
            for p in ALL:
                F = max(abs(f8.y(key, k, p)), abs(f12.y(key, k, p)))
                y8, y12 = r8.y(key, k, p), r12.y(key, k, p)
                tol = 3 * abs(da) * abs(y8) + 10 * F
                ins.append({"save": k, "t": r8.saves[k]["t"], "point": p, "Y_8": y8, "Y_12": y12, "d_alpha_c": da,
                            "floor": F, "tolerance": tol, "counted": abs(y8) >= 10 * F and abs(y12) >= 10 * F,
                            "pass": abs(y8 - y12) <= tol})
        out[key] = verdict(ins)
    vs = {out[k]["verdict"] for k in ("H", "M")}
    out["verdict"] = "FAIL" if "FAIL" in vs else ("PASS" if vs == {"PASS"} else "UNDETERMINED")
    return out


def gk1(steps, floor_steps):
    out = {}
    for stage in ("primary", "accepted"):
        ins = []
        for k in range(1, min(len(steps), len(floor_steps))):
            g, f = steps[k]["GK1"][stage], floor_steps[k]["GK1"][stage]
            for p in ALL:
                F = abs(f["residual"][p])
                ins.append({"save": k, "point": p, "pred": g["pred"][p], "residual": g["residual"][p], "floor": F,
                            "counted": abs(g["pred"][p]) >= 10 * F, "pass": abs(g["residual"][p]) <= 10 * F})
        out[stage] = verdict(ins)
    vs = {out[s]["verdict"] for s in ("primary", "accepted")}
    out["verdict"] = "FAIL" if "FAIL" in vs else ("PASS" if vs == {"PASS"} else "UNDETERMINED")
    return out


def gs2(r3, r6, r4, cfl):
    """Semi-discrete M rate at the saves. cfl: the 1a cdot floors (cdot_floor.json) or None."""
    def floor(run, k, p):
        if cfl is None:
            return None
        a = cfl.get(floor_arm(run.dr, run.dt, run.rec["r_max"]))
        return None if a is None else abs(a["saves"][k]["cM"][p])
    out = {}
    for label, other, lo, hi, scale in (("amp2", r6, 3.6, 4.4, 1.0), ("dr", r4, 4.0, math.inf, SCALE_T24)):
        ins = []
        for k in range(1, len(r3.saves)):
            for p in ALL:
                a, b = r3.saves[k]["cM"][p], other.saves[k]["cM"][p] / scale
                F3, Fo = floor(r3, k, p), floor(other, k, p)
                ratio = a / b if b != 0 else float("nan")
                cnt = F3 is not None and Fo is not None and abs(a) >= 10 * F3 and abs(other.saves[k]["cM"][p]) >= 10 * Fo
                ins.append({"save": k, "point": p, "cM_T2-3": a, f"cM_{other.name}_scaled": b, "ratio": ratio,
                            "floor_T2-3": F3, f"floor_{other.name}": Fo, "counted": cnt, "pass": lo <= ratio <= hi})
        out[label] = verdict(ins)
        if cfl is None:
            out[label]["note"] = "no cdot_floor file: nothing counted"
    vs = {out[k]["verdict"] for k in ("amp2", "dr")}
    out["verdict"] = "FAIL" if "FAIL" in vs else ("PASS" if vs == {"PASS"} else "UNDETERMINED")
    return out


def gm1(d, run):
    """60-digit scalar rho and j at each save of the arm, against the stored binary64 (A1 section 3)."""
    rows = []
    dr, dt = run.dr, run.dt
    with np.load(d / f"{run.name}.npz") as z:
        for k in range(len(run.saves)):
            n = len(z[f"t{k}_g_a"])
            g = M.Grid(n, dr)
            S, P = M.arr(z[f"t{k}_s_S"]), M.arr(z[f"t{k}_s_PS"])
            X, a = M.arr(z[f"t{k}_g_X"]), M.arr(z[f"t{k}_g_a"])
            fp = g.d(S, 1, 1)
            rho = [((0.5 * p_) * p_ + 0.5 * (((x * x) / a_) * f) * f) + 0.5 * (s_ * s_)
                   for p_, f, x, a_, s_ in zip(P, fp, X, a, S)]
            j = [-(p_ * f) for p_, f in zip(P, fp)]
            with np.load(d / f"{run.name}_save{k}.npz") as w:
                r64, j64 = np.asarray(w["m_rho_scalar"], dtype=float), np.asarray(w["m_j_scalar"], dtype=float)
            drho = np.array([float(mp.mpf(float(r64[i])) - rho[i].v) for i in range(n)])
            dj = np.array([float(mp.mpf(float(j64[i])) - j[i].v) for i in range(n)])
            sj = np.abs(8 * math.pi * np.array(fixed_points(dj, dr)))
            rate = max(abs(run.saves[k]["cM"][p]) for p in BAND)
            lim = 1e-6 * rate * dt
            rows.append({"save": k, "t": run.saves[k]["t"], "band_max_8pi_dj": float(np.max(sj[1:])),
                         "cell0_8pi_dj": float(sj[0]), "band_max_16pi_drho": float(np.max(np.abs(
                             16 * math.pi * np.array(fixed_points(drho, dr)))[1:])),
                         "limit": lim, "pass": bool(np.max(sj[1:]) <= lim)})
    gate = [r for r in rows if r["t"] > 1e-9]
    return {"verdict": "PASS" if all(r["pass"] for r in gate) else "FAIL", "saves": rows}


def simpson_split(run):
    """A2 P10: stepper part = [C(t) - C(0)] - Simpson integral of Cdot, at each save, per fixed point."""
    ps = np.array(run.rec["per_step"], dtype=float)
    dt = run.dt
    out = []
    for sv in run.saves:
        n = sv["step"]
        row = {"t": sv["t"]}
        for ci, key in ((0, "H"), (1, "M")):
            C = ps[:, 1 + 5 * ci: 6 + 5 * ci]
            Cd = ps[:, 11 + 5 * ci: 16 + 5 * ci]
            if n == 0:
                semi = np.zeros(5)
            else:
                w = np.ones(n + 1); w[1:-1:2] = 4; w[2:-1:2] = 2
                semi = dt / 3 * (w[:, None] * Cd[:n + 1]).sum(axis=0)
            row[f"semi_{key}"] = semi.tolist()
            row[f"stepper_{key}"] = (C[n] - C[0] - semi).tolist()
        out.append(row)
    return out


def gs13(r3, r5, f0, f1):
    s3, s5 = simpson_split(r3), simpson_split(r5)
    gs1, gs3 = [], []
    for k in range(1, min(len(s3), len(s5))):
        for key in ("H", "M"):
            for p in BAND:
                F0, F1 = abs(f0.y(key, k, p)), abs(f1.y(key, k, p))
                a, b = s3[k][f"semi_{key}"][p], s5[k][f"semi_{key}"][p]
                gs1.append({"save": k, "C": key, "point": p, "semi_dt0": a, "semi_dt1": b,
                            "counted": abs(a) >= 10 * F0 and abs(b) >= 10 * F1, "pass": abs(a - b) <= 0.05 * abs(a)})
                x, y = s3[k][f"stepper_{key}"][p], s5[k][f"stepper_{key}"][p]
                ratio = y / x if x != 0 else float("nan")
                gs3.append({"save": k, "C": key, "point": p, "stepper_dt0": x, "stepper_dt1": y, "ratio": ratio,
                            "counted": abs(x) >= 10 * F0 and abs(y) >= 10 * F1, "pass": 0.2 <= ratio <= 0.3})
    return {"GS1": verdict(gs1), "GS3": verdict(gs3), "split_T2-3": s3, "split_T2-5": s5}


def main(stage1, stage2, cdot_floor=None):
    s1, d = Path(stage1), Path(stage2)
    fa = s1 / "arms"
    have = lambda nm: (d / f"{nm}.json").exists()
    runs = {nm: Run(d, nm) for nm in ("T2-3", "T2-4", "T2-6", "T2-E", "T2-5") if have(nm)}
    fl = lambda r: Run(fa, floor_arm(r.dr, r.dt, r.rec["r_max"]))
    cfl = json.loads(Path(cdot_floor).read_text()) if cdot_floor else None
    rb1 = json.loads((s1 / "readback.json").read_text())["arms"]
    rb2p = d / "readback.json"
    rb2 = json.loads(rb2p.read_text()) if rb2p.exists() else {"inputs": {}, "arms": {}}
    res = {"inputs": rb2["inputs"]}
    r = runs
    if {"T2-3", "T2-6"} <= r.keys():
        res["G2.2"] = g22(r["T2-3"], r["T2-6"], fl(r["T2-3"]))
    if {"T2-3", "T2-4"} <= r.keys():
        res["G2.3"] = g23(r["T2-3"], r["T2-4"], fl(r["T2-3"]), fl(r["T2-4"]))
    res["G2.4"] = {nm: g24(r[nm], fl(r[nm])) for nm in ("T2-3", "T2-6", "T2-4") if nm in r}
    if {"T2-3", "T2-E"} <= r.keys():
        res["edge_control"] = edge(r["T2-3"], r["T2-E"], fl(r["T2-3"]), fl(r["T2-E"]))
    res["GK.1"] = {}
    for nm in ("T2-3", "T2-4", "T2-5"):
        if nm in r and nm in rb2["arms"]:
            res["GK.1"][nm] = gk1(rb2["arms"][nm]["steps"], rb1[floor_arm(r[nm].dr, r[nm].dt, r[nm].rec["r_max"])]["steps"])
    if {"T2-3", "T2-4"} <= rb2["arms"].keys():
        a, b = rb2["arms"]["T2-3"]["steps"], rb2["arms"]["T2-4"]["steps"]
        res["GK.1_scaling_reading"] = [
            {"save": k, "primary_pred_ratio_T2-3_over_T2-4_div9": [x / (y / SCALE_T24) if y else None for x, y in
                                                                   zip(a[k]["GK1"]["primary"]["pred"], b[k]["GK1"]["primary"]["pred"])]}
            for k in range(min(len(a), len(b)))]
    if {"T2-3", "T2-6", "T2-4"} <= r.keys():
        res["GS.2"] = gs2(r["T2-3"], r["T2-6"], r["T2-4"], cfl)
    if "T2-3" in r:
        res["GM.1"] = gm1(d, r["T2-3"])
        res["stepper_split_reading_T2-3"] = simpson_split(r["T2-3"])
    if {"T2-3", "T2-5"} <= r.keys():
        res["GS.1_GS.3"] = gs13(r["T2-3"], r["T2-5"], fl(r["T2-3"]), fl(r["T2-5"]))
    res["not_scored_here"] = {"G2.1": "stage 3; needs the exact Klein-Gordon reference (design 3.2)"}
    res["summary"] = {k: (v.get("verdict") if isinstance(v, dict) and "verdict" in v else
                          {kk: vv.get("verdict") for kk, vv in v.items() if isinstance(vv, dict)})
                      for k, v in res.items() if k.startswith(("G", "edge")) and isinstance(v, dict)}
    (d / "score2.json").write_text(json.dumps(res, indent=1, default=float))
    print(json.dumps(res["summary"], indent=1))


if __name__ == "__main__":
    main(*sys.argv[1:4])
