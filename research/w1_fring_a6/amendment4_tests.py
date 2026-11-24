"""Score Amendment 4 predictions P1, P2, P3, P5 on stored data (methods in
results/research/w1_fring_a6_amendment4/PRERUN.md). No simulation.

Usage: python -m research.w1_fring_a6.amendment4_tests ASTER_SP_DIR OUT.json
ASTER_SP_DIR holds the #598 replay samples {dr}-{dt}.npz.
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline

ARMS = Path("results/research/w1_fring_a6/arms")
ORIGIN = Path("results/research/w1_fring_origin_trace")
REPLAY = Path("results/research/w1_fring_aster_replay/bands_both_ways.json")
RADII = np.linspace(0.05, 0.25, 17)
BAND = [1, 2, 3, 4]
DR = 0.05
SAVES = ["0.8235", "1.6455", "3.291"]


def arm(name):
    d = json.loads((ARMS / f"{name}.json").read_text())
    ps = d["per_step"]
    t = np.array([x["t"] for x in ps])
    H = np.array([x["H"][:6] for x in ps])
    M = np.array([x["M"][:6] for x in ps])
    return d, t, H, M


def r_cells(dr, n=6):
    return (np.arange(n) + 0.5) * dr


def step_of(t, tt):
    return int(np.argmin(np.abs(t - tt)))


def dtH(t, H, n):
    if n + 1 < len(t):
        return (H[n + 1] - H[n - 1]) / (t[n + 1] - t[n - 1])
    # last step (t 3.291): second-order one-sided (amendment to PRERUN, made before reading)
    return (3 * H[n] - 4 * H[n - 1] + H[n - 2]) / (t[n] - t[n - 2])


# ---------- scorers (pure, used on data and on synthetic inputs) ----------
def score_p1(Q, sign_ok):
    Q = np.abs(np.asarray(Q)); s = np.asarray(sign_ok)
    if np.all(Q <= 0.5) and np.all(s):
        return "PASS"
    if np.mean(Q > 1) > 0.5 or np.mean(~s) > 0.5:
        return "FAIL"
    return "PARTIAL"


def near(x, lo, hi):
    return lo <= x <= hi


def score_p2(early, late):
    if near(late, 0.18, 0.35) or late > 1.25:
        return "FAIL"
    if near(early, 0.18, 0.35) and near(late, 0.8, 1.25):
        return "PASS"
    return "PARTIAL"


def score_p3(peak_c, peak_f, r_c, r_f, g_cell1, g_peak):
    s3 = peak_c == 1 and peak_f == 1 and near(g_cell1, 2, 8)
    same = abs(r_c - r_f) <= 0.025
    s2 = same and near(g_peak, 1 / 8, 1 / 2)
    both_fail = same and g_peak > 1
    return {"S3_holds": bool(s3), "S2_holds": bool(s2), "both_fail": bool(both_fail),
            "verdict": "S3" if s3 else "S2" if s2 else "BOTH FAIL" if both_fail else "NEITHER AS WRITTEN"}


def score_p5(cases):
    agree = [(gm > 1) == (gh > 1) and near(gm / gh, 0.5, 2) for gm, gh in cases]
    bad = sum(gm > 1 and gh < 1 for gm, gh in cases)
    if bad >= 2:
        return "FAIL"
    return "PASS" if all(agree) else "PARTIAL"


def synthetic():
    out = {}
    out["p1"] = [score_p1([0.1, -0.2], [True, True]), score_p1([2, 3, 0.1], [True, True, True]),
                 score_p1([0.1, 0.2, 0.3], [False, False, True])]
    out["p2"] = [score_p2(0.25, 1.0), score_p2(0.25, 0.26), score_p2(0.25, 2.0)]
    out["p3"] = [score_p3(1, 1, .075, .0375, 4, 4)["verdict"], score_p3(2, 4, .125, .1125, 1, .25)["verdict"],
                 score_p3(2, 4, .125, .1125, 1, 3)["verdict"]]
    out["p5"] = [score_p5([(2, 2.5)] * 4), score_p5([(2, .5)] * 4)]
    ok = out["p1"] == ["PASS", "FAIL", "FAIL"] and out["p2"] == ["PASS", "FAIL", "FAIL"] \
        and out["p3"] == ["S3", "S2", "BOTH FAIL"] and out["p5"] == ["PASS", "FAIL"]
    return out, ok


# ---------- data ----------
def input_checks():
    res = {}
    for a in ("A1", "A2"):
        d, t, H, M = arm(a)
        diffs = []
        for k, s in d["saves"].items():
            n = step_of(t, s["t"])
            diffs.append(float(max(np.max(np.abs(H[n] - s["row"]["H"][:6])), np.max(np.abs(M[n] - s["row"]["M"][:6])))))
        res[f"{a}_per_step_vs_saves_max_abs_diff"] = max(diffs)
    o = json.loads((ORIGIN / "A-0.05.json").read_text())["rows"][-1]
    _, t, H, M = arm("A1")
    n = step_of(t, o["t"])
    res["origin_A-0.05_vs_A1_t0.015_max_abs_diff"] = float(max(np.max(np.abs(H[n] - o["cells"]["H"][:6])),
                                                               np.max(np.abs(M[n] - o["cells"]["M"][:6]))))
    return res


def p1_p2():
    out = {"P1": {}, "P2": {}}
    r = r_cells(DR)[BAND]
    data = {a: arm(a) for a in ("A1", "A2")}
    Qs, signs = [], []
    for a, (d, t, H, M) in data.items():
        for k in SAVES:
            n = step_of(t, float(k))
            dh = dtH(t, H, n)[BAND]; m = M[n][BAND]
            Q = (4 * m / r + dh) / np.abs(dh)
            sg = np.sign(m) == np.sign(-dh)
            Qs += list(Q); signs += list(sg)
            out["P1"][f"{a}@{k}"] = {"M": m.tolist(), "d_tH": dh.tolist(), "minus_r_over_4_d_tH": (-r / 4 * dh).tolist(),
                                     "Q": Q.tolist(), "sign_M_eq_sign_minus_dtH": sg.tolist()}
        # secondary: over saves at 1.6455
        s = d["saves"]
        Hs = {k: np.array(s[k]["row"]["H"][1:5]) for k in SAVES}
        dh = (Hs["3.291"] - Hs["0.8235"]) / (s["3.291"]["t"] - s["0.8235"]["t"])
        m = np.array(s["1.6455"]["row"]["M"][1:5])
        out["P1"][f"{a}@1.6455_over_saves"] = {"d_tH": dh.tolist(), "Q": ((4 * m / r + dh) / np.abs(dh)).tolist()}
        # supplementary: all steps t >= .8235
        ns = [n for n in range(1, len(t) - 1) if t[n] >= 0.8235 - 1e-9]
        allQ = np.array([(4 * M[n][BAND] / r + dtH(t, H, n)[BAND]) / np.abs(dtH(t, H, n)[BAND]) for n in ns])
        allS = np.array([np.sign(M[n][BAND]) == np.sign(-dtH(t, H, n)[BAND]) for n in ns])
        out["P1"][f"{a}_all_steps_t_ge_0.8235"] = {
            "points": int(allQ.size), "share_absQ_le_0.5": float(np.mean(np.abs(allQ) <= 0.5)),
            "share_absQ_gt_1": float(np.mean(np.abs(allQ) > 1)), "share_sign_agree": float(np.mean(allS)),
            "median_absQ": float(np.median(np.abs(allQ)))}
    out["P1"]["max_absQ_at_saves"] = float(np.max(np.abs(Qs)))
    out["P1"]["share_absQ_gt_1_at_saves"] = float(np.mean(np.abs(np.array(Qs)) > 1))
    out["P1"]["share_sign_agree_at_saves"] = float(np.mean(signs))
    out["P1"]["score"] = score_p1(Qs, signs)
    # P2
    SH = {}
    for a, (d, t, H, M) in data.items():
        for k in ("0.015", "3.291"):
            n = step_of(t, float(k))
            SH[(a, k)] = dtH(t, H, n)[BAND] + 4 * M[n][BAND] / r
            out["P2"][f"S_H_{a}@{k}"] = SH[(a, k)].tolist()
    ratio = {}
    for k in ("0.015", "3.291"):
        ratio[k] = float(np.max(np.abs(SH[("A2", k)])) / np.max(np.abs(SH[("A1", k)])))
        out["P2"][f"per_cell_ratio_A2_over_A1@{k}"] = (np.abs(SH[("A2", k)]) / np.abs(SH[("A1", k)])).tolist()
    out["P2"]["band_max_ratio_A2_over_A1"] = ratio
    out["P2"]["score"] = score_p2(ratio["0.015"], ratio["3.291"])
    return out


def p3():
    out = {}
    _, t2, H2, _ = arm("A2")
    _, t1, H1, _ = arm("A1")
    fine = json.loads((ORIGIN / "A-0.025.json").read_text())
    assert abs(fine["dt"] - 0.00075) < 1e-15
    fr = fine["rows"]
    Hf = np.array(fr[-1]["cells"]["H"][:6]); Hf0 = np.array(fr[0]["cells"]["H"][:6])
    rc, rf = r_cells(0.05), r_cells(0.025)

    def compare(Hc, Hf_):
        pc = 1 + int(np.argmax(np.abs(Hc[1:6]))); pf = 1 + int(np.argmax(np.abs(Hf_[1:6])))
        g1 = abs(Hf_[1]) / abs(Hc[1]); gp = abs(Hf_[pf]) / abs(Hc[pc])
        return {"H_coarse_cells_0_5": Hc.tolist(), "H_fine_cells_0_5": Hf_.tolist(),
                "peak_cell_coarse": pc, "peak_r_coarse": float(rc[pc]), "peak_cell_fine": pf, "peak_r_fine": float(rf[pf]),
                "cell1_growth": float(g1), "peak_growth": float(gp),
                **score_p3(pc, pf, rc[pc], rf[pf], g1, gp)}

    n2 = step_of(t2, 0.015); n1 = step_of(t1, 0.015)
    out["matched_dt_0.00075"] = compare(H2[n2], Hf)
    out["score"] = out["matched_dt_0.00075"]["verdict"]
    out["supp_matched_dt_over_dr_A1"] = compare(H1[n1], Hf)
    out["supp_minus_t0_matched_dt"] = compare(H2[n2] - H2[0], Hf - Hf0)
    return out


def band17(R, Q):
    m = R < 0.4
    return CubicSpline(R[m], Q[m])(RADII)


def p5(sp):
    out = {"cases": {}, "replay_input_check": {}}
    replay = json.loads(REPLAY.read_text())["runs"]
    cases = []
    for dt in ("0.000375", "0.0001875"):
        z = {dr: np.load(Path(sp) / f"{dr}-{dt}.npz") for dr in ("0.025", "0.0125")}
        for dr, zz in z.items():
            got = float(np.max(np.abs(band17(zz["R_areal"][4], zz["M"][4]))))
            want = replay[f"{dr}-{dt}"]["band_17_areal_radii_cubic"]
            out["replay_input_check"][f"{dr}-{dt}"] = {"recomputed": got, "stored": want, "rel_diff": abs(got - want) / want}
        vals = {}
        for dr, zz in z.items():
            T = zz["t"]; R = zz["R_areal"]; H = zz["H"]; M = zz["M"]
            h = [band17(R[i], H[i]) for i in range(5)]
            rr = zz["r"]; cc = (rr > 0.05) & (rr < 0.25)
            vals[dr] = {
                "centred@0.01125": (float(np.max(np.abs(band17(R[3], M[3])))), float(np.max(np.abs((h[4] - h[2]) / (T[4] - T[2])))),
                                     float(np.max(np.abs(M[3][cc]))), float(np.max(np.abs((H[4][cc] - H[2][cc]) / (T[4] - T[2]))))),
                "backward@0.015": (float(np.max(np.abs(band17(R[4], M[4])))), float(np.max(np.abs((h[4] - h[3]) / (T[4] - T[3])))),
                                    float(np.max(np.abs(M[4][cc]))), float(np.max(np.abs((H[4][cc] - H[3][cc]) / (T[4] - T[3]))))),
            }
        for w in ("centred@0.01125", "backward@0.015"):
            c, f = vals["0.025"][w], vals["0.0125"][w]
            gm, gh = f[0] / c[0], f[1] / c[1]
            cases.append((gm, gh))
            out["cases"][f"dt{dt}_{w}"] = {
                "band_M_coarse": c[0], "band_M_fine": f[0], "band_dtH_coarse": c[1], "band_dtH_fine": f[1],
                "M_growth": gm, "dtH_growth": gh, "ratio_of_growths": gm / gh,
                "supp_cell_centres": {"M_growth": f[2] / c[2], "dtH_growth": f[3] / c[3]}}
    out["score"] = score_p5(cases)
    return out


def main(sp, outp):
    syn, ok = synthetic()
    res = {"synthetic": syn, "synthetic_ok": ok, "input_checks": input_checks()}
    res.update(p1_p2()); res["P3"] = p3(); res["P5"] = p5(sp)
    res["P4"] = "skipped: needs a new run with the proposed centre step (his go)"
    Path(outp).write_text(json.dumps(res, indent=1))
    print(json.dumps({"synthetic_ok": ok, "input_checks": res["input_checks"],
                      "P5_replay_check": res["P5"]["replay_input_check"],
                      "scores": {k: res[k]["score"] for k in ("P1", "P2", "P3", "P5")}}, indent=1))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
