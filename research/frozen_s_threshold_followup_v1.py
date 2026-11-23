"""Frozen-S threshold scan follow-up (predeclared design: Drive DESIGN-frozen-S-threshold-followup-Nov23-2026,
1C9wW5BRJi_1f0-obh5_ljdRcNdC1biA2EbbIQYBk7Y8; pre-run log entry 2026-11-23). Research only.
Imports run() from research/frozen_s_threshold_scan_v1.py unchanged (estimator B, CAP, floor, LAM_GROW as they are);
passes T = 10000 and the design's grid. The seed floor run() computes (np.median over t in [60, 120]) is recorded by a
wrapper around the module's numpy median; run()'s source and arithmetic are untouched.
python research/frozen_s_threshold_followup_v1.py [outdir]"""
import math, os, sys, time
from multiprocessing import Pool
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import frozen_s_threshold_scan_v1 as S

T = 10000.0; G = 2e-3
FINE = [round(5.5 + 0.1 * i, 2) for i in range(16)]


class _NP:  # numpy proxy: records the one np.median call run() makes (the seed floor); everything else passes through
    def __init__(self): self.med = []
    def __getattr__(self, k): return getattr(np, k)
    def median(self, *a, **kw):
        v = np.median(*a, **kw); self.med.append(float(v)); return v


def job(a):
    p = _NP(); S.np = p; t0 = time.time()
    o = S.run(T=T, **a); o["floor"] = p.med[-1] if p.med else float("nan"); o["sec"] = time.time() - t0
    return a, o


def K(a):
    return tuple(sorted(a.items()))


def cls(o):
    if math.isnan(o["lam"]): return "U"
    return "G" if o["lam"] > G else "N"


def bracket(rows):  # rows: [(xAR, o)] sorted by xAR -> (lo, hi, widened, notes)
    c = [cls(o) for _, o in rows]; xs = [x for x, _ in rows]; nu = [i for i, k in enumerate(c) if k != "U"]
    hi = None
    for j, i in enumerate(nu):
        if all(c[m] == "G" for m in nu[j:]): hi = i; break
    lo = None
    if hi is not None:
        n_below = [i for i in nu if i < hi and c[i] == "N"]
        lo = n_below[-1] if n_below else None
    nonmono = any(c[i] == "N" and any(c[m] == "G" for m in nu if m < i) for i in nu)
    widened = nonmono or (lo is not None and hi is not None and hi - lo > 1)
    return (xs[lo] if lo is not None else None), (xs[hi] if hi is not None else None), widened, nonmono, \
        [xs[i] for i, k in enumerate(c) if k == "U"]


def main():
    od = sys.argv[1] if len(sys.argv) > 1 else "results/frozen-s-threshold-followup-2026-11-23"
    os.makedirs(od, exist_ok=True); t0 = time.time()
    ph1 = [dict(R=12, xAR=x) for x in FINE] + [dict(R=24, xAR=x) for x in sorted(FINE + [6.72])]
    ph1 += [dict(R=12, xAR=6.0, pump=False), dict(R=24, xAR=6.5, pump=False)]
    ph1 += [dict(R=24, xAR=x, wd=0.3) for x in FINE]
    ph1 += [dict(R=12, xAR=6.5, H=0.01), dict(R=12, xAR=6.5, H=0.004)]
    with Pool(4) as p:
        D = {K(a): o for a, o in p.map(job, ph1, chunksize=1)}
    row = lambda **a: D[K(a)]
    prim = {R: [(x, row(R=R, xAR=x)) for x in (sorted(FINE + [6.72]) if R == 24 else FINE)] for R in (12, 24)}
    br = {R: bracket(prim[R]) for R in (12, 24)}
    w3 = [(x, row(R=24, xAR=x, wd=0.3)) for x in FINE]; br3 = bracket(w3)
    # phase 2: rulers on the bracket rows, chosen by the fixed rule "the two rows forming the bracket"
    ph2 = []
    for R in (12, 24):
        for x in br[R][:2]:
            if x is None: continue
            ph2.append(dict(R=R, xAR=x, dr=0.05))
            if R == 24: ph2.append(dict(R=R, xAR=x, sw=60.0))
    for x in br3[:2]:
        if x is not None: ph2.append(dict(R=24, xAR=x, wd=0.3, dr=0.05))
    with Pool(4) as p:
        D.update({K(a): o for a, o in p.map(job, ph2, chunksize=1)})
    nruns = len(ph1) + len(ph2)
    f = lambda o: ("nan" if math.isnan(o["lam"]) else f"{o['lam']:+.4f}") + f"[A {o['lamA']:+.4f}]"
    L = [f"# frozen-S threshold follow-up; T = {T:.0f}, estimator B, growth = lambda_B > 2e-3; {nruns} runs; wall {{WALL}} s",
         "# columns: lambda_B [A lambda_A] class span(e-folds) floor floor/est(0.858 R d0) t_end"]
    def line(x, o, R, d0=1e-6):
        est = 0.858 * R * d0
        return (f"  xAR {x:5.2f}  {f(o)}  {cls(o)}  span {o['span']:.2f}  floor {o['floor']:.2e}  "
                f"{o['floor'] / est:.2f}  t_end {o['tend']:.0f}")
    V = {}
    for R in (12, 24):
        L.append(f"\n## primary R {R} (dr 0.1, wd 1.0, sponge 40, H 0, seed 1e-6)")
        L += [line(x, o, R) for x, o in prim[R]]
        lo, hi, wid, nm, uu = br[R]
        L.append(f"  BRACKET ({lo}, {hi})" + (" NON-MONOTONIC" if nm else "") + (f" unresolved rows {uu}" if uu else "")
                 + (" widened" if wid else ""))
        if lo is not None and hi is not None:
            ol, oh = row(R=R, xAR=lo)["lam"], row(R=R, xAR=hi)["lam"]
            L.append(f"  linear interpolation to 2e-3 (Account): {lo + (G - ol) / (oh - ol) * (hi - lo):.3f}")
    # tests
    lam = lambda R, x, **k: row(R=R, xAR=x, **k)["lam"]
    P = {12: {}, 24: {}}
    P1 = {12: [5.5], 24: [5.5, 5.6, 5.7]}
    P2 = {12: [round(6.1 + 0.1 * i, 2) for i in range(10)], 24: sorted([round(6.3 + 0.1 * i, 2) for i in range(8)] + [6.72])}
    L.append("\n## predeclared tests")
    for R in (12, 24):
        bad1 = [x for x in P1[R] if not lam(R, x) <= G]; P[R]["P1"] = not bad1
        bad2 = [x for x in P2[R] if not (cls(row(R=R, xAR=x)) == "G" and row(R=R, xAR=x)["span"] >= 1.0)]; P[R]["P2"] = not bad2
        gr = [(x, o["lam"]) for x, o in prim[R] if cls(o) == "G"]
        bad3 = [(a[0], b[0]) for a, b in zip(gr, gr[1:]) if b[1] < 0.95 * a[1]]; P[R]["P3"] = not bad3
        L.append(f"  R {R} P1 no growth at {P1[R]}: {'PASS' if not bad1 else 'FAIL ' + str(bad1)}")
        L.append(f"  R {R} P2 growth at {P2[R][0]}-{P2[R][-1]}: {'PASS' if not bad2 else 'FAIL ' + str(bad2)}")
        L.append(f"  R {R} P3 monotonic (no fall > 5%): {'PASS' if not bad3 else 'FAIL ' + str(bad3)}")
    l672 = lam(24, 6.72); P[24]["P4"] = 0.0049 <= l672 <= 0.0091
    L.append(f"  R 24 P4 witness 6.72: {l672:.4f} vs 0.007 (0.0049-0.0091): {'PASS' if P[24]['P4'] else 'FAIL'}")
    ref5 = {12: [(6.0, 0.0048), (6.5, 0.0114), (7.0, 0.0178)], 24: [(6.5, 0.0048), (6.72, 0.0062), (7.0, 0.0079)]}
    for R in (12, 24):
        rr = [(x, lam(R, x), v, abs(lam(R, x) - v) <= 0.1 * v) for x, v in ref5[R]]; P[R]["P5"] = all(q[3] for q in rr)
        L.append(f"  R {R} P5 continuity: " + "; ".join(f"{x}: {a:.4f} vs {v} {'ok' if k else 'OUT'}" for x, a, v, k in rr)
                 + f": {'PASS' if P[R]['P5'] else 'FAIL'}")
    c6 = {12: lam(12, 6.0, pump=False), 24: lam(24, 6.5, pump=False)}
    for R in (12, 24):
        P[R]["P6"] = c6[R] <= G
        o = row(R=R, xAR=6.0 if R == 12 else 6.5, pump=False)
        L.append(f"  R {R} P6 control pump off: {f(o)} {cls(o)} floor {o['floor']:.2e}: {'PASS' if P[R]['P6'] else 'FAIL'}")
    def p7(R, x, **k):
        a, b = row(R=R, xAR=x, **{kk: vv for kk, vv in k.items() if kk == 'wd'}), row(R=R, xAR=x, **k)
        same = cls(a) == cls(b); within = abs(b["lam"] - a["lam"]) <= 0.1 * abs(a["lam"]) if not math.isnan(a["lam"] + b["lam"]) else False
        return same, within, a, b
    for R in (12, 24):
        ok = True; moved = False; txt = []
        for x in br[R][:2]:
            if x is None: continue
            for k in ([dict(dr=0.05)] + ([dict(sw=60.0)] if R == 24 else [])):
                s, w, a, b = p7(R, x, **k); ok &= s and w; moved |= not s
                txt.append(f"{x} {list(k.items())[0]}: {f(b)} {cls(b)} vs {f(a)} {cls(a)} {'ok' if s and w else 'OUT'}")
        P[R]["P7"] = ok; P[R]["P7moved"] = moved
        L.append(f"  R {R} P7 rulers: " + "; ".join(txt) + f": {'PASS' if ok else 'FAIL'}")
    pred = {12: 5.78, 24: 6.0}
    for R in (12, 24):
        lo, hi, wid, nm, uu = br[R]
        if lo is None or hi is None:
            P[R]["B1"] = False; L.append(f"  R {R} B1: bracket open ({lo}, {hi}): FAIL"); continue
        P[R]["B1"] = lo - 0.2 <= pred[R] <= hi + 0.2
        b07 = lo - 0.7 <= pred[R] <= hi + 0.7
        L.append(f"  R {R} B1 prediction {pred[R]} in ({lo}-0.2, {hi}+0.2): {'PASS' if P[R]['B1'] else 'FAIL'}"
                 f"   [continuity, not scored: widened by 0.7 {'inside' if b07 else 'outside'}]")
    L.append("\n## verdict per R")
    for R in (12, 24):
        lo, hi, wid, nm, uu = br[R]
        need = ["P1", "P2", "P3", "P6", "P7", "B1"] + (["P4"] if R == 24 else [])
        if not P[R]["P1"] or not P[R]["P6"] or not P[R]["B1"]: v = "FAIL"
        elif all(P[R][k] for k in need) and P[R]["P5"] and not wid and not uu: v = "PASS"
        else: v = "MIXED"
        V[R] = v
        L.append(f"  R {R}: {v}  (bracket ({lo}, {hi}); " + ", ".join(f"{k} {'pass' if P[R][k] else 'fail'}" for k in need + ['P5'])
                 + (", bracket widened" if wid else "") + (", ruler moved a class" if P[R]["P7moved"] else "") + ")")
    # item 3
    L.append("\n## item 3: edge width 0.3 (R 24, T 10000)")
    L += [line(x, o, 24) for x, o in w3]
    lo3, hi3, wid3, nm3, uu3 = br3
    L.append(f"  BRACKET wd 0.3 ({lo3}, {hi3})" + (" NON-MONOTONIC" if nm3 else "") + (f" unresolved {uu3}" if uu3 else ""))
    hi1 = bracket([(x, o) for x, o in prim[24] if x != 6.72])[1]
    res3 = []
    for x in (lo3, hi3):
        if x is None: continue
        s, w, a, b = p7(24, x, wd=0.3, dr=0.05)
        res3.append(s); L.append(f"  dr 0.05 on wd 0.3 row {x}: {f(b)} {cls(b)} vs {f(a)} {cls(a)}")
    if hi3 is None or hi1 is None or not all(res3):
        v3 = "UNRESOLVED"; d = None
    else:
        d = round(hi3 - hi1, 2)
        v3 = ("CONFIRMED" if 0.2 <= d <= 0.4 else "NO SHIFT (prediction ~6.3 NOT CONFIRMED)" if -0.1 <= d <= 0.1 else
              "OPPOSITE SIGN (prediction FAILS)" if d <= -0.2 else "LARGER THAN PREDICTED (prediction FAILS on size)")
    L.append(f"  first growth row wd 1.0 (0.1 grid) {hi1}; wd 0.3 {hi3}; Delta {d}: {v3}")
    if br[24][1] == 6.72:
        L.append("  note: 6.72 is the wd 1.0 first growth row; Delta on the 0.1 grid is the predeclared reading")
    # item 4
    L.append("\n## item 4 (R 12, xAR 6.5)")
    o0, oa, ob = row(R=12, xAR=6.5), row(R=12, xAR=6.5, H=0.01), row(R=12, xAR=6.5, H=0.004)
    L.append(f"  H 0: {f(o0)}; H 0.01: {f(oa)} floor {oa['floor']:.2e}; H 0.004: {f(ob)} floor {ob['floor']:.2e} span {ob['span']:.2f}")
    a4 = oa["lam"] <= G; p4b = o0["lam"] - 0.006; b4 = abs(ob["lam"] - p4b) <= 0.004
    L.append(f"  4a no growth at H 0.01 (predicted -0.0043): {'PASS' if a4 else 'FAIL'}")
    L.append(f"  4b H 0.004: {ob['lam']:.4f} vs predicted {p4b:.4f} (+/- 0.004): {'PASS' if b4 else 'FAIL'}")
    lowf = [(K(a), o["floor"]) for a, o in [(dict(kk), D[kk]) for kk in D]
            if a.get("H", 0) == 0 and o["floor"] < 0.5 * 0.858 * a["R"] * 1e-6]
    L.append(f"\n## floors below 0.5 x estimate (H 0 runs): {lowf if lowf else 'none'}")
    L.append(f"## slowest run {max(o['sec'] for o in D.values()):.0f} s")
    txt = "\n".join(L).replace("{WALL}", f"{time.time() - t0:.0f}")
    open(os.path.join(od, "summary.txt"), "w").write(txt + "\n"); print(txt)


if __name__ == "__main__":
    main()
