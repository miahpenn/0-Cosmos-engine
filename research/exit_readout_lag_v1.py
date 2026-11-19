"""RESEARCH ONLY. Exit-readout lag check (design: coordination log, Nov 18 2026). Frozen-S lattice G read at the old (rolling-max) and lag-free (energy) exit times."""
import math, sys, numpy as np
sys.path.insert(0, "research")
import pumping_spatial_toyS_v1 as P
import frozen_s_trapped_v1 as F
DT, WIN = P.DT, P.WIN
def exits(H, x0):
    sh, ps, t = math.sqrt(x0), 0.0, 0.0
    ac = lambda s, p: -3 * H * p - (s + 0.5 * s ** 3)
    A = ac(sh, ps); hist = []; told = tnew = None
    while t < 400 and (told is None or tnew is None):
        ps += 0.5 * DT * A; sh += DT * ps; t += DT; A = ac(sh, ps); ps += 0.5 * DT * A
        hist.append(sh * sh)
        if tnew is None:
            E = sh * sh / 2 + sh ** 4 / 8 + ps * ps / 2
            if -2 + math.sqrt(4 + 8 * E) < 4 / 3: tnew = t
        if told is None and len(hist) > WIN and max(hist[-WIN:]) < 4 / 3: told = t
    return told, tnew
if __name__ == "__main__":
    print("(H,x0): te_old (script) | lattice P.evolve te | te_new | lag | G(te_old) | G(te_new) | ratio")
    ok_lag = True; G = {}
    for H, x0 in ((0.0057, 2.0), (0.0057, 2.67), (0.01, 2.0), (0.01, 2.67)):
        to, tn = exits(H, x0); tl = P.evolve(H, x0, None)
        go = F.run(H, x0, 1e-8, to, "frozen")["Gk"]; gn = F.run(H, x0, 1e-8, tn, "frozen")["Gk"]; G[(H, x0)] = (go, gn)
        ok_lag &= 0.6 <= to - tn <= 3.2
        print(f"H {H} x0 {x0}: {to:.2f} | {tl:.2f} (diff {to-tl:+.2f}) | {tn:.2f} | lag {to-tn:.2f} | {go:.5g} | {gn:.5g} | {gn/go:.3f}", flush=True)
    print("LAG AS DESCRIBED (all rows in [0.6,3.2]):", ok_lag)
    for key, f in (((0.0057, 2.0), 620.0), ((0.0057, 2.67), 9e4)):
        go, gn = G[key]; eo, en = go / f - 1, gn / f - 1
        print(f"formula {f:g} at {key}: gap old {eo:+.1%} new {en:+.1%} -> {'nearer' if abs(en) < abs(eo) else 'not nearer'}")
    print("H_G10 at x0=2.0 (log-linear interpolation):")
    for lab, idx in (("old", 0), ("new", 1)):
        Hs = np.arange(0.012, 0.0301, 0.002); lg = []
        for H in Hs:
            e = exits(H, 2.0)[idx]; lg.append(math.log(F.run(H, 2.0, 1e-8, e, "frozen")["Gk"]))
        lg = np.array(lg); target = math.log(10)
        h10 = float(np.interp(target, lg[::-1], Hs[::-1])) if lg.min() < target < lg.max() else float("nan")
        print(f"  {lab}: H_G10 {h10:.4f}   (G at H=0.012..0.03: {[round(math.exp(v),1) for v in lg]})", flush=True)
