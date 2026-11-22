"""Frozen-S threshold scan (design: coordination log, Lumen 2026-11-22 + amendments A-C, Projects' entry 2026-11-22).
Linear homogeneous radial D (chi = r d, l = 0) in a prescribed top-hat S pump:
  chi_tt + (3H + gamma) chi_t = chi_rr - [M2(r) + 0.75 x_A prof cos(2 w t)] chi,  prof = (1 - tanh((r-R)/wd))/2,
  M2 = 2 - 1.5 x_A prof,  w^2 = 4 - 1.875 x_A.  Dirichlet at r = 0 and r_max = R + 100; quadratic sponge over the outer sw.
estimator A (declared): LSQ slope of ln(max|chi| per pump period pi/w) for t >= 60 until max|chi| >= 1e-3 or T = 600; estimator B (amended before the full run, see log): the same fit starting where max|chi| first reaches 3x the seed floor.
python research/frozen_s_threshold_scan_v1.py [outdir]"""
import math, os, sys
from multiprocessing import Pool
import numpy as np

T_MAX = 600.0; T_FIT = 60.0; CAP = 1e-3; LAM_GROW = 2e-3


def run(R, xAR, H=0.0, d0=1e-6, dr=0.1, wd=1.0, sw=40.0, detune=1.0, pump=True, T=T_MAX):
    xA = xAR / R; rmax = R + 100.0; N = int(round(rmax / dr)); r = np.arange(N + 1) * dr
    prof = 0.5 * (1 - np.tanh((r - R) / wd)); M2 = 2.0 - 1.5 * xA * prof; P = (0.75 * xA * prof) if pump else 0 * prof
    w = math.sqrt(4 - 1.875 * xA); wp = w * detune
    r0 = rmax - sw; gam = np.where(r > r0, ((r - r0) / sw) ** 2, 0.0)
    dt = 0.4 * dr; ns = int(round(T / dt)); a = (3 * H + gam) * dt / 2
    chi = r * d0 * np.exp(-(r / (2 * R)) ** 2); chi[0] = 0; chi[-1] = 0

    def rhs(c, t):
        out = np.zeros_like(c); out[1:-1] = (c[2:] - 2 * c[1:-1] + c[:-2]) / dr ** 2
        return out - (M2 + P * math.cos(2 * wp * t)) * c
    prev = chi + 0.5 * dt * dt * rhs(chi, 0.0)
    per = math.pi / w; t_list = []; m_list = []; cur = 0.0; tcur = per; capped = False
    for n in range(ns):
        t = n * dt
        nxt = (2 * chi - (1 - a) * prev + dt * dt * rhs(chi, t)) / (1 + a); nxt[0] = 0; nxt[-1] = 0
        prev, chi = chi, nxt
        cur = max(cur, float(np.abs(chi).max()))
        if t + dt >= tcur:
            t_list.append(tcur - 0.5 * per); m_list.append(cur); cur = 0.0; tcur += per
            if t > T_FIT and m_list[-1] >= CAP:
                break
    t = np.array(t_list); m = np.array(m_list); k = t >= T_FIT
    if k.sum() < 8:
        return dict(lam=float("nan"), lamA=float("nan"), span=0.0, ok=False, tend=float(t[-1]) if len(t) else 0.0)
    lg = np.log(m[k]); lamA = float(np.polyfit(t[k], lg, 1)[0]); spanA = float(abs(lamA) * (t[k][-1] - t[k][0]))
    # estimator B (amended before the full run): fit from the first period where max|chi| >= 3 x the seed floor (median of
    # max|chi| over t in [60, 120]) until the cap or T; needs >= 1 e-fold in that window, else 'unresolved' (= no growth).
    fl = np.median(m[(t >= T_FIT) & (t <= 2 * T_FIT)]); idx = np.nonzero((t >= T_FIT) & (m >= 3 * fl))[0]
    lam = float("nan"); span = 0.0; ok = False
    if idx.size >= 8:
        tb, mb = t[idx[0]:], m[idx[0]:]
        if len(tb) >= 8:
            lam = float(np.polyfit(tb, np.log(mb), 1)[0]); span = float(np.log(mb.max() / mb[:3].mean()))
            ok = span >= 1.0
    if not ok:
        lam = 0.0 if lamA <= LAM_GROW else float("nan")      # no emergence from the floor: no growth (A kept for the record)
    return dict(lam=lam, lamA=lamA, span=span, ok=True if lam == 0.0 else ok, tend=float(t[-1]), floor=bool(d0 <= 1e-12))


def job(a):
    return a, run(**a)


def main():
    od = sys.argv[1] if len(sys.argv) > 1 else "results/frozen-s-threshold-2026-11-22"
    os.makedirs(od, exist_ok=True)
    G = [4.0, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 10.0]; extra24 = [6.72, 8.64, 12.24]
    jobs = []
    for R in (12, 24):
        for x in G + (extra24 if R == 24 else []):
            for d0 in (1e-6, 1e-12):
                for H in (0.0, 0.01):
                    jobs.append(dict(R=R, xAR=x, H=H, d0=d0))
    for R in (12, 24):
        for x in (5.5, 6.0, 6.5):
            jobs.append(dict(R=R, xAR=x, dr=0.05)); jobs.append(dict(R=R, xAR=x, sw=60.0))
    for x in G[3:7]:
        jobs.append(dict(R=24, xAR=x, wd=0.3))
    jobs += [dict(R=24, xAR=8.64, pump=False), dict(R=24, xAR=8.64, detune=1.1)]
    for R, xs in ((24, (5.5, 6.0, 6.5, 6.72, 7.0)), (12, (5.5, 6.0, 6.5))):
        for x in xs:
            jobs.append(dict(R=R, xAR=x, T=2000.0))
    with Pool(2) as p:
        res = p.map(job, jobs, chunksize=2)
    key = lambda a: (a["R"], a["xAR"], a.get("H", 0.0), a.get("d0", 1e-6), a.get("dr", 0.1), a.get("wd", 1.0), a.get("sw", 40.0), a.get("detune", 1.0), a.get("pump", True), a.get("T", T_MAX))  # noqa
    D = {key(a): o for a, o in res}
    L = ["# frozen-S threshold scan; lambda per unit time; 'unres' = window < 1 e-fold; infinite-medium rate 0.75 x_A/(4 w) printed per row (x_A = xAR/R)"]
    def row(R, x, **kw):
        k = key(dict(R=R, xAR=x, **kw)); return D[k]
    def f(o):
        return "nan" if math.isnan(o["lam"]) else f"{o['lam']:+.4f}" + ("" if o["ok"] else "(unres)") + f"[A {o['lamA']:+.4f}]"
    L.append("\n## primary grid (dr 0.1, wd 1.0, sponge 40): lambda for seed1e-6 H0 | seed1e-12 H0 | seed1e-6 H0.01 | seed1e-12 H0.01 ; ref")
    thr = {}
    for R in (12, 24):
        L.append(f"R {R}")
        prev = None; found = None
        for x in G + (extra24 if R == 24 else []):
            o = [row(R, x, H=0.0, d0=1e-6), row(R, x, H=0.0, d0=1e-12), row(R, x, H=0.01, d0=1e-6), row(R, x, H=0.01, d0=1e-12)]
            xa = x / R; ref = 0.75 * xa / (4 * math.sqrt(4 - 1.875 * xa))
            L.append(f"  xAR {x:6.2f}  x_A {xa:.3f}   " + " | ".join(f(q) for q in o) + f"   ref {ref:.4f}")
        gs = sorted(G + (extra24 if R == 24 else []))
        for x in gs:
            o = row(R, x, H=0.0, d0=1e-6)
            if found is None and o["lam"] > LAM_GROW and o["ok"] is not None:
                found = x; break
            prev = x
        i = gs.index(found) if found else None
        lo = gs[i - 1] if i else None
        if found:
            ol, oh = row(R, lo, H=0.0, d0=1e-6)["lam"], row(R, found, H=0.0, d0=1e-6)["lam"]
            interp = lo + (LAM_GROW - ol) / (oh - ol) * (found - lo)
        else:
            interp = float("nan")
        thr[R] = (lo, found, interp)
        L.append(f"  THRESHOLD bracket ({lo}, {found}); linear interpolation (Account) {interp:.2f}")
    L.append("\n## rulers (lambda): dr 0.05 | sponge 60 | primary")
    for R in (12, 24):
        for x in (5.5, 6.0, 6.5):
            L.append(f"  R {R} xAR {x}: {f(row(R, x, dr=0.05))} | {f(row(R, x, sw=60.0))} | {f(row(R, x))}")
    L.append("\n## edge width 0.3 vs 1.0 (R 24)")
    for x in G[3:7]:
        L.append(f"  xAR {x}: wd 0.3 {f(row(24, x, wd=0.3))} | wd 1.0 {f(row(24, x))}")
    L.append("\n## ruler T = 2000 (seed 1e-6, H 0): lambda T2000 | T600")
    for R, xs in ((24, (5.5, 6.0, 6.5, 6.72, 7.0)), (12, (5.5, 6.0, 6.5))):
        for x in xs:
            L.append(f"  R {R} xAR {x}: {f(row(R, x, T=2000.0))} | {f(row(R, x))}")
    L.append("\n## controls (must NOT grow): pump off (R24, 8.64): " + f(row(24, 8.64, pump=False)) + "; detune +10%: " + f(row(24, 8.64, detune=1.1))
             + "; xAR 4.0 seeds 1e-6 / 1e-12 (R 12, R 24): " + ", ".join(f(row(R, 4.0, d0=d)) for R in (12, 24) for d in (1e-6, 1e-12)))
    L.append("\n## PASS/FAIL against the fixed predictions")
    pred = {12: 5.8, 24: 6.0}
    for R in (12, 24):
        lo, hi, ip = thr[R]
        okb = lo is not None and hi is not None and (lo - 0.7) <= pred[R] <= (hi + 0.7)
        L.append(f"  threshold R {R}: bracket ({lo}, {hi}), prediction {pred[R]}, inside bracket widened by 0.7: {'PASS' if okb else 'FAIL'}")
    for x, p in ((6.72, 0.007), (8.64, 0.018), (12.24, 0.040)):
        o = row(24, x); okg = o["ok"] and abs(o["lam"] - p) <= 0.3 * p
        L.append(f"  growth R 24 xAR {x}: lambda {f(o)} vs {p} (30%): {'PASS' if okg else 'FAIL'}")
    ng = all(not (row(R, x)["lam"] > LAM_GROW) for R in (12, 24) for x in (4.0, 5.0))
    L.append(f"  no growth at xAR 4.0 and 5.0 (both R): {'PASS' if ng else 'FAIL'}")
    sh = []
    for R in (12, 24):
        for x in G + (extra24 if R == 24 else []):
            o0, o1 = row(R, x, H=0.0, d0=1e-6), row(R, x, H=0.01, d0=1e-6)
            if o0["lam"] > LAM_GROW and o0["ok"] and o1["ok"]:
                sh.append((R, x, o0["lam"] - o1["lam"]))
    L.append("  H-shift lambda(0) - lambda(0.01), expected 0.015 +/- 0.004: " + "; ".join(f"R{R} {x}: {d:.4f}" for R, x, d in sh)
             + f"  -> {'PASS' if sh and all(abs(d - 0.015) <= 0.004 for _, _, d in sh) else 'FAIL/mixed'}")
    sd = []
    for R in (12, 24):
        for x in G + (extra24 if R == 24 else []):
            a1, a2 = row(R, x, d0=1e-6), row(R, x, d0=1e-12)
            if a1["lam"] > LAM_GROW and a1["ok"]:
                sd.append((R, x, abs(a1["lam"] - a2["lam"]) / a1["lam"]))
    L.append("  seed independence (|dlambda|/lambda, within 10%): " + "; ".join(f"R{R} {x}: {d*100:.0f}%" for R, x, d in sd)
             + f"  -> {'PASS' if sd and all(d <= 0.10 for _, _, d in sd) else 'FAIL/mixed'}")
    ctrl = [row(24, 8.64, pump=False), row(24, 8.64, detune=1.1)] + [row(R, 4.0, d0=d) for R in (12, 24) for d in (1e-6, 1e-12)]
    L.append("  controls none grow: " + ("PASS" if all(not (c["lam"] > LAM_GROW) for c in ctrl) else "FAIL"))
    open(os.path.join(od, "summary.txt"), "w").write("\n".join(L) + "\n"); print("\n".join(L))


if __name__ == "__main__":
    main()
