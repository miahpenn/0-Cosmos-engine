"""RESEARCH ONLY. The reset-to-lane handoff (Nov 15 2026, next step after PR #488). Stored states; no new runs.
Note: results/reset-handoff-2026-11-15/README.md.

1. Separation rule. D = D_slow + D_ring. D_ring's lines: D's own ring at sqrt2 (period 4.443) and the forced ring at
   twice S's clock (4, period pi/2), both in local proper time. D_slow = D averaged by a boxcar of one D period
   convolved with a boxcar of one forced period, each stretched by 1/alpha (label time). Rulers: the same without the
   second boxcar; a Hann window of two D periods; no lapse stretching.
2. Components (scaled): held step h(D_slow) = (1/2)(1 - d_slow^2/2)^2, w = -1; ringing S-D matter = the rest of the
   S-D energy above D's vacuum.
3. Signed accounting, per proper-volume element (CMC, zero shift: d ln dV/dt = -alpha K):
   Q = -d(E_h)/dt - rho_h * alpha K dV   summed over the cell  (held -> ringing; one-way if Q >= 0)
4. Handoff time ruler: the lane (S-D as dust) started at t_h = 12, 18, 24 must agree afterwards.
Usage: python research/reset_handoff_v1.py RUN.npz [...]
"""
import math
import sys

sys.path.insert(0, ".")
import numpy as np  # noqa: E402

from research.b6_frozen_test_v1 import KAP, LAM  # noqa: E402

SQ = math.sqrt(LAM)
TD = 2 * math.pi / math.sqrt(2)          # D's ring period at rest
T4 = math.pi / 2                         # the forced ring at twice S's clock (2 x 2)


def boxcar(t, y, L):
    """Exact running mean over [t - L/2, t + L/2] of the linear interpolant; NaN where the window leaves the record."""
    out = np.full_like(y, np.nan, dtype=float)
    v = ~np.isnan(y)
    tv, yv = t[v], y[v]
    c = np.concatenate([[0], np.cumsum(0.5 * (yv[1:] + yv[:-1]) * np.diff(tv))])
    lo, hi = tv - L / 2, tv + L / 2
    o = (np.interp(hi, tv, c) - np.interp(lo, tv, c)) / L
    o[(lo < tv[0]) | (hi > tv[-1])] = np.nan
    out[v] = o
    return out


def hann(t, y, L):
    out = np.full_like(y, np.nan)
    for i, tc in enumerate(t):
        if tc - L / 2 < t[0] or tc + L / 2 > t[-1]:
            continue
        m = np.abs(t - tc) <= L / 2
        w = np.cos(np.pi * (t[m] - tc) / L)**2
        out[i] = np.sum(w * y[m]) / np.sum(w)
    return out


def d_slow(d, rule):
    t = d["t"]; D = SQ * d["D"]; al = d["alpha"]
    out = np.empty_like(D)
    for j in range(D.shape[1]):
        a = 1.0 if rule == "no_lapse" else float(np.mean(al[:, j]))
        if rule in ("box2", "no_lapse"):
            out[:, j] = boxcar(t, boxcar(t, D[:, j], TD / a), T4 / a)
        elif rule == "box1":
            out[:, j] = boxcar(t, D[:, j], TD / a)
        elif rule == "hann2":
            out[:, j] = hann(t, D[:, j], 2 * TD / a)
    return out


def held(d, Ds):
    h = 0.5 * (1 - 0.5 * Ds**2)**2
    return h


def sd_energy(d, i):
    """S-D energy density above D's vacuum (scaled), per cell."""
    dr = float(d["r"][1] - d["r"][0]); grr = d["grr"][i]
    s, dd, ps, pd = (SQ * d[k][i] for k in ("S", "D", "PS", "PD"))
    g2 = (np.gradient(s, dr)**2 + np.gradient(dd, dr)**2) / grr
    V = 0.5 * s**2 - 0.5 * dd**2 + (s**4 + 6 * s**2 * dd**2 + dd**4) / 8 + 0.5
    return 0.5 * (ps**2 + pd**2) + 0.5 * g2 + V


if __name__ == "__main__":
    from research.cosmos_lane_follow_v1 import averages, lane
    from research.b6_frozen_test_v1 import P_FULL, v_sd
    for path in sys.argv[1:]:
        d = np.load(path, allow_pickle=True); t = d["t"]; dV = d["dV"]
        print(f"\n=== {path.split('/')[-1]}: x_c {float(d['x_c']):g}, r_max {d['r'][-1]+0.05:g}, closure {d['closure']}")
        rules = ("box2", "box1", "hann2", "no_lapse")
        H = {r: held(d, d_slow(d, r)) for r in rules}
        E = np.array([sd_energy(d, i) for i in range(len(t))])
        mean = lambda q, i: float(np.sum(q[i] * dV[i]) / np.sum(dV[i]))                  # noqa: E731
        print("  1-2. held step / S-D energy, by separation rule (box2 = the derived rule)   | held energy, comoving (x a^3, /S-D energy at t = 4.4)")
        print("       t     box2      box1      hann2     no_lapse   spread(rulers)/box2    |  box2")
        i0 = int(np.argmin(np.abs(t - 4.5)))
        a3 = lambda i: float(np.sum(dV[i]) / np.sum(dV[i0]))                              # noqa: E731
        for tq in (4.5, 5, 6, 7, 8, 9, 10, 12, 15, 18, 24, 30, 40, 54, 62):
            i = int(np.argmin(np.abs(t - tq)))
            v = [mean(H[r], i) for r in rules]
            if any(np.isnan(v)):
                continue
            e = mean(E, i)
            spread = (max(v) - min(v)) / max(v[0], 1e-300)
            print(f"     {t[i]:5.1f}   {v[0]/e:.5f}   {v[1]/e:.5f}   {v[2]/e:.5f}   {v[3]/e:.5f}    {spread:8.2f}                |  {v[0]*a3(i)/(mean(E,i0)):.5f}")
        # 3. signed accounting with the derived rule
        Hb = H["box2"]; K = d["K"][:, 0]
        Eh = np.array([np.nansum(Hb[i] * dV[i]) for i in range(len(t))])
        work = np.array([np.nansum(Hb[i] * d["alpha"][i] * K[i] * dV[i]) for i in range(len(t))])     # rho_h alpha K dV
        Q = -np.gradient(Eh, t) - work
        Etot = np.array([np.sum(E[i] * dV[i]) for i in range(len(t))])
        ok = ~np.isnan(Eh)
        print("  3. signed transfer held -> ringing, Q (per unit label time, cell total, / S-D energy):")
        line = []
        for tq in (5, 6, 7, 8, 9, 10, 12, 15, 20, 30, 45, 60):
            i = int(np.argmin(np.abs(t - tq)))
            if ok[i]:
                line.append(f"t={t[i]:.0f}: {Q[i]/Etot[i]:+.2e}")
        print("     " + "  ".join(line))
        m = ok & (t > 4.4) & (t < 15)
        print(f"     share of time with Q < 0 (held gaining from ringing), t = 4.4-15: {np.mean(Q[m] < 0):.2f};  t = 15-60: {np.mean(Q[ok & (t > 15) & (t < 60)] < 0):.2f}")
        Ir = np.trapezoid(Q[m], t[m]); print(f"     integrated Q over t = 4.4-15: {Ir/Eh[i0]:.3f} of the held energy at t = 4.4; held energy then left at t = 15: {Eh[int(np.argmin(np.abs(t-15)))]/Eh[i0]:.3f}")
        # 4. handoff-time ruler
        A, tau = averages(d)
        P = P_FULL; rho_vac = KAP * v_sd(0.0, math.sqrt(2 / LAM))
        print("  4. handoff-time ruler: lane (S-D as dust) from t_h; held step at t_h by the derived rule; residual against the cell")
        print("       t_h   held/S-D at t_h   at t:  dH/H(t=0)             d ln a              d rho_SD/rho_SD")
        res = {}
        for th in (8.0, 12.0, 18.0, 24.0):
            i0h = int(np.argmin(np.abs(t - th)))
            y0 = [A[i0h, 10], A[i0h, 2], A[i0h, 3], A[i0h, 4], A[i0h, 5], A[i0h, 6], A[i0h, 7], A[i0h, 8]]
            Y = lane(y0, tau[i0h:], P, rho_vac, 0.0, "raychaudhuri")
            row = []
            for tq in (30.0, 40.0, 54.0, 66.0):
                i = int(np.argmin(np.abs(t - tq))); k = i - i0h
                row.append(((Y[k, 1] - A[i, 2]) / abs(A[0, 2]), Y[k, 0] - A[i, 10], Y[k, 7] / A[i, 8] - 1))
            res[th] = row
            hf = mean(H["box2"], i0h) / mean(E, i0h)
            print(f"      {th:4.0f}    {hf:.4f}        " + "   ".join(f"t={tq:.0f}: {a:+.1e} {b:+.1e} {c:+.1e}" for tq, (a, b, c) in zip((30, 40, 54, 66), row)))
        for th in (8.0, 18.0, 24.0):
            dmax = max(abs(res[th][k][2] - res[12.0][k][2]) for k in range(4))
            print(f"      rho_SD at t = 30-66, t_h = {th:g} against t_h = 12: largest difference {dmax:.2e}")
