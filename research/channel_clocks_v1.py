"""RESEARCH ONLY. Which clock is the cosmos lane's t? (Jeremiah, Nov 14 2026: put the master clock, S's own
oscillation, to 0*.)  Stored frozen states only; no new runs.

Closures (i) and (ii) are the same slices with different time labels. For a clock c whose rate per unit label on a
slice is N_c (a lapse average), dK/d(tau_c) = <alpha Q>/N_c, and the channel per unit of that clock is
    C_c = <alpha Q>/N_c - Q(averaged fields).
N_c carries the same overall lapse scale as <alpha Q>, so C_c cannot depend on the closure for ANY clock defined on
the slice. What it does depend on is which clock. Clocks compared:
    label      N = 1                    (the closure's own label; this one is the closure-dependent one)
    wall       N = alpha at the wall
    volume     N = <alpha>, proper volume
    focusing   N = <alpha Q>/<Q>        (closure (ii)'s label)
    dust       N = <alpha> weighted by dark-matter rest mass       (matter moving with the expansion)
    S          N = <alpha> weighted by S's oscillation energy PS^2 + m^2 S^2, m^2 = 1 + 1.5 lambda D^2 (master clock)
    centre     N = alpha(0)
Also read directly: K at S's saddle crossings (crossing number as the clock), closure (i) against (ii).
Usage: python research/channel_clocks_v1.py RUN_ii.npz RUN_i.npz
"""
import sys

sys.path.insert(0, ".")
import numpy as np  # noqa: E402

from research.cosmos_channel_v1 import LAM, channel  # noqa: E402

CLOCKS = ("label", "wall", "volume", "focusing", "dust", "S", "centre")


def rates(d, i):
    dV = d["dV"][i]; al = d["alpha"][i]
    c = channel(d, i)
    wS = (d["PS"][i]**2 + (1 + 1.5 * LAM * d["D"][i]**2) * d["S"][i]**2) * dV
    wd = d["rho_dm"][i] * dV
    N = dict(label=1.0, wall=float(al[-1]), volume=float(np.sum(al * dV) / dV.sum()), focusing=c["aQ"] / c["meanQ"],
             dust=float(np.sum(al * wd) / wd.sum()), S=float(np.sum(al * wS) / wS.sum()), centre=float(al[0]))
    Qbar = c["K"]**2 / 3 + c["src"]
    return N, {k: c["aQ"] / N[k] - Qbar for k in CLOCKS}, c["fl"], c["aQ"]


def table(d, name):
    t = d["t"]; al = d["alpha"]
    tw = np.concatenate([[0], np.cumsum(0.5 * (al[1:, -1] + al[:-1, -1]) * np.diff(t))])      # wall proper time
    rows = []
    for tq in (3, 9, 15, 27, 39, 51, 63):
        sel = np.where((tw >= tq - 3) & (tw < tq + 3))[0]
        r = [rates(d, i) for i in sel]
        rows.append((tq, {k: float(np.mean([x[0][k] for x in r])) for k in CLOCKS},
                     {k: float(np.mean([x[1][k] / x[2] for x in r])) for k in CLOCKS},
                     float(np.mean([x[3] / x[2] for x in r]))))
    return rows


def crossings(d):
    ts = d["ts"]; t, s = ts[:, 0], ts[:, 1]
    i = np.where(s[:-1] * s[1:] < 0)[0]
    tc = t[i] - s[i] * (t[i + 1] - t[i]) / (s[i + 1] - s[i])
    return tc, np.interp(tc, d["t"], d["K"][:, -1]), np.interp(tc, t, ts[:, 19]), np.interp(tc, t, ts[:, 18])


if __name__ == "__main__":
    dii, di = (np.load(p, allow_pickle=True) for p in sys.argv[1:3])
    Tii, Ti = table(dii, "ii"), table(di, "i")
    print("Windows of 6 in the WALL's proper time (the same slices in both runs). C in units of the fluids' source F.")
    print("\nclock rates N per unit label, closure (ii) | closure (i)")
    print("   tau_w  " + " ".join(f"{k:>9s}" for k in CLOCKS) + "  |  " + " ".join(f"{k:>9s}" for k in CLOCKS))
    for a, b in zip(Tii, Ti):
        print(f"   {a[0]:4.0f}   " + " ".join(f"{a[1][k]:9.5f}" for k in CLOCKS) + "  |  " + " ".join(f"{b[1][k]:9.5f}" for k in CLOCKS))
    print("\nchannel per unit of each clock, C_c/F: closure (ii) | closure (i) | difference (i) - (ii)")
    for a, b in zip(Tii, Ti):
        print(f"   {a[0]:4.0f}   " + " ".join(f"{a[2][k]:+9.4f}" for k in CLOCKS) + "  |  " + " ".join(f"{b[2][k]:+9.4f}" for k in CLOCKS)
              + "  |  " + " ".join(f"{b[2][k]-a[2][k]:+9.4f}" for k in CLOCKS))
    print("\nspread of C_c/F across the six physical clocks (label left out), closure (ii); and <alpha Q>/F, the size the clock rate multiplies")
    for a in Tii:
        v = [a[2][k] for k in CLOCKS[1:]]
        print(f"   {a[0]:4.0f}   min {min(v):+.4f}  max {max(v):+.4f}  range {max(v)-min(v):.4f};   dust {a[2]['dust']:+.4f}  S {a[2]['S']:+.4f}  focusing {a[2]['focusing']:+.4f};"
              f"   <alpha Q>/F {a[3]:+.1f};   N_S - N_dust {a[1]['S']-a[1]['dust']:+.5f}")
    tc2, K2, tw2, t02 = crossings(dii); tc1, K1, tw1, t01 = crossings(di)
    n = min(len(tc1), len(tc2))
    print(f"\nS's saddle crossings as the clock: {len(tc2)} crossings in (ii), {len(tc1)} in (i)")
    print("    n    label (ii)  label (i)   K at crossing (ii)   (i)          K diff/|K0|   wall proper time (ii)  (i)")
    for j in list(range(0, min(n, 12), 2)) + list(range(12, n, 6)):
        print(f"   {j+1:3d}   {tc2[j]:8.3f}   {tc1[j]:8.3f}    {K2[j]:+.6e}  {K1[j]:+.6e}   {(K1[j]-K2[j])/abs(K2[0]):+.1e}      {tw2[j]:8.3f}  {tw1[j]:8.3f}")
    print(f"   largest |K(i) - K(ii)|/|K0| at equal crossing number, first {n} crossings: {np.max(np.abs(K1[:n]-K2[:n]))/abs(K2[0]):.1e};"
          f"  first 12: {np.max(np.abs(K1[:12]-K2[:12]))/abs(K2[0]):.1e}")
    hp = np.diff(tw2)
    for lo, hi in ((0, 13), (13, 30), (30, 55), (55, 70)):
        m = (tw2[:-1] >= lo) & (tw2[:-1] < hi)
        if m.sum() > 1:
            print(f"   S's half-period in wall proper time, crossings in {lo}-{hi}: mean {hp[m].mean():.3f}, rms scatter {hp[m].std()/hp[m].mean()*100:.1f}%  ({m.sum()} half-cycles)")
