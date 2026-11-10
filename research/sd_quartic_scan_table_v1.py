"""RESEARCH ONLY. Table for the GEAR-67 lambda scan (output of research/sd_quartic_scan_v1.py; usage: DIR N) (S0 = 0.02, Gaussian S width 5, D seed 1e-4 width 15, CMC, no wall)."""
import json, math, sys, glob, os
import numpy as np

F = sys.argv[1]
N = int(sys.argv[2]) if len(sys.argv) > 2 else 120


def load(lam, n=N):
    p = f"{F}/l{lam}_N{n}.json"
    return json.load(open(p)) if os.path.exists(p) else None


def arrays(d):
    h = d["hist"]
    t = np.array([x["t"] for x in h]); tau = np.array([x["tau"] for x in h])
    A = np.array([x["alpha"] for x in h]); S = np.array([x["S"] for x in h]); D = np.array([x["D"] for x in h])
    return np.array(d["r"]), t, tau, A, S, D


def omega2(t, S, A, i0, i1, j):
    """Proper-time omega^2 = -<S S_tt>/<S^2> / <alpha^2> at radius index j over [i0, i1)."""
    dt = t[1] - t[0]
    s = S[:, j]; stt = (s[2:] - 2 * s[1:-1] + s[:-2]) / dt**2
    sl = slice(max(i0, 1), min(i1, len(t) - 1))
    num = -np.sum(s[sl] * stt[sl.start - 1:sl.stop - 1]); den = np.sum(s[sl] ** 2)
    return num / den / np.mean(A[sl, j] ** 2)


ctrl = load(0)
r, tc, _, Ac, Sc, Dc = arrays(ctrl)
lams = sorted(int(os.path.basename(p)[1:].split("_")[0]) for p in glob.glob(f"{F}/l*_N{N}.json"))

print("lambda=0 control, centre lapse every 1.0:")
print("  " + "  ".join(f"t{tc[i]:.0f}:{Ac[i,0]:.3f}" for i in range(0, len(tc), 20)))
print("  |D|max every 2.0: " + "  ".join(f"t{tc[i]:.0f}:{np.abs(Dc[i]).max():.2e}" for i in range(0, len(tc), 40)))
print()
hdr = ("lam", "lamS0^2", "D0(t4)/ctrl", "gamma_c[0.5,4]", "held", "|D|max", "/sqrt(4/lam)", "min a0 (t)", "t a0<0.9", "t a0<0.5",
       "max a", "r*", "<D^2>r*", "w_meas", "w_pred(naive)", "w_pred(ctrl+quartic)")
print(" | ".join(hdr))
i4 = int(np.argmin(np.abs(tc - 4.0))); m = (tc >= 0.5) & (tc <= 4.0)
gc = np.polyfit(tc[m], np.log(np.abs(Dc[m, 0])), 1)[0]
for lam in lams:
    d = load(lam); r, t, tau, A, S, D = arrays(d)
    mm = (t >= 0.5) & (t <= 4.0)
    g = np.polyfit(t[mm], np.log(np.abs(D[mm, 0]) + 1e-300), 1)[0]
    ratio = D[int(np.argmin(np.abs(t - 4.0))), 0] / Dc[i4, 0]
    held = "yes" if g < 0.5 * gc else ("partial" if g < 0.9 * gc else "no")
    Dmax = np.abs(D).max(); cap = math.sqrt(4 / lam) if lam else float("nan")
    a0 = A[:, 0]; dtt = t[1] - t[0]
    below9 = dtt * np.sum(a0 < 0.9); below5 = dtt * np.sum(a0 < 0.5)
    # frequency of S where D is large: window where |D| max >= 0.7 of its peak; r* maximizes <D^2> over window
    peak_t = np.abs(D).max(axis=1)
    w = np.where(peak_t >= 0.7 * Dmax)[0]
    if lam and len(w) > 40:
        i0, i1 = w[0], len(t)
        D2 = np.mean(D[i0:i1] ** 2, axis=0); j = int(np.argmax(D2[r < 30]))
        i1 = min(i1, len(tc))
        wm = omega2(t, S, A, i0, i1, j); wc = omega2(tc, Sc, Ac, i0, i1, j) if i0 + 20 < len(tc) else float('nan')
        s2 = S[i0:i1, j] ** 2
        quart = 1.5 * lam * np.sum(D[i0:i1, j] ** 2 * s2) / np.sum(s2) + 0.5 * lam * np.sum(s2 * s2) / np.sum(s2)
        wn = math.sqrt(1 + 1.5 * lam * D2[j]); wp = math.sqrt(max(wc, 0) + quart) if wc == wc else float('nan')
        fr = f"{r[j]:.1f} | {D2[j]:.2e} | {math.sqrt(max(wm,0)):.3f} | {wn:.3f} | {wp:.3f}"
    else:
        fr = "- | - | - | - | -"
    print(f"{lam} | {lam*0.02**2:.2f} | {ratio:.3f} | {g:.3f} (ctrl {gc:.3f}) | {held} | {Dmax:.4f} | {Dmax/cap if lam else float('nan'):.3f} | "
          f"{a0.min():.3f} ({t[a0.argmin()]:.1f}) | {below9:.1f} | {below5:.1f} | {A.max():.3f} | {fr}")
