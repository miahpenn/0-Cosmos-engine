"""RESEARCH ONLY. 2c of HANDOFF Nov 15 (2), the derivation owed by amendment 1: the pocket's source spectrum against
the sea's propagating and parametric thresholds. Stored states; no new runs.
Note: results/reset-handoff-2026-11-15/SEA_ANSWER.md.

Linearised about the sea (D_s at rest, d = sqrt2, s = 0; scaled):  S waves: w^2 = 4 + k^2;  D waves: w^2 = 2 + k^2.
Couplings that can move energy into sea waves (from V):
  (3/2) s^2 d  -> pumps D at the frequencies of s^2  (parametric: D waves at half the pump, needs pump >= 2 sqrt2;
                 forced: D at the pump itself, propagates if pump >= sqrt2)
  (3/2) s d^2  -> pumps S at the frequencies of d^2 (parametric: needs pump >= 4)
  d^3 / 2      -> D's own harmonics.
Mathieu, uniform pump s = A cos(w_S t) with w_S = 2: delta'' + [2 + k^2 + (3/4)x_A (1 + cos 2w_S t)] delta = 0;
first resonance 2 + k^2 + (3/4) x_A = w_S^2, growth rate mu = (3/4) x_A / (4 w_S) per unit time (amplitude).
Usage: python research/sea_answer_v1.py RUN.npz [...]
"""
import math
import sys

sys.path.insert(0, ".")
import numpy as np  # noqa: E402

from research.reset_handoff_v1 import SQ, boxcar, d_slow, TD  # noqa: E402
from research.sd_spectrum_v1 import spectrum  # noqa: E402
from scipy.signal import find_peaks  # noqa: E402

S2 = math.sqrt(2)


def lines(tau, y, n=5):
    om, P, res, *_ = spectrum(tau, y)
    pk, _ = find_peaks(P, height=P.max() * 0.03, distance=max(1, int(0.5 * res / (om[1] - om[0]))))
    pk = sorted(pk, key=lambda i: P[i], reverse=True)[:n]
    return sorted([(om[i], P[(om > om[i] - res) & (om < om[i] + res)].sum()) for i in pk]), res


def tag(w):
    out = []
    if w >= 4: out.append("pumps S parametrically")
    if w >= 2 * S2: out.append("pumps D parametrically")
    if w >= 2: out.append("an S wave")
    if w >= S2: out.append("a D wave (forced)")
    return ", ".join(out) if out else "below every threshold"


if __name__ == "__main__":
    for path in sys.argv[1:]:
        d = np.load(path, allow_pickle=True); t = d["t"]; r = d["r"]
        print(f"\n=== {path.split('/')[-1]}: x_c {float(d['x_c']):g}, r_max {r[-1]+0.05:g}")
        Ds = d_slow(d, "box2")
        for rq, lab in ((0.0, "centre (the knot)"), (3.0, "pocket, r = 3"), (8.0, "pocket edge, r = 8"), (15.0, "the sea, r = 15")):
            j = int(np.argmin(np.abs(r - rq)))
            tau = np.concatenate([[0], np.cumsum(0.5 * (d["alpha"][1:, j] + d["alpha"][:-1, j]) * np.diff(t))])
            m = (t >= 18) & (t <= 54)
            s = SQ * d["S"][m, j]; dd = SQ * d["D"][m, j]; dr_ = dd - np.nan_to_num(Ds[m, j], nan=float(np.nanmean(Ds[m, j])))
            print(f"  {lab} (t = 18-54, local proper time)")
            for name, sig in (("S", s), ("S^2 (pump on D)", s * s), ("D's ring (D - D_slow)", dr_), ("D^2 - <D^2> (pump on S)", dd * dd - np.mean(dd * dd))):
                L, res = lines(tau[m], sig)
                print(f"     {name:24s} lines (resolution +-{res:.2f}): " + "; ".join(f"{w:.2f} ({p*100:.0f}%: {tag(w)})" for w, p in L))
        # Mathieu growth against Hubble damping in the pocket core, and D's ring energy there
        print("  Mathieu (S rings D_s into D_r) in the core r < 6: growth mu = (3/32) x_A against D_r's Hubble damping (3/2) H")
        dV = d["dV"]; core = r < 6.0
        for tq in (12, 18, 24, 30, 40, 50):
            i = int(np.argmin(np.abs(t - tq)))
            s, ps = SQ * d["S"][i], SQ * d["PS"][i]
            xA = float(np.sum(((s**2 + ps**2 / 4)) * dV[i] * core) / np.sum(dV[i] * core))
            H = -d["K"][i, 0] / 3
            print(f"     t = {t[i]:4.1f}: x_A = {xA:.3f}, mu = {3/32*xA:.4f}, (3/2)H = {1.5*H:+.4f}  ->  net {3/32*xA - 1.5*H:+.4f}")
        # D_r ring energy in the sea outside the pocket (r = 10-20), comoving
        ring = SQ * d["D"] - np.nan_to_num(Ds)
        Er = np.array([np.sum((0.5 * (SQ * d["PD"][i])**2 * 0 + ring[i]**2) * dV[i] * ((r > 10) & (r < 20))) for i in range(len(t))])
        Es = np.array([np.sum((SQ * d["S"][i])**2 * dV[i] * ((r > 10) & (r < 20))) for i in range(len(t))])
        print("  the sea, r = 10-20: D's ring (sum of ring^2 dV) and S's (sum of s^2 dV), comoving; ratio")
        print("     " + "  ".join(f"t={t[i]:.0f}: {Er[i]:.2e} / {Es[i]:.2e} = {Er[i]/max(Es[i],1e-300):.3f}" for i in [int(np.argmin(np.abs(t - q))) for q in (12, 24, 36, 48)]))
