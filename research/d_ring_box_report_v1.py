"""RESEARCH ONLY. Tables for the GEAR-105 field check (research/d_ring_box_v1.py output; HANDOFF Nov 14 2026 (3), Task A).
Usage: python research/d_ring_box_report_v1.py DIR
"""
import glob
import os
import sys

import numpy as np

sys.path.insert(0, ".")
from research.d_ring_box_v1 import DSTAR, amp_S, rates  # noqa: E402

if __name__ == "__main__":
    d0 = sys.argv[1]
    print("(a), (b) linear phase. Floquet (fine grid, research/d_ring_floquet_fine_v1.py): rate/eps^2 = 1.114, 1.167, 1.283;"
          " peak k^2/eps^2 = 3.10, 3.20, 3.15; band top k^2/eps^2 = 6.25, 6.40, 6.40 (eps = 0.05, 0.1, 0.2).")
    print("  run                              | D waves: peak rate   rate/eps^2   at k (k^2/eps^2)   band (rate > 10% of peak) | S waves: largest fitted rate   largest amplification over the fit window")
    for f in sorted(glob.glob(os.path.join(d0, "lin_*.npz")) + glob.glob(os.path.join(d0, "ctrl_*.npz"))):
        r = dict(np.load(f, allow_pickle=True)); eps = float(r["eps"]); kind = str(r["kind"])
        k, rS, rD, (t1, t2) = rates(r)
        rDm = np.where(np.isnan(rD), -np.inf, rD); rSm = np.where(np.isnan(rS), -np.inf, rS)
        j = int(np.argmax(rDm)); band = k[1:][rDm[1:] > 0.1 * rDm[j]]
        dtxt = (f"{rDm[j]:.3e}    {rDm[j]/eps**2:6.3f}     {k[j]:.4f} ({k[j]**2/eps**2:.2f})     k = {band.min():.4f}-{band.max():.4f}"
                if kind in ("D",) and eps > 0 else (f"{rDm[j]:.1e} (no ring: control)" if eps == 0 else "(S seed only: D unseeded)          "))
        stxt = f"{rSm.max():+.1e}        {amp_S(r, t1, t2):.3f}" if kind in ("S", "both") else "(D seed only: S = 0 exactly)"
        print(f"  {os.path.basename(f)[:-4]:32s} | {dtxt:75s} | {stxt}")
    print("\n(c) nonlinear end state (seed 1e-6 D* in both S and D). Shares of the total energy above the vacuum:")
    print("  run                      t/T:   uniform ring | long D waves (0<k<0.5) | short D waves (k>=0.5) | S | (sum; the rest is interaction/anharmonic)   at 0, 25, 50, 75, 100% of the run; last column: mean over the last 20%; energy drift")
    for f in sorted(glob.glob(os.path.join(d0, "end_*.npz"))):
        r = np.load(f, allow_pickle=True); E = r["E"]; t = r["t"]
        tot = E[:, 4]
        def sh(i):
            return " ".join(f"{E[i, c] / tot[i]:.3f}" for c in range(4))
        idx = [0, len(t) // 4, len(t) // 2, 3 * len(t) // 4, len(t) - 1]
        late = t >= 0.8 * t[-1]
        lm = " ".join(f"{(E[late, c] / tot[late]).mean():.3f}" for c in range(4))
        # mean wavenumber of the D waves at late times (amplitude-weighted over the recorded modes)
        k = r["k"]; w = ((r["Dtk"][late] ** 2 + (k**2 + 2) * r["Dk"][late] ** 2).mean(0))[1:]
        kbar = float(np.sum(k[1:] * w) / np.sum(w))
        print(f"  {os.path.basename(f)[:-4]:22s} T = {t[-1]:6.0f}: " + " | ".join(sh(i) for i in idx) + f" || late mean: {lm}; mean k of D waves {kbar:.3f}"
              f" ({kbar/float(r['eps']):.2f} eps); drift {(tot[-1]-tot[0])/tot[0]:+.1e}")
