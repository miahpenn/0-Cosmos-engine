"""RESEARCH ONLY. The central dark-matter spike, from existing output only (HANDOFF Nov 14 2026 (3), Task E).
No new runs. All inputs are expanding one-budget runs of research/dm_one_budget_v1.py and are PROVISIONAL (AUDIT §C:
the outer boundary does not fit an expanding background).
Per run: rho_dm at the innermost cell against r = 2 (so the pocket's general gathering is divided out) over time; the
spike's half-width (radius where the excess over r = 2 falls to half its central value); S's core size (radius where
the window mean of S^2, last 8 time units, falls to half its central value); <D^2>/D*^2 at the centre.
Usage: python research/spike_from_disk_v1.py LABEL=RUN.npz ...
"""
import sys

import numpy as np

LAM = 5000.0

if __name__ == "__main__":
    print("PROVISIONAL (expanding runs, production outer edge). Spike = rho_dm(innermost cell) / rho_dm(r = 2).")
    for arg in sys.argv[1:]:
        lab, p = arg.split("=")
        d = np.load(p); r, t = d["r"], d["t"]; dr = r[1] - r[0]; j2 = int(np.argmin(np.abs(r - 2)))
        row = []
        for tt in (24, 28, 32, 36, 40):
            if tt > t[-1] + 0.3:
                row.append("   -  "); continue
            i = int(np.argmin(np.abs(t - tt)))
            row.append(f"{d['rho_dm'][i, 0] / d['rho_dm'][i, j2]:6.2f}")
        i = len(t) - 1
        ex = d["rho_dm"][i] / d["rho_dm"][i, j2] - 1
        hw = r[np.argmax(ex < 0.5 * ex[0])] if ex[0] > 0.05 else float("nan")
        sel = t >= t[-1] - 8
        s2 = (d["S"][sel] ** 2).mean(0)
        core = r[np.argmax(s2 < 0.5 * s2[0])] if s2[0] > 0 else float("nan")
        d2 = (d["D"][sel, 0] ** 2).mean() * LAM / 2
        print(f"  {lab:28s} t_end {t[-1]:5.1f} | spike at t = 24/28/32/36/40: " + " ".join(row)
              + f" | spike half-width at t_end {hw:.2f} ({hw/dr:.1f} cells) | S core (half of <S^2>) r = {core:.2f} | <D^2>(0)/D*^2 = {d2:.2f}")
