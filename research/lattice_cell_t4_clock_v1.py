"""RESEARCH ONLY. T4 under a second ruler (HANDOFF Nov 14 2026 (4), Task 2): the same pocket in different cells,
compared in windows of the CENTRE'S PROPER TIME tau_0 = int alpha(0) dt instead of the time label, and with D's
split made on a common piece of every cell (r < 12, tapered over 10-12) instead of the whole cell.

Why: under closure (ii) the lapse is fixed by the cell average, so the central lapse differs from cell to cell at the
per-cent level (alpha(0) at t = 0: 1.070 in r_max = 12, 1.049 in r_max = 36); D grows from its seed exponentially, so a
per-cent difference of clocks is a 10-20 % difference of D's oscillation energy in a window of equal label. And the
sharp k < 0.5 filter over the whole cell has a cell-dependent set of wavenumbers, so the slow/travelling partition
moves with the cell even for identical fields.
Columns: S's action J over r < 11 (windowed sqrt(<PS^2><S^2>) dV), D at rest (r < 8, |<D>|/D*), D's slow (k < 0.5)
and travelling parts and their total (r < 8), each with the ratio to the reference (first file).
Usage: python research/lattice_cell_t4_clock_v1.py REF.npz RUN.npz ...
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, ".")
from research.sea_action_analysis_v1 import lowpass  # noqa: E402

LAM = 5000.0
DSTAR = math.sqrt(2 / LAM)
W, TR = 6.0, 12.0


def parts(p, clock):
    d = np.load(p, allow_pickle=True); r = d["r"]; dr = float(r[1] - r[0]); t = d["t"]
    al0 = d["alpha"][:, 0]
    tau = np.concatenate([[0], np.cumsum(0.5 * (al0[1:] + al0[:-1]) * np.diff(t))])
    c = tau if clock == "tau0" else t
    m = r < TR; core = (r < 8)[m]; pk = (r < 11)[m]
    w = np.where(r[m] > TR - 2, np.cos(0.25 * np.pi * (r[m] - (TR - 2)))**2, 1.0)
    out = []
    for a in np.arange(0, 36, W):
        sel = (c >= a) & (c < a + W)
        D, PD, S, PS = (d[k][sel][:, m] for k in ("D", "PD", "S", "PS"))
        v = d["dV"][sel].mean(0)[m]; grr = d["grr"][sel][:, m]
        Dm = D.mean(0); dl = D - Dm; s2 = (S**2).mean(0); m2 = np.clip(-1 + 1.5 * LAM * (Dm**2 + s2), 0, None)
        A, At = lowpass(dl * w, r[m], True, 0.5, dr), lowpass(PD * w, r[m], True, 0.5, dr); B, Bt = dl * w - A, PD * w - At
        g = lambda f: np.gradient(f, dr, axis=-1) / np.sqrt(grr)  # noqa: E731
        es = ((0.5 * At**2 + 0.5 * g(A)**2 + 0.5 * m2 * A**2).mean(0) * v)[core].sum()
        eb = ((0.5 * Bt**2 + 0.5 * g(B)**2 + 0.5 * m2 * B**2).mean(0) * v)[core].sum()
        J = (np.sqrt((PS**2).mean(0) * s2) * v)[pk].sum()
        rest = float(np.sum((np.abs(Dm) * v)[core]) / np.sum(v[core])) / DSTAR
        out.append((J, rest, es, eb, es + eb))
    return np.array(out), float(tau[-1]), float(al0[0])


if __name__ == "__main__":
    for clock in ("t", "tau0"):
        ref, _, _ = parts(sys.argv[1], clock)
        print(f"\nwindows of {W:.0f} in {'the time label t' if clock == 't' else 'the centre proper time tau_0'}; D split on r < {TR:.0f}; ratios to {os.path.basename(sys.argv[1])[:-4]}")
        for p in sys.argv[1:]:
            q, tend, a0 = parts(p, clock)
            print(f"  {os.path.basename(p)[:-4]}: alpha(0) at t = 0: {a0:.4f}; tau_0 at the end {tend:.2f}")
            for i in range(len(q)):
                print(f"     {W*i:2.0f}-{W*i+W:<2.0f}  J/J0 {q[i,0]/q[0,0]:.4f} [{q[i,0]/q[0,0]/(ref[i,0]/ref[0,0]):.4f}]  D at rest {q[i,1]:.4f} [{q[i,1]-ref[i,1]:+.4f}]"
                      f"  slow {q[i,2]:.3e} [{q[i,2]/ref[i,2]:.3f}]  travelling {q[i,3]:.3e} [{q[i,3]/ref[i,3]:.3f}]  total [{q[i,4]/ref[i,4]:.3f}]")
