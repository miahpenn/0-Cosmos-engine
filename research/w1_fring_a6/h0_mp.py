"""W1 starting-state H wiggle, part B (python3 + mpmath, no engine import).

From the slices h0_build.py stored (t = 0, bitwise the archived replays' start), per band shell and native cell:
  H64      the engine's binary64 H;
  Hmp      H in 60 digits on the same stored binary64 a, b, X, Aa, K, Lambda, rho (removes evaluation rounding);
  Hrebuilt H in 60 digits on a, b, X, Lambda rebuilt in 60 digits from the areal loop's last binary64 source
           (exact eight-node weights, exact cumulative sum, exact powers; Aa, K, rho kept as stored, undifferentiated
           in H). This removes the rounding the build put into the stored geometry.
A wiggle that survives Hrebuilt is in the source (the supplied S, PS, D, PD or the discretisation), not rounding.

Usage: python3 -m research.w1_fring_a6.h0_mp SLICES.npz OUT.json
"""
import json
import sys
from fractions import Fraction

import mpmath as mp
import numpy as np

from research.control_tests import mp60 as M

BAND = (0.05, 0.25)


def eight_node_weights():
    """Exact v55_initial._eight_node_origin_integral weights (degree-7 interpolant on each cell interval)."""
    def solve(offset, left):
        nodes = [Fraction(j) - offset + Fraction(1, 2) for j in range(8)]
        mom = [(Fraction(1, 2) ** p - Fraction(left) ** p) / p for p in range(1, 9)]
        A = [[nodes[j] ** i for j in range(8)] for i in range(8)]    # vander(nodes, increasing).T
        Ab = [row[:] + [mom[i]] for i, row in enumerate(A)]
        for c in range(8):
            piv = next(i for i in range(c, 8) if Ab[i][c] != 0)
            Ab[c], Ab[piv] = Ab[piv], Ab[c]
            for i in range(8):
                if i != c and Ab[i][c] != 0:
                    f = Ab[i][c] / Ab[c][c]
                    Ab[i] = [x - f * y for x, y in zip(Ab[i], Ab[c])]
        return [mp.mpf(Ab[i][8].numerator) / Ab[i][8].denominator / (mp.mpf(Ab[i][i].numerator) / Ab[i][i].denominator)
                for i in range(8)]
    return solve


def rebuild(src, r, dr, g):
    n = len(src)
    solve = eight_node_weights()
    dr_ = mp.mpf(float(dr))
    samples = [mp.mpf(float(x)) - 1 for x in src]   # f - 1.0 is exact in binary64 here (Sterbenz)
    cache = {}
    integ, acc = [], mp.mpf(0)
    for i in range(n - 4):                           # windows need i + 3 < n
        start = min(i - 4, n - 8)
        off = i - start
        key = (off, 0 if i == 0 else 1)
        if key not in cache:
            cache[key] = solve(off, 0 if i == 0 else Fraction(-1, 2))
        idx = [start + j for j in range(8)]
        vals = [samples[-1 - k] if k < 0 else samples[k] for k in idx]
        acc += dr_ * mp.fsum(w * v for w, v in zip(cache[key], vals))
        integ.append(r[i].v + acc)
    B = [x / r[i].v for i, x in enumerate(integ)]
    A = [1 / mp.sqrt(x) for x in B]
    a = [mp.power(x, mp.mpf(4) / 3) for x in A]
    b = [mp.power(x, mp.mpf(-2) / 3) for x in A]
    X = [mp.power(x, mp.mpf(-1) / 3) for x in A]
    m = len(a)
    gg = M.Grid(m, dr)
    Ea, Eb = [M.E(x) for x in a], [M.E(x) for x in b]
    ap, bp = gg.d(Ea, 1, 1), gg.d(Eb, 1, 1)
    Lam = [ap[i].v / (2 * a[i] * a[i]) - bp[i].v / (a[i] * b[i]) + 2 / gg.r[i].v * (1 / b[i] - 1 / a[i]) for i in range(m)]
    return {"a": a, "b": b, "X": X, "Lambda": Lam}


def H_of(g, s, rho):
    s = dict(s, alpha=M.const(1, g.n))   # geometry_terms also forms lapse terms; H does not use them
    Hv = M.hamiltonian_vac(g, s)
    return [h.v - 16 * mp.pi * mp.mpf(float(p)) for h, p in zip(Hv, rho)]


def band_read(v, r, ok):
    w = (r > BAND[0]) & (r < BAND[1]) & ok
    b = v[w]
    s = np.sign(b)
    edges = np.arange(BAND[0], BAND[1] + 1e-12, 0.05)
    sh = [float(np.sum(r[(r > lo) & (r < hi)] ** 2 * v[(r > lo) & (r < hi)]) / np.sum(r[(r > lo) & (r < hi)] ** 2))
          for lo, hi in zip(edges[:-1], edges[1:])]
    return {"cells": b.tolist(), "sign_changes": int(np.sum(s[1:] != s[:-1])), "pairs": int(len(b) - 1),
            "signs": "".join("+" if x > 0 else "-" for x in s), "shell_means": sh,
            "band_rms": float(np.sqrt(np.mean(b * b))),
            "second_difference_rms": float(np.sqrt(np.mean(np.diff(b, 2) ** 2))) if len(b) > 2 else None}


def main(slices, out):
    res = {}
    with np.load(slices) as z:
        for dr in (0.05, 0.025, 0.0125):
            t = f"dr{dr}_"
            r = np.asarray(z[t + "r"], dtype=float)
            n = len(r)
            g = M.Grid(n, dr)
            s = {k: M.arr(z[t + k]) for k in ("a", "b", "X", "Aa", "K", "Lambda")}
            H64 = np.asarray(z[t + "H"], dtype=float)
            Hmp = np.array([float(x) for x in H_of(g, s, z[t + "rho"])])
            rb = rebuild(z[t + "src_last"], g.r, dr, g)
            m = len(rb["Lambda"])
            g2 = M.Grid(m, dr)
            s2 = {k: [M.E(x) for x in rb[k]] for k in ("a", "b", "X", "Lambda")}
            s2.update({k: M.arr(z[t + k][:m]) for k in ("Aa", "K")})
            Hrb = np.full(n, np.nan)
            Hrb[:m] = [float(x) for x in H_of(g2, s2, z[t + "rho"][:m])]
            ok = np.arange(n) < m - 6                        # rebuilt-state stencil reach (edge rows of the cut)
            geo = {k: float(np.max(np.abs(np.asarray(z[t + k][:m - 6], dtype=float)
                                          - np.array([float(x) for x in rb[k][:m - 6]]))))
                   for k in ("a", "b", "X", "Lambda")}
            res[str(dr)] = {"H64": band_read(H64, r, ok), "Hmp": band_read(Hmp, r, ok),
                            "Hrebuilt": band_read(Hrb, r, ok),
                            "H64_minus_Hmp_band_rms": float(np.sqrt(np.mean((H64 - Hmp)[(r > .05) & (r < .25) & ok] ** 2))),
                            "stored_minus_rebuilt_geometry_max": geo,
                            "16pi_rho_band": float(16 * np.pi * np.max(np.abs(z[t + "rho"][(r > .05) & (r < .25)])))}
    json.dump(res, open(out, "w"), indent=1)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
