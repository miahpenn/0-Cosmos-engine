"""N vs 2N for the stone-export field run: R 24 g f 0, R 24 g f 1, R 12 8g f 0 at dr 0.2 and 0.1 (design in the coordination log).
python research/stone_export_n2n_v1.py [outdir]"""
import math, os, sys
from multiprocessing import Pool
import numpy as np
sys.path.insert(0, "research"); sys.path.insert(0, ".")
import stone_toy_v1 as S
import stone_export_field_v1 as X
import stone_export_xbin_v1 as B

CASES = [(24, 1.0, 500.0, 0.0), (24, 1.0, 500.0, 1.0), (12, 8.0, 400.0, 0.0)]
DRS = (0.2, 0.1)
if os.environ.get("N2N_FINE"):   # dr 0.1 vs 0.05 on R 24 only (design in the log)
    CASES = CASES[:2]; DRS = (0.1, 0.05)
    if os.environ.get("N2N_FINE") == "R12":
        CASES = [(12, 8.0, 400.0, 0.0), (12, 8.0, 400.0, 1.0)]


def make(R, gm, T, f, dr):
    c = X.make(R, gm, T, f); c["dr"] = dr; pr = []
    for fac in (1.5, 3.0):
        i = round(fac * R / dr); pr += [(i - 1) * dr, i * dr, (i + 1) * dr]
    c["probe"] = tuple(pr); return c


def job(a):
    (R, gm, T, f), dr = a
    c, rows = S.run(make(R, gm, T, f, dr)); o = X.analyse(c, rows)
    for j in (0, 1):
        r0, chi, chi_t, chi_r = B.pieces(c, j, 1.35, 3.6)
        rho = 0.5 * chi_t ** 2 + 0.5 * chi_r ** 2 + chi ** 2; p = 0.5 * chi_t ** 2 - chi_r ** 2 / 6 - chi ** 2
        o[("w", j)] = float(p.sum() / rho.sum())
    return a, o


def main():
    od = sys.argv[1] if len(sys.argv) > 1 else "results/stone-export-n2n-2026-11-19"
    os.makedirs(od, exist_ok=True)
    jobs = [(cs, dr) for cs in CASES for dr in DRS]
    with Pool(2) as p:
        res = dict(p.map(job, jobs, chunksize=1))
    L = []; verdict = []
    for cs in CASES:
        o2, o1 = res[(cs, DRS[0])], res[(cs, DRS[1])]; ok = True
        L.append(f"## R {cs[0]} g x{cs[1]:g} T {cs[2]:g} f {cs[3]:g}   E3 dr{DRS[0]} {o2["E3"]*100:+.1f}% dr{DRS[1]} {o1["E3"]*100:+.1f}%   drift {o2['drift']:.1e} / {o1['drift']:.1e}"
                 f"   sponge {o2['sponge']*100:.0f}% / {o1['sponge']*100:.0f}%")
        for j, pn in ((0, 1.5), (1, 3.0)):
            L.append(f"  probe {pn}R   w_tot {o2[('w', j)]:+.3f} / {o1[('w', j)]:+.3f}  (d {abs(o2[('w', j)] - o1[('w', j)]):.3f})")
            ok &= abs(o2[("w", j)] - o1[("w", j)]) <= 0.02
            for nm in X.WIN:
                a2, a1 = o2[pn][nm], o1[pn][nm]; ds = abs(a2["share"] - a1["share"]) * 100; dp = abs(a2["pk"] - a1["pk"])
                ok &= (ds <= 5 and dp <= 0.05)
                L.append(f"    {nm:8s} share {a2['share']*100:5.1f}% / {a1['share']*100:5.1f}% (d {ds:4.1f} pp)  peak {a2['pk']:.3f} / {a1['pk']:.3f} (d {dp:.3f})"
                         f"  w_s {a2['w_s']:+.3f} / {a1['w_s']:+.3f}")
        verdict.append((cs, ok))
    L.append("\nVERDICT (pass: shares within 5 pp, peaks within 0.05, total w within 0.02, at both probes): "
             + "; ".join(f"R{cs[0]} g x{cs[1]:g} f{cs[3]:g}: {'PASS' if ok else 'FAIL'}" for cs, ok in verdict))
    open(os.path.join(od, "summary.txt"), "w").write("\n".join(L) + "\n"); print("\n".join(L))


if __name__ == "__main__":
    main()
