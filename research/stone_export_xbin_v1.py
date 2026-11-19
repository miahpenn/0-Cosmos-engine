"""Offline x-binning of the stone-export field run (chat's hypothesis: the absolute-omega windows slice one pair channel
at different x(t), omega_S(x) = sqrt(4 - 1.875 x) for x < 1.07; single channel 2 omega_S for 1.07 < x < 2.67).
Reruns the same 12 configs (deterministic), keeps probe series and the centre's x(t) (rows col 9, r < 1; col 17, r < 5), and bins each
window's band-passed outward flux by x at the emission time t - r/v_g (v_g = sqrt(omega^2 - 2)/omega at the window's mean omega).
Also: total w_stress of the whole radiated field (band 1.35-3.6), per run.  python research/stone_export_xbin_v1.py [outdir]"""
import math, os, sys
from multiprocessing import Pool
import numpy as np
from scipy.signal import hilbert
sys.path.insert(0, "research"); sys.path.insert(0, ".")
import stone_toy_v1 as S
import stone_export_field_v1 as X

XE = [0, 0.33, 0.6, 0.8, 1.07, 1.5, 2.67, 4.0, 1e9]


def pieces(c, j, lo, hi):
    dt, dr, pd = c["_dt"], c["dr"], c["_probe_d"]
    r0 = (round((1.5, 3.0)[j] * c["R"] / dr) + 0.5) * dr
    a, b, cc = pd[:, 3 * j], pd[:, 3 * j + 1], pd[:, 3 * j + 2]
    chi = X.band(r0 * b, dt, lo, hi); chi_t = X.band(r0 * b, dt, lo, hi, True)
    chi_r = (X.band((r0 + dr) * cc, dt, lo, hi) - X.band((r0 - dr) * a, dt, lo, hi)) / (2 * dr)
    return r0, chi, chi_t, chi_r


def job(a):
    R, gm, T, f = a
    c, rows = S.run(X.make(R, gm, T, f)); dt = c["_dt"]
    tt = rows[:, 0]; out = {"x": {}}
    # total w of the whole radiated field
    for j in (0, 1):
        r0, chi, chi_t, chi_r = pieces(c, j, 1.35, 3.6)
        rho = 0.5 * chi_t ** 2 + 0.5 * chi_r ** 2 + chi ** 2; p = 0.5 * chi_t ** 2 - chi_r ** 2 / 6 - chi ** 2
        out[("w_tot", j)] = float(p.sum() / rho.sum())
    for xname, col in (("centre", 9), ("core", 17)):
        res = {}
        for j in (0, 1):
            for nm, (lo, hi) in X.WIN.items():
                r0, chi, chi_t, chi_r = pieces(c, j, lo, hi)
                F = -chi_t * chi_r * 4 * math.pi
                an = hilbert(chi); ph = np.unwrap(np.angle(an)); om = np.gradient(ph, dt)
                Y = np.abs(np.fft.rfft(chi)) ** 2; w = 2 * math.pi * np.fft.rfftfreq(len(chi), dt)
                wm = float((w * Y).sum() / Y.sum()); vg = math.sqrt(max(wm * wm - 2, 0)) / wm
                te = np.arange(len(chi)) * dt - r0 / max(vg, 0.05)
                xe = np.interp(te, tt, rows[:, col])
                tot = max(F.sum(), 1e-300); lst = []
                for k in range(len(XE) - 1):
                    m = (xe >= XE[k]) & (xe < XE[k + 1]) & (te > 0)
                    e = float(F[m].sum()); wt = np.abs(F[m])
                    om_m = float((om[m] * wt).sum() / wt.sum()) if wt.sum() > 0 else float("nan")
                    x_m = float((xe[m] * wt).sum() / wt.sum()) if wt.sum() > 0 else float("nan")
                    lst.append((XE[k], e, om_m, x_m))
                res[(j, nm)] = (float(F.sum() * dt), lst)
        out["x"][xname] = res
    return a, out


def main():
    od = sys.argv[1] if len(sys.argv) > 1 else "results/stone-export-xbin-2026-11-19"
    os.makedirs(od, exist_ok=True)
    jobs = [(R, gm, T, f) for R, gm, T in X.CFGS for f in X.OFFS]
    with Pool(2) as p:
        res = p.map(job, jobs, chunksize=1)
    with open(os.path.join(od, "summary.txt"), "w") as fh:
        fh.write("# x-binned window flux (emission-time x at the stone centre r<1 and core r<5). Per window and x bin: share of the window's\n"
                 "# flux (E*dt, signed), flux-weighted instantaneous omega (Hilbert) and flux-weighted x; pair prediction w_S = sqrt(4 - 1.875 x)\n"
                 "# for x < 1.07; single-channel 2 w_S for 1.07 < x < 2.67. 1.5R probe unless marked.\n")
        for a, o in res:
            R, gm, T, f = a
            fh.write(f"\n## R {R} g x{gm:g} T {T:g} f {f:g}   w_tot(1.5R) {o[('w_tot',0)]:+.3f}  w_tot(3R) {o[('w_tot',1)]:+.3f}\n")
            for xn in ("centre", "core"):
                for j, pn in ((0, "1.5R"), (1, "3R")):
                    for nm in X.WIN:
                        E, lst = o["x"][xn][(j, nm)]
                        s = "  ".join(f"[x>{x0:g}] {100*e*0.08/E if E else 0:+.0f}% w{om:.2f} x{xm:.2f}" for x0, e, om, xm in lst
                                      if not math.isnan(om) and abs(e) * 0.08 > 0.04 * abs(E))
                        fh.write(f"  {xn} {pn} {nm:8s} {s}\n")
    print("done")


if __name__ == "__main__":
    main()
