"""Stone-export field run: which channel carries energy out of a stone, and is it stable under the start offset?
Design fixed beforehand in the coordination log (2026-11-19, 'STONE-EXPORT FIELD RUN: FINAL DESIGN').

  python research/stone_export_field_v1.py [outdir]

Configs: R 24 at g (T 500), R 18 at g (T 400), R 12 at 8g (T 400); phi = 0, dr 0.2, gravity on.
Start offset f: f = 1 is #516's start (SP profile, D at the sea value + 1e-6); f < 1 is the quiet start (x(0) = 0.3, SP shape)
with D at d^2 = 2 - 1.5 (1 - f) x(r) (f = 0: D at its adiabatic rest value; Choice: f scales the D mismatch linearly).
Probes at 1.5R and 3R (cell and both neighbours). Use chi = r (d - mean d): a plane wave, so w = p/rho = k^2/(3 omega^2) exactly
(D's fluctuation mass^2 = 2).  Band-passed (rectangular FFT mask) in the log's windows; flux = -chi_t chi_r (x 4 pi) outward.
"""
import math, os, sys
from multiprocessing import Pool
import numpy as np
sys.path.insert(0, "research"); sys.path.insert(0, ".")
import stone_toy_v1 as S

WIN = {"near-gap": (1.35, 1.55), "mid": (1.55, 1.9), "broad": (2.0, 2.8), "high": (3.4, 3.6)}
M2 = 2.0
CFGS = [(24, 1.0, 500.0), (18, 1.0, 400.0), (12, 8.0, 400.0)]
OFFS = (0.0, 0.1, 0.3, 1.0)


def make(R, gm, T, f):
    base = dict(R=R, g=gm * S.G_PROD, dr=0.2, grav=True, T=T, rmax=8.0 * R)
    pr = []
    for fac in (1.5, 3.0):
        c = round(fac * R / 0.2)
        pr += [(c - 1) * 0.2, c * 0.2, (c + 1) * 0.2]      # run() records cell round(pr/dr): exactly c-1, c, c+1
    base["probe"] = tuple(pr)
    if f >= 1.0:
        return dict(base, start="sp")
    return dict(base, start="quiet", shape="sp", xc=S.XC_QUIET, seed="none", dcoef=1.5 * (1 - f))


def band(y, dt, lo, hi, deriv=False):
    n = len(y); F = np.fft.rfft(y - y.mean()); om = 2 * math.pi * np.fft.rfftfreq(n, dt)
    F = np.where((om >= lo) & (om <= hi), F, 0.0)
    if deriv:
        F = 1j * om * F
    return np.fft.irfft(F, n)


def analyse(c, rows):
    dt = c["_dt"]; pd = c["_probe_d"]; dr = c["dr"]; out = {}
    for j, fac in enumerate((1.5, 3.0)):
        a, b, cc = pd[:, 3 * j], pd[:, 3 * j + 1], pd[:, 3 * j + 2]
        r0 = (round(fac * c["R"] / dr) + 0.5) * dr
        d0 = b - b[:20].mean() if False else b
        res = {}
        tot = 0.0
        for nm, (lo, hi) in WIN.items():
            chi_t = band(r0 * b, dt, lo, hi, True)
            chi_p = band((r0 + dr) * cc, dt, lo, hi); chi_m = band((r0 - dr) * a, dt, lo, hi)
            chi_r = (chi_p - chi_m) / (2 * dr); chi = band(r0 * b, dt, lo, hi)
            F = -chi_t * chi_r * 4 * math.pi
            rho = 0.5 * chi_t ** 2 + 0.5 * chi_r ** 2 + 0.5 * M2 * chi ** 2
            p = 0.5 * chi_t ** 2 - chi_r ** 2 / 6.0 - 0.5 * M2 * chi ** 2
            Y = np.abs(np.fft.rfft(chi)) ** 2; om = 2 * math.pi * np.fft.rfftfreq(len(chi), dt)
            wm = float((om * Y).sum() / max(Y.sum(), 1e-300)); pk = float(om[np.argmax(Y)])
            E = float(F.sum() * dt); w_s = float(p.sum() / max(rho.sum(), 1e-300))
            W = len(chi) * dt; dom = 2 * math.pi / W
            w_d = (wm * wm - M2) / (3 * wm * wm); lo_, hi_ = [(x * x - M2) / (3 * x * x) for x in (wm - dom, wm + dom)]
            vg = math.sqrt(max(3 * w_s, 0.0))
            tr = np.arange(len(chi)) * dt - r0 / max(vg, 0.05)                 # retarded time
            bins = [float(F[(tr >= t0) & (tr < t0 + 40)].sum() * dt) for t0 in range(0, int(c["T"]), 40)]
            res[nm] = dict(E=E, pk=pk, wm=wm, w_s=w_s, w_d=w_d, w_dtol=(min(lo_, hi_), max(lo_, hi_)), vg=vg, bins=bins)
            tot += max(E, 0.0)
        for nm in res:
            res[nm]["share"] = max(res[nm]["E"], 0.0) / max(tot, 1e-300)
        out[fac] = res
    t = rows[:, 0]
    out["E3"] = float(rows[-1, 2] / rows[0, 2] - 1); out["drift"] = float(np.abs((rows[:, 18] + rows[:, 19]) / rows[0, 18] - 1).max())
    out["sponge"] = float(rows[-1, 19] / rows[0, 18])
    return out


def job(a):
    R, gm, T, f = a
    c, rows = S.run(make(R, gm, T, f))
    return a, analyse(c, rows)


def main():
    od = sys.argv[1] if len(sys.argv) > 1 else "results/stone-export-field-2026-11-19"
    os.makedirs(od, exist_ok=True)
    jobs = [(R, gm, T, f) for R, gm, T in CFGS for f in OFFS]
    with Pool(2) as p:
        res = p.map(job, jobs, chunksize=1)
    with open(os.path.join(od, "summary.txt"), "w") as fh:
        fh.write("# stone-export field run; windows omega: " + str(WIN) + "\n# E: time-integrated band-passed outward energy "
                 "through the probe (D only); share = positive E / sum positive E; w_s = int p / int rho (stress ruler);\n"
                 "# w_d: w from the window's power-weighted omega, [lo, hi] over omega +/- 2 pi / W; vg = sqrt(3 w_s); "
                 "peak = FFT peak omega in the window.\n")
        for a, o in res:
            R, gm, T, f = a
            fh.write(f"\n## R {R} g x{gm:g} T {T:g} f {f:g}   E3 change {o['E3']*100:+.1f}%  sponge take {o['sponge']*100:.1f}% of E0  drift {o['drift']:.1e}\n")
            for fac in (1.5, 3.0):
                fh.write(f"  probe {fac}R\n")
                for nm, v in o[fac].items():
                    fh.write(f"    {nm:8s} E {v['E']:10.4g} share {v['share']*100:5.1f}%  peak {v['pk']:.3f} mean {v['wm']:.3f}  "
                             f"w_s {v['w_s']:+.3f}  w_d {v['w_d']:+.3f} [{v['w_dtol'][0]:+.3f},{v['w_dtol'][1]:+.3f}]  vg {v['vg']:.2f}\n")
    print(open(os.path.join(od, "summary.txt")).read())


if __name__ == "__main__":
    main()
