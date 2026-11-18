"""RESEARCH ONLY. Bound stones on main (expanding) vs the frozen-box toy (Nov 17 2026, follow-up to #520).
Note: results/bound-stone-main-2026-11-17/README.md.  Nothing in engine/ or vendor/ is changed.

Same fixture, conversion, speed shim and readings as research/quiet_stone_main_v1.py (#520), with the stone's R and
x(0) as inputs. Added readings: the lapse at areal 1.5R and 3R (local proper time), the S-D energy inside the
COMOVING (coordinate) radius 3R, and S's energy half-radius, areal and coordinate (bound: the areal one does not
follow the Hubble flow).

Binding criterion (Derived; README): B = (45/64) x_w / (H^2 <r^2>) against the threshold 2.

Usage (repo root): python research/bound_stone_main_v1.py run TAG R XC DR RMAX T CFL OUTDIR     (resumable)
                   python research/bound_stone_main_v1.py toy R XC T OUTDIR
                   python research/bound_stone_main_v1.py analyse R OUTDIR
                   python research/bound_stone_main_v1.py criterion
"""
import math
import os
import sys

sys.path.insert(0, ".")
import numpy as np  # noqa: E402

import research.quiet_stone_main_v1 as q  # noqa: E402  (installs the bit-identical speed shim)
from engine.scalar_system import _d1  # noqa: E402

H_START = 0.00567          # production H/mu at the start (Read in #520)


def configure(R, xc):
    q.R0, q.XC = float(R), float(xc)
    q.SNAP_RMAX = 4.5 * q.R0
    q.SHELL = (1.5 * q.R0, 4.0 * q.R0)
    q.WINDOWS = (("b+6..b+200", (6, 200)), ("b+6..b+103", (6, 103)), ("b+103..b+200", (103, 200)),
                 ("pre-burst b-40..b-4", (-40, -4)))


_base_readings = q.readings


def readings(st):
    row = _base_readings(st)
    g, sc, grid = st.geometry, st.scalars, st.grid
    r = np.asarray(grid.centers, dtype=float); dr = float(grid.dr)
    RA = r * np.sqrt(g.b) / g.X
    grr_inv = g.X ** 2 / g.a
    Sp, Dp = _d1(grid, sc.S, 1), _d1(grid, sc.D, 1)
    Ds = q.sd_rest_value(sc.phi)
    e = q.LAM * (0.5 * (sc.PS ** 2 + sc.PD ** 2) + 0.5 * grr_inv * (Sp ** 2 + Dp ** 2)
                 + q.sd_potential(sc.S, sc.D, phi=sc.phi) - q.sd_potential(0.0 * sc.S, Ds, phi=sc.phi))
    dV = 4 * math.pi * r ** 2 * dr * np.sqrt(g.a) * g.b / g.X ** 3
    eS = q.LAM * (0.5 * sc.PS ** 2 + 0.5 * grr_inv * Sp ** 2 + 0.5 * sc.S ** 2) * dV
    inner = r < 0.75 * r[-1]
    c = np.cumsum(eS[inner]); k = int(np.searchsorted(c, 0.5 * c[-1]))
    a15, a3 = (float(np.interp(f * q.R0, RA, g.alpha)) for f in (1.5, 3.0))
    return row + [a15, a3, float(np.sum((e * dV)[r < 3 * q.R0])), float(RA[k]), float(r[k]), float(c[-1])]


q.readings = readings


def criterion():
    from research.stone_toy_v1 import quiet_profile
    r = np.linspace(0.005, 400.0, 80000); w = r ** 2
    out = []
    for R, xc in ((30, 0.3), (20, 0.3), (15, 0.3), (12, 0.3), (10, 0.3), (30, 1.2), (20, 1.0), (15, 0.6)):
        x = quiet_profile(r, R, xc, "sp")
        xw = float((x * x * w).sum() / (x * w).sum()); y = float((x * r * r * w).sum() / (x * w).sum())
        A = 45 / 64 * xw
        out.append(f"R {R:3d} x0 {xc:.1f}: x_w {xw:.3f} <r2> {y:6.0f} draw A {A:.4f} 2H^2<r2> {2 * H_START ** 2 * y:.4f} "
                   f"B = A/(H^2<r2>) {A / (H_START ** 2 * y):5.1f}  naive collapse sqrt(2<r2>/A) {math.sqrt(2 * y / A):.0f}")
    print("\n".join(out))
    return out


def toy(R, xc, T, outdir):
    R, xc = float(R), float(xc)
    from multiprocessing import Pool
    from research import stone_toy_v1 as st1
    base = dict(st1.CHAOS_BASE, R=float(R), xc=float(xc), T=float(T), g=st1.G_PROD, snap_r=4.5, snap_dt=q.SNAP_DT,
                probe=(1.5 * R,), rmax=240.0)
    cfgs = [dict(base, dr=dr, grav=False) for dr in (0.2, 0.1, 0.05)] + [dict(base, dr=0.1, grav=True)]
    with Pool(2) as p:
        res = p.map(st1.run, cfgs)
    for cfg, out in res:
        tag = f"toy_R{R:g}_x{xc:g}_dr{cfg['dr']}_{'grav' if cfg['grav'] else 'nograv'}"
        np.savez_compressed(os.path.join(outdir, tag + ".npz"), out=out, snap_t=cfg["_snap_t"],
                            snap_d=cfg["_snap_d"].astype(np.float32), rr=cfg["_rr"],
                            probe_t=np.arange(len(cfg["_probe_d"])) * cfg["_dt"],
                            probe_d=cfg["_probe_d"][:, 0].astype(np.float32), dt=cfg["_dt"])


def extra(outdir, tag):
    """Main-only readings: proper-time line, comoving-3R energy, S's half-radius."""
    z = np.load(os.path.join(outdir, tag + ".npz")); o = z["rows"]; t = o[:, 0]
    F15 = q.box(t, o[:, 5], q.W_FLUX); tb, _ = q.first_burst(t, F15)
    y = np.array([np.interp(1.5 * q.R0, z["snap_RA"][i], z["snap_d"][i]) for i in range(len(z["snap_t"]))])
    om = q.line(z["snap_t"], y, tb)
    m = (t >= tb - 16) & (t <= tb + 16)
    a15, a3 = o[m, 12].mean(), o[m, 13].mean()
    s = [f"## {tag} extra: line omega {om:.3f} per t; alpha(1.5R) {a15:.5f}, alpha(3R) {a3:.5f}, alpha(0) "
         f"{o[m, 1].mean():.5f} around t_b -> omega on proper time at 1.5R {om / a15:.3f}; alpha(1.5R) range over run "
         f"{o[:, 12].min():.5f}-{o[:, 12].max():.5f}"]
    pts = [tb + a for a in (-40, 0, 6, 128, 200, 250)] + [0.0]
    pts = sorted(p for p in pts if 0 <= p <= t[-1])
    s.append("   t | E(areal 3R)/E0 | E(comoving 3R)/E0 | S half-radius areal / coordinate | R_A/r(45-scale)")
    for p in pts:
        i = int(np.searchsorted(t, p)); i = min(i, len(t) - 1)
        s.append(f"   {t[i]:5.0f} | {o[i, 3] / o[0, 3]:.3f} | {o[i, 14] / o[0, 14]:.3f} | {o[i, 15]:.1f} / {o[i, 16]:.1f} | "
                 f"{o[i, 11]:.3f}")
    W = (tb + 6, min(tb + 200, t[-1]))
    s.append(f"   comoving-3R energy change over {W[0]:.0f}..{W[1]:.0f}: "
             f"{100 * (np.interp(W[1], t, o[:, 14]) / np.interp(W[0], t, o[:, 14]) - 1):+.1f}%")
    return s


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "run":
        tag, R, xc, dr, rmax, T, cfl, outdir = sys.argv[2:10]
        configure(R, xc)
        q.run(tag, float(dr), float(rmax), float(T), float(cfl), 0.0, outdir)
    elif mode == "toy":
        configure(sys.argv[2], sys.argv[3])
        toy(sys.argv[2], sys.argv[3], float(sys.argv[4]), sys.argv[5])
    elif mode == "analyse":
        configure(sys.argv[2], 0.0)
        outdir = sys.argv[3]
        q.analyse(outdir)
        lines = []
        for f in sorted(os.listdir(outdir)):
            if f.endswith(".npz") and not f.startswith("toy"):
                lines += extra(outdir, f[:-4])
        with open(os.path.join(outdir, "summary.txt"), "a") as fh:
            fh.write("\n".join(lines) + "\n")
        print("\n".join(lines))
    else:
        criterion()
