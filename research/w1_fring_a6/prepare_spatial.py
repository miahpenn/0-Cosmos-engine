"""Corrected W1 static slice on extra grids for the spatial replays (rebuild; Aster's prepare_spatial is not in the repo).

Same path as Aster's native.py (new HB Newton solve per aligned grid, seeded by even cubic interpolation of the
archived A1 harmonics, same gamma seed, 1e-10 residual) followed by build.py's static phase
(theta = omega_u H0 r^2 / 2 per solver node). Positive control first: dr 0.025 must reproduce the stored HB vector
(native.npz dr0.025_x) and the committed corrected slice (native_static.npz dr0.025_*) BITWISE, else STOP.

Usage: python research/w1_fring_a6/prepare_spatial.py V10_DIR NATIVE_NPZ OUT.npz DR [DR ...]
"""
import hashlib
import json
import math
import os
import sys
from pathlib import Path

v10, native_npz, out = Path(sys.argv[1]).resolve(), Path(sys.argv[2]), Path(sys.argv[3])
drs = [float(x) for x in sys.argv[4:]]
sys.path[:0] = [os.getcwd(), str(v10)]
import numpy as np  # noqa: E402
from scipy.interpolate import CubicSpline  # noqa: E402

from hbfull_u import HB  # noqa: E402
from research import b7_w2_origin_check_v1 as oc  # noqa: E402

for line in (v10 / "MANIFEST.sha256").read_text().splitlines():
    h, n = line.split()
    assert hashlib.sha256((v10 / n).read_bytes()).hexdigest() == h, n
settings = json.loads((v10 / "report.json").read_text())["settings"]
u, lam, H0, w = [settings[k] for k in ("u", "lam", "H0", "omega_hat")]
sq = math.sqrt(u)
OMEGA_U = w * sq
ns = {"np": np, "math": math}
src = (v10 / "step2_map.py").read_text()
exec(src[src.index("def d1"):src.index("def member")], ns)
d1 = ns["d1"]


def mapped(H, x, r, theta):
    s, ps, d, pd = H.snapshot(x, theta)
    N = len(r)
    sp, dp = d1(s, H.dr)[:N], d1(d, H.dr)[:N]
    return dict(S=s[:N]*sq/math.sqrt(lam), D=d[:N]*sq/math.sqrt(lam),
                PS=(ps[:N]+H0*r*sp)*u/math.sqrt(lam), PD=(pd[:N]+H0*r*dp)*u/math.sqrt(lam))


def solver(dr, R):
    return HB(w, (1, 3), (0, 2, 4), dr=dr, R=R, mu=1/u)


hb = solver(.05, 120.)
with np.load(v10 / "hb_members.npz") as z:
    old, gamma = z["x_A1"], float(z["g_A1"])
seed = CubicSpline(np.r_[-hb.r[::-1], hb.r],
                   np.concatenate([old.reshape(-1, hb.N)[:, ::-1], old.reshape(-1, hb.N)], 1), axis=1)
go, vac, _ = oc.adapter.vendor_modules()
R = round(120/(sq*.05))*(sq*.05)
stored = np.load(native_npz)
static = np.load("results/research/w1_static_slice/native_static.npz")
save, meta = {}, {}
for dr in [.025] + [d for d in drs if d != .025]:
    tag = "dr%g" % dr
    grid = go.SphericalCellGrid(round(60/dr), 60.)
    r = np.asarray(grid.centers)
    H = solver(dr*sq, R)
    assert np.max(np.abs(H.r[:len(r)]/sq - r)) < 2e-14
    x, g, ok, it = H.solve(seed(H.r).ravel(), gamma)
    res = float(max(abs(H.residual(x, g))))
    assert ok and res < 1e-10
    f = mapped(H, x, r, OMEGA_U*H0*(H.r/sq)**2/2)
    meta[tag] = dict(gamma=g, iterations=it, resid=res)
    if dr == .025:
        ctrl = bool(np.array_equal(x, stored[tag+"_x"])) and all(
            f[k].tobytes() == np.asarray(static[tag+"_"+k]).tobytes() for k in ("S", "PS", "D", "PD"))
        meta[tag]["control_bitwise"] = ctrl
        print("control dr0.025 bitwise:", ctrl, flush=True)
        if not ctrl:
            raise SystemExit("STOP: rebuild does not reproduce the committed dr 0.025 slice")
    for k in f:
        save[tag+"_"+k] = f[k]
    save[tag+"_x"] = x
np.savez(out, **save)
out.with_suffix(".json").write_text(json.dumps(meta, indent=1))
print(meta, flush=True)
