"""Corrected W1 slice in static time (W1 derivation AMENDMENT 3, Drive 1JU2EhMUs7dSio23VVu0WIPMt33iAY6u9gbD9EKE64aY,
build path per AMENDMENT 4, Drive 1J_dGcUVzAy9jceOwxByyvuQ1pyix9hHiC9owyoRNItM, section 1). No simulation.

Each grid's stored harmonic-balance solution x (native.npz "<grid>_x") goes through the producer's own mapped()
from Aster's native.py, unchanged, with the solver's snapshot phase set per node to
theta(r) = sign * omega_u * H0 * r**2 / 2 (theta0 = 0) on every solver node (engine r = r_hat / sqrt(u)).
mapped() then gives P = P_st + H0 r d1(f), so P - v f_R = P_st up to stencil error. No Newton re-solve.

Positive control (A4 1): sign 0 is the producer's own path and must equal the admitted native.npz BITWISE on both
grids; anything else is a STOP and no other variant is written.

Usage: python3 research/w1_static_slice/build.py V10_DIR NATIVE_PY NATIVE_NPZ OUTDIR
  V10_DIR    extracted ringkit_v10_chk-w1-slice.tgz (Drive 16LGCwbEof7j0E84GSr21J9HsKDMpzDKn), v10/ folder
  NATIVE_PY  native.py from Aster's W1 bundle (Drive 1OUObtOwD7B51pDMiTacKJEjbavVBqv6D)
Writes OUTDIR/static/native.npz (the corrected slice), OUTDIR/build.json, and fault variants
OUTDIR/fault_no_v2/native.npz (P = P_st + v GS_st, the v^2 P_st piece removed) and
OUTDIR/fault_wrong_sign/native.npz (theta = -omega_u H0 r^2 / 2).
"""
import hashlib
import json
import math
import os
import sys
from pathlib import Path

v10, native_py, native_npz, out = map(Path, sys.argv[1:5])
v10, out = v10.resolve(), out.resolve()
sys.path[:0] = [os.getcwd(), str(v10)]
import numpy as np  # noqa: E402

from hbfull_u import HB  # noqa: E402
from research import b7_w2_origin_check_v1 as oc  # noqa: E402


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


for line in (v10 / "MANIFEST.sha256").read_text().splitlines():
    h, n = line.split()
    assert sha(v10 / n) == h, n
NATIVE_PY_SHA = sha(native_py)
NATIVE_NPZ_SHA = sha(native_npz)
assert NATIVE_NPZ_SHA == "16151c445c9dcfcc34cf61e28e268e991ab2bfea629beaa3aa52fadc4c5a98ea"

settings = json.loads((v10 / "report.json").read_text())["settings"]
u, lam, H0, w = [settings[k] for k in ("u", "lam", "H0", "omega_hat")]
sq = math.sqrt(u)
OMEGA_U = w * sq
# The producer's own code, unchanged: d1 from step2_map.py (as native.py takes it), mapped() and solver() from native.py.
ns = {"np": np, "math": math, "oc": oc, "HB": HB, "u": u, "lam": lam, "H0": H0, "w": w, "sq": sq}
src = (v10 / "step2_map.py").read_text()
exec(src[src.index("def d1"):src.index("def member")], ns)
src = native_py.read_text()
exec(src[src.index("def mapped"):src.index("hb=solver")], ns)
mapped, solver = ns["mapped"], ns["solver"]


class Phased:
    """The solver with its snapshot phase shifted per node; everything else passes through."""

    def __init__(self, H, theta):
        self.H, self.theta = H, theta

    def __getattr__(self, k):
        return getattr(self.H, k)

    def snapshot(self, x, theta):
        return self.H.snapshot(x, theta + self.theta)


go, vac, _ = oc.adapter.vendor_modules()
R = round(120 / (sq * .05)) * (sq * .05)
z = np.load(native_npz)
FIELDS = ("S", "PS", "D", "PD")
arrays = {"static": {}, "fault_no_v2": {}, "fault_wrong_sign": {}}
build = {"v10_manifest": "verified", "native_py_sha256": NATIVE_PY_SHA, "native_npz_sha256": NATIVE_NPZ_SHA,
         "settings": {"u": u, "lam": lam, "H0": H0, "omega_hat": w, "omega_u": OMEGA_U},
         "theta": "sign * omega_u * H0 * r**2 / 2 on every solver node, r = r_hat / sqrt(u), theta0 = 0",
         "harmonics": {"S": [1, 3], "D": [0, 2, 4]}, "grids": {}}
for dr in (.05, .025):
    tag = "dr%g" % dr
    grid = go.SphericalCellGrid(round(60 / dr), 60.)
    r = np.asarray(grid.centers)
    H = solver(dr * sq, R)
    assert np.max(np.abs(H.r[:len(r)] / sq - r)) < 2e-14
    x = np.asarray(z[tag + "_x"])
    rn = H.r / sq
    th = OMEGA_U * H0 * rn**2 / 2
    f0 = mapped(Phased(H, np.zeros(H.N)), x, r)
    ctrl = {k: {"array_equal": bool(np.array_equal(f0[k], z[tag + "_" + k])),
                "bytes_equal": f0[k].tobytes() == np.asarray(z[tag + "_" + k]).tobytes()} for k in FIELDS}
    ok = all(c["bytes_equal"] for c in ctrl.values())
    build["grids"][tag] = {"control_theta0": ctrl, "control_PASS": ok, "theta_last_node": float(th[-1]),
                           "solver_nodes": int(H.N), "engine_cells": len(r)}
    if not ok:
        continue
    fs = mapped(Phased(H, th), x, r)
    fw = mapped(Phased(H, -th), x, r)
    # v^2 piece removed: P_st from the phased snapshot; P = P_st + v (d1(f) - v P_st)
    s_, ps_, d_, pd_ = H.snapshot(x, th)
    N = len(r)
    v = H0 * r
    pst = {"PS": ps_[:N] * u / math.sqrt(lam), "PD": pd_[:N] * u / math.sqrt(lam)}
    for k in FIELDS:
        arrays["static"][tag + "_" + k] = fs[k]
        arrays["fault_wrong_sign"][tag + "_" + k] = fw[k]
        arrays["fault_no_v2"][tag + "_" + k] = fs[k] - v**2 * pst[k] if k in pst else fs[k]
    for a in arrays.values():
        a[tag + "_x"] = x
build["control_PASS"] = all(g["control_PASS"] for g in build["grids"].values())
out.mkdir(parents=True, exist_ok=True)
if build["control_PASS"]:
    for name, a in arrays.items():
        (out / name).mkdir(exist_ok=True)
        np.savez(out / name / "native.npz", **a)
        build["sha256_" + name] = sha(out / name / "native.npz")
        build["sha256_" + name + "_per_grid"] = {
            g: hashlib.sha256(b"".join(np.ascontiguousarray(a[g + "_" + k]).tobytes() for k in FIELDS)).hexdigest()
            for g in ("dr0.025", "dr0.05")}
(out / "build.json").write_text(json.dumps(build, indent=1))
print("CONTROL", "PASS (bitwise, both grids)" if build["control_PASS"] else "FAIL - STOP", flush=True)
