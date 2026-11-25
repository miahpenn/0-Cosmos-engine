"""Solve the corrected W1 member on the dr=0.0125 engine grid."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.interpolate import CubicSpline

from research import b7_w2_origin_check_v1 as oc
from research.w1_fring_smoke.origin_trace import DIGEST, INPUT, ad

DR = 0.0125
FIELDS = ("S", "PS", "D", "PD")
EDGE_TOL = 1.943e-7  # Same entry threshold as the archived dr=0.025 diagnostic.


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ringkit", type=Path)
    parser.add_argument("producer", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    native = args.producer / "native.py"
    solution = args.producer / "native.npz"
    pins = {
        native: "81f8fa3a7675bcf4315a3cabed7cb3b7c56f71ce7801c95089c11b9837f46028",
        solution: "16151c445c9dcfcc34cf61e28e268e991ab2bfea629beaa3aa52fadc4c5a98ea",
        INPUT: DIGEST,
    }
    for line in (args.ringkit / "MANIFEST.sha256").read_text().splitlines():
        digest, path = line.split()
        pins[args.ringkit / path] = digest
    for path, digest in pins.items():
        if sha(path) != digest:
            raise ValueError(f"archive pin mismatch for {path}")
    sys.path.insert(0, str(args.ringkit.resolve()))
    from hbfull_u import HB

    settings = json.loads((args.ringkit / "report.json").read_text())["settings"]
    u, lam, H0, w = [settings[key] for key in ("u", "lam", "H0", "omega_hat")]
    sq = math.sqrt(u)
    ns = dict(np=np, math=math, oc=oc, HB=HB, u=u, lam=lam, H0=H0, w=w, sq=sq)
    src = (args.ringkit / "step2_map.py").read_text()
    exec(src[src.index("def d1"):src.index("def member")], ns)
    src = native.read_text()
    exec(src[src.index("def mapped"):src.index("hb=solver")], ns)
    mapped, solver = ns["mapped"], ns["solver"]

    class StaticPhase:
        def __init__(self, member):
            self.member = member

        def __getattr__(self, name):
            return getattr(self.member, name)

        def snapshot(self, x, theta):
            r = self.member.r/sq
            return self.member.snapshot(x, theta+w*sq*H0*r**2/2)

    go, _, _ = ad.vendor_modules()
    radius = round(120/(sq*.05))*sq*.05
    coarse = solver(.025*sq, radius)
    with np.load(solution, allow_pickle=False) as data:
        old = data["dr0.025_x"]
    coarse_r = go.SphericalCellGrid(2400, 60.).centers
    control = mapped(StaticPhase(coarse), old, coarse_r)
    with np.load(INPUT, allow_pickle=False) as data:
        for name in FIELDS:
            np.testing.assert_array_equal(control[name], data[f"dr0.025_{name}"])

    fine = solver(DR*sq, radius)
    profiles = old.reshape(-1, coarse.N)
    seed = CubicSpline(np.r_[-coarse.r[::-1], coarse.r],
                       np.concatenate([profiles[:, ::-1], profiles], axis=1), axis=1)
    member = json.loads((args.producer / "members.json").read_text())["dr0.025"]
    x, gamma, success, iterations = fine.solve(seed(fine.r).ravel(), member["gamma"])
    residual = float(np.max(np.abs(fine.residual(x, gamma))))
    if not success or residual >= 1e-10:
        raise ValueError(f"fine harmonic-balance solve failed: {residual}")
    r = go.SphericalCellGrid(round(60/DR), 60.).centers
    offset = float(np.max(np.abs(fine.r[:len(r)]/sq-r)))
    if offset >= 2e-14:
        raise ValueError("solver nodes are off-grid")
    fields = mapped(StaticPhase(fine), x, r)
    args.output.mkdir(parents=True, exist_ok=True)
    npz = args.output / "static_0p0125.npz"
    np.savez_compressed(npz, r=r, x=x, **fields)
    report = {"input_sha256": sha(npz), "source_sha256": {path.name: digest for path,digest in pins.items()},
              "coarse_control_bitwise": True, "n": len(r), "hb_nodes": fine.N, "dr": DR,
              "hb_iterations": iterations, "hb_residual": residual, "gamma": gamma,
              "node_offset": offset, "u": u, "H0": H0, "omega_hat": w, "edge_tol": EDGE_TOL}
    (args.output / "input.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
    print(json.dumps(report, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
