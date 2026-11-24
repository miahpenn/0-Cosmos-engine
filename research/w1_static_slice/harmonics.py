"""S harmonic split for the W1 corrected slice (W1 derivation AMENDMENT 6, section 3).

Builds S3-only and S1-only engine arrays from the stored harmonic-balance vectors (native.npz "<grid>_x")
through the producer's own mapped() from Aster's native.py, unchanged: the other S harmonic's entries in x are
zeroed. The per-node snapshot phase and the momentum recipe follow the slice's own construction
(research/w1_static_slice/build.py): sign +1 is the corrected slice, 0 the admitted slice, -1 the wrong-sign fault,
and no_v2 removes the v^2 P_st piece as the v^2-missing fault does.
"""
import math
import sys
from pathlib import Path

import numpy as np


def load(v10, native_py):
    v10 = Path(v10).resolve()
    if str(v10) not in sys.path:
        sys.path.insert(0, str(v10))
    from hbfull_u import HB
    from research import b7_w2_origin_check_v1 as oc
    import json
    settings = json.loads((v10 / "report.json").read_text())["settings"]
    u, lam, H0, w = [settings[k] for k in ("u", "lam", "H0", "omega_hat")]
    sq = math.sqrt(u)
    ns = {"np": np, "math": math, "oc": oc, "HB": HB, "u": u, "lam": lam, "H0": H0, "w": w, "sq": sq}
    src = (v10 / "step2_map.py").read_text()
    exec(src[src.index("def d1"):src.index("def member")], ns)
    src = Path(native_py).read_text()
    exec(src[src.index("def mapped"):src.index("hb=solver")], ns)
    ns["R"] = round(120 / (sq * .05)) * (sq * .05)
    return ns


class _Phased:
    def __init__(self, H, theta):
        self.H, self.theta = H, theta

    def __getattr__(self, k):
        return getattr(self.H, k)

    def snapshot(self, x, theta):
        return self.H.snapshot(x, theta + self.theta)


def split_S(ns, x, dr, r, sign=1, no_v2=False):
    """Return {S1, PS1, S3, PS3} engine arrays on cell centres r for one grid."""
    u, lam, H0, sq = ns["u"], ns["lam"], ns["H0"], ns["sq"]
    H = ns["solver"](dr * sq, ns["R"])
    assert np.max(np.abs(H.r[:len(r)] / sq - r)) < 2e-14
    th = sign * (ns["w"] * sq) * H0 * (H.r / sq) ** 2 / 2
    N, Nh = len(r), H.N
    out = {}
    for keep, drop in ((1, 3), (3, 1)):
        xs = np.array(x, dtype=np.float64, copy=True)
        for j, (f, n, p) in enumerate(H.U):
            if f == "s" and n == drop:
                xs[j * Nh:(j + 1) * Nh] = 0.0
        f = ns["mapped"](_Phased(H, th), xs, r)
        S, PS = f["S"], f["PS"]
        if no_v2:
            ps = H.snapshot(xs, th)[1]
            PS = PS - (H0 * r) ** 2 * ps[:N] * u / math.sqrt(lam)
        out["S%d" % keep], out["PS%d" % keep] = S, PS
    return out
