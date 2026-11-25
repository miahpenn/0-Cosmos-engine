"""Amendment 17 section 5: zero-step readback of the archived W1 states (no step taken, no engine import).

Reads the stored signed native-cell M and H profiles at t = 0, 0.00375, 0.0075, 0.01125, 0.015 (the five saves of the
dr 0.025 and dr 0.0125 W1 spatial replays, dt 3.75e-4 primary, 1.875e-4 supplementary) and reports, per A17 section 1:
  - native band cells r in [0.05, 0.25] (8 at .025, 16 at .0125) and the adjacent-pair sign changes (A12: >= 40%);
  - shares of the band sum of squares on A3's four coarse shells (r .075 .. .225; r^2-weighted shell means) and on
    the eight 0.025 shells, for reference (G5 also needs no shell above 40%);
  - where along the band the sign changes sit (pair index of each change), for A17's centre-spread signature;
  - band RMS of M and of M - M(t=0) against step count (A17: random walk ~ sqrt(steps), systematic ~ steps).
The signatures are A17's, declared before this read; nothing here is fitted.

Usage: python3 -m research.w1_fring_a6.a17_readback RUNS_DIR OUT.json   (RUNS_DIR/sp/<dr>-<dt>.npz)
"""
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

BAND = (0.05, 0.25)


def shells(v, r, width):
    edges = np.arange(BAND[0], BAND[1] + 1e-12, width)
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (r > lo) & (r < hi)
        out.append(float(np.sum(r[m] ** 2 * v[m]) / np.sum(r[m] ** 2)))
    return np.array(out)


def read(v, r):
    w = (r > BAND[0]) & (r < BAND[1])
    b = v[w]
    s = np.sign(b)
    ch = s[1:] != s[:-1]
    c4 = shells(v, r, 0.05)
    c8 = shells(v, r, 0.025)
    sh = lambda c: (c ** 2 / np.sum(c ** 2)).tolist()
    return {"n_cells": int(w.sum()), "pairs": int(w.sum()) - 1, "sign_changes": int(ch.sum()),
            "fraction": float(ch.mean()), "meets_40pct": bool(ch.mean() >= 0.4),
            "change_pair_index": np.nonzero(ch)[0].tolist(), "signs": "".join("+" if x > 0 else "-" for x in s),
            "coarse_shell_share": sh(c4), "coarse_max_share": float(max(sh(c4))),
            "shell8_share": sh(c8), "shell8_max_share": float(max(sh(c8))),
            "first_shell_share": float(sh(c4)[0]), "band_rms": float(np.sqrt(np.mean(b ** 2)))}


def main(runs, out):
    runs = Path(runs)
    res = {"source_files": {}, "grids": {}}
    for dr in ("0.025", "0.0125"):
        for dt in ("0.000375", "0.0001875"):
            f = runs / "sp" / f"{dr}-{dt}.npz"
            res["source_files"][f.name] = hashlib.sha256(f.read_bytes()).hexdigest()
            with np.load(f, allow_pickle=False) as z:
                meta = json.loads(str(z["metadata"]))
                t, r = z["t"], z["r"]
                rows = {}
                for q in ("M", "H"):
                    rows[q] = [dict(t=float(t[k]), step=int(round(t[k] / float(dt))), **read(z[q][k], r))
                               for k in range(len(t))]
                    w = (r > BAND[0]) & (r < BAND[1])
                    d0 = [float(np.sqrt(np.mean((z[q][k][w] - z[q][0][w]) ** 2))) for k in range(len(t))]
                    for k, x in enumerate(rows[q]):
                        x["band_rms_minus_t0"] = d0[k]
                    # log-slope of band RMS(M - M0) against steps between successive checkpoints
                    st = [x["step"] for x in rows[q]]
                    rows[q + "_growth_exponent_minus_t0"] = [
                        float(np.log(d0[k + 1] / d0[k]) / np.log(st[k + 1] / st[k])) for k in range(1, len(t) - 1)]
            res["grids"][f"{dr}/{dt}"] = {"metadata": meta, **rows}
    Path(out).write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
