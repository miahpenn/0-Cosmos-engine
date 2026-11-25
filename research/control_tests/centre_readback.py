"""Readback of #606's flat-vacuum (1a) centre growth from stored saves only; no engine stepping.

Question (Jeremiah's "Do the readback", Nov 25 2026 21:47 UTC, relayed): is the 1a centre growth at dr .025 rounding seen
through the derivative, or a centre instability amplifying rounding? Verdict rule fixed before any number was read:
  ROUNDING-LIKE: over half of ln(final/initial) by step 10 AND rho_pow <= rho_exp with p + 2SE <= 2 (window fit).
  INSTABILITY:   rho_exp < rho_pow/2, gamma - 2SE > 0 and N_e = 1.8 gamma >= 3.
  Otherwise OPEN.
#606 kept full H and M profiles only at its 8 saves (t 0, .3, ..., 2.1); steps 1, 10 and 100 and t .5, 1 were not stored.
  python -m research.control_tests.centre_readback OUT
"""
import json
import sys
from pathlib import Path

import numpy as np

ARMS = Path("results/research/control_tests_stage1/arms")
ARM = "1a_dr025_dt0_r8"
CONTEXT = ("1a_dr05_dt0_r8",)   # reading only: the grid where #606 saw H(0) level off
EPS = float(np.finfo(float).eps)
WINDOW = (0.3, 2.1)


def fit(x, y):
    A = np.vstack([np.ones_like(x), x]).T
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    r = y - A @ c
    s2 = float(r @ r) / (len(y) - 2)
    se = float(np.sqrt(s2 * np.linalg.inv(A.T @ A)[1, 1]))
    return float(c[1]), se, float(np.sqrt(np.mean(r**2)))


def arm(name):
    rec = json.loads((ARMS / f"{name}.json").read_text())
    with np.load(ARMS / f"{name}.npz") as z:
        H = [np.array(z[f"t{k}_H"]) for k in range(len(rec["saves"]))]
        M = [np.array(z[f"t{k}_M"]) for k in range(len(rec["saves"]))]
    t = np.array([s["t"] for s in rec["saves"]])
    dr = rec["dr"]
    out = {"arm": name, "dr": dr, "dt": rec["dt"], "n": rec["n"], "save_t": t.tolist(),
           "cells_0_2": {q: {f"{tt:g}": [float(abs(v)) for v in arr[:3]] for tt, arr in zip(t, src)}
                         for q, src in (("H", H), ("M", M))},
           "eps_over_dr2": EPS / dr**2, "not_stored": ["steps 1, 10, 100", "t 0.5", "t 1.0"]}
    fits = {}
    for q, src in (("H0", H), ("M0", M)):
        e = np.maximum.accumulate(np.array([abs(a[0]) for a in src]))
        w = (t >= WINDOW[0] - 1e-12) & (t <= WINDOW[1] + 1e-12) & (e > 0) & np.isfinite(e)
        y = np.log(e[w])
        p, sp, rp = fit(np.log(t[w]), y)
        g, sg, rx = fit(t[w], y)
        inst = rx < rp / 2 and g - 2 * sg > 0 and 1.8 * g >= 3
        fits[q] = {"envelope": e.tolist(), "rows": int(w.sum()), "p": p, "SE_p": sp, "rho_pow": rp,
                   "gamma": g, "SE_gamma": sg, "rho_exp": rx, "N_e": 1.8 * g,
                   "window_fit_rounding_part": bool(rp <= rx and p + 2 * sp <= 2),
                   "verdict": "INSTABILITY" if inst else "OPEN",
                   "note": "ROUNDING-LIKE needs the step-10 share, which #606 did not store; it cannot be reached"}
    out["fits"] = fits
    return out


def main(out_dir):
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    report = {"primary": arm(ARM), "context_reading_only": {n: arm(n) for n in CONTEXT}}
    (out_dir / "readback.json").write_text(json.dumps(report, indent=1))
    p = report["primary"]
    print(json.dumps({"cells_0_2": p["cells_0_2"], "eps_over_dr2": p["eps_over_dr2"],
                      "fits": {k: {kk: v for kk, v in f.items() if kk != "envelope"} for k, f in p["fits"].items()},
                      "envelopes": {k: f["envelope"] for k, f in p["fits"].items()},
                      "context": {n: {k: {kk: f[kk] for kk in ("p", "SE_p", "rho_pow", "gamma", "SE_gamma", "rho_exp", "N_e", "verdict")}
                                      for k, f in a["fits"].items()} for n, a in report["context_reading_only"].items()}},
                     indent=1))


if __name__ == "__main__":
    main(sys.argv[1])
