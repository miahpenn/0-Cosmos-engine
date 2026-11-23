"""Follow-up to the B7 W1 outer-cell trace (PR #570), research only, no engine change.

Gates fixed in Lumen2 Amendment 3 (1FD7pQV5-gKMMPn5AwQjjg_NGFDXQPI4eUBFw6BI5WGU) as amended by Amendment 4
(1gDXlmh3H9SwrisSV1VojFGT65k8p1dY6NPtX_9cyoLU); go from Jeremiah 2026-11-23 23:07 UTC. Run from a checkout of
aster/w1-profile (cfae43c) with this file and research/b7_w1_outer_cell_trace.py present:

    OPENBLAS_NUM_THREADS=1 python -m research.b7_w1_outer_cell_followup OUTDIR

G0 (first t=0.003) and T1 (every stage) first. F1: 64 steps per run (128 stage values). F2: 2400/1200 RMS ratio of the
arithmetic-only control at the outer cell in [1.3, 3.0], binding only when its relative 1-sigma (moving-block bootstrap)
is below 0.15. F4: control below 2F on both grids. F3 (report-only): control split by operation (a) stencil sum,
(b) 2(...)/r, (c) 1/b - 1/a incl. reciprocals, (e) other incl. final Lambda - Lambda_metric; (d) = departure difference
minus control.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

from research import b7_w1_outer_cell_trace as T

NSTEPS = 64
BLOCK = 8
NBOOT = 4000


def rms(x):
    x = np.asarray(x, float)
    return float(np.sqrt(np.mean(x * x)))


def block_resample(x, rng):
    x = np.asarray(x, float); n = len(x)
    nb = -(-n // BLOCK)
    starts = rng.integers(0, n - BLOCK + 1, size=nb)
    return np.concatenate([x[s:s + BLOCK] for s in starts])[:n]


def ratio_boot(fine, coarse, seed):
    rng = np.random.default_rng(seed)
    point = rms(fine) / rms(coarse)
    reps = np.array([rms(block_resample(fine, rng)) / rms(block_resample(coarse, rng)) for _ in range(NBOOT)])
    return point, float(np.std(reps) / point), [float(np.percentile(reps, 16)), float(np.percentile(reps, 84))]


def gates(runs):
    R = {(r["fixture"], r["n"], r["dt"]): r for r in runs}
    out = {}
    g0 = {f"{k[0]} {k[1]} {k[2]}": dict(measured=round(R[k]["G0_peak_e12"], 3), expected=v,
                                         exact=R[k]["G0_peak_e12"], match=round(R[k]["G0_peak_e12"], 3) == v)
          for k, v in T.G0_EXPECTED.items()}
    out["G0"] = dict(rows=g0, PASS=all(v["match"] for v in g0.values()))
    t1 = max(max(v for v in s["t1_worst_over_tol"].values()) for r in runs for s in r["steps"])
    wo = all(s["write_only_Aa_last"] for r in runs for s in r["steps"])
    out["T1"] = dict(worst_residual_over_tolerance=t1, write_touches_only_last_Aa=wo, PASS=bool(t1 <= 1.0 and wo))
    out["F1"] = dict(stage_values={f"{k[0]} {k[1]} {k[2]}": len(r["stages"]) for k, r in R.items()},
                     PASS=all(len(r["stages"]) >= 128 for r in runs))
    f2 = {}; seed = 20261123
    for fx in ("W1", "PC2"):
        for dt in (0.00075, 0.000375):
            c, f = R[(fx, 1200, dt)], R[(fx, 2400, dt)]
            xc = [s["diff_ctl"][2] for s in c["stages"]]; xf = [s["diff_ctl"][2] for s in f["stages"]]
            point, rel, band68 = ratio_boot(xf, xc, seed); seed += 1
            binds = rel < 0.15
            f2[f"{fx} {dt}"] = dict(ratio=point, rel_sigma=rel, boot_16_84=band68, binds=bool(binds),
                                    verdict=("PASS" if 1.3 <= point <= 3.0 else "FAIL") if binds else "REPORTED (scatter >= 0.15)",
                                    rms_1200=rms(xc), rms_2400=rms(xf))
    out["F2"] = dict(rows=f2, PASS=all(v["verdict"] == "PASS" for v in f2.values()),
                     any_unbound=any(not v["binds"] for v in f2.values()))
    f4 = {f"{r['fixture']} {r['n']} {r['dt']}": max(abs(s["diff_ctl"][2]) for s in r["stages"]) / r["F"] for r in runs}
    out["F4"] = dict(max_ctl_over_F=f4, PASS=all(v < 2.0 for v in f4.values()))
    f3 = {}
    for r in runs:
        st = r["stages"]
        parts = {p: rms([s["split"][p] for s in st]) for p in ("a", "b", "c", "e")}
        closure = max(abs(sum(s["split"][p] for p in "abce") - s["diff_ctl"][2]) for s in st)
        d = [s["diff_dep"][2] - s["diff_ctl"][2] for s in st]
        f3[f"{r['fixture']} {r['n']} {r['dt']}"] = dict(rms_parts=parts, rms_d=rms(d), rms_ctl=rms([s["diff_ctl"][2] for s in st]),
                                                       rms_dep_diff=rms([s["diff_dep"][2] for s in st]),
                                                       closure_max_abs=closure, F=r["F"])
    out["F3_report_only"] = f3
    out["drift_ulp_report_only"] = {f"{r['fixture']} {r['n']} {r['dt']}": r["drift_ulp"] for r in runs}
    return out


def main(outdir):
    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    runs = []
    for fx in ("W1", "PC2"):
        for n in (1200, 2400):
            for dt in (0.00075, 0.000375):
                rr = T.run(fx, n, dt, nsteps=NSTEPS, g0_t=0.003, split=True)
                runs.append(rr)
                print(fx, n, dt, "G0 peak %.4f e-12" % rr["G0_peak_e12"], flush=True)
    res = dict(gates=gates(runs), runs=runs)
    (outdir / "report.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
    print(json.dumps(res["gates"], indent=1, default=float))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "results/research/b7_w1_outer_cell_followup")
