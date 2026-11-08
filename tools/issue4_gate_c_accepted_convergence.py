"""Gate C, accepted class: locked full-Bianchi residual at N = 60, 120, 240.

Read-only driver around tools/issue4_full_bianchi_evaluator.py. For each resolution
it calls evaluate_accepted_fixture unchanged (the fixture loader checks every archive
SHA-256). The evaluator's accepted reference regression exists only for N = 60; for
N = 120 and 240 it is recorded as not executed instead of raising. Nothing else in
the evaluator is altered.

The three accepted states are separate t = 0.60 evolutions at each resolution
(dt halved with each doubling), so this is a refinement sequence of states, not a
pure spatial convergence test on one state. If the projected-Bianchi identity holds,
the tangent-minus-predictor residual should still fall under refinement at fixed
radius; a residual that levels off would indicate a missing or mismatched term.

Usage:
  python tools/issue4_gate_c_accepted_convergence.py \
      --archive 60=<zip> --archive 120=<zip> --archive 240=<zip> \
      --diagnostic-source <4aec88d0 checkout> --helper-baseline <7f47669a checkout> \
      [--output report.json]
"""
import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import tools.issue4_full_bianchi_evaluator as E  # noqa: E402

FIXED_RADII = (1.0, 2.0, 3.0, 5.0, 8.0)
TERMS = ("normal_alpha_K_M", "normal_lapse_gradient_H", "normal_H_gradient", "connection_divergence")


def regions(values, r, signal=None):
    x = np.asarray(values, float)
    inner = r <= 10.0
    annulus = (r >= 5.0) & (r <= 10.0)
    out = {
        "max_abs_r_le_10": float(np.max(np.abs(x[inner]))),
        "radius_of_max_r_le_10": float(r[inner][int(np.argmax(np.abs(x[inner])))]),
        "rms_annulus_5_10": float(np.sqrt(np.mean(x[annulus] ** 2))),
        "cells_0_1_2": [float(x[i]) for i in range(3)],
        "fixed_radius": {f"{t:g}": float(np.interp(t, r, x)) for t in FIXED_RADII},
    }
    if signal is not None:
        s = np.asarray(signal, float)
        out["relative_rms_annulus_5_10"] = out["rms_annulus_5_10"] / max(
            float(np.sqrt(np.mean(s[annulus] ** 2))), np.finfo(float).tiny)
    return out


def run(n, archive, diag, helper):
    original = E.assert_reference_regression

    def regression(result):
        if int(result["resolution"]) == 60:
            return original(result)
        return {"executed": False, "reason": f"no archived reference regression for N={result['resolution']}"}

    E.assert_reference_regression = regression
    try:
        res = E.evaluate_accepted_fixture(archive, n, diag, helper)
    finally:
        E.assert_reference_regression = original
    r = np.asarray(res["grid_centers"], float)
    out = {"N": n, "status": res["status"], "production_steps_run": res["production_steps_run"],
           "zero_step_tripwire": res["zero_step_tripwire"], "input_digests": res["input_digests"],
           "accepted_reference_regression": res["accepted_reference_regression"],
           "ricci_offset_crosscheck": res["ricci_offset_crosscheck"], "channels": {}, "terms": {}}
    st = res["signed_terms"]
    for k in TERMS + ("H_total", "M_total", "C_r"):
        out["terms"][k] = regions(st[k], r)
    for eps, chans in res["channels"].items():
        out["channels"][eps] = {}
        for name in ("pointwise", "finite_volume"):
            ch = chans[name]
            meas = np.asarray(ch["measured_tangent"], float)
            out["channels"][eps][name] = {
                "measured_tangent": regions(meas, r),
                "predictor": regions(ch["predictor"], r),
                "residual": regions(ch["residual_tangent_minus_predictor"], r, signal=meas),
            }
    eps_a, eps_b = list(res["channels"])[:2]
    out["epsilon_sensitivity_max_abs"] = {
        name: float(np.max(np.abs(np.asarray(res["channels"][eps_a][name]["residual_tangent_minus_predictor"])
                                  - np.asarray(res["channels"][eps_b][name]["residual_tangent_minus_predictor"]))))
        for name in ("pointwise", "finite_volume")}
    return out


def order(values, ns):
    return [None if a == 0 or b == 0 else round(math.log(abs(a) / abs(b)) / math.log(m / k), 2)
            for a, b, k, m in zip(values, values[1:], ns, ns[1:])]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--archive", action="append", required=True, help="N=path")
    p.add_argument("--diagnostic-source", type=Path, required=True)
    p.add_argument("--helper-baseline", type=Path, required=True)
    p.add_argument("--output", type=Path)
    a = p.parse_args()
    archives = dict(item.split("=", 1) for item in a.archive)
    ns = sorted(int(k) for k in archives)
    rows = [run(n, archives[str(n)], a.diagnostic_source, a.helper_baseline) for n in ns]
    eps = "eps1e-05"
    orders = {}
    for name in ("pointwise", "finite_volume"):
        get = lambda f: [f(row["channels"][eps][name]["residual"]) for row in rows]  # noqa: E731
        orders[name] = {
            "max_abs_r_le_10": order(get(lambda s: s["max_abs_r_le_10"]), ns),
            "rms_annulus_5_10": order(get(lambda s: s["rms_annulus_5_10"]), ns),
            "cell_0": order(get(lambda s: s["cells_0_1_2"][0]), ns),
            **{f"r_{t:g}": order(get(lambda s, t=t: s["fixed_radius"][f"{t:g}"]), ns) for t in FIXED_RADII},
        }
    rep = {"schema": "issue4_gate_c_accepted_convergence_v1", "read_only": True,
           "scientific_admission": False, "acceptance_oracle_satisfied": False,
           "resolutions": ns, "rows": rows, "residual_orders_eps1e-05": orders}
    text = json.dumps(rep, indent=1)
    if a.output:
        a.output.parent.mkdir(parents=True, exist_ok=True)
        a.output.write_text(text + "\n")
    print("GATE_C_REPORT_BEGIN")
    print(json.dumps(rep, separators=(",", ":")))
    print("GATE_C_REPORT_END")
    for row in rows:
        for name in ("pointwise", "finite_volume"):
            s = row["channels"][eps][name]["residual"]
            m = row["channels"][eps][name]["measured_tangent"]
            print(f"N={row['N']:>3} {name:13s} residual max(r<=10)={s['max_abs_r_le_10']:.3e} @r={s['radius_of_max_r_le_10']:.2f} "
                  f"annulusRMS={s['rms_annulus_5_10']:.3e} (rel {s['relative_rms_annulus_5_10']:.2e}) "
                  f"c0..2={['%+.2e' % v for v in s['cells_0_1_2']]} "
                  f"fixed={ {k: '%+.2e' % v for k, v in s['fixed_radius'].items()} } "
                  f"| measured max={m['max_abs_r_le_10']:.2e}")
    print("orders:", json.dumps(orders))


if __name__ == "__main__":
    main()
