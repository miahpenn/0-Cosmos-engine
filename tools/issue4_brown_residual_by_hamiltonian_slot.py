"""Read-only: Brown Eq. (11b) momentum residual under two Hamiltonian slots.

Predictors:
  metric_only : Brown with the pinned tool's warped-product Ricci Hamiltonian (current default)
  bssn_minus_C: Brown with H_BSSN - X^2 D_r(C), frozen in Gate A section 11
  bssn_minus_C_plus_vendor_C: bssn_minus_C plus the candidate vendor-specific term
      alpha X^2 [D_r(u) + 3 u / r],  u = (D_r a / 2a) C,
      found by fitting the converged residual of bssn_minus_C (diagnostic candidate only)

For each resolution, the measured two-sided tangent of M_r along the pinned vacuum RHS
is compared with Brown's predictor built with each slot. Cells 0, 1, the interior
(cells 2..n-5) and the outer 4 cells are reported separately, with fixed-radius samples.
No production code is modified; the tool's source and vendor pins are enforced.

Usage: COSMOS_PINNED_SOURCE_ROOT=<checkout of 7f47669a with vendor submodule> \
       python tools/issue4_brown_residual_by_hamiltonian_slot.py [output.json]
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import tools.audit_brown_vacuum_momentum_propagation as T  # noqa: E402

RES = [40, 60, 80, 120, 160, 240, 320, 480, 640, 1280]
KEYS = ("metric_only", "bssn_minus_C", "bssn_minus_C_plus_vendor_C")
R_MAX, EDGE = 40.0, 4
EPSILONS = (2.0e-5, 1.0e-5)
FIXED_RADII = (1.0, 2.0, 5.0, 10.0)


def slot_bssn_minus_c(grid, vacuum, state):
    c = vacuum.constraints(grid, state)
    X = np.asarray(state.X, float)
    return np.asarray(c["hamiltonian"], float) - X * X * T.d1(grid, np.asarray(c["connection"], float), -1)


def prediction(grid, state, constraints, H_slot):
    """Brown predictor with the H slot replaced; the tool's term formulas are reused."""
    original = T.metric_only_hamiltonian
    try:
        T.metric_only_hamiltonian = lambda _g, _s: (H_slot, {})
        terms, pred = T.brown_lower_momentum_rhs(grid, state, constraints)
    finally:
        T.metric_only_hamiltonian = original
    return terms, pred


def regions(res, meas, x, n):
    def rel(sl):
        e, m = res[sl], meas[sl]
        return float(np.sqrt(np.mean(e * e)) / max(np.sqrt(np.mean(m * m)), np.finfo(float).tiny))
    inner = slice(2, n - EDGE)
    return {
        "cell_0": float(res[0]), "cell_1": float(res[1]),
        "interior_max_abs": float(np.max(np.abs(res[inner]))),
        "interior_relative_rms": rel(inner),
        "whole_grid_relative_rms": rel(slice(0, n)),
        "outer_edge_max_abs": float(np.max(np.abs(res[n - EDGE:]))),
        "fixed_radius": {f"{t:g}": float(np.interp(t, x, res)) for t in FIXED_RADII},
    }


def run(n):
    grid_ops, vacuum, moving = T.adapter.vendor_modules()
    grid = grid_ops.SphericalCellGrid(n, R_MAX)
    state = T.manufactured_state(grid, vacuum)
    x = np.asarray(grid.centers, float)
    con = vacuum.constraints(grid, state)
    rates = T.vacuum_rhs_rates(grid, state, vacuum, moving)
    H_metric, _ = T.metric_only_hamiltonian(grid, state)
    slots = {"metric_only": H_metric, "bssn_minus_C": slot_bssn_minus_c(grid, vacuum, state)}
    preds = {k: prediction(grid, state, con, h)[1] for k, h in slots.items()}
    a = np.asarray(state.a, float)
    X = np.asarray(state.X, float)
    u = T.d1(grid, a, +1) / (2.0 * a) * np.asarray(con["connection"], float)
    vendor_c = np.asarray(state.alpha, float) * X * X * (T.d1(grid, u, -1) + 3.0 * u / x)
    preds["bssn_minus_C_plus_vendor_C"] = preds["bssn_minus_C"] + vendor_c
    row = {"N": n, "slot_difference_max_abs": float(np.max(np.abs(slots["bssn_minus_C"] - H_metric))),
           "predictor_difference_max_abs": float(np.max(np.abs(preds["bssn_minus_C"] - preds["metric_only"]))),
           "vendor_C_term_max_abs": float(np.max(np.abs(vendor_c))),
           "epsilons": {}}
    for eps in EPSILONS:
        meas = T.tangent_momentum(grid, state, rates, vacuum, eps)
        row["epsilons"][f"{eps:g}"] = {
            "measured_max_abs": float(np.max(np.abs(meas))),
            "measured_fixed_radius": {f"{t:g}": float(np.interp(t, x, meas)) for t in FIXED_RADII},
            **{k: regions(meas - p, meas, x, n) for k, p in preds.items()},
        }
    return row


def order(vals, ns):
    return [None if a == 0 or b == 0 else round(math.log(abs(a) / abs(b)) / math.log(m / k), 2)
            for a, b, k, m in zip(vals, vals[1:], ns, ns[1:])]


def main():
    rows = [run(n) for n in RES]
    eps = f"{EPSILONS[1]:g}"
    summary = {}
    for k in KEYS:
        summary[k] = {
            "order_cell_0": order([r["epsilons"][eps][k]["cell_0"] for r in rows], RES),
            "order_interior_max": order([r["epsilons"][eps][k]["interior_max_abs"] for r in rows], RES),
            **{f"order_fixed_r_{t:g}": order([r["epsilons"][eps][k]["fixed_radius"][f"{t:g}"] for r in rows], RES)
               for t in FIXED_RADII},
        }
    rep = {"schema": "issue4_brown_residual_by_hamiltonian_slot_v1", "read_only": True,
           "scientific_admission": False, "provenance": T.provenance(), "r_max": R_MAX,
           "edge_cells_excluded": EDGE, "rows": rows, "orders_at_epsilon": {eps: summary}}
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "results" / "issue4-hamiltonian-slot-decomposition-2026-11-08" / "brown_residual_by_slot.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rep, indent=1) + "\n")
    print(f"eps={eps}")
    print(f"{'N':>4} {'|meas|max':>10} | {'metric c0':>10} {'c1':>10} {'int max':>9} {'int relrms':>10} | "
          f"{'slot c0':>10} {'c1':>10} {'int max':>9} {'int relrms':>10}")
    for r in rows:
        e = r["epsilons"][eps]
        m, s = e["metric_only"], e["bssn_minus_C"]
        print(f"{r['N']:>4} {e['measured_max_abs']:10.3e} | {m['cell_0']:+10.3e} {m['cell_1']:+10.3e} "
              f"{m['interior_max_abs']:9.2e} {m['interior_relative_rms']:10.2e} | {s['cell_0']:+10.3e} "
              f"{s['cell_1']:+10.3e} {s['interior_max_abs']:9.2e} {s['interior_relative_rms']:10.2e}")
    for k, v in summary.items():
        print(k, v)
    eps_a, eps_b = (f"{e:g}" for e in EPSILONS)
    dev = max(abs(r["epsilons"][eps_a][k]["interior_max_abs"] - r["epsilons"][eps_b][k]["interior_max_abs"])
              for r in rows for k in KEYS)
    print("max change in interior residual between epsilons:", f"{dev:.2e}")
    print("fixed-radius residual / local measured signal, eps", eps)
    for r in rows:
        e = r["epsilons"][eps]
        print(f"  N={r['N']:>3} " + "  ".join(
            f"r={t:g}: slot {e['bssn_minus_C']['fixed_radius'][f'{t:g}']:+.3e} metric {e['metric_only']['fixed_radius'][f'{t:g}']:+.3e}"
            f" +vendorC {e['bssn_minus_C_plus_vendor_C']['fixed_radius'][f'{t:g}']:+.3e}"
            f" signal {e['measured_fixed_radius'][f'{t:g}']:+.2e}" for t in FIXED_RADII))


if __name__ == "__main__":
    main()
