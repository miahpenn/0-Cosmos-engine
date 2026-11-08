"""Read-only: split H_BSSN - H_metric_only into the Lambda part and a remainder.

H_BSSN is affine in Lambda with dH_BSSN = X^2 * D_r(dLambda) (the pinned tool's
discriminator checks this to ~5e-12). So, with C = Lambda - Lambda_metric from the
vendor's own constraint routine,

    H_BSSN - H_metric = X^2 * D_r(C)  +  remainder,

and the remainder vanishes in the continuum. Its behaviour under refinement tells
whether the non-Lambda difference is discretization error or a convention mismatch.
No production code is modified; the vendor and source pins are checked by the tool.

Usage: COSMOS_PINNED_SOURCE_ROOT=<checkout of 7f47669a with vendor submodule> \
       python tools/issue4_hamiltonian_slot_decomposition.py [output.json]
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import tools.audit_brown_vacuum_momentum_propagation as T  # noqa: E402

RESOLUTIONS = [40, 60, 80, 120, 160, 240, 320]
R_MAX = 40.0
EDGE = 4  # outer cells excluded from the interior region (stencil width)


def region(arr, sl):
    a = np.abs(arr[sl])
    return float(a.max()), float(np.sqrt(np.mean(arr[sl] ** 2)))


def run(n, zero_c):
    grid_ops, vacuum, _moving = T.adapter.vendor_modules()
    grid = grid_ops.SphericalCellGrid(n, R_MAX)
    state = T.manufactured_state(grid, vacuum)
    if zero_c:
        # Set Lambda to the vendor's own metric-derived value: C becomes ~0.
        c0 = np.asarray(vacuum.constraints(grid, state)["connection"])
        state.Lambda = state.Lambda - c0
    con = vacuum.constraints(grid, state)
    H_b = np.asarray(con["hamiltonian"], float)
    C = np.asarray(con["connection"], float)
    H_m, _ = T.metric_only_hamiltonian(grid, state)
    X = np.asarray(state.X, float)
    lam_part = X * X * T.d1(grid, C, -1)
    diff = H_b - H_m
    rest = diff - lam_part
    interior = slice(2, n - EDGE)
    h_scale = float(np.max(np.abs(H_b)))
    out = {"N": n, "dr": float(grid.centers[1] - grid.centers[0]),
           "max_abs_H_BSSN": h_scale,
           "max_abs_C": float(np.max(np.abs(C)))}
    for name, arr in (("diff", diff), ("lambda_part", lam_part), ("remainder", rest)):
        mx, rms = region(arr, interior)
        out[name] = {
            "interior_max_abs": mx, "interior_rms": rms,
            "interior_max_rel_to_max_H_BSSN": mx / h_scale if h_scale else None,
            "cell_0": float(arr[0]), "cell_1": float(arr[1]),
            "outer_edge_max_abs": float(np.max(np.abs(arr[n - EDGE:]))),
        }
    return out


def orders(rows, key):
    res = []
    for p, q in zip(rows, rows[1:]):
        e1, e2 = p[key]["interior_max_abs"], q[key]["interior_max_abs"]
        res.append(None if e1 <= 0 or e2 <= 0 else
                   float(np.log(e1 / e2) / np.log(q["N"] / p["N"])))
    return res


def main():
    prov = T.provenance()
    report = {"schema": "issue4_hamiltonian_slot_decomposition_v1", "read_only": True,
              "provenance": prov, "r_max": R_MAX, "edge_cells_excluded": EDGE}
    for label, zero_c in (("manufactured_C_nonzero", False), ("lambda_set_to_metric_C_zero", True)):
        rows = [run(n, zero_c) for n in RESOLUTIONS]
        report[label] = {
            "rows": rows,
            "observed_order_remainder_interior_max": orders(rows, "remainder"),
            "observed_order_diff_interior_max": orders(rows, "diff"),
        }
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "results" / "issue4-hamiltonian-slot-decomposition-2026-11-08" / "summary.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=1) + "\n")
    for label in ("manufactured_C_nonzero", "lambda_set_to_metric_C_zero"):
        print("==", label)
        print(f"{'N':>4} {'maxH':>10} {'diff':>10} {'lam_part':>10} {'remainder':>10} {'rem/maxH':>9} {'rem c0':>10} {'rem c1':>10}")
        for r in report[label]["rows"]:
            print(f"{r['N']:>4} {r['max_abs_H_BSSN']:10.3e} {r['diff']['interior_max_abs']:10.3e} "
                  f"{r['lambda_part']['interior_max_abs']:10.3e} {r['remainder']['interior_max_abs']:10.3e} "
                  f"{r['remainder']['interior_max_rel_to_max_H_BSSN']:9.2e} "
                  f"{r['remainder']['cell_0']:10.3e} {r['remainder']['cell_1']:10.3e}")
        print("order(remainder):", [None if o is None else round(o, 2) for o in report[label]["observed_order_remainder_interior_max"]])


if __name__ == "__main__":
    main()
