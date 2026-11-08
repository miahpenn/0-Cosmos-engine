"""Read-only: derived vendor-specific connection term in the vacuum M_r tangent.

Derivation (pinned vendor vacuum_bssn.py at d6052d60):
  geometry_terms: Lambda enters R^r_r as -(X^2/a)[-a D1(Lambda) - (1/2) Lambda D1(a)] and R as
  -(X^2/a)[-a D1(Lambda)].  With Lambda = Lambda_metric + C,
      alpha (R^r_r - R/3) = [C = 0 value] + alpha X^2 [ (2/3) D1(C) + u ],   u = (D1 a / 2a) C.
  primary_l2_rhs: L2^{A_a} contains alpha (R^r_r - R/3); L2^K = -lap(alpha) has no Lambda.
  momentum_constraint: M_r = D1(A_a) - (2/3) D1(K) + 6 A_a phi' + (3/2) A_a q,
      phi' = -D1(X)/(2X),  q = 2/r + D1(b)/b.
  So the A_a-rate piece delta = alpha X^2 u adds to the M_r tangent exactly
      dM = D1(delta) + (6 phi' + (3/2) q) delta,
  whose leading part is alpha X^2 [D1(u) + 3u/r]. Brown Eq. (11b) carries only the (2/3) D1(C) part.

Checks per resolution:
  - residual_after: measured tangent - Brown predictor (Gate A section 11 slot) - dM
  - rate_removal: tangent with delta removed from the A_a rate, compared with measured - dM
  - fit_form_difference: dM against the fitted form alpha X^2 [D1(u) + 3u/r]
  - locked_form: the geometric part of the locked full-Bianchi evaluator
    (tools/issue4_full_bianchi_evaluator.py: H and M at Lambda = Lambda_metric, and
    D1(alpha W) + (6 phi' + 3/r + 3 D1(b)/(2b)) alpha W with
    W = X^2 [(2/3) D1(C) + (D1 a / 2a) C]), transcribed here with rho = j = 0,
    compared with the measured vacuum tangent

Usage: COSMOS_PINNED_SOURCE_ROOT=<checkout of 7f47669a with vendor submodule> \
       python tools/issue4_vendor_connection_term_derivation.py [output.json]
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import tools.issue4_brown_residual_by_hamiltonian_slot as B  # noqa: E402

T = B.T
RES = [40, 80, 160, 320, 640, 1280]
EPS = 1.0e-5
FIXED_RADII = (1.0, 2.0, 5.0, 10.0)
EDGE = 4


def derived_term(grid, state, connection):
    x = np.asarray(grid.centers, float)
    a, b, X, alpha = (np.asarray(getattr(state, k), float) for k in ("a", "b", "X", "alpha"))
    u = T.d1(grid, a, +1) / (2.0 * a) * connection
    delta = alpha * X * X * u
    phi_p = -0.5 * T.d1(grid, X, +1) / X
    q = 2.0 / x + T.d1(grid, b, +1) / b
    dM = T.d1(grid, delta, +1) + (6.0 * phi_p + 1.5 * q) * delta
    fit = alpha * X * X * (T.d1(grid, u, -1) + 3.0 * u / x)
    return delta, dM, fit


def locked_form_prediction(grid, state, vacuum):
    """Geometric part of the locked evaluator's predictor, transcribed (vacuum)."""
    d1 = grid.cell_derivative_fourth
    r = np.asarray(grid.centers, float)
    a, b, X, alpha, K = (np.asarray(getattr(state, k), float) for k in ("a", "b", "X", "alpha", "K"))
    lam_m = d1(a, parity=1) / (2.0 * a**2) - d1(b, parity=1) / (a * b) + 2.0 * (1.0 / b - 1.0 / a) / r
    metric = state.copy()
    metric.Lambda = lam_m.copy()
    mc = vacuum.constraints(grid, metric)
    H = np.asarray(mc["hamiltonian"], float)
    M = np.asarray(mc["momentum"], float)
    chi = X**2
    chi_r = -d1(X, parity=1) / (2.0 * X)
    C = np.asarray(state.Lambda, float) - lam_m
    W = chi * ((2.0 / 3.0) * d1(C, parity=-1) + (d1(a, parity=1) / (2.0 * a)) * C)
    F = alpha * W
    conn = d1(F, parity=1) + (6.0 * chi_r + 3.0 / r + 3.0 * d1(b, parity=1) / (2.0 * b)) * F
    return (alpha * K * M - (1.0 / 3.0) * d1(alpha, parity=1) * H
            + (1.0 / 6.0) * alpha * d1(H, parity=1) + conn)


def run(n):
    grid_ops, vacuum, moving = T.adapter.vendor_modules()
    grid = grid_ops.SphericalCellGrid(n, B.R_MAX)
    state = T.manufactured_state(grid, vacuum)
    x = np.asarray(grid.centers, float)
    con = vacuum.constraints(grid, state)
    rates = T.vacuum_rhs_rates(grid, state, vacuum, moving)
    measured = T.tangent_momentum(grid, state, rates, vacuum, EPS)
    _, predicted = B.prediction(grid, state, con, B.slot_bssn_minus_c(grid, vacuum, state))
    delta, dM, fit = derived_term(grid, state, np.asarray(con["connection"], float))
    removed = dict(rates)
    removed["Aa"] = rates["Aa"] - delta
    measured_removed = T.tangent_momentum(grid, state, removed, vacuum, EPS)
    before = measured - predicted
    locked = measured - locked_form_prediction(grid, state, vacuum)
    after = before - dM
    inner = slice(2, n - EDGE)
    sample = lambda f: {f"{t:g}": float(np.interp(t, x, f)) for t in FIXED_RADII}  # noqa: E731
    return {
        "N": n,
        "derived_term_max_abs": float(np.max(np.abs(dM))),
        "fit_form_difference_max_abs": float(np.max(np.abs(dM - fit))),
        "rate_removal_closure_max_abs": float(np.max(np.abs((measured - measured_removed) - dM))),
        "residual_before": {"cell_0": float(before[0]), "fixed_radius": sample(before),
                            "interior_max_abs": float(np.max(np.abs(before[inner])))},
        "residual_after": {"cell_0": float(after[0]), "cell_1": float(after[1]), "fixed_radius": sample(after),
                           "interior_max_abs": float(np.max(np.abs(after[inner])))},
        "measured_fixed_radius": sample(measured),
        "locked_form_residual": {"cell_0": float(locked[0]), "fixed_radius": sample(locked),
                                 "interior_max_abs": float(np.max(np.abs(locked[inner])))},
        "locked_form_minus_derived_max_abs": float(np.max(np.abs(locked - after))),
    }


def order(vals, ns):
    return [None if a == 0 or b == 0 else round(math.log(abs(a) / abs(b)) / math.log(m / k), 2)
            for a, b, k, m in zip(vals, vals[1:], ns, ns[1:])]


def main():
    rows = [run(n) for n in RES]
    rep = {
        "schema": "issue4_vendor_connection_term_derivation_v1", "read_only": True,
        "scientific_admission": False, "provenance": T.provenance(), "epsilon": EPS,
        "derived_term": "dM = D1(delta) + (6 phi' + 1.5 q) delta, delta = alpha X^2 (D1 a / 2a) C",
        "rows": rows,
        "orders_after": {
            "interior_max": order([r["residual_after"]["interior_max_abs"] for r in rows], RES),
            "cell_0": order([r["residual_after"]["cell_0"] for r in rows], RES),
            **{f"fixed_r_{t:g}": order([r["residual_after"]["fixed_radius"][f"{t:g}"] for r in rows], RES)
               for t in FIXED_RADII},
        },
        "orders_locked_form": {
            "interior_max": order([r["locked_form_residual"]["interior_max_abs"] for r in rows], RES),
            **{f"fixed_r_{t:g}": order([r["locked_form_residual"]["fixed_radius"][f"{t:g}"] for r in rows], RES)
               for t in FIXED_RADII},
        },
    }
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "results" / "issue4-hamiltonian-slot-decomposition-2026-11-08" / "vendor_connection_term_derivation.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rep, indent=1) + "\n")
    for r in rows:
        print(f"N={r['N']:>4} before r=2 {r['residual_before']['fixed_radius']['2']:+.2e}  after r=2 "
              f"{r['residual_after']['fixed_radius']['2']:+.2e}  after c0 {r['residual_after']['cell_0']:+.2e}  "
              f"after int {r['residual_after']['interior_max_abs']:.2e}  |dM-fit| {r['fit_form_difference_max_abs']:.1e}  "
              f"rate-removal {r['rate_removal_closure_max_abs']:.1e}  locked r=2 "
              f"{r['locked_form_residual']['fixed_radius']['2']:+.2e} locked int {r['locked_form_residual']['interior_max_abs']:.2e}")
    print("orders after:", rep["orders_after"])


if __name__ == "__main__":
    main()
