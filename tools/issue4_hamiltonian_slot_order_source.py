"""Read-only: which discrete construction limits H_BSSN - H_metric to second order?

The manufactured state is analytic, so the exact continuum Hamiltonian H_exact and the
exact metric connection Lambda_metric_exact can be built symbolically. Each discrete
piece is then measured against them:

  E_metric = H_metric(discrete)                     - H_exact
  E_bssn   = H_BSSN(vendor, Lambda=Lambda_exact)     - H_exact
  E_conn   = C(vendor) with Lambda=Lambda_exact      (exact value 0)
  E_lam    = Lambda_metric(d1, as in the tool)       - Lambda_metric_exact
  E_slot   = [H_BSSN - X^2 D_r(C)](vendor, as manufactured) - H_exact   (candidate slot)

The candidate slot's response to the tool's Lambda-only mutation is also reported.

Usage: COSMOS_PINNED_SOURCE_ROOT=<checkout of 7f47669a with vendor submodule> \
       python tools/issue4_hamiltonian_slot_order_source.py [output.json]

Interior = cells 2..n-5; cells 0, 1 and the outer 4 cells reported separately.
"""
import json
import math
import sys
from pathlib import Path

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import tools.audit_brown_vacuum_momentum_propagation as T  # noqa: E402

RES = [40, 60, 80, 120, 160, 240, 320]
R_MAX, EDGE = 40.0, 4

r = sp.symbols("r", positive=True)
q = (r / 6) ** 2
e = sp.exp(-q)
f = sp.Rational(2, 1000) * q * e
a_s, b_s = sp.exp(f), sp.exp(-f / 2)
X_s = 1 + sp.Rational(15, 10000) * e
Aa_s = sp.Rational(2, 100000) * r**2 * e
K_s = -sp.Rational(2, 100) + sp.Rational(2, 1000) * e
Rar = r * sp.sqrt(b_s) / X_s                    # areal radius
pf = sp.sqrt(a_s) / X_s                          # dl/dr
R_l = sp.diff(Rar, r) / pf
R_ll = sp.diff(R_l, r) / pf
Ric = -2 * R_ll / Rar + 2 * (1 - R_l**2 - Rar * R_ll) / Rar**2
H_exact_s = Ric - (Aa_s**2 + 2 * (Aa_s / 2) ** 2) + sp.Rational(2, 3) * K_s**2
lam_exact_s = (sp.diff(a_s, r) / (2 * a_s**2) - sp.diff(b_s, r) / (a_s * b_s)
               + 2 / r * (1 / b_s - 1 / a_s))
H_exact = sp.lambdify(r, H_exact_s, "mpmath")
lam_exact = sp.lambdify(r, lam_exact_s, "mpmath")


def ev(fn, xs):
    import mpmath as mp
    mp.mp.dps = 40
    return np.array([float(fn(mp.mpf(float(x)))) for x in xs])


FIXED_RADII = (1.0, 2.0, 5.0, 10.0)


def summarize(err, n, x):
    inner = err[2:n - EDGE]
    i = int(np.argmax(np.abs(inner))) + 2
    return {"interior_max_abs": float(np.max(np.abs(inner))), "interior_argmax_cell": i,
            "cell_0": float(err[0]), "cell_1": float(err[1]),
            "outer_edge_max_abs": float(np.max(np.abs(err[n - EDGE:]))),
            "fixed_radius": {f"{t:g}": float(np.interp(t, x, err)) for t in FIXED_RADII}}


def slot(grid, vacuum, state):
    c = vacuum.constraints(grid, state)
    X = np.asarray(state.X, float)
    return np.asarray(c["hamiltonian"], float) - X * X * T.d1(grid, np.asarray(c["connection"], float), -1)


def run(n):
    grid_ops, vacuum, _ = T.adapter.vendor_modules()
    grid = grid_ops.SphericalCellGrid(n, R_MAX)
    state = T.manufactured_state(grid, vacuum)
    x = np.asarray(grid.centers, float)
    Hx, Lx = ev(H_exact, x), ev(lam_exact, x)
    # sanity: the tool's analytic fields equal the symbolic ones
    for name, s in (("a", a_s), ("b", b_s), ("X", X_s), ("Aa", Aa_s), ("K", K_s)):
        dev = np.max(np.abs(np.asarray(getattr(state, name)) - ev(sp.lambdify(r, s, "mpmath"), x)))
        assert dev < 1e-14, (name, dev)
    H_m, _ = T.metric_only_hamiltonian(grid, state)
    lam_tool = np.asarray(state.Lambda) - 2.5e-4 * x / 36.0 * np.exp(-(x / 6.0) ** 2)
    exact_state = state.copy()
    exact_state.Lambda = Lx.copy()
    con = vacuum.constraints(grid, exact_state)
    H_slot = slot(grid, vacuum, state)
    probe = state.copy()
    probe.Lambda = np.asarray(state.Lambda) + 1.0e-3 * np.sin((2.0 * math.pi / 8.0) * x)
    return {
        "E_slot_candidate": summarize(H_slot - Hx, n, x),
        "slot_candidate_lambda_mutation_max_abs": float(np.max(np.abs(slot(grid, vacuum, probe) - H_slot))),
        "N": n, "max_abs_H_exact": float(np.max(np.abs(Hx))),
        "E_metric": summarize(H_m - Hx, n, x),
        "E_bssn_exact_lambda": summarize(np.asarray(con["hamiltonian"]) - Hx, n, x),
        "E_conn_vendor_exact_lambda": summarize(np.asarray(con["connection"]), n, x),
        "E_lambda_metric_tool": summarize(lam_tool - Lx, n, x),
    }


def order(rows, key, field="interior_max_abs"):
    out = []
    for p, s in zip(rows, rows[1:]):
        e1, e2 = abs(p[key][field]), abs(s[key][field])
        out.append(None if e1 == 0 or e2 == 0 else round(float(np.log(e1 / e2) / np.log(s["N"] / p["N"])), 2))
    return out


def order_fixed(rows, key, t):
    out = []
    for p, s in zip(rows, rows[1:]):
        e1, e2 = abs(p[key]["fixed_radius"][t]), abs(s[key]["fixed_radius"][t])
        out.append(None if e1 == 0 or e2 == 0 else round(float(np.log(e1 / e2) / np.log(s["N"] / p["N"])), 2))
    return out


def main():
    rows = [run(n) for n in RES]
    keys = ["E_metric", "E_bssn_exact_lambda", "E_conn_vendor_exact_lambda", "E_lambda_metric_tool",
            "E_slot_candidate"]
    rep = {"schema": "issue4_hamiltonian_slot_order_source_v1", "read_only": True,
           "provenance": T.provenance(), "rows": rows,
           "orders_interior": {k: order(rows, k) for k in keys},
           "orders_cell_0": {k: order(rows, k, "cell_0") for k in keys},
           "orders_fixed_radius": {k: {f"{t:g}": order_fixed(rows, k, f"{t:g}") for t in FIXED_RADII}
                                   for k in keys}}
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "results" / "issue4-hamiltonian-slot-decomposition-2026-11-08" / "order_source.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rep, indent=1) + "\n")
    for k in keys:
        print(f"== {k}")
        for row in rows:
            s = row[k]
            print(f"  N={row['N']:>3} interior={s['interior_max_abs']:.3e} @cell {s['interior_argmax_cell']:>3}"
                  f"  c0={s['cell_0']:+.3e} c1={s['cell_1']:+.3e} outer={s['outer_edge_max_abs']:.2e}")
        print("  order interior:", rep["orders_interior"][k], " cell0:", rep["orders_cell_0"][k])
    for k in keys:
        print("fixed-radius orders", k, rep["orders_fixed_radius"][k])
    print("slot Lambda-mutation response:", [f"{row['slot_candidate_lambda_mutation_max_abs']:.1e}" for row in rows])


if __name__ == "__main__":
    main()
