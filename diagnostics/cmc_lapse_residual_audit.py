"""Same-slice diagnostic for the stage-aware CMC lapse solve.

This reproduces the exact discrete tridiagonal operator used by
engine.cmc_gauge.solve_cmc_lapse. It does not edit production equations.
The lapse is repeatedly re-solved on a frozen geometry/matter/scalar slice,
with only geometry.alpha updated between calls, so target_kdot's dependence
on the incoming lapse can be measured separately from linear-solver error.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from engine.cmc_gauge import solve_cmc_lapse, target_kdot
from engine.production_kernel import V55ProductionKernel
from engine.v55_matter import total_matter_projection

OUT = ROOT / "runs" / "cmc_lapse_residual_audit"
ITERATIONS = 8


def exact_discrete_residual(grid, geometry, scalars, matter, alpha, kdot=None):
    """Return residual of the exact discrete equation assembled in cmc_gauge.py.

    Interior equation:
      alpha'' + c alpha' - (Q/inv) alpha = -kdot/inv
    with the same centered second-order lapse stencil, fourth-order metric
    derivatives, regular-center row alpha[0]-alpha[1]=0, and alpha[-1]=1.
    """
    a = np.asarray(geometry.a, dtype=float)
    b = np.asarray(geometry.b, dtype=float)
    X = np.asarray(geometry.X, dtype=float)
    r = np.asarray(grid.centers, dtype=float)
    alpha = np.asarray(alpha, dtype=float)
    inv = X * X / a

    ap = grid.cell_derivative_fourth(a, parity=1)
    bp = grid.cell_derivative_fourth(b, parity=1)
    Xp = grid.cell_derivative_fourth(X, parity=1)
    c = -0.5 * ap / a + bp / b - Xp / X + 2.0 / r

    total = total_matter_projection(grid, geometry, scalars, matter)
    Q = (
        1.5 * geometry.Aa * geometry.Aa
        + geometry.K * geometry.K / 3.0
        + 4.0 * np.pi * (total["rho"] + total["pr"] + 2.0 * total["pt"])
    )
    m = Q / inv
    if kdot is None:
        kdot = target_kdot(grid, geometry, scalars, matter)
    rhs = -float(kdot) / inv

    h = float(grid.dr)
    residual = np.empty(grid.n, dtype=float)
    residual[0] = alpha[0] - alpha[1]
    d2 = (alpha[2:] - 2.0 * alpha[1:-1] + alpha[:-2]) / (h * h)
    d1 = (alpha[2:] - alpha[:-2]) / (2.0 * h)
    lhs_interior = d2 + c[1:-1] * d1 - m[1:-1] * alpha[1:-1]
    residual[1:-1] = lhs_interior - rhs[1:-1]
    residual[-1] = alpha[-1] - 1.0

    scale = max(
        float(np.max(np.abs(lhs_interior))),
        float(np.max(np.abs(rhs[1:-1]))),
        float(np.max(np.abs(d2))),
        float(np.max(np.abs(c[1:-1] * d1))),
        float(np.max(np.abs(m[1:-1] * alpha[1:-1]))),
        1.0,
    )
    return residual, float(kdot), scale


def rms(values):
    values = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(values * values)))


def profile_summary(grid, values):
    values = np.asarray(values, dtype=float)
    idx = int(np.argmax(np.abs(values)))
    return {
        "rms": rms(values),
        "max_abs": float(np.max(np.abs(values))),
        "center": float(values[0]),
        "outer": float(values[-1]),
        "max_abs_index": idx,
        "max_abs_radius": float(grid.centers[idx]),
        "max_abs_signed_value": float(values[idx]),
    }


def cosine_similarity(left, right):
    left = np.asarray(left, dtype=float)
    right = np.asarray(right, dtype=float)
    denom = float(np.linalg.norm(left) * np.linalg.norm(right))
    if denom == 0.0:
        return None
    return float(np.dot(left, right) / denom)


def main():
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=64,
        r_max=32.0,
        amplitude=1.0e-3,
        width=7.0,
        D_amplitude=0.0,
        include_radiation=False,
    )

    # Freeze every field other than alpha throughout this audit. A copy is
    # used at each pass to make the only changing input explicit.
    frozen = state.geometry.copy()
    previous_delta = None
    rows = []
    finite = True
    solver_reconstruction_ok = True

    for iteration in range(ITERATIONS):
        geometry_in = frozen.copy()
        if rows:
            geometry_in.alpha = rows[-1]["_alpha_out_array"].copy()

        alpha_in = np.asarray(geometry_in.alpha, dtype=float).copy()
        residual_in, kdot_in, scale_in = exact_discrete_residual(
            state.grid, geometry_in, state.scalars, state.matter, alpha_in
        )

        alpha_out, kdot_used_by_solver = solve_cmc_lapse(
            state.grid, geometry_in, state.scalars, state.matter
        )
        alpha_out = np.asarray(alpha_out, dtype=float).copy()

        # Check that our reconstructed discrete operator matches the linear
        # system actually solved, using the target computed from the input.
        geometry_out = geometry_in.copy()
        geometry_out.alpha = alpha_out.copy()
        solver_residual, _, solver_scale = exact_discrete_residual(
            state.grid, geometry_out, state.scalars, state.matter, alpha_out,
            kdot=float(kdot_used_by_solver),
        )

        # Re-evaluate the CMC target on the returned lapse. This is the
        # same-slice self-consistency residual: same non-lapse state, but the
        # target is now computed from the output lapse rather than the input.
        residual_out, kdot_out, scale_out = exact_discrete_residual(
            state.grid, geometry_out, state.scalars, state.matter, alpha_out
        )

        delta = alpha_out - alpha_in
        row = {
            "iteration": iteration,
            "kdot_input": float(kdot_in),
            "kdot_used_by_solver": float(kdot_used_by_solver),
            "kdot_recomputed_on_output": float(kdot_out),
            "kdot_output_minus_input": float(kdot_out - kdot_in),
            "input_lapse_residual": profile_summary(state.grid, residual_in),
            "input_lapse_residual_relative_max": float(
                np.max(np.abs(residual_in)) / scale_in
            ),
            "linear_solve_reconstruction_residual": profile_summary(
                state.grid, solver_residual
            ),
            "linear_solve_reconstruction_relative_max": float(
                np.max(np.abs(solver_residual)) / solver_scale
            ),
            "output_lapse_residual_using_recomputed_target": profile_summary(
                state.grid, residual_out
            ),
            "output_lapse_residual_relative_max": float(
                np.max(np.abs(residual_out)) / scale_out
            ),
            "lapse_shift": profile_summary(state.grid, delta),
            "lapse_shift_vs_previous_cosine": (
                None if previous_delta is None
                else cosine_similarity(delta, previous_delta)
            ),
            "lapse_shift_change_from_previous_rms": (
                None if previous_delta is None else rms(delta - previous_delta)
            ),
            "_alpha_out_array": alpha_out,
        }
        rows.append(row)
        previous_delta = delta.copy()

        finite = finite and all(
            np.all(np.isfinite(np.asarray(values, dtype=float)))
            for values in (alpha_in, alpha_out, residual_in, solver_residual, residual_out)
        )
        solver_reconstruction_ok = solver_reconstruction_ok and (
            row["linear_solve_reconstruction_relative_max"] <= 1.0e-10
        )

    for row in rows:
        row.pop("_alpha_out_array", None)

    shifts = [row["lapse_shift"]["rms"] for row in rows]
    shift_change = [
        row["lapse_shift_change_from_previous_rms"]
        for row in rows[1:]
    ]
    payload = {
        "test": "same-slice CMC lapse discrete residual and repeated-solve audit",
        "scope": (
            "actual V55ProductionKernel initialized at N=64, r_max=32, "
            "S amplitude=1e-3, D amplitude=0, radiation off; non-lapse "
            "geometry, scalars, and matter frozen while each returned lapse "
            "is fed into the next CMC target calculation"
        ),
        "production_equations_changed": False,
        "iterations": rows,
        "summary": {
            "lapse_shift_rms_by_iteration": shifts,
            "lapse_shift_rms_range": [
                float(min(shifts)), float(max(shifts))
            ],
            "max_shift_rms_change_between_calls": (
                float(max(shift_change)) if shift_change else None
            ),
            "max_linear_solve_reconstruction_relative_residual": float(
                max(row["linear_solve_reconstruction_relative_max"] for row in rows)
            ),
            "max_input_self_slice_residual_relative": float(
                max(row["input_lapse_residual_relative_max"] for row in rows)
            ),
            "max_output_self_slice_residual_relative": float(
                max(row["output_lapse_residual_relative_max"] for row in rows)
            ),
            "max_abs_kdot_change_per_solve": float(max(
                abs(row["kdot_output_minus_input"]) for row in rows
            )),
        },
        "checks": {
            "all_values_finite": bool(finite),
            "reconstructed_linear_system_matches_solver": bool(solver_reconstruction_ok),
        },
        "pass": bool(finite and solver_reconstruction_ok),
        "interpretation_limit": (
            "The output self-slice residual measures failure of one solve to "
            "be a fixed point when target_kdot is recomputed from the returned "
            "lapse. It does not, by itself, establish the cause of that "
            "target dependence or explain strong-D proper-time differences."
        ),
    }

    OUT.mkdir(parents=True, exist_ok=True)
    result_path = OUT / "cmc_lapse_residual_audit.json"
    result_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload, indent=2, sort_keys=True), flush=True)
    print(f"RESULT_JSON={result_path}", flush=True)
    if not payload["pass"]:
        raise SystemExit("FAIL: non-finite output or discrete solve reconstruction mismatch")
    print("PASS: residual audit completed; physical interpretation remains diagnostic-only")


if __name__ == "__main__":
    main()
