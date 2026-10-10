"""Diagnostic-only consistency audit for the legacy reference initial slice.

No trajectory is evolved, and no reference equation, source, or threshold is changed.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from . import reference_pirk_unified as ref

INNER_CELLS = 5
REFERENCE_BLOB = "a1636a2a26670a2c4dc9ba5f23985cd53c6b737d"


def _maxabs(values, start=0, stop=None):
    values = np.asarray(values, dtype=float)
    part = values[start:stop]
    if not np.all(np.isfinite(part)):
        return float("nan")
    return float(np.max(np.abs(part))) if part.size else 0.0


def _regions(values, r, n):
    inner_stop = min(INNER_CELLS, n)
    away_stop = max(inner_stop, n - 2)
    i25 = int(np.argmin(np.abs(r - 2.5)))
    return {
        "max_abs_all": _maxabs(values),
        "max_abs_first_five": _maxabs(values, 0, inner_stop),
        "max_abs_away_from_center": _maxabs(values, inner_stop, away_stop),
        "first_cell": float(values[0]),
        "at_r_near_2p5": float(values[i25]),
    }


def audit_initial_slice(resolution, r_max, amplitude=0.01, width=7.0):
    """Measure source balance and decompose the Hamiltonian mismatch."""
    if int(resolution) < 8 or not math.isfinite(float(r_max)) or r_max <= 0.0:
        raise ValueError("resolution must be >= 8 and r_max finite and positive")
    if not math.isfinite(float(amplitude)) or not math.isfinite(float(width)) or width <= 0.0:
        raise ValueError("amplitude must be finite and width finite and positive")

    grid = ref.Grid(int(resolution), float(r_max))
    state, fields = ref.make_initial(grid, Aamp=float(amplitude), width=float(width))
    r = np.asarray(grid.centers, dtype=float)
    S, PS, Df, PD, phi, Pi, rdm, rb = fields

    # Given the reference metric parametrization, b^3 = X^6 = B = A^(-2).
    B_b = np.asarray(state.b, dtype=float) ** 3
    B_X = np.asarray(state.X, dtype=float) ** 6
    B = B_b
    metric_identity = B_b - B_X

    # These analytic Gaussian gradients match those used by make_initial().
    Sp = -2.0 * r / float(width) ** 2 * S
    Dp = -2.0 * r / float(width) ** 2 * Df
    rho_source = (
        0.5 * (PS * PS + B * Sp * Sp) + 0.5 * S * S
        + 0.5 * (PD * PD + B * Dp * Dp) - 0.5 * Df * Df
        + (0.5 * Pi * Pi + ref.Vc(phi)) / ref.KAPPA
        + (rdm + rb) / ref.KAPPA
    )
    # The legacy H0 definition is intentionally repeated as an audit witness;
    # this is its own radiation-free, historically mapped reference background.
    H0 = math.sqrt((
        0.5 * 0.179055**2 + float(ref.Vc(np.array([33.8983]))[0])
        + 2.5857e-5 + 4.0306e-6
    ) / 3.0)
    Kt = np.full_like(r, -H0)
    Kr = Kt - 4.0 * math.pi * r * PD * Dp
    source = 1.0 - 8.0 * math.pi * r * r * rho_source + 2.0 * r * r * Kr * Kt + r * r * Kt * Kt
    dB = ref.D(grid, B, 1)

    # B + r B' = source is the continuum radial reconstruction identity.
    source_balance = B + r * dB - source
    R_source = 2.0 * (1.0 - source) / (r * r)
    R_polar = 2.0 * (1.0 - B) / (r * r) - 2.0 * dB / r

    # Independently evaluate the reference BSSN Ricci and constraint monitors.
    R_bssn = ref.ricci_terms(grid, state)[0]
    H_bssn, M_bssn, conn, det = ref.constraint(grid, state, fields)
    Aa = np.asarray(state.Aa)
    Ab = -0.5 * Aa
    rho_constraint = ref.matter_projection(grid, state, fields)[0]
    H_polar = R_polar - (Aa * Aa + 2.0 * Ab * Ab) + (2.0 / 3.0) * state.K * state.K - 16.0 * math.pi * rho_constraint

    # In source variables the Hamiltonian has the algebraic form
    # R_source + 4 K_r K_t + 2 K_t^2 - 16 pi rho = 0.
    H_source = R_source + 4.0 * Kr * Kt + 2.0 * Kt * Kt - 16.0 * math.pi * rho_source
    rho_mismatch = rho_constraint - rho_source
    R_polar_minus_source = R_polar - R_source
    R_bssn_minus_source = R_bssn - R_source

    # Identity closing the decomposition:
    # H_BSSN - H_source = (R_BSSN - R_source) - 16 pi (rho_constraint-rho_source).
    decomposition_closure = H_bssn - H_source - R_bssn_minus_source + 16.0 * math.pi * rho_mismatch

    arrays = (
        B, source, source_balance, R_source, R_bssn, R_polar,
        H_bssn, H_polar, H_source, rho_mismatch, decomposition_closure, M_bssn, conn, det
    )
    return {
        "resolution": int(resolution),
        "r_max": float(r_max),
        "dr": float(grid.dr),
        "amplitude": float(amplitude),
        "width": float(width),
        "H0_reference": H0,
        "metric_B_identity_max_abs_b3_minus_X6": _maxabs(metric_identity),
        "source_balance_B_plus_r_dB_minus_source": _regions(source_balance, r, int(resolution)),
        "ricci_BSSN_minus_polar": _regions(R_bssn - R_polar, r, int(resolution)),
        "ricci_polar_minus_source": _regions(R_polar_minus_source, r, int(resolution)),
        "ricci_BSSN_minus_source": _regions(R_bssn_minus_source, r, int(resolution)),
        "hamiltonian_BSSN": _regions(np.asarray(H_bssn), r, int(resolution)),
        "hamiltonian_independent_polar": _regions(np.asarray(H_polar), r, int(resolution)),
        "hamiltonian_source_algebraic_closure": _regions(np.asarray(H_source), r, int(resolution)),
        "matter_density_constraint_minus_source": _regions(rho_mismatch, r, int(resolution)),
        "hamiltonian_decomposition_closure": _regions(decomposition_closure, r, int(resolution)),
        "momentum_constraint_max_abs": _maxabs(M_bssn),
        "connection_constraint_max_abs": _maxabs(conn),
        "determinant_constraint_max_abs": _maxabs(det),
        "finite": bool(all(np.all(np.isfinite(np.asarray(v, dtype=float))) for v in arrays)),
        "interpretation": "Diagnostic witness only; this report does not approve a physical model or repair.",
    }


def run_audit():
    """Run cheap initial-slice controls: fixed domain and fixed dr=1 domain growth."""
    cases = []
    for n in (40, 80, 160, 320):
        for amp in (0.0, 0.01):
            cases.append(audit_initial_slice(n, 40.0, amplitude=amp))
    for n in (40, 80, 160, 320):
        for amp in (0.0, 0.01):
            cases.append(audit_initial_slice(n, float(n), amplitude=amp))
    return {
        "schema": "reference_initial_data_geometric_consistency_v1",
        "kind": "diagnostic_only_initial_slice_no_evolution",
        "reference_blob_sha": REFERENCE_BLOB,
        "case_count": len(cases),
        "all_finite": all(case["finite"] for case in cases),
        "cases": cases,
        "acceptance_policy": "No physical admission or repair acceptance is claimed by this diagnostic.",
    }


def main(output="runs/reference-initial-data-geometric-consistency/report.json"):
    result = run_audit()
    output_path = Path(output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(result, indent=2, allow_nan=False))
    compact = []
    for case in result["cases"]:
        compact.append({
            "N": case["resolution"], "r_max": case["r_max"], "dr": case["dr"],
            "S_amplitude": case["amplitude"],
            "source_first5": case["source_balance_B_plus_r_dB_minus_source"]["max_abs_first_five"],
            "source_away": case["source_balance_B_plus_r_dB_minus_source"]["max_abs_away_from_center"],
            "R_BSSN_minus_polar_first": case["ricci_BSSN_minus_polar"]["first_cell"],
            "R_polar_minus_source_first": case["ricci_polar_minus_source"]["first_cell"],
            "R_BSSN_minus_source_first": case["ricci_BSSN_minus_source"]["first_cell"],
            "H_source_closure_first": case["hamiltonian_source_algebraic_closure"]["first_cell"],
            "rho_density_mismatch_first": case["matter_density_constraint_minus_source"]["first_cell"],
            "H_bssn_first": case["hamiltonian_BSSN"]["first_cell"],
            "H_polar_first": case["hamiltonian_independent_polar"]["first_cell"],
            "H_bssn_r2p5": case["hamiltonian_BSSN"]["at_r_near_2p5"],
        })
    print(json.dumps({
        "schema": result["schema"], "case_count": result["case_count"],
        "all_finite": result["all_finite"], "cases": compact,
    }, indent=2), flush=True)
    return 0 if result["all_finite"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
