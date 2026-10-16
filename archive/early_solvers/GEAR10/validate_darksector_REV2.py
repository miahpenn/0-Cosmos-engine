"""Bounded archive preflight for the reconstructed GEAR-10 dark-sector solver.

This script intentionally does not claim to reproduce the complete six-rung ladder.
It checks import/RHS behavior, potential derivatives, the two corrected Euler signs,
and the general R6 momentum identity with trapezoidal convergence.
"""
from __future__ import annotations

import hashlib
import json
import platform
from pathlib import Path

import numpy as np
import scipy

import darksector_REV2 as solver

HERE = Path(__file__).resolve().parent
SOURCE = HERE / "darksector_REV2.py"
RESULTS = HERE / "validation_results_r6_sign_discriminator_2026-10-16.json"


def check_model_derivatives(model: solver.Model) -> dict:
    max_dv_error = 0.0
    max_d2v_error = 0.0
    h = 1.0e-5
    for phi in (-0.3, 0.0, 0.4):
        dV_fd = (model.V(phi + h) - model.V(phi - h)) / (2 * h)
        d2V_fd = (model.dV(phi + h) - model.dV(phi - h)) / (2 * h)
        max_dv_error = max(max_dv_error, abs(dV_fd - model.dV(phi)))
        max_d2v_error = max(max_d2v_error, abs(d2V_fd - model.d2V(phi)))
    assert max_dv_error < 2.0e-8, max_dv_error
    assert max_d2v_error < 2.0e-8, max_d2v_error
    return {
        "max_abs_dV_finite_difference_error": max_dv_error,
        "max_abs_d2V_finite_difference_error": max_d2v_error,
        "pass": True,
    }


def check_rhs_and_corrected_euler_signs(model: solver.Model) -> dict:
    y_bg = [0.3, 0.1, 0.02, 0.2, 0.1]
    background = solver.background_rhs(0.0, y_bg, model)
    assert len(background) == 5 and np.all(np.isfinite(background))
    a, phi, pp, rdm, rb = y_bg
    H = a * np.sqrt((rdm + rb + pp**2 / (2 * a**2) + model.V(phi)) / 3.0)
    y_full = y_bg + [0.01, 0.03, -0.02, 0.04, 0.005, 0.007, 0.001]
    k = 10.0
    full = solver.full_rhs(0.0, y_full, model, k)
    assert len(full) == 12 and np.all(np.isfinite(full))
    expected_tdm = (
        -H * y_full[6]
        + k * k * y_full[11]
        - model.beta * k * k * y_full[9]
        + model.beta * y_full[2] * y_full[6]
    )
    error = abs(full[6] - expected_tdm)
    assert error < 1.0e-14, error
    return {"background_rhs_length": len(background), "full_rhs_length": len(full),
            "euler_rhs_abs_error": error, "pass": True}


def _identity_metrics(out: dict, k: float, beta: float) -> dict:
    """Return the R6 balance and both historical normalization conventions."""
    a = out["a"]
    phi = out["phi"]
    theta = out["tdm"]
    delta_phi = out["dphi"]
    Phi = out["Phi"]
    tau = out["tau"]
    mass = np.exp(-beta * phi)
    lhs = float(a[-1] * mass[-1] * theta[-1] - a[0] * mass[0] * theta[0])
    integrand = a * mass * (k * k * Phi - beta * k * k * delta_phi)
    # Direct composite trapezoid sum: compatible with NumPy releases both before and after np.trapz removal.
    rhs_integral = float(np.sum(0.5 * (integrand[:-1] + integrand[1:]) * np.diff(tau)))
    abs_gap = abs(lhs - rhs_integral)
    max_side_denominator = max(abs(lhs), abs(rhs_integral), 1.0e-30)
    sum_sides_denominator = max(abs(lhs) + abs(rhs_integral), 1.0e-30)
    residual_max = abs_gap / max_side_denominator
    residual_sum = abs_gap / sum_sides_denominator
    return {
        "lhs_delta_a_m_theta": lhs,
        "rhs_integral": rhs_integral,
        "absolute_gap": abs_gap,
        "normalized_residual": residual_max,  # compatibility alias: max-side normalization
        "normalized_residual_max_side": residual_max,
        "normalized_residual_sum_sides": residual_sum,
    }


def _evolve_rung6(model: solver.Model, k: float, n: int, flip_drag_sign: bool = False) -> dict:
    """Evolve normally, or apply a test-only fault to the DM Euler drag term.

    The fault reverses only +beta*phi_prime*theta_dm in theta_dm'. It is an
    in-memory monkeypatch for this diagnostic call and never edits the archive
    solver file. The module global is always restored, including on failure.
    """
    if not flip_drag_sign:
        return solver.evolve(model, k=k, a_i=1.0e-4, a_f=1.0, n=n)

    original_rhs = solver.full_rhs

    def wrong_drag_rhs(tau, y, M, k_value):
        dy = list(original_rhs(tau, y, M, k_value))
        # Correct term is +beta*phi_prime*theta_dm; sign reversal changes it by -2x.
        dy[6] -= 2.0 * M.beta * y[2] * y[6]
        return dy

    solver.full_rhs = wrong_drag_rhs
    try:
        return solver.evolve(model, k=k, a_i=1.0e-4, a_f=1.0, n=n)
    finally:
        solver.full_rhs = original_rhs


def _compute_rung6_rows(
    model: solver.Model,
    k: float,
    sample_counts: tuple[int, ...] = (5000, 20000, 80000),
    flip_drag_sign: bool = False,
) -> list[dict]:
    rows = []
    for n in sample_counts:
        out = _evolve_rung6(model, k=k, n=n, flip_drag_sign=flip_drag_sign)
        rows.append({"n": n, **_identity_metrics(out, k=k, beta=model.beta)})

    for i, row in enumerate(rows):
        for suffix, key in (
            ("max_side", "normalized_residual_max_side"),
            ("sum_sides", "normalized_residual_sum_sides"),
        ):
            order = None
            if i > 0:
                previous = rows[i - 1]
                current_value = row[key]
                previous_value = previous[key]
                if current_value > 0.0 and previous_value > 0.0:
                    order = float(
                        np.log(previous_value / current_value)
                        / np.log((row["n"] - 1) / (previous["n"] - 1))
                    )
            row[f"observed_order_{suffix}_from_previous"] = order
        row["observed_order_from_previous"] = row["observed_order_max_side_from_previous"]
    return rows


def check_general_rung6(model: solver.Model, k: float = 1.0) -> list[dict]:
    """Backward-compatible single-case check, now reporting both normalizers."""
    rows = _compute_rung6_rows(model, k=k)
    max_side = [row["normalized_residual_max_side"] for row in rows]
    assert max_side[0] > max_side[1] > max_side[2], max_side
    orders = [row["observed_order_max_side_from_previous"] for row in rows[1:]]
    assert all(order is not None for order in orders), orders
    assert min(orders) > 1.8, orders
    return rows


def check_general_rung6_matrix() -> list[dict]:
    """Reproduce the four archived general-form R6 (k, beta) cases."""
    cases = ((1.0, 0.0), (1.0, -0.20), (100.0, -0.04), (1000.0, -0.20))
    matrix = []
    for k, beta in cases:
        model = solver.Model(beta=beta, lam=1.0, D=0.0, V0=2.07)
        rows = _compute_rung6_rows(model, k=k)
        matrix.append({
            "k": k,
            "beta": beta,
            "case": f"k={k:g}, beta={beta:g}",
            "normalization_conventions": {
                "max_side": "abs(lhs-rhs) / max(abs(lhs), abs(rhs))",
                "sum_sides": "abs(lhs-rhs) / (abs(lhs) + abs(rhs))",
            },
            "rows": rows,
        })
    return matrix


def check_drag_sign_fault_injection(matrix: list[dict]) -> list[dict]:
    """Pair corrected and wrong-drag R6 balances without touching solver source."""
    selected = ((1.0, -0.20), (100.0, -0.04))
    results = []
    for k, beta in selected:
        corrected_case = next(
            item for item in matrix if item["k"] == k and item["beta"] == beta
        )
        corrected_rows = corrected_case["rows"]
        model = solver.Model(beta=beta, lam=1.0, D=0.0, V0=2.07)
        fault_rows = _compute_rung6_rows(model, k=k, flip_drag_sign=True)
        paired = []
        for correct, faulty in zip(corrected_rows, fault_rows):
            paired.append({
                "n": correct["n"],
                "corrected": {
                    "lhs_delta_a_m_theta": correct["lhs_delta_a_m_theta"],
                    "rhs_integral": correct["rhs_integral"],
                    "normalized_residual_max_side": correct["normalized_residual_max_side"],
                    "normalized_residual_sum_sides": correct["normalized_residual_sum_sides"],
                },
                "drag_sign_flipped_test_only": {
                    "lhs_delta_a_m_theta": faulty["lhs_delta_a_m_theta"],
                    "rhs_integral": faulty["rhs_integral"],
                    "normalized_residual_max_side": faulty["normalized_residual_max_side"],
                    "normalized_residual_sum_sides": faulty["normalized_residual_sum_sides"],
                },
            })
        correct_residuals = [row["corrected"]["normalized_residual_max_side"] for row in paired]
        fault_residuals = [row["drag_sign_flipped_test_only"]["normalized_residual_max_side"] for row in paired]
        assert correct_residuals[0] > correct_residuals[1] > correct_residuals[2], correct_residuals
        assert correct_residuals[-1] < correct_residuals[0] / 100.0, correct_residuals
        assert min(fault_residuals) > 1.0e-3, fault_residuals
        plateau_ratio = fault_residuals[-1] / fault_residuals[0]
        assert 0.5 < plateau_ratio < 2.0, fault_residuals
        results.append({
            "k": k,
            "beta": beta,
            "fault_scope": "test-only reversal of beta*phi_prime*theta_dm in theta_dm'; production/transcription source unchanged",
            "rows": paired,
            "corrected_max_side_refinement_factor": correct_residuals[0] / correct_residuals[-1],
            "fault_max_side_last_to_first_ratio": plateau_ratio,
            "fault_discriminator_pass": True,
        })
    return results

def main() -> None:
    model = solver.Model(beta=-0.2, lam=1.0, D=0.0, V0=2.07)
    matrix = check_general_rung6_matrix()
    baseline_case = next(
        item for item in matrix if item["k"] == 1.0 and item["beta"] == -0.2
    )
    baseline_rows = baseline_case["rows"]
    fault_injection = check_drag_sign_fault_injection(matrix)
    payload = {
        "test_date": "2026-10-16",
        "artifact_scope": "research-only transcription from corrected REV2 PDF; bounded R6 matrix and test-only drag-sign fault injection; not full six-rung admission",
        "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "model": {
            "beta": model.beta, "lambda": model.lam, "D": model.D,
            "V0": model.V0, "k_rung6": 1.0, "a_i": 1.0e-4, "a_f": 1.0,
        },
        "potential_derivatives": check_model_derivatives(model),
        "rhs_and_corrected_euler_signs": check_rhs_and_corrected_euler_signs(model),
        "general_rung6": baseline_rows,
        "rung6_general_form_matrix": matrix,
        "rung6_drag_sign_fault_injection": fault_injection,
        "historical_normalization_note": (
            "The older verification text says max-side normalization, while the archived R6 table "
            "numerically matches abs(lhs-rhs)/(abs(lhs)+abs(rhs)) more closely. Both conventions "
            "are reported here; this provenance mismatch is preserved, not silently adjudicated."
        ),
        "full_six_rung_ladder_reproduced": False,
        "camb_rung2_independently_reproduced": False,
        "historical_pre_correction_source_bytes_recovered": False,
        "scientific_admission": False,
    }
    RESULTS.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
