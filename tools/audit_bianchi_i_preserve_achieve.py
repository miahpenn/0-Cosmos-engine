"""Research-only Bianchi-I test oracle for the preserve/achieve shear identity.

This module deliberately does not import or modify the production spherical engine.
It tests the restricted homogeneous Bianchi-I submodel with unsourced shear and a
separately conserved isotropic fluid. A passing result validates this bounded
numerical identity only; it is not a production-trajectory or scientific-admission
result.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from scipy.integrate import solve_ivp


def slope_domega_dln_a(omega_sigma: float, w: float) -> float:
    """Return d(Omega_sigma)/d(ln a) for the unsourced Bianchi-I submodel."""
    if not (math.isfinite(omega_sigma) and math.isfinite(w)):
        raise ValueError("omega_sigma and w must be finite")
    if not 0.0 <= omega_sigma <= 1.0:
        raise ValueError("omega_sigma must lie in [0, 1]")
    return 3.0 * (w - 1.0) * omega_sigma * (1.0 - omega_sigma)


def slope_domega_dcontraction_efold(omega_sigma: float, w: float) -> float:
    """Return d(Omega_sigma)/dN_c, N_c increasing during contraction."""
    return -slope_domega_dln_a(omega_sigma, w)


def slope_with_sources_dln_a(
    *,
    sigma2: float,
    rho: float,
    w: float,
    H: float,
    sigma2_source: float = 0.0,
    rho_exchange: float = 0.0,
    min_abs_H: float = 0.0,
    constraint_rtol: float = 1.0e-10,
) -> float:
    """Return the source-aware d(Omega_sigma)/d(ln a) at a valid nonzero-H sample.

    Conventions:
      3 H^2 = rho + sigma2
      dot(sigma2) = -6 H sigma2 + sigma2_source
      dot(rho) = -3 H (1+w) rho + rho_exchange

    The caller must provide source terms derived from its chosen model. This helper
    does not infer them or silently approximate a turnaround sample.
    """
    values = (sigma2, rho, w, H, sigma2_source, rho_exchange, min_abs_H, constraint_rtol)
    if not all(math.isfinite(value) for value in values):
        raise ValueError("all inputs must be finite")
    if sigma2 <= 0.0 or rho <= 0.0:
        raise ValueError("sigma2 and rho must be strictly positive for source ratios")
    if min_abs_H < 0.0 or constraint_rtol < 0.0:
        raise ValueError("min_abs_H and constraint_rtol must be non-negative")
    if abs(H) <= min_abs_H:
        raise ValueError("H is zero or inside the caller-declared turnaround exclusion band")

    friedmann_left = 3.0 * H * H
    friedmann_right = rho + sigma2
    scale = max(abs(friedmann_left), abs(friedmann_right), 1.0)
    if abs(friedmann_left - friedmann_right) > constraint_rtol * scale:
        raise ValueError("sample fails the declared 3 H^2 = rho + sigma2 constraint")

    omega_sigma = sigma2 / friedmann_left
    bracket = (
        3.0 * (w - 1.0)
        + sigma2_source / (H * sigma2)
        - rho_exchange / (H * rho)
    )
    return omega_sigma * (1.0 - omega_sigma) * bracket


def run_unsourced_case(
    *,
    w: float,
    omega_sigma_initial: float = 0.2,
    n_contraction_end: float = 2.0,
    samples: int = 801,
    rtol: float = 1.0e-11,
    atol: float = 1.0e-13,
) -> dict[str, Any]:
    """Integrate the source-free density scalings and compare to exact solutions.

    N_c = ln(a_ref / a) increases with physical time on a contracting branch.
    The independent numerical system is y' = (6 sigma2, 3(1+w) rho).
    """
    if not math.isfinite(w):
        raise ValueError("w must be finite")
    if not 0.0 < omega_sigma_initial < 1.0:
        raise ValueError("omega_sigma_initial must lie strictly between 0 and 1")
    if not math.isfinite(n_contraction_end) or n_contraction_end <= 0.0:
        raise ValueError("n_contraction_end must be finite and positive")
    if samples < 5:
        raise ValueError("at least five samples are required for interior slope checks")
    if not (math.isfinite(rtol) and math.isfinite(atol) and rtol > 0.0 and atol > 0.0):
        raise ValueError("rtol and atol must be finite and positive")

    nc = np.linspace(0.0, n_contraction_end, samples)
    sigma2_initial = omega_sigma_initial
    rho_initial = 1.0 - omega_sigma_initial

    def rhs(_nc: float, y: np.ndarray) -> tuple[float, float]:
        sigma2, rho = float(y[0]), float(y[1])
        return 6.0 * sigma2, 3.0 * (1.0 + w) * rho

    numerical = solve_ivp(
        rhs,
        (float(nc[0]), float(nc[-1])),
        (sigma2_initial, rho_initial),
        t_eval=nc,
        method="DOP853",
        rtol=rtol,
        atol=atol,
    )
    if not numerical.success:
        raise RuntimeError(f"bounded reference integration failed: {numerical.message}")

    sigma2_num, rho_num = numerical.y
    if not (np.all(np.isfinite(sigma2_num)) and np.all(np.isfinite(rho_num))):
        raise RuntimeError("reference integration produced non-finite component values")
    if np.any(sigma2_num <= 0.0) or np.any(rho_num <= 0.0):
        raise RuntimeError("reference integration left the strictly positive component domain")

    sigma2_exact = sigma2_initial * np.exp(6.0 * nc)
    rho_exact = rho_initial * np.exp(3.0 * (1.0 + w) * nc)
    omega_num = sigma2_num / (sigma2_num + rho_num)
    omega_exact = sigma2_exact / (sigma2_exact + rho_exact)

    # The numerical derivative is with respect to N_c; hence it must have the
    # opposite sign to dOmega/dln(a) during contraction.
    finite_difference = np.gradient(omega_num, nc, edge_order=2)
    analytic_contraction_slope = 3.0 * (1.0 - w) * omega_exact * (1.0 - omega_exact)
    interior = slice(2, -2)
    max_component_relative_error = max(
        float(np.max(np.abs(sigma2_num - sigma2_exact) / sigma2_exact)),
        float(np.max(np.abs(rho_num - rho_exact) / rho_exact)),
    )
    max_omega_error = float(np.max(np.abs(omega_num - omega_exact)))
    max_fd_slope_error = float(
        np.max(np.abs(finite_difference[interior] - analytic_contraction_slope[interior]))
    )
    observed_delta = float(omega_num[-1] - omega_num[0])
    if math.isclose(w, 1.0, rel_tol=0.0, abs_tol=1.0e-14):
        trend = "preserved"
    elif w < 1.0:
        trend = "increases_during_contraction"
    else:
        trend = "decreases_during_contraction"

    return {
        "w": float(w),
        "omega_sigma_initial": float(omega_sigma_initial),
        "omega_sigma_final": float(omega_num[-1]),
        "n_contraction_end": float(n_contraction_end),
        "samples": int(samples),
        "max_component_relative_error": max_component_relative_error,
        "max_omega_abs_error": max_omega_error,
        "max_interior_finite_difference_slope_abs_error": max_fd_slope_error,
        "omega_sigma_delta": observed_delta,
        "expected_contraction_trend": trend,
        "solver_success": bool(numerical.success),
    }


def run_benchmark_suite() -> dict[str, Any]:
    """Run a deterministic, bounded w/initial-fraction matrix for the submodel."""
    cases = [
        run_unsourced_case(w=w, omega_sigma_initial=omega0)
        for w in (-1.0 / 3.0, 1.0, 2.0)
        for omega0 in (0.1, 0.5, 0.9)
    ]
    limits = {
        "max_component_relative_error": 2.0e-9,
        "max_omega_abs_error": 2.0e-9,
        "max_interior_finite_difference_slope_abs_error": 1.0e-4,
    }
    passed = all(
        case["solver_success"]
        and case["max_component_relative_error"] <= limits["max_component_relative_error"]
        and case["max_omega_abs_error"] <= limits["max_omega_abs_error"]
        and case["max_interior_finite_difference_slope_abs_error"]
        <= limits["max_interior_finite_difference_slope_abs_error"]
        and (
            (case["w"] < 1.0 and case["omega_sigma_delta"] > 0.0)
            or (case["w"] == 1.0 and abs(case["omega_sigma_delta"]) <= 2.0e-9)
            or (case["w"] > 1.0 and case["omega_sigma_delta"] < 0.0)
        )
        for case in cases
    )
    return {
        "title": "A6 unsourced Bianchi-I preserve/achieve numerical benchmark",
        "model_scope": (
            "Standalone homogeneous Bianchi-I shear plus separately conserved isotropic fluid; "
            "no production-engine calls, no anisotropic stress, no energy exchange, and no H=0 samples."
        ),
        "source_equations": {
            "d_sigma2_dNc": "6 * sigma2",
            "d_rho_dNc": "3 * (1 + w) * rho",
            "Omega_sigma": "sigma2 / (sigma2 + rho) = sigma2 / (3 * H**2)",
            "dOmega_dln_a": "3 * (w - 1) * Omega_sigma * (1 - Omega_sigma)",
            "dOmega_dNc": "3 * (1 - w) * Omega_sigma * (1 - Omega_sigma)",
        },
        "scientific_admission": "NOT ADMITTED",
        "limits": limits,
        "all_checks_passed": bool(passed),
        "case_count": len(cases),
        "cases": cases,
        "interpretation_limit": (
            "A passing benchmark checks the conditional unsourced identity and its branch sign only. "
            "It does not validate a production trajectory, sourced dynamics, turnaround behavior, "
            "Bianchi IX curvature walls, global quiescence, or cyclic viability."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="optional path for the JSON report")
    args = parser.parse_args()
    report = run_benchmark_suite()
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
        print(f"WROTE {args.output}")
    else:
        print(rendered, end="")
    print(
        "A6_BENCHMARK_PASS" if report["all_checks_passed"] else "A6_BENCHMARK_FAIL",
        f"cases={report['case_count']}",
        flush=True,
    )
    return 0 if report["all_checks_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
