#!/usr/bin/env python3
"""High-accuracy continuous initializer source oracle for the origin half-cell.

This research-only observer solves the continuous radial integral equation
for the pinned initializer's analytic Gaussian scalar profiles. It compares
the continuous moment with production point-sampled and polynomial-fit
moments, and separates estimator correction from discrete-source mismatch.
It performs no production edit and takes no evolution steps.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
import sys
from typing import Any

import numpy as np
from scipy.integrate import quad, solve_ivp

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from engine import v55_pirk_adapter as adapter
from engine.production_kernel import V55ProductionKernel
from engine.scalar_system import cosmos_potential
from engine.v55_matter import RHO_B_PRESENT, RHO_DM_PRESENT, RHO_R_PRESENT
from tools.audit_current_main_centre_proper_radial_ricci import (
    R_MAX,
    RESOLUTIONS,
    SOURCE_COMMIT,
    VENDOR_COMMIT,
    _arr,
    _finite,
    _max_abs,
    _state_digest,
    _verify_source_pins,
)
from tools.audit_origin_anchored_half_cell_ricci_response import (
    FIT_RADIUS,
    LOCAL_WINDOW_CELLS,
    _make_source_function,
    fit_origin_moment,
)

OUT = Path("results/origin-half-cell-continuous-source-oracle-2026-10-15")
AMPLITUDE = 0.01
WIDTH = 7.0
D_AMPLITUDE = 1.0e-10
PHI0 = 0.179055
PI0 = 0.00341913
ODE_RTOL = 5.0e-14
ODE_ATOL = 1.0e-18
QUAD_EPSABS = 2.0e-15
QUAD_EPSREL = 2.0e-13


def _continuum_parameters(
    amplitude: float = AMPLITUDE,
    width: float = WIDTH,
    D_amplitude: float = D_AMPLITUDE,
) -> dict[str, float]:
    """Return analytic initializer constants and origin-series coefficients."""
    if not (math.isfinite(amplitude) and math.isfinite(width) and width > 0.0
            and math.isfinite(D_amplitude)):
        raise ValueError("amplitude, width and D_amplitude must be finite; width must be positive")
    rho_budget = (
        0.5 * PI0**2
        + float(cosmos_potential(np.array([PHI0]))[0])
        + float(RHO_DM_PRESENT)
        + float(RHO_B_PRESENT)
        + float(RHO_R_PRESENT)
    )
    H0 = math.sqrt(max(rho_budget / 3.0, 0.0))
    Kt = -H0
    rho_constant = rho_budget / (8.0 * math.pi)
    cancellation_residual = -8.0 * math.pi * rho_constant + 3.0 * Kt**2

    # Exact regular-centre expansion for this initializer: PD=D makes the
    # D kinetic and negative D-potential terms cancel, and
    # 3 Kt^2 = 8 pi rho_constant. Thus s(r)=1+s2*r^2+s4*r^4+s6*r^6+...
    s2 = -4.0 * math.pi * amplitude**2
    b2 = s2 / 3.0
    s4 = (
        8.0 * math.pi * amplitude**2 / width**2
        - 16.0 * math.pi * (amplitude**2 + D_amplitude**2) / width**4
        + 16.0 * math.pi * D_amplitude**2 * Kt / width**2
    )
    b4 = s4 / 5.0
    s6 = (
        -8.0 * math.pi * amplitude**2 / width**4
        + 32.0 * math.pi * (amplitude**2 + D_amplitude**2) / width**6
        - 16.0 * math.pi * b2 * (amplitude**2 + D_amplitude**2) / width**4
        - 32.0 * math.pi * D_amplitude**2 * Kt / width**4
    )
    b6 = s6 / 7.0
    return {
        "amplitude": float(amplitude),
        "width": float(width),
        "D_amplitude": float(D_amplitude),
        "H0": float(H0),
        "Kt": float(Kt),
        "rho_constant": float(rho_constant),
        "constant_cancellation_residual_float64": float(cancellation_residual),
        "s2": float(s2),
        "s4": float(s4),
        "s6": float(s6),
        "b2": float(b2),
        "b4": float(b4),
        "b6": float(b6),
    }


def continuous_source(
    radius: float | np.ndarray,
    B: float | np.ndarray,
    *,
    params: dict[str, float],
) -> float | np.ndarray:
    """Evaluate the continuous regular-centre source using analytic S,D and derivatives.

    This algebraically simplifies the production source with Kt=-H0 and
    3 H0^2=8 pi rho_constant, retaining the Gaussian scalar-gradient and
    momentum-current terms. That avoids subtracting nearly cancelling
    O(r^2) constant energy/curvature terms in floating-point arithmetic.
    """
    r = np.asarray(radius, dtype=float)
    metric_B = np.asarray(B, dtype=float)
    if np.any(~np.isfinite(r)) or np.any(r < 0.0):
        raise ValueError("radius must be finite and nonnegative")
    if np.any(~np.isfinite(metric_B)) or np.any(metric_B <= 0.0):
        raise ValueError("B must be finite and positive")
    amp, width, damp, Kt = (
        params["amplitude"], params["width"], params["D_amplitude"], params["Kt"]
    )
    S = amp * np.exp(-(r / width) ** 2)
    D = damp * np.exp(-(r / width) ** 2)
    Sp = -2.0 * r * S / width**2
    Dp = -2.0 * r * D / width**2
    j = -D * Dp
    bracket = (
        -4.0 * math.pi * S**2
        -4.0 * math.pi * metric_B * (Sp**2 + Dp**2)
        +8.0 * math.pi * r * j * Kt
    )
    values = 1.0 + r**2 * bracket
    if not np.all(np.isfinite(values)):
        raise FloatingPointError("continuous source generated non-finite values")
    return float(values) if values.ndim == 0 else values


def solve_continuous_moment(
    r0: float,
    *,
    amplitude: float = AMPLITUDE,
    width: float = WIDTH,
    D_amplitude: float = D_AMPLITUDE,
    profile_radii: np.ndarray | None = None,
) -> dict[str, Any]:
    """Solve d(r B)/dr=s(r,B), cross-check I(r0), optionally sample its profile."""
    if not math.isfinite(r0) or r0 <= 0.0:
        raise ValueError("r0 must be finite and positive")
    profile = None
    stop_radius = float(r0)
    if profile_radii is not None:
        profile = np.asarray(profile_radii, dtype=float)
        if profile.ndim != 1 or not np.all(np.isfinite(profile)) or np.any(profile < 0.0):
            raise ValueError("profile_radii must be a finite, nonnegative vector")
        if profile.size:
            stop_radius = max(stop_radius, float(np.max(profile)))
        if np.any(profile < 1.0e-12):
            raise ValueError("profile_radii must stay above the positive series-start radius")
    p = _continuum_parameters(amplitude, width, D_amplitude)
    eps = min(1.0e-3, r0 / 100.0)
    y0 = p["b2"] * eps**2 + p["b4"] * eps**4 + p["b6"] * eps**6

    def rhs(radius: float, y: np.ndarray) -> np.ndarray:
        B = 1.0 + float(y[0])
        r = float(radius)
        amp, w, damp, Kt = (
            p["amplitude"], p["width"], p["D_amplitude"], p["Kt"]
        )
        S = amp * math.exp(-(r / w) ** 2)
        D = damp * math.exp(-(r / w) ** 2)
        Sp = -2.0 * r * S / w**2
        Dp = -2.0 * r * D / w**2
        j = -D * Dp
        bracket = (
            -4.0 * math.pi * S**2
            -4.0 * math.pi * B * (Sp**2 + Dp**2)
            +8.0 * math.pi * r * j * Kt
        )
        # y=B-1 avoids subtracting two O(1) values before cancellation.
        return np.array([(-float(y[0]) + r**2 * bracket) / r])

    solution = solve_ivp(
        rhs,
        (eps, stop_radius),
        np.array([y0], dtype=float),
        method="DOP853",
        rtol=ODE_RTOL,
        atol=ODE_ATOL,
        dense_output=True,
        max_step=max((stop_radius - eps) / 128.0, 1.0e-6),
    )
    if not solution.success or solution.sol is None:
        raise RuntimeError(f"continuous-source ODE did not converge: {solution.message}")
    B_r0 = 1.0 + float(solution.sol(r0)[0])
    if not math.isfinite(B_r0) or B_r0 <= 0.0:
        raise FloatingPointError("continuous-source metric B(r0) is invalid")

    # Independent quadrature of s(r,B(r)) cross-checks the ODE identity I=r B.
    def source_on_solution(radius: float) -> float:
        B = 1.0 + float(solution.sol(radius)[0])
        return float(continuous_source(radius, B, params=p))

    quad_tail, quad_error = quad(
        source_on_solution, eps, r0,
        epsabs=QUAD_EPSABS, epsrel=QUAD_EPSREL, limit=200,
    )
    integral_initial_series = eps + p["b2"] * eps**3 + p["b4"] * eps**5 + p["b6"] * eps**7
    I_from_quadrature = float(integral_initial_series + quad_tail)
    I_from_metric = float(r0 * B_r0)
    s_cont_r0 = float(continuous_source(r0, B_r0, params=p))
    continuous_point_moment = float(r0 * s_cont_r0)
    series_moment_s6 = float(
        r0 + p["s2"] * r0**3 / 3.0 + p["s4"] * r0**5 / 5.0 + p["s6"] * r0**7 / 7.0
    )
    sampled_profile = {}
    if profile is not None:
        B_profile = 1.0 + np.asarray(solution.sol(profile)[0], dtype=float) if profile.size else np.array([], dtype=float)
        source_profile = np.asarray(continuous_source(profile, B_profile, params=p), dtype=float)
        if not _finite(B_profile) or not _finite(source_profile) or np.any(B_profile <= 0.0):
            raise FloatingPointError("continuous sampled profile is non-finite or invalid")
        sampled_profile = {
            "continuous_metric_profile_radii": [float(v) for v in profile],
            "continuous_B_profile_at_radii": [float(v) for v in B_profile],
            "continuous_source_profile_at_radii": [float(v) for v in source_profile],
        }
    return {
        "r0": float(r0),
        "continuous_ode_stop_radius": float(stop_radius),
        **sampled_profile,
        "epsilon_start": float(eps),
        "ode_success": bool(solution.success),
        "ode_message": str(solution.message),
        "ode_relative_tolerance": ODE_RTOL,
        "ode_absolute_tolerance": ODE_ATOL,
        "quad_absolute_tolerance": QUAD_EPSABS,
        "quad_relative_tolerance": QUAD_EPSREL,
        "parameters": p,
        "continuous_B_at_r0": B_r0,
        "continuous_source_at_r0": s_cont_r0,
        "continuous_point_sampled_moment": continuous_point_moment,
        "continuous_integral_from_rB": I_from_metric,
        "continuous_integral_from_independent_quad": I_from_quadrature,
        "quad_estimated_error_on_tail": float(quad_error),
        "ode_vs_quad_integral_abs_difference": abs(I_from_metric - I_from_quadrature),
        "exact_continuous_point_moment_correction": float(I_from_metric - continuous_point_moment),
        "s2_s4_s6_series_moment": series_moment_s6,
        "series_moment_minus_continuous_integral": float(series_moment_s6 - I_from_metric),
    }


def run_audit(output_dir: Path = OUT) -> dict[str, Any]:
    """Run a checksum-pinned comparison for N=160/240/320."""
    pins = _verify_source_pins()
    from engine.production_kernel import V55ProductionKernel

    _, _, _ = adapter.vendor_modules()
    reports = []
    for n in RESOLUTIONS:
        kernel = V55ProductionKernel()
        state = kernel.initialize(
            resolution=n,
            r_max=R_MAX,
            amplitude=AMPLITUDE,
            width=WIDTH,
            D_amplitude=D_AMPLITUDE,
            include_radiation=True,
        )
        state_hash = _state_digest(state)
        grid = state.grid
        r = _arr(grid.centers)
        dr = float(grid.dr)
        r0 = float(r[0])
        Bprod = _arr(state.geometry.X) ** 6
        rgrid, source_for_B = _make_source_function(grid, state)
        source_prod = _arr(source_for_B(Bprod))
        if rgrid.shape != r.shape or not np.array_equal(rgrid, r):
            raise RuntimeError(f"N={n}: production source helper grid does not match the initialized grid")
        if not _finite(source_prod) or not _finite(Bprod) or np.any(Bprod <= 0.0):
            raise FloatingPointError(f"N={n}: production initializer profiles are invalid")

        # Confirm that this continuous oracle really represents the pinned,
        # unmodified analytic scalar-profile family at the sampled coordinates.
        S_expected = AMPLITUDE * np.exp(-(r / WIDTH) ** 2)
        D_expected = D_AMPLITUDE * np.exp(-(r / WIDTH) ** 2)
        checks = {
            "S_profile_max_abs_difference": _max_abs(_arr(state.scalars.S) - S_expected),
            "D_profile_max_abs_difference": _max_abs(_arr(state.scalars.D) - D_expected),
            "PS_profile_max_abs": _max_abs(_arr(state.scalars.PS)),
            "PD_minus_D_max_abs": _max_abs(_arr(state.scalars.PD) - _arr(state.scalars.D)),
            "phi_minus_phi0_max_abs": _max_abs(_arr(state.scalars.phi) - PHI0),
            "Pi_minus_Pi0_max_abs": _max_abs(_arr(state.scalars.Pi) - PI0),
        }
        if checks["S_profile_max_abs_difference"] > 2.0e-15:
            raise RuntimeError(f"N={n}: S profile does not match continuous oracle input")
        if checks["D_profile_max_abs_difference"] > 2.0e-22:
            raise RuntimeError(f"N={n}: D profile does not match continuous oracle input")
        if checks["PS_profile_max_abs"] != 0.0 or checks["PD_minus_D_max_abs"] != 0.0:
            raise RuntimeError(f"N={n}: momentum profiles do not match continuous oracle assumptions")
        if checks["phi_minus_phi0_max_abs"] > 2.0e-15 or checks["Pi_minus_Pi0_max_abs"] > 2.0e-15:
            raise RuntimeError(f"N={n}: phi/Pi profiles do not match continuous oracle assumptions")

        cont = solve_continuous_moment(r0, profile_radii=r)
        if cont["ode_vs_quad_integral_abs_difference"] > 2.0e-12:
            raise RuntimeError(
                f"N={n}: continuous ODE and independent quadrature disagree by "
                f"{cont['ode_vs_quad_integral_abs_difference']:.6e}"
            )
        continuous_source_at_r0 = cont["continuous_source_at_r0"]
        production_point_moment = float(r0 * source_prod[0])
        continuous_point_moment = cont["continuous_point_sampled_moment"]
        point_source_difference = float(production_point_moment - continuous_point_moment)

        fit_specs = (
            ("fixed_unanchored", FIT_RADIUS, False, "constrained"),
            ("fixed_origin_anchored", FIT_RADIUS, True, "constrained"),
            ("fixed_intercept_only", FIT_RADIUS, True, "intercept_only"),
            ("local_six_origin_anchored", LOCAL_WINDOW_CELLS * dr, True, "constrained"),
        )
        fit_rows = {}
        for tag, radius, anchored, strategy in fit_specs:
            fit = fit_origin_moment(
                r, source_prod, dr, fit_radius=radius,
                degree=3, anchor_origin=anchored, anchor_strategy=strategy,
            )
            moment_estimate = float(fit["half_cell_moment_estimate"])
            correction = float(fit["moment_correction_vs_point_sample"])
            signed_error = float(moment_estimate - cont["continuous_integral_from_rB"])
            decomposition = (
                correction
                + point_source_difference
                + continuous_point_moment
                - cont["continuous_integral_from_rB"]
            )
            continuous_source_profile = np.asarray(
                cont["continuous_source_profile_at_radii"], dtype=float
            )
            continuous_fit = fit_origin_moment(
                r, continuous_source_profile, dr, fit_radius=radius,
                degree=3, anchor_origin=anchored, anchor_strategy=strategy,
            )
            continuous_fit_moment = float(continuous_fit["half_cell_moment_estimate"])
            fit_rows[tag] = {
                "fit": fit,
                "moment_estimate": moment_estimate,
                "estimator_correction_vs_production_point_sample": correction,
                "signed_moment_error_vs_continuous_oracle": signed_error,
                "production_point_minus_continuous_point_moment": point_source_difference,
                "continuous_point_minus_continuous_integral": float(
                    continuous_point_moment - cont["continuous_integral_from_rB"]
                ),
                "error_decomposition_sum": float(decomposition),
                "error_decomposition_closure_residual": float(signed_error - decomposition),
                "signed_error_relative_to_r0": float(signed_error / r0),
                "fit_on_continuous_source_samples": {
                    "fit": continuous_fit,
                    "moment_estimate": continuous_fit_moment,
                    "signed_moment_error_vs_continuous_oracle": float(
                        continuous_fit_moment - cont["continuous_integral_from_rB"]
                    ),
                    "continuous_sampled_point_moment": float(r0 * continuous_source_profile[0]),
                    "estimator_correction_vs_continuous_sampled_point": float(
                        continuous_fit["moment_correction_vs_point_sample"]
                    ),
                    "sampled_source_uses_continuum_metric_and_analytic_derivatives": True,
                },
            }
        point_estimate = production_point_moment
        point_error = point_estimate - cont["continuous_integral_from_rB"]
        reports.append({
            "resolution": int(n),
            "dr": dr,
            "r0": r0,
            "evolution_steps": 0,
            "initial_state_sha256_before": state_hash,
            "initial_state_sha256_after": _state_digest(state),
            "input_state_unchanged": state_hash == _state_digest(state),
            "profile_identity_checks": checks,
            "production_profile": {
                "B_first_cell": float(Bprod[0]),
                "source_first_cell": float(source_prod[0]),
                "point_sampled_half_cell_moment": production_point_moment,
                "point_moment_signed_error_vs_continuous_oracle": point_error,
                "source_point_difference_vs_continuum": point_source_difference,
            },
            "continuous_oracle": cont,
            "fit_estimates": fit_rows,
        })
        if state_hash != _state_digest(state):
            raise RuntimeError(f"N={n}: input state changed during read-only oracle comparison")

    output = {
        "schema": "origin_half_cell_continuous_source_oracle_v1",
        "status": "completed_diagnostic_only",
        "source_pins": pins,
        "parameters": {
            "resolution_set": list(RESOLUTIONS),
            "r_max": R_MAX,
            "amplitude": AMPLITUDE,
            "width": WIDTH,
            "D_amplitude": D_AMPLITUDE,
            "include_radiation": True,
            "evolution_steps": 0,
        },
        "continuous_model": {
            "source_formula": "s=1+r^2*(-4*pi*S^2 -4*pi*B*(Sprime^2+Dprime^2)+8*pi*r*j*Kt), j=-D*Dprime, S=A exp(-(r/w)^2), D=D0 exp(-(r/w)^2)",
            "metric_equation": "d(r*B)/dr=s(r,B), B(0)=1; equivalent to the continuous initializer's integral definition",
            "known_origin_value": 1.0,
            "constant_balance": "-8*pi*rho_constant + 3*Kt^2 = 0 analytically for initializer H0",
            "series_start": "B(epsilon)=1+b2 epsilon^2+b4 epsilon^4+b6 epsilon^6 from the analytic regular-centre series",
            "cross_check": "solve_ivp(DOP853) for the metric ODE versus independent scipy.quad integration of s(r,B(r))",
        },
        "interpretation_limits": [
            "The continuous oracle uses analytic scalar profiles and derivatives with the continuum radial integral equation; it is independent of the cell-centred source fit.",
            "The oracle is a high-accuracy numerical reference, not a formally rigorous interval enclosure. ODE/quad agreement and the regular-centre series residual are numerical checks, not a proof bound.",
            "The production point sample uses finite-difference scalar derivatives and a discrete metric. The recorded source-point difference keeps those discretization effects separate from estimator correction.",
            "The fit-minus-production-point quantity is a correction; the fit's signed error against the independent continuous oracle is separately reported.",
            "No production source, vendor code, dependency, campaign setting or initialized state is modified; no field is evolved.",
        ],
        "production_or_vendor_files_modified": False,
        "no_evolution": True,
        "scientific_admission": False,
        "errors": [],
        "resolutions": reports,
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "report.json").write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    if output["errors"]:
        raise RuntimeError("continuous-source oracle had observer errors")
    return output


def main() -> None:
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=OUT)
    args = parser.parse_args()
    report = run_audit(args.output_dir)
    condensed = {
        "schema": report["schema"],
        "status": report["status"],
        "source_pins_match": report["source_pins"]["engine_vendor_dependency_diff_empty"],
        "rows": [
            {
                "N": row["resolution"],
                "continuous_I": row["continuous_oracle"]["continuous_integral_from_rB"],
                "continuous_I_quad": row["continuous_oracle"]["continuous_integral_from_independent_quad"],
                "production_point_moment": row["production_profile"]["point_sampled_half_cell_moment"],
                "estimates": {
                    key: {
                        "moment": val["moment_estimate"],
                        "signed_error_vs_continuous_oracle": val["signed_moment_error_vs_continuous_oracle"],
                        "error_decomposition_closure_residual": val["error_decomposition_closure_residual"],
                        "continuous_source_sample_fit_error": val["fit_on_continuous_source_samples"]["signed_moment_error_vs_continuous_oracle"],
                    }
                    for key, val in row["fit_estimates"].items()
                },
            }
            for row in report["resolutions"]
        ],
        "no_evolution": True,
        "scientific_admission": False,
        "output_dir": str(args.output_dir),
    }
    print(json.dumps(condensed, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
