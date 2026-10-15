#!/usr/bin/env python3
"""No-evolution origin-anchored half-cell counterfactual and centre-Ricci test.

Initializes source-pinned production profiles, then compares copied-state
source/metric recurrences with the exact regular-centre source value s(0)=1
either inferred or imposed. Never edits production state or evolves fields.
"""
from __future__ import annotations

import copy
import json
import math
from pathlib import Path
import sys
from typing import Any

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from tools.audit_current_main_centre_proper_radial_ricci import (
    R_MAX,
    RESOLUTIONS,
    SOURCE_BLOBS,
    SOURCE_COMMIT,
    VENDOR_COMMIT,
    _arr,
    _finite,
    _max_abs,
    _state_digest,
    _verify_source_pins,
    proper_radial_ricci,
)

OUT = Path("results/origin-anchored-half-cell-ricci-response-2026-10-15")
FIT_RADIUS = 2.5
LOCAL_WINDOW_CELLS = 6
ORIGIN_SOURCE_VALUE = 1.0
MAX_ITERATIONS = 100
FIXED_POINT_TOLERANCE = 2.0e-14


def fit_origin_moment(
    r: np.ndarray,
    source: np.ndarray,
    dr: float,
    *,
    fit_radius: float,
    degree: int = 3,
    anchor_origin: bool,
    anchor_strategy: str = "constrained",
) -> dict[str, Any]:
    """Fit s(r) in z=(r/fit_radius)^2 and integrate over [0, dr/2].

    Requires at least degree+2 samples to leave residual degrees of freedom.
    Anchored strategies both impose s(0)=1. "constrained" refits the higher
    coefficients under that constraint. "intercept_only" fits without the
    constraint, then replaces only the constant coefficient; it is a pure
    coefficient-sensitivity control, not the constrained least-squares optimum.
    """
    r = np.asarray(r, dtype=float)
    source = np.asarray(source, dtype=float)
    if r.ndim != 1 or r.shape != source.shape or r.size == 0:
        raise ValueError("r and source must be matching non-empty vectors")
    if not np.all(np.isfinite(r)) or not np.all(np.isfinite(source)):
        raise ValueError("r and source must be finite")
    if not math.isfinite(dr) or dr <= 0.0:
        raise ValueError("dr must be finite and positive")
    if not math.isfinite(fit_radius) or fit_radius <= 0.0:
        raise ValueError("fit radius must be finite and positive")
    if not isinstance(degree, int) or degree < 1:
        raise ValueError("degree must be a positive integer")
    if anchor_strategy not in ("constrained", "intercept_only"):
        raise ValueError("anchor_strategy must be 'constrained' or 'intercept_only'")

    mask = r <= fit_radius + 1.0e-12 * max(1.0, fit_radius)
    x, y = r[mask], source[mask]
    nfit = int(x.size)
    if nfit < degree + 2:
        raise ValueError(
            f"degree-{degree} fit at r<={fit_radius} has {nfit} samples; "
            f"at least {degree + 2} required"
        )
    z = (x / fit_radius) ** 2
    design = np.column_stack([z**k for k in range(degree + 1)])
    coeff, _, rank, _ = np.linalg.lstsq(design, y, rcond=None)
    if rank != degree + 1:
        raise ValueError("unanchored source-fit design is rank deficient")
    if anchor_origin and anchor_strategy == "constrained":
        higher = design[:, 1:]
        higher_coeff, _, higher_rank, _ = np.linalg.lstsq(
            higher, y - ORIGIN_SOURCE_VALUE, rcond=None
        )
        if higher_rank != degree:
            raise ValueError("anchored higher-coefficient design is rank deficient")
        coeff = np.r_[ORIGIN_SOURCE_VALUE, higher_coeff]
    elif anchor_origin and anchor_strategy == "intercept_only":
        coeff = coeff.copy()
        coeff[0] = ORIGIN_SOURCE_VALUE

    fitted = design @ coeff
    residual = y - fitted
    r0 = 0.5 * dr
    moment = float(sum(
        float(coeff[k]) * r0 * (r0 / fit_radius) ** (2 * k) / (2 * k + 1)
        for k in range(degree + 1)
    ))
    point_moment = float(r0 * source[0])
    # For a saturated leave-one-out system the omitted fit has zero residual
    # degrees of freedom. Do not claim a LOO error bound in that case.
    loo_max = None
    # The standard PRESS identity applies only to an ordinary least-squares
    # projection. For a constrained refit, use the higher-power design; for
    # the intercept-only replacement, do not imply a PRESS bound at all.
    loo_design = None
    if not anchor_origin:
        loo_design = design
    elif anchor_strategy == "constrained":
        loo_design = design[:, 1:]
    if loo_design is not None and nfit >= loo_design.shape[1] + 2:
        leverages = np.diag(loo_design @ np.linalg.pinv(loo_design))
        denominator = 1.0 - leverages
        if np.all(denominator > 1.0e-10):
            loo_max = float(np.max(np.abs(residual / denominator)))

    return {
        "fit_radius": float(fit_radius),
        "degree_in_r_squared": int(degree),
        "origin_anchored": bool(anchor_origin),
        "anchor_strategy": anchor_strategy if anchor_origin else "unanchored",
        "known_s0": ORIGIN_SOURCE_VALUE if anchor_origin else None,
        "sample_count": nfit,
        "overdetermined": bool(nfit >= degree + 2),
        "coefficients_scaled_z_basis": [float(v) for v in coeff],
        "fitted_intercept": float(coeff[0]),
        "intercept_bias_vs_known_s0": float(coeff[0] - ORIGIN_SOURCE_VALUE),
        "half_cell_moment_estimate": moment,
        "production_point_sampled_moment": point_moment,
        "moment_correction_vs_point_sample": float(moment - point_moment),
        "first_cell_B_shift_vs_point_sample": float((moment - point_moment) / r0),
        "max_abs_fit_residual": _max_abs(residual),
        "rms_fit_residual": float(np.sqrt(np.mean(residual * residual))),
        "max_abs_loo_source_prediction_residual": loo_max,
        "scaled_design_condition_number": float(np.linalg.cond(design)),
    }


def _make_source_function(grid, state):
    """Reconstruct the pinned initializer's source map for a copied B profile."""
    from engine.scalar_system import cosmos_potential
    from engine.v55_matter import RHO_B_PRESENT, RHO_DM_PRESENT, RHO_R_PRESENT

    r = _arr(grid.centers)
    fields = state.scalars
    S, PS, D, PD = (_arr(fields.S), _arr(fields.PS), _arr(fields.D), _arr(fields.PD))
    phi, Pi = _arr(fields.phi), _arr(fields.Pi)
    d1 = lambda values, parity=1: grid.cell_derivative_fourth(values, parity=parity)
    Sp, Dp = d1(S, 1), d1(D, 1)
    current = -(PD * Dp)
    phi0, pi0 = 0.179055, 0.00341913
    rho0 = 0.5 * pi0**2 + float(cosmos_potential(np.array([phi0]))[0])
    rho0 += float(RHO_DM_PRESENT) + float(RHO_B_PRESENT) + float(RHO_R_PRESENT)
    H0 = math.sqrt(max(rho0 / 3.0, 0.0))
    Kt = np.full_like(r, -H0)
    Ktp = d1(Kt, 1)
    rho_matter = (
        float(RHO_DM_PRESENT) + float(RHO_B_PRESENT) + float(RHO_R_PRESENT)
    ) / (8.0 * math.pi)
    rho_phi = (0.5 * Pi * Pi + cosmos_potential(phi)) / (8.0 * math.pi)

    def source_for_B(B: np.ndarray) -> np.ndarray:
        B = np.asarray(B, dtype=float)
        rhoS = 0.5 * (PS * PS + B * Sp * Sp) + 0.5 * S * S
        rhoD = 0.5 * (PD * PD + B * Dp * Dp) - 0.5 * D * D
        rho = rhoS + rhoD + rho_phi + rho_matter
        Kr = Kt + r * Ktp + 4.0 * math.pi * r * current
        src = 1.0 - 8.0 * math.pi * r * r * rho + 2.0 * r * r * Kr * Kt + r * r * Kt * Kt
        if not _finite(src):
            raise FloatingPointError("non-finite source in copied initializer")
        return src

    return r, source_for_B


def _integrate_source(source: np.ndarray, r: np.ndarray, dr: float, first_moment: float) -> np.ndarray:
    """Integrate the custom first half-cell, retaining all later trapezoids."""
    return np.cumsum(np.r_[
        float(first_moment),
        0.5 * dr * (source[1:] + source[:-1]),
    ]) / r


def _metric_and_ricci(grid, state, vacuum, B: np.ndarray, source: np.ndarray, tag: str):
    r = _arr(grid.centers)
    if np.any(~np.isfinite(B)) or np.any(B <= 0.0):
        raise FloatingPointError(f"{tag}: invalid B profile")
    A = 1.0 / np.sqrt(B)
    a, b, X = A ** (4.0 / 3.0), A ** (-2.0 / 3.0), A ** (-1.0 / 3.0)
    d1 = lambda values, parity=1: grid.cell_derivative_fourth(values, parity=parity)
    Lambda = d1(a, 1) / (2.0 * a * a) - d1(b, 1) / (a * b) + 2.0 / r * (1.0 / b - 1.0 / a)

    geom_cf = copy.deepcopy(state.geometry)
    geom_cf.a, geom_cf.b, geom_cf.X = a.copy(), b.copy(), X.copy()
    geom_cf.Lambda = Lambda.copy()
    geom_cf.B = 0.75 * Lambda
    vendor_terms = vacuum.geometry_terms(grid, geom_cf)
    vendor_tf = _arr(vendor_terms["Rr_r"] - vendor_terms["R"] / 3.0)
    areal = r * np.sqrt(b) / X
    proper_factor = np.sqrt(a) / X
    direct = proper_radial_ricci(grid, areal, proper_factor)
    direct_tf = _arr(direct["tf_ricci"])
    gap = vendor_tf - direct_tf

    arrays = {
        f"{tag}_B": _arr(B), f"{tag}_source": _arr(source),
        f"{tag}_a": _arr(a), f"{tag}_b": _arr(b), f"{tag}_X": _arr(X),
        f"{tag}_Lambda": _arr(Lambda), f"{tag}_areal_radius": _arr(areal),
        f"{tag}_proper_factor": _arr(proper_factor),
        f"{tag}_vendor_tf_ricci": vendor_tf, f"{tag}_direct_tf_ricci": direct_tf,
        f"{tag}_signed_ricci_gap": gap, f"{tag}_signed_ricci_gap_over_r2": gap / (r * r),
    }
    invariants = {
        "a_times_b_squared_minus_one_max_abs": _max_abs(a * b * b - 1.0),
        "b_minus_X_squared_max_abs": _max_abs(b - X * X),
        "a_minus_X_minus_four_max_abs": _max_abs(a - X ** (-4.0)),
        "X6_minus_B_max_abs": _max_abs(X ** 6 - B),
        "areal_radius_minus_r_max_abs": _max_abs(areal - r),
        "all_profiles_finite": all(_finite(v) for v in arrays.values()),
    }
    return arrays, invariants, gap


def _solve_fit_mode(
    grid, state, vacuum, *, tag: str, fit_radius: float, anchored: bool,
    anchor_strategy: str = "constrained",
):
    """Run a copied source/metric recurrence, returning step 6 and fixed point."""
    before = _state_digest(state)
    r, source_for_B = _make_source_function(grid, state)
    dr = float(grid.dr)
    B_production = _arr(state.geometry.X) ** 6

    # The unmodified six-step recurrence must reproduce the pinned production X^6.
    B_baseline = np.ones_like(r)
    for _ in range(6):
        src_base = source_for_B(B_baseline)
        B_baseline = _integrate_source(src_base, r, dr, 0.5 * dr * src_base[0])
    baseline_reconstruction_error = _max_abs(B_baseline - B_production)
    if baseline_reconstruction_error > 2.0e-12 * max(1.0, _max_abs(B_production)):
        raise RuntimeError(
            f"copied point-sampled recurrence mismatch: {baseline_reconstruction_error:.6e}"
        )

    def apply_map(B: np.ndarray):
        src = source_for_B(B)
        fit = fit_origin_moment(
            r, src, dr, fit_radius=fit_radius, degree=3, anchor_origin=anchored,
            anchor_strategy=anchor_strategy,
        )
        Bnext = _integrate_source(src, r, dr, fit["half_cell_moment_estimate"])
        return Bnext, src, fit

    B = np.ones_like(r)
    history = []
    B_six = None
    fit_six = None
    source_six = None
    converged = False
    last_update = float("inf")
    for iteration in range(1, MAX_ITERATIONS + 1):
        Bnext, src, fit = apply_map(B)
        last_update = _max_abs(Bnext - B)
        history.append({
            "iteration": iteration,
            "max_abs_B_update": float(last_update),
            "half_cell_moment_estimate": fit["half_cell_moment_estimate"],
            "moment_correction_vs_point_sample": fit["moment_correction_vs_point_sample"],
            "intercept_bias_vs_known_s0": fit["intercept_bias_vs_known_s0"],
            "max_abs_fit_residual": fit["max_abs_fit_residual"],
        })
        B = Bnext
        if iteration == 6:
            B_six, fit_six, source_six = B.copy(), fit, src.copy()
        if iteration >= 6 and last_update <= FIXED_POINT_TOLERANCE * max(1.0, _max_abs(B)):
            converged = True
            break
    if B_six is None:
        raise RuntimeError("copied recurrence missed the six-step checkpoint")

    # Recompute source and map at both saved states to record fixed-point defects.
    B_six_map, src_six_fresh, fit_six_fresh = apply_map(B_six)
    B_fixed_map, source_fixed, fit_fixed = apply_map(B)
    six_defect = _max_abs(B_six_map - B_six)
    fixed_defect = _max_abs(B_fixed_map - B)
    six_arrays, six_invariants, six_gap = _metric_and_ricci(
        grid, state, vacuum, B_six, src_six_fresh, tag + "_six_step"
    )
    fixed_arrays, fixed_invariants, fixed_gap = _metric_and_ricci(
        grid, state, vacuum, B, source_fixed, tag + "_fixed_point"
    )
    r0 = float(r[0])
    report = {
        "tag": tag,
        "fit_policy": {
            "degree_in_r_squared": 3,
            "fit_radius": float(fit_radius),
            "sample_count": int(np.count_nonzero(r <= fit_radius + 1.0e-12 * max(1.0, fit_radius))),
            "anchor_origin": bool(anchored),
            "anchor_strategy": anchor_strategy if anchored else "unanchored",
            "known_s0": ORIGIN_SOURCE_VALUE if anchored else None,
        },
        "copied_recurrence": {
            "maximum_iterations": MAX_ITERATIONS,
            "tolerance": FIXED_POINT_TOLERANCE,
            "converged": converged,
            "iterations_used": len(history),
            "last_max_abs_B_update": float(last_update),
            "six_step_fixed_point_defect_max_abs": float(six_defect),
            "converged_fixed_point_defect_max_abs": float(fixed_defect),
            "history_first_six": history[:6],
            "history_tail": history[-min(5, len(history)):],
        },
        "source_moments": {
            "six_step_fit": fit_six_fresh,
            "fixed_point_fit": fit_fixed,
            "six_step_first_cell_B": float(B_six[0]),
            "fixed_point_first_cell_B": float(B[0]),
            "six_step_reconstruction_error_vs_production_X6": float(baseline_reconstruction_error),
            "input_state_sha256_before": before,
            "input_state_sha256_after": _state_digest(state),
            "input_state_unchanged": before == _state_digest(state),
        },
        "six_step_geometry_invariants": six_invariants,
        "fixed_point_geometry_invariants": fixed_invariants,
        "ricci": {
            "production_first_two_signed_gaps": [
                float(v) for v in (
                    _baseline_gap(grid, state, vacuum)[:2]
                )
            ],
            "six_step_first_two_signed_gaps": [float(v) for v in six_gap[:2]],
            "fixed_point_first_two_signed_gaps": [float(v) for v in fixed_gap[:2]],
            "six_step_first_two_signed_gap_over_r2": [float(v) for v in (six_gap[:2] / (r[:2] * r[:2]))],
            "fixed_point_first_two_signed_gap_over_r2": [float(v) for v in (fixed_gap[:2] / (r[:2] * r[:2]))],
            "production_first_six_max_abs_gap": float(_max_abs(_baseline_gap(grid, state, vacuum)[:6])),
            "six_step_first_six_max_abs_gap": float(_max_abs(six_gap[:6])),
            "fixed_point_first_six_max_abs_gap": float(_max_abs(fixed_gap[:6])),
            "fixed_point_gap_reduction_ratio_first_six": float(
                _max_abs(fixed_gap[:6]) / max(_max_abs(_baseline_gap(grid, state, vacuum)[:6]), np.finfo(float).tiny)
            ),
            "fixed_annulus_r5_to10_production_max_abs_gap": float(_max_abs(_baseline_gap(grid, state, vacuum)[(r >= 5.0) & (r <= 10.0)])),
            "fixed_annulus_r5_to10_counterfactual_max_abs_gap": float(_max_abs(fixed_gap[(r >= 5.0) & (r <= 10.0)])),
        },
    }
    arrays = {**six_arrays, **fixed_arrays}
    if not report["copied_recurrence"]["converged"]:
        raise RuntimeError(f"{tag}: corrected recurrence did not converge")
    if not report["source_moments"]["input_state_unchanged"]:
        raise RuntimeError(f"{tag}: input state was mutated")
    if not (six_invariants["all_profiles_finite"] and fixed_invariants["all_profiles_finite"]):
        raise FloatingPointError(f"{tag}: non-finite copied geometry")
    return arrays, report


def _baseline_gap(grid, state, vacuum) -> np.ndarray:
    r = _arr(grid.centers)
    g = state.geometry
    vendor = vacuum.geometry_terms(grid, g)
    vendor_tf = _arr(vendor["Rr_r"] - vendor["R"] / 3.0)
    direct = proper_radial_ricci(
        grid,
        r * np.sqrt(_arr(g.b)) / _arr(g.X),
        np.sqrt(_arr(g.a)) / _arr(g.X),
    )
    return vendor_tf - _arr(direct["tf_ricci"])


def run_audit(output_dir: Path = OUT) -> dict[str, Any]:
    """Run source-pinned baseline, free-fit, and origin-anchored counterfactuals."""
    pins = _verify_source_pins()
    from engine import v55_pirk_adapter as adapter
    from engine.production_kernel import V55ProductionKernel

    _, vacuum, _ = adapter.vendor_modules()
    reports = []
    profiles: dict[str, np.ndarray] = {}
    problems = []
    for resolution in RESOLUTIONS:
        kernel = V55ProductionKernel()
        state = kernel.initialize(
            resolution=resolution,
            r_max=R_MAX,
            amplitude=0.01,
            width=7.0,
            D_amplitude=1.0e-10,
            include_radiation=True,
        )
        baseline_state_hash = _state_digest(state)
        grid = state.grid
        r = _arr(grid.centers)
        baseline_gap = _baseline_gap(grid, state, vacuum)
        source_r, source_for_B = _make_source_function(grid, state)
        B_production = _arr(state.geometry.X) ** 6
        source_baseline = source_for_B(B_production)
        mode_specs = (
            ("fixed_unanchored", FIT_RADIUS, False, "constrained"),
            ("fixed_origin_anchored", FIT_RADIUS, True, "constrained"),
            ("fixed_intercept_only", FIT_RADIUS, True, "intercept_only"),
            ("local_six_origin_anchored", LOCAL_WINDOW_CELLS * float(grid.dr), True, "constrained"),
        )
        per_mode = {}
        mode_arrays = {}
        for tag, fit_radius, anchor, strategy in mode_specs:
            arrays, summary = _solve_fit_mode(
                grid, state, vacuum, tag=tag, fit_radius=fit_radius, anchored=anchor,
                anchor_strategy=strategy,
            )
            if not summary["source_moments"]["input_state_unchanged"]:
                problems.append(f"N={resolution}/{tag}: initialized state changed")
            if not summary["ricci"]["fixed_point_first_six_max_abs_gap"] >= 0.0:
                problems.append(f"N={resolution}/{tag}: invalid Ricci summary")
            per_mode[tag] = summary
            mode_arrays[tag] = arrays
            profiles.update({f"N{resolution}_{name}": _arr(value) for name, value in arrays.items()})

        free_fixed = per_mode["fixed_unanchored"]
        anchored_fixed = per_mode["fixed_origin_anchored"]
        intercept_only = per_mode["fixed_intercept_only"]
        anchored_local = per_mode["local_six_origin_anchored"]
        def get_delta(mode):
            return float(mode["source_moments"]["fixed_point_fit"]["moment_correction_vs_point_sample"])
        def get_gap(mode):
            return mode["ricci"]["fixed_point_first_two_signed_gaps"]
        delta_free, delta_anchor, delta_intercept, delta_local = (
            get_delta(free_fixed), get_delta(anchored_fixed), get_delta(intercept_only), get_delta(anchored_local)
        )
        gap_free, gap_anchor, gap_intercept, gap_local = (
            get_gap(free_fixed), get_gap(anchored_fixed), get_gap(intercept_only), get_gap(anchored_local)
        )

        # Same-source fits isolate estimator changes from the coupled source/metric feedback.
        same_source_contrasts = {}
        source_comparison_profiles = {
            "production_baseline_source": source_baseline,
            "unanchored_mode_fixed_point_source": mode_arrays["fixed_unanchored"]["fixed_unanchored_fixed_point_source"],
            "constrained_anchor_mode_fixed_point_source": mode_arrays["fixed_origin_anchored"]["fixed_origin_anchored_fixed_point_source"],
            "intercept_only_mode_fixed_point_source": mode_arrays["fixed_intercept_only"]["fixed_intercept_only_fixed_point_source"],
        }
        for source_name, source_samples in source_comparison_profiles.items():
            free_same = fit_origin_moment(
                r, source_samples, float(grid.dr), fit_radius=FIT_RADIUS,
                degree=3, anchor_origin=False,
            )
            constrained_same = fit_origin_moment(
                r, source_samples, float(grid.dr), fit_radius=FIT_RADIUS,
                degree=3, anchor_origin=True, anchor_strategy="constrained",
            )
            intercept_same = fit_origin_moment(
                r, source_samples, float(grid.dr), fit_radius=FIT_RADIUS,
                degree=3, anchor_origin=True, anchor_strategy="intercept_only",
            )
            same_source_contrasts[source_name] = {
                "unanchored_fit": free_same,
                "constrained_origin_fit": constrained_same,
                "intercept_only_replacement": intercept_same,
                "constrained_minus_unanchored_moment": float(
                    constrained_same["half_cell_moment_estimate"] - free_same["half_cell_moment_estimate"]
                ),
                "intercept_only_minus_unanchored_moment": float(
                    intercept_same["half_cell_moment_estimate"] - free_same["half_cell_moment_estimate"]
                ),
                "constrained_minus_intercept_only_moment": float(
                    constrained_same["half_cell_moment_estimate"] - intercept_same["half_cell_moment_estimate"]
                ),
            }
        # Frozen-source counterfactuals use one identical production source array
        # so the Ricci response isolates the direct moment/metric change before
        # any copied source/metric feedback can alter the profile.
        frozen_source_ricci_contrasts = {}
        annulus = (r >= 5.0) & (r <= 10.0)
        for tag, fit_radius, anchor, strategy in mode_specs:
            fit_same_source = fit_origin_moment(
                r, source_baseline, float(grid.dr), fit_radius=fit_radius,
                degree=3, anchor_origin=anchor, anchor_strategy=strategy,
            )
            delta_i = float(fit_same_source["moment_correction_vs_point_sample"])
            B_frozen = B_production + delta_i / r
            frozen_tag = "frozen_source_" + tag
            frozen_arrays, frozen_invariants, frozen_gap = _metric_and_ricci(
                grid, state, vacuum, B_frozen, source_baseline, frozen_tag
            )
            if not frozen_invariants["all_profiles_finite"]:
                problems.append(f"N={resolution}/{frozen_tag}: non-finite frozen-source metric")
            profiles.update({f"N{resolution}_{name}": _arr(value) for name, value in frozen_arrays.items()})
            frozen_source_ricci_contrasts[tag] = {
                "source_samples": "identical production-baseline source for all four estimator arms",
                "fit": fit_same_source,
                "delta_integral_applied": delta_i,
                "first_cell_B_production": float(B_production[0]),
                "first_cell_B_frozen_source": float(B_frozen[0]),
                "first_two_signed_ricci_gaps": [float(v) for v in frozen_gap[:2]],
                "first_two_signed_gap_change_vs_production": [
                    float(frozen_gap[i] - baseline_gap[i]) for i in range(2)
                ],
                "first_six_max_abs_ricci_gap": float(_max_abs(frozen_gap[:6])),
                "fixed_annulus_r5_to10_max_abs_ricci_gap": float(_max_abs(frozen_gap[annulus])),
                "geometry_invariants": frozen_invariants,
            }


        # Split the change from the free fixed-window fit to constrained
        # anchoring into (a) c0 replacement only and (b) all higher-coefficient
        # refit contributions. Apply each increment separately on top of the
        # *same free-fit frozen metric* so the curvature response is isolated.
        free_same = fit_origin_moment(
            r, source_baseline, float(grid.dr), fit_radius=FIT_RADIUS,
            degree=3, anchor_origin=False,
        )
        constrained_same = fit_origin_moment(
            r, source_baseline, float(grid.dr), fit_radius=FIT_RADIUS,
            degree=3, anchor_origin=True, anchor_strategy="constrained",
        )
        coeff_free = np.asarray(free_same["coefficients_scaled_z_basis"], dtype=float)
        coeff_anchor = np.asarray(constrained_same["coefficients_scaled_z_basis"], dtype=float)
        r0 = float(r[0])
        delta_c0 = float((coeff_anchor[0] - coeff_free[0]) * r0)
        delta_higher_terms = [
            float((coeff_anchor[k] - coeff_free[k]) * r0
                  * (r0 / FIT_RADIUS) ** (2 * k) / (2 * k + 1))
            for k in range(1, 4)
        ]
        delta_higher = float(sum(delta_higher_terms))
        delta_total = float(constrained_same["moment_correction_vs_point_sample"]
                            - free_same["moment_correction_vs_point_sample"])
        decomposition_roundoff = float(delta_c0 + delta_higher - delta_total)

        delta_free = float(free_same["moment_correction_vs_point_sample"])
        B_free_frozen = B_production + delta_free / r
        intervention_B = {
            "free_fit_reference": B_free_frozen,
            "c0_replacement_only": B_free_frozen + delta_c0 / r,
            "higher_coefficients_only": B_free_frozen + delta_higher / r,
            "c0_plus_higher_refit": B_free_frozen + (delta_c0 + delta_higher) / r,
        }
        coefficient_response = {}
        coefficient_gaps = {}
        for tag, Bcandidate in intervention_B.items():
            arrs, inv, gap = _metric_and_ricci(
                grid, state, vacuum, Bcandidate, source_baseline,
                "coefficient_response_" + tag,
            )
            if not inv["all_profiles_finite"]:
                problems.append(f"N={resolution}/coefficient_response_{tag}: non-finite geometry")
            profiles.update({f"N{resolution}_{name}": _arr(value) for name, value in arrs.items()})
            coefficient_gaps[tag] = gap
            coefficient_response[tag] = {
                "first_two_signed_ricci_gaps": [float(v) for v in gap[:2]],
                "first_two_gap_change_vs_free_fit": [
                    float(gap[i] - coefficient_gaps["free_fit_reference"][i])
                    for i in range(2)
                ] if "free_fit_reference" in coefficient_gaps else [0.0, 0.0],
                "first_six_max_abs_gap": float(_max_abs(gap[:6])),
                "fixed_annulus_r5_to10_max_abs_gap": float(_max_abs(gap[annulus])),
                "geometry_invariants": inv,
            }
        gap_free_c = coefficient_gaps["free_fit_reference"]
        gap_c0 = coefficient_gaps["c0_replacement_only"]
        gap_hi = coefficient_gaps["higher_coefficients_only"]
        gap_both = coefficient_gaps["c0_plus_higher_refit"]
        interaction = gap_both - gap_c0 - gap_hi + gap_free_c
        coefficient_response["decomposition_checks"] = {
            "unanchored_coefficients_scaled_z_basis": [float(v) for v in coeff_free],
            "constrained_coefficients_scaled_z_basis": [float(v) for v in coeff_anchor],
            "origin_c0": ORIGIN_SOURCE_VALUE,
            "delta_integral_from_c0_replacement_only": delta_c0,
            "delta_integral_by_higher_coefficient": delta_higher_terms,
            "delta_integral_total_higher_refit": delta_higher,
            "delta_integral_constrained_minus_free": delta_total,
            "delta_integral_decomposition_roundoff_residual": decomposition_roundoff,
            "c0_contribution_percent_of_net_moment_change": (
                100.0 * delta_c0 / delta_total if delta_total != 0.0 else None
            ),
            "higher_contribution_percent_of_net_moment_change": (
                100.0 * delta_higher / delta_total if delta_total != 0.0 else None
            ),
            "nonlinear_ricci_interaction_residual_max_abs_first_two": float(_max_abs(interaction[:2])),
            "nonlinear_ricci_interaction_residual_max_abs_first_six": float(_max_abs(interaction[:6])),
            "interpretation": "Moment decomposition is linear in polynomial coefficients; Ricci response is nonlinear in the perturbed metric, so one-at-a-time gap changes need not sum to the combined response.",
        }
        reports.append({
            "resolution": resolution,
            "r_max": R_MAX,
            "dr": float(grid.dr),
            "evolution_steps": 0,
            "input_state_sha256_before": baseline_state_hash,
            "input_state_sha256_after": _state_digest(state),
            "input_state_unchanged": baseline_state_hash == _state_digest(state),
            "production_baseline": {
                "first_two_signed_ricci_gaps": [float(v) for v in baseline_gap[:2]],
                "first_six_max_abs_ricci_gap": float(_max_abs(baseline_gap[:6])),
            },
            "counterfactual_modes": per_mode,
            "same_source_estimator_contrasts": same_source_contrasts,
            "frozen_source_ricci_contrasts": frozen_source_ricci_contrasts,
            "coefficient_response_decomposition": coefficient_response,
            "anchoring_effect_same_fixed_window": {
                "free_fit_fixed_point_moment_correction": delta_free,
                "anchored_fit_fixed_point_moment_correction": delta_anchor,
                "anchored_minus_free_moment_correction": float(delta_anchor - delta_free),
                "intercept_only_fixed_point_moment_correction": delta_intercept,
                "intercept_only_minus_free_moment_correction": float(delta_intercept - delta_free),
                "constrained_anchor_minus_intercept_only_moment_correction": float(delta_anchor - delta_intercept),
                "moment_change_relative_to_free_correction_percent": float(
                    100.0 * (delta_anchor - delta_free) / delta_free
                ),
                "free_fit_fixed_point_first_two_signed_gaps": gap_free,
                "anchored_fit_fixed_point_first_two_signed_gaps": gap_anchor,
                "anchored_minus_free_first_two_signed_gap": [
                    float(a - b) for a, b in zip(gap_anchor, gap_free)
                ],
                "intercept_only_fixed_point_first_two_signed_gaps": gap_intercept,
                "intercept_only_minus_free_first_two_signed_gap": [
                    float(a - b) for a, b in zip(gap_intercept, gap_free)
                ],
                "constrained_anchor_minus_intercept_only_first_two_signed_gap": [
                    float(a - b) for a, b in zip(gap_anchor, gap_intercept)
                ],
                "free_fit_first_six_max_abs_gap": free_fixed["ricci"]["fixed_point_first_six_max_abs_gap"],
                "anchored_fit_first_six_max_abs_gap": anchored_fixed["ricci"]["fixed_point_first_six_max_abs_gap"],
                "anchored_fit_to_free_fit_first_six_gap_ratio": float(
                    anchored_fixed["ricci"]["fixed_point_first_six_max_abs_gap"]
                    / max(free_fixed["ricci"]["fixed_point_first_six_max_abs_gap"], np.finfo(float).tiny)
                ),
                "intercept_only_first_six_max_abs_gap": intercept_only["ricci"]["fixed_point_first_six_max_abs_gap"],
                "intercept_only_to_free_first_six_gap_ratio": float(
                    intercept_only["ricci"]["fixed_point_first_six_max_abs_gap"]
                    / max(free_fixed["ricci"]["fixed_point_first_six_max_abs_gap"], np.finfo(float).tiny)
                ),
                "constrained_to_intercept_only_first_six_gap_ratio": float(
                    anchored_fixed["ricci"]["fixed_point_first_six_max_abs_gap"]
                    / max(intercept_only["ricci"]["fixed_point_first_six_max_abs_gap"], np.finfo(float).tiny)
                ),
                "local_anchored_moment_correction": delta_local,
                "local_anchored_minus_fixed_anchored_moment": float(delta_local - delta_anchor),
                "local_anchored_first_two_signed_gaps": gap_local,
                "local_anchored_first_six_max_abs_gap": anchored_local["ricci"]["fixed_point_first_six_max_abs_gap"],
            },
        })
        if baseline_state_hash != _state_digest(state):
            problems.append(f"N={resolution}: initialized state hash changed across all counterfactuals")

    report = {
        "schema": "origin_anchored_half_cell_ricci_response_v1",
        "status": "completed_diagnostic_only" if not problems else "observer_integrity_failure",
        "source_pins": pins,
        "resolutions": reports,
        "errors": problems,
        "method": {
            "origin_source_value": ORIGIN_SOURCE_VALUE,
            "regularity_basis": "s(r)=1 + O(r^2), since all non-constant source terms are multiplied by r^2",
            "fixed_fit_support": "r <= 2.5 including any cell centre exactly at r=2.5",
            "local_fit_support": "first six cell centres, radius 6*dr",
            "intercept_only_control": "fit unconstrained, then replace only c0 by 1 while retaining higher coefficients; diagnostic contrast distinct from constrained least-squares fit",
            "fixed_point_iteration": "source recomputed for current B on each copied iterate; only first half-cell moment changes; inter-centre trapezoid intervals remain",
            "ricci_discriminator": "pinned vendor trace-free Ricci minus independent proper-radial warped-product reconstruction using the same pinned grid derivative operator",
        },
        "interpretation_limits": [
            "All copied recurrences hold sampled scalar/matter fields fixed while reevaluating their energy/current contributions from the candidate metric.",
            "No fields are evolved, and no production state is modified.",
            "The exact source value at r=0 is known, but the exact continuous integral over the production half-cell remains unknown.",
            "The fitted-minus-point moment is an estimator correction, not a measured quadrature error.",
            "A smaller signed Ricci gap is localization evidence for initialized slices, not proof of continuum or evolved convergence.",
            "Free-fit, anchored-fit, and local-fit results are separate counterfactuals; anchoring is an explicit additional modeling constraint.",
        ],
        "production_or_vendor_files_modified": False,
        "no_evolution": True,
        "scientific_admission": False,
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    np.savez_compressed(output_dir / "profiles.npz", **profiles)
    if report["errors"]:
        raise RuntimeError("observer integrity failures: " + "; ".join(problems))
    if any(not row["input_state_unchanged"] for row in reports):
        raise RuntimeError("input state hash changed")
    return report


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
        "resolutions": [
            {
                "N": row["resolution"],
                "moment_and_ricci": row["anchoring_effect_same_fixed_window"],
            }
            for row in report["resolutions"]
        ],
        "evolution_steps": 0,
        "production_or_vendor_files_modified": False,
        "scientific_admission": False,
        "output_dir": str(args.output_dir),
    }
    print(json.dumps(condensed, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
