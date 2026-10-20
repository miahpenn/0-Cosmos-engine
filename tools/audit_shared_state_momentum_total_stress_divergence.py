#!/usr/bin/env python3
"""Shared accepted-state momentum / total-stress divergence witness; diagnostic-only."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from tools import audit_current_source_momentum_constraint_tangent as momentum
from tools import audit_total_stress_energy_divergence as stress
from tools import audit_brown_vacuum_momentum_propagation as brown

INITIAL_SHA256 = "a256b57f274ee120b55c0f805ed488d91382a9799647c8f9811d7df82aa337c0"
ACCEPTED_SHA256 = "c6acb3711d59a85c804a88208c9bc0356e8227a806c0684ae3a8afac72c9091a"
RESOLUTION, R_MAX, DT, EPSILON = 60, 80.0, 0.015, 1.0e-5


def summary(x):
    a = np.asarray(x, dtype=np.float64)
    if not np.all(np.isfinite(a)):
        raise FloatingPointError("non-finite value in shared-state calculation")
    i = int(np.argmax(np.abs(a)))
    return {"max_abs": float(abs(a[i])), "signed_at_max": float(a[i]),
            "index": i, "rms": float(np.sqrt(np.mean(a*a))),
            "cell_0": float(a[0]), "cell_1": float(a[1])}


def require_digest(state, expected: str, label: str) -> str:
    actual = momentum.state_digest(state)
    if actual != expected:
        raise RuntimeError(f"{label} digest mismatch: expected {expected}, got {actual}")
    return actual


def matter_source_momentum_tangent(state, projections):
    """Apply the exact M_geom functional derivative to explicit Aa/K matter sources."""
    grid, geom = state.grid, state.geometry
    r = np.asarray(grid.centers, dtype=np.float64)
    d1 = lambda values: np.asarray(grid.cell_derivative_fourth(values, parity=1), dtype=np.float64)
    x_r = d1(np.asarray(geom.X, dtype=np.float64))
    b_r = d1(np.asarray(geom.b, dtype=np.float64))
    coefficient_Aa = -3.0 * x_r / np.asarray(geom.X) + 3.0 / r + 1.5 * b_r / np.asarray(geom.b)

    def contribution(part):
        rho = np.asarray(part["rho"], dtype=np.float64)
        pr = np.asarray(part["pr"], dtype=np.float64)
        pt = np.asarray(part["pt"], dtype=np.float64)
        aa_source = -(16.0 * np.pi / 3.0) * np.asarray(geom.alpha) * (pr - pt)
        k_source = 4.0 * np.pi * np.asarray(geom.alpha) * (rho + pr + 2.0 * pt)
        return d1(aa_source) - (2.0 / 3.0) * d1(k_source) + coefficient_Aa * aa_source

    pieces = {}
    for name, part in projections["scalar_parts"].items():
        pieces[name] = contribution(part)
    for name, part in projections["fluid_parts"].items():
        pieces[name] = contribution(part)
    total = contribution({key: projections[key] for key in ("rho", "pr", "pt")})
    summed = sum(pieces.values(), np.zeros_like(total))
    operand_scale = sum((np.abs(value) for value in pieces.values()), np.zeros_like(total))
    decomposition_gap = summed - total
    # Bound cancellation from the primitive stress operands, not only the
    # already-cancelled tangent outputs. In the outer cells pr and pt can be
    # nearly equal, so forming (pr - pt) loses digits before D1 is applied.
    # Propagate the same 512*eps roundoff factor through the actual pinned D1
    # stencil's absolute weights; this changes the error estimate, not physics
    # tolerances or any scientific-admission gate.
    all_parts = tuple(projections["scalar_parts"].values()) + tuple(projections["fluid_parts"].values())
    sum_abs_rho = sum((np.abs(np.asarray(part["rho"], dtype=np.float64)) for part in all_parts),
                      np.zeros_like(total))
    sum_abs_pr = sum((np.abs(np.asarray(part["pr"], dtype=np.float64)) for part in all_parts),
                     np.zeros_like(total))
    sum_abs_pt = sum((np.abs(np.asarray(part["pt"], dtype=np.float64)) for part in all_parts),
                     np.zeros_like(total))
    alpha_abs = np.abs(np.asarray(geom.alpha, dtype=np.float64))
    aa_source_operand_scale = (16.0 * np.pi / 3.0) * alpha_abs * (
        sum_abs_pr + sum_abs_pt + np.abs(np.asarray(projections["pr"], dtype=np.float64))
        + np.abs(np.asarray(projections["pt"], dtype=np.float64))
    )
    k_source_operand_scale = 4.0 * np.pi * alpha_abs * (
        sum_abs_rho + sum_abs_pr + 2.0 * sum_abs_pt
        + np.abs(np.asarray(projections["rho"], dtype=np.float64))
        + np.abs(np.asarray(projections["pr"], dtype=np.float64))
        + 2.0 * np.abs(np.asarray(projections["pt"], dtype=np.float64))
    )
    # Recover the actual linear D1 stencil row weights by applying it to basis
    # vectors. The grid has only 60 cells here; this is a cheap diagnostic-only
    # construction and avoids assuming a generic boundary stencil norm.
    basis = np.eye(total.size, dtype=np.float64)
    d1_abs = np.abs(np.column_stack([d1(basis[:, j]) for j in range(total.size)]))
    propagated_source_scale = (
        d1_abs @ aa_source_operand_scale
        + (2.0 / 3.0) * (d1_abs @ k_source_operand_scale)
        + np.abs(coefficient_Aa) * aa_source_operand_scale
    )
    bound_scale = operand_scale + np.abs(total) + propagated_source_scale
    bound = 512.0 * np.finfo(np.float64).eps * bound_scale
    failed = np.abs(decomposition_gap) > bound
    if not np.all(~failed):
        ratio = np.divide(
            np.abs(decomposition_gap), bound,
            out=np.full_like(decomposition_gap, np.inf), where=bound > 0.0,
        )
        worst = int(np.argmax(ratio))
        # Diagnostic-only signed breakdown. Preserve the assertion and bound;
        # expose per-sector operator pieces at the worst cell and centre controls.
        def sector_terms(part):
            rho = np.asarray(part["rho"], dtype=np.float64)
            pr = np.asarray(part["pr"], dtype=np.float64)
            pt = np.asarray(part["pt"], dtype=np.float64)
            aa_source = -(16.0 * np.pi / 3.0) * np.asarray(geom.alpha) * (pr - pt)
            k_source = 4.0 * np.pi * np.asarray(geom.alpha) * (rho + pr + 2.0 * pt)
            aa_derivative = d1(aa_source)
            k_derivative = -(2.0 / 3.0) * d1(k_source)
            coefficient_term = coefficient_Aa * aa_source
            return {
                "Aa_source": aa_source,
                "K_source": k_source,
                "Aa_derivative_term": aa_derivative,
                "K_derivative_term": k_derivative,
                "geometric_coefficient_term": coefficient_term,
                "combined": aa_derivative + k_derivative + coefficient_term,
            }
        term_arrays = {name: sector_terms(part) for name, part in {
            **projections["scalar_parts"], **projections["fluid_parts"]
        }.items()}
        total_terms = sector_terms({key: projections[key] for key in ("rho", "pr", "pt")})
        diagnostic_cells = sorted(set([worst, 0, 1]))
        # Cell 58 is the second-to-last cell on the fixed 60-cell grid.
        # The pinned fourth-order D1 outer closure uses the final five cells
        # for i=n-2; record that exact support, not a window around 'worst'.
        target_cell = 58
        if total.size != 60 or target_cell != total.size - 2:
            raise RuntimeError("Cell-58 diagnostic requires the pinned 60-cell grid")
        d1_stencil_indices = list(range(total.size - 5, total.size))
        d1_method = getattr(grid.cell_derivative_fourth, "__func__", grid.cell_derivative_fourth)
        d1_implementation = (
            f"{getattr(d1_method, '__module__', type(grid).__module__)}."
            f"{getattr(d1_method, '__qualname__', type(grid).__qualname__)}"
        )
        source_neighborhood = []
        source_roundoff_factor = 512.0 * np.finfo(np.float64).eps
        for j in d1_stencil_indices:
            source_entry = {"cell": j}
            for source_name in ("K_source", "Aa_source"):
                sector_values = {
                    name: float(terms[source_name][j])
                    for name, terms in term_arrays.items()
                }
                sector_sum = float(sum(sector_values.values()))
                independent = float(total_terms[source_name][j])
                local_operand_scale = float(sum(abs(value) for value in sector_values.values()))
                local_bound = source_roundoff_factor * max(local_operand_scale, abs(independent))
                source_entry[source_name] = {
                    "sector_values": sector_values,
                    "sector_sum": sector_sum,
                    "independent_total": independent,
                    "signed_gap": sector_sum - independent,
                    "local_bound": local_bound,
                    "operand_scale": local_operand_scale,
                    "gap_over_bound": (
                        abs(sector_sum - independent) / local_bound if local_bound else None
                    ),
                }
            source_neighborhood.append(source_entry)
        diagnostic_cells = sorted(set([worst, 0, 1]))
        detail = []
        for idx in diagnostic_cells:
            sectors = {
                name: {key: float(arr[idx]) for key, arr in terms.items()}
                for name, terms in term_arrays.items()
            }
            total_detail = {key: float(arr[idx]) for key, arr in total_terms.items()}
            abs_intermediate_scale = float(sum(
                abs(sectors[name][term])
                for name in sectors
                for term in ("Aa_derivative_term", "K_derivative_term", "geometric_coefficient_term")
            ))
            signed_sector_sum = float(sum(sectors[name]["combined"] for name in sectors))
            detail.append({
                "cell": idx,
                "sector_terms": sectors,
                "total_terms": total_detail,
                "signed_sector_sum": signed_sector_sum,
                "independent_total": float(total[idx]),
                "gap": float(decomposition_gap[idx]),
                "bound": float(bound[idx]),
                "gap_over_bound": float(ratio[idx]),
                "abs_intermediate_operand_scale": abs_intermediate_scale,
                "gap_over_intermediate_scale": (
                    float(abs(decomposition_gap[idx]) / abs_intermediate_scale)
                    if abs_intermediate_scale else None
                ),
            })
        raise RuntimeError(
            "sector matter-source tangent decomposition failed: "
            f"failed_cells={int(np.count_nonzero(failed))}/{failed.size}, "
            f"max_abs_gap={float(np.max(np.abs(decomposition_gap))):.17g}, "
            f"max_bound={float(np.max(bound)):.17g}, "
            f"worst_ratio={float(ratio[worst]):.17g}, worst_index={worst}, "
            f"gap_at_worst={float(decomposition_gap[worst]):.17g}, "
            f"bound_at_worst={float(bound[worst]):.17g}, "
            f"summed_at_worst={float(summed[worst]):.17g}, "
            f"total_at_worst={float(total[worst]):.17g}, "
            f"operand_scale_at_worst={float(operand_scale[worst]):.17g}, "
            f"d1_stencil={json.dumps({'target_cell': target_cell, 'indices': d1_stencil_indices, 'implementation': d1_implementation, 'order': 1, 'parity': 1, 'closure': 'outer one-sided polynomial-moment five-node closure'}, sort_keys=True, separators=(',', ':'))}, "
            f"source_neighborhood={json.dumps(source_neighborhood, sort_keys=True, separators=(',', ':'))}, "
            f"diagnostic_cells={json.dumps(detail, sort_keys=True, separators=(',', ':'))}"
        )
    return {
        "total": total,
        "pieces": pieces,
        "sector_sum": summed,
        "decomposition_gap": decomposition_gap,
        "decomposition_bound": bound,
    }


def connection_correction(state):
    """Evaluate the source-derived off-constraint C^r correction on this exact state.

    Diagnostic decomposition only: do not fold this into Brown's predictor
    until the overlapping connection contribution is reconciled.
    """
    grid, geom = state.grid, state.geometry
    r = np.asarray(grid.centers, dtype=np.float64)
    # The pinned API's parity argument describes the INPUT profile:
    # +1 for even input, -1 for odd input.  The derivative output has opposite
    # parity; do not pass 0 or infer the argument from the output parity.
    d1_even_input = lambda values: np.asarray(
        grid.cell_derivative_fourth(values, parity=+1), dtype=np.float64
    )
    d1_odd_input = lambda values: np.asarray(
        grid.cell_derivative_fourth(values, parity=-1), dtype=np.float64
    )
    a, b, x = (np.asarray(getattr(geom, key), dtype=np.float64) for key in ("a", "b", "X"))
    alpha = np.asarray(geom.alpha, dtype=np.float64)
    lam = np.asarray(geom.Lambda, dtype=np.float64)
    da, db, dx = d1_even_input(a), d1_even_input(b), d1_even_input(x)
    lambda_metric = da / (2.0*a*a) - db/(a*b) + 2.0/r*(1.0/b - 1.0/a)
    C = lam - lambda_metric
    dC = d1_odd_input(C)
    d2C = np.asarray(grid.cell_second_derivative_fourth(C, parity=-1), dtype=np.float64)
    chi = x*x
    eta = alpha*chi*da/(2.0*a)
    deta = d1_odd_input(eta)
    etaC = eta*C
    # C and eta are odd; eta*C is even and needs the even-input stencil.
    d_etaC = d1_even_input(etaC)
    gamma_m = -3.0*dx/x + 3.0/r + 1.5*db/b
    eta_dC = eta*dC
    d_eta_plus_gamma_eta_times_C = (deta + gamma_m*eta)*C
    gamma_etaC = gamma_m*etaC
    direct = d_etaC + gamma_etaC
    expanded = eta_dC + d_eta_plus_gamma_eta_times_C
    product_remainder = direct - expanded
    # Check arithmetic assembly closure separately from the product-rule
    # remainder itself. The latter is the measured discrete operator defect
    # and must not be forced to roundoff.
    direct_assembly_remainder = direct - (d_etaC + gamma_etaC)
    expanded_assembly_remainder = expanded - (eta_dC + d_eta_plus_gamma_eta_times_C)
    eps_scale = 128.0 * np.finfo(np.float64).eps
    direct_assembly_bound = eps_scale * (
        np.abs(d_etaC) + np.abs(gamma_etaC) + np.abs(direct)
    )
    expanded_assembly_bound = eps_scale * (
        np.abs(eta_dC) + np.abs(d_eta_plus_gamma_eta_times_C) + np.abs(expanded)
    )
    arithmetic_closure = {
        "direct_assembly_remainder_max_abs": float(np.max(np.abs(direct_assembly_remainder))),
        "direct_assembly_bound_max": float(np.max(direct_assembly_bound)),
        "direct_assembly_passed": bool(np.all(np.abs(direct_assembly_remainder) <= direct_assembly_bound)),
        "expanded_assembly_remainder_max_abs": float(np.max(np.abs(expanded_assembly_remainder))),
        "expanded_assembly_bound_max": float(np.max(expanded_assembly_bound)),
        "expanded_assembly_passed": bool(np.all(np.abs(expanded_assembly_remainder) <= expanded_assembly_bound)),
        "bound_definition": "128*float64_eps*sum(abs(assembled operands)); applies only to arithmetic assembly, not direct-minus-expanded product-rule remainder",
    }
    if not arithmetic_closure["direct_assembly_passed"] or not arithmetic_closure["expanded_assembly_passed"]:
        raise ArithmeticError("connection-correction arithmetic assembly failed operand-scaled closure")
    # Keep this source-derived overlap candidate visible; its relation to Brown's scaffold still needs an equation-level reconciliation.
    brown_connection_d1c = (2.0/3.0)*alpha*chi*dC
    arrays = {
        "connection_C_signed": C,
        "connection_D1C_signed": dC,
        "connection_D2C_signed": d2C,
        "connection_eta_signed": eta,
        "connection_D1eta_signed": deta,
        "connection_etaC_signed": etaC,
        "connection_D1_etaC_signed": d_etaC,
        "connection_GammaM_signed": gamma_m,
        "connection_eta_D1C_signed": eta_dC,
        "connection_D1eta_plus_GammaM_eta_times_C_signed": d_eta_plus_gamma_eta_times_C,
        "connection_GammaM_etaC_signed": gamma_etaC,
        "connection_correction_direct_signed": direct,
        "connection_correction_expanded_signed": expanded,
        "connection_product_rule_remainder_signed": product_remainder,
        "connection_direct_assembly_remainder_signed": direct_assembly_remainder,
        "connection_direct_assembly_bound_signed": direct_assembly_bound,
        "connection_expanded_assembly_remainder_signed": expanded_assembly_remainder,
        "connection_expanded_assembly_bound_signed": expanded_assembly_bound,
        "connection_overlap_candidate_2over3_alpha_chi_D1C_signed": brown_connection_d1c,
    }
    if not all(np.all(np.isfinite(v)) for v in arrays.values()):
        raise FloatingPointError("non-finite connection-correction array")
    def cell_record(idx):
        return {"cell": idx, "radius": float(r[idx]),
                **{name: float(values[idx]) for name, values in arrays.items()}}
    cells = {str(idx): cell_record(idx) for idx in (0, 1, 2)}
    fixed_radius = {}
    for target in (2.0, 3.0, 5.0, 8.0):
        idx = int(np.argmin(np.abs(r-target)))
        fixed_radius[str(target)] = {
            "cell": idx, "grid_radius": float(r[idx]),
            **{name: float(arrays[name][idx]) for name in (
                "connection_correction_direct_signed",
                "connection_correction_expanded_signed",
                "connection_product_rule_remainder_signed",
                "connection_overlap_candidate_2over3_alpha_chi_D1C_signed")}
        }
    return {
        "arrays": arrays,
        "report": {
            "formula": "F=eta*C; direct=D1_even(F)+GammaM*F; expanded=eta*D1_odd(C)+(D1_odd(eta)+GammaM*eta)*C",
            "definitions": {
                "C": "Lambda - [D1_odd(a)/(2*a^2) - D1_odd(b)/(a*b) + 2/r*(1/b - 1/a)]",
                "eta": "alpha*X^2*D1_odd(a)/(2*a)",
                "GammaM": "-3*D1_odd(X)/X + 3/r + 3/2*D1_odd(b)/b",
                "product_remainder": "direct - expanded = D1_even(eta*C) - eta*D1_odd(C) - C*D1_odd(eta)"
            },
            "parity": {"C": "odd", "eta": "odd", "eta_times_C": "even",
                       "D1C_and_D1eta": "odd stencil", "D1_etaC": "even stencil"},
            "arithmetic_closure": arithmetic_closure,
            "signed_summary": {name: summary(arrays[name]) for name in (
                "connection_correction_direct_signed",
                "connection_correction_expanded_signed",
                "connection_product_rule_remainder_signed",
                "connection_overlap_candidate_2over3_alpha_chi_D1C_signed")},
            "cells_0_1_2": cells,
            "nearest_fixed_radii": fixed_radius,
            "overlap_policy": "report separately; compare against Brown scaffold equation before any combined prediction",
            "closure_claim": False
        }
    }


def run(output_dir: Path):
    source = momentum.provenance()
    stress_source = stress.verify_provenance()
    kernel = momentum.V55ProductionKernel()
    initial = kernel.initialize(resolution=RESOLUTION, r_max=R_MAX, include_radiation=True)
    initial_digest = require_digest(initial, INITIAL_SHA256, "initial-state")
    accepted = kernel.step(initial, DT)
    accepted_digest = require_digest(accepted, ACCEPTED_SHA256, "accepted-state")
    rates = momentum.rhs_direction(kernel, accepted)
    plus = momentum.shifted_state(kernel, accepted, rates, EPSILON)
    minus = momentum.shifted_state(kernel, accepted, rates, -EPSILON)
    base_digest = momentum.state_digest(accepted)
    plus_digest = momentum.state_digest(plus)
    minus_digest = momentum.state_digest(minus)
    require_digest(accepted, ACCEPTED_SHA256, "shared base state")
    # Negative control: the exact same state must be rejected against a wrong digest.
    try:
        require_digest(accepted, "0" * 64, "negative-control")
    except RuntimeError:
        wrong_digest_rejected = True
    else:
        wrong_digest_rejected = False
    if not wrong_digest_rejected:
        raise RuntimeError("negative control failed: mismatched state digest was accepted")
    if plus_digest == base_digest or minus_digest == base_digest:
        raise RuntimeError("directional states unexpectedly equal the accepted base state")
    if not np.all(np.asarray(accepted.geometry.beta) == 0.0):
        raise RuntimeError("observer requires zero shift on the accepted slice")

    # Both observers consume the same accepted/base/plus/minus objects.
    base_constraint = momentum.observe_constraint(accepted)
    plus_constraint = momentum.observe_constraint(plus)
    minus_constraint = momentum.observe_constraint(minus)
    mdot_total = (np.asarray(plus_constraint["M_total"]) -
                  np.asarray(minus_constraint["M_total"])) / (2*EPSILON)
    mdot_geom = (np.asarray(plus_constraint["M_geom"]) -
                 np.asarray(minus_constraint["M_geom"])) / (2*EPSILON)
    div = stress.divergence_tangent_from_states(accepted, plus, minus, EPSILON)
    sector = stress.sector_resolved_divergence_tangent_from_states(accepted, plus, minus, EPSILON)

    # Independent Brown Eq. 11b structural predictor evaluated on this same
    # accepted geometry. It is intentionally vacuum-form: the residual exposes
    # omitted matter/source and discrete product-rule contributions; it is not
    # labelled as a complete matter-coupled propagation RHS.
    _, vacuum, _ = momentum.adapter.vendor_modules()
    geometric_constraints = vacuum.constraints(accepted.grid, accepted.geometry)
    brown_terms, brown_prediction = brown.brown_lower_momentum_rhs(
        accepted.grid, accepted.geometry, geometric_constraints
    )
    brown_residual = mdot_geom - brown_prediction
    hgrad = np.asarray(brown_terms["alpha_DrH_over_6"], dtype=np.float64)
    wrong_sign_prediction = brown_prediction - 2.0*hgrad
    wrong_sign_delta = brown_prediction - wrong_sign_prediction
    operand_scale = sum((np.abs(v) for v in brown_terms.values()),
                        np.zeros_like(brown_prediction))
    negative_control_bound = 128.0 * np.finfo(np.float64).eps * operand_scale
    negative_control_detected = bool(np.any(np.abs(wrong_sign_delta) > negative_control_bound))
    if not negative_control_detected:
        raise RuntimeError("Brown H-gradient wrong-sign negative control was not detected")

    # Separate the explicit matter-source contribution to the geometric
    # constraint tangent. The projection is independently rebuilt sector by
    # sector; its summed Aa/K sources are cross-checked against the pinned
    # production L3-minus-vacuum L3 source delta before use.
    projections = stress.independent_total_projections(accepted)
    if not projections["projection_crosscheck_passed"]:
        raise RuntimeError("independent stress projection failed canonical cross-check")
    matter_source = matter_source_momentum_tangent(accepted, projections)
    source_rates = rates["geometry_stage_terms"]["primary_l3"]
    vacuum_l3 = vacuum.primary_l3_rhs(accepted.grid, accepted.geometry)
    production_aa_source = np.asarray(source_rates["Aa"]) - np.asarray(vacuum_l3["Aa"])
    production_k_source = np.asarray(source_rates["K"]) - np.asarray(vacuum_l3["K"])
    geom = accepted.geometry
    expected_aa_source = -(16.0 * np.pi / 3.0) * np.asarray(geom.alpha) * (
        np.asarray(projections["pr"]) - np.asarray(projections["pt"])
    )
    expected_k_source = 4.0 * np.pi * np.asarray(geom.alpha) * (
        np.asarray(projections["rho"]) + np.asarray(projections["pr"])
        + 2.0 * np.asarray(projections["pt"])
    )
    aa_source_delta = production_aa_source - expected_aa_source
    k_source_delta = production_k_source - expected_k_source
    source_projection_crosscheck = {
        "Aa_source_max_abs_delta": float(np.max(np.abs(aa_source_delta))),
        "K_source_max_abs_delta": float(np.max(np.abs(k_source_delta))),
        "passed": bool(
            np.allclose(production_aa_source, expected_aa_source, rtol=5.0e-12, atol=2.0e-16)
            and np.allclose(production_k_source, expected_k_source, rtol=5.0e-12, atol=2.0e-16)
        ),
    }
    if not source_projection_crosscheck["passed"]:
        raise RuntimeError("independent Aa/K matter-source cross-check failed")
    brown_plus_source = brown_prediction + matter_source["total"]
    residual_after_source = mdot_geom - brown_plus_source
    connection = connection_correction(accepted)
    r = np.asarray(accepted.grid.centers, dtype=np.float64)
    # Diagnostic-only candidate: retain the vendor-derived correction as a
    # separate signed addition. The Brown-overlap candidate is NOT subtracted
    # or folded in until its equation-level relationship is reconciled.
    connection_candidate = brown_prediction + connection["arrays"]["connection_correction_direct_signed"]
    connection_candidate_residual = mdot_geom - connection_candidate
    matter_connection_candidate = brown_plus_source + connection["arrays"]["connection_correction_direct_signed"]
    matter_connection_candidate_residual = mdot_geom - matter_connection_candidate
    connection["report"]["scaffold_comparison"] = {
        "original_brown_scaffold_signed": summary(brown_prediction),
        "measured_geometric_tangent_signed": summary(mdot_geom),
        "original_brown_residual_measured_minus_scaffold_signed": summary(brown_residual),
        "brown_plus_connection_correction_candidate_signed": summary(connection_candidate),
        "residual_after_connection_candidate_signed": summary(connection_candidate_residual),
        "brown_plus_explicit_matter_source_plus_connection_candidate_signed": summary(matter_connection_candidate),
        "residual_after_matter_and_connection_candidate_signed": summary(matter_connection_candidate_residual),
        "cells_0_1_2": {
            str(i): {
                "radius": float(accepted.grid.centers[i]),
                "measured_geometric_tangent": float(mdot_geom[i]),
                "original_brown_scaffold": float(brown_prediction[i]),
                "original_residual": float(brown_residual[i]),
                "connection_direct": float(connection["arrays"]["connection_correction_direct_signed"][i]),
                "connection_expanded": float(connection["arrays"]["connection_correction_expanded_signed"][i]),
                "connection_product_rule_remainder": float(connection["arrays"]["connection_product_rule_remainder_signed"][i]),
                "overlap_candidate": float(connection["arrays"]["connection_overlap_candidate_2over3_alpha_chi_D1C_signed"][i]),
                "brown_plus_connection_candidate": float(connection_candidate[i]),
                "residual_after_connection_candidate": float(connection_candidate_residual[i]),
                "brown_plus_matter_and_connection_candidate": float(matter_connection_candidate[i]),
                "residual_after_matter_and_connection_candidate": float(matter_connection_candidate_residual[i]),
            } for i in (0, 1, 2)
        },
        "nearest_fixed_radii": {
            str(target): {
                "cell": int(np.argmin(np.abs(r - target))),
                "grid_radius": float(r[int(np.argmin(np.abs(r - target)))]),
                "measured_geometric_tangent": float(mdot_geom[int(np.argmin(np.abs(r - target)))]),
                "original_brown_scaffold": float(brown_prediction[int(np.argmin(np.abs(r - target)))]),
                "original_residual": float(brown_residual[int(np.argmin(np.abs(r - target)))]),
                "connection_direct": float(connection["arrays"]["connection_correction_direct_signed"][int(np.argmin(np.abs(r - target)))]),
                "connection_expanded": float(connection["arrays"]["connection_correction_expanded_signed"][int(np.argmin(np.abs(r - target)))]),
                "connection_product_rule_remainder": float(connection["arrays"]["connection_product_rule_remainder_signed"][int(np.argmin(np.abs(r - target)))]),
                "overlap_candidate": float(connection["arrays"]["connection_overlap_candidate_2over3_alpha_chi_D1C_signed"][int(np.argmin(np.abs(r - target)))]),
                "residual_after_connection_candidate": float(connection_candidate_residual[int(np.argmin(np.abs(r - target)))]),
            } for target in (2.0, 3.0, 5.0, 8.0)
        },
        "overlap_policy": "overlap candidate remains separate; candidate residuals are diagnostic and are not closure or acceptance claims",
        "closure_claim": False,
    }

    after_digest = momentum.state_digest(accepted)
    if after_digest != base_digest:
        raise RuntimeError("shared accepted state mutated during observation")

    r = np.asarray(accepted.grid.centers, dtype=np.float64)
    output_dir.mkdir(parents=True, exist_ok=True)
    arrays = {
        "radius": r,
        "Mdot_total_signed": mdot_total,
        "Mdot_geometric_signed": mdot_geom,
        "brown_vacuum_structural_prediction_signed": np.asarray(brown_prediction),
        "brown_structural_residual_signed": np.asarray(brown_residual),
        "brown_H_gradient_term_signed": hgrad,
        "brown_wrong_sign_mutation_signed": wrong_sign_delta,
        "matter_source_Mdot_signed": matter_source["total"],
        "brown_plus_matter_source_prediction_signed": brown_plus_source,
        "residual_after_matter_source_signed": residual_after_source,
        "brown_plus_connection_correction_candidate_signed": connection_candidate,
        "residual_after_connection_candidate_signed": connection_candidate_residual,
        "brown_plus_matter_source_plus_connection_candidate_signed": matter_connection_candidate,
        "residual_after_matter_and_connection_candidate_signed": matter_connection_candidate_residual,
        **connection["arrays"],
        **{f"matter_source_Mdot_{name}_signed": value for name, value in matter_source["pieces"].items()},
        "div_T_t_signed": np.asarray(div["div_T_t"]),
        "div_T_r_signed": np.asarray(div["div_T_r"]),
        "div_T_t_fv_signed": np.asarray(div["div_T_t_fv"]),
        "div_T_r_fv_signed": np.asarray(div["div_T_r_fv"]),
    }
    np.savez_compressed(output_dir / "shared_state_signed_arrays.npz", **arrays)
    report = {
        "schema": "issue4_shared_state_momentum_total_stress_divergence_v3",
        "status": "completed_diagnostic_only",
        "provenance": {"momentum_observer": source, "total_stress_observer": stress_source},
        "configuration": {"resolution": RESOLUTION, "r_max": R_MAX, "dt": DT,
                          "epsilon": EPSILON, "production_steps_advanced": 1},
        "state_pairing": {"initial_sha256": initial_digest, "accepted_sha256": base_digest,
                          "plus_sha256": plus_digest, "minus_sha256": minus_digest,
                          "accepted_after_observation_sha256": after_digest,
                          "accepted_unchanged": after_digest == base_digest,
                          "all_observers_used_same_state_objects": True,
                          "wrong_digest_negative_control_rejected": wrong_digest_rejected},
        "signed_residuals": {
            "Mdot_total_constraint": summary(mdot_total),
            "Mdot_geometric_constraint": summary(mdot_geom),
            "Brown_vacuum_structural_prediction": summary(brown_prediction),
            "Brown_structural_residual_measured_minus_predicted": summary(brown_residual),
            "Brown_H_gradient_term": summary(hgrad),
            "Brown_wrong_sign_mutation": summary(wrong_sign_delta),
            "explicit_matter_source_Mdot_geometric": summary(matter_source["total"]),
            "Brown_plus_matter_source_prediction": summary(brown_plus_source),
            "residual_after_explicit_matter_source": summary(residual_after_source),
            "matter_source_Mdot_by_sector": {name: summary(value) for name, value in matter_source["pieces"].items()},
            "matter_source_sector_decomposition_gap": summary(matter_source["decomposition_gap"]),
            "div_T_t_cell_centered": summary(div["div_T_t"]),
            "div_T_r_cell_centered": summary(div["div_T_r"]),
            "div_T_t_finite_volume": summary(div["div_T_t_fv"]),
            "div_T_r_finite_volume": summary(div["div_T_r_fv"]),
            "sector_decomposition": sector["sector_decomposition"]["reconstruction"],
        },
        "brown_structural_crosscheck": {
            "formula": "beta*D_r(M)+(D_r beta)*M+alpha*K*M-(D_r alpha)*H/3"
                       "+alpha*D_r(H)/6+2*alpha*chi*D_r^2(C)/3"
                       "+[2*(D_r alpha)*chi/3-alpha*(D_r chi)/3"
                       "+alpha*D_r ln(g_theta_theta)]*D_r(C)",
            "hamiltonian_slot": "independent metric-only warped-product Ricci reconstruction",
            "wrong_sign_H_gradient_control_detected": negative_control_detected,
            "scope": "vacuum structural terms plus separately measured explicit matter source; discrete product-rule/operator residual remains",
            "full_matter_coupled_identity": False,
        },
        "connection_correction": connection["report"],
        "matter_source_crosscheck": {
            **source_projection_crosscheck,
            "independent_projection_crosscheck_passed": projections["projection_crosscheck_passed"],
            "sector_decomposition_passed": bool(np.all(
                np.abs(matter_source["decomposition_gap"]) <= matter_source["decomposition_bound"]
            )),
            "formula": "D_r(Aa_source) - 2/3 D_r(K_source) + [6 chi_r + 3/r + 3/2 D_r(b)/b]*Aa_source",
            "source_terms": "Aa_source=-(16*pi/3)*alpha*(pr-pt); K_source=4*pi*alpha*(rho+pr+2*pt)",
        },
        "raw_arrays": {"file": "shared_state_signed_arrays.npz",
                       "names": sorted(arrays), "all_finite": True},
        "gates": {"shared_state_digest_match": True,
                  "wrong_digest_negative_control_passed": wrong_digest_rejected,
                  "observer_immutability_passed": True,
                  "brown_wrong_sign_negative_control_passed": negative_control_detected,
                  "matter_source_crosscheck_passed": source_projection_crosscheck["passed"],
                  "matter_source_sector_decomposition_passed": bool(np.all(
                      np.abs(matter_source["decomposition_gap"]) <= matter_source["decomposition_bound"]
                  )),
                  "independent_full_propagation_identity_closed": False,
                  "scientific_admission": False, "acceptance_oracle_satisfied": False},
        "interpretation": (
            "Shared-state extension: signed measured geometric and total momentum tangents, "
            "independent total-stress divergence, and Brown vacuum structural terms are now "
            "evaluated on the same accepted state and directional pair. The Brown residual is "
            "not a closure claim: matter-source terms and discrete operator/product-rule terms "
            "remain to be independently separated. Signed connection-correction and product-rule "
            "remainder arrays are now emitted separately; the correction is not folded into the "
            "predictor pending Brown-overlap reconciliation."
        )
    }
    (output_dir / "shared_state_report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output-dir", type=Path,
                   default=Path("results/accepted-state-raw-momentum-witness"))
    args = p.parse_args()
    report = run(args.output_dir)
    print(json.dumps({"status": report["status"], "state_pairing": report["state_pairing"],
                      "signed_residuals": {k:v for k,v in report["signed_residuals"].items()
                                           if k != "sector_decomposition"},
                      "brown_structural_crosscheck": report["brown_structural_crosscheck"],
                      "gates": report["gates"]}, sort_keys=True), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
