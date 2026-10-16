#!/usr/bin/env python3
"""Reproducible three-arm origin-integral / five-cell Newton discriminator.

This copied-state audit runs (1) the pinned production initializer, (2) its
principled continuous-origin-half-cell-integral counterfactual, and (3) the
same principled B profile with a finite-difference Newton solve of the exact
vendor Hamiltonian residual on cells 0..4, varying B[0:5] only.

No production files are changed and no evolution steps are advanced. The report
includes signed H/M/C_Lambda/determinant and semi-discrete RHS tangent witnesses.
The independently normalized full constraint-propagation RHS is deliberately
not inferred from the measured tangent and remains the scientific acceptance
gate. The earlier supplied hybrid state lacked an exact provenance bundle; this
creates a new reproducible hybrid family, not a claim of numerical identity.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import traceback
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from engine.production_kernel import ProductionState, V55ProductionKernel
from engine.stress_energy import assemble_total_stress_energy
from engine.v55_matter import (
    RHO_B_PRESENT,
    RHO_DM_PRESENT,
    RHO_R_PRESENT,
    initialize_from_archive,
)
from engine import v55_pirk_adapter as adapter
from tools.audit_continuous_origin_half_cell_oracle import (
    analytic_half_cell_envelope,
    continuous_reference,
)
from tools.audit_current_source_momentum_constraint_tangent import (
    EXPECTED_VENDOR,
    SOURCE_BASELINE,
    GEOMETRY_RATE_NAMES,
    SCALAR_NAMES,
    provenance,
    rhs_direction,
    shifted_state,
    observe_constraint,
    tangent_for_epsilon as momentum_tangent_for_epsilon,
    compact_tangent_result,
    state_digest,
)

PI = float(math.pi)
SCHEMA = "hybrid_origin_integral_five_cell_newton_signed_tangent_v1"
DEFAULT_RESOLUTIONS = (60, 80, 120, 160)
DEFAULT_R_MAX = 80.0
AMPLITUDE = 0.01
WIDTH = 7.0
D_AMPLITUDE = 1.0e-10
PATCH_CELLS = 5
CONTINUOUS_SUBSTEPS = (256, 512, 1024, 2048, 4096, 8192)
NEWTON_MAX_ITER = 12
NEWTON_TOL = 1.0e-9
NEWTON_MAX_CONDITION = 1.0e12
NEWTON_RELATIVE_FD_EPS = 1.0e-7
DIRECTIONAL_EPSILONS = (2.0e-5, 1.0e-5)
FIXED_RADII = (1.0, 2.0, 3.0, 5.0, 8.0, 10.0)
DEFAULT_OUT = Path("artifacts/hybrid-origin-integral-patch-newton")


def finite(label: str, values: Any) -> np.ndarray:
    arr = np.asarray(values, dtype=np.float64)
    if not np.all(np.isfinite(arr)):
        raise FloatingPointError(f"non-finite values in {label}")
    return arr


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def initializer_source(state: ProductionState, B: np.ndarray) -> np.ndarray:
    """Exact copied six-pass initializer source formula with current B iterate."""
    grid, f = state.grid, state.scalars
    r = finite("source radii", grid.centers)
    B = finite("source B", B)
    S, PS = finite("S", f.S), finite("PS", f.PS)
    D, PD = finite("D", f.D), finite("PD", f.PD)
    phi, Pi = finite("phi", f.phi), finite("Pi", f.Pi)
    d1 = lambda values: grid.cell_derivative_fourth(values, parity=1)
    Sp, Dp = finite("D_r S", d1(S)), finite("D_r D", d1(D))
    current = -(PD * Dp)
    # Kt = K/3 - Aa/2 follows from the spherical initializer identities.
    Kt = finite(
        "K_t",
        np.asarray(state.geometry.K, dtype=np.float64) / 3.0
        - np.asarray(state.geometry.Aa, dtype=np.float64) / 2.0,
    )
    Ktp = finite("D_r K_t", d1(Kt))
    from engine.scalar_system import cosmos_potential
    rho_S = 0.5 * (PS * PS + B * Sp * Sp) + 0.5 * S * S
    rho_D = 0.5 * (PD * PD + B * Dp * Dp) - 0.5 * D * D
    rho_phi = (0.5 * Pi * Pi + cosmos_potential(phi)) / (8.0 * PI)
    rho_matter = np.full_like(
        r, (RHO_DM_PRESENT + RHO_B_PRESENT + RHO_R_PRESENT) / (8.0 * PI)
    )
    rho = rho_S + rho_D + rho_phi + rho_matter
    K_radial = Kt + r * Ktp + 4.0 * PI * r * current
    source = (
        1.0
        - 8.0 * PI * r * r * rho
        + 2.0 * r * r * K_radial * Kt
        + r * r * Kt * Kt
    )
    return finite("copied initializer source", source)


def map_source_to_B(
    grid: Any, source: np.ndarray, first_half_cell_moment: float | None = None
) -> np.ndarray:
    """Keep production inter-centre trapezoids; optionally replace only cell 0."""
    r = finite("B-map radii", grid.centers)
    src = finite("B-map source", source)
    if src.shape != r.shape or r.size == 0:
        raise ValueError("source shape must match the nonempty radial grid")
    dr = float(grid.dr)
    if not math.isfinite(dr) or dr <= 0.0:
        raise ValueError("grid spacing must be finite and positive")
    first = 0.5 * dr * float(src[0]) if first_half_cell_moment is None else float(first_half_cell_moment)
    if not math.isfinite(first) or first <= 0.0:
        raise ValueError("first half-cell moment must be finite and positive")
    integrated = np.cumsum(np.r_[first, 0.5 * dr * (src[1:] + src[:-1])])
    return finite("mapped metric B", integrated / r)


def source_iteration(
    state: ProductionState,
    first_half_cell_moment: float | None,
    iterations: int = 6,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Rerun the six source/B updates, changing only the origin half-cell moment."""
    B = np.ones_like(np.asarray(state.grid.centers, dtype=np.float64))
    history = []
    for it in range(1, iterations + 1):
        source = initializer_source(state, B)
        B_next = map_source_to_B(state.grid, source, first_half_cell_moment)
        if not np.all(np.isfinite(B_next)) or np.any(B_next <= 0.0):
            raise FloatingPointError(f"B iteration {it} produced nonpositive/nonfinite values")
        history.append({
            "iteration": it,
            "max_abs_B_update": float(np.max(np.abs(B_next - B))),
            "half_cell_moment_used": (
                0.5 * float(state.grid.dr) * float(source[0])
                if first_half_cell_moment is None else float(first_half_cell_moment)
            ),
            "source_cell_0": float(source[0]),
        })
        B = B_next
    return B, {
        "iterations": int(iterations),
        "history": history,
        "final_max_abs_B_update": history[-1]["max_abs_B_update"],
        "B_min": float(np.min(B)),
        "B_max": float(np.max(B)),
    }


def rebuilt_state_from_B(
    kernel: V55ProductionKernel,
    template: ProductionState,
    B_metric: np.ndarray,
) -> ProductionState:
    """Rebuild metric, archive matter and CMC lapse on a copied B profile."""
    grid = template.grid
    r = finite("rebuild radii", grid.centers)
    B_metric = finite("rebuild metric B", B_metric)
    if B_metric.shape != r.shape or np.any(B_metric <= 0.0):
        raise ValueError("metric B must be positive and match the grid")

    geom = template.geometry.copy()
    A = 1.0 / np.sqrt(B_metric)
    geom.a = finite("rebuilt a", A ** (4.0 / 3.0))
    geom.b = finite("rebuilt b", A ** (-2.0 / 3.0))
    geom.X = finite("rebuilt X", A ** (-1.0 / 3.0))
    geom.Lambda = finite(
        "rebuilt Lambda",
        grid.cell_derivative_fourth(geom.a, parity=1) / (2.0 * geom.a * geom.a)
        - grid.cell_derivative_fourth(geom.b, parity=1) / (geom.a * geom.b)
        + 2.0 / r * (1.0 / geom.b - 1.0 / geom.a),
    )
    # Match initializer ordering: make archive matter on alpha=1, beta=0,
    # then refresh the same CMC lapse on the candidate metric and matter.
    geom.alpha = np.ones_like(r)
    geom.beta.fill(0.0)
    geom.B.fill(0.0)  # gauge-driver B, not the metric integration variable
    matter = initialize_from_archive(grid, geom, include_radiation=True)
    alpha, _ = kernel._solve_lapse_for_state(
        grid, geom, template.scalars, matter, None
    )
    geom.alpha = finite("rebuilt CMC alpha", alpha).copy()
    geom.beta.fill(0.0)
    geom.B.fill(0.0)
    if np.min(geom.alpha) <= 0.0:
        raise FloatingPointError("rebuilt CMC lapse is not positive")
    return ProductionState(
        grid=grid,
        geometry=geom,
        scalars=template.scalars,
        matter=matter,
        t=template.t,
        tau=template.tau,
        e_folds=template.e_folds,
    )


def constraints(state: ProductionState, vacuum: Any) -> dict[str, np.ndarray]:
    grid, geom = state.grid, state.geometry
    r = finite("constraint radii", grid.centers)
    raw = vacuum.constraints(grid, geom)
    total = assemble_total_stress_energy(grid, geom, state.scalars, state.matter)
    H = finite("H", np.asarray(raw["hamiltonian"]) - 16.0 * PI * np.asarray(total.rho))
    M = finite("M", np.asarray(raw["momentum"]) - 8.0 * PI * np.asarray(total.j))
    d1 = lambda values: grid.cell_derivative_fourth(values, parity=1)
    a = finite("constraint a", geom.a)
    b = finite("constraint b", geom.b)
    connection_from_metric = (
        d1(a) / (2.0 * a * a)
        - d1(b) / (a * b)
        + 2.0 / r * (1.0 / b - 1.0 / a)
    )
    C = finite("C_Lambda", np.asarray(geom.Lambda) - connection_from_metric)
    determinant = finite("determinant", a * b * b - 1.0)
    return {"H": H, "M": M, "C_Lambda": C, "determinant": determinant}


def witness(values: np.ndarray, r: np.ndarray, *, index_offset: int = 0) -> dict[str, Any]:
    values = finite("witness profile", values)
    r = finite("witness radii", r)
    if values.shape != r.shape:
        raise ValueError("witness arrays must have matching shape")
    n = min(PATCH_CELLS, values.size)
    annulus = (r >= 5.0) & (r <= 10.0)
    peaks = np.flatnonzero(np.abs(values) == np.max(np.abs(values)))
    peak = int(peaks[0])
    return {
        "cell_0_to_4": [
            {"i": i + int(index_offset), "r": float(r[i]), "signed": float(values[i]), "abs": float(abs(values[i]))}
            for i in range(n)
        ],
        "fixed_radius_signed": {
            str(target): float(np.interp(target, r, values))
            for target in FIXED_RADII if r[0] <= target <= r[-1]
        },
        "full_max_abs": float(abs(values[peak])),
        "full_peak_signed": float(values[peak]),
        "full_peak_radius": float(r[peak]),
        "full_rms": float(np.sqrt(np.mean(values * values))),
        "annulus_r5_to10_max_abs": float(np.max(np.abs(values[annulus]))) if np.any(annulus) else None,
        "annulus_r5_to10_rms": float(np.sqrt(np.mean(values[annulus] ** 2))) if np.any(annulus) else None,
    }


def exterior_witness(values: np.ndarray, radii: np.ndarray) -> dict[str, Any]:
    """Witness the untouched cells with their true global cell indices."""
    out = witness(values, radii, index_offset=PATCH_CELLS)
    out["cells_5_plus_first_five"] = out.pop("cell_0_to_4")
    return out


def solve_five_cell_hybrid(
    kernel: V55ProductionKernel,
    template: ProductionState,
    B_seed: np.ndarray,
    vacuum: Any,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Newton solve H[0:5]=0 by varying B[0:5] only; fail closed."""
    full_B_seed = finite("hybrid B seed", B_seed).copy()
    if full_B_seed.size < PATCH_CELLS + 1:
        raise ValueError("grid is too small for five-cell patch and exterior")
    x = full_B_seed[:PATCH_CELLS].copy()
    history = []
    max_condition = 0.0

    def evaluate(candidate_x: np.ndarray):
        if not np.all(np.isfinite(candidate_x)) or np.any(candidate_x <= 0.0):
            raise FloatingPointError("Newton trial metric B patch must remain finite and positive")
        candidate_B = full_B_seed.copy()
        candidate_B[:PATCH_CELLS] = candidate_x
        candidate_state = rebuilt_state_from_B(kernel, template, candidate_B)
        H = constraints(candidate_state, vacuum)["H"]
        return finite("Newton H patch", H[:PATCH_CELLS]), candidate_B

    H, full_B = evaluate(x)
    for iteration in range(NEWTON_MAX_ITER + 1):
        norm0 = float(np.max(np.abs(H)))
        record = {
            "iteration": iteration,
            "max_abs_H_cells_0_to_4": norm0,
            "B_patch": [float(v) for v in x],
            "H_patch": [float(v) for v in H],
            "condition_number": None,
            "accepted_line_search_scale": 0.0,
        }
        history.append(record)
        if norm0 <= NEWTON_TOL:
            return full_B, {
                "converged": True,
                "iterations_used": iteration,
                "tolerance": NEWTON_TOL,
                "patch_cell_indices": list(range(PATCH_CELLS)),
                "unknowns": "B[0:5]; B[5:] held exactly at the principled-origin-integral profile",
                "max_condition_number_seen": max_condition,
                "history": history,
                "final_max_abs_H_patch_before_final_rebuild": norm0,
                "B_patch_final": [float(v) for v in x],
            }
        if iteration == NEWTON_MAX_ITER:
            break

        J = np.zeros((PATCH_CELLS, PATCH_CELLS), dtype=np.float64)
        for j in range(PATCH_CELLS):
            trial = x.copy()
            delta = NEWTON_RELATIVE_FD_EPS * max(abs(float(x[j])), 1.0e-8)
            trial[j] += delta
            H_trial, _ = evaluate(trial)
            J[:, j] = (H_trial - H) / delta
        if not np.all(np.isfinite(J)):
            raise FloatingPointError("Newton Jacobian contains non-finite values")
        condition = float(np.linalg.cond(J))
        record["condition_number"] = condition
        max_condition = max(max_condition, condition)
        if not math.isfinite(condition) or condition > NEWTON_MAX_CONDITION:
            raise RuntimeError(f"Newton Jacobian condition {condition:.6e} exceeds limit")
        delta_x = np.linalg.solve(J, -H)
        if not np.all(np.isfinite(delta_x)):
            raise FloatingPointError("Newton step is non-finite")

        accepted = False
        scale = 1.0
        for _ in range(24):
            trial = x + scale * delta_x
            if np.all(np.isfinite(trial)) and np.all(trial > 0.0):
                H_trial, B_trial = evaluate(trial)
                trial_norm = float(np.max(np.abs(H_trial)))
                if trial_norm < norm0:
                    x, H, full_B = trial, H_trial, B_trial
                    record["accepted_line_search_scale"] = scale
                    record["accepted_trial_max_abs_H_patch"] = trial_norm
                    accepted = True
                    break
            scale *= 0.5
        if not accepted:
            raise RuntimeError("Newton line search failed to reduce max |H[0:5]|")
    raise RuntimeError(
        f"Newton did not converge in {NEWTON_MAX_ITER} iterations; "
        f"last max |H[0:5]|={float(np.max(np.abs(H))):.6e}"
    )


def _constraint_propagation_candidate(
    state: ProductionState,
    rates: dict[str, Any],
    base_constraints: dict[str, np.ndarray],
    plus_state: ProductionState,
    minus_state: ProductionState,
    epsilon: float,
    measured_mdot: np.ndarray,
) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    """Candidate Brown-spherical BSSN M-propagation law + covariant T-divergence.

    This is deliberately an independent formula evaluation, not the finite-
    difference constraint tangent. Because the source contains matter/exchange,
    compute Q_r = nabla_mu T^mu_r from the projected total stress tensor and
    the same semi-discrete tangent state. The sign/index mapping is an explicit
    diagnostic candidate and remains non-admitted until its normalization review.
    """
    grid, geom = state.grid, state.geometry
    r = finite("propagation radii", grid.centers)
    alpha = finite("propagation alpha", geom.alpha)
    a = finite("propagation a", geom.a)
    b = finite("propagation b", geom.b)
    X = finite("propagation X", geom.X)
    K = finite("propagation K", geom.K)
    d1 = lambda value, parity=1: finite(
        "propagation D1", grid.cell_derivative_fourth(value, parity=parity)
    )
    d2 = lambda value, parity=-1: finite(
        "propagation D2", grid.cell_second_derivative_fourth(value, parity=parity)
    )
    H = finite("propagation H", base_constraints["H"])
    M = finite("propagation M", base_constraints["M"])
    C = finite("propagation C_Lambda", base_constraints["C_Lambda"])
    alpha_r = d1(alpha, 1)
    chi = X * X
    chi_r = d1(chi, 1)
    H_r = d1(H, 1)
    C_r = d1(C, -1)
    C_rr = d2(C, -1)
    conformal_angular_log_derivative = d1(b, 1) / b + 2.0 / r

    term_KM = alpha * K * M
    term_alpha_H = -(alpha_r / 3.0) * H
    term_H_gradient = (alpha / 6.0) * H_r
    term_C_second = (2.0 / 3.0) * alpha * chi * C_rr
    c_first_coeff = (
        (2.0 / 3.0) * alpha_r * chi
        - (1.0 / 3.0) * alpha * chi_r
        + alpha * chi * conformal_angular_log_derivative
    )
    term_C_first = c_first_coeff * C_r
    brown_conserved_source = (
        term_KM + term_alpha_H + term_H_gradient + term_C_second + term_C_first
    )

    total = assemble_total_stress_energy(
        grid, state.geometry, state.scalars, state.matter
    )
    plus_total = assemble_total_stress_energy(
        grid, plus_state.geometry, plus_state.scalars, plus_state.matter
    )
    minus_total = assemble_total_stress_energy(
        grid, minus_state.geometry, minus_state.scalars, minus_state.matter
    )
    j = finite("total covariant radial momentum", total.j)
    j_r_dot = finite(
        "radial momentum time tangent",
        (np.asarray(plus_total.j) - np.asarray(minus_total.j)) / (2.0 * epsilon),
    )

    a_dot = finite("a RHS for stress divergence", rates["geometry"]["a"])
    b_dot = finite("b RHS for stress divergence", rates["geometry"]["b"])
    X_dot = finite("X RHS for stress divergence", rates["geometry"]["X"])
    gamma_rr_log_dot = a_dot / a - 2.0 * X_dot / X
    areal_radius_log_dot = 0.5 * b_dot / b - X_dot / X
    areal_radius_log_r = 1.0 / r + 0.5 * d1(b, 1) / b - d1(X, 1) / X

    # For beta=0 and ds^2=-alpha^2 dt^2 + gamma_rr dr^2 + R^2 dOmega^2,
    # Q_r = nabla_mu T^mu_r. This uses T^t_r=j_r/alpha, T^r_r=p_r.
    Q_r = (
        j_r_dot / alpha
        + (j / alpha) * (
            0.5 * gamma_rr_log_dot + 2.0 * areal_radius_log_dot
        )
        + d1(np.asarray(total.pr), 1)
        + (alpha_r / alpha) * (np.asarray(total.rho) + np.asarray(total.pr))
        + 2.0 * areal_radius_log_r * (np.asarray(total.pr) - np.asarray(total.pt))
    )
    Q_r = finite("total stress-energy radial divergence", Q_r)
    matter_source = -8.0 * math.pi * alpha * Q_r
    predicted = finite(
        "candidate predicted M tangent",
        brown_conserved_source + matter_source,
    )
    residual = finite(
        "measured-minus-predicted M tangent",
        np.asarray(measured_mdot, dtype=np.float64) - predicted,
    )
    scale = max(
        float(np.sqrt(np.mean(np.asarray(measured_mdot, dtype=np.float64) ** 2))),
        float(np.sqrt(np.mean(predicted * predicted))),
        np.finfo(np.float64).tiny,
    )
    report = {
        "formula": (
            "Brown spherical BSSN momentum-constraint propagation (beta=0) "
            "with H/M replaced by the repository total-stress constraints, "
            "C by C_Lambda, chi=X^2, and an explicitly computed "
            "-8*pi*alpha*nabla_mu(T^mu_r) source."
        ),
        "normalization_reviewed_for_scientific_admission": False,
        "matter_conservation_source_sign_and_index_mapping": (
            "candidate -8*pi*alpha*Q_r; output is a discriminator, not an accepted identity"
        ),
        "brown_conserved_source": witness(brown_conserved_source, r),
        "matter_source_minus_8pi_alpha_Qr": witness(matter_source, r),
        "Q_r_total_stress_divergence": witness(Q_r, r),
        "predicted_Mdot": witness(predicted, r),
        "measured_Mdot": witness(measured_mdot, r),
        "measured_minus_predicted": witness(residual, r),
        "relative_rms_error_full_domain": float(
            np.sqrt(np.mean(residual * residual)) / scale
        ),
        "relative_rms_error_r_ge_5": float(
            np.sqrt(np.mean(residual[r >= 5.0] ** 2))
            / max(
                float(np.sqrt(np.mean(np.asarray(measured_mdot)[r >= 5.0] ** 2))),
                float(np.sqrt(np.mean(predicted[r >= 5.0] ** 2))),
                np.finfo(np.float64).tiny,
            )
        ),
        "finite": bool(
            np.all(np.isfinite(predicted))
            and np.all(np.isfinite(residual))
            and np.all(np.isfinite(Q_r))
        ),
        "acceptance_oracle_satisfied": False,
    }
    arrays = {
        "propagation_brown_conserved_source": brown_conserved_source,
        "propagation_Qr_total_stress_divergence": Q_r,
        "propagation_matter_source_minus_8pi_alpha_Qr": matter_source,
        "propagation_Mdot_predicted": predicted,
        "propagation_Mdot_measured_minus_predicted": residual,
        "propagation_term_KM": term_KM,
        "propagation_term_alpha_H": term_alpha_H,
        "propagation_term_H_gradient": term_H_gradient,
        "propagation_term_connection_second": term_C_second,
        "propagation_term_connection_first": term_C_first,
    }
    return report, arrays


def directional_tangent(
    kernel: V55ProductionKernel, state: ProductionState, vacuum: Any
) -> tuple[dict[str, Any], dict[float, dict[str, np.ndarray]]]:
    """Measure H/M/C_Lambda/det tangents with a fresh CMC lapse per endpoint."""
    original = state_digest(state)
    rates = rhs_direction(kernel, state)
    r = np.asarray(state.grid.centers, dtype=np.float64)
    baseline = constraints(state, vacuum)
    rows = []
    raw_rows: dict[float, dict[str, np.ndarray]] = {}
    raw_budget_profiles: dict[float, dict[str, np.ndarray]] = {}
    for eps in DIRECTIONAL_EPSILONS:
        plus_state = shifted_state(kernel, state, rates, eps)
        minus_state = shifted_state(kernel, state, rates, -eps)
        plus = constraints(plus_state, vacuum)
        minus = constraints(minus_state, vacuum)
        dot = {
            name: finite(f"{name} tangent", (plus[name] - minus[name]) / (2.0 * eps))
            for name in baseline
        }
        raw_rows[eps] = dot
        # Use the shared exact-pin momentum observer for a signed decomposition
        # into geometric terms and all total-current sectors. This observer
        # independently refreshes CMC alpha at both directional endpoints.
        momentum_budget = momentum_tangent_for_epsilon(kernel, state, rates, eps)
        propagation_report, propagation_arrays = _constraint_propagation_candidate(
            state, rates, baseline, plus_state, minus_state, eps, dot["M"]
        )
        raw_budget_profiles[eps] = {
            **{
                name: finite(f"saved {name}", values).copy()
                for name, values in propagation_arrays.items()
            },
            "Mdot": finite("signed Mdot profile", momentum_budget["Mdot"]).copy(),
            "Mdot_geom": finite("signed geometric Mdot profile", momentum_budget["Mdot_geom"]).copy(),
            **{
                f"Mdot_geom_term_{name}": finite(f"signed geometric term {name}", values).copy()
                for name, values in momentum_budget["Mdot_geom_terms"].items()
            },
            **{
                f"Mdot_sector_{name}": finite(f"signed current sector {name}", values).copy()
                for name, values in momentum_budget["Mdot_sector"].items()
            },
        }
        rows.append({
            "epsilon": eps,
            "tangent_witnesses": {name: witness(arr, r) for name, arr in dot.items()},
            "signed_momentum_budget": compact_tangent_result(momentum_budget, r),
            "independent_constraint_propagation_candidate": propagation_report,
        })
    if state_digest(state) != original:
        raise RuntimeError("tangent observer mutated the frozen state")
    coarse, fine = raw_rows[DIRECTIONAL_EPSILONS[0]], raw_rows[DIRECTIONAL_EPSILONS[1]]
    report = {
        "frozen_state_sha256": original,
        "baseline_constraints": {name: witness(arr, r) for name, arr in baseline.items()},
        "tangent_by_epsilon": rows,
        "epsilon_difference_coarse_minus_fine": {
            name: witness(coarse[name] - fine[name], r) for name in baseline
        },
        "state_unchanged": state_digest(state) == original,
        "tangent_scope": (
            "Measured directional tangent of the implemented semi-discrete RHS; each +/-epsilon endpoint rebuilds shifted geometry/scalars/matter and re-solves the production CMC lapse. "
            "The signed M budget reports geometric terms and all six total-current sectors, with raw signed profiles retained in the NPZ. "
            "This remains a measured tangent, not the independent normalized full constraint-propagation prediction."
        ),
    }
    return report, raw_budget_profiles

def save_state_profiles(
    profile_store: dict[str, np.ndarray],
    n: int,
    arm: str,
    state: ProductionState,
    B_metric: np.ndarray,
    tangent_report: dict[str, Any],
    raw_budget_profiles: dict[float, dict[str, np.ndarray]],
) -> None:
    prefix = f"N{n}_{arm}"
    profile_store[f"{prefix}_r"] = np.asarray(state.grid.centers, dtype=np.float64)
    profile_store[f"{prefix}_metric_B_X6"] = np.asarray(B_metric, dtype=np.float64)
    for name in ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B"):
        profile_store[f"{prefix}_geometry_{name}"] = np.asarray(getattr(state.geometry, name), dtype=np.float64)
    for name in SCALAR_NAMES:
        profile_store[f"{prefix}_scalar_{name}"] = np.asarray(getattr(state.scalars, name), dtype=np.float64)
    for species_name in ("dark_matter", "baryons", "radiation"):
        species = getattr(state.matter, species_name)
        for name in ("rest", "energy_t", "momentum_r"):
            profile_store[f"{prefix}_matter_{species_name}_{name}"] = np.asarray(getattr(species, name), dtype=np.float64)
    for item in tangent_report["tangent_by_epsilon"]:
        eps_tag = format(float(item["epsilon"]), ".0e").replace("+", "")
        eps = float(item["epsilon"])
        for name, values in raw_budget_profiles[eps].items():
            profile_store[f"{prefix}_tangent_epsilon_{eps_tag}_{name}"] = np.asarray(values, dtype=np.float64)
        for key in ("H", "M", "C_Lambda", "determinant"):
            profile_store[f"{prefix}_tangent_epsilon_{eps_tag}_{key}_cell_0_to_4"] = np.asarray(
                [cell["signed"] for cell in item["tangent_witnesses"][key]["cell_0_to_4"]],
                dtype=np.float64,
            )


def run(args: argparse.Namespace) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    pin_info = provenance()
    if pin_info["production_source_baseline_anchor"] != SOURCE_BASELINE:
        raise RuntimeError("production source anchor differs from the declared baseline")
    if pin_info["vendor_commit"] != EXPECTED_VENDOR:
        raise RuntimeError("vendor revision differs from the declared pin")
    if len(pin_info["production_source_blobs"]) != 11:
        raise RuntimeError("expected exactly eleven pinned production engine blobs")

    kernel = V55ProductionKernel()
    _, vacuum, _ = adapter.vendor_modules()
    records = []
    profiles: dict[str, np.ndarray] = {}
    for n in args.resolutions:
        template = kernel.initialize(
            resolution=int(n),
            r_max=args.r_max,
            amplitude=AMPLITUDE,
            width=WIDTH,
            D_amplitude=D_AMPLITUDE,
            include_radiation=True,
        )
        template_digest = state_digest(template)
        r = np.asarray(template.grid.centers, dtype=np.float64)
        r0 = float(r[0])

        # Independently replay the standard six-pass B recurrence and require
        # its X^6 profile to match the production initializer.
        baseline_B_replay, baseline_iter = source_iteration(template, None, iterations=6)
        production_B = finite("production metric B", np.asarray(template.geometry.X) ** 6)
        baseline_error = float(np.max(np.abs(baseline_B_replay - production_B)))
        if baseline_error > 3.0e-12:
            raise RuntimeError(f"N={n}: replayed baseline B differs from production by {baseline_error:.6e}")

        reference = continuous_reference(r0, n_levels=CONTINUOUS_SUBSTEPS)
        envelope = analytic_half_cell_envelope(r0)
        moment = float(reference["reference_moment"])
        uncertainty = float(reference["estimated_numeric_uncertainty"])
        pad = 128.0 * np.finfo(np.float64).eps * max(1.0, abs(moment))
        low = float(envelope["continuous_integral_lower_bound"])
        high = float(envelope["continuous_integral_upper_bound"])
        if not low - pad <= moment <= high + pad:
            raise RuntimeError(f"N={n}: Richardson moment lies outside analytic enclosure")
        if uncertainty > 5.0e-13 * max(1.0, abs(moment)):
            raise RuntimeError(f"N={n}: continuous origin moment fails numerical-uncertainty gate")

        B_principled, principled_iter = source_iteration(
            template, first_half_cell_moment=moment, iterations=6
        )
        principal_state = rebuilt_state_from_B(kernel, template, B_principled)
        if state_digest(template) != template_digest:
            raise RuntimeError(f"N={n}: principled counterfactual mutated production input")

        B_hybrid, newton = solve_five_cell_hybrid(
            kernel, template, B_principled, vacuum
        )
        hybrid_state = rebuilt_state_from_B(kernel, template, B_hybrid)
        if state_digest(template) != template_digest:
            raise RuntimeError(f"N={n}: hybrid counterfactual mutated production input")

        all_states = {
            "baseline_production_initializer": (template, production_B),
            "principled_origin_integral": (principal_state, B_principled),
            "hybrid_five_cell_patch_newton": (hybrid_state, B_hybrid),
        }
        arms: dict[str, Any] = {}
        constraint_arrays = {arm: constraints(st, vacuum) for arm, (st, _) in all_states.items()}
        final_patch_max = float(np.max(np.abs(constraint_arrays["hybrid_five_cell_patch_newton"]["H"][:PATCH_CELLS])))
        if not newton["converged"] or final_patch_max > NEWTON_TOL:
            raise RuntimeError(f"N={n}: final rebuilt H patch {final_patch_max:.6e} exceeds Newton tolerance")
        if not np.array_equal(B_hybrid[PATCH_CELLS:], B_principled[PATCH_CELLS:]):
            raise RuntimeError(f"N={n}: hybrid changed B outside the five-cell patch")

        for arm, (st, metric_B) in all_states.items():
            start_digest = state_digest(st)
            tang, raw_budget_profiles = directional_tangent(kernel, st, vacuum)
            arms[arm] = {
                "state_sha256": start_digest,
                "r_cells_0_to_4": [float(v) for v in r[:PATCH_CELLS]],
                "metric_B_min": float(np.min(metric_B)),
                "metric_B_max": float(np.max(metric_B)),
                "B_positive": bool(np.all(metric_B > 0.0)),
                "metric_identity_max_abs_X6_minus_B": float(np.max(np.abs(np.asarray(st.geometry.X) ** 6 - metric_B))),
                "determinant_identity_max_abs": float(np.max(np.abs(np.asarray(st.geometry.a) * np.asarray(st.geometry.b) ** 2 - 1.0))),
                "lapse_min": float(np.min(st.geometry.alpha)),
                "lapse_max": float(np.max(st.geometry.alpha)),
                "constraints": {name: witness(arr, r) for name, arr in constraint_arrays[arm].items()},
                "measured_rhs_tangent": tang,
            }
            save_state_profiles(profiles, int(n), arm, st, metric_B, tang, raw_budget_profiles)
            if state_digest(st) != start_digest:
                raise RuntimeError(f"N={n}/{arm}: observer mutated the input state")

        delta_rows = {}
        for name in ("H", "M", "C_Lambda", "determinant"):
            delta_rows[name] = {
                "principled_minus_baseline": witness(
                    constraint_arrays["principled_origin_integral"][name] - constraint_arrays["baseline_production_initializer"][name], r
                ),
                "hybrid_minus_principled": witness(
                    constraint_arrays["hybrid_five_cell_patch_newton"][name] - constraint_arrays["principled_origin_integral"][name], r
                ),
                "hybrid_minus_baseline": witness(
                    constraint_arrays["hybrid_five_cell_patch_newton"][name] - constraint_arrays["baseline_production_initializer"][name], r
                ),
            }
        records.append({
            "resolution": int(n),
            "r_max": float(args.r_max),
            "dr": float(template.grid.dr),
            "r_cells_0_to_4": [float(v) for v in r[:PATCH_CELLS]],
            "template_state_sha256": template_digest,
            "baseline_reconstruction": {
                "source_iteration": baseline_iter,
                "max_abs_difference_replayed_B_vs_production_X6": baseline_error,
                "tolerance": 3.0e-12,
            },
            "continuous_origin_integral": {
                "radius": r0,
                "moment": moment,
                "equivalent_B0_from_integral": moment / r0,
                "estimated_numeric_uncertainty": uncertainty,
                "change_between_final_extrapolants": float(reference["change_between_final_extrapolants"]),
                "analytic_lower_bound": low,
                "analytic_upper_bound": high,
                "inside_analytic_enclosure": True,
                "method": reference["method"],
                "limitation": reference["limit"],
            },
            "principled_origin_iteration": principled_iter,
            "hybrid_newton": {
                **newton,
                "final_fully_rebuilt_H_cells_0_to_4_max_abs": final_patch_max,
                "outside_patch_B_unchanged": bool(np.array_equal(B_hybrid[PATCH_CELLS:], B_principled[PATCH_CELLS:])),
                "max_abs_B_change_inside_patch": float(np.max(np.abs(B_hybrid[:PATCH_CELLS] - B_principled[:PATCH_CELLS]))),
                "H_residual_outside_patch": exterior_witness(
                    constraint_arrays["hybrid_five_cell_patch_newton"]["H"][PATCH_CELLS:],
                    r[PATCH_CELLS:],
                ) if len(r) > PATCH_CELLS else None,
            },
            "arms": arms,
            "signed_constraint_deltas": delta_rows,
            "input_mutation_check": {
                "baseline_template_unchanged": state_digest(template) == template_digest,
                "all_arm_tangent_inputs_unchanged": all(arm["measured_rhs_tangent"]["state_unchanged"] for arm in arms.values()),
            },
            "independent_full_propagation_prediction": {
                "status": "NOT_COMPUTED",
                "acceptance_oracle_satisfied": False,
                "required_terms": [
                    "repository-normalized Hamiltonian-constraint gradient",
                    "reference-metric radial derivative/coupling of the connection constraint",
                    "total stress-energy conservation/exchange residual, including scalar-DM exchange",
                ],
                "reason": "This report measures the implemented semi-discrete tangent but does not claim to implement the independently normalized full propagation law.",
            },
            "scientific_admission": False,
        })
    report = {
        "schema": SCHEMA,
        "status": "completed_diagnostic_only",
        "scientific_admission": False,
        "production_or_vendor_files_modified": False,
        "evolution_steps_advanced": 0,
        "source_provenance": pin_info,
        "configuration": {
            "resolutions": [int(n) for n in args.resolutions],
            "r_max": float(args.r_max),
            "amplitude": AMPLITUDE,
            "width": WIDTH,
            "D_amplitude": D_AMPLITUDE,
            "radiation_enabled": True,
            "continuous_substeps": list(CONTINUOUS_SUBSTEPS),
            "hybrid_patch_cells": list(range(PATCH_CELLS)),
            "newton": {
                "max_iterations": NEWTON_MAX_ITER,
                "tolerance": NEWTON_TOL,
                "max_condition_number": NEWTON_MAX_CONDITION,
                "relative_finite_difference_epsilon": NEWTON_RELATIVE_FD_EPS,
                "backtracking": "halve step up to 24 times; accept only if max_abs_H_patch strictly decreases",
            },
            "directional_epsilons": list(DIRECTIONAL_EPSILONS),
            "fixed_radii": list(FIXED_RADII),
            "arm_order": [
                "baseline_production_initializer",
                "principled_origin_integral",
                "hybrid_five_cell_patch_newton",
            ],
        },
        "experiment_scope": {
            "hybrid_constructed_reproducibly": True,
            "hybrid_is_exact_reproduction_of_prior_supplied_state": False,
            "state_arrays_saved": True,
            "measured_signed_rhs_tangents_reported": True,
            "independent_full_propagation_prediction_computed": False,
            "reason": "The former supplied hybrid result lacked an exact method/configuration/state artifact. This creates a new pinned hybrid state family; it does not claim bitwise or numerical identity with that earlier state.",
        },
        "interpretation_limits": [
            "The principled arm changes only the origin half-cell moment to the converged continuous-source reference; inter-centre trapezoids remain the production ones.",
            "The hybrid arm solves the exact vendor H residual in cells 0..4 by varying B[0:5] only. B[5:] remains identical to the principled arm; dependent metric, archive matter and CMC lapse are rebuilt at each residual evaluation.",
            "The copied initializer's six-step source/B iteration is a fixed finite iteration, not a proof of a converged coupled density-metric fixed point.",
            "A small initial H patch and an exact accounting closure do not imply constraint propagation. The independent full propagation RHS remains uncomputed and is the acceptance gate.",
            "No evolution steps, production edits, equation substitutions, field projections or scientific-admission claims are made.",
        ],
        "resolutions": records,
    }
    return report, profiles


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resolutions", nargs="+", type=int, default=list(DEFAULT_RESOLUTIONS))
    parser.add_argument("--r-max", type=float, default=DEFAULT_R_MAX)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    if not args.resolutions or any(n < PATCH_CELLS + 2 for n in args.resolutions):
        raise ValueError("each resolution must be at least PATCH_CELLS + 2")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    report_path = args.output_dir / "report.json"
    profile_path = args.output_dir / "state_profiles.npz"
    try:
        report, profile_arrays = run(args)
        np.savez_compressed(profile_path, **profile_arrays)
        report["artifacts"] = {
            "report_json": "report.json",
            "state_profiles_npz": "state_profiles.npz",
            "profile_key_count": int(len(profile_arrays)),
        }
        report["artifact_hashes"] = {
            "state_profiles_npz_sha256": sha256_file(profile_path)
        }
        report_path.write_text(
            json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        (args.output_dir / "SHA256SUMS.txt").write_text(
            f"{sha256_file(report_path)}  report.json\n"
            f"{sha256_file(profile_path)}  state_profiles.npz\n",
            encoding="utf-8",
        )
        code = 0
    except Exception as exc:
        report = {
            "schema": SCHEMA,
            "status": "failed_closed",
            "scientific_admission": False,
            "production_or_vendor_files_modified": False,
            "error_type": type(exc).__name__,
            "error": str(exc),
            "traceback": traceback.format_exc(),
        }
        report_path.write_text(
            json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        code = 1

    print(json.dumps({
        "schema": report.get("schema"),
        "status": report.get("status"),
        "source_commit": report.get("source_provenance", {}).get("repository_commit"),
        "vendor_commit": report.get("source_provenance", {}).get("vendor_commit"),
        "resolutions": [
            {
                "N": row.get("resolution"),
                "newton_converged": row.get("hybrid_newton", {}).get("converged"),
                "final_H_patch_max": row.get("hybrid_newton", {}).get("final_fully_rebuilt_H_cells_0_to_4_max_abs"),
                "Mdot_cell0": {
                    arm: armrow.get("measured_rhs_tangent", {}).get("tangent_by_epsilon", [{}])[0]
                        .get("tangent_witnesses", {}).get("M", {}).get("cell_0_to_4", [{}])[0].get("signed")
                    for arm, armrow in row.get("arms", {}).items()
                } if "arms" in row else {},
            } for row in report.get("resolutions", [])
        ],
        "experiment_scope": report.get("experiment_scope"),
        "error_type": report.get("error_type"),
        "error": report.get("error"),
        "scientific_admission": False,
    }, indent=2, sort_keys=True), flush=True)
    return code


if __name__ == "__main__":
    raise SystemExit(main())
