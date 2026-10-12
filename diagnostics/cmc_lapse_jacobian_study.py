"""Finite-difference CMC lapse-map Jacobian and input-slice residual audit.

Diagnostic only: freezes all state except alpha; does not edit production
physics. It measures G(alpha)=solve_cmc_lapse(alpha), the local directional
Jacobian of G, and the discrete CMC residual on the incoming slice.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import engine.cmc_gauge as cmc
from engine.production_kernel import V55ProductionKernel
from diagnostics.cmc_lapse_residual_audit import exact_discrete_residual

OUT = ROOT / "runs" / "cmc_lapse_jacobian_study"
EPSILONS = (1.0e-4, 1.0e-5, 1.0e-6)


def rms(x):
    x = np.asarray(x, dtype=float)
    return float(np.sqrt(np.mean(x * x)))


def max_abs(x):
    return float(np.max(np.abs(np.asarray(x, dtype=float))))


def cosine_similarity(left, right):
    left = np.asarray(left, dtype=float)
    right = np.asarray(right, dtype=float)
    denom = float(np.linalg.norm(left) * np.linalg.norm(right))
    return None if denom == 0.0 else float(np.dot(left, right) / denom)


def solve_map(state, frozen_geometry, alpha):
    geometry = frozen_geometry.copy()
    geometry.alpha = np.asarray(alpha, dtype=float).copy()
    output, target_used = cmc.solve_cmc_lapse(
        state.grid, geometry, state.scalars, state.matter
    )
    return np.asarray(output, dtype=float), float(target_used)


def build_modes(grid):
    """Build smooth, independent radial perturbations respecting both ends."""
    r = np.asarray(grid.centers, dtype=float)
    rmax = float(r[-1])
    x = r / rmax
    raw = [
        ("central_gaussian", np.exp(-(r / 2.5) ** 2)),
        ("broad_gaussian", np.exp(-(r / 10.0) ** 2)),
        ("off_center_shell", np.exp(-0.5 * ((r - 12.0) / 4.0) ** 2)),
        ("signed_radial_mode", np.cos(2.0 * np.pi * x) * (1.0 - x) ** 2),
    ]
    modes = []
    for name, values in raw:
        values = np.asarray(values, dtype=float)
        values[0] = values[1]  # satisfy the discrete regular-center row
        values[-1] = 0.0       # preserve the outer lapse normalization
        scale = rms(values)
        if scale == 0.0:
            raise RuntimeError(f"degenerate perturbation mode: {name}")
        modes.append((name, values / scale))
    return modes


def profile_summary(values):
    values = np.asarray(values, dtype=float)
    return {
        "rms": rms(values),
        "max_abs": max_abs(values),
        "center": float(values[0]),
        "outer": float(values[-1]),
    }


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
    frozen = state.geometry.copy()
    alpha0 = np.asarray(frozen.alpha, dtype=float).copy()

    # F(alpha0) = G(alpha0) - alpha0, where G is the actual production
    # solve/target pair. All non-lapse state remains frozen.
    alpha_base, kdot_used = solve_map(state, frozen, alpha0)
    base_shift = alpha_base - alpha0
    geometry_base = frozen.copy()
    geometry_base.alpha = alpha_base.copy()
    kdot_recomputed = float(
        cmc.target_kdot(state.grid, geometry_base, state.scalars, state.matter)
    )

    residual_input, kdot_input, input_scale = exact_discrete_residual(
        state.grid, frozen, state.scalars, state.matter, alpha0
    )
    # First, evaluate the mapped output against the same frozen target. This
    # isolates linear-solver closure from target self-consistency.
    residual_base_fixed_target, _, fixed_scale = exact_discrete_residual(
        state.grid, geometry_base, state.scalars, state.matter, alpha_base,
        kdot=kdot_input,
    )
    residual_base_recomputed, _, recomputed_scale = exact_discrete_residual(
        state.grid, geometry_base, state.scalars, state.matter, alpha_base
    )

    # Candidate profiles are all RMS-normalized. Central differences estimate
    # J_G v = dG(alpha0)[v] over a short amplitude ladder to expose truncation
    # and roundoff sensitivity.
    modes = build_modes(state.grid)
    directional = []
    derivative_vectors = []
    eps_details = {}
    finite = True
    for name, mode in modes:
        by_epsilon = []
        for eps in EPSILONS:
            plus_alpha = alpha0 + eps * mode
            minus_alpha = alpha0 - eps * mode
            plus, plus_target = solve_map(state, frozen, plus_alpha)
            minus, minus_target = solve_map(state, frozen, minus_alpha)
            jv = (plus - minus) / (2.0 * eps)
            dkdot = (plus_target - minus_target) / (2.0 * eps)
            item = {
                "epsilon": eps,
                "map_directional_derivative": profile_summary(jv),
                "directional_response_rms_ratio": rms(jv) / max(rms(mode), 1.0e-300),
                "identity_response_error_rms_ratio": rms(jv - mode) / max(rms(mode), 1.0e-300),
                "rayleigh_quotient": float(np.dot(mode, jv) / np.dot(mode, mode)),
                "response_cosine_with_input_mode": cosine_similarity(mode, jv),
                "dkdot_d_direction": float(dkdot),
                "plus_target": plus_target,
                "minus_target": minus_target,
                "all_finite": bool(
                    np.all(np.isfinite(plus_alpha))
                    and np.all(np.isfinite(minus_alpha))
                    and np.all(np.isfinite(plus))
                    and np.all(np.isfinite(minus))
                    and np.isfinite(dkdot)
                ),
            }
            by_epsilon.append(item)
            finite = finite and item["all_finite"]
            if eps == EPSILONS[1]:
                derivative_vectors.append(jv)
        directional.append({"mode": name, "input_mode": profile_summary(mode), "by_epsilon": by_epsilon})
        eps_details[name] = by_epsilon

    # Estimate the Jacobian on an orthonormal radial subspace spanning the
    # candidate modes. The full output response is also checked for components
    # outside that subspace; this is not a claim about unsampled directions.
    raw_matrix = np.column_stack([m for _, m in modes])
    q, _ = np.linalg.qr(raw_matrix / np.sqrt(state.grid.n))
    basis = q * np.sqrt(state.grid.n)  # columns have RMS norm 1
    reduced_jacobian_columns = []
    off_subspace = []
    eps_matrix = EPSILONS[1]
    for j in range(basis.shape[1]):
        v = basis[:, j]
        gp, _ = solve_map(state, frozen, alpha0 + eps_matrix * v)
        gm, _ = solve_map(state, frozen, alpha0 - eps_matrix * v)
        jv = (gp - gm) / (2.0 * eps_matrix)
        coeff = (basis.T @ jv) / state.grid.n
        reduced_jacobian_columns.append(coeff)
        off_subspace.append(rms(jv - basis @ coeff))
        finite = finite and bool(np.all(np.isfinite(jv)) and np.all(np.isfinite(coeff)))

    reduced_jacobian = np.column_stack(reduced_jacobian_columns)
    residual_jacobian = reduced_jacobian - np.eye(reduced_jacobian.shape[0])
    jacobian_eigenvalues = np.linalg.eigvals(reduced_jacobian)
    defect_eigenvalues = np.linalg.eigvals(residual_jacobian)
    singular_values = np.linalg.svd(residual_jacobian, compute_uv=False)

    # The lapse shift is an inverse-elliptic response to the input CMC
    # residual, not the residual itself (they have different units/operators).
    # Subtracting same-target residuals checks L(G(alpha0)-alpha0) = -R(alpha0).
    correction_closure = (
        residual_base_fixed_target - residual_input + residual_input
    )

    payload = {
        "test": "CMC lapse-map finite-difference Jacobian and input-slice residual",
        "scope": (
            "V55ProductionKernel initial slice; N=64, r_max=32, localized S "
            "amplitude=1e-3, D amplitude=0, radiation off; alpha alone is "
            "perturbed while all other state fields remain frozen"
        ),
        "production_equations_changed": False,
        "fixed_point_map": "G(alpha) = solve_cmc_lapse(grid, geometry(alpha), scalars, matter); F(alpha) = G(alpha) - alpha",
        "base_state": {
            "input_alpha": profile_summary(alpha0),
            "map_offset_F_rms": rms(base_shift),
            "map_offset_F_max_abs": max_abs(base_shift),
            "map_offset_F_profile": profile_summary(base_shift),
            "kdot_used_by_map": float(kdot_used),
            "kdot_recomputed_on_mapped_output": kdot_recomputed,
            "kdot_output_minus_input": float(kdot_recomputed - kdot_input),
            "input_target_kdot": float(kdot_input),
            "input_slice_discrete_cmc_residual": profile_summary(residual_input),
            "input_slice_cmc_residual_relative_max": max_abs(residual_input) / input_scale,
            "mapped_output_residual_with_input_target": profile_summary(residual_base_fixed_target),
            "mapped_output_fixed_target_residual_relative_max": max_abs(residual_base_fixed_target) / fixed_scale,
            "mapped_output_residual_with_recomputed_target": profile_summary(residual_base_recomputed),
            "mapped_output_recomputed_target_residual_relative_max": max_abs(residual_base_recomputed) / recomputed_scale,
            "shift_residual_cosine_descriptive_only": cosine_similarity(base_shift, residual_input),
            "shift_residual_note": (
                "Do not equate shift and residual directly: the shift is the "
                "elliptic inverse-operator response to the residual. The same-target "
                "mapped solve is separately checked for residual closure."
            ),
            "same_target_residual_closure_max_abs": max_abs(correction_closure),
        },
        "finite_difference_epsilons": list(EPSILONS),
        "directional_modes": directional,
        "reduced_jacobian": {
            "basis": "orthonormal span of four smooth center-regular, outer-zero radial profiles",
            "finite_difference_epsilon": eps_matrix,
            "J_G_matrix": reduced_jacobian.tolist(),
            "J_G_eigenvalues": [
                {"real": float(np.real(z)), "imag": float(np.imag(z))}
                for z in jacobian_eigenvalues
            ],
            "J_F_equals_J_G_minus_I_matrix": residual_jacobian.tolist(),
            "J_F_eigenvalues": [
                {"real": float(np.real(z)), "imag": float(np.imag(z))}
                for z in defect_eigenvalues
            ],
            "J_F_singular_values": [float(v) for v in singular_values],
            "min_J_F_singular_value": float(np.min(singular_values)),
            "off_subspace_response_rms_by_basis_direction": [float(v) for v in off_subspace],
            "max_off_subspace_response_rms": float(max(off_subspace)),
            "interpretation_limit": (
                "Finite differences characterize only the local response on the "
                "sampled lapse slice and four-dimensional radial subspace. They "
                "cannot prove global existence/nonexistence of a self-consistent "
                "lapse or rule out behavior in unsampled directions."
            ),
        },
        "checks": {
            "all_finite": bool(finite),
            "fixed_target_linear_solve_residual_small": bool(
                max_abs(residual_base_fixed_target) / fixed_scale < 1.0e-10
            ),
            "input_and_output_diagnostic_residuals_finite": bool(
                np.all(np.isfinite(residual_input))
                and np.all(np.isfinite(residual_base_fixed_target))
                and np.all(np.isfinite(residual_base_recomputed))
            ),
        },
        "pass": bool(
            finite
            and max_abs(residual_base_fixed_target) / fixed_scale < 1.0e-10
            and np.all(np.isfinite(residual_input))
            and np.all(np.isfinite(residual_base_recomputed))
        ),
        "interpretation_limit": (
            "Weak-field initial slice only. A nonsingular sampled J_F is local "
            "evidence that root solving is not degenerate in the tested subspace; "
            "it does not guarantee a global root or explain strong-D proper time."
        ),
    }

    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / "cmc_lapse_jacobian_study.json"
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload, indent=2, sort_keys=True), flush=True)
    print(f"RESULT_JSON={path}", flush=True)
    if not payload["pass"]:
        raise SystemExit("FAIL: non-finite result or fixed-target solve residual mismatch")
    print("PASS: lapse Jacobian and input-slice residual characterized; diagnostic only")


if __name__ == "__main__":
    main()
