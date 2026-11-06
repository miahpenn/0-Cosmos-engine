"""Separate affine-Jacobian construction error from float64 root-solve error.

This is a new Issue #4 diagnostic. It sweeps the finite-difference step used
by the historical affine-root construction and compares float64 perturbed-root
differences with a 70-digit reference. No production code is modified.
"""
from __future__ import annotations

import json

import mpmath as mp
import numpy as np

from test_issue4_split_derivative_clean_rerun import _case, _perturbed_mp_state
from test_issue4_cmc_lapse_tangent import _d1, _d2, _fixed_point_map_output
from issue4_split_derivative_mp import as_mp_state, exact_fixed_point
from tools.issue4_cmc_lapse_tangent import solve_cmc_lapse_fixed_point


def _float_root(args):
    return solve_cmc_lapse_fixed_point(**args, d1=_d1, d2=_d2).alpha_in


def _affine_fd_root(args, step):
    n = len(args["r"])
    reference = np.ones(n)
    fref = _fixed_point_map_output(args, reference).alpha
    jac = np.empty((n, n))
    eye = np.eye(n)
    for j in range(n):
        plus = _fixed_point_map_output(args, reference + step * eye[j]).alpha
        minus = _fixed_point_map_output(args, reference - step * eye[j]).alpha
        jac[:, j] = (plus - minus) / (2.0 * step)
    system = np.eye(n) - jac
    root = np.linalg.solve(system, fref - jac @ reference)
    singular_values = np.linalg.svd(system, compute_uv=False)
    return root, jac, float(singular_values[-1]), float(np.linalg.cond(system))


def _affine_basis_root(args):
    """Recover the affine Jacobian from a finite-amplitude exact secant."""
    n = len(args["r"])
    anchor = np.ones(n)
    amplitude = 1.0e-2
    fanchor = _fixed_point_map_output(args, anchor).alpha
    jac = np.empty((n, n))
    for j in range(n):
        basis = anchor.copy()
        basis[j] += amplitude
        f_basis = _fixed_point_map_output(args, basis).alpha
        # The map is affine in incoming lapse, so this secant has no
        # truncation error; the finite amplitude keeps every input physical.
        jac[:, j] = (f_basis - fanchor) / amplitude
    system = np.eye(n) - jac
    intercept = fanchor - jac @ anchor
    root = np.linalg.solve(system, intercept)
    singular_values = np.linalg.svd(system, compute_uv=False)
    return root, jac, intercept, float(singular_values[-1]), float(np.linalg.cond(system))


def _perturbed_float_args(args, rates, epsilon):
    shifted = dict(args)
    for key in ("a", "b", "X", "Aa", "K", "rho", "pr", "pt"):
        shifted[key] = args[key] + epsilon * rates[key + "_dot"]
    return shifted


def _norm(x):
    return float(np.linalg.norm(x))


def test_report_jacobian_and_root_error_separation():
    _, _, args, fixed, rates, _, implicit = _case()
    mp_state = as_mp_state(args)
    mp_root_mp, mp_jac_mp = exact_fixed_point(mp_state, dps=70)
    mp_root = np.asarray([float(x) for x in mp_root_mp])
    mp_jac = np.asarray([[float(mp_jac_mp[i, j]) for j in range(len(mp_root))]
                         for i in range(len(mp_root))])

    basis_root, basis_jac, basis_intercept, basis_sigma_min, basis_condition = _affine_basis_root(args)
    basis_result = {
        "root_error_l2_vs_mp": _norm(basis_root - mp_root),
        "root_error_max_abs_vs_mp": float(np.max(np.abs(basis_root - mp_root))),
        "jacobian_error_l2_vs_mp": _norm(basis_jac - mp_jac),
        "system_sigma_min": basis_sigma_min,
        "system_condition_number": basis_condition,
    }

    # Solve the float64 affine system again in 70-digit arithmetic using
    # exactly the same binary64 coefficients. This separates coefficient
    # construction error from forward error in the float64 linear solve.
    basis_system = np.eye(len(mp_root)) - basis_jac
    with mp.workdps(70):
        mp_system64 = mp.matrix([
            [mp.mpf(float(basis_system[i, j])) for j in range(len(mp_root))]
            for i in range(len(mp_root))
        ])
        mp_intercept64 = mp.matrix([
            mp.mpf(float(value)) for value in basis_intercept
        ])
        mp_root_from_float64_coefficients = mp.lu_solve(
            mp_system64, mp_intercept64)
        mp_root_from_float64_coefficients_float = np.asarray([
            float(value) for value in mp_root_from_float64_coefficients
        ])
        mp_root_error_of_float64_solution = mp.matrix([
            mp.mpf(float(value)) for value in basis_root
        ]) - mp_root_from_float64_coefficients
        mp_system_residual = (
            mp_system64 * mp.matrix([mp.mpf(float(value)) for value in basis_root])
            - mp_intercept64
        )
        mp_intercept_reference = mp.matrix(mp_root_mp) - mp_jac_mp * mp.matrix(mp_root_mp)
        mp_system_reference = mp.eye(len(mp_root)) - mp_jac_mp
        solve_error_decomposition = {
            "float64_basis_root_error_l2_vs_mp": _norm(basis_root - mp_root),
            "mp_root_from_float64_coefficients_error_l2_vs_mp": _norm(
                mp_root_from_float64_coefficients_float - mp_root),
            "float64_solve_roundoff_error_l2_vs_mp_solve_same_coefficients": float(
                mp.sqrt(sum(value**2 for value in mp_root_error_of_float64_solution))),
            "float64_intercept_error_l2_vs_mp": _norm(
                basis_intercept - np.asarray([float(x) for x in mp_intercept_reference])),
            "float64_system_matrix_error_l2_vs_mp": _norm(
                basis_system - np.asarray([[float(mp_system_reference[i, j])
                                             for j in range(len(mp_root))]
                                            for i in range(len(mp_root))])),
            "float64_solve_residual_l2_high_precision": float(
                mp.sqrt(sum(value**2 for value in mp_system_residual))),
            "float64_solve_residual_linf_high_precision": float(
                max(abs(value) for value in mp_system_residual)),
        }

    affine_sweep = []
    for step in (1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8):
        root, jac, sigma_min, condition = _affine_fd_root(args, step)
        affine_sweep.append({
            "jacobian_step": step,
            "jacobian_error_l2": _norm(jac - mp_jac),
            "root_error_l2_vs_mp": _norm(root - mp_root),
            "root_error_max_abs_vs_mp": float(np.max(np.abs(root - mp_root))),
            "system_sigma_min": sigma_min,
            "system_condition_number": condition,
        })

    root_fd_sweep = []
    for epsilon in (2e-3, 1e-3, 5e-4, 2e-4, 1e-4, 5e-5, 2e-5, 1e-5, 5e-6):
        plus_float = _float_root(_perturbed_float_args(args, rates, epsilon))
        minus_float = _float_root(_perturbed_float_args(args, rates, -epsilon))
        float_fd = (plus_float - minus_float) / (2.0 * epsilon)

        plus_mp, _ = exact_fixed_point(
            _perturbed_mp_state(args, rates, epsilon), dps=70)
        minus_mp, _ = exact_fixed_point(
            _perturbed_mp_state(args, rates, -epsilon), dps=70)
        mp_fd = np.asarray([
            float((plus_mp[i] - minus_mp[i]) /
                  (2 * mp.mpf(str(epsilon))))
            for i in range(len(plus_mp))
        ])
        plus_mp_float = np.asarray([float(value) for value in plus_mp])
        minus_mp_float = np.asarray([float(value) for value in minus_mp])
        root_fd_sweep.append({
            "float64_plus_root_error_l2_vs_mp": _norm(plus_float - plus_mp_float),
            "float64_minus_root_error_l2_vs_mp": _norm(minus_float - minus_mp_float),
            "root_error_difference_amplified_by_1_over_2epsilon": _norm(
                ((plus_float - plus_mp_float) - (minus_float - minus_mp_float))
                / (2.0 * epsilon)),
            "epsilon": epsilon,
            "float64_fd_vs_mp_exact_fd_l2": _norm(float_fd - mp_fd),
            "float64_fd_vs_implicit_tangent_l2": _norm(float_fd - implicit.alpha_dot),
            "mp_exact_fd_vs_implicit_tangent_l2": _norm(mp_fd - implicit.alpha_dot),
            "float64_plus_root_residual": float(
                solve_cmc_lapse_fixed_point(
                    **_perturbed_float_args(args, rates, epsilon),
                    d1=_d1, d2=_d2).fixed_point_residual),
            "float64_minus_root_residual": float(
                solve_cmc_lapse_fixed_point(
                    **_perturbed_float_args(args, rates, -epsilon),
                    d1=_d1, d2=_d2).fixed_point_residual),
        })

    result = {
        "schema": "issue4.jacobian_root_error_separation.v1",
        "reference": "70-decimal-digit mpmath affine fixed point and Jacobian",
        "baseline_float64_root_error_l2": _norm(fixed.alpha_in - mp_root),
        "affine_basis_jacobian_exact_secant": basis_result,
        "affine_root_error_decomposition": solve_error_decomposition,
        "affine_jacobian_finite_difference_sweep": affine_sweep,
        "perturbed_float64_root_finite_difference_sweep": root_fd_sweep,
        "interpretation_limits": [
            "The two sweeps isolate numerical contributions; neither changes a solver tolerance.",
            "The affine-root sweep reproduces the old central-difference construction with an explicit step sweep.",
            "The perturbed-root sweep compares float64 solver roots against exact high-precision roots at matching state perturbations.",
            "No cause is assigned until the measured curves support it."
        ]
    }
    print("ISSUE4_JACOBIAN_ROOT_SEPARATION=" + json.dumps(result, sort_keys=True))
    assert all(np.isfinite(row["root_error_l2_vs_mp"])
               for row in affine_sweep)
    assert all(np.isfinite(row["float64_fd_vs_mp_exact_fd_l2"])
               for row in root_fd_sweep)
