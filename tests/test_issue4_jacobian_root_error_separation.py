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
    mp_root, mp_jac = exact_fixed_point(mp_state, dps=70)
    mp_root = np.asarray([float(x) for x in mp_root])
    mp_jac = np.asarray([[float(mp_jac[i, j]) for j in range(len(mp_root))]
                         for i in range(len(mp_root))])

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
        root_fd_sweep.append({
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
