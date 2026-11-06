"""Independent numerical sensitivity study for the Issue #4 affine-root fixture.

This protocol is new and explicit. It does not claim to reproduce Lumen's
unavailable perturbation schedule or decision rule. Historical reported
boundaries are comparison points only.
"""
from __future__ import annotations

import json

import mpmath as mp
import numpy as np

from test_issue4_split_derivative_clean_rerun import _case
from test_issue4_cmc_lapse_tangent import _stage_state, _d1, _d2
from issue4_split_derivative_mp import as_mp_state, exact_fixed_point
from tools.issue4_cmc_lapse_tangent import solve_cmc_lapse_fixed_point


def _mp_perturb(state, field, amplitude, scale, shape):
    out = {key: (list(value) if isinstance(value, list) else value)
           for key, value in state.items()}
    step = mp.mpf(amplitude)
    if field == "__uniform__":
        for key in ("a", "b", "X", "Aa", "K", "rho", "pr", "pt"):
            field_scale = max(abs(x) for x in state[key])
            out[key] = [
                state[key][i] + step * field_scale * shape[i]
                for i in range(len(shape))
            ]
    else:
        out[field] = list(state[field])
        out[field][len(shape) // 2] += step
    return out


def _float_perturb(args, field, amplitude, shape):
    out = dict(args)
    if field == "__uniform__":
        for key in ("a", "b", "X", "Aa", "K", "rho", "pr", "pt"):
            values = np.asarray(args[key])
            scale = float(np.max(np.abs(values)))
            out[key] = values + float(amplitude) * scale * shape
    else:
        values = np.asarray(args[field]).copy()
        values[len(shape) // 2] += float(amplitude)
        out[field] = values
    return out


def _root_float(args):
    return solve_cmc_lapse_fixed_point(**args, d1=_d1, d2=_d2).alpha_in


def _root_mp(state):
    root, _ = exact_fixed_point(state, dps=70)
    return np.asarray([float(value) for value in root])


def _measure(field, amplitudes, args, mp_base, float_base, shape, baseline_error):
    rows = []
    for amplitude in amplitudes:
        mp_state = _mp_perturb(
            mp_base, field, amplitude, None,
            [mp.mpf(str(float(x))) for x in shape])
        float_args = _float_perturb(args, field, amplitude, shape)
        mp_delta = _root_mp(mp_state) - mp_base_root
        float_delta = _root_float(float_args) - float_base
        signal = float(np.linalg.norm(mp_delta))
        discrepancy = float(np.linalg.norm(float_delta - mp_delta))
        resolved = (
            signal > 10.0 * baseline_error
            and discrepancy <= 0.1 * signal
        )
        rows.append({
            "amplitude": float(amplitude),
            "root_shift_l2": signal,
            "float64_vs_mp_delta_l2": discrepancy,
            "resolved": bool(resolved),
        })
    resolved_amplitudes = [row["amplitude"] for row in rows if row["resolved"]]
    return {
        "perturbation": "uniform_relative_all_state_fields" if field == "__uniform__"
                        else f"single_cell_absolute_{field}[{len(shape)//2}]",
        "smallest_resolved_amplitude": min(resolved_amplitudes)
                                      if resolved_amplitudes else None,
        "rows": rows,
    }


def test_report_new_affine_root_sensitivity_protocol(case=None):
    global mp_base_root
    _, _, args, fixed, _, _, _ = _case()
    mp_base = as_mp_state(args)
    mp_root_values, mp_jacobian = exact_fixed_point(mp_base, dps=70)
    mp_base_root = np.asarray([float(value) for value in mp_root_values])
    float_base = np.asarray(fixed.alpha_in)
    baseline_error = float(np.linalg.norm(float_base - mp_base_root))
    n = len(float_base)
    system = np.eye(n) - np.asarray(
        [[float(mp_jacobian[i, j]) for j in range(n)] for i in range(n)])
    sigma_min = float(np.linalg.svd(system, compute_uv=False)[-1])
    conditioning_bound = np.sqrt(n) * 2.0e-12 / sigma_min
    shape = np.linspace(0.75, 1.25, n)

    uniform_amplitudes = ["1e-11", "1e-10", "1e-9", "1e-8",
                          "1e-7", "1e-6", "1e-5"]
    single_amplitudes = ["1e-13", "1e-12", "1e-11", "1e-10",
                         "1e-9", "1e-8", "1e-7"]

    results = {
        "schema": "issue4.affine_root_sensitivity_study.v1",
        "protocol": {
            "reference": "70-decimal-digit mpmath affine-root solution",
            "float64_solver": "solve_cmc_lapse_fixed_point",
            "detection_rule": (
                "root-shift L2 norm exceeds 10 times baseline float64-vs-mpmath "
                "root error, and float64 root-shift discrepancy is at most "
                "10 percent of the high-precision root shift"
            ),
            "uniform_perturbation": (
                "all eight state arrays a,b,X,Aa,K,rho,pr,pt; each receives "
                "the same positive smooth shape times its own max-absolute scale"
            ),
            "single_component_perturbation": (
                "one absolute perturbation to the midpoint cell of a; this is "
                "one explicit component choice, not a claim about Lumen's choice"
            ),
            "historical_comparison_only": {
                "uniform_normwise_relative_reported": 1e-7,
                "single_component_absolute_reported": 2e-10,
                "original_schedule_and_decision_rule_available": False,
            },
        },
        "baseline_float64_vs_mpmath_root_error_l2": baseline_error,
        "implicit_system_sigma_min": sigma_min,
        "root_error_conditioning_bound_from_2e_12_residual": conditioning_bound,
        "baseline_error_over_conditioning_bound": baseline_error / conditioning_bound,
        "uniform": _measure("__uniform__", uniform_amplitudes, args, mp_base,
                            float_base, shape, baseline_error),
        "single_component": _measure("a", single_amplitudes, args, mp_base,
                                     float_base, shape, baseline_error),
    }
    print("ISSUE4_SENSITIVITY_RESULT=" + json.dumps(results, sort_keys=True))
    assert baseline_error >= 0.0
    assert all(np.isfinite(row["root_shift_l2"])
               for section in (results["uniform"], results["single_component"])
               for row in section["rows"])
    assert any(row["resolved"] for section in
               (results["uniform"], results["single_component"])
               for row in section["rows"])
