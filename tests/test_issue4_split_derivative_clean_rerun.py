"""Clean, blocking Issue #4 split-derivative contracts.

All checks use the independent research implementation and synthetic fixture.
The original strict tests remain untouched as visible historical witnesses.
"""
from __future__ import annotations

import numpy as np
import pytest

from test_issue4_cmc_lapse_tangent import (
    _stage_state,
    _d1,
    _d2,
    _physical_six_sector_split,
    _fixed_point_map_block_tangent,
    _fixed_point_map_output,
)
from tools.issue4_cmc_lapse_tangent import (
    solve_cmc_lapse_fixed_point,
    differentiate_cmc_lapse_fixed_point,
)
from tests.issue4_split_derivative_mp import as_mp_state, exact_fixed_point


EPSILONS = (2.0e-5, 1.0e-5)
ROOT_RESIDUAL_GATE = 2.0e-12
TANGENT_RESIDUAL_GATE = 2.0e-10
ABSOLUTE_ERROR_GATE = 2.0e-8
NORMWISE_RELATIVE_GATE = 1.0e-7


def _case():
    base, stage_rates = _stage_state()
    sectors, directions, summed, _ = _physical_six_sector_split(
        r=base["r"], dr=base["dr"], a=base["a"], b=base["b"],
        X=base["X"], Aa=base["Aa"], K=base["K"],
        alpha=base["alpha_in"])
    args = {key: base[key] for key in
            ("r", "dr", "a", "b", "X", "Aa", "K")}
    args.update({key: summed["base"][key] for key in ("rho", "pr", "pt")})
    fixed = solve_cmc_lapse_fixed_point(**args, d1=_d1, d2=_d2)

    # Recompute the physical directions at the fixed-point incoming lapse.
    sectors, directions, summed, physical_rates = _physical_six_sector_split(
        r=base["r"], dr=base["dr"], a=base["a"], b=base["b"],
        X=base["X"], Aa=base["Aa"], K=base["K"],
        alpha=fixed.alpha_in)
    state_rates = {
        key + "_dot": physical_rates["geometry"][key]
        for key in ("a", "b", "X", "Aa", "K")
    }
    for key in ("rho", "pr", "pt"):
        state_rates[key + "_dot"] = summed["direction"][key + "_dot"]
    zeros = {key: np.zeros_like(base["r"]) for key in state_rates}
    implicit = differentiate_cmc_lapse_fixed_point(
        fixed_point=fixed, **args,
        a_dot=state_rates["a_dot"], b_dot=state_rates["b_dot"],
        X_dot=state_rates["X_dot"], Aa_dot=state_rates["Aa_dot"],
        K_dot=state_rates["K_dot"], rho_dot=state_rates["rho_dot"],
        pr_dot=state_rates["pr_dot"], pt_dot=state_rates["pt_dot"],
        d1=_d1, d2=_d2)
    return base, stage_rates, args, fixed, state_rates, zeros, implicit


@pytest.fixture(scope="module")
def case():
    return _case()


def _perturbed_mp_state(args, state_rates, epsilon):
    perturbed = {"r": args["r"], "dr": args["dr"]}
    for key in ("a", "b", "X", "Aa", "K", "rho", "pr", "pt"):
        perturbed[key] = args[key] + epsilon * state_rates[key + "_dot"]
    return as_mp_state(perturbed)


def _normwise_relative(error, reference):
    denominator = float(np.linalg.norm(reference))
    if denominator == 0.0:
        return float("inf")
    return float(np.linalg.norm(error) / denominator)


def _within_contract(candidate, reference):
    error = candidate - reference
    return (
        float(np.max(np.abs(error))) <= ABSOLUTE_ERROR_GATE
        and _normwise_relative(error, reference) <= NORMWISE_RELATIVE_GATE
    )


def test_fixed_point_root_residual_is_within_solver_gate(case):
    _, _, _, fixed, _, _, _ = case
    assert fixed.fixed_point_residual <= ROOT_RESIDUAL_GATE
    assert fixed.linearization_evaluations == len(fixed.alpha_in) + 1
    assert np.all(np.isfinite(fixed.alpha_in))
    assert np.min(fixed.alpha_in) > 0.0


def test_affine_root_matches_extended_precision_conditioning_bound(case):
    _, _, args, fixed, _, _, _ = case
    mp_state = as_mp_state(args)
    root, jacobian = exact_fixed_point(mp_state, dps=70)
    reference = np.asarray([float(value) for value in root])
    system = np.eye(len(root)) - np.asarray(
        [[float(jacobian[i, j]) for j in range(len(root))]
         for i in range(len(root))])
    sigma_min = float(np.linalg.svd(system, compute_uv=False)[-1])
    assert np.isfinite(sigma_min) and sigma_min > 0.0
    bound = np.sqrt(len(root)) * ROOT_RESIDUAL_GATE / sigma_min
    error = float(np.max(np.abs(fixed.alpha_in - reference)))
    assert error <= bound


def test_mpmath_port_matches_float64_tool_map(case):
    _, _, args, fixed, _, _, _ = case
    mp_state = as_mp_state(args)
    mapped = exact_fixed_point  # Keep the reference port import explicit below.
    from tests.issue4_split_derivative_mp import _solve_map
    mp_output = _solve_map(mp_state, [__import__("mpmath").mpf(float(x))
                                      for x in fixed.alpha_in], dps=70)
    mp_output = np.asarray([float(value) for value in mp_output])
    float_output = _fixed_point_map_output(args, fixed.alpha_in).alpha
    np.testing.assert_allclose(mp_output, float_output, rtol=2.0e-10, atol=2.0e-12)


def test_incoming_lapse_only_direction_matches_two_central_differences(case):
    base, stage_rates, args, fixed, _, zeros, _ = case
    alpha_direction = stage_rates["alpha_in_dot"]
    _, tangent = _fixed_point_map_block_tangent(
        args, fixed.alpha_in, zeros, alpha_direction)
    assert np.max(np.abs(tangent.differentiated_residual)) < TANGENT_RESIDUAL_GATE
    for epsilon in EPSILONS:
        plus = _fixed_point_map_output(
            args, fixed.alpha_in + epsilon * alpha_direction)
        minus = _fixed_point_map_output(
            args, fixed.alpha_in - epsilon * alpha_direction)
        finite_difference = (plus.alpha - minus.alpha) / (2.0 * epsilon)
        np.testing.assert_allclose(
            finite_difference, tangent.alpha_dot,
            rtol=4.0e-5, atol=2.0e-8)


def test_geometry_and_six_sector_state_only_direction_matches_two_differences(case):
    _, _, args, fixed, state_rates, zeros, _ = case
    _, tangent = _fixed_point_map_block_tangent(
        args, fixed.alpha_in, state_rates, np.zeros_like(fixed.alpha_in))
    assert np.max(np.abs(tangent.differentiated_residual)) < TANGENT_RESIDUAL_GATE
    for epsilon in EPSILONS:
        delta_plus = {key: epsilon * value for key, value in state_rates.items()}
        delta_minus = {key: -epsilon * value for key, value in state_rates.items()}
        plus = _fixed_point_map_output(args, fixed.alpha_in, delta_plus)
        minus = _fixed_point_map_output(args, fixed.alpha_in, delta_minus)
        finite_difference = (plus.alpha - minus.alpha) / (2.0 * epsilon)
        np.testing.assert_allclose(
            finite_difference, tangent.alpha_dot,
            rtol=4.0e-5, atol=2.0e-8)


def test_combined_implicit_direction_matches_exact_roots_at_both_epsilons(case):
    _, _, args, fixed, state_rates, zeros, implicit = case
    _, state_tangent = _fixed_point_map_block_tangent(
        args, fixed.alpha_in, state_rates, np.zeros_like(fixed.alpha_in))
    _, lapse_tangent = _fixed_point_map_block_tangent(
        args, fixed.alpha_in, zeros, implicit.alpha_dot)
    _, combined_tangent = _fixed_point_map_block_tangent(
        args, fixed.alpha_in, state_rates, implicit.alpha_dot)
    np.testing.assert_allclose(
        combined_tangent.alpha_dot,
        state_tangent.alpha_dot + lapse_tangent.alpha_dot,
        rtol=2.0e-10, atol=2.0e-12)
    np.testing.assert_allclose(
        combined_tangent.alpha_dot, implicit.alpha_dot,
        rtol=2.0e-10, atol=2.0e-12)
    assert implicit.implicit_residual_max < TANGENT_RESIDUAL_GATE

    for epsilon in EPSILONS:
        plus_root, _ = exact_fixed_point(
            _perturbed_mp_state(args, state_rates, epsilon), dps=70)
        minus_root, _ = exact_fixed_point(
            _perturbed_mp_state(args, state_rates, -epsilon), dps=70)
        exact_fd = np.asarray([
            float((plus_root[i] - minus_root[i]) / (2 * __import__("mpmath").mp.mpf(epsilon)))
            for i in range(len(plus_root))
        ])
        assert _within_contract(implicit.alpha_dot, exact_fd), (
            f"epsilon={epsilon:g}; max_abs="
            f"{np.max(np.abs(implicit.alpha_dot-exact_fd)):.9e}; "
            f"normwise_relative={_normwise_relative(implicit.alpha_dot-exact_fd, exact_fd):.9e}"
        )


def test_normwise_relative_gate_rejects_one_e_minus_six_negative_control(case):
    _, _, _, _, _, _, implicit = case
    reference = implicit.alpha_dot
    norm = float(np.linalg.norm(reference))
    assert norm > 0.0
    direction = np.ones_like(reference) / np.sqrt(reference.size)
    negative_control = reference + 1.0e-6 * norm * direction
    assert _normwise_relative(negative_control - reference, reference) > NORMWISE_RELATIVE_GATE
    assert not _within_contract(negative_control, reference)
