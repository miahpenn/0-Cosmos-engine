"""Issue #4 synthetic implicit-tangent diagnostic (separate from the original finite-difference witness).

Owner decision: the original test test_self_consistent_cmc_lapse_fixed_point_and_implicit_direction_two_epsilons
(tests/test_issue4_cmc_lapse_tangent.py) is kept unchanged as a visible witness. Its finite-difference
comparison is built from float64 solver roots, so it inherits root error divided by 2*epsilon
(measured: 1.81248e-4 for this fixture). This file tests the implicit tangent against an exact-root
reference instead.

Assertions:
  1. Solver/tangent contract: tangent.implicit_residual_max < 2e-10 (as in the original test).
  2. For epsilon in {2e-5, 1e-5}, both criteria must pass against the exact-root central difference:
       (i)  absolute:  max_i |T_i - T_i^ref| <= 2e-8                  (predeclared; original absolute tolerance)
       (ii) normwise relative: ||T - T^ref||_inf / ||T^ref||_inf <= 1e-7  (predeclared; owner decision)
     The normwise relative form avoids unstable componentwise relative errors near zero components.
  2b. Committed detection check: a uniform 1e-6 relative error (absolute error ~1e-9, below the absolute
     threshold) must be rejected by the relative criterion.
  3. Exact-root central differences at the two epsilons agree to <= 2e-8 per component.

Reference: tests/data/issue4_implicit_derivative/issue4_implicit_derivative_reference.json,
generated under pinned mpmath 1.3.0 (cross-checked identical under mpmath 1.4.1).
Scope: synthetic _stage_state fixture only. Not a production-state claim.
"""
import importlib.util
import json
import pathlib

import numpy as np
import pytest

from tools.issue4_cmc_lapse_tangent import (
    differentiate_cmc_lapse_fixed_point,
    solve_cmc_lapse_fixed_point,
)

HERE = pathlib.Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location(
    "issue4_historical", HERE / "test_issue4_cmc_lapse_tangent.py")
H = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(H)

REFERENCE = HERE / "data" / "issue4_implicit_derivative" / "issue4_implicit_derivative_reference.json"
IMPLICIT_RESIDUAL_GATE = 2.0e-10     # as in the original test
ABS_THRESHOLD = 2.0e-8               # predeclared absolute threshold (owner decision)
REL_THRESHOLD = 1.0e-7               # predeclared normwise relative threshold (owner decision)


def _criteria(tangent_alpha_dot, reference):
    """Return (absolute_max, normwise_relative) for the two predeclared criteria."""
    t = np.asarray(tangent_alpha_dot, dtype=float)
    ref = np.asarray(reference, dtype=float)
    absolute = float(np.max(np.abs(t - ref)))
    normwise_rel = float(np.max(np.abs(t - ref)) / np.max(np.abs(ref)))
    return absolute, normwise_rel


def _tangent_and_args():
    base, rates = H._stage_state()
    args = {k: base[k] for k in ("r", "dr", "a", "b", "X", "Aa", "K", "rho", "pr", "pt")}
    fixed = solve_cmc_lapse_fixed_point(**args, d1=H._d1, d2=H._d2)
    tangent = differentiate_cmc_lapse_fixed_point(
        fixed_point=fixed, **args,
        a_dot=rates["a_dot"], b_dot=rates["b_dot"], X_dot=rates["X_dot"],
        Aa_dot=rates["Aa_dot"], K_dot=rates["K_dot"], rho_dot=rates["rho_dot"],
        pr_dot=rates["pr_dot"], pt_dot=rates["pt_dot"], d1=H._d1, d2=H._d2)
    return tangent


def _reference():
    data = json.loads(REFERENCE.read_text())
    assert int(data["dps"]) >= 30
    return {e: np.array([float(x) for x in data["derivative"][e]]) for e in data["epsilons"]}, data


def test_implicit_tangent_residual_contract():
    """Solver/tangent contract: the implicit tangent satisfies its residual gate."""
    tangent = _tangent_and_args()
    assert tangent.implicit_residual_max < IMPLICIT_RESIDUAL_GATE
    assert np.max(np.abs(tangent.alpha_dot)) > 0.0


@pytest.mark.parametrize("eps_key", ["2e-05", "1e-05"])
def test_tangent_matches_exact_root_directional_derivative(eps_key):
    """Tangent accuracy: both predeclared criteria, each at this epsilon."""
    ref, data = _reference()
    tangent = _tangent_and_args()
    fd = ref[eps_key]
    absolute, normwise_rel = _criteria(tangent.alpha_dot, fd)
    comp_rel = np.abs(np.asarray(tangent.alpha_dot) - fd) / np.maximum(np.abs(fd), 1e-300)
    print(f"\neps={eps_key}: absolute max = {absolute:.3e} (threshold {ABS_THRESHOLD:.0e}); "
          f"normwise relative = {normwise_rel:.3e} (threshold {REL_THRESHOLD:.0e}); "
          f"componentwise relative max = {comp_rel.max():.3e} (reported only); mpmath={data['provenance']['mpmath']}")
    assert absolute <= ABS_THRESHOLD, f"absolute criterion failed at eps={eps_key}"
    assert normwise_rel <= REL_THRESHOLD, f"normwise relative criterion failed at eps={eps_key}"


def test_relative_criterion_rejects_uniform_relative_error_absolute_accepts():
    """Committed detection check: a 1e-6 relative error must be rejected by the normwise criterion
    even though its absolute error is below the absolute threshold."""
    ref, _ = _reference()
    tangent = _tangent_and_args()
    injected = np.asarray(tangent.alpha_dot) * (1.0 + 1.0e-6)
    absolute, normwise_rel = _criteria(injected, ref["2e-05"])
    assert absolute <= ABS_THRESHOLD          # the absolute criterion alone would accept it
    assert normwise_rel > REL_THRESHOLD       # the normwise criterion must reject it


def test_exact_root_derivative_is_epsilon_independent():
    """Reference self-consistency: the exact-root central differences agree at both epsilons,
    under both predeclared criteria."""
    ref, _ = _reference()
    absolute, normwise_rel = _criteria(ref["2e-05"], ref["1e-05"])
    print(f"\nexact-root FD: absolute = {absolute:.3e}; normwise relative = {normwise_rel:.3e}")
    assert absolute <= ABS_THRESHOLD
    assert normwise_rel <= REL_THRESHOLD
