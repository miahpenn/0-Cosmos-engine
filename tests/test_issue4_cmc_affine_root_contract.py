"""Issue #4 synthetic CMC affine-root: separated contract tests.

Owner decision (test-only): the solver's declared contract and the root-accuracy question are
asserted separately.

* test_solver_residual_contract_holds: the fixed-point solver satisfies its declared residual
  gate (2e-12). This is the solver's contract.
* test_root_accuracy_within_contract_bound: the solver root agrees with an extended-precision
  reference to a bound that follows from the declared residual gate and the smallest singular
  value of the fixed-point system. No fitted constant is used.

The original strict root assertion (absolute 2e-10 agreement) is retained, unchanged, in
tests/test_issue4_cmc_lapse_tangent.py::test_fixed_point_root_matches_independent_affine_map_and_reports_conditioning
as the historical failure witness. This file does not erase that record.

Scope: synthetic _stage_state fixture only. Says nothing about production-derived states.
The reference root is a synthetic-state reference, not a production reference.
"""
import importlib.util
import json
import math
import pathlib

import numpy as np

from tools.issue4_cmc_lapse_tangent import (
    differentiate_cmc_discrete,
    physical_k_rhs_and_direction,
    solve_cmc_discrete,
    solve_cmc_lapse_fixed_point,
)

HERE = pathlib.Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location(
    "issue4_historical", HERE / "test_issue4_cmc_lapse_tangent.py")
H = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(H)

REFERENCE = HERE / "data" / "issue4_affine_root" / "issue4_affine_root_reference.json"
RESIDUAL_GATE = 2.0e-12  # solver's declared fixed-point residual gate (solve_cmc_lapse_fixed_point)


def _fixture():
    base, rates = H._stage_state()
    args = {k: base[k] for k in ("r", "dr", "a", "b", "X", "Aa", "K", "rho", "pr", "pt")}
    return base, args


def _float64_analytic_jacobian(base):
    """Replicates the solver's analytic Jacobian loop, with the solver's own tool functions."""
    n = base["r"].size
    z = np.zeros(n)
    A = {k: base[k] for k in ("a", "b", "X", "Aa", "K", "rho", "pr", "pt")}
    ref = np.ones(n)
    kr0, _ = physical_k_rhs_and_direction(
        r=base["r"], dr=base["dr"], **A, alpha_in=ref, a_dot=z, b_dot=z, X_dot=z,
        alpha_in_dot=z, Aa_dot=z, K_dot=z, rho_dot=z, pr_dot=z, pt_dot=z,
        d1=H._d1, d2=H._d2)
    rb = solve_cmc_discrete(r=base["r"], dr=base["dr"], **A, k_rhs=kr0, d1=H._d1)
    J = np.empty((n, n))
    for j in range(n):
        e = np.zeros(n); e[j] = 1.0
        _, krd = physical_k_rhs_and_direction(
            r=base["r"], dr=base["dr"], **A, alpha_in=ref, a_dot=z, b_dot=z, X_dot=z,
            alpha_in_dot=e, Aa_dot=z, K_dot=z, rho_dot=z, pr_dot=z, pt_dot=z,
            d1=H._d1, d2=H._d2)
        t = differentiate_cmc_discrete(
            base=rb, r=base["r"], dr=base["dr"], **A, k_rhs=kr0, a_dot=z, b_dot=z, X_dot=z,
            Aa_dot=z, K_dot=z, rho_dot=z, pr_dot=z, pt_dot=z, k_rhs_dot=krd, d1=H._d1)
        J[:, j] = t.alpha_dot
    return J


def test_solver_residual_contract_holds():
    """Solver contract: the fixed point satisfies its declared residual gate."""
    base, args = _fixture()
    fixed = solve_cmc_lapse_fixed_point(**args, d1=H._d1, d2=H._d2)
    # The solver raises if its gate fails; re-check the reported residual against the contract.
    assert math.isfinite(fixed.fixed_point_residual)
    assert fixed.fixed_point_residual <= RESIDUAL_GATE
    # Independent re-evaluation of the map residual at the returned root.
    out = H._fixed_point_map_output(base, fixed.alpha_in)
    map_residual = float(np.max(np.abs(out.alpha - fixed.alpha_in)))
    assert map_residual <= RESIDUAL_GATE


def test_root_accuracy_within_contract_bound():
    """Root accuracy: against the extended-precision reference, with a conditioning-aware bound.

    For the affine map F(x) = Jx + c with S = I - J, the exact root x* satisfies
    F(x) - x = -S (x - x*). Therefore
        ||x - x*||_inf <= ||x - x*||_2 <= ||S^{-1}||_2 ||F(x) - x||_2 = ||F(x) - x||_2 / sigma_min(S)
                       <= sqrt(n) * gate / sigma_min(S),
    using ||F(x) - x||_inf <= gate, which is the solver's declared contract.
    Assumptions behind the bound (each is required; none is fitted):
      1. The map is affine in the incoming lapse, F(x) = Jx + c. Checked for this fixture: the
         centred Jacobian is step-independent to ~7e-25 between h = 1e-5 and 2e-5, and the
         extended-precision and float64 Jacobians agree to ~2e-15 relative.
      2. The Jacobian J and sigma_min(I - J) refer to the same state (this fixture, same reference
         lapse = 1) and the same norm: sigma_min is the 2-norm smallest singular value, and the
         2-norm bound is converted to the infinity norm used by the gate via the factor sqrt(n).
      3. The residual gate is a pointwise (max-abs) bound, so ||F(x) - x||_inf <= gate is the
         hypothesis; the factor sqrt(n) converts the 2-norm bound back to the pointwise scale.
      4. The reference root and the float64 map use consistent boundary conventions: the inner
         row is A[0,0] = 1, A[0,1] = -1 (regularity), and the outer row pins alpha = 1
         (rhs[-1] = 1). The reference root's outer component is exactly 1.0, and the extended
         map agrees with the float64 map to ~7e-15 on the same inputs.
      5. The float64 map residual matches the exact-arithmetic residual to well within the gate.
         Measured for this fixture: 2.9e-14 float64 vs 2.9e-14 extended.
    Scope: synthetic _stage_state fixture only. Not a bound for production-derived states.
    """
    base, args = _fixture()
    fixed = solve_cmc_lapse_fixed_point(**args, d1=H._d1, d2=H._d2)
    n = base["r"].size
    J = _float64_analytic_jacobian(base)
    S = np.eye(n) - J
    sigma_min = float(np.linalg.svd(S, compute_uv=False)[-1])
    bound = math.sqrt(n) * RESIDUAL_GATE / sigma_min

    data = json.loads(REFERENCE.read_text())
    assert int(data["dps"]) >= 30
    reference = np.array([float(x) for x in data["root"]])
    err = np.abs(fixed.alpha_in - reference)
    max_err = float(np.max(err))
    cond2 = float(np.linalg.cond(S))

    # Reported regardless of pass/fail, so the conditioning and precision are documented.
    print(f"\nroot accuracy: max|x - x_ref| = {max_err:.3e}; cell0 signed = "
          f"{fixed.alpha_in[0] - reference[0]:+.3e}")
    print(f"sigma_min(I-J) = {sigma_min:.3e}; cond2 = {cond2:.3e}; contract bound = {bound:.3e}")
    print(f"reference: dps={data['dps']}, mpmath={data['provenance']['mpmath']}, "
          f"numpy={data['provenance']['numpy']}")

    assert max_err <= bound
