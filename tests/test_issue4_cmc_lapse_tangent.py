from __future__ import annotations

import numpy as np
import pytest

from tools.issue4_cmc_lapse_tangent import (
    solve_cmc_discrete,
    differentiate_cmc_discrete,
    central_difference_witness,
    fourth_order_cell_derivative,
    solve_stage_aware_cmc_tangent,
)


def _d1(values, parity=1):
    return fourth_order_cell_derivative(values, 0.5, parity=parity, order=1)


def _state():
    n=20; dr=0.5; r=(np.arange(n)+0.5)*dr
    a=1.0+0.015*np.exp(-(r/3.0)**2)
    b=1.0-0.008*np.exp(-(r/4.0)**2)
    X=1.0+0.006*np.exp(-(r/5.0)**2)
    Aa=1.0e-4*np.exp(-(r/4.0)**2)
    K=np.full(n,-2.0e-4)
    rho=2.0e-5*np.exp(-(r/6.0)**2)
    pr=0.2*rho
    pt=0.15*rho
    kr=np.full(n,2.0e-5)
    rates={
        "a_dot": 1.0e-4*np.exp(-(r/3.5)**2),
        "b_dot": -0.7e-4*np.exp(-(r/4.5)**2),
        "X_dot": 0.5e-4*np.exp(-(r/5.5)**2),
        "Aa_dot": 0.3e-4*np.exp(-(r/4.0)**2),
        "K_dot": -0.2e-4*np.exp(-(r/4.0)**2),
        "rho_dot": 0.1*rho,
        "pr_dot": 0.2*pr,
        "pt_dot": -0.1*pt,
        "k_rhs_dot": 0.15*kr,
    }
    base=dict(r=r,dr=dr,a=a,b=b,X=X,Aa=Aa,K=K,rho=rho,pr=pr,pt=pt,k_rhs=kr,d1=_d1)
    return base,rates


def _solve_perturbed(base, rates, eps):
    values={}
    for name in ("a","b","X","Aa","K","rho","pr","pt","k_rhs"):
        values[name]=base[name]+eps*rates[name+"_dot"]
    return solve_cmc_discrete(**{**base,**values})


def test_independent_discrete_cmc_solve_satisfies_all_rows_and_positive_lapse():
    base,_=_state()
    result=solve_cmc_discrete(**base)
    assert np.min(result.alpha)>0.0
    assert result.alpha[0]==pytest.approx(result.alpha[1],abs=1e-13)
    assert result.alpha[-1]==pytest.approx(1.0,abs=1e-13)
    assert np.max(np.abs(result.residual))<1e-12
    assert result.operator.shape==(base["r"].size,base["r"].size)


def test_differentiated_cmc_equation_matches_independent_plus_minus_solves():
    base,rates=_state()
    result=solve_cmc_discrete(**base)
    tangent=differentiate_cmc_discrete(
        base=result, r=base["r"], dr=base["dr"],
        **{k:base[k] for k in ("a","b","X","Aa","K","rho","pr","pt","k_rhs")},
        **rates, d1=base["d1"],
    )
    assert np.max(np.abs(tangent.differentiated_residual))<1e-11
    assert tangent.alpha_dot[0]==pytest.approx(tangent.alpha_dot[1],abs=1e-12)
    assert tangent.alpha_dot[-1]==pytest.approx(0.0,abs=1e-12)

    eps=2.0e-5
    plus=_solve_perturbed(base,rates,eps)
    minus=_solve_perturbed(base,rates,-eps)
    witness=central_difference_witness(plus,minus,eps,tangent.alpha_dot)
    assert witness["max_abs_difference"]<2e-8


def test_cmc_tangent_includes_proper_volume_projection_response():
    base,rates=_state()
    result=solve_cmc_discrete(**base)
    tangent=differentiate_cmc_discrete(
        base=result,r=base["r"],dr=base["dr"],
        **{k:base[k] for k in ("a","b","X","Aa","K","rho","pr","pt","k_rhs")},
        **rates,d1=base["d1"],
    )
    assert np.max(np.abs(tangent.weight_dot))>0.0
    assert tangent.proper_volume_total_dot!=0.0
    assert tangent.kdot_cmc_dot!=0.0


def test_cmc_solver_rejects_nonpositive_metric_and_nonfinite_projection():
    base,_=_state()
    bad=dict(base); bad["a"]=base["a"].copy(); bad["a"][3]=-1.0
    with pytest.raises(FloatingPointError,match="metric coefficients must be positive"):
        solve_cmc_discrete(**bad)
    bad=dict(base); bad["k_rhs"]=base["k_rhs"].copy(); bad["k_rhs"][4]=np.nan
    with pytest.raises(FloatingPointError,match="K RHS contains non-finite"):
        solve_cmc_discrete(**bad)

# CI gate: independent CMC discrete solve and differentiated operator contract.
# Regression rerun after preserving the mass-coefficient and potential tangents.


def test_independent_derivative_transcription_matches_even_and_odd_polynomial_stencils():
    dr=0.5; r=(np.arange(12)+0.5)*dr
    even=fourth_order_cell_derivative(r**2,dr,parity=1,order=1)
    odd=fourth_order_cell_derivative(r,dr,parity=-1,order=1)
    np.testing.assert_allclose(even,2.0*r,rtol=0.0,atol=2e-12)
    np.testing.assert_allclose(odd,np.ones_like(r),rtol=0.0,atol=2e-12)


def _d2(values, parity=1):
    return fourth_order_cell_derivative(values, 0.5, parity=parity, order=2)


def _stage_state():
    base,rates=_state()
    r=base["r"]
    base["alpha_in"]=1.0-0.04*np.exp(-(r/4.0)**2)
    rates["alpha_in_dot"]=0.2e-4*np.exp(-(r/3.0)**2)
    return base,rates


def _stage_solve(base,rates,eps=0.0):
    values={k:base[k]+eps*rates[k+"_dot"] for k in
        ("a","b","X","Aa","K","rho","pr","pt","alpha_in")}
    return solve_stage_aware_cmc_tangent(
        r=base["r"],dr=base["dr"],**values,
        **{k:v for k,v in rates.items() if k != 'k_rhs_dot'},
        d1=_d1,d2=_d2)


def test_stage_aware_physical_k_rhs_and_lapse_tangent_match_plus_minus_solves():
    from tools.issue4_cmc_lapse_tangent import solve_stage_aware_cmc_tangent
    base,rates=_stage_state()
    result=_stage_solve(base,rates)
    assert np.min(result.base.alpha)>0.0
    assert np.max(np.abs(result.base.residual))<1e-11
    assert np.max(np.abs(result.tangent.differentiated_residual))<1e-10
    assert np.all(np.isfinite(result.k_rhs_in))
    assert np.all(np.isfinite(result.k_rhs_in_dot))
    eps=2.0e-5
    plus=_stage_solve(base,rates,eps)
    minus=_stage_solve(base,rates,-eps)
    fd=(plus.base.alpha-minus.base.alpha)/(2.0*eps)
    np.testing.assert_allclose(fd,result.tangent.alpha_dot,rtol=2e-6,atol=2e-9)
    target_fd=(plus.base.kdot_cmc-minus.base.kdot_cmc)/(2.0*eps)
    assert target_fd==pytest.approx(result.tangent.kdot_cmc_dot,rel=2e-6,abs=2e-10)


def test_stage_aware_physical_k_rhs_direction_is_not_an_external_input():
    import inspect
    from tools.issue4_cmc_lapse_tangent import (
        physical_k_rhs_and_direction, solve_stage_aware_cmc_tangent)
    assert "k_rhs_dot" not in inspect.signature(physical_k_rhs_and_direction).parameters
    assert "k_rhs_dot" not in inspect.signature(solve_stage_aware_cmc_tangent).parameters
    base,rates=_stage_state()
    result=_stage_solve(base,rates)
    assert np.max(np.abs(result.k_rhs_in_dot))>0.0
    assert np.max(np.abs(result.tangent.alpha_dot))>0.0
