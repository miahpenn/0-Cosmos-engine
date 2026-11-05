from __future__ import annotations

import numpy as np
import pytest

from tools.issue4_cmc_lapse_tangent import (
    solve_cmc_discrete,
    differentiate_cmc_discrete,
    central_difference_witness,
)


def _d1(values, parity=1):
    # A fixed, linear, second-order operator for the synthetic operator test.
    return np.gradient(np.asarray(values, dtype=float), 0.5, edge_order=2)


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
# Regression rerun after preserving the mass-coefficient tangent witness.
