from __future__ import annotations

import numpy as np
import pytest

from tools.issue4_cmc_lapse_tangent import (
    solve_cmc_discrete,
    differentiate_cmc_discrete,
    central_difference_witness,
    fourth_order_cell_derivative,
    solve_stage_aware_cmc_tangent,
    solve_stage_aware_cmc_tangent_six_sector,
    sum_six_sector_projections,
    project_six_sector_stress_from_arrays,
    independent_scalar_and_geometry_rates,
    six_sector_projection_directions,
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
    for eps in (2.0e-5,1.0e-5):
        plus=_stage_solve(base,rates,eps)
        minus=_stage_solve(base,rates,-eps)
        fd=(plus.base.alpha-minus.base.alpha)/(2.0*eps)
        np.testing.assert_allclose(fd,result.tangent.alpha_dot,rtol=3e-6,atol=2e-9)
        target_fd=(plus.base.kdot_cmc-minus.base.kdot_cmc)/(2.0*eps)
        assert target_fd==pytest.approx(result.tangent.kdot_cmc_dot,rel=3e-6,abs=2e-10)


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


def test_stage_aware_k_rhs_uses_pinned_fieldwise_l2_lapse_stencil():
    base,rates=_stage_state()
    result=_stage_solve(base,rates)
    r=base["r"]; a=base["a"]; b=base["b"]; X=base["X"]
    alpha=base["alpha_in"]; Aa=base["Aa"]; K=base["K"]
    ap=_d1(a,1); bp=_d1(b,1); Xp=_d1(X,1)
    chip=-0.5*Xp/X
    lap=(X**2/a)*(_d2(alpha,1)-_d1(alpha,1)*(
        0.5*ap/a-bp/b-2.0*chip-2.0/r))
    Q=1.5*Aa**2+K**2/3.0+4.0*np.pi*(base["rho"]+base["pr"]+2.0*base["pt"])
    np.testing.assert_allclose(result.k_rhs_in,-lap+alpha*Q,rtol=0.0,atol=2e-15)



def test_six_sector_stage_aware_binding_requires_all_named_base_and_directional_profiles():
    base,rates=_stage_state()
    names=("dark_matter","baryons","radiation","S","D","phi")
    weights=(0.20,0.10,0.05,0.30,0.05,0.30)
    sector_base={}; sector_dot={}
    for name,w in zip(names,weights):
        sector_base[name]={k:w*base[k] for k in ("rho","pr","pt")}
        sector_base[name]["j"]=np.zeros_like(base["rho"])
        sector_dot[name]={k+"_dot":w*rates[k+"_dot"] for k in ("rho","pr","pt")}
        sector_dot[name]["j_dot"]=np.zeros_like(base["rho"])
    summed=sum_six_sector_projections(
        sector_projections=sector_base,sector_projection_directions=sector_dot)
    np.testing.assert_allclose(summed["base"]["rho"],base["rho"],rtol=0,atol=1e-20)
    np.testing.assert_allclose(summed["direction"]["rho_dot"],rates["rho_dot"],rtol=0,atol=1e-20)
    stage_args={k:v for k,v in base.items() if k not in ("rho","pr","pt","k_rhs")}
    stage_args["d2"]=_d2
    stage_args.update({k:v for k,v in rates.items() if k not in
        ("rho_dot","pr_dot","pt_dot","k_rhs_dot")})
    six=solve_stage_aware_cmc_tangent_six_sector(
        sector_projections=sector_base,
        sector_projection_directions=sector_dot,**stage_args)
    direct=_stage_solve(base,rates)
    np.testing.assert_allclose(six.k_rhs_in,direct.k_rhs_in,rtol=0,atol=2e-15)
    bad=dict(sector_base); bad.pop("D")
    with pytest.raises(ValueError,match="exactly the six locked sectors"):
        sum_six_sector_projections(
            sector_projections=bad,sector_projection_directions=sector_dot)


def test_independent_six_sector_projector_reconstructs_all_locked_sectors():
    n=20; dr=0.5; r=(np.arange(n)+0.5)*dr
    a=1.0+0.01*np.exp(-(r/4.0)**2)
    b=1.0-0.005*np.exp(-(r/5.0)**2)
    X=1.0+0.003*np.exp(-(r/6.0)**2)
    sqrtg=np.sqrt(a)*b*r*r/X**3
    scalars={
        "S":0.01*np.exp(-(r/3.0)**2),
        "PS":np.full(n,2.0e-4),
        "D":0.002*np.exp(-(r/4.0)**2),
        "PD":np.full(n,1.0e-5),
        "phi":np.full(n,0.18),
        "Pi":np.full(n,3.0e-3),
    }
    fluids={}
    for name,density in (("dark_matter",1.0e-5),("baryons",2.0e-6)):
        fluids[name]={
            "rest":sqrtg*np.full(n,density),
            "energy_t":sqrtg*np.full(n,density),
            "momentum_r":np.zeros(n),
        }
    fluids["radiation"]={
        "rest":np.zeros(n),
        "energy_t":sqrtg*np.full(n,1.0e-8),
        "momentum_r":np.zeros(n),
    }
    sectors=project_six_sector_stress_from_arrays(
        r=r,a=a,b=b,X=X,scalars=scalars,fluids=fluids,
        d1=lambda x,parity:fourth_order_cell_derivative(x,dr,parity=parity,order=1))
    assert set(sectors)=={"dark_matter","baryons","radiation","S","D","phi"}
    for sector in sectors.values():
        assert set(sector)=={"rho","pr","pt","j"}
        assert all(np.all(np.isfinite(v)) for v in sector.values())


def test_pinned_instantaneous_six_sector_rates_and_projection_directions_two_epsilons():
    n=20; dr=0.5; r=(np.arange(n)+0.5)*dr
    a=1.0+0.012*np.exp(-(r/3.5)**2)
    b=1.0-0.006*np.exp(-(r/4.5)**2)
    X=1.0+0.004*np.exp(-(r/5.5)**2)
    alpha=1.0-0.025*np.exp(-(r/4.0)**2)
    beta=np.zeros(n)
    Aa=1.0e-4*np.exp(-(r/4.0)**2); K=np.full(n,-1.5e-4)
    Lambda=1.0e-5*r*np.exp(-(r/5.0)**2)
    scalars={
        "S":0.008*np.exp(-(r/3.0)**2),"PS":np.full(n,1.5e-4),
        "D":0.0015*np.exp(-(r/4.0)**2),"PD":np.full(n,8.0e-6),
        "phi":0.16+0.01*np.exp(-(r/5.0)**2),"Pi":np.full(n,2.0e-3),
    }
    sqrtg=np.sqrt(a)*b*r*r/X**3
    fluids={
        "dark_matter":{"rest":sqrtg*np.full(n,2.5e-5),
            "energy_t":-sqrtg*alpha*np.full(n,2.5e-5),
            "momentum_r":sqrtg*np.full(n,2.0e-7)},
        "baryons":{"rest":sqrtg*np.full(n,4.0e-6),
            "energy_t":-sqrtg*alpha*np.full(n,4.0e-6),
            "momentum_r":sqrtg*np.full(n,-1.0e-7)},
        "radiation":{"rest":np.zeros(n),
            "energy_t":sqrtg*np.full(n,8.0e-9),
            "momentum_r":sqrtg*np.full(n,1.0e-10)},
    }
    d1=lambda x,parity=1:fourth_order_cell_derivative(x,dr,parity=parity,order=1)
    d2=lambda x,parity=1:fourth_order_cell_derivative(x,dr,parity=parity,order=2)
    rates=independent_scalar_and_geometry_rates(
        r=r,dr=dr,a=a,b=b,X=X,alpha=alpha,beta=beta,Aa=Aa,K=K,
        Lambda=Lambda,scalars=scalars,fluids=fluids,d1=d1,d2=d2)
    assert set(rates["matter"])=={"dark_matter","baryons","radiation"}
    assert set(rates["sectors"])=={"dark_matter","baryons","radiation","S","D","phi"}
    dots=six_sector_projection_directions(
        r=r,a=a,b=b,X=X,scalars=scalars,fluids=fluids,
        geometry_rates=rates["geometry"],scalar_rates=rates["scalars"],
        matter_rates=rates["matter"],d1=d1)
    for eps in (2.0e-5,1.0e-5):
        ap=a+eps*rates["geometry"]["a"]; am=a-eps*rates["geometry"]["a"]
        bp=b+eps*rates["geometry"]["b"]; bm=b-eps*rates["geometry"]["b"]
        Xp=X+eps*rates["geometry"]["X"]; Xm=X-eps*rates["geometry"]["X"]
        sp={k:scalars[k]+eps*rates["scalars"][k] for k in scalars}
        sm={k:scalars[k]-eps*rates["scalars"][k] for k in scalars}
        fp={name:{k:fluids[name][k]+eps*rates["matter"][name][k]
                  for k in fluids[name]} for name in fluids}
        fm={name:{k:fluids[name][k]-eps*rates["matter"][name][k]
                  for k in fluids[name]} for name in fluids}
        pp=project_six_sector_stress_from_arrays(r=r,a=ap,b=bp,X=Xp,scalars=sp,fluids=fp,d1=d1)
        pm=project_six_sector_stress_from_arrays(r=r,a=am,b=bm,X=Xm,scalars=sm,fluids=fm,d1=d1)
        for sector in dots:
            for key,basekey in (("rho_dot","rho"),("pr_dot","pr"),("pt_dot","pt"),("j_dot","j")):
                fd=(pp[sector][basekey]-pm[sector][basekey])/(2.0*eps)
                np.testing.assert_allclose(fd,dots[sector][key],rtol=2.0e-4,atol=2.0e-11)
