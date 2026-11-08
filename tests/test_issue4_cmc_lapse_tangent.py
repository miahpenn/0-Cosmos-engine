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



def test_self_consistent_cmc_lapse_fixed_point_and_implicit_direction_two_epsilons():
    from tools.issue4_cmc_lapse_tangent import (
        solve_cmc_lapse_fixed_point,
        differentiate_cmc_lapse_fixed_point,
    )
    base,rates=_stage_state()
    args={k:base[k] for k in
        ("r","dr","a","b","X","Aa","K","rho","pr","pt")}
    fixed=solve_cmc_lapse_fixed_point(**args,d1=_d1,d2=_d2)
    assert fixed.linearization_evaluations==base["r"].size+1
    assert fixed.fixed_point_residual<2.0e-13
    np.testing.assert_allclose(fixed.alpha_in,fixed.base.alpha,rtol=0.0,atol=2.0e-13)

    tangent=differentiate_cmc_lapse_fixed_point(
        fixed_point=fixed,**args,
        a_dot=rates["a_dot"],b_dot=rates["b_dot"],X_dot=rates["X_dot"],
        Aa_dot=rates["Aa_dot"],K_dot=rates["K_dot"],
        rho_dot=rates["rho_dot"],pr_dot=rates["pr_dot"],pt_dot=rates["pt_dot"],
        d1=_d1,d2=_d2)
    assert tangent.implicit_residual_max<2.0e-10
    assert np.max(np.abs(tangent.alpha_dot))>0.0

    for eps in (2.0e-5,1.0e-5):
        perturbed={}
        for key in ("a","b","X","Aa","K","rho","pr","pt"):
            perturbed[key]=base[key]+eps*rates[key+"_dot"]
        plus=solve_cmc_lapse_fixed_point(
            r=base["r"],dr=base["dr"],**perturbed,d1=_d1,d2=_d2)
        perturbed={}
        for key in ("a","b","X","Aa","K","rho","pr","pt"):
            perturbed[key]=base[key]-eps*rates[key+"_dot"]
        minus=solve_cmc_lapse_fixed_point(
            r=base["r"],dr=base["dr"],**perturbed,d1=_d1,d2=_d2)
        fd=(plus.alpha_in-minus.alpha_in)/(2.0*eps)
        np.testing.assert_allclose(fd,tangent.alpha_dot,rtol=4.0e-5,atol=2.0e-8)


def test_independent_derivative_transcription_matches_even_and_odd_polynomial_stencils():
    dr=0.5; r=(np.arange(12)+0.5)*dr
    even=fourth_order_cell_derivative(r**2,dr,parity=1,order=1)
    odd=fourth_order_cell_derivative(r,dr,parity=-1,order=1)
    np.testing.assert_allclose(even,2.0*r,rtol=0.0,atol=2e-12)
    np.testing.assert_allclose(odd,np.ones_like(r),rtol=0.0,atol=2e-12)



def test_fixed_point_map_combined_direction_is_sum_of_independent_blocks_two_epsilons():
    """Check lapse/state additivity before interpreting the implicit fixed-point mismatch."""
    from tools.issue4_cmc_lapse_tangent import solve_cmc_lapse_fixed_point
    base,rates=_stage_state()
    r=base["r"]
    _,_,summed,physical_rates=_physical_six_sector_split(
        r=r,dr=base["dr"],a=base["a"],b=base["b"],X=base["X"],
        Aa=base["Aa"],K=base["K"],alpha=base["alpha_in"])
    args={key:base[key] for key in ("r","dr","a","b","X","Aa","K")}
    args.update({key:summed["base"][key] for key in ("rho","pr","pt")})
    fixed=solve_cmc_lapse_fixed_point(**args,d1=_d1,d2=_d2)

    # Recompute the physical instantaneous state direction at the fixed-point
    # incoming lapse, keeping the independently projected base state fixed.
    _,_,summed,physical_rates=_physical_six_sector_split(
        r=r,dr=base["dr"],a=base["a"],b=base["b"],X=base["X"],
        Aa=base["Aa"],K=base["K"],alpha=fixed.alpha_in)
    state_rates={key+"_dot":physical_rates["geometry"][key]
                 for key in ("a","b","X","Aa","K")}
    for key in ("rho","pr","pt"):
        state_rates[key+"_dot"]=summed["direction"][key+"_dot"]
    zeros={key:np.zeros_like(r) for key in
           ("a_dot","b_dot","X_dot","Aa_dot","K_dot",
            "rho_dot","pr_dot","pt_dot")}
    alpha_direction=rates["alpha_in_dot"]

    _,lapse_tangent=_fixed_point_map_block_tangent(
        args,fixed.alpha_in,zeros,alpha_direction)
    _,state_tangent=_fixed_point_map_block_tangent(
        args,fixed.alpha_in,state_rates,np.zeros_like(r))
    combined_rates={key:state_rates[key] for key in state_rates}
    _,combined_tangent=_fixed_point_map_block_tangent(
        args,fixed.alpha_in,combined_rates,alpha_direction)

    # The derivative of this map must be additive in independently varied inputs.
    np.testing.assert_allclose(
        combined_tangent.alpha_dot,
        lapse_tangent.alpha_dot+state_tangent.alpha_dot,
        rtol=2.0e-10,atol=2.0e-12)

    for eps in (2.0e-5,1.0e-5):
        plus=_fixed_point_map_output(
            args,fixed.alpha_in+eps*alpha_direction,
            {name.removesuffix("_dot"):eps*value
             for name,value in combined_rates.items()})
        minus=_fixed_point_map_output(
            args,fixed.alpha_in-eps*alpha_direction,
            {name.removesuffix("_dot"):-eps*value
             for name,value in combined_rates.items()})
        fd=(plus.alpha-minus.alpha)/(2.0*eps)
        np.testing.assert_allclose(
            fd,combined_tangent.alpha_dot,rtol=4.0e-5,atol=2.0e-8)


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


def _physical_six_sector_split(*, r, dr, a, b, X, Aa, K, alpha):
    """Build base and directional profiles through the independent six-sector projector."""
    n=r.size
    Lambda=1.0e-5*r*np.exp(-(r/5.0)**2)
    beta=np.zeros(n)
    scalars={
        "S":0.008*np.exp(-(r/3.0)**2),"PS":np.full(n,1.5e-4),
        "D":0.0015*np.exp(-(r/4.0)**2),"PD":np.full(n,8.0e-6),
        "phi":0.16+0.01*np.exp(-(r/5.0)**2),"Pi":np.full(n,2.0e-3),
    }
    sqrtg=np.sqrt(a)*b*r*r/(X**3)
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
    sectors=project_six_sector_stress_from_arrays(
        r=r,a=a,b=b,X=X,scalars=scalars,fluids=fluids,d1=_d1)
    rates=independent_scalar_and_geometry_rates(
        r=r,dr=dr,a=a,b=b,X=X,alpha=alpha,beta=beta,Aa=Aa,K=K,
        Lambda=Lambda,scalars=scalars,fluids=fluids,d1=_d1,d2=_d2)
    directions=six_sector_projection_directions(
        r=r,a=a,b=b,X=X,scalars=scalars,fluids=fluids,
        geometry_rates=rates["geometry"],scalar_rates=rates["scalars"],
        matter_rates=rates["matter"],d1=_d1)
    summed=sum_six_sector_projections(
        sector_projections=sectors,sector_projection_directions=directions)
    return sectors,directions,summed,rates


def _d2(values, parity=1):
    return fourth_order_cell_derivative(values, 0.5, parity=parity, order=2)

def _fixed_point_map_output(base, alpha_in, state_delta=None):
    """Evaluate the independent input-lapse -> output-lapse map."""
    from tools.issue4_cmc_lapse_tangent import (
        physical_k_rhs_and_direction, solve_cmc_discrete)
    state_delta={} if state_delta is None else state_delta
    fields={key:base[key]+state_delta.get(key,np.zeros_like(base[key]))
            for key in ("a","b","X","Aa","K","rho","pr","pt")}
    zeros=np.zeros_like(base["r"])
    kr,_=physical_k_rhs_and_direction(
        r=base["r"],dr=base["dr"],a=fields["a"],b=fields["b"],X=fields["X"],
        alpha_in=alpha_in,Aa=fields["Aa"],K=fields["K"],
        rho=fields["rho"],pr=fields["pr"],pt=fields["pt"],
        a_dot=zeros,b_dot=zeros,X_dot=zeros,alpha_in_dot=zeros,
        Aa_dot=zeros,K_dot=zeros,rho_dot=zeros,pr_dot=zeros,pt_dot=zeros,
        d1=_d1,d2=_d2)
    return solve_cmc_discrete(
        r=base["r"],dr=base["dr"],a=fields["a"],b=fields["b"],X=fields["X"],
        Aa=fields["Aa"],K=fields["K"],rho=fields["rho"],pr=fields["pr"],
        pt=fields["pt"],k_rhs=kr,d1=_d1)


def _fixed_point_map_block_tangent(base, alpha_in, state_rates, alpha_in_dot):
    """Differentiate one independently selected input block of the output map."""
    from tools.issue4_cmc_lapse_tangent import (
        physical_k_rhs_and_direction, differentiate_cmc_discrete)
    fields={key:base[key] for key in ("a","b","X","Aa","K","rho","pr","pt")}
    zeros=np.zeros_like(base["r"])
    kr,kr_dot=physical_k_rhs_and_direction(
        r=base["r"],dr=base["dr"],a=fields["a"],b=fields["b"],X=fields["X"],
        alpha_in=alpha_in,Aa=fields["Aa"],K=fields["K"],
        rho=fields["rho"],pr=fields["pr"],pt=fields["pt"],
        a_dot=state_rates["a_dot"],b_dot=state_rates["b_dot"],
        X_dot=state_rates["X_dot"],alpha_in_dot=alpha_in_dot,
        Aa_dot=state_rates["Aa_dot"],K_dot=state_rates["K_dot"],
        rho_dot=state_rates["rho_dot"],pr_dot=state_rates["pr_dot"],
        pt_dot=state_rates["pt_dot"],d1=_d1,d2=_d2)
    output=_fixed_point_map_output(base,alpha_in)
    tangent=differentiate_cmc_discrete(
        base=output,r=base["r"],dr=base["dr"],**fields,k_rhs=kr,
        a_dot=state_rates["a_dot"],b_dot=state_rates["b_dot"],
        X_dot=state_rates["X_dot"],Aa_dot=state_rates["Aa_dot"],
        K_dot=state_rates["K_dot"],rho_dot=state_rates["rho_dot"],
        pr_dot=state_rates["pr_dot"],pt_dot=state_rates["pt_dot"],
        k_rhs_dot=kr_dot,d1=_d1)
    return output,tangent




def test_stage_fixture_fixed_point_map_state_and_combined_blocks_two_epsilons():
    """Cross-check the original stage fixture's state block and implicit combined direction."""
    from tools.issue4_cmc_lapse_tangent import (
        solve_cmc_lapse_fixed_point,differentiate_cmc_lapse_fixed_point)
    base,rates=_stage_state()
    args={key:base[key] for key in
          ("r","dr","a","b","X","Aa","K","rho","pr","pt")}
    fixed=solve_cmc_lapse_fixed_point(**args,d1=_d1,d2=_d2)
    state_rates={key+"_dot":rates[key+"_dot"] for key in
                 ("a","b","X","Aa","K","rho","pr","pt")}
    zeros=np.zeros_like(base["r"])
    _,state_tangent=_fixed_point_map_block_tangent(
        args,fixed.alpha_in,state_rates,zeros)
    implicit=differentiate_cmc_lapse_fixed_point(
        fixed_point=fixed,**args,
        a_dot=state_rates["a_dot"],b_dot=state_rates["b_dot"],
        X_dot=state_rates["X_dot"],Aa_dot=state_rates["Aa_dot"],
        K_dot=state_rates["K_dot"],rho_dot=state_rates["rho_dot"],
        pr_dot=state_rates["pr_dot"],pt_dot=state_rates["pt_dot"],
        d1=_d1,d2=_d2)
    _,combined_tangent=_fixed_point_map_block_tangent(
        args,fixed.alpha_in,state_rates,implicit.alpha_dot)
    np.testing.assert_allclose(
        combined_tangent.alpha_dot,implicit.alpha_dot,
        rtol=2.0e-10,atol=2.0e-12)
    for eps in (2.0e-5,1.0e-5):
        plus_state={key:eps*state_rates[key+"_dot"] for key in
                    ("a","b","X","Aa","K","rho","pr","pt")}
        minus_state={key:-eps*state_rates[key+"_dot"] for key in
                     ("a","b","X","Aa","K","rho","pr","pt")}
        plus=_fixed_point_map_output(args,fixed.alpha_in,plus_state)
        minus=_fixed_point_map_output(args,fixed.alpha_in,minus_state)
        state_fd=(plus.alpha-minus.alpha)/(2.0*eps)
        np.testing.assert_allclose(
            state_fd,state_tangent.alpha_dot,rtol=4.0e-5,atol=2.0e-8)
        plus=_fixed_point_map_output(
            args,fixed.alpha_in+eps*implicit.alpha_dot,plus_state)
        minus=_fixed_point_map_output(
            args,fixed.alpha_in-eps*implicit.alpha_dot,minus_state)
        combined_fd=(plus.alpha-minus.alpha)/(2.0*eps)
        np.testing.assert_allclose(
            combined_fd,implicit.alpha_dot,rtol=4.0e-5,atol=2.0e-8)


def test_fixed_point_map_full_incoming_lapse_jacobian_basis_two_epsilons():
    """Compare every analytic incoming-lapse basis column with direct map differences."""
    from tools.issue4_cmc_lapse_tangent import solve_cmc_lapse_fixed_point
    base,_=_stage_state()
    args={key:base[key] for key in
          ("r","dr","a","b","X","Aa","K","rho","pr","pt")}
    fixed=solve_cmc_lapse_fixed_point(**args,d1=_d1,d2=_d2)
    zeros={key:np.zeros_like(base["r"]) for key in
           ("a_dot","b_dot","X_dot","Aa_dot","K_dot",
            "rho_dot","pr_dot","pt_dot")}
    n=base["r"].size
    for eps in (2.0e-5,1.0e-5):
        for j in range(n):
            basis=np.zeros(n); basis[j]=1.0
            _,analytic=_fixed_point_map_block_tangent(
                args,fixed.alpha_in,zeros,basis)
            plus=_fixed_point_map_output(args,fixed.alpha_in+eps*basis)
            minus=_fixed_point_map_output(args,fixed.alpha_in-eps*basis)
            fd=(plus.alpha-minus.alpha)/(2.0*eps)
            np.testing.assert_allclose(
                fd,analytic.alpha_dot,rtol=4.0e-5,atol=2.0e-8)


def test_fixed_point_map_incoming_lapse_block_only_two_epsilons():
    """Incoming-lapse derivative: geometry and all six matter sectors fixed."""
    from tools.issue4_cmc_lapse_tangent import solve_cmc_lapse_fixed_point
    base,rates=_stage_state()
    args={key:base[key] for key in
          ("r","dr","a","b","X","Aa","K","rho","pr","pt")}
    fixed=solve_cmc_lapse_fixed_point(**args,d1=_d1,d2=_d2)
    zeros={key:np.zeros_like(base["r"]) for key in
           ("a_dot","b_dot","X_dot","Aa_dot","K_dot",
            "rho_dot","pr_dot","pt_dot")}
    alpha_direction=rates["alpha_in_dot"]
    output,tangent=_fixed_point_map_block_tangent(
        base,fixed.alpha_in,zeros,alpha_direction)
    assert np.max(np.abs(tangent.differentiated_residual))<1.0e-10
    for eps in (2.0e-5,1.0e-5):
        plus=_fixed_point_map_output(
            base,fixed.alpha_in+eps*alpha_direction)
        minus=_fixed_point_map_output(
            base,fixed.alpha_in-eps*alpha_direction)
        fd=(plus.alpha-minus.alpha)/(2.0*eps)
        np.testing.assert_allclose(
            fd,tangent.alpha_dot,rtol=4.0e-5,atol=2.0e-8)


def test_fixed_point_map_geometry_and_six_sector_state_block_only_two_epsilons():
    """State derivative uses independently projected six-sector base and tangent profiles."""
    from tools.issue4_cmc_lapse_tangent import solve_cmc_lapse_fixed_point
    base,rates=_stage_state()
    r=base["r"]
    sectors,directions,summed,_=_physical_six_sector_split(
        r=r,dr=base["dr"],a=base["a"],b=base["b"],X=base["X"],
        Aa=base["Aa"],K=base["K"],alpha=base["alpha_in"])
    args={key:base[key] for key in
          ("r","dr","a","b","X","Aa","K")}
    args.update({key:summed["base"][key] for key in ("rho","pr","pt")})
    fixed=solve_cmc_lapse_fixed_point(**args,d1=_d1,d2=_d2)
    # Recompute the independent instantaneous directions at the fixed-point
    # incoming lapse; the projection base remains the same physical six-sector state.
    sectors,directions,summed,rates=_physical_six_sector_split(
        r=r,dr=base["dr"],a=base["a"],b=base["b"],X=base["X"],
        Aa=base["Aa"],K=base["K"],alpha=fixed.alpha_in)
    state_rates={key+"_dot":rates["geometry"][key] for key in
                 ("a","b","X","Aa","K")}
    for key in ("rho","pr","pt"):
        state_rates[key+"_dot"]=summed["direction"][key+"_dot"]
    zeros=np.zeros_like(r)
    output,tangent=_fixed_point_map_block_tangent(
        args,fixed.alpha_in,state_rates,zeros)
    assert np.max(np.abs(tangent.differentiated_residual))<1.0e-10
    delta={key:state_rates[key+"_dot"] for key in
           ("a","b","X","Aa","K","rho","pr","pt")}
    for eps in (2.0e-5,1.0e-5):
        plus=_fixed_point_map_output(
            args,fixed.alpha_in,
            {key:eps*value for key,value in delta.items()})
        minus=_fixed_point_map_output(
            args,fixed.alpha_in,
            {key:-eps*value for key,value in delta.items()})
        fd=(plus.alpha-minus.alpha)/(2.0*eps)
        np.testing.assert_allclose(
            fd,tangent.alpha_dot,rtol=4.0e-5,atol=2.0e-8)


def test_fixed_point_root_matches_independent_affine_map_and_reports_conditioning():
    """Separate affine-root construction from implicit tangent conditioning."""
    from tools.issue4_cmc_lapse_tangent import solve_cmc_lapse_fixed_point
    base,rates=_stage_state()
    keys=("a","b","X","Aa","K","rho","pr","pt")
    args={key:base[key] for key in
          ("r","dr","a","b","X","Aa","K","rho","pr","pt")}
    fixed=solve_cmc_lapse_fixed_point(**args,d1=_d1,d2=_d2)
    n=base["r"].size
    reference=np.ones(n)
    eye=np.eye(n)
    epsilons=(2.0e-5,1.0e-5)

    def affine_root(state_delta):
        shifted=dict(args)
        for key in keys:
            shifted[key]=args[key]+state_delta.get(key,np.zeros(n))
        fref=_fixed_point_map_output(shifted,reference).alpha
        jac=np.empty((n,n))
        for j in range(n):
            step=1.0e-5
            basis=step*eye[j]
            plus=_fixed_point_map_output(shifted,reference+basis).alpha
            minus=_fixed_point_map_output(shifted,reference-basis).alpha
            jac[:,j]=(plus-minus)/(2.0*step)
        system=np.eye(n)-jac
        root=np.linalg.solve(system,fref-jac@reference)
        return root,float(np.linalg.cond(system))

    base_root,base_condition=affine_root({})
    np.testing.assert_allclose(base_root,fixed.alpha_in,rtol=2.0e-10,atol=2.0e-12)
    print("fixed-point affine-system condition number:",base_condition)

    for eps in epsilons:
        plus_delta={key:eps*rates[key+"_dot"] for key in keys}
        minus_delta={key:-eps*rates[key+"_dot"] for key in keys}
        plus_root,plus_condition=affine_root(plus_delta)
        minus_root,minus_condition=affine_root(minus_delta)
        plus_args={key:args[key]+plus_delta.get(key,np.zeros(n)) for key in keys}
        minus_args={key:args[key]+minus_delta.get(key,np.zeros(n)) for key in keys}
        plus=solve_cmc_lapse_fixed_point(
            r=base["r"],dr=base["dr"],**plus_args,d1=_d1,d2=_d2)
        minus=solve_cmc_lapse_fixed_point(
            r=base["r"],dr=base["dr"],**minus_args,d1=_d1,d2=_d2)
        np.testing.assert_allclose(plus.alpha_in,plus_root,rtol=2.0e-10,atol=2.0e-12)
        np.testing.assert_allclose(minus.alpha_in,minus_root,rtol=2.0e-10,atol=2.0e-12)
        fd=(plus.alpha_in-minus.alpha_in)/(2.0*eps)
        print(
            f"epsilon={eps:g}; plus/minus affine-system condition numbers="
            f"{plus_condition:.9e}/{minus_condition:.9e}; "
            f"root-FD max={np.max(np.abs(fd)):.9e}")
