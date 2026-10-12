#!/usr/bin/env python3
"""Diagnostic A/B test of the regular-origin half-cell integral in v55_initial.

Baseline exactly reproduces the pinned initializer's first integral entry:
(Delta r/2)*source[0]. The counterfactual replaces only that half-cell integral
with the three-node Lagrange integral in r^2 through source(r0), source(r1),
source(r2); weights are derived analytically as 401/720, -31/480, 11/1440.
That quadrature is exact for even source polynomials through r^4. All subsequent
center-to-center trapezoids, six nonlinear iterations, field definitions,
matter initialization, CMC lapse solve, and constraints are held unchanged.

No production file is edited and no state is evolved. This is a numerical
initial-data counterfactual, not yet a proposed production repair.
"""
from __future__ import annotations
import json
import math
from pathlib import Path
import sys
import numpy as np

sys.path.insert(0, str(Path("source").resolve()))
from engine import v55_pirk_adapter as adapter
from engine.production_kernel import V55ProductionKernel, ProductionState, LAMBDA_M
from engine.scalar_system import ScalarFields, cosmos_potential
from engine.v55_matter import (
    RHO_DM_PRESENT, RHO_B_PRESENT, RHO_R_PRESENT,
    initialize_from_archive,
)
from engine.stress_energy import assemble_total_stress_energy

SOURCE_COMMIT = "6d8b3f2790644591326df57c55908863f994cb0c"
RESOLUTIONS = (40, 60, 80, 120)
R_MAX = 80.0
AMPLITUDE = 0.01
WIDTH = 7.0
D_AMPLITUDE = 1.0e-10
ANNULUS = (2.0, 5.0)
CENTRE_CELLS = (0, 1, 2, 3, 4)
FIXED_RADII = (2.0, 2.5, 3.0, 4.0, 5.0)
ORIGIN_WEIGHTS = (401.0/720.0, -31.0/480.0, 11.0/1440.0)


def d1(grid, x, parity=+1):
    return np.asarray(grid.cell_derivative_fourth(x, parity=parity), dtype=float)


def counterfactual_build_initial_data(grid, vacuum_state_factory):
    """Exact copy of pinned build_initial_data except the origin half-cell rule."""
    r = np.asarray(grid.centers, dtype=float)
    n = len(r)
    dr = float(grid.dr)
    rho_r = RHO_R_PRESENT

    phi0 = 0.179055
    pi0 = 0.00341913
    dm0 = RHO_DM_PRESENT
    b0 = RHO_B_PRESENT

    rho0 = 0.5 * pi0**2 + float(cosmos_potential(np.array([phi0]))[0])
    rho0 += dm0 + b0 + rho_r
    H0 = math.sqrt(max(rho0 / 3.0, 0.0))

    S = AMPLITUDE * np.exp(-(r / WIDTH) ** 2)
    D = D_AMPLITUDE * np.exp(-(r / WIDTH) ** 2)
    PS = np.zeros_like(r)
    PD = D.copy()
    phi = np.full(n, phi0)
    Pi = np.full(n, pi0)
    Kt = np.full(n, -H0)

    Sp = d1(grid, S, +1)
    Dp = d1(grid, D, +1)
    Ktp = d1(grid, Kt, +1)
    j = -(PD * Dp)
    B = np.ones_like(r)

    for _ in range(6):
        rhoS = 0.5 * (PS**2 + B * Sp**2) + 0.5 * S**2
        rhoD = 0.5 * (PD**2 + B * Dp**2) - 0.5 * D**2
        rhoPhi = (0.5 * Pi**2 + cosmos_potential(phi)) / (8.0 * math.pi)
        rhoMatter = (dm0 + b0 + rho_r) / (8.0 * math.pi)
        rho = rhoS + rhoD + rhoPhi + rhoMatter
        Kr = Kt + r * Ktp + 4.0 * math.pi * r * j
        src = (
            1.0 - 8.0 * math.pi * r**2 * rho
            + 2.0 * r**2 * Kr * Kt
            + r**2 * Kt**2
        )

        if n < 3:
            raise ValueError("origin quadrature requires at least three radial cells")
        half_cell = dr * sum(ORIGIN_WEIGHTS[k] * src[k] for k in range(3))
        increments = 0.5 * dr * (src[1:] + src[:-1])
        integral = np.cumsum(np.r_[half_cell, increments])
        B = integral / r

    if np.any(~np.isfinite(B)) or np.any(B <= 0.0):
        raise FloatingPointError("counterfactual quadrature produced invalid B")

    A = 1.0 / np.sqrt(B)
    Kr = Kt + r * Ktp + 4.0 * math.pi * r * j
    K = Kr + 2.0 * Kt
    Aa = Kr - K / 3.0
    a = A ** (4.0 / 3.0)
    b = A ** (-2.0 / 3.0)
    X = A ** (-1.0 / 3.0)
    Lambda = (
        d1(grid, a, +1) / (2.0*a*a)
        - d1(grid, b, +1) / (a*b)
        + 2.0/r * (1.0/b - 1.0/a)
    )

    geom = vacuum_state_factory(grid)
    geom.a[:] = a
    geom.b[:] = b
    geom.X[:] = X
    geom.alpha[:] = 1.0
    geom.beta[:] = 0.0
    geom.Aa[:] = Aa
    geom.K[:] = K
    geom.Lambda[:] = Lambda
    geom.B[:] = 0.75*Lambda

    scalars = ScalarFields(S, PS, D, PD, phi, Pi)
    matter = initialize_from_archive(grid, geom, include_radiation=True)
    return geom, scalars, matter, B, src


def independent_metric_ricci(grid, g):
    r = np.asarray(grid.centers, dtype=float)
    areal = r*np.sqrt(np.asarray(g.b))/np.asarray(g.X)
    d_ell_dr = np.asarray(g.X)/np.sqrt(np.asarray(g.a))
    areal_ell = d_ell_dr*d1(grid, areal, -1)
    areal_ell_ell = d_ell_dr*d1(grid, areal_ell, +1)
    rr = -2.0*areal_ell_ell/areal
    angular = (1.0-areal_ell**2-areal*areal_ell_ell)/areal**2
    scalar = rr+2.0*angular
    return {"Rr_r":rr,"Rtheta_theta":angular,"R":scalar,"R_tf":rr-scalar/3.0}


def maxabs(x):
    return float(np.max(np.abs(np.asarray(x,dtype=float))))


def rms(x):
    x=np.asarray(x,dtype=float)
    return float(np.sqrt(np.mean(x*x)))


def evaluate(state, label, B_profile, src_profile):
    grid,g=state.grid,state.geometry
    r=np.asarray(grid.centers,dtype=float)
    _,vacuum,_=adapter.vendor_modules()
    metric_lam=(d1(grid,g.a,+1)/(2*np.asarray(g.a)**2)
         -d1(grid,g.b,+1)/(np.asarray(g.a)*np.asarray(g.b))
         +2.0/r*(1.0/np.asarray(g.b)-1.0/np.asarray(g.a)))
    metric_geom=g.copy()
    metric_geom.Lambda=metric_lam.copy()
    terms=vacuum.geometry_terms(grid,g)
    metric_terms=vacuum.geometry_terms(grid,metric_geom)
    metric_ricci=independent_metric_ricci(grid,g)
    tf_vendor=np.asarray(terms["Rr_r"])-np.asarray(terms["R"])/3.0
    tf_vendor_metric_lam=np.asarray(metric_terms["Rr_r"])-np.asarray(metric_terms["R"])/3.0
    tf_metric=np.asarray(metric_ricci["R_tf"])
    conn_piece=tf_vendor-tf_vendor_metric_lam
    route_piece=tf_vendor_metric_lam-tf_metric
    constraints=vacuum.constraints(grid,g)
    total=assemble_total_stress_energy(grid,g,state.scalars,state.matter)
    H=np.asarray(constraints["hamiltonian"])-16.0*math.pi*np.asarray(total.rho)
    M=np.asarray(constraints["momentum"])-8.0*math.pi*np.asarray(total.j)
    C=np.asarray(g.Lambda)-metric_lam

    stages=adapter.geometry_stage_terms(grid,g,state.scalars,state.matter,lambda_m=LAMBDA_M)
    aa_lapse=-(np.asarray(terms["DrDr_alpha"])-np.asarray(terms["lap_alpha"])/3.0)
    aa_tf=np.asarray(g.alpha)*tf_vendor
    aa_rhs=np.asarray(stages["primary_l2"]["Aa"])+np.asarray(stages["primary_l3"]["Aa"])
    aa_rebuilt=aa_lapse+aa_tf+np.asarray(vacuum.primary_l3_rhs(grid,g)["Aa"])+(
        np.asarray(stages["primary_l3"]["Aa"])-np.asarray(vacuum.primary_l3_rhs(grid,g)["Aa"])
    )

    ann=(r>=ANNULUS[0])&(r<ANNULUS[1])
    first=[]
    for i in CENTRE_CELLS:
        first.append({
            "i":i,"r":float(r[i]),"B":float(B_profile[i]),"source":float(src_profile[i]),
            "a":float(g.a[i]),"b":float(g.b[i]),"X":float(g.X[i]),"Aa":float(g.Aa[i]),
            "a_over_b_minus1_over_r2":float((g.a[i]/g.b[i]-1.0)/r[i]**2),
            "Lambda_over_r":float(g.Lambda[i]/r[i]),
            "R_tf_vendor":float(tf_vendor[i]),"R_tf_vendor_over_r2":float(tf_vendor[i]/r[i]**2),
            "R_tf_vendor_metric_lambda":float(tf_vendor_metric_lam[i]),
            "R_tf_metric":float(tf_metric[i]),"R_tf_metric_over_r2":float(tf_metric[i]/r[i]**2),
            "R_tf_connection_piece":float(conn_piece[i]),"R_tf_route_piece":float(route_piece[i]),
            "Aa_rhs_lapse":float(aa_lapse[i]),"Aa_rhs_tf_Ricci":float(aa_tf[i]),
            "Aa_rhs_total":float(aa_rhs[i]),"H":float(H[i]),"M":float(M[i]),"C_lambda":float(C[i]),
        })
    fixed={}
    for rr in FIXED_RADII:
        i=int(np.argmin(np.abs(r-rr)))
        fixed[str(rr)]={
            "i":i,"r":float(r[i]),"R_tf_vendor":float(tf_vendor[i]),"R_tf_metric":float(tf_metric[i]),
            "R_tf_route_piece":float(route_piece[i]),"H":float(H[i]),"M":float(M[i]),"C_lambda":float(C[i]),
        }
    return {
        "arm":label,"resolution":int(grid.n),"dr":float(grid.dr),"t":float(state.t),"tau":float(state.tau),
        "first_cells":first,"fixed_radius_nearest_cells":fixed,
        "fixed_annulus":{
            "cell_count":int(np.count_nonzero(ann)),
            "R_tf_vendor_rms":rms(tf_vendor[ann]),"R_tf_metric_rms":rms(tf_metric[ann]),
            "R_tf_route_error_rms":rms(route_piece[ann]),"R_tf_connection_piece_rms":rms(conn_piece[ann]),
            "H_rms":rms(H[ann]),"H_max_abs":maxabs(H[ann]),
            "M_rms":rms(M[ann]),"M_max_abs":maxabs(M[ann]),
            "C_lambda_max_abs":maxabs(C[ann]),
        },
        "full_domain":{
            "H_max_abs":maxabs(H),"M_max_abs":maxabs(M),"C_lambda_max_abs":maxabs(C),
            "determinant_max_abs":maxabs(np.asarray(g.a)*np.asarray(g.b)**2-1.0),
            "Aa_rhs_reconstruction_max_abs":maxabs(aa_rhs-aa_rebuilt),
        },
        "source_commit":SOURCE_COMMIT,
    }


def main():
    kernel=V55ProductionKernel()
    _,vacuum,_=adapter.vendor_modules()
    rows=[]
    for n in RESOLUTIONS:
        base=kernel.initialize(resolution=n,r_max=R_MAX,amplitude=AMPLITUDE,width=WIDTH,
             D_amplitude=D_AMPLITUDE,include_radiation=True)
        rows.append(evaluate(base,"pinned_production_initializer",np.nan*np.ones(n),np.zeros(n)))

        geom,scalars,matter,B,src=counterfactual_build_initial_data(
            base.grid,vacuum.flat_state
        )
        geom.alpha=kernel._solve_lapse(base.grid,geom,scalars,matter)[0]
        geom.beta.fill(0.0)
        geom.B.fill(0.0)
        candidate=ProductionState(
            grid=base.grid,geometry=geom,scalars=scalars,matter=matter,
            t=0.0,tau=0.0,e_folds=0.0,
        )
        rows.append(evaluate(candidate,"even_polynomial_origin_quadrature",B,src))

    comparisons=[]
    for n in RESOLUTIONS:
        a=next(x for x in rows if x["resolution"]==n and x["arm"]=="pinned_production_initializer")
        b=next(x for x in rows if x["resolution"]==n and x["arm"]=="even_polynomial_origin_quadrature")
        comparisons.append({
            "N":n,"dr":a["dr"],
            "baseline":a,"counterfactual":b,
            "delta_center_i0":{
                k:b["first_cells"][0][k]-a["first_cells"][0][k]
                for k in ("R_tf_vendor","R_tf_metric","R_tf_vendor_over_r2","R_tf_metric_over_r2",
                          "Aa_rhs_tf_Ricci","Aa_rhs_total","H","M","C_lambda")
            },
        })
    report={
        "schema":"v55_origin_half_cell_quadrature_regularity_ab_v1",
        "status":"completed","diagnostic_only":True,"source_commit":SOURCE_COMMIT,
        "quadrature_definition":{
            "baseline":"0.5*dr*source[0]",
            "counterfactual":"dr*(401/720*source[0] - 31/480*source[1] + 11/1440*source[2])",
            "derivation":"Integrate the unique quadratic interpolant in x=r^2 through the first three cell-centred source samples over [0,dr/2]; exact for source=1,r^2,r^4.",
            "all_other_updates":"Pinned six-iteration source loop and inter-centre trapezoids remain unchanged."
        },
        "resolutions":list(RESOLUTIONS),"r_max":R_MAX,"reports":rows,"paired_comparisons":comparisons,
        "guardrails":[
            "Only the origin-to-first-cell half integral differs between arms.",
            "The counterfactual is not a production edit or a fitted correction.",
            "No evolution step is taken.",
            "An improved initial regularity witness is necessary but not sufficient for scientific admission.",
        ],
    }
    out=Path("results/origin-half-cell-quadrature-regularity")
    out.mkdir(parents=True,exist_ok=True)
    (out/"report.json").write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n")
    summary=[]
    for p in comparisons:
        for arm in ("baseline","counterfactual"):
            x=p[arm]
            c=x["first_cells"][0]
            summary.append({
                "N":p["N"],"arm":arm,
                "Rtf_vendor_i0":c["R_tf_vendor"],"Rtf_metric_i0":c["R_tf_metric"],
                "Rtf_vendor_over_r2_i0":c["R_tf_vendor_over_r2"],
                "Rtf_metric_over_r2_i0":c["R_tf_metric_over_r2"],
                "Aa_rhs_total_i0":c["Aa_rhs_total"],"Aa_rhs_tf_i0":c["Aa_rhs_tf_Ricci"],
                "H_i0":c["H"],"M_i0":c["M"],"C_lambda_i0":c["C_lambda"],
                "H_annulus_max":x["fixed_annulus"]["H_max_abs"],
                "M_annulus_max":x["fixed_annulus"]["M_max_abs"],
            })
    print("ORIGIN_HALF_CELL_QUADRATURE_REGULARITY_AB "+json.dumps(summary,sort_keys=True))


if __name__=="__main__":
    main()
