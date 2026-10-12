#!/usr/bin/env python3
"""Operator-level centre audit using frozen production states only.

1. Verify pinned fourth-order cell derivatives against manufactured regular
   even/odd polynomial profiles at the first cells.
2. Compare the production momentum/connection witnesses with an independent
   local even-in-r polynomial reconstruction over a fixed physical centre
   window. This is a sensitivity probe only, never an evolution repair.

States: production initial slice and saved t=2.16 checkpoints at N=40,60,80,120.
No evolution steps are taken. The independent fit changes no production field
or source; its output is kept in a separate diagnostic profile.
"""
from __future__ import annotations
import json
from pathlib import Path
import sys
import numpy as np

sys.path.insert(0, str(Path("source").resolve()))
from engine.production_kernel import V55ProductionKernel
from engine.scalar_system import ScalarFields, scalar_projection
from engine.matter_system import ConservedSpecies, Species, project_species
from engine.v55_matter import V55MatterState, metric_slice_from_q
from engine.stress_energy import assemble_total_stress_energy
from engine import v55_pirk_adapter as adapter

_, vacuum, _ = adapter.vendor_modules()
RESOLUTIONS=(40,60,80,120)
R_MAX=80.0
TARGET_T=2.16
SPATIAL_RUN_ID=38212879394
SOURCE_COMMIT="6d8b3f2790644591326df57c55908863f994cb0c"
SCALAR_NAMES=("S","PS","D","PD","phi","Pi")
SPECIES_PREFIX={"dark_matter":"dm","baryons":"baryon","radiation":"radiation"}
CONSERVE=("rest","energy_t","momentum_r")
GEOMETRY_NAMES=("a","b","X","alpha","beta","Aa","K","Lambda","B")
FIT_WINDOW_R=6.0
REPORT_WINDOW_R=5.0
SAMPLE_RADII=(2.0,3.0,4.0,5.0)


def d1(grid, arr, parity=+1):
    return np.asarray(grid.cell_derivative_fourth(arr,parity=parity),dtype=float)


def d2(grid, arr, parity=+1):
    return np.asarray(grid.cell_second_derivative_fourth(arr,parity=parity),dtype=float)


def peak(arr,r,mask):
    ids=np.flatnonzero(mask)
    if not ids.size: return None
    i=int(ids[np.argmax(np.abs(arr[ids]))])
    return {"i":i,"r":float(r[i]),"signed":float(arr[i]),"abs":float(abs(arr[i]))}


def poly_even_derivative(r, values, fit_mask):
    """Derivative of least-squares f(r)=c0+c2*r^2+c4*r^4 on fixed r<=6 window."""
    r=np.asarray(r,dtype=float); values=np.asarray(values,dtype=float)
    idx=np.flatnonzero(fit_mask)
    if idx.size<3:
        raise ValueError(f"need >=3 nodes in centre fit window, got {idx.size}")
    x=r[idx]**2
    A=np.column_stack([np.ones_like(x),x,x*x])
    coeff,_,rank,_=np.linalg.lstsq(A,values[idx],rcond=None)
    if rank<3: raise ValueError("rank-deficient even polynomial reconstruction")
    deriv=2.0*coeff[1]*r+4.0*coeff[2]*r**3
    return deriv,{"fit_cell_indices":idx.tolist(),"fit_node_count":int(idx.size),
                  "coefficients_c0_c2_c4":[float(v) for v in coeff],
                  "max_fit_residual":float(np.max(np.abs(A@coeff-values[idx])))}


def manufactured_test(grid):
    r=np.asarray(grid.centers,dtype=float)
    even=1.0+0.2*r**2-0.001*r**4
    even_d=0.4*r-0.004*r**3
    even_d2=0.4-0.012*r**2
    odd=0.3*r+0.01*r**3
    odd_d=0.3+0.03*r**2
    odd_d2=0.06*r
    ix=np.arange(min(5,len(r)))
    rows={}
    for name,got,expected in (
        ("even_d1",d1(grid,even,+1),even_d),
        ("even_d2",d2(grid,even,+1),even_d2),
        ("odd_d1",d1(grid,odd,-1),odd_d),
        ("odd_d2",d2(grid,odd,-1),odd_d2),
    ):
        err=np.asarray(got)-np.asarray(expected)
        rows[name]={"max_abs_error_first_five":float(np.max(np.abs(err[ix]))),
                    "signed_errors_first_five":[float(v) for v in err[ix]],
                    "tolerance":1e-10,"pass":bool(np.max(np.abs(err[ix]))<1e-10)}
    return rows


def restore(path,kernel):
    with np.load(path,allow_pickle=False) as z: data={k:z[k].copy() for k in z.files}
    n=int(data["r"].size); dr=float(data["r"][1]-data["r"][0])
    state=kernel.initialize(resolution=n,r_max=n*dr,amplitude=.01,width=7.0,D_amplitude=1.e-10,include_radiation=True)
    if not np.allclose(state.grid.centers,data["r"],rtol=0.0,atol=1.e-12*max(1.0,n*dr)):
        raise ValueError(f"{path}: grid mismatch")
    for name in GEOMETRY_NAMES: setattr(state.geometry,name,np.asarray(data[name],dtype=float).copy())
    state.scalars=ScalarFields(*(np.asarray(data[name],dtype=float).copy() for name in SCALAR_NAMES))
    species={}
    for name,prefix in SPECIES_PREFIX.items():
        species[name]=ConservedSpecies(*(np.asarray(data[f"{prefix}_{field}"],dtype=float).copy() for field in CONSERVE))
    state.matter=V55MatterState(dark_matter=species["dark_matter"],baryons=species["baryons"],radiation=species["radiation"])
    state.t=float(data["t"]); state.tau=float(data["tau"]); state.e_folds=float(data["e_folds"])
    state.history=[]; state.handoffs=[]
    return state


def audit_state(state,label,scheduled_t=None):
    grid=state.grid; r=np.asarray(grid.centers,dtype=float); g=state.geometry; s=state.scalars; matter=state.matter
    total=assemble_total_stress_energy(grid,g,s,matter)
    raw=vacuum.constraints(grid,g)
    base_geom=np.asarray(raw["momentum"],dtype=float)
    base_source=8.0*np.pi*np.asarray(total.j,dtype=float)
    base_M=base_geom-base_source
    base_C=np.asarray(g.Lambda,dtype=float)-(
        d1(grid,g.a,+1)/(2*np.asarray(g.a)**2)
        -d1(grid,g.b,+1)/(np.asarray(g.a)*np.asarray(g.b))
        +2.0/r*(1.0/np.asarray(g.b)-1.0/np.asarray(g.a))
    )
    h_base=np.asarray(raw["hamiltonian"],dtype=float)-16*np.pi*np.asarray(total.rho,dtype=float)

    fit_mask=r<=FIT_WINDOW_R
    report_mask=r<=REPORT_WINDOW_R
    derivatives={}
    fits={}
    for name in ("a","b","X","Aa","K","S","D","phi"):
        values=np.asarray(getattr(g,name),dtype=float) if name in GEOMETRY_NAMES else np.asarray(getattr(s,name),dtype=float)
        deriv,fit=poly_even_derivative(r,values,fit_mask)
        derivatives[name]=deriv
        fits[name]=fit

    Aa=np.asarray(g.Aa,dtype=float); Ab=-0.5*Aa
    Mgeom_poly=(
        derivatives["Aa"]-(2.0/3.0)*derivatives["K"]
        +6.0*Aa*(-0.5*derivatives["X"]/np.asarray(g.X,dtype=float))
        +(Aa-Ab)*(2.0/r+derivatives["b"]/np.asarray(g.b,dtype=float))
    )

    # Same conservative fluids, but independently differentiate the three scalar
    # fields in the local current projection.
    jS_poly=-np.asarray(s.PS,dtype=float)*derivatives["S"]
    jD_poly=-np.asarray(s.PD,dtype=float)*derivatives["D"]
    jphi_poly=-np.asarray(s.Pi,dtype=float)*derivatives["phi"]/(8.0*np.pi)
    metric=metric_slice_from_q(grid,g)
    jfluid=sum((np.asarray(project_species(metric,getattr(matter,name),species)["j"],dtype=float)
                for name,species in (("dark_matter",Species.DARK_MATTER),
                                     ("baryons",Species.BARYON),
                                     ("radiation",Species.RADIATION))),np.zeros_like(r))
    Msource_poly=8.0*np.pi*(jS_poly+jD_poly+jphi_poly+jfluid)
    M_poly=Mgeom_poly-Msource_poly

    C_poly=np.asarray(g.Lambda,dtype=float)-(
        derivatives["a"]/(2.0*np.asarray(g.a,dtype=float)**2)
        -derivatives["b"]/(np.asarray(g.a,dtype=float)*np.asarray(g.b,dtype=float))
        +2.0/r*(1.0/np.asarray(g.b,dtype=float)-1.0/np.asarray(g.a,dtype=float))
    )

    sample={}
    for rr in SAMPLE_RADII:
        i=int(np.argmin(np.abs(r-rr)))
        sample[str(rr)]={"target_r":rr,"grid_i":i,"grid_r":float(r[i]),
            "H_production":float(h_base[i]),
            "M_geom_production":float(base_geom[i]),"M_source_production":float(base_source[i]),
            "M_residual_production":float(base_M[i]),"C_lambda_production":float(base_C[i]),
            "M_geom_poly":float(Mgeom_poly[i]),"M_source_poly":float(Msource_poly[i]),
            "M_residual_poly":float(M_poly[i]),"C_lambda_poly":float(C_poly[i]),
            "delta_M_residual_poly_minus_prod":float(M_poly[i]-base_M[i]),
            "delta_C_lambda_poly_minus_prod":float(C_poly[i]-base_C[i])}
    first_cells={}
    for i in range(min(6,len(r))):
        first_cells[str(i)]={"r":float(r[i]),"M_geom_production":float(base_geom[i]),
            "M_source_production":float(base_source[i]),"M_residual_production":float(base_M[i]),
            "M_geom_poly":float(Mgeom_poly[i]),"M_source_poly":float(Msource_poly[i]),
            "M_residual_poly":float(M_poly[i]),"delta_M":float(M_poly[i]-base_M[i]),
            "C_lambda_production":float(base_C[i]),"C_lambda_poly":float(C_poly[i]),
            "delta_C":float(C_poly[i]-base_C[i]),"H_production":float(h_base[i])}

    mask_summaries={}
    for name,mask in (("r_le_5",report_mask),("full_domain",np.ones_like(r,dtype=bool))):
        mask_summaries[name]={
            "M_production_peak":peak(base_M,r,mask),
            "M_polynomial_probe_peak":peak(M_poly,r,mask),
            "max_abs_delta_M":float(np.max(np.abs((M_poly-base_M)[mask]))),
            "C_production_peak":peak(base_C,r,mask),
            "C_polynomial_probe_peak":peak(C_poly,r,mask),
            "max_abs_delta_C":float(np.max(np.abs((C_poly-base_C)[mask]))),
        }

    return {
        "label":label,"resolution":int(len(r)),"actual_t":float(state.t),"scheduled_t":scheduled_t,
        "r_max":float(len(r)*grid.dr),"dr":float(grid.dr),"fit_window_radius":FIT_WINDOW_R,
        "fit_node_count":int(np.count_nonzero(fit_mask)),"fit":fits,
        "first_cells":first_cells,"fixed_radius_samples":sample,"regions":mask_summaries,
        "max_fit_residual_over_fields":max(f["max_fit_residual"] for f in fits.values()),
        "interpretation_limit":"The polynomial derivative is a frozen-state sensitivity probe only; it does not change production fields, evolution, or stress sources.",
    }


def main():
    kernel=V55ProductionKernel()
    all_test={str(n):manufactured_test(kernel.initialize(resolution=n,r_max=R_MAX,include_radiation=True).grid) for n in RESOLUTIONS}
    cases=[]
    root=Path("spatial-artifact/runs/lambda-spatial-refinement")
    for n in RESOLUTIONS:
        initial=kernel.initialize(resolution=n,r_max=R_MAX,amplitude=.01,width=7.0,D_amplitude=1.e-10,include_radiation=True)
        cases.append(audit_state(initial,f"N{n}_INITIAL",0.0))
        path=root/f"N{n}"/f"checkpoint_t{TARGET_T:012.6f}.npz"
        if not path.is_file(): raise FileNotFoundError(path)
        saved=restore(path,kernel)
        cases.append(audit_state(saved,f"N{n}_T2P16",TARGET_T))
    report={
      "schema":"v55_center_operator_parity_audit_v1","status":"completed","diagnostic_only":True,
      "source_commit":SOURCE_COMMIT,"spatial_refinement_run_id":SPATIAL_RUN_ID,
      "resolutions":list(RESOLUTIONS),"state_count":len(cases),
      "manufactured_regular_polynomial_tests":all_test,"states":cases,
      "guardrails":[
        "Manufactured polynomial tests check the pinned operator's parity implementation only.",
        "Local even-in-r polynomial derivatives are an alternate diagnostic reconstruction, not a proposed production stencil.",
        "All evolution fields and sources remain unchanged; no trajectory is advanced.",
        "A change in residual under the alternate derivative measures stencil sensitivity, not physical causality."
      ]
    }
    out=Path("results/center-operator-parity-audit");out.mkdir(parents=True,exist_ok=True)
    (out/"report.json").write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("CENTER_OPERATOR_PARITY_AUDIT "+json.dumps({
      "status":report["status"],
      "manufactured_tests_pass":all(t["pass"] for byn in all_test.values() for t in byn.values()),
      "max_manufactured_error":max(v["max_abs_error_first_five"] for byn in all_test.values() for v in byn.values()),
      "case_count":len(cases),
      "rows":[{"label":q["label"],"N":q["resolution"],"t":q["actual_t"],
               "fit_nodes":q["fit_node_count"],"max_abs_delta_M_rle5":q["regions"]["r_le_5"]["max_abs_delta_M"],
               "samples":q["fixed_radius_samples"]}
              for q in cases]
    },sort_keys=True))


if __name__=="__main__":
    main()
