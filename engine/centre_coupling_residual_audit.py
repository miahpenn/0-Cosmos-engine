"""Diagnostic-only central Hamiltonian, projection-stage, and scalar/DM exchange audit.
No production equation, default, boundary condition, or parameter is changed.
"""
import json, math, os
from pathlib import Path
import numpy as np
from engine.discrete_consistent_initial import discrete_consistent_state
from engine.hamiltonian_evolution_increment_decomposition import hamiltonian_residual
from engine.production_kernel import V55ProductionKernel, BETA_DM
from engine.stress_energy import assemble_total_stress_energy
from engine.scalar_system import cosmos_potential
from engine.v55_matter import dm_density, metric_slice_from_q
from engine.matter_system import Species, project_species
from engine.valencia import dark_matter_covector_source, spherical_metric_from_bssn
from engine import v55_pirk_adapter as adapter

SET = dict(N=160, r_max=40.0, amplitude=0.01, width=7.0, D_amplitude=1e-10,
           radiation=True, cfl=0.03, final_time=3.0,
           probes=[0.0, 0.75, 1.5, 2.25, 3.0])
NC = 5
TOL = 2.0e-10
STAGES = ["step_start_projection", "explicit_predictor_projection",
          "primary_predictor_projection", "final_explicit_block_projection",
          "completed_slice_projection"]

def first(x):
    return [float(v) for v in np.asarray(x)[:NC]]

def ricci_parts(grid, g):
    r = np.asarray(grid.centers); a, b, X = g.a, g.b, g.X
    d1, d2 = grid.cell_derivative_fourth, grid.cell_second_derivative_fourth
    ap, bp, app, bpp = d1(a, parity=1), d1(b, parity=1), d2(a, parity=1), d2(b, parity=1)
    Xp = d1(X, parity=1); chip = -0.5*Xp/X
    chipp = -0.5*(d2(X, parity=1)/X-(Xp/X)**2)
    Lp = d1(g.Lambda, parity=-1); ml = (1.0-a/b)/r**2; inv = X**2/a
    br = {
      "a_second": 0.5*app/a, "b_second": bpp/b, "Lambda_gradient": -a*Lp,
      "a_gradient_squared": -(ap/a)**2, "b_gradient_squared": 0.5*(bp/b)**2,
      "spherical_b_gradient": 2.0*(3.0-a/b)*bp/(r*b), "metric_lambda": 4.0*ml,
      "chi_second_and_square": 8.0*(chipp+chip**2),
      "chi_gradient_mixed": -8.0*chip*(0.5*ap/a-bp/b-2.0/r)
    }
    parts = {k: -inv*v for k,v in br.items()}
    _, vac, _ = adapter.vendor_modules()
    R = vac.geometry_terms(grid, g)["R"]
    err = float(np.max(np.abs(sum(parts.values(), np.zeros_like(r))-R)))
    return parts, err

def geometry_snap(grid, g):
    _, vac, _ = adapter.vendor_modules()
    c, t = vac.constraints(grid,g), vac.geometry_terms(grid,g)
    rp, re = ricci_parts(grid,g); Ab=-0.5*np.asarray(g.Aa)
    return {
      "R": first(t["R"]), "H_geometry_only": first(c["hamiltonian"]),
      "a": first(g.a), "b": first(g.b), "X": first(g.X),
      "Aa": first(g.Aa), "K": first(g.K), "Lambda": first(g.Lambda),
      "connection_constraint": first(c["connection"]),
      "determinant_constraint": first(c["determinant"]),
      "ricci_terms": {k:first(v) for k,v in rp.items()},
      "ricci_reconstruction_error": re
    }

def snap(kernel, s, label):
    grid,g,sf,matter=s.grid,s.geometry,s.scalars,s.matter
    _,vac,_=adapter.vendor_modules()
    c,t=vac.constraints(grid,g),vac.geometry_terms(grid,g)
    met=metric_slice_from_q(grid,g)
    total=assemble_total_stress_energy(grid,g,sf,matter)
    H=np.asarray(c["hamiltonian"])-16*math.pi*np.asarray(total.rho)
    inv=np.asarray(g.X)**2/np.asarray(g.a)
    srho={}
    for name,f,p,pot,scale in (
      ("S",sf.S,sf.PS,0.5*sf.S**2,1.0),
      ("D",sf.D,sf.PD,-0.5*sf.D**2,1.0),
      ("phi",sf.phi,sf.Pi,cosmos_potential(sf.phi),8*math.pi)):
        fp=grid.cell_derivative_fourth(f,parity=1)
        srho[name]=(0.5*p**2+0.5*inv*fp**2+pot)/scale
    frho={}
    for name,sp in (("dark_matter",Species.DARK_MATTER),("baryons",Species.BARYON),("radiation",Species.RADIATION)):
        frho[name]=project_species(met,getattr(matter,name),sp)["rho"]
    Ab=-0.5*np.asarray(g.Aa)
    hparts={"Ricci_R":np.asarray(t["R"]),
            "extrinsic_A_negative":-(np.asarray(g.Aa)**2+2*Ab**2),
            "trace_K_positive":(2/3)*np.asarray(g.K)**2}
    for k,v in {**srho,**frho}.items(): hparts["matter_"+k]=-16*math.pi*np.asarray(v)
    hsum=sum(hparts.values(),np.zeros_like(H))
    rhosum=sum(srho.values(),np.zeros_like(H))+sum(frho.values(),np.zeros_like(H))
    dmrho=dm_density(met,matter); srhs,_,_=kernel._rhs(s)
    phir=grid.cell_derivative_fourth(sf.phi,parity=1)
    qt,qr=dark_matter_covector_source(BETA_DM,float(dmrho[0]),float(srhs.phi[0]),float(phir[0]))
    sm=spherical_metric_from_bssn(float(grid.centers[0]),float(g.a[0]),float(g.b[0]),float(g.X[0]),float(g.alpha[0]),float(g.beta[0]))
    qe=-sm.sqrt_gamma*(qt-sm.beta*qr); qs=sm.sqrt_gamma*qr
    live=BETA_DM*float(dmrho[0]); comp=float(g.alpha[0])*live
    dPi=comp-live; direct=-2*float(sf.Pi[0])*dPi
    rp,re=ricci_parts(grid,g)
    return {
      "label":label,"t":float(s.t),"tau":float(s.tau),
      "central_H":float(H[0]),"max_abs_H_all_cells":float(np.max(np.abs(H))),
      "H_first_five":first(H),"geometry_H_first_five":first(c["hamiltonian"]),
      "H_terms_first_five":{k:first(v) for k,v in hparts.items()},
      "H_component_reconstruction_error":float(np.max(np.abs(hsum-H))),
      "rho_component_reconstruction_error":float(np.max(np.abs(rhosum-np.asarray(total.rho)))),
      "ricci_terms_first_five":{k:first(v) for k,v in rp.items()},
      "ricci_reconstruction_error":re,
      "exchange_cell0":{
        "alpha":float(g.alpha[0]),"shift":float(g.beta[0]),"Pi":float(sf.Pi[0]),
        "rho_dm":float(dmrho[0]),"beta_dm":float(BETA_DM),
        "scalar_Pi_source_live":live,"alpha_weighted_comparator_not_applied":comp,
        "comparator_minus_live":dPi,
        "relative_difference_if_lapse_weighting_is_required":abs(dPi)/max(abs(live),abs(comp),1e-300),
        "coordinate_phi_t_live":float(srhs.phi[0]),"DM_Q_t":float(qt),"DM_Q_r":float(qr),
        "DM_energy_source":float(qe),"DM_momentum_source":float(qs),
        "direct_H_rate_comparator_minus_live_source_only":float(direct),
        "direct_H_rate_live_source_only":float(-2*float(sf.Pi[0])*live),
        "direct_H_rate_comparator_source_only":float(-2*float(sf.Pi[0])*comp),
        "warning":"Comparator only, not an applied fix. Confirm sign and normalization from the frozen covariant convention."
      }
    }

def run():
    kernel=V55ProductionKernel()
    state,B,hist,sol=discrete_consistent_state(
      kernel,resolution=SET["N"],r_max=SET["r_max"],amplitude=SET["amplitude"],
      width=SET["width"],D_amplitude=SET["D_amplitude"],
      include_radiation=SET["radiation"],max_iter=12,tol=1e-13,
      max_cond=1e12,return_info=True)
    state.geometry.alpha=np.asarray(kernel._solve_lapse(state.grid,state.geometry,state.scalars,state.matter)[0]).copy()
    state.geometry.beta.fill(0.0); state.geometry.B.fill(0.0)
    report={"schema":"centre_coupling_residual_audit_v1","kind":"diagnostic_only_source_and_stage_trace",
      "source_commit":os.environ.get("GITHUB_SHA","unavailable"),"settings":SET,
      "admission":"NOT_ADMITTED_DIAGNOSTIC_ONLY","production_equations_changed":False,
      "production_defaults_changed":False,"physical_parameters_changed":False,
      "solver":{"iterations":int(sol.get("iterations",-1)),"final_max_residual":float(sol.get("final_max_residual",math.nan)),
        "max_condition_number":float(sol.get("max_condition_number",math.nan)),"B_min":float(np.min(B)),"B_max":float(np.max(B)),
        "residual_history":[float(x) for x in hist]},
      "samples":[],"probe_steps":[],"status":"in_progress","failure":None}
    orig=kernel._enforce_center_regularity
    tracing=False; trace_t=None; pidx=0; precords=[]
    def traced(grid,g):
        nonlocal pidx
        if not tracing: return orig(grid,g)
        before=geometry_snap(grid,g); fixed=orig(grid,g); after=geometry_snap(grid,fixed)
        precords.append({"stage":STAGES[pidx] if pidx<len(STAGES) else "unexpected_extra_projection",
          "index":pidx,"before":before,"after":after,
          "delta_R_center":after["R"][0]-before["R"][0],
          "delta_H_geometry_center":after["H_geometry_only"][0]-before["H_geometry_only"][0],
          "delta_connection_center":after["connection_constraint"][0]-before["connection_constraint"][0],
          "delta_determinant_center":after["determinant_constraint"][0]-before["determinant_constraint"][0]})
        pidx+=1
        return fixed
    kernel._enforce_center_regularity=traced
    dt0=SET["cfl"]*(SET["r_max"]/SET["N"])
    tol=32*np.finfo(float).eps*max(SET["final_time"],dt0)
    targets=SET["probes"]; nxt=0; steps=0
    try:
      while SET["final_time"]-state.t>tol:
        at_probe=nxt<len(targets) and abs(state.t-targets[nxt])<=max(TOL,tol)
        probe=None
        if at_probe:
          probe=snap(kernel,state,"t=%g"%targets[nxt]); report["samples"].append(probe)
          tracing=targets[nxt]<SET["final_time"]-tol; trace_t=float(state.t); pidx=0; precords=[]
          nxt+=1
        dt=min(dt0,SET["final_time"]-state.t)
        Hb=hamiltonian_residual(state)
        after=kernel.step(state,dt); Ha=hamiltonian_residual(after); steps+=1
        if at_probe:
          e=probe["exchange_cell0"]; dh=float(Ha[0]-Hb[0]); d_est=float(dt*e["direct_H_rate_comparator_minus_live_source_only"])
          report["probe_steps"].append({"t_before":float(state.t),"t_after":float(after.t),"dt":float(dt),
            "H_before":float(Hb[0]),"H_after":float(Ha[0]),"observed_delta_H":dh,
            "fixed_other_fields_estimated_delta_H_from_comparator":d_est,
            "estimated_fraction_of_observed_delta":abs(d_est)/abs(dh) if abs(dh)>1e-300 else None,
            "projections":list(precords),"projection_count":len(precords),"expected_projection_count":len(STAGES)})
        state=after; tracing=False
      if nxt<len(targets) and abs(state.t-targets[nxt])<=max(TOL,tol):
        report["samples"].append(snap(kernel,state,"t=%g"%targets[nxt])); nxt+=1
      h=hamiltonian_residual(state)
      report["accepted_steps"]=steps
      report["final"]={"t":float(state.t),"tau":float(state.tau),"central_H":float(h[0]),
        "max_abs_H_all_cells":float(np.max(np.abs(h))),"cycle_event_count":len(state.cycle.events),"handoff_count":len(state.handoffs)}
      report["status"]="completed" if abs(state.t-SET["final_time"])<=max(TOL,tol) and np.all(np.isfinite(h)) else "numerical_failure"
    except Exception as exc:
      report["status"]="numerical_failure"; report["failure"]={"t":float(state.t),"type":type(exc).__name__,"message":str(exc)}
    finally:
      kernel._enforce_center_regularity=orig
    projections=[z for p in report["probe_steps"] for z in p["projections"]]
    report["summary_checks"]={
      "expected_steps":int(round(SET["final_time"]/dt0)),"actual_steps":steps,
      "max_H_reconstruction_error":max((x["H_component_reconstruction_error"] for x in report["samples"]),default=None),
      "max_rho_reconstruction_error":max((x["rho_component_reconstruction_error"] for x in report["samples"]),default=None),
      "max_Ricci_reconstruction_error":max([x["ricci_reconstruction_error"] for x in report["samples"]]+[z["before"]["ricci_reconstruction_error"] for z in projections]+[z["after"]["ricci_reconstruction_error"] for z in projections],default=None),
      "all_sample_H_finite":all(math.isfinite(x["central_H"]) and math.isfinite(x["max_abs_H_all_cells"]) for x in report["samples"]),
      "admission":"NOT_ADMITTED_DIAGNOSTIC_ONLY"}
    out=Path("runs/centre-coupling-residual-audit"); out.mkdir(parents=True,exist_ok=True)
    (out/"report.json").write_text(json.dumps(report,indent=2,allow_nan=False))
    lines=["# Centre residual and coupling audit","","Status: %s"%report["status"],
      "Source commit: %s"%report["source_commit"],"Accepted steps: %s"%steps,
      "Final central H: %.12e"%report["final"]["central_H"],
      "Maximum H reconstruction error: %.6e"%report["summary_checks"]["max_H_reconstruction_error"],
      "Maximum Ricci-term reconstruction error: %.6e"%report["summary_checks"]["max_Ricci_reconstruction_error"],
      "Admission: NOT_ADMITTED_DIAGNOSTIC_ONLY","","No equation or parameter was changed. The lapse-weighted source is a comparator only.","",
      "| t | central H | max abs H | alpha(0) | live Pi source | alpha comparator | direct H-rate difference |",
      "|---:|---:|---:|---:|---:|---:|---:|"]
    for x in report["samples"]:
      e=x["exchange_cell0"]
      lines.append("| %.6g | %.8e | %.8e | %.8e | %.8e | %.8e | %.8e |"%(
        x["t"],x["central_H"],x["max_abs_H_all_cells"],e["alpha"],e["scalar_Pi_source_live"],
        e["alpha_weighted_comparator_not_applied"],e["direct_H_rate_comparator_minus_live_source_only"]))
    (out/"report.md").write_text("\n".join(lines)+"\n")
    print(json.dumps({"status":report["status"],"source_commit":report["source_commit"],
      "accepted_steps":steps,"final":report.get("final"),"samples":report["samples"],
      "probe_steps":report["probe_steps"],"summary_checks":report["summary_checks"],"failure":report["failure"]},indent=2,allow_nan=False))
    if report["status"]!="completed": raise SystemExit(2)
    return report

if __name__=="__main__": run()
