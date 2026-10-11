"""Fixed-domain strong-D Hamiltonian-ridge refinement; diagnostic only.

Reuses the admitted N=160, R=80 production capture and runs one N=320, R=80
production trajectory. Physics, production equations/defaults, gauge, projection,
and boundary conditions remain unchanged. Fixed CFL means both dr and dt halve.
"""
from __future__ import annotations
import json, math, pathlib, sys, traceback
import numpy as np
ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from engine import v55_pirk_adapter as adapter
from engine.production_kernel import V55ProductionKernel
from engine.stress_energy import assemble_total_stress_energy
from diagnostics.curvature_consistency.independent_metric_ricci import finite_array, metric_connection
_, vacuum, _ = adapter.vendor_modules()
N, RMAX, CFL, TFINAL, SAMPLE = 320, 80.0, 0.0075, 12.0, 0.25
AMP, WIDTH, DAMP, RADIATION = 0.01, 7.0, 1.0e-4, True
TARGETS = tuple(float(x) for x in np.arange(SAMPLE, TFINAL + SAMPLE/2, SAMPLE))
OUT = ROOT / "runs" / "strong-D-ridge-fixed-domain-refinement"
EPS = np.finfo(float).eps

def peaks(v, r, lo=4.0, hi=20.0, nmax=6):
    q = np.abs(np.asarray(v, float)); mask = (r >= lo) & (r <= hi)
    ids = np.flatnonzero(mask)
    if not len(ids): return []
    a, b = int(ids[0]), int(ids[-1]); local = [i for i in range(max(2,a), min(len(q)-2,b)+1) if q[i]>=q[i-1] and q[i]>q[i+1]]
    local.append(a + int(np.argmax(q[a:b+1]))); local = sorted(set(local), key=lambda i:q[i], reverse=True)
    return [{"radius":float(r[i]),"cell":int(i),"H":float(v[i]),"abs_H":float(q[i])} for i in local[:nmax]]

def budget(i, r, H, R, A, K, S, F):
    m = S + F; den = abs(R[i])+abs(A[i])+abs(K[i])+abs(m[i]); md=abs(S[i])+abs(F[i])
    return {"cell":int(i),"radius":float(r[i]),"H":float(H[i]),"abs_H":abs(float(H[i])),
      "R":float(R[i]),"term_A":float(A[i]),"term_K":float(K[i]),"term_scalar":float(S[i]),"term_fluid":float(F[i]),"term_matter":float(m[i]),
      "fraction_abs_R":abs(float(R[i]))/den if den else None,"fraction_abs_A":abs(float(A[i]))/den if den else None,
      "fraction_abs_K":abs(float(K[i]))/den if den else None,"fraction_abs_matter":abs(float(m[i]))/den if den else None,
      "scalar_fraction_abs_matter":abs(float(S[i]))/md if md else None,"fluid_fraction_abs_matter":abs(float(F[i]))/md if md else None,
      "cancellation_ratio_abs_H_over_piece_sum":abs(float(H[i]))/den if den else None}

def sample(state, label, requested):
    g, geo = state.grid, state.geometry; r=np.asarray(g.centers,float)
    raw=vacuum.constraints(g,geo); gt=vacuum.geometry_terms(g,geo); total=assemble_total_stress_energy(g,geo,state.scalars,state.matter)
    R=finite_array("R",np.asarray(gt["R"],float)); aa=finite_array("Aa",np.asarray(geo.Aa,float)); ab=-.5*aa; k=finite_array("K",np.asarray(geo.K,float))
    A=finite_array("term_A",-(aa**2+2*ab**2)); K=finite_array("term_K",(2/3)*k**2)
    rho_s=finite_array("rho_scalar",np.asarray(total.scalar_rho,float)); rho_f=finite_array("rho_fluid",np.asarray(total.fluid_rho,float))
    S=finite_array("term_scalar",-16*math.pi*rho_s); F=finite_array("term_fluid",-16*math.pi*rho_f); m=S+F
    H=finite_array("H_vendor",np.asarray(raw["hamiltonian"],float)+m); rebuilt=R+A+K+m
    closure=float(np.max(np.abs(H-rebuilt))); tol=2048*EPS*max(1.,float(np.max(np.abs(R))),float(np.max(np.abs(A))),float(np.max(np.abs(K))),float(np.max(np.abs(m))))
    alpha=finite_array("alpha",np.asarray(geo.alpha,float)); cl=finite_array("C_lambda",np.asarray(geo.Lambda,float)-metric_connection(g,geo))
    off=2+int(np.argmax(np.abs(H[2:]))); shell=np.flatnonzero((r>=4)&(r<=20)); rid=int(shell[np.argmax(np.abs(H[shell]))]); amin=int(np.argmin(alpha))
    return {"label":label,"requested_time":float(requested),"t":float(state.t),"tau":float(state.tau),"e_folds":float(state.e_folds),
      "max_abs_H_all":float(np.max(np.abs(H))),"max_abs_H_offcentre_cells_2plus":float(np.max(np.abs(H[2:]))),"max_abs_H_offcentre_radius":float(r[off]),
      "ridge_window":{"r_min":4.,"r_max":20.,"max_abs_H":float(abs(H[rid])),"radius":float(r[rid]),"budget":budget(rid,r,H,R,A,K,S,F),"local_peaks":peaks(H,r)},
      "top_offcentre_peaks":peaks(H,r,lo=float(r[2]),hi=float(r[-2])),
      "central":{"H":float(H[0]),"R":float(R[0]),"term_A":float(A[0]),"term_K":float(K[0]),"term_scalar":float(S[0]),"term_fluid":float(F[0]),"alpha":float(alpha[0]),"D":float(state.scalars.D[0]),"S":float(state.scalars.S[0])},
      "lapse":{"alpha_center":float(alpha[0]),"alpha_min":float(alpha[amin]),"alpha_min_radius":float(r[amin])},
      "constraints":{"max_abs_C_lambda":float(np.max(np.abs(cl))),"term_reconstruction_closure":closure,"tolerance":tol,"pass":bool(closure<=tol)},
      "profiles":{"radius":r.tolist(),"H_vendor":H.tolist(),"R_vendor":R.tolist(),"term_A":A.tolist(),"term_K":K.tolist(),"term_scalar":S.tolist(),"term_fluid":F.tolist(),"term_matter":m.tolist(),"rho_scalar":rho_s.tolist(),"rho_fluid":rho_f.tolist(),"alpha":alpha.tolist(),"D":np.asarray(state.scalars.D,float).tolist(),"S":np.asarray(state.scalars.S,float).tolist(),"PD":np.asarray(state.scalars.PD,float).tolist(),"PS":np.asarray(state.scalars.PS,float).tolist()}}

def load_reference(path):
    rows=[json.loads(x) for x in path.read_text().splitlines() if x.strip()]; rows.sort(key=lambda x:float(x["t"]))
    if len(rows)!=49: raise ValueError(f"expected 49 reference slices, got {len(rows)}")
    d=rows[0].get("domain",{})
    if (int(d.get("N",-1)),float(d.get("r_max",-1)))!=(160,80.0): raise ValueError(f"wrong reference domain: {d}")
    if any(len(x.get("H_vendor",[]))!=160 or len(x.get("radius",[]))!=160 for x in rows): raise ValueError("reference profile arrays have unexpected length")
    return rows

def interpolate_history(rows,t,field):
    ts=np.asarray([float(x["t"]) for x in rows]); j=int(np.searchsorted(ts,t,side="left"))
    if j<=0: lo=hi=0
    elif j>=len(rows): lo=hi=len(rows)-1
    else: lo,hi=j-1,j
    a,b=rows[lo],rows[hi]; va=np.asarray(a[field],float); vb=np.asarray(b[field],float)
    w=0. if hi==lo or float(b["t"])==float(a["t"]) else (t-float(a["t"]))/(float(b["t"])-float(a["t"]))
    return np.asarray(a["radius"],float),(1-w)*va+w*vb,[float(a["t"]),float(b["t"])]

def compare(records, ref):
    out=[]
    for x in records:
        p=x["profiles"]; fr=np.asarray(p["radius"],float); tr=float(x["t"])
        cr,ch,tbr=interpolate_history(ref,tr,"H_vendor"); _,ca,_=interpolate_history(ref,tr,"alpha_profile")
        valid=(fr>=cr[0]-1e-12)&(fr<=cr[-1]+1e-12); r=fr[valid]; fh=np.asarray(p["H_vendor"],float)[valid]; falpha=np.asarray(p["alpha"],float)[valid]
        ch_interp=np.interp(r,cr,ch); ca_interp=np.interp(r,cr,ca); delta=fh-ch_interp
        sf=(r>=4)&(r<=20); sc=(cr>=4)&(cr<=20); off=(r>=1.25)
        fi=int(np.flatnonzero(sf)[np.argmax(np.abs(fh[sf]))]); ci=int(np.flatnonzero(sc)[np.argmax(np.abs(ch[sc]))])
        ratio=float(abs(fh[fi])/abs(ch[ci])) if abs(ch[ci])>0 else None
        shell_delta=delta[sf]; shell_ref=ch_interp[sf]; off_delta=delta[off]; off_ref=ch_interp[off]
        out.append({"label":x["label"],"t":tr,"reference_time_bracket":tbr,"matched_domain_r_max":80.0,
          "N_reference":160,"N_refined":320,"dr_reference":.5,"dr_refined":.25,"CFL_both":CFL,
          "max_abs_profile_delta_H_including_center":float(np.max(np.abs(delta))),
          "offcentre_max_abs_profile_delta_H_r_ge_1p25":float(np.max(np.abs(off_delta))),
          "offcentre_relative_l2_profile_delta_H_r_ge_1p25":float(np.linalg.norm(off_delta)/max(np.linalg.norm(off_ref),np.finfo(float).tiny)),
          "shell_max_abs_profile_delta_H_r4_to20":float(np.max(np.abs(shell_delta))),
          "shell_relative_l2_profile_delta_H_r4_to20":float(np.linalg.norm(shell_delta)/max(np.linalg.norm(shell_ref),np.finfo(float).tiny)),
          "offcentre_rms_profile_delta_H":float(np.sqrt(np.mean(off_delta**2))),
          "max_abs_profile_delta_alpha":float(np.max(np.abs(falpha-ca_interp))),
          "reference_ridge_window_peak":{"radius":float(cr[ci]),"abs_H":float(abs(ch[ci]))},
          "refined_ridge_window_peak":{"radius":float(r[fi]),"abs_H":float(abs(fh[fi]))},
          "ridge_peak_amplitude_ratio_refined_over_reference":ratio,
          "reference_ridge_local_peaks":peaks(ch,cr),"refined_ridge_local_peaks":peaks(fh,r),
          "alpha_center_reference":float(ca[0]),"alpha_center_refined":float(p["alpha"][0]),
          "central_H_reference":float(ch[0]),"central_H_refined":float(p["H_vendor"][0])})
    late=[float(x["ridge_peak_amplitude_ratio_refined_over_reference"]) for x in out if float(x["t"])>=8 and x["ridge_peak_amplitude_ratio_refined_over_reference"] is not None]
    if late and all(v<.9 for v in late): decision="ridge-window peaks shrink by more than 10% at all sampled t>=8; resolution-sensitive, not automatically pure discretization error"
    elif late and all(.9<=v<=1.1 for v in late): decision="late ridge-window peak amplitudes remain within 10%; if local-peak identities also match, proceed to D/S/momentum/gradient scalar-source split"
    elif late and all(v<1 for v in late): decision="late ridge-window peaks are consistently lower at N=320; quantify magnitude and peak identity before classification"
    else: decision="mixed/intermediate response or peak switching; do not force a binary conclusion; inspect local peaks and full-profile norms"
    return {"reference_run_id":38150925743,"reference_artifact_id":11706797173,"reference_zip_sha256":"8b118d1f81e01f2a5d291d160b8bac42abc6a11051d8453efd9c19d0826d6d8b","time_interpolation":"linear between actual N=160 sample times","space_interpolation":"linear onto N=320 cell centers; no extrapolation","decision_rule":"Consistent decline in late ridge amplitude and profile mismatch supports resolution sensitivity; near-stability with matched peak identities supports scalar-source split. Two resolutions do not prove continuum convergence.","late_time_peak_ratios":late,"preliminary_decision":decision,"samples":out}

def main():
    OUT.mkdir(parents=True,exist_ok=True); refpath=pathlib.Path(sys.argv[1]) if len(sys.argv)>1 else None
    rows=load_reference(refpath) if refpath else None
    prof=OUT/"rmax080_N320_production_ridge_profiles.jsonl"; prof.write_text("")
    kernel=V55ProductionKernel(); state=kernel.initialize(resolution=N,r_max=RMAX,amplitude=AMP,width=WIDTH,D_amplitude=DAMP,include_radiation=RADIATION)
    if abs(float(state.grid.dr)-RMAX/N)>1e-13: raise ValueError(f"unexpected dr={state.grid.dr}")
    dt0=CFL*float(state.grid.dr); records=[]; steps=0; failure=None
    def save(t,label):
        rec=sample(state,label,t); records.append(rec)
        with prof.open("a") as f: f.write(json.dumps(rec,sort_keys=True,allow_nan=False)+"\n"); f.flush()
        print("[RIDGE_N320] SAMPLE "+json.dumps({"label":label,"t":state.t,"steps":steps,"ridge":rec["ridge_window"],"closure":rec["constraints"]},sort_keys=True),flush=True)
    save(0.0,"t=0"); j=0
    while state.t<TFINAL-1e-13:
        dt=min(dt0,TFINAL-state.t)
        try:
            nxt=kernel.step(state,dt); nxt.geometry.assert_finite_positive()
            for name,arr in (("S",nxt.scalars.S),("D",nxt.scalars.D),("PS",nxt.scalars.PS),("PD",nxt.scalars.PD),("a",nxt.geometry.a),("b",nxt.geometry.b),("X",nxt.geometry.X),("alpha",nxt.geometry.alpha),("Lambda",nxt.geometry.Lambda)): finite_array(name,np.asarray(arr,float))
            state=nxt; steps+=1
        except Exception as e:
            failure={"stage":"evolution_step","attempt":steps+1,"t_before":float(state.t),"type":type(e).__name__,"error":str(e),"traceback":traceback.format_exc(limit=12)}; break
        if steps%256==0: print("[RIDGE_N320] HEARTBEAT "+json.dumps({"t":state.t,"steps":steps,"slices":len(records)}),flush=True)
        while j<len(TARGETS) and (state.t+.5*dt>=TARGETS[j] or abs(state.t-TARGETS[j])<1e-12): save(TARGETS[j],f"target={TARGETS[j]:g}"); j+=1
    status="completed" if failure is None and state.t>=TFINAL-1e-10 and len(records)==49 else "numerical_or_sampling_failure"
    report={"schema":"strong_D_ridge_fixed_domain_refinement_v1","status":status,"diagnostic_only":True,
      "configuration":{"N":N,"r_max":RMAX,"dr":RMAX/N,"CFL":CFL,"nominal_dt":dt0,"final_time":TFINAL,"sample_interval":SAMPLE,"amplitude":AMP,"width":WIDTH,"D_amplitude":DAMP,"radiation":RADIATION,"production_physics_changed":False,"production_defaults_changed":False,"gauge_changed":False,"projection_changed":False,"boundary_conditions_changed":False},
      "accepted_steps":steps,"sample_count":len(records),"closure_all_samples_pass":all(x["constraints"]["pass"] for x in records),"failure":failure,
      "samples":[{k:x[k] for k in ("label","t","tau","max_abs_H_all","max_abs_H_offcentre_cells_2plus","max_abs_H_offcentre_radius","ridge_window","central","lapse","constraints")} for x in records],
      "interpretation_guardrail":"Fixed-domain fixed-CFL refinement halves both dx and dt. The completed half-CFL control had <0.2% effects at late sampled targets, but two spatial resolutions do not prove continuum convergence; peak switching and lapse history must be inspected."}
    if rows is not None: report["reference_comparison"]=compare(records,rows)
    (OUT/"report.json").write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n")
    md=["# Strong-D Hamiltonian ridge — fixed-domain N=160/320 comparison","",f"Status: {status}; diagnostic only.","",f"Configuration: N=320, r_max=80, dr=0.25, CFL=0.0075, D=1e-4, radiation on, t=12; {steps} accepted steps and {len(records)} slices.","","| t | N=160 peak | N=320 peak | ratio 320/160 | shell max ΔH |","|---:|---:|---:|---:|---:|"]
    if rows is not None:
        for c in report["reference_comparison"]["samples"]: md.append(f"| {c['t']:.5f} | {c['reference_ridge_window_peak']['abs_H']:.6e} @ {c['reference_ridge_window_peak']['radius']:.3f} | {c['refined_ridge_window_peak']['abs_H']:.6e} @ {c['refined_ridge_window_peak']['radius']:.3f} | {c['ridge_peak_amplitude_ratio_refined_over_reference']:.4f} | {c['shell_max_abs_profile_delta_H_r4_to20']:.6e} |")
        md += ["", "**Preliminary decision:** "+report["reference_comparison"]["preliminary_decision"]]
    md += ["",f"Hamiltonian reconstruction closure passed every slice: **{report['closure_all_samples_pass']}**.","",report["interpretation_guardrail"],""]
    (OUT/"report.md").write_text("\n".join(md))
    print("[RIDGE_N320] FINAL "+json.dumps({"status":status,"steps":steps,"slices":len(records),"closure_all_samples_pass":report["closure_all_samples_pass"],"decision":report.get("reference_comparison",{}).get("preliminary_decision")},sort_keys=True),flush=True)
    if status!="completed": raise RuntimeError("N=320 ridge refinement did not complete; checkpoint profiles preserved")
    if not report["closure_all_samples_pass"]: raise RuntimeError("Hamiltonian term reconstruction failed")

if __name__=="__main__": main()
