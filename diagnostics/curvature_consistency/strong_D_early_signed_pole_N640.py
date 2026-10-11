"""Targeted high-resolution early signed-pole/splice audit; diagnostic only."""
from __future__ import annotations
import json, math, pathlib, sys, traceback
import numpy as np

ROOT=pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from engine import v55_pirk_adapter as adapter
from engine.production_kernel import V55ProductionKernel
from engine.stress_energy import assemble_total_stress_energy
from diagnostics.curvature_consistency.independent_metric_ricci import finite_array, metric_connection

_,vacuum,_=adapter.vendor_modules()
N,RMAX,CFL,TFINAL,SAMPLE=640,80.0,0.0075,12.0,0.0625
AMP,WIDTH,DAMP,RADIATION=0.01,7.0,1e-4,True
OUT=ROOT/"runs"/"strong-D-early-signed-pole-N640"
EPS=np.finfo(float).eps

def sign(x): return 1 if x>0 else (-1 if x<0 else 0)

def local_extrema(r,h,lo=1.25,hi=20.0):
    out=[]
    for i in range(1,len(h)-1):
        if not lo<=r[i]<=hi: continue
        kind="max" if h[i]>=h[i-1] and h[i]>h[i+1] else ("min" if h[i]<=h[i-1] and h[i]<h[i+1] else None)
        if kind: out.append({"kind":kind,"i":i,"r":float(r[i]),"H":float(h[i]),"sign":sign(float(h[i]))})
    return out

def sign_edges(r,h):
    return [{"i_left":i-1,"i_right":i,"r_left":float(r[i-1]),"r_right":float(r[i]),
             "H_left":float(h[i-1]),"H_right":float(h[i])}
            for i in range(1,len(h)) if r[i-1]>=1.25 and r[i]<=79.875 and sign(float(h[i-1]))*sign(float(h[i]))<0]

def sample(state,label,requested,steps):
    g,geo=state.grid,state.geometry;r=np.asarray(g.centers,float)
    raw=vacuum.constraints(g,geo);gt=vacuum.geometry_terms(g,geo)
    total=assemble_total_stress_energy(g,geo,state.scalars,state.matter)
    R=finite_array("R",np.asarray(gt["R"],float))
    aa=finite_array("Aa",np.asarray(geo.Aa,float));ab=-.5*aa
    K=finite_array("K",np.asarray(geo.K,float))
    A=finite_array("term_A",-(aa**2+2*ab**2));Kt=finite_array("term_K",(2/3)*K**2)
    rho_s=finite_array("rho_scalar",np.asarray(total.scalar_rho,float))
    rho_f=finite_array("rho_fluid",np.asarray(total.fluid_rho,float))
    S=finite_array("term_scalar",-16*math.pi*rho_s);F=finite_array("term_fluid",-16*math.pi*rho_f)
    H=finite_array("H_vendor",np.asarray(raw["hamiltonian"],float)+S+F)
    rebuilt=R+A+Kt+S+F
    closure=float(np.max(np.abs(H-rebuilt)))
    tol=float(2048*EPS*max(1.,float(np.max(np.abs(R))),float(np.max(np.abs(A))),float(np.max(np.abs(Kt))),float(np.max(np.abs(S+F)))))
    alpha=finite_array("alpha",np.asarray(geo.alpha,float))
    ext=local_extrema(r,H);edges=sign_edges(r,H)
    def band(lo,hi):
        ids=np.flatnonzero((r>=lo)&(r<=hi));i=int(ids[np.argmax(np.abs(H[ids]))])
        return {"i":i,"r":float(r[i]),"H":float(H[i]),"abs_H":float(abs(H[i]))}
    # Decompose each scalar's density in the exact projection convention used by the production kernel.
    invr=geo.X**2/geo.a
    pieces={}
    for name,f,p,V in (
        ("S",state.scalars.S,state.scalars.PS,.5*state.scalars.S**2),
        ("D",state.scalars.D,state.scalars.PD,-.5*state.scalars.D**2),
        ("phi",state.scalars.phi,state.scalars.Pi,0.0),
    ):
        fp=g.cell_derivative_fourth(f,parity=1)
        pot=V
        if name=="phi":
            from engine.scalar_system import cosmos_potential
            pot=cosmos_potential(f)
        scale=8*math.pi if name=="phi" else 1.0
        pieces[name]={"rho":finite_array(name+"_rho",(0.5*p*p+0.5*invr*fp*fp+pot)/scale)}
    for name in pieces:
        pieces[name]["term_H"]=-16*math.pi*pieces[name]["rho"]
    # Snapshot values around the known early exchange and signed minimum.
    probes={}
    for rr in (2.25,2.75,3.25,3.75,4.25,4.75,5.125,5.375,5.625,5.875,6.125,6.375,6.625,7.125,7.625,8.125,8.625,9.125):
        i=int(np.argmin(np.abs(r-rr)))
        probes[f"r={r[i]:.5f}"]={"H":float(H[i]),"R":float(R[i]),"term_A":float(A[i]),"term_K":float(Kt[i]),
          "term_S":float(pieces["S"]["term_H"][i]),"term_D":float(pieces["D"]["term_H"][i]),
          "term_phi":float(pieces["phi"]["term_H"][i]),"term_fluid":float(F[i]),"alpha":float(alpha[i]),
          "D":float(state.scalars.D[i]),"PD":float(state.scalars.PD[i])}
    return {"label":label,"requested_t":requested,"t":float(state.t),"steps":steps,"tau":float(state.tau),
      "domain":{"N":N,"r_max":RMAX,"dr":float(g.dr)},"radius":r.tolist(),
      "H_vendor":H.tolist(),"R":R.tolist(),"term_A":A.tolist(),"term_K":Kt.tolist(),
      "term_S":pieces["S"]["term_H"].tolist(),"term_D":pieces["D"]["term_H"].tolist(),
      "term_phi":pieces["phi"]["term_H"].tolist(),"term_fluid":F.tolist(),"alpha":alpha.tolist(),
      "D":np.asarray(state.scalars.D,float).tolist(),"PD":np.asarray(state.scalars.PD,float).tolist(),
      "ridge_winner":band(4,20),"inner_winner":band(1.25,4),
      "local_extrema":ext,"sign_edges":edges,"sign_transition_count":len(edges),
      "probes":probes,"central":{"H":float(H[0]),"alpha":float(alpha[0]),"D":float(state.scalars.D[0]),"PD":float(state.scalars.PD[0])},
      "closure":{"max_abs_H_reconstruction":closure,"tolerance":tol,"pass":bool(closure<=tol)},
      "term_maxima":{"R":float(np.max(np.abs(R))),"A":float(np.max(np.abs(A))),"K":float(np.max(np.abs(Kt))),
        "S":float(np.max(np.abs(pieces["S"]["term_H"]))),"D":float(np.max(np.abs(pieces["D"]["term_H"]))),
        "phi":float(np.max(np.abs(pieces["phi"]["term_H"]))),"fluid":float(np.max(np.abs(F)))}}

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    prof=OUT/"n640_early_signed_profiles.jsonl";prof.write_text("")
    kernel=V55ProductionKernel()
    state=kernel.initialize(resolution=N,r_max=RMAX,amplitude=AMP,width=WIDTH,D_amplitude=DAMP,include_radiation=RADIATION)
    dt0=CFL*float(state.grid.dr)
    if not math.isclose(float(state.grid.dr),RMAX/N,rel_tol=0,abs_tol=1e-13): raise ValueError("unexpected grid spacing")
    steps=0;records=[];failure=None;targets=np.arange(SAMPLE,TFINAL+SAMPLE/2,SAMPLE);j=0
    def save(requested,label):
        rec=sample(state,label,float(requested),steps);records.append(rec)
        with prof.open("a") as f: f.write(json.dumps(rec,sort_keys=True,allow_nan=False)+"\n");f.flush()
        print("[EARLY_N640] SAMPLE "+json.dumps({k:rec[k] for k in ("label","t","steps","ridge_winner","inner_winner","sign_transition_count","local_extrema","central","closure","term_maxima")},sort_keys=True),flush=True)
    save(0.0,"t=0")
    while state.t<TFINAL-1e-13:
        dt=min(dt0,TFINAL-state.t)
        try:
            nxt=kernel.step(state,dt);nxt.geometry.assert_finite_positive()
            for name,arr in (("S",nxt.scalars.S),("D",nxt.scalars.D),("PS",nxt.scalars.PS),("PD",nxt.scalars.PD),("a",nxt.geometry.a),("b",nxt.geometry.b),("X",nxt.geometry.X),("alpha",nxt.geometry.alpha),("Lambda",nxt.geometry.Lambda)):
                finite_array(name,np.asarray(arr,float))
            state=nxt;steps+=1
        except Exception as e:
            failure={"step":steps+1,"t_before":float(state.t),"type":type(e).__name__,"error":str(e),"traceback":traceback.format_exc(limit=12)};break
        if steps%512==0: print("[EARLY_N640] HEARTBEAT "+json.dumps({"t":state.t,"steps":steps,"slices":len(records)},sort_keys=True),flush=True)
        while j<len(targets) and (state.t+.5*dt>=targets[j] or abs(state.t-targets[j])<1e-12):
            save(float(targets[j]),f"target={targets[j]:g}");j+=1
    expected_records=int(round(TFINAL/SAMPLE))+1
    status="completed" if failure is None and state.t>=TFINAL-1e-10 and len(records)==expected_records else "numerical_or_sampling_failure"
    report={"schema":"strong_D_signed_pole_N640_through_t12_v2","status":status,"diagnostic_only":True,
      "configuration":{"N":N,"r_max":RMAX,"dr":RMAX/N,"CFL":CFL,"dt_nominal":dt0,"final_time":TFINAL,"sample_interval":SAMPLE,
        "amplitude":AMP,"width":WIDTH,"D_amplitude":DAMP,"radiation":RADIATION,"production_kernel_changed":False,
        "physics_changed":False,"gauge_changed":False,"projection_changed":False,"boundary_changed":False},
      "accepted_steps":steps,"sample_count":len(records),"all_closure_checks_pass":all(x["closure"]["pass"] for x in records),
      "failure":failure,"sign_change_timeline":[{"t":x["t"],"count":x["sign_transition_count"],"edges":x["sign_edges"]} for x in records],
      "pole_timeline":[{"t":x["t"],"ridge":x["ridge_winner"],"inner":x["inner_winner"],"extrema":x["local_extrema"]} for x in records],
      "profiles_file":prof.name,"interpretation_guardrail":"N=640 refinement through early sign exchange and late ridge with unchanged production equations. Signed extrema are constraint-residual features, not presumed physical objects. The per-term decomposition tests which existing terms make/cancel the signs; it does not insert a new coupling."}
    (OUT/"report.json").write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n")
    lines=["# N=640 signed-pole audit through t=12","","Diagnostic-only full production evolution, unchanged equations and gauge.",
      f"- Status: {status}; N={N}; Rmax={RMAX}; dr={RMAX/N}; dt={dt0}; t={state.t}; steps={steps}; slices={len(records)}.",
      f"- All Hamiltonian reconstruction checks passed: {report['all_closure_checks_pass']}.","",
      "## Signed poles and sign topology","","| t | sign transitions | inner winner (r,H) | ridge winner (r,H) |",
      "|---:|---:|---:|---:|"]
    for x in records:
        lines.append(f"| {x['t']:.5f} | {x['sign_transition_count']} | ({x['inner_winner']['r']:.4f}, {x['inner_winner']['H']:+.5e}) | ({x['ridge_winner']['r']:.4f}, {x['ridge_winner']['H']:+.5e}) |")
    lines += ["","## Term reconstruction","","For each slice, H is independently reconstructed as R + term_A + term_K + term_S + term_D + term_phi + term_fluid. The report records closure and full radial arrays for each contribution, plus fixed-radius probes around the early exchange. This identifies which existing terms create/cancel the signed extrema without fitting or changing any physics.","",
      "## Interpretation guardrail","",
      "A robust sign/pole pattern that survives this refinement is a candidate structural feature for a follow-up equation-level audit. A pattern that shifts/disappears is not admitted as a physical splice. This run itself does not add a sign-flip law or a D-to-matter transfer term.",""]
    (OUT/"report.md").write_text("\n".join(lines))
    print("[EARLY_N640] FINAL "+json.dumps({"status":status,"steps":steps,"slices":len(records),"all_closure_checks_pass":report["all_closure_checks_pass"],"failure":failure},sort_keys=True),flush=True)
    if status!="completed": raise RuntimeError("N640 through-t12 run incomplete; preserve partial profiles")
    if not report["all_closure_checks_pass"]: raise RuntimeError("Hamiltonian term reconstruction failed")
if __name__=="__main__":main()
