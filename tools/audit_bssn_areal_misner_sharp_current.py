#!/usr/bin/env python3
"""Research-only transformed Misner-Sharp current identity probe.

For each saved full-continuous-profile accepted t=0.60 state, take one copied
production-kernel step at dt and dt/2. Build the areal radius R=r*sqrt(b)/X,
compute m_MS from the invariant norm of grad(R), transform T_t^R into the
(t,R) chart, and compare the fixed-R finite-difference mass rate to the
time-averaged transformed current. The saved accepted state is never mutated.
The current production branch is zero-shift; nonzero shift fails closed.
"""
from __future__ import annotations
import hashlib,json,os,platform,sys
from pathlib import Path
from typing import Any
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
os.environ.setdefault("COSMOS_PINNED_SOURCE_ROOT",str(ROOT))
from engine.v55_matter import total_matter_projection
from engine import production_kernel as pkmod
from engine.production_kernel import V55ProductionKernel
from tools import audit_full_profile_propagation_candidate as saved
from tools import audit_origin_cell_joint_shadow_one_step as origin_shadow

SOURCE="4aec88d0e3ff028e0a65a8ac36b30948625df792"
VENDOR="d6052d605673ce9d82cdc99b0df79855fa2ea215"
N=int(os.environ.get("AREAL_CURRENT_RESOLUTION","120"))
FULL_DIR=Path(os.environ.get("FULL_PROFILE_ARTIFACT_DIR",f"reference/full/N{N}"))
OUT=Path(os.environ.get("AREAL_CURRENT_OUT",f"results/areal-current-identity-N{N}-2026-10-25"))
TARGET_TIME=0.60
CFL=0.03
FOUR_PI=4.0*np.pi

def sha(path:Path)->str:
 h=hashlib.sha256()
 with path.open("rb") as f:
  for chunk in iter(lambda:f.read(1<<20),b""):h.update(chunk)
 return h.hexdigest()

def finite(name,arr):
 a=np.asarray(arr,dtype=np.float64)
 if not np.all(np.isfinite(a)):raise RuntimeError(f"{name} contains non-finite values")
 return a

def clone_state(s):
 return type(s)(
  grid=s.grid,geometry=s.geometry.copy(),scalars=s.scalars.copy(),matter=s.matter.copy(),
  electromagnetic=None,t=float(s.t),tau=float(s.tau),e_folds=float(s.e_folds),
  cycle=getattr(s,"cycle",0),history=list(getattr(s,"history",[])),
  handoffs=list(getattr(s,"handoffs",[])))

def refreshed_state(kernel,s):
 q=clone_state(s)
 q.geometry.alpha=kernel._solve_lapse_for_state(
  q.grid,q.geometry,q.scalars,q.matter,q.electromagnetic
 )[0]
 q.geometry.beta.fill(0.0)
 q.geometry.B.fill(0.0)
 if np.max(np.abs(q.geometry.beta))>1e-14:
  raise RuntimeError("areal-current probe requires the pinned zero-shift branch")
 return q

def geometry_rate_and_mass(kernel,s):
 """Return areal radius, R_t, R_r, Misner-Sharp mass and transformed flux."""
 q=refreshed_state(kernel,s)
 g=q.geometry;grid=q.grid
 r=finite("grid centres",grid.centers)
 a=finite("a",g.a);b=finite("b",g.b);X=finite("X",g.X)
 alpha=finite("alpha",g.alpha);beta=finite("beta",g.beta)
 if np.any(a<=0) or np.any(b<=0) or np.any(X<=0) or np.any(alpha<=0):
  raise RuntimeError("non-positive metric field in areal-current probe")
 terms=kernel._geometry_stage_terms(grid,g,q.scalars,q.matter,q.electromagnetic)
 explicit=terms["explicit"]
 adot=finite("a RHS",explicit["a"]);bdot=finite("b RHS",explicit["b"]);Xdot=finite("X RHS",explicit["X"])
 R=finite("areal radius",r*np.sqrt(b)/X)
 Rr=finite("areal-radius radial derivative",grid.cell_derivative_fourth(R,parity=-1))
 Rt=finite("areal-radius time derivative",R*(0.5*bdot/b-Xdot/X))
 if np.any(np.diff(R)<=0) or np.any(Rr<=0):
  raise RuntimeError("areal-radius map is not strictly monotone")
 grr=a/(X*X)
 chi_spatial=finite("spatial-gradient contribution to chi",(X*X/a)*Rr*Rr)
 chi_time=finite("time/shift contribution to chi",-((Rt-beta*Rr)/alpha)**2)
 chi=finite("invariant radial norm",chi_spatial+chi_time)
 mass=finite("Misner-Sharp mass",0.5*R*(1.0-chi))
 proj=total_matter_projection(grid,g,q.scalars,q.matter)
 E=finite("total Eulerian energy",proj["rho"])
 j=finite("total covariant radial current",proj["j"])
 pr=finite("total radial pressure",proj["pr"])
 j_up=j/grr
 # Mixed components in the pinned zero-shift (t,r) chart:
 # T_t^t=-E, T_t^r=-alpha*j^r, T_r^t=j/alpha, T_r^r=p_r.
 T_t_t=-E
 T_t_r=-alpha*j_up
 T_r_t=j/alpha
 T_r_r=pr
 # Tensor transformation to y=(t,R(t,r)), x=(t,r), split by source component.
 TtR_energy=finite("T_t^R energy-motion term",T_t_t*Rt)
 TtR_radial_current=finite("T_t^R radial-current term",T_t_r*Rr)
 TtR_mixed_current=finite("T_t^R time-index current term",-(Rt/Rr)*(T_r_t*Rt))
 TtR_pressure=finite("T_t^R radial-pressure motion term",-(Rt/Rr)*(T_r_r*Rr))
 T_t_R=finite("transformed mixed current T_t^R",
  TtR_energy+TtR_radial_current+TtR_mixed_current+TtR_pressure)
 rhs_components={
  "energy_motion":finite("current RHS energy-motion",FOUR_PI*R*R*TtR_energy),
  "radial_current":finite("current RHS radial-current",FOUR_PI*R*R*TtR_radial_current),
  "mixed_current":finite("current RHS mixed-current",FOUR_PI*R*R*TtR_mixed_current),
  "pressure_motion":finite("current RHS pressure-motion",FOUR_PI*R*R*TtR_pressure),
 }
 rhs=finite("Misner-Sharp current RHS",FOUR_PI*R*R*T_t_R)
 return {"state":q,"R":R,"Rt":Rt,"Rr":Rr,"chi":chi,"chi_spatial":chi_spatial,
         "chi_time":chi_time,"mass":mass,"T_t_R":T_t_R,"rhs":rhs,
         "rhs_components":rhs_components,"max_abs_beta":float(np.max(np.abs(beta)))}

def compare_step(kernel,base,dt):
 before=saved.tangent.state_digest(base)
 start=geometry_rate_and_mass(kernel,base)
 advanced=kernel.step(clone_state(base),dt)
 end=geometry_rate_and_mass(kernel,advanced)
 after=saved.tangent.state_digest(base)
 if before!=after:raise RuntimeError("saved accepted input state mutated by copied step")
 lo=max(float(start["R"][0]),float(end["R"][0]))
 hi=min(float(start["R"][-1]),float(end["R"][-1]),10.0)
 if not hi>lo:raise RuntimeError("no common areal-radius interval")
 Rcommon=np.linspace(lo,hi, max(300,base.grid.n*2))
 m0=np.interp(Rcommon,start["R"],start["mass"])
 m1=np.interp(Rcommon,end["R"],end["mass"])
 q0=np.interp(Rcommon,start["R"],start["rhs"])
 q1=np.interp(Rcommon,end["R"],end["rhs"])
 chi0=np.interp(Rcommon,start["R"],start["chi"])
 chi1=np.interp(Rcommon,end["R"],end["chi"])
 chi_spatial0=np.interp(Rcommon,start["R"],start["chi_spatial"])
 chi_spatial1=np.interp(Rcommon,end["R"],end["chi_spatial"])
 chi_time0=np.interp(Rcommon,start["R"],start["chi_time"])
 chi_time1=np.interp(Rcommon,end["R"],end["chi_time"])
 # At fixed R, m=R(1-chi)/2, so compute the rate from interpolated chi.
 # Keep the separately interpolated-m rate to quantify its extra interpolation error.
 mdot_interpolated_mass=(m1-m0)/dt
 mdot=-0.5*Rcommon*(chi1-chi0)/dt
 rhs_mid=0.5*(q0+q1)
 mdot_spatial=-0.5*Rcommon*(chi_spatial1-chi_spatial0)/dt
 mdot_time=-0.5*Rcommon*(chi_time1-chi_time0)/dt
 mass_interpolation_gap=mdot_interpolated_mass-mdot
 rhs_comp_mid={}
 for name in start["rhs_components"]:
  s0=np.interp(Rcommon,start["R"],start["rhs_components"][name])
  s1=np.interp(Rcommon,end["R"],end["rhs_components"][name])
  rhs_comp_mid[name]=0.5*(s0+s1)
 residual=finite("fixed-R mass-current residual",mdot-rhs_mid)
 residual_interpolated_mass=finite("interpolated-mass fixed-R residual",mdot_interpolated_mass-rhs_mid)
 residual_decomp=finite("decomposed mass-current residual",
  mdot_spatial+mdot_time-sum(rhs_comp_mid.values()))
 decomp_gap=finite("decomposition closure gap",residual-residual_decomp)
 inner=Rcommon<=3.0
 annulus=(Rcommon>=2.0)&(Rcommon<=10.0)
 def endpoint_at_radius(radius,profile):
  chi_r=float(np.interp(radius,profile["R"],profile["chi"]))
  mass_r=0.5*radius*(1.0-chi_r)
  rhs_r=float(np.interp(radius,profile["R"],profile["rhs"]))
  comp_r={name:float(np.interp(radius,profile["R"],values)) for name,values in profile["rhs_components"].items()}
  return {"chi":chi_r,"mass":mass_r,"rhs":rhs_r,"components":comp_r}
 fixed={}
 radii=(0.5,1.0,2.0,3.0,5.0,8.0,10.0)
 endpoint_fixed={}
 for radius in radii:
  if lo<=radius<=hi:
   sR=endpoint_at_radius(radius,start);eR=endpoint_at_radius(radius,end)
   mdot_r=(eR["mass"]-sR["mass"])/dt
   mdot_chi_r=-0.5*radius*(eR["chi"]-sR["chi"])/dt
   rhs_mid_r=0.5*(sR["rhs"]+eR["rhs"])
   row={
    "mass_rate":mdot_chi_r,
    "mass_rate_interpolated_mass":(float(np.interp(radius,end["R"],end["mass"]))-float(np.interp(radius,start["R"],start["mass"])))/dt,
    "mass_interpolation_gap":((float(np.interp(radius,end["R"],end["mass"]))-float(np.interp(radius,start["R"],start["mass"])))/dt)-mdot_chi_r,
    "mass_rate_spatial_metric":float(np.interp(radius,Rcommon,mdot_spatial)),
    "mass_rate_time_shift_norm":float(np.interp(radius,Rcommon,mdot_time)),
    "current_rhs_midpoint":rhs_mid_r,
    "residual":mdot_chi_r-rhs_mid_r}
   row.update({f"current_{name}":0.5*(sR["components"][name]+eR["components"][name]) for name in sR["components"]})
   fixed[str(radius)]=row
   endpoint_fixed[radius]=(sR,eR)
 shell_balances=[]
 shell_bounds=tuple(r for r in (1.0,2.0,3.0,5.0,8.0,10.0) if r in endpoint_fixed)
 for rin,rout in zip(shell_bounds[:-1],shell_bounds[1:]):
  sin,ein=endpoint_fixed[rin];sout,eout=endpoint_fixed[rout]
  shell_mass_rate=((eout["mass"]-ein["mass"])-(sout["mass"]-sin["mass"]))/dt
  net_current=0.5*((sout["rhs"]-sin["rhs"])+(eout["rhs"]-ein["rhs"]))
  shell_residual=shell_mass_rate-net_current
  shell_balances.append({
   "R_in":rin,"R_out":rout,"shell_mass_rate":shell_mass_rate,
   "net_current_midpoint":net_current,"residual":shell_residual,
   "relative_abs_residual":abs(shell_residual)/max(abs(net_current),1e-300),
   "telescoping_residual_closure":shell_residual-((fixed[str(rout)]["residual"])-(fixed[str(rin)]["residual"]))
  })
 return {
  "dt":float(dt),"accepted_time":float(base.t),"probe_end_time":float(advanced.t),
  "steps_advanced_from_saved_input":1,"saved_input_immutable":True,
  "common_areal_radius_interval":[lo,hi],"common_samples":len(Rcommon),
  "max_abs_residual_R_le_3":float(np.max(np.abs(residual[inner]))),
  "rms_residual_R_2_to_10":float(np.sqrt(np.mean(residual[annulus]**2))),
  "max_abs_residual_R_2_to_10":float(np.max(np.abs(residual[annulus]))),
  "max_abs_decomposition_closure":float(np.max(np.abs(decomp_gap))),
  "max_abs_mass_interpolation_gap_R_le_3":float(np.max(np.abs(mass_interpolation_gap[inner]))),
  "rms_mass_interpolation_gap_R_2_to_10":float(np.sqrt(np.mean(mass_interpolation_gap[annulus]**2))),
  "max_abs_residual_using_interpolated_mass_R_le_3":float(np.max(np.abs(residual_interpolated_mass[inner]))),
  "rms_residual_using_interpolated_mass_R_2_to_10":float(np.sqrt(np.mean(residual_interpolated_mass[annulus]**2))),
  "rms_mass_rate_spatial_metric_R_2_to_10":float(np.sqrt(np.mean(mdot_spatial[annulus]**2))),
  "rms_mass_rate_time_shift_norm_R_2_to_10":float(np.sqrt(np.mean(mdot_time[annulus]**2))),
  "rms_current_components_R_2_to_10":{name:float(np.sqrt(np.mean(values[annulus]**2))) for name,values in rhs_comp_mid.items()},
  "fixed_areal_radius":fixed,
  "shell_balances":shell_balances,
  "max_abs_shell_balance_residual":float(max((abs(x["residual"]) for x in shell_balances),default=0.0)),
  "max_abs_shell_telescoping_closure":float(max((abs(x["telescoping_residual_closure"]) for x in shell_balances),default=0.0)),
  "start":{"R_cell0":float(start["R"][0]),"Rr_cell0":float(start["Rr"][0]),
           "Rt_cell0":float(start["Rt"][0]),"chi_cell0":float(start["chi"][0]),
           "mass_cell0":float(start["mass"][0]),"T_t_R_cell0":float(start["T_t_R"][0])},
  "end":{"R_cell0":float(end["R"][0]),"Rr_cell0":float(end["Rr"][0]),
         "Rt_cell0":float(end["Rt"][0]),"chi_cell0":float(end["chi"][0]),
         "mass_cell0":float(end["mass"][0]),"T_t_R_cell0":float(end["T_t_R"][0])},
  "max_abs_beta_start_end":max(start["max_abs_beta"],end["max_abs_beta"]),
  "input_state_immutable":before==after,
 }

def main()->int:
 if N not in (60,120,240):raise ValueError("resolution must be 60, 120 or 240")
 pin=origin_shadow.verify_source_pins()
 if pin.get("source_commit")!=SOURCE or pin.get("vendor_commit")!=VENDOR:
  raise RuntimeError("source/vendor pin mismatch")
 if pin.get("protected_engine_vendor_dependency_diff_empty") is not True:
  raise RuntimeError("protected source diff not empty")
 manifest=saved.check_artifact(FULL_DIR,N,full_profile=True)
 arrays=saved.load_profile_arrays(FULL_DIR)
 prefix=f"N{N}_full_continuous_metric_profile_final_"
 kernel,state=saved.reconstruct_state(N,arrays,prefix)
 origin_shadow._install_regular_row(kernel)
 if abs(float(state.t)-TARGET_TIME)>1e-12:raise RuntimeError("wrong accepted input time")
 dt=CFL*float(state.grid.dr)
 reports=[compare_step(kernel,state,dt),compare_step(kernel,state,dt/2.0)]
 report={
  "schema":"bssn-areal-misner-sharp-current-probe/v1","status":"completed_diagnostic_only",
  "source_pins":{"source_commit":SOURCE,"vendor_commit":VENDOR,"live_pin_check":pin,"full_profile_artifact":manifest},
  "runtime":{"python":sys.version,"numpy":np.__version__,"platform":platform.platform()},
  "configuration":{"resolution":N,"accepted_input_time":TARGET_TIME,"cfl_dt":dt,"dt_arms":[dt,dt/2.0],
   "areal_radius":"R=r*sqrt(b)/X","sqrt_gamma":"sqrt(a)*b*r^2/X^3",
   "chi":"(X^2/a)*R_r^2 - (R_t-beta*R_r)^2/alpha^2",
   "mass":"m_MS=R*(1-chi)/2 (G=1)",
   "current_transform":"T_t^R=(dx^mu/dy^t)(dy^R/dx^nu)T_mu^nu for y=(t,R(t,r)); zero shift enforced",
   "balance":"d m_MS/dt at fixed R compared with time-averaged 4*pi*R^2*T_t^R"},
  "probes":reports,
  "acceptance":{"source_pins_match":True,"accepted_artifact_sha_verified":True,
   "input_state_immutable":all(x["input_state_immutable"] for x in reports),
   "zero_shift_precondition":all(x["max_abs_beta_start_end"]<=1e-14 for x in reports),
   "all_profiles_finite":True,"integrated_current_identity_closed":False,
   "full_matter_coupled_propagation_identity_closed":False,
   "total_stress_energy_conservation_closed":False,
   "production_source_or_vendor_changed":False,"scientific_admission":False},
  "limitations":[
   "Finite-difference time rates use one copied step at dt and dt/2; temporal truncation remains and must be judged by the paired timestep comparison.",
   "The source identity is compared after transforming the mixed stress tensor and fixed-R derivative; this diagnostic does not authorize production changes.",
   "The archive's areal-coordinate identity is a source target; no source normalization or physical admission is implied."
  ],
  "scientific_admission":False}
 OUT.mkdir(parents=True,exist_ok=True)
 rp=OUT/"report.json";rp.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
 (OUT/"SHA256SUMS.txt").write_text(f"{sha(rp)}  report.json\n")
 print(json.dumps({"status":report["status"],"resolution":N,"dt_results":[
  {"dt":x["dt"],"max_abs_R_le_3":x["max_abs_residual_R_le_3"],"rms_R_2_to_10":x["rms_residual_R_2_to_10"],"fixed_R":x["fixed_areal_radius"]} for x in reports],
  "input_immutable":report["acceptance"]["input_state_immutable"],"scientific_admission":False},indent=2))
 return 0
if __name__=="__main__":raise SystemExit(main())
