"""Initial-slice radial-quadrature control; no evolution or production edits."""
import csv, inspect, json, math, os
from pathlib import Path
import numpy as np
from scipy.integrate import cumulative_simpson
import engine.production_kernel as pk
import engine.v55_initial as initial_module
from engine.stress_energy import assemble_total_stress_energy
from engine import v55_pirk_adapter as adapter
from engine.discrete_consistent_initial import discrete_consistent_state

OUT=Path("runs/radial_quadrature_control"); OUT.mkdir(parents=True,exist_ok=True)
N_VALUES=(40,80,160,320); RMAX=40.0

def maxabs(x):
    x=np.asarray(x,dtype=float)
    return float(np.max(np.abs(x))) if x.size else 0.0

def high_order_builder():
    source=inspect.getsource(initial_module.build_initial_data)
    old='''integ = np.cumsum(
            np.r_[0.5 * dr * src[0],
                  0.5 * dr * (src[1:] + src[:-1])]
        )'''
    new='''# The radial source is src(r) = 1 + r^2*F(r), so src(0) is exactly 1.
        xq = np.r_[0.0, r]
        yq = np.r_[1.0, src]
        integ = cumulative_simpson(yq, x=xq, initial=0.0)[1:]'''
    if old not in source: raise RuntimeError("radial integral source changed; refusing to patch")
    ns=dict(initial_module.__dict__); ns["cumulative_simpson"]=cumulative_simpson
    exec(compile(source.replace(old,new),"<quadrature_control_initializer>","exec"),ns)
    return ns["build_initial_data"]

def measure(state,mode,n):
    _,vacuum,_=adapter.vendor_modules()
    raw=vacuum.constraints(state.grid,state.geometry)
    total=assemble_total_stress_energy(state.grid,state.geometry,state.scalars,state.matter)
    Hraw=np.asarray(raw["hamiltonian"],dtype=float)
    source=16.0*math.pi*np.asarray(total.rho,dtype=float)
    H=Hraw-source; r=np.asarray(state.grid.centers)
    off_mask=(r >= 2.0)
    if not np.any(off_mask): raise RuntimeError("no grid cells at r >= 2 for off-centre audit")
    off_indices=np.flatnonzero(off_mask)
    oi=int(off_indices[np.argmax(np.abs(H[off_mask]))])
    oi_cells=int(np.argmax(np.abs(H[2:]))+2)
    row={"mode":mode,"N":n,"R":RMAX,"dr":float(state.grid.dr),"t":float(state.t),
      "H_center":float(H[0]),"H_center_abs":float(abs(H[0])),
      "H_offmax_r_ge_2":maxabs(H[off_mask]),"H_offmax_r_ge_2_r":float(r[oi]),
      "H_offmax_excluding_first_two":maxabs(H[2:]),"H_offmax_excluding_first_two_r":float(r[oi_cells]),
      "raw_constraint_center":float(Hraw[0]),"16pi_rho_center":float(source[0]),
      "scalar_rho_center":float(total.scalar_rho[0]),"fluid_rho_center":float(total.fluid_rho[0]),
      "H_max_all":maxabs(H),"M_center":float(np.asarray(raw["momentum"])[0]-8*math.pi*np.asarray(total.j)[0]),
      "finite":bool(np.all(np.isfinite(H)))}
    if not row["finite"]: raise RuntimeError(f"nonfinite residual {mode} N={n}")
    print("RESULT "+json.dumps(row,sort_keys=True),flush=True)
    return row

def checkpoint(rows):
    OUT.mkdir(parents=True,exist_ok=True)
    payload={"schema":"initial_slice_radial_quadrature_control_partial_v1","R":RMAX,
      "resolutions":list(N_VALUES),"evolution_performed":False,"rows":rows}
    (OUT/"partial_rows.json").write_text(json.dumps(payload,indent=2,sort_keys=True)+"\\n")

def run():
    rows=[]; kernel=pk.V55ProductionKernel(); original=pk.build_initial_data
    for n in N_VALUES:
        pk.build_initial_data=original
        state=kernel.initialize(resolution=n,r_max=RMAX,amplitude=0.01,width=7.0,D_amplitude=1e-4,include_radiation=True)
        rows.append(measure(state,"production_trapezoid",n)); checkpoint(rows)
        pk.build_initial_data=high_order_builder()
        state=kernel.initialize(resolution=n,r_max=RMAX,amplitude=0.01,width=7.0,D_amplitude=1e-4,include_radiation=True)
        rows.append(measure(state,"high_order_cumulative_simpson",n)); checkpoint(rows)
        pk.build_initial_data=original
        state,B,hist,info=discrete_consistent_state(kernel,resolution=n,r_max=RMAX,amplitude=0.01,width=7.0,D_amplitude=1e-4,include_radiation=True,max_iter=16,tol=2e-13,max_cond=1e12,return_info=True)
        state.geometry.alpha=np.asarray(kernel._solve_lapse(state.grid,state.geometry,state.scalars,state.matter)[0],dtype=float).copy()
        state.geometry.beta.fill(0.0); state.geometry.B.fill(0.0)
        row=measure(state,"discrete_consistent",n)
        row.update({"solver_iterations":info["iterations"],"solver_final_max_residual":info["final_max_residual"],"solver_max_condition_number":info["max_condition_number"],"solver_residual_history":[float(x) for x in hist]})
        rows.append(row); checkpoint(rows)
    pk.build_initial_data=original
    by={}
    for x in rows: by.setdefault(x["N"],{})[x["mode"]]=x
    comparisons=[]
    for n,g in sorted(by.items()):
        comparisons.append({"N":n,"production_H_center":g["production_trapezoid"]["H_center"],
          "high_order_H_center":g["high_order_cumulative_simpson"]["H_center"],
          "high_order_minus_production_H_center":g["high_order_cumulative_simpson"]["H_center"]-g["production_trapezoid"]["H_center"],
          "discrete_consistent_H_center":g["discrete_consistent"]["H_center"],
          "production_offmax_r_ge_2":g["production_trapezoid"]["H_offmax_r_ge_2"],
          "high_order_offmax_r_ge_2":g["high_order_cumulative_simpson"]["H_offmax_r_ge_2"],
          "discrete_consistent_H_max_all":g["discrete_consistent"]["H_max_all"]})
    report={"schema":"initial_slice_radial_quadrature_control_v1","source_commit":os.environ.get("GITHUB_SHA","unavailable"),
      "R":RMAX,"resolutions":list(N_VALUES),"off_centre_radius_threshold":2.0,
      "discrete_solver_max_iter":16,"discrete_solver_tolerance":2e-13,
      "evolution_performed":False,
      "constraint_operator_unchanged":True,"hamiltonian_definition_unchanged":True,
      "quadrature":"cumulative Simpson on [0, cell centers], using even central source value at origin",
      "term_breakdown":"H_raw - 16*pi*(scalar_rho + fluid_rho); raw, scalar, and fluid terms recorded separately",
      "rows":rows,"comparisons":comparisons,
      "interpretation_guard":"No automatic attribution or physics admission."}
    (OUT/"report.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    with (OUT/"report.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    print("SUMMARY "+json.dumps(comparisons,sort_keys=True),flush=True)
    print("PASS: all initial-slice scans complete; no evolution",flush=True)

if __name__=="__main__": run()
