"""Diagnostic decomposition of the existing D-field equation.

No evolution change. Records the already-used scalar RHS terms:
  D_t = alpha * PD
  PD_t = alpha*lap(D) + alpha*K*PD + X^-2*alpha_r*D_r + alpha*D
"""
from __future__ import annotations
import json, os
from pathlib import Path
import numpy as np
from engine.production_kernel import ProductionState, V55ProductionKernel
from engine.v55_initial import build_initial_data
from engine import v55_pirk_adapter as adapter

N=int(os.environ["RESOLUTION"]); RMAX=float(os.environ["RMAX"]); CFL=float(os.environ["CFL"])
TARGET=float(os.environ.get("TARGET","22.5")); DAMP=float(os.environ["DAMP"])
RAD=os.environ["RADIATION"]=="1"; CASE=os.environ["CASE_LABEL"]

class Probe(V55ProductionKernel):
    def initialize(self,resolution=N,r_max=RMAX,**kw):
        grid_ops,vac,_=adapter.vendor_modules()
        grid=grid_ops.SphericalCellGrid(resolution,r_max)
        init=build_initial_data(grid,vac.flat_state,D_amplitude=DAMP,include_radiation=RAD)
        init.geometry.alpha=self._solve_lapse(grid,init.geometry,init.scalars,init.matter)[0]
        init.geometry.beta.fill(0.); init.geometry.B.fill(0.)
        return ProductionState(grid=grid,geometry=init.geometry,scalars=init.scalars,matter=init.matter)

def rhs_terms(st):
    g=st.grid; q=st.geometry; D=st.scalars.D; PD=st.scalars.PD
    r=np.asarray(g.centers,float)
    d=lambda y,p=1: np.asarray(g.cell_derivative_fourth(np.asarray(y),parity=p),float)
    invr=q.X**2/q.a
    ap=d(q.a); bp=d(q.b); Xp=d(q.X); ar=d(q.alpha); Dr=d(D)
    Drr=np.asarray(g.cell_second_derivative_fourth(D,parity=1),float)
    lap=invr*(Drr-Dr*(.5*ap/q.a-bp/q.b+Xp/q.X-2./r))
    return r,D,PD,Dr,q.alpha*lap,q.alpha*q.K*PD,q.X**-2*ar*Dr,q.alpha*D

kernel=Probe(); state=kernel.initialize(); rows=[]
while state.t < TARGET:
    r,D,PD,Dr,lap,k,lapse,pot=rhs_terms(state)
    mask=(r>=1)&(r<=60); i=int(np.argmax(np.where(mask,np.abs(Dr),-np.inf)))
    den=max(abs(float(lap[i]))+abs(float(k[i]))+abs(float(lapse[i]))+abs(float(pot[i])),1e-300)
    w=mask
    def wl2(x): return float(np.sqrt(np.mean(np.asarray(x)[w]**2)))
    rows.append({
      "t":float(state.t),"tau":float(state.tau),"r_front":float(r[i]),
      "D_front":float(D[i]),"PD_front":float(PD[i]),"alpha_front":float(state.geometry.alpha[i]),
      "lap_D_front":float(lap[i]),"KPD_front":float(k[i]),
      "lapse_grad_front":float(lapse[i]),"potential_D_front":float(pot[i]),
      "signed_sum_front":float(lap[i]+k[i]+lapse[i]+pot[i]),
      "abs_fraction_lap":float(abs(lap[i])/den),
      "abs_fraction_K":float(abs(k[i])/den),
      "abs_fraction_lapse_grad":float(abs(lapse[i])/den),
      "abs_fraction_potential":float(abs(pot[i])/den),
      "L2_lap":wl2(lap),"L2_KPD":wl2(k),
      "L2_lapse_grad":wl2(lapse),"L2_potential":wl2(pot),
      "D_t_front":float(state.geometry.alpha[i]*PD[i]),
    })
    state=kernel.step(state,min(CFL*state.grid.dr,TARGET-state.t))

t=np.array([x["t"] for x in rows]); rf=np.array([x["r_front"] for x in rows])
m=t>=max(7.5,TARGET/2)
late=rows[int(np.argmax(t>=max(7.5,TARGET/2))):]
summary={
 "front_speed_late":float(np.polyfit(t[m],rf[m],1)[0]),
 "final_front":float(rf[-1]),"final_tau":float(state.tau),
 "final_term_values":{k:rows[-1][k] for k in ("lap_D_front","KPD_front","lapse_grad_front","potential_D_front","signed_sum_front")},
 "late_mean_abs_fractions":{k:float(np.mean([row[k] for row in late])) for k in ("abs_fraction_lap","abs_fraction_K","abs_fraction_lapse_grad","abs_fraction_potential")}
}
out=Path("runs/0star-d-field-budget")/CASE; out.mkdir(parents=True,exist_ok=True)
(out/"report.json").write_text(json.dumps({"case":CASE,"D_amplitude":DAMP,"radiation":RAD,"resolution":N,"rows":rows,"summary":summary},indent=2,default=float))
print(json.dumps(summary,indent=2,default=float))
