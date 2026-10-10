"""One-step constraint propagation: H(step(dt)) - H(0) vs dt.
Conservative flow at a constraint-satisfying state: change is O(dt^2) or smaller (rate -> 0).
First-order (rate ~ const) change means the RHS does not conserve H at that place."""
import sys, math, numpy as np
import pathlib as _p; sys.path.insert(0, str(_p.Path(__file__).resolve().parents[2]))
import engine.production_kernel as pk
from engine.discrete_consistent_initial import discrete_consistent_state
from engine.stress_energy import assemble_total_stress_energy
from engine import v55_pirk_adapter as ad
grid_ops, vacuum, _ = ad.vendor_modules()
def Hprof(st):
    raw=vacuum.constraints(st.grid,st.geometry); tot=assemble_total_stress_energy(st.grid,st.geometry,st.scalars,st.matter)
    return np.asarray(raw["hamiltonian"])-16*math.pi*np.asarray(tot.rho)
def mprof(st):
    raw=vacuum.constraints(st.grid,st.geometry); return np.asarray(raw["momentum"])
def test(label, init, amp, D):
    k=pk.V55ProductionKernel()
    if init=="candidate":
        st0,_,_=discrete_consistent_state(k,40,40.0,amplitude=amp,D_amplitude=D)
    else:
        st0=k.initialize(resolution=40,r_max=40.0,amplitude=amp,width=7.0,D_amplitude=D,include_radiation=True)
    H0=Hprof(st0); M0=mprof(st0); dr=st0.grid.dr
    out=[]
    for f in (0.03, 0.015, 0.0075):
        dt=f*dr
        st1=k.step(st0,dt)
        dH=Hprof(st1)-H0
        out.append((dt, np.max(np.abs(dH)), np.max(np.abs(dH[:5]))/dt, np.max(np.abs(dH))/dt))
    print(f"{label}: |M0| max={np.max(np.abs(M0[2:])):.2e}  |H0| max(cells>=2)={np.max(np.abs(H0[2:])):.2e}")
    for dt,dmax,rate0,rate in out:
        print(f"   dt={dt:.5f}  max|dH|={dmax:.3e}  max|dH|/dt={rate:.3e}  cells0-4 max|dH|/dt={rate0:.3e}")
test("candidate A=0 (background)", "candidate", 0.0, 0.0)
test("candidate A=0.01 D=1e-10 (full S/D)", "candidate", 0.01, 1e-10)
test("candidate A=0.01 D=0 (S only)", "candidate", 0.01, 0.0)
test("baseline   A=0.01 D=1e-10", "baseline", 0.01, 1e-10)
