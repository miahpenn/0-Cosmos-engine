"""Diagnostic: regularity-projection change of the candidate initial state, and
successive one-step Hamiltonian changes at fixed dt. Measurement only."""
import math, copy, sys, pathlib
import numpy as np
ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from engine.production_kernel import V55ProductionKernel
from engine.discrete_consistent_initial import discrete_consistent_state
from engine.stress_energy import assemble_total_stress_energy
from engine import v55_pirk_adapter as ad
grid_ops, vacuum, _ = ad.vendor_modules()
def Hprof(st):
    raw = vacuum.constraints(st.grid, st.geometry)
    tot = assemble_total_stress_energy(st.grid, st.geometry, st.scalars, st.matter)
    return np.asarray(raw["hamiltonian"]) - 16*math.pi*np.asarray(tot.rho)
k = V55ProductionKernel()
st0, _, _ = discrete_consistent_state(k, 40, 40.0)
dt = 0.03*st0.grid.dr
G = st0.geometry
Q = vacuum.enforce_algebraic_regularity(st0.grid, copy.deepcopy(G))
print("projection change of initial state: "
      f"max|da|={np.max(np.abs(np.asarray(Q.a)-np.asarray(G.a))):.3e} "
      f"max|db|={np.max(np.abs(np.asarray(Q.b)-np.asarray(G.b))):.3e} "
      f"max|dX|={np.max(np.abs(np.asarray(Q.X)-np.asarray(G.X))):.3e} "
      f"max|dK|={np.max(np.abs(np.asarray(Q.K)-np.asarray(G.K))):.3e}")
st = st0; Hp = Hprof(st)
for n in range(1, 5):
    st = k.step(st, dt); H = Hprof(st)
    print(f"step {n}: max|dH| cells0-4={np.max(np.abs((H-Hp)[:5])):.3e} all={np.max(np.abs(H-Hp)):.3e}")
    Hp = H
