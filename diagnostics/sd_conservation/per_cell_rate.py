"""Diagnostic: per-cell Hamiltonian rate after the first (projection) step, from a
post-projection state. Measurement only."""
import math, sys, pathlib
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
dr = st0.grid.dr
st1 = k.step(st0, 0.03*dr)            # first step applies the projection
r = np.asarray(st1.grid.centers); dt = 0.0075*dr
Ha = Hprof(st1); st2 = k.step(st1, dt); Hb = Hprof(st2)
rate = (Hb - Ha)/dt
for i in list(range(0, 6)) + [8, 12, 16, 20, 30, 39]:
    print(f"cell {i:2d} r={r[i]:6.2f} rate={rate[i]:+.3e}")
print(f"max|rate| at cell {int(np.argmax(np.abs(rate)))}; sum|rate| cells0-4={np.sum(np.abs(rate[:5])):.3e}; all={np.sum(np.abs(rate)):.3e}")
