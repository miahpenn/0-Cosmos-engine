"""Vacuum check: PunctureGaugeKernel vs the vendor damped moving-puncture step on Schwarzschild."""
import sys; sys.path.insert(0, ".")
import numpy as np, dataclasses
from engine import v55_pirk_adapter as adapter
from engine.matter_system import _metric_arrays, initialize_dust, initialize_radiation
from engine.v55_matter import V55MatterState, metric_slice_from_q
from research.puncture_gauge_kernel_v1 import PunctureGaugeKernel
_, vac, moving = adapter.vendor_modules()
N, RMAX, M, ETA, KO = int(sys.argv[1]), 30.0, 1.0, 2.0, float(sys.argv[2])
k = PunctureGaugeKernel(eta=ETA, ko_epsilon=KO)
st = k.initialize(resolution=N, r_max=RMAX, include_radiation=False)
g = vac.flat_state(st.grid); r = np.asarray(st.grid.centers)
g.X = (1 + M / (2 * r)) ** -2; g.alpha.fill(1.0)
z = np.zeros_like(r); tiny = np.full_like(r, 1e-30)
mets = _metric_arrays(metric_slice_from_q(st.grid, g))
st = dataclasses.replace(st, geometry=g, scalars=dataclasses.replace(st.scalars, S=z, PS=z, D=z, PD=z, phi=z, Pi=z),
    matter=V55MatterState(initialize_dust(mets, tiny), initialize_dust(mets, tiny), initialize_radiation(mets, tiny), None))
import engine.scalar_system as ss; ss.V0 = 0.0
ref = g.copy(); dt = 0.1 * float(st.grid.dr); steps = int(sys.argv[3])
for n in range(steps):
    st = k.step(st, dt)
    ref = moving.damped_moving_puncture_pirk2_step(st.grid, ref, dt, eta=ETA, ko_epsilon=KO)
diffs = {f: float(np.max(np.abs(getattr(st.geometry, f) - getattr(ref, f)))) for f in ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B")}
print(f"N={N} ko={KO} t={st.t:.3f}: max |kernel - vendor| =", {f: f"{v:.1e}" for f, v in diffs.items()})
print("min alpha", float(st.geometry.alpha.min()), "max |beta|", float(np.max(abs(st.geometry.beta))))
