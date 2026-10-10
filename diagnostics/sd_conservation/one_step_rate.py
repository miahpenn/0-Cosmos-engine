"""One-step Hamiltonian rate test, measured from the POST-PROJECTION state.

The production step applies regularity projection internally, so the first step
from the raw initial state contains a dt-independent projection jump that inflates
dH/dt at small dt. This script takes one projecting step first, then measures the
rate from that state, so the jump is excluded. Controls: background (A=0) and D off.
Measurement only; no production change."""
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


def post_projection_state(k, amp, D, init):
    if init == "candidate":
        st, _, _ = discrete_consistent_state(
            k, 40, 40.0, amplitude=amp, D_amplitude=D
        )
    else:
        st = k.initialize(
            resolution=40, r_max=40.0, amplitude=amp, width=7.0,
            D_amplitude=D, include_radiation=True
        )
    dr = st.grid.dr
    st = k.step(st, 0.03*dr)  # projecting step, excluded from the rate
    return st


def rate_test(label, init, amp, D):
    k = V55ProductionKernel()
    base = post_projection_state(k, amp, D, init)
    H0 = Hprof(base)
    dr = base.grid.dr
    print(
        f"{label}: max|H| after projection (cells>=2) = "
        f"{np.max(np.abs(H0[2:])):.2e}"
    )
    for f in (0.03, 0.015, 0.0075):
        dt = f*dr
        st1 = k.step(base, dt)
        d = np.max(np.abs(Hprof(st1) - H0))
        print(f"   dt={dt:.5f}  max|dH|={d:.3e}  rate={d/dt:.3e}")


if __name__ == "__main__":
    rate_test("candidate A=0 (background)", "candidate", 0.0, 0.0)
    rate_test("candidate A=0.01 D=1e-10 (full S/D)", "candidate", 0.01, 1e-10)
    rate_test("candidate A=0.01 D=0 (S only)", "candidate", 0.01, 0.0)
    rate_test("baseline   A=0.01 D=1e-10", "baseline", 0.01, 1e-10)
