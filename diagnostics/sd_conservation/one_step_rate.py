"""One-step Hamiltonian scaling from a projection-only base.

The regularity projection is applied as a state map; no physical evolution step
is used to construct the base. The projection's one-time H offset is reported
separately. Measurement only; production equations are untouched.
"""
import copy
import math
import pathlib
import sys
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from engine.production_kernel import V55ProductionKernel
from engine.discrete_consistent_initial import discrete_consistent_state
from engine.stress_energy import assemble_total_stress_energy
from engine import v55_pirk_adapter as ad

_, vacuum, _ = ad.vendor_modules()


def Hprof(st):
    raw = vacuum.constraints(st.grid, st.geometry)
    tot = assemble_total_stress_energy(st.grid, st.geometry, st.scalars, st.matter)
    return np.asarray(raw["hamiltonian"]) - 16.0 * math.pi * np.asarray(tot.rho)


def projection_only_base(k, amp, D, init):
    if init == "candidate":
        st, _, _ = discrete_consistent_state(
            k, 40, 40.0, amplitude=amp, D_amplitude=D
        )
    else:
        st = k.initialize(
            resolution=40, r_max=40.0, amplitude=amp, width=7.0,
            D_amplitude=D, include_radiation=True
        )
    Hraw = Hprof(st)
    projected = copy.copy(st)
    projected.geometry = vacuum.enforce_algebraic_regularity(
        st.grid, st.geometry.copy()
    )
    Hproj = Hprof(projected)
    print(
        f"   projection-only offset: cell0={Hproj[0]-Hraw[0]:+.3e}; "
        f"max|dH|={np.max(np.abs(Hproj-Hraw)):.3e}"
    )
    return projected, Hraw, Hproj


def rate_test(label, init, amp, D):
    k = V55ProductionKernel()
    base, _, H0 = projection_only_base(k, amp, D, init)
    dr = base.grid.dr
    print(
        f"{label}: projection-only base max|H| cells>=2="
        f"{np.max(np.abs(H0[2:])):.3e}; cell0={H0[0]:+.3e}"
    )
    dts, changes = [], []
    for f in (0.03, 0.015, 0.0075):
        dt = f * dr
        st1 = k.step(base, dt)
        delta = Hprof(st1) - H0
        d = float(np.max(np.abs(delta)))
        dts.append(dt)
        changes.append(d)
        print(f"   dt={dt:.5f}  max|dH|={d:.3e}  rate={d/dt:.3e}")
    if min(changes) > 0.0:
        slope = float(np.polyfit(np.log(dts), np.log(changes), 1)[0])
        print(f"   log-log slope of max|dH| vs dt = {slope:.4f} (2 is second-order scaling)")


if __name__ == "__main__":
    rate_test("candidate A=0 (background)", "candidate", 0.0, 0.0)
    rate_test("candidate A=0.01 D=1e-10 (full S/D)", "candidate", 0.01, 1e-10)
    rate_test("candidate A=0.01 D=0 (S only)", "candidate", 0.01, 0.0)
    rate_test("baseline A=0.01 D=1e-10 (non-admissible H base)", "baseline", 0.01, 1e-10)
