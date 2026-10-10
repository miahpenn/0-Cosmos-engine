"""Short-run ledger of regularity-projection displacement versus Hamiltonian drift.

Diagnostic only. The pinned projection routine is wrapped to measure its actual
per-call state displacement during evolution; no production source is changed.
Uses resolution=40, r_max=40, dt=0.03*dr, target t=3.
"""
import copy
import math
import pathlib
import sys
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from engine.production_kernel import V55ProductionKernel
from engine.stress_energy import assemble_total_stress_energy
from engine import v55_pirk_adapter as ad

_, vacuum, _ = ad.vendor_modules()
FIELDS = ("a", "b", "X", "Aa", "K", "Lambda")


def Hprof(st):
    raw = vacuum.constraints(st.grid, st.geometry)
    tot = assemble_total_stress_energy(st.grid, st.geometry, st.scalars, st.matter)
    return np.asarray(raw["hamiltonian"]) - 16.0 * math.pi * np.asarray(tot.rho)


def run_case(label, amplitude, D_amplitude):
    kernel = V55ProductionKernel()
    original_projection = vacuum.enforce_algebraic_regularity
    ledger = {
        name: {"signed": 0.0, "absolute": 0.0, "max_abs_single": 0.0}
        for name in FIELDS
    }
    calls = 0

    def measured_projection(grid, geometry):
        nonlocal calls
        projected = original_projection(grid, geometry)
        calls += 1
        for name in FIELDS:
            before = np.asarray(getattr(geometry, name), dtype=float)
            after = np.asarray(getattr(projected, name), dtype=float)
            delta = after - before
            ledger[name]["signed"] += float(delta[0])
            ledger[name]["absolute"] += float(abs(delta[0]))
            ledger[name]["max_abs_single"] = max(
                ledger[name]["max_abs_single"], float(np.max(np.abs(delta)))
            )
        return projected

    vacuum.enforce_algebraic_regularity = measured_projection
    try:
        state = kernel.initialize(
            resolution=40, r_max=40.0, amplitude=amplitude, width=7.0,
            D_amplitude=D_amplitude, include_radiation=True
        )
        # Exclude initializer activity; ledger is specifically for evolution.
        for values in ledger.values():
            values.update(signed=0.0, absolute=0.0, max_abs_single=0.0)
        calls = 0
        H_initial = Hprof(state)
        dt = 0.03 * state.grid.dr
        target_t = 3.0
        nsteps = int(round(target_t / dt))
        for _ in range(nsteps):
            state = kernel.step(state, dt)
        H_final = Hprof(state)
        drift = H_final - H_initial
        print(f"{label}: steps={nsteps} dt={dt:.6g} t_final={state.t:.6g} projection_calls={calls}")
        print(
            f"   H drift cell0={drift[0]:+.6e}; "
            f"max|drift| cells0-4={np.max(np.abs(drift[:5])):.6e}; "
            f"max|drift| all={np.max(np.abs(drift)):.6e}"
        )
        for name in FIELDS:
            x = ledger[name]
            print(
                f"   {name}: cell0 signed_sum={x['signed']:+.6e} "
                f"abs_sum={x['absolute']:.6e} "
                f"max_single_any_cell={x['max_abs_single']:.6e}"
            )
    finally:
        vacuum.enforce_algebraic_regularity = original_projection


if __name__ == "__main__":
    run_case("candidate S/D", 0.01, 1.0e-10)
    run_case("background", 0.0, 0.0)
