"""Projection-displacement ledger with matched-time and scalar-amplitude controls.

Diagnostic-only: wraps the pinned regularity projection and measures its actual
per-call state displacement. No production source/equations are changed.
CFL comparison: resolution=40, r_max=40, t=3, dt/dr=(0.03, 0.015, 0.0075),
candidate S/D and zero-amplitude background. Matched-time snapshots are emitted
near t=(0.75, 1.5, 2.25, 3.0). A separate S-only amplitude sweep holds D=0 and
dt/dr=0.015 fixed at A=(0.005, 0.01, 0.02), isolating scalar-amplitude scaling.
Correlation between projection displacement and H drift is not proof of cause.
"""
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
DT_FACTORS = (0.03, 0.015, 0.0075)
TARGET_T = 3.0
SNAPSHOT_TIMES = (0.75, 1.5, 2.25, 3.0)


def Hprof(st):
    raw = vacuum.constraints(st.grid, st.geometry)
    tot = assemble_total_stress_energy(st.grid, st.geometry, st.scalars, st.matter)
    return np.asarray(raw["hamiltonian"]) - 16.0 * math.pi * np.asarray(tot.rho)


def new_ledger():
    return {
        name: {"signed": 0.0, "absolute": 0.0, "max_abs_single": 0.0}
        for name in FIELDS
    }


def run_case(label, amplitude, D_amplitude, dt_factor, matched_times=False):
    kernel = V55ProductionKernel()
    original_projection = vacuum.enforce_algebraic_regularity
    ledger = new_ledger()
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
        # Exclude initializer projection activity; measure evolution only.
        ledger = new_ledger()
        calls = 0
        H_initial = Hprof(state)
        dt = dt_factor * state.grid.dr
        nsteps = int(round(TARGET_T / dt))
        snapshot_steps = {
            int(round(t / dt)): t for t in SNAPSHOT_TIMES
        } if matched_times else {}
        snapshots = []
        if matched_times and 0 in snapshot_steps:
            snapshots.append((state.t, Hprof(state) - H_initial, {k: v.copy() for k, v in ledger.items()}))
        for step in range(1, nsteps + 1):
            state = kernel.step(state, dt)
            if step in snapshot_steps:
                snapshots.append((
                    state.t, Hprof(state) - H_initial,
                    {k: v.copy() for k, v in ledger.items()}
                ))
        H_final = Hprof(state)
        drift = H_final - H_initial
        print(
            f"{label}: A={amplitude:g} D={D_amplitude:g} dt/dr={dt_factor:.5g} "
            f"steps={nsteps} dt={dt:.6g} t_final={state.t:.6g} projection_calls={calls}"
        )
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
        for actual_t, h_drift, snap_ledger in snapshots:
            print(
                f"   MATCHED_TIME requested-near={actual_t:.6g} "
                f"actual_t={actual_t:.6g} H_cell0={h_drift[0]:+.6e} "
                f"max|H_drift|all={np.max(np.abs(h_drift)):.6e} "
                f"a_cell0_signed={snap_ledger['a']['signed']:+.6e} "
                f"a_cell0_abs={snap_ledger['a']['absolute']:.6e} "
                f"b_cell0_signed={snap_ledger['b']['signed']:+.6e} "
                f"b_cell0_abs={snap_ledger['b']['absolute']:.6e}"
            )
    finally:
        vacuum.enforce_algebraic_regularity = original_projection


def run_frequency_case(label, dt_factor, projection_mode):
    """Causal projection-frequency test; all evolution equations remain fixed."""
    kernel = V55ProductionKernel()
    original_projection = vacuum.enforce_algebraic_regularity
    ledger = new_ledger()
    total_calls = 0
    applied_calls = 0
    skipped_calls = 0
    evolution_started = False

    def frequency_projection(grid, geometry):
        nonlocal total_calls, applied_calls, skipped_calls
        if not evolution_started:
            return original_projection(grid, geometry)
        call_index = total_calls
        total_calls += 1
        # The current step implementation invokes the projection five times
        # per step (g0, predictor explicit block, primary block, final
        # explicit block, completed block). Keep/drop whole five-call groups.
        step_index = call_index // 5 + 1
        if projection_mode == "every-step":
            apply_now = True
        elif projection_mode == "every-second-step":
            apply_now = (step_index % 2 == 1)
        elif projection_mode == "initial-only":
            apply_now = False
        else:
            raise ValueError(f"unknown projection mode: {projection_mode}")
        if not apply_now:
            skipped_calls += 1
            return geometry
        projected = original_projection(grid, geometry)
        applied_calls += 1
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

    vacuum.enforce_algebraic_regularity = frequency_projection
    try:
        state = kernel.initialize(
            resolution=40, r_max=40.0, amplitude=0.01, width=7.0,
            D_amplitude=1.0e-10, include_radiation=True
        )
        ledger = new_ledger()
        total_calls = applied_calls = skipped_calls = 0
        H_initial = Hprof(state)
        dt = dt_factor * state.grid.dr
        nsteps = int(round(TARGET_T / dt))
        snapshots = []
        evolution_started = True
        for step in range(1, nsteps + 1):
            state = kernel.step(state, dt)
            if step in {int(round(t / dt)) for t in SNAPSHOT_TIMES}:
                h = Hprof(state) - H_initial
                snapshots.append((state.t, h.copy(), {
                    k: v.copy() for k, v in ledger.items()
                }))
        drift = Hprof(state) - H_initial
        print(
            f"FREQUENCY {label}: mode={projection_mode} A=0.01 D=1e-10 "
            f"dt/dr={dt_factor:.5g} steps={nsteps} t_final={state.t:.6g} "
            f"projection_calls={total_calls} applied={applied_calls} skipped={skipped_calls}"
        )
        print(
            f"   H drift cell0={drift[0]:+.6e}; "
            f"max|drift| cells0-4={np.max(np.abs(drift[:5])):.6e}; "
            f"max|drift| all={np.max(np.abs(drift)):.6e}"
        )
        for name in ("a", "b"):
            x = ledger[name]
            print(
                f"   {name}: cell0 signed_sum={x['signed']:+.6e} "
                f"abs_sum={x['absolute']:.6e} "
                f"max_single_any_cell={x['max_abs_single']:.6e}"
            )
        for actual_t, h_drift, snap_ledger in snapshots:
            print(
                f"   MATCHED_TIME t={actual_t:.6g} H_cell0={h_drift[0]:+.6e} "
                f"max|H_drift|all={np.max(np.abs(h_drift)):.6e} "
                f"a_signed={snap_ledger['a']['signed']:+.6e} "
                f"a_abs={snap_ledger['a']['absolute']:.6e} "
                f"b_signed={snap_ledger['b']['signed']:+.6e} "
                f"b_abs={snap_ledger['b']['absolute']:.6e}"
            )
    finally:
        vacuum.enforce_algebraic_regularity = original_projection


def constraint_budget(state):
    """Split H = H_geometry - 16*pi*(rho_scalar + rho_fluid)."""
    raw = vacuum.constraints(state.grid, state.geometry)
    total = assemble_total_stress_energy(
        state.grid, state.geometry, state.scalars, state.matter
    )
    return {
        "H_geometry": np.asarray(raw["hamiltonian"], dtype=float),
        "rho_scalar": np.asarray(total.scalar_rho, dtype=float),
        "rho_fluid": np.asarray(total.fluid_rho, dtype=float),
        "H_total": np.asarray(raw["hamiltonian"], dtype=float)
        - 16.0 * math.pi * np.asarray(total.rho, dtype=float),
    }


def run_constraint_budget():
    """Decompose the measured constraint drift without changing the evolution."""
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=40, r_max=40.0, amplitude=0.01, width=7.0,
        D_amplitude=1.0e-10, include_radiation=True
    )
    initial = constraint_budget(state)
    dt = 0.015 * state.grid.dr
    nsteps = int(round(TARGET_T / dt))
    snapshot_steps = {int(round(t / dt)): t for t in SNAPSHOT_TIMES}
    print("=== Hamiltonian drift budget: geometry vs scalar/fluid energy ===")
    print("Convention: Delta H = Delta H_geometry - 16*pi*Delta rho_scalar - 16*pi*Delta rho_fluid")
    for step in range(1, nsteps + 1):
        state = kernel.step(state, dt)
        if step not in snapshot_steps:
            continue
        now = constraint_budget(state)
        components = {
            "dH_geometry": now["H_geometry"] - initial["H_geometry"],
            "scalar_term": -16.0 * math.pi * (now["rho_scalar"] - initial["rho_scalar"]),
            "fluid_term": -16.0 * math.pi * (now["rho_fluid"] - initial["rho_fluid"]),
            "dH_total": now["H_total"] - initial["H_total"],
        }
        i = 0
        print(
            f"CONSTRAINT_BUDGET t={state.t:.6g} "
            f"cell0 dHgeom={components['dH_geometry'][i]:+.6e} "
            f"scalar={components['scalar_term'][i]:+.6e} "
            f"fluid={components['fluid_term'][i]:+.6e} "
            f"sum={components['dH_total'][i]:+.6e}"
        )
        for name, values in components.items():
            print(
                f"   {name}: maxabs_cells0-4={np.max(np.abs(values[:5])):.6e} "
                f"maxabs_all={np.max(np.abs(values)):.6e}"
            )


if __name__ == "__main__":
    print("=== CFL comparison: candidate S/D and background; matched-time ledger ===")
    for factor in DT_FACTORS:
        run_case("candidate S/D", 0.01, 1.0e-10, factor, matched_times=True)
        run_case("background", 0.0, 0.0, factor, matched_times=True)
    print("=== Scalar-amplitude scaling: S only (D=0), fixed dt/dr=0.015 ===")
    for amplitude in (0.005, 0.01, 0.02):
        run_case("S-only amplitude", amplitude, 0.0, 0.015, matched_times=True)
    run_constraint_budget()
    print("=== Projection-frequency causal test: fixed equations, S/D candidate ===")
    for factor in (0.03, 0.015):
        for mode in ("every-step", "every-second-step", "initial-only"):
            run_frequency_case("candidate S/D", factor, mode)
