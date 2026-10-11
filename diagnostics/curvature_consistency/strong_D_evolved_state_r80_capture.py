"""Full-state R=80 trajectory capture for a matched-domain history comparison.

This is diagnostic-only. It reproduces the original production parameters at
N=160, r_max=80, dr=0.5 and saves all evolved state arrays every 0.25 time
units. Compare these snapshots with the already archived N=320, r_max=160
trajectory restricted to its first 160 cells. No physics is changed.
"""
from __future__ import annotations
import json
import math
import pathlib
import sys
import traceback
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from engine import v55_pirk_adapter as adapter
from engine.production_kernel import V55ProductionKernel
from engine.stress_energy import assemble_total_stress_energy
from diagnostics.curvature_consistency.independent_metric_ricci import finite_array, metric_connection

grid_ops, vacuum, _ = adapter.vendor_modules()

RESOLUTION = 160
R_MAX = 80.0
AMPLITUDE = 0.01
WIDTH = 7.0
D_AMPLITUDE = 1.0e-4
INCLUDE_RADIATION = True
CFL = 0.0075
FINAL_TIME = 12.0
SAMPLE_INTERVAL = 0.25
TARGET_TIMES = tuple(float(x) for x in np.arange(
    SAMPLE_INTERVAL, FINAL_TIME + 0.5 * SAMPLE_INTERVAL, SAMPLE_INTERVAL
))
OUT_DIR = ROOT / "runs" / "strong-D-evolved-state-r80-capture"
GEOMETRY_ARRAYS = ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B")
SCALAR_ARRAYS = ("S", "PS", "D", "PD", "phi", "Pi")
MATTER_NAMES = ("dark_matter", "baryons", "radiation")
CONSERVED_ARRAYS = ("rest", "energy_t", "momentum_r")

def state_snapshot(state):
    geometry = {
        name: [float(x) for x in np.asarray(getattr(state.geometry, name), dtype=float)]
        for name in GEOMETRY_ARRAYS
    }
    scalars = {
        name: [float(x) for x in np.asarray(getattr(state.scalars, name), dtype=float)]
        for name in SCALAR_ARRAYS
    }
    matter = {}
    for species_name in MATTER_NAMES:
        species = getattr(state.matter, species_name)
        matter[species_name] = {
            name: [float(x) for x in np.asarray(getattr(species, name), dtype=float)]
            for name in CONSERVED_ARRAYS
        }
    return {"geometry": geometry, "scalars": scalars, "matter": matter}

def sample_state(state, requested_time, label):
    grid, geom = state.grid, state.geometry
    radius = np.asarray(grid.centers, dtype=float)
    raw = vacuum.constraints(grid, geom)
    total = assemble_total_stress_energy(grid, geom, state.scalars, state.matter)
    rho_scalar = finite_array("scalar rho", np.asarray(total.scalar_rho, dtype=float))
    rho_fluid = finite_array("fluid rho", np.asarray(total.fluid_rho, dtype=float))
    H = finite_array(
        "vendor H", np.asarray(raw["hamiltonian"], dtype=float) - 16.0 * math.pi * (rho_scalar + rho_fluid)
    )
    alpha = finite_array("alpha", np.asarray(geom.alpha, dtype=float))
    C_lambda = finite_array(
        "connection constraint", np.asarray(geom.Lambda, dtype=float) - metric_connection(grid, geom)
    )
    radius_idx = int(np.argmax(np.abs(H)))
    off = np.abs(H[2:])
    off_idx = int(np.argmax(off)) + 2 if len(off) else 0
    amin_idx = int(np.argmin(alpha))
    snapshot = state_snapshot(state)
    for group, arrays in snapshot.items():
        for name, values in arrays.items():
            if group == "matter" and isinstance(values, dict):
                for field, data in values.items():
                    finite_array(f"{group}.{name}.{field}", np.asarray(data, dtype=float))
            elif isinstance(values, list):
                finite_array(f"{group}.{name}", np.asarray(values, dtype=float))
    return {
        "schema": "strong_D_evolved_state_r80_capture_v1",
        "diagnostic_only": True,
        "domain": {"N": int(grid.n), "r_max": float(R_MAX), "dr": float(grid.dr)},
        "label": label,
        "requested_time": float(requested_time),
        "t": float(state.t),
        "tau": float(state.tau),
        "e_folds": float(state.e_folds),
        "cycle_event_count": len(state.cycle.events),
        "handoff_count": len(state.handoffs),
        "central_H_vendor": float(H[0]),
        "max_abs_H_all": float(np.max(np.abs(H))),
        "max_abs_H_all_radius": float(radius[radius_idx]),
        "max_abs_H_offcentre_cells_2plus": float(np.max(off)) if len(off) else 0.0,
        "max_abs_H_offcentre_radius": float(radius[off_idx]),
        "C_lambda_center": float(C_lambda[0]),
        "max_abs_C_lambda": float(np.max(np.abs(C_lambda))),
        "min_a": float(np.min(geom.a)),
        "min_b": float(np.min(geom.b)),
        "min_X": float(np.min(geom.X)),
        "min_alpha": float(alpha[amin_idx]),
        "min_alpha_radius": float(radius[amin_idx]),
        "radius": [float(x) for x in radius],
        "H_vendor": [float(x) for x in H],
        "alpha_profile": [float(x) for x in alpha],
        "state_snapshot": snapshot,
    }

def write_progress(state, accepted_steps, records, status, failure=None):
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema": "strong_D_evolved_state_r80_capture_v1",
        "status": status,
        "diagnostic_only": True,
        "configuration": {
            "N": RESOLUTION, "r_max": R_MAX, "dr": 0.5, "CFL": CFL,
            "final_time": FINAL_TIME, "sample_interval": SAMPLE_INTERVAL,
            "amplitude": AMPLITUDE, "width": WIDTH, "D_amplitude": D_AMPLITUDE,
            "include_radiation": INCLUDE_RADIATION,
            "production_physics_changed": False, "production_defaults_changed": False,
            "gauge_changed": False, "projection_changed": False,
            "evolution_boundary_conditions_changed": False,
        },
        "accepted_steps": accepted_steps,
        "sample_count": len(records),
        "sample_times": [float(row["t"]) for row in records],
        "failure": failure,
    }
    (OUT_DIR / "campaign_progress.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\\n"
    )

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    profiles_path = OUT_DIR / "rmax080_N160_production_full_state_profiles.jsonl"
    profiles_path.write_text("")
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=RESOLUTION, r_max=R_MAX, amplitude=AMPLITUDE, width=WIDTH,
        D_amplitude=D_AMPLITUDE, include_radiation=INCLUDE_RADIATION,
    )
    if abs(float(state.grid.dr) - 0.5) > 1.0e-13:
        raise ValueError(f"unexpected grid spacing: {state.grid.dr}")
    dt_nominal = CFL * float(state.grid.dr)
    records = []
    accepted_steps = 0
    failure = None

    initial = sample_state(state, 0.0, "t=0")
    records.append(initial)
    with profiles_path.open("a") as stream:
        stream.write(json.dumps(initial, sort_keys=True, allow_nan=False) + "\\n")
        stream.flush()
    write_progress(state, accepted_steps, records, "in_progress")

    target_index = 0
    while state.t < FINAL_TIME - 1.0e-13:
        dt = min(dt_nominal, FINAL_TIME - state.t)
        try:
            next_state = kernel.step(state, dt)
            next_state.geometry.assert_finite_positive()
            for name, array in (
                ("S", next_state.scalars.S), ("D", next_state.scalars.D),
                ("phi", next_state.scalars.phi), ("a", next_state.geometry.a),
                ("b", next_state.geometry.b), ("X", next_state.geometry.X),
                ("Lambda", next_state.geometry.Lambda), ("alpha", next_state.geometry.alpha),
            ):
                finite_array(name, np.asarray(array, dtype=float))
            state = next_state
            accepted_steps += 1
            if accepted_steps % 128 == 0:
                print("[STRONG_D_R80_STATE_CAPTURE] HEARTBEAT " + json.dumps({
                    "t": float(state.t), "accepted_steps": accepted_steps,
                }, sort_keys=True), flush=True)
        except Exception as exc:
            failure = {
                "stage": "evolution_step", "step_attempt": accepted_steps + 1,
                "t_before": float(state.t), "error_type": type(exc).__name__,
                "error": str(exc), "traceback": traceback.format_exc(limit=12),
            }
            break

        while target_index < len(TARGET_TIMES) and (
            state.t + 0.5 * dt >= TARGET_TIMES[target_index]
            or abs(state.t - TARGET_TIMES[target_index]) < 1.0e-12
        ):
            target = TARGET_TIMES[target_index]
            sample = sample_state(state, target, f"target={target:g}")
            records.append(sample)
            with profiles_path.open("a") as stream:
                stream.write(json.dumps(sample, sort_keys=True, allow_nan=False) + "\\n")
                stream.flush()
            target_index += 1
            write_progress(state, accepted_steps, records, "in_progress", failure)
            print("[STRONG_D_R80_STATE_CAPTURE] CHECKPOINT " + json.dumps({
                "t": float(state.t), "accepted_steps": accepted_steps,
                "sample_count": len(records), "central_H_vendor": sample["central_H_vendor"],
                "max_abs_H": sample["max_abs_H_all"],
            }, sort_keys=True), flush=True)

    completed = failure is None and state.t >= FINAL_TIME - 1.0e-10
    status = "completed" if completed and len(records) == len(TARGET_TIMES) + 1 else "numerical_or_sampling_failure"
    write_progress(state, accepted_steps, records, status, failure)
    report = {
        "schema": "strong_D_evolved_state_r80_capture_v1",
        "status": status, "diagnostic_only": True,
        "accepted_steps": accepted_steps, "sample_count": len(records),
        "final_t": float(state.t), "final_tau": float(state.tau),
        "final_e_folds": float(state.e_folds),
        "cycle_event_count": len(state.cycle.events), "handoff_count": len(state.handoffs),
        "failure": failure, "profiles_file": profiles_path.name,
        "configuration": {
            "N": RESOLUTION, "r_max": R_MAX, "dr": float(state.grid.dr),
            "CFL": CFL, "final_time": FINAL_TIME, "sample_interval": SAMPLE_INTERVAL,
            "amplitude": AMPLITUDE, "width": WIDTH, "D_amplitude": D_AMPLITUDE,
            "include_radiation": INCLUDE_RADIATION,
            "production_physics_changed": False, "production_defaults_changed": False,
            "gauge_changed": False, "projection_changed": False,
            "evolution_boundary_conditions_changed": False,
        },
    }
    (OUT_DIR / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\\n")
    print("[STRONG_D_R80_STATE_CAPTURE] FINAL " + json.dumps(report, sort_keys=True), flush=True)
    if status != "completed":
        raise RuntimeError("R=80 full-state diagnostic capture did not complete cleanly.")

if __name__ == "__main__":
    main()
