"""Diagnostic-only frozen-state crop audit of the native CMC lapse solve.

One production trajectory is evolved on the largest domain. At each saved slice,
the exact current state is cropped to the smaller radii and the existing CMC
lapse equation is solved again. This separates direct domain/truncation effects
in the elliptic solve from differences accumulated by separately evolved runs.

No production equation, gauge, boundary condition, projection, or evolution step
is changed. Alternative RHS/boundary values are diagnostic counterfactual solves
only, explicitly labelled in the output.
"""
from __future__ import annotations

import copy
import json
import math
import pathlib
import sys
import traceback

import numpy as np
from scipy.linalg import solve_banded

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from engine import v55_pirk_adapter as adapter
from engine.cmc_gauge import _cmc_operator_banded, solve_cmc_lapse
from engine.matter_system import ConservedSpecies
from engine.production_kernel import V55ProductionKernel
from engine.scalar_system import ScalarFields
from engine.v55_matter import V55MatterState, total_matter_projection

grid_ops, vacuum, _ = adapter.vendor_modules()

AMPLITUDE = 0.01
WIDTH = 7.0
D_AMPLITUDE = 1.0e-4
INCLUDE_RADIATION = True
CFL = 0.0075
FINAL_TIME = 12.0
SAMPLE_INTERVAL = 0.25
FULL_N = 320
FULL_R_MAX = 160.0
CROPS = (
    {"N": 160, "r_max": 80.0},
    {"N": 240, "r_max": 120.0},
    {"N": 320, "r_max": 160.0},
)
TARGET_TIMES = tuple(float(x) for x in np.arange(
    SAMPLE_INTERVAL, FINAL_TIME + 0.5 * SAMPLE_INTERVAL, SAMPLE_INTERVAL
))
OUT_DIR = ROOT / "runs" / "cmc-lapse-frozen-state-crop-audit"
GEOMETRY_ARRAYS = ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B")
SCALAR_ARRAYS = ("S", "PS", "D", "PD", "phi", "Pi")
MATTER_NAMES = ("dark_matter", "baryons", "radiation")
CONSERVED_ARRAYS = ("rest", "energy_t", "momentum_r")


def _crop_geometry(geometry, full_n: int, n: int):
    out = geometry.copy()
    values = vars(out)
    cropped = 0
    for name, value in list(values.items()):
        if isinstance(value, np.ndarray) and value.ndim == 1 and value.shape[0] == full_n:
            setattr(out, name, value[:n].copy())
            cropped += 1
    missing = [name for name in GEOMETRY_ARRAYS
               if not isinstance(getattr(out, name, None), np.ndarray)
               or getattr(out, name).shape != (n,)]
    if missing:
        raise ValueError(f"geometry crop missing expected length-{n} arrays: {missing}")
    if cropped < len(GEOMETRY_ARRAYS):
        raise ValueError(f"only cropped {cropped} geometry arrays; expected at least {len(GEOMETRY_ARRAYS)}")
    return out


def _crop_scalars(fields: ScalarFields, n: int) -> ScalarFields:
    return ScalarFields(*(np.asarray(getattr(fields, name))[:n].copy()
                          for name in SCALAR_ARRAYS))


def _crop_matter(matter: V55MatterState, n: int) -> V55MatterState:
    out = matter.copy()
    for name in MATTER_NAMES:
        src = getattr(out, name)
        setattr(out, name, ConservedSpecies(
            np.asarray(src.rest)[:n].copy(),
            np.asarray(src.energy_t)[:n].copy(),
            np.asarray(src.momentum_r)[:n].copy(),
        ))
    return out


def _cmc_coefficients(grid, geometry, scalars, matter):
    total = total_matter_projection(grid, geometry, scalars, matter)
    a = np.asarray(geometry.a, dtype=float)
    b = np.asarray(geometry.b, dtype=float)
    X = np.asarray(geometry.X, dtype=float)
    r = np.asarray(grid.centers, dtype=float)
    inv = X * X / a
    ap = np.asarray(grid.cell_derivative_fourth(a, parity=1), dtype=float)
    bp = np.asarray(grid.cell_derivative_fourth(b, parity=1), dtype=float)
    Xp = np.asarray(grid.cell_derivative_fourth(X, parity=1), dtype=float)
    c = -0.5 * ap / a + bp / b - Xp / X + 2.0 / r
    Q = (
        1.5 * np.asarray(geometry.Aa, dtype=float) ** 2
        + np.asarray(geometry.K, dtype=float) ** 2 / 3.0
        + 4.0 * math.pi * (
            np.asarray(total["rho"], dtype=float)
            + np.asarray(total["pr"], dtype=float)
            + 2.0 * np.asarray(total["pt"], dtype=float)
        )
    )
    m = Q / inv
    return c, Q, inv, m


def _banded_matvec(band: np.ndarray, values: np.ndarray, lower=4, upper=4):
    n = len(values)
    out = np.zeros(n, dtype=float)
    for j in range(n):
        lo = max(0, j - upper)
        hi = min(n - 1, j + lower)
        for i in range(lo, hi + 1):
            out[i] += band[upper + i - j, j] * values[j]
    return out


def _solve_with_diagnostic_rhs(grid, geometry, scalars, matter,
                               kdot: float, outer_alpha: float):
    """Use the native discrete operator with explicitly diagnostic RHS/BC overrides."""
    c, Q, inv, m = _cmc_coefficients(grid, geometry, scalars, matter)
    band = _cmc_operator_banded(grid, c, m)
    rhs = np.zeros(grid.n, dtype=float)
    rhs[:-1] = -float(kdot) / inv[:-1]
    rhs[-1] = float(outer_alpha)
    alpha = solve_banded((4, 4), band, rhs, check_finite=False)
    if not np.all(np.isfinite(alpha)):
        raise FloatingPointError("counterfactual CMC solve returned non-finite lapse")
    residual = _banded_matvec(band, alpha) - rhs
    return alpha, {
        "max_abs_discrete_operator_residual": float(np.max(np.abs(residual))),
        "max_abs_interior_operator_residual": float(np.max(np.abs(residual[:-1]))),
        "boundary_row_residual": float(residual[-1]),
        "c": [float(x) for x in c],
        "Q": [float(x) for x in Q],
        "inv": [float(x) for x in inv],
        "m": [float(x) for x in m],
    }


def _profile_stats(radius, lhs, rhs):
    delta = np.asarray(lhs, dtype=float) - np.asarray(rhs, dtype=float)
    radius = np.asarray(radius, dtype=float)
    core = radius <= min(20.0, float(radius[-1]) - 5.0)
    return {
        "max_abs_difference_shared_profile": float(np.max(np.abs(delta))),
        "rms_difference_shared_profile": float(np.sqrt(np.mean(delta * delta))),
        "max_abs_difference_r_le_20": float(np.max(np.abs(delta[core]))),
        "rms_difference_r_le_20": float(np.sqrt(np.mean(delta[core] * delta[core]))),
    }


def _state_snapshot(state):
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


def analyze_slice(kernel, state):
    full_grid = state.grid
    full_n = int(full_grid.n)
    seed_geometry = state.geometry.copy()
    # Native full-domain solve on this frozen slice; seed alpha is held fixed
    # while target_Kdot is projected, matching the native solver's call contract.
    alpha_full, kdot_full = solve_cmc_lapse(
        full_grid, seed_geometry, state.scalars, state.matter
    )
    full_radius = np.asarray(full_grid.centers, dtype=float)
    full_native_diff = alpha_full - np.asarray(state.geometry.alpha, dtype=float)

    crops_out = []
    for crop in CROPS:
        n = int(crop["N"])
        r_max = float(crop["r_max"])
        if n > full_n or abs(r_max / n - float(full_grid.dr)) > 1.0e-13:
            raise ValueError(f"invalid fixed-spacing crop specification: {crop}")
        grid = grid_ops.SphericalCellGrid(n, r_max)
        geom = _crop_geometry(state.geometry, full_n, n)
        scalars = _crop_scalars(state.scalars, n)
        matter = _crop_matter(state.matter, n)
        radius = np.asarray(grid.centers, dtype=float)

        # Native crop solve: crop-specific target_Kdot and alpha(R)=1.
        alpha_local, kdot_local = solve_cmc_lapse(grid, geom.copy(), scalars, matter)

        # Counterfactual A: hold target_Kdot at the full-domain projection,
        # while retaining the crop operator and native alpha(R)=1 boundary.
        alpha_full_kdot, diag_full_kdot = _solve_with_diagnostic_rhs(
            grid, geom.copy(), scalars, matter,
            kdot=float(kdot_full), outer_alpha=1.0,
        )

        # Counterfactual B: keep the same full-domain target_Kdot and the same
        # crop operator, but match alpha at the crop edge to the full-domain
        # solved value there. Comparing A and B isolates the Dirichlet-value
        # effect within one fixed cropped operator.
        boundary_match = float(alpha_full[n - 1])
        alpha_matched_boundary, diag_matched_boundary = _solve_with_diagnostic_rhs(
            grid, geom.copy(), scalars, matter,
            kdot=float(kdot_full), outer_alpha=boundary_match,
        )

        radius_common = full_radius[:n]
        local_vs_full = _profile_stats(radius, alpha_local, alpha_full[:n])
        fullkdot_vs_full = _profile_stats(radius, alpha_full_kdot, alpha_full[:n])
        boundary_vs_fullkdot = _profile_stats(
            radius, alpha_matched_boundary, alpha_full_kdot
        )
        local_vs_fullkdot = _profile_stats(
            radius, alpha_local, alpha_full_kdot
        )

        crops_out.append({
            "N": n,
            "r_max": r_max,
            "dr": float(grid.dr),
            "target_Kdot_full_domain": float(kdot_full),
            "target_Kdot_crop_native": float(kdot_local),
            "target_Kdot_crop_minus_full": float(kdot_local - kdot_full),
            "target_Kdot_crop_relative_difference": (
                float((kdot_local - kdot_full) / kdot_full)
                if abs(kdot_full) > 1.0e-300 else None
            ),
            "full_domain_solved_alpha_at_crop_edge": boundary_match,
            "native_crop_vs_full_resolve": local_vs_full,
            "fixed_full_Kdot_vs_full_resolve": fullkdot_vs_full,
            "boundary_value_effect_fixed_crop_operator": boundary_vs_fullkdot,
            "native_crop_vs_fixed_full_Kdot": local_vs_fullkdot,
            "native_crop_operator_residual": {
                "max_abs_discrete_operator_residual": float(
                    _solve_with_diagnostic_rhs(
                        grid, geom.copy(), scalars, matter,
                        kdot=float(kdot_local), outer_alpha=1.0,
                    )[1]["max_abs_discrete_operator_residual"]
                ),
            },
            "fixed_full_Kdot_operator": diag_full_kdot,
            "fixed_full_Kdot_boundary_matched_operator": diag_matched_boundary,
            "radius": [float(x) for x in radius],
            "alpha_native_crop": [float(x) for x in alpha_local],
            "alpha_fixed_full_Kdot_outer_alpha_1": [float(x) for x in alpha_full_kdot],
            "alpha_fixed_full_Kdot_boundary_matched": [float(x) for x in alpha_matched_boundary],
            "alpha_full_domain_resolved_restricted": [float(x) for x in alpha_full[:n]],
            "alpha_input_seed_restricted": [float(x) for x in state.geometry.alpha[:n]],
            "state_crop_uses_same_interior_arrays": True,
        })

    return {
        "t": float(state.t),
        "tau": float(state.tau),
        "e_folds": float(state.e_folds),
        "full_N": full_n,
        "full_r_max": float(FULL_R_MAX),
        "dr": float(full_grid.dr),
        "accepted_state_cycle_event_count": len(state.cycle.events),
        "accepted_state_handoff_count": len(state.handoffs),
        "target_Kdot_full_domain": float(kdot_full),
        "full_resolved_vs_stored_alpha": _profile_stats(
            full_radius, alpha_full, np.asarray(state.geometry.alpha, dtype=float)
        ),
        "radius_full": [float(x) for x in full_radius],
        "alpha_full_resolved": [float(x) for x in alpha_full],
        "alpha_full_stored": [float(x) for x in state.geometry.alpha],
        "crops": crops_out,
        "state_snapshot": _state_snapshot(state),
    }


def write_progress(state, accepted_steps, records, status, failure=None):
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema": "cmc_lapse_frozen_state_crop_audit_v1",
        "status": status,
        "diagnostic_only": True,
        "configuration": {
            "full_N": FULL_N,
            "full_r_max": FULL_R_MAX,
            "crop_domains": list(CROPS),
            "dr": 0.5,
            "CFL": CFL,
            "final_time": FINAL_TIME,
            "sample_interval": SAMPLE_INTERVAL,
            "amplitude": AMPLITUDE,
            "width": WIDTH,
            "D_amplitude": D_AMPLITUDE,
            "include_radiation": INCLUDE_RADIATION,
            "production_physics_changed": False,
            "production_defaults_changed": False,
            "gauge_changed": False,
            "projection_changed": False,
            "evolution_boundary_conditions_changed": False,
        },
        "accepted_steps": accepted_steps,
        "sample_count": len(records),
        "completed_sample_times": [float(row["t"]) for row in records],
        "failure": failure,
    }
    (OUT_DIR / "campaign_progress.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    profiles_path = OUT_DIR / "frozen_state_crop_profiles.jsonl"
    profiles_path.write_text("")
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=FULL_N,
        r_max=FULL_R_MAX,
        amplitude=AMPLITUDE,
        width=WIDTH,
        D_amplitude=D_AMPLITUDE,
        include_radiation=INCLUDE_RADIATION,
    )
    dr = float(state.grid.dr)
    dt_nominal = CFL * dr
    if abs(dr - 0.5) > 1.0e-13:
        raise ValueError(f"full-domain spacing changed unexpectedly: dr={dr}")

    records = []
    accepted_steps = 0
    failure = None
    initial = analyze_slice(kernel, state)
    records.append(initial)
    with profiles_path.open("a") as stream:
        stream.write(json.dumps(initial, sort_keys=True, allow_nan=False) + "\n")
        stream.flush()
    write_progress(state, accepted_steps, records, "in_progress")

    target_index = 0
    while state.t < FINAL_TIME - 1.0e-13:
        dt = min(dt_nominal, FINAL_TIME - state.t)
        try:
            next_state = kernel.step(state, dt)
            next_state.geometry.assert_finite_positive()
            for name, array in (
                ("S", next_state.scalars.S),
                ("D", next_state.scalars.D),
                ("phi", next_state.scalars.phi),
                ("a", next_state.geometry.a),
                ("b", next_state.geometry.b),
                ("X", next_state.geometry.X),
                ("Lambda", next_state.geometry.Lambda),
                ("alpha", next_state.geometry.alpha),
            ):
                arr = np.asarray(array, dtype=float)
                if not np.all(np.isfinite(arr)):
                    raise FloatingPointError(f"non-finite array in accepted state: {name}")
            state = next_state
            accepted_steps += 1
        except Exception as exc:
            failure = {
                "stage": "evolution_step",
                "step_attempt": accepted_steps + 1,
                "t_before": float(state.t),
                "error_type": type(exc).__name__,
                "error": str(exc),
                "traceback": traceback.format_exc(limit=12),
            }
            break

        while target_index < len(TARGET_TIMES) and (
            state.t + 0.5 * dt >= TARGET_TIMES[target_index]
            or abs(state.t - TARGET_TIMES[target_index]) < 1.0e-12
        ):
            row = analyze_slice(kernel, state)
            records.append(row)
            with profiles_path.open("a") as stream:
                stream.write(json.dumps(row, sort_keys=True, allow_nan=False) + "\n")
                stream.flush()
            target_index += 1
            write_progress(state, accepted_steps, records, "in_progress", failure)
            print(
                "[CMC_LAPSE_FROZEN_STATE_CROP_AUDIT] CHECKPOINT "
                + json.dumps({
                    "t": float(state.t),
                    "accepted_steps": accepted_steps,
                    "sample_count": len(records),
                    "target_Kdot_full_domain": float(row["target_Kdot_full_domain"]),
                    "kdot_crops": [
                        {
                            "N": c["N"],
                            "kdot": c["target_Kdot_crop_native"],
                            "max_alpha_diff": c["native_crop_vs_full_resolve"]["max_abs_difference_shared_profile"],
                        }
                        for c in row["crops"]
                    ],
                }, sort_keys=True),
                flush=True,
            )

    completed = failure is None and state.t >= FINAL_TIME - 1.0e-10
    status = "completed" if completed and len(records) == len(TARGET_TIMES) + 1 else "numerical_or_sampling_failure"
    write_progress(state, accepted_steps, records, status, failure)
    summary = {
        "schema": "cmc_lapse_frozen_state_crop_audit_v1",
        "status": status,
        "diagnostic_only": True,
        "accepted_steps": accepted_steps,
        "sample_count": len(records),
        "final_t": float(state.t),
        "final_tau": float(state.tau),
        "final_e_folds": float(state.e_folds),
        "cycle_event_count": len(state.cycle.events),
        "handoff_count": len(state.handoffs),
        "failure": failure,
        "profiles_file": profiles_path.name,
        "configuration": {
            "full_N": FULL_N,
            "full_r_max": FULL_R_MAX,
            "crop_domains": list(CROPS),
            "dr": dr,
            "CFL": CFL,
            "final_time": FINAL_TIME,
            "sample_interval": SAMPLE_INTERVAL,
            "production_physics_changed": False,
            "production_defaults_changed": False,
            "gauge_changed": False,
            "projection_changed": False,
            "evolution_boundary_conditions_changed": False,
        },
    }
    (OUT_DIR / "report.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    print("[CMC_LAPSE_FROZEN_STATE_CROP_AUDIT] FINAL "
          + json.dumps(summary, sort_keys=True), flush=True)
    if status != "completed":
        raise RuntimeError("Frozen-state CMC lapse crop audit did not complete cleanly.")


if __name__ == "__main__":
    main()
