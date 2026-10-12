"""Paired provisional-vs-self-consistent S-D bond diagnostic.

Both models use the same input parameters and scalar profile. Each is initialized
under its own stress prescription. The provisional mode adds gSD only to scalar
momentum RHS (intentionally missing bond stress); the self-consistent mode also
adds its stress projection to geometry/CMC and constraint consumers. On every
state, both H-residual definitions are postprocessed using that exact same slice.
This is diagnostic only, not a physics claim.
"""
import csv
import json
import math
import os
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import engine.production_kernel as pk
import engine.scalar_system as scalar_system
import engine.stress_energy as stress_energy
import engine.v55_matter as v55_matter
from engine import v55_pirk_adapter as adapter

OUT = Path(os.environ.get("BOND_SWEEP_OUT", "runs/bond_sweep"))
OUT.mkdir(parents=True, exist_ok=True)

ORIGINAL_RHS = pk.scalar_rhs_arrays
ORIGINAL_PROJECTION = scalar_system.scalar_projection
G_BOND = 0.0
STRESS_BOND_ENABLED = False
CLOSURE_TOL = 5.0e-12


def maxabs(*arrays):
    values = []
    for arr in arrays:
        a = np.asarray(arr)
        if not np.all(np.isfinite(a)):
            return math.inf
        values.append(float(np.max(np.abs(a))) if a.size else 0.0)
    return max(values, default=0.0)


def rhs_with_bond(grid, geometry, fields, *args, **kwargs):
    out = ORIGINAL_RHS(grid, geometry, fields, *args, **kwargs)
    if G_BOND == 0.0:
        return out
    import dataclasses
    return dataclasses.replace(
        out,
        PS=out.PS - geometry.alpha * G_BOND * fields.D,
        PD=out.PD - geometry.alpha * G_BOND * fields.S,
    )


def projection_with_bond(grid, geometry, fields):
    rho, pr, pt, j = ORIGINAL_PROJECTION(grid, geometry, fields)
    if G_BOND == 0.0 or not STRESS_BOND_ENABLED:
        return rho, pr, pt, j
    Vbond = G_BOND * np.asarray(fields.S) * np.asarray(fields.D)
    return rho + Vbond, pr - Vbond, pt - Vbond, j


# These are the two independently imported projection aliases used by the
# geometry/CMC source and the constraint-diagnostic stress-energy assembler.
pk.scalar_rhs_arrays = rhs_with_bond
v55_matter.scalar_projection = projection_with_bond
stress_energy.scalar_projection = projection_with_bond


def measure(state, model):
    _, vacuum, _ = adapter.vendor_modules()
    raw = vacuum.constraints(state.grid, state.geometry)
    total = stress_energy.assemble_total_stress_energy(
        state.grid, state.geometry, state.scalars, state.matter
    )
    geometry_total = v55_matter.total_matter_projection(
        state.grid, state.geometry, state.scalars, state.matter
    )
    base_rho, base_pr, base_pt, base_j = ORIGINAL_PROJECTION(
        state.grid, state.geometry, state.scalars
    )
    S = np.asarray(state.scalars.S)
    D = np.asarray(state.scalars.D)
    full_V = G_BOND * S * D
    applied_V = full_V if STRESS_BOND_ENABLED else np.zeros_like(full_V)

    # Expected potential deltas and agreement between both consumers.
    projection_error = maxabs(
        geometry_total["scalar"]["rho"] - (base_rho + applied_V),
        geometry_total["scalar"]["pr"] - (base_pr - applied_V),
        geometry_total["scalar"]["pt"] - (base_pt - applied_V),
        geometry_total["scalar"]["j"] - base_j,
        total.scalar_rho - (base_rho + applied_V),
        total.scalar_pr - (base_pr - applied_V),
        total.scalar_pt - (base_pt - applied_V),
        total.scalar_j - base_j,
    )
    consumer_error = maxabs(
        geometry_total["rho"] - total.rho,
        geometry_total["pr"] - total.pr,
        geometry_total["pt"] - total.pt,
        geometry_total["j"] - total.j,
    )
    if projection_error > CLOSURE_TOL or consumer_error > CLOSURE_TOL:
        raise RuntimeError(
            f"projection closure failed model={model} t={state.t:.12g}: "
            f"projection_error={projection_error:.4e}, "
            f"consumer_error={consumer_error:.4e}"
        )

    # Same-slice comparison: compute both H definitions on this exact state.
    # Remove whatever bond term is active to recover the common baseline rho,
    # then reapply it algebraically. This separates definition effects from
    # trajectory/evolution effects.
    rho_base = np.asarray(total.rho) - applied_V
    H_raw = np.asarray(raw["hamiltonian"])
    H_no_bond = H_raw - 16.0 * math.pi * rho_base
    H_with_bond = H_raw - 16.0 * math.pi * (rho_base + full_V)
    H_selected = H_with_bond if STRESS_BOND_ENABLED else H_no_bond
    M = np.asarray(raw["momentum"]) - 8.0 * math.pi * np.asarray(total.j)

    r = np.asarray(state.grid.centers)
    alpha = np.asarray(state.geometry.alpha)
    mask = r < 10.0
    weights = r * r * float(state.grid.dr)
    idx325 = int(np.argmin(np.abs(r - 3.25)))
    off_idx = int(np.argmax(np.abs(H_selected[2:])) + 2)
    row = {
        "model": model,
        "g": float(G_BOND),
        "t": float(state.t),
        "H_center": float(H_selected[0]),
        "H_center_no_bond_same_state": float(H_no_bond[0]),
        "H_center_with_bond_same_state": float(H_with_bond[0]),
        "H_center_bond_correction_same_state": float(H_with_bond[0] - H_no_bond[0]),
        "H_offmax_excluding_first_two": float(np.max(np.abs(H_selected[2:]))),
        "H_offmax_r": float(r[off_idx]),
        "H_offmax_no_bond_same_state": float(np.max(np.abs(H_no_bond[2:]))),
        "H_offmax_with_bond_same_state": float(np.max(np.abs(H_with_bond[2:]))),
        "M_offmax_excluding_first_two": float(np.max(np.abs(M[2:]))),
        "lapse_center": float(alpha[0]),
        "one_minus_lapse_r3p25": float(1.0 - alpha[idx325]),
        "D_center": float(D[0]),
        "S_center": float(S[0]),
        "bond_V_center": float(full_V[0]),
        "stress_bond_enabled": bool(STRESS_BOND_ENABLED),
        "lapse_min": float(np.min(alpha)),
        "projection_closure_max": projection_error,
        "geometry_vs_constraint_consumer_max": consumer_error,
        "E_S_inside_r10": float(np.sum(
            (0.5 * (np.asarray(state.scalars.PS)**2 + S**2))[mask] * weights[mask]
        )),
        "E_D_kinetic_inside_r10": float(np.sum(
            (0.5 * np.asarray(state.scalars.PD)**2)[mask] * weights[mask]
        )),
    }
    return row


def log_row(row, stage):
    print(
        f"{stage} model={row['model']} g={row['g']:+.3f} t={row['t']:.2f} "
        f"H0={row['H_center']:+.5e} "
        f"H0_noV_same={row['H_center_no_bond_same_state']:+.5e} "
        f"H0_withV_same={row['H_center_with_bond_same_state']:+.5e} "
        f"Hmax={row['H_offmax_excluding_first_two']:.5e} "
        f"Mmax={row['M_offmax_excluding_first_two']:.5e} "
        f"alpha0={row['lapse_center']:.7f} D0={row['D_center']:+.5e} "
        f"V0={row['bond_V_center']:+.5e} "
        f"Perr={row['projection_closure_max']:.2e} "
        f"Cerr={row['geometry_vs_constraint_consumer_max']:.2e}",
        flush=True,
    )


def run_model(g, model):
    global G_BOND, STRESS_BOND_ENABLED
    G_BOND = float(g)
    STRESS_BOND_ENABLED = model == "self_consistent"
    kernel = pk.V55ProductionKernel()
    # Initialize in the selected mode so the self-consistent lapse/geometry
    # source includes the potential from the first slice. Input scalar profiles
    # and numerical settings remain the same across modes.
    state = kernel.initialize(
        resolution=40, r_max=40.0, amplitude=0.01, width=7.0,
        D_amplitude=1.0e-4, include_radiation=True,
    )
    dt_base = 0.0075 * float(state.grid.dr)
    rows = []
    initial = measure(state, model)
    rows.append(initial)
    log_row(initial, "INITIAL")
    for target in (2.0, 3.0, 4.0, 5.0):
        while state.t < target - 1e-12:
            state = kernel.step(state, min(dt_base, target - state.t))
        row = measure(state, model)
        rows.append(row)
        log_row(row, "SAMPLE")
    return rows


def main():
    global G_BOND, STRESS_BOND_ENABLED
    rows = []
    couplings = (0.0, 0.05, -0.05, 0.10)
    for g in couplings:
        for model in ("rhs_only", "self_consistent"):
            rows.extend(run_model(g, model))

    # Pair trajectories by coupling/time; residual definitions on each same
    # slice are separately provided in the main rows above.
    keyed = {(r["g"], r["model"], round(r["t"], 8)): r for r in rows}
    deltas = []
    for g in couplings:
        for t in (0.0, 2.0, 3.0, 4.0, 5.0):
            p = keyed[(float(g), "rhs_only", round(float(t), 8))]
            c = keyed[(float(g), "self_consistent", round(float(t), 8))]
            deltas.append({
                "g": float(g), "t": float(t),
                "delta_D_center_consistent_minus_rhs_only": c["D_center"] - p["D_center"],
                "delta_lapse_center_consistent_minus_rhs_only": c["lapse_center"] - p["lapse_center"],
                "delta_H_selected_center_consistent_minus_rhs_only": c["H_center"] - p["H_center"],
                "rhs_only_H0_no_bond": p["H_center_no_bond_same_state"],
                "rhs_only_H0_with_bond_posthoc": p["H_center_with_bond_same_state"],
                "self_consistent_H0_no_bond_posthoc": c["H_center_no_bond_same_state"],
                "self_consistent_H0_with_bond": c["H_center_with_bond_same_state"],
                "rhs_only_projection_closure_max": p["projection_closure_max"],
                "self_consistent_projection_closure_max": c["projection_closure_max"],
            })

    for name, data in (
        ("bond_sweep_mode_pair", rows),
        ("bond_sweep_mode_pair_deltas", deltas),
    ):
        (OUT / f"{name}.json").write_text(
            json.dumps(data, indent=2, sort_keys=True) + "\n"
        )
        with (OUT / f"{name}.csv").open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(data[0]))
            writer.writeheader()
            writer.writerows(data)
    print(f"PASS: projection closure passed for {len(rows)} rows; "
          f"same-state H definitions and paired trajectories recorded")
    print(f"outputs={OUT} rows={len(rows)} pairs={len(deltas)}")


if __name__ == "__main__":
    main()
