"""Self-consistent diagnostic-only S-D bond sweep with consumer-closure audit.

No production file is modified. The interaction V_bond=g*S*D is a candidate
diagnostic term, not an archive-derived physical law.
"""
import csv
import json
import math
import os
from pathlib import Path
import sys

import numpy as np

# Executing this file by path puts its directory, not the repository root,
# at sys.path[0]. Make the project package import deterministic in CI.
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
CLOSURE_TOL = 5.0e-12


def maxabs(*arrays):
    return max((float(np.max(np.abs(np.asarray(a)))) for a in arrays), default=0.0)


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
    if G_BOND == 0.0:
        return rho, pr, pt, j
    Vbond = G_BOND * np.asarray(fields.S) * np.asarray(fields.D)
    # A potential contributes +V to rho and -V to each principal pressure.
    # It has no spatial gradient and therefore no momentum density.
    return rho + Vbond, pr - Vbond, pt - Vbond, j


# Geometry/lapse consumers use v55_matter.total_matter_projection; the
# Hamiltonian/momentum diagnostic uses stress_energy.assemble_total_stress_energy.
# Patch both module-level aliases so both paths see the same interaction.
pk.scalar_rhs_arrays = rhs_with_bond
v55_matter.scalar_projection = projection_with_bond
stress_energy.scalar_projection = projection_with_bond


def measure(state):
    _, vacuum, _ = adapter.vendor_modules()
    raw = vacuum.constraints(state.grid, state.geometry)
    total = stress_energy.assemble_total_stress_energy(
        state.grid, state.geometry, state.scalars, state.matter
    )
    through_geometry = v55_matter.total_matter_projection(
        state.grid, state.geometry, state.scalars, state.matter
    )

    # Independently compute the unbonded scalar projection and check that
    # both actual consumers receive exactly the expected +/- g*S*D term.
    base_rho, base_pr, base_pt, base_j = ORIGINAL_PROJECTION(
        state.grid, state.geometry, state.scalars
    )
    Vbond = G_BOND * np.asarray(state.scalars.S) * np.asarray(state.scalars.D)
    projection_error = maxabs(
        through_geometry["scalar"]["rho"] - (base_rho + Vbond),
        through_geometry["scalar"]["pr"] - (base_pr - Vbond),
        through_geometry["scalar"]["pt"] - (base_pt - Vbond),
        through_geometry["scalar"]["j"] - base_j,
        total.scalar_rho - (base_rho + Vbond),
        total.scalar_pr - (base_pr - Vbond),
        total.scalar_pt - (base_pt - Vbond),
        total.scalar_j - base_j,
    )
    consumer_error = maxabs(
        through_geometry["rho"] - total.rho,
        through_geometry["pr"] - total.pr,
        through_geometry["pt"] - total.pt,
        through_geometry["j"] - total.j,
    )
    if projection_error > CLOSURE_TOL or consumer_error > CLOSURE_TOL:
        raise RuntimeError(
            f"bond projection closure failed at t={state.t:.12g}: "
            f"projection_error={projection_error:.4e}, "
            f"consumer_error={consumer_error:.4e}"
        )

    H = np.asarray(raw["hamiltonian"]) - 16.0 * math.pi * np.asarray(total.rho)
    M = np.asarray(raw["momentum"]) - 8.0 * math.pi * np.asarray(total.j)
    r = np.asarray(state.grid.centers)
    alpha = np.asarray(state.geometry.alpha)
    S = np.asarray(state.scalars.S)
    D = np.asarray(state.scalars.D)
    mask = r < 10.0
    weights = r * r * float(state.grid.dr)
    return {
        "t": float(state.t),
        "H_center": float(H[0]),
        "H_offmax_excluding_first_two": float(np.max(np.abs(H[2:]))),
        "H_offmax_r": float(r[int(np.argmax(np.abs(H[2:])) + 2)]),
        "M_offmax_excluding_first_two": float(np.max(np.abs(M[2:]))),
        "lapse_center": float(alpha[0]),
        "one_minus_lapse_r3p25": float(1.0 - alpha[int(np.argmin(np.abs(r - 3.25)))]),
        "D_center": float(D[0]),
        "S_center": float(S[0]),
        "E_S_inside_r10": float(np.sum(
            (0.5 * (np.asarray(state.scalars.PS)**2 + S**2))[mask] * weights[mask]
        )),
        "E_D_kinetic_inside_r10": float(np.sum(
            (0.5 * np.asarray(state.scalars.PD)**2)[mask] * weights[mask]
        )),
        "bond_V_center": float(Vbond[0]),
        "lapse_min": float(np.min(alpha)),
        "projection_closure_max": projection_error,
        "geometry_vs_constraint_consumer_max": consumer_error,
    }


def log_row(row, stage):
    print(
        f"{stage} g={row['g']:+.3f} t={row['t']:.2f} "
        f"H0={row['H_center']:+.5e} "
        f"Hmax={row['H_offmax_excluding_first_two']:.5e} "
        f"Mmax={row['M_offmax_excluding_first_two']:.5e} "
        f"alpha0={row['lapse_center']:.7f} "
        f"1-alpha(r=3.25)={row['one_minus_lapse_r3p25']:+.5e} "
        f"D0={row['D_center']:+.5e} Vbond0={row['bond_V_center']:+.5e} "
        f"Perr={row['projection_closure_max']:.2e} "
        f"Cerr={row['geometry_vs_constraint_consumer_max']:.2e}",
        flush=True,
    )


def run_one(g):
    global G_BOND
    G_BOND = float(g)
    kernel = pk.V55ProductionKernel()
    state = kernel.initialize(
        resolution=40, r_max=40.0, amplitude=0.01, width=7.0,
        D_amplitude=1.0e-4, include_radiation=True,
    )
    dt_base = 0.0075 * float(state.grid.dr)
    rows = []
    initial = {"g": G_BOND, **measure(state)}
    rows.append(initial)
    log_row(initial, "INITIAL")
    for target in (2.0, 3.0, 4.0, 5.0):
        while state.t < target - 1e-12:
            state = kernel.step(state, min(dt_base, target - state.t))
        row = {"g": G_BOND, **measure(state)}
        rows.append(row)
        log_row(row, "SAMPLE")
    return rows


def main():
    all_rows = []
    for g in (0.0, 0.01, -0.01, 0.05, -0.05, 0.10):
        all_rows.extend(run_one(g))
    json_path = OUT / "bond_sweep_self_consistent.json"
    csv_path = OUT / "bond_sweep_self_consistent.csv"
    json_path.write_text(json.dumps(all_rows, indent=2, sort_keys=True) + "\n")
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(all_rows[0]))
        writer.writeheader()
        writer.writerows(all_rows)
    print(f"JSON={json_path} rows={len(all_rows)}")
    print(f"CSV={csv_path} rows={len(all_rows)}")
    print("PASS: all scalar-projection and cross-consumer closure checks passed")


if __name__ == "__main__":
    main()
