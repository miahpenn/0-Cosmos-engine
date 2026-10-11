"""Self-consistent diagnostic-only S-D bond sweep.

This does not edit production physics. In-process wrappers add V_bond=g*S*D
to scalar momentum RHS and its potential stress projection. Run at N=40,
Rmax=40, CFL=.0075, t=5 for g=0,+/-.01,+/-.05,+.10.
"""
import csv
import json
import math
import os
from pathlib import Path
import sys

import numpy as np

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


# Patch every consumer of scalar_projection in this process. v55_matter's
# total_matter_projection is used by the geometry RHS and CMC lapse target;
# stress_energy's assembler is used by the Hamiltonian/momentum diagnostics.
pk.scalar_rhs_arrays = rhs_with_bond
v55_matter.scalar_projection = projection_with_bond
stress_energy.scalar_projection = projection_with_bond


def measure(state):
    _, vacuum, _ = adapter.vendor_modules()
    raw = vacuum.constraints(state.grid, state.geometry)
    total = stress_energy.assemble_total_stress_energy(
        state.grid, state.geometry, state.scalars, state.matter
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
        "bond_V_center": float(G_BOND * S[0] * D[0]),
        "lapse_min": float(np.min(alpha)),
    }


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
    # Include the initial slice and common target times.
    rows.append({"g": G_BOND, **measure(state)})
    for target in (2.0, 3.0, 4.0, 5.0):
        while state.t < target - 1e-12:
            state = kernel.step(state, min(dt_base, target - state.t))
        row = {"g": G_BOND, **measure(state)}
        rows.append(row)
        print(
            f"g={G_BOND:+.3f} t={row['t']:.2f} "
            f"H0={row['H_center']:+.5e} Hmax={row['H_offmax_excluding_first_two']:.5e} "
            f"Mmax={row['M_offmax_excluding_first_two']:.5e} "
            f"alpha0={row['lapse_center']:.7f} "
            f"1-alpha(r=3.25)={row['one_minus_lapse_r3p25']:+.5e} "
            f"D0={row['D_center']:+.5e} Vbond0={row['bond_V_center']:+.5e}",
            flush=True,
        )
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


if __name__ == "__main__":
    main()
