"""Initial-slice central Hamiltonian residual audit (diagnostic only).

Matched N=40/N=80 initial data; no evolution and no physics changes.
The audit reports the selected H and momentum residuals using the same
definitions as diagnostics/bond_sweep/bond_sweep.py, with g=0.
"""
import csv
import json
import math
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import engine.production_kernel as pk
import engine.stress_energy as stress_energy
from engine import v55_pirk_adapter as adapter

OUT = Path("runs/central_residual_audit")
OUT.mkdir(parents=True, exist_ok=True)


def maxabs(a):
    a = np.asarray(a)
    if not np.all(np.isfinite(a)):
        return math.inf
    return float(np.max(np.abs(a))) if a.size else 0.0


def audit_resolution(n):
    kernel = pk.V55ProductionKernel()
    state = kernel.initialize(
        resolution=n,
        r_max=40.0,
        amplitude=0.01,
        width=7.0,
        D_amplitude=1.0e-4,
        include_radiation=True,
    )
    _, vacuum, _ = adapter.vendor_modules()
    raw = vacuum.constraints(state.grid, state.geometry)
    total = stress_energy.assemble_total_stress_energy(
        state.grid, state.geometry, state.scalars, state.matter
    )
    H = np.asarray(raw["hamiltonian"]) - 16.0 * math.pi * np.asarray(total.rho)
    M = np.asarray(raw["momentum"]) - 8.0 * math.pi * np.asarray(total.j)
    r = np.asarray(state.grid.centers)
    alpha = np.asarray(state.geometry.alpha)
    row = {
        "resolution": int(n),
        "r_max": 40.0,
        "t": float(state.t),
        "dr": float(state.grid.dr),
        "central_H_residual": float(H[0]),
        "central_H_abs": float(abs(H[0])),
        "H_max_excluding_first_two": maxabs(H[2:]),
        "H_max_excluding_first_two_r": float(r[int(np.argmax(np.abs(H[2:])) + 2)]),
        "central_M_residual": float(M[0]),
        "M_max_excluding_first_two": maxabs(M[2:]),
        "central_lapse": float(alpha[0]),
        "central_D": float(np.asarray(state.scalars.D)[0]),
        "central_S": float(np.asarray(state.scalars.S)[0]),
        "finite": bool(np.all(np.isfinite(H)) and np.all(np.isfinite(M))),
    }
    if not row["finite"]:
        raise RuntimeError(f"non-finite constraint residual at N={n}")
    print(
        f"AUDIT N={n} t={state.t:.2f} dr={state.grid.dr:.8g} "
        f"H0={H[0]:+.8e} Hmax_off={row['H_max_excluding_first_two']:.8e} "
        f"M0={M[0]:+.8e} Mmax_off={row['M_max_excluding_first_two']:.8e} "
        f"alpha0={alpha[0]:.8f}",
        flush=True,
    )
    return row


def main():
    rows = [audit_resolution(n) for n in (40, 80)]
    baseline = rows[0]
    refined = rows[1]
    comparison = {
        "N40_to_N80_central_H_abs_ratio": (
            refined["central_H_abs"] / baseline["central_H_abs"]
            if baseline["central_H_abs"] else None
        ),
        "N40_to_N80_offcenter_Hmax_ratio": (
            refined["H_max_excluding_first_two"] / baseline["H_max_excluding_first_two"]
            if baseline["H_max_excluding_first_two"] else None
        ),
        "N40_to_N80_central_H_sign_change": (
            np.sign(refined["central_H_residual"]) != np.sign(baseline["central_H_residual"])
        ),
        "interpretation": (
            "Two-resolution initial-slice comparison only; this is not a "
            "convergence claim and does not identify the defect source."
        ),
    }
    (OUT / "central_residual_audit.json").write_text(
        json.dumps({"rows": rows, "comparison": comparison}, indent=2, sort_keys=True) + "\n"
    )
    with (OUT / "central_residual_audit.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print("PASS: finite initial-slice residuals recorded for N=40 and N=80")
    print("COMPARISON " + json.dumps(comparison, sort_keys=True))


if __name__ == "__main__":
    main()
