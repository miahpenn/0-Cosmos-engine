#!/usr/bin/env python3
"""Step 1 Hamiltonian-term recorder.

Diagnostic-only. This records the four already-defined Hamiltonian terms on
cells 0..4 and the corresponding constraint residual. It does not alter the
production evolution or add a source, repair, or fitted coefficient.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from engine.production_kernel import V55ProductionKernel
from engine.stress_energy import assemble_total_stress_energy
from engine import v55_pirk_adapter as adapter


class NoRepairKernel(V55ProductionKernel):
    """A/B diagnostic kernel: bypass only the existing center projection."""

    @staticmethod
    def _enforce_center_regularity(grid, geometry):
        return geometry


def hamiltonian_terms(state):
    """Return the existing four Hamiltonian terms without changing state."""
    grid = state.grid
    geom = state.geometry
    _, vacuum, _ = adapter.vendor_modules()
    total = assemble_total_stress_energy(
        grid, geom, state.scalars, state.matter
    )
    raw = vacuum.constraints(grid, geom)

    # These are exactly the four terms used by the admitted decomposition:
    # H = R - (Aa^2 + 2 Ab^2) + 2 K^2/3 - 16 pi rho.
    Ab = -0.5 * geom.Aa
    geometry = vacuum.geometry_terms(grid, geom)
    curvature = np.asarray(geometry["R"], dtype=float)
    extrinsic_A = -(geom.Aa**2 + 2.0 * Ab**2)
    extrinsic_K = (2.0 / 3.0) * geom.K**2
    matter_source = -16.0 * math.pi * total.rho

    reconstructed = curvature + extrinsic_A + extrinsic_K + matter_source
    H = np.asarray(raw["hamiltonian"], dtype=float) - 16.0 * math.pi * total.rho
    decomposition_error = reconstructed - H
    return {
        "r": np.asarray(grid.centers, dtype=float),
        "curvature": curvature,
        "extrinsic_A": np.asarray(extrinsic_A, dtype=float),
        "extrinsic_K": np.asarray(extrinsic_K, dtype=float),
        "matter_source": np.asarray(matter_source, dtype=float),
        "reconstructed": np.asarray(reconstructed, dtype=float),
        "vendor_H": np.asarray(H, dtype=float),
        "decomposition_error": np.asarray(decomposition_error, dtype=float),
    }


def record_sample(state, sample_index):
    terms = hamiltonian_terms(state)
    r = terms["r"]
    H = terms["vendor_H"]
    cells = []
    for i in range(5):
        cells.append({
            "cell": i,
            "r": float(r[i]),
            "curvature": float(terms["curvature"][i]),
            "abs_curvature": float(abs(terms["curvature"][i])),
            "extrinsic_A": float(terms["extrinsic_A"][i]),
            "abs_extrinsic_A": float(abs(terms["extrinsic_A"][i])),
            "extrinsic_K": float(terms["extrinsic_K"][i]),
            "abs_extrinsic_K": float(abs(terms["extrinsic_K"][i])),
            "matter_source": float(terms["matter_source"][i]),
            "abs_matter_source": float(abs(terms["matter_source"][i])),
            "reconstructed": float(terms["reconstructed"][i]),
            "vendor_H": float(H[i]),
            "decomposition_error": float(terms["decomposition_error"][i]),
        })

    masked = H[2:]
    i_global = int(np.argmax(np.abs(H)))
    i_masked = 2 + int(np.argmax(np.abs(masked)))
    return {
        "sample_index": int(sample_index),
        "t": float(state.t),
        "cells_0_to_4": cells,
        "global_max_abs_H": float(np.max(np.abs(H))),
        "global_max_cell": i_global,
        "global_max_r": float(r[i_global]),
        "masked_max_abs_H": float(np.max(np.abs(masked))),
        "masked_max_cell": i_masked,
        "masked_max_r": float(r[i_masked]),
        "decomposition_error_max": float(
            np.max(np.abs(terms["decomposition_error"]))
        ),
    }


def run_case(
    *,
    D=1.0e-4,
    repair=True,
    N=160,
    r_max=80.0,
    cfl=0.0075,
    final_time=24.0,
    sample_times=(0.0, 6.0, 12.0, 18.0, 24.0),
):
    """Run one repair state and record only at the registered sample times."""
    if N < 5:
        raise ValueError("N must be at least 5")
    Kernel = V55ProductionKernel if repair else NoRepairKernel
    kernel = Kernel()
    state = kernel.initialize(
        resolution=N,
        r_max=r_max,
        D_amplitude=D,
        include_radiation=True,
    )
    dt_nominal = cfl * state.grid.dr
    targets = tuple(float(x) for x in sample_times)
    if not targets or targets[0] != 0.0:
        raise ValueError("sample_times must begin with 0.0")

    samples = [record_sample(state, 0)]
    target_index = 1
    failed = None

    while state.t < final_time - 1.0e-14:
        dt = min(dt_nominal, final_time - state.t)
        try:
            state = kernel.step(state, dt)
            state.geometry.assert_finite_positive()
        except (FloatingPointError, ValueError) as exc:
            failed = {"t": float(state.t), "error": str(exc)}
            break

        if (
            target_index < len(targets)
            and state.t + 0.5 * dt >= targets[target_index]
        ):
            samples.append(record_sample(state, target_index))
            target_index += 1

    return {
        "status": "numerical_failure" if failed else "completed",
        "failure": failed,
        "repair": bool(repair),
        "D_amplitude": float(D),
        "resolution": int(N),
        "r_max": float(r_max),
        "cfl": float(cfl),
        "requested_final_time": float(final_time),
        "final_time": float(state.t),
        "samples": samples,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repair", choices=("on", "off"), default="on")
    parser.add_argument("--D", type=float, default=1.0e-4)
    parser.add_argument("--N", type=int, default=160)
    parser.add_argument("--r-max", type=float, default=80.0)
    parser.add_argument("--cfl", type=float, default=0.0075)
    parser.add_argument("--final-time", type=float, default=24.0)
    parser.add_argument(
        "--output",
        default="runs/0star-step1/hamiltonian_terms.json",
    )
    args = parser.parse_args()

    result = run_case(
        D=args.D,
        repair=args.repair == "on",
        N=args.N,
        r_max=args.r_max,
        cfl=args.cfl,
        final_time=args.final_time,
    )
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
