#!/usr/bin/env python3
"""Recompute signed constraint profiles from archived production checkpoints.

This is a post-processing diagnostic. It does not evolve the state or change
any production physics. The formulas are the production source assembly and
pinned vendor constraint operators applied to saved state arrays.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np


GEOMETRY_FIELDS = ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B")
SCALAR_FIELDS = ("S", "PS", "D", "PD", "phi", "Pi")
MATTER_FIELDS = {
    "dark_matter": "dm",
    "baryons": "baryon",
    "radiation": "radiation",
}
CONSERVED_FIELDS = ("rest", "energy_t", "momentum_r")
R_CUT = 5.0
SAMPLE_RADII = (0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 15.0, 20.0, 40.0, 60.0, 75.0)


def _peak(values: np.ndarray, radii: np.ndarray, mask: np.ndarray) -> dict:
    indices = np.flatnonzero(mask)
    if indices.size == 0:
        return {"max_abs": None, "radius": None, "signed_value": None}
    idx = int(indices[np.argmax(np.abs(values[indices]))])
    return {
        "max_abs": float(abs(values[idx])),
        "radius": float(radii[idx]),
        "signed_value": float(values[idx]),
        "index": idx,
    }


def _restore_checkpoint(path: Path, kernel):
    from engine.production_kernel import ProductionState

    with np.load(path, allow_pickle=False) as archive:
        data = {name: archive[name].copy() for name in archive.files}

    radii_saved = np.asarray(data["r"], dtype=float)
    n = int(radii_saved.size)
    if n < 4:
        raise ValueError(f"{path}: expected at least 4 grid cells, got {n}")
    dr = float(radii_saved[1] - radii_saved[0])
    r_max = float(n * dr)
    state = kernel.initialize(
        resolution=n,
        r_max=r_max,
        include_radiation=True,
    )
    radii = np.asarray(state.grid.centers, dtype=float)
    if radii.shape != radii_saved.shape or not np.allclose(
        radii, radii_saved, rtol=0.0, atol=2.0e-13 * max(1.0, r_max)
    ):
        raise ValueError(f"{path}: archived grid does not match reconstructed grid")

    for name in GEOMETRY_FIELDS:
        setattr(state.geometry, name, np.asarray(data[name], dtype=float).copy())
    for name in SCALAR_FIELDS:
        setattr(state.scalars, name, np.asarray(data[name], dtype=float).copy())
    for species_name, prefix in MATTER_FIELDS.items():
        species_state = getattr(state.matter, species_name)
        for field_name in CONSERVED_FIELDS:
            setattr(
                species_state,
                field_name,
                np.asarray(data[f"{prefix}_{field_name}"], dtype=float).copy(),
            )

    state.t = float(data["t"])
    state.tau = float(data["tau"])
    state.e_folds = float(data["e_folds"])
    state.history = []
    state.handoffs = []
    return state


def audit_checkpoint(path: Path, kernel, vacuum):
    from engine.stress_energy import assemble_total_stress_energy

    state = _restore_checkpoint(path, kernel)
    radii = np.asarray(state.grid.centers, dtype=float)
    total = assemble_total_stress_energy(
        state.grid, state.geometry, state.scalars, state.matter
    )
    raw = vacuum.constraints(state.grid, state.geometry)
    profiles = {
        "hamiltonian": np.asarray(raw["hamiltonian"], dtype=float)
        - 16.0 * np.pi * np.asarray(total.rho, dtype=float),
        "momentum": np.asarray(raw["momentum"], dtype=float)
        - 8.0 * np.pi * np.asarray(total.j, dtype=float),
        "connection": np.asarray(raw["connection"], dtype=float),
        "determinant_error": (
            np.asarray(state.geometry.a, dtype=float)
            * np.asarray(state.geometry.b, dtype=float) ** 2 - 1.0
        ),
    }
    for name, values in profiles.items():
        if values.shape != radii.shape:
            raise ValueError(f"{path}: {name} profile has unexpected shape {values.shape}")
        if not np.all(np.isfinite(values)):
            raise FloatingPointError(f"{path}: non-finite values in {name} profile")

    inner = radii < R_CUT
    exterior = radii >= R_CUT
    summary = {
        "resolution": int(radii.size),
        "t": float(state.t),
        "tau": float(state.tau),
        "r_cut": R_CUT,
        "peaks": {
            name: {
                "inner_r_lt_5": _peak(values, radii, inner),
                "exterior_r_ge_5": _peak(values, radii, exterior),
                "full_domain": _peak(values, radii, np.ones_like(radii, dtype=bool)),
            }
            for name, values in profiles.items()
        },
        "sampled_signed_profiles": {},
    }
    for target_radius in SAMPLE_RADII:
        idx = int(np.argmin(np.abs(radii - target_radius)))
        summary["sampled_signed_profiles"][str(target_radius)] = {
            "grid_radius": float(radii[idx]),
            "values": {
                name: float(values[idx]) for name, values in profiles.items()
            },
        }
    key = f"N{radii.size}_t{state.t:06.2f}".replace(".", "p")
    arrays = {"r": radii, "t": np.asarray(state.t), "tau": np.asarray(state.tau)}
    arrays.update(profiles)
    return summary, key, arrays


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--run-id", type=int, required=True)
    parser.add_argument("--source-sha", required=True)
    args = parser.parse_args()

    sys.path.insert(0, str(args.source.resolve()))
    from engine.production_kernel import V55ProductionKernel
    from engine import v55_pirk_adapter as adapter

    _, vacuum, _ = adapter.vendor_modules()
    kernel = V55ProductionKernel()
    checkpoints = sorted(args.input.glob("N*/checkpoint_t*.npz"))
    if not checkpoints:
        raise FileNotFoundError(f"no checkpoint files found under {args.input}")

    summaries = []
    profile_arrays = {}
    for path in checkpoints:
        summary, key, arrays = audit_checkpoint(path, kernel, vacuum)
        summaries.append(summary)
        for name, values in arrays.items():
            profile_arrays[f"{key}_{name}"] = values

    # Ensure the archived campaign is the expected synchronized checkpoint set.
    expected = {(n, float(t)) for n in (40, 60, 80) for t in (2, 4, 6, 8, 10, 12)}
    observed = {(int(item["resolution"]), float(item["t"])) for item in summaries}
    if observed != expected:
        raise ValueError(
            f"checkpoint set mismatch; missing={sorted(expected-observed)}, "
            f"unexpected={sorted(observed-expected)}"
        )

    args.output.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.output / "signed_constraint_profiles.npz", **profile_arrays)
    report = {
        "status": "PROFILE_AUDIT_COMPLETE",
        "classification": "postprocessing_only_not_physical_admission",
        "production_run_id": args.run_id,
        "production_source_sha": args.source_sha,
        "r_cut": R_CUT,
        "checkpoint_count": len(summaries),
        "resolution_count": len({int(item["resolution"]) for item in summaries}),
        "summaries": sorted(summaries, key=lambda item: (item["resolution"], item["t"])),
        "limitations": [
            "Profiles are recomputed with the production stress assembly and pinned vendor constraint operators.",
            "This audits spatial localization but is not an independent re-derivation of the continuum Einstein constraints.",
            "No evolution is performed and no physics or boundary condition is modified.",
        ],
    }
    (args.output / "profile_audit.json").write_text(json.dumps(report, indent=2))
    print(json.dumps({
        "status": report["status"],
        "production_run_id": args.run_id,
        "checkpoint_count": len(summaries),
        "output": str(args.output),
    }, indent=2))


if __name__ == "__main__":
    main()
