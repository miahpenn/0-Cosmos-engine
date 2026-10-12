#!/usr/bin/env python3
"""Offline Hamiltonian/momentum source budget for an admitted V5.5 run.

This script does not call kernel.step(). It evaluates the t=0 initialization and
reads the already-saved production checkpoints, then compares the geometric
constraint terms with the total stress-energy source on each exact slice.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from engine import v55_pirk_adapter as adapter
from engine.matter_system import ConservedSpecies
from engine.production_kernel import V55ProductionKernel
from engine.scalar_system import ScalarFields
from engine.stress_energy import assemble_total_stress_energy
from engine.v55_matter import V55MatterState


GEOMETRY_FIELDS = ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B")
SCALAR_FIELDS = ("S", "PS", "D", "PD", "phi", "Pi")
MATTER_FIELDS = {
    "dark_matter": ("dm_rest", "dm_energy_t", "dm_momentum_r"),
    "baryons": ("baryon_rest", "baryon_energy_t", "baryon_momentum_r"),
    "radiation": ("radiation_rest", "radiation_energy_t", "radiation_momentum_r"),
}
SNAPSHOT_TIMES = (0, 2, 4, 6, 8, 10, 12)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def scalar_from_npz(data, name: str) -> float:
    value = np.asarray(data[name])
    if value.size != 1:
        raise ValueError(f"checkpoint field {name!r} is not scalar")
    return float(value.reshape(()))


def restore_checkpoint(kernel: V55ProductionKernel, resolution: int, r_max: float, path: Path):
    with np.load(path, allow_pickle=False) as archive:
        data = {key: archive[key].copy() for key in archive.files}

    state = kernel.initialize(
        resolution=resolution,
        r_max=r_max,
        amplitude=0.01,
        width=7.0,
        D_amplitude=1.0e-10,
        include_radiation=True,
    )
    archived_r = np.asarray(data["r"], dtype=float)
    actual_r = np.asarray(state.grid.centers, dtype=float)
    if archived_r.shape != actual_r.shape or not np.array_equal(archived_r, actual_r):
        raise ValueError(f"{path}: checkpoint grid does not match N={resolution}, r_max={r_max}")

    for field in GEOMETRY_FIELDS:
        setattr(state.geometry, field, np.asarray(data[field], dtype=float).copy())
    state.scalars = ScalarFields(*(np.asarray(data[field], dtype=float).copy() for field in SCALAR_FIELDS))

    species = {}
    for name, keys in MATTER_FIELDS.items():
        species[name] = ConservedSpecies(
            *(np.asarray(data[key], dtype=float).copy() for key in keys)
        )
    state.matter = V55MatterState(
        dark_matter=species["dark_matter"],
        baryons=species["baryons"],
        radiation=species["radiation"],
    )
    state.t = scalar_from_npz(data, "t")
    state.tau = scalar_from_npz(data, "tau")
    state.e_folds = scalar_from_npz(data, "e_folds")
    return state


def witness(values: np.ndarray, radii: np.ndarray) -> dict:
    index = int(np.argmax(np.abs(values)))
    return {
        "absmax": float(abs(values[index])),
        "signed_value": float(values[index]),
        "radius": float(radii[index]),
        "index": index,
    }


def witness_mask(values: np.ndarray, radii: np.ndarray, mask: np.ndarray) -> dict:
    selected = np.flatnonzero(np.asarray(mask, dtype=bool))
    if selected.size == 0:
        raise ValueError("physical-radius mask selects no grid points")
    index = int(selected[np.argmax(np.abs(values[selected]))])
    return {
        "absmax": float(abs(values[index])),
        "signed_value": float(values[index]),
        "radius": float(radii[index]),
        "index": index,
    }


def finite_arrays(**arrays: np.ndarray) -> None:
    for name, values in arrays.items():
        if not np.all(np.isfinite(values)):
            raise FloatingPointError(f"non-finite profile data in {name}")


def audit_state(state, vacuum) -> tuple[dict, dict]:
    grid, geometry = state.grid, state.geometry
    total = assemble_total_stress_energy(grid, geometry, state.scalars, state.matter)
    raw = vacuum.constraints(grid, geometry)
    r = np.asarray(grid.centers, dtype=float)

    h_geometry = np.asarray(raw["hamiltonian"], dtype=float)
    h_source_scalar = 16.0 * math.pi * np.asarray(total.scalar_rho, dtype=float)
    h_source_fluid = 16.0 * math.pi * np.asarray(total.fluid_rho, dtype=float)
    h_source = h_source_scalar + h_source_fluid
    h_residual = h_geometry - h_source

    m_geometry = np.asarray(raw["momentum"], dtype=float)
    m_source_scalar = 8.0 * math.pi * np.asarray(total.scalar_j, dtype=float)
    m_source_fluid = 8.0 * math.pi * np.asarray(total.fluid_j, dtype=float)
    m_source = m_source_scalar + m_source_fluid
    m_residual = m_geometry - m_source

    arrays = {
        "h_geometry": h_geometry,
        "h_source_scalar": h_source_scalar,
        "h_source_fluid": h_source_fluid,
        "h_source": h_source,
        "h_residual": h_residual,
        "m_geometry": m_geometry,
        "m_source_scalar": m_source_scalar,
        "m_source_fluid": m_source_fluid,
        "m_source": m_source,
        "m_residual": m_residual,
    }
    finite_arrays(**arrays)

    inner_mask = r < 5.0
    exterior_mask = r >= 5.0
    h_full = witness(h_residual, r)
    h_inner = witness_mask(h_residual, r, inner_mask)
    h_exterior = witness_mask(h_residual, r, exterior_mask)
    m_full = witness(m_residual, r)
    m_inner = witness_mask(m_residual, r, inner_mask)
    m_exterior = witness_mask(m_residual, r, exterior_mask)
    summary = {
        "resolution": int(grid.n),
        "t": float(state.t),
        "tau": float(state.tau),
        "dr": float(grid.dr),
        "fixed_inner_exclusion_radius": 5.0,
        "finite": True,
        "hamiltonian": {
            "full": h_full,
            "inside_fixed_radius_r_lt_5": h_inner,
            "outside_fixed_radius_r_ge_5": h_exterior,
            "geometry_at_full_residual_peak": float(h_geometry[h_full["index"]]),
            "source_at_full_residual_peak": float(h_source[h_full["index"]]),
            "source_scalar_at_full_residual_peak": float(h_source_scalar[h_full["index"]]),
            "source_fluid_at_full_residual_peak": float(h_source_fluid[h_full["index"]]),
            "geometry_absmax": witness(h_geometry, r),
            "source_absmax": witness(h_source, r),
            "scalar_source_absmax": witness(h_source_scalar, r),
            "fluid_source_absmax": witness(h_source_fluid, r),
            "full_rms": float(np.sqrt(np.mean(h_residual**2))),
            "inside_fixed_radius_rms": float(np.sqrt(np.mean(h_residual[inner_mask]**2))),
            "outside_fixed_radius_rms": float(np.sqrt(np.mean(h_residual[exterior_mask]**2))),
        },
        "momentum": {
            "full": m_full,
            "inside_fixed_radius_r_lt_5": m_inner,
            "outside_fixed_radius_r_ge_5": m_exterior,
            "geometry_at_full_residual_peak": float(m_geometry[m_full["index"]]),
            "source_at_full_residual_peak": float(m_source[m_full["index"]]),
            "source_scalar_at_full_residual_peak": float(m_source_scalar[m_full["index"]]),
            "source_fluid_at_full_residual_peak": float(m_source_fluid[m_full["index"]]),
            "geometry_absmax": witness(m_geometry, r),
            "source_absmax": witness(m_source, r),
            "scalar_source_absmax": witness(m_source_scalar, r),
            "fluid_source_absmax": witness(m_source_fluid, r),
            "full_rms": float(np.sqrt(np.mean(m_residual**2))),
            "inside_fixed_radius_rms": float(np.sqrt(np.mean(m_residual[inner_mask]**2))),
            "outside_fixed_radius_rms": float(np.sqrt(np.mean(m_residual[exterior_mask]**2))),
        },
    }
    profile = {
        "resolution": int(grid.n),
        "t": float(state.t),
        "tau": float(state.tau),
        "r": r.tolist(),
        "H_geometry": h_geometry.tolist(),
        "H_source_scalar": h_source_scalar.tolist(),
        "H_source_fluid": h_source_fluid.tolist(),
        "H_source_total": h_source.tolist(),
        "H_residual": h_residual.tolist(),
        "M_geometry": m_geometry.tolist(),
        "M_source_scalar": m_source_scalar.tolist(),
        "M_source_fluid": m_source_fluid.tolist(),
        "M_source_total": m_source.tolist(),
        "M_residual": m_residual.tolist(),
    }
    return summary, profile


def common_radius_comparison(summaries_by_n_time, profiles_by_n_time) -> list[dict]:
    rows = []
    base_n = 40
    common_times = sorted({float(t) for t in SNAPSHOT_TIMES})
    for t in common_times:
        base = profiles_by_n_time[(base_n, t)]
        base_r = np.asarray(base["r"], dtype=float)
        fields = ("H_geometry", "H_source_total", "H_residual", "M_geometry", "M_source_total", "M_residual")
        row = {"t": t, "reference_resolution": base_n}
        for n in (60, 80):
            other = profiles_by_n_time[(n, t)]
            other_r = np.asarray(other["r"], dtype=float)
            key = f"N{n}_minus_N{base_n}"
            row[key] = {}
            for field in fields:
                b = np.asarray(base[field], dtype=float)
                o = np.interp(base_r, other_r, np.asarray(other[field], dtype=float))
                delta = o - b
                shell = base_r >= 5.0
                row[key][field] = {
                    "rms_on_N40_centers_r1_to_r79": float(np.sqrt(np.mean(delta**2))),
                    "max_abs_on_N40_centers_r1_to_r79": float(np.max(np.abs(delta))),
                    "rms_at_r_ge_5": float(np.sqrt(np.mean(delta[shell]**2))),
                    "max_abs_at_r_ge_5": float(np.max(np.abs(delta[shell]))),
                }
        rows.append(row)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifact-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--source-run-id", type=int, required=True)
    parser.add_argument("--source-code-sha", required=True)
    args = parser.parse_args()
    artifact_dir = args.artifact_dir
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    campaign_path = artifact_dir / "campaign_summary.json"
    if not campaign_path.is_file():
        raise FileNotFoundError(f"missing {campaign_path}")
    campaign = json.loads(campaign_path.read_text())
    config = campaign.get("campaign", {})
    expected = {"resolutions": [40, 60, 80], "r_max": 80.0, "cfl": 0.03, "final_time": 12.0}
    for key, value in expected.items():
        if config.get(key) != value:
            raise ValueError(f"unexpected input campaign {key}: expected {value!r}, got {config.get(key)!r}")
    if campaign.get("admission", {}).get("passed") is not True:
        raise ValueError("source campaign did not pass process-level admission")
    for result in campaign.get("results", []):
        if result.get("status") != "completed" or float(result.get("t_final", -1)) != 12.0:
            raise ValueError(f"source campaign has incomplete resolution result: {result}")

    kernel = V55ProductionKernel()
    _, vacuum, _ = adapter.vendor_modules()
    summaries = []
    profiles_by_n_time = {}
    for n in (40, 60, 80):
        # t=0 is initialized at the same production defaults; no evolution step is taken.
        initial = kernel.initialize(
            resolution=n,
            r_max=80.0,
            amplitude=0.01,
            width=7.0,
            D_amplitude=1.0e-10,
            include_radiation=True,
        )
        summary, profile = audit_state(initial, vacuum)
        summaries.append(summary)
        profiles_by_n_time[(n, 0.0)] = profile

        for t in (2, 4, 6, 8, 10, 12):
            checkpoint = artifact_dir / f"N{n}" / f"checkpoint_t{t:012.6f}.npz"
            if not checkpoint.is_file():
                raise FileNotFoundError(f"missing saved checkpoint {checkpoint}")
            state = restore_checkpoint(kernel, n, 80.0, checkpoint)
            if abs(state.t - t) > 1.0e-10:
                raise ValueError(f"{checkpoint} contains t={state.t}, expected t={t}")
            summary, profile = audit_state(state, vacuum)
            summary["input_checkpoint"] = checkpoint.relative_to(artifact_dir).as_posix()
            summary["input_checkpoint_sha256"] = sha256_file(checkpoint)
            summaries.append(summary)
            profiles_by_n_time[(n, float(t))] = profile

    comparisons = common_radius_comparison(summaries, profiles_by_n_time)
    with (output_dir / "constraint_source_profiles.jsonl").open("w") as handle:
        for n in (40, 60, 80):
            for t in SNAPSHOT_TIMES:
                handle.write(json.dumps(profiles_by_n_time[(n, float(t))], separators=(",", ":")) + "\n")

    output = {
        "schema": "v55_hamiltonian_momentum_source_budget_v1",
        "diagnostic_only": True,
        "source_run_id": args.source_run_id,
        "source_code_commit": args.source_code_sha,
        "source_campaign_summary_sha256": sha256_file(campaign_path),
        "source_campaign_config": config,
        "source_campaign_admission": campaign.get("admission"),
        "input_artifact_root": str(artifact_dir),
        "snapshots_per_resolution": list(SNAPSHOT_TIMES),
        "snapshot_count": len(summaries),
        "method": {
            "hamiltonian_residual": "raw vacuum Hamiltonian constraint - 16*pi*(scalar_rho + fluid_rho)",
            "momentum_residual": "raw vacuum momentum constraint - 8*pi*(scalar_j + fluid_j)",
            "geometry_and_source_terms": "evaluated from the same restored slice using the pinned production operator and stress-energy assembly",
            "common_radius_comparison": "linear interpolation of N=60 and N=80 profiles onto N=40 cell centers r=1..79; fixed-radius exterior comparison uses r>=5",
            "production_steps_run": 0,
        },
        "snapshots": summaries,
        "resolution_comparisons": comparisons,
    }
    (output_dir / "constraint_source_budget_summary.json").write_text(json.dumps(output, indent=2))
    print(json.dumps({
        "schema": output["schema"],
        "source_run_id": args.source_run_id,
        "snapshot_count": output["snapshot_count"],
        "profiles": str(output_dir / "constraint_source_profiles.jsonl"),
        "summary": str(output_dir / "constraint_source_budget_summary.json"),
        "production_steps_run": 0,
        "finite": all(item["finite"] for item in summaries),
    }, indent=2))


if __name__ == "__main__":
    main()
