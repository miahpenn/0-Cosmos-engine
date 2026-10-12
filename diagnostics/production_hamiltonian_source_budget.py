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
from engine.production_kernel import ProductionState, V55ProductionKernel
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


H_GROUPS = (
    # Split the former spatial-metric group into its evolved variables. This
    # identifies which variable leads the observed H drift without an evolution
    # step; it remains a Shapley attribution on saved endpoint states.
    ("a", "geometry", ("a",)),
    ("b", "geometry", ("b",)),
    ("X", "geometry", ("X",)),
    ("lapse_shift_gauge", "geometry", ("alpha", "beta", "B")),
    ("A_a", "geometry", ("Aa",)),
    ("K", "geometry", ("K",)),
    ("Lambda", "geometry", ("Lambda",)),
    ("matter", "matter", ("scalars", "matter")),
)
H_GROUP_COUNT = len(H_GROUPS)
H_FULL_MASK = (1 << H_GROUP_COUNT) - 1


def clone_production_state(state: ProductionState) -> ProductionState:
    """Copy an evolved state while sharing only the immutable spatial grid."""
    return ProductionState(
        grid=state.grid,
        geometry=state.geometry.copy(),
        scalars=state.scalars.copy(),
        matter=state.matter.copy(),
        t=float(state.t),
        tau=float(state.tau),
        e_folds=float(state.e_folds),
    )


def hamiltonian_for_state(state: ProductionState, vacuum) -> np.ndarray:
    total = assemble_total_stress_energy(
        state.grid, state.geometry, state.scalars, state.matter
    )
    raw = vacuum.constraints(state.grid, state.geometry)
    return (
        np.asarray(raw["hamiltonian"], dtype=float)
        - 16.0 * math.pi * np.asarray(total.rho, dtype=float)
    )


def hybrid_production_state(
    before: ProductionState, after: ProductionState, mask: int
) -> ProductionState:
    """Create a fresh pre/post hybrid state for the selected groups."""
    hybrid = clone_production_state(before)
    for index, (_, owner, fields) in enumerate(H_GROUPS):
        if not (mask & (1 << index)):
            continue
        if owner == "geometry":
            for field in fields:
                setattr(
                    hybrid.geometry,
                    field,
                    np.array(getattr(after.geometry, field), dtype=float, copy=True),
                )
        else:
            hybrid.scalars = after.scalars.copy()
            hybrid.matter = after.matter.copy()
    return hybrid


def decompose_saved_interval(
    before: ProductionState, after: ProductionState, vacuum
) -> tuple[dict, dict]:
    """Shapley-attribute H(after)-H(before) across six saved-state groups."""
    if before.grid.n != after.grid.n or not np.array_equal(
        before.grid.centers, after.grid.centers
    ):
        raise ValueError("saved interval endpoints must use identical spatial grids")

    values = {}
    for mask in range(H_FULL_MASK + 1):
        hybrid = hybrid_production_state(before, after, mask)
        values[mask] = hamiltonian_for_state(hybrid, vacuum)

    denominator = math.factorial(H_GROUP_COUNT)
    component_profiles = {}
    for index, (label, _, _) in enumerate(H_GROUPS):
        bit = 1 << index
        contribution = np.zeros_like(values[0], dtype=float)
        for mask in range(H_FULL_MASK + 1):
            if mask & bit:
                continue
            size = mask.bit_count()
            weight = (
                math.factorial(size)
                * math.factorial(H_GROUP_COUNT - size - 1)
                / denominator
            )
            contribution += weight * (values[mask | bit] - values[mask])
        component_profiles[label] = contribution

    h_before = values[0]
    h_after = values[H_FULL_MASK]
    delta_h = h_after - h_before
    component_sum = sum(
        component_profiles.values(), np.zeros_like(delta_h, dtype=float)
    )
    closure = delta_h - component_sum
    radius = np.asarray(before.grid.centers, dtype=float)
    inner = radius < 5.0
    exterior = radius >= 5.0
    finite = all(
        np.all(np.isfinite(array))
        for array in (
            *values.values(), delta_h, closure, *component_profiles.values()
        )
    )

    component_summary = {}
    for label, profile in component_profiles.items():
        inner_index = int(np.flatnonzero(inner)[np.argmax(np.abs(profile[inner]))])
        exterior_index = int(np.flatnonzero(exterior)[np.argmax(np.abs(profile[exterior]))])
        component_summary[label] = {
            "cell0_signed_delta_H": float(profile[0]),
            "max_abs_full": float(np.max(np.abs(profile))),
            "max_abs_inside_r_lt_5": float(abs(profile[inner_index])),
            "radius_of_inner_max": float(radius[inner_index]),
            "max_abs_exterior_r_ge_5": float(abs(profile[exterior_index])),
            "radius_of_exterior_max": float(radius[exterior_index]),
            "full_rms": float(np.sqrt(np.mean(profile**2))),
            "inside_r_lt_5_rms": float(np.sqrt(np.mean(profile[inner]**2))),
            "exterior_r_ge_5_rms": float(np.sqrt(np.mean(profile[exterior]**2))),
        }

    cells = []
    for i in range(min(5, len(radius))):
        cells.append({
            "cell": i,
            "r": float(radius[i]),
            "H_before": float(h_before[i]),
            "H_after": float(h_after[i]),
            "delta_H": float(delta_h[i]),
            "component_delta_H": {
                label: float(profile[i])
                for label, profile in component_profiles.items()
            },
            "sum_component_delta_H": float(component_sum[i]),
            "closure_error": float(closure[i]),
        })

    summary = {
        "finite": bool(finite),
        "max_abs_H_before_full": float(np.max(np.abs(h_before))),
        "max_abs_H_after_full": float(np.max(np.abs(h_after))),
        "max_abs_delta_H_full": float(np.max(np.abs(delta_h))),
        "max_abs_delta_H_inside_r_lt_5": float(np.max(np.abs(delta_h[inner]))),
        "max_abs_delta_H_exterior_r_ge_5": float(np.max(np.abs(delta_h[exterior]))),
        "closure_max_abs_full": float(np.max(np.abs(closure))),
        "components": component_summary,
        "cells_0_to_4": cells,
    }
    profile = {
        "radius": radius.tolist(),
        "H_before": h_before.tolist(),
        "H_after": h_after.tolist(),
        "delta_H": delta_h.tolist(),
        "component_delta_H": {
            label: values.tolist() for label, values in component_profiles.items()
        },
        "closure_error": closure.tolist(),
    }
    return summary, profile


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
    states_by_n_time = {}
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
        states_by_n_time[(n, 0.0)] = initial

        for t in (2, 4, 6, 8, 10, 12):
            checkpoint = artifact_dir / f"N{n}" / f"checkpoint_t{t:012.6f}.npz"
            if not checkpoint.is_file():
                raise FileNotFoundError(f"missing saved checkpoint {checkpoint}")
            state = restore_checkpoint(kernel, n, 80.0, checkpoint)
            # Checkpoints are stored on the nearest accepted step surrounding
            # the scheduled target; keep the actual timestamp and permit at most
            # half a nominal production step of offset.
            checkpoint_time_tolerance = 0.5 * 0.03 * float(state.grid.dr) + 1.0e-10
            if abs(state.t - t) > checkpoint_time_tolerance:
                raise ValueError(
                    f"{checkpoint} contains t={state.t}, outside half-step "
                    f"tolerance {checkpoint_time_tolerance} around scheduled t={t}"
                )
            summary, profile = audit_state(state, vacuum)
            summary["scheduled_checkpoint_time"] = float(t)
            summary["checkpoint_time_offset"] = float(state.t - t)
            summary["input_checkpoint"] = checkpoint.relative_to(artifact_dir).as_posix()
            summary["input_checkpoint_sha256"] = sha256_file(checkpoint)
            summaries.append(summary)
            profile["scheduled_checkpoint_time"] = float(t)
            profile["checkpoint_time_offset"] = float(state.t - t)
            profiles_by_n_time[(n, float(t))] = profile
            states_by_n_time[(n, float(t))] = state

    comparisons = common_radius_comparison(summaries, profiles_by_n_time)
    increment_attributions = []
    increment_profiles = []
    for n in (40, 60, 80):
        for t_before, t_after in zip(SNAPSHOT_TIMES[:-1], SNAPSHOT_TIMES[1:]):
            before = states_by_n_time[(n, float(t_before))]
            after = states_by_n_time[(n, float(t_after))]
            interval_summary, interval_profile = decompose_saved_interval(
                before, after, vacuum
            )
            interval_summary.update({
                "resolution": int(n),
                "scheduled_t_before": float(t_before),
                "scheduled_t_after": float(t_after),
                "actual_t_before": float(before.t),
                "actual_t_after": float(after.t),
                "actual_dt": float(after.t - before.t),
                "tau_before": float(before.tau),
                "tau_after": float(after.tau),
            })
            increment_attributions.append(interval_summary)
            interval_profile.update({
                "resolution": int(n),
                "scheduled_t_before": float(t_before),
                "scheduled_t_after": float(t_after),
                "actual_t_before": float(before.t),
                "actual_t_after": float(after.t),
                "actual_dt": float(after.t - before.t),
            })
            increment_profiles.append(interval_profile)

    (output_dir / "checkpoint_hamiltonian_increment_attribution.json").write_text(
        json.dumps({
            "schema": "v55_checkpoint_H_increment_shapley_by_variable_v2",
            "diagnostic_only": True,
            "source_run_id": args.source_run_id,
            "source_code_commit": args.source_code_sha,
            "source_campaign_summary_sha256": sha256_file(campaign_path),
            "interval_count": len(increment_attributions),
            "method": "Shapley attribution of H(after)-H(before) across separate a, b, X, lapse/shift gauge (alpha,beta,B), A_a, K, Lambda, and matter (all scalars and fluids); all 256 pre/post hybrids evaluated on the same saved grid; no evolution steps",
            "intervals": increment_attributions,
        }, indent=2)
    )
    with (output_dir / "checkpoint_hamiltonian_increment_profiles.jsonl").open("w") as handle:
        for profile in increment_profiles:
            handle.write(json.dumps(profile, separators=(",", ":")) + "\n")

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
        "checkpoint_increment_interval_count": len(increment_attributions),
        "method": {
            "hamiltonian_residual": "raw vacuum Hamiltonian constraint - 16*pi*(scalar_rho + fluid_rho)",
            "momentum_residual": "raw vacuum momentum constraint - 8*pi*(scalar_j + fluid_j)",
            "geometry_and_source_terms": "evaluated from the same restored slice using the pinned production operator and stress-energy assembly",
            "common_radius_comparison": "linear interpolation of N=60 and N=80 profiles onto N=40 cell centers r=1..79 grouped by scheduled checkpoint; actual checkpoint-time offsets are retained and at most half a nominal step; fixed-radius exterior comparison uses r>=5",
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
        "checkpoint_increment_attribution": str(output_dir / "checkpoint_hamiltonian_increment_attribution.json"),
        "checkpoint_increment_interval_count": len(increment_attributions),
        "production_steps_run": 0,
        "finite": all(item["finite"] for item in summaries),
    }, indent=2))


if __name__ == "__main__":
    main()
