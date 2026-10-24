#!/usr/bin/env python3
"""Source-pinned production evolution audit of regular-centre A_a/K behaviour.

Research-only observer. Runs the unchanged production kernel at N=60/120/240
to t=0.6, captures the actual CMC/PIRK geometry stages, decomposes the signed
A_a L2 source, and measures even-series origin extrapolants on fixed physical
support. It neither modifies production equations nor rejects a state because
of the measured regularity residual.
"""
from __future__ import annotations

import hashlib
import json
import os
import platform
import subprocess
from pathlib import Path

import numpy as np

SOURCE_COMMIT = "4aec88d0e3ff028e0a65a8ac36b30948625df792"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
SOURCE_BLOBS = {
    "engine/production_kernel.py": "1720615628a98a8d7cfc9fcadcc984fe28c88be1",
    "engine/v55_pirk_adapter.py": "be28ceea63fa1b59efc5cc5c460b9df19ae1cbc3",
    "engine/v55_initial.py": "7db2e1a10cf16fb40e4221b0df504fdf39aae611",
    "engine/v55_matter.py": "0ee3c51fac7dc46c63feba2747a750542bcee3d0",
    "engine/cmc_gauge.py": "7b49b5751abebd4faf0130208032cd2080b71f64",
    "engine/matter_system.py": "0969f4d996e4f8e0fc5b70a7dff147f024907302",
    "engine/matter_rhs.py": "0153c232ee15b5a9d3dbac37ea858e4d2d730011",
    "engine/scalar_system.py": "e2c2dae9a016b4d6454cf74a9eae8b0167d5733d",
    "engine/stress_energy.py": "ac3e7518c62bb72e2218a61489fa486bb787809f",
    "engine/worldtube.py": "aea19090cb0fd2ebe5481a8b7af0d26f119a9ad6",
    "engine/native_geometry/__init__.py": "796d338d7cf2ed87bcb6ae63808dcf08110aeff8",
    "engine/native_geometry/evolution.py": "7b47d700680b12ab3c82cb6a556d40580b103ce9",
    "engine/native_geometry/operators.py": "2748be304b58328e2620f0a7fde49b01c0794ae2",
    "engine/native_geometry/geometry.py": "1bf12c1139572ff4597eaec0570aa43f2563fd91",
    "engine/native_geometry/constraints.py": "124d1f237e54d082b4d76213f9ebe144e02b02ba",
    "engine/handoff.py": "1e0b8a8a126c67dafaa5d370babdd519cbedf974",
    "engine/production_contract.py": "0a3b6f0a4223f93a5e33ab30706dc803f69643ca",
    "engine/invariant_diagnostics.py": "4bbcc6174c978b07c4f2175d1dcfeccb2d0c12b1",
    "engine/cosmology_observables.py": "4d6f39308685bae5e5bad2e17baa9c4043b67b4e",
    "engine/atomic_io.py": "484020bf0e2b8b1a37fe85a0f54b27d4cec5d2dd",
    "engine/radial_maxwell.py": "81269f1d1334c46c4c2f44ba816d792722c0a38e",
    "requirements.txt": "1e2a5bcbfd897f6465081a8794e14966f6e5db09",
    ".gitmodules": "502d49645859dc22c848efee0f654ed2170e18b9",
}

RESOLUTION = int(os.environ.get("REGULARITY_RESOLUTION", "120"))
R_MAX = float(os.environ.get("REGULARITY_R_MAX", "80.0"))
TARGET_TIME = float(os.environ.get("REGULARITY_TARGET_TIME", "0.6"))
CFL = float(os.environ.get("REGULARITY_CFL", "0.03"))
FIT_R_MAX = float(os.environ.get("REGULARITY_FIT_R_MAX", "8.0"))
FIT_DEGREE = int(os.environ.get("REGULARITY_FIT_DEGREE", "3"))
OUT = Path(os.environ.get(
    "REGULARITY_OUT",
    f"results/evolved-origin-aa-regularity-N{RESOLUTION}-t{TARGET_TIME:g}-2026-10-24",
))
GEOMETRY_STAGE_NAMES = (
    "start_lapse_refreshed",
    "explicit_predictor",
    "primary_predictor",
)


def _git(*args: str) -> str:
    return subprocess.check_output(["git", *args], text=True).strip()


def verify_source_pins() -> dict:
    base_vendor = _git("rev-parse", f"{SOURCE_COMMIT}:vendor/bb-palatini-unified-r0")
    head_vendor = _git("rev-parse", "HEAD:vendor/bb-palatini-unified-r0")
    if base_vendor != VENDOR_COMMIT or head_vendor != VENDOR_COMMIT:
        raise RuntimeError(
            f"vendor gitlink mismatch: base={base_vendor}, head={head_vendor}, expected={VENDOR_COMMIT}"
        )
    pin_rows = {}
    for path, expected in SOURCE_BLOBS.items():
        base_blob = _git("rev-parse", f"{SOURCE_COMMIT}:{path}")
        head_blob = _git("rev-parse", f"HEAD:{path}")
        matched = base_blob == expected and head_blob == expected
        pin_rows[path] = {
            "expected_blob": expected,
            "base_blob": base_blob,
            "head_blob": head_blob,
            "matches": matched,
        }
        if not matched:
            raise RuntimeError(f"source pin mismatch for {path}: {pin_rows[path]}")
    changed = _git(
        "diff", "--name-only", SOURCE_COMMIT, "HEAD", "--",
        "engine", "vendor", ".gitmodules", "requirements.txt",
    )
    if changed:
        raise RuntimeError(
            "Production/source/vendor/dependency paths changed from the baseline: " + changed
        )
    return {
        "baseline_commit": SOURCE_COMMIT,
        "audit_head": _git("rev-parse", "HEAD"),
        "vendor_commit": head_vendor,
        "production_vendor_dependency_diff_empty": True,
        "blobs": pin_rows,
    }


def _arr(values) -> np.ndarray:
    return np.asarray(values, dtype=float).copy()


def _max_abs(values) -> float:
    values = _arr(values)
    return float(np.max(np.abs(values))) if values.size else 0.0


def _state_digest(state) -> str:
    arrays = {}
    for name in ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B"):
        arrays["geometry." + name] = getattr(state.geometry, name)
    for name in ("S", "PS", "D", "PD", "phi", "Pi"):
        arrays["scalars." + name] = getattr(state.scalars, name)
    for sector in ("dark_matter", "baryons", "radiation", "charged_dust"):
        item = getattr(state.matter, sector, None)
        if item is None:
            continue
        for name in ("rest", "energy_t", "momentum_r"):
            arrays["matter." + sector + "." + name] = getattr(item, name)
    digest = hashlib.sha256()
    for name in sorted(arrays):
        values = np.ascontiguousarray(np.asarray(arrays[name]))
        digest.update(name.encode("utf-8"))
        digest.update(str(values.dtype).encode("ascii"))
        digest.update(repr(values.shape).encode("ascii"))
        digest.update(values.tobytes())
    if state.electromagnetic is not None:
        for name in ("enclosed_charge_faces", "electric_field_faces"):
            if hasattr(state.electromagnetic, name):
                values = np.ascontiguousarray(np.asarray(getattr(state.electromagnetic, name)))
                digest.update(("electromagnetic." + name).encode("utf-8"))
                digest.update(str(values.dtype).encode("ascii"))
                digest.update(repr(values.shape).encode("ascii"))
                digest.update(values.tobytes())
    return digest.hexdigest()


def _even_origin_fit(r: np.ndarray, values: np.ndarray) -> dict:
    """Fit c0+c2*r^2+... on fixed physical support, excluding cell zero."""
    r = _arr(r)
    y = _arr(values)
    if r.shape != y.shape or r.ndim != 1:
        raise ValueError("even-origin fit requires matching one-dimensional arrays")
    mask = (r > 0.0) & (r <= FIT_R_MAX + 1.0e-12)
    mask[0] = False
    count = int(np.count_nonzero(mask))
    degree = min(FIT_DEGREE, count - 1)
    if degree < 1:
        raise ValueError(f"too few cells for origin fit: Nfit={count}")
    x = (r[mask] / FIT_R_MAX) ** 2
    design = np.column_stack([x**order for order in range(degree + 1)])
    coeff_scaled, *_ = np.linalg.lstsq(design, y[mask], rcond=None)
    fitted = design @ coeff_scaled
    residual = y[mask] - fitted
    physical_coeff = [
        float(coeff_scaled[order] / (FIT_R_MAX ** (2 * order)))
        for order in range(degree + 1)
    ]
    return {
        "origin_intercept": physical_coeff[0],
        "coefficients_c0_c2_etc": physical_coeff,
        "degree": degree,
        "fit_cells": count,
        "fit_r_max": FIT_R_MAX,
        "fit_max_abs_residual": _max_abs(residual),
        "fit_rms_residual": float(np.sqrt(np.mean(residual * residual))),
        "profile_max_abs_on_fit_support": _max_abs(y[mask]),
        "cell0": float(y[0]),
        "cell0_minus_fit_intercept": float(y[0] - physical_coeff[0]),
        "cell0_cell1_cell2_cell3": [float(v) for v in y[:4]],
    }


def _capture_arrays(
    *,
    label: str,
    step_index: int,
    t_start: float,
    grid,
    geometry,
    vacuum,
    terms: dict | None = None,
    l2_only: dict | None = None,
) -> tuple[dict, dict[str, np.ndarray]]:
    r = _arr(grid.centers)
    g = geometry
    field_names = ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda")
    arrays: dict[str, np.ndarray] = {
        "r": r.copy(),
        **{name: _arr(getattr(g, name)) for name in field_names},
    }
    geom_terms = vacuum.geometry_terms(grid, g)
    alpha = _arr(g.alpha)
    aa_lapse = -(
        _arr(geom_terms["DrDr_alpha"]) - _arr(geom_terms["lap_alpha"]) / 3.0
    )
    aa_ricci = alpha * (
        _arr(geom_terms["Rr_r"]) - _arr(geom_terms["R"]) / 3.0
    )
    arrays["Aa_lapse_hessian_tf"] = aa_lapse
    arrays["Aa_tracefree_ricci"] = aa_ricci
    l2_aa = None
    l2_k = None
    l3_aa = None
    l3_k = None
    if terms is not None:
        l2_aa = _arr(terms["primary_l2"]["Aa"])
        l2_k = _arr(terms["primary_l2"]["K"])
        l3_aa = _arr(terms["primary_l3"]["Aa"])
        l3_k = _arr(terms["primary_l3"]["K"])
        arrays["Aa_l3_total"] = l3_aa
        arrays["K_l3_total"] = l3_k
        arrays["a_rhs"] = _arr(terms["explicit"]["a"])
        arrays["b_rhs"] = _arr(terms["explicit"]["b"])
        arrays["X_rhs"] = _arr(terms["explicit"]["X"])
        vacuum_l3 = vacuum.primary_l3_rhs(grid, g)
        arrays["Aa_l3_vacuum"] = _arr(vacuum_l3["Aa"])
        arrays["Aa_l3_matter_increment"] = l3_aa - arrays["Aa_l3_vacuum"]
    elif l2_only is not None:
        l2_aa = _arr(l2_only["Aa"])
        l2_k = _arr(l2_only["K"])
    if l2_aa is not None:
        arrays["Aa_l2_total"] = l2_aa
    if l2_k is not None:
        arrays["K_l2_total"] = l2_k
    if l2_aa is not None:
        arrays["Aa_l2_closure"] = l2_aa - aa_lapse - aa_ricci
    if l2_aa is not None and l3_aa is not None:
        arrays["Aa_rhs_l2_plus_l3"] = l2_aa + l3_aa
        arrays["K_rhs_l2_plus_l3"] = l2_k + l3_k

    def fit_if_present(name: str) -> dict | None:
        return _even_origin_fit(r, arrays[name]) if name in arrays else None

    even_fields = {
        name: _even_origin_fit(r, arrays[name])
        for name in ("a", "b", "X", "alpha", "Aa", "K")
    }
    source_fits = {
        name: fit_if_present(name)
        for name in (
            "Aa_lapse_hessian_tf", "Aa_tracefree_ricci",
            "Aa_l2_total", "Aa_l3_total", "Aa_l3_vacuum",
            "Aa_l3_matter_increment", "Aa_rhs_l2_plus_l3", "K_rhs_l2_plus_l3",
        )
    }
    source_fits = {key: value for key, value in source_fits.items() if value is not None}

    identities = {}
    if terms is not None:
        a, b, X = _arr(g.a), _arr(g.b), _arr(g.X)
        Aa, K = _arr(g.Aa), _arr(g.K)
        identities = {
            "max_abs_a_rhs_plus_2_alpha_a_Aa": _max_abs(
                arrays["a_rhs"] + 2.0 * alpha * a * Aa
            ),
            "max_abs_b_rhs_minus_alpha_b_Aa": _max_abs(
                arrays["b_rhs"] - alpha * b * Aa
            ),
            "max_abs_X_rhs_minus_alpha_X_K_over3": _max_abs(
                arrays["X_rhs"] - alpha * X * K / 3.0
            ),
            "beta_max_abs": _max_abs(g.beta),
        }
    closure = _max_abs(arrays["Aa_l2_closure"]) if "Aa_l2_closure" in arrays else None
    record = {
        "label": label,
        "step_index": int(step_index),
        "accepted_time_at_step_start": float(t_start),
        "regularity_fits_even_fields": even_fields,
        "regularity_fits_Aa_source_channels": source_fits,
        "metric_rhs_identities": identities,
        "Aa_L2_signed_decomposition": {
            "formula": "Aa_L2 = -(DrDr_alpha - lap_alpha/3) + alpha*(Rr_r - R/3)",
            "cell0_lapse_hessian": float(aa_lapse[0]),
            "cell0_tracefree_ricci": float(aa_ricci[0]),
            "cell0_primary_l2": float(l2_aa[0]) if l2_aa is not None else None,
            "cell0_signed_closure": (
                float(arrays["Aa_l2_closure"][0]) if "Aa_l2_closure" in arrays else None
            ),
            "max_abs_closure_full_grid": closure,
        },
        "available_profiles": sorted(k for k in arrays if k != "r"),
    }
    return record, arrays


def main() -> int:
    if RESOLUTION < 16 or R_MAX <= 0.0 or TARGET_TIME <= 0.0 or CFL <= 0.0:
        raise ValueError("invalid resolution / domain / time / CFL")
    if FIT_R_MAX <= 0.0:
        raise ValueError("fit radius must be positive")
    provenance = verify_source_pins()

    from engine.production_kernel import V55ProductionKernel
    from engine.v55_pirk_adapter import vendor_modules

    _, vacuum, _ = vendor_modules()
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=RESOLUTION,
        r_max=R_MAX,
        amplitude=0.01,
        width=7.0,
        D_amplitude=1.0e-10,
        include_radiation=True,
    )
    dt = CFL * R_MAX / RESOLUTION
    steps = int(round(TARGET_TIME / dt))
    if steps < 1 or not np.isclose(steps * dt, TARGET_TIME, rtol=0.0, atol=1.0e-12):
        raise ValueError(
            f"target time {TARGET_TIME} is not an integer number of steps at dt={dt}"
        )

    out_profiles: dict[str, np.ndarray] = {"r": _arr(state.grid.centers)}
    records: list[dict] = []
    step_budgets: list[dict] = []
    active = {
        "enabled": False,
        "inside_geometry_stage": False,
        "step_index": None,
        "t_start": None,
        "geometry_calls": 0,
        "stage_records": [],
        "final_l2_record": None,
    }
    original_geometry_stage = kernel._geometry_stage_terms
    original_primary_l2 = vacuum.primary_l2_rhs

    def capture(label, geometry, terms=None, l2_only=None):
        rec, arrays = _capture_arrays(
            label=label,
            step_index=active["step_index"] if active["step_index"] is not None else -1,
            t_start=active["t_start"] if active["t_start"] is not None else float(state.t),
            grid=state.grid,
            geometry=geometry,
            vacuum=vacuum,
            terms=terms,
            l2_only=l2_only,
        )
        prefix = f"{len(records):04d}_{label}"
        for name, values in arrays.items():
            if name == "r":
                continue
            out_profiles[prefix + "__" + name] = _arr(values)
        rec["profile_prefix"] = prefix
        rec["capture_index"] = len(records)
        records.append(rec)
        return rec, arrays

    def wrapped_geometry_stage(grid, geometry, scalars, matter, electromagnetic):
        if not active["enabled"]:
            return original_geometry_stage(grid, geometry, scalars, matter, electromagnetic)
        active["inside_geometry_stage"] = True
        try:
            terms = original_geometry_stage(grid, geometry, scalars, matter, electromagnetic)
        finally:
            active["inside_geometry_stage"] = False
        call_index = active["geometry_calls"]
        if call_index >= len(GEOMETRY_STAGE_NAMES):
            raise RuntimeError(f"unexpected geometry stage call {call_index}")
        label = GEOMETRY_STAGE_NAMES[call_index]
        active["geometry_calls"] += 1
        rec, arrays = capture("step%03d_%s" % (active["step_index"], label), geometry, terms=terms)
        active["stage_records"].append((label, rec, arrays))
        return terms

    def wrapped_primary_l2(grid, geometry):
        result = original_primary_l2(grid, geometry)
        if active["enabled"] and not active["inside_geometry_stage"]:
            if active["final_l2_record"] is not None:
                raise RuntimeError("more than one out-of-stage primary_l2 call in a production step")
            rec, arrays = capture(
                "step%03d_final_u_old_v" % active["step_index"],
                geometry,
                l2_only=result,
            )
            active["final_l2_record"] = (rec, arrays)
        return result

    kernel._geometry_stage_terms = wrapped_geometry_stage
    vacuum.primary_l2_rhs = wrapped_primary_l2
    try:
        # Accepted-state observer: same grid and operators, no evolution, then
        # verify that observation leaves every saved state array byte-identical.
        for step_index in range(steps + 1):
            before = _state_digest(state)
            terms = original_geometry_stage(
                state.grid, state.geometry, state.scalars, state.matter, state.electromagnetic
            )
            rec, arrays = capture("accepted_%03d" % step_index, state.geometry, terms=terms)
            after = _state_digest(state)
            if before != after:
                raise RuntimeError(f"accepted-state observer mutated state at index {step_index}")
            rec["state_sha256_before_observation"] = before
            rec["state_sha256_after_observation"] = after
            rec["observer_left_state_unchanged"] = True

            if step_index == steps:
                break

            old_state = state
            active.update({
                "enabled": True,
                "inside_geometry_stage": False,
                "step_index": step_index,
                "t_start": float(old_state.t),
                "geometry_calls": 0,
                "stage_records": [],
                "final_l2_record": None,
            })
            state = kernel.step(old_state, dt)
            active["enabled"] = False
            if active["geometry_calls"] != 3 or active["final_l2_record"] is None:
                raise RuntimeError(
                    "production step capture incomplete: "
                    f"geometry calls={active['geometry_calls']}, final L2={active['final_l2_record'] is not None}"
                )
            stages = {name: arrays for name, _, arrays in active["stage_records"]}
            start = stages["start_lapse_refreshed"]
            primary = stages["primary_predictor"]
            final_l2 = active["final_l2_record"][1]
            r = _arr(state.grid.centers)
            interior = r <= 10.0
            actual_Aa_delta = _arr(state.geometry.Aa) - _arr(old_state.geometry.Aa)
            actual_K_delta = _arr(state.geometry.K) - _arr(old_state.geometry.K)
            predicted_Aa_delta = 0.5 * dt * (
                _arr(start["Aa_l2_total"]) + _arr(final_l2["Aa_l2_total"])
                + _arr(start["Aa_l3_total"]) + _arr(primary["Aa_l3_total"])
            )
            predicted_K_delta = 0.5 * dt * (
                _arr(start["K_l2_total"]) + _arr(final_l2["K_l2_total"])
                + _arr(start["K_l3_total"]) + _arr(primary["K_l3_total"])
            )
            closure_Aa = actual_Aa_delta - predicted_Aa_delta
            closure_K = actual_K_delta - predicted_K_delta
            budget = {
                "step_index": step_index,
                "time_start": float(old_state.t),
                "time_end": float(state.t),
                "dt": dt,
                "Aa_cell0_start": float(old_state.geometry.Aa[0]),
                "Aa_cell0_end": float(state.geometry.Aa[0]),
                "Aa_cell0_actual_delta": float(actual_Aa_delta[0]),
                "Aa_cell0_predicted_delta": float(predicted_Aa_delta[0]),
                "Aa_cell0_signed_update_closure": float(closure_Aa[0]),
                "Aa_update_closure_max_abs_r_le_10": _max_abs(closure_Aa[interior]),
                "Aa_update_closure_max_abs_full_grid": _max_abs(closure_Aa),
                "K_cell0_start": float(old_state.geometry.K[0]),
                "K_cell0_end": float(state.geometry.K[0]),
                "K_cell0_actual_delta": float(actual_K_delta[0]),
                "K_cell0_predicted_delta": float(predicted_K_delta[0]),
                "K_cell0_signed_update_closure": float(closure_K[0]),
                "K_update_closure_max_abs_r_le_10": _max_abs(closure_K[interior]),
                "K_update_closure_max_abs_full_grid": _max_abs(closure_K),
                "Aa_regular_origin_fit_start": start.get("Aa", None),
                "Aa_regular_origin_fit_primary_predictor": primary.get("Aa", None),
            }
            # Preserve the signed pointwise residuals in the compact report;
            # full profiles below carry every cell's closure and source channels.
            for name, values in (
                ("Aa_actual_delta", actual_Aa_delta),
                ("Aa_predicted_delta", predicted_Aa_delta),
                ("Aa_update_closure", closure_Aa),
                ("K_actual_delta", actual_K_delta),
                ("K_predicted_delta", predicted_K_delta),
                ("K_update_closure", closure_K),
            ):
                prefix = f"step{step_index:03d}_budget__{name}"
                out_profiles[prefix] = _arr(values)
            step_budgets.append(budget)
            if budget["Aa_update_closure_max_abs_r_le_10"] > 2.0e-12:
                raise RuntimeError(
                    "production A_a update did not close at r<=10 after exact captured-stage reconstruction: "
                    f"step {step_index}, residual={budget['Aa_update_closure_max_abs_r_le_10']:.6e}"
                )
            if budget["K_update_closure_max_abs_r_le_10"] > 2.0e-12:
                raise RuntimeError(
                    "production K update did not close at r<=10 after exact captured-stage reconstruction: "
                    f"step {step_index}, residual={budget['K_update_closure_max_abs_r_le_10']:.6e}"
                )
    finally:
        kernel._geometry_stage_terms = original_geometry_stage
        vacuum.primary_l2_rhs = original_primary_l2

    OUT.mkdir(parents=True, exist_ok=True)
    report = {
        "schema": "source_pinned_evolved_origin_Aa_regularity_v1",
        "status": "completed_bounded_production_origin_regularity_audit",
        "scientific_admission": False,
        "production_or_vendor_files_modified": False,
        "source_pins": provenance,
        "configuration": {
            "resolution": RESOLUTION,
            "r_max": R_MAX,
            "dr": float(R_MAX / RESOLUTION),
            "first_cell_center": float(R_MAX / RESOLUTION / 2.0),
            "cfl": CFL,
            "dt": dt,
            "target_time": TARGET_TIME,
            "accepted_steps": steps,
            "fit": {
                "parity_assumption": "even Taylor series in r: c0+c2*r^2+c4*r^4+...",
                "physical_support": f"0 < r <= {FIT_R_MAX}",
                "cell_zero_excluded_from_fit": True,
                "degree": FIT_DEGREE,
                "purpose": "extrapolation diagnostic only; no regularity acceptance threshold",
            },
        },
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "platform": platform.platform(),
        },
        "initial_state_sha256": None,
        "final_state_sha256": _state_digest(state),
        "accepted_state_observer_mutation_check_passed": all(
            row.get("observer_left_state_unchanged", False)
            for row in records if row["label"].startswith("accepted_")
        ),
        "records": records,
        "accepted_step_Aa_K_update_budgets": step_budgets,
        "interpretation_limits": [
            "The fit window and polynomial degree are explicit diagnostics, not a physical regularity bound.",
            "A nonzero extrapolated Aa intercept flags incompatibility with the expected even, trace-free centre series for this discrete profile; it does not identify the continuum equation responsible.",
            "The A_a L2 split is the pinned production source decomposition; its algebraic closure alone is not a continuum regularity proof.",
            "Stage histories distinguish the start-lapse-refreshed slice, explicit predictor, primary predictor, final implicit L2 operand, and accepted state. This does not supply internal sub-stages outside the documented production PIRK sequence.",
            "The accepted A_a/K update closure is checked only on r<=10 to avoid conflating the interior source budget with the finite outer-boundary reconstruction.",
            "No production physics, equation, boundary, coefficient, threshold, or scientific-admission flag is changed by this diagnostic.",
        ],
    }
    # Record the first digest separately after all source-preserving runs.
    init_record = next(row for row in records if row["label"] == "accepted_000")
    report["initial_state_sha256"] = init_record["state_sha256_before_observation"]
    report_path = OUT / "report.json"
    OUT.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    npz_path = OUT / "stage_profiles.npz"
    np.savez_compressed(npz_path, **out_profiles)
    checkpoint_path = OUT / "final_state.npz"
    kernel.checkpoint(state, checkpoint_path)
    sums = {}
    for path in (report_path, npz_path, checkpoint_path):
        sums[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
    (OUT / "SHA256SUMS.json").write_text(json.dumps(sums, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    final_record = next(row for row in records if row["label"] == f"accepted_{steps:03d}")
    print(json.dumps({
        "status": report["status"],
        "resolution": RESOLUTION,
        "target_time": TARGET_TIME,
        "dt": dt,
        "accepted_steps": steps,
        "initial_Aa_cell0": init_record["regularity_fits_even_fields"]["Aa"]["cell0"],
        "initial_Aa_even_fit_origin": init_record["regularity_fits_even_fields"]["Aa"]["origin_intercept"],
        "final_Aa_cell0": final_record["regularity_fits_even_fields"]["Aa"]["cell0"],
        "final_Aa_even_fit_origin": final_record["regularity_fits_even_fields"]["Aa"]["origin_intercept"],
        "final_Aa_cell0_minus_fit": final_record["regularity_fits_even_fields"]["Aa"]["cell0_minus_fit_intercept"],
        "final_K_cell0_minus_fit": final_record["regularity_fits_even_fields"]["K"]["cell0_minus_fit_intercept"],
        "max_Aa_update_closure_r_le_10": max(x["Aa_update_closure_max_abs_r_le_10"] for x in step_budgets),
        "max_K_update_closure_r_le_10": max(x["K_update_closure_max_abs_r_le_10"] for x in step_budgets),
        "artifact_checksums": sums,
        "scientific_admission": False,
    }, indent=2, allow_nan=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
