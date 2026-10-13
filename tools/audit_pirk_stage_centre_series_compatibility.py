#!/usr/bin/env python3
"""Capture the real V5.5 PIRK stage states and centre-series witnesses.

Diagnostic-only: wraps the production split-term function without replacing its
output, performs one normal production step at each listed resolution, records
the actual intermediate states, and checks the implemented final weighted
increments against the same stage terms. It does not project or modify physics.
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

# Running a script by path puts tools/ on sys.path, not the repository root.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

SOURCE_COMMIT = "6d8b3f2790644591326df57c55908863f994cb0c"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
RESOLUTIONS = (40, 60, 80, 120, 160)
R_MAX = 80.0
CFL = 0.03
FIXED_RADII = (2.0, 3.0, 4.0, 5.0, 8.0)
CENTRE_CELLS = tuple(range(6))
OUT = Path("results/pirk-stage-centre-series-compatibility")
GEOMETRY_FIELDS = ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B")
SCALAR_FIELDS = ("S", "PS", "D", "PD", "phi", "Pi")
MATTER_SECTORS = ("dark_matter", "baryons", "radiation")
MATTER_FIELDS = ("rest", "energy_t", "momentum_r")
QUOTIENT_FIELDS = (
    "Aa_over_r2",
    "Lambda_over_r",
    "a_over_b_minus_1_over_r2",
    "Dr_a_over_r",
    "Dr_b_over_r",
    "Dr_X_over_r",
    "Dr_alpha_over_r",
    "Dr_K_over_r",
    "Dr_S_over_r",
    "Dr_PS_over_r",
    "Dr_D_over_r",
    "Dr_PD_over_r",
    "Dr_phi_over_r",
    "Dr_Pi_over_r",
    "scalar_j_over_r",
    "fluid_j_over_r",
    "total_j_over_r",
    "Aa_rhs_L2_over_r2",
    "Aa_rhs_L3_over_r2",
    "Aa_rhs_total_over_r2",
    "Dr_K_rhs_L2_over_r",
    "Dr_K_rhs_L3_over_r",
    "Dr_K_rhs_total_over_r",
    "Lambda_rhs_L2_over_r",
    "Lambda_rhs_L3_over_r",
    "Lambda_rhs_total_over_r",
    "Dr_a_rhs_over_r",
    "Dr_b_rhs_over_r",
    "Dr_X_rhs_over_r",
)


def _git_output(*args: str) -> str:
    return subprocess.run(
        ["git", *args], check=True, capture_output=True, text=True
    ).stdout.strip()


def _array(value) -> np.ndarray:
    return np.asarray(value, dtype=float).copy()


def _max_abs(values) -> float:
    values = np.asarray(values, dtype=float)
    if values.size == 0:
        return 0.0
    return float(np.max(np.abs(values)))


def _finite(values) -> bool:
    return bool(np.all(np.isfinite(np.asarray(values, dtype=float))))


def _state_digest(state) -> str:
    """Hash every input-state numerical array to detect unintended mutation."""
    h = hashlib.sha256()
    h.update(float(state.t).hex().encode("ascii"))
    h.update(float(state.tau).hex().encode("ascii"))
    h.update(float(state.e_folds).hex().encode("ascii"))
    named = {}
    for name in GEOMETRY_FIELDS:
        named[f"geometry.{name}"] = getattr(state.geometry, name)
    for name in SCALAR_FIELDS:
        named[f"scalars.{name}"] = getattr(state.scalars, name)
    for sector in MATTER_SECTORS:
        bundle = getattr(state.matter, sector)
        for name in MATTER_FIELDS:
            named[f"matter.{sector}.{name}"] = getattr(bundle, name)
    for name in sorted(named):
        values = np.ascontiguousarray(np.asarray(named[name]))
        h.update(name.encode("utf-8"))
        h.update(str(values.dtype).encode("ascii"))
        h.update(repr(values.shape).encode("ascii"))
        h.update(values.tobytes())
    return h.hexdigest()


def _d1(grid, values, parity=+1) -> np.ndarray:
    return _array(grid.cell_derivative_fourth(np.asarray(values, dtype=float), parity=parity))


def _terms_to_rates(terms: dict, grid) -> dict[str, np.ndarray]:
    """Return signed L2/L3 and full rate sums from the exact split result."""
    rates: dict[str, np.ndarray] = {}
    if "primary_l2" in terms:
        for field in ("Aa", "K"):
            rates[f"{field}_rhs_L2"] = _array(terms["primary_l2"][field])
    if "primary_l3" in terms:
        for field in ("Aa", "K"):
            rates[f"{field}_rhs_L3"] = _array(terms["primary_l3"][field])
    for field in ("Aa", "K"):
        l2 = rates.get(f"{field}_rhs_L2")
        l3 = rates.get(f"{field}_rhs_L3")
        if l2 is not None and l3 is not None:
            rates[f"{field}_rhs_total"] = l2 + l3

    if "lambda_l2" in terms:
        rates["Lambda_rhs_L2"] = _array(terms["lambda_l2"])
    if "lambda_l3" in terms:
        rates["Lambda_rhs_L3"] = _array(terms["lambda_l3"])
    if "Lambda_rhs_L2" in rates and "Lambda_rhs_L3" in rates:
        rates["Lambda_rhs_total"] = rates["Lambda_rhs_L2"] + rates["Lambda_rhs_L3"]

    if "explicit" in terms:
        for field in ("a", "b", "X"):
            if field in terms["explicit"]:
                rates[f"{field}_rhs_explicit"] = _array(terms["explicit"][field])
    return rates


def _quotients(grid, geometry, scalars, sources, rates) -> dict[str, np.ndarray]:
    r = _array(grid.centers)
    q: dict[str, np.ndarray] = {
        "Aa_over_r2": _array(geometry.Aa) / r**2,
        "Lambda_over_r": _array(geometry.Lambda) / r,
        "a_over_b_minus_1_over_r2": (_array(geometry.a) / _array(geometry.b) - 1.0) / r**2,
    }
    for name in ("a", "b", "X", "alpha", "K"):
        q[f"Dr_{name}_over_r"] = _d1(grid, getattr(geometry, name), +1) / r
    for name in SCALAR_FIELDS:
        q[f"Dr_{name}_over_r"] = _d1(grid, getattr(scalars, name), +1) / r

    if sources is not None:
        q["scalar_j_over_r"] = _array(sources.scalar_j) / r
        q["fluid_j_over_r"] = _array(sources.fluid_j) / r
        q["total_j_over_r"] = _array(sources.j) / r

    for field in ("Aa", "K", "Lambda"):
        l2_key = f"{field}_rhs_L2" if field != "Lambda" else "Lambda_rhs_L2"
        l3_key = f"{field}_rhs_L3" if field != "Lambda" else "Lambda_rhs_L3"
        total_key = f"{field}_rhs_total" if field != "Lambda" else "Lambda_rhs_total"
        if l2_key in rates:
            rate = rates[l2_key]
            if field == "Aa":
                q["Aa_rhs_L2_over_r2"] = rate / r**2
            elif field == "K":
                q["Dr_K_rhs_L2_over_r"] = _d1(grid, rate, +1) / r
            else:
                q["Lambda_rhs_L2_over_r"] = rate / r
        if l3_key in rates:
            rate = rates[l3_key]
            if field == "Aa":
                q["Aa_rhs_L3_over_r2"] = rate / r**2
            elif field == "K":
                q["Dr_K_rhs_L3_over_r"] = _d1(grid, rate, +1) / r
            else:
                q["Lambda_rhs_L3_over_r"] = rate / r
        if total_key in rates:
            rate = rates[total_key]
            if field == "Aa":
                q["Aa_rhs_total_over_r2"] = rate / r**2
            elif field == "K":
                q["Dr_K_rhs_total_over_r"] = _d1(grid, rate, +1) / r
            else:
                q["Lambda_rhs_total_over_r"] = rate / r

    for field in ("a", "b", "X"):
        key = f"{field}_rhs_explicit"
        if key in rates:
            q[f"Dr_{field}_rhs_over_r"] = _d1(grid, rates[key], +1) / r
    return q


def _series_summary(values, radii) -> dict:
    values = np.asarray(values, dtype=float)
    finite_mask = np.isfinite(values)
    if not bool(np.all(finite_mask)):
        return {"all_finite": False}
    peak = int(np.argmax(np.abs(values)))
    centres = [
        {"i": int(i), "r": float(radii[i]), "value": float(values[i])}
        for i in CENTRE_CELLS if i < values.size
    ]
    fixed = []
    for target in FIXED_RADII:
        i = int(np.argmin(np.abs(radii - target)))
        fixed.append(
            {"target_r": float(target), "i": i, "r": float(radii[i]), "value": float(values[i])}
        )
    return {
        "all_finite": True,
        "max_abs": _max_abs(values),
        "peak_i": peak,
        "peak_r": float(radii[peak]),
        "peak_signed": float(values[peak]),
        "rms": float(np.sqrt(np.mean(values**2))),
        "centre_cells": centres,
        "fixed_radius_samples": fixed,
    }


def _snapshot(
    *,
    name: str,
    grid,
    geometry,
    scalars=None,
    matter=None,
    terms=None,
    t: float,
    tau: float = 0.0,
    e_folds: float = 0.0,
) -> tuple[dict, dict[str, np.ndarray]]:
    r = _array(grid.centers)
    n = int(r.size)
    arrays: dict[str, np.ndarray] = {"r": r}
    for field in GEOMETRY_FIELDS:
        arrays[f"g_{field}"] = _array(getattr(geometry, field))
    if scalars is not None:
        for field in SCALAR_FIELDS:
            arrays[f"s_{field}"] = _array(getattr(scalars, field))
    if matter is not None:
        for sector in MATTER_SECTORS:
            bundle = getattr(matter, sector)
            for field in MATTER_FIELDS:
                arrays[f"m_{sector}_{field}"] = _array(getattr(bundle, field))

    rates = _terms_to_rates(terms, grid) if terms is not None else {}
    arrays.update({f"rhs_{name}": values for name, values in rates.items()})
    sources = None
    if scalars is not None and matter is not None:
        from engine.stress_energy import assemble_total_stress_energy
        sources = assemble_total_stress_energy(grid, geometry, scalars, matter)
        for field in ("rho", "pr", "pt", "j", "scalar_j", "fluid_j"):
            arrays[f"source_{field}"] = _array(getattr(sources, field))

    quotients = _quotients(grid, geometry, scalars, sources, rates) if scalars is not None else {}
    arrays.update({f"q_{key}": value for key, value in quotients.items()})
    quotient_summary = {
        key: _series_summary(value, r) for key, value in quotients.items()
    }

    closure = {}
    if "Aa_rhs_L2" in rates and "Aa_rhs_L3" in rates:
        closure["Aa_L2_plus_L3"] = _max_abs(
            rates["Aa_rhs_total"] - rates["Aa_rhs_L2"] - rates["Aa_rhs_L3"]
        )
    if "K_rhs_L2" in rates and "K_rhs_L3" in rates:
        closure["K_L2_plus_L3"] = _max_abs(
            rates["K_rhs_total"] - rates["K_rhs_L2"] - rates["K_rhs_L3"]
        )
    if "Lambda_rhs_L2" in rates and "Lambda_rhs_L3" in rates:
        closure["Lambda_L2_plus_L3"] = _max_abs(
            rates["Lambda_rhs_total"] - rates["Lambda_rhs_L2"] - rates["Lambda_rhs_L3"]
        )

    all_arrays_finite = all(_finite(values) for values in arrays.values())
    record = {
        "stage": name,
        "N": n,
        "r_max": float(n * float(grid.dr)),
        "dr": float(grid.dr),
        "t": float(t),
        "tau": float(tau),
        "e_folds": float(e_folds),
        "all_arrays_finite": bool(all_arrays_finite),
        "series_quotients": quotient_summary,
        "signed_rate_reconstruction_closure_max_abs": closure,
    }
    return record, arrays


def main() -> None:
    # Source provenance is checked by the workflow; these runtime checks ensure
    # that this report carries the same source/vendor identity in its metadata.
    head = _git_output("rev-parse", "HEAD")
    source_blob = _git_output("rev-parse", f"{SOURCE_COMMIT}:engine/production_kernel.py")
    head_blob = _git_output("rev-parse", "HEAD:engine/production_kernel.py")
    vendor = _git_output("rev-parse", "HEAD:vendor/bb-palatini-unified-r0")
    if source_blob != head_blob:
        raise RuntimeError(
            f"production_kernel.py is not byte-identical to frozen source {SOURCE_COMMIT}: "
            f"source blob={source_blob}, checked-out blob={head_blob}"
        )
    if vendor != VENDOR_COMMIT:
        raise RuntimeError(f"pinned vendor mismatch: expected {VENDOR_COMMIT}, found {vendor}")

    from engine import v55_pirk_adapter as adapter
    from engine.production_kernel import ProductionState, V55ProductionKernel

    grid_ops, vacuum, _ = adapter.vendor_modules()
    all_reports = []
    all_profiles: dict[str, np.ndarray] = {}
    implementation_closures = []

    for resolution in RESOLUTIONS:
        kernel = V55ProductionKernel()
        state = kernel.initialize(
            resolution=resolution,
            r_max=R_MAX,
            amplitude=0.01,
            width=7.0,
            D_amplitude=1.0e-10,
            include_radiation=True,
        )
        dt = CFL * float(state.grid.dr)
        input_digest = _state_digest(state)

        stage_terms_records = []
        final_l2_captures = []
        final_lambda_captures = []
        inside_geometry_stage_terms = False
        original_terms = adapter.geometry_stage_terms
        original_primary_l2 = vacuum.primary_l2_rhs
        original_lambda_l2 = vacuum.lambda_l2_rhs

        stage_names = (
            "rhs0_initial_slice",
            "rhs_explicit_predictor",
            "rhs_primary_predictor",
        )

        def wrapped_geometry_stage_terms(grid, geometry, fields, matter, lambda_m=2.0):
            nonlocal inside_geometry_stage_terms
            if len(stage_terms_records) >= len(stage_names):
                raise RuntimeError("unexpected additional production geometry-stage evaluation")
            inside_geometry_stage_terms = True
            try:
                terms = original_terms(grid, geometry, fields, matter, lambda_m=lambda_m)
            finally:
                inside_geometry_stage_terms = False
            label = stage_names[len(stage_terms_records)]
            record, profiles = _snapshot(
                name=label,
                grid=grid,
                geometry=geometry,
                scalars=fields,
                matter=matter,
                terms=terms,
                t=float(state.t),
                tau=float(state.tau),
                e_folds=float(state.e_folds),
            )
            stage_terms_records.append({"record": record, "profiles": profiles, "terms": terms})
            all_profiles.update(
                {f"N{resolution}_{label}_{key}": value for key, value in profiles.items()}
            )
            return terms

        def wrapped_primary_l2(grid, geometry):
            result = original_primary_l2(grid, geometry)
            if not inside_geometry_stage_terms:
                record, profiles = _snapshot(
                    name="final_l2_explicit_state",
                    grid=grid,
                    geometry=geometry,
                    t=float(state.t),
                    tau=float(state.tau),
                    e_folds=float(state.e_folds),
                )
                profiles["rhs_Aa_L2_final"] = _array(result["Aa"])
                profiles["rhs_K_L2_final"] = _array(result["K"])
                # The generic snapshot is geometry-only for this final PIRK
                # input, which is exactly the L2 operation's argument.
                record["l2_arrays_finite"] = bool(
                    _finite(result["Aa"]) and _finite(result["K"])
                )
                record["series_quotients"]["Aa_rhs_L2_over_r2"] = _series_summary(
                    _array(result["Aa"]) / _array(grid.centers)**2,
                    _array(grid.centers),
                )
                record["series_quotients"]["Dr_K_rhs_L2_over_r"] = _series_summary(
                    _d1(grid, result["K"], +1) / _array(grid.centers),
                    _array(grid.centers),
                )
                final_l2_captures.append({"record": record, "profiles": profiles, "rhs": result})
            return result

        def wrapped_lambda_l2(grid, geometry, lambda_m=2.0):
            result = original_lambda_l2(grid, geometry, lambda_m=lambda_m)
            if not inside_geometry_stage_terms:
                record, profiles = _snapshot(
                    name="final_lambda_l2_primary_state",
                    grid=grid,
                    geometry=geometry,
                    t=float(state.t),
                    tau=float(state.tau),
                    e_folds=float(state.e_folds),
                )
                profiles["rhs_Lambda_L2_final"] = _array(result)
                record["l2_arrays_finite"] = bool(_finite(result))
                record["series_quotients"]["Lambda_rhs_L2_over_r"] = _series_summary(
                    _array(result) / _array(grid.centers),
                    _array(grid.centers),
                )
                final_lambda_captures.append({"record": record, "profiles": profiles, "rhs": result})
            return result

        adapter.geometry_stage_terms = wrapped_geometry_stage_terms
        vacuum.primary_l2_rhs = wrapped_primary_l2
        vacuum.lambda_l2_rhs = wrapped_lambda_l2
        try:
            accepted = kernel.step(state, dt)
        finally:
            adapter.geometry_stage_terms = original_terms
            vacuum.primary_l2_rhs = original_primary_l2
            vacuum.lambda_l2_rhs = original_lambda_l2

        output_digest = _state_digest(state)
        if input_digest != output_digest:
            raise AssertionError(
                f"N={resolution}: step mutated the input state in place "
                f"(before={input_digest}, after={output_digest})"
            )
        if len(stage_terms_records) != 3:
            raise AssertionError(
                f"N={resolution}: expected three actual geometry split stages, "
                f"captured {len(stage_terms_records)}"
            )
        if len(final_l2_captures) != 1 or len(final_lambda_captures) != 1:
            raise AssertionError(
                f"N={resolution}: expected one direct final Aa/K L2 and Lambda L2 call; "
                f"captured {len(final_l2_captures)} and {len(final_lambda_captures)}"
            )
        if not (np.isfinite(accepted.t) and accepted.t > state.t):
            raise AssertionError(f"N={resolution}: normal production step did not advance time")
        # The geometry contract raises on failure and returns None on success.
        accepted.geometry.assert_finite_positive()

        accepted_record, accepted_profiles = _snapshot(
            name="accepted_final_state",
            grid=accepted.grid,
            geometry=accepted.geometry,
            scalars=accepted.scalars,
            matter=accepted.matter,
            t=float(accepted.t),
            tau=float(accepted.tau),
            e_folds=float(accepted.e_folds),
        )
        all_profiles.update(
            {f"N{resolution}_accepted_final_state_{key}": value for key, value in accepted_profiles.items()}
        )

        # Recompute the accepted state's next-step RHS on a copy, without
        # applying a second evolution step or changing the accepted state.
        final_geometry = accepted.geometry.copy()
        final_alpha, _ = kernel._solve_lapse(
            accepted.grid, final_geometry, accepted.scalars, accepted.matter
        )
        final_geometry.alpha = _array(final_alpha)
        final_geometry.beta.fill(0.0)
        final_geometry.B.fill(0.0)
        final_stage = ProductionState(
            grid=accepted.grid,
            geometry=final_geometry,
            scalars=accepted.scalars,
            matter=accepted.matter,
            t=accepted.t,
            tau=accepted.tau,
            e_folds=accepted.e_folds,
        )
        kernel._rhs(final_stage)  # Force the same scalar/fluid RHS path as a real next step.
        post_terms = adapter.geometry_stage_terms(
            accepted.grid, final_geometry, accepted.scalars, accepted.matter
        )
        post_record, post_profiles = _snapshot(
            name="rhs_poststep_accepted_state",
            grid=accepted.grid,
            geometry=final_geometry,
            scalars=accepted.scalars,
            matter=accepted.matter,
            terms=post_terms,
            t=float(accepted.t),
            tau=float(accepted.tau),
            e_folds=float(accepted.e_folds),
        )
        all_profiles.update(
            {f"N{resolution}_rhs_poststep_accepted_state_{key}": value for key, value in post_profiles.items()}
        )

        # Verify exact code-path weighted increments for a,b,X,Aa,K,Lambda.
        t0 = stage_terms_records[0]["terms"]
        t_primary = stage_terms_records[2]["terms"]
        final_l2 = final_l2_captures[0]["rhs"]
        final_lambda_l2 = final_lambda_captures[0]["rhs"]
        checks = {}
        explicit_fields = ("a", "b", "X")
        for field in explicit_fields:
            expected = 0.5 * (
                _array(t0["explicit"][field]) + _array(t_primary["explicit"][field])
            )
            measured = (
                _array(getattr(accepted.geometry, field)) - _array(getattr(state.geometry, field))
            ) / dt
            checks[f"{field}_weighted_increment_minus_rhs_max_abs"] = _max_abs(measured - expected)
        for field in ("Aa", "K"):
            l2_final = _array(final_l2[field])
            expected = 0.5 * (
                _array(t0["primary_l2"][field])
                + l2_final
                + _array(t0["primary_l3"][field])
                + _array(t_primary["primary_l3"][field])
            )
            measured = (
                _array(getattr(accepted.geometry, field)) - _array(getattr(state.geometry, field))
            ) / dt
            checks[f"{field}_weighted_increment_minus_rhs_max_abs"] = _max_abs(measured - expected)
        expected_lambda = 0.5 * (
            _array(t0["lambda_l2"]) + _array(final_lambda_l2)
            + _array(t0["lambda_l3"]) + _array(t_primary["lambda_l3"])
        )
        measured_lambda = (
            _array(accepted.geometry.Lambda) - _array(state.geometry.Lambda)
        ) / dt
        checks["Lambda_weighted_increment_minus_rhs_max_abs"] = _max_abs(
            measured_lambda - expected_lambda
        )

        for field in ("a", "b", "X", "Aa", "K", "Lambda"):
            value = checks[f"{field}_weighted_increment_minus_rhs_max_abs"]
            implementation_closures.append({"N": resolution, "field": field, "max_abs": value})

        if not bool(np.all(np.isfinite(np.asarray(list(checks.values()), dtype=float)))):
            raise FloatingPointError(f"N={resolution}: non-finite weighted increment closure")

        stage_records = [x["record"] for x in stage_terms_records]
        stage_records.extend(x["record"] for x in final_l2_captures)
        stage_records.extend(x["record"] for x in final_lambda_captures)
        stage_records.extend([accepted_record, post_record])
        if not all(record["all_arrays_finite"] for record in stage_records):
            raise FloatingPointError(
                f"N={resolution}: non-finite captured field or quotient in at least one stage"
            )

        all_profiles.update(
            {f"N{resolution}_final_l2_explicit_state_{key}": value
             for key, value in final_l2_captures[0]["profiles"].items()}
        )
        all_profiles.update(
            {f"N{resolution}_final_lambda_l2_primary_state_{key}": value
             for key, value in final_lambda_captures[0]["profiles"].items()}
        )

        all_reports.append({
            "resolution": resolution,
            "r_max": R_MAX,
            "dr": float(state.grid.dr),
            "dt": dt,
            "cfl": CFL,
            "initial_t": float(state.t),
            "accepted_t": float(accepted.t),
            "input_state_sha256_before": input_digest,
            "input_state_sha256_after": output_digest,
            "input_state_unchanged": input_digest == output_digest,
            "actual_geometry_stage_call_order": [x["record"]["stage"] for x in stage_terms_records],
            "captured_stages": stage_records,
            "weighted_increment_reconstruction": checks,
            "raw_constraint_or_conservation_admission": "NOT_EVALUATED_BY_THIS_DIAGNOSTIC",
            "disposition": (
                "one ordinary production step completed; stage traces localize numerical "
                "compatibility but do not constitute continuum regularity or admission"
            ),
        })

    OUT.mkdir(parents=True, exist_ok=True)
    report = {
        "schema": "v55_pirk_stage_centre_series_compatibility_v1",
        "status": "completed_diagnostic_only",
        "source_commit": SOURCE_COMMIT,
        "source_production_kernel_blob": source_blob,
        "diagnostic_head_commit": head,
        "vendor_commit": vendor,
        "resolutions": list(RESOLUTIONS),
        "r_max": R_MAX,
        "cfl": CFL,
        "fixed_radius_targets": list(FIXED_RADII),
        "centre_cells_recorded": list(CENTRE_CELLS),
        "no_projection_no_physics_edit": True,
        "full_step_attempts_per_resolution": 1,
        "stage_capture_attempts": 3,
        "poststep_rhs_is_read_only": True,
        "weighted_increment_reconstruction_max_abs": implementation_closures,
        "resolution_reports": all_reports,
        "guardrails": [
            "All stages are captured from the exact production kernel implementation at the pinned source commit.",
            "adapter.geometry_stage_terms is wrapped only to record arguments and returned terms; the original result is returned unchanged.",
            "The vendor L2 wrappers record only the direct final-stage calls; calls nested inside geometry_stage_terms are not double-recorded.",
            "The accepted step uses CFL=0.03 and the grid's own dr; no artificial timestep or threshold is imposed.",
            "This diagnostic does not pass/fail a physics tolerance and does not authorize a production repair."
        ],
    }
    report_path = OUT / "report.json"
    report_path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    npz_path = OUT / "profiles.npz"
    np.savez_compressed(npz_path, **all_profiles)

    # Tabular view of every captured profile array, retaining full radial data.
    csv_path = OUT / "profiles.csv"
    rows_written = 0
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        columns = ("profile_key", "i", "r", "value")
        writer.writerow(columns)
        for key in sorted(all_profiles):
            values = np.asarray(all_profiles[key])
            if values.ndim != 1:
                continue
            if key.endswith("_r"):
                continue
            # Associate every vector with its radial grid from the same N/stage.
            prefix = key.rsplit("_", 1)[0]
            # Array names have form N<resolution>_<stage>_<field>. Stage keys can
            # contain underscores; locate the matching r profile by prefix scan.
            parts = key.split("_")
            n_label = parts[0]
            matching = f"{n_label}_"
            candidate_rs = [
                v for k, v in all_profiles.items()
                if k.startswith(matching) and k.endswith("_r")
            ]
            # Canonical r arrays are inserted as stage records' 'r'; their key is
            # reconstructed from the profile key after stripping its final field.
            stage_prefixes = []
            for stage_record in all_reports:
                n_tag = f"N{stage_record['resolution']}_"
                for stage in stage_record["captured_stages"]:
                    if key.startswith(n_tag + stage["stage"] + "_"):
                        stage_prefixes.append(n_tag + stage["stage"] + "_")
                for stage_name in ("accepted_final_state", "rhs_poststep_accepted_state",
                                   "final_l2_explicit_state", "final_lambda_l2_primary_state"):
                    if key.startswith(n_tag + stage_name + "_"):
                        stage_prefixes.append(n_tag + stage_name + "_")
            stage_prefix = max(stage_prefixes, key=len) if stage_prefixes else None
            r_key = f"{stage_prefix}r" if stage_prefix else None
            radii = all_profiles.get(r_key) if r_key else None
            if radii is None or np.asarray(radii).shape != values.shape:
                # Non-profile metadata arrays are excluded from the radial CSV.
                continue
            for i, (radius, value) in enumerate(zip(radii, values)):
                writer.writerow((key, i, float(radius), float(value)))
                rows_written += 1

    (OUT / "SOURCE_COMMIT.txt").write_text(SOURCE_COMMIT + "\n", encoding="utf-8")
    (OUT / "VENDOR_COMMIT.txt").write_text(VENDOR_COMMIT + "\n", encoding="utf-8")
    (OUT / "DIAGNOSTIC_COMMIT.txt").write_text(head + "\n", encoding="utf-8")
    checksum_files = (report_path, npz_path, csv_path, OUT / "SOURCE_COMMIT.txt",
                      OUT / "VENDOR_COMMIT.txt", OUT / "DIAGNOSTIC_COMMIT.txt")
    sums = []
    for path in checksum_files:
        sums.append(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.name}")
    (OUT / "SHA256SUMS.txt").write_text("\n".join(sums) + "\n", encoding="utf-8")

    print("PIRK_STAGE_CENTRE_SERIES_AUDIT_PASS")
    print(f"source={SOURCE_COMMIT} vendor={vendor} diagnostic_head={head}")
    print(f"resolutions={list(RESOLUTIONS)} single_step=True")
    print(f"stage_records={sum(len(x['captured_stages']) for x in all_reports)} csv_profile_rows={rows_written}")
    print(f"report_sha256={hashlib.sha256(report_path.read_bytes()).hexdigest()}")
    print(f"profiles_npz_sha256={hashlib.sha256(npz_path.read_bytes()).hexdigest()}")
    # Do not fail on large or non-convergent physical quotients: those are the
    # findings being measured, not implementation-integrity failures.


if __name__ == "__main__":
    main()
