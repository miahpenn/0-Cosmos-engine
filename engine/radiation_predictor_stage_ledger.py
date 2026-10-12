"""Diagnostic-only ledger separating raw geometry RHS, Euler predictor, and projection.

This module wraps existing production calls; it does not alter the evolved state,
physics, admissibility rules, or step ordering. Run only on an isolated diagnostic
branch after deterministic tests pass.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path

import numpy as np

from . import production_kernel as production_kernel_module
from . import v55_pirk_adapter as adapter
from .radiation_geometry_stage_trace import RadiationStageTrace, _trace_provenance

CELLS = 5
CAPTURE_FROM_T = 44.0


def _grid_size(grid):
    """Return cell count through the pinned grid's public coordinate array."""
    return len(grid.centers)

def _attach_radiation_state(values, radiation):
    """Attach the exact predictor conservative state and recompute cone metrics."""
    for cell in values:
        i = int(cell["i"])
        U_E = float(radiation.energy_t[i])
        U_r = float(radiation.momentum_r[i])
        cell["U_E"] = U_E
        cell["U_r"] = U_r
        cell.update(reconstruct_metric_cone(cell["a"], cell["X"], U_E, U_r))
    return values


def _capture_first_geometry_terms_call(step):
    """Capture only gterms0, not later predictor geometry evaluations."""
    if float(step.get("t0", 0.0)) < CAPTURE_FROM_T:
        return False
    count = int(step.get("ledger_geometry_terms_count", 0))
    step["ledger_geometry_terms_count"] = count + 1
    return count == 0 and int(step.get("rhs_count", 0)) == 1


def reconstruct_metric_cone(a, X, U_E, U_r):
    """Reconstruct the radial inverse metric and radiation cone invariant."""
    a, X, U_E, U_r = map(float, (a, X, U_E, U_r))
    if not all(math.isfinite(v) for v in (a, X, U_E, U_r)):
        raise ValueError("metric/conservative inputs must be finite")
    if a <= 0.0 or X <= 0.0:
        raise ValueError("radial metric inputs a and X must be positive")
    if U_E <= 0.0:
        raise ValueError("radiation conservative energy must be positive")
    gamma_rr_inv = X * X / a
    margin = U_E - abs(U_r) * math.sqrt(gamma_rr_inv)
    ratio = abs(U_r) * math.sqrt(gamma_rr_inv) / U_E
    return {
        "gamma_rr_inv_reconstructed": gamma_rr_inv,
        "cone_margin_reconstructed": margin,
        "ratio_absS_over_E_reconstructed": ratio,
    }


def validate_ledger_steps(steps, expected_indices):
    """Machine-check that every recorded step contains a closed stage ledger."""
    expected = [int(i) for i in expected_indices]
    errors = []

    if not steps:
        return {
            "ok": False,
            "step_count": 0,
            "expected_cell_indices": expected,
            "errors": ["no predictor-stage ledger steps were captured"],
        }

    stages = (
        ("accepted_geometry_and_rhs", "accepted_geometry"),
        ("unprojected_euler_predictor", None),
        ("post_regularity_projection", None),
        ("first_predictor_cmc_inputs", None),
    )
    metric_fields = ("a", "b", "X", "alpha", "beta")

    def rows_have_expected_indices(rows, where):
        if not isinstance(rows, list):
            errors.append(f"{where}: expected a list of cell records")
            return False
        observed = [int(row.get("i", -1)) for row in rows]
        if observed != expected:
            errors.append(f"{where}: cell indices {observed!r} != {expected!r}")
            return False
        return True

    for step in steps:
        step_id = step.get("step_index", "?")
        prefix = f"step {step_id}"
        accepted_block = step.get("accepted_geometry_and_rhs")
        accepted = (
            accepted_block.get("accepted_geometry")
            if isinstance(accepted_block, dict) else None
        )
        stage_values = {
            "accepted_geometry_and_rhs": accepted,
            "unprojected_euler_predictor": step.get("unprojected_euler_predictor"),
            "post_regularity_projection": step.get("post_regularity_projection"),
            "first_predictor_cmc_inputs": step.get("first_predictor_cmc_inputs"),
        }

        for stage_name, _ in stages:
            values = stage_values[stage_name]
            if not values:
                errors.append(f"{prefix}: missing stage {stage_name}")
                continue
            rows_have_expected_indices(values, f"{prefix}/{stage_name}")
            for cell in values:
                for field in metric_fields:
                    value = cell.get(field)
                    if not isinstance(value, (int, float)) or not math.isfinite(float(value)):
                        errors.append(f"{prefix}/{stage_name}/cell {cell.get('i')}: invalid {field}")
                    elif field in ("a", "b", "X") and float(value) <= 0.0:
                        errors.append(f"{prefix}/{stage_name}/cell {cell.get('i')}: non-positive {field}")
                for field in (
                    "U_E", "U_r", "gamma_rr_inv_reconstructed",
                    "cone_margin_reconstructed", "ratio_absS_over_E_reconstructed",
                ):
                    value = cell.get(field)
                    if not isinstance(value, (int, float)) or not math.isfinite(float(value)):
                        errors.append(f"{prefix}/{stage_name}/cell {cell.get('i')}: invalid or missing {field}")
                if all(isinstance(cell.get(k), (int, float)) and math.isfinite(float(cell[k]))
                       for k in ("a", "X", "U_E", "U_r")) and float(cell["a"]) > 0 and float(cell["X"]) > 0 and float(cell["U_E"]) > 0:
                    try:
                        reconstructed = reconstruct_metric_cone(
                            cell["a"], cell["X"], cell["U_E"], cell["U_r"]
                        )
                        for field, expected_value in reconstructed.items():
                            if not math.isclose(
                                float(cell[field]), expected_value,
                                rel_tol=1.0e-12, abs_tol=1.0e-15
                            ):
                                errors.append(
                                    f"{prefix}/{stage_name}/cell {cell.get('i')}: {field} does not reconstruct"
                                )
                    except (ValueError, TypeError, KeyError):
                        errors.append(f"{prefix}/{stage_name}/cell {cell.get('i')}: reconstruction failed")

        rhs = accepted_block.get("raw_explicit_geometry_rhs") if isinstance(accepted_block, dict) else None
        if not isinstance(rhs, dict):
            errors.append(f"{prefix}: missing raw explicit geometry RHS")
            continue
        rhs_by_var = {}
        for variable in ("a", "b", "X"):
            records = rhs.get(variable)
            if not rows_have_expected_indices(records, f"{prefix}/raw_rhs/{variable}"):
                continue
            rhs_by_var[variable] = {}
            for record in records:
                value = record.get("value")
                if not isinstance(value, (int, float)) or not math.isfinite(float(value)):
                    errors.append(f"{prefix}/raw_rhs/{variable}/cell {record.get('i')}: invalid value")
                else:
                    rhs_by_var[variable][int(record["i"])] = float(value)

        if accepted and step.get("unprojected_euler_predictor") and rhs_by_var.keys() == {"a", "b", "X"}:
            unprojected = {int(c["i"]): c for c in step["unprojected_euler_predictor"]}
            accepted_map = {int(c["i"]): c for c in accepted}
            dt = float(step.get("dt", math.nan))
            if not math.isfinite(dt) or dt <= 0:
                errors.append(f"{prefix}: invalid dt for Euler-predictor reconstruction")
            else:
                for i in expected:
                    for variable in ("a", "b", "X"):
                        predicted = float(accepted_map[i][variable]) + dt * rhs_by_var[variable][i]
                        if not math.isclose(
                            float(unprojected[i][variable]), predicted,
                            rel_tol=2.0e-12, abs_tol=2.0e-14
                        ):
                            errors.append(f"{prefix}/cell {i}: unprojected {variable} is not accepted + dt*RHS")

        projected = stage_values["post_regularity_projection"]
        if projected:
            for cell in projected:
                if all(isinstance(cell.get(k), (int, float)) and math.isfinite(float(cell[k]))
                       for k in ("a", "b")) and not math.isclose(
                           float(cell["a"]) * float(cell["b"]) ** 2,
                           1.0, rel_tol=1.0e-10, abs_tol=1.0e-12
                       ):
                    errors.append(f"{prefix}/cell {cell.get('i')}: post-projection a*b^2 != 1")

        first_cmc = stage_values["first_predictor_cmc_inputs"]
        if projected and first_cmc:
            projected_map = {int(c["i"]): c for c in projected}
            cmc_map = {int(c["i"]): c for c in first_cmc}
            for i in expected:
                for variable in ("a", "b", "X"):
                    if not math.isclose(
                        float(projected_map[i][variable]), float(cmc_map[i][variable]),
                        rel_tol=0.0, abs_tol=0.0
                    ):
                        errors.append(f"{prefix}/cell {i}: first CMC input {variable} differs from post-projection metric")

    return {
        "ok": not errors,
        "step_count": len(steps),
        "expected_cell_indices": expected,
        "errors": errors,
    }




def validate_primary_source_capture(steps, expected_indices, failure_step_index=None):
    """Fail closed on missing/misaligned Aa/K source-term captures."""
    errors = []
    max_split_error = 0.0
    max_transition_error = 0.0
    expected_indices = list(expected_indices)
    labels = (
        "accepted_start_gterms0",
        "explicit_predictor_gterms_pred",
        "primary_predictor_gterms1",
    )
    previous = None
    for step in steps:
        step_index = int(step.get("step_index", -1))
        is_failure = failure_step_index is not None and step_index == int(failure_step_index)
        stages = step.get("primary_source_stages")
        expected_labels = [labels[0]] if is_failure else list(labels)
        if not isinstance(stages, list):
            errors.append(f"step {step_index}: primary_source_stages missing")
            stages = []
        got_labels = [str(s.get("stage_label")) for s in stages]
        if got_labels != expected_labels:
            errors.append(f"step {step_index}: stage labels {got_labels!r} != {expected_labels!r}")
        for stage in stages:
            rows = stage.get("cells", [])
            got_indices = [int(r.get("i", -1)) for r in rows]
            if got_indices != expected_indices:
                errors.append(f"step {step_index}/{stage.get('stage_label')}: cells {got_indices!r} != {expected_indices!r}")
            for cell in rows:
                for block in ("primary_l2", "primary_l3_total", "primary_l3_vacuum",
                              "primary_l3_matter_delta", "explicit_rhs"):
                    values = cell.get(block, {})
                    required = ("Aa", "K") if block != "explicit_rhs" else ("a", "b", "X")
                    if any(k not in values or not math.isfinite(float(values[k])) for k in required):
                        errors.append(f"step {step_index}/{stage.get('stage_label')}/cell {cell.get('i')}: invalid {block}")
                for key in ("Aa", "K"):
                    total = float(cell["primary_l3_total"][key])
                    vacuum = float(cell["primary_l3_vacuum"][key])
                    matter_delta = float(cell["primary_l3_matter_delta"][key])
                    max_split_error = max(max_split_error, abs(total - vacuum - matter_delta))
        if not is_failure:
            final_l2 = step.get("final_explicit_state_primary_l2")
            if not isinstance(final_l2, dict):
                errors.append(f"step {step_index}: final_explicit_state_primary_l2 missing")
            else:
                for section in ("geometry",):
                    got = [int(r.get("i", -1)) for r in final_l2.get(section, [])]
                    if got != expected_indices:
                        errors.append(f"step {step_index}/final_explicit_state_primary_l2/{section}: cells {got!r} != {expected_indices!r}")
                for key in ("Aa", "K"):
                    got = [int(r.get("i", -1)) for r in final_l2.get("rhs", {}).get(key, [])]
                    if got != expected_indices:
                        errors.append(f"step {step_index}/final_explicit_state_primary_l2/{key}: cells {got!r} != {expected_indices!r}")
            for name in ("primary_predictor_after_outer_boundary", "completed_state_after_outer_boundary"):
                rows = step.get(name)
                if not isinstance(rows, list) or [int(r.get("i", -1)) for r in rows] != expected_indices:
                    errors.append(f"step {step_index}: {name} missing or cell mismatch")
        if previous is not None:
            previous_step, previous_rows = previous
            if previous_step + 1 == step_index and not is_failure:
                new_accepted = {
                    int(r["i"]): r for r in step.get("accepted_geometry_and_rhs", {}).get("accepted_geometry", [])
                }
                for old in previous_rows:
                    i = int(old["i"])
                    if i not in new_accepted:
                        continue
                    for field in ("a", "b", "X", "Aa", "K"):
                        delta = abs(float(old[field]) - float(new_accepted[i][field]))
                        max_transition_error = max(max_transition_error, delta)
                        if delta > 1.0e-10:
                            errors.append(f"step {previous_step}->{step_index}/cell {i}/{field}: completed-to-accepted mismatch {delta:.6g}")
        completed = step.get("completed_state_after_outer_boundary")
        if completed:
            previous = (step_index, completed)
        else:
            previous = (step_index, [])
    return {
        "ok": not errors and max_split_error <= 1.0e-12 and max_transition_error <= 1.0e-10,
        "step_count": len(steps),
        "expected_cell_indices": expected_indices,
        "max_abs_primary_l3_split_error": max_split_error,
        "max_abs_completed_to_next_accepted_error": max_transition_error,
        "errors": errors,
    }

def _primary_source_term_record(grid, geometry, terms, vacuum_l3, call_index):
    """Serialize existing Aa/K L2/L3 terms for the outer capture cells."""
    n = len(grid.centers)
    lo = max(0, n - CELLS)
    labels = (
        "accepted_start_gterms0",
        "explicit_predictor_gterms_pred",
        "primary_predictor_gterms1",
    )
    record = {
        "call_index": int(call_index),
        "stage_label": labels[call_index] if call_index < len(labels)
            else f"geometry_terms_call_{call_index + 1}",
        "cells": [],
    }
    for i in range(lo, n):
        l2 = terms["primary_l2"]
        l3 = terms["primary_l3"]
        row = {
            "i": i,
            "r": float(grid.centers[i]),
            "geometry": {
                key: float(getattr(geometry, key)[i])
                for key in ("a", "b", "X", "alpha", "Aa", "K", "Lambda")
                if hasattr(geometry, key)
            },
            "primary_l2": {
                key: float(np.asarray(l2[key])[i]) for key in ("Aa", "K")
            },
            "primary_l3_total": {
                key: float(np.asarray(l3[key])[i]) for key in ("Aa", "K")
            },
            "primary_l3_vacuum": {
                key: float(np.asarray(vacuum_l3[key])[i]) for key in ("Aa", "K")
            },
            "primary_l3_matter_delta": {
                key: float(np.asarray(l3[key])[i] - np.asarray(vacuum_l3[key])[i])
                for key in ("Aa", "K")
            },
            "explicit_rhs": {
                key: float(np.asarray(terms["explicit"][key])[i])
                for key in ("a", "b", "X")
            },
        }
        record["cells"].append(row)
    return record

class RadiationPredictorStageLedger(RadiationStageTrace):
    """Extend the existing source trace with the missing pre/post-projection stages."""

    def __init__(self, output_dir: Path):
        self.ledger_steps: dict[str, dict] = {}
        super().__init__(output_dir)
        kernel_cls = production_kernel_module.V55ProductionKernel
        self._ledger_original_project = kernel_cls._enforce_center_regularity
        self._ledger_original_solve = kernel_cls._solve_lapse
        self._ledger_original_geometry_terms = adapter.geometry_stage_terms
        _, vacuum, _ = adapter.vendor_modules()
        self._ledger_vacuum_module = vacuum
        self._ledger_original_primary_l2 = vacuum.primary_l2_rhs
        self._ledger_original_primary_l3 = vacuum.primary_l3_rhs
        ledger = self

        def capture_project(grid, geometry):
            step = ledger.current_step
            if step is None:
                return ledger._ledger_original_project(grid, geometry)
            n = int(step.setdefault("ledger_project_count", 0))
            step["ledger_project_count"] = n + 1
            if n != 1 or float(step.get("t0", 0.0)) < CAPTURE_FROM_T:
                return ledger._ledger_original_project(grid, geometry)
            raw = geometry.copy()
            projected = ledger._ledger_original_project(grid, geometry)
            row = ledger._row_for_step(step)
            row["unprojected_euler_predictor"] = ledger._geometry_values(
                grid, raw, None
            )
            row["post_regularity_projection"] = ledger._geometry_values(
                grid, projected, None
            )
            return projected

        def capture_solve(grid, geometry, scalars, matter):
            step = ledger.current_step
            solve_index = int(step.get("solve_count", 0)) if step is not None else -1
            if step is not None and float(step.get("t0", 0.0)) >= CAPTURE_FROM_T:
                row = ledger._row_for_step(step)
                if solve_index == 1:
                    row["first_predictor_cmc_inputs"] = ledger._geometry_values(
                        grid, geometry, matter.radiation
                    )
                    raw_values = row.get("unprojected_euler_predictor", [])
                    projected_values = row.get("post_regularity_projection", [])
                    accepted_values = row.get("accepted_geometry_and_rhs", {}).get(
                        "accepted_geometry", []
                    )
                    # Keep the actual predictor conservatives paired with every metric.
                    for values in (accepted_values, raw_values, projected_values):
                        if values:
                            _attach_radiation_state(values, matter.radiation)
                elif solve_index == 2:
                    row["primary_predictor_after_outer_boundary"] = (
                        ledger._geometry_values(grid, geometry, matter.radiation)
                    )
                elif solve_index == 3:
                    row["completed_state_after_outer_boundary"] = (
                        ledger._geometry_values(grid, geometry, matter.radiation)
                    )
            return ledger._ledger_original_solve(grid, geometry, scalars, matter)

        def capture_primary_l2(grid, geometry):
            result = ledger._ledger_original_primary_l2(grid, geometry)
            step = ledger.current_step
            if (step is not None
                    and not step.get("ledger_inside_geometry_terms", False)
                    and float(step.get("t0", 0.0)) >= CAPTURE_FROM_T):
                row = ledger._row_for_step(step)
                n = len(grid.centers)
                lo = max(0, n - CELLS)
                row["final_explicit_state_primary_l2"] = {
                    "geometry": [
                        {
                            "i": i, "r": float(grid.centers[i]),
                            **{
                                key: float(getattr(geometry, key)[i])
                                for key in ("a", "b", "X", "alpha", "Aa", "K", "Lambda")
                                if hasattr(geometry, key)
                            },
                        }
                        for i in range(lo, n)
                    ],
                    "rhs": {
                        key: [
                            {"i": i, "r": float(grid.centers[i]),
                             "value": float(np.asarray(result[key])[i])}
                            for i in range(lo, n)
                        ]
                        for key in ("Aa", "K")
                    },
                }
            return result

        def capture_geometry_terms(grid, geometry, fields, matter, lambda_m=2.0):
            step = ledger.current_step
            active = step is not None and float(step.get("t0", 0.0)) >= CAPTURE_FROM_T
            if active:
                step["ledger_inside_geometry_terms"] = True
            try:
                terms = ledger._ledger_original_geometry_terms(
                    grid, geometry, fields, matter, lambda_m=lambda_m
                )
            finally:
                if active:
                    step["ledger_inside_geometry_terms"] = False
            if active:
                row = ledger._row_for_step(step)
                call_index = int(step.get("ledger_primary_terms_call_count", 0))
                step["ledger_primary_terms_call_count"] = call_index + 1
                vacuum_l3 = ledger._ledger_original_primary_l3(grid, geometry)
                row.setdefault("primary_source_stages", []).append(
                    _primary_source_term_record(
                        grid, geometry, terms, vacuum_l3, call_index
                    )
                )
            if step is not None and _capture_first_geometry_terms_call(step):
                row = ledger._row_for_step(step)
                row["accepted_geometry_and_rhs"] = {
                    "accepted_geometry": ledger._geometry_values(
                        grid, geometry, None
                    ),
                    "raw_explicit_geometry_rhs": {
                        name: [
                            {"i": i, "r": float(grid.centers[i]),
                             "value": float(np.asarray(terms["explicit"][name])[i])}
                            for i in range(max(0, len(grid.centers) - CELLS), len(grid.centers))
                        ]
                        for name in ("a", "b", "X")
                    },
                }
            return terms

        kernel_cls._enforce_center_regularity = staticmethod(capture_project)
        kernel_cls._solve_lapse = staticmethod(capture_solve)
        adapter.geometry_stage_terms = capture_geometry_terms
        vacuum.primary_l2_rhs = capture_primary_l2

    def _row_for_step(self, step):
        key = str(int(step["step_index"]))
        row = self.ledger_steps.setdefault(key, {
            "step_index": int(step["step_index"]),
            "t_accepted": float(step["t0"]),
            "t_attempted": float(step["t0"] + step["dt"]),
            "dt": float(step["dt"]),
            "cells": list(range(max(0, _grid_size(self.latest_state.grid) - CELLS),
                                _grid_size(self.latest_state.grid)))
                if self.latest_state is not None else [],
        })
        return row

    @staticmethod
    def _geometry_values(grid, geometry, radiation):
        n = len(grid.centers)
        lo = max(0, n - CELLS)
        rows = []
        for i in range(lo, n):
            row = {
                "i": i, "r": float(grid.centers[i]),
                "a": float(geometry.a[i]), "b": float(geometry.b[i]),
                "X": float(geometry.X[i]), "alpha": float(geometry.alpha[i]),
                "beta": float(geometry.beta[i]),
            }
            for field in ("Aa", "K", "Lambda"):
                if hasattr(geometry, field):
                    row[field] = float(getattr(geometry, field)[i])
            if radiation is not None:
                U_E, U_r = float(radiation.energy_t[i]), float(radiation.momentum_r[i])
                row.update({"U_E": U_E, "U_r": U_r})
                row.update(reconstruct_metric_cone(row["a"], row["X"], U_E, U_r))
            rows.append(row)
        return rows

    def write_outputs(self, state, status, dt_nominal):
        result = super().write_outputs(state, status, dt_nominal)
        trace_path = self.output_dir / "trace.json"
        trace_payload = json.loads(trace_path.read_text()) if trace_path.is_file() else {}
        trace_provenance = trace_payload.get("provenance", {})
        configuration = trace_provenance.get("configuration", {})
        provenance = _trace_provenance()
        hashes = {}
        for rel in (
            "engine/production_kernel.py",
            "engine/cmc_gauge.py",
            "engine/matter_system.py",
            "engine/matter_rhs.py",
            "engine/v55_matter.py",
            "engine/v55_pirk_adapter.py",
            "engine/radiation_geometry_stage_trace.py",
            "engine/radiation_predictor_stage_ledger.py",
            ".github/workflows/zero_star_radiation_stage_trace.yml",
            ".github/workflows/radiation-predictor-metric-stage-ledger.yml",
            "tests/test_radiation_predictor_stage_ledger.py",
        ):
            path = Path(rel)
            hashes[rel] = hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None

        centers = np.asarray(state.grid.centers, dtype="<f8") if state is not None else np.asarray([], dtype="<f8")
        grid_centers_sha256 = hashlib.sha256(centers.tobytes()).hexdigest()
        config_payload = {
            "trace_configuration": configuration,
            "resolution": int(_grid_size(state.grid)) if state is not None else None,
            "dr": float(state.grid.dr) if state is not None else None,
            "dt_nominal": float(dt_nominal) if math.isfinite(float(dt_nominal)) else None,
            "capture_from_t": CAPTURE_FROM_T,
            "captured_cell_count": CELLS,
            "grid_centers_sha256": grid_centers_sha256,
        }
        canonical_config = json.dumps(
            config_payload, sort_keys=True, separators=(",", ":"), allow_nan=False
        ).encode()
        configuration_sha256 = hashlib.sha256(canonical_config).hexdigest()
        input_fingerprint = {
            "configuration_sha256": configuration_sha256,
            "grid_centers_sha256": grid_centers_sha256,
            "source_file_sha256": hashes,
            "vendor_submodule_commit": provenance.get("vendor_submodule_commit"),
            "run_trigger_commit": os.environ.get("GITHUB_SHA"),
        }
        input_sha256 = hashlib.sha256(json.dumps(
            input_fingerprint, sort_keys=True, separators=(",", ":"), allow_nan=False
        ).encode()).hexdigest()

        expected_indices = (
            list(range(max(0, int(_grid_size(state.grid)) - CELLS), int(_grid_size(state.grid))))
            if state is not None else []
        )
        steps = [self.ledger_steps[k] for k in sorted(self.ledger_steps, key=int)]
        validation = validate_ledger_steps(steps, expected_indices)
        failure_step_index = trace_payload.get("failure", {}).get("step_index")
        primary_source_validation = validate_primary_source_capture(
            steps, expected_indices, failure_step_index=failure_step_index
        )
        payload = {
            "schema_name": "0star_radiation_predictor_metric_stage_ledger_v2",
            "purpose": "diagnostic-only raw-RHS / Euler predictor / projection attribution",
            "status": status,
            "branch": os.environ.get("GITHUB_REF_NAME", "local"),
            "run_trigger_commit": os.environ.get("GITHUB_SHA"),
            "source_commit_before_trace_workflow": provenance.get("source_commit_before_trace_workflow"),
            "instrumentation_commit": provenance.get("instrumentation_commit"),
            "workflow_commit": provenance.get("workflow_commit"),
            "source_file_sha256": hashes,
            "vendor_submodule_commit": provenance.get("vendor_submodule_commit"),
            "configuration": config_payload,
            "configuration_sha256": configuration_sha256,
            "input_fingerprint_sha256": input_sha256,
            "input_fingerprint": input_fingerprint,
            "capture_from_t": CAPTURE_FROM_T,
            "cells": expected_indices,
            "steps": steps,
            "validation": validation,
            "primary_source_capture_validation": primary_source_validation,
            "interpretation_guardrail": "A recorded stage attribution is numerical evidence, not a physical-law conclusion.",
            "trace_sha256": hashlib.sha256(trace_path.read_bytes()).hexdigest() if trace_path.is_file() else None,
        }
        ledger_path = self.output_dir / "stage_ledger.json"
        ledger_path.write_text(json.dumps(payload, indent=2, allow_nan=False))
        summary = {
            "schema_name": payload["schema_name"],
            "status": status,
            "step_count": len(steps),
            "captured_stage_counts": {
                "raw_rhs": sum("accepted_geometry_and_rhs" in s for s in steps),
                "unprojected_predictor": sum("unprojected_euler_predictor" in s for s in steps),
                "post_projection": sum("post_regularity_projection" in s for s in steps),
                "first_cmc_inputs": sum("first_predictor_cmc_inputs" in s for s in steps),
            },
            "source_file_sha256": hashes,
            "configuration_sha256": configuration_sha256,
            "input_fingerprint_sha256": input_sha256,
            "trace_sha256": payload["trace_sha256"],
            "stage_ledger_sha256": hashlib.sha256(ledger_path.read_bytes()).hexdigest(),
            "admission": validation,
            "primary_source_capture_admission": primary_source_validation,
            "all_steps_have_expected_order": validation["ok"] and primary_source_validation["ok"],
        }
        summary_path = self.output_dir / "stage_ledger_summary.json"
        summary_path.write_text(json.dumps(summary, indent=2, allow_nan=False))
        checksum_path = self.output_dir / "SHA256SUMS.txt"
        with checksum_path.open("a") as stream:
            stream.write(
                hashlib.sha256(ledger_path.read_bytes()).hexdigest() + "  stage_ledger.json\n"
                + hashlib.sha256(summary_path.read_bytes()).hexdigest() + "  stage_ledger_summary.json\n"
            )
        print(json.dumps({"stage_ledger_summary": summary}, indent=2), flush=True)
        return result if validation["ok"] and primary_source_validation["ok"] else 1


def main():
    trace = RadiationPredictorStageLedger(Path("runs/radiation-predictor-stage-ledger"))
    return trace.run()


if __name__ == "__main__":
    raise SystemExit(main())
