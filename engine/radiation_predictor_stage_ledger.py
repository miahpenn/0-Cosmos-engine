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


class RadiationPredictorStageLedger(RadiationStageTrace):
    """Extend the existing source trace with the missing pre/post-projection stages."""

    def __init__(self, output_dir: Path):
        self.ledger_steps: dict[str, dict] = {}
        super().__init__(output_dir)
        kernel_cls = production_kernel_module.V55ProductionKernel
        self._ledger_original_project = kernel_cls._enforce_center_regularity
        self._ledger_original_solve = kernel_cls._solve_lapse
        self._ledger_original_geometry_terms = adapter.geometry_stage_terms
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
            if (step is not None
                    and int(step.get("solve_count", 0)) == 1
                    and float(step.get("t0", 0.0)) >= CAPTURE_FROM_T):
                row = ledger._row_for_step(step)
                row["first_predictor_cmc_inputs"] = ledger._geometry_values(
                    grid, geometry, matter.radiation
                )
                raw_values = row.get("unprojected_euler_predictor", [])
                projected_values = row.get("post_regularity_projection", [])
                accepted_values = row.get("accepted_geometry_and_rhs", {}).get(
                    "accepted_geometry", []
                )
                # Independently recompute the invariant from primitive a/X and
                # conservative U_E/U_r for all recorded intermediate metrics.
                for stage_name, values in (
                    ("accepted_metric", accepted_values),
                    ("unprojected_metric", raw_values),
                    ("post_projection_metric", projected_values),
                ):
                    if not values:
                        continue
                    for cell in values:
                        i = cell["i"]
                        rad = matter.radiation
                        cell.update(reconstruct_metric_cone(
                            cell["a"], cell["X"],
                            rad.energy_t[i], rad.momentum_r[i],
                        ))
            return ledger._ledger_original_solve(grid, geometry, scalars, matter)

        def capture_geometry_terms(grid, geometry, fields, matter, lambda_m=2.0):
            terms = ledger._ledger_original_geometry_terms(
                grid, geometry, fields, matter, lambda_m=lambda_m
            )
            step = ledger.current_step
            if (step is not None
                    and int(step.get("rhs_count", 0)) == 1
                    and float(step.get("t0", 0.0)) >= CAPTURE_FROM_T):
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

    def _row_for_step(self, step):
        key = str(int(step["step_index"]))
        row = self.ledger_steps.setdefault(key, {
            "step_index": int(step["step_index"]),
            "t_accepted": float(step["t0"]),
            "t_attempted": float(step["t0"] + step["dt"]),
            "dt": float(step["dt"]),
            "cells": list(range(max(0, self.latest_state.grid.resolution - CELLS),
                                self.latest_state.grid.resolution))
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
            ".github/workflows/zero_star_radiation_predictor_stage_ledger.yml",
            "tests/test_radiation_predictor_stage_ledger.py",
        ):
            path = Path(rel)
            hashes[rel] = hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None

        centers = np.asarray(state.grid.centers, dtype="<f8") if state is not None else np.asarray([], dtype="<f8")
        grid_centers_sha256 = hashlib.sha256(centers.tobytes()).hexdigest()
        config_payload = {
            "trace_configuration": configuration,
            "resolution": int(state.grid.resolution) if state is not None else None,
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
            list(range(max(0, int(state.grid.resolution) - CELLS), int(state.grid.resolution)))
            if state is not None else []
        )
        steps = [self.ledger_steps[k] for k in sorted(self.ledger_steps, key=int)]
        validation = validate_ledger_steps(steps, expected_indices)
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
            "all_steps_have_expected_order": validation["ok"],
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
        return result if validation["ok"] else 1


def main():
    trace = RadiationPredictorStageLedger(Path("runs/radiation-predictor-stage-ledger"))
    return trace.run()


if __name__ == "__main__":
    raise SystemExit(main())
