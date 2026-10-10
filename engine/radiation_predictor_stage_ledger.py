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
from .radiation_geometry_stage_trace import RadiationStageTrace, _safe_float

CELLS = 5
CAPTURE_FROM_T = 44.0


def reconstruct_metric_cone(a, X, U_E, U_r):
    """Independent algebraic reconstruction used by ledger and regressions."""
    a, X, U_E, U_r = map(float, (a, X, U_E, U_r))
    if not all(math.isfinite(v) for v in (a, X, U_E, U_r)) or a <= 0.0:
        raise ValueError("metric/conservative inputs must be finite and a positive")
    gamma_rr_inv = X * X / a
    margin = U_E - abs(U_r) * math.sqrt(gamma_rr_inv)
    ratio = abs(U_r) * math.sqrt(gamma_rr_inv) / U_E if U_E > 0.0 else None
    return {
        "gamma_rr_inv_reconstructed": gamma_rr_inv,
        "cone_margin_reconstructed": margin,
        "ratio_absS_over_E_reconstructed": ratio,
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
                    and int(step.get("rhs_count", 0)) == 0
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
        hashes = {}
        for rel in (
            "engine/production_kernel.py", "engine/v55_pirk_adapter.py",
            "engine/radiation_geometry_stage_trace.py",
            "engine/radiation_predictor_stage_ledger.py",
        ):
            path = Path(rel)
            hashes[rel] = hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None
        payload = {
            "schema_name": "0star_radiation_predictor_metric_stage_ledger_v1",
            "purpose": "diagnostic-only raw-RHS / Euler predictor / projection attribution",
            "status": status,
            "branch": os.environ.get("GITHUB_REF_NAME", "local"),
            "run_trigger_commit": os.environ.get("GITHUB_SHA"),
            "source_file_sha256": hashes,
            "vendor_submodule_commit": self._trace_provenance()["vendor_submodule_commit"]
                if hasattr(self, "_trace_provenance") else None,
            "capture_from_t": CAPTURE_FROM_T,
            "cells": "outermost five cells (75-79 for N=80)",
            "steps": [self.ledger_steps[k] for k in sorted(self.ledger_steps, key=int)],
            "interpretation_guardrail": "A recorded stage attribution is numerical evidence, not a physical-law conclusion.",
        }
        ledger_path = self.output_dir / "stage_ledger.json"
        ledger_path.write_text(json.dumps(payload, indent=2, allow_nan=False))
        summary = {
            "schema_name": payload["schema_name"],
            "status": status,
            "step_count": len(payload["steps"]),
            "captured_stage_counts": {
                "raw_rhs": sum("accepted_geometry_and_rhs" in s for s in payload["steps"]),
                "unprojected_predictor": sum("unprojected_euler_predictor" in s for s in payload["steps"]),
                "post_projection": sum("post_regularity_projection" in s for s in payload["steps"]),
                "first_cmc_inputs": sum("first_predictor_cmc_inputs" in s for s in payload["steps"]),
            },
            "source_file_sha256": hashes,
            "stage_ledger_sha256": hashlib.sha256(ledger_path.read_bytes()).hexdigest(),
            "all_steps_have_expected_order": all(
                ("accepted_geometry_and_rhs" in s
                 and "unprojected_euler_predictor" in s
                 and "post_regularity_projection" in s)
                for s in payload["steps"]
            ),
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
        return result


def main():
    trace = RadiationPredictorStageLedger(Path("runs/radiation-predictor-stage-ledger"))
    return trace.run()


if __name__ == "__main__":
    raise SystemExit(main())
