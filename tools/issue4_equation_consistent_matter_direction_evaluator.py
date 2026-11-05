#!/usr/bin/env python3
"""Run the locked Issue #4 predictor on the separate equation-consistent control."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from tools.issue4_equation_consistent_matter_direction import STATE_ARRAY_NAMES
from tools.issue4_equation_consistent_matter_direction_fixture import load_manufactured_fixture
from tools.issue4_zero_step_full_bianchi_observer import verify_provenance
from tools.issue4_full_bianchi_evaluator import _evaluate_loaded_fixture, _jsonable


def evaluate_equation_consistent_fixture(archive_path, resolution, diagnostic_root, helper_baseline_root):
    if resolution != 60:
        raise ValueError("equation-consistent matter-direction fixture is pinned to N=60")
    diagnostic_root = Path(diagnostic_root).resolve()
    helper_baseline_root = Path(helper_baseline_root).resolve()
    provenance = verify_provenance(diagnostic_root, helper_baseline_root)
    loaded = load_manufactured_fixture(archive_path)
    report, arrays = loaded["report"], loaded["arrays"]
    state_arrays = {
        name: np.asarray(arrays[f"state.base.{name}"], dtype=np.float64)
        for name in STATE_ARRAY_NAMES
    }
    records = report["raw_arrays"]["array_records"]
    paths = {record["path"] for record in records}
    required = {
        "frozen_input_state": {f"state.base.{name}" for name in STATE_ARRAY_NAMES},
        "analytic_family_tangent": {"tangent_eps2e-05.Mdot", "tangent_eps1e-05.Mdot"},
        "pointwise_total_stress_divergence": {
            "total_T_eps2e-05.total.div_T_r", "total_T_eps1e-05.total.div_T_r",
        },
        "finite_volume_total_stress_divergence": {
            "total_T_eps2e-05.total.div_T_r_fv", "total_T_eps1e-05.total.div_T_r_fv",
        },
    }
    missing = {key: sorted(value - paths) for key, value in required.items() if value - paths}
    if missing:
        raise RuntimeError(f"equation-consistent fixture witnesses are incomplete: {missing}")
    structural_inventory = {
        "lens": "equation_consistent_matter_direction_structural_inventory",
        "stored_operand_families": {key: sorted(value & paths) for key, value in required.items()},
        "missing_required_witnesses": missing,
        "interpretation": (
            "Equation-consistent matter-direction states, six-sector projections, "
            "central-difference tangent, and independent pointwise/FV divergence "
            "channels are stored in the verified archive. Geometry remains frozen."
        ),
    }
    fixture = {
        "resolution": resolution,
        "accepted_time": 0.0,
        "state_arrays": state_arrays,
        "arrays": arrays,
        "report": report,
        "helper_state_digest": report["state_digests"]["base"],
        "archive_state_digest": report["state_digests"]["base"],
        "structural_inventory": structural_inventory,
        "state_class": "equation_consistent_matter_direction",
        "tangent_semantics": report["tangent_semantics"],
    }
    return _evaluate_loaded_fixture(
        fixture, resolution, diagnostic_root, helper_baseline_root, provenance,
        state_class="equation_consistent_matter_direction",
        tangent_semantics=report["tangent_semantics"],
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--resolution", type=int, required=True)
    parser.add_argument("--diagnostic-source", type=Path, required=True)
    parser.add_argument("--helper-baseline", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = evaluate_equation_consistent_fixture(
        args.archive, args.resolution, args.diagnostic_source, args.helper_baseline
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(_jsonable(result), indent=2, sort_keys=True, allow_nan=False) + "\n")
    print(json.dumps({
        "status": result["status"],
        "state_class": result["state_class"],
        "tangent_semantics": result["tangent_semantics"],
        "production_steps_run": result["production_steps_run"],
        "state_advanced": result["state_advanced"],
        "zero_step_tripwire": result["zero_step_tripwire"],\n        "summaries": _jsonable(result["summaries"]),
        "output": str(args.output),
    }, indent=2))


if __name__ == "__main__":
    main()
