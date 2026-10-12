"""Known-physics check: [flat-vacuum-geometry-smoke] production geometry operators on exact flat vacuum.

This probes the production adapter's geometry RHS blocks, not the full
PIRK stepper. Exact flat spherical BSSN data with all stress-energy zero must
produce zero explicit, primary-L2/L3, and Lambda-L2/L3 RHS to roundoff.
[flat-vacuum-geometry-smoke]
"""
from __future__ import annotations

import json
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from engine import v55_pirk_adapter as adapter
from engine.scalar_system import ScalarFields
from engine.matter_system import ConservedSpecies
from engine.v55_matter import V55MatterState

OUT = ROOT / "runs" / "flat_vacuum_geometry_smoke"


def main():
    grid_ops, vacuum, _ = adapter.vendor_modules()
    grid = grid_ops.SphericalCellGrid(64, 32.0)
    geom = vacuum.flat_state(grid)
    for name, value in (
        ("alpha", 1.0), ("beta", 0.0), ("B", 0.0),
        ("a", 1.0), ("b", 1.0), ("X", 1.0),
        ("Aa", 0.0), ("K", 0.0), ("Lambda", 0.0),
    ):
        getattr(geom, name)[:] = value

    zeros = np.zeros(grid.n, dtype=float)
    fields = ScalarFields(*(zeros.copy() for _ in range(6)))
    # The current Valencia dust recovery rejects exactly zero rest-current.
    # Use a positive 1e-30 admissible seed so the adapter can evaluate its
    # matter projections; its source contribution is negligible at the
    # stated 1e-12 RHS tolerance and is explicitly recorded as near-vacuum.
    tiny = np.full(grid.n, 1.0e-30, dtype=float)
    fluid = ConservedSpecies(tiny.copy(), tiny.copy(), zeros.copy())
    matter = V55MatterState(fluid.copy(), fluid.copy(), fluid.copy())
    terms = adapter.geometry_stage_terms(grid, geom, fields, matter)

    groups = {
        "explicit": terms["explicit"],
        "primary_l2": terms["primary_l2"],
        "primary_l3": terms["primary_l3"],
        "lambda_l2": {"Lambda": terms["lambda_l2"]},
        "lambda_l3": {"Lambda": terms["lambda_l3"]},
    }
    maxima = {}
    for group, values in groups.items():
        maxima[group] = {
            name: float(np.max(np.abs(np.asarray(value, dtype=float))))
            for name, value in values.items()
        }
    overall = max(v for group in maxima.values() for v in group.values())
    finite = all(np.isfinite(v) for group in maxima.values() for v in group.values())
    # This is a roundoff-level static-state test, not a tolerance fit.
    tolerance = 1.0e-12
    passed = finite and overall <= tolerance
    payload = {
        "test": "near-vacuum flat state through production geometry RHS adapter",
        "scope": "geometry RHS blocks only; tiny positive fluid seed required by current matter API; does not exercise full production stepper",
        "grid": {"resolution": grid.n, "r_max": 32.0, "dr": float(grid.dr)},
        "state": "a=b=X=alpha=1; beta=B=Aa=K=Lambda=0; scalar fields zero; each fluid species has rest=energy_t=1e-30 and zero momentum (near-vacuum seed required by dust recovery)",
        "max_abs_rhs_by_block": maxima,
        "overall_max_abs_rhs": overall,
        "tolerance": tolerance,
        "finite": finite,
        "pass": passed,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "flat_vacuum_geometry_smoke.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(payload, indent=2, sort_keys=True), flush=True)
    if not passed:
        raise SystemExit("FAIL: flat vacuum is not stationary under production geometry RHS blocks")
    print("PASS: exact flat vacuum is stationary under production geometry RHS blocks")


if __name__ == "__main__":
    main()
