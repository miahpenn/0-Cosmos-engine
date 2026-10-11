#!/usr/bin/env python3
"""Threshold-free spatial/time audit of an existing repair ON/OFF state-array JSON.

Usage: python state_array_spatial_evolution.py state-array-ab.json [output.json]
No simulation is run. No physical or state-difference threshold is applied.
"""
import json
import math
import sys
from pathlib import Path

REGIONS = [
    ("inner_r_le_20", 0.0, 20.0),
    ("middle_20_lt_r_lt_40", 20.0, 40.0),
    ("middle_40_le_r_lt_64", 40.0, 64.0),
    ("outer_r_ge_64", 64.0, math.inf),
]
FIELDS = [
    "geometry.Lambda", "geometry.Aa", "geometry.a", "geometry.b",
    "geometry.alpha", "scalars.D", "matter.dark_matter.energy_t",
]


def summarize(source):
    data = json.loads(Path(source).read_text())
    on_samples, off_samples = data["repair_on"]["samples"], data["repair_off"]["samples"]
    if len(on_samples) != len(off_samples):
        raise ValueError("ON/OFF sample count mismatch")
    output = {
        "schema": "0star-state-array-spatial-evolution-v1",
        "diagnostic_only": True,
        "source": {
            "input_filename": Path(source).name,
            "source_artifact": "state-array-ab-comparison; workflow run 37859756564",
            "source_code_commit": "7b9bac434f402df597641cfe7b24cf1a126a02e1",
            "method": "absolute OFF-minus-ON difference; regional max and unweighted RMS per recorded array",
        },
        "configuration": data["configuration"],
        "thresholds": {"state_difference_threshold_applied": False, "physical_significance_threshold_applied": False},
        "regions": [{"name": name, "r_min_inclusive": lo, "r_max_exclusive": hi if math.isfinite(hi) else None} for name, lo, hi in REGIONS],
        "samples": [],
        "interpretation_limits": [
            "Nonzero differences and their amplitudes are reported without a threshold for physical significance.",
            "Regional growth in absolute differences does not by itself establish propagation, causal influence, or a propagation speed.",
            "The ON/OFF states already differ at t=0; the comparison is not an identical-initial-state experiment.",
            "This is post-processing of an existing artifact, not a new simulation or a physics modification.",
        ],
    }
    for index, (on, off) in enumerate(zip(on_samples, off_samples)):
        if abs(float(on["actual_t"]) - float(off["actual_t"])) > 1e-12:
            raise ValueError(f"physical sample times differ at sample {index}")
        centers = [float(x) for x in on["arrays"]["grid.centers"]]
        if centers != [float(x) for x in off["arrays"]["grid.centers"]]:
            raise ValueError(f"grid centers differ at sample {index}")
        sample = {
            "sample_index": index,
            "actual_t": float(on["actual_t"]),
            "tau_on": on["tau"],
            "tau_off": off["tau"],
            "tau_off_minus_on": float(off["tau"]) - float(on["tau"]),
            "fields": {},
        }
        for field in FIELDS:
            a, b = on["arrays"][field], off["arrays"][field]
            if len(a) != len(centers) or len(b) != len(centers):
                raise ValueError(f"array length mismatch at sample {index}, field {field}")
            diff = [abs(float(y) - float(x)) for x, y in zip(a, b)]
            regions = {}
            for name, lo, hi in REGIONS:
                vals = [diff[j] for j, r in enumerate(centers) if r >= lo and r < hi]
                if not vals:
                    raise ValueError(f"empty region {name}")
                regions[name] = {
                    "max_abs_difference": max(vals),
                    "rms_abs_difference": math.sqrt(sum(v*v for v in vals) / len(vals)),
                    "cell_count": len(vals),
                }
            sample["fields"][field] = regions
        output["samples"].append(sample)
    return output


if __name__ == "__main__":
    if len(sys.argv) not in (2, 3):
        raise SystemExit(__doc__)
    result = summarize(sys.argv[1])
    text = json.dumps(result, indent=2, allow_nan=False) + "\n"
    if len(sys.argv) == 3:
        Path(sys.argv[2]).write_text(text)
    else:
        print(text, end="")
