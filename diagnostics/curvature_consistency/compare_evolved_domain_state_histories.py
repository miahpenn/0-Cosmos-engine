"""Compare two archived full-state trajectories on a shared radial grid.

Intended use: R=80 production full-state capture vs the existing R=160 full-
state CMC crop artifact (optionally plus the prior R=160 H_vendor profile).
This is post-processing only; it does not initialize or evolve the engine.
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import numpy as np

GEOMETRY_FIELDS = ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B")
SCALAR_FIELDS = ("S", "PS", "D", "PD", "phi", "Pi")
MATTER_FIELDS = ("rest", "energy_t", "momentum_r")
SPECIES = ("dark_matter", "baryons", "radiation")

def read_jsonl(path: Path):
    rows = [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
    if not rows:
        raise ValueError(f"no JSONL records found: {path}")
    return rows

def field_pairs(small_snapshot, reference_snapshot):
    for key in GEOMETRY_FIELDS:
        yield f"geometry.{key}", small_snapshot["geometry"][key], reference_snapshot["geometry"][key]
    for key in SCALAR_FIELDS:
        yield f"scalars.{key}", small_snapshot["scalars"][key], reference_snapshot["scalars"][key]
    for species in SPECIES:
        for key in MATTER_FIELDS:
            yield f"matter.{species}.{key}", small_snapshot["matter"][species][key], reference_snapshot["matter"][species][key]

def metrics(small, reference, radius):
    s = np.asarray(small, dtype=float)
    r = np.asarray(reference, dtype=float)
    if s.shape != r.shape or s.shape != radius.shape:
        raise ValueError(f"shared-grid shape mismatch: small={s.shape}, reference={r.shape}, radius={radius.shape}")
    delta = s - r
    abs_delta = np.abs(delta)
    i = int(np.argmax(abs_delta))
    core = radius <= min(20.0, float(radius[-1]))
    core_indices = np.where(core)[0]
    ic = int(core_indices[int(np.argmax(abs_delta[core]))])
    norm_r = float(np.linalg.norm(r))
    norm_d = float(np.linalg.norm(delta))
    relative_l2 = norm_d / norm_r if norm_r > 0.0 else (0.0 if norm_d == 0.0 else None)
    scale = max(float(np.max(np.abs(s))), float(np.max(np.abs(r))))
    numeric_threshold = max(1.0e-12, 1.0e-10 * scale)
    return {
        "max_abs": float(abs_delta[i]),
        "max_abs_radius": float(radius[i]),
        "signed_delta_at_max": float(delta[i]),
        "rms_abs": float(np.sqrt(np.mean(delta * delta))),
        "relative_l2_to_reference": relative_l2,
        "core_max_abs_r_le_20": float(abs_delta[core].max()),
        "core_max_radius_r_le_20": float(radius[ic]),
        "core_rms_abs_r_le_20": float(np.sqrt(np.mean(delta[core] * delta[core]))),
        "numeric_threshold_for_onset": numeric_threshold,
        "above_numeric_threshold": bool(abs_delta[i] > numeric_threshold),
    }

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--small", required=True, type=Path, help="R=80 full-state JSONL from the new capture")
    parser.add_argument("--reference", required=True, type=Path, help="R=160 full-state JSONL from the CMC crop artifact")
    parser.add_argument("--reference-h", type=Path, help="Optional R=160 legacy profile JSONL containing profile.H_vendor")
    parser.add_argument("--out", required=True, type=Path, help="Output directory")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    small_rows = read_jsonl(args.small)
    reference_rows = read_jsonl(args.reference)
    if len(small_rows) != len(reference_rows):
        raise ValueError(f"sample count mismatch: small={len(small_rows)}, reference={len(reference_rows)}")
    h_rows = read_jsonl(args.reference_h) if args.reference_h else None
    if h_rows is not None and len(h_rows) != len(small_rows):
        raise ValueError(f"H-profile sample count mismatch: small={len(small_rows)}, reference_h={len(h_rows)}")

    records = []
    for idx, (small_row, reference_row) in enumerate(zip(small_rows, reference_rows)):
        ts = float(small_row["t"])
        tr = float(reference_row["t"])
        if abs(ts - tr) > 1.0e-8:
            raise ValueError(f"sample times do not align at index {idx}: small={ts}, reference={tr}")
        radius = np.asarray(small_row.get("radius", small_row.get("profile", {}).get("radius")), dtype=float)
        ref_radius = np.asarray(reference_row.get("radius_full", reference_row.get("radius", reference_row.get("profile", {}).get("radius"))), dtype=float)
        if radius.ndim != 1 or ref_radius.ndim != 1 or len(radius) > len(ref_radius):
            raise ValueError(f"invalid radius arrays at sample {idx}")
        if not np.allclose(radius, ref_radius[:len(radius)], rtol=0.0, atol=1.0e-12):
            raise ValueError(f"radial grids do not share identical cell centers at sample {idx}")
        s_snapshot = small_row["state_snapshot"]
        r_snapshot = reference_row["state_snapshot"]
        sample_record = {"index": idx, "t": ts, "tau_small": small_row.get("tau"), "tau_reference": reference_row.get("tau"), "fields": {}}
        for field_name, s_values, r_values in field_pairs(s_snapshot, r_snapshot):
            m = metrics(s_values, np.asarray(r_values, dtype=float)[:len(radius)], radius)
            sample_record["fields"][field_name] = m
            records.append({"index": idx, "t": ts, "field": field_name, **m})
        # Compare H_vendor when the matching large-domain diagnostic profile is supplied.
        if h_rows is not None:
            hrow = h_rows[idx]
            if abs(float(hrow["t"]) - ts) > 1.0e-8:
                raise ValueError(f"H-profile times do not align at sample {idx}")
            hs = np.asarray(small_row["H_vendor"], dtype=float)
            hr = np.asarray(hrow["profile"]["H_vendor"], dtype=float)[:len(radius)]
            hm = metrics(hs, hr, radius)
            sample_record["fields"]["diagnostic.H_vendor"] = hm
            records.append({"index": idx, "t": ts, "field": "diagnostic.H_vendor", **hm})
        # Alpha profile duplicate is intentionally checked against the state snapshot.
        alpha_s = np.asarray(s_snapshot["geometry"]["alpha"], dtype=float)
        alpha_r = np.asarray(r_snapshot["geometry"]["alpha"], dtype=float)[:len(radius)]
        if not np.array_equal(alpha_s, np.asarray(small_row["alpha_profile"], dtype=float)):
            raise ValueError(f"R=80 alpha profile and state snapshot disagree at sample {idx}")
        sample_record["alpha_duplicate_check"] = metrics(alpha_s, alpha_r, radius)
    
    # Field-level extrema identify where the independently evolved states differ most.
    field_names = sorted({row["field"] for row in records})
    summaries = []
    for name in field_names:
        field_rows = [row for row in records if row["field"] == name]
        bymax = max(field_rows, key=lambda row: row["max_abs"])
        byrms = max(field_rows, key=lambda row: row["rms_abs"])
        t0 = min(field_rows, key=lambda row: abs(row["t"]))
        t4 = min(field_rows, key=lambda row: abs(row["t"] - 4.0))
        t8 = min(field_rows, key=lambda row: abs(row["t"] - 8.0))
        t12 = min(field_rows, key=lambda row: abs(row["t"] - 12.0))
        onset = next((row for row in field_rows if row["above_numeric_threshold"]), None)
        summaries.append({
            "field": name,
            "max_abs_over_time": bymax["max_abs"],
            "time_of_max_abs": bymax["t"],
            "radius_of_max_abs": bymax["max_abs_radius"],
            "signed_delta_at_max_abs": bymax["signed_delta_at_max"],
            "max_rms_over_time": byrms["rms_abs"],
            "time_of_max_rms": byrms["t"],
            "t0_max_abs": t0["max_abs"],
            "t4_max_abs": t4["max_abs"],
            "t8_max_abs": t8["max_abs"],
            "t12_max_abs": t12["max_abs"],
            "t12_rms_abs": t12["rms_abs"],
            "t12_relative_l2_to_reference": t12["relative_l2_to_reference"],
            "t12_core_max_abs_r_le_20": t12["core_max_abs_r_le_20"],
            "t12_core_rms_abs_r_le_20": t12["core_rms_abs_r_le_20"],
            "first_time_above_numeric_threshold": onset["t"] if onset else None,
            "numeric_threshold_at_onset": onset["numeric_threshold_for_onset"] if onset else None,
        })
    summaries.sort(key=lambda row: row["max_abs_over_time"], reverse=True)

    with (args.out / "field_history.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0].keys()))
        writer.writeheader()
        writer.writerows(records)
    report = {
        "schema": "evolved_state_domain_pair_comparison_v1",
        "status": "completed",
        "diagnostic_only": True,
        "small_input": str(args.small),
        "reference_input": str(args.reference),
        "reference_h_input": str(args.reference_h) if args.reference_h else None,
        "sample_count": len(small_rows),
        "time_start": float(small_rows[0]["t"]),
        "time_end": float(small_rows[-1]["t"]),
        "small_domain": small_rows[0].get("domain"),
        "reference_domain": {"N": reference_rows[0].get("full_N"), "r_max": reference_rows[0].get("full_r_max")},
        "onset_criterion": "first saved sample where max absolute field difference exceeds max(1e-12, 1e-10 * max field magnitude on the shared slice); diagnostic reporting threshold only, not a physics-admission gate",
        "field_summaries_descending_by_max_abs": summaries,
    }
    (args.out / "comparison_report.json").write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n")
    md = [
        "# Evolved-state domain comparison",
        "",
        f"- Samples: {len(small_rows)}; t={small_rows[0]['t']:.6g} to {small_rows[-1]['t']:.6g}",
        "- Comparison: R=80 production state vs R=160 production state restricted to the same radial cell centers.",
        "- Threshold above is only to label the first clearly non-roundoff saved slice; it is not a physics pass/fail criterion.",
        "",
        "| Field | First above numeric threshold | Max difference (time, radius) | RMS at t=12 | Relative L2 at t=12 | Core max at t=12 |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for s in summaries:
        rel = "n/a" if s["t12_relative_l2_to_reference"] is None else f"{s['t12_relative_l2_to_reference']:.4g}"
        onset = "none" if s["first_time_above_numeric_threshold"] is None else f"{s['first_time_above_numeric_threshold']:.4f}"
        md.append(f"| {s['field']} | {onset} | {s['max_abs_over_time']:.4g} ({s['time_of_max_abs']:.4f}, {s['radius_of_max_abs']:.2f}) | {s['t12_rms_abs']:.4g} | {rel} | {s['t12_core_max_abs_r_le_20']:.4g} |")
    md += ["", "## First interpretation rule", "", "Read each field separately. Initial-slice differences are initialization/domain effects; growth from t=0 shows later evolution divergence. Do not infer a physical cause solely from correlations or from the largest field difference.", ""]
    (args.out / "comparison_report.md").write_text("\n".join(md) + "\n")
    print(json.dumps({"status": report["status"], "sample_count": len(small_rows), "fields_compared": len(field_names), "output": str(args.out), "largest_final_time_fields": summaries[:5]}, indent=2))

if __name__ == "__main__":
    main()
