#!/usr/bin/env python3
"""Audit kinematic consistency in a pinned radiation replay ZIP; no evolution is run."""
from __future__ import annotations
import argparse
import hashlib
import io
import json
from pathlib import Path
import sys
import zipfile
import numpy as np

EXPECTED_ZIP_SHA256 = "5c3dee4cbe62fb410ed7bcdf97707385019e00bae4f925ab51d4242b42933d09"
ROOT = "artifacts/current-main-radiation-frozen-state-replay"
SCHEMA = "radiation_explicit_predictor_stage_capture_v5_retry_identity_and_recomputable_independent_l3"

def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()

def key(row: dict, ordinal_name: str) -> tuple:
    return (float(row["input_t"]), int(row["attempt_id"]),
            float(row["trial_dt"]), int(row[ordinal_name]))

def maxabs(x: np.ndarray) -> float:
    return float(np.max(np.abs(x))) if x.size else 0.0

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact_zip", type=Path)
    parser.add_argument("--allow-other-zip-hash", action="store_true")
    args = parser.parse_args()
    blob = args.artifact_zip.read_bytes()
    zip_hash = digest(blob)
    if not args.allow_other_zip_hash and zip_hash != EXPECTED_ZIP_SHA256:
        raise ValueError("exact artifact ZIP hash mismatch: " + zip_hash)

    with zipfile.ZipFile(io.BytesIO(blob)) as z:
        names = z.namelist()
        capture_name = ROOT + "/stage_capture.json"
        profiles_name = ROOT + "/stage_capture_profiles.npz"
        manifest_name = ROOT + "/SHA256SUMS.txt"
        for name in (capture_name, profiles_name, manifest_name):
            if name not in names:
                raise ValueError("missing artifact member: " + name)
        manifest_files = []
        for line in z.read(manifest_name).decode().splitlines():
            if not line.strip():
                continue
            expected, filename = line.strip().split(maxsplit=1)
            filename = filename.lstrip("* ")
            matches = [n for n in names if n.rsplit("/", 1)[-1] == filename]
            if len(matches) != 1 or digest(z.read(matches[0])) != expected:
                raise ValueError("internal SHA256SUMS failure: " + filename)
            manifest_files.append(filename)
        capture = json.loads(z.read(capture_name))
        profiles = np.load(io.BytesIO(z.read(profiles_name)), allow_pickle=False)

    if capture.get("schema") != SCHEMA:
        raise ValueError("unexpected stage-capture schema")
    geos = capture["geometry_source_term_capture_records"]
    rates = capture["geometry_stage_explicit_rates"]
    rhss = capture["rhs_component_decomposition"]
    if (len(geos), len(rates), len(rhss)) != (703, 703, 677):
        raise ValueError(f"unexpected record counts: {len(geos)}, {len(rates)}, {len(rhss)}")

    geo_map = {key(x, "ordinal"): x for x in geos}
    rate_map = {key(x, "ordinal"): x for x in rates}
    rhs_map = {key(x, "rhs_ordinal"): x for x in rhss}
    if len(geo_map) != 703 or len(rate_map) != 703 or len(rhs_map) != 677:
        raise ValueError("duplicate stage/RHS identity key")
    if geo_map.keys() != rate_map.keys():
        raise ValueError("geometry source/rate key mismatch")

    id_abs, id_rel, outer_abs = [], [], []
    formula_errors = {"adot": [], "bdot": [], "Xdot": []}
    for k, g in geo_map.items():
        p = g["prefix"]
        if g["stage_state_sha256"] != rate_map[k]["stage_state_sha256"]:
            raise ValueError("geometry source/rate state hash mismatch: " + p)
        def arr(field):
            return np.asarray(profiles[f"{p}__{field}"], dtype=float)
        a, b, X = arr("metric__a"), arr("metric__b"), arr("metric__X")
        Aa, K, alpha, beta = arr("metric__Aa"), arr("metric__K"), arr("metric__alpha"), arr("metric__beta")
        adot, bdot, Xdot = arr("explicit_rate__a"), arr("explicit_rate__b"), arr("explicit_rate__X")
        if np.max(np.abs(beta)) != 0.0:
            raise ValueError("zero-shift assumption failed in stage " + p)
        formula_errors["adot"].append(maxabs(adot + 2.0*alpha*a*Aa))
        formula_errors["bdot"].append(maxabs(bdot - alpha*b*Aa))
        formula_errors["Xdot"].append(maxabs(Xdot - alpha*X*K/3.0))
        lhs = 2.0*Xdot/X - adot/a
        rhs = 2.0*alpha*(K/3.0 + Aa)
        resid = lhs-rhs
        denom = np.abs(lhs)+np.abs(rhs)
        rel = np.divide(np.abs(resid), denom, out=np.zeros_like(resid), where=denom > 0)
        id_abs.append(maxabs(resid))
        id_rel.append(maxabs(rel))
        outer_abs.append(abs(float(resid[-1])))

    cone_abs, cone_rel = [], []
    matching_count = 0
    max_state_diff = 0.0
    outer_rows = []
    for k, r in rhs_map.items():
        if k[3] != 0:
            continue
        if k not in geo_map or k not in rate_map:
            raise ValueError("ordinal-0 RHS has no matching same-time/attempt/dt geometry")
        g = geo_map[k]
        p, rp = g["prefix"], r["prefix"]
        # Confirm exact stage-state profile identity before combining arrays.
        diffs = []
        for field in ("a", "X", "Aa", "K", "alpha"):
            x = np.asarray(profiles[f"{p}__metric__{field}"], dtype=float)
            y = np.asarray(profiles[f"{rp}__stage_geometry__{field}"], dtype=float)
            diffs.append(maxabs(x-y))
        state_diff = max(diffs)
        max_state_diff = max(max_state_diff, state_diff)
        if state_diff != 0.0:
            raise ValueError("ordinal-0 RHS and geometry stage profiles differ")
        a = np.asarray(profiles[f"{p}__metric__a"], dtype=float)
        X = np.asarray(profiles[f"{p}__metric__X"], dtype=float)
        Aa = np.asarray(profiles[f"{p}__metric__Aa"], dtype=float)
        K = np.asarray(profiles[f"{p}__metric__K"], dtype=float)
        alpha = np.asarray(profiles[f"{p}__metric__alpha"], dtype=float)
        adot = np.asarray(profiles[f"{p}__explicit_rate__a"], dtype=float)
        Xdot = np.asarray(profiles[f"{p}__explicit_rate__X"], dtype=float)
        Ur = np.asarray(profiles[f"{rp}__radiation__momentum_r"], dtype=float)
        UE = np.asarray(profiles[f"{rp}__radiation__energy_t"], dtype=float)
        Urdot = np.asarray(profiles[f"{rp}__rhs__radiation__momentum_r"], dtype=float)
        UEdot = np.asarray(profiles[f"{rp}__rhs__radiation__energy_t"], dtype=float)
        qmetric = X/np.sqrt(a)
        direct = -qmetric*np.abs(Ur)*(Xdot/X-0.5*adot/a)
        reduced = -qmetric*np.abs(Ur)*alpha*(K/3.0+Aa)
        err = direct-reduced
        denom = np.abs(direct)+np.abs(reduced)
        rel = np.divide(np.abs(err), denom, out=np.zeros_like(err), where=denom > 0)
        cone_abs.append(maxabs(err))
        cone_rel.append(maxabs(rel))
        matter = UEdot-qmetric*np.sign(Ur)*Urdot
        C = UE-qmetric*np.abs(Ur)
        outer_rows.append({
            "t": k[0], "attempt": k[1], "dt": k[2],
            "C": float(C[-1]), "matter_dotC": float(matter[-1]),
            "metric_dotC_direct": float(direct[-1]),
            "metric_dotC_reduced": float(reduced[-1]),
            "net_dotC": float(matter[-1]+direct[-1]),
            "metric_identity_abs_error": float(err[-1]),
            "alpha_times_K_over_3_plus_Aa": float(alpha[-1]*(K[-1]/3.0+Aa[-1])),
        })
        matching_count += 1

    first = next(x for x in outer_rows if x["t"] == 47.130000000001026 and x["attempt"] == 1 and x["dt"] == 0.03)
    last = next(x for x in outer_rows if x["t"] == 47.15300258518951 and x["dt"] == 6.36411169718682e-15)
    print(json.dumps({
        "artifact_zip_sha256": zip_hash,
        "internal_manifest_verified": sorted(manifest_files),
        "source_baseline": "29d828b4d72b017c12a20b6dc90549f2511f025c",
        "vendor_pin": "d6052d605673ce9d82cdc99b0df79855fa2ea215",
        "counts": {"geometry_stages": len(geos), "geometry_cells": len(geos)*80,
                   "rhs_records": len(rhss), "ordinal0_exact_matches": matching_count,
                   "max_rhs_geometry_profile_difference": max_state_diff},
        "max_abs_explicit_rate_formula_residual": {k: max(v) for k,v in formula_errors.items()},
        "kinematic_identity": {
            "equation": "2*Xdot/X - adot/a = 2*alpha*(K/3 + Aa)",
            "max_abs_residual_all_stages_cells": max(id_abs),
            "max_relative_residual_all_stages_cells": max(id_rel),
            "max_abs_residual_outer_cell": max(outer_abs),
        },
        "cone_metric_term": {
            "max_abs_direct_vs_reduced_difference": max(cone_abs),
            "max_relative_direct_vs_reduced_difference": max(cone_rel),
            "first_rejected_predictor_outer_cell": first,
            "terminal_machine_resolution_attempt_outer_cell": last,
        },
        "admission_passed": False,
        "interpretation": "Kinematic/RHS ledger consistency only; not physical/scientific admission.",
    }, indent=2, sort_keys=True))
    return 0

if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, KeyError, ValueError, zipfile.BadZipFile, json.JSONDecodeError) as exc:
        print("AUDIT FAILED CLOSED: " + str(exc), file=sys.stderr)
        raise SystemExit(2)
