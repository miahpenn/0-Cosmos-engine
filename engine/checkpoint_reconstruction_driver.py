"""Replay the six preregistered N320 Dbase checkpoints; never evolve the state."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

from engine.checkpoint_reconstruction import (
    ARTIFACT_ID, ARTIFACT_SHA256, DECOMP_TOL, L2_FIELDS, SOURCE_COMMIT,
    TARGETS, MATCHED, load_checkpoint_npz,
    match_ledger_row_by_checkpoint_time, reconstruct_l2, compare_field,
)

TARGET_TIME_TOL = 0.001875 / 2.0 + 1.0e-9
EXPECTED_LEDGER_ROWS = 12001


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def checkpoint_label(target: float) -> str:
    # The target is used only to locate the intended checkpoint filename.
    # Ledger pairing always uses the checkpoint's embedded stored timestamp.
    return f"checkpoint_t{target:012.6f}.npz"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifact-root", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--source-commit", default=SOURCE_COMMIT)
    parser.add_argument("--artifact-zip-sha256", default=ARTIFACT_SHA256)
    args = parser.parse_args()

    if args.artifact_zip_sha256 != ARTIFACT_SHA256:
        raise SystemExit(
            f"artifact ZIP SHA-256 mismatch: {args.artifact_zip_sha256}"
        )
    ledger_files = list(args.artifact_root.glob("**/ledgers.json"))
    if len(ledger_files) != 1:
        raise SystemExit(f"expected exactly one ledgers.json, found {len(ledger_files)}")

    ledger_path = ledger_files[0]
    raw_bytes = ledger_path.read_bytes()
    raw = json.loads(raw_bytes)
    rows = raw["history"]
    if len(rows) != EXPECTED_LEDGER_ROWS:
        raise SystemExit(
            f"expected the archived 12001-row history, found {len(rows)}"
        )
    if rows[-1].get("t") != 22.5:
        raise SystemExit(f"expected final ledger row t=22.5, got {rows[-1].get('t')!r}")
    last_interval = float(rows[-1]["t"]) - float(rows[-2]["t"])
    if not (0.0 < last_interval < 1.0e-9):
        raise SystemExit(f"unexpected final ledger interval: {last_interval!r}")

    case_root = ledger_path.parent
    out_results = []
    integrity_pass = True
    total_exact = 0
    total_fields = 0

    for target in TARGETS:
        checkpoint_path = case_root / checkpoint_label(target)
        if not checkpoint_path.is_file():
            raise SystemExit(f"missing intended checkpoint: {checkpoint_path}")
        checkpoint_sha = sha256_file(checkpoint_path)
        state = load_checkpoint_npz(checkpoint_path, r_max=80.0)
        checkpoint_t = float(state.t)
        target_offset = checkpoint_t - target
        if abs(target_offset) > TARGET_TIME_TOL:
            raise SystemExit(
                f"checkpoint {checkpoint_path.name} outside frozen target tolerance: "
                f"stored t={checkpoint_t:.17g}, target={target:.17g}, "
                f"offset={target_offset:.17g}"
            )

        pairing = match_ledger_row_by_checkpoint_time(checkpoint_t, rows)
        if pairing["status"] != MATCHED:
            raise SystemExit(
                f"target {target}: same-step pairing {pairing['status']}; "
                f"candidate timestamps={pairing['candidates']}"
            )
        recorded_row = pairing["row"]
        reconstructed = reconstruct_l2(state)

        field_results = {}
        fields_valid = True
        for field in L2_FIELDS:
            value = recorded_row.get(field)
            result = compare_field(reconstructed[field], value)
            field_results[field] = {
                "reconstructed": reconstructed[field],
                "recorded": None if value is None else float(value),
                **result,
            }
            total_fields += 1
            total_exact += int(result["float64_bitwise_equal"])
            fields_valid = fields_valid and result["outcome"] in ("BITWISE_MATCH", "NUMERIC_DIFFERENCE")

        decomposition_error = reconstructed["hamiltonian_decomposition_error_max"]
        decomposition_ok = decomposition_error <= DECOMP_TOL
        target_integrity_ok = fields_valid and decomposition_ok
        integrity_pass = integrity_pass and target_integrity_ok

        out_results.append({
            "target": target,
            "checkpoint_filename": checkpoint_path.name,
            "checkpoint_sha256": checkpoint_sha,
            "checkpoint_t": checkpoint_t,
            "target_offset": target_offset,
            "target_time_tolerance": TARGET_TIME_TOL,
            "ledger_row_index": int(pairing["row_index"]),
            "ledger_t": float(recorded_row["t"]),
            "ledger_offset": float(pairing["offset"]),
            "ledger_candidate_count": len(pairing["candidates"]),
            "l2_fields": field_results,
            "hamiltonian_decomposition_error_max": decomposition_error,
            "decomposition_gate_tolerance": DECOMP_TOL,
            "decomposition_gate_ok": bool(decomposition_ok),
            "integrity_status": "PASS" if target_integrity_ok else "FAIL",
        })

    payload = {
        "protocol": (
            "Archived N320 Dbase checkpoint reconstruction; exact stored-time "
            "pairing applies only to this dataset; diagnostic replay only"
        ),
        "dataset_scope": "Actions artifact 11638114258 only",
        "source_commit_for_archived_run": SOURCE_COMMIT,
        "audit_code_commit": args.source_commit,
        "artifact_id": ARTIFACT_ID,
        "artifact_zip_sha256": args.artifact_zip_sha256,
        "ledger_path_within_artifact": str(ledger_path.relative_to(args.artifact_root)),
        "ledger_file_sha256": sha256_file(ledger_path),
        "ledger_rows": len(rows),
        "final_two_ledger_times": [float(rows[-2]["t"]), float(rows[-1]["t"])],
        "final_ledger_interval": last_interval,
        "targets": list(TARGETS),
        "l2_fields": list(L2_FIELDS),
        "comparison_rule": (
            "Report exact float64 equality, absolute errors, and relative errors. "
            "Do not assign numerical acceptance PASS/FAIL: the frozen reconstruction "
            "acceptance tolerance was not found in available records. Exact timestamp "
            "pairing is Rev 4 dataset-scoped only."
        ),
        "formal_comparison_acceptance": "UNRESOLVED_NO_FROZEN_NUMERICAL_TOLERANCE_LOCATED",
        "decomposition_tolerance": DECOMP_TOL,
        "exact_float64_matches": total_exact,
        "field_comparisons": total_fields,
        "audit_status": "COMPLETED" if integrity_pass else "INTEGRITY_FAILURE",
        "integrity_checks_pass": bool(integrity_pass),
        "numerical_acceptance_decision": "UNRESOLVED",
        "results": out_results,
        "interpretation_limits": (
            "This checks reproducibility of stored diagnostic values from saved checkpoint states. "
            "The production diagnostic method, state dataclasses, and pinned vendor operators are shared "
            "dependencies. It does not independently validate the residual equations or physical model. "
            "No state evolution or physics campaign is executed."
        ),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )

    for item in out_results:
        exact = sum(v["float64_bitwise_equal"] for v in item["l2_fields"].values())
        print(
            f"target={item['target']:>4} stored_t={item['checkpoint_t']:.15f} "
            f"ledger_row={item['ledger_row_index']} offset={item['ledger_offset']:.1e} "
            f"integrity={'PASS' if item['integrity_status'] == 'PASS' else 'FAIL'} "
            f"bitwise={exact}/4 decomposition={item['hamiltonian_decomposition_error_max']:.3e}"
        )
    print(f"exact_float64_matches={total_exact}/{total_fields}")
    print(f"integrity_checks={'PASS' if integrity_pass else 'FAIL'}")
    print("formal_numerical_acceptance=UNRESOLVED")
    print(f"results_json={args.out}")
    return 0 if integrity_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
