#!/usr/bin/env python3
"""Frozen-state finite-dt comparison of baseline vs direct matter-RHS update.

This diagnostic evolves one production finite step from independently initialized
states at each resolution. It compares arithmetic implementation output only.
"""
from __future__ import annotations
import argparse, hashlib, json, math, subprocess
from pathlib import Path
import numpy as np

BASELINE_SHA = "3a87679faf676414d81cd30f0c17e0ae33358a56"
DIRECT_RHS_SHA = "a4f7319384a8c060de2ceeb1981ee312c9b98bd9"
VENDOR_SHA = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
RESOLUTIONS = (40, 60, 80, 120, 160, 240, 320)
R_MAX, CFL, BETA_DM = 80.0, 0.03, -0.04
COMMON_SOURCE_PATHS = (
    "engine/production_kernel.py", "engine/v55_initial.py", "engine/v55_matter.py",
    "engine/scalar_system.py", "engine/cmc_gauge.py", "engine/valencia.py",
    "engine/stress_energy.py",
)
CHANGED_PATHS = ("engine/matter_system.py", "engine/matter_rhs.py")
MATTER_CASES = (("dark_matter", "DARK_MATTER"), ("baryons", "BARYON"), ("radiation", "RADIATION"))
STATE_FIELDS = ("rest", "energy_t", "momentum_r")


def _git(*args):
    return subprocess.check_output(["git", *args], text=True, stderr=subprocess.STDOUT).strip()


def _float_list(values):
    arr = np.asarray(values, dtype=float).reshape(-1)
    return [float(x) if math.isfinite(float(x)) else None for x in arr]


def _array_digest_update(h, label, values):
    arr = np.ascontiguousarray(np.asarray(values, dtype=np.float64))
    h.update(label.encode("utf-8")); h.update(b"\0")
    h.update(str(arr.shape).encode("ascii")); h.update(b"\0")
    h.update(arr.tobytes(order="C"))


def _state_digest(state):
    h = hashlib.sha256()
    for name in ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B"):
        _array_digest_update(h, f"geometry.{name}", getattr(state.geometry, name))
    for name in ("S", "PS", "D", "PD", "phi", "Pi"):
        _array_digest_update(h, f"scalars.{name}", getattr(state.scalars, name))
    for species_name, _ in MATTER_CASES:
        species = getattr(state.matter, species_name)
        for name in STATE_FIELDS:
            _array_digest_update(h, f"matter.{species_name}.{name}", getattr(species, name))
    return h.hexdigest()


def _metadata():
    paths = (*COMMON_SOURCE_PATHS, *CHANGED_PATHS)
    blobs = {}
    for path in paths:
        try: blobs[path] = _git("rev-parse", f"HEAD:{path}")
        except subprocess.CalledProcessError: blobs[path] = None
    try: vendor = _git("rev-parse", "HEAD:vendor/bb-palatini-unified-r0")
    except subprocess.CalledProcessError: vendor = None
    return {"commit": _git("rev-parse", "HEAD"), "vendor_commit": vendor, "source_blobs": blobs}


def _radiation_cone(metric_slice, radiation):
    from engine.matter_system import _metric_arrays
    metrics = _metric_arrays(metric_slice)
    sg = np.asarray([m.sqrt_gamma for m in metrics], dtype=float)
    gir = np.asarray([m.gamma_rr_inv for m in metrics], dtype=float)
    energy = np.asarray(radiation.energy_t, dtype=float) / sg
    momentum = np.asarray(radiation.momentum_r, dtype=float) / sg
    abs_s = np.sqrt(np.maximum(gir * momentum**2, 0.0))
    margin = energy - abs_s
    ratio = abs_s / np.maximum(np.abs(energy), 1e-300)
    finite = bool(np.all(np.isfinite(energy)) and np.all(np.isfinite(abs_s))
                  and np.all(np.isfinite(margin)) and np.all(np.isfinite(ratio)))
    imin, imax = int(np.argmin(margin)), int(np.argmax(ratio))
    bad = np.flatnonzero((energy < 0.0) | (abs_s > energy))
    return {
        "finite": finite, "admissible_all_cells": bool(finite and bad.size == 0),
        "bad_cell_count": int(bad.size), "first_bad_cell": int(bad[0]) if bad.size else None,
        "min_margin": float(margin[imin]) if math.isfinite(float(margin[imin])) else None,
        "min_margin_cell": imin,
        "max_ratio": float(ratio[imax]) if math.isfinite(float(ratio[imax])) else None,
        "max_ratio_cell": imax,
    }


def _rates_context(kernel, state):
    from engine.scalar_system import scalar_rhs_arrays
    from engine.v55_matter import dm_density, metric_slice_from_q
    from engine.matter_system import geometry_metric_derivatives
    grid = state.grid
    metric = metric_slice_from_q(grid, state.geometry)
    rho_dm = dm_density(metric, state.matter)
    scalar_rhs = scalar_rhs_arrays(grid, state.geometry, state.scalars, beta_dm=BETA_DM, rho_dm=rho_dm)
    derivatives = geometry_metric_derivatives(metric, lambda values, parity: kernel._d1(grid, values, parity))
    return metric, derivatives, scalar_rhs.phi, kernel._d1(grid, state.scalars.phi, 1)


def capture(out_path):
    from engine.production_kernel import V55ProductionKernel
    from engine.matter_system import Species, evolve_species
    kernel = V55ProductionKernel()
    rows = []
    for n in RESOLUTIONS:
        state = kernel.initialize(resolution=n, r_max=R_MAX)
        metric, derivatives, dphi_t, dphi_r = _rates_context(kernel, state)
        row = {
            "N": n, "r_max": R_MAX, "dr": float(state.grid.dr),
            "dt": CFL * float(state.grid.dr), "t_initial": float(state.t),
            "initial_state_sha256": _state_digest(state),
            "initial_radiation_cone": _radiation_cone(metric, state.matter.radiation),
            "species": {},
        }
        for species_name, enum_name in MATTER_CASES:
            enum = getattr(Species, enum_name)
            before = getattr(state.matter, species_name)
            try:
                after = evolve_species(
                    metric, derivatives, before, enum, row["dt"],
                    dphi_t=dphi_t if species_name == "dark_matter" else None,
                    dphi_r=dphi_r if species_name == "dark_matter" else None,
                    beta_dm=BETA_DM,
                )
                case = {
                    "status": "completed",
                    "fields": {name: _float_list(getattr(after, name)) for name in STATE_FIELDS},
                    "nonfinite_counts": {name: int(np.count_nonzero(~np.isfinite(np.asarray(getattr(after, name), dtype=float)))) for name in STATE_FIELDS},
                }
                if species_name == "radiation":
                    case["radiation_cone_after"] = _radiation_cone(metric, after)
            except Exception as exc:
                case = {"status": "exception", "exception_type": type(exc).__name__, "exception_message": str(exc)}
            row["species"][species_name] = case
        rows.append(row)
    return {
        "schema": "matter_rhs_finite_dt_update_order_capture_v1",
        "kind": "single_real_finite_dt_update_from_default_initialized_state",
        "metadata": _metadata(),
        "expected_source": {
            "baseline_commit": BASELINE_SHA, "direct_rhs_commit": DIRECT_RHS_SHA,
            "vendor_commit": VENDOR_SHA, "resolutions": list(RESOLUTIONS),
            "r_max": R_MAX, "dt_rule": "dt=0.03*dr",
        },
        "rows": rows,
        "interpretation_limit": "One actual finite-dt matter update from default initialization. This is arithmetic/contract evidence only, not physical admission, convergence proof, or long-horizon validation.",
    }


def _json_array(values):
    return np.asarray([np.nan if x is None else float(x) for x in values], dtype=np.float64)


def _ordered_float_ints(arr):
    bits = np.asarray(arr, dtype=np.float64).view(np.uint64)
    sign = np.uint64(1) << np.uint64(63)
    return np.where((bits & sign) == 0, bits | sign, ~bits)


def _field_delta(base_values, direct_values):
    base, direct = _json_array(base_values), _json_array(direct_values)
    if base.shape != direct.shape:
        return {"shape_match": False, "baseline_shape": list(base.shape), "direct_shape": list(direct.shape)}
    fb, fd = np.isfinite(base), np.isfinite(direct)
    both = fb & fd
    result = {
        "shape_match": True, "length": int(base.size),
        "all_finite_in_both": bool(np.all(fb) and np.all(fd)),
        "baseline_nonfinite_count": int(np.count_nonzero(~fb)),
        "direct_nonfinite_count": int(np.count_nonzero(~fd)),
        "compared_finite_cell_count": int(np.count_nonzero(both)),
    }
    if not np.any(both):
        return {**result, "exactly_equal_on_finite_cells": False,
                "differing_finite_cell_count": None, "max_abs_difference": None,
                "max_relative_difference_symmetric": None, "max_ulp_difference": None}
    b, d = base[both], direct[both]
    delta = np.abs(d - b)
    scale = np.maximum(np.maximum(np.abs(b), np.abs(d)), 1e-300)
    ob, od = _ordered_float_ints(b).tolist(), _ordered_float_ints(d).tolist()
    ulps = [abs(int(x) - int(y)) for x, y in zip(ob, od)]
    result.update({
        "exactly_equal_on_finite_cells": bool(np.all(d == b)),
        "differing_finite_cell_count": int(np.count_nonzero(d != b)),
        "max_abs_difference": float(np.max(delta)),
        "max_relative_difference_symmetric": float(np.max(delta / scale)),
        "max_ulp_difference": max(ulps, default=0),
    })
    return result


def compare(base_path, direct_path, out_path, md_path):
    base, direct = json.loads(base_path.read_text()), json.loads(direct_path.read_text())
    problems = []
    if base.get("metadata", {}).get("commit") != BASELINE_SHA:
        problems.append("baseline commit does not match pinned main")
    if direct.get("metadata", {}).get("commit") != DIRECT_RHS_SHA:
        problems.append("direct-RHS commit does not match pinned PR #8 head")
    for label, doc in (("baseline", base), ("direct_rhs", direct)):
        if doc.get("metadata", {}).get("vendor_commit") != VENDOR_SHA:
            problems.append(f"{label} vendor gitlink differs from required pin")
        for path in CHANGED_PATHS:
            if not doc.get("metadata", {}).get("source_blobs", {}).get(path):
                problems.append(f"{label} capture missing source blob for {path}")
    for path in COMMON_SOURCE_PATHS:
        if base.get("metadata", {}).get("source_blobs", {}).get(path) != direct.get("metadata", {}).get("source_blobs", {}).get(path):
            problems.append(f"common production source drift: {path}")
    rb_map = {int(x["N"]): x for x in base.get("rows", [])}
    rd_map = {int(x["N"]): x for x in direct.get("rows", [])}
    initial_identity, cases = {}, []
    for n in RESOLUTIONS:
        rb, rd = rb_map.get(n), rd_map.get(n)
        if rb is None or rd is None:
            problems.append(f"missing N={n} in a capture")
            continue
        same_state = rb.get("initial_state_sha256") == rd.get("initial_state_sha256")
        same_dt = rb.get("dt") == rd.get("dt")
        initial_identity[str(n)] = {
            "same_state_hash": same_state, "baseline_hash": rb.get("initial_state_sha256"),
            "direct_hash": rd.get("initial_state_sha256"), "same_dt": same_dt,
            "baseline_dt": rb.get("dt"), "direct_dt": rd.get("dt"),
        }
        if not same_state: problems.append(f"initial state mismatch at N={n}")
        if not same_dt: problems.append(f"dt mismatch at N={n}")
        for species_name, _ in MATTER_CASES:
            sb, sd = rb.get("species", {}).get(species_name, {}), rd.get("species", {}).get(species_name, {})
            item = {
                "N": n, "dt": rb.get("dt"), "species": species_name,
                "baseline_status": sb.get("status"), "direct_rhs_status": sd.get("status"),
                "baseline_exception": None if sb.get("status") == "completed" else {"type": sb.get("exception_type"), "message": sb.get("exception_message")},
                "direct_rhs_exception": None if sd.get("status") == "completed" else {"type": sd.get("exception_type"), "message": sd.get("exception_message")},
                "matched_input_state": bool(same_state and same_dt),
            }
            if sb.get("status") == sd.get("status") == "completed" and same_state and same_dt:
                item["fields"] = {field: _field_delta(sb["fields"][field], sd["fields"][field]) for field in STATE_FIELDS}
                if species_name == "radiation":
                    cb, ab, ad = rb.get("initial_radiation_cone", {}), sb.get("radiation_cone_after", {}), sd.get("radiation_cone_after", {})
                    item["initial_cone_admissible"] = cb.get("admissible_all_cells")
                    item["baseline_cone_after"], item["direct_rhs_cone_after"] = ab, ad
                    item["radiation_cone_delta"] = {
                        key: {"baseline": ab.get(key), "direct_rhs": ad.get(key),
                              "direct_minus_baseline": None if ab.get(key) is None or ad.get(key) is None else float(ad[key] - ab[key])}
                        for key in ("min_margin", "max_ratio")
                    }
            else:
                item["fields"] = None
            cases.append(item)
    completed_pairs = sum(1 for x in cases if x.get("matched_input_state") and x.get("baseline_status") == x.get("direct_rhs_status") == "completed")
    field_deltas = [d for row in cases for d in (row.get("fields") or {}).values() if d.get("shape_match")]
    report = {
        "schema": "matter_rhs_finite_dt_update_order_comparison_v1",
        "kind": "baseline_componentwise_update_vs_direct_rhs_update_on_matched_states",
        "baseline_capture": str(base_path), "direct_rhs_capture": str(direct_path),
        "pins": {"baseline": BASELINE_SHA, "direct_rhs": DIRECT_RHS_SHA, "vendor": VENDOR_SHA},
        "common_source_blobs": {p: base.get("metadata", {}).get("source_blobs", {}).get(p) for p in COMMON_SOURCE_PATHS},
        "changed_source_blobs": {label: {p: doc.get("metadata", {}).get("source_blobs", {}).get(p) for p in CHANGED_PATHS} for label, doc in (("baseline", base), ("direct_rhs", direct))},
        "initial_identity_by_resolution": initial_identity,
        "quality_gate": {"comparison_pins_and_inputs_valid": not problems, "problems": problems,
                         "matched_species_update_pairs_completed": completed_pairs,
                         "expected_species_update_pairs": len(RESOLUTIONS) * len(MATTER_CASES)},
        "cases": cases,
        "overall_max_abs_difference": max((x.get("max_abs_difference") or 0.0 for x in field_deltas), default=0.0),
        "overall_max_relative_difference_symmetric": max((x.get("max_relative_difference_symmetric") or 0.0 for x in field_deltas), default=0.0),
        "overall_max_ulp_difference": max((x.get("max_ulp_difference") or 0 for x in field_deltas), default=0),
        "interpretation_limit": "Matched one-step numerical comparison only. It does not select an implementation by itself, validate radiation source physics, close centre/constraint or Misner-Sharp gates, or admit a long-horizon campaign.",
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    md_path.parent.mkdir(parents=True, exist_ok=True)
    md_path.write_text(_markdown(report))
    return report


def _markdown(report):
    q = report["quality_gate"]
    lines = [
        "# Matter RHS update-order A/B", "",
        "One real finite-dt update from the default production initializer. Arithmetic/contract evidence only; not scientific admission.", "",
        f"- Baseline commit: {report['pins']['baseline']}",
        f"- Direct-RHS commit: {report['pins']['direct_rhs']}",
        f"- Vendor commit: {report['pins']['vendor']}",
        f"- Input/pin gate: **{'PASS' if q['comparison_pins_and_inputs_valid'] else 'FAIL'}**",
        f"- Matched successful species updates: {q['matched_species_update_pairs_completed']} / {q['expected_species_update_pairs']}",
        f"- Overall max absolute delta: {report['overall_max_abs_difference']:.9e}",
        f"- Overall max symmetric relative delta: {report['overall_max_relative_difference_symmetric']:.9e}",
        f"- Overall max ULP delta: {report['overall_max_ulp_difference']}", "",
        "| N | dt | Species | Field | Differing cells | Max abs delta | Max rel delta | Max ULP |",
        "|---:|---:|---|---|---:|---:|---:|---:|",
    ]
    for row in report["cases"]:
        for field, d in (row.get("fields") or {}).items():
            lines.append(f"| {row['N']} | {row['dt']:.9e} | {row['species']} | {field} | {d.get('differing_finite_cell_count', 'n/a')} | {d.get('max_abs_difference', 'n/a')} | {d.get('max_relative_difference_symmetric', 'n/a')} | {d.get('max_ulp_difference', 'n/a')} |")
        if row["baseline_status"] != "completed" or row["direct_rhs_status"] != "completed":
            lines.append(f"| {row['N']} | {row['dt']:.9e} | {row['species']} | status | n/a | baseline={row['baseline_status']} / direct={row['direct_rhs_status']} | n/a | n/a |")
        if row["species"] == "radiation" and row.get("radiation_cone_delta"):
            m, r = row["radiation_cone_delta"]["min_margin"], row["radiation_cone_delta"]["max_ratio"]
            lines.append(f"| {row['N']} | {row['dt']:.9e} | radiation | cone min-margin delta / max-ratio delta | n/a | {m['direct_minus_baseline']} | {r['direct_minus_baseline']} | n/a |")
    lines.extend(["", "## Gate problems", ""])
    lines.extend(f"- {p}" for p in q["problems"]) if q["problems"] else lines.append("None.")
    lines.extend(["", "## Interpretation limit", "", report["interpretation_limit"], ""])
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="mode", required=True)
    p = sub.add_parser("capture")
    p.add_argument("--out", type=Path, required=True)
    p = sub.add_parser("compare")
    p.add_argument("--baseline", type=Path, required=True)
    p.add_argument("--direct", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--markdown", type=Path, required=True)
    args = parser.parse_args()
    if args.mode == "capture":
        report = capture(args.out)
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
        print(json.dumps({"status": "capture_complete", "commit": report["metadata"]["commit"], "vendor_commit": report["metadata"]["vendor_commit"],
                          "resolutions": [x["N"] for x in report["rows"]],
                          "species_statuses": {str(x["N"]): {k: v["status"] for k, v in x["species"].items()} for x in report["rows"]},
                          "output": str(args.out)}, allow_nan=False))
    else:
        report = compare(args.baseline, args.direct, args.out, args.markdown)
        print(json.dumps({"status": "comparison_complete", "quality_gate": report["quality_gate"],
                          "overall_max_abs_difference": report["overall_max_abs_difference"],
                          "overall_max_relative_difference_symmetric": report["overall_max_relative_difference_symmetric"],
                          "overall_max_ulp_difference": report["overall_max_ulp_difference"],
                          "output": str(args.out), "markdown": str(args.markdown)}, allow_nan=False))
        if not report["quality_gate"]["comparison_pins_and_inputs_valid"]:
            raise SystemExit(2)


if __name__ == "__main__":
    main()
