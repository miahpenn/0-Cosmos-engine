#!/usr/bin/env python3
"""Diagnostic-only initial-state trace; compares candidate and pinned historical engines.

This is a newly authored runner, not the hash-mismatched reconstructed Slack file.
It never changes expected digests or production code. A digest mismatch is reported
and returns nonzero only after the JSON evidence has been written.
"""
from __future__ import annotations
import argparse, hashlib, json, os, platform, subprocess, sys, traceback
from pathlib import Path
import numpy as np

EXPECTED_INITIAL = "a256b57f274ee120b55c0f805ed488d91382a9799647c8f9811d7df82aa337c0"
EXPECTED_ACCEPTED = "c6acb3711d59a85c804a88208c9bc0356e8227a806c0684ae3a8afac72c9091a"
FIELD_SPECS = [
    ("grid.centers", "grid", "centers"),
    *[(f"geometry.{n}", "geometry", n) for n in ("a","b","X","Aa","K","Lambda","alpha","beta","B")],
    *[(f"scalars.{n}", "scalars", n) for n in ("S","PS","D","PD","phi","Pi")],
    *[(f"matter.{s}.{n}", "matter."+s, n) for s in ("dark_matter","baryons","radiation") for n in ("rest","energy_t","momentum_r")],
]

def sha256_array(value):
    original = np.asarray(value)
    arr = np.ascontiguousarray(np.asarray(value, dtype=np.float64))
    return {"shape": list(original.shape), "source_dtype": str(original.dtype),
            "source_byteorder": original.dtype.byteorder,
            "source_c_contiguous": bool(original.flags.c_contiguous),
            "canonical_dtype": str(arr.dtype), "canonical_byteorder": arr.dtype.byteorder,
            "canonical_sha256": hashlib.sha256(arr.tobytes(order="C")).hexdigest(),
            "source_raw_sha256": hashlib.sha256(np.ascontiguousarray(original).tobytes(order="C")).hexdigest(),
            "min": float(np.min(arr)) if arr.size else None,
            "max": float(np.max(arr)) if arr.size else None,
            "max_abs": float(np.max(np.abs(arr))) if arr.size else None,
            "sum": float(np.sum(arr, dtype=np.float64)) if arr.size else 0.0,
            "first_5": arr.reshape(-1)[:5].tolist(),
            "last_5": arr.reshape(-1)[-5:].tolist()}

def resolve_field(state, path):
    obj = state
    for part in path.split("."):
        obj = getattr(obj, part)
    return obj

def collect(engine_root, output_path, role):
    root = Path(engine_root).resolve()
    sys.path.insert(0, str(root))
    from engine.production_kernel import V55ProductionKernel
    state = V55ProductionKernel().initialize(resolution=60, r_max=80.0, include_radiation=True)
    arrays = {name: np.ascontiguousarray(np.asarray(resolve_field(state, ".".join(rest)), dtype=np.float64))
              for name, *rest in FIELD_SPECS}
    np.savez_compressed(output_path.with_suffix(".npz"), **arrays)
    fields = {name: sha256_array(arrays[name]) for name, *_ in FIELD_SPECS}
    h = hashlib.sha256()
    for name, *_ in FIELD_SPECS:
        spec = fields[name]
        h.update(name.encode("utf-8")); h.update(str(tuple(spec["shape"])).encode("ascii"))
        # Rebuild the historical state digest's canonical bytes by resolving the
        # explicit object path below (the named field key is already stable).
        path = next(".".join(rest) for field_name, *rest in FIELD_SPECS if field_name == name)
        arr = np.ascontiguousarray(np.asarray(resolve_field(state, path), dtype=np.float64))
        h.update(arr.tobytes(order="C"))
    trace = {
        "role": role, "engine_root": str(root),
        "module_paths": {name: str(Path(sys.modules[name].__file__).resolve()) for name in
                         ("engine.production_kernel", "engine.v55_initial", "engine.v55_pirk_adapter")},
        "initial_digest": h.hexdigest(),
        "expected_historical_initial_digest": EXPECTED_INITIAL,
        "expected_accepted_digest_not_checked_by_init_trace": EXPECTED_ACCEPTED,
        "fields": fields,
        "boundary_values": {
            
            "geometry.Aa_outer_3": np.asarray(state.geometry.Aa, dtype=float)[-3:].tolist(),
            "geometry.Lambda_outer_3": np.asarray(state.geometry.Lambda, dtype=float)[-3:].tolist(),
            "geometry.alpha_outer_3": np.asarray(state.geometry.alpha, dtype=float)[-3:].tolist(),
            "geometry.beta_outer_3": np.asarray(state.geometry.beta, dtype=float)[-3:].tolist(),
        },
    }
    output_path.write_text(json.dumps(trace, indent=2, sort_keys=True, allow_nan=False)+"\n", encoding="utf-8")
    return trace

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--engine-root", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--role", default="candidate")
    p.add_argument("--historical-root")
    args = p.parse_args()
    output = Path(args.output); output.parent.mkdir(parents=True, exist_ok=True)
    report_path = output.parent / "init_trace_report.json"
    report = {
        "schema": "initial_state_identity_trace_v1",
        "diagnostic_only": True,
        "expected_initial_digest_immutable": EXPECTED_INITIAL,
        "accepted_digest_not_evaluated": EXPECTED_ACCEPTED,
        "environment": {
            "python": sys.version, "numpy": np.__version__, "platform": platform.platform(),
            "machine": platform.machine(), "processor": platform.processor(),
            "cwd": os.getcwd(), "python_executable": sys.executable,
            "cpuinfo_model_lines": [line.strip() for line in Path("/proc/cpuinfo").read_text(errors="replace").splitlines()
                                    if line.lower().startswith(("model name", "flags", "microcode"))][:12]
                                    if Path("/proc/cpuinfo").exists() else [],
            "env": {k:v for k,v in os.environ.items() if any(s in k.upper() for s in
                    ("BLAS","OPENBLAS","OMP","MKL","NUMEXPR","VECLIB","PYTHONHASHSEED"))},
        },
        "candidate": None, "historical": None, "first_divergent_field": None,
        "status": "trace_started", "errors": [],
        "gates": {"state_identity_passed": False, "paired_rhs_allowed": False,
                  "scientific_admission": False, "acceptance_oracle_satisfied": False},
    }
    try:
        report["candidate"] = collect(args.engine_root, output, args.role)
    except BaseException as exc:
        report["errors"].append({"role":"candidate", "error":repr(exc), "traceback":traceback.format_exc()})
    if args.historical_root:
        try:
            hist_path = output.parent / "historical_init_trace.json"
            child = subprocess.run([sys.executable, str(Path(__file__).resolve()), "--engine-root", args.historical_root,
                                    "--output", str(hist_path), "--role", "historical-child"],
                                   check=False, capture_output=True, text=True, timeout=600)
            if hist_path.exists():
                report["historical"] = json.loads(hist_path.read_text(encoding="utf-8"))
                report["historical_child_exit_code"] = child.returncode
                report["historical_child_stdout"] = child.stdout[-4000:]
                report["historical_child_stderr"] = child.stderr[-4000:]
            else:
                raise RuntimeError("historical trace child produced no trace JSON: " + child.stderr[-4000:])
        except BaseException as exc:
            report["errors"].append({"role":"historical", "error":repr(exc), "traceback":traceback.format_exc()})
    cand, hist = report["candidate"], report["historical"]
    if cand and hist and cand.get("fields") and hist.get("fields"):
        for name, *_ in FIELD_SPECS:
            if cand["fields"].get(name, {}).get("canonical_sha256") != hist["fields"].get(name, {}).get("canonical_sha256"):
                report["first_divergent_field"] = {
                    "field": name, "candidate": cand["fields"][name], "historical": hist["fields"][name]
                }
                break
    report["status"] = "trace_complete" if cand and (hist or not args.historical_root) else "trace_incomplete"
    report["digest_comparison"] = {
        "candidate_initial_matches_historical_reference": bool(cand and cand.get("initial_digest")==EXPECTED_INITIAL),
        "historical_initial_matches_historical_reference": bool(hist and hist.get("initial_digest")==EXPECTED_INITIAL),
    }
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False)+"\n", encoding="utf-8")
    print(json.dumps({"status":report["status"], "candidate_digest":(cand or {}).get("initial_digest"),
                      "historical_digest":(hist or {}).get("initial_digest"),
                      "first_divergent_field":report["first_divergent_field"],
                      "errors":report["errors"], "report":str(report_path)}, sort_keys=True), flush=True)
    # Fail closed after report is durable. The workflow uploads with if: always().
    if report["errors"] or not cand or (args.historical_root and not hist):
        return 2
    return 0 if cand.get("initial_digest")==EXPECTED_INITIAL and (not hist or hist.get("initial_digest")==EXPECTED_INITIAL) else 1

if __name__ == "__main__":
    raise SystemExit(main())
