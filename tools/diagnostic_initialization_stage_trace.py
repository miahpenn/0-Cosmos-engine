#!/usr/bin/env python3
"""Diagnostic-only stage trace for initialization identity; never changes production arithmetic."""
from __future__ import annotations
import argparse, hashlib, json, linecache, os, subprocess, sys, traceback
from pathlib import Path
import numpy as np

EXPECTED = "a256b57f274ee120b55c0f805ed488d91382a9799647c8f9811d7df82aa337c0"

def describe(value):
    a = np.asarray(value)
    c = np.ascontiguousarray(a)
    return {
        "shape": list(a.shape), "dtype": str(a.dtype),
        "sha256": hashlib.sha256(c.tobytes(order="C")).hexdigest(),
        "min": float(np.min(a)) if a.size else None,
        "max": float(np.max(a)) if a.size else None,
        "first3": a.reshape(-1)[:3].tolist(),
        "last3": a.reshape(-1)[-3:].tolist(),
    }

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--engine-root", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--role", default="candidate")
    p.add_argument("--historical-root")
    a = p.parse_args()
    root, out = Path(a.engine_root).resolve(), Path(a.output).resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, str(root))
    report = {"schema":"init_stage_trace_v1", "diagnostic_only":True, "role":a.role,
              "engine_root":str(root), "expected_digest_immutable":EXPECTED,
              "stages":[], "errors":[], "first_divergent_stage":None}
    arrays = {}
    pass_counter = {"n":0}

    def snap(stage, values):
        item = {"stage":stage, "arrays":{}}
        for key in ("S","D","PS","PD","phi","Pi","Sp","Dp","Kt","Ktp","j",
                    "rhoS","rhoD","rhoPhi","rhoMatter","rho_em","rho","Kr","src","integ","B","A","K","Aa","a","b","X","Lambda"):
            v = values.get(key)
            if isinstance(v, np.ndarray) and v.ndim > 0:
                item["arrays"][key] = describe(v)
                arrays[f"{stage}__{key}"] = np.asarray(v)
        report["stages"].append(item)

    def tracer(frame, event, arg):
        if frame.f_code.co_name != "build_initial_data" or "v55_initial.py" not in frame.f_code.co_filename:
            return tracer
        loc = frame.f_locals
        source = linecache.getline(frame.f_code.co_filename, frame.f_lineno).strip()
        if event == "line" and source.startswith("rhoS ="):
            if "src" in loc:
                pass_counter["n"] += 1
                snap(f"reconstruction_pass_{pass_counter['n']}", loc)
        elif event == "line" and source.startswith("if np.any(~np.isfinite(B))"):
            if "src" in loc:
                pass_counter["n"] = 6
                snap("reconstruction_pass_6", loc)
        elif event == "return":
            state = loc.get("state")
            if state is not None:
                snap("build_initial_data_return", {
                    **loc, "a":state.a, "b":state.b, "X":state.X, "alpha":state.alpha,
                    "beta":state.beta, "Aa":state.Aa, "K":state.K, "Lambda":state.Lambda, "B":state.B
                })
        return tracer

    try:
        from engine.production_kernel import V55ProductionKernel
        kernel = V55ProductionKernel()
        original_solve = kernel._solve_lapse_for_state
        solve_counter = {"n":0}
        def wrapped_solve(grid, geometry, scalars, matter, electromagnetic):
            solve_counter["n"] += 1
            n = solve_counter["n"]
            before = {"a":geometry.a,"b":geometry.b,"X":geometry.X,"alpha":geometry.alpha,
                      "Aa":geometry.Aa,"K":geometry.K,"Lambda":geometry.Lambda,"B":geometry.B,
                      "S":scalars.S,"PS":scalars.PS,"D":scalars.D,"PD":scalars.PD,
                      "phi":scalars.phi,"Pi":scalars.Pi}
            snap(f"lapse_solve_{n}_before", before)
            result = original_solve(grid, geometry, scalars, matter, electromagnetic)
            alpha = result[0] if isinstance(result, tuple) else result
            snap(f"lapse_solve_{n}_after", {**before, "alpha":np.asarray(alpha)})
            return result
        kernel._solve_lapse_for_state = wrapped_solve
        sys.settrace(tracer)
        try:
            state = kernel.initialize(resolution=60, r_max=80.0, include_radiation=True)
        finally:
            sys.settrace(None)
        snap("initialize_return", {"a":state.geometry.a,"b":state.geometry.b,"X":state.geometry.X,
              "alpha":state.geometry.alpha,"beta":state.geometry.beta,"Aa":state.geometry.Aa,
              "K":state.geometry.K,"Lambda":state.geometry.Lambda,"B":state.geometry.B,
              "S":state.scalars.S,"PS":state.scalars.PS,"D":state.scalars.D,"PD":state.scalars.PD,
              "phi":state.scalars.phi,"Pi":state.scalars.Pi})
        from tools.audit_current_source_momentum_constraint_tangent import state_digest
        digest = state_digest(state)
        report["initial_digest"] = digest
        report["matches_immutable_digest"] = digest == EXPECTED
        report["alpha_uniform_one_at_build_return"] = None
        for st in report["stages"]:
            if st["stage"] == "build_initial_data_return" and "alpha" in st["arrays"]:
                report["alpha_uniform_one_at_build_return"] = (
                    st["arrays"]["alpha"]["min"] == 1.0 and st["arrays"]["alpha"]["max"] == 1.0
                )
        np.savez_compressed(out.with_suffix(".npz"), **arrays)
    except BaseException as exc:
        report["errors"].append({"error":repr(exc),"traceback":traceback.format_exc()})
    if a.historical_root:
        child_out = out.parent / "historical_stages.json"
        child = subprocess.run([sys.executable, str(Path(__file__).resolve()), "--engine-root",
             str(Path(a.historical_root).resolve()), "--output", str(child_out), "--role","historical-child"],
             capture_output=True, text=True, check=False, timeout=600)
        report["historical_child_exit_code"] = child.returncode
        if child_out.exists():
            hist = json.loads(child_out.read_text())
            report["historical_digest"] = hist.get("initial_digest")
            report["historical_stages"] = hist.get("stages", [])
            report["historical_errors"] = hist.get("errors", [])
            candidate_by_name = {s["stage"]:s for s in report["stages"]}
            for hs in report["historical_stages"]:
                cs = candidate_by_name.get(hs["stage"])
                if cs is None:
                    report["first_divergent_stage"] = {"stage":hs["stage"],"reason":"missing_candidate_stage"}
                    break
                keys = sorted(set(cs["arrays"]) | set(hs["arrays"]))
                diff = next((k for k in keys if cs["arrays"].get(k,{}).get("sha256") != hs["arrays"].get(k,{}).get("sha256")), None)
                if diff:
                    report["first_divergent_stage"] = {"stage":hs["stage"],"array":diff,
                        "candidate":cs["arrays"].get(diff),"historical":hs["arrays"].get(diff)}
                    break
        else:
            report["errors"].append({"error":"historical child did not produce trace","stderr":child.stderr[-4000:]})
    out.write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print(json.dumps({"role":a.role,"digest":report.get("initial_digest"),
      "matches_immutable_digest":report.get("matches_immutable_digest"),
      "first_divergent_stage":report.get("first_divergent_stage"),
      "stage_count":len(report["stages"]),"errors":report["errors"],"output":str(out)},sort_keys=True),flush=True)
    return 2 if report["errors"] else (0 if report.get("matches_immutable_digest") else 1)

if __name__ == "__main__":
    raise SystemExit(main())
