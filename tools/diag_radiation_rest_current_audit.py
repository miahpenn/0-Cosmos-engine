"""Diagnostic-only audit of radiation's otherwise-unused rest-current channel.

Compare the existing generic finite-volume transport with the existing
radiation-specific flux/conserved helpers on the same one-step production state.
No production module is changed by this runner; the second arm monkeypatches
only the diagnostic process and restores all functions before saving output.
"""
from __future__ import annotations

import hashlib
import json
import math
import traceback
from pathlib import Path

import numpy as np

from engine import matter_system as ms
from engine.production_kernel import V55ProductionKernel
from engine.v55_matter import metric_slice_from_q

RESOLUTIONS = (80, 160)
CFL = 0.03
VENDOR = "d6052d605673ce9d82cdc99b0df79855fa2ea215"


def number(value):
    value = float(value)
    return value if math.isfinite(value) else ("nan" if math.isnan(value) else ("inf" if value > 0 else "-inf"))


def radiation_state_profile(state):
    metric = metric_slice_from_q(state.grid, state.geometry)
    metrics = ms._metric_arrays(metric)
    rad = state.matter.radiation
    sqrt_gamma = np.asarray([m.sqrt_gamma for m in metrics], dtype=float)
    gamma_rr_inv = np.asarray([m.gamma_rr_inv for m in metrics], dtype=float)
    E = np.asarray(rad.energy_t, dtype=float) / sqrt_gamma
    S = np.asarray(rad.momentum_r, dtype=float) / sqrt_gamma
    abs_S = np.sqrt(np.maximum(gamma_rr_inv * S * S, 0.0))
    ratio = abs_S / np.maximum(np.abs(E), 1.0e-300)
    cone_margin = E - abs_S
    r = np.asarray(state.grid.centers, dtype=float)
    finite = np.isfinite(E) & np.isfinite(abs_S) & np.isfinite(ratio)
    inadmissible = (~finite) | (E < 0.0) | (abs_S > E)
    rest = np.asarray(rad.rest, dtype=float)
    max_ratio_i = int(np.argmax(np.where(finite, ratio, np.inf)))
    max_abs_rest_i = int(np.argmax(np.abs(rest)))
    rows = [
        {
            "i": int(i), "r": number(r[i]),
            "rest": number(rest[i]),
            "energy_t": number(rad.energy_t[i]),
            "momentum_r": number(rad.momentum_r[i]),
            "E": number(E[i]),
            "abs_S": number(abs_S[i]),
            "abs_S_over_E": number(ratio[i]),
            "ratio_margin": number(1.0 - ratio[i]),
            "cone_margin": number(cone_margin[i]),
        }
        for i in range(len(r))
    ]
    return {
        "summary": {
            "cell_count": int(len(r)),
            "finite_cell_count": int(np.count_nonzero(finite)),
            "inadmissible_cell_count": int(np.count_nonzero(inadmissible)),
            "rest_nonzero_cell_count": int(np.count_nonzero(rest != 0.0)),
            "rest_max_abs": number(np.max(np.abs(rest))),
            "rest_max_abs_cell": max_abs_rest_i,
            "rest_max_abs_radius": number(r[max_abs_rest_i]),
            "rest_min": number(np.min(rest)),
            "rest_max": number(np.max(rest)),
            "max_abs_S_over_E": number(ratio[max_ratio_i]),
            "max_ratio_cell": max_ratio_i,
            "max_ratio_radius": number(r[max_ratio_i]),
            "min_ratio_margin": number(np.min(1.0 - ratio)),
            "min_cone_margin": number(np.min(cone_margin)),
            "centre_abs_S_over_E": number(ratio[0]),
            "outer_abs_S_over_E": number(ratio[-1]),
        },
        "profile": rows,
    }


def geometry_profile(state):
    names = ("a", "b", "X", "alpha", "Aa", "K", "Lambda")
    return {
        name: [number(x) for x in np.asarray(getattr(state.geometry, name), dtype=float)]
        for name in names if hasattr(state.geometry, name)
    }


def diagnostic_summary(kernel, state):
    obs = kernel.diagnostics(state, profiles=True)
    return {
        "t": number(state.t), "tau": number(state.tau),
        "H_eff": number(obs.get("H_eff", np.nan)),
        "max_abs_H_admitted": number(obs.get("hamiltonian_max_admitted", obs.get("hamiltonian_max", np.nan))),
        "max_abs_M_admitted": number(obs.get("momentum_max_admitted", obs.get("momentum_max", np.nan))),
        "max_abs_C_Lambda": number(obs.get("connection_max", np.nan)),
        "MS_current_max": number(obs.get("misner_sharp_current_max", np.nan)),
        "MS_current_rms": number(obs.get("misner_sharp_current_rms", np.nan)),
    }


def install_radiation_channel_dispatch():
    """Diagnostic override: use existing radiation-specific helpers for Species.RADIATION."""
    active = {"species": None}
    originals = {
        "hll": ms._hll_flux,
        "valencia_flux": ms._valencia_flux,
        "valencia_conserved": ms._valencia_conserved,
    }

    def hll(metric, left, right, species):
        active["species"] = species
        return originals["hll"](metric, left, right, species)

    def valencia_flux(metric, q):
        if active["species"] is ms.Species.RADIATION:
            return ms._radiation_flux(metric, q)
        return originals["valencia_flux"](metric, q)

    def valencia_conserved(metric, q):
        if active["species"] is ms.Species.RADIATION:
            return ms._radiation_conserved(metric, q)
        return originals["valencia_conserved"](metric, q)

    ms._hll_flux = hll
    ms._valencia_flux = valencia_flux
    ms._valencia_conserved = valencia_conserved

    def restore():
        ms._hll_flux = originals["hll"]
        ms._valencia_flux = originals["valencia_flux"]
        ms._valencia_conserved = originals["valencia_conserved"]
    return restore


def run_case(resolution, mode):
    result = {
        "resolution": resolution, "mode": mode, "status": "failed",
        "dt": None, "pre_summary": None, "post_summary": None,
        "pre_radiation": None, "post_radiation": None,
        "pre_geometry": None, "post_geometry": None,
        "diagnostics_pre": None, "diagnostics_post": None,
        "error_type": None, "error": None, "traceback": None,
    }
    restore = None
    try:
        kernel = V55ProductionKernel()
        state = kernel.initialize(resolution=resolution, r_max=40.0, include_radiation=True)
        result["dt"] = CFL * float(state.grid.centers[1] - state.grid.centers[0])
        result["pre_radiation"] = radiation_state_profile(state)
        result["pre_geometry"] = geometry_profile(state)
        result["diagnostics_pre"] = diagnostic_summary(kernel, state)
        if mode == "radiation_specific_dispatch":
            restore = install_radiation_channel_dispatch()
        advanced = kernel.step(state, result["dt"])
        result["post_radiation"] = radiation_state_profile(advanced)
        result["post_geometry"] = geometry_profile(advanced)
        result["diagnostics_post"] = diagnostic_summary(kernel, advanced)
        result["status"] = "passed"
    except Exception as exc:
        result["error_type"] = type(exc).__name__
        result["error"] = str(exc)
        result["traceback"] = traceback.format_exc()
        print(f"CASE FAILED N={resolution} mode={mode}: {type(exc).__name__}: {exc}", flush=True)
    finally:
        if restore is not None:
            restore()
    return result


def max_delta(left, right, field):
    rows_a = left["post_radiation"]["profile"]
    rows_b = right["post_radiation"]["profile"]
    diffs = [
        (abs(float(b[field]) - float(a[field])), float(b[field]) - float(a[field]),
         int(a["i"]), float(a["r"]), float(a[field]), float(b[field]))
        for a, b in zip(rows_a, rows_b)
    ]
    peak = max(diffs, key=lambda x: x[0])
    return {
        "max_abs_delta": number(peak[0]),
        "signed_delta_specific_minus_generic": number(peak[1]),
        "cell": peak[2], "radius": number(peak[3]),
        "generic_value": number(peak[4]), "specific_value": number(peak[5]),
    }


def compare_pair(generic, specific):
    pair = {
        "resolution": generic["resolution"],
        "generic_status": generic["status"],
        "radiation_specific_status": specific["status"],
        "comparison_available": False,
        "radiation_profile_max_abs_delta": {},
        "geometry_profile_max_abs_delta": {},
        "diagnostics_specific_minus_generic": {},
    }
    if generic["status"] != "passed" or specific["status"] != "passed":
        return pair
    pair["comparison_available"] = True
    for field in ("rest", "energy_t", "momentum_r", "E", "abs_S", "abs_S_over_E", "ratio_margin", "cone_margin"):
        pair["radiation_profile_max_abs_delta"][field] = max_delta(generic, specific, field)
    ga, gb = generic["post_geometry"], specific["post_geometry"]
    for field in ga:
        vals = [
            (abs(float(b) - float(a)), float(b) - float(a), i, float(a), float(b))
            for i, (a, b) in enumerate(zip(ga[field], gb[field]))
        ]
        peak = max(vals, key=lambda x: x[0])
        pair["geometry_profile_max_abs_delta"][field] = {
            "max_abs_delta": number(peak[0]),
            "signed_delta_specific_minus_generic": number(peak[1]),
            "cell": peak[2], "generic_value": number(peak[3]), "specific_value": number(peak[4]),
        }
    for field in ("H_eff", "max_abs_H_admitted", "max_abs_M_admitted", "max_abs_C_Lambda"):
        pair["diagnostics_specific_minus_generic"][field] = number(
            float(specific["diagnostics_post"][field]) - float(generic["diagnostics_post"][field])
        )
    pair["expected_invariant_check"] = {
        "specific_rest_max_abs": specific["post_radiation"]["summary"]["rest_max_abs"],
        "generic_rest_max_abs": generic["post_radiation"]["summary"]["rest_max_abs"],
        "all_energy_momentum_channels_unchanged_exactly": all(
            float(pair["radiation_profile_max_abs_delta"][field]["max_abs_delta"]) == 0.0
            for field in ("energy_t", "momentum_r", "E", "abs_S", "abs_S_over_E", "ratio_margin", "cone_margin")
        ),
        "all_geometry_fields_unchanged_exactly": all(
            float(pair["geometry_profile_max_abs_delta"][field]["max_abs_delta"]) == 0.0
            for field in pair["geometry_profile_max_abs_delta"]
        ),
    }
    return pair


def main():
    cases = []
    pairs = []
    for resolution in (80, 160):
        for mode in ("generic_all_species_dispatch", "radiation_specific_dispatch"):
            print(f"Running N={resolution} mode={mode}", flush=True)
            cases.append(run_case(resolution, mode))
        generic = next(c for c in cases if c["resolution"] == resolution and c["mode"] == "generic_all_species_dispatch")
        specific = next(c for c in cases if c["resolution"] == resolution and c["mode"] == "radiation_specific_dispatch")
        pairs.append(compare_pair(generic, specific))
    payload_obj = {
        "schema": "radiation_rest_current_channel_audit_v1",
        "diagnostic_only": True,
        "production_physics_changed": False,
        "vendor_pin": VENDOR,
        "cases": cases,
        "matched_pairs": pairs,
        "interpretation_limit": (
            "This isolates the passive first/rest channel. Radiation recovery and stress-energy projection use E/S; "
            "the run does not establish causality for the late-time radiation-cone failure."
        ),
    }
    out = Path("artifacts/radiation_rest_current_channel_audit.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(payload_obj, indent=2, sort_keys=True, allow_nan=False) + "\n"
    out.write_text(text, encoding="utf-8")
    digest = hashlib.sha256(text.encode("utf-8")).hexdigest()
    print(f"JSON_SHA256={digest}", flush=True)
    failures = [case for case in cases if case["status"] != "passed"]
    invariant_failures = [
        pair for pair in pairs if pair["comparison_available"] and (
            not pair["expected_invariant_check"]["all_energy_momentum_channels_unchanged_exactly"]
            or not pair["expected_invariant_check"]["all_geometry_fields_unchanged_exactly"]
            or float(pair["expected_invariant_check"]["specific_rest_max_abs"]) != 0.0
        )
    ]
    print(f"Cases: {len(cases)-len(failures)} passed, {len(failures)} failed", flush=True)
    print("PAIRS=" + json.dumps(pairs, sort_keys=True), flush=True)
    if failures or invariant_failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
