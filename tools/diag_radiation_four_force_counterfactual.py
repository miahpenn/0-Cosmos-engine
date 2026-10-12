#!/usr/bin/env python3
"""Diagnostic-only radiation four-force counterfactual with centre-state capture.

Four cases: N=80/160 x beta_rad=0/-0.04. One ordinary production step only.
Monkeypatch wrappers are observational: they call the original functions and
return their results unchanged. This is not a physics validation campaign.
"""
from __future__ import annotations

import hashlib
import json
import math
import traceback
from pathlib import Path

import numpy as np

import engine.matter_system as ms
from engine.production_kernel import V55ProductionKernel
from engine.v55_matter import metric_slice_from_q

RESOLUTIONS = (80, 160)
BETAS = (0.0, -0.04)
CFL = 0.03
CELLS = (0, 1, 2)


def _number(value):
    try:
        x = float(value)
        return x if math.isfinite(x) else str(x)
    except Exception:
        return str(value)


def _primitive(q, metric):
    W = q.lorentz()
    E = metric.sqrt_gamma
    # Conserved-state ratio is added by the caller where raw U is available.
    return {
        "rho": _number(q.rho), "pressure": _number(q.pressure),
        "v_r": _number(q.v_r), "W": _number(W),
        "h": _number(q.rho + q.pressure),
        "vhat": _number(math.sqrt(metric.gamma_rr) * q.v_r),
    }


def _metric_capture(state, metric_slice, metrics):
    out = {"grid_centres": [float(x) for x in state.grid.centers],
           "dr": float(state.grid.centers[1] - state.grid.centers[0]),
           "cells": {}}
    U = state.matter.radiation
    prims = ms.primitives(metrics, U, ms.Species.RADIATION)
    for i in CELLS:
        m = metrics[i]
        E = float(U.energy_t[i] / m.sqrt_gamma)
        S = float(U.momentum_r[i] / m.sqrt_gamma)
        ratio = math.sqrt(max(m.gamma_rr_inv * S * S, 0.0)) / max(abs(E), 1e-300)
        out["cells"][str(i)] = {
            "conserved": {"rest": _number(U.rest[i]), "energy_t": _number(U.energy_t[i]),
                          "momentum_r": _number(U.momentum_r[i])},
            "metric": {"alpha": _number(m.alpha), "beta": _number(m.beta),
                       "a": _number(metric_slice.a[i]), "b": _number(metric_slice.b[i]),
                       "X": _number(metric_slice.X[i]), "K": _number(metric_slice.K[i]),
                       "Aa": _number(metric_slice.Aa[i]),
                       "gamma_rr": _number(m.gamma_rr), "gamma_rr_inv": _number(m.gamma_rr_inv),
                       "sqrt_gamma": _number(m.sqrt_gamma)},
            "primitive": {**_primitive(prims[i], m), "abs_S_over_E": _number(ratio)},
        }
    return out


def _install_capture(capture):
    originals = {name: getattr(ms, name) for name in
                 ("primitives", "_reconstructed_primitive", "_hll_flux", "_valencia_source")}
    source_count = {"n": 0}
    def primitives(metrics, U, species):
        result = originals["primitives"](metrics, U, species)
        if species is ms.Species.RADIATION and capture["rhs_entry"] is None:
            cap = {"cells": {}}
            for i in CELLS:
                m, q = metrics[i], result[i]
                E = float(U.energy_t[i] / m.sqrt_gamma)
                S = float(U.momentum_r[i] / m.sqrt_gamma)
                ratio = math.sqrt(max(m.gamma_rr_inv * S*S, 0.0)) / max(abs(E), 1e-300)
                cap["cells"][str(i)] = {
                    "conserved": {"rest": _number(U.rest[i]), "energy_t": _number(U.energy_t[i]),
                                  "momentum_r": _number(U.momentum_r[i])},
                    "primitive": {**_primitive(q, m), "abs_S_over_E": _number(ratio)},
                    "metric": {"alpha": _number(m.alpha), "beta": _number(m.beta),
                               "gamma_rr": _number(m.gamma_rr), "gamma_rr_inv": _number(m.gamma_rr_inv),
                               "sqrt_gamma": _number(m.sqrt_gamma)},
                }
            capture["rhs_entry"] = cap
        return result

    def reconstructed(prim, i, side, gamma_rr):
        result = originals["_reconstructed_primitive"](prim, i, side, gamma_rr)
        if i in (0, 1) and len(capture["faces"]) < 2:
            capture["faces"].append({"cell_i": i, "side": side, "gamma_rr": _number(gamma_rr),
                                     "primitive": {"rho": _number(result.rho),
                                                   "pressure": _number(result.pressure),
                                                   "v_r": _number(result.v_r),
                                                   "W": _number(result.lorentz())}})
        return result

    def hll(metric, left, right, species):
        result = originals["_hll_flux"](metric, left, right, species)
        if species is ms.Species.RADIATION and capture["first_face"] is None:
            cs = 1.0 / math.sqrt(3.0)
            def speed(q, sign):
                vh = math.sqrt(metric.gamma_rr) * q.v_r
                return metric.alpha * (q.v_r + sign*cs/math.sqrt(metric.gamma_rr)) / (1.0 + sign*vh*cs) - metric.beta
            sm = min(0.0, speed(left, -1.0), speed(right, -1.0))
            sp = max(0.0, speed(left, 1.0), speed(right, 1.0))
            capture["first_face"] = {
                "face_between_cells": [0, 1],
                "metric": {k: _number(getattr(metric, k)) for k in
                           ("alpha", "beta", "gamma_rr", "gamma_rr_inv", "sqrt_gamma")},
                "left": {"rho": _number(left.rho), "pressure": _number(left.pressure),
                         "v_r": _number(left.v_r), "W": _number(left.lorentz())},
                "right": {"rho": _number(right.rho), "pressure": _number(right.pressure),
                          "v_r": _number(right.v_r), "W": _number(right.lorentz())},
                "characteristic_speeds": {"s_minus": _number(sm), "s_plus": _number(sp)},
                "hll_flux": [_number(x) for x in result],
            }
        return result

    def source(metric, q, K, Aa, dr_alpha, dr_beta, dr_grr, dr_gthth):
        result = originals["_valencia_source"](metric, q, K, Aa, dr_alpha, dr_beta, dr_grr, dr_gthth)
        source_count["n"] += 1
        # evolve_species is called in DM, baryon, radiation order; each species
        # calls sources in cell order. Capture the centre source on radiation.
        n = len(capture["grid_centres"])
        if source_count["n"] == 2*n + 1:
            capture["centre_source"] = {
                "energy": _number(result[0]), "momentum": _number(result[1]),
                "K": _number(K), "Aa": _number(Aa),
                "dr_alpha": _number(dr_alpha), "dr_beta": _number(dr_beta),
                "dr_gamma_rr": _number(dr_grr), "dr_gamma_thth": _number(dr_gthth),
                "primitive": {"rho": _number(q.rho), "pressure": _number(q.pressure),
                              "v_r": _number(q.v_r), "W": _number(q.lorentz())},
            }
        return result

    ms.primitives = primitives
    ms._reconstructed_primitive = reconstructed
    ms._hll_flux = hll
    ms._valencia_source = source
    def restore():
        for name, value in originals.items():
            setattr(ms, name, value)
    return restore


def _snapshot(kernel, state, label):
    obs = kernel.diagnostics(state, profiles=True)
    return {
        "label": label, "t": float(state.t), "tau": float(state.tau),
        "H_eff": _number(obs.get("H_eff", np.nan)),
        "max_abs_H_admitted": _number(obs.get("hamiltonian_max_admitted", obs.get("hamiltonian_max", np.nan))),
        "max_abs_M_admitted": _number(obs.get("momentum_max_admitted", obs.get("momentum_max", np.nan))),
        "max_abs_C_Lambda": _number(obs.get("connection_max", np.nan)),
        "MS_current_max": _number(obs.get("misner_sharp_current_max", np.nan)),
        "MS_current_rms": _number(obs.get("misner_sharp_current_rms", np.nan)),
        "rad_margin_i0": _number(obs.get("radiation_margin_i0", np.nan)),
        "rad_margin_admitted": _number(obs.get("radiation_margin_admitted", np.nan)),
    }


def run_one(resolution, beta_rad):
    result = {"resolution": resolution, "beta_rad": beta_rad, "status": "failed",
              "dt": None, "pre": None, "post": None, "capture_A": None,
              "capture_B": {"rhs_entry": None, "faces": [], "first_face": None,
                            "centre_source": None, "grid_centres": []},
              "error_type": None, "error": None, "traceback": None}
    kernel = None
    restore = None
    try:
        kernel = V55ProductionKernel()
        kernel.set_diagnostic_beta_rad(beta_rad)
        state = kernel.initialize(resolution=resolution, r_max=40.0, include_radiation=True)
        metric_slice = metric_slice_from_q(state.grid, state.geometry)
        metrics = ms._metric_arrays(metric_slice)
        result["dt"] = float(CFL * (state.grid.centers[1] - state.grid.centers[0]))
        result["capture_A"] = _metric_capture(state, metric_slice, metrics)
        result["pre"] = _snapshot(kernel, state, "t0")
        result["capture_B"]["grid_centres"] = [float(x) for x in state.grid.centers]
        restore = _install_capture(result["capture_B"])
        state = kernel.step(state, result["dt"])
        result["post"] = _snapshot(kernel, state, "t0+dt")
        result["status"] = "passed"
    except Exception as exc:
        result["error_type"] = type(exc).__name__
        result["error"] = str(exc)
        result["traceback"] = traceback.format_exc()
        print(f"CASE FAILED N={resolution} beta_rad={beta_rad}: {type(exc).__name__}: {exc}")
    finally:
        if restore:
            restore()
    return result


def main():
    results = []
    for resolution in RESOLUTIONS:
        for beta in BETAS:
            print(f"Running N={resolution} beta_rad={beta}")
            results.append(run_one(resolution, beta))
    out = Path("artifacts/radiation_four_force_counterfactual.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(results, indent=2, sort_keys=True, allow_nan=False)
    out.write_text(payload + "\n", encoding="utf-8")
    digest = hashlib.sha256((payload + "\n").encode("utf-8")).hexdigest()
    print(f"Wrote {out} SHA-256={digest}")
    failures = [r for r in results if r["status"] != "passed"]
    print(f"Counterfactual cases: {len(results)-len(failures)} passed, {len(failures)} failed")
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
