#!/usr/bin/env python3
"""Research-only operand-wise fixed-R interpolation discriminator.

Compare two ways of evaluating the same pinned cell-centred operands at a
common fixed areal radius:
  (1) interpolate the already assembled chi and F=4*pi*R**2*T_t^R;
  (2) interpolate the underlying geometry/stress operands, then reconstruct
      chi and the transformed current at the query radius.

Run both linear (low-order control) and four-point cubic interpolation, at
dt and dt/2, on copied one-step states from the same accepted t=0.60 input.
This isolates interpolation-before-versus-after-nonlinear-reassembly. It does
not recompute R_r from interpolated a,b,X, nor does it claim physical closure.
Production equations and pinned inputs remain untouched.
"""
from __future__ import annotations

import hashlib
import json
import os
import platform
import sys
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
os.environ.setdefault("COSMOS_PINNED_SOURCE_ROOT", str(ROOT))

from engine.v55_matter import total_matter_projection
from tools import audit_bssn_areal_misner_sharp_current as areal
from tools import audit_full_profile_propagation_candidate as saved
from tools import audit_origin_cell_joint_shadow_one_step as origin_shadow

SOURCE = "4aec88d0e3ff028e0a65a8ac36b30948625df792"
VENDOR = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
N = int(os.environ.get("AREAL_OPERAND_RESOLUTION", "120"))
FULL_DIR = Path(os.environ.get("FULL_PROFILE_ARTIFACT_DIR", f"reference/full/N{N}"))
OUT = Path(os.environ.get("AREAL_OPERAND_OUT", f"results/bssn-areal-operand-interpolation-N{N}-2026-10-25"))
TARGET_TIME = 0.60
CFL = 0.03
FOUR_PI = 4.0 * np.pi
FIXED_RADII = (0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 10.0)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def finite(name: str, values: Any) -> np.ndarray:
    arr = np.asarray(values, dtype=np.float64)
    if not np.all(np.isfinite(arr)):
        raise RuntimeError(f"{name} contains non-finite values")
    return arr


def scalar_finite(name: str, value: Any) -> float:
    out = float(value)
    if not np.isfinite(out):
        raise RuntimeError(f"{name} is non-finite")
    return out


def operand_arrays(profile: dict[str, Any]) -> dict[str, np.ndarray]:
    """Collect raw fields used by chi and the transformed mixed current."""
    q = profile["state"]
    g = q.geometry
    proj = total_matter_projection(q.grid, g, q.scalars, q.matter)
    arrays = {
        "R": finite("R", profile["R"]),
        "a": finite("a", g.a),
        "X": finite("X", g.X),
        "alpha": finite("alpha", g.alpha),
        "Rr": finite("Rr", profile["Rr"]),
        "Rt": finite("Rt", profile["Rt"]),
        "E": finite("Eulerian energy", proj["rho"]),
        "j": finite("covariant radial current", proj["j"]),
        "pr": finite("radial pressure", proj["pr"]),
    }
    for name, vals in arrays.items():
        if vals.ndim != 1:
            raise RuntimeError(f"{name} is not a one-dimensional cell profile")
    if len({len(v) for v in arrays.values()}) != 1:
        raise RuntimeError("raw operands do not share the same cell grid")
    if np.any(arrays["a"] <= 0) or np.any(arrays["X"] <= 0) or np.any(arrays["alpha"] <= 0) or np.any(arrays["Rr"] <= 0):
        raise RuntimeError("positive geometry precondition failed")
    return arrays


def recombine_at_query(fields: dict[str, float], radius: float) -> tuple[float, float]:
    a = fields["a"]
    X = fields["X"]
    alpha = fields["alpha"]
    Rr = fields["Rr"]
    Rt = fields["Rt"]
    E = fields["E"]
    j = fields["j"]
    pr = fields["pr"]
    grr = a / (X * X)
    j_up = j / grr
    chi = (X * X / a) * Rr * Rr - (Rt / alpha) ** 2
    T_t_R = -E * Rt - alpha * j_up * Rr - (Rt / Rr) * (j / alpha * Rt + pr * Rr)
    F = FOUR_PI * radius * radius * T_t_R
    return scalar_finite("operand-reassembled chi", chi), scalar_finite("operand-reassembled current", F)


def query_operand(profile_arrays: dict[str, np.ndarray], radius: float, method: str) -> dict[str, float]:
    sample = {}
    for name in ("a", "X", "alpha", "Rr", "Rt", "E", "j", "pr"):
        if method == "linear":
            value = np.interp(radius, profile_arrays["R"], profile_arrays[name])
        elif method == "cubic":
            value = areal.interp4(radius, profile_arrays["R"], profile_arrays[name])
        else:
            raise ValueError(f"unsupported interpolation method {method}")
        sample[name] = scalar_finite(f"{method}.{name}@R={radius}", value)
    chi, current = recombine_at_query(sample, radius)
    return {"chi": chi, "current": current, **sample}


def query_composite(profile: dict[str, Any], radius: float, method: str) -> dict[str, float]:
    if method == "linear":
        chi = np.interp(radius, profile["R"], profile["chi"])
        current = np.interp(radius, profile["R"], profile["rhs"])
    elif method == "cubic":
        chi = areal.interp4(radius, profile["R"], profile["chi"])
        current = areal.interp4(radius, profile["R"], profile["rhs"])
    else:
        raise ValueError(f"unsupported interpolation method {method}")
    return {
        "chi": scalar_finite(f"{method}.composite_chi@R={radius}", chi),
        "current": scalar_finite(f"{method}.composite_current@R={radius}", current),
    }


def verify_native_recombination(profile: dict[str, Any], arrays: dict[str, np.ndarray]) -> dict[str, float]:
    raw = {}
    for name in ("a", "X", "alpha", "Rr", "Rt", "E", "j", "pr"):
        raw[name] = arrays[name]
    a, X, alpha, Rr, Rt = (raw[x] for x in ("a", "X", "alpha", "Rr", "Rt"))
    E, j, pr = (raw[x] for x in ("E", "j", "pr"))
    grr = a / (X * X)
    chi = (X * X / a) * Rr * Rr - (Rt / alpha) ** 2
    current = FOUR_PI * arrays["R"] ** 2 * (
        -E * Rt - alpha * (j / grr) * Rr - (Rt / Rr) * (j / alpha * Rt + pr * Rr)
    )
    chi_gap = finite("native chi recombination gap", chi - finite("observer chi", profile["chi"]))
    current_gap = finite("native current recombination gap", current - finite("observer current", profile["rhs"]))
    out = {
        "max_abs_chi_reassembly_gap_on_native_grid": float(np.max(np.abs(chi_gap))),
        "max_abs_current_reassembly_gap_on_native_grid": float(np.max(np.abs(current_gap))),
    }
    if out["max_abs_chi_reassembly_gap_on_native_grid"] > 1e-10:
        raise RuntimeError(f"native chi reassembly failed: {out}")
    if out["max_abs_current_reassembly_gap_on_native_grid"] > 1e-10:
        raise RuntimeError(f"native current reassembly failed: {out}")
    return out


def residual_summary(residual: np.ndarray, radii: np.ndarray, current: np.ndarray) -> dict[str, float]:
    residual = finite("residual", residual)
    current = finite("current-average scale", current)
    inner = radii <= 3.0
    annulus = (radii >= 2.0) & (radii <= 10.0)
    if not inner.any() or not annulus.any():
        raise RuntimeError("common sample does not cover inner and annular regions")
    rms = float(np.sqrt(np.mean(residual[annulus] ** 2)))
    scale = float(np.sqrt(np.mean(current[annulus] ** 2)))
    return {
        "max_abs_R_le_3": float(np.max(np.abs(residual[inner]))),
        "rms_R_2_to_10": rms,
        "max_abs_R_2_to_10": float(np.max(np.abs(residual[annulus]))),
        "rms_relative_to_current_R_2_to_10": rms / max(scale, 1e-300),
    }


def one_arm(kernel, accepted, h: float) -> dict[str, Any]:
    accepted_digest = saved.tangent.state_digest(accepted)
    start = areal.geometry_rate_and_mass(kernel, accepted)
    state1 = kernel.step(areal.clone_state(accepted), h)
    state1_digest = saved.tangent.state_digest(state1)
    end = areal.geometry_rate_and_mass(kernel, state1)
    if saved.tangent.state_digest(accepted) != accepted_digest:
        raise RuntimeError("accepted input was mutated by the one-step observer")
    if saved.tangent.state_digest(state1) != state1_digest:
        raise RuntimeError("copied evolved state was mutated by observer")
    if abs(float(state1.t) - (float(accepted.t) + h)) > 1e-10:
        raise RuntimeError("shadow step did not advance exactly by dt")
    if max(start["max_abs_beta"], end["max_abs_beta"]) > 1e-14:
        raise RuntimeError("zero-shift precondition failed")

    start_arrays = operand_arrays(start)
    end_arrays = operand_arrays(end)
    native_checks = {
        "start": verify_native_recombination(start, start_arrays),
        "end": verify_native_recombination(end, end_arrays),
    }
    lo = max(float(start["R"][0]), float(end["R"][0]))
    hi = min(float(start["R"][-1]), float(end["R"][-1]), 10.0)
    if not hi > lo:
        raise RuntimeError("no common fixed-areal-radius interval")
    radii = np.linspace(lo, hi, max(300, int(accepted.grid.n) * 2))

    methods = (
        ("composite_linear", "linear", "composite"),
        ("operand_linear", "linear", "operand"),
        ("composite_cubic", "cubic", "composite"),
        ("operand_cubic", "cubic", "operand"),
    )
    sampled: dict[str, dict[str, np.ndarray]] = {}
    for label, order, source in methods:
        chi0 = np.empty_like(radii)
        chi1 = np.empty_like(radii)
        f0 = np.empty_like(radii)
        f1 = np.empty_like(radii)
        for i, radius in enumerate(radii):
            if source == "composite":
                q0 = query_composite(start, float(radius), order)
                q1 = query_composite(end, float(radius), order)
            else:
                q0 = query_operand(start_arrays, float(radius), order)
                q1 = query_operand(end_arrays, float(radius), order)
            chi0[i], chi1[i] = q0["chi"], q1["chi"]
            f0[i], f1[i] = q0["current"], q1["current"]
        rate = finite(f"{label}.mass_rate", -0.5 * radii * (chi1 - chi0) / h)
        current_average = finite(f"{label}.current_average", 0.5 * (f0 + f1))
        residual = finite(f"{label}.residual", rate - current_average)
        sampled[label] = {
            "chi_start": chi0,
            "chi_end": chi1,
            "current_start": f0,
            "current_end": f1,
            "mass_rate": rate,
            "current_average": current_average,
            "residual": residual,
        }

    summary = {}
    for label, values in sampled.items():
        summary[label] = {
            "residual": residual_summary(values["residual"], radii, values["current_average"]),
            "fixed_R": {},
        }
        for R in FIXED_RADII:
            if not lo <= R <= hi:
                continue
            idx = int(np.argmin(np.abs(radii - R)))
            summary[label]["fixed_R"][str(R)] = {
                "sample_radius": float(radii[idx]),
                "mass_rate": float(values["mass_rate"][idx]),
                "current_average": float(values["current_average"][idx]),
                "residual": float(values["residual"][idx]),
            }

    paired = {}
    for order in ("linear", "cubic"):
        composite = sampled[f"composite_{order}"]
        operand = sampled[f"operand_{order}"]
        for quantity in ("mass_rate", "current_average", "residual"):
            delta = finite(f"{order}.{quantity}.operand_minus_composite", operand[quantity] - composite[quantity])
            paired[f"{order}_{quantity}_operand_minus_composite"] = {
                "rms_R_2_to_10": float(np.sqrt(np.mean(delta[(radii >= 2.0) & (radii <= 10.0)] ** 2))),
                "max_abs_R_2_to_10": float(np.max(np.abs(delta[(radii >= 2.0) & (radii <= 10.0)]))),
            }

    return {
        "step_h": float(h),
        "base_time": float(accepted.t),
        "end_time": float(state1.t),
        "accepted_input_immutable": True,
        "evolved_copy_immutable_during_observation": True,
        "common_areal_radius_interval": [float(lo), float(hi)],
        "common_samples": int(radii.size),
        "native_operand_reassembly_checks": native_checks,
        "methods": summary,
        "operand_minus_composite": paired,
    }


def main() -> int:
    if N not in (60, 120, 240):
        raise ValueError("resolution must be 60, 120 or 240")
    pins = origin_shadow.verify_source_pins()
    if pins.get("source_commit") != SOURCE or pins.get("vendor_commit") != VENDOR:
        raise RuntimeError("source/vendor pin mismatch")
    if pins.get("protected_engine_vendor_dependency_diff_empty") is not True:
        raise RuntimeError("protected source/vendor/dependency diff is not empty")

    manifest = saved.check_artifact(FULL_DIR, N, full_profile=True)
    arrays = saved.load_profile_arrays(FULL_DIR)
    prefix = f"N{N}_full_continuous_metric_profile_final_"
    kernel, accepted = saved.reconstruct_state(N, arrays, prefix)
    origin_shadow._install_regular_row(kernel)
    if abs(float(accepted.t) - TARGET_TIME) > 1e-12:
        raise RuntimeError("accepted input is not at t=0.60")
    digest = saved.tangent.state_digest(accepted)
    dt = CFL * float(accepted.grid.dr)
    arms = [{"arm": label, **one_arm(kernel, accepted, h)} for label, h in (("dt", dt), ("dt_half", dt / 2.0))]
    if saved.tangent.state_digest(accepted) != digest:
        raise RuntimeError("accepted state digest changed across interpolation arms")

    report = {
        "schema": "bssn-areal-operand-interpolation-discriminator/v1",
        "status": "completed_diagnostic_only",
        "source_pins": {
            "source_commit": SOURCE,
            "vendor_commit": VENDOR,
            "live_pin_check": pins,
            "accepted_artifact_manifest": manifest,
            "accepted_state_digest_before_after_match": True,
        },
        "runtime": {"python": sys.version, "numpy": np.__version__, "platform": platform.platform()},
        "configuration": {
            "resolution": N,
            "accepted_time": TARGET_TIME,
            "dt_rule": "dt=0.03*dr; compare dt and dt/2",
            "same accepted artifact reused": True,
            "time evolution": "one copied forward step per arm; no backward stepping",
            "composite route": "interpolate preassembled chi and F=4*pi*R^2*T_t^R",
            "operand route": "interpolate a, X, alpha, R_r, R_t, Eulerian energy E, radial current j and radial pressure p_r; then recompute chi and F at fixed R",
            "limitation": "R_r and R_t are the existing pinned observer fields; this test does not recompute R_r from interpolated primitive metric fields",
        },
        "arms": arms,
        "interpretation_guard": [
            "Native-grid reassembly is checked against the existing observer before interpolation is varied.",
            "Linear interpolation is retained as a low-order control; cubic interpolation is compared both before and after nonlinear reassembly.",
            "A smaller residual under an interpolation method is diagnostic evidence only and does not prove continuum or total stress-energy conservation.",
            "No production equation, source normalization, initializer, vendor or dependency was changed.",
        ],
        "acceptance": {
            "source_pins_match": True,
            "accepted_artifact_sha_verified": True,
            "accepted_input_immutable": True,
            "native_operand_reassembly_checked": True,
            "linear_low_order_control_included": True,
            "cubic_composite_vs_operand_route_compared": True,
            "integrated_current_identity_closed": False,
            "full_matter_coupled_propagation_identity_closed": False,
            "total_stress_energy_conservation_closed": False,
            "scientific_admission": False,
        },
        "scientific_admission": False,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    report_path = OUT / "report.json"
    report_path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    (OUT / "SHA256SUMS.txt").write_text(f"{sha256(report_path)}  report.json\n", encoding="utf-8")
    print(json.dumps({
        "status": report["status"],
        "resolution": N,
        "source_vendor_pins_match": True,
        "accepted_input_immutable": True,
        "scientific_admission": False,
        "arms": [{
            "arm": arm["arm"],
            "methods": {name: item["residual"] for name, item in arm["methods"].items()},
            "operand_minus_composite": arm["operand_minus_composite"],
        } for arm in arms],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
