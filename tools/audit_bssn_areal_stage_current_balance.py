#!/usr/bin/env python3
"""Research-only comparison of endpoint current quadrature with PIRK stage current.

The production kernel advances a/b/X and the scalar/fluid state with a
predictor/trapezoidal corrector. An endpoint trapezoid, (F(t)+F(t+dt))/2,
is not the same quadrature as the kernel's actual start/predictor-stage pairing.
This diagnostic captures the exact stage geometry/scalar/matter operands passed
to _geometry_stage_terms during an unchanged shadow step, reconstructs the
fixed-R current on those stage states, and compares:

  mass difference over accepted t -> accepted t+dt
  against endpoint-averaged current and kernel-stage-averaged current.

The diagnostic is read-only with respect to production source and the saved
accepted state. It does not use a fitted source or claim that a successful
stage balance alone proves the continuum identity.
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

from tools import audit_bssn_areal_misner_sharp_current as areal
from tools import audit_full_profile_propagation_candidate as saved
from tools import audit_origin_cell_joint_shadow_one_step as origin_shadow

SOURCE = "4aec88d0e3ff028e0a65a8ac36b30948625df792"
VENDOR = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
N = int(os.environ.get("AREAL_STAGE_CURRENT_RESOLUTION", "120"))
FULL_DIR = Path(os.environ.get("FULL_PROFILE_ARTIFACT_DIR", f"reference/full/N{N}"))
OUT = Path(os.environ.get("AREAL_STAGE_CURRENT_OUT", f"results/bssn-areal-stage-current-N{N}-2026-10-25"))
TARGET_TIME = 0.60
CFL = 0.03
FIXED_RADII = (1.0, 2.0, 3.0, 5.0, 8.0, 10.0)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def finite(name: str, values: Any) -> np.ndarray:
    out = np.asarray(values, dtype=np.float64)
    if not np.all(np.isfinite(out)):
        raise RuntimeError(f"{name} contains non-finite values")
    return out


def scalar(name: str, value: Any) -> float:
    out = float(value)
    if not np.isfinite(out):
        raise RuntimeError(f"{name} is non-finite")
    return out


def sample(profile: dict[str, Any], radii: np.ndarray, order: str, field: str) -> np.ndarray:
    values = finite(f"profile {field}", profile[field])
    if order == "linear":
        out = np.interp(radii, finite("profile R", profile["R"]), values)
    elif order == "cubic":
        out = np.asarray([areal.interp4(float(R), profile["R"], values) for R in radii], dtype=np.float64)
    else:
        raise ValueError(f"unknown interpolation order {order}")
    return finite(f"{order} samples of {field}", out)


def summary(values: np.ndarray, radii: np.ndarray) -> dict[str, float]:
    values = finite("summary residual", values)
    inner = radii <= 3.0
    annulus = (radii >= 2.0) & (radii <= 10.0)
    if not inner.any() or not annulus.any():
        raise RuntimeError("common fixed-R sample fails to cover the required regions")
    return {
        "max_abs_R_le_3": float(np.max(np.abs(values[inner]))),
        "rms_R_2_to_10": float(np.sqrt(np.mean(values[annulus] ** 2))),
        "max_abs_R_2_to_10": float(np.max(np.abs(values[annulus]))),
    }


def capture_stage_states(kernel, accepted, dt: float):
    if accepted.electromagnetic is not None:
        raise RuntimeError("stage-current probe is pinned to the current EM-free accepted-state lane")
    accepted_digest = saved.tangent.state_digest(accepted)
    input_state = areal.clone_state(accepted)
    original_method = kernel._geometry_stage_terms
    captured = []

    def wrapped(grid, geometry, scalars, matter, electromagnetic):
        if electromagnetic is not None:
            raise RuntimeError("unexpected electromagnetic stage in the pinned current lane")
        snap = areal.clone_state(accepted)
        snap.grid = grid
        snap.geometry = geometry.copy()
        snap.scalars = scalars.copy()
        snap.matter = matter.copy()
        snap.electromagnetic = None
        captured.append(snap)
        return original_method(grid, geometry, scalars, matter, electromagnetic)

    kernel._geometry_stage_terms = wrapped
    try:
        candidate = kernel.step(input_state, dt)
    finally:
        kernel._geometry_stage_terms = original_method

    if saved.tangent.state_digest(accepted) != accepted_digest:
        raise RuntimeError("saved accepted state was mutated")
    if saved.tangent.state_digest(input_state) != saved.tangent.state_digest(areal.clone_state(accepted)):
        # The kernel's public step contract is pure on the supplied input; fail
        # closed if a future implementation begins mutating the passed state.
        raise RuntimeError("shadow input state was mutated during step")
    if abs(float(candidate.t) - (float(accepted.t) + dt)) > 1e-10:
        raise RuntimeError("accepted shadow candidate did not advance by dt")
    if len(captured) != 3:
        raise RuntimeError(f"expected three explicit-geometry calls (start, predictor, corrected stage), observed {len(captured)}")

    stage0, predictor, stage1 = captured
    if max(np.max(np.abs(q.geometry.beta)) for q in (stage0, predictor, stage1, candidate)) > 1e-14:
        raise RuntimeError("zero-shift precondition failed")
    digests_before = [saved.tangent.state_digest(q) for q in captured]
    return candidate, captured, digests_before


def compare_arm(kernel, accepted, dt: float) -> dict[str, Any]:
    candidate, captured, stage_digests = capture_stage_states(kernel, accepted, dt)
    stage0, predictor, stage1 = captured

    # Reconstruct the same independent observable on accepted, both actual
    # integrator stages, predictor (negative/control), and accepted output.
    profiles = {
        "accepted_start": areal.geometry_rate_and_mass(kernel, accepted),
        "kernel_start_stage": areal.geometry_rate_and_mass(kernel, stage0),
        "forward_euler_predictor": areal.geometry_rate_and_mass(kernel, predictor),
        "kernel_corrected_stage": areal.geometry_rate_and_mass(kernel, stage1),
        "accepted_output": areal.geometry_rate_and_mass(kernel, candidate),
    }
    for name, state in zip(("kernel_start_stage", "forward_euler_predictor", "kernel_corrected_stage"), captured):
        if saved.tangent.state_digest(state) != stage_digests[("kernel_start_stage", "forward_euler_predictor", "kernel_corrected_stage").index(name)]:
            raise RuntimeError(f"{name} was mutated while the observable was reconstructed")

    lo = max(float(p["R"][0]) for p in profiles.values())
    hi = min(float(p["R"][-1]) for p in profiles.values())
    hi = min(hi, 10.0)
    if not hi > lo:
        raise RuntimeError("no common areal-radius range across accepted and stage states")
    radii = np.linspace(lo, hi, max(300, int(accepted.grid.n) * 2))

    method_reports: dict[str, Any] = {}
    for order in ("linear", "cubic"):
        sampled = {
            name: {
                "chi": sample(profile, radii, order, "chi"),
                "current": sample(profile, radii, order, "rhs"),
            }
            for name, profile in profiles.items()
        }
        start = sampled["accepted_start"]
        end = sampled["accepted_output"]
        start_stage = sampled["kernel_start_stage"]
        predictor_values = sampled["forward_euler_predictor"]
        corrected_stage = sampled["kernel_corrected_stage"]

        # Mass rate is measured on the accepted endpoints. The two current
        # quadratures differ only in which current samples represent this step:
        # actual accepted endpoints vs the kernel's start/corrected stage pair.
        mass_rate = finite(f"{order}.accepted-endpoint mass rate",
            -0.5 * radii * (end["chi"] - start["chi"]) / dt)
        endpoint_average = finite(f"{order}.endpoint current average",
            0.5 * (start["current"] + end["current"]))
        stage_average = finite(f"{order}.stage current average",
            0.5 * (start_stage["current"] + corrected_stage["current"]))
        predictor_average = finite(f"{order}.predictor current average",
            0.5 * (start_stage["current"] + predictor_values["current"]))
        endpoint_residual = finite(f"{order}.endpoint residual", mass_rate - endpoint_average)
        stage_residual = finite(f"{order}.stage residual", mass_rate - stage_average)
        predictor_residual = finite(f"{order}.predictor residual", mass_rate - predictor_average)
        quadrature_delta = finite(f"{order}.stage minus endpoint quadrature", stage_average - endpoint_average)

        method_reports[order] = {
            "endpoint_trapezoid_residual": summary(endpoint_residual, radii),
            "kernel_stage_trapezoid_residual": summary(stage_residual, radii),
            "predictor_stage_trapezoid_residual_control": summary(predictor_residual, radii),
            "stage_minus_endpoint_current_quadrature": summary(quadrature_delta, radii),
            "stage_quadrature_residual_delta_vs_endpoint": summary(stage_residual - endpoint_residual, radii),
            "fixed_R": {},
        }
        for R in FIXED_RADII:
            if not lo <= R <= hi:
                continue
            q = {
                name: {
                    field: scalar(f"{order}.{name}.{field}@R={R}", sample(p, np.asarray([R]), order, field)[0])
                    for field in ("chi", "rhs")
                }
                for name, p in profiles.items()
            }
            mass_rate_R = -0.5 * R * (q["accepted_output"]["chi"] - q["accepted_start"]["chi"]) / dt
            end_avg_R = 0.5 * (q["accepted_start"]["rhs"] + q["accepted_output"]["rhs"])
            stage_avg_R = 0.5 * (q["kernel_start_stage"]["rhs"] + q["kernel_corrected_stage"]["rhs"])
            predictor_avg_R = 0.5 * (q["kernel_start_stage"]["rhs"] + q["forward_euler_predictor"]["rhs"])
            method_reports[order]["fixed_R"][str(R)] = {
                "mass_rate": scalar("fixed-R mass rate", mass_rate_R),
                "endpoint_average_current": scalar("fixed-R endpoint-average current", end_avg_R),
                "kernel_stage_average_current": scalar("fixed-R kernel-stage-average current", stage_avg_R),
                "predictor_control_average_current": scalar("fixed-R predictor control current", predictor_avg_R),
                "endpoint_residual": scalar("fixed-R endpoint residual", mass_rate_R - end_avg_R),
                "kernel_stage_residual": scalar("fixed-R stage residual", mass_rate_R - stage_avg_R),
                "predictor_control_residual": scalar("fixed-R predictor control residual", mass_rate_R - predictor_avg_R),
            }

    return {
        "dt": float(dt),
        "accepted_start_time": float(accepted.t),
        "accepted_output_time": float(candidate.t),
        "captured_stage_times": [float(stage0.t), float(predictor.t), float(stage1.t)],
        "captured_stage_labels": ["kernel_start_stage", "forward_euler_predictor", "kernel_corrected_stage"],
        "stage_input_digests_unchanged_after_observation": True,
        "candidate_input_state_immutable": True,
        "common_areal_radius_interval": [float(lo), float(hi)],
        "common_samples": int(radii.size),
        "methods": method_reports,
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
    accepted_digest = saved.tangent.state_digest(accepted)
    dt = CFL * float(accepted.grid.dr)
    arms = [{"arm": label, **compare_arm(kernel, accepted, h)} for label, h in (("dt", dt), ("dt_half", dt / 2.0))]
    if saved.tangent.state_digest(accepted) != accepted_digest:
        raise RuntimeError("accepted state digest changed across both temporal arms")

    report = {
        "schema": "bssn-areal-stage-current-balance/v1",
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
            "dt_rule": "dt=0.03*dr; dt and dt/2",
            "actual_kernel_update": "TRUE-PIRK explicit a/b/X and scalar/fluid variables use the captured start and corrected-stage trapezoidal pairing; internal predictor is retained as a control",
            "endpoint_quadrature": "(F(accepted start)+F(accepted output))/2",
            "kernel_stage_quadrature": "(F(start stage)+F(corrected predictor stage))/2",
            "fixed_R_mass_rate": "-R*(chi(accepted output)-chi(accepted start))/(2*dt)",
            "stage_capture": "temporarily wrap the pure _geometry_stage_terms observer during a shadow step; copy stage operands immediately; restore the original bound method before diagnostic reconstruction",
            "interpolation": "linear control and four-point cubic, same fixed areal-radius sample grid",
        },
        "arms": arms,
        "interpretation_guard": [
            "The captured stage current is a quadrature consistency probe; it does not by itself establish the exact continuum matter-coupled Einstein identity.",
            "The actual accepted output is compared to both stage and endpoint current quadratures; the predictor current is a labeled control only.",
            "Changing the time quadrature is not a physical source adjustment. No source term is fitted.",
            "No production source file, vendor, dependency, initializer, or scientific-admission gate is changed.",
        ],
        "acceptance": {
            "source_pins_match": True,
            "accepted_artifact_sha_verified": True,
            "accepted_input_immutable": True,
            "captured_stage_states_remain_immutable": True,
            "three_internal_geometry_stages_captured": True,
            "endpoint_vs_kernel_stage_current_quadrature_compared": True,
            "linear_control_included": True,
            "integrated_current_identity_closed": False,
            "full_matter_coupled_propagation_identity_closed": False,
            "total_stress_energy_conservation_closed": False,
            "scientific_admission": False,
        },
        "scientific_admission": False,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / "report.json"
    path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    (OUT / "SHA256SUMS.txt").write_text(f"{sha256(path)}  report.json\n", encoding="utf-8")
    print(json.dumps({
        "status": report["status"],
        "resolution": N,
        "scientific_admission": False,
        "arms": [{
            "arm": arm["arm"],
            "methods": {
                order: {
                    "endpoint": result["endpoint_trapezoid_residual"],
                    "kernel_stage": result["kernel_stage_trapezoid_residual"],
                    "predictor_control": result["predictor_stage_trapezoid_residual_control"],
                    "stage_minus_endpoint": result["stage_minus_endpoint_current_quadrature"],
                } for order, result in arm["methods"].items()
            },
            "fixed_R_cubic": arm["methods"]["cubic"]["fixed_R"],
        } for arm in arms],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
