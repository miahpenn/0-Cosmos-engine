#!/usr/bin/env python3
"""Research-only kinematic crosswalk: metric-derived R_t versus K/Aa-derived R_t.

At the same accepted state and captured TRUE-PIRK stages, compare

  R_t(metric RHS) = R*(0.5*b_dot/b - X_dot/X)
  R_t(extrinsic)  = -alpha*R*(K/3 - Aa/2)

which agree in the continuum for zero shift when the metric evolution and
extrinsic-curvature definitions are compatible. The diagnostic then rebuilds
chi and T_t^R with each rate and measures how much each rate choice moves the
same fixed-R mass/current balance. This is a crosswalk, not an authority
selection and not a production patch.
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
from tools import audit_bssn_areal_stage_current_balance as stage_balance
from tools import audit_full_profile_propagation_candidate as saved
from tools import audit_origin_cell_joint_shadow_one_step as origin_shadow
from engine.v55_matter import total_matter_projection

SOURCE = "4aec88d0e3ff028e0a65a8ac36b30948625df792"
VENDOR = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
N = int(os.environ.get("AREAL_RT_KINEMATIC_RESOLUTION", "120"))
FULL_DIR = Path(os.environ.get("FULL_PROFILE_ARTIFACT_DIR", f"reference/full/N{N}"))
OUT = Path(os.environ.get("AREAL_RT_KINEMATIC_OUT", f"results/bssn-areal-Rt-kinematic-crosswalk-N{N}-2026-10-25"))
TARGET_TIME = 0.60
CFL = 0.03
FOUR_PI = 4.0 * np.pi
FIXED_RADII = (1.0, 2.0, 3.0, 5.0, 8.0, 10.0)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def finite(name: str, values: Any) -> np.ndarray:
    result = np.asarray(values, dtype=np.float64)
    if not np.all(np.isfinite(result)):
        raise RuntimeError(f"{name} contains non-finite values")
    return result


def scalar(name: str, value: Any) -> float:
    result = float(value)
    if not np.isfinite(result):
        raise RuntimeError(f"{name} is non-finite")
    return result


def interpolate(profile: dict[str, Any], field: str, radii: np.ndarray, order: str) -> np.ndarray:
    x = finite("profile areal radii", profile["R"])
    y = finite(f"profile {field}", profile[field])
    if order == "linear":
        out = np.interp(radii, x, y)
    elif order == "cubic":
        out = np.asarray([areal.interp4(float(R), x, y) for R in radii], dtype=np.float64)
    else:
        raise ValueError(f"unknown interpolation order {order}")
    return finite(f"{order} interpolation of {field}", out)


def summary(values: np.ndarray, radii: np.ndarray) -> dict[str, float]:
    values = finite("summary profile", values)
    inner = radii <= 3.0
    annulus = (radii >= 2.0) & (radii <= 10.0)
    if not inner.any() or not annulus.any():
        raise RuntimeError("fixed-R interval does not cover required summary regions")
    return {
        "max_abs_R_le_3": float(np.max(np.abs(values[inner]))),
        "rms_R_2_to_10": float(np.sqrt(np.mean(values[annulus] ** 2))),
        "max_abs_R_2_to_10": float(np.max(np.abs(values[annulus]))),
    }


def enrich_kinematic_profile(kernel, profile: dict[str, Any]) -> dict[str, Any]:
    q = profile["state"]
    g = q.geometry
    R = finite("areal radius", profile["R"])
    Rr = finite("metric-operator Rr", profile["Rr"])
    Rt_metric = finite("metric-RHS Rt", profile["Rt"])
    alpha = finite("lapse", g.alpha)
    K = finite("trace extrinsic curvature", g.K)
    Aa = finite("trace-free extrinsic curvature Aa", g.Aa)
    Ktheta = K / 3.0 - Aa / 2.0
    Rt_extrinsic = finite("extrinsic-curvature Rt", -alpha * R * Ktheta)
    delta_Rt = finite("metric-minus-extrinsic Rt", Rt_metric - Rt_extrinsic)
    alpha_safe = np.maximum(alpha, 1e-300)
    delta_normal_rate = finite("normal-rate difference", delta_Rt / alpha_safe)

    projection = total_matter_projection(q.grid, g, q.scalars, q.matter)
    E = finite("total Eulerian energy", projection["rho"])
    j = finite("total covariant radial current", projection["j"])
    pr = finite("total radial pressure", projection["pr"])
    a = finite("a metric", g.a)
    X = finite("X metric", g.X)
    grr = a / (X * X)
    j_up = j / grr

    def current_for_Rt(Rt: np.ndarray, name: str) -> tuple[np.ndarray, np.ndarray]:
        chi = finite(f"{name} chi", (X * X / a) * Rr * Rr - (Rt / alpha) ** 2)
        T_t_R = finite(
            f"{name} transformed T_t^R",
            -E * Rt - alpha * j_up * Rr - (Rt / Rr) * (j / alpha * Rt + pr * Rr),
        )
        F = finite(f"{name} current", FOUR_PI * R * R * T_t_R)
        return chi, F

    chi_metric, F_metric = current_for_Rt(Rt_metric, "metric-RHS rate")
    chi_extrinsic, F_extrinsic = current_for_Rt(Rt_extrinsic, "K/Aa rate")
    if not np.allclose(chi_metric, profile["chi"], rtol=1e-10, atol=1e-10):
        raise RuntimeError("metric-RHS chi reconstruction failed native-grid cross-check")
    if not np.allclose(F_metric, profile["rhs"], rtol=1e-10, atol=1e-10):
        raise RuntimeError("metric-RHS current reconstruction failed native-grid cross-check")
    return {
        "profile": profile,
        "R_t_metric": Rt_metric,
        "R_t_extrinsic": Rt_extrinsic,
        "delta_R_t": delta_Rt,
        "delta_normal_rate": delta_normal_rate,
        "chi_metric": chi_metric,
        "chi_extrinsic": chi_extrinsic,
        "current_metric": F_metric,
        "current_extrinsic": F_extrinsic,
    }


def run_arm(kernel, accepted, dt: float) -> dict[str, Any]:
    accepted_digest = saved.tangent.state_digest(accepted)
    candidate, stage_states, stage_digests = stage_balance.capture_stage_states(kernel, accepted, dt)
    raw_states = [accepted, stage_states[0], stage_states[1], stage_states[2], candidate]
    labels = ["accepted_start", "kernel_start_stage", "forward_euler_predictor", "kernel_corrected_stage", "accepted_output"]
    enriched = {
        label: enrich_kinematic_profile(kernel, areal.geometry_rate_and_mass(kernel, state))
        for label, state in zip(labels, raw_states)
    }
    if saved.tangent.state_digest(accepted) != accepted_digest:
        raise RuntimeError("accepted state mutated during kinematic crosswalk")
    for label, state, digest in zip(labels[1:4], stage_states, stage_digests):
        if saved.tangent.state_digest(state) != digest:
            raise RuntimeError(f"{label} was mutated during kinematic observation")
    if max(float(np.max(np.abs(item["profile"]["state"].geometry.beta))) for item in enriched.values()) > 1e-14:
        raise RuntimeError("zero-shift precondition failed")

    lo = max(float(item["profile"]["R"][0]) for item in enriched.values())
    hi = min(float(item["profile"]["R"][-1]) for item in enriched.values())
    hi = min(hi, 10.0)
    if not hi > lo:
        raise RuntimeError("no common areal-radius range across accepted and stage slices")
    radii = np.linspace(lo, hi, max(300, int(accepted.grid.n) * 2))

    kinematic_gap_report = {}
    for label, item in enriched.items():
        p = item["profile"]
        # Native coordinate-grid rates and the geometric state are sampled on
        # the same slice; fixed-R interpolation is separately summarized below.
        rmask = p["R"] <= 3.0
        amask = (p["R"] >= 2.0) & (p["R"] <= 10.0)
        kinematic_gap_report[label] = {
            "native_R_t_gap_max_abs_R_le_3": float(np.max(np.abs(item["delta_R_t"][rmask]))),
            "native_R_t_gap_rms_R_2_to_10": float(np.sqrt(np.mean(item["delta_R_t"][amask] ** 2))),
            "native_normal_rate_gap_rms_R_2_to_10": float(np.sqrt(np.mean(item["delta_normal_rate"][amask] ** 2))),
            "native_chi_gap_rms_R_2_to_10": float(np.sqrt(np.mean((item["chi_metric"] - item["chi_extrinsic"])[amask] ** 2))),
            "native_current_gap_rms_R_2_to_10": float(np.sqrt(np.mean((item["current_metric"] - item["current_extrinsic"])[amask] ** 2))),
        }

    method_reports: dict[str, Any] = {}
    per_method_samples: dict[str, dict[str, dict[str, np.ndarray]]] = {}
    for order in ("linear", "cubic"):
        method_series = {}
        for label, item in enriched.items():
            method_series[label] = {
                "chi_metric": interpolate(item, "chi_metric", radii, order),
                "chi_extrinsic": interpolate(item, "chi_extrinsic", radii, order),
                "current_metric": interpolate(item, "current_metric", radii, order),
                "current_extrinsic": interpolate(item, "current_extrinsic", radii, order),
                "R_t_metric": interpolate(item, "R_t_metric", radii, order),
                "R_t_extrinsic": interpolate(item, "R_t_extrinsic", radii, order),
            }
        start, end = method_series["accepted_start"], method_series["accepted_output"]
        stage0, stage1 = method_series["kernel_start_stage"], method_series["kernel_corrected_stage"]
        balance = {}
        for lane, chi_key, current_key in (
            ("metric_rate", "chi_metric", "current_metric"),
            ("extrinsic_rate", "chi_extrinsic", "current_extrinsic"),
        ):
            mdot = finite(f"{order}.{lane} fixed-R mass rate",
                -0.5 * radii * (end[chi_key] - start[chi_key]) / dt)
            endpoint_average = finite(f"{order}.{lane} endpoint current average",
                0.5 * (start[current_key] + end[current_key]))
            stage_average = finite(f"{order}.{lane} stage current average",
                0.5 * (stage0[current_key] + stage1[current_key]))
            balance[lane] = {
                "endpoint_quadrature": summary(mdot - endpoint_average, radii),
                "kernel_stage_quadrature": summary(mdot - stage_average, radii),
            }
        delta_chi_endpoint = finite(f"{order}.chi endpoint path delta",
            (end["chi_extrinsic"] - start["chi_extrinsic"]) - (end["chi_metric"] - start["chi_metric"]))
        delta_mdot = finite(f"{order}.mass rate extrinsic minus metric", -0.5 * radii * delta_chi_endpoint / dt)
        delta_current_endpoint = finite(f"{order}.endpoint current delta",
            0.5 * ((start["current_extrinsic"] - start["current_metric"]) +
                   (end["current_extrinsic"] - end["current_metric"])))
        delta_current_stage = finite(f"{order}.stage current delta",
            0.5 * ((stage0["current_extrinsic"] - stage0["current_metric"]) +
                   (stage1["current_extrinsic"] - stage1["current_metric"])))
        balance["kinematic_swap_delta"] = {
            "delta_mass_rate_extrinsic_minus_metric": summary(delta_mdot, radii),
            "delta_endpoint_current_extrinsic_minus_metric": summary(delta_current_endpoint, radii),
            "delta_stage_current_extrinsic_minus_metric": summary(delta_current_stage, radii),
            "delta_endpoint_residual_extrinsic_minus_metric": summary(delta_mdot - delta_current_endpoint, radii),
            "delta_stage_residual_extrinsic_minus_metric": summary(delta_mdot - delta_current_stage, radii),
        }
        method_reports[order] = balance
        per_method_samples[order] = method_series

    fixed_R = {}
    for R in FIXED_RADII:
        if not lo <= R <= hi:
            continue
        fixed_R[str(R)] = {}
        for order, series in per_method_samples.items():
            start, end = series["accepted_start"], series["accepted_output"]
            stage0, stage1 = series["kernel_start_stage"], series["kernel_corrected_stage"]
            fixed_R[str(R)][order] = {}
            for lane, chi_key, current_key in (
                ("metric_rate", "chi_metric", "current_metric"),
                ("extrinsic_rate", "chi_extrinsic", "current_extrinsic"),
            ):
                mdot = -0.5 * R * (end[chi_key] - start[chi_key]) / dt
                endpoint_current = 0.5 * (start[current_key] + end[current_key])
                stage_current = 0.5 * (stage0[current_key] + stage1[current_key])
                fixed_R[str(R)][order][lane] = {
                    "mass_rate": scalar("fixed-R mass rate", mdot),
                    "endpoint_current_average": scalar("fixed-R endpoint current", endpoint_current),
                    "stage_current_average": scalar("fixed-R stage current", stage_current),
                    "endpoint_residual": scalar("fixed-R endpoint residual", mdot - endpoint_current),
                    "stage_residual": scalar("fixed-R stage residual", mdot - stage_current),
                }

    return {
        "dt": float(dt),
        "accepted_start_time": float(accepted.t),
        "accepted_output_time": float(candidate.t),
        "captured_stages": labels[1:4],
        "accepted_input_immutable": True,
        "captured_stages_immutable": True,
        "common_areal_radius_interval": [float(lo), float(hi)],
        "common_samples": int(radii.size),
        "kinematic_gap_by_slice": kinematic_gap_report,
        "fixed_R_rate_and_balance": method_reports,
        "fixed_R": fixed_R,
    }


def main() -> int:
    if N not in (60, 120, 240):
        raise ValueError("resolution must be 60, 120 or 240")
    pins = origin_shadow.verify_source_pins()
    if pins.get("source_commit") != SOURCE or pins.get("vendor_commit") != VENDOR:
        raise RuntimeError("source/vendor pin mismatch")
    if pins.get("protected_engine_vendor_dependency_diff_empty") is not True:
        raise RuntimeError("protected production/vendor/dependency diff is not empty")
    manifest = saved.check_artifact(FULL_DIR, N, full_profile=True)
    arrays = saved.load_profile_arrays(FULL_DIR)
    prefix = f"N{N}_full_continuous_metric_profile_final_"
    kernel, accepted = saved.reconstruct_state(N, arrays, prefix)
    origin_shadow._install_regular_row(kernel)
    if abs(float(accepted.t) - TARGET_TIME) > 1e-12:
        raise RuntimeError("accepted input is not at t=0.60")
    digest = saved.tangent.state_digest(accepted)
    dt = CFL * float(accepted.grid.dr)
    arms = [{"arm": label, **run_arm(kernel, accepted, h)} for label, h in (("dt", dt), ("dt_half", dt / 2.0))]
    if saved.tangent.state_digest(accepted) != digest:
        raise RuntimeError("accepted state changed across kinematic crosswalk")

    report = {
        "schema": "bssn-areal-Rt-kinematic-crosswalk/v1",
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
            "metric_rate": "R_t=R*(0.5*b_dot/b-X_dot/X) from the pinned explicit metric RHS",
            "extrinsic_rate": "R_t=-alpha*R*(K/3-Aa/2), from the same slice's kinematic BSSN relation",
            "current_rebuild": "recompute chi and transformed total-matter T_t^R using each R_t and otherwise identical native R_r, a/X/alpha, E/j/p_r",
            "actual_stages": "accepted start, captured start stage, Euler predictor control, corrected stage, and accepted output",
        },
        "arms": arms,
        "interpretation_guard": [
            "The metric and extrinsic-curvature rates coincide in the continuum kinematic relation for zero shift; this test measures their discrete compatibility and its impact on the fixed-R balance.",
            "A smaller residual from either rate is not a method-selection rule. The field definitions and pinned equations remain authoritative.",
            "The predictor state is a control; only accepted endpoints and the kernel's start/corrected-stage pair represent the production temporal quadrature.",
            "No production source, initializer, vendor, dependency, source normalization, or scientific-admission gate was changed.",
        ],
        "acceptance": {
            "source_pins_match": True,
            "accepted_artifact_sha_verified": True,
            "accepted_input_immutable": True,
            "captured_stages_immutable": True,
            "native_metric_Rt_reconstruction_checked": True,
            "metric_vs_K_Aa_Rt_crosswalk_completed": True,
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
            "kinematic_gap": arm["kinematic_gap_by_slice"],
            "balance": arm["fixed_R_rate_and_balance"],
        } for arm in arms],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
