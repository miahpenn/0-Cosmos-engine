#!/usr/bin/env python3
"""Evaluate existing Issue #4 propagation candidates on saved accepted states.

Research diagnostic only. Loads SHA-verified t=0.60 profiles from the prior
production-baseline/joint-origin run and the full-continuous-metric shadow.
Runs the existing RHS tangent and provisional Brown/source-divergence candidate without redefining its equation convention
on those exact states. Lapse perturbations use the production CMC solver for
the baseline and the copied regular-origin CMC row for the two shadow arms.

This does not assert the candidate formula is the accepted propagation oracle:
the matter source mapping/index normalization and complete sector conservation
identity remain provisional. No evolution is repeated and no production source
is changed.
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
# The existing sector-resolved divergence module requires a source-root variable
# at import time; actual archived source/vendor provenance is verified separately
# by verify_source_pins() and the input artifact manifests below.
os.environ.setdefault("COSMOS_PINNED_SOURCE_ROOT", str(ROOT))

from engine.production_kernel import V55ProductionKernel
from engine.scalar_system import ScalarFields
from engine.matter_system import ConservedSpecies
from engine.v55_matter import V55MatterState
from engine import v55_pirk_adapter as adapter
from tools import audit_current_source_momentum_constraint_tangent as tangent
from tools import audit_hybrid_origin_integral_patch_newton as propagation
from tools import audit_total_stress_energy_divergence as stress_divergence
from tools.audit_origin_cell_joint_shadow_one_step import (
    SOURCE_COMMIT,
    VENDOR_COMMIT,
    regular_origin_cmc_lapse,
    verify_source_pins,
)

RESOLUTION = int(os.environ.get("PROPAGATION_RESOLUTION", "120"))
JOINT_DIR = Path(os.environ.get("JOINT_ARTIFACT_DIR", f"reference/joint/N{RESOLUTION}"))
FULL_DIR = Path(os.environ.get("FULL_PROFILE_ARTIFACT_DIR", f"reference/full/N{RESOLUTION}"))
OUT = Path(os.environ.get("PROPAGATION_OUT", f"results/full-profile-propagation-N{RESOLUTION}-2026-10-25"))
EXPECTED_SOURCE = "4aec88d0e3ff028e0a65a8ac36b30948625df792"
EXPECTED_VENDOR = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
TARGET_TIME = 0.60
ARMS = (
    ("production_baseline", "production_baseline", False),
    ("partial_joint", "joint_integral_and_cmc_row", True),
    ("full_continuous_profile", "full_continuous_metric_profile", True),
)
SCALAR_NAMES = ("S", "PS", "D", "PD", "phi", "Pi")
GEOMETRY_NAMES = ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B")
MATTER_FIELDS = ("rest", "energy_t", "momentum_r")
EPSILONS = (2.0e-5, 1.0e-5)
FULL_METRIC_LAPSE_SHADOW = False
SHIFTED_STATES: dict[float, Any] = {}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def check_artifact(path: Path, resolution: int, *, full_profile: bool) -> dict[str, Any]:
    if not path.is_dir():
        raise FileNotFoundError(f"artifact directory does not exist: {path}")
    manifest_path = path / "SHA256SUMS.txt"
    if not manifest_path.is_file():
        raise FileNotFoundError(f"artifact checksum manifest missing: {manifest_path}")
    manifest = {}
    for line in manifest_path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        parts = line.split(None, 1)
        if len(parts) != 2:
            raise RuntimeError(f"malformed checksum row in {manifest_path}: {line!r}")
        expected, name = parts
        name = name.strip().lstrip("*")
        candidate = (path / name).resolve()
        if path.resolve() not in candidate.parents:
            raise RuntimeError(f"manifest path escapes artifact directory: {name}")
        if not candidate.is_file() or sha256(candidate) != expected:
            raise RuntimeError(f"artifact hash mismatch: {candidate}")
        manifest[name] = expected
    for required in ("report.json", "profiles.npz"):
        if required not in manifest:
            raise RuntimeError(f"artifact manifest omits {required}: {path}")
    report = json.loads((path / "report.json").read_text(encoding="utf-8"))
    cfg = report.get("configuration", {})
    if int(cfg.get("resolution", -1)) != resolution:
        raise RuntimeError(f"{path}: resolution mismatch in report")
    if abs(float(cfg.get("target_time", -1.0)) - TARGET_TIME) > 1e-12:
        raise RuntimeError(f"{path}: target-time mismatch in report")
    if float(cfg.get("r_max", -1.0)) != 80.0:
        raise RuntimeError(f"{path}: domain mismatch in report")
    pins = report.get("source_pins", {}).get("source_pin_check", {})
    if pins.get("source_commit") != EXPECTED_SOURCE or pins.get("vendor_commit") != EXPECTED_VENDOR:
        raise RuntimeError(f"{path}: source/vendor provenance mismatch")
    if pins.get("protected_engine_vendor_dependency_diff_empty") is not True:
        raise RuntimeError(f"{path}: protected source diff was not empty")
    if report.get("scientific_admission") is not False:
        raise RuntimeError(f"{path}: unexpectedly changed scientific_admission")
    if full_profile and cfg.get("arms") != ["full_continuous_metric_profile"]:
        raise RuntimeError(f"{path}: not a single full-profile artifact")
    return {
        "directory": str(path),
        "manifest_sha256": sha256(manifest_path),
        "manifest": manifest,
        "report_sha256": sha256(path / "report.json"),
        "profiles_sha256": sha256(path / "profiles.npz"),
        "status": report.get("status"),
        "source_commit": pins["source_commit"],
        "vendor_commit": pins["vendor_commit"],
        "configuration": cfg,
    }


def load_profile_arrays(path: Path) -> dict[str, np.ndarray]:
    with np.load(path / "profiles.npz", allow_pickle=False) as data:
        return {key: np.asarray(data[key]).copy() for key in data.files}


def scalar(data: dict[str, np.ndarray], key: str) -> float:
    return float(np.asarray(data[key]).reshape(()))


def reconstruct_state(resolution: int, data: dict[str, np.ndarray], prefix: str):
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=resolution,
        r_max=80.0,
        amplitude=0.01,
        width=7.0,
        D_amplitude=1.0e-10,
        include_radiation=True,
    )
    r_saved = np.asarray(data[prefix + "r"], dtype=np.float64)
    r_live = np.asarray(state.grid.centers, dtype=np.float64)
    if r_saved.shape != r_live.shape or not np.array_equal(r_saved, r_live):
        difference = float(np.max(np.abs(r_saved - r_live))) if r_saved.shape == r_live.shape else None
        raise RuntimeError(f"{prefix}: stored grid does not exactly match the pinned initializer grid; max_delta={difference}")

    geom = state.geometry.copy()
    for name in GEOMETRY_NAMES:
        value = np.asarray(data[prefix + name], dtype=np.float64)
        if value.shape != r_live.shape or not np.all(np.isfinite(value)):
            raise RuntimeError(f"{prefix}: invalid geometry field {name}")
        setattr(geom, name, value.copy())
    state.geometry = geom

    scalar_arrays = []
    for name in SCALAR_NAMES:
        value = np.asarray(data[prefix + name], dtype=np.float64)
        if value.shape != r_live.shape or not np.all(np.isfinite(value)):
            raise RuntimeError(f"{prefix}: invalid scalar field {name}")
        scalar_arrays.append(value.copy())
    state.scalars = ScalarFields(*scalar_arrays)

    def species(short: str) -> ConservedSpecies:
        values = []
        for name in MATTER_FIELDS:
            value = np.asarray(data[prefix + short + "_" + name], dtype=np.float64)
            if value.shape != r_live.shape or not np.all(np.isfinite(value)):
                raise RuntimeError(f"{prefix}: invalid matter field {short}.{name}")
            values.append(value.copy())
        return ConservedSpecies(*values)

    state.matter = V55MatterState(
        species("dm"), species("baryon"), species("radiation")
    )
    state.t = scalar(data, prefix + "t")
    state.tau = scalar(data, prefix + "tau")
    state.e_folds = scalar(data, prefix + "e_folds")
    if abs(state.t - TARGET_TIME) > 1e-12:
        raise RuntimeError(f"{prefix}: saved state is not accepted at t={TARGET_TIME}")
    if getattr(state, "electromagnetic", None) is not None:
        raise RuntimeError("this six-sector ledger must fail closed on an electromagnetic state")
    if getattr(state.matter, "charged_dust", None) is not None:
        raise RuntimeError("this six-sector ledger must fail closed on charged dust")
    return kernel, state


def shifted_state_with_matched_lapse(kernel, state, rates, delta):
    """Reuse the pinned tangent displacement, but match each arm's CMC row."""
    geom = state.geometry.copy()
    for name in tangent.GEOMETRY_RATE_NAMES:
        updated = (
            np.asarray(getattr(state.geometry, name), dtype=np.float64)
            + delta * np.asarray(rates["geometry"][name], dtype=np.float64)
        )
        setattr(geom, name, updated.copy())

    scalars = ScalarFields(*(
        np.asarray(getattr(state.scalars, name), dtype=np.float64)
        + delta * np.asarray(rates["scalars"][name], dtype=np.float64)
        for name in SCALAR_NAMES
    ))
    matter_parts = {}
    for species_name, short in (
        ("dark_matter", "dm"), ("baryons", "baryon"), ("radiation", "radiation")
    ):
        base = getattr(state.matter, species_name)
        rate = rates["matter"][species_name]
        matter_parts[species_name] = ConservedSpecies(*(
            np.asarray(getattr(base, name), dtype=np.float64)
            + delta * np.asarray(rate[name], dtype=np.float64)
            for name in MATTER_FIELDS
        ))
    matter = V55MatterState(
        matter_parts["dark_matter"], matter_parts["baryons"], matter_parts["radiation"]
    )
    geom.beta = np.zeros_like(geom.beta)
    geom.B = np.zeros_like(geom.B)
    if FULL_METRIC_LAPSE_SHADOW:
        alpha, _ = regular_origin_cmc_lapse(state.grid, geom, scalars, matter)
    else:
        alpha = kernel._solve_lapse(state.grid, geom, scalars, matter)[0]
    geom.alpha = np.asarray(alpha, dtype=np.float64).copy()
    geom.beta = np.zeros_like(geom.beta)
    geom.B = np.zeros_like(geom.B)
    shifted = type(state)(
        grid=state.grid,
        geometry=geom,
        scalars=scalars,
        matter=matter,
        t=state.t,
        tau=state.tau,
        e_folds=state.e_folds,
    )
    SHIFTED_STATES[float(delta)] = shifted
    return shifted


def summarize_profile(values: Any, radii: np.ndarray) -> dict[str, Any]:
    arr = np.asarray(values, dtype=np.float64)
    if arr.shape != radii.shape or not np.all(np.isfinite(arr)):
        raise RuntimeError("invalid profile supplied to residual summary")
    centre = radii <= 3.0
    annulus = (radii >= 2.0) & (radii < 10.0)
    outer = radii >= 5.0
    return {
        "cell_0_to_2": [float(x) for x in arr[:3]],
        "cell_0_to_4": [float(x) for x in arr[:5]],
        "max_abs_r_le_3": float(np.max(np.abs(arr[centre]))) if np.any(centre) else None,
        "rms_r_2_to_10": float(np.sqrt(np.mean(arr[annulus] ** 2))) if np.any(annulus) else None,
        "max_abs_r_ge_5": float(np.max(np.abs(arr[outer]))) if np.any(outer) else None,
        "rms_full_domain": float(np.sqrt(np.mean(arr * arr))),
        "fixed_radius": {
            str(q): float(np.interp(q, radii, arr))
            for q in (1.0, 2.0, 3.0, 5.0, 8.0, 10.0)
            if radii[0] <= q <= radii[-1]
        },
    }


def compact_candidate(item: dict[str, Any], raw: dict[float, dict[str, np.ndarray]], radii: np.ndarray) -> dict[str, Any]:
    eps = float(item["epsilon"])
    prop = item["independent_constraint_propagation_candidate"]
    direct = prop["published_covariant_M_r_candidate"]
    geometric = prop["geometric_constraint_propagation_candidate"]
    direct_residual = direct["measured_minus_predicted_covariant"]
    source_map = geometric["matter_source_composite_minus_stress_divergence"]
    brown_residual = geometric["measured_total_minus_Brown_plus_actual_source"]
    raw_item = raw[eps]

    return {
        "epsilon": eps,
        "direct_published_covariant_Mr_candidate": {
            "relative_rms_error_full_domain": direct["relative_rms_error_full_domain"],
            "relative_rms_error_r_ge_5": direct["relative_rms_error_r_ge_5"],
            "measured_minus_predicted": direct_residual,
        },
        "geometric_Brown_plus_actual_source_candidate": {
            "relative_rms": geometric["full_domain_relative_rms"]["total_Brown_plus_actual_matter_composite"],
            "measured_minus_predicted": brown_residual,
            "matter_source_composite_minus_minus8pi_alpha_Qr": source_map,
            "source_mapping_relative_rms": geometric["full_domain_relative_rms"]["actual_matter_composite_minus_minus8pi_alpha_Qr"],
            "nonmatter_geometric_minus_Brown_vacuum_relative_rms": geometric["full_domain_relative_rms"]["nonmatter_geometric_minus_Brown_vacuum"],
        },
        "explicit_propagation_profiles": {
            name: summarize_profile(raw_item[name], radii)
            for name in (
                "propagation_Mtotal_measured_minus_Brown_plus_actual_matter_composite",
                "propagation_matter_source_mapping_residual",
                "propagation_exact_chain_rule_residual",
                "propagation_Mtotal_decomposition_residual",
                "propagation_term_H_gradient",
                "propagation_term_connection_second",
                "propagation_term_connection_first",
                "propagation_published_Mr_term_KM",
                "propagation_published_Mr_term_alpha_H",
                "propagation_published_Mr_matter_source",
                "propagation_Qr_total_stress_divergence",
            )
            if name in raw_item
        },
    }


def compact_sector_divergence(
    state,
    epsilon: float,
    candidate_raw: dict[str, np.ndarray],
) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    """Reuse the existing additive sector-resolved total-T divergence observer."""
    try:
        plus_state = SHIFTED_STATES[float(epsilon)]
        minus_state = SHIFTED_STATES[-float(epsilon)]
    except KeyError as exc:
        raise RuntimeError(f"matched tangent endpoint not captured for epsilon={epsilon}") from exc
    witnessed = stress_divergence.sector_resolved_divergence_tangent_from_states(
        state, plus_state, minus_state, float(epsilon)
    )
    total = witnessed["total"]
    decomposition = witnessed["sector_decomposition"]
    radii = np.asarray(state.grid.centers, dtype=np.float64)
    total_r = np.asarray(total["div_T_r"], dtype=np.float64)
    total_r_fv = np.asarray(total["div_T_r_fv"], dtype=np.float64)
    candidate_qr = np.asarray(
        candidate_raw["propagation_Qr_total_stress_divergence"], dtype=np.float64
    )
    q_compare = total_r - candidate_qr
    recon = decomposition["reconstruction"]
    channels = recon["channels"]
    sector_rows = {}
    profile_store: dict[str, np.ndarray] = {
        "sector_resolved_total_div_T_r": total_r,
        "sector_resolved_total_div_T_r_fv": total_r_fv,
        "sector_resolved_total_minus_candidate_Qr": q_compare,
        "sector_resolved_sector_sum_minus_total_div_T_r": np.asarray(
            channels["div_T_r"]["sum_minus_total"], dtype=np.float64
        ),
        "sector_resolved_unexplained_gap_after_time_density_div_T_r": np.asarray(
            channels["div_T_r"]["unexplained_gap_after_time_density_accounting"],
            dtype=np.float64,
        ),
    }
    for name, data in decomposition["sectors"].items():
        raw = data["raw_signed_arrays"]
        sector_r = np.asarray(raw["div_T_r"], dtype=np.float64)
        sector_r_fv = np.asarray(raw["div_T_r_fv"], dtype=np.float64)
        profile_store[f"sector_{name}_div_T_r"] = sector_r
        profile_store[f"sector_{name}_div_T_r_fv"] = sector_r_fv
        sector_rows[name] = {
            "div_T_r": summarize_profile(sector_r, radii),
            "div_T_r_finite_volume": summarize_profile(sector_r_fv, radii),
            "cell0_cell1": [float(sector_r[0]), float(sector_r[1])],
        }
        # For scalar sectors, also retain the existing observer's signed split
        # into time-density, radial-flux, and metric-source components.
        operator = data.get("operator_component_decomposition")
        if operator:
            pointwise = operator["pointwise"]["raw_signed_arrays"]
            names = (
                "time_density_div_T_r",
                "radial_flux_div_T_r",
                "metric_source_div_T_r",
                "component_sum_div_T_r",
                "component_sum_gap_div_T_r",
            )
            sector_rows[name]["pointwise_operator_terms"] = {
                term: summarize_profile(np.asarray(pointwise[term], dtype=np.float64), radii)
                for term in names if term in pointwise
            }
            for term in names:
                if term in pointwise:
                    profile_store[f"sector_{name}_{term}"] = np.asarray(
                        pointwise[term], dtype=np.float64
                    )
    compact = {
        "epsilon": float(epsilon),
        "observer": "existing audit_total_stress_energy_divergence.sector_resolved_divergence_tangent_from_states",
        "scope": "same saved base accepted state and the exact matched-lapse plus/minus endpoint states used by the existing propagation candidate",
        "total_div_T_r": summarize_profile(total_r, radii),
        "total_div_T_r_finite_volume": summarize_profile(total_r_fv, radii),
        "independent_total_divergence_minus_candidate_Qr": summarize_profile(q_compare, radii),
        "sector_contributions": sector_rows,
        "additive_sector_reconstruction": {
            "all_sectors_reconstruct_total": bool(recon["all_sectors_reconstruct_total"]),
            "total_div_T_r_sum_passed": bool(channels["div_T_r"]["passed"]),
            "total_div_T_r_sum_max_abs_gap": float(channels["div_T_r"]["max_abs_gap"]),
            "total_div_T_r_component_corrected_passed": bool(
                channels["div_T_r"]["component_corrected_passed"]
            ),
            "total_div_T_r_unexplained_max_abs_gap": float(
                channels["div_T_r"]["max_abs_unexplained_gap"]
            ),
            "time_density_component_accounting_passed": bool(
                recon["time_density_component_accounting_passed"]
            ),
            "time_density_gap_div_T_r": summarize_profile(
                np.asarray(
                    decomposition["time_density_reconstruction"]["raw_signed_arrays"]["time_density_gap_r"],
                    dtype=np.float64,
                ),
                radii,
            ),
        },
        "sector_pairing": decomposition["state_pairing"],
        "scientific_admission": False,
    }
    return compact, profile_store


def main() -> int:
    if RESOLUTION not in (60, 120, 240):
        raise ValueError(f"unsupported resolution {RESOLUTION}")
    if SOURCE_COMMIT != EXPECTED_SOURCE or VENDOR_COMMIT != EXPECTED_VENDOR:
        raise RuntimeError("helper source/vendor constants do not match the archived runs")
    pin = verify_source_pins()
    if pin.get("source_commit") != EXPECTED_SOURCE or pin.get("vendor_commit") != EXPECTED_VENDOR:
        raise RuntimeError("live source/vendor pin verification failed")

    joint_manifest = check_artifact(JOINT_DIR, RESOLUTION, full_profile=False)
    full_manifest = check_artifact(FULL_DIR, RESOLUTION, full_profile=True)
    joint_arrays = load_profile_arrays(JOINT_DIR)
    full_arrays = load_profile_arrays(FULL_DIR)

    _, vacuum, _ = adapter.vendor_modules()
    OUT.mkdir(parents=True, exist_ok=True)
    summaries = []
    saved_profiles: dict[str, np.ndarray] = {}
    state_hashes = {}
    global FULL_METRIC_LAPSE_SHADOW

    for label, arm_name, uses_regular_lapse in ARMS:
        if uses_regular_lapse:
            prefix_data = full_arrays if arm_name == "full_continuous_metric_profile" else joint_arrays
            prefix = f"N{RESOLUTION}_{arm_name}_final_"
        else:
            prefix_data = joint_arrays
            prefix = f"N{RESOLUTION}_{arm_name}_final_"

        FULL_METRIC_LAPSE_SHADOW = bool(uses_regular_lapse)
        kernel, state = reconstruct_state(RESOLUTION, prefix_data, prefix)
        initial_digest = tangent.state_digest(state)
        tangent.shifted_state = shifted_state_with_matched_lapse
        propagation.shifted_state = shifted_state_with_matched_lapse

        SHIFTED_STATES.clear()
        report, raw = propagation.directional_tangent(kernel, state, vacuum)
        if report.get("state_unchanged") is not True:
            raise RuntimeError(f"{label}: propagation observer mutated its frozen accepted state")
        final_digest = tangent.state_digest(state)
        if initial_digest != final_digest:
            raise RuntimeError(f"{label}: accepted state digest changed during observation")
        state_hashes[label] = initial_digest

        radii = np.asarray(state.grid.centers, dtype=np.float64)
        arm_rows = {
            "arm": label,
            "source_profile_arm": arm_name,
            "lapse_treatment": "copied regular-origin CMC row" if uses_regular_lapse else "production CMC row",
            "accepted_time": float(state.t),
            "state_digest_sha256": initial_digest,
            "constraint_baseline": report["baseline_constraints"],
            "epsilon_results": [],
        }
        for item in report["tangent_by_epsilon"]:
            eps = float(item["epsilon"])
            compact = compact_candidate(item, raw, radii)
            compact["signed_momentum_budget"] = item["signed_momentum_budget"]
            compact["exact_discrete_chain_rule"] = item["independent_constraint_propagation_candidate"]["exact_discrete_chain_rule"]
            sector_report, sector_profiles = compact_sector_divergence(
                state, eps, raw[eps]
            )
            compact["sector_resolved_stress_energy_divergence"] = sector_report
            arm_rows["epsilon_results"].append(compact)
            for name, values in raw[eps].items():
                saved_profiles[f"N{RESOLUTION}_{label}_eps{eps:.0e}_{name}"] = np.asarray(values, dtype=np.float64).copy()
            for name, values in sector_profiles.items():
                saved_profiles[
                    f"N{RESOLUTION}_{label}_eps{eps:.0e}_{name}"
                ] = np.asarray(values, dtype=np.float64).copy()
        summaries.append(arm_rows)

    report_path = OUT / "report.json"
    result = {
        "schema": "full-profile-accepted-state-propagation-candidate/v2",
        "status": "completed_diagnostic_only",
        "purpose": "Apply the existing convention/index-sensitive propagation candidate and its existing additive sector-resolved total-stress-energy divergence observer to source-pinned accepted t=0.60 states from baseline, partial-joint, and full-continuous-metric arms.",
        "source_pins": {
            "source_commit": EXPECTED_SOURCE,
            "vendor_commit": EXPECTED_VENDOR,
            "live_pin_check": pin,
            "protected_engine_vendor_dependency_diff_empty": pin.get("protected_engine_vendor_dependency_diff_empty", False),
            "joint_artifact": joint_manifest,
            "full_profile_artifact": full_manifest,
        },
        "runtime": {"python": sys.version, "numpy": np.__version__, "platform": platform.platform()},
        "configuration": {
            "resolution": RESOLUTION,
            "r_max": 80.0,
            "accepted_time": TARGET_TIME,
            "arms": [x[0] for x in ARMS],
            "directional_epsilons": list(EPSILONS),
            "lapse_policy": {
                "production_baseline": "the production CMC lapse solve at perturbed endpoints",
                "partial_joint": "the copied regular-origin CMC row at perturbed endpoints",
                "full_continuous_profile": "the copied regular-origin CMC row at perturbed endpoints",
            },
        },
        "state_sha256": state_hashes,
        "arms": summaries,
        "acceptance": {
            "all_source_pins_match": True,
            "input_artifacts_sha256_verified": True,
            "all_frozen_state_observers_immutable": True,
            "all_reported_profiles_finite": True,
            "full_matter_coupled_propagation_identity_closed": False,
            "matter_source_sign_and_normalization_reviewed": False,
            "reference_metric_connection_coupling_convention_frozen": False,
            "total_stress_energy_sector_conservation_closed": False,
            "production_source_or_vendor_changed": False,
            "scientific_admission": False,
        },
        "interpretation_limits": [
            "This evaluates the existing Brown Eq. (11b)-based candidates, not a newly asserted acceptance oracle.",
            "The direct published-covariant-M_r candidate uses matter-inclusive H/M/C with a provisional -8*pi*alpha*Q_r source.",
            "The geometric-vacuum comparator uses H_geom/M_geom and separately compares the actual matter source composite with the independently reconstructed stress-divergence source.",
            "The source-mapping, scalar spatial-operator, variable-mass exchange, and full sector conservation identities remain incomplete/unreviewed.",
            "Sector-resolved div_T_r is an existing component-accounting observer on the same directional endpoint states; its reconstruction result is measured, not treated as an acceptance threshold.",
            "Centre cells 0/1 are reported separately from fixed physical radii and annuli; no cell-index comparison is presented as a matched-radius comparison.",
            "No evolution is advanced and no initializer, engine, vendor, dependency, or production equation is modified.",
        ],
    }
    report_path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    profile_path = OUT / "profiles.npz"
    np.savez_compressed(profile_path, **saved_profiles)
    (OUT / "SHA256SUMS.txt").write_text(
        f"{sha256(report_path)}  report.json\n{sha256(profile_path)}  profiles.npz\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "status": result["status"],
        "resolution": RESOLUTION,
        "output_dir": str(OUT),
        "state_sha256": state_hashes,
        "direct_Mr_candidate_rel_rms": {
            row["arm"]: [
                item["direct_published_covariant_Mr_candidate"]["relative_rms_error_full_domain"]
                for item in row["epsilon_results"]
            ] for row in summaries
        },
        "matter_source_mapping_rel_rms": {
            row["arm"]: [
                item["geometric_Brown_plus_actual_source_candidate"]["source_mapping_relative_rms"]
                for item in row["epsilon_results"]
            ] for row in summaries
        },
        "scientific_admission": False,
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
