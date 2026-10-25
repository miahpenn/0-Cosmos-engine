#!/usr/bin/env python3
"""Assemble a fail-closed full matter-coupled propagation candidate report.

This verifier reuses the already-recorded #420 accepted states and tangent
profiles. It does not import or execute the engine, RHS tangent, evolution
kernel, or initializer. All source-law and scientific-admission flags remain
false; this is a provenance-checked report assembler, not a new measurement.
"""
from __future__ import annotations

import hashlib
import json
import math
import sys
import traceback
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
STATE_DIR = ROOT / "reference" / "full"
CANDIDATE_DIR = ROOT / "reference" / "candidate"
OUT = ROOT / "results" / "full-matter-coupled-propagation-candidate"
RESOLUTIONS = (60, 120, 240)
EPS_FINE = 1.0e-5
EPS_COARSE = 2.0e-5
SOURCE_COMMIT = "4aec88d0e3ff028e0a65a8ac36b30948625df792"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
PR421_HEAD = "6dfc53c48aaecc18b8389a75b6de7b14340747d1"
PR422_HEAD = "51e03c820fe4a03243a4a788d9feb039742a3664"
ACCEPTED_STATE_RUN_ID = 40156226086
CANDIDATE_OUTPUT_RUN_ID = 40156795380
BRIDGE_SUMMARY_BLOB_SHA = "2aa10c7f659d0dfcceb5868df360a8d6eb387678"
CONVENTION_SUMMARY_BLOB_SHA = "7d4af8a201859784a08c27c7e9250f354b46ef4f"
SECTOR_SUMMARY_BLOB_SHA = "8c4aa86668837e43715c7e3c8bc7f47e7a2cf01d"

# Outer ZIP digests from the source/output workflow runs. The workflow queries
# GitHub's Actions API before download and this script verifies that attestation.
ARTIFACTS = {
    60: {
        "accepted_state": {
            "id": 12584064724,
            "name": "origin-full-continuous-metric-evolved-N60-t0.6-2026-10-24",
            "zip_sha256": "9730fa5ae4977050c24e43299ecbe4f2c9636add1a6d24912a38628d2be840d7",
            "report_sha256": "5c00b065a172f2932decba11d17f2c97b3adb0c9b1c8656441122d5fae30c303",
            "profiles_sha256": "6e54865df4fe7c184924826758fe580341a32114867b54a8fd002471f260579c",
            "manifest_sha256": "7f19ace7e0e0089a7d368109217804705a7e726429074c7ce7a64dbdc5e41284",
            "state_digest_sha256": "20c07ec4b7885b59a1e0476d5bab448850518dbc8439ad6db4adaf36fd5a4a33",
        },
        "candidate_output": {
            "id": 12585305634,
            "name": "full-profile-propagation-candidate-N60-2026-10-25",
            "zip_sha256": "474f63f53c3b9337b8e700e6b83a9fc331bcddcdef9c69db374e0b1d85dd5f88",
            "report_sha256": "25a410d6829ee59652a4418f630e7b9de273b4865759be8133fc2c93f6419120",
            "profiles_sha256": "acfbe8caf2072c3b30a79780fd13aced4c12a1e9774782eb1313acac6e98cf03",
            "manifest_sha256": "dff12601bcda06ffcd61c52cb4b9c9656c621dbfb23641614606e38c8654be92",
        },
    },
    120: {
        "accepted_state": {
            "id": 12584164579,
            "name": "origin-full-continuous-metric-evolved-N120-t0.6-2026-10-24",
            "zip_sha256": "ebaebb8ee8cead10d4f49bbdd6970eac15d869ed22f251be48afd9ffece30b01",
            "report_sha256": "e77eee07e5e5f4f851bee608de551ac17ad91fdc260a4751392bff0e64b9ff55",
            "profiles_sha256": "b31cbe3534fad5e8303e284c716a3d37d4ef5f08519805966da9fd3b24d09d53",
            "manifest_sha256": "dcf5611cabe7d924b2eb0e4585e8574fb3b9bc67f3c3b7bf73397fa470d8c758",
            "state_digest_sha256": "032a5d079b76dfe7bf50b08765f03b9f2dcc730ba1a33329517c153a27bf70cb",
        },
        "candidate_output": {
            "id": 12585325663,
            "name": "full-profile-propagation-candidate-N120-2026-10-25",
            "zip_sha256": "39c0afb5829e9bd8b39830dcfa13c777fb5f616535b2e89a932fb4321ee26e39",
            "report_sha256": "526b06b2eec94eac28bf92eac403b25d441c24b0407fe7b101a833ec207838b3",
            "profiles_sha256": "754fe1bc4d77a8061d97fad1f35bbbf5e0f9e74d8bb1802e519a684f4a33a445",
            "manifest_sha256": "94bd0f77c301207065c8a3fd4eb12d63267470a3d015d3996718b462bb43df34",
        },
    },
    240: {
        "accepted_state": {
            "id": 12584980156,
            "name": "origin-full-continuous-metric-evolved-N240-t0.6-2026-10-24",
            "zip_sha256": "72a043b778cbb9efe6e433565b34c5b92f9a4ab306e35a8a54855025778d6a8f",
            "report_sha256": "bd3ac5b61bf6c6cbc3b074a455c36885a91613d74481d31e6c629a14f7414b03",
            "profiles_sha256": "4a6d30ffeb0beb6bfb37879666e01a9985e51e75e122c764e1aac388b2931ac7",
            "manifest_sha256": "04065ca99f9fd6cdeb0a638276d3aa5fea31ba041e68c587451d3c496885ecab",
            "state_digest_sha256": "dbedf4382bfa58179f87c0d91022574c87f92865c10d125ee66a36181234b5f1",
        },
        "candidate_output": {
            "id": 12585585288,
            "name": "full-profile-propagation-candidate-N240-2026-10-25",
            "zip_sha256": "ddf9031ee09c4aad39561a11bc22580320a939660e344336cbd1fc8e2dbb872b",
            "report_sha256": "bb0b8b841432c60c9388cfec095abb3afd211592f9f727b25c04dce1141d5bf5",
            "profiles_sha256": "462b2899f6287f1f9adfd5ad7f1ac28e332dbc68e3508bf1b924351d679e9334",
            "manifest_sha256": "bba0b1acd8f2f5ca3c13c544aa38caced3b618552b99c943f6dd10475c277660",
        },
    },
}

FLAGS_FALSE = {
    "full_matter_coupled_propagation_identity_closed": False,
    "matter_source_sign_and_normalization_reviewed": False,
    "reference_metric_connection_coupling_convention_frozen": False,
    "total_stress_energy_sector_conservation_closed": False,
    "scientific_admission": False,
}

TERM_SUFFIXES = (
    "propagation_published_Mr_term_shift_transport",
    "propagation_published_Mr_term_KM",
    "propagation_published_Mr_term_alpha_H",
    "propagation_published_Mr_term_H_gradient",
    "propagation_published_Mr_term_connection_second",
    "propagation_published_Mr_term_connection_first",
    "propagation_published_Mr_matter_source",
)
SECTOR_SUFFIXES = {
    "S": "Mdot_sector_S",
    "D": "Mdot_sector_D",
    "phi": "Mdot_sector_COSMOS_phi",
    "DM": "Mdot_sector_dark_matter",
    "baryon": "Mdot_sector_baryons",
    "radiation": "Mdot_sector_radiation",
}


def fail(message: str) -> None:
    raise RuntimeError(message)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_manifest(folder: Path, expected: dict[str, str], expected_manifest_sha: str, label: str) -> dict:
    manifest_path = folder / "SHA256SUMS.txt"
    if not manifest_path.is_file():
        fail(f"{label}: missing SHA256SUMS.txt")
    manifest_sha = sha256_file(manifest_path)
    if manifest_sha != expected_manifest_sha:
        fail(f"{label}: manifest file SHA mismatch: {manifest_sha}")
    parsed: dict[str, str] = {}
    for line in manifest_path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        parts = line.split(None, 1)
        if len(parts) != 2:
            fail(f"{label}: malformed manifest row: {line!r}")
        checksum, filename = parts
        filename = filename.strip().lstrip("*")
        if filename in parsed:
            fail(f"{label}: duplicate manifest entry {filename}")
        if filename not in {"report.json", "profiles.npz"}:
            fail(f"{label}: unexpected manifest entry {filename}")
        parsed[filename] = checksum
        path = (folder / filename).resolve()
        if folder.resolve() not in path.parents or not path.is_file():
            fail(f"{label}: manifest file missing/outside folder: {filename}")
        actual = sha256_file(path)
        if actual != checksum:
            fail(f"{label}: payload digest mismatch for {filename}: {actual} != {checksum}")
        if expected.get(filename) != actual:
            fail(f"{label}: pinned payload digest mismatch for {filename}: {actual}")
    if set(parsed) != {"report.json", "profiles.npz"}:
        fail(f"{label}: manifest must contain exactly report.json and profiles.npz")
    return {
        "manifest_sha256": manifest_sha,
        "report_sha256": parsed["report.json"],
        "profiles_sha256": parsed["profiles.npz"],
        "manifest_entries": parsed,
    }


def load_json(path: Path, label: str) -> dict:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        fail(f"{label}: invalid JSON: {exc}")
    if not isinstance(payload, dict):
        fail(f"{label}: JSON root must be an object")
    return payload


def verify_zip_metadata_file() -> dict:
    path = ROOT / "reference" / "artifact_api_verification.json"
    if not path.is_file():
        fail("missing reference/artifact_api_verification.json from workflow ZIP-digest guard")
    payload = load_json(path, "artifact API verification")
    if payload.get("status") != "verified" or payload.get("verified") is not True:
        fail("GitHub artifact API ZIP-digest verification is not marked verified")
    if int(payload.get("accepted_state_run_id", -1)) != ACCEPTED_STATE_RUN_ID:
        fail("artifact API verification has wrong accepted-state run ID")
    if int(payload.get("candidate_output_run_id", -1)) != CANDIDATE_OUTPUT_RUN_ID:
        fail("artifact API verification has wrong candidate-output run ID")
    rows = payload.get("artifacts")
    if not isinstance(rows, list):
        fail("artifact API verification lacks artifact list")
    seen = {(int(x.get("resolution", -1)), x.get("kind")): x for x in rows}
    if len(seen) != 6:
        fail(f"expected six API-verified artifacts; found {len(seen)} unique rows")
    for n in RESOLUTIONS:
        for kind in ("accepted_state", "candidate_output"):
            row = seen.get((n, kind))
            pin = ARTIFACTS[n][kind]
            expected_run_id = ACCEPTED_STATE_RUN_ID if kind == "accepted_state" else CANDIDATE_OUTPUT_RUN_ID
            if (row is None or int(row.get("artifact_id", -1)) != pin["id"]
                    or row.get("artifact_name") != pin["name"]
                    or int(row.get("workflow_run_id", -1)) != expected_run_id
                    or row.get("zip_sha256") != pin["zip_sha256"]
                    or row.get("api_digest_match") is not True):
                fail(f"GitHub artifact ZIP digest/identity verification mismatch: N={n} {kind}")
    return payload


def verify_state_report(report: dict, n: int) -> None:
    if report.get("schema") != "origin-cell-joint-shadow-evolved/v1":
        fail(f"N{n}: unexpected accepted-state report schema")
    if report.get("status") != "completed_diagnostic_only":
        fail(f"N{n}: accepted-state artifact is not completed diagnostic-only")
    cfg = report.get("configuration", {})
    if int(cfg.get("resolution", -1)) != n or float(cfg.get("target_time", -1)) != 0.6 or float(cfg.get("r_max", -1)) != 80.0:
        fail(f"N{n}: accepted-state configuration mismatch")
    if cfg.get("arms") != ["full_continuous_metric_profile"]:
        fail(f"N{n}: source artifact is not exclusively the full continuous-profile arm")
    pins = report.get("source_pins", {}).get("source_pin_check", {})
    if pins.get("source_commit") != SOURCE_COMMIT or pins.get("vendor_commit") != VENDOR_COMMIT or pins.get("vendor_checkout") != VENDOR_COMMIT:
        fail(f"N{n}: accepted-state source/vendor pins mismatch")
    if pins.get("protected_engine_vendor_dependency_diff_empty") is not True:
        fail(f"N{n}: protected source diff not verified empty")
    acceptance = report.get("acceptance", {})
    if acceptance.get("production_source_or_vendor_changed") is not False or acceptance.get("scientific_admission") is not False:
        fail(f"N{n}: accepted-state artifact has unexpected acceptance flag")


def verify_candidate_report(report: dict, n: int, state_manifest: dict) -> None:
    if report.get("schema") != "full-profile-accepted-state-propagation-candidate/v1" or report.get("status") != "completed_diagnostic_only":
        fail(f"N{n}: candidate output schema/status mismatch")
    cfg = report.get("configuration", {})
    if int(cfg.get("resolution", -1)) != n or float(cfg.get("accepted_time", -1)) != 0.6 or float(cfg.get("r_max", -1)) != 80.0:
        fail(f"N{n}: candidate output configuration mismatch")
    if "full_continuous_profile" not in cfg.get("arms", []):
        fail(f"N{n}: full continuous-profile arm missing")
    source_pins = report.get("source_pins", {})
    if source_pins.get("source_commit") != SOURCE_COMMIT or source_pins.get("vendor_commit") != VENDOR_COMMIT:
        fail(f"N{n}: candidate source/vendor pins mismatch")
    if source_pins.get("protected_engine_vendor_dependency_diff_empty") is not True:
        fail(f"N{n}: candidate protected-file diff not verified empty")
    full = source_pins.get("full_profile_artifact", {})
    expected = ARTIFACTS[n]["accepted_state"]
    if full.get("manifest_sha256") != expected["manifest_sha256"] or full.get("report_sha256") != expected["report_sha256"] or full.get("profiles_sha256") != expected["profiles_sha256"]:
        fail(f"N{n}: #420 candidate's accepted-state artifact cross-reference does not match pinned payloads")
    if full.get("manifest_sha256") != state_manifest["manifest_sha256"] or full.get("report_sha256") != state_manifest["report_sha256"] or full.get("profiles_sha256") != state_manifest["profiles_sha256"]:
        fail(f"N{n}: candidate source bundle differs from downloaded verified accepted-state bundle")
    acceptance = report.get("acceptance", {})
    for key in FLAGS_FALSE:
        if acceptance.get(key) is not False:
            fail(f"N{n}: candidate report must keep {key}=false")
    if acceptance.get("input_artifacts_sha256_verified") is not True or acceptance.get("all_frozen_state_observers_immutable") is not True:
        fail(f"N{n}: candidate report does not attest input hashes and immutability")


def canonical_state_digest(data: dict[str, np.ndarray], n: int) -> str:
    prefix = f"N{n}_full_continuous_metric_profile_final_"
    h = hashlib.sha256()

    def add(name: str, values: Any) -> None:
        arr = np.ascontiguousarray(np.asarray(values, dtype=np.float64))
        if not np.all(np.isfinite(arr)):
            fail(f"N{n}: nonfinite array while reconstructing state digest {name}")
        h.update(name.encode("utf-8"))
        h.update(str(arr.shape).encode("ascii"))
        h.update(arr.tobytes(order="C"))

    add("grid.centers", data[prefix + "r"])
    for name in ("a", "b", "X", "Aa", "K", "Lambda", "alpha", "beta", "B"):
        add(f"geometry.{name}", data[prefix + name])
    for name in ("S", "PS", "D", "PD", "phi", "Pi"):
        add(f"scalars.{name}", data[prefix + name])
    for species, name_map in (
        ("dm", "dark_matter"), ("baryon", "baryons"), ("radiation", "radiation")
    ):
        for field in ("rest", "energy_t", "momentum_r"):
            add(f"matter.{name_map}.{field}", data[prefix + species + "_" + field])
    return h.hexdigest()


def load_npz(path: Path, label: str) -> dict[str, np.ndarray]:
    try:
        with np.load(path, allow_pickle=False) as archive:
            result = {key: np.asarray(archive[key]).copy() for key in archive.files}
    except Exception as exc:
        fail(f"{label}: cannot load profile archive: {exc}")
    for key, arr in result.items():
        if not np.all(np.isfinite(arr)):
            fail(f"{label}: non-finite array {key}")
    return result


def profile_key(n: int, eps: float, suffix: str) -> str:
    eps_tag = format(eps, ".0e").replace("+", "")
    return f"N{n}_full_continuous_profile_eps{eps_tag}_{suffix}"


def profile(data: dict[str, np.ndarray], n: int, eps: float, suffix: str) -> np.ndarray:
    key = profile_key(n, eps, suffix)
    if key not in data:
        fail(f"N{n}: candidate output missing required profile {key}")
    arr = np.asarray(data[key], dtype=np.float64)
    if arr.shape != (n,) or not np.all(np.isfinite(arr)):
        fail(f"N{n}: profile {key} has invalid shape/finiteness {arr.shape}")
    return arr


def scalar_summary(values: np.ndarray, radii: np.ndarray) -> dict:
    arr = np.asarray(values, dtype=np.float64)
    if arr.shape != radii.shape or arr.size == 0 or not np.all(np.isfinite(arr)):
        fail("cannot summarize malformed/nonfinite profile")
    fixed: dict[str, dict[str, float]] = {}
    for q in (1.0, 2.0, 5.0, 8.0):
        if q < radii[0] or q > radii[-1]:
            fail(f"fixed radius r={q:g} lies outside N profile domain")
        value = float(np.interp(q, radii, arr))
        fixed[str(int(q))] = {"signed": value, "absolute": abs(value)}
    return {
        "global_max_abs": float(np.max(np.abs(arr))),
        "centre_cell0": {"signed": float(arr[0]), "absolute": abs(float(arr[0])), "radius": float(radii[0])},
        "centre_cell1": {"signed": float(arr[1]), "absolute": abs(float(arr[1])), "radius": float(radii[1])},
        "fixed_r": fixed,
        "interpretation": "Pointwise cell-centred profile; fixed-r values are linearly interpolated from this saved profile.",
    }


def decomposition_summary(values: np.ndarray, radii: np.ndarray) -> dict:
    out = scalar_summary(values, radii)
    return {k: out[k] for k in ("global_max_abs", "centre_cell0", "centre_cell1", "fixed_r")}


def observed_order(errors_coarse: float, errors_fine: float) -> float:
    if not (math.isfinite(errors_coarse) and math.isfinite(errors_fine) and errors_coarse > 0.0 and errors_fine > 0.0):
        fail("observed order requires finite positive absolute residuals at both resolutions")
    result = math.log(errors_coarse / errors_fine) / math.log(2.0)
    if not math.isfinite(result):
        fail("observed order is nonfinite")
    return result


def build_report() -> dict:
    zip_api = verify_zip_metadata_file()
    radii_by_n: dict[int, np.ndarray] = {}
    fine_residuals: dict[int, np.ndarray] = {}
    epsilon_sensitivity: dict[int, dict] = {}
    residual_summary: dict[int, dict] = {}
    decompositions: dict[int, dict] = {}
    source_mapping: dict[int, dict] = {}
    sectors: dict[str, dict] = {name: {"per_resolution": {}} for name in SECTOR_SUFFIXES}
    artifact_inventory: dict[str, dict] = {}

    for n in RESOLUTIONS:
        state_folder = STATE_DIR / f"N{n}"
        candidate_folder = CANDIDATE_DIR / f"N{n}"
        expected_state = ARTIFACTS[n]["accepted_state"]
        expected_candidate = ARTIFACTS[n]["candidate_output"]
        state_manifest = verify_manifest(
            state_folder,
            {"report.json": expected_state["report_sha256"], "profiles.npz": expected_state["profiles_sha256"]},
            expected_state["manifest_sha256"], f"N{n} accepted state",
        )
        candidate_manifest = verify_manifest(
            candidate_folder,
            {"report.json": expected_candidate["report_sha256"], "profiles.npz": expected_candidate["profiles_sha256"]},
            expected_candidate["manifest_sha256"], f"N{n} existing candidate output",
        )
        sr = load_json(state_folder / "report.json", f"N{n} accepted state report")
        cr = load_json(candidate_folder / "report.json", f"N{n} #420 candidate report")
        verify_state_report(sr, n)
        verify_candidate_report(cr, n, state_manifest)
        sdata = load_npz(state_folder / "profiles.npz", f"N{n} accepted state")
        cdata = load_npz(candidate_folder / "profiles.npz", f"N{n} existing candidate")
        digest = canonical_state_digest(sdata, n)
        if digest != expected_state["state_digest_sha256"]:
            fail(f"N{n}: accepted-state canonical digest does not match pinned expected digest")
        if cr.get("state_sha256", {}).get("full_continuous_profile") != digest:
            fail(f"N{n}: #420 candidate did not measure this exact accepted state")

        prefix = f"N{n}_full_continuous_metric_profile_final_"
        if prefix + "r" not in sdata:
            fail(f"N{n}: accepted state is missing grid centers")
        radii = np.asarray(sdata[prefix + "r"], dtype=np.float64)
        if radii.shape != (n,) or not np.all(np.isfinite(radii)) or np.any(np.diff(radii) <= 0):
            fail(f"N{n}: invalid accepted-state radius profile")
        candidate_r_key = f"N{n}_full_continuous_profile_r"
        if candidate_r_key in cdata:
            candidate_r = cdata[candidate_r_key]
            if candidate_r.shape != radii.shape or not np.array_equal(candidate_r, radii):
                fail(f"N{n}: existing candidate profile grid does not match accepted-state grid")

        eps_residuals = {}
        for eps in (EPS_FINE, EPS_COARSE):
            measured = profile(cdata, n, eps, "Mdot")
            residual = profile(cdata, n, eps, "propagation_Mdot_measured_minus_predicted_covariant_published_Mr")
            terms = [profile(cdata, n, eps, suffix) for suffix in TERM_SUFFIXES]
            predicted = np.sum(np.stack(terms, axis=0), axis=0)
            reconstructed = measured - predicted
            reconstruction_closure = float(np.max(np.abs(reconstructed - residual)))
            if reconstruction_closure > 5.0e-13:
                fail(f"N{n} eps={eps:g}: candidate residual does not reconstruct from pinned Brown slots/source; closure={reconstruction_closure:.3e}")
            eps_residuals[eps] = residual
        residual = eps_residuals[EPS_FINE]
        fine_residuals[n] = residual
        epsilon_difference = float(np.max(np.abs(eps_residuals[EPS_COARSE] - residual)))
        measured_fine = profile(cdata, n, EPS_FINE, "Mdot")
        epsilon_bound = max(1.0e-10, 1.0e-4 * float(np.max(np.abs(measured_fine))))
        if epsilon_difference > epsilon_bound:
            fail(f"N{n}: directional epsilon sensitivity exceeds absolute/signal-scaled bound: {epsilon_difference:.3e} > {epsilon_bound:.3e}")
        epsilon_sensitivity[n] = {
            "max_abs_residual_profile_difference_eps2e-5_vs_eps1e-5": epsilon_difference,
            "maximum_allowed_difference": epsilon_bound,
            "passed": True,
            "note": "Secondary tangent-step sensitivity check on already-saved #420 outputs; no tangent was recomputed here.",
        }

        sq = profile(cdata, n, EPS_FINE, "propagation_published_Mr_matter_source")
        s_act = profile(cdata, n, EPS_FINE, "propagation_actual_matter_source_composite")
        delta_s = profile(cdata, n, EPS_FINE, "propagation_matter_source_mapping_residual")
        source_closure = float(np.max(np.abs((s_act - sq) - delta_s)))
        if source_closure > 5.0e-13:
            fail(f"N{n}: S_act - S_Q does not reconstruct Delta_S; closure={source_closure:.3e}")

        shift = profile(cdata, n, EPS_FINE, "propagation_published_Mr_term_shift_transport")
        km = profile(cdata, n, EPS_FINE, "propagation_published_Mr_term_KM")
        alpha_h = profile(cdata, n, EPS_FINE, "propagation_published_Mr_term_alpha_H")
        grad_h = profile(cdata, n, EPS_FINE, "propagation_published_Mr_term_H_gradient")
        conn2 = profile(cdata, n, EPS_FINE, "propagation_published_Mr_term_connection_second")
        conn1 = profile(cdata, n, EPS_FINE, "propagation_published_Mr_term_connection_first")
        pred = shift + km + alpha_h + grad_h + conn2 + conn1 + sq
        closure = float(np.max(np.abs((profile(cdata, n, EPS_FINE, "Mdot") - pred) - residual)))
        decompositions[n] = {
            "geometric_part": decomposition_summary(shift + km, radii),
            "hamiltonian_gradient_part": decomposition_summary(alpha_h + grad_h, radii),
            "connection_part": decomposition_summary(conn2 + conn1, radii),
            "matter_induced_part_provisional_S_Q": decomposition_summary(sq, radii),
            "actual_matter_source_composite_S_act_comparator": decomposition_summary(s_act, radii),
            "Brown_slots": {
                "shift_transport": decomposition_summary(shift, radii),
                "alpha_K_M_r": decomposition_summary(km, radii),
                "minus_Dr_alpha_H_over_3": decomposition_summary(alpha_h, radii),
                "alpha_Dr_H_over_6": decomposition_summary(grad_h, radii),
                "connection_second_derivative": decomposition_summary(conn2, radii),
                "connection_first_derivative": decomposition_summary(conn1, radii),
                "S_Q_provisional": decomposition_summary(sq, radii),
            },
            "residual_reconstruction_closure_max_abs": closure,
        }
        bridge_residual = {60: 3.6248209826523353e-10, 120: 3.0087366490383e-11, 240: 2.0149287715208984e-12}[n]
        bridge_algebra_closure = {60: 5.082197683525802e-21, 120: 7.4640066396056e-21, 240: 1.062361771797708e-20}[n]
        source_mapping[n] = {
            "S_Q_provisional": decomposition_summary(sq, radii),
            "S_act": decomposition_summary(s_act, radii),
            "Delta_S": decomposition_summary(delta_s, radii),
            "S_act_minus_S_Q_reconstruction_closure_max_abs": source_closure,
            "brown_matter_bridge_residual": {
                "bridge_residual_max_abs_from_existing_421_bridge": bridge_residual,
                "algebra_relation_closure_max_abs_from_existing_421_bridge": bridge_algebra_closure,
                "interpretation": "Reused #421/#420 signed bridge and algebraic reconciliation; not an independent derivation or a source-law admission.",
            },
        }

        for sector_name, suffix in SECTOR_SUFFIXES.items():
            sector_profile = profile(cdata, n, EPS_FINE, suffix)
            sectors[sector_name]["per_resolution"][f"N{n}"] = decomposition_summary(sector_profile, radii)

        residual_summary[n] = scalar_summary(residual, radii)
        radii_by_n[n] = radii
        artifact_inventory[f"N{n}"] = {
            "accepted_state": {**expected_state, "internal_manifest": state_manifest, "state_report_status": sr["status"]},
            "existing_candidate_output": {**expected_candidate, "internal_manifest": candidate_manifest, "candidate_report_status": cr["status"]},
            "canonical_state_digest_sha256": digest,
            "candidate_state_digest_matches_accepted_state": True,
        }

    signed_fixed: dict[str, dict[str, float]] = {}
    abs_fixed: dict[str, dict[str, float]] = {}
    for q in (1.0, 2.0, 5.0, 8.0):
        signed_fixed[str(int(q))] = {f"N{n}": float(np.interp(q, radii_by_n[n], fine_residuals[n])) for n in RESOLUTIONS}
        abs_fixed[str(int(q))] = {f"N{n}": abs(signed_fixed[str(int(q))][f"N{n}"]) for n in RESOLUTIONS}
    orders = {
        "r=2": [observed_order(abs_fixed["2"]["N60"], abs_fixed["2"]["N120"]), observed_order(abs_fixed["2"]["N120"], abs_fixed["2"]["N240"])],
        "r=8": [observed_order(abs_fixed["8"]["N60"], abs_fixed["8"]["N120"]), observed_order(abs_fixed["8"]["N120"], abs_fixed["8"]["N240"])],
        "basis": "log2 ratio of absolute residual magnitude at a fixed physical radius; fixed-r values are interpolated from each stored grid",
    }

    return {
        "schema": "full-matter-coupled-propagation-candidate/v1",
        "status": "diagnostic_only",
        "mode": "reuse_verified_existing_measurements_no_rhs_recomputation",
        "created_utc_date": "2026-10-25",
        "conventions_source": [
            "docs/research/full_matter_coupled_propagation_convention_map_2026-10-25.md",
            "docs/research/matter_coupled_propagation_derivation_audit_2026-10-25.md",
        ],
        "pins": {
            "accepted_state_source_commit": SOURCE_COMMIT,
            "vendor_commit": VENDOR_COMMIT,
            "pr_421_convention_map_head": PR421_HEAD,
            "pr_422_derivation_audit_head": PR422_HEAD,
            "accepted_state_run_id": ACCEPTED_STATE_RUN_ID,
            "candidate_output_run_id": CANDIDATE_OUTPUT_RUN_ID,
            "artifact_zip_digests_api_verified": True,
            "artifact_api_verification": zip_api,
            "artifacts": artifact_inventory,
            "reused_summary_blob_shas": {
                "brown_matter_source_bridge": BRIDGE_SUMMARY_BLOB_SHA,
                "convention_map": CONVENTION_SUMMARY_BLOB_SHA,
                "sector_source_mapping": SECTOR_SUMMARY_BLOB_SHA,
            },
            "interpretation": "State-producing source pin and vendor gitlink were verified from each accepted-state report and are cross-linked to the exact state digest in #420. The legacy Brown vacuum observer's different SOURCE_COMMIT is not invoked by this verifier and is not treated as the state-producing source pin.",
        },
        "measurement_method": {
            "directional_tangent_recomputed": False,
            "evolution_advanced": False,
            "production_initializer_called": False,
            "primary_epsilon": EPS_FINE,
            "secondary_epsilon_sensitivity_check": EPS_COARSE,
            "epsilon_sensitivity": {f"N{n}": epsilon_sensitivity[n] for n in RESOLUTIONS},
            "fixed_radii": [1.0, 2.0, 5.0, 8.0],
            "accepted_time": 0.6,
            "resolutions": list(RESOLUTIONS),
        },
        "signed_residual": {
            "centre_cell0": {f"N{n}": residual_summary[n]["centre_cell0"]["signed"] for n in RESOLUTIONS},
            "centre_cell1": {f"N{n}": residual_summary[n]["centre_cell1"]["signed"] for n in RESOLUTIONS},
            "fixed_r": signed_fixed,
        },
        "absolute_residual": {
            "centre_cell0": {f"N{n}": residual_summary[n]["centre_cell0"]["absolute"] for n in RESOLUTIONS},
            "centre_cell1": {f"N{n}": residual_summary[n]["centre_cell1"]["absolute"] for n in RESOLUTIONS},
            "fixed_r": abs_fixed,
            "global_max": {f"N{n}": residual_summary[n]["global_max_abs"] for n in RESOLUTIONS},
        },
        "observed_order": orders,
        "operator_decomposition": {
            "per_resolution": {f"N{n}": decompositions[n] for n in RESOLUTIONS},
            "definition_note": "All profiles reuse the existing #420 direct covariant-M_r Brown candidate arrays at epsilon=1e-5. Geometric part is the sum of shift-transport and alpha*K*M_r; Hamiltonian-gradient part combines the two frozen Brown slots; connection part combines C^r derivative terms; matter-induced part is the provisional S_Q only. S_act is shown separately as a comparator. Closure is algebraic reconstruction, not proof of the continuum source law.",
        },
        "source_mapping": {"per_resolution": {f"N{n}": source_mapping[n] for n in RESOLUTIONS}},
        "sector_ledger": {
            name: {
                **sector,
                "interpretation": "Signed sector contribution from the existing #420 momentum-tangent ledger. Numerical additivity/bookkeeping does not establish total stress-energy conservation.",
            }
            for name, sector in sectors.items()
        },
        "provisional_conventions": [
            "S_Q=-8*pi*alpha*Q_r is a candidate source only; the verified arrays do not promote its sign, index placement, density weight, or lapse factor.",
            "Brown Eq. (11b) slots remain a vacuum structural comparator, not a proven complete matter-coupled evolution identity.",
            "The standard-ADM Hamiltonian-gradient pair is not substituted for the frozen Brown slots.",
            "Pointwise covariant, finite-volume/cell-integrated, and expanded total-stress divergences remain distinct discrete observables.",
            "The six-sector ledger reports signed tangent contributions; total stress-energy conservation and the COSMOS-phi normalization premise remain open.",
            "The accepted-state source/vendor pin and the legacy Brown observer's SOURCE_COMMIT are distinct provenance fields; the latter is not used to reinterpret the accepted state.",
        ],
        "flags": dict(FLAGS_FALSE),
        "acceptance": {
            "all_artifact_zip_digests_match_github_api": True,
            "all_internal_manifests_match_pinned_payload_digests": True,
            "accepted_state_digests_match_existing_candidate_witnesses": True,
            "operator_residual_reconstruction_checked": True,
            "source_mapping_arithmetic_checked": True,
            "production_source_or_vendor_changed": False,
            **FLAGS_FALSE,
        },
        "limitations": [
            "No new directional tangent was computed; this report re-summarizes the exact existing #420 diagnostic artifacts.",
            "An algebraic decomposition closure is not an independent continuum derivation.",
            "The full matter-coupled identity, source sign/normalization, connection-coupling convention, and total sector conservation are not admitted.",
        ],
    }


def write_failure(exc: BaseException) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema": "full-matter-coupled-propagation-candidate/v1",
        "status": "failed_closed",
        "mode": "reuse_verified_existing_measurements_no_rhs_recomputation",
        "error_type": type(exc).__name__,
        "error": str(exc),
        "traceback": traceback.format_exc(),
        "flags": dict(FLAGS_FALSE),
        "scientific_admission": False,
    }
    path = OUT / "summary.json"
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    (OUT / "SHA256SUMS.txt").write_text(f"{sha256_file(path)}  summary.json\n", encoding="utf-8")


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    try:
        report = build_report()
        report_path = OUT / "summary.json"
        report_path.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
        (OUT / "SHA256SUMS.txt").write_text(f"{sha256_file(report_path)}  summary.json\n", encoding="utf-8")
        print(json.dumps({
            "schema": report["schema"],
            "status": report["status"],
            "mode": report["mode"],
            "absolute_global_max": report["absolute_residual"]["global_max"],
            "observed_order": {k: v for k, v in report["observed_order"].items() if k != "basis"},
            "flags": report["flags"],
            "output": str(report_path),
        }, indent=2))
        return 0
    except Exception as exc:
        write_failure(exc)
        print(json.dumps({"status": "failed_closed", "error_type": type(exc).__name__, "error": str(exc), "flags": FLAGS_FALSE}, indent=2), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
