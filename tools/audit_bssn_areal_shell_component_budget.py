#!/usr/bin/env python3
"""Research-only shell mass-rate/current reconstruction component budget.

Reuses the existing pinned fixed-R observer on the exact accepted t=0.60
full-profile artifacts. For each shell, compare the cell-centred baseline with
parity-aware face and hybrid reconstruction variants. Decompose

    residual_variant - residual_cell
      = (mass_rate_variant - mass_rate_cell)
        - (net_current_variant - net_current_cell)

so the face-minus-cell difference is attributed to the two numerical sides of
the balance. This is operator accounting, not a choice of physical authority.
No production source is changed.
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
N = int(os.environ.get("AREAL_SHELL_BUDGET_RESOLUTION", "120"))
FULL_DIR = Path(os.environ.get("FULL_PROFILE_ARTIFACT_DIR", f"reference/full/N{N}"))
OUT = Path(os.environ.get("AREAL_SHELL_BUDGET_OUT", f"results/bssn-areal-shell-component-budget-N{N}-2026-10-25"))
TARGET_TIME = 0.60
CFL = 0.03

VARIANTS = {
    "cell_centered_linear": ("shell_balances", "net_current_midpoint"),
    "parity_face_primitive_current": ("face_reconstructed_shell_balances", "net_face_current_midpoint"),
    "cell_derivative_face_geometry_current": ("cell_derivative_face_shell_balances", "net_face_current_midpoint"),
    "metric_face_chi_combined_cell_current": ("metric_chi_cell_flux_shell_balances", "net_current_midpoint"),
    "cell_chi_face_primitive_current": ("cell_chi_primitive_flux_shell_balances", "net_current_midpoint"),
    "cell_chi_combined_cell_current": ("cell_chi_cell_flux_shell_balances", "net_current_midpoint"),
    "cell_centered_cubic_interpolation": ("cell_shell_balances_cubic_interpolation", "net_current_midpoint"),
    "parity_face_cubic_interpolation": ("face_shell_balances_cubic_interpolation", "net_current_midpoint"),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def finite_float(value: Any, label: str) -> float:
    result = float(value)
    if not np.isfinite(result):
        raise RuntimeError(f"{label} is non-finite")
    return result


def shell_rows(probe: dict[str, Any], key: str, current_key: str) -> dict[str, dict[str, float]]:
    raw = probe.get(key)
    if not isinstance(raw, list):
        raise RuntimeError(f"probe does not contain shell variant {key}")
    rows: dict[str, dict[str, float]] = {}
    for item in raw:
        rin = finite_float(item["R_in"], f"{key}.R_in")
        rout = finite_float(item["R_out"], f"{key}.R_out")
        shell = f"{rin:g}-{rout:g}"
        mass_rate = finite_float(item["shell_mass_rate"], f"{key}.{shell}.shell_mass_rate")
        current_value = item.get(current_key)
        current_field = current_key
        if current_value is None:
            for fallback_key in ("net_face_current_midpoint", "net_current_midpoint", "net_current"):
                if fallback_key in item:
                    current_value = item[fallback_key]
                    current_field = fallback_key
                    break
        if current_value is None:
            raise RuntimeError(f"{key}.{shell}: no recognized net-current field; found={sorted(item.keys())}")
        net_current = finite_float(current_value, f"{key}.{shell}.{current_field}")
        residual = finite_float(item["residual"], f"{key}.{shell}.residual")
        closure = residual - (mass_rate - net_current)
        if abs(closure) > 5e-12 * max(1.0, abs(mass_rate), abs(net_current), abs(residual)):
            raise RuntimeError(f"{key}.{shell}: residual != mass-rate minus net-current")
        rows[shell] = {
            "R_in": rin,
            "R_out": rout,
            "shell_mass_rate": mass_rate,
            "net_current": net_current,
            "residual": residual,
            "residual_decomposition_closure": closure,
        }
    return rows


def decompose_probe(probe: dict[str, Any]) -> dict[str, Any]:
    variants = {
        label: shell_rows(probe, fields[0], fields[1])
        for label, fields in VARIANTS.items()
    }
    baseline = variants["cell_centered_linear"]
    output: dict[str, Any] = {"shells": {}, "variant_maxima": {}}

    for label, rows in variants.items():
        output["variant_maxima"][label] = {
            "max_abs_shell_residual": max((abs(row["residual"]) for row in rows.values()), default=0.0),
            "max_abs_decomposition_closure": max((abs(row["residual_decomposition_closure"]) for row in rows.values()), default=0.0),
        }

    for shell, base in baseline.items():
        row: dict[str, Any] = {"cell_centered_baseline": base, "variants": {}}
        for label, rows in variants.items():
            if shell not in rows:
                continue
            variant = rows[shell]
            delta_mass = variant["shell_mass_rate"] - base["shell_mass_rate"]
            delta_current = variant["net_current"] - base["net_current"]
            delta_residual = variant["residual"] - base["residual"]
            closure = delta_residual - (delta_mass - delta_current)
            row["variants"][label] = {
                "shell_mass_rate": variant["shell_mass_rate"],
                "net_current": variant["net_current"],
                "residual": variant["residual"],
                "delta_mass_rate_vs_cell": delta_mass,
                "delta_net_current_vs_cell": delta_current,
                "delta_residual_vs_cell": delta_residual,
                "delta_budget_closure": closure,
            }
            if abs(closure) > 5e-12 * max(1.0, abs(delta_residual), abs(delta_mass), abs(delta_current)):
                raise RuntimeError(f"{label}.{shell}: crossed-component difference failed closure")
        output["shells"][shell] = row

    baseline_max = output["variant_maxima"]["cell_centered_linear"]["max_abs_shell_residual"]
    face_max = output["variant_maxima"]["parity_face_primitive_current"]["max_abs_shell_residual"]
    output["headline"] = {
        "cell_centered_max_abs_shell_residual": baseline_max,
        "parity_face_max_abs_shell_residual": face_max,
        "parity_face_minus_cell_max_abs_delta_residual": max(
            (abs(row["variants"]["parity_face_primitive_current"]["delta_residual_vs_cell"])
             for row in output["shells"].values()
             if "parity_face_primitive_current" in row["variants"]),
            default=0.0,
        ),
        "operator_choice_closes_physical_identity": False,
    }
    return output


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
        raise RuntimeError("accepted state is not at t=0.60")
    input_digest = saved.tangent.state_digest(accepted)

    dt = CFL * float(accepted.grid.dr)
    reports = []
    for label, h in (("dt", dt), ("dt_half", dt / 2.0)):
        before = saved.tangent.state_digest(accepted)
        probe = areal.compare_step(kernel, accepted, h)
        after = saved.tangent.state_digest(accepted)
        if before != after:
            raise RuntimeError(f"accepted input mutated in {label} arm")
        if probe.get("input_state_immutable") is not True:
            raise RuntimeError(f"underlying observer did not affirm immutable input in {label} arm")
        if probe.get("max_abs_beta_start_end", 0.0) > 1e-14:
            raise RuntimeError(f"pinned zero-shift precondition failed in {label} arm")
        reports.append({
            "arm": label,
            "dt": float(h),
            "steps_advanced_from_accepted_input": 1,
            "accepted_input_immutable": True,
            "probe_shell_component_budget": decompose_probe(probe),
            "probe_accepted_state_bookkeeping": {
                "input_state_immutable": probe["input_state_immutable"],
                "max_abs_shell_telescoping_closure": probe["max_abs_shell_telescoping_closure"],
                "face_reconstruction_max_abs_difference_from_cell_shell_residual":
                    probe["face_reconstruction_max_abs_difference_from_cell_shell_residual"],
                "max_abs_Rr_face_metric_vs_cell_gap_r_le10":
                    probe["max_abs_Rr_face_metric_vs_cell_gap_r_le10"],
                "max_abs_Rt_face_metric_vs_cell_gap_r_le10":
                    probe["max_abs_Rt_face_metric_vs_cell_gap_r_le10"],
            },
        })
    if saved.tangent.state_digest(accepted) != input_digest:
        raise RuntimeError("accepted state digest changed across component-budget arms")

    report = {
        "schema": "bssn-areal-shell-component-budget/v1",
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
            "areal_radius_map": "R=r*sqrt(b)/X",
            "fixed_R_mass": "m=R*(1-chi)/2; chi interpolated to common R before differencing",
            "baseline": "cell-centred shell mass-rate and endpoint-time-averaged net current",
            "difference_identity": "Delta(residual) = Delta(shell mass-rate) - Delta(net current)",
        },
        "probes": reports,
        "interpretation_guard": [
            "Each variant is a separately named discrete reconstruction; no variant is promoted as physical authority.",
            "Small mass/current component difference does not by itself establish the continuum law; component cancellation and refinement must both be inspected.",
            "The exact shell telescoping and component-difference closures are bookkeeping checks, not total stress-energy conservation.",
            "No source normalization, stress formula, production equation, initializer, vendor, or dependency was changed.",
        ],
        "acceptance": {
            "source_pins_match": True,
            "accepted_artifact_sha_verified": True,
            "accepted_input_immutable": True,
            "zero_shift_precondition": True,
            "shell_component_decomposition_closes": True,
            "face_vs_cell_residual_component_budget_recorded": True,
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
    checksum_path = OUT / "SHA256SUMS.txt"
    checksum_path.write_text(f"{sha256(report_path)}  report.json\n", encoding="utf-8")

    print(json.dumps({
        "status": report["status"],
        "resolution": N,
        "source_vendor_pins_match": True,
        "accepted_input_immutable": True,
        "scientific_admission": False,
        "arms": [{
            "arm": item["arm"],
            **item["probe_shell_component_budget"]["headline"],
            "shells": item["probe_shell_component_budget"]["shells"],
        } for item in reports],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
