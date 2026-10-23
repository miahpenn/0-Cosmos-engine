#!/usr/bin/env python3
"""One-step, current-source fixed-areal-radius ledger witness; research only.

This uses the existing fixed-R algebra and passive production-stage observer,
but has its own current-source pins and a six-component matter-source observer.
It never calls the historical audit's old SOURCE_PINS verifier and never treats
a small physical balance remainder as a pass condition.
"""
from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import tempfile

import numpy as np

# Running this file directly puts tools/ first on sys.path, not the repository root.
_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from engine import matter_rhs as matter_rhs_module
from engine import scalar_system
from engine.campaign import CampaignConfig
from engine.run_production import run_campaign
from engine import production_kernel as production_kernel_module
from engine.v55_matter import dm_density, metric_slice_from_q
from tools import audit_em_total_energy_budget as fixed_r

# Blob pins for the current production equation/geometry path at base
# 4aec88d0e3ff028e0a65a8ac36b30948625df792. These are intentionally separate
# from the historical audit's SOURCE_PINS, which must remain frozen.
SOURCE_BASE = "4aec88d0e3ff028e0a65a8ac36b30948625df792"
CURRENT_SOURCE_PINS = {
    "engine/campaign.py": "119942d9debf63b99997cda96c541754e7e8aded",
    "engine/cmc_gauge.py": "7b49b5751abebd4faf0130208032cd2080b71f64",
    "engine/cosmology_observables.py": "4d6f39308685bae5e5bad2e17baa9c4043b67b4e",
    "engine/handoff.py": "1e0b8a8a126c67dafaa5d370babdd519cbedf974",
    "engine/invariant_diagnostics.py": "4bbcc6174c978b07c4f2175d1dcfeccb2d0c12b1",
    "engine/matter_rhs.py": "0153c232ee15b5a9d3dbac37ea858e4d2d730011",
    "engine/matter_system.py": "0969f4d996e4f8e0fc5b70a7dff147f024907302",
    "engine/native_geometry/__init__.py": "796d338d7cf2ed87bcb6ae63808dcf08110aeff8",
    "engine/native_geometry/constraints.py": "124d1f237e54d082b4d76213f9ebe144e02b02ba",
    "engine/native_geometry/evolution.py": "7b47d700680b12ab3c82cb6a556d40580b103ce9",
    "engine/native_geometry/geometry.py": "1bf12c1139572ff4597eaec0570aa43f2563fd91",
    "engine/native_geometry/operators.py": "2748be304b58328e2620f0a7fde49b01c0794ae2",
    "engine/production_contract.py": "0a3b6f0a4223f93a5e33ab30706dc803f69643ca",
    "engine/production_kernel.py": "1720615628a98a8d7cfc9fcadcc984fe28c88be1",
    "engine/radial_maxwell.py": "81269f1d1334c46c4c2f44ba816d792722c0a38e",
    "engine/run_production.py": "726f3fd09d46142862574d348aa0755f83bc0612",
    "engine/scalar_system.py": "e2c2dae9a016b4d6454cf74a9eae8b0167d5733d",
    "engine/stress_energy.py": "ac3e7518c62bb72e2218a61489fa486bb787809f",
    "engine/v55_initial.py": "7db2e1a10cf16fb40e4221b0df504fdf39aae611",
    "engine/v55_matter.py": "0ee3c51fac7dc46c63feba2747a750542bcee3d0",
    "engine/valencia.py": "03ed86398f424c03d9d6a3e646da9222a34fac5c",
    "engine/v55_pirk_adapter.py": "be28ceea63fa1b59efc5cc5c460b9df19ae1cbc3",
    "engine/worldtube.py": "aea19090cb0fd2ebe5481a8b7af0d26f119a9ad6",
    "tools/audit_em_total_energy_budget.py": "b8feca93ca3231bcdc49cbe565a206013daeb9e8",
    "tools/fixed_r_discriminators.py": "22f9b0333487487ddcf3443939266418f311e1fc",
}
VENDOR_PATH = "vendor/bb-palatini-unified-r0"
VENDOR_PIN = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
RESOLUTION = 24
R_MAX = 24.0
CFL = 0.01
FINAL_TIME = 0.01
CHECKPOINT_INTERVAL = 0.01
ROUND_OFF_MULTIPLIER = 128.0


def _git(root: Path, *args: str) -> str:
    result = subprocess.run(
        ["git", *args], cwd=root, check=True, capture_output=True, text=True
    )
    return result.stdout.strip()


def verify_current_source_pins(root: Path) -> dict:
    """Fail closed on changed current-source blobs and the vendor gitlink."""
    blobs = {}
    for path, expected in CURRENT_SOURCE_PINS.items():
        actual = _git(root, "rev-parse", f"HEAD:{path}")
        blobs[path] = {"expected": expected, "observed": actual, "match": actual == expected}
        if actual != expected:
            raise RuntimeError(
                f"current-source pin mismatch for {path}: expected {expected}, got {actual}"
            )
    vendor = _git(root, "rev-parse", f"HEAD:{VENDOR_PATH}")
    if vendor != VENDOR_PIN:
        raise RuntimeError(
            f"vendor gitlink mismatch: expected {VENDOR_PIN}, got {vendor}"
        )
    checkout = _git(root, "rev-parse", "HEAD")
    github_sha = os.environ.get("GITHUB_SHA", "").strip()
    if github_sha and github_sha != checkout:
        raise RuntimeError(
            f"GITHUB_SHA={github_sha} does not match checkout HEAD={checkout}"
        )
    return {
        "source_base_anchor": SOURCE_BASE,
        "checkout_sha": checkout,
        "github_sha": github_sha or None,
        "pr_source_head_sha": os.environ.get("PR_SOURCE_HEAD_SHA") or None,
        "production_blob_pins": blobs,
        "vendor_gitlink": {"expected": VENDOR_PIN, "observed": vendor, "match": vendor == VENDOR_PIN},
    }


def _arr(value) -> np.ndarray:
    result = np.asarray(value, dtype=np.float64)
    if not np.all(np.isfinite(result)):
        raise RuntimeError("captured operand contains a non-finite value")
    return result.copy()


def _profile_check(actual, terms: list[np.ndarray], *, label: str) -> dict:
    actual_arr = _arr(actual)
    term_arrays = [_arr(term) for term in terms]
    if any(term.shape != actual_arr.shape for term in term_arrays):
        raise RuntimeError(f"{label}: term shapes do not match the observed RHS")
    expected = term_arrays[0].copy()
    for term in term_arrays[1:]:
        expected = expected + term
    residual = actual_arr - expected
    scale = np.abs(actual_arr)
    for term in term_arrays:
        scale = scale + np.abs(term)
    bound = ROUND_OFF_MULTIPLIER * np.finfo(np.float64).eps * np.maximum(
        scale, np.finfo(np.float64).tiny
    )
    ratio = np.abs(residual) / bound
    return {
        "label": label,
        "max_abs_residual": float(np.max(np.abs(residual))),
        "max_roundoff_bound": float(np.max(bound)),
        "max_bound_ratio": float(np.max(ratio)),
        "within_roundoff_bound": bool(np.all(np.abs(residual) <= bound)),
        "residual_profile": residual,
        "roundoff_bound_profile": bound,
    }


class CurrentSourceLedgerKernel(fixed_r.AuditedKernel):
    """Reuse stage observer, additionally retain actual scalar RHS profiles."""

    def _rhs(self, state):
        scalar_rhs, matter_rhs, metric_derivatives = super()._rhs(state)
        if self._active_step is not None and self._active_step.get("rhs_stages"):
            stage = self._active_step["rhs_stages"][-1]
            stage["scalar_rhs_profiles"] = {
                name: _arr(getattr(scalar_rhs, name))
                for name in ("S", "PS", "D", "PD", "phi", "Pi")
            }
        return scalar_rhs, matter_rhs, metric_derivatives


from contextlib import contextmanager


@contextmanager
def instrument_current_species_rhs(kernel: CurrentSourceLedgerKernel):
    """Capture the current six returned species components without modifying physics."""
    original = matter_rhs_module._species_rhs_components

    def wrapper(*args, **kwargs):
        result = original(*args, **kwargs)
        active = kernel._active_rhs
        if active is None:
            return result
        if len(result) != 6:
            raise RuntimeError(
                "current matter RHS component contract changed: expected "
                "(face_flux, source_energy, source_momentum, q_energy, q_momentum, q_rest)"
            )
        stage_state = active["stage_state"]
        species_state = args[2]
        name = fixed_r._species_name(stage_state, species_state)
        face_flux, source_energy, source_momentum, q_energy, q_momentum, q_rest = result
        active["components"][name] = {
            "face_flux": _arr(face_flux),
            "source_energy": _arr(source_energy),
            "source_momentum": _arr(source_momentum),
            "dm_exchange_energy": _arr(q_energy),
            "dm_exchange_momentum": _arr(q_momentum),
            "rest_current_source": _arr(q_rest),
        }
        return result

    matter_rhs_module._species_rhs_components = wrapper
    try:
        yield
    finally:
        matter_rhs_module._species_rhs_components = original


def _matter_component_ledger(stage: dict, dr: float) -> dict:
    """Reconstruct all three current conservative RHS channels from named operands."""
    rows = {}
    for species, components in stage["components"].items():
        rhs = stage["matter_rhs"][species]
        flux = _arr(components["face_flux"])
        if flux.ndim != 2 or flux.shape[1] != 3:
            raise RuntimeError(f"{species}: face flux must have shape (N+1, 3)")
        inv_dr = 1.0 / float(dr)
        rest_flux_div = -inv_dr * (flux[1:, 0] - flux[:-1, 0])
        energy_flux_div = -inv_dr * (flux[1:, 1] - flux[:-1, 1])
        momentum_flux_div = -inv_dr * (flux[1:, 2] - flux[:-1, 2])
        checks = {
            "rest": _profile_check(
                rhs.rest,
                [rest_flux_div, _arr(components["rest_current_source"])],
                label=f"{species}.rest = -D(face J) + q_D",
            ),
            "energy_t": _profile_check(
                rhs.energy_t,
                [
                    energy_flux_div,
                    _arr(components["source_energy"]),
                    _arr(components["dm_exchange_energy"]),
                ],
                label=f"{species}.energy_t = -D(face T) + geometric source + q_E",
            ),
            "momentum_r": _profile_check(
                rhs.momentum_r,
                [
                    momentum_flux_div,
                    _arr(components["source_momentum"]),
                    _arr(components["dm_exchange_momentum"]),
                ],
                label=f"{species}.momentum_r = -D(face T) + geometric source + q_Sr",
            ),
        }
        profiles = {
            "face_flux": flux,
            "face_flux_divergences": {
                "rest": rest_flux_div,
                "energy_t": energy_flux_div,
                "momentum_r": momentum_flux_div,
            },
            "geometric_sources": {
                "energy": _arr(components["source_energy"]),
                "momentum": _arr(components["source_momentum"]),
            },
            "conditional_exchange_sources": {
                "energy": _arr(components["dm_exchange_energy"]),
                "momentum": _arr(components["dm_exchange_momentum"]),
                "rest_current_q_D": _arr(components["rest_current_source"]),
            },
            "observed_rhs": {
                "rest": _arr(rhs.rest),
                "energy_t": _arr(rhs.energy_t),
                "momentum_r": _arr(rhs.momentum_r),
            },
        }
        rows[species] = {
            "operator_reconstruction_checks": {
                name: {
                    key: value
                    for key, value in check.items()
                    if key not in ("residual_profile", "roundoff_bound_profile")
                }
                for name, check in checks.items()
            },
            "profiles": profiles,
        }
    return rows


def _scalar_exchange_contract(stage: dict, state) -> dict:
    """Compare live Pi RHS against same-state no-exchange evaluation."""
    metric = metric_slice_from_q(state.grid, state.geometry)
    rho_dm = _arr(dm_density(metric, state.matter))
    beta_dm = float(production_kernel_module.BETA_DM)
    baseline = scalar_system.scalar_rhs_arrays(
        state.grid, state.geometry, state.scalars,
        beta_dm=beta_dm, rho_dm=None,
    )
    actual_pi = _arr(stage["scalar_rhs_profiles"]["Pi"])
    baseline_pi = _arr(baseline.Pi)
    expected = _arr(state.geometry.alpha) * beta_dm * rho_dm
    observed = actual_pi - baseline_pi
    residual = observed - expected
    # Subtracting two full Pi RHS profiles can lose precision when the
    # exchange is much smaller than the uncoupled scalar evolution. Bound that
    # reconstruction using the actual operands, not only their small difference.
    scale = np.abs(actual_pi) + np.abs(baseline_pi) + np.abs(observed) + np.abs(expected)
    bound = ROUND_OFF_MULTIPLIER * np.finfo(np.float64).eps * np.maximum(
        scale, np.finfo(np.float64).tiny
    )
    ratio = np.abs(residual) / bound
    check = {
        "label": "live scalar Pi RHS minus same-state no-exchange RHS = alpha beta_DM rho_DM",
        "max_abs_residual": float(np.max(np.abs(residual))),
        "max_roundoff_bound": float(np.max(bound)),
        "max_bound_ratio": float(np.max(ratio)),
        "within_roundoff_bound": bool(np.all(np.abs(residual) <= bound)),
    }
    return {
        "status": "conditional source algebra only; action premise and sector normalization remain unresolved",
        "beta_dm": beta_dm,
        "rho_dm_profile": rho_dm,
        "actual_live_pi_rhs": actual_pi,
        "same_state_pi_rhs_without_exchange": baseline_pi,
        "observed_exchange_difference": observed,
        "expected_lapse_weighted_exchange_source": expected,
        "reconstruction_check": check,
        "difference_from_expected_profile": residual,
    }


def _assert_stage_geometry_matches(stage: dict, geometry_stage: dict) -> dict:
    geometry = stage["stage_state"].geometry
    captured = geometry_stage["geometry"]
    comparisons = {}
    for name in ("a", "b", "X", "alpha", "beta", "Aa", "K"):
        live = _arr(getattr(geometry, name))
        observed = _arr(captured[name])
        exact = bool(np.array_equal(live, observed))
        comparisons[name] = {
            "exact_array_match": exact,
            "max_abs_difference": float(np.max(np.abs(live - observed))),
        }
        if not exact:
            raise RuntimeError(
                f"stage-state/captured geometry mismatch for {name}; refusing mixed-stage ledger"
            )
    return comparisons


def run_current_source_witness(root: Path, output_path: Path | None = None) -> dict:
    """Run exactly one accepted N=24 production step and serialize the live start-stage ledger."""
    root = root.resolve()
    provenance = verify_current_source_pins(root)
    with tempfile.TemporaryDirectory(prefix="current-source-fixed-r-ledger-") as directory:
        config = CampaignConfig(
            resolutions=(RESOLUTION,),
            r_max=R_MAX,
            cfl=CFL,
            final_time=FINAL_TIME,
            checkpoint_interval=CHECKPOINT_INTERVAL,
            output_dir=directory,
        )
        kernel = CurrentSourceLedgerKernel()
        with instrument_current_species_rhs(kernel):
            results = run_campaign(kernel=kernel, config=config)

        if len(results) != 1 or results[0].status != "completed":
            raise RuntimeError(
                "bounded current-source witness did not complete: "
                + repr([(item.status, item.steps, item.t_final) for item in results])
            )
        result = results[0]
        if result.steps != 1 or len(kernel.audit_runs) != 1:
            raise RuntimeError(
                f"expected exactly one accepted production step; got steps={result.steps}, "
                f"audit_runs={len(kernel.audit_runs)}"
            )
        run = kernel.audit_runs[0]
        if len(run["steps"]) != 1:
            raise RuntimeError(f"expected one captured step, got {len(run['steps'])}")
        step = run["steps"][0]
        if len(step["rhs_stages"]) < 1 or len(step["geometry_stages"]) != 3:
            raise RuntimeError("observer did not retain the start RHS and all PIRK geometry stages")
        stage = step["rhs_stages"][0]
        geometry_stage = step["geometry_stages"][0]
        stage_state = stage["stage_state"]
        geometry_match = _assert_stage_geometry_matches(stage, geometry_stage)

        if not math.isclose(float(stage_state.t), float(step["t_before"]), rel_tol=0.0, abs_tol=0.0):
            raise RuntimeError("captured start RHS time differs from the accepted-step start time")

        component_ledger = _matter_component_ledger(stage, float(stage_state.grid.dr))
        all_component_checks = [
            item
            for species_row in component_ledger.values()
            for item in species_row["operator_reconstruction_checks"].values()
        ]
        scalar_exchange = _scalar_exchange_contract(stage, stage_state)
        if not all(item["within_roundoff_bound"] for item in all_component_checks):
            raise RuntimeError("one or more live matter RHS component reconstructions failed")
        if not scalar_exchange["reconstruction_check"]["within_roundoff_bound"]:
            raise RuntimeError("live scalar exchange reconstruction failed")

        sample = fixed_r._worldtube_sample(
            kernel, stage_state, geometry_stage_capture=geometry_stage
        )
        fixed_R = fixed_r._fixed_R_rate_snapshot(sample)
        raw_balance = sample.get("tube_fixed_R_balance_raw_witness")
        if not isinstance(raw_balance, dict):
            raise RuntimeError("fixed-R observer did not produce its signed raw balance witness")

        retries_path = Path(directory) / f"N{RESOLUTION}" / "radiation_step_retries.json"
        retries = json.loads(retries_path.read_text(encoding="utf-8")) if retries_path.exists() else []
        if retries:
            raise RuntimeError(f"bounded one-step witness needed radiation retries: {retries!r}")

        report = {
            "classification": (
                "one-step current-source fixed-areal-radius ledger witness; research only; "
                "not total stress-energy conservation proof"
            ),
            "scientific_admission": False,
            "acceptance_oracle_satisfied": False,
            "physical_fixed_R_residual_threshold_applied": False,
            "run_settings": {
                "resolution": RESOLUTION,
                "r_max": R_MAX,
                "cfl": CFL,
                "final_time": FINAL_TIME,
                "checkpoint_interval": CHECKPOINT_INTERVAL,
                "accepted_steps": int(result.steps),
                "t_final": float(result.t_final),
                "stage_used_for_joint_ledger": "start RHS at accepted-step t_before",
                "stage_time": float(stage_state.t),
                "actual_shift_absmax": float(np.max(np.abs(stage_state.geometry.beta))),
                "fixed_areal_radius": float(fixed_r.FIXED_AREAL_RADIUS),
                "radial_maxwell_enabled": False,
                "radiation_retries": retries,
            },
            "provenance": {
                **provenance,
                "environment": {
                    "python": sys.version,
                    "platform": platform.platform(),
                    "numpy": np.__version__,
                    "OPENBLAS_CORETYPE": os.environ.get("OPENBLAS_CORETYPE"),
                    "NPY_DISABLE_CPU_FEATURES": os.environ.get("NPY_DISABLE_CPU_FEATURES"),
                },
            },
            "same_stage_geometry_match": geometry_match,
            "live_matter_rhs_component_ledger": component_ledger,
            "live_scalar_exchange_contract": scalar_exchange,
            "geometry_rhs_stage_operands": {
                name: _arr(values)
                for name, values in geometry_stage["explicit"].items()
            },
            "geometry_primary_rhs_stage_operands": {
                label: {
                    name: _arr(values)
                    for name, values in geometry_stage[label].items()
                }
                for label in ("primary_l2", "primary_l3", "vacuum_primary_l3", "einstein_source_l3")
            },
            "fixed_R_signed_rate_operators": fixed_R,
            "fixed_R_raw_mass_gradient_and_stress_witness": raw_balance,
            "total_stress_projection_profiles_by_sector": raw_balance.get("stress_sectors", {}),
            "interpretation_limits": [
                "Matter RHS checks reconstruct the live semi-discrete component assembly; they are not an independent integrated conservation proof.",
                "The scalar exchange check validates the conditional source implementation on the same state, not the action premise or physical coupling value.",
                "Both factorized and composite-interpolated fixed-R flux/work values are retained; neither is chosen by residual size.",
                "The fixed-R physical balance residuals are reported signed and unthresholded.",
                "Issue #210 normalization provenance, full evolving total-T closure, and Issue #4 remain open.",
            ],
        }

    safe = fixed_r._json_safe(report)
    if output_path is not None:
        output_path = output_path.resolve()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(json.dumps(safe, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return safe


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = run_current_source_witness(args.repo_root, args.output)
    fixed = report["fixed_R_signed_rate_operators"]
    scalar = report["live_scalar_exchange_contract"]["reconstruction_check"]
    print(json.dumps({
        "classification": report["classification"],
        "checkout_sha": report["provenance"]["checkout_sha"],
        "resolution": report["run_settings"]["resolution"],
        "accepted_steps": report["run_settings"]["accepted_steps"],
        "stage_time": report["run_settings"]["stage_time"],
        "matter_component_checks": sum(
            len(row["operator_reconstruction_checks"])
            for row in report["live_matter_rhs_component_ledger"].values()
        ),
        "all_matter_component_checks_passed": all(
            check["within_roundoff_bound"]
            for row in report["live_matter_rhs_component_ledger"].values()
            for check in row["operator_reconstruction_checks"].values()
        ),
        "scalar_exchange_check_passed": scalar["within_roundoff_bound"],
        "factorized_fixed_R_residual": fixed["tube_chain_rule_fixed_R_residual"],
        "composite_fixed_R_residual": fixed["tube_chain_rule_flux_work_operator_comparison"][
            "fixed_R_residual_local_composite_interpolated"
        ],
        "physical_residual_threshold_applied": False,
        "output": str(args.output.resolve()),
    }, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
