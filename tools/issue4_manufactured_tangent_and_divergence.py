#!/usr/bin/env python3
"""Diagnostic-only RHS tangent and independent covariant-divergence channels."""
from __future__ import annotations

import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import zipfile
from types import SimpleNamespace

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DIAG = Path(os.environ["ISSUE4_DIAGNOSTIC_ROOT"]).resolve()
BASELINE = Path(os.environ["ISSUE4_HELPER_BASELINE_ROOT"]).resolve()
OUT = Path(os.environ.get("ISSUE4_CHANNEL_OUT", "issue4-channel-results")).resolve()
OUT.mkdir(parents=True, exist_ok=True)
EPSILONS = (2.0e-5, 1.0e-5)
N, R_MAX = 60, 40.0
EXPECTED_DIAG = "4aec88d0e3ff028e0a65a8ac36b30948625df792"
EXPECTED_BASELINE = "7f47669a1c0edfafef9fc5f27b2db80840d3f035"
EXPECTED_VENDOR = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
EXPECTED_TANGENT_BLOB = "3f00a4e0a2063b86d2e007c5a825f9edfff61a4c"
EXPECTED_TOTAL_BLOB = "c10f9e8b9e654d5195fd55c493c3fe06ba5c7082"


def git(root: Path, *args: str) -> str:
    return subprocess.check_output(["git", "-C", str(root), *args], text=True).strip()


def digest_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def load_pinned_module(name: str, path: Path, env_name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load pinned helper: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    old = os.environ.get(env_name)
    os.environ[env_name] = str(BASELINE)
    try:
        spec.loader.exec_module(module)
    finally:
        if old is None:
            os.environ.pop(env_name, None)
        else:
            os.environ[env_name] = old
    return module


def construct_production_state(arrays, r0, vacuum, tangent):
    from engine.production_kernel import ProductionState
    from engine.scalar_system import ScalarFields
    from engine.matter_system import ConservedSpecies
    from engine.v55_matter import V55MatterState

    r = np.asarray(arrays["grid.centers"], dtype=np.float64)
    grid = r0.SphericalCellGrid(r.size, float(r[0] * 2.0 * r.size))
    if not np.array_equal(np.asarray(grid.centers), r):
        raise RuntimeError("pinned grid differs from constructed state centers")
    geometry = vacuum.VacuumState(**{
        key: np.asarray(arrays[f"geometry.{key}"], dtype=np.float64).copy()
        for key in ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B")
    })
    scalars = ScalarFields(*(
        np.asarray(arrays[f"scalars.{key}"], dtype=np.float64).copy()
        for key in ("S", "PS", "D", "PD", "phi", "Pi")
    ))
    species = {}
    for name in ("dark_matter", "baryons", "radiation"):
        species[name] = ConservedSpecies(*(
            np.asarray(arrays[f"matter.{name}.{field}"], dtype=np.float64).copy()
            for field in ("rest", "energy_t", "momentum_r")
        ))
    matter = V55MatterState(species["dark_matter"], species["baryons"], species["radiation"])
    return ProductionState(
        grid=grid, geometry=geometry, scalars=scalars, matter=matter,
        t=0.0, tau=0.0, e_folds=0.0,
    )


def add_projection_arrays(out, prefix, projection):
    for key in ("rho", "pr", "pt", "j"):
        out[f"{prefix}.total.{key}"] = np.asarray(projection[key], dtype=np.float64)
    for sector, values in projection["scalar_parts"].items():
        for key, value in values.items():
            out[f"{prefix}.scalar.{sector}.{key}"] = np.asarray(value, dtype=np.float64)
    for sector, values in projection["fluid_parts"].items():
        for key, value in values.items():
            out[f"{prefix}.fluid.{sector}.{key}"] = np.asarray(value, dtype=np.float64)


def add_nested_arrays(out, prefix, values):
    if isinstance(values, dict):
        for key, value in values.items():
            add_nested_arrays(out, f"{prefix}.{key}", value)
    elif isinstance(values, np.ndarray):
        out[prefix] = np.asarray(values, dtype=np.float64)


def main():
    diag_commit = git(DIAG, "rev-parse", "HEAD")
    baseline_commit = git(BASELINE, "rev-parse", "HEAD")
    vendor_root = DIAG / "vendor" / "bb-palatini-unified-r0"
    vendor_commit = git(vendor_root, "rev-parse", "HEAD")
    tangent_blob = git(DIAG, "rev-parse", "HEAD:tools/audit_current_source_momentum_constraint_tangent.py")
    total_blob = git(DIAG, "rev-parse", "HEAD:tools/audit_total_stress_energy_divergence.py")
    if (diag_commit, baseline_commit, vendor_commit, tangent_blob, total_blob) != (
        EXPECTED_DIAG, EXPECTED_BASELINE, EXPECTED_VENDOR, EXPECTED_TANGENT_BLOB, EXPECTED_TOTAL_BLOB
    ):
        raise RuntimeError("one or more diagnostic/helper/vendor pins differ from the locked identities")

    sys.path.insert(0, str(ROOT))
    sys.path.insert(0, str(vendor_root))
    sys.path.insert(0, str(DIAG))
    from tools.issue4_manufactured_regular_family import construct_regular_family
    import r0_operators
    import vacuum_bssn

    tangent = load_pinned_module(
        "_issue4_manufactured_signed_tangent",
        DIAG / "tools" / "audit_current_source_momentum_constraint_tangent.py",
        "COSMOS_PINNED_SOURCE_ROOT",
    )
    total = load_pinned_module(
        "_issue4_manufactured_total_stress",
        DIAG / "tools" / "audit_total_stress_energy_divergence.py",
        "COSMOS_PINNED_SOURCE_ROOT",
    )
    tangent_provenance = tangent.provenance()
    total_provenance = total.verify_provenance()

    r = (np.arange(N, dtype=np.float64) + 0.5) * (R_MAX / N)
    arrays = construct_regular_family(r, 0.0)
    state = construct_production_state(arrays, r0_operators, vacuum_bssn, tangent)
    initial_digest = tangent.state_digest(state)
    kernel = tangent.V55ProductionKernel()

    # Hard tripwire: no production step is permitted in this observer.
    def forbidden_step(*args, **kwargs):
        raise RuntimeError("production step forbidden in manufactured tangent/divergence observer")
    kernel.step = forbidden_step

    rates = tangent.rhs_direction(kernel, state)
    base_projection = total.independent_total_projections(state)
    if not base_projection["projection_crosscheck_passed"]:
        raise RuntimeError(f"base projection cross-check failed: {base_projection['projection_crosscheck']}")

    retained: dict[str, np.ndarray] = {}
    for name, value in arrays.items():
        retained[f"state.base.{name}"] = np.asarray(value, dtype=np.float64)
    add_projection_arrays(retained, "projection.base", base_projection)
    family_projection_reports = {"base": {
        "state_digest": initial_digest,
        "projection_crosscheck_passed": True,
        "projection_crosscheck": base_projection["projection_crosscheck"],
    }}
    # Preserve all four analytically constructed +/- family members in this
    # same complete fixture archive, not only in the earlier projection run.
    for family_eps in EPSILONS:
        family_label = f"eps{family_eps:.0e}".replace("+", "")
        for sign, signed_eps in (("plus", family_eps), ("minus", -family_eps)):
            member = construct_regular_family(r, signed_eps)
            member_state = construct_production_state(member, r0_operators, vacuum_bssn, tangent)
            member_digest = tangent.state_digest(member_state)
            member_projection = total.independent_total_projections(member_state)
            if not member_projection["projection_crosscheck_passed"]:
                raise RuntimeError(
                    f"family projection cross-check failed for {family_label}_{sign}: "
                    f"{member_projection['projection_crosscheck']}"
                )
            member_key = f"{family_label}_{sign}"
            for key, value in member.items():
                retained[f"state.{member_key}.{key}"] = np.asarray(value, dtype=np.float64)
            add_projection_arrays(retained, f"projection.{member_key}", member_projection)
            if tangent.state_digest(member_state) != member_digest:
                raise RuntimeError(f"family state mutated during projection: {member_key}")
            family_projection_reports[member_key] = {
                "state_digest": member_digest,
                "projection_crosscheck_passed": True,
                "projection_crosscheck": member_projection["projection_crosscheck"],
            }
    tangent_reports = {}
    divergence_reports = {}
    epsilon_sensitivity = {}

    for eps in EPSILONS:
        label = f"eps{eps:.0e}".replace("+", "")
        tangent_item = tangent.tangent_for_epsilon(kernel, state, rates, eps)
        tangent_reports[label] = {
            "epsilon": eps,
            "closure_ok": bool(tangent_item["closure_ok"]),
            "closure_max_abs": float(tangent_item["closure_max_abs"]),
            "closure_max_bound": float(tangent_item["closure_max_bound"]),
            "Mdot_geom_tangent_closure_ok": bool(tangent_item["Mdot_geom_tangent_closure_ok"]),
            "Mdot_geom_tangent_closure_max_abs": float(tangent_item["Mdot_geom_tangent_closure_max_abs"]),
            "Mdot_geom_tangent_closure_max_bound": float(tangent_item["Mdot_geom_tangent_closure_max_bound"]),
            "plus_minus_current_projection_closure_ok": bool(tangent_item["plus_minus_current_projection_closure_ok"]),
        }
        for key in ("Mdot", "Mdot_geom", "Mdot_reconstructed", "closure_abs", "closure_bound",
                    "Mdot_geom_tangent_closure_abs", "Mdot_geom_tangent_closure_bound"):
            retained[f"tangent.{label}.{key}"] = np.asarray(tangent_item[key], dtype=np.float64)
        add_nested_arrays(retained, f"tangent.{label}.sector", tangent_item["Mdot_sector"])
        add_nested_arrays(retained, f"tangent.{label}.geometry_terms", tangent_item["Mdot_geom_terms"])

        plus_state = tangent.shifted_state(kernel, state, rates, eps)
        minus_state = tangent.shifted_state(kernel, state, rates, -eps)
        plus_projection = total.independent_total_projections(plus_state)
        minus_projection = total.independent_total_projections(minus_state)
        if not plus_projection["projection_crosscheck_passed"] or not minus_projection["projection_crosscheck_passed"]:
            raise RuntimeError(f"directional-state projection cross-check failed at {label}")
        add_projection_arrays(retained, f"projection.{label}.plus", plus_projection)
        add_projection_arrays(retained, f"projection.{label}.minus", minus_projection)

        plus_obj = total._covariant_objects_from_components(plus_state, plus_projection)
        minus_obj = total._covariant_objects_from_components(minus_state, minus_projection)
        denom = 2.0 * eps
        dt_density_t = (plus_obj["density_t"] - minus_obj["density_t"]) / denom
        dt_density_r = (plus_obj["density_r"] - minus_obj["density_r"]) / denom
        dt_gtt = (plus_obj["gtt"] - minus_obj["gtt"]) / denom
        plus_grr = plus_state.geometry.a / plus_state.geometry.X**2
        minus_grr = minus_state.geometry.a / minus_state.geometry.X**2
        plus_gth = plus_state.geometry.b * r**2 / plus_state.geometry.X**2
        minus_gth = minus_state.geometry.b * r**2 / minus_state.geometry.X**2
        dt_grr = (plus_grr - minus_grr) / denom
        dt_gth = (plus_gth - minus_gth) / denom
        alpha = np.asarray(state.geometry.alpha, dtype=np.float64)
        gamma_rr = np.asarray(state.geometry.a / state.geometry.X**2, dtype=np.float64)
        gamma_thth = np.asarray(state.geometry.b * r**2 / state.geometry.X**2, dtype=np.float64)
        d1 = lambda values, parity: state.grid.cell_derivative_fourth(values, parity=parity)
        args = (
            r, d1, alpha, gamma_rr, gamma_thth,
            base_projection["rho"], base_projection["j"],
            base_projection["pr"], base_projection["pt"],
            dt_density_t, dt_density_r, dt_gtt, dt_grr, dt_gth,
        )
        point = total.covariant_divergence_terms_from_profiles(*args)
        fv = total.covariant_divergence_fv_terms_from_profiles(*args)
        for channel_name, result in (("pointwise", point), ("finite_volume", fv)):
            retained[f"divergence.{label}.{channel_name}.div_T_t"] = np.asarray(result["div_T_t"], dtype=np.float64)
            retained[f"divergence.{label}.{channel_name}.div_T_r"] = np.asarray(result["div_T_r"], dtype=np.float64)
            add_nested_arrays(retained, f"divergence.{label}.{channel_name}.terms", result["terms"])
        divergence_reports[label] = {
            "plus_projection_crosscheck_passed": bool(plus_projection["projection_crosscheck_passed"]),
            "minus_projection_crosscheck_passed": bool(minus_projection["projection_crosscheck_passed"]),
            "pointwise_div_T_r_max_abs": float(np.max(np.abs(point["div_T_r"]))),
            "finite_volume_div_T_r_max_abs": float(np.max(np.abs(fv["div_T_r"]))),
            "pointwise_component_sum_gap_max_abs": float(np.max(np.abs(point["terms"]["component_sum_gap_div_T_r"]))),
            "finite_volume_component_sum_gap_max_abs": float(np.max(np.abs(fv["terms"]["component_sum_gap_div_T_r"]))),
        }

    labels = list(tangent_reports)
    epsilon_sensitivity["max_abs_Mdot_difference"] = float(np.max(np.abs(
        retained[f"tangent.{labels[0]}.Mdot"] - retained[f"tangent.{labels[1]}.Mdot"]
    )))
    if tangent.state_digest(state) != initial_digest:
        raise RuntimeError("base manufactured state mutated during tangent/divergence observer")

    npz = OUT / "arrays.npz"
    np.savez_compressed(npz, **retained)
    manifest = {
        "schema": "issue4_manufactured_complete_fixture_v1",
        "diagnostic_only": True,
        "state_class": "manufactured_regular_family",
        "resolution": N,
        "r_max": R_MAX,
        "directional_epsilons": list(EPSILONS),
        "complete_constructed_state_classes": [
            "base", "eps2e-05_plus", "eps2e-05_minus",
            "eps1e-05_plus", "eps1e-05_minus",
        ],
        "family_projection_reports": family_projection_reports,
        "production_rhs_evaluated": True,
        "production_steps": 0,
        "scientific_admission": False,
        "acceptance_oracle_satisfied": False,
        "base_state_digest": initial_digest,
        "diagnostic_source_commit": diag_commit,
        "helper_baseline_commit": baseline_commit,
        "vendor_commit": vendor_commit,
        "tangent_helper_blob": tangent_blob,
        "total_stress_helper_blob": total_blob,
        "tangent_provenance": tangent_provenance,
        "total_stress_provenance": total_provenance,
        "tangent_reports": tangent_reports,
        "divergence_reports": divergence_reports,
        "epsilon_sensitivity": epsilon_sensitivity,
        "array_count": len(retained),
        "array_names": sorted(retained),
        "arrays_member": npz.name,
        "arrays_sha256": digest_file(npz),
        "archive_members": ["SHA256SUMS.txt", "report.json", "arrays.npz"],
        "interpretation_limits": [
            "The tangent is the pinned semi-discrete production RHS directional derivative from the manufactured base state; no production step or trajectory is taken.",
            "Both epsilon channels use the same frozen base state and the pinned two-sided RHS direction with the helper's CMC lapse re-solve.",
            "Pointwise and finite-volume covariant divergences are independent spatial-flux treatments using the same signed time-density and metric derivatives from the two-sided RHS direction.",
            "This arithmetic diagnostic does not establish physical conservation closure or scientific admission.",
        ],
    }
    mp = OUT / "report.json"
    mp.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    sums = OUT / "SHA256SUMS.txt"
    sums.write_text(
        f"{digest_file(mp)}  {mp.name}\n{digest_file(npz)}  {npz.name}\n"
    )
    fixture_zip = OUT / "issue4_manufactured_regular_family_N60.zip"
    with zipfile.ZipFile(fixture_zip, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        archive.write(sums, arcname="SHA256SUMS.txt")
        archive.write(mp, arcname="report.json")
        archive.write(npz, arcname="arrays.npz")
    print(json.dumps({
        "status": "PASS",
        "base_state_digest": initial_digest,
        "family_projection_reports": family_projection_reports,
        "tangent": tangent_reports,
        "divergence": divergence_reports,
        "epsilon_sensitivity": epsilon_sensitivity,
        "array_count": len(retained),
        "arrays_sha256": digest_file(npz),
        "fixture_zip": fixture_zip.name,
        "fixture_zip_sha256": digest_file(fixture_zip),
        "production_steps": 0,
    }, indent=2))


if __name__ == "__main__":
    main()
