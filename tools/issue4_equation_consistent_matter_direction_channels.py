#!/usr/bin/env python3
"""Build the equation-consistent matter-direction Issue #4 fixture.

No production initializer, production RHS, kernel step, or vendor evolution
step is called. The tangent is a central difference of the declared analytic
family U(s), and both divergence channels use the same frozen triplets.
"""
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
CODE_PATHS = (
    "tools/issue4_full_bianchi_evaluator.py",
    "tools/issue4_zero_step_full_bianchi_observer.py",
    "tools/issue4_accepted_state_fixture.py",
    "tools/issue4_manufactured_regular_family.py",
    "tools/issue4_equation_consistent_matter_direction.py",
    "tools/issue4_equation_consistent_matter_direction_fixture.py",
    "tools/issue4_equation_consistent_matter_direction_channels.py",
    "tools/issue4_equation_consistent_matter_direction_evaluator.py",
    "tests/test_issue4_equation_consistent_matter_direction.py",
    "tests/test_issue4_equation_consistent_matter_direction_fixture.py",
    "tests/test_issue4_accepted_state_fixture.py",
    "tests/test_issue4_bianchi_reference_regression.py",
    "tests/test_issue4_full_bianchi_evaluator_contract.py",
)


def git(root: Path, *args: str) -> str:
    return subprocess.check_output(["git", "-C", str(root), *args], text=True).strip()


def digest_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def digest_file(path: Path) -> str:
    return digest_bytes(path.read_bytes())


def load_pinned_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load pinned helper: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    old = os.environ.get("COSMOS_PINNED_SOURCE_ROOT")
    os.environ["COSMOS_PINNED_SOURCE_ROOT"] = str(BASELINE)
    try:
        spec.loader.exec_module(module)
    finally:
        if old is None:
            os.environ.pop("COSMOS_PINNED_SOURCE_ROOT", None)
        else:
            os.environ["COSMOS_PINNED_SOURCE_ROOT"] = old
    return module


def construct_production_state(arrays, r0, vacuum):
    from engine.production_kernel import ProductionState
    from engine.scalar_system import ScalarFields
    from engine.matter_system import ConservedSpecies
    from engine.v55_matter import V55MatterState

    r = np.asarray(arrays["grid.centers"], dtype=np.float64)
    grid = r0.SphericalCellGrid(r.size, float(r[0] * 2.0 * r.size))
    if not np.array_equal(np.asarray(grid.centers), r):
        raise RuntimeError("pinned vendor grid differs from constructed state centers")
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


def add_list_arrays(out, prefix, values):
    for key, value in values.items():
        out[f"{prefix}.{key}"] = np.asarray(value, dtype=np.float64)


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
        raise RuntimeError("one or more diagnostic/helper/vendor pins differ from locked identities")

    sys.path.insert(0, str(ROOT))
    sys.path.insert(0, str(vendor_root))
    sys.path.insert(0, str(DIAG))
    from tools.issue4_equation_consistent_matter_direction import (
        PROFILE, STATE_ARRAY_NAMES, analytic_matter_rates, construct_equation_consistent_family,
    )
    import r0_operators
    import vacuum_bssn

    tangent = load_pinned_module(
        "_issue4_manufactured_constraint_observer",
        DIAG / "tools" / "audit_current_source_momentum_constraint_tangent.py",
    )
    total = load_pinned_module(
        "_issue4_manufactured_total_stress",
        DIAG / "tools" / "audit_total_stress_energy_divergence.py",
    )
    tangent_provenance = tangent.provenance()
    total_provenance = total.verify_provenance()

    pinned_grid = r0_operators.SphericalCellGrid(N, R_MAX)
    r = np.asarray(pinned_grid.centers, dtype=np.float64)
    rates = analytic_matter_rates(
        r, d1=pinned_grid.cell_derivative_fourth,
        d2=pinned_grid.cell_second_derivative_fourth,
    )
    family_arrays = {"base": construct_equation_consistent_family(r, 0.0, rates=rates)}
    for eps in EPSILONS:
        label = f"eps{eps:.0e}".replace("+", "")
        family_arrays[f"{label}_plus"] = construct_equation_consistent_family(r, eps, rates=rates)
        family_arrays[f"{label}_minus"] = construct_equation_consistent_family(r, -eps, rates=rates)

    states = {}
    state_digests = {}
    projections = {}
    retained: dict[str, np.ndarray] = {}
    for rate_name, rate_values in analytic_matter_rates(r).items():
        retained[f"analytic_rate.{rate_name}"] = np.asarray(rate_values, dtype=np.float64)
    family_reports = {}
    for label, arrays in family_arrays.items():
        state = construct_production_state(arrays, r0_operators, vacuum_bssn)
        before = tangent.state_digest(state)
        projection = total.independent_total_projections(state)
        observed = tangent.observe_constraint(state)
        after = tangent.state_digest(state)
        if before != after:
            raise RuntimeError(f"state {label} mutated during zero-step projection")
        if not projection["projection_crosscheck_passed"]:
            raise RuntimeError(f"six-sector projection cross-check failed for {label}")
        if not observed["M_geom_closure_ok"] or not observed["J_closure_ok"]:
            raise RuntimeError(f"constraint/current projection accounting failed for {label}")
        states[label], state_digests[label], projections[label] = state, before, projection
        for key in STATE_ARRAY_NAMES:
            retained[f"state.{label}.{key}"] = np.asarray(arrays[key], dtype=np.float64)
        add_projection_arrays(retained, f"projection.{label}", projection)
        family_reports[label] = {
            "state_digest": before,
            "projection_crosscheck": projection["projection_crosscheck"],
            "projection_crosscheck_passed": True,
            "geometric_constraint_term_accounting_passed": True,
            "six_sector_current_accounting_passed": True,
        }

    tangent_reports = {}
    divergence_reports = {}
    eps_labels = []
    machine_eps = float(np.finfo(np.float64).eps)
    for eps in EPSILONS:
        label = f"eps{eps:.0e}".replace("+", "")
        eps_labels.append(label)
        plus_label, minus_label = f"{label}_plus", f"{label}_minus"
        base_state, plus_state, minus_state = (
            states["base"], states[plus_label], states[minus_label]
        )
        base_obs = tangent.observe_constraint(base_state)
        plus_obs = tangent.observe_constraint(plus_state)
        minus_obs = tangent.observe_constraint(minus_state)
        denom = 2.0 * eps
        mdot = (plus_obs["M_total"] - minus_obs["M_total"]) / denom
        mdot_geom = (plus_obs["M_geom"] - minus_obs["M_geom"]) / denom
        mdot_geom_terms = {
            key: (plus_obs["M_geom_pieces"][key] - minus_obs["M_geom_pieces"][key]) / denom
            for key in tangent.GEOMETRIC_M_NAMES
        }
        mdot_sector = {
            key: -8.0 * np.pi * (
                plus_obs["J_pieces"][key] - minus_obs["J_pieces"][key]
            ) / denom
            for key in tangent.SECTOR_NAMES
        }
        geom_sum = np.zeros_like(mdot_geom)
        for value in mdot_geom_terms.values():
            geom_sum += value
        reconstructed = mdot_geom.copy()
        for value in mdot_sector.values():
            reconstructed += value
        closure = mdot - reconstructed
        operand_scale = (
            np.abs(plus_obs["M_total"]) + np.abs(minus_obs["M_total"])
            + np.abs(plus_obs["M_geom"]) + np.abs(minus_obs["M_geom"])
        )
        for key in tangent.SECTOR_NAMES:
            operand_scale += 8.0 * np.pi * (
                np.abs(plus_obs["J_pieces"][key]) + np.abs(minus_obs["J_pieces"][key])
            )
        decomposition_roundoff_bound = (128.0 * machine_eps / denom) * operand_scale
        geom_operand_scale = np.zeros_like(mdot_geom)
        for key in tangent.GEOMETRIC_M_NAMES:
            geom_operand_scale += (
                np.abs(plus_obs["M_geom_pieces"][key])
                + np.abs(minus_obs["M_geom_pieces"][key])
            )
        geom_bound = (128.0 * machine_eps / denom) * geom_operand_scale
        geom_closure = mdot_geom - geom_sum
        current_closure_ok = bool(plus_obs["J_closure_ok"] and minus_obs["J_closure_ok"])
        mdot_closure_ok = bool(np.all(np.abs(closure) <= decomposition_roundoff_bound))
        geom_closure_ok = bool(np.all(np.abs(geom_closure) <= geom_bound))

        retained[f"tangent_{label}.Mdot"] = np.asarray(mdot, dtype=np.float64)
        retained[f"tangent_{label}.Mdot_geom"] = np.asarray(mdot_geom, dtype=np.float64)
        retained[f"tangent_{label}.Mdot_reconstructed"] = np.asarray(reconstructed, dtype=np.float64)
        retained[f"tangent_{label}.decomposition_gap_abs"] = np.abs(closure)
        retained[f"tangent_{label}.decomposition_roundoff_bound"] = decomposition_roundoff_bound
        retained[f"tangent_{label}.Mdot_geom_accounting_gap_abs"] = np.abs(geom_closure)
        retained[f"tangent_{label}.Mdot_geom_accounting_roundoff_bound"] = geom_bound
        for key, value in mdot_sector.items():
            retained[f"tangent_{label}.Mdot_sector.{key}"] = np.asarray(value, dtype=np.float64)
        for key, value in mdot_geom_terms.items():
            retained[f"tangent_{label}.Mdot_geom_terms.{key}"] = np.asarray(value, dtype=np.float64)
        tangent_reports[label] = {
            "epsilon": eps,
            "semantics": "analytic_equation_consistent_matter_direction_central_difference",
            "Mdot_decomposition_accounting_passed": mdot_closure_ok,
            "Mdot_decomposition_gap_max_abs": float(np.max(np.abs(closure))),
            "Mdot_decomposition_roundoff_bound_max": float(np.max(decomposition_roundoff_bound)),
            "Mdot_geom_term_accounting_passed": geom_closure_ok,
            "Mdot_geom_term_accounting_gap_max_abs": float(np.max(np.abs(geom_closure))),
            "Mdot_geom_term_accounting_roundoff_bound_max": float(np.max(geom_bound)),
            "plus_minus_six_sector_current_accounting_passed": current_closure_ok,
        }
        if not (mdot_closure_ok and geom_closure_ok and current_closure_ok):
            raise RuntimeError(f"analytic-family signed tangent accounting failed for {label}")

        sector_result = total.sector_resolved_divergence_tangent_from_states(
            base_state, plus_state, minus_state, eps
        )
        total_result = sector_result["total"]
        decomposition = sector_result["sector_decomposition"]
        for key in ("div_T_t", "div_T_r", "div_T_t_fv", "div_T_r_fv"):
            retained[f"total_T_{label}.total.{key}"] = np.asarray(total_result[key], dtype=np.float64)
        for sector, record in decomposition["sectors"].items():
            add_list_arrays(
                retained, f"total_T_{label}.sectors.{sector}",
                record["raw_signed_arrays"],
            )
            operator = record.get("operator_component_decomposition")
            if operator:
                for channel in ("pointwise", "finite_volume"):
                    for key, value in operator[channel]["raw_signed_arrays"].items():
                        retained[f"total_T_{label}.sectors.{sector}.{channel}_operator.{key}"] = np.asarray(value, dtype=np.float64)
        reconstruction = decomposition["reconstruction"]
        for channel, record in reconstruction["channels"].items():
            for key in ("sum_of_sectors", "sum_minus_total",
                        "time_density_component_predicted_gap",
                        "unexplained_gap_after_time_density_accounting"):
                retained[f"total_T_{label}.reconstruction.{channel}.{key}"] = np.asarray(record[key], dtype=np.float64)
        for key, value in decomposition["time_density_reconstruction"]["raw_signed_arrays"].items():
            retained[f"total_T_{label}.time_density_reconstruction.{key}"] = np.asarray(value, dtype=np.float64)

        divergence_reports[label] = {
            "pointwise_div_T_r_max_abs": float(np.max(np.abs(total_result["div_T_r"]))),
            "finite_volume_div_T_r_max_abs": float(np.max(np.abs(total_result["div_T_r_fv"]))),
            "sector_names": decomposition["sector_names"],
            "raw_six_sector_reconstruction_passed": reconstruction["all_sectors_reconstruct_total"],
            "time_density_component_accounting_passed": reconstruction["time_density_component_accounting_passed"],
            "six_sector_reconstruction_after_time_density_accounting": reconstruction["all_sectors_reconstruct_total_after_time_density_accounting"],
            "raw_sum_max_gap": reconstruction["raw_sum_max_gap"],
            "unexplained_max_gap_after_time_density_accounting": reconstruction["unexplained_max_gap"],
            "channel_accounting": {
                key: {
                    "max_abs_gap": value["max_abs_gap"],
                    "max_abs_unexplained_gap": value["max_abs_unexplained_gap"],
                    "passed": value["passed"],
                    "component_corrected_passed": value["component_corrected_passed"],
                }
                for key, value in reconstruction["channels"].items()
            },
            "pointwise_operator_component_reconstruction_passed": all(
                decomposition["sectors"][name].get("operator_component_decomposition", {}).get(
                    "component_reconstruction_passed", True
                ) for name in ("S", "D", "phi")
            ),
            "state_pairing": decomposition["state_pairing"],
        }

    base_digest_after = tangent.state_digest(states["base"])
    if base_digest_after != state_digests["base"]:
        raise RuntimeError("base frozen state mutated during manufactured fixture construction")

    # Preserve the full state and witness inventory in accepted-fixture style.
    records = []
    for key in sorted(retained):
        value = np.ascontiguousarray(np.asarray(retained[key], dtype=np.float64))
        if not np.all(np.isfinite(value)):
            raise FloatingPointError(f"non-finite fixture array {key}")
        records.append({
            "path": key,
            "key": key,
            "shape": list(value.shape),
            "dtype": str(value.dtype),
            "sha256": digest_bytes(value.tobytes(order="C")),
        })

    code_blobs = {path: git(ROOT, "rev-parse", f"HEAD:{path}") for path in CODE_PATHS}
    report = {
        "schema": "issue4_equation_consistent_matter_direction_fixture_v1",
        "status": "completed_manufactured_fixture_construction",
        "state_class": "equation_consistent_matter_direction",
        "diagnostic_only": True,
        "resolution": N,
        "r_max": R_MAX,
        "dr": R_MAX / N,
        "grid_centers": "r_i=(i+1/2)*dr",
        "profile": PROFILE,
        "family_parameterization": {
            "direction": "instantaneous matter-equation direction at s=0",
            "geometry_tangent": "frozen_flat_background_not_Einstein_evolution",
            "scalar_spatial_operator": "pinned_vendor_fourth_order_cell_derivative_and_second_derivative",
            "parameter": "s",
            "U(0)": "base",
            "U(+epsilon)": "plus member independently constructed by analytic formulas",
            "U(-epsilon)": "minus member independently constructed by analytic formulas",
            "held_fixed": [
                "geometry.*", "scalars.S", "scalars.D", "scalars.phi",
                "primitive fluid rho_species and pressure_species profiles",
            ],
            "changed": [
                "scalars.PS, PD, Pi from pinned flat-slice scalar accelerations",
                "dark-matter velocity from Q_DM,r=beta_DM*rho_DM*partial_r(phi)",
                "radiation velocity from -partial_r(p)/(rho+p)",
                "baryon velocity remains zero",
                "conserved fluid arrays reconstructed from these primitives",
            ],
            "construction_method": "initializer_free_analytic_profile_constructor",
        },
        "parity": {
            "centre_regularity": {
                "even_fields": ["a", "b", "X", "alpha", "Aa", "K", "Lambda", "B",
                                "S", "D", "phi", "PS", "PD", "Pi", "rho_species",
                                "pressure_species", "conserved rest", "conserved energy_t"],
                "odd_fields": ["fluid radial velocity", "conserved momentum_r", "C^r"],
                "inner_face_flux": "zero for regular odd radial fluxes",
            },
            "geometry": "even",
            "scalar_fields": "even",
            "scalar_conjugate_momenta": "even profiles with signed family parameter",
            "fluid_rest_density": "even",
            "fluid_radial_velocity": "odd regular r*exp(-r^2/w^2)",
            "connection_difference_Cr": "odd; C^r=Lambda-Lambda_metric=0 on exact flat geometry",
        },
        "electromagnetic_carrier": None,
        "em_on_content": False,
        "sector_names": ["dark_matter", "baryons", "radiation", "S", "D", "phi"],
        "sector_projection_fields": ["rho", "pr", "pt", "j"],
        "directional_epsilons": list(EPSILONS),
        "tangent_semantics": "analytic_equation_consistent_matter_direction_central_difference",
        "tangent_definition": "Mdot_manufactured(epsilon)=[M_total(U(+epsilon))-M_total(U(-epsilon))]/(2*epsilon)",
        "tangent_is_production_rhs_tangent": False,
        "production_steps_run": 0,
        "state_advanced": False,
        "production_initializer_called": False,
        "production_rhs_called": False,
        "vendor_step_called": False,
        "scientific_admission": False,
        "acceptance_oracle_satisfied": False,
        "conservation_closure_admitted": False,
        "physical_exchange_law_admitted": False,
        "admission_flags": {
            "scientific_admission": False,
            "acceptance_oracle_satisfied": False,
            "conservation_closure_admitted": False,
            "physical_exchange_law_admitted": False,
        },
        "state_digests": state_digests,
        "state_unchanged_after_construction": True,
        "family_projection_reports": family_reports,
        "tangent_reports": tangent_reports,
        "divergence_reports": divergence_reports,
        "pins": {
            "diagnostic_source_anchor": diag_commit,
            "helper_baseline": baseline_commit,
            "vendor": vendor_commit,
            "tangent_helper_blob": tangent_blob,
            "total_stress_helper_blob": total_blob,
            "reviewed_code_blobs": code_blobs,
        },
        "raw_arrays": {"array_count": len(records), "array_records": records},
        "array_names": sorted(retained),
        "arithmetic_witnesses": {
            "all_family_projections_match_canonical": all(
                row["projection_crosscheck_passed"] for row in family_reports.values()
            ),
            "all_tangent_accounting_passed": all(
                row["Mdot_decomposition_accounting_passed"] and row["Mdot_geom_term_accounting_passed"]
                and row["plus_minus_six_sector_current_accounting_passed"]
                for row in tangent_reports.values()
            ),
            "all_pointwise_operator_component_reconstructions_passed": all(
                row["pointwise_operator_component_reconstruction_passed"]
                for row in divergence_reports.values()
            ),
        },
        "analytic_direction": {
            "rates": "analytic_rate.* arrays in the archive",
            "fluid_continuity": "all primitive density time derivatives vanish at the zero-velocity base",
            "dark_matter_force": "beta_DM*partial_r(phi), beta_DM=-0.04",
            "radiation_force": "-partial_r(rho/3)/(4*rho/3)",
            "geometry": "held fixed; no Einstein geometric tangent is supplied",
            "scalar_operator": "pinned vendor fourth-order cell derivative and second derivative",
            "matter_gradient_operator": "pinned vendor fourth-order cell derivative",
        },
        "interpretation_limits": [
            "The matter-direction tangent is an analytic-family central difference, not a production-RHS tangent.",
            "The constructor calls no production initializer, production RHS, kernel step, or vendor evolution step. Geometry is intentionally frozen and the scalar direction is continuum-analytic.",
            "Pointwise and finite-volume total-stress divergences are independently evaluated from the same frozen analytic family triplets by the pinned helper.",
            "Six-sector raw sums and time-density component accounting are diagnostic witnesses, not physical conservation or scientific admission.",
        ],
    }

    arrays_path = OUT / "tangent_arrays.npz"
    np.savez_compressed(arrays_path, **retained)
    report_path = OUT / "report.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    sums_path = OUT / "SHA256SUMS.txt"
    sums_path.write_text(
        f"{digest_file(report_path)}  report.json\n{digest_file(arrays_path)}  tangent_arrays.npz\n"
    )
    fixture_zip = OUT / "issue4_equation_consistent_matter_direction_N60.zip"
    with zipfile.ZipFile(fixture_zip, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        archive.write(sums_path, arcname="SHA256SUMS.txt")
        archive.write(report_path, arcname="report.json")
        archive.write(arrays_path, arcname="tangent_arrays.npz")
    print(json.dumps({
        "status": "PASS",
        "state_class": report["state_class"],
        "tangent_semantics": report["tangent_semantics"],
        "family_projection_count": len(family_reports),
        "tangent": tangent_reports,
        "divergence": divergence_reports,
        "array_count": len(retained),
        "fixture_zip": fixture_zip.name,
        "fixture_zip_sha256": digest_file(fixture_zip),
        "production_initializer_called": False,
        "production_rhs_called": False,
        "vendor_step_called": False,
        "production_steps_run": 0,
    }, indent=2))


if __name__ == "__main__":
    main()
