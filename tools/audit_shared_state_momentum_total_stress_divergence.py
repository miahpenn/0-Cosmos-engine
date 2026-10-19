#!/usr/bin/env python3
"""Shared accepted-state momentum / total-stress divergence witness; diagnostic-only."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from tools import audit_current_source_momentum_constraint_tangent as momentum
from tools import audit_total_stress_energy_divergence as stress
from tools import audit_brown_vacuum_momentum_propagation as brown

INITIAL_SHA256 = "a256b57f274ee120b55c0f805ed488d91382a9799647c8f9811d7df82aa337c0"
ACCEPTED_SHA256 = "c6acb3711d59a85c804a88208c9bc0356e8227a806c0684ae3a8afac72c9091a"
RESOLUTION, R_MAX, DT, EPSILON = 60, 80.0, 0.015, 1.0e-5


def summary(x):
    a = np.asarray(x, dtype=np.float64)
    if not np.all(np.isfinite(a)):
        raise FloatingPointError("non-finite value in shared-state calculation")
    i = int(np.argmax(np.abs(a)))
    return {"max_abs": float(abs(a[i])), "signed_at_max": float(a[i]),
            "index": i, "rms": float(np.sqrt(np.mean(a*a))),
            "cell_0": float(a[0]), "cell_1": float(a[1])}


def require_digest(state, expected: str, label: str) -> str:
    actual = momentum.state_digest(state)
    if actual != expected:
        raise RuntimeError(f"{label} digest mismatch: expected {expected}, got {actual}")
    return actual


def run(output_dir: Path):
    source = momentum.provenance()
    stress_source = stress.verify_provenance()
    kernel = momentum.V55ProductionKernel()
    initial = kernel.initialize(resolution=RESOLUTION, r_max=R_MAX, include_radiation=True)
    initial_digest = require_digest(initial, INITIAL_SHA256, "initial-state")
    accepted = kernel.step(initial, DT)
    accepted_digest = require_digest(accepted, ACCEPTED_SHA256, "accepted-state")
    rates = momentum.rhs_direction(kernel, accepted)
    plus = momentum.shifted_state(kernel, accepted, rates, EPSILON)
    minus = momentum.shifted_state(kernel, accepted, rates, -EPSILON)
    base_digest = momentum.state_digest(accepted)
    plus_digest = momentum.state_digest(plus)
    minus_digest = momentum.state_digest(minus)
    require_digest(accepted, ACCEPTED_SHA256, "shared base state")
    # Negative control: the exact same state must be rejected against a wrong digest.
    try:
        require_digest(accepted, "0" * 64, "negative-control")
    except RuntimeError:
        wrong_digest_rejected = True
    else:
        wrong_digest_rejected = False
    if not wrong_digest_rejected:
        raise RuntimeError("negative control failed: mismatched state digest was accepted")
    if plus_digest == base_digest or minus_digest == base_digest:
        raise RuntimeError("directional states unexpectedly equal the accepted base state")
    if not np.all(np.asarray(accepted.geometry.beta) == 0.0):
        raise RuntimeError("observer requires zero shift on the accepted slice")

    # Both observers consume the same accepted/base/plus/minus objects.
    base_constraint = momentum.observe_constraint(accepted)
    plus_constraint = momentum.observe_constraint(plus)
    minus_constraint = momentum.observe_constraint(minus)
    mdot_total = (np.asarray(plus_constraint["M_total"]) -
                  np.asarray(minus_constraint["M_total"])) / (2*EPSILON)
    mdot_geom = (np.asarray(plus_constraint["M_geom"]) -
                 np.asarray(minus_constraint["M_geom"])) / (2*EPSILON)
    div = stress.divergence_tangent_from_states(accepted, plus, minus, EPSILON)
    sector = stress.sector_resolved_divergence_tangent_from_states(accepted, plus, minus, EPSILON)

    # Independent Brown Eq. 11b structural predictor evaluated on this same
    # accepted geometry. It is intentionally vacuum-form: the residual exposes
    # omitted matter/source and discrete product-rule contributions; it is not
    # labelled as a complete matter-coupled propagation RHS.
    _, vacuum, _ = momentum.adapter.vendor_modules()
    geometric_constraints = vacuum.constraints(accepted.grid, accepted.geometry)
    brown_terms, brown_prediction = brown.brown_lower_momentum_rhs(
        accepted.grid, accepted.geometry, geometric_constraints
    )
    brown_residual = mdot_geom - brown_prediction
    hgrad = np.asarray(brown_terms["alpha_DrH_over_6"], dtype=np.float64)
    wrong_sign_prediction = brown_prediction - 2.0*hgrad
    wrong_sign_delta = brown_prediction - wrong_sign_prediction
    operand_scale = sum((np.abs(v) for v in brown_terms.values()),
                        np.zeros_like(brown_prediction))
    negative_control_bound = 128.0 * np.finfo(np.float64).eps * operand_scale
    negative_control_detected = bool(np.any(np.abs(wrong_sign_delta) > negative_control_bound))
    if not negative_control_detected:
        raise RuntimeError("Brown H-gradient wrong-sign negative control was not detected")

    after_digest = momentum.state_digest(accepted)
    if after_digest != base_digest:
        raise RuntimeError("shared accepted state mutated during observation")

    r = np.asarray(accepted.grid.centers, dtype=np.float64)
    output_dir.mkdir(parents=True, exist_ok=True)
    arrays = {
        "radius": r,
        "Mdot_total_signed": mdot_total,
        "Mdot_geometric_signed": mdot_geom,
        "brown_vacuum_structural_prediction_signed": np.asarray(brown_prediction),
        "brown_structural_residual_signed": np.asarray(brown_residual),
        "brown_H_gradient_term_signed": hgrad,
        "brown_wrong_sign_mutation_signed": wrong_sign_delta,
        "div_T_t_signed": np.asarray(div["div_T_t"]),
        "div_T_r_signed": np.asarray(div["div_T_r"]),
        "div_T_t_fv_signed": np.asarray(div["div_T_t_fv"]),
        "div_T_r_fv_signed": np.asarray(div["div_T_r_fv"]),
    }
    np.savez_compressed(output_dir / "shared_state_signed_arrays.npz", **arrays)
    report = {
        "schema": "issue4_shared_state_momentum_total_stress_divergence_v2",
        "status": "completed_diagnostic_only",
        "provenance": {"momentum_observer": source, "total_stress_observer": stress_source},
        "configuration": {"resolution": RESOLUTION, "r_max": R_MAX, "dt": DT,
                          "epsilon": EPSILON, "production_steps_advanced": 1},
        "state_pairing": {"initial_sha256": initial_digest, "accepted_sha256": base_digest,
                          "plus_sha256": plus_digest, "minus_sha256": minus_digest,
                          "accepted_after_observation_sha256": after_digest,
                          "accepted_unchanged": after_digest == base_digest,
                          "all_observers_used_same_state_objects": True,
                          "wrong_digest_negative_control_rejected": wrong_digest_rejected},
        "signed_residuals": {
            "Mdot_total_constraint": summary(mdot_total),
            "Mdot_geometric_constraint": summary(mdot_geom),
            "Brown_vacuum_structural_prediction": summary(brown_prediction),
            "Brown_structural_residual_measured_minus_predicted": summary(brown_residual),
            "Brown_H_gradient_term": summary(hgrad),
            "Brown_wrong_sign_mutation": summary(wrong_sign_delta),
            "div_T_t_cell_centered": summary(div["div_T_t"]),
            "div_T_r_cell_centered": summary(div["div_T_r"]),
            "div_T_t_finite_volume": summary(div["div_T_t_fv"]),
            "div_T_r_finite_volume": summary(div["div_T_r_fv"]),
            "sector_decomposition": sector["sector_decomposition"]["reconstruction"],
        },
        "brown_structural_crosscheck": {
            "formula": "beta*D_r(M)+(D_r beta)*M+alpha*K*M-(D_r alpha)*H/3"
                       "+alpha*D_r(H)/6+2*alpha*chi*D_r^2(C)/3"
                       "+[2*(D_r alpha)*chi/3-alpha*(D_r chi)/3"
                       "+alpha*D_r ln(g_theta_theta)]*D_r(C)",
            "hamiltonian_slot": "independent metric-only warped-product Ricci reconstruction",
            "wrong_sign_H_gradient_control_detected": negative_control_detected,
            "scope": "vacuum structural terms only; matter-source and discrete residual remain explicit",
            "full_matter_coupled_identity": False,
        },
        "raw_arrays": {"file": "shared_state_signed_arrays.npz",
                       "names": sorted(arrays), "all_finite": True},
        "gates": {"shared_state_digest_match": True,
                  "wrong_digest_negative_control_passed": wrong_digest_rejected,
                  "observer_immutability_passed": True,
                  "brown_wrong_sign_negative_control_passed": negative_control_detected,
                  "independent_full_propagation_identity_closed": False,
                  "scientific_admission": False, "acceptance_oracle_satisfied": False},
        "interpretation": (
            "Shared-state extension: signed measured geometric and total momentum tangents, "
            "independent total-stress divergence, and Brown vacuum structural terms are now "
            "evaluated on the same accepted state and directional pair. The Brown residual is "
            "not a closure claim: matter-source terms and discrete operator/product-rule terms "
            "remain to be independently separated."
        )
    }
    (output_dir / "shared_state_report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output-dir", type=Path,
                   default=Path("results/accepted-state-raw-momentum-witness"))
    args = p.parse_args()
    report = run(args.output_dir)
    print(json.dumps({"status": report["status"], "state_pairing": report["state_pairing"],
                      "signed_residuals": {k:v for k,v in report["signed_residuals"].items()
                                           if k != "sector_decomposition"},
                      "brown_structural_crosscheck": report["brown_structural_crosscheck"],
                      "gates": report["gates"]}, sort_keys=True), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
