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

def run(output_dir: Path):
    source = momentum.provenance()
    stress_source = stress.verify_provenance()
    kernel = momentum.V55ProductionKernel()
    initial = kernel.initialize(resolution=RESOLUTION, r_max=R_MAX, include_radiation=True)
    initial_digest = momentum.state_digest(initial)
    if initial_digest != INITIAL_SHA256:
        raise RuntimeError(f"initial-state digest mismatch: {initial_digest}")
    accepted = kernel.step(initial, DT)
    accepted_digest = momentum.state_digest(accepted)
    if accepted_digest != ACCEPTED_SHA256:
        raise RuntimeError(f"accepted-state digest mismatch: {accepted_digest}")
    rates = momentum.rhs_direction(kernel, accepted)
    plus = momentum.shifted_state(kernel, accepted, rates, EPSILON)
    minus = momentum.shifted_state(kernel, accepted, rates, -EPSILON)
    base_digest = momentum.state_digest(accepted)
    plus_digest = momentum.state_digest(plus)
    minus_digest = momentum.state_digest(minus)
    if base_digest != ACCEPTED_SHA256:
        raise RuntimeError("shared base state digest mismatch")
    if plus_digest == base_digest or minus_digest == base_digest:
        raise RuntimeError("directional states unexpectedly equal the accepted base state")
    if not np.all(np.asarray(accepted.geometry.beta) == 0.0):
        raise RuntimeError("observer requires zero shift on the accepted slice")

    # Both independent observers consume these exact objects; neither may reinitialize.
    base_constraint = momentum.observe_constraint(accepted)
    plus_constraint = momentum.observe_constraint(plus)
    minus_constraint = momentum.observe_constraint(minus)
    mdot = (np.asarray(plus_constraint["M_total"]) - np.asarray(minus_constraint["M_total"])) / (2*EPSILON)
    div = stress.divergence_tangent_from_states(accepted, plus, minus, EPSILON)
    sector = stress.sector_resolved_divergence_tangent_from_states(accepted, plus, minus, EPSILON)
    after_digest = momentum.state_digest(accepted)
    if after_digest != base_digest:
        raise RuntimeError("shared accepted state mutated during observation")

    r = np.asarray(accepted.grid.centers, dtype=np.float64)
    output_dir.mkdir(parents=True, exist_ok=True)
    arrays = {
        "radius": r,
        "Mdot_signed": mdot,
        "div_T_t_signed": np.asarray(div["div_T_t"]),
        "div_T_r_signed": np.asarray(div["div_T_r"]),
        "div_T_t_fv_signed": np.asarray(div["div_T_t_fv"]),
        "div_T_r_fv_signed": np.asarray(div["div_T_r_fv"]),
    }
    np.savez_compressed(output_dir / "shared_state_signed_arrays.npz", **arrays)
    report = {
        "schema": "issue4_shared_state_momentum_total_stress_divergence_v1",
        "status": "completed_diagnostic_only",
        "provenance": {"momentum_observer": source, "total_stress_observer": stress_source},
        "configuration": {"resolution": RESOLUTION, "r_max": R_MAX, "dt": DT,
                          "epsilon": EPSILON, "production_steps_advanced": 1},
        "state_pairing": {"initial_sha256": initial_digest, "accepted_sha256": base_digest,
                          "plus_sha256": plus_digest, "minus_sha256": minus_digest,
                          "accepted_after_observation_sha256": after_digest,
                          "accepted_unchanged": after_digest == base_digest,
                          "all_observers_used_same_state_objects": True},
        "signed_residuals": {
            "Mdot_total_constraint": summary(mdot),
            "div_T_t_cell_centered": summary(div["div_T_t"]),
            "div_T_r_cell_centered": summary(div["div_T_r"]),
            "div_T_t_finite_volume": summary(div["div_T_t_fv"]),
            "div_T_r_finite_volume": summary(div["div_T_r_fv"]),
            "sector_decomposition": sector["sector_decomposition"]["reconstruction"],
        },
        "raw_arrays": {"file": "shared_state_signed_arrays.npz",
                       "names": sorted(arrays), "all_finite": True},
        "gates": {"shared_state_digest_match": True,
                  "observer_immutability_passed": True,
                  "independent_full_propagation_identity_closed": False,
                  "scientific_admission": False, "acceptance_oracle_satisfied": False},
        "interpretation": ("First same-state calculation: the accepted production state and one "
            "shared plus/minus RHS-direction pair are passed to both momentum-constraint and "
            "independent total-stress divergence observers. Signed profiles are recorded, not "
            "asserted to close. This is diagnostic evidence, not a propagation-identity pass.")
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
                      "gates": report["gates"]}, sort_keys=True), flush=True)
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
