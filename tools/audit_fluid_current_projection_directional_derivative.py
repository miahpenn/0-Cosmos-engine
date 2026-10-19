#!/usr/bin/env python3
"""Projection-map directional derivative diagnostic; not a propagation closure oracle."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from tools import audit_current_source_momentum_constraint_tangent as audit

EXPECTED_INITIAL = "a256b57f274ee120b55c0f805ed488d91382a9799647c8f9811d7df82aa337c0"
EXPECTED_ACCEPTED = "c6acb3711d59a85c804a88208c9bc0356e8227a806c0684ae3a8afac72c9091a"
RESOLUTION = 60
R_MAX = 80.0
DT = 0.015
EPSILONS = (8.0e-5, 4.0e-5, 2.0e-5, 1.0e-5)
SPECIES = (
    ("dark_matter", audit.Species.DARK_MATTER),
    ("baryons", audit.Species.BARYON),
    ("radiation", audit.Species.RADIATION),
)


def project_fluid_currents(state) -> dict[str, np.ndarray]:
    """Evaluate only the per-species fluid projection map, not the constraint."""
    metric = audit.metric_slice_from_q(state.grid, state.geometry)
    out = {}
    for name, species in SPECIES:
        projected = audit.project_species(metric, getattr(state.matter, name), species)
        value = np.asarray(projected["j"], dtype=np.float64)
        if value.shape != np.asarray(state.grid.centers).shape or not np.all(np.isfinite(value)):
            raise FloatingPointError(f"invalid projected current for {name}")
        out[name] = value
    return out


def perturbed_partial_state(state, rates, delta: float, *, metric_direction: bool):
    """Build a partial projection-map direction without evaluating M_r.

    metric_direction=True varies only spatial metric fields a,b,X. Otherwise,
    only the conservative fluid arrays vary. This is a partial-Jacobian control;
    neither control is mislabeled as the full accepted RHS direction.
    """
    geom = state.geometry.copy()
    matter = state.matter.copy()
    if metric_direction:
        for name in ("a", "b", "X"):
            setattr(
                geom,
                name,
                np.asarray(getattr(state.geometry, name), dtype=np.float64)
                + delta * np.asarray(rates["geometry"][name], dtype=np.float64),
            )
    else:
        for species_name, _species in SPECIES:
            base = getattr(state.matter, species_name)
            rate = rates["matter"][species_name]
            updated = type(base)(
                *(
                    np.asarray(getattr(base, field), dtype=np.float64)
                    + delta * np.asarray(rate[field], dtype=np.float64)
                    for field in ("rest", "energy_t", "momentum_r")
                )
            )
            setattr(matter, species_name, updated)
    return audit.ProductionState(
        grid=state.grid,
        geometry=geom,
        scalars=state.scalars,
        matter=matter,
        t=state.t,
        tau=state.tau,
        e_folds=state.e_folds,
    )


def derivative_sequence(state, rates, kernel, direction: str) -> dict:
    rows = []
    for eps in EPSILONS:
        if direction == "full_rhs_with_lapse_resolve":
            plus_state = audit.shifted_state(kernel, state, rates, eps)
            minus_state = audit.shifted_state(kernel, state, rates, -eps)
        else:
            is_metric = direction == "metric_only"
            plus_state = perturbed_partial_state(state, rates, eps, metric_direction=is_metric)
            minus_state = perturbed_partial_state(state, rates, -eps, metric_direction=is_metric)
        plus = project_fluid_currents(plus_state)
        minus = project_fluid_currents(minus_state)
        rows.append({
            "epsilon": eps,
            "by_species": {
                name: (plus[name] - minus[name]) / (2.0 * eps)
                for name, _species in SPECIES
            },
        })

    summaries = {}
    for name, _species in SPECIES:
        summaries[name] = []
        for row in rows:
            arr = row["by_species"][name]
            summaries[name].append({
                "epsilon": row["epsilon"],
                "cell_0": float(arr[0]),
                "cell_1": float(arr[1]),
                "max_abs": float(np.max(np.abs(arr))),
                "rms": float(np.sqrt(np.mean(arr * arr))),
                "all_finite": bool(np.all(np.isfinite(arr))),
            })
    convergence = {}
    for name, _species in SPECIES:
        coarse = rows[-2]["by_species"][name]
        fine = rows[-1]["by_species"][name]
        scale = max(float(np.max(np.abs(fine))), float(np.max(np.abs(coarse))), 1.0e-300)
        convergence[name] = {
            "last_two_epsilon_max_abs_difference": float(np.max(np.abs(fine - coarse))),
            "last_two_epsilon_relative_difference": float(np.max(np.abs(fine - coarse)) / scale),
        }
    return {
        "epsilon_sequence": list(EPSILONS),
        "per_species": summaries,
        "last_two_step_convergence": convergence,
        "all_finite": all(
            item["all_finite"]
            for per_species in summaries.values()
            for item in per_species
        ),
    }


def run() -> dict:
    provenance = audit.provenance()
    kernel = audit.V55ProductionKernel()
    initial = kernel.initialize(
        resolution=RESOLUTION, r_max=R_MAX, include_radiation=True
    )
    initial_digest = audit.state_digest(initial)
    if initial_digest != EXPECTED_INITIAL:
        raise RuntimeError(f"initial digest mismatch: {initial_digest}")
    accepted = kernel.step(initial, DT)
    accepted_digest_before = audit.state_digest(accepted)
    if accepted_digest_before != EXPECTED_ACCEPTED:
        raise RuntimeError(f"accepted digest mismatch: {accepted_digest_before}")
    rates = audit.rhs_direction(kernel, accepted)
    before = audit.state_digest(accepted)

    results = {
        "metric_only": derivative_sequence(accepted, rates, kernel, "metric_only"),
        "conserved_state_only": derivative_sequence(accepted, rates, kernel, "conserved_state_only"),
        "full_rhs_with_lapse_resolve": derivative_sequence(
            accepted, rates, kernel, "full_rhs_with_lapse_resolve"
        ),
    }
    after = audit.state_digest(accepted)
    if before != after:
        raise RuntimeError("projection diagnostic mutated the accepted state")
    all_finite = all(item["all_finite"] for item in results.values())
    return {
        "schema": "issue4_fluid_projection_directional_derivative_v1",
        "status": "projection_diagnostic_completed" if all_finite else "projection_diagnostic_nonfinite",
        "provenance": provenance,
        "configuration": {
            "resolution": RESOLUTION,
            "r_max": R_MAX,
            "dt": DT,
            "epsilon_sequence": list(EPSILONS),
            "production_steps_advanced": 1,
        },
        "state_identity": {
            "initial_sha256": initial_digest,
            "accepted_expected_sha256": EXPECTED_ACCEPTED,
            "accepted_before_sha256": before,
            "accepted_after_sha256": after,
            "accepted_state_unchanged": before == after,
        },
        "directions": results,
        "gates": {
            "projection_arrays_finite": all_finite,
            "projection_step_size_plateau_established": False,
            "independent_full_propagation_identity_closed": False,
            "scientific_admission": False,
            "acceptance_oracle_satisfied": False,
        },
        "interpretation_limit": (
            "This is a projection-map derivative diagnostic. The last-two-step "
            "difference is reported, not treated as a passing plateau without a "
            "predeclared numerical bound. It does not establish full propagation closure."
        ),
    }


def main() -> int:
    result = run()
    out = Path("results/accepted-state-raw-momentum-witness/fluid-current-projection-directional-derivative")
    out.mkdir(parents=True, exist_ok=True)
    (out / "report.json").write_text(
        json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "status": result["status"],
        "state_identity": result["state_identity"],
        "gates": result["gates"],
        "last_two_step_convergence": {
            direction: item["last_two_step_convergence"]
            for direction, item in result["directions"].items()
        },
    }, sort_keys=True), flush=True)
    return 0 if result["status"] == "projection_diagnostic_completed" else 2


if __name__ == "__main__":
    raise SystemExit(main())
