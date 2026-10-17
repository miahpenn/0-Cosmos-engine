"""Observation-only explicit-RHS shadow on CMC-aligned copied stage slices.

The stage-transition ordering is pinned by test_native_cmc_stage_contract_preflight.py.
This witness builds only stage 0, the explicit predictor, and the primary
predictor from copies using the existing production-owned CMC/RHS/boundary
helpers. It never calls V55ProductionKernel.step. The comparison is restricted
to the explicit metric fields a/b/X; alpha/beta RHSs belong to the moving-
puncture gauge and are not compared against production CMC gauge evolution.

PR #340 remains the initialized frozen-slice/source decomposition baseline.
This is one fixed-grid stage observation, not a repeated resolution matrix.
"""
from __future__ import annotations

import hashlib
import inspect
import json

import numpy as np

from engine.matter_system import Species, validate_species_admissibility
from engine.native_geometry import (
    NativeSphericalCellGrid,
    native_bssn_geometry_terms,
    native_explicit_rhs,
)
from engine.production_kernel import ProductionState, V55ProductionKernel
from engine.scalar_system import ScalarFields
from engine.v55_matter import V55MatterState, metric_slice_from_q


_GEOMETRY_FIELDS = ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B")
_SCALAR_FIELDS = ("S", "PS", "D", "PD", "phi", "Pi")
_MATTER_SPECIES = ("dark_matter", "baryons", "radiation", "charged_dust")
_MATTER_FIELDS = ("rest", "energy_t", "momentum_r")
_EXPLICIT_METRIC_FIELDS = ("a", "b", "X")
_DT = 1.0e-3


def _slice_sha256(geometry, scalars, matter):
    digest = hashlib.sha256()
    for label, obj, fields in (
        ("geometry", geometry, _GEOMETRY_FIELDS),
        ("scalars", scalars, _SCALAR_FIELDS),
    ):
        for name in fields:
            values = np.ascontiguousarray(np.asarray(getattr(obj, name)))
            digest.update(label.encode("ascii") + b"\0" + name.encode("ascii") + b"\0")
            digest.update(values.dtype.str.encode("ascii") + b"\0")
            digest.update(repr(values.shape).encode("ascii") + b"\0")
            digest.update(values.tobytes(order="C"))

    for species_name in _MATTER_SPECIES:
        species = getattr(matter, species_name, None)
        if species is None:
            continue
        for name in _MATTER_FIELDS:
            values = np.ascontiguousarray(np.asarray(getattr(species, name)))
            digest.update(b"matter\0" + species_name.encode("ascii") + b"\0")
            digest.update(name.encode("ascii") + b"\0")
            digest.update(values.dtype.str.encode("ascii") + b"\0")
            digest.update(repr(values.shape).encode("ascii") + b"\0")
            digest.update(values.tobytes(order="C"))
    return digest.hexdigest()


def _regional_rhs_residuals(native_rhs, adapter_rhs, geometry, native_grid):
    areal = np.asarray(
        native_bssn_geometry_terms(native_grid, geometry)["areal_radius"],
        dtype=float,
    )
    masks = {
        "centre_R_areal_lt_2": areal < 2.0,
        "fixed_annulus_2_le_R_areal_le_6": (areal >= 2.0) & (areal <= 6.0),
        "outer_R_areal_gt_6": areal > 6.0,
    }
    summaries = {}
    for region, mask in masks.items():
        selected_count = int(np.count_nonzero(mask))
        assert selected_count >= 4, (
            f"{region}: expected at least four descriptive cells, got {selected_count}"
        )
        fields = {}
        for name in _EXPLICIT_METRIC_FIELDS:
            native = np.asarray(native_rhs[name], dtype=float)
            reference = np.asarray(adapter_rhs[name], dtype=float)
            assert native.shape == reference.shape == areal.shape, (
                f"{region}/{name}: incompatible profile shapes"
            )
            residual = native - reference
            assert np.all(np.isfinite(residual)), (
                f"{region}/{name}: non-finite native-minus-adapter residual"
            )
            selected = residual[mask]
            fields[name] = {
                "signed_mean_native_minus_adapter": float(np.mean(selected)),
                "rms_native_minus_adapter": float(np.sqrt(np.mean(selected**2))),
                "max_abs_native_minus_adapter": float(np.max(np.abs(selected))),
                "bitwise_equal_on_region": bool(np.array_equal(native[mask], reference[mask])),
            }
        summaries[region] = {"cells": selected_count, "fields": fields}
    return summaries


def _observe_stage(kernel, native_grid, stage, *, evaluate_coupled_rhs):
    geometry = stage.geometry
    scalars = stage.scalars
    matter = stage.matter
    digest_before = _slice_sha256(geometry, scalars, matter)

    scalar_rhs = None
    matter_rhs = None
    if evaluate_coupled_rhs:
        scalar_rhs, matter_rhs, _ = kernel._rhs(stage)

    blocks = kernel._geometry_stage_terms(
        stage.grid, geometry, scalars, matter, stage.electromagnetic
    )
    native_rhs = native_explicit_rhs(native_grid, geometry)
    digest_after = _slice_sha256(geometry, scalars, matter)
    assert digest_after == digest_before, (
        f"{stage.geometry!r}: production helper or RHS evaluation mutated its "
        "copied input slice"
    )

    record = {
        "stage": stage.stage_name if hasattr(stage, "stage_name") else "unnamed",
        "geometry_input_sha256_before": digest_before,
        "geometry_input_sha256_after": digest_after,
        "native_minus_adapter_explicit_metric_rhs": _regional_rhs_residuals(
            native_rhs, blocks["explicit"], geometry, native_grid
        ),
        "compared_fields": list(_EXPLICIT_METRIC_FIELDS),
        "excluded_gauge_rhs_fields": {
            "alpha": "production CMC lapse is re-solved, not advanced by moving-puncture 1+log RHS",
            "beta": "production CMC branch pins beta/B to zero instead of advancing moving-puncture shift",
        },
    }
    return blocks, scalar_rhs, matter_rhs, record


def test_native_explicit_rhs_on_cmc_stage_slices_without_calling_full_step(capsys):
    kernel = V55ProductionKernel()
    accepted = kernel.initialize(
        resolution=80, r_max=40.0, include_radiation=True
    )
    accepted_hash_before = _slice_sha256(
        accepted.geometry, accepted.scalars, accepted.matter
    )

    grid = accepted.grid
    native_grid = NativeSphericalCellGrid(grid.n, grid.r_max)
    np.testing.assert_allclose(
        native_grid.centers, grid.centers, rtol=0.0, atol=1.0e-14
    )

    # Stage 0: only copies are refreshed onto the current CMC gauge.
    g0 = accepted.geometry.copy()
    s0 = accepted.scalars.copy()
    m0 = accepted.matter.copy()
    em0 = accepted.electromagnetic
    validate_species_admissibility(
        metric_slice_from_q(grid, g0), m0.radiation, Species.RADIATION,
        stage="initial",
    )
    g0.alpha = kernel._solve_lapse_for_state(
        grid, g0, s0, m0, em0
    )[0]
    g0.beta.fill(0.0)
    g0.B.fill(0.0)
    validate_species_admissibility(
        metric_slice_from_q(grid, g0), m0.radiation, Species.RADIATION,
        stage="initial",
    )

    stage0 = ProductionState(
        grid=grid,
        geometry=g0,
        scalars=s0,
        matter=m0,
        electromagnetic=em0,
        t=accepted.t,
        tau=accepted.tau,
        e_folds=accepted.e_folds,
    )
    blocks0, srhs0, mrhs0, record0 = _observe_stage(
        kernel, native_grid, stage0, evaluate_coupled_rhs=True
    )
    record0["stage"] = "stage0_cmc"

    # Explicit predictor: copy production's source-locked a/b/X update and
    # scalar/matter predictor, then refresh the lapse on that CMC slice.
    g_explicit1 = g0.copy()
    for name in ("a", "b", "X"):
        setattr(
            g_explicit1,
            name,
            getattr(g0, name) + _DT * blocks0["explicit"][name],
        )
    g_explicit1.beta.fill(0.0)
    g_explicit1.B.fill(0.0)

    s1 = ScalarFields(*(
        getattr(s0, name) + _DT * getattr(srhs0, name)
        for name in _SCALAR_FIELDS
    ))
    m1_base = tuple(
        kernel._add_matter(getattr(m0, name), mrhs0[name], _DT)
        for name in ("dark_matter", "baryons", "radiation")
    )
    m1_charged = (
        None if m0.charged_dust is None
        else kernel._add_matter(m0.charged_dust, mrhs0["charged_dust"], _DT)
    )
    m1 = V55MatterState(*m1_base, charged_dust=m1_charged)

    validate_species_admissibility(
        metric_slice_from_q(grid, g_explicit1), m1.radiation,
        Species.RADIATION, stage="explicit-predictor",
    )
    g_explicit1.alpha = kernel._solve_lapse_for_state(
        grid, g_explicit1, s1, m1, em0
    )[0]
    validate_species_admissibility(
        metric_slice_from_q(grid, g_explicit1), m1.radiation,
        Species.RADIATION, stage="explicit-predictor",
    )
    predictor = ProductionState(
        grid=grid,
        geometry=g_explicit1,
        scalars=s1,
        matter=m1,
        electromagnetic=em0,
        t=accepted.t,
        tau=accepted.tau,
        e_folds=accepted.e_folds,
    )
    blocks_pred, _, _, record_pred = _observe_stage(
        kernel, native_grid, predictor, evaluate_coupled_rhs=False
    )
    record_pred["stage"] = "explicit_predictor_cmc"

    # Primary predictor uses the existing CMC L2/L3 split and source channels.
    # Apply the established outer boundary only after Lambda is current, then
    # refresh CMC lapse and pin beta/B before observing the next RHS slice.
    g1 = g_explicit1.copy()
    g1.Aa = g0.Aa + _DT * (
        0.5 * blocks0["primary_l2"]["Aa"]
        + 0.5 * blocks_pred["primary_l2"]["Aa"]
        + blocks0["primary_l3"]["Aa"]
    )
    g1.K = g0.K + _DT * (
        0.5 * blocks0["primary_l2"]["K"]
        + 0.5 * blocks_pred["primary_l2"]["K"]
        + blocks0["primary_l3"]["K"]
    )
    g1.Lambda = g0.Lambda + _DT * (
        0.5 * blocks0["lambda_l2"]
        + 0.5 * blocks_pred["lambda_l2"]
        + blocks0["lambda_l3"]
    )
    validate_species_admissibility(
        metric_slice_from_q(grid, g1), m1.radiation,
        Species.RADIATION, stage="primary-predictor",
    )
    kernel._apply_outer_boundary_for_state(grid, g1, s1, m1, em0)
    g1.alpha = kernel._solve_lapse_for_state(
        grid, g1, s1, m1, em0
    )[0]
    g1.beta.fill(0.0)
    g1.B.fill(0.0)
    g1.assert_finite_positive()
    validate_species_admissibility(
        metric_slice_from_q(grid, g1), m1.radiation,
        Species.RADIATION, stage="primary-predictor",
    )

    primary_predictor = ProductionState(
        grid=grid,
        geometry=g1,
        scalars=s1,
        matter=m1,
        electromagnetic=em0,
        t=accepted.t,
        tau=accepted.tau,
        e_folds=accepted.e_folds,
    )
    blocks1, _, _, record1 = _observe_stage(
        kernel, native_grid, primary_predictor, evaluate_coupled_rhs=True
    )
    record1["stage"] = "primary_predictor_cmc"

    accepted_hash_after = _slice_sha256(
        accepted.geometry, accepted.scalars, accepted.matter
    )
    assert accepted_hash_after == accepted_hash_before, (
        "accepted production input changed while observing copied CMC stages"
    )

    source = inspect.getsource(V55ProductionKernel.step).encode("utf-8")
    report = {
        "schema": "native-cmc-explicit-rhs-stage-shadow/v1",
        "production_step_called": False,
        "stage_build": (
            "copied-state reconstruction using current production-owned lapse/RHS/"
            "boundary helpers and the merged source-level stage-order contract"
        ),
        "production_step_source_sha256": hashlib.sha256(source).hexdigest(),
        "resolution": int(grid.n),
        "r_max": float(grid.r_max),
        "dt": _DT,
        "native_and_production_grid_centres_match": True,
        "accepted_state_sha256_before": accepted_hash_before,
        "accepted_state_sha256_after": accepted_hash_after,
        "stage_observations": [record0, record_pred, record1],
        "reuse": [
            "PR #340 initialized frozen CMC RHS/source decomposition",
            "PR #351 CMC stage-order/boundary/reset source contract",
        ],
        "new_numeric_thresholds": False,
        "production_runtime_selector_added": False,
        "scientific_admission": False,
        "interpretation_limit": (
            "This compares a/b/X explicit RHS on three source-reconstructed CMC "
            "stage slices; it is not a full production-step capture, full-step "
            "equivalence result, or scientific-admission witness."
        ),
    }
    with capsys.disabled():
        print(
            "NATIVE_CMC_EXPLICIT_RHS_STAGE_SHADOW "
            + json.dumps(report, sort_keys=True, separators=(",", ":"), allow_nan=False),
            flush=True,
        )
