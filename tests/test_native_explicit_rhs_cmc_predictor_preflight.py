"""Stage-local parity check for the native explicit RHS on a CMC predictor.

The test constructs only the copied stage-0 and explicit-predictor inputs using
the existing production stage partition. It never calls V55ProductionKernel.step.
The first initialized CMC slice is already covered by #340; this follow-on asks
only whether native and pinned-adapter a/b/X RHS agree on the later, CMC-resolved
explicit-predictor slice.
"""
from __future__ import annotations

import hashlib
import inspect
import json

import numpy as np

from engine.matter_system import Species, validate_species_admissibility
from engine.native_geometry import NativeSphericalCellGrid, native_explicit_rhs
from engine.production_kernel import ProductionState, V55ProductionKernel
from engine.scalar_system import ScalarFields
from engine.v55_matter import V55MatterState, metric_slice_from_q


_GEOMETRY_FIELDS = (
    "a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B",
)
_SCALAR_FIELDS = ("S", "PS", "D", "PD", "phi", "Pi")
_MATTER_FIELDS = ("rest", "energy_t", "momentum_r")
_MATTER_SPECIES = ("dark_matter", "baryons", "radiation", "charged_dust")
_EXPLICIT_FIELDS = ("a", "b", "X")


def _slice_sha256(geometry, scalars, matter):
    digest = hashlib.sha256()
    for label, obj, names in (
        ("geometry", geometry, _GEOMETRY_FIELDS),
        ("scalars", scalars, _SCALAR_FIELDS),
    ):
        for name in names:
            values = np.ascontiguousarray(np.asarray(getattr(obj, name)))
            digest.update(label.encode("ascii") + b"\0" + name.encode("ascii") + b"\0")
            digest.update(values.dtype.str.encode("ascii") + b"\0")
            digest.update(repr(values.shape).encode("ascii") + b"\0")
            digest.update(values.tobytes(order="C"))
    for species in _MATTER_SPECIES:
        obj = getattr(matter, species, None)
        if obj is None:
            digest.update(b"matter\0" + species.encode("ascii") + b"\0none\0")
            continue
        for name in _MATTER_FIELDS:
            values = np.ascontiguousarray(np.asarray(getattr(obj, name)))
            digest.update(b"matter\0" + species.encode("ascii") + b"\0")
            digest.update(name.encode("ascii") + b"\0")
            digest.update(values.dtype.str.encode("ascii") + b"\0")
            digest.update(repr(values.shape).encode("ascii") + b"\0")
            digest.update(values.tobytes(order="C"))
    return digest.hexdigest()


def _finite_profile(values, name, expected_shape):
    values = np.asarray(values, dtype=float)
    assert values.shape == expected_shape, f"{name}: expected {expected_shape}, got {values.shape}"
    assert np.all(np.isfinite(values)), f"{name}: non-finite values"
    return values


def test_native_explicit_rhs_matches_adapter_on_cmc_resolved_predictor(
    monkeypatch, capsys
):
    """Compare same-input kinematic RHS on a copied later CMC stage only."""
    kernel = V55ProductionKernel()
    accepted = kernel.initialize(
        resolution=48, r_max=12.0, include_radiation=True
    )
    assert accepted.electromagnetic is None, (
        "this bounded parity lane covers the established default no-EM CMC path"
    )

    accepted_hash_before = _slice_sha256(
        accepted.geometry, accepted.scalars, accepted.matter
    )

    def forbidden_step(*_args, **_kwargs):
        raise AssertionError(
            "explicit-RHS CMC predictor preflight must never call production .step"
        )

    monkeypatch.setattr(kernel, "step", forbidden_step)

    grid = accepted.grid
    native_grid = NativeSphericalCellGrid(grid.n, grid.r_max)
    np.testing.assert_allclose(
        native_grid.centers, grid.centers, rtol=0.0, atol=1.0e-14
    )

    # Stage 0: use copies and the same production CMC lapse/gauge ownership.
    g0 = accepted.geometry.copy()
    s0 = accepted.scalars.copy()
    m0 = accepted.matter.copy()
    g0.alpha = kernel._solve_lapse_for_state(grid, g0, s0, m0, None)[0]
    g0.beta.fill(0.0)
    g0.B.fill(0.0)
    g0.assert_finite_positive()
    validate_species_admissibility(
        metric_slice_from_q(grid, g0), m0.radiation, Species.RADIATION,
        stage="native-explicit-cmc-stage0-preflight",
    )

    dt = 1.0e-5
    stage0_hash_before = _slice_sha256(g0, s0, m0)
    stage0 = ProductionState(
        grid=grid,
        geometry=g0,
        scalars=s0,
        matter=m0,
        electromagnetic=None,
        t=accepted.t,
        tau=accepted.tau,
        e_folds=accepted.e_folds,
    )
    scalar_rhs0, matter_rhs0, _ = kernel._rhs(stage0)
    terms0 = kernel._geometry_stage_terms(grid, g0, s0, m0, None)
    stage0_hash_after = _slice_sha256(g0, s0, m0)
    assert stage0_hash_after == stage0_hash_before, (
        "stage-0 RHS evaluation mutated its copied inputs"
    )

    # Construct only the explicit predictor partition, matching production's
    # ownership: advance a/b/X, advance scalar/matter inputs, and re-solve alpha.
    predictor = g0.copy()
    for name in _EXPLICIT_FIELDS:
        setattr(
            predictor,
            name,
            getattr(g0, name) + dt * np.asarray(terms0["explicit"][name], dtype=float),
        )
    predictor.beta.fill(0.0)
    predictor.B.fill(0.0)

    scalars1 = ScalarFields(*(
        getattr(s0, name) + dt * getattr(scalar_rhs0, name)
        for name in _SCALAR_FIELDS
    ))
    matter1_base = tuple(
        kernel._add_matter(getattr(m0, name), matter_rhs0[name], dt)
        for name in ("dark_matter", "baryons", "radiation")
    )
    matter1 = V55MatterState(*matter1_base, charged_dust=None)

    validate_species_admissibility(
        metric_slice_from_q(grid, predictor), matter1.radiation,
        Species.RADIATION, stage="native-explicit-cmc-predictor-preflight",
    )
    predictor.alpha = kernel._solve_lapse_for_state(
        grid, predictor, scalars1, matter1, None
    )[0]
    predictor.beta.fill(0.0)
    predictor.B.fill(0.0)
    predictor.assert_finite_positive()
    validate_species_admissibility(
        metric_slice_from_q(grid, predictor), matter1.radiation,
        Species.RADIATION, stage="native-explicit-cmc-predictor-preflight",
    )

    predictor_hash_before = _slice_sha256(predictor, scalars1, matter1)
    native = native_explicit_rhs(native_grid, predictor)
    reference_terms = kernel._geometry_stage_terms(
        grid, predictor, scalars1, matter1, None
    )
    reference = reference_terms["explicit"]
    predictor_hash_after = _slice_sha256(predictor, scalars1, matter1)
    assert predictor_hash_after == predictor_hash_before, (
        "same-input explicit RHS evaluations mutated the CMC predictor"
    )

    residuals = {}
    for name in _EXPLICIT_FIELDS:
        native_values = _finite_profile(native[name], f"native {name}", grid.centers.shape)
        reference_values = _finite_profile(
            reference[name], f"adapter {name}", grid.centers.shape
        )
        residual = native_values - reference_values
        assert np.all(np.isfinite(residual)), f"{name}: non-finite residual"
        # #340 found exact same-slice explicit parity. The predictor-stage gate
        # intentionally asks for exact equality too; it adds no fitted tolerance.
        np.testing.assert_array_equal(
            native_values,
            reference_values,
            err_msg=f"native/adapter explicit RHS differs for {name} on CMC predictor",
        )
        residuals[name] = {
            "array_equal": bool(np.array_equal(native_values, reference_values)),
            "max_abs_native_minus_adapter": float(np.max(np.abs(residual))),
            "rms_native_minus_adapter": float(np.sqrt(np.mean(residual**2))),
        }

    accepted_hash_after = _slice_sha256(
        accepted.geometry, accepted.scalars, accepted.matter
    )
    assert accepted_hash_after == accepted_hash_before, (
        "preflight mutated caller-owned accepted state"
    )

    record = {
        "schema": "native-explicit-rhs-cmc-predictor-preflight/v1",
        "resolution": int(grid.n),
        "r_max": float(grid.r_max),
        "dt_for_copied_predictor": dt,
        "comparison_stage": "explicit predictor after CMC lapse refresh",
        "compared_rhs_fields": list(_EXPLICIT_FIELDS),
        "excluded_rhs_fields": {
            "alpha": "CMC-solved gauge variable; not substituted from moving-puncture RHS",
            "beta": "pinned to zero in production CMC path",
            "B": "pinned to zero in production CMC path",
            "primary_L2_L3_and_Lambda": "outside this one-block explicit-RHS lane",
        },
        "native_and_production_grid_centres_match": True,
        "native_explicit_rhs_source_sha256": hashlib.sha256(
            inspect.getsource(native_explicit_rhs).encode("utf-8")
        ).hexdigest(),
        "production_geometry_stage_source_sha256": hashlib.sha256(
            inspect.getsource(kernel._geometry_stage_terms).encode("utf-8")
        ).hexdigest(),
        "production_step_called": False,
        "accepted_state_sha256_before": accepted_hash_before,
        "accepted_state_sha256_after": accepted_hash_after,
        "accepted_state_unchanged": True,
        "stage0_sha256_before_rhs": stage0_hash_before,
        "stage0_sha256_after_rhs": stage0_hash_after,
        "stage0_inputs_unchanged": True,
        "predictor_sha256_before_rhs": predictor_hash_before,
        "predictor_sha256_after_rhs": predictor_hash_after,
        "predictor_inputs_unchanged": True,
        "native_minus_adapter_explicit_rhs": residuals,
        "new_residual_thresholds": False,
        "runtime_selector_added": False,
        "scientific_admission": False,
        "interpretation_limit": (
            "single-resolution, same-input explicit RHS implementation parity only; "
            "not a full-step/runtime substitution or scientific-admission result"
        ),
    }
    with capsys.disabled():
        print(
            "NATIVE_EXPLICIT_RHS_CMC_PREDICTOR_PREFLIGHT "
            + json.dumps(record, sort_keys=True, separators=(",", ":"), allow_nan=False),
            flush=True,
        )
