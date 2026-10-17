"""Contract for the first native production RHS seam.

Only the explicit metric RHS is migrated here. CMC alpha is still solved by
the production lapse routine, beta/B are still reset by the CMC stage contract,
and all primary/Lambda/source blocks remain supplied by the pinned adapter.
"""
from __future__ import annotations

import hashlib

import numpy as np

from engine import v55_pirk_adapter as adapter
from engine.native_geometry import NativeSphericalCellGrid, native_explicit_rhs
from engine.production_kernel import LAMBDA_M, V55ProductionKernel


_GEOMETRY_FIELDS = ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B")
_SCALAR_FIELDS = ("S", "PS", "D", "PD", "phi", "Pi")
_MATTER_FIELDS = ("rest", "energy_t", "momentum_r")
_MATTER_SPECIES = ("dark_matter", "baryons", "radiation", "charged_dust")
_NATIVE_EXPLICIT_FIELDS = ("a", "b", "X", "alpha", "beta")
_NON_EXPLICIT_BLOCKS = ("primary_l2", "primary_l3", "lambda_l2", "lambda_l3")


def _state_sha256(state):
    digest = hashlib.sha256()
    for label, obj, names in (
        ("geometry", state.geometry, _GEOMETRY_FIELDS),
        ("scalars", state.scalars, _SCALAR_FIELDS),
    ):
        for name in names:
            values = np.ascontiguousarray(np.asarray(getattr(obj, name)))
            digest.update(label.encode("ascii") + b"\0" + name.encode("ascii") + b"\0")
            digest.update(values.dtype.str.encode("ascii") + b"\0")
            digest.update(repr(values.shape).encode("ascii") + b"\0")
            digest.update(values.tobytes(order="C"))
    for species_name in _MATTER_SPECIES:
        species = getattr(state.matter, species_name, None)
        if species is None:
            digest.update(b"matter\0" + species_name.encode("ascii") + b"\0none\0")
            continue
        for name in _MATTER_FIELDS:
            values = np.ascontiguousarray(np.asarray(getattr(species, name)))
            digest.update(b"matter\0" + species_name.encode("ascii") + b"\0")
            digest.update(name.encode("ascii") + b"\0")
            digest.update(values.dtype.str.encode("ascii") + b"\0")
            digest.update(repr(values.shape).encode("ascii") + b"\0")
            digest.update(values.tobytes(order="C"))
    return digest.hexdigest()


def _assert_tree_array_equal(actual, expected, prefix):
    assert actual.keys() == expected.keys(), (
        f"{prefix}: block keys differ: {actual.keys()} != {expected.keys()}"
    )
    for key in actual:
        path = f"{prefix}.{key}"
        if isinstance(actual[key], dict):
            assert isinstance(expected[key], dict), f"{path}: expected nested block"
            _assert_tree_array_equal(actual[key], expected[key], path)
        else:
            left = np.asarray(actual[key])
            right = np.asarray(expected[key])
            assert left.shape == right.shape, f"{path}: shape mismatch"
            assert np.array_equal(left, right), f"{path}: unchanged adapter block differs"


def test_production_geometry_terms_use_native_explicit_block_only(monkeypatch):
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=48, r_max=12.0, include_radiation=True
    )
    assert state.electromagnetic is None

    # Preserve an independent pinned-adapter oracle before guarding its explicit
    # routine. All later comparisons use the exact same copied input slice.
    reference = adapter.geometry_stage_terms(
        state.grid,
        state.geometry,
        state.scalars,
        state.matter,
        lambda_m=LAMBDA_M,
    )
    reference_explicit = reference["explicit"]
    native_grid = NativeSphericalCellGrid(state.grid.n, state.grid.r_max)
    np.testing.assert_allclose(
        native_grid.centers, state.grid.centers, rtol=0.0, atol=1.0e-14
    )
    native = native_explicit_rhs(native_grid, state.geometry)
    for name in ("a", "b", "X"):
        np.testing.assert_array_equal(
            native[name], reference_explicit[name],
            err_msg=f"native/adapter explicit RHS {name} on the accepted CMC slice",
        )
    accepted_hash_before = _state_sha256(state)

    _, _, moving = adapter.vendor_modules()

    def forbidden_borrowed_explicit(*_args, **_kwargs):
        raise AssertionError(
            "production CMC stage terms must not invoke the borrowed "
            "moving-puncture explicit RHS after native wiring"
        )

    monkeypatch.setattr(
        moving, "moving_puncture_explicit_rhs", forbidden_borrowed_explicit
    )

    actual = kernel._geometry_stage_terms(
        state.grid, state.geometry, state.scalars, state.matter,
        state.electromagnetic,
    )

    # Native module owns the explicit dictionary in this production path.
    # The CMC step itself consumes only a/b/X; alpha remains elliptically
    # resolved and beta/B remain pinned by the existing stage contract.
    assert actual["explicit"].keys() == native.keys()
    for name in _NATIVE_EXPLICIT_FIELDS:
        np.testing.assert_array_equal(
            actual["explicit"][name], native[name], err_msg=f"native explicit {name}"
        )
    # The unchanged split remains the adapter's: L2/L3 and matter-source
    # contributions must equal the separately captured pinned reference.
    for name in _NON_EXPLICIT_BLOCKS:
        _assert_tree_array_equal(actual[name], reference[name], name)

    accepted_hash_after = _state_sha256(state)
    assert accepted_hash_after == accepted_hash_before, (
        "native production RHS wiring mutated caller-owned accepted state"
    )
