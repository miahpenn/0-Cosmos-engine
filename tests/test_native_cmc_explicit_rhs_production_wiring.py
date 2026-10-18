"""Integration contract for the bounded native CMC explicit-metric seam."""
from __future__ import annotations

import hashlib
import numpy as np

from engine import v55_pirk_adapter as adapter
from engine.native_geometry import NativeSphericalCellGrid, native_cmc_explicit_metric_rhs
from engine.production_kernel import LAMBDA_M, V55ProductionKernel, _copy_final_cmc_state


_GEOMETRY_FIELDS = ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B")
_SCALAR_FIELDS = ("S", "PS", "D", "PD", "phi", "Pi")


def _digest(geometry, scalars):
    h = hashlib.sha256()
    for label, obj, names in (
        ("geometry", geometry, _GEOMETRY_FIELDS),
        ("scalars", scalars, _SCALAR_FIELDS),
    ):
        for name in names:
            arr = np.ascontiguousarray(np.asarray(getattr(obj, name)))
            h.update(label.encode() + b"\0" + name.encode() + b"\0")
            h.update(arr.dtype.str.encode() + b"\0" + repr(arr.shape).encode() + b"\0")
            h.update(arr.tobytes())
    return h.hexdigest()


def _tree_equal(actual, expected, path="block"):
    assert actual.keys() == expected.keys(), f"{path}: keys differ"
    for key in actual:
        if isinstance(expected[key], dict):
            _tree_equal(actual[key], expected[key], f"{path}.{key}")
        else:
            np.testing.assert_array_equal(actual[key], expected[key], err_msg=f"{path}.{key}")


def test_production_cmc_terms_use_native_a_b_x_without_calling_vendor_explicit(monkeypatch):
    kernel = V55ProductionKernel()
    state = kernel.initialize(resolution=48, r_max=12.0, include_radiation=True)
    assert state.electromagnetic is None
    before = _digest(state.geometry, state.scalars)

    reference = adapter.geometry_stage_terms(
        state.grid, state.geometry, state.scalars, state.matter, lambda_m=LAMBDA_M
    )
    native_grid = NativeSphericalCellGrid(state.grid.n, state.grid.r_max)
    np.testing.assert_allclose(
        native_grid.centers, state.grid.centers, rtol=0.0, atol=1.0e-14
    )
    native = native_cmc_explicit_metric_rhs(native_grid, state.geometry)
    for name in ("a", "b", "X"):
        np.testing.assert_array_equal(native[name], reference["explicit"][name])

    _, _, moving = adapter.vendor_modules()

    def forbidden(*_args, **_kwargs):
        raise AssertionError("CMC production path called pinned explicit RHS")

    monkeypatch.setattr(moving, "moving_puncture_explicit_rhs", forbidden)
    actual = kernel._geometry_stage_terms(
        state.grid, state.geometry, state.scalars, state.matter, state.electromagnetic
    )

    assert set(actual["explicit"]) == {"a", "b", "X"}
    for name in ("a", "b", "X"):
        np.testing.assert_array_equal(actual["explicit"][name], native[name])
    for name in ("primary_l2", "primary_l3", "lambda_l2", "lambda_l3"):
        if isinstance(reference[name], dict):
            _tree_equal(actual[name], reference[name], name)
        else:
            np.testing.assert_array_equal(actual[name], reference[name])
    assert _digest(state.geometry, state.scalars) == before


def test_final_cmc_state_copy_keeps_gauge_fields_outside_rhs_key_contract():
    kernel = V55ProductionKernel()
    state = kernel.initialize(resolution=24, r_max=8.0, include_radiation=True)
    target = state.geometry.copy()
    source = state.geometry.copy()
    untouched_Aa = target.Aa.copy()

    for offset, name in enumerate(("a", "b", "X", "alpha", "beta"), start=1):
        setattr(source, name, np.asarray(getattr(source, name), dtype=float) + offset * 1.0e-4)

    _copy_final_cmc_state(target, source)
    for name in ("a", "b", "X", "alpha", "beta"):
        np.testing.assert_array_equal(getattr(target, name), getattr(source, name))
    np.testing.assert_array_equal(target.Aa, untouched_Aa)

    source.alpha[:] = -123.0
    assert not np.array_equal(target.alpha, source.alpha)


def test_native_cmc_override_preserves_non_none_electromagnetic_projection(monkeypatch):
    kernel = V55ProductionKernel()
    state = kernel.initialize(resolution=32, r_max=10.0, include_radiation=True)
    projection = {
        "rho": np.full(state.grid.n, 2.0e-5),
        "pr": np.full(state.grid.n, 3.0e-6),
        "pt": np.full(state.grid.n, 1.0e-6),
        "j": np.linspace(-2.0e-7, 2.0e-7, state.grid.n),
    }

    class ProjectedEM:
        def stress_projection(self, grid, geometry):
            assert grid is state.grid
            return {key: value.copy() for key, value in projection.items()}

    native_grid = NativeSphericalCellGrid(state.grid.n, state.grid.r_max)
    native = native_cmc_explicit_metric_rhs(native_grid, state.geometry)
    reference = adapter.geometry_stage_terms(
        state.grid,
        state.geometry,
        state.scalars,
        state.matter,
        lambda_m=LAMBDA_M,
        electromagnetic_projection=projection,
        explicit_rhs_override=native,
    )

    _, _, moving = adapter.vendor_modules()

    def forbidden(*_args, **_kwargs):
        raise AssertionError("CMC production path called pinned explicit RHS")

    monkeypatch.setattr(moving, "moving_puncture_explicit_rhs", forbidden)
    actual = kernel._geometry_stage_terms(
        state.grid, state.geometry, state.scalars, state.matter, ProjectedEM()
    )

    assert set(actual["explicit"]) == {"a", "b", "X"}
    for name in ("a", "b", "X"):
        np.testing.assert_array_equal(actual["explicit"][name], native[name])
    _tree_equal(actual["primary_l3"], reference["primary_l3"], "primary_l3")
    _tree_equal(actual["lambda_l3"], reference["lambda_l3"], "lambda_l3")
