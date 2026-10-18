"""Integration contract for the bounded native CMC explicit-metric seam."""
from __future__ import annotations

import hashlib
import numpy as np

from engine import v55_pirk_adapter as adapter
from engine.native_geometry import NativeSphericalCellGrid, native_cmc_explicit_metric_rhs
from engine.production_kernel import LAMBDA_M, V55ProductionKernel


def _digest(*objects):
    h = hashlib.sha256()
    for obj in objects:
        for name in getattr(obj, "__dataclass_fields__", {}):
            value = getattr(obj, name)
            if isinstance(value, np.ndarray):
                arr = np.ascontiguousarray(value)
                h.update(name.encode() + b"\0" + arr.dtype.str.encode() + b"\0")
                h.update(repr(arr.shape).encode() + b"\0" + arr.tobytes())
    return h.hexdigest()


def _tree_equal(actual, expected, path="block"):
    assert actual.keys() == expected.keys(), f"{path}: keys differ"
    for key in actual:
        if isinstance(actual[key], dict):
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
