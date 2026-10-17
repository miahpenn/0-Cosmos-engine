"""CMC-specific native explicit metric interface contract."""
from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from engine.native_geometry import (
    NativeSphericalCellGrid,
    native_cmc_explicit_metric_rhs,
    native_explicit_rhs,
)


def _cmc_metric_state(grid):
    r = np.asarray(grid.centers, dtype=float)
    envelope = np.exp(-(r / 8.0) ** 2)
    # A smooth, positive CMC slice with beta pinned exactly to zero. This
    # state intentionally has no B-driver field.
    return SimpleNamespace(
        a=1.0 + 0.012 * envelope,
        b=1.0 + 0.009 * envelope,
        X=1.0 + 0.006 * envelope,
        alpha=0.92 + 0.08 * (1.0 - envelope),
        beta=np.zeros_like(r),
        Aa=2.0e-4 * envelope,
        K=1.0e-3 * envelope,
    )


def test_cmc_metric_api_returns_only_shared_exact_native_a_b_x_rhs():
    grid = NativeSphericalCellGrid(80, 40.0)
    state = _cmc_metric_state(grid)
    before = {
        name: np.asarray(getattr(state, name)).copy()
        for name in ("a", "b", "X", "alpha", "beta", "Aa", "K")
    }

    cmc_rhs = native_cmc_explicit_metric_rhs(grid, state)
    assert set(cmc_rhs) == {"a", "b", "X"}
    for name, values in cmc_rhs.items():
        assert values.shape == grid.centers.shape
        assert np.all(np.isfinite(values))

    # The CMC-specific interface and the general native explicit/gauge API
    # must share the exact metric formula implementation.
    full_state = SimpleNamespace(**vars(state), B=np.zeros_like(grid.centers))
    full_rhs = native_explicit_rhs(grid, full_state)
    for name in ("a", "b", "X"):
        np.testing.assert_array_equal(cmc_rhs[name], full_rhs[name])

    # The CMC interface neither requires B nor mutates any input profile.
    assert not hasattr(state, "B")
    for name, values in before.items():
        np.testing.assert_array_equal(getattr(state, name), values)


def test_cmc_metric_api_rejects_a_nonzero_shift_instead_of_using_gauge_rhs():
    grid = NativeSphericalCellGrid(80, 40.0)
    state = _cmc_metric_state(grid)
    shifted = SimpleNamespace(**vars(state))
    shifted.beta = 1.0e-3 * grid.centers * np.exp(-(grid.centers / 8.0) ** 2)

    with pytest.raises(ValueError, match="requires beta=0"):
        native_cmc_explicit_metric_rhs(grid, shifted)
