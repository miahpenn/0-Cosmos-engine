"""Stage-order tests for the undamped, source-backed native PIRK2 candidate."""

import hashlib
import json
from types import SimpleNamespace

import numpy as np
import pytest

from engine.native_geometry import (
    NativeSphericalCellGrid,
    metric_connection_lambda,
    native_bssn_constraints,
    native_bssn_geometry_terms,
    native_pirk2_step,
)
from engine.v55_pirk_adapter import vendor_modules
import engine.native_geometry.evolution as evolution_module


_FIELDS = ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B")


def _state(grid, *, curved=False):
    r = grid.centers
    zeros = np.zeros(grid.n, dtype=float)
    if curved:
        f = 4.0e-4 * r**2 + 6.0e-6 * r**4
        a = np.exp(2.0 * f)
        b = np.exp(-f)
        X = 1.0 + 2.0e-4 * r**2 / (1.0 + r**2)
        alpha = 1.0 + 0.01 * r**2 / (1.0 + r**2)
        beta = 0.009 * r / (1.0 + 0.18 * r**2)
        Aa = 1.8e-4 * r**2 / (1.0 + 0.12 * r**2)
        K = 0.018 + 3.5e-4 * r**2 / (1.0 + r**2)
    else:
        a = b = X = alpha = np.ones(grid.n, dtype=float)
        beta = Aa = K = zeros.copy()
    state = SimpleNamespace(
        a=np.asarray(a, dtype=float).copy(),
        b=np.asarray(b, dtype=float).copy(),
        X=np.asarray(X, dtype=float).copy(),
        alpha=np.asarray(alpha, dtype=float).copy(),
        beta=np.asarray(beta, dtype=float).copy(),
        Aa=np.asarray(Aa, dtype=float).copy(),
        K=np.asarray(K, dtype=float).copy(),
        Lambda=zeros.copy(),
        B=zeros.copy(),
    )
    if curved:
        state.Lambda = metric_connection_lambda(grid, state)
        state.B = 0.75 * state.Lambda
    # Preserve the exact flat control instead of importing roundoff from a
    # reconstructed connection diagnostic into an otherwise-zero state.
    return state


def _vendor_state(vacuum, state):
    return vacuum.VacuumState(**{
        name: np.asarray(getattr(state, name), dtype=float).copy()
        for name in _FIELDS
    })


def test_flat_zero_state_is_fixed_by_native_pirk2_step():
    grid = NativeSphericalCellGrid(64, 12.0)
    state = _state(grid)
    result = native_pirk2_step(grid, state, 1.0e-3)

    for name in _FIELDS:
        np.testing.assert_allclose(
            getattr(result, name), getattr(state, name), rtol=0.0, atol=1.0e-16,
            err_msg=f"{name} should remain fixed to floating-point roundoff",
        )
        np.testing.assert_array_equal(getattr(state, name), getattr(_state(grid), name))


def test_native_pirk2_step_is_read_only_and_rejects_invalid_step_contracts():
    grid = NativeSphericalCellGrid(64, 12.0)
    state = _state(grid, curved=True)
    before = {name: getattr(state, name).copy() for name in _FIELDS}

    result = native_pirk2_step(grid, state, 1.0e-3, lambda_m=1.7)
    # Whole-step output must preserve the independently sourced undamped
    # gauge law, not merely the standalone native B RHS.
    np.testing.assert_allclose(
        result.B,
        state.B + 0.75 * (result.Lambda - state.Lambda),
        rtol=2.0e-14,
        atol=2.0e-14,
    )
    for name in _FIELDS:
        np.testing.assert_array_equal(getattr(state, name), before[name])
        assert np.all(np.isfinite(getattr(result, name))), name

    with pytest.raises(ValueError, match="dt must be positive"):
        native_pirk2_step(grid, state, 0.0)
    with pytest.raises(ValueError, match="dt must be finite"):
        native_pirk2_step(grid, state, np.nan)
    with pytest.raises(ValueError, match="dt must be positive"):
        native_pirk2_step(grid, state, -1.0e-3)
    with pytest.raises(ValueError, match="dt must be finite"):
        native_pirk2_step(grid, state, np.inf)
    with pytest.raises(ValueError, match="dt must be a finite real scalar"):
        native_pirk2_step(grid, state, True)

    bad_shape = SimpleNamespace(**vars(state))
    bad_shape.Lambda = np.zeros(grid.n - 1)
    with pytest.raises(ValueError, match="Lambda must have shape"):
        native_pirk2_step(grid, bad_shape, 1.0e-3)

    bad_finite = SimpleNamespace(**vars(state))
    bad_finite.B = np.full(grid.n, np.inf)
    with pytest.raises(ValueError, match="B must be finite"):
        native_pirk2_step(grid, bad_finite, 1.0e-3)

    bad_lambda_m = SimpleNamespace(**vars(state))
    with pytest.raises(ValueError, match="lambda_m must be finite"):
        native_pirk2_step(grid, bad_lambda_m, 1.0e-3, lambda_m=np.nan)


def _state_fingerprint(state):
    """Hash all named input profiles with dtype, shape, and raw bytes."""
    digest = hashlib.sha256()
    for name in _FIELDS:
        values = np.ascontiguousarray(np.asarray(getattr(state, name)))
        digest.update(name.encode("ascii") + b"\0")
        digest.update(values.dtype.str.encode("ascii") + b"\0")
        digest.update(repr(values.shape).encode("ascii") + b"\0")
        digest.update(values.tobytes(order="C"))
    return digest.hexdigest()


def _region_residual_summary(lhs, rhs, masks, names):
    """Summarize signed profile differences region-by-region, after differencing."""
    summary = {}
    for region, mask in masks.items():
        count = int(np.count_nonzero(mask))
        assert count >= 4, f"{region}: insufficient cells for a region summary ({count})"
        region_fields = {}
        for name in names:
            signed = np.asarray(getattr(lhs, name), dtype=float) - np.asarray(
                getattr(rhs, name), dtype=float
            )
            assert np.all(np.isfinite(signed)), f"{region}/{name}: non-finite signed residual"
            selected = signed[mask]
            region_fields[name] = {
                "rms": float(np.sqrt(np.mean(selected**2))),
                "max_abs": float(np.max(np.abs(selected))),
            }
        summary[region] = {"cells": count, "fields": region_fields}
    return summary


def _constraint_residual_summary(lhs, rhs, masks):
    """Apply one source-free native constraint convention to both output states."""
    names = ("H", "M", "connection", "determinant")
    result = {}
    for region, mask in masks.items():
        fields = {}
        for name in names:
            # Difference the full signed profile before restricting/reducing.
            signed = np.asarray(lhs[name], dtype=float) - np.asarray(rhs[name], dtype=float)
            assert np.all(np.isfinite(signed)), f"{region}/{name}: non-finite constraint residual"
            selected = signed[mask]
            fields[name] = {
                "rms": float(np.sqrt(np.mean(selected**2))),
                "max_abs": float(np.max(np.abs(selected))),
            }
        result[region] = {"fields": fields}
    return result


def test_native_pirk2_step_matches_pinned_stage_order_on_fixed_areal_annulus(capsys):
    operators, vacuum, moving = vendor_modules()
    rms_by_field = {name: [] for name in _FIELDS}
    shadow_records = []

    for n in (48, 96, 192):
        grid = NativeSphericalCellGrid(n, 12.0)
        vendor_grid = operators.SphericalCellGrid(n, 12.0)
        state = _state(grid, curved=True)
        vendor_state = _vendor_state(vacuum, state)
        input_hashes_before = {
            "native": _state_fingerprint(state),
            "legacy_copy": _state_fingerprint(vendor_state),
        }
        dt = 1.0e-3
        # Only the source-backed undamped Gamma-driver convention is in scope.
        native = native_pirk2_step(grid, state, dt, lambda_m=2.0)
        reference = moving.damped_moving_puncture_pirk2_step(
            vendor_grid, vendor_state, dt, eta=0.0, ko_epsilon=0.0, lambda_m=2.0
        )
        input_hashes_after = {
            "native": _state_fingerprint(state),
            "legacy_copy": _state_fingerprint(vendor_state),
        }
        assert input_hashes_after == input_hashes_before, (
            f"native/legacy paired comparison mutated an input at N={n}: "
            f"before={input_hashes_before}, after={input_hashes_after}"
        )

        areal = native_bssn_geometry_terms(grid, state)["areal_radius"]
        masks = {
            "centre_R_areal_lt_2": areal < 2.0,
            "fixed_annulus_2_le_R_areal_le_6": (areal >= 2.0) & (areal <= 6.0),
            "outer_R_areal_gt_6": areal > 6.0,
        }
        assert int(np.count_nonzero(masks["fixed_annulus_2_le_R_areal_le_6"])) >= 8

        field_regions = _region_residual_summary(native, reference, masks, _FIELDS)
        zero_sources = np.zeros(grid.n)
        constraints_native = native_bssn_constraints(
            grid, native, rho=zero_sources, radial_momentum_source=zero_sources
        )
        constraints_reference = native_bssn_constraints(
            grid, reference, rho=zero_sources, radial_momentum_source=zero_sources
        )
        constraint_regions = _constraint_residual_summary(
            constraints_native, constraints_reference, masks
        )
        for name in _FIELDS:
            rms_by_field[name].append(
                field_regions["fixed_annulus_2_le_R_areal_le_6"]["fields"][name]["rms"]
            )

        physical_regions = {}
        for region, region_summary in field_regions.items():
            mask = masks[region]
            physical_regions[region] = {
                "cells": region_summary["cells"],
                "r_areal_min": float(np.min(areal[mask])),
                "r_areal_max": float(np.max(areal[mask])),
                "field_signed_residuals": region_summary["fields"],
            }
        shadow_records.append({
            "resolution": n,
            "dt": dt,
            "physical_regions": physical_regions,
            "source_free_native_constraint_residuals": constraint_regions,
            "input_sha256_before": input_hashes_before,
            "input_sha256_after": input_hashes_after,
        })

    # Emit a compact machine-readable witness to the Actions job log even when
    # pytest passes. This is diagnostic output; only the existing annulus
    # comparison below is an acceptance threshold.
    with capsys.disabled():
        print(
            "NATIVE_PIRK2_SHADOW_WITNESS "
            + json.dumps(
                {"schema": "native-pirk2-shadow-witness/v1", "records": shadow_records},
                sort_keys=True,
                separators=(",", ":"),
                allow_nan=False,
            ),
            flush=True,
        )

    # Preserve the established annulus bound exactly; centre and outer metrics
    # are descriptive only until their error envelopes are independently justified.
    for name, errors in rms_by_field.items():
        assert len(errors) == 3 and np.all(np.isfinite(errors)), name
        assert max(errors) < 1.0e-4, f"{name} native/vendor PIRK2 RMS: {errors}"

def test_native_pirk2_stage_schedule_uses_declared_time_levels(monkeypatch):
    """Independent scalar witness that discriminates PIRK2 evaluation stages."""
    grid = NativeSphericalCellGrid(8, 8.0)
    state = SimpleNamespace(
        a=np.ones(grid.n),
        b=np.ones(grid.n),
        X=np.ones(grid.n),
        alpha=np.ones(grid.n),
        beta=np.zeros(grid.n),
        Aa=np.full(grid.n, 0.1),
        K=np.full(grid.n, 0.2),
        Lambda=np.full(grid.n, 0.3),
        B=np.full(grid.n, 0.225),
    )
    n = grid.n

    def explicit(g, s):
        return {
            "a": np.full(n, 0.01 * s.a[0]),
            "b": np.full(n, 0.02 * s.b[0]),
            "X": np.full(n, 0.03 * s.X[0]),
            "alpha": np.full(n, 0.04 * s.alpha[0]),
            "beta": s.B.copy(),
        }

    def primary_l2(g, s):
        return {
            "Aa": np.full(n, 0.1 * s.a[0]),
            "K": np.full(n, 0.2 * s.alpha[0]),
        }

    def primary_l3(g, s):
        return {
            "Aa": np.full(n, 0.3 * s.beta[0] + 0.04),
            "K": np.full(n, 0.5 * s.Aa[0] + 0.01),
        }

    def lambda_l2(g, s, *, lambda_m=2.0):
        return np.full(n, 0.6 * s.Aa[0] + 0.1 * s.a[0])

    def lambda_l3(g, s):
        return np.full(n, 0.25 * s.K[0] + 0.02)

    monkeypatch.setattr(evolution_module, "native_explicit_rhs", explicit)
    monkeypatch.setattr(evolution_module, "native_primary_l2_rhs", primary_l2)
    monkeypatch.setattr(evolution_module, "native_primary_l3_rhs", primary_l3)
    monkeypatch.setattr(evolution_module, "native_lambda_l2_rhs", lambda_l2)
    monkeypatch.setattr(evolution_module, "native_lambda_l3_rhs", lambda_l3)

    result = native_pirk2_step(grid, state, 0.1)

    # Independently evaluated scalar result of the declared u -> v -> Lambda
    # predictor/corrector. These values discriminate the L2 time levels and
    # the use of the same Lambda stage rate for the undamped B driver.
    expected = {
        "a": 1.0010005,
        "b": 1.002002,
        "X": 1.0030045,
        "alpha": 1.004008,
        "beta": 0.023378443125,
        "Aa": 0.1143425025,
        "K": 0.226390205,
        "Lambda": 0.323760777575,
        "B": 0.24282058318125,
    }
    for name, value in expected.items():
        np.testing.assert_allclose(
            getattr(result, name), np.full(n, value),
            rtol=0.0, atol=5.0e-14, err_msg=name,
        )




def _snapshot_native_stage(state):
    """Freeze one real stage state observed at an RHS call, without changing the stepper."""
    return SimpleNamespace(**{
        name: np.asarray(getattr(state, name)).copy()
        for name in _FIELDS
    })


def test_native_pirk2_actual_stage_states_shadow_against_pinned_vendor(
    monkeypatch, capsys
):
    """Capture real predictor/corrector stage inputs through the existing RHS calls."""
    operators, vacuum, moving = vendor_modules()
    native_names = (
        "native_explicit_rhs",
        "native_primary_l2_rhs",
        "native_primary_l3_rhs",
        "native_lambda_l2_rhs",
        "native_lambda_l3_rhs",
    )
    vendor_names = (
        "moving_puncture_explicit_rhs",
        "primary_l2_rhs",
        "primary_l3_rhs",
        "lambda_l2_rhs",
        "lambda_l3_rhs",
    )
    native_originals = {
        name: getattr(evolution_module, name) for name in native_names
    }
    vendor_originals = {name: getattr(moving, name) for name in vendor_names}
    stage_records = []

    def instrument(original, bucket, key):
        def wrapped(grid, state, *args, **kwargs):
            bucket[key].append(_snapshot_native_stage(state))
            return original(grid, state, *args, **kwargs)
        return wrapped

    for n in (48, 96, 192):
        grid = NativeSphericalCellGrid(n, 12.0)
        vendor_grid = operators.SphericalCellGrid(n, 12.0)
        state = _state(grid, curved=True)
        vendor_state = _vendor_state(vacuum, state)
        hashes_before = {
            "native": _state_fingerprint(state),
            "legacy_copy": _state_fingerprint(vendor_state),
        }
        native_calls = {
            "explicit": [], "primary_l2": [], "primary_l3": [],
            "lambda_l2": [], "lambda_l3": [],
        }
        vendor_calls = {
            "explicit": [], "primary_l2": [], "primary_l3": [],
            "lambda_l2": [], "lambda_l3": [],
        }

        native_aliases = {
            "native_explicit_rhs": "explicit",
            "native_primary_l2_rhs": "primary_l2",
            "native_primary_l3_rhs": "primary_l3",
            "native_lambda_l2_rhs": "lambda_l2",
            "native_lambda_l3_rhs": "lambda_l3",
        }
        vendor_aliases = {
            "moving_puncture_explicit_rhs": "explicit",
            "primary_l2_rhs": "primary_l2",
            "primary_l3_rhs": "primary_l3",
            "lambda_l2_rhs": "lambda_l2",
            "lambda_l3_rhs": "lambda_l3",
        }
        with monkeypatch.context() as patch:
            for name, key in native_aliases.items():
                patch.setattr(
                    evolution_module, name,
                    instrument(native_originals[name], native_calls, key),
                )
            for name, key in vendor_aliases.items():
                patch.setattr(
                    moving, name,
                    instrument(vendor_originals[name], vendor_calls, key),
                )
            native_result = native_pirk2_step(grid, state, 1.0e-3, lambda_m=2.0)
            vendor_result = moving.damped_moving_puncture_pirk2_step(
                vendor_grid, vendor_state, 1.0e-3,
                eta=0.0, ko_epsilon=0.0, lambda_m=2.0,
            )

        hashes_after = {
            "native": _state_fingerprint(state),
            "legacy_copy": _state_fingerprint(vendor_state),
        }
        assert hashes_after == hashes_before, (
            f"stage-instrumented step mutated an input at N={n}: "
            f"before={hashes_before}, after={hashes_after}"
        )

        # Both steppers invoke these real RHS stages in the same semantic order.
        # Snapshots are copied inside the wrappers, so later Lambda/B assignments
        # cannot retroactively alter earlier pre-update stage evidence.
        stage_pairs = {
            "U1_predictor": (native_calls["primary_l2"][1], vendor_calls["primary_l2"][1]),
            "V1_before_Lambda_update": (
                native_calls["lambda_l2"][1], vendor_calls["lambda_l2"][1]
            ),
            "stage1_after_Lambda_B_update": (
                native_calls["explicit"][1], vendor_calls["explicit"][1]
            ),
            "U_corrector_before_V_update": (
                native_calls["primary_l2"][2], vendor_calls["primary_l2"][2]
            ),
            "V_corrector_before_final_Lambda_B_update": (
                native_calls["lambda_l2"][2], vendor_calls["lambda_l2"][2]
            ),
        }
        expected_counts = {
            "explicit": 2, "primary_l2": 3, "primary_l3": 2,
            "lambda_l2": 3, "lambda_l3": 2,
        }
        for key, expected_count in expected_counts.items():
            assert len(native_calls[key]) == expected_count, (
                f"native {key} call count at N={n}: {len(native_calls[key])}"
            )
            assert len(vendor_calls[key]) == expected_count, (
                f"pinned vendor {key} call count at N={n}: {len(vendor_calls[key])}"
            )

        areal = native_bssn_geometry_terms(grid, state)["areal_radius"]
        masks = {
            "centre_R_areal_lt_2": areal < 2.0,
            "fixed_annulus_2_le_R_areal_le_6": (areal >= 2.0) & (areal <= 6.0),
            "outer_R_areal_gt_6": areal > 6.0,
        }
        per_stage = {}
        for stage_name, (native_stage, vendor_stage) in stage_pairs.items():
            per_stage[stage_name] = _region_residual_summary(
                native_stage, vendor_stage, masks, _FIELDS
            )
        # The return values are included as a stage mapping cross-check. Their
        # existing acceptance gate lives in the prior shadow-witness test.
        per_stage["final_output_crosscheck"] = _region_residual_summary(
            native_result, vendor_result, masks, _FIELDS
        )
        stage_records.append({
            "resolution": n,
            "physical_regions": list(masks.keys()),
            "stages": per_stage,
            "input_sha256_before": hashes_before,
            "input_sha256_after": hashes_after,
        })

    with capsys.disabled():
        print(
            "NATIVE_PIRK2_ACTUAL_STAGE_WITNESS "
            + json.dumps(
                {"schema": "native-pirk2-actual-stage-witness/v1", "records": stage_records},
                sort_keys=True,
                separators=(",", ":"),
                allow_nan=False,
            ),
            flush=True,
        )
