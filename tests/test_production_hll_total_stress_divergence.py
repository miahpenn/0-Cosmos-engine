"""Synthetic contracts for the pinned production HLL total-T observer."""

from __future__ import annotations

import importlib.util
import os
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("COSMOS_PINNED_SOURCE_ROOT", str(ROOT))
_SPEC = importlib.util.spec_from_file_location(
    "_issue4_production_hll_total_stress_divergence",
    ROOT / "tools" / "audit_production_hll_total_stress_divergence.py",
)
if _SPEC is None or _SPEC.loader is None:
    raise ImportError("cannot load production HLL divergence observer")
AUDIT_HLL = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(AUDIT_HLL)


def test_native_component_decomposition_matches_signed_flux_source_identity():
    face_flux = np.asarray([
        [0.0, 0.0, 0.0],
        [1.0, 2.0, 3.0],
        [3.0, 6.0, 9.0],
        [6.0, 12.0, 18.0],
    ])
    source_e = np.asarray([10.0, 20.0, 30.0])
    source_s = np.asarray([2.0, 4.0, 6.0])
    q_e = np.asarray([0.5, 1.0, 1.5])
    q_s = np.asarray([1.0, 2.0, 3.0])

    result = AUDIT_HLL.rhs_from_native_components(
        face_flux, source_e, source_s, q_e, q_s, dr=1.0
    )

    np.testing.assert_array_equal(result["rest"], [-1.0, -2.0, -3.0])
    np.testing.assert_array_equal(result["energy_t"], [8.5, 17.0, 25.5])
    np.testing.assert_array_equal(result["momentum_r"], [0.0, 0.0, 0.0])


def test_negative_controls_retain_sign_flip_and_omitted_geometric_source_effects():
    face_flux = np.asarray([
        [0.0, 0.0, 0.0],
        [1.0, 2.0, 3.0],
        [3.0, 6.0, 9.0],
        [6.0, 12.0, 18.0],
    ])
    source_e = np.asarray([10.0, 20.0, 30.0])
    source_s = np.asarray([2.0, 4.0, 6.0])
    q_e = np.asarray([0.5, 1.0, 1.5])
    q_s = np.asarray([1.0, 2.0, 3.0])

    controls = AUDIT_HLL.negative_control_rhs(
        face_flux, source_e, source_s, q_e, q_s, dr=1.0
    )

    assert controls["sign_flipped_face_flux"]["delta_vs_native"]["energy_t"]["max_abs"] > 0.0
    assert controls["sign_flipped_face_flux"]["delta_vs_native"]["momentum_r"]["max_abs"] > 0.0
    assert controls["omitted_geometric_source"]["delta_vs_native"]["energy_t"]["max_abs"] == 30.0
    assert controls["omitted_geometric_source"]["delta_vs_native"]["momentum_r"]["max_abs"] == 6.0
    # The geometric-source omission must not fabricate a rest-current source.
    assert controls["omitted_geometric_source"]["delta_vs_native"]["rest"]["max_abs"] == 0.0


def test_custom_face_fv_matches_merged_projection_reconstruction():
    n = 64
    dr = 0.125
    r = (np.arange(n, dtype=np.float64) + 0.5) * dr
    alpha = 1.0 + 0.001 * r**2
    gamma_rr = 1.0 + 0.0001 * r**2
    gamma_thth = r**2 * (1.0 + 0.0002 * r**2)
    rho = 0.5 + 0.01 * r**2
    j = 2.0e-5 * r * (1.0 + 0.01 * r**2)
    pr = 0.2 + 0.003 * r**2
    pt = 0.18 + 0.002 * r**2
    dt_density_t = -3.0e-5 * r**2 * (1.0 + 0.01 * r**2)
    dt_density_r = 1.0e-6 * r**3 * (1.0 + 0.02 * r**2)
    dt_gtt = -2.0e-4 * r**2
    dt_gamma_rr = 1.0e-4 * r**2
    dt_gamma_thth = 3.0e-4 * r**2

    def d1(values, parity):
        del parity
        return np.gradient(np.asarray(values, dtype=np.float64), dr, edge_order=2)

    flux_t_cells = -(alpha**2) * np.sqrt(gamma_rr) * gamma_thth * j / gamma_rr
    flux_r_cells = alpha * np.sqrt(gamma_rr) * gamma_thth * pr
    flux_t_faces = AUDIT_HLL.AUDIT.spherical_flux_face_values(
        flux_t_cells, dr, origin_power=3
    )
    flux_r_faces = AUDIT_HLL.AUDIT.spherical_flux_face_values(
        flux_r_cells, dr, origin_power=2
    )
    args = (
        r, d1, alpha, gamma_rr, gamma_thth, rho, j, pr, pt,
        dt_density_t, dt_density_r, dt_gtt, dt_gamma_rr, dt_gamma_thth,
    )

    expected_t, expected_r = AUDIT_HLL.AUDIT.covariant_divergence_fv_from_profiles(*args)
    actual_t, actual_r = AUDIT_HLL.covariant_fv_with_custom_flux_faces(
        *args, flux_t_faces, flux_r_faces
    )

    np.testing.assert_allclose(actual_t, expected_t, rtol=5.0e-12, atol=2.0e-18)
    np.testing.assert_allclose(actual_r, expected_r, rtol=5.0e-12, atol=2.0e-18)


def test_component_contract_rejects_bad_face_shape_and_nonfinite_source():
    face_flux = np.zeros((4, 3), dtype=np.float64)
    zeros = np.zeros(3, dtype=np.float64)
    with pytest.raises(ValueError):
        AUDIT_HLL.rhs_from_native_components(np.zeros((3, 2)), zeros, zeros, zeros, zeros, 1.0)
    invalid = zeros.copy()
    invalid[1] = np.nan
    with pytest.raises(FloatingPointError):
        AUDIT_HLL.rhs_from_native_components(face_flux, invalid, zeros, zeros, zeros, 1.0)
