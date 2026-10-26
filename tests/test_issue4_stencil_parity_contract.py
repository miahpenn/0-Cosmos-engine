"""Fast operator/parity controls for the Issue #4 manufactured-profile audit."""
from types import SimpleNamespace
import json
import sys

import numpy as np

from engine.native_geometry import NativeSphericalCellGrid
from tools import audit_issue4_rhs_stencil_parity as audit


def test_vendor_and_native_derivative_operators_match_on_declared_parity_profiles():
    r0, _, _ = audit.adapter.vendor_modules()
    vendor_grid = r0.SphericalCellGrid(64, 40.0)
    native_grid = NativeSphericalCellGrid(64, 40.0)
    r = vendor_grid.centers

    even_profile = 1.0 + 0.05 * (r / 40.0) ** 2 + 0.01 * (r / 40.0) ** 4
    odd_profile = (r / 40.0) + 0.05 * (r / 40.0) ** 3
    profiles = ((even_profile, +1), (odd_profile, -1))
    for values, parity in profiles:
        np.testing.assert_allclose(
            audit.d1(vendor_grid, values, parity),
            audit.d1(native_grid, values, parity),
            rtol=0.0,
            atol=2.0e-14,
        )
        np.testing.assert_allclose(
            audit.d2(vendor_grid, values, parity),
            audit.d2(native_grid, values, parity),
            rtol=0.0,
            atol=2.0e-14,
        )


def test_production_and_native_momentum_formulas_match_in_flat_regular_chart():
    r0, _, _ = audit.adapter.vendor_modules()
    vendor_grid = r0.SphericalCellGrid(80, 8.0)
    native_grid = NativeSphericalCellGrid(80, 8.0)
    r = vendor_grid.centers

    c = 0.025
    k0 = 0.04
    k = 0.003
    state = SimpleNamespace(
        a=np.ones_like(r),
        b=np.ones_like(r),
        X=np.ones_like(r),
        alpha=np.ones_like(r),
        Lambda=np.zeros_like(r),
        Aa=c * r**2,
        K=k0 + k * r**2,
    )
    production = sum(
        audit.momentum_pieces_production(vendor_grid, state).values(),
        np.zeros_like(r),
    )
    native = sum(
        audit.momentum_pieces_native(native_grid, state).values(),
        np.zeros_like(r),
    )
    mask = (r >= 1.0) & (r <= 6.0)
    expected = (5.0 * c - (4.0 / 3.0) * k) * r
    np.testing.assert_allclose(production[mask], expected[mask], rtol=1.0e-12, atol=1.0e-12)
    np.testing.assert_allclose(native[mask], expected[mask], rtol=1.0e-12, atol=1.0e-12)
    np.testing.assert_allclose(production[mask], native[mask], rtol=1.0e-12, atol=1.0e-12)


def test_negative_operator_parity_is_saved_as_a_diagnostic_result(monkeypatch, tmp_path, capsys):
    report = {
        "schema": "issue4_manufactured_profile_rhs_stencil_parity_v1",
        "status": "diagnostic_only_not_scientific_admission",
        "provenance": {
            "production_source_commit": "source-sha",
            "production_source_vendor_commit": "vendor-sha",
        },
        "input_state": {"immutable": True},
        "decision": {
            "shared_operator_parity_check": {
                "passed": False,
                "checks": {"D1_Aa_even": {"passed": False}},
            },
            "production_constraint_decomposition": {"passed": True},
            "native_constraint_decomposition": {"passed": True},
            "current_projection_reconstruction": {"passed": True},
            "primary_l3_stress_source_reconstruction": {
                "Aa": {"passed": True},
                "K": {"passed": True},
            },
        },
    }
    arrays = {
        "r": np.array([0.5, 1.5]),
        "signed_witness": np.array([1.25e-6, -2.5e-7]),
    }
    monkeypatch.setattr(audit, "run_audit", lambda resolution, r_max: (report, arrays))
    monkeypatch.setattr(
        sys, "argv", ["audit_issue4_rhs_stencil_parity.py", "--output-dir", str(tmp_path)]
    )

    assert audit.main() == 0
    persisted = json.loads((tmp_path / "report.json").read_text(encoding="utf-8"))
    assert persisted["decision"]["shared_operator_parity_check"]["passed"] is False
    assert persisted["diagnostic_integrity"]["passed"] is True
    assert persisted["diagnostic_integrity"]["failed_checks"] == []
    with np.load(tmp_path / "signed_arrays.npz") as saved:
        np.testing.assert_array_equal(saved["signed_witness"], arrays["signed_witness"])
    manifest = json.loads((tmp_path / "SHA256SUMS.json").read_text(encoding="utf-8"))
    assert manifest["report.json"]
    assert manifest["signed_arrays.npz"] == persisted["arrays"]["sha256"]
    assert "shared_operator_parity_pass\": false" in capsys.readouterr().out

