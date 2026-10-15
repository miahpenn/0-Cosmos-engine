"""Regression tests for strict-JSON failure evidence in the PIRK L2 observer."""
import json

import numpy as np
import pytest

from tools.audit_tracefree_ricci_proper_radial import _write_strict_json_report, _independent_metric_crosscheck


def test_nonfinite_diagnostic_failure_still_writes_strict_json(tmp_path):
    report = {
        "status": "diagnostic_crosscheck_failure",
        "errors": ["N=80: nonfinite observer input"],
        "resolution_reports": [
            {
                "resolution": 80,
                "all_checks_pass": False,
                "stage_summaries": [
                    {
                        "stage": "final_implicit_l2_rhs",
                        "finite": False,
                        "independent_metric_crosscheck": {
                            "ricci_gap_max_abs_full": float("nan"),
                            "ricci_gap_peak_full": float("inf"),
                            "direct_d2_hessian_gap_max_abs_full": np.float64("-inf"),
                            "centre_witness": np.array([0.125, np.nan, -np.inf]),
                        },
                    }
                ],
            }
        ],
    }
    output = tmp_path / "report.json"

    _write_strict_json_report(output, report)

    raw = output.read_text(encoding="utf-8")
    decoded = json.loads(
        raw,
        parse_constant=lambda token: pytest.fail(f"non-standard JSON constant: {token}"),
    )
    stage = decoded["resolution_reports"][0]["stage_summaries"][0]
    assert decoded["status"] == "diagnostic_crosscheck_failure"
    assert decoded["errors"] == ["N=80: nonfinite observer input"]
    assert stage["finite"] is False
    metric = stage["independent_metric_crosscheck"]
    assert metric["ricci_gap_max_abs_full"] is None
    assert metric["ricci_gap_peak_full"] is None
    assert metric["direct_d2_hessian_gap_max_abs_full"] is None
    assert metric["centre_witness"] == [0.125, None, None]
    assert "NaN" not in raw
    assert "Infinity" not in raw


def test_manufactured_constant_curvature_regular_profile_is_negative_control():
    from engine.v55_pirk_adapter import vendor_modules
    grid_ops, vacuum, _ = vendor_modules()
    grid = grid_ops.SphericalCellGrid(160, 8.0)
    r = np.asarray(grid.centers, dtype=float)
    k = 1.0e-3
    f = 1.0 / np.sqrt(1.0-k*r*r)
    geometry = vacuum.flat_state(grid)
    geometry.a[:] = f**(4.0/3.0)
    geometry.b[:] = f**(-2.0/3.0)
    geometry.X[:] = f**(-1.0/3.0)
    geometry.alpha[:] = 1.0
    geometry.beta[:] = 0.0
    d1 = lambda values, parity: grid.cell_derivative_fourth(values, parity=parity)
    geometry.Lambda[:] = d1(geometry.a,1)/(2.0*geometry.a*geometry.a)-d1(geometry.b,1)/(geometry.a*geometry.b)+2.0/r*(1.0/geometry.b-1.0/geometry.a)
    terms = vacuum.geometry_terms(grid,geometry)
    capture = {"r":r,"a":geometry.a,"b":geometry.b,"X":geometry.X,"alpha":geometry.alpha,"Lambda":geometry.Lambda,"R":terms["R"],"Rr_r":terms["Rr_r"],"DrDr_alpha":terms["DrDr_alpha"],"lap_alpha":terms["lap_alpha"]}
    arrays, summary = _independent_metric_crosscheck(capture,float(grid.dr))
    interior = r <= 2.0
    assert np.max(np.abs(arrays["ricci_tracefree_metric"][interior])) < 5.0e-6
    assert np.max(np.abs(arrays["ricci_tracefree_vendor"][interior])) < 5.0e-6
    assert np.max(np.abs(arrays["ricci_tracefree_gap"][interior])) < 5.0e-6
    assert summary["ricci_gap_max_abs_centre_i0_i1"] == pytest.approx(np.max(np.abs(arrays["ricci_tracefree_gap"][:2])))
    annulus = summary["ricci_gap_fixed_areal_annulus_R5_R8"]
    assert annulus["count"] > 0 and annulus["max_abs"] is not None
