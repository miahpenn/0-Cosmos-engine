"""Regression tests for strict-JSON failure evidence in the PIRK L2 observer."""
import json

import numpy as np
import pytest

from tools.audit_current_main_pirk_l2_metric_crosscheck import _write_strict_json_report


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
