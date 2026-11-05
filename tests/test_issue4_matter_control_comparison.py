import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from tools.issue4_matter_control_comparison import compare


def report(state_class, residual_offset=0.0):
    grid = np.arange(4, dtype=float)
    channels = {}
    summaries = {}
    for epsilon in ("eps2e-05", "eps1e-05"):
        channels[epsilon] = {
            "measured_tangent": grid + residual_offset,
            "pointwise": {
                "predictor": grid,
                "stress_source": grid,
                "residual_tangent_minus_predictor": grid + residual_offset,
            },
            "finite_volume": {
                "predictor": grid,
                "stress_source": grid,
                "residual_tangent_minus_predictor": grid + residual_offset,
            },
        }
        summaries[epsilon] = {}
        for name in ("pointwise", "finite_volume"):
            summaries[epsilon][name] = {
                "predictor": {"max_abs_r_le_10": 3.0, "radius_at_max_r_le_10": 2.0},
                "residual_tangent_minus_predictor": {
                    "max_abs_r_le_10": 3.0 + residual_offset,
                    "radius_at_max_r_le_10": 2.0,
                    "max_abs_annulus_5_to_10": 0.0,
                    "rms_annulus_5_to_10": 0.0,
                },
            }
    return {
        "status": "diagnostic_evaluation_complete",
        "state_class": state_class,
        "tangent_semantics": "test",
        "production_steps_run": 0,
        "state_advanced": False,
        "zero_step_tripwire": "armed_and_not_triggered",
        "channels": channels,
        "summaries": summaries,
    }


class MatterControlComparisonTests(unittest.TestCase):
    def test_three_class_comparison_and_directional_deltas(self):
        with tempfile.TemporaryDirectory() as d:
            paths = [Path(d) / f"{i}.json" for i in range(3)]
            for path, value in zip(paths, (
                report("manufactured_regular_family"),
                report("equation_consistent_matter_direction", 1e-12),
                report("accepted_t_0p60"),
            )):
                path.write_text(json.dumps(value))
            result = compare(*paths)
            self.assertEqual(result["status"], "comparison_complete_diagnostic_only")
            self.assertEqual(len(result["classes"]), 3)
            self.assertAlmostEqual(
                result["manufactured_direction_deltas"]["eps1e-05"]["pointwise"]
                ["residual_tangent_minus_predictor_max_abs_equation_minus_arbitrary"],
                1e-12,
            )

    def test_rejects_nonzero_production_steps(self):
        with tempfile.TemporaryDirectory() as d:
            paths = [Path(d) / f"{i}.json" for i in range(3)]
            values = [report("a"), report("b"), report("c")]
            values[1]["production_steps_run"] = 1
            for path, value in zip(paths, values):
                path.write_text(json.dumps(value))
            with self.assertRaisesRegex(RuntimeError, "zero-step"):
                compare(*paths)


if __name__ == "__main__":
    unittest.main()
