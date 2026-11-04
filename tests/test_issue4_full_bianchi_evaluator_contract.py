"""Static zero-step contract tests for the canonical Bianchi evaluator."""
from __future__ import annotations

import ast
from pathlib import Path
import unittest

SOURCE = Path(__file__).resolve().parents[1] / "tools" / "issue4_full_bianchi_evaluator.py"


class FullBianchiEvaluatorContractTests(unittest.TestCase):
    def test_no_step_or_pirk_advance_call_is_present(self):
        tree = ast.parse(SOURCE.read_text(encoding="utf-8"))
        forbidden = []
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            func = node.func
            if isinstance(func, ast.Attribute) and func.attr in {
                "step", "damped_moving_puncture_pirk2_step"
            }:
                forbidden.append((func.attr, node.lineno))
        self.assertEqual(forbidden, [])

    def test_runtime_engine_import_is_deferred_until_after_provenance_gate(self):
        source = SOURCE.read_text(encoding="utf-8")
        self.assertNotIn("from engine.", source)
        self.assertNotIn("import engine.", source)
        self.assertLess(
            source.index("provenance = verify_provenance"),
            source.index("fixture = load_accepted_fixture"),
        )

    def test_canonical_predictor_keeps_connection_parity_and_chi_distinctions(self):
        source = SOURCE.read_text(encoding="utf-8")
        self.assertIn("d1(C_r, parity=-1)", source)
        self.assertIn("chi = X**2", source)
        self.assertIn("-d1(X, parity=1) / (2.0 * X)", source)
        self.assertIn("d1(F_C, parity=1)", source)
        self.assertIn("vacuum_bssn.geometry_terms(grid, geom)", source)
        self.assertIn("ricci_offset_crosscheck", source)


if __name__ == "__main__":
    unittest.main()
