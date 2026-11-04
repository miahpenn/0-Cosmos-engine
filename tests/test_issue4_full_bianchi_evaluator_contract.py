"""Static zero-step contract tests for the canonical Bianchi evaluator."""
from __future__ import annotations

import ast
import sys
import types
from pathlib import Path
from unittest.mock import patch
import unittest

SOURCE = Path(__file__).resolve().parents[1] / "tools" / "issue4_full_bianchi_evaluator.py"

import importlib.util
SPEC = importlib.util.spec_from_file_location("issue4_full_bianchi_evaluator_under_test", SOURCE)
assert SPEC is not None and SPEC.loader is not None
evaluator = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = evaluator
SPEC.loader.exec_module(evaluator)


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

    def test_vendor_loader_rejects_preloaded_module_outside_pinned_checkout(self):
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            vendor = Path(tmp) / "vendor" / "bb-palatini-unified-r0"
            vendor.mkdir(parents=True)
            poisoned = types.ModuleType("r0_operators")
            poisoned.__file__ = str(Path(tmp) / "untrusted" / "r0_operators.py")
            with patch.dict(sys.modules, {"r0_operators": poisoned}):
                with self.assertRaisesRegex(RuntimeError, "vendor module origin mismatch"):
                    evaluator._load_vendor(Path(tmp))

    def test_zero_step_guard_restores_methods_after_success_and_exception(self):
        class Kernel:
            def step(self):
                return "original-kernel"

        class Vendor:
            def damped_moving_puncture_pirk2_step(self):
                return "original-vendor"

        original_kernel = Kernel.step
        original_vendor = Vendor.damped_moving_puncture_pirk2_step
        calls = []
        with evaluator._zero_step_guard(Kernel, Vendor, calls):
            self.assertIsNot(Kernel.step, original_kernel)
            self.assertIsNot(Vendor.damped_moving_puncture_pirk2_step, original_vendor)
            with self.assertRaisesRegex(RuntimeError, "forbids ProductionKernel.step"):
                Kernel().step()
        self.assertIs(Kernel.step, original_kernel)
        self.assertIs(Vendor.damped_moving_puncture_pirk2_step, original_vendor)
        self.assertEqual(calls, ["ProductionKernel.step"])

        with self.assertRaisesRegex(ValueError, "injected"):
            with evaluator._zero_step_guard(Kernel, Vendor, calls):
                raise ValueError("injected")
        self.assertIs(Kernel.step, original_kernel)
        self.assertIs(Vendor.damped_moving_puncture_pirk2_step, original_vendor)

    def test_output_serialization_uses_real_newline(self):
        source = SOURCE.read_text(encoding="utf-8")
        self.assertEqual(source.count('write_text(payload + "\\\\n"'), 0)
        self.assertEqual(source.count('write_text(payload + "\\n"'), 2)


if __name__ == "__main__":
    unittest.main()
