"""Unit tests for the zero-step observer's non-evolving foundation."""
from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import unittest

MODULE_PATH = Path(__file__).resolve().parents[1] / "tools" / "issue4_zero_step_full_bianchi_observer.py"
SPEC = importlib.util.spec_from_file_location("issue4_zero_step_observer", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
observer = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = observer
SPEC.loader.exec_module(observer)


class ObserverFoundationTests(unittest.TestCase):
    def test_digest_is_order_independent(self):
        import numpy as np
        a = {"geometry.a": np.array([1.0, 2.0]), "geometry.X": np.array([3.0])}
        b = {"geometry.X": np.array([3.0]), "geometry.a": np.array([1.0, 2.0])}
        self.assertEqual(observer.digest_named_arrays(a), observer.digest_named_arrays(b))

    def test_digest_changes_when_state_changes(self):
        import numpy as np
        before = observer.digest_named_arrays({"geometry.a": np.array([1.0, 2.0])})
        after = observer.digest_named_arrays({"geometry.a": np.array([1.0, 2.0001])})
        self.assertNotEqual(before, after)

    def test_digest_rejects_nonfinite_state(self):
        import numpy as np
        with self.assertRaises(ValueError):
            observer.digest_named_arrays({"geometry.a": np.array([1.0, float("nan")])})

    def test_unchanged_state_guard_accepts_identical_digest(self):
        observer.require_unchanged("abc", "abc")

    def test_unchanged_state_guard_rejects_mutation(self):
        with self.assertRaises(RuntimeError):
            observer.require_unchanged("before", "after")

    def test_no_production_kernel_import_at_module_load(self):
        self.assertNotIn("engine.production_kernel", sys.modules)


if __name__ == "__main__":
    unittest.main()
