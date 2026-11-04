"""Unit tests for the Issue #4 zero-step observer foundation."""
from __future__ import annotations

import ast
import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

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
        tree = ast.parse(MODULE_PATH.read_text(encoding="utf-8"))
        imported = []
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom) and (node.module or "").startswith("engine"):
                imported.append(node.module)
            elif isinstance(node, ast.Import):
                imported.extend(alias.name for alias in node.names if alias.name.startswith("engine"))
        self.assertEqual(imported, [])

    def test_provenance_accepts_only_the_two_distinct_source_checkouts(self):
        with tempfile.TemporaryDirectory() as tmp:
            diagnostic = Path(tmp) / "diagnostic"
            baseline = Path(tmp) / "baseline"
            vendor_root = diagnostic / "vendor" / "bb-palatini-unified-r0"
            vendor_root.mkdir(parents=True)
            baseline.mkdir()
            expected = {
                ("HEAD",): observer.ACCEPTED_SOURCE,
                ("HEAD:tools/audit_current_source_momentum_constraint_tangent.py",): observer.TANGENT_HELPER_BLOB,
                ("HEAD:tools/audit_total_stress_energy_divergence.py",): observer.TOTAL_STRESS_HELPER_BLOB,
            }
            expected.update({
                (f"HEAD:{path}",): blob
                for path, blob in observer.REQUIRED_SOURCE_BLOBS.items()
            })

            def fake_git(root, *args):
                if root == baseline and args == ("rev-parse", "HEAD"):
                    return observer.HELPER_BASELINE
                if root == diagnostic and args == ("rev-parse", "HEAD"):
                    return observer.ACCEPTED_SOURCE
                if root == vendor_root and args == ("rev-parse", "HEAD"):
                    return observer.VENDOR
                if root == diagnostic and args[:1] == ("rev-parse",):
                    return expected.get((args[1],), "unexpected")
                raise AssertionError(f"unexpected git probe: {root} {args}")

            with patch.object(observer, "_git", side_effect=fake_git):
                manifest = observer.verify_provenance(diagnostic, baseline)
            self.assertEqual(manifest["accepted_diagnostic_source_commit"], observer.ACCEPTED_SOURCE)
            self.assertEqual(manifest["protected_helper_baseline_commit"], observer.HELPER_BASELINE)
            self.assertEqual(manifest["vendor_commit"], observer.VENDOR)
            self.assertEqual(manifest["accepted_source_blobs"], observer.REQUIRED_SOURCE_BLOBS)

    def test_each_pinned_blob_and_vendor_mismatch_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            diagnostic = Path(tmp) / "diagnostic"
            baseline = Path(tmp) / "baseline"
            vendor_root = diagnostic / "vendor" / "bb-palatini-unified-r0"
            vendor_root.mkdir(parents=True)
            baseline.mkdir()

            probes = [
                (diagnostic, ("rev-parse", f"HEAD:{path}"))
                for path in observer.REQUIRED_SOURCE_BLOBS
            ]
            probes.extend([
                (diagnostic, ("rev-parse", "HEAD:tools/audit_current_source_momentum_constraint_tangent.py")),
                (diagnostic, ("rev-parse", "HEAD:tools/audit_total_stress_energy_divergence.py")),
                (vendor_root, ("rev-parse", "HEAD")),
            ])

            for bad_root, bad_args in probes:
                with self.subTest(root=str(bad_root), args=bad_args):
                    def fake_git(root, *args):
                        if root == baseline and args == ("rev-parse", "HEAD"):
                            return observer.HELPER_BASELINE
                        if root == diagnostic and args == ("rev-parse", "HEAD"):
                            return observer.ACCEPTED_SOURCE
                        if root == vendor_root and args == ("rev-parse", "HEAD"):
                            return "wrong-vendor" if (root, args) == (bad_root, bad_args) else observer.VENDOR
                        if root == diagnostic and args[:1] == ("rev-parse",):
                            key = (root, args)
                            if key == (bad_root, bad_args):
                                return "wrong-blob"
                            name = args[1]
                            if name.startswith("HEAD:tools/audit_current_source_momentum_constraint_tangent.py"):
                                return observer.TANGENT_HELPER_BLOB
                            if name.startswith("HEAD:tools/audit_total_stress_energy_divergence.py"):
                                return observer.TOTAL_STRESS_HELPER_BLOB
                            source_path = name.removeprefix("HEAD:")
                            return observer.REQUIRED_SOURCE_BLOBS.get(source_path, "unexpected")
                        raise AssertionError(f"unexpected git probe: {root} {args}")

                    with patch.object(observer, "_git", side_effect=fake_git):
                        with self.assertRaises(RuntimeError):
                            observer.verify_provenance(diagnostic, baseline)

    def test_historical_baseline_cannot_substitute_for_accepted_source(self):
        with tempfile.TemporaryDirectory() as tmp:
            diagnostic = Path(tmp) / "diagnostic"
            baseline = Path(tmp) / "baseline"
            diagnostic.mkdir()
            baseline.mkdir()

            def fake_git(root, *args):
                if root == diagnostic and args == ("rev-parse", "HEAD"):
                    return observer.HELPER_BASELINE
                if root == baseline and args == ("rev-parse", "HEAD"):
                    return observer.HELPER_BASELINE
                raise AssertionError(f"unexpected git probe: {root} {args}")

            with patch.object(observer, "_git", side_effect=fake_git):
                with self.assertRaisesRegex(RuntimeError, "accepted diagnostic source commit mismatch"):
                    observer.verify_provenance(diagnostic, baseline)


if __name__ == "__main__":
    unittest.main()
