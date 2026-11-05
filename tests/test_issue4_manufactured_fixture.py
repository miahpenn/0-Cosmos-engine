"""Fail-closed schema tests for the manufactured Issue #4 fixture loader."""
from __future__ import annotations

import tempfile
import unittest
import zipfile
from pathlib import Path

from tools.issue4_manufactured_fixture import load_manufactured_fixture


class ManufacturedFixtureLoaderTests(unittest.TestCase):
    def test_rejects_missing_archive_members_before_state_loading(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "incomplete.zip"
            with zipfile.ZipFile(path, "w") as archive:
                archive.writestr("report.json", "{}")
            with self.assertRaisesRegex(RuntimeError, "member set mismatch"):
                load_manufactured_fixture(path)


if __name__ == "__main__":
    unittest.main()
