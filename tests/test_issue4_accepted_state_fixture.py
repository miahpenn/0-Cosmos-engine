"""Tests for the checksum-verified Issue #4 accepted-state fixture loader."""
from __future__ import annotations

import hashlib
import tempfile
import unittest
from pathlib import Path

import numpy as np

from tools.issue4_accepted_state_fixture import (
    INPUT_ARRAY_KEYS,
    STATE_ARRAY_NAMES,
    helper_state_digest,
    load_accepted_fixture,
    structural_inventory,
)


class AcceptedStateFixtureTests(unittest.TestCase):
    def test_helper_digest_matches_pinned_byte_convention(self):
        arrays = {
            name: np.asarray([float(i), float(i + 1)], dtype=np.float64)
            for i, name in enumerate(STATE_ARRAY_NAMES)
        }
        expected = hashlib.sha256()
        for name in STATE_ARRAY_NAMES:
            value = np.ascontiguousarray(arrays[name], dtype=np.float64)
            expected.update(name.encode("utf-8"))
            expected.update(str(value.shape).encode("ascii"))
            expected.update(value.tobytes(order="C"))
        self.assertEqual(helper_state_digest(arrays), expected.hexdigest())

    def test_helper_digest_refuses_missing_and_nonfinite_fields(self):
        with self.assertRaises(KeyError):
            helper_state_digest({})
        arrays = {
            name: np.asarray([1.0], dtype=np.float64)
            for name in STATE_ARRAY_NAMES
        }
        arrays[STATE_ARRAY_NAMES[0]] = np.asarray([float("nan")])
        with self.assertRaises(ValueError):
            helper_state_digest(arrays)

    def test_0star_inventory_keeps_tangent_and_divergence_channels_separate(self):
        paths = set(INPUT_ARRAY_KEYS.values())
        paths.update({
            "tangent_eps2e-05.Mdot",
            "tangent_eps1e-05.Mdot",
            "tangent_eps2e-05.Mdot_geom",
            "tangent_eps1e-05.Mdot_geom",
            "total_T_eps2e-05.total.div_T_r",
            "total_T_eps1e-05.total.div_T_r",
            "total_T_eps2e-05.total.div_T_r_fv",
            "total_T_eps1e-05.total.div_T_r_fv",
        })
        report = {"raw_arrays": {"array_records": [
            {"path": path} for path in sorted(paths)
        ]}}
        inventory = structural_inventory({}, report)
        self.assertFalse(any(inventory["missing_required_witnesses"].values()))
        self.assertEqual(
            inventory["stored_operand_families"]["pointwise_total_stress_divergence"],
            ["total_T_eps1e-05.total.div_T_r", "total_T_eps2e-05.total.div_T_r"],
        )
        self.assertEqual(
            inventory["stored_operand_families"]["finite_volume_total_stress_divergence"],
            ["total_T_eps1e-05.total.div_T_r_fv", "total_T_eps2e-05.total.div_T_r_fv"],
        )
        self.assertEqual(
            set(inventory["predictor_terms_requiring_reconstruction"]),
            {
                "shift_transport",
                "normal_alpha_K_M",
                "normal_lapse_gradient_H",
                "normal_H_gradient",
                "evolved_connection_tracefree_divergence",
                "total_stress_bianchi_source",
            },
        )

    def test_archive_loader_rejects_nonmatching_archive_digest_before_parsing(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "not-the-accepted-fixture.zip"
            path.write_bytes(b"not an accepted archive")
            with self.assertRaisesRegex(RuntimeError, "archive SHA-256 mismatch"):
                load_accepted_fixture(path, 60)


if __name__ == "__main__":
    unittest.main()
