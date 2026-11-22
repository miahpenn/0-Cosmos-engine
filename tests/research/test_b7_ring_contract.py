"""RESEARCH ONLY. B7 ring contract tests: deterministic positive and negative controls (see research/b7_ring_*).

Synthetic tests need no archive. The archive tests need B7_V7_TGZ, B7_V5_TGZ and B7_MANIFEST_C and are skipped otherwise.
Expected failing gates were fixed before the run (coordination log, 2026-11-22).
"""
import copy
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, ".")
from research import b7_ring_constraint_v1 as con  # noqa: E402
from research import b7_ring_mapping_v1 as mp  # noqa: E402
from research import b7_ring_source_v1 as src  # noqa: E402

DR, N = 0.05, 2400


def man():
    return {"slice": {"time": 300, "radial_grid": {"n_cells": N, "dr": DR, "r_max": 120}}}


def synth():
    r = (np.arange(N) + 0.5) * DR
    s = 0.3 * np.exp(-r ** 2 / 9)
    d = np.sqrt(2) - 0.5 * np.exp(-r ** 2 / 9)
    return {"r": r, "s": s, "d": d, "ps": 0.5 * s, "pd": 0.1 * np.sin(r)}


def gates(g):
    return {c["check"][:3]: c["result"] for c in g["checks"]}


def test_g1_passes_clean_synthetic_without_ledger():
    g = src.g1_snapshot(synth(), man(), None)
    r = gates(g)
    assert r["G1a"] == r["G1b"] == r["G1c"] == src.PASS and r["G1d"] == src.UNTESTED and g["result"] == src.UNTESTED


def test_g1_rejects_nan_nonuniform_missing_and_bare():
    f = synth(); f["s"][5] = np.nan
    assert gates(src.g1_snapshot(f, man()))["G1a"] == src.FAIL
    f = synth(); f["r"][10] += 1e-6
    assert gates(src.g1_snapshot(f, man()))["G1b"] == src.FAIL
    f = synth(); f["ps"][:] = 0
    assert gates(src.g1_snapshot(f, man()))["G1c"] == src.FAIL
    f = synth(); f["d"] = f["d"][:-1]
    assert gates(src.g1_snapshot(f, man()))["G1a"] == src.FAIL


def test_energy_ledger_identity_and_cross_term_sensitivity():
    f = synth()
    e = src.energy_ledger(f["s"], f["d"], f["ps"], f["pd"], DR)
    assert mp.g2_mapping(f, DR, e)["result"] == src.UNTESTED          # G2a PASS, G2b UNTESTED
    assert gates(mp.g2_mapping(f, DR, e))["G2a"] == src.PASS
    assert gates(mp.g2_mapping(f, DR, e, mapping_cross=False))["G2a"] == src.FAIL


def test_momentum_density_needs_ps_term():
    f = synth()
    assert gates(con.g4_gr_initial_data(f, DR))["G4a"] == src.PASS
    assert gates(con.g4_gr_initial_data(f, DR, con.momentum_density(f, DR, include_ps=False)))["G4a"] == src.FAIL
    assert con.g4_gr_initial_data(f, DR)["result"] != src.PASS          # never PASS without a solver path


def test_boundary_and_population_controls():
    f = synth()
    g = mp.g5_boundary(f, DR, {"open_edge": {"kind": "lattice_wall", "reflecting": True}}, True)
    assert gates(g)["G5b"] == src.FAIL
    assert mp.g3_population(None, True)["result"] == src.FAIL
    assert mp.g3_population(None, False)["result"] == src.UNTESTED


def test_adapter_not_ready_without_destination():
    a = mp.adapter_status(None)
    assert a["adapter_ready"] is False and set(a["missing_inputs"]) == set(mp.REQUIRED_DESTINATION)


ARCH = all(k in os.environ for k in ("B7_V7_TGZ", "B7_V5_TGZ", "B7_MANIFEST_C"))


@pytest.mark.skipif(not ARCH, reason="source archives not supplied (B7_V7_TGZ, B7_V5_TGZ, B7_MANIFEST_C)")
def test_real_sources_expected_outcomes(tmp_path):
    from research import b7_ring_run_v1 as run
    res = run.main(os.environ["B7_V7_TGZ"], os.environ["B7_V5_TGZ"], os.environ["B7_MANIFEST_C"], str(tmp_path))
    c = res["sources"]["C"]["gates"]
    assert c["G0"]["result"] == src.FAIL and [x for x in c["G0"]["checks"] if x["result"] == src.FAIL][0]["check"].startswith("G0e")
    assert res["sources"]["C"]["proposed_corrected_manifest_G0"]["result"] == src.PASS
    for k in "ABC":
        g = res["sources"][k]["gates"]
        assert g["G1"]["result"] == src.PASS and g["G4"]["result"] != src.PASS and not res["sources"][k]["adapter"]["adapter_ready"]
    assert all(x["as_expected"] for x in res["controls"])
