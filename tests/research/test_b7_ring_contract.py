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


def test_h1_keys_alone_never_unlock_adapter():
    keys_none = {k: None for k in mp.REQUIRED_DESTINATION}
    keys_empty = {k: [] for k in mp.REQUIRED_DESTINATION}
    full = {k: [1.0] for k in mp.REQUIRED_DESTINATION}
    allpass = {g: {"result": src.PASS} for g in ("G0", "G1", "G2", "G3", "G4", "G5", "G6")}
    untested = dict(allpass, G4={"result": src.UNTESTED})
    assert not mp.adapter_status(keys_none, allpass)["adapter_ready"]
    assert not mp.adapter_status(keys_empty, allpass)["adapter_ready"]
    assert not mp.adapter_status(full, untested)["adapter_ready"]
    assert not mp.adapter_status(full, None)["adapter_ready"]
    assert mp.adapter_status(full, allpass)["adapter_ready"]            # synthetic only


def _g0(tmp_path, code_ref, files=None, schema=None):
    root = tmp_path / "m"
    root.mkdir(parents=True)
    (root / "MANIFEST.sha256").write_text("")
    arch = tmp_path / "a.tgz"
    arch.write_bytes(b"x")
    npz = root / "d.npz"
    np.savez(npz, r=np.zeros(1), t300_s=np.zeros(1), t300_ps=np.zeros(1), t300_d=np.zeros(1), t300_pd=np.zeros(1))
    for n, c in (files or {}).items():
        (root / n).write_text(c)
    m = {"source": {"archive_sha256": src.sha256_file(arch), "data_member": "d.npz", "data_sha256": src.sha256_file(npz),
                    "code_ref": code_ref, "environment": {"python": "3", "numpy": "2", "scipy": "1"}},
         "slice": {"field_keys": {"r": "r", "s": "t300_s", "ps": "t300_ps", "d": "t300_d", "pd": "t300_pd"}}}
    g = src.g0_provenance(m, str(arch), str(root), schema_path=schema, expected_code_files=("a.py",))
    return {c["check"][:3]: c["result"] for c in g["checks"]}


def test_h2_code_ref_is_fail_closed(tmp_path):
    h = src.sha256_bytes(b"print(1)")
    ok = _g0(tmp_path / "1", "a.py sha256=" + h, {"a.py": "print(1)"})
    assert ok["G0e"] == src.PASS
    assert _g0(tmp_path / "2", "", {"a.py": "print(1)"})["G0e"] == src.FAIL                       # zero refs
    assert _g0(tmp_path / "3", "a.py sha256=" + h, {})["G0e"] == src.FAIL                          # member missing
    assert _g0(tmp_path / "4", "a.py sha256=" + h + "; a.py sha256=" + h, {"a.py": "print(1)"})["G0e"] == src.FAIL   # duplicate
    assert _g0(tmp_path / "5", "a.py sha256=" + h[:12], {"a.py": "print(1)"})["G0e"] == src.FAIL   # short digest
    assert _g0(tmp_path / "6", "a.py sha256=" + h + "; dep sha256=" + h, {"a.py": "print(1)"})["G0e"] == src.FAIL  # unparsed ref
    assert _g0(tmp_path / "7", "b.py sha256=" + h, {"b.py": "print(1)"})["G0e"] == src.FAIL        # expected a.py not cited
    assert _g0(tmp_path / "8", "../a.py sha256=" + h, {"a.py": "print(1)"})["G0e"] == src.FAIL     # unsafe name


def test_h3_schema_enforced_and_unsupplied_is_untested(tmp_path):
    h = src.sha256_bytes(b"print(1)")
    sch = os.path.join("research", "b7_ring_schemas", "b7_ring_source_v1.schema.json")
    r = _g0(tmp_path / "1", "a.py sha256=" + h, {"a.py": "print(1)"}, schema=sch)
    assert r["G0g"] == src.FAIL                                    # minimal manifest lacks required schema fields
    assert _g0(tmp_path / "2", "a.py sha256=" + h, {"a.py": "print(1)"})["G0g"] == src.UNTESTED


def _h5(tmp_path, data_member, in_manifest=True, symlink_to=None, bad_manifest_name=False):
    root = tmp_path / "m"
    root.mkdir(parents=True)
    arch = tmp_path / "a.tgz"
    arch.write_bytes(b"x")
    outside = tmp_path / "outside.npz"
    z = dict(r=np.zeros(1), t300_s=np.zeros(1), t300_ps=np.zeros(1), t300_d=np.zeros(1), t300_pd=np.zeros(1))
    np.savez(outside, **z)
    inner = root / "d.npz"
    np.savez(inner, **z)
    if symlink_to:
        (root / "link.npz").symlink_to(tmp_path / symlink_to)
    lines = ["%s  d.npz" % src.sha256_file(inner)] if in_manifest else []
    if bad_manifest_name:
        lines.append("%s  ../outside.npz" % src.sha256_file(outside))
    (root / "MANIFEST.sha256").write_text("\n".join(lines) + "\n")
    ref = data_member if not os.path.isabs(data_member) else data_member
    target = (root / "d.npz") if data_member in ("d.npz", "link.npz") else outside
    m = {"source": {"archive_sha256": src.sha256_file(arch), "data_member": ref, "data_sha256": src.sha256_file(target),
                    "code_ref": "", "environment": {"python": "3", "numpy": "2", "scipy": "1"}},
         "slice": {"field_keys": {"r": "r", "s": "t300_s", "ps": "t300_ps", "d": "t300_d", "pd": "t300_pd"}}}
    sch = os.path.join("research", "b7_ring_schemas", "b7_ring_source_v1.schema.h5.json")
    g = src.g0_provenance(m, str(arch), str(root), schema_path=sch)
    r = {c["check"][:3]: c["result"] for c in g["checks"]}
    return g["result"], r, str(outside)


def test_h5_traversal_absolute_and_nested_names_fail_g0(tmp_path):
    for i, name in enumerate(["../outside.npz", "nested/../d.npz", "nested/d.npz", "..", ".hidden.npz"]):
        res, r, _ = _h5(tmp_path / str(i), name)
        assert res == src.FAIL and r["G0c"] == src.FAIL and r["G0h"] == src.FAIL, (name, r)
        assert r["G0g"] == src.FAIL, (name, r)
    p = str(tmp_path / "x" / "outside.npz")
    res, r, out = _h5(tmp_path / "x", p)
    assert res == src.FAIL and r["G0c"] == src.FAIL and r["G0h"] == src.FAIL


def test_h5_member_must_be_in_internal_manifest(tmp_path):
    res, r, _ = _h5(tmp_path / "a", "d.npz", in_manifest=True)
    assert r["G0c"] == src.PASS and r["G0h"] == src.PASS
    res, r, _ = _h5(tmp_path / "b", "d.npz", in_manifest=False)
    assert r["G0c"] == src.PASS and r["G0h"] == src.FAIL and res == src.FAIL


def test_h5_symlink_escaping_the_folder_fails(tmp_path):
    res, r, _ = _h5(tmp_path / "s", "link.npz", symlink_to="outside.npz")
    assert r["G0c"] == src.FAIL and res == src.FAIL


def test_h5_unsafe_name_inside_manifest_fails_g0b(tmp_path):
    res, r, _ = _h5(tmp_path / "u", "d.npz", bad_manifest_name=True)
    assert r["G0b"] == src.FAIL and res == src.FAIL


def test_h5_outside_file_is_never_read(tmp_path, monkeypatch):
    seen = []
    real = np.load
    monkeypatch.setattr(np, "load", lambda p, *a, **k: (seen.append(str(p)), real(p, *a, **k))[1])
    res, r, outside = _h5(tmp_path, "../outside.npz")
    assert res == src.FAIL and not any(x.endswith("outside.npz") for x in seen)


ARCH = all(k in os.environ for k in ("B7_V7_TGZ", "B7_V5_TGZ", "B7_MANIFEST_C_R0", "B7_MANIFEST_C_R1"))


@pytest.mark.skipif(not ARCH, reason="archives/manifests not supplied (B7_V7_TGZ, B7_V5_TGZ, B7_MANIFEST_C_R0, B7_MANIFEST_C_R1)")
@pytest.mark.parametrize("rev", ["R0", "R1"])
def test_real_sources_expected_outcomes(tmp_path, rev):
    from research import b7_ring_run_v1 as run
    res = run.main(os.environ["B7_V7_TGZ"], os.environ["B7_V5_TGZ"], os.environ["B7_MANIFEST_C_" + rev], str(tmp_path), rev.lower())
    c = res["sources"]["C"]["gates"]
    if rev == "R0":                                     # historical: only the code_ref check fails
        assert c["G0"]["result"] == src.FAIL
        assert [x["check"][:3] for x in c["G0"]["checks"] if x["result"] == src.FAIL] == ["G0e"]
        assert c["G2"]["result"] == src.UNTESTED and "blocked" in c["G2"]["checks"][0]["check"]
    else:
        assert c["G0"]["result"] == src.PASS and c["G1"]["result"] == src.PASS
    for k in "ABC":
        g = res["sources"][k]["gates"]
        assert [x for x in g["G0"]["checks"] if x["check"].startswith("G0h")][0]["result"] == src.PASS
        assert g["G4"]["result"] != src.PASS and not res["sources"][k]["adapter"]["adapter_ready"]
    for k in "AB":
        assert res["sources"][k]["gates"]["G0"]["result"] == src.FAIL and "G2" in res["sources"][k]["gates"]
        assert "blocked" in res["sources"][k]["gates"]["G2"]["checks"][0]["check"]
    assert all(x["as_expected"] for x in res["controls"])
