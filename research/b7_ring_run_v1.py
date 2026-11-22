"""RESEARCH ONLY. B7 harness driver: runs G0-G7 per source and the pre-declared negative controls; writes results.

Usage: python research/b7_ring_run_v1.py V7.tgz V5.tgz MANIFEST_C.json OUTDIR
Read-only on its inputs; extracts into temporary directories; writes only under OUTDIR. No engine run, no G7.
"""
import copy
import json
import os
import sys
import tempfile

sys.path.insert(0, ".")
import numpy as np  # noqa: E402

from research import b7_ring_constraint_v1 as con  # noqa: E402
from research import b7_ring_mapping_v1 as mp  # noqa: E402
from research import b7_ring_source_v1 as src  # noqa: E402

V5_SHA = "dd4f468edbd5438229e02af189cc89f9265cbb2477b71af163ef8c4b0374c9f4"
# Projects-built manifests for A and B from B7 design section 2 (declared values only; no source environment recorded)
V5_SOURCES = {
    "A": ("run_w1917_dr0025.npz", "4b3155a5ad423c18d3cf622f34882f96f6cf7fd536b19053aa497bcd76b93df9", 4800, 0.025),
    "B": ("run_w1877_dr005.npz", "af0b765302a57fb36f2d1e3d3e4e0bf53283a565267e66647b9e7615ef13ffc6", 2400, 0.05),
}


def v5_manifest(member, sha, n, dr):
    return {"schema_version": "b7-ring-source-v1", "claim": "candidate_not_admitted",
            "source": {"archive_name": "ringkit_v5.tgz", "archive_sha256": V5_SHA, "manifest_member": "MANIFEST.sha256",
                       "data_member": member, "data_sha256": sha, "run_label": "Projects-built from B7 design section 2",
                       "code_ref": "", "environment": "unknown"},
            "slice": {"time": 300, "time_unit": "mu^-1", "phase_basis": "stored PDE slice t300",
                      "radial_grid": {"n_cells": n, "dr": dr, "r_max": 120, "radius_unit": "mu^-1", "centering": "cell-centred"},
                      "field_keys": {"r": "r", "s": "t300_s", "ps": "t300_ps", "d": "t300_d", "pd": "t300_pd"}}}


# Ledger column of E(<15) in rows: v7 evo.py:31 comment (t, s0, d@6 probes, E..., so col 8); v5 ringrun.py rows
# (t, s0, d@6 probes, s@2 probes, E...) per ringkit README.txt:3 and analyze_runs.py:1, so col 10.
LEDGER_COL = {"C": 8, "A": 10, "B": 10}


def run_source(label, manifest, archive, root):
    keys = manifest["slice"]["field_keys"]
    rep = {"source": label, "data_member": manifest["source"]["data_member"], "data_sha256": manifest["source"]["data_sha256"]}
    g0 = src.g0_provenance(manifest, archive, root)
    dp = os.path.join(root, manifest["source"]["data_member"])
    f, d = src.load_slice(dp, keys)
    dr = manifest["slice"]["radial_grid"]["dr"]
    row = src.ledger_row_at(d["rows"], manifest["slice"]["time"], lo=LEDGER_COL[label])
    gates = {"G0": g0, "G1": src.g1_snapshot(f, manifest, row), "G2": mp.g2_mapping(f, dr, row),
             "G3": mp.g3_population(), "G4": con.g4_gr_initial_data(f, dr), "G5": mp.g5_boundary(f, dr),
             "G6": {"gate": "G6", "result": src.UNTESTED, "checks": [{"check": "same-source dressed BVP at measured omega", "result": src.UNTESTED,
                    "detail": "not implemented; dressed energies 176.38 vs 181.69 (E-matched omega) recorded as not closed by Aster"}]},
             "G7": {"gate": "G7", "result": src.UNTESTED, "checks": [{"check": "dynamics", "result": src.UNTESTED, "detail": "not authorised"}]}}
    rep["gates"] = gates
    rep["adapter"] = mp.adapter_status()
    rep["fingerprint_sha256"] = src.fingerprint(f)
    rep["ledger_row_t300"] = row
    return rep, f, d, row, dr


def controls(manifest, archive, root, f, row, dr):
    out = []

    def rec(name, expect, got, ok, note=""):
        out.append({"control": name, "expected_failing": expect, "observed": got, "as_expected": ok, "note": note})

    m1 = copy.deepcopy(manifest)
    s = m1["source"]["archive_sha256"]
    m1["source"]["archive_sha256"] = s[:-1] + ("0" if s[-1] != "0" else "1")
    g = src.g0_provenance(m1, archive, root)
    a = [c for c in g["checks"] if c["check"].startswith("G0a")][0]["result"]
    rec("c1 altered source sha", "G0a", {"G0a": a}, a == src.FAIL)
    f2 = {k: v.copy() for k, v in f.items()}
    f2["pd"][:] = 0
    g = src.g1_snapshot(f2, manifest, row)
    r = {c["check"][:3]: c["result"] for c in g["checks"]}
    rec("c2 pd stripped", "G1c and G1d", r, r["G1c"] == src.FAIL and r["G1d"] == src.FAIL)
    f3 = {k: v.copy() for k, v in f.items()}
    f3["pd"] = -f3["pd"]
    g1 = src.g1_snapshot(f3, manifest, row)
    g5 = mp.g5_boundary(f3, dr)
    a5 = [c for c in g5["checks"] if c["check"].startswith("G5a")][0]["result"]
    rec("c3 D phase reversed", "G5a (G1 stays PASS)", {"G1": g1["result"], "G5a": a5}, a5 == src.FAIL and g1["result"] == src.PASS,
        "time reversal is valid physics: not a Hamiltonian failure")
    g2 = mp.g2_mapping(f, dr, row, mapping_cross=False)
    a2 = [c for c in g2["checks"] if c["check"].startswith("G2a")][0]["result"]
    rec("c4 S-D cross term omitted", "G2a", {"G2a": a2}, a2 == src.FAIL)
    g5 = mp.g5_boundary(f, dr, destination={"open_edge": {"kind": "lattice_wall", "reflecting": True}}, claim_isolated_escape=True)
    b5 = [c for c in g5["checks"] if c["check"].startswith("G5b")][0]["result"]
    rec("c5 reflecting wall with isolated-escape claim", "G5b", {"G5b": b5}, b5 == src.FAIL)
    g3 = mp.g3_population(destination=None, dust_declared_overlap=True)
    rec("c6 dust reused without partition", "G3", {"G3": g3["result"]}, g3["result"] == src.FAIL)
    g4 = con.g4_gr_initial_data(f, dr, j_supplied=con.momentum_density(f, dr, include_ps=False))
    a4 = [c for c in g4["checks"] if c["check"].startswith("G4a")][0]["result"]
    rec("c7 PS nonzero, momentum source without -PS S_r", "G4a", {"G4a": a4}, a4 == src.FAIL)
    return out


def main(v7, v5, mc, outdir):
    manifestC = json.load(open(mc))
    results = {"runtime": src.environment_report(), "inputs": {"v7": src.sha256_file(v7), "v5": src.sha256_file(v5),
               "manifestC": src.sha256_file(mc)}, "sources": {}, "controls": {}}
    with tempfile.TemporaryDirectory() as t7, tempfile.TemporaryDirectory() as t5:
        r7 = os.path.join(src.safe_extract(v7, t7), "ringkit_v7")
        r5 = os.path.join(src.safe_extract(v5, t5), "ringkit")
        runs = {"C": (manifestC, v7, r7)}
        for k, (mem, sha, n, dr) in V5_SOURCES.items():
            runs[k] = (v5_manifest(mem, sha, n, dr), v5, r5)
        for k, (m, arch, root) in runs.items():
            rep, f, d, row, dr = run_source(k, m, arch, root)
            if k == "C":
                # Projects-PROPOSED corrected code_ref (file hash of hbfull.py; v6 archive sha is unverified here)
                mc2 = copy.deepcopy(m)
                hb = src.sha256_file(os.path.join(root, "hbfull.py"))
                mc2["source"]["code_ref"] = m["source"]["code_ref"].replace("108cf5d52b9937b7e623323993e261eb58583cb97c43e78fecc0a2ca42f3cd1e", hb)
                g = src.g0_provenance(mc2, arch, root)
                rep["proposed_corrected_manifest_G0"] = {"result": g["result"], "note": "Projects-proposed only; hbfull.py file sha " + hb +
                                                         "; ringkit_v6.tgz sha 108cf5d5... listed by Lumen, not verified by Projects (archive not held)"}
                results["controls"] = controls(m, arch, root, f, row, dr)
                results["controls_ran_on"] = "C"
                if row is not None:
                    rep["ledger_row_t300"] = row
            results["sources"][k] = rep
    pre = {"C": "f3675f71", "A": "4b3155a5", "B": "af0b7653"}
    for k, rep in results["sources"].items():
        od = os.path.join(outdir, "b7_ring_import_" + pre[k])
        os.makedirs(od, exist_ok=True)
        json.dump({"runtime": results["runtime"], "inputs": results["inputs"], **rep,
                   **({"controls": results["controls"]} if k == "C" else {})}, open(os.path.join(od, "report.json"), "w"), indent=1, default=str)
    json.dump(results, open(os.path.join(outdir, "summary.json"), "w"), indent=1, default=str)
    return results


if __name__ == "__main__":
    main(*sys.argv[1:5])
