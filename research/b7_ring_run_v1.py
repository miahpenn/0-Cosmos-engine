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


def blocked(g, why):
    return {"gate": g, "result": src.UNTESTED, "checks": [{"check": g + " blocked", "result": src.UNTESTED, "detail": why}]}


def run_source(label, manifest, archive, root, schema=None, expected=None):
    keys = manifest["slice"]["field_keys"]
    rep = {"source": label, "data_member": manifest["source"]["data_member"], "data_sha256": manifest["source"]["data_sha256"]}
    g0 = src.g0_provenance(manifest, archive, root, schema_path=schema, expected_code_files=expected)
    dp = src.safe_member_path(root, manifest["source"]["data_member"])
    if dp is None or not os.path.isfile(dp):
        rep["gates"] = {"G0": g0, **{g: blocked(g, "blocked: data member path unsafe or missing; nothing loaded") for g in
                                     ("G1", "G2", "G3", "G4", "G5", "G6", "G7")}}
        rep["adapter"] = mp.adapter_status(None, rep["gates"])
        return rep, None, None, None, None
    f, d = src.load_slice(dp, keys)
    dr = manifest["slice"]["radial_grid"]["dr"]
    row = src.ledger_row_at(d["rows"], manifest["slice"]["time"], lo=LEDGER_COL[label])
    g1 = src.g1_snapshot(f, manifest, row)
    gates = {"G0": g0, "G1": g1}
    if g0["result"] == src.FAIL:
        rep["forensic_only"] = "G0 failed: the G1 numbers are forensic diagnostics, not evidence for mapping or admission"
    if g0["result"] == src.FAIL or g1["result"] == src.FAIL:
        why = "blocked: G0/G1 failed; no mapping computed"
        for g in ("G2", "G3", "G4", "G5"):
            gates[g] = blocked(g, why)
    else:
        gates.update({"G2": mp.g2_mapping(f, dr, row), "G3": mp.g3_population(), "G4": con.g4_gr_initial_data(f, dr),
                      "G5": mp.g5_boundary(f, dr)})
    gates["G6"] = {"gate": "G6", "result": src.UNTESTED, "checks": [{"check": "same-source dressed BVP at measured omega",
                   "result": src.UNTESTED, "detail": "not implemented in this harness; the v5 dressed2 energy comparison is superseded "
                   "by v6/v7; nothing inferred"}]}
    gates["G7"] = {"gate": "G7", "result": src.UNTESTED, "checks": [{"check": "dynamics", "result": src.UNTESTED, "detail": "not authorised"}]}
    rep["gates"] = gates
    rep["adapter"] = mp.adapter_status(None, gates)
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


def git_ids():
    import subprocess
    try:
        c = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
        t = subprocess.run(["git", "rev-parse", "HEAD^{tree}"], capture_output=True, text=True).stdout.strip()
        return {"commit": c, "tree": t}
    except OSError:
        return {"commit": "unknown", "tree": "unknown"}


SCHEMA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "b7_ring_schemas", "b7_ring_source_v1.schema.h5.json")
C_EXPECTED = ("dressed_run.py", "evo.py", "hbfull.py")


def main(v7, v5, mc, outdir, tag=""):
    manifestC = json.load(open(mc))
    results = {"runtime": src.environment_report(), "harness": git_ids(), "tag": tag,
               "inputs": {"v7": src.sha256_file(v7), "v5": src.sha256_file(v5), "manifestC": src.sha256_file(mc),
                          "schema": src.sha256_file(SCHEMA), "schema_original": src.sha256_file(os.path.join(os.path.dirname(SCHEMA), "b7_ring_source_v1.schema.json"))}, "sources": {}, "controls": {}}
    with tempfile.TemporaryDirectory() as t7, tempfile.TemporaryDirectory() as t5:
        r7 = os.path.join(src.safe_extract(v7, t7), "ringkit_v7")
        r5 = os.path.join(src.safe_extract(v5, t5), "ringkit")
        runs = {"C": (manifestC, v7, r7, C_EXPECTED)}
        for k, (mem, sha, n, dr) in V5_SOURCES.items():
            runs[k] = (v5_manifest(mem, sha, n, dr), v5, r5, None)
        for k, (m, arch, root, exp) in runs.items():
            rep, f, d, row, dr = run_source(k, m, arch, root, SCHEMA, exp)
            if k == "C" and f is not None:
                results["controls"] = controls(m, arch, root, f, row, dr)
                results["controls_ran_on"] = "C" + ("" if rep["gates"]["G0"]["result"] == src.PASS else " (forensic: G0 failed)")
            results["sources"][k] = rep
    pre = {"C": "f3675f71", "A": "4b3155a5", "B": "af0b7653"}
    sfx = ("_" + tag) if tag else ""
    for k, rep in results["sources"].items():
        od = os.path.join(outdir, "b7_ring_import_" + pre[k] + sfx)
        os.makedirs(od, exist_ok=True)
        json.dump({"runtime": results["runtime"], "harness": results["harness"], "inputs": results["inputs"], **rep,
                   **({"controls": results["controls"], "controls_ran_on": results["controls_ran_on"]} if k == "C" else {})},
                  open(os.path.join(od, "report.json"), "w"), indent=1, default=str)
    json.dump(results, open(os.path.join(outdir, "summary%s.json" % sfx), "w"), indent=1, default=str)
    return results


if __name__ == "__main__":
    main(*sys.argv[1:6])
