"""RESEARCH ONLY. B7 ring source intake: fail-closed, read-only (no engine, no vendor, no production archive).

Design: B7 build design (Drive 1MrTDWNnMajevIAzCXfR6vegBvWktsjovqRz1KB3Ucfg) sections 3, 5, 6; Projects' pre-run
entry in the coordination log (2026-11-22). Gates G0 (provenance) and G1 (valid snapshot) live here, with sub-checks.
Structural checks admit no tolerance. The one numeric check, the energy ledger identity (G1d), is fixed at relative
1e-9 before any run (same formula, same float64 arrays). Nothing here admits a source: every manifest says
candidate_not_admitted and this code never sets anything else.
"""
import hashlib
import json
import os
import platform
import sys
import tarfile

import numpy as np

PASS, FAIL, UNTESTED = "PASS", "FAIL", "UNTESTED"
E_LEDGER_RTOL = 1e-9          # fixed before the run (design entry): same formula, float64
FLUX_RADII = (15, 25, 40)
KEYS_T300 = ("r", "t300_s", "t300_ps", "t300_d", "t300_pd")


def sha256_bytes(b):
    return hashlib.sha256(b).hexdigest()


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for blk in iter(lambda: f.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def runtime():
    import scipy
    return {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__}


def check(name, ok, detail=""):
    return {"check": name, "result": PASS if ok else FAIL, "detail": detail}


def gate_result(checks):
    if any(c["result"] == FAIL for c in checks):
        return FAIL
    if any(c["result"] == UNTESTED for c in checks):
        return UNTESTED
    return PASS


def safe_member_path(root, name):
    """H5: basename-only member name resolved strictly inside root (symlinks included); None if unsafe."""
    if not isinstance(name, str) or not name or "\x00" in name or "/" in name or "\\" in name:
        return None
    if name in (".", "..") or name.startswith(".") or os.path.isabs(name):
        return None
    rr = os.path.realpath(root)
    p = os.path.realpath(os.path.join(root, name))
    if os.path.dirname(p) != rr:
        return None
    return os.path.join(root, name)


def safe_relative_path(root, name):
    """H5: MANIFEST entries may sit in sub-folders but never escape: no absolute path, no '.', '..' or empty component,
    no backslash/NUL, realpath (symlinks resolved) inside realpath(root). None if unsafe."""
    if not isinstance(name, str) or not name or "\x00" in name or "\\" in name or os.path.isabs(name):
        return None
    if any(part in ("", ".", "..") for part in name.split("/")):
        return None
    rr = os.path.realpath(root)
    p = os.path.realpath(os.path.join(root, name))
    if os.path.commonpath([rr, p]) != rr or p == rr:
        return None
    return os.path.join(root, name)


def manifest_digests(root):
    """name -> sha256 from root/MANIFEST.sha256 (names as written; no path is opened here)."""
    d = {}
    with open(os.path.join(root, "MANIFEST.sha256")) as f:
        for line in f:
            line = line.strip()
            if line:
                sha, name = line.split(None, 1)
                d[name.strip().lstrip("*")] = sha
    return d


def verify_internal_manifest(root):
    """Every line of root/MANIFEST.sha256 ('<sha>  <name>') against the file in root. Returns (ok, n, bad)."""
    path = os.path.join(root, "MANIFEST.sha256")
    bad, n = [], 0
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            sha, name = line.split(None, 1)
            name = name.strip().lstrip("*")
            n += 1
            p = safe_relative_path(root, name)
            if p is None or not os.path.isfile(p) or sha256_file(p) != sha:
                bad.append(name)
    return (not bad), n, bad


def safe_extract(tgz, dest):
    """Extract into an empty directory; refuse absolute paths, '..' and links."""
    os.makedirs(dest, exist_ok=True)
    with tarfile.open(tgz, "r:gz") as t:
        for m in t.getmembers():
            n = os.path.normpath(m.name)
            if n.startswith("..") or os.path.isabs(n) or m.issym() or m.islnk():
                raise ValueError("unsafe tar member: " + m.name)
        t.extractall(dest, filter="data")
    return dest


def g0_provenance(manifest, archive_path, member_root, env_required=True, schema_path=None, expected_code_files=None):
    """G0. manifest = parsed source manifest (b7-ring-source-v1); member_root = extracted folder holding the members."""
    src = manifest["source"]
    cs = []
    got = sha256_file(archive_path)
    cs.append(check("G0a archive sha256 equals declared", got == src["archive_sha256"], got))
    ok, n, bad = verify_internal_manifest(member_root)
    cs.append(check("G0b internal MANIFEST.sha256 verifies", ok, "%d entries, bad: %s" % (n, bad)))
    dp = safe_member_path(member_root, src.get("data_member"))
    dgot = sha256_file(dp) if (dp and os.path.isfile(dp)) else ("unsafe path" if dp is None else "missing")
    cs.append(check("G0c data member path safe, sha256 equals declared", dp is not None and dgot == src["data_sha256"], dgot))
    try:
        listed = manifest_digests(member_root).get(src.get("data_member"))
    except OSError:
        listed = None
    cs.append(check("G0h data member listed in internal MANIFEST.sha256 with the declared digest",
                    dp is not None and listed is not None and listed == src["data_sha256"], str(listed)))
    keys = manifest["slice"]["field_keys"]
    want = [keys[k] for k in ("r", "s", "ps", "d", "pd")]
    try:
        have = list(np.load(dp).keys()) if (dp and os.path.isfile(dp)) else []
    except Exception:  # noqa: BLE001 - unreadable data member is a FAIL
        have = []
    miss = [k for k in want if k not in have]
    cs.append(check("G0d declared field keys present", not miss, "missing: %s" % miss))
    # G0e (H2): code_ref must cite a non-empty, exact, duplicate-free list of '<file>.py ... sha256=<64 hex>' refs. Every
    # 'sha256=' in code_ref must parse, each file must exist in the member folder with a safe name and match in full, and
    # every expected file (if given) must be cited. Zero refs is FAIL.
    import re
    cr = src.get("code_ref", "") or ""
    refs = re.findall(r"([^\s;()]+\.py)\b[^;]*?sha256=([0-9a-fA-F]+)", cr)
    n_sha = len(re.findall(r"sha256=", cr))
    probs = []
    if not refs:
        probs.append("no code refs parsed")
    if len(refs) != n_sha:
        probs.append("%d 'sha256=' but %d parsed refs" % (n_sha, len(refs)))
    names = [x[0] for x in refs]
    if len(set(names)) != len(names):
        probs.append("duplicate refs")
    for fname, sha in refs:
        if safe_member_path(member_root, fname) is None:
            probs.append("unsafe member name " + fname)
            continue
        if len(sha) != 64:
            probs.append("%s: digest is %d hex, not 64" % (fname, len(sha)))
            continue
        p = os.path.join(member_root, fname)
        if not os.path.isfile(p):
            probs.append(fname + ": member file missing")
        elif sha256_file(p) != sha.lower():
            probs.append("%s: cited %s... actual %s..." % (fname, sha[:8], sha256_file(p)[:8]))
    for fname in (expected_code_files or ()):
        if fname not in names:
            probs.append("expected file not cited: " + fname)
    cs.append(check("G0e code_ref: non-empty exact file hashes, all members match", not probs, "; ".join(probs) or "%d refs" % len(refs)))
    env = src.get("environment", "unknown")
    env_ok = isinstance(env, dict) and all(v not in (None, "", "unknown") for v in env.values())
    cs.append(check("G0f source environment recorded", env_ok or not env_required, str(env)))
    # G0g (H3): the manifest must validate against the b7-ring-source-v1 schema (Draft 2020-12). No schema = UNTESTED.
    if schema_path is None:
        cs.append({"check": "G0g manifest validates against b7-ring-source-v1 schema", "result": UNTESTED, "detail": "no schema supplied"})
    else:
        try:
            import jsonschema
            jsonschema.Draft202012Validator(json.load(open(schema_path))).validate(manifest)
            cs.append(check("G0g manifest validates against b7-ring-source-v1 schema", True, "schema sha256 " + sha256_file(schema_path)[:16]))
        except Exception as e:  # noqa: BLE001 - any failure to validate is a FAIL
            cs.append(check("G0g manifest validates against b7-ring-source-v1 schema", False, str(e).splitlines()[0][:200]))
    return {"gate": "G0", "result": gate_result(cs), "checks": cs}


def load_slice(npz_path, keys):
    d = np.load(npz_path)
    return {k: np.asarray(d[v]) for k, v in keys.items()}, d


def energy_ledger(s, d, ps, pd, dr, js=None):
    """E(<R) exactly as in ringkit evo.py energy(): cell energy + face gradient energy, J = round(R/dr) cells."""
    n = len(s)
    r = (np.arange(n) + 0.5) * dr
    rf = np.arange(n + 1) * dr
    vol = 4 * np.pi * r ** 2 * dr
    V = 0.5 * s * s - 0.5 * d * d + (s ** 4 + 6 * s * s * d * d + d ** 4) / 8 + 0.5
    ce = (0.5 * ps ** 2 + 0.5 * pd ** 2 + V) * vol
    gs, gd = np.diff(s) / dr, np.diff(d) / dr
    ge = 0.5 * (gs ** 2 + gd ** 2) * 4 * np.pi * rf[1:n] ** 2 * dr
    out = []
    for R in (js or FLUX_RADII):
        j = int(round(R / dr))
        out.append(float(ce[:j].sum() + ge[:j - 1].sum()))
    return out


def g1_snapshot(f, manifest, ledger_row=None):
    """G1 on a slice dict f (r, s, ps, d, pd). ledger_row = stored (E<15, E<25, E<40) at the snapshot time or None."""
    cs = []
    arrs = [f[k] for k in ("r", "s", "ps", "d", "pd")]
    cs.append(check("G1a float64, equal length, finite",
                    all(a.dtype == np.float64 for a in arrs) and len({a.size for a in arrs}) == 1
                    and all(np.isfinite(a).all() for a in arrs), "n=%d" % arrs[0].size))
    r = f["r"]
    rg = manifest["slice"]["radial_grid"]
    dr0 = rg["dr"]
    dd = np.diff(r)
    uni = bool(dd.size and (dd > 0).all() and np.max(np.abs(dd - dr0)) <= 1e-12)
    cs.append(check("G1b r uniform, increasing, cell-centred, declared grid",
                    uni and abs(r[0] - dr0 / 2) <= 1e-12 and r.size == rg["n_cells"]
                    and abs(r[-1] + dr0 / 2 - rg["r_max"]) <= 1e-9 if r.size else False,
                    "n=%d dr0=%g r0=%g" % (r.size, dr0, r[0] if r.size else float("nan"))))
    nz = bool(np.any(f["ps"] != 0) and np.any(f["pd"] != 0))
    cs.append(check("G1c ps and pd nonzero (bare profile rejected)", nz))
    if ledger_row is None:
        cs.append({"check": "G1d energy ledger identity", "result": UNTESTED, "detail": "no ledger row at this time"})
    else:
        e = energy_ledger(f["s"], f["d"], f["ps"], f["pd"], dr0)
        rel = [abs(a - b) / abs(b) for a, b in zip(e, ledger_row)]
        cs.append(check("G1d energy E(<15,<25,<40) equals stored ledger (rel %.0e)" % E_LEDGER_RTOL,
                        max(rel) <= E_LEDGER_RTOL, "E=%s ledger=%s maxrel=%.2e" % (
                            ["%.9g" % x for x in e], ["%.9g" % x for x in ledger_row], max(rel))))
    return {"gate": "G1", "result": gate_result(cs), "checks": cs}


def ledger_row_at(rows, t, lo=8):
    """rows layout (ringkit evo.py): col 0 = t, cols 8..10 = E(<15,<25,<40). None if t is not sampled."""
    rows = np.asarray(rows)
    i = np.where(np.abs(rows[:, 0] - t) < 1e-6)[0]
    return None if i.size == 0 or rows.shape[1] < lo + 3 else [float(x) for x in rows[i[0], lo:lo + 3]]


def fingerprint(f):
    h = hashlib.sha256()
    for k in ("r", "s", "ps", "d", "pd"):
        h.update(np.ascontiguousarray(f[k], dtype="<f8").tobytes())
    return h.hexdigest()


def load_manifest(path, schema_path=None):
    m = json.load(open(path))
    if schema_path:
        import jsonschema
        jsonschema.Draft202012Validator(json.load(open(schema_path))).validate(m)
    assert m["schema_version"] == "b7-ring-source-v1" and m["claim"] == "candidate_not_admitted"
    return m


def environment_report():
    r = runtime()
    r["platform"] = platform.platform()
    r["argv0"] = os.path.basename(sys.argv[0])
    return r
