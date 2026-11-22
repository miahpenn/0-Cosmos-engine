"""RESEARCH ONLY. B7 ring mapping: a pure diagnostic proposal, never a production initial_structure option.

Gates G2 (mapping), G3 (population) and G5 (boundary) live here. A mapping needs a separately sourced destination
operating point (phi, Pi, units chart / lambda, P1 partition, open edge). With none, adapter_ready is False, the
missing inputs are listed, and G2b/G3/G5 are UNTESTED: that is the correct outcome, not a failure of the build.
Only what can be decided from the source slice itself is decided: G2a (stress-source identity) and G5a (outgoing sign).
Fixed before the run (design entry 2026-11-22). Nothing in engine/ or vendor/ is read or changed.
"""
import numpy as np

from research.b7_ring_source_v1 import (FAIL, PASS, UNTESTED, E_LEDGER_RTOL, check, gate_result, energy_ledger)

REQUIRED_DESTINATION = ("phi", "Pi", "units_chart_lambda", "p1_partition", "open_edge", "gr_cmc_solver_for_ringing_slice")


def potential_toy(s, d, cross=True):
    """Toy S-D potential of the source (lambda = 1, static sea): cross = the 6 s^2 d^2/8 coupling term."""
    return 0.5 * s * s - 0.5 * d * d + (s ** 4 + (6 * s * s * d * d if cross else 0.0) + d ** 4) / 8 + 0.5


def energy_with_potential(f, dr, cross=True, js=(15, 25, 40)):
    s, d, ps, pd = f["s"], f["d"], f["ps"], f["pd"]
    n = s.size
    r = (np.arange(n) + 0.5) * dr
    rf = np.arange(n + 1) * dr
    vol = 4 * np.pi * r ** 2 * dr
    ce = (0.5 * ps ** 2 + 0.5 * pd ** 2 + potential_toy(s, d, cross)) * vol
    gs, gd = np.diff(s) / dr, np.diff(d) / dr
    ge = 0.5 * (gs ** 2 + gd ** 2) * 4 * np.pi * rf[1:n] ** 2 * dr
    out = []
    for R in js:
        j = int(round(R / dr))
        out.append(float(ce[:j].sum() + ge[:j - 1].sum()))
    return out


def g2_mapping(f, dr, ledger_row, destination=None, mapping_cross=True):
    cs = []
    if ledger_row is None:
        cs.append({"check": "G2a stress-source identity", "result": UNTESTED, "detail": "no ledger row"})
    else:
        e = energy_with_potential(f, dr, cross=mapping_cross)
        rel = max(abs(a - b) / abs(b) for a, b in zip(e, ledger_row))
        cs.append(check("G2a mapping potential reproduces ledger energy (rel %.0e)" % E_LEDGER_RTOL, rel <= E_LEDGER_RTOL,
                        "maxrel=%.2e cross_term=%s" % (rel, mapping_cross)))
    # G2b: field/momentum conversion is only decidable in a declared chart. Flat chart (alpha 1, beta 0) round trip is
    # machinery, not physics, so the real sources stay UNTESTED without a destination.
    if destination is None:
        cs.append({"check": "G2b units chart / lambda / momentum conversion", "result": UNTESTED,
                   "detail": "no destination operating point"})
    else:
        cs.append({"check": "G2b units chart / lambda / momentum conversion", "result": UNTESTED,
                   "detail": "destination given but chart conversion not implemented in this harness"})
    return {"gate": "G2", "result": gate_result(cs), "checks": cs}


def g3_population(destination=None, dust_declared_overlap=False):
    if dust_declared_overlap and not (destination or {}).get("p1_partition"):
        cs = [check("G3 no double count: overlapping dust declared without a signed P1 partition", False)]
    elif not (destination or {}).get("p1_partition"):
        cs = [{"check": "G3 signed P1 partition", "result": UNTESTED, "detail": "absent"}]
    else:
        cs = [{"check": "G3 signed P1 partition", "result": UNTESTED, "detail": "present but accounting not implemented"}]
    return {"gate": "G3", "result": gate_result(cs), "checks": cs}


def outgoing_sign(f, dr, lo=15.0, hi=40.0):
    """r^2-weighted mean of -(pd d_r) over lo<=r<=hi: positive = D energy flows outward (one snapshot, ~14 wavelengths)."""
    r, d, pd = f["r"], f["d"], f["pd"]
    dd = np.gradient(d, dr)
    m = (r >= lo) & (r <= hi)
    w = r[m] ** 2
    return float(np.sum(-pd[m] * dd[m] * w) / np.sum(w))


def g5_boundary(f, dr, destination=None, claim_isolated_escape=False):
    cs = []
    v = outgoing_sign(f, dr)
    cs.append(check("G5a outgoing D flux sign over r 15..40 is positive", v > 0, "mean(-pd d_r) r^2-weighted = %.3e" % v))
    edge = (destination or {}).get("open_edge")
    if edge is None:
        cs.append({"check": "G5b destination edge non-reflecting for an isolated-escape claim", "result": UNTESTED,
                   "detail": "no destination edge"})
    else:
        ok = not (edge.get("reflecting", True) and claim_isolated_escape)
        cs.append(check("G5b destination edge non-reflecting for an isolated-escape claim", ok, str(edge)))
    cs.append({"check": "G5c worldtube signed closure on destination", "result": UNTESTED, "detail": "no destination"})
    return {"gate": "G5", "result": gate_result(cs), "checks": cs}


def adapter_status(destination=None):
    have = set((destination or {}).keys())
    miss = [k for k in REQUIRED_DESTINATION if k not in have]
    return {"adapter_ready": not miss, "missing_inputs": miss}
