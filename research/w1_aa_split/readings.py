"""Closure gate and predeclared readings P1-P4 for the W1 A_a rate term split.

DESIGN-W1-A_a-rate-term-split-predeclared-Nov24-2026 (Drive 1tt1m0E8bkTvUGbsGsWG06DIiK2w3a5YEDEP-9XQRpTI), sections 2 and 4.
These functions are fixed before any split number is seen; synthetic.py exercises them.
"""
import numpy as np

EPS = np.finfo(float).eps  # binary64 machine epsilon, the "eps" of the closure gate
PAIRS = ("1-2", "2-4")     # nominal-half, half-quarter


def closure(terms, production):
    """Section 2: |sum of tagged terms - production rate| <= 10 eps * sum |tagged terms|, per cell.

    ``terms`` is a sequence of arrays (one per tag), ``production`` the production A_a rate on the same cells.
    Returns (passed, worst ratio of the left side to the right side)."""
    stacked = np.asarray(terms, dtype=np.longdouble)
    gap = np.abs(stacked.sum(axis=0) - np.asarray(production, dtype=np.longdouble))
    bound = 10*EPS*np.abs(stacked).sum(axis=0)
    passed = bool(np.all(gap <= bound))
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(bound > 0, gap/bound, np.where(gap > 0, np.inf, 0))
    return passed, float(np.max(ratio))


def resolution(total):
    """Readback resolution: one binary64 spacing of the signed total as read back."""
    return float(np.spacing(abs(float(total))))


def _conditions(value, others, data):
    """The three P1 conditions for one contribution (a tag, or the origin sum) against ``others``."""
    half = {p: data[p]["total"] != 0 and value[p]/data[p]["total"] >= .5 for p in PAIRS}
    dominant = {p: abs(value[p]) >= 2*max((abs(o[p]) for o in others), default=0.) for p in PAIRS}
    fall = value["2-4"] != 0 and value["1-2"]/value["2-4"] >= 2
    return {"half_of_total": half, "twice_next_largest": dominant, "falls_2x": bool(fall),
            "met": all(half.values()) and all(dominant.values()) and bool(fall)}


def classify(data, origin):
    """Section 4. ``data[pair]`` = {"tags": {name: signed contribution}, "total": signed tracefree-rate
    difference, "c_split": closure term}; ``origin`` = set of origin-tagged names.
    Precedence is the design's order P1, P2, P3, P4; every condition is reported."""
    names = list(data[PAIRS[0]]["tags"])
    value = {n: {p: data[p]["tags"][n] for p in PAIRS} for n in names}
    out = {"floor": {}, "unresolved": {}, "absolute_sum_over_total": {}, "largest_opposing_pair": {}}
    for p in PAIRS:
        floor = 10*max(abs(data[p]["c_split"]), resolution(data[p]["total"]))
        out["floor"][p] = floor
        out["unresolved"][p] = [n for n in names if abs(value[n][p]) < floor]
        absolute = sum(abs(value[n][p]) for n in names)
        out["absolute_sum_over_total"][p] = absolute/abs(data[p]["total"]) if data[p]["total"] else float("inf")
        positive = max((n for n in names if value[n][p] > 0), key=lambda n: value[n][p], default=None)
        negative = min((n for n in names if value[n][p] < 0), key=lambda n: value[n][p], default=None)
        out["largest_opposing_pair"][p] = [positive, negative]

    p1 = {n: _conditions(value[n], [value[m] for m in names if m != n], data) for n in names}
    out["P1_by_tag"] = p1
    pinned = [n for n in names if p1[n]["met"]]
    group = {p: sum(value[n][p] for n in names if n in origin) for p in PAIRS}
    p2 = _conditions(group, [value[n] for n in names if n not in origin], data)
    out["origin_sum"] = group
    out["P2"] = p2
    p3 = any(out["absolute_sum_over_total"][p] > 10 for p in PAIRS)
    out["P3_condition"] = p3
    if pinned:
        out["reading"], out["pinned"] = "P1", pinned[0]
    elif p2["met"]:
        out["reading"], out["pinned"] = "P2", "origin group"
    elif p3:
        out["reading"], out["pinned"] = "P3", None
    else:
        out["reading"], out["pinned"] = "P4", None
    return out
