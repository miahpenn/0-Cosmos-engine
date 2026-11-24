"""Section-order step 1: show that the closure gate and readings P1-P4 can discriminate on synthetic tag sets.

No production number enters here. Exits non-zero if any case gives other than its expected verdict."""
import json

import numpy as np

from research.w1_aa_split import readings

TOTAL = {"1-2": 6.5e-11, "2-4": 6.5e-11/3.5}
C = 1e-16


def pairs(tags, totals=TOTAL):
    return {p: {"tags": {n: v[i] for n, v in tags.items()}, "total": totals[p], "c_split": C}
            for i, p in enumerate(readings.PAIRS)}


def scaled(**fractions):
    """Tags as fractions of each pair's total."""
    return {n: (f*TOTAL["1-2"], f*TOTAL["2-4"]) for n, f in fractions.items()}


def cases():
    flat = {f"t{i}": (TOTAL["1-2"]/8, TOTAL["2-4"]/8) for i in range(8)}
    no_fall = {p: 6.5e-11 for p in readings.PAIRS}
    return {
        "one dominant tag": (pairs(scaled(A=.9, B=.06, C=.04)), set(), ("P1", "A")),
        "dominant tag that does not fall 2x": (
            pairs({"A": (.9*6.5e-11, .9*6.5e-11), "B": (.06*6.5e-11,)*2, "C": (.04*6.5e-11,)*2}, no_fall),
            set(), ("P4", None)),
        "dominant tag under half the total": (pairs(scaled(A=.45, B=-.2, C=.2, D=.2, E=.2, F=.15)), set(), ("P4", None)),
        "origin group, no single tag": (pairs(scaled(o1=.3, o2=.3, o3=.3, b1=.05, b2=.05)), {"o1", "o2", "o3"},
                                        ("P2", "origin group")),
        "cancelling pair": (pairs(scaled(A=40., B=-39., C=.02)), set(), ("P3", None)),
        "flat split": (pairs(flat), set(), ("P4", None)),
    }


def gate_cases():
    rng = np.random.default_rng(588)
    terms = rng.normal(size=(40, 6))*np.logspace(0, 4, 40)[:, None]
    exact = np.sum(terms.astype(np.longdouble), axis=0)
    binary64 = np.zeros(6)
    for t in terms:            # a production-style sequential binary64 sum
        binary64 = binary64 + t
    bound = 10*readings.EPS*np.abs(terms).sum(axis=0)
    return {
        "exact sum": (terms, exact, True),
        "binary64 production sum": (terms, binary64, True),
        "closure fault (production off by 20 eps * sum|terms| at one cell)": (
            terms, exact + np.where(np.arange(6) == 3, 2*bound, 0), False),
        "dropped term": (terms[1:], exact, False),
    }


def main():
    report, failed = {"readings": {}, "closure_gate": {}}, False
    for name, (data, origin, expected) in cases().items():
        got = readings.classify(data, origin)
        ok = (got["reading"], got["pinned"]) == expected
        failed |= not ok
        report["readings"][name] = {"expected": list(expected), "got": [got["reading"], got["pinned"]],
                                    "pass": ok, "absolute_sum_over_total": got["absolute_sum_over_total"],
                                    "largest_opposing_pair": got["largest_opposing_pair"]}
    for name, (terms, production, expected) in gate_cases().items():
        passed, ratio = readings.closure(terms, production)
        ok = passed == expected
        failed |= not ok
        report["closure_gate"][name] = {"expected_pass": expected, "gate_passed": passed,
                                        "worst_gap_over_bound": ratio, "pass": ok}
    report["all_discriminate"] = not failed
    print(json.dumps(report, indent=2, allow_nan=False))
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
