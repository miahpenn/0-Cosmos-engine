"""Detection-floor study for the implicit-derivative criteria (research record, not a CI test).

Injects known errors into the solver's implicit tangent and records which predeclared criterion
rejects each. Thresholds are the committed predeclared values and are NOT adjusted here.
Run from repository root: python tests/data/issue4_implicit_derivative/detection_floor_study.py <out.json>
"""
import importlib.util
import json
import pathlib
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO))
_spec = importlib.util.spec_from_file_location("impl_test", REPO / "tests" / "test_issue4_implicit_derivative_contract.py")
D = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(D)

LEVELS = [1e-9, 1e-8, 1e-7, 2e-7, 5e-7, 1e-6, 2e-6, 5e-6, 1e-5, 2e-5, 5e-5, 1e-4]
ABS_LEVELS = [1e-12, 1e-11, 1e-10, 2e-10, 5e-10, 1e-9, 2e-9, 5e-9, 1e-8, 2e-8, 5e-8, 1e-7]


def decision(absolute, normwise_rel):
    a_fail = absolute > D.ABS_THRESHOLD
    r_fail = normwise_rel > D.REL_THRESHOLD
    if a_fail and r_fail:
        return "both"
    if a_fail:
        return "absolute"
    if r_fail:
        return "relative"
    return "none"


def main(out):
    ref, data = D._reference()
    base = np.asarray(D._tangent_and_args().alpha_dot, dtype=float)
    fd = ref["2e-05"]
    rows = {"uniform_relative": [], "single_component_absolute": []}
    for d in LEVELS:
        inj = base * (1.0 + d)
        a, r = D._criteria(inj, fd)
        rows["uniform_relative"].append(dict(level=d, absolute=a, normwise_rel=r, rejected_by=decision(a, r)))
    for e in ABS_LEVELS:
        inj = base.copy(); inj[0] += e
        a, r = D._criteria(inj, fd)
        rows["single_component_absolute"].append(dict(level=e, absolute=a, normwise_rel=r, rejected_by=decision(a, r)))

    def boundary(kind):
        rej = [x for x in rows[kind] if x["rejected_by"] != "none"]
        return min(rej, key=lambda x: x["level"]) if rej else None

    summary = dict(
        thresholds=dict(absolute=D.ABS_THRESHOLD, normwise_relative=D.REL_THRESHOLD),
        max_abs_tangent=float(np.max(np.abs(base))),
        boundary_uniform_relative=boundary("uniform_relative"),
        boundary_single_component_absolute=boundary("single_component_absolute"),
        note="Thresholds are the committed predeclared values; they were not adjusted after observing this study.",
        rows=rows, provenance=dict(mpmath=data["provenance"]["mpmath"], reference_dps=data["dps"]))
    json.dump(summary, open(out, "w"), indent=1)
    print(json.dumps({k: summary[k] for k in ("boundary_uniform_relative", "boundary_single_component_absolute")}, indent=1))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else str(HERE / "detection_floor_study.json"))
