"""Non-destructive field-by-field comparison for Issue #4 reference regeneration.

Writes a JSON report. Never updates committed reference fixtures.
"""
import hashlib
import json
import os
import platform
import sys
from decimal import Decimal
from pathlib import Path

import numpy as np
import mpmath as mp

ROOT = Path(__file__).resolve().parents[3]


def sha(path):
    p = ROOT / path
    return hashlib.sha256(p.read_bytes()).hexdigest() if p.is_file() else None


def dec(value):
    return Decimal(str(value))


def compare_vector(label, committed, regenerated):
    if len(committed) != len(regenerated):
        raise ValueError(f"{label}: vector length mismatch {len(committed)} != {len(regenerated)}")
    deltas = [dec(b) - dec(a) for a, b in zip(committed, regenerated)]
    abs_deltas = [abs(x) for x in deltas]
    return {
        "length": len(deltas),
        "exact_match": all(x == 0 for x in deltas),
        "mismatch_count": sum(x != 0 for x in deltas),
        "max_abs_delta": str(max(abs_deltas, default=Decimal(0))),
        "max_abs_delta_index": (abs_deltas.index(max(abs_deltas)) if abs_deltas else None),
        "signed_deltas_regenerated_minus_committed": [str(x) for x in deltas],
    }


def main():
    root_path = Path(sys.argv[1])
    deriv_path = Path(sys.argv[2])
    report_path = Path(sys.argv[3])
    committed_root_path = ROOT / "tests/data/issue4_affine_root/issue4_affine_root_reference.json"
    committed_deriv_path = ROOT / "tests/data/issue4_implicit_derivative/issue4_implicit_derivative_reference.json"
    a = json.loads(committed_root_path.read_text())
    b = json.loads(root_path.read_text())
    c = json.loads(committed_deriv_path.read_text())
    d = json.loads(deriv_path.read_text())

    root_vector = compare_vector("root", a["root"], b["root"])
    root_cond = {
        "committed": a["cond_inf"],
        "regenerated": b["cond_inf"],
        "exact_match": dec(a["cond_inf"]) == dec(b["cond_inf"]),
        "signed_delta_regenerated_minus_committed": str(dec(b["cond_inf"]) - dec(a["cond_inf"])),
    }
    derivatives = {}
    for epsilon in c["epsilons"]:
        derivatives[epsilon] = compare_vector(
            f"derivative[{epsilon}]", c["derivative"][epsilon], d["derivative"][epsilon]
        )

    # Metadata is reported independently; it is not part of the CI numeric equality predicate.
    root_metadata_keys = sorted(set(a) | set(b)) 
    deriv_metadata_keys = sorted(set(c) | set(d))
    root_metadata = {k: {"committed": a.get(k), "regenerated": b.get(k), "equal": a.get(k) == b.get(k)}
                     for k in root_metadata_keys if k not in ("root", "cond_inf") and a.get(k) != b.get(k)}
    deriv_metadata = {k: {"committed": c.get(k), "regenerated": d.get(k), "equal": c.get(k) == d.get(k)}
                      for k in deriv_metadata_keys if k != "derivative" and c.get(k) != d.get(k)}

    inputs = [
        "tests/data/issue4_affine_root/generate_reference.py",
        "tests/data/issue4_affine_root/issue4_affine_root_reference.json",
        "tests/data/issue4_implicit_derivative/generate_derivative_reference.py",
        "tests/data/issue4_implicit_derivative/issue4_implicit_derivative_reference.json",
        "tests/data/issue4_research/mp_fixture_helper.py",
        "tests/test_issue4_cmc_lapse_tangent.py",
        "tests/test_issue4_research_port_validation.py",
        "tools/issue4_cmc_lapse_tangent.py",
        "tests/data/issue4_research/mp_fixture_port.py",
    ]
    provenance = {
        "python": platform.python_version(),
        "python_executable": sys.executable,
        "numpy": np.__version__,
        "mpmath": mp.__version__,
        "platform": platform.platform(),
        "runner_temp": os.environ.get("RUNNER_TEMP"),
        "repo_commit": os.popen("git rev-parse HEAD").read().strip(),
        "sha256": {p: sha(p) for p in inputs},
        "missing_paths": [p for p in inputs if sha(p) is None],
        "generator_reported_provenance": {
            "committed_root": a.get("provenance"),
            "regenerated_root": b.get("provenance"),
            "committed_derivative": c.get("provenance"),
            "regenerated_derivative": d.get("provenance"),
        },
    }
    result = {
        "schema": "issue4-reference-regeneration-comparison-v1",
        "scope": "synthetic research fixtures only; no production physics claim",
        "numeric_comparison": {
            "affine_root_vector": root_vector,
            "affine_root_cond_inf": root_cond,
            "implicit_derivative_by_epsilon": derivatives,
        },
        "metadata_differences_not_in_numeric_predicate": {
            "affine_root": root_metadata,
            "implicit_derivative": deriv_metadata,
        },
        "provenance": provenance,
    }
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "affine_root_exact": root_vector["exact_match"] and root_cond["exact_match"],
        "derivative_exact": all(v["exact_match"] for v in derivatives.values()),
        "root_max_abs_delta": root_vector["max_abs_delta"],
        "derivative_max_abs_delta": {k: v["max_abs_delta"] for k, v in derivatives.items()},
        "report": str(report_path),
        "provenance_sha256": provenance["sha256"],
        "missing_paths": provenance["missing_paths"],
    }, indent=2))
    if not (root_vector["exact_match"] and root_cond["exact_match"]
            and all(v["exact_match"] for v in derivatives.values())):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
