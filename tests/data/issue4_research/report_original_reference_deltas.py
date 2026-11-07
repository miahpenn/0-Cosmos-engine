"""Emit Decimal deltas for the original Issue #4 reference comparison outputs.

This does not regenerate or alter fixtures. The workflow invokes the original
reference generators first; this script only compares their temporary outputs
against the committed JSON using the same fields as the original predicate.
"""
import json
import os
from decimal import Decimal
from pathlib import Path

root = Path(os.environ["RUNNER_TEMP"])
committed_root = json.loads(Path("tests/data/issue4_affine_root/issue4_affine_root_reference.json").read_text())
generated_root = json.loads((root / "root.json").read_text())
committed_deriv = json.loads(Path("tests/data/issue4_implicit_derivative/issue4_implicit_derivative_reference.json").read_text())
generated_deriv = json.loads((root / "deriv.json").read_text())

def compare(label, expected, actual):
    deltas = [Decimal(str(a)) - Decimal(str(e)) for e, a in zip(expected, actual)]
    mismatches = [i for i, d in enumerate(deltas) if d != 0]
    max_abs = max((abs(d) for d in deltas), default=Decimal(0))
    index = max(range(len(deltas)), key=lambda i: abs(deltas[i])) if deltas else None
    return {
        "label": label, "length_expected": len(expected), "length_actual": len(actual),
        "mismatch_count": len(mismatches), "mismatch_indices": mismatches,
        "max_abs_delta": str(max_abs), "max_abs_delta_index": index,
        "signed_deltas_actual_minus_committed": [str(d) for d in deltas],
        "lengths_equal": len(expected) == len(actual),
    }

report = {
    "scope": "original generator entry points and original workflow numeric fields; diagnostic instrumentation only",
    "root": compare("affine_root.root", committed_root["root"], generated_root["root"]),
    "cond_inf": {
        "committed": str(committed_root["cond_inf"]),
        "generated": str(generated_root["cond_inf"]),
        "signed_delta_generated_minus_committed": str(Decimal(str(generated_root["cond_inf"])) - Decimal(str(committed_root["cond_inf"]))),
        "exact_match": Decimal(str(generated_root["cond_inf"])) == Decimal(str(committed_root["cond_inf"])),
    },
    "derivatives": {},
}
for eps in committed_deriv["epsilons"]:
    report["derivatives"][eps] = compare(
        f"implicit_derivative[{eps}]",
        committed_deriv["derivative"][eps],
        generated_deriv["derivative"][eps],
    )
report["original_predicate"] = {
    "root_equal": committed_root["root"] == generated_root["root"],
    "cond_inf_equal": committed_root["cond_inf"] == generated_root["cond_inf"],
    "derivative_equal": all(committed_deriv["derivative"][e] == generated_deriv["derivative"][e]
                            for e in committed_deriv["epsilons"]),
}
(root / "original-path-delta-report.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report, indent=2))
