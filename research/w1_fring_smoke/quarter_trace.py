"""Compare three timestep sizes on the fine W1 ring grid."""
import hashlib
import json
import math
from pathlib import Path

from research.w1_fring_smoke import areal_trace, origin_trace

PREVIOUS = Path("results/research/w1_fring_areal_trace")
END_TIME = 0.0075
RESIDUALS = ("H_production", "H_metric_7", "H_metric_9")


def main():
    files = [PREVIOUS / f"A-0.025{suffix}.json" for suffix in ("", "-half")]
    prior = [json.loads(path.read_text()) for path in files]
    quarter_rows = []

    def record(state, vacuum, row):
        quarter_rows.append(areal_trace.observe(state, vacuum, row))

    replay = origin_trace.run(0.025, "A", 4, end_time=END_TIME, accepted_observer=record)
    for div, data in enumerate(prior, start=1):
        div = 2**(div-1)
        meta = data["metadata"]
        for key in ("source", "vendor", "input_sha256", "arm", "dr", "sample_cells", "sample_r"):
            if meta[key] != replay[key]:
                raise ValueError(f"mismatched baseline metadata: {key}")
        if meta["dt"] != replay["dt"] * (4 // div) or data["rows"][0] != quarter_rows[0]:
            raise ValueError("mismatched timestep or initializer")

    nominal, half = (data["rows"] for data in prior)
    comparisons = []
    for i in range(1, len(quarter_rows) // 4 + 1):
        rows = (nominal[i], half[2*i], quarter_rows[4*i])
        time = rows[0]["t"]
        if any(not math.isclose(row["t"], time, abs_tol=1e-17, rel_tol=0) for row in rows):
            raise ValueError("time mismatch")
        residuals = {}
        for key in RESIDUALS:
            coarse, medium, fine = (row["cells"][key][0] for row in rows)
            delta_coarse = coarse - medium
            delta_fine = medium - fine
            residuals[key] = {
                "H": [coarse, medium, fine],
                "deltas": [delta_coarse, delta_fine],
                "ratio": delta_coarse / delta_fine if delta_fine else None,
            }
        comparisons.append({"t": time, "cell": 0, "residuals": residuals})
    metadata = {key: value for key, value in replay.items() if key not in ("rows", "first_step_stages")}
    paths = [Path(__file__), Path(areal_trace.__file__), *files]
    metadata["sha256"] = {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}
    print(json.dumps({"metadata": metadata, "quarter_rows": quarter_rows,
                      "quarter_first_step": replay["first_step_stages"],
                      "matched_centre": comparisons}, allow_nan=False))


if __name__ == "__main__":
    main()
