"""Telescope the direct Aa/K write errors and momentum over pinned W1 runs."""
import argparse
from decimal import Decimal, localcontext
import json
from pathlib import Path
from unittest.mock import patch

import numpy as np

from research.w1_fring_smoke import stage_momentum as replay
from research.w1_fring_smoke.live_lapse_audit import signs
from research.w1_fring_smoke.momentum_precision import derivative, wide

FIELDS = ("Aa", "K")
PARTS = ("computed", "increment_arithmetic", "write")
CELLS = replay.LOCAL_CELLS


def decimals(values):
    return list(map(wide, values[:CELLS]))


class MomentumHistory(replay.Observe):
    def __init__(self, vacuum):
        super().__init__(vacuum)
        self.rows = []
        self.totals = {field: {part: [Decimal(0)] * CELLS for part in PARTS}
                       for field in FIELDS}

    def record(self, name, grid, geometry, scalars, matter):
        if name in ("initial_lapse", "accepted_start", "accepted_end"):
            super().record(name, grid, geometry, scalars, matter)
        if name == "accepted_start" and not self.rows:
            self.initial = {field: decimals(getattr(geometry, field)) for field in FIELDS}
            self.initial_c = self.shear_coefficient(grid, geometry)
            self.initial_rows = self.stages[name]
        elif name == "corrector_before_boundary":
            self.pre_boundary = {field: getattr(geometry, field)[:CELLS].copy()
                                 for field in FIELDS}
        elif name == "accepted_end":
            self.append_row(grid, geometry)

    def shear_coefficient(self, grid, geometry):
        dx, db = (self.vacuum._d(grid, field, 1)[:CELLS] for field in (geometry.X, geometry.b))
        return 3 / grid.centers[:CELLS] - 3 * dx / geometry.X[:CELLS] + 1.5 * db / geometry.b[:CELLS]

    def _geometry_stage_terms(self, grid, geometry, scalars, matter, electromagnetic):
        idx = self.geometry_calls
        terms = super()._geometry_stage_terms(grid, geometry, scalars, matter, electromagnetic)
        if idx in (0, 2):
            self.rhs[idx] = {block: {key: terms[block][key][:CELLS].copy() for key in FIELDS}
                             for block in ("primary_l2", "primary_l3")}
        return terms

    def step(self, state, dt):
        self.rhs = {}
        final_l2 = []
        original = self.vacuum.primary_l2_rhs

        def capture(grid, geometry):
            terms = original(grid, geometry)
            if 2 in self.rhs and self.boundary_calls == 1:
                final_l2.append({key: terms[key][:CELLS].copy() for key in FIELDS})
            return terms

        with patch.object(self.vacuum, "primary_l2_rhs", capture):
            result = super().step(state, dt)
        if len(final_l2) != 1:
            raise ValueError("expected one corrector L2 evaluation")
        self.t = result.t
        self.latest = {}
        with localcontext() as ctx:
            ctx.prec = 60
            for key in FIELDS:
                a = self.rhs[0]["primary_l2"][key]
                b = final_l2[0][key]
                c = self.rhs[0]["primary_l3"][key]
                d = self.rhs[2]["primary_l3"][key]
                increment = .5 * dt * (a + b + c + d)
                initial = getattr(state.geometry, key)[:CELLS]
                stored = self.pre_boundary[key]
                np.testing.assert_array_equal(stored, initial + increment)
                np.testing.assert_array_equal(stored, getattr(result.geometry, key)[:CELLS])
                pieces = {part: [] for part in PARTS}
                for operands, old, new, performed in zip(zip(a, b, c, d), initial, stored, increment):
                    computed = wide(dt) * sum(map(wide, operands)) / 2
                    pieces["computed"].append(computed)
                    pieces["increment_arithmetic"].append(wide(performed) - computed)
                    pieces["write"].append(wide(new) - wide(old) - wide(performed))
                self.latest[key] = pieces
                for part in PARTS:
                    self.totals[key][part] = [x + y for x, y in zip(self.totals[key][part], pieces[part])]
        return result

    def append_row(self, grid, geometry):
        n = CELLS - 2
        r = grid.centers[:n]
        coefficient = self.shear_coefficient(grid, geometry)
        row = self.stages["accepted_end"]
        initial = self.initial_rows
        with localcontext() as ctx:
            ctx.prec = 60
            for field in FIELDS:
                endpoint = decimals(getattr(geometry, field))
                totals = self.totals[field]
                for idx, values in enumerate(zip(*(totals[part] for part in PARTS))):
                    if abs(sum(values) - (endpoint[idx] - self.initial[field][idx])) > Decimal("1e-54"):
                        raise ValueError("stored field updates do not telescope")

            def response(field, values):
                dr = wide(grid.dr)
                d = [derivative(values, i, dr) for i in range(n)]
                if field == "Aa":
                    return np.array([float(v + wide(c) * w)
                                     for v, c, w in zip(d, coefficient, values)])
                return np.array([float(-2 * v / 3) for v in d])

            cumulative = {field: {part: response(field, self.totals[field][part])
                                  for part in PARTS} for field in FIELDS}
            last_write = {field: response(field, self.latest[field]["write"]) for field in FIELDS}

        # Use the current shear coefficient for all Aa increments, and account
        # for the response of the initial Aa separately.
        coefficient_term = (coefficient - self.initial_c)[:n] * np.asarray(self.initial["Aa"][:n], dtype=float)
        sources = {key: (row[key] - initial[key])[:n] for key in replay.PARTS[4:]}
        readback = row["M"][:n]
        actual = readback - initial["M"][:n]
        joint_write = sum(cumulative[field]["write"] for field in FIELDS)
        predicted = (sum(sum(parts.values()) for parts in cumulative.values())
                     + coefficient_term + sum(sources.values()))
        closure = actual - predicted
        # The ordinary derivative evaluates on the full background K;
        # its subtraction error scales with |K|/dr, not with the small delta.
        evaluation_bound = 4 * np.finfo(float).eps * np.max(abs(row["K"])) / grid.dr
        if np.max(abs(closure)) > evaluation_bound:
            raise ValueError("accumulated momentum budget does not close")

        def stats(values):
            return {"signed": values.tolist(), "band_rms": replay.local_rms(r, values)}

        self.rows.append({
            "t": self.t, "M": stats(readback), "delta_M": stats(actual),
            "M_signs": signs(r, readback),
            "cumulative": {field: {part: stats(v) for part, v in pieces.items()}
                           for field, pieces in cumulative.items()},
            "last_write": {field: stats(v) for field, v in last_write.items()},
            "joint_write": stats(joint_write), "M_minus_direct_write": stats(readback - joint_write),
            "shear_coefficient_change": stats(coefficient_term),
            "matter_change": {key: stats(v) for key, v in sources.items()},
            "signed_budget_closure": stats(closure),
        })


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dr", type=float, choices=(.025, .0125))
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    meta = replay.provenance()
    go, vacuum, _ = replay.pins.ad.vendor_modules()
    if vacuum.KO_EPSILON != 0:
        raise ValueError("expected KO off")
    observer = MomentumHistory(vacuum)
    advance = replay.advance
    with (patch.object(replay, "Observe", lambda _: observer),
          patch.object(replay, "advance", lambda kernel, state, dt, watch: advance(kernel, state, dt, True))):
        report = replay.fine_replay(args.dr, go, vacuum, meta)
    archive = json.loads((replay.FOLDER / "comparison.json").read_text())
    samples = archive["runs"][str(args.dr)]["0.000375"]["rows"]
    for i, target in enumerate(samples[1:], 1):
        np.testing.assert_equal(observer.rows[10*i-1]["M"]["signed"][0], target["cell_0"]["M"])
    if len(observer.rows) != 40:
        raise ValueError("expected 40 accepted steps")
    report.pop("observed_steps")
    report["r"] = report["r"][:CELLS-2]
    report.update(rows=observer.rows, initial_M=observer.initial_rows["M"][:CELLS-2].tolist(),
                  engine=replay.pins.SOURCE, vendor=replay.pins.VENDOR,
                  script_sha256=replay.sha(Path(__file__)))
    args.output.write_text(json.dumps(report, allow_nan=False, separators=(",", ":")) + "\n")


if __name__ == "__main__":
    main()
