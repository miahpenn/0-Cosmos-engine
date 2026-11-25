"""Capture the operands of the last W1 corrector on a pinned replay."""
import argparse
from pathlib import Path
from unittest.mock import patch

import numpy as np

from engine.cmc_gauge import total_matter_projection
from research.w1_fring_smoke import stage_momentum as replay


FINAL_STAGES = {"0.025": .015, "0.0125": .015, "A1": .0165}
FIELDS = ("a", "b", "X", "Aa", "K", "Lambda", "alpha", "beta", "B")
PRIMARY_FIELDS = ("K", "Aa")


class Capture(replay.Observe):
    def __init__(self, vacuum, case):
        super().__init__(vacuum)
        self.end = FINAL_STAGES[case]
        self.capturing = False
        self.data = {}

    def _solve_lapse_for_state(self, grid, geometry, scalars, matter, electromagnetic):
        alpha, target = super()._solve_lapse_for_state(
            grid, geometry, scalars, matter, electromagnetic
        )
        if self.capturing:
            self.target = target
        return alpha, target

    def _geometry_stage_terms(self, grid, geometry, scalars, matter, electromagnetic):
        call = self.geometry_calls
        terms = super()._geometry_stage_terms(grid, geometry, scalars, matter, electromagnetic)
        if self.capturing and call in (0, 2):
            if electromagnetic is not None or np.any(geometry.beta):
                raise ValueError("expected the zero-shift, EM-free witness")
            name = {0: "start", 2: "primary"}[call]
            self.save_geometry(name, geometry)
            self.data[f"{name}_target"] = self.target
            projection = total_matter_projection(grid, geometry, scalars, matter)
            for key in ("rho", "pr", "pt"):
                self.data[f"{name}_{key}"] = projection[key].copy()
            for block in ("primary_l2", "primary_l3"):
                for field in PRIMARY_FIELDS:
                    self.data[f"{name}_{block}_{field}"] = terms[block][field].copy()
        return terms

    def save_geometry(self, name, geometry):
        for field in FIELDS:
            self.data[f"{name}_{field}"] = getattr(geometry, field).copy()

    def record(self, name, grid, geometry, scalars, matter):
        super().record(name, grid, geometry, scalars, matter)
        if self.capturing and name == "corrector_before_boundary":
            self.data["corrector_before_boundary_Aa"] = geometry.Aa[:replay.LOCAL_CELLS].copy()

    def step(self, state, dt):
        self.capturing = self.active and np.isclose(state.t + dt, self.end, atol=2e-17, rtol=0)
        if not self.capturing:
            return super().step(state, dt)

        calls = 0
        primary_l2_rhs = self.vacuum.primary_l2_rhs

        def capture_corrector(grid, geometry):
            nonlocal calls
            terms = primary_l2_rhs(grid, geometry)
            if "primary_target" in self.data and self.boundary_calls == 1:
                calls += 1
                self.save_geometry("corrector_copy", geometry)
                for field in PRIMARY_FIELDS:
                    self.data[f"corrector_copy_primary_l2_{field}"] = terms[field].copy()
            return terms

        with patch.object(self.vacuum, "primary_l2_rhs", capture_corrector):
            result = super().step(state, dt)
        if calls != 1:
            raise ValueError("expected a single corrector L2 evaluation")
        for name in ("initial_lapse", "corrector_before_boundary"):
            for key, value in self.stages[name].items():
                self.data[f"{name}_{key}"] = value
        self.data.update(dt=dt, dr=state.grid.dr, t_start=state.t, t_end=result.t)
        self.capturing = False
        return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=tuple(FINAL_STAGES))
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    meta = replay.provenance()
    go, vacuum, _ = replay.pins.ad.vendor_modules()
    if vacuum.KO_EPSILON != 0:
        raise ValueError("expected KO off")
    observer = Capture(vacuum, args.case)
    with patch.object(replay, "Observe", lambda _: observer):
        if args.case == "A1":
            report = replay.archived_replay(args.case, vacuum)
        else:
            report = replay.fine_replay(float(args.case), go, vacuum, meta)
    if not observer.data:
        raise ValueError("no corrector snapshot collected")
    digest = report.get("input_sha256", report.get("state_sha256"))
    np.savez_compressed(args.output, **observer.data, input_sha256=digest,
                        engine=replay.pins.SOURCE, vendor=replay.pins.VENDOR,
                        capture_sha256=replay.sha(Path(__file__)))


if __name__ == "__main__":
    main()
