"""Branch one late W1 production step with alternate research-only lapse solves."""
import argparse
from copy import deepcopy
from dataclasses import fields
import hashlib
import json
from pathlib import Path
from unittest.mock import patch

import numpy as np

from engine.cmc_gauge import total_matter_projection
from engine.stress_energy import assemble_total_stress_energy
from research.w1_fring_smoke import matched_lapse, stage_momentum as replay
from research.w1_fring_smoke.capture_lapse_stages import FINAL_STAGES


SOLVE_STAGES = ("start", "explicit_predictor", "primary_predictor", "accepted_end")


def state_arrays(state):
    for name in ("geometry", "scalars"):
        group = getattr(state, name)
        for field in fields(group):
            yield f"{name}.{field.name}", getattr(group, field.name)
    for species in ("dark_matter", "baryons", "radiation"):
        group = getattr(state.matter, species)
        for name in ("rest", "energy_t", "momentum_r"):
            yield f"matter.{species}.{name}", getattr(group, name)
    yield "clock", np.array([state.t, state.tau, state.e_folds])


def state_digest(state):
    digest = hashlib.sha256()
    for name, values in state_arrays(state):
        digest.update(name.encode())
        digest.update(np.asarray(values, dtype="<f8").tobytes())
    return digest.hexdigest()


class AlternateLapse(replay.Observe):
    def __init__(self, vacuum, order):
        super().__init__(vacuum)
        self.order = order
        self.solves = {}

    def _solve_lapse_for_state(self, grid, geometry, scalars, matter, electromagnetic):
        if electromagnetic is not None or np.any(geometry.beta):
            raise ValueError("expected the zero-shift, EM-free witness")
        alpha, target = super()._solve_lapse_for_state(
            grid, geometry, scalars, matter, electromagnetic
        )
        name = SOLVE_STAGES[len(self.solves)]
        scratch = geometry.copy()
        scratch.alpha = alpha
        stress = total_matter_projection(grid, geometry, scalars, matter)
        trace = stress["rho"] + stress["pr"] + 2.0 * stress["pt"]
        lapse, info = matched_lapse.lapse_variant(
            grid, scratch, (trace, None), target, self.order
        )
        self.solves[name] = {"projected_Kdot": target, **info}
        return lapse, target


def compare_stages(r, baseline, other):
    changes = {}
    for name, reference in baseline.items():
        signed = {key: np.asarray(other[name]["rows"][key]) - reference["rows"][key]
                  for key in (*replay.PARTS, "M")}
        geometry = sum(signed[key] for key in replay.PARTS[:4])
        matter = sum(signed[key] for key in replay.PARTS[4:])
        changes[name] = {
            "rms": {key: replay.local_rms(r, value) for key, value in signed.items()},
            "signed": {key: value.tolist() for key, value in signed.items()},
            "geometry_rms": replay.local_rms(r, geometry),
            "matter_rms": replay.local_rms(r, matter),
            "closure_max": float(np.max(abs(signed["M"] - geometry - matter))),
        }
    return changes


def end_constraints(vacuum, state):
    geometry, grid = state.geometry, state.grid
    raw = vacuum.constraints(grid, geometry)
    total = assemble_total_stress_energy(grid, geometry, state.scalars, state.matter)
    H = raw["hamiltonian"] - 16 * np.pi * total.rho
    n = replay.LOCAL_CELLS
    M = raw["momentum"] - 8 * np.pi * total.j
    return {
        key: {"centre": float(value[0]), "band_rms": replay.local_rms(grid.centers[:n], value[:n]),
              "domain_max_abs": float(np.max(abs(value)))}
        for key, value in (("H", H), ("M", M))
    }


class BranchLateStep(replay.Observe):
    def __init__(self, vacuum, case):
        super().__init__(vacuum)
        self.end_time = FINAL_STAGES[case]
        self.report = None

    def step(self, state, dt):
        if not self.active or not np.isclose(state.t + dt, self.end_time, atol=2e-17, rtol=0):
            return super().step(state, dt)

        digest = state_digest(state)
        baseline = super().step(state, dt)
        self.record("accepted_end", baseline.grid, baseline.geometry, baseline.scalars, baseline.matter)
        r = state.grid.centers[:replay.LOCAL_CELLS]
        baseline_stages = replay.summarize(r, self.stages)
        self.report = {
            "start_sha256": digest,
            "production": {"stages": baseline_stages, "end_constraints": end_constraints(self.vacuum, baseline)},
        }
        baseline_arrays = dict(state_arrays(baseline))
        for name, order in (("three_point_control", 2), ("fourth_order_match", 4)):
            shadow = AlternateLapse(self.vacuum, order)
            start = deepcopy(state)
            if state_digest(start) != digest:
                raise ValueError("branch starts from a different checkpoint")
            evolved, sampled = replay.advance(shadow, start, dt, True)
            if evolved.t != baseline.t or tuple(shadow.solves) != SOLVE_STAGES:
                raise ValueError("branch did not follow the complete step")
            self.report[name] = {
                "stages": sampled["stages"],
                "changes_vs_production": compare_stages(r, baseline_stages, sampled["stages"]),
                "end_constraints": end_constraints(self.vacuum, evolved),
                "end_field_max_delta": {
                    key: float(np.max(abs(values - baseline_arrays[key])))
                    for key, values in state_arrays(evolved)
                },
                "lapse_solves": shadow.solves,
            }

        if state_digest(state) != digest:
            raise ValueError("live replay modified the original checkpoint")
        return baseline


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=tuple(FINAL_STAGES))
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    meta = replay.provenance()
    go, vacuum, _ = replay.pins.ad.vendor_modules()
    if vacuum.KO_EPSILON != 0:
        raise ValueError("expected KO off")
    observer = BranchLateStep(vacuum, args.case)
    with patch.object(replay, "Observe", lambda _: observer):
        raw = (replay.archived_replay(args.case, vacuum) if args.case == "A1"
               else replay.fine_replay(float(args.case), go, vacuum, meta))

    report = observer.report
    if report is None:
        raise ValueError("expected one late checkpoint")
    sampled = raw.pop("observed_steps")[-1]
    if sampled["stages"] != report["production"]["stages"]:
        raise ValueError("branched production step differs from ordinary replay")
    report.update(raw, t_start=sampled["t_start"], t_end=sampled["t_end"],
                  engine=replay.pins.SOURCE, vendor=replay.pins.VENDOR,
                  script_sha256=replay.sha(Path(__file__)),
                  solver_sha256=replay.sha(Path(matched_lapse.__file__)))
    args.output.write_text(json.dumps(report, allow_nan=False, separators=(",", ":")) + "\n")


if __name__ == "__main__":
    main()
