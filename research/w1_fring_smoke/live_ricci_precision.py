"""Replay the first W1 step with stable metric second derivatives in primary L2."""
import argparse
from copy import deepcopy
import json
from pathlib import Path
from unittest.mock import patch

import numpy as np

from research.w1_fring_smoke import live_lapse, stage_momentum as replay
from research.w1_fring_smoke.first_step_rate_precision import FirstStepComplete
from research.w1_fring_smoke.frozen_centre_derivative import shell_rms
from research.w1_fring_smoke.live_lapse_audit import signs

ARMS = {"production": (), "stable_X": ("X",), "stable_metric": ("a", "b", "X")}


def centered_second_rows(values, dr):
    """Even reflection at the centre and centred rows through n-3."""
    padded = np.r_[values[:2][::-1], values]
    left2, left1, centre, right1, right2 = (padded[i:i+len(values)-2] for i in range(5))
    return (16*((left1-centre) + (right1-centre))
            - ((left2-centre) + (right2-centre))) / (12*dr**2)


class StableRicci(replay.Observe):
    def __init__(self, vacuum, fields):
        super().__init__(vacuum)
        self.fields = fields
        self.calls = 0

    def step(self, state, dt):
        original_l2 = self.vacuum.primary_l2_rhs
        original_d2 = self.vacuum._d2

        def l2(grid, geometry):
            targets = {id(getattr(geometry, name)) for name in self.fields}
            selected = set()

            def d2(local_grid, values, parity):
                result = original_d2(local_grid, values, parity)
                if id(values) in targets:
                    if parity != 1:
                        raise ValueError("expected even metric reflection")
                    selected.add(id(values))
                    result[:-2] = centered_second_rows(values, local_grid.dr)
                return result

            reference = original_l2(grid, geometry)
            with patch.object(self.vacuum, "_d2", d2):
                result = original_l2(grid, geometry)
            if selected != targets:
                raise ValueError("expected each selected Ricci stencil to be evaluated")
            np.testing.assert_array_equal(reference["K"], result["K"])
            self.calls += 1
            return result

        with patch.object(self.vacuum, "primary_l2_rhs", l2):
            return super().step(state, dt)


def stats(r, profile):
    values = shell_rms(r, profile)
    return {"band_rms": values["rms"], "rough_rms": values["rough_rms"],
            "centre": float(profile[0]), "signs": signs(r, profile)}


def branches(case, go, vacuum, meta):
    original = replay.advance
    reports = {}

    def run_first(kernel, state, dt, watch):
        grid = state.grid
        r = grid.centers[:replay.LOCAL_CELLS]
        digest = live_lapse.state_digest(state)
        for name, fields in ARMS.items():
            observer = StableRicci(vacuum, fields) if fields else kernel
            start = deepcopy(state)
            if live_lapse.state_digest(start) != digest:
                raise ValueError("branch inputs differ")
            accepted, row = original(observer, start, dt, watch)
            reports[name] = {**row, "end_constraints": live_lapse.end_constraints(vacuum, accepted),
                             "end_sha256": live_lapse.state_digest(accepted)}
            if isinstance(observer, StableRicci):
                if observer.calls < 4:
                    raise ValueError("missing primary L2 stage evaluations")
                reports[name]["changed_stencil_evaluations"] = observer.calls
            for stage, data in row["stages"].items():
                data["M_statistics"] = stats(r, np.asarray(data["rows"]["M"]))
        if live_lapse.state_digest(state) != digest:
            raise ValueError("shared checkpoint changed")
        reports["start_sha256"] = digest
        reports["r"] = r.tolist()
        raise FirstStepComplete

    with patch.object(replay, "advance", run_first):
        try:
            replay.fine_replay(float(case), go, vacuum, meta)
        except FirstStepComplete:
            pass
    return reports


def validate_and_compare(report, archive, readback):
    stages = report["production"]["stages"]
    r = np.asarray(report["r"])
    archived = archive["observed_steps"][0]["stages"]
    for stage, data in stages.items():
        reference = {key:value for key,value in data.items() if key != "M_statistics"}
        if reference != archived[stage]:
            raise ValueError("production branch differs from archived first-step stages")

    for name in ("stable_X", "stable_metric"):
        arm = report[name]
        arm["changes_vs_production"] = live_lapse.compare_stages(r, stages, arm["stages"])
        comparisons = {}
        for stage in ("primary_before_boundary", "corrector_before_boundary"):
            frozen = readback["stages"][stage]
            signed_change = np.asarray(arm["changes_vs_production"][stage]["signed"]["M"])
            if name == "stable_X":
                expected = np.asarray(frozen["changes"]["ricci_X_second"]["sum"]["signed"])
            else:
                expected = np.asarray(frozen["centered_second_control"]["signed"])
            count = len(expected)
            difference = signed_change[:count] - expected
            comparisons[stage] = {
                "frozen_signed_change": expected.tolist(),
                "live_minus_frozen_signed": difference.tolist(),
                "live_minus_frozen_band_rms": replay.local_rms(r[:count], difference),
            }
        arm["frozen_comparison"] = comparisons


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=("0.025", "0.0125"))
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    meta = replay.provenance()
    go, vacuum, _ = replay.pins.ad.vendor_modules()
    if vacuum.KO_EPSILON != 0:
        raise ValueError("expected KO off")
    root = Path("results/research")
    archive = root/"w1_fring_stage_momentum"/f"{args.case}.json"
    frozen = root/"w1_fring_first_rate_precision"/f"{args.case}.json"
    reference = json.loads(archive.read_text())
    report = branches(args.case, go, vacuum, meta)
    validate_and_compare(report, reference, json.loads(frozen.read_text()))
    report.update(engine=replay.pins.SOURCE, vendor=replay.pins.VENDOR,
                  input_sha256=reference["input_sha256"], dt=reference["dt"], dr=reference["dr"],
                  script_sha256=replay.sha(Path(__file__)),
                  archive_sha256=replay.sha(archive), frozen_sha256=replay.sha(frozen))
    args.output.write_text(json.dumps(report, allow_nan=False, separators=(",", ":"))+"\n")


if __name__ == "__main__":
    main()
