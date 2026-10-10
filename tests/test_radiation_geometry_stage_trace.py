"""Unit tests for radiation-trace provenance and failure preservation."""
import json
import time

import engine.radiation_geometry_stage_trace as trace_module


def test_push_trace_provenance_separates_source_instrumentation_and_trigger(
    monkeypatch,
):
    source_commit = "1" * 40
    instrumentation_commit = "2" * 40
    trigger_commit = "3" * 40

    monkeypatch.setenv("GITHUB_EVENT_NAME", "push")
    monkeypatch.setenv("GITHUB_REF_NAME", "physics/spatial-beta-covariant-source-closure")
    monkeypatch.setenv("GITHUB_SHA", trigger_commit)

    def fake_git(command):
        key = " ".join(command)
        return {
            "rev-parse HEAD^": instrumentation_commit,
            "rev-parse HEAD~2": source_commit,
        }.get(key, "unavailable")

    monkeypatch.setattr(trace_module, "_git", fake_git)
    provenance = trace_module._trace_provenance()

    assert provenance["branch"] == "physics/spatial-beta-covariant-source-closure"
    assert provenance["run_trigger_commit"] == trigger_commit
    assert provenance["instrumentation_commit"] == instrumentation_commit
    assert provenance["workflow_commit"] == instrumentation_commit
    assert provenance["source_commit_before_trace_workflow"] == source_commit


def test_trace_classifies_only_explicit_radiation_admissibility_as_radiation():
    radiation_failure = {
        "exception_message":
            "radiation conservative state violates E>=|S|: ratio=1.0001"
    }
    cmc_failure = {
        "exception_message":
            "CMC lapse solve has no positive solution on this stage"
    }

    assert (
        trace_module._failure_status(radiation_failure)
        == "radiation_admissibility_failure_captured"
    )
    assert trace_module._failure_status(cmc_failure) == "numerical_failure_captured"


def test_trace_writes_hashed_failure_artifact_when_initialization_has_no_state(
    tmp_path, monkeypatch
):
    monkeypatch.setenv("GITHUB_EVENT_NAME", "workflow_dispatch")
    monkeypatch.setenv("GITHUB_REF_NAME", "physics/spatial-beta-covariant-source-closure")
    monkeypatch.setenv("GITHUB_SHA", "4" * 40)

    trace = trace_module.RadiationStageTrace.__new__(
        trace_module.RadiationStageTrace
    )
    trace.output_dir = tmp_path / "trace"
    trace.output_dir.mkdir(parents=True)
    trace.started_wall = time.time() - 0.25
    trace.step_index = 0
    trace.failure = {
        "stage": "trace_harness",
        "exception_class": "ValueError",
        "exception_message": "test initialization failure",
    }
    trace.global_geometry_minima = {}
    trace.progress_samples = []
    trace.stage_snapshots = []
    trace.predictor_budgets = []
    trace.instrumentation_errors = []

    result = trace.write_outputs(
        None, "trace_harness_failure", float("nan")
    )
    assert result == 0

    payload = json.loads((trace.output_dir / "trace.json").read_text())
    summary = json.loads((trace.output_dir / "summary.json").read_text())
    checksums = (trace.output_dir / "SHA256SUMS.txt").read_text()

    assert payload["final_state"]["t"] is None
    assert payload["final_state"]["tau"] is None
    assert payload["provenance"]["configuration"]["dt_nominal"] is None
    assert payload["provenance"]["branch"] == (
        "physics/spatial-beta-covariant-source-closure"
    )
    assert payload["failure"]["exception_message"] == "test initialization failure"
    assert summary["trace_sha256"]
    assert "trace.json" in checksums and "summary.json" in checksums
