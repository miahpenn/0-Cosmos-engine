"""Canonical end-to-end campaign runner.

This module is deliberately not invoked from import-time code or CI controls.
It is the final numerical entry point after construction is complete.

A finite final_time is a campaign boundary, not a physical stop condition.
The evolution equations themselves contain no reset, bounce, branch flip,
ejection threshold, or phenomenological feedback source.
"""
from dataclasses import asdict, dataclass, is_dataclass, replace
from pathlib import Path
import json
import math
import numpy as np

from .atomic_io import atomic_write_text
from .adaptive_step import (
    AdaptiveStepOutcome,
    RadiationStepSizeUnderflow,
    advance_with_radiation_admissibility_retries,
)
from .campaign import CampaignConfig, prepare_campaign
from .production_kernel import V55ProductionKernel
from .worldtube import current_residual


@dataclass
class ResolutionResult:
    resolution: int
    status: str
    t_final: float
    tau_final: float
    steps: int
    cycle_events: int
    handoffs: int
    checkpoint_count: int
    turnaround_events: int
    reexpansion_events: int
    trapped_root_count_max: int
    misner_sharp_current_max: float | None
    misner_sharp_current_rms: float | None
    last_diagnostics: dict
    endpoint_roundoff_correction: float = 0.0


def _due_checkpoint_times(
    accepted_time: float,
    next_checkpoint: float,
    checkpoint_interval: float,
) -> tuple[list[float], float]:
    """Return checkpoint ticks reached by an accepted endpoint and the next tick.

    Checkpoint files serialize accepted states, not interpolated states at their
    filename times. A scheduled tick is therefore due only after the accepted
    endpoint reaches it.
    """
    accepted = float(accepted_time)
    next_tick = float(next_checkpoint)
    interval = float(checkpoint_interval)
    due = []
    while accepted >= next_tick:
        due.append(next_tick)
        next_tick += interval
    return due, next_tick


def _assert_finite_payload(value, label: str) -> None:
    """Fail immediately when an evolved state or saved diagnostic is non-finite."""
    if isinstance(value, dict):
        for key, item in value.items():
            _assert_finite_payload(item, f"{label}.{key}")
        return
    if isinstance(value, (list, tuple)):
        for index, item in enumerate(value):
            _assert_finite_payload(item, f"{label}[{index}]")
        return
    if isinstance(value, (np.ndarray, np.number, int, float, complex)):
        array = np.asarray(value)
        if array.dtype.kind in "fciub" and not bool(np.all(np.isfinite(array))):
            raise FloatingPointError(f"non-finite value detected in {label}")
        return
    if hasattr(value, "__dict__"):
        for key, item in vars(value).items():
            _assert_finite_payload(item, f"{label}.{key}")


def assert_finite_state(state) -> None:
    """Check all numeric geometry, scalar and conservative-matter state arrays."""
    for label in ("geometry", "scalars", "matter"):
        _assert_finite_payload(getattr(state, label), f"state.{label}")


def _canonicalize_endpoint_roundoff(state, target_time: float, tolerance: float) -> float:
    """Snap an accepted endpoint and every ledger record for that same slice.

    This does not take an additional evolution step. It normalizes the
    coordinate-time label after the nominal last step and keeps the matching
    history, handoff, and cycle-event timestamps synchronized. Older ledger
    entries are left untouched.
    """
    previous_time = float(state.t)
    correction = float(target_time) - previous_time
    if abs(correction) <= tolerance:
        state.t = float(target_time)
        if state.history:
            state.history[-1]["t"] = float(target_time)

        # The accepted endpoint is recorded independently in immutable
        # HandoffPoint objects and CycleEvent objects. Snap only a terminal
        # record whose time exactly matches the pre-correction state time;
        # never relabel a prior event/handoff.
        record_lists = (
            getattr(state, "handoffs", None),
            getattr(getattr(state, "cycle", None), "events", None),
        )
        for records in record_lists:
            if not records:
                continue
            record = records[-1]
            try:
                record_time = float(record.t)
            except (AttributeError, TypeError, ValueError, OverflowError):
                continue
            if record_time != previous_time:
                continue
            if is_dataclass(record):
                records[-1] = replace(record, t=float(target_time))
            else:
                setattr(record, "t", float(target_time))
        return correction
    return 0.0


def campaign_admission_errors(results, config) -> list[str]:
    """Return fail-closed admission errors for a requested campaign.

    A trajectory counts as completed only if every requested resolution is
    present, each run reports completed, and each final coordinate time reaches
    the requested campaign horizon within floating-point endpoint tolerance.
    """
    errors: list[str] = []
    expected = tuple(int(n) for n in config.resolutions)
    observed = tuple(int(result.resolution) for result in results)
    if observed != expected:
        errors.append(
            f"resolution results do not match request: expected {expected}, got {observed}"
        )

    target_time = float(config.final_time)
    if not math.isfinite(target_time) or target_time <= 0.0:
        errors.append(f"invalid requested final_time={target_time!r}")
        time_tolerance = 0.0
    else:
        time_tolerance = 1.0e-10 * max(1.0, abs(target_time))

    for result in results:
        label = f"N={result.resolution}"
        if result.status != "completed":
            errors.append(f"{label} status is {result.status!r}")
        final_time = float(result.t_final)
        if not math.isfinite(final_time):
            errors.append(f"{label} final time is not finite: {final_time!r}")
        elif math.isfinite(target_time) and abs(final_time - target_time) > time_tolerance:
            errors.append(
                f"{label} ended at t={final_time:.17g}, "
                f"short of requested t={target_time:.17g}"
            )

        # None means the current residual is undefined because too few time
        # samples exist. A numeric value, when supplied, must be finite.
        for field in (
            "misner_sharp_current_max",
            "misner_sharp_current_rms",
        ):
            value = getattr(result, field, None)
            if value is None:
                continue
            try:
                metric = float(value)
            except (TypeError, ValueError, OverflowError):
                errors.append(f"{label} {field} is not a finite number: {value!r}")
                continue
            if not math.isfinite(metric):
                errors.append(f"{label} {field} is not finite: {value!r}")
    return errors


def _json_safe(value):
    if isinstance(value, dict):
        return {str(k): _json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(v) for v in value]
    if isinstance(value, np.ndarray):
        # Recurse through tolist() so non-finite array leaves are normalized too.
        return _json_safe(value.item() if value.ndim == 0 else value.tolist())
    if isinstance(value, np.generic):
        return _json_safe(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def _write_ledgers(
    state,
    run_dir: Path,
    result: ResolutionResult,
    radiation_step_retries: list[dict] | None = None,
) -> None:
    run_dir.mkdir(parents=True, exist_ok=True)
    # An initialization failure has no physical state to serialize, but it
    # still needs a durable fail-closed ledger and retry log.
    events = [] if state is None else state.cycle.events
    handoffs = [] if state is None else state.handoffs
    history = [] if state is None else state.history
    payload = {
        "result": _json_safe(asdict(result)),
        "cycle_events": [_json_safe(asdict(e)) for e in events],
        "handoffs": [_json_safe(asdict(h)) for h in handoffs],
        "history": _json_safe(history),
        "radiation_step_retries": _json_safe(radiation_step_retries or []),
    }
    atomic_write_text(
        run_dir / "ledgers.json", json.dumps(payload, indent=2, allow_nan=False)
    )
    atomic_write_text(
        run_dir / "radiation_step_retries.json",
        json.dumps(_json_safe(radiation_step_retries or []), indent=2, allow_nan=False),
    )


def run_campaign(
    kernel: V55ProductionKernel | None = None,
    config: CampaignConfig = CampaignConfig(),
) -> list[ResolutionResult]:
    """Run all requested resolutions after validating the production contract."""
    if kernel is None:
        kernel = V55ProductionKernel()
    # Validate the production capability contract and campaign inputs before
    # any expensive evolution begins.
    prepare_campaign(kernel, config)

    results: list[ResolutionResult] = []
    interruption: KeyboardInterrupt | None = None
    for resolution in config.resolutions:
        try:
            state = kernel.initialize(
                resolution=resolution,
                r_max=config.r_max,
                include_radiation=True,
            )
        except KeyboardInterrupt as exc:
            # No state is available to continue from; record the interrupted
            # resolution and stop before beginning another initialization.
            result = ResolutionResult(
                resolution=resolution,
                status="interrupted:initialization",
                t_final=0.0,
                tau_final=0.0,
                steps=0,
                cycle_events=0,
                handoffs=0,
                checkpoint_count=0,
                turnaround_events=0,
                reexpansion_events=0,
                trapped_root_count_max=0,
                misner_sharp_current_max=None,
                misner_sharp_current_rms=None,
                last_diagnostics={
                    "interruption_class": type(exc).__name__,
                    "interruption_message": str(exc),
                },
            )
            run_dir = Path(config.output_dir) / f"N{resolution}"
            _write_ledgers(None, run_dir, result, [])
            results.append(result)
            interruption = exc
            break
        except (FloatingPointError, ValueError, TypeError, RuntimeError) as exc:
            # Initialization can reject an inadmissible radiation seed before
            # a ProductionState exists. Persist that rejection as a failed
            # resolution so the campaign summary is still written fail-closed.
            result = ResolutionResult(
                resolution=resolution,
                status="numerical_failure:initialization:" + str(exc),
                t_final=0.0,
                tau_final=0.0,
                steps=0,
                cycle_events=0,
                handoffs=0,
                checkpoint_count=0,
                turnaround_events=0,
                reexpansion_events=0,
                trapped_root_count_max=0,
                misner_sharp_current_max=None,
                misner_sharp_current_rms=None,
                last_diagnostics={"initialization_error": str(exc)},
            )
            run_dir = Path(config.output_dir) / f"N{resolution}"
            _write_ledgers(None, run_dir, result, [])
            results.append(result)
            continue

        target_time = float(config.final_time)
        endpoint_tolerance = 1.0e-12 * max(1.0, abs(target_time))
        endpoint_correction = 0.0
        dt_nominal = config.cfl * state.grid.dr
        next_checkpoint = config.checkpoint_interval
        checkpoint_count = 0
        radiation_step_retries: list[dict] = []
        status = "completed"

        # Include the initial slice so the evolution ledger has an actual t=0
        # baseline rather than starting at the first post-step sample.
        try:
            assert_finite_state(state)
            initial = kernel.diagnostics(state)
            if isinstance(initial, dict):
                initial = {
                    key: value for key, value in initial.items()
                    if key not in (
                        "total_rho", "total_pr", "total_pt", "total_j",
                        "areal_radius", "chi", "misner_sharp"
                    )
                }
            else:
                raise TypeError("kernel diagnostics must return a mapping")
            initial["t"] = float(state.t)
            initial["tau"] = float(state.tau)
            initial["cycle_event"] = None
            _assert_finite_payload(initial, f"initial_diagnostics_N{resolution}")
            state.history.append(initial)
        except KeyboardInterrupt as exc:
            status = "interrupted:initial_diagnostics"
            interruption = exc
        except (FloatingPointError, ValueError, TypeError, RuntimeError) as exc:
            status = "numerical_failure:initial_state:" + str(exc)

        while status == "completed" and state.t < target_time:
            remaining = target_time - state.t
            if remaining <= endpoint_tolerance:
                endpoint_correction = _canonicalize_endpoint_roundoff(
                    state, target_time, endpoint_tolerance
                )
                break
            requested_dt = min(dt_nominal, remaining)
            state_before_time = float(state.t)
            try:
                outcome: AdaptiveStepOutcome = advance_with_radiation_admissibility_retries(
                    kernel, state, requested_dt
                )
                candidate_state = outcome.state
                # A successful cone retry is only a step candidate. Retain the
                # strict whole-state and diagnostic checks before acceptance.
                assert_finite_state(candidate_state)
                candidate_state.geometry.assert_finite_positive()
                if not candidate_state.history:
                    raise FloatingPointError("kernel step produced no diagnostic history")
                _assert_finite_payload(
                    candidate_state.history[-1],
                    f"diagnostics_N{resolution}_t{candidate_state.t:.17g}",
                )
                state = candidate_state
                dt = outcome.accepted_dt
                if outcome.retry_count:
                    radiation_step_retries.append({
                        "t_before": state_before_time,
                        "requested_dt": outcome.requested_dt,
                        "accepted_dt": outcome.accepted_dt,
                        "retry_count": outcome.retry_count,
                        "attempts": list(outcome.attempts),
                        "reason": "typed_radiation_admissibility_rejection",
                    })
                correction = _canonicalize_endpoint_roundoff(
                    state, target_time, endpoint_tolerance
                )
                if correction:
                    endpoint_correction = correction
            except KeyboardInterrupt as exc:
                # The candidate is not accepted until the adaptive step call
                # returns. Keep the last accepted state and persist its ledger
                # before propagating the interruption to the caller.
                status = "interrupted:step"
                interruption = exc
                break
            except (FloatingPointError, ValueError, TypeError, RuntimeError) as exc:
                if isinstance(exc, RadiationStepSizeUnderflow):
                    radiation_step_retries.append({
                        "t_before": state_before_time,
                        "requested_dt": requested_dt,
                        "accepted_dt": None,
                        "retry_count": len(exc.attempts),
                        "attempts": list(exc.attempts),
                        "status": "failed_closed_at_time_precision",
                        "exception_class": type(exc).__name__,
                        "exception_message": str(exc),
                    })
                status = "numerical_failure:" + str(exc)
                break

            due_checkpoints, next_checkpoint = _due_checkpoint_times(
                float(state.t), next_checkpoint, config.checkpoint_interval
            )
            for checkpoint_time in due_checkpoints:
                path = (
                    Path(config.output_dir)
                    / f"N{resolution}"
                    / f"checkpoint_t{checkpoint_time:012.6f}.npz"
                )
                try:
                    kernel.checkpoint(state, path)
                except KeyboardInterrupt as exc:
                    # The atomic checkpoint writer preserves the previous
                    # canonical file if staging is interrupted. Still record
                    # the last accepted state and a fail-closed campaign result.
                    status = "interrupted:checkpoint"
                    interruption = exc
                    break
                checkpoint_count += 1

        if state.history:
            last = state.history[-1]
        elif status != "completed":
            # The returned initial state or its initial diagnostics already
            # failed validation. Do not recompute diagnostics on that invalid
            # state: it could raise again or produce non-finite values before
            # the fail-closed ledger/summary is written.
            last = {
                "diagnostics_available": False,
                "failure_status": status,
            }
        else:
            try:
                last = kernel.diagnostics(state)
                if not isinstance(last, dict):
                    raise TypeError("kernel diagnostics must return a mapping")
                _assert_finite_payload(last, f"final_diagnostics_N{resolution}")
            except (FloatingPointError, ValueError, TypeError, RuntimeError) as exc:
                status = "numerical_failure:final_diagnostics:" + str(exc)
                last = {
                    "diagnostics_available": False,
                    "diagnostics_error": str(exc),
                }
        if len(state.history) >= 3:
            times = np.asarray([row["t"] for row in state.history])
            masses = np.asarray([row["M_MS"] for row in state.history])
            rhs = np.asarray([
                row["flux_T"] + row["work_pR"] for row in state.history
            ])
            ms_residual = np.asarray(
                current_residual(times, masses, rhs), dtype=float
            )
            if ms_residual.size == 0 or not bool(np.all(np.isfinite(ms_residual))):
                # A computed but invalid residual is a numerical failure, not
                # an undefined metric and never a passing campaign result.
                status = "numerical_failure:invalid_misner_sharp_current_residual"
                ms_max = None
                ms_rms = None
            else:
                ms_max = float(np.max(np.abs(ms_residual)))
                ms_rms = float(np.sqrt(np.mean(ms_residual**2)))
        else:
            # The discrete current residual needs three distinct time samples.
            # Use None so output is standard JSON null, not the non-standard NaN.
            ms_max = None
            ms_rms = None
        turnaround = sum(e.kind == "turnaround" for e in state.cycle.events)
        reexpansion = sum(e.kind == "re_expansion_crossing" for e in state.cycle.events)
        max_roots = max(
            (len(row.get("trapped_roots", [])) for row in state.history),
            default=0,
        )
        result = ResolutionResult(
            resolution=resolution,
            status=status,
            t_final=float(state.t),
            tau_final=float(state.tau),
            steps=max(0, len(state.history) - 1),
            endpoint_roundoff_correction=float(endpoint_correction),
            cycle_events=len(state.cycle.events),
            handoffs=len(state.handoffs),
            checkpoint_count=checkpoint_count,
            turnaround_events=turnaround,
            reexpansion_events=reexpansion,
            trapped_root_count_max=max_roots,
            misner_sharp_current_max=ms_max,
            misner_sharp_current_rms=ms_rms,
            last_diagnostics=_json_safe(last),
        )
        run_dir = Path(config.output_dir) / f"N{resolution}"
        _write_ledgers(state, run_dir, result, radiation_step_retries)
        results.append(result)
        if interruption is not None:
            # Do not begin additional resolutions once an operator interrupts.
            break

    admission_errors = campaign_admission_errors(results, config)
    summary = {
        "campaign": _json_safe(asdict(config)),
        "results": _json_safe([asdict(r) for r in results]),
        "admission": {
            "passed": not admission_errors,
            "errors": admission_errors,
        },
    }
    out = Path(config.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    atomic_write_text(
        out / "campaign_summary.json", json.dumps(summary, indent=2, allow_nan=False)
    )
    if interruption is not None:
        # The fail-closed summary/ledger has now been written. Preserve normal
        # Ctrl-C semantics rather than disguising an interruption as a numerical
        # admission exception.
        raise interruption
    if admission_errors:
        raise RuntimeError(
            "Production campaign admission failed: " + "; ".join(admission_errors)
        )
    return results


if __name__ == "__main__":
    raise SystemExit(
        "Campaign runner is intentionally not auto-executed during build. "
        "Invoke run_campaign from the final campaign job."
    )
