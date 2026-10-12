"""Canonical end-to-end campaign runner.

This module is deliberately not invoked from import-time code or CI controls.
It is the final numerical entry point after construction is complete.

A finite final_time is a campaign boundary, not a physical stop condition.
The evolution equations themselves contain no reset, bounce, branch flip,
ejection threshold, or phenomenological feedback source.
"""
from dataclasses import asdict, dataclass
from pathlib import Path
import json
import math
import numpy as np

from .campaign import CampaignConfig, prepare_campaign
from .production_contract import require_production_capabilities
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
    misner_sharp_current_max: float
    misner_sharp_current_rms: float
    last_diagnostics: dict


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
    return errors


def _json_safe(value):
    if isinstance(value, dict):
        return {str(k): _json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(v) for v in value]
    if isinstance(value, np.ndarray):
        if value.ndim == 0:
            return float(value)
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    return value


def _write_ledgers(
    state,
    run_dir: Path,
    result: ResolutionResult,
) -> None:
    run_dir.mkdir(parents=True, exist_ok=True)
    payload = {
        "result": _json_safe(asdict(result)),
        "cycle_events": [_json_safe(asdict(e)) for e in state.cycle.events],
        "handoffs": [_json_safe(asdict(h)) for h in state.handoffs],
        "history": _json_safe(state.history),
    }
    (run_dir / "ledgers.json").write_text(
        json.dumps(payload, indent=2)
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
    for resolution in config.resolutions:
        state = kernel.initialize(
            resolution=resolution,
            r_max=config.r_max,
            include_radiation=True,
        )
        dt_nominal = config.cfl * state.grid.dr
        next_checkpoint = config.checkpoint_interval
        checkpoint_count = 0
        status = "completed"

        while state.t < config.final_time:
            dt = min(dt_nominal, config.final_time - state.t)
            try:
                state = kernel.step(state, dt)
                # This is a numerical validity check, not a physical event.
                state.geometry.assert_finite_positive()
            except (FloatingPointError, ValueError) as exc:
                status = "numerical_failure:" + str(exc)
                break

            while state.t + 0.5 * dt >= next_checkpoint:
                path = (
                    Path(config.output_dir)
                    / f"N{resolution}"
                    / f"checkpoint_t{next_checkpoint:012.6f}.npz"
                )
                kernel.checkpoint(state, path)
                checkpoint_count += 1
                next_checkpoint += config.checkpoint_interval

        last = state.history[-1] if state.history else kernel.diagnostics(state)
        if len(state.history) >= 3:
            times = np.asarray([row["t"] for row in state.history])
            masses = np.asarray([row["M_MS"] for row in state.history])
            rhs = np.asarray([
                row["flux_T"] + row["work_pR"] for row in state.history
            ])
            ms_residual = current_residual(times, masses, rhs)
            ms_max = float(np.max(np.abs(ms_residual)))
            ms_rms = float(np.sqrt(np.mean(ms_residual**2)))
        else:
            ms_max = float("nan")
            ms_rms = float("nan")
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
            steps=len(state.history),
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
        _write_ledgers(state, run_dir, result)
        results.append(result)

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
    (out / "campaign_summary.json").write_text(
        json.dumps(summary, indent=2)
    )
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
