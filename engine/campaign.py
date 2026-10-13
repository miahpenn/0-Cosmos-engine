"""Final numerical campaign orchestration contract.

This module refuses unsafe campaign settings before an evolution loop can
begin. It also blocks the incomplete CMC/reference driver from being promoted
silently into the production campaign path.
"""
from collections.abc import Sequence
from dataclasses import dataclass
from numbers import Integral, Real
from pathlib import Path
import json
import math
import os

from .atomic_io import atomic_write_text
from .production_contract import require_production_capabilities


def _positive_finite(name: str, value) -> float:
    """Validate a strictly positive finite real-valued campaign parameter."""
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite positive real number")
    try:
        numeric = float(value)
    except (OverflowError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite positive real number") from exc
    if not math.isfinite(numeric) or numeric <= 0.0:
        raise ValueError(f"{name} must be a finite positive real number")
    return numeric


def _validated_resolutions(values) -> tuple[int, ...]:
    """Require a non-empty, duplicate-free sequence of positive integers."""
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence):
        raise ValueError("resolutions must be a sequence of positive integers")
    if not values:
        raise ValueError("at least one resolution is required")

    resolutions = []
    for value in values:
        if isinstance(value, bool) or not isinstance(value, Integral):
            raise ValueError("resolutions must contain only positive integers")
        if value <= 0:
            raise ValueError("resolutions must contain only positive integers")
        resolutions.append(int(value))

    if len(set(resolutions)) != len(resolutions):
        raise ValueError("resolutions must not contain duplicates")
    return tuple(resolutions)


def _validate_config_values(config) -> dict:
    """Validate configuration values without running or changing physics."""
    resolutions = _validated_resolutions(config.resolutions)
    r_max = _positive_finite("r_max", config.r_max)
    cfl = _positive_finite("cfl", config.cfl)
    final_time = _positive_finite("final_time", config.final_time)
    checkpoint_interval = _positive_finite(
        "checkpoint_interval", config.checkpoint_interval
    )

    # The pinned SphericalCellGrid defines dr = r_max / N. Validate the
    # smallest and largest nominal steps across the requested resolutions
    # before the runner can enter a time/checkpoint loop.
    try:
        min_nominal_dt = cfl * (r_max / max(resolutions))
        max_nominal_dt = cfl * (r_max / min(resolutions))
    except (OverflowError, ZeroDivisionError) as exc:
        raise ValueError(
            "cfl and r_max/resolutions produce an unrepresentable nominal timestep"
        ) from exc
    if not math.isfinite(min_nominal_dt) or min_nominal_dt <= 0.0:
        raise ValueError(
            "cfl and r_max/resolutions produce a zero or non-finite nominal timestep"
        )
    if not math.isfinite(max_nominal_dt) or max_nominal_dt <= 0.0:
        raise ValueError(
            "cfl and r_max/resolutions produce a zero or non-finite nominal timestep"
        )

    # The largest float spacing strictly below final_time is the ULP of
    # its predecessor (important when final_time is an exact power of two).
    # A step no larger than half that spacing can round back to the same t,
    # depending on ties-to-even. Require representable forward progress, not
    # a physics-based CFL floor.
    previous_time = math.nextafter(final_time, -math.inf)
    min_progress_step = math.ulp(previous_time) / 2.0
    if min_nominal_dt <= min_progress_step:
        raise ValueError(
            "cfl and resolution produce a nominal timestep too small to guarantee "
            "floating-point time progress through final_time"
        )

    # The current runner writes the accepted state for each checkpoint tick
    # crossed; it does not interpolate states at checkpoint times. Keep the
    # interval at least as large as the largest nominal step so one accepted
    # step cannot trigger an unbounded catch-up loop of duplicate snapshots.
    if checkpoint_interval < max_nominal_dt:
        raise ValueError(
            "checkpoint_interval must be at least the largest nominal timestep "
            "across resolutions; finer checkpoint times cannot produce distinct "
            "accepted states"
        )

    # The production runner currently formats checkpoint filenames with
    # six decimal places. A finer cadence can make distinct scheduled times
    # resolve to the same path, silently overwriting an earlier snapshot.
    # This is an output-integrity constraint, not a physical CFL threshold.
    if checkpoint_interval < 1.0e-6:
        raise ValueError(
            "checkpoint_interval must be at least 1e-6 because checkpoint "
            "filenames preserve only six decimal places"
        )

    try:
        output_dir = os.fspath(config.output_dir)
    except TypeError as exc:
        raise ValueError("output_dir must be a non-empty string or path") from exc
    if not isinstance(output_dir, str) or not output_dir.strip():
        raise ValueError("output_dir must be a non-empty string or path")
    return {
        "resolutions": resolutions,
        "r_max": r_max,
        "cfl": cfl,
        "final_time": final_time,
        "checkpoint_interval": checkpoint_interval,
        "output_dir": output_dir,
    }


@dataclass(frozen=True)
class CampaignConfig:
    # Preserve the original dr values while moving the finite-radius
    # worldtube outward so the strong outgoing D geometry remains interior.
    resolutions: tuple = (80, 120, 160)
    r_max: float = 80.0
    cfl: float = 0.03
    final_time: float = 160.0
    checkpoint_interval: float = 2.0
    output_dir: str = "runs"

    def __post_init__(self) -> None:
        # Fail at construction time for settings that could make the evolution
        # or checkpoint loop non-advancing/non-terminating.
        values = _validate_config_values(self)
        # Freeze normalized values too, so a list of resolutions or a Path
        # cannot mutate/serialize differently after the preflight check.
        for name, value in values.items():
            object.__setattr__(self, name, value)


def write_checkpoint(path: Path, metadata: dict, trajectory: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {"metadata": metadata, "trajectory": trajectory}
    atomic_write_text(path, json.dumps(payload, indent=2, default=float))


def prepare_campaign(kernel, config: CampaignConfig = CampaignConfig()):
    """Validate inputs/capabilities and return a campaign plan, without running."""
    values = _validate_config_values(config)
    require_production_capabilities(kernel)
    return {
        **values,
        "status": "prepared_not_executed",
    }
