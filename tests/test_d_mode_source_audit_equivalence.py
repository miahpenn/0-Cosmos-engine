"""Bitwise evolution equivalence: archived inline kernel vs extracted module."""
import math
import subprocess
import textwrap

import numpy as np

from engine.production_kernel import ProductionState, V55ProductionKernel
from engine.v55_initial import build_initial_data
from engine.v55_matter import metric_slice_from_q, project_species
from engine.scalar_system import KAPPA
from engine import v55_pirk_adapter as adapter
from engine.d_mode_source_audit import SourceAuditKernel


BASE_COMMIT = "042d4fd7dd4815cea881b48455c52e8e89e9d3eb"
WORKFLOW = ".github/workflows/zero_star_d_mode_time_shape.yml"
CAMPAIGN_AMPLITUDES = (0.0, 5.0e-11, 1.0e-10, 2.0e-10)
EVOLUTION_STEPS = 20


def _legacy_inline_kernel():
    workflow = subprocess.check_output(
        ["git", "show", f"{BASE_COMMIT}:{WORKFLOW}"], text=True
    )
    start = workflow.index("          class SourceAuditKernel(")
    end = workflow.index("\n          cases = {", start)
    source = textwrap.dedent(workflow[start:end]).replace(
        "class SourceAuditKernel(", "class LegacyInlineKernel(", 1
    )
    namespace = {
        "math": math,
        "np": np,
        "ProductionState": ProductionState,
        "V55ProductionKernel": V55ProductionKernel,
        "build_initial_data": build_initial_data,
        "metric_slice_from_q": metric_slice_from_q,
        "project_species": project_species,
        "KAPPA": KAPPA,
        "adapter": adapter,
    }
    exec(compile(source, "<archived-inline-SourceAuditKernel>", "exec"), namespace)
    return namespace["LegacyInlineKernel"]


def _assert_nested_arrays_identical(left, right, path):
    if isinstance(left, np.ndarray):
        assert isinstance(right, np.ndarray), path
        assert left.dtype == right.dtype, path
        assert left.shape == right.shape, path
        assert left.tobytes() == right.tobytes(), f"bitwise mismatch: {path}"
        return
    if isinstance(left, dict):
        assert isinstance(right, dict), path
        assert left.keys() == right.keys(), path
        for key in left:
            _assert_nested_arrays_identical(left[key], right[key], f"{path}.{key}")
        return
    if isinstance(left, (list, tuple)):
        assert type(left) is type(right) and len(left) == len(right), path
        for i, (a, b) in enumerate(zip(left, right)):
            _assert_nested_arrays_identical(a, b, f"{path}[{i}]")
        return
    if hasattr(left, "__dict__") and not isinstance(left, type):
        assert type(left) is type(right), path
        _assert_nested_arrays_identical(vars(left), vars(right), path)
        return
    if isinstance(left, (float, int, bool, str, type(None), np.generic)):
        assert type(left) is type(right) or (
            isinstance(left, (float, int, np.generic))
            and isinstance(right, (float, int, np.generic))
        ), path
        if isinstance(left, (float, np.floating)):
            assert np.float64(left).tobytes() == np.float64(right).tobytes(), path
        else:
            assert left == right, path


def _assert_evolved_state_identical(left, right, step_index, amplitude):
    prefix = f"D={amplitude:g}.step{step_index}"
    for name in ("geometry", "scalars", "matter"):
        _assert_nested_arrays_identical(
            vars(getattr(left, name)), vars(getattr(right, name)),
            f"{prefix}.{name}",
        )
    for name in ("t", "tau", "e_folds"):
        a, b = getattr(left, name), getattr(right, name)
        assert np.float64(a).tobytes() == np.float64(b).tobytes(), (
            f"{prefix}.{name}: {a!r} != {b!r}"
        )
    _assert_nested_arrays_identical(
        np.asarray(left.grid.centers), np.asarray(right.grid.centers),
        f"{prefix}.grid.centers",
    )


def test_extracted_kernel_matches_archived_inline_evolution_bitwise():
    LegacyInlineKernel = _legacy_inline_kernel()

    for amplitude in CAMPAIGN_AMPLITUDES:
        old = LegacyInlineKernel(amplitude)
        new = SourceAuditKernel(amplitude)
        old_state = old.initialize(resolution=160, r_max=160.0)
        new_state = new.initialize(resolution=160, r_max=160.0)
        _assert_evolved_state_identical(old_state, new_state, "initial", amplitude)

        dt = 0.0075 * old_state.grid.dr
        for step_index in range(1, EVOLUTION_STEPS + 1):
            old_state = old.step(old_state, dt)
            new_state = new.step(new_state, dt)
            _assert_evolved_state_identical(
                old_state, new_state, step_index, amplitude
            )
