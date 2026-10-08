"""Non-interference test for the Step 1 Hamiltonian recorder."""

from __future__ import annotations

import hashlib

import numpy as np

from engine.production_kernel import V55ProductionKernel
from engine.step1_hamiltonian_terms import record_sample


def checksum(state):
    h = hashlib.sha256()
    for label, obj in (
        ("geometry", state.geometry),
        ("scalars", state.scalars),
        ("matter", state.matter),
    ):
        h.update(label.encode())
        for name, value in vars(obj).items():
            if isinstance(value, np.ndarray):
                a = np.ascontiguousarray(value)
                h.update(name.encode())
                h.update(str(a.dtype).encode())
                h.update(str(a.shape).encode())
                h.update(a.tobytes(order="C"))
    for name in ("t", "tau", "e_folds"):
        h.update(name.encode())
        h.update(
            np.asarray(
                [float(getattr(state, name))], dtype=np.float64
            ).tobytes()
        )
    return h.hexdigest()


def evolve(recording):
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=32,
        r_max=16.0,
        width=7.0,
        D_amplitude=1.0e-4,
        include_radiation=True,
    )
    checksums = [checksum(state)]
    dt = 0.0075 * state.grid.dr

    for step in range(1, 4):
        state = kernel.step(state, dt)
        if recording:
            record_sample(state, step)
        checksums.append(checksum(state))
    return checksums


def test_step1_recording_is_bitwise_non_interfering():
    off = evolve(recording=False)
    on = evolve(recording=True)
    assert on == off
