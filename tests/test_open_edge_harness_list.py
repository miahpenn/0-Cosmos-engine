"""The open-edge harness list is explicit and guarded (Jeremiah's review of PR #507, Nov 16 2026).

The eleven contract modules written against the open edge's internals get outer_edge="open_edge" from tests/conftest.py,
centrally, so their files stay byte-identical (some are path triggers of source-pinned workflows). This module pins that
list, and fails if a test module that touches the open edge's internals and builds a production state is in neither
harness list, nor asks for its edge explicitly, nor is reviewed below as safe on the wall. A future test written against
the open edge's internals therefore cannot land silently on the wall.
"""
import pathlib

import conftest

TESTS = pathlib.Path(__file__).resolve().parent

EXPECTED_OPEN_EDGE_MODULES = {
    "tests/test_radiation_rhs_admissibility_contract.py",
    "tests/test_strict_radiation_rhs_admissibility_contract.py",
    "tests/test_native_cmc_explicit_rhs_production_wiring.py",
    "tests/test_constraint_source_normalization_contract.py",
    "tests/test_native_cmc_explicit_rhs_stage_shadow.py",
    "tests/test_native_constraints.py",
    "tests/test_v55_frozen_slice_rhs_shadow.py",
    "tests/test_v55_lambda_l3_matter_contract.py",
    "tests/test_research_puncture_gauge_kernel.py",
    "tests/test_em_consumer_gain_contract.py",
    "tests/test_em_preflight_baseline_projection_parity.py",
}

# Names whose use marks a test as written against the open edge's implementation seams.
OPEN_EDGE_SEAMS = ("_solve_lapse(", "species_rhs", "_apply_outer_light_boundary", "solve_cmc_lapse", "target_kdot",
                   "SphericalCellGrid")

# Reviewed: touch a seam name but are correct on whatever edge their start gets.
REVIEWED_WALL_SAFE = {
    "tests/test_production_derived_clock_quadrature.py":
        "calls the kernel's _solve_lapse, which dispatches on the run's edge (wall: K-dot from the wall's clock)",
    "tests/test_em_lorentz_rhs_injection.py":
        "radial-Maxwell (EM) start: the edge default is the open edge for EM runs (lattice_wall.default_edge)",
    "tests/test_lattice_wall_kernel.py": "the wall's own contracts; asks for each edge explicitly",
    "tests/test_open_edge_harness_list.py": "this guard",
}


def test_open_edge_list_is_exactly_the_eleven_and_they_exist():
    assert set(conftest.OPEN_EDGE_REFERENCE_MODULES) == EXPECTED_OPEN_EDGE_MODULES
    assert len(EXPECTED_OPEN_EDGE_MODULES) == 11
    for rel in EXPECTED_OPEN_EDGE_MODULES:
        assert (TESTS.parent / rel).is_file(), rel
        assert conftest.OPEN_EDGE_REFERENCE_MODULES[rel], f"{rel}: the reason must be stated"


def test_no_module_touching_open_edge_seams_is_unaccounted_for():
    covered = set(conftest.OPEN_EDGE_REFERENCE_MODULES) | set(conftest.HAND_FLOOR_REFERENCE_MODULES) | set(REVIEWED_WALL_SAFE)
    unaccounted = []
    for path in sorted(TESTS.glob("test_*.py")):
        rel = f"tests/{path.name}"
        text = path.read_text(encoding="utf-8")
        if ".initialize(" not in text or rel in covered or "outer_edge" in text:
            continue
        hits = [s for s in OPEN_EDGE_SEAMS if s in text]
        if hits:
            unaccounted.append((rel, hits))
    assert not unaccounted, (
        "test modules touch the open edge's internals, build a production state and have no explicit edge: add them to "
        f"tests/conftest.py OPEN_EDGE_REFERENCE_MODULES (with a reason), pass outer_edge, or review them here: {unaccounted}"
    )
