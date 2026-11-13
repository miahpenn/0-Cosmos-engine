"""DECISION-map-row-5 (Nov 11 2026): the GEAR-63/65 reduced lane is provenance only; the GEAR-67 quartic is the one
arrest closure. The module is kept, unchanged in behaviour, and no production kernel path may import it."""
import math
import pathlib
import re

import engine.local as local
from engine.run_controls import run_controls

PRODUCTION_KERNEL_PATH = (
    "production_kernel.py", "scalar_system.py", "v55_initial.py", "v55_matter.py", "v55_pirk_adapter.py",
    "matter_system.py", "matter_rhs.py", "cmc_gauge.py", "stress_energy.py", "valencia.py", "worldtube.py",
    "coupled_engine.py", "bidirectional_production.py",
)


def test_module_is_marked_provenance_only():
    assert local.__doc__.startswith("PROVENANCE ONLY.")
    assert "DECISION-map-row-5" in local.__doc__


def test_no_production_kernel_module_imports_the_reduced_lane():
    root = pathlib.Path(local.__file__).parent
    pat = re.compile(r"^\s*(from\s+(\.|engine\.)local\s+import|import\s+engine\.local|from\s+(\.|engine)\s+import\s+local\b)", re.M)
    for name in PRODUCTION_KERNEL_PATH:
        assert not pat.search((root / name).read_text()), name


def test_behaviour_unchanged_for_provenance_comparisons():
    p = local.LocalParams()
    y = (0.02, 0.0, 1.0e-8, 0.0)
    Sd, Sdd, Dd, Ddd = local.local_rhs(y, p, H=0.0)
    assert Sd == 0.0 and Dd == 0.0
    assert Sdd == -(p.kappa + p.g) * 0.02
    assert Ddd == (p.g - p.kappa) * 1.0e-8 - (2.0 * math.pi / 3.0) * p.G * p.alpha_D**2 * 1.0e-24


def test_cubic_is_the_quartic_at_lambda_four_pi_over_three():
    # -(2 pi/3) G a^2 D^3 = -d/dD [(lambda/8) D^4] with lambda = (4 pi/3) G a^2; D*^2 = 2/lambda = 3/(2 pi)
    lam = 4.0 * math.pi / 3.0
    D = 0.37
    assert math.isclose((2.0 * math.pi / 3.0) * D**3, (lam / 2.0) * D**3, rel_tol=1e-15)
    assert math.isclose(math.sqrt(2.0 / lam), 0.690988, rel_tol=1e-6)


def test_control_runner_labels_the_lane():
    assert run_controls()["local"]["status"] == "PROVENANCE_ONLY_GEAR63_65_RETIRED"
