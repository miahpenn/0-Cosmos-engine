"""The shared model constants exist once (engine/model_constants.py). Gap #2, HANDOFF Nov 14 2026 (5).

The COSMOS lane, the local scalar lane and the production kernel read that one copy. A few files are exact-blob pinned
by the Issue #4 audits and are not edited; their own literals are held equal to the one copy here, so a divergence
fails this test instead of passing silently."""
import inspect
import math

from engine import cosmos, matter_rhs, matter_system, model_constants as mc, production_kernel, scalar_system, valencia


def test_values_are_unchanged():
    assert mc.KAPPA == 8.0 * math.pi
    assert mc.V0 == 8.242522415500654e-5
    assert mc.D_POT == 0.10
    assert mc.LAMBDA_PHI == 1.0
    assert mc.BETA_DM == -0.04


def test_the_lanes_read_the_one_copy():
    p = cosmos.CosmosParams()
    assert (p.V0, p.lam, p.D, p.beta) == (mc.V0, mc.LAMBDA_PHI, mc.D_POT, mc.BETA_DM)
    assert (scalar_system.V0, scalar_system.LAMBDA, scalar_system.D_POT, scalar_system.KAPPA) == (
        mc.V0, mc.LAMBDA_PHI, mc.D_POT, mc.KAPPA)
    assert production_kernel.BETA_DM == mc.BETA_DM
    assert inspect.signature(scalar_system.scalar_rhs_arrays).parameters["beta_dm"].default == mc.BETA_DM
    # no second literal left in the three files that now import the copy
    for module, literals in ((cosmos, ("8.242522415500654e-5", "-0.04")), (scalar_system, ("8.242522415500654e-5", "-0.04")),
                             (production_kernel, ("BETA_DM = -0.04",))):
        source = inspect.getsource(module)
        for literal in literals:
            assert literal not in source, (module.__name__, literal)


def test_pinned_mirrors_equal_the_one_copy():
    assert matter_system.KAPPA == mc.KAPPA and valencia.KAPPA == mc.KAPPA
    assert inspect.signature(matter_rhs.species_rhs).parameters["beta_dm"].default == mc.BETA_DM
    defaults = [param.default for _, fn in inspect.getmembers(matter_system, inspect.isfunction)
                for name, param in inspect.signature(fn).parameters.items() if name == "beta_dm"
                and param.default is not inspect.Parameter.empty]
    assert defaults and all(d == mc.BETA_DM for d in defaults)


def test_cosmos_potential_agrees_between_lanes():
    p = cosmos.CosmosParams()
    for phi in (-0.7, 0.0, 0.4, 2.302585092994046):
        assert scalar_system.cosmos_potential(phi) == p.V0 * (math.exp(-p.lam * phi) - p.D)
