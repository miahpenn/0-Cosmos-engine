"""The model constants shared by the COSMOS lane and the local lane: one copy.

Before this module the same numbers were written twice, in ``engine.cosmos.CosmosParams`` and in
``engine.scalar_system`` / ``engine.production_kernel``, and nothing kept the copies equal. The lane conversion
(results/kappa-channel-2026-11-14/KAPPA_CONVERSION_DECISION.md) holds only if they are equal.

Labels (docs/GEAR_TO_CODE_MAP.md): V0 and BETA_DM are inputs not derived from the gears (BETA_DM borrowed);
LAMBDA_PHI is the potential's slope.

KAPPA is 8 pi G with G = c = 1: the factor between a physical density and the Einstein-normalised density the
COSMOS lane carries (rho_E = KAPPA * rho).

Hand-floor removal, step 2 (Nov 16 2026; step 1 research: PRs #499-#502, results/floor-swap-2026-11-15/)
--------------------------------------------------------------------------------------------------------
The hand floor (V = V0 (e^-phi - D_POT), D_POT = 0.10 put in by hand) is replaced by the floor the S-D sector makes
itself: D at rest in the sea, D_s^2 = 2 e^(-2 beta_D phi) / lambda, has the vacuum energy

    V_SD(0, D_s; phi) = -e^(-4 beta_D phi) / (2 lambda)                (G = 1; Einstein-normalised: times KAPPA)

with the GEAR-91 coupling of S and D to phi,

    V_SD(S, D; phi) = e^(-2 beta_S phi) S^2/2 - e^(-2 beta_D phi) D^2/2 + (lambda/8)(S^4 + 6 S^2 D^2 + D^4),

and phi's equation gains 8 pi dV_SD/dphi (phi at 8 pi G = 1; S, D at G = 1).

Labels carried with these constants (and in the PR, the code and the tests):

- D_POT = 0: the hand floor is off. Its old value is kept only as HAND_FLOOR_REFERENCE (the labelled reference path)
  and as FLOOR_FRACTION_TODAY, the ruler the new lambda is matched against.
- Floor: MECHANISM DERIVED (D's step, felt through phi by GEAR-91); SIZE SET THROUGH lambda.
- SD_QUARTIC_LAMBDA = 1.57e6 is A RULER, not a derived number: it is matched so that the sea's floor today equals the
  old hand floor, KAPPA e^(-4 beta_D phi_today) / (2 lambda) = FLOOR_FRACTION_TODAY * V0, through the BORROWED "today"
  (V0 and phi_today, GEAR-92). C6's reset criterion gives lambda* ~ 10-1000 (conditional on C6); the gap of three to
  five orders of magnitude is THE RELOCATED HIERARCHY, recorded here, not hidden. lambda is computed below from the
  constants, not written as digits.
- BETA_D: FORM DERIVED (GEAR-91, confirmed by 0*, witnessed externally); VALUE NOT DERIVED (it is the locked dark-matter
  coupling, BETA_DM, borrowed). BETA_S = 0: S is tangible matter (AMENDMENT-3-GEAR-91: a coupled S meets the Cassini
  bound ~100x over); a switch, not derived.
- All times read from the lane are CHART-DEPENDENT through PHI_TODAY (borrowed).
- Held-step readings with the coupling on are made against D's phi-dependent rest value (engine.scalar_system.
  sd_rest_value), not against a fixed D*.
- Unpaid inputs: H/m_S (the one unpaid number under the one-budget rule); V0 and phi_today (borrowed "today");
  FLOOR_FRACTION_TODAY = 0.10 (forced: the old hand floor's size, kept as the matching ruler); BETA_D's value.
"""
import math

KAPPA = 8.0 * math.pi
V0 = 8.242522415500654e-5
LAMBDA_PHI = 1.0
BETA_DM = -0.04

# Borrowed chart: phi today (the archive-reconstructed present-day state, GEAR-92). Same value as
# engine.cosmos.construct_present_day and engine.v55_initial.
PHI_TODAY = 0.179055

# The hand floor is off in production (step 2). The old value survives only as the labelled reference path and as the
# ruler that lambda is matched against.
D_POT = 0.0
HAND_FLOOR_REFERENCE = 0.10
FLOOR_FRACTION_TODAY = HAND_FLOOR_REFERENCE   # forced: the floor's size today, kept as a ruler (not derived)

# GEAR-91 couplings of S and D to phi. Form derived; BETA_D's value is the borrowed dark-matter coupling.
BETA_S = 0.0
BETA_D = BETA_DM

# GEAR-67 on-site quartic, matched TODAY (a ruler, ~1.5689e6): the sea's floor equals FLOOR_FRACTION_TODAY * V0 at
# PHI_TODAY, so today's reconstructed budget is unchanged. C6: lambda* ~ 10-1000 (the relocated hierarchy).
SD_QUARTIC_LAMBDA = (KAPPA / (2.0 * FLOOR_FRACTION_TODAY * V0)) * math.exp(-4.0 * BETA_D * PHI_TODAY)
