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

The dark-matter replacement (Nov 16 2026; evidence results/sd-dark-matter-2026-11-16/)
------------------------------------------------------------------------------------
The COSMOS lane's imported, fully coupled dust is replaced by S-D matter in two parts: the S part (beta_S = 0,
uncoupled) and the D part (drained D_r plus D's share in the stones), which carries beta = beta_D (GEAR-103).

- The split's FORM is DERIVED (GEAR-91): with beta_S = 0 the S-D sector's only phi-dependent term is
  -e^(-2 beta_D phi) D^2/2, so S-D matter's source in phi's equation is beta_D times f_D rho_SD with
  f_D = e^(-2 beta_D phi) (D_s^2 - <D^2>) / rho_SD. Checks: a free D ring (D_r) has f = 1 (GEAR-103); S ringing in
  the sea has f = 3u/(1 + 3u) = 3/4 at u = e^(-2 beta_D phi) = 1 (S's mass^2 = 1 + 3u: 1 is S's own, uncoupled;
  3u comes from D's sea: D's share in the stone); the held step has f = 4.
- SD_MATTER_COUPLED_FRACTION = 0.76 is READ from the stored cells (five cells and two fine-cadence cells, x_c 1.5 and
  2.0, r_max 16-36, both closures; whole cell and cores r < 6, 8; t = 13-30 and 30-54): the mean of the averaged series
  is 0.756-0.771 in every case; the averaged series itself ranges over SD_MATTER_COUPLED_FRACTION_BOUNDS. Derived
  anchor: 3/4 (1 - r) + r with D_r's share r = 0.06-0.08 (AMENDMENT-1) gives 0.765-0.770.
- The two parts are each held at a fixed coupling. The stone's share 3u/(1 + 3u) drifts with phi (0.753 today, 0.786
  at phi = 2.5): an approximation, recorded (Open).
- The QUANTITY is not closed: the total today is the old imported amount, matched through the borrowed "today". Only the
  KIND is replaced. Baryons and radiation stay imported. The production kernel still carries fully coupled dust
  (not changed here).
- Unpaid inputs of the Read value: the cells' lambda = 5000 and their H/m_S, one stone per cell, the cells' x_c, the
  averaging windows (D's ring period, S's period).
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

# The dark-matter replacement: S-D matter's coupled (D) fraction today. READ from the stored cells (form derived above);
# the bounds are the averaged series' full range over cells, radii and windows (results at both ends in the PR).
SD_MATTER_COUPLED_FRACTION = 0.76
SD_MATTER_COUPLED_FRACTION_BOUNDS = (0.672, 0.853)


def stone_coupled_share(phi):
    """D's share of a small S ring's mass in the sea (Derived, GEAR-91 with beta_S = 0): 3u / (1 + 3u)."""
    u = math.exp(-2.0 * BETA_D * phi)
    return 3.0 * u / (1.0 + 3.0 * u)
