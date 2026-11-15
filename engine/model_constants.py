"""The model constants shared by the COSMOS lane and the local lane: one copy.

Before this module the same numbers were written twice, in ``engine.cosmos.CosmosParams`` and in
``engine.scalar_system`` / ``engine.production_kernel``, and nothing kept the copies equal. The lane conversion
(results/kappa-channel-2026-11-14/KAPPA_CONVERSION_DECISION.md) holds only if they are equal.

Values unchanged; production is bit-identical. Labels (docs/GEAR_TO_CODE_MAP.md): V0 and BETA_DM are inputs not
derived from the gears (BETA_DM borrowed); D_POT is the hand floor; LAMBDA_PHI the potential's slope.

KAPPA is 8 pi G with G = c = 1: the factor between a physical density and the Einstein-normalised density the
COSMOS lane carries (rho_E = KAPPA * rho).
"""
import math

KAPPA = 8.0 * math.pi
V0 = 8.242522415500654e-5
D_POT = 0.10
LAMBDA_PHI = 1.0
BETA_DM = -0.04
