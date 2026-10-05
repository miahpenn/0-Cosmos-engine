# V5.5 production architecture locks

The production spherical branch is one spacetime and one total stress-energy tensor.

### Matter
Local S:
    V_S = +S^2/2

Local inverted D:
    V_D = -D^2/2

COSMOS scalar:
    V(phi) = V0 [exp(-phi) - 0.10]
with corrected V0 = 8.242522415500654e-5.

Dark matter:
    nabla_mu T_DM^{mu nu} = + beta rho_DM nabla^nu phi
    beta = -0.04

Scalar receives the opposite exchange.

Baryons are pressureless dust.

Radiation is a perfect fluid with:
    p_r = rho_r/3

The archive's Bianchi-I shear has:
    p_shear = rho_shear
and is excluded from the exact spherical source until the geometry is lifted.

### Geometry
Reference-metric spherical BSSN and sequential PIRK2 are supplied by the pinned vendor/bb-palatini-unified-r0 commit documented in docs/PRODUCTION_KERNEL.md.

The external matter equations are not imported.

### Coupling
There is no shell equation of state, phenomenological interface source, fitted feedback coefficient, D-to-matter identification law, clock-conversion formula, bounce rule, reset, branch flip, lapse floor, ejection threshold, or physical stop condition.

The former interface is a set of derived observables:
stress-energy, geometric flux, Misner-Sharp mass/current, proper time, and cycle events.

### Numerical boundaries
A finite campaign horizon is a host-computation boundary, not physics.
Non-finite or nonpositive evolution variables constitute numerical failure,
not a physical interpretation.
