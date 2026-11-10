# GEAR-91 coupling, research setting (HANDOFF-to-projects Nov 10, Task 2)

`research/gear91_coupling_v1.py`, tests `tests/test_research_gear91_coupling_v1.py`. No production file changes:
`coupled(beta_S, beta_D, alpha0)` wraps the two production entry points the kernel uses
(`engine.production_kernel.scalar_rhs_arrays`, `engine.v55_matter.scalar_projection`) for the duration of a block.

    V_SD(S, D; phi) = e^(-2 beta_S phi) S^2/2 - e^(-2 beta_D phi) alpha0^2 D^2/2 + (lambda/8)(S^4 + 6 S^2 D^2 + D^4)

- Switches: beta_S (default **0**, S uncoupled), beta_D (default **beta = -0.04**); universal = both beta (a switch,
  not a default). alpha0 = 1. phi_ref = 0 (GEAR-78 convention; the +0.7% mass shift at phi = 0.179 is accepted).
- phi's equation gains **8 pi dV_SD/dphi** (phi at 8 pi G_N = 1, S/D at G_N = 1); S/D see the phi-dependent masses;
  the S/D stress carries the full V_SD; V_SD(0, 0; phi) = 0 (the saddle's vacuum is 0*, supplied by phi).

## Tests (all pass)

| Test | Result |
|---|---|
| (a) switches off (beta_S = beta_D = 0), one kernel step | bit-identical to production (all scalars and geometry) |
| (b) exchange, cell by cell, inhomogeneous state (S, D, phi, Pi, lapse non-uniform): S/D energy gain alpha sqrt(gamma) (d e/d phi) Pi from the **projection** vs phi's loss sqrt(gamma) Pi dPi/dt / 8 pi from the **kernel RHS** | balanced to 3e-9 of the transfer (finite-difference limit) for beta_D only, beta_S only, and universal; momentum likewise |
| (b) negative controls | missing 8 pi: residual 0.96 of the transfer; reversed sign: 2.0. Both fail as required |
| S/D RHS vs projection | dPS/dt and dPD/dt change by -alpha d(Delta V)/dS, d(Delta V)/dD of the projected Delta V |
| (c) universal switch | Omega_S^2 / Gamma_D^2 = 1/alpha0^2 for phi in [-5, 5] (alpha0 = 1 and 0.8) |
| potential identities | V_SD(0,0;phi) = 0 for all phi; even in D; analytic gradient = finite difference (S, D, phi) |

## Logged, not interpreted

`g_kappa(phi)`: with beta_S = 0, K_S = 1, K_D = e^(-2 beta_D phi); g = (1 + K_D)/2, kappa = (1 - K_D)/2, g/kappa.
For beta_D phi < 0 (production: beta_D = -0.04, phi > 0) kappa < 0. Its meaning for the saddle classification is not
derived (AMENDMENT-3-GEAR-91 §2): flagged for chat.

Unpaid inputs for any run using it: lambda; V0 and phi_today (runs with phi on are "at production V0"); the
beta_S / beta_D choice.
