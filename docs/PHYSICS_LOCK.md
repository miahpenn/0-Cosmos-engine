# Physics Lock

## Local GEAR-03
With kappa=0 and g=1:
- omega_S^2 = kappa + g = 1
- alpha_D^2 = g - kappa = 1
- V_S = +S^2/2
- V_D = -D^2/2

## GEAR-63/65 (retired to provenance, Nov 11 2026)
RETIRED: by DECISION-map-row-5 (Jeremiah, Nov 11 2026; `docs/GEAR_TO_CODE_MAP.md` row 5) the GEAR-67 S-D quartic is the
one arrest closure. The cubic barrier and the archive damping below are provenance only; `engine/local.py` keeps them,
unchanged in behaviour, for archive comparisons and the reduced control runner; no production kernel path imports it.
The cubic force is that of the quartic (pi/6) G alpha_D^2 D^4, i.e. GEAR-67's D^4 term at lambda = (4 pi/3) G alpha_D^2
= 4.19 (D* = sqrt(3/(2 pi)) = 0.690988), with no S-D cross term: the only place the archive fixed lambda from
something deeper, kept on record should D -> matter conversion reopen.

The reduced local sector (`engine/local.py`) retains the cubic geometric barrier -(2 pi/3) G alpha_D^2 D^3.

Provenance: the coefficient is derived from a closure, not fitted. It is the conformal force -(R/6) D with the curvature sourced by a pressureless dust remnant of density (1/2) alpha_D^2 D^2, R = 8 pi G (1/2 alpha_D^2 D^2) = 4 pi G alpha_D^2 D^2 (the GEAR-67E converted/dust branch).

Status: conditional on D -> matter conversion, which is open. With D's own stress-energy as the source, the same conformal coupling gives R = -8 pi G alpha_D^2 D^2 and a force of +(4 pi/3) G alpha_D^2 D^3, which is anti-restoring (GEAR-63, Sep 30 rewrite, section 3(b); `docs/CONDITIONAL_CANONICAL_D_TRACE_AUDIT_2026-10-15.md`). With no curvature coupling there is no cubic force; a restoring cubic from the mode's own curvature needs xi < 0, a new parameter. The code applies the barrier at all times, i.e. the post-conversion force also acts on the pre-conversion mode; this is part of the conditional model.

Not adopted: the Clock-Bridge coefficient Lambda_u(R) (Oct 3 report). It scales as Lambda_u = 0.039823 R^2, about 1.0008 R^2/(8 pi), as dimensional scaling of the cavity projection requires, so it replaces the closure with the hand-placed cavity radius R.

The archive damping eta_G D^2 Ddot (`research/gear_archive_machine_v1.py`, `archive_docx` mode only) is a separate imposed term with no derived origin. `engine/local.py` declares `eta_G` but does not use it.

Recorded in GEAR-73 and its amendment (Nov 8, 2026).

## COSMOS
The corrected archive operating point uses the reconstructed potential normalization rather than the legacy hard-coded V0=1 control. Archived corrected V0 is approximately 8.2425e-5.

## Coupling rules
Only quantities already present in the machine may cross the interface:
1. local stress-energy;
2. local geometric flux;
3. local energy/transport;
4. COSMOS component/source-energy ledger quantities.

No phenomenological efficiency, EOS, pressure closure, time-scale conversion, bounce law, or matter-identification law may be introduced merely to make the coupled system close.

## Promotion rule
A candidate coupling must satisfy algebraic consistency, donor-plus-recipient conservation accounting, radial and temporal resolution checks, constraint monitoring, and no hidden source term.
