# Checks behind the REV20 brief answers (Nov 10 2026)

`research/rev20_brief_checks_v1.py` -> `output.txt`. Research only; no production change. Unpaid inputs: lambda, V0,
phi_today (B and F are at production V0); beta_D = -0.04 in B. G_N = 1; G_reset is not computed here.

| Check | Result | Status |
|---|---|---|
| A | V0 is the origin of phi: (V0, D, phi) -> (V0 e^-c, D e^c, phi - c) leaves the lane identical to 10 digits (c = 1, -2, 30). The only invariant constant of phi's potential is the floor F = D V0 = 8.2425e-6. | derived (algebra) + run |
| B | If the floor is D's vacuum step (GEAR-91, D slaved to its well, beta_D = -0.04, matched today), turnaround comes at 1.739 e-folds (a_turn/a_today 5.69, H0 t_turn 4.26, phi_turn 2.785) instead of 1.988 (7.30, 5.42, 3.100). D's well depth is not periodic in the lane: x1.60 today -> turnaround, x3.5e5 to the proton clock-stop, x1.4e8 to the Planck endpoint. | run, at production V0 |
| C | In the w = 1 contraction S is exactly J0/Y0 of mu*tau. Last zero crossing at mu*tau = 0.894-2.405, i.e. \|H\| = 0.14-0.37 mu; the constant-H criterion 3\|H\| = 2 mu (0.667 mu) is 0.19-0.52 e-folds too deep. | derived + run |
| D | Field-space curvature of the GEAR-97 scalar-tensor form is K = -6 lambda/(2 + 3 lambda D^2)^2, never zero; the machine's is zero. No field redefinition connects them. K = 0 forces polar coordinates (radius = gravity's scalar), and then the potential must be U(rho) + rho^4 W(theta): the GEAR-67 quartic is exactly rho^4 (lambda/8)(1 + sin^2 2 theta); the mass terms -rho^2 cos(2 theta)/2 are not of that form. | derived (symbolic) |
| E | The uniform D mode is held for x = 4/3 to 2.45, unstable for 2.45 < x < ~50 (peak 0.62 at x = 10), and held again above x ~ 60. Single-harmonic S; unreliable at large x, where S is strongly anharmonic. | run |
| F | Under floor = step the D-holding pocket's density is 1.8-3.1 floors at the centre, 0.39-0.57 on average; phi's vacuum today is 7.36 floors. A static pocket binds only if rho > 2 rho_vac, i.e. phi > 1.4-1.7 (centre) or > 2.05-2.12 (mean), against phi_today = 0.179. Mean density / critical density today = 0.033-0.049 (the machine's Omega_b is 0.042). | estimate, at production V0 |
| G | In the exact scalar-tensor theory with the same potential, D inside matter is driven by matter's trace, not by S^2. Same time average, different pump: D's uniform mode grows at 0.77-1.35 for x = 4/3 to 2.45, where the machine holds it (rate 0). | derived + run (linear, single harmonic, no gravity) |

Also run for the audit: a one-solar-mass black hole fills the horizon (r_s H = 1) at 17.40 e-folds below today by a
pure w = 1 conversion (the archive's number) and at 15.86 below today (17.85 below turnaround) when integrated through
the lane, the convention used for the 30.6 and 45.6 beside it in GEAR-94.
