# Discrete-consistent initial state: multi-constraint RHS tangent (2026-10-12)

## Run and provenance

- Diagnostic workflow: [38290633971](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38290633971)
- Commit: `5f7fd471f4d3afe46320535de9438ed71a4dd524`
- Input paired initial-data campaign: [38284401548](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38284401548)
- Artifact: `discrete-initial-constraint-tangent`, ID `11766405316`
- Artifact ZIP SHA-256: `2011162ba822bea1b52d2b8c2099284f97f7a0ca5e76e2d2bdc1788c324a7186`
- Frozen production source: `6d8b3f2790644591326df57c55908863f994cb0c`
- Pinned vendor: `d6052d605673ce9d82cdc99b0df79855fa2ea215`
- Resolutions: N=60, 80, 120; (r_{max}=80); radiation enabled; (dt=0.03\,dr) for the one-step comparator.
- Six frozen t=0 checkpoints were tested (baseline/candidate × three resolutions). Each state has two symmetric RHS-direction perturbations, (epsilon=2\times10^{-5}) and (10^{-5}), with the CMC lapse re-solved. A separate copy takes one ordinary production step for a local comparison.
- No long trajectory was run and no production equation, physical coefficient, gauge, or boundary condition was changed. `main` remains untouched.

## Main result

**The discrete-consistent initial metric removes the initial Hamiltonian residual by construction, but it does not remove the momentum-constraint RHS tendency.** The actual one-step (Delta M/dt) closely matches the semi-discrete tangent at the first admitted cell for all tested resolutions.

| N | First admitted cell (r) | Baseline (dot M) | Candidate (dot M) | Candidate one-step (Delta M/dt) |
|---:|---:|---:|---:|---:|
| 60 | 3.333 | (+1.2113\times10^{-4}) | (-8.7553\times10^{-5}) | (-8.6967\times10^{-5}) |
| 80 | 2.500 | (+1.4873\times10^{-4}) | (-1.2354\times10^{-4}) | (-1.2295\times10^{-4}) |
| 120 | 1.667 | (+2.1261\times10^{-4}) | (-1.9063\times10^{-4}) | (-1.8993\times10^{-4}) |

The one-step candidate rate differs from the tangent by below 0.7% at N=60/80/120 in these cells. The candidate does not make (dot M) vanish; its sign flips and its magnitude remains resolution-dependent. This shows the mismatch is in the coupled semi-discrete RHS, not merely a finite-difference reporting artifact. It does not identify a single incorrect equation.

## What improves, and what remains

At t=0, the candidate's maximum H residual is (7.65\times10^{-16}), (1.29\times10^{-15}), and (3.16\times10^{-15}) at N=60/80/120, while M and (C_\Lambda) remain at roundoff. Candidate admitted-inner instantaneous H-residual tangents are about (7\times10^{-9}) for all three resolutions, much smaller than the baseline H tangent. This improvement is expected to be distinguished from time-step stage effects: the normal one-step H quotient is larger than that frozen-state tangent, so it should not be described as proof of H constraint propagation.

The frozen-state (C_\Lambda) tangent is at roundoff in all three candidate inputs. The one-step quotient is small (order (10^{-6})); the current dominant witness is still momentum rather than an independently faulty Lambda RHS.

## Signed momentum budget — first admitted cell

The residual is (M=M_{geom}-8\pi J). In the candidate t=0 states, the derivative of the source side is dominated by local S scalar current, but the sum of geometry terms fails to cancel it at the first admitted cell:

| N | (dM_{geom}/dt) | (d(-8\pi j_S)/dt) | Net (dot M) |
|---:|---:|---:|---:|
| 60 | (-3.2966\times10^{-4}) | (+2.4214\times10^{-4}) | (-8.7553\times10^{-5}) |
| 80 | (-3.4710\times10^{-4}) | (+2.2359\times10^{-4}) | (-1.2354\times10^{-4}) |
| 120 | (-3.6351\times10^{-4}) | (+1.7290\times10^{-4}) | (-1.9063\times10^{-4}) |

For the geometry derivative, the important signed terms are (D_r\dot A_a), (-\frac23 D_r\dot K), and the radial-b contribution. For example, at N=120/i=2 they are (-4.6093\times10^{-4}), (+1.5690\times10^{-4}), and (-5.9462\times10^{-5}); together with the negligible (6A_a\chi_r) term they produce the listed (dM_{geom}/dt). No one term has been demonstrated to be the unique defect. The local scalar-source rate barely changes between the baseline and candidate; most of the A/B change is in the geometric response.

## Fixed-radius check: localization, not a convergence claim

The first admitted cell moves inward as N increases, so its values are not a like-for-like spatial comparison. At fixed radii the results are non-monotone near r=2, while the candidate is close to geometric/source cancellation at r=3–5 for N=120:

| Radius | N=60 candidate (dot M) | N=80 candidate (dot M) | N=120 candidate (dot M) |
|---:|---:|---:|---:|
| 2 | (+1.6085\times10^{-4}) | (+4.0873\times10^{-5}) | (-7.1734\times10^{-5}) |
| 3 | (-2.5453\times10^{-5}) | (-4.5373\times10^{-5}) | (+1.2286\times10^{-7}) |
| 4 | (-3.2010\times10^{-5}) | (+1.6493\times10^{-5}) | (+3.2882\times10^{-7}) |
| 5 | (+1.7011\times10^{-5}) | (-1.5194\times10^{-7}) | (-1.7416\times10^{-8}) |

These data support a central spatial/derivative compatibility layer whose influence is much stronger near r=2 than at fixed r=3–5 on the finest tested grid. Do not fit a formal order to this short, sign-changing sequence, and do not equate the N=120 cancellation at selected radii with proof of central regularity or long-time constraint preservation.

## Closure checks

The signed decompositions close at roundoff: base H reconstruction max error (8.67\times10^{-19}), base M reconstruction (0), tangent M reconstruction at most (8.67\times10^{-19}), and tangent H reconstruction at most (4.37\times10^{-14}). The latter is still a diagnostic accounting closure, not a proof of the continuum Einstein-matter identity.

## Decision and next bounded test

- Do not promote the discrete-consistent initial-metric solve as a production repair: its initial H closure is built into the Newton solve and it does not eliminate the M tangent.
- Do not change the S source, A/K source coefficients, CMC lapse or Lambda equation based on this report alone.
- Next, compare the **signed (\dot A_a) and (\dot K) components themselves** (lapse Hessian, Ricci, anisotropic-stress source, nonlinear/extrinsic terms) against their radial series identities on the paired t=0 states, then repeat on the already saved t=2.16 states. Focus on i=0,1,2 and fixed r=2,3,4,5. Check exact parity/derivative expressions and term closure; no time evolution is needed for this step.
- In parallel, keep the actual production source/current continuity identities separate from the geometric budget; a current tangent by itself is not an independent covariant stress-energy conservation test.
- Only if one discrete RHS identity demonstrably fails should a single isolated numerical-layer repair be proposed and re-tested. No long campaign until that admission gate passes.
