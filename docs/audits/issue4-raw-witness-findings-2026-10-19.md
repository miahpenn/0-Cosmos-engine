# Issue #4 — accepted-state raw witness findings (2026-10-19)

## Scope and provenance

This is a diagnostic-only interpretation of the uploaded artifact from [workflow run #39187860872](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/39187860872), job `exact-replay`, PR #408. The workflow completed successfully and uploaded artifact `accepted-state-raw-momentum-witness`.

- Artifact ZIP SHA-256: `2de2375d53b74b918d61ba9dc9e06f7f2324aa313889059c60b7455477d0bb35`
- Diagnostic branch head at run: `21e39211f91693107dc684ffc76d9176fb063ea8`
- Historical production source: `7f47669a1c0edfafef9fc5f27b2db80840d3f035`
- Vendor pin: `d6052d605673ce9d82cdc99b0df79855fa2ea215`
- Initial-state digest: `a256b57f274ee120b55c0f805ed488d91382a9799647c8f9811d7df82aa337c0`
- Accepted-state digest: `c6acb3711d59a85c804a88208c9bc0356e8227a806c0684ae3a8afac72c9091a`
- Replay configuration: resolution 60, (r_{\max}=80), (\Delta t=0.015), one production step.
- Raw witness: 131 finite arrays. The observer did not mutate the accepted state.

A green workflow means the diagnostic scripts executed and their provenance assertions passed. It does **not** mean every diagnostic acceptance gate passed.

## Result 1 — current identity

The tested candidate identity is
`project_species(metric, conserved, species)["j"] = conserved.momentum_r / (KAPPA * sqrt_gamma)`.

The script deliberately uses a strict roundoff-only bound and does not soften it for radiation.

| State | Species | Max absolute difference | Max current magnitude | Relative difference | Strict roundoff check |
|---|---|---:|---:|---:|---|
| Accepted | Dark matter | (3.2312\times10^{-27}) | (1.8410\times10^{-11}) | (1.76\times10^{-16}) | Pass |
| Accepted | Baryons | (4.0390\times10^{-28}) | (2.8697\times10^{-12}) | (1.41\times10^{-16}) | Pass |
| Accepted | Radiation | (4.7395\times10^{-27}) | (8.4903\times10^{-15}) | (5.58\times10^{-13}) | **Fail** |
| Manufactured | Dark matter | (6.6174\times10^{-24}) | (2.4049\times10^{-8}) | (2.75\times10^{-16}) | Pass |
| Manufactured | Baryons | (1.6544\times10^{-24}) | (3.9528\times10^{-9}) | (4.19\times10^{-16}) | Pass |
| Manufactured | Radiation | (3.8128\times10^{-24}) | (5.1069\times10^{-12}) | (7.47\times10^{-13}) | **Fail** |

The radiation relative discrepancy is approximately (5.6\times10^{-13}) on the accepted state and (7.5\times10^{-13}) on the manufactured state. The pinned radiation primitive recovery uses an iterative solve with relative residual tolerance (10^{-12}), making recovery error a plausible explanation, but this is not yet proven as the sole cause. Do not relax the bound without a separately justified recovery-aware error budget and a controlled test.

The check is algebraic only. It does not differentiate the projection along the RHS or establish the momentum-constraint propagation identity.

## Result 2 — projection directional derivative

The projection diagnostic ran four central-difference scales: (8\times10^{-5}, 4\times10^{-5}, 2\times10^{-5}, 10^{-5}). All reported arrays were finite and the accepted-state digest remained unchanged.

Relative max-norm difference between the two finest scales:

| Direction | Dark matter | Baryons | Radiation |
|---|---:|---:|---:|
| Conservative-state only | (3.95\times10^{-13}) | (3.17\times10^{-13}) | (5.12\times10^{-10}) |
| Metric only | (2.04\times10^{-9}) | (1.23\times10^{-9}) | (1.35\times10^{-9}) |
| Full RHS with lapse resolve | (3.29\times10^{-13}) | (3.17\times10^{-13}) | (5.12\times10^{-10}) |

These numbers are reported convergence indicators, **not a passing plateau**: the diagnostic explicitly has `projection_step_size_plateau_established=false`, because a numerical acceptance bound has not been predeclared. Radiation's conservative/full-direction differences are visibly larger than the dust-species differences and warrant a recovery-aware derivative error analysis.

## Required next steps

1. Preserve the current strict identity result as a failure for radiation; do not tune thresholds to force acceptance.
2. Separate radiation primitive-recovery error from projection algebra using an independently verified manufactured-state recovery round trip. Record primitive residuals and current error together, over a declared tolerance/conditioning sweep.
3. Predeclare derivative plateau criteria (absolute and relative scales, including a roundoff floor) before interpreting the epsilon sequence. Retain raw signed arrays.
4. Once the projection derivative is independently controlled, assemble the matter-current contribution with independently computed geometric and scalar-current derivatives. Compare the signed momentum-constraint tangent at cells 0/1 and fixed-coordinate annuli.
5. Keep all gates unchanged until evidence satisfies them: `projection_step_size_plateau_established=false`, `independent_full_propagation_identity_closed=false`, `scientific_admission=false`, and `acceptance_oracle_satisfied=false`.

No production equations, source pins, or scientific-admission rules are changed by this note.
