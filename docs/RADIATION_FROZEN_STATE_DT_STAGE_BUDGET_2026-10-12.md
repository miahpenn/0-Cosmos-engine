# Frozen-State Radiation Timestep Budget — Signed Stage Evidence

**Date:** 12 October 2026 UTC  
**Classification:** Diagnostic-only, one-step A/B from a checksum-pinned accepted state.  
**Branch:** `diag/radiation-frozen-state-stage-budget-2026-10-12`  
**Workflow:** [38324794759](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38324794759) — success.  
**Artifact:** [11779854674 — radiation-frozen-state-dt-stage-budget](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38324794759/artifacts/11779854674)  
**Artifact ZIP SHA-256:** `94bdbafbf3b988a85be84bc6779d0a74fdd7d8192679beb820858daf99c86e28`  
**Embedded report SHA-256:** `6983779ab1f5a0bf85e03c439f828c00e6e02d64515762db709a5f1751201d34`  
**Instrumented branch head:** `c355a72a2ac0d147dbe2ec29e8db7ae54edaacdc`  
**Frozen source baseline:** `cdcec8282bff0443f2fe4674b15929822ff38736`  
**Pinned numerical vendor:** `d6052d605673ce9d82cdc99b0df79855fa2ea215`  
**Accepted-state NPZ SHA-256:** `ccd1fa0fbabe65c51c3dd7ac2c9227d6d25aa680c169eef1afd3fca4326294cf`  
**Accepted-state metadata SHA-256:** `0e6251ec58b843d5d5f0dfc8911afea0dbc5da27b409f34d3435628c61676858`

## Question and method

Does the observed radiation admissibility rejection at \(t=47.16\) persist when one production step is taken with smaller \(\Delta t\), starting from exactly the same saved accepted state at \(t=47.13\)?

The audit downloads the accepted full numerical state from the already completed [boundary/CMC audit run 38322337250](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38322337250). It verifies the frozen NPZ SHA-256 and source/vendor metadata, reconstructs the numerical arrays, and performs independent production-kernel step attempts at \(\Delta t=0.03,\ 0.015,\ 0.0075,\ 0.00375\). It does **not** replay the 1,571-step lead-in. Each arm begins from a deep copy of the same frozen arrays; all physics and the exact \(E\geq|S|\) guard remain unchanged.

Workflow gates: pinned source/vendor scope check passed; 1 targeted metric identity test passed; **91 repository-wide tests passed**; all four arms were recorded; no missing radiation source cells or instrumentation errors. This success validates the diagnostic recorder, not a completed multi-step trajectory.

## Decision table

| \(\Delta t\) | Real production predictor | Maximum \(|S|/E\) in predictor | Minimum \(C\) over full grid |
|---:|---|---:|---:|
| 0.03000 | Rejected at explicit predictor, cell 79 | 1.0004012501659771 | \(-9.96703374185\times10^{-8}\) |
| 0.01500 | Completed one step | 0.9996133976881819 | \(+2.03884377471\times10^{-9}\) |
| 0.00750 | Completed one step | 0.9992315204412235 | \(+2.03846434337\times10^{-9}\) |
| 0.00375 | Completed one step | 0.9990436008271782 | \(+2.03827460084\times10^{-9}\) |

Here \(C=U_E-\sqrt{\gamma^{rr}}|U_r|\). The full-grid minimum is at the regular-centre index \(i=0\) in the successful smaller-step arms; the near-saturated outer-cell budget is separately given below. Those are different witnesses and should not be conflated.

## Signed outer-cell budget at \(i=79,\ r=79.5\)

| \(C\) witness | \(\Delta t=0.03\), rejected | \(\Delta t=0.015\) | \(\Delta t=0.0075\) | \(\Delta t=0.00375\) |
|---|---:|---:|---:|---:|
| Accepted metric | \(+2.80644382492\times10^{-7}\) | same | same | same |
| After flux, accepted metric | \(+1.04420641674\times10^{-6}\) | \(+6.51858239046\times10^{-7}\) | \(+4.68897761730\times10^{-7}\) | \(+3.75426706844\times10^{-7}\) |
| After geometric source, accepted metric | \(+1.68264528290\times10^{-6}\) | \(+9.81644832695\times10^{-7}\) | \(+6.31144607594\times10^{-7}\) | \(+4.55894495043\times10^{-7}\) |
| Same predicted conservative state on explicit predictor metric | \(-9.96703374185\times10^{-8}\) | \(+9.55066906386\times10^{-8}\) | \(+1.89324177500\times10^{-7}\) | \(+2.35295660932\times10^{-7}\) |

The actual predicted conservative state matches the reconstructed RHS state exactly at stored precision in all four arms: maximum \(U_E\) and \(U_r\) reconstruction errors are both zero, and the signed split closure error is zero. No term is added to rescue the cone; the order of accounting is fixed as accepted metric → flux → source → explicit predictor metric.

At the failing nominal \(\Delta t=0.03\), flux and source contributions add \(+7.63562034246\times10^{-7}\) and \(+6.38438866160\times10^{-7}\) to the outer-cell margin when evaluated on the accepted metric. The same conservative state then evaluated on the explicit predictor metric loses \(1.78231562032\times10^{-6}\), making \(C=-9.9670\times10^{-8}\). The stage budget therefore localizes the crossing to the **explicit metric conversion of the near-cone-boundary radiation state**, not the flux/source update alone.

## Metric identity closure

At \(i=79\), the predictor metric is \(\gamma^{rr}=X^2/a\). For the nominal step:

| Quantity | Value |
|---|---:|
| Accepted \(\gamma^{rr}\) | 1.5315453355337358 |
| After \(a\) update only | 1.5454793873641697 |
| After \(X\) update | 1.5537534421946575 |
| \(\Delta\gamma^{rr}\) from \(a\) | \(+0.013934051830433836\) |
| \(\Delta\gamma^{rr}\) from \(X\) | \(+0.00827405483048782\) |
| Regularity projection / other contribution | 0 |
| Decomposition closure error | 0 |

The implemented zero-shift identities \(\dot a=-2\alpha a A_a\) and \(\dot X=\alpha X K/3\) agree with the recorded explicit RHS to zero / \(1.39\times10^{-17}\) maximum residual, respectively. This shows arithmetic consistency between the selected RHS and predictor. It does not by itself establish whether the discrete splitting/timestep or a continuum/boundary-conditioned trajectory is responsible for the long-run sensitivity.

## Constraint and Misner–Sharp witnesses

The archived accepted state begins with full sourced witnesses:
- \(H_{\max}=0.084934723579255\)
- \(M_{\max}=0.018343263859204154\)
- connection maximum \(=0.05829447235243363\)
- determinant-constraint maximum \(=1.638927844904714\times10^{-5}\)
- \(M_{\rm MS}(R_\sigma)=4.501726366496262\), \(R_\sigma=10.494223796016815\)
- \(\chi_\sigma=0.1420563437564698\)
- \(H_{\rm eff}=-0.09091019278724598\)
- outer CMC residual maximum \(=2.36009193898\times10^{-4}\); last-interior residual \(=2.25537742968\times10^{-5}\).

For completed smaller-step arms, those full sourced diagnostics remain nonzero and change with the step:
- \(dt=0.015\): \(H_{\max}=0.08505320892648259,\ M_{\max}=0.01840765888000545,\ M_{\rm MS}(R_\sigma)=4.591636438124435,\ H_{\rm eff}=-0.09180154849558375\).
- \(dt=0.0075\): \(H_{\max}=0.08499383097452469,\ M_{\max}=0.018375423347974665,\ M_{\rm MS}(R_\sigma)=4.546439645940305,\ H_{\rm eff}=-0.09135439499189411\).
- \(dt=0.00375\): \(H_{\max}=0.08496424355808607,\ M_{\max}=0.01835933409732787,\ M_{\rm MS}(R_\sigma)=4.524022903773679,\ H_{\rm eff}=-0.0911319262400953\).

The raw geometry-only constraints at the rejected predictor are not sourced \(H/M\) and must not be presented as such: their largest values are different because they omit stress-energy subtraction and were computed on a rejected stage. Full sourced \(H/M\) and \(M_{\rm MS}\) for successful arms are evaluated on their admissible completed state. The data do not show production admission: constraint and CMC residuals remain nonzero and one-step evidence is not convergence or conservation closure over time.

## Supported conclusion and guardrails

**Supported:** nominal \(dt=0.03\) crosses the exact cone at the existing explicit-predictor guard, while the same accepted numerical arrays make admissible one-step predictors at \(dt/2,\ dt/4,\ dt/8\). The local signed budget closes and attributes the negative margin to the predictor metric conversion, with both \(a\) and \(X\) raising \(\gamma^{rr}\).

**Not supported:** this does not yet prove a time-step retry yields stable evolution over many subsequent steps, that the finite-radius boundary is correct/incorrect, or that all Einstein constraints and Misner–Sharp balance pass. No physics change is justified by this result alone.

Preserve the exact radiation admissibility condition. No floors, clamps, cone relaxation, damping, fitted coefficients, source changes, boundary changes, or forced bounce are introduced. The real adaptive-retry integration is documented separately on the candidate branch and is independently gated. Production `main`, the pinned vendor, and canonical GEAR archive were not changed by this timestep-budget audit.
