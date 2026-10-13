# Origin-half-cell-only matched-time A/B — 13 October 2026

**Disposition:** diagnostic evidence only; no production admission.  
**Draft PR:** [#79](https://github.com/miahpenn/0-Cosmos-engine/pull/79)  
**Exact-head run:** [38343269247](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38343269247)  
**Artifact:** [11787518609](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38343269247/artifacts/11787518609)

## Exact test and provenance

- Current-main source baseline: \`c01f2467347e8aebd57d220374fd7bc4e93045bc\`.
- Baseline initializer blob: \`4cfbfe4be6ac3ce6ed492f4e435186997d94116a\`.
- Candidate changes only \`engine/v55_initial.py\`; production kernel blob remains \`2954e7a79656849b5f088053ea0f4ac97d4bb443\`.
- Vendor submodule: \`d6052d605673ce9d82cdc99b0df79855fa2ea215\`.
- No vendor/dependency changes and no other engine files differ.
- Full repository suite: **233 passed**. Source pin and candidate-scope gates passed.
- Both arms reached the same coordinate-time horizon \(t=0.12\), with radiation enabled and the default machine initialization, using \(r_{\max}=80\), CFL 0.03, and N=80/120/160 (4/6/8 steps).
- No production source, physical equation, CMC gauge, matter model, boundary, evolution operator, or campaign parameter is changed.

The candidate integrates only the origin half-cell with
\[
I_{1/2}=dr\left(\frac{401}{720}s_0-\frac{31}{480}s_1+\frac{11}{1440}s_2\right),
\]
where the three cell-centred samples are those of the same iterated initializer source. These coefficients exactly integrate the moments \(1,r^2,r^4\) for an even source. Every inter-centre interval remains the original trapezoid.

Artifact ZIP SHA-256: \`45a0cbc274ea646104f841a983d1fc1da4b1df86c2fbd191414d1d00563a2e93\`. The report JSON SHA-256 printed by the observer is \`7760364e8eef76847a3a1d5c39c0d732ed6914d97ea1321ee0bd123337f06116\`.

## Central/constraint result

| N | H(cell 0), t=0 baseline → origin-only | H(cell 0), t=.12 baseline → origin-only | Full max \(|H|\), t=.12 baseline → origin-only | Full max \(|M|\), t=.12 baseline → origin-only | Full max \(|C_\Lambda|\), t=.12 baseline → origin-only |
|---:|---:|---:|---:|---:|---:|
| 80 | 3.98859e-3 → 8.17222e-5 | 3.69186e-3 → 1.26977e-4 | 3.69186e-3 → 2.60191e-4 | 5.99673e-4 → 7.92436e-5 | 7.43726e-5 → 9.87486e-6 |
| 120 | 4.02811e-3 → 6.59952e-5 | 3.38685e-3 → 1.79491e-4 | 3.38685e-3 → 2.93906e-4 | 8.51211e-4 → 1.30423e-4 | 1.08323e-4 → 1.67199e-5 |
| 160 | 4.04156e-3 → 6.50608e-5 | 2.94708e-3 → 2.66657e-4 | 2.94708e-3 → 2.82793e-4 | 1.05399e-3 → 1.66432e-4 | 1.39283e-4 → 2.22571e-5 |

At \(t=0.12\), the full-grid Hamiltonian maximum is reduced by **90.4–93.0%**, momentum by **84.2–86.8%**, and connection by **84.0–86.7%**, across N=80/120/160. The central \(\dot A_a/r_0^2\) witness changes from \(+3.89\times10^{-3},+7.53\times10^{-3},+1.056\times10^{-2}\) to \(-2.06\times10^{-4},-3.77\times10^{-4},-1.78\times10^{-4}\). This indicates that the origin-half-cell correction suppresses most of the observed central source defect, but the signed value remains nonzero and changes sign; this is not a regularity/convergence proof. The even-lapse fit residual barely changes, as expected because the CMC row itself is not modified.

## Mass/current/interval ledger

The conservation check uses the production-defined \(dM_{\rm MS}/dt-(F_T+W_{pR})\) at the nearest coordinate cell to \(r=10\) (the cell is \(r=9.5\) at N=80 and \(r=9.75\) at N=120/160). It is **not** the separate fixed-areal-radius worldtube test.

| N | Current residual max, baseline → origin-only | Change | Interval mass-defect max, baseline → origin-only | Change |
|---:|---:|---:|---:|---:|
| 80 | 1.12259009e-6 → 1.12346111e-6 | +0.0776% | 2.90750137e-8 → 2.91033820e-8 | +0.0976% |
| 120 | 9.27015148e-8 → 9.27341979e-8 | +0.0353% | 1.04709548e-9 → 1.04783678e-9 | +0.0708% |
| 160 | 1.05742216e-7 → 1.05737064e-7 | -0.0049% | 9.96210487e-10 → 9.96139921e-10 | -0.0071% |

These differences are tiny compared with the existing residuals, unlike the 4-arm all-interval quadrature candidate which worsened current residuals by factors at some resolutions. Under this short matched horizon the origin-only change does not produce a material worsening of the coordinate-worldtube current or interval defect.

The candidate's absolute \(M_{\rm MS}(r\approx10)\) is lower already at \(t=0\): by \(-5.12320\times10^{-5}\), \(-1.53938\times10^{-5}\), and \(-6.51796\times10^{-6}\) at N=80/120/160. That is an initial-data change, not a time-integrated loss. Over \(0\le t\le0.12\), the \(M_{\rm MS}\) increments agree to about \(10^{-11}\) at each N, so the evolution increment and current budget are essentially unchanged over this window. A production-ready disposition still requires the fixed-areal-radius ledger and longer-horizon/continuum tests; neither is closed here.

Outer radiation witnesses track closely between arms through \(t=0.12\). This short result does not test the archived late-time radiation continuation failure.

## Conclusion and next gate

This is the narrowest useful attribution so far:

- The actual-profile centre-Ricci mismatch and the initial Hamiltonian peak largely disappear when only the origin half-cell integral is corrected.
- The remaining inter-centre trapezoid is untouched.
- The short-horizon coordinate-worldtube current residual and interval-integrated mass defect are effectively unchanged.
- The central \(\dot A_a/r_0^2\) is reduced but not proven convergent; momentum/connection later develop central errors.
- Fixed-areal-radius conservation, longer-time stability, matched-time spatial convergence, and the archived outer-radiation continuation gate remain open.

Keep PR #79 draft. The result supports testing this *specific origin-only numerical ingredient* further; it does not justify merging it yet or replacing the broader PR #18 candidate wholesale. Next, compare this run to the existing fixed-areal-radius A/B infrastructure and extend only if that can be done with the same pinned source, reporting both centre constraints and the full signed mass/flux/work ledger. Do not alter physical equations or fit quadrature weights to current residuals.
