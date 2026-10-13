# Origin-half-cell-only matched-time A/B — 13 October 2026

**Disposition:** diagnostic evidence only; no production admission.  
**Draft PR:** [#79](https://github.com/miahpenn/0-Cosmos-engine/pull/79)  
**Latest exact-head run (fixed-R ledger included):** [38343526931](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38343526931)  
**Latest artifact:** [11787991716](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38343526931/artifacts/11787991716)  
**Initial coordinate-ledger-only run:** [38343269247](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38343269247)

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

The conservation check uses the production-defined \(dM_{\rm MS}/dt-(F_T+W_{pR})\) at the nearest coordinate cell to \(r=10\): the selected coordinates are \(r=9.5\) at N=80, \(r=10.333333333333332\) at N=120, and \(r=9.75\) at N=160. At N=120, \(r=10\) is exactly midway between two cell centres; the existing floating-point `argmin` tie-break selects the outer centre. Baseline and candidate use the same selected cell at each N, but this is **not** a matched fixed-coordinate-radius or fixed-areal-radius worldtube test.

| N | Current residual max, baseline → origin-only | Change | Interval mass-defect max, baseline → origin-only | Change |
|---:|---:|---:|---:|---:|
| 80 | 1.12259009e-6 → 1.12346111e-6 | +0.0776% | 2.90750137e-8 → 2.91033820e-8 | +0.0976% |
| 120 | 9.27015148e-8 → 9.27341979e-8 | +0.0353% | 1.04709548e-9 → 1.04783678e-9 | +0.0708% |
| 160 | 1.05742216e-7 → 1.05737064e-7 | -0.0049% | 9.96210487e-10 → 9.96139921e-10 | -0.0071% |

These differences are tiny compared with the existing residuals, unlike the 4-arm all-interval quadrature candidate which worsened current residuals by factors at some resolutions. Under this short matched horizon the origin-only change does not produce a material worsening of the coordinate-worldtube current or interval defect.

The candidate's absolute \(M_{\rm MS}(r\approx10)\) is lower already at \(t=0\): by \(-5.12320\times10^{-5}\), \(-1.53938\times10^{-5}\), and \(-6.51796\times10^{-6}\) at N=80/120/160. That is an initial-data change, not a time-integrated loss. Over \(0\le t\le0.12\), the \(M_{\rm MS}\) increments agree to about \(10^{-11}\) at each N, so the evolution increment and current budget are essentially unchanged over this window. The later exact-head run also includes the separate moving fixed-areal-radius ledger; results are below. Longer-horizon/continuum tests remain open.

Outer radiation witnesses track closely between arms through \(t=0.12\). This short result does not test the archived late-time radiation continuation failure.


## Fixed-areal-radius \(R=10\) ledger — exact-head follow-up

The first run sampled the nearest coordinate cell to \(r=10\), which does not follow a physical-radius tube. The follow-up exact-head run [38343526931](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38343526931), artifact [11787991716](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38343526931/artifacts/11787991716), preserves the same source pins and adds a separate \(R=10\) ledger at every sample over the same matched \(t=0.12\) runs. The full repository regression suite passed again (233 tests); both baseline and candidate reached the horizon at every resolution.

At fixed areal radius the tube's coordinate location moves to approximately \(r=9.99317\) by \(t=0.12\). The predicted rate includes the sweep contribution
\[
\left.\frac{dM_{\rm MS}}{dt}\right|_{R=10}
\approx F_T+W_{pR}+v_r\,M_{,r}^{\rm stress},
\qquad
v_r=-\frac{\dot R}{R_{,r}},
\]
and compares that against the finite-time derivative of the interpolated \(M_{\rm MS}(R=10,t)\). The independent radial-gradient identity is also measured as \(D_r M_{\rm MS}-4\pi R^2(\rho R_{,r}+j\,\dot R/\alpha)\), using the repository's regular odd-parity \(R_{,r}\) operator. This is still a short finite-difference diagnostic, not an instantaneous conservation proof.

| N | Fixed-\(R\) current residual max, baseline → origin-only | Change | Fixed-\(R\) interval defect max, baseline → origin-only | Change | \(|\mathrm{radial\ identity\ error}|\) at final sample, baseline → origin-only |
|---:|---:|---:|---:|---:|---:|
| 80 | 4.24933594e-7 → 4.25639337e-7 | +0.1661% | 1.15505752e-8 → 1.15735285e-8 | +0.1987% | 2.74473654e-5 → 2.74474261e-5 |
| 120 | 3.45276112e-7 → 3.45256696e-7 | -0.0056% | 6.29092999e-9 → 6.29049948e-9 | -0.0068% | 1.17045455e-5 → 1.17045507e-5 |
| 160 | 2.53691396e-7 → 2.53688119e-7 | -0.0013% | 3.58901279e-9 → 3.58896532e-9 | -0.0013% | 6.49281626e-6 → 6.49281729e-6 |

The fixed-\(R\) ledger therefore tells the same bounded story as the coordinate-cell ledger: the origin-only candidate produces no material change in outer current/interval balance through \(t=0.12\). The residual changes by less than 0.2% at every resolution. The radial mass-gradient identity at \(R=10\) is also unchanged to \(6.1\times10^{-11}\), \(5.2\times10^{-12}\), and \(1.1\times10^{-12}\) absolute at N=80/120/160. In contrast, the broader all-interval quadrature candidate in PR #18 had larger fixed-\(R\) residual and interval defects at every tested resolution. That contrast supports keeping the origin-only counterfactual separate from PR #18's full radial integration rule.

The fixed-\(R\) residual itself is not zero, and the radial mass-gradient identity error remains around \(2.7\times10^{-5}\), \(1.2\times10^{-5}\), and \(6.5\times10^{-6}\) at the sampled tube. Thus this diagnostic does **not** close the conservation/admission gate. It establishes only that, in a short source-pinned A/B, the central improvement from the origin-only term is not accompanied by the outer fixed-\(R\) balance penalty previously observed for the broader quadrature candidate.

**Updated disposition:** retain draft/unmerged PR #79 as a narrow numerical candidate. Next gates are centre/constraint convergence under matched resolution, more temporal samples / longer bounded horizon, and the existing archived radiation continuation test. Do not promote the broader PR #18 quadrature by association and do not change physics equations or fit weights to balance residuals.

## Higher-resolution matched-time check — \(N=160,240,320\)

The higher-resolution sibling run [38343781776](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38343781776) completed successfully on the same production-source baseline, initializer candidate, and vendor pins. It passed the complete repository suite (**233 passed**) and all six baseline/candidate arms reached \(t=0.12\) with radiation enabled. Artifact [11788575091](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38343781776/artifacts/11788575091), artifact ZIP SHA-256 \`96d03a2f09abb54ca4ed06b852adbabf13952f1964600adaf73dd05697e6d671\`, report JSON SHA-256 \`b678d42b40c234b595206adfe171de8a898538b9891314e0fa9516387e8a7626\`.

| N | Final full-grid \(|H|\) max, baseline → origin-only | Final full-grid \(|M|\) max, baseline → origin-only | Final full-grid \(|C_\Lambda|\) max, baseline → origin-only | Final \( \dot A_a/r_0^2\), baseline → origin-only |
|---:|---:|---:|---:|---:|
| 160 | 2.94708e-3 → 2.82793e-4 | 1.05399e-3 → 1.66432e-4 | 1.39283e-4 → 2.22571e-5 | +1.05604e-2 → -1.78030e-4 |
| 240 | 1.83935e-3 → 4.85100e-4 | 1.26318e-3 → 1.98920e-4 | 1.88580e-4 → 3.04936e-5 | +7.62152e-3 → +2.61484e-3 |
| 320 | 6.80833e-4 → 7.18254e-4 | 1.16397e-3 → 1.77971e-4 | 2.16770e-4 → 3.48008e-5 | -1.83637e-2 → +1.02814e-2 |

The candidate continues to improve the final momentum and connection maxima by about 84% at N=240/320. However, its Hamiltonian advantage narrows with resolution and reverses at N=320: the candidate final full-grid H maximum is **5.50% higher** than baseline. The central \( \dot A_a/r_0^2\) witness is also not convergent: it changes sign and grows in magnitude across the candidate's refinement sequence (N=160/240/320: \(-1.78\times10^{-4},+2.61\times10^{-3},+1.03\times10^{-2}\)). Initial \(H(r_0)\) remains greatly reduced (\(\approx6.5\text{–}6.7\times10^{-5}\), versus \(\approx4.05\times10^{-3}\) baseline), so the origin-only integration changes the initial constraint but does not by itself produce a regular, stable central evolution.

The fixed-\(R=10\) balance remains nearly identical between arms over the short horizon: candidate-minus-baseline changes in current-residual max are \(-0.00129\%\), \(-0.00510\%\), and \(+0.0128\%\) at N=160/240/320; interval-defect max changes are \(-0.00132\%\), \(-0.00314\%\), and \(+0.00483\%\). This again separates the origin-only candidate from the conservation deterioration of the broader all-interval rule. It does not make the finite residual vanish or certify conservation.

**High-resolution disposition:** the origin-only half-cell correction is a verified way to suppress the *initial* centre constraint/Ricci anomaly, and it does not materially worsen the short-horizon fixed-R ledger. But the residual central RHS and evolved Hamiltonian do not show reliable refinement behavior. Keep it diagnostic-only: the result does not support merging this initializer patch, and the next engineering layer to diagnose is the regular-centre discretization/RHS propagation after initialization, not further adjustment of the quadrature coefficients. No physical equation or coefficient should be changed from these data alone.


## Conclusion and next gate

This is the narrowest useful attribution so far:

- The actual-profile centre-Ricci mismatch and the initial Hamiltonian peak largely disappear when only the origin half-cell integral is corrected.
- The remaining inter-centre trapezoid is untouched.
- The short-horizon coordinate-worldtube current residual and interval-integrated mass defect are effectively unchanged.
- The central \(\dot A_a/r_0^2\) is reduced on the first slices but is explicitly non-convergent after evolution at N=160/240/320; the Hamiltonian advantage reverses at N=320 even though momentum/connection maxima remain lower.
- The fixed-areal-radius ledger has now been compared and is effectively unchanged in the origin-only arm over this short horizon; longer-time stability, matched-time spatial/temporal convergence, and archived outer-radiation continuation remain open.

Keep PR #79 draft. The result supports testing this *specific origin-only numerical ingredient* further; it does not justify merging it yet or replacing the broader PR #18 candidate wholesale. Next, reconcile the nearly unchanged fixed-R ledger against the existing PR #35 component-resolved current analysis, then choose a short matched-resolution/time refinement that records both centre constraints and the full signed mass/flux/work ledger. Do not alter physical equations or fit quadrature weights to current residuals.
