# 0-Cosmos V5.5 — Current Continuation Handoff: Centre-Series RHS, Quadrature, and Constraint Propagation

**Updated:** 2026-10-12 (UTC)  
**Repository:** [miahpenn/0-Cosmos-engine](https://github.com/miahpenn/0-Cosmos-engine)  
**Production branch status:** \`main\` at \`809228c380d80af83858412557c78b6329ffa5c8\` when this handoff was assembled; no diagnostic result below has been merged into main.  
**Frozen production source used in the numerical audits:** \`6d8b3f2790644591326df57c55908863f994cb0c\`  
**Pinned geometry/PIRK vendor:** \`d6052d605673ce9d82cdc99b0df79855fa2ea215\`  
**Active quadrature/centre-series branch:** \`diag/discrete-consistent-initial-production-test-2026-10-12\`  
**Separate frozen-tangent report branch:** \`diag/constraint-tangent-discrete-initial-2026-10-12\`  
**Admission:** Execution through the t=12 horizon is proven for N=40/60/80; scientific production admission remains OPEN.

## 1. Executive result

The strongest new result is that **radial quadrature accounts for a large share of the initial Hamiltonian residual, but it does not solve the centre-series RHS problem or close the full constraint-propagation gate**.

A diagnostic-only initializer uses a fourth-order cubic Lagrange integral on each centre-to-centre radial interval, alongside the already tested r²-exact origin half-cell integral. The initial first-cell Hamiltonian residual falls from roughly \(4.05\times10^{-3}\) in the pinned initializer to \(4.99\times10^{-7}\) at N=320. The formal-stencil H maximum falls to \(1.60\times10^{-7}\), and the fixed annulus \(2\le r<5\) H maximum reaches \(7.95\times10^{-9}\).

However, the first-cell production \(\dot A_a\) remains approximately \(4.0\)–\(4.14\times10^{-5}\) from N=80 through N=320. Since \(r_0=\Delta r/2\), \(\dot A_a/r_0^2\) grows with refinement instead of approaching a bounded centre limit. Therefore the better initial H residual **does not show that the regular spherical trace-free series is preserved**.

The matched one-step test supports that distinction:
- All three arms (pinned baseline, origin-half-cell-only, and origin-plus-uniform-fourth-order intervals) take one ordinary production step at N=40/60/80/120/160.
- The uniform-fourth-order arm retains substantially smaller H errors after that step, with post-step full-domain H maxima falling from \(3.85\times10^{-4}\) at N=40 to \(5.93\times10^{-6}\) at N=160; post-step formal-stencil H maxima fall to \(1.43\times10^{-7}\) at N=160.
- Momentum and connection residuals are nonzero after the normal step, and the centre \(\dot A_a/r^2\) scaling criterion still fails.
- At N=240 and N=320, **all three arms** fail at the first radiation cell during the RHS/update with an inadmissible state (\(|S|/E>1\)). The workflow consequently ended with case failures, while preserving its artifact. This is a high-resolution matter-admissibility blocker, not evidence that one quadrature arm caused the failure.

The independent regular-centre Ricci/\(A_a\) identity test converges at approximately fourth order on smooth manufactured profiles. That negative control argues against a generic sign/parity defect in the vendor operator; it does not validate the actual initialized centre profile. A discrete-consistent initial metric also fails to make the momentum-constraint tangent vanish. No single production equation has been proved to be the unique cause.

**Decision:** keep all quadrature candidates diagnostic-only. Do not merge a production repair, change physics/gauge/source coefficients, impose a centre projection, or start t=160.

## 2. Exact artifacts and reproducibility ledger

All reported numerical audits use the frozen source and pinned vendor commits stated above unless explicitly identified as a separate invariant-diagnostic branch.

| Audit | Run | Artifact ID | SHA-256 of Actions artifact ZIP |
|---|---|---:|---|
| Baseline production campaign through t=12, N=40/60/80 | [38208253687](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38208253687) | as archived in existing handoff | \`dbbaa3be94242cdeecd437d0a41c00d4cf3f2d625e1f3aaaad9e1b45833f2047\` |
| Discrete-consistent initial-metric A/B | [38284401548](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38284401548) | linked in run | The full evidence is retained in the paired checkpoint artifact |
| Discrete-initial ordinary one-step and saved time series | [38285057887](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38285057887), [38287772729](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38287772729) | linked in run | See the full A/B handoff on the active diagnostic branch |
| Discrete-consistent initial signed K/CMC-target budget | [38288404354](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38288404354) | in run artifacts | Signed term budget closes arithmetically; this is not a propagation proof |
| Manufactured momentum-constraint/operator identity | [38288419513](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38288419513) | \`11763654184\` | \`83c738421ca014b781d5a07d1218d3de116ddc43bea5bc4d762233bc1dfb2d77\` |
| Discrete-initial centre-series RHS rate audit | [38289053751](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38289053751) | \`11763734669\` | \`c29969262e7bd79a8dde8ace357945685a5b5a418db3344f5aa08fda94b52021\` |
| Production-initial trace-free Ricci cross-check | [38289275758](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38289275758) | \`11764109840\` | \`d0525ab5494c7230c428bcaba1503aab3f12e7720316aefd4c50bbb604b82718\` |
| Frozen \(\Lambda\)/Ricci sensitivity | [38289657512](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38289657512) | \`11764557943\` | \`5cb556cc9a9311fb6ba1a86b1daf398e4c01c9ffeb7f458960ff45e767bb47a2\` |
| Independent manufactured regular-centre Ricci/\(A_a\) RHS identity | [38290099699](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38290099699) | \`11765776185\` | \`345f91141bb19920051935b9c5c516bf1ddc798046012a36583426d9e5e34e97\` |
| Origin half-cell quadrature A/B, frozen initial data | [38290342829](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38290342829) | \`11766265011\` | \`f0b217ee820a1d414740099822a912151aced36faa67c1876c09083d962f3a46\` |
| Discrete-initial multi-constraint RHS tangent | [38290633971](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38290633971) | \`11766405316\` | \`2011162ba822bea1b52d2b8c2099284f97f7a0ca5e76e2d2bdc1788c324a7186\` |
| Origin-only matched ordinary one-step propagation | [38290654079](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38290654079) | \`11766655192\` | \`a511f67db4d9a1a12ce27652c82a200dbb4b876266fa3045bf8f995ddc087987\` |
| Origin quadrature refinement, N=40–320, no evolution | [38291486881](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38291486881) | \`11767325140\` | \`42d898b6f195bdd823fba35a240cdb58bb557146c828c87f8b3a8b2a58636638\` |
| Uniform fourth-order radial quadrature, initial profiles only | [38291938761](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38291938761) | \`11767500772\` | \`937a139ffd360ccae753f06572da34a400c26183dbdd0f308f1f8d67ef2bc918\` |
| Matched one-step comparison of all three quadrature arms | [38292762118](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38292762118) | \`11766972176\` | \`0faff832630fe189143e2aa178013967e79d30de5c0b1d0276447691fee7993e\` |

The uniform-quadrature one-step action result is **failure by design of the admission report** because the N=240/320 cases produce an inadmissible radiation state. Treat the artifact as valuable partial diagnostic data, not as a passing one-step gate.

## 3. Initial-data quadrature: exact numerical findings

### 3.1 What was changed in the diagnostic arm

No continuum equation or matter source changes. Only the radial integral used when reconstructing the initial metric differs.

The origin half-cell candidate integrates the first cell with a polynomial in \(r^2\) exact for even source polynomials through \(r^4\), using the weights already recorded in the origin-quadrature handoff. The uniform arm then replaces each centre-to-centre trapezoid with a four-point cubic Lagrange integral:
- Interior interval weights at nodes \(i-2,\ldots,i+1\): \((-1,13,13,-1)/24\), multiplied by \(\Delta r\).
- First interval: even-reflection ghost convention, with weights \((12,13,-1)/24\) on the first three stored samples.
- Final interval: one-sided weights \((1,-5,19,9)/24\).
- The existing six fixed-point source iterations, matter initialization, and later CMC lapse solve are retained.

This is an analytic interpolation/integration rule, not a fitted coefficient. The variants use the exact pinned production code for the subsequent RHS/step.

### 3.2 First-cell H reduction

Signed H is defined with the same production diagnostic at every resolution and arm.

| N | Pinned baseline \(|H_0|\) | Origin half-cell only \(|H_0|\) | Uniform fourth-order intervals \(|H_0|\) |
|---:|---:|---:|---:|
| 40 | \(3.7692\times10^{-3}\) | \(3.9264\times10^{-4}\) | \(3.8648\times10^{-4}\) |
| 60 | \(3.9316\times10^{-3}\) | \(1.3137\times10^{-4}\) | \(1.0031\times10^{-4}\) |
| 80 | \(3.9886\times10^{-3}\) | \(8.1722\times10^{-5}\) | \(3.6766\times10^{-5}\) |
| 120 | \(4.0281\times10^{-3}\) | \(6.5995\times10^{-5}\) | \(1.0569\times10^{-5}\) |
| 160 | \(4.0416\times10^{-3}\) | \(6.5061\times10^{-5}\) | \(6.8805\times10^{-6}\) |
| 240 | \(4.0510\times10^{-3}\) | \(6.5985\times10^{-5}\) | \(1.0699\times10^{-6}\) |
| 320 | \(4.0543\times10^{-3}\) | \(6.6634\times10^{-5}\) | \(4.9875\times10^{-7}\) |

The formal-stencil maximum over \(i\ge2\) with uniform quadrature is:
- N=40: \(1.495\times10^{-5}\)
- N=60: \(6.041\times10^{-6}\)
- N=80: \(4.273\times10^{-6}\)
- N=120: \(1.309\times10^{-6}\)
- N=160: \(3.805\times10^{-7}\)
- N=240: \(2.208\times10^{-7}\)
- N=320: \(1.596\times10^{-7}\)

At N=320, the fixed-annulus \(2\le r<5\) H maximum is \(7.95\times10^{-9}\) on the uniform arm, compared with \(1.23\times10^{-6}\) for baseline. The origin-only arm by itself plateaus in the first-cell H residual near \(6.5\)–\(6.7\times10^{-5}\), which is why the additional interval integration was tested.

**Interpretation:** the initial radial reconstruction is a real contributor to the initial H residual. But H improvement is only one gate; it does not establish the regularity or propagation identities of the evolved RHS.

### 3.3 The centre-series \(A_a\) RHS does not pass

For a regular spherical trace-free tensor, \(A_a=O(r^2)\), so the RHS should preserve that class: \(\dot A_a=O(r^2)\). In the uniform arm:

| N | First-cell r | First-cell \(\dot A_a\) | First-cell \(\dot A_a/r^2\) |
|---:|---:|---:|---:|
| 40 | 1.0000 | \(+1.4570\times10^{-4}\) | \(+1.4570\times10^{-4}\) |
| 60 | 0.6667 | \(+5.1839\times10^{-5}\) | \(+1.1664\times10^{-4}\) |
| 80 | 0.5000 | \(+3.8368\times10^{-5}\) | \(+1.5347\times10^{-4}\) |
| 120 | 0.3333 | \(+3.7292\times10^{-5}\) | \(+3.3563\times10^{-4}\) |
| 160 | 0.2500 | \(+3.8843\times10^{-5}\) | \(+6.2149\times10^{-4}\) |
| 240 | 0.1667 | \(+4.0613\times10^{-5}\) | \(+1.4621\times10^{-3}\) |
| 320 | 0.1250 | \(+4.1368\times10^{-5}\) | \(+2.6475\times10^{-3}\) |

The raw rate tends toward a nonzero value around \(4.1\times10^{-5}\), rather than scaling as \(r^2\). The quotient therefore rises strongly with refinement. Origin-only quadrature is no better on this criterion: its raw first-cell RHS approaches about \(-9.5\times10^{-5}\), while the quotient magnitude rises. This establishes that the quadrature repair is incomplete; it does not establish which curvature/stage expression must change.

## 4. Matched one-step propagation: latest outcome

### 4.1 Runs and coverage

The matched test used \(dt=0.03\Delta r\), identical to the prior bounded campaign's CFL convention, with one ordinary production step per arm/resolution. No field was projected, no damping/clamp was added, and no production RHS was changed.

At N=40/60/80/120/160, every arm took a normal step. The uniform arm retained the H improvement:
 
| N | Uniform-arm first-cell \(|H|\) before step | after step | Uniform-arm formal-stencil max \(|H|\) after step |
|---:|---:|---:|---:|
| 40 | \(3.865\times10^{-4}\) | \(3.847\times10^{-4}\) | \(1.542\times10^{-5}\) |
| 60 | \(1.003\times10^{-4}\) | \(1.009\times10^{-4}\) | \(5.912\times10^{-6}\) |
| 80 | \(3.677\times10^{-5}\) | \(3.776\times10^{-5}\) | \(3.754\times10^{-6}\) |
| 120 | \(1.057\times10^{-5}\) | \(1.019\times10^{-5}\) | \(7.868\times10^{-7}\) |
| 160 | \(6.881\times10^{-6}\) | \(5.926\times10^{-6}\) | \(1.426\times10^{-7}\) |

The full-domain post-step \(|M|\) / \(|C_\Lambda|\) maxima for the uniform arm remain finite but nonzero:
 
| N | post-step max \(|M|\) | post-step max \(|C_\Lambda|\) |
|---:|---:|---:|
| 40 | \(1.530\times10^{-5}\) | \(9.321\times10^{-7}\) |
| 60 | \(5.850\times10^{-6}\) | \(2.390\times10^{-7}\) |
| 80 | \(6.123\times10^{-6}\) | \(1.875\times10^{-7}\) |
| 120 | \(6.185\times10^{-6}\) | \(1.262\times10^{-7}\) |
| 160 | \(6.166\times10^{-6}\) | \(9.437\times10^{-8}\) |

At N=40–160, the uniform arm generally improves H and often M/C relative to baseline, but full-domain norms include centre cells and cannot replace the required regularity and matched fixed-radius tests. The centre-series \(\dot A_a/r^2\) test still fails.

### 4.2 N=240 and N=320 fail all initializer arms on radiation admissibility

For N=240 and N=320, the pinned baseline, origin-only arm, and uniform arm all fail during the unchanged radiation RHS/one-step path at cell \(i=0\). The error says the radiation RHS has produced an inadmissible state with \(|S|/E>1\):
- N=240: \(E\approx8.81887\times10^{-9}\), \(|S|\approx8.9976\times10^{-9}\), ratio about \(1.0203\).
- N=320: \(E\approx8.81887\times10^{-9}\), \(|S|\approx1.1996\times10^{-8}\), ratio about \(1.3603\).

This outcome is common to all three initialization variants and should be tracked as a separate first-cell radiation admissibility failure. It is not cured by these metric quadrature alternatives. Since those cases fail, there is no admitted one-step result for N=240/320; do not silently include them as successful evolution.

## 5. Constraint RHS tangent: a constraint-consistent initial slice still regenerates momentum residual

A separate tangent audit evaluates the full production semi-discrete RHS on six frozen t=0 checkpoints (baseline/candidate × N=60/80/120). It perturbs state sectors by \(\pm\epsilon\) with \(\epsilon=2\times10^{-5}\) and \(10^{-5}\), re-solves the CMC lapse, finite-differences signed H/M/\(C_\Lambda\) terms, and compares one normal production step from a separate copy. Source: frozen commit \`6d8b3f2790644591326df57c55908863f994cb0c\`; vendor \`d6052d605673ce9d82cdc99b0df79855fa2ea215\`.

The discrete-consistent candidate has initial \(\max|H|\) near roundoff (N=60/80/120: \(7.65\times10^{-16}\), \(1.29\times10^{-15}\), \(3.16\times10^{-15}\)); M and \(C_\Lambda\) are also at roundoff. Yet at the first admitted cell:
 
| N | first admitted r | candidate tangent \(\dot M\) | one-step \(\Delta M/\Delta t\) |
|---:|---:|---:|---:|
| 60 | 3.333 | \(-8.7553\times10^{-5}\) | \(-8.6967\times10^{-5}\) |
| 80 | 2.500 | \(-1.2354\times10^{-4}\) | \(-1.2295\times10^{-4}\) |
| 120 | 1.667 | \(-1.9063\times10^{-4}\) | \(-1.8993\times10^{-4}\) |

The tangent/step agreement is very close at those witnesses. It places the leading issue in the coupled semi-discrete RHS rather than only in a reporting finite difference or timestep quotient. Since the first-admitted cell moves inward with resolution, the table alone is not fixed-radius convergence.

At the candidate N=120 first admitted cell, the signed budget gives approximately:
- geometric momentum tendency: \(-3.6351\times10^{-4}\);
- S-current term contribution to the residual tendency: \(+1.7290\times10^{-4}\);
- net \(\dot M\): \(-1.9063\times10^{-4}\).

The geometry terms \(\mathrm{D}_r\dot A_a\), \(-\frac23\mathrm{D}_r\dot K\), radial-b, and \(6A_a\chi_r\) reconstruct the geometric budget at roundoff; source-sector assembly also closes at roundoff. But the geometry/current tendencies do not cancel. The scalar-side current tendency is dominated by local S and changes relatively little between the paired baseline/candidate states; most of the A/B shift lies on the geometry side.

**Limit:** this is a semi-discrete diagnostic, not an independent proof of covariant stress-energy conservation and not proof that any one term is wrong.

## 6. Independent Ricci/Aa controls and what they rule out

### 6.1 Smooth manufactured profiles pass the pinned operator convergence control

Run [38290099699](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38290099699) uses smooth analytic even metric/lapse profiles with \(ab^2=1\), analytic regular connection data, and an independently computed curvature route from proper radial distance and areal radius. It takes no evolution step. The flat control is roundoff-sized; on the fixed physical annulus, the descriptive convergence orders are about:
- trace-free Ricci: 4.017 (analytic Lambda) and 4.031 (discrete-constraint Lambda);
- trace-free lapse Hessian: 4.006;
- complete \(A_a\) RHS: 4.017 and 4.030.

This is a negative control against a generic vendor Ricci/parity defect on smooth regular data. It does not validate the production-initialized centre profile.

### 6.2 Production-initialized profile is the outlier

On the exact pinned production initial profile, the first-cell raw \(\dot A_a\) is about \(+1.09\times10^{-3}\) in the baseline over N=40–320, yielding \(\dot A_a/r^2\) that rises from \(1.10\times10^{-3}\) at N=40 to \(6.96\times10^{-2}\) at N=320. The half-cell quadrature changes raw first-cell \(\dot A_a\), but it changes sign and does not restore \(O(r^2)\) scaling. The uniform arm lowers the raw rate to about \(4.1\times10^{-5}\), but that too does not tend to zero at the required rate.

The current evidence therefore separates:
1. a genuine initial metric integral/quadrature contribution to H;
2. a persistent initialized-centre \(A_a\) RHS series inconsistency;
3. later Lambda/connection sensitivity to momentum-constraint growth;
4. an independent high-resolution radiation admissibility failure.

Do not conflate these into a single proven defect or replace the production Ricci operator with the independent route.

## 7. Rejected / demoted leads

- **Generic parity/sign error in the shared derivative:** disfavoured by manufactured smooth-profile operator tests with descriptive near-fourth-order convergence and roundoff signed-budget reconstruction.
- **CMC centre row as primary propagation cause:** demoted. Even-parity centre-row and production-matched fourth-order lapse counterfactuals reduce the algebraic first-cell K-RHS mismatch but increase admitted-inner M after one normal step by about \(3\times10^{-6}\); neither justifies a CMC patch.
- **Lambda-only projection/reconstruction:** not justified. At t=0 the connection constraint is already roundoff while the production \(A_a\) RHS mismatch exists; at t=2.16 Lambda sensitivity is resolution/sign dependent.
- **Scalar lapse-gradient inverse metric candidate:** the declared \(\gamma_{rr}=a/X^2\) implies \(\gamma^{rr}=X^2/a\), while archived scalar RHS has \(X^{-2}\) for the lapse-gradient contraction. This remains a legitimate covariance audit item, but the matched candidate's trajectory changes are tiny and constraint effects mixed; it does not explain the main defect by itself. Keep separate.
- **Origin half-cell alone as complete fix:** rejected. It reduces first-cell H but redistributes error to adjacent/first-admitted cells and does not fix the \(A_a\) RHS series.
- **Discrete-consistent H initializer as complete fix:** rejected. H=0 is imposed by the Newton solve at t=0; the unchanged RHS immediately develops a nonzero M tangent.
- **Areal-radius parity correction:** a separate invariant-diagnostic correction confined to cells \(i=0,1\); it fixes a false central trapping indicator and does not repair H/M/\(C_\Lambda\) propagation.
- **Failed timeline report join:** run 38285163101 failed with a \`KeyError: (40, 0.72, 'baseline')\` in report-join logic. It is a reporting/coverage bug, not a production physics result.

## 8. Current status and next bounded gate

### What has been established
1. The N=40/60/80 campaign reached the t=12 horizon with finite state/diagnostics. This remains execution completion only; H/M/\(C_\Lambda\), Misner-Sharp balance and resolution convergence are not admitted.
2. Origin quadrature significantly contributes to the initial H error.
3. The uniform fourth-order interval integral yields the best tested initial H profile, and the H improvement survives one normal step through N=160.
4. A constraint-consistent initial H/M/\(C_\Lambda\) slice still develops a nonzero momentum RHS tendency under the unchanged production vector field.
5. Smooth manufactured data do not show a generic vendor Ricci/trace-free operator defect.
6. The actual initialized profile does not demonstrate the required first-cell \(O(r^2)\) \(\dot A_a\) series.
7. The radiation primitive/admissibility gate fails in the first cell for one-step N=240/320 cases, independent of quadrature arm.

### Next steps — read-only, no long trajectory
1. Decompose the actual production \(A_a\) RHS at the **baseline and uniform-initialized profiles**, both t=0 and already saved t=2.16: trace-free Ricci, trace-free lapse Hessian, anisotropic stress, nonlinear/extrinsic term, and numerical dissipation if active. Preserve signs and report \(r^2\dot A_a\) and \(\dot A_a/r^2\) at cells i=0,1,2 plus fixed r=2,3,4,5.
2. In the same frozen states, split \(\dot K\) into lapse Laplacian and nonlinear/extrinsic-curvature plus matter terms, and check the series for \(\mathrm{D}_r\dot K/r\). Verify that both the analytic reference and discrete expression use the same metric/conformal conventions.
3. Independently audit first-cell radiation source/flux tendencies with the same source state at N=160/240/320; report \(E-|S|\) before/after each RHS stage and term closure. Do not loosen the radiation cone/admissibility criterion or clip/project the state.
4. Keep the uniform quadrature only as an experiment. If a regular-centre identity is shown to fail in a specific term, run one bounded one-step counterfactual changing only that numerical layer; compare H, M, \(C_\Lambda\), determinant, fixed radii, matter admissibility and Misner-Sharp balance.
5. Only after those gates pass should an evolution campaign be resumed. No t=160 campaign yet.

## 9. Invariants and project guardrails

- Production physics, EOS, scalar potentials, DM coupling, matter definitions, gauge, Lambda equation and boundary conditions remain unchanged by these diagnostic arms.
- No hand-tuned centre correction, fitted coefficient, arbitrary damping, floor/clamp, centre projection, made-up source, imposed bounce/turnaround, reset, D-to-matter law or artificial stop.
- Do not claim covariant source conservation from the source budget's arithmetic closure alone.
- Do not fit formal convergence order to sign-changing centre-cell values or compare “first admitted cell” across N as if it were a fixed physical radius.
- Distinguish workflow success, numerical stage completion, scientific admission and physical interpretation.
- The authoritative REV16/GEAR archive and the full bidirectional-machine objective are preserved. The current “time-tube” picture remains an interpretive guide, not an extra equation to impose.


---

## Latest result — signed production Aa/K RHS component split (2026-10-12)

- Run: [38294736014](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38294736014)
- Artifact: \`production-Aa-K-centre-source-tangent\`, ID \`11768515284\`
- Artifact ZIP SHA-256: \`2128028ad4fcd80c4d603107d481a78a1add233d18b35eaa3d1d19a43cae2a91\`
- Detailed report: [PRODUCTION_AA_K_RHS_CENTRE_SOURCE_SPLIT_2026-10-12.md](https://github.com/miahpenn/0-Cosmos-engine/blob/diag/discrete-consistent-initial-production-test-2026-10-12/docs/PRODUCTION_AA_K_RHS_CENTRE_SOURCE_SPLIT_2026-10-12.md)
- Frozen production source and vendor remain \`6d8b3f2790644591326df57c55908863f994cb0c\` / \`d6052d605673ce9d82cdc99b0df79855fa2ea215\`.

### What the signed terms show

The audit covers 17 frozen states: baseline and discrete-consistent-initial-metric states at t=0 and saved t=2.16 for N=60/80/120, plus the uniform-fourth-order quadrature t=0 profiles for N=40/60/80/120/160. No state is evolved.

1. In the pinned baseline initial profile at N=60/80/120, first-cell \(\dot A_a\approx1.09\times10^{-3}\) is dominated by the trace-free Ricci contribution \(\alpha(R^r{}_r-R/3)\approx1.05\times10^{-3}\); the trace-free lapse Hessian is only \(3.6\)–\(4.1\times10^{-5}\).
2. The discrete-consistent metric and uniform-quadrature initial profiles largely remove the first-cell Ricci contribution at these resolutions. Yet the uniform arm's trace-free lapse-Hessian contribution remains \(3.87\)–\(4.15\times10^{-5}\) from N=80–160, and its total \(\dot A_a\) remains \(3.84\)–\(3.88\times10^{-5}\). The quotient \(\dot A_a/r_0^2\) grows from \(1.53\times10^{-4}\) at N=80 to \(6.21\times10^{-4}\) at N=160. The centre-series gate therefore remains open.
3. On the uniform arm, \(\dot K(r_0)\) is about \(-1.38\) to \(-1.52\times10^{-3}\), while the next cell is about \(+1.0\)–\(1.45\times10^{-4}\). The lapse-Laplacian term changes from negative at i=0 to \(+0.96\)–\(1.29\times10^{-3}\) at i=1, while the matter-source term changes much less. \(D_r\dot K/r\) at i=0 rises from \(5.11\times10^{-4}\) (N=40) to \(8.68\times10^{-3}\) (N=160).
4. Signed RHS reconstruction closes to roundoff: maximum base Aa/K reconstruction error \(2.17\times10^{-19}\), momentum geometry closure 0, per-sector momentum-current closure \(4.07\times10^{-20}\). The two epsilon tangents agree closely (largest full Aa tangent difference \(2.84\times10^{-10}\); full K tangent difference \(8.54\times10^{-10}\)). This validates accounting of the implemented RHS, not continuum conservation.

### Current interpretation and next gate

The old baseline's dominant initial central term is Ricci; high-order metric quadrature largely suppresses that. The persistent residual is now most directly visible in the trace-free lapse-Hessian term. The \(K\)-RHS central-to-neighbour jump also localizes attention on the CMC lapse solution / central-row compatibility. That is a **lead, not yet a proved faulty equation**: manufactured smooth-profile tests pass the pinned trace-free lapse-Hessian and Ricci operators at approximately fourth order on a fixed annulus, and prior CMC-row counterfactuals improved a local K witness but worsened an admitted-inner momentum witness.

Next, do a read-only audit of the exact pinned \`_solve_lapse\` centre row: record \(\alpha_i\), even-series defect, row equation residual, and the resulting \(-\Delta\alpha\) for i=0,1,2 and fixed r=2,3,4,5 across N=40–160; compare to independent smooth even lapse profiles. In the same frozen states, reconstruct regular-series identities for trace-free Hessian, trace-free Ricci, \(\dot A_a\) and \(\dot K\). Do not modify the row or project the centre until a discrete identity is shown to fail with term-level closure.

The N=240/320 first-cell radiation inadmissibility discovered by the matched quadrature one-step run remains a separate open blocker. No t=160 campaign; no production physics/source/gauge/boundary changes; \`main\` remains untouched.


---

## Latest CMC closure and centre-row checks — 2026-10-12

### Target-closure rank probe

- Run: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38294955148
- Artifact ZIP SHA-256: 64cc1666a2010d73b462cac48c78706ac48be1a1d02fea7de343bd09172c4e24
- Report: https://github.com/miahpenn/0-Cosmos-engine/blob/diag/independent-metric-ricci-audit/docs/CMC_TARGET_CLOSURE_RANK_AUDIT_2026-10-12.md

Nine fixed Kdot targets over 6.12947e-5 to 6.30289e-5 yield a target-output slope of 0.999999999757 and closure residual slope -2.42686e-10. The H residual remains around -8.58e-14; interior PDE rows close to 3.71e-15 max. The outer normalization cell is not a PDE row but its K RHS enters the target average; its residual is around -4.58e-12 and it carries 1.8629% of the R80 projection weight. H equals the weighted mean pointwise residual to about 1.14e-20. This indicates near-redundancy of the current target closure on the frozen slice, not proof of a continuum family or authority to change the target definition.

### Read-only CMC lapse centre-row/operator audit

- Run: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38296556473
- Artifact ID: 11768278320
- Artifact ZIP SHA-256: e02c2c0f384d20692fc02bb3ee771409d0743a2479dc39276ee3a14b1156101f
- Report: https://github.com/miahpenn/0-Cosmos-engine/blob/diag/independent-metric-ricci-audit/docs/CMC_LAPSE_CENTRE_ROW_AUDIT_2026-10-12.md

On frozen archived t=0 R80/R160 states, both lapse profiles remain positive, alpha fits an even polynomial in r over the first eight cells with max residual 4.49e-8, and the production K-RHS / CMC elliptic rows close to roundoff (max 6.04e-15 / 7.32e-15). Most importantly, the production lapse-geometric K term and the independently reconstructed CMC lapse term agree to around 2e-15 at the first cells and fixed radii r≈2,3,4,5. This weakens the hypothesis that the archived-state CMC centre row is itself faulty. It does not test the uniform-quadrature initializer arm and is not a refinement-order proof.

### Combined interpretation / next gate

The centre-series issue is still open, but the present evidence does not support patching the CMC centre row. Keep the target-closure rank issue separate: compare its discrete derivation with the authoritative CMC condition, including the outer normalization cell. Next run the same CMC-versus-production lapse-operator comparison on the exact uniform-fourth-order quadrature initial profiles N=40/60/80/120/160, then compare trace-free Hessian, Ricci and K source terms on identical frozen states. If those operators agree there too, shift focus to initializer geometry/matter construction and the regular-centre identities for Aa/K/X/a/b. Radiation realizability failures at N=240/320 remain a separate admission blocker. No physics/source/gauge changes, centre projection, damping, floors/clamps, or long campaign until the specific defect and all admission gates are demonstrated.

All work above is diagnostic-only. Production main is unchanged.


## Main-branch reconciliation (2026-10-12)

Current main HEAD is 87e71c3793d8afd5519ca1099dc69ebff707546b. Compared with the recorded production baseline 809228c380d80af83858412557c78b6329ffa5c8, the only transient difference was the diagnostic workflow .github/workflows/cmc-target-closure-rank-probe.yml, which was added in commit 6d3c508 and removed in the current revert commit. The current main tree has no production source, physics, or kernel changes relative to that recorded baseline. Diagnostic scripts/reports remain on diagnostic branches.
