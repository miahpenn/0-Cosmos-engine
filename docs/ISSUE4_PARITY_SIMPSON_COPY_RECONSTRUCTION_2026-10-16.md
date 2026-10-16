# Issue #4 copied-initializer Simpson discriminator result

**Current result status:** diagnostic-only, applicable CI gates passed at head
3f584c866e716ec09400eeed6f53d1539273a5de. **Scientific admission remains NOT ADMITTED.**

## Current exact-head validation and constraint interpretation

The latest specialist workflow and the four standard validation workflows all passed on this head:
- [Issue #4 specialist](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38853886332)
- [Lightweight validation + one-step smoke](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38853886259)
- [Strict JSON/residual contract](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38853886194)
- [Campaign/radiation-retry integration](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38853886242)
- [Current-main momentum/Brown diagnostic bundle](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38853886192)
- [Bounded total Valencia stage-budget diagnostic](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38853886362)

The specialist run passed **820 repository tests** and **14 focused mathematical controls**. Source/vendor pins, copied production recurrence reconstruction, fixed-point convergence, finite/positive profiles, metric identities, matter-residual finiteness, and original-state digests all passed. Artifact [12026039922](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38853886332/artifacts/12026039922), ZIP SHA-256 ac27ec017e0e684ba556a17bc53a1ff2ac8df9932b387a666237caff4fd05aab.

Important interpretation correction: the earlier phrase “Hamiltonian-constraint maximum” described the raw vacuum.constraints Hamiltonian metric-only geometric side, not the full matter-subtracted Hamiltonian residual. The updated report now records the geometric side, its canonical matter source, and the subtracted residual separately. With the repository convention, H_total = H_geo − 16πρ_total and M_total = M_geo − 8πj_total. These are same-path initial constraint witnesses, not independent stress-energy conservation or S/D source-normalization proofs.

The separate Brown diagnostic job succeeded as software, but its internal acceptance flag is still acceptance_oracle_satisfied=false; the matter-coupled Brown/Bianchi propagation gate remains open.

## Actual initialized-source comparison

The diagnostic adds a \`parity_composite_simpson\` arm to the *copied* fixed-point recurrence for the standard, uncharged initializer. The source, scalar state, vendor, engine, and production initializer are source-pinned and never edited; zero production evolution steps are run. The Simpson arm uses the regular even extension and composite Simpson 1/3 plus a final Simpson 3/8 panel, with \(s(0)=1\).

| N | Fixed-s2 fit: signed centre Ricci gap | Copied parity-Simpson: signed centre Ricci gap | Simpson fixed-annulus maximum gap | Simpson metric-only H-side maximum | Fixed-s2 matter-subtracted H residual maximum | Simpson matter-subtracted H residual maximum |
|---:|---:|---:|---:|---:|
| 160 | -1.88494e-4 | +4.77069e-6 | 8.84763e-8 | 2.69618e-3 | 3.55752e-4 | 9.67149e-6 |
| 240 | -1.93165e-4 | +2.11393e-6 | 1.88714e-8 | 2.70174e-3 | 3.70172e-4 | 4.30277e-6 |
| 320 | -1.94790e-4 | +1.18783e-6 | 6.18941e-9 | 2.70366e-3 | 3.75276e-4 | 2.41976e-6 |

The parity-Simpson centre gap decreases with order approximately \(p=2.00\) when resolution is refined (roughly factor four when \(N\) doubles). This is a significant change from the original/fitted arms, whose centre gap tends toward a nonzero value near \(-1.95\times10^{-4}\). It is consistent with PR #284's independent manufactured calibration, where the metric-only route error and inter-route gap converge at about second order at the centre while the BSSN Ricci error converges at about fourth order. The Simpson fixed-annulus gap retains the expected near-fourth-order decrease.

The values in the “metric-only H-side” column are not full Hamiltonian residuals. On this latest head, the maximum absolute matter-subtracted H residual in the fixed-s2 fit arm is 3.55752e-4, 3.70172e-4, and 3.75276e-4 at N=160, 240, 320; for parity-Simpson it is 9.67149e-6, 4.30277e-6, and 2.41976e-6. That is a reduction of about 97.3%, 98.8%, and 99.4%. The Simpson first-cell/second-cell signed H residuals are [-3.77706e-6, -9.67149e-6], [-1.62219e-6, -4.30277e-6], and [-9.00992e-7, -2.41976e-6], respectively.

The matter-subtracted momentum residual remains at order 1e-17 for these initialized-state witnesses. The connection and metric-mapping identities pass under the declared representation. These residuals are computed by subtracting the canonical local stress-energy projection from the copied metric-side constraint; this is not independent evidence for stress-energy conservation or the archive S/D source normalization. Every copied variant remains finite and positive, converges within the declared fixed-point threshold, and leaves the original initialized-state digest unchanged. The report retains the signed profiles and raw source operands.

## Interpretation and next scientific gate

The strongest current explanation is now the point-sampled/cell-centre cumulative source-integration rule: fixing only the first half-cell did not remove the centre plateau, whereas replacing the integration in a copied actual-source recurrence removed that plateau and restored second-order centre-gap convergence. The test is a discriminator, not a proof that this is the only defect and not authorization to apply a production patch. The remaining O(h^2) centre Ricci-route gap and nonzero matter-subtracted Hamiltonian residual still require source/formula review and convergence assessment.

Do **not** close Issue #4 or claim scientific admission from these curvature comparisons. The independent same-stage total-stress-energy conservation/balance and Lumen's matter-coupled Brown/Bianchi constraint-propagation oracle remain separate mandatory gates. The PR remains draft, with no engine/vendor/dependency edits and \`scientific_admission=false\`.
