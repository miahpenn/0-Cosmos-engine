# Issue #4 copied-initializer Simpson discriminator result

**Result status:** diagnostic-only, all CI gates passed at head
\`43605353cca054e4e69e911f308df44fdfcc398e\`. **Scientific admission remains NOT ADMITTED.**

## Exact-head validation

The specialist run [38853110175](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38853110175) and companion checks all succeeded on this head:

- [Issue #4 specialist](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38853110175)
- [Lightweight validation and one-step smoke](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38853114068)
- [Strict JSON/residual contract](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38853114099)
- [Campaign/radiation-retry integration](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38853114054)

The specialist run passed 818 repository tests and 12 focused mathematical tests. Source/vendor pins, the six-pass production-recurrence reproduction, copied fixed-point convergence, finite profile checks, metric identities, and original-state digests all passed. The retained artifact is [12025674853](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38853110175/artifacts/12025674853), SHA-256 \`e3cc98edc8a583ed95705e541f43886e2f2907b260300dc994cad52c36e5307d\`.

## Actual initialized-source comparison

The diagnostic adds a \`parity_composite_simpson\` arm to the *copied* fixed-point recurrence for the standard, uncharged initializer. The source, scalar state, vendor, engine, and production initializer are source-pinned and never edited; zero production evolution steps are run. The Simpson arm uses the regular even extension and composite Simpson 1/3 plus a final Simpson 3/8 panel, with \(s(0)=1\).

| Resolution | Existing fixed-\(s_2\) polynomial arm: signed centre gap | Copied parity-Simpson arm: signed centre gap | Simpson fixed-annulus maximum | Simpson Hamiltonian-constraint max |
|---:|---:|---:|---:|---:|
| 160 | \(-1.88494\times10^{-4}\) | \(+4.77069\times10^{-6}\) | \(8.84763\times10^{-8}\) | \(2.69618\times10^{-3}\) |
| 240 | \(-1.93165\times10^{-4}\) | \(+2.11393\times10^{-6}\) | \(1.88714\times10^{-8}\) | \(2.70174\times10^{-3}\) |
| 320 | \(-1.94790\times10^{-4}\) | \(+1.18783\times10^{-6}\) | \(6.18941\times10^{-9}\) | \(2.70366\times10^{-3}\) |

The parity-Simpson centre gap decreases with order approximately \(p=2.00\) when resolution is refined (roughly factor four when \(N\) doubles). This is a significant change from the original/fitted arms, whose centre gap tends toward a nonzero value near \(-1.95\times10^{-4}\). It is consistent with PR #284's independent manufactured calibration, where the metric-only route error and inter-route gap converge at about second order at the centre while the BSSN Ricci error converges at about fourth order. The Simpson fixed-annulus gap retains the expected near-fourth-order decrease.

The global Hamiltonian-constraint maximum improves only modestly (about 10–12%) relative to the fixed-\(s_2\) arm and does not vanish. The copied momentum-constraint maximum is about \(10^{-17}\); the connection constraint is zero under the declared representation. Every copied variant remains finite and positive, converges within the declared fixed-point threshold, and leaves the original initialized state digest unchanged. The report retains the full profiles and signed errors.

## Interpretation and next scientific gate

The strongest current explanation is now the point-sampled/cell-centre cumulative source-integration rule: fixing only the first half-cell did not remove the centre plateau, whereas replacing the integration in a copied actual-source recurrence removed that plateau and restored second-order centre-gap convergence. The test is a discriminator, not a proof that this is the only defect and not authorization to apply a production patch. The remaining \(O(h^2)\) centre gap and nonzero Hamiltonian constraint still require source/formula review and convergence assessment.

Do **not** close Issue #4 or claim scientific admission from these curvature comparisons. The independent same-stage total-stress-energy conservation/balance and Lumen's matter-coupled Brown/Bianchi constraint-propagation oracle remain separate mandatory gates. The PR remains draft, with no engine/vendor/dependency edits and \`scientific_admission=false\`.
