# Matched-radius evolved-state comparison: full continuous-source profile — 2026-10-25

**Status:** post-processing of existing accepted-state artifacts; no evolution rerun and no production changes. All arms reached the same coordinate time \(t=0.60\). The full continuous-source arm is the primary research comparator; the partial-joint arm is retained only as a control.

## Provenance

- Partial-joint/full-evolution artifacts: [run 40153734946](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40153734946), artifact IDs N60=12584100185, N120=12583184844, N240=12583826100.
- Full continuous-source artifacts: [run 40156226086](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40156226086), artifact IDs N60=12584064724, N120=12584164579, N240=12584980156.
- Source pin: \`4aec88d0e3ff028e0a65a8ac36b30948625df792\`; vendor pin: \`d6052d605673ce9d82cdc99b0df79855fa2ea215\`.
- The underlying workflow reports verify protected engine/vendor/dependency files were unchanged and \`scientific_admission=false\`.

## Method

The archived \`profiles.npz\` arrays were read from the six original workflow artifacts. Profiles were linearly interpolated from each resolution's cell-centre radii to a common 1,000-point radius set spanning \(r=0.5\) to \(10\), then max-absolute, RMS, and cross-resolution differences were computed. Fixed-radius samples below are also interpolated, never compared by cell index. This is an offline comparison of the already accepted \(t=0.60\) states.

## Main result: full-profile Hamiltonian residual converges strongly at matched physical radii

For the full continuous-source arm:

| N | max \(|H|\), \(0.5\le r\le10\) | RMS \(H\), \(0.5\le r\le10\) | RMS \(H\), \(2\le r\le10\) |
|---:|---:|---:|---:|
| 60 | 2.1104e-5 | 7.4423e-6 | 4.6921e-6 |
| 120 | 1.6685e-6 | 5.1426e-7 | 3.3699e-7 |
| 240 | 6.4989e-8 | 2.7801e-8 | 2.1795e-8 |

The max-absolute Hamiltonian residual falls by factors of about 12.6 and 25.7 on successive refinements. At fixed physical radius:

| Radius | H (N=60) | H (N=120) | H (N=240) |
|---:|---:|---:|---:|
| 1 | -1.7395e-5 | -1.1115e-6 | -5.5241e-8 |
| 2 | -6.2698e-6 | -4.0290e-7 | -2.1123e-8 |
| 3 | +2.6738e-6 | +2.5025e-7 | +1.7838e-8 |
| 5 | +7.9285e-6 | +5.6591e-7 | +3.5245e-8 |
| 8 | -1.9497e-6 | -1.7406e-7 | -1.2362e-8 |

The full-profile momentum residual \(M\) and Lambda-identity residual \(C_\Lambda\) also improve across these matched-radius profiles:

| Witness | N=60 max abs | N=120 max abs | N=240 max abs |
|---|---:|---:|---:|
| \(M\) | 6.3169e-6 | 4.3960e-7 | 2.9939e-8 |
| \(C_\Lambda\) | 3.9073e-6 | 2.8226e-7 | 2.1467e-8 |

These maxima use the same interpolated \(r=0.5\ldots10\) window.

## Control comparison: partial-joint arm does not show the same convergence

For the partial-joint arm, the max-absolute \(H\) over the same common-radius window is \(2.8871\times10^{-4}\), \(3.8333\times10^{-4}\), and \(3.8767\times10^{-4}\) for N=60/120/240. Its RMS over \(r=2\ldots10\) does decrease, but its near-origin max remains non-monotone and roughly two orders of magnitude larger than the full-profile N=240 maximum. This is consistent with keeping it as a control rather than treating the partial correction as a candidate remedy.


## Central regularity witness: full-profile arm converges, partial control does not

The normalized centre witness is \(A_a(r_0)/r_0^2\). Comparing it with the immediately adjacent cell's \(A_a(r_1)/r_1^2\) is a useful local regularity consistency check because it does not compare unlike cell indices across resolutions.

| N | Full-profile centre \(A_a/r_0^2\) | Full-profile adjacent \(A_a/r_1^2\) | Absolute gap |
|---:|---:|---:|---:|
| 60 | -8.21965e-6 | -8.47322e-6 | 2.53564e-7 |
| 120 | -9.88439e-6 | -9.81941e-6 | 6.49810e-8 |
| 240 | -1.00810e-5 | -1.00752e-5 | 5.80029e-9 |

The centre coefficient approaches approximately \(-1.01\times10^{-5}\), while the centre-to-adjacent gap shrinks from 2.54e-7 to 5.80e-9. The difference between successive centre values falls by a factor 8.47 and then the observed order is about 3.08 under each doubling of N. This is evidence of convergence to a finite regular \(A_a/r^2\) coefficient, not evidence that the coefficient itself must vanish.

For the partial-joint control, the accepted centre values are -7.23741e-5 (N=60), +1.62028e-4 (N=120), and -5.48849e-4 (N=240); the sign and magnitude do not converge. The full-profile arm therefore improves both the matched-radius Hamiltonian residual and the local normalized centre regularity witness in the same accepted t=0.60 states.

## Interpretation and limits

This is strong evidence that the full continuous-source metric profile produces a much more resolution-consistent accepted-state geometry in the near-origin-to-inner-profile window than the partial first-half-cell correction. The result now holds in both the stage-level \(C_\Lambda\)-rate probe and the evolved \(t=0.60\) Hamiltonian/momentum/\(C_\Lambda\) profiles at matched physical radii.

It does **not** close the complete matter-coupled momentum-constraint propagation identity or total stress-energy conservation. It also does not independently validate the archive's physical source normalization/sign, determine the authoritative pointwise-versus-integrated scalar spatial operator, or authorize promotion of this diagnostic initializer to production. Keep all source, initializer, vendor and dependency code unchanged; keep \`scientific_admission=false\`.
