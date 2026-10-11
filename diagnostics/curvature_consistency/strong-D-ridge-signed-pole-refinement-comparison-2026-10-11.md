# Signed Hamiltonian-residual pole comparison: N=160 vs N=320

Date: 2026-10-11. Diagnostic-only post-processing; no new trajectory and no production code changes.

## Reproducible source runs

- N=160, Rmax=80: Actions run 38150925743, artifact 11706797173, SHA-256 8b118d1f81e01f2a5d291d160b8bac42abc6a11051d8453efd9c19d0826d6d8b.
- N=320, Rmax=80: Actions run 38165315602, artifact 11712817995, ZIP SHA-256 d266c9aedf259222e720cd52067b97ab5cc44f7556f98177db5e00e656ee449f.
- Both are production trajectories with D=1e-4, scalar amplitude 0.01, width 7, radiation ON, CFL=0.0075, t=0..12. Spatial spacing halves from dr=0.5 to 0.25 while the fixed CFL also halves dt. The N=320 workflow verified Hamiltonian term reconstruction closure at every saved slice. Two resolutions are not a proof of continuum convergence.

## Early signed-pole pattern is resolution-sensitive

Radial sign-transition counts for H_vendor, measured over 1.25 <= r <= 79.75, differ materially in the early low-amplitude interval:

| Requested t | Actual t, N=160 | Actual t, N=320 | Sign transitions N=160 / N=320 |
|---:|---:|---:|---:|
| 5.50 | 5.50125 | 5.49938 | 2 / 4 |
| 5.75 | 5.74875 | 5.75063 | 2 / 0 |
| 6.00 | 6.00000 | 6.00000 | 2 / 0 |
| 6.25 | 6.25125 | 6.24938 | 0 / 0 |

The local positive maxima also differ:

| Requested t | Positive local maxima at N=160 (r, H) | Positive local maxima at N=320 (r, H) |
|---:|---|---|
| 5.50 | (9.75, +2.50e-7) | (4.875, +1.40e-8), (8.625, +1.60e-7) |
| 5.75 | (9.25, +3.21e-7) | (5.125, +4.68e-7), (7.125, +3.63e-7) |
| 6.00 | (5.75, +1.38e-7), (8.25, +5.60e-7) | (5.375, +1.30e-6) |
| 6.25 | (5.75, +1.76e-6) | (5.375, +2.84e-6) |

At N=160, two distinct positive maxima coexist near r=5.75 and r=8.25 at t=6.00; the r=8.25 maximum is gone by t=6.25 while the r=5.75 maximum wins. At N=320 this exact pairing is not reproduced: the maxima and sign transitions change earlier, and by t=6.00 only one positive local maximum remains in 4 <= r <= 20. Thus the existence of a local sign change and pole/peak switching is confirmed in the N=160 slice data, but the precise early sign topology and pole identities are not spatial-resolution-stable. Do not interpret the N=160 zero-time estimate as physical.

## Late outer residual ridge is much more stable

| Requested t | Ridge peak N=160 (r, H) | Ridge peak N=320 (r, H) | Amplitude change N=320 vs N=160 |
|---:|---|---|---:|
| 7.50 | (5.75, +3.34049e-5) | (5.875, +3.42468e-5) | +2.520% |
| 8.00 | (6.25, +5.87044e-5) | (6.125, +5.92557e-5) | +0.939% |
| 8.50 | (6.75, +8.82380e-5) | (6.625, +8.83675e-5) | +0.147% |
| 8.75 | (6.75, +1.03498e-4) | (6.875, +1.03819e-4) | +0.310% |
| 9.00 | (7.25, +1.19173e-4) | (7.125, +1.19350e-4) | +0.148% |
| 10.00 | (8.25, +1.80231e-4) | (8.125, +1.80434e-4) | +0.113% |
| 12.00 | (9.75, +2.82918e-4) | (9.875, +2.82763e-4) | -0.055% |

The late outer ridge closely agrees across the two resolutions in amplitude and radius. Its measured movement is in the Hamiltonian constraint residual; that robustness does not turn it into a material shell or prove a physical propagation law.

## Verdict

**Early sign/pole exchange:** observed in the N=160 profile, but not sufficiently stable under N=320 refinement for physical interpretation. Treat it as resolution-sensitive until further signed-profile analysis establishes otherwise.

**Late outward residual ridge:** survives this fixed-domain refinement with nearly identical late amplitudes and neighboring peak radii. It remains a constraint-residual feature, not evidence by itself of D-to-matter conversion, a physical front, bounce, turnaround, or cycle.

The refinement reduces no uncertainty about the missing GEAR couplings: D-to-matter and local-to-cosmic transfer remain open gates, not assumed terms.
