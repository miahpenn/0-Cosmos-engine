# Strong-D domain-size residual audit — full-profile findings

**Date:** 2026-10-11  
**Branch:** `diag/independent-metric-ricci-audit`  
**Workflow:** https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38108719482  
**Run commit:** `9087dbea4c8507894088cc2603c6608b6cfc8623`  
**Artifact:** `strong-D-domain-size-residual-audit`, ID `11692468908`  
**Artifact SHA-256:** `98395c9934e950dda46a70e5a011bcc27778641b8ca02b670880bc367e4204b4`  
**Status:** completed successfully; diagnostic-only.

## Campaign integrity

Four trajectories completed to coordinate time (t=12), each with 3,200 accepted steps and 49 full radial profiles (initial slice plus quarter-time samples). All use (\Delta r=0.5), nominal (\Delta t=0.00375), (D=10^{-4}), scalar amplitude (0.01), width (7), and radiation enabled:

- Production initializer: ((N,r_{max})=(160,80)) and ((320,160)).
- Regular-(F) cubic/Gauss shadow initializer: the same two domains. The shadow changes only the initial radial (B) quadrature and then uses the unchanged production evolution kernel.

No evolution failures occurred. The recorded cycle-event and handoff counts remained zero. Production physics, defaults, projection, and gauge setup were not changed. The run holds radial spacing fixed; it is a domain-size test, **not** a spatial-resolution convergence test.

## Finding 1 — the global maximum changes identity

The full profiles resolve the ambiguity left by the half-CFL summaries. The maximum of (|H|) over off-centre cells can jump from an inner peak to a separate, growing outer peak. The crossover depends strongly on the initializer:

| Domain and initializer | First sampled time outer peak exceeds inner peak |
|---|---:|
| (r_{max}=80), production | (t=8.74875) |
| (r_{max}=160), production | (t=9.00000) |
| (r_{max}=80), regular-(F) shadow | (t=6.25125) |
| (r_{max}=160), regular-(F) shadow | (t=6.25125) |

Thus the jump in the **global argmax radius** is an inner/outer peak switch, not a reliable trajectory measurement by itself.

Separately, the outer (|H|) ridge has a coherent outward-moving maximum in the sampled profiles. At late times its peak radius progresses stepwise from approximately (r=6.25) at (t=8), to (r=7.25) around (t=9), and to (r=9.25)–(9.75) at (t=12). The two initializer variants closely agree on this outer peak within a fixed domain by (t\approx9). For example, at (t=12):

| Domain | Production outer peak | Shadow outer peak |
|---|---:|---:|
| (r_{max}=80) | (2.8291823\times10^{-4}) at (r=9.75) | (2.8281186\times10^{-4}) at (r=9.75) |
| (r_{max}=160) | (2.5979611\times10^{-4}) at (r=9.25) | (2.5969952\times10^{-4}) at (r=9.25) |

That is strong evidence for a shared **outer constraint-residual feature** in these four trajectories. It is not proof of a physical front, matter/D propagation, or domain-independent dynamics.

## Finding 2 — the inner trajectory remains initializer-dependent

The near-equality of the late outer peak does not mean the full trajectories converge.

At (t=12), for (r_{max}=80), the production inner-sector maximum is approximately (1.008\times10^{-4}) at (r=1.25); the shadow inner-sector maximum is (4.334\times10^{-5}) at (r=2.75). The central vendor Hamiltonian residual is (+1.0933\times10^{-3}) in production versus (-7.5884\times10^{-6}) in the shadow. The (r_{max}=160) pair shows the same split.

The shadow's outer peak overtakes its inner peak at the first sampled time (t=6.25125); production's inner peak persists as the global off-centre maximum much later. The initializer therefore remains material to the inner residual and to the timing of the global-maximum switch, even while the outer ridge later becomes similar.

## Finding 3 — domain size changes the lapse profile substantially

The domain effect is not confined to the first-cell lapse. On the common radial grid, the maximum absolute difference between the two production lapse profiles grows from (5.67\times10^{-4}) at (t=0), to (0.04741) at (t=8), (0.08259) at (t=10), and (0.11614) at (t=12). At (t=12), the largest difference is at (r=39.75): (\alpha_{80}=0.9538122) and (\alpha_{160}=0.8376731). The RMS lapse-profile difference over the shared grid at that time is approximately (0.07828). The shadow pair gives essentially the same domain-dependent lapse profile; initializer changes alter the lapse only at the few-(10^{-6}) level or less in these samples.

The first-cell lapse also differs with domain. In production it is (0.0165559) at (r_{max}=80) and (0.0172775) at (r_{max}=160) at (t=12), a (4.36\%) relative difference. Its absolute domain difference peaks earlier, at about (-0.00653) near (t=6.25125).

Meanwhile, at (t=12), the off-centre (|H|) maximum is about (8.17\%) smaller on the larger domain and lies one cell ((0.5) in radius) inward. Across the late-time samples the outer-peak amplitude difference reaches about (8.43\%). These are material domain effects at fixed (\Delta r), not evidence of spatial-resolution error.

The lapse results are consistent with an important global-domain/gauge dependency, but this campaign alone does not prove a unique causal mechanism. Do not treat the outer ridge as domain-independent until this dependence is better bounded.

## Health and interpretation gates

- Every case reaches (t=12) with 3,200 accepted steps and no recorded numerical failure.
- Metric variables remain finite and positive. Across the sampled trajectories, minima stay approximately (a\ge0.952), (b\ge0.99739), and (X\ge0.9353).
- The production connection-constraint maximum remains larger than the shadow's, even though the late outer (H) peaks agree closely; the latter is not a general statement that the initial-data defect has been removed.
- No cycle event, handoff, bounce, or turnaround is registered in these runs.

## Decision and next diagnostic

**No physics patch.** The present result is: the global-max jump is a switch between distinct inner and outer peaks; a coherent outward-moving maximum exists in the outer Hamiltonian-residual profile; its physical meaning is unresolved; and its late amplitude/radius and the full lapse profile depend materially on domain size.

Next, add an intermediate domain ((N,r_{max})=(240,120)) at the same (\Delta r=0.5), with production and regular-(F) shadow initializers, identical physics/gauge/configuration, quarter-time profiles, and the same fail-closed numerical checks. This tests whether the observed domain effect varies smoothly between 80 and 160; it does not substitute for a spatial-resolution study. Keep this a diagnostic-only campaign.
