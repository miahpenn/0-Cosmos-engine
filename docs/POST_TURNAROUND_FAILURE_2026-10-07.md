# 0-Cosmos V5.5 post-turnaround failure adjudication — 2026-10-07

## Campaign
True-CMC long post-turnaround gate:
- N=80,120,160
- R_max=80
- CFL=0.06
- target t=80
- native V55TrueCMCPIRKKernel
- no bounce/reset/branch flip/lapse floor/physical stop

## Result
All three resolutions failed in the same late post-turnaround window:
- N=80: t=46.920000000000265
- N=120: t=46.91999999999918
- N=160: t=46.98000000000102

The immediate exception was a radiation conservative-state realizability failure in the outermost cell:
E < |S| by approximately 0.03–0.14 percent. This is not accepted as a physical term or repaired by clipping.

## Cross-resolution geometry signature
At failure the runs simultaneously show:
- one turnaround and one handoff;
- no re-expansion;
- two trapped roots, with the outer root near r=65–68;
- strong radial metric stretching concentrated around r~70–74;
- rapidly growing Hamiltonian residuals in the outer region;
- CMC residuals much larger than the clean t=40 gate;
- radiation outer characteristic flow is fully outgoing.

The radiation-off ablation independently reached later but still failed through geometric/metric breakdown (nonpositive conformal metric ratio) at t=48.84, with very large outer-boundary geometry derivatives and Hamiltonian residuals.

## Domain extension / radiation-off evidence
The matched-spacing R=120 and R=160 campaigns established that enlarging the domain moves the failure later but does not remove it. The R=160 radiation-off run failed at t=72.78 and r=119.25, well inside R_max=160, with invalid spatial volume element from a<0. This keeps radiation from being the fundamental continuation blocker and shifts attention to strong-field geometry/CMC continuation.

## Boundary-control discriminator correction
The completed R160 no-light-CPBC experiment must not be treated as a valid light-boundary discriminator.

Reason: V55TrueCMCPIRKKernel.step() does not invoke the outer-light-boundary hook. That hook exists in V55ProductionKernel.step(), but the True-CMC production override bypasses it. Therefore the no-light subclass used by that campaign changed no active operation; its result is effectively another True-CMC geometry continuation run.

The run remains useful as an independent reproduction of the same R160 radiation-off geometry failure, but it provides no causal evidence about the pinned light-constraint boundary operator.

## Current next gate
The next experiment is the actual coupled V5.5 machine plus an explicit CMC-operator diagnostic, not another boundary-removal test.

The coupled diagnostic:
- uses V55ProductionKernel, the current unified production state graph;
- keeps the full local S/D + COSMOS + DM + baryon + radiation stress-energy system active;
- preserves the locked DM/phi exchange;
- varies only the initial D amplitude across 1e-10, 1e-9, 1e-8;
- records D and S phase, COSMOS state, proper time, H_eff, lapse, Misner-Sharp/worldtube quantities, stress-energy flux/work, constraints, and trapped-root counts;
- introduces no reciprocal phenomenological source.

The independent CMC probe measures the same CMC elliptic operator's singular-value conditioning, solve residual, lapse profile, and kdot alongside H_eff through the post-turnaround window.

This is the correct reconnection point: measure whether the local↔COSMOS dynamics alter the continuation signature while simultaneously determining whether the CMC operator itself is becoming ill-conditioned.

## Standing rule
No clipping of E or S, no radiation floor, no fitted damping, no manufactured bounce/continuation rule, and no promotion of an experimental gauge without a controlled discriminator. Repair only the numerical layer actually isolated by the campaign evidence.
