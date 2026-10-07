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

## Adjudication
The radiation inversion is therefore a valid numerical witness, but not the sole fundamental blocker. The common failure window and radiation-off continuation identify the finite-radius strong-field/outer-worldtube geometry as the next primary diagnostic target.

The next experiment is a domain-extension test using the existing True-CMC R120 matched-spacing workflow:
- R_max=120
- N=120,180,240
- CFL=0.06
- long horizon
- require H_eff=0

This changes only the computational domain/resolution relationship, not the physical equations.

## Standing rule
No clipping of E or S, no radiation floor, no fitted damping, and no manufactured bounce/continuation rule. Promote a radiation repair only if the discretization defect is isolated independently of the geometric boundary failure.
