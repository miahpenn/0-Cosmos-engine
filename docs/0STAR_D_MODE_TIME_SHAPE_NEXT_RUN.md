# 0* D-mode time-shape investigation — next-run proposal

## Current result

The completed D-mode interior causality campaign (Actions run 37512951152) ran four N=160, Rmax=160 trajectories with the same low coordinate-CFL control and only the initial inverted-D amplitude varied: D0=0, Dhalf=5e-11, Dbase=1e-10, Ddouble=2e-10.

All four completed t=22.5. No turnaround or reexpansion occurred.

| case | tau(t=22.5) | H_center | lapse_min | max |alpha_r| | D_center | D_active_center |
|---|---:|---:|---:|---:|---:|---:|
| D0 | 22.3974 | +9.34e-4 | 0.9900 | 0.00127 | 0 | 0 |
| Dhalf | 21.4338 | -1.001e-2 | 0.2671 | 0.0376 | 0.03735 | 0.00403 |
| Dbase | 20.8697 | -1.401e-2 | 0.1574 | 0.0408 | 0.04338 | 0.00546 |
| Ddouble | 20.2487 | -1.718e-2 | 0.0954 | 0.0422 | 0.04768 | 0.00660 |

The alpha-gradient maximum moves outward with increasing D amplitude, approximately r=9.5, 11.5, 12.5 for Dhalf, Dbase, and Ddouble.

## Interpretation status

The data support a D-amplitude-dependent interior dynamical time-shape. Increasing the initial inverted-D amplitude advances the development of central contraction and lapse deformation, while the D=0 control does not show the same transition.

This is not yet proof of one-way causality. D and geometry are coupled. The next experiment must resolve source ordering and state-space alignment rather than merely compare endpoints.

No physical term, fitted coefficient, boundary condition, or bounce rule is justified by this result.

## Next-run objective

Determine whether clock deformation is organized by the instantaneous D state/source rather than coordinate time or initial amplitude, and identify the existing stress-energy/geometry channel through which D feeds the CMC response.

### 1. Keep physics unchanged
Use N=160, Rmax=160, the same pinned vendor kernel, the same CMC gauge, the same low coordinate-CFL rule, and the same four D amplitudes. No new source, damping, floor, boundary condition, or coefficient.

### 2. Increase temporal resolution
Record diagnostics at about Delta-t=0.25 through t=22.5. The goal is to resolve onset ordering.

### 3. Add source decomposition
At the center and at the instantaneous maximum-lapse-gradient location, record D, PD, D growth rate, existing active-D source 2 PD^2 + D^2, S-mode source, COSMOS scalar contribution, existing matter contribution, exact sector contributions to L3(K) where the existing adapter permits them, CMC Kdot, lapse/alpha_r, central H, H_eff, and central proper-time rate.

These are diagnostics of existing equations only; no reconstructed source is to be invented.

### 4. Perform two alignments
Time alignment: compare observables against coordinate t.
State alignment: compare the same observables against instantaneous D_active rather than t.
The decisive question is whether different initial-D runs approximately collapse onto a common relation D_active -> geometric source -> alpha response despite different coordinate-time histories.

### 5. Establish ordering
For each run determine the onset ordering: D growth -> D stress -> matter/geometry L3 source -> CMC Kdot -> alpha/alpha_r -> clock-rate change.
Use thresholds only as post-run analysis markers, never as forces in the engine. Report ordering in both coordinate time and D-state.

### 6. Keep D0 mandatory
The D0 trajectory remains the control for whether the sequence exists without the inverted mode.

### 7. Follow with resolution testing
If D-state alignment is strong at N=160, repeat only the most discriminating cases at N=80 and N=120. The goal is convergence of the mechanism, not merely convergence of a failure time.

## Decision tree
1. D-state collapse + D-source precedes clock response: strong support for a D-driven dynamical time-shape; trace the exact archived field/geometry equation carrying the source.
2. D-state collapse but ordering is bidirectional: treat it as coupled feedback; map the loop before changing anything.
3. No D-state collapse: initial amplitude is controlling the trajectory more directly; investigate intermediate geometry/matter variables.
4. D0 develops the same structure: D is not the unique driver; return to the full source decomposition.
5. Strong resolution dependence: stop physical interpretation and isolate the numerical layer.

## 0* gate
0* remains the diagnostic/scoring layer, not a construction rule. After the dynamical mechanism is measured, apply t: turnaround {H=0, dot H<0}; b: bounce {H=0, dot H>0}; w >= 1 anti-BKL during contraction; w > 1 stronger smoothing; G as the separate global timescale condition.

The immediate question is narrower: What interior dynamics actually produces the observed time-shape?