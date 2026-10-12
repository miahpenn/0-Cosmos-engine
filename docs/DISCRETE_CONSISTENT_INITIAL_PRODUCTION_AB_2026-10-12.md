# Discrete-consistent initial metric — matched production A/B result

**Date:** 2026-10-12  
**Workflow:** [38284401548](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38284401548)  
**Branch:** \`diag/discrete-consistent-initial-production-test-2026-10-12\`  
**Artifact:** \`discrete-consistent-initial-production-ab-t2p16\` (ID \`11762503161\`)

## What was tested

The baseline and candidate both used the same strict-admission production kernel, source equations, matter state, CMC gauge, PIRK stepping, outer domain and CFL. The candidate changed only the *initial metric construction*: its opt-in Newton solve varied the radial metric variable \(B\) until the discrete Hamiltonian residual used by the production constraint diagnostic was zero to tolerance. The baseline continued to use the existing cumulative-integral initializer.

Both arms ran at N=40,60,80,120, \(r_{\max}=80\), CFL \(=0.03\), radiation on, through \(t=2.16\). Each requested run reached its horizon and passed the strict execution gate. The test suite and explicit initializer regression passed. The baseline histories/checkpoints from this matched run agree with the archived exact-source t=2.16 profiles to roundoff (maximum norm differences typically \(10^{-14}\) or smaller).

No evolution equation, stress source, physical coefficient, gauge equation, or production default was changed.

## Initial-slice result

All four candidate solves converged in **two Newton iterations**. Maximum discrete Hamiltonian residual at \(t=0\) was:

| Resolution | Baseline max \(|H|\) | Candidate max \(|H|\) | Maximum Jacobian condition number |
|---:|---:|---:|---:|
| N=40 | \(3.77\times10^{-3}\) | \(4.27\times10^{-15}\) | \(2.73\times10^3\) |
| N=60 | \(3.93\times10^{-3}\) | \(3.34\times10^{-15}\) | \(6.15\times10^3\) |
| N=80 | \(3.99\times10^{-3}\) | \(2.31\times10^{-15}\) | \(1.09\times10^4\) |
| N=120 | \(4.03\times10^{-3}\) | \(6.11\times10^{-15}\) | \(2.46\times10^4\) |

The candidate retains \(K\), \(A_a\), all scalar fields and every conservative matter field exactly from the baseline initializer. The initial momentum residual is unchanged at roundoff; \(C_\Lambda\) remains at roundoff because Lambda is reconstructed from the metric; and \(a b^2=1\) remains at machine precision. Thus this is a genuinely constraint-satisfied discrete initial slice for H/M/C—not an H-only tweak that immediately injects an initial momentum error.

## Results at t=2.16

The large full-domain central residuals improve substantially at N=60 and N=80, but **the result is not uniform under refinement**. The finest tested grid develops a larger first-cell H and M residual than its own baseline.

### Full-domain maximum norms

| N | \(\max |H|\), baseline → candidate | \(\max |M|\), baseline → candidate | \(\max |C_\Lambda|\), baseline → candidate |
|---:|---:|---:|---:|
| 40 | \(7.41\!\times10^{-4}\to7.26\!\times10^{-4}\) | \(1.78\!\times10^{-3}\to9.17\!\times10^{-5}\) | \(1.64\!\times10^{-3}\to2.80\!\times10^{-4}\) |
| 60 | \(9.27\!\times10^{-3}\to9.93\!\times10^{-4}\) | \(1.11\!\times10^{-3}\to2.39\!\times10^{-4}\) | \(1.99\!\times10^{-3}\to2.67\!\times10^{-4}\) |
| 80 | \(6.97\!\times10^{-3}\to1.60\!\times10^{-3}\) | \(2.54\!\times10^{-3}\to5.04\!\times10^{-4}\) | \(9.73\!\times10^{-4}\to5.96\!\times10^{-4}\) |
| 120 | \(4.74\!\times10^{-3}\to9.85\!\times10^{-3}\) | \(2.43\!\times10^{-3}\to3.66\!\times10^{-3}\) | \(1.20\!\times10^{-3}\to1.40\!\times10^{-3}\) |

The coarse N=40 constraint peak is itself at the centre, and N=40 contains no cell satisfying both \(i\ge2\) and \(r<5\). Use the admitted-inner/fixed-radius measures below before interpreting any resolution trend.

### Admitted inner region (\(i\ge2,\ r<5\))

| N | \(\max |H|\), baseline → candidate | \(\max |M|\), baseline → candidate | \(\max |C_\Lambda|\), baseline → candidate |
|---:|---:|---:|---:|
| N=40 | No cells in this mask | No cells in this mask | No cells in this mask |
| N=60 | \(7.51\times10^{-4}\to6.19\times10^{-4}\) | \(5.00\times10^{-5}\to3.40\times10^{-5}\) | \(4.74\times10^{-4}\to2.61\times10^{-4}\) |
| N=80 | \(1.07\times10^{-3}\to1.02\times10^{-3}\) | \(3.04\times10^{-4}\to2.27\times10^{-4}\) | \(3.24\times10^{-4}\to2.38\times10^{-4}\) |
| N=120 | \(1.54\times10^{-3}\to1.67\times10^{-3}\) | \(9.36\times10^{-4}\to9.09\times10^{-4}\) | \(2.56\times10^{-4}\to1.84\times10^{-4}\) |

The fixed exterior \(r\ge5\) improves at every resolution. For instance, N=120 exterior H falls from \(5.46\times10^{-6}\) to \(1.50\times10^{-7}\), while N=80 falls from \(8.92\times10^{-6}\) to \(1.09\times10^{-6}\). The discrepancy remains concentrated at the centre/first retained cells.

The proper-time changes at t=2.16 are small (about \(-2.10\times10^{-5}\), \(-5.92\times10^{-6}\), \(-2.42\times10^{-6}\), and \(-5.99\times10^{-7}\) for N=40/60/80/120). Misner–Sharp current RMS also improves at N=40/60/80 but worsens at N=120, so the balance witness is not uniformly repaired either.

## Interpretation

This experiment separates two facts:

1. **The old initializer contributes a genuine discrete H-constraint defect.** The candidate removes it at t=0 without spoiling the initial M or connection constraint.
2. **That initial defect is not the whole production failure.** Starting from the constraint-satisfied slice, the existing evolution regenerates central H/M/C residuals by t=2.16. The residual response is strongly resolution-dependent: N=60/80 improves, whereas at N=120 the central H/M maxima are worse than baseline. That is consistent with a central discrete compatibility layer that an H-only initial solve does not cure.

The candidate's \(t=0\) H closure is built into its Newton solve, so it is a solver check, not independent validation. The short evolution and independent same-source baseline comparison are the meaningful tests. No long campaign was run.

## Decision and next discriminator

**Do not replace the production default initializer, merge this candidate, or authorize t=160 from this evidence.** Keep the solver opt-in.

The next test is a one-step comparison from the exact saved baseline and candidate \(t=0\) states at all four resolutions. Record \(\Delta H/\Delta t\), \(\Delta M/\Delta t\), and \(\Delta C_\Lambda/\Delta t\) after one unmodified production step, including the first five cells and fixed physical radii. This will locate how much of the regenerated defect is already in the first discrete update versus emerging over the first 2.16 units. If the constraint-satisfied initial data immediately develops the same signed defect, the next repair must be in evolution/stage regularity or derivative compatibility rather than further initial B fitting.

## Reference context

Constraint propagation studies emphasize that satisfying the constraints initially is necessary but does not guarantee a discrete evolution preserves them; the continuum BSSN constraint subsystem and its zero-speed modes make the evolution/constraint interface important. See the Living Reviews treatment of Einstein initial-boundary value problems and BSSN constraint propagation ([link](https://link.springer.com/article/10.12942/lrr-2012-9)) and Bernuzzi & Hilditch's comparison of BSSNOK and Z4 constraint propagation in spherical symmetry ([Phys. Rev. D 81, 084003 (2010)](https://doi.org/10.1103/PhysRevD.81.084003)). Neither paper validates this particular coupled CMC/PIRK kernel; they support the decision to test constraint propagation directly rather than infer success from initial residual closure.
