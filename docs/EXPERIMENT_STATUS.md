# Experimental status — 0star central clock

**Branch:** `0star-central-clock`  
**Production:** `main`  
**PR:** #1, draft

## Completed evidence

- D-dependent interior evolution is reproducible.
- Rephasing by `D_active` produces strong trajectory collapse in D, momentum, and geometry proxies.
- The early-phase long campaign reached coordinate time 40 but accumulated only about 0.78 proper-time units after synchronization because the central lapse collapsed.
- The pre-repair resolution campaign showed an inward-moving Hamiltonian-constraint maximum tied to the first retained radial cell.
- CMC elliptic residuals were small while the full Hamiltonian constraint was not, so the old lapse-collapse result cannot yet be treated as clean physical evidence.

## Repair now under validation

The production kernel reapplies the existing spherical-center algebraic regularity projection during the evolution stages. This is a numerical repair, not a new physical law.

## Required gates

- Smoke/regression suite.
- N=80/120/160 center-regularity resolution audit.
- N=160 repaired Dhalf/Dbase/Ddouble phase campaign.
- CFL/domain audit.
- Only after these pass: extended turnaround/bounce campaign.

## Interpretation rule

A persistent lapse collapse after the repair must be tested against regional/normalized constraint diagnostics and resolution/CFL/domain convergence before being called a physical clock-stalling or singular behavior.
