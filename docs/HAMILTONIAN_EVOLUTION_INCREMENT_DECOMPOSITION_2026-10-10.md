# Accepted-step Hamiltonian increment decomposition — 2026-10-10

**Branch:** `diag/hamiltonian-evolution-increment-decomposition`  
**Base:** `0star-central-clock`  
**Scope:** diagnostic-only; no evolution equation, source, gauge choice, initializer default, or production file is changed  
**Status:** open numerical attribution question; not admitted as a physical result

## Purpose

The previous initial-slice audit and Lumen's local discrete-consistent initial-data prototype point to two potentially distinct effects: an initial Hamiltonian mismatch, and a residual that reappears during evolution after the proposed initializer removes the initial mismatch.

This branch isolates the second question first. It runs a short, explicitly bounded baseline trajectory and measures the Hamiltonian residual before and after every accepted step. Each observed increment is attributed symmetrically to six state groups:

1. Spatial metric: (a,b,X).
2. Lapse/shift gauge state: (alpha,eta,B).
3. Trace-free curvature variable: (A_a).
4. Trace curvature: (K).
5. Connection variable: (Lambda).
6. Matter: all scalar fields and all conservative fluid species.

The Hamiltonian convention is the production Step 1 convention:
`H = R - (A_a^2 + 2 A_b^2) + (2/3) K^2 - 16 pi rho`, with `A_b = -A_a/2`. The implementation obtains the residual through the pinned vendor constraint operator and subtracts the total matter density using the existing stress-energy assembly.

## Why the attribution is symmetric

The Hamiltonian is nonlinear. A sequential “change metric, then change K, then matter” decomposition would assign interaction terms according to the chosen ordering. Instead, this diagnostic evaluates all (2^6=64) pre/post hybrid states for a step and computes the Shapley attribution for each group. Thus:

- component contributions sum to the measured (H_{n+1}-H_n), to floating-point closure;
- nonlinear interactions are distributed symmetrically over the participating groups;
- the reported closure error checks bookkeeping only. It is not a physics threshold or constraint-admission gate.

The hybrid states are copies. They are never evolved and are never written back to the trajectory.

## Registered baseline settings

The workflow's initial run uses:
- (N=40), (r_{max}=40), (Delta r=1);
- amplitude (S_0=0.01), width (7);
- (D_0=10^{-10}), radiation enabled;
- CFL (=0.0075), requested final coordinate time (t=3).

At this spacing the nominal step is (0.0075), giving about 400 accepted steps. The small trajectory is intended to resolve the initial residual increment, not to test a turnaround or a long-time outcome. The report includes the exact source commit, initial/final residual summaries, per-step changes at cells 0–4, all six component attributions, and full-grid component-sum closure.

## Critical limitation: candidate initializer not yet included

The branch on GitHub does **not** contain Lumen's reported local-only `ac7bdf6` discrete-consistent initializer candidate. The first run therefore characterizes the currently checked-in baseline initializer only. It is not the requested A/B comparison and it must not be described as confirming the prototype's numbers. Once the candidate's exact code is published or supplied, it should be added as a separate opt-in initializer and the baseline/candidate comparison rerun under identical settings. No production default should be changed.

## Interpretation guardrails

- Reappearance of (H) does not on its own identify an incorrect field equation; it identifies the component contributions requiring further operator-level review.
- (Lambda) is a connection variable, not an independent physical source.
- A large Shapley contribution is an attribution of a discrete constraint change under the chosen grouping; it is not proof that the attributed equation is wrong.
- No bounce, turnaround, completed cycle, or physical validity claim follows.
- This branch does not alter `main` or `0star-central-clock`; it does not alter production equations, add damping, add floors, or apply a residual repair.

## Reproduction

```bash
python -m pytest -q tests/test_hamiltonian_evolution_increment_decomposition.py
python -m engine.hamiltonian_evolution_increment_decomposition
```

The workflow uploads `runs/hamiltonian-evolution-increment-decomposition/report.json`. Retain it with the exact source commit. A completed run with finite outputs and closed component bookkeeping is only a successful diagnostic execution, not trajectory admission.
