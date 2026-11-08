# Issue #4 Gate C, accepted class: locked full-Bianchi residual under refinement (2026-11-08)

Read-only. `scientific_admission=false`, `acceptance_oracle_satisfied=false`. No production code,
state, tolerance or evaluator logic changed; zero production steps at every resolution.

## What was run

`tools/issue4_gate_c_accepted_convergence.py` calls `tools/issue4_full_bianchi_evaluator.py`
(`evaluate_accepted_fixture`, unchanged) on the archived accepted t = 0.60 states at N = 60, 120, 240
(artifacts 12615136565, 12615555594, 12615550675 from run 40236402303; archive SHA-256 checked by the
workflow and by the fixture loader). The evaluator's accepted reference regression exists only for
N = 60; it ran and **passed** there, and is recorded as not executed for N = 120 and 240.

- Workflow run: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/42505508956
  (head `0c8c550848d3f34b0c99e674a38bb6660f24a6ee`; diagnostic source `4aec88d0`, helper baseline
  `7f47669a`, vendor `d6052d60`; Python 3.12.15, numpy 2.5.4, scipy 1.18.1).
- Artifact `issue4-gate-c-accepted-convergence-42505508956` (id 13627112984), published digest
  `sha256:5786032c5dad70643e89a3bc2eb6a224b20d0703aa42a326654135156f152302`; the JSON inside has
  SHA-256 `535b0aa46d3beed58bec6cdf2df39aeaf74f4d021526bc0db1a46a94c2658c75` (printed in the run log).
- `gate_c_accepted_convergence_from_run_log.json` here is the same report, printed compactly between
  the `GATE_C_REPORT_BEGIN/END` markers of that run's log (SHA-256
  `a4950907b0e59d7749fd8a4bb3549194c03a614c1111eb59b015bb11d5cb5ab0`). The artifact could not be
  downloaded from the session that wrote this note.

The three states are separate evolutions to t = 0.60 with dt halved at each doubling, so this is a
refinement sequence of states, not a spatial-only test on one state.

## Results (ε = 1e-5; ε sensitivity ≤ 1.3e-12 at every N)

Residual = measured M_r tangent − locked predictor.

| N | pointwise max, r ≤ 10 | FV max, r ≤ 10 | pointwise r = 2 | pointwise r = 8 | annulus 5–10 RMS (pw) | measured tangent max |
|---|---|---|---|---|---|---|
| 60 | 1.84e-5 (r = 0.67) | 9.64e-6 | +1.40e-5 | −1.35e-6 | 1.45e-6 | 9.93e-6 |
| 120 | 1.53e-6 (r = 0.33) | 1.28e-6 | +9.16e-7 | −8.24e-8 | 9.96e-8 | 7.24e-7 |
| 240 | 3.99e-7 (r = 0.17) | 4.08e-8 | +6.03e-8 | −5.11e-9 | 6.60e-9 | 8.79e-8 |

Observed orders (60→120, 120→240):

- Pointwise: r = 2: 3.94, 3.93; r = 3: 3.82, 3.97; r = 8: 4.03, 4.01; annulus RMS: 3.86, 3.92.
- Finite volume: r = 2: 4.01, 3.88; r = 3: 3.75, 3.98; r = 8: 4.01, 4.01; annulus RMS: 3.84, 3.94.
- The pointwise maximum sits in the first cells and converges at 3.59, then 1.94; centre cells are
  reported separately and not used for the fixed-radius reading.

## Reading

At every fixed radius and in the 5–10 annulus, in both divergence channels, the residual falls at
about fourth order, the stencil order. There is no plateau, so on the accepted states the locked
projected-Bianchi predictor leaves no residual that survives refinement.

**Limit of this test.** On these states every operand is itself truncation-sized and falls at about
fourth order: max |H_total| 2.4e-5 → 2.2e-6 → 8.4e-8, |M_total| 6.3e-6 → 4.4e-7 → 3.7e-8,
|C^r| 3.9e-6 → 2.8e-7 → 2.2e-8, and the measured tangent itself 9.9e-6 → 7.2e-7 → 8.8e-8 (the
residual is 0.76–0.94 of the local tangent RMS in the annulus at every N). So this rules out a missing
term that is nonzero on a constraint-satisfying, conserving state, but it tests terms proportional to
constraint or conservation violations only weakly, because those vanish with the operands. The
geometric constraint-proportional terms were tested with O(1) violations on the vacuum manufactured
state (`results/issue4-hamiltonian-slot-decomposition-2026-11-08/`). The matter-side counterparts
(the ΔB_sub conversion terms and −8πα·D_tot with nonzero D_tot) still need a state with nonzero matter
constraint violation under the actual coupled RHS.
