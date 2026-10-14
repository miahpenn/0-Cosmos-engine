# GEAR Cycle Operator v1 — high-precision reconstruction

**Status:** research-only numerical utility; not a physical bounce mechanism  
**Source branch:** \`research/gear-cycle-operator-v1-2026-10-14\`  
**Archive references:** \`gear_cyclic_27_cycle_operator_reconstruction_2026-10-02.py\` and \`GEAR_cyclic_operator_integration_2026-10-02.md\`.

## Operator

One reconstructed pole swap acts on \((a,H,\dot H)\) as

\[
(a,H,\dot H)\mapsto \left(\frac{a_T^2}{a},-H,-\dot H\right).
\]

A full numerical cycle is defined as two pole swaps. The archived regression uses decimal-string seeds, 80-digit arithmetic, 27 full cycles, and therefore 54 swaps. The current implementation uses \`mpmath.workdps\` to isolate the requested precision instead of permanently mutating the caller's global precision.

Run the isolated contract tests with:

\`\`\`bash
pytest -q tests/test_gear_cycle_operator_v1.py
\`\`\`

## Evidence and limits

The preserved seed is \(a_T=9.80859460229\), \(a_0=3.2718281828459045235360287471352662497757\), \(H_0=0.123456789012345678901234567890123456789\), and \(\dot H_0=-0.1131778\). At the declared 80-digit arithmetic, the archived seed returns exactly after 27 cycles / 54 swaps.

This verifies the reconstructed map's arithmetic closure for the stated input, operation order, and precision. It does **not** establish:
- a first-principles physical bounce law;
- that a cosmological evolution physically completed 27 cycles;
- conservation or stability of an evolving COSMOS trajectory;
- equivalence to any unrecovered original notebook implementation.

The original archive labels this source a **reconstruction, not recovered original notebook source**. Keep the operator isolated from \`engine/cosmos.py\`, the production timestepper, and the physical bounce gate. Any new physical interpretation requires its own equation-derived argument and independent validation.
