# Matched baseline decomposition trigger

Diagnostic-only workflow trigger record.

The workflow must run `tests/test_discrete_consistent_initial.py` before the matched N=40 baseline Hamiltonian-drift decomposition. The baseline mode must use `V55ProductionKernel.initialize(...)` unchanged. Settings are pinned in `engine/discrete_consistent_initial_decomposition.py`.

Admission: NOT_ADMITTED_DIAGNOSTIC_ONLY. No production changes or physical interpretation are authorized by this run.
