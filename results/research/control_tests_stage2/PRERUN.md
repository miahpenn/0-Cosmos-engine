# Stage 2 pre-run note (written 2026-11-25T08:43:49Z, before launch)

Go: Jeremiah, 08:42:29 UTC, relayed: "Go on stage 2" (answering the coordinator's question after A23). Code: main dd746f6e (#607, #608, #609). Overrides V0 = 0, SD_QUARTIC_LAMBDA = 0.
Arms: T2-3 (1e-2, .05, dt0, r8, Cdot per step, R-stage), T2-4 (3e-2, .025, dt0, r8, R-stage), T2-6 (5e-3, .05, dt0, r8), T2-E (1e-2, .05, dt0, r12); t_end 2.1.
Gates as predeclared (design, A1-A5, A20, A21, A22, A23), none changed. Order: input checks (test2 inputs + mp_readback); any failure stops the arms; then arms; then cdot_floor, mp_readback, score2.
An earlier input phase (08:33, under A22) stopped on the A22 item 3 P_S RHS term; A23 replaced that term. No arm ran before this note.
Stage 3 (T2-1, T2-2, T2-5) is not covered by this go.
