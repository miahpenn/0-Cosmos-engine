# The floor run: GEAR-91 on, hand floor off, φ feels D's floor (research only)

Approved Nov 15. One-budget lattice cell (λ = 5000, V0 and densities × c², hand floor off, D's ejected vacuum the only
floor in geometry). New: β_D = −0.04 through the GEAR-91 switch (β_S = 0), so φ also feels the floor (8π ∂V_SD/∂φ).
Three short runs (r_max 24, N 240, to t = 70): background cells (no matter) at β_D = 0 and −0.04, and the pocket
cell x_c 2.0 at β_D = −0.04 (its β_D = 0 twin is the stored u_x20_cold).
`research/lattice_cell_run_coupled_v1.py` (driver), `research/phi_floor_v1.py` → `phi_floor.txt`, `held_coupled.txt`.

## Summary

- **The floor is counted once, and the machine matches the derivation.** The coupled background cell follows the
  derived background to ≤ 1.5e-4 in φ, Π and H up to t = 50: V(φ) = V0 e^(−φ), vacuum −8π e^(0.16φ)/(2λ) counted
  once, and φ pushed by its slope. D tracks the moving sea rest D* e^(−β_D φ) to 1e-3. Control: the uncoupled cell
  matches its prediction to 1e-7.
- **The floor moves with φ, and that ends the research universe sooner.** As φ rolls, D's step deepens (e^(0.16φ)). φ
  rolls 3% further by t = 20 and 10% by t = 40. It crosses its node earlier (1.985 instead of 2.303), and the growing
  negative vacuum turns the cell around at t ≈ 42 instead of ≈ 55. It crunches at t ≈ 58.8, where the runs stop: the
  geometry becomes singular, for the background and the pocket cell alike. (Read; the derivation reproduces it to
  1.3e-2 at t = 58.5.)
- **The reset itself barely changes.** The held step (measured against the φ-dependent sea rest) peaks at 0.23 at
  t = 3.5, against 0.20 at t = 4.5 uncoupled, and tracks the uncoupled values within about 10% afterwards (t = 6–40).
- **φ now feels matter through the floor.** Where D is held near 0 there is no floor push (∂V_SD/∂φ ∝ D²). φ in the
  pocket falls behind the sea by an extra −0.006 to −0.02 against the uncoupled cell (pocket − sea: −0.029 → −0.041,
  against −0.023 → −0.023). The effect is small and has the derived sign.
- **No second reset before the crunch.** As the cell collapses, the core's x rises from 0.03 (t = 40) to 0.39 (t = 58).
  D in the core falls to 0.87 of the sea's value. Both stay far from the hold band (x = 1.27) when the code stops. Chat's
  gravity-driven second reset would need the run to get closer to the crunch than this code can. Here it would be
  cosmos-wide (the whole cell collapsing), not a local stone. **Open.**

## Background, cell − derived (β_D = −0.04)

| t | φ − φ0 | dφ | dΠ/Π | H | dH/H | D/D* | e^(−β_D φ) |
|---|---|---|---|---|---|---|---|
| 10 | 0.6584 | −1e-4 | −2e-4 | 0.0608 | +2e-4 | 1.0333 | 1.0341 |
| 30 | 1.8474 | −1e-4 | +8e-5 | 0.0224 | −8e-5 | 1.0843 | 1.0844 |
| 50 | 3.2072 | +1e-4 | +4e-4 | −0.0261 | +8e-4 | 1.1450 | 1.1451 |

## What this says about the floor question (Read)

- The machine pays the floor from D's step through φ, with no hand floor and no double count. That part of the floor
  question is now done in the research cell.
- The derived consequence is a **moving floor**. Today's dark energy then depends on how far φ has rolled. A
  universe whose φ reaches its node turns around. At production's φ_today = 0.18 the node (1.985) is far off. When φ
  gets there is a lane question (production H, the floor at λ = 1.5e6), not a cell question. **Next, if wanted:** run
  the production lane with the coupled floor and find when φ reaches 1.985, against the hand floor's 2.303.
- A caution for any held-step read-out with the coupling on: the fixed-vacuum measure of the handoff note reads the
  φ-driven shift of D_s as held step (0.24 at t = 12, 0.8 at t = 24 here). It must use the sea rest e^(−2β_D φ), as
  `held_coupled()` does.

Unpaid inputs: λ = 5000 (one-budget c² = 304.9); β_D = −0.04 (borrowed); β_S = 0; φ_today = 0.179 (imported); x_c.
