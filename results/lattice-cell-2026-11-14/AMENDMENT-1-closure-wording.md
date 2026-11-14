# Amendment 1 to the lattice-cell README (14 Nov 2026): what closure (ii) is

The README (Task 1, W4, and the map-row-7 proposal) writes closure (ii) as "K̇ = ⟨Q⟩". Stated properly:

- dK/dt = ⟨αQ⟩ (proper-volume average) holds for **any** closure. It follows from the lapse equation and the zero
  lapse gradient at the wall. **Derived**; read to 2e-6 in every stored run.
- Closure (ii) is the choice of clock in which ⟨αQ⟩ = ⟨Q⟩: the Q-weighted mean lapse is 1.
- Closure (i) is the clock with α = 1 at the wall.

No number in the README changes. Derivation: `results/kappa-channel-2026-11-14/CHANNEL_NOTE.md`.
