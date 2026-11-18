# AMENDMENT-2 (C6_NOTE): the share s was inverted (18 Nov 2026; docs only, no run)

Applies to `results/sd-equation-of-state-2026-11-15/C6_NOTE.md` (§3, variant B and its table). Found by chat's hierarchy pass and
confirmed by Duet (coordination log, 18 Nov). Not to be confused with the Drive documents also called "Amendment 2"; this file is cited by path.

## What was wrong (Derived)
With s = S–D matter's share of the total density, ρ_tot = ρ_SD / s, and H² = (κ/λ) ρ_tot / 3 gives

    H² = (κ/λ) ρ̃ / (3 s)        (C6_NOTE §3 wrote H² = (κ/λ) s ρ̃ / 3: s inverted).

Variant B's balance therefore fixes **κ*/λ* / s**, not κ*/λ* · s: κ*/λ* = (column value) × s.

## Table fix (values unchanged, labels changed)
- Column "B: κ*/λ*·s" (0.043–0.31 at x_A = 1.28) is 3H²/ρ̃ with H = the archive band 0.139–0.373. Read it as **κ*/λ* ÷ s**.
- Column "B: λ*/s" (80–590) is 8π / that value. Read it as **λ*·s**. So λ* = (80–590) / s.
- With production's grip κ/λ = 8π/1.57×10⁶ = 1.60×10⁻⁵: **s = 5×10⁻⁵ to 3.7×10⁻⁴** (x_A = 1.28); production's λ sits inside (80–590)/s for that s.

## The share s at the lane's x = 4/3 crossing (Derived arithmetic from stored #509 numbers; no run)
From `results/sd-dark-matter-2026-11-16/POKE2_RECONCILIATION.md`: at x = 4/3 the lane has |H|/m_S = 0.30, so |H| = 0.60 (m_S = 2).
ρ̃(4/3) = x/2 + x²/8 + ½ = 1.389. Then s = (κ/λ) ρ̃ / (3H²) = **2.1×10⁻⁵** (1.99×10⁻⁵ at x_A = 1.28).
- Chat's "f_S ≈ 3×10⁻⁵" is not stored in #509's files; this recomputation gives about 2×10⁻⁵ (same order). Label: arithmetic on stored numbers.
- s ∝ 1/H² at a fixed grip, so this is **not independent of the C6 requirement**: the lane's crossing has |H| = 0.60, which is
  ×1.6 above the clock-stop band's upper edge (0.373) and ×4.3 above its lower edge (0.139); that is exactly why s is ×2.5–×18
  below the 5×10⁻⁵–3.7×10⁻⁴ needed. The real comparison is |H| against the band, not s.

## Labels
Carrier = φ's kinetic energy on the contracting branch (Read, on a constructed branch: the lane's crunch follows φ's constructed
positive term, AUDIT-4; not the account's design target, `docs/STATUS.md`). Whether the same s holds for a local reset in an
expanding background is Open. A radiation component for the bar (CMB, light elements) is still needed; its seed (Ω_r = 9.2×10⁻⁵) is borrowed.
