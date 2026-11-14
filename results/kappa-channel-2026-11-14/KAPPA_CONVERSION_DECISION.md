# Gap #2: the density / kappa convention. Decision document

HANDOFF Nov 14 2026 (5), Task A. RESEARCH ONLY: no file in `engine/` or `vendor/` is changed. Production λ = 0.

Marks: **Derived** = follows from the two lanes' equations; **Read** = read off a run; **Open** = not settled here.

Files: `research/lane_conversion_v1.py` (the conversion as one function pair, the chart states, the two arms),
`tests/test_lane_conversion.py` (A5), `two_lane_comparison.txt` (A3 raw output).

## Decision (A4)

**(a) The exact conversion, as one function pair.** A1 is clean. `cosmos_to_local` and `local_to_cosmos` in
`research/lane_conversion_v1.py` are the single statement of it; both arms of the comparison and both tests go through
them.

For Jeremiah's review (production, not made): move the two functions into `engine/`, and make the model's shared
constants one copy (below, "Two copies"). Until then the functions refuse to convert if the copies differ.

## A1. The conversion, derived

Both lanes use c = 1, G = 1. κ = 8π.

**Local lane** (BSSN, as coded: `engine/cmc_gauge.py:191`–`:196`, the Hamiltonian diagnostic with 16πρ):

    R + (2/3) K² − A_ij A^ij = 16π ρ                           (Hamiltonian constraint)
    ∂t K = −∇²α + α [ K_ij K^ij + 4π (ρ + S) ],   K_ij K^ij = K²/3 + (3/2) Aa²

On a uniform flat slice with α = 1, Aa = 0, K = −3H, S = 3p:

    H² = (8π/3) ρ            (Friedmann)
    dH/dt = −4π (ρ + p)       (Raychaudhuri, after using Friedmann)

**COSMOS lane** (`engine/cosmos.py:136`–`:146`):

    H² = ρ_E / 3
    dH/dt = −(ρ_E + p_E)/2

**Derived:** the two pairs are the same equations if and only if

    ρ_E = κ ρ,   p_E = κ p,   t_cosmos = t_local,   H = −K/3.

The pair alone would allow a common rescaling of time (t → s t, ρ_E → ρ_E/s²). That freedom is closed by the scalar:
both lanes carry the same number V0 in the same potential V0(e^{−λφ} − D) with a dimensionless φ, and the local
lane's S and D masses are 1 in the same time unit. So s = 1: **time and length units are the same in both lanes.**

Species by species:

| COSMOS | local lane | factor | how it is stored locally | mark |
|---|---|---|---|---|
| H | −K/3 | 1 | K | Derived |
| a | 1/X (a = b = 1) | 1, ratios only | conformal factor X | Derived; the absolute value of a is a label |
| φ, π | φ, Π | 1 | the local φ is already √κ × the canonical field | Derived |
| π²/2 + V(φ) | (Π²/2 + (∇φ)²/2 + V)/κ in the projection | 1/κ on the source | `scalar_system.py:142` | Derived |
| ρ_dm, ρ_b | the dust states' rest density | 1 in the state, 1/κ in the projection | `matter_system.py:395` | Derived |
| ρ_r (p = ρ/3) | the radiation state's density | same | same | Derived |
| β | `BETA_DM` | **1** | see below | Derived |
| ρ_shear (p = ρ) | (3/4) Aa² pointwise | — | cannot be uniform in spherical symmetry | Derived; untestable here |
| — | S, D | κ × their projected ρ, p | no COSMOS variable | Derived |

**β picks up no factor.** Write the canonical field χ with ρ = χ̇²/2 + U(χ) and coupling b: χ̈ + 3Hχ̇ + U′ = b ρ_phys,
ρ̇_phys + 3Hρ_phys = −b ρ_phys χ̇. Multiply through by κ with φ = √κ χ, V = κU, ρ_E = κ ρ_phys:

    φ̈ + 3Hφ̇ + V′(φ) = (b/√κ) ρ_E,     ρ̇_E + 3Hρ_E = −(b/√κ) ρ_E φ̇.

So β = b/√κ in *both* lanes, because both lanes evolve φ = √κ χ and both feed the Einstein-normalised ρ_dm into the
exchange (`production_kernel.py:288`–`:296` passes the raw rest density, not the projected one). A factor would
appear only if one lane switched to χ or to the projected density. In canonical terms b = √κ β = −0.2005. (β's value
stays borrowed; that is a separate label.)

**Shear.** From the Hamiltonian constraint with Aa ≠ 0, H² = (8π/3)ρ + Aa²/4, so ρ_shear = (3/4)Aa², and it falls as
a⁻⁶ in a homogeneous anisotropic space. About a regular centre Aa(0) = 0, so the spherical local lane cannot hold a
uniform shear. The conversion refuses a state with ρ_shear ≠ 0. This is a limit of spherical symmetry, not of the
convention.

**λ does not enter.** λ appears only in V_SD(S, D), and S and D have no COSMOS variable. With λ > 0 under the
one-budget rule D at rest in its vacuum adds the constant κ V_SD(0, D*) = −4π/λ to ρ_E; that is the rule already used
and tested in the lattice T1 run. Task A was therefore not repeated at λ = 5000.

**Two copies (Read).** The constants the conversion relies on exist twice: `CosmosParams` (V0, lam, D, beta) and
`scalar_system.V0 / LAMBDA / D_POT` with `production_kernel.BETA_DM`. They are equal today. Nothing in the code keeps
them equal. `check_shared_parameters` refuses to convert if they differ.

## A2. A homogeneous state from the model's own constants

No φ_today, no Ω, no present-day H. Chart labels (`chart_state`):

- **chart "V0":** φ = 0 at rest; ρ_dm = ρ_b = ρ_r = V0; H the expanding root of H² = ρ_E/3 (H = 1.035e-2).
- **chart "node"** (second chart, as a ruler): φ at the potential's zero, φ = ln(1/D)/λ; π = √V0; ρ_dm = 2V0,
  ρ_b = V0/2, ρ_r = V0/4 (H = 9.45e-3). This one has φ moving, so the β exchange is live from the first step.

These are chart labels, not claims about the universe.

## A3. The same state in both lanes

- COSMOS arm: `cosmos_rhs`, integrated to 1e-13.
- Local arm 1: the production kernel on a uniform flat slice (production edge), read at the centre cell.
- Local arm 2: the same state behind the research wall of GEAR-107, closure (i). Labelled: research, not production.

**Frozen state (no stepping), both charts, both local arms:** source density 2e-16; pressure 4e-17; Friedmann 1.5e-16;
K̇ against −3 dH/dt 2e-12; dΠ/dt against dπ/dt 5e-11 or better; round trip exact.

**Few steps, to t = 4** (largest |local − COSMOS|/value over a, H, φ, π, ρ_dm, ρ_b, ρ_r):

| chart | dt = 0.1 | 0.05 | 0.025 | ratio per halving |
|---|---|---|---|---|
| V0 (all but φ) | 5.6e-7 | 1.4e-7 | 3.5e-8 | 4.00 |
| V0, φ (starts from rest at zero) | 2.3e-5 | 5.7e-6 | 1.4e-6 | 4.00 |
| node | 2.4e-7 | 6.0e-8 | 1.5e-8 | 4.00 |

- **Read:** agreement. The residual is the local lane's time-stepping error: second order in dt in every variable,
  the same to all printed digits on the production edge and behind the wall.
- **Read:** it does not depend on the grid: n = 8, 16, 32 and r_max = 8, 16 give the same residuals to four digits
  (the state has no gradients).
- **Read:** the fields stay uniform to 1e-13 over these steps. The production edge leaves |Aa|/|K| of 1e-11 at the
  centre and up to 6e-10 at the last cell by t = 4. That is the outer-boundary gap (#3), small here, known to grow
  (1.3% by t = 40, PR #477).
- Longer evidence, labelled as such: the lattice T1 uniform cell (research wall, one-budget, λ = 5000) follows the same
  ODE to 7e-6 of H0 over t = 70, through the turnaround, second order in dt.

## A5. Regression tests (`tests/test_lane_conversion.py`)

- COSMOS lane: each chart obeys H² = ρ_E/3 and dH = −(ρ_E + p_E)/2, and the conversion hands the local lane K = −3H,
  the same φ, π and fluid densities, and sources ρ_E/κ, p_E/κ.
- Local lane: the production stress assembly, lapse solve and right-hand side reproduce those values on the frozen
  state; a few steps agree with the COSMOS lane with the residual falling fourfold per halving of dt.
- Negative control: β × √κ (the canonical coupling) is refused by the parameter check and misses dΠ/dt by more
  than 1%.
- Refusals: ρ_shear ≠ 0; S or D present.

## A6. 0* quantification

| item | verdict | why |
|---|---|---|
| the factor κ = 8π between the lanes' densities | **COORDINATE** | a normalisation; it changes no equation and no observable; removed exactly by one function |
| time and length units | **COORDINATE**, and equal | fixed by the shared V0 and the unit S, D mass |
| the absolute value of a | **COORDINATE** | only ratios of a have a local counterpart |
| β across the lanes | no factor (**Derived**) | both lanes use φ = √κ χ and ρ_E |
| the comparison's residual | **RULER-DEPENDENT in dt only**, → 0 as dt² | independent of grid, cell size, chart, and edge against wall |
| the chart states | chart labels | chosen, not derived; two charts agree |
| ρ_shear | **Open (structural)** | no uniform counterpart in spherical symmetry |
| two copies of the constants | **Open (production)** | equal today by hand; proposal above |

**Verdict:** the density / kappa convention is a COORDINATE matter and is tied down. Numbers from the two lanes can be
compared through `cosmos_to_local` / `local_to_cosmos`. Nothing borrowed was needed. What was already written in
`docs/UNITS_CONVENTION.md` is confirmed and extended (φ, β, radiation, time units, evolution; that note covered the
fluid density on a frozen state only).

Unpaid inputs: V0, D = 0.10 (the hand floor), λ_φ = 1, β = −0.04 (borrowed); the chart states. Forced: spherical
symmetry; α = 1 uniform slices; for arm 2, the wall's identical neighbours.
