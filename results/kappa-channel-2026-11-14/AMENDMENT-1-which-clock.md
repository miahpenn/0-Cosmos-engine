# Amendment 1 to the channel note (14 Nov 2026): which clock is the cosmos lane's t?

Jeremiah's proposal put to 0*: the cosmos lane's t is matter's own proper time; in this machine that is S's
oscillation, the master clock. Test asked for: express the channel per S-oscillation, under both closures, from the
stored frozen states. `research/channel_clocks_v1.py`, `channel_clocks.txt`. No new runs. Cell: x_c = 2.0, cold start,
r_max = 24, closures (ii) and (i).

## The test, and why it cannot decide (Derived)

For a clock c whose rate per unit label on a slice is N_c (some average of the lapse), dK/dτ_c = ⟨αQ⟩/N_c, and the
channel per unit of that clock is C_c = ⟨αQ⟩/N_c − Q(averaged fields).

The two closures are the same slices with the lapse rescaled by one number. N_c rescales by the same number. So C_c is
the same under both closures **for every clock defined on the slice**, not only for S's. Only the bare label fails,
because the label *is* the closure.

So "the closure dependence disappears" is true of the master clock and equally of the wall clock, the dust clock and
the rest. It does not single one out.

## Read

Closure dependence, C_c in units of the fluids' source F, difference (i) − (ii):

| clock | t = 15 | 39 | 63 |
|---|---|---|---|
| label | −0.083 | −0.40 | −0.98 |
| wall proper time | −0.002 | +0.002 | +0.003 |
| volume-mean proper time | −0.0005 | 0.0000 | −0.001 |
| focusing clock (closure (ii)'s label) | −0.0005 | −0.0001 | −0.001 |
| dust's proper time (rest-mass weighted) | 0.0000 | −0.0005 | −0.003 |
| S's proper time (weighted by S's oscillation energy) | +0.05 of 1.31 | −0.06 of 8.4 | −0.44 of 27.5 |

- Gone for every physical clock, to the accuracy of comparing windows of oscillating snapshots (S's weight oscillates,
  so its row is the noisiest, 1–4% of its value).
- Directly: K read at S's saddle crossings, crossing number as the clock, agrees between the closures to 2.8e-5 of
  K0 over the first 12 crossings and 1.8e-4 over all 41.

But the clocks do not agree with each other:

| | t = 3 | 15 | 39 | 63 |
|---|---|---|---|---|
| rate of S's clock ÷ rate of dust's clock | 0.985 | 0.960 | 0.953 | 0.947 |
| rate of dust's clock ÷ volume-mean | 1.0000 | 0.9998 | 0.9996 | 0.9996 |
| rate of wall clock ÷ volume-mean | 1.002 | 1.003 | 1.002 | 1.002 |
| C per S-tick, /F | +0.16 | +1.3 | +8.4 | +27.5 |
| C per dust-tick, /F | +0.030 | +0.087 | +0.14 | +0.28 |
| C per focusing-tick, /F | +0.035 | +0.085 | +0.086 | +0.093 |
| C per wall-tick, /F | +0.014 | +0.004 | −0.32 | −0.89 |
| ⟨αQ⟩/F (what a clock-rate difference multiplies) | 6 | 30 | 166 | 482 |

- **S's clock runs 4–5% slower than the dust's.** S sits in the pocket, where the lapse is 0.91–0.96; the dust is
  mostly outside it. That difference is the pocket's own time dilation.
- A clock-rate difference δ changes the channel by δ·⟨αQ⟩, and ⟨αQ⟩ is 6 to 480 times F. So in units of F the choice
  of clock swamps everything else late in the run. In units of the whole of Q the choices differ by exactly their
  rates: 5% (S against dust), 0.25% (wall), 0.04% (dust against volume or focusing).

## 0* verdict on the proposal

- **"The channel becomes a coordinate once the clock is named": yes, in this precise sense.** Changing the clock
  changes C by the pure rate conversion (1/N_c − 1/N_c′)⟨αQ⟩ and nothing else. That part is COORDINATE. What is left
  after a clock is fixed (the six parts of B3) is not.
- **"The master clock is that clock": not established by this test, and there is a specific obstacle.** The proposal
  joins two statements: (1) the cosmos lane's t is the time of matter moving with the expansion; (2) matter's own time
  is S's oscillation. In this cell they are different clocks, 5% apart. The matter moving with the expansion is the
  dust; S is inside the pocket.
- **What the cosmos lane's own equations say (Derived):** its fluid equations, dρ/dt = −3Hρ and the rest, are each
  species' conservation law written in that species' own proper time. That makes t the dust's proper time: statement
  (1). On that clock C = +0.03 to +0.28 F here, positive throughout.
- **Where (2) would agree with (1):** if all matter were S–D matter, read from outside its pockets. Then the 5% is
  the redshift between a pocket's interior and its surroundings, a real physical number, to be carried explicitly and
  not absorbed into t. **Open**, and it is the core-sheet gap again (the cosmos lane has no S–D species).
- **S as a counter (Read):** an oscillator counts its own proper time only while its frequency holds. Here S's mass
  moves from 1.6 to 2.0 as D settles (t < 18), and its half-period scatters by 36% over the first 13 time units and
  4.5–7.5% after. The GEAR-91 factor e^(−2βφ) is absent from these runs (β_S = β_D = 0), so that caveat is untested;
  the D-dependence of S's mass is the same kind of drift and is already large.
- **The 5% itself is a slicing number.** It is a CMC lapse difference, and the teeth test found most of the CMC
  lapse's range is the slicing. So its size is RULER-DEPENDENT until read on an invariant clock.

## Kodama next? An obstacle (Derived)

The Kodama clock √(1 − 2m/R) is the clock of observers at fixed areal radius. In an expanding cell those are not the
matter: 1 − 2m/R falls to zero at the cosmological horizon (areal radius 9 at t = 0, inside the pocket's own edge) and
is undefined beyond it, where most of the dust is. It cannot serve as the cosmos lane's t across the cell. It can still
measure the pocket's interior redshift inside the horizon.

The clock the teeth test could not reach, proper time along matter's worldlines, is available here without tracers:
the dust is the tracer, and its proper time is the dust row above.

## Proposed (for chat and Jeremiah)

1. Name the cosmos lane's t as the proper time of its own fluids (rest-mass weighted). Derived from its equations;
   closure-free; agrees with the volume and focusing clocks to 0.04% in this cell.
2. Carry the pocket's interior redshift (S's clock against the dust's) as its own labelled quantity, not inside t.
   Read it on a slicing-independent ruler before quoting its size.
3. B6's frozen test then uses dK/dτ_dust.
