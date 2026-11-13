"""RESEARCH ONLY. Outer boundaries and the return-flow test (HANDOFF to projects Nov 11 (2), Tasks 2 and 3).

Same physics as research/sd_bag_controls_v1.py (PR #470): flat space, no gravity, GEAR-67 quartic at lambda = 5000,
S and D both evolved (full back-reaction), cold start (D first relaxed with friction against the time-averaged S, so
the energy released where D falls into its well is removed before the test), smooth edge (width 1), seed 1e-10 inside.

    S'' = lap S - S - (lambda/2) S^3 - (3 lambda/2) S D^2
    D'' = lap D + D - (lambda/2) D^3 - (3 lambda/2) S^2 D

Geometries and boundaries:
  sph  absorb   spherical, sponge for r > 40 (outgoing energy is removed: PR #470's open domain)
  sph  reflect  spherical, no sponge, Neumann wall at r_max = 70 (everything that leaves comes back)
  row  periodic planar slab pocket in a periodic box of length L (a row of identical pockets, spacing L)
  row  absorb   the same slab with sponges at both ends (control for the row)
Energy bookkeeping (vacuum zero at D's ejected state, so the energy density is non-negative far from matter):
  e = S_t^2/2 + S'^2/2 + D_t^2/2 + D'^2/2 + V(S, D) + 1/(2 lambda)
  E_S(pocket)  = int over r < R + 3 of S_t^2/2 + S'^2/2 + S^2/2 + (lambda/8) S^4
  D's gradient + kinetic energy split by wavenumber band (FFT of D_t and D' on the grid): k < 1, 1-3, 3-10, > 10, and the
  top quarter of the grid's wavenumbers (grid scale).
Task 2 reports energy leaving (absorbed by the sponge) and returning (reflecting wall).
Task 3: does S's energy in the pocket settle (balance) or keep draining, and does D's energy pile up at the grid scale?

Unpaid input: lambda = 5000. No V0, phi, beta_D, gravity (labelled "no phi / no vacuum"). G_reset not involved.
Usage: python research/sd_return_flow_v1.py [T] [dr]       (all cases, x_in = 2.0 and 3.266)
"""
import math
import sys

import numpy as np

LAM = 5000.0
DSTAR = math.sqrt(2.0 / LAM)
VAC = 1.0 / (2 * LAM)


class Setup:
    def __init__(self, geom, bc, x_in, R=8.0, w_edge=1.0, dr=0.04, rmax=70.0, L=40.0, r_s=40.0, g0=2.0):
        self.geom, self.bc, self.x_in, self.R, self.dr = geom, bc, x_in, R, dr
        if geom == "sph":
            self.r = np.arange(1, int(round(rmax / dr)) + 1) * dr
            dist = self.r
            w = 4 * math.pi * self.r**2
        else:   # planar slab centred at L/2 in a box of length L (periodic) or 2L (absorbing ends)
            Lbox = L if bc == "periodic" else 3 * L
            n = int(round(Lbox / dr))
            self.r = (np.arange(n) + 0.5) * dr
            dist = np.abs(self.r - Lbox / 2)
            w = np.ones_like(self.r)
            self.Lbox = Lbox
        self.dist, self.w = dist, w
        x = x_in * 0.5 * (1 - np.tanh((dist - R) / w_edge))
        self.A = np.sqrt(x / LAM)
        if bc == "absorb":
            far = dist.max()
            start = r_s if geom == "sph" else 1.5 * L - 10.0
            self.gam = np.where(dist > start, g0 * ((dist - start) / (far - start)) ** 2, 0.0)
        else:
            self.gam = np.zeros_like(self.r)
        self.in_pocket = dist < R + 3.0

    # fields are stored as u = r f (sph) or f (row)
    def to_field(self, u):
        return u / self.r if self.geom == "sph" else u

    def lap(self, u):
        dr = self.dr
        if self.geom == "row" and self.bc == "periodic":
            return (np.roll(u, -1) - 2 * u + np.roll(u, 1)) / dr**2
        o = np.empty_like(u)
        o[1:-1] = (u[2:] - 2 * u[1:-1] + u[:-2]) / dr**2
        if self.geom == "sph":
            o[0] = (u[1] - 2 * u[0]) / dr**2
            g = u[-2] + 2 * dr * u[-1] / self.r[-1]          # Neumann on the field at r_max
        else:
            o[0] = (u[1] - u[0]) / dr**2                      # Neumann at both ends of the slab box
            g = u[-2]
        o[-1] = (g - 2 * u[-1] + u[-2]) / dr**2 if self.geom == "sph" else (u[-2] - u[-1]) / dr**2
        return o

    def grad(self, f):
        if self.geom == "row" and self.bc == "periodic":
            return (np.roll(f, -1) - np.roll(f, 1)) / (2 * self.dr)
        return np.gradient(f, self.dr)


def cold(st, T=200.0, gamma=0.3, seed=1e-10):
    r, dr = st.r, st.dr
    x = LAM * st.A**2
    D = np.sqrt(np.clip((2 / LAM) * (1 - 0.75 * x), 0, None))
    mult = st.r if st.geom == "sph" else 1.0
    w = mult * D; wp = w.copy()
    S2 = 0.5 * st.A**2
    dt = 0.4 * dr
    for _ in range(int(round(T / dt))):
        Dv = st.to_field(w)
        aw = st.lap(w) + w * (1 - 0.5 * LAM * Dv * Dv - 1.5 * LAM * S2)
        c = 0.5 * gamma * dt
        wn = (2 * w - (1 - c) * wp + dt * dt * aw) / (1 + c)
        wp, w = w, wn
    D = st.to_field(w) + seed * np.exp(-(st.dist / (0.5 * st.R)) ** 2)
    return D


def bands(st, Dt, Dr):
    """D's kinetic + gradient energy by wavenumber band (FFT; spherical: of r D_t and r D', a proxy)."""
    m = st.r if st.geom == "sph" else 1.0
    n = len(st.r); k = 2 * math.pi * np.fft.rfftfreq(n, st.dr)
    P = (np.abs(np.fft.rfft(m * Dt)) ** 2 + np.abs(np.fft.rfft(m * Dr)) ** 2)
    tot = P.sum() + 1e-300
    kq = 0.75 * k[-1]
    return [P[(k < 1)].sum() / tot, P[(k >= 1) & (k < 3)].sum() / tot, P[(k >= 3) & (k < 10)].sum() / tot,
            P[k >= 10].sum() / tot, P[k >= kq].sum() / tot]


def run(geom, bc, x_in, T=200.0, dr=0.04, every=1.0):
    st = Setup(geom, bc, x_in, dr=dr)
    D = cold(st)
    mult = st.r if geom == "sph" else 1.0
    u, w = mult * st.A, mult * D
    up, wp = u.copy(), w.copy()
    dt = 0.4 * dr

    def acc(u, w):
        S, Dv = st.to_field(u), st.to_field(w)
        return (st.lap(u) - u * (1 + 0.5 * LAM * S * S + 1.5 * LAM * Dv * Dv),
                st.lap(w) + w * (1 - 0.5 * LAM * Dv * Dv - 1.5 * LAM * S * S))
    au, aw = acc(u, w)
    u, w = u + 0.5 * dt * dt * au, w + 0.5 * dt * dt * aw
    vol = st.w * dr
    absorbed = 0.0
    rec = []
    rec_every = int(round(every / dt))
    E0 = None
    for n in range(1, int(round(T / dt)) + 1):
        au, aw = acc(u, w)
        c = 0.5 * st.gam * dt
        un = (2 * u - (1 - c) * up + dt * dt * au) / (1 + c)
        wn = (2 * w - (1 - c) * wp + dt * dt * aw) / (1 + c)
        St, Dt = st.to_field(un - up) / (2 * dt), st.to_field(wn - wp) / (2 * dt)
        absorbed += float(np.sum(st.gam * (St**2 + Dt**2) * vol)) * dt
        if n % rec_every == 0 or n == 1:
            S, Dv = st.to_field(u), st.to_field(w)
            Sr, Dr = st.grad(S), st.grad(Dv)
            e = (0.5 * St**2 + 0.5 * Sr**2 + 0.5 * Dt**2 + 0.5 * Dr**2 + 0.5 * S**2 - 0.5 * Dv**2
                 + LAM / 8 * (S**4 + 6 * S**2 * Dv**2 + Dv**4) + VAC)
            eS = 0.5 * St**2 + 0.5 * Sr**2 + 0.5 * S**2 + LAM / 8 * S**4
            Etot = float(np.sum(e * vol)); ES = float(np.sum((eS * vol)[st.in_pocket]))
            if E0 is None:
                E0 = Etot
            rec.append([n * dt, Etot, absorbed, ES, float(np.abs(Dv[st.in_pocket]).max()) / DSTAR] + bands(st, Dt, Dr))
        up, u, wp, w = u, un, w, wn
    return np.array(rec), E0


def report(T=200.0, dr=0.04):
    print(f"lambda = {LAM:g}; flat space, no gravity, no phi / no vacuum; cold start; R = 8, smooth edge; dr = {dr}; T = {T}")
    print("E_S = S's energy in the pocket (r < R + 3) relative to its start; D bands = D kinetic+gradient energy fractions")
    for x_in in (2.0, 3.266):
        for geom, bc in (("sph", "absorb"), ("sph", "reflect"), ("row", "absorb"), ("row", "periodic")):
            rec, E0 = run(geom, bc, x_in, T, dr)
            t, Et, Ab, ES = rec[:, 0], rec[:, 1], rec[:, 2], rec[:, 3]
            print(f"\nx_in = {x_in}  {geom:3s} {bc:8s}: total energy start {E0:.4e}; at end in box {Et[-1]:.4e}, absorbed {Ab[-1]:.4e}, "
                  f"box + absorbed - start = {(Et[-1] + Ab[-1] - E0)/E0:+.2e} of start")
            print("     t     E_S/E_S(0)   max|D|/D* in pocket   D bands: k<1   1-3    3-10   >10    grid-scale(top quarter)")
            for tt in (0, 10, 20, 40, 60, 80, 100, 130, 160, 200):
                if tt > t[-1] + 1e-9:
                    continue
                i = int(np.argmin(np.abs(t - tt)))
                b = rec[i, 5:]
                print(f"  {t[i]:6.1f}   {ES[i]/ES[0]:8.3f}       {rec[i,4]:8.3f}          {b[0]:.3f}  {b[1]:.3f}  {b[2]:.3f}  {b[3]:.3f}  {b[4]:.2e}")


if __name__ == "__main__":
    T = float(sys.argv[1]) if len(sys.argv) > 1 else 200.0
    dr = float(sys.argv[2]) if len(sys.argv) > 2 else 0.04
    report(T, dr)
