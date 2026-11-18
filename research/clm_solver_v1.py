import numpy as np
def make(N, nu):
    k = np.fft.fftfreq(N, 1.0 / N); deal = (np.abs(k) <= N / 3).astype(float)
    return k, deal
def rhs_nl(wh, k, deal, sgn):
    # H symbol: -i sgn(k) (sgn=+1) ; omega*H(omega) in physical space, dealiased
    Hh = -1j * sgn * np.sign(k) * wh
    w = np.fft.ifft(wh).real; H = np.fft.ifft(Hh).real
    return np.fft.fft(w * H) * deal
def run(N, nu, ic, sgn=1, T=35.0, dt=5e-4, gamma=0.0, band=(10, 20), fext=None, blow=1e3, rec=False, every=1):
    """fext(t) -> f_hat array or None. gamma: closed-loop feedback forcing gamma*nu*k^2*w_hat in band."""
    k, deal = make(N, nu); x = 2 * np.pi * np.arange(N) / N
    wh = np.fft.fft(ic(x)) * deal
    inb = ((np.abs(k) >= band[0]) & (np.abs(k) <= band[1])).astype(float)
    def F(wh, t):
        f = gamma * nu * k**2 * wh * inb if gamma else 0.0
        if fext is not None: f = f + fext(t)
        return rhs_nl(wh, k, deal, sgn) - nu * k**2 * wh + f
    n = int(round(T / dt)); out = []; frec = []; t = 0.0
    for i in range(n):
        if rec and i % every == 0: frec.append((gamma * nu * k**2 * wh * inb)[inb > 0].copy())
        k1 = F(wh, t); k2 = F(wh + 0.5 * dt * k1, t + 0.5 * dt); k3 = F(wh + 0.5 * dt * k2, t + 0.5 * dt); k4 = F(wh + dt * k3, t + dt)
        wh = wh + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4); t += dt
        if i % 200 == 0:
            m = np.abs(np.fft.ifft(wh).real).max(); out.append((t, m))
            if m > blow or not np.isfinite(m): return dict(blew=True, t=t, out=out, frec=frec, inb=inb)
    return dict(blew=False, t=t, out=out, frec=frec, inb=inb, wmax=np.abs(np.fft.ifft(wh).real).max())
