"""Verifier v_reach: fragmentation growth rates of the homogeneous phi != 0 rotating condensate, two independent routes
(neither uses the charge_phase cubic):
 (A) 6x6 first-order linear system in the co-rotating frame (re-derived by hand from the field equations),
     max Re(eigenvalue) over k;
 (B) direct 1D nonlinear lattice evolution (spectral derivatives, leapfrog) of p, q (X = v(p+iq)/sqrt2) and S,
     background p+iq = sqrt(x) e^{i Omega t}, S = S0, seeded at the fastest k; growth rate fitted.
 (C) small-x closed form: NR attractive condensate, gamma_max = (a^2-b) x/2, k_max = sqrt((a^2-b) x).
Units: sqrt(lam_phi) v for rates, ell^-1 for k. a = 1, b = 0.1.
"""
import json, sys
import numpy as np
from scipy.optimize import minimize_scalar

a, b = 1.0, 0.1


def bg(x):
    S0 = np.sqrt(1 - a * x); Om = np.sqrt(a * S0**2 + b * x)
    return S0, Om


def mat(x, k):
    S0, Om = bg(x)
    ms2 = 3 * S0**2 - 1 + a * x
    c = 2 * a * S0 * np.sqrt(x)
    M = np.zeros((6, 6))
    M[0, 3] = M[1, 4] = M[2, 5] = 1
    # u'' = 2 Om w' - (k^2 + 2 b x) u - c s
    M[3, 4] = 2 * Om; M[3, 0] = -(k**2 + 2 * b * x); M[3, 2] = -c
    # w'' = -2 Om u' - k^2 w
    M[4, 3] = -2 * Om; M[4, 1] = -k**2
    # s'' = -(k^2 + ms2) s - c u
    M[5, 2] = -(k**2 + ms2); M[5, 0] = -c
    return M


def growth(x, k):
    return float(np.max(np.linalg.eigvals(mat(x, k)).real))


def gmax(x):
    ks = np.sqrt(x) * np.linspace(0.05, 2.5, 400)
    g = np.array([growth(x, k) for k in ks])
    i = int(np.argmax(g))
    r = minimize_scalar(lambda k: -growth(x, k), bracket=(ks[max(i - 1, 0)], ks[i], ks[min(i + 1, len(ks) - 1)]),
                        tol=1e-12)
    return -r.fun, r.x


def lattice_growth(x, k, nper=1, N=32, dt=None, tmax=None):
    """1D periodic box of one wavelength 2 pi/k; returns fitted growth rate of the k-mode of |psi|^2."""
    S0, Om = bg(x)
    L = 2 * np.pi / k * nper
    z = np.arange(N) * L / N
    kk = 2 * np.pi * np.fft.fftfreq(N, d=L / N)
    lap = lambda f: np.fft.ifft(-(kk**2) * np.fft.fft(f)).real
    gam_est = (a * a - b) * x / 2
    if dt is None:
        dt = min(0.02, 0.2 / max(kk.max(), 1.0))
    if tmax is None:
        tmax = 22.0 / gam_est
    eps = 1e-9
    p = np.sqrt(x) * (1 + eps * np.cos(k * z)) + 0 * z
    q = np.zeros(N)
    S = S0 * np.ones(N)
    pd = np.zeros(N); qd = np.sqrt(x) * Om * (1 + eps * np.cos(k * z)); Sd = np.zeros(N)

    def acc(p, q, S):
        F2 = p * p + q * q
        m2 = a * S * S + b * F2
        return lap(p) - m2 * p, lap(q) - m2 * q, lap(S) - S * (S * S - 1) - a * S * F2

    ap, aq, aS = acc(p, q, S)
    nsteps = int(tmax / dt)
    ts, amps = [], []
    for n in range(nsteps):
        pd += 0.5 * dt * ap; qd += 0.5 * dt * aq; Sd += 0.5 * dt * aS
        p += dt * pd; q += dt * qd; S += dt * Sd
        ap, aq, aS = acc(p, q, S)
        pd += 0.5 * dt * ap; qd += 0.5 * dt * aq; Sd += 0.5 * dt * aS
        if n % 20 == 0:
            F2 = p * p + q * q
            amp = np.abs(np.fft.fft(F2)[nper]) / N
            ts.append(n * dt); amps.append(amp)
    ts = np.array(ts); amps = np.array(amps)
    la = np.log(amps)
    # fit in the window where the mode amplitude is between 1e-7 x and 1e-3 x (linear, past transients)
    sel = (amps > 1e-7 * x) & (amps < 1e-3 * x)
    if sel.sum() < 10:
        return None
    A = np.vstack([ts[sel], np.ones(sel.sum())]).T
    slope = np.linalg.lstsq(A, la[sel], rcond=None)[0][0]
    return float(slope)


if __name__ == '__main__':
    out = dict(rows=[])
    for x in [1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 0.1, 0.25, 0.5]:
        g, k = gmax(x)
        row = dict(x=x, gamma_max_6x6=g, k_max_6x6=k, gamma_over_x=g / x, k_over_sqrtx=k / np.sqrt(x),
                   closed_form_gamma=(a * a - b) * x / 2, closed_form_k=np.sqrt((a * a - b) * x))
        if x >= 0.01:
            row['gamma_lattice_1D'] = lattice_growth(x, k)
        out['rows'].append(row)
        print({kk: (round(v, 6) if isinstance(v, float) else v) for kk, v in row.items()}, flush=True)
    json.dump(out, open(sys.argv[1], 'w'), indent=1)
