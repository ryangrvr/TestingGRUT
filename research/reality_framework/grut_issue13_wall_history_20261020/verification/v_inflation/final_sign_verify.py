"""p_- after post-inflation trapping (verifier v_inflation; independent of module code).

Forward 1D FP from a delta at y0 = v/H (fresh start) on the full line, by eigen-decomposition of the
symmetrized central-difference generator (reflecting ends).  Each end value phi_end = y H is mapped to its
trapped sign with the RD ODE of postinf_lame.crossings, started at the end of inflation
(eta0 = sqrt(sqrt(lam) v / H_inf), conformal units a = eta, H = 1/eta^2 in units sqrt(lam) v).
lam = 0.01 (tree; RD era, xi term absent since R = 0).
"""
import json, math, sys
import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import brentq
import fp1d_verify as F
from postinf_lame import crossings

DC = F.DC


def end_distribution(H, N, h=0.004):
    w, wp, a2, lam = F.wfun(H)
    y0 = F.V_GEV / H
    k = max(4, int(round(y0 / h))); h = y0 / k
    wmin = w(math.sqrt(max(-a2, 0) / lam)) if a2 < 0 else 0.0
    Y = max(1.05 * y0, 1.0)
    while (8 * math.pi ** 2 / 3) * (w(Y) - wmin) < 60:
        Y *= 1.05
    m = int(math.ceil(Y / h))
    y = h * np.arange(-m, m + 1)
    A = -wp(y) / 3
    up = DC / h ** 2 + A / (2 * h); lo = DC / h ** 2 - A / (2 * h)
    up[-1] = 0.0; lo[-1] = 2 * DC / h ** 2
    lo[0] = 0.0; up[0] = 2 * DC / h ** 2
    assert np.all(up[:-1] > 0) and np.all(lo[1:] > 0)
    lpi = np.zeros(len(y))
    for i in range(len(y) - 1):
        lpi[i + 1] = lpi[i] + math.log(up[i]) - math.log(lo[i + 1])
    lpi -= lpi.max()
    d = up + lo
    e = -np.sqrt(up[:-1] * lo[1:])
    lam_, V = eigh_tridiagonal(d, e)
    j0 = m + k
    p = np.exp(0.5 * (lpi - lpi[j0])) * (V @ (np.exp(-lam_ * N) * V[j0, :]))
    return y, p


_sign_cache = {}


def final_sign(Phi, eta0):
    key = (round(Phi, 6), round(eta0, 6))
    if key not in _sign_cache:
        if Phi < 1.0:
            _sign_cache[key] = 1
        else:
            nz, ta, sg = crossings(Phi, eta0=eta0, rtol=1e-9)
            _sign_cache[key] = sg
    return _sign_cache[key]


def p_after(H, N):
    y, p = end_distribution(H, N)
    eta0 = math.sqrt(0.1 * F.V_GEV / H)
    pend = p[y < 0].sum() + 0.5 * p[y == 0].sum()
    neg = 0.0
    for yi, pi in zip(y, p):
        if pi < 1e-14 or yi == 0.0:
            continue
        Phi = abs(yi) * H / F.V_GEV
        s = final_sign(Phi, eta0) * (1 if yi > 0 else -1)
        if s < 0:
            neg += pi
    neg += 0.5 * p[y == 0].sum()
    return dict(p_end=float(pend), p_after=float(neg), norm=float(p.sum()))


if __name__ == '__main__':
    out = {}
    for H in (1000.0, 1400.0, 2000.0, 3000.0):
        out['%g' % H] = {('N%d' % N): p_after(H, float(N)) for N in (50, 60)}
        print(H, out['%g' % H], flush=True)
    thr = {}
    for N in (50, 60):
        thr['N%d' % N] = brentq(lambda H: p_after(H, float(N))['p_after'] - 0.3116077, 1100.0, 1800.0, xtol=1.0)
        print(N, thr['N%d' % N], flush=True)
    out['H_perc_after_trapping'] = thr
    json.dump(out, open('results_final_sign.json', 'w'), indent=1)
