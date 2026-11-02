#!/usr/bin/env python3
"""inf_potential (Claude, own code): calibrate each potential family to A_s, H at the slice and (where the family has a
free shape parameter) n_s, with N_* consistent with the stated post-inflation history.

Slice conventions (all reported; a = 1 at the slice when handed to the engine):
  'S1'   : 1 e-fold before eps_H = 1 (the end of the slow-roll hold; |eta_V| ~ 1 there for exponential plateaus).
           This is the PRIMARY convention: H_i = H(S1).
  'Send' : eps_H = 1.
History conventions for N_*:
  ('gamma', g) : homogeneous oscillation with constant width Gamma_I = g H_i (the engine's mode-0 convention; g = 0.5 is
                 M1r's supplied value), integrated until radiation holds all but 1e-7 of the energy; then entropy
                 conservation with g_* = 106.75.
  ('wT', w, T) : constant w from the end of inflation to T_RH = T (GeV), then radiation.
Labels: [C] given the model choices.
"""
import math, sys
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
sys.path.insert(0, str(Path(__file__).resolve().parent))
import inf_common as ic


def slice_N(conv):
    return {'S1': 1.0, 'Send': 0.0, 'S2': 2.0, 'S05': 0.5}[conv]


def a_end_over_a0_gamma(m, bg, V0, gam_over_Hi, h_slice, tmax_HiUnits=400.0):
    """a_end/a_0 for a constant perturbative width (homogeneous); also returns K = a^2 H (H_i units, a = 1 at the end) and
    the instantaneous-equality time."""
    g = gam_over_Hi * h_slice               # width in sqrt(V0) units
    tend = max(40.0 / g, 20.0) if g > 0 else 50.0
    sol = ic.post_inflation(m, bg, g, tend, rtol=1e-9, atol=1e-13)
    lna, f, fd, rr = sol.y[:, -1]
    rho_tot = 0.5 * fd * fd + m.v(f) + rr
    if rr / rho_tot < 1 - 1e-6:
        raise RuntimeError('not radiation dominated at t_end: %g' % (rr / rho_tot))
    rho_GeV = rr * V0 * ic.MPL_GEV ** 4
    Tx = (30 * rho_GeV / (math.pi ** 2 * ic.G_RH)) ** 0.25
    a0x = (ic.GS0 / ic.G_RH) ** (1 / 3) * ic.T0_GEV / Tx
    a_end_a0 = a0x / math.exp(lna)
    H = math.sqrt(rho_tot / 3)
    K_end = math.exp(2 * lna) * H / h_slice     # a^2 H in H_i units with a = 1 at the END of inflation
    return a_end_a0, dict(K_end_norm=K_end, T_last=Tx, sol=sol)


def nstar(m, bg, V0, hist, h_slice, cache=None):
    """N_* for k = 0.05/Mpc: N_* = ln(H_*/k_*) + ln(a_end/a_0), solved by fixed point."""
    f, fN = ic.bg_at(bg, 0.0)
    h_end, _ = ic.h_of(m, f, fN)
    rho_end = 3 * h_end ** 2 * V0 * ic.MPL_GEV ** 4
    if hist[0] == 'gamma':
        if cache is not None and 'aea0' in cache:
            aea0 = cache['aea0']
        else:
            aea0, info = a_end_over_a0_gamma(m, bg, V0, hist[1], h_slice)
            if cache is not None:
                cache['aea0'] = aea0; cache['info'] = info
    elif hist[0] == 'wT':
        aea0, Nre, Trh = ic.a_end_over_a0_wT(rho_end, hist[2], hist[1])
    elif hist[0] == 'inst':
        aea0, Nre, Trh = ic.a_end_over_a0_wT(rho_end, 1e30, 1 / 3)
    N = 53.0
    for it in range(50):
        f, fN = ic.bg_at(bg, N)
        h, _ = ic.h_of(m, f, fN)
        Hs = math.sqrt(V0) * h * ic.MPL_GEV
        Nn = ic.Nstar_from_history(Hs, aea0)
        if abs(Nn - N) < 1e-7:
            break
        N = Nn
    return N, dict(rho_end=rho_end, T_inst=ic.T_instant(rho_end), a_end_over_a0=aea0)


def evaluate(m, Hi_GeV, hist, conv='S1', bg=None, need_spec=True):
    """Given a shape, set V0 so that H(slice) = Hi, find N_*, compute the pivot spectrum."""
    if bg is None:
        bg = ic.background(m, 75.0)
    Ns = slice_N(conv)
    f, fN = ic.bg_at(bg, Ns)
    hs, _ = ic.h_of(m, f, fN)
    V0 = (Hi_GeV / ic.MPL_GEV / hs) ** 2
    cache = {}
    Nst, ninfo = nstar(m, bg, V0, hist, hs, cache)
    out = dict(V0=V0, h_slice=hs, Nstar=Nst, bg=bg, ninfo=ninfo, cache=cache)
    if need_spec:
        sp = ic.spectrum_at_pivot(m, bg, Nst)
        out.update(As=V0 * sp['Pz'], ns=sp['ns'], r=sp['r'], alpha_s=sp['alpha_s'], n_t=sp['n_t'])
        f, fN = ic.bg_at(bg, Nst)
        h, _ = ic.h_of(m, f, fN)
        out['Hstar_GeV'] = math.sqrt(V0) * h * ic.MPL_GEV
    return out


def calibrate_alpha(kind, Hi_GeV, hist, conv='S1', a_guess=None):
    """T/E: find alpha with A_s = target at H(slice) = Hi."""
    if a_guess is None:
        a_guess = 6e-5 * (Hi_GeV / 1.29e11) ** 2
    def F(la):
        m = ic.Model(kind, alpha=math.exp(la))
        r = evaluate(m, Hi_GeV, hist, conv)
        return math.log(r['As'] / ic.AS_TARGET), r, m
    la0 = math.log(a_guess)
    # secant in ln alpha (A_s ~ 1/alpha at fixed H)
    x0, x1 = la0, la0 + 0.05
    f0 = F(x0)[0]
    for it in range(30):
        f1, r1, m1 = F(x1)
        if abs(f1) < 2e-6:
            return m1, r1
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        x0, f0, x1 = x1, f1, x2
    raise RuntimeError('alpha calibration failed')


def calibrate_hilltop(Hi_GeV, hist, p, beta=0.0, conv='S1', vev_guess=0.4):
    """HT at fixed (p, beta): vev with A_s = target."""
    def F(lv):
        m = ic.Model('HT', vev=math.exp(lv), p=p, beta=beta)
        r = evaluate(m, Hi_GeV, hist, conv)
        return math.log(r['As'] / ic.AS_TARGET), r, m
    x0, x1 = math.log(vev_guess), math.log(vev_guess) + 0.02
    f0 = F(x0)[0]
    for it in range(40):
        f1, r1, m1 = F(x1)
        if abs(f1) < 2e-6:
            return m1, r1
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        x0, f0, x1 = x1, f1, x2
    raise RuntimeError('hilltop calibration failed')


def calibrate_hilltop_ns(Hi_GeV, hist, p, ns_target=ic.NS_TARGET, conv='S1', x0=(0.3785, 1e-5), tol=2e-6, verbose=False):
    """HT with the linear tilt: (vev, beta) with A_s = target and n_s = ns_target (2D Newton in ln vev, ln beta)."""
    def F(z):
        m = ic.Model('HT', vev=math.exp(z[0]), p=p, beta=math.exp(z[1]))
        r = evaluate(m, Hi_GeV, hist, conv)
        return np.array([math.log(r['As'] / ic.AS_TARGET), (r['ns'] - ns_target) * 100.0]), r, m
    z = np.array([math.log(x0[0]), math.log(x0[1])])
    for it in range(25):
        f, r, m = F(z)
        if verbose:
            print('  it', it, 'vev %.6f beta %.4e As %.5e ns %.5f' % (m.vev, m.beta, r['As'], r['ns']), flush=True)
        if np.max(np.abs(f)) < tol:
            return m, r
        J = np.zeros((2, 2)); h = 0.01
        for j in range(2):
            dz = np.zeros(2); dz[j] = h
            J[:, j] = (F(z + dz)[0] - f) / h
        step = -np.linalg.solve(J, f)
        lim = 0.3
        if np.max(np.abs(step)) > lim:
            step *= lim / np.max(np.abs(step))
        z = z + step
    raise RuntimeError('hilltop (A_s, n_s) calibration failed')


def calibrate_param(kind, pname, guess, Hi_GeV, hist, conv='S1', fixed=None, tol=2e-6):
    """generic 1D calibration of one shape parameter (secant in ln param) to A_s."""
    fixed = dict(fixed or {})
    def F(lx):
        kw = dict(fixed); kw[pname] = math.exp(lx)
        m = ic.Model(kind, **kw)
        r = evaluate(m, Hi_GeV, hist, conv)
        return math.log(r['As'] / ic.AS_TARGET), r, m
    x0, x1 = math.log(guess), math.log(guess) + 0.05
    f0 = F(x0)[0]
    for it in range(40):
        f1, r1, m1 = F(x1)
        if abs(f1) < tol:
            return m1, r1
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        x0, f0, x1 = x1, f1, x2
    raise RuntimeError('calibration failed')
