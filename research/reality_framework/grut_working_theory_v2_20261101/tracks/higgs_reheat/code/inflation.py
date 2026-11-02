"""Item 2 of track higgs_reheat: the Higgs during inflation in M1r (H = H_i constant, R = 12 H^2, xi_H R h^2/2).

Own code.  Writes ../results/inflation.json.  Run: python3 -B inflation.py  (~3 min, one thread).

Quantities [C]:
  * xi(mu) along the 3-loop trajectory for xi(M_t) = 0.1 (and for xi(H_i) = 0.1 as the alternative convention).
  * Delta V_dS(h): de Sitter one-loop RG-improved potential (smpot.SM.V_dS, 1804.02020 Tables 1-2), its barrier h_b and
    height; the Hawking-Moss exponent B = 8 pi^2 Delta V_b/(3 H^4) (Hawking-Moss; equivalently the Fokker-Planck
    equilibrium P ~ exp(-8 pi^2 V/3H^4), Starobinsky-Yokoyama).
  * Expected number of AdS bubbles in our past light cone, <N> = (4 pi/9) e^{3N} e^{-B} (constant H; the form used by
    Mantziris-Markkanen-Rajantie 2011.03763 reduces to this for constant H [A]); survival needs <N> < 1.
  * xi_min: the smallest xi(M_t) (and xi(H_i)) with <N> = 1, N = 50, 55, 60.
  * Variances: stochastic light-field <phi_i^2> = 3H^4/(8 pi^2 m^2) per real component; the exact de Sitter
    super-horizon (k < eps a H) variance from Bunch-Davies Hankel modes, P(x) = (H^2/8 pi) x^3 |H_nu(x)|^2, over N
    e-folds (BD normalisation e^{i pi nu/2}, so |.|^2 = e^{-pi Im nu} for a heavy field); the curvature mass m^2 = 2 Delta V/h^2 at h -> 0 (one loop).
  * Per-patch probability that the 4-component Higgs exceeds the FLAT barrier h_max (2-loop) and the de Sitter barrier
    h_b: chi^2_4 tail e^{-x}(1+x), x = h_c^2/(2 sigma^2); and the Fokker-Planck (non-Gaussian) tail with the full
    potential and radial measure h^3.  Expected counts among e^{3N} Hubble patches.
"""
import json
import os
import time

import mpmath as mp
import numpy as np
from scipy.optimize import brentq

import smpot as S

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, '..', 'results')
HS = {'H_free-X 1.29e11': 1.29e11, 'H_M1 1.3428e11': 1.342842735e11, 'H_Q-matter 1.51e11': 1.51e11}
MTS = (172.33, 172.60, 172.87)
NS = (50, 55, 60)
LN49 = np.log(4 * np.pi / 9)


def B_HM(sm, H, xi0, mu0=None, cR=1.0, loop=True):
    hb, Vb, _, _ = S.barrier(sm, H, xi0, cR=cR, mu0=mu0, loop=loop)
    if hb is None:
        return None, None, 0.0
    return hb, Vb, 8 * np.pi ** 2 * Vb / (3 * H ** 4)


def xi_min(sm, H, N, mu0=None, cR=1.0, loop=True, lo=-0.2, hi=0.3):
    target = 3 * N + LN49

    def f(x):
        return B_HM(sm, H, x, mu0=mu0, cR=cR, loop=loop)[2] - target
    xs = np.linspace(lo, hi, 51)
    fs = [f(x) for x in xs]
    for i in range(len(xs) - 1):
        if fs[i] < 0 <= fs[i + 1]:
            return float(brentq(f, xs[i], xs[i + 1], xtol=1e-6))
    return None


def mass2_curv(sm, H, xi0, mu0=None, cR=1.0):
    h = 1e-3 * H
    return float(2 * sm.V_dS(np.array([h]), H, xi0, cR=cR, mu0=mu0)[0] / h ** 2)


def ds_variance(m2_over_H2, eps=1.0, N=60.0):
    """Super-horizon variance per real component (units H^2) of a Bunch-Davies scalar with mass m, modes with
    eps e^{-N} < k/(aH) < eps.  P(x) = x^3 |H_nu(x)|^2 / (8 pi)."""
    nu = mp.sqrt(mp.mpf(9) / 4 - m2_over_H2)

    def P(lx):
        x = mp.e ** lx
        return x ** 3 * abs(mp.hankel1(nu, x)) ** 2 * mp.e ** (-mp.pi * mp.im(nu)) / (8 * mp.pi)
    return float(mp.quad(P, [mp.log(eps) - N, mp.log(eps) - 10, mp.log(eps) - 3, mp.log(eps)]))


def tail_chi4(hc, sig2):
    x = hc ** 2 / (2 * sig2)
    return x, float(mp.e ** (-x) * (1 + x))


def fp_tail(sm, H, xi0, hc, hb, n=4000):
    """Fokker-Planck equilibrium (radial, 4 components): P(h) ~ h^3 exp(-8pi^2 V/3H^4) on [0, h_b];
    returns log10 of the probability of h > hc."""
    hs = np.linspace(1e-6 * H, hb, n)
    V = sm.V_dS(hs, H, xi0)
    w = 3 * np.log(hs) - 8 * np.pi ** 2 * V / (3 * H ** 4)
    wm = w.max()
    p = np.exp(w - wm)
    tot = np.trapezoid(p, hs)
    m = hs >= hc
    tail = np.trapezoid(p[m], hs[m]) if m.sum() > 1 else 0.0
    if tail <= 0:   # analytic continuation beyond float range
        j = np.searchsorted(hs, hc)
        return float((w[j] - wm + np.log(H) - np.log(tot)) / np.log(10))
    return float(np.log10(tail / tot))


def main():
    t0 = time.time()
    out = {'conventions': {
        'V': 'Delta V = V(h) - V(0) at the same mu, 1804.02020 de Sitter one loop + 3-loop running; mu^2 = h^2 + R',
        'xi_EW': 'xi(M_t) = 0.1 is the main reading of "xi_H = 0.1"; xi(H_i) = 0.1 is the alternative',
        'survival': '<N_bubbles> = (4pi/9) e^{3N - B} < 1; HM prefactor taken as 1 (H^4)',
        'N': 'N = e-folds between the exit of our present horizon and the end of inflation (50-60)'}}
    flat = {}
    for Mt in MTS:
        sm = S.SM(Mt)
        flat[Mt] = sm.flat_scales()
        rows = {}
        for hl, H in HS.items():
            r = {'H': H, 'flat': flat[Mt], 'h_max_over_H': flat[Mt]['h_max'] / H}
            r['xi_at_mu'] = {f'{mu:.3g}': float(sm.xi(mu, 0.1)) for mu in (S.MT0, 1e10, H, np.sqrt(12) * H, 1e12, 1e13)}
            for lab, mu0 in (('xiEW=0.1', None), ('xi(H_i)=0.1', H)):
                d = {}
                for cR in (1.0, 0.0):
                    for loop in (True, False):
                        hb, Vb, B = B_HM(sm, H, 0.1, mu0=mu0, cR=cR, loop=loop)
                        d[f'cR={cR},loop={loop}'] = dict(h_b_GeV=hb, h_b_over_H=(hb / H if hb else None),
                                                         dV_over_H4=(Vb / H ** 4 if Vb else None), B_HM=B,
                                                         log10_N_bubbles={N: float((LN49 + 3 * N - B) / np.log(10))
                                                                          for N in NS})
                m2 = mass2_curv(sm, H, 0.1, mu0=mu0)
                d['m2_curv_over_H2_one_loop'] = m2 / H ** 2
                d['m2_tree_12xi(mu=sqrtR)'] = 12 * float(sm.xi(np.sqrt(12) * H, 0.1, mu0))
                rows[lab] = d
            # xi_min
            r['xi_min'] = {}
            for N in NS:
                r['xi_min'][N] = {
                    'xi(M_t), one loop dS': xi_min(sm, H, N),
                    'xi(M_t), tree RGI': xi_min(sm, H, N, loop=False),
                    'xi(M_t), one loop, mu^2=h^2': xi_min(sm, H, N, cR=0.0),
                    'xi(H_i), one loop dS': xi_min(sm, H, N, mu0=H),
                    'xi(H_i), tree RGI': xi_min(sm, H, N, mu0=H, loop=False)}
            # variances and tails at xi_EW = 0.1
            m2 = rows['xiEW=0.1']['m2_curv_over_H2_one_loop']
            var = {'m2_over_H2': m2, 'stochastic_light_3H^4/(8pi^2 m^2)': 3 / (8 * np.pi ** 2 * m2)}
            for eps in (1.0, 0.3, 0.1):
                var[f'dS_superhorizon_eps={eps}'] = ds_variance(m2, eps)
            var['dS_superhorizon_eps=1_N=50'] = ds_variance(m2, 1.0, N=50)
            hb = rows['xiEW=0.1']['cR=1.0,loop=True']['h_b_GeV']
            tails = {}
            for vl in ('stochastic_light_3H^4/(8pi^2 m^2)', 'dS_superhorizon_eps=1.0'):
                s2 = var[vl] * H ** 2
                for cl, hc in (('flat h_max', flat[Mt]['h_max']), ('dS barrier h_b', hb)):
                    x, p = tail_chi4(hc, s2)
                    tails[f'{vl} | {cl}'] = dict(x=x, log10_P_patch=float(mp.log10(p)),
                                                 log10_expected_patches={N: float(mp.log10(p)) + 3 * N / np.log(10)
                                                                         for N in NS})
            tails['Fokker-Planck full V | flat h_max'] = {'log10_P_patch': fp_tail(sm, H, 0.1, flat[Mt]['h_max'], hb)}
            lp = tails['Fokker-Planck full V | flat h_max']['log10_P_patch']
            tails['Fokker-Planck full V | flat h_max']['log10_expected_patches'] = {N: lp + 3 * N / np.log(10)
                                                                                     for N in NS}
            # xi needed so that even the flat h_max is not reached by any of e^{3N} patches (stochastic light
            # formula, tree mass 12 xi H^2), and the exact heavy-field statement
            need = {}
            for N in NS:
                need[N] = float((3 * N + 5) * H ** 2 / (16 * np.pi ** 2 * flat[Mt]['h_max'] ** 2))
            r['variance_units_H2'] = var
            r['tails_at_end_of_inflation'] = tails
            r['xi(H) for no patch above flat h_max (light formula, x = 3N+5)'] = need
            r.update(rows)
            rows_all = out.setdefault(f'Mt={Mt}', {})
            rows_all[hl] = r
            print(Mt, hl, 'B', rows['xiEW=0.1']['cR=1.0,loop=True']['B_HM'], 'ximin', r['xi_min'][60], time.time() - t0,
                  flush=True)
    # heavy-field variance for xi above 3/16
    out['heavy_check'] = {f'm2={m2}': ds_variance(m2, 1.0) for m2 in (2.0, 2.2, 2.4, 3.0, 4.0)}
    out['seconds'] = time.time() - t0
    json.dump(out, open(os.path.join(RES, 'inflation.json'), 'w'), indent=1, default=float)


if __name__ == '__main__':
    main()
