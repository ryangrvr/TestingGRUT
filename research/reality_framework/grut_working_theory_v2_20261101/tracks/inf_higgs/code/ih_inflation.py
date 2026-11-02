#!/usr/bin/env python3
"""Track inf_higgs: the Higgs during inflation with the benchmark's own H(N) and R(N) (own code).

higgs_reheat/code/inflation.py assumed H = H_i and R = 12 H^2 over the last 60 e-folds.  On the benchmarks H(N) is
constant to 3e-4 until the last e-fold, but R/H^2 = 6(2 - eps_H) falls from 12 to 6 at eps_H = 1, and H falls by 1.3%
(T1, E1) or 4.3% (HT6t).  The curvature mass, the de Sitter barrier and the Hawking-Moss exponent fall with it.

Criterion (the literature's: Mantziris-Markkanen-Rajantie arXiv:2011.03763 [Q abstract], who find decay "most likely a
few e-foldings before the end of inflation"): expected number of AdS bubbles in our past light cone
    <N> = (4 pi/3) int_0^{N_0} dN_b e^{3(N_0 - N_b)} (H(N_b)/H(N_0))^3 e^{-B(N_b)},   B = 8 pi^2 DeltaV_b/(3 H^4),
N_b = e-folds before eps_H = 1, N_0 = exit of today's horizon (N_obs + 1).  For constant B it is (4 pi/9) e^{3 N_0 - B},
higgs_reheat's form.  DeltaV_b is the barrier of
    tree:     V = xi(mu) R phi^2/2 + lam(mu) phi^4/4, phi = e^{Gamma(mu)} h, mu^2 = h^2 + R   (RG-improved tree);
    one loop: smpot.SM.V_dS (1804.02020) evaluated at the quasi-de Sitter H_eff = (R/12)^{1/2}  [model choice: the
              de Sitter one-loop potential at the instantaneous curvature; the noise uses the true H].
Extra Higgs mass terms (model variants): R-odd A phi_I |H|^2 (m^2 += A phi_I), R-even portal (lam_IH/2) phi_I^2 |H|^2
(m^2 += lam_IH phi_I^2/2).  Writes ../results/inflation_HN.json.
Also the non-equilibrium Gaussian tail against the instantaneous de Sitter barrier, from the mode runs (ih_run cache).
"""
import sys
sys.dont_write_bytecode = True
import json
import math
import time

import numpy as np
from scipy.optimize import brentq

import ih_common as C
import ih_run as IR

S = C.S
IC = C.IC
LN43 = math.log(4 * math.pi / 3)
_SMS = {}


def SMfull(Mt, Mh=C.PDG24['Mh']):
    k = (round(Mt, 4), Mh)
    if k not in _SMS:
        _SMS[k] = S.SM(Mt, Mh=Mh, als=C.PDG24['als'])
    return _SMS[k]


def history(name, iH=0, Nmax=4.0, n=60):
    """H(N_b), R(N_b), phi_I(N_b) in GeV units for N_b in [0, Nmax] (dense near the end) and the plateau value."""
    B = C.load_bench(name, iH)
    m, A = B['model'], B['A']
    bg = IC.background(m, Nwant=75.0)
    Nb = np.unique(np.concatenate([np.linspace(0, 0.3, 16), np.linspace(0.3, 1.2, 13), np.linspace(1.2, Nmax, 8)]))
    rows = []
    for x in Nb:
        f, fN = IC.bg_at(bg, x)
        eps = 0.5 * fN * fN
        H = math.sqrt(A * m.v(f) / (3 - eps))
        rows.append((x, H * B['H_GeV'], 6 * H * H * (2 - eps) * B['H_GeV'] ** 2, eps, f - B['phimin']))
    f, fN = IC.bg_at(bg, 50.0)
    eps = 0.5 * fN * fN
    Hp = math.sqrt(A * m.v(f) / (3 - eps)) * B['H_GeV']
    return dict(rows=np.array(rows), H_plateau=Hp, phiI_plateau=f - B['phimin'], H_GeV=B['H_GeV'])


def barrier_B(sm, H, R, xi_EW, loop, m2x=0.0, n=900):
    """Hawking-Moss exponent at (H, R) (GeV); m2x (GeV^2) extra Higgs mass^2 (canonical)."""
    hs = np.exp(np.linspace(math.log(1e-3 * H), math.log(1e17), n))

    def V(h):
        if loop:
            Heff = math.sqrt(R / 12)
            v = sm.V_dS(h, Heff, xi_EW)
        else:
            mu = np.sqrt(h ** 2 + R)
            c = sm.cpl(mu)
            xi = sm.xi(mu, xi_EW)
            phi = np.exp(c['Gam']) * h
            v = 0.5 * xi * R * phi ** 2 + 0.25 * c['lam'] * phi ** 4
        if m2x:
            c = sm.cpl(np.sqrt(h ** 2 + R))
            v = v + 0.5 * m2x * (np.exp(c['Gam']) * h) ** 2
        return v
    Vv = V(hs)
    k = int(np.argmax(Vv))
    if k == 0 or Vv[k] <= 0:
        return 0.0, None
    if k == n - 1:
        return np.inf, None
    hh = np.exp(np.linspace(math.log(hs[k - 1]), math.log(hs[k + 1]), 201))
    VV = V(hh)
    j = int(np.argmax(VV))
    return 8 * math.pi ** 2 * float(VV[j]) / (3 * H ** 4), float(hh[j])


def log10_bubbles(hist, sm, xi_EW, N0, loop, A=0.0, lamIH=0.0, Nmax_int=None):
    """log10 <N> with the history; A (GeV) R-odd trilinear, lamIH portal.  Beyond the tabulated range B is the
    plateau value (constant to 3e-4)."""
    r = hist['rows']
    Nb, H, R, eps, fI = r[:, 0], r[:, 1], r[:, 2], r[:, 3], r[:, 4]
    MP = C.MPL
    Bs = np.empty(len(Nb))
    for i in range(len(Nb)):
        m2x = A * fI[i] * MP + lamIH * (fI[i] * MP) ** 2 / 2
        Bs[i] = barrier_B(sm, H[i], R[i], xi_EW, loop, m2x)[0]
    Hp = hist['H_plateau']
    m2p = A * hist['phiI_plateau'] * MP + lamIH * (hist['phiI_plateau'] * MP) ** 2 / 2
    Bp = barrier_B(sm, Hp, 12 * Hp ** 2, xi_EW, loop, m2p)[0]
    # integrand in ln: 3(N0 - Nb) + 3 ln(H/Hp) - B
    lnI = 3 * (N0 - Nb) + 3 * np.log(H / Hp) - Bs
    m = np.max(lnI)
    part = m + math.log(np.trapezoid(np.exp(lnI - m), Nb))
    # tail Nb > Nmax at the plateau value: int_{Nmax}^{N0} e^{3(N0-Nb) - Bp} = e^{3(N0-Nmax) - Bp}/3
    tail = 3 * (N0 - Nb[-1]) - Bp - math.log(3.0)
    tot = np.logaddexp(part, tail) + LN43
    imax = int(np.argmax(lnI))
    return float(tot / math.log(10)), dict(B_end=float(Bs[0]), B_plateau=float(Bp), B_at_max=float(Bs[imax]),
                                           Nb_at_max=float(Nb[imax]), B_min=float(np.min(Bs)), Nb_Bmin=float(Nb[int(np.argmin(Bs))]))


def log10_bubbles_const(sm, H, xi_EW, N0, loop):
    """higgs_reheat's constant-H de Sitter form (4 pi/9) e^{3 N0 - B}"""
    Bc = barrier_B(sm, H, 12 * H * H, xi_EW, loop)[0]
    return float((math.log(4 * math.pi / 9) + 3 * N0 - Bc) / math.log(10)), Bc


def xi_min(f, lo=0.0, hi=0.3):
    """smallest xi_EW with f(xi) <= 0 (f = log10 <N>, decreasing in xi)"""
    xs = np.linspace(lo, hi, 11)
    fs = [f(x) for x in xs]
    for i in range(len(xs) - 1):
        if fs[i] > 0 >= fs[i + 1]:
            return float(brentq(f, xs[i], xs[i + 1], xtol=2e-5))
    return None if fs[-1] > 0 else f'< {lo}'


def gaussian_tail(name, hist_key, iH, sm, xi_EW, N0):
    """non-equilibrium Gaussian count of Hubble patches beyond the instantaneous tree de Sitter barrier during the last
    e-folds, from the mode run: log10[e^{3(N0 - Nb)} (1 + x) e^{-x}], x = phi_b^2 / (2 sigma^2(k < aH))."""
    H_GeV = C.load_bench(name, iH)['H_GeV']
    xi = IR.smtab(172.57).xi(H_GeV, 0.1)
    R = IR.run_cached(name, hist_key, xi, iH)
    rows, spec, ks, wq = R['rows'], R['spec'], R['ks'], R['wq']
    te = R['res']['tau_endinf']
    out = []
    best = -np.inf
    for j in range(len(rows)):
        t, a, H, Rc = rows[j, 0], rows[j, 1], rows[j, 2], rows[j, 3]
        if t > te:
            break
        Nb = math.log(math.exp(1.0) / a) if True else 0.0      # a(S1) = 1, end at a = e (to 1e-3)
        m = ks < a * H
        s2 = float(np.sum(wq[m] * spec[m, j]))
        Bv, hb = barrier_B(sm, H * H_GeV, Rc * H_GeV ** 2, xi_EW, False)
        if hb is None or not np.isfinite(Bv):
            continue
        phib = math.exp(float(sm.cpl(hb)['Gam'])) * hb / H_GeV
        x = phib ** 2 / (2 * s2)
        lp = (3 * (N0 - Nb) + math.log1p(x) - x) / math.log(10)
        out.append((t, Nb, s2, phib, x, lp))
        best = max(best, lp)
    return dict(best_log10_patches=best, rows=out[::4])


def main():
    t0 = time.time()
    path = C.RES / 'inflation_HN.json'
    try:
        out = json.load(open(path))
    except Exception:
        out = {}
    for name in C.BENCH:
        for iH in (0, 1):
            key = f'{name}|iH{iH}'
            if key in out:
                continue
            hist = history(name, iH)
            H = hist['H_GeV']
            r = hist['rows']
            o = dict(H_i=H, H_end_over_Hi=float(r[0, 1] / H), H_plateau_over_Hi=hist['H_plateau'] / H,
                     R_over_H2_end=float(r[0, 2] / r[0, 1] ** 2), phiI_end=float(r[0, 4]), phiI_plateau=hist['phiI_plateau'])
            N0s = dict(MD=IR.N_obs(name, 'MD', iH) + 1, G05=IR.N_obs(name, 'G05', iH) + 1, conv60=60.0)
            o['N0'] = N0s
            for Mt in IR.MT_PDG:
                sm = SMfull(Mt)
                row = {}
                for loop in (False, True):
                    lab = 'loop' if loop else 'tree'
                    for nk, N0 in N0s.items():
                        lb, info = log10_bubbles(hist, sm, 0.1, N0, loop)
                        lc, Bc = log10_bubbles_const(sm, H, 0.1, N0, loop)
                        row[f'{lab} {nk}'] = dict(log10N_HN=lb, log10N_constH=lc, B_constH=Bc, **info)
                    for nk in ('MD', 'conv60'):
                        N0 = N0s[nk]
                        row[f'{lab} {nk} xi_min_HN'] = xi_min(lambda x: log10_bubbles(hist, sm, x, N0, loop)[0])
                        row[f'{lab} {nk} xi_min_constH'] = xi_min(lambda x: log10_bubbles_const(sm, H, x, N0, loop)[0])
                o[f'{Mt}'] = row
                print(name, iH, Mt, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in row.items() if 'xi_min' in k},
                      'B_end/B_const (tree MD) %.0f/%.0f' % (row['tree MD']['B_end'], row['tree MD']['B_constH']),
                      '%.0fs' % (time.time() - t0), flush=True)
            out[key] = o
            C.jdump(out, path)


def cmd_Aterm():
    """bound on an R-odd A phi_I |H|^2 (either sign) during inflation at xi_EW = 0.1, tree and loop, MD N0;
    and a check that the R-even portal only stabilises."""
    path = C.RES / 'inflation_Aterm.json'
    out = {}
    for name in C.BENCH:
        hist = history(name, 0)
        H = hist['H_GeV']
        N0 = IR.N_obs(name, 'MD', 0) + 1
        o = {}
        for Mt in IR.MT_PDG:
            sm = SMfull(Mt)
            row = {}
            for loop in (False, True):
                lab = 'loop' if loop else 'tree'
                # destabilising sign: A phi_I < 0
                sgn = -1.0 if hist['phiI_plateau'] > 0 else 1.0
                f = lambda lgA: log10_bubbles(hist, sm, 0.1, N0, loop, A=sgn * 10 ** lgA)[0]
                lo, hi = -3.0, 8.0
                if f(lo) > 0:
                    row[f'{lab} A_max_GeV'] = 0.0
                elif f(hi) <= 0:
                    row[f'{lab} A_max_GeV'] = f'> 1e{hi}'
                else:
                    lg = brentq(f, lo, hi, xtol=1e-3)
                    row[f'{lab} A_max_GeV'] = 10 ** lg
                    row[f'{lab} A_max_over_Hi'] = 10 ** lg / H
                    row[f'{lab} A_max_phiI_end_over_Hi2'] = 10 ** lg * abs(hist['rows'][0, 4]) * C.MPL / H ** 2
                    row[f'{lab} A_max_phiI_plateau_over_Hi2'] = 10 ** lg * abs(hist['phiI_plateau']) * C.MPL / H ** 2
                row[f'{lab} lamIH=1e-10 log10N'] = log10_bubbles(hist, sm, 0.1, N0, loop, lamIH=1e-10)[0]
                row[f'{lab} lamIH=0 log10N'] = log10_bubbles(hist, sm, 0.1, N0, loop)[0]
            o[f'{Mt}'] = row
            print(name, Mt, row, flush=True)
        out[name] = o
        C.jdump(out, path)


def cmd_tail():
    path = C.RES / 'inflation_gaussian_tail.json'
    out = {}
    for name in C.BENCH:
        hk = 'MD' if name == 'HT6t' else 'MDf'
        N0 = IR.N_obs(name, 'MD', 0) + 1
        o = {}
        for Mt in IR.MT_PDG:
            sm = SMfull(Mt)
            o[f'{Mt}'] = gaussian_tail(name, hk, 0, sm, 0.1, N0)
            print(name, Mt, o[f'{Mt}']['best_log10_patches'], flush=True)
        out[name] = o
        C.jdump(out, path)


if __name__ == '__main__':
    c = sys.argv[1] if len(sys.argv) > 1 else 'main'
    {'main': main, 'Aterm': cmd_Aterm, 'tail': cmd_tail}[c]()
