#!/usr/bin/env python3
"""Skeptic of inf_higgs: the inflationary Hawking-Moss bubble count with the benchmark's own H(N), R(N), re-done with
own background (sk_common.initial_state, LSODA) and own barrier search on the SM splines, and extended in M_t
(the track evaluated it only at 172.28 / 172.57 / 172.86 but claims survival to 174.6 GeV with a plasma; a plasma
cannot act during inflation).
Criterion (the track's, Mantziris-Markkanen-Rajantie 2011.03763 convention):
    <N> = (4 pi/3) int dN_b e^{3(N0 - N_b)} (H/H0)^3 e^{-B(N_b)},  B = 8 pi^2 DeltaV_b / (3 H^4),
tree RG-improved V = xi(mu) R phi^2/2 + m_x^2 phi^2/2 + lam(mu) phi^4/4, phi = e^{Gamma(mu)} h, mu^2 = h^2 + R.
  python3 -B sk_infl.py main          -> results/sk_infl.json (xi_EW,min vs M_t, B_end, B_const)
  python3 -B sk_infl.py induced       -> results/sk_induced.json (S1 loop-induced Higgs mass along the trajectory)
"""
import sys
sys.dont_write_bytecode = True
import json
import math

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

import sk_common as C

MTS = [172.28, 172.57, 172.86, 173.15, 173.45, 173.75, 174.05, 174.6]
H0C_LN = math.log(0.05 / 2.248e-4)


def history(name, iH=0):
    pot = C.Pot(name, iH)
    hist = pot.rec['history']
    row = max(hist, key=lambda r: r['N_before_end'] if r['N_before_end'] <= 14 else -1)
    f0 = row['phi_MP']
    H0 = math.sqrt(pot.V(f0) / 3)
    fd0 = -pot.dV(f0) / (3 * H0)

    def rhs(t, y):
        H = math.sqrt((0.5 * y[2] ** 2 + pot.V(y[1])) / 3)
        return [H, y[2], -3 * H * y[2] - pot.dV(y[1])]

    def ev(t, y):
        H2 = (0.5 * y[2] ** 2 + pot.V(y[1])) / 3
        return 0.5 * y[2] ** 2 / H2 - 1
    ev.terminal = True; ev.direction = 1
    s = solve_ivp(rhs, (0, 1e4), [0, f0, fd0], method='LSODA', rtol=1e-11, atol=1e-14, dense_output=True, events=ev)
    te = s.t_events[0][0]
    lae = s.y_events[0][0][0]
    # sample in N_b = lae - lna
    Nb = np.unique(np.concatenate([np.linspace(0, 0.3, 31), np.linspace(0.3, 1.5, 25), np.linspace(1.5, lae - 0.5, 20)]))
    out = []
    for x in Nb:
        t = brentq(lambda tt: s.sol(tt)[0] - (lae - x), 0, te) if x > 0 else te
        la, f, fd = s.sol(t)
        H = math.sqrt((0.5 * fd * fd + pot.V(f)) / 3)
        eps = 0.5 * fd * fd / H ** 2
        out.append((x, H, 6 * H * H * (2 - eps), eps, f - pot.phimin))
    return pot, np.array(out)


def B_tree(smi, H, R, xi_EW, m2x=0.0):
    """Hawking-Moss exponent (H, R in GeV, GeV^2; m2x canonical extra mass^2 in GeV^2)"""
    hs = np.exp(np.linspace(math.log(1e-3 * H), math.log(1e17), 1500))
    lmu = 0.5 * np.log(hs ** 2 + R)
    lam = smi.sp['lam'](lmu)
    G = smi.sp['Gam'](lmu)
    Ixi = smi.sp['Ixi'](lmu)
    xi = 1 / 6 + (xi_EW - 1 / 6) * np.exp(Ixi - smi.d['Ixi_Mt'])
    phi = np.exp(G) * hs
    V = 0.5 * (xi * R + m2x) * phi ** 2 + 0.25 * lam * phi ** 4
    k = int(np.argmax(V))
    if k == 0 or V[k] <= 0:
        return 0.0
    if k == len(hs) - 1:
        return np.inf
    # refine
    a, b = hs[k - 1], hs[k + 1]
    hh = np.exp(np.linspace(math.log(a), math.log(b), 200))
    lmu = 0.5 * np.log(hh ** 2 + R)
    xi = 1 / 6 + (xi_EW - 1 / 6) * np.exp(smi.sp['Ixi'](lmu) - smi.d['Ixi_Mt'])
    phi = np.exp(smi.sp['Gam'](lmu)) * hh
    VV = 0.5 * (xi * R + m2x) * phi ** 2 + 0.25 * smi.sp['lam'](lmu) * phi ** 4
    return 8 * math.pi ** 2 * float(VV.max()) / (3 * H ** 4)


def log10N(smi, pot, hist, xi_EW, N0, m2x_fn=None):
    Hg = pot.H_GeV
    Nb, H, R, eps, fI = hist.T
    Bs = np.array([B_tree(smi, H[i] * Hg, R[i] * Hg ** 2, xi_EW,
                          (m2x_fn(fI[i]) * Hg ** 2 if m2x_fn else 0.0)) for i in range(len(Nb))])
    Hp = H[-1]
    lnI = 3 * (N0 - Nb) + 3 * np.log(H / Hp) - Bs
    m = lnI.max()
    part = m + math.log(np.trapezoid(np.exp(lnI - m), Nb))
    tail = 3 * (N0 - Nb[-1]) - Bs[-1] - math.log(3.0)
    tot = np.logaddexp(part, tail) + math.log(4 * math.pi / 3)
    return float(tot / C.LN10), Bs


def N0_MD(name, iH=0, T_RH=1e11):
    r = json.load(open(C.TRACKS / 'inf_potential' / 'results' / f'bench_{name}.json'))['by_Hi'][iH]
    inst = [x for x in r['Nstar_alternatives'] if x['history'] == ['inst']][0]['Nstar']
    return inst - math.log(r['T_inst_GeV'] / T_RH) / 3 + H0C_LN


def xi_min(f, lo=0.0, hi=0.4):
    xs = np.linspace(lo, hi, 9)
    fs = [f(x) for x in xs]
    for i in range(len(xs) - 1):
        if fs[i] > 0 >= fs[i + 1]:
            return float(brentq(f, xs[i], xs[i + 1], xtol=2e-4))
    return None


def cmd_main():
    out = {}
    for name in ('HT6t', 'T1'):
        pot, hist = history(name)
        N0 = N0_MD(name)
        o = dict(N0=N0, R_over_H2_end=float(hist[0, 2] / hist[0, 1] ** 2), H_end=float(hist[0, 1]))
        for Mt in MTS:
            smi = C.sm_in(Mt)
            l, Bs = log10N(smi, pot, hist, 0.1, N0)
            Bc = B_tree(smi, pot.H_GeV, 12 * pot.H_GeV ** 2, 0.1)
            xm = xi_min(lambda x: log10N(smi, pot, hist, x, N0)[0])
            xc = xi_min(lambda x: (math.log(4 * math.pi / 9) + 3 * N0 - B_tree(smi, pot.H_GeV, 12 * pot.H_GeV ** 2, x)) / C.LN10)
            o[str(Mt)] = dict(log10N_HN=l, B_end=float(Bs[0]), B_const=Bc, xi_EW_min_HN=xm, xi_EW_min_constH=xc,
                              h_max=smi.h_max)
            print(name, Mt, 'log10N %.1f B_end %.0f B_const %.0f xi_min H(N) %s const %s' % (l, Bs[0], Bc, xm, xc), flush=True)
        out[name] = o
    C.jdump(out, C.RES / 'sk_infl.json')


if __name__ == '__main__':
    if sys.argv[1] == 'main':
        cmd_main()


def B_loop(sm, H, R, xi_EW):
    """one loop: smpot.SM.V_dS (1804.02020) at H_eff = (R/12)^{1/2}, the track's choice; noise with the true H"""
    Heff = math.sqrt(R / 12)
    hs = np.exp(np.linspace(math.log(1e-3 * H), math.log(1e17), 900))
    V = sm.V_dS(hs, Heff, xi_EW)
    k = int(np.argmax(V))
    if k == 0 or V[k] <= 0:
        return 0.0
    if k == len(hs) - 1:
        return np.inf
    hh = np.exp(np.linspace(math.log(hs[k - 1]), math.log(hs[k + 1]), 120))
    return 8 * math.pi ** 2 * float(np.max(sm.V_dS(hh, Heff, xi_EW))) / (3 * H ** 4)


def cmd_loop():
    out = json.load(open(C.RES / 'sk_infl.json'))
    name = 'HT6t'
    pot, hist = history(name)
    N0 = N0_MD(name)
    Hg = pot.H_GeV
    sel = hist[(hist[:, 0] < 1.6) | (np.arange(len(hist)) == len(hist) - 1)]
    for Mt in (172.86, 173.45, 173.75, 174.05, 174.6):
        sm = C.S.SM(Mt, Mh=125.20, als=0.1180)

        def l10(x):
            Bs = np.array([B_loop(sm, r[1] * Hg, r[2] * Hg ** 2, x) for r in sel])
            lnI = 3 * (N0 - sel[:, 0]) + 3 * np.log(sel[:, 1] / sel[-1, 1]) - Bs
            m = lnI[:-1].max()
            part = m + math.log(np.trapezoid(np.exp(lnI[:-1] - m), sel[:-1, 0]))
            tail = 3 * (N0 - sel[-2, 0]) - Bs[-1] - math.log(3.0)
            return float((np.logaddexp(part, tail) + math.log(4 * math.pi / 3)) / C.LN10), Bs
        l, Bs = l10(0.1)
        xm = xi_min(lambda x: l10(x)[0], 0.0, 0.2)
        out[name][str(Mt)]['loop'] = dict(log10N_HN=l, B_end=float(Bs[0]), xi_EW_min_HN=xm)
        print(name, Mt, 'loop log10N %.1f B_end %.0f xi_min %s' % (l, Bs[0], xm), flush=True)
        C.jdump(out, C.RES / 'sk_infl.json')


if __name__ == '__main__' and sys.argv[1] == 'loop':
    cmd_loop()


def cmd_Abound():
    """inflation-only Hawking-Moss bound on the destabilising A (track R10), own code, tree"""
    out = json.load(open(C.RES / 'sk_infl.json'))
    for name in ('HT6t', 'T1', 'E1'):
        pot, hist = history(name)
        N0 = N0_MD(name)
        sgn = -1.0 if hist[-1, 4] > 0 else 1.0
        for Mt in (172.57, 172.86):
            smi = C.sm_in(Mt)
            f = lambda lgA: log10N(smi, pot, hist, 0.1, N0,
                                   m2x_fn=lambda fI, lgA=lgA: sgn * 10 ** lgA * fI * C.MPL / pot.H_GeV ** 2)[0]
            lg = brentq(f, 1.0, 6.0, xtol=1e-3)
            out.setdefault(name, {}).setdefault(str(Mt), {})['A_max_inflation_TeV'] = 10 ** lg / 1e3
            print(name, Mt, 'A_max(inflation, tree) %.1f TeV' % (10 ** lg / 1e3), flush=True)
    C.jdump(out, C.RES / 'sk_infl.json')


if __name__ == '__main__' and sys.argv[1] == 'Abound':
    cmd_Abound()
