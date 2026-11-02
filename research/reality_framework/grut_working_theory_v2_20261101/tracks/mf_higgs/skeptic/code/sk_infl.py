#!/usr/bin/env python3
"""mf_higgs SKEPTIC: the Higgs during inflation on HFa, recomputed with my own barrier code and four RG/potential
schemes (smpot couplings read-only).  Criterion choices (each a model choice):
  * constant-H Hawking-Moss count <N> = (4 pi/9) e^{3 N0 - B}, B = 8 pi^2 DeltaV_b / 3 H^4 (track's baseline), and
    the same with the past-light-cone radius eta_0 = 3.2/(a0 H0) instead of 1/(a0 H0) (x 3.2^3);
  * instantaneous HM through the waterfall to eps_H = 1 on the bench H(t), R(t) (pessimistic bracket).
Schemes for DeltaV_b:
  A  tree, mu^2 = h^2 + R, V = xi(mu) R phi^2/2 + lam(mu) phi^4/4, phi = e^{Gamma(mu)} h   (the track's choice)
  B  flat 2-loop lambda_eff(h) (Landau gauge, mu = h) + xi(h) R e^{2 Gamma(h)} h^2 / 2
  C  one-loop de Sitter potential of 1804.02020 (smpot.V_dS), mu^2 = h^2 + R
  D  same, mu = h (cR = 0)
Usage: python3 -B sk_infl.py [main|cap|inst]
"""
import sys
sys.dont_write_bytecode = True
import math
import time

import numpy as np
from scipy.optimize import brentq

import sk_common as C

LN10 = math.log(10)
H_I = 1.30002426e11
N0_EXP = 52.644453 + C.LN_K_H0 + 0.192656   # N_obs(K = 2.048) + t(eps_H = 1)
N0_PAULI = N0_EXP - 0.5 * math.log(7.054 / 2.048)


def Vfun(s, H, R, xi_EW, scheme):
    def V(h):
        h = np.atleast_1d(h)
        if scheme == 'A':
            mu = np.sqrt(h ** 2 + R)
            c = s.cpl(mu)
            xi = s.xi(mu, xi_EW)
            phi = np.exp(c['Gam']) * h
            return 0.5 * xi * R * phi ** 2 + 0.25 * c['lam'] * phi ** 4
        if scheme == 'B':
            c = s.cpl(h)
            xi = s.xi(h, xi_EW)
            le = s.lam_eff_flat(h)
            return 0.25 * le * h ** 4 + 0.5 * xi * R * np.exp(2 * c['Gam']) * h ** 2
        if scheme in ('C', 'D'):
            return s.V_dS(h, math.sqrt(R / 12), xi_EW, cR=1.0 if scheme == 'C' else 0.0)
    return V


def Bexp(s, H, R, xi_EW, scheme, n=500):
    V = Vfun(s, H, R, xi_EW, scheme)
    hs = np.exp(np.linspace(math.log(0.05 * H), math.log(3e16), n))
    v = V(hs)
    k = int(np.argmax(v))
    if k == 0 or v[k] <= 0:
        return 0.0, None
    if k == n - 1:
        return np.inf, None
    hh = np.exp(np.linspace(math.log(hs[k - 1]), math.log(hs[k + 1]), 121))
    vv = V(hh)
    j = int(np.argmax(vv))
    return 8 * math.pi ** 2 * float(vv[j]) / (3 * H ** 4), float(hh[j])


def log10N_const(Mt, xi_EW, scheme, N0=N0_EXP, H=H_I, lightcone=False, als=C.PDG['als']):
    s = C.sm(Mt, als=als)
    B, _ = Bexp(s, H, 12 * H ** 2, xi_EW, scheme)
    extra = 3 * math.log(3.2) if lightcone else 0.0
    return (math.log(4 * math.pi / 9) + 3 * N0 - B + extra) / LN10


def xi_min(f, lo=0.02, hi=0.25):
    a, b = f(lo), f(hi)
    if a > 0 > b:
        return float(brentq(f, lo, hi, xtol=2e-4))
    return None if a <= 0 else f'> {hi}'


def cmd_main():
    out = dict(N0_EXP=N0_EXP, N0_PAULI=N0_PAULI, H_i=H_I, rows={})
    t0 = time.time()
    for Mt in (172.57, 172.86, 173.45, 174.05, 174.6):
        r = {}
        for sc in ('A', 'B', 'C', 'D'):
            r[f'{sc} log10N(0.1)'] = log10N_const(Mt, 0.1, sc)
            r[f'{sc} xi_min'] = xi_min(lambda x: log10N_const(Mt, x, sc))
        r['A xi_min lightcone'] = xi_min(lambda x: log10N_const(Mt, x, 'A', lightcone=True))
        r['A xi_min N0_PAULI'] = xi_min(lambda x: log10N_const(Mt, x, 'A', N0=N0_PAULI))
        out['rows'][f'{Mt}'] = r
        print(Mt, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()}, '%.0fs' % (time.time() - t0),
              flush=True)
        C.jdump(out, C.RES / 'sk_infl.json')


def cmd_cap():
    """M_t at which log10<N>(xi_EW = 0.1) = 0, constant H, schemes A-D; alpha_s +- 0.0009 (scheme A).
    M_t > 174.6 GeV is an extrapolation of the SMDR fit grid (quadratic in M_t)."""
    import json
    out = json.load(open(C.RES / 'sk_infl.json'))
    caps = {}
    for sc in ('A', 'B', 'C', 'D'):
        f = lambda m: log10N_const(round(m, 3), 0.1, sc)
        a, b = f(174.6), f(177.0)
        caps[sc] = dict(log10N_174p6=a, log10N_177=b,
                        cap=(float(brentq(f, 174.6, 177.0, xtol=5e-3)) if a < 0 < b else ('> 177' if b < 0 else '< 174.6')))
        print(sc, caps[sc], flush=True)
    for als in (0.1171, 0.1189):
        f = lambda m: log10N_const(round(m, 3), 0.1, 'A', als=als)
        a, b = f(173.5), f(177.0)
        caps[f'A als={als}'] = dict(log10N_173p5=a, log10N_177=b,
                                    cap=(float(brentq(f, 173.5, 177.0, xtol=5e-3)) if a < 0 < b else ('> 177' if b < 0 else '< 173.5')))
        print(als, caps[f'A als={als}'], flush=True)
    out['caps_const_H'] = caps
    C.jdump(out, C.RES / 'sk_infl.json')


def inst_log10N(Mt, sc='A', name='HFa', als=C.PDG['als'], endpoint=True, stride=1):
    """instantaneous HM through the waterfall to eps_H = 1 (bench H(t), R(t)); endpoint=True adds the exact eps_H = 1
    row (R = 6 H^2), which the row table only approaches to within 1e-3 H^-1."""
    d, h = C.load_hist(C.TRACKS / 'mf_hybrid' / 'results' / f'bench_{name}.json')
    te = d['waterfall']['t_eps1']
    m = (h['t'] >= 0.12) & (h['t'] <= te + 1e-9)
    t, N, Hh, R = h['t'][m][::stride], h['N'][m][::stride], h['H'][m][::stride], h['R'][m][::stride]
    He = d['waterfall']['H_eps1_over_Hi']
    Ne = d['waterfall']['N_eps1']
    if endpoint:
        t = np.append(t, te); N = np.append(N, Ne); Hh = np.append(Hh, He); R = np.append(R, 6 * He ** 2)
    s = C.sm(Mt, als=als)
    Bs = np.array([Bexp(s, hh * H_I, r * H_I ** 2, 0.1, sc, n=300)[0] for hh, r in zip(Hh, R)])
    B12 = Bexp(s, H_I, 12 * H_I ** 2, 0.1, sc)[0]
    Nb = Ne - N
    lnI = 3 * (N0_EXP - Nb) + 3 * np.log(Hh) - Bs
    o = np.argsort(Nb)
    x, y = Nb[o], lnI[o]
    mx = y.max()
    part = mx + math.log(np.trapezoid(np.exp(y - mx), x))
    tail = 3 * (N0_EXP - x[-1]) - B12 - math.log(3)
    tot = (np.logaddexp(part, tail) + math.log(4 * math.pi / 3)) / LN10
    return float(tot), dict(B_last=float(Bs[-1]), R_last_over_H2=float(R[-1] / Hh[-1] ** 2), B_const=float(B12),
                            t_last=float(t[-1]))


def cmd_inst():
    """instantaneous bracket (scheme A, the track's): with and without the eps_H = 1 endpoint; M_t cap; alpha_s; HFb."""
    import json
    out = json.load(open(C.RES / 'sk_infl.json'))
    res = {}
    for Mt in (172.57, 172.86, 173.45, 174.05):
        for ep in (True, False):
            v, info = inst_log10N(Mt, endpoint=ep, stride=1 if ep else 4)
            res[f'{Mt} endpoint={ep}'] = dict(log10N=v, **info)
            print(Mt, ep, round(v, 2), info, flush=True)
    caps = {}
    for name, als in (('HFa', 0.1180), ('HFa', 0.1171), ('HFa', 0.1189), ('HFb', 0.1180)):
        f = lambda m: inst_log10N(round(m, 3), name=name, als=als, stride=2)[0]
        lo, hi = 172.9, 174.6
        a, b = f(lo), f(hi)
        caps[f'{name} als={als}'] = float(brentq(f, lo, hi, xtol=1e-2)) if a < 0 < b else (a, b)
        print(name, als, caps[f'{name} als={als}'], flush=True)
        out['instantaneous'] = dict(rows=res, caps_with_endpoint=caps)
        C.jdump(out, C.RES / 'sk_infl.json')


if __name__ == '__main__':
    dict(main=cmd_main, cap=cmd_cap, inst=cmd_inst)[sys.argv[1] if len(sys.argv) > 1 else 'main']()
