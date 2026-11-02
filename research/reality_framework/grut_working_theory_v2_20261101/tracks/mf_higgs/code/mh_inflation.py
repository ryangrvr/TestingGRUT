#!/usr/bin/env python3
"""Track mf_higgs: the Higgs during inflation on the hybrid benchmarks (own code; barrier from inf_higgs's
ih_inflation.barrier_B, imported read-only, on smpot.SM).

The hybrid is V0-dominated: over the whole observable inflation H = H_i to 2e-6 (H_* = 1.0000023 H_i) and R = 12 H^2 to
1e-7, until sigma_c.  The waterfall then ends inflation within 0.19 H^-1: R/H^2 falls from 12 to 6 (eps_H = 1) in the last
~0.02 H^-1 and swings negative afterwards.

Criteria (each a model choice):
  * HM-const [A, the baseline]: Hawking-Moss bubble count with constant H and R = 12 H^2 over N0 e-folds,
      <N> = (4 pi/9) e^{3 N0 - B},  B = 8 pi^2 DeltaV_b / (3 H^4),
    N0 = e-folds from the exit of today's horizon to eps_H = 1 (N_obs + 0.193).  HM needs the barrier to persist for
    the Higgs relaxation time ~ 3 H/m^2 ~ 3 H^-1; the sub-Hubble drop of R at the waterfall is outside its validity.
  * HM-inst [pessimistic bracket]: the literature's instantaneous form (Mantziris-Markkanen-Rajantie 2011.03763) taken
    literally through the waterfall to eps_H = 1 with the bench's H(t), R(t):
      <N> = (4 pi/3) int dN_b e^{3(N0 - N_b)} (H/H_0)^3 e^{-B(N_b)}.
  * Gaussian tail [non-equilibrium]: from the mode runs (mh_modes), log10 of the expected number of Hubble patches with
    the Higgs beyond the instantaneous tree de Sitter barrier, e^{3(N0 - N_b)} (1 + x) e^{-x}, x = phi_b^2/(2 sigma^2(k < aH)).
Potentials: tree RG-improved (mu^2 = h^2 + R) or the one-loop de Sitter potential of 1804.02020 (smpot.V_dS) at the
instantaneous H_eff = (R/12)^{1/2}.
"""
import sys
sys.dont_write_bytecode = True
import json
import math
import time

import numpy as np
from scipy.optimize import brentq

import mh_common as C
sys.path.insert(0, str(C.TRACKS / 'inf_higgs' / 'code'))
import ih_inflation as IHI        # noqa: E402  barrier_B (read-only)

LN10 = math.log(10.0)
_SM = {}


def SMfull(Mt, Mh=C.PDG24['Mh'], als=C.PDG24['als']):
    k = (round(Mt, 4), round(Mh, 3), round(als, 5))
    if k not in _SM:
        _SM[k] = C.S.SM(Mt, Mh=Mh, als=als)
    return _SM[k]


def B_of(sm, H, R, xi_EW, loop, m2x=0.0):
    return IHI.barrier_B(sm, H, R, xi_EW, loop, m2x)[0]


def log10N_const(sm, H_GeV, xi_EW, N0, loop, m2x=0.0):
    B = B_of(sm, H_GeV, 12 * H_GeV ** 2, xi_EW, loop, m2x)
    return (math.log(4 * math.pi / 9) + 3 * N0 - B) / LN10, B


def waterfall_table(name, n=40):
    """bench H(t), R(t) from t = 0 (sigma_c) to eps_H = 1; N_b = e-folds before eps_H = 1."""
    d = C.load_bench(name)
    h = d['h']
    te = d['t_eps1']
    m = (h['t'] >= 0.10) & (h['t'] <= te)
    t, N, H, R = h['t'][m], h['N'][m], h['H'][m], h['R'][m]
    Ne = float(np.interp(te, h['t'], h['N']))
    return dict(t=t, Nb=Ne - N, H=H, R=R, N_end=Ne)


def log10N_inst(sm, name, H_GeV, xi_EW, N0, loop, stride=4):
    """instantaneous HM through the waterfall (pessimistic bracket)."""
    W = waterfall_table(name)
    Nb, H, R = W['Nb'][::stride], W['H'][::stride], W['R'][::stride]
    Bs = np.array([B_of(sm, h * H_GeV, r * H_GeV ** 2, xi_EW, loop) if r > 0 else 0.0 for h, r in zip(H, R)])
    B12 = B_of(sm, H_GeV, 12 * H_GeV ** 2, xi_EW, loop)
    lnI = 3 * (N0 - Nb) + 3 * np.log(H) - Bs
    order = np.argsort(Nb)
    x, y = Nb[order], lnI[order]
    mx = np.max(y)
    part = mx + math.log(max(np.trapezoid(np.exp(y - mx), x), 1e-300))
    tail = 3 * (N0 - x[-1]) - B12 - math.log(3.0)
    tot = np.logaddexp(part, tail) + math.log(4 * math.pi / 3)
    return float(tot / LN10), dict(B12=B12, B_min=float(np.min(Bs)), R_at_Bmin=float(R[int(np.argmin(Bs))]))


def xi_min(f, lo=0.0, hi=0.4, n=17):
    xs = np.linspace(lo, hi, n)
    fs = [f(x) for x in xs]
    for i in range(len(xs) - 1):
        if fs[i] > 0 >= fs[i + 1]:
            return float(brentq(f, xs[i], xs[i + 1], xtol=2e-5))
    return None if fs[-1] > 0 else f'< {lo}'


def Mt_cap(f, lo=172.0, hi=175.5):
    """M_t at which log10<N>(xi_EW = 0.1) crosses 0 (f increasing in M_t)"""
    a, b = f(lo), f(hi)
    if a > 0:
        return f'< {lo}'
    if b <= 0:
        return f'> {hi}'
    return float(brentq(f, lo, hi, xtol=2e-3))


MT_LIST = (172.28, 172.57, 172.86, 173.15, 173.45, 173.75, 174.05, 174.35, 174.6)


def cmd_main():
    """xi_min (constant H, tree and loop; instantaneous bracket tree) at the M_t list, the M_t cap at xi_EW = 0.1,
    for HFa and HFb, with N0 for the exported K and the Pauli K = 7.05; alpha_s +- 0.0009; H_i dependence."""
    path = C.RES / 'inflation.json'
    try:
        out = json.load(open(path))
    except Exception:
        out = {}
    t0 = time.time()
    names = sys.argv[2].split(',') if len(sys.argv) > 2 else ('HFa', 'HFb')
    for name in names:
        d = C.load_bench(name)
        H = d['H_i_GeV']
        Hst = d['cmb']['ms']['H_star']
        W = waterfall_table(name)
        o = out.get(name) or dict(H_i=H, H_star_over_Hi=Hst, R_over_H2_pivot=12.0, eps_star=d['cmb']['ms']['eps_star'],
                                  N_eps1_after_slice=W['N_end'], rows={})
        N0s = dict(EXP=C.N_obs(name, d['reheating']['K']) + W['N_end'], PAULI0=C.N_obs(name, 7.054) + W['N_end'],
                   conv60=60.0)
        o['N0'] = N0s
        for Mt in MT_LIST:
            if f'{Mt}' in o['rows']:
                continue
            sm = SMfull(Mt)
            r = {}
            for lab, loop in (('tree', False), ('loop', True)):
                for nk in ('EXP', 'PAULI0', 'conv60'):
                    N0 = N0s[nk]
                    r[f'{lab} {nk} log10N xi0.1'] = log10N_const(sm, H, 0.1, N0, loop)[0]
                    r[f'{lab} {nk} xi_min'] = xi_min(lambda x: log10N_const(sm, H, x, N0, loop)[0])
                if not loop:
                    N0 = N0s['EXP']
                    r['tree EXP inst log10N xi0.1'] = log10N_inst(sm, name, H, 0.1, N0, False)[0]
                    r['tree EXP inst xi_min'] = xi_min(lambda x: log10N_inst(sm, name, H, x, N0, False)[0], hi=0.3, n=13)
            o['rows'][f'{Mt}'] = r
            out[name] = o
            C.jdump(out, path)
            print(name, Mt, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items() if 'xi_min' in k},
                  '%.0fs' % (time.time() - t0), flush=True)
        # M_t caps at xi_EW = 0.1
        if 'Mt_cap_xi0.1' in o:
            continue
        caps = {}
        for lab, loop in (('tree', False), ('loop', True)):
            N0 = N0s['EXP']
            caps[f'{lab} const'] = Mt_cap(lambda m: log10N_const(SMfull(round(m, 3)), H, 0.1, N0, loop)[0])
        caps['tree inst'] = Mt_cap(lambda m: log10N_inst(SMfull(round(m, 3)), name, H, 0.1, N0s['EXP'], False)[0],
                                   lo=172.5, hi=175.0)
        o['Mt_cap_xi0.1'] = caps
        print(name, 'caps', caps, flush=True)
        out[name] = o
        C.jdump(out, path)


def cmd_sens():
    """alpha_s +- 0.0009 and H_i dependence of the inflationary cap and xi_min (HFa, constant H; tree and loop)."""
    path = C.RES / 'inflation_sens.json'
    name = 'HFa'
    d = C.load_bench(name)
    W = waterfall_table(name)
    N0 = C.N_obs(name, d['reheating']['K']) + W['N_end']
    out = dict(N0=N0, alpha_s={}, H_i={})
    for als in (0.1171, 0.1180, 0.1189):
        r = {}
        for lab, loop in (('tree', False), ('loop', True)):
            r[f'{lab} cap'] = Mt_cap(lambda m: log10N_const(SMfull(round(m, 3), als=als), d['H_i_GeV'], 0.1, N0, loop)[0],
                                     lo=171.5, hi=176.0)
            for Mt in (172.57, 172.86):
                r[f'{lab} xi_min {Mt}'] = xi_min(lambda x: log10N_const(SMfull(Mt, als=als), d['H_i_GeV'], x, N0, loop)[0])
        out['alpha_s'][f'{als}'] = r
        print('als', als, r, flush=True)
        C.jdump(out, path)
    for Hi in (6e10, 1.0e11, 1.3e11, 2.1e11, 3e11, 5e11, 1e12):
        r = {}
        # N0 at another H_i (same K): N_* shifts by ln(H_i/H_ref)/2 ... use the bench bookkeeping: N ~ const + (1/2) ln H_i
        N0h = N0 + 0.5 * math.log(Hi / d['H_i_GeV'])
        for lab, loop in (('tree', False), ('loop', True)):
            r[f'{lab} cap'] = Mt_cap(lambda m: log10N_const(SMfull(round(m, 3)), Hi, 0.1, N0h, loop)[0], lo=171.0, hi=176.0)
            r[f'{lab} xi_min 172.57'] = xi_min(lambda x: log10N_const(SMfull(172.57), Hi, x, N0h, loop)[0])
            r[f'{lab} xi_min 172.86'] = xi_min(lambda x: log10N_const(SMfull(172.86), Hi, x, N0h, loop)[0])
        r['N0'] = N0h
        out['H_i'][f'{Hi:.3g}'] = r
        print('H_i', Hi, r, flush=True)
        C.jdump(out, path)


if __name__ == '__main__':
    {'main': cmd_main, 'sens': cmd_sens}[sys.argv[1] if len(sys.argv) > 1 else 'main']()
