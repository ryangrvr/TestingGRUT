#!/usr/bin/env python3
"""e2e_higgs: Higgs during inflation at the final H_i (mf_higgs mh_inflation, read-only), N0 from e2e_kinetic's K(H_i).

  * constant-H Hawking-Moss count (tree mu^2 = h^2 + R, and one-loop dS), xi_min and the M_t cap at xi_EW = 0.1;
  * the instantaneous count through the waterfall WITH the exact eps_H = 1 endpoint (mf_higgs skeptic N1): all bench
    rows (stride 1) from t = 0.10 to t(eps_H = 1), plus the exact end row (H(t_e), R = 6 H^2).  As the skeptic showed,
    this reading is not a bracket (continued 0.003 H^-1 further it diverges); it is reported only as the number the
    track's "pessimistic bracket" becomes with the endpoint fixed.
The waterfall up to eps_H = 1 is the homogeneous bench (the lattice reproduces it to 3e-4 with fluctuations off and
its fluctuations are ~1e-5 of the field there), so it is common to every final history; only H_i and N0 change.
Output: results/inflation.json.   [C] on the HM criterion [E]; model choices throughout.
"""
import sys
sys.dont_write_bytecode = True
import json
import math

import numpy as np

import eh_common as E
import mh_inflation as MI

MC = E.MC
LN10 = math.log(10)
HIS = (5e10, 1.3e11, 2.4e11)


def inst_exact(sm, H_GeV, xi_EW, N0, loop=False):
    W = MI.waterfall_table('HFa')
    b = MC.load_bench('HFa')
    te = b['t_eps1']
    h = b['h']
    He = float(np.interp(te, h['t'], h['H']))
    Nb = np.concatenate([W['Nb'], [0.0]])
    H = np.concatenate([W['H'], [He]])
    R = np.concatenate([W['R'], [6 * He ** 2]])
    Bs = np.array([MI.B_of(sm, hh * H_GeV, r * H_GeV ** 2, xi_EW, loop) if r > 0 else 0.0 for hh, r in zip(H, R)])
    B12 = MI.B_of(sm, H_GeV, 12 * H_GeV ** 2, xi_EW, loop)
    lnI = 3 * (N0 - Nb) + 3 * np.log(H) - Bs
    o = np.argsort(Nb)
    x, y = Nb[o], lnI[o]
    mx = np.max(y)
    part = mx + math.log(max(np.trapezoid(np.exp(y - mx), x), 1e-300))
    tail = 3 * (N0 - x[-1]) - B12 - math.log(3.0)
    tot = np.logaddexp(part, tail) + math.log(4 * math.pi / 3)
    return float(tot / LN10), float(Bs[-1])


def main():
    path = E.RES / 'inflation.json'
    out = {}
    W = MI.waterfall_table('HFa')
    for Hi in HIS:
        K = E.K_KIN[Hi]
        N0 = MC.N_obs('HFa', K) + 0.5 * math.log(Hi / E.H_BENCH) + W['N_end']
        r = dict(H_i=Hi, K=K, N0=N0, rows={})
        for Mt in (172.57, 172.86, 173.45, 174.05, 174.6):
            sm = MI.SMfull(Mt)
            row = {}
            for lab, loop in (('tree', False), ('loop', True)):
                row[f'{lab} const log10N xi0.1'] = MI.log10N_const(sm, Hi, 0.1, N0, loop)[0]
                row[f'{lab} const xi_min'] = MI.xi_min(lambda x: MI.log10N_const(sm, Hi, x, N0, loop)[0])
            l10, Bend = inst_exact(sm, Hi, 0.1, N0)
            row['tree inst_exact log10N xi0.1'] = l10
            row['tree inst_exact B_end'] = Bend
            if Mt in (172.57, 172.86):
                row['tree inst_exact xi_min'] = MI.xi_min(lambda x: inst_exact(sm, Hi, x, N0)[0], hi=0.3, n=13)
            r['rows'][f'{Mt}'] = row
            print(Hi, Mt, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
        caps = {}
        for lab, loop in (('tree', False), ('loop', True)):
            caps[f'{lab} const'] = MI.Mt_cap(lambda m: MI.log10N_const(MI.SMfull(round(m, 3)), Hi, 0.1, N0, loop)[0],
                                             lo=172.0, hi=176.0)
        caps['tree inst_exact'] = MI.Mt_cap(lambda m: inst_exact(MI.SMfull(round(m, 3)), Hi, 0.1, N0)[0], lo=172.5, hi=175.5)
        for als in (0.1171, 0.1189):
            caps[f'tree inst_exact als{als}'] = MI.Mt_cap(
                lambda m: inst_exact(MI.SMfull(round(m, 3), als=als), Hi, 0.1, N0)[0], lo=172.0, hi=175.5)
        r['Mt_cap_xi0.1'] = caps
        print(Hi, 'caps', caps, flush=True)
        out[f'{Hi:.3g}'] = r
        E.jdump(out, path)


if __name__ == '__main__':
    main()
