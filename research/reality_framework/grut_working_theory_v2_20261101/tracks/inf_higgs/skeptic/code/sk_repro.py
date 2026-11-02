#!/usr/bin/env python3
"""Skeptic of inf_higgs: reproduction of R4 (R < 0 windows), R7 (deadlines, fine tau_th grid), R8 (lower xi edge) with
own code.  -> results/sk_repro.json"""
import sys
sys.dont_write_bytecode = True
import math

import numpy as np

import sk_common as C
import sk_bench as B


def windows(ev, n=6):
    F = ev['fine']
    t, R, H = F[:, 0], F[:, 1], F[:, 4]
    te = ev['tau_end_inf']
    m = t > te - 0.05
    t, R, H = t[m], R[m], H[m]
    neg = R < 0
    out = []
    i = 0
    while i < len(t) and len(out) < n:
        if neg[i]:
            j = i
            while j < len(t) and neg[j]:
                j += 1
            k = i + int(np.argmin(R[i:j]))
            out.append(dict(tau=float(t[k]), R_min_over_H2=float(R[k] / H[k] ** 2), duration_H=float((t[j - 1] - t[i]) * H[k]),
                            exponent=float(np.trapezoid(np.sqrt(np.maximum(-ev['xi'] * R[i:j], 0)), t[i:j]))))
            i = j
        else:
            i += 1
    post = (t > te) & (t < te + 3)
    return dict(first=out, frac_R_neg_first3=float(np.mean(R[post] < 0)))


def deadlines(ev, Mts=(172.90, 172.95, 173.05, 173.15, 173.45)):
    ev['Rbar'] = C.period_mean(ev, 1)
    taus = np.arange(2.0, 40.01, 0.25)
    out = {}
    for Mt in Mts:
        smi = C.sm_in(Mt)
        Ds = np.array([C.danger(ev, smi, ev['Rbar'], Ns=(60.0,), tau_th=tt)[60.0]['D'] for tt in taus])
        Df = np.where(np.isfinite(Ds), Ds, -1e9)
        i = np.where(Df >= 0)[0]
        if len(i) and i[0] > 0:
            a, b = i[0] - 1, i[0]
            d = taus[a] + (0 - Df[a]) * (taus[b] - taus[a]) / (Df[b] - Df[a])
        else:
            d = None
        out[str(Mt)] = d
        print('deadline', Mt, d, flush=True)
    return out


if __name__ == '__main__':
    o = {}
    for name, hist in (('HT6t', 'MD'), ('T1', 'MD'), ('E1', 'MD'), ('M1r', 'G05')):
        ev = B.load(B.key(name, hist, 0))
        w = windows(ev)
        o[f'windows {name}'] = w
        print(name, 'first windows', [(round(x['R_min_over_H2'], 2), round(x['duration_H'], 4), round(x['exponent'], 4)) for x in w['first'][:3]],
              'frac', round(w['frac_R_neg_first3'], 3), flush=True)
    for name, hist in (('HT6t', 'MD'), ('HT6t', 'G05')):
        o[f'deadline {name} {hist}'] = deadlines(B.load(B.key(name, hist, 0)))
    # lower edge: xi_EW = 0.095 at 172.86 and 0.067 at 172.57 (track's MD edges)
    for xe, Mt in ((0.095, 172.86), (0.067, 172.57)):
        ev = B.run('HT6t', 'MD', 0, xi_EW=xe, tag=f'_xe{xe}')
        D = C.danger(ev, C.sm_in(Mt), ev['Rbar'], Ns=(60.0,))[60.0]['D']
        o[f'lower edge xi_EW={xe} Mt={Mt}'] = D
        print('xi_EW', xe, 'Mt', Mt, 'D60', D, flush=True)
    C.jdump(o, C.RES / 'sk_repro.json')
