#!/usr/bin/env python3
"""e2e_kinetic: single half-period transition probabilities pp(k; A) = p_-(k) + p_+(k) of the axial channels [C].

Exact two-level evolution (Bloch rotation, ek_modes.c 'table' mode) of both helicities from the free vacuum over one
half-period [0, pi] of phi = S cos(m t) (mu5 = -A sin t, units m = 1), no expansion: the per-period transition
probability of each helicity in the DEPHASED model (Stueckelberg phases between half-periods randomised).
Regression: the h = -1 half-period energy reproduces mf_hybrid's axial_check.json to 5 digits (ek_modes table mode).
Smoothing: moving average over w = min(2, 0.25 A + 0.3) m (several Stueckelberg fringes; the fringes drift with A and
k between half-periods), sampled on x = k / ((1.25 A + 3 M/m + 2) m) in [0, 1], 400 points.
Usage: python3 -B ek_ptable.py   -> $SP/ptable_<Mh>.npz, ../results/ptable_summary.json
"""
import math, subprocess, sys
import numpy as np
import ek_common as C

NX = 400
AGRID = np.round(0.25 * 2.0 ** (np.arange(0, 20) / 2.0), 6)          # 0.25 ... 181


def read_table(fn):
    out = {}
    L = open(fn).read().split('\n')
    i = 0
    while i < len(L):
        if L[i].startswith('A '):
            _, A, nk = L[i].split(); nk = int(nk)
            out[float(A)] = np.array([list(map(float, x.split())) for x in L[i + 1:i + 1 + nk]])
            i += nk + 1
        else:
            i += 1
    return out


def smooth(k, y, w):
    dk = k[1] - k[0]
    n = max(1, int(round(w / dk)))
    if n < 2:
        return y
    ker = np.ones(n) / n
    return np.convolve(y, ker, mode='same')


def build(Mh):
    C.compile_engine()
    fn = C.SP / ('tab_%.5f.txt' % Mh)
    subprocess.run([str(C.BIN), 'table', str(fn), '%.10g' % Mh] + ['%.10g' % A for A in AGRID], check=True)
    T = read_table(fn)
    lA, pp, summ = [], [], []
    for A in sorted(T):
        a = T[A]
        k, pm, ppl = a[:, 0], a[:, 1], a[:, 2]
        tot = pm + ppl
        w = min(2.0, 0.25 * A + 0.3)
        sm = smooth(k, tot, w)
        ks = 1.25 * A + 3 * Mh + 2
        x = np.linspace(0, 1, NX)
        pp.append(np.interp(x * ks, k, sm, left=sm[0], right=0.0))
        lA.append(math.log(A))
        om = np.sqrt(k ** 2 + Mh ** 2)
        # per-period energy (both helicities, one species) and band capacity at n = 1/2 (dephased saturation)
        e_raw = np.trapezoid(k ** 2 * tot * om, k) / (2 * math.pi ** 2)
        e_sm = np.trapezoid(k ** 2 * sm * om, k) / (2 * math.pi ** 2)
        summ.append(dict(A=A, eps_per_period=e_raw, eps_smoothed=e_sm, h_plus_share=float(np.trapezoid(k ** 2 * ppl * om, k) /
                    max(np.trapezoid(k ** 2 * tot * om, k), 1e-300)), pmax=float(tot.max())))
    np.savez(C.SP / ('ptable_%.5f.npz' % Mh), lA=np.array(lA), pp=np.array(pp), A=np.array(sorted(T)), Mh=Mh)
    return summ


if __name__ == '__main__':
    Mhs = [float(x) for x in sys.argv[1:]] or [250 / 600., 400 / 250., 200 / 250., 120 / 250.]
    out = dict(label='C', note=__doc__, tables={})
    for Mh in Mhs:
        s = build(Mh)
        out['tables']['%.5f' % Mh] = s
        print('Mh %.4f' % Mh, ['A %.3g e %.3g' % (r['A'], r['eps_per_period']) for r in s[::3]], flush=True)
    import json, os
    fn = C.RES / 'ptable_summary.json'
    old = json.load(open(fn)) if fn.exists() else dict(tables={})
    old.update(dict(label=out['label'], note=out['note']))
    old['tables'].update(out['tables'])
    C.dump(fn, old)
