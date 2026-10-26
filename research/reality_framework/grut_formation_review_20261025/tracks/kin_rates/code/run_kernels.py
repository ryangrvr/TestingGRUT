#!/usr/bin/env python3
"""kin_rates: production tabulation of the n->2 kernels Gt(s) and the 2->n cross sections sigma(s) (MC, own code).
usage: python3 run_kernels.py <out_dir> [--quick]
Writes <out_dir>/kernels_n2.json and <out_dir>/xsec_2n.json.  Grid: sqrt s = sqrt s_th (1 + delta), delta geometric.
Two regulators (f_reg = 0.02 primary, 0.04 secondary) for the n->2 channels; the difference is the regulator systematic.
"""
import json
import sys
import time
from multiprocessing import Pool
import numpy as np
import model as M
import kernels_n2 as K

QUICK = '--quick' in sys.argv
OUT = sys.argv[1]
N6 = 60000 if QUICK else 300000
N5 = 60000 if QUICK else 200000
NPTS = 10 if QUICK else 26
FREGS = (0.02,) if QUICK else (0.02, 0.04)


def grid(sq_th, top_over_m=400.0, npts=NPTS):
    dmax = top_over_m * M.MX / sq_th - 1
    return sq_th * (1 + np.geomspace(2e-3, dmax, npts))


def run_n2(ch):
    spec = K.CHANNELS_N2[ch]
    legs = spec['legs']
    m_in = [K.MASS[K.species(l)] for l in legs if not K.is_out(l)]
    m_out = [K.MASS[K.species(l)] for l in legs if K.is_out(l)]
    sq_th = max(sum(m_in), sum(m_out))
    N = N6 if len(legs) == 6 else N5
    rows = []
    t0 = time.time()
    for freg in FREGS:
        rng = np.random.default_rng(abs(hash((ch, freg))) % (2**32))
        for sq in grid(sq_th):
            r = K.kernel_n2(ch, sq, N, rng, freg=freg)
            r['sqrt_s_over_mX'] = sq / M.MX
            rows.append(r)
    return ch, dict(legs=[list(l) if isinstance(l, tuple) else l for l in legs], Sin=spec['Sin'], Sout=spec['Sout'],
                    m_in=m_in, m_out=m_out, sqrt_s_th=sq_th, rows=rows, seconds=time.time() - t0)


def run_2n(ch):
    spec = K.CHANNELS_2N[ch]
    legs = spec['legs']
    m_in = [K.MASS[K.species(l)] for l in legs if not K.is_out(l)]
    m_out = [K.MASS[K.species(l)] for l in legs if K.is_out(l)]
    sq_th = max(sum(m_in), sum(m_out))
    N = N6 if len(legs) == 6 else N5
    rows = []
    t0 = time.time()
    rng = np.random.default_rng(abs(hash(ch)) % (2**32))
    for sq in sq_th * (1 + np.geomspace(1e-2, 1e4 * M.MX / sq_th, 10 if QUICK else 22)):
        r = K.sigma_2n(ch, sq, N, rng, freg=0.02)
        r['sqrt_s_over_mX'] = sq / M.MX
        rows.append(r)
    return ch, dict(legs=[list(l) if isinstance(l, tuple) else l for l in legs], Sin=spec['Sin'], Sout=spec['Sout'],
                    m_in=m_in, m_out=m_out, sqrt_s_th=sq_th, rows=rows, seconds=time.time() - t0)


def main():
    n2 = [c for c in K.CHANNELS_N2 if not c.startswith('ctrl')]
    n2 = sorted(n2, key=lambda c: -len(K.CHANNELS_N2[c]['legs']))   # 6-point first (longest)
    twon = list(K.CHANNELS_2N)
    t0 = time.time()
    with Pool(4) as pool:
        res_n2 = dict(pool.map(run_n2, n2, chunksize=1))
    json.dump(dict(N6=N6, N5=N5, fregs=FREGS, channels=res_n2, seconds=time.time() - t0),
              open(f'{OUT}/kernels_n2.json', 'w'), indent=1)
    print('n->2 done', time.time() - t0, flush=True)
    t1 = time.time()
    with Pool(4) as pool:
        res_2n = dict(pool.map(run_2n, twon, chunksize=1))
    json.dump(dict(N6=N6, N5=N5, channels=res_2n, seconds=time.time() - t1),
              open(f'{OUT}/xsec_2n.json', 'w'), indent=1)
    print('2->n done', time.time() - t1, flush=True)


if __name__ == '__main__':
    main()
