#!/usr/bin/env python3
"""kin_rates: production tabulation of the n->2 kernels Gt(s) and the 2->n cross sections sigma(s) (MC, own code).
usage: OMP_NUM_THREADS=1 python3 -B run_kernels.py <out_dir> [--quick]
Writes <out_dir>/kernels_n2.json and <out_dir>/xsec_2n.json.
Grid: sqrt s = sqrt s_th (1 + delta), delta geometric from 2e-3 to the top (400 m_X for n->2, 1e4 m_X for 2->n).
Regulators (Gamma_reg,s/m_s, Gamma_reg,X/m_X): primary (0.02, 0.05) on the full grid; secondary (0.04, 0.10) on
every other grid point with half the samples: the difference is the regulator systematic of the finite part.
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
N6 = 20000 if QUICK else 60000
N5 = 20000 if QUICK else 100000
NPTS = 6 if QUICK else 16
FREGS = ((0.02, 0.05, 1.0, 1), (0.04, 0.10, 0.5, 2))     # (freg_s, freg_X, sample fraction, grid stride)


def grid(sq_th, top_over_m, npts=NPTS):
    dmax = top_over_m * M.MX / sq_th - 1
    return sq_th * (1 + np.geomspace(2e-3, dmax, npts))


def _masses(legs):
    m_in = [K.MASS[K.species(l)] for l in legs if not K.is_out(l)]
    m_out = [K.MASS[K.species(l)] for l in legs if K.is_out(l)]
    return m_in, m_out


def run_n2(ch):
    spec = K.CHANNELS_N2[ch]
    legs = spec['legs']
    m_in, m_out = _masses(legs)
    sq_th = max(sum(m_in), sum(m_out))
    N = N6 if len(legs) == 6 else N5
    rows = []
    t0 = time.time()
    for freg, fregX, frac, stride in FREGS:
        rng = np.random.default_rng(abs(hash((ch, freg))) % (2**32))
        for sq in grid(sq_th, 400.0)[::stride]:
            r = K.kernel_n2(ch, sq, int(N * frac), rng, freg=freg, fregX=fregX)
            r['sqrt_s_over_mX'] = sq / M.MX
            rows.append(r)
            print(ch, f'{sq / M.MX:.3f}', f"{r['G']:.3e}", f"{r['err']:.1e}", f"full {r['G_full']:.3e}", flush=True)
    return ch, dict(legs=[list(l) if isinstance(l, tuple) else l for l in legs], Sin=spec['Sin'], Sout=spec['Sout'],
                    m_in=m_in, m_out=m_out, sqrt_s_th=sq_th, rows=rows, seconds=time.time() - t0)


def run_2n(ch):
    spec = K.CHANNELS_2N[ch]
    legs = spec['legs']
    m_in, m_out = _masses(legs)
    sq_th = max(sum(m_in), sum(m_out))
    N = N6 if len(legs) == 6 else N5
    rows = []
    t0 = time.time()
    for freg, fregX, frac, stride in FREGS:
        rng = np.random.default_rng(abs(hash((ch, freg))) % (2**32))
        for sq in grid(sq_th, 1e4)[::stride]:
            r = K.sigma_2n(ch, sq, int(N * frac), rng, freg=freg, fregX=fregX)
            r['sqrt_s_over_mX'] = sq / M.MX
            r['fregX'] = fregX
            rows.append(r)
            print(ch, f'{sq / M.MX:.3f}', f"{r['sigma']:.3e}", f"{r['err']:.1e}", f"full {r['sigma_full']:.3e}", flush=True)
    return ch, dict(legs=[list(l) if isinstance(l, tuple) else l for l in legs], Sin=spec['Sin'], Sout=spec['Sout'],
                    m_in=m_in, m_out=m_out, sqrt_s_th=sq_th, rows=rows, seconds=time.time() - t0)


def main():
    n2 = [c for c in K.CHANNELS_N2 if not c.startswith('ctrl')]
    n2 = sorted(n2, key=lambda c: -len(K.CHANNELS_N2[c]['legs']))   # 6-point first (longest)
    twon = sorted(K.CHANNELS_2N, key=lambda c: -len(K.CHANNELS_2N[c]['legs']))
    t0 = time.time()
    with Pool(4) as pool:
        res_n2 = dict(pool.imap_unordered(run_n2, n2, chunksize=1))
    json.dump(dict(N6=N6, N5=N5, fregs=FREGS, npts=NPTS, channels=res_n2, seconds=time.time() - t0,
                   units='Gt(s) in GeV^(2n-6)... see NOTES: dPhi_n dPhi_2 |M|^2 with |M|^2 dimensionless for n=4 (GeV^-4 x GeV^0): '
                         'Gt has dimension GeV^(2(n-2)-4) = GeV^0 for n=4, GeV^-2 for n=3; sigma in GeV^-2'),
              open(f'{OUT}/kernels_n2.json', 'w'), indent=1)
    print('n->2 done', time.time() - t0, flush=True)
    t1 = time.time()
    with Pool(4) as pool:
        res_2n = dict(pool.imap_unordered(run_2n, twon, chunksize=1))
    json.dump(dict(N6=N6, N5=N5, fregs=FREGS, npts=NPTS, channels=res_2n, seconds=time.time() - t1),
              open(f'{OUT}/xsec_2n.json', 'w'), indent=1)
    print('2->n done', time.time() - t1, flush=True)


if __name__ == '__main__':
    main()
