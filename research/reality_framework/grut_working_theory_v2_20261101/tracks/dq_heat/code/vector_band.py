#!/usr/bin/env python3
"""dq_heat (Claude, own code): resolve the release-era (piece B) resonance band of the longitudinal mode.

The X-orbit modulation of the Higgs norm, a^2 F^2 = 2 chi^2(eta) + 72 |Sigma|^2 a^2 + ..., pumps the relativistic
longitudinal (Goldstone-regime) mode near half the radial frequency of the eccentric X orbit; vector_modes.py samples it
with one point (k = 4.11).  Here: a dense grid over k in [1, 30] (variants 'full' and 'noX'), integrated in ln k, replacing
the coarse k > 1 part of vector_modes.json.  In the Goldstone regime n_k is e-independent (checked at e = 1e-6 and 1e-5).
Usage: python3 -B vector_band.py   (DQ_BG=planck for the other renormalisation condition)
"""
import json, sys, time
from multiprocessing import Pool
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import vector_modes as V

RES = Path(__file__).resolve().parent.parent / 'results'


def main():
    t0 = time.time()
    V.load_bg()
    kb = np.unique(np.concatenate([np.geomspace(1.0, 2.0, 8), np.linspace(2.0, 8.0, 121), np.geomspace(8.0, 30.0, 16)]))
    jobs = [(k, 1e-6, 'full') for k in kb] + [(k, 1e-6, 'noX') for k in kb] + [(k, 1e-5, 'full') for k in kb[::8]]
    if V.BGTAG != 'm1r':
        jobs = [(k, 1e-6, 'full') for k in kb]
    with Pool(2) as pool:
        rows = pool.map(V.job, jobs, chunksize=4)
    out = dict(method=__doc__, bg=V.BGTAG, kgrid=kb.tolist(), sets={})
    for var, e in (('full', 1e-6), ('noX', 1e-6), ('full', 1e-5)):
        rr = sorted([r for r in rows if r['variant'] == var and r['e'] == e], key=lambda r: r['k'])
        if not rr:
            continue
        k = np.array([r['k'] for r in rr]); n = np.array([r['n_k'] for r in rr])
        integ = k ** 3 * n / (2 * np.pi ** 2)
        N = float(np.trapezoid(integ, np.log(k)))
        ip = int(np.argmax(n))
        # resonance band: where n > 10 x the smooth floor (median of the outer points)
        floor = float(np.median(n[(k < 2) | (k > 8)]))
        band = n > 10 * floor
        out['sets']['%s_%.0e' % (var, e)] = dict(variant=var, e=e, N_A_k_ge_1=N, k_peak=float(k[ip]), n_peak=float(n[ip]),
                                                 band_k=[float(k[band].min()), float(k[band].max())] if band.any() else None,
                                                 N_band=float(np.trapezoid(np.where(band, integ, 0), np.log(k))),
                                                 spectrum=[[float(a), float(b)] for a, b in zip(k, n)])
        print('%-5s e=%.0e N_A(k>=1)=%.5f peak k=%.3f n=%.3e band %s N_band=%.5f' % (var, e, N, k[ip], n[ip],
              out['sets']['%s_%.0e' % (var, e)]['band_k'], out['sets']['%s_%.0e' % (var, e)]['N_band']), flush=True)
    out['seconds'] = time.time() - t0
    json.dump(out, open(RES / ('vector_band.json' if V.BGTAG == 'm1r' else 'vector_band_planck.json'), 'w'), indent=1)
    print('seconds', out['seconds'])


if __name__ == '__main__':
    main()
