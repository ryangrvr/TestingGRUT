"""Lower edge of the gravitational-only window (N = 60) in the canonical normalisation for criterion variants:
region wavenumber k_c = f a m_c (f = 0.5, 1, 2) and a Gaussian window.  Uses sk_modes.  Writes
../results/sk_edges_band.json."""
import json
import os
import sys
import time

import numpy as np
import sk_modes as K

HERE = os.path.dirname(os.path.abspath(__file__))
sel = [int(v) for v in sys.argv[1:]] if len(sys.argv) > 1 else range(5)
out_path = os.path.join(HERE, '..', 'results', 'sk_edges_band_%s.json' % '_'.join(str(i) for i in sel))
out = {}
t0 = time.time()
for Mt, Hg in [K.corners()[i] for i in sel]:
    smi = K.SMIn(Mt)
    cache = {}

    def ev(x):
        if x not in cache:
            cache[x] = K.evolve(x)
        return cache[x]
    res = {}
    for lab, kw in (('f=0.5', dict(kc_factor=0.5)), ('f=2', dict(kc_factor=2.0)), ('gauss R=1/k_c', dict(window='gauss'))):
        D = lambda x: K.danger(ev(x), smi, Hg, norm='canon', **kw)['D']
        lo, hi = 0.0, 0.25
        dlo, dhi = D(lo + 1e-4), D(hi)
        if not (dlo > 0 > dhi):
            res[lab] = dict(note='no sign change', D_lo=dlo, D_hi=dhi)
            continue
        lo = 1e-4
        for _ in range(12):
            mid = 0.5 * (lo + hi)
            if D(mid) > 0:
                lo = mid
            else:
                hi = mid
        xe = 0.5 * (lo + hi)
        res[lab] = dict(xi_H_lower=xe, xi_EW_lower=float(smi.xi_EW(Hg, xe)))
    out[f'Mt={Mt} H={Hg:.3g}'] = res
    print(Mt, Hg, res, '%.0fs' % (time.time() - t0), flush=True)
json.dump(out, open(out_path, 'w'), indent=1, default=float)
