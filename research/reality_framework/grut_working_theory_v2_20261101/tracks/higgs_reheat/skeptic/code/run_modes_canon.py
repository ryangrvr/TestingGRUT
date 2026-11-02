"""Benchmark xi_EW = 0.1 at the five corners with the track's modes code, normalisation-fixed copy (modes_canon.py).
Writes ../results/modes_canon_benchmark.json."""
import json
import os
import time

import modes_canon as M

HERE = os.path.dirname(os.path.abspath(__file__))
t0 = time.time()
out = {}
for Mt, H in ((172.6, 1.29e11), (172.6, 1.51e11), (172.87, 1.29e11), (172.87, 1.51e11), (172.33, 1.51e11)):
    lam = M.Lam(Mt)
    xb = float(lam.sm.xi(H, 0.1))
    r, _ = M.run(xb, Mt=Mt, H_GeV=H, lam=lam, nsave=800, cHs=(0.5, 2.0))
    b = r['danger']['best'][60]
    out[f'Mt={Mt} H={H:.3g}'] = dict(xi_H=xb, D60=b['D'], tau_c=b['tau_c'], sigma2_c=b['sigma2_c'], x=b['x'],
                                     D60_cH05=r['danger_cH']['cH=0.5'][60]['D'], D60_cH2=r['danger_cH']['cH=2.0'][60]['D'])
    print(Mt, H, out[f'Mt={Mt} H={H:.3g}'], '%.0fs' % (time.time() - t0), flush=True)
json.dump(out, open(os.path.join(HERE, '..', 'results', 'modes_canon_benchmark.json'), 'w'), indent=1, default=float)
