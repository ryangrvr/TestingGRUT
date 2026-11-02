"""Numerical controls for modes.py (k-grid, K_max, tolerance, save cadence) at the benchmark and at the two edges of
the gravitational-only window, M_t = 172.6, H_i = 1.29e11 GeV.  Writes ../results/controls.json."""
import json
import os

import modes as M

HERE = os.path.dirname(os.path.abspath(__file__))
lam = M.Lam(172.60)
out = {}
for xi in (0.0567, 0.08775, 30.8):
    rows = {}
    for lab, kw in (('base n=260 rtol=1e-8 nsave=800', dict(nsave=800)),
                    ('n=400', dict(n=400, nsave=800)),
                    ('Kmax x1.5', dict(Kmax=1.5 * (6.0 + 4.0 * (12 * xi) ** 0.5), nsave=800)),
                    ('rtol=1e-10', dict(rtol=1e-10, atol=1e-14, nsave=800)),
                    ('nsave=1600', dict(nsave=1600)),
                    ('no Hartree', dict(hartree=False, nsave=800))):
        r, _ = M.run(xi, lam=lam, **kw)
        b = r['danger']['best'][60]
        rows[lab] = dict(D60=b['D'], tau_c=b.get('tau_c'), sigma2_c=b.get('sigma2_c'), max_sigma2=r['max_sigma2'],
                         sigma2_end=r['sigma2_end'], wronskian_err=r['wronskian_err'])
        print(xi, lab, rows[lab], flush=True)
    out[f'xi_H={xi}'] = rows
json.dump(out, open(os.path.join(HERE, '..', 'results', 'controls.json'), 'w'), indent=1, default=float)
