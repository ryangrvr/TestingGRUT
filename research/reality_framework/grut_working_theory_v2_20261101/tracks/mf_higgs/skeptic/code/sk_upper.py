#!/usr/bin/env python3
"""mf_higgs SKEPTIC: the upper xi edge (mean-field 2 sigma_prod = phi_max, inf_higgs convention: produced occupation at
t_end mapped back adiabatically, maximised over t > t(eps_H = 1) + 0.5) with my own modes on the bench export.
[C on the argued criterion]  Output: results/sk_upper.json"""
import sys
sys.dont_write_bytecode = True
import math
import numpy as np
import sk_common as C
import sk_modes as M

out = {}
for xi in (10.0, 30.0, 100.0, 300.0):
    Kmax = 6 + 4 * math.sqrt(12 * xi)
    o = M.run('EXP', xi, t0=-1.0, t_end=20.0, n=320, Kmax=Kmax, rtol=1e-7)
    nk = np.maximum(o['nk_end'], 0)
    t, bg = o['t'], o['bg']
    late = t > o['te'] + 0.5
    sp = []
    for i in np.where(late)[0]:
        a, H, Rb = bg[i, 0], bg[i, 1], bg[i, 3]
        om = np.sqrt(o['ks'] ** 2 / a ** 2 + max(xi * Rb, 0) + H ** 2)
        sp.append(float(np.sum(o['w'] * nk / (a ** 3 * om))))
    s2 = max(sp)
    row = dict(xi_H=xi, xi_EW=C.Flat(172.57).xi_EW(o['Hi'], xi), Kmax=Kmax, max_nk=float(nk.max()),
               k_at_max_nk=float(o['ks'][int(np.argmax(nk))]), sigma2_prod=s2, seconds=o['seconds'], nfev=o['nfev'])
    for Mt in (172.28, 172.57, 172.86):
        fl = C.Flat(Mt)
        phim = math.exp(fl.Gam(o['Hi'])) * fl.h_max / o['Hi']
        row[f'2sigma_over_phimax {Mt}'] = 2 * math.sqrt(s2) / phim
    out[f'{xi:g}'] = row
    print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
    C.jdump(out, C.RES / 'sk_upper.json')
