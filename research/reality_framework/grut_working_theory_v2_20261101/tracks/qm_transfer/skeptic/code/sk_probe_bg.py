"""qm_transfer SKEPTIC: post-freeze-out vapour conditions (n_v, T_v, T_o) from the track's own qt_hidden.Vapour
(read-only import) on the physical line S/Q = 1.235, s = 1, 'thermo' and 'sudden' -> input of sk_object_mc.py.
usage: sk_probe_bg.py OUT.json"""
import sys, math, json, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..', '..', 'code'))
from qt_bg import Cosmo
from qt_hidden import Vapour
c = Cosmo(recomb=False)
out = {}
for z in (1e-2, 3e-3, 1e-3, 3e-4, 1e-4):
    for mode in ('thermo', 'sudden'):
        v = Vapour(c, SQ=1.235, Q=z**-3, s=1.0, mode=mode, Nq=24, Lv=10)
        af = v.a_f
        rows = []
        for fac in (1, 2, 5, 20, 100):
            j = int(np.argmin(abs(v.k_lna - math.log(af * fac))))
            N = v.k_N[j]; Nv = float(v.wq @ N); q2 = float(v.wq @ (N * v.q**2)) / Nv
            aa = math.exp(v.k_lna[j])
            rows.append(dict(fac=fac, a=aa, T_o=float(v.k_T[j]), Nv=Nv, n_v=Nv * v.NQ0 / aa**3, Tv_rms=q2 / 3 / aa**2, H=v.H_mX(aa)))
        out['z%g_%s' % (z, mode)] = dict(Q=z**-3, a_f=af, a_f_eq=v.a_f_eq, T_f_eq=v.T_f_eq, sigma_p=v.sigma_p, NQ0=v.NQ0,
                                         Delta=v.Delta, f_v_charge=v.f_v_charge, rows=rows)
json.dump(out, open(sys.argv[1], 'w'), indent=1)
