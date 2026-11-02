"""qm_transfer SKEPTIC: is the own-code thermal relic 'colder than CLASS' a radiation-streaming-switch artefact?
Re-runs scan family W (FD thermal relic, own code, read-only import) at K_rsa = 200 (production) and 800 around the
half-mode scale and compares k_hm with the Vogel & Abazajian 2022 CLASS-based fit (arXiv:2210.10753, Table 2:
a = 0.0437, b = -1.188, nu = 1.049, theta = 2.012, eta = 0.2463, alpha in Mpc/h).  [C]
usage: sk_wdm_K.py OUT.json m_keV"""
import json, os, sys, math, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'code'))
from qt_bg import Cosmo
from qt_boltz import Mode, CDM
import scan
out_fn = sys.argv[1]; m = float(sys.argv[2])
res = json.load(open(out_fn)) if os.path.exists(out_fn) else {}
c = Cosmo()
fr, extra = scan.build(c, dict(fam='W', m_keV=m))
cW = extra['cos']
h = 0.6736
al = 0.0437 * m**-1.188
T_va = lambda k: (1 + (al * k / h)**(2 * 1.049))**(-5 / 1.049)
row = res.setdefault('m%g' % m, {})
for K in (200.0, 800.0):
    for k in (34.0, 48.0, 68.0):
        t0 = time.time()
        kw = dict(K_rsa=K, K_end=2 * K, rtol=1e-7)
        P = (Mode(cW, fr, k, **kw).run()[1.0] / Mode(c, CDM(), k, **kw).run()[1.0])**2
        row['K%g_k%g' % (K, k)] = dict(P=P, P_VA22=T_va(k)**2)
        print(m, K, k, 'P_own=%.4f  P_VA22=%.4f  (%.0fs)' % (P, T_va(k)**2, time.time() - t0), flush=True)
    ks = [34.0, 48.0, 68.0]; Ps = [row['K%g_k%g' % (K, k)]['P'] for k in ks]
    row['k_hm_K%g' % K] = float(math.exp(np.interp(0.25, Ps[::-1], np.log(ks)[::-1]))) if min(Ps) < 0.25 < max(Ps) else None
x = (2**(1.049 / 5) - 1)**(1 / (2 * 1.049)); row['k_hm_VA22'] = x / al * h
print(row)
json.dump(res, open(out_fn, 'w'), indent=1)
