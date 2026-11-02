"""qm_transfer validation 0 (own): LCDM transfer shape of the Boltzmann code vs the Eisenstein & Hu 1998 no-wiggle fit
(astro-ph/9709112) and vs the RSA switch (K_rsa = 200 production, 400, 800).  T(k) ~ delta_m(k, z=0)/k^2 normalised at
k = 0.002/Mpc.  [C]   usage: val_lcdm.py OUT.json"""
import json, os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from qt_bg import Cosmo
from qt_boltz import Mode, CDM
import constraints as C
c = Cosmo(); out = {}
d0 = Mode(c, CDM(), 0.002, rtol=1e-7).run()[1.0]
for K in (200.0, 400.0, 800.0):
    row = {}
    for k in (0.05, 0.2, 1.0, 3.0, 10.0, 30.0, 100.0):
        d = Mode(c, CDM(), k, K_rsa=K, K_end=2 * K, rtol=1e-7).run()[1.0]
        T = (d / k**2) / (d0 / 0.002**2)
        row[str(k)] = dict(T=T, T_EHnw=float(C.eh_nw(k) / C.eh_nw(0.002)), ratio=T / float(C.eh_nw(k) / C.eh_nw(0.002)))
    out['K_rsa_%g' % K] = row
    print(K, {k: round(v['ratio'], 4) for k, v in row.items()}, flush=True)
json.dump(out, open(sys.argv[1], 'w'), indent=1)
