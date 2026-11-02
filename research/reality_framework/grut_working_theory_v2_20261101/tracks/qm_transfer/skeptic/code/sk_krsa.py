"""qm_transfer SKEPTIC: convergence of P_mix/P_LCDM in the radiation-streaming switch K_rsa (own driver; imports the
track's code read-only).  The track tested only K_rsa = 200 vs 400; here 200 / 400 / 800 / 1600 with K_end = 2 K_rsa.
usage: sk_krsa.py OUT.json zeta SQ mode k1,k2,... K1,K2,...   (appends to OUT.json)"""
import json, os, sys, time, math
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'code'))
from qt_bg import Cosmo
from qt_hidden import Vapour
from qt_boltz import Mode, CDM

out_fn = sys.argv[1]; z = float(sys.argv[2]); SQ = float(sys.argv[3]); mode = sys.argv[4]
ks = [float(x) for x in sys.argv[5].split(',')]; Ks = [float(x) for x in sys.argv[6].split(',')]
res = json.load(open(out_fn)) if os.path.exists(out_fn) else {}
c = Cosmo()
v = Vapour(c, SQ=SQ, Q=z**-3, s=1.0, mode=mode, Nq=24, Lv=10)
key = 'SQ%g_z%g_%s' % (SQ, z, mode)
res.setdefault(key, {})
for K in Ks:
    for k in ks:
        t0 = time.time()
        kw = dict(K_rsa=K, K_end=2 * K, rtol=1e-7)
        dm = Mode(c, v, k, **kw).run()[1.0]
        dc = Mode(c, CDM(), k, **kw).run()[1.0]
        res[key]['K%g_k%g' % (K, k)] = dict(P=(dm / dc)**2, d_cdm=dc, d_dm=dm, t=time.time() - t0)
        print(key, K, k, (dm / dc)**2, dc, '%.0fs' % (time.time() - t0), flush=True)
        json.dump(res, open(out_fn, 'w'), indent=1)
