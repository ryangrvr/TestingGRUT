"""qm_transfer SKEPTIC: the kinetic-vs-two-fluid comparison (track R3) at a converged radiation-streaming switch.
Same construction as the track's val_twofluid.py (fx_observables history_B, imported read-only) but K_rsa = 800,
K_end = 1600 for both the vapour and the CDM reference.  [C]
usage: sk_twofluid_K.py OUT.json Q k1,k2,..."""
import importlib.util, json, math, os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
TR = os.path.join(HERE, '..', '..', 'code'); sys.path.insert(0, TR)
from qt_bg import Cosmo
from qt_hidden import Vapour
from qt_boltz import Mode, CDM
FXO = os.path.abspath(os.path.join(HERE, '..', '..', '..', 'fx_observables'))
spec = importlib.util.spec_from_file_location('tfc', os.path.join(FXO, 'skeptic', 'code', 'two_fluid_check.py'))
tfc = importlib.util.module_from_spec(spec); spec.loader.exec_module(tfc)
vt = tfc.vt
out_fn = sys.argv[1]; Q = float(sys.argv[2]); ks = [float(x) for x in sys.argv[3].split(',')]
SQ, s = 1.235, 1.0
res = json.load(open(out_fn)) if os.path.exists(out_fn) else {}
c = Cosmo(h=0.674, omb=0.0224, omc=0.12, Neff=3.044)
Hf = tfc.history_full(Q, SQ, s, sigma_m=None)
p = vt.obj_props(Q); R = p['R']; sig_m = math.pi * R * R
fxh = dict(a=Hf['a'], Nv=Hf['Nv'], Tv=Hf['Tv'], af=Hf['af'], am=Hf['track']['a_m'], sig_c=s * math.pi * R * R, sig_m=sig_m, Q=Q, NQ0=vt.NQ0)
v = Vapour(c, SQ=SQ, Q=Q, s=s, mode='fx', fx_hist=fxh, NQ0=vt.NQ0, Delta=1 - p['Omega'], Nq=28, Lv=10,
           sig_a=s * math.pi * R * R, sig_el=sig_m - s * math.pi * R * R)
row = res.setdefault('%g' % Q, {})
for k in ks:
    for K in (200.0, 800.0):
        t0 = time.time()
        kw = dict(K_rsa=K, K_end=2 * K, rtol=1e-7)
        P = (Mode(c, v, k, **kw).run()[1.0] / Mode(c, CDM(), k, **kw).run()[1.0])**2
        row['K%g_k%g' % (K, k)] = P
        print('Q=%.0e k=%g K=%g kinetic P=%.4f (%.0fs)' % (Q, k, K, P, time.time() - t0), flush=True)
        json.dump(res, open(out_fn, 'w'), indent=1)
