"""qm_transfer SKEPTIC: re-run selected production models on the full k grid at a higher radiation-streaming switch
K_rsa (K_end = 2 K_rsa), with its own CDM reference at the same switch, so that the constraint metrics can be
recomputed (imports the track's code read-only).  Restartable, time budget.
usage: sk_grid_K.py OUT.json K_rsa budget_s name1,name2,...   names like P_SQ1.235_z3e-03 or S_SQ1.235_z1e-03"""
import json, os, sys, time, math
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'code'))
from qt_bg import Cosmo
from qt_hidden import Vapour
from qt_boltz import Mode, CDM
KS = [1.5, 3.0, 5.0, 8.0, 12.0, 17.0, 24.0, 34.0, 48.0, 68.0, 96.0, 136.0]
out_fn = sys.argv[1]; K = float(sys.argv[2]); budget = float(sys.argv[3]); names = sys.argv[4].split(',')
t_start = time.time()
res = json.load(open(out_fn)) if os.path.exists(out_fn) else {}
c = Cosmo()
kw = dict(K_rsa=K, K_end=2 * K, rtol=1e-7)
ref = res.setdefault('cdm_K%g' % K, {})
def save():
    json.dump(res, open(out_fn + '.tmp', 'w'), indent=1); os.replace(out_fn + '.tmp', out_fn)
for name in names:
    fam, sq, zz = name.split('_'); SQ = float(sq[2:]); z = float(zz[1:])
    mode = {'P': 'thermo', 'S': 'sudden'}[fam]
    v = None
    row = res.setdefault('%s_K%g' % (name, K), {})
    for k in KS:
        if str(k) in row:
            continue
        if time.time() - t_start > budget:
            print('budget'); save(); sys.exit(0)
        if v is None:
            v = Vapour(c, SQ=SQ, Q=z**-3, s=1.0, mode=mode, Nq=24, Lv=10)
        t0 = time.time()
        if str(k) not in ref:
            ref[str(k)] = Mode(c, CDM(), k, **kw).run()[1.0]
        dm = Mode(c, v, k, **kw).run()[1.0]
        row[str(k)] = (dm / ref[str(k)])**2
        print(name, K, k, row[str(k)], '%.0fs' % (time.time() - t0), flush=True)
        save()
save()
