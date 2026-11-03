"""Self-consistency of the calibrated closure: on the background that the calibrated closure itself produced (KFmid),
compare the line fermion energy (lfline) with the passive calibrated closure (kappa_j(t) table, local blocking, no
decays) on the same lines, and with P0 / GM (lf3d tracers).  Claude, lf_fermion (issue #13 context).
Usage: python3 -B check_kf.py RUN VAR -> ../results/check_<RUN>.json"""
import json, math, sys
import numpy as np
sys.dont_write_bytecode = True
import lines_io, lfcommon as LC
import calib_seq as CS
import kappa_table as KT

SPI = {'N1': 0, 'N2': 1, 'N3': 2, 'N4': 3}
run = sys.argv[1] if len(sys.argv) > 1 else 'KFmid'
var = sys.argv[2] if len(sys.argv) > 2 else 'mid'
rd = LC.RUNS / run
D = lines_io.read_lines(rd)
cmpd = json.loads((LC.TRACK / 'results' / ('cmp_%s.json' % run)).read_text())
tg = np.array(cmpd['t'])
T = KT.table('REC_P0')[var]
out = dict(run=run, var=var, label='C', note=__doc__, t=list(tg), species={})
for nm, d in cmpd['species'].items():
    j = SPI[nm]
    t, H, ch = CS.prep(rd, D, d['lines'], 1.1)
    kh, kf = T[nm]['k_hom'], T[nm]['k_frag']
    B = np.zeros(ch[j]['G'].shape[1:]); TE = np.zeros_like(B)
    te = np.zeros(len(t))
    for n in range(len(t) - 1):
        tt = t[n]
        kap = kh if tt <= 0.45 else (kf if tt >= 0.55 else kh + (kf - kh) * (tt - 0.45) / 0.1)
        B, TE = CS.advance(ch[j], t, H, n, n + 1, kap, B, TE)
        te[n + 1] = TE.mean()
    kfc = np.interp(tg, t, te)
    q = np.array(d['Q_mean'])
    out['species'][nm] = dict(Q=list(q), Q_sem=d['Q_sem'], closure_kappa=list(kfc), P0=d['tracer_mean']['0.0'], GM=d['tracer_mean']['1.0'],
                              ratio_Q_closure=list(q / kfc), lambda_out=d['lambda_out'])
    print(nm, ' '.join('t%.2f Q/KF %.2f (Q/P0 %.2f Q/GM %.2f)' % (tg[i], q[i] / kfc[i], d['ratio_Q_P0'][i], d['ratio_Q_GM'][i]) for i in (5, 7, 9, 13, 17, 21)))
LC.dump(LC.TRACK / 'results' / ('check_%s.json' % run), out)
