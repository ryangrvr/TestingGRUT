#!/usr/bin/env python3
"""hw_verify (Claude, own code): reproduction of ls_scaling's heavy ladder with a different integrator.

Ladder: M_Y/H = 10, 15, 20, 30, 40 at r = 0.45, h = 40 H, edge fraction 0.5 (kappa = 1.153e-4 m), c_H = 10,
Y-relaxed preparation, delta = pi/2 and pi/4, H_i = 4.3487e11 GeV; instant radiation release; a_end = 100,
window a = 30-100 (ls_scaling: a_end = 40, window 30-40).  Controls at the spot point (m = 30, pi/2): Yoshida
c = 0.1/0.05/0.025, DP5(4) rtol 1e-11/1e-13, a_end = 300, the yield-normalised GeV point, Codex preparation.
Usage: python3 -B ladder.py ../results/ladder.json     (~3 min on 4 cores)
"""
import json, sys, time
from multiprocessing import Pool
from pathlib import Path
import numpy as np
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import hw

H = hw.H_CAND
S_NORM = 5.889090163799063**(-2 / 3)      # scale factor that brings the raw yield at H_CAND to 1.00 Y_Q


def job(spec):
    name, cp, Hs, delta, prep, kw = spec
    t0 = time.time()
    z0, info = (hw.prep_relaxed if prep == 'relaxed' else hw.prep_codex)(cp, Hs, hw.CH, delta)
    a_end = kw.pop('a_end', 100.0); win = kw.pop('window', '30-100')
    r = hw.run(cp, Hs, z0, a_end=a_end, windows=(1, 10, 30, 100, 300), **kw)
    w = r['windows'][win]
    out = hw.summarize(cp, Hs, w, prep_info=info,
                       extra=dict(name=name, delta=delta, prep=prep, window=win, a_end=a_end, controls=r['end'],
                                  method=kw.get('method', 8), cpar=kw.get('cpar', 0.05),
                                  windows_all={k: {q: v[q] for q in ('eps', 'NX', 'EX', 'NYH', 'NYL', 'Nsig', 'Gauss')}
                                               for k, v in r['windows'].items()},
                                  codex_envelope=hw.codex_envelope_bound(cp, Hs, hw.CH, info['x']),
                                  seconds=time.time() - t0))
    print(name, delta, prep, 'eps %.6f R %.4g RL %s %.0fs' % (out['eps'], out['R_central'], out['R_L'], out['seconds']),
          flush=True)
    return out


def specs():
    sp = []
    for m in (10.0, 15.0, 20.0, 30.0, 40.0):
        for d in (np.pi / 2, np.pi / 4):
            sp.append((f'm{m:g}', hw.heavy(m=m), H, d, 'relaxed', {}))
    cp = hw.heavy()
    sp.append(('spot_c0.1', cp, H, np.pi / 2, 'relaxed', dict(cpar=0.1)))
    sp.append(('spot_c0.025', cp, H, np.pi / 2, 'relaxed', dict(cpar=0.025)))
    sp.append(('spot_dp5_1e-11', cp, H, np.pi / 2, 'relaxed', dict(method=5, cpar=1e-11)))
    sp.append(('spot_dp5_1e-13', cp, H, np.pi / 2, 'relaxed', dict(method=5, cpar=1e-13)))
    sp.append(('spot_a300', cp, H, np.pi / 2, 'relaxed', dict(a_end=300.0, window='100-300')))
    sp.append(('spot_codexprep', cp, H, np.pi / 2, 'codex', {}))
    sp.append(('spot_codexprep_a300', cp, H, np.pi / 2, 'codex', dict(a_end=300.0, window='100-300')))
    sp.append(('spot_Hnorm', hw.heavy(Hs=S_NORM * H), S_NORM * H, np.pi / 2, 'relaxed', {}))
    return sp


if __name__ == '__main__':
    out = Path(sys.argv[1])
    sp = specs()
    sp.sort(key=lambda s: -s[5].get('a_end', 100.0)**2 * s[1]['M_GeV'] / s[2])
    rows = []
    with Pool(4) as p:
        for r in p.imap_unordered(job, sp):
            rows.append(r)
            out.write_text(json.dumps(dict(method=__doc__, H_GeV=H, cH=hw.CH, s_norm=S_NORM, rows=rows), indent=1,
                                      default=float))
    print('rows', len(rows))
