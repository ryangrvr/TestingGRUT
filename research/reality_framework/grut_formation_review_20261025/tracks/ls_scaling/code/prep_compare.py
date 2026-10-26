#!/usr/bin/env python3
"""ls_scaling (Claude, own code): preparation dependence of the stated-coupling debris at Codex's candidates.

Codex preparation (b448f25 / 373d37e): aligned radial inflationary minimum, then Sigma -> -|Sigma| e^{i delta}
with Y left at its delta = 0 value.  Y-relaxed preparation (our variant): conditional minimum over |X|, complex Y
and |Sigma| at the same fixed phases (a heavy Y tracks it during inflation).  Nominal phases only, a_end = 300,
window a = 100-300 (control: candidate-10 R within 0.03% of the a = 1000 value).
Usage: python3 -B prep_compare.py <out.json>      (~5 min on 4 cores)
"""
import json, sys
from multiprocessing import Pool
from pathlib import Path
import numpy as np

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import scal as S  # noqa: E402

H, CH = S.H_CAND, 10.0
CANDS = {10: np.pi / 6, 17: 5 * np.pi / 24}


def job(a):
    cand, prep = a
    cp = dict(S.STATED); d = CANDS[cand]
    if prep == 'codex':
        z0, m, ms = S.initial(cp, H, CH, d); info = dict(x=m['x'], y=m['y'], s=m['s'], n_minima=len(ms))
    else:
        z0, info = S.relaxed_initial(cp, H, CH, d)
    rr = S.run(cp, H, z0, a_end=300.0, rtol=1e-11, atol=1e-14, nsamp=60)
    w = rr['windows']['100-300']
    L = S.ledger(cp, H, w['NYH'], w['NYL'], w['Nsig'], w['EX'])
    return dict(candidate=cand, delta=d, prep=prep, eps=w['eps'], NX=w['NX'], EX=w['EX'], NYH=w['NYH'], NYL=w['NYL'],
                Nsig=w['Nsig'], R_Y=L['R_Y'], R_Sigma=L['R_Sigma'], R_cascade=L['R_cascade'], R_lower=L['R_lower'],
                R_central=L['R_central'], R_L_normalised=S.R_L(w['eps']), yield_over_YQ=S.yield_over_YQ(w['NX'], H),
                prep_info=info, controls=rr['end'],
                history=[dict(a=s['a'], eps=s['eps'], NYH=s['NYH'], NYL=s['NYL'], Nsig=s['Nsig']) for s in rr['samples']])


if __name__ == '__main__':
    with Pool(4) as p:
        rows = p.map(job, [(c, pr) for c in CANDS for pr in ('codex', 'relaxed')])
    Path(sys.argv[1]).write_text(json.dumps(dict(method=__doc__, rows=rows), indent=1, default=float))
    for r in rows:
        print(r['candidate'], r['prep'], 'eps %.5f NYH %.3f NYL %.3f Nsig %.3f R_lower %.2f R_central %.2f R_L %.1f'
              % (r['eps'], r['NYH'], r['NYL'], r['Nsig'], r['R_lower'], r['R_central'], r['R_L_normalised']))
