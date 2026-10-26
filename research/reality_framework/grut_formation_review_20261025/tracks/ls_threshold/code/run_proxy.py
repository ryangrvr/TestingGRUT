#!/usr/bin/env python3
"""ls_threshold driver, own code: benchmark (proxy-mode) thresholds, ledger relics and s-channel probe.

Usage: OMP_NUM_THREADS=1 python3 -B run_proxy.py <out.json>     (~30-40 min on 4 cores)
Ledger R values and yields are quoted [Q] from b448f25 DECAY_RESULTS.json / VECTOR_RESULTS.json.
"""
import json, sys
from multiprocessing import Pool
from pathlib import Path
import lsrelic as LR

C = {10: dict(eps=0.15741283451037602, ym=5.69009662564405, R_Y=113.0689908402781, R_lower=199.4042010902021,
              R_central=228.40581445512555 + 0.5 * (0.010909139591541852 + 0.011258527826447283),
              R_upper=266.7472531530517),
     17: dict(eps=0.1470461788202921, ym=5.309416605350421, R_Y=122.80645189533328, R_lower=249.69026892092234,
              R_central=280.5199866572872 + 0.5 * (0.010925443805502822 + 0.011275354217293932),
              R_upper=322.8328805439777)}


def th_job(a):
    cand, ym = a
    r = LR.threshold(C[cand]['eps'], ym, 40.0, 200.0, 11, N=400)
    lo, hi = r['bracket']
    # check the bracket edges at doubled grid resolution
    r['edge_check_N800'] = [dict(R=x, f=LR.proxy_run(C[cand]['eps'], x, ym, N=800)['f']) for x in (lo, hi)]
    r['candidate'] = cand
    return r


def pt_job(a):
    cand, name, ym, s_ch, N = a
    r = LR.proxy_run(C[cand]['eps'], C[cand][name], ym, N=N, s_channel=s_ch, keep=(name == 'R_lower'))
    r.update(candidate=cand, ledger=name)
    return r


if __name__ == '__main__':
    out = dict(method=__doc__)
    pts = [(k, n, ym, False, 800) for k in C for n in ('R_Y', 'R_lower', 'R_central', 'R_upper') for ym in (1.0, C[k]['ym'])]
    pts += [(k, 'R_lower', C[k]['ym'], False, 1600) for k in C]                     # grid control
    pts += [(k, n, ym, True, 800) for k in C for n in ('R_Y', 'R_lower', 'R_central') for ym in (1.0, C[k]['ym'])]
    ths = [(k, ym) for k in C for ym in (1.0, C[k]['ym'])]
    with Pool(4) as p:
        out['points'] = p.map(pt_job, pts)
        Path(sys.argv[1]).write_text(json.dumps(out, indent=1, default=float))
        out['thresholds'] = p.map(th_job, ths)
    Path(sys.argv[1]).write_text(json.dumps(out, indent=1, default=float))
    for r in out['points']:
        print(r['candidate'], r['ledger'], round(r['ymult'], 4), r['s_channel'], r['N'], 'f=%.5g' % r['f'], 'sigma=%.4f' % r['sigma'])
    for r in out['thresholds']:
        print('R_L', r['candidate'], round(r['ymult'], 4), r['bracket'], r['edge_check_N800'])
