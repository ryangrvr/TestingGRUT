#!/usr/bin/env python3
"""ls_threshold driver, own code: where the f = 0.01 crossing lies when the hidden state keeps the actual
X-bar number (ke mode), candidate 10 at its raw yield.  Coarse scan, not a bisection.

Usage: OMP_NUM_THREADS=1 python3 -B run_ke_scan.py <out.json>
"""
import json, sys
from multiprocessing import Pool
from pathlib import Path
import lsrelic as LR
from run_proxy import C


def job(R):
    try:
        r = LR.ke_run(C[10]['eps'], R, C[10]['ym'], N=500, antilumps=True)
    except Exception as exc:              # recorded, not hidden
        return dict(R=R, f=float('nan'), sigma0=float('nan'), status=repr(exc))
    r.pop('history')
    return r


if __name__ == '__main__':
    Rs = [0.0, 10.0, 20.0, 35.0, 50.0, 70.0, 90.0, 103.0]
    with Pool(4) as p:
        rows = p.map(job, Rs)
    Path(sys.argv[1]).write_text(json.dumps(dict(method=__doc__, rows=rows), indent=1, default=float))
    for r in rows:
        print('R', r['R'], 'f=%.4g' % r['f'], 'sigma0 %.3g' % r['sigma0'])
