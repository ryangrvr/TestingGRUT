#!/usr/bin/env python3
"""ls_critic (Claude): locate the ke (energy-conserving) case-L crossing at the ls_scaling heavy point more finely,
with a step-size control on the steep edge and the on-shell s route.  Same inputs as heavy_threshold.py
(eps = 0.14571744, yield normalised); imports tracks/ls_threshold/code/lsrelic.py unchanged.
Usage: OMP_NUM_THREADS=1 python3 -B ke_refine.py <out.json>     (~15 min on 4 cores)
"""
import json, sys
from multiprocessing import Pool
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
from heavy_threshold import ke, W, A  # noqa: E402

if __name__ == '__main__':
    jobs = [(R, W, False, 500) for R in (6.0, 7.0, 8.0)]
    jobs += [(5.0, W, False, 1000), (8.9, W, False, 1000), (7.0, W, False, 1000)]
    jobs += [(R, W, True, 500) for R in (20.0, 25.0)]
    jobs += [(5.0, A, False, 500)]
    rows = []
    with Pool(4, maxtasksperchild=1) as p:
        res = [(a, p.apply_async(ke, (a,))) for a in jobs]
        for a, r in res:
            try:
                rows.append(r.get())
            except Exception as exc:
                rows.append(dict(kind='failed', args=list(a), status=repr(exc)))
            Path(sys.argv[1]).write_text(json.dumps(dict(method=__doc__, rows=rows), indent=1, default=float))
            x = rows[-1]
            if x['kind'] == 'ke':
                print('ke R=%.4g mu_c=%.4f s=%s N=%d f=%.4g fallbacks=%d' % (x['R'], x['mu_c'], x['s_channel'], x['N'],
                                                                         x['f'], x['root_fallbacks']), flush=True)
            else:
                print(x, flush=True)
