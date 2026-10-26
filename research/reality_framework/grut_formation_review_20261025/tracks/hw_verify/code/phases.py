#!/usr/bin/env python3
"""hw_verify (Claude, own code): robustness of the heavy spot point across the supplied phase and the edge fraction.

Spot point (m = 30, r = 0.45, h = 40 H, frac = 0.5, c_H = 10, H_i = 4.3487e11 GeV), instant release, a_end = 100,
window 30-100.  (a) delta = k pi/36, k = 1..35, Y-relaxed preparation; delta = 3 pi/2 as the CP mirror of pi/2;
(b) delta = k pi/12 with the Codex preparation; (c) edge fraction 0.1 ... 1.5 at delta = pi/2 (relaxed).
Usage: python3 -B phases.py ../results/phases.json      (~3 min on 4 cores)
"""
import json, sys, time
from multiprocessing import Pool
from pathlib import Path
import numpy as np
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import hw

H = hw.H_CAND


def job(spec):
    name, cp, delta, prep, frac = spec
    t0 = time.time()
    try:
        z0, info = (hw.prep_relaxed if prep == 'relaxed' else hw.prep_codex)(cp, H, hw.CH, delta)
        r = hw.run(cp, H, z0, a_end=100.0, windows=(1, 10, 30, 100))
    except Exception as exc:
        return dict(name=name, delta=delta, prep=prep, frac=frac, status=repr(exc))
    w = r['windows']['30-100']
    out = hw.summarize(cp, H, w, prep_info=info, extra=dict(name=name, delta=delta, prep=prep, frac=frac,
                                                             controls=r['end'], seconds=time.time() - t0,
                                                             eps_10_30=r['windows']['10-30']['eps']))
    print(name, '%.4f' % delta, prep, frac, 'eps %.5f R %.4g RL %s' % (out['eps'], out['R_central'], out['R_L']), flush=True)
    return out


if __name__ == '__main__':
    out = Path(sys.argv[1])
    cp = hw.heavy()
    sp = [('phase', cp, k * np.pi / 36, 'relaxed', 0.5) for k in range(1, 36)]
    sp.append(('mirror', cp, 3 * np.pi / 2, 'relaxed', 0.5))
    sp += [('phase_codex', cp, k * np.pi / 12, 'codex', 0.5) for k in range(1, 12)]
    sp += [('frac', hw.heavy(frac=f), np.pi / 2, 'relaxed', f) for f in (0.1, 0.25, 0.5, 0.75, 0.9, 1.0, 1.2, 1.5)]
    rows = []
    with Pool(4) as p:
        for r in p.imap_unordered(job, sp):
            rows.append(r)
            out.write_text(json.dumps(dict(method=__doc__, H_GeV=H, rows=rows), indent=1, default=float))
    print('rows', len(rows))
