#!/usr/bin/env python3
"""st_response track (Claude, own code): finely sampled first crossings (spacing 0.002 in tau inside windows that
cover each sampled bracket and everything after the previous coarse sample below 0.8 of the level), to test whether
the 0.5-spaced samples miss an earlier narrow excursion.  Usage: python3 -B fine.py <case> ; writes scratch json.
Cases: q3 = (3, 0.5) k = 2.13 ; q5 = (5, 1) k = 2.0 and 2.13 ; q3e, q5e = whole range 0-84 at 0.002 (10% level).
"""
import json, sys
from pathlib import Path
import numpy as np
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import st_tangent as S

OUT = Path('/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/review/st_response/fine')
CASES = {'q3': dict(b=3.0, gam=0.5, ks=[2.13], windows=[(60, 140), (500, 780)], tmax=780.0),
         'q5': dict(b=5.0, gam=1.0, ks=[2.0, 2.13], windows=[(60, 150), (540, 790)], tmax=790.0),
         # early 10% warning level: Codex's sample grid below tau = 80 is sparse (MODE_PROFILES grid)
         'q3e': dict(b=3.0, gam=0.5, ks=[0.02, 2.13, 2.4], windows=[(0, 84)], tmax=84.0),
         'q5e': dict(b=5.0, gam=1.0, ks=[0.02, 2.0, 2.13], windows=[(0, 84)], tmax=84.0)}


def first(t, v, level):
    w = np.flatnonzero(np.abs(v) >= level)
    return None if not len(w) else float(t[w[0]])


def main(case):
    c = CASES[case]
    coarse = np.arange(0.0, c['tmax'] + 1e-9, 0.5)
    fine = np.concatenate([np.arange(lo, hi + 1e-9, 0.002) for lo, hi in c['windows']])
    times = np.unique(np.round(np.r_[coarse, fine], 9))
    model, sol, info, sec = S.run(c['b'], c['gam'], c['ks'], times)
    ob = S.observables(model, sol.y.T)
    t = sol.t
    norm = abs(np.trapezoid(ob['NX'][(t >= 64) & (t <= 80)], t[(t >= 64) & (t <= 80)]) / 16.0)
    if case.endswith('e'):   # also record where the sparse archived grid would have sampled (for the bracket)
        C = np.load(Path('/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/b73/research/'
                         'reality_framework/grut_source_timing_20261026') / ('quiet3.npz' if c['b'] == 3.0 else 'quiet5.npz'))
        tc = C['times'][C['times'] <= c['tmax']]
        vals = {}
        for j, k in enumerate(model.ks):
            v = 0.03 * np.abs(ob['amp'][:, j]); i = int(np.flatnonzero(v >= 0.1)[0])
            tc_before = tc[tc <= t[i]]
            vals['k=%g' % k] = dict(first_fine_0p1=float(t[i]), archived_grid_samples_around=[float(tc_before[-1]) if len(tc_before) else None,
                                    float(tc[tc > t[i]][0]) if np.any(tc > t[i]) else None],
                                    max_0p03A_before_80=float(np.max(v[t <= 80])), tau_of_that_max=float(t[t <= 80][np.argmax(v[t <= 80])]))
        print('EARLY', json.dumps(vals), flush=True)
    res = dict(case=case, spec={k: v for k, v in c.items()}, seconds=sec, norm_fine_grid=norm, rows=[])
    onc = np.isin(np.round(t, 9), np.round(coarse, 9))
    for j, k in enumerate(model.ks):
        for seed in (0.03, 1e-6):
            for name, v in (('amplitude', seed * ob['amp'][:, j]), ('local_charge', seed * ob['dN'][:, j] / norm)):
                for level in (0.1, 1.0):
                    if name == 'local_charge' and level != 1.0:
                        continue
                    res['rows'].append(dict(k=float(k), seed=seed, observable=name, level=level,
                                            first_fine=first(t, v, level), first_coarse_0p5=first(t[onc], v[onc], level)))
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / (case + '.json')).write_text(json.dumps(res, indent=1))
    print(json.dumps(res), flush=True)


if __name__ == '__main__':
    for c in sys.argv[1:]:
        main(c)
